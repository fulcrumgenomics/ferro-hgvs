//! Unit 2: build the DAG of every minimal-cost alignment.
//!
//! Global alignment (both ends anchored) under **unit-cost Levenshtein**. Both
//! choices are deliberate:
//!
//! - *Global*, because the caller trims common flanks before calling in, so
//!   both ends correspond by construction. That trimming is a precondition,
//!   not a convenience — see [`AlignmentDag::build`].
//! - *Unit cost, not affine*. The criterion is not "which alignment is best"
//!   but "which steps appear in **all** minimal alignments" — the cost model's
//!   job is to define that set, not to pick a winner. Affine penalties are
//!   tuned parameters, and tuning is where the previous principled attempts at
//!   this died. Unit cost has no knob to tune.
//!
//! There is **no shifting here**. The DAG holds all minimal alignments at once,
//! so it has no 5'/3' bias to choose; that is exactly why it is confluent.
//! Shifting is the renderer's job.
//!
//! There is also **no seeding from the input's own alignment** — that would
//! reintroduce input-relativity by a new route.
//!
//! # The DP is banded, and the band is exact (#1928)
//!
//! A cell `(i, j)` sits on diagonal `d = i - j`. Reaching it from the source
//! costs at least `|d|`, because every step off the main diagonal is an indel;
//! reaching the sink from it costs at least `|delta - d|`, where
//! `delta = ref_len - alt_len`. So every cell on an alignment of cost `k`
//! satisfies `|d| + |delta - d| <= k`.
//!
//! `on_optimal_path` **is** the set of cells lying on some minimal alignment, so
//! every cell this module retains already satisfies that bound for `k` = the
//! edit distance. Computing only the band is therefore a **cost** change and not
//! a semantic one: the retained cells, the out-edge masks and
//! [`AlignmentDag::dominators`] come out bit-identical to what the full
//! `(n + 1) x (m + 1)` sweep produces. That is not left to the argument: the
//! five `the_band_agrees_with_the_full_grid_*` tests keep the pre-band
//! implementation verbatim and diff the two builds field by field.
//!
//! **They are not redundant with the oracles already here.** Deliberately
//! widening the band's acceptance test by one — so a band one doubling too
//! narrow is accepted — is caught *only* by those tests:
//! `dominators_agree_with_brute_force_on_exhaustive_small_inputs` passes, because
//! the brute-force probe is computed from the DAG itself and a systematically
//! narrowed DAG is internally consistent. An oracle that reads the artefact
//! cannot see the artefact being too small.
//!
//! The band is found by doubling `k` until the band-confined distance no longer
//! exceeds its own band, which costs at most twice the final pass. The **DP** is
//! then `Θ((n + m) · k)` rather than `Θ(n · m)`.
//!
//! # `build` as a whole is not, and the measurement says so (#1928)
//!
//! **Storage is unchanged** — the grids keep their full `(n + 1) x (m + 1)`
//! shape and their `i * width + j` addressing. So `build` still *allocates and
//! zeroes* `Θ(n · m)` even though it now *fills* `Θ((n + m) · k)`, and past a
//! certain width that allocation is the whole cost.
//!
//! Measured on a quiet dedicated c7i.4xlarge (x86_64, 4 kB pages), criterion
//! against the immediately preceding commit, so banding is the only variable.
//! **Ratios only** — absolute wall-clock is not comparable across hosts, and the
//! raw figures are in `campaigns/2026-08-14-1928-banded-dp-benchmark-protocol.md`
//! with the run that produced them:
//!
//! | n (k = 2) | speedup | cell-count ratio |
//! |---:|---:|---:|
//! | 256 | 37.0x | 74x |
//! | 1024 | 6.6x | 293x |
//! | 4096 | 7.9x | 1171x |
//!
//! The cell-count column is **derived, not measured**, and is reproducible from
//! this file: at `n = 4096, k = 2` the full grids take `2 * 4097^2 = 33,570,818`
//! cell-writes against the band's 28,675 — counting *both* the prefix and suffix
//! grids and the doubling ladder's discarded `k = 1` pass. So geometry
//! **over-predicts the measured speedup by ~148x**, and the gap widens with `n`
//! (only ~2x at `n = 256`). At `n = 4096` `build` allocates ~101 MB of grids, and
//! 101 MB / 20 ms is ~5 GB/s — i.e. the banded time is essentially all
//! allocation and none of it DP.
//!
//! Two things follow, and neither is obvious from the code:
//!
//! - **Do not quote a single speedup.** `by_length_fixed_k/len1024_k2` and
//!   `by_divergence/k2` are the *same* computation and measured **5.8x apart**,
//!   because the latter runs after the `n = 4096` cases have already faulted in a
//!   large arena. Banded timings are sensitive to allocator state in a way the
//!   unbanded ones were not (their baselines differ by only 15%).
//! - **Narrowing the storage to the band is now the lever**, not a tidy-up. It is
//!   a separate change with a real refactor behind it; see [`AlignmentDag::build`].
//!
//! Whether that unfilled allocation costs resident *memory* is **not measured
//! here** — criterion has no memory axis, and the peak-RSS harness was not run
//! against this commit. Do not read the time win as a memory win.
//!
//! At high divergence the band is the whole grid and the doubling search is
//! overhead paid for nothing: at `n = 1024`, `k = 614` (~0.6n, the uniform-random
//! regime) this is **11.2% slower**, having still been 1.56x faster at `k = 512`.
//! Real variants sit at single-digit `k`, which is the regime the corpus presents.

/// One step of an alignment, as an edge in the DAG.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Step {
    /// Reference and alternate bases agree; consumes one of each.
    Match,
    /// Reference and alternate bases differ; consumes one of each.
    Sub,
    /// Consumes a reference base with no alternate base.
    Del,
    /// Consumes an alternate base with no reference base.
    Ins,
}

/// The DAG of every minimal-cost alignment of a reference block to an alternate
/// block.
///
/// Nodes are grid cells `(i, j)`, meaning "`i` reference bases and `j`
/// alternate bases consumed". A cell is present only if some minimal alignment
/// passes through it. Edges are the steps a minimal alignment may take.
///
/// `pub` (rather than `pub(crate)`) so the `seqfirst_align` criterion
/// benchmark can name the type; actually reachable from outside the crate
/// only when the enclosing `seqfirst`/`align` modules are also `pub`, which
/// they are only under the `dev` feature — see `seqfirst/mod.rs`.
pub struct AlignmentDag {
    ref_len: u32,
    alt_len: u32,
    /// Read only through [`AlignmentDag::edit_distance`], whose sole non-test
    /// caller is the `dev`-gated `seqfirst_align` benchmark — so without that
    /// feature the field is genuinely unread. Scoped on the feature rather than
    /// allowed unconditionally, so it starts warning again if the benchmark is
    /// ever the only thing keeping it alive under `dev` too.
    #[cfg_attr(not(feature = "dev"), allow(dead_code))]
    total: u32,
    /// Lower bound (inclusive) on `i - j` for the diagonals [`AlignmentDag::build`]
    /// examined. Every on-path cell lies inside `band_lo ..= band_hi`; a cell
    /// outside was never written and so reads `false` from `on_optimal_path`,
    /// which is what makes [`AlignmentDag::cells`] safe to clip to the same
    /// range.
    band_lo: i64,
    /// Upper bound (inclusive) on `i - j`; see [`AlignmentDag::band_lo`].
    band_hi: i64,
    /// Row-major `(ref_len + 1) x (alt_len + 1)`; `true` when the cell lies on
    /// at least one minimal alignment.
    on_optimal_path: Vec<bool>,
    /// Row-major, parallel to `on_optimal_path`. One bit per possible out-edge:
    /// [`OUT_MATCH`], [`OUT_SUB`], [`OUT_DEL`], [`OUT_INS`].
    ///
    /// A flat bitmask rather than a `Vec<Vec<_>>`. A cell has at most three
    /// out-edges (diagonal, down, right) and the diagonal's kind is fixed by
    /// whether the two bases are equal, so four bits describe a cell exactly.
    /// The nested form cost 24 bytes of `Vec` header per grid cell whether or
    /// not the cell was on-path — about 400 MB before a single edge was stored
    /// for the 4096 x 4096 block `benches/seqfirst_align.rs` measures — plus one
    /// allocator call per on-path cell. This keeps the same `Θ(n·m)` space with
    /// one allocation and an order-of-magnitude smaller constant.
    out: Vec<u8>,
}

/// Diagonal out-edge whose two bases are equal (`Step::Match`).
const OUT_MATCH: u8 = 1 << 0;
/// Diagonal out-edge whose two bases differ (`Step::Sub`).
const OUT_SUB: u8 = 1 << 1;
/// Down out-edge: a reference base is deleted (`Step::Del`).
const OUT_DEL: u8 = 1 << 2;
/// Right out-edge: an alternate base is inserted (`Step::Ins`).
const OUT_INS: u8 = 1 << 3;

/// Largest `n + m` (reference length plus alternate length) whose alignment
/// [`AlignmentDag::build`] will compute.
///
/// # Why `n + m`, and why it is exact rather than conservative
///
/// The cost grids are `u16`. A prefix cost `prefix[i][j]` is the edit distance
/// of `ref[..i]` to `alt[..j]`, so it is bounded by `max(i, j) <= i + j`; a
/// suffix cost `suffix(i, j)` is likewise bounded by `(n - i) + (m - j)`. Every
/// value the DP writes is one of those, and the largest sum `build` ever forms
/// is `on_optimal_path`'s `prefix[i][j] + suffix(i, j) <= n + m`. The three
/// `here + 1` / `here + cost` comparisons are formed only where `i < n` or
/// `j < m`, so `here <= n + m - 1` there. `n + m <= u16::MAX` therefore covers
/// every arithmetic site in `build` with nothing left over — the bound is the
/// exact one, not a margin.
///
/// # Why it is a refusal and not a `debug_assert`
///
/// `[profile.release]` sets neither `debug-assertions` nor `overflow-checks`,
/// so in the shipped library and in the Python wheel a `debug_assert` is a
/// comment: past the bound the grid fill's boundary row would truncate `j as
/// u16` and the sums would wrap, and `on_optimal_path` would be computed from
/// wrong costs. That is a silently wrong partition and so a silently wrong
/// emitted description — strictly worse than declining. Every gate in front of
/// `build`
/// bounds **cells**, i.e. `(n + 1) * (m + 1)`, and a degenerate `1 x 70,000`
/// grid clears the widest of those by four orders of magnitude while `n + m` is
/// past `u16::MAX`; a long literal `ins` payload reaches exactly that shape.
/// See #1970.
pub(crate) const MAX_ALIGNMENT_SPAN: usize = u16::MAX as usize;
/// Cost of a cell no path *confined to the band* reaches.
///
/// It never survives to be read as a real cost: every cell the band admits is
/// reachable by a staircase that stays inside the band, so this exists only to
/// keep the `min` total at the band's edges, where one or two of the three
/// predecessors lie outside it.
const UNREACHED: u16 = u16::MAX;

/// `ceil(x / 2)`, for a possibly negative `x`.
fn ceil_half(x: i64) -> i64 {
    (x + 1).div_euclid(2)
}

/// `floor(x / 2)`, for a possibly negative `x`.
fn floor_half(x: i64) -> i64 {
    x.div_euclid(2)
}

/// The diagonals an alignment of cost at most `k` can occupy, as inclusive
/// bounds on `i - j`.
///
/// From `|d| + |delta - d| <= k` (see the module header), with
/// `delta = ref_len - alt_len`:
///
/// ```text
/// ceil((delta - k) / 2)  <=  d  <=  floor((delta + k) / 2)
/// ```
///
/// `lo + hi == delta` exactly, so the interval is symmetric about `delta / 2`.
/// That is not decoration: reversing both blocks maps diagonal `d` to
/// `delta - d`, so the *same* interval bands the suffix grid and no second
/// derivation is needed.
///
/// `k` must be at least `|delta|` — no alignment is cheaper than the length
/// difference, and below that the interval excludes the source or the sink.
fn band_for(ref_len: usize, alt_len: usize, k: usize) -> (i64, i64) {
    let delta = ref_len as i64 - alt_len as i64;
    let k = k as i64;
    debug_assert!(
        k >= delta.abs(),
        "band k = {k} is below |ref_len - alt_len| = {}",
        delta.abs()
    );
    let (lo, hi) = (ceil_half(delta - k), floor_half(delta + k));
    debug_assert_eq!(lo + hi, delta, "the band must be symmetric about delta/2");
    (lo, hi)
}

/// The columns row `i` occupies inside the band `lo ..= hi`, clamped to the
/// grid. Never empty for a `k >= |delta|` band.
fn row_span(i: usize, lo: i64, hi: i64, alt_len: usize) -> (usize, usize) {
    debug_assert!(lo <= 0 && hi >= 0, "the band must contain diagonal 0");
    let row = i as i64;
    let first = (row - hi).max(0) as usize;
    let last = (row - lo).min(alt_len as i64) as usize;
    debug_assert!(first <= last, "row {i} has an empty band span");
    (first, last)
}

/// Levenshtein cost grid restricted to diagonals `lo ..= hi`, written into
/// `grid` (row-major, `alt.len() + 1` columns — the *full* grid's shape).
///
/// Cells outside the band are **left untouched**, which is the whole point and
/// is why nothing may read them. Each predecessor is guarded by an explicit
/// in-band test rather than by pre-filling the grid with a sentinel: that
/// pre-fill would put back the `Θ(n·m)` pass the band exists to remove.
///
/// Safe to call repeatedly on one buffer with a widening band, which is what the
/// doubling search in [`AlignmentDag::build`] does — a wider band is a superset,
/// and every cell of it is recomputed here in dependency order.
fn fill_edit_grid_banded(reference: &[u8], alt: &[u8], lo: i64, hi: i64, grid: &mut [u16]) {
    let n = reference.len();
    let m = alt.len();
    let width = m + 1;
    for i in 0..=n {
        let (first, last) = row_span(i, lo, hi, m);
        let row = i as i64;
        for j in first..=last {
            if i == 0 && j == 0 {
                grid[0] = 0;
                continue;
            }
            let col = j as i64;
            let mut best = UNREACHED;
            // Diagonal: same `d`, so in band whenever this cell is.
            if i > 0 && j > 0 {
                let cost = u16::from(reference[i - 1] != alt[j - 1]);
                best = best.min(grid[(i - 1) * width + (j - 1)].saturating_add(cost));
            }
            // Up (a reference base deleted). Its diagonal is `(i - 1) - j`,
            // one below this cell's, so it is in band iff `row - col > lo`.
            if i > 0 && row - col > lo {
                best = best.min(grid[(i - 1) * width + j].saturating_add(1));
            }
            // Left (an alternate base inserted). Its diagonal is `i - (j - 1)`,
            // one above this cell's, so it is in band iff `row - col < hi`.
            if j > 0 && row - col < hi {
                best = best.min(grid[i * width + (j - 1)].saturating_add(1));
            }
            grid[i * width + j] = best;
        }
    }
}

impl AlignmentDag {
    /// Build the DAG for `ref_block` -> `alt_block`.
    ///
    /// Both slices are raw bases; either may be empty.
    ///
    /// **The caller must trim common flanks first.** The partition is a
    /// function of this exact pair, and an untrimmed flank base can be absorbed
    /// into a neighbouring member rather than sitting outside it: `AT -> AAT`
    /// yields one member at ref `0..1`, while the flank-trimmed `T -> AT`
    /// yields a pure insertion at ref `1`. Trimming is therefore a
    /// precondition, not a convenience.
    ///
    /// `None` when `ref_block.len() + alt_block.len()` exceeds
    /// [`MAX_ALIGNMENT_SPAN`], which is the width of the `u16` cost grid. Every
    /// caller in the normalizer already has a decline path, and declining is the
    /// only safe answer: past the bound the costs wrap rather than panic in the
    /// shipped profile. See [`MAX_ALIGNMENT_SPAN`].
    ///
    /// The refusal is checked **before** the band is derived, so it bounds the
    /// whole computation and not merely the fill: `k` itself is bounded by
    /// `n + m`, and every value the banded DP writes is a prefix or suffix cost
    /// under the same bound as the full sweep. Banding narrows which cells are
    /// visited, never how large a cost can get.
    ///
    /// # Cost
    ///
    /// The **DP** is `Θ((n + m) · k)` for `k` the edit distance — it visits only
    /// the band of diagonals a minimal alignment can occupy (module header), and
    /// the doubling search for `k` costs at most twice its final pass. At
    /// `k ≈ 0.6 n` (a uniform-random pair) the band is essentially the whole
    /// grid and this degenerates to the previous `Θ(n·m)`, plus the doubling's
    /// constant — measured at `n = 1024`, that is 1.56x *faster* still at
    /// `k = 512` and 11.2% *slower* at `k = 614`.
    ///
    /// **`build` as a whole is not `Θ((n + m) · k)`, because space is still
    /// `Θ(n · m)`** — deliberately, at this step: the grids keep their full
    /// `(n + 1) x (m + 1)` shape and their `i * width + j` addressing, so
    /// `index`, `contains_cell`, `out_edges` and `dominators`'s in-degree table
    /// are untouched. So the allocation is unnarrowed while the fill is narrowed,
    /// and past a certain width the allocation is the entire cost: the measured
    /// speedup at `n = 4096, k = 2` is **7.9x against a ~1171x cell-count
    /// ratio**. Narrowing the *storage* to the band is a separate change with a
    /// real refactor behind it, and the module header's table is the argument for
    /// doing it.
    ///
    /// # What this does to resident memory is UNMEASURED — do not infer it
    ///
    /// The untouched part of a zero-initialised `Vec` is never resident, so a
    /// banded write pattern could in principle shrink peak RSS by itself. Whether
    /// it does **depends on the page size**, and the answer differs across the
    /// two platforms this project builds on:
    ///
    /// - **16 kB pages** (Apple silicon): a grid row is 7.7 kB at the canonical
    ///   window, so a row is *under* one page and consecutive rows share pages.
    ///   The band crosses every row, so every page is faulted in regardless and
    ///   banding the computation should save nothing.
    /// - **4 kB pages** (the Linux host the timings above come from): a row spans
    ///   ~2 pages, and the band is diagonal — row `i` touches columns near `i`,
    ///   hence only *one* of that row's pages. So a substantial fraction of each
    ///   grid may never be faulted.
    ///
    /// Both are predictions from geometry. **Neither has been measured**:
    /// criterion has no memory axis, and the peak-RSS harness was not run against
    /// this commit on either page size. Do not cite the time win as a memory win,
    /// and do not carry the 16 kB conclusion onto a 4 kB host.
    pub fn build(ref_block: &[u8], alt_block: &[u8]) -> Option<Self> {
        let n = ref_block.len();
        let m = alt_block.len();
        if n.checked_add(m)? > MAX_ALIGNMENT_SPAN {
            return None;
        }
        let width = m + 1;
        let cell_count = (n + 1) * width;
        let delta = (n as i64 - m as i64).unsigned_abs() as usize;
        // The band at which every diagonal is in play. Reaching it terminates
        // the loop on an exact answer whatever the doubling does, so the
        // convergence test below is an optimisation and not the only exit.
        let widest = (n + m).max(1);

        let mut prefix = vec![0u16; cell_count];
        // No alignment costs less than the length difference, so that is the
        // narrowest band worth trying; `max(1)` keeps a diagonal of width one
        // for the equal-length case.
        let mut k = delta.max(1).min(widest);
        let (mut lo, mut hi) = band_for(n, m, k);
        loop {
            fill_edit_grid_banded(ref_block, alt_block, lo, hi, &mut prefix);
            // Confining the DP can only overestimate, so a distance the band
            // itself admits is the true one: every cell of every minimal
            // alignment then lies inside `lo ..= hi`, and the band-confined
            // value at each of them is its true value.
            if prefix[n * width + m] as usize <= k || k >= widest {
                break;
            }
            k = (2 * k).min(widest);
            (lo, hi) = band_for(n, m, k);
        }
        let total: u16 = prefix[n * width + m];

        // Suffix distances, computed as prefix distances of the reversed
        // blocks: `suffix(i, j)` is the cost of aligning `ref[i..]` to
        // `alt[j..]`. Reversal maps diagonal `d` to `delta - d`, and the band is
        // symmetric about `delta / 2`, so it carries over unchanged.
        let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
        let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
        let mut suffix_grid = vec![0u16; cell_count];
        fill_edit_grid_banded(&ref_rev, &alt_rev, lo, hi, &mut suffix_grid);
        let suffix = |i: usize, j: usize| suffix_grid[(n - i) * width + (m - j)];

        let mut on_optimal_path = vec![false; cell_count];
        for i in 0..=n {
            let (first, last) = row_span(i, lo, hi, m);
            for j in first..=last {
                on_optimal_path[i * width + j] =
                    prefix[i * width + j].saturating_add(suffix(i, j)) == total;
            }
        }

        let mut out = vec![0u8; cell_count];
        for i in 0..=n {
            let (first, last) = row_span(i, lo, hi, m);
            for j in first..=last {
                if !on_optimal_path[i * width + j] {
                    continue;
                }
                let here = prefix[i * width + j];
                // `on_optimal_path` is tested *first* in each pair below. A
                // neighbour outside the band holds an unwritten `prefix` cell,
                // and its `false` flag is what keeps that cell from being read.
                // Diagonal: match or substitution.
                if i < n && j < m {
                    let cost = u16::from(ref_block[i] != alt_block[j]);
                    let next = (i + 1) * width + (j + 1);
                    if on_optimal_path[next] && here + cost == prefix[next] {
                        out[i * width + j] |= if cost == 0 { OUT_MATCH } else { OUT_SUB };
                    }
                }
                // Down: deletion of a reference base.
                if i < n {
                    let next = (i + 1) * width + j;
                    if on_optimal_path[next] && here + 1 == prefix[next] {
                        out[i * width + j] |= OUT_DEL;
                    }
                }
                // Right: insertion of an alternate base.
                if j < m {
                    let next = i * width + (j + 1);
                    if on_optimal_path[next] && here + 1 == prefix[next] {
                        out[i * width + j] |= OUT_INS;
                    }
                }
            }
        }

        Some(Self {
            ref_len: n as u32,
            alt_len: m as u32,
            total: u32::from(total),
            band_lo: lo,
            band_hi: hi,
            on_optimal_path,
            out,
        })
    }

    /// Minimal edit distance between the two blocks.
    ///
    /// Called by the `dev`-gated `seqfirst_align` benchmark and by this module's
    /// tests; the shadow comparison in `normalize_allele` does not need it, hence
    /// the feature-scoped allowance.
    #[cfg_attr(not(feature = "dev"), allow(dead_code))]
    pub fn edit_distance(&self) -> u32 {
        self.total
    }

    /// Length of the reference block.
    pub(crate) fn ref_len(&self) -> u32 {
        self.ref_len
    }

    /// Length of the alternate block.
    pub(crate) fn alt_len(&self) -> u32 {
        self.alt_len
    }

    fn index(&self, i: u32, j: u32) -> usize {
        i as usize * (self.alt_len as usize + 1) + j as usize
    }

    /// Whether the cell lies on at least one minimal alignment.
    pub(crate) fn contains_cell(&self, i: u32, j: u32) -> bool {
        i <= self.ref_len && j <= self.alt_len && self.on_optimal_path[self.index(i, j)]
    }

    /// Every cell on a minimal alignment, in topological order (`i + j`
    /// ascending, then `i` ascending). Every edge strictly increases `i + j`,
    /// so this ordering is a valid topological sort.
    ///
    /// Emitted by walking anti-diagonals (`d = i + j`) directly, rather than
    /// collecting every on-path cell and sorting it — diagonal order is
    /// already topological, so no sort is needed.
    ///
    /// Each anti-diagonal is further clipped to the band `build` examined: with
    /// `j = d - i`, the constraint `band_lo <= i - j <= band_hi` is
    /// `2i - d` in that range, i.e. `i` in
    /// `ceil((d + band_lo) / 2) ..= floor((d + band_hi) / 2)`. Nothing is lost —
    /// a cell outside the band was never written and so is `false` in
    /// `on_optimal_path` — the clip removes scanning, not cells. That is what
    /// takes this walk from `Θ(n·m)` to `Θ((n + m) · k)`.
    pub(crate) fn cells(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        let (lo, hi) = (self.band_lo, self.band_hi);
        (0..=self.ref_len + self.alt_len)
            .flat_map(move |d| {
                let diag = i64::from(d);
                let i_min = i64::from(d.saturating_sub(self.alt_len)).max(ceil_half(diag + lo));
                let i_max = i64::from(d.min(self.ref_len)).min(floor_half(diag + hi));
                (i_min..=i_max).map(move |i| (i as u32, (diag - i) as u32))
            })
            .filter(move |&(i, j)| self.contains_cell(i, j))
    }

    /// Out-edges of a cell, as `(next_i, next_j, step)`.
    ///
    /// `(i, j)` must be within the grid (`i <= ref_len`, `j <= alt_len`).
    /// Off-grid is not uniformly caught in release: with `i < ref_len` an
    /// out-of-range `j` wraps into a later row and silently returns another
    /// cell's edges; in the last row (`i == ref_len`), or for any
    /// `i > ref_len`, the flat index runs past the backing `Vec` and panics.
    /// The `debug_assert` below catches both in debug builds.
    pub(crate) fn out_edges(&self, i: u32, j: u32) -> impl Iterator<Item = (u32, u32, Step)> + '_ {
        debug_assert!(
            i <= self.ref_len && j <= self.alt_len,
            "out_edges({i}, {j}) is off-grid ({} x {}) and would alias another cell",
            self.ref_len,
            self.alt_len
        );
        // Reconstructed from the mask and `(i, j)`. The order matches what the
        // nested-`Vec` form pushed — diagonal, then down, then right — so any
        // caller that depended on edge order is unaffected.
        let mask = self.out[self.index(i, j)];
        [
            (mask & OUT_MATCH != 0).then_some((i + 1, j + 1, Step::Match)),
            (mask & OUT_SUB != 0).then_some((i + 1, j + 1, Step::Sub)),
            (mask & OUT_DEL != 0).then_some((i + 1, j, Step::Del)),
            (mask & OUT_INS != 0).then_some((i, j + 1, Step::Ins)),
        ]
        .into_iter()
        .flatten()
    }
}

/// Steps that every minimal alignment takes.
///
/// `matched` holds `(reference offset, alternate offset)` for each match edge
/// present in every minimal alignment — the positions that are genuinely
/// unchanged, together with where they sit in the alternate block. Sorted by
/// reference offset. `forced_ins_junctions` holds junctions (reference offsets,
/// with junction `k` sitting between reference bases `k - 1` and `k`) where
/// every minimal alignment inserts.
///
/// The alternate offset is not decoration: it is the only reliable way to
/// recover a member's alternate span. See `partition::member_alt_spans`.
///
/// # Why node agreement is not enough
///
/// A reference position can lie on every minimal alignment's path (node
/// agreement) without any single match *edge* into it being common to all of
/// them — and only a common edge pins a single alternate coordinate.
///
/// Worked example: `AT -> AAT`. Reference offset 0 (`A`) is `Match`-stepped on
/// both minimal alignments, but at different alternate offsets: one path
/// matches it to alternate offset 0, the other to alternate offset 1 (the
/// inserted `A` sits before it on one path and after it on the other). No
/// single match edge is common to both, so reference offset 0 is *not* a
/// dominator match — it counts as changed. `matched = [(1, 2)]`, `changed =
/// [0]`.
///
/// (`AAAAAAA -> AACACAAAA` is not a usable example here: all seven of its
/// reference positions are dominator-matched, so its `changed` set is empty.)
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct Dominators {
    pub(crate) matched: Vec<(u32, u32)>,
    pub(crate) forced_ins_junctions: Vec<u32>,
}

impl Dominators {
    /// Reference offsets that every minimal alignment matches.
    pub(crate) fn matched_ref(&self) -> Vec<u32> {
        self.matched.iter().map(|&(r, _)| r).collect()
    }

    /// The alternate offset pinned by the dominator match at `ref_pos`, if that
    /// reference position is a dominator match.
    pub(crate) fn alt_at(&self, ref_pos: u32) -> Option<u32> {
        self.matched
            .binary_search_by_key(&ref_pos, |&(r, _)| r)
            .ok()
            .map(|idx| self.matched[idx].1)
    }
}

impl AlignmentDag {
    /// Edges present in **every** minimal alignment.
    ///
    /// An edge lies on every source-to-sink path exactly when removing it
    /// disconnects them. In a DAG this is found with one topological sweep:
    /// track how many edges cross from the processed set to the unprocessed
    /// set, and whenever that count falls to one, that single edge is on every
    /// path.
    ///
    /// The sweep is `O(V + E)`, and since #1928 `cells()` walks only the band —
    /// `Θ((n + m) · k)` — rather than filtering the whole grid. What is still
    /// `Θ(n·m)` here is the `in_degree` table, which is sized to the full grid
    /// because `index` addresses the full grid; it is a zeroed allocation of
    /// which only band cells are ever touched. Deliberately not the per-edge
    /// reachability probe the calibration prototype used — that is `O(E)` per
    /// edge, which is fine for a 52 nt block and not for production. The probe
    /// survives as the test oracle.
    pub(crate) fn dominators(&self) -> Dominators {
        let mut result = Dominators::default();
        let sink = (self.ref_len, self.alt_len);
        let order: Vec<(u32, u32)> = self.cells().collect();

        // In-degree within the DAG, needed to know when a node's incoming edges
        // stop crossing the frontier.
        let mut in_degree = vec![0u32; (self.ref_len as usize + 1) * (self.alt_len as usize + 1)];
        for &(i, j) in &order {
            for (ni, nj, _) in self.out_edges(i, j) {
                in_degree[self.index(ni, nj)] += 1;
            }
        }

        let mut crossing: i64 = 0;
        for &(i, j) in &order {
            crossing -= i64::from(in_degree[self.index(i, j)]);
            // Counted, not collected. The edge itself is only needed on the
            // `crossing == 1` branch below, and that is rare; collecting on every
            // cell put a heap allocation back on the sweep that the flat
            // adjacency exists to keep off it.
            let out_degree = self.out_edges(i, j).count();
            crossing += out_degree as i64;

            // If exactly one edge crosses the frontier now, it is on every
            // path. Any node with two or more out-edges contributes two or more
            // crossings itself, so a count of one means this node has exactly
            // one out-edge and nothing older is still crossing.
            if crossing == 1 && (i, j) != sink {
                debug_assert_eq!(
                    out_degree, 1,
                    "a single crossing must be this cell's only out-edge"
                );
                let Some((_, _, step)) = self.out_edges(i, j).next() else {
                    continue;
                };
                match step {
                    // `(i, j)` is the cell the edge leaves, so it pins the
                    // alternate coordinate as well as the reference one.
                    Step::Match => result.matched.push((i, j)),
                    Step::Ins => result.forced_ins_junctions.push(i),
                    Step::Sub | Step::Del => {}
                }
            }
        }

        result.matched.sort_unstable();
        result.matched.dedup();
        result.forced_ins_junctions.sort_unstable();
        result.forced_ins_junctions.dedup();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// [`AlignmentDag::build`] on blocks small enough that the refusal cannot
    /// fire. Every block in this module is a handful of bases, so an `Option`
    /// here would only be noise — but the refusal is exercised on purpose by
    /// [`tests::the_span_bound_is_enforced_from_both_sides`].
    fn dag(reference: &[u8], alt: &[u8]) -> AlignmentDag {
        AlignmentDag::build(reference, alt).expect("test blocks are far below MAX_ALIGNMENT_SPAN")
    }

    #[test]
    fn edit_distance_matches_known_values() {
        for (r, a, expected) in [
            (&b""[..], &b""[..], 0),
            (&b"ACGT"[..], &b"ACGT"[..], 0),
            (&b"ACGT"[..], &b"AGGT"[..], 1),
            (&b"ACGT"[..], &b"ACT"[..], 1),
            (&b"ACGT"[..], &b"ACGGT"[..], 1),
            (&b"AAAAAAA"[..], &b"AACACAAAA"[..], 2),
            (&b"AAAAAAA"[..], &b"ACAAAA"[..], 2),
        ] {
            let dag = dag(r, a);
            assert_eq!(
                dag.edit_distance(),
                expected,
                "d({}, {})",
                String::from_utf8_lossy(r),
                String::from_utf8_lossy(a)
            );
        }
    }

    #[test]
    fn every_cell_lies_on_a_minimal_path() {
        // A cell is in the DAG only if some minimal alignment passes through
        // it: prefix(i,j) + suffix(i,j) == total.
        let dag = dag(b"ACGTACGT", b"ATGTTCGT");
        assert!(dag.contains_cell(0, 0), "source must be present");
        assert!(
            dag.contains_cell(dag.ref_len(), dag.alt_len()),
            "sink must be present"
        );
        for (i, j) in dag.cells() {
            assert!(
                dag.out_edges(i, j).next().is_some() || (i, j) == (dag.ref_len(), dag.alt_len()),
                "cell ({i},{j}) is a dead end but is not the sink"
            );
        }
    }

    #[test]
    fn every_edge_advances_the_topological_key() {
        // i + j strictly increases along every edge, which is what makes
        // sorting cells by i + j a valid topological order.
        let dag = dag(b"ACGTACGTAC", b"ACGTTACGTGAC");
        for (i, j) in dag.cells() {
            for (ni, nj, _step) in dag.out_edges(i, j) {
                assert!(
                    ni + nj > i + j,
                    "edge ({i},{j}) -> ({ni},{nj}) does not advance i+j"
                );
            }
        }
    }

    /// Every word of length `len` over a two-letter alphabet.
    fn words(len: u32) -> Vec<Vec<u8>> {
        if len == 0 {
            return vec![Vec::new()];
        }
        let mut out = Vec::new();
        for shorter in words(len - 1) {
            for base in *b"AC" {
                let mut w = shorter.clone();
                w.push(base);
                out.push(w);
            }
        }
        out
    }

    /// Brute-force oracle: an edge is a dominator iff deleting it disconnects
    /// the source from the sink. `O(E)` per edge — correct but far too slow for
    /// production, which is exactly why it belongs in the tests.
    fn dominators_by_brute_force(dag: &AlignmentDag) -> (Vec<(u32, u32)>, Vec<u32>) {
        let sink = (dag.ref_len(), dag.alt_len());
        let mut matched = Vec::new();
        let mut forced_ins = Vec::new();
        for (i, j) in dag.cells().collect::<Vec<_>>() {
            for (ni, nj, step) in dag.out_edges(i, j).collect::<Vec<_>>() {
                if matches!(step, Step::Sub | Step::Del) {
                    continue;
                }
                let banned = (i, j, ni, nj);
                let mut seen = vec![(0u32, 0u32)];
                let mut stack = vec![(0u32, 0u32)];
                let mut reached = false;
                while let Some(v) = stack.pop() {
                    if v == sink {
                        reached = true;
                        break;
                    }
                    for (wi, wj, _) in dag.out_edges(v.0, v.1).collect::<Vec<_>>() {
                        if (v.0, v.1, wi, wj) == banned || seen.contains(&(wi, wj)) {
                            continue;
                        }
                        seen.push((wi, wj));
                        stack.push((wi, wj));
                    }
                }
                if !reached {
                    match step {
                        Step::Match => matched.push((i, j)),
                        Step::Ins => forced_ins.push(i),
                        _ => unreachable!(),
                    }
                }
            }
        }
        matched.sort_unstable();
        matched.dedup();
        forced_ins.sort_unstable();
        forced_ins.dedup();
        (matched, forced_ins)
    }

    #[test]
    fn lrg199_has_exactly_one_dominator_match() {
        // The spec's own worked example, `delins.md:44` / LRG_199t1:c.850_901.
        // Real sequence, pulled from LRG_199.fasta; CDS offset verified against
        // the spec's own `c.145_147 = CGC = Arg49`.
        let dag = dag(
            b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT",
            b"TTCCTCGATGCCTG",
        );
        let dom = dag.dominators();
        assert_eq!(
            dom.matched_ref(),
            vec![48],
            "LRG_199 must have exactly one dominator match, at reference offset 48"
        );
        assert_eq!(
            dom.alt_at(48),
            Some(13),
            "the dominator must pin alternate offset 13, not merely some coordinate"
        );
    }

    #[test]
    fn a_pure_insertion_is_seen_as_a_forced_junction() {
        // #1260: two `insC` at adjacent gaps in a poly-A run. A node model
        // reports this as zero members because no reference position changes;
        // the junctions are the whole signal.
        let dag = dag(b"AAAAAAA", b"AACACAAAA");
        let dom = dag.dominators();
        assert_eq!(
            dom.matched_ref().len(),
            7,
            "every reference position is matched in every minimal alignment"
        );
        assert_eq!(
            dom.forced_ins_junctions,
            vec![2, 3],
            "the insertions must show up as forced junctions at exactly 2 and 3, or they vanish"
        );
    }

    #[test]
    fn dominators_agree_with_brute_force_on_exhaustive_small_inputs() {
        // Exhaustive over a 2-letter alphabet, lengths 0..=5 both sides. Small
        // alphabet and short blocks maximise alternative minimal alignments,
        // which is where a wrong dominator rule shows up.
        let all: Vec<Vec<u8>> = (0..=5u32).flat_map(words).collect();
        let mut checked = 0usize;
        for r in &all {
            for a in &all {
                let dag = dag(r, a);
                let got = dag.dominators();
                let (want_matched, want_ins) = dominators_by_brute_force(&dag);
                assert_eq!(
                    (got.matched, got.forced_ins_junctions),
                    (want_matched, want_ins),
                    "disagreement on {} -> {}",
                    String::from_utf8_lossy(r),
                    String::from_utf8_lossy(a)
                );
                checked += 1;
            }
        }
        assert!(
            checked > 3000,
            "expected a large exhaustive sweep, ran {checked}"
        );
    }

    /// #1970. The `u16` narrowing's bound must be a refusal, pinned from **both**
    /// sides so neither a widened nor a tightened bound passes silently.
    ///
    /// The over-limit arm is not a "does it panic" test: it asserts `is_none()`,
    /// which is a value the shipped profile produces too. A `debug_assert` would
    /// satisfy a panic-shaped test and still wrap in release, which is the whole
    /// defect.
    #[test]
    fn the_span_bound_is_enforced_from_both_sides() {
        // Degenerate on purpose: `1 x m` is the shape that clears every cells
        // gate in front of `build` while `n + m` runs away. A long literal `ins`
        // payload arrives here as exactly this.
        let one = vec![b'A'; 1];

        // n + m == MAX_ALIGNMENT_SPAN: the last accepted pair.
        let just_under = vec![b'C'; MAX_ALIGNMENT_SPAN - 1];
        let built = AlignmentDag::build(&one, &just_under)
            .expect("n + m == MAX_ALIGNMENT_SPAN must still be computed");
        assert_eq!(
            built.edit_distance(),
            MAX_ALIGNMENT_SPAN as u32 - 1,
            "one base against {} differing bases is a substitution plus {} insertions",
            MAX_ALIGNMENT_SPAN - 1,
            MAX_ALIGNMENT_SPAN - 2
        );

        // n + m == MAX_ALIGNMENT_SPAN + 1: the first refused pair. One base
        // wider than the pair above, so nothing but the bound separates them.
        let just_over = vec![b'C'; MAX_ALIGNMENT_SPAN];
        assert!(
            AlignmentDag::build(&one, &just_over).is_none(),
            "n + m == {} must be refused, not wrapped",
            MAX_ALIGNMENT_SPAN + 1
        );

        // And far past it, which is the size #1970 reported.
        assert!(
            AlignmentDag::build(&one, &vec![b'C'; 70_000]).is_none(),
            "a 70,000-base alternate must be refused"
        );
    }

    // ------------------------------------------------------------------
    // #1928: the band is a COST change. These guard that claim directly.
    // ------------------------------------------------------------------

    /// The **pre-band** build, kept verbatim as the differential oracle: a full
    /// `(n + 1) x (m + 1)` Levenshtein grid, its reverse, the
    /// `prefix + suffix == total` membership test and the out-edge masks,
    /// exactly as this module computed them before the band was introduced.
    ///
    /// Deliberately a second implementation rather than "the banded code run
    /// with a full band" — running one implementation with two parameter values
    /// tests the parameters, not the algorithm, and every band-specific
    /// off-by-one would cancel.
    fn full_grid_build(ref_block: &[u8], alt_block: &[u8]) -> AlignmentDag {
        fn edit_grid(reference: &[u8], alt: &[u8]) -> Vec<u16> {
            let n = reference.len();
            let m = alt.len();
            let width = m + 1;
            let mut d = vec![0u16; (n + 1) * width];
            for (i, x) in d.iter_mut().step_by(width).enumerate() {
                *x = i as u16;
            }
            for (j, x) in d.iter_mut().take(width).enumerate() {
                *x = j as u16;
            }
            for i in 1..=n {
                for j in 1..=m {
                    let cost = u16::from(reference[i - 1] != alt[j - 1]);
                    d[i * width + j] = (d[(i - 1) * width + j] + 1)
                        .min(d[i * width + j - 1] + 1)
                        .min(d[(i - 1) * width + j - 1] + cost);
                }
            }
            d
        }

        let n = ref_block.len();
        let m = alt_block.len();
        let width = m + 1;
        let prefix = edit_grid(ref_block, alt_block);
        let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
        let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
        let suffix_grid = edit_grid(&ref_rev, &alt_rev);
        let suffix = |i: usize, j: usize| suffix_grid[(n - i) * width + (m - j)];
        let total: u16 = prefix[n * width + m];

        let mut on_optimal_path = vec![false; (n + 1) * width];
        for i in 0..=n {
            for j in 0..=m {
                on_optimal_path[i * width + j] = prefix[i * width + j] + suffix(i, j) == total;
            }
        }
        let mut out = vec![0u8; (n + 1) * width];
        for i in 0..=n {
            for j in 0..=m {
                if !on_optimal_path[i * width + j] {
                    continue;
                }
                let here = prefix[i * width + j];
                if i < n && j < m {
                    let cost = u16::from(ref_block[i] != alt_block[j]);
                    let next = (i + 1) * width + (j + 1);
                    if here + cost == prefix[next] && on_optimal_path[next] {
                        out[i * width + j] |= if cost == 0 { OUT_MATCH } else { OUT_SUB };
                    }
                }
                if i < n {
                    let next = (i + 1) * width + j;
                    if here + 1 == prefix[next] && on_optimal_path[next] {
                        out[i * width + j] |= OUT_DEL;
                    }
                }
                if j < m {
                    let next = i * width + (j + 1);
                    if here + 1 == prefix[next] && on_optimal_path[next] {
                        out[i * width + j] |= OUT_INS;
                    }
                }
            }
        }
        AlignmentDag {
            ref_len: n as u32,
            alt_len: m as u32,
            total: u32::from(total),
            // The whole grid, expressed as a band, so `cells()` scans all of it
            // and the two walks are compared like for like.
            band_lo: -(m as i64),
            band_hi: n as i64,
            on_optimal_path,
            out,
        }
    }

    /// Assert the banded build and the pre-band build agree on everything the
    /// DAG exposes: the distance, the retained cell set, every out-edge mask,
    /// the topological walk, and `dominators()`. Returns the edit distance so a
    /// caller can assert which `k` it actually exercised.
    fn assert_band_matches_full_grid(r: &[u8], a: &[u8]) -> u32 {
        // `build` refuses past `MAX_ALIGNMENT_SPAN` (#1970); every fixture here is
        // far below it, so a `None` is a bug in the fixture, not a decline to test.
        let banded = AlignmentDag::build(r, a)
            .expect("differential fixtures are far below MAX_ALIGNMENT_SPAN");
        let full = full_grid_build(r, a);
        let label = format!(
            "{:?} -> {:?}",
            String::from_utf8_lossy(r),
            String::from_utf8_lossy(a)
        );
        assert_eq!(banded.total, full.total, "edit distance differs on {label}");
        assert_eq!(
            banded.on_optimal_path, full.on_optimal_path,
            "retained cell set differs on {label}"
        );
        assert_eq!(banded.out, full.out, "out-edge masks differ on {label}");
        assert_eq!(
            banded.cells().collect::<Vec<_>>(),
            full.cells().collect::<Vec<_>>(),
            "topological cell walk differs on {label}"
        );
        assert_eq!(
            banded.dominators(),
            full.dominators(),
            "dominators differ on {label}"
        );
        banded.total
    }

    /// Deterministic pseudo-random stream; no RNG dependency, so any failure
    /// reproduces from the seed named at the call site.
    fn next(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 33
    }

    /// `alt`, mutated by `edits` single-base changes drawn from `state` — a mix
    /// of substitutions, deletions and insertions, so the alternate is not
    /// confined to the main diagonal.
    fn apply_random_edits(reference: &[u8], edits: usize, state: &mut u64) -> Vec<u8> {
        let mut alt = reference.to_vec();
        for _ in 0..edits {
            if alt.is_empty() {
                alt.push(b"ACGT"[(next(state) as usize) % 4]);
                continue;
            }
            let at = (next(state) as usize) % alt.len();
            match next(state) % 3 {
                0 => alt[at] = b"ACGT"[(next(state) as usize) % 4],
                1 => {
                    alt.remove(at);
                }
                _ => alt.insert(at, b"ACGT"[(next(state) as usize) % 4]),
            }
        }
        alt
    }

    #[test]
    fn the_band_agrees_with_the_full_grid_on_exhaustive_small_inputs() {
        // The same exhaustive corpus the brute-force dominator oracle uses:
        // every word over a 2-letter alphabet up to length 5, both sides. It
        // covers the empty block on either side, equal and unequal lengths, and
        // `k` all the way to `n`, since a 2-letter alphabet makes short blocks
        // maximally ambiguous.
        let all: Vec<Vec<u8>> = (0..=5u32).flat_map(words).collect();
        let mut checked = 0usize;
        let mut widest_k = 0u32;
        for r in &all {
            for a in &all {
                widest_k = widest_k.max(assert_band_matches_full_grid(r, a));
                checked += 1;
            }
        }
        assert!(
            checked > 3000,
            "expected a large exhaustive sweep, ran {checked}"
        );
        assert_eq!(
            widest_k, 5,
            "the sweep must reach k == n (a block wholly replaced), saw {widest_k}"
        );
    }

    #[test]
    fn the_band_agrees_with_the_full_grid_at_the_shape_boundaries() {
        // The cases the doubling search is most likely to get wrong, named one
        // by one rather than left to a random draw: both blocks empty, either
        // block empty, a pure insertion, a pure deletion, a length ratio far
        // from one (so `|delta|` and not the divergence sets the initial band),
        // and equal-length blocks sharing no base at all.
        for (r, a, what) in [
            (&b""[..], &b""[..], "both empty"),
            (&b""[..], &b"ACGT"[..], "reference empty — a pure insertion"),
            (&b"ACGT"[..], &b""[..], "alternate empty — a pure deletion"),
            (&b"A"[..], &b""[..], "one base deleted"),
            (&b""[..], &b"A"[..], "one base inserted"),
            (&b"ACGTACGT"[..], &b"ACGTACGT"[..], "k = 0"),
            (&b"ACGTACGT"[..], &b"ACGTACGA"[..], "k = 1, a substitution"),
            (&b"ACGTACGT"[..], &b"ACGTACG"[..], "k = 1, a deletion"),
            (&b"ACGTACGT"[..], &b"ACGTACGTA"[..], "k = 1, an insertion"),
            (
                &b"AAAAAAAAAAAAAAAAAAAA"[..],
                &b"AAA"[..],
                "delta = 17 with k = delta",
            ),
            (
                &b"AAA"[..],
                &b"AAAAAAAAAAAAAAAAAAAA"[..],
                "delta = -17 with k = delta",
            ),
            (
                &b"AAAAAAAAAAAAAAAAAAAA"[..],
                &b"CCC"[..],
                "delta = 17 with k > delta",
            ),
            (
                &b"AAAAAAAA"[..],
                &b"CCCCCCCC"[..],
                "equal length, k = n, no shared base",
            ),
            (
                &b"ACGTACGTACGT"[..],
                &b"TGCATGCATGCA"[..],
                "equal length, maximal substitution divergence",
            ),
        ] {
            let k = assert_band_matches_full_grid(r, a);
            // The assertion above is the guard; this only records that the case
            // exercised the `k` its name claims, so a rewritten fixture cannot
            // silently stop covering the boundary.
            let expected = match what {
                "k = 0" => Some(0),
                "k = 1, a substitution" | "k = 1, a deletion" | "k = 1, an insertion" => Some(1),
                "one base deleted" | "one base inserted" => Some(1),
                "both empty" => Some(0),
                _ => None,
            };
            if let Some(expected) = expected {
                assert_eq!(k, expected, "case {what:?} no longer has k = {expected}");
            }
        }
    }

    #[test]
    fn the_band_agrees_with_the_full_grid_at_every_doubling_boundary() {
        // The doubling schedule is 1, 2, 4, 8, 16, ... so the interesting `k`
        // are the powers of two and their neighbours: one below converges a
        // round early, one above forces the extra widening. Substitutions are
        // spaced so the distance is exactly the number of them, which is
        // asserted rather than assumed — the point of the test is that a
        // *named* `k` was exercised.
        let mut state = 0x1928_u64 | 1;
        let reference: Vec<u8> = (0..256)
            .map(|_| b"ACGT"[(next(&mut state) as usize) % 4])
            .collect();
        for k in [0usize, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33] {
            let mut alt = reference.clone();
            for slot in 0..k {
                // Spaced far apart, so no cheaper indel-bearing alignment
                // exists and the distance really is `k`.
                let at = slot * 7 + 3;
                alt[at] = match alt[at] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
            let seen = assert_band_matches_full_grid(&reference, &alt);
            assert_eq!(
                seen as usize, k,
                "the fixture for k = {k} no longer has that edit distance"
            );
        }
    }

    #[test]
    fn the_band_agrees_with_the_full_grid_across_the_whole_divergence_range() {
        // Randomised, and deliberately spanning the range the band cannot
        // narrow: from `k = 0` (identical blocks, one diagonal) through single
        // digits (a real variant) to a uniform-random alternate of unrelated
        // length, where `k` approaches `n` and the band is the whole grid.
        let mut state = 0xBA4D_u64 | 1;
        let mut checked = 0usize;
        let mut saw_zero = false;
        let mut saw_near_n = false;
        for &n in &[0usize, 1, 2, 5, 9, 16, 33, 64] {
            let reference: Vec<u8> = (0..n)
                .map(|_| b"ACGT"[(next(&mut state) as usize) % 4])
                .collect();
            for trial in 0..24usize {
                let alt: Vec<u8> = if trial >= 20 {
                    // An unrelated alternate of unrelated length: the
                    // `k ≈ 0.6n` regime, plus a `delta` the band must centre on.
                    let len = if n == 0 {
                        0
                    } else {
                        1 + (next(&mut state) as usize) % (2 * n)
                    };
                    (0..len)
                        .map(|_| b"ACGT"[(next(&mut state) as usize) % 4])
                        .collect()
                } else {
                    apply_random_edits(&reference, trial.min(n), &mut state)
                };
                let k = assert_band_matches_full_grid(&reference, &alt);
                saw_zero |= k == 0;
                saw_near_n |= n > 0 && k as usize >= n;
                checked += 1;
            }
        }
        assert_eq!(checked, 192, "expected the full sweep, ran {checked}");
        assert!(saw_zero, "the sweep never exercised k = 0");
        assert!(saw_near_n, "the sweep never exercised k >= n");
    }

    #[test]
    fn the_band_agrees_with_the_full_grid_at_a_realistic_block_size() {
        // The regime the change exists for, and the one small blocks cannot
        // show: hundreds of bases with a single-digit edit distance, where the
        // band is a sliver of the grid rather than most of it. Unequal lengths
        // are included because the indel edits change the block length.
        let mut state = 0xB4E1_u64 | 1;
        let reference: Vec<u8> = (0..512)
            .map(|_| b"ACGT"[(next(&mut state) as usize) % 4])
            .collect();
        for k in [0usize, 1, 2, 4, 8, 16] {
            let alt = apply_random_edits(&reference, k, &mut state);
            assert_band_matches_full_grid(&reference, &alt);
        }
    }

    /// The band must be NARROW, not merely correct.
    ///
    /// Every `the_band_agrees_with_the_full_grid_*` test above asserts that the
    /// banded build **agrees** with the full-grid one — and a full-width band
    /// satisfies all of them trivially. So a change that quietly stopped
    /// narrowing (a doubling seeded at `widest`, a sign slip in `ceil_half`, a
    /// "simplification" of `band_for`) would leave the entire suite green while
    /// deleting the optimization outright. Nothing else in the crate would
    /// notice: the only other signal is a criterion run, and CI does not do one.
    ///
    /// This is the mirror of the hazard the module header records in the other
    /// direction. A band one doubling too *narrow* is caught by the differential
    /// tests; a band too *wide* is caught by nothing, and that is the direction a
    /// well-meaning refactor actually takes.
    ///
    /// Checked against a deliberate sabotage rather than assumed: forcing
    /// `band_for` to return the full-width band (`k = n + m`) leaves the module
    /// at **95 passed, 1 failed** — this test being the one failure, with all
    /// five differential tests and the brute-force dominator oracle still green.
    #[test]
    fn the_band_is_narrow_at_low_divergence_so_a_no_op_regression_is_caught() {
        // `band_for` directly: at a real-variant edit distance the band is a
        // handful of diagonals regardless of how long the blocks are.
        for &n in &[256usize, 1024, 4096] {
            let (lo, hi) = band_for(n, n, 2);
            let width = hi - lo + 1;
            assert_eq!(
                width, 3,
                "band is {width} diagonals at n={n}, k=2; expected a 3-diagonal sliver"
            );
        }
        // Unequal lengths: the band tracks `delta`, so it stays narrow while
        // being centred away from the main diagonal.
        let (lo, hi) = band_for(1024, 1020, 8);
        assert_eq!((lo, hi), (-2, 6), "the band must centre on delta/2 = 2");

        // And the walk the band drives must visit `O((n + m) * k)` cells, not
        // `O(n * m)` — the property the cell-count column in the module header
        // claims, checked rather than asserted in prose.
        let mut state = 0x1928_u64 | 1;
        let reference: Vec<u8> = (0..1024)
            .map(|_| b"ACGT"[(next(&mut state) as usize) % 4])
            .collect();
        let alt = apply_random_edits(&reference, 2, &mut state);
        let dag = AlignmentDag::build(&reference, &alt)
            .expect("1024 + 1024 is far below MAX_ALIGNMENT_SPAN");
        let visited = dag.cells().count();
        let grid = (reference.len() + 1) * (alt.len() + 1);
        assert!(
            visited * 100 < grid,
            "cells() visited {visited} of {grid} cells — that is grid-scale, so \
             the band is not narrowing the walk"
        );
    }
}
