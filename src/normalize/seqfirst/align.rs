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
//! # Storage is banded too (#1988)
//!
//! Every grid is band-major: `(n + 1)` rows of `min(hi - lo, m) + 1` cells
//! (`hi - lo` is `k` up to a parity adjustment — see [`BandLayout`]), addressed
//! by [`BandLayout`], so `build` allocates and zeroes `Θ((n + m) · k)` rather
//! than `Θ(n · m)`. The stride's `min` is what makes that a narrowing in every
//! shape rather than only in the low-divergence one — see [`BandLayout`] for
//! why a diagonal-count stride would be a 500x *pessimisation* on a `1000 x 1`
//! block.
//!
//! This is a **storage** change and not a semantic one, on the same argument as
//! the band itself: a cell outside the band lies on no minimal alignment, so
//! declining to reserve a slot for it removes nothing that was ever `true`. It
//! is not left to the argument either — the five
//! `the_band_agrees_with_the_full_grid_*` tests diff the banded build against
//! the pre-band one over **every cell of the full grid**, so a narrowed cell
//! that should have been retained shows up as a disagreement at that cell.
//!
//! # Why it was worth doing: the shortfall the band left (#1928)
//!
//! Before #1988, storage kept its full `(n + 1) x (m + 1)` shape while the fill
//! was banded, so `build` still allocated and zeroed `Θ(n · m)`, and past a
//! certain width that allocation was the whole cost.
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
//! - **Narrowing the storage to the band was the lever that remained**, and it
//!   is what the section above does (#1988). At `n = 4096, k = 2` the two grids
//!   `build` retains go from `2 * 4097^2 = 33,570,818` cells to
//!   `2 * 4097 * 3 = 24,582`, and the two cost grids with them.
//!
//! **No speedup is claimed here for #1988, and the two figures above must not be
//! read as one.** The 7.9x is #1928's, measured on that host with criterion; no
//! equivalent run has been made for the storage change, and this file's own
//! first bullet is the reason not to guess one from a local timing.
//!
//! What *was* measured for #1988 is **peak RSS, on 16 kB pages only**:
//! **101.8 MiB -> 5.73 MiB** for one `build` at `n = m = 4096, k = 2`. So on
//! that platform the old allocation really was resident in full, and the win is
//! real. The 4 kB-page host is **still unmeasured** and is where the geometry
//! predicts the least to reclaim — see [`AlignmentDag::build`] for both, and do
//! not carry either figure onto the other platform, or read a memory win as a
//! time win.
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
    /// The diagonals [`AlignmentDag::build`] examined, and how their cells are
    /// addressed. Every on-path cell lies inside the band; a cell outside it has
    /// no slot, which is why every read goes through
    /// [`BandLayout::contains`] first and what makes [`AlignmentDag::cells`]
    /// safe to clip to the same range.
    band: BandLayout,
    /// Band-major, `band.len()` cells; `true` when the cell lies on at least one
    /// minimal alignment.
    on_optimal_path: Vec<bool>,
    /// Band-major, parallel to `on_optimal_path`. One bit per possible out-edge:
    /// [`OUT_MATCH`], [`OUT_SUB`], [`OUT_DEL`], [`OUT_INS`].
    ///
    /// A flat bitmask rather than a `Vec<Vec<_>>`. A cell has at most three
    /// out-edges (diagonal, down, right) and the diagonal's kind is fixed by
    /// whether the two bases are equal, so four bits describe a cell exactly.
    /// The nested form cost 24 bytes of `Vec` header per grid cell whether or
    /// not the cell was on-path — about 400 MB before a single edge was stored
    /// for the 4096 x 4096 block `benches/seqfirst_align.rs` measures — plus one
    /// allocator call per on-path cell. This keeps one allocation and an
    /// order-of-magnitude smaller constant, over `Θ((n + m) · k)` cells since
    /// #1988.
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

/// Where a band's cells live in memory (#1988).
///
/// Every grid this module allocates — the two cost grids, `on_optimal_path`,
/// `out`, `dominators`'s in-degree table and `CanonicalAlignment`'s `to_go`
/// planes — is addressed through one of these, so a change here narrows all of
/// them at once.
///
/// **That coupling is structural, not tested, and the distinction matters.** All
/// six read `band.len()` today, and a seventh consumer written the same way
/// inherits the narrowing for free — but nothing *enforces* it. The guard,
/// `the_storage_is_narrow_at_low_divergence_so_a_no_op_regression_is_caught`,
/// pins an exact size on the two grids the DAG **retains** (`on_optimal_path`,
/// `out`) and on the layout itself: **three** of the six. A regression that put
/// `dominators`'s `in_degree` back to `(ref_len + 1) * (alt_len + 1)`, or sized
/// `of`'s `to_go` planes from the grid, passes it untouched — and those are the
/// two largest transients in the pipeline, `to_go` being `2 x size`. They are
/// unguarded because both are locals of methods that drop them before returning,
/// so asserting on them needs a test hook this change does not add. Read the
/// guard as covering the retained grids only; the other three are carried by
/// review.
///
/// # The addressing
///
/// Row `i` occupies columns `first(i) ..= last(i)`, with
/// `first(i) = max(i - hi, 0)` and `last(i) = min(i - lo, m)`. Storing the rows
/// back to back at a fixed `stride` gives
///
/// ```text
/// index(i, j) = i * stride + (j - first(i))
/// ```
///
/// which is the full grid's `i * (m + 1) + j` with the out-of-band prefix of
/// each row squeezed out.
///
/// # Why the stride is `min(hi - lo, m) + 1` and not `hi - lo + 1`
///
/// A band's diagonal count is `hi - lo + 1`, which is `k + 1` only when
/// `delta + k` is **even**: [`band_for`] rounds the two ends toward each other,
/// so `hi - lo` is `k` on even parity and `k - 1` on odd, i.e.
/// `hi - lo ∈ {k - 1, k}`. That is not a nicety — `delta = 0, k = 1` gives
/// `lo = hi = 0`, a **one**-diagonal band, which is exactly the first rung of
/// the doubling ladder `the_storage_is_narrow_at_low_divergence_…` pins.
///
/// Using the diagonal count as the stride is a **pessimisation** whenever it
/// exceeds `m`: a `1000 x 1` block has `delta = 999` and so `k >= 999`, and a
/// diagonal-count stride would store ~1000 columns per row where the grid itself
/// has 2. Each row's span is clipped to the grid, and that clip
/// bounds it by `min(hi - lo, m) + 1` (proved by cases on which of the two
/// clamps bind — see the `debug_assert` in `row_span`). Taking the
/// minimum therefore makes this layout **never larger than the full grid**, and
/// exactly the full grid when the band covers it.
///
/// The rows are stored at a fixed stride rather than packed at their exact
/// spans: a packed layout needs a per-row prefix-sum table to address, which is
/// another `Θ(n)` allocation and a dependent load on every access, to save at
/// most `k` cells per row near the two corners.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BandLayout {
    ref_len: usize,
    alt_len: usize,
    /// Lower bound (inclusive) on `i - j`.
    lo: i64,
    /// Upper bound (inclusive) on `i - j`.
    hi: i64,
    /// Cells reserved per row.
    stride: usize,
}

impl BandLayout {
    /// The layout for the `(ref_len + 1) x (alt_len + 1)` grid restricted to
    /// diagonals `lo ..= hi`.
    ///
    /// The band must contain diagonal 0 — every band [`band_for`] produces for a
    /// `k >= |delta|` does, which is the only `k` the doubling search ever uses.
    fn new(ref_len: usize, alt_len: usize, lo: i64, hi: i64) -> Self {
        debug_assert!(lo <= 0 && hi >= 0, "the band must contain diagonal 0");
        let stride = ((hi - lo) as usize).min(alt_len) + 1;
        Self {
            ref_len,
            alt_len,
            lo,
            hi,
            stride,
        }
    }

    /// Cells to allocate for one grid under this layout.
    pub(crate) fn len(&self) -> usize {
        (self.ref_len + 1) * self.stride
    }

    /// The columns row `i` occupies inside the band, clamped to the grid. Never
    /// empty for a `k >= |delta|` band.
    fn row_span(&self, i: usize) -> (usize, usize) {
        let row = i as i64;
        let first = (row - self.hi).max(0) as usize;
        let last = (row - self.lo).min(self.alt_len as i64) as usize;
        debug_assert!(first <= last, "row {i} has an empty band span");
        debug_assert!(
            last - first < self.stride,
            "row {i} spans {} cells but the stride is {}",
            last - first + 1,
            self.stride
        );
        (first, last)
    }

    /// `index(i, j) - j`, hoisted out of a row's inner loop.
    fn row_base(&self, i: usize) -> usize {
        // `first(i) <= i <= i * stride` because `stride >= 1`, so this cannot
        // underflow.
        i * self.stride - (i as i64 - self.hi).max(0) as usize
    }

    /// Whether `(i, j)` is on the grid **and** inside the band. A cell outside
    /// the band has no slot at all under this layout, so this is the guard every
    /// read must pass first.
    pub(crate) fn contains(&self, i: usize, j: usize) -> bool {
        let d = i as i64 - j as i64;
        i <= self.ref_len && j <= self.alt_len && self.lo <= d && d <= self.hi
    }

    /// Slot of `(i, j)`. The caller must have established
    /// [`contains`](Self::contains); there is no in-range answer otherwise.
    pub(crate) fn index(&self, i: usize, j: usize) -> usize {
        debug_assert!(
            self.contains(i, j),
            "({i}, {j}) is outside the band {}..={} of a {} x {} grid",
            self.lo,
            self.hi,
            self.ref_len,
            self.alt_len
        );
        self.row_base(i) + j
    }
}

/// Levenshtein cost grid restricted to `band`'s diagonals, written into `grid`
/// (band-major, `band.len()` cells — see [`BandLayout`]).
///
/// Cells outside the band are not merely unwritten, they are **unallocated**, so
/// nothing may read one. Each predecessor is guarded by an explicit in-band test
/// rather than by pre-filling the grid with a sentinel: that pre-fill would put
/// back the `Θ(n·m)` pass the band exists to remove.
///
/// Every cell of `band` is written here in dependency order, so the buffer's
/// prior contents are irrelevant — which is what lets [`AlignmentDag::build`]
/// hand it a freshly sized buffer on each round of the doubling search.
fn fill_edit_grid_banded(reference: &[u8], alt: &[u8], band: &BandLayout, grid: &mut [u16]) {
    let n = reference.len();
    for i in 0..=n {
        let (first, last) = band.row_span(i);
        let base = band.row_base(i);
        let row = i as i64;
        for j in first..=last {
            if i == 0 && j == 0 {
                grid[base] = 0;
                continue;
            }
            let col = j as i64;
            let mut best = UNREACHED;
            // Diagonal: same `d`, so in band whenever this cell is.
            if i > 0 && j > 0 {
                let cost = u16::from(reference[i - 1] != alt[j - 1]);
                best = best.min(grid[band.index(i - 1, j - 1)].saturating_add(cost));
            }
            // Up (a reference base deleted). Its diagonal is `(i - 1) - j`,
            // one below this cell's, so it is in band iff `row - col > lo`.
            if i > 0 && row - col > band.lo {
                best = best.min(grid[band.index(i - 1, j)].saturating_add(1));
            }
            // Left (an alternate base inserted). Its diagonal is `i - (j - 1)`,
            // one above this cell's, so it is in band iff `row - col < hi`.
            if j > 0 && row - col < band.hi {
                best = best.min(grid[band.index(i, j - 1)].saturating_add(1));
            }
            grid[base + j] = best;
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
    /// **Storage is banded too, since #1988.** Every grid is `Θ((n + m) · k)`
    /// rather than `Θ(n · m)`, addressed through [`BandLayout`] — see there for
    /// the index arithmetic and for why the stride is `min(hi - lo, m) + 1`
    /// (`hi - lo` being `k` up to a parity adjustment, never simply `k`), which
    /// is what keeps this layout from ever being *larger* than the full grid. The
    /// doubling search reallocates as the band widens; the strides double, so
    /// the whole ladder costs less than twice its final allocation.
    ///
    /// Before that, the fill was banded and the allocation was not, and past a
    /// certain width the allocation was the entire cost: the speedup measured at
    /// `n = 4096, k = 2` was **7.9x against a ~1171x cell-count ratio**, i.e.
    /// ~101 MB of zeroing at ~5 GB/s and almost no DP. The module header carries
    /// that table and the argument it made.
    ///
    /// # Resident memory: MEASURED on 16 kB pages, still unmeasured on 4 kB
    ///
    /// Narrowing the *allocation* is not the same claim as lowering peak RSS —
    /// the untouched part of a zero-initialised `Vec` is never resident, so how
    /// much of the old full-size allocation was ever faulted in **depends on the
    /// page size**, and the geometry predicts opposite answers on the two
    /// platforms this project builds on. Only one of the two has been measured,
    /// so read the two bullets differently:
    ///
    /// - **16 kB pages** (Apple silicon) — **measured, and the prediction
    ///   held**. A grid row is 7.7 kB at the canonical window, so a row sits
    ///   *under* one page and consecutive rows share pages; the band crosses
    ///   every row, so every page was faulted in regardless and the whole
    ///   allocation was resident. Peak RSS for one `build` at
    ///   `n = m = 4096, k = 2`, release profile, four reps either side:
    ///   **106,758,144 B (101.8 MiB) before, 6,012,928 B (5.73 MiB) after** — a
    ///   17.8x reduction, ~96 MiB. Rep spread was 0.015% / 0.14%, and the
    ///   "before" figure is within 1% of the ~101 MB of grids the arithmetic
    ///   above predicts, which is what says they really were all resident. The
    ///   5.73 MiB that remains is the process floor — binary, runtime, the two
    ///   4 KiB blocks — not the DAG. **Say which grids, because the two answers
    ///   differ by 3x:** at this shape the layout is `4097 × 3 = 12,291` cells,
    ///   so the two grids the DAG *retains* (`on_optimal_path` and `out`, one
    ///   byte per cell each) are **24 KiB** together, while the four live
    ///   simultaneously inside `build` — those two plus `prefix` and
    ///   `suffix_grid` at two bytes per cell — are **72 KiB**. Both are noise
    ///   beside 5.73 MiB, which is the point being made; neither is 36 KiB,
    ///   which is what this sentence used to say and is not either quantity.
    /// - **4 kB pages** (the Linux host the timings above come from) — **still
    ///   unmeasured, and the win is expected to be smaller**. A row spans ~2
    ///   pages and the band is diagonal, so row `i` touches only *one* of them
    ///   and much of each grid may never have been faulted in the first place.
    ///   Do **not** carry the 16 kB figure onto that host: it is the platform
    ///   where the geometry predicts the least to reclaim.
    ///
    /// The measurement is peak RSS of a process doing one `build` and nothing
    /// else, so it isolates these grids and does not describe a normalizer run.
    /// It says nothing about **time**: no criterion run has been made for this
    /// change, and the module header's first bullet is the reason not to guess
    /// one from a local timing.
    pub fn build(ref_block: &[u8], alt_block: &[u8]) -> Option<Self> {
        let n = ref_block.len();
        let m = alt_block.len();
        if n.checked_add(m)? > MAX_ALIGNMENT_SPAN {
            return None;
        }
        let delta = (n as i64 - m as i64).unsigned_abs() as usize;
        // The band at which every diagonal is in play. Reaching it terminates
        // the loop on an exact answer whatever the doubling does, so the
        // convergence test below is an optimisation and not the only exit.
        let widest = (n + m).max(1);

        // No alignment costs less than the length difference, so that is the
        // narrowest band worth trying; `max(1)` keeps a diagonal of width one
        // for the equal-length case.
        let mut k = delta.max(1).min(widest);
        let (lo, hi) = band_for(n, m, k);
        let mut band = BandLayout::new(n, m, lo, hi);
        let mut prefix = vec![0u16; band.len()];
        loop {
            fill_edit_grid_banded(ref_block, alt_block, &band, &mut prefix);
            // Confining the DP can only overestimate, so a distance the band
            // itself admits is the true one: every cell of every minimal
            // alignment then lies inside the band, and the band-confined value
            // at each of them is its true value.
            if prefix[band.index(n, m)] as usize <= k || k >= widest {
                break;
            }
            k = (2 * k).min(widest);
            let (lo, hi) = band_for(n, m, k);
            band = BandLayout::new(n, m, lo, hi);
            // The band widened, so the buffer must too. Reallocated rather than
            // resized: the next fill writes every cell of the wider band in
            // dependency order, so nothing carried over would be read.
            prefix = vec![0u16; band.len()];
        }
        let total: u16 = prefix[band.index(n, m)];

        // Suffix distances, computed as prefix distances of the reversed
        // blocks: `suffix(i, j)` is the cost of aligning `ref[i..]` to
        // `alt[j..]`. Reversal maps diagonal `d` to `delta - d`, and the band is
        // symmetric about `delta / 2`, so both the band and its layout carry
        // over unchanged — `(n - i, m - j)` sits on diagonal `delta - d`, which
        // is in band exactly when `d` is.
        let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
        let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
        let mut suffix_grid = vec![0u16; band.len()];
        fill_edit_grid_banded(&ref_rev, &alt_rev, &band, &mut suffix_grid);
        let suffix = |i: usize, j: usize| suffix_grid[band.index(n - i, m - j)];

        let mut on_optimal_path = vec![false; band.len()];
        for i in 0..=n {
            let (first, last) = band.row_span(i);
            let base = band.row_base(i);
            for j in first..=last {
                on_optimal_path[base + j] = prefix[base + j].saturating_add(suffix(i, j)) == total;
            }
        }

        let mut out = vec![0u8; band.len()];
        // `needless_range_loop` fires on `ref_block[i]` below and its suggestion
        // is wrong: the loop runs `0..=n` over the grid's `n + 1` rows while
        // `ref_block` holds `n` bases, so the offered
        // `ref_block.iter().enumerate().take(n + 1)` silently drops the last row
        // — the one holding the sink. `i` indexes three things of two different
        // lengths here and is a grid coordinate, not a slice cursor.
        #[allow(clippy::needless_range_loop)]
        for i in 0..=n {
            let (first, last) = band.row_span(i);
            let base = band.row_base(i);
            for j in first..=last {
                if !on_optimal_path[base + j] {
                    continue;
                }
                let here = prefix[base + j];
                // Each neighbour is tested for band membership *first*. A cell
                // outside the band has no slot at all under this layout, so —
                // unlike the full-size grids this replaced, where an out-of-band
                // neighbour was merely an unwritten cell holding `false` — the
                // check is what keeps the index in range, not merely what keeps
                // an unwritten value from being read.
                // Diagonal: match or substitution. Same diagonal as this cell,
                // so in band whenever it is on the grid.
                if i < n && j < m {
                    let cost = u16::from(ref_block[i] != alt_block[j]);
                    let next = band.index(i + 1, j + 1);
                    if on_optimal_path[next] && here + cost == prefix[next] {
                        out[base + j] |= if cost == 0 { OUT_MATCH } else { OUT_SUB };
                    }
                }
                // Down: deletion of a reference base. One diagonal above.
                if i < n && band.contains(i + 1, j) {
                    let next = band.index(i + 1, j);
                    if on_optimal_path[next] && here + 1 == prefix[next] {
                        out[base + j] |= OUT_DEL;
                    }
                }
                // Right: insertion of an alternate base. One diagonal below.
                if j < m && band.contains(i, j + 1) {
                    let next = band.index(i, j + 1);
                    if on_optimal_path[next] && here + 1 == prefix[next] {
                        out[base + j] |= OUT_INS;
                    }
                }
            }
        }

        Some(Self {
            ref_len: n as u32,
            alt_len: m as u32,
            total: u32::from(total),
            band,
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

    /// How this DAG's cells are addressed, so a caller allocating a plane over
    /// the same cells can be narrowed by the same band rather than sized to the
    /// full grid. `CanonicalAlignment::of`'s `to_go` planes are the one such
    /// caller.
    pub(crate) fn band(&self) -> BandLayout {
        self.band
    }

    fn index(&self, i: u32, j: u32) -> usize {
        self.band.index(i as usize, j as usize)
    }

    /// Whether the cell lies on at least one minimal alignment.
    ///
    /// A cell outside the band is not on any minimal alignment (module header),
    /// and has no slot under this DAG's layout — so the band test is both the
    /// answer and the bounds check.
    pub(crate) fn contains_cell(&self, i: u32, j: u32) -> bool {
        self.band.contains(i as usize, j as usize) && self.on_optimal_path[self.index(i, j)]
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
        let (lo, hi) = (self.band.lo, self.band.hi);
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
    /// `(i, j)` should be within the grid (`i <= ref_len`, `j <= alt_len`); the
    /// `debug_assert` below catches a caller that strays. Off-grid and
    /// out-of-band both yield **no edges** in release rather than aliasing
    /// another cell — under the banded layout (#1988) every read is gated on
    /// [`BandLayout::contains`], and a cell outside the band is on no minimal
    /// alignment and so genuinely has no out-edges.
    pub(crate) fn out_edges(&self, i: u32, j: u32) -> impl Iterator<Item = (u32, u32, Step)> + '_ {
        debug_assert!(
            i <= self.ref_len && j <= self.alt_len,
            "out_edges({i}, {j}) is off-grid ({} x {})",
            self.ref_len,
            self.alt_len
        );
        // Reconstructed from the mask and `(i, j)`. The order matches what the
        // nested-`Vec` form pushed — diagonal, then down, then right — so any
        // caller that depended on edge order is unaffected.
        let mask = if self.band.contains(i as usize, j as usize) {
            self.out[self.index(i, j)]
        } else {
            0
        };
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
    /// `Θ((n + m) · k)` — rather than filtering the whole grid. The `in_degree`
    /// table is banded too since #1988: it is sized from `on_optimal_path`
    /// rather than from the grid, so it cannot regress to `Θ(n·m)` on its own.
    /// Deliberately not the per-edge reachability probe the calibration
    /// prototype used — that is `O(E)` per edge, which is fine for a 52 nt block
    /// and not for production. The probe survives as the test oracle.
    pub(crate) fn dominators(&self) -> Dominators {
        let mut result = Dominators::default();
        let sink = (self.ref_len, self.alt_len);
        let order: Vec<(u32, u32)> = self.cells().collect();

        // In-degree within the DAG, needed to know when a node's incoming edges
        // stop crossing the frontier.
        let mut in_degree = vec![0u32; self.on_optimal_path.len()];
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
    ///
    /// **What it establishes, and what it does not (#1959).** Every term here —
    /// `dag.cells()`, `dag.out_edges(..)`, the reachability walk — is read from
    /// the DAG; nothing is re-derived from the two sequences. So it proves only
    /// that [`AlignmentDag::dominators`] agrees with edge-banned reachability
    /// *over this graph*. It **cannot** detect a DAG that is uniformly too small
    /// — one missing cells that every minimal alignment of the sequences
    /// actually occupies — because a shrunken graph is internally consistent and
    /// both sides read the same missing cells. It catches asymmetric or
    /// inconsistent damage; it is blind to uniform shrinkage. What does catch
    /// that is the independent [`full_grid_build`], diffed against the shipped
    /// `build` by the `the_band_agrees_with_the_full_grid_*` tests. See
    /// [`the_brute_force_oracle_is_blind_to_a_uniformly_too_small_dag`], which
    /// pins both halves.
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

    /// Build the DAG confined to the band for a caller-supplied `k`, with **no**
    /// doubling search — the one thing [`AlignmentDag::build`] does that would
    /// correct a too-narrow band. When `k` is smaller than the true edit
    /// distance the band cannot hold every minimal alignment, so the returned
    /// DAG is *uniformly too small*: it is internally consistent (source reaches
    /// sink, `on_optimal_path` and `out` agree with its own `total`), yet it is
    /// missing cells the true DAG retains.
    ///
    /// This exists only to manufacture that artefact for
    /// [`the_brute_force_oracle_is_blind_to_a_uniformly_too_small_dag`],
    /// reproducing the #1959 sabotage exactly ("narrowing the band by one
    /// doubling step, accepting a band-confined distance where only a smaller
    /// one is sound"). It deliberately mirrors `build`'s banded fill so the
    /// shrunk DAG is a faithful narrowing rather than a hand-corrupted one; it
    /// is NOT an independent oracle — that role is [`full_grid_build`]'s.
    ///
    /// "Mirrors" is a claim, and a hand copy has no mechanism of its own to keep
    /// it true, so it is pinned rather than asserted:
    /// [`the_hand_copied_band_confined_build_still_matches_build`] requires that
    /// at a band wide enough to hold every minimal alignment this reproduces
    /// `build` exactly.
    fn build_confined_to_band(ref_block: &[u8], alt_block: &[u8], k: usize) -> AlignmentDag {
        let n = ref_block.len();
        let m = alt_block.len();
        let (lo, hi) = band_for(n, m, k);
        let band = BandLayout::new(n, m, lo, hi);
        let mut prefix = vec![0u16; band.len()];
        fill_edit_grid_banded(ref_block, alt_block, &band, &mut prefix);
        let total = prefix[band.index(n, m)];
        let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
        let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
        let mut suffix_grid = vec![0u16; band.len()];
        fill_edit_grid_banded(&ref_rev, &alt_rev, &band, &mut suffix_grid);
        let suffix = |i: usize, j: usize| suffix_grid[band.index(n - i, m - j)];
        let mut on_optimal_path = vec![false; band.len()];
        for i in 0..=n {
            let (first, last) = band.row_span(i);
            let base = band.row_base(i);
            for j in first..=last {
                on_optimal_path[base + j] = prefix[base + j].saturating_add(suffix(i, j)) == total;
            }
        }
        let mut out = vec![0u8; band.len()];
        // `i` is a grid coordinate here, not a slice cursor — the loop runs the
        // grid's `n + 1` rows while `ref_block` holds `n` bases — so the
        // `needless_range_loop` suggestion silently drops the sink row, exactly
        // as it does on the same loop in `build`.
        #[allow(clippy::needless_range_loop)]
        for i in 0..=n {
            let (first, last) = band.row_span(i);
            let base = band.row_base(i);
            for j in first..=last {
                if !on_optimal_path[base + j] {
                    continue;
                }
                let here = prefix[base + j];
                if i < n && j < m {
                    let cost = u16::from(ref_block[i] != alt_block[j]);
                    let next = band.index(i + 1, j + 1);
                    if on_optimal_path[next] && here + cost == prefix[next] {
                        out[base + j] |= if cost == 0 { OUT_MATCH } else { OUT_SUB };
                    }
                }
                if i < n && band.contains(i + 1, j) {
                    let next = band.index(i + 1, j);
                    if on_optimal_path[next] && here + 1 == prefix[next] {
                        out[base + j] |= OUT_DEL;
                    }
                }
                if j < m && band.contains(i, j + 1) {
                    let next = band.index(i, j + 1);
                    if on_optimal_path[next] && here + 1 == prefix[next] {
                        out[base + j] |= OUT_INS;
                    }
                }
            }
        }
        AlignmentDag {
            ref_len: n as u32,
            alt_len: m as u32,
            total: u32::from(total),
            band,
            on_optimal_path,
            out,
        }
    }

    /// #1959. The brute-force oracle
    /// ([`dominators_by_brute_force`]) is computed **from the DAG it is handed**
    /// — its cell set, its edges, its reachability — and re-derives nothing from
    /// the two sequences. So it establishes only that
    /// [`AlignmentDag::dominators`] agrees with edge-banned reachability *over
    /// that graph*, and it cannot see the graph itself being uniformly too
    /// small: a shrunken DAG is internally consistent, so both sides read the
    /// same missing cells and the comparison is satisfied.
    ///
    /// This pins the limitation as an executable fact and pins what closes it:
    ///
    /// 1. **The blind spot is real.** On the issue's own `AC -> CA`, a band one
    ///    doubling too narrow retains 3 of the true DAG's 7 cells, and the
    ///    oracle still agrees with `dominators()` on it — the same agreement
    ///    [`dominators_agree_with_brute_force_on_exhaustive_small_inputs`]
    ///    treats as a pass. The shrinkage can even leave the dominator output
    ///    *correct* (`AC -> CA`) or make it *wrong* (`CAG -> AGA` loses every
    ///    dominator); the oracle cannot tell either from the truth.
    /// 2. **What closes it is the INDEPENDENT implementation, not the oracle.**
    ///    [`full_grid_build`] recomputes the DAG from the sequences, and the
    ///    `the_band_agrees_with_the_full_grid_*` tests diff the shipped `build`
    ///    against it cell by cell. That comparison detects the shrinkage the
    ///    oracle misses — which is why a DAG-construction change that lands
    ///    without such a second implementation to diff against would not be
    ///    protected (issue suggestion 3).
    #[test]
    fn the_brute_force_oracle_is_blind_to_a_uniformly_too_small_dag() {
        // (reference, alt, too-narrow `k`, the shrunk DAG is a strict SUBSET of
        // the truth, `dominators()` on it still MATCHES the truth, the edit
        // DISTANCE on it still matches the truth).
        for (r, a, k, subset_of_truth, dominators_stay_correct, distance_stays_correct) in [
            // The issue's verbatim example. At k = 1 the band holds only the
            // two-substitution diagonal, a strict subset of the seven cells the
            // true (k = 2) DAG spans — and with the SAME distance (2) and the
            // SAME (empty) dominators, so nothing derived from the DAG alone,
            // not even `dominators()`, reveals the shrinkage. The subtlest form.
            (&b"AC"[..], &b"CA"[..], 1usize, true, true, true),
            // The dangerous form: at k = 1 the band collapses to the
            // all-substitution diagonal (distance 3, not the true 2), so
            // `dominators()` comes back EMPTY where the truth is
            // {(1,0),(2,1)} + forced junction 3 — a wrong answer the oracle
            // nonetheless blesses. Not a subset (its cells lie off every
            // minimal alignment), so only its size is compared below.
            (&b"CAG"[..], &b"AGA"[..], 1, false, false, false),
            // The two rows above demonstrate the oracle's agreement as
            // `empty == empty`: `dominators()` is `{}` on both sides of the
            // first, and `{}` on the shrunk side of the second. A `dominators()`
            // that returned nothing at all, for any reason, would satisfy both.
            // These two rows carry NON-EMPTY dominator sets, so the agreement
            // has content.
            //
            // The first is `AC -> CA`'s shape at a size that has dominators:
            // 12 of the truth's 15 cells, same distance (4), and the same
            // non-empty `{(0,0)}` on both sides — every derived quantity agrees
            // and the graph is still too small.
            (&b"AAAC"[..], &b"ACCCAA"[..], 2, true, true, true),
            // The second is `CAG -> AGA`'s shape with a non-empty WRONG answer,
            // and the only fixture here that exercises `forced_ins_junctions`
            // at all: 8 of the truth's 14 cells, and `dominators()` invents both
            // a matched dominator `{(2,2)}` and a forced insertion junction
            // `{4}` where the truth has NEITHER. Unlike `CAG -> AGA` the
            // distance does not give it away — 4 on both sides — so this is a
            // wrong non-empty answer with no tell anywhere in the DAG.
            (&b"AACA"[..], &b"CCCAAC"[..], 2, true, false, true),
        ] {
            let shrunk = build_confined_to_band(r, a, k);
            let truth = full_grid_build(r, a);
            let label = label_for(r, a);

            // (1) It really is too small: fewer retained cells than the truth,
            // and — for the issue's own case — a strict subset of them.
            let shrunk_cells: Vec<_> = shrunk.cells().collect();
            let true_cells: Vec<_> = truth.cells().collect();
            assert!(
                shrunk_cells.len() < true_cells.len(),
                "{label}: the manufactured DAG must retain fewer cells than the truth \
                 ({} vs {})",
                shrunk_cells.len(),
                true_cells.len()
            );
            if subset_of_truth {
                assert!(
                    shrunk_cells.iter().all(|c| true_cells.contains(c)),
                    "{label}: the shrunk cells must be a strict subset of the true ones"
                );
            }

            // (2) THE BLIND SPOT. The brute-force oracle, reading only the
            // shrunk DAG, agrees with `dominators()` on it: both are pure
            // functions of the DAG, so the oracle is SATISFIED by a graph it
            // should reject. This is exactly the check
            // `dominators_agree_with_brute_force_on_exhaustive_small_inputs`
            // runs, reproduced on a deliberately shrunken graph.
            let (oracle_matched, oracle_ins) = dominators_by_brute_force(&shrunk);
            let produced = shrunk.dominators();
            assert_eq!(
                (
                    produced.matched.clone(),
                    produced.forced_ins_junctions.clone()
                ),
                (oracle_matched, oracle_ins),
                "{label}: the brute-force oracle must AGREE with dominators() on the \
                 shrunk DAG — that agreement IS the blind spot"
            );

            // (3) And yet the shrunk dominator answer may be correct (AC -> CA)
            // or wrong (CAG -> AGA). Either way (2) held, so the oracle could
            // not distinguish it from the truth.
            assert_eq!(
                shrunk.dominators() == truth.dominators(),
                dominators_stay_correct,
                "{label}: unexpected dominator (dis)agreement with the truth"
            );

            // (4) WHAT CLOSES THE GAP — and, more usefully, WHICH ARM of it
            // does the closing.
            //
            // The comparison that catches uniform shrinkage is the shipped
            // `build` against the independent [`full_grid_build`], which
            // `the_band_agrees_with_the_full_grid_*` make. Run it here, on these
            // very inputs: it is green, so the closing comparison is genuinely
            // available at the locus where the oracle above is blind.
            assert_band_matches_full_grid(r, a);

            // That comparison has three arms — the edit distance, the retained
            // cell set with its out-edges, and `dominators()` — and the rows
            // above are chosen so that on some of them two arms are blind. (3)
            // has already pinned the dominator arm as blind wherever
            // `dominators_stay_correct`; this pins the distance arm the same
            // way, which the prose above claims ("the SAME distance (2)") and
            // nothing asserted.
            //
            // On the two rows where both flags are `true` — `AC -> CA` and
            // `AAAC -> ACCCAA` — that leaves the retained-cell/out-edge arm as
            // the ONLY one that can see the narrowing, and (1) has already
            // shown it does. Which is why a DAG-construction change landing
            // without a second implementation to diff against would be
            // unprotected (issue suggestion 3).
            //
            // This replaces an `assert_ne!(shrunk_cells, true_cells)` that could
            // never fail: (1) asserts the two lengths differ, and `Vec`'s
            // equality compares lengths first, so the inequality was already
            // established by the assertion above it.
            assert_eq!(
                shrunk.edit_distance() == truth.edit_distance(),
                distance_stays_correct,
                "{label}: unexpected edit-distance (dis)agreement with the truth ({} vs {})",
                shrunk.edit_distance(),
                truth.edit_distance()
            );
        }
    }

    /// [`build_confined_to_band`] is a hand copy of [`AlignmentDag::build`]'s
    /// post-loop body, kept so the shrunk DAG in
    /// [`the_brute_force_oracle_is_blind_to_a_uniformly_too_small_dag`] is a
    /// faithful narrowing rather than a hand-corrupted one. **Nothing about that
    /// role makes the copy track the original**: it asserts properties of the
    /// *shrunk* DAG — fewer cells than the truth, particular dominators — and
    /// never agreement with `build`, so if `build`'s fill, its `out`-bit
    /// semantics or `BandLayout` change, whether that test notices is incidental.
    ///
    /// **Do not read that as "it would stay green".** Three deliberate drifts
    /// were applied to the copy — dropping its `OUT_INS` arm, truncating its
    /// `out` loop from `0..=n` to `0..n`, and swapping its `OUT_DEL`/`OUT_INS`
    /// bits — and it reddened on all three. What it did *not* do on any of them
    /// is say why: it reports `a single crossing must be this cell's only
    /// out-edge`, `out_edges(5, 4) is off-grid (4 x 6)` or a dominator
    /// disagreement, every one of which reads as `dominators()` or the DAG being
    /// wrong rather than as the copy having drifted. This pin answers with
    /// `out-edges differ at (1, 2)`, which is the actual fault.
    ///
    /// So the value here is attribution plus a guarantee the other test does not
    /// offer — not coverage of a hole measured open.
    ///
    /// The guarantee: handed a band wide enough to contain every minimal
    /// alignment, the copy must produce exactly what the shipped `build`
    /// produces. `n + m` is that band — it is `build`'s own `widest`, the width
    /// at which its doubling loop terminates unconditionally — so the two are
    /// being asked the same question and only a drift in the shared post-loop
    /// body can separate them.
    ///
    /// Deliberately compared against `build` and not against
    /// [`full_grid_build`]: `full_grid_build` is the *pre-band* implementation,
    /// so a copy that agreed with it could still have drifted from every
    /// band-specific detail the copy exists to reproduce.
    #[test]
    fn the_hand_copied_band_confined_build_still_matches_build() {
        for (r, a) in [
            // Every fixture the blind-spot test above narrows, so a drift is
            // caught at exactly the inputs that test depends on.
            (&b"AC"[..], &b"CA"[..]),
            (&b"CAG"[..], &b"AGA"[..]),
            (&b"AAAC"[..], &b"ACCCAA"[..]),
            (&b"AACA"[..], &b"CCCAAC"[..]),
            // Plus the shapes a band is most likely to get wrong, since the
            // copy carries its own `band.contains` / `row_span` arithmetic:
            // both blocks empty, either block empty, and a length ratio far
            // from one so `|delta|` and not the divergence sets the band.
            (&b""[..], &b""[..]),
            (&b""[..], &b"ACGT"[..]),
            (&b"ACGT"[..], &b""[..]),
            (&b"A"[..], &b"ACGTACGT"[..]),
            (&b"ACGTACGT"[..], &b"A"[..]),
            (&b"AAAAA"[..], &b"AAAAA"[..]),
        ] {
            let widest = (r.len() + a.len()).max(1);
            let copy = build_confined_to_band(r, a, widest);
            let shipped =
                AlignmentDag::build(r, a).expect("these fixtures are far below MAX_ALIGNMENT_SPAN");
            assert_dags_agree(&copy, &shipped, &label_for(r, a));
        }
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
            // and the two walks are compared like for like. `BandLayout`'s
            // stride collapses to `m + 1` for a band this wide, so the flat
            // `i * width + j` vectors built above are addressed unchanged —
            // which is what keeps this oracle the *pre-band* implementation
            // rather than a re-expression of the banded one.
            band: BandLayout::new(n, m, -(m as i64), n as i64),
            on_optimal_path,
            out,
        }
    }

    /// Assert the banded build and the pre-band build agree on everything the
    /// DAG exposes: the distance, the retained cell set, every out-edge mask,
    /// the topological walk, and `dominators()`. Returns the edit distance so a
    /// caller can assert which `k` it actually exercised.
    ///
    /// The cell set and the edge masks are compared **over every cell of the
    /// full `(n + 1) x (m + 1)` grid**, through the public accessors, rather
    /// than by diffing the backing `Vec`s. Since #1988 the two builds do not
    /// share a layout — that is the change — so a `Vec` diff would compare
    /// storage rather than content and would report a difference for a
    /// narrowing that lost nothing. Reading every grid cell is the stronger
    /// question anyway: it additionally pins that a cell the band excludes
    /// answers `false` and yields no edges, which a `Vec` diff never asked.
    fn assert_band_matches_full_grid(r: &[u8], a: &[u8]) -> u32 {
        // `build` refuses past `MAX_ALIGNMENT_SPAN` (#1970); every fixture here is
        // far below it, so a `None` is a bug in the fixture, not a decline to test.
        let banded = AlignmentDag::build(r, a)
            .expect("differential fixtures are far below MAX_ALIGNMENT_SPAN");
        let full = full_grid_build(r, a);
        assert_dags_agree(&banded, &full, &label_for(r, a));
        banded.total
    }

    /// `"<reference>" -> "<alternate>"`, the label every differential assertion
    /// in this module reports a failure against.
    fn label_for(r: &[u8], a: &[u8]) -> String {
        format!(
            "{:?} -> {:?}",
            String::from_utf8_lossy(r),
            String::from_utf8_lossy(a)
        )
    }

    /// Assert two DAGs agree on everything the type exposes: the distance, the
    /// retained cell set, every out-edge mask, the topological walk, and
    /// `dominators()`.
    ///
    /// Extracted from [`assert_band_matches_full_grid`] so the same comparison
    /// can be aimed at a second pair of builds
    /// ([`the_hand_copied_band_confined_build_still_matches_build`]) without a
    /// second copy of the sweep — a copied comparison drifting from the one it
    /// was copied from is the very failure that test exists to prevent.
    ///
    /// Everything is read **through the public accessors, over every cell of the
    /// full `(n + 1) x (m + 1)` grid**, rather than by diffing the backing
    /// `Vec`s. The two sides need not share a layout — since #1988 the banded
    /// and full-grid builds do not — so a `Vec` diff would compare storage
    /// rather than content and would report a difference for a narrowing that
    /// lost nothing. Reading every grid cell is the stronger question anyway: it
    /// additionally pins that a cell either band excludes answers `false` and
    /// yields no edges, which a `Vec` diff never asked.
    fn assert_dags_agree(left: &AlignmentDag, right: &AlignmentDag, label: &str) {
        assert_eq!(left.total, right.total, "edit distance differs on {label}");
        assert_eq!(
            left.ref_len(),
            right.ref_len(),
            "ref_len differs on {label}"
        );
        assert_eq!(
            left.alt_len(),
            right.alt_len(),
            "alt_len differs on {label}"
        );
        for i in 0..=left.ref_len() {
            for j in 0..=left.alt_len() {
                assert_eq!(
                    left.contains_cell(i, j),
                    right.contains_cell(i, j),
                    "retained cell set differs at ({i}, {j}) on {label}"
                );
                assert_eq!(
                    left.out_edges(i, j).collect::<Vec<_>>(),
                    right.out_edges(i, j).collect::<Vec<_>>(),
                    "out-edges differ at ({i}, {j}) on {label}"
                );
            }
        }
        assert_eq!(
            left.cells().collect::<Vec<_>>(),
            right.cells().collect::<Vec<_>>(),
            "topological cell walk differs on {label}"
        );
        assert_eq!(
            left.dominators(),
            right.dominators(),
            "dominators differ on {label}"
        );
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

    /// The **storage** must be narrow, not merely the fill (#1988).
    ///
    /// This is the storage analogue of the guard directly above, and it exists
    /// for the same reason: every `the_band_agrees_with_the_full_grid_*` test
    /// asserts that the banded build *agrees* with the pre-band one, and a build
    /// that allocated full-size grids and addressed them at
    /// `i * (m + 1) + j` satisfies all of them trivially — that is exactly what
    /// this module did before #1988, and it was correct then. So a regression to
    /// full-size storage is invisible to the entire correctness suite, and the
    /// only other signal is a criterion run, which CI does not do.
    ///
    /// Checked against a deliberate sabotage rather than assumed: pinning
    /// `BandLayout::new`'s stride to `alt_len + 1` — which is behaviour-
    /// preserving, since a wider stride only leaves unused gaps between rows —
    /// leaves the `seqfirst` modules at **42 passed, 1 failed**, this test being
    /// the one failure (`on_optimal_path holds 1050625 cells ... expected
    /// 3075`), with all five differential tests, the brute-force dominator
    /// oracle and the band-narrowness guard above still green.
    #[test]
    fn the_storage_is_narrow_at_low_divergence_so_a_no_op_regression_is_caught() {
        // Two substitutions, spaced far apart, so the edit distance is exactly 2
        // and no cheaper indel-bearing alignment exists. That fixes the whole
        // doubling ladder — `delta = 0` seeds `k = 1`, whose one-diagonal band
        // reports distance 2, so `k` doubles once to 2 and converges — and hence
        // fixes the final band at three diagonals. The expected size is
        // therefore an exact figure rather than a bound.
        let mut state = 0x1988_u64 | 1;
        let reference: Vec<u8> = (0..1024)
            .map(|_| b"ACGT"[(next(&mut state) as usize) % 4])
            .collect();
        let mut alt = reference.clone();
        for at in [3usize, 700] {
            alt[at] = match alt[at] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        let dag = AlignmentDag::build(&reference, &alt)
            .expect("1024 + 1024 is far below MAX_ALIGNMENT_SPAN");
        assert_eq!(dag.edit_distance(), 2, "the fixture no longer has k = 2");

        let rows = reference.len() + 1;
        let grid = rows * (alt.len() + 1);
        let expected = rows * 3;
        for (what, stored) in [
            ("on_optimal_path", dag.on_optimal_path.len()),
            ("out", dag.out.len()),
            ("the layout the sweeps allocate from", dag.band().len()),
        ] {
            assert_eq!(
                stored, expected,
                "{what} holds {stored} cells for a 3-diagonal band over {rows} rows; \
                 expected {expected}, and the full grid would be {grid}"
            );
        }

        // The same property stated as the one the issue is about, so a future
        // change to the exact stride cannot quietly restore grid-scale storage
        // while keeping the equalities above passing by moving `expected`.
        assert!(
            expected * 100 < grid,
            "storage is {expected} of {grid} cells — that is grid-scale"
        );
    }

    /// Storage must never exceed the full grid either — the direction a
    /// diagonal-count stride gets wrong.
    ///
    /// `hi - lo + 1` is the obvious stride and is a **pessimisation** wherever
    /// `k > m`: a `1000 x 1` block has `delta = 999`, so the band is 1000
    /// diagonals wide while the grid has 2 columns, and storing one row per
    /// diagonal would be ~500x *larger* than the thing it replaces. Degenerate
    /// shapes are not hypothetical here — `MAX_ALIGNMENT_SPAN` exists because a
    /// long literal `ins` payload arrives as exactly this (#1970).
    #[test]
    fn storage_never_exceeds_the_full_grid_on_a_degenerate_shape() {
        for (n, m) in [(1000usize, 1usize), (1, 1000), (512, 0), (0, 512), (0, 0)] {
            let reference = vec![b'A'; n];
            let alt = vec![b'C'; m];
            let dag = AlignmentDag::build(&reference, &alt).expect("far below MAX_ALIGNMENT_SPAN");
            let grid = (n + 1) * (m + 1);
            assert!(
                dag.band().len() <= grid,
                "a {n} x {m} block stores {} cells, more than the full grid's {grid}",
                dag.band().len()
            );
        }
    }
}
