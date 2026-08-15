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
/// comment: past the bound `edit_grid`'s boundary row would truncate `j as u16`
/// and the sums would wrap, and `on_optimal_path` would be computed from wrong
/// costs. That is a silently wrong partition and so a silently wrong emitted
/// description — strictly worse than declining. Every gate in front of `build`
/// bounds **cells**, i.e. `(n + 1) * (m + 1)`, and a degenerate `1 x 70,000`
/// grid clears the widest of those by four orders of magnitude while `n + m` is
/// past `u16::MAX`; a long literal `ins` payload reaches exactly that shape.
/// See #1970.
pub(crate) const MAX_ALIGNMENT_SPAN: usize = u16::MAX as usize;

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
    /// Cost is `O(ref_len * alt_len)` in time and space.
    ///
    /// `None` when `ref_block.len() + alt_block.len()` exceeds
    /// [`MAX_ALIGNMENT_SPAN`], which is the width of the `u16` cost grid. Every
    /// caller in the normalizer already has a decline path, and declining is the
    /// only safe answer: past the bound the costs wrap rather than panic in the
    /// shipped profile. See [`MAX_ALIGNMENT_SPAN`].
    pub fn build(ref_block: &[u8], alt_block: &[u8]) -> Option<Self> {
        let n = ref_block.len();
        let m = alt_block.len();
        if n.checked_add(m)? > MAX_ALIGNMENT_SPAN {
            return None;
        }
        let prefix = edit_grid(ref_block, alt_block);

        // Suffix distances, computed as prefix distances of the reversed
        // blocks: `suffix(i, j)` is the cost of aligning `ref[i..]` to
        // `alt[j..]`.
        let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
        let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
        let suffix_grid = edit_grid(&ref_rev, &alt_rev);
        let suffix = |i: usize, j: usize| suffix_grid[(n - i) * (m + 1) + (m - j)];

        let total: u16 = prefix[n * (m + 1) + m];
        let width = m + 1;

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
                // Diagonal: match or substitution.
                if i < n && j < m {
                    let cost = u16::from(ref_block[i] != alt_block[j]);
                    let next = (i + 1) * width + (j + 1);
                    if here + cost == prefix[next] && on_optimal_path[next] {
                        out[i * width + j] |= if cost == 0 { OUT_MATCH } else { OUT_SUB };
                    }
                }
                // Down: deletion of a reference base.
                if i < n {
                    let next = (i + 1) * width + j;
                    if here + 1 == prefix[next] && on_optimal_path[next] {
                        out[i * width + j] |= OUT_DEL;
                    }
                }
                // Right: insertion of an alternate base.
                if j < m {
                    let next = i * width + (j + 1);
                    if here + 1 == prefix[next] && on_optimal_path[next] {
                        out[i * width + j] |= OUT_INS;
                    }
                }
            }
        }

        Some(Self {
            ref_len: n as u32,
            alt_len: m as u32,
            total: u32::from(total),
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
    pub(crate) fn cells(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        (0..=self.ref_len + self.alt_len)
            .flat_map(move |d| {
                let i_min = d.saturating_sub(self.alt_len);
                let i_max = d.min(self.ref_len);
                (i_min..=i_max).map(move |i| (i, d - i))
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
    /// `Θ(n·m)` — no worse than `build()`'s own cost. `cells()` walks the full
    /// `(ref_len + 1) x (alt_len + 1)` grid once to filter the on-path cells
    /// (`V` of them) in topological order, and the in-degree table below is
    /// sized to that same grid, so this pass cannot beat `Θ(n·m)` regardless of
    /// how few cells are actually on-path; on top of that it is `O(V + E)` to
    /// sweep. Deliberately not the per-edge reachability probe the calibration
    /// prototype used — that is `O(E)` per edge, which is fine for a 52 nt
    /// block and not for production. The probe survives as the test oracle.
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
            // cell put a heap allocation back on the `Θ(n·m)` sweep that the flat
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

/// Full Levenshtein distance grid, row-major with `alt.len() + 1` columns.
fn edit_grid(reference: &[u8], alt: &[u8]) -> Vec<u16> {
    let n = reference.len();
    let m = alt.len();
    // A cost never exceeds `n + m`, so the `u16` grid holds every value exactly
    // when `n + m <= MAX_ALIGNMENT_SPAN`. This is a statement of the
    // precondition `AlignmentDag::build` already enforces by refusing, not a
    // second gate — `edit_grid` is private and `build` is its only caller. Kept
    // so a future caller added here is caught in a debug build; the refusal, not
    // this line, is what makes the narrowing safe in the shipped profile.
    debug_assert!(
        n + m <= MAX_ALIGNMENT_SPAN,
        "edit_grid narrowed to u16: n+m = {} exceeds {MAX_ALIGNMENT_SPAN}",
        n + m
    );
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
}
