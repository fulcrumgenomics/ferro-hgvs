//! Unit 2: build the DAG of every minimal-cost alignment.
//!
//! Global alignment (both ends anchored) under **unit-cost Levenshtein**. Both
//! choices are deliberate:
//!
//! - *Global*, because the flanks were trimmed by the caller and are identical
//!   by construction, so there is no "where do these correspond" question and
//!   every base must be accounted for.
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
    total: u32,
    /// Row-major `(ref_len + 1) x (alt_len + 1)`; `true` when the cell lies on
    /// at least one minimal alignment.
    on_optimal_path: Vec<bool>,
    /// Row-major, parallel to `on_optimal_path`. Out-edges of each cell.
    out: Vec<Vec<(u32, u32, Step)>>,
}

impl AlignmentDag {
    /// Build the DAG for `ref_block` -> `alt_block`.
    ///
    /// Both slices are raw bases with common flanks already trimmed. Either may
    /// be empty. Cost is `O(ref_len * alt_len)` in time and space.
    pub fn build(ref_block: &[u8], alt_block: &[u8]) -> Self {
        let n = ref_block.len();
        let m = alt_block.len();
        let prefix = edit_grid(ref_block, alt_block);

        // Suffix distances, computed as prefix distances of the reversed
        // blocks: `suffix(i, j)` is the cost of aligning `ref[i..]` to
        // `alt[j..]`.
        let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
        let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
        let suffix_grid = edit_grid(&ref_rev, &alt_rev);
        let suffix = |i: usize, j: usize| suffix_grid[(n - i) * (m + 1) + (m - j)];

        let total = prefix[n * (m + 1) + m];
        let width = m + 1;

        let mut on_optimal_path = vec![false; (n + 1) * width];
        for i in 0..=n {
            for j in 0..=m {
                on_optimal_path[i * width + j] = prefix[i * width + j] + suffix(i, j) == total;
            }
        }

        let mut out: Vec<Vec<(u32, u32, Step)>> = vec![Vec::new(); (n + 1) * width];
        for i in 0..=n {
            for j in 0..=m {
                if !on_optimal_path[i * width + j] {
                    continue;
                }
                let here = prefix[i * width + j];
                // Diagonal: match or substitution.
                if i < n && j < m {
                    let cost = u32::from(ref_block[i] != alt_block[j]);
                    let next = (i + 1) * width + (j + 1);
                    if here + cost == prefix[next] && on_optimal_path[next] {
                        let step = if cost == 0 { Step::Match } else { Step::Sub };
                        out[i * width + j].push((i as u32 + 1, j as u32 + 1, step));
                    }
                }
                // Down: deletion of a reference base.
                if i < n {
                    let next = (i + 1) * width + j;
                    if here + 1 == prefix[next] && on_optimal_path[next] {
                        out[i * width + j].push((i as u32 + 1, j as u32, Step::Del));
                    }
                }
                // Right: insertion of an alternate base.
                if j < m {
                    let next = i * width + (j + 1);
                    if here + 1 == prefix[next] && on_optimal_path[next] {
                        out[i * width + j].push((i as u32, j as u32 + 1, Step::Ins));
                    }
                }
            }
        }

        Self {
            ref_len: n as u32,
            alt_len: m as u32,
            total,
            on_optimal_path,
            out,
        }
    }

    /// Minimal edit distance between the two blocks.
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
    /// ascending). Every edge strictly increases `i + j`, so this ordering is a
    /// valid topological sort.
    pub(crate) fn cells(&self) -> impl Iterator<Item = (u32, u32)> + '_ {
        let mut cells: Vec<(u32, u32)> = (0..=self.ref_len)
            .flat_map(move |i| (0..=self.alt_len).map(move |j| (i, j)))
            .filter(|&(i, j)| self.contains_cell(i, j))
            .collect();
        cells.sort_by_key(|&(i, j)| (i + j, i));
        cells.into_iter()
    }

    /// Out-edges of a cell, as `(next_i, next_j, step)`.
    pub(crate) fn out_edges(&self, i: u32, j: u32) -> impl Iterator<Item = (u32, u32, Step)> + '_ {
        self.out[self.index(i, j)].iter().copied()
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
    /// `O(V + E)`. Deliberately not the per-edge reachability probe the
    /// calibration prototype used — that is `O(E)` per edge, which is fine for
    /// a 52 nt block and not for production. The probe survives as the test
    /// oracle.
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
            let outs: Vec<(u32, u32, Step)> = self.out_edges(i, j).collect();
            crossing += outs.len() as i64;

            // If exactly one edge crosses the frontier now, it is on every
            // path. Any node with two or more out-edges contributes two or more
            // crossings itself, so a count of one means this node has exactly
            // one out-edge and nothing older is still crossing.
            if crossing == 1 && (i, j) != sink {
                debug_assert_eq!(
                    outs.len(),
                    1,
                    "a single crossing must be this cell's only out-edge"
                );
                let (_, _, step) = outs[0];
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
fn edit_grid(reference: &[u8], alt: &[u8]) -> Vec<u32> {
    let n = reference.len();
    let m = alt.len();
    let width = m + 1;
    let mut d = vec![0u32; (n + 1) * width];
    for (i, x) in d.iter_mut().step_by(width).enumerate() {
        *x = i as u32;
    }
    for (j, x) in d.iter_mut().take(width).enumerate() {
        *x = j as u32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let cost = u32::from(reference[i - 1] != alt[j - 1]);
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
            let dag = AlignmentDag::build(r, a);
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
        let dag = AlignmentDag::build(b"ACGTACGT", b"ATGTTCGT");
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
        let dag = AlignmentDag::build(b"ACGTACGTAC", b"ACGTTACGTGAC");
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
        let dag = AlignmentDag::build(
            b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT",
            b"TTCCTCGATGCCTG",
        );
        let dom = dag.dominators();
        assert_eq!(
            dom.matched_ref(),
            vec![48],
            "LRG_199 must have exactly one dominator match, at reference offset 48"
        );
        assert!(
            dom.alt_at(48).is_some(),
            "the dominator must pin an alternate coordinate too"
        );
    }

    #[test]
    fn a_pure_insertion_is_seen_as_a_forced_junction() {
        // #1260: two `insC` at adjacent gaps in a poly-A run. A node model
        // reports this as zero members because no reference position changes;
        // the junctions are the whole signal.
        let dag = AlignmentDag::build(b"AAAAAAA", b"AACACAAAA");
        let dom = dag.dominators();
        assert_eq!(
            dom.matched_ref().len(),
            7,
            "every reference position is matched in every minimal alignment"
        );
        assert!(
            !dom.forced_ins_junctions.is_empty(),
            "the insertions must show up as forced junctions, or they vanish"
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
                let dag = AlignmentDag::build(r, a);
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
}
