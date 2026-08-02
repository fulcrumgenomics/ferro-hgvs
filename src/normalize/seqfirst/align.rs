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
pub(crate) struct AlignmentDag {
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
    pub(crate) fn build(ref_block: &[u8], alt_block: &[u8]) -> Self {
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
    pub(crate) fn edit_distance(&self) -> u32 {
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
}
