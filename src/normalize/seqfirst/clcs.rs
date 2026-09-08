//! The 2025 cLCS-graph canonical-variant rule, as published, plus one named
//! extension.
//!
//! Santcroos, Kosters, Lefter, Laros & Vis, *"A Graph-based Approach to Variant
//! Description Extraction from Sequences"*, NAR Genomics and Bioinformatics
//! 7(4), Dec 2025 (arXiv 2503.18472v2); reference implementation
//! `github.com/mutalyzer/algebra` v1.5.2.
//!
//! # What the paper's rule is
//!
//! §3.3, verbatim:
//!
//! > We consider the subgraph of the cLCS-graph consisting of only the shortest
//! > paths. All edges have weight 1 except for unlabeled (λ) edges which have
//! > weight 0. … If there is no unique shortest path, we employ the extraction
//! > method from Section 3.2 on the shortest paths subgraph.
//!
//! §3.2 is the **local-supremal** collapse: the nodes lying on *every*
//! source→sink path (post-dominators) partition the graph, and each partition
//! becomes one coarser replacement.
//!
//! So the pipeline is three stages, and all three are implemented here:
//!
//! 1. the graph of all minimal alignments under the **simple edit distance**;
//! 2. restrict to the shortest-path subgraph under "1 per labeled replacement";
//! 3. collapse that subgraph at its post-dominators.
//!
//! # The one thing that is NOT ferro's existing DAG
//!
//! Ferro's [`super::align::AlignmentDag`] is **unit-cost Levenshtein**: a
//! substitution is one edge costing 1. The paper's metric is the **LCS/indel**
//! metric, in which there is no substitution primitive at all — a substitution
//! is a delete plus an insert, costing 2. That single difference is the whole
//! reason this module exists rather than reusing `AlignmentDag`, and it is also
//! the reason the rule behaves as it does on inversions (see below).
//!
//! # Why this cannot describe an inversion, measured rather than argued
//!
//! A whole-block inversion is the alignment that matches **nothing** — reference
//! `X` replaced wholesale by `revcomp(X)`. Under the simple edit distance that
//! alignment costs `|X| + |revcomp(X)| = 2n`, and
//! `d(X, Y) = |X| + |Y| - 2·LCS(X, Y)`. So the whole-block replacement is a
//! *minimal* alignment — and therefore present in the cLCS-graph at all — **if
//! and only if `LCS(X, revcomp(X)) == 0`**.
//!
//! For DNA it never is. Measured on real GRCh37 chr13 spans, `LCS(X,
//! revcomp(X)) / n` sits at 0.58–0.69 across spans 8–1024 — the Chvátal–Sankoff
//! constant for a 4-letter alphabet (γ₄ ≈ 0.654), i.e. exactly what two
//! *unrelated* strings score. `NC_000013.10:g.100809575_100810031inv`, a perfect
//! 457-base reverse complement, has `LCS = 294` and minimal distance **326**
//! against **914** for the single replacement.
//!
//! The consequence is structural, not a corner case: **the alignment that spells
//! the inversion is not in the graph**, so no rule defined over the graph can
//! select it. §3.4's inversion clause — "when the shortest path to describe an
//! inversion is a single replacement, we attempt to see if that replacement can
//! also be expressed as an inversion" — therefore cannot fire for any inversion
//! longer than a handful of bases. The paper concedes the failure mode in the
//! same paragraph ("fails when the inversion is better expressed as an allele");
//! what the measurement adds is that it is the *generic* case for DNA, not an
//! edge case.
//!
//! # The extension, stated rather than smuggled
//!
//! [`Rule::ClcsInv`] adds inversion to the algebra, which is the only thing that
//! can fix the above. An **inversion edge** is a *relabeling* of an already
//! minimal sub-path: where `ref[a..a+L]` is the reverse complement of
//! `alt[b..b+L]`, the region may be described by one labeled replacement instead
//! of by whatever the character-level alignment does inside it. It changes no
//! distance — it changes only what counts as *one* replacement — so it
//! introduces no weight to tune. The single constant, `MIN_INV_LEN = 2`, is read
//! off `DNA/inversion.md:5` ("**more than one nucleotide**"), not fitted.
//!
//! The candidate is the **longest reverse-complement common substring** of the
//! two blocks, found in `Θ(n·m)` by [`longest_revcomp_match`]. The block is then
//! costed both ways under the paper's own metric — `1 + members(prefix) +
//! members(suffix)` against `members(whole)` — and the cheaper wins, with ties
//! going to the un-inverted form. No threshold, no density heuristic.
//!
//! **What is deliberately NOT built:** the general case of *many* inversion
//! edges. Admitting every reverse-complement sub-block requires, for each cell,
//! a scan over the co-ordinate the reversal couples, which is `Θ(n·m²)` — 3.8e8
//! cell-visits on the 723-base block of `NG_009874.2` alone. One edge (the
//! longest) reaches both the whole-block inversion and the flanked-inversion
//! shape, which is what the two known failures are; nested inversions are out of
//! reach and are reported as such.

use std::collections::BTreeSet;

/// Which of the two rules to apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Rule {
    /// The published rule, faithfully: no inversion in the algebra.
    Clcs,
    /// The published rule plus one inversion edge — see the module doc.
    ClcsInv,
}

/// Shortest reverse complement admitted as an inversion.
///
/// `DNA/inversion.md:5` defines an inversion as "a sequence change where,
/// compared to a reference sequence, **more than one nucleotide** replacing the
/// original sequence is the reverse complement of the original sequence". Two,
/// therefore — read off the clause, not fitted to anything.
const MIN_INV_LEN: usize = 2;

/// Largest block this rule will attempt.
///
/// The rule is `Θ(n·m)` in both time and space, so a block of `n` reference and
/// `m` alternate bases costs one `u32` grid plus two path-count grids. This cap
/// matches `merge::MAX_SPLIT_BLOCK` so the arm declines exactly where the live
/// rule already stops splitting, rather than inventing a second size regime.
const MAX_BLOCK: usize = 1024;

/// One member of the partition: a half-open reference span and the half-open
/// alternate span replacing it, both relative to the start of the trimmed block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ClcsBlock {
    pub(crate) ref_start: u32,
    pub(crate) ref_end: u32,
    pub(crate) alt_start: u32,
    pub(crate) alt_end: u32,
}

/// Two moduli, so a path-count collision needs both to fail at once.
const MOD_A: u64 = 1_000_000_007;
const MOD_B: u64 = 998_244_353;

/// A path count, tracked modulo two primes.
type Count = (u64, u64);

const ZERO: Count = (0, 0);
const ONE: Count = (1, 1);

fn add(x: Count, y: Count) -> Count {
    ((x.0 + y.0) % MOD_A, (x.1 + y.1) % MOD_B)
}

fn mul(x: Count, y: Count) -> Count {
    ((x.0 * y.0) % MOD_A, (x.1 * y.1) % MOD_B)
}

/// IUPAC-blind complement, matching `merge::reverse_complement_bytes`'s
/// treatment of the four unambiguous bases; anything else complements to itself
/// and so can never satisfy an inversion test by accident.
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Partition `ref_block` -> `alt_block` by the cLCS rule.
///
/// Both blocks must have had their common flanks trimmed, for the same reason
/// [`super::align::AlignmentDag::build`] states.
///
/// Returns `None` when the block exceeds [`MAX_BLOCK`], so the caller can fall
/// back rather than pay a quadratic cost on a megabase block.
pub(crate) fn partition_clcs(
    ref_block: &[u8],
    alt_block: &[u8],
    rule: Rule,
) -> Option<Vec<ClcsBlock>> {
    if ref_block.len() > MAX_BLOCK || alt_block.len() > MAX_BLOCK {
        return None;
    }
    if ref_block.is_empty() && alt_block.is_empty() {
        return Some(Vec::new());
    }

    let plain = partition_faithful(ref_block, alt_block);
    if rule == Rule::Clcs {
        return Some(plain);
    }

    // --- the extension: cost the block with one inversion edge, and compare ---
    let Some((ref_at, alt_at, len)) = longest_revcomp_match(ref_block, alt_block) else {
        return Some(plain);
    };
    if len < MIN_INV_LEN {
        return Some(plain);
    }

    let prefix = partition_faithful(&ref_block[..ref_at], &alt_block[..alt_at]);
    let suffix = partition_faithful(&ref_block[ref_at + len..], &alt_block[alt_at + len..]);

    // The paper's own metric: one labeled replacement per member.  The inversion
    // edge is one replacement; the flanks cost whatever they cost.
    let with_inversion = 1 + prefix.len() + suffix.len();
    if with_inversion >= plain.len() {
        // Ties go to the un-inverted form.  That is not a tuned preference: it
        // keeps the rule from *introducing* an inversion where the character
        // alignment already describes the block just as economically.
        return Some(plain);
    }

    let mut out = Vec::with_capacity(with_inversion);
    out.extend(prefix.iter().copied());
    out.push(ClcsBlock {
        ref_start: ref_at as u32,
        ref_end: (ref_at + len) as u32,
        alt_start: alt_at as u32,
        alt_end: (alt_at + len) as u32,
    });
    let ref_off = (ref_at + len) as u32;
    let alt_off = (alt_at + len) as u32;
    out.extend(suffix.iter().map(|b| ClcsBlock {
        ref_start: b.ref_start + ref_off,
        ref_end: b.ref_end + ref_off,
        alt_start: b.alt_start + alt_off,
        alt_end: b.alt_end + alt_off,
    }));
    Some(out)
}

/// The published rule with nothing added: stages 1–3 of the module doc.
fn partition_faithful(ref_block: &[u8], alt_block: &[u8]) -> Vec<ClcsBlock> {
    let n = ref_block.len();
    let m = alt_block.len();
    if n == 0 && m == 0 {
        return Vec::new();
    }
    if n == 0 || m == 0 {
        // A pure insertion or a pure deletion is one replacement either way.
        return vec![ClcsBlock {
            ref_start: 0,
            ref_end: n as u32,
            alt_start: 0,
            alt_end: m as u32,
        }];
    }

    let width = m + 1;
    let idx = |i: usize, j: usize| i * width + j;

    // --- stage 1: the graph of all minimal alignments (simple edit distance) ---
    let prefix = indel_grid(ref_block, alt_block);
    let ref_rev: Vec<u8> = ref_block.iter().rev().copied().collect();
    let alt_rev: Vec<u8> = alt_block.iter().rev().copied().collect();
    let suffix_grid = indel_grid(&ref_rev, &alt_rev);
    let suffix = |i: usize, j: usize| suffix_grid[(n - i) * width + (m - j)];
    let total = prefix[idx(n, m)];

    let mut on_path = vec![false; (n + 1) * width];
    for i in 0..=n {
        for j in 0..=m {
            on_path[idx(i, j)] = prefix[idx(i, j)] + suffix(i, j) == total;
        }
    }

    // Admitted out-edges.  There is deliberately no diagonal on a mismatch: the
    // simple edit distance has no substitution primitive.
    let out_edges = |i: usize, j: usize| -> [Option<(usize, usize, bool)>; 3] {
        let here = prefix[idx(i, j)];
        let mut e: [Option<(usize, usize, bool)>; 3] = [None, None, None];
        if i < n && j < m && ref_block[i] == alt_block[j] {
            let next = idx(i + 1, j + 1);
            if here == prefix[next] && on_path[next] {
                e[0] = Some((i + 1, j + 1, true)); // match: a λ edge, weight 0
            }
        }
        if i < n {
            let next = idx(i + 1, j);
            if here + 1 == prefix[next] && on_path[next] {
                e[1] = Some((i + 1, j, false)); // delete
            }
        }
        if j < m {
            let next = idx(i, j + 1);
            if here + 1 == prefix[next] && on_path[next] {
                e[2] = Some((i, j + 1, false)); // insert
            }
        }
        e
    };

    // --- stage 2: shortest path, weight 1 per labeled replacement -------------
    // A maximal run of consecutive non-match steps is ONE labeled edge, so the
    // run state is part of the search state.
    const CLOSED: usize = 0;
    const OPEN: usize = 1;
    const INF: u32 = u32::MAX;
    let size = (n + 1) * width;

    // `back[s][c]` = fewest labeled edges from cell c in run state s to the sink.
    let mut back = vec![[INF, INF]; size];
    back[idx(n, m)] = [0, 0];
    for i in (0..=n).rev() {
        for j in (0..=m).rev() {
            if i == n && j == m {
                continue;
            }
            if !on_path[idx(i, j)] {
                continue;
            }
            let e = out_edges(i, j);
            for state in [CLOSED, OPEN] {
                let mut best = INF;
                for edge in e.iter().flatten() {
                    let (ni, nj, is_match) = *edge;
                    let (cost, next_state) = if is_match {
                        (0u32, CLOSED)
                    } else {
                        (u32::from(state == CLOSED), OPEN)
                    };
                    let v = back[idx(ni, nj)][next_state];
                    if v != INF {
                        best = best.min(v + cost);
                    }
                }
                back[idx(i, j)][state] = best;
            }
        }
    }
    let opt = back[idx(0, 0)][CLOSED];
    if opt == INF {
        // Cannot happen for a well-formed grid, but declining beats a panic.
        return vec![ClcsBlock {
            ref_start: 0,
            ref_end: n as u32,
            alt_start: 0,
            alt_end: m as u32,
        }];
    }

    // Forward distances and path counts over the shortest-path subgraph.
    let mut fwd = vec![[INF, INF]; size];
    let mut fcnt = vec![[ZERO, ZERO]; size];
    fwd[idx(0, 0)][CLOSED] = 0;
    fcnt[idx(0, 0)][CLOSED] = ONE;
    for i in 0..=n {
        for j in 0..=m {
            if !on_path[idx(i, j)] {
                continue;
            }
            let e = out_edges(i, j);
            for state in [CLOSED, OPEN] {
                let base = fwd[idx(i, j)][state];
                if base == INF {
                    continue;
                }
                let c = fcnt[idx(i, j)][state];
                for edge in e.iter().flatten() {
                    let (ni, nj, is_match) = *edge;
                    let (cost, next_state) = if is_match {
                        (0u32, CLOSED)
                    } else {
                        (u32::from(state == CLOSED), OPEN)
                    };
                    let nc = base + cost;
                    let slot = &mut fwd[idx(ni, nj)][next_state];
                    if nc < *slot {
                        *slot = nc;
                        fcnt[idx(ni, nj)][next_state] = c;
                    } else if nc == *slot {
                        let prev = fcnt[idx(ni, nj)][next_state];
                        fcnt[idx(ni, nj)][next_state] = add(prev, c);
                    }
                }
            }
        }
    }

    // Backward path counts, over optimal completions only.
    let mut bcnt = vec![[ZERO, ZERO]; size];
    bcnt[idx(n, m)] = [ONE, ONE];
    for i in (0..=n).rev() {
        for j in (0..=m).rev() {
            if i == n && j == m {
                continue;
            }
            if !on_path[idx(i, j)] {
                continue;
            }
            let e = out_edges(i, j);
            for state in [CLOSED, OPEN] {
                let here = back[idx(i, j)][state];
                if here == INF {
                    continue;
                }
                let mut acc = ZERO;
                for edge in e.iter().flatten() {
                    let (ni, nj, is_match) = *edge;
                    let (cost, next_state) = if is_match {
                        (0u32, CLOSED)
                    } else {
                        (u32::from(state == CLOSED), OPEN)
                    };
                    let v = back[idx(ni, nj)][next_state];
                    if v != INF && v + cost == here {
                        acc = add(acc, bcnt[idx(ni, nj)][next_state]);
                    }
                }
                bcnt[idx(i, j)][state] = acc;
            }
        }
    }

    let mut total_paths = ZERO;
    for state in [CLOSED, OPEN] {
        if fwd[idx(n, m)][state] == opt {
            total_paths = add(total_paths, fcnt[idx(n, m)][state]);
        }
    }

    // --- stage 3: local supremal — the nodes on EVERY shortest path -----------
    let mut cuts: BTreeSet<(usize, usize)> = BTreeSet::new();
    for i in 0..=n {
        for j in 0..=m {
            if !on_path[idx(i, j)] {
                continue;
            }
            let mut through = ZERO;
            for state in [CLOSED, OPEN] {
                let f = fwd[idx(i, j)][state];
                let b = back[idx(i, j)][state];
                if f == INF || b == INF || f + b != opt {
                    continue;
                }
                through = add(through, mul(fcnt[idx(i, j)][state], bcnt[idx(i, j)][state]));
            }
            if through == total_paths && through != ZERO {
                cuts.insert((i, j));
            }
        }
    }
    cuts.insert((0, 0));
    cuts.insert((n, m));

    // Is the match step OUT of `(i, j)` on a shortest labeled path?  A match is a
    // weight-0 edge to `(i+1, j+1)` that resets the run state to CLOSED, so it is
    // optimal iff `fwd(here) + 0 + back(next, CLOSED) == opt` for some incoming
    // run state.  This is the one signal that distinguishes a run boundary from a
    // cell interior to a run.
    let has_opt_match_out = |i: usize, j: usize| -> bool {
        if i >= n || j >= m || ref_block[i] != alt_block[j] {
            return false;
        }
        let head = idx(i + 1, j + 1);
        if !on_path[head] {
            return false;
        }
        let b = back[head][CLOSED];
        if b == INF {
            return false;
        }
        [CLOSED, OPEN]
            .into_iter()
            .any(|state| fwd[idx(i, j)][state] != INF && fwd[idx(i, j)][state] + b == opt)
    };

    // Drop post-dominators that sit STRICTLY INSIDE a labeled run.  In the paper's
    // compressed graph a maximal run of non-match steps is one edge with no
    // interior nodes, so such a cell is not a partition boundary; keeping it splits
    // one forced insertion/deletion run into several adjacent members — the stage-3
    // over-split.  A genuine boundary always has an optimal match step incident, on
    // its out-side (a match starting the next labeled edge) or its in-side (a match
    // ending the previous one); the source and sink are boundaries by fiat.
    // Dropping an interior cut only merges the two non-match windows either side of
    // it, so the partition stays disjoint, ordered and covering.
    let ordered: Vec<(usize, usize)> = cuts
        .into_iter()
        .filter(|&(i, j)| {
            (i, j) == (0, 0)
                || (i, j) == (n, m)
                || has_opt_match_out(i, j)
                || (i > 0 && j > 0 && has_opt_match_out(i - 1, j - 1))
        })
        .collect();
    let mut members = Vec::new();
    for pair in ordered.windows(2) {
        let (i0, j0) = pair[0];
        let (i1, j1) = pair[1];
        if i1 - i0 == 1 && j1 - j0 == 1 && ref_block[i0] == alt_block[j0] {
            continue;
        }
        members.push(ClcsBlock {
            ref_start: i0 as u32,
            ref_end: i1 as u32,
            alt_start: j0 as u32,
            alt_end: j1 as u32,
        });
    }
    members
}

/// Edit-distance grid under the **simple edit distance**: insertions and
/// deletions only, each costing 1, with no substitution primitive.
///
/// This is the one line that separates the paper's metric from ferro's
/// [`super::align::AlignmentDag`], which charges 1 for a substitution.
fn indel_grid(reference: &[u8], alt: &[u8]) -> Vec<u32> {
    let n = reference.len();
    let m = alt.len();
    let width = m + 1;
    let mut grid = vec![0u32; (n + 1) * width];
    for (j, slot) in grid[..width].iter_mut().enumerate() {
        *slot = j as u32;
    }
    for i in 1..=n {
        grid[i * width] = i as u32;
        for j in 1..=m {
            grid[i * width + j] = if reference[i - 1] == alt[j - 1] {
                grid[(i - 1) * width + (j - 1)]
            } else {
                1 + grid[(i - 1) * width + j].min(grid[i * width + (j - 1)])
            };
        }
    }
    grid
}

/// The longest substring of `reference` that is the reverse complement of a
/// substring of `alt`.
///
/// Returns `(ref_offset, alt_offset, length)`, or `None` when no two bases pair.
///
/// Works by reverse-complementing `alt` once and running the standard
/// longest-common-substring DP, which is `Θ(n·m)` time and `Θ(m)` space. The
/// co-ordinate flip is the whole trick: `ref[a..a+L] == revcomp(alt[b..b+L])` is
/// the same statement as `ref[a..a+L] == rc_alt[m-b-L..m-b]`, so a common
/// substring of `reference` and `rc_alt` at `(a, u)` of length `L` is an
/// inversion edge at `b = m - u - L`.
fn longest_revcomp_match(reference: &[u8], alt: &[u8]) -> Option<(usize, usize, usize)> {
    let n = reference.len();
    let m = alt.len();
    if n == 0 || m == 0 {
        return None;
    }
    let rc_alt: Vec<u8> = alt.iter().rev().map(|b| complement(*b)).collect();

    let mut prev = vec![0u32; m + 1];
    let mut cur = vec![0u32; m + 1];
    let mut best = (0usize, 0usize, 0usize);
    for i in 1..=n {
        for u in 1..=m {
            cur[u] = if reference[i - 1].eq_ignore_ascii_case(&rc_alt[u - 1]) {
                prev[u - 1] + 1
            } else {
                0
            };
            let len = cur[u] as usize;
            if len > best.2 {
                // The run ends at ref[i-1] and rc_alt[u-1], so it starts at
                // ref[i-len] and rc_alt[u-len].
                best = (i - len, u - len, len);
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    if best.2 == 0 {
        return None;
    }
    let (ref_at, rc_at, len) = best;
    Some((ref_at, m - rc_at - len, len))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|b| complement(*b)).collect()
    }

    /// The measurement the module doc rests on, as an executable check: the
    /// whole-block replacement is minimal only when the two blocks share no
    /// common subsequence at all, and a reverse complement always shares one.
    #[test]
    fn a_whole_block_inversion_is_not_a_minimal_alignment() {
        let reference = b"ACGTGCATGCATGCAA";
        let alt = revcomp(reference);
        let grid = indel_grid(reference, &alt);
        let distance = grid[reference.len() * (alt.len() + 1) + alt.len()];
        let single_replacement = (reference.len() + alt.len()) as u32;
        assert!(
            distance < single_replacement,
            "the single replacement ({single_replacement}) must not be minimal ({distance}); \
             if it were, the published rule could describe the inversion"
        );
    }

    /// And therefore the published rule shreds it.  This is the negative result
    /// the prototype exists to establish, pinned so it cannot be quietly lost.
    #[test]
    fn the_published_rule_shreds_a_clean_inversion() {
        let reference = b"ACGTGCATGCATGCAA";
        let alt = revcomp(reference);
        let members = partition_clcs(reference, &alt, Rule::Clcs).expect("within MAX_BLOCK");
        assert!(
            members.len() > 1,
            "expected the published rule to shred; got {members:?}"
        );
    }

    /// The extension recovers it — one member spanning the whole block.
    #[test]
    fn the_inversion_edge_recovers_a_clean_inversion() {
        let reference = b"ACGTGCATGCATGCAA";
        let alt = revcomp(reference);
        let members = partition_clcs(reference, &alt, Rule::ClcsInv).expect("within MAX_BLOCK");
        assert_eq!(members.len(), 1, "got {members:?}");
        assert_eq!(members[0].ref_start, 0);
        assert_eq!(members[0].ref_end, reference.len() as u32);
    }

    /// The flanked shape: a reverse-complement core with plain deletions either
    /// side.  One inversion edge is enough to find it.
    #[test]
    fn the_inversion_edge_finds_a_flanked_core() {
        let core = b"ACGTGCATGCATGCAAGGTTCA";
        let mut reference = b"TTTTTT".to_vec();
        reference.extend_from_slice(core);
        reference.extend_from_slice(b"GGGG");
        let alt = revcomp(core);
        let members = partition_clcs(&reference, &alt, Rule::ClcsInv).expect("within MAX_BLOCK");
        assert!(
            members.iter().any(|b| {
                let r = &reference[b.ref_start as usize..b.ref_end as usize];
                let a = &alt[b.alt_start as usize..b.alt_end as usize];
                r.len() > 1 && r == revcomp(a).as_slice()
            }),
            "no member is a reverse complement: {members:?}"
        );
    }

    /// `longest_revcomp_match`'s co-ordinate flip is the easiest thing here to
    /// get subtly wrong, so it is checked against the definition directly.
    #[test]
    fn longest_revcomp_match_returns_a_real_reverse_complement() {
        let reference = b"AAAACGTGCATTTT";
        let alt = b"GGGGATGCACGTCC";
        let (r, a, len) = longest_revcomp_match(reference, alt).expect("a match exists");
        assert!(len >= 2, "len {len}");
        assert_eq!(
            &reference[r..r + len],
            revcomp(&alt[a..a + len]).as_slice(),
            "reported span is not a reverse complement"
        );
    }

    /// Ferro's committed `AAGCTA -> TAGCTT` guard, which `general.md:56` ranks as
    /// two substitutions rather than one inversion.  The published rule agrees
    /// with the guard; the extension does NOT, and that collision is recorded
    /// here rather than re-blessed.
    #[test]
    fn the_1040_control_splits_under_the_published_rule_and_merges_under_the_extension() {
        let reference = b"AAGCTA";
        let alt = b"TAGCTT";
        let published = partition_clcs(reference, alt, Rule::Clcs).expect("within MAX_BLOCK");
        assert_eq!(
            published.len(),
            2,
            "the published rule should keep the two substitutions apart: {published:?}"
        );
        let extended = partition_clcs(reference, alt, Rule::ClcsInv).expect("within MAX_BLOCK");
        assert_eq!(
            extended.len(),
            1,
            "KNOWN COLLISION: the inversion edge costs 1 and the two substitutions cost 2, so \
             the extension merges what `inversion-vs-two-delins-76-83` and #1040 require split. \
             This is reported as a blocker for the extension, not re-blessed: {extended:?}"
        );
    }

    /// A forced insertion run is ONE labeled edge, not one member per base.  This
    /// is the stage-3 over-split: every interior cell of the run `AC -> AGGC` is a
    /// post-dominator, so the uncompressed grid used to cut it into `[insG; insG]`.
    /// The compressed-graph post-dominators — the cells with a match step incident
    /// — collapse it to a single `insGG`.
    #[test]
    fn a_forced_insertion_run_is_one_member_not_one_per_base() {
        let reference = b"AC";
        let alt = b"AGGC";
        let members = partition_clcs(reference, alt, Rule::Clcs).expect("within MAX_BLOCK");
        assert_eq!(
            members.len(),
            1,
            "the forced run of two insertions must be one member: {members:?}"
        );
        let only = members[0];
        assert_eq!(
            only.ref_start, only.ref_end,
            "a pure insertion consumes no reference"
        );
        assert_eq!(only.ref_start, 1);
        assert_eq!(
            (only.alt_start, only.alt_end),
            (1, 3),
            "the inserted span is the whole GG run"
        );
    }

    /// The same, one length up and as a deletion-bearing mix, so the fix is not
    /// pinned to the exact two-base shape: `ATG -> AGGGTG` inserts a forced `GGG`
    /// run between two matched walls and must stay one member.
    #[test]
    fn a_longer_forced_run_still_collapses() {
        let reference = b"ATG";
        let alt = b"AGGGTG";
        let members = partition_clcs(reference, alt, Rule::Clcs).expect("within MAX_BLOCK");
        let insertions: Vec<_> = members
            .iter()
            .filter(|b| b.ref_start == b.ref_end && b.alt_end > b.alt_start)
            .collect();
        assert_eq!(
            insertions.len(),
            1,
            "the forced GGG run must not be split into separate insertions: {members:?}"
        );
    }

    /// Members must partition the block: ascending, disjoint, and covering every
    /// changed column.  Checked exhaustively over a small alphabet so a coverage
    /// gap cannot hide behind a hand-picked case.
    #[test]
    fn members_partition_the_block_exhaustively() {
        let alphabet = *b"ACG";
        let mut checked = 0usize;
        for rn in 0..=3usize {
            for an in 0..=3usize {
                for rcode in 0..alphabet.len().pow(rn as u32) {
                    let reference: Vec<u8> = (0..rn)
                        .map(|k| alphabet[(rcode / alphabet.len().pow(k as u32)) % alphabet.len()])
                        .collect();
                    for acode in 0..alphabet.len().pow(an as u32) {
                        let alt: Vec<u8> = (0..an)
                            .map(|k| {
                                alphabet[(acode / alphabet.len().pow(k as u32)) % alphabet.len()]
                            })
                            .collect();
                        for rule in [Rule::Clcs, Rule::ClcsInv] {
                            let members = partition_clcs(&reference, &alt, rule).expect("small");
                            let mut last_ref = 0u32;
                            let mut last_alt = 0u32;
                            for b in &members {
                                assert!(b.ref_start >= last_ref, "{reference:?} {alt:?} {b:?}");
                                assert!(b.alt_start >= last_alt, "{reference:?} {alt:?} {b:?}");
                                assert!(b.ref_end >= b.ref_start);
                                assert!(b.alt_end >= b.alt_start);
                                assert!(b.ref_end <= reference.len() as u32);
                                assert!(b.alt_end <= alt.len() as u32);
                                // The uncovered gap before this member must be an
                                // unchanged match: equal length on both axes, and
                                // byte-identical column by column.
                                let (ref_gap, alt_gap) =
                                    (b.ref_start - last_ref, b.alt_start - last_alt);
                                assert_eq!(
                                    ref_gap, alt_gap,
                                    "uncovered gap must advance both axes equally: \
                                     {reference:?} {alt:?} {members:?}"
                                );
                                assert_eq!(
                                    &reference[last_ref as usize..b.ref_start as usize],
                                    &alt[last_alt as usize..b.alt_start as usize],
                                    "uncovered gap is not an unchanged match: \
                                     {reference:?} {alt:?} {members:?}"
                                );
                                last_ref = b.ref_end;
                                last_alt = b.alt_end;
                            }
                            assert!(last_ref <= reference.len() as u32);
                            assert!(last_alt <= alt.len() as u32);
                            // And so must the uncovered tail after the last member.
                            assert_eq!(
                                &reference[last_ref as usize..],
                                &alt[last_alt as usize..],
                                "uncovered tail is not an unchanged match: \
                                 {reference:?} {alt:?} {members:?}"
                            );
                            checked += 1;
                        }
                    }
                }
            }
        }
        assert!(
            checked > 500,
            "denominator too small to mean anything: {checked}"
        );
    }
}
