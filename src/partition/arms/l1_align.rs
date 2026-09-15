//! Alignment-based L1 cutters that complete the locked L1 match-definition set
//! (design §4.2): single-min-edit (#2), affine-gap (#3), position-wise (#4), and
//! union/LCS-maximal (#5). Each aligns `(reference, resulting)` under a different
//! cost model and cuts on the columns the alignment leaves matched; the walls are
//! the unchanged bases and the segments are the maximal changed runs between them.
//!
//! Every cutter first trims the shared flanks and declines (falls back to one
//! whole-span segment) on a changed block wider than `MAX_BLOCK`, so the O(n·m)
//! DPs stay bounded on the large genomic windows. Tracebacks use a fixed,
//! deterministic tie-break — the property the confluence gate needs — biased to
//! defer matches (a 3′-leaning placement of indels).

use crate::partition::arm::L1Segmenter;
use crate::partition::block_ctx::{Molecule, Provenance};
use crate::partition::output::Segment;

/// Largest changed block these O(n·m) rules will align; matches the cLCS cap.
const MAX_BLOCK: usize = 1024;

/// One aligned column.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Op {
    /// ref[i] == alt[j], both consumed — an unchanged wall.
    Match,
    /// ref[i] != alt[j], both consumed.
    Sub,
    /// ref[i] consumed, no alt (deletion).
    Del,
    /// alt[j] consumed, no ref (insertion).
    Ins,
}

fn shared_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}
fn shared_suffix(a: &[u8], b: &[u8], prefix: usize) -> usize {
    let max = a.len().min(b.len()) - prefix;
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take(max)
        .take_while(|(x, y)| x == y)
        .count()
}

/// Whole-span fallback segment (used when a block exceeds `MAX_BLOCK`).
fn whole_span(reference: &[u8], resulting: &[u8]) -> Vec<Segment> {
    let p = shared_prefix(reference, resulting);
    let s = shared_suffix(reference, resulting, p);
    vec![Segment {
        ref_start: p,
        ref_end: reference.len() - s,
        res_start: p,
        res_end: resulting.len() - s,
    }]
}

/// Turn an op sequence over the trimmed block into absolute segments — the
/// maximal runs of non-`Match` ops, shifted back by the prefix length.
fn ops_to_segments(ops: &[Op], prefix: usize) -> Vec<Segment> {
    let mut segs = Vec::new();
    let (mut ri, mut ai) = (0usize, 0usize);
    let mut open: Option<(usize, usize)> = None; // (ref_start, res_start) of the open run
    let close =
        |open: &mut Option<(usize, usize)>, ri: usize, ai: usize, segs: &mut Vec<Segment>| {
            if let Some((rs, as_)) = open.take() {
                segs.push(Segment {
                    ref_start: prefix + rs,
                    ref_end: prefix + ri,
                    res_start: prefix + as_,
                    res_end: prefix + ai,
                });
            }
        };
    for op in ops {
        match op {
            Op::Match => {
                close(&mut open, ri, ai, &mut segs);
                ri += 1;
                ai += 1;
            }
            Op::Sub => {
                open.get_or_insert((ri, ai));
                ri += 1;
                ai += 1;
            }
            Op::Del => {
                open.get_or_insert((ri, ai));
                ri += 1;
            }
            Op::Ins => {
                open.get_or_insert((ri, ai));
                ai += 1;
            }
        }
    }
    close(&mut open, ri, ai, &mut segs);
    segs
}

/// Run `align` over the trimmed, size-capped block, else fall back to one span.
fn segment_with(
    reference: &[u8],
    resulting: &[u8],
    align: impl Fn(&[u8], &[u8]) -> Vec<Op>,
) -> Vec<Segment> {
    let p = shared_prefix(reference, resulting);
    let s = shared_suffix(reference, resulting, p);
    let ref_block = &reference[p..reference.len() - s];
    let alt_block = &resulting[p..resulting.len() - s];
    if ref_block.len() > MAX_BLOCK || alt_block.len() > MAX_BLOCK {
        return whole_span(reference, resulting);
    }
    let ops = align(ref_block, alt_block);
    ops_to_segments(&ops, p)
}

// --- #2 single min-edit path (unit-cost Levenshtein), one deterministic path ---

fn levenshtein_ops(r: &[u8], a: &[u8]) -> Vec<Op> {
    let (n, m) = (r.len(), a.len());
    let mut d = vec![vec![0u32; m + 1]; n + 1];
    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i as u32;
    }
    for (j, cell) in d[0].iter_mut().enumerate() {
        *cell = j as u32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let sub = d[i - 1][j - 1] + if r[i - 1] == a[j - 1] { 0 } else { 1 };
            d[i][j] = sub.min(d[i - 1][j] + 1).min(d[i][j - 1] + 1);
        }
    }
    // Traceback, deferring the diagonal so indels sit toward the 3' end.
    let (mut i, mut j) = (n, m);
    let mut ops = Vec::new();
    while i > 0 || j > 0 {
        let here = d[i][j];
        if i > 0 && here == d[i - 1][j] + 1 {
            ops.push(Op::Del);
            i -= 1;
        } else if j > 0 && here == d[i][j - 1] + 1 {
            ops.push(Op::Ins);
            j -= 1;
        } else {
            ops.push(if r[i - 1] == a[j - 1] {
                Op::Match
            } else {
                Op::Sub
            });
            i -= 1;
            j -= 1;
        }
    }
    ops.reverse();
    ops
}

/// One ATOMIC segment per non-`Match` op — the finest sound segmentation: every
/// substituted column, every deleted base, every inserted base is its own
/// one-wide segment. Used by `MaximalSplitL1`.
fn atomic_segments(ops: &[Op], prefix: usize) -> Vec<Segment> {
    let mut segs = Vec::new();
    let (mut ri, mut ai) = (0usize, 0usize);
    for op in ops {
        match op {
            Op::Match => {
                ri += 1;
                ai += 1;
            }
            Op::Sub => {
                segs.push(Segment {
                    ref_start: prefix + ri,
                    ref_end: prefix + ri + 1,
                    res_start: prefix + ai,
                    res_end: prefix + ai + 1,
                });
                ri += 1;
                ai += 1;
            }
            Op::Del => {
                segs.push(Segment {
                    ref_start: prefix + ri,
                    ref_end: prefix + ri + 1,
                    res_start: prefix + ai,
                    res_end: prefix + ai,
                });
                ri += 1;
            }
            Op::Ins => {
                segs.push(Segment {
                    ref_start: prefix + ri,
                    ref_end: prefix + ri,
                    res_start: prefix + ai,
                    res_end: prefix + ai + 1,
                });
                ai += 1;
            }
        }
    }
    segs
}

/// L1: the MAXIMAL (finest) sound split — one atomic segment per changed column /
/// indel base, from a single min-edit alignment. This is the extreme the
/// merge-only architecture starts from: every legal coarser description is a
/// coarsening of it, so a canonicaliser can reach any target by merging alone,
/// never splitting. (Design experiment for the max-split + merge-only idea.)
pub struct MaximalSplitL1;
impl L1Segmenter for MaximalSplitL1 {
    fn name(&self) -> &str {
        "maximal-split"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        _molecule: Molecule,
        _provenance: &Provenance,
    ) -> Vec<Segment> {
        let p = shared_prefix(reference, resulting);
        let s = shared_suffix(reference, resulting, p);
        let ref_block = &reference[p..reference.len() - s];
        let alt_block = &resulting[p..resulting.len() - s];
        if ref_block.len() > MAX_BLOCK || alt_block.len() > MAX_BLOCK {
            return whole_span(reference, resulting);
        }
        atomic_segments(&levenshtein_ops(ref_block, alt_block), p)
    }
}

/// L1 #2: one minimal-edit (unit-cost Levenshtein) alignment. Where #1 keeps only
/// the bases matched in *every* minimal alignment, this commits to a single path,
/// so it cuts on more walls.
pub struct SingleMinEditL1;
impl L1Segmenter for SingleMinEditL1 {
    fn name(&self) -> &str {
        "single-min-edit"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        _molecule: Molecule,
        _provenance: &Provenance,
    ) -> Vec<Segment> {
        segment_with(reference, resulting, levenshtein_ops)
    }
}

// --- #3 affine-gap (Gotoh): one gap-open plus per-base extension ---

fn affine_ops(r: &[u8], a: &[u8]) -> Vec<Op> {
    const SUB: u32 = 1;
    const OPEN: u32 = 2; // first gap base
    const EXTEND: u32 = 1; // each subsequent gap base
    let inf = u32::MAX / 4;
    let (n, m) = (r.len(), a.len());
    // m_[i][j]: best ending in match/sub; dx: gap in alt (del); dy: gap in ref (ins).
    let mut m_ = vec![vec![inf; m + 1]; n + 1];
    let mut dx = vec![vec![inf; m + 1]; n + 1];
    let mut dy = vec![vec![inf; m + 1]; n + 1];
    m_[0][0] = 0;
    for (i, row) in dx.iter_mut().enumerate().skip(1) {
        row[0] = OPEN + (i as u32 - 1) * EXTEND;
    }
    for (j, cell) in dy[0].iter_mut().enumerate().skip(1) {
        *cell = OPEN + (j as u32 - 1) * EXTEND;
    }
    for i in 1..=n {
        for j in 1..=m {
            let diag = m_[i - 1][j - 1].min(dx[i - 1][j - 1]).min(dy[i - 1][j - 1]);
            m_[i][j] = diag + if r[i - 1] == a[j - 1] { 0 } else { SUB };
            dx[i][j] = (m_[i - 1][j] + OPEN).min(dx[i - 1][j] + EXTEND);
            dy[i][j] = (m_[i][j - 1] + OPEN).min(dy[i][j - 1] + EXTEND);
        }
    }
    // Traceback across the three planes; defer the diagonal (gaps first).
    #[derive(Clone, Copy)]
    enum Plane {
        M,
        Dx,
        Dy,
    }
    let (mut i, mut j) = (n, m);
    let mut plane = {
        let (fm, fx, fy) = (m_[n][m], dx[n][m], dy[n][m]);
        if fx <= fm && fx <= fy {
            Plane::Dx
        } else if fy <= fm {
            Plane::Dy
        } else {
            Plane::M
        }
    };
    let mut ops = Vec::new();
    while i > 0 || j > 0 {
        match plane {
            Plane::Dx => {
                // `dx[0][j]` is INF and `Plane::M` redirects at `i == 0`, so a
                // deletion plane is never entered at the top boundary.
                debug_assert!(i >= 1, "Dx plane entered at i == 0");
                ops.push(Op::Del);
                let from_open = dx[i][j] == m_[i - 1][j] + OPEN;
                i -= 1;
                if from_open {
                    plane = Plane::M;
                }
            }
            Plane::Dy => {
                debug_assert!(j >= 1, "Dy plane entered at j == 0");
                ops.push(Op::Ins);
                let from_open = dy[i][j] == m_[i][j - 1] + OPEN;
                j -= 1;
                if from_open {
                    plane = Plane::M;
                }
            }
            Plane::M => {
                if i == 0 {
                    plane = Plane::Dy;
                    continue;
                }
                if j == 0 {
                    plane = Plane::Dx;
                    continue;
                }
                ops.push(if r[i - 1] == a[j - 1] {
                    Op::Match
                } else {
                    Op::Sub
                });
                let diag = m_[i - 1][j - 1].min(dx[i - 1][j - 1]).min(dy[i - 1][j - 1]);
                i -= 1;
                j -= 1;
                // Which plane the optimal diagonal predecessor came from.
                plane = if diag == dx[i][j] {
                    Plane::Dx
                } else if diag == dy[i][j] {
                    Plane::Dy
                } else {
                    Plane::M
                };
            }
        }
    }
    ops.reverse();
    ops
}

/// L1 #3: affine-gap alignment. A gap-open penalty makes it prefer one long indel
/// over several scattered ones, so it coalesces indels where the unit-cost rules
/// would spread them.
pub struct AffineGapL1;
impl L1Segmenter for AffineGapL1 {
    fn name(&self) -> &str {
        "affine-gap"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        _molecule: Molecule,
        _provenance: &Provenance,
    ) -> Vec<Segment> {
        segment_with(reference, resulting, affine_ops)
    }
}

// --- #5 union / LCS-maximal: no substitutions, maximize matched bases ---

fn lcs_ops(r: &[u8], a: &[u8]) -> Vec<Op> {
    let (n, m) = (r.len(), a.len());
    let mut l = vec![vec![0u32; m + 1]; n + 1];
    for i in 1..=n {
        for j in 1..=m {
            l[i][j] = if r[i - 1] == a[j - 1] {
                l[i - 1][j - 1] + 1
            } else {
                l[i - 1][j].max(l[i][j - 1])
            };
        }
    }
    let (mut i, mut j) = (n, m);
    let mut ops = Vec::new();
    while i > 0 || j > 0 {
        if i > 0 && j > 0 && r[i - 1] == a[j - 1] {
            ops.push(Op::Match);
            i -= 1;
            j -= 1;
        } else if i > 0 && (j == 0 || l[i - 1][j] >= l[i][j - 1]) {
            ops.push(Op::Del);
            i -= 1;
        } else {
            ops.push(Op::Ins);
            j -= 1;
        }
    }
    ops.reverse();
    ops
}

/// L1 #5: cut on a maximal common subsequence — the permissive extreme. It admits
/// no substitutions, so it keeps the most bases matched and fragments the most;
/// the diagnostic upper bracket on how far segmentation can split a variant.
pub struct UnionLcsL1;
impl L1Segmenter for UnionLcsL1 {
    fn name(&self) -> &str {
        "union-lcs"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        _molecule: Molecule,
        _provenance: &Provenance,
    ) -> Vec<Segment> {
        segment_with(reference, resulting, lcs_ops)
    }
}

// --- #4 position-wise: column i of ref pairs with column i of alt (no indels) ---

/// L1 #4: position-wise cutting. Pairs reference column `i` with resulting column
/// `i` — blind to indels by construction — so an unchanged column is a wall and
/// each maximal run of differing columns is one segment. Undefined when the
/// lengths differ, so it falls back to one whole-span segment there.
pub struct PositionWiseL1;
impl L1Segmenter for PositionWiseL1 {
    fn name(&self) -> &str {
        "position-wise-cut"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        _molecule: Molecule,
        _provenance: &Provenance,
    ) -> Vec<Segment> {
        if reference.len() != resulting.len() {
            return whole_span(reference, resulting);
        }
        let ops: Vec<Op> = reference
            .iter()
            .zip(resulting)
            .map(|(r, a)| if r == a { Op::Match } else { Op::Sub })
            .collect();
        ops_to_segments(&ops, 0)
    }
}
