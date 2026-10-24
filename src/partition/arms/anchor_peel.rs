//! `AnchorPeelL1` — the Fable round-1 independent design for the peel frontier.
//!
//! Thesis: the member-beside-a-member peel is an **L1 (wall) question, not a typer
//! question** — structured members (inv, dup) corrupt the plain edit metric that
//! defines walls, so they must be peeled BEFORE wall-finding, not labeled after.
//! Walls come in two classes:
//!
//! * **Positional walls** (zero relative displacement): maximal runs of
//!   position-wise identical columns in an equal-length residual. Genuine at any
//!   length (`general.md:33`; B's measured 100% on Delins/multi-Sub buckets).
//! * **Displaced walls** (identity or revcomp anchors at nonzero displacement):
//!   manufactured by gap placement unless long enough that a chance match is
//!   excluded (`unequal-length-block-a-placed-gap-is-not-a-separation` widened to
//!   a composition floor — the same house-policy shape as the inverted-dup
//!   coincidence floor).
//!
//! The recursion, per trimmed block:
//! 1. **W0** — whole-span revcomp (>=2) is ONE segment, unconditionally
//!    (`whole-span-reverse-complement-types-as-inv`; a ruling, not a chance call).
//! 2. **Dup peel** — a net insertion wholly consumed by copies of its immediate
//!    5' reference is peeled as zero-width dup segments; the equal-length residual
//!    recurses (`duplication-must-ranks-the-label-not-the-partition`: the label is
//!    applied to the derived piece; the peel is what derives the piece).
//! 3. **Anchor split** — the longest identity/revcomp substring match passing the
//!    chance floor pins a correspondence; flanks recurse independently. A revcomp
//!    anchor emits an inversion segment; an identity anchor emits nothing (it is
//!    unchanged content, possibly displaced).
//! 4. **Leaf** — equal-length: positional walls (position-wise runs — the
//!    coalesce-structural reading #2174 ratified). Unequal-length: forced-unchanged
//!    walls (`unchanged-is-read-over-every-minimal-alignment`), leaving the
//!    manufactured single-base separations to the downstream molecule/Dial-B
//!    carve-outs (W1 / C2 / C3), where the ledger already scopes them.
//!
//! Molecule-only by construction: nothing here reads the frame; the DNA
//! coincidence collapse stays in the `AxisAware` wrapper / Dial-B.

use crate::partition::arm::L1Segmenter;
use crate::partition::arms::arms::AllAlignmentSplitL1;
use crate::partition::block_ctx::{Molecule, Provenance};
use crate::partition::output::Segment;

/// Largest block the O(n·m) anchor DP will process; over it, one whole-span segment.
const MAX_BLOCK: usize = 1024;
/// Recursion guard (anchors strictly shrink the block, so this is belt-and-braces).
const MAX_DEPTH: usize = 64;

fn shared_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}
fn shared_suffix(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| x == y)
        .count()
}

fn revcomp(b: &[u8]) -> Vec<u8> {
    b.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b'u',
            b'u' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
}

/// Longest common substring of `r` and `a` as `(r_off, a_off, len)`, deterministic
/// (first maximal match in (r_off, a_off) scan order), `None` if no byte matches.
fn lcs_substring(r: &[u8], a: &[u8]) -> Option<(usize, usize, usize)> {
    if r.is_empty() || a.is_empty() {
        return None;
    }
    let mut best = (0usize, 0usize, 0usize); // (len, r_off, a_off)
    let mut prev = vec![0usize; a.len() + 1];
    let mut cur = vec![0usize; a.len() + 1];
    for i in 1..=r.len() {
        cur[0] = 0;
        for j in 1..=a.len() {
            cur[j] = if r[i - 1] == a[j - 1] {
                prev[j - 1] + 1
            } else {
                0
            };
            if cur[j] > best.0 {
                best = (cur[j], i - cur[j], j - cur[j]);
            }
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    if best.0 == 0 {
        None
    } else {
        Some((best.1, best.2, best.0))
    }
}

/// The chance floor: accept an anchor of length `n` between blocks of length
/// `rl`/`al` iff the expected count of chance n-mer matches (`positions × p^n`)
/// is at most `cap`. `p` is the per-column chance-match probability: 0.25 for the
/// uniform floor, or the window-composition estimate (the inverted-dup record's
/// composition-aware shape) for the composition variant.
fn anchor_significant(n: usize, rl: usize, al: usize, p: f64, cap: f64) -> bool {
    if n == 0 {
        return false;
    }
    let positions = ((rl - n + 1) as f64) * ((al - n + 1) as f64);
    positions * p.powi(n as i32) <= cap
}

/// Per-column chance that a reference base drawn from `rb`'s composition matches
/// an alt base drawn from `ab`'s, under `map` (identity, or complement for the
/// revcomp/inversion hypothesis). In uniform-random content this is 0.25; in a
/// low-complexity window it approaches 1, which is what makes a short "match"
/// there meaningless (`inverted-duplication-is-derived-as-ins-range-inv`'s
/// measured lesson: "in a skewed window a short match is far likelier than
/// 4^-n").
fn column_match_p(rb: &[u8], ab: &[u8], map: impl Fn(u8) -> u8) -> f64 {
    let mut fr = [0f64; 256];
    let mut fa = [0f64; 256];
    for &b in rb {
        fr[map(b) as usize] += 1.0;
    }
    for &b in ab {
        fa[b as usize] += 1.0;
    }
    let (nr, na) = (rb.len() as f64, ab.len() as f64);
    if nr == 0.0 || na == 0.0 {
        return 0.25;
    }
    let mut p = 0.0;
    for i in 0..256 {
        p += (fr[i] / nr) * (fa[i] / na);
    }
    p.clamp(0.25, 1.0)
}

/// The anchor-peel L1. `cap` is the chance floor for displaced anchors
/// (identity and revcomp); W0 and the positional walls are unconditional.
/// `composition` switches the floor's per-column probability from the uniform
/// 0.25 to the window-composition estimate.
pub struct AnchorPeelL1 {
    cap: f64,
    composition: bool,
    name: String,
}

impl AnchorPeelL1 {
    pub fn new(cap: f64, name: &str) -> Self {
        Self {
            cap,
            composition: false,
            name: name.to_string(),
        }
    }
    pub fn new_composition(cap: f64, name: &str) -> Self {
        Self {
            cap,
            composition: true,
            name: name.to_string(),
        }
    }
}

impl L1Segmenter for AnchorPeelL1 {
    fn name(&self) -> &str {
        &self.name
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        molecule: Molecule,
        provenance: &Provenance,
    ) -> Vec<Segment> {
        let mut out = Vec::new();
        peel(
            reference,
            resulting,
            0,
            reference.len(),
            0,
            resulting.len(),
            molecule,
            provenance,
            Floor {
                cap: self.cap,
                composition: self.composition,
            },
            0,
            &mut out,
        );
        out
    }
}

/// The displaced-anchor chance floor: expected-count cap plus whether the
/// per-column probability is composition-derived.
#[derive(Clone, Copy)]
struct Floor {
    cap: f64,
    composition: bool,
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b'u',
        b'u' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        other => other,
    }
}

#[allow(clippy::too_many_arguments)]
fn peel(
    reference: &[u8],
    resulting: &[u8],
    orig_r0: usize,
    orig_r1: usize,
    orig_a0: usize,
    orig_a1: usize,
    molecule: Molecule,
    provenance: &Provenance,
    floor: Floor,
    depth: usize,
    out: &mut Vec<Segment>,
) {
    let (mut r0, mut r1, mut a0, mut a1) = (orig_r0, orig_r1, orig_a0, orig_a1);
    // Trim the shared flanks of this sub-block.
    let p = shared_prefix(&reference[r0..r1], &resulting[a0..a1]);
    r0 += p;
    a0 += p;
    let s = shared_suffix(&reference[r0..r1], &resulting[a0..a1]);
    r1 -= s;
    a1 -= s;
    let (rl, al) = (r1 - r0, a1 - a0);
    if rl == 0 && al == 0 {
        return;
    }
    let rb = &reference[r0..r1];
    let ab = &resulting[a0..a1];
    let whole = Segment {
        ref_start: r0,
        ref_end: r1,
        res_start: a0,
        res_end: a1,
    };
    if rl > MAX_BLOCK || al > MAX_BLOCK || depth >= MAX_DEPTH {
        out.push(whole);
        return;
    }
    // W0: whole-span reverse complement is one (inverted) segment, unconditionally.
    if rl >= 2 && rl == al && ab == revcomp(rb).as_slice() {
        out.push(whole);
        return;
    }
    // Dup peel: a net insertion wholly consumed by copies of the immediate 5'
    // reference, leaving an equal-length (or empty) residual.
    if al > rl {
        let mut units: Vec<usize> = Vec::new();
        let mut cursor = a0;
        let mut net = al - rl;
        while net > 0 {
            let maxu = net.min(r0);
            let mut found = 0usize;
            for u in (1..=maxu).rev() {
                if reference[r0 - u..r0] == resulting[cursor..cursor + u] {
                    found = u;
                    break;
                }
            }
            if found == 0 {
                break;
            }
            units.push(found);
            cursor += found;
            net -= found;
        }
        if net == 0 && !units.is_empty() {
            let mut c = a0;
            for u in units {
                out.push(Segment {
                    ref_start: r0,
                    ref_end: r0,
                    res_start: c,
                    res_end: c + u,
                });
                c += u;
            }
            peel(
                reference,
                resulting,
                r0,
                r1,
                c,
                a1,
                molecule,
                provenance,
                floor,
                depth + 1,
                out,
            );
            return;
        }
        // Mirror: 3'-edge dup peel. An alt-block suffix that duplicates the
        // reference immediately FOLLOWING the block 3'-shifts through it and is a
        // duplication of those bases (the trim absorbed the copy's source into
        // the suffix). Same completion condition: the net insertion must be
        // wholly consumed, leaving an equal-length residual.
        let mut units3: Vec<usize> = Vec::new();
        let mut end = a1;
        let mut net3 = al - rl;
        while net3 > 0 {
            let maxu = net3.min(reference.len() - r1).min(end - a0);
            let mut found = 0usize;
            for u in (1..=maxu).rev() {
                if reference[r1..r1 + u] == resulting[end - u..end] {
                    found = u;
                    break;
                }
            }
            if found == 0 {
                break;
            }
            units3.push(found);
            end -= found;
            net3 -= found;
        }
        // Emit the shifted single-unit form: the duplication's canonical position
        // is AFTER the reference bases it copies (it 3'-shifts through them), so
        // the segment sits at ref `r1 + u` with its payload read from the matched
        // suffix (`resulting[a1..a1+u]` == the unit, by the trim property) — which
        // is what lets the plain ladder's 5'-copy check label it `dup`.
        if net3 == 0 && units3.len() == 1 {
            let u = units3[0];
            if a1 + u <= orig_a1 {
                peel(
                    reference,
                    resulting,
                    r0,
                    r1,
                    a0,
                    end,
                    molecule,
                    provenance,
                    floor,
                    depth + 1,
                    out,
                );
                out.push(Segment {
                    ref_start: r1 + u,
                    ref_end: r1 + u,
                    res_start: a1,
                    res_end: a1 + u,
                });
                return;
            }
        }
    }
    // Anchor split: the longest identity/revcomp match passing the chance floor.
    //
    // Identity anchors are searched only in UNEQUAL-length blocks: in an
    // equal-length (net-neutral) block a displaced identity match is by
    // construction the product of a compensating indel pair — the rotation shape
    // #2174 ratifies as one spanning delins — so it is a shift artifact at any
    // length, while at zero displacement it is already a positional wall the leaf
    // finds. Revcomp anchors are searched everywhere (an inversion is a property
    // of the bases, not of the length balance). Identity preferred on a tie.
    let id_anchor = if rl != al {
        let p_id = if floor.composition {
            column_match_p(rb, ab, |b| b)
        } else {
            0.25
        };
        lcs_substring(rb, ab).filter(|&(_, _, n)| anchor_significant(n, rl, al, p_id, floor.cap))
    } else {
        None
    };
    let rc = revcomp(rb);
    let inv_anchor = lcs_substring(&rc, ab)
        .map(|(x, ai, n)| (rl - x - n, ai, n))
        .filter(|&(_, _, n)| {
            let p_inv = if floor.composition {
                column_match_p(rb, ab, comp)
            } else {
                0.25
            };
            anchor_significant(n, rl, al, p_inv, floor.cap)
        })
        // Palindrome guard: an anchor whose reference content is its own reverse
        // complement explains nothing an identity does not — "inverting" it is a
        // no-op, so labeling it `inv` mints a spurious member over unchanged
        // content (measured: Sub+Sub blocks in A/T-rich windows grew an Inv(6)
        // over an unchanged ATATAT run).
        .filter(|&(ri, _, n)| {
            let content = &rb[ri..ri + n];
            revcomp(content) != content
        });
    let best: Option<(bool, usize, usize, usize)> = match (id_anchor, inv_anchor) {
        (Some((_, _, n)), Some((iri, iai, inn))) if inn > n => Some((true, iri, iai, inn)),
        (Some((ri, ai, n)), _) => Some((false, ri, ai, n)),
        (None, Some((ri, ai, n))) => Some((true, ri, ai, n)),
        (None, None) => None,
    };
    if let Some((is_inv, ri, ai, n)) = best {
        if is_inv {
            // Extend the inversion outward while the complement holds — into the
            // trimmed flanks if need be. Greedy flank trimming can nibble an
            // inversion's edge (a flank base coinciding with the revcomp edge is
            // trimmed as identity), which mis-attributes one inverted base to the
            // flank and strands its partner in a residual; extension recovers the
            // maximal inversion span the ceiling derives.
            // Guard: a step must keep at least one of its two bases inside the
            // trimmed block — a step whose bases BOTH come from identity flanks
            // would convert unchanged content into inverted content (and mint a
            // compensating member elsewhere), which is how a chance complement
            // pair at the seam over-extends.
            let (mut x, mut y) = (r0 + ri, r0 + ri + n);
            let (mut u, mut v) = (a0 + ai, a0 + ai + n);
            while x > orig_r0
                && v < orig_a1
                && (x > r0 || v < a1)
                && resulting[v] == comp(reference[x - 1])
            {
                x -= 1;
                v += 1;
            }
            while y < orig_r1
                && u > orig_a0
                && (y < r1 || u > a0)
                && resulting[u - 1] == comp(reference[y])
            {
                y += 1;
                u -= 1;
            }
            peel(
                reference,
                resulting,
                orig_r0,
                x,
                orig_a0,
                u,
                molecule,
                provenance,
                floor,
                depth + 1,
                out,
            );
            out.push(Segment {
                ref_start: x,
                ref_end: y,
                res_start: u,
                res_end: v,
            });
            peel(
                reference,
                resulting,
                y,
                orig_r1,
                v,
                orig_a1,
                molecule,
                provenance,
                floor,
                depth + 1,
                out,
            );
        } else {
            peel(
                reference,
                resulting,
                r0,
                r0 + ri,
                a0,
                a0 + ai,
                molecule,
                provenance,
                floor,
                depth + 1,
                out,
            );
            peel(
                reference,
                resulting,
                r0 + ri + n,
                r1,
                a0 + ai + n,
                a1,
                molecule,
                provenance,
                floor,
                depth + 1,
                out,
            );
        }
        return;
    }
    // Leaf.
    if rl == al {
        // Positional walls: maximal runs of position-wise differing columns.
        let mut i = 0usize;
        while i < rl {
            if rb[i] == ab[i] {
                i += 1;
                continue;
            }
            let st = i;
            while i < rl && rb[i] != ab[i] {
                i += 1;
            }
            out.push(Segment {
                ref_start: r0 + st,
                ref_end: r0 + i,
                res_start: a0 + st,
                res_end: a0 + i,
            });
        }
    } else {
        // Forced-unchanged walls over the unequal residual. Each segment is then
        // re-trimmed (an alignment-pinned segment can carry a shared flank, e.g.
        // `A -> AT` where the ceiling's 3'-shifted form is a bare insertion), and
        // on a DNA molecule maximal runs separated by a SINGLE unchanged base are
        // merged back into one spanning segment: a lone forced base inside a net
        // indel is an alignment-placed coincidence, not a separation
        // (`unequal-length-block-a-placed-gap-is-not-a-separation`; the ceiling's
        // own region isolation walls only at runs >= 2). RNA keeps every wall
        // (no `delins.md:47` counterpart off the DNA axis).
        let mut segs: Vec<Segment> = Vec::new();
        for s in AllAlignmentSplitL1.segment(rb, ab, molecule, provenance) {
            let (mut rs, mut re) = (s.ref_start, s.ref_end);
            let (mut ss, mut se) = (s.res_start, s.res_end);
            let p2 = shared_prefix(&rb[rs..re], &ab[ss..se]);
            rs += p2;
            ss += p2;
            let s2 = shared_suffix(&rb[rs..re], &ab[ss..se]);
            re -= s2;
            se -= s2;
            if rs == re && ss == se {
                continue;
            }
            segs.push(Segment {
                ref_start: r0 + rs,
                ref_end: r0 + re,
                res_start: a0 + ss,
                res_end: a0 + se,
            });
        }
        // Gated the same way the shared `coincidence` core gates its callers — an
        // individuation claim suppresses this leaf collapse too, even on DNA. Not
        // routed through the core itself: this collapse is deliberately
        // direction-unbounded (see `collapse_licensed_gap1_runs`'s doc), a real,
        // documented divergence from the shared predicates. The suppression matches
        // the core's BY POLICY, not by shared code; if a third gate copy ever
        // appears, route all three through the core instead.
        if molecule == Molecule::Dna && !provenance.individuation.suppresses_coincidence_collapse()
        {
            out.extend(collapse_licensed_gap1_runs(segs));
        } else {
            out.extend(segs);
        }
    }
}

/// Collapse maximal runs of segments separated by exactly one unchanged base (in
/// BOTH coordinates) into one spanning segment, where licensed:
/// * some segment is gap-bearing — supplies bases while consuming a different
///   count (the `delins.md:46` coincidence mechanism; direction-unbounded here,
///   because the measured ceiling merges net-insertion `[ins;1;ins]` runs too); or
/// * the C3 shape — a net-deletion run whose segments all consume reference, at
///   most one a pure deletion, some supplying bases
///   (`unequal-length-block-a-placed-gap-is-not-a-separation`).
///
/// An unlicensed run (e.g. position-anchored subs beside a deletion) keeps its
/// walls — the unconditional version of this merge measurably destroyed
/// Sub+Sub+Del agreement.
fn collapse_licensed_gap1_runs(segs: Vec<Segment>) -> Vec<Segment> {
    if segs.len() < 2 {
        return segs;
    }
    let gap_bearing = |s: &Segment| {
        let (rl, al) = (s.ref_end - s.ref_start, s.res_end - s.res_start);
        al > 0 && rl != al
    };
    let mut out: Vec<Segment> = Vec::new();
    let mut run_start = 0usize;
    for i in 0..segs.len() {
        let run_ends = i + 1 == segs.len()
            || segs[i + 1].ref_start != segs[i].ref_end + 1
            || segs[i + 1].res_start != segs[i].res_end + 1;
        if !run_ends {
            continue;
        }
        let run = &segs[run_start..=i];
        run_start = i + 1;
        if run.len() < 2 {
            out.extend_from_slice(run);
            continue;
        }
        let ref_span = run[run.len() - 1].ref_end - run[0].ref_start;
        let res_span = run[run.len() - 1].res_end - run[0].res_start;
        let c3_shape = ref_span > res_span
            && run.iter().all(|s| s.ref_end > s.ref_start)
            && run.iter().any(|s| s.res_end > s.res_start)
            && run.iter().filter(|s| s.res_end == s.res_start).count() <= 1;
        if run.iter().any(gap_bearing) || c3_shape {
            out.push(Segment {
                ref_start: run[0].ref_start,
                ref_end: run[run.len() - 1].ref_end,
                res_start: run[0].res_start,
                res_end: run[run.len() - 1].res_end,
            });
        } else {
            out.extend_from_slice(run);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segs(reference: &[u8], resulting: &[u8], cap: f64) -> Vec<(usize, usize, usize, usize)> {
        AnchorPeelL1::new(cap, "peel")
            .segment(reference, resulting, Molecule::Dna, &Provenance::none())
            .iter()
            .map(|s| (s.ref_start, s.ref_end, s.res_start, s.res_end))
            .collect()
    }

    /// dup-ac-sub (#2175): the CA dup peels as a zero-width segment, the residual
    /// is the A>C sub — the ratified `[13_14dup;15A>C]` geometry.
    #[test]
    fn dup_peel_reaches_the_2175_dup_beside_a_sub() {
        let reference = b"ACGTTCAGGTCACAATTAGCTAGCTAG";
        let resulting = b"ACGTTCAGGTCACACACTTAGCTAGCTAG";
        let got = segs(reference, resulting, 0.5);
        // Zero-width dup at ref 14 with res span of 2 (the CA copy), then the sub.
        assert_eq!(got, vec![(14, 14, 14, 16), (14, 15, 16, 17)]);
    }

    /// rot-actg (#2174): GACT->ACTG is an equal-length rotation; the interior ACT
    /// is a displaced (chance-level) anchor and must NOT split — one segment.
    #[test]
    fn rotation_stays_one_segment() {
        let reference = b"ACGTTCAGGTGACTTTAGCTAGCTAG";
        let resulting = b"ACGTTCAGGTACTGTTAGCTAGCTAG";
        let got = segs(reference, resulting, 0.5);
        assert_eq!(got, vec![(10, 14, 10, 14)]);
    }

    /// Whole-span revcomp is one segment (W0), unconditionally.
    #[test]
    fn whole_span_revcomp_is_one_segment() {
        let got = segs(
            b"ACGTTCAGGTAAGCTATTAGCTAGCTAG",
            b"ACGTTCAGGTTAGCTTTTAGCTAGCTAG",
            0.01,
        );
        assert_eq!(got, vec![(10, 16, 10, 16)]);
    }

    /// A long inversion beside a deletion: the revcomp anchor pins the inversion,
    /// the flank leaf carries the deletion.
    #[test]
    fn inv_anchor_peels_a_deletion_neighbour() {
        // ref = P + DEL(3) + INV(10) + Q ; res = P + revcomp(INV) + Q
        let p = b"ACGTACGTAC".to_vec();
        let del = b"TTT".to_vec();
        let inv = b"AACCGGTTCA".to_vec();
        let q = b"GCATGCATGC".to_vec();
        let mut reference = p.clone();
        reference.extend(&del);
        reference.extend(&inv);
        reference.extend(&q);
        let mut resulting = p.clone();
        resulting.extend(revcomp(&inv));
        resulting.extend(&q);
        let got = segs(&reference, &resulting, 0.5);
        // A deletion segment (res-empty) and an equal-length revcomp segment.
        assert_eq!(got.len(), 2, "got {got:?}");
        assert!(got.iter().any(|&(rs, re, ss, se)| re > rs && se == ss));
        assert!(got
            .iter()
            .any(|&(rs, re, ss, se)| re - rs == 10 && se - ss == 10));
    }

    /// Two subs separated by three unchanged bases stay two segments
    /// (positional walls; `general.md:33`).
    #[test]
    fn separated_subs_stay_separate() {
        let got = segs(
            b"ACGTTCAGGTGTACTTAGCTAGCTAG",
            b"ACGTTCAGGTATACATAGCTAGCTAG",
            0.01,
        );
        assert_eq!(got, vec![(10, 11, 10, 11), (14, 15, 14, 15)]);
    }

    /// The #1610 unequal residual: on DNA the single-base placed-gap separation is
    /// disbelieved at the leaf (one spanning segment,
    /// `unequal-length-block-a-placed-gap-is-not-a-separation`); RNA keeps the
    /// [delins; del] split (no `delins.md:47` counterpart off the DNA axis).
    #[test]
    fn unequal_leaf_collapses_single_base_walls_on_dna_only() {
        let got = segs(b"CGCG", b"AAC", 0.01);
        assert_eq!(got, vec![(0, 4, 0, 3)]);
        let rna: Vec<(usize, usize, usize, usize)> = AnchorPeelL1::new(0.01, "peel")
            .segment(b"CGCG", b"AAC", Molecule::Rna, &Provenance::none())
            .iter()
            .map(|s| (s.ref_start, s.ref_end, s.res_start, s.res_end))
            .collect();
        assert_eq!(rna, vec![(0, 2, 0, 2), (3, 4, 3, 3)]);
    }

    /// An individuation claim suppresses the leaf collapse even on DNA — mirrors
    /// `unequal_leaf_collapses_single_base_walls_on_dna_only`.
    #[test]
    fn keep_separate_provenance_keeps_the_unequal_leaf_split_on_dna() {
        use crate::partition::block_ctx::{IndividuationPolicy, Provenance};
        let keep = Provenance {
            individuation: IndividuationPolicy::KeepSeparate,
            examined: Vec::new(),
        };
        let got: Vec<_> = AnchorPeelL1::new(0.01, "peel")
            .segment(b"CGCG", b"AAC", Molecule::Dna, &keep)
            .iter()
            .map(|s| (s.ref_start, s.ref_end, s.res_start, s.res_end))
            .collect();
        assert_eq!(
            got,
            vec![(0, 2, 0, 2), (3, 4, 3, 3)],
            "leaf collapse suppressed"
        );
    }
}
