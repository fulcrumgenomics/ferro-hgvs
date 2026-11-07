//! Control arms: the null/trim segmenter (#10) and always-one-delins typer (#8).

use crate::normalize::seqfirst::align::{AlignmentDag, Step};
use crate::normalize::seqfirst::clcs::{partition_clcs, Rule};
use crate::partition::arm::{DialBConfig, L1Segmenter, L2Typer, MergeContext};
use crate::partition::arms::dial_b::apply_dial_b;
use crate::partition::block_ctx::{FrameContext, IndividuationPolicy, Molecule, Provenance};
use crate::partition::coincidence::{self, CoincidenceParams, CoincidencePiece};
use crate::partition::output::{EditKind, Member, Segment};

/// Length of the shared prefix of two byte slices.
fn shared_prefix(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

/// Length of the shared suffix that does not overlap the already-counted prefix.
fn shared_suffix(a: &[u8], b: &[u8], prefix: usize) -> usize {
    let max = a.len().min(b.len()) - prefix;
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .take(max)
        .take_while(|(x, y)| x == y)
        .count()
}

/// L1 control (#10): no interior cuts. Trims shared ends, returns one segment
/// spanning the changed middle. Pure fn of `(reference, resulting)`.
pub struct TrimOnlyL1;
impl L1Segmenter for TrimOnlyL1 {
    fn name(&self) -> &str {
        "trim-only"
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
        vec![Segment {
            ref_start: p,
            ref_end: reference.len() - s,
            res_start: p,
            res_end: resulting.len() - s,
        }]
    }
}

/// L1 #1 (design §4.2): segment on bases matched in EVERY minimal alignment.
/// Uses the crate's `AlignmentDag`/`Dominators` (adjacent unit-cost metric — the
/// structure the `unchanged-is-read-over-every-minimal-alignment` ledger record
/// discusses). Cut-rule k=1: any single dominator-matched base is a wall, so the
/// segments are the maximal changed runs between walls. Falls back to one
/// whole-span segment when the block is too large for the DAG.
pub struct MinAlignmentIntersectionL1;

impl L1Segmenter for MinAlignmentIntersectionL1 {
    fn name(&self) -> &str {
        "min-alignment-intersection"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        molecule: Molecule,
        provenance: &Provenance,
    ) -> Vec<Segment> {
        // Trim the shared flanks before building the DAG. The flanks are matched
        // in every alignment (so they contribute no changed segment either way),
        // and trimming keeps the banded DAG off the large unchanged padding of a
        // genomic window — the segments produced are identical.
        let p = shared_prefix(reference, resulting);
        let s = shared_suffix(reference, resulting, p);
        let ref_block = &reference[p..reference.len() - s];
        let alt_block = &resulting[p..resulting.len() - s];
        let Some(dag) = AlignmentDag::build(ref_block, alt_block) else {
            return TrimOnlyL1.segment(reference, resulting, molecule, provenance);
        };
        // (ref, alt) offsets matched on every minimal alignment, sorted —
        // relative to the trimmed block, so shift back by the prefix length.
        let matched = dag.dominators().matched;
        let mut segs = Vec::new();
        let (mut prev_ref, mut prev_alt) = (0usize, 0usize);
        let block_ref_len = ref_block.len();
        let block_alt_len = alt_block.len();
        for &(r, a) in &matched {
            let (r, a) = (r as usize, a as usize);
            debug_assert!(
                r >= prev_ref && a >= prev_alt,
                "dominator-matched offsets must be coordinate-monotonic"
            );
            if r > prev_ref || a > prev_alt {
                segs.push(Segment {
                    ref_start: p + prev_ref,
                    ref_end: p + r,
                    res_start: p + prev_alt,
                    res_end: p + a,
                });
            }
            prev_ref = r + 1;
            prev_alt = a + 1;
        }
        if block_ref_len > prev_ref || block_alt_len > prev_alt {
            segs.push(Segment {
                ref_start: p + prev_ref,
                ref_end: p + block_ref_len,
                res_start: p + prev_alt,
                res_end: p + block_alt_len,
            });
        }
        segs
    }
}

/// L1: split at bases matched in EVERY minimal alignment — the node criterion of
/// the ledger ruling `unchanged-is-read-over-every-minimal-alignment`. Where
/// `MaximalSplitL1` commits to one Levenshtein path (so every match on THAT path
/// is a wall), this walls only at columns no minimal alignment can change, which
/// is strictly coarser and alignment-canonical. It is the principled "all minimal
/// alignments" foundation the paper and the ledger both point at — the
/// counterpart, in the merge-only architecture, to a single-path atomic split.
///
/// A column `i` is forced-unchanged iff no cell in row `i` leaves by a `Sub` or
/// `Del` edge on any minimal path; the alt position each forced column matches is
/// read off a canonical walk (prefer diagonal, then deletion, then insertion —
/// the order `out_edges` yields), which is a single minimal alignment and so pins
/// the res coordinate deterministically.
pub struct AllAlignmentSplitL1;
impl L1Segmenter for AllAlignmentSplitL1 {
    fn name(&self) -> &str {
        "all-alignment-split"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        molecule: Molecule,
        provenance: &Provenance,
    ) -> Vec<Segment> {
        // Axis-blind: wall at every forced-unchanged column, so a coincidence
        // separation stands on every molecule. This IS the RNA behaviour; it is the
        // baseline the axis-aware arm collapses down from on a DNA axis.
        forced_unchanged_segments(reference, resulting)
            .unwrap_or_else(|| TrimOnlyL1.segment(reference, resulting, molecule, provenance))
    }
}

/// Segment a `(reference, resulting)` block at its **forced-unchanged columns** —
/// the ledger's `unchanged-is-read-over-every-minimal-alignment` reading. A
/// reference column `i` is a wall iff no minimal alignment leaves row `i` by a
/// `Sub` or `Del` edge (column-based, so it is unchanged even where different
/// minimal alignments match it to *different* alt offsets — the `GACA→AGAT` case
/// the cell-based `Dominators::matched` misses). The segments are the maximal
/// changed runs between walls. `None` when the block is too large for the DAG (the
/// caller falls back to one whole-span segment).
///
/// This is the axis-agnostic wall foundation. The DNA coincidence carve-out is a
/// separate pass ([`disbelieve_dna_coincidences`]) applied on top of these walls;
/// RNA keeps them verbatim.
fn forced_unchanged_segments(reference: &[u8], resulting: &[u8]) -> Option<Vec<Segment>> {
    let p = shared_prefix(reference, resulting);
    let s = shared_suffix(reference, resulting, p);
    let ref_block = &reference[p..reference.len() - s];
    let alt_block = &resulting[p..resulting.len() - s];
    let dag = AlignmentDag::build(ref_block, alt_block)?;
    let n = ref_block.len();
    // forced[i]: no minimal path changes reference column i.
    let mut forced = vec![true; n];
    for (i, j) in dag.cells() {
        if (i as usize) < n
            && dag
                .out_edges(i, j)
                .any(|(_, _, step)| matches!(step, Step::Sub | Step::Del))
        {
            forced[i as usize] = false;
        }
    }
    // Canonical walk to pin the alt offset each forced column matches.
    let mut matched_alt = vec![None; n];
    let (mut i, mut j) = (0u32, 0u32);
    while (i, j) != (dag.ref_len(), dag.alt_len()) {
        let Some((ni, nj, step)) = dag.out_edges(i, j).next() else {
            break;
        };
        if matches!(step, Step::Match) {
            matched_alt[i as usize] = Some(j);
        }
        (i, j) = (ni, nj);
    }
    // Walls are the forced columns; segments are the changed regions between
    // them (block-relative, shifted back by the prefix).
    let mut segs = Vec::new();
    let (mut prev_ref, mut prev_alt) = (0usize, 0usize);
    for col in 0..n {
        if forced[col] {
            if let Some(a) = matched_alt[col] {
                let a = a as usize;
                if col > prev_ref || a > prev_alt {
                    segs.push(Segment {
                        ref_start: p + prev_ref,
                        ref_end: p + col,
                        res_start: p + prev_alt,
                        res_end: p + a,
                    });
                }
                prev_ref = col + 1;
                prev_alt = a + 1;
            }
        }
    }
    let alt_len = alt_block.len();
    if n > prev_ref || alt_len > prev_alt {
        segs.push(Segment {
            ref_start: p + prev_ref,
            ref_end: p + n,
            res_start: p + prev_alt,
            res_end: p + alt_len,
        });
    }
    Some(segs)
}

/// Reverse complement over both alphabets — DNA (`ACGT`) and RNA (`acgu`) — so the
/// whole-span inversion check is correct on either molecule (RNA pairs `a`↔`u`). The
/// bakeoff compares byte structure, so a mixed-case input is left as-is on the bases
/// it does not recognise.
fn revcomp_bytes(b: &[u8]) -> Vec<u8> {
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

/// If the whole common-flank-trimmed block is an exact reverse complement of its
/// reference span (equal length, ≥2 bases), return it as ONE segment. This is
/// `whole-span-reverse-complement-types-as-inv`: the inversion is a property of the
/// entire span, so it must be seen before the forced-unchanged detector can split it
/// at an interior base that coincidentally equals its own complement.
fn whole_span_inversion(reference: &[u8], resulting: &[u8]) -> Option<Segment> {
    let p = shared_prefix(reference, resulting);
    let s = shared_suffix(reference, resulting, p);
    let ref_block = &reference[p..reference.len() - s];
    let alt_block = &resulting[p..resulting.len() - s];
    if ref_block.len() >= 2
        && ref_block.len() == alt_block.len()
        && alt_block == revcomp_bytes(ref_block).as_slice()
    {
        Some(Segment {
            ref_start: p,
            ref_end: reference.len() - s,
            res_start: p,
            res_end: resulting.len() - s,
        })
    } else {
        None
    }
}

/// A composable set of independently-toggleable wall-policy carve-outs — the L1
/// analogue of [`DialBConfig`]. Each carve-out is keyed to a ledger record and can
/// be measured in isolation, so a `base × policy` grid cell is attributable to a
/// specific `(wall-finder, wall-policy)` pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WallPolicy {
    /// W0 (preemption): a whole-span reverse complement is typed as one `inv`,
    /// short-circuiting the base entirely (`whole-span-reverse-complement-types-as-inv`).
    pub whole_span_inv: bool,
    /// W1 (post-transform, DNA only): disbelieve the payload-coincidence separations
    /// `DNA/delins.md:44-47` licenses collapsing, incl. the #1610 unequal-length
    /// residual (`unequal-length-block-a-placed-gap-is-not-a-separation`).
    pub coincidence_collapse: bool,
    /// The bounded-gap reach and embed budget W1's net-deletion pass uses (M4 commit
    /// 2c), in the shared [`CoincidenceParams`] so ONE sweep pair drives this seam and
    /// Dial-B's C2 (`DialBConfig::coincidence`) identically. Only the net-deletion
    /// pass reads it; W1's equal-length (#2174) and placed-gap (#1610) passes stay
    /// single-base by their own antecedents. Not part of [`WallPolicy::name`] — the
    /// label names the enabled carve-outs, not their numeric tuning.
    pub coincidence: CoincidenceParams,
}

impl WallPolicy {
    /// Every carve-out enabled — the configuration the shipped `AxisAware` arm uses.
    pub fn all_on() -> Self {
        Self {
            whole_span_inv: true,
            coincidence_collapse: true,
            coincidence: CoincidenceParams::default(),
        }
    }
    /// No carve-out enabled — the wrapper is a pass-through around its base.
    pub fn all_off() -> Self {
        Self {
            whole_span_inv: false,
            coincidence_collapse: false,
            coincidence: CoincidenceParams::default(),
        }
    }
    /// A short, stable label naming the enabled carve-outs, for the arm name and the
    /// JSONL record: `"W0+W1"`, `"W0"`, `"W1"`, or `"off"`.
    pub fn name(&self) -> &'static str {
        match (self.whole_span_inv, self.coincidence_collapse) {
            (true, true) => "W0+W1",
            (true, false) => "W0",
            (false, true) => "W1",
            (false, false) => "off",
        }
    }
}

/// **Axis-aware wall policy as a reusable wrapper** around any base [`L1Segmenter`]
/// (the pivot). It is itself an `L1Segmenter`, so it drops into the grid unchanged,
/// and it decides WHERE walls go — a wall decision, not a typing (Dial-A) or
/// merge-back (Dial-B) decision. The typer downstream just types whatever segments
/// it is handed, with no need to re-derive structure.
///
/// Precedence, driven by [`WallPolicy`]:
/// * **PRE** — W0 whole-span inversion, highest precedence, short-circuits the base
///   entirely: a span whose whole content is its exact reverse complement is one
///   inverted segment regardless of interior coincidence, uniformly across axes
///   (`whole-span-reverse-complement-types-as-inv`). Without this the base's walls
///   fragment an inversion wherever an interior base matches its own complement.
/// * **base** — the wrapped wall-finder produces the forced-unchanged walls (the RNA
///   reading: every separation stands).
/// * **POST** — W1 coincidence collapse, on a **DNA** molecule only. Disbelieves the
///   single-base payload-coincidence separations `DNA/delins.md:44-47` licenses
///   collapsing ([`disbelieve_dna_coincidences`]); RNA/Protein keep the separations.
///
/// Grounding: RNA — `RNA/delins.md:17` governs unqualified (no `:47` counterpart on
/// the RNA axis), so a coincidence separation is genuine: fragment. DNA — the
/// `delins-payload-coincidence-carve-out-is-coding-dna-scoped` ruling, superseded to
/// all DNA axes by #2155, disbelieves it.
///
/// The shipped arm is `AxisAware::new(Box::new(AllAlignmentSplitL1), WallPolicy::all_on())`
/// — subsuming the old bundled `AxisAwareL1`.
pub struct AxisAware {
    base: Box<dyn L1Segmenter + Send + Sync>,
    policy: WallPolicy,
    /// Precomputed so [`L1Segmenter::name`] can return a `&str`: e.g.
    /// `"axis-aware(all-alignment-split;W0+W1)"`, attributable to a specific
    /// `(wall-finder, wall-policy)` pair.
    name: String,
}

impl AxisAware {
    /// Wrap `base` under `policy`, precomputing the attributable arm name.
    pub fn new(base: Box<dyn L1Segmenter + Send + Sync>, policy: WallPolicy) -> Self {
        let name = format!("axis-aware({};{})", base.name(), policy.name());
        Self { base, policy, name }
    }
}

impl L1Segmenter for AxisAware {
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
        // PRE: whole-span inversion preempts the base — its property is of the whole
        // span and has no term for interior columns, so it must be seen before any
        // wall-finder splits at a coincidental interior match.
        if self.policy.whole_span_inv {
            if let Some(seg) = whole_span_inversion(reference, resulting) {
                return vec![seg];
            }
        }
        let segs = self
            .base
            .segment(reference, resulting, molecule, provenance);
        // POST: DNA-only coincidence collapse. RNA/Protein keep the base's walls
        // verbatim (there is no `:47` counterpart off the DNA axis).
        if self.policy.coincidence_collapse && molecule == Molecule::Dna {
            disbelieve_dna_coincidences(
                segs,
                reference,
                resulting,
                self.policy.coincidence,
                provenance.individuation,
            )
        } else {
            segs
        }
    }
}

/// Whether a segment **reaches** the `DNA/delins.md:44-47` carve-out: it supplies
/// bases *while consuming a different number of reference bases* — i.e. it is an
/// indel-bearing member, the thing that creates an alignment shift and so
/// manufactures the coincidence walls beside it.
///
/// This is `delins-recommendation-reach-when-the-input-arrives-split` verbatim:
/// `:47` reaches a payload-coincidence split "only where some derived member
/// supplies bases while consuming a **different** number of reference bases". Two
/// consequences the too-weak "supplies any base" predicate got wrong:
/// * a pure **substitution** (`ref_len == alt_len`) consumes the SAME count, so it
///   is position-anchored — two subs separated by an unchanged base are genuinely
///   individual (`general.md:33`, `substitution.md`) and must NOT collapse on any
///   axis; and
/// * a pure **deletion** (empty alt) supplies no bases, so a `[3_4del;9del]` split
///   inserts nothing to have coincided and stays split.
///
/// The gap-bearing test itself now lives in `coincidence::CoincidencePiece::
/// is_gap_bearing` (shared with Dial-B), and the licensing decision in
/// `coincidence::{base_carve_out, placed_gap_extension}`.
///
/// Collapse the coincidence separations `DNA/delins.md:44-47` licenses merging, on a
/// DNA axis, in **two passes** — the L1 mirror of Dial-B's C2 (bounded net-deletion) +
/// C3/#2174 (single-base) composition:
///
/// * **Pass 1 — bounded-gap net-deletion (M4 commit 2c).** Groups segments into
///   maximal runs separated by `1..=params.max_separation` forced-unchanged bases (the
///   shared [`coincidence::bounded_gap_runs_over`]) and merges each run the embed-gated
///   net-deletion predicate [`coincidence::base_carve_out`] accepts — byte-aware, so a
///   wider wall must actually spell the payload within `params.mismatch_budget`. Same
///   grouper, same predicate, and the SAME [`CoincidenceParams`] as Dial-B's C2, so a
///   run W1 merges here is byte-for-byte the run C2 would merge at Dial-B (the W1/C2
///   parity the sweep and the equivalence acceptance depend on).
/// * **Pass 2 — single-base equal-length (#2174) + placed-gap (#1610).** The residual
///   is grouped at single-base separation only and collapsed by
///   [`run_collapse_is_licensed`] (the byte-blind `base_carve_out_or_equal_length ||
///   placed_gap_extension`). Both are single-base by their own antecedents and are left
///   exactly as before 2c; a bounded run pass 1 declined falls through here and its
///   single-base sub-runs are handled unchanged.
///
/// At `max_separation == 1` pass 1's net-deletion merges are a subset of pass 2's
/// byte-blind ones, so the two-pass result reduces to the pre-2c single-base behaviour
/// — 2c's movement is exactly the new `2..=max` net-deletion merges. A separation of
/// two or more unchanged bases that is NOT a net-deletion coincidence (an equal-length
/// #2174 block, two position-anchored subs) is still kept split: pass 1's net-deletion
/// predicate declines it and pass 2 never groups across a >1 gap.
fn disbelieve_dna_coincidences(
    segs: Vec<Segment>,
    reference: &[u8],
    resulting: &[u8],
    params: CoincidenceParams,
    policy: IndividuationPolicy,
) -> Vec<Segment> {
    let segs = collapse_bounded_net_del(segs, reference, resulting, params, policy);
    collapse_single_base_coincidences(segs, resulting, policy)
}

/// Pass 1: merge each bounded-gap (`1..=max_separation`) run the embed-gated
/// net-deletion carve-out [`coincidence::base_carve_out`] accepts. The L1-seam mirror
/// of Dial-B's `pass_c2_coincidence`. `pub(crate)` so the W1/C2 parity pin
/// (`dial_b` tests) can drive this seam and C2's directly from one run.
pub(crate) fn collapse_bounded_net_del(
    segs: Vec<Segment>,
    reference: &[u8],
    resulting: &[u8],
    params: CoincidenceParams,
    policy: IndividuationPolicy,
) -> Vec<Segment> {
    if segs.len() < 2 {
        return segs;
    }
    let spans: Vec<(usize, usize)> = segs.iter().map(|s| (s.ref_start, s.ref_end)).collect();
    let mut out: Vec<Segment> = Vec::new();
    let mut next = 0usize;
    for (start, end) in coincidence::bounded_gap_runs_over(&spans, params.max_separation) {
        out.extend_from_slice(&segs[next..start]);
        next = end;
        let run = &segs[start..end];
        let pieces = project_segment_pieces(run, resulting);
        if coincidence::base_carve_out(&pieces, reference, params.mismatch_budget, policy) {
            out.push(merge_run(run));
        } else {
            out.extend_from_slice(run);
        }
    }
    out.extend_from_slice(&segs[next..]);
    out
}

/// Pass 2: the pre-2c single-base collapse — group at single-base separation and merge
/// each run [`run_collapse_is_licensed`] accepts (equal-length #2174 + placed-gap
/// #1610). Unchanged from before 2c; runs pass 1 already merged appear here as lone
/// segments and are not re-grouped (their neighbours are `> 1` or `0` away).
fn collapse_single_base_coincidences(
    segs: Vec<Segment>,
    resulting: &[u8],
    policy: IndividuationPolicy,
) -> Vec<Segment> {
    if segs.len() < 2 {
        return segs;
    }
    let mut out: Vec<Segment> = Vec::new();
    let mut run_start = 0usize;
    for i in 0..segs.len() {
        // A run ends here unless the NEXT segment is separated from this one by
        // exactly one forced-unchanged reference base (a single coincidence wall).
        // A separation of two or more unchanged bases is genuine and breaks the run.
        let run_ends_here = i + 1 == segs.len() || segs[i + 1].ref_start != segs[i].ref_end + 1;
        if !run_ends_here {
            continue;
        }
        let run = &segs[run_start..=i];
        if run.len() >= 2 && run_collapse_is_licensed(run, resulting, policy) {
            out.push(merge_run(run));
        } else {
            out.extend_from_slice(run);
        }
        run_start = i + 1;
    }
    out
}

/// Project a run of segments into the shared molecule-blind [`CoincidencePiece`] view
/// the coincidence predicates decide from — `alt` is the resulting-sequence slice each
/// segment supplies. Shared by both passes so their projection cannot drift.
fn project_segment_pieces(run: &[Segment], resulting: &[u8]) -> Vec<CoincidencePiece> {
    run.iter()
        .map(|s| CoincidencePiece {
            ref_start: s.ref_start,
            ref_end: s.ref_end,
            alt: resulting[s.res_start..s.res_end].to_vec(),
        })
        .collect()
}

/// Whether a single-base-separated run may be collapsed on a DNA axis. Run-based
/// (evaluated over the maximal single-base-separated run, not per candidate wall) —
/// equivalent for these carve-outs, since every separation in the run is already one
/// unchanged base, and simpler than the spec's recommended per-wall disbelief
/// predicate. Two licensing paths, both keyed to ledger records:
///
/// * **Base rule** (`delins-recommendation-reach-when-the-input-arrives-split` +
///   `delins-merge-vs-individual-gap-two-or-more`): the run contains an indel-bearing
///   member (so its interior walls are alignment coincidences, not position-anchored)
///   AND the merged span is not a net insertion.
/// * **#1610 extension** (`unequal-length-block-a-placed-gap-is-not-a-separation`,
///   widened to all DNA by #2155): an *extension* of the base rule, not a second
///   carve-out. A run with no indel-bearing member still collapses when all four of
///   the record's conditions hold — (a) the merged span is a net deletion, (b) every
///   separation is exactly one unchanged base (guaranteed by the run grouping), (c)
///   at least one member supplies bases, and (d) every member renders as a `delins`
///   or a pure deletion, with at most one pure-deletion residual (the "one member
///   wider" gap the aligner placed). This admits `CGCG->AAC` = `[CG->AA; G->del]`,
///   which the base rule wrongly declines because neither member is indel-bearing.
///
/// See [`disbelieve_dna_coincidences`] for how the run grouping bounds these.
fn run_collapse_is_licensed(
    run: &[Segment],
    resulting: &[u8],
    policy: IndividuationPolicy,
) -> bool {
    // Decide from the shared core (decision #1: one predicate, called here and by
    // Dial-B's C3, so the two cannot drift). This pass-2 predicate applies both the
    // equal-length-inclusive base carve-out (#2174) and the #1610 placed-gap
    // extension; both are byte-blind (they read only length/emptiness). Pass 1's
    // net-deletion carve-out is the byte-aware one.
    let pieces = project_segment_pieces(run, resulting);
    coincidence::base_carve_out_or_equal_length(&pieces, policy)
        || coincidence::placed_gap_extension(&pieces, policy)
}

/// The single spanning segment covering an entire run (its interior coincidence
/// walls included) — the `delins` form the DNA carve-out prefers.
fn merge_run(run: &[Segment]) -> Segment {
    let first = &run[0];
    let last = &run[run.len() - 1];
    Segment {
        ref_start: first.ref_start,
        ref_end: last.ref_end,
        res_start: first.res_start,
        res_end: last.res_end,
    }
}

/// Segment a variant with the paper-faithful cLCS rule (design §4.4 / M4). The
/// blocks come back relative to the common-flank-trimmed block, so they are
/// shifted back into absolute coordinates. Falls back to one whole-span segment
/// when the block exceeds the rule's size cap.
fn clcs_segments(
    reference: &[u8],
    resulting: &[u8],
    rule: Rule,
    molecule: Molecule,
    provenance: &Provenance,
) -> Vec<Segment> {
    let p = shared_prefix(reference, resulting);
    let s = shared_suffix(reference, resulting, p);
    let ref_block = &reference[p..reference.len() - s];
    let alt_block = &resulting[p..resulting.len() - s];
    match partition_clcs(ref_block, alt_block, rule) {
        Some(blocks) => blocks
            .iter()
            .map(|b| Segment {
                ref_start: p + b.ref_start as usize,
                ref_end: p + b.ref_end as usize,
                res_start: p + b.alt_start as usize,
                res_end: p + b.alt_end as usize,
            })
            .collect(),
        None => TrimOnlyL1.segment(reference, resulting, molecule, provenance),
    }
}

/// L1 arm: the paper-faithful cLCS-graph rule (LCS/indel metric), no inversion
/// edge. The LCS/indel-metric counterpart to `MinAlignmentIntersectionL1`'s
/// unit-cost Levenshtein.
pub struct ClcsL1;
impl L1Segmenter for ClcsL1 {
    fn name(&self) -> &str {
        "clcs"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        molecule: Molecule,
        provenance: &Provenance,
    ) -> Vec<Segment> {
        clcs_segments(reference, resulting, Rule::Clcs, molecule, provenance)
    }
}

/// L1 arm: the cLCS rule plus one inversion edge — a whole-span reverse
/// complement is segmented as a single inverted block (which the ladder then
/// types as `Inv`), where the plain rule would fragment it into substitutions.
pub struct ClcsInvL1;
impl L1Segmenter for ClcsInvL1 {
    fn name(&self) -> &str {
        "clcs-inv"
    }
    fn segment(
        &self,
        reference: &[u8],
        resulting: &[u8],
        molecule: Molecule,
        provenance: &Provenance,
    ) -> Vec<Segment> {
        clcs_segments(reference, resulting, Rule::ClcsInv, molecule, provenance)
    }
}

/// L2 control (#8): always one delins over the segment (Identity if empty).
pub struct AlwaysOneDelinsL2;
impl L2Typer for AlwaysOneDelinsL2 {
    fn name(&self) -> &str {
        "always-one-delins"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        _r: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let ref_empty = seg.ref_start == seg.ref_end;
        let res_empty = seg.res_start == seg.res_end;
        if ref_empty && res_empty {
            return vec![Member {
                kind: EditKind::Identity,
                ref_start: seg.ref_start,
                ref_end: seg.ref_end,
                inserted: Vec::new(),
            }];
        }
        vec![Member {
            kind: EditKind::Delins,
            ref_start: seg.ref_start,
            ref_end: seg.ref_end,
            inserted: resulting[seg.res_start..seg.res_end].to_vec(),
        }]
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        _r: &[u8],
        _s: &[u8],
        _ctx: &MergeContext,
        _c: &DialBConfig,
    ) -> Vec<Member> {
        m
    }
}

fn revcomp(b: &[u8]) -> Vec<u8> {
    b.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            other => *other,
        })
        .collect()
}

/// L2 Dial A #1 (design §5.2/§5.3): the recommended-form ladder. Types one
/// segment by spec-stated precedence, deriving every choice from the sequences:
///
/// 1. empty↔empty → `Identity`
/// 2. pure insertion (no reference consumed) → `Dup` if the inserted bases copy
///    the immediately-5' reference bases (`duplication-must-ranks-the-label-not-
///    the-partition`, ≥2 bases), else `Ins`
/// 3. pure deletion (nothing inserted) → `Del`
/// 4. equal length, whole span is the reverse complement (≥2) → `Inv`
///    (`whole-span-reverse-complement-types-as-inv`)
/// 5. single reference base replaced by a single base → `Sub`
/// 6. otherwise → one spanning `Delins`
///
/// Merge-back (Dial B) is identity here; the carve-outs are a separate concern.
pub struct LadderL2;

impl LadderL2 {
    fn type_one(seg: &Segment, reference: &[u8], resulting: &[u8]) -> Member {
        let ref_len = seg.ref_end - seg.ref_start;
        let alt = resulting[seg.res_start..seg.res_end].to_vec();
        let member = |kind: EditKind, inserted: Vec<u8>| Member {
            kind,
            ref_start: seg.ref_start,
            ref_end: seg.ref_end,
            inserted,
        };
        if ref_len == 0 && alt.is_empty() {
            return member(EditKind::Identity, Vec::new());
        }
        if ref_len == 0 {
            // A pure insertion. Dup iff it copies the immediately-5' bases. n >= 1: a
            // single-nucleotide 5'-copy insertion is a dup too (`duplication.md:5`
            // "one or more", the `c.20dup` example, `c.19_20insT` marked not allowed —
            // operator ruling 2026-08-22, Side A). Kept in step with `type_run` in
            // `l2_typers.rs` and `conformance.rs::copies_immediate_5prime`.
            let n = alt.len();
            if n >= 1
                && seg.ref_start >= n
                && reference[seg.ref_start - n..seg.ref_start] == alt[..]
            {
                return member(EditKind::Dup, alt);
            }
            return member(EditKind::Ins, alt);
        }
        if alt.is_empty() {
            return member(EditKind::Del, Vec::new());
        }
        let ref_slice = &reference[seg.ref_start..seg.ref_end];
        if alt.len() == ref_len && ref_len >= 2 && alt == revcomp(ref_slice) {
            return member(EditKind::Inv, Vec::new());
        }
        if ref_len == 1 && alt.len() == 1 {
            return member(EditKind::Sub, alt);
        }
        member(EditKind::Delins, alt)
    }
}

/// L2 Dial A #6 (design §5.2): position-wise typing — the "live" back-half
/// behavior. An equal-length segment becomes one `Sub` per differing column
/// (unchanged columns emit nothing); an unequal-length segment falls back to one
/// spanning `Delins`. Blind to dup/inv by construction — useful as the contrast
/// that isolates what the ladder's typing buys over naive position pairing.
pub struct PositionWiseL2;

impl L2Typer for PositionWiseL2 {
    fn name(&self) -> &str {
        "position-wise"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        let ref_len = seg.ref_end - seg.ref_start;
        let alt_len = seg.res_end - seg.res_start;
        if ref_len != alt_len {
            // Undefined position-wise; fall back to one spanning delins.
            return vec![Member {
                kind: EditKind::Delins,
                ref_start: seg.ref_start,
                ref_end: seg.ref_end,
                inserted: resulting[seg.res_start..seg.res_end].to_vec(),
            }];
        }
        // One Sub per column that actually differs; unchanged columns emit no
        // member (the applier copies those reference bases between members).
        let mut members = Vec::new();
        for i in 0..ref_len {
            let ref_base = reference[seg.ref_start + i];
            let alt_base = resulting[seg.res_start + i];
            if ref_base != alt_base {
                members.push(Member {
                    kind: EditKind::Sub,
                    ref_start: seg.ref_start + i,
                    ref_end: seg.ref_start + i + 1,
                    inserted: vec![alt_base],
                });
            }
        }
        members
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        apply_dial_b(m, reference, resulting, ctx, cfg)
    }
}

impl L2Typer for LadderL2 {
    fn name(&self) -> &str {
        "recommended-form-ladder"
    }
    fn type_segment(
        &self,
        seg: &Segment,
        reference: &[u8],
        resulting: &[u8],
        _f: &FrameContext,
    ) -> Vec<Member> {
        vec![Self::type_one(seg, reference, resulting)]
    }
    fn merge_back(
        &self,
        m: Vec<Member>,
        reference: &[u8],
        resulting: &[u8],
        ctx: &MergeContext,
        cfg: &DialBConfig,
    ) -> Vec<Member> {
        apply_dial_b(m, reference, resulting, ctx, cfg)
    }
}

#[cfg(test)]
mod axis_aware_tests {
    use super::*;

    /// A segment as `(ref_start, ref_end, res_start, res_end)` for terse asserts.
    fn tup(s: &Segment) -> (usize, usize, usize, usize) {
        (s.ref_start, s.ref_end, s.res_start, s.res_end)
    }
    fn tups(segs: &[Segment]) -> Vec<(usize, usize, usize, usize)> {
        segs.iter().map(tup).collect()
    }

    /// The shipped axis-aware arm: `all-alignment-split` under an all-on wall policy.
    /// This is the configuration the deleted bundled `AxisAwareL1` was, so the guards
    /// below exercise it in that arm's place.
    fn axis_aware() -> AxisAware {
        AxisAware::new(Box::new(AllAlignmentSplitL1), WallPolicy::all_on())
    }

    /// The #2155 block `CTTAGTTA -> AAACAAAC` (equal length, edit distance 7, two
    /// interior A's forced-unchanged in every minimal alignment). The forced-unchanged
    /// walls fragment it into three pieces — the RNA answer and the input to the DNA
    /// carve-out.
    #[test]
    fn forced_unchanged_walls_fragment_the_2155_block() {
        let segs = forced_unchanged_segments(b"CTTAGTTA", b"AAACAAAC").expect("dag builds");
        // seg1 CTT->AA, wall A@3, seg2 GTT->CAA, wall A@7, seg3 ins C.
        assert_eq!(tups(&segs), vec![(0, 3, 0, 2), (4, 7, 3, 6), (8, 8, 7, 8)]);
    }

    /// RNA keeps every forced-unchanged wall (no `:47` counterpart) — fragments.
    #[test]
    fn axis_aware_rna_fragments_the_2155_block() {
        let segs =
            axis_aware().segment(b"CTTAGTTA", b"AAACAAAC", Molecule::Rna, &Provenance::none());
        assert_eq!(segs.len(), 3, "RNA must keep the three members");
        assert_eq!(tups(&segs), vec![(0, 3, 0, 2), (4, 7, 3, 6), (8, 8, 7, 8)]);
    }

    /// DNA disbelieves the single-base coincidence separations (`delins.md:44-47`)
    /// and collapses to one spanning delins over the whole block.
    #[test]
    fn axis_aware_dna_collapses_the_2155_block() {
        let segs =
            axis_aware().segment(b"CTTAGTTA", b"AAACAAAC", Molecule::Dna, &Provenance::none());
        assert_eq!(segs.len(), 1, "DNA must collapse to one spanning member");
        assert_eq!(tup(&segs[0]), (0, 8, 0, 8));
    }

    /// An individuation claim suppresses the DNA coincidence carve-out even though
    /// the molecule is DNA — the changes were reported/observed individually, so
    /// the separation is not an alignment artifact.
    #[test]
    fn keep_separate_provenance_keeps_the_2155_block_fragmented_on_dna() {
        use crate::partition::block_ctx::{IndividuationPolicy, Provenance};
        let keep = Provenance {
            individuation: IndividuationPolicy::KeepSeparate,
            examined: Vec::new(),
        };
        let segs = axis_aware().segment(b"CTTAGTTA", b"AAACAAAC", Molecule::Dna, &keep);
        assert_eq!(
            segs.len(),
            3,
            "an individuation claim suppresses W1 even on DNA"
        );
        let none = Provenance::none();
        let dna = axis_aware().segment(b"CTTAGTTA", b"AAACAAAC", Molecule::Dna, &none);
        assert_eq!(dna.len(), 1, "∅ provenance still collapses");
    }

    /// Two substitutions separated by ONE unchanged base stay individual even on a
    /// DNA axis (`CAC -> GAT`): the interior A is position-anchored, not an alignment
    /// coincidence — no member is indel-bearing, so `general.md:33` / `substitution.md`
    /// govern. This is the case the too-weak "supplies any base" predicate collapsed
    /// wrongly; `CoincidencePiece::is_gap_bearing` requires a DIFFERENT ref/alt count.
    #[test]
    fn axis_aware_dna_keeps_separated_substitutions() {
        let dna = axis_aware().segment(b"CAC", b"GAT", Molecule::Dna, &Provenance::none());
        let rna = axis_aware().segment(b"CAC", b"GAT", Molecule::Rna, &Provenance::none());
        assert_eq!(
            dna, rna,
            "separated subs are molecule-independent (no coincidence)"
        );
        assert_eq!(
            tups(&dna),
            vec![(0, 1, 0, 1), (2, 3, 2, 3)],
            "two individual subs"
        );
    }

    /// A pure-deletion split (`ACACA -> AAA`, delete both C's) stays split even on a
    /// DNA axis: nothing was re-aligned, so `general.md:33` governs unqualified
    /// (`delins-recommendation-reach-when-the-input-arrives-split`).
    #[test]
    fn axis_aware_dna_keeps_a_pure_deletion_split() {
        let dna = axis_aware().segment(b"ACACA", b"AAA", Molecule::Dna, &Provenance::none());
        let rna = axis_aware().segment(b"ACACA", b"AAA", Molecule::Rna, &Provenance::none());
        assert_eq!(dna, rna, "a pure-deletion split is molecule-independent");
        assert_eq!(tups(&dna), vec![(1, 2, 1, 1), (3, 4, 2, 2)]);
    }

    /// The DNA carve-out is direction-scoped: a net insertion is not collapsed even
    /// when its separations are single-base coincidences
    /// (`delins-merge-vs-individual-gap-two-or-more`; `duplication.md:90-92` keeps a
    /// net insertion split). Tested on the predicate directly with a net-insertion run
    /// (merged ref span 3, alt span 5).
    #[test]
    fn dna_carve_out_leaves_a_net_insertion_split() {
        let net_ins = vec![
            Segment {
                ref_start: 0,
                ref_end: 1,
                res_start: 0,
                res_end: 2,
            },
            Segment {
                ref_start: 2,
                ref_end: 3,
                res_start: 3,
                res_end: 5,
            },
        ];
        // Single-base gap (2 == 1 + 1) but net insertion (5 > 3) -> not collapsed.
        // `resulting` must be long enough for the segments' res indices (the projection
        // now reads the real alt bytes); its content is irrelevant to the net-insertion
        // verdict, which is byte-blind.
        let out = disbelieve_dna_coincidences(
            net_ins.clone(),
            b"",
            b"AAGCC",
            CoincidenceParams::default(),
            IndividuationPolicy::Unspecified,
        );
        assert_eq!(out, net_ins, "net-insertion coincidence run stays split");
    }

    /// A whole-span reverse complement is ONE segment, uniformly across axes, even
    /// when an interior base coincidentally equals its own complement (which the
    /// forced-unchanged detector would otherwise wall). `AAGCTA -> TAGCTT` is an exact
    /// revcomp; without the whole-span check it fragments into subs
    /// (`whole-span-reverse-complement-types-as-inv`).
    #[test]
    fn axis_aware_types_a_whole_span_revcomp_as_one_segment() {
        for mol in [Molecule::Dna, Molecule::Rna] {
            let segs = axis_aware().segment(b"AAGCTA", b"TAGCTT", mol, &Provenance::none());
            assert_eq!(segs.len(), 1, "{mol:?}: whole-span revcomp is one segment");
            assert_eq!(tup(&segs[0]), (0, 6, 0, 6));
        }
    }

    /// The whole-span check fires only on an EXACT revcomp; a near-revcomp that is not
    /// one is segmented normally (no false inversion).
    #[test]
    fn whole_span_check_ignores_a_non_revcomp() {
        // TAGCTT with one base changed is no longer revcomp(AAGCTA); segment normally.
        assert!(whole_span_inversion(b"AAGCTA", b"TAGCTA").is_none());
        assert!(whole_span_inversion(b"AAGCTA", b"TAGCTT").is_some());
    }

    /// A two-base-separated run that is a NET INSERTION stays split even under the M4
    /// widening — `base_carve_out` declines a net insertion before the separation
    /// bound ever matters (`DNA/duplication.md:90-92`). Renamed from the pre-2c
    /// `dna_carve_out_keeps_a_two_base_separation`: that test passed for the wrong
    /// reason (it read the split as coming from the 2-base wall, which 2c overturns for
    /// net deletions — the run stays split because it is a net INSERTION).
    #[test]
    fn w1_keeps_a_two_base_net_insertion_split() {
        // Two delins runs separated by a 2-base gap (ref_end 1, next ref_start 3),
        // net insertion (alt 5 > ref 4). `reference` is unused: net insertion is
        // declined before the span slice.
        let two_gap = vec![
            Segment {
                ref_start: 0,
                ref_end: 1,
                res_start: 0,
                res_end: 2,
            },
            Segment {
                ref_start: 3,
                ref_end: 4,
                res_start: 4,
                res_end: 5,
            },
        ];
        let out = disbelieve_dna_coincidences(
            two_gap.clone(),
            b"",
            b"AAGCC",
            CoincidenceParams::default(),
            IndividuationPolicy::Unspecified,
        );
        assert_eq!(
            out, two_gap,
            "a net-insertion run stays split at any separation"
        );
    }

    /// M4 (2c): W1's bounded net-deletion pass merges a coincidence separated by TWO
    /// unchanged bases — the widening past single-base. `CGGAG -> CGA` is
    /// `[delins CG->C; 2-base wall "GA"; del G]`: net-deletion, payload "CGA" embeds in
    /// the span "CGGAG" with 0 substitutions, so `base_carve_out` accepts it. At
    /// `max_separation == 1` the SAME run stays split (the pre-2c behaviour), so the
    /// merge is attributable to the widening and nothing else.
    #[test]
    fn w1_bounded_pass_merges_a_two_base_net_deletion_coincidence() {
        let two_base_net_del = vec![
            Segment {
                ref_start: 0,
                ref_end: 2,
                res_start: 0,
                res_end: 1,
            }, // CG -> C (gap-bearing)
            Segment {
                ref_start: 4,
                ref_end: 5,
                res_start: 3,
                res_end: 3,
            }, // G -> del
        ];
        let widened = disbelieve_dna_coincidences(
            two_base_net_del.clone(),
            b"CGGAG",
            b"CGA",
            CoincidenceParams::default(), // max_separation 8
            IndividuationPolicy::Unspecified,
        );
        assert_eq!(
            widened.len(),
            1,
            "2c: a 2-base net-deletion coincidence merges to one spanning delins"
        );
        assert_eq!(tup(&widened[0]), (0, 5, 0, 3));

        let single_base = disbelieve_dna_coincidences(
            two_base_net_del.clone(),
            b"CGGAG",
            b"CGA",
            CoincidenceParams {
                max_separation: 1,
                mismatch_budget: 1,
            },
            IndividuationPolicy::Unspecified,
        );
        assert_eq!(
            single_base, two_base_net_del,
            "max_separation 1 reproduces the pre-2c split (attribution)"
        );
    }

    /// Pass-2 fall-through (M4 2c, Q2 interaction): a bounded run pass 1 DECLINES on
    /// the embed test still lets pass 2 collapse a single-base sub-run inside it.
    /// `AAAAAAAAA -> XYAAAA` = `[delins AAA->XY; gap1; del A; gap2; del A]`, trailing A
    /// unchanged. Pass 1 groups the whole `[A,B,C]` bounded run (gaps 1 and 2) and
    /// declines — payload "XYAAA" needs 2 substitutions to embed in "AAAAAAAA" (budget
    /// 1). Pass 2 then groups only the single-base sub-run `[A,B]` (the B–C gap is 2)
    /// and collapses it byte-blind (gap-bearing, net-deletion), leaving C. Result: 2
    /// segments — not 1 (pass 1 did not over-merge) and not 3 (pass 2 did collapse the
    /// residual sub-run).
    #[test]
    fn w1_pass1_declines_on_embed_then_pass2_collapses_a_single_base_subrun() {
        let run = vec![
            Segment {
                ref_start: 0,
                ref_end: 3,
                res_start: 0,
                res_end: 2,
            }, // AAA -> XY
            Segment {
                ref_start: 4,
                ref_end: 5,
                res_start: 3,
                res_end: 3,
            }, // A -> del
            Segment {
                ref_start: 7,
                ref_end: 8,
                res_start: 5,
                res_end: 5,
            }, // A -> del
        ];
        let out = disbelieve_dna_coincidences(
            run,
            b"AAAAAAAAA",
            b"XYAAAA",
            CoincidenceParams::default(),
            IndividuationPolicy::Unspecified,
        );
        assert_eq!(
            tups(&out),
            vec![(0, 5, 0, 3), (7, 8, 5, 5)],
            "pass 1 declines on embed, pass 2 collapses the single-base sub-run [A,B]"
        );
    }

    /// The #1610 shape (`unequal-length-block-a-placed-gap-is-not-a-separation`,
    /// widened to all DNA by #2155): `CGCG -> AAC` segments to `[CG->AA (equal-length
    /// delins); G->del (pure deletion)]`, one unchanged base between. Neither member
    /// is indel-bearing, so the base rule declines it — but it is a net deletion whose
    /// members are a delins plus one pure-deletion residual, so the #1610 extension
    /// collapses it to ONE spanning delins on a DNA axis. RNA keeps the two members
    /// (no `:47` counterpart off the DNA axis).
    #[test]
    fn axis_aware_dna_collapses_the_1610_unequal_length_block() {
        let dna = axis_aware().segment(b"CGCG", b"AAC", Molecule::Dna, &Provenance::none());
        assert_eq!(
            dna.len(),
            1,
            "DNA #1610 must collapse to one spanning delins"
        );
        assert_eq!(tup(&dna[0]), (0, 4, 0, 3));

        let rna = axis_aware().segment(b"CGCG", b"AAC", Molecule::Rna, &Provenance::none());
        assert_eq!(rna.len(), 2, "RNA keeps the #1610 block split");
    }

    /// Policy toggle: with W1 (coincidence collapse) OFF but W0 (whole-span inversion)
    /// ON, the octamer must NOT collapse on a DNA axis — the base walls stand.
    #[test]
    fn wall_policy_w1_off_does_not_collapse_the_octamer() {
        let arm = AxisAware::new(
            Box::new(AllAlignmentSplitL1),
            WallPolicy {
                coincidence_collapse: false,
                whole_span_inv: true,
                coincidence: CoincidenceParams::default(),
            },
        );
        let segs = arm.segment(b"CTTAGTTA", b"AAACAAAC", Molecule::Dna, &Provenance::none());
        assert_eq!(
            segs.len(),
            3,
            "W1 off: DNA keeps the base's three walls (no collapse)"
        );
    }

    /// Policy toggle: with W0 (whole-span inversion) OFF but W1 ON, a whole-span
    /// revcomp is NOT preempted — it fragments through the base wall-finder instead of
    /// returning one inverted segment.
    #[test]
    fn wall_policy_w0_off_does_not_preempt_a_whole_span_revcomp() {
        let arm = AxisAware::new(
            Box::new(AllAlignmentSplitL1),
            WallPolicy {
                whole_span_inv: false,
                coincidence_collapse: true,
                coincidence: CoincidenceParams::default(),
            },
        );
        let segs = arm.segment(b"AAGCTA", b"TAGCTT", Molecule::Dna, &Provenance::none());
        assert!(
            segs.len() > 1,
            "W0 off: whole-span revcomp is not preempted, so it fragments (got {} segments)",
            segs.len()
        );
    }
}
