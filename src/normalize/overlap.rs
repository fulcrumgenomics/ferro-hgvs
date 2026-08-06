//! Detect overlapping cis-allele sub-variants.
//!
//! See `docs/superpowers/specs/2026-05-05-allele-overlap-same-position-design.md`.
//!
//! Two cis-phase, same-accession sub-variants overlap when they describe edits
//! on the same reference territory — a case with no canonical HGVS form.
//! Rather than silently picking a winner, ferro preserves the input verbatim
//! and emits one [`NormalizationWarning::OverlapConflict`] per overlap (strict
//! mode promotes it to an error).
//!
//! Two detectors run at different points in the allele pipeline:
//!
//! - [`detect_overlap_conflicts`] — *post-shift* coincident bounds: span edits
//!   (`sub`/`del`/`delins`/`dup`/`inv`/`repeat`) whose `(region, start, end)`
//!   keys are identical. Insertions are excluded here (they anchor at a
//!   boundary, not a single-base span).
//! - [`detect_insertion_overlaps`] — *pre-merge* junction overlaps: two
//!   junction-anchored edits at one junction, or one interior to a span edit
//!   (mutalyzer `EOVERLAP`, #486). Must run before the normalizer's merge step
//!   collapses overlapping cis edits into one combined edit. "Junction-
//!   anchored" covers a true `ins` plus the `dup` and `repeat` spellings (which
//!   land their extra copies at the junction after their own span), so the
//!   detector sees the same conflict whichever way the pipeline has spelled it.

use std::collections::BTreeMap;

use crate::hgvs::edit::NaEdit;
use crate::hgvs::interval::{Interval, UncertainBoundary};
use crate::hgvs::location::{CdsPos, GenomePos, RnaPos, TxPos};
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, HgvsVariant, LocEdit};
use crate::normalize::merge::Region;
use crate::normalize::NormalizationWarning;

/// Detect coincident-bounds groups in a cis allele.
///
/// Returns one warning per group of `>= 2` same-accession sub-variants whose
/// `(coordinate-system region, signed start, signed end)` keys are identical.
/// Insertions are excluded — they anchor at boundaries (`[end, start]`) and
/// have no single-base location to coincide on. Trans / mosaic / chimeric /
/// unknown phases short-circuit (the warning only applies to a same-haplotype
/// allele, where the conflict is real).
///
/// Input `variants` are the post-normalization sub-variants of an
/// [`AlleleVariant`]. The pass is purely observational: it does not mutate
/// or reorder its input. Warnings are emitted in deterministic key order
/// (BTreeMap iteration), so two equivalent inputs yield identical warning
/// sequences regardless of source-line ordering of the conflicting edits.
pub(crate) fn detect_overlap_conflicts(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<NormalizationWarning> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return Vec::new();
    }

    // Group sub-variants by (accession, coord_system, region, start, end). The
    // `Vec<usize>` collects 0-based indices into `variants` in input order so
    // each group's `edit_kinds` reflect the source ordering of the conflict.
    let mut groups: BTreeMap<GroupKey, Vec<usize>> = BTreeMap::new();
    for (idx, variant) in variants.iter().enumerate() {
        let Some(key) = group_key(variant) else {
            continue;
        };
        groups.entry(key).or_default().push(idx);
    }

    let mut warnings = Vec::new();
    for (key, indices) in &groups {
        if indices.len() < 2 {
            continue;
        }
        let edit_kinds: Vec<String> = indices
            .iter()
            .filter_map(|&i| edit_kind(&variants[i]).map(|s| s.to_string()))
            .collect();
        // group_key already filtered to variants with a known edit kind, so
        // edit_kinds.len() == indices.len() in practice. Preserved as a
        // debug assert against future drift.
        debug_assert_eq!(edit_kinds.len(), indices.len());
        let location_str = location_for_variant(&variants[indices[0]])
            .expect("group_key established a renderable location");

        warnings.push(NormalizationWarning::OverlapConflict {
            accession: key.accession.clone(),
            coordinate_system: key.coord_system.to_string(),
            location: location_str,
            edit_kinds,
        });
    }
    warnings
}

/// Detect overlaps that involve at least one junction-anchored edit within a
/// cis allele.
///
/// An insertion `a_(a+1)ins…` occupies the zero-width junction between
/// reference positions `a` and `a+1`. It overlaps:
///
/// - **another insertion** at the *same* junction (`[4_5insT;4_5insA]`); and
/// - **a span edit** (`del`/`delins`/`dup`/`inv`/`sub`/`repeat`) whose
///   reference range *strictly encloses* the junction, i.e.
///   `range.start <= a` and `a + 1 <= range.end`
///   (`[274_275delinsT;274_275insA]`).
///
/// An insertion abutting a span's edge — e.g. `100_101ins` next to a single-
/// base `100` substitution — is *not* interior, so it does not overlap. This
/// keeps the spec-valid `[273_274insT;274G>T;274_275insA]` accepted.
///
/// A **`dup` and a `repeat` also occupy a junction** (see
/// [`junction_writing_kind`]): each writes its extra copies directly 3' of its
/// own span, i.e. at the junction after `end`. Both are therefore registered as
/// junction occupants *in addition to* being spans, so `[5_9inv;5_6dup]` and
/// `[1005_1009inv;1005_1006A[4]]` are each flagged just as the `ins`-spelled
/// input they normalize from is.
///
/// Without this the detector is spelling-sensitive and normalization is not
/// idempotent on a conflicting allele. The first pass respells the interior
/// `insAA` — as a `dup` on the CDS axis, as a `repeat` on the genomic axis —
/// the second pass no longer recognises either as a conflict, and so it
/// reorders members that the first pass deliberately left in authored order
/// (`#395`). Covering only one of the two respellings fixes only one axis.
///
/// An occupant's own span never encloses its own junction (`gap == end`, and
/// the interior test is `gap < end`), so a lone `dup`/`repeat` cannot conflict
/// with itself.
///
/// One warning is emitted per same-junction occupant group and one per span
/// edit that encloses ≥1 occupant. Iteration is in deterministic order
/// (BTreeMap junction key, then input index) so equivalent inputs yield
/// identical warning sequences.
///
/// This must run on the *pre-merge* allele members: the normalizer collapses
/// overlapping cis edits into a single combined edit before
/// [`detect_overlap_conflicts`] sees the post-shift list, so by then the
/// overlap is no longer observable as two sub-variants.
pub(crate) fn detect_insertion_overlaps(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<NormalizationWarning> {
    junction_overlaps(variants, phase, true)
}

/// Only the *interior-junction* half of [`detect_insertion_overlaps`]: a
/// junction strictly inside a span edit that keeps the bases it spans.
///
/// This is the subset that is genuinely ambiguous about what the allele
/// denotes, and it is the only subset the #1406 conflict-preservation gates may
/// act on. The same-junction half must stay out of them: two insertions at one
/// junction are a two-member **spelling** of a variant that also has a
/// single-member spelling, and merging them is what makes the two converge.
/// Treating that as unresolvable made `g.[263_264insAC;264_265insAA]` settle
/// apart from `g.264_265insCAAA` — measured, and caught by
/// `cis_spelling_confluence_gap::converged_pairs_stay_converged`, which reports
/// it as "#1301 regressed". Losing confluence to gain a strict-mode property is
/// the wrong trade: confluence is what downstream consumers key on.
pub(crate) fn detect_interior_junction_conflicts(
    variants: &[HgvsVariant],
    phase: AllelePhase,
) -> Vec<NormalizationWarning> {
    junction_overlaps(variants, phase, false)
}

/// Shared body. `include_same_junction` selects branch (a) below.
fn junction_overlaps(
    variants: &[HgvsVariant],
    phase: AllelePhase,
    include_same_junction: bool,
) -> Vec<NormalizationWarning> {
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return Vec::new();
    }
    /// An edit that anchors at a zero-width junction: a true `ins`, or the
    /// copies a `dup`/`repeat` lands 3' of its own span. `kind` is the label the
    /// emitted warning reports for this member (`"ins"` / `"dup"` / `"repeat"`).
    struct Insertion {
        idx: usize,
        accession: String,
        coord_system: &'static str,
        region: Region,
        gap: i64,
        kind: &'static str,
    }
    struct Span {
        idx: usize,
        accession: String,
        coord_system: &'static str,
        region: Region,
        start: i64,
        end: i64,
        /// Whether the edit **removes** every base it spans, leaving nothing
        /// for an interior junction to be positioned against. See branch (b).
        removes_its_bases: bool,
    }

    let mut insertions: Vec<Insertion> = Vec::new();
    let mut spans: Vec<Span> = Vec::new();
    for (idx, variant) in variants.iter().enumerate() {
        let Some((coord_system, region, start, end)) = simple_span(variant) else {
            continue;
        };
        let Some(accession) = variant.accession().map(|a| a.to_string()) else {
            continue;
        };
        let Some(edit) = inner_edit(variant) else {
            continue;
        };
        if matches!(edit, NaEdit::Insertion { .. }) {
            // A canonical insertion anchors at two adjacent positions; anything
            // else (e.g. a malformed single-position insertion) has no junction.
            if end == start + 1 {
                insertions.push(Insertion {
                    idx,
                    accession,
                    coord_system,
                    region,
                    gap: start,
                    kind: "ins",
                });
            }
        } else if is_overlap_edit(edit) {
            // A `dup` and a `repeat` are each both a span (the run they read)
            // and a junction occupant (the copies they write, landing
            // immediately 3' of that run). Registering only the span would miss
            // the conflict whenever the per-member pipeline has respelled an
            // interior `ins` as one of these — which it does routinely, and
            // which axis it picks decides *which* spelling: the CDS axis
            // rewrites `c.5_6insAA` to `c.5_6dup`, the genomic axis rewrites
            // the same shape to `g.1005_1006A[4]`.
            if let Some(kind) = junction_writing_kind(edit) {
                insertions.push(Insertion {
                    idx,
                    accession: accession.clone(),
                    coord_system,
                    region,
                    gap: end,
                    kind,
                });
            }
            spans.push(Span {
                idx,
                accession,
                coord_system,
                region,
                start,
                end,
                removes_its_bases: matches!(edit, NaEdit::Deletion { .. }),
            });
        }
    }

    let mut warnings = Vec::new();

    // (a) Two or more junction-anchored edits sharing a junction. Dups and
    // repeats count alongside true insertions: `[3_6dup;5_6dup]` writes both
    // copies into the slot after position 6 with no defined order between them,
    // exactly as `[5_6insA;5_6insT]` does. Branch (b) does not cover that pair —
    // the inner dup's junction sits at the outer dup's *edge*, not interior to
    // it — so omitting them here would leave the shape unreported even though
    // the sequence-first guard in `merge::…` refuses it (one `after[slot]`, two
    // writers). The same holds for the `repeat` spelling.
    // Each group carries `(index, kind)` so the warning can name every member by
    // its own spelling rather than labelling the whole group `ins`.
    /// `(accession, coordinate system, region, gap)` — one zero-width junction
    /// on one molecule.
    type JunctionKey = (String, &'static str, Region, i64);
    /// The members writing at a junction, as `(index into `variants`, kind)`.
    type Occupants = Vec<(usize, &'static str)>;

    let mut by_junction: BTreeMap<JunctionKey, Occupants> = BTreeMap::new();
    for ins in &insertions {
        by_junction
            .entry((ins.accession.clone(), ins.coord_system, ins.region, ins.gap))
            .or_default()
            .push((ins.idx, ins.kind));
    }
    for ((accession, coord_system, _region, _gap), occupants) in &by_junction {
        if !include_same_junction || occupants.len() < 2 {
            continue;
        }
        // Render the junction via the occupant's HGVS Display (like branch (b)
        // and `detect_overlap_conflicts`) so region prefixes (`*`/`-`) survive;
        // the raw signed `gap` drops them (e.g. 3'UTR `*1_*2` → `1_2`).
        let location = location_for_variant(&variants[occupants[0].0])
            .expect("same-junction occupant has a renderable location");
        warnings.push(NormalizationWarning::OverlapConflict {
            accession: accession.clone(),
            coordinate_system: coord_system.to_string(),
            location,
            edit_kinds: occupants.iter().map(|(_, k)| k.to_string()).collect(),
        });
    }

    // (b) An insertion junction strictly interior to a span edit's range.
    //
    // **Except a pure deletion (#1406).** The conflict this branch reports is
    // that the junction's position within the span is meaningful and the
    // combination therefore has no single answer. A deletion removes every base
    // it spans, so an interior junction has nothing left to be positioned
    // against: `g.[2_3del;2_3insAA]` denotes `AA` in place of `2_3` whichever
    // order the two members are applied in, and whichever interior junction the
    // insertion had named. It composes uniquely, so it was never a conflict —
    // the same argument #1411 used to stop rejecting `g.[24dup;24C>G]`, whose
    // members likewise write to disjoint places.
    //
    // Reporting it anyway had a concrete cost beyond the false verdict: strict
    // mode rejected the input while accepting its own lenient output
    // (`g.2_3delinsAA`), which is the laundering #1406 row 3 is about, and the
    // merged form is what `delins.md:86-89` asks for outright.
    //
    // Deletion only, deliberately. A `delins` does **not** qualify: its payload
    // survives, so an interior insertion has a position relative to it —
    // `g.[2_3delinsGG;2_3insA]` genuinely does not say whether the `A` lands
    // before, inside or after the `GG`. Nor does `inv`, `dup` or `repeat`,
    // each of which keeps the spanned bases and so keeps the interior junction
    // meaningful. Those remain conflicts.
    for span in &spans {
        if span.removes_its_bases {
            continue;
        }
        let interior = insertions.iter().filter(|ins| {
            // A `dup`/`repeat` is registered as both span and occupant, so it
            // meets itself here. The interior test already excludes it — its
            // junction is its own `end`, and `gap < end` is strict — but state
            // the invariant locally rather than leave it resting on that
            // coincidence.
            ins.idx != span.idx
                && ins.accession == span.accession
                && ins.coord_system == span.coord_system
                && ins.region == span.region
                && span.start <= ins.gap
                // `gap + 1 <= end` (junction interior); `gap < end` for ints.
                && ins.gap < span.end
        });
        let mut edit_kinds = vec![edit_kind(&variants[span.idx])
            .expect("span edit has a known kind")
            .to_string()];
        edit_kinds.extend(interior.map(|ins| ins.kind.to_string()));
        if edit_kinds.len() < 2 {
            continue;
        }
        warnings.push(NormalizationWarning::OverlapConflict {
            accession: span.accession.clone(),
            coordinate_system: span.coord_system.to_string(),
            location: location_for_variant(&variants[span.idx])
                .expect("span edit has a renderable location"),
            edit_kinds,
        });
    }

    warnings
}

/// Group key for the coincident-bounds detector.
///
/// `accession` is rendered as the canonical `Accession::Display` string so
/// equality is value-based (two variants pointing at distinct `Arc<str>`
/// instances of the same accession compare equal). `coord_system` is the
/// HGVS prefix character (`g`/`c`/`n`/`r`/`m`) — note `m.` differs from
/// `g.` even though both share `Region::Genome`, so we key on it explicitly
/// rather than collapsing.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct GroupKey {
    accession: String,
    coord_system: &'static str,
    region: Region,
    start: i64,
    end: i64,
}

/// Build the group key for a sub-variant, or `None` if it can't participate
/// in coincident-bounds detection.
///
/// Returns `None` for:
/// - non-NaEdit variants (`Protein`, `RnaFusion`, `Allele`, `NullAllele`,
///   `UnknownAllele`, `Circular`)
/// - uncertain edits (`Mu::Uncertain` / `Mu::Unknown`)
/// - `NaEdit` kinds with no fixed reference span (`Insertion`, `DupIns`,
///   `Identity`, `Conversion`, `Unknown`, `Methylation`, `CopyNumber`,
///   `MultiRepeat`)
/// - positions with intronic offsets, `?` sentinels, or special anchors
///   (`pter`/`qter`/`cen`)
fn group_key(variant: &HgvsVariant) -> Option<GroupKey> {
    let (coord_system, region, start, end) = simple_range(variant)?;
    edit_kind(variant)?;
    let accession = variant.accession()?.to_string();
    Some(GroupKey {
        accession,
        coord_system,
        region,
        start,
        end,
    })
}

/// Extract `(coord_system, region, start, end)` for an overlap-eligible
/// sub-variant. Mirrors [`super::merge::simple_range_for_variant`] but is
/// more permissive on edit kind: `Duplication`, `Inversion`, and
/// (single-base or range) `Repeat` all have a definite reference span and
/// can collide with other edits, so they're included here. `Insertion` is
/// excluded — its anchor is `[end, start]` (zero-width at a boundary) and
/// the spec excludes insertions from this rule.
fn simple_range(variant: &HgvsVariant) -> Option<(&'static str, Region, i64, i64)> {
    match variant {
        HgvsVariant::Genome(g) => {
            na_range(&g.loc_edit, genome_range).map(|(r, s, e)| ("g", r, s, e))
        }
        HgvsVariant::Cds(c) => na_range(&c.loc_edit, cds_range).map(|(r, s, e)| ("c", r, s, e)),
        HgvsVariant::Tx(t) => na_range(&t.loc_edit, tx_range).map(|(r, s, e)| ("n", r, s, e)),
        HgvsVariant::Rna(r) => na_range(&r.loc_edit, rna_range).map(|(rg, s, e)| ("r", rg, s, e)),
        HgvsVariant::Mt(m) => na_range(&m.loc_edit, genome_range).map(|(r, s, e)| ("m", r, s, e)),
        _ => None,
    }
}

fn na_range<L>(
    loc_edit: &LocEdit<Interval<L>, NaEdit>,
    range_fn: impl Fn(&Interval<L>) -> Option<(Region, i64, i64)>,
) -> Option<(Region, i64, i64)> {
    if !loc_edit.edit.is_certain() {
        return None;
    }
    let edit = loc_edit.edit.inner()?;
    if !is_overlap_edit(edit) {
        return None;
    }
    if writes_only_at_a_junction(edit) {
        return None;
    }
    range_fn(&loc_edit.location)
}

/// Whether this edit's *write* lands at a junction rather than over the span it
/// names — in which case coincident bounds are not a conflict.
///
/// Coincident-bounds detection asks whether two members claim the same bases.
/// That is a question about what they **write**, and for most edit kinds the
/// two coincide: a substitution, deletion, delins or inversion replaces the span
/// it names. A `dup` does not. `duplication.md:5` places the copy "directly 3'
/// of the original copy", so `g.23dup` *reads* base 23 and *writes* at the
/// junction 23|24. Against `g.23A>G`, which writes base 23, the two write
/// footprints are disjoint and the composition is unique — `g.22_23insG` — with
/// no ordering to choose. Grouping them by the read span reported a conflict
/// that is not one, and made the verdict depend on the spelling: the same
/// variant written `g.[23A>G;23_24insA]` was accepted.
///
/// The spec says merge rather than reject for this shape. `delins.md:86-89`
/// marks `NM_007294.3:c.[2077G>A;2077_2078insTA]` invalid and gives
/// `c.2077delinsATA` as the correct description, and `general.md:56`
/// prioritisation then requires that form be spelled as a `dup` where it is one.
///
/// **`Repeat` is deliberately not here**, though it also writes at a junction
/// when it grows. Whether it *only* does so depends on whether its unit tiles
/// the reference tract — `A[8]`, `G[6]` and `TA[3]` are indistinguishable at
/// this layer, and `g.[4_9G[6];4_9inv]` rewrites all six bases while satisfying
/// any purely syntactic "does not shorten" test. This function has no reference
/// provider and cannot tell them apart, so a repeat keeps its span. That is the
/// principled line: a `dup`'s write footprint is reference-*independent*, a
/// repeat's is not, which is also why a repeat is registered in *both* detectors
/// and a `dup` in only one.
///
/// An uncertain-extent `dup` keeps its span too — its footprint is not known to
/// be that junction. Same condition `merge::has_same_gap_insertions` uses.
fn writes_only_at_a_junction(edit: &NaEdit) -> bool {
    matches!(
        edit,
        NaEdit::Duplication {
            uncertain_extent: None,
            ..
        }
    )
}

/// The certain inner [`NaEdit`] of an NaEdit-bearing variant, or `None` for
/// non-NaEdit variants and uncertain edits.
fn inner_edit(variant: &HgvsVariant) -> Option<&NaEdit> {
    match variant {
        HgvsVariant::Genome(g) => g.loc_edit.edit.inner(),
        HgvsVariant::Cds(c) => c.loc_edit.edit.inner(),
        HgvsVariant::Tx(t) => t.loc_edit.edit.inner(),
        HgvsVariant::Rna(r) => r.loc_edit.edit.inner(),
        HgvsVariant::Mt(m) => m.loc_edit.edit.inner(),
        _ => None,
    }
}

/// Extract `(coord_system, region, start, end)` for a certain NaEdit-bearing
/// variant *regardless of edit kind*.
///
/// This is [`simple_range`] without the `is_overlap_edit` gate: insertion-
/// overlap detection needs the flanking-position span of an `Insertion`,
/// which `simple_range` deliberately drops. Position-validity gates (special
/// anchors, intronic offsets, `?` sentinels, region splits) still apply via
/// the per-coordinate-system range helpers.
fn simple_span(variant: &HgvsVariant) -> Option<(&'static str, Region, i64, i64)> {
    fn na_span<L>(
        loc_edit: &LocEdit<Interval<L>, NaEdit>,
        range_fn: impl Fn(&Interval<L>) -> Option<(Region, i64, i64)>,
    ) -> Option<(Region, i64, i64)> {
        if !loc_edit.edit.is_certain() {
            return None;
        }
        range_fn(&loc_edit.location)
    }
    match variant {
        HgvsVariant::Genome(g) => {
            na_span(&g.loc_edit, genome_range).map(|(r, s, e)| ("g", r, s, e))
        }
        HgvsVariant::Cds(c) => na_span(&c.loc_edit, cds_range).map(|(r, s, e)| ("c", r, s, e)),
        HgvsVariant::Tx(t) => na_span(&t.loc_edit, tx_range).map(|(r, s, e)| ("n", r, s, e)),
        HgvsVariant::Rna(r) => na_span(&r.loc_edit, rna_range).map(|(rg, s, e)| ("r", rg, s, e)),
        HgvsVariant::Mt(m) => na_span(&m.loc_edit, genome_range).map(|(r, s, e)| ("m", r, s, e)),
        _ => None,
    }
}

/// The occupant label for a span edit that also *writes* at the junction 3' of
/// its own span, or `None` for one that only rewrites the span in place.
///
/// A `dup` copies its run to the slot after its last base; a `repeat` rewrites
/// a tract to a different copy count, and any expansion lands at that same
/// slot. Both are therefore junction occupants as well as spans, which is what
/// lets the detector recognise an interior insertion after the per-member
/// pipeline has respelled it as one of them. Everything else
/// (`sub`/`del`/`delins`/`inv`) edits its span in place and writes no junction.
fn junction_writing_kind(edit: &NaEdit) -> Option<&'static str> {
    match edit {
        NaEdit::Duplication { .. } => Some("dup"),
        NaEdit::Repeat { .. } => Some("repeat"),
        _ => None,
    }
}

/// Edit kinds that have a single, definite reference span. Insertions are
/// deliberately excluded: their anchor is between two bases, so the
/// "coincident bounds" notion does not apply.
fn is_overlap_edit(edit: &NaEdit) -> bool {
    match edit {
        NaEdit::Substitution { .. }
        | NaEdit::SubstitutionNoRef { .. }
        | NaEdit::Deletion { .. }
        | NaEdit::Delins { .. }
        | NaEdit::Duplication { .. }
        | NaEdit::Inversion { .. }
        | NaEdit::Repeat { .. } => true,
        NaEdit::Insertion { .. }
        | NaEdit::BreakpointInsertion { .. }
        | NaEdit::DupIns { .. }
        | NaEdit::MultiRepeat { .. }
        | NaEdit::Identity { .. }
        | NaEdit::Conversion { .. }
        | NaEdit::Unknown { .. }
        | NaEdit::Methylation { .. }
        | NaEdit::CopyNumber { .. }
        | NaEdit::Splice { .. }
        // N-padded deletions sit over an uncertain `(start_end)` range with no
        // definite reference span, so they are not overlap edits.
        | NaEdit::NPaddedDeletion { .. }
        | NaEdit::NoProduct
        | NaEdit::PositionOnly => false,
    }
}

/// Short tag for the edit kind reported in the warning.
fn edit_kind(variant: &HgvsVariant) -> Option<&'static str> {
    let inner: Option<&NaEdit> = match variant {
        HgvsVariant::Genome(g) => g.loc_edit.edit.inner(),
        HgvsVariant::Cds(c) => c.loc_edit.edit.inner(),
        HgvsVariant::Tx(t) => t.loc_edit.edit.inner(),
        HgvsVariant::Rna(r) => r.loc_edit.edit.inner(),
        HgvsVariant::Mt(m) => m.loc_edit.edit.inner(),
        _ => None,
    };
    let edit = inner?;
    Some(match edit {
        NaEdit::Substitution { .. } | NaEdit::SubstitutionNoRef { .. } => "sub",
        NaEdit::Deletion { .. } => "del",
        NaEdit::Delins { .. } => "delins",
        NaEdit::Duplication { .. } => "dup",
        NaEdit::Inversion { .. } => "inv",
        NaEdit::Repeat { .. } => "repeat",
        // Excluded by `is_overlap_edit`; should never reach here for variants
        // that survived `group_key`. Returning None keeps `edit_kind` honest
        // when called from a future caller without the same prefilter.
        _ => return None,
    })
}

/// Render the canonical span text for the warning's `location` field.
///
/// Matches HGVS Display: a point uses bare `100` / `*1` / `-3`, a range
/// uses `start_end`. Returns `None` for non-NaEdit variants.
fn location_for_variant(variant: &HgvsVariant) -> Option<String> {
    match variant {
        HgvsVariant::Genome(g) => Some(format_interval(&g.loc_edit.location)),
        HgvsVariant::Cds(c) => Some(format_interval(&c.loc_edit.location)),
        HgvsVariant::Tx(t) => Some(format_interval(&t.loc_edit.location)),
        HgvsVariant::Rna(r) => Some(format_interval(&r.loc_edit.location)),
        HgvsVariant::Mt(m) => Some(format_interval(&m.loc_edit.location)),
        _ => None,
    }
}

/// Render an interval using the underlying position type's Display impl,
/// collapsing `start == end` to a single point.
fn format_interval<P: std::fmt::Display + PartialEq>(interval: &Interval<P>) -> String {
    let start = render_boundary(&interval.start);
    let end = render_boundary(&interval.end);
    if start == end {
        start
    } else {
        format!("{}_{}", start, end)
    }
}

/// Render a single side of an interval. `simple_range` rejects any
/// variant whose boundary isn't `Single(Certain(_))`, so render_boundary
/// is only ever called on certain boundaries — the other arms exist to
/// document the invariant and panic loudly if it's ever violated.
fn render_boundary<P: std::fmt::Display>(boundary: &UncertainBoundary<P>) -> String {
    match boundary {
        UncertainBoundary::Single(Mu::Certain(p)) => p.to_string(),
        UncertainBoundary::Single(Mu::Uncertain(_))
        | UncertainBoundary::Single(Mu::Unknown)
        | UncertainBoundary::Range { .. } => {
            unreachable!("simple_range gates these out")
        }
    }
}

// --- Per-coordinate-system range extraction ----------------------------------
//
// These helpers mirror the equivalents in `merge.rs`. They are duplicated
// (rather than shared) because `merge.rs`'s versions are intentionally
// scoped to that module — re-exporting them all just to share four
// six-line helpers would widen merge.rs's surface area for no real
// benefit. The shared piece is `Region`, which is `pub(crate)`.

fn genome_range(interval: &Interval<GenomePos>) -> Option<(Region, i64, i64)> {
    let s = simple_genome(interval.start.as_single()?)?;
    let e = simple_genome(interval.end.as_single()?)?;
    Some((Region::Genome, s, e))
}

fn simple_genome(mu: &Mu<GenomePos>) -> Option<i64> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_special() || pos.offset.is_some() {
        return None;
    }
    i64::try_from(pos.base).ok()
}

fn cds_range(interval: &Interval<CdsPos>) -> Option<(Region, i64, i64)> {
    let (rs, s) = simple_cds(interval.start.as_single()?)?;
    let (re, e) = simple_cds(interval.end.as_single()?)?;
    if rs != re {
        return None;
    }
    Some((rs, s, e))
}

fn simple_cds(mu: &Mu<CdsPos>) -> Option<(Region, i64)> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_unknown() || pos.is_intronic() {
        return None;
    }
    if pos.is_3utr() {
        return (pos.base >= 1).then_some((Region::ThreePrimeUtr, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::FivePrimeUtr, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Cds, pos.base));
    }
    None
}

fn tx_range(interval: &Interval<TxPos>) -> Option<(Region, i64, i64)> {
    let (rs, s) = simple_tx(interval.start.as_single()?)?;
    let (re, e) = simple_tx(interval.end.as_single()?)?;
    if rs != re {
        return None;
    }
    Some((rs, s, e))
}

fn simple_tx(mu: &Mu<TxPos>) -> Option<(Region, i64)> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_intronic() {
        return None;
    }
    if pos.is_downstream() {
        return (pos.base >= 1).then_some((Region::TxDownstream, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::TxUpstream, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Tx, pos.base));
    }
    None
}

fn rna_range(interval: &Interval<RnaPos>) -> Option<(Region, i64, i64)> {
    let (rs, s) = simple_rna(interval.start.as_single()?)?;
    let (re, e) = simple_rna(interval.end.as_single()?)?;
    if rs != re {
        return None;
    }
    Some((rs, s, e))
}

fn simple_rna(mu: &Mu<RnaPos>) -> Option<(Region, i64)> {
    let pos = match mu {
        Mu::Certain(p) => p,
        _ => return None,
    };
    if pos.is_intronic() {
        return None;
    }
    if pos.is_3utr() {
        return (pos.base >= 1).then_some((Region::ThreePrimeUtr, pos.base));
    }
    if pos.base < 0 {
        return Some((Region::FivePrimeUtr, pos.base));
    }
    if pos.base > 0 {
        return Some((Region::Cds, pos.base));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

    fn parse_allele(s: &str) -> (Vec<HgvsVariant>, AllelePhase) {
        let v = parse_hgvs(s).expect("parse");
        match v {
            HgvsVariant::Allele(a) => (a.variants, a.phase),
            other => panic!("expected allele in test, got {:?}", other),
        }
    }

    #[test]
    fn same_position_two_subs_emits_one_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100G>A;100A>C]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(
            warnings.len(),
            1,
            "expected one warning, got {:?}",
            warnings
        );
        let NormalizationWarning::OverlapConflict {
            accession,
            coordinate_system,
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(accession, "NC_000001.11");
        assert_eq!(coordinate_system, "g");
        assert_eq!(location, "100");
        assert_eq!(edit_kinds, &vec!["sub".to_string(), "sub".to_string()]);
    }

    #[test]
    fn same_position_sub_plus_del_emits_one_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100del;100A>C]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(warnings.len(), 1);
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!();
        };
        assert!(edit_kinds.contains(&"del".to_string()));
        assert!(edit_kinds.contains(&"sub".to_string()));
    }

    #[test]
    fn coincident_range_del_inv_emits_one_warning() {
        let (variants, phase) = parse_allele("NM_TEST.1:c.[100_103del;100_103inv]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(warnings.len(), 1);
        let NormalizationWarning::OverlapConflict {
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!();
        };
        assert_eq!(location, "100_103");
        assert!(edit_kinds.contains(&"del".to_string()));
        assert!(edit_kinds.contains(&"inv".to_string()));
    }

    #[test]
    fn three_subs_at_one_base_emit_one_group_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100A>C;100A>G;100A>T]");
        let warnings = detect_overlap_conflicts(&variants, phase);
        assert_eq!(warnings.len(), 1, "groups, not pairs");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!();
        };
        assert_eq!(edit_kinds.len(), 3);
    }

    #[test]
    fn adjacent_subs_no_warning() {
        let (variants, phase) = parse_allele("NC_000001.11:g.[100A>C;101A>G]");
        assert!(detect_overlap_conflicts(&variants, phase).is_empty());
    }

    #[test]
    fn multi_accession_no_warning() {
        let (variants, phase) = parse_allele("[NC_000001.11:g.100A>C;NC_000002.11:g.100A>G]");
        assert!(detect_overlap_conflicts(&variants, phase).is_empty());
    }

    #[test]
    fn insertion_at_boundary_no_warning() {
        // Insertions anchor at p_p+1 — no single-base location to coincide.
        let (variants, phase) = parse_allele("NC_000001.11:g.[100A>C;100_101insT]");
        assert!(detect_overlap_conflicts(&variants, phase).is_empty());
    }

    // --- Insertion overlap detection (#486 EOVERLAP) -------------------------

    #[test]
    fn two_insertions_at_same_junction_emit_one_warning() {
        // Two insertions at the same interspace `4_5` overlap (mutalyzer
        // EOVERLAP). The inserted sequence is irrelevant — the junction is
        // shared.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_5insT;4_5insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict {
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(location, "4_5");
        assert_eq!(edit_kinds, &vec!["ins".to_string(), "ins".to_string()]);
    }

    #[test]
    fn same_junction_insertions_render_utr_location_with_star_prefix() {
        // Regression: the same-junction branch must render the location via
        // HGVS Display, not the raw signed base. A 3'UTR junction `*1_*2`
        // would otherwise drop the `*` prefix and print `1_2`.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[*1_*2insT;*1_*2insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { location, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(location, "*1_*2");
    }

    #[test]
    fn insertion_interior_to_delins_emits_one_warning() {
        // An insertion `274_275ins` whose junction sits strictly inside a
        // `274_275delins` range overlaps it (mutalyzer EOVERLAP).
        let (variants, phase) = parse_allele("NM_TEST.1:c.[274_275delinsT;274_275insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert!(edit_kinds.contains(&"delins".to_string()));
        assert!(edit_kinds.contains(&"ins".to_string()));
    }

    #[test]
    fn insertion_interior_to_deletion_is_not_a_conflict() {
        // Insertion junction `5_6` is strictly interior to the deleted range
        // `4_7`, but a deletion removes every base it spans, so the insertion
        // has nothing left to be positioned against: the pair denotes the same
        // bases whichever order it is applied in, and whichever interior
        // junction the insertion had named. Not a conflict (#1406).
        //
        // `insertion_interior_to_delins_emits_one_warning` is the
        // discriminating sibling: a `delins` keeps a payload, so an
        // interior insertion *does* have a position relative to it. Only a
        // pure deletion is exempt.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_7del;5_6insAA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert!(
            warnings.is_empty(),
            "an insertion interior to a pure deletion composes uniquely, so it \
             must not be reported; got {warnings:?}"
        );
    }

    #[test]
    fn insertions_sharing_junction_and_interior_to_span_emit_two_warnings() {
        // `5_6insA` and `5_6insT` share junction `5` (branch (a)) *and* both
        // sit strictly interior to the `4_7delinsGG` range positions 4..=7
        // (branch (b)), so the pass emits two `OverlapConflict` warnings for
        // the one overlap cluster. The combined-shape outcome is still a
        // correct rejection; this pins the per-branch warning count so a
        // future dedup refactor doesn't silently change it.
        //
        // The span is a `delins` rather than the `4_7del` this used to use:
        // a pure deletion no longer contributes a branch-(b) warning (#1406),
        // so with `del` here only branch (a) would fire and the test would no
        // longer exercise the two-branch shape it exists for.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_7delinsGG;5_6insA;5_6insT]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 2, "expected two warnings, got {warnings:?}");
        assert!(
            warnings
                .iter()
                .all(|w| matches!(w, NormalizationWarning::OverlapConflict { .. })),
            "all warnings must be OverlapConflict: {warnings:?}"
        );
    }

    #[test]
    fn duplication_interior_to_inversion_emits_one_warning() {
        // A `dup` lands its copy at the junction after its own last base, so
        // `5_6dup` occupies junction 6, strictly interior to `5_9inv`
        // (positions 5..=9). This is the `ins`-spelled conflict
        // `[4_10inv;5_6insAA]` after the per-member pipeline has respelled it,
        // and it must be recognised as the same conflict — otherwise
        // normalization is not idempotent on the conflicting allele.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[5_9inv;5_6dup]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict {
            location,
            edit_kinds,
            ..
        } = &warnings[0]
        else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(location, "5_9", "the warning names the enclosing span");
        assert_eq!(edit_kinds, &vec!["inv".to_string(), "dup".to_string()]);
    }

    #[test]
    fn duplication_conflict_is_detected_in_either_member_order() {
        // Detection must not depend on authored order: the sort that renders
        // cis members in genomic order is gated *off* by this very warning, so
        // both spellings have to be reported or the two orders normalize to
        // each other and idempotency breaks.
        for allele in ["NM_TEST.1:c.[5_9inv;5_6dup]", "NM_TEST.1:c.[5_6dup;5_9inv]"] {
            let (variants, phase) = parse_allele(allele);
            assert_eq!(
                detect_insertion_overlaps(&variants, phase).len(),
                1,
                "expected one warning for {allele}"
            );
        }
    }

    #[test]
    fn lone_duplication_does_not_conflict_with_itself() {
        // A `dup` is registered as both a span and a junction occupant. Its own
        // junction (`gap == end`) is not strictly interior to its own span
        // (`gap < end` fails), and the self-pairing is skipped explicitly, so a
        // duplication alongside an unrelated edit stays clean.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[5_6dup;20A>G]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn duplication_abutting_a_span_edge_no_warning() {
        // `1_4dup` writes at junction 4, which is the *edge* of `5_9inv`
        // (interior junctions are 5..=8), not interior to it. Mirrors
        // `insertions_flanking_a_sub_no_warning` for the dup spelling.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[1_4dup;5_9inv]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn two_duplications_at_one_junction_emit_one_warning() {
        // Both dups write their copy into the slot after position 6, with no
        // defined order between the two payloads — the same ambiguity branch
        // (a) reports for `[5_6insA;5_6insT]`. Branch (b) does *not* cover this
        // pair: the inner dup's junction (6) is the outer dup's end, i.e. its
        // edge, and the interior test is strict (`gap < end`).
        let (variants, phase) = parse_allele("NM_TEST.1:c.[3_6dup;5_6dup]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(edit_kinds, &vec!["dup".to_string(), "dup".to_string()]);
    }

    #[test]
    fn repeat_interior_to_inversion_emits_one_warning() {
        // The genomic axis respells an interior insertion as a repeat
        // expansion rather than a dup, so the `repeat` spelling has to be
        // recognised as a junction occupant for the same reason: `1005_1006A[4]`
        // writes its extra copies at junction 1006, interior to `1005_1009inv`.
        let (variants, phase) = parse_allele("NC_000001.11:g.[1005_1009inv;1005_1006A[4]]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(edit_kinds, &vec!["inv".to_string(), "repeat".to_string()]);
    }

    #[test]
    fn duplication_and_insertion_at_one_junction_emit_one_warning() {
        // Mixed spellings compete for the same slot too, and the warning must
        // name each member by its own kind rather than labelling both `ins`.
        let (variants, phase) = parse_allele("NM_TEST.1:c.[5_6dup;6_7insA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
        let NormalizationWarning::OverlapConflict { edit_kinds, .. } = &warnings[0] else {
            panic!("expected OverlapConflict, got {:?}", warnings[0]);
        };
        assert_eq!(edit_kinds, &vec!["dup".to_string(), "ins".to_string()]);
    }

    #[test]
    fn two_insertions_at_different_junctions_no_warning() {
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_5insT;8_9insA]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn insertions_flanking_a_sub_no_warning() {
        // `273_274ins` and `274_275ins` are at distinct junctions either side
        // of the substitution at `274`; none overlaps the single-base sub.
        // Mutalyzer normalizes this to `c.274delinsTTA` (no EOVERLAP).
        let (variants, phase) = parse_allele("NM_TEST.1:c.[273_274insT;274G>T;274_275insA]");
        assert!(
            detect_insertion_overlaps(&variants, phase).is_empty(),
            "non-overlapping flanking insertions must not warn: {:?}",
            detect_insertion_overlaps(&variants, phase)
        );
    }

    #[test]
    fn insertion_overlap_only_in_cis() {
        // Trans phase: the edits are on different haplotypes, so coincident
        // junctions are not a conflict.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_5insT];[4_5insA]");
        assert!(detect_insertion_overlaps(&variants, phase).is_empty());
    }

    #[test]
    fn end_to_end_normalize_emits_warning() {
        use crate::normalize::Normalizer;
        use crate::reference::mock::MockProvider;

        let normalizer = Normalizer::new(MockProvider::new());
        let v = parse_hgvs("NC_000001.11:g.[100G>A;100A>C]").expect("parse");
        let result = normalizer
            .normalize_with_diagnostics(&v)
            .expect("normalize");
        assert!(
            result
                .warnings
                .iter()
                .any(|w| w.code() == "OVERLAP_CONFLICTING_EDITS"),
            "expected OVERLAP_CONFLICTING_EDITS in warnings, got {:?}",
            result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
        );
        let out = result.result.to_string();
        assert!(
            out.contains("100G>A") && out.contains("100A>C"),
            "expected pass-through, got {out}"
        );
    }
}
