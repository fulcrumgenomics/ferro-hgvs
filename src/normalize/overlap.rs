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
//! - [`detect_insertion_overlaps`] — *pre-merge* insertion overlaps: two
//!   insertions at one junction, or an insertion interior to a span edit
//!   (mutalyzer `EOVERLAP`, #486). Must run before the normalizer's merge step
//!   collapses overlapping cis edits into one combined edit.

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

/// Detect overlaps that involve at least one insertion within a cis allele.
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
/// One warning is emitted per same-junction insertion group and one per span
/// edit that encloses ≥1 insertion. Iteration is in deterministic order
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
    if phase != AllelePhase::Cis || variants.len() < 2 {
        return Vec::new();
    }
    struct Insertion {
        idx: usize,
        accession: String,
        coord_system: &'static str,
        region: Region,
        gap: i64,
    }
    struct Span {
        idx: usize,
        accession: String,
        coord_system: &'static str,
        region: Region,
        start: i64,
        end: i64,
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
                });
            }
        } else if is_overlap_edit(edit) {
            spans.push(Span {
                idx,
                accession,
                coord_system,
                region,
                start,
                end,
            });
        }
    }

    let mut warnings = Vec::new();

    // (a) Two or more insertions sharing a junction.
    let mut by_junction: BTreeMap<(String, &'static str, Region, i64), Vec<usize>> =
        BTreeMap::new();
    for ins in &insertions {
        by_junction
            .entry((ins.accession.clone(), ins.coord_system, ins.region, ins.gap))
            .or_default()
            .push(ins.idx);
    }
    for ((accession, coord_system, _region, _gap), indices) in &by_junction {
        if indices.len() < 2 {
            continue;
        }
        // Render the junction via the insertion's HGVS Display (like branch (b)
        // and `detect_overlap_conflicts`) so region prefixes (`*`/`-`) survive;
        // the raw signed `gap` drops them (e.g. 3'UTR `*1_*2` → `1_2`).
        let location = location_for_variant(&variants[indices[0]])
            .expect("same-junction insertion has a renderable location");
        warnings.push(NormalizationWarning::OverlapConflict {
            accession: accession.clone(),
            coordinate_system: coord_system.to_string(),
            location,
            edit_kinds: indices.iter().map(|_| "ins".to_string()).collect(),
        });
    }

    // (b) An insertion junction strictly interior to a span edit's range.
    for span in &spans {
        let interior = insertions.iter().filter(|ins| {
            ins.accession == span.accession
                && ins.coord_system == span.coord_system
                && ins.region == span.region
                && span.start <= ins.gap
                // `gap + 1 <= end` (junction interior); `gap < end` for ints.
                && ins.gap < span.end
        });
        let mut edit_kinds = vec![edit_kind(&variants[span.idx])
            .expect("span edit has a known kind")
            .to_string()];
        edit_kinds.extend(interior.map(|_| "ins".to_string()));
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
    range_fn(&loc_edit.location)
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
    fn insertion_interior_to_deletion_emits_one_warning() {
        // Insertion junction `5_6` is strictly interior to the deleted range
        // `4_7` (positions 4..=7), so they overlap.
        let (variants, phase) = parse_allele("NG_012337.1:g.[4_7del;5_6insAA]");
        let warnings = detect_insertion_overlaps(&variants, phase);
        assert_eq!(warnings.len(), 1, "expected one warning, got {warnings:?}");
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
