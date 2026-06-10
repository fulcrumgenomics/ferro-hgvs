//! Detect cis-allele sub-variants whose reference bounds coincide.
//!
//! See `docs/superpowers/specs/2026-05-05-allele-overlap-same-position-design.md`.
//!
//! Two or more cis-phase, same-accession sub-variants whose `(coord-system
//! region, signed start, signed end)` keys are identical have no canonical
//! HGVS form. Rather than silently picking a winner, ferro preserves the
//! input verbatim and emits one [`NormalizationWarning::OverlapConflict`]
//! per coincident-bounds group. Insertions are excluded (they anchor at a
//! boundary `[end, start]`, not a single-base location).

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
