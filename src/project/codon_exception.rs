//! The coding one-amino-acid exception, evaluated per rendered axis (#1664).
//!
//! `DNA/delins.md:18` carves an exception out of the separation rule:
//!
//! > **exception**: two variants separated by one nucleotide, together affecting
//! > one amino acid, should be described as a "delins".
//!
//! The reading is settled for this project —
//! `rulings[delins-codon-carve-out-gap-one]` is `decided`, `delins.md:18`
//! governing — and the normalizer implements it, in
//! `merge::apply_coding_codon_exception`, gated on `AxisFrame::reading_frame`.
//! That gate asks whether *the axis a description is written on* carries a
//! reading frame. For a bare description that is the whole question. For a
//! projection it is half of one, because the axes are not independent: a genomic
//! input is partitioned by a frameless normalization, and the `c.`/`n.`/`r.`/`p.`
//! axes derived from it inherit a partition the frame would not have chosen. One
//! variant then renders as one member or as two on the *coding* axis depending
//! only on which spelling the caller happened to hand in.
//!
//! This module holds the shape test the projector needs to re-decide that
//! against the transcript, plus the member count it compares either side of the
//! attempt, and nothing else: the merge itself is the normalizer's, performed by
//! re-normalizing the coding axis, which has a frame already.
//!
//! # What is deliberately *not* changed: any `g.` axis, bare or derived
//!
//! Closed issue #79 scoped the exception out of the bare genomic axis, on the
//! grounds that a `g.` description names no transcript and so offers no frame
//! to consult. Nothing here reopens that, and
//! `rulings[projection-codon-exception-is-decided-by-the-rendered-axis]`
//! (`decided`) gives a derived genomic axis the same answer for a reason of its
//! own: `:18` is a conditional whose second conjunct — "together affecting one
//! amino acid" — cannot be stated on a genomic reference, so the exception does
//! not fire there and `DNA/delins.md:17` governs unopposed, asking for the
//! individual descriptions. The predicates below are only ever reached from
//! [`crate::project::VariantProjector`], which always has a transcript, and they
//! re-decide the transcript axes only. A genomic input's own genomic axis is
//! still returned as that input normalizes on its own, which is #79's answer.

use crate::hgvs::edit::{Base, NaEdit};
use crate::hgvs::interval::UncertainBoundary;
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{AllelePhase, HgvsVariant};
use crate::normalize::merge::same_codon;
use crate::project::VariantProjection;

/// The number of top-level members a description carries: 1 for a single
/// variant, the member count for an allele.
pub(crate) fn member_count(variant: &HgvsVariant) -> usize {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.len(),
        _ => 1,
    }
}

/// Carry the parts of the pre-merge projection that the re-projection must not
/// re-derive onto `reprojected` — and **only** those parts.
///
/// The merge re-projects from the merged *coding* form, so every axis the
/// re-projection derives (`c.`, `n.`, `r.`, `p.`) is the merged form's own
/// answer and must be reported as such, decline reason included. Two things are
/// not the merged form's to answer:
///
/// * **the genomic axis.** For a genomic input it is that input's own
///   normalization — #79's answer, and the one
///   `rulings[projection-codon-exception-is-decided-by-the-rendered-axis]`
///   keeps — so both the axis and its decline reason come from before the
///   merge.
/// * **the normalization warnings**, which are the caller's input's (an
///   auto-corrected reference mismatch, say). Re-normalizing an
///   already-corrected description does not re-emit them.
///
/// **A decline reason belongs to the axis whose value it explains.** Assigning
/// the whole [`AxisDeclineReasons`] across — which is what this function
/// replaced — leaves `coding`/`noncoding`/`protein`/`rna` explaining values
/// that no longer exist, so a caller can read "this axis declined because …"
/// beside a rendered axis, or a stale reason beside a genuinely absent one.
/// `AxisDeclineReasons` is per-axis precisely so the two halves cannot be
/// swapped wholesale; pinned by
/// `only_the_genomic_decline_reason_survives_the_re_projection`.
///
/// [`AxisDeclineReasons`]: crate::project::AxisDeclineReasons
pub(crate) fn carry_pre_merge_state(
    reprojected: &mut VariantProjection,
    pre_merge: &mut VariantProjection,
) {
    reprojected.genomic = pre_merge.genomic.take();
    reprojected.axis_decline_reasons.genomic = pre_merge.axis_decline_reasons.genomic.take();
    reprojected.normalization_warnings = std::mem::take(&mut pre_merge.normalization_warnings);
}

/// The first adjacent member pair on `coding` that `delins.md:18` keeps
/// together, as its two CDS positions.
///
/// The shape is the one `merge::apply_coding_codon_exception` matches — two
/// lone substitutions with exactly one unchanged base between them — with the
/// codon test run on the CDS axis, where `same_codon` is meaningful. The
/// `codon-carve-out-shape-restriction` ruling widens the exception past that
/// shape in principle but records that the widening is *not implemented*; this
/// predicate deliberately matches what ships rather than what is licensed, so
/// the two halves of the rule cannot drift apart here.
///
/// **The caller still has to check that the pair lies inside one exon.** A codon
/// may be split across an exon/exon junction, and then the pair is "separated by
/// one nucleotide" on the `c.` axis while being separated by a whole intron on
/// the sequence. Merging those is a different question from the one
/// `delins.md:18` settles, and the repository already answers it the other way
/// for the same shape — `project_cis_intron_split_codon_combines_to_single_missense`
/// pins `c.[4G>C;6T>A]` staying split across a ~100 bp intron and combining only
/// at the protein level — so that decision is left where it already is.
///
/// The pair is returned rather than a `bool` so the caller can run that test
/// only when there is something to test: a transcript lookup on every projection
/// is a cost this rule has no right to impose on descriptions it never touches.
///
/// Answering `Some` says only that the coding axis has a partition the frame
/// would not have chosen. The projector re-normalizes the coding allele to get
/// the merged form, so the merge itself stays in one place.
pub(crate) fn coding_codon_pair(coding: &HgvsVariant) -> Option<(i64, i64)> {
    let HgvsVariant::Allele(allele) = coding else {
        return None;
    };
    if allele.uncertain || allele.phase != AllelePhase::Cis {
        return None;
    }
    allele.variants.windows(2).find_map(|pair| {
        match (cds_substitution(&pair[0]), cds_substitution(&pair[1])) {
            (Some((left, _)), Some((right, _))) if right == left + 2 && same_codon(left, right) => {
                Some((left, right))
            }
            _ => None,
        }
    })
}

/// The position a boundary names, **only** when it is stated with certainty.
///
/// `UncertainBoundary::inner` deliberately answers for `(pos)` as well as for
/// `pos`, which is the wrong reading here: everything this module does ends in
/// rewriting two members as one certain `delins`, so admitting an uncertain
/// endpoint would silently promote `(4)` to `4`. Uncertainty is information the
/// author supplied, and the codon exception is not a licence to drop it — a
/// description that carries any is left exactly as normalization produced it.
fn certain_boundary<T>(boundary: &UncertainBoundary<T>) -> Option<&T> {
    match boundary {
        UncertainBoundary::Single(Mu::Certain(position)) => Some(position),
        _ => None,
    }
}

/// The edit a `Mu` names, only when it is stated with certainty. The companion
/// of [`certain_boundary`], and for the same reason: `Mu::inner` answers for
/// `Mu::Uncertain` too.
fn certain_edit(edit: &Mu<NaEdit>) -> Option<&NaEdit> {
    match edit {
        Mu::Certain(edit) => Some(edit),
        _ => None,
    }
}

/// `(cds base, alternative base)` for a lone single-base `c.` substitution
/// sitting plainly in the CDS.
///
/// `None` for an intronic offset, a `*`-numbered 3'UTR position, a `pter`-class
/// marker or a 5'UTR (negative) base: `same_codon` is only meaningful on the
/// positive CDS body, and the exception is scoped to the coding sequence.
fn cds_substitution(variant: &HgvsVariant) -> Option<(i64, Base)> {
    let HgvsVariant::Cds(cds) = variant else {
        return None;
    };
    let start = certain_boundary(&cds.loc_edit.location.start)?;
    let end = certain_boundary(&cds.loc_edit.location.end)?;
    if start != end {
        return None;
    }
    if start.offset.is_some() || start.utr3 || start.special.is_some() || start.base < 1 {
        return None;
    }
    match certain_edit(&cds.loc_edit.edit)? {
        NaEdit::Substitution { alternative, .. } => Some((start.base, *alternative)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::NormalizationWarning;
    use crate::parse_hgvs;
    use crate::project::AxisDeclineReasons;

    /// Replace every member's edit `Mu` in an allele, so an uncertain spelling
    /// can be built without depending on the parser accepting one.
    fn make_edits_uncertain(variant: &HgvsVariant) -> HgvsVariant {
        let mut variant = variant.clone();
        let HgvsVariant::Allele(allele) = &mut variant else {
            panic!("not an allele");
        };
        for member in &mut allele.variants {
            let HgvsVariant::Cds(cds) = member else {
                panic!("not a c. member");
            };
            let edit = cds.loc_edit.edit.inner().expect("member states an edit");
            cds.loc_edit.edit = Mu::Uncertain(edit.clone());
        }
        variant
    }

    fn pair(description: &str) -> Option<(i64, i64)> {
        coding_codon_pair(&parse_hgvs(description).expect("parses"))
    }

    #[test]
    fn a_gap_of_one_inside_one_codon_is_the_shape_the_exception_keeps_together() {
        // `c.4` and `c.6` are both codon 2, separated by the single unchanged
        // `c.5` — `DNA/delins.md:18`'s "two variants separated by one
        // nucleotide, together affecting one amino acid".
        assert_eq!(pair("NM_000000.1:c.[4C>A;6G>A]"), Some((4, 6)));
    }

    #[test]
    fn a_gap_of_two_is_left_split_because_general_md_34_asks_for_it() {
        // `general.md:34` — "two variants separated by one or more nucleotides
        // should be described individually" — governs everything the exception
        // does not carve out, and the carve-out is a gap of exactly one.
        assert_eq!(pair("NM_000000.1:c.[4C>A;7G>A]"), None);
    }

    #[test]
    fn a_gap_of_one_that_straddles_two_codons_is_not_one_amino_acid() {
        // `c.5` and `c.7` are codon 2 and codon 3, so the pair affects two amino
        // acids and the exception's precondition fails even though the
        // separation is right.
        assert_eq!(pair("NM_000000.1:c.[5C>A;7G>A]"), None);
    }

    #[test]
    fn an_intronic_or_utr_member_is_not_on_the_coding_body() {
        // `same_codon` is only meaningful on the positive CDS body, so an
        // intronic offset, a `*`-numbered 3'UTR position and a negative 5'UTR
        // position each decline rather than being counted into a codon.
        assert_eq!(pair("NM_000000.1:c.[4+1C>A;6G>A]"), None);
        assert_eq!(pair("NM_000000.1:c.[*4C>A;*6G>A]"), None);
        assert_eq!(pair("NM_000000.1:c.[-6C>A;-4G>A]"), None);
    }

    #[test]
    fn a_trans_allele_is_not_two_variants_on_one_sequence() {
        // The exception merges two changes that lie on the same molecule; `;`
        // in a trans allele asserts they do not.
        assert_eq!(pair("NM_000000.1:c.[4C>A];[6G>A]"), None);
    }

    /// A projection whose every axis carries a decline reason naming `tag`, so
    /// that "which side did this field come from" is answerable field by field
    /// and no assertion below can pass by two sides happening to agree.
    fn projection_tagged(tag: &str, genomic: Option<HgvsVariant>) -> VariantProjection {
        VariantProjection {
            genomic,
            coding: None,
            noncoding: None,
            protein: None,
            rna: None,
            transcript_id: "LRG_199t1".to_string(),
            gene_symbol: None,
            is_frameshift: false,
            is_intronic: false,
            is_utr: false,
            affects_init: false,
            normalization_warnings: vec![NormalizationWarning::OverlapConflict {
                accession: tag.to_string(),
                coordinate_system: "c".to_string(),
                location: "4_6".to_string(),
                edit_kinds: vec!["sub".to_string(), "sub".to_string()],
            }],
            axis_decline_reasons: AxisDeclineReasons {
                genomic: Some(format!("{tag} g.")),
                coding: Some(format!("{tag} c.")),
                noncoding: Some(format!("{tag} n.")),
                protein: Some(format!("{tag} p.")),
                rna: Some(format!("{tag} r.")),
            },
        }
    }

    /// A decline reason explains one axis's absence, so it may only travel with
    /// that axis. The merge carries the genomic axis across from before the
    /// merge and lets the re-projection answer for `c.`/`n.`/`r.`/`p.`; carrying
    /// the whole `AxisDeclineReasons` across — which is what this pins against —
    /// leaves those four explaining values that are no longer there.
    #[test]
    fn only_the_genomic_decline_reason_survives_the_re_projection() {
        let genomic = parse_hgvs("LRG_199:g.[499798A>T;499800G>T]").expect("parses");
        let mut pre_merge = projection_tagged("pre-merge", Some(genomic.clone()));
        let mut reprojected = projection_tagged("re-projected", None);

        carry_pre_merge_state(&mut reprojected, &mut pre_merge);

        // The genomic axis is the pre-merge one, so its explanation is too.
        assert_eq!(
            reprojected.genomic.as_ref().map(HgvsVariant::to_string),
            Some(genomic.to_string())
        );
        assert_eq!(
            reprojected.axis_decline_reasons.genomic.as_deref(),
            Some("pre-merge g.")
        );
        // The other four axes are the merged form's own, and so are theirs.
        assert_eq!(
            [
                reprojected.axis_decline_reasons.coding.as_deref(),
                reprojected.axis_decline_reasons.noncoding.as_deref(),
                reprojected.axis_decline_reasons.protein.as_deref(),
                reprojected.axis_decline_reasons.rna.as_deref(),
            ],
            [
                Some("re-projected c."),
                Some("re-projected n."),
                Some("re-projected p."),
                Some("re-projected r."),
            ],
            "a re-derived axis must be explained by the re-projection that derived it"
        );
        // The warnings are the caller's input's, not the merged form's.
        assert!(
            format!("{:?}", reprojected.normalization_warnings).contains("pre-merge"),
            "the input's normalization warnings must survive the re-projection"
        );
    }

    #[test]
    fn an_uncertain_edit_is_left_exactly_as_the_author_wrote_it() {
        // The merge always produces a *certain* `delins`, so admitting an
        // uncertain member would promote `(4C>A)` to `4C>A` — dropping
        // information the author supplied. `Mu::inner` answers for both, which
        // is why these predicates use `certain_edit` instead.
        let certain = parse_hgvs("NM_000000.1:c.[4C>A;6G>A]").expect("parses");
        assert_eq!(coding_codon_pair(&certain), Some((4, 6)));
        assert_eq!(coding_codon_pair(&make_edits_uncertain(&certain)), None);
    }
}
