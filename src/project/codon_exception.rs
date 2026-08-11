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
//! reading frame, which is the right question for a bare description and the
//! wrong one for a projection: a projection names a transcript, so the frame is
//! in hand on every axis it renders. Deciding the partition once, on whichever
//! axis the caller happened to author, and letting the other axes inherit it is
//! what makes `LRG_199t1:c.145_147delinsTGG` render as
//! `LRG_199:g.[494841C>T;494843C>G]` — the exact string `DNA/delins.md:42`
//! renders `class="invalid"` and calls "not correct".
//!
//! This module holds the two shape tests the projector needs to re-decide it
//! against the transcript, and nothing else: the merge itself is either the
//! normalizer's (on the coding axis, which has a frame already) or a
//! three-base rewrite of a pair the normalizer has just split (on the genomic
//! axis, which does not).
//!
//! # What is deliberately *not* changed: the bare `g.` axis
//!
//! Closed issue #79 scoped the exception out of the bare genomic axis, on the
//! grounds that a `g.` description names no transcript and so offers no frame
//! to consult. Nothing here reopens that. The predicates below are only ever
//! reached from [`crate::project::VariantProjector`], which always has a
//! transcript, and the projector applies the genomic half only to a genomic
//! axis it *derived* — a genomic input's own axis is still returned as that
//! input normalizes on its own, which is #79's answer.

use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::UncertainBoundary;
use crate::hgvs::uncertainty::Mu;
use crate::hgvs::variant::{Accession, AllelePhase, HgvsVariant};
use crate::normalize::merge::same_codon;

/// The number of top-level members a description carries: 1 for a single
/// variant, the member count for an allele.
pub(crate) fn member_count(variant: &HgvsVariant) -> usize {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.len(),
        _ => 1,
    }
}

/// Is `variant` a genomic description — a `g.`/`m.` variant, or a cis allele of
/// them?
///
/// Used by the projector to tell a *derived* genomic axis from one that is the
/// input's own normalization. Only the former is re-partitioned; see the module
/// doc on #79.
pub(crate) fn is_genomic_description(variant: &HgvsVariant) -> bool {
    match variant {
        HgvsVariant::Genome(_) | HgvsVariant::Mt(_) => true,
        HgvsVariant::Allele(allele) => allele
            .variants
            .first()
            .is_some_and(|first| matches!(first, HgvsVariant::Genome(_) | HgvsVariant::Mt(_))),
        _ => false,
    }
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

/// Re-merge a genomic axis the frameless re-derivation split, when the coding
/// axis — which does carry the frame — kept the same variant whole.
///
/// Returns `None` whenever the exception does not apply, which is the common
/// case; the caller then leaves the genomic axis exactly as normalization
/// produced it.
///
/// The precondition is read off `coding`: a lone `delins` occupying exactly one
/// codon *is* "two variants separated by one nucleotide, together affecting one
/// amino acid" once the genomic axis has split it into that pair. Reading it
/// there rather than re-mapping the genomic positions back through cdot is what
/// keeps this correct in an `NG_`/`LRG_` frame, whose coordinates are not the
/// chromosome coordinates cdot is placed against.
///
/// # Why the payload is not taken from the coding axis
///
/// The merged genomic edit replaces three bases with `[alt, unchanged, alt]`,
/// and the two flanking bases are the substitutions themselves. Only the middle
/// one has to be looked up, and `middle_base` supplies it from the genomic
/// reference. Transliterating the coding payload instead would need the
/// transcript's orientation *relative to the accession the genomic axis is
/// written on*, which is not `Transcript::strand`: `LRG_199t1` is a minus-strand
/// transcript on `NC_000023.11` and yet ascends on `LRG_199`, because an LRG
/// record is stored in its own gene's orientation. Reading the one base that is
/// actually unknown sidesteps the question.
pub(crate) fn merge_genomic_codon_split(
    genomic: &HgvsVariant,
    coding: &HgvsVariant,
    middle_base: impl Fn(&Accession, u64) -> Option<Base>,
) -> Option<HgvsVariant> {
    codon_delins_replacement(coding)?;
    let HgvsVariant::Allele(allele) = genomic else {
        return None;
    };
    if allele.uncertain || allele.phase != AllelePhase::Cis || allele.variants.len() != 2 {
        return None;
    }
    let (left_pos, left_alt) = genome_substitution(&allele.variants[0])?;
    let (right_pos, right_alt) = genome_substitution(&allele.variants[1])?;
    if right_pos != left_pos + 2 {
        return None;
    }
    let HgvsVariant::Genome(left) = &allele.variants[0] else {
        return None;
    };
    let HgvsVariant::Genome(right) = &allele.variants[1] else {
        return None;
    };
    if left.accession != right.accession {
        return None;
    }
    let middle = middle_base(&left.accession, left_pos + 1)?;
    let bases = vec![left_alt, middle, right_alt];

    let mut merged = left.clone();
    merged.loc_edit.location.start = UncertainBoundary::certain(position(left_pos));
    merged.loc_edit.location.end = UncertainBoundary::certain(position(right_pos));
    merged.loc_edit.edit = Mu::Certain(NaEdit::Delins {
        sequence: InsertedSequence::Literal(Sequence::new(bases)),
        deleted: None,
        deleted_length: None,
        substitution_reference: None,
    });
    Some(HgvsVariant::Genome(merged))
}

/// A plain genomic position with no offset and no `pter`/`qter`/`cen` marker.
fn position(base: u64) -> crate::hgvs::location::GenomePos {
    crate::hgvs::location::GenomePos {
        base,
        special: None,
        offset: None,
    }
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

/// `(position, alternative base)` for a lone single-base `g.` substitution.
fn genome_substitution(variant: &HgvsVariant) -> Option<(u64, Base)> {
    let HgvsVariant::Genome(genome) = variant else {
        return None;
    };
    let start = certain_boundary(&genome.loc_edit.location.start)?;
    let end = certain_boundary(&genome.loc_edit.location.end)?;
    if start != end || start.offset.is_some() || start.special.is_some() {
        return None;
    }
    match certain_edit(&genome.loc_edit.edit)? {
        NaEdit::Substitution { alternative, .. } => Some((start.base, *alternative)),
        _ => None,
    }
}

/// Does `coding` name a lone `c.` `delins` occupying exactly one codon?
///
/// This is the precondition of `delins.md:18` read on the axis that can answer
/// it: a three-position span inside one codon "affects one amino acid", and the
/// genomic split of it is "two variants separated by one nucleotide".
fn codon_delins_replacement(coding: &HgvsVariant) -> Option<()> {
    let HgvsVariant::Cds(cds) = coding else {
        return None;
    };
    let start = certain_boundary(&cds.loc_edit.location.start)?;
    let end = certain_boundary(&cds.loc_edit.location.end)?;
    if start.offset.is_some() || start.utr3 || start.special.is_some() || start.base < 1 {
        return None;
    }
    if end.offset.is_some() || end.utr3 || end.special.is_some() {
        return None;
    }
    if end.base != start.base + 2 || !same_codon(start.base, end.base) {
        return None;
    }
    match certain_edit(&cds.loc_edit.edit)? {
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(sequence),
            ..
        } if sequence.len() == 3 => Some(()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;

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

    #[test]
    fn a_lone_codon_delins_is_the_genomic_halfs_precondition() {
        let whole = parse_hgvs("NM_000000.1:c.4_6delinsATA").expect("parses");
        assert_eq!(codon_delins_replacement(&whole), Some(()));
        // Two codons wide: three bases replaced is not the test, one codon is.
        let wide = parse_hgvs("NM_000000.1:c.5_7delinsATA").expect("parses");
        assert_eq!(codon_delins_replacement(&wide), None);
        // A payload that is not three bases does not restore the triplet the
        // genomic split came from.
        let unbalanced = parse_hgvs("NM_000000.1:c.4_6delinsAT").expect("parses");
        assert_eq!(codon_delins_replacement(&unbalanced), None);
    }

    #[test]
    fn a_genomic_description_is_told_apart_from_a_transcript_one() {
        // The genomic half applies only to an axis the projector derived; this
        // is what keeps a genomic *input*'s own normalization at #79's answer.
        for genomic in ["NC_000001.11:g.4C>A", "NC_000001.11:g.[4C>A;6G>A]"] {
            assert!(is_genomic_description(
                &parse_hgvs(genomic).expect("parses")
            ));
        }
        for transcript in ["NM_000000.1:c.4C>A", "NM_000000.1:c.[4C>A;6G>A]"] {
            assert!(!is_genomic_description(
                &parse_hgvs(transcript).expect("parses")
            ));
        }
    }
}
