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
///   already-corrected description does not re-emit them. **This carry is
///   defensive today and moves nothing** — see below.
///
/// # The warnings carry is defensive, not live coverage
///
/// On every production path both sides are empty at this call site, so the
/// third line of the body is a no-op. The warnings a projection reports are
/// attached by the *outer* entry points — `project_bare_genomic_parent`,
/// `project_variant`, and both arms of `project_variant_all` — and every one of
/// them assigns after `project_variant_inner` (or `project_normalized_all_inner`)
/// has already returned. `apply_codon_frame_exception` runs *inside*
/// `project_variant_inner`, so nothing has written warnings onto either
/// projection by the time it is reached, and the outer assignment overwrites
/// whatever it did. `project_normalized` assigns none at all.
///
/// It is kept rather than deleted because the invariant it encodes — the
/// warnings describe the **caller's input**, not the merged form the
/// re-projection re-normalized — is the one a future pre-merge warning source
/// would need, and the line is what makes this function's contract complete
/// rather than accidentally correct. Its unit coverage is stated as defensive
/// for the same reason: read that assertion as pinning the contract, not as
/// evidence any shipped path exercises it.
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

/// **Every** adjacent member pair on `coding` that `delins.md:18` keeps
/// together, each as its two CDS positions, in member order.
///
/// # This is the TRIGGER, not the authorization
///
/// A non-empty answer means the coding axis holds at least one partition the
/// frame would not have chosen, so the re-normalization is worth attempting.
/// It does **not** enumerate everything that re-normalization will merge, and
/// nothing may treat it as if it did — see the section below, and
/// [`merged_member_spans`], which is what actually authorizes the result.
///
/// The pairs are returned rather than a `bool` so the caller can run its exon
/// test only when there is something to test: a transcript lookup on every
/// projection is a cost this rule has no right to impose on descriptions it
/// never touches.
///
/// # Why every pair, and not just the first
///
/// The merge this feeds is `merge::apply_coding_codon_exception`, and that
/// function loops (`while index < pieces.len()`) — it merges *each* qualifying
/// pair it walks past, with no exon awareness of its own, `same_codon` being
/// pure CDS arithmetic. So a predicate that answered for the first pair alone
/// would guard one pair while licensing the rest. Pinned by
/// `a_second_candidate_pair_crossing_an_exon_junction_declines_the_whole_merge`.
///
/// # The shape, and where it is narrower than what the re-normalization merges
///
/// A pair here is two **lone single-base** `c.` substitutions with exactly one
/// unchanged base between them, codon-tested on the CDS axis where `same_codon`
/// is meaningful. What the re-normalization actually merges is wider in two
/// separate ways, and **both were measured**, not inferred:
///
/// * `merge::apply_coding_codon_exception` is **asymmetric**: it requires the
///   *left* piece to be a lone substitution and the *right* only to **begin**
///   with one — its `starts_with_substitution` test is
///   `ref_end - ref_start == alt.len()`, i.e. any same-length replacement block.
///   Measured on `NM_SPLIT2.1`: `c.[4G>C;6_7delinsAG]` has **no** pair here and
///   normalizes to `c.4_7delinsCAAG`.
/// * The re-normalization is a whole `normalize` call, so **every** merge rule
///   runs, not only this exception. `delins-adjacent-members-when-both-consume-reference`
///   merges at separation **zero**, which this predicate does not look at at
///   all.
///
/// **So a whole-allele `all()` over these pairs is not a safe authorization.**
/// It reads as one — "every candidate cleared the exon test" — while saying
/// nothing about the merges it never enumerated. Measured end-to-end through
/// `project_variant` on `NM_SPLIT2.1`, same geometry, opposite answers:
///
/// ```text
/// g.[1009G>C;1100_1103delinsCAAC]                  -> c.[10G>C;11_14delinsCAAC]
/// g.[1003G>C;1005T>A;1009G>C;1100_1103delinsCAAC]  -> c.[4_6delinsCAA;10_14delinsCCAAC]
/// ```
///
/// `c.10` ends exon 1 and `c.11` opens exon 2, so the second line merged across
/// the junction — licensed only by an *unrelated* pair, `(4, 6)`, happening to
/// be enumerable and to sit inside exon 1. The paragraph this replaces argued
/// the narrowness was safe because "being narrower can only decline a merge the
/// normalizer would have made". That holds for a **per-pair** gate; it is false
/// for a whole-allele one, which is what this feeds. Pinned by
/// `an_unenumerated_merge_crossing_an_exon_junction_declines_the_whole_merge`.
///
/// **The caller still has to check that each pair lies inside one exon.** A
/// codon may be split across an exon/exon junction, and then the pair is
/// "separated by one nucleotide" on the `c.` axis while being separated by a
/// whole intron on the sequence. Merging those is a different question from the
/// one `delins.md:18` settles, and the repository already answers it the other
/// way for the same shape — `project_cis_intron_split_codon_combines_to_single_missense`
/// pins `c.[4G>C;6T>A]` staying split across a ~100 bp intron and combining only
/// at the protein level — so that decision is left where it already is.
///
/// A non-empty answer says only that the coding axis has a partition the frame
/// would not have chosen. The projector re-normalizes the coding allele to get
/// the merged form, so the merge itself stays in one place.
pub(crate) fn coding_codon_pairs(coding: &HgvsVariant) -> Vec<(i64, i64)> {
    let HgvsVariant::Allele(allele) = coding else {
        return Vec::new();
    };
    if allele.uncertain || allele.phase != AllelePhase::Cis {
        return Vec::new();
    }
    allele
        .variants
        .windows(2)
        .filter_map(
            |pair| match (cds_substitution(&pair[0]), cds_substitution(&pair[1])) {
                (Some((left, _)), Some((right, _)))
                    if right == left + 2 && same_codon(left, right) =>
                {
                    Some((left, right))
                }
                _ => None,
            },
        )
        .collect()
}

/// The CDS spans of the members the re-normalization **produced** — every
/// member of `merged` that was not already a member of `coding`, as
/// `(first, last)` CDS bases.
///
/// # Why the authorization is here and not on [`coding_codon_pairs`]
///
/// The caller has to answer "may this merge stand?", and the honest way to
/// answer it is to look at what the merge did, not to predict it. Predicting
/// requires re-deriving the normalizer's whole rule set — its piece
/// decomposition, `apply_coding_codon_exception`'s asymmetric right-hand shape,
/// the separation-zero adjacency merge, and whatever is added next — in a
/// second place, and a prediction that falls behind by one rule silently
/// becomes a licence. That is exactly how the junction-crossing merge recorded
/// on [`coding_codon_pairs`] got through.
///
/// Reading the result costs one pass over the members and cannot fall behind,
/// because it is stated over the output rather than over the rules.
///
/// # `None` means "cannot verify", and the caller must decline on it
///
/// A produced member whose span is not plainly on the positive CDS body — an
/// intronic offset, a `*`-numbered or negative position, an uncertain endpoint,
/// a non-`c.` member — cannot be exon-tested by
/// `codon_pair_is_within_one_exon`, which is CDS arithmetic. There is then no
/// basis to claim it stays inside one exon, which is the same direction
/// `codon_pair_is_within_one_exon` already takes for an unservable transcript.
///
/// Members are compared by their rendered form. Two members that render
/// identically denote the same thing on the same axis, which is the property
/// this needs: anything the merge left untouched is not the merge's to justify.
pub(crate) fn merged_member_spans(
    coding: &HgvsVariant,
    merged: &HgvsVariant,
) -> Option<Vec<(i64, i64)>> {
    let before: std::collections::BTreeSet<String> =
        members(coding).map(ToString::to_string).collect();
    members(merged)
        .filter(|member| !before.contains(&member.to_string()))
        .map(cds_span)
        .collect()
}

/// The top-level members of a description: an allele's own, or the description
/// itself when it is a lone variant. The iterator companion of
/// [`member_count`].
fn members(variant: &HgvsVariant) -> impl Iterator<Item = &HgvsVariant> {
    match variant {
        HgvsVariant::Allele(allele) => allele.variants.iter(),
        single => std::slice::from_ref(single).iter(),
    }
}

/// `(first, last)` CDS bases of a `c.` member sitting plainly in the CDS.
///
/// `None` on the same grounds as [`cds_substitution`]: an intronic offset, a
/// `*`-numbered 3'UTR position, a `pter`-class marker or a 5'UTR (negative)
/// base is not a position the exon arithmetic can place. Unlike
/// `cds_substitution` this admits any edit and any span — it answers *where* a
/// member sits, not *what shape* it is.
fn cds_span(variant: &HgvsVariant) -> Option<(i64, i64)> {
    let HgvsVariant::Cds(cds) = variant else {
        return None;
    };
    let start = certain_boundary(&cds.loc_edit.location.start)?;
    let end = certain_boundary(&cds.loc_edit.location.end)?;
    for position in [start, end] {
        if position.offset.is_some()
            || position.utr3
            || position.special.is_some()
            || position.base < 1
        {
            return None;
        }
    }
    Some((start.base, end.base))
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

    fn pairs(description: &str) -> Vec<(i64, i64)> {
        coding_codon_pairs(&parse_hgvs(description).expect("parses"))
    }

    /// No pair matched, spelled once so the assertions below read as a claim
    /// about the predicate rather than as a bare empty literal.
    fn no_pairs() -> Vec<(i64, i64)> {
        Vec::new()
    }

    #[test]
    fn a_gap_of_one_inside_one_codon_is_the_shape_the_exception_keeps_together() {
        // `c.4` and `c.6` are both codon 2, separated by the single unchanged
        // `c.5` — `DNA/delins.md:18`'s "two variants separated by one
        // nucleotide, together affecting one amino acid".
        assert_eq!(pairs("NM_000000.1:c.[4C>A;6G>A]"), vec![(4, 6)]);
    }

    #[test]
    fn every_qualifying_pair_is_surfaced_not_only_the_first() {
        // `merge::apply_coding_codon_exception` merges each qualifying pair it
        // walks past, so the caller's exon test has to be able to see all of
        // them. `c.4`/`c.6` are codon 2 and `c.10`/`c.12` are codon 4; answering
        // with the first alone would leave the second merged unchecked.
        assert_eq!(
            pairs("NM_000000.1:c.[4C>A;6G>A;10C>A;12G>A]"),
            vec![(4, 6), (10, 12)]
        );
    }

    fn spans(coding: &str, merged: &str) -> Option<Vec<(i64, i64)>> {
        merged_member_spans(
            &parse_hgvs(coding).expect("parses"),
            &parse_hgvs(merged).expect("parses"),
        )
    }

    /// The authorization sees what the enumeration missed, which is the whole
    /// reason it is stated over the result.
    ///
    /// The pair `(10, 12)` is invisible to [`coding_codon_pairs`] here —
    /// `c.11_14delinsCAAC` is not a lone substitution — yet `c.10_14delinsCCAAC`
    /// is a member the merge produced, and its span is exactly what the exon
    /// test has to be handed. Measured on `NM_SPLIT2.1`, where `c.10` ends exon
    /// 1 and `c.11` opens exon 2, so this span is the junction-crossing merge.
    #[test]
    fn the_result_check_reports_a_span_the_pair_enumeration_never_saw() {
        let coding = "NM_SPLIT2.1:c.[4C>A;6G>A;10G>C;11_14delinsCAAC]";
        assert_eq!(pairs(coding), vec![(4, 6)]);
        assert_eq!(
            spans(coding, "NM_SPLIT2.1:c.[4_6delinsCAA;10_14delinsCCAAC]"),
            Some(vec![(4, 6), (10, 14)]),
            "both produced members must be reported, including the one no pair \
             enumerated"
        );
    }

    /// A member the merge left alone is not the merge's to justify.
    #[test]
    fn an_untouched_member_is_not_reported_as_produced() {
        assert_eq!(
            spans(
                "NM_SPLIT2.1:c.[4C>A;6G>A;20G>C]",
                "NM_SPLIT2.1:c.[4_6delinsAGA;20G>C]"
            ),
            Some(vec![(4, 6)])
        );
    }

    /// A produced member the CDS arithmetic cannot place refuses to answer, so
    /// the caller declines rather than assuming it stays inside an exon.
    #[test]
    fn a_produced_member_off_the_cds_body_cannot_be_verified() {
        for merged in [
            "NM_SPLIT2.1:c.[*4_*6delinsAGA]",
            "NM_SPLIT2.1:c.[-6_-4delinsAGA]",
            "NM_SPLIT2.1:c.[4+1_4+3delinsAGA]",
        ] {
            assert_eq!(
                spans("NM_SPLIT2.1:c.[4C>A;6G>A]", merged),
                None,
                "{merged} is not placeable on the CDS body, so it must not be \
                 reported as verified"
            );
        }
    }

    #[test]
    fn a_gap_of_two_is_left_split_because_general_md_34_asks_for_it() {
        // `general.md:34` — "two variants separated by one or more nucleotides
        // should be described individually" — governs everything the exception
        // does not carve out, and the carve-out is a gap of exactly one.
        assert_eq!(pairs("NM_000000.1:c.[4C>A;7G>A]"), no_pairs());
    }

    #[test]
    fn a_gap_of_one_that_straddles_two_codons_is_not_one_amino_acid() {
        // `c.5` and `c.7` are codon 2 and codon 3, so the pair affects two amino
        // acids and the exception's precondition fails even though the
        // separation is right.
        assert_eq!(pairs("NM_000000.1:c.[5C>A;7G>A]"), no_pairs());
    }

    #[test]
    fn an_intronic_or_utr_member_is_not_on_the_coding_body() {
        // `same_codon` is only meaningful on the positive CDS body, so an
        // intronic offset, a `*`-numbered 3'UTR position and a negative 5'UTR
        // position each decline rather than being counted into a codon.
        assert_eq!(pairs("NM_000000.1:c.[4+1C>A;6G>A]"), no_pairs());
        assert_eq!(pairs("NM_000000.1:c.[*4C>A;*6G>A]"), no_pairs());
        assert_eq!(pairs("NM_000000.1:c.[-6C>A;-4G>A]"), no_pairs());
    }

    #[test]
    fn a_trans_allele_is_not_two_variants_on_one_sequence() {
        // The exception merges two changes that lie on the same molecule; `;`
        // in a trans allele asserts they do not.
        assert_eq!(pairs("NM_000000.1:c.[4C>A];[6G>A]"), no_pairs());
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
        // DEFENSIVE, not production coverage: on every shipped path both sides
        // are empty here, because the warnings are attached by the outer entry
        // points after `project_variant_inner` returns. This constructs the
        // pre-merge warnings by hand to pin the contract — that they describe
        // the caller's input and so must not be replaced by the merged form's —
        // and cannot fail against any current caller. See the note on
        // `carry_pre_merge_state`.
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
        assert_eq!(coding_codon_pairs(&certain), vec![(4, 6)]);
        assert_eq!(
            coding_codon_pairs(&make_edits_uncertain(&certain)),
            no_pairs()
        );
    }
}
