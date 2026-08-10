//! Issue #1264 — ferro must not build descriptions it cannot re-parse.
//!
//! The issue is filed as a *parser leniency* bug: "the parser accepted
//! non-adjacent `ins` anchors on the way in and rejects them on the way back
//! out". That diagnosis is wrong, and measuring it is what found the real
//! defects. `validate_no_point_insertion`
//! (`src/hgvs/parser/variant.rs`, #446/#450) already rejects both the
//! single-position and the non-adjacent anchor inbound on all six nucleotide
//! axes — `LRG_199:g.699646_700297insNNN…`, the exact string the issue quotes,
//! does not parse.
//!
//! The unparseable descriptions in the issue were never produced by
//! `parse_hgvs`. They were built by the **projector**, and they reached the
//! `FERRO_ASSERT_REPARSE` oracle only as its *input*. The oracle's exemption —
//! "skip when the input does not itself re-parse" — is a blanket wide enough to
//! swallow every one of them, which is precisely why #1264 exists ("an
//! exemption with no tracked follow-up quietly becomes permanent").
//!
//! Running the suite with that exemption instrumented to report rather than
//! return silently gives the census (9005 tests, 0 failures, 18 hits):
//!
//! | rendered shape | hits | cause |
//! |---|---|---|
//! | `LRG_199:g.699646_700297insNNN…` | 8 | insertion at a splice junction |
//! | `NC_000023.11:g.spl`, `g.(spl)` | 8 | the RNA-only `spl` edit on a `g.` axis |
//! | `[]` | 2 | an `AlleleVariant` built with no members |
//!
//! Each gets a test below, plus one for a genuine parser hole the issue did not
//! find: a `GenomeRing` segment is never anchor-checked at all.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::hgvs::variant::{AllelePhase, AlleleVariant};
use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, HgvsVariant, VariantProjector};

// ---------------------------------------------------------------------------
// The parser hole the issue's framing points at — real, but in `GenomeRing`.
// ---------------------------------------------------------------------------

/// A ring segment must obey the same two insertion-anchor rules as any other
/// insertion: `DNA/insertion.md:15` adjacency and `DNA/insertion.md:95-101`
/// single-position, both MUST-level.
///
/// A ring's `segments` are `LocEdit<GenomeInterval, NaEdit>` hanging off
/// `GenomeRing` rather than nested variants, so no validator reaches them by
/// matching on variant kind — they are validated only because the parse
/// driver's `for_each_leaf` walk re-wraps each segment as an ordinary `g.`
/// leaf (and unwraps the `sup` marker) before running the leaf rules
/// (#1264, #1578). This test pins that delivery for the insertion-anchor
/// rules, on both ring spellings.
///
/// A parser test is the only guard that can: an accepted escape round-trips
/// — ferro parses it *and* renders it back — so the re-parse oracle never
/// sees it.
#[test]
fn a_ring_segment_insertion_anchor_must_be_flanking() {
    for input in [
        // Non-adjacent anchors: 100 and 102 are two apart.
        "NC_000022.11:g.100_102insATG::200_201insG",
        // Single-position anchor.
        "NC_000022.11:g.100insATG::200_201insG",
        // Same, reached through the supernumerary wrapper.
        "NC_000022.11:g.[100_102insATG::200_201insG]sup",
    ] {
        assert!(
            parse_hgvs(input).is_err(),
            "`{input}` has an insertion anchor that is not a flanking pair; the same \
             spelling outside a `::`-join is rejected, so the ring must reject it too",
        );
    }
}

/// The control for the test above: a ring whose insertion anchors *are* flanking
/// is not rejected **by the anchor rule**.
///
/// Without this, tightening the ring path could pass by rejecting every ring
/// insertion for the wrong reason.
///
/// **The observable changed when `validate_ring_segments_are_wellformed` landed,
/// and the guarantee did not.** This used to assert the input parses. It no longer
/// does: an insertion segment designates no break point junction, so a ring may
/// not contain one at all (`DNA/complex.md:39`, and `:60-64` on the committee
/// withdrawing `::` as a general join operator). Every ring insertion is now
/// rejected — deliberately — which makes "it parses" unavailable as evidence.
///
/// So the control asserts the *reason* instead: this input is rejected by the
/// **ring** rule and not by the anchor rule, which is exactly the original claim
/// — a flanking pair is not something the anchor rule refuses. Asserting the
/// reason is in fact stronger than asserting acceptance was, because acceptance
/// could have come from the anchor rule silently not running.
///
/// Note the sibling test above still fails for its own reason: the per-leaf
/// validators run *before* the ring rules, so a non-flanking anchor is caught by
/// the anchor rule and never reaches the junction check.
#[test]
fn a_ring_segment_with_flanking_anchors_is_not_rejected_by_the_anchor_rule() {
    let input = "NC_000022.11:g.100_101insATG::200_201insG";
    let rendered = parse_hgvs(input)
        .expect_err("an insertion segment designates no break junction")
        .to_string();
    assert!(
        rendered.contains("break point junction"),
        "`{input}`'s anchors are a flanking pair, so it must be the ring rule that \
         refuses it, not the insertion-anchor rule; got: {rendered}"
    );
    assert!(
        !rendered.contains("single-position insertion") && !rendered.contains("flanking"),
        "the anchor rule must not be what rejects a flanking pair; got: {rendered}"
    );
}

// ---------------------------------------------------------------------------
// The projector defects that actually produced the issue's strings.
// ---------------------------------------------------------------------------

/// Two exons of 30 transcript bases each, separated by a 170-base intron.
///
/// The intron is the whole point: transcript positions 30 and 31 are adjacent
/// in the spliced RNA but 171 bases apart in the genome, so an insertion
/// anchored on that pair sits exactly on the exon–exon junction.
///
/// The sequence is deliberately non-repetitive. An earlier draft used `CGCG…`
/// exons and an all-`N` intron; against that reference the junction insertion
/// normalized into `g.1030_1199N[173]` — a repeat spanning the whole intron,
/// which *does* re-parse, so the test passed while proving nothing. A
/// homopolymer is not a representative intron.
const TX_SEQUENCE: &str = concat!(
    "ATGCTAGCATCGATCGTACGATCAGCTAGT", //  1..30  exon 1  (genomic 1000..1029)
    "GGACCTTGCAATCGGATCCAAGCTTGATAA", // 31..60  exon 2  (genomic 1200..1229)
);

/// A 170-base intron with no long repeat and no homopolymer run, generated
/// deterministically so the fixture stays reproducible.
///
/// The period-11 stride over a 4-letter alphabet keeps successive bases
/// differing and avoids the tandem structure that would let a junction
/// insertion be re-spelled as a repeat.
fn intron_sequence(len: usize) -> String {
    (0..len)
        .map(|i| ["A", "C", "G", "T"][(i * 7 + i / 11) % 4])
        .collect()
}

/// Exon 1 occupies genomic `[1000, 1030)`, exon 2 `[1200, 1230)` — so the
/// intron is 170 bases and transcript positions 30/31 map 171 apart.
const EXON1_GENOMIC_START: u64 = 1000;
const EXON1_GENOMIC_END: u64 = 1030;
const EXON2_GENOMIC_START: u64 = 1200;
const EXON2_GENOMIC_END: u64 = 1230;
const CONTIG: &str = "NC_000001.11";

/// The cdot alignment is built explicitly rather than via
/// `CdotMapper::from_transcripts`, because the genomic axis is only reachable
/// when a real exon alignment has been ingested. Without one the projector
/// reports `genomic = None` for every input and the tests below would pass
/// **vacuously** — declining to answer looks identical to declining correctly.
fn projector() -> VariantProjector<JsonProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("MYGENE".to_string()),
            contig: CONTIG.to_string(),
            strand: Strand::Plus,
            // [genomic_start, genomic_end, tx_start, tx_end)
            exons: vec![
                [EXON1_GENOMIC_START, EXON1_GENOMIC_END, 0, 30],
                [EXON2_GENOMIC_START, EXON2_GENOMIC_END, 30, 60],
            ],
            cds_start: Some(0),
            cds_end: Some(60),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );

    let mut provider = JsonProvider::new();
    provider.add_transcript(
        Transcript::new(
            "NM_TEST.1".to_string(),
            Some("MYGENE".to_string()),
            TxStrand::Plus,
            TX_SEQUENCE.to_string(),
            Some(1),
            Some(60),
            vec![Exon::new(1, 1, 30), Exon::new(2, 31, 60)],
            Some(CONTIG.to_string()),
            Some(EXON1_GENOMIC_START),
            Some(EXON2_GENOMIC_END - 1),
            Default::default(),
            ManeStatus::default(),
            None,
            None,
        )
        .with_protein_id(Some("NP_TEST.1".to_string())),
    );
    // Exon 1 at 1000, a 170-base intron, exon 2 at 1200 (all 1-based here).
    let prefix = intron_sequence((EXON1_GENOMIC_START - 1) as usize);
    let intron = intron_sequence((EXON2_GENOMIC_START - EXON1_GENOMIC_END) as usize);
    let suffix = intron_sequence(100);
    provider.add_genomic_sequence(
        CONTIG,
        format!(
            "{prefix}{}{intron}{}{suffix}",
            &TX_SEQUENCE[..30],
            &TX_SEQUENCE[30..],
        ),
    );

    VariantProjector::new(Projector::new(cdot), provider)
}

/// Every rendered axis of a projection of `input`.
fn projected_axes(input: &str) -> Vec<(&'static str, String)> {
    let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    let projection = projector()
        .project_variant(&variant, "NM_TEST.1")
        .unwrap_or_else(|e| panic!("project {input}: {e}"));
    [
        ("g", projection.genomic.as_ref().map(ToString::to_string)),
        ("c", projection.coding.as_ref().map(ToString::to_string)),
        ("n", projection.noncoding.as_ref().map(ToString::to_string)),
        ("r", projection.rna.as_ref().map(ToString::to_string)),
        ("p", projection.protein.as_ref().map(ToString::to_string)),
    ]
    .into_iter()
    .filter_map(|(axis, rendered)| rendered.map(|r| (axis, r)))
    .collect()
}

/// **The invariant.** Every axis a projection reports must be a description
/// ferro can read back.
///
/// Stated as a property over the projection rather than as an expected string,
/// so it survives a change to *which* axes a given input can reach. A projector
/// that declines an axis satisfies it vacuously — declining is always a legal
/// answer; emitting an unreadable description never is.
fn assert_every_projected_axis_re_parses(input: &str) {
    for (axis, rendered) in projected_axes(input) {
        parse_hgvs(&rendered).unwrap_or_else(|e| {
            panic!(
                "projecting `{input}` produced `{rendered}` on the {axis}. axis, which \
                 ferro cannot re-parse: {e}",
            )
        });
    }
}

/// `spl` is an RNA-only effect. Carried onto another axis it renders `g.spl` /
/// `c.1spl`, neither of which is HGVS.
///
/// `RNA/splicing.md` defines `r.spl` as "splicing affected" — a statement about
/// the transcript, with no genomic or coding counterpart. `hgvs_to_vcf` already
/// declines it for that reason (`vcf/from_hgvs.rs:790`), but the projector has
/// no `NaEdit::Splice` handling anywhere, so it converts the position and keeps
/// the edit token verbatim.
///
/// Both spellings are covered: `r.spl` (splicing very likely affected) and
/// `r.spl?` (might be affected), plus the predicted wrapper `r.(spl)`.
/// Both accession forms are covered because they take different code paths: a
/// bare `NM_` routes through `project_noncoding_direct` (which never reports a
/// genomic axis), while the gene-selector form
/// `NC_000001.11(NM_TEST.1)` goes through `project_to_genomic_nc`. The
/// enumeration hit `g.spl` only on the second and `c.1spl` on both.
#[test]
fn the_rna_only_splice_edit_has_no_representation_on_another_axis() {
    for input in [
        "NM_TEST.1:r.spl",
        "NM_TEST.1:r.spl?",
        "NM_TEST.1:r.(spl)",
        "NC_000001.11(NM_TEST.1):r.spl",
        "NC_000001.11(NM_TEST.1):r.spl?",
    ] {
        assert_every_projected_axis_re_parses(input);
    }
}

/// The same rule on the `n.` axis, on the route the direct path does not cover.
///
/// `0` ("no product") is not confined to `r.`: it is a valid spec form for a
/// non-coding transcript, and `n.0` parses to `HgvsVariant::Tx` because the `n.`
/// whole-entity path probes the same `parse_rna_no_product`. There are two
/// routes to a projection, and only one of them saw that.
/// `a_no_product_transcript_variant_reports_the_mismatch_before_the_rna_reason`
/// covers the direct route (bare `NM_` → `project_noncoding_direct`), whose
/// guard reads the extracted edit and so already declined. The genome-pivot
/// route (`NC_000001.11(NM_TEST.1)` → `project_to_genomic_nc`) guarded on
/// `if let HgvsVariant::Rna(v)` instead — a test on the variant *kind*, which a
/// `Tx` input passes straight through. It reported the genomic axis as
/// `NC_000001.11:g.0`: the same defect as `g.spl`, reached by the other axis.
///
/// The enumeration corpus carries no `n.0` input, so the census that found the
/// `r.` half could not see this one — which is why the fix is to ask the
/// question of the extracted edit on both routes rather than to add a second
/// variant-kind arm.
#[test]
fn the_rna_only_no_product_edit_is_declined_on_the_transcript_axis_too() {
    for input in [
        "NM_TEST.1:n.0",
        "NM_TEST.1:n.(0)",
        "NC_000001.11(NM_TEST.1):n.0",
        "NC_000001.11(NM_TEST.1):n.(0)",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        // Erroring the whole projection and declining per axis are both legal
        // answers, and the two accession forms currently give different ones —
        // a bare `NM_` routes through the direct path, the gene-selector form
        // through `project_to_genomic_nc`. What is never legal is reporting the
        // statement on an axis that cannot carry it. `n.` is excluded from the
        // check because it is the input's *own* axis: were a future change to
        // claim it the way `project_rna_only_effect` claims `r.`, that would be
        // this rule being followed, not broken.
        if let Ok(projection) = projector().project_variant(&variant, "NM_TEST.1") {
            let carried: Vec<String> = [
                ("g", projection.genomic.as_ref()),
                ("c", projection.coding.as_ref()),
                ("p", projection.protein.as_ref()),
            ]
            .into_iter()
            .filter_map(|(axis, v)| v.map(|v| format!("{axis}. = `{v}`")))
            .collect();
            assert!(
                carried.is_empty(),
                "`{input}` states that no product is made, which has no g./c./p. \
                 counterpart, yet the projection reported {}",
                carried.join(", "),
            );
        }
    }
}

/// An insertion whose transcript anchors straddle a splice junction has no
/// genomic spelling, so the genomic axis must be declined rather than invented.
///
/// `r.30_31insNNN` is well formed: 30 and 31 are adjacent in the spliced RNA.
/// Their genomic images are 1029 and 1200 — 171 bases apart, because the intron
/// lies between them. Rendering that pair as `g.1029_1200ins…` asserts an
/// insertion between two bases that are not neighbours, which is the very thing
/// `DNA/insertion.md:15` forbids and ferro's own parser rejects.
///
/// There is no repair available. The insertion could be attached to the 3' end
/// of exon 1 or the 5' end of exon 2, and those are different genomic variants;
/// HGVS has no way to say "at the junction". Declining is the honest answer,
/// and it is the same answer `spl` gets above.
#[test]
fn an_insertion_across_a_splice_junction_has_no_genomic_spelling() {
    // The gene-selector form, so the genomic axis is actually attempted — a
    // bare `NM_` routes through `project_noncoding_direct` and reports no
    // genomic axis at all, which would make this vacuous.
    assert_every_projected_axis_re_parses("NC_000001.11(NM_TEST.1):r.30_31insnnn");
}

/// Declining the genomic axis must not cost the axes that are still correct.
///
/// The junction only breaks the *genomic* spelling: `c.30_31insNNN` and
/// `n.30_31insNNN` are perfectly good descriptions of the same variant. An
/// earlier attempt guarded the genome pivot itself, which failed the whole
/// projection and threw those away — a worse answer than the bug, since the
/// caller went from "four right axes and one wrong one" to nothing at all.
/// The guard therefore sits on the reported axis, not on the pivot that also
/// feeds the c./n./p. derivation.
#[test]
fn a_declined_genomic_axis_does_not_take_the_other_axes_with_it() {
    let axes = projected_axes("NC_000001.11(NM_TEST.1):r.30_31insnnn");
    let reported: Vec<&str> = axes.iter().map(|(axis, _)| *axis).collect();
    assert!(
        !reported.contains(&"g"),
        "the genomic axis straddles the junction and must be withheld; got {axes:?}",
    );
    for axis in ["c", "n", "r"] {
        assert!(
            reported.contains(&axis),
            "the {axis}. axis is unaffected by the junction and must still be \
             reported; got {axes:?}",
        );
    }
}

/// The control, and the anti-vacuity guard for every projector test above.
///
/// An insertion *inside* an exon has genomically adjacent anchors, so it must
/// still reach the genomic axis. Without this, the fixture could report
/// `genomic = None` for everything — no alignment ingested, say — and the tests
/// above would pass while proving nothing, because "declined" and "never
/// attempted" are indistinguishable from the outside.
#[test]
fn an_insertion_inside_an_exon_still_projects_to_the_genomic_axis() {
    let axes = projected_axes("NC_000001.11(NM_TEST.1):r.10_11insnnn");
    let genomic = axes
        .iter()
        .find(|(axis, _)| *axis == "g")
        .map(|(_, rendered)| rendered.as_str());
    assert!(
        genomic.is_some(),
        "an intra-exonic insertion has adjacent genomic anchors and must still project — \
         if this is None the fixture reaches no genomic axis at all and every other \
         projector test here is vacuous; got {axes:?}",
    );
    let genomic = genomic.expect("checked above");
    parse_hgvs(genomic).unwrap_or_else(|e| panic!("`{genomic}` must re-parse: {e}"));
    assert!(
        genomic.contains("ins"),
        "the projected genomic form should still be an insertion, got `{genomic}`",
    );
}

/// Declining an axis must not cost the *reason*, even when the decline is
/// inherited from an allele member.
///
/// The allele's `.genomic` is all-or-nothing: one member straddling the junction
/// drops the whole genomic allele. The explanation for that drop already exists
/// — it is sitting in the member's own `axis_decline_reasons.genomic` — so
/// reporting the aggregate as an unexplained `None` throws away the only answer
/// anyone has. `.noncoding` has propagated its member's reason since #1198; this
/// pins the genomic axis to the same contract, which matters most for exactly
/// the mixed allele this issue introduces (one junction-straddling member, one
/// ordinary one).
#[test]
fn an_allele_inherits_the_genomic_decline_reason_of_the_member_that_caused_it() {
    let variant = parse_hgvs("NC_000001.11(NM_TEST.1):r.[30_31insnnn;40a>g]")
        .expect("the mixed allele must parse");
    let projection = projector()
        .project_variant(&variant, "NM_TEST.1")
        .expect("the allele still projects; only its genomic axis declines");

    assert!(
        projection.genomic.is_none(),
        "one member straddles the junction, so the all-or-nothing genomic allele \
         must be withheld; got {:?}",
        projection.genomic.as_ref().map(ToString::to_string),
    );
    let reason = projection
        .axis_decline_reasons
        .genomic
        .as_deref()
        .expect("the withheld genomic axis must say why, as `.noncoding` does");
    assert!(
        reason.contains("allele member"),
        "the reason must attribute the decline to a member rather than the whole \
         allele, since it describes one member's coordinates; got `{reason}`",
    );
    assert!(
        reason.contains("non-adjacent genomic positions"),
        "the member's own explanation must survive the trip up, not be replaced by \
         generic wording; got `{reason}`",
    );
}

// ---------------------------------------------------------------------------
// A declined axis must not swallow the transcript_id contract.
// ---------------------------------------------------------------------------

/// An RNA-only effect answered up front must still respect the transcript it is
/// being projected against.
///
/// `project_rna_only_effect` short-circuits ahead of the normal dispatch, which
/// is where every other transcript-coordinate input gets its `transcript_id`
/// checked. Without its own check the short-circuit hands back a projection
/// stamped with a transcript the variant is not on — and because the fan-out
/// path calls the projector once per overlapping transcript, a single `r.spl`
/// would be reported against every one of them, each row claiming a different
/// `transcript_id` while carrying the same `r.` value.
///
/// Both accession shapes are covered because they reach the check by different
/// routes: a bare `NM_` falls through to `project_noncoding_direct`, the
/// gene-selector form to the genome-pivot path.
#[test]
fn an_rna_only_effect_is_still_refused_against_a_foreign_transcript() {
    for input in [
        "NM_TEST.1:r.spl",
        "NM_TEST.1:r.spl?",
        "NM_TEST.1:r.0",
        "NC_000001.11(NM_TEST.1):r.spl",
    ] {
        let variant = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
        let err = projector()
            .project_variant(&variant, "NM_FOREIGN.9")
            .expect_err("an RNA-only effect on another transcript must be refused");
        let rendered = err.to_string();
        assert!(
            rendered.contains("transcript_id mismatch"),
            "projecting `{input}` against a foreign transcript must report the \
             mismatch, not the RNA-level reason; got `{rendered}`",
        );
    }
}

/// The transcript check must also come first for `n.0`, which the up-front
/// short-circuit does not own.
///
/// `project_rna_only_effect` only claims `r.` inputs, so a `Tx` no-product
/// (`n.0` — a valid spec form) still reaches `project_noncoding_direct`. There
/// the RNA-only guard and the mismatch guard are both applicable, and the
/// mismatch is the more fundamental objection: it holds whatever the edit is, so
/// answering "this edit has no other representation" for a variant that was
/// never on this transcript describes the wrong problem.
#[test]
fn a_no_product_transcript_variant_reports_the_mismatch_before_the_rna_reason() {
    let variant = parse_hgvs("NM_TEST.1:n.0").expect("`n.0` is a valid spec form");

    let foreign = projector()
        .project_variant(&variant, "NM_FOREIGN.9")
        .expect_err("`n.0` on another transcript must be refused");
    assert!(
        foreign.to_string().contains("transcript_id mismatch"),
        "the mismatch outranks the RNA-level reason; got `{foreign}`",
    );

    // The control: on its *own* transcript the RNA-level reason is still what
    // `n.0` gets, so the reordering above narrowed nothing.
    let own = projector()
        .project_variant(&variant, "NM_TEST.1")
        .expect_err("`n.0` has no representation on another axis");
    assert!(
        own.to_string().contains("RNA-level effect"),
        "on its own transcript `n.0` must still be declined as an RNA-level \
         effect; got `{own}`",
    );
}

// ---------------------------------------------------------------------------
// The empty allele.
// ---------------------------------------------------------------------------

/// The hazard `AlleleVariant::checked` exists to prevent, characterized on the
/// **unchecked** constructor so the test is not merely restating `checked`.
///
/// `AlleleVariant::new(vec![], _)` renders `[]`, and `[]` is not HGVS — so the
/// value does not round-trip and every consumer that chains a second call fails
/// on it. That asymmetry is the empty-allele half of #1264.
///
/// The parser refuses empty brackets in every position that accepts them
/// (`empty_brackets_do_not_parse` below), so the shape is reachable only by
/// direct construction.
#[test]
fn an_allele_with_no_members_renders_something_that_does_not_parse() {
    let empty = HgvsVariant::Allele(AlleleVariant::new(Vec::new(), AllelePhase::Cis));
    assert_eq!(empty.to_string(), "[]");
    assert!(
        parse_hgvs(&empty.to_string()).is_err(),
        "`[]` is not HGVS, so an empty allele does not round-trip — which is why a member \
         list assembled from data must go through `AlleleVariant::checked`",
    );
}

/// The checked constructor refuses exactly that shape, and nothing else.
///
/// `new` is deliberately left infallible — it is called throughout the parser
/// and normalizer on lists that are non-empty by construction. `checked` is for
/// the boundaries where a list is assembled from data (a projection, a filter,
/// an external payload), mirroring `GenomeRing::new`, which refuses a ring of
/// fewer than two segments for the same round-trip reason.
#[test]
fn the_checked_allele_constructor_refuses_an_empty_member_list() {
    assert!(
        AlleleVariant::checked(Vec::new(), AllelePhase::Cis).is_none(),
        "an allele with no members has no HGVS rendering, so it must not be constructible",
    );
    let single = AlleleVariant::checked(
        vec![parse_hgvs("NC_000022.11:g.100A>G").expect("parses")],
        AllelePhase::Cis,
    );
    assert!(single.is_some(), "a one-member allele is still valid");
}

/// The parser's half of the same rule, pinned so it cannot regress.
#[test]
fn empty_brackets_do_not_parse() {
    for input in ["[]", "NC_000022.11:g.[]", "NC_000022.11:g.[];[]"] {
        assert!(parse_hgvs(input).is_err(), "`{input}` must not parse");
    }
}

/// Guard against re-introducing the shape by hand: rendering any allele ferro
/// builds must round-trip.
#[test]
fn a_constructed_allele_round_trips() {
    let members: Vec<HgvsVariant> = ["NC_000022.11:g.100A>G", "NC_000022.11:g.200del"]
        .iter()
        .map(|s| parse_hgvs(s).expect("parses"))
        .collect();
    let allele = AlleleVariant::checked(members, AllelePhase::Cis).expect("two members is valid");
    let rendered = HgvsVariant::Allele(allele).to_string();
    parse_hgvs(&rendered).unwrap_or_else(|e| panic!("`{rendered}` must re-parse: {e}"));
}

/// `checked` at the boundary its own doc comment names, exercised through that
/// boundary rather than directly.
///
/// `AlleleVariant::checked` is only worth adding if something calls it, and
/// `project_to_genomic`'s allele branch is exactly the case it describes — a
/// member list assembled from data rather than fixed by the grammar. Asserting
/// on the constructor alone (`the_checked_allele_constructor_refuses_an_empty_member_list`)
/// leaves the wiring untested, so this drives the same rule through the
/// projector and pins that it declines with an error instead of handing back an
/// allele that renders `[]`.
#[test]
fn projecting_an_empty_allele_to_the_genomic_axis_is_declined() {
    let empty = HgvsVariant::Allele(AlleleVariant::new(Vec::new(), AllelePhase::Cis));
    let err = projector()
        .project_to_genomic(&empty)
        .expect_err("an empty allele has no genomic rendering, so projection must decline");
    let message = err.to_string();
    assert!(
        message.contains("no members"),
        "the decline must say what was wrong with the input; got `{message}`",
    );
}
