//! End-to-end tests for the opt-in confluence self-check (#1892).
//!
//! The self-check answers the question a consumer actually asks — *will my
//! variants change, and are any of them non-confluent on my data?* — by
//! composing two decided pieces of machinery: the equivalence relation
//! (`EquivalenceChecker`) and the normalizer. It groups the inputs into
//! equivalence classes under a caller-chosen relation and flags any class whose
//! members reach more than one distinct normalized output. It reports; it never
//! decides which relation gates a release (that is the deferred adjudication in
//! #1890), and it emits no pass/fail verdict.
//!
//! # The hermetic violation witnesses
//!
//! This file used to say there could be no hermetic "real violation" test —
//! that the decided non-confluences all need a real reference window. That was
//! wrong, and it left the feature's headline behaviour (*report the
//! non-confluence*) with no end-to-end guard at all: replacing
//! `violations.push(group)` in `check_confluence` with a no-op left every test
//! in the repository green, because every one of them asserted
//! `violations == []`.
//!
//! Two witnesses close it, and both are `MockProvider`-only. Each is a pair that
//! the cross-axis relation groups into **one** class — apply-equal on every
//! determined axis, established by `compare_denotations`, which never
//! normalizes — whose members the normalizer carries to **two** distinct
//! strings. See `a_class_that_normalizes_apart_is_reported` and
//! `two_no_ops_at_different_positions_are_a_witness`. If either pair ever
//! converges, that is a representation change and the test must be re-pointed at
//! another divergent pair, not deleted: what it guards is that the report
//! *fires*, and a green run over a corpus with nothing to find does not guard
//! that.

use ferro_hgvs::equivalence::{
    ConfluenceRelation, ConfluenceSkipKind, DeclineSite, EquivalenceChecker, TripleDecline,
};
use ferro_hgvs::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MockProvider};

/// The 40-base contig `spdi::apply`'s own tests use.
fn genomic_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_KEY.1", "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT");
    provider
}

fn variants(descriptors: &[&str]) -> Vec<HgvsVariant> {
    descriptors
        .iter()
        .map(|d| parse_hgvs(d).expect("fixture must parse"))
        .collect()
}

/// A two-exon transcript whose junction splits a two-base `T` run, so
/// duplicating either base yields one transcript sequence and two different
/// genomic sequences — the `c.3921dup`/`c.3922dup` shape, reduced. This is the
/// only geometry where the cross-axis relation and the (own-axis) SPDI relation
/// group differently.
fn two_exon_provider() -> MockProvider {
    const TX: &str = "NM_TWOEXON.1";
    const EXON1_TX: &str = "ACGCACGCAT";
    const EXON2_TX: &str = "TATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCGTACGTA";
    let tx_sequence = format!("{EXON1_TX}{EXON2_TX}");
    let exon2_end_tx = 11 + EXON2_TX.len() as u64 - 1;
    let exon2_end_genomic = 2001 + EXON2_TX.len() as u64 - 1;
    let mut genomic = String::new();
    genomic.push_str(&"A".repeat(1000));
    genomic.push_str(EXON1_TX);
    genomic.push_str(&"G".repeat(990));
    genomic.push_str(EXON2_TX);
    genomic.push_str(&"A".repeat(91));
    let transcript = Transcript::new(
        TX.to_string(),
        Some("TWOEXON".to_string()),
        Strand::Plus,
        tx_sequence,
        Some(1),
        Some(78),
        vec![
            Exon::with_genomic(1, 1, 10, 1001, 1010),
            Exon::with_genomic(2, 11, exon2_end_tx, 2001, exon2_end_genomic),
        ],
        Some("chr_two_exon".to_string()),
        Some(1001),
        Some(exon2_end_genomic),
        GenomeBuild::GRCh38,
        ManeStatus::None,
        None,
        None,
    );
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("chr_two_exon", genomic);
    provider.add_transcript(transcript);
    provider
}

/// Several spellings of one variant that all normalize to one string form a
/// single equivalence class with a single output, so the corpus is confluent
/// and no violation is reported.
#[test]
fn a_confluent_corpus_reports_no_violations() {
    let checker = EquivalenceChecker::new(genomic_provider());
    // A literal insertion and the same bases named by reference range: cross-axis
    // apply-equal, and the normalizer converges both.
    let corpus = variants(&["NC_KEY.1:g.5_6insGGCATTAGCCT", "NC_KEY.1:g.5_6ins9_19"]);

    let report = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);

    assert!(report.is_confluent(), "one variant, one output: {report:?}");
    assert!(report.violations().is_empty());
    assert!(report.skipped().is_empty());
    assert_eq!(
        report.classes_checked(),
        1,
        "the two spellings are one equivalence class"
    );
    assert_eq!(report.undecided_pairs(), 0, "the pair was decided");
    assert!(
        report.is_complete(),
        "nothing was skipped and nothing was undecidable"
    );
    assert_eq!(
        report.relation(),
        ConfluenceRelation::CrossAxisSequenceMatch
    );
}

/// Two genuinely different variants land in separate classes, and two spellings
/// of one variant land together.
///
/// The corpus is deliberately three inputs reducing to **two** classes rather
/// than two inputs reducing to two. The earlier two-input form asserted
/// `classes_checked == 2`, which is also exactly what a grouping that never
/// merges anything produces — so it pinned the count without pinning the
/// grouping, and survived "make `group_by_cross_axis` never merge" untouched.
/// Three inputs and two classes is a number only a working grouping can reach.
#[test]
fn distinct_variants_form_separate_classes() {
    let checker = EquivalenceChecker::new(genomic_provider());
    // `g.3delinsG` is `g.3A>G` respelled: apply-equal, so it must join it.
    // `g.3A>C` denotes a different base, so it must not.
    let corpus = variants(&["NC_KEY.1:g.3A>G", "NC_KEY.1:g.3A>C", "NC_KEY.1:g.3delinsG"]);

    let report = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);

    assert!(report.is_confluent());
    assert!(report.violations().is_empty());
    assert_eq!(
        report.classes_checked(),
        2,
        "three inputs, two variants — a grouping that never merges would say 3: {report:?}"
    );
    assert!(report.is_complete(), "{report:?}");
}

/// An input the chosen relation cannot place — a non-literal `insN[5]` payload
/// has no SPDI key — is reported as skipped, with a reason, rather than silently
/// dropped or falsely grouped. A skip that read as a pass is the failure this
/// self-check exists to avoid.
#[test]
fn inputs_the_relation_cannot_place_are_reported_as_skipped() {
    let checker = EquivalenceChecker::new(genomic_provider());
    // An unspecified inserted length denotes no resolvable sequence, so it has
    // no SPDI key (`spdi_key`'s documented refusals).
    let corpus = variants(&["NC_KEY.1:g.10_11ins(10)", "NC_KEY.1:g.5_6ins(6)"]);

    for relation in [
        ConfluenceRelation::Spdi,
        ConfluenceRelation::CrossAxisSequenceMatch,
    ] {
        let report = checker.check_confluence(&corpus, relation);

        assert!(report.violations().is_empty());
        assert_eq!(
            report.skipped().len(),
            2,
            "neither non-literal payload denotes a sequence under {relation:?}: {report:?}"
        );
        assert!(
            report
                .skipped()
                .iter()
                .all(|s| s.kind == ConfluenceSkipKind::Unplaceable),
            "the relation never placed either: {report:?}"
        );
        assert!(
            report.skipped().iter().all(|s| !s.reason.is_empty()),
            "every skip must carry a reason"
        );
        assert_eq!(
            report.classes_checked(),
            0,
            "an unplaceable input forms no class, under either relation: {report:?}"
        );
        assert!(
            !report.is_complete(),
            "a corpus that was entirely skipped is not a complete run: {report:?}"
        );
    }
}

/// The chosen relation changes the grouping. Two junction-straddling
/// duplications denote the same transcript bases (one SPDI key) but different
/// genomic bases, so the SPDI relation groups them into one class while the
/// cross-axis relation keeps them in two. Both are confluent here; what is being
/// pinned is that the relation parameter is honored.
///
/// The corpus carries a **third** input, `c.11delinsTT`, which is `c.11dup`
/// respelled and which the cross-axis relation therefore must merge with it.
/// Without it the cross-axis assertion read `classes_checked == 2` over a
/// two-input corpus — the same number a grouping that never merges anything
/// produces — so only the SPDI half of this test was load-bearing and the
/// mutation "make `group_by_cross_axis` never merge" left it green. Three inputs
/// and two cross-axis classes cannot be reached without merging.
#[test]
fn the_relation_parameter_changes_the_grouping() {
    let checker = EquivalenceChecker::new(two_exon_provider());
    let corpus = variants(&[
        "NM_TWOEXON.1:c.10dup",
        "NM_TWOEXON.1:c.11dup",
        "NM_TWOEXON.1:c.11delinsTT",
    ]);

    let spdi = checker.check_confluence(&corpus, ConfluenceRelation::Spdi);
    assert_eq!(
        spdi.classes_checked(),
        1,
        "same transcript bases: one SPDI class over all three: {spdi:?}"
    );

    let cross = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);
    assert_eq!(
        cross.classes_checked(),
        2,
        "`c.11dup` and `c.11delinsTT` share a genomic denotation and must merge; \
         `c.10dup` must not join them — a grouping that never merges would say 3: {cross:?}"
    );

    assert!(spdi.is_confluent() && cross.is_confluent());
    assert!(spdi.is_complete() && cross.is_complete());
}

/// An empty corpus is vacuously confluent.
#[test]
fn an_empty_corpus_is_confluent() {
    let checker = EquivalenceChecker::new(genomic_provider());
    let report = checker.check_confluence(&[], ConfluenceRelation::CrossAxisSequenceMatch);
    assert!(report.is_confluent());
    assert_eq!(report.classes_checked(), 0);
    assert!(report.skipped().is_empty());
    assert_eq!(report.undecided_pairs(), 0);
    assert!(report.is_complete(), "nothing was handed over to miss");
}

/// **The headline behaviour, end to end.** A class the relation groups whose
/// members the normalizer carries to two distinct strings is reported as a
/// violation, carrying both inputs and both outputs.
///
/// The witness is one variant spelled two ways: `g.[3A>G;10=]` states an
/// unchanged base alongside the substitution, `g.3A>G` does not. Applying either
/// to the reference produces the same bases — an `=` member changes nothing — so
/// `compare_denotations` reaches `CrossAxisSequenceMatch` and the two are one
/// class. The normalizer preserves the `=` member, so the class reaches two
/// outputs.
///
/// This is a real non-confluence in ferro, not a contrivance: a caller that
/// records the bases it checked is spelling the same variant. If it is ever
/// converged, re-point this test at another divergent pair — see the module
/// docs.
#[test]
fn a_class_that_normalizes_apart_is_reported() {
    let checker = EquivalenceChecker::new(genomic_provider());
    let corpus = variants(&["NC_KEY.1:g.[3A>G;10=]", "NC_KEY.1:g.3A>G"]);

    let report = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);

    assert_eq!(
        report.classes_checked(),
        1,
        "apply-equal on every determined axis, so one class: {report:?}"
    );
    assert!(
        report.is_complete(),
        "nothing here is skipped or undecidable, so the finding is about the corpus: {report:?}"
    );
    assert!(
        !report.is_confluent(),
        "one class, two outputs, is a non-confluence witness: {report:?}"
    );

    let violations = report.violations();
    assert_eq!(violations.len(), 1, "{report:?}");
    assert_eq!(
        violations[0].inputs,
        vec![
            "NC_KEY.1:g.[3A>G;10=]".to_string(),
            "NC_KEY.1:g.3A>G".to_string()
        ],
        "the inputs are recorded in input order so a consumer sees which collided"
    );
    assert_eq!(
        violations[0].outputs,
        vec![
            "NC_KEY.1:g.3A>G".to_string(),
            "NC_KEY.1:g.[3A>G;10=]".to_string()
        ],
        "the distinct outputs the class reached, sorted"
    );
}

/// A second, independent violation witness, on the shape `spdi_key`'s own docs
/// single out: two descriptions that change nothing.
///
/// `g.10=` and `g.20=` both denote the unchanged reference, so they are
/// apply-equal on every determined axis and the cross-axis relation groups them;
/// the normalizer keeps each at its own position, so the class reaches two
/// outputs. This one is deliberately stable — `src/equivalence/key.rs` records
/// that giving two no-ops one answer "means choosing one, and nothing here needs
/// it", so converging them is a decision nobody has made rather than a defect
/// waiting to be fixed. Under `Spdi` they key at their own positions and do not
/// group, which is the same fact seen from the weaker relation.
#[test]
fn two_no_ops_at_different_positions_are_a_witness() {
    let checker = EquivalenceChecker::new(genomic_provider());
    let corpus = variants(&["NC_KEY.1:g.10=", "NC_KEY.1:g.20="]);

    let cross = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);
    assert_eq!(cross.classes_checked(), 1, "{cross:?}");
    assert!(cross.is_complete(), "{cross:?}");
    assert!(!cross.is_confluent(), "{cross:?}");
    assert_eq!(
        cross.violations()[0].outputs,
        vec!["NC_KEY.1:g.10=".to_string(), "NC_KEY.1:g.20=".to_string()],
        "{cross:?}"
    );

    let spdi = checker.check_confluence(&corpus, ConfluenceRelation::Spdi);
    assert_eq!(
        spdi.classes_checked(),
        2,
        "the weaker relation keys each at its own position, so it finds nothing: {spdi:?}"
    );
    assert!(spdi.is_confluent(), "{spdi:?}");
}

/// A pair the checker could not examine must not be reported as a clean corpus.
///
/// The two spellings are one cis allele written in two member orders, 150 kb
/// apart — past `MAX_SEQUENCE_COMPARE_WINDOW`, so `compare_denotations` answers
/// `Indeterminate`. That does not merge, which is byte-identical to what a pair
/// decided *apart* does, so before this was counted the run reported
/// `is_confluent: true`, `classes_checked: 2`, `skipped: []` — indistinguishable
/// from a corpus that was read and found clean.
#[test]
fn a_pair_the_checker_could_not_examine_is_counted_not_assumed() {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_WIDE.1", "ACGT".repeat(50_000));
    let checker = EquivalenceChecker::new(provider);
    let corpus = variants(&[
        "NC_WIDE.1:g.[10del;150010del]",
        "NC_WIDE.1:g.[150010del;10del]",
    ]);

    let report = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);

    assert!(
        report.skipped().is_empty(),
        "both inputs denote a sequence — nothing was dropped: {report:?}"
    );
    assert_eq!(
        report.classes_checked(),
        2,
        "the undecidable pair split a class it could not examine: {report:?}"
    );
    assert_eq!(
        report.undecided_pairs(),
        1,
        "the one comparison came back undecidable: {report:?}"
    );
    assert!(
        !report.is_complete(),
        "a run that could not examine a pair has not covered its corpus: {report:?}"
    );
}

/// A reference that cannot serve the bases is reported as *that*, and never as
/// a clean corpus.
///
/// This is the consumer failure the report exists to prevent: a partly
/// provisioned reference — one contig absent — over which every unreadable
/// description was split into a singleton and every non-confluence among them
/// made invisible, while the report read `is_confluent: true`, `skipped: []`.
/// Both spellings here are well formed and denote one variant, and both
/// normalize to `g.2del`.
#[test]
fn an_unreadable_reference_is_named_rather_than_read_as_clean() {
    // Serves no bases at all: `NC_ABSENT.1` is not in it.
    let checker = EquivalenceChecker::new(MockProvider::new());
    let corpus = variants(&["NC_ABSENT.1:g.2del", "NC_ABSENT.1:g.2delC"]);

    for relation in [
        ConfluenceRelation::CrossAxisSequenceMatch,
        ConfluenceRelation::Spdi,
    ] {
        let report = checker.check_confluence(&corpus, relation);

        assert!(
            !report.is_complete(),
            "the reference could not be read, so the run covered nothing it claims: {report:?}"
        );
        assert!(
            !report.skipped().is_empty(),
            "at least one input has no computable denotation here: {report:?}"
        );
        assert!(
            report
                .skipped()
                .iter()
                .all(|s| s.kind == ConfluenceSkipKind::Unplaceable),
            "{report:?}"
        );
        // Fix the misdiagnosis: the descriptions are well formed and it is the
        // provider that declined, so the skip must say so.
        //
        // **This keys on the classifier, never on the sentence.** Asserting
        // `reason.contains("reference")` — which this test did until the typed
        // decline existed — passes whenever that word happens to appear and
        // breaks silently the moment someone rewords the message, so it tests
        // the prose rather than the classification. Matching the variant makes
        // the assertion immune to any rewording and sensitive to exactly the
        // thing it is about: which class `edit_triples` put the refusal in.
        assert!(
            report
                .skipped()
                .iter()
                .any(|s| matches!(s.decline, Some(TripleDecline::ReferenceUnavailable(_)))),
            "the refusal must be classified against the reference, not the description: {report:?}"
        );
    }
}

/// Every reachable `(class, site)` pair renders a **distinct** sentence, and the
/// count is what it is rather than what is convenient to claim.
///
/// Six refusal sites exist in `edit_triples`. Only `MemberConversion` can be
/// either class — it is the only site holding a `ConversionError`, and only a
/// conversion error can be "want of reference data" — so the reachable pairs
/// number **seven**, not twelve and not six. Enumerated here so that a site
/// added without its own wording (or two sites whose prose collapses together)
/// fails rather than quietly reducing what a report can tell a consumer apart.
#[test]
fn every_reachable_decline_renders_a_distinct_message() {
    use ferro_hgvs::spdi::ConversionError;

    let subject = parse_hgvs("NC_KEY.1:g.3A>G").unwrap();
    let conversion = || ConversionError::MissingReferenceData {
        description: "probe".to_string(),
    };

    // The five sites that only ever classify as `Unrepresentable`, plus the one
    // site that reaches both classes — once under each.
    let reachable = [
        TripleDecline::Unrepresentable(DeclineSite::MultiMoleculeAllele {
            phase: ferro_hgvs::hgvs::variant::AllelePhase::Trans,
        }),
        TripleDecline::Unrepresentable(DeclineSite::NullOrUnknownAllele),
        TripleDecline::Unrepresentable(DeclineSite::EmptyAllele),
        TripleDecline::Unrepresentable(DeclineSite::MemberConversion {
            error: conversion(),
        }),
        TripleDecline::Unrepresentable(DeclineSite::CrossAccession {
            first: "NC_A.1".to_string(),
            second: "NC_B.1".to_string(),
        }),
        TripleDecline::Unrepresentable(DeclineSite::UnresolvedAccession),
        TripleDecline::ReferenceUnavailable(DeclineSite::MemberConversion {
            error: conversion(),
        }),
    ];

    // The counts below are the *properties*, not a restated literal: a seventh
    // pair that renders like an existing one, or two sites sharing a name, fails
    // here. What guarantees a newly *added* site reaches this list at all is not
    // a count but the wildcard-free matches in `into_reason` and
    // `DeclineSite::as_str`, which make an unhandled site a compile error.
    let rendered: std::collections::BTreeSet<String> = reachable
        .iter()
        .map(|d| d.clone().into_reason(&subject))
        .collect();
    assert_eq!(
        rendered.len(),
        7,
        "two decline pairs render the same sentence, so a report cannot tell them apart: {rendered:#?}"
    );

    // The six site names are distinct too — that is the machine-readable half,
    // and it is what a consumer matches on.
    let sites: std::collections::BTreeSet<&str> =
        reachable.iter().map(|d| d.site().as_str()).collect();
    assert_eq!(sites.len(), 6, "{sites:?}");
}

/// A skip is a skip because the relation could not place the input — never
/// because the input is malformed or the normalizer refused it.
///
/// `g.pterdel` is the sharpest form of that: it normalizes perfectly well (to
/// `g.2del`), and it still has no computable denotation, because `pter` is a
/// landmark of the assembled chromosome rather than a coordinate on the
/// sequence. So it is `Unplaceable`, forms no class, and is excluded from
/// `classes_checked` — which is the half of the accounting the two
/// `ConfluenceSkipKind`s exist to keep apart. A
/// `NormalizationDeclined` skip is excluded the *other* way: its class is
/// formed and counted, and only its output is missing. That arm is stated as
/// contract on `ConfluenceReport::classes_checked` rather than reproduced here,
/// because no description this provider can place also fails to normalize.
#[test]
fn an_unplaceable_input_forms_no_class_even_when_it_normalizes() {
    let checker = EquivalenceChecker::new(genomic_provider());
    let corpus = variants(&["NC_KEY.1:g.3A>G", "NC_KEY.1:g.pterdel"]);

    let report = checker.check_confluence(&corpus, ConfluenceRelation::CrossAxisSequenceMatch);

    assert_eq!(report.skipped().len(), 1, "{report:?}");
    assert_eq!(
        report.skipped()[0].input,
        "NC_KEY.1:g.pterdel",
        "{report:?}"
    );
    assert_eq!(
        report.skipped()[0].kind,
        ConfluenceSkipKind::Unplaceable,
        "it was never placed — this is not a normalization failure: {report:?}"
    );
    assert_eq!(
        report.classes_checked(),
        1,
        "only the placeable input formed a class: {report:?}"
    );
    assert!(!report.is_complete(), "{report:?}");
}
