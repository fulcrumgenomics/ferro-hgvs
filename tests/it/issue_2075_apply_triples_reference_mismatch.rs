//! #2075: an `apply_triples` failing **after a successful fetch** must report
//! [`EquivalenceLevel::Indeterminate`], not a decided negative.
//!
//! This is the last of the three reference-level obstacles that used to reach a
//! decided `NotEquivalent`. Its two siblings are closed:
//!
//! * #1989/#2053 — the sequence rung's **window fetch** failing after both sides
//!   convert to SPDI ([`compare_triples`]'s missing reference bases);
//! * #2056/#2069 — a member's own **HGVS→SPDI conversion** needing the reference
//!   and failing (`edit_triples`).
//!
//! This one is one step further in: the window fetch **succeeds** and both sides
//! convert, but `apply_triples` still declines — most commonly because a stated
//! reference base contradicts the base the reference actually carries. That is
//! `compare_triples`'s catch-all arm, and before this fix it fell through to the
//! ladder's `NotEquivalent` tail: a decided negative asserted for a pair whose
//! resulting sequences were never reconstructed.
//!
//! # The measured wrong answer
//!
//! `SyntheticBuilder::genomic` pads its core with 256 `ACGT` bases on each side,
//! so `NC_TEST.1` is **524** bases and positions 20 and 21 are *served*, not
//! out-of-bounds. Position 20 reads `T` and position 21 reads `A`, so
//! `g.20A>G` and `g.21C>T` each state a reference base the reference contradicts.
//! Both SPDIs convert without consulting the provider (a substitution states its
//! own bases), so neither takes the `edit_triples` route (#2056) — the fetch
//! succeeds and the only thing left to decline is the apply. Before this fix the
//! checker answered `NotEquivalent`; the honest verdict is `Indeterminate`.
//!
//! # Scope
//!
//! This moves **only** the reference-mismatch decline. A genuinely different
//! served pair, and a genuinely equivalent served pair, both keep their decided
//! verdicts, and so does the one decline `apply_triples` raises that really is a
//! fact about the description — an allele whose members overlap.
//!
//! # #2104: the other half of the split, and why it was unpinned
//!
//! `ApplyDecline` began as a two-way split, with everything that is not a
//! reference mismatch collapsed into one `Unappliable`. Measured over the whole
//! suite, a `panic!()` in the arm that routed it was never reached
//! (`11405 tests run: 11405 passed`), so both mutations that matter were green:
//! over-mapping it onto the reference-obstacle verdict, and reversing the
//! precedence between the two.
//!
//! The reason is not that the shape is unreachable — it is that only one entry
//! point was ever driven with it. `EquivalenceChecker::check` **normalizes
//! before the sequence rung**, and normalization repairs the overlap: measured,
//! `NC_TEST.1:g.[267_270AT[4];270_276del]` normalizes to the single-member
//! `NC_TEST.1:g.270_272del`, and one member cannot overlap anything. Every one
//! of `issue_1244_equivalence_overlap_panic`'s tests uses `check`, which is why
//! this PR's original claim that they preserved the route was wrong — they pass
//! with the arm over-mapped, and they never reach it.
//!
//! `compare_denotations` compares the pair **as written**, so it does reach it.
//! `an_overlapping_cis_allele_is_a_decided_negative_not_a_reference_obstacle`
//! below is that pin, verified by the same `panic!()` probe: with the arm armed,
//! exactly one test in the suite fires, and it is that one.

use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::parse_hgvs;

use crate::common::synthetic::SyntheticBuilder;

/// The core `SyntheticBuilder::genomic` wraps in 256 `ACGT` pad bases per side,
/// so `NC_TEST.1` is 524 bases long. All positions probed below are served.
const CORE: &str = "ACGTACGTACGT";

fn checker() -> EquivalenceChecker<ferro_hgvs::JsonProvider> {
    EquivalenceChecker::new(SyntheticBuilder::genomic(CORE).build())
}

/// The primary reproducer. Both variants state a reference base the served
/// reference contradicts, so `apply_triples` declines both after the fetch
/// succeeds. Nothing was reconstructed, so the verdict must not be decided.
#[test]
fn both_sides_contradict_the_served_reference_and_are_indeterminate() {
    let checker = checker();

    // g.20 reads T (stated A) and g.21 reads A (stated C): both contradicted.
    let a = parse_hgvs("NC_TEST.1:g.20A>G").unwrap();
    let b = parse_hgvs("NC_TEST.1:g.21C>T").unwrap();

    for (label, result) in [
        (
            "compare_denotations",
            checker.compare_denotations(&a, &b).unwrap(),
        ),
        ("check", checker.check(&a, &b).unwrap()),
    ] {
        assert_eq!(
            result.level,
            EquivalenceLevel::Indeterminate,
            "{label}: a stated-ref-base contradiction is 'could not tell', not a decided \
             negative; got {:?}",
            result.level,
        );
        assert!(
            !result.level.is_decided(),
            "{label}: an unreconstructable comparison is not a decided verdict",
        );
        assert!(
            !result.is_equivalent(),
            "{label}: and it is not a positive one either",
        );
    }
}

/// One side contradicts the reference while the other applies cleanly — a
/// different fact from neither applying, but still one where a resulting
/// sequence could not be built, so still `Indeterminate`.
#[test]
fn one_side_contradicting_the_reference_is_indeterminate() {
    let checker = checker();

    // g.20A>G contradicts (g.20 reads T); g.21A>G is faithful (g.21 reads A).
    let contradicts = parse_hgvs("NC_TEST.1:g.20A>G").unwrap();
    let faithful = parse_hgvs("NC_TEST.1:g.21A>G").unwrap();

    let result = checker
        .compare_denotations(&contradicts, &faithful)
        .unwrap();
    assert_eq!(
        result.level,
        EquivalenceLevel::Indeterminate,
        "one side could not be reconstructed, so the verdict is not decided; got {:?}",
        result.level,
    );
    assert!(!result.level.is_decided());
}

/// The discriminating control, and why this fix does not over-map: when both
/// stated reference bases are faithful, two variants that produce genuinely
/// different sequences are reconstructed and compared, and that decided negative
/// must survive unchanged.
#[test]
fn a_faithful_genuinely_different_pair_stays_a_decided_negative() {
    let checker = checker();

    // g.21 reads A, g.22 reads C — both faithful, and the two edits differ.
    let a = parse_hgvs("NC_TEST.1:g.21A>G").unwrap();
    let b = parse_hgvs("NC_TEST.1:g.22C>G").unwrap();

    assert_eq!(
        checker.compare_denotations(&a, &b).unwrap().level,
        EquivalenceLevel::NotEquivalent,
        "a reconstructed, genuinely different pair is still `NotEquivalent`",
    );
    assert_eq!(
        checker.check(&a, &b).unwrap().level,
        EquivalenceLevel::NotEquivalent,
    );
}

/// The other half of the control: two faithful spellings of one edit still reach
/// a decided **positive**. The fix moves only the case the reference could not
/// carry, never a case it could.
#[test]
fn a_faithful_equivalent_pair_stays_a_decided_positive() {
    let checker = checker();

    // g.21 reads A: both spellings are the same faithful edit A->G.
    let a = parse_hgvs("NC_TEST.1:g.21A>G").unwrap();
    let b = parse_hgvs("NC_TEST.1:g.21delAinsG").unwrap();

    let result = checker.compare_denotations(&a, &b).unwrap();
    assert!(
        result.is_equivalent(),
        "a reconstructed, equivalent pair still reaches a decided positive; got {:?}",
        result.level,
    );
}

// ---------------------------------------------------------------------------
// #2104: the other half of the same split.
// ---------------------------------------------------------------------------

/// The 30-base core #1244 uses, whose core positions 11..=20 are an `(AT)`
/// tandem repeat. Under `SyntheticBuilder::genomic` core `n` sits at HGVS
/// `256 + n`, so the repeat spans `g.267`..=`g.276`.
const OVERLAP_CORE: &str = "CAGCTAGCTGATATATATATGCGCGCGCGC";

/// A repeat expansion over core 11..=14 followed by a deletion over core
/// 14..=20: both members claim core position 14, so the allele's resulting
/// sequence is undefined. #1244's own reproducer.
const OVERLAPPING_ALLELE: &str = "NC_TEST.1:g.[267_270AT[4];270_276del]";

/// The description-level decline, driven through the ladder.
///
/// `compare_denotations` is load-bearing here and `check` cannot substitute for
/// it: `check` normalizes first and normalization collapses this allele to one
/// member (see the module docs), which is why the arm had never been reached.
#[test]
fn an_overlapping_cis_allele_is_a_decided_negative_not_a_reference_obstacle() {
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic(OVERLAP_CORE).build());
    let overlapping = parse_hgvs(OVERLAPPING_ALLELE).unwrap();
    let elsewhere = parse_hgvs("NC_TEST.1:g.257C>A").unwrap();

    let result = checker
        .compare_denotations(&overlapping, &elsewhere)
        .unwrap();
    assert_eq!(
        result.level,
        EquivalenceLevel::NotEquivalent,
        "an allele whose members overlap denotes no single sequence, which the \
         description itself says — so this stays a decided negative and must not be \
         swept up with the reference-level declines this PR moves; got {:?}",
        result.level,
    );
}
