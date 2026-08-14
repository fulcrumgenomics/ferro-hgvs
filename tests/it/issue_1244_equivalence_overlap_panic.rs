//! Issue #1244: `EquivalenceChecker::check` panicked (slice index out of
//! bounds) instead of returning a verdict, for a cis allele whose members
//! overlap.
//!
//! # The bug
//!
//! `same_resulting_sequence` rebuilds each variant's edited window by splicing
//! its SPDI triples into a copy of the reference. `apply_triples` applies them
//! in **descending** position order, which lets it validate each stated
//! deletion against the *original* reference bases — every already-applied
//! splice sits strictly 3' of the next one, so the span it is about to read is
//! still pristine.
//!
//! That reasoning holds only while the triples are disjoint. Its bounds check
//! compared against `ref_bytes.len()` — the *original* window length — while
//! splicing into `bytes`, whose length shrinks as deletions are applied. Two
//! overlapping triples therefore let a range pass the check and then index past
//! the end of the mutated buffer:
//!
//! ```text
//! window          10 bases
//! triple applied   rel 5..10, insertion ""   -> bytes.len() becomes 5
//! triple next      rel 0..8                  -> 8 <= ref_bytes.len() (10), passes
//!                                              bytes.splice(0..8, ..) on len 5 -> PANIC
//! ```
//!
//! A cis allele of a repeat expansion immediately followed by an adjacent
//! deletion produces exactly that shape: `g.[11_14AT[4];14_20del]` has two
//! members both claiming position 14.
//!
//! # Why a panic is worse than a wrong answer
//!
//! `check` is the documented fallback when normalized strings differ, and it
//! returns `Result<_, FerroError>` — callers reasonably handle failure by
//! matching the error. A panic bypasses that contract entirely, and through the
//! Python bindings it surfaces as `pyo3_runtime.PanicException`, which
//! subclasses `BaseException` and so escapes `except Exception`.
//!
//! # The fix
//!
//! Overlapping triples have no single well-defined resulting sequence, so
//! `apply_triples` declines them (`None`) rather than splicing. `None` already
//! means "cannot decide" to the only caller, which uses the comparison solely
//! to *upgrade* a `NotEquivalent` verdict — so declining is always safe and
//! never invents an equivalence.

use crate::common::synthetic::SyntheticBuilder;
use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::Normalizer;

/// The issue's synthetic reference: 30 nt whose core positions 11..=20 are an
/// `(AT)` tandem repeat. Under `SyntheticBuilder::genomic` core position `n`
/// sits at HGVS `256 + n`, so the repeat spans `g.267`..=`g.276`.
const CORE: &str = "CAGCTAGCTGATATATATATGCGCGCGCGC";

/// `[267_270AT[4];270_276del]` — a repeat expansion over core 11..=14 followed
/// by a deletion over core 14..=20. Both members claim core position 14, so the
/// allele's members overlap and its resulting sequence is undefined.
const OVERLAPPING_ALLELE: &str = "NC_TEST.1:g.[267_270AT[4];270_276del]";

/// A substitution at core position 1, i.e. a different locus entirely.
const ELSEWHERE: &str = "NC_TEST.1:g.257C>A";

fn checker() -> EquivalenceChecker<ferro_hgvs::JsonProvider> {
    EquivalenceChecker::new(SyntheticBuilder::genomic(CORE).build())
}

// ---------------------------------------------------------------------------
// The reported crash.
// ---------------------------------------------------------------------------

#[test]
fn an_overlapping_cis_allele_returns_a_verdict_instead_of_panicking() {
    let result = checker()
        .check(
            &parse_hgvs(OVERLAPPING_ALLELE).unwrap(),
            &parse_hgvs(ELSEWHERE).unwrap(),
        )
        .expect("check must return a verdict, not fail");

    // The two describe different loci, so the honest answer is NotEquivalent.
    // What matters is that a verdict is produced at all.
    assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
    assert!(!result.is_equivalent());
}

#[test]
fn the_crash_is_symmetric_in_argument_order() {
    let result = checker()
        .check(
            &parse_hgvs(ELSEWHERE).unwrap(),
            &parse_hgvs(OVERLAPPING_ALLELE).unwrap(),
        )
        .expect("check must return a verdict, not fail");
    assert_eq!(result.level, EquivalenceLevel::NotEquivalent);
}

#[test]
fn an_overlapping_allele_compared_with_itself_is_still_decided() {
    // Both sides project to the same overlapping triples. The string-equality
    // rung settles this before any sequence reconstruction is attempted, so it
    // must not reach the splice at all.
    let result = checker()
        .check(
            &parse_hgvs(OVERLAPPING_ALLELE).unwrap(),
            &parse_hgvs(OVERLAPPING_ALLELE).unwrap(),
        )
        .expect("check must return a verdict, not fail");
    assert!(result.is_equivalent());
}

// ---------------------------------------------------------------------------
// The decline must be narrow: disjoint alleles still compare equal.
//
// These two assert `is_equivalent()` rather than a specific rung on purpose.
// Which rung settles a pair is the normalizer's business and it moves as
// normalization improves — both of these used to reach the resulting-sequence
// rung and now settle one rung earlier, at `NormalizedMatch`, because the two
// spellings normalize to the same string. Pinning the rung would turn that
// improvement into a failure. What the overlap guard owes them is only that it
// does not decline them into `NotEquivalent`.
//
// The sequence rung itself is still exercised with a disjoint allele — see the
// coincident-insertion tests below, which reach `SequenceMatch`.
// ---------------------------------------------------------------------------

#[test]
fn a_disjoint_cis_allele_is_not_declined() {
    // #1158's case A: a delins against the cis allele that spells out the same
    // edit member by member. If the overlap guard were too eager it would
    // decline this and the verdict would fall back to NotEquivalent.
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic("AGTCAGT").build());
    let result = checker
        .check(
            &parse_hgvs("NC_TEST.1:g.257_263delinsGATTA").unwrap(),
            &parse_hgvs("NC_TEST.1:g.[257A>G;258G>A;260C>T;262_263del]").unwrap(),
        )
        .unwrap();
    assert!(
        result.is_equivalent(),
        "a disjoint allele must stay comparable, got {:?}",
        result.level
    );
}

#[test]
fn members_that_merely_abut_are_not_treated_as_overlapping() {
    // `257_258del` ends where `259T>A` begins — adjacent, not overlapping. The
    // guard keys on a genuine positional clash, so this must still compare
    // equal rather than decline.
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic("AGTCAGT").build());
    let result = checker
        .check(
            &parse_hgvs("NC_TEST.1:g.[257_258del;259T>A]").unwrap(),
            &parse_hgvs("NC_TEST.1:g.257_259delinsA").unwrap(),
        )
        .unwrap();
    assert!(
        result.is_equivalent(),
        "abutting members describe one edit and must stay comparable, got {:?}",
        result.level
    );
}

#[test]
fn pure_insertions_at_one_position_are_not_treated_as_overlapping() {
    // An insertion deletes nothing, so it claims no reference base and cannot
    // clash with anything — not even another insertion at the same interbase
    // point. `apply_triples`' overlap watermark says so explicitly; this pins
    // it, because an off-by-one there (`>=` for `>`) would decline the whole family while
    // every other test in this file still passed.
    //
    // The assertion is deliberately order-agnostic. Coincident insertions are
    // spliced in an order the guard does not define, so which concatenation
    // comes out is an implementation artifact, not the contract. What *is* the
    // contract: `apply_triples` returns a sequence rather than declining — and
    // a decline is observable, because it would leave both orderings at
    // `NotEquivalent`.
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic(CORE).build());
    let allele = parse_hgvs("NC_TEST.1:g.[277_278insA;277_278insT]").unwrap();

    let matched = sequence_rung_matches(
        &checker,
        &allele,
        &["NC_TEST.1:g.277_278insAT", "NC_TEST.1:g.277_278insTA"],
    );

    assert_eq!(
        matched, 1,
        "two insertions at one position must reach the sequence rung and match \
         exactly one concatenation; 0 means `apply_triples` declined them"
    );
}

#[test]
fn three_insertions_at_one_position_are_still_not_overlapping() {
    // "any number of pure insertions at one interbase position" — the running
    // 3' reach never advances past the shared position, so adding members
    // cannot make the guard fire.
    let checker = EquivalenceChecker::new(SyntheticBuilder::genomic(CORE).build());
    let allele = parse_hgvs("NC_TEST.1:g.[277_278insA;277_278insC;277_278insT]").unwrap();

    let matched = sequence_rung_matches(
        &checker,
        &allele,
        &[
            "NC_TEST.1:g.277_278insACT",
            "NC_TEST.1:g.277_278insATC",
            "NC_TEST.1:g.277_278insCAT",
            "NC_TEST.1:g.277_278insCTA",
            "NC_TEST.1:g.277_278insTAC",
            "NC_TEST.1:g.277_278insTCA",
        ],
    );

    assert_eq!(
        matched, 1,
        "three coincident insertions must still be spliced, not declined"
    );
}

/// How many of `spellings` the checker settles against `allele` at the
/// *resulting-sequence* rung.
///
/// The rung matters here, unlike in the disjoint-allele tests above: only a
/// sequence rung runs `apply_triples`, so reaching one is direct evidence the
/// overlap guard let the triples through. Coincident insertions cannot
/// short-circuit at `NormalizedMatch` or `Identical` — a two-member allele and
/// a single `ins` are neither the same string nor normalize to one — so nothing
/// weaker would prove it, and neither of the two textual rungs can satisfy the
/// floor below by accident.
///
/// These are genomic descriptions, which determine exactly one axis, so
/// apply-equality on it is the whole equivalence relation and the verdict is
/// `CrossAxisSequenceMatch`. Written as a floor on the denotational order
/// rather than as an equality: what this helper is counting is "a sequence rung
/// fired", and pinning the exact rung would make it fail again the next time
/// the ladder gains one.
fn sequence_rung_matches(
    checker: &EquivalenceChecker<ferro_hgvs::JsonProvider>,
    allele: &ferro_hgvs::HgvsVariant,
    spellings: &[&str],
) -> usize {
    spellings
        .iter()
        .filter(|spelling| {
            checker
                .check(allele, &parse_hgvs(spelling).unwrap())
                .expect("check must return a verdict")
                .level
                .is_at_least(EquivalenceLevel::SequenceMatch)
        })
        .count()
}

// ---------------------------------------------------------------------------
// Sibling surface: the same allele must not crash the normalizer either.
// ---------------------------------------------------------------------------

#[test]
fn normalizing_the_same_overlapping_allele_does_not_panic() {
    let normalizer = Normalizer::new(SyntheticBuilder::genomic(CORE).build());
    let parsed = parse_hgvs(OVERLAPPING_ALLELE).unwrap();

    // Whether the normalizer accepts or rejects this allele is its own choice,
    // but it must make one — and if it accepts, the result has to satisfy the
    // same contract as any other: `norm(norm(x)) == norm(x)`. Discarding the
    // result would let a successful-but-non-idempotent normalization through.
    let normalized = normalizer
        .normalize(&parsed)
        .expect("the overlapping allele normalizes; a rejection would be a behaviour change")
        .to_string();

    let renormalized = normalizer
        .normalize(&parse_hgvs(&normalized).expect("normalized output must re-parse"))
        .expect("re-normalizing normalized output must succeed")
        .to_string();

    assert_eq!(
        normalized, renormalized,
        "normalization must be idempotent even for an allele whose members overlap"
    );
}
