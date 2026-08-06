//! #1448 — strict mode accepted span edits whose spans *intersect* without
//! coinciding, and those descriptions denote no sequence at all.
//!
//! `detect_overlap_conflicts` grouped sub-variants by **exact**
//! `(accession, coord_system, region, start, end)` equality. Two members that
//! overlap only partially — or that nest — land in different groups, so nothing
//! reported them and strict mode accepted them. The independent applier declines
//! them, in either member order:
//!
//! ```text
//! g.[10_14del;12_16del]        ACCEPTED, denotes nothing
//! g.[10_14inv;12_16inv]        ACCEPTED, denotes nothing
//! g.[10_14del;10_16del]        ACCEPTED, denotes nothing  (nested)
//! ```
//!
//! Measured over two-member span-edit alleles on four sequences at widths 1-3
//! and separations 0-4: **10,499 of the 25,848 strict mode accepted — 41% —
//! denoted nothing.** With the intersection test, zero.
//!
//! None of the three seam oracles could see this. `FERRO_ASSERT_REPARSE` passes
//! (both spellings parse), `FERRO_ASSERT_IN_BOUNDS` passes (every coordinate is
//! in range), and `FERRO_ASSERT_IDEMPOTENT` passes whenever the output is a
//! fixed point — which most of these are. The applier is not wired in as an
//! oracle, so "denotes no sequence" was invisible.
//!
//! The test restricted to edits whose write footprint **is** their span. That is
//! #1411's principle, and it is what makes the `dup` behaviour fall out instead
//! of needing a special case: a `dup`/`repeat` reads its span and writes at the
//! junction 3' of it, so two overlapping `dup` spans genuinely compose — and no
//! `dup` family appeared anywhere in the measured census.
//!
//! The `dup` exemption is asserted where it can be isolated — as a unit test
//! on `detect_overlap_conflicts` in `src/normalize/overlap.rs`. End to end it
//! cannot be: `detect_insertion_overlaps` reports an overlapping `dup` pair
//! for its own reasons on `main` (#1437), so an integration assertion here
//! would be measuring that detector rather than this one.

use crate::common::cis_apply_oracle::apply;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// 29 bases, with a `TTT` run at 9-11 so a member has somewhere to shift.
const SEQUENCE: &str = "ACGTACGTTTTACGTACGTGGGCCCACGT";

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", SEQUENCE.to_string());
    provider
}

fn strict_accepts(input: &str) -> bool {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::with_config(provider(), NormalizeConfig::strict())
        .normalize(&variant)
        .is_ok()
}

fn lenient(input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::new(provider())
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// The headline, stated as the invariant rather than as four verdicts: strict
/// mode must not accept a description that denotes no sequence.
///
/// Written this way on purpose. Pinning "these four are rejected" would be
/// satisfied by a change that rejected them for an unrelated reason while
/// letting a fifth shape through; tying the verdict to the applier's answer is
/// what makes it a property.
#[test]
fn strict_never_accepts_a_description_that_denotes_no_sequence() {
    for input in [
        // Partial overlap, every combination of the span-writing kinds.
        "TEMPLATE:g.[10_14del;12_16del]",
        "TEMPLATE:g.[10_14delinsAA;12_16del]",
        // Both members `delins`, which the row above does not reach: a `delins`
        // is the only span-writing kind that carries a payload, so a pair of
        // them is the case where two *different* payloads compete for the same
        // bases rather than one payload competing with a removal.
        "TEMPLATE:g.[10_14delinsAA;12_16delinsCC]",
        "TEMPLATE:g.[10_14inv;12_16inv]",
        "TEMPLATE:g.[10_14del;12_16inv]",
        "TEMPLATE:g.[10_14inv;12_16delinsAA]",
        // Nested, which shares no bound at all with its sibling.
        "TEMPLATE:g.[10_14del;10_16del]",
        // A single-base member interior to a span — caught before this fix, and
        // it must stay caught.
        "TEMPLATE:g.[10_14del;12C>G]",
        "TEMPLATE:g.[10_14inv;12C>G]",
    ] {
        assert!(
            apply(SEQUENCE, input).is_none(),
            "fixture error: `{input}` was chosen because it denotes no sequence, \
             but the applier accepted it"
        );
        assert!(
            !strict_accepts(input),
            "`{input}` denotes no sequence, so strict mode must reject it"
        );
    }
}

/// A conflicting allele is returned as authored, not shifted.
///
/// Before the fix, two of these were *shifted* while left overlapping —
/// `g.[10_14del;12_16del]` came back as `g.[11_15del;12_16del]` — which is the
/// pipeline repairing a member of an allele it should have been preserving
/// (#395). Reporting the conflict routes them to the verbatim path.
#[test]
fn a_conflicting_allele_is_returned_as_authored() {
    for input in [
        "TEMPLATE:g.[10_14del;12_16del]",
        "TEMPLATE:g.[10_14del;10_16del]",
    ] {
        assert_eq!(
            lenient(input),
            input,
            "`{input}` conflicts, so it must be preserved rather than shifted"
        );
    }
}

/// Disjoint members are untouched. The intersection test must not widen into
/// members that merely sit near each other.
#[test]
fn disjoint_span_members_are_still_accepted() {
    for input in [
        // Adjacent, sharing no base.
        "TEMPLATE:g.[10_14del;15_16del]",
        // Separated.
        "TEMPLATE:g.[10_14del;20_22del]",
        "TEMPLATE:g.[10_14inv;20_22inv]",
    ] {
        assert!(
            strict_accepts(input),
            "`{input}`'s members share no base, so it must still be accepted"
        );
        assert!(
            apply(SEQUENCE, input).is_some(),
            "`{input}` must denote a sequence"
        );
    }
}

/// Exactly-coincident members are reported once, not twice.
///
/// The pre-existing grouped loop names every member of a coincident group; the
/// intersection pass skips equal-bound pairs so it does not re-report them.
#[test]
fn a_coincident_pair_is_reported_once() {
    use ferro_hgvs::normalize::NormalizationWarning;
    let variant = parse_hgvs("TEMPLATE:g.[10_14del;10_14inv]").expect("parse");
    let outcome = Normalizer::new(provider())
        .normalize_with_diagnostics(&variant)
        .expect("lenient normalization must not reject");
    let warnings = outcome.warnings;
    let conflicts = warnings
        .iter()
        .filter(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }))
        .count();
    assert_eq!(
        conflicts, 1,
        "one coincident pair must yield one conflict warning, got {conflicts}: {warnings:?}"
    );
}

/// …and so is an intersecting pair, which the coincident test cannot show.
///
/// The pass walks unordered pairs (`for i`, then `for j` skipping `i + 1`), so
/// each pair is visited once — but nothing pinned that, and a loop that visited
/// both orderings would report every partial overlap twice while still passing
/// every other test in this file. The count is the only thing that notices.
///
/// Both a partial overlap and a nested pair, because they reach the test
/// through different branches of the interval comparison.
#[test]
fn an_intersecting_pair_is_reported_once_too() {
    use ferro_hgvs::normalize::NormalizationWarning;
    for input in [
        "TEMPLATE:g.[10_14del;12_16del]",
        "TEMPLATE:g.[10_14del;10_16del]",
    ] {
        let variant = parse_hgvs(input).expect("parse");
        let outcome = Normalizer::new(provider())
            .normalize_with_diagnostics(&variant)
            .expect("lenient normalization must not reject");
        let warnings = outcome.warnings;
        let conflicts = warnings
            .iter()
            .filter(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }))
            .count();
        assert_eq!(
            conflicts, 1,
            "`{input}` is one intersecting pair and must yield one conflict warning, \
             got {conflicts}: {warnings:?}"
        );
    }
}
