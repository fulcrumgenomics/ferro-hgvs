//! The conformance census's **run-backed** properties (#2063, carved out #2094).
//!
//! These are the tests from `conformance_census_instrument` that need a *real*
//! census run — each one normalizes the whole spec corpus, which costs tens of
//! seconds to minutes rather than the milliseconds a pure test costs. They live
//! in their own file so the census carve-out can address them by module name:
//! `ci.yml`'s `CENSUS_FILTER` negates this module from the `test` / `test-oracle`
//! shards and the `censuses` job runs it off the optimized archive, exactly as it
//! already does for `spec_conformance_axis` and the other census modules. The
//! pure and subprocess tests stay in `conformance_census_instrument`, which is
//! cheap and stays in the shards — and which keeps the `#[cfg(feature =
//! "benchmark")]` refusal test with the `ferro-benchmark` binary it needs, a
//! binary the soak driver does not build.
//!
//! **Why the split, and why it is a move rather than an optimization (#2094).**
//! Before the carve-out both 60s+ census tests hashed into one `Test` shard, and
//! that shard ran 277s while its other 1,885 tests finished in ~72s — three
//! quarters of the shard's wall was two tests. `--partition hash:K/N` balances
//! shard *counts* and is blind to cost, so re-sharding cannot fix it; a shard's
//! wall can never fall below its slowest single test. Moving the run-backed tests
//! to the job that owns the census does. The tests are otherwise unchanged: they
//! call `run_census` / `CensusReport::build` directly, as they always did.
//! Deduplicating the censuses across tests was tried and dropped — nextest runs
//! each test in its own process, so a shared `LazyLock` would never be reused, and
//! `a_looser_relation_converges_at_least_as_many_families` genuinely needs both a
//! `sequence` and a `partition` census in one process regardless.
//!
//! This module is named in `ci.yml`'s `ORACLE_EXCLUDE` for the same reason
//! `conformance_census_instrument` and `spec_conformance_axis` are: it measures
//! the census, and an armed seam oracle makes a census read better than the truth
//! (a panicking row contributes no output, so its family silently converges).
//! Since the refusal moved inside `run_census` it would be a red `OracleArmed`
//! error here instead, which is the louder half of the same fact.

use std::collections::BTreeSet;

use ferro_hgvs::conformance::census::{run_census, CensusReport, CorpusSource, Equivalence};
use ferro_hgvs::ShuffleDirection;

/// The instrument produces a machine-readable report whose corpus completeness is
/// accounted, whose confluence decision under a real relation flags the outputs it
/// could not key, and which round-trips through JSON.
///
/// The `CaptureCounts` stamped beside the census is the discipline from
/// `CONTRIBUTING.md`'s "a generator must account for what it dropped": `attempted ==
/// succeeded + dropped`, and the corpus really does drop designs, so `dropped` is
/// positive and every drop carries a reason. SPDI has no offset notation, so the
/// corpus's intronic outputs have no `sequence` key — a family carrying one cannot
/// be certified converged, so it is `underdetermined` and the unkeyable output is
/// counted rather than swept into a bucket.
#[test]
fn a_spec_census_run_is_accounted_and_flags_unkeyable_outputs() {
    let run = run_census(
        CorpusSource::Spec,
        ShuffleDirection::ThreePrime,
        Equivalence::Sequence,
    )
    .expect("the hermetic spec census runs with no reference");
    let report = &run.report;

    assert_eq!(report.corpus, "spec");
    assert_eq!(report.direction, "3prime");
    assert_eq!(report.equivalence, "sequence");

    // The corpus shape travels beside the census (the #1460 hazard), and a zero
    // would be a claim about the corpus, so the denominators are non-zero.
    assert!(report.shape.rows > 0, "VACUOUS: the corpus built no rows");
    assert!(report.shape.family_rows > 0, "VACUOUS: no family rows");
    assert!(
        report.census.outputs > 0,
        "VACUOUS: no outputs were produced"
    );

    // Completeness is a claim carried in the artifact, and it is self-consistent.
    assert!(
        report.capture.is_self_consistent(),
        "the capture counts do not add up: {:?}",
        report.capture
    );
    assert_eq!(
        report.capture.succeeded, report.shape.rows,
        "every surviving row must be a captured success"
    );
    assert!(
        report.capture.dropped > 0,
        "the corpus drops designs by construction, so a zero drop count means the ledger \
         was not actually routing the population"
    );
    assert!(
        !report.capture.dropped_by_reason.is_empty(),
        "a drop with no reason is exactly the silent drop the ledger exists to prevent"
    );

    // The allowance is a CEILING, not decoration. It was `at_most(attempted)`,
    // which with `empty_pass_permitted: false` refuses exactly one thing — a run
    // that attempted nothing — so every drop count from 0 to 100% passed while the
    // command's docs called a drop past the allowance a refusal. A bound that
    // cannot be exceeded is a bound nobody is standing behind.
    let allowance = report
        .capture
        .allowance
        .as_ref()
        .expect("the run declares the allowance it accepted its shortfall under");
    assert!(
        allowance.max_dropped < report.capture.attempted,
        "an allowance of {} against {} attempted designs permits dropping EVERYTHING, so it \
         refuses nothing a run could actually do",
        allowance.max_dropped,
        report.capture.attempted
    );
    assert!(
        report.capture.dropped <= allowance.max_dropped,
        "the corpus dropped {} of {} designs, past its own ceiling of {} — either the corpus \
         changed shape or the ceiling is wrong; both need a person",
        report.capture.dropped,
        report.capture.attempted,
        allowance.max_dropped
    );

    // A real relation cannot key every output, and the unkeyable ones are counted
    // rather than folded into a converged family.
    assert!(
        report.confluence_unkeyable_outputs > 0,
        "the corpus produces intronic outputs SPDI cannot key, so a real relation must report \
         some unkeyable outputs; zero means they were silently swept into a bucket"
    );
    assert!(
        report.census.underdetermined > 0,
        "families carrying an unkeyable output cannot be certified converged, so some must be \
         undecided"
    );

    // Every family lands in exactly one confluence bucket — nothing is lost.
    let census = &report.census;
    assert_eq!(
        census.converged
            + census.split_two
            + census.split_three
            + census.split_more
            + census.underdetermined,
        report.shape.family_rows,
        "the confluence buckets must partition the family rows"
    );

    // The artifact round-trips, which is the whole point of a machine-readable
    // report: a consumer stores it and diffs a later one.
    let json = serde_json::to_string(report).expect("serializes");
    let back: CensusReport = serde_json::from_str(&json).expect("deserializes");
    assert_eq!(*report, back);
}

/// The summary's divergence cap, restated here as the *vacuity floor* the test
/// below needs rather than as a value under test.
///
/// It is deliberately not imported: `MAX_DIVERGENCES` is private to the library,
/// and the assertion it guards is "the export is not capped", which needs only
/// *some* number the corpus is known to exceed. If the summary's cap ever rises
/// above this, the worst outcome is that the vacuity guard is looser than it could
/// be — never that a truncated export passes, because the export is compared
/// against the census's own bucket counts and not against this literal.
const MAX_SUMMARY_DIVERGENCES: usize = 12;

/// **Every** divergent family reaches the export, not the first twelve.
///
/// This is the property `--findings` rests on. The cap was sized for a panic
/// message and was carried into a machine-readable export, so two releases' files
/// held the same corpus-order-first twelve however many families moved behind
/// them, and #1891's release-to-release diff read "no change" on a release in
/// which hundreds changed.
///
/// Measured under `partition` deliberately, and that choice is the whole test:
/// `sequence` finds **2** divergent families on this corpus, which is under the
/// cap, so the same assertions there would pass whether or not the export
/// truncates — a green test proving nothing. The vacuity guard below is what
/// stops that happening silently again if the corpus changes shape.
#[test]
fn every_divergent_family_reaches_the_export() {
    let run = run_census(
        CorpusSource::Spec,
        ShuffleDirection::ThreePrime,
        Equivalence::Partition,
    )
    .expect("the hermetic spec census runs with no reference");

    // The three split buckets ARE the divergent families, so this is not a
    // restated literal — it is the same fact read off the census.
    let census = &run.report.census;
    let divergent = census.split_two + census.split_three + census.split_more;
    assert!(
        divergent > MAX_SUMMARY_DIVERGENCES,
        "VACUOUS: the corpus produced only {divergent} divergent families under `partition`, \
         at or under the summary's cap of {MAX_SUMMARY_DIVERGENCES} — so this test cannot \
         tell a complete export from a truncated one and is not checking what it claims"
    );
    assert_eq!(
        run.measured.divergences.len(),
        divergent,
        "the export must carry EVERY divergent family, not the first \
         {MAX_SUMMARY_DIVERGENCES}: a consumer diffing two releases' findings files cannot \
         tell a truncated file from a complete one"
    );

    let ids: BTreeSet<&str> = run
        .measured
        .divergences
        .iter()
        .map(|d| d.id.as_str())
        .collect();
    assert_eq!(
        ids.len(),
        run.measured.divergences.len(),
        "a family must be exported once"
    );
    assert!(
        run.measured.divergences.iter().all(|d| d.outputs.len() > 1),
        "a family with one output is not divergent and must not be exported as one"
    );
}

/// A looser relation can only converge MORE families than a stricter one — never
/// fewer.
///
/// `sequence` (same bases) is looser than `partition` (same bases AND same member
/// partition), so under `sequence` at least as many families reach one bucket.
/// This is the property that makes "confluence under relation R" meaningful; if it
/// failed, the relations would not be ordered by strictness.
#[test]
fn a_looser_relation_converges_at_least_as_many_families() {
    let sequence = CensusReport::build(
        CorpusSource::Spec,
        ShuffleDirection::ThreePrime,
        Equivalence::Sequence,
    )
    .expect("sequence census runs");
    let partition = CensusReport::build(
        CorpusSource::Spec,
        ShuffleDirection::ThreePrime,
        Equivalence::Partition,
    )
    .expect("partition census runs");

    assert!(
        sequence.census.converged >= partition.census.converged,
        "sequence (looser) converged {} families but partition (stricter) converged {} — \
         a looser relation cannot converge fewer",
        sequence.census.converged,
        partition.census.converged,
    );
}

/// The `spdi` relation runs end-to-end and its report round-trips. Its member-key
/// grouping is the one the issue names as already having a backbone
/// (`group_by_spdi_key`), so it is exercised over the real corpus rather than only
/// in the ordering test.
#[test]
fn the_spdi_relation_runs_and_round_trips() {
    let report = CensusReport::build(
        CorpusSource::Spec,
        ShuffleDirection::ThreePrime,
        Equivalence::Spdi,
    )
    .expect("spdi census runs");
    assert_eq!(report.equivalence, "spdi");
    assert!(report.census.outputs > 0, "VACUOUS: no outputs");

    let json = serde_json::to_string(&report).expect("serializes");
    let back: CensusReport = serde_json::from_str(&json).expect("deserializes");
    assert_eq!(report, back);
}
