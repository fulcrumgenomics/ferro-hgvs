//! The runnable conformance census, exercised as an instrument (#1890).
//!
//! `tests/it/spec_conformance_axis.rs` is the *gate* — it pins the census over
//! the spec corpus and fails when a number moves. This module tests the
//! **instrument** the gate and the two companion pieces (#1891 stability harness,
//! #1892 consumer self-check) all consume: the census run as a library call that
//! emits structured, machine-readable counts a consumer can diff without editing
//! a test.
//!
//! The gate consumes the same `measure` this module does — it passes
//! [`Equivalence::ExactString`], preserving its pins — so a separate cross-check
//! that the two agree is unnecessary: the gate's own `assert_eq!(census, pinned)`
//! is that cross-check.
//!
//! Each full census run normalizes the whole corpus (~seconds), so the tests that
//! need a real run are kept few and the properties that are pure — the burn-down
//! classification and the seam-oracle prefix rule — are tested against hand-built
//! inputs so they cost nothing and do not race the process environment.
//!
//! The one property that needs a real *process* rather than a real run is the
//! seam-oracle **refusal**, since the thing that can be wrong is reading the live
//! environment. It gets a subprocess with a controlled environment, which costs a
//! process start and not a census — the refusal happens before any measurement,
//! and the test asserts that too.
//!
//! This module is named in `ci.yml`'s `ORACLE_EXCLUDE` for the same reason
//! `spec_conformance_axis` is: it measures the census, and an armed oracle makes a
//! census read better than the truth. Since the refusal moved inside `run_census`
//! it would simply be red there instead, which is the louder half of the same
//! fact.

use std::collections::{BTreeMap, BTreeSet};

use ferro_hgvs::conformance::census::{
    compare_reports, oracle_armed, run_census, Census, CensusReport, CorpusShape, CorpusSource,
    DirectionOfVirtue, Equivalence, Movement,
};
use ferro_hgvs::conformance::completeness::CaptureCounts;
use ferro_hgvs::ShuffleDirection;

/// A hand-built report for the pure tests, so they neither run a full census nor
/// touch the environment.
fn synthetic_report(census: Census) -> CensusReport {
    CensusReport {
        corpus: "spec".to_string(),
        direction: "3prime".to_string(),
        equivalence: "sequence".to_string(),
        shape: CorpusShape::default(),
        capture: CaptureCounts {
            attempted: 1,
            succeeded: 1,
            dropped: 0,
            dropped_by_reason: BTreeMap::new(),
            allowance: None,
        },
        census,
        confluence_unkeyable_outputs: 0,
    }
}

/// The instrument produces a machine-readable report whose corpus completeness is
/// accounted, whose confluence decision under a real relation flags the outputs it
/// could not key, and which round-trips through JSON.
///
/// The `CaptureCounts` stamped beside the census is the discipline from
/// `CLAUDE.md`'s "a generator must account for what it dropped": `attempted ==
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
///
/// It costs one more full census run than the module strictly needs. That is the
/// price of asserting a completeness property against a population large enough
/// for completeness to mean anything.
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

/// The **prefix rule** [`oracle_armed`] applies, over hand-written key lists.
///
/// This is a string filter and nothing more, which is why it is no longer named
/// after the refusal: the test that was called
/// `a_seam_oracle_in_the_environment_is_refused` asserted only on this helper, so
/// its name claimed a refusal it never reached — neither `oracle_armed_in_env`,
/// which is the half that reads a real environment, nor the non-zero exit. That
/// refusal now has its own test below.
///
/// Kept pure over an iterator of names so it costs nothing and does not mutate the
/// process environment, which races every other test in the binary under
/// `cargo test` (`tests/it/common/cis_apply_oracle.rs` records that lesson at
/// length).
#[test]
fn the_oracle_detector_matches_the_ferro_assert_prefix_and_nothing_else() {
    let armed = [
        "PATH".to_string(),
        "FERRO_ASSERT_IDEMPOTENT".to_string(),
        "HOME".to_string(),
    ];
    assert_eq!(
        oracle_armed(armed),
        Some("FERRO_ASSERT_IDEMPOTENT".to_string()),
        "an armed FERRO_ASSERT_* variable must be detected and named"
    );

    let clean = [
        "PATH".to_string(),
        "HOME".to_string(),
        "FERRO_PARTITION".to_string(),
    ];
    assert_eq!(
        oracle_armed(clean),
        None,
        "no FERRO_ASSERT_* variable means the census may run"
    );

    // The lexicographic tie-break, so the message names a stable variable when
    // several are armed.
    assert_eq!(
        oracle_armed([
            "FERRO_ASSERT_REPARSE".to_string(),
            "FERRO_ASSERT_IN_BOUNDS".to_string(),
        ]),
        Some("FERRO_ASSERT_IN_BOUNDS".to_string()),
        "with several armed, the lexicographically first is named"
    );
}

/// The census **refuses** to run with a seam oracle armed, exercised end to end:
/// a real process, a real environment, the real entry point, and the exit status.
///
/// A panicking row contributes no output, so its family silently converges and the
/// census reads BETTER than the truth. That is the failure this refusal exists to
/// prevent, and none of it is observable from the pure prefix filter above — the
/// parts that can actually be wrong are `oracle_armed_in_env` reading the live
/// environment, `run_census` consulting it before doing any work, and the command
/// exiting non-zero instead of writing a flattering artifact.
///
/// A **subprocess**, rather than `set_var` in this process: the environment is
/// process-global and `cargo test` runs the suite as threads in one process, so
/// arming an oracle here would arm it for every concurrently-executing test. The
/// child gets a controlled environment for free.
///
/// It also asserts the refusal happened *before* measuring — no census summary is
/// printed — because a refusal that runs the whole census first would still write
/// the flattering numbers to stderr.
#[cfg(feature = "benchmark")]
#[test]
fn a_seam_oracle_in_the_environment_is_refused() {
    let output = std::process::Command::new(env!("CARGO_BIN_EXE_ferro-benchmark"))
        .args(["conformance", "census", "--equivalence", "sequence"])
        .env("FERRO_ASSERT_IDEMPOTENT", "1")
        .output()
        .expect("the ferro-benchmark binary runs");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        output.status.code(),
        Some(2),
        "an armed oracle must be a non-zero refusal, not a census. stderr was:\n{stderr}"
    );
    assert!(
        stderr.contains("FERRO_ASSERT_IDEMPOTENT"),
        "the refusal must NAME the variable, so the reader knows what to unset:\n{stderr}"
    );
    assert!(
        stderr.contains("seam oracle armed"),
        "the refusal must say why it refused:\n{stderr}"
    );
    assert!(
        !stderr.contains("spec conformance census"),
        "the refusal must come BEFORE the measurement — a census printed and then \
         refused has already published the flattering numbers:\n{stderr}"
    );
}

/// `--compare` reports per-counter movement with the direction of virtue attached,
/// and never collapses to a single verdict — some counters should fall and some
/// should rise, so an average would report a correctness improvement as a
/// regression.
///
/// Pure: it perturbs a hand-built census, so it exercises the classification
/// without a full run.
#[test]
fn compare_reports_movement_with_the_direction_of_virtue() {
    let before = synthetic_report(Census {
        outputs_denoting_no_sequence: 3,
        converged: 100,
        coding_axis_separation_two_or_more_merges: 5,
        ..Census::default()
    });
    let after = synthetic_report(Census {
        outputs_denoting_no_sequence: 1, // a failure counter falls: IMPROVED
        converged: 101,                  // converged rises: IMPROVED
        coding_axis_separation_two_or_more_merges: 7, // an instrument moves: NEUTRAL
        ..Census::default()
    });

    let movements = compare_reports(&before, &after);

    let denoting = movements
        .iter()
        .find(|m| m.name == "outputs_denoting_no_sequence")
        .expect("the perturbed failure counter is reported");
    assert_eq!(denoting.virtue, DirectionOfVirtue::LowerIsBetter);
    assert_eq!(denoting.movement, Movement::Improved);

    let converged = movements
        .iter()
        .find(|m| m.name == "converged")
        .expect("converged is reported");
    assert_eq!(converged.virtue, DirectionOfVirtue::HigherIsBetter);
    assert_eq!(converged.movement, Movement::Improved);

    // The instrument counter is neutral: a movement there is a population change,
    // not a verdict, so it must never read as WORSE — but it must not read as
    // `Unchanged` either, which is what it did until #2063's review. 5 -> 7 is a
    // movement; `Unchanged` describes the value, and the value changed.
    let instrument = movements
        .iter()
        .find(|m| m.name == "coding_axis_separation_two_or_more_merges")
        .expect("the coding-axis instrument is reported");
    assert_eq!(instrument.virtue, DirectionOfVirtue::Neutral);
    assert_eq!(
        instrument.movement,
        Movement::Moved,
        "a neutral counter that moved 5 -> 7 must report a movement, not `Unchanged`"
    );

    // A counter nobody perturbed really did not move.
    let still = movements
        .iter()
        .find(|m| m.name == "guard_violations")
        .expect("guard_violations is reported");
    assert_eq!(still.movement, Movement::Unchanged);
}

/// The burn-down reports the two **completeness** figures — the ones that exist so
/// a partial run cannot read as a complete one — and not only the census counters.
///
/// Without them a release that doubled its unkeyable outputs, or whose corpus
/// enumeration halved, showed up only indirectly (via `underdetermined`) or not at
/// all, and the comparison path is the one #1891 consumes.
#[test]
fn the_burn_down_reports_completeness_as_well_as_conformance() {
    let mut before = synthetic_report(Census::default());
    before.capture = CaptureCounts {
        attempted: 100,
        succeeded: 90,
        dropped: 10,
        dropped_by_reason: BTreeMap::new(),
        allowance: None,
    };
    before.confluence_unkeyable_outputs = 4;

    let mut after = synthetic_report(Census::default());
    after.capture = CaptureCounts {
        attempted: 100,
        succeeded: 60,
        dropped: 40,
        dropped_by_reason: BTreeMap::new(),
        allowance: None,
    };
    after.confluence_unkeyable_outputs = 9;

    let movements = compare_reports(&before, &after);
    let named = |name: &str| {
        movements
            .iter()
            .find(|m| m.name == name)
            .unwrap_or_else(|| panic!("the burn-down must report `{name}`"))
    };

    let unkeyable = named("confluence_unkeyable_outputs");
    assert_eq!((unkeyable.before, unkeyable.after), (4, 9));
    assert_eq!(
        unkeyable.movement,
        Movement::Worse,
        "more outputs the relation cannot key is a LESS complete confluence decision"
    );

    let dropped = named("capture_dropped");
    assert_eq!((dropped.before, dropped.after), (10, 40));
    assert_eq!(
        dropped.movement,
        Movement::Moved,
        "a drop count is a population, so it moves without being graded"
    );
    assert_eq!(named("capture_succeeded").movement, Movement::Moved);
    assert_eq!(
        named("capture_attempted").movement,
        Movement::Unchanged,
        "the enumeration did not change size in this pair"
    );
}

/// Comparing two runs decided under **different** relations grades nothing about
/// confluence.
///
/// Confluence is a question asked under a relation, so `sequence` against
/// `partition` is a change of question. The CLI already printed a NOTE saying so
/// and then printed `converged 11750 -> 11071 WORSE` underneath it, which is the
/// half a reader carries away. The movement stays visible; the verdict goes.
#[test]
fn a_relation_change_leaves_the_confluence_counters_ungraded() {
    let mut before = synthetic_report(Census {
        converged: 11_750,
        split_two: 10,
        declined: 5,
        ..Census::default()
    });
    before.equivalence = "sequence".to_string();
    let mut after = synthetic_report(Census {
        converged: 11_071,
        split_two: 40,
        declined: 9,
        ..Census::default()
    });
    after.equivalence = "partition".to_string();

    let movements = compare_reports(&before, &after);
    let named = |name: &str| {
        movements
            .iter()
            .find(|m| m.name == name)
            .unwrap_or_else(|| panic!("the burn-down must report `{name}`"))
    };

    for counter in [
        "converged",
        "split_two",
        "underdetermined",
        "confluence_unkeyable_outputs",
    ] {
        assert_eq!(
            named(counter).virtue,
            DirectionOfVirtue::Neutral,
            "`{counter}` is decided under the relation, so a relation change gives it no virtue"
        );
    }
    assert_eq!(
        named("converged").movement,
        Movement::Moved,
        "the movement must still be visible — only the grade is withheld"
    );

    // A counter that does NOT depend on the relation keeps its verdict.
    assert_eq!(named("declined").virtue, DirectionOfVirtue::LowerIsBetter);
    assert_eq!(named("declined").movement, Movement::Worse);

    // The same pair under ONE relation is graded, so the neutralization is the
    // relation change and not something about these numbers.
    let mut same = after.clone();
    same.equivalence = "sequence".to_string();
    let graded = compare_reports(&before, &same);
    assert_eq!(
        graded
            .iter()
            .find(|m| m.name == "converged")
            .expect("converged is reported")
            .movement,
        Movement::Worse
    );
}
