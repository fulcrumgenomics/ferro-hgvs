//! Driver for the exhaustive, non-redundant HGVS spec test enumeration.
//!
//! **What a green run here does and does not prove.** Every assertion in this
//! file is either a *rejection/repair* assertion (input side) or a
//! *well-formedness* assertion over the string ferro emits. None of it is a
//! correctness oracle: "is this legal HGVS?" is not "is this the right answer".
//! The correctness oracles are the differential corpora
//! (`mutalyzer_normalize_tests.rs`, `biocommons_normalize_tests.rs`,
//! `hgvs_rs_projection_tests.rs`). A green enumeration means ferro's output
//! parses and obeys the MUST-level shape rules — nothing more.
//!
//! The fixture is generated and gitignored; regenerate it with
//! `cargo run --features dev --example generate_spec_enumeration`.
//!
//! Rows whose recorded behaviour diverges from the spec are **not** failures
//! here. They are classified (`repair-diverges`, `false-acceptance`,
//! `invariant-violation-must`, …) and counted against a committed budget, so
//! the suite stays green while the divergence is recorded and any regression
//! or improvement shows up as a budget mismatch.

use std::collections::BTreeMap;
use std::path::PathBuf;

use ferro_hgvs::conformance::reference_window::WindowProvider;
use ferro_hgvs::conformance::spec_projection;
use ferro_hgvs::conformance::{Expectation, NormativeLevel, Status};
use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;

/// Committed hermetic reference slice backing the `project-*` dimensions.
/// Anchored on `CARGO_MANIFEST_DIR` (via the shared helper) so the path
/// resolves regardless of the test's working directory. Unlike the generated
/// fixtures this slice is committed, so it is returned directly — no
/// `ensure_generated_fixture` regeneration step.
fn projection_windows_path() -> PathBuf {
    crate::common::fixture_gen::fixture_path("tests/fixtures/grammar/spec_enumeration_windows.json")
}

#[derive(Debug, Deserialize)]
struct Enumeration {
    rows: Vec<Row>,
}

#[derive(Debug, Deserialize)]
struct Row {
    id: String,
    dimension: String,
    operation: String,
    error_mode: String,
    target: String,
    expectation: Expectation,
    normative_level: NormativeLevel,
    #[serde(default)]
    expected: Option<String>,
    observed: String,
    status: Status,
    spec_citation: String,
}

fn load() -> Enumeration {
    crate::common::spec_enumeration::ensure_spec_enumeration();
    let path = crate::common::spec_enumeration::spec_enumeration_path();
    let text = std::fs::read_to_string(&path).expect("read enumeration fixture");
    serde_json::from_str(&text).expect("parse enumeration fixture")
}

/// Everything a replay needs: the hermetic normalizer, plus the projector built
/// from the committed reference slice (absent only on a checkout that has not
/// been given the slice, in which case projection rows are skipped rather than
/// failed).
struct Replayer {
    normalizer: Normalizer<MockProvider>,
    projector: Option<VariantProjector<WindowProvider>>,
}

impl Replayer {
    fn new() -> Self {
        // The slice is committed, so its absence means only a checkout that has
        // not materialised it — skip the `project-*` rows then. But if the file
        // *is* present, any load failure (malformed JSON, a missing sidecar
        // FASTA) is a corrupt fixture: fail loudly rather than silently drop all
        // projection coverage, which a blanket `.ok()` would have done.
        let windows_path = projection_windows_path();
        let projector = if windows_path.exists() {
            let windows = spec_projection::load_slice(&windows_path).unwrap_or_else(|e| {
                panic!(
                    "committed projection slice {} is present but failed to load \
                     (corrupt fixture?): {e}",
                    windows_path.display()
                )
            });
            let cdot = CdotMapper::from_transcripts(windows.transcripts.iter());
            Some(VariantProjector::new(
                Projector::new(cdot),
                windows.to_provider(),
            ))
        } else {
            None
        };
        Replayer {
            normalizer: Normalizer::new(MockProvider::new()),
            projector,
        }
    }
}

/// Re-run a row's operation against live ferro.
///
/// Returns `None` for rows whose operation cannot be replayed from the fixture
/// alone (`invariant-check` needs the generator's catalog, which lives with the
/// generator so the rules and their spec citations stay in one place).
fn replay(row: &Row, ctx: &Replayer) -> Option<String> {
    let normalizer = &ctx.normalizer;
    if let Some(axis) = row.operation.strip_prefix("project-") {
        let code = axis.chars().next()?;
        let projector = ctx.projector.as_ref()?;
        let variant = parse_hgvs(&row.target).ok()?;
        let pass = spec_projection::project_all_axes(projector, &variant);
        return pass.axes.get(&code).map(|r| r.as_observed());
    }
    match row.operation.as_str() {
        "reject" | "normalize" => Some(match parse_hgvs(&row.target) {
            Err(e) => format!("parse error: {e}"),
            Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
                Err(e) => format!("normalize error: {e}"),
                Ok(n) => format!("{}", n.result),
            },
        }),
        "parse" if row.dimension == "grammar-form" => Some(match parse_hgvs(&row.target) {
            Ok(v) => match v.coordinate_axis() {
                Some(ax) => format!("axis={}", ax.code()),
                None => "axis=none".to_string(),
            },
            Err(e) => format!("parse error: {e}"),
        }),
        "parse" => {
            let cfg = match row.error_mode.as_str() {
                "strict" => ErrorConfig::strict(),
                "lenient" => ErrorConfig::lenient(),
                "silent" => ErrorConfig::silent(),
                _ => return None,
            };
            Some(match parse_hgvs_with_config(&row.target, cfg) {
                Err(e) => format!("parse error: {e}"),
                Ok(r) => {
                    let mut codes: Vec<String> = r
                        .warnings
                        .iter()
                        .map(|w| format!("{:?}", w.error_type))
                        .collect();
                    codes.sort();
                    codes.dedup();
                    if codes.is_empty() {
                        format!("{}", r.result)
                    } else {
                        format!("{} warnings={}", r.result, codes.join(","))
                    }
                }
            })
        }
        _ => None,
    }
}

/// Every replayable row must reproduce its recorded behaviour exactly. This is
/// the drift lock: it catches a change in ferro that the generator recorded but
/// nothing else asserts.
#[test]
fn enumeration_replays_recorded_behavior() {
    let fx = load();
    let ctx = Replayer::new();
    let mut diffs = Vec::new();
    let mut replayed = 0usize;

    for row in &fx.rows {
        let Some(observed) = replay(row, &ctx) else {
            continue;
        };
        replayed += 1;
        if observed != row.observed {
            diffs.push(format!(
                "  id       : {}\n    target   : {}\n    recorded : {}\n    observed : {}",
                row.id, row.target, row.observed, observed
            ));
        }
    }

    eprintln!(
        "spec enumeration: replayed {replayed} of {} rows",
        fx.rows.len()
    );
    assert!(replayed > 300, "replayed too few rows ({replayed})");
    assert!(
        diffs.is_empty(),
        "{} enumeration row(s) drifted. Regenerate if intentional:\n  \
         cargo run --features dev --example generate_spec_enumeration\n\n{}",
        diffs.len(),
        diffs
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// Every row must carry usable provenance: a spec citation pinned to the
/// submodule SHA. The `normative_level` and `expectation` vocabularies are now
/// enforced structurally — an unknown value fails deserialization by naming
/// itself — so the only cross-field invariant left to assert is that a pinned
/// baseline never masquerades as a MUST-level spec expectation.
#[test]
fn every_row_carries_provenance() {
    let fx = load();
    for row in &fx.rows {
        assert!(
            row.spec_citation.contains('@'),
            "row {} has no SHA-pinned citation: {:?}",
            row.id,
            row.spec_citation
        );
        // A pinned baseline must never masquerade as a spec expectation.
        if row.expectation == Expectation::PinnedBaseline {
            assert_ne!(
                row.normative_level,
                NormativeLevel::Must,
                "row {} pins current behaviour but claims MUST-level force",
                row.id
            );
        }
    }
}

/// Divergence budget.
///
/// Each entry is the number of rows currently in a status that records ferro
/// diverging from the spec. The suite stays green, but any change — a
/// regression *or* a fix — trips this and must be re-blessed deliberately.
///
/// Regenerate the numbers with:
///   `cargo run --features dev --example generate_spec_enumeration -- --census`
/// and read `by_status` in the generated fixture.
const DIVERGENCE_BUDGET: &[(Status, usize)] = &[
    // Spec-forbidden strings ferro parses and renders anyway.
    (Status::FalseAcceptance, 6),
    // Spec names a canonical replacement; ferro accepts the bad string and
    // renders something else. These are the genuine violations.
    (Status::RepairDiverges, 4),
    // Spec names a canonical replacement; ferro rejects the bad string instead
    // of repairing it. Spec-conformant — rejection is always permitted.
    (Status::RejectedNotRepaired, 13),
    // Repairs that need real reference bases; not assertable hermetically.
    (Status::RequiresReference, 10),
    // MUST-level output invariants violated by a string ferro emits.
    (Status::InvariantViolationMust, 0),
    // SHOULD-level (advisory) output invariants. Never a hard failure.
    (Status::InvariantViolationShould, 2),
    // syntax.yaml examples that do not parse into their declared axis.
    (Status::FormAxisDiverges, 3),
    // Projections whose rendered form differs from the one the spec states.
    (Status::ProjectionDiverges, 0),
];

/// Does a projected string match the form the spec states? Mirrors
/// `generate_spec_enumeration`'s rule: a fully-qualified spec value is compared
/// verbatim, an accession-less one against the coordinate part only.
fn projection_matches(observed: &str, expected: &str) -> bool {
    if expected.contains(':') {
        observed == expected
    } else {
        observed.rsplit(':').next() == Some(expected)
    }
}

#[test]
fn divergence_budget_is_unchanged() {
    let fx = load();
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for row in &fx.rows {
        *counts.entry(row.status.as_str()).or_default() += 1;
    }
    let mut mismatches = Vec::new();
    for (status, budget) in DIVERGENCE_BUDGET {
        let actual = counts.get(status.as_str()).copied().unwrap_or(0);
        if actual != *budget {
            mismatches.push(format!("  {status}: budget {budget}, actual {actual}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "spec-divergence budget changed. If ferro improved, lower the budget; if it \
         regressed, this is a real defect. Update DIVERGENCE_BUDGET deliberately.\n{}",
        mismatches.join("\n")
    );
}

/// Known passing/expected statuses that are deliberately **not** part of the
/// divergence budget: every row in one of these is a conformant outcome, not a
/// recorded spec divergence, so `divergence_budget_is_unchanged` rightly ignores
/// them. Keeping the list explicit is the point of issue #1107 — a *new* status
/// the generator starts emitting is then neither budgeted nor allowlisted, so
/// `every_observed_status_is_accounted_for` fails and names it, turning "a new
/// status category appeared" into a deliberate, reviewed decision instead of a
/// silent pass. Regenerate the roster of observed statuses with:
///   `cargo run --features dev --example generate_spec_enumeration -- --census`
/// and read `by_status` in the generated fixture.
const KNOWN_PASSING_STATUSES: &[Status] = &[
    // Spec-forbidden string that ferro rejected outright (parse error) — the
    // conformant counterpart of the budgeted `false-acceptance`.
    Status::CorrectlyRejected,
    // Bad string repaired to the canonical form the spec names — the ideal
    // outcome of a repair.
    Status::Repaired,
    // Grammar example that parses into the coordinate axis the spec declares.
    Status::FormAxisOk,
    // Grammar example whose axis the spec does not state; ferro's parsed axis is
    // pinned as a baseline (not a spec matter). Emittable but currently 0 rows.
    Status::FormAxisPinned,
    // Projection that matches the form the spec states — the conformant
    // counterpart of the budgeted `projection-diverges`.
    Status::Preserved,
    // Projection with no spec-stated expectation, pinned as a baseline: ferro
    // rendered a form, declined/unavailable, or errored. None is a spec
    // divergence — the spec mandates no expectation for these.
    Status::ProjectionPinned,
    Status::ProjectionUnavailablePinned,
    Status::ProjectionErrorPinned,
    // Error-mode outcome pinned as ferro policy: the spec says nothing about
    // ferro's error modes, so a per-mode divergence is never a spec expectation.
    Status::ModeDivergencePinned,
];

/// A status is accounted for if it is either budgeted (a tracked divergence in
/// `DIVERGENCE_BUDGET`) or allowlisted (a known passing outcome in
/// `KNOWN_PASSING_STATUSES`).
fn status_is_accounted_for(status: &Status) -> bool {
    DIVERGENCE_BUDGET.iter().any(|(s, _)| s == status) || KNOWN_PASSING_STATUSES.contains(status)
}

/// Returns the distinct observed statuses that are neither budgeted (a tracked
/// spec divergence in `DIVERGENCE_BUDGET`) nor allowlisted (a known passing
/// outcome in `KNOWN_PASSING_STATUSES`), sorted for a stable message. Such a
/// status is a new, unreviewed category — `divergence_budget_is_unchanged` is
/// blind to it, so this is what `every_observed_status_is_accounted_for` gates
/// on.
fn unaccounted_statuses<'a>(statuses: impl IntoIterator<Item = &'a Status>) -> Vec<&'a Status> {
    let mut unknown: Vec<&Status> = statuses
        .into_iter()
        .filter(|status| !status_is_accounted_for(status))
        .collect();
    unknown.sort_unstable();
    unknown.dedup();
    unknown
}

#[test]
fn unaccounted_statuses_flags_a_new_status() {
    // A budgeted status is accounted for; a status neither budgeted nor
    // allowlisted (e.g. a future `projection panicked` defect indicator) must
    // be surfaced by name so it becomes a deliberate decision, not a silent
    // pass. This holds independent of the allowlist's contents.
    let observed = [
        Status::FalseAcceptance,
        Status::from_wire("projection panicked"),
    ];
    let expected = Status::Unknown("projection panicked".to_string());
    assert_eq!(unaccounted_statuses(&observed), vec![&expected]);
}

#[test]
fn budget_and_allowlist_are_disjoint() {
    // `DIVERGENCE_BUDGET` = tracked divergences, `KNOWN_PASSING_STATUSES` =
    // conformant outcomes: the two carve the status vocabulary into disjoint
    // halves. A status in both would blur that intent (and be double-classified),
    // so guard the split explicitly.
    let overlap: Vec<&Status> = KNOWN_PASSING_STATUSES
        .iter()
        .filter(|passing| {
            DIVERGENCE_BUDGET
                .iter()
                .any(|(budgeted, _)| budgeted == *passing)
        })
        .collect();
    assert!(
        overlap.is_empty(),
        "a status is in both DIVERGENCE_BUDGET and KNOWN_PASSING_STATUSES: {overlap:?}"
    );
}

#[test]
fn every_observed_status_is_accounted_for() {
    let fx = load();
    let unknown = unaccounted_statuses(fx.rows.iter().map(|row| &row.status));
    assert!(
        unknown.is_empty(),
        "the enumeration emitted status(es) that are neither in DIVERGENCE_BUDGET \
         nor KNOWN_PASSING_STATUSES: {unknown:?}. A new status category appeared — \
         classify it deliberately: budget it in DIVERGENCE_BUDGET if it records a \
         spec divergence, or add it to KNOWN_PASSING_STATUSES if it is a conformant \
         outcome. (`divergence_budget_is_unchanged` alone would not have caught this.)"
    );
}

/// The acceptance gate from the conformance design: the output-invariant
/// catalog must produce **zero** MUST-level violations over the outputs the
/// repo already blesses. Every MUST violation must trace back to a
/// `false-acceptance` row — a string the spec forbids that ferro accepted —
/// never to a blessed `preserved` output. A violation on a blessed output means
/// the *rule* is wrong, not ferro.
#[test]
fn invariant_catalog_has_no_false_positives_on_blessed_output() {
    let fx = load();
    let offenders: Vec<&Row> = fx
        .rows
        .iter()
        .filter(|r| r.status == Status::InvariantViolationMust)
        .filter(|r| {
            // The generator records the source row's pinned status in `note`;
            // a blessed output is one whose source row is `preserved`.
            parse_hgvs(&r.target).is_ok() && r.expected.as_deref() == Some("no violation")
        })
        .filter(|r| !r.id.contains("output-invariant/"))
        .collect();
    assert!(
        offenders.is_empty(),
        "invariant catalog fired outside the output-invariant dimension: {:?}",
        offenders.iter().map(|r| &r.id).collect::<Vec<_>>()
    );

    // Cross-check the same claim from the other side: no MUST-level invariant
    // row may cite a rule the design lists as advisory.
    for row in fx.rows.iter().filter(|r| r.dimension == "output-invariant") {
        if row.id.contains("/A1-") || row.id.contains("/A2-") {
            assert_eq!(
                row.normative_level,
                NormativeLevel::Should,
                "advisory rule {} must never be MUST-level",
                row.id
            );
        }
    }
}

/// Spec-mandated MUST rows that currently pass must keep passing. This is the
/// ratchet that turns the enumeration into a regression net rather than a
/// snapshot.
#[test]
fn passing_spec_mandated_musts_stay_passing() {
    let fx = load();
    let ctx = Replayer::new();
    let mut regressions = Vec::new();
    let mut checked = 0usize;

    for row in &fx.rows {
        if row.expectation != Expectation::SpecMandated
            || row.normative_level != NormativeLevel::Must
        {
            continue;
        }
        let passing = matches!(
            row.status,
            Status::CorrectlyRejected
                | Status::Repaired
                | Status::RejectedNotRepaired
                | Status::FormAxisOk
                // A projection that already matches the form the spec states.
                | Status::Preserved
        );
        if !passing {
            continue;
        }
        let Some(observed) = replay(row, &ctx) else {
            continue;
        };
        checked += 1;
        let ok = match row.status {
            Status::CorrectlyRejected | Status::RejectedNotRepaired => {
                observed.starts_with("parse error") || observed.starts_with("normalize error")
            }
            Status::Repaired => Some(&observed) == row.expected.as_ref(),
            Status::FormAxisOk => Some(&observed) == row.expected.as_ref(),
            // The spec writes a projected form with or without its accession;
            // an accession-less statement is compared against the coordinate
            // part only, mirroring the generator (see `projection_matches`).
            Status::Preserved => row
                .expected
                .as_ref()
                .is_some_and(|e| projection_matches(&observed, e)),
            _ => true,
        };
        if !ok {
            regressions.push(format!(
                "  {} ({}): expected {:?}, observed {:?}",
                row.id, row.spec_citation, row.expected, observed
            ));
        }
    }

    eprintln!("spec enumeration: {checked} passing spec-mandated MUST rows re-checked");
    assert!(
        regressions.is_empty(),
        "{} spec-mandated MUST assertion(s) regressed:\n{}",
        regressions.len(),
        regressions.join("\n")
    );
}
