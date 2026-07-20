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
use std::path::Path;

use ferro_hgvs::conformance::reference_window::WindowProvider;
use ferro_hgvs::conformance::spec_projection;
use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;

/// Committed hermetic reference slice backing the `project-*` dimensions.
const PROJECTION_WINDOWS: &str = "tests/fixtures/grammar/spec_enumeration_windows.json";

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
    expectation: String,
    normative_level: String,
    #[serde(default)]
    expected: Option<String>,
    observed: String,
    status: String,
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
        let projector = spec_projection::load_slice(Path::new(PROJECTION_WINDOWS))
            .ok()
            .map(|windows| {
                let cdot = CdotMapper::from_transcripts(windows.transcripts.iter());
                VariantProjector::new(Projector::new(cdot), windows.to_provider())
            });
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
/// submodule SHA, a dimension, and an RFC 2119 level.
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
        assert!(
            matches!(row.normative_level.as_str(), "must" | "should" | "n/a"),
            "row {} has an unknown normative level {:?}",
            row.id,
            row.normative_level
        );
        assert!(
            matches!(
                row.expectation.as_str(),
                "spec-mandated" | "pinned-baseline"
            ),
            "row {} has an unknown expectation kind {:?}",
            row.id,
            row.expectation
        );
        // A pinned baseline must never masquerade as a spec expectation.
        if row.expectation == "pinned-baseline" {
            assert_ne!(
                row.normative_level, "must",
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
const DIVERGENCE_BUDGET: &[(&str, usize)] = &[
    // Spec-forbidden strings ferro parses and renders anyway.
    ("false-acceptance", 6),
    // Spec names a canonical replacement; ferro accepts the bad string and
    // renders something else. These are the genuine violations.
    ("repair-diverges", 4),
    // Spec names a canonical replacement; ferro rejects the bad string instead
    // of repairing it. Spec-conformant — rejection is always permitted.
    ("rejected-not-repaired", 13),
    // Repairs that need real reference bases; not assertable hermetically.
    ("requires-reference", 10),
    // MUST-level output invariants violated by a string ferro emits.
    ("invariant-violation-must", 0),
    // SHOULD-level (advisory) output invariants. Never a hard failure.
    ("invariant-violation-should", 2),
    // syntax.yaml examples that do not parse into their declared axis.
    ("form-axis-diverges", 3),
    // Projections whose rendered form differs from the one the spec states.
    ("projection-diverges", 0),
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
        let actual = counts.get(status).copied().unwrap_or(0);
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
        .filter(|r| r.status == "invariant-violation-must")
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
                row.normative_level, "should",
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
        if row.expectation != "spec-mandated" || row.normative_level != "must" {
            continue;
        }
        let passing = matches!(
            row.status.as_str(),
            "correctly-rejected" | "repaired" | "rejected-not-repaired" | "form-axis-ok"
                // A projection that already matches the form the spec states.
                | "preserved"
        );
        if !passing {
            continue;
        }
        let Some(observed) = replay(row, &ctx) else {
            continue;
        };
        checked += 1;
        let ok = match row.status.as_str() {
            "correctly-rejected" | "rejected-not-repaired" => {
                observed.starts_with("parse error") || observed.starts_with("normalize error")
            }
            "repaired" => Some(&observed) == row.expected.as_ref(),
            "form-axis-ok" => Some(&observed) == row.expected.as_ref(),
            // The spec writes a projected form with or without its accession;
            // an accession-less statement is compared against the coordinate
            // part only, mirroring the generator (see `projection_matches`).
            "preserved" => row
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
