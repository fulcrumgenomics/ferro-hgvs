//! Hermetic gate for the generated conformance summary docs (#509).
//!
//! Each corpus's `failure-patterns.md` is a derived view over its `cases.json`
//! disposition annotations + cluster taxonomy. This test re-renders it and
//! fails if the committed file has drifted — the same check the CI step
//! `generate_conformance_summary -- --check` runs, but inside `cargo nextest`
//! so a stale doc (or a dangling cluster reference) is caught without invoking
//! the example binary. Runs without a reference manifest.

use std::path::Path;

use ferro_hgvs::conformance::{biocommons, mutalyzer, summary};

#[test]
fn mutalyzer_failure_patterns_is_current() {
    let content = std::fs::read_to_string("tests/fixtures/mutalyzer-normalize/cases.json")
        .expect("read mutalyzer cases.json");
    let fixture: mutalyzer::Fixture =
        serde_json::from_str(&content).expect("parse mutalyzer cases.json");
    fixture
        .validate_clusters()
        .expect("every mutalyzer cluster reference must resolve");
    summary::check(
        Path::new("tests/fixtures/mutalyzer-normalize/failure-patterns.md"),
        &fixture.to_summary(),
    )
    .expect("regenerate via `cargo run --features dev --example generate_conformance_summary`");
}

#[test]
fn biocommons_failure_patterns_is_current() {
    let content = std::fs::read_to_string("tests/fixtures/biocommons-normalize/cases.json")
        .expect("read biocommons cases.json");
    let fixture: biocommons::Fixture =
        serde_json::from_str(&content).expect("parse biocommons cases.json");
    fixture
        .validate_clusters()
        .expect("every biocommons cluster reference must resolve");
    summary::check(
        Path::new("tests/fixtures/biocommons-normalize/failure-patterns.md"),
        &fixture.to_summary(),
    )
    .expect("regenerate via `cargo run --features dev --example generate_conformance_summary`");
}
