//! HGVS v21.0 spec normalization regression test.
//!
//! Companion to issue #84 / #83. The fixture pins ferro's current normalize()
//! output for every variant string in the v21.0 spec. Rows whose `current`
//! diverges from `spec_expected` carry a `todo` link to the #83 audit.
//!
//! Regenerate the fixture: `cargo run --features dev --example generate_spec_fixture`.
//! Verify byte-identical regen:    `cargo run --features dev --example generate_spec_fixture -- --check`.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;
use std::path::PathBuf;

#[derive(Debug, Deserialize)]
struct Fixture {
    rows: Vec<Row>,
}

#[allow(dead_code)] // many fields are present for downstream tooling, not the assertion
#[derive(Debug, Deserialize)]
struct Row {
    input: String,
    current: String,
    spec_expected: String,
    status: String,
    coordinate_system: String,
    source_kind: String,
    source_paths: Vec<String>,
    #[serde(default)]
    working_group: Option<String>,
    #[serde(default)]
    todo: Option<String>,
}

fn fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/fixtures/grammar/hgvs_spec_normalization.json");
    p
}

fn observe(normalizer: &Normalizer<MockProvider>, input: &str) -> String {
    match parse_hgvs(input) {
        Err(e) => format!("parse error: {e}"),
        Ok(v) => match normalizer.normalize(&v) {
            Err(e) => format!("normalize error: {e}"),
            Ok(n) => format!("{n}"),
        },
    }
}

#[test]
fn pinned_v21_normalization_behavior() {
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let normalizer = Normalizer::new(MockProvider::new());
    let mut diffs: Vec<String> = Vec::new();
    let mut skipped_needs_ref = 0usize;
    let mut tested = 0usize;

    for row in &fx.rows {
        if row.status == "needs-reference" {
            skipped_needs_ref += 1;
            continue;
        }
        let observed = observe(&normalizer, &row.input);
        if observed != row.current {
            diffs.push(format!(
                "  input        : {}\n    expected   : {}\n    observed   : {}\n    status     : {}",
                row.input, row.current, observed, row.status,
            ));
        }
        tested += 1;
    }

    eprintln!(
        "hgvs_spec_normalization: tested {tested}, skipped(needs-reference) {skipped_needs_ref}"
    );

    assert!(
        tested > 0,
        "pinned_v21_normalization_behavior exercised no cases (fixture empty or all needs-reference)"
    );

    if !diffs.is_empty() {
        let preview = diffs
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n");
        panic!(
            "{} row(s) drifted from the pinned current behavior. First {} shown.\n\n{}\n\n\
             If this drift is intentional, regenerate the fixture:\n  \
             cargo run --features dev --example generate_spec_fixture",
            diffs.len(),
            diffs.len().min(20),
            preview,
        );
    }
}
