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
    /// Default-prefixed form for bare fragments (#4). When `Some`, this is
    /// what ferro is run against — `input` stays as the verbatim spec text.
    #[serde(default)]
    input_prefixed: Option<String>,
    current: String,
    /// `None` means the spec rejects this input (sentinel for #2).
    spec_expected: Option<String>,
    status: String,
    coordinate_system: String,
    source_kind: String,
    source_paths: Vec<String>,
    #[serde(default)]
    working_group: Option<String>,
    #[serde(default)]
    todo: Option<String>,
    /// `Some(true)` means the row needs a real reference sequence to
    /// evaluate (e.g. §2.1 3'-rule shifting). Skipped until #82 lands.
    #[serde(default)]
    requires_reference: Option<bool>,
    /// Warning codes ferro currently emits for this row, sorted alphabetically.
    /// Defaults to `vec![]` for rows with no warnings.
    #[serde(default)]
    expected_warnings: Vec<String>,
}

fn fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/fixtures/grammar/hgvs_spec_normalization.json");
    p
}

fn observe(normalizer: &Normalizer<MockProvider>, input: &str) -> (String, Vec<String>) {
    match parse_hgvs(input) {
        Err(e) => (format!("parse error: {e}"), Vec::new()),
        Ok(v) => match normalizer.normalize_with_diagnostics(&v) {
            Err(e) => (format!("normalize error: {e}"), Vec::new()),
            Ok(n) => {
                let mut codes: Vec<String> =
                    n.warnings.iter().map(|w| w.code().to_string()).collect();
                codes.sort();
                codes.dedup();
                (format!("{}", n.result), codes)
            }
        },
    }
}

/// True when `s` is a bare coordinate reference (coord prefix + only position
/// characters, no edit) — mirrors the harvester's drop predicate so this test
/// pins the *outcome* (a clean fixture) independent of the generator internals.
fn is_bare_coordinate(s: &str) -> bool {
    let core = s.rsplit(':').next().unwrap_or(s);
    match core.get(0..2) {
        Some("c." | "g." | "m." | "n." | "o." | "p." | "r.") => core[2..]
            .chars()
            .all(|ch| matches!(ch, '0'..='9' | '_' | '+' | '-' | '*' | '?')),
        _ => false,
    }
}

/// #955 harvester cleanup. The generated fixture must contain only real variant
/// strings: spec-prose position references (`c.5690`), `<code class="invalid">`
/// HTML remnants, and syntax-template fragments are dropped, while variants the
/// spec split across `<code class="spotN">` highlights are reassembled intact.
#[test]
fn harvester_yields_only_real_variants() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    // No HTML remnant leaks through — no real HGVS string contains `<` or `"`.
    for row in &fx.rows {
        assert!(
            !row.input.contains('<') && !row.input.contains('"'),
            "HTML remnant leaked into fixture: {:?}",
            row.input
        );
    }

    // Outcome check over *every* parse-error row: an edit-less coordinate ferro
    // rejects is spec prose, not a variant, and is dropped. This must scan all
    // rows, not fixed examples — whether a bare-shaped coordinate is dropped
    // depends on ferro's parse result (`!parse_ok`), not shape alone: e.g.
    // `g.12345678` is bare-shaped but PARSES, so it legitimately survives. The
    // predicate below only pins the outcome; the predicate's char-set is
    // unit-tested independently in the generator (`is_edit_less_coordinate`
    // tests). The generator's drop is also gated on `ov.is_none()` — an explicit
    // override deliberately keeps its row; no current override is a bare
    // coordinate, so this asserts over every parse-error row; if one is ever
    // added, exempt overridden inputs here.
    for row in &fx.rows {
        if row.status == "parse-error" {
            assert!(
                !is_bare_coordinate(&row.input),
                "bare coordinate reference should have been dropped, not left as parse-error: {:?}",
                row.input
            );
        }
    }

    // Variants the spec split across a `<code class="spotN">` highlight are
    // reassembled and parse cleanly (previously lost as severed fragments).
    let by_input: std::collections::HashMap<&str, &Row> =
        fx.rows.iter().map(|r| (r.input.as_str(), r)).collect();
    for recovered in [
        "LRG_199t1:c.11T>G",
        "LRG_199p1:p.(Val25Gly)",
        "g.123_456dup",
    ] {
        let row = by_input
            .get(recovered)
            .unwrap_or_else(|| panic!("reassembled variant missing from fixture: {recovered}"));
        assert_eq!(
            row.status, "preserved",
            "reassembled {recovered} should round-trip (preserved)"
        );
    }
}

#[test]
fn pinned_v21_normalization_behavior() {
    crate::common::spec_fixture::ensure_spec_fixture();
    let text = std::fs::read_to_string(fixture_path()).expect("read fixture");
    let fx: Fixture = serde_json::from_str(&text).expect("parse fixture");

    let normalizer = Normalizer::new(MockProvider::new());
    let mut diffs: Vec<String> = Vec::new();
    let mut skipped_needs_ref = 0usize;
    let mut skipped_requires_ref = 0usize;
    let mut tested = 0usize;

    for row in &fx.rows {
        if row.status == "needs-reference" {
            skipped_needs_ref += 1;
            continue;
        }
        if row.requires_reference == Some(true) {
            skipped_requires_ref += 1;
            continue;
        }
        let target = row.input_prefixed.as_deref().unwrap_or(&row.input);
        let (observed, mut observed_warnings) = observe(&normalizer, target);
        observed_warnings.sort();
        observed_warnings.dedup();
        let mut expected_warnings = row.expected_warnings.clone();
        expected_warnings.sort();
        expected_warnings.dedup();
        if observed != row.current || observed_warnings != expected_warnings {
            diffs.push(format!(
                "  input            : {}\n    target           : {}\n    expected         : {}\n    observed         : {}\n    expected_warnings: {:?}\n    observed_warnings: {:?}\n    status           : {}",
                row.input, target, row.current, observed, expected_warnings, observed_warnings, row.status,
            ));
        }
        tested += 1;
    }

    eprintln!(
        "hgvs_spec_normalization: tested {tested}, skipped(needs-reference) {skipped_needs_ref}, skipped(requires-reference) {skipped_requires_ref}"
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
