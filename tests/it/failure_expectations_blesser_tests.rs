//! Merge semantics of the failure-expectations blesser
//! (`UPDATE_FAILURE_EXPECTATIONS=1`).
//!
//! The blesser rewrites a committed `*_failure_expectations.json` snapshot
//! from the current parse output. Those snapshots carry hand-curated
//! content — named kind ids, their `category`, their prose `reason` citing
//! the spec, and the input -> kind mapping. None of that is recoverable
//! from a parser error string, so the binding property is: **a bless must
//! never lose human-authored content.** These tests pin that property.

use crate::common::failure_expectations::{bless, FailureExpectations, FixtureCheck};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// A curated snapshot: two named kinds with spec-citing reasons, three
/// inputs mapped onto them.
const CURATED: &str = r##"{
  "kinds": {
    "spec_self_cancelling_allele": {
      "category": "rejected_by_design",
      "reason": "Overlapping del+dup in the same cis-allele that partially undo each other.",
      "tracking": "#115"
    },
    "single_position_insertion": {
      "category": "rejected_by_design",
      "reason": "Insertion must name the two flanking positions."
    }
  },
  "expected_failures": {
    "NM_1.1:c.[1del;1_3dup]": "spec_self_cancelling_allele",
    "NM_1.1:c.[4del;4_6dup]": "spec_self_cancelling_allele",
    "NM_1.1:c.5insA": "single_position_insertion"
  }
}
"##;

/// Scratch directory for this test module's temp snapshots.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("ferro-blesser-tests");
    std::fs::create_dir_all(&dir).expect("create scratch dir");
    let path = dir.join(name);
    let _ = std::fs::remove_file(&path);
    path
}

fn write_curated(name: &str) -> PathBuf {
    let path = scratch(name);
    std::fs::write(&path, CURATED).expect("write curated snapshot");
    path
}

fn read_back(path: &Path) -> FailureExpectations {
    let raw = std::fs::read_to_string(path).expect("read blessed snapshot");
    serde_json::from_str(&raw).expect("parse blessed snapshot")
}

fn check<'a>(failures: &[(&'a str, &str)], total_inputs: usize) -> FixtureCheck<'a> {
    let map: BTreeMap<&'a str, String> = failures
        .iter()
        .map(|(input, err)| (*input, (*err).to_string()))
        .collect();
    FixtureCheck {
        total_inputs,
        failures: map,
    }
}

/// The core case: a bless where some curated inputs still fail, one new
/// input fails, and one curated input now passes. Curated kind ids,
/// categories, reasons, tracking, and the input -> kind mapping for the
/// still-failing inputs must all survive verbatim.
#[test]
fn bless_preserves_curated_kinds_and_mappings() {
    let path = write_curated("mixed.json");
    let report = bless(
        &path,
        &check(
            &[
                // (a) still-failing curated inputs
                ("NM_1.1:c.[1del;1_3dup]", "E3006 SelfCancellingAllele"),
                ("NM_1.1:c.5insA", "E1234 single position insertion"),
                // (b) a new, previously-unseen failure
                (
                    "NM_1.1:c.9zzz",
                    r#"Parse error at position 7: input: "zzz""#,
                ),
                // (c) "NM_1.1:c.[4del;4_6dup]" now passes -> absent here
            ],
            100,
        ),
    );
    let after = read_back(&path);

    // Curated kinds survive byte-for-byte.
    let cancelling = after
        .kinds
        .get("spec_self_cancelling_allele")
        .expect("curated kind `spec_self_cancelling_allele` was destroyed by the bless");
    assert_eq!(
        cancelling.reason,
        "Overlapping del+dup in the same cis-allele that partially undo each other."
    );
    assert_eq!(cancelling.tracking.as_deref(), Some("#115"));
    assert!(
        after.kinds.contains_key("single_position_insertion"),
        "curated kind `single_position_insertion` was destroyed by the bless"
    );

    // Curated mappings for still-failing inputs survive.
    assert_eq!(
        after
            .expected_failures
            .get("NM_1.1:c.[1del;1_3dup]")
            .map(String::as_str),
        Some("spec_self_cancelling_allele")
    );
    assert_eq!(
        after
            .expected_failures
            .get("NM_1.1:c.5insA")
            .map(String::as_str),
        Some("single_position_insertion")
    );

    // The new failure is recorded under a clearly-auto kind needing triage.
    let new_kind = after
        .expected_failures
        .get("NM_1.1:c.9zzz")
        .expect("new failure was not added");
    assert!(
        new_kind.starts_with("auto:"),
        "new kind id `{}` is not distinguishable from curated ids",
        new_kind
    );
    assert_eq!(report.new_inputs, vec!["NM_1.1:c.9zzz".to_string()]);

    // The now-passing input is removed from the mapping but reported...
    assert!(!after
        .expected_failures
        .contains_key("NM_1.1:c.[4del;4_6dup]"));
    assert_eq!(
        report.resolved_inputs,
        vec!["NM_1.1:c.[4del;4_6dup]".to_string()]
    );
    // ...and its kind is still referenced by the other input, so not orphaned.
    assert!(
        report.orphaned_kinds.is_empty(),
        "{:?}",
        report.orphaned_kinds
    );
}

/// A kind whose every input now passes becomes unreferenced. Its prose is
/// human-authored, so the blesser keeps the definition and reports it
/// rather than deleting it.
#[test]
fn bless_retains_and_reports_orphaned_kinds() {
    let path = write_curated("orphan.json");
    let report = bless(
        &path,
        &check(
            &[("NM_1.1:c.5insA", "E1234 single position insertion")],
            100,
        ),
    );
    let after = read_back(&path);

    assert!(
        after.kinds.contains_key("spec_self_cancelling_allele"),
        "orphaned curated kind was deleted, losing its curated reason"
    );
    assert_eq!(
        report.orphaned_kinds,
        vec!["spec_self_cancelling_allele".to_string()]
    );
}

/// Blessing twice with the same input set must be a no-op on disk.
#[test]
fn bless_is_idempotent() {
    let path = write_curated("idempotent.json");
    let failures = [
        ("NM_1.1:c.[1del;1_3dup]", "E3006 SelfCancellingAllele"),
        ("NM_1.1:c.[4del;4_6dup]", "E3006 SelfCancellingAllele"),
        ("NM_1.1:c.5insA", "E1234 single position insertion"),
        (
            "NM_1.1:c.9zzz",
            r#"Parse error at position 7: input: "zzz""#,
        ),
    ];
    bless(&path, &check(&failures, 100));
    let first = std::fs::read_to_string(&path).expect("read first bless");
    let report = bless(&path, &check(&failures, 100));
    let second = std::fs::read_to_string(&path).expect("read second bless");
    assert_eq!(first, second, "second bless mutated the snapshot");
    assert!(report.new_inputs.is_empty());
    assert!(report.resolved_inputs.is_empty());
}

/// If the fixture iteration produced no inputs at all (dataset missing or
/// skipped), a bless would wipe the whole snapshot. Refuse instead.
#[test]
#[should_panic(expected = "refusing to bless")]
fn bless_refuses_when_fixture_produced_no_inputs() {
    let path = write_curated("empty.json");
    bless(&path, &check(&[], 0));
}

/// An unknown field in a curated snapshot is human-authored content the
/// blesser cannot round-trip. Fail loudly rather than silently drop it.
#[test]
#[should_panic(expected = "unknown field")]
fn bless_refuses_snapshot_with_unknown_fields() {
    let path = scratch("unknown_field.json");
    std::fs::write(
        &path,
        r#"{"kinds":{"k":{"category":"malformed_input","reason":"r","note":"curated note"}},
            "expected_failures":{"x":"k"}}"#,
    )
    .expect("write snapshot");
    bless(&path, &check(&[("x", "boom")], 1));
}

/// An input mapped to a kind id that is not defined gets a placeholder
/// definition (so the enforcing test can explain itself) and is reported.
#[test]
fn bless_synthesizes_placeholder_for_undefined_kind() {
    let path = scratch("undefined_kind.json");
    std::fs::write(
        &path,
        r#"{"kinds":{},"expected_failures":{"NM_1.1:c.5insA":"typo_kind_id"}}"#,
    )
    .expect("write snapshot");
    let report = bless(
        &path,
        &check(&[("NM_1.1:c.5insA", "E1234 single position insertion")], 1),
    );
    let after = read_back(&path);
    assert!(after.kinds.contains_key("typo_kind_id"));
    assert_eq!(report.synthesized_kinds, vec!["typo_kind_id".to_string()]);
}

/// Round-trip proof against a copy of the real, hand-curated ClinVar
/// unique snapshot: blessing it with exactly its own expected failures
/// must leave every kind id, category, reason, and tracking untouched.
#[test]
fn bless_roundtrips_real_clinvar_unique_snapshot() {
    // Anchor on CARGO_MANIFEST_DIR so the path resolves regardless of the
    // test's working directory (repo convention). The fixture is a committed
    // file, so it must be present — a silent skip here would turn the most
    // valuable test into a no-op that still reports green.
    let source = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bulk/clinvar_hgvs_unique_failure_expectations.json");
    let raw = std::fs::read_to_string(&source)
        .unwrap_or_else(|e| panic!("read committed snapshot {}: {e}", source.display()));
    let before: FailureExpectations = serde_json::from_str(&raw).expect("parse real snapshot");
    assert!(before.kinds.len() > 50, "fixture unexpectedly small");

    // Work on a copy — never touch the committed fixture.
    let path = scratch("clinvar_unique_copy.json");
    std::fs::write(&path, &raw).expect("write copy");

    // Simulate the real bless: every curated input still fails, with an
    // error string that bears no relation to its curated kind id.
    let inputs: Vec<String> = before.expected_failures.keys().cloned().collect();
    let failures: BTreeMap<&str, String> = inputs
        .iter()
        .map(|i| {
            (
                i.as_str(),
                "Parse error at position 3: some error".to_string(),
            )
        })
        .collect();
    let report = bless(
        &path,
        &FixtureCheck {
            total_inputs: 4_000_000,
            failures,
        },
    );

    let after = read_back(&path);
    assert_eq!(
        after.kinds.len(),
        before.kinds.len(),
        "kind count changed during bless"
    );
    for (id, kind) in &before.kinds {
        let got = after
            .kinds
            .get(id)
            .unwrap_or_else(|| panic!("curated kind `{}` lost", id));
        assert_eq!(got.category, kind.category, "category changed for `{}`", id);
        assert_eq!(got.reason, kind.reason, "reason changed for `{}`", id);
        assert_eq!(got.tracking, kind.tracking, "tracking changed for `{}`", id);
    }
    assert_eq!(after.expected_failures, before.expected_failures);
    assert!(report.new_inputs.is_empty());
    assert!(report.resolved_inputs.is_empty());
    assert!(report.orphaned_kinds.is_empty());

    // Every curated reason string survives verbatim in the rewritten text
    // (the committed file is not in canonical key order, so the bytes may
    // be reordered; the human-authored content may not change).
    let rewritten = std::fs::read_to_string(&path).expect("read blessed copy");
    for kind in before.kinds.values() {
        let escaped = serde_json::to_string(&kind.reason).expect("escape reason");
        assert!(
            rewritten.contains(&escaped),
            "curated reason not found verbatim after bless: {}",
            kind.reason
        );
    }

    // A second bless over the same inputs is byte-identical.
    bless(
        &path,
        &FixtureCheck {
            total_inputs: 4_000_000,
            failures: inputs
                .iter()
                .map(|i| {
                    (
                        i.as_str(),
                        "Parse error at position 3: some error".to_string(),
                    )
                })
                .collect(),
        },
    );
    assert_eq!(
        rewritten,
        std::fs::read_to_string(&path).expect("read second bless"),
        "second bless mutated the snapshot"
    );
}
