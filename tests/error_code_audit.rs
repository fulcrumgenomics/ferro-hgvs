//! Sync check between `src/error_handling/registry.rs` and the audit fixture
//! `tests/fixtures/error_code_audit.json` (companion to
//! `docs/spec/error_code_audit.md`).
//!
//! The audit document maps every error/warning code to the HGVS spec section
//! it enforces. To prevent the document from drifting silently as codes are
//! added, removed, or renamed, this test asserts that:
//!
//!   1. every code in the runtime registry has a matching audit entry, and
//!   2. every audit entry refers to a code that exists in the runtime registry
//!      (no stale rows).
//!
//! Tracking issue: fulcrumgenomics/ferro-hgvs#81 item L1.

use ferro_hgvs::error_handling::list_all_codes;
use serde::Deserialize;
use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;

/// One row from the audit fixture. Only fields that are validated here are
/// modeled; other fields (descriptions, spec URLs, notes) are tolerated as
/// pass-through metadata so the fixture stays human-friendly.
#[derive(Debug, Deserialize)]
struct AuditEntry {
    code: String,
    status: String,
    #[serde(default)]
    spec_section: Option<String>,
}

#[derive(Debug, Deserialize)]
struct AuditFile {
    codes: Vec<AuditEntry>,
}

/// Allowed values of the `status` field. Kept in lockstep with the legend in
/// `docs/spec/error_code_audit.md`.
const VALID_STATUSES: &[&str] = &["enforced", "partial", "registered", "missing", "infra"];

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/error_code_audit.json")
}

/// Source-tree subset where `ErrorType::<Variant>` or
/// `NormalizationWarning::<Variant>` references *should* count as "this
/// variant is emitted as a warning at runtime."
///
/// Excluded files are bookkeeping shims that map every variant by definition
/// (`config.rs`, `python.rs`, `bin/ferro.rs`, `error_handling/types.rs`,
/// `error_handling/codes.rs`, `error_handling/registry.rs`) and the override
/// surfaces (`error_handling/mod.rs`, `normalize/config.rs`) that only set
/// per-code policy without raising the warning. The remaining set captures
/// the real emission path: preprocessor invocations of correction functions
/// that return `ErrorType::*`, post-normalization advisory emissions through
/// `NormalizationWarning::*`, and any future emission site.
const EMISSION_SCAN_PATHS: &[&str] = &[
    "src/error_handling/preprocessor.rs",
    "src/hgvs",
    "src/normalize/mod.rs",
    "src/normalize/overlap.rs",
    "src/normalize/rules.rs",
    "src/normalize/shuffle.rs",
    "src/convert",
    "src/reference",
    "src/spdi",
    "src/vcf",
];

/// Map from a `NormalizationWarning::<Variant>` name to the registry's
/// `name` field (the `ErrorType::*` variant name the audit row identifies
/// the code by). Use when the emission variant and the registry name diverge.
const NORMALIZATION_WARNING_TO_REGISTRY_NAME: &[(&str, &str)] = &[
    // W5001 — both names match, so no mapping entry needed.
    ("OverlapConflict", "OverlapConflictingEdits"),
    // W5003: the warning describes the procedural effect
    // (canonical-split skipped), the registry/error name describes the
    // underlying spec violation (variant exceeds the reference). Same
    // condition, two names.
    ("CanonicalSplitSkipped", "VariantExceedsReference"),
];

/// Scan [`EMISSION_SCAN_PATHS`] under `CARGO_MANIFEST_DIR` for distinct
/// variant names actually used to drive warning emission. Recognizes both
/// `ErrorType::<Variant>` (preprocessor-style emission) and
/// `NormalizationWarning::<Variant>` (post-normalization emission), mapping
/// the latter through [`NORMALIZATION_WARNING_TO_REGISTRY_NAME`] so the
/// returned set keys on the registry's `name` field.
fn emitted_error_type_variants() -> BTreeSet<String> {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut found = BTreeSet::new();
    for rel in EMISSION_SCAN_PATHS {
        scan_for_variants(&manifest.join(rel), &mut found);
    }
    found
}

fn scan_for_variants(path: &std::path::Path, out: &mut BTreeSet<String>) {
    if !path.exists() {
        return;
    }
    if path.is_file() {
        if path.extension().and_then(|e| e.to_str()) == Some("rs") {
            if let Ok(text) = std::fs::read_to_string(path) {
                extract_error_type_variants(&text, out);
                extract_normalization_warning_variants(&text, out);
            }
        }
        return;
    }
    if let Ok(rd) = std::fs::read_dir(path) {
        for entry in rd.flatten() {
            scan_for_variants(&entry.path(), out);
        }
    }
}

fn extract_normalization_warning_variants(text: &str, out: &mut BTreeSet<String>) {
    let needle = "NormalizationWarning::";
    let mut i = 0;
    while let Some(pos) = text[i..].find(needle) {
        let start = i + pos + needle.len();
        let rest = &text[start..];
        let end = rest
            .find(|c: char| !c.is_ascii_alphanumeric() && c != '_')
            .unwrap_or(rest.len());
        if end > 0 {
            let name = &rest[..end];
            if name.chars().next().is_some_and(|c| c.is_ascii_uppercase()) {
                let mapped = NORMALIZATION_WARNING_TO_REGISTRY_NAME
                    .iter()
                    .find_map(|(emit, reg)| (*emit == name).then_some(*reg))
                    .unwrap_or(name);
                out.insert(mapped.to_string());
            }
        }
        i = start + end.max(1);
    }
}

fn extract_error_type_variants(text: &str, out: &mut BTreeSet<String>) {
    // Match `ErrorType::<Variant>` where `<Variant>` starts with an ASCII
    // uppercase letter and continues with alphanumerics or underscore.
    let needle = "ErrorType::";
    let mut i = 0;
    while let Some(pos) = text[i..].find(needle) {
        let start = i + pos + needle.len();
        let rest = &text[start..];
        let end = rest
            .find(|c: char| !c.is_ascii_alphanumeric() && c != '_')
            .unwrap_or(rest.len());
        if end > 0 {
            let name = &rest[..end];
            if name.chars().next().is_some_and(|c| c.is_ascii_uppercase()) {
                out.insert(name.to_string());
            }
        }
        i = start + end.max(1);
    }
}

fn load_fixture() -> AuditFile {
    let path = fixture_path();
    let raw = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read audit fixture at {}: {}", path.display(), e));
    serde_json::from_str(&raw).unwrap_or_else(|e| {
        panic!(
            "audit fixture at {} is not valid JSON matching the expected schema: {}",
            path.display(),
            e
        )
    })
}

#[test]
fn audit_fixture_covers_every_registry_code() {
    let registry_codes: BTreeSet<String> = list_all_codes()
        .iter()
        .map(|c| c.code.to_string())
        .collect();
    let audit = load_fixture();
    let audit_codes: BTreeSet<String> = audit.codes.iter().map(|e| e.code.clone()).collect();

    let missing_in_audit: Vec<&String> = registry_codes.difference(&audit_codes).collect();
    assert!(
        missing_in_audit.is_empty(),
        "registry codes have no entry in tests/fixtures/error_code_audit.json (and therefore no \
         row in docs/spec/error_code_audit.md). Add an audit row for each code listed below, \
         then re-run this test:\n  {:?}",
        missing_in_audit
    );
}

#[test]
fn audit_fixture_has_no_stale_codes() {
    let registry_codes: BTreeSet<String> = list_all_codes()
        .iter()
        .map(|c| c.code.to_string())
        .collect();
    let audit = load_fixture();
    let audit_codes: BTreeSet<String> = audit.codes.iter().map(|e| e.code.clone()).collect();

    let stale_in_audit: Vec<&String> = audit_codes.difference(&registry_codes).collect();
    assert!(
        stale_in_audit.is_empty(),
        "tests/fixtures/error_code_audit.json refers to codes that no longer exist in the \
         registry (src/error_handling/registry.rs). Remove the stale rows from the fixture and \
         from docs/spec/error_code_audit.md:\n  {:?}",
        stale_in_audit
    );
}

#[test]
fn audit_fixture_codes_are_unique() {
    let audit = load_fixture();
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    let mut duplicates: Vec<&str> = Vec::new();
    for entry in &audit.codes {
        if !seen.insert(entry.code.as_str()) {
            duplicates.push(entry.code.as_str());
        }
    }
    assert!(
        duplicates.is_empty(),
        "tests/fixtures/error_code_audit.json contains duplicate code entries: {:?}",
        duplicates
    );
}

#[test]
fn audit_fixture_statuses_are_valid() {
    let audit = load_fixture();
    let mut bad: Vec<(String, String)> = Vec::new();
    for entry in &audit.codes {
        if !VALID_STATUSES.contains(&entry.status.as_str()) {
            bad.push((entry.code.clone(), entry.status.clone()));
        }
    }
    assert!(
        bad.is_empty(),
        "tests/fixtures/error_code_audit.json contains entries with invalid `status` values \
         (allowed: {:?}). Update either the entry or the legend in \
         docs/spec/error_code_audit.md:\n  {:?}",
        VALID_STATUSES,
        bad
    );
}

#[test]
fn audit_fixture_non_infra_entries_cite_a_spec_section() {
    // `infra` rows are allowed to omit `spec_section` (they document I/O,
    // missing reference data, and similar internal failures with no spec
    // mapping). Every other status must cite a spec section so reviewers can
    // verify the mapping without re-reading the registry.
    let audit = load_fixture();
    let mut missing_citation: Vec<String> = Vec::new();
    for entry in &audit.codes {
        if entry.status == "infra" {
            continue;
        }
        let cited = entry
            .spec_section
            .as_ref()
            .map(|s| !s.trim().is_empty())
            .unwrap_or(false);
        if !cited {
            missing_citation.push(entry.code.clone());
        }
    }
    assert!(
        missing_citation.is_empty(),
        "audit entries with status enforced/partial/missing must cite a `spec_section`. \
         Codes lacking a citation:\n  {:?}",
        missing_citation
    );
}

#[test]
fn audit_enforced_warning_rows_have_emission_sites() {
    // For every fixture row marked `enforced` for a W-code (warning), look up
    // the `ErrorType::*` variant name from the registry's `name` field and
    // assert that the variant is referenced from one of the emission-relevant
    // files under `src/`. A failure here means either:
    //
    //   (a) an audit row was reclassified to `enforced` without wiring an
    //       emission site, or
    //   (b) an emission site was deleted but the audit row was not updated.
    //
    // E-codes use a different code path (`ErrorCode::*` raised via
    // `Diagnostic`) and are not covered by this test; they are spot-checked
    // by reviewer judgement and the existing parser/preprocessor unit tests.
    let registry: HashMap<String, &'static str> = list_all_codes()
        .iter()
        .map(|c| (c.code.to_string(), c.name))
        .collect();
    let audit = load_fixture();
    let emitted = emitted_error_type_variants();

    let mut not_emitted: Vec<(String, &'static str)> = Vec::new();
    for entry in &audit.codes {
        if !entry.code.starts_with('W') {
            continue;
        }
        if entry.status != "enforced" {
            continue;
        }
        let Some(name) = registry.get(&entry.code) else {
            // covered by audit_fixture_has_no_stale_codes
            continue;
        };
        if !emitted.contains(*name) {
            not_emitted.push((entry.code.clone(), *name));
        }
    }

    assert!(
        not_emitted.is_empty(),
        "audit rows classified `enforced` for W-codes whose `ErrorType::<Variant>` \
         is not referenced from any emission-relevant file under src/ \
         (preprocessor, hgvs/, normalize/{{rules,shuffle}}.rs, convert/, \
         reference/, spdi/, vcf/). Either wire an emission site, or reclassify \
         the audit row to `registered` (no emission) or `partial` (some \
         emission with a documented gap).\n  {:?}",
        not_emitted
    );
}

#[test]
fn audit_registered_warning_rows_have_no_emission_site() {
    // The companion to `audit_enforced_warning_rows_have_emission_sites`. A
    // row classified `registered` declares "no preprocessor or parser site
    // emits this code." If a future PR wires up emission, the row's status
    // must move to `enforced` (or `partial`) — this test catches the case
    // where someone wires emission but forgets to update the audit.
    let registry: HashMap<String, &'static str> = list_all_codes()
        .iter()
        .map(|c| (c.code.to_string(), c.name))
        .collect();
    let audit = load_fixture();
    let emitted = emitted_error_type_variants();

    let mut should_upgrade: Vec<(String, &'static str)> = Vec::new();
    for entry in &audit.codes {
        if !entry.code.starts_with('W') {
            continue;
        }
        if entry.status != "registered" {
            continue;
        }
        let Some(name) = registry.get(&entry.code) else {
            continue;
        };
        if emitted.contains(*name) {
            should_upgrade.push((entry.code.clone(), *name));
        }
    }

    assert!(
        should_upgrade.is_empty(),
        "audit rows classified `registered` for W-codes that *do* have an \
         `ErrorType::<Variant>` reference in an emission-relevant file under \
         src/. Either the code is now wired (promote to `enforced` / \
         `partial` and remove the `registered` status), or the reference is a \
         test/import — in which case extend EMISSION_SCAN_PATHS' exclusions \
         instead.\n  {:?}",
        should_upgrade
    );
}
