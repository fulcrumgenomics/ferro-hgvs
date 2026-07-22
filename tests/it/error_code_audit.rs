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
    // The deprecated-form correction functions (`correct_deprecated_protein_forms`
    // etc.) construct the `ErrorType::*` values the preprocessor emits; this is
    // their real production emission site (#1123).
    "src/error_handling/corrections.rs",
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
    // W4005: the warning is generic over telomere/centromere markers; the
    // registry/error name is specific to the only unresolvable one (cen).
    ("UnresolvableSpecialPosition", "UnresolvableCentromere"),
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

/// Remove `#[cfg(test)] mod … { … }` blocks from Rust source so `ErrorType::*`
/// references inside unit-test modules are not counted as production emission
/// sites (#1123). Only `mod` blocks are stripped; a non-`mod` `#[cfg(test)]`
/// item (e.g. `#[cfg(test)] use …;`) is left untouched (there are none in the
/// scanned tree, and skipping them avoids over-stripping).
///
/// The closing brace is found by a brace-matcher that skips over string, char,
/// raw-string, and comment content, so braces inside literals/comments (e.g.
/// format strings like `"{x}"`, or `r#"…{…"#`) do not unbalance the count.
fn strip_cfg_test_modules(text: &str) -> String {
    const ATTR: &str = "#[cfg(test)]";
    let bytes = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    let mut copied = 0usize; // next byte index not yet copied to `out`
    let mut search = 0usize; // where to look for the next attribute
    while let Some(rel) = text[search..].find(ATTR) {
        let attr = search + rel;
        let after = attr + ATTR.len();
        // The attribute must be the first non-whitespace on its line. A real
        // `#[cfg(test)]` is emitted at line start (rustfmt), whereas a mention
        // embedded in a doc comment (`/// … #[cfg(test)] …`) or a string is
        // preceded by other text on the line. Without this guard the raw-text
        // search below could read such a mention as a module attribute and,
        // finding `mod ` before the next `;`, strip real production code after
        // it (a false FAIL). The brace-matcher is already literal/comment-aware;
        // this guards the search that *locates* the block.
        let line_start = text[..attr].rfind('\n').map_or(0, |n| n + 1);
        if text[line_start..attr]
            .bytes()
            .any(|c| c != b' ' && c != b'\t')
        {
            search = after;
            continue;
        }
        // Find the first `{` and first `;` after the attribute. A `{` reached
        // before any `;`, with `mod ` in between, marks a test-module block.
        let brace = text[after..].find('{').map(|b| after + b);
        let semi = text[after..].find(';').map(|s| after + s);
        let is_mod_block = matches!((brace, semi),
            (Some(b), s) if s.is_none_or(|s| b < s) && text[after..b].contains("mod "));
        if !is_mod_block {
            search = after;
            continue;
        }
        let open = brace.unwrap();
        let end = match_block_end(bytes, open);
        out.push_str(&text[copied..attr]);
        copied = end;
        search = end;
    }
    out.push_str(&text[copied..]);
    out
}

/// Given `bytes[open] == b'{'`, return the index just past its matching `}`,
/// skipping string / char / raw-string literals and `//` / `/* */` comments so
/// braces inside them are not counted. Best-effort lexing sufficient for the
/// well-formed Rust sources scanned here.
fn match_block_end(bytes: &[u8], open: usize) -> usize {
    let mut depth = 0usize;
    let mut i = open;
    while i < bytes.len() {
        let b = bytes[i];
        match b {
            b'{' => depth += 1,
            b'}' => {
                depth -= 1;
                if depth == 0 {
                    return i + 1;
                }
            }
            b'/' if i + 1 < bytes.len() && bytes[i + 1] == b'/' => {
                i += 2;
                while i < bytes.len() && bytes[i] != b'\n' {
                    i += 1;
                }
                continue;
            }
            b'/' if i + 1 < bytes.len() && bytes[i + 1] == b'*' => {
                // Rust block comments nest (`/* /* */ */`). Track depth so an
                // inner `/* … */` does not end the outer comment early and leak
                // its (possibly brace-bearing) tail into the brace count.
                let mut comment_depth = 1usize;
                i += 2;
                while i + 1 < bytes.len() && comment_depth > 0 {
                    if bytes[i] == b'/' && bytes[i + 1] == b'*' {
                        comment_depth += 1;
                        i += 2;
                    } else if bytes[i] == b'*' && bytes[i + 1] == b'/' {
                        comment_depth -= 1;
                        i += 2;
                    } else {
                        i += 1;
                    }
                }
                continue;
            }
            b'"' => {
                i += 1;
                while i < bytes.len() {
                    match bytes[i] {
                        b'\\' => i += 2,
                        b'"' => {
                            i += 1;
                            break;
                        }
                        _ => i += 1,
                    }
                }
                continue;
            }
            b'\'' => {
                // Char literal (with escapes) or a lifetime (`'a`, no closer).
                // Scan a short window for a closing quote; if none, advance one.
                // 6 bytes covers every char literal up to `'\xNN'`; the only
                // longer form, a `'\u{…}'` escape, carries a *balanced* `{…}`
                // pair, so overshooting it cannot unbalance the brace count.
                let mut j = i + 1;
                let mut closed = false;
                let mut steps = 0;
                while j < bytes.len() && steps < 6 {
                    match bytes[j] {
                        b'\\' => {
                            j += 2;
                            steps += 2;
                        }
                        b'\'' => {
                            closed = true;
                            j += 1;
                            break;
                        }
                        _ => {
                            j += 1;
                            steps += 1;
                        }
                    }
                }
                i = if closed { j } else { i + 1 };
                continue;
            }
            b'r' | b'b' => {
                // Possible raw string: r"…", r#…"…"…#, br"…", br#…"…"…#. Only
                // when `r`/`b` does not continue an identifier.
                let prev_ident =
                    i > 0 && (bytes[i - 1].is_ascii_alphanumeric() || bytes[i - 1] == b'_');
                if !prev_ident {
                    if let Some(next) = raw_string_end(bytes, i) {
                        i = next;
                        continue;
                    }
                }
            }
            _ => {}
        }
        i += 1;
    }
    bytes.len()
}

/// If a raw string literal starts at `start` (`r`/`br` then `#`* then `"`),
/// return the index just past its closing `"#*`; otherwise `None`.
fn raw_string_end(bytes: &[u8], start: usize) -> Option<usize> {
    let mut j = start;
    if bytes[j] == b'b' {
        j += 1;
    }
    if j >= bytes.len() || bytes[j] != b'r' {
        return None;
    }
    j += 1;
    let hashes = j;
    while j < bytes.len() && bytes[j] == b'#' {
        j += 1;
    }
    if j >= bytes.len() || bytes[j] != b'"' {
        return None;
    }
    let n = j - hashes;
    j += 1;
    while j < bytes.len() {
        if bytes[j] == b'"' {
            let mut k = j + 1;
            let mut cnt = 0;
            while k < bytes.len() && cnt < n && bytes[k] == b'#' {
                k += 1;
                cnt += 1;
            }
            if cnt == n {
                return Some(k);
            }
        }
        j += 1;
    }
    Some(bytes.len())
}

fn scan_for_variants(path: &std::path::Path, out: &mut BTreeSet<String>) {
    if !path.exists() {
        return;
    }
    if path.is_file() {
        if path.extension().and_then(|e| e.to_str()) == Some("rs") {
            if let Ok(text) = std::fs::read_to_string(path) {
                // Strip `#[cfg(test)]` modules so `ErrorType::*` references
                // inside unit tests do not count as production emission sites
                // (#1123). Otherwise a test-only reference masks a deleted
                // production emitter and the emission-site audit passes falsely.
                let production = strip_cfg_test_modules(&text);
                extract_error_type_variants(&production, out);
                extract_normalization_warning_variants(&production, out);
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

/// A variant referenced only inside a `#[cfg(test)]` module must NOT count as a
/// production emission site (#1123). Otherwise a test-only reference can mask a
/// deleted production emitter, defeating `audit_enforced_warning_rows_have_emission_sites`.
#[test]
fn emission_scan_excludes_cfg_test_module_refs() {
    let src = r#"
        fn emit() {
            let _ = ErrorType::RealEmitter;
        }

        #[cfg(test)]
        mod tests {
            use super::*;

            #[test]
            fn t() {
                let _ = ErrorType::TestOnly;
            }
        }
    "#;
    let mut found = BTreeSet::new();
    extract_error_type_variants(&strip_cfg_test_modules(src), &mut found);
    assert!(
        found.contains("RealEmitter"),
        "production ErrorType reference must be kept: {found:?}",
    );
    assert!(
        !found.contains("TestOnly"),
        "ErrorType reference inside a #[cfg(test)] module must be excluded: {found:?}",
    );
}

/// The `#[cfg(test)]` stripper must skip braces inside string/raw-string/
/// comment content so an unbalanced-looking brace in a literal does not end the
/// block early and leak a test-only reference (#1123).
#[test]
fn strip_handles_braces_in_strings_raw_strings_and_comments() {
    let src = r####"
        fn emit() { let _ = ErrorType::Real; }

        #[cfg(test)]
        mod tests {
            fn t() {
                let _a = "unbalanced { brace only in string";
                let _b = r#"raw string with } and { braces"#;
                // a line comment with a } brace
                let _ = ErrorType::TestOnly;
            }
        }
    "####;
    let mut found = BTreeSet::new();
    extract_error_type_variants(&strip_cfg_test_modules(src), &mut found);
    assert!(found.contains("Real"), "production ref kept: {found:?}");
    assert!(
        !found.contains("TestOnly"),
        "test-module ref excluded despite braces in literals/comments: {found:?}",
    );
}

/// The stripper must end a `#[cfg(test)]` block at exactly its closing brace, so
/// production code AFTER the module still counts. Placing a production reference
/// after the test module pins the block-end index — a matcher that over-consumes
/// (e.g. a literal/comment brace miscount) would drop `AfterMod` and produce a
/// spurious audit failure (#1123).
#[test]
fn strip_keeps_production_ref_after_cfg_test_module() {
    let src = r####"
        fn emit() { let _ = ErrorType::Before; }

        #[cfg(test)]
        mod tests {
            fn t() {
                let _a = "trailing brace in string }";
                let _b = r#"and a raw one }"#;
                /* block comment with a } brace */
                let _ = ErrorType::TestOnly;
            }
        }

        fn emit_more() { let _ = ErrorType::AfterMod; }
    "####;
    let mut found = BTreeSet::new();
    extract_error_type_variants(&strip_cfg_test_modules(src), &mut found);
    assert!(
        found.contains("Before"),
        "ref before the module kept: {found:?}"
    );
    assert!(
        found.contains("AfterMod"),
        "ref AFTER the module must survive — the block must end at its own brace: {found:?}",
    );
    assert!(
        !found.contains("TestOnly"),
        "test-module ref excluded: {found:?}"
    );
}

/// A `#[cfg(test)]` mention inside a doc comment or string is not a module
/// attribute. The raw-text search must not treat it as one and strip the
/// production item that follows (#1123). Both mentions here are load-bearing:
/// each is followed by `mod ` and a `{` before the next `;`, so without the
/// line-start guard the raw-text search would classify it as a module block and
/// strip the real code after it.
#[test]
fn strip_ignores_cfg_test_mention_in_non_attribute_position() {
    let src = r#"
        /// Mirrors the pattern used by #[cfg(test)] mod tests { below.
        fn documented() { let _ = ErrorType::Documented; }

        fn note() {
            let _msg = "guard #[cfg(test)] mod names { in a message";
            let _ = ErrorType::Noted;
        }
    "#;
    let mut found = BTreeSet::new();
    extract_error_type_variants(&strip_cfg_test_modules(src), &mut found);
    assert!(
        found.contains("Documented") && found.contains("Noted"),
        "a #[cfg(test)] mention in a comment/string must not strip real code: {found:?}",
    );
}

/// Rust block comments nest; the brace-matcher must too. An inner `/* … */`
/// whose tail carries an unbalanced brace must not end the outer comment early
/// and corrupt the depth count, which could leak a test ref (false PASS) or drop
/// production code after the module (false FAIL) (#1123).
#[test]
fn strip_handles_nested_block_comments() {
    let src = r####"
        fn emit() { let _ = ErrorType::Real; }

        #[cfg(test)]
        mod tests {
            fn t() {
                /* outer /* inner } */ still in outer { */
                let _ = ErrorType::TestOnly;
            }
        }

        fn after() { let _ = ErrorType::AfterMod; }
    "####;
    let mut found = BTreeSet::new();
    extract_error_type_variants(&strip_cfg_test_modules(src), &mut found);
    assert!(
        found.contains("Real"),
        "production ref before kept: {found:?}"
    );
    assert!(
        found.contains("AfterMod"),
        "production ref after module kept: {found:?}"
    );
    assert!(
        !found.contains("TestOnly"),
        "test-module ref excluded despite nested block comments: {found:?}",
    );
}
