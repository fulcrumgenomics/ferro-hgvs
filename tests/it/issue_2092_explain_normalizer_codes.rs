//! `ferro explain <CODE>` must resolve every code the tool can emit, not just
//! the preprocessor's `W`/`E`-numbered half (#2092).
//!
//! `ferro normalize` prints two code namespaces to stderr in the same
//! `warning[<CODE>]: <message>` shape: the preprocessor's `W`-numbered SVA
//! codes and the normalizer's SCREAMING_SNAKE codes (`src/bin/ferro.rs`
//! documents the split at the emission site). `explain` resolved only the
//! first, so a user who saw `warning[MEMBERS_COALESCED_FROM_REPORTED_FORM]:`
//! and asked about it got a bare "Unknown code" back.
//!
//! # Why the code set is scanned out of `src/normalize/mod.rs`
//!
//! A hardcoded list of normalizer codes here would be blind to the next
//! variant somebody adds — which is the failure mode the issue asks the guard
//! to survive. So [`normalizer_codes_from_source`] reads the `code()` match
//! arms, and [`normalizer_warning_variants_from_source`] independently reads
//! the enum *declarations*. The two must agree
//! (`scanned_code_arms_match_the_declared_variants`), so a variant whose arm
//! the scan cannot parse fails loudly instead of quietly shrinking the set the
//! other tests then check. Without that pairing every test in this file would
//! pass vacuously the day the scan stops matching.

use ferro_hgvs::error_handling::{get_code_info, CodeInfo};
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::process::Command;

/// The file both scans read. The two enums and their `code()` impls live here.
const NORMALIZE_MOD: &str = "src/normalize/mod.rs";

/// The two diagnostic enums the normalizer reports through
/// `NormalizeResult` — warnings and info-grade signals. Each owns a
/// `pub fn code(&self) -> &'static str`.
const DIAGNOSTIC_ENUMS: &[&str] = &["NormalizationWarning", "NormalizationInfo"];

fn normalize_mod_source() -> String {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(NORMALIZE_MOD);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("failed to read {}: {e}", path.display()))
}

/// Every `Self::<Variant> { .. } => "<CODE>",` arm in the file, as
/// `variant -> code`.
///
/// That shape occurs only inside the two `code()` impls: the `Display` impls
/// match the same variants but bind fields and end in `write!`, so they cannot
/// collide with this pattern.
fn normalizer_codes_from_source() -> BTreeMap<String, String> {
    let mut found = BTreeMap::new();
    for line in normalize_mod_source().lines() {
        let line = line.trim();
        let Some(rest) = line.strip_prefix("Self::") else {
            continue;
        };
        let Some((variant, tail)) = rest.split_once(" { .. } => \"") else {
            continue;
        };
        let Some(code) = tail.strip_suffix("\",") else {
            continue;
        };
        if !variant.chars().all(|c| c.is_ascii_alphanumeric()) {
            continue;
        }
        if code.is_empty() || !code.chars().all(|c| c.is_ascii_uppercase() || c == '_') {
            continue;
        }
        if let Some(prev) = found.insert(variant.to_string(), code.to_string()) {
            panic!("variant {variant} has two code() arms: {prev} and {code}");
        }
    }
    found
}

/// Every variant declared by the enums in [`DIAGNOSTIC_ENUMS`], read from the
/// declaration site rather than from `code()`.
fn normalizer_warning_variants_from_source() -> BTreeSet<String> {
    let source = normalize_mod_source();
    let mut found = BTreeSet::new();
    for enum_name in DIAGNOSTIC_ENUMS {
        let header = format!("pub enum {enum_name} {{");
        let start = source
            .find(&header)
            .unwrap_or_else(|| panic!("`{header}` not found in {NORMALIZE_MOD}"));
        let body = &source[start + header.len()..];
        // Both enums are top-level items, so the first line that is exactly a
        // closing brace at column 0 ends the declaration.
        let end = body
            .find("\n}")
            .unwrap_or_else(|| panic!("unterminated `{enum_name}` declaration"));
        for line in body[..end].lines() {
            let Some(variant) = line.strip_prefix("    ").and_then(|l| l.strip_suffix(" {")) else {
                continue;
            };
            if variant.starts_with(|c: char| c.is_ascii_uppercase())
                && variant.chars().all(|c| c.is_ascii_alphanumeric())
            {
                assert!(
                    found.insert(variant.to_string()),
                    "variant {variant} declared twice"
                );
            }
        }
    }
    found
}

/// The scan is only as good as its agreement with the declarations. If a new
/// variant is spelled in a way the arm pattern misses, this fails here rather
/// than silently narrowing what every other test in the file checks.
#[test]
fn scanned_code_arms_match_the_declared_variants() {
    let by_arm: BTreeSet<String> = normalizer_codes_from_source().into_keys().collect();
    let declared = normalizer_warning_variants_from_source();

    assert!(
        !declared.is_empty(),
        "declaration scan found no variants in {NORMALIZE_MOD} — the scan is broken, not the code"
    );
    assert_eq!(
        by_arm,
        declared,
        "the `code()` arms scanned out of {NORMALIZE_MOD} do not match the variants declared \
         there.\n  arms but not declared: {:?}\n  declared but no arm scanned: {:?}",
        by_arm.difference(&declared).collect::<Vec<_>>(),
        declared.difference(&by_arm).collect::<Vec<_>>(),
    );
}

/// The issue's headline defect: the codes `ferro normalize` prints are not the
/// codes `ferro explain` knows.
#[test]
fn every_normalizer_code_resolves_through_get_code_info() {
    let codes = normalizer_codes_from_source();
    assert!(
        !codes.is_empty(),
        "code() scan found nothing in {NORMALIZE_MOD} — the scan is broken, not the code"
    );

    let unresolved: Vec<String> = codes
        .iter()
        .filter(|(_, code)| get_code_info(code).is_none())
        .map(|(variant, code)| format!("{code} (NormalizationWarning::{variant})"))
        .collect();

    assert!(
        unresolved.is_empty(),
        "{} of {} normalizer codes do not resolve via `ferro explain`:\n  {}",
        unresolved.len(),
        codes.len(),
        unresolved.join("\n  "),
    );
}

/// Resolving is not enough: an entry that resolves to empty prose is the same
/// dead end wearing a different hat.
///
/// An unresolvable code is reported here too, deliberately: skipping it would
/// make this test pass vacuously on exactly the state #2092 describes, where
/// nothing resolves at all.
#[test]
fn every_normalizer_code_resolves_to_substantive_content() {
    let mut thin: Vec<String> = Vec::new();
    for (variant, code) in normalizer_codes_from_source() {
        match get_code_info(&code) {
            None => thin.push(format!("{code} ({variant}): resolves to nothing")),
            Some(info) if info.summary.trim().len() < 20 || info.explanation.trim().len() < 80 => {
                thin.push(format!(
                    "{code} ({variant}): summary {} chars, explanation {} chars",
                    info.summary.trim().len(),
                    info.explanation.trim().len(),
                ));
            }
            Some(_) => {}
        }
    }
    assert!(
        thin.is_empty(),
        "normalizer codes with no real explanatory content behind them:\n  {}",
        thin.join("\n  "),
    );
}

/// The three codes the issue names, pinned to the exact entry each must reach.
/// Asserting the entry's `name` (not merely that *something* came back) is what
/// stops a wrong-but-populated alias from passing.
#[test]
fn the_codes_named_in_issue_2092_resolve_to_their_own_entries() {
    let expected: &[(&str, &str)] = &[
        (
            "MEMBERS_COALESCED_FROM_REPORTED_FORM",
            "MembersCoalescedFromReportedForm",
        ),
        ("INSERTED_SEQUENCE_EXPANDED", "InsertedSequenceExpanded"),
        ("REFSEQ_MISMATCH", "RefSeqMismatch"),
    ];
    for (code, name) in expected {
        let info: &CodeInfo = get_code_info(code)
            .unwrap_or_else(|| panic!("`ferro explain {code}` resolves nothing"));
        assert_eq!(
            info.name, *name,
            "{code} resolved to the wrong registry entry ({} / {})",
            info.code, info.name
        );
    }
}

/// Lowercase input already works for `W`-codes (`get_code_info` uppercases), so
/// it must work for the normalizer namespace too.
#[test]
fn normalizer_codes_resolve_case_insensitively() {
    let upper = get_code_info("REFSEQ_MISMATCH").expect("REFSEQ_MISMATCH must resolve");
    let lower = get_code_info("refseq_mismatch").expect("refseq_mismatch must resolve");
    assert_eq!(upper.code, lower.code);
}

/// `W5005` is minted by `ErrorType::MembersCoalescedFromReportedForm.code()`
/// but was absent from the registry, so `ferro explain W5005` failed for a code
/// the tool documents as its own. Sibling gap to the headline defect, same
/// root: the registry did not cover everything that can reach a user.
#[test]
fn w5005_resolves_under_its_minted_w_code() {
    let info = get_code_info("W5005").expect("W5005 must be registered");
    assert_eq!(info.name, "MembersCoalescedFromReportedForm");
}

fn explain(code: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
        .args(["explain", code])
        .output()
        .expect("failed to run `ferro explain`")
}

/// End-to-end through the binary: the exact command the issue says returns
/// nothing useful.
#[test]
fn explain_cli_resolves_a_normalizer_code() {
    let out = explain("MEMBERS_COALESCED_FROM_REPORTED_FORM");
    assert!(
        out.status.success(),
        "`ferro explain MEMBERS_COALESCED_FROM_REPORTED_FORM` exited {:?}: {}",
        out.status.code(),
        String::from_utf8_lossy(&out.stderr),
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("MembersCoalescedFromReportedForm"),
        "explain output does not name the code's registry entry:\n{stdout}"
    );
}

/// The defect's own shape, guarded against: a code that cannot be resolved must
/// say so in a way the user can tell apart from a code that resolved. Both
/// halves are asserted in one test so neither can be satisfied by the other's
/// behaviour.
#[test]
fn explain_cli_distinguishes_an_unresolvable_code_from_a_resolved_one() {
    let resolved = explain("REFSEQ_MISMATCH");
    let unknown = explain("NOT_A_FERRO_CODE_AT_ALL");

    assert!(resolved.status.success(), "a registered code must succeed");
    assert!(
        !unknown.status.success(),
        "an unregistered code must not exit 0 — a silent success is exactly the failure #2092 \
         reports"
    );

    let unknown_msg = format!(
        "{}{}",
        String::from_utf8_lossy(&unknown.stdout),
        String::from_utf8_lossy(&unknown.stderr)
    );
    assert!(
        unknown_msg.contains("NOT_A_FERRO_CODE_AT_ALL"),
        "the failure must echo the code the user typed:\n{unknown_msg}"
    );
    // The user has just been handed a SCREAMING_SNAKE code by `ferro normalize`
    // and told it is unknown. Naming both namespaces is what turns that from a
    // dead end into a next step.
    assert!(
        unknown_msg.contains("REFSEQ_MISMATCH"),
        "the failure for a normalizer-shaped code must name a real normalizer code so the user \
         can tell 'no such code' from 'this namespace is not covered':\n{unknown_msg}"
    );
}

/// `explain --list` must enumerate the normalizer namespace too — a listing
/// that omits half the codes reproduces the defect for anyone who goes looking
/// rather than asking about a specific code.
#[test]
fn explain_list_enumerates_the_normalizer_namespace() {
    let out = Command::new(env!("CARGO_BIN_EXE_ferro"))
        .args(["explain", "--list"])
        .output()
        .expect("failed to run `ferro explain --list`");
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);

    let missing: Vec<String> = normalizer_codes_from_source()
        .into_values()
        .filter(|code| !stdout.contains(code.as_str()))
        .collect();
    assert!(
        missing.is_empty(),
        "`ferro explain --list` does not list these normalizer codes: {missing:?}"
    );
}
