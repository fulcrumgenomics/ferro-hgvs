//! Per-input failure-expectations framework for the bulk parser fixtures.
//!
//! Replaces the coarse `pass_rate > 99%` thresholds with a committed
//! snapshot of every input expected to fail and why. Tracking issue: #174.
//!
//! ## Snapshot shape
//!
//! Each fixture has a sibling JSON file:
//!
//! ```json
//! {
//!   "kinds": {
//!     "<kind_id>": {
//!       "category": "rejected_by_design" | "spec_not_yet_supported" |
//!                   "malformed_input" | "parser_bug" | "needs_triage",
//!       "reason": "human-readable summary",
//!       "tracking": "#NNN"   // optional
//!     }
//!   },
//!   "expected_failures": {
//!     "<input_string>": "<kind_id>"
//!   }
//! }
//! ```
//!
//! ## Test contract
//!
//! Given a fixture iteration that produces `(input, parse_error_or_none)`
//! pairs, the framework asserts:
//!
//! 1. Every input not in `expected_failures` parsed successfully — a fail
//!    for any such input is a regression and the test panics with the
//!    actual error.
//! 2. Every input in `expected_failures` failed to parse — a previously-
//!    failing input that now passes is an "improvement"; the test panics
//!    with instructions to remove the entry from the snapshot.
//! 3. Every kind referenced in `expected_failures` is defined in `kinds`.
//! 4. No expected failure is categorized `parser_bug` or `needs_triage` —
//!    both fail the build. `parser_bug` indicates a tracked real bug;
//!    `needs_triage` is the placeholder assigned by snapshot
//!    regeneration to newly-observed kinds and **must** be moved to a
//!    real category before merge (Phase 2 of #174 enforced this; new
//!    kinds discovered after that point need triage in the same PR).
//!
//! Set `UPDATE_FAILURE_EXPECTATIONS=1` to regenerate the snapshot from the
//! current parse output. Newly-observed kinds are written with category
//! `needs_triage` and the parser's verbatim error string as `reason`;
//! existing kinds' categories are preserved so manual triage is sticky.
//! Any new kinds force a triage decision before the test passes.

use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::OnceLock;

/// Why a given input is allowed to fail.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Category {
    /// We never want to accept this input; rejection is the design.
    RejectedByDesign,
    /// Spec-compliant input; parser support not yet implemented.
    SpecNotYetSupported,
    /// The input itself is malformed (data-quality, not a parser issue).
    MalformedInput,
    /// Real parser bug; should be fixed. Always fails the build.
    ParserBug,
    /// Not yet categorized. Tracked in #174 phase 2. Permitted by the
    /// framework in phase 1 so the framework can land before triage.
    NeedsTriage,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Kind {
    pub category: Category,
    pub reason: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tracking: Option<String>,
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct FailureExpectations {
    pub kinds: BTreeMap<String, Kind>,
    pub expected_failures: BTreeMap<String, String>,
}

/// One fixture's parse outcome handed to [`enforce`].
///
/// `failures[input] = parser_error_string`. Inputs that parsed successfully
/// are simply absent. `total_inputs` is the total number of fixture entries
/// (used only for diagnostic eprintln output).
pub struct FixtureCheck<'a> {
    pub total_inputs: usize,
    pub failures: BTreeMap<&'a str, String>,
}

/// Compare a fixture's parse outcome against the expectations snapshot at
/// `expectations_path`, or regenerate the snapshot if `update_env_var` is
/// set in the environment.
///
/// Panics on regression, improvement, undefined kind, or `parser_bug` hit
/// — with a message that explains how to fix it (regenerate the snapshot,
/// update a category, etc.).
pub fn enforce(expectations_path: &Path, update_env_var: &str, check: FixtureCheck<'_>) {
    if std::env::var_os(update_env_var).is_some() {
        regenerate_snapshot(expectations_path, &check);
        return;
    }

    let snapshot = read_snapshot(expectations_path);
    let mut errors: Vec<String> = Vec::new();

    let mut regressions: Vec<(&str, &str)> = Vec::new();
    let mut parser_bug_hits: Vec<(&str, &str)> = Vec::new();
    let mut needs_triage_hits: Vec<(&str, &str)> = Vec::new();
    let mut undefined_kinds: Vec<(&str, &str)> = Vec::new();

    for (input, err) in &check.failures {
        match snapshot.expected_failures.get(*input) {
            None => regressions.push((input, err.as_str())),
            Some(kind_id) => match snapshot.kinds.get(kind_id) {
                None => undefined_kinds.push((input, kind_id.as_str())),
                Some(kind) => match kind.category {
                    Category::ParserBug => parser_bug_hits.push((input, kind_id.as_str())),
                    Category::NeedsTriage => needs_triage_hits.push((input, kind_id.as_str())),
                    _ => { /* expected failure with permitted category */ }
                },
            },
        }
    }

    let mut improvements: Vec<&str> = Vec::new();
    for input in snapshot.expected_failures.keys() {
        if !check.failures.contains_key(input.as_str()) {
            improvements.push(input.as_str());
        }
    }

    if !regressions.is_empty() {
        let mut msg = format!(
            "{} fixture input(s) regressed (used to parse, now fail). Showing up to 20:\n",
            regressions.len()
        );
        for (input, err) in regressions.iter().take(20) {
            msg.push_str(&format!("  {} -> {}\n", input, err));
        }
        msg.push_str(&format!(
            "\nIf the new failure is intended (parser change rejecting more inputs), \
             regenerate the snapshot:\n  {}=1 cargo nextest run --features dev <test-name>\n",
            update_env_var
        ));
        errors.push(msg);
    }

    if !improvements.is_empty() {
        let mut msg = format!(
            "{} fixture input(s) used to fail but now parse — \"improvement\". \
             Update the snapshot to remove these entries (regenerate):\n  {}=1 cargo nextest run \
             --features dev <test-name>\n\nShowing up to 20:\n",
            improvements.len(),
            update_env_var
        );
        for input in improvements.iter().take(20) {
            msg.push_str(&format!("  {}\n", input));
        }
        errors.push(msg);
    }

    if !undefined_kinds.is_empty() {
        let mut msg = format!(
            "{} expected_failures entry/entries reference an undefined kind. \
             Either define the kind under \"kinds\" or regenerate the snapshot. \
             Showing up to 20:\n",
            undefined_kinds.len()
        );
        for (input, kind_id) in undefined_kinds.iter().take(20) {
            msg.push_str(&format!("  {} -> kind={}\n", input, kind_id));
        }
        errors.push(msg);
    }

    if !parser_bug_hits.is_empty() {
        let mut msg = format!(
            "{} expected failure(s) are categorized `parser_bug`. These are tracked \
             real bugs and the build should not pass while they fail. Either fix the \
             bug or re-categorize the kind.\n\nShowing up to 20:\n",
            parser_bug_hits.len()
        );
        for (input, kind_id) in parser_bug_hits.iter().take(20) {
            msg.push_str(&format!("  {} -> kind={}\n", input, kind_id));
        }
        errors.push(msg);
    }

    if !needs_triage_hits.is_empty() {
        let mut msg = format!(
            "{} expected failure(s) reference a kind still tagged `needs_triage`. \
             Phase 2 of #174 closed: every kind must be assigned a non-failing \
             category (`rejected_by_design` / `spec_not_yet_supported` / \
             `malformed_input`) before merge. (`parser_bug` is also a real \
             category but still fails the build until the bug is fixed.) Edit \
             the snapshot's `kinds.<kind>.category` field directly.\n\n\
             Showing up to 20:\n",
            needs_triage_hits.len()
        );
        for (input, kind_id) in needs_triage_hits.iter().take(20) {
            msg.push_str(&format!("  {} -> kind={}\n", input, kind_id));
        }
        errors.push(msg);
    }

    if !errors.is_empty() {
        let header = format!(
            "Failure-expectations check failed for {} ({} of {} inputs failed):\n\n",
            expectations_path.display(),
            check.failures.len(),
            check.total_inputs
        );
        panic!("{}{}", header, errors.join("\n"));
    }
}

fn read_snapshot(path: &Path) -> FailureExpectations {
    if !path.exists() {
        panic!(
            "Failure-expectations snapshot not found at {}.\n\
             Generate it from current parse output:\n\
             \n  UPDATE_FAILURE_EXPECTATIONS=1 cargo nextest run --features dev <test-name>\n",
            path.display()
        );
    }
    let raw =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {}", path.display(), e));
    serde_json::from_str(&raw).unwrap_or_else(|e| panic!("parse {}: {}", path.display(), e))
}

fn regenerate_snapshot(path: &Path, check: &FixtureCheck<'_>) {
    // Preserve any existing categories/tracking so manual triage is sticky:
    // re-read the current snapshot (if any) and only insert *new* kinds.
    let mut existing: FailureExpectations = if path.exists() {
        let raw = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("read {}: {}", path.display(), e));
        serde_json::from_str(&raw)
            .unwrap_or_else(|e| panic!("parse existing {}: {}", path.display(), e))
    } else {
        FailureExpectations::default()
    };

    // Group failing inputs by *normalized* error string so per-input variability
    // (position numbers and the unparsed-remainder substring inside the nom
    // error) collapses into a small set of kinds. Phase 2 triage will rename
    // these to nicer ids and assign categories.
    let mut next_expected: BTreeMap<String, String> = BTreeMap::new();
    let mut new_kinds_seen: BTreeMap<String, ()> = BTreeMap::new();
    for (input, err) in &check.failures {
        let kind_id = normalize_error_for_kind(err);
        next_expected.insert((*input).to_string(), kind_id.clone());
        new_kinds_seen.insert(kind_id, ());
    }

    // Add any newly-observed kinds with NeedsTriage; preserve existing ones.
    for kind_id in new_kinds_seen.keys() {
        existing
            .kinds
            .entry(kind_id.clone())
            .or_insert_with(|| Kind {
                category: Category::NeedsTriage,
                reason: kind_id.clone(),
                tracking: None,
            });
    }

    // Drop kinds no longer referenced by any expected failure.
    existing
        .kinds
        .retain(|kind_id, _| next_expected.values().any(|v| v == kind_id));

    existing.expected_failures = next_expected;

    let json = serde_json::to_string_pretty(&existing).expect("serialize FailureExpectations");
    std::fs::write(path, json + "\n").unwrap_or_else(|e| panic!("write {}: {}", path.display(), e));
    eprintln!(
        "Updated failure-expectations snapshot: {} ({} expected failures, {} kinds)",
        path.display(),
        existing.expected_failures.len(),
        existing.kinds.len()
    );
}

/// Normalize a parser error string to use as a kind identifier so that
/// failures sharing the same parser-issue type (but differing in position
/// number or unparsed-input remainder) collapse to one kind. Two
/// replacements:
///
/// - `position 17` -> `position <N>`: per-input position drops out.
/// - `input: "<anything>"` -> `input: "<remainder>"`: the nom-error
///   substring containing the unparsed tail of the input drops out.
fn normalize_error_for_kind(err: &str) -> String {
    static POSITION_RE: OnceLock<Regex> = OnceLock::new();
    static INPUT_RE: OnceLock<Regex> = OnceLock::new();
    let position_re = POSITION_RE.get_or_init(|| Regex::new(r"position \d+").unwrap());
    let input_re = INPUT_RE.get_or_init(|| Regex::new(r#"input: "[^"]*""#).unwrap());
    let after_pos = position_re.replace_all(err, "position <N>");
    input_re
        .replace_all(&after_pos, r#"input: "<remainder>""#)
        .into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_collapses_position_and_input_remainder() {
        let a = r#"Parse error at position 10: Failed to parse variant: Error(Error { input: "(1_?)del", code: Digit })"#;
        let b = r#"Parse error at position 12: Failed to parse variant: Error(Error { input: "X_Y_Z", code: Digit })"#;
        assert_eq!(normalize_error_for_kind(a), normalize_error_for_kind(b));
    }

    #[test]
    fn normalize_keeps_distinct_codes_distinct() {
        let digit = r#"Parse error at position 10: Failed to parse variant: Error(Error { input: "x", code: Digit })"#;
        let tag = r#"Parse error at position 10: Failed to parse variant: Error(Error { input: "x", code: Tag })"#;
        assert_ne!(
            normalize_error_for_kind(digit),
            normalize_error_for_kind(tag)
        );
    }
}
