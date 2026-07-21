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
//! ## Re-blessing (`UPDATE_FAILURE_EXPECTATIONS=1`)
//!
//! Set `UPDATE_FAILURE_EXPECTATIONS=1` to update the snapshot from the
//! current parse output. The snapshots are hand-curated — the kind ids are
//! human-chosen names, and each `reason` is prose citing the HGVS spec.
//! None of that is recoverable from a parser error string, so the blesser
//! is a **merge**, never a regeneration:
//!
//! - An input that is already in `expected_failures` and still fails keeps
//!   its curated kind id verbatim. Its kind's `category`, `reason`, and
//!   `tracking` are never rewritten.
//! - An input that fails and is *not* in `expected_failures` is added under
//!   a synthetic kind id prefixed [`AUTO_KIND_PREFIX`] (derived from the
//!   normalized error string) with category `needs_triage`. The prefix
//!   makes machine-written ids obvious next to curated ones, and
//!   `needs_triage` fails the build until a human names and categorizes it.
//! - An input that no longer fails is removed from `expected_failures` —
//!   leaving it would permanently fail the "improvement" check — but its
//!   *kind* is kept, so the curated prose survives, and the removal is
//!   printed in the bless report.
//! - Kinds are **never deleted**, even when nothing references them any
//!   more. Unreferenced kinds are reported so a human can delete them
//!   deliberately; the blesser will not throw prose away on its own.
//! - The blesser refuses to run (panics) when it cannot guarantee the
//!   above: an unparseable snapshot, a snapshot carrying fields this
//!   framework does not model (`deny_unknown_fields`), or a fixture run
//!   that produced zero inputs (dataset missing/skipped — blessing would
//!   wipe the file).
//!
//! Merge behaviour is covered directly by
//! `tests/it/failure_expectations_blesser_tests.rs`.

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
#[serde(deny_unknown_fields)]
pub struct Kind {
    pub category: Category,
    pub reason: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tracking: Option<String>,
}

#[derive(Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FailureExpectations {
    pub kinds: BTreeMap<String, Kind>,
    pub expected_failures: BTreeMap<String, String>,
}

/// One fixture's parse outcome handed to [`enforce`].
///
/// `failures[input] = parser_error_string`. Inputs that parsed successfully
/// are simply absent. `total_inputs` is the total number of fixture entries;
/// besides diagnostic output it drives the zero-inputs safety guard in
/// [`bless`] (a run of 0 inputs must not empty a curated snapshot), so it must
/// reflect the real corpus size, not just the failing subset.
pub struct FixtureCheck<'a> {
    pub total_inputs: usize,
    pub failures: BTreeMap<&'a str, String>,
}

/// Compare a fixture's parse outcome against the expectations snapshot at
/// `expectations_path`, or bless (non-destructively merge into) the snapshot if
/// `update_env_var` is set in the environment.
///
/// Panics on regression, improvement, undefined kind, or `parser_bug` hit
/// — with a message that explains how to fix it (regenerate the snapshot,
/// update a category, etc.).
pub fn enforce(expectations_path: &Path, update_env_var: &str, check: FixtureCheck<'_>) {
    if std::env::var_os(update_env_var).is_some() {
        bless(expectations_path, &check);
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
             bless the snapshot (a non-destructive merge that keeps curated kinds):\n  \
             {}=1 cargo nextest run --features dev <test-name>\n",
            update_env_var
        ));
        errors.push(msg);
    }

    if !improvements.is_empty() {
        let mut msg = format!(
            "{} fixture input(s) used to fail but now parse — \"improvement\". \
             Bless the snapshot to drop these entries (a merge; curated kinds are kept):\n  \
             {}=1 cargo nextest run --features dev <test-name>\n\nShowing up to 20:\n",
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
             Either define the kind under \"kinds\" or bless the snapshot (which writes a \
             needs_triage placeholder). Showing up to 20:\n",
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

/// Prefix on kind ids the blesser invents for newly-observed failures, so
/// machine-written ids are never mistaken for hand-curated ones.
pub const AUTO_KIND_PREFIX: &str = "auto:";

/// What a bless changed, returned so callers (and tests) can inspect it
/// instead of scraping stderr.
#[derive(Debug, Default)]
pub struct BlessReport {
    /// Failing inputs whose curated kind mapping was carried over as-is.
    pub preserved_inputs: usize,
    /// Newly-failing inputs added to `expected_failures`.
    pub new_inputs: Vec<String>,
    /// Kind ids invented for those new inputs (all `needs_triage`).
    pub new_kinds: Vec<String>,
    /// Inputs that no longer fail; removed from `expected_failures`.
    pub resolved_inputs: Vec<String>,
    /// Kinds no longer referenced by any expected failure. Kept on disk.
    pub orphaned_kinds: Vec<String>,
    /// Kinds referenced by an expected failure but never defined; given a
    /// `needs_triage` placeholder so the enforcing test can explain itself.
    pub synthesized_kinds: Vec<String>,
}

/// Merge the current parse outcome into the snapshot at `path`, preserving
/// every piece of human-authored content (see the module docs for the exact
/// semantics), write it back, and print a report of what changed.
///
/// Panics rather than proceeding whenever curated content would be at risk:
/// an unreadable/unparseable snapshot, a snapshot with fields this framework
/// does not model, or a fixture run that produced no inputs at all.
pub fn bless(path: &Path, check: &FixtureCheck<'_>) -> BlessReport {
    let mut snapshot: FailureExpectations = if path.exists() {
        let raw = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("read {}: {}", path.display(), e));
        serde_json::from_str(&raw).unwrap_or_else(|e| {
            panic!(
                "refusing to bless {}: the existing snapshot could not be parsed \
                 ({e}). Rewriting it would destroy hand-curated kinds and reasons; \
                 fix the file (or the framework's schema) first.",
                path.display()
            )
        })
    } else {
        FailureExpectations::default()
    };

    // A fixture that yielded zero inputs means the dataset was missing or
    // skipped, not that everything now parses. Blessing would empty the file.
    if check.total_inputs == 0 && !snapshot.expected_failures.is_empty() {
        panic!(
            "refusing to bless {}: the fixture produced 0 inputs, so this run \
             carries no evidence about the {} curated expected failure(s). \
             Check that the fixture data is present (Git LFS) and re-run.",
            path.display(),
            snapshot.expected_failures.len()
        );
    }

    let report = merge_into(&mut snapshot, check);

    let json = serde_json::to_string_pretty(&snapshot).expect("serialize FailureExpectations");
    std::fs::write(path, json + "\n").unwrap_or_else(|e| panic!("write {}: {}", path.display(), e));

    print_bless_report(path, &snapshot, &report);
    report
}

/// Pure merge step of [`bless`]: updates `snapshot` in place and describes
/// what changed. Kept separate from all I/O so it is directly testable.
/// Merge this run's failures into `snapshot`, returning a [`BlessReport`].
///
/// Contract (what "merge" preserves): the hand-curated **kinds** (their
/// category / reason / tracking prose) are never deleted — that is the whole
/// point over the old destructive regenerate. `expected_failures`, by contrast,
/// is rebuilt from the current run: a curated mapping whose input still fails is
/// carried over verbatim, a newly-failing input gets an `auto:` needs-triage
/// kind, and an input **absent from this run is treated as resolved and its
/// mapping dropped** (its kind is kept). This assumes the run covers the full
/// corpus; a wholly-missing dataset is caught by [`bless`]'s zero-inputs guard,
/// but a *partial* run would silently drop the absent inputs' mappings.
fn merge_into(snapshot: &mut FailureExpectations, check: &FixtureCheck<'_>) -> BlessReport {
    let mut report = BlessReport::default();
    let mut next_expected: BTreeMap<String, String> = BTreeMap::new();

    for (input, err) in &check.failures {
        match snapshot.expected_failures.get(*input) {
            // Curated (or previously blessed) mapping: carry it over verbatim.
            Some(kind_id) => {
                next_expected.insert((*input).to_string(), kind_id.clone());
                report.preserved_inputs += 1;
            }
            // Genuinely new failure: invent an obviously-machine-written kind
            // id from the normalized error, awaiting human triage.
            None => {
                let normalized = normalize_error_for_kind(err);
                let kind_id = format!("{}{}", AUTO_KIND_PREFIX, normalized);
                if !snapshot.kinds.contains_key(&kind_id) {
                    snapshot.kinds.insert(
                        kind_id.clone(),
                        Kind {
                            category: Category::NeedsTriage,
                            reason: normalized,
                            tracking: None,
                        },
                    );
                    report.new_kinds.push(kind_id.clone());
                }
                next_expected.insert((*input).to_string(), kind_id);
                report.new_inputs.push((*input).to_string());
            }
        }
    }

    // Inputs that no longer fail: drop the mapping (keeping it would fail the
    // "improvement" check forever) but keep the kind, so the curated reason
    // survives and the change is visible in the report.
    for input in snapshot.expected_failures.keys() {
        if !next_expected.contains_key(input) {
            report.resolved_inputs.push(input.clone());
        }
    }
    snapshot.expected_failures = next_expected;

    // Dangling kind references (a typo, or a hand-edit that removed a kind):
    // give them a placeholder so the enforcing test reports `needs_triage`
    // rather than a bare "undefined kind".
    let dangling: Vec<String> = snapshot
        .expected_failures
        .values()
        .filter(|kind_id| !snapshot.kinds.contains_key(*kind_id))
        .cloned()
        .collect();
    for kind_id in dangling {
        if snapshot.kinds.contains_key(&kind_id) {
            continue;
        }
        snapshot.kinds.insert(
            kind_id.clone(),
            Kind {
                category: Category::NeedsTriage,
                reason: format!(
                    "Referenced by expected_failures but never defined; needs a category and reason ({kind_id})."
                ),
                tracking: None,
            },
        );
        report.synthesized_kinds.push(kind_id);
    }

    // Report — never delete — kinds nothing references any more.
    report.orphaned_kinds = snapshot
        .kinds
        .keys()
        .filter(|kind_id| !snapshot.expected_failures.values().any(|v| v == *kind_id))
        .cloned()
        .collect();

    report
}

/// Print a human-readable summary of a bless to stderr.
fn print_bless_report(path: &Path, snapshot: &FailureExpectations, report: &BlessReport) {
    eprintln!(
        "Blessed failure-expectations snapshot: {} ({} expected failures, {} kinds; \
         {} mapping(s) preserved)",
        path.display(),
        snapshot.expected_failures.len(),
        snapshot.kinds.len(),
        report.preserved_inputs
    );
    let sections: [(&str, &Vec<String>); 5] = [
        (
            "NEW failure(s) added (kind id prefixed `auto:`, category `needs_triage` — \
             rename and categorize before merge)",
            &report.new_inputs,
        ),
        (
            "NEW kind id(s) synthesized for those failures (the `auto:`-prefixed ids to \
             rename and give a real category/tracking before merge)",
            &report.new_kinds,
        ),
        (
            "input(s) NO LONGER failing — removed from expected_failures (their kind \
             definitions were kept)",
            &report.resolved_inputs,
        ),
        (
            "kind(s) now UNREFERENCED — kept so the curated reason is not lost; delete \
             by hand if truly obsolete",
            &report.orphaned_kinds,
        ),
        (
            "kind(s) referenced but UNDEFINED — a `needs_triage` placeholder was written",
            &report.synthesized_kinds,
        ),
    ];
    for (label, items) in sections {
        if items.is_empty() {
            continue;
        }
        eprintln!("  {} {}:", items.len(), label);
        for item in items.iter().take(20) {
            eprintln!("    {}", item);
        }
        if items.len() > 20 {
            eprintln!("    ... and {} more", items.len() - 20);
        }
    }
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
