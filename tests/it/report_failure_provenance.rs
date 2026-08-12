//! The failure reporter must not assert that its callers cannot fail on a pull
//! request.
//!
//! `report-failure.yml` appends a provenance sentence to every tracking issue it
//! opens. That sentence used to be the fixed text *"This workflow does not run on
//! pull requests, so this failure could not have been caught before merge"* — a
//! claim about the CALLER, hardcoded in the CALLEE, and therefore unable to
//! notice when a caller's triggers change.
//!
//! It went stale immediately. #1225 both created this reporter (because
//! `Coverage` stayed red for 15 commits while every PR showed a green board) and
//! gave `coverage.yml` a `pull_request` trigger. From that commit the sentence
//! was false for `Coverage` and true for the other five callers, and five of the
//! recurrences filed on the rolling coverage issue carried it.
//!
//! **The direction of the error is what makes it worth a guard.** It excused a
//! failure that *was* visible on the pull request, and pointed the reader away
//! from the one place the answer was sitting. A wrong claim that discourages
//! checking is worse than no claim.
//!
//! The fix is to derive the sentence from `github.event_name`, which in a
//! `workflow_call` describes the caller's own triggering event. This test pins
//! that: the reporter must key on the event, and must not carry a blanket claim
//! about pull-request coverage while any caller runs on `pull_request`.
//!
//! Modelled on `sweep_filter_invariant.rs` and `coderabbit_config_paths.rs`,
//! which make the same class of silent config rot loud.

use std::path::{Path, PathBuf};

/// The reusable reporter under test.
const REPORTER: &str = "report-failure.yml";

/// Phrases that assert, as a fixed fact, that the failure could not have been
/// seen before a merge. Any of these in the reporter is the defect this file
/// exists to catch — they can only be true of *some* callers.
///
/// Kept as fragments rather than the whole sentence so a light rewording cannot
/// slip past: the defect is the CLAIM, not its exact punctuation.
const BLANKET_CLAIMS: &[&str] = &[
    "does not run on pull requests",
    "could not have been caught before merge",
    "cannot fail on a pull request",
];

fn workflows_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(".github")
        .join("workflows")
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
}

/// The reporter with comment lines removed.
///
/// The guard below asks what the workflow **emits**, and a comment emits
/// nothing. Scanning the raw file conflates the two, and not in a theoretical
/// way: the fix for this very defect quotes the offending sentence in a comment
/// to explain why it was removed, so a raw scan would fail the file for
/// documenting its own correction. That is the `SELF`-exclusion problem
/// `sweep_filter_invariant.rs` records, in a different costume.
///
/// One rule covers both layers: a line whose first non-whitespace character is
/// `#` is a comment in YAML and in the `run:` block's bash alike.
fn emitted_text(path: &Path) -> String {
    read(path)
        .lines()
        .filter(|line| !line.trim_start().starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Every workflow that calls the reporter, paired with whether it can run on a
/// pull request.
///
/// The `on:` scan is deliberately crude — it reads the block between `on:` and
/// the next top-level key — because the alternative is a YAML dependency for one
/// boolean. It is a floor: a caller whose trigger it cannot parse counts as
/// non-PR, which biases toward *not* firing this test. That is the safe
/// direction for a guard whose failure message accuses the reporter.
fn callers() -> Vec<(String, bool)> {
    let dir = workflows_dir();
    let mut out = Vec::new();
    for entry in std::fs::read_dir(&dir).expect("read .github/workflows") {
        let path = entry.expect("dir entry").path();
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_string();
        if name == REPORTER || !name.ends_with(".yml") {
            continue;
        }
        let body = read(&path);
        if !body.contains(REPORTER) {
            continue;
        }
        let mut in_on = false;
        let mut runs_on_pr = false;
        for line in body.lines() {
            if line.starts_with("on:") {
                in_on = true;
                continue;
            }
            // A new top-level key ends the `on:` block.
            if in_on
                && !line.starts_with(char::is_whitespace)
                && !line.trim().is_empty()
                && !line.trim_start().starts_with('#')
            {
                in_on = false;
            }
            if in_on && line.contains("pull_request") {
                runs_on_pr = true;
            }
        }
        out.push((name, runs_on_pr));
    }
    out.sort();
    out
}

#[test]
fn the_reporter_has_callers_to_speak_for() {
    let callers = callers();
    assert!(
        !callers.is_empty(),
        "no workflow calls {REPORTER}. Either the reporter is dead — delete it and this \
         test — or the scan stopped matching, in which case every assertion below is \
         vacuous and this test is proving nothing."
    );
}

#[test]
fn at_least_one_caller_can_fail_on_a_pull_request() {
    let callers = callers();
    let pr_callers: Vec<&String> = callers
        .iter()
        .filter(|(_, pr)| *pr)
        .map(|(name, _)| name)
        .collect();
    assert!(
        !pr_callers.is_empty(),
        "no caller of {REPORTER} runs on `pull_request`, so the blanket claim this test \
         forbids would now be TRUE. That is a legitimate state — but it means the guard \
         below has nothing to guard, so decide deliberately rather than leaving a test \
         that passes for the wrong reason. Callers: {callers:?}"
    );
}

#[test]
fn the_reporter_does_not_assert_that_failures_were_uncatchable() {
    let body = emitted_text(&workflows_dir().join(REPORTER));
    let callers = callers();
    let pr_callers: Vec<&str> = callers
        .iter()
        .filter(|(_, pr)| *pr)
        .map(|(name, _)| name.as_str())
        .collect();

    for claim in BLANKET_CLAIMS {
        assert!(
            !body.contains(claim),
            "{REPORTER} contains the blanket claim {claim:?}, but {n} of its callers run \
             on `pull_request` ({pr_callers:?}). For those, a failure IS visible before \
             merge, and the sentence sends the reader away from the pull request that \
             carries the answer.\n\n\
             Derive the sentence from `github.event_name` instead — in a `workflow_call` \
             that is the caller's own triggering event, so it cannot go stale when a \
             caller gains or loses a trigger.",
            n = pr_callers.len(),
        );
    }
}

#[test]
fn the_reporter_derives_provenance_from_the_triggering_event() {
    let body = read(&workflows_dir().join(REPORTER));
    assert!(
        body.contains("github.event_name"),
        "{REPORTER} no longer reads `github.event_name`. The provenance sentence in the \
         issue body must be derived from this run's triggering event rather than asserted \
         about the workflow — that is the whole correction, and without it the text can \
         silently become false again the next time a caller changes its triggers."
    );
    assert!(
        body.contains("pull_request"),
        "{REPORTER} reads `github.event_name` but never compares it against \
         `pull_request`, so both branches of the provenance sentence cannot be reachable."
    );
}
