//! The failure reporter must file something even when its lookup fails.
//!
//! `report-failure.yml` runs under `set -euo pipefail` and opens by searching
//! for the rolling issue to comment on. That search is a network call to the
//! GitHub search API, and it is the *first* thing the job does — so before this
//! was guarded, any transient failure of it (rate limit, API blip, a title the
//! search chokes on) aborted the step before the `gh issue create` fallback was
//! ever reached. The alerting path would then fail silently underneath the very
//! failure it exists to announce.
//!
//! That is the same shape as the defect this workflow was already carrying: a
//! path nobody exercises, whose broken state is indistinguishable from its
//! working state. Asserting the `GH_REPO` line is present does not touch it,
//! because the question is not what the file says but what the script *does*.
//!
//! So this test runs the script. It extracts the `run:` block from the workflow
//! and executes it under `bash` with a stub `gh` on `PATH`, then checks which
//! command the script reached. Three scenarios, all of which must file
//! something:
//!
//! | lookup | expected |
//! |---|---|
//! | fails outright | warn, then **create** — a duplicate beats silence |
//! | succeeds, no match | **create** |
//! | succeeds, match | **comment** on the issue it found |
//!
//! Deliberately NOT asserted: that `gh issue create` or `gh issue comment`
//! failing still reddens the step. It does, and it should — at that point
//! nothing was recorded anywhere and the step's own red X is the last remaining
//! signal — but a test that pinned it would be pinning `set -e`, not this
//! script.

#![cfg(unix)]

use serde_yaml::Value;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The `run:` script of the reporter's only step, read out of the workflow.
///
/// Read from the workflow rather than duplicated here: a copy would keep
/// passing after the real script changed, which is the failure mode this whole
/// file is about.
fn reporter_script() -> String {
    let path = repo_root().join(".github/workflows/report-failure.yml");
    let text = std::fs::read_to_string(&path).expect("report-failure.yml is readable");
    let doc: Value = serde_yaml::from_str(&text).expect("report-failure.yml parses");
    doc.get("jobs")
        .and_then(|j| j.get("report"))
        .and_then(|j| j.get("steps"))
        .and_then(Value::as_sequence)
        .and_then(|steps| steps.iter().find_map(|s| s.get("run")))
        .and_then(Value::as_str)
        .expect("the `report` job has a `run` step")
        .to_string()
}

/// A `gh` that records which subcommand it was asked for instead of calling
/// GitHub, and can be told to fail the lookup.
fn write_stub_gh(dir: &Path) {
    let stub = dir.join("gh");
    std::fs::write(
        &stub,
        r#"#!/bin/bash
case "$1 $2" in
  "issue list")
      if [ "${GH_STUB_LOOKUP:-ok}" = "fail" ]; then
        echo "gh: API rate limit exceeded" >&2
        exit 1
      fi
      echo "${GH_STUB_FOUND:-}"
      ;;
  "issue comment") echo "STUB_COMMENTED $3" ;;
  "issue create")  echo "STUB_CREATED" ;;
  *) echo "STUB: unexpected: $*" >&2; exit 99 ;;
esac
"#,
    )
    .expect("write stub gh");
    let mut perms = std::fs::metadata(&stub).expect("stat stub").permissions();
    std::os::unix::fs::PermissionsExt::set_mode(&mut perms, 0o755);
    std::fs::set_permissions(&stub, perms).expect("chmod stub");
}

struct Outcome {
    status: Option<i32>,
    output: String,
}

fn run_reporter(lookup: &str, found: &str) -> Outcome {
    let dir = tempfile::tempdir().expect("tempdir");
    write_stub_gh(dir.path());
    let script = dir.path().join("step.sh");
    std::fs::write(&script, reporter_script()).expect("write script");

    let path = format!(
        "{}:{}",
        dir.path().display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let out = Command::new("bash")
        .arg(&script)
        .env("PATH", path)
        .env("GH_STUB_LOOKUP", lookup)
        .env("GH_STUB_FOUND", found)
        // The values the workflow would supply. `set -u` aborts on any unset
        // one, so every variable the script reads has to be here.
        .env("GH_TOKEN", "stub-token")
        .env("GH_REPO", "example/repo")
        .env("TITLE", "a stable failure-mode title")
        .env("SUMMARY", "something failed")
        .env("HINT", "reproduce with: cargo test")
        .env("RUN_URL", "https://example.invalid/run/1")
        .env("WORKFLOW", "Coverage")
        .env("SHA", "0123456789abcdef")
        .output()
        .expect("run the reporter script under bash");

    Outcome {
        status: out.status.code(),
        output: format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        ),
    }
}

#[test]
fn a_failed_lookup_still_files_an_issue() {
    let run = run_reporter("fail", "");
    assert_eq!(
        run.status,
        Some(0),
        "a failed rolling-issue lookup aborted the reporter instead of falling back. \
         The step dies before it can file anything, so the failure it was reporting \
         goes unrecorded — the defect this workflow exists to prevent.\n{}",
        run.output
    );
    assert!(
        run.output.contains("STUB_CREATED"),
        "the reporter exited cleanly but filed nothing after a failed lookup\n{}",
        run.output
    );
    assert!(
        run.output.contains("::warning::"),
        "the reporter degraded silently. A run that could not check for the rolling \
         issue may have filed a duplicate, and the log has to say so\n{}",
        run.output
    );
}

#[test]
fn no_existing_issue_is_created_fresh() {
    let run = run_reporter("ok", "");
    assert_eq!(run.status, Some(0), "{}", run.output);
    assert!(
        run.output.contains("STUB_CREATED"),
        "with no rolling issue found, the reporter must create one\n{}",
        run.output
    );
}

#[test]
fn an_existing_issue_is_commented_on_rather_than_duplicated() {
    let run = run_reporter("ok", "4242");
    assert_eq!(run.status, Some(0), "{}", run.output);
    assert!(
        run.output.contains("STUB_COMMENTED 4242"),
        "the reporter must comment on the rolling issue it found, not open a second \
         one — one issue per failure mode is the whole point of the lookup\n{}",
        run.output
    );
    assert!(
        !run.output.contains("STUB_CREATED"),
        "the reporter both commented and created\n{}",
        run.output
    );
}
