//! CLI surface tests for `ferro bug-report` (Task 12).
//!
//! `ferro bug-report` consumes a completed `ferro arbitrate --format json`
//! bundle; no live reference, network, or browser launch is exercised here.
//! Every test pipes a hand-built `Arbitration` JSON into the spawned
//! binary's stdin via `--from-arbitration -` and passes `--no-open` so the
//! command never attempts to spawn a browser opener (`should_offer_to_open_
//! browser` in `src/bin/ferro.rs` also gates on `$BROWSER`/a TTY, but
//! `--no-open` alone is sufficient and doesn't depend on the test
//! environment).

use std::io::Write;
use std::process::{Command, Stdio};

use ferro_hgvs::arbitrate::spec_citations::SpecCitation;
use ferro_hgvs::arbitrate::{Arbitration, ArbitrationCategory, Compliance, OtherResult, Verdict};

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// A genuine ferro-is-wrong arbitration: mirrors `bug_report::tests::sample`
/// (the other tool's spelling is spec-compliant, ferro's is not).
fn sample_arbitration() -> Arbitration {
    Arbitration {
        input: "NC_000001.11:g.51dup".to_string(),
        verdict: Verdict::Different,
        compliance: Compliance::Other,
        category: ArbitrationCategory::MutalyzerCorrect,
        ferro_output: Some("NC_000001.11:g.51dup".to_string()),
        other: OtherResult {
            tool: "mutalyzer".to_string(),
            status: "ok".to_string(),
            output: Some("NC_000001.11:g.52dup".to_string()),
        },
        spec_citations: vec![SpecCitation {
            spec_version: "21.0.4".to_string(),
            file: "duplication.md".to_string(),
            heading: "Duplication".to_string(),
            excerpt: "A duplication must be described using the shortest possible 3' \
                      representation."
                .to_string(),
        }],
        ferro_spdi: Some("NC_000001.11:50:A:AA".to_string()),
        other_spdi: Some("NC_000001.11:51:A:AA".to_string()),
    }
}

/// Run `ferro bug-report` with the given extra args, piping `stdin` in and
/// capturing combined stdout/stderr as separate strings.
fn run_bug_report(stdin: &str, extra_args: &[&str]) -> (bool, String, String) {
    let mut child = ferro()
        .arg("bug-report")
        .args(extra_args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("failed to spawn ferro bug-report");
    child
        .stdin
        .take()
        .expect("child stdin")
        .write_all(stdin.as_bytes())
        .expect("failed to write to child stdin");
    let out = child.wait_with_output().expect("failed to wait on child");
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).to_string(),
        String::from_utf8_lossy(&out.stderr).to_string(),
    )
}

#[test]
fn bug_report_is_listed_in_help() {
    let out = ferro().arg("--help").output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("bug-report"),
        "top-level --help should list the bug-report subcommand: {stdout}"
    );
}

#[test]
fn bug_report_help_documents_its_flags() {
    let out = ferro().args(["bug-report", "--help"]).output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    for flag in [
        "--from-arbitration",
        "--category",
        "--notes",
        "--no-open",
        "--include-environment",
    ] {
        assert!(
            stdout.contains(flag),
            "bug-report --help missing {flag}: {stdout}"
        );
    }
}

#[test]
fn bug_report_requires_from_arbitration() {
    let out = ferro().args(["bug-report", "--no-open"]).output().unwrap();
    assert!(
        !out.status.success(),
        "omitting --from-arbitration must fail"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--from-arbitration"),
        "stderr should explain the missing bundle: {stderr}"
    );
}

#[test]
fn bug_report_bad_json_exits_nonzero_not_a_panic() {
    let (success, _stdout, stderr) =
        run_bug_report("not valid json", &["--from-arbitration", "-", "--no-open"]);
    assert!(!success, "malformed bundle JSON must fail cleanly");
    assert!(!stderr.contains("panicked"), "stderr: {stderr}");
}

#[test]
fn bug_report_missing_file_exits_nonzero_not_a_panic() {
    // A `--from-arbitration <path>` that does not exist on disk (as opposed
    // to `-` for stdin, exercised by the other tests here) must fail
    // cleanly via the `std::fs::read_to_string` error mapping, not panic.
    let (success, _stdout, stderr) = run_bug_report(
        "",
        &[
            "--from-arbitration",
            "/nonexistent/path/does-not-exist.json",
            "--no-open",
        ],
    );
    assert!(!success, "a nonexistent --from-arbitration path must fail");
    assert!(!stderr.contains("panicked"), "stderr: {stderr}");
    assert!(
        stderr.contains("failed to read arbitration bundle"),
        "stderr should explain the missing file: {stderr}"
    );
}

#[test]
fn bug_report_no_open_prints_prefilled_url_without_labels() {
    let json = serde_json::to_string(&sample_arbitration()).unwrap();
    let (success, stdout, stderr) =
        run_bug_report(&json, &["--from-arbitration", "-", "--no-open"]);
    assert!(success, "stderr={stderr}\nstdout={stdout}");

    assert!(
        stdout.contains("https://github.com/fulcrumgenomics/ferro-hgvs/issues/new?"),
        "stdout should contain the prefilled new-issue URL: {stdout}"
    );
    assert!(stdout.contains("title="), "stdout: {stdout}");
    assert!(stdout.contains("body="), "stdout: {stdout}");
    assert!(
        !stdout.contains("labels="),
        "labels= 404s for read-only external users and must never appear: {stdout}"
    );
    assert!(
        !stdout.contains("assignees="),
        "assignees= must never appear either: {stdout}"
    );

    // --no-open must not require any confirmation input: this test writes
    // only the JSON to stdin (no trailing "y"/"n" line), yet the command
    // still succeeds and prints the URL.
    assert!(
        stdout.contains("arbitration: NC_000001.11:g.51dup"),
        "the rendered issue title should appear in stdout: {stdout}"
    );
}

#[test]
fn bug_report_body_includes_verdict_and_omits_environment_by_default() {
    let json = serde_json::to_string(&sample_arbitration()).unwrap();
    let (success, stdout, stderr) =
        run_bug_report(&json, &["--from-arbitration", "-", "--no-open"]);
    assert!(success, "stderr={stderr}\nstdout={stdout}");
    assert!(stdout.contains("mutalyzer"), "stdout: {stdout}");
    assert!(
        !stdout.to_lowercase().contains("os:"),
        "environment details must be opt-in: {stdout}"
    );
}

#[test]
fn bug_report_include_environment_adds_environment_section() {
    let json = serde_json::to_string(&sample_arbitration()).unwrap();
    let (success, stdout, stderr) = run_bug_report(
        &json,
        &[
            "--from-arbitration",
            "-",
            "--no-open",
            "--include-environment",
        ],
    );
    assert!(success, "stderr={stderr}\nstdout={stdout}");
    assert!(
        stdout.to_lowercase().contains("os:"),
        "--include-environment should add the environment section: {stdout}"
    );
}

#[test]
fn bug_report_notes_are_rendered_verbatim() {
    let json = serde_json::to_string(&sample_arbitration()).unwrap();
    let (success, stdout, stderr) = run_bug_report(
        &json,
        &[
            "--from-arbitration",
            "-",
            "--no-open",
            "--notes",
            "seen while normalizing a ClinVar batch",
        ],
    );
    assert!(success, "stderr={stderr}\nstdout={stdout}");
    assert!(
        stdout.contains("seen while normalizing a ClinVar batch"),
        "stdout: {stdout}"
    );
}

#[test]
fn bug_report_mismatched_category_warns_but_still_succeeds() {
    let json = serde_json::to_string(&sample_arbitration()).unwrap();
    let (success, stdout, stderr) = run_bug_report(
        &json,
        &[
            "--from-arbitration",
            "-",
            "--no-open",
            "--category",
            "ferro_correct",
        ],
    );
    assert!(success, "stderr={stderr}\nstdout={stdout}");
    assert!(
        stdout.contains("https://github.com/fulcrumgenomics/ferro-hgvs/issues/new?"),
        "a category mismatch should only warn, not block the report: {stdout}"
    );
    assert!(
        stderr.contains("--category"),
        "a mismatched --category should produce a warning: {stderr}"
    );
}
