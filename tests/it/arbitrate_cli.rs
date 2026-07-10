//! CLI surface tests for `ferro arbitrate` (Task 10).
//!
//! `ferro arbitrate` always requires a prepared reference directory (the
//! normalize/projection step it runs ferro's own output through needs one),
//! so a full end-to-end run cannot be exercised without a manifest. Two kinds
//! of test here:
//!
//! - spawned-binary surface tests (clap validation, exit codes) that need no
//!   reference data at all;
//! - a manifest-gated smoke test (skipped, not failed, when no manifest is
//!   configured) that exercises the real `--other-output` + `--format json`
//!   / `--format text` path end-to-end against a live prepared reference.
//!
//! No live Mutalyzer calls are made anywhere in this file — every
//! manifest-gated case uses `--other-output` to avoid the network.

use std::path::PathBuf;
use std::process::Command;

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

// ===== spawned-binary surface tests (no reference data needed) =====

#[test]
fn arbitrate_is_listed_in_help() {
    let out = ferro().arg("--help").output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("arbitrate"),
        "top-level --help should list the arbitrate subcommand: {stdout}"
    );
}

#[test]
fn arbitrate_help_documents_its_flags() {
    let out = ferro().args(["arbitrate", "--help"]).output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    for flag in [
        "--reference",
        "--other-output",
        "--other-tool",
        "--mutalyzer-url",
        "--format",
        "--bug-report",
        "--no-open",
        "--include-environment",
        "--notes",
    ] {
        assert!(
            stdout.contains(flag),
            "arbitrate --help missing {flag}: {stdout}"
        );
    }
}

#[test]
fn arbitrate_accepts_bug_report_flags_at_clap_level() {
    // Clap-level acceptance only (no reference data needed): --bug-report,
    // --no-open, --include-environment, and --notes must parse without
    // clap rejecting them, even though this specific invocation fails later
    // for an unrelated reason (bad --reference dir). If clap didn't know
    // these flags it would fail immediately with an "unexpected argument"
    // error instead.
    let out = ferro()
        .args([
            "arbitrate",
            "NM_003002.4:c.274G>T",
            "--reference",
            "/nonexistent-xyz-ferro-arbitrate-bug-report",
            "--other-output",
            "NM_003002.4:c.274G>T",
            "--bug-report",
            "--no-open",
            "--include-environment",
            "--notes",
            "seen while triaging",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("unexpected argument") && !stderr.contains("unrecognized"),
        "clap should accept the bug-report flags: stderr={stderr}"
    );
}

#[test]
fn arbitrate_requires_reference() {
    // `--reference` is `required = true` in clap, so omitting it must fail
    // before any reference I/O is attempted.
    let out = ferro()
        .args([
            "arbitrate",
            "NM_003002.4:c.274G>T",
            "--other-output",
            "NM_003002.4:c.274G>T",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success(), "missing --reference must fail");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("reference"), "stderr: {stderr}");
}

#[test]
fn arbitrate_bad_reference_dir_exits_nonzero() {
    // A reference directory with no manifest.json is a hard error (nonzero
    // exit), not a panic.
    let out = ferro()
        .args([
            "arbitrate",
            "NM_003002.4:c.274G>T",
            "--reference",
            "/nonexistent-xyz-ferro-arbitrate",
            "--other-output",
            "NM_003002.4:c.274G>T",
        ])
        .output()
        .unwrap();
    assert!(!out.status.success());
}

// ===== manifest-gated end-to-end smoke test =====

/// `FERRO_MANIFEST`, when set, is authoritative — no fallback to
/// well-known paths. Same convention as the conformance test suites (see
/// e.g. `tests/it/mutalyzer_normalize_tests.rs`), so `FERRO_MANIFEST` can
/// explicitly disable this test even on a host with a reference mounted.
fn manifest_path() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("FERRO_MANIFEST") {
        let p = PathBuf::from(path);
        return if p.exists() { Some(p) } else { None };
    }
    let p = PathBuf::from("benchmark-output/manifest.json");
    if p.exists() {
        return Some(p);
    }
    None
}

/// The prepared reference *directory* `ferro arbitrate --reference` expects
/// (it joins `manifest.json` itself), derived from `manifest_path()`.
fn reference_dir() -> Option<PathBuf> {
    manifest_path().and_then(|p| p.parent().map(|d| d.to_path_buf()))
}

#[test]
fn arbitrate_json_reports_equivalent_for_identical_other_output() {
    let Some(reference) = reference_dir() else {
        println!(
            "arbitrate_cli: skipping — no manifest at FERRO_MANIFEST or benchmark-output/manifest.json"
        );
        return;
    };

    // NM_003002.4:c.274G>T is a plain substitution: it normalizes to itself,
    // so feeding it back as `--other-output` must yield an Equivalent verdict
    // regardless of the specific reference's transcript/genome content, as
    // long as the accession resolves (which the worktree CLAUDE.md documents
    // it does for the blessed/backup prepared reference).
    let variant = "NM_003002.4:c.274G>T";
    let out = ferro()
        .args(["arbitrate", variant, "--reference"])
        .arg(&reference)
        .args(["--other-output", variant, "--format", "json"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "arbitrate should succeed: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let json: serde_json::Value =
        serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not valid JSON: {e}\n{stdout}"));
    assert_eq!(json["verdict"], "equivalent", "full output: {json}");
    assert_eq!(json["compliance"], "not_applicable", "full output: {json}");
    assert_eq!(json["ferro_output"], variant, "full output: {json}");
    assert_eq!(json["other"]["output"], variant, "full output: {json}");
}

#[test]
fn arbitrate_text_verdict_line_leads_the_output() {
    let Some(reference) = reference_dir() else {
        println!(
            "arbitrate_cli: skipping — no manifest at FERRO_MANIFEST or benchmark-output/manifest.json"
        );
        return;
    };

    let variant = "NM_003002.4:c.274G>T";
    let out = ferro()
        .args(["arbitrate", variant, "--reference"])
        .arg(&reference)
        .args(["--other-output", variant, "--format", "text"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    let first_line = stdout.lines().next().unwrap_or_default();
    assert!(
        first_line.starts_with("VERDICT: Equivalent"),
        "verdict line must lead the text output: {stdout}"
    );
}

#[test]
fn arbitrate_basis_mismatch_on_version_difference_suggests_aligning_versions() {
    let Some(reference) = reference_dir() else {
        println!(
            "arbitrate_cli: skipping — no manifest at FERRO_MANIFEST or benchmark-output/manifest.json"
        );
        return;
    };

    // Same transcript, two different versions (NM_003002.2 vs .4): both are
    // real, resolvable accessions in the prepared reference (see the
    // worktree CLAUDE.md's NM_003002.2 note), so this is a genuine
    // version-basis mismatch rather than a bug-report-worthy disagreement.
    let out = ferro()
        .args(["arbitrate", "NM_003002.4:c.274G>T", "--reference"])
        .arg(&reference)
        .args(["--other-output", "NM_003002.2:c.274G>T", "--format", "text"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("VERDICT: BasisMismatch") && stdout.contains("align versions"),
        "expected a version-alignment hint, got: {stdout}"
    );
}

// ===== --bug-report gating (manifest-gated) =====
//
// A reliably ferro-implicated case (compliance/category naming the other
// tool, or a ferro parse error) requires a live disagreement against the
// real prepared reference; the committed `baseline-failures/*.txt` live-FAIL
// snapshot is currently empty (see its README — the snapshot only exists
// when a manifest run has found one), so there is no fixed input known to
// reproduce one deterministically here. `ferro_is_implicated`'s full
// true/false matrix is covered directly (and reference-free) by the unit
// tests in `src/bin/ferro.rs`'s `arbitrate_bug_report_gate_tests` module;
// this file sticks to what a real reference run can assert reliably: the
// not-implicated path, end to end through the CLI.

#[test]
fn arbitrate_bug_report_on_equivalent_case_declines_to_file() {
    let Some(reference) = reference_dir() else {
        println!(
            "arbitrate_cli: skipping — no manifest at FERRO_MANIFEST or benchmark-output/manifest.json"
        );
        return;
    };

    // Same case as arbitrate_json_reports_equivalent_for_identical_other_output:
    // feeding a variant back as its own --other-output always yields
    // Verdict::Equivalent / Compliance::NotApplicable, which
    // ferro_is_implicated must treat as not-implicated.
    let variant = "NM_000255.4:c.446dup";
    let out = ferro()
        .args(["arbitrate", variant, "--reference"])
        .arg(&reference)
        .args(["--other-output", variant, "--bug-report", "--no-open"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "arbitrate --bug-report should still succeed on a not-implicated verdict: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("No ferro bug indicated"),
        "expected the not-implicated note on stderr: stderr={stderr}"
    );
    assert!(
        !stdout.contains("github.com/fulcrumgenomics/ferro-hgvs/issues/new")
            && !stderr.contains("github.com/fulcrumgenomics/ferro-hgvs/issues/new"),
        "not-implicated case must not file: stdout={stdout} stderr={stderr}"
    );
}

#[test]
fn arbitrate_bug_report_json_mode_keeps_stdout_pipeable() {
    let Some(reference) = reference_dir() else {
        println!(
            "arbitrate_cli: skipping — no manifest at FERRO_MANIFEST or benchmark-output/manifest.json"
        );
        return;
    };

    let variant = "NM_000255.4:c.446dup";
    let out = ferro()
        .args(["arbitrate", variant, "--reference"])
        .arg(&reference)
        .args([
            "--other-output",
            variant,
            "--format",
            "json",
            "--bug-report",
            "--no-open",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    // stdout must be exactly the arbitration JSON — nothing else — so a
    // `--format json | ferro bug-report --from-arbitration -` pipe still
    // works when combined with --bug-report.
    let json: serde_json::Value =
        serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not valid JSON: {e}\n{stdout}"));
    assert_eq!(json["verdict"], "equivalent", "full output: {json}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("No ferro bug indicated"),
        "expected the not-implicated note on stderr: stderr={stderr}"
    );
}
