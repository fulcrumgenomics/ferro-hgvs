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
fn arbitrate_help_documents_the_exit_status() {
    // A script branching on the exit status has to be able to learn the
    // contract from --help: which verdicts exit 0, and which numbers the
    // two "no verdict" results exit with.
    let out = ferro().args(["arbitrate", "--help"]).output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8_lossy(&out.stdout);
    for needle in [
        "Exit status:",
        "3  Inconclusive",
        "4  FerroParseError",
        "OtherUnparseable",
    ] {
        assert!(
            stdout.contains(needle),
            "arbitrate --help missing {needle:?}: {stdout}"
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
    // A runtime error (the bad reference), not clap's usage error (2).
    assert_eq!(
        out.status.code(),
        Some(1),
        "stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("unexpected argument") && !stderr.contains("unrecognized"),
        "clap should accept the bug-report flags: stderr={stderr}"
    );
}

#[test]
fn bug_report_only_flags_are_refused_without_bug_report() {
    // `--no-open`, `--include-environment` and `--notes` are read only on the
    // `--bug-report` path. Without `requires`, `--notes X` alone would be
    // silently dropped; clap must refuse each one before any reference I/O.
    for extra in [
        vec!["--no-open"],
        vec!["--include-environment"],
        vec!["--notes", "seen while triaging"],
    ] {
        let mut args = vec![
            "arbitrate",
            "NM_003002.4:c.274G>T",
            "--reference",
            "/nonexistent-xyz-ferro-arbitrate-bug-report",
            "--other-output",
            "NM_003002.4:c.274G>T",
        ];
        args.extend(extra.iter().copied());
        let out = ferro().args(&args).output().unwrap();
        assert_eq!(
            out.status.code(),
            Some(2),
            "{extra:?} without --bug-report must be a usage error"
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains("--bug-report"),
            "{extra:?}: clap should name the required --bug-report flag: {stderr}"
        );
    }
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
    assert_eq!(
        out.status.code(),
        Some(2),
        "missing --reference must be a usage error"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("reference"), "stderr: {stderr}");
}

#[test]
fn arbitrate_bad_reference_dir_exits_with_a_runtime_error() {
    // A reference directory with no manifest.json is a runtime error (exit 1),
    // not a panic and not a usage error.
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
    assert_eq!(
        out.status.code(),
        Some(1),
        "stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
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

/// Say why a manifest-gated test is skipping and what still covers the exit
/// status without a reference, so a skip on PR CI is not read as coverage.
fn skip_without_manifest() {
    eprintln!(
        "arbitrate_cli: skipping — no manifest at FERRO_MANIFEST or \
         benchmark-output/manifest.json; the reference-free gates are the \
         clap-level tests in this file plus `arbitrate_exit_code_tests` and \
         `arbitrate_bug_report_gate_tests` in src/bin/ferro.rs"
    );
}

#[test]
fn arbitrate_json_reports_equivalent_for_identical_other_output() {
    let Some(reference) = reference_dir() else {
        skip_without_manifest();
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
        skip_without_manifest();
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
        skip_without_manifest();
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

// ===== exit status (manifest-gated) =====
//
// The exit status says whether ferro delivered a verdict, not what the
// verdict was: every delivered verdict exits 0 and the script reads
// `.verdict`; `Inconclusive` exits 3 and `FerroParseError` 4, each after the
// full result has been printed. `FerroParseError` means ferro's own output
// did not parse, which no input reproduces on purpose, so its code is pinned
// by the unit tests next to `arbitrate_exit_code` in `src/bin/ferro.rs`;
// the two cases below are the ones a real reference run can produce.

#[test]
fn arbitrate_inconclusive_exits_3_after_printing_the_result() {
    let Some(reference) = reference_dir() else {
        skip_without_manifest();
        return;
    };

    // `c.=` asserts the whole reference is unchanged and names no interval,
    // so the oracle cannot reduce it to SPDI: a comparison error, caught and
    // reported as Inconclusive with the reason. Deterministic on any
    // reference that resolves the accession, and independent of the
    // substitution on the ferro side.
    let variant = "NM_003002.4:c.274G>T";
    for format in ["json", "text"] {
        let out = ferro()
            .args(["arbitrate", variant, "--reference"])
            .arg(&reference)
            .args(["--other-output", "NM_003002.4:c.=", "--format", format])
            .output()
            .unwrap();
        assert_eq!(
            out.status.code(),
            Some(3),
            "--format {format}: Inconclusive must exit 3, got {:?}; stderr={}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        if format == "json" {
            let json: serde_json::Value = serde_json::from_str(&stdout)
                .unwrap_or_else(|e| panic!("not valid JSON: {e}\n{stdout}"));
            assert_eq!(json["verdict"], "inconclusive", "full output: {json}");
            assert!(
                json["reason"].as_str().is_some_and(|r| !r.is_empty()),
                "an Inconclusive result must carry its reason: {json}"
            );
        } else {
            let first_line = stdout.lines().next().unwrap_or_default();
            assert!(
                first_line.starts_with("INCONCLUSIVE:"),
                "the result must be printed in full before the non-zero exit: {stdout}"
            );
        }
    }
}

#[test]
fn arbitrate_other_unparseable_is_a_verdict_and_exits_0() {
    let Some(reference) = reference_dir() else {
        skip_without_manifest();
        return;
    };

    // The other side's output not being HGVS is something ferro found, not
    // something it failed to do: a delivered verdict, exit 0, read it from
    // the output like any other.
    let variant = "NM_003002.4:c.274G>T";
    let out = ferro()
        .args(["arbitrate", variant, "--reference"])
        .arg(&reference)
        .args(["--other-output", "not hgvs at all", "--format", "json"])
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(0),
        "OtherUnparseable is a delivered verdict: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let json: serde_json::Value =
        serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not valid JSON: {e}\n{stdout}"));
    assert_eq!(json["verdict"], "other_unparseable", "full output: {json}");
}

#[test]
fn arbitrate_unreachable_mutalyzer_is_an_error_not_a_verdict() {
    let Some(reference) = reference_dir() else {
        skip_without_manifest();
        return;
    };

    // Port 1 on the loopback refuses immediately: no network, no live
    // Mutalyzer call. A fetch that never reached the other tool used to
    // come back as `status: "unavailable"` with no output, which
    // `arbitrate()` reports as OtherUnparseable -- exit 0 and "could not be
    // parsed as HGVS" for a tool that was never asked. It is a failure to
    // run the arbitration: exit 1, nothing on stdout, and the error names
    // the offline fallback.
    let out = ferro()
        .args(["arbitrate", "NM_003002.4:c.274G>T", "--reference"])
        .arg(&reference)
        .args(["--mutalyzer-url", "http://127.0.0.1:1", "--format", "json"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(1),
        "an unreachable Mutalyzer must be a runtime error, got {:?}; stderr={stderr}",
        out.status
    );
    assert!(
        out.stdout.is_empty(),
        "no arbitration result must be printed: {}",
        String::from_utf8_lossy(&out.stdout)
    );
    assert!(
        stderr.contains("--other-output") && stderr.contains("127.0.0.1:1"),
        "the error must name the fallback flag and the URL it could not reach: {stderr}"
    );
}

#[test]
fn arbitrate_bug_report_on_inconclusive_declines_and_still_exits_3() {
    let Some(reference) = reference_dir() else {
        skip_without_manifest();
        return;
    };

    // `--bug-report` neither files on an Inconclusive result (#886) nor
    // changes its exit status: the status reports the arbitration, not the
    // side effect, so a script reads it the same way with or without the
    // flag. The JSON stays alone on stdout and the note goes to stderr.
    let variant = "NM_003002.4:c.274G>T";
    let out = ferro()
        .args(["arbitrate", variant, "--reference"])
        .arg(&reference)
        .args([
            "--other-output",
            "NM_003002.4:c.=",
            "--format",
            "json",
            "--bug-report",
            "--no-open",
        ])
        .output()
        .unwrap();
    assert_eq!(
        out.status.code(),
        Some(3),
        "Inconclusive must exit 3 with --bug-report too, got {:?}; stderr={}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let json: serde_json::Value =
        serde_json::from_str(&stdout).unwrap_or_else(|e| panic!("not valid JSON: {e}\n{stdout}"));
    assert_eq!(json["verdict"], "inconclusive", "full output: {json}");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("No ferro bug indicated"),
        "expected the not-implicated note on stderr: stderr={stderr}"
    );
    assert!(
        !stderr.contains("github.com/fulcrumgenomics/ferro-hgvs/issues/new"),
        "an Inconclusive result must not be filed: stderr={stderr}"
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
        skip_without_manifest();
        return;
    };

    // Feeding a variant back as its own --other-output reduces both sides to
    // the same SPDI edit, so the verdict is always Equivalent. This uses the
    // plain substitution from
    // arbitrate_json_reports_equivalent_for_identical_other_output rather than
    // a dup: a duplication also runs the dup-versus-insertion compliance
    // predicate, which can name a tool and would make this a ferro-implicated
    // case for reasons unrelated to the gating logic under test.
    let variant = "NM_003002.4:c.274G>T";
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
        skip_without_manifest();
        return;
    };

    // A substitution, for the same reason as the test above.
    let variant = "NM_003002.4:c.274G>T";
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
