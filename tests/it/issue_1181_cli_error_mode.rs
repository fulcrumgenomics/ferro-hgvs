//! CLI surface test for `ferro normalize --error-mode` (#1181).
//!
//! `--error-mode` reached the *preprocessor* but never the *normalizer*, so
//! every normalizer-level rejection was unreachable from the CLI: strict, lenient
//! and silent all normalized leniently. The same `ErrorConfig` also carries the
//! `--ignore` / `--reject` overrides, so those were inert for normalization too.
//!
//! **These tests must drive the binary, not the library.** Strict mode is already
//! well covered at the library level (`error_mode_tests.rs`,
//! `issue_336_position_past_end.rs`, `issue_486_position_out_of_bounds.rs`,
//! `issue_393_mt_w4004.rs`) and was correct there the whole time — the defect
//! lived purely in the CLI seam. A library-level regression test would have
//! passed against the bug. Before this fix, exactly one test in the suite passed
//! `--error-mode` to `normalize` (`cli_parallel_normalize.rs`), and it passed
//! `silent` only to keep output quiet, asserting nothing about mode behavior.
//!
//! For the same reason the fixture is a real on-disk reference rather than the
//! mock provider: under the mock, `normalize` echoes its input and no
//! reference-base validation runs, so all modes are trivially identical and the
//! test would prove nothing. A genomic (`g.`) input needs only an indexed genome
//! FASTA — no cdot — which keeps the fixture small enough for CI with no external
//! data.
//!
//! Scope note: `RefSeqMismatch` (W5001) is the one normalizer-level category
//! reachable from a genome-only fixture. `PositionPastEnd` (W4004) is gated by
//! `check_cds_pos_past_end`, which is CDS-specific and does not fire for a `g.`
//! coordinate past a contig end, so it is not asserted here — it is covered at
//! the library level. That is sufficient, because the CLI seam is a *single*
//! wiring point: one category flowing through plus both override directions
//! proves the whole `ErrorConfig` arrives.

use std::path::Path;
use std::process::Command;

use ferro_hgvs::prepare::manifest::ReferenceManifest;

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// The single contig in the fixture: 60 bases cycling `ACGT`.
const CONTIG: &str = "NC_000001.11";
const CONTIG_LEN: usize = 60;

/// The base at 1-based position `pos` of the fixture contig. Computed rather than
/// hardcoded so the sequence and the expectations cannot drift apart.
fn base_at(pos: usize) -> char {
    ['A', 'C', 'G', 'T'][(pos - 1) % 4]
}

/// Write a minimal reference dir: an **indexed** genome FASTA with one small
/// contig, plus a manifest pointing at it. The `.fai` is required —
/// `MultiFastaProvider` refuses a reference with no indexed FASTA.
fn minimal_genome_reference(dir: &Path) {
    let seq: String = (1..=CONTIG_LEN).map(base_at).collect();
    let header = format!(">{CONTIG}\n");
    std::fs::write(dir.join("genome.fa"), format!("{header}{seq}\n")).unwrap();
    // `name\tlength\toffset\tlinebases\tlinewidth` — the sequence is a single
    // line, so linebases == length and linewidth == length + 1 (the newline).
    std::fs::write(
        dir.join("genome.fa.fai"),
        format!(
            "{CONTIG}\t{CONTIG_LEN}\t{}\t{CONTIG_LEN}\t{}\n",
            header.len(),
            CONTIG_LEN + 1
        ),
    )
    .unwrap();

    let mut manifest = ReferenceManifest {
        reference_dir: dir.to_path_buf(),
        genome_fasta: Some(std::path::PathBuf::from("genome.fa")),
        transcript_count: 0,
        available_prefixes: vec!["NC_".to_string()],
        ..Default::default()
    };
    manifest.save().unwrap();
}

/// A substitution whose *stated* reference base contradicts the actual base at
/// position 10, so it trips `RefSeqMismatch` in the normalizer.
fn mismatching_substitution() -> String {
    let actual = base_at(10);
    let wrong = if actual == 'A' { 'C' } else { 'A' };
    let alt = if actual == 'T' { 'G' } else { 'T' };
    format!("{CONTIG}:g.10{wrong}>{alt}")
}

/// Run `ferro normalize` with the given extra args; returns (success, combined
/// output).
fn normalize_with(reference: &str, variant: &str, extra: &[&str]) -> (bool, String) {
    let mut cmd = ferro();
    cmd.args(["normalize", "--reference", reference]);
    cmd.args(extra);
    cmd.arg(variant);
    let out = cmd.output().unwrap();
    let combined = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    (out.status.success(), combined)
}

/// Set up the fixture and return (tempdir guard, reference path, variant).
fn fixture() -> (tempfile::TempDir, String, String) {
    let dir = tempfile::tempdir().unwrap();
    minimal_genome_reference(dir.path());
    let reference = dir.path().to_str().unwrap().to_string();
    (dir, reference, mismatching_substitution())
}

/// The heart of the regression: a reference-base mismatch must be *rejected*
/// under `--error-mode strict` and *accepted* under `lenient`. Before the fix
/// both printed the same successful output.
#[test]
fn strict_rejects_a_reference_mismatch_that_lenient_accepts() {
    let (_dir, reference, variant) = fixture();

    let (lenient_ok, lenient_out) =
        normalize_with(&reference, &variant, &["--error-mode", "lenient"]);
    assert!(
        lenient_ok,
        "lenient must accept a reference mismatch; output: {lenient_out}"
    );

    let (strict_ok, strict_out) = normalize_with(&reference, &variant, &["--error-mode", "strict"]);
    assert!(
        !strict_ok,
        "strict must reject a reference mismatch, but it exited 0 — this is the \
         #1181 bug: --error-mode never reached the normalizer. output: {strict_out}"
    );
    assert!(
        strict_out.contains("Reference mismatch"),
        "the strict rejection should name the reference mismatch; got: {strict_out}"
    );
}

/// `silent` must also accept rather than reject. Pinning all three modes guards
/// against a "fix" that hardwires strict everywhere.
///
/// Note what is and is not distinct here. `strict` is distinguishable — it
/// rejects, and names the mismatch. `lenient` and `silent` are **not**: their
/// combined stdout+stderr is byte-identical. That equality is asserted rather
/// than glossed, because the tempting assertion — "silent omits the diagnostic
/// that lenient prints" — would keep passing if the diagnostic disappeared
/// everywhere.
///
/// **What that equality now contains has changed.** It used to hold because
/// neither mode printed a normalizer diagnostic at all: `ferro normalize` was on
/// `Normalizer::normalize`, the exit that discards warnings. Since the CLI moved
/// to `normalize_with_diagnostics`, both modes print
/// `warning[REFSEQ_MISMATCH]`, so the equality is now "both disclose it
/// identically" rather than "neither discloses anything". The assertion below
/// pins the disclosure explicitly, so the equality can never again be satisfied
/// by two silences.
///
/// That `silent` prints it is a genuine open question, recorded rather than
/// decided here. It used to be stated against `ErrorMode::emits_warnings()`,
/// which "says only `Lenient` should" — **#1629 removed that predicate**, since
/// it had no call site in `src/` and could not have acquired one without
/// overruling the per-code `--ignore`/`--reject` overrides. So the question
/// survives, in a sharper form: `REFSEQ_MISMATCH` (W5001) reaches the CLI's
/// diagnostic printer unconditionally rather than through
/// `NormalizeConfig::should_warn_ref_mismatch()` — which *is* the per-code
/// authority, *does* answer `false` under silent mode, and is itself unwired
/// (the sole `warnings.push(NormalizationWarning::RefSeqMismatch { … })` in
/// `src/normalize/mod.rs` sits in `normalize_na_edit`'s `if !validation.valid`
/// arm, guarded by nothing but the mismatch itself; only
/// `should_reject_ref_mismatch` is consulted, at the strict gate — cited by
/// expression rather than by line, because line numbers here drift).
/// `ferro project`
/// (#1182) prints normalizer warnings in every mode for the same reason. The
/// two sibling commands are consistent with each other, which is the property
/// this change preserved.
#[test]
fn strict_is_distinct_from_lenient_and_silent_at_the_cli() {
    let (_dir, reference, variant) = fixture();

    let (strict_ok, strict_out) = normalize_with(&reference, &variant, &["--error-mode", "strict"]);
    let (lenient_ok, lenient_out) =
        normalize_with(&reference, &variant, &["--error-mode", "lenient"]);
    let (silent_ok, silent_out) = normalize_with(&reference, &variant, &["--error-mode", "silent"]);

    assert!(!strict_ok, "strict must reject");
    assert!(lenient_ok, "lenient must accept");
    assert!(silent_ok, "silent must accept");
    // Both accepting modes surface the same normalizer diagnostic, so pin the
    // equality as the fact it is. If they ever diverge, this fails and the
    // distinction gets asserted properly instead of assumed.
    assert_eq!(
        lenient_out, silent_out,
        "lenient and silent are indistinguishable on this input; if that changes, \
         assert the actual difference rather than relaxing this"
    );
    // …and pin what they agree *on*, so the equality above cannot degenerate
    // back into two silences without failing here.
    assert!(
        lenient_out.contains("warning[REFSEQ_MISMATCH]"),
        "an accepted reference mismatch must be disclosed, not silently taken; got: {lenient_out}"
    );
    assert!(
        strict_out.contains("Reference mismatch"),
        "only strict is distinct, and it must say why; got: {strict_out}"
    );
}

/// The documented default is `strict` (`#[arg(long, default_value = "strict")]`),
/// so omitting the flag must be identical to passing it. This is the half that
/// silently behaved as lenient since the initial commit.
#[test]
fn the_default_mode_is_strict() {
    let (_dir, reference, variant) = fixture();

    let (default_ok, default_out) = normalize_with(&reference, &variant, &[]);
    let (strict_ok, strict_out) = normalize_with(&reference, &variant, &["--error-mode", "strict"]);

    assert!(
        !default_ok,
        "the default must reject like strict; output: {default_out}"
    );
    assert_eq!(
        (default_ok, default_out),
        (strict_ok, strict_out),
        "omitting --error-mode must be identical to --error-mode strict",
    );
}

/// `--reject` is folded into the same `ErrorConfig`, so it was inert for
/// normalization too. Promoting a warning to a rejection must work even from the
/// lenient base mode.
#[test]
fn reject_override_promotes_a_warning_from_lenient_base_mode() {
    let (_dir, reference, variant) = fixture();

    let (plain_ok, _) = normalize_with(&reference, &variant, &["--error-mode", "lenient"]);
    assert!(plain_ok, "premise: plain lenient must accept");

    let (rejected_ok, rejected_out) = normalize_with(
        &reference,
        &variant,
        &["--error-mode", "lenient", "--reject", "W5001"],
    );
    assert!(
        !rejected_ok,
        "--reject W5001 must reject from a lenient base mode; output: {rejected_out}"
    );
    // Pin *why* it failed: the rejection must come from the promoted
    // reference-mismatch rule, not from the command erroring for some unrelated
    // reason, which a bare `!rejected_ok` cannot distinguish.
    assert!(
        rejected_out.contains("Reference mismatch"),
        "the rejection must be the promoted W5001 mismatch; got: {rejected_out}"
    );
}

/// The mirror direction: `--ignore` must demote a category that strict would
/// otherwise reject.
#[test]
fn ignore_override_demotes_a_rejection_from_strict_base_mode() {
    let (_dir, reference, variant) = fixture();

    let (plain_ok, _) = normalize_with(&reference, &variant, &["--error-mode", "strict"]);
    assert!(!plain_ok, "premise: plain strict must reject");

    let (ignored_ok, ignored_out) = normalize_with(
        &reference,
        &variant,
        &["--error-mode", "strict", "--ignore", "W5001"],
    );
    assert!(
        ignored_ok,
        "--ignore W5001 must let a strict run succeed; output: {ignored_out}"
    );
}
