//! The normalizer's warning channel, from the normalizer out to an ordinary caller.
//!
//! `Normalizer::normalize` returns `Result<HgvsVariant, FerroError>` — a return
//! type with nowhere to put a warning. Its body is
//! `Ok(self.normalize_core_checked(variant)?.0)`, and the `.0` discards the
//! `Vec<NormalizationWarning>` the core hands back. Every warning the normalizer
//! raises therefore existed and was thrown away before any caller on that exit
//! could see it.
//!
//! That is not a cosmetic loss. Several of the repairs it reports are **lossy in
//! a way the output string does not record**, so a caller cannot recover the
//! fact from the result:
//!
//! - `MEMBERS_COALESCED_FROM_REPORTED_FORM` — the input described N cis members
//!   and the canonical form describes fewer. `DNA/delins.md:79-84` names exactly
//!   why that matters ("the two variants may have been reported (or might occur)
//!   individually"), and the provenance is unrecoverable from the normalized
//!   string.
//! - `INSERTED_SEQUENCE_EXPANDED` — a `ins[120_130]` reference-range payload was
//!   replaced by the literal bases it denotes. The output no longer says it came
//!   from a range.
//! - `REFSEQ_MISMATCH` — the caller's stated reference base disagreed with the
//!   reference. Under a non-strict error mode the description is accepted anyway.
//!
//! `ferro project` was given this channel by #1182 (`src/bin/ferro.rs`'s
//! `run_project` prints `warning[CODE]: message` for the normalization warnings
//! its projection carried). `ferro normalize` was not: its only `warning[...]`
//! prints iterated `preprocess_result.warnings`, which is the **preprocessor's**
//! population — a different phase entirely. So the two sibling commands
//! disagreed about whether a normalizer diagnostic is worth showing.
//!
//! These tests drive the **binary**, not the library, because the library was
//! never the broken half: `normalize_with_diagnostics` has always returned these
//! warnings and `coalesced_members_diagnostic.rs` has always asserted so. The
//! defect lived purely in which exit the callers were on, and a library-level
//! test passes against the bug.
//!
//! The last test is the other half of the contract, and the one that makes this
//! a *disclosure* change rather than a behaviour change: the normalized strings
//! are pinned to the bytes the pre-fix binary produced.

use std::path::Path;
use std::process::Command;

use ferro_hgvs::prepare::manifest::ReferenceManifest;

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// The single contig in the fixture: 400 bases cycling `ACGT`.
const CONTIG: &str = "NC_000001.11";
const CONTIG_LEN: usize = 400;

/// The base at 1-based position `pos` of the fixture contig.
fn base_at(pos: usize) -> char {
    ['A', 'C', 'G', 'T'][(pos - 1) % 4]
}

/// Write a minimal reference dir: an **indexed** genome FASTA with one contig,
/// plus a manifest pointing at it. Modelled on `issue_1181_cli_error_mode.rs`,
/// which needs the same thing — a real on-disk reference, because under the mock
/// provider no reference-base work runs and the diagnostics under test never
/// fire.
fn minimal_genome_reference(dir: &Path) {
    let seq: String = (1..=CONTIG_LEN).map(base_at).collect();
    let header = format!(">{CONTIG}\n");
    std::fs::write(dir.join("genome.fa"), format!("{header}{seq}\n")).unwrap();
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

/// `(exit ok, stdout, stderr)` from one `ferro normalize` invocation.
///
/// stdout and stderr are returned **separately**, deliberately: the whole point
/// of the fix is that the result stream stays exactly what it was while the
/// diagnostics go to stderr. Merging them would make the byte-identity
/// assertion below unable to see the difference.
fn normalize(reference: &str, variant: &str, extra: &[&str]) -> (bool, String, String) {
    let mut cmd = ferro();
    cmd.args(["normalize", "--reference", reference]);
    cmd.args(extra);
    cmd.arg(variant);
    let out = cmd.output().unwrap();
    (
        out.status.success(),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

/// Set up the fixture; returns (tempdir guard, reference path).
fn fixture() -> (tempfile::TempDir, String) {
    let dir = tempfile::tempdir().unwrap();
    minimal_genome_reference(dir.path());
    let reference = dir.path().to_str().unwrap().to_string();
    (dir, reference)
}

/// Two adjacent substitutions the normalizer merges into one `delins`. Bases at
/// 201/202 are `A`/`C`, so `[201A>T;202C>G]` states the reference correctly and
/// the only repair is the coalescing itself.
const COALESCING_INPUT: &str = "NC_000001.11:g.[201A>T;202C>G]";

/// A bracketed reference-range `ins` payload. Positions 120..=130 spell
/// `TACGTACGTAC` on the fixture contig.
const RANGE_INS_INPUT: &str = "NC_000001.11:g.100_101ins[120_130]";

/// A substitution stating `C` where the contig actually carries `A`.
const MISMATCHING_INPUT: &str = "NC_000001.11:g.201C>T";

/// The fixture bases the three inputs above assume. Asserted rather than
/// trusted, so a change to `CONTIG_LEN` or `base_at` fails here — where the
/// cause is obvious — instead of somewhere downstream as a mystery.
#[test]
fn the_fixture_carries_the_bases_the_inputs_assume() {
    assert_eq!(base_at(201), 'A');
    assert_eq!(base_at(202), 'C');
    let payload: String = (120..=130).map(base_at).collect();
    assert_eq!(payload, "TACGTACGTAC");
}

/// W5005 reaches the CLI. Before the fix this produced the same stdout line with
/// **no stderr diagnostic at all** — the coalescing was completely silent.
#[test]
fn the_coalesce_provenance_warning_reaches_the_cli() {
    let (_dir, reference) = fixture();
    let (ok, stdout, stderr) = normalize(&reference, COALESCING_INPUT, &[]);

    assert!(ok, "must normalize successfully; stderr: {stderr}");
    assert!(
        stdout.contains("NC_000001.11:g.[201A>T;202C>G] -> NC_000001.11:g.201_202delinsTG"),
        "precondition: the input must actually coalesce; stdout: {stdout}"
    );
    // The whole line, verbatim: a caller has to be able to select on the code
    // AND read the reason, so both are pinned.
    assert!(
        stderr.contains(
            "warning[MEMBERS_COALESCED_FROM_REPORTED_FORM]: NC_000001.11:g. — input described 2 \
             cis members, normalized form describes 1; the individually reported form is not \
             recoverable from the normalized string (DNA/delins.md:79-84)"
        ),
        "the coalesce must be disclosed with its code and its reason; stderr: {stderr}"
    );
}

/// `INSERTED_SEQUENCE_EXPANDED` reaches the CLI, and the message names both the
/// payload that was thrown away and the literal that replaced it — which is the
/// only place the discarded `[120_130]` provenance now survives.
#[test]
fn the_inserted_sequence_expansion_warning_reaches_the_cli() {
    let (_dir, reference) = fixture();
    let (ok, stdout, stderr) = normalize(&reference, RANGE_INS_INPUT, &[]);

    assert!(ok, "must normalize successfully; stderr: {stderr}");
    assert!(
        stdout
            .contains("NC_000001.11:g.100_101ins[120_130] -> NC_000001.11:g.100_101insTACGTACGTAC"),
        "precondition: the payload must actually be expanded; stdout: {stdout}"
    );
    assert!(
        stderr.contains(
            "warning[INSERTED_SEQUENCE_EXPANDED]: NC_000001.11: ins payload [120_130] expanded \
             to literal TACGTACGTAC"
        ),
        "the expansion must be disclosed, naming the original payload; stderr: {stderr}"
    );
}

/// `REFSEQ_MISMATCH` reaches the CLI under a non-strict error mode — the case
/// where the description is *accepted* despite contradicting the reference, and
/// so the one where silence is most costly.
///
/// The rendered position used to be left unpinned here. On this input the
/// message read `at 100-100` for a variant at `g.201` — a window-relative
/// offset leaking into a user-facing coordinate — and freezing it would have
/// made the eventual repair fail this test for the wrong reason.
///
/// That finding was filed as **#1612** and is now **fixed**, so the coordinate
/// is pinned: the warning is built from the fetched window's `start`/`end`, and
/// those are converted back through `MismatchFrame` before they are formatted.
/// `tests/it/issue_1612_refseq_mismatch_coordinate.rs` owns the full guard —
/// every axis, the strict-mode `FerroError::ReferenceMismatch` `location` (the
/// same string), and the agreement with the cis producer that the
/// `position`-keyed dedupe depends on. What is pinned *here* is only that the
/// repaired coordinate reaches the CLI, which is this file's own subject.
#[test]
fn the_reference_mismatch_warning_reaches_the_cli() {
    let (_dir, reference) = fixture();
    let (ok, _stdout, stderr) =
        normalize(&reference, MISMATCHING_INPUT, &["--error-mode", "lenient"]);

    assert!(ok, "lenient must accept the mismatch; stderr: {stderr}");
    assert!(
        stderr.contains("warning[REFSEQ_MISMATCH]: reference sequence mismatch at 201-201:"),
        "the mismatch must be disclosed under its code, at the coordinate the \
         input names (#1612); stderr: {stderr}"
    );
    assert!(
        stderr.contains(r#"stated "C", actual "A""#),
        "the disclosure must name both bases, or a caller cannot act on it; stderr: {stderr}"
    );
}

/// Strict mode is not a substitute for reading the warnings, which is the reason
/// this fix is needed at all rather than "just use `--error-mode strict`".
///
/// For the coalescing input, strict and lenient produce **byte-identical stdout
/// and byte-identical stderr**: the strict rejection ladder does not cover
/// W5005, so strict rejects nothing and reports exactly what lenient does. A
/// caller who picked strict specifically to be told about repairs was told
/// nothing.
#[test]
fn strict_mode_mitigates_nothing_for_a_coalesced_allele() {
    let (_dir, reference) = fixture();
    let (strict_ok, strict_out, strict_err) =
        normalize(&reference, COALESCING_INPUT, &["--error-mode", "strict"]);
    let (lenient_ok, lenient_out, lenient_err) =
        normalize(&reference, COALESCING_INPUT, &["--error-mode", "lenient"]);

    assert!(strict_ok && lenient_ok, "both modes accept this input");
    assert_eq!(
        strict_out, lenient_out,
        "strict does not change the normalized output for a coalesced allele"
    );
    assert_eq!(
        strict_err, lenient_err,
        "nor the diagnostics — strict promotes no rung that covers W5005"
    );
    assert!(
        strict_err.contains("warning[MEMBERS_COALESCED_FROM_REPORTED_FORM]"),
        "and both must disclose it; strict stderr: {strict_err}"
    );
}

/// `--format json` carries the warnings in the record, under their own key.
///
/// `corrections` is the preprocessor's population and stays exactly what it was;
/// `warnings` is the normalizer's, and is always present (`[]` included) so a
/// consumer can index it unconditionally — the same contract
/// `ferro project --format json` adopted in #1182.
#[test]
fn json_records_carry_the_warnings_under_their_own_key() {
    let (_dir, reference) = fixture();

    let (ok, stdout, stderr) = normalize(&reference, COALESCING_INPUT, &["--format", "json"]);
    assert!(ok, "must normalize successfully; stderr: {stderr}");
    assert!(
        stdout.contains(r#""corrections": []"#),
        "the preprocessor population must stay separate and stay empty here; stdout: {stdout}"
    );
    assert!(
        stdout.contains(r#""warnings": [{"code":"MEMBERS_COALESCED_FROM_REPORTED_FORM""#),
        "the normalizer warning must be in the JSON record; stdout: {stdout}"
    );

    // The key is unconditional, so a clean input still carries an empty array
    // rather than omitting it.
    let (clean_ok, clean_stdout, _) =
        normalize(&reference, "NC_000001.11:g.201A>T", &["--format", "json"]);
    assert!(clean_ok);
    assert!(
        clean_stdout.contains(r#""warnings": []"#),
        "a clean record must still carry the key; stdout: {clean_stdout}"
    );
}

/// `--format tsv` correlates the warning to the row that caused it, via the
/// `detail` column — the same column the preprocessor's corrections already use,
/// and for the same reason ("greppable where it always was and correlatable to
/// the row it explains").
#[test]
fn tsv_rows_carry_the_warning_in_the_detail_column() {
    let (_dir, reference) = fixture();
    let (ok, stdout, stderr) = normalize(&reference, COALESCING_INPUT, &["--format", "tsv"]);

    assert!(ok, "must normalize successfully; stderr: {stderr}");
    let row = stdout
        .lines()
        .find(|l| l.contains("NC_000001.11:g.[201A>T;202C>G]"))
        .unwrap_or_else(|| panic!("no data row in tsv output: {stdout}"));
    let detail = row
        .split('\t')
        .next_back()
        .expect("a tsv row always has a last field");
    assert!(
        detail.starts_with("MEMBERS_COALESCED_FROM_REPORTED_FORM: "),
        "the detail column must carry the warning; row: {row}"
    );
    // stderr keeps its copy too, so an existing `2>&1 | grep warning` pipeline
    // sees it in the shape it always saw preprocessor warnings in.
    assert!(
        stderr.contains("warning[MEMBERS_COALESCED_FROM_REPORTED_FORM]"),
        "tsv must also report on stderr, as the text format does; stderr: {stderr}"
    );
}

/// **The disclosure moved no string.**
///
/// Every normalized form below is the byte-for-byte stdout the *pre-fix* binary
/// produced for that input — captured by running `ferro normalize` at the
/// campaign base (`35de96c8`) against this same fixture, before any of the
/// changes in this test's PR existed. Routing the CLI from `normalize` to
/// `normalize_with_diagnostics` cannot move the result, because both exits go
/// through the same `normalize_core_checked` and return the same variant; this
/// pins that claim to observed bytes rather than to the argument for it.
///
/// It covers all three formats, because each renders the result through a
/// different path and only `text` shares code with the other two.
#[test]
fn the_disclosure_does_not_move_the_normalized_string() {
    let (_dir, reference) = fixture();

    // (input, exact stdout line the pre-fix binary emitted for `--format text`)
    let pinned: [(&str, &str); 3] = [
        (
            COALESCING_INPUT,
            "NC_000001.11:g.[201A>T;202C>G] -> NC_000001.11:g.201_202delinsTG",
        ),
        (
            RANGE_INS_INPUT,
            "NC_000001.11:g.100_101ins[120_130] -> NC_000001.11:g.100_101insTACGTACGTAC",
        ),
        // Unchanged: lenient accepts the contradicting base and preserves it.
        (MISMATCHING_INPUT, "NC_000001.11:g.201C>T"),
    ];

    for (input, expected) in pinned {
        let (ok, stdout, stderr) = normalize(&reference, input, &["--error-mode", "lenient"]);
        assert!(ok, "{input} must normalize; stderr: {stderr}");
        assert_eq!(
            stdout.trim_end(),
            expected,
            "the result stream moved for {input} — this PR is a disclosure change \
             and must not touch the normalized form"
        );

        // The same string, through the JSON renderer.
        let normalized = expected.rsplit(" -> ").next().unwrap();
        let (json_ok, json_stdout, _) = normalize(
            &reference,
            input,
            &["--error-mode", "lenient", "--format", "json"],
        );
        assert!(json_ok, "{input} must normalize under --format json");
        assert!(
            json_stdout.contains(&format!(r#""output": "{normalized}""#)),
            "json output moved for {input}: {json_stdout}"
        );

        // And through the TSV renderer, whose `normalized` column is field 3.
        let (tsv_ok, tsv_stdout, _) = normalize(
            &reference,
            input,
            &["--error-mode", "lenient", "--format", "tsv"],
        );
        assert!(tsv_ok, "{input} must normalize under --format tsv");
        let row = tsv_stdout
            .lines()
            .find(|l| l.starts_with(&format!("\t{input}\t")))
            .unwrap_or_else(|| panic!("no tsv row for {input}: {tsv_stdout}"));
        assert_eq!(
            row.split('\t').nth(2),
            Some(normalized),
            "tsv normalized column moved for {input}"
        );
    }
}
