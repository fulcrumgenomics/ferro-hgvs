//! `ferro hgvs-to-vcf` must report a declined description and keep going (#1764).
//!
//! The batch loop used to propagate a per-record conversion error out of the
//! whole command, so one description the converter declines discarded every
//! remaining line — leaving the caller a *truncated* VCF on stdout plus a
//! non-zero exit, with no indication of which input line was responsible.
//!
//! The convention these tests pin is not new. `run_normalize`, `run_parse` and
//! `run_project` in the same binary already all: continue past a per-record
//! failure, format the diagnostic to **stderr** through
//! `output_error_with_context` with the 1-based input line number, count the
//! failures, and end the run with `Err("<n> variant(s) failed to …")` so the
//! exit status is non-zero whenever any record was declined. `hgvs-to-vcf` was
//! the odd one out; these tests hold it to the same contract.
//!
//! Deliberately spawned-binary tests rather than in-process ones: the behaviour
//! under test *is* the command's contract — stream completeness, stream
//! separation and exit status — and none of those three is observable from a
//! library call.
//!
//! Every pinned string below was **measured** from a run of the fixed binary,
//! not inferred from reading the formatter.

use std::io::Write;
use std::process::{Command, Stdio};

use tempfile::NamedTempFile;

fn ferro() -> Command {
    Command::new(env!("CARGO_BIN_EXE_ferro"))
}

/// The three-line reproducer from #1764, plus a fourth row: a convertible row,
/// a row the converter declines (#1734 correctly refuses an offset on a genomic
/// position), then two more convertible rows.
const DECLINE_IN_THE_MIDDLE: &str = "NC_000001.11:g.10C>T\n\
                                     NC_000001.11:g.10+2del\n\
                                     NC_000001.11:g.12C>A\n\
                                     NC_000001.11:g.20A>G\n";

/// The exact VCF rows the three convertible descriptions above produce.
const ROW_10: &str = "chr1\t10\t.\tC\tT\t.\tPASS\t.";
const ROW_12: &str = "chr1\t12\t.\tC\tA\t.\tPASS\t.";
const ROW_20: &str = "chr1\t20\t.\tA\tG\t.\tPASS\t.";

/// The declined row's diagnostic, as far as the wording is this command's
/// contract: the `ERROR (line N): <input> - ` shape plus the start of the
/// converter's own refusal. The rest of the converter's message is #1734's to
/// word, so pinning it here would make this test fail on an unrelated reword.
const DECLINE_REPORT: &str =
    "ERROR (line 2): NC_000001.11:g.10+2del - Invalid coordinates: genomic position 10+2";

/// Write `lines` to a temp file and hand back the handle. The handle is
/// returned (not just its path) so the caller keeps it alive for the run — a
/// dropped `NamedTempFile` deletes the file.
fn input_file(lines: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().expect("create temp input");
    f.write_all(lines.as_bytes()).expect("write temp input");
    f.flush().expect("flush temp input");
    f
}

/// Run `hgvs-to-vcf -i <temp file> -f <format>`, returning
/// `(stdout, stderr, success)`.
fn run_on_file(lines: &str, format: &str) -> (String, String, bool) {
    let input = input_file(lines);
    let out = ferro()
        .args(["hgvs-to-vcf", "-i"])
        .arg(input.path())
        .args(["-f", format])
        .output()
        .expect("spawn ferro");
    (
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
        out.status.success(),
    )
}

/// The exact reproducer from #1764. All three convertible rows must reach
/// stdout — the two *after* the declined one are what the abort discarded.
#[test]
fn a_declined_description_does_not_truncate_the_stream() {
    let (stdout, stderr, ok) = run_on_file(DECLINE_IN_THE_MIDDLE, "vcf");

    let rows: Vec<&str> = stdout.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        vec![ROW_10, ROW_12, ROW_20],
        "every convertible row must be emitted, in input order;\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );

    // The declined row is reported on stderr, naming its 1-based input line and
    // the description — the `ERROR (line N): <input> - <error>` shape the
    // sibling drivers emit.
    assert!(
        stderr.contains(DECLINE_REPORT),
        "stderr must carry {DECLINE_REPORT:?}; got:\n{stderr}"
    );

    // A run in which some records were declined is a failure, per the sibling
    // convention (`run_normalize`/`run_parse`/`run_project` all end with
    // `Err("<n> variant(s) failed to …")`). Silence here would make a partial
    // conversion indistinguishable from a complete one.
    assert!(
        !ok,
        "a run with declined records must exit non-zero; stderr:\n{stderr}"
    );
    assert!(
        stderr.contains("1 variant(s) failed to convert"),
        "the run must say how many records it declined:\n{stderr}"
    );
}

/// A non-genomic description is a per-record decline like any other: it belongs
/// on stderr and it must make the run fail. Before #1764 it was written to
/// **stdout**, inside the VCF body — which corrupts a downstream parse — and
/// the run still exited 0, so the failure was invisible.
#[test]
fn a_non_genomic_description_is_reported_on_stderr_and_fails_the_run() {
    let (stdout, stderr, ok) = run_on_file(
        "NC_000001.11:g.10C>T\n\
         NM_000088.3:c.459A>G\n\
         NC_000001.11:g.12C>A\n",
        "vcf",
    );

    let rows: Vec<&str> = stdout.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        vec![ROW_10, ROW_12],
        "both genomic rows must convert, and nothing else may enter the VCF body;\nstdout:\n{stdout}"
    );

    assert!(
        stderr.contains(
            "ERROR (line 2): NM_000088.3:c.459A>G - Unsupported variant type: \
             only genomic (g.) descriptions convert to VCF"
        ),
        "the non-genomic description must be reported on stderr with its line:\n{stderr}"
    );
    assert!(
        !ok,
        "a non-genomic record is a failure, not a silent skip; stderr:\n{stderr}"
    );
}

/// The failure counter must accumulate, and every declined line must be named.
///
/// Every other test in this module declines exactly **one** record, so all of
/// them pass against a run that hardcoded `1 variant(s) failed to convert` — a
/// `error_count = 1` where the code means `error_count += 1` is invisible to
/// them. Two declines, of two different kinds (a refused genomic offset and a
/// non-genomic description), on non-adjacent lines, so the count, the plural
/// wording and the per-line attribution are all separated from the single-
/// failure case.
#[test]
fn two_declines_are_counted_and_reported_line_by_line() {
    let (stdout, stderr, ok) = run_on_file(
        "NC_000001.11:g.10C>T\n\
         NC_000001.11:g.10+2del\n\
         NC_000001.11:g.12C>A\n\
         NM_000088.3:c.459A>G\n\
         NC_000001.11:g.20A>G\n",
        "vcf",
    );

    let rows: Vec<&str> = stdout.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        vec![ROW_10, ROW_12, ROW_20],
        "both declines must be stepped over, not one;\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );

    assert!(
        stderr.contains(DECLINE_REPORT),
        "the line-2 decline must be reported against line 2:\n{stderr}"
    );
    assert!(
        stderr.contains("ERROR (line 4): NM_000088.3:c.459A>G - "),
        "the line-4 decline must be reported against line 4:\n{stderr}"
    );
    assert!(
        stderr.contains("2 variant(s) failed to convert"),
        "the terminal report must count both declines, not just the last:\n{stderr}"
    );
    assert!(!ok, "stderr:\n{stderr}");
}

/// The all-convertible case must stay exactly as it was: every row on stdout,
/// nothing on stderr, exit 0. This is the discriminating test against "report
/// everything as a failure" — an exit status that is always non-zero separates
/// nothing, which is precisely what the CLI contract asks it to separate.
#[test]
fn a_fully_convertible_run_still_exits_zero_with_a_quiet_stderr() {
    let (stdout, stderr, ok) = run_on_file("NC_000001.11:g.10C>T\nNC_000001.11:g.12C>A\n", "vcf");

    assert!(ok, "an all-convertible run must exit 0; stderr:\n{stderr}");
    let rows: Vec<&str> = stdout.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(rows, vec![ROW_10, ROW_12], "stdout:\n{stdout}");
    assert!(
        stderr.is_empty(),
        "an all-convertible run must not report anything:\n{stderr}"
    );
}

/// `--format json` keeps the result stream parseable: the decline is a JSON
/// object on **stderr** carrying the `line`, and the success objects on stdout
/// are uninterrupted. This is the format-awareness the sibling drivers get from
/// routing declines through `output_error_with_context` rather than a bare
/// `writeln!`.
#[test]
fn json_format_reports_the_decline_as_json_on_stderr() {
    let (stdout, stderr, ok) = run_on_file(DECLINE_IN_THE_MIDDLE, "json");

    let rows: Vec<&str> = stdout.lines().filter(|l| !l.trim().is_empty()).collect();
    assert_eq!(
        rows,
        vec![
            r#"{"hgvs": "NC_000001.11:g.10C>T", "chrom": "chr1", "pos": 10, "ref": "C", "alt": "T"}"#,
            r#"{"hgvs": "NC_000001.11:g.12C>A", "chrom": "chr1", "pos": 12, "ref": "C", "alt": "A"}"#,
            r#"{"hgvs": "NC_000001.11:g.20A>G", "chrom": "chr1", "pos": 20, "ref": "A", "alt": "G"}"#,
        ],
        "the JSON result stream must carry every convertible row and nothing else;\nstdout:\n{stdout}"
    );

    assert!(
        stderr.contains(r#""input": "NC_000001.11:g.10+2del""#)
            && stderr.contains(r#""line": 2"#)
            && stderr.contains(r#""status": "error""#),
        "the decline must be a JSON status:error object naming its input line:\n{stderr}"
    );
    assert!(!ok, "stderr:\n{stderr}");
}

/// The stdin path is a third copy of the batch loop, alongside `--input` and the
/// single-variant path. Pinned separately because a fix applied to only one of
/// the three would leave the spellings disagreeing — and stdin is how the
/// command is used in a pipeline, which is exactly where a truncated stream
/// does the most damage.
#[test]
fn the_stdin_path_continues_past_a_decline_too() {
    let mut child = ferro()
        .args(["hgvs-to-vcf", "-f", "vcf"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn ferro");
    child
        .stdin
        .take()
        .expect("stdin piped")
        .write_all(DECLINE_IN_THE_MIDDLE.as_bytes())
        .expect("write stdin");
    let out = child.wait_with_output().expect("wait for ferro");

    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);

    let rows: Vec<&str> = stdout.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        rows,
        vec![ROW_10, ROW_12, ROW_20],
        "the stdin path must emit every convertible row;\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    assert!(
        stderr.contains(DECLINE_REPORT),
        "the stdin path must report the decline with its line number:\n{stderr}"
    );
    assert!(!out.status.success(), "stderr:\n{stderr}");
}

/// The single-variant path (`ferro hgvs-to-vcf <variant>`) declines the same
/// way, minus the line number — there is no file to number lines in. Pinned
/// because `run_normalize` and `run_project` both special-case that path.
#[test]
fn the_single_variant_path_reports_on_stderr_without_a_line_number() {
    let out = ferro()
        .args(["hgvs-to-vcf", "NC_000001.11:g.10+2del", "-f", "vcf"])
        .output()
        .expect("spawn ferro");

    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert!(!out.status.success(), "stderr:\n{stderr}");
    assert!(
        stderr
            .contains("ERROR: NC_000001.11:g.10+2del - Invalid coordinates: genomic position 10+2"),
        "stderr must name the description, with no line number:\n{stderr}"
    );
    assert!(
        !stderr.contains("line 1"),
        "there is no input line to name on the single-variant path:\n{stderr}"
    );
    // Only the header reached stdout — no data row, and no diagnostic.
    assert_eq!(
        stdout.lines().filter(|l| !l.starts_with('#')).count(),
        0,
        "nothing but the header may reach stdout:\n{stdout}"
    );
}
