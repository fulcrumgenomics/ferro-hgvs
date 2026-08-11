//! Integration tests for `ferro normalize --format tsv`.
//!
//! The format exists to answer "which of these variants does ferro change, and
//! to what?" over a list, so the properties worth pinning are: a header is
//! always present, every input yields exactly one row (a failure included), the
//! `changed` column agrees with both the `input`/`normalized` columns and the
//! stderr summary, each row names the input line it came from, and nothing but
//! rows reaches the output stream.
//!
//! All tests run under the mock provider (no `--reference`), so they need no
//! external data. `--error-mode silent` is used where a preprocessing
//! correction is wanted as the source of a change: under the default strict mode
//! those same inputs are rejected instead, which the strict-mode test covers.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// The header the table must always start with.
const HEADER: &str = "line\tinput\tnormalized\tchanged\tstatus\tdetail";

/// Number of columns in the table.
const COLUMNS: usize = 6;

/// A variant the mock provider normalizes back to itself.
const UNCHANGED: &str = "NM_000088.3:c.459A>G";
/// A variant whose redundant deleted sequence is dropped by the *preprocessor*,
/// so it normalizes to a different string.
const CHANGED_IN: &str = "NM_000088.3:c.459delA";
const CHANGED_OUT: &str = "NM_000088.3:c.459del";
/// A deletion inside the `CCCC` homopolymer at c.16-19 of the mock `NM_000088.3`
/// sequence: the *normalizer* 3'-shifts it, so this is a change produced by the
/// feature's headline code path rather than by a preprocessing rewrite.
const SHIFTED_IN: &str = "NM_000088.3:c.16del";
const SHIFTED_OUT: &str = "NM_000088.3:c.19del";

/// Mirror of `BATCH_CHUNK_LINES` in `src/bin/ferro.rs`. Integration tests cannot
/// import items from a `bin` target, so this is duplicated here; the runtime
/// assertion in the parallel test verifies the chosen input size still spans
/// multiple chunks, so drift in the real constant fails loudly.
const BATCH_CHUNK_LINES: usize = 8192;

/// One completed `ferro normalize` run.
struct Run {
    stdout: String,
    stderr: String,
    success: bool,
}

impl Run {
    /// The output stream split into lines, with the header line removed after
    /// asserting it is present, first, and the only one.
    fn rows(&self) -> Vec<&str> {
        let mut lines = self.stdout.lines();
        assert_eq!(
            lines.next(),
            Some(HEADER),
            "first output line must be the TSV header; got:\n{}",
            self.stdout
        );
        let rows: Vec<&str> = lines.collect();
        assert!(
            !rows.contains(&HEADER),
            "the header must appear exactly once:\n{}",
            self.stdout
        );
        rows
    }

    /// The single row for a one-variant run, split into its fields.
    fn only_row_fields(&self) -> Vec<&str> {
        let rows = self.rows();
        assert_eq!(rows.len(), 1, "expected exactly one row, got {rows:?}");
        fields(rows[0])
    }

    /// The `summary:` line from stderr.
    fn summary(&self) -> &str {
        self.stderr
            .lines()
            .find(|l| l.starts_with("summary:"))
            .unwrap_or_else(|| panic!("no summary line on stderr:\n{}", self.stderr))
    }
}

/// Split one row into its fields, asserting the column count.
fn fields(row: &str) -> Vec<&str> {
    let f: Vec<&str> = row.split('\t').collect();
    assert_eq!(
        f.len(),
        COLUMNS,
        "row does not have {COLUMNS} fields: {row:?}"
    );
    f
}

/// Run `ferro normalize` with the given extra arguments.
fn run(args: &[&str]) -> Run {
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .arg("normalize")
        .args(args)
        .output()
        .expect("run ferro normalize");
    Run {
        stdout: String::from_utf8(out.stdout).expect("stdout is UTF-8"),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        success: out.status.success(),
    }
}

/// Write `variants` to a temporary input file, one per line. The handle is
/// returned so the file outlives the run.
fn input_file(variants: &[&str]) -> NamedTempFile {
    let mut tf = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    for v in variants {
        writeln!(tf, "{v}").unwrap();
    }
    tf.flush().unwrap();
    tf
}

/// Run `ferro normalize --format tsv` over an input file of `variants`.
fn run_tsv(variants: &[&str]) -> Run {
    let tf = input_file(variants);
    run(&[
        "--error-mode",
        "silent",
        "-f",
        "tsv",
        "-i",
        tf.path().to_str().unwrap(),
    ])
}

#[test]
fn emits_a_header_even_for_empty_input() {
    // A consumer reading the table with a header-aware parser must not have to
    // special-case "no variants matched".
    let run = run_tsv(&[]);
    assert!(run.success, "empty input should succeed: {}", run.stderr);
    assert!(run.rows().is_empty());
}

#[test]
fn reports_an_unchanged_variant_as_unchanged() {
    let run = run_tsv(&[UNCHANGED]);
    assert!(run.success, "{}", run.stderr);
    assert_eq!(
        run.only_row_fields(),
        vec!["1", UNCHANGED, UNCHANGED, "false", "ok", ""]
    );
}

#[test]
fn reports_a_changed_variant_with_its_new_form() {
    let run = run_tsv(&[CHANGED_IN]);
    assert!(run.success, "{}", run.stderr);
    assert_eq!(
        run.only_row_fields(),
        vec!["1", CHANGED_IN, CHANGED_OUT, "true", "ok", ""]
    );
}

#[test]
fn reports_a_change_produced_by_the_normalizer_itself() {
    // The headline case: no preprocessing rewrite is involved, the normalizer
    // 3'-shifts the deletion through a homopolymer. Runs under the default
    // (strict) error mode precisely to show the change is not a correction.
    let tf = input_file(&[SHIFTED_IN]);
    let run = run(&["-f", "tsv", "-i", tf.path().to_str().unwrap()]);
    assert!(run.success, "{}", run.stderr);
    assert_eq!(
        run.only_row_fields(),
        vec![
            "1",
            SHIFTED_IN,
            SHIFTED_OUT,
            "true",
            "ok",
            // No correction was applied, so nothing explains the change but the
            // normalization itself.
            ""
        ]
    );
    assert!(
        run.summary().contains("changed=1"),
        "a normalizer-driven change must be counted: {}",
        run.summary()
    );
}

#[test]
fn a_malformed_variant_becomes_a_row_and_does_not_abort_the_run() {
    // The whole point of the format: one bad line must not cost the caller the
    // rows for the good lines that follow it.
    let run = run_tsv(&[UNCHANGED, "not-a-variant", CHANGED_IN]);
    let rows = run.rows();
    assert_eq!(rows.len(), 3, "every input must yield a row: {rows:?}");

    let bad = fields(rows[1]);
    assert_eq!(bad[0], "2");
    assert_eq!(bad[1], "not-a-variant");
    assert_eq!(bad[2], "", "a failed variant has no normalized form");
    assert_eq!(bad[3], "", "a failed variant has no changed verdict");
    assert_eq!(bad[4], "parse_error");
    assert!(!bad[5].is_empty(), "a failure must carry a detail message");

    // The rows around it are the same as they would be on their own.
    assert_eq!(fields(rows[0])[4], "ok");
    assert_eq!(
        fields(rows[2]),
        vec!["3", CHANGED_IN, CHANGED_OUT, "true", "ok", ""]
    );

    // Exit code behavior is unchanged: any failed variant still fails the run.
    assert!(!run.success, "a failed variant must still exit non-zero");
}

#[test]
fn the_changed_column_agrees_with_the_summary() {
    // The summary excludes failures from its `unchanged` count, so the column has
    // to leave `changed` *empty* on a failure — otherwise `changed=false` rows and
    // the summary's `unchanged` count are two different numbers for the same
    // question.
    let run = run_tsv(&[UNCHANGED, CHANGED_IN, "not-a-variant", UNCHANGED]);
    let rows = run.rows();
    let count = |want: &str| rows.iter().filter(|r| fields(r)[3] == want).count();
    let (changed, unchanged, failed) = (count("true"), count("false"), count(""));
    assert_eq!((changed, unchanged, failed), (1, 2, 1), "rows: {rows:?}");
    assert_eq!(
        run.summary(),
        format!(
            "summary: total={} changed={changed} unchanged={unchanged} failed={failed}",
            rows.len()
        )
    );
}

#[test]
fn distinguishes_the_phase_a_variant_failed_in() {
    // Under the default strict mode a redundant deleted sequence is rejected by
    // the preprocessor rather than corrected, and an out-of-range position gets
    // past the parser and fails in the normalizer. The tokens must tell those
    // apart from a grammar failure.
    let tf = input_file(&[CHANGED_IN, "NM_000088.3:c.99999A>G", "not-a-variant"]);
    let run = run(&[
        "--error-mode",
        "strict",
        "-f",
        "tsv",
        "-i",
        tf.path().to_str().unwrap(),
    ]);
    let rows = run.rows();
    assert_eq!(rows.len(), 3, "{rows:?}");
    assert_eq!(fields(rows[0])[4], "preprocess_error");
    assert_eq!(fields(rows[1])[4], "normalize_error");
    assert_eq!(fields(rows[2])[4], "parse_error");
    for row in &rows {
        let f = fields(row);
        assert_eq!(f[2], "", "a failed variant has no normalized form: {row:?}");
        assert_eq!(f[3], "", "a failed variant has no changed verdict: {row:?}");
        assert!(!f[5].is_empty(), "a failure must carry a detail: {row:?}");
    }
}

#[test]
fn a_successful_row_explains_a_correction_in_its_detail() {
    // Under lenient mode the preprocessor rewrites the input and the row is
    // `changed` because of it. The reason has to be *in the row*: the stderr
    // warning carries no line or row identifier, so it cannot be correlated back.
    let tf = input_file(&[UNCHANGED, CHANGED_IN]);
    let run = run(&[
        "--error-mode",
        "lenient",
        "-f",
        "tsv",
        "-i",
        tf.path().to_str().unwrap(),
    ]);
    assert!(run.success, "{}", run.stderr);
    let rows = run.rows();
    // Row 0 is uncorrected *by the preprocessor*, which is the contrast this
    // test draws — so what must be absent from its detail is the W3025 code,
    // not all content.
    //
    // Its detail is no longer empty. `UNCHANGED` is `c.459` on the mock
    // `NM_000088.3`, whose CDS is 60 bases, so the **normalizer** raises W4004
    // for it — and always did. Until `ferro normalize` moved off
    // `Normalizer::normalize` (the exit that discards warnings) that diagnostic
    // was produced and dropped, which is why this row read as clean. It is
    // pinned rather than dodged with a different fixture: an out-of-CDS position
    // silently accepted is exactly the thing worth showing.
    let uncorrected = fields(rows[0]);
    assert!(
        !uncorrected[5].contains("W3025"),
        "row 0 must carry no preprocessor correction: {uncorrected:?}"
    );
    assert_eq!(
        uncorrected[5], "POSITION_PAST_END: NM_000088.3:c.459 lies past the cds-end (bound 60)",
        "the normalizer's own diagnostic belongs in the row that caused it"
    );

    let corrected = fields(rows[1]);
    assert_eq!(corrected[3], "true");
    assert_eq!(corrected[4], "ok");
    assert!(
        corrected[5].starts_with("W3025: "),
        "a corrected row must name the correction code and message: {corrected:?}"
    );

    // The warning also still goes to stderr, exactly as it did before.
    assert!(
        run.stderr.contains("warning[W3025]: "),
        "the stderr warning must be kept: {}",
        run.stderr
    );
}

#[test]
fn line_numbers_survive_skipped_lines() {
    // The row index is not a substitute for the line number: blank and comment
    // lines are skipped, so the two diverge from the first skip onwards.
    let mut tf = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    writeln!(tf, "# a comment").unwrap();
    writeln!(tf).unwrap();
    writeln!(tf, "{UNCHANGED}").unwrap();
    writeln!(tf).unwrap();
    writeln!(tf, "not-a-variant").unwrap();
    tf.flush().unwrap();

    let run = run(&[
        "--error-mode",
        "silent",
        "-f",
        "tsv",
        "-i",
        tf.path().to_str().unwrap(),
    ]);
    let rows = run.rows();
    assert_eq!(rows.len(), 2, "{rows:?}");
    assert_eq!(fields(rows[0])[0], "3");
    assert_eq!(fields(rows[1])[0], "5");
}

#[test]
fn a_failure_detail_cannot_break_the_table() {
    // A parse error message quotes the offending input, so an input containing
    // control characters is the CLI-level route to a would-be multi-line detail
    // field. Every physical line of the output must still be one full row, and no
    // field may carry anything a line-splitter would break on — including the
    // separators `str::lines()` ignores but Python's `splitlines()` honours.
    let exotic = "NM_000088.3:c.459A>G\tx\u{000B}y\u{2028}z\u{0085}w\u{001B}[31m";
    let tf = input_file(&[exotic]);
    let run = run(&[
        "--error-mode",
        "silent",
        "-f",
        "tsv",
        "-i",
        tf.path().to_str().unwrap(),
    ]);
    let rows = run.rows();
    assert_eq!(rows.len(), 1, "the input must stay one row: {rows:?}");
    for row in &rows {
        for f in fields(row) {
            assert!(
                !f.chars()
                    .any(|c| c.is_control() || matches!(c, '\u{0085}' | '\u{2028}' | '\u{2029}')),
                "no field may carry a row-breaking character: {row:?}"
            );
        }
    }
}

#[test]
fn the_summary_goes_to_stderr_and_never_to_the_table() {
    let run = run_tsv(&[UNCHANGED, CHANGED_IN, "not-a-variant"]);
    assert_eq!(
        run.summary(),
        "summary: total=3 changed=1 unchanged=1 failed=1"
    );
    assert!(
        !run.stdout.contains("summary:"),
        "the summary must not pollute the table:\n{}",
        run.stdout
    );
    // Every output line is a row (plus the header), so a naive `tail -n +2`
    // consumer never sees a non-row line.
    for row in &run.rows() {
        fields(row);
    }
}

#[test]
fn honors_the_output_flag() {
    let tf = input_file(&[UNCHANGED, CHANGED_IN]);
    let out = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
    let run = run(&[
        "--error-mode",
        "silent",
        "-f",
        "tsv",
        "-i",
        tf.path().to_str().unwrap(),
        "-o",
        out.path().to_str().unwrap(),
    ]);
    assert!(run.success, "{}", run.stderr);
    assert!(
        run.stdout.is_empty(),
        "with --output nothing should reach stdout: {}",
        run.stdout
    );
    let written = std::fs::read_to_string(out.path()).unwrap();
    let lines: Vec<&str> = written.lines().collect();
    assert_eq!(lines.len(), 3, "header plus two rows: {lines:?}");
    assert_eq!(lines[0], HEADER);
    assert_eq!(fields(lines[1])[3], "false");
    assert_eq!(fields(lines[2])[3], "true");
}

#[test]
fn a_missing_input_file_writes_no_table() {
    // Writing the header while opening the output stream made this produce a
    // header-only artifact byte-identical to a *successful* empty-input run, with
    // no summary to tell the two apart.
    let out = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
    let missing = out.path().with_extension("does-not-exist.txt");
    assert!(!missing.exists());
    let run = run(&[
        "--error-mode",
        "silent",
        "-f",
        "tsv",
        "-i",
        missing.to_str().unwrap(),
        "-o",
        out.path().to_str().unwrap(),
    ]);
    assert!(!run.success, "a missing input file must fail the run");
    assert_eq!(
        std::fs::read_to_string(out.path()).unwrap(),
        "",
        "a failed open must leave no table behind"
    );
}

#[test]
fn works_for_a_single_positional_variant() {
    // `--format tsv` must not be silently ignored outside `--input` mode. The
    // `line` column is empty there: there is no file whose lines to number.
    let ok = run(&["--error-mode", "silent", "-f", "tsv", CHANGED_IN]);
    assert!(ok.success, "{}", ok.stderr);
    assert_eq!(
        ok.only_row_fields(),
        vec!["", CHANGED_IN, CHANGED_OUT, "true", "ok", ""]
    );

    let bad = run(&["--error-mode", "silent", "-f", "tsv", "not-a-variant"]);
    let f = bad.only_row_fields();
    assert_eq!(f[0], "");
    assert_eq!(f[1], "not-a-variant");
    assert_eq!(f[4], "parse_error");
    assert!(!bad.success, "a failed variant must still exit non-zero");
}

#[test]
fn a_single_positional_variant_honors_the_output_flag() {
    let out = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
    let run = run(&[
        "--error-mode",
        "silent",
        "-f",
        "tsv",
        "-o",
        out.path().to_str().unwrap(),
        CHANGED_IN,
    ]);
    assert!(run.success, "{}", run.stderr);
    assert!(
        run.stdout.is_empty(),
        "with --output nothing should reach stdout: {}",
        run.stdout
    );
    let written = std::fs::read_to_string(out.path()).unwrap();
    assert_eq!(
        written,
        format!("{HEADER}\n\t{CHANGED_IN}\t{CHANGED_OUT}\ttrue\tok\t\n")
    );
}

#[test]
fn parallel_tsv_output_is_byte_identical_to_serial() {
    // The batch driver may fan out across a rayon pool, so the table must be
    // invariant under worker count: one header, the same bytes, input order
    // preserved across chunk boundaries.
    let n = 20_000;
    assert!(
        n > 2 * BATCH_CHUNK_LINES,
        "n={n} no longer spans multiple chunks (BATCH_CHUNK_LINES={BATCH_CHUNK_LINES}); \
         raise n so ordering is exercised across chunk boundaries"
    );
    let mut tf = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    for i in 1..=n {
        // Comments and blank lines are skipped, so the `line` column and the row
        // index diverge — which is exactly what makes the column worth checking.
        if i % 100 == 0 {
            writeln!(tf, "# comment at {i}").unwrap();
        }
        if i % 50 == 0 {
            writeln!(tf).unwrap();
        }
        writeln!(tf, "NM_000088.3:c.{i}A>G").unwrap();
    }
    tf.flush().unwrap();
    let path = tf.path().to_str().unwrap();

    let serial = run(&["--error-mode", "silent", "-f", "tsv", "-j", "1", "-i", path]);
    assert!(serial.success, "{}", serial.stderr);
    let rows = serial.rows();
    assert_eq!(rows.len(), n, "every variant must yield a row");

    // Input order preserved: the inputs come back in the order written, and the
    // line numbers are strictly increasing.
    let mut last_line = 0usize;
    for (i, row) in rows.iter().enumerate() {
        let f = fields(row);
        assert_eq!(
            f[1],
            format!("NM_000088.3:c.{}A>G", i + 1),
            "row {i} out of order"
        );
        let line: usize = f[0].parse().expect("line column is a number");
        assert!(line > last_line, "line numbers must increase: {row:?}");
        last_line = line;
    }

    for workers in [2usize, 4, 8] {
        let parallel = run(&[
            "--error-mode",
            "silent",
            "-f",
            "tsv",
            "-j",
            &workers.to_string(),
            "-i",
            path,
        ]);
        assert!(parallel.success, "{}", parallel.stderr);
        assert_eq!(
            serial.stdout, parallel.stdout,
            "the table for -j{workers} differs from serial (-j1)"
        );
        assert_eq!(
            serial.summary(),
            parallel.summary(),
            "the summary for -j{workers} differs from serial (-j1)"
        );
    }
}

#[test]
fn text_and_json_output_is_untouched() {
    // The `tsv` addition must be purely additive: these are the exact byte
    // streams the two pre-existing formats produced before it existed.
    //
    // The json record has since gained one key — `warnings`, the normalizer's
    // population, separate from the preprocessor's `corrections` — when
    // `ferro normalize` was moved onto the diagnostics exit. It is present
    // unconditionally so a consumer can index it without branching, which is why
    // it shows here as `[]`: under `--error-mode silent` W4004 is suppressed, so
    // these two rows raise nothing. (Under `lenient` the first row carries
    // `POSITION_PAST_END` — see
    // `a_successful_row_explains_a_correction_in_its_detail`.)
    let variants = [UNCHANGED, CHANGED_IN, "not-a-variant"];
    let tf = input_file(&variants);
    let path = tf.path().to_str().unwrap();

    let text = run(&["--error-mode", "silent", "-f", "text", "-i", path]);
    assert_eq!(
        text.stdout,
        format!("{UNCHANGED}\n{CHANGED_IN} -> {CHANGED_OUT}\n")
    );
    assert!(
        !text.stderr.contains("summary:"),
        "the TSV summary must not appear in text mode:\n{}",
        text.stderr
    );

    let json = run(&["--error-mode", "silent", "-f", "json", "-i", path]);
    assert_eq!(
        json.stdout,
        format!(
            "{{\"input\": \"{UNCHANGED}\", \"output\": \"{UNCHANGED}\", \
             \"status\": \"ok\", \"corrections\": [], \"warnings\": []}}\n\
             {{\"input\": \"{CHANGED_IN}\", \"output\": \"{CHANGED_OUT}\", \
             \"status\": \"ok\", \"corrections\": [], \"warnings\": []}}\n"
        )
    );
    assert!(
        !json.stderr.contains("summary:"),
        "the TSV summary must not appear in json mode:\n{}",
        json.stderr
    );

    // The default format is still `text`, not the new one.
    let default = run(&["--error-mode", "silent", "-i", path]);
    assert_eq!(default.stdout, text.stdout);
}
