//! Output formatting utilities for CLI operations
//!
//! # `ferro normalize --format tsv`
//!
//! The TSV writer ([`NORMALIZE_TSV_HEADER`], [`normalize_tsv_row`],
//! [`normalize_tsv_summary`]) renders a *report*: one row per input variant,
//! failures included, so the table answers "which of my variants changed, and to
//! what?" over a whole file. Its columns are:
//!
//! | column | meaning |
//! |---|---|
//! | `line` | 1-based line number in the `--input` file; empty for a single positional variant, which has no line context |
//! | `input` | the variant exactly as read, sanitized (see below) |
//! | `normalized` | the normalized HGVS string; empty for a failure |
//! | `changed` | `true`/`false` for a success; **empty** for a failure, which has no normalized string to compare against — this is what keeps the column consistent with [`normalize_tsv_summary`], whose `unchanged` count also excludes failures |
//! | `status` | `ok`, or a [`NormalizeTsvFailure`] token |
//! | `detail` | the failure message, or — on an `ok` row — the preprocessing corrections that were applied, `; `-joined |
//!
//! Two properties a consumer can rely on:
//!
//! - **`status` names the stage that rejected the input under the active error
//!   mode**, not an intrinsic property of the input. See
//!   [`NormalizeTsvFailure`].
//! - **Every physical line of the output is either the header or one row of
//!   exactly six tab-separated fields.** Fields are sanitized (see
//!   [`normalize_tsv_row`]) so no control character can split or shift a row.
//!   Skip exactly the first line to reach the rows; an input line that happens
//!   to be the literal text `input` cannot be mistaken for a second header,
//!   because its row is preceded by the `line` column.

use crate::error::FerroError;
use std::io::{self, Write};
use std::str::FromStr;

/// Output format for CLI results
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    /// Plain text format (default)
    #[default]
    Text,
    /// JSON format
    Json,
    /// VCF format
    Vcf,
}

impl FromStr for OutputFormat {
    type Err = std::convert::Infallible;

    /// Parse an output format from a string
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::cli::OutputFormat;
    /// use std::str::FromStr;
    ///
    /// assert!(matches!(OutputFormat::from_str("json").unwrap(), OutputFormat::Json));
    /// assert!(matches!(OutputFormat::from_str("text").unwrap(), OutputFormat::Text));
    /// assert!(matches!(OutputFormat::from_str("vcf").unwrap(), OutputFormat::Vcf));
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s.to_lowercase().as_str() {
            "json" => OutputFormat::Json,
            "vcf" => OutputFormat::Vcf,
            _ => OutputFormat::Text,
        })
    }
}

/// Write a successful result to the output
///
/// # Arguments
///
/// * `writer` - The output writer (can be stdout, file, or buffer for testing)
/// * `input` - The original input string
/// * `output` - The processed output string
/// * `format` - The output format
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::{output_result, OutputFormat};
/// use std::io::Cursor;
///
/// let mut buffer = Cursor::new(Vec::new());
/// output_result(&mut buffer, "input", "output", OutputFormat::Text).unwrap();
/// let result = String::from_utf8(buffer.into_inner()).unwrap();
/// assert!(result.contains("input -> output"));
/// ```
pub fn output_result<W: Write>(
    writer: &mut W,
    input: &str,
    output: &str,
    format: OutputFormat,
) -> io::Result<()> {
    match format {
        OutputFormat::Json => {
            writeln!(
                writer,
                r#"{{"input": "{}", "output": "{}", "status": "ok"}}"#,
                escape_json(input),
                escape_json(output)
            )
        }
        OutputFormat::Text | OutputFormat::Vcf => {
            if input == output {
                writeln!(writer, "{}", output)
            } else {
                writeln!(writer, "{} -> {}", input, output)
            }
        }
    }
}

/// Write one `ferro project` record.
///
/// - text: a rendered axis is the bare HGVS string (pipe-able); an unavailable
///   axis is written as a `#`-prefixed comment line on the same stream so it
///   never pollutes the HGVS output while remaining visible to humans.
/// - json: one object per record carrying axis/transcript/output/status.
///
/// # Arguments
///
/// * `writer` - The output writer
/// * `input` - The original input string
/// * `axis` - The selected output axis
/// * `outcome` - The axis selection outcome (rendered or unavailable)
/// * `format` - The output format
/// * `line_number` - Optional line number in the input file
pub fn output_projection(
    writer: &mut dyn Write,
    input: &str,
    axis: crate::cli::project::Axis,
    outcome: &crate::cli::project::AxisOutcome,
    format: OutputFormat,
    line_number: Option<usize>,
) -> io::Result<()> {
    use crate::cli::project::AxisOutcome;
    match (format, outcome) {
        (
            OutputFormat::Json,
            AxisOutcome::Rendered {
                transcript_id,
                output,
                warnings,
            },
        ) => writeln!(
            writer,
            r#"{{"input": "{}", "axis": "{}", "transcript": "{}", "output": "{}", "status": "ok", "warnings": {}}}"#,
            escape_json(input),
            axis.code(),
            escape_json(transcript_id),
            escape_json(output),
            warnings_json(warnings)
        ),
        (
            OutputFormat::Json,
            AxisOutcome::Unavailable {
                transcript_id,
                reason,
                warnings,
            },
        ) => writeln!(
            writer,
            r#"{{"input": "{}", "axis": "{}", "transcript": {}, "output": null, "status": "unavailable", "reason": "{}", "warnings": {}}}"#,
            escape_json(input),
            axis.code(),
            match transcript_id {
                Some(t) => format!(r#""{}""#, escape_json(t)),
                None => "null".to_string(),
            },
            escape_json(reason),
            warnings_json(warnings)
        ),
        (_, AxisOutcome::Rendered { output, .. }) => writeln!(writer, "{}", output),
        (_, AxisOutcome::Unavailable { reason, .. }) => match line_number {
            Some(n) => writeln!(
                writer,
                "# {} (line {}) [{}]: unavailable: {}",
                input,
                n,
                axis.code(),
                reason
            ),
            None => writeln!(
                writer,
                "# {} [{}]: unavailable: {}",
                input,
                axis.code(),
                reason
            ),
        },
    }
}

/// Write a `ferro project` hard-failure record (the exit-1 class: parse error,
/// transcript not found, `--transcript` mismatch, ambiguous bare-g, IO) to
/// `writer`, honoring the output format so `--format json` stays parseable.
///
/// Mirrors [`output_error_with_context`]'s JSON shape but takes a plain message
/// string, because a `ferro project` hard failure is carried as a
/// `ProjectCliError`, not a `FerroError`.
pub fn output_project_error(
    writer: &mut dyn Write,
    input: &str,
    message: &str,
    format: OutputFormat,
    line_number: Option<usize>,
) -> io::Result<()> {
    match format {
        OutputFormat::Json => match line_number {
            Some(n) => writeln!(
                writer,
                r#"{{"input": "{}", "error": "{}", "line": {}, "status": "error"}}"#,
                escape_json(input),
                escape_json(message),
                n
            ),
            None => writeln!(
                writer,
                r#"{{"input": "{}", "error": "{}", "status": "error"}}"#,
                escape_json(input),
                escape_json(message)
            ),
        },
        OutputFormat::Text | OutputFormat::Vcf => match line_number {
            Some(n) => writeln!(writer, "ERROR (line {}): {} - {}", n, input, message),
            None => writeln!(writer, "ERROR: {} - {}", input, message),
        },
    }
}

/// Write an error to the output
///
/// # Arguments
///
/// * `writer` - The output writer (can be stderr, file, or buffer for testing)
/// * `input` - The original input string that caused the error
/// * `error` - The error that occurred
/// * `format` - The output format
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::{output_error, OutputFormat};
/// use ferro_hgvs::FerroError;
/// use std::io::Cursor;
///
/// let mut buffer = Cursor::new(Vec::new());
/// let error = FerroError::Parse { msg: "test error".to_string(), pos: 0, diagnostic: None };
/// output_error(&mut buffer, "input", &error, OutputFormat::Text).unwrap();
/// let result = String::from_utf8(buffer.into_inner()).unwrap();
/// assert!(result.contains("ERROR: input"));
/// ```
pub fn output_error<W: Write>(
    writer: &mut W,
    input: &str,
    error: &FerroError,
    format: OutputFormat,
) -> io::Result<()> {
    output_error_with_context(writer, input, error, format, None)
}

/// Write an error to the output with optional line number context
///
/// # Arguments
///
/// * `writer` - The output writer (can be stderr, file, or buffer for testing)
/// * `input` - The original input string that caused the error
/// * `error` - The error that occurred
/// * `format` - The output format
/// * `line_number` - Optional line number in the input file
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::{output_error_with_context, OutputFormat};
/// use ferro_hgvs::FerroError;
/// use std::io::Cursor;
///
/// let mut buffer = Cursor::new(Vec::new());
/// let error = FerroError::Parse { msg: "test error".to_string(), pos: 0, diagnostic: None };
/// output_error_with_context(&mut buffer, "input", &error, OutputFormat::Text, Some(42)).unwrap();
/// let result = String::from_utf8(buffer.into_inner()).unwrap();
/// assert!(result.contains("line 42"));
/// ```
pub fn output_error_with_context<W: Write>(
    writer: &mut W,
    input: &str,
    error: &FerroError,
    format: OutputFormat,
    line_number: Option<usize>,
) -> io::Result<()> {
    match format {
        OutputFormat::Json => {
            if let Some(line) = line_number {
                writeln!(
                    writer,
                    r#"{{"input": "{}", "error": "{}", "line": {}, "status": "error"}}"#,
                    escape_json(input),
                    escape_json(&error.to_string()),
                    line
                )
            } else {
                writeln!(
                    writer,
                    r#"{{"input": "{}", "error": "{}", "status": "error"}}"#,
                    escape_json(input),
                    escape_json(&error.to_string())
                )
            }
        }
        OutputFormat::Text | OutputFormat::Vcf => {
            if let Some(line) = line_number {
                writeln!(writer, "ERROR (line {}): {} - {}", line, input, error)
            } else {
                writeln!(writer, "ERROR: {} - {}", input, error)
            }
        }
    }
}

/// Format a VCF-to-HGVS result
///
/// # Arguments
///
/// * `writer` - The output writer
/// * `chrom` - Chromosome name
/// * `pos` - Position
/// * `ref_allele` - Reference allele
/// * `alt_allele` - Alternate allele
/// * `hgvs` - The HGVS string
/// * `format` - The output format
pub fn output_vcf_to_hgvs<W: Write>(
    writer: &mut W,
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    hgvs: &str,
    format: OutputFormat,
) -> io::Result<()> {
    match format {
        OutputFormat::Json => {
            writeln!(
                writer,
                r#"{{"chrom": "{}", "pos": {}, "ref": "{}", "alt": "{}", "hgvs": "{}"}}"#,
                escape_json(chrom),
                pos,
                escape_json(ref_allele),
                escape_json(alt_allele),
                escape_json(hgvs)
            )
        }
        OutputFormat::Text | OutputFormat::Vcf => {
            writeln!(
                writer,
                "{}:{} {}/{} -> {}",
                chrom, pos, ref_allele, alt_allele, hgvs
            )
        }
    }
}

/// Format an HGVS-to-VCF result
///
/// # Arguments
///
/// * `writer` - The output writer
/// * `hgvs` - The original HGVS string
/// * `chrom` - Chromosome name
/// * `pos` - Position
/// * `ref_allele` - Reference allele
/// * `alt_allele` - Alternate allele
/// * `format` - The output format
pub fn output_hgvs_to_vcf<W: Write>(
    writer: &mut W,
    hgvs: &str,
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    format: OutputFormat,
) -> io::Result<()> {
    match format {
        OutputFormat::Json => {
            writeln!(
                writer,
                r#"{{"hgvs": "{}", "chrom": "{}", "pos": {}, "ref": "{}", "alt": "{}"}}"#,
                escape_json(hgvs),
                escape_json(chrom),
                pos,
                escape_json(ref_allele),
                escape_json(alt_allele)
            )
        }
        OutputFormat::Vcf => {
            writeln!(
                writer,
                "{}\t{}\t.\t{}\t{}",
                chrom, pos, ref_allele, alt_allele
            )
        }
        OutputFormat::Text => {
            writeln!(
                writer,
                "{} -> {}:{}:{}/{}",
                hgvs, chrom, pos, ref_allele, alt_allele
            )
        }
    }
}

/// Format a transcript annotation result
///
/// # Arguments
///
/// * `writer` - The output writer
/// * `transcript` - Transcript accession
/// * `gene` - Gene symbol
/// * `hgvs` - The HGVS string
/// * `format` - The output format
pub fn output_transcript_annotation<W: Write>(
    writer: &mut W,
    transcript: &str,
    gene: &str,
    hgvs: &str,
    format: OutputFormat,
) -> io::Result<()> {
    match format {
        OutputFormat::Json => {
            writeln!(
                writer,
                r#"{{"transcript": "{}", "gene": "{}", "hgvs": "{}"}}"#,
                escape_json(transcript),
                escape_json(gene),
                escape_json(hgvs)
            )
        }
        OutputFormat::Text | OutputFormat::Vcf => {
            writeln!(writer, "  {} ({}) -> {}", transcript, gene, hgvs)
        }
    }
}

/// Header row emitted once at the top of `ferro normalize --format tsv` output.
///
/// Kept next to [`normalize_tsv_row`] so the column order can only be changed in
/// one place; a consumer that pins column positions therefore cannot be broken
/// by an edit to only one of the two.
pub const NORMALIZE_TSV_HEADER: &str = "line\tinput\tnormalized\tchanged\tstatus\tdetail";

/// The phase in which a variant failed, as the short token written to the
/// `status` column of a `normalize --format tsv` row.
///
/// The phases are distinguished because they tell a caller different things: a
/// preprocessing rejection means the input was refused before parsing (usually
/// under `--error-mode strict`), a parse failure means the HGVS grammar rejected
/// it, and a normalization failure means it parsed but could not be normalized
/// (e.g. the reference lacks the accession).
///
/// # `status` is a property of `(input, error mode)`, not of the input alone
///
/// The token names **which stage rejected the input under the active error
/// mode**, and stages run in order: preprocess, then parse, then normalize. The
/// preprocessor is error-mode-driven, so widening or narrowing `--error-mode`
/// moves the boundary between `preprocess_error` and the later tokens *for the
/// same input*. A grammar violation the preprocessor happens to inspect first is
/// therefore reported as `preprocess_error` when that mode rejects it — for
/// example `NM_000088.3:c.(1_2)ins` under `--error-mode strict` yields
/// `preprocess_error` with a parse-shaped `detail`, while a mode that lets it
/// through reaches the parser and yields `parse_error`.
///
/// A consumer that wants to bucket "malformed input" against "policy rejection"
/// must therefore key off the error mode it passed, not off the token alone;
/// treating the token as an intrinsic classification of the variant string will
/// mis-bucket exactly these cases. What the token *is* good for is telling a
/// caller how far a given run got, and — held at a fixed error mode — it is
/// stable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizeTsvFailure {
    /// Input rejected by the error-mode preprocessor, before parsing.
    Preprocess,
    /// Input could not be parsed as HGVS.
    Parse,
    /// Input parsed but could not be normalized.
    Normalize,
}

impl NormalizeTsvFailure {
    /// The stable, machine-readable token for this phase.
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::cli::NormalizeTsvFailure;
    ///
    /// assert_eq!(NormalizeTsvFailure::Parse.token(), "parse_error");
    /// ```
    pub fn token(self) -> &'static str {
        match self {
            NormalizeTsvFailure::Preprocess => "preprocess_error",
            NormalizeTsvFailure::Parse => "parse_error",
            NormalizeTsvFailure::Normalize => "normalize_error",
        }
    }
}

/// What became of one input variant, as rendered by [`normalize_tsv_row`].
///
/// Modelled as a sum type so a row can never claim both a normalized string and
/// a failure status: the `changed` and `status` columns are derived from this,
/// never passed in alongside it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalizeTsvOutcome<'a> {
    /// Normalization succeeded.
    Normalized {
        /// The normalized HGVS string.
        normalized: &'a str,
        /// Preprocessing corrections applied on the way, as `(code, message)`
        /// pairs. Rendered into the `detail` column so a row that is `changed`
        /// only because the input was corrected says so *in the table*, rather
        /// than only in an un-correlatable stderr warning.
        corrections: &'a [(&'a str, &'a str)],
    },
    /// Normalization did not happen, carrying the failing phase and its message.
    Failed(NormalizeTsvFailure, &'a str),
}

/// Whether `c` would break the "one row, one physical line, six fields" shape.
///
/// [`char::is_control`] covers the obvious tab/LF/CR plus the rest of C0/C1 —
/// notably VT (`U+000B`), FF (`U+000C`), the file/record/group/unit separators
/// (`U+001C`–`U+001F`) and ESC (`U+001B`), which would inject ANSI escapes into
/// a terminal that `cat`s the file. The three explicit cases are the non-control
/// separators: NEL (`U+0085`) is C1 (and so already a control char, listed here
/// for symmetry), and LS/PS (`U+2028`/`U+2029`) are not control characters at
/// all. All of them are line breaks to Python's `str.splitlines()` — the most
/// likely way a pipeline reads this table — so leaving any of them through turns
/// one row into two.
fn breaks_tsv_row(c: char) -> bool {
    c.is_control() || matches!(c, '\u{0085}' | '\u{2028}' | '\u{2029}')
}

/// Replace every run of row-breaking characters with a single space.
///
/// An error message is arbitrary text and may legitimately contain a tab or a
/// line break; emitting it raw would split one logical row into several physical
/// ones (or shift the columns), silently corrupting the table for every
/// downstream consumer. Substituting a space keeps the row readable while making
/// the "one line, six fields" shape unconditional. Runs collapse to one space so
/// a `\r\n` (or an indented multi-line diagnostic) does not leave a trail of
/// whitespace behind.
fn sanitize_tsv_field(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut in_run = false;
    for c in s.chars() {
        if breaks_tsv_row(c) {
            if !in_run {
                out.push(' ');
                in_run = true;
            }
        } else {
            out.push(c);
            in_run = false;
        }
    }
    out
}

/// One rendered `normalize --format tsv` row plus the `changed` verdict it was
/// built from.
///
/// The flag is handed back rather than recomputed by the caller so the run
/// summary's `changed` count and the row's `changed` column cannot drift apart:
/// there is exactly one comparison, here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizeTsvRow {
    /// The rendered row, without a trailing newline.
    pub row: String,
    /// Whether the normalized string differed from the input. Always `false` for
    /// a failure, which has nothing to compare — the row's `changed` column is
    /// left *empty* in that case (see [`normalize_tsv_row`]), so a caller
    /// counting changes and a caller filtering the column agree.
    pub changed: bool,
}

/// Render one `ferro normalize --format tsv` row (no trailing newline).
///
/// Columns are `line`, `input`, `normalized`, `changed`, `status`, `detail`,
/// matching [`NORMALIZE_TSV_HEADER`] and documented in the [module
/// docs](self#ferro-normalize---format-tsv). `line`, `normalized`, `changed` and
/// `detail` are empty for the cases where they do not apply, and every text field
/// is sanitized so the row cannot contain a character that would split or shift
/// it.
///
/// # Arguments
///
/// * `line` - 1-based input line number, or `None` for a single positional variant
/// * `input` - The variant exactly as read from the input
/// * `outcome` - What normalization produced for it
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::{normalize_tsv_row, NormalizeTsvFailure, NormalizeTsvOutcome};
///
/// let changed = normalize_tsv_row(
///     Some(4),
///     "NM_000088.3:c.459delA",
///     NormalizeTsvOutcome::Normalized {
///         normalized: "NM_000088.3:c.459del",
///         corrections: &[],
///     },
/// );
/// assert!(changed.changed);
/// assert_eq!(changed.row, "4\tNM_000088.3:c.459delA\tNM_000088.3:c.459del\ttrue\tok\t");
///
/// // A failure leaves both `normalized` and `changed` empty.
/// let failed = normalize_tsv_row(
///     None,
///     "bogus",
///     NormalizeTsvOutcome::Failed(NormalizeTsvFailure::Parse, "not an accession"),
/// );
/// assert!(!failed.changed);
/// assert_eq!(failed.row, "\tbogus\t\t\tparse_error\tnot an accession");
/// ```
pub fn normalize_tsv_row(
    line: Option<usize>,
    input: &str,
    outcome: NormalizeTsvOutcome<'_>,
) -> NormalizeTsvRow {
    // One match yields every derived piece, so no two columns can disagree about
    // which arm they came from. The comparison is on the *unsanitized* strings:
    // sanitization is a rendering concern, and an HGVS string never contains a
    // character that sanitization would touch anyway.
    //
    // `changed` renders as the empty string on failure rather than `false`: a
    // failed row has no normalized string, so claiming it is "unchanged" would be
    // a false statement about it — and would contradict `normalize_tsv_summary`,
    // whose `unchanged` count excludes failures.
    let (normalized, changed, changed_col, status, detail) = match outcome {
        NormalizeTsvOutcome::Normalized {
            normalized,
            corrections,
        } => {
            let changed = normalized != input;
            (
                normalized,
                changed,
                if changed { "true" } else { "false" },
                "ok",
                join_tsv_corrections(corrections),
            )
        }
        NormalizeTsvOutcome::Failed(phase, detail) => {
            ("", false, "", phase.token(), detail.to_string())
        }
    };
    let row = format!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        line.map(|n| n.to_string()).unwrap_or_default(),
        sanitize_tsv_field(input),
        sanitize_tsv_field(normalized),
        changed_col,
        status,
        sanitize_tsv_field(&detail)
    );
    NormalizeTsvRow { row, changed }
}

/// Render preprocessing corrections as one `detail` value: `CODE: message`,
/// `; `-joined.
///
/// The separator is `; ` because a correction message may itself contain a comma
/// or a colon, so a more common separator would be ambiguous to split on.
fn join_tsv_corrections(corrections: &[(&str, &str)]) -> String {
    corrections
        .iter()
        .map(|(code, message)| format!("{code}: {message}"))
        .collect::<Vec<_>>()
        .join("; ")
}

/// Render the one-line `normalize --format tsv` run summary (no trailing newline).
///
/// Written to stderr, never to the output stream, so stdout stays a table a
/// consumer can parse without having to skip a trailing non-row line.
///
/// # Arguments
///
/// * `total` - Rows emitted (every non-skipped input variant)
/// * `changed` - Rows whose normalized string differed from the input
/// * `failed` - Rows that carry a failure status
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::normalize_tsv_summary;
///
/// assert_eq!(
///     normalize_tsv_summary(10, 3, 1),
///     "summary: total=10 changed=3 unchanged=6 failed=1"
/// );
/// ```
pub fn normalize_tsv_summary(total: usize, changed: usize, failed: usize) -> String {
    // Unchanged counts only *successful* rows: a failed row has no normalized
    // string, so calling it "unchanged" would be a false claim about it.
    let unchanged = total.saturating_sub(changed).saturating_sub(failed);
    format!("summary: total={total} changed={changed} unchanged={unchanged} failed={failed}")
}

/// Render projection warnings as a JSON array of `{code, message}` objects.
///
/// Always emits an array — `[]` when there are none — so consumers can index
/// `.warnings` unconditionally rather than branching on the key's presence.
/// Hand-built for consistency with the surrounding hand-built JSON in this
/// module (which is why `escape_json` exists at all). #1182.
fn warnings_json(warnings: &[crate::cli::project::ProjectionWarning]) -> String {
    let mut out = String::from("[");
    for (i, w) in warnings.iter().enumerate() {
        if i > 0 {
            out.push_str(", ");
        }
        out.push_str(&format!(
            r#"{{"code": "{}", "message": "{}"}}"#,
            escape_json(&w.code),
            escape_json(&w.message)
        ));
    }
    out.push(']');
    out
}

/// Escape special characters for JSON strings
///
/// Escapes backslashes, quotes, and control characters.
fn escape_json(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => result.push_str("\\\""),
            '\\' => result.push_str("\\\\"),
            '\n' => result.push_str("\\n"),
            '\r' => result.push_str("\\r"),
            '\t' => result.push_str("\\t"),
            c if c.is_control() => {
                result.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => result.push(c),
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    // ===== OutputFormat Tests =====

    #[test]
    fn test_output_format_from_str() {
        use std::str::FromStr;
        assert_eq!(OutputFormat::from_str("json").unwrap(), OutputFormat::Json);
        assert_eq!(OutputFormat::from_str("JSON").unwrap(), OutputFormat::Json);
        assert_eq!(OutputFormat::from_str("Json").unwrap(), OutputFormat::Json);
        assert_eq!(OutputFormat::from_str("text").unwrap(), OutputFormat::Text);
        assert_eq!(OutputFormat::from_str("TEXT").unwrap(), OutputFormat::Text);
        assert_eq!(OutputFormat::from_str("vcf").unwrap(), OutputFormat::Vcf);
        assert_eq!(OutputFormat::from_str("VCF").unwrap(), OutputFormat::Vcf);
    }

    #[test]
    fn test_output_format_default() {
        use std::str::FromStr;
        assert_eq!(
            OutputFormat::from_str("unknown").unwrap(),
            OutputFormat::Text
        );
        assert_eq!(OutputFormat::from_str("").unwrap(), OutputFormat::Text);
        assert_eq!(OutputFormat::default(), OutputFormat::Text);
    }

    // ===== output_result Tests =====

    #[test]
    fn test_output_result_text_changed() {
        let mut buffer = Cursor::new(Vec::new());
        output_result(&mut buffer, "input", "output", OutputFormat::Text).unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert_eq!(result, "input -> output\n");
    }

    #[test]
    fn test_output_result_text_unchanged() {
        let mut buffer = Cursor::new(Vec::new());
        output_result(&mut buffer, "same", "same", OutputFormat::Text).unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert_eq!(result, "same\n");
    }

    #[test]
    fn test_output_result_json() {
        let mut buffer = Cursor::new(Vec::new());
        output_result(&mut buffer, "input", "output", OutputFormat::Json).unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.contains(r#""input": "input""#));
        assert!(result.contains(r#""output": "output""#));
        assert!(result.contains(r#""status": "ok""#));
    }

    #[test]
    fn test_output_result_json_escaping() {
        let mut buffer = Cursor::new(Vec::new());
        output_result(
            &mut buffer,
            "with\"quote",
            "with\\slash",
            OutputFormat::Json,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.contains(r#"with\"quote"#));
        assert!(result.contains(r#"with\\slash"#));
    }

    // ===== output_error Tests =====

    #[test]
    fn test_output_error_text() {
        let mut buffer = Cursor::new(Vec::new());
        let error = FerroError::Parse {
            msg: "test error".to_string(),
            pos: 0,
            diagnostic: None,
        };
        output_error(&mut buffer, "input", &error, OutputFormat::Text).unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.starts_with("ERROR: input"));
        assert!(result.contains("test error"));
    }

    #[test]
    fn test_output_error_json() {
        let mut buffer = Cursor::new(Vec::new());
        let error = FerroError::Parse {
            msg: "test error".to_string(),
            pos: 0,
            diagnostic: None,
        };
        output_error(&mut buffer, "input", &error, OutputFormat::Json).unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.contains(r#""input": "input""#));
        assert!(result.contains(r#""status": "error""#));
    }

    // ===== output_vcf_to_hgvs Tests =====

    #[test]
    fn test_output_vcf_to_hgvs_text() {
        let mut buffer = Cursor::new(Vec::new());
        output_vcf_to_hgvs(
            &mut buffer,
            "chr1",
            12345,
            "A",
            "G",
            "g.12345A>G",
            OutputFormat::Text,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert_eq!(result, "chr1:12345 A/G -> g.12345A>G\n");
    }

    #[test]
    fn test_output_vcf_to_hgvs_json() {
        let mut buffer = Cursor::new(Vec::new());
        output_vcf_to_hgvs(
            &mut buffer,
            "chr1",
            12345,
            "A",
            "G",
            "g.12345A>G",
            OutputFormat::Json,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.contains(r#""chrom": "chr1""#));
        assert!(result.contains(r#""pos": 12345"#));
        assert!(result.contains(r#""ref": "A""#));
        assert!(result.contains(r#""alt": "G""#));
        assert!(result.contains(r#""hgvs": "g.12345A>G""#));
    }

    // ===== output_hgvs_to_vcf Tests =====

    #[test]
    fn test_output_hgvs_to_vcf_text() {
        let mut buffer = Cursor::new(Vec::new());
        output_hgvs_to_vcf(
            &mut buffer,
            "NC_000001.11:g.12345A>G",
            "chr1",
            12345,
            "A",
            "G",
            OutputFormat::Text,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert_eq!(result, "NC_000001.11:g.12345A>G -> chr1:12345:A/G\n");
    }

    #[test]
    fn test_output_hgvs_to_vcf_vcf_format() {
        let mut buffer = Cursor::new(Vec::new());
        output_hgvs_to_vcf(
            &mut buffer,
            "NC_000001.11:g.12345A>G",
            "chr1",
            12345,
            "A",
            "G",
            OutputFormat::Vcf,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert_eq!(result, "chr1\t12345\t.\tA\tG\n");
    }

    #[test]
    fn test_output_hgvs_to_vcf_json() {
        let mut buffer = Cursor::new(Vec::new());
        output_hgvs_to_vcf(
            &mut buffer,
            "NC_000001.11:g.12345A>G",
            "chr1",
            12345,
            "A",
            "G",
            OutputFormat::Json,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.contains(r#""hgvs": "NC_000001.11:g.12345A>G""#));
        assert!(result.contains(r#""chrom": "chr1""#));
    }

    // ===== output_transcript_annotation Tests =====

    #[test]
    fn test_output_transcript_annotation_text() {
        let mut buffer = Cursor::new(Vec::new());
        output_transcript_annotation(
            &mut buffer,
            "NM_000088.3",
            "BRCA1",
            "c.100A>G",
            OutputFormat::Text,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert_eq!(result, "  NM_000088.3 (BRCA1) -> c.100A>G\n");
    }

    #[test]
    fn test_output_transcript_annotation_json() {
        let mut buffer = Cursor::new(Vec::new());
        output_transcript_annotation(
            &mut buffer,
            "NM_000088.3",
            "BRCA1",
            "c.100A>G",
            OutputFormat::Json,
        )
        .unwrap();
        let result = String::from_utf8(buffer.into_inner()).unwrap();
        assert!(result.contains(r#""transcript": "NM_000088.3""#));
        assert!(result.contains(r#""gene": "BRCA1""#));
        assert!(result.contains(r#""hgvs": "c.100A>G""#));
    }

    // ===== output_projection Tests =====

    #[test]
    fn output_projection_text_rendered_is_bare_hgvs() {
        use crate::cli::project::{Axis, AxisOutcome};
        let mut buf = Cursor::new(Vec::new());
        let outcome = AxisOutcome::Rendered {
            transcript_id: "NM_000088.3".to_string(),
            output: "NP_000079.2:p.(Gly197Cys)".to_string(),
            warnings: Vec::new(),
        };
        output_projection(
            &mut buf,
            "NM_000088.3:c.589G>T",
            Axis::Protein,
            &outcome,
            OutputFormat::Text,
            None,
        )
        .unwrap();
        let out = String::from_utf8(buf.into_inner()).unwrap();
        assert_eq!(out, "NP_000079.2:p.(Gly197Cys)\n");
    }

    #[test]
    fn output_projection_json_rendered_has_fields() {
        use crate::cli::project::{Axis, AxisOutcome, ProjectionWarning};
        let mut buf = Cursor::new(Vec::new());
        // Non-empty on purpose: surfacing warnings is the point of #1182, and an
        // empty vector cannot distinguish "serialised correctly" from "silently
        // dropped" — `warnings_json` renders `[]` either way.
        let outcome = AxisOutcome::Rendered {
            transcript_id: "NM_000088.3".to_string(),
            output: "NP_000079.2:p.(Gly197Cys)".to_string(),
            warnings: vec![ProjectionWarning {
                code: "POSITION_PAST_END".to_string(),
                message: "position 999 is past the end".to_string(),
            }],
        };
        output_projection(
            &mut buf,
            "NM_000088.3:c.589G>T",
            Axis::Protein,
            &outcome,
            OutputFormat::Json,
            None,
        )
        .unwrap();
        let out = String::from_utf8(buf.into_inner()).unwrap();
        assert!(out.contains(r#""axis": "p""#));
        assert!(out.contains(r#""transcript": "NM_000088.3""#));
        assert!(out.contains(r#""output": "NP_000079.2:p.(Gly197Cys)""#));
        assert!(out.contains(r#""status": "ok""#));
        assert!(
            out.contains(r#""warnings": [{"code": "POSITION_PAST_END", "message": "position 999 is past the end"}]"#),
            "a rendered projection must serialise its warnings; got: {out}"
        );
    }

    #[test]
    fn output_projection_json_unavailable_status() {
        use crate::cli::project::{Axis, AxisOutcome, ProjectionWarning};
        let mut buf = Cursor::new(Vec::new());
        let outcome = AxisOutcome::Unavailable {
            transcript_id: Some("NM_000088.3".to_string()),
            reason: "no p. representation for this variant".to_string(),
            // See the note in the `rendered` test: the unavailable branch builds
            // its JSON separately, so it needs its own non-empty case.
            warnings: vec![ProjectionWarning {
                code: "INTRONIC".to_string(),
                message: "intronic offset dropped".to_string(),
            }],
        };
        output_projection(
            &mut buf,
            "NM_000088.3:c.100-5A>G",
            Axis::Protein,
            &outcome,
            OutputFormat::Json,
            None,
        )
        .unwrap();
        let out = String::from_utf8(buf.into_inner()).unwrap();
        assert!(out.contains(r#""status": "unavailable""#));
        assert!(out.contains(r#""reason""#));
        assert!(
            out.contains(
                r#""warnings": [{"code": "INTRONIC", "message": "intronic offset dropped"}]"#
            ),
            "an unavailable projection must serialise its warnings too; got: {out}"
        );
    }

    #[test]
    fn output_projection_text_unavailable_is_comment() {
        use crate::cli::project::{Axis, AxisOutcome};
        let mut buf = Cursor::new(Vec::new());
        let outcome = AxisOutcome::Unavailable {
            transcript_id: Some("NM_000088.3".to_string()),
            reason: "no p. representation for this variant".to_string(),
            warnings: Vec::new(),
        };
        output_projection(
            &mut buf,
            "NM_000088.3:c.100-5A>G",
            Axis::Protein,
            &outcome,
            OutputFormat::Text,
            None,
        )
        .unwrap();
        let out = String::from_utf8(buf.into_inner()).unwrap();
        assert!(
            out.starts_with("# NM_000088.3:c.100-5A>G [p]: unavailable:"),
            "got {out}"
        );
    }

    #[test]
    fn output_project_error_json_is_parseable_object() {
        let mut buf = Cursor::new(Vec::new());
        output_project_error(
            &mut buf,
            "NC_000001.11:g.1003C>A",
            "a genomic (g.) input requires --transcript",
            OutputFormat::Json,
            Some(3),
        )
        .unwrap();
        let out = String::from_utf8(buf.into_inner()).unwrap();
        assert!(out.contains(r#""status": "error""#));
        assert!(out.contains(r#""line": 3"#));
        assert!(out.contains(r#""error": "a genomic (g.) input requires --transcript""#));
    }

    #[test]
    fn output_project_error_text_is_plain() {
        let mut buf = Cursor::new(Vec::new());
        output_project_error(&mut buf, "x", "boom", OutputFormat::Text, None).unwrap();
        let out = String::from_utf8(buf.into_inner()).unwrap();
        assert_eq!(out, "ERROR: x - boom\n");
    }

    // ===== escape_json Tests =====

    #[test]
    fn test_escape_json_quotes() {
        assert_eq!(escape_json(r#"hello"world"#), r#"hello\"world"#);
    }

    #[test]
    fn test_escape_json_backslash() {
        assert_eq!(escape_json(r"hello\world"), r"hello\\world");
    }

    #[test]
    fn test_escape_json_newline() {
        assert_eq!(escape_json("hello\nworld"), r"hello\nworld");
    }

    #[test]
    fn test_escape_json_tab() {
        assert_eq!(escape_json("hello\tworld"), r"hello\tworld");
    }

    #[test]
    fn test_escape_json_carriage_return() {
        assert_eq!(escape_json("hello\rworld"), r"hello\rworld");
    }

    #[test]
    fn test_escape_json_no_escaping() {
        assert_eq!(escape_json("hello world"), "hello world");
    }

    #[test]
    fn test_escape_json_empty() {
        assert_eq!(escape_json(""), "");
    }

    #[test]
    fn test_escape_json_combined() {
        assert_eq!(
            escape_json("line1\nline2\t\"quoted\"\r\\"),
            r#"line1\nline2\t\"quoted\"\r\\"#
        );
    }

    // ===== normalize --format tsv Tests =====

    /// Build a `Normalized` outcome with no preprocessing corrections.
    fn normalized(s: &str) -> NormalizeTsvOutcome<'_> {
        NormalizeTsvOutcome::Normalized {
            normalized: s,
            corrections: &[],
        }
    }

    /// Split a rendered row into its fields, asserting the row is one physical
    /// line with exactly as many fields as the header declares, and that no field
    /// carries a character that any reasonable line-splitter would break on.
    ///
    /// The "one physical line" check deliberately does *not* go through
    /// `str::lines()`: that splits on exactly `\n`/`\r\n`, i.e. it shares the
    /// blind spot an incomplete sanitizer would have. Inspecting for
    /// [`char::is_control`] plus the three non-control separators instead means
    /// the assertion is independent of the implementation's notion of a break.
    fn tsv_fields(row: &str) -> Vec<&str> {
        let expected = NORMALIZE_TSV_HEADER.split('\t').count();
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(
            fields.len(),
            expected,
            "row does not have {expected} fields: {row:?}"
        );
        for c in row.chars() {
            assert!(
                !(c.is_control() || matches!(c, '\u{0085}' | '\u{2028}' | '\u{2029}')) || c == '\t',
                "row carries {c:?} (U+{:04X}), which splits or shifts it: {row:?}",
                c as u32
            );
        }
        fields
    }

    #[test]
    fn tsv_header_names_the_six_documented_columns() {
        assert_eq!(
            tsv_fields(NORMALIZE_TSV_HEADER),
            vec!["line", "input", "normalized", "changed", "status", "detail"]
        );
    }

    #[test]
    fn tsv_row_marks_an_unchanged_variant() {
        let row = normalize_tsv_row(
            Some(7),
            "NM_000088.3:c.459A>G",
            normalized("NM_000088.3:c.459A>G"),
        );
        assert!(!row.changed);
        assert_eq!(
            tsv_fields(&row.row),
            vec![
                "7",
                "NM_000088.3:c.459A>G",
                "NM_000088.3:c.459A>G",
                "false",
                "ok",
                ""
            ]
        );
    }

    #[test]
    fn tsv_row_marks_a_changed_variant() {
        let row = normalize_tsv_row(
            Some(1),
            "NM_000088.3:c.459delA",
            normalized("NM_000088.3:c.459del"),
        );
        assert!(row.changed);
        assert_eq!(
            tsv_fields(&row.row),
            vec![
                "1",
                "NM_000088.3:c.459delA",
                "NM_000088.3:c.459del",
                "true",
                "ok",
                ""
            ]
        );
    }

    #[test]
    fn tsv_row_omits_the_line_column_without_line_context() {
        // The single-positional path has no line number to report; the column is
        // empty rather than a fabricated `1`.
        let row = normalize_tsv_row(None, "x", normalized("x"));
        assert_eq!(tsv_fields(&row.row)[0], "");
    }

    #[test]
    fn tsv_row_reports_preprocessing_corrections_as_the_detail() {
        // The most common reason a row is `changed` is a preprocessing
        // correction; the table has to say which one, or the change is
        // unexplained.
        let row = normalize_tsv_row(
            Some(2),
            "NM_000088.3:c.459delA",
            NormalizeTsvOutcome::Normalized {
                normalized: "NM_000088.3:c.459del",
                corrections: &[
                    ("W3025", "deletion with explicit deleted sequence"),
                    ("W2003", "extra whitespace removed"),
                ],
            },
        );
        assert_eq!(
            tsv_fields(&row.row)[5],
            "W3025: deletion with explicit deleted sequence; W2003: extra whitespace removed"
        );
    }

    #[test]
    fn tsv_row_leaves_normalized_and_changed_empty_on_failure() {
        // `changed` must be empty, not `false`: `normalize_tsv_summary` excludes
        // failures from its `unchanged` count, so a `false` here would make
        // `awk -F'\t' '$4=="false"'` and the summary disagree.
        for (phase, token) in [
            (NormalizeTsvFailure::Preprocess, "preprocess_error"),
            (NormalizeTsvFailure::Parse, "parse_error"),
            (NormalizeTsvFailure::Normalize, "normalize_error"),
        ] {
            let row =
                normalize_tsv_row(Some(3), "bogus", NormalizeTsvOutcome::Failed(phase, "why"));
            assert!(!row.changed);
            assert_eq!(
                tsv_fields(&row.row),
                vec!["3", "bogus", "", "", token, "why"]
            );
        }
    }

    #[test]
    fn tsv_row_sanitizes_every_row_breaking_character() {
        // The motivating case: an error message containing a tab and a newline
        // would otherwise split one row into three and shift the columns of the
        // first, corrupting the whole table from that point on. The exotic
        // separators matter just as much — Python's `str.splitlines()`, the
        // likeliest reader, breaks on all of them.
        let row = normalize_tsv_row(
            Some(1),
            "in\tput\nwith\u{000B}breaks\u{2028}here",
            NormalizeTsvOutcome::Failed(
                NormalizeTsvFailure::Parse,
                "one\u{000C}two\u{001C}three\u{0085}four\u{2029}five\u{001B}[31m\r\n",
            ),
        );
        // `tsv_fields` itself asserts the row carries no breaking character.
        assert_eq!(
            tsv_fields(&row.row),
            vec![
                "1",
                "in put with breaks here",
                "",
                "",
                "parse_error",
                // Runs collapse, so the trailing `\r\n` leaves one space.
                "one two three four five [31m ",
            ]
        );
    }

    #[test]
    fn tsv_row_sanitizes_a_normalized_string_too() {
        // `normalized` comes from `Display` on a parsed variant and so should
        // never contain a tab — sanitizing it anyway means no future `Display`
        // change can corrupt the table.
        let row = normalize_tsv_row(Some(1), "in", normalized("a\tb"));
        assert_eq!(
            tsv_fields(&row.row),
            vec!["1", "in", "a b", "true", "ok", ""]
        );
    }

    #[test]
    fn tsv_summary_derives_unchanged_from_total() {
        assert_eq!(
            normalize_tsv_summary(10, 3, 1),
            "summary: total=10 changed=3 unchanged=6 failed=1"
        );
        assert_eq!(
            normalize_tsv_summary(0, 0, 0),
            "summary: total=0 changed=0 unchanged=0 failed=0"
        );
    }

    #[test]
    fn tsv_summary_never_underflows() {
        // Defensive: the counters are maintained by the CLI, and a saturating
        // subtraction is preferable to a panic in a summary line.
        assert_eq!(
            normalize_tsv_summary(1, 5, 5),
            "summary: total=1 changed=5 unchanged=0 failed=5"
        );
    }
}
