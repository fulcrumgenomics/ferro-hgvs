//! Output formatting utilities for CLI operations

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
}
