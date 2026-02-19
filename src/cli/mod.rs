//! CLI utilities for ferro-hgvs
//!
//! This module provides testable functions used by the CLI binary.
//! By extracting pure functions and I/O-abstracted functions to the library,
//! we enable comprehensive unit testing without requiring end-to-end CLI tests.

pub mod format;
pub mod parse;
pub mod sequence;

// Re-export commonly used items
pub use format::{output_error, output_error_with_context, output_result, OutputFormat};
pub use parse::{
    parse_genome_build, parse_shuffle_direction, parse_vcf_line, parse_vcf_line_with_build,
};
pub use sequence::reverse_complement;

/// UTF-8 BOM (Byte Order Mark) constant
const UTF8_BOM: &str = "\u{feff}";

/// Strip UTF-8 BOM from the beginning of a string if present.
///
/// This is common when files are exported from Windows applications or Excel.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::strip_bom;
///
/// assert_eq!(strip_bom("\u{feff}NM_000088.3:c.459del"), "NM_000088.3:c.459del");
/// assert_eq!(strip_bom("NM_000088.3:c.459del"), "NM_000088.3:c.459del");
/// ```
pub fn strip_bom(s: &str) -> &str {
    s.strip_prefix(UTF8_BOM).unwrap_or(s)
}

/// Strip inline comments from a variant line.
///
/// Comments start with `#` and extend to the end of the line.
/// Leading/trailing whitespace is also trimmed.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::strip_inline_comment;
///
/// assert_eq!(strip_inline_comment("NM_000088.3:c.459del  # my note"), "NM_000088.3:c.459del");
/// assert_eq!(strip_inline_comment("NM_000088.3:c.459del"), "NM_000088.3:c.459del");
/// assert_eq!(strip_inline_comment("# full line comment"), "");
/// ```
pub fn strip_inline_comment(s: &str) -> &str {
    match s.find('#') {
        Some(pos) => s[..pos].trim(),
        None => s.trim(),
    }
}

/// Process an input line: trim whitespace, strip BOM (for first line), and strip inline comments.
///
/// Returns None if the line is empty or a comment-only line.
///
/// # Arguments
///
/// * `line` - The input line to process
/// * `is_first_line` - Whether this is the first line of input (for BOM handling).
///   Note: UTF-8 BOM only appears at the beginning of a file, never on subsequent lines,
///   so we intentionally only check for it on the first line.
///
/// # Examples
///
/// ```
/// use ferro_hgvs::cli::process_input_line;
///
/// // Normal line
/// assert_eq!(process_input_line("NM_000088.3:c.459del", false), Some("NM_000088.3:c.459del"));
///
/// // Line with inline comment
/// assert_eq!(process_input_line("NM_000088.3:c.459del  # note", false), Some("NM_000088.3:c.459del"));
///
/// // First line with BOM
/// assert_eq!(process_input_line("\u{feff}NM_000088.3:c.459del", true), Some("NM_000088.3:c.459del"));
///
/// // Empty line
/// assert_eq!(process_input_line("", false), None);
///
/// // Comment-only line
/// assert_eq!(process_input_line("# comment", false), None);
/// ```
pub fn process_input_line(line: &str, is_first_line: bool) -> Option<&str> {
    let line = line.trim();
    let line = if is_first_line { strip_bom(line) } else { line };
    let line = strip_inline_comment(line);

    if line.is_empty() {
        None
    } else {
        Some(line)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strip_bom() {
        assert_eq!(strip_bom("\u{feff}test"), "test");
        assert_eq!(strip_bom("test"), "test");
        assert_eq!(strip_bom("\u{feff}"), "");
        assert_eq!(strip_bom(""), "");
    }

    #[test]
    fn test_strip_inline_comment() {
        assert_eq!(strip_inline_comment("variant  # comment"), "variant");
        assert_eq!(strip_inline_comment("variant#comment"), "variant");
        assert_eq!(strip_inline_comment("# full comment"), "");
        assert_eq!(strip_inline_comment("variant"), "variant");
        assert_eq!(strip_inline_comment("  variant  "), "variant");
    }

    #[test]
    fn test_process_input_line() {
        // Normal cases
        assert_eq!(process_input_line("variant", false), Some("variant"));
        assert_eq!(
            process_input_line("variant  # comment", false),
            Some("variant")
        );

        // BOM handling on first line
        assert_eq!(process_input_line("\u{feff}variant", true), Some("variant"));
        // BOM on non-first line is preserved (edge case - user has malformed input)
        // The BOM character is not whitespace so trim() won't remove it
        assert_eq!(
            process_input_line("\u{feff}variant", false),
            Some("\u{feff}variant")
        );

        // Empty/comment lines
        assert_eq!(process_input_line("", false), None);
        assert_eq!(process_input_line("   ", false), None);
        assert_eq!(process_input_line("# comment", false), None);
    }
}
