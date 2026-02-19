//! Correction logic for each error type.
//!
//! This module provides functions to detect and correct specific error types
//! in HGVS input strings.

use super::types::{ErrorType, ResolvedAction};
use std::cmp::min;

/// A detected correction with its details.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DetectedCorrection {
    /// The type of error detected.
    pub error_type: ErrorType,
    /// The original text that was corrected.
    pub original: String,
    /// The corrected text.
    pub corrected: String,
    /// Start position (byte offset) in the original input.
    pub start: usize,
    /// End position (byte offset) in the original input.
    pub end: usize,
}

impl DetectedCorrection {
    /// Create a new detected correction.
    pub fn new(
        error_type: ErrorType,
        original: impl Into<String>,
        corrected: impl Into<String>,
        start: usize,
        end: usize,
    ) -> Self {
        Self {
            error_type,
            original: original.into(),
            corrected: corrected.into(),
            start,
            end,
        }
    }

    /// Format a warning message for this correction.
    pub fn warning_message(&self) -> String {
        format!(
            "{} at position {}: '{}' → '{}'",
            self.error_type, self.start, self.original, self.corrected
        )
    }
}

/// Normalize wrong dash characters (en-dash, em-dash) to hyphen.
///
/// Returns the corrected string and a list of corrections made.
pub fn correct_dash_characters(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut result = String::with_capacity(input.len());
    let mut corrections = Vec::new();

    for (i, c) in input.char_indices() {
        match c {
            // En-dash (U+2013)
            '\u{2013}' => {
                corrections.push(DetectedCorrection::new(
                    ErrorType::WrongDashCharacter,
                    c.to_string(),
                    "-",
                    i,
                    i + c.len_utf8(),
                ));
                result.push('-');
            }
            // Em-dash (U+2014)
            '\u{2014}' => {
                corrections.push(DetectedCorrection::new(
                    ErrorType::WrongDashCharacter,
                    c.to_string(),
                    "-",
                    i,
                    i + c.len_utf8(),
                ));
                result.push('-');
            }
            // Minus sign (U+2212)
            '\u{2212}' => {
                corrections.push(DetectedCorrection::new(
                    ErrorType::WrongDashCharacter,
                    c.to_string(),
                    "-",
                    i,
                    i + c.len_utf8(),
                ));
                result.push('-');
            }
            _ => result.push(c),
        }
    }

    (result, corrections)
}

/// Normalize wrong quote characters (smart quotes) to regular quotes.
///
/// Returns the corrected string and a list of corrections made.
pub fn correct_quote_characters(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut result = String::with_capacity(input.len());
    let mut corrections = Vec::new();

    for (i, c) in input.char_indices() {
        match c {
            // Left single quote (U+2018), right single quote (U+2019)
            '\u{2018}' | '\u{2019}' => {
                corrections.push(DetectedCorrection::new(
                    ErrorType::WrongQuoteCharacter,
                    c.to_string(),
                    "'",
                    i,
                    i + c.len_utf8(),
                ));
                result.push('\'');
            }
            // Left double quote (U+201C), right double quote (U+201D)
            '\u{201C}' | '\u{201D}' => {
                corrections.push(DetectedCorrection::new(
                    ErrorType::WrongQuoteCharacter,
                    c.to_string(),
                    "\"",
                    i,
                    i + c.len_utf8(),
                ));
                result.push('"');
            }
            _ => result.push(c),
        }
    }

    (result, corrections)
}

/// Normalize extra whitespace in the input.
///
/// This removes leading/trailing whitespace and collapses internal
/// whitespace around significant characters.
///
/// Returns the corrected string and a list of corrections made.
pub fn correct_whitespace(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();
    let trimmed = input.trim();

    // Track if we trimmed leading/trailing whitespace
    if input != trimmed {
        let leading = input.len() - input.trim_start().len();
        let trailing = input.trim_end().len();
        if leading > 0 || trailing < input.len() {
            corrections.push(DetectedCorrection::new(
                ErrorType::ExtraWhitespace,
                input.to_string(),
                trimmed.to_string(),
                0,
                input.len(),
            ));
        }
    }

    // Remove spaces around significant HGVS characters
    let mut result = String::with_capacity(trimmed.len());
    let mut prev_space = false;
    let mut chars = trimmed.char_indices().peekable();

    while let Some((i, c)) = chars.next() {
        if c.is_whitespace() {
            // Check if next char is a significant HGVS character
            if let Some(&(_, next_c)) = chars.peek() {
                // Don't add space before these characters
                if matches!(
                    next_c,
                    '.' | ':' | '>' | '<' | '+' | '-' | '_' | '(' | ')' | '[' | ']'
                ) {
                    if corrections.is_empty()
                        || corrections.last().unwrap().error_type != ErrorType::ExtraWhitespace
                    {
                        corrections.push(DetectedCorrection::new(
                            ErrorType::ExtraWhitespace,
                            " ".to_string(),
                            "".to_string(),
                            i,
                            i + 1,
                        ));
                    }
                    prev_space = false;
                    continue;
                }
            }
            prev_space = true;
        } else {
            // Check if previous char was space and current is significant
            if prev_space
                && matches!(
                    c,
                    '.' | ':' | '>' | '<' | '+' | '-' | '_' | '(' | ')' | '[' | ']'
                )
            {
                // Don't add space before this character
                prev_space = false;
                result.push(c);
                continue;
            }

            if prev_space {
                result.push(' ');
                prev_space = false;
            }
            result.push(c);
        }
    }

    (result, corrections)
}

/// Check if the input contains position zero, which is never valid.
///
/// Returns Some with the position if found, None otherwise.
pub fn detect_position_zero(input: &str) -> Option<usize> {
    // Look for patterns like "c.0", "g.0", "p.0", etc.
    // Also check for ".0_" or ".0+" or ".0-" patterns
    let mut chars = input.char_indices().peekable();

    while let Some((i, c)) = chars.next() {
        if c == '.' {
            if let Some(&(_, next_c)) = chars.peek() {
                if next_c == '0' {
                    // Check if followed by non-digit or end of string
                    chars.next();
                    if let Some(&(_, after_zero)) = chars.peek() {
                        if !after_zero.is_ascii_digit() {
                            return Some(i + 1); // Position of the '0'
                        }
                    } else {
                        return Some(i + 1);
                    }
                }
            }
        }
    }
    None
}

/// Correct protein substitution arrow syntax.
///
/// Converts `p.Val600>Glu` to `p.Val600Glu`.
pub fn correct_protein_arrow(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();

    // Simple pattern: look for ">Xxx" where Xxx is an amino acid
    // This is a simplified implementation; full implementation would use the parser
    if let Some(arrow_pos) = input.find('>') {
        // Check if this looks like a protein variant (starts with p.)
        if input.starts_with("p.") || input.contains(":p.") {
            // Check if after '>' there's an amino acid-like pattern
            let after_arrow = &input[arrow_pos + 1..];
            if after_arrow.len() >= 3 {
                let potential_aa = &after_arrow[..3];
                if is_amino_acid_like(potential_aa) {
                    let corrected = format!("{}{}", &input[..arrow_pos], &input[arrow_pos + 1..]);
                    corrections.push(DetectedCorrection::new(
                        ErrorType::ProteinSubstitutionArrow,
                        ">",
                        "",
                        arrow_pos,
                        arrow_pos + 1,
                    ));
                    return (corrected, corrections);
                }
            }
        }
    }

    (input.to_string(), corrections)
}

/// Check if a string looks like a 3-letter amino acid code.
fn is_amino_acid_like(s: &str) -> bool {
    if s.len() < 3 {
        return false;
    }
    let bytes = s.as_bytes();
    bytes[0].is_ascii_alphabetic()
        && bytes[1].is_ascii_alphabetic()
        && bytes[2].is_ascii_alphabetic()
}

/// Correct lowercase amino acid codes to proper case.
///
/// Converts `val` to `Val`, `GLU` to `Glu`, etc.
///
/// Note: This function will be used by the smart error correction feature.
#[allow(dead_code)]
pub fn correct_amino_acid_case(token: &str) -> Option<(String, ErrorType)> {
    // Three-letter amino acid codes (case-insensitive check)
    let amino_acids = [
        ("ala", "Ala"),
        ("arg", "Arg"),
        ("asn", "Asn"),
        ("asp", "Asp"),
        ("cys", "Cys"),
        ("gln", "Gln"),
        ("glu", "Glu"),
        ("gly", "Gly"),
        ("his", "His"),
        ("ile", "Ile"),
        ("leu", "Leu"),
        ("lys", "Lys"),
        ("met", "Met"),
        ("phe", "Phe"),
        ("pro", "Pro"),
        ("sec", "Sec"),
        ("ser", "Ser"),
        ("thr", "Thr"),
        ("trp", "Trp"),
        ("tyr", "Tyr"),
        ("val", "Val"),
        ("ter", "Ter"),
        ("xaa", "Xaa"),
    ];

    let lower = token.to_lowercase();
    for (pattern, correct) in amino_acids {
        if lower == pattern && token != correct {
            return Some((correct.to_string(), ErrorType::LowercaseAminoAcid));
        }
    }

    None
}

/// Convert single-letter amino acid code to three-letter code.
///
/// Note: This function will be used by the smart error correction feature.
#[allow(dead_code)]
pub fn single_to_three_letter_aa(single: char) -> Option<&'static str> {
    match single.to_ascii_uppercase() {
        'A' => Some("Ala"),
        'R' => Some("Arg"),
        'N' => Some("Asn"),
        'D' => Some("Asp"),
        'C' => Some("Cys"),
        'Q' => Some("Gln"),
        'E' => Some("Glu"),
        'G' => Some("Gly"),
        'H' => Some("His"),
        'I' => Some("Ile"),
        'L' => Some("Leu"),
        'K' => Some("Lys"),
        'M' => Some("Met"),
        'F' => Some("Phe"),
        'P' => Some("Pro"),
        'U' => Some("Sec"),
        'S' => Some("Ser"),
        'T' => Some("Thr"),
        'W' => Some("Trp"),
        'Y' => Some("Tyr"),
        'V' => Some("Val"),
        '*' => Some("Ter"),
        'X' => Some("Xaa"),
        _ => None,
    }
}

/// Correct lowercase accession prefix.
///
/// Converts `nm_000088.3` to `NM_000088.3`.
pub fn correct_accession_prefix_case(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();

    // Known prefixes that should be uppercase
    let prefixes = [
        ("nm_", "NM_"),
        ("np_", "NP_"),
        ("nc_", "NC_"),
        ("ng_", "NG_"),
        ("nr_", "NR_"),
        ("xm_", "XM_"),
        ("xp_", "XP_"),
        ("xr_", "XR_"),
        ("enst", "ENST"),
        ("ensp", "ENSP"),
        ("lrg_", "LRG_"),
    ];

    for (lower, correct) in prefixes {
        if input.to_lowercase().starts_with(lower) && !input.starts_with(correct) {
            let corrected = format!("{}{}", correct, &input[lower.len()..]);
            corrections.push(DetectedCorrection::new(
                ErrorType::LowercaseAccessionPrefix,
                input[..lower.len()].to_string(),
                correct.to_string(),
                0,
                lower.len(),
            ));
            return (corrected, corrections);
        }
    }

    (input.to_string(), corrections)
}

/// Correct mixed case edit types.
///
/// Converts `Del`, `INS`, `DELINS` to `del`, `ins`, `delins`.
///
/// Note: This function will be used by the smart error correction feature.
#[allow(dead_code)]
pub fn correct_edit_type_case(token: &str) -> Option<(String, ErrorType)> {
    let edit_types = [
        ("del", "del"),
        ("ins", "ins"),
        ("dup", "dup"),
        ("inv", "inv"),
        ("con", "con"),
        ("delins", "delins"),
    ];

    let lower = token.to_lowercase();
    for (pattern, correct) in edit_types {
        if lower == pattern && token != correct {
            return Some((correct.to_string(), ErrorType::MixedCaseEditType));
        }
    }

    None
}

/// Detect missing accession version.
///
/// Looks for accession patterns like `NM_000088` without a version suffix.
/// Returns Some with position and the accession if found.
pub fn detect_missing_version(input: &str) -> Option<(usize, String)> {
    // Match patterns like NM_000088 (without .N version)
    // Common prefixes: NM_, NP_, NC_, NG_, NR_, XM_, XP_, XR_, ENST, ENSP, LRG_
    let prefixes = [
        ("NM_", false),
        ("NP_", false),
        ("NC_", false),
        ("NG_", false),
        ("NR_", false),
        ("XM_", false),
        ("XP_", false),
        ("XR_", false),
        ("ENST", true), // Ensembl doesn't always have versions
        ("ENSP", true),
        ("LRG_", false),
    ];

    for (prefix, optional) in prefixes {
        if let Some(start) = input.find(prefix) {
            // Find the end of the accession (colon or end of string)
            let after_prefix = &input[start + prefix.len()..];
            let end = after_prefix.find(':').unwrap_or(after_prefix.len());
            let accession_body = &after_prefix[..end];

            // Check if it has a version (contains a period followed by digits)
            let has_version = accession_body.contains('.')
                && accession_body
                    .rsplit('.')
                    .next()
                    .map(|v| v.chars().all(|c| c.is_ascii_digit()))
                    .unwrap_or(false);

            // If no version and not optional, report it
            if !has_version && !optional {
                let full_accession = format!("{}{}", prefix, accession_body);
                return Some((start, full_accession));
            }
        }
    }
    None
}

/// Detect swapped interval positions.
///
/// Looks for patterns like `c.200_100del` where start > end.
/// Returns Some with the detected swap information and a suggested correction.
pub fn detect_swapped_positions(input: &str) -> Option<DetectedCorrection> {
    // Look for interval patterns: number_number
    // Pattern: .123_456 (position after coordinate type marker)

    // Find coordinate marker
    let coord_markers = [".c.", ".g.", ".n.", ".m.", ".o.", ".r."];
    let mut search_start = 0;

    for marker in &coord_markers {
        if let Some(pos) = input.find(marker) {
            search_start = pos + marker.len();
            break;
        }
    }

    // Also check without accession: "c.100_200del"
    if search_start == 0 {
        if let Some(pos) = input.find(['c', 'g', 'n', 'm', 'r']) {
            if input.get(pos + 1..pos + 2) == Some(".") {
                search_start = pos + 2;
            }
        }
    }

    if search_start == 0 {
        return None;
    }

    let after_marker = &input[search_start..];

    // Parse first number (may have offset like 100+5)
    let mut chars = after_marker.char_indices().peekable();
    let mut first_num_str = String::new();

    // Handle negative positions (like c.-100)
    if let Some(&(_, '-')) = chars.peek() {
        first_num_str.push('-');
        chars.next();
    }

    // Collect digits for first position
    while let Some(&(_, c)) = chars.peek() {
        if c.is_ascii_digit() {
            first_num_str.push(c);
            chars.next();
        } else {
            break;
        }
    }

    if first_num_str.is_empty() || first_num_str == "-" {
        return None;
    }

    // Skip any offset (like +5 or -10)
    while let Some(&(_, c)) = chars.peek() {
        if c == '+' || c == '-' || c.is_ascii_digit() {
            chars.next();
        } else {
            break;
        }
    }

    // Check for underscore (interval marker)
    if chars.next().map(|(_, c)| c) != Some('_') {
        return None;
    }

    let _second_start = chars.peek().map(|(i, _)| *i).unwrap_or(0);
    let mut second_num_str = String::new();

    // Handle negative positions for second number
    if let Some(&(_, '-')) = chars.peek() {
        second_num_str.push('-');
        chars.next();
    }

    // Collect digits for second position
    while let Some(&(_, c)) = chars.peek() {
        if c.is_ascii_digit() {
            second_num_str.push(c);
            chars.next();
        } else {
            break;
        }
    }

    if second_num_str.is_empty() || second_num_str == "-" {
        return None;
    }

    // Parse the numbers
    let first_num: i64 = first_num_str.parse().ok()?;
    let second_num: i64 = second_num_str.parse().ok()?;

    // Check if swapped (start > end)
    if first_num > second_num {
        // Construct corrected string
        let prefix = &input[..search_start];
        let suffix_start =
            search_start + chars.peek().map(|(i, _)| *i).unwrap_or(after_marker.len());
        let suffix = &input[suffix_start..];

        let _corrected = format!("{}{}_{}{}", prefix, second_num_str, first_num_str, suffix);

        return Some(DetectedCorrection::new(
            ErrorType::SwappedPositions,
            format!("{}_{}", first_num_str, second_num_str),
            format!("{}_{}", second_num_str, first_num_str),
            search_start,
            suffix_start,
        ));
    }

    None
}

/// Detect and strip trailing protein annotation.
///
/// ClinVar and other databases often append protein consequence annotations
/// to HGVS expressions, like `NM_003467.3:c.708G>A (p.Lys236=)`.
/// This function strips the trailing annotation.
///
/// Returns the stripped string and a list of corrections made.
pub fn strip_trailing_annotation(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();

    // Look for patterns like " (p.XXX)" or "(p.XXX)" at the end
    // Pattern: optional whitespace, opening paren, 'p.', content, closing paren, optional whitespace, end
    let trimmed = input.trim_end();

    // Find the last occurrence of '(p.' to handle cases with parentheses elsewhere
    if let Some(start_idx) = trimmed.rfind("(p.") {
        // Check that this is followed by a closing paren
        if let Some(close_idx) = trimmed[start_idx..].find(')') {
            let abs_close_idx = start_idx + close_idx;

            // Verify there's nothing significant after the closing paren
            let after_close = trimmed[abs_close_idx + 1..].trim();
            if after_close.is_empty() {
                // Check for leading whitespace before the annotation
                let before_annotation = &trimmed[..start_idx];
                let stripped = before_annotation.trim_end();

                // The annotation includes any whitespace before it
                let annotation_start = stripped.len();
                let annotation = &input[annotation_start..];

                corrections.push(DetectedCorrection::new(
                    ErrorType::TrailingAnnotation,
                    annotation.trim(),
                    "",
                    annotation_start,
                    input.len(),
                ));

                return (stripped.to_string(), corrections);
            }
        }
    }

    (input.to_string(), corrections)
}

/// Infer and add missing coordinate type prefix for genomic accessions.
///
/// Converts patterns like `NC_000001.11:12345A>G` to `NC_000001.11:g.12345A>G`.
/// Also handles uncertain boundaries like `NC_000001.11:(?_pos)_(pos_?)del`.
///
/// Returns the corrected string and a list of corrections made.
pub fn correct_missing_coordinate_prefix(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();

    // Find the colon that separates accession from variant
    let colon_pos = match input.find(':') {
        Some(pos) => pos,
        None => return (input.to_string(), corrections),
    };

    let accession = &input[..colon_pos];
    let after_colon = &input[colon_pos + 1..];

    // Check if this is a genomic accession (NC_, NG_, NT_, NW_, or LRG_ without t suffix)
    let is_genomic_accession = accession.starts_with("NC_")
        || accession.starts_with("NG_")
        || accession.starts_with("NT_")
        || accession.starts_with("NW_")
        || (accession.starts_with("LRG_") && !accession.contains('t'));

    if !is_genomic_accession {
        return (input.to_string(), corrections);
    }

    // Check if already has a coordinate prefix (g., c., p., n., r., m.)
    let has_prefix = after_colon.starts_with("g.")
        || after_colon.starts_with("c.")
        || after_colon.starts_with("p.")
        || after_colon.starts_with("n.")
        || after_colon.starts_with("r.")
        || after_colon.starts_with("m.");

    if has_prefix {
        return (input.to_string(), corrections);
    }

    // Check if this looks like a variant (starts with digit, paren, or ?)
    let first_char = after_colon.chars().next().unwrap_or(' ');
    if !first_char.is_ascii_digit() && first_char != '(' && first_char != '?' && first_char != '[' {
        return (input.to_string(), corrections);
    }

    // Add g. prefix
    let corrected = format!("{}:g.{}", accession, after_colon);
    corrections.push(DetectedCorrection::new(
        ErrorType::MissingCoordinatePrefix,
        format!("{}:", accession),
        format!("{}:g.", accession),
        0,
        colon_pos + 1,
    ));

    (corrected, corrections)
}

/// Correct swapped interval positions.
///
/// Swaps `c.200_100del` to `c.100_200del`.
pub fn correct_swapped_positions(input: &str) -> (String, Vec<DetectedCorrection>) {
    if let Some(correction) = detect_swapped_positions(input) {
        let corrected = format!(
            "{}{}{}",
            &input[..correction.start],
            correction.corrected,
            &input[correction.end..]
        );
        (corrected, vec![correction])
    } else {
        (input.to_string(), vec![])
    }
}

/// Correct old/deprecated allele format to current HGVS standard.
///
/// Converts patterns like `NM_000088.3:[c.100A>G;c.200C>T]` (coordinate type inside brackets)
/// to `NM_000088.3:c.[100A>G;200C>T]` (coordinate type before brackets).
///
/// Returns the corrected string and a list of corrections made.
pub fn correct_old_allele_format(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();

    // Pattern to detect: ref:[x.edit1;x.edit2;...] where x is the coordinate type
    // We need to find :[c. or :[g. or :[n. or :[r. or :[m. or :[p. followed by variants with same prefix

    // Find the colon followed by opening bracket
    let colon_bracket = match input.find(":[") {
        Some(pos) => pos,
        None => return (input.to_string(), corrections),
    };

    // Check what's after the bracket
    let after_bracket = &input[colon_bracket + 2..];

    // Determine the coordinate type by looking at the first variant
    let coord_type = if after_bracket.starts_with("c.") {
        "c."
    } else if after_bracket.starts_with("g.") {
        "g."
    } else if after_bracket.starts_with("n.") {
        "n."
    } else if after_bracket.starts_with("r.") {
        "r."
    } else if after_bracket.starts_with("m.") {
        "m."
    } else if after_bracket.starts_with("p.") {
        "p."
    } else {
        // Not an old allele format we can correct
        return (input.to_string(), corrections);
    };

    // Find the closing bracket
    let close_bracket = match input[colon_bracket..].find(']') {
        Some(pos) => colon_bracket + pos,
        None => return (input.to_string(), corrections),
    };

    // Extract the content between brackets
    let content = &input[colon_bracket + 2..close_bracket];

    // Check if all variants have the same coordinate prefix
    let parts: Vec<&str> = content.split(';').collect();
    if parts.is_empty() {
        return (input.to_string(), corrections);
    }

    // Verify all parts start with the same coordinate type
    let all_same_prefix = parts.iter().all(|p| p.trim().starts_with(coord_type));
    if !all_same_prefix {
        return (input.to_string(), corrections);
    }

    // Strip the coordinate prefix from each part
    let stripped_parts: Vec<String> = parts
        .iter()
        .map(|p| {
            let trimmed = p.trim();
            if let Some(stripped) = trimmed.strip_prefix(coord_type) {
                stripped.to_string()
            } else {
                trimmed.to_string()
            }
        })
        .collect();

    // Build the corrected string: ref:x.[edit1;edit2;...]
    let accession = &input[..colon_bracket];
    let new_content = stripped_parts.join(";");
    let remaining = &input[close_bracket + 1..];
    let corrected = format!("{}:{}[{}]{}", accession, coord_type, new_content, remaining);

    let original_pattern = &input[colon_bracket..close_bracket + 1];
    let corrected_pattern = format!(":{}[{}]", coord_type, new_content);

    corrections.push(DetectedCorrection::new(
        ErrorType::OldAlleleFormat,
        original_pattern,
        &corrected_pattern,
        colon_bracket,
        close_bracket + 1,
    ));

    (corrected, corrections)
}

/// Apply a correction based on the resolved action.
///
/// Returns the corrected input and whether a warning should be emitted.
///
/// Note: This function will be used by the smart error correction feature.
#[allow(dead_code)]
pub fn apply_correction(
    original: &str,
    correction: &DetectedCorrection,
    action: ResolvedAction,
) -> (String, bool) {
    match action {
        ResolvedAction::Reject => (original.to_string(), false),
        ResolvedAction::WarnCorrect => {
            let corrected = format!(
                "{}{}{}",
                &original[..correction.start],
                correction.corrected,
                &original[correction.end..]
            );
            (corrected, true)
        }
        ResolvedAction::SilentCorrect => {
            let corrected = format!(
                "{}{}{}",
                &original[..correction.start],
                correction.corrected,
                &original[correction.end..]
            );
            (corrected, false)
        }
        ResolvedAction::Accept => (original.to_string(), false),
    }
}

// =============================================================================
// Fuzzy Matching (Levenshtein Distance)
// =============================================================================

/// Calculate the Levenshtein (edit) distance between two strings.
///
/// The Levenshtein distance is the minimum number of single-character edits
/// (insertions, deletions, or substitutions) required to change one string
/// into the other.
///
/// # Examples
/// ```
/// use ferro_hgvs::error_handling::corrections::levenshtein_distance;
///
/// assert_eq!(levenshtein_distance("kitten", "sitting"), 3);
/// assert_eq!(levenshtein_distance("del", "del"), 0);
/// assert_eq!(levenshtein_distance("del", "dek"), 1);
/// ```
pub fn levenshtein_distance(a: &str, b: &str) -> usize {
    let a_chars: Vec<char> = a.chars().collect();
    let b_chars: Vec<char> = b.chars().collect();
    let a_len = a_chars.len();
    let b_len = b_chars.len();

    // Early returns for edge cases
    if a_len == 0 {
        return b_len;
    }
    if b_len == 0 {
        return a_len;
    }

    // Use a single row for space efficiency
    let mut prev_row: Vec<usize> = (0..=b_len).collect();
    let mut curr_row: Vec<usize> = vec![0; b_len + 1];

    for (i, a_char) in a_chars.iter().enumerate() {
        curr_row[0] = i + 1;

        for (j, b_char) in b_chars.iter().enumerate() {
            let cost = if a_char == b_char { 0 } else { 1 };
            curr_row[j + 1] = min(
                min(
                    prev_row[j + 1] + 1, // deletion
                    curr_row[j] + 1,     // insertion
                ),
                prev_row[j] + cost, // substitution
            );
        }

        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row[b_len]
}

/// A fuzzy match result with the matched string and its edit distance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FuzzyMatch {
    /// The matched string.
    pub matched: String,
    /// The Levenshtein distance from the input.
    pub distance: usize,
}

impl FuzzyMatch {
    /// Create a new fuzzy match.
    pub fn new(matched: impl Into<String>, distance: usize) -> Self {
        Self {
            matched: matched.into(),
            distance,
        }
    }
}

/// Known accession prefixes for fuzzy matching.
const ACCESSION_PREFIXES: &[&str] = &[
    "NM_", "NP_", "NC_", "NG_", "NR_", "NW_", "NT_", "XM_", "XP_", "XR_", "ENST", "ENSG", "ENSP",
    "ENSE", "LRG_",
];

/// Known edit types for fuzzy matching.
const EDIT_TYPES: &[&str] = &["del", "ins", "dup", "inv", "con", "delins", "fs", "ext"];

/// Known amino acid codes (three-letter) for fuzzy matching.
const AMINO_ACIDS: &[&str] = &[
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
    "Phe", "Pro", "Sec", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter", "Xaa",
];

/// Find the closest match in a list of candidates.
///
/// Returns the best match if the distance is within the threshold.
///
/// # Arguments
/// * `input` - The input string to match
/// * `candidates` - List of valid strings to match against
/// * `max_distance` - Maximum edit distance for a valid match
/// * `case_sensitive` - Whether to perform case-sensitive matching
///
/// # Examples
/// ```
/// use ferro_hgvs::error_handling::corrections::find_closest_match;
///
/// let candidates = vec!["del", "ins", "dup"];
/// let result = find_closest_match("dek", &candidates, 1, true);
/// assert!(result.is_some());
/// assert_eq!(result.unwrap().matched, "del");
/// ```
pub fn find_closest_match(
    input: &str,
    candidates: &[&str],
    max_distance: usize,
    case_sensitive: bool,
) -> Option<FuzzyMatch> {
    let input_normalized = if case_sensitive {
        input.to_string()
    } else {
        input.to_lowercase()
    };

    let mut best_match: Option<FuzzyMatch> = None;

    for &candidate in candidates {
        let candidate_normalized = if case_sensitive {
            candidate.to_string()
        } else {
            candidate.to_lowercase()
        };

        let distance = levenshtein_distance(&input_normalized, &candidate_normalized);

        if distance <= max_distance
            && (best_match.is_none() || distance < best_match.as_ref().unwrap().distance)
        {
            best_match = Some(FuzzyMatch::new(candidate, distance));
        }
    }

    best_match
}

/// Detect a typo in an accession prefix.
///
/// Returns the suggested correction if a close match is found.
///
/// # Examples
/// ```
/// use ferro_hgvs::error_handling::corrections::detect_accession_typo;
///
/// // Transposed letters
/// let result = detect_accession_typo("MN_000088.3");
/// assert!(result.is_some());
/// assert_eq!(result.unwrap().matched, "NM_");
///
/// // Correct prefix - no suggestion
/// let result = detect_accession_typo("NM_000088.3");
/// assert!(result.is_none());
/// ```
pub fn detect_accession_typo(input: &str) -> Option<FuzzyMatch> {
    // Check if the input starts with any known prefix (case-insensitive)
    let upper = input.to_uppercase();
    for prefix in ACCESSION_PREFIXES {
        if upper.starts_with(prefix) {
            return None; // Valid prefix, no typo
        }
    }

    // Extract potential prefix (first 2-4 characters + optional underscore)
    let prefix_end = input
        .char_indices()
        .take(5)
        .find(|(_, c)| c.is_ascii_digit())
        .map(|(i, _)| i)
        .unwrap_or_else(|| input.len().min(5));

    if prefix_end < 2 {
        return None;
    }

    let potential_prefix = &input[..prefix_end];

    // Find closest match
    find_closest_match(potential_prefix, ACCESSION_PREFIXES, 2, false).filter(|m| m.distance > 0)
    // Only return if there's actually a difference
}

/// Detect a typo in an edit type keyword.
///
/// Returns the suggested correction if a close match is found.
///
/// # Examples
/// ```
/// use ferro_hgvs::error_handling::corrections::detect_edit_type_typo;
///
/// // Typo: 'dek' instead of 'del'
/// let result = detect_edit_type_typo("dek");
/// assert!(result.is_some());
/// assert_eq!(result.unwrap().matched, "del");
///
/// // Correct - no suggestion
/// let result = detect_edit_type_typo("del");
/// assert!(result.is_none());
/// ```
pub fn detect_edit_type_typo(input: &str) -> Option<FuzzyMatch> {
    // Check if already a valid edit type
    let lower = input.to_lowercase();
    if EDIT_TYPES.contains(&lower.as_str()) {
        return None;
    }

    // Find closest match (max distance 1 for short keywords)
    find_closest_match(input, EDIT_TYPES, 1, false)
}

/// Detect a typo in an amino acid code.
///
/// Returns the suggested correction if a close match is found.
///
/// # Examples
/// ```
/// use ferro_hgvs::error_handling::corrections::detect_amino_acid_typo;
///
/// // Typo: 'Vasl' instead of 'Val'
/// let result = detect_amino_acid_typo("Vasl");
/// assert!(result.is_some());
/// assert_eq!(result.unwrap().matched, "Val");
///
/// // Correct - no suggestion
/// let result = detect_amino_acid_typo("Val");
/// assert!(result.is_none());
/// ```
pub fn detect_amino_acid_typo(input: &str) -> Option<FuzzyMatch> {
    // Check if already a valid amino acid (case-insensitive)
    for aa in AMINO_ACIDS {
        if input.eq_ignore_ascii_case(aa) {
            return None;
        }
    }

    // Find closest match (max distance 1 for short codes)
    find_closest_match(input, AMINO_ACIDS, 1, false)
}

/// A detected typo with its context and suggested correction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TypoSuggestion {
    /// The original text that appears to be a typo.
    pub original: String,
    /// The suggested correction.
    pub suggestion: String,
    /// The type of token (prefix, edit, amino acid).
    pub token_type: TypoTokenType,
    /// Start position in the input string.
    pub start: usize,
    /// End position in the input string.
    pub end: usize,
    /// Edit distance from original to suggestion.
    pub distance: usize,
}

/// Type of token that contains a typo.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypoTokenType {
    /// Accession prefix (NM_, NP_, etc.)
    AccessionPrefix,
    /// Edit type keyword (del, ins, dup, etc.)
    EditType,
    /// Amino acid code (Val, Glu, etc.)
    AminoAcid,
}

impl std::fmt::Display for TypoTokenType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TypoTokenType::AccessionPrefix => write!(f, "accession prefix"),
            TypoTokenType::EditType => write!(f, "edit type"),
            TypoTokenType::AminoAcid => write!(f, "amino acid"),
        }
    }
}

impl TypoSuggestion {
    /// Format a suggestion message for this typo.
    pub fn suggestion_message(&self) -> String {
        format!(
            "possible typo in {}: '{}' → did you mean '{}'?",
            self.token_type, self.original, self.suggestion
        )
    }
}

/// Detect typos in an HGVS input string.
///
/// This function scans the input for potential typos in:
/// - Accession prefixes (NM_, NP_, etc.)
/// - Edit type keywords (del, ins, dup, etc.)
/// - Amino acid codes (Val, Glu, etc.)
///
/// # Examples
/// ```
/// use ferro_hgvs::error_handling::corrections::detect_typos;
///
/// // Typo in accession prefix
/// let typos = detect_typos("MN_000088.3:c.100A>G");
/// assert!(!typos.is_empty());
/// assert_eq!(typos[0].suggestion, "NM_");
///
/// // No typos
/// let typos = detect_typos("NM_000088.3:c.100A>G");
/// assert!(typos.is_empty());
/// ```
pub fn detect_typos(input: &str) -> Vec<TypoSuggestion> {
    let mut suggestions = Vec::new();

    // Check for accession prefix typos
    if let Some(fuzzy) = detect_accession_typo(input) {
        // Find where the prefix ends
        let prefix_end = input
            .char_indices()
            .take(10)
            .find(|(_, c)| c.is_ascii_digit())
            .map(|(i, _)| i)
            .unwrap_or(fuzzy.matched.len());

        suggestions.push(TypoSuggestion {
            original: input[..prefix_end].to_string(),
            suggestion: fuzzy.matched,
            token_type: TypoTokenType::AccessionPrefix,
            start: 0,
            end: prefix_end,
            distance: fuzzy.distance,
        });
    }

    // Check for edit type typos (look for keywords after position numbers)
    // Pattern: after digits or after '_', look for 3-7 letter sequences
    let edit_pattern = regex_lite_find_edit_types(input);
    for (start, end, token) in edit_pattern {
        if let Some(fuzzy) = detect_edit_type_typo(&token) {
            suggestions.push(TypoSuggestion {
                original: token,
                suggestion: fuzzy.matched,
                token_type: TypoTokenType::EditType,
                start,
                end,
                distance: fuzzy.distance,
            });
        }
    }

    // Check for amino acid typos (look for patterns after p. or in protein variants)
    if input.contains("p.") || input.contains(":p.") {
        let aa_pattern = find_potential_amino_acids(input);
        for (start, end, token) in aa_pattern {
            if let Some(fuzzy) = detect_amino_acid_typo(&token) {
                suggestions.push(TypoSuggestion {
                    original: token,
                    suggestion: fuzzy.matched,
                    token_type: TypoTokenType::AminoAcid,
                    start,
                    end,
                    distance: fuzzy.distance,
                });
            }
        }
    }

    suggestions
}

/// Find potential edit type keywords in the input (simple pattern matching).
fn regex_lite_find_edit_types(input: &str) -> Vec<(usize, usize, String)> {
    let mut results = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        // Look for positions after digits where edit types might appear
        if chars[i].is_ascii_digit() {
            // Skip to end of number
            while i < len && (chars[i].is_ascii_digit() || chars[i] == '+' || chars[i] == '-') {
                i += 1;
            }

            // Collect alphabetic sequence
            let start = i;
            let mut token = String::new();
            while i < len && chars[i].is_ascii_alphabetic() {
                token.push(chars[i]);
                i += 1;
            }

            // Check if it's a potential edit type (2-7 letters)
            if token.len() >= 2 && token.len() <= 7 {
                // Skip known valid patterns
                let lower = token.to_lowercase();
                if !lower.starts_with("ins")
                    && !EDIT_TYPES.contains(&lower.as_str())
                    && !AMINO_ACIDS.iter().any(|aa| aa.eq_ignore_ascii_case(&token))
                {
                    results.push((start, i, token));
                }
            }
        } else {
            i += 1;
        }
    }

    results
}

/// Find potential amino acid codes in protein variants.
fn find_potential_amino_acids(input: &str) -> Vec<(usize, usize, String)> {
    let mut results = Vec::new();

    // Find the p. prefix
    let p_pos = if let Some(pos) = input.find(":p.") {
        pos + 3
    } else if let Some(pos) = input.find("p.") {
        pos + 2
    } else {
        return results;
    };

    let chars: Vec<char> = input[p_pos..].chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        // Look for sequences that start with uppercase letter
        if chars[i].is_ascii_uppercase() {
            let start = p_pos + i;
            let mut token = String::new();
            token.push(chars[i]);
            i += 1;

            // Collect remaining letters (lowercase typically for amino acids)
            while i < len && chars[i].is_ascii_alphabetic() && !chars[i].is_ascii_uppercase() {
                token.push(chars[i]);
                i += 1;
            }

            // Check if it looks like an amino acid (3-4 letters)
            if token.len() >= 3 && token.len() <= 4 {
                let end = p_pos + i;
                results.push((start, end, token));
            }
        } else {
            i += 1;
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    // Dash character correction tests
    #[test]
    fn test_correct_dash_en_dash() {
        let (corrected, corrections) = correct_dash_characters("c.100\u{2013}200del");
        assert_eq!(corrected, "c.100-200del");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::WrongDashCharacter);
    }

    #[test]
    fn test_correct_dash_em_dash() {
        let (corrected, corrections) = correct_dash_characters("c.100\u{2014}200del");
        assert_eq!(corrected, "c.100-200del");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_dash_minus_sign() {
        let (corrected, corrections) = correct_dash_characters("c.100\u{2212}200del");
        assert_eq!(corrected, "c.100-200del");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_dash_no_change() {
        let (corrected, corrections) = correct_dash_characters("c.100-200del");
        assert_eq!(corrected, "c.100-200del");
        assert!(corrections.is_empty());
    }

    // Quote character correction tests
    #[test]
    fn test_correct_smart_single_quotes() {
        let (corrected, corrections) = correct_quote_characters("c.100ins\u{2018}A\u{2019}");
        assert_eq!(corrected, "c.100ins'A'");
        assert_eq!(corrections.len(), 2);
    }

    #[test]
    fn test_correct_smart_double_quotes() {
        let (corrected, corrections) = correct_quote_characters("c.100ins\u{201C}ATG\u{201D}");
        assert_eq!(corrected, "c.100ins\"ATG\"");
        assert_eq!(corrections.len(), 2);
    }

    // Whitespace correction tests
    #[test]
    fn test_correct_whitespace_trim() {
        let (corrected, corrections) = correct_whitespace("  c.100A>G  ");
        assert_eq!(corrected, "c.100A>G");
        assert!(!corrections.is_empty());
    }

    #[test]
    fn test_correct_whitespace_no_change() {
        let (corrected, _corrections) = correct_whitespace("c.100A>G");
        assert_eq!(corrected, "c.100A>G");
        // May have corrections detected but result should be same
        assert_eq!(corrected, "c.100A>G");
    }

    // Position zero detection tests
    #[test]
    fn test_detect_position_zero() {
        assert!(detect_position_zero("c.0A>G").is_some());
        assert!(detect_position_zero("g.0del").is_some());
        assert!(detect_position_zero("c.10A>G").is_none());
        assert!(detect_position_zero("c.100A>G").is_none());
    }

    // Protein arrow correction tests
    #[test]
    fn test_correct_protein_arrow() {
        let (corrected, corrections) = correct_protein_arrow("p.Val600>Glu");
        assert_eq!(corrected, "p.Val600Glu");
        assert_eq!(corrections.len(), 1);
        assert_eq!(
            corrections[0].error_type,
            ErrorType::ProteinSubstitutionArrow
        );
    }

    #[test]
    fn test_correct_protein_arrow_no_change() {
        let (corrected, corrections) = correct_protein_arrow("p.Val600Glu");
        assert_eq!(corrected, "p.Val600Glu");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_protein_arrow_not_protein() {
        // DNA variant with > should not be modified
        let (corrected, corrections) = correct_protein_arrow("c.100A>G");
        assert_eq!(corrected, "c.100A>G");
        assert!(corrections.is_empty());
    }

    // Amino acid case correction tests
    #[test]
    fn test_correct_amino_acid_case_lowercase() {
        let result = correct_amino_acid_case("val");
        assert!(result.is_some());
        let (corrected, error_type) = result.unwrap();
        assert_eq!(corrected, "Val");
        assert_eq!(error_type, ErrorType::LowercaseAminoAcid);
    }

    #[test]
    fn test_correct_amino_acid_case_uppercase() {
        let result = correct_amino_acid_case("VAL");
        assert!(result.is_some());
        let (corrected, _) = result.unwrap();
        assert_eq!(corrected, "Val");
    }

    #[test]
    fn test_correct_amino_acid_case_correct() {
        let result = correct_amino_acid_case("Val");
        assert!(result.is_none());
    }

    // Single letter amino acid tests
    #[test]
    fn test_single_to_three_letter_aa() {
        assert_eq!(single_to_three_letter_aa('V'), Some("Val"));
        assert_eq!(single_to_three_letter_aa('E'), Some("Glu"));
        assert_eq!(single_to_three_letter_aa('*'), Some("Ter"));
        assert_eq!(single_to_three_letter_aa('Z'), None);
    }

    // Accession prefix case tests
    #[test]
    fn test_correct_accession_prefix_lowercase() {
        let (corrected, corrections) = correct_accession_prefix_case("nm_000088.3");
        assert_eq!(corrected, "NM_000088.3");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_accession_prefix_correct() {
        let (corrected, corrections) = correct_accession_prefix_case("NM_000088.3");
        assert_eq!(corrected, "NM_000088.3");
        assert!(corrections.is_empty());
    }

    // Edit type case tests
    #[test]
    fn test_correct_edit_type_uppercase() {
        let result = correct_edit_type_case("DEL");
        assert!(result.is_some());
        let (corrected, _) = result.unwrap();
        assert_eq!(corrected, "del");
    }

    #[test]
    fn test_correct_edit_type_mixed() {
        let result = correct_edit_type_case("Del");
        assert!(result.is_some());
        let (corrected, _) = result.unwrap();
        assert_eq!(corrected, "del");
    }

    #[test]
    fn test_correct_edit_type_correct() {
        let result = correct_edit_type_case("del");
        assert!(result.is_none());
    }

    // DetectedCorrection tests
    #[test]
    fn test_detected_correction_warning_message() {
        let correction =
            DetectedCorrection::new(ErrorType::WrongDashCharacter, "\u{2013}", "-", 5, 8);
        let msg = correction.warning_message();
        assert!(msg.contains("wrong dash character"));
        assert!(msg.contains("position 5"));
    }

    // Missing version detection tests
    #[test]
    fn test_detect_missing_version_no_version() {
        let result = detect_missing_version("NM_000088:c.100A>G");
        assert!(result.is_some());
        let (pos, acc) = result.unwrap();
        assert_eq!(pos, 0);
        assert_eq!(acc, "NM_000088");
    }

    #[test]
    fn test_detect_missing_version_with_version() {
        let result = detect_missing_version("NM_000088.3:c.100A>G");
        assert!(result.is_none());
    }

    #[test]
    fn test_detect_missing_version_ensembl() {
        // Ensembl accessions don't require versions
        let result = detect_missing_version("ENST00000123456:c.100A>G");
        assert!(result.is_none());
    }

    #[test]
    fn test_detect_missing_version_nc() {
        let result = detect_missing_version("NC_000001:g.12345A>G");
        assert!(result.is_some());
    }

    // Swapped position detection tests
    #[test]
    fn test_detect_swapped_positions_swapped() {
        let result = detect_swapped_positions("c.200_100del");
        assert!(result.is_some());
        let correction = result.unwrap();
        assert_eq!(correction.error_type, ErrorType::SwappedPositions);
        assert_eq!(correction.original, "200_100");
        assert_eq!(correction.corrected, "100_200");
    }

    #[test]
    fn test_detect_swapped_positions_correct() {
        let result = detect_swapped_positions("c.100_200del");
        assert!(result.is_none());
    }

    #[test]
    fn test_detect_swapped_positions_with_accession() {
        let result = detect_swapped_positions("NM_000088.3:c.500_100del");
        assert!(result.is_some());
        let correction = result.unwrap();
        assert_eq!(correction.original, "500_100");
        assert_eq!(correction.corrected, "100_500");
    }

    #[test]
    fn test_detect_swapped_positions_genomic() {
        let result = detect_swapped_positions("g.2000_1000del");
        assert!(result.is_some());
    }

    #[test]
    fn test_detect_swapped_positions_negative() {
        // Negative positions are valid in c. coordinates (5' UTR)
        let result = detect_swapped_positions("c.-10_-50del");
        // -10 > -50, so this should be detected as swapped
        assert!(result.is_some());
    }

    #[test]
    fn test_correct_swapped_positions() {
        let (corrected, corrections) = correct_swapped_positions("c.200_100del");
        assert_eq!(corrected, "c.100_200del");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_swapped_positions_no_change() {
        let (corrected, corrections) = correct_swapped_positions("c.100_200del");
        assert_eq!(corrected, "c.100_200del");
        assert!(corrections.is_empty());
    }

    // ==========================================================================
    // Fuzzy Matching (Levenshtein Distance) Tests
    // ==========================================================================

    #[test]
    fn test_levenshtein_distance_identical() {
        assert_eq!(levenshtein_distance("del", "del"), 0);
        assert_eq!(levenshtein_distance("", ""), 0);
    }

    #[test]
    fn test_levenshtein_distance_one_empty() {
        assert_eq!(levenshtein_distance("", "abc"), 3);
        assert_eq!(levenshtein_distance("abc", ""), 3);
    }

    #[test]
    fn test_levenshtein_distance_substitution() {
        assert_eq!(levenshtein_distance("del", "dek"), 1);
        assert_eq!(levenshtein_distance("cat", "bat"), 1);
    }

    #[test]
    fn test_levenshtein_distance_insertion() {
        assert_eq!(levenshtein_distance("del", "delx"), 1);
        assert_eq!(levenshtein_distance("abc", "abcd"), 1);
    }

    #[test]
    fn test_levenshtein_distance_deletion() {
        assert_eq!(levenshtein_distance("delx", "del"), 1);
        assert_eq!(levenshtein_distance("abcd", "abc"), 1);
    }

    #[test]
    fn test_levenshtein_distance_transposition() {
        // Transposition is 2 operations (delete + insert)
        assert_eq!(levenshtein_distance("ab", "ba"), 2);
        assert_eq!(levenshtein_distance("NM_", "MN_"), 2);
    }

    #[test]
    fn test_levenshtein_distance_classic_example() {
        // Classic example: kitten -> sitting
        assert_eq!(levenshtein_distance("kitten", "sitting"), 3);
    }

    #[test]
    fn test_find_closest_match_exact() {
        let candidates = vec!["del", "ins", "dup"];
        let result = find_closest_match("del", &candidates, 1, true);
        // Exact matches return with distance 0
        assert!(result.is_some());
        assert_eq!(result.unwrap().distance, 0);
    }

    #[test]
    fn test_find_closest_match_typo() {
        let candidates = vec!["del", "ins", "dup"];
        let result = find_closest_match("dek", &candidates, 1, true);
        assert!(result.is_some());
        let m = result.unwrap();
        assert_eq!(m.matched, "del");
        assert_eq!(m.distance, 1);
    }

    #[test]
    fn test_find_closest_match_case_insensitive() {
        let candidates = vec!["Val", "Glu", "Ala"];
        let result = find_closest_match("VAL", &candidates, 0, false);
        assert!(result.is_some());
        assert_eq!(result.unwrap().matched, "Val");
    }

    #[test]
    fn test_find_closest_match_no_match() {
        let candidates = vec!["del", "ins", "dup"];
        let result = find_closest_match("xyz", &candidates, 1, true);
        assert!(result.is_none());
    }

    #[test]
    fn test_detect_accession_typo_transposed() {
        // NM_ -> MN_ (transposition)
        let result = detect_accession_typo("MN_000088.3");
        assert!(result.is_some());
        let m = result.unwrap();
        assert_eq!(m.matched, "NM_");
    }

    #[test]
    fn test_detect_accession_typo_correct() {
        let result = detect_accession_typo("NM_000088.3");
        assert!(result.is_none());
    }

    #[test]
    fn test_detect_accession_typo_nc() {
        // CN_ has distance 2 from both NM_ and NC_ (transposition)
        // The function returns the first match found, which is NM_
        let result = detect_accession_typo("CN_000001.10");
        assert!(result.is_some());
        // Both NM_ and NC_ are distance 2, returns first found
        let matched = result.unwrap().matched;
        assert!(matched == "NM_" || matched == "NC_");
    }

    #[test]
    fn test_detect_accession_typo_ensembl() {
        let result = detect_accession_typo("ESNT00000123456");
        assert!(result.is_some());
        assert_eq!(result.unwrap().matched, "ENST");
    }

    #[test]
    fn test_detect_edit_type_typo_dek() {
        let result = detect_edit_type_typo("dek");
        assert!(result.is_some());
        assert_eq!(result.unwrap().matched, "del");
    }

    #[test]
    fn test_detect_edit_type_typo_inx() {
        let result = detect_edit_type_typo("inx");
        assert!(result.is_some());
        assert_eq!(result.unwrap().matched, "ins");
    }

    #[test]
    fn test_detect_edit_type_typo_correct() {
        assert!(detect_edit_type_typo("del").is_none());
        assert!(detect_edit_type_typo("ins").is_none());
        assert!(detect_edit_type_typo("dup").is_none());
    }

    #[test]
    fn test_detect_amino_acid_typo_vasl() {
        let result = detect_amino_acid_typo("Vasl");
        assert!(result.is_some());
        assert_eq!(result.unwrap().matched, "Val");
    }

    #[test]
    fn test_detect_amino_acid_typo_gul() {
        // "Gul" to "Glu" is distance 2 (transposition), which exceeds max distance of 1
        // So this should NOT match - the threshold is intentionally low to avoid false positives
        let result = detect_amino_acid_typo("Gul");
        assert!(result.is_none());
    }

    #[test]
    fn test_detect_amino_acid_typo_single_char_error() {
        // "Vak" to "Val" is distance 1 (single substitution)
        let result = detect_amino_acid_typo("Vak");
        assert!(result.is_some());
        assert_eq!(result.unwrap().matched, "Val");
    }

    #[test]
    fn test_detect_amino_acid_typo_correct() {
        assert!(detect_amino_acid_typo("Val").is_none());
        assert!(detect_amino_acid_typo("Glu").is_none());
        assert!(detect_amino_acid_typo("Ala").is_none());
    }

    #[test]
    fn test_detect_typos_accession_prefix() {
        let typos = detect_typos("MN_000088.3:c.100A>G");
        assert!(!typos.is_empty());
        assert_eq!(typos[0].token_type, TypoTokenType::AccessionPrefix);
        assert_eq!(typos[0].suggestion, "NM_");
    }

    #[test]
    fn test_detect_typos_no_typos() {
        let typos = detect_typos("NM_000088.3:c.100A>G");
        assert!(typos.is_empty());
    }

    #[test]
    fn test_detect_typos_edit_type() {
        // Use a full accession to avoid accession prefix detection
        let typos = detect_typos("NM_000088.3:c.100dek");
        assert!(!typos.is_empty());
        // Find the edit type suggestion
        let edit_typo = typos
            .iter()
            .find(|t| t.token_type == TypoTokenType::EditType);
        assert!(edit_typo.is_some());
        assert_eq!(edit_typo.unwrap().suggestion, "del");
    }

    #[test]
    fn test_detect_typos_amino_acid() {
        let typos = detect_typos("p.Vasl600Glu");
        assert!(!typos.is_empty());
        assert_eq!(typos[0].token_type, TypoTokenType::AminoAcid);
        assert_eq!(typos[0].suggestion, "Val");
    }

    #[test]
    fn test_typo_suggestion_message() {
        let suggestion = TypoSuggestion {
            original: "MN_".to_string(),
            suggestion: "NM_".to_string(),
            token_type: TypoTokenType::AccessionPrefix,
            start: 0,
            end: 3,
            distance: 2,
        };
        let msg = suggestion.suggestion_message();
        assert!(msg.contains("accession prefix"));
        assert!(msg.contains("MN_"));
        assert!(msg.contains("NM_"));
    }

    #[test]
    fn test_typo_token_type_display() {
        assert_eq!(
            format!("{}", TypoTokenType::AccessionPrefix),
            "accession prefix"
        );
        assert_eq!(format!("{}", TypoTokenType::EditType), "edit type");
        assert_eq!(format!("{}", TypoTokenType::AminoAcid), "amino acid");
    }

    #[test]
    fn test_fuzzy_match_struct() {
        let m = FuzzyMatch::new("del", 1);
        assert_eq!(m.matched, "del");
        assert_eq!(m.distance, 1);
    }

    // ==========================================================================
    // Trailing Annotation Tests
    // ==========================================================================

    #[test]
    fn test_strip_trailing_annotation_synonymous() {
        let (corrected, corrections) =
            strip_trailing_annotation("NM_003467.3(CXCR4):c.708G>A (p.Lys236=)");
        assert_eq!(corrected, "NM_003467.3(CXCR4):c.708G>A");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::TrailingAnnotation);
    }

    #[test]
    fn test_strip_trailing_annotation_missense() {
        let (corrected, corrections) =
            strip_trailing_annotation("NM_021831.6(AGBL5):c.2083G>A (p.Val695Ile)");
        assert_eq!(corrected, "NM_021831.6(AGBL5):c.2083G>A");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_strip_trailing_annotation_frameshift() {
        let (corrected, corrections) =
            strip_trailing_annotation("NM_178127.5(ANGPTL5):c.1097dup (p.Asn366fs)");
        assert_eq!(corrected, "NM_178127.5(ANGPTL5):c.1097dup");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_strip_trailing_annotation_no_space() {
        let (corrected, corrections) = strip_trailing_annotation("NM_000088.3:c.459A>G(p.Lys153=)");
        assert_eq!(corrected, "NM_000088.3:c.459A>G");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_strip_trailing_annotation_no_annotation() {
        let (corrected, corrections) = strip_trailing_annotation("NM_000088.3:c.459A>G");
        assert_eq!(corrected, "NM_000088.3:c.459A>G");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_strip_trailing_annotation_valid_uncertain() {
        // This has parentheses but they're part of valid HGVS syntax, not annotation
        // The p. inside is valid predicted protein change notation
        let (corrected, corrections) = strip_trailing_annotation("NP_000079.2:p.(Val600Glu)");
        // This should NOT be stripped because it's valid HGVS
        assert_eq!(corrected, "NP_000079.2:p.(Val600Glu)");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_strip_trailing_annotation_protein_variant() {
        // Pure protein variant should not have its content stripped
        let (corrected, corrections) = strip_trailing_annotation("NP_000079.2:p.Val600Glu");
        assert_eq!(corrected, "NP_000079.2:p.Val600Glu");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_strip_trailing_annotation_all_parse_gaps() {
        // Test all 17 parse gap patterns from parse_gaps.json
        let patterns = [
            (
                "NM_003467.3(CXCR4):c.708G>A (p.Lys236=)",
                "NM_003467.3(CXCR4):c.708G>A",
            ),
            (
                "NM_001267550.2(TTN):c.30570C>A (p.Thr10190=)",
                "NM_001267550.2(TTN):c.30570C>A",
            ),
            (
                "NM_000199.5(SGSH):c.1428C>T (p.His476=)",
                "NM_000199.5(SGSH):c.1428C>T",
            ),
            (
                "NM_001366385.1(CARD14):c.27C>T (p.Ser9=)",
                "NM_001366385.1(CARD14):c.27C>T",
            ),
            (
                "NM_001267550.2(TTN):c.49476T>C (p.Pro16492=)",
                "NM_001267550.2(TTN):c.49476T>C",
            ),
            (
                "NM_198253.3(TERT):c.603G>A (p.Arg201=)",
                "NM_198253.3(TERT):c.603G>A",
            ),
            (
                "NM_178127.5(ANGPTL5):c.1097dup (p.Asn366fs)",
                "NM_178127.5(ANGPTL5):c.1097dup",
            ),
            (
                "NM_032303.5(HSDL2):c.894A>G (p.Lys298=)",
                "NM_032303.5(HSDL2):c.894A>G",
            ),
            (
                "NM_021831.6(AGBL5):c.2083G>A (p.Val695Ile)",
                "NM_021831.6(AGBL5):c.2083G>A",
            ),
            (
                "NM_201555.2(FHL2):c.507C>T (p.Ile169=)",
                "NM_201555.2(FHL2):c.507C>T",
            ),
            (
                "NM_025137.4(SPG11):c.1821A>G (p.Ser607=)",
                "NM_025137.4(SPG11):c.1821A>G",
            ),
            (
                "NM_025137.4(SPG11):c.6892A>G (p.Ile2298Val)",
                "NM_025137.4(SPG11):c.6892A>G",
            ),
            (
                "NM_001148.6(ANK2):c.8484T>C (p.Asp2828=)",
                "NM_001148.6(ANK2):c.8484T>C",
            ),
            (
                "NM_001148.6(ANK2):c.6492T>G (p.Leu2164=)",
                "NM_001148.6(ANK2):c.6492T>G",
            ),
            (
                "NM_001148.6(ANK2):c.231G>A (p.Val77=)",
                "NM_001148.6(ANK2):c.231G>A",
            ),
            (
                "NM_000059.3(BRCA2):c.3570G>C (p.Arg1190=)",
                "NM_000059.3(BRCA2):c.3570G>C",
            ),
            (
                "NM_030773.4(TUBB1):c.1045G>A (p.Val349Ile)",
                "NM_030773.4(TUBB1):c.1045G>A",
            ),
        ];

        for (input, expected) in patterns {
            let (corrected, corrections) = strip_trailing_annotation(input);
            assert_eq!(corrected, expected, "Failed for input: {}", input);
            assert_eq!(corrections.len(), 1, "No correction for: {}", input);
        }
    }

    // ==========================================================================
    // Missing Coordinate Prefix Tests
    // ==========================================================================

    #[test]
    fn test_correct_missing_coordinate_prefix_nc() {
        let (corrected, corrections) = correct_missing_coordinate_prefix("NC_000017.11:12345A>G");
        assert_eq!(corrected, "NC_000017.11:g.12345A>G");
        assert_eq!(corrections.len(), 1);
        assert_eq!(
            corrections[0].error_type,
            ErrorType::MissingCoordinatePrefix
        );
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_uncertain_range() {
        let (corrected, corrections) =
            correct_missing_coordinate_prefix("NC_000017.11:(?_31094927)_(31377677_?)del");
        assert_eq!(corrected, "NC_000017.11:g.(?_31094927)_(31377677_?)del");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_ng() {
        let (corrected, corrections) = correct_missing_coordinate_prefix("NG_007489.1:100_200del");
        assert_eq!(corrected, "NG_007489.1:g.100_200del");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_already_present() {
        let (corrected, corrections) = correct_missing_coordinate_prefix("NC_000017.11:g.12345A>G");
        assert_eq!(corrected, "NC_000017.11:g.12345A>G");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_nm() {
        // NM_ is transcript, not genomic, so should not auto-add g.
        let (corrected, corrections) = correct_missing_coordinate_prefix("NM_000088.3:459A>G");
        assert_eq!(corrected, "NM_000088.3:459A>G");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_lrg() {
        let (corrected, corrections) = correct_missing_coordinate_prefix("LRG_292:100_200del");
        assert_eq!(corrected, "LRG_292:g.100_200del");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_lrg_transcript() {
        // LRG transcript (LRG_*t*) should not auto-add g.
        let (corrected, corrections) = correct_missing_coordinate_prefix("LRG_292t1:100_200del");
        assert_eq!(corrected, "LRG_292t1:100_200del");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_missing_coordinate_prefix_bracket() {
        let (corrected, corrections) =
            correct_missing_coordinate_prefix("NC_000004.12:[144539078A>G]");
        assert_eq!(corrected, "NC_000004.12:g.[144539078A>G]");
        assert_eq!(corrections.len(), 1);
    }
}
