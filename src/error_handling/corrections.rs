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

/// Return true if the character is whitespace per Rust's `char::is_whitespace`,
/// or one of the common zero-width invisible characters that users frequently
/// paste from PDFs and rich-text sources:
/// - U+200B ZERO WIDTH SPACE
/// - U+200C ZERO WIDTH NON-JOINER
/// - U+200D ZERO WIDTH JOINER
/// - U+FEFF ZERO WIDTH NO-BREAK SPACE (BOM)
///
/// HGVS expressions never contain any of these; we treat them all as
/// `ExtraWhitespace` (W2003) and strip them under the same soft-warn contract.
#[inline]
fn is_invisible_whitespace(c: char) -> bool {
    c.is_whitespace() || matches!(c, '\u{200B}' | '\u{200C}' | '\u{200D}' | '\u{FEFF}')
}

/// Normalize extra whitespace in the input.
///
/// HGVS expressions are not permitted to contain whitespace (the spec is
/// explicit that the format `reference:description` has spaces "added for
/// clarity only"). This function therefore removes all whitespace from the
/// input — leading, trailing, and embedded — along with the zero-width
/// invisible characters U+200B/U+200C/U+200D/U+FEFF that are functionally
/// indistinguishable from whitespace. Each contiguous run of stripped
/// characters is recorded as a single [`DetectedCorrection`] of type
/// [`ErrorType::ExtraWhitespace`].
///
/// In strict mode the preprocessor will turn those corrections into a hard
/// rejection; in lenient/silent modes they become warnings (or are applied
/// silently). The function is idempotent on its own output. See issue #128.
///
/// Returns the corrected string and a list of corrections made.
pub fn correct_whitespace(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();
    let mut result = String::with_capacity(input.len());
    let mut run_start: Option<usize> = None;

    for (i, c) in input.char_indices() {
        if is_invisible_whitespace(c) {
            if run_start.is_none() {
                run_start = Some(i);
            }
        } else {
            if let Some(start) = run_start.take() {
                corrections.push(DetectedCorrection::new(
                    ErrorType::ExtraWhitespace,
                    input[start..i].to_string(),
                    String::new(),
                    start,
                    i,
                ));
            }
            result.push(c);
        }
    }

    // Trailing whitespace run, if any.
    if let Some(start) = run_start {
        corrections.push(DetectedCorrection::new(
            ErrorType::ExtraWhitespace,
            input[start..].to_string(),
            String::new(),
            start,
            input.len(),
        ));
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

/// Correct deprecated stop-codon and frameshift forms in protein descriptions.
///
/// Detects four deprecated forms within `p.` segments of the input and rewrites
/// them to their canonical `Ter`-based equivalents:
///
/// - `fsXN` (digits) → `fsTerN` — emits `DeprecatedFrameshiftX` (W3010).
/// - `fs*N` (digits) → `fsTerN` — emits `DeprecatedFrameshiftStar` (W3009).
/// - position-then-`X` at edit boundary → `Ter` — emits `DeprecatedStopCodonX` (W3008).
///   Distinguished from `Xaa` by requiring `X` not followed by a letter.
/// - position-then-`*` at edit boundary → `Ter` — emits `DeprecatedStopCodonStar` (W3007).
///   Distinguished from `fs*N` (already handled), from `c.*N` UTR positions
///   (digits do NOT precede `*` there), and from `delext*N` (preceding `t` is `t`).
///
/// The function preserves byte offsets relative to the *original* input by tracking
/// the input cursor, even when replacements change the output length. Only operates
/// when `p.` appears in the input. Idempotent: a second pass over canonical input
/// yields no corrections.
pub fn correct_deprecated_protein_forms(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();

    // Gate: only operate on inputs that contain a protein variant.
    if !input.starts_with("p.") && !input.contains(":p.") && !input.contains("[p.") {
        return (input.to_string(), corrections);
    }

    let bytes = input.as_bytes();
    let mut out = String::with_capacity(input.len());
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];

        // Detect "fs*N" or "fsXN" (frameshift termination notation).
        // Must begin at a literal "fs" boundary; the byte before the `f` is not
        // a lowercase letter (so we don't catch e.g. "ffs" — though that's
        // implausible in HGVS).
        if c == b'f'
            && i + 1 < bytes.len()
            && bytes[i + 1] == b's'
            && i + 2 < bytes.len()
            && (bytes[i + 2] == b'*' || bytes[i + 2] == b'X')
            && i + 3 < bytes.len()
            && bytes[i + 3].is_ascii_digit()
            && in_protein_segment(input, i + 2)
        {
            // Walk the digit run.
            let star_or_x = bytes[i + 2];
            let digits_start = i + 3;
            let mut digits_end = digits_start;
            while digits_end < bytes.len() && bytes[digits_end].is_ascii_digit() {
                digits_end += 1;
            }
            let digits = &input[digits_start..digits_end];
            let original = &input[i + 2..digits_end]; // "*23" or "X23"
            let corrected = format!("Ter{}", digits);
            let error_type = if star_or_x == b'*' {
                ErrorType::DeprecatedFrameshiftStar
            } else {
                ErrorType::DeprecatedFrameshiftX
            };
            corrections.push(DetectedCorrection::new(
                error_type,
                original.to_string(),
                corrected.clone(),
                i + 2,
                digits_end,
            ));
            // Emit "fs" + "Ter" + digits to output.
            out.push_str("fs");
            out.push_str(&corrected);
            i = digits_end;
            continue;
        }

        // Detect a stop-codon `*` or `X` immediately following a digit at the
        // edit boundary. The `*` or `X` must:
        // - be preceded by an ASCII digit (the position number),
        // - not be followed by a digit (else it's a position offset like c.*5),
        // - for `X`, not be followed by an ASCII letter (else it's `Xaa`),
        // - not be part of a `fs*N` / `fsXN` run (handled above).
        //
        // We also gate on protein context: only when, looking back, we find
        // `p.` before a top-level boundary (`:` or `[` or start) without an
        // intervening boundary that would put us back outside protein context.
        if (c == b'*' || c == b'X') && i > 0 && bytes[i - 1].is_ascii_digit() {
            let next_is_digit = i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit();
            let next_is_alpha = i + 1 < bytes.len() && bytes[i + 1].is_ascii_alphabetic();
            // For `X`, also reject when followed by an ASCII letter (Xaa, Xxx).
            let invalid_for_x = c == b'X' && next_is_alpha;
            // For `*`, `c.*5` UTR position has `*` BEFORE digits, not after, so
            // we don't need to filter that here; but `*N` followed by another
            // digit is rare in protein context — still, safer to skip.
            if !next_is_digit && !invalid_for_x && in_protein_segment(input, i) {
                let error_type = if c == b'*' {
                    ErrorType::DeprecatedStopCodonStar
                } else {
                    ErrorType::DeprecatedStopCodonX
                };
                corrections.push(DetectedCorrection::new(
                    error_type,
                    (c as char).to_string(),
                    "Ter".to_string(),
                    i,
                    i + 1,
                ));
                out.push_str("Ter");
                i += 1;
                continue;
            }
        }

        // Default: copy one UTF-8 char (1–4 bytes). The pattern branches above
        // only fire on ASCII bytes, so any non-ASCII byte at `i` is the start
        // of a multi-byte sequence — push the whole char to preserve UTF-8.
        if c.is_ascii() {
            out.push(c as char);
            i += 1;
        } else {
            let ch = input[i..]
                .chars()
                .next()
                .expect("input is valid UTF-8 and i is at a char boundary");
            let len = ch.len_utf8();
            out.push_str(&input[i..i + len]);
            i += len;
        }
    }

    (out, corrections)
}

/// Returns true if byte index `pos` in `input` is inside a `p.` (protein) segment.
///
/// Walks backward from `pos` looking for either a top-level boundary that
/// introduces a coordinate-type prefix (`:` after an accession; or start of
/// input), or for an embedded `p.` token. Brackets/parens/semicolons are
/// transparent — they are *inside* the protein segment, not boundaries to the
/// surrounding accession.
///
/// In practice we just look for the nearest preceding `p.` token before any
/// `:` (which would mark the end of the protein description on the left),
/// and confirm the input contains `p.` at all.
fn in_protein_segment(input: &str, pos: usize) -> bool {
    let bytes = input.as_bytes();
    // Walk back from pos, looking for "p." — but stop if we cross a ':'
    // (which would put us BEFORE the protein description, in the accession).
    let mut j = pos;
    while j >= 2 {
        if bytes[j - 2] == b'p' && bytes[j - 1] == b'.' {
            return true;
        }
        if bytes[j - 1] == b':' {
            // Crossed a colon — `pos` is on the left side of `:` (accession),
            // not inside a protein description.
            return false;
        }
        j -= 1;
    }
    // No "p." found before pos and no colon crossed; only true if input itself
    // begins with "p.".
    input.starts_with("p.")
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
/// Used by `correct_amino_acid_case_in_protein` to classify individual
/// three-letter tokens.
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

/// Scan only the `p.` segment of an HGVS expression and rewrite three-letter
/// amino-acid tokens to canonical capitalization (e.g. `val` -> `Val`,
/// `VAL` -> `Val`).
///
/// Records one `DetectedCorrection` per token rewritten. Inputs without a
/// `:p.` (or leading `p.`) variant-type marker are returned unchanged. This
/// scoping prevents false matches against accession bodies, gene-symbol
/// selectors, or coding-variant nucleotide tokens.
pub fn correct_amino_acid_case_in_protein(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();
    let Some(start) = find_protein_segment_start(input) else {
        return (input.to_string(), corrections);
    };
    let prefix = &input[..start];
    let body = &input[start..];

    let mut out = String::with_capacity(input.len());
    out.push_str(prefix);

    // Operate on byte offsets within `body`. AA tokens are pure ASCII, so
    // 3-byte windows align with 3-character windows whenever the leading
    // byte is ASCII. We treat non-ASCII bytes as opaque single chars and
    // copy them through verbatim using `char_indices`.
    let bytes = body.as_bytes();
    let len = bytes.len();
    let mut i = 0usize;
    while i < len {
        let c = bytes[i];
        if c.is_ascii_alphabetic() {
            // Greedy: collect the run of ASCII letters, then probe a 3-letter
            // window from the current offset against the AA table.
            let mut run_end = i;
            while run_end < len && bytes[run_end].is_ascii_alphabetic() {
                run_end += 1;
            }
            let mut j = i;
            while j + 3 <= run_end {
                let token = &body[j..j + 3];
                if let Some((canonical, error_type)) = correct_amino_acid_case(token) {
                    let abs_start = start + j;
                    corrections.push(DetectedCorrection::new(
                        error_type,
                        token.to_string(),
                        canonical.clone(),
                        abs_start,
                        abs_start + 3,
                    ));
                    out.push_str(&canonical);
                    j += 3;
                } else {
                    // Push one ASCII letter and advance.
                    out.push(bytes[j] as char);
                    j += 1;
                }
            }
            while j < run_end {
                out.push(bytes[j] as char);
                j += 1;
            }
            i = run_end;
        } else if c.is_ascii() {
            out.push(c as char);
            i += 1;
        } else {
            // Non-ASCII char: copy the full UTF-8 sequence verbatim.
            let ch_len = utf8_char_len(c);
            out.push_str(&body[i..i + ch_len]);
            i += ch_len;
        }
    }

    (out, corrections)
}

/// UTF-8 leading-byte → byte count.
///
/// Returns 1 for ASCII bytes (`< 0x80`) and for invalid continuation/lead
/// bytes (`< 0xC0`); consuming one byte at a time avoids stalling on
/// malformed input. Returns 2/3/4 for valid 2/3/4-byte UTF-8 sequences.
fn utf8_char_len(b: u8) -> usize {
    if b < 0xC0 {
        1
    } else if b < 0xE0 {
        2
    } else if b < 0xF0 {
        3
    } else {
        4
    }
}

/// Scan only the `p.` segment and expand uppercase one-letter amino-acid
/// codes (e.g. `V` -> `Val`) to their three-letter canonical forms.
///
/// Only **uppercase** ASCII letters are treated as one-letter AA candidates:
/// HGVS one-letter AA codes are uppercase, and lowercase letters inside the
/// `p.` segment belong to edit-type keywords (`del`, `ins`, `dup`, `inv`,
/// `con`, `delins`, `fs`, `ext`). The leading character of a canonical
/// three-letter token (`Ala`, `Arg`, ..., `Val`, `Ter`, `Xaa`) is also
/// skipped — by the time this phase runs (after W1001), every recognised
/// three-letter run already has canonical capitalisation.
///
/// One `DetectedCorrection` is emitted per replaced one-letter code.
pub fn correct_single_letter_aa_in_protein(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut corrections = Vec::new();
    let Some(start) = find_protein_segment_start(input) else {
        return (input.to_string(), corrections);
    };
    let prefix = &input[..start];
    let body = &input[start..];

    let mut out = String::with_capacity(input.len());
    out.push_str(prefix);

    let bytes = body.as_bytes();
    let len = bytes.len();
    let mut i = 0usize;
    while i < len {
        let c = bytes[i];

        // Skip three-letter canonical AA tokens already in place.
        if c.is_ascii_uppercase() && i + 3 <= len {
            let token = &body[i..i + 3];
            if is_canonical_three_letter_aa(token) {
                out.push_str(token);
                i += 3;
                continue;
            }
        }

        // Only uppercase ASCII letters are one-letter AA candidates.
        if c.is_ascii_uppercase() {
            if let Some(three) = single_to_three_letter_aa(c as char) {
                let abs_start = start + i;
                corrections.push(DetectedCorrection::new(
                    ErrorType::SingleLetterAminoAcid,
                    (c as char).to_string(),
                    three.to_string(),
                    abs_start,
                    abs_start + 1,
                ));
                out.push_str(three);
                i += 1;
                continue;
            }
        }

        if c.is_ascii() {
            out.push(c as char);
            i += 1;
        } else {
            let ch_len = utf8_char_len(c);
            out.push_str(&body[i..i + ch_len]);
            i += ch_len;
        }
    }

    (out, corrections)
}

/// Find the byte offset *within `input`* at which the protein description
/// body begins (one byte past the `p.` marker), or `None` if the input does
/// not have a protein variant type.
fn find_protein_segment_start(input: &str) -> Option<usize> {
    if let Some(idx) = input.find(":p.") {
        return Some(idx + 3);
    }
    if input.starts_with("p.") {
        return Some(2);
    }
    None
}

/// Returns true if `token` is a canonical three-letter amino-acid code.
fn is_canonical_three_letter_aa(token: &str) -> bool {
    matches!(
        token,
        "Ala"
            | "Arg"
            | "Asn"
            | "Asp"
            | "Cys"
            | "Gln"
            | "Glu"
            | "Gly"
            | "His"
            | "Ile"
            | "Leu"
            | "Lys"
            | "Met"
            | "Phe"
            | "Pro"
            | "Sec"
            | "Ser"
            | "Thr"
            | "Trp"
            | "Tyr"
            | "Val"
            | "Ter"
            | "Xaa"
    )
}

/// Convert single-letter amino acid code to three-letter code.
///
/// Used by `correct_single_letter_aa_in_protein` to expand W1002 candidates.
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

/// Find *all* RefSeq-style accessions in the input that lack a `.<version>`
/// suffix. Returns one `DetectedCorrection` per occurrence, with `original`
/// set to the unversioned accession and `corrected` left empty (W3001 is not
/// auto-correctable: we have no way to know which version was intended).
///
/// Walks the input character-by-character, identifying contiguous
/// `[A-Za-z_]+\d+` accession bodies, then classifying each by prefix.
/// Accessions whose prefix is in the optional-version set (Ensembl
/// ENST/ENSP/ENSG/ENSR) are skipped.
pub fn detect_missing_versions(input: &str) -> Vec<DetectedCorrection> {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if !bytes[i].is_ascii_uppercase() {
            i += 1;
            continue;
        }
        let prefix_start = i;
        let mut j = i;
        while j < bytes.len() && (bytes[j].is_ascii_alphabetic() || bytes[j] == b'_') {
            j += 1;
        }
        let alpha_end = j;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        let digit_end = j;

        if alpha_end == prefix_start || digit_end == alpha_end {
            i = (alpha_end + 1).max(i + 1);
            continue;
        }

        let accession_alpha = &input[prefix_start..alpha_end];
        let accession_full = &input[prefix_start..digit_end];

        let (recognised, optional) = classify_accession_prefix(accession_alpha);
        if !recognised {
            i = digit_end;
            continue;
        }

        let has_version = digit_end < bytes.len()
            && bytes[digit_end] == b'.'
            && bytes
                .get(digit_end + 1)
                .copied()
                .map(|b| b.is_ascii_digit())
                .unwrap_or(false);

        if !has_version && !optional {
            // LRG transcript (`LRG_NNNtM`) and protein (`LRG_NNNpM`) refs
            // encode the version as an alphabetic suffix immediately after
            // the numeric ID, not as `.<version>`. Skip them rather than
            // flagging a false W3001.
            let is_lrg_suffix_ref = accession_alpha == "LRG_"
                && digit_end < bytes.len()
                && bytes[digit_end].is_ascii_alphabetic();
            if !is_lrg_suffix_ref {
                hits.push(DetectedCorrection::new(
                    ErrorType::MissingVersion,
                    accession_full.to_string(),
                    String::new(),
                    prefix_start,
                    digit_end,
                ));
            }
        }

        i = digit_end;
    }
    hits
}

/// Classify an accession's leading alphabetic prefix.
///
/// Returns `(recognised, optional_version)`:
/// - `recognised = true` when the prefix matches a known RefSeq/Ensembl/LRG
///   accession family.
/// - `optional_version = true` when that family routinely appears in the
///   wild without a version suffix (Ensembl).
fn classify_accession_prefix(prefix: &str) -> (bool, bool) {
    match prefix {
        "NM_" | "NP_" | "NC_" | "NG_" | "NR_" | "NT_" | "NW_" | "XM_" | "XP_" | "XR_" | "LRG_" => {
            (true, false)
        }
        "ENST" | "ENSP" | "ENSG" | "ENSR" => (true, true),
        _ => (false, false),
    }
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

// =============================================================================
// Issue #127 — non-canonical input forms (SVA-007 / 008 / 009 / 010 / 014 / 027)
// =============================================================================

/// Returns true if `input` contains a non-protein coordinate-type prefix
/// (`g.`, `c.`, `n.`, `r.`, `m.`, `o.`).
///
/// Used to scope DNA/RNA-only checks and exclude protein descriptions, where
/// the same surface syntax (`del32`, `_<n>_<n>del`) has different semantics
/// and a different spec section.
fn has_non_protein_description(bytes: &[u8]) -> bool {
    let mut i = 0usize;
    while i + 1 < bytes.len() {
        let coord = bytes[i];
        let dot = bytes[i + 1];
        if dot == b'.'
            && (coord == b'g'
                || coord == b'c'
                || coord == b'n'
                || coord == b'r'
                || coord == b'm'
                || coord == b'o')
        {
            // Require the byte before `coord` to be a colon, an open paren /
            // bracket, a semicolon, or the start of input. This avoids
            // matching `del.` / interior letters within a token.
            let prev_ok = i == 0 || matches!(bytes[i - 1], b':' | b'(' | b'[' | b';');
            if prev_ok {
                return true;
            }
        }
        i += 1;
    }
    false
}

/// Returns the byte index immediately after the UTF-8 character starting at
/// `i` in `bytes`. Falls back to `i + 1` for malformed input.
#[inline]
fn next_char_end(bytes: &[u8], i: usize) -> usize {
    if i >= bytes.len() {
        return i;
    }
    let b = bytes[i];
    let extra = if b < 0x80 {
        0
    } else if b & 0xE0 == 0xC0 {
        1
    } else if b & 0xF0 == 0xE0 {
        2
    } else if b & 0xF8 == 0xF0 {
        3
    } else {
        0
    };
    (i + 1 + extra).min(bytes.len())
}

/// Detect deletions described with a size-count suffix instead of a position
/// range (W3011, SVA-007).
///
/// Matches `del<digits>` occurrences inside non-protein descriptions (after
/// `g.`, `c.`, `n.`, `r.`, `m.`, or `o.`), excluding the `delins` keyword.
/// Returns one `DetectedCorrection` per occurrence with `corrected` empty:
/// W3011 is not auto-correctable because synthesizing the end position
/// requires intronic / UTR-aware coordinate arithmetic that the preprocessor
/// does not have.
pub fn detect_del_size_suffix(input: &str) -> Vec<DetectedCorrection> {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();
    if !has_non_protein_description(bytes) {
        return hits;
    }

    let mut i = 0usize;
    while i + 3 <= bytes.len() {
        if &bytes[i..i + 3] != b"del" {
            i += 1;
            continue;
        }
        // Skip `delins`.
        if i + 6 <= bytes.len() && &bytes[i + 3..i + 6] == b"ins" {
            i += 6;
            continue;
        }
        // Look for one-or-more digits after `del`.
        let digit_start = i + 3;
        let mut j = digit_start;
        while j < bytes.len() && bytes[j].is_ascii_digit() {
            j += 1;
        }
        if j == digit_start {
            i += 3;
            continue;
        }
        // The token must end here — followed by end-of-string or one of the
        // HGVS terminator bytes (`)`, `;`, `]`, whitespace).
        let end_byte = bytes.get(j).copied();
        let is_terminator =
            end_byte.is_none_or(|b| b == b')' || b == b';' || b == b']' || b.is_ascii_whitespace());
        if !is_terminator {
            i = j;
            continue;
        }
        hits.push(DetectedCorrection::new(
            ErrorType::DelSizeSuffix,
            &input[i..j],
            String::new(),
            i,
            j,
        ));
        i = j;
    }
    hits
}

/// Detect and correct deletion-insertions whose inserted sequence is empty
/// (W3012, SVA-010).
///
/// Matches `delins` followed by a terminator (end-of-input, `)`, `;`, `]`,
/// or whitespace). Rewrites to `del` and emits one `DetectedCorrection`
/// per occurrence. Pre-existing `delinsATG` / `delins[…]` / `delinsN[12]` /
/// `delins(10_20)` / `delins10` forms are left untouched.
pub fn correct_empty_delins(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();
    if !has_non_protein_description(bytes) {
        return (input.to_string(), hits);
    }

    let mut result = String::with_capacity(input.len());
    let mut i = 0usize;
    while i < bytes.len() {
        if i + 6 <= bytes.len() && &bytes[i..i + 6] == b"delins" {
            let after = bytes.get(i + 6).copied();
            let is_empty = match after {
                None => true,
                Some(b) => b == b')' || b == b';' || b == b']' || b.is_ascii_whitespace(),
            };
            if is_empty {
                hits.push(DetectedCorrection::new(
                    ErrorType::EmptyDelinsInsert,
                    "delins",
                    "del",
                    i,
                    i + 6,
                ));
                result.push_str("del");
                i += 6;
                continue;
            }
        }
        let ch_end = next_char_end(bytes, i);
        result.push_str(&input[i..ch_end]);
        i = ch_end;
    }
    (result, hits)
}

/// If `bytes` starting at `i` matches `<sign?><digits>_<sign?><digits>` and
/// the two numeric tokens (with sign) are equal, return `(pair_end,
/// canonical_single_position_text)`; otherwise `None`.
fn match_equal_position_pair(bytes: &[u8], i: usize, input: &str) -> Option<(usize, String)> {
    let mut j = i;
    let first_start = j;
    if j < bytes.len() && bytes[j] == b'-' {
        j += 1;
    }
    let num1_digit_start = j;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    if j == num1_digit_start {
        return None;
    }
    let first_end = j;
    if j >= bytes.len() || bytes[j] != b'_' {
        return None;
    }
    j += 1;
    let second_start = j;
    if j < bytes.len() && bytes[j] == b'-' {
        j += 1;
    }
    let num2_digit_start = j;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    if j == num2_digit_start {
        return None;
    }
    let second_end = j;

    let first = &input[first_start..first_end];
    let second = &input[second_start..second_end];
    if first != second {
        return None;
    }
    Some((j, first.to_string()))
}

/// Returns true when `bytes` at `i` starts with one of `del` (not `delins`),
/// `dup`, or `inv`.
fn matches_single_pos_keyword(bytes: &[u8], i: usize) -> bool {
    if i + 3 > bytes.len() {
        return false;
    }
    let kw = &bytes[i..i + 3];
    if kw == b"del" {
        // Reject `delins`.
        return !(i + 6 <= bytes.len() && &bytes[i + 3..i + 6] == b"ins");
    }
    kw == b"dup" || kw == b"inv"
}

/// Detect and collapse single-position ranges in `del` / `dup` / `inv`
/// descriptions (W4003, SVA-008 / 009 / 014).
///
/// Matches `<sign?><digits>_<sign?><digits>` where the two integers are
/// equal, followed by `del` (but not `delins`), `dup`, or `inv`. Rewrites to
/// the single-position form and emits one `DetectedCorrection` per
/// occurrence. The check is scoped to non-protein descriptions; protein has
/// its own L2 issue track.
pub fn correct_single_position_range(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();
    if !has_non_protein_description(bytes) {
        return (input.to_string(), hits);
    }

    let mut result = String::with_capacity(input.len());
    let mut p = 0usize;
    while p < bytes.len() {
        if let Some((pair_end, value_text)) = match_equal_position_pair(bytes, p, input) {
            if matches_single_pos_keyword(bytes, pair_end) {
                hits.push(DetectedCorrection::new(
                    ErrorType::SinglePositionRange,
                    &input[p..pair_end],
                    value_text.clone(),
                    p,
                    pair_end,
                ));
                result.push_str(&value_text);
                p = pair_end;
                continue;
            }
        }
        let ch_end = next_char_end(bytes, p);
        result.push_str(&input[p..ch_end]);
        p = ch_end;
    }
    (result, hits)
}

/// If `bytes` starting at `i` matches `<sign?><digits>_<sign?><digits>`,
/// return the byte index just past the second number; otherwise `None`. Does
/// NOT require the two numbers to be equal — used to find candidate
/// positions-pair shapes for the redundant-repeat-label detector.
fn match_position_pair(bytes: &[u8], i: usize) -> Option<usize> {
    let mut j = i;
    if j < bytes.len() && bytes[j] == b'-' {
        j += 1;
    }
    let num1_start = j;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    if j == num1_start {
        return None;
    }
    if j >= bytes.len() || bytes[j] != b'_' {
        return None;
    }
    j += 1;
    if j < bytes.len() && bytes[j] == b'-' {
        j += 1;
    }
    let num2_start = j;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    if j == num2_start {
        return None;
    }
    Some(j)
}

/// Detect and strip redundant base labels in RNA repeat descriptions (W3013,
/// SVA-027).
///
/// Matches `r.<num1>_<num2><rna-bases>[<digits>...]` where `<rna-bases>` is
/// a run of one or more lowercase a/c/g/u characters between the second
/// position and the `[` count. Rewrites to drop the base label.
///
/// Scoped to `r.` (RNA) descriptions: the HGVS
/// `recommendations/RNA/repeated.md` page is the only one that calls this
/// form non-canonical; DNA repeats keep the base label by convention.
pub fn correct_redundant_repeat_label(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();

    // Locate the `r.` description start.
    let mut desc_start = None;
    let mut i = 0usize;
    while i + 1 < bytes.len() {
        if bytes[i] == b'r'
            && bytes[i + 1] == b'.'
            && (i == 0 || matches!(bytes[i - 1], b':' | b'(' | b'[' | b';'))
        {
            desc_start = Some(i + 2);
            break;
        }
        i += 1;
    }
    let Some(start) = desc_start else {
        return (input.to_string(), hits);
    };

    let mut result = String::with_capacity(input.len());
    result.push_str(&input[..start]);
    let mut p = start;
    while p < bytes.len() {
        // Stop scanning once we cross into a non-RNA description so we don't
        // strip lowercase repeats that happen to live inside a later
        // `c.`/`g.`/`m.`/`n.`/`o.`/`p.` description (e.g.
        // `r.100_102cug[4];c.50_52acg[3]`).
        if matches!(bytes[p], b'c' | b'g' | b'm' | b'n' | b'o' | b'p')
            && bytes.get(p + 1).copied() == Some(b'.')
            && p > start
            && matches!(bytes[p - 1], b':' | b'(' | b'[' | b';' | b'|')
        {
            break;
        }

        let pair_start = p;
        let pair = match_position_pair(bytes, p);
        if let Some(after_pair) = pair {
            let bases_start = after_pair;
            let mut bases_end = bases_start;
            while bases_end < bytes.len() && matches!(bytes[bases_end], b'a' | b'c' | b'g' | b'u') {
                bases_end += 1;
            }
            if bases_end > bases_start && bytes.get(bases_end).copied() == Some(b'[') {
                result.push_str(&input[pair_start..bases_start]);
                hits.push(DetectedCorrection::new(
                    ErrorType::RedundantRepeatLabel,
                    &input[bases_start..bases_end],
                    String::new(),
                    bases_start,
                    bases_end,
                ));
                p = bases_end;
                continue;
            }
        }
        let ch_end = next_char_end(bytes, p);
        result.push_str(&input[p..ch_end]);
        p = ch_end;
    }
    result.push_str(&input[p..]);
    (result, hits)
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

    #[test]
    fn test_correct_whitespace_strips_zero_width_chars() {
        // U+200B ZERO WIDTH SPACE, U+FEFF ZWNBSP/BOM, U+200C ZWNJ, U+200D ZWJ
        // are functionally invisible; users frequently paste them from PDFs/web.
        // Treat them like whitespace under W2003.
        let input = "NM_000088.3:c.100\u{200B}A>G";
        let (corrected, corrections) = correct_whitespace(input);
        assert_eq!(corrected, "NM_000088.3:c.100A>G");
        assert_eq!(corrections.len(), 1, "should record one correction");
        assert_eq!(corrections[0].error_type, ErrorType::ExtraWhitespace);

        let input = "\u{FEFF}NM_000088.3:c.100A>G";
        let (corrected, _) = correct_whitespace(input);
        assert_eq!(corrected, "NM_000088.3:c.100A>G");

        let input = "c.100\u{200C}A>G";
        let (corrected, _) = correct_whitespace(input);
        assert_eq!(corrected, "c.100A>G");

        let input = "c.100\u{200D}A>G";
        let (corrected, _) = correct_whitespace(input);
        assert_eq!(corrected, "c.100A>G");
    }

    #[test]
    fn test_correct_whitespace_embedded_single_run_one_warning() {
        // Multiple whitespace chars in one contiguous run = ONE correction.
        let (corrected, corrections) = correct_whitespace("c.100   A>G");
        assert_eq!(corrected, "c.100A>G");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].original, "   ");
        assert_eq!(corrections[0].start, 5);
        assert_eq!(corrections[0].end, 8);
    }

    #[test]
    fn test_correct_whitespace_multiple_runs_count() {
        // Three separate whitespace runs = THREE corrections.
        let (corrected, corrections) = correct_whitespace(" c.100 A> G");
        assert_eq!(corrected, "c.100A>G");
        assert_eq!(corrections.len(), 3);
    }

    #[test]
    fn test_correct_whitespace_idempotent() {
        // Re-normalizing a corrected output produces zero corrections.
        let (first, _) = correct_whitespace("  c. 100 A>G  ");
        let (second, corrections) = correct_whitespace(&first);
        assert_eq!(first, second);
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_whitespace_mixed_unicode_whitespace() {
        // Tab, NBSP, em-space, vertical tab — all stripped, one run = one warning.
        let input = "c.100\t\u{00A0}\u{2003}\u{000B}A>G";
        let (corrected, corrections) = correct_whitespace(input);
        assert_eq!(corrected, "c.100A>G");
        assert_eq!(corrections.len(), 1);
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

    // Protein-scoped AA case correction (W1001) tests
    #[test]
    fn test_correct_amino_acid_case_in_protein_lowercase() {
        let (corrected, corrections) =
            correct_amino_acid_case_in_protein("NP_000079.2:p.val600glu");
        assert_eq!(corrected, "NP_000079.2:p.Val600Glu");
        assert_eq!(corrections.len(), 2);
        assert!(corrections
            .iter()
            .all(|c| c.error_type == ErrorType::LowercaseAminoAcid));
        assert_eq!(corrections[0].original, "val");
        assert_eq!(corrections[0].corrected, "Val");
        assert_eq!(corrections[1].original, "glu");
        assert_eq!(corrections[1].corrected, "Glu");
    }

    #[test]
    fn test_correct_amino_acid_case_in_protein_uppercase() {
        let (corrected, corrections) =
            correct_amino_acid_case_in_protein("NP_000079.2:p.VAL600GLU");
        assert_eq!(corrected, "NP_000079.2:p.Val600Glu");
        assert_eq!(corrections.len(), 2);
    }

    #[test]
    fn test_correct_amino_acid_case_in_protein_canonical_no_warning() {
        let (corrected, corrections) =
            correct_amino_acid_case_in_protein("NP_000079.2:p.Val600Glu");
        assert_eq!(corrected, "NP_000079.2:p.Val600Glu");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_amino_acid_case_in_protein_only_runs_on_protein() {
        let (corrected, corrections) = correct_amino_acid_case_in_protein("NM_000088.3:c.100A>G");
        assert_eq!(corrected, "NM_000088.3:c.100A>G");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_amino_acid_case_in_protein_predicted_parens() {
        let (corrected, corrections) =
            correct_amino_acid_case_in_protein("NP_000079.2:p.(val600glu)");
        assert_eq!(corrected, "NP_000079.2:p.(Val600Glu)");
        assert_eq!(corrections.len(), 2);
    }

    #[test]
    fn test_correct_amino_acid_case_in_protein_idempotent() {
        let (once, _) = correct_amino_acid_case_in_protein("NP_000079.2:p.val600glu");
        let (twice, corrections) = correct_amino_acid_case_in_protein(&once);
        assert_eq!(twice, once);
        assert!(corrections.is_empty());
    }

    // Protein-scoped single-letter AA expansion (W1002) tests
    #[test]
    fn test_correct_single_letter_aa_substitution() {
        let (corrected, corrections) = correct_single_letter_aa_in_protein("NP_000079.2:p.V600E");
        assert_eq!(corrected, "NP_000079.2:p.Val600Glu");
        assert_eq!(corrections.len(), 2);
        assert!(corrections
            .iter()
            .all(|c| c.error_type == ErrorType::SingleLetterAminoAcid));
    }

    #[test]
    fn test_correct_single_letter_aa_canonical_no_warning() {
        let (corrected, corrections) =
            correct_single_letter_aa_in_protein("NP_000079.2:p.Val600Glu");
        assert_eq!(corrected, "NP_000079.2:p.Val600Glu");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_single_letter_aa_no_p_segment() {
        let (corrected, corrections) = correct_single_letter_aa_in_protein("NM_000088.3:c.459A>G");
        assert_eq!(corrected, "NM_000088.3:c.459A>G");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_single_letter_aa_predicted() {
        let (corrected, corrections) = correct_single_letter_aa_in_protein("NP_000079.2:p.(V600E)");
        assert_eq!(corrected, "NP_000079.2:p.(Val600Glu)");
        assert_eq!(corrections.len(), 2);
    }

    #[test]
    fn test_correct_single_letter_aa_delins() {
        let (corrected, corrections) =
            correct_single_letter_aa_in_protein("NP_000079.2:p.V600_E601delinsK");
        assert_eq!(corrected, "NP_000079.2:p.Val600_Glu601delinsLys");
        assert_eq!(corrections.len(), 3);
    }

    #[test]
    fn test_correct_single_letter_aa_idempotent() {
        let (once, _) = correct_single_letter_aa_in_protein("NP_000079.2:p.V600E");
        let (twice, corrections) = correct_single_letter_aa_in_protein(&once);
        assert_eq!(twice, once);
        assert!(corrections.is_empty());
    }

    // Missing-version (W3001) detection tests
    #[test]
    fn test_detect_missing_versions_single() {
        let hits = detect_missing_versions("NM_000088:c.100A>G");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].original, "NM_000088");
        assert_eq!(hits[0].error_type, ErrorType::MissingVersion);
        assert_eq!(hits[0].start, 0);
        assert_eq!(hits[0].end, "NM_000088".len());
    }

    #[test]
    fn test_detect_missing_versions_with_version_no_hit() {
        let hits = detect_missing_versions("NM_000088.3:c.100A>G");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_detect_missing_versions_inner_accession() {
        let hits = detect_missing_versions("NG_012232(NM_004006.2):c.93+1G>T");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].original, "NG_012232");
    }

    #[test]
    fn test_detect_missing_versions_both_missing() {
        let hits = detect_missing_versions("NG_012232(NM_004006):c.93+1G>T");
        assert_eq!(hits.len(), 2);
        let originals: Vec<&str> = hits.iter().map(|h| h.original.as_str()).collect();
        assert!(originals.contains(&"NG_012232"));
        assert!(originals.contains(&"NM_004006"));
    }

    #[test]
    fn test_detect_missing_versions_ensembl_optional() {
        let hits = detect_missing_versions("ENST00000380152:c.100A>G");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_detect_missing_versions_no_accession() {
        let hits = detect_missing_versions("c.100A>G");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_detect_missing_versions_lrg() {
        let hits = detect_missing_versions("LRG_199:c.100A>G");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].original, "LRG_199");
    }

    #[test]
    fn test_detect_missing_versions_lrg_transcript_no_hit() {
        // LRG transcript references (`LRG_NNNtM`) carry the version as a
        // `t<digit>` suffix directly after the numeric ID, not as `.<digit>`.
        let hits = detect_missing_versions("LRG_292t1:c.100A>G");
        assert!(
            hits.is_empty(),
            "LRG transcript suffix should not trigger W3001, got {hits:?}",
        );
    }

    #[test]
    fn test_detect_missing_versions_lrg_protein_no_hit() {
        // Same for LRG protein references (`LRG_NNNpM`).
        let hits = detect_missing_versions("LRG_292p1:p.Ser68Arg");
        assert!(
            hits.is_empty(),
            "LRG protein suffix should not trigger W3001, got {hits:?}",
        );
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

    // ----------------------------------------------------------------------
    // correct_deprecated_protein_forms — SVA-003..SVA-006 (issue #125)
    // ----------------------------------------------------------------------

    #[test]
    fn test_deprecated_stop_star_substitution() {
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.Arg97*");
        assert_eq!(corrected, "NP_000079.2:p.Arg97Ter");
        assert_eq!(corrections.len(), 1);
        assert_eq!(
            corrections[0].error_type,
            ErrorType::DeprecatedStopCodonStar
        );
        assert_eq!(corrections[0].original, "*");
        assert_eq!(corrections[0].corrected, "Ter");
    }

    #[test]
    fn test_deprecated_stop_x_substitution() {
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.Arg97X");
        assert_eq!(corrected, "NP_000079.2:p.Arg97Ter");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::DeprecatedStopCodonX);
    }

    #[test]
    fn test_deprecated_frameshift_star() {
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.Arg97fs*23");
        assert_eq!(corrected, "NP_000079.2:p.Arg97fsTer23");
        assert_eq!(corrections.len(), 1);
        assert_eq!(
            corrections[0].error_type,
            ErrorType::DeprecatedFrameshiftStar
        );
        assert_eq!(corrections[0].original, "*23");
        assert_eq!(corrections[0].corrected, "Ter23");
    }

    #[test]
    fn test_deprecated_frameshift_x() {
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.Arg97fsX23");
        assert_eq!(corrected, "NP_000079.2:p.Arg97fsTer23");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::DeprecatedFrameshiftX);
    }

    #[test]
    fn test_deprecated_frameshift_with_new_aa() {
        // p.Arg97ProfsTer23 with deprecated Star or X.
        let (corrected, corrections) =
            correct_deprecated_protein_forms("NP_000079.2:p.Arg97Profs*23");
        assert_eq!(corrected, "NP_000079.2:p.Arg97ProfsTer23");
        assert_eq!(corrections.len(), 1);
        assert_eq!(
            corrections[0].error_type,
            ErrorType::DeprecatedFrameshiftStar
        );

        let (corrected, corrections) =
            correct_deprecated_protein_forms("NP_000079.2:p.Arg97ProfsX23");
        assert_eq!(corrected, "NP_000079.2:p.Arg97ProfsTer23");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::DeprecatedFrameshiftX);
    }

    #[test]
    fn test_deprecated_canonical_input_no_corrections() {
        // Canonical inputs must not trigger any warning.
        for input in [
            "NP_000079.2:p.Arg97Ter",
            "NP_000079.2:p.Arg97ProfsTer23",
            "NP_000079.2:p.Tyr180fs",
            "NP_000079.2:p.Val600Glu",
            "NP_000079.2:p.Arg782Xaa",
            "NP_000079.2:p.Met1ext-5",
            "NP_001166937.1:p.Ter514LeuextTer?",
        ] {
            let (corrected, corrections) = correct_deprecated_protein_forms(input);
            assert_eq!(corrected, input, "input changed: {}", input);
            assert!(
                corrections.is_empty(),
                "expected no corrections for {}, got {:?}",
                input,
                corrections
            );
        }
    }

    // ========================================================================
    // Issue #127 — non-canonical input forms (W3011 / W3012 / W3013 / W4003)
    // ========================================================================

    // --- W3011 DelSizeSuffix ---

    #[test]
    fn test_detect_del_size_suffix_basic() {
        let hits = detect_del_size_suffix("NG_012232.1:g.123del6");
        assert_eq!(hits.len(), 1, "expected one hit, got {:?}", hits);
        let h = &hits[0];
        assert_eq!(h.error_type, ErrorType::DelSizeSuffix);
        assert_eq!(h.original, "del6");
    }

    #[test]
    fn test_detect_del_size_suffix_canonical_no_hit() {
        let hits = detect_del_size_suffix("NG_012232.1:g.123_128del");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_detect_del_size_suffix_plain_del_no_hit() {
        let hits = detect_del_size_suffix("NM_000088.3:c.123del");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_detect_del_size_suffix_delins_no_hit() {
        // `delinsNN` and `delins<digit><N>` patterns must not match — `delins`
        // is a different keyword. Empty `delins` is handled separately.
        let hits = detect_del_size_suffix("NC_000001.11:g.100_102delins10");
        assert!(hits.is_empty());
        let hits = detect_del_size_suffix("NC_000001.11:g.100_102delinsATG");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_detect_del_size_suffix_with_range_still_hits() {
        // `g.100_120del6` — range form with size suffix is still a violation
        // (the size is redundant when both endpoints are given).
        let hits = detect_del_size_suffix("NG_012232.1:g.100_120del6");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_detect_del_size_suffix_protein_no_hit() {
        // Protein `delN` is a deletion of N residues but the spec section
        // (DNA/deletion) does not apply; protein has its own L2 track.
        let hits = detect_del_size_suffix("NP_000079.2:p.Lys100del32");
        assert!(hits.is_empty());
    }

    // --- W3012 EmptyDelinsInsert ---

    #[test]
    fn test_correct_empty_delins_basic() {
        let (out, hits) = correct_empty_delins("NC_000001.11:g.100_102delins");
        assert_eq!(out, "NC_000001.11:g.100_102del");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].error_type, ErrorType::EmptyDelinsInsert);
        assert_eq!(hits[0].original, "delins");
        assert_eq!(hits[0].corrected, "del");
    }

    #[test]
    fn test_correct_empty_delins_with_payload_no_change() {
        for input in [
            "NC_000001.11:g.100_102delinsATG",
            "NC_000001.11:g.100_102delinsN[12]",
            "NC_000001.11:g.100_102delins[A;G]",
            "NC_000001.11:g.100_102delins(10_20)",
            "NC_000001.11:g.100_102delins10",
        ] {
            let (out, hits) = correct_empty_delins(input);
            assert_eq!(out, input, "input {} should not change", input);
            assert!(
                hits.is_empty(),
                "input {} should not warn, got {:?}",
                input,
                hits
            );
        }
    }

    #[test]
    fn test_deprecated_idempotent() {
        // Re-running over already-corrected text must be a no-op.
        let once = correct_deprecated_protein_forms("NP_000079.2:p.Arg97fs*23").0;
        let twice = correct_deprecated_protein_forms(&once);
        assert_eq!(twice.0, once);
        assert!(twice.1.is_empty());
    }

    #[test]
    fn test_deprecated_no_protein_context_no_corrections() {
        // CDS coords with a literal '*' (e.g. 3'UTR position c.*5A>G) must NOT
        // be touched — only protein contexts are eligible.
        let (corrected, corrections) = correct_deprecated_protein_forms("NM_000088.3:c.123*");
        assert_eq!(corrected, "NM_000088.3:c.123*");
        assert!(corrections.is_empty());

        let (corrected, corrections) = correct_deprecated_protein_forms("NM_000088.3:c.123X");
        assert_eq!(corrected, "NM_000088.3:c.123X");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_deprecated_multiple_in_compound_allele() {
        // Compound protein allele with two deprecated forms should emit two
        // corrections in order.
        let (corrected, corrections) =
            correct_deprecated_protein_forms("NP_000079.2:p.[Arg97*;Arg100X]");
        assert_eq!(corrected, "NP_000079.2:p.[Arg97Ter;Arg100Ter]");
        assert_eq!(corrections.len(), 2);
        assert_eq!(
            corrections[0].error_type,
            ErrorType::DeprecatedStopCodonStar
        );
        assert_eq!(corrections[1].error_type, ErrorType::DeprecatedStopCodonX);
    }

    #[test]
    fn test_deprecated_does_not_touch_xaa_or_extension() {
        // 'Xaa' (any amino acid) must NOT be flagged as deprecated stop X.
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.Arg782Xaa");
        assert_eq!(corrected, "NP_000079.2:p.Arg782Xaa");
        assert!(corrections.is_empty());

        // 'extTer?' / 'ext*?' / 'extTer17' / 'ext*17' are extension notations.
        // The '*?' / '*17' here is canonical-equivalent to 'Ter?'/'Ter17'; but
        // by this corrector's rule the '*' is preceded by 't' (from "ext"), not
        // a digit — so it is not modified.
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.Met1ext*-5");
        assert_eq!(corrected, "NP_000079.2:p.Met1ext*-5");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_deprecated_predicted_paren_form() {
        // Predicted protein consequences in parentheses still match.
        let (corrected, corrections) = correct_deprecated_protein_forms("NP_000079.2:p.(Arg97*)");
        assert_eq!(corrected, "NP_000079.2:p.(Arg97Ter)");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_deprecated_byte_offsets_track_original_input() {
        // The DetectedCorrection's start/end refer to the ORIGINAL input bytes,
        // not the rewritten output, so downstream span reporting is accurate.
        let input = "NP_000079.2:p.Arg97fs*23";
        let (_, corrections) = correct_deprecated_protein_forms(input);
        assert_eq!(corrections.len(), 1);
        // "*23" begins at byte index of the '*' in the input.
        let star_pos = input.find('*').unwrap();
        assert_eq!(corrections[0].start, star_pos);
        assert_eq!(corrections[0].end, input.len()); // up through "23"
    }

    #[test]
    fn test_deprecated_preserves_non_ascii_utf8() {
        // Non-ASCII bytes elsewhere in the input must pass through unchanged
        // (the byte-walking loop must not split multi-byte UTF-8 into stray
        // Latin-1 codepoints). Input is not valid HGVS but the corrector
        // should be UTF-8 correct regardless.
        let input = "NP_000079.2:p.Arg97* αβ";
        let (corrected, corrections) = correct_deprecated_protein_forms(input);
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrected, "NP_000079.2:p.Arg97Ter αβ");
        // Round-trip the rewritten output through valid UTF-8 char count.
        assert_eq!(corrected.chars().filter(|c| !c.is_ascii()).count(), 2);
    }

    #[test]
    fn test_deprecated_fs_branch_gated_to_protein_segment() {
        // 'fs*N' outside the protein segment (here, in the accession-side text
        // before the ':p.' boundary) must not be rewritten, even when a `p.`
        // segment exists elsewhere in the input. Mirrors the per-position
        // gate the stop-codon branch already enforces.
        let input = "fs*23:p.Arg97Ter";
        let (corrected, corrections) = correct_deprecated_protein_forms(input);
        assert_eq!(corrected, input);
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_empty_delins_in_compound_allele() {
        let (out, hits) = correct_empty_delins("NM_000088.3:c.[100_102delins;200T>G]");
        assert_eq!(out, "NM_000088.3:c.[100_102del;200T>G]");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_empty_delins_idempotent() {
        let (out1, hits1) = correct_empty_delins("NC_000001.11:g.100_102delins");
        let (out2, hits2) = correct_empty_delins(&out1);
        assert_eq!(out2, "NC_000001.11:g.100_102del");
        assert_eq!(hits1.len(), 1);
        assert!(hits2.is_empty());
    }

    // --- W3013 RedundantRepeatLabel ---

    #[test]
    fn test_correct_redundant_repeat_label_basic() {
        let (out, hits) = correct_redundant_repeat_label("NM_000088.3:r.100_102cug[4]");
        assert_eq!(out, "NM_000088.3:r.100_102[4]");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].error_type, ErrorType::RedundantRepeatLabel);
        assert_eq!(hits[0].original, "cug");
        assert_eq!(hits[0].corrected, "");
    }

    #[test]
    fn test_correct_redundant_repeat_label_negative_positions() {
        // The exemplar from the spec: `r.-125_-123cug[4]`
        let (out, hits) = correct_redundant_repeat_label("NM_000088.3:r.-125_-123cug[4]");
        assert_eq!(out, "NM_000088.3:r.-125_-123[4]");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_redundant_repeat_label_canonical_no_change() {
        let (out, hits) = correct_redundant_repeat_label("NM_000088.3:r.100_102[4]");
        assert_eq!(out, "NM_000088.3:r.100_102[4]");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_correct_redundant_repeat_label_dna_no_change() {
        // The spec note is RNA-specific; DNA `c.100_102CAG[4]` keeps its label.
        let (out, hits) = correct_redundant_repeat_label("NM_000088.3:c.100_102CAG[4]");
        assert_eq!(out, "NM_000088.3:c.100_102CAG[4]");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_correct_redundant_repeat_label_range_count() {
        let (out, hits) = correct_redundant_repeat_label("NM_000088.3:r.100_102cug[4_8]");
        assert_eq!(out, "NM_000088.3:r.100_102[4_8]");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_redundant_repeat_label_idempotent() {
        let input = "NM_000088.3:r.100_102cug[4]";
        let (out1, _) = correct_redundant_repeat_label(input);
        let (out2, hits2) = correct_redundant_repeat_label(&out1);
        assert_eq!(out1, out2);
        assert!(hits2.is_empty());
    }

    #[test]
    fn test_correct_redundant_repeat_label_does_not_strip_non_rna_repeat_after_r_description() {
        // Defensive regression: even though mixed `r.;c.` descriptions are not
        // canonical HGVS, a lowercase `acg[` run inside a later non-r. segment
        // must not be treated as a redundant RNA base label.
        let input = "r.100_102cug[4];c.50_52acg[3]";
        let (out, hits) = correct_redundant_repeat_label(input);
        assert_eq!(out, "r.100_102[4];c.50_52acg[3]");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].original, "cug");
    }

    // --- W4003 SinglePositionRange ---

    #[test]
    fn test_correct_single_position_range_del() {
        let (out, hits) = correct_single_position_range("NM_000088.3:c.123_123del");
        assert_eq!(out, "NM_000088.3:c.123del");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].error_type, ErrorType::SinglePositionRange);
        assert_eq!(hits[0].original, "123_123");
        assert_eq!(hits[0].corrected, "123");
    }

    #[test]
    fn test_correct_single_position_range_dup() {
        let (out, hits) = correct_single_position_range("NM_000088.3:c.123_123dup");
        assert_eq!(out, "NM_000088.3:c.123dup");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_single_position_range_inv() {
        let (out, hits) = correct_single_position_range("NM_000088.3:c.100_100inv");
        assert_eq!(out, "NM_000088.3:c.100inv");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_single_position_range_distinct_positions_no_hit() {
        for input in [
            "NM_000088.3:c.123_126del",
            "NM_000088.3:c.123_126dup",
            "NM_000088.3:c.100_102inv",
        ] {
            let (out, hits) = correct_single_position_range(input);
            assert_eq!(out, input);
            assert!(hits.is_empty(), "input {} should not warn", input);
        }
    }

    #[test]
    fn test_correct_single_position_range_negative_positions() {
        let (out, hits) = correct_single_position_range("NM_000088.3:c.-50_-50del");
        assert_eq!(out, "NM_000088.3:c.-50del");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_single_position_range_multiple_in_allele() {
        let (out, hits) = correct_single_position_range("NM_000088.3:c.[100_100del;200_200dup]");
        assert_eq!(out, "NM_000088.3:c.[100del;200dup]");
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn test_correct_single_position_range_idempotent() {
        let input = "NM_000088.3:c.123_123del";
        let (out1, _) = correct_single_position_range(input);
        let (out2, hits2) = correct_single_position_range(&out1);
        assert_eq!(out1, out2);
        assert!(hits2.is_empty());
    }

    #[test]
    fn test_correct_single_position_range_does_not_touch_ins() {
        // `ins` requires a range of two adjacent positions per spec; don't collapse.
        let (out, hits) = correct_single_position_range("NM_000088.3:c.100_100insATG");
        assert_eq!(out, "NM_000088.3:c.100_100insATG");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_correct_single_position_range_does_not_touch_substitution() {
        let (out, hits) = correct_single_position_range("NM_000088.3:c.100A>G");
        assert_eq!(out, "NM_000088.3:c.100A>G");
        assert!(hits.is_empty());
    }
}
