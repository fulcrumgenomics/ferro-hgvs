//! Correction logic for each error type.
//!
//! This module provides functions to detect and correct specific error types
//! in HGVS input strings.

use super::types::{ErrorType, ResolvedAction};
use crate::convert::mapper::CoordinateMapper;
use crate::hgvs::location::CdsPos;
use crate::reference::provider::ReferenceProvider;
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

/// Check if the input contains position zero on a numeric axis.
///
/// Returns `Some(byte_offset)` pointing at the `0` if found, `None`
/// otherwise.
///
/// Covered axes: `c.`, `g.`, `n.`, `r.`, `m.`, `o.` (1-based numbering
/// per HGVS `background/numbering.md` line 16 — "numbering starts with
/// c.1 at the A of the ATG translation initiation codon" — and the
/// analogous rule for the other coordinate systems).
///
/// **Excluded**: the protein axis `p.0` is the HGVS-defined "no protein
/// product" form (translation abolished; see
/// `recommendations/protein/description.md`); `p.0?` is the explicit
/// "possibly no product" variant. Both are spec-valid descriptions and
/// must NOT be flagged by this detector. Closes the `p.0` over-match
/// gap that issue #269 surfaced while resolving the W4002 audit row.
pub fn detect_position_zero(input: &str) -> Option<usize> {
    let bytes = input.as_bytes();
    let mut i = 0usize;
    while i + 1 < bytes.len() {
        if bytes[i] != b'.' {
            i += 1;
            continue;
        }
        // Only numeric axes; the byte immediately before `.` identifies
        // the axis marker (`c`/`g`/`n`/`r`/`m`/`o`). `p` is excluded
        // because `p.0` is the spec-defined "no protein product" form.
        //
        // Require an HGVS boundary before the axis letter so identifier-
        // embedded bytes (e.g. the `c` in `abc.0`) do not misfire: the
        // axis must sit at string start, or the byte before it must not
        // be an ASCII alphanumeric or `_` — that rules out accession-tail
        // letters/digits while still admitting real grammar separators
        // (`:`, `[`, `(`, `;`, `,`, `|`, whitespace).
        let axis = if i == 0 { 0 } else { bytes[i - 1] };
        let axis_has_boundary =
            i <= 1 || !(bytes[i - 2].is_ascii_alphanumeric() || bytes[i - 2] == b'_');
        let is_numeric_axis =
            axis_has_boundary && matches!(axis, b'c' | b'g' | b'n' | b'r' | b'm' | b'o');
        if !is_numeric_axis {
            i += 1;
            continue;
        }
        let zero_pos = i + 1;
        if bytes[zero_pos] != b'0' {
            i += 1;
            continue;
        }
        // Reject multi-digit `c.10`, `c.100`, … (the `0` is part of a
        // larger number, not a bare position zero).
        let next_is_digit = bytes
            .get(zero_pos + 1)
            .map(|b| b.is_ascii_digit())
            .unwrap_or(false);
        if next_is_digit {
            i += 1;
            continue;
        }
        return Some(zero_pos);
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
/// - `fsXN` (digits) or `fsX?` (uncertain) → `fsTerN` / `fsTer?` — emits
///   `DeprecatedFrameshiftX` (W3010).
/// - `fs*N` (digits) or `fs*?` (uncertain) → `fsTerN` / `fsTer?` — emits
///   `DeprecatedFrameshiftStar` (W3009).
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

        // Detect "fs*N", "fsXN", "fs*?", or "fsX?" (frameshift termination
        // notation, both known-position and uncertain-position forms).
        // Must begin at a literal "fs" boundary; the byte before the `f` is not
        // a lowercase letter (so we don't catch e.g. "ffs" — though that's
        // implausible in HGVS).
        if c == b'f'
            && i + 1 < bytes.len()
            && bytes[i + 1] == b's'
            && i + 2 < bytes.len()
            && (bytes[i + 2] == b'*' || bytes[i + 2] == b'X')
            && i + 3 < bytes.len()
            && (bytes[i + 3].is_ascii_digit() || bytes[i + 3] == b'?')
            && in_protein_segment(input, i + 2)
        {
            let star_or_x = bytes[i + 2];
            // Walk the trailing token: either a digit run ("23") or a single
            // "?" (uncertain termination position, VEP notation).
            let trail_start = i + 3;
            let (corrected_trail, trail_end) = if bytes[trail_start] == b'?' {
                ("?".to_string(), trail_start + 1)
            } else {
                let mut digits_end = trail_start;
                while digits_end < bytes.len() && bytes[digits_end].is_ascii_digit() {
                    digits_end += 1;
                }
                (input[trail_start..digits_end].to_string(), digits_end)
            };
            let original = &input[i + 2..trail_end]; // "*23", "X23", "*?", or "X?"
            let corrected = format!("Ter{}", corrected_trail);
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
                trail_end,
            ));
            // Emit "fs" + "Ter" + trail to output.
            out.push_str("fs");
            out.push_str(&corrected);
            i = trail_end;
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
/// Full-string corrector for mixed-case edit-type tokens.
///
/// Per HGVS recommendations/general.md lines 87-104, edit-type tokens
/// (`del`, `ins`, `dup`, `inv`, `delins`, `con`) are spelled in
/// lowercase. Inputs like `c.100Del`, `c.100_101INSATG`, or
/// `c.100_200CONNM_001:c.5_105` are rewritten by this corrector to the
/// canonical lowercase form and tagged with W1004.
///
/// Detection anchors the token after a digit, `]`, `)`, or `?` (the
/// shapes that legitimately precede an edit type in HGVS syntax) and
/// requires at least one uppercase letter in the token to fire.
/// Already-lowercase tokens are passed through. `delins` is checked
/// before `del` (longer prefix).
pub fn correct_edit_type_case_full(input: &str) -> (String, Vec<DetectedCorrection>) {
    let bytes = input.as_bytes();
    let mut corrections = Vec::new();
    let mut result = String::with_capacity(input.len());
    let mut i = 0usize;
    while i < bytes.len() {
        let prev_ok = i > 0
            && (matches!(bytes[i - 1], b']' | b')' | b'?' | b'*') || bytes[i - 1].is_ascii_digit());
        if prev_ok {
            // Try longest tokens first: `delins` (6), then 3-char ones.
            let three = bytes.get(i..i + 3);
            let six = bytes.get(i..i + 6);
            let mut matched: Option<(&'static str, usize)> = None;
            if let Some(s) = six {
                if eq_ignore_case(s, b"delins") {
                    matched = Some(("delins", 6));
                }
            }
            if matched.is_none() {
                if let Some(s) = three {
                    if eq_ignore_case(s, b"del") {
                        matched = Some(("del", 3));
                    } else if eq_ignore_case(s, b"ins") {
                        matched = Some(("ins", 3));
                    } else if eq_ignore_case(s, b"dup") {
                        matched = Some(("dup", 3));
                    } else if eq_ignore_case(s, b"inv") {
                        matched = Some(("inv", 3));
                    } else if eq_ignore_case(s, b"con") {
                        matched = Some(("con", 3));
                    }
                }
            }
            if let Some((canonical, len)) = matched {
                let raw = &input[i..i + len];
                if raw != canonical {
                    corrections.push(DetectedCorrection::new(
                        ErrorType::MixedCaseEditType,
                        raw,
                        canonical.to_string(),
                        i,
                        i + len,
                    ));
                }
                result.push_str(canonical);
                i += len;
                continue;
            }
        }
        let ch_end = next_char_end(bytes, i);
        result.push_str(&input[i..ch_end]);
        i = ch_end;
    }
    (result, corrections)
}

#[inline]
fn eq_ignore_case(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.eq_ignore_ascii_case(y))
}

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

/// One endpoint of a range. Captures the raw substring (for splicing on
/// rewrite) and a comparable sort key.
///
/// Layout: a position is `[-|*]?<digits>([+-]<digits>)?`. The main axis
/// value is signed: `-N` is negative (5'UTR / upstream), bare `N` is the
/// CDS integer, and `*N` is mapped to `CDS_END_SENTINEL + N` so that any
/// `*` position sorts strictly after any non-`*` position on the same
/// transcript. Offsets are signed integers that refine the position.
#[derive(Debug, Clone)]
struct EndpointToken {
    /// Original substring of the input from `start`..`end`.
    raw: String,
    /// Byte offset of the start of `raw` in the original input.
    start: usize,
    /// Byte offset of the end of `raw` (exclusive) in the original input.
    end: usize,
    /// Comparable sort key: `(main_axis, offset)`.
    sort_key: (i64, i64),
}

/// Sentinel added to `*N` positions to ensure they sort after any
/// non-`*` position. The HGVS coding axis tops out around `c.4548` on
/// COL1A1 (~5 kb genes go up to ~10 kb CDS); 10^12 is comfortably above
/// any conceivable real position.
const STAR_MARKER_SENTINEL: i64 = 1_000_000_000_000;

/// Parse one endpoint starting at `chars`. On success returns
/// `EndpointToken` and advances `chars` past the consumed bytes.
/// On any malformed input returns `None`.
fn parse_endpoint_token(
    after_marker: &str,
    base_offset: usize,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
) -> Option<EndpointToken> {
    let start_idx = chars.peek().map(|(i, _)| *i)?;
    let mut raw = String::new();
    let mut main_str = String::new();
    let mut is_star = false;
    let mut is_neg = false;

    // Optional class prefix: `-` (5'UTR / upstream) or `*` (3'UTR /
    // downstream). Bare digit means CDS body / transcript body.
    match chars.peek() {
        Some(&(_, '-')) => {
            is_neg = true;
            raw.push('-');
            chars.next();
        }
        Some(&(_, '*')) => {
            is_star = true;
            raw.push('*');
            chars.next();
        }
        _ => {}
    }

    // Main axis digits.
    while let Some(&(_, c)) = chars.peek() {
        if c.is_ascii_digit() {
            raw.push(c);
            main_str.push(c);
            chars.next();
        } else {
            break;
        }
    }
    if main_str.is_empty() {
        return None;
    }
    let main_val: i64 = main_str.parse().ok()?;
    let signed_main = if is_neg { -main_val } else { main_val };
    let class_main = if is_star {
        STAR_MARKER_SENTINEL + main_val
    } else {
        signed_main
    };

    // Optional offset: `+N` or `-N`.
    let mut offset: i64 = 0;
    if let Some(&(_, c)) = chars.peek() {
        if c == '+' || c == '-' {
            let sign = c;
            raw.push(c);
            chars.next();
            let mut off_str = String::new();
            while let Some(&(_, oc)) = chars.peek() {
                if oc.is_ascii_digit() {
                    raw.push(oc);
                    off_str.push(oc);
                    chars.next();
                } else {
                    break;
                }
            }
            if off_str.is_empty() {
                // Trailing `+` / `-` with no digits — malformed; bail.
                return None;
            }
            let off_val: i64 = off_str.parse().ok()?;
            offset = if sign == '-' { -off_val } else { off_val };
        }
    }

    let end_idx = chars.peek().map(|(i, _)| *i).unwrap_or(after_marker.len());
    Some(EndpointToken {
        raw,
        start: base_offset + start_idx,
        end: base_offset + end_idx,
        sort_key: (class_main, offset),
    })
}

/// Detect swapped interval positions.
///
/// Looks for patterns like `c.200_100del` where start > end, including
/// offset-bearing forms (`c.100+5_99+3del`) and 3'UTR `*N` markers
/// (`c.*5_*1del`). Returns Some with the detected swap information and
/// the swapped endpoint pair (with offsets and markers preserved).
///
/// Axis-aware: the `m.` (mitochondrial) and `o.` (circular) axes are skipped
/// entirely. On a circular contig a range with `start > end` is a valid
/// origin-crossing wraparound (SVD-WG006), not a transposition error, so
/// rewriting it would destroy the wraparound signal. See issue #467.
pub fn detect_swapped_positions(input: &str) -> Option<DetectedCorrection> {
    // Find the coordinate marker and remember which axis it names so circular
    // axes (m./o.) can be skipped.
    let coord_markers = [
        (".c.", 'c'),
        (".g.", 'g'),
        (".n.", 'n'),
        (".m.", 'm'),
        (".o.", 'o'),
        (".r.", 'r'),
    ];
    let mut search_start = 0;
    let mut axis = None;
    // Pick the earliest (leftmost) marker in the input, not the first marker
    // kind in `coord_markers`. A later embedded coord prefix (e.g. a source
    // HGVS appended after the main variant) must not override the axis of the
    // leading expression, which would parse the wrong range. See issue #467.
    if let Some((pos, marker_len, marker_axis)) = coord_markers
        .iter()
        .filter_map(|(marker, marker_axis)| {
            input
                .find(marker)
                .map(|pos| (pos, marker.len(), *marker_axis))
        })
        .min_by_key(|(pos, _, _)| *pos)
    {
        search_start = pos + marker_len;
        axis = Some(marker_axis);
    }
    if search_start == 0 {
        if let Some(pos) = input.find(['c', 'g', 'n', 'm', 'o', 'r']) {
            if input.get(pos + 1..pos + 2) == Some(".") {
                search_start = pos + 2;
                axis = input[pos..].chars().next();
            }
        }
    }
    if search_start == 0 {
        return None;
    }

    // Skip circular axes: a reversed range here is a valid origin-crossing
    // wraparound, not a swapped-position error.
    if matches!(axis, Some('m') | Some('o')) {
        return None;
    }

    let after_marker = &input[search_start..];
    let mut chars = after_marker.char_indices().peekable();
    let first = parse_endpoint_token(after_marker, search_start, &mut chars)?;
    if chars.next().map(|(_, c)| c) != Some('_') {
        return None;
    }
    let second = parse_endpoint_token(after_marker, search_start, &mut chars)?;

    if first.sort_key > second.sort_key {
        return Some(DetectedCorrection::new(
            ErrorType::SwappedPositions,
            format!("{}_{}", first.raw, second.raw),
            format!("{}_{}", second.raw, first.raw),
            first.start,
            second.end,
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

/// Detect and correct the deprecated multi-base substitution syntax
/// (`c.79_80GC>TT`, `c.79GC>TT`, `c.100_102>ATG`, and the lowercase RNA
/// equivalents `r.79_80gc>uu` etc.) by rewriting it to `delins` form, per
/// HGVS spec recommendations/DNA/substitution.md line 26.
///
/// Single-base substitutions (`c.100A>G`, `r.100a>g`) are canonical and
/// not touched.
///
/// The detector recognises three patterns inside the input string
/// (IUPAC matching is case-insensitive, so the same patterns apply to
/// uppercase DNA and lowercase RNA inputs):
///
/// 1. `<pos>_<pos>[IUPAC]+>[IUPAC]+` — explicit ref bases, range form
/// 2. `<pos>[IUPAC]{2,}>[IUPAC]+`    — single position, multi-base ref
/// 3. `<pos>(_<pos>)?>[IUPAC]+`      — no ref bases
///
/// Case is preserved on output: the rewrite slices bases from the
/// original input rather than any normalised buffer.
///
/// Returns the (possibly rewritten) string and the list of corrections.
pub fn correct_old_substitution_syntax(input: &str) -> (String, Vec<DetectedCorrection>) {
    let bytes = input.as_bytes();
    let mut result = String::with_capacity(input.len());
    let mut corrections = Vec::new();
    let mut i = 0usize;

    while i < bytes.len() {
        // The previous character must NOT be an ASCII alphanumeric or `_`,
        // so we don't match digits embedded in identifiers like `NM_000088`.
        let prev_ok = i == 0
            || !{
                let p = bytes[i - 1];
                (p as char).is_ascii_alphanumeric() || p == b'_'
            };

        // Walk a position token: optional `-`/`*`, then digits.
        let pos1_start = i;
        let mut j = i;
        if j < bytes.len() && matches!(bytes[j] as char, '-' | '*') {
            j += 1;
        }
        let digits_start = j;
        while j < bytes.len() && (bytes[j] as char).is_ascii_digit() {
            j += 1;
        }
        if j == digits_start || !prev_ok {
            // Not a real position at this offset; emit one byte and advance.
            result.push(bytes[i] as char);
            i += 1;
            continue;
        }
        let pos1_end = j;

        // Optional `_<pos2>` for a range.
        let has_range = bytes.get(j).copied() == Some(b'_');
        let mut pos2_end = j;
        if has_range {
            let mut k = j + 1;
            if k < bytes.len() && matches!(bytes[k] as char, '-' | '*') {
                k += 1;
            }
            let d2 = k;
            while k < bytes.len() && (bytes[k] as char).is_ascii_digit() {
                k += 1;
            }
            if k > d2 {
                pos2_end = k;
            }
        }
        let positions_end = if has_range && pos2_end > j {
            pos2_end
        } else {
            j
        };

        // Now look for `[ACGTU...]*>[ACGTU...]+`.
        let mut k = positions_end;
        let refs_start = k;
        while k < bytes.len() && is_iupac_base(bytes[k] as char) {
            k += 1;
        }
        let refs_end = k;
        let ref_count = refs_end - refs_start;

        if k >= bytes.len() || bytes[k] != b'>' {
            // Not a substitution form — emit the position token verbatim and
            // advance past it.
            result.push_str(&input[pos1_start..positions_end]);
            i = positions_end;
            continue;
        }
        let arrow = k;

        // RHS: 1+ IUPAC bases.
        let mut m = arrow + 1;
        while m < bytes.len() && is_iupac_base(bytes[m] as char) {
            m += 1;
        }
        let rhs_count = m - (arrow + 1);
        if rhs_count == 0 {
            result.push_str(&input[pos1_start..positions_end]);
            i = positions_end;
            continue;
        }

        // Canonical sub: single-pos with exactly 1 ref base and 1 RHS base.
        let is_canonical_sub = !has_range && ref_count == 1 && rhs_count == 1;
        if is_canonical_sub {
            result.push_str(&input[pos1_start..m]);
            i = m;
            continue;
        }

        // It's the deprecated form. Build the rewrite.
        let pos_str = &input[pos1_start..positions_end];
        let rhs_str = &input[arrow + 1..m];
        let original = &input[pos1_start..m];

        let new_pos = if has_range {
            pos_str.to_string()
        } else if ref_count > 1 {
            // Single-pos with multi-base ref: extend to a range. Only safe
            // when pos1 is a simple positive integer; otherwise leave as-is.
            let pos1_text = &input[pos1_start..pos1_end];
            let pos1_num: Option<i64> = pos1_text.parse().ok();
            match pos1_num {
                Some(n) if n >= 1 => {
                    format!("{}_{}", n, n + ref_count as i64 - 1)
                }
                _ => pos1_text.to_string(),
            }
        } else {
            pos_str.to_string()
        };

        let corrected = format!("{}delins{}", new_pos, rhs_str);
        corrections.push(DetectedCorrection::new(
            ErrorType::OldSubstitutionSyntax,
            original.to_string(),
            corrected.clone(),
            pos1_start,
            m,
        ));
        result.push_str(&corrected);
        i = m;
    }

    (result, corrections)
}

/// IUPAC nucleotide base (incl. ambiguity codes used in HGVS).
///
/// Case-insensitive. HGVS uses uppercase letters in DNA descriptions
/// (`g.`/`c.`/`n.`/`m.`) and lowercase letters in RNA descriptions
/// (`r.`); both denote the same underlying nucleotide alphabet, so this
/// predicate accepts either form. Callers that need to preserve case
/// in their output must slice from the original input rather than
/// transforming through this predicate.
fn is_iupac_base(c: char) -> bool {
    matches!(
        c.to_ascii_uppercase(),
        'A' | 'C'
            | 'G'
            | 'T'
            | 'U'
            | 'R'
            | 'Y'
            | 'S'
            | 'W'
            | 'K'
            | 'M'
            | 'B'
            | 'D'
            | 'H'
            | 'V'
            | 'N'
    )
}

/// Detect retracted `c.IVS<n>(+|-)<offset>...` / `n.IVS<n>...` /
/// `r.IVS<n>...` intronic notation. Per HGVS spec
/// background/numbering.md line 32 the form is retracted; ferro cannot
/// auto-rewrite it without genomic intron metadata so the detector returns
/// the positions but does not mutate the input.
///
/// Returns a list of `DetectedCorrection`s with `original = "IVSn"` and an
/// empty `corrected` field (rewrite is not possible). Callers that warn /
/// reject should consume `error_type` + `start`/`end` for diagnostic spans.
pub fn detect_deprecated_ivs(input: &str) -> Vec<DetectedCorrection> {
    let mut out = Vec::new();
    let bytes = input.as_bytes();
    // Axis prefixes (`c.`, `n.`, `r.`) are matched literal (HGVS axis
    // markers are always lowercase), but the retracted `IVS` token is
    // matched case-insensitively so that lowercase RNA inputs like
    // `r.ivs5+1g>a` are detected on parity with the uppercase form.
    let axes: &[&[u8]] = &[b"c.", b"n.", b"r."];
    for axis in axes {
        let mut search_start = 0usize;
        while let Some(rel) = find_subslice(&bytes[search_start..], axis) {
            let prefix_start = search_start + rel;
            let ivs_start = prefix_start + 2;
            // Need at least "IVS<digit>" past the axis marker.
            if ivs_start + 4 > bytes.len() {
                break;
            }
            // Validate left boundary: char before the axis must NOT be
            // an ASCII alphanumeric (avoid matching embedded identifiers).
            let ok_left = prefix_start == 0 || !bytes[prefix_start - 1].is_ascii_alphanumeric();
            let ok_ivs = bytes[ivs_start..ivs_start + 3].eq_ignore_ascii_case(b"IVS");
            let mut j = ivs_start + 3;
            let digits_start = j;
            while j < bytes.len() && (bytes[j] as char).is_ascii_digit() {
                j += 1;
            }
            if ok_left && ok_ivs && j > digits_start {
                out.push(DetectedCorrection::new(
                    ErrorType::DeprecatedIvsNotation,
                    &input[ivs_start..j],
                    "",
                    ivs_start,
                    j,
                ));
                search_start = j;
            } else {
                search_start = prefix_start + 1;
            }
        }
    }
    // Sort by start so multi-axis scans return ordered results.
    out.sort_by_key(|c| c.start);
    out
}

/// Find the first occurrence of `needle` inside `haystack`, byte-wise.
fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}

/// Detect a protein-edit insertion whose inserted sequence is given as a
/// bracketed amino-acid list (e.g. `p.Arg97_Trp98ins[Ala;Pro]`).
///
/// HGVS v21's protein insertion notation concatenates 3-letter codes
/// without separators (`p.Arg97_Trp98insAlaPro`). The bracketed
/// `[Ala;Pro]` form is not in the spec — brackets at the edit level
/// would collide with the variant-level allele syntax. ferro rejects
/// the form with W3021 and points at the canonical concatenated shape.
///
/// Scope: only the protein description (text after `p.`) is scanned, so
/// allele/repeat brackets in `c.`, `g.`, `n.`, `r.`, `m.` descriptions
/// are not falsely triggered. Multi-variant inputs (compound alleles)
/// scan each `p.` slice independently.
///
/// Returns one `DetectedCorrection` per offending `ins[...]` span. The
/// `corrected` field carries the best-effort canonical rewrite when the
/// bracket contents parse as `<Aa3>;<Aa3>;...`; otherwise it is empty
/// (and the diagnostic still names the canonical shape via the
/// preprocessor hint). The `start`/`end` span covers the `ins[...]`
/// substring (inclusive of the closing `]`).
pub fn detect_protein_bracketed_aa_insertion(input: &str) -> Vec<DetectedCorrection> {
    let mut out = Vec::new();
    let bytes = input.as_bytes();

    // Bare-prefix form `p.…` with no accession (e.g. `p.Arg97_Trp98ins[Ala;Pro]`).
    // The full `<acc>:p.` loop below only matches `:p.`, so we'd otherwise
    // miss this shape. Scoped tightly: only when `p.` is the very first
    // token, to avoid colliding with substrings like `op.` mid-input.
    if bytes.starts_with(b"p.") {
        scan_protein_body(input, 2, &mut out);
    }

    let mut search_start = 0usize;
    while let Some(rel) = find_subslice(&bytes[search_start..], b":p.") {
        let p_dot_start = search_start + rel;
        // Body of the protein description starts after ":p."
        let body_start = p_dot_start + 3;
        let body_end = scan_protein_body(input, body_start, &mut out);
        search_start = body_end;
    }

    out
}

/// Scan the protein-edit body that begins at `body_start` for the
/// `ins[…]` shape, appending one `DetectedCorrection` per offending
/// span. Returns the body end (next whitespace / comma / EOF), so the
/// caller can advance its outer search cursor past the body without
/// rescanning. Shared between the `:p.` and bare-`p.` entry points.
fn scan_protein_body(input: &str, body_start: usize, out: &mut Vec<DetectedCorrection>) -> usize {
    let bytes = input.as_bytes();
    // Body ends at next whitespace, comma (multi-variant separator), or
    // end of input. This bounds the scan so we don't bleed into a sibling
    // description in compound inputs. Note: `;` is intentionally NOT a
    // terminator here — it is the intra-bracket AA separator inside
    // `ins[Ala;Pro]`, which is exactly the shape we are looking for.
    // Treating `;` as a body terminator would clip the bracket and miss
    // the detection.
    let mut body_end = body_start;
    while body_end < bytes.len()
        && !(bytes[body_end] as char).is_whitespace()
        && bytes[body_end] != b','
    {
        body_end += 1;
    }

    let mut j = body_start;
    while j + 4 <= body_end {
        if &bytes[j..j + 4] == b"ins["
            // `ins[` is the bracketed-AA form only when `ins` is a
            // protein edit token, not part of a longer word. Reject
            // if preceded by an alphabetic byte (so `delins[…]`
            // remains untouched — that's a different edit type
            // outside this issue's scope).
            && (j == body_start || !bytes[j - 1].is_ascii_alphabetic())
        {
            // Find the matching close `]` within the body.
            let open = j + 3; // index of '['
            let close_search_start = open + 1;
            let mut k = close_search_start;
            let mut close = None;
            while k < body_end {
                if bytes[k] == b']' {
                    close = Some(k);
                    break;
                }
                k += 1;
            }
            if let Some(close_idx) = close {
                let inner = &input[close_search_start..close_idx];
                let canonical = canonicalize_protein_bracketed_aa_list(inner);
                let original = &input[j..close_idx + 1];
                let corrected = match canonical {
                    Some(ref cat) => format!("ins{}", cat),
                    None => String::new(),
                };
                out.push(DetectedCorrection::new(
                    ErrorType::ProteinBracketedAaInsertion,
                    original,
                    corrected,
                    j,
                    close_idx + 1,
                ));
                j = close_idx + 1;
                continue;
            }
        }
        j += 1;
    }

    body_end
}

/// Try to canonicalize the contents of a `[...]` AA list to a
/// concatenated 3-letter sequence. Returns `Some("AlaPro")` for a
/// well-formed `Ala;Pro` body; `None` if the body is empty, has a
/// leading/trailing separator, an empty segment, or any segment that is
/// not a recognized 3-letter amino acid (case-sensitive: HGVS canonical
/// 3-letter codes start with an uppercase letter).
fn canonicalize_protein_bracketed_aa_list(inner: &str) -> Option<String> {
    if inner.is_empty() {
        return None;
    }
    if inner.starts_with(';') || inner.ends_with(';') {
        return None;
    }
    let mut out = String::with_capacity(inner.len());
    for segment in inner.split(';') {
        if segment.is_empty() {
            return None;
        }
        if !is_canonical_three_letter_aa(segment) {
            return None;
        }
        out.push_str(segment);
    }
    Some(out)
}

/// Detect the deprecated `con` (sequence conversion) edit syntax and rewrite
/// it as `delins`. Per HGVS spec recommendations/DNA/delins.md line 19:
/// "The previous format 'con' is no longer used".
///
/// Recognises `<pos>_<pos>con<source>` where `<pos>` is an HGVS position
/// (digits, optional leading `-` or `*`, optional `+`/`-` offset). The
/// `<source>` portion runs to end of input or whitespace; the rewrite is
/// `<pos>_<pos>delins<source>`.
pub fn correct_deprecated_con(input: &str) -> (String, Vec<DetectedCorrection>) {
    let bytes = input.as_bytes();
    let mut result = String::with_capacity(input.len() + 4);
    let mut corrections = Vec::new();
    let mut i = 0usize;

    while i < bytes.len() {
        if i + 3 <= bytes.len() && &bytes[i..i + 3] == b"con" {
            let con_start = i;
            // Must be immediately preceded by an ASCII digit (end of pos2).
            let preceded_by_digit =
                con_start > 0 && (bytes[con_start - 1] as char).is_ascii_digit();
            if preceded_by_digit && has_range_position_before(bytes, con_start) {
                // Right-hand source: take until whitespace or end, but only
                // when the first byte starts a valid HGVS source. Without
                // this guard, `c.100_200conditional` would be misread as
                // `delinsditional` (everything until whitespace).
                let src_start = con_start + 3;
                if !has_valid_con_source_start(bytes, src_start) {
                    result.push(bytes[i] as char);
                    i += 1;
                    continue;
                }
                let mut j = src_start;
                while j < bytes.len() && !(bytes[j] as char).is_whitespace() {
                    j += 1;
                }
                let src_end = j;
                if src_end > src_start {
                    let src = &input[src_start..src_end];
                    let original = &input[con_start..src_end];
                    let corrected = format!("delins{}", src);
                    corrections.push(DetectedCorrection::new(
                        ErrorType::DeprecatedConSyntax,
                        original.to_string(),
                        corrected.clone(),
                        con_start,
                        src_end,
                    ));
                    result.push_str(&corrected);
                    i = src_end;
                    continue;
                }
            }
        }
        result.push(bytes[i] as char);
        i += 1;
    }

    (result, corrections)
}

/// True if the bytes preceding `end` end in an HGVS-position range of the
/// form `<pos>_<pos>` (with optional `-`/`*` prefixes and `+`/`-` offsets).
/// Conservatively requires the `_` separator — single-position `con` is not
/// a recognised pattern.
fn has_range_position_before(bytes: &[u8], end: usize) -> bool {
    let mut i = end;
    // digits at the right end (the second pos)
    let mut saw_digit = false;
    while i > 0 && (bytes[i - 1] as char).is_ascii_digit() {
        i -= 1;
        saw_digit = true;
    }
    if !saw_digit {
        return false;
    }
    // optional sign for the right-hand pos (`+`/`-` offset, or `-`/`*` UTR)
    if i > 0 && matches!(bytes[i - 1] as char, '+' | '-' | '*') {
        i -= 1;
    }
    // optional more digits (offset start)
    while i > 0 && (bytes[i - 1] as char).is_ascii_digit() {
        i -= 1;
    }
    if i > 0 && matches!(bytes[i - 1] as char, '-' | '*') {
        i -= 1;
    }
    // require `_`
    if i == 0 || bytes[i - 1] != b'_' {
        return false;
    }
    i -= 1;
    // left-hand pos: digits
    let mut saw_digit_l = false;
    while i > 0 && (bytes[i - 1] as char).is_ascii_digit() {
        i -= 1;
        saw_digit_l = true;
    }
    if !saw_digit_l {
        return false;
    }
    if i > 0 && matches!(bytes[i - 1] as char, '+' | '-' | '*') {
        i -= 1;
    }
    while i > 0 && (bytes[i - 1] as char).is_ascii_digit() {
        i -= 1;
    }
    if i > 0 && matches!(bytes[i - 1] as char, '-' | '*') {
        i -= 1;
    }
    let _ = i;
    true
}

/// True if the byte at `start` is a plausible first character of a `delins`
/// source. Allowed starters mirror what HGVS lets follow `delins`:
/// - ASCII digit (bare-position source like `100_200`)
/// - `-`/`*`/`(` (UTR offsets and uncertainty)
/// - An accession start: ASCII uppercase letter (e.g. `NM_001`, `NC_000022`)
/// - A coordinate-prefix shorthand: lowercase `g`/`c`/`n`/`r`/`m` immediately
///   followed by `.`
///
/// Rejecting other starters prevents `c.100_200conditional` and similar
/// embedded-`con` words from being rewritten to `delinsditional`.
fn has_valid_con_source_start(bytes: &[u8], start: usize) -> bool {
    if start >= bytes.len() {
        return false;
    }
    let c = bytes[start] as char;
    if c.is_ascii_digit() || matches!(c, '-' | '*' | '(') || c.is_ascii_uppercase() {
        return true;
    }
    // Lowercase coordinate prefix must be one of g/c/n/r/m followed by `.`.
    if matches!(c, 'g' | 'c' | 'n' | 'r' | 'm')
        && start + 1 < bytes.len()
        && bytes[start + 1] == b'.'
    {
        return true;
    }
    false
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
/// Detect mismatch between the position-range length and the explicit
/// reference-sequence length on `del` / `dup` / `inv` / `delins<ref>ins<alt>`.
///
/// HGVS range syntax implies `len(ref) == end - start + 1`. ferro
/// silently accepts mismatches today (`g.100_110delAAAATTTGCC` has
/// range 11 but 10 bases). This detector flags them so the preprocessor
/// can emit W3016.
///
/// Only handles ranges with simple integer endpoints (no offsets,
/// no `-N` / `*N` markers, no `?`). Skips inputs where any endpoint
/// involves an offset or marker — computing length requires a provider
/// in those cases.
pub fn detect_length_mismatch(input: &str) -> Vec<DetectedCorrection> {
    detect_length_mismatch_inner(input, None)
}

/// Provider-aware variant of [`detect_length_mismatch`]. Handles the
/// mixed-shape intronic endpoint cases that the no-provider scan
/// deliberately skips:
///
/// - different anchor bases (`c.100+5_200-3del...`),
/// - opposite-sign offsets at the same anchor (`c.100+5_100-3del...`,
///   intron-internal range across the anchor),
/// - one exonic + one intronic endpoint (`c.100_200+5del...`).
///
/// Resolution path: parse the accession prefix of each ranged segment
/// (`<ACC>:c.X_Ydel...`), look up the transcript via the provider, and
/// resolve each endpoint to a genomic position via
/// [`crate::convert::mapper::CoordinateMapper::cds_to_genomic_with_intron`].
/// The reference-frame span is then `|g_end - g_start| + 1` (absolute
/// value handles minus-strand transcripts where genomic order is
/// reversed relative to transcript order).
///
/// Falls back to the no-provider behavior when:
/// - the input has no accession prefix (compound-allele inner members
///   without a fresh accession),
/// - the provider lookup fails,
/// - either endpoint can't be resolved to a genomic position
///   (transcript missing genomic coords, offset outside any intron),
/// - the axis is anything other than `c.` (only c. has the
///   exon/intron coordinate system that benefits from this path; g./m./o.
///   ranges already use the exonic same-axis branch).
///
/// #429 (follow-up to #390 item 5).
pub fn detect_length_mismatch_with_provider(
    input: &str,
    provider: &dyn ReferenceProvider,
) -> Vec<DetectedCorrection> {
    detect_length_mismatch_inner(input, Some(provider))
}

fn detect_length_mismatch_inner(
    input: &str,
    provider: Option<&dyn ReferenceProvider>,
) -> Vec<DetectedCorrection> {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();
    if !has_non_protein_description(bytes) {
        return hits;
    }

    // Walk the input looking for any of the coord markers (`c.`, `g.`,
    // `n.`, `m.`, `o.`, `r.`) preceded by `:`, `[`, `(`, `;`, or start of
    // string. Compound-allele inner members (digits immediately after
    // `[` or `;` — they inherit the outer coord axis, no second `c.`
    // prefix appears) are also scanned so `c.[100A>G;200_205delAAAA]`
    // catches the inner `200_205delAAAA` (#390 item 4). Multiple ranges
    // per input are handled by continuing past each match.
    //
    // `last_accession_span` tracks the most recent accession prefix
    // seen (the `ACC:` before a coord marker); `last_axis` tracks the
    // axis byte at that marker. Compound-allele inner members inherit
    // both. Used for the provider-aware path (#429) — without correct
    // axis tracking, a `g.[X;Y]` compound's inner `Y` would be
    // erroneously mis-routed through the c.-axis provider path.
    let mut i = 0usize;
    let mut last_accession_span: Option<(usize, usize)> = None;
    let mut last_axis: Option<u8> = None;
    while i + 1 < bytes.len() {
        let prev_ok = i == 0 || matches!(bytes[i - 1], b':' | b'(' | b'[' | b';');
        let is_coord =
            matches!(bytes[i], b'c' | b'g' | b'n' | b'm' | b'o' | b'r') && bytes[i + 1] == b'.';
        if prev_ok && is_coord {
            // If the coord marker was preceded by `:`, capture the
            // accession span ending at that `:`. Walk back over
            // non-separator bytes to find the accession start.
            if i > 0 && bytes[i - 1] == b':' {
                let acc_end = i - 1;
                let mut acc_start = acc_end;
                while acc_start > 0
                    && !matches!(
                        bytes[acc_start - 1],
                        b':' | b'(' | b'[' | b';' | b' ' | b'\t' | b'\n'
                    )
                {
                    acc_start -= 1;
                }
                if acc_start < acc_end {
                    last_accession_span = Some((acc_start, acc_end));
                }
            }
            let axis = bytes[i];
            last_axis = Some(axis);
            let start_idx = i + 2;
            let accession = last_accession_span.map_or("", |(s, e)| &input[s..e]);
            if let Some(hit) =
                try_detect_length_mismatch_at(bytes, input, start_idx, axis, provider, accession)
            {
                hits.push(hit);
            }
            i = start_idx;
            continue;
        }
        // Compound-allele inner-member fast path: a digit right after
        // `[` or `;` is the start of a range that inherits the outer
        // axis + accession (no fresh `<axis>.` marker on the inner
        // member). `try_detect_length_mismatch_at` returns `None` for
        // anything that isn't a `<int>_<int><edit><refseq>` pattern,
        // so over-triggering is safe. The `i > 0` guard prevents an
        // underflow on the `bytes[i - 1]` index at the very first byte.
        //
        // The axis is the outer axis captured above (`last_axis`).
        // Falling back to a hard-coded `b'c'` would mis-route a
        // `g.[X;Y]` compound's inner `Y` into the c.-axis provider
        // path. The axis-guard inside `try_detect_length_mismatch_at`
        // ensures provider resolution only fires on the spec-supported
        // c. axis; non-`c.` outer axes therefore correctly skip the
        // provider path on inner members too.
        let prev_is_compound_separator = i > 0 && matches!(bytes[i - 1], b'[' | b';');
        if prev_is_compound_separator && bytes[i].is_ascii_digit() {
            let axis = last_axis.unwrap_or(b'c');
            let accession = last_accession_span.map_or("", |(s, e)| &input[s..e]);
            if let Some(hit) =
                try_detect_length_mismatch_at(bytes, input, i, axis, provider, accession)
            {
                hits.push(hit);
            }
        }
        i += 1;
    }

    hits
}

/// Try to detect a length-mismatch at `pos` (right after a coord marker
/// like `.c.`). Returns the detection on success.
///
/// `axis` is the coord-system byte (`b'c'`, `b'g'`, etc.) so the
/// provider-aware path knows whether to invoke `c.`-specific
/// resolution. `accession` is the prefix accession (e.g. `NM_000088.3`)
/// — empty when none is present; the provider-aware path only fires
/// when both `provider` and a non-empty `accession` are available.
fn try_detect_length_mismatch_at(
    bytes: &[u8],
    input: &str,
    pos: usize,
    axis: u8,
    provider: Option<&dyn ReferenceProvider>,
    accession: &str,
) -> Option<DetectedCorrection> {
    // Endpoint 1: bare integer with optional `+offset` / `-offset`.
    let (start_base, start_off, mut j) = parse_endpoint_with_offset(bytes, pos)?;
    // Underscore.
    if j >= bytes.len() || bytes[j] != b'_' {
        return None;
    }
    j += 1;
    // Endpoint 2: same shape.
    let (end_base, end_off, mut k) = parse_endpoint_with_offset(bytes, j)?;

    // Compute the range length on the supported endpoint shapes:
    //   (a) Both exonic (no offset on either side): `end - start + 1`.
    //   (b) Both intronic with the same anchor base and same-sign
    //       offsets (`c.100+5_100+10` or `c.100-5_100-2`): the span
    //       is linear inside the intron, length = `|off2 - off1| + 1`.
    //   (c) Mixed-shape (different anchors / opposite signs / one
    //       exonic + one intronic): provider-aware resolution via the
    //       transcript's c→g mapper. Only fires on the `c.` axis when
    //       `provider` and `accession` are both available. #429
    //       (follow-up to #390 item 5).
    // Without a provider, mixed-shape inputs are skipped (no false
    // positives, no panic).
    let range_len = match (start_off, end_off) {
        (0, 0) => {
            if end_base < start_base {
                return None; // W4001's job.
            }
            (end_base - start_base + 1) as usize
        }
        (so, eo) if start_base == end_base && so.signum() == eo.signum() => {
            if eo < so {
                return None; // W4001's job.
            }
            (eo - so + 1) as usize
        }
        _ => {
            // Mixed-shape: try provider-aware resolution. Only the
            // `c.` axis carries the exon/intron coordinate system that
            // the c→g mapper resolves; other axes fall through to the
            // skip (per the issue's scoping).
            if axis != b'c' || accession.is_empty() {
                return None;
            }
            let provider = provider?;
            let g_start =
                resolve_cds_endpoint_to_genomic(provider, accession, start_base, start_off)?;
            let g_end = resolve_cds_endpoint_to_genomic(provider, accession, end_base, end_off)?;
            // Absolute difference handles both plus-strand (g_end >
            // g_start) and minus-strand (g_start > g_end) transcripts;
            // the W3016 length check only cares about the reference-frame
            // span magnitude.
            (g_start.max(g_end) - g_start.min(g_end) + 1) as usize
        }
    };

    // Identify the edit keyword. Order matters: `delins` is a prefix of
    // `del`, so check the longer first.
    let edit_kind = match () {
        _ if bytes.len() >= k + 6 && &bytes[k..k + 6] == b"delins" => EditKind::Delins,
        _ if bytes.len() >= k + 3 && &bytes[k..k + 3] == b"del" => EditKind::Del,
        _ if bytes.len() >= k + 3 && &bytes[k..k + 3] == b"dup" => EditKind::Dup,
        _ if bytes.len() >= k + 3 && &bytes[k..k + 3] == b"inv" => EditKind::Inv,
        _ => return None,
    };
    let keyword_end = k + match edit_kind {
        EditKind::Delins => 6,
        _ => 3,
    };
    k = keyword_end;

    // For del/dup/inv: the ref sequence (if any) is alphabetic bases
    // immediately following the keyword. For delins: the ref sequence
    // (if any) appears between `del` and `ins`, i.e. for
    // `delins<ref>` only the inserted seq is given. To compare lengths
    // for delins we need the deleted-ref form `del<ref>ins<alt>`. Look
    // it up by re-scanning: this branch only fires when the input is
    // the explicit form `del<ref>ins<alt>` (legacy/clinvar).
    let ref_seq = match edit_kind {
        EditKind::Del | EditKind::Dup | EditKind::Inv => take_ref_seq_run(bytes, k),
        EditKind::Delins => {
            // The simple `delins<ins>` form has no `del<ref>` segment to
            // measure; we'd be measuring the *inserted* seq against the
            // range length, which is the wrong rule. Skip.
            return None;
        }
    };
    if !ref_seq.is_empty() {
        let end_of_seq = k + ref_seq.len();
        if range_len == ref_seq.len() {
            return None;
        }
        // Mismatch: emit a detection covering the position range + edit
        // keyword + ref seq.
        let span_text = &input[pos..end_of_seq];
        return Some(DetectedCorrection::new(
            ErrorType::LengthMismatch,
            span_text,
            String::new(),
            pos,
            end_of_seq,
        ));
    }

    // No alphabetic ref seq run — check for a numeric length suffix
    // (e.g. `c.100_102del3`, `c.100_101dup3`, `c.100_104inv3`). #439
    // (follow-up to #427): when the explicit `length` field is
    // populated by the parser but disagrees with the position-interval
    // span, surface W3016 just like a disagreeing explicit ref-seq
    // would. The pre-#439 detector silently let these through because
    // `take_ref_seq_run` (alphabetic-only) returned empty.
    //
    // Scoped to `del`/`dup`/`inv` here; `delins` is handled in #394
    // item 2 by the validator that rejects `delins` with both
    // `deleted` and `deleted_length` set. `detect_del_size_suffix`
    // (W3011) only flags the deprecated-syntax aspect of `del<N>`
    // (not `dup<N>` / `inv<N>` — see `detect_del_size_suffix`'s
    // `b"del"`-only match); the W3016 numeric-suffix length check
    // here applies uniformly across all three.
    let digits_start = k;
    let mut digits_end = digits_start;
    while digits_end < bytes.len() && bytes[digits_end].is_ascii_digit() {
        digits_end += 1;
    }
    if digits_end == digits_start {
        // Neither alphabetic ref seq nor numeric length — nothing to
        // compare against `range_len`.
        return None;
    }
    let declared_len: usize = std::str::from_utf8(&bytes[digits_start..digits_end])
        .ok()?
        .parse()
        .ok()?;
    if range_len == declared_len {
        return None;
    }
    let span_text = &input[pos..digits_end];
    Some(DetectedCorrection::new(
        ErrorType::LengthMismatch,
        span_text,
        String::new(),
        pos,
        digits_end,
    ))
}

/// Variant: scan also covers explicit `del<ref>ins<alt>` (legacy
/// ClinVar shape). Public so the preprocessor can call it once across
/// the input.
#[derive(Copy, Clone)]
enum EditKind {
    Del,
    Dup,
    Inv,
    Delins,
}

/// Parse a coord endpoint at `pos` as `<base>[+<offset>|-<offset>]`,
/// returning `(base, offset, end_idx)`. `offset` is `0` when no
/// `+`/`-` suffix follows. Refuses inputs that start with `-`, `*`,
/// `+`, or anything other than a digit (#390 item 5).
fn parse_endpoint_with_offset(bytes: &[u8], pos: usize) -> Option<(i64, i64, usize)> {
    if pos >= bytes.len() || !bytes[pos].is_ascii_digit() {
        return None;
    }
    let mut j = pos;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    let base: i64 = std::str::from_utf8(&bytes[pos..j]).ok()?.parse().ok()?;
    // HGVS positions are 1-based; reject `0` (invalid HGVS) so the
    // W3016 scan doesn't compute a `0+N..0+M` range against an
    // explicit ref seq when the input is malformed.
    if base == 0 {
        return None;
    }
    let mut offset: i64 = 0;
    if j < bytes.len() && (bytes[j] == b'+' || bytes[j] == b'-') {
        let sign: i64 = if bytes[j] == b'+' { 1 } else { -1 };
        let off_start = j + 1;
        let mut o = off_start;
        while o < bytes.len() && bytes[o].is_ascii_digit() {
            o += 1;
        }
        if o == off_start {
            // `+` or `-` not followed by digits — not a valid offset
            // shape; surface no detection rather than guess.
            return None;
        }
        let mag: i64 = std::str::from_utf8(&bytes[off_start..o])
            .ok()?
            .parse()
            .ok()?;
        offset = sign * mag;
        j = o;
    }
    Some((base, offset, j))
}

/// Provider-aware resolver: convert a `c.` endpoint (`base`, `offset`)
/// to a genomic position via the transcript's c→g mapper. Returns
/// `None` when the transcript is missing genomic coords, the offset
/// can't be located in any intron, or `cds_to_genomic_with_intron`
/// otherwise fails. Used by the mixed-shape W3016 path (#429).
fn resolve_cds_endpoint_to_genomic(
    provider: &dyn ReferenceProvider,
    accession: &str,
    base: i64,
    offset: i64,
) -> Option<u64> {
    // The backward accession scan stops at `(` / `[` / `;` separators but
    // does not strip a trailing `)` left by a `NG_...(NM_...):c.` selector,
    // so the captured span can be `NM_004006.2)`. Trim enclosing
    // parentheses / whitespace before the provider lookup so the
    // transcript id resolves cleanly (mixed-shape W3016, #429).
    let accession = accession.trim_matches(|c| matches!(c, '(' | ')' | ' ' | '\t' | '\n'));
    let transcript = provider.get_transcript(accession).ok()?;
    let mapper = CoordinateMapper::new(&transcript);
    let cds_pos = CdsPos {
        base,
        offset: if offset == 0 { None } else { Some(offset) },
        utr3: false,
    };
    mapper.cds_to_genomic_with_intron(&cds_pos).ok()
}

/// Consume an IUPAC nucleotide run at `pos`. Returns the substring; empty
/// if no IUPAC base follows. Using IUPAC (rather than any alphabetic byte)
/// ensures the scan stops at the `ins` boundary when reading the ref-seq
/// of a legacy `del<ref>ins<alt>` form (the lowercase `i` is not a valid
/// IUPAC base, so consumption halts there).
fn take_ref_seq_run(bytes: &[u8], pos: usize) -> &str {
    let mut j = pos;
    while j < bytes.len() && is_iupac_base(bytes[j] as char) {
        j += 1;
    }
    std::str::from_utf8(&bytes[pos..j]).unwrap_or("")
}

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

/// Canonicalize thymine (`t`/`T`) to `u` inside `r.` (RNA) descriptions
/// (W3020, closes #282).
///
/// Per HGVS v21.0 RNA nomenclature, the RNA alphabet is `a/c/g/u`; `t` is
/// non-canonical input. PR #293 (issue #276) already canonicalizes `t/T` to
/// `u` on `Display`. This function closes the loop by canonicalizing at
/// preprocess time so the parsed `Base` enum holds `U` from the first parse,
/// and so a soft-validation warning fires per occurrence.
///
/// Scoping rules:
/// - The replacement is confined to the byte span of each `r.` description
///   (start = byte after `r.`, where the `r.` is preceded by `:`, `(`,
///   `[`, `;`, `|`, or appears at start of input; end = either the next
///   observed coord-prefix (`c.`/`g.`/`n.`/`m.`/`o.`/`p.`) under the same
///   prev-byte rule, or a top-level (bracket-depth 0) `;`/`|` separator,
///   or end of input). Bracket/paren bytes themselves do not end the span
///   so a bracketed compound like `r.[1a>t;2t>g]` keeps both halves under
///   the r. context. Top-level `;`/`|` DO end the span so that a malformed
///   unbracketed multi-description input such as
///   `r.1a>t;NT_000001.1:c.5A>T` doesn't rewrite the `T` in the accession
///   `NT_` before the `c.` boundary is reached.
/// - HGVS edit-type keywords inside an `r.` description (`del`, `ins`,
///   `dup`, `inv`, `delins`, `con`, `spl`) never contain `t` or `T`, so any
///   `t`/`T` byte inside an `r.` description span is unambiguously a base.
/// - Each `t`/`T` byte produces one `DetectedCorrection` and is rewritten
///   to lowercase `u` (the canonical RNA letter casing).
pub fn correct_rna_thymine(input: &str) -> (String, Vec<DetectedCorrection>) {
    let mut hits = Vec::new();
    let bytes = input.as_bytes();

    let mut result = String::with_capacity(input.len());
    let mut i = 0usize;
    let mut in_rna = false;
    let mut bracket_depth: u32 = 0;
    while i < bytes.len() {
        let b = bytes[i];

        // Track bracket/paren nesting so we can distinguish a top-level
        // `;`/`|` separator (which ends a description) from a `;`/`|`
        // inside a bracketed compound allele (which does not).
        match b {
            b'[' | b'(' => bracket_depth = bracket_depth.saturating_add(1),
            b']' | b')' => bracket_depth = bracket_depth.saturating_sub(1),
            b';' | b'|' if bracket_depth == 0 => in_rna = false,
            _ => {}
        }

        // Detect entry into a new coordinate-type description. An `r.`
        // start (preceded by `:` / `(` / `[` / `;` / `|`, or at index 0)
        // begins an RNA region; any other coord prefix ends it.
        if i + 1 < bytes.len() && bytes[i + 1] == b'.' {
            let coord = b;
            let prev_ok = i == 0 || matches!(bytes[i - 1], b':' | b'(' | b'[' | b';' | b'|');
            if prev_ok {
                match coord {
                    b'r' => in_rna = true,
                    b'c' | b'g' | b'n' | b'm' | b'o' | b'p' => in_rna = false,
                    _ => {}
                }
            }
        }

        if in_rna && (b == b't' || b == b'T') {
            hits.push(DetectedCorrection::new(
                ErrorType::RnaThymineCanonicalized,
                std::str::from_utf8(&bytes[i..i + 1]).unwrap_or("t"),
                "u",
                i,
                i + 1,
            ));
            result.push('u');
            i += 1;
            continue;
        }
        let ch_end = next_char_end(bytes, i);
        result.push_str(&input[i..ch_end]);
        i = ch_end;
    }
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

    /// Issue #269: `p.0` is the spec-defined "no protein product" form
    /// (HGVS `recommendations/protein/description.md`) — translation is
    /// abolished. The detector must NOT flag it as an invalid position.
    /// The full `p.0?` ("possibly no product") form is also spec-valid.
    #[test]
    fn test_detect_position_zero_skips_protein_axis() {
        assert!(
            detect_position_zero("NP_003997.2:p.0").is_none(),
            "p.0 is spec-valid 'no protein product'"
        );
        assert!(
            detect_position_zero("NP_003997.2:p.0?").is_none(),
            "p.0? is the explicit 'possibly no product' form"
        );
        assert!(
            detect_position_zero("p.0").is_none(),
            "bare p.0 (no accession) must also be accepted"
        );
    }

    /// Issue #269: position zero on every numeric axis is still detected.
    #[test]
    fn test_detect_position_zero_covers_all_numeric_axes() {
        for input in ["c.0A>G", "g.0del", "n.0A>G", "r.0a>g", "m.0A>G", "o.0del"] {
            assert!(
                detect_position_zero(input).is_some(),
                "{input} must be detected as position zero"
            );
        }
    }

    /// CodeRabbit PR #370 review: the axis check used to match any `.0`
    /// whose preceding byte happened to be `c/g/n/r/m/o`, so identifier-
    /// embedded suffixes like `abc.0` were misclassified as position-zero
    /// errors. The detector must require an HGVS boundary before the axis
    /// letter (either string start, or a token separator that the HGVS
    /// grammar can produce between the axis and the previous token).
    #[test]
    fn test_detect_position_zero_requires_axis_boundary() {
        // Identifier-embedded axis byte must NOT trigger.
        for input in ["abc.0", "xyzc.0del", "_g.0", "tagn.0A>G", "FOO123c.0"] {
            assert!(
                detect_position_zero(input).is_none(),
                "{input} has no HGVS axis boundary; must not be flagged"
            );
        }
        // Real HGVS boundaries (accession-separator colon, group openers)
        // before the axis letter must still trigger.
        for input in [
            "NM_000088.3:c.0A>G",
            "g.0del",   // string start
            "(c.0A>G)", // paren-wrapped allele
            "[c.0A>G]", // bracketed allele
            "c.1A>G;c.0del",
        ] {
            assert!(
                detect_position_zero(input).is_some(),
                "{input} has a real HGVS axis boundary; must be flagged"
            );
        }
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
    fn test_detect_swapped_positions_picks_leftmost_marker() {
        // The leftmost coordinate marker is a circular axis (`m.`), so a reversed
        // range there is a valid origin-crossing wraparound and must NOT be flagged.
        // A later embedded `c.` HGVS (e.g. a source expression) with a swapped
        // range must not cause the detector to lock onto that later marker and
        // emit a spurious correction. See issue #467 / PR #474.
        // The dotted accession-with-version form (`...920.1.m.`) exercises the
        // primary marker scan. The leftmost marker is `.m.`; a later `.c.` is
        // earlier in the marker array, so a naive "first marker kind" scan would
        // lock onto the embedded `.c.` range (`500_100`, reversed) and wrongly
        // report a swap. Picking the leftmost marker selects the circular `m.`
        // axis instead, which is skipped (valid origin-crossing wraparound).
        let result =
            detect_swapped_positions("NC_012920.1.m.16100_16000del NM_000088.3.c.500_100del");
        assert!(
            result.is_none(),
            "leftmost marker is circular (m.), so no swap should be reported even though a \
             later embedded c. range is reversed; got {result:?}"
        );
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

    // ----- Issue #115: deprecated multi-base substitution syntax (W3003) -----

    #[test]
    fn test_correct_old_substitution_with_refs() {
        let (out, corrections) = correct_old_substitution_syntax("NM_000088.3:c.79_80GC>TT");
        assert_eq!(out, "NM_000088.3:c.79_80delinsTT");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::OldSubstitutionSyntax);
        assert_eq!(corrections[0].original, "79_80GC>TT");
        assert_eq!(corrections[0].corrected, "79_80delinsTT");
    }

    #[test]
    fn test_correct_old_substitution_no_refs() {
        let (out, corrections) = correct_old_substitution_syntax("NM_000088.3:c.100_102>ATG");
        assert_eq!(out, "NM_000088.3:c.100_102delinsATG");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_old_substitution_single_pos_multi_ref() {
        let (out, corrections) = correct_old_substitution_syntax("NM_000088.3:c.79GC>TT");
        assert_eq!(out, "NM_000088.3:c.79_80delinsTT");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_old_substitution_canonical_no_op() {
        let (out, corrections) = correct_old_substitution_syntax("NM_000088.3:c.100A>G");
        assert_eq!(out, "NM_000088.3:c.100A>G");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_old_substitution_idempotent() {
        let (out1, _) = correct_old_substitution_syntax("c.79_80GC>TT");
        let (out2, c2) = correct_old_substitution_syntax(&out1);
        assert_eq!(out1, out2);
        assert!(c2.is_empty());
    }

    #[test]
    fn test_correct_old_substitution_negative_position() {
        let (out, corrections) = correct_old_substitution_syntax("c.-10_-8>ATG");
        assert_eq!(out, "c.-10_-8delinsATG");
        assert_eq!(corrections.len(), 1);
    }

    // ----- Issue #170: deprecated multi-base substitution on r. (lowercase IUPAC) -----

    #[test]
    fn test_correct_old_substitution_r_lowercase_with_refs() {
        // RNA is canonically lowercase per HGVS; the detector must recognise
        // lowercase IUPAC bases and preserve case when rewriting.
        let (out, corrections) = correct_old_substitution_syntax("NR_001234.1:r.79_80gc>uu");
        assert_eq!(out, "NR_001234.1:r.79_80delinsuu");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::OldSubstitutionSyntax);
        assert_eq!(corrections[0].original, "79_80gc>uu");
        assert_eq!(corrections[0].corrected, "79_80delinsuu");
    }

    #[test]
    fn test_correct_old_substitution_r_lowercase_no_refs() {
        let (out, corrections) = correct_old_substitution_syntax("NR_001234.1:r.100_102>aug");
        assert_eq!(out, "NR_001234.1:r.100_102delinsaug");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_old_substitution_r_lowercase_single_pos_multi_ref() {
        let (out, corrections) = correct_old_substitution_syntax("NR_001234.1:r.79gc>uu");
        assert_eq!(out, "NR_001234.1:r.79_80delinsuu");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_old_substitution_r_lowercase_canonical_no_op() {
        let (out, corrections) = correct_old_substitution_syntax("NR_001234.1:r.100a>g");
        assert_eq!(out, "NR_001234.1:r.100a>g");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_old_substitution_r_lowercase_idempotent() {
        let (out1, _) = correct_old_substitution_syntax("r.79_80gc>uu");
        let (out2, c2) = correct_old_substitution_syntax(&out1);
        assert_eq!(out1, out2);
        assert!(c2.is_empty());
    }

    #[test]
    fn test_correct_old_substitution_mixed_case_preserves_input_case() {
        // IUPAC matching is case-blind but the rewrite must echo the
        // submitter's original case verbatim — no implicit upcase/downcase.
        let (out, corrections) = correct_old_substitution_syntax("r.79_80Gc>uU");
        assert_eq!(out, "r.79_80delinsuU");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].original, "79_80Gc>uU");
        assert_eq!(corrections[0].corrected, "79_80delinsuU");
    }

    // ----- Issue #170: latent W3016 LengthMismatch coverage on lowercase r. -----
    //
    // The same case-blindness of `is_iupac_base` unblocks
    // `take_ref_seq_run`, which `detect_length_mismatch` (W3016) calls to
    // scan the ref-seq of `del<ref>` / `dup<ref>` / `inv<ref>`. Before
    // this PR a lowercase `r.10_15dela` consumed zero ref bytes and the
    // detector returned None; now the lowercase `a` is consumed and the
    // length mismatch is reported as expected.

    #[test]
    fn test_detect_length_mismatch_r_lowercase_del() {
        let hits = detect_length_mismatch("NR_001234.1:r.10_15dela");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].error_type, ErrorType::LengthMismatch);
        assert_eq!(hits[0].original, "10_15dela");
    }

    #[test]
    fn test_detect_length_mismatch_r_lowercase_dup() {
        let hits = detect_length_mismatch("NR_001234.1:r.10_15dupacg");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].error_type, ErrorType::LengthMismatch);
    }

    /// Compound allele inner members inherit the outer coord-system
    /// axis; the W3016 scan must walk them too. Pre-#390 only ranges
    /// preceded by a literal coord marker (`c.`, `g.`, …) — or by `:`
    /// at the variant-string start — triggered detection, so inner
    /// members after `;`/`,` inside `c.[…]` were silently skipped.
    #[test]
    fn test_detect_length_mismatch_compound_bracket_inner_member() {
        // `200_205delAAAA`: 4-base ref vs 6-base range → mismatch.
        let hits = detect_length_mismatch("NM_x:c.[100A>G;200_205delAAAA]");
        assert_eq!(
            hits.len(),
            1,
            "inner member after ';' inside c.[…] should be scanned, got: {:?}",
            hits.iter().map(|h| &h.original).collect::<Vec<_>>()
        );
        assert_eq!(hits[0].error_type, ErrorType::LengthMismatch);
        assert_eq!(hits[0].original, "200_205delAAAA");
    }

    /// First member of a compound allele (right after `[`): no `c.`
    /// prefix immediately precedes it, but it must still be scanned.
    #[test]
    fn test_detect_length_mismatch_compound_bracket_first_member() {
        // `100_103delAA`: 2-base ref vs 4-base range → mismatch.
        let hits = detect_length_mismatch("NM_x:c.[100_103delAA;200A>G]");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].original, "100_103delAA");
    }

    /// Intronic / offset-bearing endpoints with explicit ref
    /// sequences must be length-checked (#390 item 5). Pre-fix the
    /// W3016 endpoint parser bailed on any endpoint followed by
    /// `+` / `-`, so `c.100+5_100+10delAAAA` (intronic span of
    /// 6 bases, 4-base ref) silently passed.
    #[test]
    fn test_detect_length_mismatch_intronic_same_base_positive_offsets() {
        let hits = detect_length_mismatch("NM_x:c.100+5_100+10delAAAA");
        assert_eq!(hits.len(), 1, "got: {:?}", hits);
        assert_eq!(hits[0].original, "100+5_100+10delAAAA");
    }

    /// Same-base, negative-offset endpoints (`100-5_100-2`) — also
    /// linear within the intron, span = `(-2) - (-5) + 1 = 4` bases.
    #[test]
    fn test_detect_length_mismatch_intronic_same_base_negative_offsets() {
        // 4-base range, 2-base ref → mismatch.
        let hits = detect_length_mismatch("NM_x:c.100-5_100-2delAA");
        assert_eq!(hits.len(), 1, "got: {:?}", hits);
        assert_eq!(hits[0].original, "100-5_100-2delAA");
    }

    /// Differing bases or mixed-sign offsets need the intron length
    /// to compute the range, which requires a provider. The
    /// length-mismatch scan must skip these (no false positives).
    #[test]
    fn test_detect_length_mismatch_intronic_mixed_endpoints_skipped() {
        // Different anchor bases — span needs the intron length.
        let hits = detect_length_mismatch("NM_x:c.100+5_101-3delAA");
        assert!(
            hits.is_empty(),
            "mixed-anchor intronic ranges need a provider, got: {:?}",
            hits
        );
    }

    /// Multiple mismatches inside one compound allele are all
    /// surfaced.
    #[test]
    fn test_detect_length_mismatch_compound_bracket_multiple() {
        let hits = detect_length_mismatch("NM_x:c.[100_103delA;200_205delAAAA]");
        assert_eq!(hits.len(), 2);
        let originals: Vec<&str> = hits.iter().map(|h| h.original.as_str()).collect();
        assert!(originals.contains(&"100_103delA"));
        assert!(originals.contains(&"200_205delAAAA"));
    }

    // ----- Issue #115: retracted c.IVS notation (W3014) -----

    #[test]
    fn test_detect_deprecated_ivs_basic() {
        let detections = detect_deprecated_ivs("NM_000088.3:c.IVS2+2T>G");
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].error_type, ErrorType::DeprecatedIvsNotation);
        assert_eq!(detections[0].original, "IVS2");
    }

    #[test]
    fn test_detect_deprecated_ivs_minus_offset() {
        let detections = detect_deprecated_ivs("NM_000088.3:c.IVS5-1G>T");
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].original, "IVS5");
    }

    #[test]
    fn test_detect_deprecated_ivs_canonical_no_op() {
        let detections = detect_deprecated_ivs("NM_000088.3:c.88+2T>G");
        assert!(detections.is_empty());
    }

    #[test]
    fn test_detect_deprecated_ivs_does_not_match_in_accession() {
        let detections = detect_deprecated_ivs("NM_IVS_TEST.3:c.100A>G");
        assert!(detections.is_empty());
    }

    #[test]
    fn test_detect_deprecated_ivs_n_prefix() {
        let detections = detect_deprecated_ivs("NR_001234.1:n.IVS3+1G>A");
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].original, "IVS3");
    }

    #[test]
    fn test_detect_deprecated_ivs_r_prefix() {
        let detections = detect_deprecated_ivs("NR_001234.1:r.IVS3+1g>a");
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].error_type, ErrorType::DeprecatedIvsNotation);
        assert_eq!(detections[0].original, "IVS3");
    }

    #[test]
    fn test_detect_deprecated_ivs_r_lowercase() {
        // Lowercase `ivs` on an `r.` input must be detected on parity
        // with the uppercase form. Original-case is preserved in the
        // diagnostic span text.
        let detections = detect_deprecated_ivs("NR_001234.1:r.ivs3+1g>a");
        assert_eq!(detections.len(), 1);
        assert_eq!(detections[0].error_type, ErrorType::DeprecatedIvsNotation);
        assert_eq!(detections[0].original, "ivs3");
    }

    // ----- Issue #115: deprecated `con` syntax (W3015) -----

    #[test]
    fn test_correct_deprecated_con_basic() {
        let (out, corrections) = correct_deprecated_con("NM_004006.2:c.100_200conNM_001.1:c.5_105");
        assert_eq!(out, "NM_004006.2:c.100_200delinsNM_001.1:c.5_105");
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].error_type, ErrorType::DeprecatedConSyntax);
    }

    #[test]
    fn test_correct_deprecated_con_canonical_no_op() {
        let (out, corrections) =
            correct_deprecated_con("NM_004006.2:c.100_200delinsNM_001.1:c.5_105");
        assert_eq!(out, "NM_004006.2:c.100_200delinsNM_001.1:c.5_105");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_deprecated_con_genomic() {
        let (out, corrections) = correct_deprecated_con("g.1000_2000conNC_000022.10:g.5_1005");
        assert_eq!(out, "g.1000_2000delinsNC_000022.10:g.5_1005");
        assert_eq!(corrections.len(), 1);
    }

    #[test]
    fn test_correct_deprecated_con_idempotent() {
        let (out1, _) = correct_deprecated_con("c.100_200conNM_001.1:c.5_105");
        let (out2, c2) = correct_deprecated_con(&out1);
        assert_eq!(out1, out2);
        assert!(c2.is_empty());
    }

    #[test]
    fn test_correct_deprecated_con_does_not_match_inside_word() {
        let (out, corrections) = correct_deprecated_con("concept text");
        assert_eq!(out, "concept text");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_deprecated_con_does_not_match_no_separator_word() {
        // `con` followed by lowercase letters with no whitespace must not be
        // misread as `delins`-style source. The function only rewrites when
        // the byte after `con` starts a valid HGVS source: digit, `-`, `*`,
        // `(`, or an accession-prefix letter / coordinate-prefix `g.`/`c.`/
        // `n.`/`r.`/`m.`.
        let (out, corrections) = correct_deprecated_con("c.100_200conditional");
        assert_eq!(out, "c.100_200conditional");
        assert!(corrections.is_empty());
    }

    #[test]
    fn test_correct_deprecated_con_bare_position_source() {
        // Source starts with a digit (no accession prefix) — still valid.
        let (out, corrections) = correct_deprecated_con("c.100_200con5_105");
        assert_eq!(out, "c.100_200delins5_105");
        assert_eq!(corrections.len(), 1);
    }

    // -------------------------------------------------------------------------
    // correct_rna_thymine — RNA → non-RNA boundary
    // -------------------------------------------------------------------------

    #[test]
    fn test_correct_rna_thymine_rewrites_only_inside_rna_span() {
        // A compound input that flips from `r.` to `c.` mid-string: the `t`
        // in the `r.` half must be canonicalized to `u`; the `T` in the `c.`
        // half (canonical DNA) must be left alone. This pins the
        // span-handoff documented on `correct_rna_thymine` — the boundary
        // is `c.` (a new coord prefix), not a bracket/paren/semicolon byte.
        let input = "NM_000088.3:r.1a>t;NM_000088.3:c.5A>T";
        let (out, hits) = correct_rna_thymine(input);
        assert_eq!(out, "NM_000088.3:r.1a>u;NM_000088.3:c.5A>T");
        assert_eq!(hits.len(), 1, "exactly one W3020 hit on r. half");
        assert_eq!(hits[0].error_type, ErrorType::RnaThymineCanonicalized);
        assert_eq!(hits[0].original, "t");
        assert_eq!(hits[0].corrected, "u");
    }

    #[test]
    fn test_correct_rna_thymine_at_start_of_input() {
        // Verify the `at start of input` doc claim: an `r.` at byte 0 (no
        // preceding `:` / `(` / `[` / `;` / `|`) still opens an RNA span.
        let (out, hits) = correct_rna_thymine("r.123a>t");
        assert_eq!(out, "r.123a>u");
        assert_eq!(hits.len(), 1);
    }

    #[test]
    fn test_correct_rna_thymine_leaves_non_rna_alone() {
        let (out, hits) = correct_rna_thymine("NM_000088.3:c.123A>T");
        assert_eq!(out, "NM_000088.3:c.123A>T");
        assert!(hits.is_empty());
    }

    #[test]
    fn test_correct_rna_thymine_top_level_semicolon_ends_span() {
        // Top-level `;` ends the r. description span so a `T` in a
        // following accession (e.g. `NT_000001.1`) is not rewritten.
        // Pins CodeRabbit feedback: in_rna previously leaked across a
        // top-level `;`, corrupting `NT_` to `Nu_`.
        let input = "r.1a>t;NT_000001.1:c.5A>T";
        let (out, hits) = correct_rna_thymine(input);
        assert_eq!(out, "r.1a>u;NT_000001.1:c.5A>T");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].error_type, ErrorType::RnaThymineCanonicalized);
        assert_eq!(hits[0].original, "t");
    }

    #[test]
    fn test_correct_rna_thymine_top_level_pipe_ends_span() {
        // Top-level `|` (allele-group separator) behaves like `;`.
        let input = "r.1a>t|NT_000001.1:c.5A>T";
        let (out, _hits) = correct_rna_thymine(input);
        assert_eq!(out, "r.1a>u|NT_000001.1:c.5A>T");
    }

    #[test]
    fn test_correct_rna_thymine_bracketed_compound_keeps_span() {
        // Inside a `[...]` compound allele the `;` does NOT end the r.
        // span — both halves remain under the r. coord and `t`/`T` bytes
        // in either half get rewritten.
        let input = "r.[1a>t;2t>g]";
        let (out, hits) = correct_rna_thymine(input);
        assert_eq!(out, "r.[1a>u;2u>g]");
        assert_eq!(hits.len(), 2);
    }
}
