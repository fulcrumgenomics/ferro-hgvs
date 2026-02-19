//! Edit parsing
//!
//! Parses the edit portion of HGVS variants (substitution, deletion, insertion, etc.)

use crate::hgvs::edit::{
    AminoAcidSeq, Base, ExtDirection, InsertedPart, InsertedSequence, MethylationStatus, NaEdit,
    ProteinEdit, RepeatCount, RepeatUnit, Sequence,
};
use crate::hgvs::location::AminoAcid;
use crate::hgvs::parser::position::{parse_amino_acid, parse_amino_acid_one_letter};
use nom::{
    branch::alt,
    bytes::complete::{tag, take_while1},
    character::complete::{char, digit1},
    combinator::{map, opt},
    multi::many1,
    sequence::preceded,
    IResult, Parser,
};
use std::borrow::Cow;
use std::str::FromStr;

/// Lookup table for valid IUPAC nucleotide characters (O(1) instead of O(n) linear search)
/// Valid: A, C, G, T, U, N, R, Y, S, W, K, M, B, D, H, V (and lowercase)
const fn is_iupac_base(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'N'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'a'
            | b'c'
            | b'g'
            | b't'
            | b'u'
            | b'n'
            | b'r'
            | b'y'
            | b's'
            | b'w'
            | b'k'
            | b'm'
            | b'b'
            | b'd'
            | b'h'
            | b'v'
    )
}

/// Parse a single nucleotide base (including IUPAC ambiguity codes)
#[inline]
fn parse_base(input: &str) -> IResult<&str, Base> {
    let bytes = input.as_bytes();
    if bytes.is_empty() || !is_iupac_base(bytes[0]) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::OneOf,
        )));
    }
    let c = bytes[0] as char;
    Ok((&input[1..], Base::from_char(c).unwrap()))
}

/// Parse a nucleotide sequence (including IUPAC ambiguity codes)
#[inline]
fn parse_sequence(input: &str) -> IResult<&str, Sequence> {
    // Use byte-based lookup table instead of string contains (O(1) vs O(n))
    let bytes = input.as_bytes();
    let mut end = 0;
    while end < bytes.len() && is_iupac_base(bytes[end]) {
        end += 1;
    }

    if end == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::TakeWhile1,
        )));
    }

    let s = &input[..end];
    // Only allocate if not already uppercase (optimization)
    // Use to_ascii_uppercase which is faster for ASCII-only data
    let normalized: Cow<'_, str> = if s.bytes().all(|b| b.is_ascii_uppercase()) {
        Cow::Borrowed(s)
    } else {
        Cow::Owned(s.to_ascii_uppercase())
    };
    Ok((&input[end..], Sequence::from_str(&normalized).unwrap()))
}

/// Parse optional sequence (for deletions that may or may not specify the deleted bases)
fn parse_opt_sequence(input: &str) -> IResult<&str, Option<Sequence>> {
    opt(parse_sequence).parse(input)
}

/// Parse a substitution (e.g., A>G)
/// Only matches single-base substitutions - multi-base patterns like G>AA are handled by parse_multibase_substitution
/// Optimized with direct byte matching to avoid combinator overhead
#[inline]
fn parse_substitution(input: &str) -> IResult<&str, NaEdit> {
    let bytes = input.as_bytes();
    // Need at least 3 bytes: base + '>' + base
    if bytes.len() < 3 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Check pattern: IUPAC_BASE > IUPAC_BASE
    if !is_iupac_base(bytes[0]) || bytes[1] != b'>' || !is_iupac_base(bytes[2]) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Reject if followed by another IUPAC base (this is a multi-base substitution like G>AA)
    if bytes.len() > 3 && is_iupac_base(bytes[3]) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    // Safe: we verified these are valid IUPAC bases
    let reference = Base::from_char(bytes[0] as char).unwrap();
    let alternative = Base::from_char(bytes[2] as char).unwrap();

    Ok((
        &input[3..],
        NaEdit::Substitution {
            reference,
            alternative,
        },
    ))
}

/// Parse a substitution without reference base (e.g., >A, >G)
/// This non-standard notation specifies only the alternative allele
#[inline]
fn parse_substitution_no_ref(input: &str) -> IResult<&str, NaEdit> {
    let bytes = input.as_bytes();
    // Need at least 2 bytes: '>' + base
    if bytes.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Check pattern: '>' + IUPAC_BASE
    if bytes[0] != b'>' || !is_iupac_base(bytes[1]) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Reject if followed by another IUPAC base (this would be >AA which is unusual)
    if bytes.len() > 2 && is_iupac_base(bytes[2]) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    let alternative = Base::from_char(bytes[1] as char).unwrap();

    Ok((&input[2..], NaEdit::SubstitutionNoRef { alternative }))
}

/// Parse a multi-base substitution (e.g., GG>G, ATT>T, G>AA) which is converted to delins
/// This is non-standard HGVS notation but appears in some databases.
/// Handles both multi-ref (GG>G) and single-ref multi-alt (G>AA) patterns.
fn parse_multibase_substitution(input: &str) -> IResult<&str, NaEdit> {
    let (input, ref_seq) = parse_sequence(input)?;
    let (input, _) = tag(">").parse(input)?;
    let (input, alt_seq) = parse_sequence(input)?;

    // Match if either reference or alternative has more than 1 base
    // This catches both GG>G and G>AA patterns
    if ref_seq.len() <= 1 && alt_seq.len() <= 1 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    // Convert to delins
    Ok((
        input,
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(alt_seq),
        },
    ))
}

/// Parse a deletion (e.g., del, delA, delATG, del101)
/// Optimized with byte-based dispatch
#[inline]
fn parse_deletion(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("del").parse(input)?;

    let bytes = input.as_bytes();
    if bytes.is_empty() {
        // Plain "del" with no sequence or length
        return Ok((
            input,
            NaEdit::Deletion {
                sequence: None,
                length: None,
            },
        ));
    }

    match bytes[0] {
        b'0'..=b'9' => {
            // Length: del101
            let (remaining, len_str) = digit1::<&str, nom::error::Error<&str>>(input)?;
            Ok((
                remaining,
                NaEdit::Deletion {
                    sequence: None,
                    length: Some(len_str.parse().unwrap_or(0)),
                },
            ))
        }
        c if is_iupac_base(c) => {
            // Sequence: delATG
            let (remaining, seq) = parse_sequence(input)?;
            Ok((
                remaining,
                NaEdit::Deletion {
                    sequence: Some(seq),
                    length: None,
                },
            ))
        }
        _ => {
            // Neither sequence nor length (e.g., followed by other chars)
            Ok((
                input,
                NaEdit::Deletion {
                    sequence: None,
                    length: None,
                },
            ))
        }
    }
}

/// Parse an insertion (e.g., insA, insATG, ins10, ins(10), ins(10_20), insA[10], ins[A[10];T])
fn parse_insertion(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("ins").parse(input)?;
    let (input, sequence) = parse_inserted_sequence(input)?;
    Ok((input, NaEdit::Insertion { sequence }))
}

/// Parse an inserted sequence in any format
/// Optimized with byte-based dispatch to avoid trying multiple alternatives
#[inline]
fn parse_inserted_sequence(input: &str) -> IResult<&str, InsertedSequence> {
    let bytes = input.as_bytes();

    if bytes.is_empty() {
        // Empty sequence
        return Ok((input, InsertedSequence::Empty));
    }

    // Check for reference accession (e.g., NC_012920.1:m.12435_12527, AF118569:g.14094_14382)
    // This must come before IUPAC base detection since NC, NG, etc. start with IUPAC letters
    // Also handles GenBank accessions like AF118569, AC010542.7
    if let Some(reference) = parse_reference_location(input) {
        let ref_len = reference.len();
        return Ok((&input[ref_len..], InsertedSequence::Reference(reference)));
    }

    match bytes[0] {
        b'[' => {
            // Complex: [A[10];T] or [(10_20)]
            parse_bracketed_inserted_sequence(input)
        }
        b'(' => {
            // Parenthesized range: (10_20)
            parse_parenthesized_count(input)
        }
        b'0'..=b'9' => {
            // Could be CDS position range (423-1933_423-1687inv) or simple count (10)
            // Try CDS position range first if it looks like one (has _ with possible +/- offsets)
            if let Ok((remaining, part)) = parse_cds_position_range(input) {
                // Convert InsertedPart to InsertedSequence
                use crate::hgvs::edit::InsertedPart;
                let seq = match part {
                    InsertedPart::PositionRange { start, end } => {
                        InsertedSequence::Complex(vec![InsertedPart::PositionRange { start, end }])
                    }
                    InsertedPart::PositionRangeInv { start, end } => {
                        InsertedSequence::Complex(vec![InsertedPart::PositionRangeInv {
                            start,
                            end,
                        }])
                    }
                    InsertedPart::CdsPositionRange(s) => {
                        InsertedSequence::Complex(vec![InsertedPart::CdsPositionRange(s)])
                    }
                    _ => InsertedSequence::Empty,
                };
                Ok((remaining, seq))
            } else {
                // Simple count: 10
                parse_simple_count(input)
            }
        }
        c if is_iupac_base(c) => {
            // Could be repeated base (A[10]), literal sequence (ATG), or named element (AluYb8, ALU, LINE1)
            // Check if followed by '[' for repeat notation
            if bytes.len() > 1 && bytes[1] == b'[' {
                parse_repeated_base_insertion(input)
            } else {
                // Check if this is a pure IUPAC sequence or a named element (like AluYb8, ALU, LINE1)
                // Named elements can contain:
                // - lowercase letters (IUPAC is uppercase only)
                // - digits mixed with letters
                // - uppercase letters that are NOT valid IUPAC bases (e.g., L in ALU)
                let mut has_non_iupac = false;
                let mut i = 0;
                while i < bytes.len() {
                    let b = bytes[i];
                    if is_iupac_base(b) {
                        i += 1;
                    } else if b.is_ascii_lowercase() {
                        // Lowercase letter - this is a named element (IUPAC uses uppercase)
                        has_non_iupac = true;
                        i += 1;
                    } else if b.is_ascii_digit() && i > 0 {
                        // Digit after letters - could be named element
                        has_non_iupac = true;
                        i += 1;
                    } else if b.is_ascii_uppercase() {
                        // Uppercase letter that's NOT a valid IUPAC base (e.g., L in ALU)
                        // This indicates a named element
                        has_non_iupac = true;
                        i += 1;
                    } else {
                        break;
                    }
                }

                if has_non_iupac && i > 0 {
                    // Named element like AluYb8, ALU, LINE1, L1
                    let name = &input[..i];
                    return Ok((&input[i..], InsertedSequence::Named(name.to_string())));
                }

                // Literal IUPAC sequence - but check for repeat count suffix
                let (remaining, seq) = parse_sequence(input)?;
                // Check if followed by [count] for sequence repeat
                if remaining.starts_with('[') {
                    if let Ok((remaining2, count)) = parse_repeat_count(remaining) {
                        return Ok((
                            remaining2,
                            InsertedSequence::SequenceRepeat {
                                sequence: seq,
                                count,
                            },
                        ));
                    }
                }
                Ok((remaining, InsertedSequence::Literal(seq)))
            }
        }
        c if c.is_ascii_uppercase() => {
            // Non-IUPAC uppercase letter - could be a named element (e.g., LINE1, L1)
            // Collect all alphanumeric characters
            let mut i = 0;
            while i < bytes.len() {
                let b = bytes[i];
                if b.is_ascii_alphanumeric() {
                    i += 1;
                } else {
                    break;
                }
            }
            if i > 0 {
                let name = &input[..i];
                Ok((&input[i..], InsertedSequence::Named(name.to_string())))
            } else {
                Ok((input, InsertedSequence::Empty))
            }
        }
        _ => {
            // Empty (unrecognized character means end of insertion)
            Ok((input, InsertedSequence::Empty))
        }
    }
}

/// Parse a bracketed inserted sequence: [A[10];T] or [(10_20)] or [NC_000022.11:g.100_200]
fn parse_bracketed_inserted_sequence(input: &str) -> IResult<&str, InsertedSequence> {
    let (input, _) = char('[').parse(input)?;

    // Check for reference location inside brackets (e.g., [NC_000022.11:g.100_200])
    // Reference accession prefixes: NC_, NG_, NM_, NP_, NR_, NT_, NW_, XM_, XP_, XR_, ENST, ENSP, LRG_
    if is_reference_accession_prefix(input) {
        // Find the closing bracket
        if let Some(close_pos) = input.find(']') {
            let reference = &input[..close_pos];
            let remaining = &input[close_pos + 1..];
            return Ok((
                remaining,
                InsertedSequence::Reference(reference.to_string()),
            ));
        }
    }

    // Check for (min_max) range inside brackets
    if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>('(').parse(input) {
        let (remaining, min) = digit1.parse(remaining)?;
        let (remaining, _) = char('_').parse(remaining)?;
        let (remaining, max) = digit1.parse(remaining)?;
        let (remaining, _) = char(')').parse(remaining)?;
        let (remaining, _) = char(']').parse(remaining)?;
        return Ok((
            remaining,
            InsertedSequence::Range(min.parse().unwrap_or(0), max.parse().unwrap_or(0)),
        ));
    }

    // Parse parts separated by semicolon
    // Pre-allocate for typical case (most complex insertions have 1-3 parts)
    let mut parts = Vec::with_capacity(4);
    let mut remaining_input = input;

    while let Ok((remaining, part)) = parse_inserted_part(remaining_input) {
        parts.push(part);
        remaining_input = remaining;

        // Check for semicolon separator
        if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>(';').parse(remaining_input) {
            remaining_input = remaining;
        } else {
            break;
        }
    }

    let (input, _) = char(']').parse(remaining_input)?;

    if parts.is_empty() {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )))
    } else {
        Ok((input, InsertedSequence::Complex(parts)))
    }
}

/// Check if the input starts with a reference accession prefix
///
/// Supports:
/// - RefSeq: NC_, NG_, NM_, NP_, NR_, NT_, NW_, XM_, XP_, XR_
/// - Ensembl: ENST, ENSP, ENSG
/// - LRG: LRG_
/// - GenBank: Two-letter prefixes like PQ, PX, MT, MF, KT, KY, DQ, AC, AB, PP, etc.
fn is_reference_accession_prefix(input: &str) -> bool {
    // Check for known multi-character prefixes first
    let known_prefixes = [
        "NC_", "NG_", "NM_", "NP_", "NR_", "NT_", "NW_", "XM_", "XP_", "XR_", "ENST", "ENSP",
        "ENSG", "LRG_",
    ];
    for prefix in known_prefixes {
        if input.starts_with(prefix) {
            return true;
        }
    }

    // Check for GenBank-style accessions: 1-2 letters followed by digits
    // Pattern: [A-Z]{1,2}[0-9]+ (e.g., PQ998981, MT113356, U12345)
    let bytes = input.as_bytes();
    if bytes.len() >= 3 {
        // Must start with uppercase letter
        if bytes[0].is_ascii_uppercase() {
            // Check for 1 or 2 letter prefix followed by digit
            if bytes[1].is_ascii_digit() {
                // Single letter prefix (e.g., U12345)
                return true;
            } else if bytes[1].is_ascii_uppercase() && bytes.len() >= 4 && bytes[2].is_ascii_digit()
            {
                // Two letter prefix (e.g., PQ998981)
                return true;
            }
        }
    }

    false
}

/// Parse a reference location (e.g., NC_012920.1:m.12435_12527 or AF118569:g.14094_14382)
/// Returns the full reference string if valid, None otherwise
fn parse_reference_location(input: &str) -> Option<String> {
    // Must contain a colon followed by coordinate type
    let colon_pos = input.find(':')?;

    // Check that we have at least "X:Y." pattern
    let after_colon = &input[colon_pos + 1..];
    if after_colon.len() < 2 {
        return None;
    }

    // Must have coordinate type (g., c., m., n., r., p., o.)
    let coord_type = after_colon.as_bytes()[0];
    let dot = after_colon.as_bytes()[1];
    if dot != b'.' || !matches!(coord_type, b'g' | b'c' | b'm' | b'n' | b'r' | b'p' | b'o') {
        return None;
    }

    // The accession part must look like an accession (alphanumeric, dots, underscores)
    let accession = &input[..colon_pos];
    if accession.is_empty() {
        return None;
    }

    // Check that accession starts with a letter and contains valid accession characters
    let acc_bytes = accession.as_bytes();
    if !acc_bytes[0].is_ascii_alphabetic() {
        return None;
    }

    for &b in acc_bytes {
        if !b.is_ascii_alphanumeric() && b != b'_' && b != b'.' {
            return None;
        }
    }

    // Now find where the coordinates end
    let mut end_pos = colon_pos + 3; // Skip past ":X."
    let pos_part = &input[end_pos..];

    for c in pos_part.chars() {
        if c.is_ascii_digit() || c == '_' || c == '+' || c == '-' || c == '?' {
            end_pos += c.len_utf8();
        } else {
            break;
        }
    }

    // Must have at least one position character
    if end_pos <= colon_pos + 3 {
        return None;
    }

    Some(input[..end_pos].to_string())
}

/// Parse an external sequence reference as an inserted part
fn parse_external_ref_part(input: &str) -> IResult<&str, crate::hgvs::edit::InsertedPart> {
    use crate::hgvs::edit::InsertedPart;

    if !is_reference_accession_prefix(input) {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Find the end of the reference (up to ; or ])
    let end_pos = input.find([';', ']']).unwrap_or(input.len());

    if end_pos == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let reference = &input[..end_pos];
    // Validate it contains a colon (accession:position format)
    if !reference.contains(':') {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    Ok((
        &input[end_pos..],
        InsertedPart::ExternalRef(reference.to_string()),
    ))
}

/// Parse a CDS-style position range (e.g., 401_419, 244-8_249, 100+5_200-10, with optional inv suffix)
fn parse_cds_position_range(input: &str) -> IResult<&str, crate::hgvs::edit::InsertedPart> {
    use crate::hgvs::edit::InsertedPart;

    // Try to parse CDS-style position (may have offsets like -8 or +5)
    // Format: base[+-offset]_base[+-offset][inv]
    let mut end_pos = 0;
    let bytes = input.as_bytes();

    // Parse first position (digits with optional +/- offset)
    while end_pos < bytes.len() && bytes[end_pos].is_ascii_digit() {
        end_pos += 1;
    }
    if end_pos == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Digit,
        )));
    }
    // Check for offset (+N or -N)
    if end_pos < bytes.len() && (bytes[end_pos] == b'+' || bytes[end_pos] == b'-') {
        end_pos += 1;
        while end_pos < bytes.len() && bytes[end_pos].is_ascii_digit() {
            end_pos += 1;
        }
    }

    // Must have underscore separator
    if end_pos >= bytes.len() || bytes[end_pos] != b'_' {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }
    end_pos += 1; // skip underscore

    // Parse second position
    let second_start = end_pos;
    while end_pos < bytes.len() && bytes[end_pos].is_ascii_digit() {
        end_pos += 1;
    }
    if end_pos == second_start {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Digit,
        )));
    }
    // Check for offset on second position
    if end_pos < bytes.len() && (bytes[end_pos] == b'+' || bytes[end_pos] == b'-') {
        end_pos += 1;
        while end_pos < bytes.len() && bytes[end_pos].is_ascii_digit() {
            end_pos += 1;
        }
    }

    let range_str = &input[..end_pos];
    let mut remaining = &input[end_pos..];
    let has_offset = range_str.contains('+') || range_str.contains('-');

    // Check for inversion suffix
    let has_inv = if remaining.starts_with("inv") {
        remaining = &remaining[3..];
        true
    } else {
        false
    };

    // Determine which variant to use based on whether there are offsets
    if has_offset {
        let range_with_inv = if has_inv {
            format!("{}inv", range_str)
        } else {
            range_str.to_string()
        };
        Ok((remaining, InsertedPart::CdsPositionRange(range_with_inv)))
    } else {
        // Parse as simple numbers
        let parts: Vec<&str> = range_str.split('_').collect();
        if parts.len() == 2 {
            let start: u64 = parts[0].parse().unwrap_or(0);
            let end: u64 = parts[1].parse().unwrap_or(0);
            if has_inv {
                Ok((remaining, InsertedPart::PositionRangeInv { start, end }))
            } else {
                Ok((remaining, InsertedPart::PositionRange { start, end }))
            }
        } else {
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
    }
}

/// Parse a single part of a complex insertion
fn parse_inserted_part(input: &str) -> IResult<&str, crate::hgvs::edit::InsertedPart> {
    use crate::hgvs::edit::InsertedPart;

    alt((
        // External sequence reference: KT192064.1:1_310 or NC_000020.11:g.2823027_2826302
        parse_external_ref_part,
        // Repeated base: A[10], N[15_30]
        |input| {
            let (remaining, base) = parse_base(input)?;
            let (remaining, count) = parse_repeat_count(remaining)?;
            Ok((remaining, InsertedPart::Repeat { base, count }))
        },
        // Position range reference with optional inversion: 401_419 or 45043709_45310738inv
        // Also handles CDS-style positions with offsets: 244-8_249 or 100+5_200-10
        parse_cds_position_range,
        // Single base or short sequence
        map(parse_sequence, InsertedPart::Literal),
    ))
    .parse(input)
}

/// Parse a parenthesized count: (10) or (10_20)
fn parse_parenthesized_count(input: &str) -> IResult<&str, InsertedSequence> {
    let (input, _) = char('(').parse(input)?;

    // Check for uncertain: (?)
    if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("?)").parse(input) {
        return Ok((remaining, InsertedSequence::Uncertain));
    }

    // Check for uncertain range: (min_?)
    if let Ok((remaining, min_str)) = digit1::<_, nom::error::Error<&str>>.parse(input) {
        if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("_?)").parse(remaining) {
            return Ok((
                remaining,
                InsertedSequence::Range(min_str.parse().unwrap_or(0), u64::MAX),
            ));
        }

        // Check for range: (min_max)
        if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>('_').parse(remaining) {
            let (remaining, max_str) = digit1.parse(remaining)?;
            let (remaining, _) = char(')').parse(remaining)?;

            let start1: u64 = min_str.parse().unwrap_or(0);
            let end1: u64 = max_str.parse().unwrap_or(0);

            // Check if followed by _(pos_pos) for complex position ranges
            // e.g., ins(23632682_23625413)_(23625324_23619334)
            if let Some(rest) = remaining.strip_prefix("_(") {
                if let Ok((rest, start2_str)) = digit1::<_, nom::error::Error<&str>>.parse(rest) {
                    if let Ok((rest, _)) = char::<_, nom::error::Error<&str>>('_').parse(rest) {
                        if let Ok((rest, end2_str)) =
                            digit1::<_, nom::error::Error<&str>>.parse(rest)
                        {
                            if let Ok((rest, _)) =
                                char::<_, nom::error::Error<&str>>(')').parse(rest)
                            {
                                let start2: u64 = start2_str.parse().unwrap_or(0);
                                let end2: u64 = end2_str.parse().unwrap_or(0);
                                return Ok((
                                    rest,
                                    InsertedSequence::Complex(vec![
                                        InsertedPart::PositionRange {
                                            start: start1,
                                            end: end1,
                                        },
                                        InsertedPart::PositionRange {
                                            start: start2,
                                            end: end2,
                                        },
                                    ]),
                                ));
                            }
                        }
                    }
                }
            }

            return Ok((remaining, InsertedSequence::Range(start1, end1)));
        }

        let (remaining, _) = char(')').parse(remaining)?;
        return Ok((
            remaining,
            InsertedSequence::Count(min_str.parse().unwrap_or(0)),
        ));
    }

    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Digit,
    )))
}

/// Parse a repeated base insertion: A[10], N[15], N[15_30]
fn parse_repeated_base_insertion(input: &str) -> IResult<&str, InsertedSequence> {
    let (input, base) = parse_base(input)?;
    let (input, count) = parse_repeat_count(input)?;

    Ok((input, InsertedSequence::Repeat { base, count }))
}

/// Parse a simple count: 10 (just digits)
#[inline]
fn parse_simple_count(input: &str) -> IResult<&str, InsertedSequence> {
    let (remaining, count_str) = digit1.parse(input)?;
    let start: u64 = count_str.parse().unwrap_or(0);

    // Check for position range: start_end or start_endinv
    if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>('_').parse(remaining) {
        let (remaining, end_str) = digit1.parse(remaining)?;
        let end: u64 = end_str.parse().unwrap_or(0);

        // Check for inversion suffix
        if let Ok((remaining, _)) = tag::<_, _, nom::error::Error<&str>>("inv").parse(remaining) {
            return Ok((remaining, InsertedSequence::PositionRangeInv { start, end }));
        }

        return Ok((remaining, InsertedSequence::PositionRange { start, end }));
    }

    // Make sure it's not followed by more bases (which would be a sequence)
    // Use byte access instead of char iteration (optimization)
    if remaining
        .as_bytes()
        .first()
        .is_some_and(|&b| b.is_ascii_alphabetic())
    {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    Ok((remaining, InsertedSequence::Count(start)))
}

/// Parse a delins (e.g., delinsA, delinsATG, delins10, delins(10_20), delinsA[10])
fn parse_delins(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("delins").parse(input)?;
    let (input, sequence) = parse_inserted_sequence(input)?;
    Ok((input, NaEdit::Delins { sequence }))
}

/// Parse a delins with explicit deleted sequence (e.g., delAinsT, delATGinsCCC)
/// This is a legacy/redundant format but appears in some databases.
/// The deleted sequence is parsed and then discarded (normalized to delinsXXX format).
fn parse_delins_with_deleted_seq(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("del").parse(input)?;
    let (input, _deleted_seq) = parse_sequence(input)?;
    let (input, _) = tag("ins").parse(input)?;
    let (input, inserted_seq) = parse_inserted_sequence(input)?;

    Ok((
        input,
        NaEdit::Delins {
            sequence: inserted_seq,
        },
    ))
}

/// Parse a delins with explicit deleted count (e.g., del35insTA, del17insTA)
/// This is a legacy/redundant format but appears in some databases.
/// The deleted count is parsed and then discarded (normalized to delinsXXX format).
fn parse_delins_with_deleted_count(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("del").parse(input)?;
    let (input, _deleted_count) = digit1.parse(input)?;
    let (input, _) = tag("ins").parse(input)?;
    let (input, inserted_seq) = parse_inserted_sequence(input)?;

    Ok((
        input,
        NaEdit::Delins {
            sequence: inserted_seq,
        },
    ))
}

/// Parse a dupins (duplication followed by insertion, e.g., dupinsCTCA)
/// This is a non-standard but commonly used notation in databases like ClinVar
#[inline]
fn parse_dupins(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("dupins").parse(input)?;
    let (input, sequence) = parse_inserted_sequence(input)?;
    Ok((input, NaEdit::DupIns { sequence }))
}

/// Parse a duplication (e.g., dup, dupA, dupATG, dup101, dup(731_741), dup?, dup(?))
/// Optimized with byte-based dispatch
#[inline]
fn parse_duplication(input: &str) -> IResult<&str, NaEdit> {
    use crate::hgvs::edit::UncertainDupExtent;

    let (input, _) = tag("dup").parse(input)?;

    let bytes = input.as_bytes();

    // Fast path: empty input means plain "dup"
    if bytes.is_empty() {
        return Ok((
            input,
            NaEdit::Duplication {
                sequence: None,
                length: None,
                uncertain_extent: None,
            },
        ));
    }

    // Parse sequence or length based on first byte
    let (input, sequence, length) = match bytes[0] {
        c if is_iupac_base(c) => {
            let (remaining, seq) = parse_sequence(input)?;
            (remaining, Some(seq), None)
        }
        b'0'..=b'9' => {
            let (remaining, len_str) = digit1::<&str, nom::error::Error<&str>>(input)?;
            (remaining, None, Some(len_str.parse().unwrap_or(0)))
        }
        _ => (input, None, None),
    };

    // Parse uncertain extent: ?, (?), (start_end)
    let bytes = input.as_bytes();
    let (input, uncertain_extent) = if bytes.is_empty() {
        (input, None)
    } else {
        match bytes[0] {
            b'?' => (&input[1..], Some(UncertainDupExtent::Unknown)),
            b'(' => {
                // dup(?) or dup(start_end)
                let remaining = &input[1..];
                if let Some(stripped) = remaining.strip_prefix("?)") {
                    (stripped, Some(UncertainDupExtent::Unknown))
                } else if let Ok((remaining, start_str)) =
                    digit1::<&str, nom::error::Error<&str>>(remaining)
                {
                    if let Ok((remaining, _)) =
                        char::<_, nom::error::Error<&str>>('_').parse(remaining)
                    {
                        if let Ok((remaining, end_str)) =
                            digit1::<&str, nom::error::Error<&str>>(remaining)
                        {
                            if let Ok((remaining, _)) =
                                char::<_, nom::error::Error<&str>>(')').parse(remaining)
                            {
                                let start = start_str.parse().unwrap_or(0);
                                let end = end_str.parse().unwrap_or(0);
                                (remaining, Some(UncertainDupExtent::Range(start, end)))
                            } else {
                                (input, None)
                            }
                        } else {
                            (input, None)
                        }
                    } else {
                        (input, None)
                    }
                } else {
                    (input, None)
                }
            }
            _ => (input, None),
        }
    };

    Ok((
        input,
        NaEdit::Duplication {
            sequence,
            length,
            uncertain_extent,
        },
    ))
}

/// Parse an inversion (inv, inv3, invATG)
fn parse_inversion(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("inv").parse(input)?;

    // Try to parse sequence first, then length if no sequence
    if let Ok((remaining, seq)) = parse_sequence(input) {
        return Ok((
            remaining,
            NaEdit::Inversion {
                sequence: Some(seq),
                length: None,
            },
        ));
    }

    // Try to parse length (digits only)
    if let Ok((remaining, len_str)) = digit1::<&str, nom::error::Error<&str>>(input) {
        if let Ok(len) = len_str.parse::<u64>() {
            return Ok((
                remaining,
                NaEdit::Inversion {
                    sequence: None,
                    length: Some(len),
                },
            ));
        }
    }

    // Neither sequence nor length
    Ok((
        input,
        NaEdit::Inversion {
            sequence: None,
            length: None,
        },
    ))
}

/// Scan ASCII digits from bytes, returning (parsed_value, bytes_consumed)
#[inline]
fn scan_digits(bytes: &[u8]) -> (u64, usize) {
    let mut end = 0;
    let mut value = 0u64;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        value = value * 10 + (bytes[end] - b'0') as u64;
        end += 1;
    }
    (value, end)
}

/// Parse a repeat count (e.g., [12], [10_15], [?], [?_?], [10_?], [?_10])
/// Optimized with byte-based dispatch to avoid alt() backtracking
#[inline]
fn parse_repeat_count(input: &str) -> IResult<&str, RepeatCount> {
    let bytes = input.as_bytes();

    // Must start with '['
    if bytes.is_empty() || bytes[0] != b'[' {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let inner = &bytes[1..];
    if inner.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Dispatch based on first character after '['
    match inner[0] {
        b'(' => {
            // Parenthesized range: [(7_28)] or [(200_?)] or [(?_20)] - range in parens within brackets
            // Skip the opening paren
            let after_paren = &inner[1..];

            // Check for [(?_num)] pattern
            if !after_paren.is_empty() && after_paren[0] == b'?' {
                // [(?_num)] or [(?_?)]
                if after_paren.len() > 1 && after_paren[1] == b'_' {
                    let after_underscore = &after_paren[2..];
                    if !after_underscore.is_empty() && after_underscore[0] == b'?' {
                        // [(?_?)]
                        if after_underscore.len() > 1
                            && after_underscore[1] == b')'
                            && after_underscore.len() > 2
                            && after_underscore[2] == b']'
                        {
                            return Ok((&input[7..], RepeatCount::Unknown)); // [(?_?)]
                        }
                    } else {
                        // [(?_num)]
                        let (num, consumed) = scan_digits(after_underscore);
                        if consumed > 0 {
                            let check_pos = consumed;
                            if after_underscore.len() > check_pos + 1
                                && after_underscore[check_pos] == b')'
                                && after_underscore[check_pos + 1] == b']'
                            {
                                return Ok((
                                    &input[6 + consumed..],
                                    RepeatCount::MaxUncertain(num),
                                ));
                            }
                        }
                    }
                }
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Tag,
                )));
            }

            let (num1, consumed1) = scan_digits(after_paren);
            if consumed1 == 0 {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Digit,
                )));
            }
            let after_num1 = &after_paren[consumed1..];
            // Expect _
            if after_num1.is_empty() || after_num1[0] != b'_' {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Tag,
                )));
            }

            // Check for [(num_?)] pattern
            let after_underscore = &after_num1[1..];
            if !after_underscore.is_empty() && after_underscore[0] == b'?' {
                // [(num_?)]
                if after_underscore.len() > 1
                    && after_underscore[1] == b')'
                    && after_underscore.len() > 2
                    && after_underscore[2] == b']'
                {
                    // Total consumed: [ + ( + num1 + _ + ? + ) + ] = 1 + 1 + consumed1 + 1 + 1 + 1 + 1
                    let total_consumed = 6 + consumed1;
                    return Ok((&input[total_consumed..], RepeatCount::MinUncertain(num1)));
                }
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Tag,
                )));
            }

            // Parse second number for [(num1_num2)]
            let (num2, consumed2) = scan_digits(after_underscore);
            if consumed2 == 0 {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Digit,
                )));
            }
            let after_num2 = &after_underscore[consumed2..];
            // Expect )]
            if after_num2.len() < 2 || after_num2[0] != b')' || after_num2[1] != b']' {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Tag,
                )));
            }
            // Total consumed: [ + ( + num1 + _ + num2 + ) + ] = 1 + 1 + consumed1 + 1 + consumed2 + 1 + 1
            let total_consumed = 5 + consumed1 + consumed2;
            Ok((&input[total_consumed..], RepeatCount::Range(num1, num2)))
        }
        b'?' => {
            // Could be: [?], [?_?], [?_20]
            if inner.len() > 1 && inner[1] == b'_' {
                if inner.len() > 2 && inner[2] == b'?' {
                    // [?_?] - fully uncertain
                    if inner.len() > 3 && inner[3] == b']' {
                        return Ok((&input[5..], RepeatCount::Unknown));
                    }
                } else {
                    // [?_<digits>] - max uncertain
                    let (value, consumed) = scan_digits(&inner[2..]);
                    if consumed > 0 {
                        let bracket_pos = 2 + consumed;
                        if inner.len() > bracket_pos && inner[bracket_pos] == b']' {
                            return Ok((
                                &input[bracket_pos + 2..],
                                RepeatCount::MaxUncertain(value),
                            ));
                        }
                    }
                }
            } else if inner.len() > 1 && inner[1] == b']' {
                // [?] - unknown
                return Ok((&input[3..], RepeatCount::Unknown));
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        b'0'..=b'9' => {
            // Could be: [12], [10_?], [10_15]
            let (num1, consumed1) = scan_digits(inner);
            if consumed1 == 0 {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Digit,
                )));
            }

            let after_num1 = &inner[consumed1..];
            if after_num1.is_empty() {
                return Err(nom::Err::Error(nom::error::Error::new(
                    input,
                    nom::error::ErrorKind::Tag,
                )));
            }

            if after_num1[0] == b']' {
                // [12] - exact count
                return Ok((&input[consumed1 + 2..], RepeatCount::Exact(num1)));
            } else if after_num1[0] == b'_' {
                // [10_?] or [10_15]
                if after_num1.len() > 1 && after_num1[1] == b'?' {
                    // [10_?] - min uncertain
                    if after_num1.len() > 2 && after_num1[2] == b']' {
                        return Ok((&input[consumed1 + 4..], RepeatCount::MinUncertain(num1)));
                    }
                } else {
                    // [10_15] - range
                    let (num2, consumed2) = scan_digits(&after_num1[1..]);
                    if consumed2 > 0 {
                        let bracket_pos = consumed1 + 1 + consumed2;
                        if inner.len() > bracket_pos && inner[bracket_pos] == b']' {
                            return Ok((&input[bracket_pos + 2..], RepeatCount::Range(num1, num2)));
                        }
                    }
                }
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        _ => Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        ))),
    }
}

/// Parse a repeat (e.g., [12], CAG[12], [10_15], [3]T, A[6][1], GT[2]GC[2])
/// VEP sometimes outputs repeats with trailing nucleotides (e.g., c.212-18[3]T)
/// Genotype notation uses multiple counts (e.g., A[6][1] for heterozygous)
/// Multi-repeat notation has multiple sequence[count] pairs (e.g., GT[2]GC[2]GTGCATGAGTGTGCG[1])
#[inline]
fn parse_repeat(input: &str) -> IResult<&str, NaEdit> {
    let (input, first_sequence) = parse_opt_sequence(input)?;
    let (input, first_count) = parse_repeat_count(input)?;

    // Look ahead to see what comes next
    // - If it's [count] (no sequence), it's genotype notation like A[6][1]
    // - If it's sequence[count], it's multi-repeat like GT[2]GC[2]
    // - If it's just sequence with no [count], it's trailing sequence like [3]T

    let mut remaining = input;

    // First, check if the next character starts a sequence (letter) followed by [count]
    // This distinguishes multi-repeat from genotype notation
    if let Ok((after_seq, next_seq)) = parse_sequence(remaining) {
        // We got a sequence. Check if it's followed by [count]
        if parse_repeat_count(after_seq).is_ok() {
            // This is multi-repeat notation: sequence[count]sequence[count]...
            let mut units = vec![RepeatUnit {
                sequence: first_sequence.unwrap_or_else(|| Sequence::new(vec![])),
                count: first_count,
            }];

            // Parse the next unit we already peeked at
            let (after_count, next_count) = parse_repeat_count(after_seq)?;
            units.push(RepeatUnit {
                sequence: next_seq,
                count: next_count,
            });
            remaining = after_count;

            // Continue parsing more sequence[count] pairs
            while let Ok((after_seq, seq)) = parse_sequence(remaining) {
                if let Ok((after_count, count)) = parse_repeat_count(after_seq) {
                    units.push(RepeatUnit {
                        sequence: seq,
                        count,
                    });
                    remaining = after_count;
                } else {
                    break;
                }
            }

            return Ok((remaining, NaEdit::MultiRepeat { units }));
        }
    }

    // Not multi-repeat - try genotype notation (additional counts without sequence)
    let mut additional_counts = Vec::new();
    while let Ok((next_input, additional)) = parse_repeat_count(remaining) {
        additional_counts.push(additional);
        remaining = next_input;
    }

    // Parse optional trailing sequence
    let (remaining, trailing) = parse_opt_sequence(remaining)?;

    Ok((
        remaining,
        NaEdit::Repeat {
            sequence: first_sequence,
            count: first_count,
            additional_counts,
            trailing,
        },
    ))
}

/// Parse identity (=, =ATG, or G=)
/// HGVS uses G= format for identity with reference base (e.g., c.6317G=)
/// Optimized with byte-based dispatch
#[inline]
fn parse_identity(input: &str) -> IResult<&str, NaEdit> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    if bytes[0] == b'=' {
        // Format: =ATG or just =
        let (remaining, sequence) = parse_opt_sequence(&input[1..])?;
        Ok((
            remaining,
            NaEdit::Identity {
                sequence,
                whole_entity: false,
            },
        ))
    } else if is_iupac_base(bytes[0]) {
        // Format: G= (sequence before equals sign)
        let (remaining, sequence) = parse_sequence(input)?;
        let remaining_bytes = remaining.as_bytes();
        if remaining_bytes.is_empty() || remaining_bytes[0] != b'=' {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )));
        }
        Ok((
            &remaining[1..],
            NaEdit::Identity {
                sequence: Some(sequence),
                whole_entity: false,
            },
        ))
    } else {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )))
    }
}

/// Parse conversion (conNM_000001.1:c.100_200)
fn parse_conversion(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("con").parse(input)?;
    // Take everything until end or whitespace as the source
    let (input, source) = take_while1(|c: char| !c.is_whitespace()).parse(input)?;
    Ok((
        input,
        NaEdit::Conversion {
            source: source.to_string(),
        },
    ))
}

/// Parse methylation status (|gom, |lom, |met=)
fn parse_methylation(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = char('|').parse(input)?;
    alt((
        map(tag("gom"), |_| NaEdit::Methylation {
            status: MethylationStatus::GainOfMethylation,
        }),
        map(tag("lom"), |_| NaEdit::Methylation {
            status: MethylationStatus::LossOfMethylation,
        }),
        map(tag("met="), |_| NaEdit::Methylation {
            status: MethylationStatus::Unchanged,
        }),
    ))
    .parse(input)
}

/// Parse copy number (copy2, copy4, etc.)
fn parse_copy_number(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = tag("copy").parse(input)?;
    let (input, count_str) = digit1.parse(input)?;
    let count = count_str.parse::<u64>().unwrap_or(2);
    Ok((input, NaEdit::CopyNumber { count }))
}

/// Parse a parenthesized repeat count at the end of an interval (e.g., `g.100_200(30_125)`)
/// This represents a repeat with uncertain count.
fn parse_parenthesized_repeat(input: &str) -> IResult<&str, NaEdit> {
    let (input, _) = char('(').parse(input)?;
    let (input, min_str) = digit1.parse(input)?;

    // Check for range: (min_max)
    if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>('_').parse(input) {
        // Check for uncertain max: (min_?)
        if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>('?').parse(remaining) {
            let (remaining, _) = char(')').parse(remaining)?;
            return Ok((
                remaining,
                NaEdit::Repeat {
                    sequence: None,
                    count: RepeatCount::MinUncertain(min_str.parse().unwrap_or(0)),
                    additional_counts: Vec::new(),
                    trailing: None,
                },
            ));
        }

        let (remaining, max_str) = digit1.parse(remaining)?;
        let (remaining, _) = char(')').parse(remaining)?;
        return Ok((
            remaining,
            NaEdit::Repeat {
                sequence: None,
                count: RepeatCount::Range(
                    min_str.parse().unwrap_or(0),
                    max_str.parse().unwrap_or(0),
                ),
                additional_counts: Vec::new(),
                trailing: None,
            },
        ));
    }

    // Just a single count: (min)
    let (input, _) = char(')').parse(input)?;
    Ok((
        input,
        NaEdit::Repeat {
            sequence: None,
            count: RepeatCount::Exact(min_str.parse().unwrap_or(0)),
            additional_counts: Vec::new(),
            trailing: None,
        },
    ))
}

/// Parse a bare reference sequence (e.g., TCACA after interval means c.100_104TCACA)
/// This is used by ClinVar to assert what sequence exists at a position
/// Mutalyzer interprets this as a repeat type with the sequence
#[inline]
fn parse_reference_sequence(input: &str) -> IResult<&str, NaEdit> {
    // Must be at least 2 bases to distinguish from other patterns
    // Single base would be ambiguous with substitution reference
    let (remaining, sequence) = parse_sequence(input)?;
    if sequence.len() < 2 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }
    // Treat as identity with the given sequence
    Ok((
        remaining,
        NaEdit::Identity {
            sequence: Some(sequence),
            whole_entity: false,
        },
    ))
}

/// Parse any nucleic acid edit
pub fn parse_na_edit(input: &str) -> IResult<&str, NaEdit> {
    alt((
        parse_substitution,
        parse_substitution_no_ref,    // Substitution without ref (e.g., >A)
        parse_multibase_substitution, // Must come after substitution (e.g., GG>G)
        parse_delins_with_deleted_seq, // Must come before delins and deletion (e.g., delAinsT)
        parse_delins_with_deleted_count, // Must come before delins and deletion (e.g., del35insTA)
        parse_delins,                 // Must come before deletion
        parse_deletion,
        parse_insertion,
        parse_dupins, // Must come before duplication (dupins starts with dup)
        parse_duplication,
        parse_inversion,
        parse_conversion,  // Gene conversion
        parse_methylation, // Methylation (epigenetic)
        parse_copy_number, // Copy number (copy2, copy4)
        parse_repeat,
        parse_parenthesized_repeat, // (30_125) - repeat count in parentheses
        parse_identity,
        parse_reference_sequence, // Must be last - catches bare sequence as reference statement
    ))
    .parse(input)
}

/// Parse a protein substitution (e.g., Val600Glu after the position)
fn parse_protein_substitution(input: &str) -> IResult<&str, ProteinEdit> {
    let (remaining, alternative) = parse_amino_acid.parse(input)?;
    // Allow optional trailing number (e.g., p.Arg725Trp5) - appears in ClinVar data
    // This number is typically an annotation suffix and is ignored
    let remaining = if let Ok((rest, _)) = digit1::<_, nom::error::Error<&str>>.parse(remaining) {
        rest
    } else {
        remaining
    };
    Ok((
        remaining,
        ProteinEdit::Substitution {
            reference: AminoAcid::Xaa, // Will be filled in from position
            alternative,
        },
    ))
}

/// Parse protein deletion (e.g., del, del32, delQRGEPG, delLys) or deletion-extension (e.g., delextTer?)
/// The delext pattern is used when a stop codon is deleted and the protein extends
fn parse_protein_deletion(input: &str) -> IResult<&str, ProteinEdit> {
    // First try delext (deletion followed by extension)
    if let Ok((remaining, ext)) = parse_protein_del_extension(input) {
        return Ok((remaining, ext));
    }

    // Parse "del" followed by optional count or amino acid sequence
    let (remaining, _) = tag("del").parse(input)?;

    // Try to parse a count first (e.g., del32)
    if let Ok((remaining2, count_str)) = digit1::<_, nom::error::Error<&str>>.parse(remaining) {
        if let Ok(count) = count_str.parse::<u64>() {
            return Ok((
                remaining2,
                ProteinEdit::Deletion {
                    sequence: None,
                    count: Some(count),
                },
            ));
        }
    }

    // Try to parse an amino acid sequence (e.g., delQRGEPG or delLys)
    if let Ok((remaining2, seq)) = parse_amino_acid_seq(remaining) {
        return Ok((
            remaining2,
            ProteinEdit::Deletion {
                sequence: Some(seq),
                count: None,
            },
        ));
    }

    // Plain deletion without count or sequence
    Ok((
        remaining,
        ProteinEdit::Deletion {
            sequence: None,
            count: None,
        },
    ))
}

/// Parse protein deletion-extension (e.g., delextTer?, delext*17)
/// This is used when a stop codon is deleted and the protein extends
fn parse_protein_del_extension(input: &str) -> IResult<&str, ProteinEdit> {
    let (input, _) = tag("delext").parse(input)?;

    // Parse the extension details (same as parse_protein_extension but without the 'ext' prefix)
    let (input, result) = alt((
        // delextTer? or delext*? (extension with unknown stop)
        map(tag("Ter?"), |_| (ExtDirection::CTerminal, None)),
        map(tag("*?"), |_| (ExtDirection::CTerminal, None)),
        // delextTer17 (extension with stop position)
        map(preceded(tag("Ter"), digit1), |n: &str| {
            (ExtDirection::CTerminal, Some(n.parse::<i64>().unwrap_or(0)))
        }),
        // delext*17 (extension with stop position, alternate notation)
        map(preceded(tag("*"), digit1), |n: &str| {
            (ExtDirection::CTerminal, Some(n.parse::<i64>().unwrap_or(0)))
        }),
        // delext? (unknown extension)
        map(tag("?"), |_| (ExtDirection::CTerminal, None)),
        // Just delext (defaults to C-terminal, no count)
        map(tag(""), |_| (ExtDirection::CTerminal, None)),
    ))
    .parse(input)?;

    Ok((
        input,
        ProteinEdit::Extension {
            new_aa: None, // The deleted amino acid is implicit
            direction: result.0,
            count: result.1,
        },
    ))
}

/// Parse protein duplication
fn parse_protein_duplication(input: &str) -> IResult<&str, ProteinEdit> {
    map(tag("dup"), |_| ProteinEdit::Duplication).parse(input)
}

/// Parse protein frameshift (e.g., fs, fsTer12, fs*, fsTer?, fs*?, Profs, ProfsTer23)
fn parse_protein_frameshift(input: &str) -> IResult<&str, ProteinEdit> {
    // Try to parse optional amino acid before "fs"
    let (input, new_aa) = opt(parse_amino_acid).parse(input)?;
    let (input, _) = tag("fs").parse(input)?;
    // Parse optional termination position: Ter12 or *12 or Ter? or *?
    let (input, ter_pos) = opt(alt((
        // Ter? or *? - uncertain termination position (VEP notation)
        map(tag("Ter?"), |_| Some(None)),
        map(tag("*?"), |_| Some(None)),
        // Ter12 or *12 - specific termination position
        preceded(
            tag("Ter"),
            map(digit1, |s: &str| Some(s.parse::<u64>().ok())),
        ),
        preceded(tag("*"), map(digit1, |s: &str| Some(s.parse::<u64>().ok()))),
    )))
    .parse(input)?;

    Ok((
        input,
        ProteinEdit::Frameshift {
            new_aa,
            ter_pos: ter_pos.flatten().flatten(),
        },
    ))
}

/// Parse protein identity (= or (=) for predicted)
/// Note: whole_protein flag will be set by the variant parser based on context
fn parse_protein_identity(input: &str) -> IResult<&str, ProteinEdit> {
    alt((
        // Predicted identity: (=)
        map(tag("(=)"), |_| ProteinEdit::Identity {
            predicted: true,
            whole_protein: false, // Will be updated by variant parser if no position
        }),
        // Certain identity: =
        map(tag("="), |_| ProteinEdit::Identity {
            predicted: false,
            whole_protein: false, // Will be updated by variant parser if no position
        }),
    ))
    .parse(input)
}

/// Parse protein unknown (?)
fn parse_protein_unknown(input: &str) -> IResult<&str, ProteinEdit> {
    map(tag("?"), |_| ProteinEdit::position_unknown()).parse(input)
}

/// Parse no protein production (0 or 0? for predicted)
fn parse_protein_no_protein(input: &str) -> IResult<&str, ProteinEdit> {
    alt((
        // Predicted no-protein: 0?
        map(tag("0?"), |_| ProteinEdit::NoProtein { predicted: true }),
        // Certain no-protein: 0
        map(tag("0"), |_| ProteinEdit::NoProtein { predicted: false }),
    ))
    .parse(input)
}

/// Parse a sequence of amino acids (e.g., GlnProArg)
fn parse_amino_acid_seq(input: &str) -> IResult<&str, AminoAcidSeq> {
    map(many1(parse_amino_acid), AminoAcidSeq::new).parse(input)
}

/// Parse protein insertion (e.g., insGln, insGlnProArg)
fn parse_protein_insertion(input: &str) -> IResult<&str, ProteinEdit> {
    map(preceded(tag("ins"), parse_amino_acid_seq), |sequence| {
        ProteinEdit::Insertion { sequence }
    })
    .parse(input)
}

/// Parse protein delins (e.g., delinsGlu, delinsTrpVal)
fn parse_protein_delins(input: &str) -> IResult<&str, ProteinEdit> {
    map(preceded(tag("delins"), parse_amino_acid_seq), |sequence| {
        ProteinEdit::Delins { sequence }
    })
    .parse(input)
}

/// Parse protein extension (e.g., ext-5, extTer17, extTer?, ext*?, Glnext*17)
fn parse_protein_extension(input: &str) -> IResult<&str, ProteinEdit> {
    // Try to parse optional amino acid before "ext"
    let (input, new_aa) = opt(parse_amino_acid).parse(input)?;
    let (input, _) = tag("ext").parse(input)?;

    // Try to parse the extension details
    // N-terminal: ext-5 (negative number, upstream of Met1)
    // C-terminal: extTer17 or ext*17 or extTer? or ext*? (downstream of stop)
    let (input, result) = alt((
        // ext-5 (N-terminal extension)
        map(preceded(char('-'), digit1), |n: &str| {
            (
                ExtDirection::NTerminal,
                Some(-(n.parse::<i64>().unwrap_or(0))),
            )
        }),
        // extTer? or ext*? (C-terminal extension with unknown stop - VEP notation)
        map(tag("Ter?"), |_| (ExtDirection::CTerminal, None)),
        map(tag("*?"), |_| (ExtDirection::CTerminal, None)),
        // extTer17 (C-terminal extension with stop position)
        map(preceded(tag("Ter"), digit1), |n: &str| {
            (ExtDirection::CTerminal, Some(n.parse::<i64>().unwrap_or(0)))
        }),
        // ext*17 (C-terminal extension with stop position, alternate notation)
        map(preceded(tag("*"), digit1), |n: &str| {
            (ExtDirection::CTerminal, Some(n.parse::<i64>().unwrap_or(0)))
        }),
        // ext? (unknown extension)
        map(tag("?"), |_| (ExtDirection::CTerminal, None)),
        // Just ext (defaults to C-terminal, no count)
        map(tag(""), |_| (ExtDirection::CTerminal, None)),
    ))
    .parse(input)?;

    Ok((
        input,
        ProteinEdit::Extension {
            new_aa,
            direction: result.0,
            count: result.1,
        },
    ))
}

/// Parse uncertain extension annotation pattern (e.g., (41_?), (?_41))
/// This is used in non-standard notation to indicate uncertain C-terminal extension
/// Pattern: (N_?) or (?_N) where N is a position number
fn parse_uncertain_extension_annotation(input: &str) -> IResult<&str, ProteinEdit> {
    let bytes = input.as_bytes();

    // Must start with '('
    if bytes.is_empty() || bytes[0] != b'(' {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    let inner = &bytes[1..];
    if inner.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Try (N_?) pattern - position followed by unknown end
    if inner[0].is_ascii_digit() {
        let (num, consumed) = scan_digits(inner);
        if consumed > 0 {
            let after_num = &inner[consumed..];
            // Check for _?)
            if after_num.len() >= 3
                && after_num[0] == b'_'
                && after_num[1] == b'?'
                && after_num[2] == b')'
            {
                let total_consumed = 1 + consumed + 3; // ( + num + _?)
                return Ok((
                    &input[total_consumed..],
                    ProteinEdit::Extension {
                        new_aa: None,
                        direction: ExtDirection::CTerminal,
                        count: Some(num as i64), // Start position of extension
                    },
                ));
            }
        }
    }

    // Try (?_N) pattern - unknown start followed by position
    if inner.len() >= 2 && inner[0] == b'?' && inner[1] == b'_' {
        let (num, consumed) = scan_digits(&inner[2..]);
        if consumed > 0 {
            let after_num = &inner[2 + consumed..];
            // Check for )
            if !after_num.is_empty() && after_num[0] == b')' {
                let total_consumed = 1 + 2 + consumed + 1; // ( + ?_ + num + )
                return Ok((
                    &input[total_consumed..],
                    ProteinEdit::Extension {
                        new_aa: None,
                        direction: ExtDirection::CTerminal,
                        count: Some(num as i64), // End position of extension
                    },
                ));
            }
        }
    }

    Err(nom::Err::Error(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Tag,
    )))
}

/// Parse a single protein repeat unit (e.g., PA[3], QAV[1], SGGG[3])
/// Returns the sequence and count
#[inline]
fn parse_single_protein_repeat_unit(input: &str) -> IResult<&str, (AminoAcidSeq, RepeatCount)> {
    // Parse one or more single-letter amino acids
    let mut aas = Vec::new();
    let mut remaining = input;

    while !remaining.is_empty() {
        let bytes = remaining.as_bytes();
        if bytes[0] == b'[' {
            // Stop at repeat count bracket
            break;
        }
        // Try to parse a single-letter amino acid
        if let Ok((rest, aa)) = parse_amino_acid_one_letter(remaining) {
            aas.push(aa);
            remaining = rest;
        } else {
            break;
        }
    }

    if aas.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Parse repeat count
    let (remaining, count) = parse_repeat_count(remaining)?;

    Ok((remaining, (AminoAcidSeq::new(aas), count)))
}

/// Parse a protein repeat (e.g., PA[3], QAV[1], SGGG[3])
/// Also handles multi-repeat (e.g., A[4]SAAAA[1])
/// This is used for short tandem repeat expansions with single-letter amino acid codes
#[inline]
fn parse_protein_repeat(input: &str) -> IResult<&str, ProteinEdit> {
    // Parse the first repeat unit
    let (mut remaining, first_unit) = parse_single_protein_repeat_unit(input)?;

    // Check if there are more repeat units
    let mut units = vec![first_unit];
    while !remaining.is_empty() {
        // Check if next character could start a new repeat unit (amino acid letter)
        let bytes = remaining.as_bytes();
        if !bytes[0].is_ascii_uppercase() {
            break;
        }

        // Try to parse another repeat unit
        if let Ok((rest, unit)) = parse_single_protein_repeat_unit(remaining) {
            units.push(unit);
            remaining = rest;
        } else {
            break;
        }
    }

    if units.len() == 1 {
        // Single repeat unit
        let (sequence, count) = units.pop().unwrap();
        Ok((remaining, ProteinEdit::Repeat { sequence, count }))
    } else {
        // Multiple repeat units
        Ok((remaining, ProteinEdit::MultiRepeat { units }))
    }
}

/// Parse any protein edit
/// Optimized with byte-based dispatch to avoid trying all 10 alternatives
#[inline]
pub fn parse_protein_edit(input: &str) -> IResult<&str, ProteinEdit> {
    let bytes = input.as_bytes();
    if bytes.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    // Dispatch based on first character
    match bytes[0] {
        b'f' => {
            // fs (frameshift without leading amino acid)
            if bytes.len() >= 2 && bytes[1] == b's' {
                return parse_protein_frameshift(input);
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        b'd' => {
            // delins, del, or dup - try longer matches first
            if bytes.len() >= 6 && &bytes[..6] == b"delins" {
                return parse_protein_delins(input);
            }
            if bytes.len() >= 3 && &bytes[..3] == b"del" {
                return parse_protein_deletion(input);
            }
            if bytes.len() >= 3 && &bytes[..3] == b"dup" {
                return parse_protein_duplication(input);
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        b'i' => {
            // ins (insertion)
            if bytes.len() >= 3 && &bytes[..3] == b"ins" {
                return parse_protein_insertion(input);
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        b'e' => {
            // ext (extension without leading amino acid)
            if bytes.len() >= 3 && &bytes[..3] == b"ext" {
                return parse_protein_extension(input);
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        b'=' => parse_protein_identity(input),
        b'(' => {
            // (=) predicted identity
            if bytes.len() >= 3 && bytes[1] == b'=' && bytes[2] == b')' {
                return parse_protein_identity(input);
            }
            // (N_?) uncertain extension annotation pattern (e.g., (41_?))
            // This indicates an uncertain C-terminal extension position
            if let Ok(result) = parse_uncertain_extension_annotation(input) {
                return Ok(result);
            }
            Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
        b'0' => parse_protein_no_protein(input),
        b'?' => parse_protein_unknown(input),
        b'[' => {
            // Bracket-only repeat count: [14] without preceding amino acid sequence
            // This is used when the repeat unit is inferred from the position (e.g., p.Asp38[14])
            let (remaining, count) = parse_repeat_count(input)?;
            Ok((
                remaining,
                ProteinEdit::Repeat {
                    sequence: AminoAcidSeq::new(Vec::new()),
                    count,
                },
            ))
        }
        _ => {
            // Could be amino acid leading to: frameshift (Profs), extension (Glnext),
            // duplication (dup), repeat (AA[n]), or substitution (just amino acid)
            // Try frameshift first (amino acid + fs), then extension (amino acid + ext),
            // then dup, then repeat, then substitution
            if let Ok(result) = parse_protein_frameshift(input) {
                return Ok(result);
            }
            if let Ok(result) = parse_protein_extension(input) {
                return Ok(result);
            }
            if bytes.len() >= 3 && &bytes[..3] == b"dup" {
                return parse_protein_duplication(input);
            }
            // Try protein repeat (single-letter AA sequence followed by [count])
            if let Ok(result) = parse_protein_repeat(input) {
                return Ok(result);
            }
            // Fall back to substitution (just amino acid)
            parse_protein_substitution(input)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::InsertedPart;

    #[test]
    fn test_parse_substitution() {
        let (remaining, edit) = parse_na_edit("A>G").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(edit, NaEdit::Substitution { .. }));
    }

    #[test]
    fn test_parse_substitution_iupac_ambiguity() {
        // Test IUPAC ambiguity codes in substitutions
        // H = A, C, or T (not G)
        let (remaining, edit) = parse_na_edit("G>H").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Substitution {
            reference,
            alternative,
        } = edit
        {
            assert_eq!(reference, Base::G);
            assert_eq!(alternative, Base::H);
        } else {
            panic!("Expected substitution");
        }

        // R = A or G (purine)
        let (remaining, edit) = parse_na_edit("A>R").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Substitution { alternative, .. } = edit {
            assert_eq!(alternative, Base::R);
        }

        // Y = C or T (pyrimidine)
        let (remaining, edit) = parse_na_edit("C>Y").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Substitution { alternative, .. } = edit {
            assert_eq!(alternative, Base::Y);
        }
    }

    #[test]
    fn test_parse_deletion_simple() {
        let (remaining, edit) = parse_na_edit("del").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(
            edit,
            NaEdit::Deletion {
                sequence: None,
                length: None
            }
        ));
    }

    #[test]
    fn test_parse_deletion_with_seq() {
        let (remaining, edit) = parse_na_edit("delATG").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Deletion {
            sequence: Some(seq),
            length: None,
        } = edit
        {
            assert_eq!(seq.len(), 3);
        } else {
            panic!("Expected deletion with sequence");
        }
    }

    #[test]
    fn test_parse_deletion_with_length() {
        let (remaining, edit) = parse_na_edit("del101").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Deletion {
            sequence: None,
            length: Some(len),
        } = edit
        {
            assert_eq!(len, 101);
        } else {
            panic!("Expected deletion with length");
        }
    }

    #[test]
    fn test_parse_insertion() {
        let (remaining, edit) = parse_na_edit("insATG").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            assert_eq!(sequence.len(), Some(3));
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_delins() {
        let (remaining, edit) = parse_na_edit("delinsATG").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(edit, NaEdit::Delins { .. }));
    }

    #[test]
    fn test_parse_duplication() {
        let (remaining, edit) = parse_na_edit("dup").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(
            edit,
            NaEdit::Duplication {
                sequence: None,
                length: None,
                ..
            }
        ));
    }

    #[test]
    fn test_parse_duplication_with_length() {
        let (remaining, edit) = parse_na_edit("dup101").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Duplication {
            sequence: None,
            length: Some(len),
            ..
        } = edit
        {
            assert_eq!(len, 101);
        } else {
            panic!("Expected duplication with length");
        }
    }

    #[test]
    fn test_parse_inversion() {
        let (remaining, edit) = parse_na_edit("inv").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(
            edit,
            NaEdit::Inversion {
                sequence: None,
                length: None
            }
        ));
    }

    #[test]
    fn test_parse_inversion_with_length() {
        let (remaining, edit) = parse_na_edit("inv3").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Inversion {
            sequence: None,
            length: Some(len),
        } = edit
        {
            assert_eq!(len, 3);
        } else {
            panic!("Expected inversion with length");
        }
    }

    #[test]
    fn test_parse_inversion_with_sequence() {
        let (remaining, edit) = parse_na_edit("invATG").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Inversion {
            sequence: Some(seq),
            length: None,
        } = edit
        {
            assert_eq!(seq.len(), 3);
        } else {
            panic!("Expected inversion with sequence");
        }
    }

    #[test]
    fn test_parse_repeat_exact() {
        let (remaining, edit) = parse_na_edit("CAG[12]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Repeat {
            count,
            additional_counts,
            ..
        } = edit
        {
            assert!(matches!(count, RepeatCount::Exact(12)));
            assert!(additional_counts.is_empty());
        } else {
            panic!("Expected repeat");
        }
    }

    #[test]
    fn test_parse_repeat_range() {
        let (remaining, edit) = parse_na_edit("[10_15]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Repeat {
            count,
            additional_counts,
            ..
        } = edit
        {
            assert!(matches!(count, RepeatCount::Range(10, 15)));
            assert!(additional_counts.is_empty());
        } else {
            panic!("Expected repeat");
        }
    }

    #[test]
    fn test_parse_repeat_genotype() {
        // Nested repeats for genotype notation: A[6][1]
        let (remaining, edit) = parse_na_edit("A[6][1]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Repeat {
            sequence,
            count,
            additional_counts,
            trailing,
        } = edit
        {
            assert!(sequence.is_some());
            assert!(matches!(count, RepeatCount::Exact(6)));
            assert_eq!(additional_counts.len(), 1);
            assert!(matches!(additional_counts[0], RepeatCount::Exact(1)));
            assert!(trailing.is_none());
        } else {
            panic!("Expected repeat");
        }

        // Three repeat counts: [4][5][6]
        let (remaining, edit) = parse_na_edit("[4][5][6]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Repeat {
            count,
            additional_counts,
            ..
        } = edit
        {
            assert!(matches!(count, RepeatCount::Exact(4)));
            assert_eq!(additional_counts.len(), 2);
            assert!(matches!(additional_counts[0], RepeatCount::Exact(5)));
            assert!(matches!(additional_counts[1], RepeatCount::Exact(6)));
        } else {
            panic!("Expected repeat");
        }
    }

    #[test]
    fn test_parse_identity() {
        let (remaining, edit) = parse_na_edit("=").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(edit, NaEdit::Identity { .. }));
    }

    #[test]
    fn test_parse_identity_with_ref_base() {
        // HGVS format: G= (reference base before equals sign)
        let (remaining, edit) = parse_na_edit("G=").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Identity {
            sequence: Some(seq),
            whole_entity: false,
        } = edit
        {
            assert_eq!(seq.len(), 1);
        } else {
            panic!("Expected identity with sequence");
        }

        // Multiple bases before equals
        let (remaining, edit) = parse_na_edit("ATG=").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Identity {
            sequence: Some(seq),
            whole_entity: false,
        } = edit
        {
            assert_eq!(seq.len(), 3);
        } else {
            panic!("Expected identity with sequence");
        }
    }

    #[test]
    fn test_parse_reference_sequence() {
        // ClinVar format: bare sequence after interval (e.g., c.100_104TCACA)
        // This represents a reference statement (the sequence at this position is X)
        let (remaining, edit) = parse_na_edit("TCACA").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Identity {
            sequence: Some(seq),
            whole_entity: false,
        } = edit
        {
            assert_eq!(seq.len(), 5);
            assert_eq!(seq.to_string(), "TCACA");
        } else {
            panic!("Expected identity with sequence, got {:?}", edit);
        }

        // Longer sequence
        let (remaining, edit) = parse_na_edit("ATGCATGC").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Identity {
            sequence: Some(seq),
            ..
        } = edit
        {
            assert_eq!(seq.len(), 8);
        } else {
            panic!("Expected identity with sequence");
        }
    }

    #[test]
    fn test_parse_protein_substitution() {
        let (remaining, edit) = parse_protein_edit("Glu").unwrap();
        assert_eq!(remaining, "");
        assert!(matches!(edit, ProteinEdit::Substitution { .. }));
    }

    #[test]
    fn test_parse_protein_frameshift() {
        let (remaining, edit) = parse_protein_edit("fsTer12").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Frameshift { new_aa, ter_pos } = edit {
            assert_eq!(new_aa, None);
            assert_eq!(ter_pos, Some(12));
        } else {
            panic!("Expected frameshift");
        }

        // Frameshift with new amino acid
        let (remaining, edit) = parse_protein_edit("ProfsTer23").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Frameshift { new_aa, ter_pos } = edit {
            assert_eq!(new_aa, Some(crate::hgvs::location::AminoAcid::Pro));
            assert_eq!(ter_pos, Some(23));
        } else {
            panic!("Expected frameshift");
        }
    }

    #[test]
    fn test_parse_protein_insertion() {
        // Single amino acid
        let (remaining, edit) = parse_protein_edit("insGln").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Insertion { sequence } = edit {
            assert_eq!(sequence.len(), 1);
        } else {
            panic!("Expected insertion");
        }

        // Multiple amino acids
        let (remaining, edit) = parse_protein_edit("insGlyPro").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Insertion { sequence } = edit {
            assert_eq!(sequence.len(), 2);
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_protein_delins() {
        let (remaining, edit) = parse_protein_edit("delinsTrpVal").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Delins { sequence } = edit {
            assert_eq!(sequence.len(), 2);
        } else {
            panic!("Expected delins");
        }
    }

    #[test]
    fn test_parse_protein_extension() {
        // N-terminal extension
        let (remaining, edit) = parse_protein_edit("ext-5").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Extension {
            new_aa,
            direction,
            count,
        } = edit
        {
            assert_eq!(new_aa, None);
            assert_eq!(direction, ExtDirection::NTerminal);
            assert_eq!(count, Some(-5));
        } else {
            panic!("Expected extension");
        }

        // C-terminal extension with Ter
        let (remaining, edit) = parse_protein_edit("extTer17").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Extension {
            new_aa,
            direction,
            count,
        } = edit
        {
            assert_eq!(new_aa, None);
            assert_eq!(direction, ExtDirection::CTerminal);
            assert_eq!(count, Some(17));
        } else {
            panic!("Expected extension");
        }

        // C-terminal extension with *
        let (remaining, edit) = parse_protein_edit("ext*?").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Extension {
            new_aa,
            direction,
            count,
        } = edit
        {
            assert_eq!(new_aa, None);
            assert_eq!(direction, ExtDirection::CTerminal);
            assert_eq!(count, None);
        } else {
            panic!("Expected extension");
        }

        // C-terminal extension with new amino acid
        let (remaining, edit) = parse_protein_edit("Glnext*17").unwrap();
        assert_eq!(remaining, "");
        if let ProteinEdit::Extension {
            new_aa,
            direction,
            count,
        } = edit
        {
            assert_eq!(new_aa, Some(crate::hgvs::location::AminoAcid::Gln));
            assert_eq!(direction, ExtDirection::CTerminal);
            assert_eq!(count, Some(17));
        } else {
            panic!("Expected extension");
        }
    }

    #[test]
    fn test_parse_insertion_count() {
        // Simple count: ins10
        let (remaining, edit) = parse_na_edit("ins10").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            assert!(matches!(sequence, InsertedSequence::Count(10)));
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_insertion_parenthesized_count() {
        // Parenthesized count: ins(10)
        let (remaining, edit) = parse_na_edit("ins(10)").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            assert!(matches!(sequence, InsertedSequence::Count(10)));
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_insertion_range() {
        // Range: ins(10_20)
        let (remaining, edit) = parse_na_edit("ins(10_20)").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            assert!(matches!(sequence, InsertedSequence::Range(10, 20)));
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_insertion_repeat() {
        // Repeated base: insA[10]
        let (remaining, edit) = parse_na_edit("insA[10]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Repeat { base, count } = sequence {
                assert_eq!(base, Base::A);
                assert!(matches!(count, RepeatCount::Exact(10)));
            } else {
                panic!("Expected repeat insertion");
            }
        } else {
            panic!("Expected insertion");
        }

        // N repeat: insN[15]
        let (remaining, edit) = parse_na_edit("insN[15]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Repeat { base, count } = sequence {
                assert_eq!(base, Base::N);
                assert!(matches!(count, RepeatCount::Exact(15)));
            } else {
                panic!("Expected repeat insertion");
            }
        } else {
            panic!("Expected insertion");
        }

        // N repeat with range: insN[15_30]
        let (remaining, edit) = parse_na_edit("insN[15_30]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Repeat { base, count } = sequence {
                assert_eq!(base, Base::N);
                assert!(matches!(count, RepeatCount::Range(15, 30)));
            } else {
                panic!("Expected repeat insertion");
            }
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_insertion_complex() {
        // Complex: ins[A[10];T]
        let (remaining, edit) = parse_na_edit("ins[A[10];T]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Complex(parts) = sequence {
                assert_eq!(parts.len(), 2);
            } else {
                panic!("Expected complex insertion");
            }
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_delins_count() {
        // Delins with count: delins10
        let (remaining, edit) = parse_na_edit("delins10").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Delins { sequence } = edit {
            assert!(matches!(sequence, InsertedSequence::Count(10)));
        } else {
            panic!("Expected delins");
        }
    }

    #[test]
    fn test_parse_delins_n_repeat() {
        // Delins with N[count]: delinsN[12]
        let (remaining, edit) = parse_na_edit("delinsN[12]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Delins { sequence } = edit {
            if let InsertedSequence::Repeat { base, count } = sequence {
                assert_eq!(base, Base::N);
                assert!(matches!(count, RepeatCount::Exact(12)));
            } else {
                panic!("Expected repeat delins");
            }
        } else {
            panic!("Expected delins");
        }
    }

    #[test]
    fn test_parse_methylation_gom() {
        let (remaining, edit) = parse_na_edit("|gom").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Methylation { status } = edit {
            assert_eq!(status, MethylationStatus::GainOfMethylation);
        } else {
            panic!("Expected methylation");
        }
    }

    #[test]
    fn test_parse_methylation_lom() {
        let (remaining, edit) = parse_na_edit("|lom").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Methylation { status } = edit {
            assert_eq!(status, MethylationStatus::LossOfMethylation);
        } else {
            panic!("Expected methylation");
        }
    }

    #[test]
    fn test_parse_methylation_unchanged() {
        let (remaining, edit) = parse_na_edit("|met=").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Methylation { status } = edit {
            assert_eq!(status, MethylationStatus::Unchanged);
        } else {
            panic!("Expected methylation");
        }
    }

    #[test]
    fn test_parse_insertion_external_sequence_reference() {
        // GenBank accession with genomic position range: ins[PQ998981.1:g.1_6057]
        let (remaining, edit) = parse_na_edit("ins[PQ998981.1:g.1_6057]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Reference(ref_str) = sequence {
                assert_eq!(ref_str, "PQ998981.1:g.1_6057");
            } else {
                panic!("Expected reference insertion, got {:?}", sequence);
            }
        } else {
            panic!("Expected insertion");
        }

        // Two-letter GenBank prefix: ins[MF045863.1:g.1_36978]
        let (remaining, edit) = parse_na_edit("ins[MF045863.1:g.1_36978]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Reference(ref_str) = sequence {
                assert_eq!(ref_str, "MF045863.1:g.1_36978");
            } else {
                panic!("Expected reference insertion, got {:?}", sequence);
            }
        } else {
            panic!("Expected insertion");
        }

        // MT prefix (mitochondrial): ins[MT113356.1:g.1_2409]
        let (remaining, edit) = parse_na_edit("ins[MT113356.1:g.1_2409]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Reference(ref_str) = sequence {
                assert_eq!(ref_str, "MT113356.1:g.1_2409");
            } else {
                panic!("Expected reference insertion, got {:?}", sequence);
            }
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_is_reference_accession_prefix() {
        // RefSeq prefixes
        assert!(is_reference_accession_prefix("NC_000001.11:g.123"));
        assert!(is_reference_accession_prefix("NM_000001.1:c.123"));
        assert!(is_reference_accession_prefix("NG_000001.1:g.123"));

        // Ensembl prefixes
        assert!(is_reference_accession_prefix("ENST00000123456:c.123"));
        assert!(is_reference_accession_prefix("ENSP00000123456:p.123"));

        // LRG
        assert!(is_reference_accession_prefix("LRG_1:g.123"));

        // GenBank two-letter prefixes
        assert!(is_reference_accession_prefix("PQ998981.1:g.1"));
        assert!(is_reference_accession_prefix("MF045863.1:g.1"));
        assert!(is_reference_accession_prefix("MT113356.1:g.1"));
        assert!(is_reference_accession_prefix("KT192064.1:1"));
        assert!(is_reference_accession_prefix("KY923049.1:g.1"));
        assert!(is_reference_accession_prefix("DQ831669.1:1"));
        assert!(is_reference_accession_prefix("AC010542.7:g.1"));
        assert!(is_reference_accession_prefix("AB191243.1:g.1"));
        assert!(is_reference_accession_prefix("PP887427.1:g.1"));

        // GenBank single-letter prefixes
        assert!(is_reference_accession_prefix("U12345.1:g.1"));

        // Not valid prefixes
        assert!(!is_reference_accession_prefix("123ABC:g.1"));
        assert!(!is_reference_accession_prefix(":g.1"));
    }

    #[test]
    fn test_parse_complex_insertion_with_external_ref() {
        // Complex insertion with literal and external reference: ins[TCTT;KT192064.1:1_310]
        let (remaining, edit) = parse_na_edit("ins[TCTT;KT192064.1:1_310]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Complex(parts) = sequence {
                assert_eq!(parts.len(), 2);
                assert!(matches!(&parts[0], InsertedPart::Literal(_)));
                if let InsertedPart::ExternalRef(ref_str) = &parts[1] {
                    assert_eq!(ref_str, "KT192064.1:1_310");
                } else {
                    panic!("Expected external reference, got {:?}", parts[1]);
                }
            } else {
                panic!("Expected complex insertion, got {:?}", sequence);
            }
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_complex_delins_with_inversion() {
        // Complex delins with position ranges and inversion
        let (remaining, edit) =
            parse_na_edit("delins[AGAAGGAAATTT;45310743_46521014;45043709_45310738inv]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Delins { sequence } = edit {
            if let InsertedSequence::Complex(parts) = sequence {
                assert_eq!(parts.len(), 3);
                assert!(matches!(&parts[0], InsertedPart::Literal(_)));
                assert!(matches!(
                    &parts[1],
                    InsertedPart::PositionRange {
                        start: 45310743,
                        end: 46521014
                    }
                ));
                assert!(matches!(
                    &parts[2],
                    InsertedPart::PositionRangeInv {
                        start: 45043709,
                        end: 45310738
                    }
                ));
            } else {
                panic!("Expected complex delins, got {:?}", sequence);
            }
        } else {
            panic!("Expected delins");
        }
    }

    #[test]
    fn test_parse_insertion_with_multiple_external_refs() {
        // ins[80114172_80114186;NC_000020.11:g.2823027_2826302;AAA]
        let (remaining, edit) =
            parse_na_edit("ins[80114172_80114186;NC_000020.11:g.2823027_2826302;AAA]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::Complex(parts) = sequence {
                assert_eq!(parts.len(), 3);
                assert!(matches!(
                    &parts[0],
                    InsertedPart::PositionRange {
                        start: 80114172,
                        end: 80114186
                    }
                ));
                if let InsertedPart::ExternalRef(ref_str) = &parts[1] {
                    assert_eq!(ref_str, "NC_000020.11:g.2823027_2826302");
                } else {
                    panic!("Expected external reference, got {:?}", parts[1]);
                }
                assert!(matches!(&parts[2], InsertedPart::Literal(_)));
            } else {
                panic!("Expected complex insertion, got {:?}", sequence);
            }
        } else {
            panic!("Expected insertion");
        }
    }

    #[test]
    fn test_parse_substitution_no_ref() {
        // Substitution without reference base: >A
        let (remaining, edit) = parse_na_edit(">A").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::SubstitutionNoRef { alternative } = edit {
            assert_eq!(alternative, Base::A);
        } else {
            panic!("Expected substitution without ref, got {:?}", edit);
        }

        // Other bases
        let (remaining, edit) = parse_na_edit(">G").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::SubstitutionNoRef { alternative } = edit {
            assert_eq!(alternative, Base::G);
        } else {
            panic!("Expected substitution without ref");
        }

        let (remaining, edit) = parse_na_edit(">C").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::SubstitutionNoRef { alternative } = edit {
            assert_eq!(alternative, Base::C);
        } else {
            panic!("Expected substitution without ref");
        }

        let (remaining, edit) = parse_na_edit(">T").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::SubstitutionNoRef { alternative } = edit {
            assert_eq!(alternative, Base::T);
        } else {
            panic!("Expected substitution without ref");
        }
    }

    #[test]
    fn test_parse_delins_sequence_repeat() {
        // Delins with sequence repeat count: delinsTCGGCAGCGGCACAGCGAGG[13]
        let (remaining, edit) = parse_na_edit("delinsTCGGCAGCGGCACAGCGAGG[13]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Delins { sequence } = edit {
            if let InsertedSequence::SequenceRepeat {
                sequence: seq,
                count,
            } = sequence
            {
                assert_eq!(seq.len(), 20);
                assert!(matches!(count, RepeatCount::Exact(13)));
            } else {
                panic!("Expected sequence repeat delins, got {:?}", sequence);
            }
        } else {
            panic!("Expected delins");
        }

        // Delins with sequence repeat range: delinsAAAGG[400_2000]
        let (remaining, edit) = parse_na_edit("delinsAAAGG[400_2000]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Delins { sequence } = edit {
            if let InsertedSequence::SequenceRepeat {
                sequence: seq,
                count,
            } = sequence
            {
                assert_eq!(seq.len(), 5);
                assert!(matches!(count, RepeatCount::Range(400, 2000)));
            } else {
                panic!("Expected sequence repeat delins, got {:?}", sequence);
            }
        } else {
            panic!("Expected delins");
        }
    }

    #[test]
    fn test_parse_insertion_sequence_repeat() {
        // Insertion with sequence repeat count: insTCGGCAGCGGCACAGCGAGG[13]
        let (remaining, edit) = parse_na_edit("insTCGGCAGCGGCACAGCGAGG[13]").unwrap();
        assert_eq!(remaining, "");
        if let NaEdit::Insertion { sequence } = edit {
            if let InsertedSequence::SequenceRepeat {
                sequence: seq,
                count,
            } = sequence
            {
                assert_eq!(seq.len(), 20);
                assert!(matches!(count, RepeatCount::Exact(13)));
            } else {
                panic!("Expected sequence repeat insertion, got {:?}", sequence);
            }
        } else {
            panic!("Expected insertion");
        }
    }
}
