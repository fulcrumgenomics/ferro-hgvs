//! Position parsing
//!
//! Parses positions in different coordinate systems (genomic, CDS, transcript).

use crate::hgvs::location::{AminoAcid, CdsPos, GenomePos, ProtPos, RnaPos, TxPos};
use nom::{
    branch::alt,
    character::complete::{char, digit1},
    combinator::opt,
    IResult, Parser,
};

/// Parse a genomic position (unsigned integer >= 1, or special marker like pter/qter/cen)
///
/// Also supports optional offsets for uncertain intronic-style notation (e.g., g.12345-? or g.67890+?).
///
/// # Validation
///
/// Position must be >= 1 (HGVS positions are 1-based). Position 0 is rejected with an error.
///
/// **Note on offsets**: When an offset is specified (e.g., `g.12345+10` or `g.100-5`), the parser
/// stores the base position and offset separately. Semantic validation (e.g., whether `position + offset`
/// resolves to a valid genomic coordinate) is NOT performed at parse time because:
/// 1. Genomic offsets are primarily used for uncertain intronic-style notation
/// 2. Validation requires reference sequence context to determine valid ranges
/// 3. The offset may represent uncertainty (`+?` or `-?`) rather than a concrete value
///
/// Callers that need to validate resolved positions should do so after parsing.
#[inline]
pub fn parse_genome_pos(input: &str) -> IResult<&str, GenomePos> {
    // Check for special position markers first
    let bytes = input.as_bytes();
    if bytes.len() >= 4 {
        if bytes[0..4] == *b"pter" {
            return Ok((&input[4..], GenomePos::pter()));
        }
        if bytes[0..4] == *b"qter" {
            return Ok((&input[4..], GenomePos::qter()));
        }
    }
    if bytes.len() >= 3 && bytes[0..3] == *b"cen" {
        return Ok((&input[3..], GenomePos::cen()));
    }

    // Regular numeric position
    let (remaining, s) = digit1.parse(input)?;
    // Use checked parsing to detect overflow (returns error instead of silent 0)
    let base: u64 = s.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    // HGVS positions are 1-based, position 0 is invalid
    if base == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }

    // Try to parse optional offset (for intronic-style notation in genomic coordinates)
    let (remaining, offset) = opt(parse_offset).parse(remaining)?;

    match offset {
        Some(off) => Ok((remaining, GenomePos::with_offset(base, off))),
        None => Ok((remaining, GenomePos::new(base))),
    }
}

/// Parse an intronic offset (+5, -10, etc.)
///
/// # Sentinel Values for Unknown Offsets
///
/// HGVS supports uncertain offsets written as `+?` or `-?`. Since offsets use `i64`,
/// we use sentinel values to represent these uncertain offsets:
///
/// - [`OFFSET_UNKNOWN_POSITIVE`] (`i64::MAX`): Represents `+?` (unknown positive offset)
/// - [`OFFSET_UNKNOWN_NEGATIVE`] (`i64::MIN`): Represents `-?` (unknown negative offset)
///
/// These sentinel values were chosen because:
/// 1. Real genomic offsets never approach `i64::MAX` or `i64::MIN` in practice
/// 2. Using the type's extreme values avoids the need for an Option wrapper
/// 3. Callers can check for unknown offsets by comparing against these constants
///
/// # Example
///
/// ```rust,ignore
/// if offset == OFFSET_UNKNOWN_POSITIVE {
///     // Handle unknown positive offset (+?)
/// }
/// ```
/// Sentinel value for unknown positive offset (`+?` in HGVS notation)
pub const OFFSET_UNKNOWN_POSITIVE: i64 = i64::MAX;
/// Sentinel value for unknown negative offset (`-?` in HGVS notation)
pub const OFFSET_UNKNOWN_NEGATIVE: i64 = i64::MIN;

#[inline]
fn parse_offset(input: &str) -> IResult<&str, i64> {
    let (input, sign) = alt((char('+'), char('-'))).parse(input)?;

    // Check for uncertain offset: +? or -?
    if let Ok((remaining, _)) = char::<_, nom::error::Error<&str>>('?').parse(input) {
        let offset = if sign == '+' {
            OFFSET_UNKNOWN_POSITIVE
        } else {
            OFFSET_UNKNOWN_NEGATIVE
        };
        return Ok((remaining, offset));
    }

    let (input, num) = digit1.parse(input)?;
    // Use checked parsing to detect overflow (returns error instead of silent 0)
    let value: i64 = num.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    let signed_value = if sign == '-' { -value } else { value };

    Ok((input, signed_value))
}

/// Parse a CDS position (can include offset, negative positions, * notation)
/// Note: Position 0 is invalid in HGVS CDS notation (we go from -1 to 1)
/// Also supports ? for unknown positions (e.g., c.?-232 for unknown position with offset)
#[inline]
pub fn parse_cds_pos(input: &str) -> IResult<&str, CdsPos> {
    // ? or ?+5 or ?-232 (unknown position with optional offset)
    if let Some(rest) = input.strip_prefix('?') {
        let (remaining, offset) = opt(parse_offset).parse(rest)?;
        return Ok((remaining, CdsPos::unknown(offset)));
    }

    // *123 or *123+5 (3' UTR positions) - position 0 is invalid here too
    if let Some(rest) = input.strip_prefix('*') {
        let (remaining, s) = digit1.parse(rest)?;
        // Use checked parsing to detect overflow
        let base: i64 = s.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
        })?;
        if base == 0 {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (remaining, offset) = opt(parse_offset).parse(remaining)?;
        return Ok((
            remaining,
            CdsPos {
                base,
                offset,
                utr3: true,
            },
        ));
    }

    // -123 or -123+5 (5' UTR positions) - negative, so 0 not possible from this format
    if let Some(rest) = input.strip_prefix('-') {
        let (remaining, s) = digit1.parse(rest)?;
        // Use checked parsing to detect overflow
        let parsed: i64 = s.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
        })?;
        let base: i64 = -parsed;
        // -0 would be parsed but is semantically 0, which is invalid
        if base == 0 {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (remaining, offset) = opt(parse_offset).parse(remaining)?;
        return Ok((
            remaining,
            CdsPos {
                base,
                offset,
                utr3: false,
            },
        ));
    }

    // 123 or 123+5 (normal CDS positions) - position 0 is invalid
    let (remaining, s) = digit1.parse(input)?;
    // Use checked parsing to detect overflow
    let base: i64 = s.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    if base == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }
    let (remaining, offset) = opt(parse_offset).parse(remaining)?;
    Ok((
        remaining,
        CdsPos {
            base,
            offset,
            utr3: false,
        },
    ))
}

/// Parse a transcript position (n. coordinates)
///
/// Supports:
/// - Positive positions: n.100, n.100+5
/// - Negative positions (upstream of transcript): n.-30, n.-30+5
/// - Downstream positions: n.*5, n.*5+10 (after transcript end)
///
/// Note: Position 0 is invalid in HGVS notation (goes from -1 to 1)
#[inline]
pub fn parse_tx_pos(input: &str) -> IResult<&str, TxPos> {
    // Handle downstream positions (after transcript end)
    if let Some(rest) = input.strip_prefix('*') {
        let (remaining, s) = digit1.parse(rest)?;
        // Use checked parsing to detect overflow
        let base: i64 = s.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
        })?;
        if base == 0 {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (remaining, offset) = opt(parse_offset).parse(remaining)?;
        return Ok((
            remaining,
            TxPos {
                base,
                offset,
                downstream: true,
            },
        ));
    }

    // Handle negative positions (upstream of transcript start)
    if let Some(rest) = input.strip_prefix('-') {
        let (remaining, s) = digit1.parse(rest)?;
        // Use checked parsing to detect overflow
        let parsed: i64 = s.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
        })?;
        let base: i64 = -parsed;
        // -0 would be parsed but is semantically 0, which is invalid
        if base == 0 {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (remaining, offset) = opt(parse_offset).parse(remaining)?;
        return Ok((
            remaining,
            TxPos {
                base,
                offset,
                downstream: false,
            },
        ));
    }

    // Positive positions
    let (remaining, s) = digit1.parse(input)?;
    // Use checked parsing to detect overflow
    let base: i64 = s.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    // HGVS positions are 1-based, position 0 is invalid
    if base == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }
    let (remaining, offset) = opt(parse_offset).parse(remaining)?;
    Ok((
        remaining,
        TxPos {
            base,
            offset,
            downstream: false,
        },
    ))
}

/// Parse an RNA position (r. coordinates)
///
/// RNA positions support UTR notation like CDS positions:
/// - r.-14 for 5' UTR positions
/// - r.*41 for 3' UTR positions
/// - r.100 for normal coding positions
///
/// Note: Position 0 is invalid in HGVS RNA notation
#[inline]
pub fn parse_rna_pos(input: &str) -> IResult<&str, RnaPos> {
    // *123 or *123+5 (3' UTR positions) - position 0 is invalid
    if let Some(rest) = input.strip_prefix('*') {
        let (remaining, s) = digit1.parse(rest)?;
        // Use checked parsing to detect overflow
        let base: i64 = s.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
        })?;
        if base == 0 {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (remaining, offset) = opt(parse_offset).parse(remaining)?;
        return Ok((
            remaining,
            RnaPos {
                base,
                offset,
                utr3: true,
            },
        ));
    }

    // -123 or -123+5 (5' UTR positions)
    if let Some(rest) = input.strip_prefix('-') {
        let (remaining, s) = digit1.parse(rest)?;
        // Use checked parsing to detect overflow
        let parsed: i64 = s.parse().map_err(|_| {
            nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
        })?;
        let base: i64 = -parsed;
        if base == 0 {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
        let (remaining, offset) = opt(parse_offset).parse(remaining)?;
        return Ok((
            remaining,
            RnaPos {
                base,
                offset,
                utr3: false,
            },
        ));
    }

    // 123 or 123+5 (normal RNA positions) - position 0 is invalid
    let (remaining, s) = digit1.parse(input)?;
    // Use checked parsing to detect overflow
    let base: i64 = s.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    if base == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }
    let (remaining, offset) = opt(parse_offset).parse(remaining)?;
    Ok((
        remaining,
        RnaPos {
            base,
            offset,
            utr3: false,
        },
    ))
}

/// Parse a 3-letter amino acid code (case-insensitive)
///
/// Accepts any case combination: "Val", "VAL", "val", "vAl" all parse to Val.
#[inline]
fn parse_amino_acid_three_letter(input: &str) -> IResult<&str, AminoAcid> {
    // Check if we have at least 3 ASCII characters
    let bytes = input.as_bytes();
    if bytes.len() < 3 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )));
    }

    // Normalize to title case using a fixed-size array (no allocation)
    let normalized: [u8; 3] = [
        bytes[0].to_ascii_uppercase(),
        bytes[1].to_ascii_lowercase(),
        bytes[2].to_ascii_lowercase(),
    ];

    // Match statement - compiler optimizes to jump table or decision tree
    let aa = match &normalized {
        b"Ala" => AminoAcid::Ala,
        b"Arg" => AminoAcid::Arg,
        b"Asn" => AminoAcid::Asn,
        b"Asp" => AminoAcid::Asp,
        b"Cys" => AminoAcid::Cys,
        b"Gln" => AminoAcid::Gln,
        b"Glu" => AminoAcid::Glu,
        b"Gly" => AminoAcid::Gly,
        b"His" => AminoAcid::His,
        b"Ile" => AminoAcid::Ile,
        b"Leu" => AminoAcid::Leu,
        b"Lys" => AminoAcid::Lys,
        b"Met" => AminoAcid::Met,
        b"Phe" => AminoAcid::Phe,
        b"Pro" => AminoAcid::Pro,
        b"Pyl" => AminoAcid::Pyl,
        b"Sec" => AminoAcid::Sec,
        b"Ser" => AminoAcid::Ser,
        b"Thr" => AminoAcid::Thr,
        b"Trp" => AminoAcid::Trp,
        b"Tyr" => AminoAcid::Tyr,
        b"Val" => AminoAcid::Val,
        b"Ter" => AminoAcid::Ter,
        b"Xaa" => AminoAcid::Xaa,
        _ => {
            return Err(nom::Err::Error(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Tag,
            )))
        }
    };

    Ok((&input[3..], aa))
}

/// Parse a 1-letter amino acid code
#[inline]
pub fn parse_amino_acid_one_letter(input: &str) -> IResult<&str, AminoAcid> {
    if input.is_empty() {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Eof,
        )));
    }

    let c = input.chars().next().unwrap();
    if let Some(aa) = AminoAcid::from_one_letter(c) {
        Ok((&input[1..], aa))
    } else {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Tag,
        )))
    }
}

/// Parse an amino acid code (3-letter or 1-letter)
///
/// Tries 3-letter codes first, then falls back to 1-letter codes.
/// This ensures unambiguous parsing (e.g., "Val" is parsed as Val, not V+al).
/// Optimized: tries 3-letter first only if we have enough characters.
#[inline]
pub fn parse_amino_acid(input: &str) -> IResult<&str, AminoAcid> {
    // Try 3-letter first if we have enough characters
    if input.len() >= 3 {
        if let Ok(result) = parse_amino_acid_three_letter(input) {
            return Ok(result);
        }
    }
    // Fall back to 1-letter
    parse_amino_acid_one_letter(input)
}

/// Parse a protein position (e.g., Met1, Val600)
/// Note: Position 0 is invalid in HGVS protein notation
#[inline]
pub fn parse_prot_pos(input: &str) -> IResult<&str, ProtPos> {
    let (remaining, aa) = parse_amino_acid.parse(input)?;
    let (remaining, s) = digit1.parse(remaining)?;
    // Use checked parsing to detect overflow
    let num: u64 = s.parse().map_err(|_| {
        nom::Err::Error(nom::error::Error::new(input, nom::error::ErrorKind::Digit))
    })?;
    // HGVS positions are 1-based, position 0 is invalid
    if num == 0 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Verify,
        )));
    }
    Ok((remaining, ProtPos::new(aa, num)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genome_pos() {
        let (remaining, pos) = parse_genome_pos("12345del").unwrap();
        assert_eq!(remaining, "del");
        assert_eq!(pos.base, 12345);
    }

    #[test]
    fn test_parse_cds_pos_simple() {
        let (remaining, pos) = parse_cds_pos("100del").unwrap();
        assert_eq!(remaining, "del");
        assert_eq!(pos.base, 100);
        assert_eq!(pos.offset, None);
        assert!(!pos.utr3);
    }

    #[test]
    fn test_parse_cds_pos_with_positive_offset() {
        let (remaining, pos) = parse_cds_pos("100+5G>A").unwrap();
        assert_eq!(remaining, "G>A");
        assert_eq!(pos.base, 100);
        assert_eq!(pos.offset, Some(5));
    }

    #[test]
    fn test_parse_cds_pos_with_negative_offset() {
        let (remaining, pos) = parse_cds_pos("100-10G>A").unwrap();
        assert_eq!(remaining, "G>A");
        assert_eq!(pos.base, 100);
        assert_eq!(pos.offset, Some(-10));
    }

    #[test]
    fn test_parse_cds_pos_5utr() {
        let (remaining, pos) = parse_cds_pos("-20G>A").unwrap();
        assert_eq!(remaining, "G>A");
        assert_eq!(pos.base, -20);
        assert!(!pos.utr3);
    }

    #[test]
    fn test_parse_cds_pos_3utr() {
        let (remaining, pos) = parse_cds_pos("*50G>A").unwrap();
        assert_eq!(remaining, "G>A");
        assert_eq!(pos.base, 50);
        assert!(pos.utr3);
    }

    #[test]
    fn test_parse_amino_acid() {
        let (remaining, aa) = parse_amino_acid("Met1").unwrap();
        assert_eq!(remaining, "1");
        assert_eq!(aa, AminoAcid::Met);
    }

    #[test]
    fn test_parse_prot_pos() {
        let (remaining, pos) = parse_prot_pos("Val600Glu").unwrap();
        assert_eq!(remaining, "Glu");
        assert_eq!(pos.aa, AminoAcid::Val);
        assert_eq!(pos.number, 600);
    }

    #[test]
    fn test_parse_amino_acid_single_letter() {
        // Single letter codes
        let (remaining, aa) = parse_amino_acid("V600").unwrap();
        assert_eq!(remaining, "600");
        assert_eq!(aa, AminoAcid::Val);

        let (remaining, aa) = parse_amino_acid("E").unwrap();
        assert_eq!(remaining, "");
        assert_eq!(aa, AminoAcid::Glu);

        let (remaining, aa) = parse_amino_acid("M1").unwrap();
        assert_eq!(remaining, "1");
        assert_eq!(aa, AminoAcid::Met);

        let (remaining, aa) = parse_amino_acid("*").unwrap();
        assert_eq!(remaining, "");
        assert_eq!(aa, AminoAcid::Ter);
    }

    #[test]
    fn test_parse_amino_acid_prefers_three_letter() {
        // Should parse as "Val" not "V" + "al"
        let (remaining, aa) = parse_amino_acid("Val600").unwrap();
        assert_eq!(remaining, "600");
        assert_eq!(aa, AminoAcid::Val);

        // Should parse as "Met" not "M" + "et"
        let (remaining, aa) = parse_amino_acid("Met1").unwrap();
        assert_eq!(remaining, "1");
        assert_eq!(aa, AminoAcid::Met);
    }

    #[test]
    fn test_parse_prot_pos_single_letter() {
        // Single letter format: V600E
        let (remaining, pos) = parse_prot_pos("V600E").unwrap();
        assert_eq!(remaining, "E");
        assert_eq!(pos.aa, AminoAcid::Val);
        assert_eq!(pos.number, 600);

        let (remaining, pos) = parse_prot_pos("M1I").unwrap();
        assert_eq!(remaining, "I");
        assert_eq!(pos.aa, AminoAcid::Met);
        assert_eq!(pos.number, 1);
    }

    #[test]
    fn test_parse_cds_pos_unknown() {
        use crate::hgvs::location::CDS_BASE_UNKNOWN;

        // Simple unknown position: ?
        let (remaining, pos) = parse_cds_pos("?dup").unwrap();
        assert_eq!(remaining, "dup");
        assert_eq!(pos.base, CDS_BASE_UNKNOWN);
        assert!(pos.is_unknown());
        assert_eq!(pos.offset, None);

        // Unknown position with negative offset: ?-232
        let (remaining, pos) = parse_cds_pos("?-232_4484+?del").unwrap();
        assert_eq!(remaining, "_4484+?del");
        assert_eq!(pos.base, CDS_BASE_UNKNOWN);
        assert!(pos.is_unknown());
        assert_eq!(pos.offset, Some(-232));

        // Unknown position with positive offset: ?+10
        let (remaining, pos) = parse_cds_pos("?+10del").unwrap();
        assert_eq!(remaining, "del");
        assert_eq!(pos.base, CDS_BASE_UNKNOWN);
        assert!(pos.is_unknown());
        assert_eq!(pos.offset, Some(10));

        // Unknown offset: 148-?
        let (remaining, pos) = parse_cds_pos("148-?_228+?").unwrap();
        assert_eq!(remaining, "_228+?");
        assert_eq!(pos.base, 148);
        assert!(!pos.is_unknown());
        assert_eq!(pos.offset, Some(OFFSET_UNKNOWN_NEGATIVE));

        // Unknown offset: 228+?
        let (remaining, pos) = parse_cds_pos("228+?").unwrap();
        assert_eq!(remaining, "");
        assert_eq!(pos.base, 228);
        assert!(!pos.is_unknown());
        assert_eq!(pos.offset, Some(OFFSET_UNKNOWN_POSITIVE));
    }

    #[test]
    fn test_cds_pos_unknown_display() {
        // Simple unknown position
        let pos = CdsPos::unknown(None);
        assert_eq!(pos.to_string(), "?");

        // Unknown position with negative offset
        let pos = CdsPos::unknown(Some(-232));
        assert_eq!(pos.to_string(), "?-232");

        // Unknown position with positive offset
        let pos = CdsPos::unknown(Some(10));
        assert_eq!(pos.to_string(), "?+10");
    }

    #[test]
    fn test_parse_tx_pos_simple() {
        let (remaining, pos) = parse_tx_pos("100del").unwrap();
        assert_eq!(remaining, "del");
        assert_eq!(pos.base, 100);
        assert_eq!(pos.offset, None);
        assert!(!pos.downstream);
    }

    #[test]
    fn test_parse_tx_pos_upstream() {
        let (remaining, pos) = parse_tx_pos("-30C>G").unwrap();
        assert_eq!(remaining, "C>G");
        assert_eq!(pos.base, -30);
        assert!(!pos.downstream);
    }

    #[test]
    fn test_parse_tx_pos_downstream() {
        // n.*5 - downstream position
        let (remaining, pos) = parse_tx_pos("*5C>G").unwrap();
        assert_eq!(remaining, "C>G");
        assert_eq!(pos.base, 5);
        assert!(pos.downstream);
        assert_eq!(pos.to_string(), "*5");
    }

    #[test]
    fn test_parse_tx_pos_downstream_with_offset() {
        // n.*5+10 - downstream position with offset
        let (remaining, pos) = parse_tx_pos("*5+10C>G").unwrap();
        assert_eq!(remaining, "C>G");
        assert_eq!(pos.base, 5);
        assert_eq!(pos.offset, Some(10));
        assert!(pos.downstream);
        assert_eq!(pos.to_string(), "*5+10");
    }

    #[test]
    fn test_parse_tx_pos_with_offset() {
        let (remaining, pos) = parse_tx_pos("100+5G>A").unwrap();
        assert_eq!(remaining, "G>A");
        assert_eq!(pos.base, 100);
        assert_eq!(pos.offset, Some(5));
        assert!(!pos.downstream);
    }
}
