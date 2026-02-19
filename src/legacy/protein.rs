//! Legacy protein notation parser.
//!
//! Handles deprecated protein notation formats:
//! - Arrow notation: `V600>E` → `Val600Glu`
//! - Number-first: `600V>E` → `Val600Glu`

use super::LegacyFormat;
use crate::hgvs::location::AminoAcid;

/// Legacy protein format detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LegacyProteinFormat {
    /// Arrow between amino acids: `V600>E`
    ArrowSubstitution,
    /// Number before amino acids: `600V>E`
    NumberFirst,
}

/// Parse a legacy protein notation and convert to modern format.
///
/// Returns the legacy format detected and the converted string.
#[allow(dead_code)]
pub fn parse_legacy_protein(input: &str) -> Option<(LegacyProteinFormat, String)> {
    convert_legacy_protein(input).map(|(format, s)| {
        let prot_format = match format {
            super::LegacyFormat::ArrowProteinSubstitution => LegacyProteinFormat::ArrowSubstitution,
            super::LegacyFormat::NumberFirstProtein => LegacyProteinFormat::NumberFirst,
            _ => LegacyProteinFormat::ArrowSubstitution, // Fallback
        };
        (prot_format, s)
    })
}

/// Convert legacy protein notation to modern format.
///
/// Handles:
/// - `V600>E` → `Val600Glu` (arrow notation)
/// - `600V>E` → `Val600Glu` (number-first)
/// - `V600>*` → `Val600Ter` (stop codon with arrow)
pub fn convert_legacy_protein(input: &str) -> Option<(LegacyFormat, String)> {
    // Try number-first format: 600V>E
    if let Some(result) = try_number_first(input) {
        return Some((LegacyFormat::NumberFirstProtein, result));
    }

    // Try arrow format: V600>E
    if let Some(result) = try_arrow_format(input) {
        return Some((LegacyFormat::ArrowProteinSubstitution, result));
    }

    None
}

/// Try to parse number-first format: 600V>E
fn try_number_first(input: &str) -> Option<String> {
    // Must start with digit
    if !input.starts_with(|c: char| c.is_ascii_digit()) {
        return None;
    }

    // Find where number ends
    let num_end = input.find(|c: char| !c.is_ascii_digit())?;
    let position: u64 = input[..num_end].parse().ok()?;

    let rest = &input[num_end..];

    // Must have at least 2 chars for ref AA and something else
    if rest.len() < 2 {
        return None;
    }

    // Check for arrow
    if !rest.contains('>') {
        return None;
    }

    // Parse ref amino acid (1-letter)
    let ref_char = rest.chars().next()?;
    let ref_aa = AminoAcid::from_one_letter(ref_char)?;

    // Find arrow
    let arrow_pos = rest.find('>')?;

    // Parse alt amino acid
    let alt_part = &rest[arrow_pos + 1..];
    let (alt_aa, suffix) = parse_amino_acid_with_suffix(alt_part)?;

    // Build modern notation
    Some(format!("{}{}{}{}", ref_aa, position, alt_aa, suffix))
}

/// Try to parse arrow format: V600>E
fn try_arrow_format(input: &str) -> Option<String> {
    // Must start with letter (amino acid)
    if !input.starts_with(|c: char| c.is_ascii_alphabetic()) {
        // Also check for * (Ter/stop)
        if !input.starts_with('*') {
            return None;
        }
    }

    // Must contain arrow
    if !input.contains('>') {
        return None;
    }

    // Try to parse as: AA<number>>AA
    let mut chars = input.chars().peekable();

    // Parse ref amino acid (1 or 3 letter)
    let ref_aa = if let Some(aa) = try_parse_three_letter(input) {
        // Consume 3 chars
        for _ in 0..3 {
            chars.next();
        }
        aa
    } else {
        let c = chars.next()?;
        AminoAcid::from_one_letter(c)?
    };

    // Collect position number
    let mut position = String::new();
    while let Some(&c) = chars.peek() {
        if c.is_ascii_digit() {
            position.push(chars.next().unwrap());
        } else {
            break;
        }
    }

    if position.is_empty() {
        return None;
    }

    // Check for arrow
    if chars.next() != Some('>') {
        return None;
    }

    // Parse alt amino acid
    let remaining: String = chars.collect();
    let (alt_aa, suffix) = parse_amino_acid_with_suffix(&remaining)?;

    // Build modern notation
    Some(format!("{}{}{}{}", ref_aa, position, alt_aa, suffix))
}

/// Try to parse a 3-letter amino acid at the start of input.
fn try_parse_three_letter(input: &str) -> Option<AminoAcid> {
    if input.len() < 3 {
        return None;
    }

    let codes = [
        ("Ala", AminoAcid::Ala),
        ("Arg", AminoAcid::Arg),
        ("Asn", AminoAcid::Asn),
        ("Asp", AminoAcid::Asp),
        ("Cys", AminoAcid::Cys),
        ("Gln", AminoAcid::Gln),
        ("Glu", AminoAcid::Glu),
        ("Gly", AminoAcid::Gly),
        ("His", AminoAcid::His),
        ("Ile", AminoAcid::Ile),
        ("Leu", AminoAcid::Leu),
        ("Lys", AminoAcid::Lys),
        ("Met", AminoAcid::Met),
        ("Phe", AminoAcid::Phe),
        ("Pro", AminoAcid::Pro),
        ("Sec", AminoAcid::Sec),
        ("Ser", AminoAcid::Ser),
        ("Thr", AminoAcid::Thr),
        ("Trp", AminoAcid::Trp),
        ("Tyr", AminoAcid::Tyr),
        ("Val", AminoAcid::Val),
        ("Ter", AminoAcid::Ter),
        ("Xaa", AminoAcid::Xaa),
    ];

    let prefix = &input[..3];
    // Case-insensitive comparison
    let normalized: String = prefix
        .chars()
        .enumerate()
        .map(|(i, c)| {
            if i == 0 {
                c.to_ascii_uppercase()
            } else {
                c.to_ascii_lowercase()
            }
        })
        .collect();

    for (code, aa) in codes {
        if normalized == code {
            return Some(aa);
        }
    }

    None
}

/// Parse amino acid with potential suffix (like trailing characters after the AA).
fn parse_amino_acid_with_suffix(input: &str) -> Option<(AminoAcid, &str)> {
    // Try 3-letter first
    if let Some(aa) = try_parse_three_letter(input) {
        return Some((aa, &input[3..]));
    }

    // Try 1-letter
    if input.is_empty() {
        return None;
    }

    let c = input.chars().next()?;
    let aa = AminoAcid::from_one_letter(c)?;
    Some((aa, &input[1..]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arrow_substitution() {
        let (format, converted) = convert_legacy_protein("V600>E").unwrap();
        assert_eq!(format, LegacyFormat::ArrowProteinSubstitution);
        assert_eq!(converted, "Val600Glu");
    }

    #[test]
    fn test_arrow_substitution_single_letter() {
        let (format, converted) = convert_legacy_protein("M1>I").unwrap();
        assert_eq!(format, LegacyFormat::ArrowProteinSubstitution);
        assert_eq!(converted, "Met1Ile");
    }

    #[test]
    fn test_arrow_to_stop() {
        let (format, converted) = convert_legacy_protein("W288>*").unwrap();
        assert_eq!(format, LegacyFormat::ArrowProteinSubstitution);
        assert_eq!(converted, "Trp288Ter");
    }

    #[test]
    fn test_number_first() {
        let (format, converted) = convert_legacy_protein("600V>E").unwrap();
        assert_eq!(format, LegacyFormat::NumberFirstProtein);
        assert_eq!(converted, "Val600Glu");
    }

    #[test]
    fn test_number_first_to_stop() {
        let (format, converted) = convert_legacy_protein("288W>*").unwrap();
        assert_eq!(format, LegacyFormat::NumberFirstProtein);
        assert_eq!(converted, "Trp288Ter");
    }

    #[test]
    fn test_three_letter_arrow() {
        let (format, converted) = convert_legacy_protein("Val600>Glu").unwrap();
        assert_eq!(format, LegacyFormat::ArrowProteinSubstitution);
        assert_eq!(converted, "Val600Glu");
    }

    #[test]
    fn test_no_arrow_returns_none() {
        assert!(convert_legacy_protein("V600E").is_none());
    }

    #[test]
    fn test_normal_hgvs_returns_none() {
        // This is already valid modern HGVS, should return None
        assert!(convert_legacy_protein("Val600Glu").is_none());
    }

    #[test]
    fn test_try_parse_three_letter() {
        assert_eq!(try_parse_three_letter("Val600"), Some(AminoAcid::Val));
        assert_eq!(try_parse_three_letter("val600"), Some(AminoAcid::Val));
        assert_eq!(try_parse_three_letter("VAL600"), Some(AminoAcid::Val));
        assert_eq!(try_parse_three_letter("Ter"), Some(AminoAcid::Ter));
        assert_eq!(try_parse_three_letter("Xy"), None); // Too short
        assert_eq!(try_parse_three_letter("Xyz"), None); // Invalid
    }
}
