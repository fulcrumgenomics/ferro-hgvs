//! IVS (Intervening Sequence) notation parser.
//!
//! IVS notation is a deprecated way to describe intronic variants:
//! - `c.IVS4+1G>A` means 1 base into intron 4 (from the 5' donor side)
//! - `c.IVS4-2A>G` means 2 bases before exon 5 (from the 3' acceptor side)

use super::ExonData;
use crate::error::FerroError;

/// Parsed IVS position.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IvsPosition {
    /// Intron number (1-based).
    pub intron: u32,
    /// Offset from exon boundary (positive = donor side, negative = acceptor side).
    pub offset: i64,
}

/// Parse IVS notation from a string.
///
/// Accepts formats like:
/// - `IVS4+1` (1 base into intron 4)
/// - `IVS4-2` (2 bases before exon 5)
/// - `IVS12+100` (100 bases into intron 12)
pub fn parse_ivs(input: &str) -> Result<IvsPosition, FerroError> {
    let input = input.trim();

    // Must start with IVS
    if !input.starts_with("IVS") {
        return Err(FerroError::Parse {
            pos: 0,
            msg: format!("Expected IVS prefix: {}", input),
            diagnostic: None,
        });
    }

    let rest = &input[3..];

    // Find the sign (+ or -)
    let sign_pos = rest.find('+').or_else(|| rest.find('-'));

    match sign_pos {
        Some(pos) => {
            // Parse intron number
            let intron_str = &rest[..pos];
            let intron: u32 = intron_str.parse().map_err(|_| FerroError::Parse {
                pos: 3,
                msg: format!("Invalid intron number: {}", intron_str),
                diagnostic: None,
            })?;

            // Parse offset (including sign)
            let offset_str = &rest[pos..];
            let offset: i64 = offset_str.parse().map_err(|_| FerroError::Parse {
                pos: 3 + pos,
                msg: format!("Invalid offset: {}", offset_str),
                diagnostic: None,
            })?;

            if offset == 0 {
                return Err(FerroError::InvalidCoordinates {
                    msg: "IVS offset cannot be 0".to_string(),
                });
            }

            Ok(IvsPosition { intron, offset })
        }
        None => Err(FerroError::Parse {
            pos: 0,
            msg: format!("IVS notation requires +/- offset: {}", input),
            diagnostic: None,
        }),
    }
}

/// Convert IVS notation in a variant string to modern format.
///
/// Takes the full variant string and exon data, returns the converted string.
///
/// # Example
///
/// ```ignore
/// let exon_data = ExonData::new(vec![(1, 100), (201, 300), (401, 500), (601, 700)]);
/// let converted = convert_ivs_notation("NM_000088.3:c.IVS3+1G>A", &exon_data)?;
/// // Returns "NM_000088.3:c.500+1G>A"
/// ```
pub fn convert_ivs_notation(
    input: &str,
    exon_data: &ExonData,
) -> Result<Option<String>, FerroError> {
    // Find IVS in the input
    let ivs_start = match input.find("IVS") {
        Some(pos) => pos,
        None => return Ok(None),
    };

    // Find the end of the IVS notation (look for the nucleotide after offset)
    let after_ivs = &input[ivs_start..];

    // Parse to find where IVS notation ends
    // Pattern: IVS<num>+/-<offset><rest>
    let mut end_pos = 3; // Skip "IVS"

    // Skip digits (intron number)
    while end_pos < after_ivs.len()
        && after_ivs
            .as_bytes()
            .get(end_pos)
            .is_some_and(|&b| b.is_ascii_digit())
    {
        end_pos += 1;
    }

    // Must have + or -
    if end_pos >= after_ivs.len() {
        return Ok(None);
    }

    let sign_char = after_ivs.as_bytes()[end_pos];
    if sign_char != b'+' && sign_char != b'-' {
        return Ok(None);
    }
    end_pos += 1;

    // Skip offset digits
    while end_pos < after_ivs.len()
        && after_ivs
            .as_bytes()
            .get(end_pos)
            .is_some_and(|&b| b.is_ascii_digit())
    {
        end_pos += 1;
    }

    // Extract the IVS part
    let ivs_str = &after_ivs[..end_pos];
    let ivs = parse_ivs(ivs_str)?;

    // Convert to CDS position
    let (base, offset) = exon_data.ivs_to_cds(&ivs)?;

    // Build the modern position string
    let modern_pos = match offset {
        Some(off) if off > 0 => format!("{}+{}", base, off),
        Some(off) => format!("{}{}", base, off), // Negative already has sign
        None => format!("{}", base),
    };

    // Reconstruct the variant string
    let prefix = &input[..ivs_start];
    let suffix = &after_ivs[end_pos..];

    Ok(Some(format!("{}{}{}", prefix, modern_pos, suffix)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_ivs_donor() {
        let ivs = parse_ivs("IVS4+1").unwrap();
        assert_eq!(ivs.intron, 4);
        assert_eq!(ivs.offset, 1);
    }

    #[test]
    fn test_parse_ivs_acceptor() {
        let ivs = parse_ivs("IVS4-2").unwrap();
        assert_eq!(ivs.intron, 4);
        assert_eq!(ivs.offset, -2);
    }

    #[test]
    fn test_parse_ivs_large_offset() {
        let ivs = parse_ivs("IVS12+100").unwrap();
        assert_eq!(ivs.intron, 12);
        assert_eq!(ivs.offset, 100);
    }

    #[test]
    fn test_parse_ivs_negative_large() {
        let ivs = parse_ivs("IVS1-500").unwrap();
        assert_eq!(ivs.intron, 1);
        assert_eq!(ivs.offset, -500);
    }

    #[test]
    fn test_parse_ivs_invalid_no_offset() {
        assert!(parse_ivs("IVS4").is_err());
    }

    #[test]
    fn test_parse_ivs_invalid_zero_offset() {
        assert!(parse_ivs("IVS4+0").is_err());
    }

    #[test]
    fn test_parse_ivs_invalid_no_ivs() {
        assert!(parse_ivs("c.100+5").is_err());
    }

    #[test]
    fn test_convert_ivs_notation() {
        // Exons: 1-100, 201-300, 401-500, 601-700
        let exon_data = ExonData::new(vec![(1, 100), (201, 300), (401, 500), (601, 700)]);

        // IVS1+5: 5 bases after exon 1 (which ends at 100)
        let result = convert_ivs_notation("NM_000088.3:c.IVS1+5G>A", &exon_data).unwrap();
        assert_eq!(result, Some("NM_000088.3:c.100+5G>A".to_string()));

        // IVS2-10: 10 bases before exon 3 (which starts at 401)
        let result = convert_ivs_notation("NM_000088.3:c.IVS2-10A>G", &exon_data).unwrap();
        assert_eq!(result, Some("NM_000088.3:c.401-10A>G".to_string()));
    }

    #[test]
    fn test_convert_ivs_notation_interval() {
        // Test interval like IVS7+5_IVS7+10del
        let exon_data = ExonData::new(vec![
            (1, 100),
            (201, 300),
            (401, 500),
            (601, 700),
            (801, 900),
            (1001, 1100),
            (1201, 1300),
        ]);

        // IVS7+5 (after exon 7 which ends at 1300)
        let result = convert_ivs_notation("c.IVS7+5del", &exon_data).unwrap();
        assert_eq!(result, Some("c.1300+5del".to_string()));
    }

    #[test]
    fn test_convert_ivs_no_ivs() {
        let exon_data = ExonData::new(vec![(1, 100)]);
        let result = convert_ivs_notation("NM_000088.3:c.100+5G>A", &exon_data).unwrap();
        assert_eq!(result, None);
    }
}
