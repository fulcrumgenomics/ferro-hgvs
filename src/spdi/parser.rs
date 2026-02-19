//! SPDI parser.
//!
//! Parses SPDI format strings into [`SpdiVariant`] instances.

use super::SpdiVariant;
use std::fmt;

/// Error type for SPDI parsing failures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpdiParseError {
    /// Not enough colons in the input (expected exactly 3).
    NotEnoughParts {
        /// Number of parts found.
        found: usize,
    },
    /// Too many colons in the input (expected exactly 3).
    TooManyParts {
        /// Number of parts found.
        found: usize,
    },
    /// Empty sequence identifier.
    EmptySequence,
    /// Invalid position (not a valid unsigned integer).
    InvalidPosition {
        /// The invalid position string.
        value: String,
    },
    /// Position is out of range.
    PositionOutOfRange {
        /// The position value that was too large.
        value: String,
    },
}

impl fmt::Display for SpdiParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SpdiParseError::NotEnoughParts { found } => {
                write!(
                    f,
                    "invalid SPDI format: expected 4 colon-separated parts, found {}",
                    found
                )
            }
            SpdiParseError::TooManyParts { found } => {
                write!(
                    f,
                    "invalid SPDI format: expected 4 colon-separated parts, found {}",
                    found
                )
            }
            SpdiParseError::EmptySequence => {
                write!(
                    f,
                    "invalid SPDI format: sequence identifier cannot be empty"
                )
            }
            SpdiParseError::InvalidPosition { value } => {
                write!(
                    f,
                    "invalid SPDI format: position '{}' is not a valid number",
                    value
                )
            }
            SpdiParseError::PositionOutOfRange { value } => {
                write!(
                    f,
                    "invalid SPDI format: position '{}' is out of range",
                    value
                )
            }
        }
    }
}

impl std::error::Error for SpdiParseError {}

/// Parse an SPDI format string.
///
/// The SPDI format is: `sequence:position:deletion:insertion`
///
/// # Arguments
///
/// * `input` - The SPDI string to parse
///
/// # Returns
///
/// * `Ok(SpdiVariant)` - Successfully parsed variant
/// * `Err(SpdiParseError)` - Parse error with details
///
/// # Examples
///
/// ```
/// use ferro_hgvs::spdi::parse_spdi;
///
/// // Substitution
/// let spdi = parse_spdi("NC_000001.11:12344:A:G").unwrap();
/// assert_eq!(spdi.sequence, "NC_000001.11");
/// assert_eq!(spdi.position, 12344);
///
/// // Deletion
/// let spdi = parse_spdi("NC_000001.11:99:ATG:").unwrap();
/// assert!(spdi.is_deletion());
///
/// // Insertion
/// let spdi = parse_spdi("NC_000001.11:100::ATG").unwrap();
/// assert!(spdi.is_insertion());
/// ```
pub fn parse_spdi(input: &str) -> Result<SpdiVariant, SpdiParseError> {
    let input = input.trim();

    // Split on colons - we expect exactly 4 parts
    // Note: SPDI allows empty deletion/insertion, so we need to handle that
    let parts: Vec<&str> = input.splitn(4, ':').collect();

    if parts.len() < 4 {
        // Check if we have exactly 3 colons by counting
        let colon_count = input.chars().filter(|&c| c == ':').count();
        if colon_count < 3 {
            return Err(SpdiParseError::NotEnoughParts {
                found: colon_count + 1,
            });
        }
    }

    // Now split properly - splitn(4, ':') should give us exactly 4 parts
    // if there are at least 3 colons
    let parts: Vec<&str> = input.splitn(4, ':').collect();

    if parts.len() != 4 {
        return Err(SpdiParseError::NotEnoughParts { found: parts.len() });
    }

    // Check for too many colons in the last part (which would indicate >4 parts)
    // The sequence and position should not contain colons
    // Deletion and insertion might contain colons in theory (though unusual)
    // For now, we'll be strict about the format

    let sequence = parts[0];
    let position_str = parts[1];
    let deletion = parts[2];
    let insertion = parts[3];

    // Validate sequence
    if sequence.is_empty() {
        return Err(SpdiParseError::EmptySequence);
    }

    // Parse position
    let position: u64 = position_str
        .parse()
        .map_err(|_| SpdiParseError::InvalidPosition {
            value: position_str.to_string(),
        })?;

    Ok(SpdiVariant {
        sequence: sequence.to_string(),
        position,
        deletion: deletion.to_string(),
        insertion: insertion.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_substitution() {
        let spdi = parse_spdi("NC_000001.11:12344:A:G").unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 12344);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_parse_deletion() {
        let spdi = parse_spdi("NC_000001.11:99:ATG:").unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_parse_insertion() {
        let spdi = parse_spdi("NC_000001.11:100::ATG").unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 100);
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_parse_delins() {
        let spdi = parse_spdi("NC_000001.11:99:ATG:TTCC").unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");
    }

    #[test]
    fn test_parse_identity() {
        let spdi = parse_spdi("NC_000001.11:100:A:A").unwrap();
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "A");
    }

    #[test]
    fn test_parse_with_whitespace() {
        let spdi = parse_spdi("  NC_000001.11:12344:A:G  ").unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
    }

    #[test]
    fn test_parse_empty_both() {
        // Both empty is valid (represents identity at a position)
        let spdi = parse_spdi("NC_000001.11:100::").unwrap();
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_parse_ensembl_sequence() {
        let spdi = parse_spdi("ENST00000357654.9:100:A:G").unwrap();
        assert_eq!(spdi.sequence, "ENST00000357654.9");
    }

    #[test]
    fn test_parse_position_zero() {
        let spdi = parse_spdi("NC_000001.11:0:A:G").unwrap();
        assert_eq!(spdi.position, 0);
    }

    #[test]
    fn test_parse_large_position() {
        let spdi = parse_spdi("NC_000001.11:249000000:A:G").unwrap();
        assert_eq!(spdi.position, 249_000_000);
    }

    #[test]
    fn test_parse_long_sequence() {
        let long_del = "A".repeat(100);
        let input = format!("NC_000001.11:100:{}:G", long_del);
        let spdi = parse_spdi(&input).unwrap();
        assert_eq!(spdi.deletion.len(), 100);
    }

    // Error cases

    #[test]
    fn test_parse_not_enough_parts() {
        let result = parse_spdi("NC_000001.11:12344:A");
        assert!(matches!(
            result,
            Err(SpdiParseError::NotEnoughParts { found: 3 })
        ));
    }

    #[test]
    fn test_parse_too_few_colons() {
        let result = parse_spdi("NC_000001.11:12344");
        assert!(matches!(result, Err(SpdiParseError::NotEnoughParts { .. })));
    }

    #[test]
    fn test_parse_empty_sequence() {
        let result = parse_spdi(":12344:A:G");
        assert!(matches!(result, Err(SpdiParseError::EmptySequence)));
    }

    #[test]
    fn test_parse_invalid_position() {
        let result = parse_spdi("NC_000001.11:abc:A:G");
        assert!(matches!(
            result,
            Err(SpdiParseError::InvalidPosition { .. })
        ));
    }

    #[test]
    fn test_parse_negative_position() {
        let result = parse_spdi("NC_000001.11:-1:A:G");
        assert!(matches!(
            result,
            Err(SpdiParseError::InvalidPosition { .. })
        ));
    }

    #[test]
    fn test_parse_float_position() {
        let result = parse_spdi("NC_000001.11:12.5:A:G");
        assert!(matches!(
            result,
            Err(SpdiParseError::InvalidPosition { .. })
        ));
    }

    #[test]
    fn test_parse_empty_string() {
        let result = parse_spdi("");
        assert!(matches!(result, Err(SpdiParseError::NotEnoughParts { .. })));
    }

    #[test]
    fn test_error_display() {
        let err = SpdiParseError::NotEnoughParts { found: 2 };
        assert!(err.to_string().contains("expected 4"));
        assert!(err.to_string().contains("found 2"));

        let err = SpdiParseError::EmptySequence;
        assert!(err.to_string().contains("empty"));

        let err = SpdiParseError::InvalidPosition {
            value: "abc".to_string(),
        };
        assert!(err.to_string().contains("abc"));
    }

    #[test]
    fn test_roundtrip() {
        let input = "NC_000001.11:12344:A:G";
        let spdi = parse_spdi(input).unwrap();
        assert_eq!(spdi.to_string(), input);
    }

    #[test]
    fn test_roundtrip_deletion() {
        let input = "NC_000001.11:99:ATG:";
        let spdi = parse_spdi(input).unwrap();
        assert_eq!(spdi.to_string(), input);
    }

    #[test]
    fn test_roundtrip_insertion() {
        let input = "NC_000001.11:100::ATG";
        let spdi = parse_spdi(input).unwrap();
        assert_eq!(spdi.to_string(), input);
    }
}
