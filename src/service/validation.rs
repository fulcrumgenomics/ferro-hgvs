//! Input validation for web service security
//!
//! Provides comprehensive validation for HGVS strings, request parameters,
//! and other user inputs to prevent injection attacks and ensure data integrity.

use serde::{Deserialize, Serialize};

use crate::hgvs::parser::parse_hgvs_lenient;

/// Maximum allowed length for HGVS variant strings
const MAX_HGVS_LENGTH: usize = 1000;

/// Maximum allowed timeout in seconds
const MAX_TIMEOUT_SECONDS: u32 = 300;

/// Minimum allowed timeout in seconds
const MIN_TIMEOUT_SECONDS: u32 = 1;

/// Maximum allowed variants in batch request
const MAX_BATCH_SIZE: usize = 1000;

/// Validation errors for user input
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ValidationError {
    /// Input string is empty
    Empty,
    /// Input string is too long
    TooLong { max: usize, actual: usize },
    /// Input contains non-ASCII characters
    NonAscii,
    /// Input doesn't match HGVS format pattern
    InvalidFormat,
    /// Input contains potentially dangerous characters
    DangerousCharacters,
    /// Timeout value is out of valid range
    InvalidTimeout { min: u32, max: u32, actual: u32 },
    /// Batch size exceeds maximum allowed
    BatchTooLarge { max: usize, actual: usize },
}

impl std::fmt::Display for ValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ValidationError::Empty => write!(f, "Input cannot be empty"),
            ValidationError::TooLong { max, actual } => {
                write!(f, "Input too long: {} characters (max: {})", actual, max)
            }
            ValidationError::NonAscii => write!(f, "Input must contain only ASCII characters"),
            ValidationError::InvalidFormat => write!(f, "Input does not match valid HGVS format"),
            ValidationError::DangerousCharacters => {
                write!(f, "Input contains potentially dangerous characters")
            }
            ValidationError::InvalidTimeout { min, max, actual } => {
                write!(
                    f,
                    "Timeout {} is out of valid range ({}-{})",
                    actual, min, max
                )
            }
            ValidationError::BatchTooLarge { max, actual } => {
                write!(f, "Batch size {} exceeds maximum allowed ({})", actual, max)
            }
        }
    }
}

impl std::error::Error for ValidationError {}

/// Validate an HGVS variant string for security and format compliance
pub fn validate_hgvs(input: &str) -> Result<(), ValidationError> {
    // Check for empty input
    if input.is_empty() {
        return Err(ValidationError::Empty);
    }

    // Check length limits to prevent memory exhaustion
    if input.len() > MAX_HGVS_LENGTH {
        return Err(ValidationError::TooLong {
            max: MAX_HGVS_LENGTH,
            actual: input.len(),
        });
    }

    // Security: ensure only ASCII characters (prevents encoding attacks)
    if !input.is_ascii() {
        return Err(ValidationError::NonAscii);
    }

    // Security: check for dangerous characters that could be used in injection attacks.
    // Note: > is allowed (HGVS substitutions like A>G).
    // [] and () are allowed (repeats like C[8], predicted effects like p.(Val600Glu),
    // uncertain positions, allele notation, compound references like NG_x(NM_y), and
    // assembly-prefixed accessions).
    // ; is allowed (allele phase separator like [var1;var2]).
    // This runs before the format gate so injection attempts surface as
    // `DangerousCharacters` rather than a generic format error.
    if input.chars().any(|c| "<|&`${}\\".contains(c)) {
        return Err(ValidationError::DangerousCharacters);
    }

    // Format validation: defer to the HGVS parser itself as the single source of truth
    // for "valid HGVS format", rather than a second, hand-rolled grammar approximation
    // that inevitably drifts and rejects valid inputs (e.g. the compound-reference form
    // `NG_x(NM_y):c.…`; see issue #1102). The parser is a pure, I/O-free, fuzz-hardened
    // function, and the length/ASCII/dangerous-character guards above have already
    // bounded and sanitized the input, so running it here adds no execution risk. The
    // lenient entry point is used so the gate is at least as permissive as the handlers,
    // which parse in lenient mode; the caller's chosen error mode still applies downstream.
    if parse_hgvs_lenient(input).is_err() {
        return Err(ValidationError::InvalidFormat);
    }

    Ok(())
}

/// Validate a list of HGVS variant strings for batch processing
pub fn validate_hgvs_batch(variants: &[String]) -> Result<(), ValidationError> {
    // Check batch size limits
    if variants.len() > MAX_BATCH_SIZE {
        return Err(ValidationError::BatchTooLarge {
            max: MAX_BATCH_SIZE,
            actual: variants.len(),
        });
    }

    // Validate each variant in the batch
    for variant in variants {
        validate_hgvs(variant)?;
    }

    Ok(())
}

/// Validate timeout value is within acceptable bounds
pub fn validate_timeout(timeout_seconds: u32) -> Result<(), ValidationError> {
    if !(MIN_TIMEOUT_SECONDS..=MAX_TIMEOUT_SECONDS).contains(&timeout_seconds) {
        return Err(ValidationError::InvalidTimeout {
            min: MIN_TIMEOUT_SECONDS,
            max: MAX_TIMEOUT_SECONDS,
            actual: timeout_seconds,
        });
    }
    Ok(())
}

/// Validate an optional timeout value
pub fn validate_optional_timeout(timeout_seconds: Option<u32>) -> Result<(), ValidationError> {
    if let Some(timeout) = timeout_seconds {
        validate_timeout(timeout)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_hgvs_valid_inputs() {
        // Valid HGVS strings should pass
        assert!(validate_hgvs("NM_000001.1:c.123A>G").is_ok());
        assert!(validate_hgvs("NC_000001.11:g.123456del").is_ok());
        assert!(validate_hgvs("NP_000001.1:p.Arg123Cys").is_ok());
        assert!(validate_hgvs("NM_000001.1:c.123+1G>A").is_ok());
        assert!(validate_hgvs("NM_000001.1:c.123-10_123-5del").is_ok());
    }

    #[test]
    fn test_validate_hgvs_empty_input() {
        let result = validate_hgvs("");
        assert!(matches!(result, Err(ValidationError::Empty)));
    }

    #[test]
    fn test_validate_hgvs_too_long() {
        let long_input = "A".repeat(MAX_HGVS_LENGTH + 1);
        let result = validate_hgvs(&long_input);
        assert!(matches!(result, Err(ValidationError::TooLong { .. })));
    }

    #[test]
    fn test_validate_hgvs_non_ascii() {
        let result = validate_hgvs("NM_000001.1:c.123A>🦀");
        assert!(matches!(result, Err(ValidationError::NonAscii)));
    }

    #[test]
    fn test_validate_hgvs_invalid_format() {
        // Not HGVS format
        assert!(matches!(
            validate_hgvs("not-hgvs-format"),
            Err(ValidationError::InvalidFormat)
        ));

        // Missing accession
        assert!(matches!(
            validate_hgvs(":c.123A>G"),
            Err(ValidationError::InvalidFormat)
        ));

        // Invalid variant type
        assert!(matches!(
            validate_hgvs("NM_000001.1:x.123A>G"),
            Err(ValidationError::InvalidFormat)
        ));

        // Trailing garbage after an otherwise-valid variant. The parser is
        // anchored — it must consume the WHOLE string — so a valid prefix
        // followed by junk is rejected (the security-relevant property this
        // gate relies on now that the regex is gone).
        assert!(matches!(
            validate_hgvs("NM_000001.1:c.123A>G junk"),
            Err(ValidationError::InvalidFormat)
        ));
        assert!(matches!(
            validate_hgvs("NM_000001.1:c.123A>Gxyzzy"),
            Err(ValidationError::InvalidFormat)
        ));
    }

    #[test]
    fn test_validate_hgvs_accepts_surrounding_whitespace() {
        // The parser-based gate trims surrounding whitespace (lenient mode),
        // matching how the downstream `ferro` tool parses the same input, so a
        // padded-but-valid variant is accepted. The old regex gate (anchored
        // `^…$`) rejected it — pin the intended behavior of the swap.
        assert!(validate_hgvs("  NM_000001.1:c.123A>G  ").is_ok());
    }

    #[test]
    fn test_validate_hgvs_compound_reference_selector() {
        // Regression (issue #1102): a gene/genomic reference carrying a parenthesized
        // transcript or gene selector — the HGVS "reference with coordinates from an
        // aligned transcript" form — is valid HGVS the parser accepts, and must pass
        // the security gate. The web portal's spec example is the first case here.
        assert!(validate_hgvs("NG_012232.1(NM_004006.2):c.93+1G>T").is_ok());
        assert!(validate_hgvs("NC_000001.10(NM_002074.3):c.58del").is_ok());
        assert!(validate_hgvs("ENSG00000184937.16(ENST00000452863.10):c.10del").is_ok());
        assert!(validate_hgvs("NG_007485.1(CDKN2A):n.204_205insATC").is_ok());
    }

    #[test]
    fn test_validate_hgvs_lrg_and_gene_symbol_references() {
        // LRG transcript/protein references (t1/p1 suffix) and bare gene/protein-symbol
        // references are also valid HGVS the parser accepts; the gate must not reject them.
        assert!(validate_hgvs("LRG_199t1:c.100A>G").is_ok());
        assert!(validate_hgvs("LRG_199p1:p.Val100Glu").is_ok());
        assert!(validate_hgvs("BRAF:p.Val600Glu").is_ok());
    }

    #[test]
    fn test_validate_hgvs_repeats() {
        // Repeat notation uses square brackets — must be accepted
        assert!(validate_hgvs("NM_003820.4:c.495_500C[8]").is_ok());
        assert!(validate_hgvs("NC_000001.11:g.2560658_2560663C[8]").is_ok());
        // Multi-repeat
        assert!(validate_hgvs("NM_000001.1:c.100_105CTG[9]TTG[1]CTG[13]").is_ok());
        // Uncertain repeat count
        assert!(validate_hgvs("NM_000001.1:c.100_105CAG[10_15]").is_ok());
    }

    #[test]
    fn test_validate_hgvs_predicted_effects() {
        // Predicted protein effects use parentheses — must be accepted
        assert!(validate_hgvs("NP_000001.1:p.(Val600Glu)").is_ok());
        assert!(validate_hgvs("NP_000001.1:p.(=)").is_ok());
        assert!(validate_hgvs("NP_000001.1:p.(?)").is_ok());
    }

    #[test]
    fn test_validate_hgvs_uncertain_positions() {
        // Uncertain positions use parentheses — must be accepted
        assert!(validate_hgvs("NM_000001.1:c.(100_200)del").is_ok());
        assert!(validate_hgvs("NC_000001.11:g.(?_100)_(200_?)del").is_ok());
    }

    #[test]
    fn test_validate_hgvs_allele_notation() {
        // Allele notation uses square brackets and semicolons
        assert!(validate_hgvs("NM_000001.1:c.[123A>G;456C>T]").is_ok());
        assert!(validate_hgvs("NM_000001.1:c.[123A>G];[456C>T]").is_ok());
    }

    #[test]
    fn test_validate_hgvs_assembly_prefixed_accessions() {
        // Assembly-prefixed genomic accessions (GRCh and hg short names)
        assert!(validate_hgvs("GRCh38(chr1):g.2560658_2560663C[8]").is_ok());
        assert!(validate_hgvs("GRCh37(chr1):g.12345A>G").is_ok());
        assert!(validate_hgvs("GRCh36(chr1):g.100del").is_ok());
        assert!(validate_hgvs("hg38(chr1):g.12345A>G").is_ok());
        assert!(validate_hgvs("hg19(chr1):g.12345A>G").is_ok());
        assert!(validate_hgvs("hg18(chr1):g.12345A>G").is_ok());
    }

    #[test]
    fn test_validate_hgvs_ensembl_accessions() {
        assert!(validate_hgvs("ENST00000123456.1:c.123A>G").is_ok());
        assert!(validate_hgvs("ENSG00000123456.1:g.100A>G").is_ok());
        assert!(validate_hgvs("ENSP00000123456.1:p.Val600Glu").is_ok());
        assert!(validate_hgvs("ENSE00000123456.1:g.100A>G").is_ok());
        assert!(validate_hgvs("ENSR00000123456.1:g.100A>G").is_ok());
    }

    #[test]
    fn test_validate_hgvs_lrg_accessions() {
        assert!(validate_hgvs("LRG_1:g.12345A>G").is_ok());
    }

    #[test]
    fn test_validate_hgvs_genbank_accessions() {
        assert!(validate_hgvs("U12345.1:g.100A>G").is_ok());
    }

    #[test]
    fn test_validate_hgvs_dangerous_characters() {
        // Test various dangerous characters that could be used in attacks
        // Note: > is allowed as it's used in HGVS substitutions
        let dangerous_inputs = vec![
            "NM_000001.1:c.123A<G",
            "NM_000001.1:c.123del|rm -rf /",
            "NM_000001.1:c.123del$(whoami)",
            "NM_000001.1:c.123del{backdoor}",
            "NM_000001.1:c.123del\\evil",
            "NM_000001.1:c.123del`whoami`",
            "NM_000001.1:c.123del&rm",
        ];

        for input in dangerous_inputs {
            let result = validate_hgvs(input);
            assert!(
                matches!(result, Err(ValidationError::DangerousCharacters)),
                "Expected DangerousCharacters for input: {input}"
            );
        }
    }

    #[test]
    fn test_validate_hgvs_batch() {
        // Valid batch should pass
        let valid_batch = vec![
            "NM_000001.1:c.123A>G".to_string(),
            "NM_000002.2:c.456C>T".to_string(),
        ];
        assert!(validate_hgvs_batch(&valid_batch).is_ok());

        // Batch with invalid variant should fail
        let invalid_batch = vec![
            "NM_000001.1:c.123A>G".to_string(),
            "invalid-variant".to_string(),
        ];
        assert!(validate_hgvs_batch(&invalid_batch).is_err());
    }

    #[test]
    fn test_validate_hgvs_batch_too_large() {
        // Create batch larger than allowed
        let large_batch = vec!["NM_000001.1:c.123A>G".to_string(); MAX_BATCH_SIZE + 1];
        let result = validate_hgvs_batch(&large_batch);
        assert!(matches!(result, Err(ValidationError::BatchTooLarge { .. })));
    }

    #[test]
    fn test_validate_timeout() {
        // Valid timeouts should pass
        assert!(validate_timeout(30).is_ok());
        assert!(validate_timeout(MIN_TIMEOUT_SECONDS).is_ok());
        assert!(validate_timeout(MAX_TIMEOUT_SECONDS).is_ok());

        // Invalid timeouts should fail
        assert!(matches!(
            validate_timeout(0),
            Err(ValidationError::InvalidTimeout { .. })
        ));
        assert!(matches!(
            validate_timeout(MAX_TIMEOUT_SECONDS + 1),
            Err(ValidationError::InvalidTimeout { .. })
        ));
    }

    #[test]
    fn test_validate_optional_timeout() {
        // None should pass
        assert!(validate_optional_timeout(None).is_ok());

        // Valid timeout should pass
        assert!(validate_optional_timeout(Some(30)).is_ok());

        // Invalid timeout should fail
        assert!(validate_optional_timeout(Some(0)).is_err());
    }
}
