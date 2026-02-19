//! Input validation for web service security
//!
//! Provides comprehensive validation for HGVS strings, request parameters,
//! and other user inputs to prevent injection attacks and ensure data integrity.

use once_cell::sync::Lazy;
use regex::Regex;
use serde::{Deserialize, Serialize};

/// Maximum allowed length for HGVS variant strings
const MAX_HGVS_LENGTH: usize = 1000;

/// Maximum allowed timeout in seconds
const MAX_TIMEOUT_SECONDS: u32 = 300;

/// Minimum allowed timeout in seconds
const MIN_TIMEOUT_SECONDS: u32 = 1;

/// Maximum allowed variants in batch request
const MAX_BATCH_SIZE: usize = 1000;

/// Comprehensive HGVS pattern for basic format validation
/// Matches patterns like: NM_000001.1:c.123A>G, NC_000001.11:g.123del, etc.
/// Note: > is allowed as it's used in substitutions like A>G
/// The regex excludes dangerous shell characters but allows valid HGVS characters
static HGVS_PATTERN: Lazy<Regex> = Lazy::new(|| {
    // Use a simpler pattern that validates basic HGVS format
    // Dangerous characters are checked separately in validate_hgvs
    Regex::new(r"^[A-Z]{2}_\d+\.\d+:[cgmnpr]\..+$").unwrap()
});

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

    // Basic HGVS format validation using regex
    if !HGVS_PATTERN.is_match(input) {
        return Err(ValidationError::InvalidFormat);
    }

    // Security: check for dangerous characters that could be used in injection attacks
    // Note: > is allowed as it's used in HGVS substitutions like A>G
    if input.chars().any(|c| "<|&;`$(){}[]\\".contains(c)) {
        return Err(ValidationError::DangerousCharacters);
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
    }

    #[test]
    fn test_validate_hgvs_dangerous_characters() {
        // Test various dangerous characters that could be used in attacks
        // Note: > is allowed as it's used in HGVS substitutions
        let dangerous_inputs = vec![
            "NM_000001.1:c.123A<G", // < is not allowed
            "NM_000001.1:c.123del|rm -rf /",
            "NM_000001.1:c.123del$(whoami)",
            "NM_000001.1:c.123del{backdoor}",
            "NM_000001.1:c.123del[injection]",
            "NM_000001.1:c.123del\\evil",
        ];

        for input in dangerous_inputs {
            let result = validate_hgvs(input);
            assert!(matches!(result, Err(ValidationError::DangerousCharacters)));
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
