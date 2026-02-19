//! HGVS parser using nom
//!
//! This module provides a complete parser for HGVS variant nomenclature.
//!
//! # Error Handling Modes
//!
//! The parser supports three error handling modes:
//!
//! - **Strict** (default): Reject all non-standard input
//! - **Lenient**: Auto-correct common errors with warnings
//! - **Silent**: Auto-correct common errors without warnings
//!
//! Use [`parse_hgvs_with_config`] to specify the error handling mode.

pub mod accession;
pub mod edit;
pub mod fast_path;
pub mod position;
pub mod variant;

use crate::error::FerroError;
use crate::error_handling::{ErrorConfig, InputPreprocessor, ParseResultWithWarnings};
use crate::hgvs::HgvsVariant;

/// Parse an HGVS string into a variant
///
/// Uses strict error handling mode by default. For configurable error handling,
/// use [`parse_hgvs_with_config`] instead.
///
/// # Example
///
/// ```
/// use ferro_hgvs::parse_hgvs;
///
/// let variant = parse_hgvs("NM_000088.3:c.459del").unwrap();
/// println!("Parsed: {}", variant);
/// ```
pub fn parse_hgvs(input: &str) -> Result<HgvsVariant, FerroError> {
    variant::parse_variant(input)
}

/// Parse an HGVS string with fast-path optimization for common patterns
///
/// This function attempts to use specialized fast-path parsers for the most
/// common HGVS patterns (RefSeq, Ensembl, LRG, Assembly substitutions).
/// Falls back to the standard parser for complex or unusual patterns.
///
/// # Performance
///
/// For simple substitution patterns (which represent ~87% of ClinVar variants),
/// this provides **45-58% speedup**:
///
/// | Pattern Type | Example | Speedup |
/// |--------------|---------|---------|
/// | RefSeq genomic | `NC_000001.11:g.12345A>G` | ~50% |
/// | RefSeq coding | `NM_000088.3:c.459A>G` | ~57% |
/// | Ensembl | `ENST00000357033.8:c.100A>G` | ~54% |
/// | Assembly | `GRCh38(chr1):g.12345A>G` | ~47% |
///
/// # Tradeoffs
///
/// The fast-path adds a small overhead (~3-6%) for patterns it cannot optimize:
/// - Intronic variants (`c.100+5G>A`)
/// - UTR variants (`c.*100A>G`, `c.-50A>G`)
/// - Non-coding RNA (`n.100A>G`)
/// - RNA variants (`r.100a>g`)
/// - All deletions, insertions, duplications, etc.
///
/// # When to Use
///
/// **Recommended for:**
/// - Clinical variant databases (ClinVar, gnomAD)
/// - Batch processing of SNV-heavy datasets
/// - Performance-critical pipelines
///
/// **Use standard [`parse_hgvs`] instead when:**
/// - Data contains many complex variants (indels, intronic)
/// - Consistent performance across all variant types is needed
/// - Data composition is unknown
///
/// # Example
///
/// ```
/// use ferro_hgvs::parse_hgvs_fast;
///
/// // Fast path for RefSeq substitution (~2x faster)
/// let variant = parse_hgvs_fast("NC_000001.11:g.12345A>G").unwrap();
///
/// // Falls back to standard parser for complex patterns
/// let variant = parse_hgvs_fast("NM_000088.3:c.100+5A>G").unwrap();
/// ```
#[inline]
pub fn parse_hgvs_fast(input: &str) -> Result<HgvsVariant, FerroError> {
    let trimmed = input.trim();

    // Try fast path first
    match fast_path::try_fast_path(trimmed) {
        fast_path::FastPathResult::Success(variant) => Ok(variant),
        fast_path::FastPathResult::Fallback => variant::parse_variant(input),
    }
}

/// Parse an HGVS string with configurable error handling.
///
/// This function applies preprocessing based on the error configuration,
/// then parses the (potentially corrected) input.
///
/// # Example
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
/// use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode};
///
/// // Lenient mode: auto-correct common errors with warnings
/// let config = ErrorConfig::lenient();
/// let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
/// assert!(result.is_ok());
///
/// let parsed = result.unwrap();
/// assert!(parsed.had_corrections()); // Whitespace was trimmed
/// ```
pub fn parse_hgvs_with_config(
    input: &str,
    config: ErrorConfig,
) -> Result<ParseResultWithWarnings<HgvsVariant>, FerroError> {
    // Create preprocessor and preprocess input
    let preprocessor = InputPreprocessor::new(config);
    let preprocess_result = preprocessor.preprocess(input);

    // Check if preprocessing failed
    if !preprocess_result.success {
        return Err(preprocess_result
            .error
            .unwrap_or_else(|| FerroError::Parse {
                pos: 0,
                msg: "Preprocessing failed without error details".to_string(),
                diagnostic: None,
            }));
    }

    // Parse the preprocessed input
    let variant = variant::parse_variant(&preprocess_result.preprocessed)?;

    // Return result with warnings
    Ok(ParseResultWithWarnings::new(
        variant,
        preprocess_result.warnings,
        preprocess_result.original,
        preprocess_result.preprocessed,
    ))
}

/// Parse an HGVS string with lenient error handling.
///
/// This is a convenience function that uses lenient mode, which auto-corrects
/// common errors and returns warnings.
///
/// # Example
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // This will auto-correct the en-dash to hyphen
/// let result = parse_hgvs_lenient("NM_000088.3:c.100\u{2013}200del");
/// assert!(result.is_ok());
/// ```
pub fn parse_hgvs_lenient(input: &str) -> Result<ParseResultWithWarnings<HgvsVariant>, FerroError> {
    parse_hgvs_with_config(input, ErrorConfig::lenient())
}

/// Parse an HGVS string with silent error handling.
///
/// This is a convenience function that uses silent mode, which auto-corrects
/// common errors without generating warnings.
///
/// # Example
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_silent;
///
/// // This will silently auto-correct the en-dash to hyphen
/// let result = parse_hgvs_silent("NM_000088.3:c.100\u{2013}200del");
/// assert!(result.is_ok());
/// assert!(!result.unwrap().has_warnings()); // No warnings in silent mode
/// ```
pub fn parse_hgvs_silent(input: &str) -> Result<ParseResultWithWarnings<HgvsVariant>, FerroError> {
    parse_hgvs_with_config(input, ErrorConfig::silent())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_handling::{ErrorOverride, ErrorType};

    #[test]
    fn test_parse_simple_substitution() {
        let result = parse_hgvs("NC_000001.11:g.12345A>G");
        assert!(result.is_ok());
    }

    #[test]
    fn test_parse_deletion() {
        let result = parse_hgvs("NM_000088.3:c.459del");
        assert!(result.is_ok());
    }

    // Error handling mode tests
    #[test]
    fn test_parse_with_config_strict_valid() {
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("NM_000088.3:c.459del", config);
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(!parsed.had_corrections());
        assert!(!parsed.has_warnings());
    }

    #[test]
    fn test_parse_with_config_strict_rejects_en_dash() {
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("NM_000088.3:c.100\u{2013}200del", config);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_with_config_strict_rejects_whitespace() {
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_with_config_lenient_corrects_whitespace() {
        let config = ErrorConfig::lenient();
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert!(parsed.has_warnings());
    }

    #[test]
    fn test_parse_with_config_silent_no_warnings() {
        let config = ErrorConfig::silent();
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert!(!parsed.has_warnings()); // Silent mode = no warnings
    }

    #[test]
    fn test_parse_with_config_override() {
        // Lenient mode but override whitespace to reject
        let config =
            ErrorConfig::lenient().with_override(ErrorType::ExtraWhitespace, ErrorOverride::Reject);
        let result = parse_hgvs_with_config("  NM_000088.3:c.459del  ", config);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_lenient() {
        let result = parse_hgvs_lenient("  NM_000088.3:c.459del  ");
        assert!(result.is_ok());
        assert!(result.unwrap().had_corrections());
    }

    #[test]
    fn test_parse_silent() {
        let result = parse_hgvs_silent("  NM_000088.3:c.459del  ");
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert!(!parsed.has_warnings());
    }

    #[test]
    fn test_parse_lowercase_accession_lenient() {
        let result = parse_hgvs_lenient("nm_000088.3:c.459del");
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(parsed.had_corrections());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.459del");
    }
}
