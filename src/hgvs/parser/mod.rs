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
use crate::error_handling::{
    CorrectionWarning, ErrorConfig, ErrorType, InputPreprocessor, ParseResultWithWarnings,
    ResolvedAction,
};
use crate::hgvs::HgvsVariant;

/// Parse an HGVS string into a variant
///
/// Uses strict error handling mode by default. For configurable error handling,
/// use [`parse_hgvs_with_config`] instead.
///
/// # Performance
///
/// Common substitutions and plain deletions/duplications (RefSeq / Ensembl /
/// LRG / Assembly `g.`/`c.` — ~72% of real ClinVar) take a specialized fast
/// path that is ~1.7x
/// faster end-to-end over the full ClinVar corpus; everything else falls back
/// to the full [`variant::parse_variant`] parser. The fast path is
/// observationally identical to the generic parser — both accept the same
/// inputs and produce the same [`HgvsVariant`] (guarded by the differential
/// test in `tests/fast_path_differential.rs`) — so it is transparent to
/// callers. To force the generic parser, call [`variant::parse_variant`]
/// directly.
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
    let trimmed = input.trim();
    match fast_path::try_fast_path(trimmed) {
        fast_path::FastPathResult::Success(variant) => Ok(variant),
        // Fall back on the original (untrimmed) input so the generic parser
        // sees exactly what it did before the fast path existed.
        fast_path::FastPathResult::Fallback => variant::parse_variant(input),
    }
}

/// Parse an HGVS string using the fast path for common patterns.
///
/// **Equivalent to [`parse_hgvs`].** The fast path is now the default, so this
/// is a thin synonym retained for backward compatibility; new code should call
/// [`parse_hgvs`]. See [`parse_hgvs`] for the performance characteristics and
/// the identity guarantee versus the generic parser.
///
/// # Example
///
/// ```
/// use ferro_hgvs::parse_hgvs_fast;
///
/// let variant = parse_hgvs_fast("NC_000001.11:g.12345A>G").unwrap();
/// // Complex patterns transparently fall back to the full parser.
/// let variant = parse_hgvs_fast("NM_000088.3:c.100+5A>G").unwrap();
/// ```
#[inline]
pub fn parse_hgvs_fast(input: &str) -> Result<HgvsVariant, FerroError> {
    parse_hgvs(input)
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
    // Resolve the bracket-cardinality action before `config` is consumed by the
    // preprocessor; the rule itself is applied post-parse on the AST below.
    let cardinality_action = config.action_for(ErrorType::NonConformantBracketCardinality);

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

    // Apply the bracket/allele cardinality conformance rule on the parsed AST
    // (#493). This is a structural rule, identical across all coordinate
    // systems, so it is enforced once here rather than per-axis in the parser.
    let mut warnings = preprocess_result.warnings;
    let variant = apply_bracket_cardinality_rule(
        variant,
        cardinality_action,
        &preprocess_result.preprocessed,
        &mut warnings,
    )?;

    // Return result with warnings
    Ok(ParseResultWithWarnings::new(
        variant,
        warnings,
        preprocess_result.original,
        preprocess_result.preprocessed,
    ))
}

/// Enforce the HGVS bracket/allele cardinality rule on a parsed variant (#493).
///
/// `[ ]` is allele syntax. The spec admits only two conformant shapes, identically
/// across every coordinate system (`c/g/n/m/o/r/p`): one bracket group with two or
/// more cis members (`c.[76A>C;88G>T]`), or two or more trans groups
/// (`c.[76A>C];[88G>T]`). A standalone single-member bracket — `c.[76A>C]`,
/// `g.[1000G>A]`, `p.[=]`, `c.[?]`, `p.[(?)]` — is non-conformant
/// (`DNA/alleles.md`). It parses to an [`HgvsVariant::Allele`] with exactly one
/// member, which is the uniform signal this rule keys off.
///
/// The canonical repair drops the redundant brackets, unwrapping the singleton to
/// its sole member ([`AlleleVariant::Display`](crate::hgvs::variant::AlleleVariant)
/// already renders a singleton bracket-free, so the corrected output is identical).
/// `Reject` (strict) returns a `W3026` parse error; `WarnCorrect` (lenient) unwraps
/// and pushes a `W3026` warning; `SilentCorrect` (silent) unwraps quietly; `Accept`
/// leaves the wrapper intact. Conformant multi-member brackets pass through
/// untouched.
fn apply_bracket_cardinality_rule(
    variant: HgvsVariant,
    action: ResolvedAction,
    source: &str,
    warnings: &mut Vec<CorrectionWarning>,
) -> Result<HgvsVariant, FerroError> {
    // Only a standalone single-member allele bracket is non-conformant. Multi-member
    // cis groups and >=2 trans groups carry more than one member and pass through.
    let is_singleton = matches!(&variant, HgvsVariant::Allele(av) if av.variants.len() == 1);
    if !is_singleton {
        return Ok(variant);
    }

    match action {
        ResolvedAction::Accept => Ok(variant),
        ResolvedAction::Reject => {
            let canonical = variant.to_string();
            Err(FerroError::Parse {
                pos: 0,
                msg: format!(
                    "[{}] standalone single-member allele bracket is not HGVS-conformant: \
                     allele brackets require two or more cis members (`c.[a;b]`) or two or \
                     more trans groups (`c.[a];[b]`). Drop the brackets (canonical \
                     `{}`) or pair with a second allele (e.g. `[a];[=]`).",
                    ErrorType::NonConformantBracketCardinality.code(),
                    canonical,
                ),
                diagnostic: None,
            })
        }
        ResolvedAction::WarnCorrect | ResolvedAction::SilentCorrect => {
            // Canonical repair: unwrap the singleton wrapper to its sole member.
            let member = match variant {
                HgvsVariant::Allele(av) => av
                    .variants
                    .into_iter()
                    .next()
                    .expect("singleton allele has exactly one member"),
                // Unreachable: `is_singleton` already established this is an Allele.
                other => return Ok(other),
            };
            if action == ResolvedAction::WarnCorrect {
                let canonical = member.to_string();
                warnings.push(CorrectionWarning::new(
                    ErrorType::NonConformantBracketCardinality,
                    format!(
                        "standalone single-member allele bracket unwrapped to `{canonical}`; \
                         allele brackets require >=2 cis members or >=2 trans groups"
                    ),
                    None,
                    source.to_string(),
                    canonical,
                ));
            }
            Ok(member)
        }
    }
}

/// Parse an HGVS string with lenient error handling.
///
/// This is a convenience function that uses lenient mode, which auto-corrects
/// common errors and returns warnings.
///
/// # Examples
///
/// Auto-corrects an en-dash to a hyphen:
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let result = parse_hgvs_lenient("NM_000088.3:c.100\u{2013}200del");
/// assert!(result.is_ok());
/// ```
///
/// Soft-validation warnings emitted in lenient mode include W1001 for
/// lowercase amino-acid codes:
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let parsed = parse_hgvs_lenient("NP_000079.2:p.val600glu").unwrap();
/// assert_eq!(parsed.preprocessed_input, "NP_000079.2:p.Val600Glu");
/// assert!(parsed
///     .warnings
///     .iter()
///     .any(|w| w.error_type.code() == "W1001"));
/// ```
///
/// W1002 for one-letter amino-acid abbreviations:
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let parsed = parse_hgvs_lenient("NP_000079.2:p.V600E").unwrap();
/// assert_eq!(parsed.preprocessed_input, "NP_000079.2:p.Val600Glu");
/// assert!(parsed
///     .warnings
///     .iter()
///     .any(|w| w.error_type.code() == "W1002"));
/// ```
///
/// W3001 for accessions missing a `.<version>` suffix (the warning fires
/// but the input is not auto-corrected — the version cannot be synthesised):
///
/// ```
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// let parsed = parse_hgvs_lenient("NM_000088:c.100A>G").unwrap();
/// assert_eq!(parsed.preprocessed_input, "NM_000088:c.100A>G");
/// assert!(parsed
///     .warnings
///     .iter()
///     .any(|w| w.error_type.code() == "W3001"));
/// ```
///
/// ```
/// use ferro_hgvs::error_handling::ErrorType;
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // SVA-008: single-position range is collapsed (W4003)
/// let result = parse_hgvs_lenient("NM_000088.3:c.123_123del").unwrap();
/// assert_eq!(result.preprocessed_input, "NM_000088.3:c.123del");
/// assert!(result
///     .warnings
///     .iter()
///     .any(|w| w.error_type == ErrorType::SinglePositionRange));
/// ```
///
/// ```
/// use ferro_hgvs::error_handling::ErrorType;
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // SVA-010: empty delins is rewritten to del (W3012)
/// let result = parse_hgvs_lenient("NC_000001.11:g.100_102delins").unwrap();
/// assert_eq!(result.preprocessed_input, "NC_000001.11:g.100_102del");
/// assert!(result
///     .warnings
///     .iter()
///     .any(|w| w.error_type == ErrorType::EmptyDelinsInsert));
/// ```
///
/// ```
/// use ferro_hgvs::error_handling::ErrorType;
/// use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
///
/// // SVA-007: deletion with size-count suffix warns but is not rewritten (W3011)
/// let result = parse_hgvs_lenient("NG_012232.1:g.123del6").unwrap();
/// assert!(result
///     .warnings
///     .iter()
///     .any(|w| w.error_type == ErrorType::DelSizeSuffix));
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
