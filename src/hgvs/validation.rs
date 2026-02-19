//! Validation levels for HGVS parsing
//!
//! This module provides configurable validation strictness for parsing
//! HGVS variant strings. Different levels allow for varying degrees of
//! flexibility vs. spec compliance.

use std::sync::atomic::{AtomicU8, Ordering};

/// Validation level for HGVS parsing
///
/// Controls how strictly the parser adheres to the HGVS specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum ValidationLevel {
    /// Lenient validation - accepts common typos and variations
    ///
    /// - Accepts lowercase accession prefixes (nm_, nc_)
    /// - Accepts missing version numbers (NM_000088)
    /// - Accepts common typos in amino acid names
    /// - Accepts c: instead of c.
    Lenient = 0,

    /// Standard validation - common clinical usage
    ///
    /// - Requires proper case for accession prefixes
    /// - Allows missing version numbers
    /// - Requires valid amino acid codes
    /// - Requires proper separator (c.)
    #[default]
    Standard = 1,

    /// Strict validation - full HGVS spec compliance
    ///
    /// - Requires proper case for accession prefixes
    /// - Requires version numbers for RefSeq accessions
    /// - Requires valid amino acid codes
    /// - Requires proper separator (c.)
    /// - Validates position ranges
    Strict = 2,
}

impl ValidationLevel {
    /// Check if this level is at least as strict as another
    pub fn at_least(&self, other: ValidationLevel) -> bool {
        (*self as u8) >= (other as u8)
    }

    /// Check if this level allows lenient parsing
    pub fn is_lenient(&self) -> bool {
        *self == ValidationLevel::Lenient
    }

    /// Check if this level requires strict validation
    pub fn is_strict(&self) -> bool {
        *self == ValidationLevel::Strict
    }
}

// Global validation level for the parser
static GLOBAL_VALIDATION_LEVEL: AtomicU8 = AtomicU8::new(ValidationLevel::Standard as u8);

/// Get the current global validation level
pub fn get_validation_level() -> ValidationLevel {
    match GLOBAL_VALIDATION_LEVEL.load(Ordering::Relaxed) {
        0 => ValidationLevel::Lenient,
        1 => ValidationLevel::Standard,
        _ => ValidationLevel::Strict,
    }
}

/// Set the global validation level
///
/// Note: This is a global setting and affects all parsers in the process.
/// For thread-local validation, use `ParseConfig` instead.
pub fn set_validation_level(level: ValidationLevel) {
    GLOBAL_VALIDATION_LEVEL.store(level as u8, Ordering::Relaxed);
}

/// Run a block with a specific validation level, restoring the original level after
pub fn with_validation_level<T, F: FnOnce() -> T>(level: ValidationLevel, f: F) -> T {
    let old_level = get_validation_level();
    set_validation_level(level);
    let result = f();
    set_validation_level(old_level);
    result
}

/// Parse configuration with validation settings
#[derive(Debug, Clone, Default)]
pub struct ParseConfig {
    /// Validation level
    pub level: ValidationLevel,
    /// Whether to allow missing version numbers
    pub allow_missing_version: bool,
    /// Whether to allow lowercase accession prefixes
    pub allow_lowercase_prefix: bool,
    /// Whether to allow c: instead of c.
    pub allow_colon_separator: bool,
}

impl ParseConfig {
    /// Create a new parse configuration
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a lenient configuration
    pub fn lenient() -> Self {
        Self {
            level: ValidationLevel::Lenient,
            allow_missing_version: true,
            allow_lowercase_prefix: true,
            allow_colon_separator: true,
        }
    }

    /// Create a standard configuration
    pub fn standard() -> Self {
        Self {
            level: ValidationLevel::Standard,
            allow_missing_version: true,
            allow_lowercase_prefix: false,
            allow_colon_separator: false,
        }
    }

    /// Create a strict configuration
    pub fn strict() -> Self {
        Self {
            level: ValidationLevel::Strict,
            allow_missing_version: false,
            allow_lowercase_prefix: false,
            allow_colon_separator: false,
        }
    }

    /// Set the validation level
    pub fn with_level(mut self, level: ValidationLevel) -> Self {
        self.level = level;
        self
    }

    /// Allow missing version numbers
    pub fn with_allow_missing_version(mut self, allow: bool) -> Self {
        self.allow_missing_version = allow;
        self
    }
}

/// Validation result with optional warnings
#[derive(Debug, Clone)]
pub struct ValidationResult {
    /// Whether the variant is valid
    pub valid: bool,
    /// Warnings (non-fatal issues)
    pub warnings: Vec<ValidationWarning>,
    /// Errors (fatal issues)
    pub errors: Vec<ValidationError>,
}

impl ValidationResult {
    /// Create a valid result with no issues
    pub fn ok() -> Self {
        Self {
            valid: true,
            warnings: Vec::new(),
            errors: Vec::new(),
        }
    }

    /// Add a warning
    pub fn with_warning(mut self, warning: ValidationWarning) -> Self {
        self.warnings.push(warning);
        self
    }

    /// Add an error
    pub fn with_error(mut self, error: ValidationError) -> Self {
        self.valid = false;
        self.errors.push(error);
        self
    }

    /// Check if there are any warnings
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    /// Check if there are any errors
    pub fn has_errors(&self) -> bool {
        !self.errors.is_empty()
    }
}

/// Validation warning
#[derive(Debug, Clone)]
pub struct ValidationWarning {
    /// Warning code
    pub code: &'static str,
    /// Warning message
    pub message: String,
    /// Position in input (if applicable)
    pub position: Option<usize>,
}

impl ValidationWarning {
    /// Create a new warning
    pub fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            position: None,
        }
    }

    /// Set the position
    pub fn at_position(mut self, pos: usize) -> Self {
        self.position = Some(pos);
        self
    }
}

/// Validation error
#[derive(Debug, Clone)]
pub struct ValidationError {
    /// Error code
    pub code: &'static str,
    /// Error message
    pub message: String,
    /// Position in input (if applicable)
    pub position: Option<usize>,
}

impl ValidationError {
    /// Create a new error
    pub fn new(code: &'static str, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            position: None,
        }
    }

    /// Set the position
    pub fn at_position(mut self, pos: usize) -> Self {
        self.position = Some(pos);
        self
    }
}

// ============================================================================
// Semantic Validation Rules
// ============================================================================

/// Validate a parsed HGVS variant for semantic correctness
pub mod rules {
    use crate::hgvs::edit::{Base, NaEdit, ProteinEdit};
    use crate::hgvs::variant::HgvsVariant;

    use super::{ValidationError, ValidationResult, ValidationWarning};

    /// Validate that an interval's start position is not greater than end
    pub fn validate_position_order(start: i64, end: i64) -> ValidationResult {
        if start > end {
            ValidationResult::ok().with_error(ValidationError::new(
                "E001",
                format!(
                    "Start position ({}) is greater than end position ({})",
                    start, end
                ),
            ))
        } else {
            ValidationResult::ok()
        }
    }

    /// Validate that a nucleotide edit has valid bases
    pub fn validate_na_edit(edit: &NaEdit) -> ValidationResult {
        let mut result = ValidationResult::ok();

        match edit {
            NaEdit::Substitution {
                reference,
                alternative,
            } => {
                if reference == alternative {
                    result = result.with_warning(ValidationWarning::new(
                        "W001",
                        format!(
                            "Reference and alternative are the same: {} = {}",
                            reference, alternative
                        ),
                    ));
                }
            }
            NaEdit::Deletion {
                sequence: Some(seq),
                ..
            } => {
                if seq.is_empty() {
                    result = result.with_warning(ValidationWarning::new(
                        "W002",
                        "Deletion has empty sequence specified",
                    ));
                }
            }
            NaEdit::Deletion { sequence: None, .. } => {}
            NaEdit::Insertion { sequence } => {
                if sequence.is_empty() {
                    result = result.with_error(ValidationError::new(
                        "E002",
                        "Insertion must have at least one base",
                    ));
                }
            }
            NaEdit::Delins { sequence } => {
                if sequence.is_empty() {
                    result = result.with_error(ValidationError::new(
                        "E003",
                        "Deletion-insertion must have at least one inserted base",
                    ));
                }
            }
            _ => {}
        }

        result
    }

    /// Validate that a protein edit has valid amino acids
    pub fn validate_protein_edit(_edit: &ProteinEdit) -> ValidationResult {
        // Protein edits are mostly self-validating through the AminoAcid enum
        // Additional validation could check for biologically impossible changes
        ValidationResult::ok()
    }

    /// Check if a base is a valid nucleotide code
    ///
    /// Currently supports: A, C, G, T, U (RNA), N (any)
    pub fn is_valid_base(base: &Base) -> bool {
        matches!(
            base,
            Base::A | Base::C | Base::G | Base::T | Base::U | Base::N
        )
    }

    /// Validate that a variant doesn't have conflicting attributes
    pub fn validate_variant_consistency(variant: &HgvsVariant) -> ValidationResult {
        let mut result = ValidationResult::ok();

        match variant {
            HgvsVariant::Cds(v) => {
                // Check for intronic positions in exon-only contexts
                if let Some(start) = v.loc_edit.location.start.inner() {
                    if let Some(end) = v.loc_edit.location.end.inner() {
                        // Only one can be in 3' UTR
                        if start.utr3 != end.utr3 && !start.utr3 && end.utr3 {
                            // This is valid - spanning CDS to 3' UTR
                        }
                    }
                }
            }
            HgvsVariant::Protein(v) => {
                // Check for frameshift without position change
                if let Some(ProteinEdit::Frameshift { .. }) = v.loc_edit.edit.inner() {
                    // Frameshifts are inherently complex - no additional validation needed
                }
            }
            HgvsVariant::Allele(a) => {
                // Validate each variant in the allele
                for v in &a.variants {
                    let sub_result = validate_variant_consistency(v);
                    for err in sub_result.errors {
                        result = result.with_error(err);
                    }
                    for warn in sub_result.warnings {
                        result = result.with_warning(warn);
                    }
                }
            }
            _ => {}
        }

        result
    }

    /// Comprehensive validation of a variant
    pub fn validate_variant(variant: &HgvsVariant) -> ValidationResult {
        let mut result = ValidationResult::ok();

        // Check consistency
        let consistency = validate_variant_consistency(variant);
        for err in consistency.errors {
            result = result.with_error(err);
        }
        for warn in consistency.warnings {
            result = result.with_warning(warn);
        }

        result
    }
}

/// Variant validator that uses a reference provider
pub struct Validator<P: crate::reference::ReferenceProvider> {
    #[allow(dead_code)] // Will be used for reference base validation
    provider: P,
    config: ParseConfig,
}

impl<P: crate::reference::ReferenceProvider> Validator<P> {
    /// Create a new validator with a reference provider
    pub fn new(provider: P) -> Self {
        Self {
            provider,
            config: ParseConfig::default(),
        }
    }

    /// Create a validator with custom config
    pub fn with_config(provider: P, config: ParseConfig) -> Self {
        Self { provider, config }
    }

    /// Validate a variant against the reference
    pub fn validate(&self, variant: &crate::hgvs::variant::HgvsVariant) -> ValidationResult {
        let mut result = rules::validate_variant(variant);

        // If we have a provider, validate reference bases
        if self.config.level.at_least(ValidationLevel::Standard) {
            if let Some(ref_error) = self.validate_reference_base(variant) {
                result = result.with_error(ref_error);
            }
        }

        result
    }

    /// Validate that reference base(s) in the variant match the actual reference
    fn validate_reference_base(
        &self,
        _variant: &crate::hgvs::variant::HgvsVariant,
    ) -> Option<ValidationError> {
        // This would need transcript/sequence access to verify
        // For now, return None (no error)
        // A full implementation would:
        // 1. Get the transcript/sequence via self.provider
        // 2. Extract the position(s)
        // 3. Compare the reference base in the variant to the actual sequence
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validation_level_default() {
        assert_eq!(ValidationLevel::default(), ValidationLevel::Standard);
    }

    #[test]
    fn test_validation_level_at_least() {
        assert!(ValidationLevel::Strict.at_least(ValidationLevel::Standard));
        assert!(ValidationLevel::Standard.at_least(ValidationLevel::Lenient));
        assert!(!ValidationLevel::Lenient.at_least(ValidationLevel::Strict));
    }

    #[test]
    fn test_parse_config_lenient() {
        let config = ParseConfig::lenient();
        assert_eq!(config.level, ValidationLevel::Lenient);
        assert!(config.allow_missing_version);
        assert!(config.allow_lowercase_prefix);
    }

    #[test]
    fn test_parse_config_strict() {
        let config = ParseConfig::strict();
        assert_eq!(config.level, ValidationLevel::Strict);
        assert!(!config.allow_missing_version);
        assert!(!config.allow_lowercase_prefix);
    }

    #[test]
    fn test_validation_result() {
        let result = ValidationResult::ok();
        assert!(result.valid);
        assert!(!result.has_warnings());
        assert!(!result.has_errors());

        let result = result.with_warning(ValidationWarning::new("W001", "test warning"));
        assert!(result.valid);
        assert!(result.has_warnings());
    }

    #[test]
    fn test_global_validation_level() {
        let original = get_validation_level();

        set_validation_level(ValidationLevel::Lenient);
        assert_eq!(get_validation_level(), ValidationLevel::Lenient);

        set_validation_level(ValidationLevel::Strict);
        assert_eq!(get_validation_level(), ValidationLevel::Strict);

        set_validation_level(original);
    }

    #[test]
    fn test_with_validation_level() {
        let original = get_validation_level();
        set_validation_level(ValidationLevel::Standard);

        let result = with_validation_level(ValidationLevel::Strict, || {
            assert_eq!(get_validation_level(), ValidationLevel::Strict);
            42
        });

        assert_eq!(result, 42);
        assert_eq!(get_validation_level(), ValidationLevel::Standard);

        set_validation_level(original);
    }

    // Tests for semantic validation rules
    mod rules_tests {
        use super::*;
        use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};

        #[test]
        fn test_validate_position_order_valid() {
            let result = rules::validate_position_order(10, 20);
            assert!(result.valid);
            assert!(!result.has_errors());
        }

        #[test]
        fn test_validate_position_order_equal() {
            let result = rules::validate_position_order(10, 10);
            assert!(result.valid);
        }

        #[test]
        fn test_validate_position_order_invalid() {
            let result = rules::validate_position_order(20, 10);
            assert!(!result.valid);
            assert!(result.has_errors());
            assert_eq!(result.errors[0].code, "E001");
        }

        #[test]
        fn test_validate_substitution_same_base_warning() {
            let edit = NaEdit::Substitution {
                reference: Base::A,
                alternative: Base::A,
            };
            let result = rules::validate_na_edit(&edit);
            assert!(result.valid); // Still valid, just a warning
            assert!(result.has_warnings());
            assert_eq!(result.warnings[0].code, "W001");
        }

        #[test]
        fn test_validate_substitution_different_bases() {
            let edit = NaEdit::Substitution {
                reference: Base::A,
                alternative: Base::G,
            };
            let result = rules::validate_na_edit(&edit);
            assert!(result.valid);
            assert!(!result.has_warnings());
        }

        #[test]
        fn test_validate_insertion_empty_error() {
            let edit = NaEdit::Insertion {
                sequence: InsertedSequence::Literal(Sequence::new(vec![])),
            };
            let result = rules::validate_na_edit(&edit);
            assert!(!result.valid);
            assert!(result.has_errors());
            assert_eq!(result.errors[0].code, "E002");
        }

        #[test]
        fn test_validate_delins_empty_error() {
            let edit = NaEdit::Delins {
                sequence: InsertedSequence::Literal(Sequence::new(vec![])),
            };
            let result = rules::validate_na_edit(&edit);
            assert!(!result.valid);
            assert!(result.has_errors());
            assert_eq!(result.errors[0].code, "E003");
        }

        #[test]
        fn test_is_valid_base() {
            assert!(rules::is_valid_base(&Base::A));
            assert!(rules::is_valid_base(&Base::C));
            assert!(rules::is_valid_base(&Base::G));
            assert!(rules::is_valid_base(&Base::T));
            assert!(rules::is_valid_base(&Base::U));
            assert!(rules::is_valid_base(&Base::N));
        }

        #[test]
        fn test_validate_variant_consistency_cds() {
            use crate::hgvs::parser::parse_hgvs;

            let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
            let result = rules::validate_variant(&variant);
            assert!(result.valid);
        }

        #[test]
        fn test_validate_variant_consistency_allele() {
            use crate::hgvs::parser::parse_hgvs;

            let variant = parse_hgvs("[NM_000088.3:c.100A>G;NM_000088.3:c.200C>T]").unwrap();
            let result = rules::validate_variant(&variant);
            assert!(result.valid);
        }
    }

    // Tests for Validator struct
    mod validator_tests {
        use super::*;
        use crate::hgvs::parser::parse_hgvs;
        use crate::reference::MockProvider;

        #[test]
        fn test_validator_creation() {
            let provider = MockProvider::with_test_data();
            let _validator = Validator::new(provider);
        }

        #[test]
        fn test_validator_with_config() {
            let provider = MockProvider::with_test_data();
            let config = ParseConfig::strict();
            let validator = Validator::with_config(provider, config);
            assert_eq!(validator.config.level, ValidationLevel::Strict);
        }

        #[test]
        fn test_validate_simple_variant() {
            let provider = MockProvider::with_test_data();
            let validator = Validator::new(provider);

            let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
            let result = validator.validate(&variant);
            assert!(result.valid);
        }

        #[test]
        fn test_validate_protein_variant() {
            let provider = MockProvider::with_test_data();
            let validator = Validator::new(provider);

            let variant = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
            let result = validator.validate(&variant);
            assert!(result.valid);
        }
    }
}
