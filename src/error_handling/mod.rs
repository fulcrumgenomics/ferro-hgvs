//! Configurable error handling for HGVS parsing.
//!
//! This module provides three error handling modes (Strict, Lenient, Silent)
//! plus per-error-type configuration overrides. This allows users to control
//! how the parser handles common input errors like wrong dash characters,
//! lowercase amino acids, and extra whitespace.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::error_handling::{ErrorConfig, ErrorMode, ErrorType, ErrorOverride};
//!
//! // Strict mode (default): reject all non-standard input
//! let config = ErrorConfig::strict();
//!
//! // Lenient mode: auto-correct with warnings
//! let config = ErrorConfig::lenient();
//!
//! // Silent mode: auto-correct without warnings
//! let config = ErrorConfig::silent();
//!
//! // Custom: lenient mode but reject lowercase amino acids
//! let config = ErrorConfig::lenient()
//!     .with_override(ErrorType::LowercaseAminoAcid, ErrorOverride::Reject);
//! ```
//!
//! # Error Types
//!
//! The following error types can be individually configured:
//!
//! | Error Type | Example | Description |
//! |------------|---------|-------------|
//! | `LowercaseAminoAcid` | `val` → `Val` | Lowercase 3-letter amino acid codes |
//! | `MissingVersion` | `NM_000088` → `NM_000088.3` | Missing accession version |
//! | `WrongDashCharacter` | `–` → `-` | En-dash/em-dash instead of hyphen |
//! | `ExtraWhitespace` | `c.100 A>G` → `c.100A>G` | Extra spaces in description |
//! | `ProteinSubstitutionArrow` | `Val600>Glu` → `Val600Glu` | Arrow in protein substitution |
//! | `PositionZero` | `c.0A>G` | Invalid position zero (always rejected) |
//!
//! # Override Behaviors
//!
//! Each error type can have one of these override behaviors:
//!
//! - `Default`: Use the base mode's behavior
//! - `Reject`: Always return an error
//! - `WarnCorrect`: Auto-correct and emit a warning
//! - `SilentCorrect`: Auto-correct without warning
//! - `Accept`: Accept the input as-is without correction

pub mod codes;
pub mod corrections;
mod preprocessor;
pub mod registry;
mod types;

pub use codes::{CodeCategory, CodeInfo, ModeAction, ModeBehavior};
pub use corrections::{
    detect_accession_typo, detect_amino_acid_typo, detect_edit_type_typo, detect_missing_version,
    detect_swapped_positions, detect_typos, find_closest_match, levenshtein_distance,
    DetectedCorrection, FuzzyMatch, TypoSuggestion, TypoTokenType,
};
pub use preprocessor::{CorrectionWarning, InputPreprocessor, PreprocessResult};
pub use registry::{get_code_info, list_all_codes, list_error_codes, list_warning_codes};
pub use types::{ErrorMode, ErrorOverride, ErrorType, ResolvedAction};

use std::collections::HashMap;

/// Error handling configuration.
///
/// Controls how the parser handles common input errors. The configuration
/// consists of a base mode plus optional per-error-type overrides.
#[derive(Debug, Clone)]
pub struct ErrorConfig {
    /// Base error handling mode.
    pub mode: ErrorMode,
    /// Per-error-type overrides.
    pub overrides: HashMap<ErrorType, ErrorOverride>,
}

impl ErrorConfig {
    /// Create a new configuration with the given mode.
    pub fn new(mode: ErrorMode) -> Self {
        Self {
            mode,
            overrides: HashMap::new(),
        }
    }

    /// Create a strict configuration.
    ///
    /// In strict mode, all non-standard input is rejected.
    pub fn strict() -> Self {
        Self::new(ErrorMode::Strict)
    }

    /// Create a lenient configuration.
    ///
    /// In lenient mode, common errors are auto-corrected with warnings.
    pub fn lenient() -> Self {
        Self::new(ErrorMode::Lenient)
    }

    /// Create a silent configuration.
    ///
    /// In silent mode, common errors are auto-corrected without warnings.
    pub fn silent() -> Self {
        Self::new(ErrorMode::Silent)
    }

    /// Add an override for a specific error type.
    ///
    /// # Example
    ///
    /// ```
    /// use ferro_hgvs::error_handling::{ErrorConfig, ErrorType, ErrorOverride};
    ///
    /// let config = ErrorConfig::lenient()
    ///     .with_override(ErrorType::LowercaseAminoAcid, ErrorOverride::Reject)
    ///     .with_override(ErrorType::ExtraWhitespace, ErrorOverride::SilentCorrect);
    /// ```
    pub fn with_override(mut self, error_type: ErrorType, override_: ErrorOverride) -> Self {
        self.overrides.insert(error_type, override_);
        self
    }

    /// Set an override for a specific error type.
    ///
    /// This is the mutable version of `with_override`.
    pub fn set_override(&mut self, error_type: ErrorType, override_: ErrorOverride) {
        self.overrides.insert(error_type, override_);
    }

    /// Remove an override for a specific error type.
    pub fn remove_override(&mut self, error_type: ErrorType) {
        self.overrides.remove(&error_type);
    }

    /// Get the resolved action for an error type.
    ///
    /// This applies the override if one exists, otherwise uses the base mode.
    pub fn action_for(&self, error_type: ErrorType) -> ResolvedAction {
        let override_ = self.overrides.get(&error_type).copied().unwrap_or_default();
        override_.resolve(self.mode)
    }

    /// Returns true if the given error type should be rejected.
    pub fn should_reject(&self, error_type: ErrorType) -> bool {
        self.action_for(error_type).should_reject()
    }

    /// Returns true if the given error type should be corrected.
    pub fn should_correct(&self, error_type: ErrorType) -> bool {
        self.action_for(error_type).should_correct()
    }

    /// Returns true if the given error type should emit a warning.
    pub fn should_warn(&self, error_type: ErrorType) -> bool {
        self.action_for(error_type).should_warn()
    }

    /// Create a preprocessor with this configuration.
    pub fn preprocessor(&self) -> InputPreprocessor {
        InputPreprocessor::new(self.clone())
    }
}

impl Default for ErrorConfig {
    fn default() -> Self {
        Self::strict()
    }
}

/// Parse result with warnings.
///
/// This struct wraps a parsed variant along with any warnings
/// generated during preprocessing or parsing.
#[derive(Debug, Clone)]
pub struct ParseResultWithWarnings<T> {
    /// The parsed result.
    pub result: T,
    /// Warnings generated during parsing.
    pub warnings: Vec<CorrectionWarning>,
    /// The original input.
    pub original_input: String,
    /// The preprocessed input (may be same as original).
    pub preprocessed_input: String,
}

impl<T> ParseResultWithWarnings<T> {
    /// Create a new parse result with warnings.
    pub fn new(
        result: T,
        warnings: Vec<CorrectionWarning>,
        original_input: String,
        preprocessed_input: String,
    ) -> Self {
        Self {
            result,
            warnings,
            original_input,
            preprocessed_input,
        }
    }

    /// Create a parse result without warnings.
    pub fn without_warnings(result: T, input: String) -> Self {
        Self {
            result,
            warnings: Vec::new(),
            original_input: input.clone(),
            preprocessed_input: input,
        }
    }

    /// Returns true if there were any corrections made.
    pub fn had_corrections(&self) -> bool {
        self.original_input != self.preprocessed_input
    }

    /// Returns true if there are any warnings.
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }

    /// Map the result to a new type.
    pub fn map<U, F: FnOnce(T) -> U>(self, f: F) -> ParseResultWithWarnings<U> {
        ParseResultWithWarnings {
            result: f(self.result),
            warnings: self.warnings,
            original_input: self.original_input,
            preprocessed_input: self.preprocessed_input,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ErrorConfig tests
    #[test]
    fn test_error_config_default() {
        let config = ErrorConfig::default();
        assert_eq!(config.mode, ErrorMode::Strict);
        assert!(config.overrides.is_empty());
    }

    #[test]
    fn test_error_config_strict() {
        let config = ErrorConfig::strict();
        assert_eq!(config.mode, ErrorMode::Strict);
        assert!(config.should_reject(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_error_config_lenient() {
        let config = ErrorConfig::lenient();
        assert_eq!(config.mode, ErrorMode::Lenient);
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
        assert!(config.should_warn(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_error_config_silent() {
        let config = ErrorConfig::silent();
        assert_eq!(config.mode, ErrorMode::Silent);
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
        assert!(!config.should_warn(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_error_config_with_override() {
        let config = ErrorConfig::lenient()
            .with_override(ErrorType::LowercaseAminoAcid, ErrorOverride::Reject);

        // Overridden error type should reject
        assert!(config.should_reject(ErrorType::LowercaseAminoAcid));

        // Non-overridden error types should use base mode
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_error_config_set_override() {
        let mut config = ErrorConfig::strict();
        config.set_override(ErrorType::WrongDashCharacter, ErrorOverride::SilentCorrect);
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
        assert!(!config.should_warn(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_error_config_remove_override() {
        let mut config = ErrorConfig::lenient()
            .with_override(ErrorType::WrongDashCharacter, ErrorOverride::Reject);

        // With override
        assert!(config.should_reject(ErrorType::WrongDashCharacter));

        // After removing override
        config.remove_override(ErrorType::WrongDashCharacter);
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_error_config_action_for() {
        let config =
            ErrorConfig::lenient().with_override(ErrorType::PositionZero, ErrorOverride::Reject);

        assert_eq!(
            config.action_for(ErrorType::WrongDashCharacter),
            ResolvedAction::WarnCorrect
        );
        assert_eq!(
            config.action_for(ErrorType::PositionZero),
            ResolvedAction::Reject
        );
    }

    #[test]
    fn test_error_config_preprocessor() {
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();

        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(result.success);
        assert_eq!(result.preprocessed, "c.100-200del");
    }

    // ParseResultWithWarnings tests
    #[test]
    fn test_parse_result_with_warnings_new() {
        let result = ParseResultWithWarnings::new(
            42,
            vec![CorrectionWarning::new(
                ErrorType::WrongDashCharacter,
                "test",
                None,
                "",
                "",
            )],
            "original".to_string(),
            "preprocessed".to_string(),
        );

        assert_eq!(result.result, 42);
        assert!(result.has_warnings());
        assert!(result.had_corrections());
    }

    #[test]
    fn test_parse_result_with_warnings_without_warnings() {
        let result = ParseResultWithWarnings::without_warnings(42, "input".to_string());

        assert_eq!(result.result, 42);
        assert!(!result.has_warnings());
        assert!(!result.had_corrections());
    }

    #[test]
    fn test_parse_result_with_warnings_map() {
        let result = ParseResultWithWarnings::without_warnings(42, "input".to_string());
        let mapped = result.map(|x| x.to_string());

        assert_eq!(mapped.result, "42");
    }
}
