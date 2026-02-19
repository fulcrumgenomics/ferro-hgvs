//! Input preprocessor for normalizing HGVS input strings.
//!
//! The preprocessor applies corrections to common input errors before
//! parsing, based on the configured error handling mode.

use super::corrections::{
    correct_accession_prefix_case, correct_dash_characters, correct_missing_coordinate_prefix,
    correct_old_allele_format, correct_protein_arrow, correct_quote_characters, correct_whitespace,
    detect_position_zero, strip_trailing_annotation, DetectedCorrection,
};
use super::types::{ErrorType, ResolvedAction};
use super::ErrorConfig;
use crate::error::{Diagnostic, ErrorCode, FerroError, SourceSpan};

/// Warning about a correction made during preprocessing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CorrectionWarning {
    /// The type of error that was corrected.
    pub error_type: ErrorType,
    /// Human-readable message about the correction.
    pub message: String,
    /// Position in the original input (start, end).
    pub span: Option<(usize, usize)>,
    /// The original value that was corrected.
    pub original: String,
    /// The corrected value.
    pub corrected: String,
}

impl CorrectionWarning {
    /// Create a new correction warning.
    pub fn new(
        error_type: ErrorType,
        message: impl Into<String>,
        span: Option<(usize, usize)>,
        original: impl Into<String>,
        corrected: impl Into<String>,
    ) -> Self {
        Self {
            error_type,
            message: message.into(),
            span,
            original: original.into(),
            corrected: corrected.into(),
        }
    }

    /// Create from a DetectedCorrection.
    pub fn from_correction(correction: &DetectedCorrection) -> Self {
        Self {
            error_type: correction.error_type,
            message: correction.warning_message(),
            span: Some((correction.start, correction.end)),
            original: correction.original.clone(),
            corrected: correction.corrected.clone(),
        }
    }
}

/// Result of preprocessing an input string.
#[derive(Debug, Clone)]
pub struct PreprocessResult {
    /// The original input.
    pub original: String,
    /// The preprocessed input (may be same as original if no corrections).
    pub preprocessed: String,
    /// Warnings generated during preprocessing.
    pub warnings: Vec<CorrectionWarning>,
    /// Whether preprocessing was successful (no rejected errors).
    pub success: bool,
    /// Error if preprocessing failed due to a rejected error.
    pub error: Option<FerroError>,
}

impl PreprocessResult {
    /// Create a successful result with no changes.
    pub fn unchanged(input: String) -> Self {
        Self {
            original: input.clone(),
            preprocessed: input,
            warnings: Vec::new(),
            success: true,
            error: None,
        }
    }

    /// Create a successful result with corrections.
    pub fn corrected(
        original: String,
        preprocessed: String,
        warnings: Vec<CorrectionWarning>,
    ) -> Self {
        Self {
            original,
            preprocessed,
            warnings,
            success: true,
            error: None,
        }
    }

    /// Create a failed result.
    pub fn failed(original: String, error: FerroError) -> Self {
        Self {
            original: original.clone(),
            preprocessed: original,
            warnings: Vec::new(),
            success: false,
            error: Some(error),
        }
    }

    /// Returns true if there were any corrections made.
    pub fn has_corrections(&self) -> bool {
        self.original != self.preprocessed
    }

    /// Returns true if there were any warnings.
    pub fn has_warnings(&self) -> bool {
        !self.warnings.is_empty()
    }
}

/// Input preprocessor that normalizes HGVS input strings.
#[derive(Debug, Clone)]
pub struct InputPreprocessor {
    /// The error handling configuration.
    config: ErrorConfig,
}

impl InputPreprocessor {
    /// Create a new preprocessor with the given configuration.
    pub fn new(config: ErrorConfig) -> Self {
        Self { config }
    }

    /// Create a preprocessor with strict mode.
    pub fn strict() -> Self {
        Self::new(ErrorConfig::strict())
    }

    /// Create a preprocessor with lenient mode.
    pub fn lenient() -> Self {
        Self::new(ErrorConfig::lenient())
    }

    /// Create a preprocessor with silent mode.
    pub fn silent() -> Self {
        Self::new(ErrorConfig::silent())
    }

    /// Get the resolved action for an error type.
    fn action_for(&self, error_type: ErrorType) -> ResolvedAction {
        self.config.action_for(error_type)
    }

    /// Preprocess the input string.
    ///
    /// Applies corrections based on the configured error handling mode.
    pub fn preprocess(&self, input: &str) -> PreprocessResult {
        // Start with the original input
        let mut current = input.to_string();
        let mut all_warnings = Vec::new();

        // Phase 1: Check for position zero (never correctable, always rejected)
        // Position zero is a fundamental HGVS error that cannot be auto-corrected,
        // so we always reject it regardless of the error handling mode.
        if let Some(pos) = detect_position_zero(&current) {
            return PreprocessResult::failed(
                input.to_string(),
                FerroError::parse_with_diagnostic(
                    pos,
                    "Position 0 is not valid in HGVS notation",
                    Diagnostic::new()
                        .with_code(ErrorCode::InvalidPosition)
                        .with_span(SourceSpan::new(pos, pos + 1))
                        .with_source(input)
                        .with_hint("HGVS positions start at 1, not 0"),
                ),
            );
        }

        // Phase 2: Normalize dash characters
        let (corrected, corrections) = correct_dash_characters(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::WrongDashCharacter);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!("Invalid dash character '{}', expected '-'", first.original),
                            Diagnostic::new()
                                .with_code(ErrorCode::UnexpectedChar)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone()),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {
                    // Keep original
                }
            }
        }

        // Phase 3: Normalize quote characters
        let (corrected, corrections) = correct_quote_characters(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::WrongQuoteCharacter);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Invalid quote character '{}', expected regular quotes",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::UnexpectedChar)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone()),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Phase 4: Normalize whitespace
        let (corrected, corrections) = correct_whitespace(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::ExtraWhitespace);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Extra whitespace in HGVS description",
                            Diagnostic::new()
                                .with_code(ErrorCode::UnexpectedChar)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone()),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Phase 5: Correct accession prefix case
        let (corrected, corrections) = correct_accession_prefix_case(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::LowercaseAccessionPrefix);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Lowercase accession prefix '{}', expected uppercase",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidAccession)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone()),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Phase 6: Correct protein arrow syntax
        let (corrected, corrections) = correct_protein_arrow(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::ProteinSubstitutionArrow);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Arrow '>' in protein substitution is not standard HGVS",
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint("Use p.Val600Glu instead of p.Val600>Glu"),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Phase 7: Infer missing coordinate prefix for genomic accessions
        let (corrected, corrections) = correct_missing_coordinate_prefix(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::MissingCoordinatePrefix);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Missing coordinate type prefix (e.g., 'g.' for genomic)",
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidAccession)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint("For genomic accessions (NC_, NG_), add 'g.' before the position"),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Phase 8: Strip trailing protein annotations (e.g., "(p.Lys236=)")
        let (corrected, corrections) = strip_trailing_annotation(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::TrailingAnnotation);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Trailing annotation '{}' is not valid HGVS syntax",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::UnexpectedChar)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint("Protein consequence annotations should be separate from the HGVS expression"),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Phase 9: Correct old allele format (e.g., ":[c.100A>G;c.200C>T]" → ":c.[100A>G;200C>T]")
        let (corrected, corrections) = correct_old_allele_format(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::OldAlleleFormat);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Old/deprecated allele format with coordinate type inside brackets",
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "Use c.[edit1;edit2] format instead of [c.edit1;c.edit2]",
                                ),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &corrections {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                    current = corrected;
                }
                ResolvedAction::SilentCorrect => {
                    current = corrected;
                }
                ResolvedAction::Accept => {}
            }
        }

        // Return result
        if current == input && all_warnings.is_empty() {
            PreprocessResult::unchanged(input.to_string())
        } else {
            PreprocessResult::corrected(input.to_string(), current, all_warnings)
        }
    }
}

impl Default for InputPreprocessor {
    fn default() -> Self {
        Self::strict()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_handling::ErrorOverride;

    // PreprocessResult tests
    #[test]
    fn test_preprocess_result_unchanged() {
        let result = PreprocessResult::unchanged("c.100A>G".to_string());
        assert!(result.success);
        assert!(!result.has_corrections());
        assert!(!result.has_warnings());
        assert_eq!(result.original, "c.100A>G");
        assert_eq!(result.preprocessed, "c.100A>G");
    }

    #[test]
    fn test_preprocess_result_corrected() {
        let result = PreprocessResult::corrected(
            "c.100\u{2013}200del".to_string(),
            "c.100-200del".to_string(),
            vec![CorrectionWarning::new(
                ErrorType::WrongDashCharacter,
                "test warning",
                Some((5, 8)),
                "\u{2013}",
                "-",
            )],
        );
        assert!(result.success);
        assert!(result.has_corrections());
        assert!(result.has_warnings());
    }

    // InputPreprocessor strict mode tests
    #[test]
    fn test_preprocessor_strict_valid_input() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("c.100A>G");
        assert!(result.success);
        assert!(!result.has_corrections());
    }

    #[test]
    fn test_preprocessor_strict_rejects_en_dash() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_preprocessor_strict_rejects_whitespace() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("  c.100A>G  ");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_strict_rejects_position_zero() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("c.0A>G");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    // InputPreprocessor lenient mode tests
    #[test]
    fn test_preprocessor_lenient_corrects_en_dash() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(result.success);
        assert_eq!(result.preprocessed, "c.100-200del");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_lenient_corrects_whitespace() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("  c.100A>G  ");
        assert!(result.success);
        assert_eq!(result.preprocessed, "c.100A>G");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_lenient_corrects_protein_arrow() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("p.Val600>Glu");
        assert!(result.success);
        assert_eq!(result.preprocessed, "p.Val600Glu");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_lenient_rejects_position_zero() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("c.0A>G");
        // Position zero is always rejected
        assert!(!result.success);
    }

    // InputPreprocessor silent mode tests
    #[test]
    fn test_preprocessor_silent_corrects_without_warnings() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(result.success);
        assert_eq!(result.preprocessed, "c.100-200del");
        assert!(!result.has_warnings());
    }

    #[test]
    fn test_preprocessor_silent_corrects_multiple() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("  nm_000088.3:c.100\u{2013}200del  ");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.100-200del");
        assert!(!result.has_warnings());
    }

    // Override tests
    #[test]
    fn test_preprocessor_override_reject_in_lenient() {
        let config = ErrorConfig::lenient()
            .with_override(ErrorType::WrongDashCharacter, ErrorOverride::Reject);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_override_silent_in_lenient() {
        let config = ErrorConfig::lenient()
            .with_override(ErrorType::WrongDashCharacter, ErrorOverride::SilentCorrect);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(result.success);
        assert!(!result.has_warnings()); // Silent = no warnings
    }

    #[test]
    fn test_preprocessor_override_correct_in_strict() {
        let config = ErrorConfig::strict()
            .with_override(ErrorType::WrongDashCharacter, ErrorOverride::WarnCorrect);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("c.100\u{2013}200del");
        assert!(result.success);
        assert!(result.has_warnings());
        assert_eq!(result.preprocessed, "c.100-200del");
    }

    // CorrectionWarning tests
    #[test]
    fn test_correction_warning_from_correction() {
        let correction =
            DetectedCorrection::new(ErrorType::WrongDashCharacter, "\u{2013}", "-", 5, 8);
        let warning = CorrectionWarning::from_correction(&correction);
        assert_eq!(warning.error_type, ErrorType::WrongDashCharacter);
        assert!(warning.message.contains("dash"));
        assert_eq!(warning.span, Some((5, 8)));
    }

    // Trailing annotation tests
    #[test]
    fn test_preprocessor_strict_rejects_trailing_annotation() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NM_000088.3:c.459A>G (p.Lys153=)");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_preprocessor_lenient_strips_trailing_annotation() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_000088.3:c.459A>G (p.Lys153=)");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.459A>G");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_silent_strips_trailing_annotation() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("NM_000088.3:c.459A>G (p.Lys153=)");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.459A>G");
        assert!(!result.has_warnings());
    }

    #[test]
    fn test_preprocessor_lenient_clinvar_pattern() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_003467.3(CXCR4):c.708G>A (p.Lys236=)");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_003467.3(CXCR4):c.708G>A");
    }

    #[test]
    fn test_preprocessor_override_accept_trailing_annotation() {
        // Also need to override whitespace handling since there's a space before the annotation
        let config = ErrorConfig::strict()
            .with_override(ErrorType::ExtraWhitespace, ErrorOverride::SilentCorrect)
            .with_override(ErrorType::TrailingAnnotation, ErrorOverride::WarnCorrect);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("NM_000088.3:c.459A>G (p.Lys153=)");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.459A>G");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_override_trailing_annotation_no_space() {
        // Test without the space - only TrailingAnnotation override needed
        let config = ErrorConfig::strict()
            .with_override(ErrorType::TrailingAnnotation, ErrorOverride::WarnCorrect);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("NM_000088.3:c.459A>G(p.Lys153=)");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.459A>G");
        assert!(result.has_warnings());
    }

    // Missing coordinate prefix tests
    #[test]
    fn test_preprocessor_strict_rejects_missing_prefix() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NC_000017.11:12345A>G");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_preprocessor_lenient_adds_missing_prefix() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NC_000017.11:12345A>G");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NC_000017.11:g.12345A>G");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_lenient_adds_missing_prefix_uncertain() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NC_000017.11:(?_31094927)_(31377677_?)del");
        assert!(result.success);
        assert_eq!(
            result.preprocessed,
            "NC_000017.11:g.(?_31094927)_(31377677_?)del"
        );
    }

    #[test]
    fn test_preprocessor_silent_adds_missing_prefix() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("NC_000017.11:12345A>G");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NC_000017.11:g.12345A>G");
        assert!(!result.has_warnings());
    }
}
