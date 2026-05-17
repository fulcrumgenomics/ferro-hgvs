//! Input preprocessor for normalizing HGVS input strings.
//!
//! The preprocessor applies corrections to common input errors before
//! parsing, based on the configured error handling mode.

use super::corrections::{
    correct_accession_prefix_case, correct_amino_acid_case_in_protein, correct_dash_characters,
    correct_deprecated_con, correct_deprecated_protein_forms, correct_edit_type_case_full,
    correct_empty_delins, correct_missing_coordinate_prefix, correct_old_allele_format,
    correct_old_substitution_syntax, correct_protein_arrow, correct_quote_characters,
    correct_redundant_repeat_label, correct_single_letter_aa_in_protein,
    correct_single_position_range, correct_swapped_positions, correct_whitespace,
    detect_del_size_suffix, detect_deprecated_ivs, detect_missing_versions, detect_position_zero,
    strip_trailing_annotation, DetectedCorrection,
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

        // Phase 5a: Flag accessions that lack a `.<version>` suffix (W3001).
        // This is `warn_accept` per the registry — there is no canonical
        // version to inject, so lenient mode warns without mutating
        // `current`, silent mode accepts silently, and strict mode rejects.
        let missing_versions = detect_missing_versions(&current);
        if !missing_versions.is_empty() {
            let action = self.action_for(ErrorType::MissingVersion);
            match action {
                ResolvedAction::Reject => {
                    let first = &missing_versions[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Accession '{}' is missing a version suffix",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidAccession)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_hint(
                                    "RefSeq accessions require a `.<version>` suffix (e.g. NM_000088.3, NC_000023.11)",
                                ),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    // No correction is possible (we don't know the intended
                    // version); emit one warning per occurrence and leave
                    // `current` unchanged.
                    for c in &missing_versions {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                }
                ResolvedAction::SilentCorrect | ResolvedAction::Accept => {}
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

        // Phase 6a: Correct lowercase / mis-cased three-letter amino-acid
        // tokens within the protein description (W1001).
        let (corrected, corrections) = correct_amino_acid_case_in_protein(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::LowercaseAminoAcid);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Lowercase or mis-cased amino-acid code '{}', expected '{}'",
                                first.original, first.corrected
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidAminoAcid)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS prefers three-letter amino-acid codes with a leading capital (e.g. Val, Glu, Ala)",
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

        // Phase 6b: Rewrite deprecated stop-codon and frameshift forms in protein
        // descriptions (SVA-003..SVA-006, issue #125):
        //   p.Arg97*       → p.Arg97Ter         (W3007 DeprecatedStopCodonStar)
        //   p.Arg97X       → p.Arg97Ter         (W3008 DeprecatedStopCodonX)
        //   p.Arg97fs*23   → p.Arg97fsTer23     (W3009 DeprecatedFrameshiftStar)
        //   p.Arg97fsX23   → p.Arg97fsTer23     (W3010 DeprecatedFrameshiftX)
        //
        // Must run BEFORE Phase 6c (single-letter expansion): otherwise the
        // `X` in `p.Arg97X` would be expanded to `Xaa` ("any amino acid") and
        // the deprecated-stop-codon signal would be lost.
        //
        // Each detection's action is resolved independently and applied
        // per-correction, so mixing `Accept` with `WarnCorrect`/`SilentCorrect`
        // overrides across the four W-codes yields a partial rewrite —
        // Accept-marked tokens are preserved even when sibling detections in
        // the same input are rewritten.
        let (corrected_full, corrections) = correct_deprecated_protein_forms(&current);
        if !corrections.is_empty() {
            // Reject is sticky: any `Reject` detection fails the whole input,
            // and the suggestion shows the fully-rewritten canonical form.
            for c in &corrections {
                if matches!(self.action_for(c.error_type), ResolvedAction::Reject) {
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            c.start,
                            format!(
                                "Deprecated protein notation '{}', use '{}'",
                                c.original, c.corrected
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(c.start, c.end))
                                .with_source(input)
                                .with_suggestion(corrected_full.clone())
                                .with_hint(
                                    "HGVS uses 'Ter' (or 'fsTerN') for translation termination; \
                                    'X' and '*' are deprecated alternatives.",
                                ),
                        ),
                    );
                }
            }

            // Walk the original input and apply each non-Accept correction in
            // place. `corrections` is already in left-to-right byte order.
            let mut rebuilt = String::with_capacity(current.len());
            let mut cursor = 0usize;
            for c in &corrections {
                rebuilt.push_str(&current[cursor..c.start]);
                match self.action_for(c.error_type) {
                    ResolvedAction::Accept => {
                        rebuilt.push_str(&current[c.start..c.end]);
                    }
                    ResolvedAction::WarnCorrect => {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                        rebuilt.push_str(&c.corrected);
                    }
                    ResolvedAction::SilentCorrect => {
                        rebuilt.push_str(&c.corrected);
                    }
                    ResolvedAction::Reject => unreachable!("rejected above"),
                }
                cursor = c.end;
            }
            rebuilt.push_str(&current[cursor..]);
            current = rebuilt;
        }

        // Phase 6b1: Lowercase mixed-case edit-type tokens (W1004 — HGVS
        // spec recommendations/general.md lines 87-104). Must run BEFORE
        // Phase 6c (single-letter AA expansion) so a token like `DEL` in
        // a protein descriptor (`p.Arg8_Lys10DEL`) is lowercased to
        // `del` rather than mis-expanded to `AspGluLeu`. Also runs
        // before the edit-type-dependent NA phases (W4003 single-position
        // range, W3011 del-size, W3012 empty delins, W3013 redundant
        // repeat, W3015 con) so they see canonical lowercase tokens.
        let (corrected, corrections) = correct_edit_type_case_full(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::MixedCaseEditType);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!("Edit-type token '{}' must be lowercase", first.original),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS edit tokens (del, ins, dup, inv, delins, con) are spelled lowercase",
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

        // Phase 6c: Expand single-letter amino-acid codes to canonical
        // three-letter form within the protein description (W1002).
        let (corrected, corrections) = correct_single_letter_aa_in_protein(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::SingleLetterAminoAcid);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Single-letter amino-acid code '{}', expected three-letter '{}'",
                                first.original, first.corrected
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidAminoAcid)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS recommends three-letter amino-acid codes (e.g. Val, Glu, Ala) over one-letter abbreviations",
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

        // Phase 10: Collapse single-position ranges in del/dup/inv (W4003).
        // Per HGVS spec, `c.123_123del`, `c.123_123dup`, and `c.100_100inv`
        // should be written as the single-position form.
        let (corrected, corrections) = correct_single_position_range(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::SinglePositionRange);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Single-position range '{}' is non-canonical, expected '{}'",
                                first.original, first.corrected
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidPosition)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS recommends a single position for del/dup/inv when the range collapses to one base",
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

        // Phase 11: Rewrite `delins` with empty inserted sequence as `del` (W3012).
        let (corrected, corrections) = correct_empty_delins(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::EmptyDelinsInsert);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Deletion-insertion has empty inserted sequence",
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "An empty `delins` is semantically equivalent to a plain deletion (`del`)",
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

        // Phase 12: Strip redundant base labels in RNA repeat descriptions (W3013).
        let (corrected, corrections) = correct_redundant_repeat_label(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::RedundantRepeatLabel);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Repeat description has redundant base label '{}'",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "RNA repeat descriptions should omit the base label when positions already define the unit",
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

        // Phase 13: Flag deletions described with a size-count suffix (W3011).
        // `warn_accept` semantics: lenient warns without rewriting the input
        // (we cannot synthesise the end position safely; offset/intronic
        // semantics defeat naive `start + length - 1`), silent accepts
        // silently, strict rejects.
        let del_size_hits = detect_del_size_suffix(&current);
        if !del_size_hits.is_empty() {
            let action = self.action_for(ErrorType::DelSizeSuffix);
            match action {
                ResolvedAction::Reject => {
                    let first = &del_size_hits[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Deletion '{}' uses size-count suffix; canonical form names both endpoints",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_hint(
                                    "Write `g.<start>_<end>del` instead of `g.<start>del<size>`",
                                ),
                        ),
                    );
                }
                ResolvedAction::WarnCorrect => {
                    for c in &del_size_hits {
                        all_warnings.push(CorrectionWarning::from_correction(c));
                    }
                }
                ResolvedAction::SilentCorrect | ResolvedAction::Accept => {}
            }
        }

        // Phase 14: Rewrite deprecated multi-base substitution syntax to
        // delins (W3003 — HGVS spec recommendations/DNA/substitution.md line 26).
        let (corrected, corrections) = correct_old_substitution_syntax(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::OldSubstitutionSyntax);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Deprecated multi-base substitution syntax",
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS reserves '>' for single-base substitutions; use 'delins' for multi-base changes",
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

        // Phase 15: Rewrite deprecated `con` (sequence conversion) syntax to
        // delins (W3015 — HGVS spec recommendations/DNA/delins.md line 19).
        let (corrected, corrections) = correct_deprecated_con(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::DeprecatedConSyntax);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Deprecated 'con' (conversion) edit syntax",
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidEdit)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS retired the 'con' edit type; describe conversions as 'delins<source>'",
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

        // Phase 15a: Swap inverted-order interval positions (W4001 — HGVS spec
        // recommendations/general.md, range semantics start ≤ end). Covers
        // integer pairs, offset-bearing forms (`c.100+5_99+3del`), and 3'UTR
        // star positions (`c.*5_*1del`).
        let (corrected, corrections) = correct_swapped_positions(&current);
        if !corrections.is_empty() {
            let action = self.action_for(ErrorType::SwappedPositions);
            match action {
                ResolvedAction::Reject => {
                    let first = &corrections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            format!(
                                "Interval positions are swapped ('{}'); expected start ≤ end",
                                first.original
                            ),
                            Diagnostic::new()
                                .with_code(ErrorCode::InvalidPosition)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_suggestion(corrected.clone())
                                .with_hint(
                                    "HGVS range syntax requires start ≤ end; rewrite to the ascending form",
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

        // Phase 16: Reject retracted c.IVS intronic notation (W3014 — HGVS spec
        // background/numbering.md line 32). Cannot be auto-corrected without
        // genomic intron metadata, so the W3014 mode behavior is
        // always_reject; only ErrorOverride::Accept lets the input pass.
        let detections = detect_deprecated_ivs(&current);
        if !detections.is_empty() {
            let action = self.action_for(ErrorType::DeprecatedIvsNotation);
            match action {
                ResolvedAction::Reject
                | ResolvedAction::WarnCorrect
                | ResolvedAction::SilentCorrect => {
                    let first = &detections[0];
                    return PreprocessResult::failed(
                        input.to_string(),
                        FerroError::parse_with_diagnostic(
                            first.start,
                            "Retracted c.IVS intronic notation",
                            Diagnostic::new()
                                .with_code(ErrorCode::UnexpectedChar)
                                .with_span(SourceSpan::new(first.start, first.end))
                                .with_source(input)
                                .with_hint(
                                    "Use the canonical intronic-offset form (e.g. c.88+2T>G) — IVS notation has been retracted by HGVS and is ambiguous without genomic context",
                                ),
                        ),
                    );
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

    // ----------------------------------------------------------------------
    // Deprecated stop-codon and frameshift forms (issue #125, SVA-003..006).
    // ----------------------------------------------------------------------

    #[test]
    fn test_preprocessor_strict_rejects_deprecated_stop_star() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97*");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_preprocessor_strict_rejects_deprecated_stop_x() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97X");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_strict_rejects_deprecated_frameshift_star() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fs*23");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_strict_rejects_deprecated_frameshift_x() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fsX23");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_lenient_corrects_deprecated_stop_star() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97*");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97Ter");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(
            result.warnings[0].error_type,
            ErrorType::DeprecatedStopCodonStar
        );
    }

    #[test]
    fn test_preprocessor_lenient_corrects_deprecated_stop_x() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97X");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97Ter");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(
            result.warnings[0].error_type,
            ErrorType::DeprecatedStopCodonX
        );
    }

    #[test]
    fn test_preprocessor_lenient_corrects_deprecated_frameshift_star() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fs*23");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97fsTer23");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(
            result.warnings[0].error_type,
            ErrorType::DeprecatedFrameshiftStar
        );
    }

    #[test]
    fn test_preprocessor_lenient_corrects_deprecated_frameshift_x() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fsX23");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97fsTer23");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(
            result.warnings[0].error_type,
            ErrorType::DeprecatedFrameshiftX
        );
    }

    #[test]
    fn test_preprocessor_silent_corrects_deprecated_without_warnings() {
        let preprocessor = InputPreprocessor::silent();
        for input in [
            "NP_000079.2:p.Arg97*",
            "NP_000079.2:p.Arg97X",
            "NP_000079.2:p.Arg97fs*23",
            "NP_000079.2:p.Arg97fsX23",
        ] {
            let result = preprocessor.preprocess(input);
            assert!(result.success, "expected success for {}", input);
            assert!(!result.has_warnings(), "expected no warnings for {}", input);
            assert!(
                result.preprocessed.contains("Ter"),
                "expected Ter in {}",
                result.preprocessed
            );
        }
    }

    #[test]
    fn test_preprocessor_lenient_canonical_no_warnings() {
        let preprocessor = InputPreprocessor::lenient();
        for input in [
            "NP_000079.2:p.Arg97Ter",
            "NP_000079.2:p.Arg97ProfsTer23",
            "NP_000079.2:p.Tyr180fs",
            "NP_000079.2:p.Val600Glu",
            "NP_000079.2:p.Arg782Xaa",
        ] {
            let result = preprocessor.preprocess(input);
            assert!(result.success, "expected success for {}", input);
            assert_eq!(
                result.preprocessed, input,
                "expected unchanged for {}",
                input
            );
            assert!(
                !result.has_warnings(),
                "expected no warnings for {}, got {:?}",
                input,
                result.warnings
            );
        }
    }

    #[test]
    fn test_preprocessor_lenient_compound_protein_allele_two_warnings() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NP_000079.2:p.[Arg97*;Arg100X]");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.[Arg97Ter;Arg100Ter]");
        assert_eq!(result.warnings.len(), 2);
        assert_eq!(
            result.warnings[0].error_type,
            ErrorType::DeprecatedStopCodonStar
        );
        assert_eq!(
            result.warnings[1].error_type,
            ErrorType::DeprecatedStopCodonX
        );
    }

    #[test]
    fn test_preprocessor_lenient_idempotent_on_corrected_output() {
        // Re-running the preprocessor on its own output yields no further
        // deprecated-form warnings.
        let preprocessor = InputPreprocessor::lenient();
        let first = preprocessor.preprocess("NP_000079.2:p.Arg97fs*23");
        let second = preprocessor.preprocess(&first.preprocessed);
        assert_eq!(second.preprocessed, first.preprocessed);
        assert!(!second.has_warnings());
    }

    #[test]
    fn test_preprocessor_override_accept_keeps_deprecated_form() {
        // When the override is Accept, the deprecated form passes through
        // unchanged with no warning.
        let config = ErrorConfig::lenient()
            .with_override(ErrorType::DeprecatedStopCodonStar, ErrorOverride::Accept);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97*");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97*");
        assert!(!result.has_warnings());
    }

    #[test]
    fn test_preprocessor_override_silent_in_strict_mode() {
        // SilentCorrect override in strict mode rewrites without warning.
        let config = ErrorConfig::strict().with_override(
            ErrorType::DeprecatedFrameshiftStar,
            ErrorOverride::SilentCorrect,
        );
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fs*23");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97fsTer23");
        assert!(!result.has_warnings());
    }

    #[test]
    fn test_preprocessor_lenient_does_not_affect_cds_utr_position() {
        // c.*5A>G is a 3'UTR position, not a deprecated stop codon. Must NOT
        // be rewritten.
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_000088.3:c.*5A>G");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.*5A>G");
        assert!(!result.has_warnings());
    }

    #[test]
    fn test_preprocessor_partial_accept_only_rewrites_non_accept_codes() {
        // With DeprecatedStopCodonStar = Accept and DeprecatedStopCodonX left
        // at the lenient default (WarnCorrect), the `*` must remain literal
        // while the `X` is rewritten to `Ter`. Per-correction action
        // resolution — without it, the Accept override is silently ignored
        // when any sibling detection is non-Accept.
        let config = ErrorConfig::lenient()
            .with_override(ErrorType::DeprecatedStopCodonStar, ErrorOverride::Accept);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("NP_000079.2:p.[Arg97*;Arg100X]");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.[Arg97*;Arg100Ter]");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(
            result.warnings[0].error_type,
            ErrorType::DeprecatedStopCodonX
        );
    }

    // ===== Issue #115: deprecated multi-base substitution (W3003) =====

    #[test]
    fn test_preprocessor_strict_rejects_old_substitution_with_refs() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NM_000088.3:c.79_80GC>TT");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_preprocessor_lenient_corrects_old_substitution_with_refs() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_000088.3:c.79_80GC>TT");
        assert!(result.success, "lenient should rewrite to delins");
        assert_eq!(result.preprocessed, "NM_000088.3:c.79_80delinsTT");
        assert!(result.has_warnings());
        assert!(result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::OldSubstitutionSyntax));
    }

    #[test]
    fn test_preprocessor_lenient_corrects_old_substitution_no_refs() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_000088.3:c.100_102>ATG");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.100_102delinsATG");
        assert!(result.has_warnings());
    }

    #[test]
    fn test_preprocessor_silent_rewrites_old_substitution_no_warning() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("NM_000088.3:c.79_80GC>TT");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.79_80delinsTT");
        assert!(!result.has_warnings());
    }

    #[test]
    fn test_preprocessor_canonical_substitution_unchanged() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NM_000088.3:c.100A>G");
        assert!(result.success);
        assert!(!result.has_corrections());
    }

    // ===== Issue #115: deprecated `con` syntax (W3015) =====

    #[test]
    fn test_preprocessor_strict_rejects_con_syntax() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NM_004006.2:c.100_200conNM_001.1:c.5_105");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_lenient_corrects_con_to_delins() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_004006.2:c.100_200conNM_001.1:c.5_105");
        assert!(result.success);
        assert_eq!(
            result.preprocessed,
            "NM_004006.2:c.100_200delinsNM_001.1:c.5_105"
        );
        assert!(result.has_warnings());
        assert!(result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::DeprecatedConSyntax));
    }

    #[test]
    fn test_preprocessor_silent_corrects_con_no_warning() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("NM_004006.2:c.100_200conNM_001.1:c.5_105");
        assert!(result.success);
        assert_eq!(
            result.preprocessed,
            "NM_004006.2:c.100_200delinsNM_001.1:c.5_105"
        );
        assert!(!result.has_warnings());
    }

    // ===== Issue #115: retracted c.IVS notation (W3014) =====

    #[test]
    fn test_preprocessor_strict_rejects_ivs_notation() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NM_000088.3:c.IVS2+2T>G");
        assert!(!result.success);
        assert!(result.error.is_some());
    }

    #[test]
    fn test_preprocessor_lenient_rejects_ivs_notation() {
        let preprocessor = InputPreprocessor::lenient();
        let result = preprocessor.preprocess("NM_000088.3:c.IVS2+2T>G");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_silent_rejects_ivs_notation() {
        let preprocessor = InputPreprocessor::silent();
        let result = preprocessor.preprocess("NM_000088.3:c.IVS2+2T>G");
        assert!(!result.success);
    }

    #[test]
    fn test_preprocessor_override_accept_keeps_ivs_notation() {
        let config = ErrorConfig::strict()
            .with_override(ErrorType::DeprecatedIvsNotation, ErrorOverride::Accept);
        let preprocessor = InputPreprocessor::new(config);
        let result = preprocessor.preprocess("NM_000088.3:c.IVS2+2T>G");
        // With Accept override, the preprocessor leaves the input alone; the
        // downstream parser will of course still fail — but that's the
        // user's choice.
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.IVS2+2T>G");
    }

    #[test]
    fn test_preprocessor_canonical_intronic_unchanged() {
        let preprocessor = InputPreprocessor::strict();
        let result = preprocessor.preprocess("NM_000088.3:c.88+2T>G");
        assert!(result.success);
        assert!(!result.has_corrections());
    }
}
