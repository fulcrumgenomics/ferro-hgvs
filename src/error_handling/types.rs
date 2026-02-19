//! Core types for error handling modes.
//!
//! This module defines the error mode, error type, and override enums
//! used for configurable error handling.

use std::fmt;

/// Error handling mode.
///
/// Controls how the parser handles common input errors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ErrorMode {
    /// Reject all non-standard input.
    ///
    /// Any deviation from strict HGVS syntax will result in an error.
    /// This is the default mode for maximum compliance.
    #[default]
    Strict,

    /// Auto-correct common errors with warnings.
    ///
    /// Common errors like wrong dash characters or lowercase amino acids
    /// will be automatically corrected, but warnings will be generated.
    Lenient,

    /// Auto-correct silently without warnings.
    ///
    /// Common errors are automatically corrected without generating
    /// any warnings. Use this for batch processing where you want
    /// maximum tolerance.
    Silent,
}

impl ErrorMode {
    /// Returns true if this mode should reject non-standard input.
    pub fn is_strict(&self) -> bool {
        matches!(self, ErrorMode::Strict)
    }

    /// Returns true if this mode allows auto-correction.
    pub fn allows_correction(&self) -> bool {
        matches!(self, ErrorMode::Lenient | ErrorMode::Silent)
    }

    /// Returns true if this mode should emit warnings.
    pub fn emits_warnings(&self) -> bool {
        matches!(self, ErrorMode::Lenient)
    }
}

impl fmt::Display for ErrorMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ErrorMode::Strict => write!(f, "strict"),
            ErrorMode::Lenient => write!(f, "lenient"),
            ErrorMode::Silent => write!(f, "silent"),
        }
    }
}

/// Individual error type that can be configured.
///
/// Each variant represents a specific type of input error that can
/// have custom handling behavior.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ErrorType {
    /// Lowercase amino acid codes (e.g., `val` instead of `Val`).
    LowercaseAminoAcid,

    /// Missing accession version (e.g., `NM_000088` instead of `NM_000088.3`).
    MissingVersion,

    /// Wrong dash character (en-dash/em-dash instead of hyphen).
    WrongDashCharacter,

    /// Extra whitespace in the description.
    ExtraWhitespace,

    /// Arrow in protein substitution (e.g., `Val600>Glu` instead of `Val600Glu`).
    ProteinSubstitutionArrow,

    /// Invalid position zero (position 0 is never valid in HGVS).
    PositionZero,

    /// Single-letter amino acid codes (e.g., `V` instead of `Val`).
    SingleLetterAminoAcid,

    /// Wrong quote characters (smart quotes instead of regular quotes).
    WrongQuoteCharacter,

    /// Lowercase reference sequence type (e.g., `nm_` instead of `NM_`).
    LowercaseAccessionPrefix,

    /// Mixed case in edit type (e.g., `Del` instead of `del`).
    MixedCaseEditType,

    /// Old-style substitution syntax (`>` for multi-base substitution).
    OldSubstitutionSyntax,

    /// Invalid Unicode characters that should be ASCII.
    InvalidUnicodeCharacter,

    /// Swapped interval positions (e.g., `c.200_100del` instead of `c.100_200del`).
    SwappedPositions,

    /// Trailing protein annotation (e.g., `(p.Lys236=)` after a coding variant).
    ///
    /// ClinVar and other databases often append protein consequence annotations
    /// to HGVS expressions. These are not valid HGVS syntax but are common in
    /// real-world data.
    TrailingAnnotation,

    /// Missing coordinate type prefix (e.g., `NC_000001.11:12345A>G` instead of `NC_000001.11:g.12345A>G`).
    ///
    /// Some sources omit the coordinate type prefix (g., c., etc.). For chromosomal
    /// accessions (NC_), the coordinate type can be inferred as genomic (g.).
    MissingCoordinatePrefix,

    /// Old/deprecated allele format with coordinate type inside brackets.
    ///
    /// Old format: `NM_000088.3:[c.100A>G;c.200C>T]` (coordinate type repeated inside brackets)
    /// Correct format: `NM_000088.3:c.[100A>G;200C>T]` (coordinate type before brackets)
    ///
    /// This format was used in older databases but is not valid per current HGVS spec.
    OldAlleleFormat,

    // === Normalization Error Types ===
    /// Reference sequence mismatch during normalization.
    ///
    /// The reference base(s) stated in the HGVS expression do not match
    /// the actual reference sequence. For example, `c.100G>A` but the
    /// reference at position 100 is actually T.
    ///
    /// In strict mode, this is rejected. In lenient/silent mode, normalization
    /// proceeds using the actual reference sequence.
    RefSeqMismatch,
}

impl ErrorType {
    /// Returns the warning code (W-code) for this error type.
    pub fn code(&self) -> &'static str {
        match self {
            ErrorType::LowercaseAminoAcid => "W1001",
            ErrorType::SingleLetterAminoAcid => "W1002",
            ErrorType::LowercaseAccessionPrefix => "W1003",
            ErrorType::MixedCaseEditType => "W1004",
            ErrorType::WrongDashCharacter => "W2001",
            ErrorType::WrongQuoteCharacter => "W2002",
            ErrorType::ExtraWhitespace => "W2003",
            ErrorType::InvalidUnicodeCharacter => "W2004",
            ErrorType::MissingVersion => "W3001",
            ErrorType::ProteinSubstitutionArrow => "W3002",
            ErrorType::OldSubstitutionSyntax => "W3003",
            ErrorType::OldAlleleFormat => "W3004",
            ErrorType::TrailingAnnotation => "W3005",
            ErrorType::MissingCoordinatePrefix => "W3006",
            ErrorType::SwappedPositions => "W4001",
            ErrorType::PositionZero => "W4002",
            ErrorType::RefSeqMismatch => "W5001",
        }
    }

    /// Returns a human-readable description of this error type.
    pub fn description(&self) -> &'static str {
        match self {
            ErrorType::LowercaseAminoAcid => "lowercase amino acid code",
            ErrorType::MissingVersion => "missing accession version",
            ErrorType::WrongDashCharacter => "wrong dash character (en-dash or em-dash)",
            ErrorType::ExtraWhitespace => "extra whitespace in description",
            ErrorType::ProteinSubstitutionArrow => "arrow in protein substitution",
            ErrorType::PositionZero => "invalid position zero",
            ErrorType::SingleLetterAminoAcid => "single-letter amino acid code",
            ErrorType::WrongQuoteCharacter => "wrong quote character (smart quotes)",
            ErrorType::LowercaseAccessionPrefix => "lowercase accession prefix",
            ErrorType::MixedCaseEditType => "mixed case in edit type",
            ErrorType::OldSubstitutionSyntax => "old-style substitution syntax",
            ErrorType::InvalidUnicodeCharacter => "invalid Unicode character",
            ErrorType::SwappedPositions => "swapped interval positions",
            ErrorType::TrailingAnnotation => "trailing protein annotation",
            ErrorType::MissingCoordinatePrefix => "missing coordinate type prefix",
            ErrorType::OldAlleleFormat => "old/deprecated allele format",
            ErrorType::RefSeqMismatch => "reference sequence mismatch",
        }
    }

    /// Returns true if this error type can be auto-corrected.
    pub fn is_correctable(&self) -> bool {
        match self {
            // Position zero is never valid and cannot be auto-corrected
            ErrorType::PositionZero => false,
            // RefSeqMismatch can be "corrected" by using the actual reference
            ErrorType::RefSeqMismatch => true,
            // All other error types can potentially be auto-corrected
            _ => true,
        }
    }

    /// Returns an example of this error type for documentation.
    pub fn example(&self) -> (&'static str, &'static str) {
        match self {
            ErrorType::LowercaseAminoAcid => ("val600Glu", "Val600Glu"),
            ErrorType::MissingVersion => ("NM_000088:c.100A>G", "NM_000088.3:c.100A>G"),
            ErrorType::WrongDashCharacter => ("c.100–200del", "c.100-200del"),
            ErrorType::ExtraWhitespace => ("c.100 A>G", "c.100A>G"),
            ErrorType::ProteinSubstitutionArrow => ("p.Val600>Glu", "p.Val600Glu"),
            ErrorType::PositionZero => ("c.0A>G", "(invalid)"),
            ErrorType::SingleLetterAminoAcid => ("p.V600E", "p.Val600Glu"),
            ErrorType::WrongQuoteCharacter => {
                ("c.100_101ins\u{201C}ATG\u{201D}", "c.100_101ins\"ATG\"")
            }
            ErrorType::LowercaseAccessionPrefix => ("nm_000088.3", "NM_000088.3"),
            ErrorType::MixedCaseEditType => ("c.100Del", "c.100del"),
            ErrorType::OldSubstitutionSyntax => ("c.100_102>ATG", "c.100_102delinsATG"),
            ErrorType::InvalidUnicodeCharacter => ("c.100\u{2019}200del", "c.100-200del"),
            ErrorType::SwappedPositions => ("c.200_100del", "c.100_200del"),
            ErrorType::TrailingAnnotation => {
                ("NM_000088.3:c.459A>G (p.Lys153=)", "NM_000088.3:c.459A>G")
            }
            ErrorType::MissingCoordinatePrefix => {
                ("NC_000001.11:12345A>G", "NC_000001.11:g.12345A>G")
            }
            ErrorType::OldAlleleFormat => (
                "NM_000088.3:[c.100A>G;c.200C>T]",
                "NM_000088.3:c.[100A>G;200C>T]",
            ),
            ErrorType::RefSeqMismatch => ("c.100G>A (ref is T)", "c.100T>A (corrected)"),
        }
    }
}

impl fmt::Display for ErrorType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.description())
    }
}

/// Per-error-type override behavior.
///
/// Allows overriding the default mode behavior for specific error types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum ErrorOverride {
    /// Use the default mode behavior.
    #[default]
    Default,

    /// Always reject this error type.
    Reject,

    /// Always warn and correct.
    WarnCorrect,

    /// Always correct silently.
    SilentCorrect,

    /// Always accept as-is without correction.
    Accept,
}

impl ErrorOverride {
    /// Resolve this override to an effective action given the base mode.
    pub fn resolve(&self, mode: ErrorMode) -> ResolvedAction {
        match self {
            ErrorOverride::Default => match mode {
                ErrorMode::Strict => ResolvedAction::Reject,
                ErrorMode::Lenient => ResolvedAction::WarnCorrect,
                ErrorMode::Silent => ResolvedAction::SilentCorrect,
            },
            ErrorOverride::Reject => ResolvedAction::Reject,
            ErrorOverride::WarnCorrect => ResolvedAction::WarnCorrect,
            ErrorOverride::SilentCorrect => ResolvedAction::SilentCorrect,
            ErrorOverride::Accept => ResolvedAction::Accept,
        }
    }
}

impl fmt::Display for ErrorOverride {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ErrorOverride::Default => write!(f, "default"),
            ErrorOverride::Reject => write!(f, "reject"),
            ErrorOverride::WarnCorrect => write!(f, "warn+correct"),
            ErrorOverride::SilentCorrect => write!(f, "silent correct"),
            ErrorOverride::Accept => write!(f, "accept"),
        }
    }
}

/// Resolved action after applying mode and override.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolvedAction {
    /// Reject the input with an error.
    Reject,
    /// Correct the input and emit a warning.
    WarnCorrect,
    /// Correct the input silently.
    SilentCorrect,
    /// Accept the input as-is.
    Accept,
}

impl ResolvedAction {
    /// Returns true if this action should reject the input.
    pub fn should_reject(&self) -> bool {
        matches!(self, ResolvedAction::Reject)
    }

    /// Returns true if this action should correct the input.
    pub fn should_correct(&self) -> bool {
        matches!(
            self,
            ResolvedAction::WarnCorrect | ResolvedAction::SilentCorrect
        )
    }

    /// Returns true if this action should emit a warning.
    pub fn should_warn(&self) -> bool {
        matches!(self, ResolvedAction::WarnCorrect)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ErrorMode tests
    #[test]
    fn test_error_mode_default() {
        assert_eq!(ErrorMode::default(), ErrorMode::Strict);
    }

    #[test]
    fn test_error_mode_is_strict() {
        assert!(ErrorMode::Strict.is_strict());
        assert!(!ErrorMode::Lenient.is_strict());
        assert!(!ErrorMode::Silent.is_strict());
    }

    #[test]
    fn test_error_mode_allows_correction() {
        assert!(!ErrorMode::Strict.allows_correction());
        assert!(ErrorMode::Lenient.allows_correction());
        assert!(ErrorMode::Silent.allows_correction());
    }

    #[test]
    fn test_error_mode_emits_warnings() {
        assert!(!ErrorMode::Strict.emits_warnings());
        assert!(ErrorMode::Lenient.emits_warnings());
        assert!(!ErrorMode::Silent.emits_warnings());
    }

    #[test]
    fn test_error_mode_display() {
        assert_eq!(format!("{}", ErrorMode::Strict), "strict");
        assert_eq!(format!("{}", ErrorMode::Lenient), "lenient");
        assert_eq!(format!("{}", ErrorMode::Silent), "silent");
    }

    // ErrorType tests
    #[test]
    fn test_error_type_description() {
        assert_eq!(
            ErrorType::LowercaseAminoAcid.description(),
            "lowercase amino acid code"
        );
        assert_eq!(
            ErrorType::MissingVersion.description(),
            "missing accession version"
        );
        assert_eq!(
            ErrorType::WrongDashCharacter.description(),
            "wrong dash character (en-dash or em-dash)"
        );
    }

    #[test]
    fn test_error_type_is_correctable() {
        assert!(ErrorType::LowercaseAminoAcid.is_correctable());
        assert!(ErrorType::WrongDashCharacter.is_correctable());
        assert!(!ErrorType::PositionZero.is_correctable());
        // RefSeqMismatch is correctable (use actual reference)
        assert!(ErrorType::RefSeqMismatch.is_correctable());
    }

    #[test]
    fn test_ref_seq_mismatch_error_type() {
        let err = ErrorType::RefSeqMismatch;
        assert_eq!(err.description(), "reference sequence mismatch");
        assert!(err.is_correctable());
        let (bad, good) = err.example();
        assert!(bad.contains("ref is T"));
        assert!(good.contains("corrected"));
    }

    #[test]
    fn test_error_type_example() {
        let (bad, good) = ErrorType::LowercaseAminoAcid.example();
        assert_eq!(bad, "val600Glu");
        assert_eq!(good, "Val600Glu");
    }

    #[test]
    fn test_error_type_display() {
        assert_eq!(
            format!("{}", ErrorType::LowercaseAminoAcid),
            "lowercase amino acid code"
        );
    }

    #[test]
    fn test_error_type_code() {
        assert_eq!(ErrorType::LowercaseAminoAcid.code(), "W1001");
        assert_eq!(ErrorType::SingleLetterAminoAcid.code(), "W1002");
        assert_eq!(ErrorType::WrongDashCharacter.code(), "W2001");
        assert_eq!(ErrorType::MissingVersion.code(), "W3001");
        assert_eq!(ErrorType::SwappedPositions.code(), "W4001");
        assert_eq!(ErrorType::PositionZero.code(), "W4002");
        assert_eq!(ErrorType::RefSeqMismatch.code(), "W5001");
    }

    // ErrorOverride tests
    #[test]
    fn test_error_override_default() {
        assert_eq!(ErrorOverride::default(), ErrorOverride::Default);
    }

    #[test]
    fn test_error_override_resolve_default() {
        // Default with Strict mode
        assert_eq!(
            ErrorOverride::Default.resolve(ErrorMode::Strict),
            ResolvedAction::Reject
        );

        // Default with Lenient mode
        assert_eq!(
            ErrorOverride::Default.resolve(ErrorMode::Lenient),
            ResolvedAction::WarnCorrect
        );

        // Default with Silent mode
        assert_eq!(
            ErrorOverride::Default.resolve(ErrorMode::Silent),
            ResolvedAction::SilentCorrect
        );
    }

    #[test]
    fn test_error_override_resolve_explicit() {
        // Explicit overrides should ignore the mode
        assert_eq!(
            ErrorOverride::Reject.resolve(ErrorMode::Silent),
            ResolvedAction::Reject
        );

        assert_eq!(
            ErrorOverride::WarnCorrect.resolve(ErrorMode::Strict),
            ResolvedAction::WarnCorrect
        );

        assert_eq!(
            ErrorOverride::SilentCorrect.resolve(ErrorMode::Strict),
            ResolvedAction::SilentCorrect
        );

        assert_eq!(
            ErrorOverride::Accept.resolve(ErrorMode::Strict),
            ResolvedAction::Accept
        );
    }

    #[test]
    fn test_error_override_display() {
        assert_eq!(format!("{}", ErrorOverride::Default), "default");
        assert_eq!(format!("{}", ErrorOverride::Reject), "reject");
        assert_eq!(format!("{}", ErrorOverride::WarnCorrect), "warn+correct");
        assert_eq!(
            format!("{}", ErrorOverride::SilentCorrect),
            "silent correct"
        );
        assert_eq!(format!("{}", ErrorOverride::Accept), "accept");
    }

    // ResolvedAction tests
    #[test]
    fn test_resolved_action_should_reject() {
        assert!(ResolvedAction::Reject.should_reject());
        assert!(!ResolvedAction::WarnCorrect.should_reject());
        assert!(!ResolvedAction::SilentCorrect.should_reject());
        assert!(!ResolvedAction::Accept.should_reject());
    }

    #[test]
    fn test_resolved_action_should_correct() {
        assert!(!ResolvedAction::Reject.should_correct());
        assert!(ResolvedAction::WarnCorrect.should_correct());
        assert!(ResolvedAction::SilentCorrect.should_correct());
        assert!(!ResolvedAction::Accept.should_correct());
    }

    #[test]
    fn test_resolved_action_should_warn() {
        assert!(!ResolvedAction::Reject.should_warn());
        assert!(ResolvedAction::WarnCorrect.should_warn());
        assert!(!ResolvedAction::SilentCorrect.should_warn());
        assert!(!ResolvedAction::Accept.should_warn());
    }
}
