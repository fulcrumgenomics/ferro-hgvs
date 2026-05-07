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

    /// Deprecated `*` for stop codon in protein substitution position.
    ///
    /// Per HGVS checklist (`recommendations/checklist.md`), `Ter` and `*` may
    /// both indicate a translation stop codon, but the three-letter `Ter` is
    /// preferred. Example: `p.Arg97*` (corrected to `p.Arg97Ter`).
    DeprecatedStopCodonStar,

    /// Deprecated `X` for stop codon in protein substitution position.
    ///
    /// Per HGVS checklist (`recommendations/checklist.md`), "the X should not
    /// be used" to indicate a translation stop codon. The single-letter `X`
    /// is also reserved for the "any amino acid" symbol Xaa
    /// (`recommendations/protein/substitution.md`), making `p.Arg97X` ambiguous.
    /// Corrected to `p.Arg97Ter`.
    DeprecatedStopCodonX,

    /// Deprecated `fs*N` frameshift termination notation.
    ///
    /// Per HGVS frameshift recommendation (`recommendations/protein/frameshift.md`),
    /// the canonical form is `fsTerN` (e.g. `p.Arg123LysfsTer34`). The `fs*N`
    /// form is permitted as an alternative but `fsTerN` is preferred.
    /// Example: `p.Arg97fs*23` (corrected to `p.Arg97fsTer23`).
    DeprecatedFrameshiftStar,

    /// Deprecated `fsXN` frameshift termination notation.
    ///
    /// Per HGVS frameshift recommendation, neither `X` nor `Xaa` is used in
    /// frameshift termination; the canonical form is `fsTerN`. Example:
    /// `p.Arg97fsX23` (corrected to `p.Arg97fsTer23`).
    DeprecatedFrameshiftX,

    /// Deletion described with a size-count suffix instead of a position range
    /// (e.g. `g.123del6` instead of `g.123_128del`).
    ///
    /// Per HGVS spec (`recommendations/DNA/deletion.md`): "a deletion of more
    /// than one residue should mention the first and last residue deleted,
    /// e.g. `NG_012232.1:g.123_128del` and not `NG_012232.1:g.123del6`."
    DelSizeSuffix,

    /// Deletion-insertion with an empty inserted sequence (e.g.
    /// `g.100_102delins` instead of `g.100_102del`).
    ///
    /// A `delins` whose inserted sequence is absent collapses semantically to
    /// a plain deletion. Per `recommendations/DNA/delins.md` and the #81 A3
    /// canonicalization rule, the canonical form is `del`.
    EmptyDelinsInsert,

    /// Repeat described with a redundant base label inside the brackets
    /// (e.g. `r.100_102cug[4]` instead of `r.100_102[4]`).
    ///
    /// Per HGVS spec (`recommendations/RNA/repeated.md`): "the format
    /// `r.-125_-123cug[4]`, should not be used; it contains redundant
    /// information ('-125_-123' and 'cug')."
    RedundantRepeatLabel,

    /// Single-position range used where a single position is canonical (e.g.
    /// `c.123_123del`, `c.123_123dup`, or `c.100_100inv`).
    ///
    /// Per HGVS spec (`recommendations/DNA/deletion.md`,
    /// `recommendations/DNA/duplication.md`,
    /// `recommendations/DNA/inversion.md`): "position(s)_deleted should
    /// contain two different positions, e.g. 123_126 not 123_123" — the same
    /// principle applies to duplications and inversions.
    SinglePositionRange,

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

    /// Retracted `c.IVS` intronic notation (e.g., `c.IVS2+2T>G`).
    ///
    /// The `c.IVSn+offset` / `c.IVSn-offset` form was retracted by HGVS;
    /// canonical intronic offsets `c.<exon-pos>+<offset>` should be used
    /// instead. The retracted form cannot be auto-rewritten without genomic
    /// intron metadata (the IVS number is not unique across transcripts), so
    /// all modes reject it with an actionable hint.
    DeprecatedIvsNotation,

    /// Deprecated `con` (sequence conversion) edit syntax.
    ///
    /// The HGVS spec retired `con` in favour of `delins`. ferro can rewrite
    /// `c.100_200conNM_001:c.5_105` → `c.100_200delinsNM_001:c.5_105` in
    /// lenient/silent modes; strict mode rejects the deprecated form.
    DeprecatedConSyntax,
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
            ErrorType::DeprecatedStopCodonStar => "W3007",
            ErrorType::DeprecatedStopCodonX => "W3008",
            ErrorType::DeprecatedFrameshiftStar => "W3009",
            ErrorType::DeprecatedFrameshiftX => "W3010",
            ErrorType::DelSizeSuffix => "W3011",
            ErrorType::EmptyDelinsInsert => "W3012",
            ErrorType::RedundantRepeatLabel => "W3013",
            ErrorType::DeprecatedIvsNotation => "W3014",
            ErrorType::DeprecatedConSyntax => "W3015",
            ErrorType::SwappedPositions => "W4001",
            ErrorType::PositionZero => "W4002",
            ErrorType::SinglePositionRange => "W4003",
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
            ErrorType::DeprecatedStopCodonStar => "deprecated '*' for stop codon (use 'Ter')",
            ErrorType::DeprecatedStopCodonX => "deprecated 'X' for stop codon (use 'Ter')",
            ErrorType::DeprecatedFrameshiftStar => {
                "deprecated 'fs*N' frameshift notation (use 'fsTerN')"
            }
            ErrorType::DeprecatedFrameshiftX => {
                "deprecated 'fsXN' frameshift notation (use 'fsTerN')"
            }
            ErrorType::DelSizeSuffix => "deletion described with a size-count suffix",
            ErrorType::EmptyDelinsInsert => "delins with empty inserted sequence",
            ErrorType::RedundantRepeatLabel => "repeat description with redundant base label",
            ErrorType::SinglePositionRange => {
                "single-position range used where a single position is canonical"
            }
            ErrorType::RefSeqMismatch => "reference sequence mismatch",
            ErrorType::DeprecatedIvsNotation => "retracted c.IVS intronic notation",
            ErrorType::DeprecatedConSyntax => "deprecated con (conversion) edit syntax",
        }
    }

    /// Returns true if this error type can be auto-corrected.
    ///
    /// Explicit arms per variant — adding a new `ErrorType` triggers an
    /// exhaustiveness error here so the correctability decision is forced
    /// at the source rather than silently falling through to `true`.
    pub fn is_correctable(&self) -> bool {
        match self {
            ErrorType::LowercaseAminoAcid => true,
            ErrorType::MissingVersion => true,
            ErrorType::WrongDashCharacter => true,
            ErrorType::ExtraWhitespace => true,
            ErrorType::ProteinSubstitutionArrow => true,
            // Position zero is never valid and cannot be auto-corrected
            ErrorType::PositionZero => false,
            ErrorType::SingleLetterAminoAcid => true,
            ErrorType::WrongQuoteCharacter => true,
            ErrorType::LowercaseAccessionPrefix => true,
            ErrorType::MixedCaseEditType => true,
            ErrorType::OldSubstitutionSyntax => true,
            ErrorType::InvalidUnicodeCharacter => true,
            ErrorType::SwappedPositions => true,
            ErrorType::TrailingAnnotation => true,
            ErrorType::MissingCoordinatePrefix => true,
            ErrorType::OldAlleleFormat => true,
            // RefSeqMismatch can be "corrected" by using the actual reference
            ErrorType::RefSeqMismatch => true,
            ErrorType::DeprecatedStopCodonStar => true,
            ErrorType::DeprecatedStopCodonX => true,
            ErrorType::DeprecatedFrameshiftStar => true,
            ErrorType::DeprecatedFrameshiftX => true,
            // Lenient mode warns without rewriting (warn_accept); strict rejects.
            ErrorType::DelSizeSuffix => false,
            ErrorType::EmptyDelinsInsert => true,
            ErrorType::RedundantRepeatLabel => true,
            ErrorType::SinglePositionRange => true,
            // IVS notation cannot be auto-rewritten without genomic metadata
            ErrorType::DeprecatedIvsNotation => false,
            // `con` is rewritten to `delins` per SVD-WG009
            ErrorType::DeprecatedConSyntax => true,
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
            ErrorType::DeprecatedStopCodonStar => ("p.Arg97*", "p.Arg97Ter"),
            ErrorType::DeprecatedStopCodonX => ("p.Arg97X", "p.Arg97Ter"),
            ErrorType::DeprecatedFrameshiftStar => ("p.Arg97fs*23", "p.Arg97fsTer23"),
            ErrorType::DeprecatedFrameshiftX => ("p.Arg97fsX23", "p.Arg97fsTer23"),
            ErrorType::DelSizeSuffix => ("g.123del6", "g.123_128del"),
            ErrorType::EmptyDelinsInsert => ("g.100_102delins", "g.100_102del"),
            ErrorType::RedundantRepeatLabel => ("r.100_102cug[4]", "r.100_102[4]"),
            ErrorType::SinglePositionRange => ("c.123_123del", "c.123del"),
            ErrorType::RefSeqMismatch => ("c.100G>A (ref is T)", "c.100T>A (corrected)"),
            ErrorType::DeprecatedIvsNotation => ("c.IVS2+2T>G", "c.88+2T>G"),
            ErrorType::DeprecatedConSyntax => (
                "c.100_200conNM_001.1:c.5_105",
                "c.100_200delinsNM_001.1:c.5_105",
            ),
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
        // Pin the correctability decision per variant. Adding a new variant
        // forces an explicit choice in is_correctable() and an entry here.
        assert!(ErrorType::LowercaseAminoAcid.is_correctable());
        assert!(ErrorType::MissingVersion.is_correctable());
        assert!(ErrorType::WrongDashCharacter.is_correctable());
        assert!(ErrorType::ExtraWhitespace.is_correctable());
        assert!(ErrorType::ProteinSubstitutionArrow.is_correctable());
        assert!(!ErrorType::PositionZero.is_correctable());
        assert!(ErrorType::SingleLetterAminoAcid.is_correctable());
        assert!(ErrorType::WrongQuoteCharacter.is_correctable());
        assert!(ErrorType::LowercaseAccessionPrefix.is_correctable());
        assert!(ErrorType::MixedCaseEditType.is_correctable());
        assert!(ErrorType::OldSubstitutionSyntax.is_correctable());
        assert!(ErrorType::InvalidUnicodeCharacter.is_correctable());
        assert!(ErrorType::SwappedPositions.is_correctable());
        assert!(ErrorType::TrailingAnnotation.is_correctable());
        assert!(ErrorType::MissingCoordinatePrefix.is_correctable());
        assert!(ErrorType::OldAlleleFormat.is_correctable());
        assert!(ErrorType::RefSeqMismatch.is_correctable());
        assert!(!ErrorType::DeprecatedIvsNotation.is_correctable());
        assert!(ErrorType::DeprecatedConSyntax.is_correctable());
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
        assert_eq!(ErrorType::DeprecatedStopCodonStar.code(), "W3007");
        assert_eq!(ErrorType::DeprecatedStopCodonX.code(), "W3008");
        assert_eq!(ErrorType::DeprecatedFrameshiftStar.code(), "W3009");
        assert_eq!(ErrorType::DeprecatedFrameshiftX.code(), "W3010");
        assert_eq!(ErrorType::DelSizeSuffix.code(), "W3011");
        assert_eq!(ErrorType::EmptyDelinsInsert.code(), "W3012");
        assert_eq!(ErrorType::RedundantRepeatLabel.code(), "W3013");
        assert_eq!(ErrorType::DeprecatedIvsNotation.code(), "W3014");
        assert_eq!(ErrorType::DeprecatedConSyntax.code(), "W3015");
        assert_eq!(ErrorType::SwappedPositions.code(), "W4001");
        assert_eq!(ErrorType::PositionZero.code(), "W4002");
        assert_eq!(ErrorType::SinglePositionRange.code(), "W4003");
        assert_eq!(ErrorType::RefSeqMismatch.code(), "W5001");
    }

    #[test]
    fn test_error_type_is_correctable_new_variants() {
        // IVS notation cannot be auto-corrected without metadata.
        assert!(!ErrorType::DeprecatedIvsNotation.is_correctable());
        // Con syntax CAN be rewritten to delins.
        assert!(ErrorType::DeprecatedConSyntax.is_correctable());
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
