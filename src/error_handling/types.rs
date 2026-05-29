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

    // Note: the W4002 `PositionZero` variant was retired in issue #269.
    // `c.0…` inputs are rejected by `detect_position_zero` as
    // `ErrorCode::InvalidPosition` (E1003) — a hard parse error, not a
    // soft-validation warning — which is the canonical (and only) code
    // for the symptom now.
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

    /// Length mismatch between an explicit reference sequence and the
    /// position range it spans.
    ///
    /// For `del` / `dup` / `inv` / `delins` with an explicit reference
    /// sequence, the sequence length must equal the range length
    /// (`end - start + 1`). `g.100_110delAAAATTTGCC` has range 11 but
    /// sequence length 10. Lenient mode emits W3016 without rewriting
    /// (the principled correction depends on which endpoint is wrong,
    /// which the parser cannot decide); strict mode rejects.
    LengthMismatch,

    /// Allele-fraction / heteroplasmy annotation appended to an HGVS
    /// expression (e.g. `m.3243A>G[level=70%]`, `m.3243A>G(80%)`).
    ///
    /// The HGVS spec does not encode allele fraction inside the variant
    /// string; fraction belongs in accompanying metadata (VCF FORMAT/AF,
    /// ClinVar's heteroplasmy field, etc.). All modes reject this
    /// annotation with a targeted W3017 diagnostic so downstream tooling
    /// can recognize the intent rather than seeing a generic
    /// trailing-character parse error.
    AlleleFractionAnnotation,

    /// ClinVar prose multi-allelic shorthand
    /// `<acc>:m.<pos><ref>><alt>/<alt2>` (e.g. `m.3243A>G/T`) where the
    /// RHS of the slash is a bare nucleotide letter — no edit operator,
    /// no accession, no `=` reference marker.
    ///
    /// This shape is not a HGVS spec form: the spec's mosaic `/`
    /// separator requires either two fully-qualified variants on each
    /// side, or the spec compact mosaic `<pos>=/<alt-edit>` form (`=`
    /// stands in for the reference, edit names the alternative). All
    /// modes reject this annotation with a targeted W3018 diagnostic
    /// pointing at the three spec-supported alternatives.
    ClinVarProseMultiAllelic,

    /// Thymine (`t`/`T`) used as a base in an `r.` (RNA) description.
    ///
    /// Per HGVS v21.0 RNA nomenclature the RNA alphabet is `a/c/g/u`;
    /// `t` is non-canonical input. ferro canonicalizes `t` → `u` inside
    /// any `r.` description in lenient/silent modes; strict mode rejects.
    /// Issue #282 (closes #232 follow-up).
    RnaThymineCanonicalized,

    /// Non-spec mosaic / chimeric form: nested `/` + `//` at the same
    /// nesting level (e.g. `m.[3243A>G/T]//[3243A>C/G]`,
    /// chimeric-of-mosaic or mosaic-of-chimeric) **and** the
    /// bracketed `[a/b]` / `[a//b]` mosaic group (a form HGVS does
    /// not define). One code, two detectors — both surface the same
    /// W3019 family so downstream tools can key off a single
    /// diagnostic. Always-reject: no canonical alternative exists, so
    /// auto-correction would be a guess.
    NonSpecMosaicForm,

    /// Bracketed amino-acid list inside a protein insertion edit
    /// (e.g., `p.Arg97_Trp98ins[Ala;Pro]`).
    ///
    /// HGVS v21's protein insertion notation concatenates 3-letter codes
    /// without separators (`p.Arg97_Trp98insAlaPro`); brackets are reserved
    /// for alleles at the variant level, not for amino-acid lists inside an
    /// edit. Cannot be auto-rewritten without ambiguity (single-letter vs
    /// three-letter mixing), so all modes reject with an actionable hint
    /// pointing at the canonical `insAlaPro` form.
    ProteinBracketedAaInsertion,

    /// The canonical `p.dup` form produced by ins→dup or delins→dup
    /// canonicalization includes the initiator methionine (position 1).
    ///
    /// HGVS Prioritization (`general.md §3`) requires `dup` to be preferred
    /// over `ins` when both describe the same change; the rule is
    /// unconditional, so ferro rewrites to `dup` even when the interval
    /// includes p.1. The spec uses Met1-inclusive ranges elsewhere
    /// (`deletion.md §63-65`: `p.(Met1_Leu46del)`), so the form is
    /// permitted. However, the predicted protein-level consequence may
    /// also be described per the substitution rule for start-codon variants
    /// (`substitution.md §45-65`: `p.0`, `p.0?`, or `p.(Met1?)`). This
    /// warning fires so downstream consumers can apply that policy.
    /// Closes-after: #92. Code: `W3022`.
    InitiatorMetCanonicalization,

    /// Duplication described with a size-count suffix instead of a position
    /// range (e.g. `g.123dup6` instead of `g.123_128dup`).
    ///
    /// Per HGVS spec (`recommendations/DNA/duplication.md:140-143` Q&A):
    /// "a duplication of more than one nucleotide should give the position
    /// of the first and last nucleotide duplicated, separated using the
    /// range symbol (`_`, underscore), e.g., `g.123_128dup`." Position
    /// ambiguity ("at" vs "after") prevents safe rewrite.
    DupSizeSuffix,

    /// Duplication described with an explicit duplicated sequence (e.g.
    /// `c.20_23dupTAGA` instead of `c.20_23dup`).
    ///
    /// Per HGVS spec (`recommendations/DNA/duplication.md:35-36, 50-51`):
    /// "the recommendation is not to describe the variant as
    /// `c.20_23dupTAGA`. This description is longer, it contains redundant
    /// information, and chances to make an error increases." Sequence is
    /// safely droppable since the position range fully determines it.
    DupExplicitSeq,

    /// Deletion described with an explicit deleted sequence (e.g.
    /// `g.33344590_33344592delGAT` instead of `g.33344590_33344592del`).
    ///
    /// Per HGVS spec (`recommendations/DNA/deletion.md:30-31, 45-46`):
    /// "the recommendation is not to describe the variant as
    /// `NC_000023.11:g.33344590_33344592delGAT`. This description is
    /// longer, it contains redundant information, and chances to make an
    /// error increases." Sequence is safely droppable.
    DelExplicitSeq,

    /// Variant position range exceeds the reference sequence the
    /// provider returned (`fetch_ref_for_canonical_split` got fewer
    /// bytes than the HGVS interval span).
    ///
    /// Per HGVS spec `recommendations/background/refseq.md` §43: "the
    /// entirety of the variant sequence **must** be encompassed by the
    /// selected reference sequence." Biocommons hgvs raises
    /// `HGVSInvalidVariantError` for this shape. Strict mode rejects
    /// with `FerroError::VariantExceedsReference`; lenient mode emits
    /// the existing `CanonicalSplitSkipped` warning and preserves the
    /// input; silent mode preserves the input without warning.
    /// Closes-after: #355.
    VariantExceedsReference,

    /// A `c.`, `c.*N`, `c.-N`, or `n.` position lies past the bound of
    /// its coordinate sub-axis:
    ///   - `c.<N>` past `cds_end - cds_start + 1`
    ///   - `c.*<N>` past `sequence_length - cds_end`
    ///   - `c.-<N>` past `cds_start - 1`
    ///   - `n.<N>` past `sequence_length`
    ///
    /// Examples: `c.946G>C` on a transcript whose CDS is only 945 bases;
    /// `c.*9G>C` on a transcript whose 3'UTR is only 8 bases;
    /// `c.-6G>C` on a transcript whose 5'UTR is only 5 bases;
    /// `n.21G>C` on a 20-base transcript. Strict mode rejects; lenient
    /// emits W4004 and short-circuits normalize() to the canonical
    /// variant; silent mode skips the check entirely. Intronic offsets
    /// remain out of scope (their bounds depend on intron size).
    PositionPastEnd,

    /// Two or more cis-allele edits share identical reference bounds
    /// (e.g. `g.[100A>C;100A>G]` — both substitutions at position 100).
    /// The HGVS spec does not define a canonical form for this case;
    /// ferro preserves the input verbatim and emits this warning.
    /// Strict mode rejects with `FerroError::InvalidCoordinates`;
    /// lenient and silent modes both accept and preserve the input,
    /// with the warning still surfaced by `normalize_with_diagnostics`
    /// (the emit site at `src/normalize/overlap.rs:88` is
    /// unconditional — only strict-mode promotion is gated by mode).
    /// Closes #395 item 6 — previously the warning was emitted but
    /// strict mode did not reject, bypassing the registry's
    /// `always_warn_if_not_rejected` policy table that declared
    /// Strict→Reject.
    OverlapConflictingEdits,
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
            ErrorType::LengthMismatch => "W3016",
            ErrorType::AlleleFractionAnnotation => "W3017",
            ErrorType::ClinVarProseMultiAllelic => "W3018",
            ErrorType::NonSpecMosaicForm => "W3019",
            ErrorType::RnaThymineCanonicalized => "W3020",
            ErrorType::ProteinBracketedAaInsertion => "W3021",
            ErrorType::InitiatorMetCanonicalization => "W3022",
            ErrorType::DupSizeSuffix => "W3023",
            ErrorType::DupExplicitSeq => "W3024",
            ErrorType::DelExplicitSeq => "W3025",
            ErrorType::SwappedPositions => "W4001",
            ErrorType::SinglePositionRange => "W4003",
            ErrorType::PositionPastEnd => "W4004",
            ErrorType::RefSeqMismatch => "W5001",
            ErrorType::OverlapConflictingEdits => "W5002",
            ErrorType::VariantExceedsReference => "W5003",
        }
    }

    /// Every `ErrorType` variant, in W-code order.
    ///
    /// This is the authoritative list used by `from_code` and by config
    /// parsing to map a user-facing W-code back to its `ErrorType`. Adding
    /// a new variant requires extending both `code()` (compile-time
    /// exhaustiveness) and this array (the
    /// `test_from_code_round_trips_all_variants` test asserts that every
    /// code returned by `code()` resolves via `from_code`).
    pub const ALL: &'static [ErrorType] = &[
        ErrorType::LowercaseAminoAcid,
        ErrorType::SingleLetterAminoAcid,
        ErrorType::LowercaseAccessionPrefix,
        ErrorType::MixedCaseEditType,
        ErrorType::WrongDashCharacter,
        ErrorType::WrongQuoteCharacter,
        ErrorType::ExtraWhitespace,
        ErrorType::InvalidUnicodeCharacter,
        ErrorType::MissingVersion,
        ErrorType::ProteinSubstitutionArrow,
        ErrorType::OldSubstitutionSyntax,
        ErrorType::OldAlleleFormat,
        ErrorType::TrailingAnnotation,
        ErrorType::MissingCoordinatePrefix,
        ErrorType::DeprecatedStopCodonStar,
        ErrorType::DeprecatedStopCodonX,
        ErrorType::DeprecatedFrameshiftStar,
        ErrorType::DeprecatedFrameshiftX,
        ErrorType::DelSizeSuffix,
        ErrorType::EmptyDelinsInsert,
        ErrorType::RedundantRepeatLabel,
        ErrorType::DeprecatedIvsNotation,
        ErrorType::DeprecatedConSyntax,
        ErrorType::LengthMismatch,
        ErrorType::AlleleFractionAnnotation,
        ErrorType::ClinVarProseMultiAllelic,
        ErrorType::NonSpecMosaicForm,
        ErrorType::RnaThymineCanonicalized,
        ErrorType::ProteinBracketedAaInsertion,
        ErrorType::InitiatorMetCanonicalization,
        ErrorType::DupSizeSuffix,
        ErrorType::DupExplicitSeq,
        ErrorType::DelExplicitSeq,
        ErrorType::SwappedPositions,
        ErrorType::SinglePositionRange,
        ErrorType::PositionPastEnd,
        ErrorType::RefSeqMismatch,
        ErrorType::OverlapConflictingEdits,
        ErrorType::VariantExceedsReference,
    ];

    /// Parse a warning code (e.g. `"W3007"`) into its `ErrorType`.
    ///
    /// Returns `None` when the code is not (or is no longer) a registered
    /// soft-validation warning. The retired W4002 `PositionZero` code
    /// (issue #269) is one such case — `c.0…` inputs are now rejected at
    /// parse time as `E1003 InvalidPosition`. Comparison is
    /// case-insensitive so `.ferro.toml` entries like `"w3003"` resolve.
    ///
    /// The mapping is derived from `ErrorType::ALL` and
    /// `ErrorType::code()` so it cannot drift away from the variant set;
    /// see `test_from_code_round_trips_all_variants`.
    pub fn from_code(code: &str) -> Option<ErrorType> {
        let upper = code.to_uppercase();
        ErrorType::ALL.iter().copied().find(|t| t.code() == upper)
    }

    /// Returns a human-readable description of this error type.
    pub fn description(&self) -> &'static str {
        match self {
            ErrorType::LowercaseAminoAcid => "lowercase amino acid code",
            ErrorType::MissingVersion => "missing accession version",
            ErrorType::WrongDashCharacter => "wrong dash character (en-dash or em-dash)",
            ErrorType::ExtraWhitespace => "extra whitespace in description",
            ErrorType::ProteinSubstitutionArrow => "arrow in protein substitution",
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
            ErrorType::LengthMismatch => {
                "explicit reference sequence length does not match position range"
            }
            ErrorType::AlleleFractionAnnotation => "allele-fraction / heteroplasmy annotation",
            ErrorType::ClinVarProseMultiAllelic => {
                "ClinVar prose multi-allelic shorthand m.<pos><ref>><alt>/<alt2>"
            }
            ErrorType::RnaThymineCanonicalized => {
                "thymine (t) used in r. RNA description; canonicalized to u"
            }
            ErrorType::NonSpecMosaicForm => {
                "non-spec mosaic/chimeric form (nested / + // or bracketed [a/b])"
            }
            ErrorType::ProteinBracketedAaInsertion => {
                "bracketed amino-acid list inside protein insertion edit"
            }
            ErrorType::InitiatorMetCanonicalization => {
                "canonical form `p.Met1dup` (or analogous) includes the initiator methionine"
            }
            ErrorType::DupSizeSuffix => "duplication described with a size-count suffix",
            ErrorType::DupExplicitSeq => "duplication with explicit duplicated sequence",
            ErrorType::DelExplicitSeq => "deletion with explicit deleted sequence",
            ErrorType::VariantExceedsReference => {
                "variant position range exceeds reference sequence"
            }
            ErrorType::PositionPastEnd => "position lies past CDS-end or transcript-end",
            ErrorType::OverlapConflictingEdits => {
                "two or more cis-allele edits share identical reference bounds"
            }
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
            // Length-mismatch input has no safe auto-correction — the user
            // could have meant either endpoint or a different ref seq.
            ErrorType::LengthMismatch => false,
            // Allele-fraction annotations carry data that does not belong
            // in an HGVS string; no auto-correction is meaningful.
            ErrorType::AlleleFractionAnnotation => false,
            // The malformed ClinVar prose shape carries no shared
            // structure with the spec-supported alternatives; ferro
            // cannot pick one for the caller, so this is never
            // auto-corrected.
            ErrorType::ClinVarProseMultiAllelic => false,
            // RNA thymine is rewritten to `u`
            ErrorType::RnaThymineCanonicalized => true,
            // No canonical alternative for nested mosaic/chimeric or
            // bracketed [a/b] — see W3019 explanation in registry.
            ErrorType::NonSpecMosaicForm => false,
            // Bracketed AA list inside an insertion has no spec-defined
            // rewrite: mixing 3-letter and 1-letter inside `[...]` is
            // ambiguous, so all modes reject with a hint.
            ErrorType::ProteinBracketedAaInsertion => false,
            // The dup form is correct per HGVS Prioritization; the warning
            // is purely informational. We do not rewrite the canonical output
            // back to ins or to p.0?/p.(Met1?) — that is the consumer's call.
            ErrorType::InitiatorMetCanonicalization => false,
            // Lenient mode warns without rewriting (warn_accept); strict rejects.
            // Position ambiguity ("at" vs "after") prevents safe synthesis.
            ErrorType::DupSizeSuffix => false,
            // Sequence is redundant given the position range; lenient mode
            // drops it (standard_correctable).
            ErrorType::DupExplicitSeq => true,
            ErrorType::DelExplicitSeq => true,
            // Variant exceeds reference: no safe auto-correction — the
            // variant references positions past the provider's window
            // and we cannot conjure missing bases.
            ErrorType::VariantExceedsReference => false,
            // Past-end positions cannot be auto-corrected — there's no safe
            // way to guess the intended canonical position.
            ErrorType::PositionPastEnd => false,
            // Coincident-bounds cis-allele edits have no spec-defined
            // canonical form; ferro preserves the input verbatim.
            ErrorType::OverlapConflictingEdits => false,
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
            ErrorType::LengthMismatch => ("g.100_110delAAAATTTGCC", "(no auto-correct)"),
            ErrorType::AlleleFractionAnnotation => {
                ("NC_012920.1:m.3243A>G[level=70%]", "NC_012920.1:m.3243A>G")
            }
            ErrorType::ClinVarProseMultiAllelic => {
                ("NC_012920.1:m.3243A>G/T", "NC_012920.1:m.[3243A>G;3243A>T]")
            }
            ErrorType::RnaThymineCanonicalized => ("r.123a>t", "r.123a>u"),
            ErrorType::NonSpecMosaicForm => (
                "NM_000088.3:c.[100A>G/200T>C]",
                "NM_000088.3:c.[100A>G;200T>C] or NM_000088.3:c.100A>G/NM_000088.3:c.200T>C",
            ),
            ErrorType::ProteinBracketedAaInsertion => {
                ("p.Arg97_Trp98ins[Ala;Pro]", "p.Arg97_Trp98insAlaPro")
            }
            ErrorType::InitiatorMetCanonicalization => ("p.Met1_Lys2insMet", "p.Met1dup"),
            ErrorType::DupSizeSuffix => ("g.123dup6", "g.123_128dup"),
            ErrorType::DupExplicitSeq => ("c.20_23dupTAGA", "c.20_23dup"),
            ErrorType::DelExplicitSeq => ("g.33344590_33344592delGAT", "g.33344590_33344592del"),
            ErrorType::VariantExceedsReference => (
                "NG_032871.1:g.32476_53457delinsAATTAAGGTATA (ref shorter than 53457)",
                "(no auto-correct; spec rejects per refseq.md \u{00A7}43)",
            ),
            ErrorType::PositionPastEnd => ("c.946G>C (CDS length 945)", "(no auto-correct)"),
            ErrorType::OverlapConflictingEdits => (
                "g.[100A>C;100A>G] (two cis edits at the same bound)",
                "(no auto-correct)",
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
        assert!(!ErrorType::LengthMismatch.is_correctable());
        assert!(ErrorType::RnaThymineCanonicalized.is_correctable());
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
        assert_eq!(ErrorType::LengthMismatch.code(), "W3016");
        assert_eq!(ErrorType::NonSpecMosaicForm.code(), "W3019");
        assert_eq!(ErrorType::RnaThymineCanonicalized.code(), "W3020");
        assert_eq!(ErrorType::ProteinBracketedAaInsertion.code(), "W3021");
        assert_eq!(ErrorType::InitiatorMetCanonicalization.code(), "W3022");
        assert_eq!(ErrorType::SwappedPositions.code(), "W4001");
        assert_eq!(ErrorType::SinglePositionRange.code(), "W4003");
        assert_eq!(ErrorType::PositionPastEnd.code(), "W4004");
        assert_eq!(ErrorType::RefSeqMismatch.code(), "W5001");
        assert_eq!(ErrorType::VariantExceedsReference.code(), "W5003");
    }

    #[test]
    fn test_error_type_is_correctable_new_variants() {
        // IVS notation cannot be auto-corrected without metadata.
        assert!(!ErrorType::DeprecatedIvsNotation.is_correctable());
        // Con syntax CAN be rewritten to delins.
        assert!(ErrorType::DeprecatedConSyntax.is_correctable());
        // Non-spec mosaic / chimeric forms have no canonical
        // alternative and cannot be auto-corrected.
        assert!(!ErrorType::NonSpecMosaicForm.is_correctable());
        // Past-end positions could mean a different position or a different
        // reference — no safe auto-correction.
        assert!(!ErrorType::PositionPastEnd.is_correctable());
        // Variants exceeding the reference bounds have no safe correction.
        assert!(!ErrorType::VariantExceedsReference.is_correctable());
        // InitiatorMetCanonicalization is informational only — the dup
        // form is kept; no rewrite back to ins or to p.0?/p.(Met1?).
        assert!(!ErrorType::InitiatorMetCanonicalization.is_correctable());
        // DupSizeSuffix (W3023) uses warn_accept — the input is preserved as-is,
        // so there is no rewrite to apply.
        assert!(!ErrorType::DupSizeSuffix.is_correctable());
        // DupExplicitSeq (W3024) and DelExplicitSeq (W3025) use
        // standard_correctable — the redundant sequence suffix is stripped.
        assert!(ErrorType::DupExplicitSeq.is_correctable());
        assert!(ErrorType::DelExplicitSeq.is_correctable());
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

    // ErrorType::ALL / from_code tests
    #[test]
    fn test_all_codes_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for variant in ErrorType::ALL {
            assert!(
                seen.insert(variant.code()),
                "duplicate code {} in ErrorType::ALL",
                variant.code()
            );
        }
    }

    #[test]
    fn test_from_code_round_trips_all_variants() {
        // Every variant in `ALL` must be reachable via `from_code(code())`.
        for variant in ErrorType::ALL {
            let code = variant.code();
            assert_eq!(
                ErrorType::from_code(code),
                Some(*variant),
                "ErrorType::from_code({code:?}) did not round-trip"
            );
        }
    }

    #[test]
    fn test_all_contains_every_error_type_variant() {
        // Compile-time + runtime exhaustiveness: every `ErrorType` variant
        // must appear in `ErrorType::ALL` so `from_code` can resolve it.
        //
        // The `match` arms below enumerate every variant; adding a new
        // variant to `ErrorType` produces a `non-exhaustive patterns`
        // compile error here, forcing the author to also extend `ALL`.
        // The runtime `contains(&v)` check catches the case where the
        // arm is added but the entry in `ALL` is forgotten.
        let every_variant: [ErrorType; ErrorType::ALL.len()] = [
            ErrorType::LowercaseAminoAcid,
            ErrorType::SingleLetterAminoAcid,
            ErrorType::LowercaseAccessionPrefix,
            ErrorType::MixedCaseEditType,
            ErrorType::WrongDashCharacter,
            ErrorType::WrongQuoteCharacter,
            ErrorType::ExtraWhitespace,
            ErrorType::InvalidUnicodeCharacter,
            ErrorType::MissingVersion,
            ErrorType::ProteinSubstitutionArrow,
            ErrorType::OldSubstitutionSyntax,
            ErrorType::OldAlleleFormat,
            ErrorType::TrailingAnnotation,
            ErrorType::MissingCoordinatePrefix,
            ErrorType::DeprecatedStopCodonStar,
            ErrorType::DeprecatedStopCodonX,
            ErrorType::DeprecatedFrameshiftStar,
            ErrorType::DeprecatedFrameshiftX,
            ErrorType::DelSizeSuffix,
            ErrorType::EmptyDelinsInsert,
            ErrorType::RedundantRepeatLabel,
            ErrorType::DeprecatedIvsNotation,
            ErrorType::DeprecatedConSyntax,
            ErrorType::LengthMismatch,
            ErrorType::AlleleFractionAnnotation,
            ErrorType::ClinVarProseMultiAllelic,
            ErrorType::NonSpecMosaicForm,
            ErrorType::RnaThymineCanonicalized,
            ErrorType::ProteinBracketedAaInsertion,
            ErrorType::InitiatorMetCanonicalization,
            ErrorType::DupSizeSuffix,
            ErrorType::DupExplicitSeq,
            ErrorType::DelExplicitSeq,
            ErrorType::SwappedPositions,
            ErrorType::SinglePositionRange,
            ErrorType::PositionPastEnd,
            ErrorType::RefSeqMismatch,
            ErrorType::OverlapConflictingEdits,
            ErrorType::VariantExceedsReference,
        ];
        // Exhaustiveness probe: if a new variant is added to the enum,
        // this match fails to compile until `every_variant` is extended.
        // The match must be on *every* possible value of `ErrorType`, so
        // we run it once per slot above.
        for v in every_variant {
            match v {
                ErrorType::LowercaseAminoAcid
                | ErrorType::SingleLetterAminoAcid
                | ErrorType::LowercaseAccessionPrefix
                | ErrorType::MixedCaseEditType
                | ErrorType::WrongDashCharacter
                | ErrorType::WrongQuoteCharacter
                | ErrorType::ExtraWhitespace
                | ErrorType::InvalidUnicodeCharacter
                | ErrorType::MissingVersion
                | ErrorType::ProteinSubstitutionArrow
                | ErrorType::OldSubstitutionSyntax
                | ErrorType::OldAlleleFormat
                | ErrorType::TrailingAnnotation
                | ErrorType::MissingCoordinatePrefix
                | ErrorType::DeprecatedStopCodonStar
                | ErrorType::DeprecatedStopCodonX
                | ErrorType::DeprecatedFrameshiftStar
                | ErrorType::DeprecatedFrameshiftX
                | ErrorType::DelSizeSuffix
                | ErrorType::EmptyDelinsInsert
                | ErrorType::RedundantRepeatLabel
                | ErrorType::DeprecatedIvsNotation
                | ErrorType::DeprecatedConSyntax
                | ErrorType::LengthMismatch
                | ErrorType::AlleleFractionAnnotation
                | ErrorType::ClinVarProseMultiAllelic
                | ErrorType::NonSpecMosaicForm
                | ErrorType::RnaThymineCanonicalized
                | ErrorType::ProteinBracketedAaInsertion
                | ErrorType::InitiatorMetCanonicalization
                | ErrorType::DupSizeSuffix
                | ErrorType::DupExplicitSeq
                | ErrorType::DelExplicitSeq
                | ErrorType::SwappedPositions
                | ErrorType::SinglePositionRange
                | ErrorType::PositionPastEnd
                | ErrorType::RefSeqMismatch
                | ErrorType::OverlapConflictingEdits
                | ErrorType::VariantExceedsReference => {}
            }
            assert!(
                ErrorType::ALL.contains(&v),
                "ErrorType::{v:?} is missing from ErrorType::ALL"
            );
        }
        // Belt-and-braces: lengths must agree.
        assert_eq!(every_variant.len(), ErrorType::ALL.len());
    }

    #[test]
    fn test_from_code_is_case_insensitive() {
        assert_eq!(
            ErrorType::from_code("w3003"),
            Some(ErrorType::OldSubstitutionSyntax)
        );
        assert_eq!(
            ErrorType::from_code("W3003"),
            Some(ErrorType::OldSubstitutionSyntax)
        );
    }

    #[test]
    fn test_from_code_covers_codes_previously_missing_from_config() {
        // Regression for CodeRabbit PR #370 review: these registered
        // warning codes were missing from the old hard-coded
        // `code_to_error_type` table in `src/config.rs`, so
        // `.ferro.toml` overrides targeting them silently no-op'd.
        let cases = [
            ("W3007", ErrorType::DeprecatedStopCodonStar),
            ("W3008", ErrorType::DeprecatedStopCodonX),
            ("W3009", ErrorType::DeprecatedFrameshiftStar),
            ("W3010", ErrorType::DeprecatedFrameshiftX),
            ("W3011", ErrorType::DelSizeSuffix),
            ("W3012", ErrorType::EmptyDelinsInsert),
            ("W3013", ErrorType::RedundantRepeatLabel),
            ("W3014", ErrorType::DeprecatedIvsNotation),
            ("W3015", ErrorType::DeprecatedConSyntax),
            ("W3017", ErrorType::AlleleFractionAnnotation),
            ("W3018", ErrorType::ClinVarProseMultiAllelic),
            ("W3020", ErrorType::RnaThymineCanonicalized),
            ("W3021", ErrorType::ProteinBracketedAaInsertion),
            ("W3022", ErrorType::InitiatorMetCanonicalization),
            ("W4003", ErrorType::SinglePositionRange),
        ];
        for (code, expected) in cases {
            assert_eq!(
                ErrorType::from_code(code),
                Some(expected),
                "ErrorType::from_code({code:?}) regressed"
            );
        }
    }

    #[test]
    fn test_from_code_rejects_retired_and_unknown() {
        // W4002 was retired in issue #269.
        assert_eq!(ErrorType::from_code("W4002"), None);
        // Bare nonsense.
        assert_eq!(ErrorType::from_code("W9999"), None);
        assert_eq!(ErrorType::from_code(""), None);
        assert_eq!(ErrorType::from_code("not-a-code"), None);
    }
}
