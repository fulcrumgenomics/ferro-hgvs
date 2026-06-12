//! Mutalyzer ↔ ferro error-code mapping table.
//!
//! Provides translation between mutalyzer/mutalyzer's E-code taxonomy and
//! ferro-hgvs's `FerroError` variants for cross-tool diagnostics, parity
//! audits, and the mutalyzer corpus runner (`tests/mutalyzer_normalize_tests.rs`).
//!
//! The forward map (`mutalyzer_to_ferro`) is keyed by the mutalyzer error code
//! string (e.g. `"ESEQUENCEMISMATCH"`) and yields a [`FerroErrorTag`] whose
//! `debug_tag()` returns the `Debug` representation of the corresponding
//! `FerroError` variant (e.g. `"ReferenceMismatch"`). The reverse map
//! (`ferro_to_mutalyzer`) yields every mutalyzer code that *may* correspond to
//! a given `FerroError` — there is no one-to-one bijection because mutalyzer's
//! taxonomy is finer-grained.
//!
//! # Coverage
//!
//! Every E-code emitted by mutalyzer/mutalyzer's `errors.py` at the pinned
//! upstream commit (`00045d6`, the version PR #323's corpus targets) is
//! classified — either to a concrete ferro variant or explicitly to "no ferro
//! equivalent" via the [`NO_FERRO_EQUIV`] allowlist. The exhaustiveness test
//! [`tests::forward_all_upstream_codes_classified`] guarantees no code drifts
//! into an unclassified state.
//!
//! # Lossy mappings
//!
//! Several mutalyzer codes lack an exact ferro counterpart and are
//! deliberately mapped to the closest structural fit with `is_lossy() == true`
//! to record the approximation. Examples:
//!
//! - `EOVERLAP` → `InvalidCoordinates` (ferro's `ErrorCode::SelfCancellingAllele`
//!   is structurally close but lives in a different taxonomy)
//! - `ESEQUENCEMISMATCH`, `ELENGTHMISMATCH` → `ReferenceMismatch`
//! - `ECOORDINATESYSTEMMISMATCH` → `Parse` (remapped in #543; was `UnsupportedVariant`)
//! - `EOUTOFBOUNDARY` → `InvalidCoordinates`

use crate::error::FerroError;

/// Strongly-typed mutalyzer error code.
///
/// Closed enum covering every E-code emitted by mutalyzer/mutalyzer at the
/// pinned upstream commit (`00045d6`). New codes added upstream will fail to
/// parse via [`MutalyzerCode::parse`] until added here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MutalyzerCode {
    EAminoAcidMismatch,
    ECdsSlices,
    ECoordinateSystemInvalid,
    ECoordinateSystemMismatch,
    EInsertedLength,
    EInsertionLocation,
    EInsertionRange,
    EIntronic,
    EIntronicRna,
    EInvalidInput,
    ELengthMismatch,
    ELengthsDifference,
    ELocationSlice,
    EMismatch,
    EMissingParameter,
    ENoCds,
    ENoCoordinateSystem,
    ENoDna,
    ENoInputs,
    ENoInputsOther,
    ENoRna,
    ENoSelectorFound,
    ENoToSelector,
    EOffset,
    EOutOfBoundary,
    EOutsideCds,
    EOverlap,
    EPositionInvalid,
    EPositionSyntax,
    ERangeReversed,
    ERepeatMismatch,
    ERepeatReferenceLength,
    ERepeatUnsupported,
    ERetr,
    ESelectorOptions,
    ESequenceLength,
    ESequenceMismatch,
    ESliceOption,
    ESpliceSite,
    ESyntaxNested,
    ESyntaxUc,
    ESyntaxUeof,
    EUncertain,
    EVariantNotSupported,
}

impl MutalyzerCode {
    /// Parse a mutalyzer code string (e.g. `"ESEQUENCEMISMATCH"`).
    ///
    /// Returns `None` for unknown codes.
    pub fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "EAMINOACIDMISMATCH" => Self::EAminoAcidMismatch,
            "ECDSSLICES" => Self::ECdsSlices,
            "ECOORDINATESYSTEMINVALID" => Self::ECoordinateSystemInvalid,
            "ECOORDINATESYSTEMMISMATCH" => Self::ECoordinateSystemMismatch,
            "EINSERTEDLENGTH" => Self::EInsertedLength,
            "EINSERTIONLOCATION" => Self::EInsertionLocation,
            "EINSERTIONRANGE" => Self::EInsertionRange,
            "EINTRONIC" => Self::EIntronic,
            "EINTRONICRNA" => Self::EIntronicRna,
            "EINVALIDINPUT" => Self::EInvalidInput,
            "ELENGTHMISMATCH" => Self::ELengthMismatch,
            "ELENGTHSDIFFERENCE" => Self::ELengthsDifference,
            "ELOCATIONSLICE" => Self::ELocationSlice,
            "EMISMATCH" => Self::EMismatch,
            "EMISSINGPARAMETER" => Self::EMissingParameter,
            "ENOCDS" => Self::ENoCds,
            "ENOCOORDINATESYSTEM" => Self::ENoCoordinateSystem,
            "ENODNA" => Self::ENoDna,
            "ENOINPUTS" => Self::ENoInputs,
            "ENOINPUTSOTHER" => Self::ENoInputsOther,
            "ENORNA" => Self::ENoRna,
            "ENOSELECTORFOUND" => Self::ENoSelectorFound,
            "ENOTOSELECTOR" => Self::ENoToSelector,
            "EOFFSET" => Self::EOffset,
            "EOUTOFBOUNDARY" => Self::EOutOfBoundary,
            "EOUTSIDECDS" => Self::EOutsideCds,
            "EOVERLAP" => Self::EOverlap,
            "EPOSITIONINVALID" => Self::EPositionInvalid,
            "EPOSITIONSYNTAX" => Self::EPositionSyntax,
            "ERANGEREVERSED" => Self::ERangeReversed,
            "EREPEATMISMATCH" => Self::ERepeatMismatch,
            "EREPEATREFERENCELENGTH" => Self::ERepeatReferenceLength,
            "EREPEATUNSUPPORTED" => Self::ERepeatUnsupported,
            "ERETR" => Self::ERetr,
            "ESELECTOROPTIONS" => Self::ESelectorOptions,
            "ESEQUENCELENGTH" => Self::ESequenceLength,
            "ESEQUENCEMISMATCH" => Self::ESequenceMismatch,
            "ESLICEOPTION" => Self::ESliceOption,
            "ESPLICESITE" => Self::ESpliceSite,
            "ESYNTAXNESTED" => Self::ESyntaxNested,
            "ESYNTAXUC" => Self::ESyntaxUc,
            "ESYNTAXUEOF" => Self::ESyntaxUeof,
            "EUNCERTAIN" => Self::EUncertain,
            "EVARIANTNOTSUPPORTED" => Self::EVariantNotSupported,
            _ => return None,
        })
    }

    /// Get the canonical upper-case mutalyzer string for this code.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::EAminoAcidMismatch => "EAMINOACIDMISMATCH",
            Self::ECdsSlices => "ECDSSLICES",
            Self::ECoordinateSystemInvalid => "ECOORDINATESYSTEMINVALID",
            Self::ECoordinateSystemMismatch => "ECOORDINATESYSTEMMISMATCH",
            Self::EInsertedLength => "EINSERTEDLENGTH",
            Self::EInsertionLocation => "EINSERTIONLOCATION",
            Self::EInsertionRange => "EINSERTIONRANGE",
            Self::EIntronic => "EINTRONIC",
            Self::EIntronicRna => "EINTRONICRNA",
            Self::EInvalidInput => "EINVALIDINPUT",
            Self::ELengthMismatch => "ELENGTHMISMATCH",
            Self::ELengthsDifference => "ELENGTHSDIFFERENCE",
            Self::ELocationSlice => "ELOCATIONSLICE",
            Self::EMismatch => "EMISMATCH",
            Self::EMissingParameter => "EMISSINGPARAMETER",
            Self::ENoCds => "ENOCDS",
            Self::ENoCoordinateSystem => "ENOCOORDINATESYSTEM",
            Self::ENoDna => "ENODNA",
            Self::ENoInputs => "ENOINPUTS",
            Self::ENoInputsOther => "ENOINPUTSOTHER",
            Self::ENoRna => "ENORNA",
            Self::ENoSelectorFound => "ENOSELECTORFOUND",
            Self::ENoToSelector => "ENOTOSELECTOR",
            Self::EOffset => "EOFFSET",
            Self::EOutOfBoundary => "EOUTOFBOUNDARY",
            Self::EOutsideCds => "EOUTSIDECDS",
            Self::EOverlap => "EOVERLAP",
            Self::EPositionInvalid => "EPOSITIONINVALID",
            Self::EPositionSyntax => "EPOSITIONSYNTAX",
            Self::ERangeReversed => "ERANGEREVERSED",
            Self::ERepeatMismatch => "EREPEATMISMATCH",
            Self::ERepeatReferenceLength => "EREPEATREFERENCELENGTH",
            Self::ERepeatUnsupported => "EREPEATUNSUPPORTED",
            Self::ERetr => "ERETR",
            Self::ESelectorOptions => "ESELECTOROPTIONS",
            Self::ESequenceLength => "ESEQUENCELENGTH",
            Self::ESequenceMismatch => "ESEQUENCEMISMATCH",
            Self::ESliceOption => "ESLICEOPTION",
            Self::ESpliceSite => "ESPLICESITE",
            Self::ESyntaxNested => "ESYNTAXNESTED",
            Self::ESyntaxUc => "ESYNTAXUC",
            Self::ESyntaxUeof => "ESYNTAXUEOF",
            Self::EUncertain => "EUNCERTAIN",
            Self::EVariantNotSupported => "EVARIANTNOTSUPPORTED",
        }
    }
}

/// Tag identifying which `FerroError` variant a mutalyzer code maps to.
///
/// `debug_tag()` returns the variant name as it appears in `format!("{e:?}")`
/// for an instance of the corresponding `FerroError` (e.g. `"ReferenceMismatch"`).
/// This is the substring the mutalyzer corpus runner asserts against when it
/// checks that a ferro error matches an expected mutalyzer code.
///
/// `is_lossy()` records that the mapping is approximate — the mutalyzer code
/// captures a finer-grained distinction than ferro's `FerroError` taxonomy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FerroErrorTag {
    name: &'static str,
    lossy: bool,
}

impl FerroErrorTag {
    /// The `FerroError` variant name as it appears in `Debug` output.
    pub fn debug_tag(self) -> &'static str {
        self.name
    }

    /// Returns true if this mapping is approximate (no exact 1:1 with the
    /// mutalyzer code's semantics).
    pub fn is_lossy(self) -> bool {
        self.lossy
    }
}

/// Mutalyzer codes that have no ferro equivalent.
///
/// Per-entry rationale follows. Codes split between mutalyzer's mapper /
/// batch-API layer (which ferro does not implement) and normalizer-surface
/// codes that ferro does not currently emit.
pub const NO_FERRO_EQUIV: &[&str] = &[
    // mapper API: CDS-slice description; no ferro normalizer path.
    "ECDSSLICES",
    // mapper API: malformed batch input; ferro has no batch-input surface.
    "EINVALIDINPUT",
    // mapper API: differential-slice length mismatch; ferro has no slice op.
    "ELENGTHSDIFFERENCE",
    // mapper API: location-slice request; ferro has no slice op.
    "ELOCATIONSLICE",
    // mapper API: zero inputs for a batch operation.
    "ENOINPUTS",
    // mapper API: insufficient inputs for a batch operation.
    "ENOINPUTSOTHER",
    // mapper API: sequence-length request; ferro has no equivalent.
    "ESEQUENCELENGTH",
    // mapper API: slice-option parameter error.
    "ESLICEOPTION",
    // normalizer-surface: splice-site predicate; ferro does not yet emit
    // splice diagnostics structurally.
    "ESPLICESITE",
    // normalizer-surface: `(?)` uncertain-position semantic concept ferro
    // does not model as a distinct error.
    "EUNCERTAIN",
];

/// `FerroError` variants with no mutalyzer counterpart.
///
/// These variants surface ferro-specific concerns (IO, JSON, projection
/// shape, ferro's split between transcript and protein reference lookups)
/// that mutalyzer's error taxonomy does not model.
pub const FERRO_NO_MUTALYZER_EQUIV: &[&str] = &[
    "AlignmentGap",
    "ConversionError",
    "ExonIntronBoundary",
    "GenomicReferenceNotAvailable",
    "Io",
    "Json",
    "ProteinReferenceNotAvailable",
    "ProteinSequenceUnavailable",
    "TranscriptNotOverlapping",
    "TranscriptVersionNotExact",
    "UnsupportedProjection",
    "UtrCdsBoundary",
    "VariantExceedsReference",
];

/// Forward map: mutalyzer code → ferro variant tag.
///
/// Returns `None` for codes in [`NO_FERRO_EQUIV`].
pub fn mutalyzer_to_ferro(code: &str) -> Option<FerroErrorTag> {
    let parsed = MutalyzerCode::parse(code)?;
    forward_map(parsed)
}

fn forward_map(code: MutalyzerCode) -> Option<FerroErrorTag> {
    use MutalyzerCode::*;
    Some(match code {
        // Exact-ish 1:1 mappings.
        EAminoAcidMismatch => exact("AminoAcidMismatch"),
        EMismatch => exact("ReferenceMismatch"),
        ERetr => exact("ReferenceNotFound"),
        EIntronic => exact("IntronicVariant"),
        ERangeReversed => exact("InvalidCoordinates"),
        ERepeatUnsupported => exact("UnsupportedVariant"),
        EVariantNotSupported => exact("UnsupportedVariant"),
        ESyntaxUc => exact("Parse"),
        ESyntaxUeof => exact("Parse"),
        ESyntaxNested => exact("Parse"),
        EPositionSyntax => exact("Parse"),
        EPositionInvalid => exact("Parse"),

        // Lossy: closest structural fit, semantics narrower or broader.
        ELengthMismatch => lossy("ReferenceMismatch"),
        ESequenceMismatch => lossy("ReferenceMismatch"),
        ERepeatMismatch => lossy("ReferenceMismatch"),
        ENoSelectorFound => lossy("ReferenceNotFound"),
        ESelectorOptions => lossy("ReferenceNotFound"),
        ENoToSelector => lossy("ReferenceNotFound"),
        ECoordinateSystemInvalid => lossy("Parse"),
        ENoCoordinateSystem => lossy("Parse"),
        // #543 added a parse-time `CoordinateSystemMismatch` rejection for
        // `c./g./m.` on a non-coding RNA (NR_/XR_), surfaced as a
        // `FerroError::Parse`. Map to "Parse" alongside the sibling
        // coordinate-system codes above (was the pre-#543 generic
        // "UnsupportedVariant", which no longer matches ferro's error).
        ECoordinateSystemMismatch => lossy("Parse"),
        EOffset => lossy("Parse"),
        EOutsideCds => lossy("InvalidCoordinates"),
        EIntronicRna => lossy("IntronicVariant"),
        // EOUTOFBOUNDARY → InvalidCoordinates (lossy). PositionOutOfBounds is
        // an ErrorCode-level distinction; the FerroError variant surface
        // collapses it into InvalidCoordinates.
        EOutOfBoundary => lossy("InvalidCoordinates"),
        EInsertionRange => lossy("InvalidCoordinates"),
        EInsertionLocation => lossy("InvalidCoordinates"),
        ERepeatReferenceLength => lossy("InvalidCoordinates"),
        EInsertedLength => lossy("InvalidCoordinates"),
        // EOVERLAP → InvalidCoordinates (lossy). ferro's
        // ErrorCode::SelfCancellingAllele (E3006) is structurally close but
        // lives in a different taxonomy; matching it would require touching
        // the FerroError enum.
        EOverlap => lossy("InvalidCoordinates"),
        ENoDna => lossy("Parse"),
        ENoRna => lossy("Parse"),
        EMissingParameter => lossy("Parse"),
        ENoCds => lossy("UnsupportedVariant"),

        // No ferro equivalent.
        ECdsSlices | EInvalidInput | ELengthsDifference | ELocationSlice | ENoInputs
        | ENoInputsOther | ESequenceLength | ESliceOption | ESpliceSite | EUncertain => {
            return None
        }
    })
}

fn exact(name: &'static str) -> FerroErrorTag {
    FerroErrorTag { name, lossy: false }
}

fn lossy(name: &'static str) -> FerroErrorTag {
    FerroErrorTag { name, lossy: true }
}

/// Reverse map: ferro error variant → all mutalyzer codes that may correspond.
///
/// Empty slice means there is no mutalyzer equivalent (see
/// [`FERRO_NO_MUTALYZER_EQUIV`]).
pub fn ferro_to_mutalyzer(err: &FerroError) -> &'static [&'static str] {
    match err {
        FerroError::Parse { .. } => &[
            "ECOORDINATESYSTEMINVALID",
            "ECOORDINATESYSTEMMISMATCH",
            "EMISSINGPARAMETER",
            "ENOCOORDINATESYSTEM",
            "ENODNA",
            "ENORNA",
            "EOFFSET",
            "EPOSITIONINVALID",
            "EPOSITIONSYNTAX",
            "ESYNTAXNESTED",
            "ESYNTAXUC",
            "ESYNTAXUEOF",
        ],
        FerroError::ReferenceNotFound { .. } => &[
            "ENOSELECTORFOUND",
            "ENOTOSELECTOR",
            "ERETR",
            "ESELECTOROPTIONS",
        ],
        FerroError::InvalidCoordinates { .. } => &[
            "EINSERTEDLENGTH",
            "EINSERTIONLOCATION",
            "EINSERTIONRANGE",
            "EOUTOFBOUNDARY",
            "EOUTSIDECDS",
            "EOVERLAP",
            "ERANGEREVERSED",
            "EREPEATREFERENCELENGTH",
        ],
        FerroError::UnsupportedVariant { .. } => {
            &["ENOCDS", "EREPEATUNSUPPORTED", "EVARIANTNOTSUPPORTED"]
        }
        FerroError::IntronicVariant { .. } => &["EINTRONIC", "EINTRONICRNA"],
        FerroError::ReferenceMismatch { .. } => &[
            "ELENGTHMISMATCH",
            "EMISMATCH",
            "EREPEATMISMATCH",
            "ESEQUENCEMISMATCH",
        ],
        FerroError::AminoAcidMismatch { .. } => &["EAMINOACIDMISMATCH"],

        // No mutalyzer equivalent — see FERRO_NO_MUTALYZER_EQUIV.
        FerroError::AlignmentGap { .. }
        | FerroError::ExonIntronBoundary { .. }
        | FerroError::UtrCdsBoundary { .. }
        | FerroError::GenomicReferenceNotAvailable { .. }
        | FerroError::ProteinReferenceNotAvailable { .. }
        | FerroError::ProteinSequenceUnavailable { .. }
        | FerroError::TranscriptNotOverlapping { .. }
        | FerroError::UnsupportedProjection { .. }
        | FerroError::ConversionError { .. }
        | FerroError::VariantExceedsReference { .. }
        | FerroError::TranscriptVersionNotExact { .. }
        | FerroError::Io { .. }
        | FerroError::Json { .. } => &[],
    }
}

/// Every known mutalyzer code (canonical upper-case form), for exhaustiveness
/// tests in downstream callers.
pub fn all_mutalyzer_codes() -> &'static [&'static str] {
    ALL_CODES
}

const ALL_CODES: &[&str] = &[
    "EAMINOACIDMISMATCH",
    "ECDSSLICES",
    "ECOORDINATESYSTEMINVALID",
    "ECOORDINATESYSTEMMISMATCH",
    "EINSERTEDLENGTH",
    "EINSERTIONLOCATION",
    "EINSERTIONRANGE",
    "EINTRONIC",
    "EINTRONICRNA",
    "EINVALIDINPUT",
    "ELENGTHMISMATCH",
    "ELENGTHSDIFFERENCE",
    "ELOCATIONSLICE",
    "EMISMATCH",
    "EMISSINGPARAMETER",
    "ENOCDS",
    "ENOCOORDINATESYSTEM",
    "ENODNA",
    "ENOINPUTS",
    "ENOINPUTSOTHER",
    "ENORNA",
    "ENOSELECTORFOUND",
    "ENOTOSELECTOR",
    "EOFFSET",
    "EOUTOFBOUNDARY",
    "EOUTSIDECDS",
    "EOVERLAP",
    "EPOSITIONINVALID",
    "EPOSITIONSYNTAX",
    "ERANGEREVERSED",
    "EREPEATMISMATCH",
    "EREPEATREFERENCELENGTH",
    "EREPEATUNSUPPORTED",
    "ERETR",
    "ESELECTOROPTIONS",
    "ESEQUENCELENGTH",
    "ESEQUENCEMISMATCH",
    "ESLICEOPTION",
    "ESPLICESITE",
    "ESYNTAXNESTED",
    "ESYNTAXUC",
    "ESYNTAXUEOF",
    "EUNCERTAIN",
    "EVARIANTNOTSUPPORTED",
];

#[cfg(test)]
mod tests {
    use super::*;

    /// Codes that appear in the upstream mutalyzer-normalize corpus's
    /// `errors` axis at the pinned commit. Sourced from grepping
    /// `tests/fixtures/mutalyzer-normalize/cases.json` and from the runner's
    /// recorded diagnostics for the `errors` axis. Used by
    /// [`forward_roundtrip_all_corpus_codes`] to guarantee that every code
    /// the runner sees in practice resolves to a non-`None` mapping.
    const CORPUS_ERRORS_AXIS_CODES: &[&str] = &[
        "ECOORDINATESYSTEMINVALID",
        "ECOORDINATESYSTEMMISMATCH",
        "EINSERTIONRANGE",
        "EINTRONIC",
        "EINTRONICRNA",
        "ELENGTHMISMATCH",
        "ENOCDS",
        "ENODNA",
        "ENOSELECTORFOUND",
        "EOUTOFBOUNDARY",
        "EOVERLAP",
        "EREPEATMISMATCH",
        "EREPEATUNSUPPORTED",
        "ERETR",
        "ESELECTOROPTIONS",
        "ESEQUENCEMISMATCH",
        "ESYNTAXNESTED",
        "ESYNTAXUEOF",
        "EVARIANTNOTSUPPORTED",
    ];

    /// Sample `FerroError` values, one per discriminant, for use by the
    /// reverse-map and Debug-substring tests.
    fn sample_ferro_errors() -> Vec<FerroError> {
        vec![
            FerroError::parse(0, "sample"),
            FerroError::AlignmentGap {
                msg: "sample".to_string(),
            },
            FerroError::ReferenceNotFound {
                id: "NM_TEST.1".to_string(),
            },
            FerroError::TranscriptVersionNotExact {
                requested: "NM_000088.4".to_string(),
            },
            FerroError::ExonIntronBoundary {
                exon: 1,
                variant: "c.1del".to_string(),
            },
            FerroError::UtrCdsBoundary {
                variant: "c.1del".to_string(),
            },
            FerroError::InvalidCoordinates {
                msg: "sample".to_string(),
            },
            FerroError::UnsupportedVariant {
                variant_type: "x".to_string(),
            },
            FerroError::IntronicVariant {
                variant: "c.1+1del".to_string(),
            },
            FerroError::GenomicReferenceNotAvailable {
                contig: "chr1".to_string(),
                start: 1,
                end: 2,
            },
            FerroError::ProteinReferenceNotAvailable {
                accession: "NP_TEST.1".to_string(),
                start: 1,
                end: 2,
            },
            FerroError::AminoAcidMismatch {
                accession: "NP_TEST.1".to_string(),
                position: 1,
                expected: "Ala".to_string(),
                found: "Gly".to_string(),
            },
            FerroError::ReferenceMismatch {
                location: "c.1".to_string(),
                expected: "A".to_string(),
                found: "G".to_string(),
            },
            FerroError::ConversionError {
                msg: "sample".to_string(),
            },
            FerroError::TranscriptNotOverlapping {
                variant: "g.1A>G".to_string(),
                transcript_id: "NM_TEST.1".to_string(),
            },
            FerroError::ProteinSequenceUnavailable {
                accession: "NP_TEST.1".to_string(),
            },
            FerroError::UnsupportedProjection {
                reason: "sample".to_string(),
            },
            FerroError::VariantExceedsReference {
                accession: "NC_000001.11".to_string(),
                hgvs_start: 1,
                hgvs_end: 10,
                expected_span: 10,
                actual_bytes: 5,
            },
            FerroError::Io {
                msg: "sample".to_string(),
            },
            FerroError::Json {
                msg: "sample".to_string(),
            },
        ]
    }

    /// The `Debug` representation of a `FerroError` always starts with the
    /// variant name followed by a delimiter (`{` for struct-variants, `(` for
    /// tuple variants, ` ` or end-of-string for unit variants). The corpus
    /// runner asserts that the tag substring appears in the `{:?}` output,
    /// which is satisfied so long as the tag is the variant name itself.
    fn debug_variant_name(err: &FerroError) -> &'static str {
        match err {
            FerroError::Parse { .. } => "Parse",
            FerroError::AlignmentGap { .. } => "AlignmentGap",
            FerroError::ReferenceNotFound { .. } => "ReferenceNotFound",
            FerroError::TranscriptVersionNotExact { .. } => "TranscriptVersionNotExact",
            FerroError::ExonIntronBoundary { .. } => "ExonIntronBoundary",
            FerroError::UtrCdsBoundary { .. } => "UtrCdsBoundary",
            FerroError::InvalidCoordinates { .. } => "InvalidCoordinates",
            FerroError::UnsupportedVariant { .. } => "UnsupportedVariant",
            FerroError::IntronicVariant { .. } => "IntronicVariant",
            FerroError::GenomicReferenceNotAvailable { .. } => "GenomicReferenceNotAvailable",
            FerroError::ProteinReferenceNotAvailable { .. } => "ProteinReferenceNotAvailable",
            FerroError::AminoAcidMismatch { .. } => "AminoAcidMismatch",
            FerroError::ReferenceMismatch { .. } => "ReferenceMismatch",
            FerroError::ConversionError { .. } => "ConversionError",
            FerroError::TranscriptNotOverlapping { .. } => "TranscriptNotOverlapping",
            FerroError::ProteinSequenceUnavailable { .. } => "ProteinSequenceUnavailable",
            FerroError::UnsupportedProjection { .. } => "UnsupportedProjection",
            FerroError::VariantExceedsReference { .. } => "VariantExceedsReference",
            FerroError::Io { .. } => "Io",
            FerroError::Json { .. } => "Json",
        }
    }

    /// Forward-map every code the upstream corpus's `errors` axis emits, and
    /// require a non-`None` mapping. If a new code shows up in a future corpus
    /// refresh, add it to `CORPUS_ERRORS_AXIS_CODES`. If the new code legitimately
    /// has no ferro equivalent, also list it in `NO_FERRO_EQUIV` and drop it
    /// from `CORPUS_ERRORS_AXIS_CODES`.
    #[test]
    fn forward_roundtrip_all_corpus_codes() {
        for code in CORPUS_ERRORS_AXIS_CODES {
            assert!(
                mutalyzer_to_ferro(code).is_some(),
                "corpus code {code:?} has no ferro mapping; either add a mapping or move it to NO_FERRO_EQUIV"
            );
        }
    }

    /// Every code in `all_mutalyzer_codes()` is either mapped to a ferro
    /// variant or appears in the explicit `NO_FERRO_EQUIV` allowlist.
    #[test]
    fn forward_all_upstream_codes_classified() {
        for code in all_mutalyzer_codes() {
            let mapped = mutalyzer_to_ferro(code).is_some();
            let allowlisted = NO_FERRO_EQUIV.contains(code);
            assert!(
                mapped ^ allowlisted,
                "code {code:?} must be either mapped or in NO_FERRO_EQUIV (exclusive); mapped={mapped}, allowlisted={allowlisted}"
            );
        }
    }

    /// Every `FerroError` discriminant either reverse-maps to one or more
    /// mutalyzer codes or appears in `FERRO_NO_MUTALYZER_EQUIV`.
    #[test]
    fn reverse_every_ferro_variant_has_codes_or_sentinel() {
        for err in sample_ferro_errors() {
            let codes = ferro_to_mutalyzer(&err);
            let name = debug_variant_name(&err);
            let allowlisted = FERRO_NO_MUTALYZER_EQUIV.contains(&name);
            assert!(
                !codes.is_empty() || allowlisted,
                "ferro variant {name:?} has no reverse mapping and is not in FERRO_NO_MUTALYZER_EQUIV"
            );
            // Codes returned must round-trip through the forward map.
            for c in codes {
                assert!(
                    mutalyzer_to_ferro(c).is_some(),
                    "reverse map for {name:?} returned {c:?} which has no forward mapping"
                );
            }
        }
    }

    /// `MutalyzerCode::parse` round-trips through `as_str()` for every variant
    /// listed in `all_mutalyzer_codes()`.
    #[test]
    fn mutalyzer_code_parse_roundtrip() {
        for code in all_mutalyzer_codes() {
            let parsed = MutalyzerCode::parse(code)
                .unwrap_or_else(|| panic!("known code {code:?} failed to parse"));
            assert_eq!(
                parsed.as_str(),
                *code,
                "round-trip failed for {code:?}: parsed.as_str() = {:?}",
                parsed.as_str()
            );
        }
    }

    /// For every (mutalyzer code → ferro tag) mapping, format the sample
    /// `FerroError` of that variant with `{:?}` and assert the tag is a
    /// substring. This is load-bearing for the corpus runner: its `axis_errors`
    /// test asserts `format!("{e:?}").contains(tag.debug_tag())`.
    #[test]
    fn tag_appears_in_ferro_debug_repr() {
        // Build a lookup: variant-name → sample FerroError.
        let samples = sample_ferro_errors();
        for code in all_mutalyzer_codes() {
            let Some(tag) = mutalyzer_to_ferro(code) else {
                continue;
            };
            let want = tag.debug_tag();
            let sample = samples
                .iter()
                .find(|e| debug_variant_name(e) == want)
                .unwrap_or_else(|| panic!("no sample FerroError for tag {want:?}"));
            let rendered = format!("{sample:?}");
            assert!(
                rendered.contains(want),
                "tag {want:?} for code {code:?} not found in {rendered:?}"
            );
        }
    }

    /// Unknown codes return `None` rather than panicking.
    #[test]
    fn unknown_code_returns_none() {
        assert!(mutalyzer_to_ferro("EBOGUS").is_none());
        assert!(mutalyzer_to_ferro("").is_none());
        assert!(mutalyzer_to_ferro("not-a-code").is_none());
    }

    /// Documented lossy mappings flip `is_lossy()`.
    #[test]
    fn lossy_flag_set_for_documented_lossy_mappings() {
        for code in &[
            "ECOORDINATESYSTEMMISMATCH",
            "ELENGTHMISMATCH",
            "EOUTOFBOUNDARY",
            "EOVERLAP",
            "ESEQUENCEMISMATCH",
        ] {
            let tag =
                mutalyzer_to_ferro(code).unwrap_or_else(|| panic!("expected mapping for {code:?}"));
            assert!(tag.is_lossy(), "expected {code:?} to be marked lossy");
        }
        // Spot-check: at least one exact mapping is NOT lossy.
        let exact = mutalyzer_to_ferro("EAMINOACIDMISMATCH").expect("EAMINOACIDMISMATCH maps");
        assert!(
            !exact.is_lossy(),
            "EAMINOACIDMISMATCH is an exact mapping, should not be lossy"
        );
    }
}
