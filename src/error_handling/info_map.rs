//! Mutalyzer ↔ ferro info-code mapping table.
//!
//! Companion to [`crate::error_handling::mutalyzer_map`] for the `infos`
//! axis. Mutalyzer's normalize output emits structured info codes (`I*`)
//! alongside the `errors` array; this module translates between those
//! upstream codes and ferro's [`crate::normalize::NormalizationInfo`]
//! variants for cross-tool diagnostics and the mutalyzer corpus runner
//! (`tests/mutalyzer_normalize_tests.rs`).
//!
//! Coverage policy mirrors `mutalyzer_map`: every `I*` code emitted by
//! mutalyzer/mutalyzer at the pinned upstream commit is classified — either
//! to a concrete `NormalizationInfo` variant tag or explicitly to "no ferro
//! equivalent" via [`NO_FERRO_INFO_EQUIV`]. The exhaustiveness test
//! [`tests::forward_all_upstream_info_codes_classified`] guarantees no
//! upstream code drifts into an unclassified state.

/// Strongly-typed mutalyzer info code.
///
/// Closed enum covering every `I*` code emitted by mutalyzer/mutalyzer at
/// the pinned upstream commit. New codes added upstream will fail to parse
/// via [`MutalyzerInfoCode::parse`] until added here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MutalyzerInfoCode {
    ICorrectedCoordinateSystem,
    ICorrectedLrgReference,
    ICorrectedPoint,
    ICorrectedSelectorId,
    ICorrectedVariantType,
    ILrgWarning,
    IMrnaGenomicTip,
    ISortedVariants,
    IWholeTranscriptExon,
}

impl MutalyzerInfoCode {
    /// Parse a mutalyzer info code string (e.g. `"ICORRECTEDPOINT"`).
    ///
    /// Returns `None` for unknown codes.
    pub fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "ICORRECTEDCOORDINATESYSTEM" => Self::ICorrectedCoordinateSystem,
            "ICORRECTEDLRGREFERENCE" => Self::ICorrectedLrgReference,
            "ICORRECTEDPOINT" => Self::ICorrectedPoint,
            "ICORRECTEDSELECTORID" => Self::ICorrectedSelectorId,
            "ICORRECTEDVARIANTTYPE" => Self::ICorrectedVariantType,
            "ILRGWARNING" => Self::ILrgWarning,
            "IMRNAGENOMICTIP" => Self::IMrnaGenomicTip,
            "ISORTEDVARIANTS" => Self::ISortedVariants,
            "IWHOLETRANSCRIPTEXON" => Self::IWholeTranscriptExon,
            _ => return None,
        })
    }

    /// Get the canonical upper-case mutalyzer string for this code.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ICorrectedCoordinateSystem => "ICORRECTEDCOORDINATESYSTEM",
            Self::ICorrectedLrgReference => "ICORRECTEDLRGREFERENCE",
            Self::ICorrectedPoint => "ICORRECTEDPOINT",
            Self::ICorrectedSelectorId => "ICORRECTEDSELECTORID",
            Self::ICorrectedVariantType => "ICORRECTEDVARIANTTYPE",
            Self::ILrgWarning => "ILRGWARNING",
            Self::IMrnaGenomicTip => "IMRNAGENOMICTIP",
            Self::ISortedVariants => "ISORTEDVARIANTS",
            Self::IWholeTranscriptExon => "IWHOLETRANSCRIPTEXON",
        }
    }
}

/// Tag identifying which [`crate::normalize::NormalizationInfo`] variant a
/// mutalyzer info code maps to.
///
/// `code()` returns the stable identifier string returned by
/// [`crate::normalize::NormalizationInfo::code`] for the corresponding
/// variant (e.g. `"SHUFFLE_APPLIED"`). The corpus runner uses
/// this for cross-tool equivalence checks: an expected mutalyzer info maps
/// to a ferro info code; the runner asserts the same code appears in the
/// emitted `infos` list.
///
/// `is_lossy()` records that the mapping is approximate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FerroInfoTag {
    code: &'static str,
    lossy: bool,
}

impl FerroInfoTag {
    /// The `NormalizationInfo` code string this tag identifies.
    pub fn code(self) -> &'static str {
        self.code
    }

    /// Returns true if this mapping is approximate (no exact 1:1 with the
    /// mutalyzer code's semantics).
    pub fn is_lossy(self) -> bool {
        self.lossy
    }
}

/// Mutalyzer info codes that have no ferro equivalent.
///
/// Per-entry rationale follows. ferro does not currently model these
/// signals as structured info codes; some are mutalyzer-specific
/// behaviours (LRG handling, auto coordinate-system inference) that ferro
/// intentionally treats as errors or no-ops instead.
pub const NO_FERRO_INFO_EQUIV: &[&str] = &[
    // mutalyzer auto-infers a missing coordinate system prefix
    // (e.g. `NM_xxx:274del` → `NM_xxx:c.274del`). ferro errors on missing
    // prefixes rather than silently correcting, so there is no info to
    // emit.
    "ICORRECTEDCOORDINATESYSTEM",
    // LRG-specific reference correction; ferro does not implement the LRG
    // reference system.
    "ICORRECTEDLRGREFERENCE",
    // mutalyzer corrects a missing/wrong selector ID on the input. ferro
    // rewrites a coding-transcript selector to its protein accession for a
    // `p.` coordinate on a genomic reference (NM_→NP_, #502) but emits no
    // structured info for the rewrite; for other selector classes it requires
    // explicit selectors. Either way there is no ferro info to map.
    "ICORRECTEDSELECTORID",
    // mutalyzer corrects the variant type (e.g. del→delins) when the
    // stated bases force the change. ferro canonicalizes via
    // `canonicalize_edit` but does not emit a structured info for the
    // rewrite.
    "ICORRECTEDVARIANTTYPE",
    // LRG-specific warning surface; see `ICORRECTEDLRGREFERENCE`.
    "ILRGWARNING",
    // mutalyzer recommends adding a genomic tip ("g." accession) when a
    // c./n. variant is reported alone; ferro has no equivalent
    // recommendation surface.
    "IMRNAGENOMICTIP",
    // mutalyzer reorders cis-allele members by position and emits this
    // info. ferro's `normalize::merge` layer combines positionally-
    // adjacent edits but does not sort members; if a future change does
    // start sorting, this rationale (and the surface) will need
    // revisiting.
    "ISORTEDVARIANTS",
    // mutalyzer flags transcripts where the entire transcript is a single
    // exon. ferro has no exon-classification info surface.
    "IWHOLETRANSCRIPTEXON",
];

/// Forward map: mutalyzer info code → ferro info-variant tag.
///
/// Returns `None` for codes in [`NO_FERRO_INFO_EQUIV`].
pub fn mutalyzer_info_to_ferro(code: &str) -> Option<FerroInfoTag> {
    let parsed = MutalyzerInfoCode::parse(code)?;
    forward_info_map(parsed)
}

fn forward_info_map(code: MutalyzerInfoCode) -> Option<FerroInfoTag> {
    use MutalyzerInfoCode::*;
    Some(match code {
        // Exact-ish 1:1 mapping: mutalyzer's `ICORRECTEDPOINT` fires when
        // a position is moved by the 3'-rule shuffle. ferro emits the
        // direction-agnostic `SHUFFLE_APPLIED` code carrying the actual
        // direction in the variant payload — under the default config
        // (`ShuffleDirection::ThreePrime`) the two signals coincide.
        // Mapping is lossy because ferro may also emit `SHUFFLE_APPLIED`
        // for `FivePrime` configurations that mutalyzer never produces.
        ICorrectedPoint => lossy("SHUFFLE_APPLIED"),

        // No ferro equivalent — see NO_FERRO_INFO_EQUIV.
        ICorrectedCoordinateSystem
        | ICorrectedLrgReference
        | ICorrectedSelectorId
        | ICorrectedVariantType
        | ILrgWarning
        | IMrnaGenomicTip
        | ISortedVariants
        | IWholeTranscriptExon => return None,
    })
}

#[allow(dead_code)]
fn exact(code: &'static str) -> FerroInfoTag {
    FerroInfoTag { code, lossy: false }
}

fn lossy(code: &'static str) -> FerroInfoTag {
    FerroInfoTag { code, lossy: true }
}

/// Every known mutalyzer info code (canonical upper-case form), for
/// exhaustiveness tests in downstream callers.
pub fn all_mutalyzer_info_codes() -> &'static [&'static str] {
    ALL_INFO_CODES
}

const ALL_INFO_CODES: &[&str] = &[
    "ICORRECTEDCOORDINATESYSTEM",
    "ICORRECTEDLRGREFERENCE",
    "ICORRECTEDPOINT",
    "ICORRECTEDSELECTORID",
    "ICORRECTEDVARIANTTYPE",
    "ILRGWARNING",
    "IMRNAGENOMICTIP",
    "ISORTEDVARIANTS",
    "IWHOLETRANSCRIPTEXON",
];

#[cfg(test)]
mod tests {
    use super::*;

    /// Every code listed in `all_mutalyzer_info_codes()` is either mapped
    /// to a ferro info variant or appears in the `NO_FERRO_INFO_EQUIV`
    /// allowlist (exclusive). This is the load-bearing exhaustiveness
    /// guarantee — a future corpus refresh that adds a new `I*` code will
    /// trip this test until it is classified.
    #[test]
    fn forward_all_upstream_info_codes_classified() {
        for code in all_mutalyzer_info_codes() {
            let mapped = mutalyzer_info_to_ferro(code).is_some();
            let allowlisted = NO_FERRO_INFO_EQUIV.contains(code);
            assert!(
                mapped ^ allowlisted,
                "info code {code:?} must be either mapped or in NO_FERRO_INFO_EQUIV (exclusive); mapped={mapped}, allowlisted={allowlisted}"
            );
        }
    }

    /// `MutalyzerInfoCode::parse` round-trips through `as_str()` for every
    /// variant listed in `all_mutalyzer_info_codes()`.
    #[test]
    fn mutalyzer_info_code_parse_roundtrip() {
        for code in all_mutalyzer_info_codes() {
            let parsed = MutalyzerInfoCode::parse(code)
                .unwrap_or_else(|| panic!("known info code {code:?} failed to parse"));
            assert_eq!(
                parsed.as_str(),
                *code,
                "round-trip failed for {code:?}: parsed.as_str() = {:?}",
                parsed.as_str()
            );
        }
    }

    /// Unknown info codes return `None` rather than panicking.
    #[test]
    fn unknown_info_code_returns_none() {
        assert!(mutalyzer_info_to_ferro("IBOGUS").is_none());
        assert!(mutalyzer_info_to_ferro("").is_none());
        assert!(mutalyzer_info_to_ferro("ECORRECTEDPOINT").is_none()); // wrong prefix
    }

    /// `ICORRECTEDPOINT` maps to the ferro shuffle-applied info code that
    /// `NormalizationInfo::ShuffleApplied::code()` returns. Pin the exact
    /// string so the corpus runner's check is stable. Marked lossy
    /// because ferro's signal carries an explicit direction
    /// (3'/5') while mutalyzer's only ever fires for 3'-shift.
    #[test]
    fn icorrectedpoint_maps_to_shuffle_applied() {
        let tag = mutalyzer_info_to_ferro("ICORRECTEDPOINT")
            .expect("ICORRECTEDPOINT should map to a ferro info variant");
        assert_eq!(tag.code(), "SHUFFLE_APPLIED");
        assert!(
            tag.is_lossy(),
            "ICORRECTEDPOINT → SHUFFLE_APPLIED is lossy (direction-agnostic)"
        );
    }
}
