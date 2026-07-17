//! Normalization configuration options

use crate::error_handling::{ErrorConfig, ErrorMode, ErrorOverride, ErrorType, ResolvedAction};
use serde::{Deserialize, Serialize};

/// Direction for variant shuffling during normalization
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ShuffleDirection {
    /// Shuffle towards 3' end (default for HGVS)
    #[default]
    ThreePrime,
    /// Shuffle towards 5' end (for VCF compatibility)
    FivePrime,
}

impl std::fmt::Display for ShuffleDirection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ShuffleDirection::ThreePrime => write!(f, "3prime"),
            ShuffleDirection::FivePrime => write!(f, "5prime"),
        }
    }
}

impl std::str::FromStr for ShuffleDirection {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "3prime" | "3'" | "three_prime" => Ok(ShuffleDirection::ThreePrime),
            "5prime" | "5'" | "five_prime" => Ok(ShuffleDirection::FivePrime),
            _ => Err(format!("Invalid shuffle direction: {}", s)),
        }
    }
}

/// Configuration for variant normalization
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NormalizeConfig {
    /// Direction to shuffle variants (default: 3')
    pub shuffle_direction: ShuffleDirection,

    /// Whether to allow crossing exon-intron boundaries
    pub cross_boundaries: bool,

    /// Error handling configuration (controls reference validation behavior)
    #[serde(skip)]
    pub error_config: ErrorConfig,

    /// Window size for reference sequence fetching
    pub window_size: u64,

    /// Vestigial overlap-prevention flag retained for source compatibility.
    ///
    /// Overlap is now prevented structurally: `normalize_allele` merges
    /// adjacent sub-variants before the per-variant pipeline runs, and
    /// `merge_consecutive_edits` uses a strict `prev.end + 1 == next.start`
    /// adjacency rule, so the normalizer cannot emit overlapping ranges
    /// from non-overlapping inputs. Detection of input-time overlap is
    /// handled unconditionally by `detect_overlap_conflicts`. This field
    /// has no effect and is preserved only so existing callers (notably
    /// `with_overlap_prevention`) keep compiling.
    pub prevent_overlap: bool,
}

impl Default for NormalizeConfig {
    fn default() -> Self {
        Self {
            shuffle_direction: ShuffleDirection::ThreePrime,
            cross_boundaries: false,
            // Default to Lenient mode for backwards compatibility
            // (previous behavior was to not validate at all)
            error_config: ErrorConfig::lenient(),
            window_size: 100,
            prevent_overlap: true,
        }
    }
}

impl PartialEq for NormalizeConfig {
    fn eq(&self, other: &Self) -> bool {
        self.shuffle_direction == other.shuffle_direction
            && self.cross_boundaries == other.cross_boundaries
            && self.window_size == other.window_size
            && self.prevent_overlap == other.prevent_overlap
    }
}

impl Eq for NormalizeConfig {}

impl NormalizeConfig {
    /// Create a new config with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a config with strict error handling (reject reference mismatches)
    pub fn strict() -> Self {
        Self {
            error_config: ErrorConfig::strict(),
            ..Default::default()
        }
    }

    /// Create a config with lenient error handling (warn on reference mismatches)
    pub fn lenient() -> Self {
        Self {
            error_config: ErrorConfig::lenient(),
            ..Default::default()
        }
    }

    /// Create a config with silent error handling (ignore reference mismatches)
    pub fn silent() -> Self {
        Self {
            error_config: ErrorConfig::silent(),
            ..Default::default()
        }
    }

    /// Set shuffle direction
    pub fn with_direction(mut self, direction: ShuffleDirection) -> Self {
        self.shuffle_direction = direction;
        self
    }

    /// Allow crossing boundaries
    pub fn allow_crossing_boundaries(mut self) -> Self {
        self.cross_boundaries = true;
        self
    }

    /// Set error handling mode
    pub fn with_error_mode(mut self, mode: ErrorMode) -> Self {
        self.error_config = ErrorConfig::new(mode);
        self
    }

    /// Replace the entire error-handling configuration.
    pub fn with_error_config(mut self, error_config: ErrorConfig) -> Self {
        self.error_config = error_config;
        self
    }

    /// Set a specific error type override
    pub fn with_error_override(mut self, error_type: ErrorType, action: ErrorOverride) -> Self {
        self.error_config = self.error_config.with_override(error_type, action);
        self
    }

    /// Disable reference validation (sets RefSeqMismatch to SilentCorrect)
    #[deprecated(
        since = "0.2.0",
        note = "Use with_error_mode(ErrorMode::Silent) instead"
    )]
    pub fn skip_validation(mut self) -> Self {
        self.error_config = self
            .error_config
            .with_override(ErrorType::RefSeqMismatch, ErrorOverride::SilentCorrect);
        self
    }

    /// Set the vestigial `prevent_overlap` flag.
    ///
    /// This builder is retained only for source compatibility with existing
    /// callers. The flag has no runtime effect: overlap prevention is handled
    /// structurally by `normalize_allele`, `merge_consecutive_edits`, and
    /// `detect_overlap_conflicts`. See [`NormalizeConfig::prevent_overlap`]
    /// for details.
    pub fn with_overlap_prevention(mut self, prevent: bool) -> Self {
        self.prevent_overlap = prevent;
        self
    }

    /// Get the resolved action for reference sequence mismatch
    pub fn ref_mismatch_action(&self) -> ResolvedAction {
        self.error_config.action_for(ErrorType::RefSeqMismatch)
    }

    /// Returns true if reference mismatches should be rejected
    pub fn should_reject_ref_mismatch(&self) -> bool {
        self.ref_mismatch_action().should_reject()
    }

    /// Returns true if reference mismatches should emit warnings
    pub fn should_warn_ref_mismatch(&self) -> bool {
        self.ref_mismatch_action().should_warn()
    }

    /// Get the resolved action for `VariantExceedsReference` (W5003) —
    /// fires when the provider returns fewer bytes than the HGVS
    /// interval span (the input violates HGVS spec refseq.md §43).
    /// Closes-after: #355.
    pub fn variant_exceeds_reference_action(&self) -> ResolvedAction {
        self.error_config
            .action_for(ErrorType::VariantExceedsReference)
    }

    /// Returns true if `VariantExceedsReference` should be rejected
    /// (strict mode default).
    pub fn should_reject_variant_exceeds_reference(&self) -> bool {
        self.variant_exceeds_reference_action().should_reject()
    }

    /// Returns true if `VariantExceedsReference` should emit a warning
    /// (lenient mode default; silent mode suppresses).
    pub fn should_warn_variant_exceeds_reference(&self) -> bool {
        self.variant_exceeds_reference_action().should_warn()
    }

    /// Get the resolved action for `PositionPastEnd` (W4004).
    pub fn position_past_end_action(&self) -> ResolvedAction {
        self.error_config.action_for(ErrorType::PositionPastEnd)
    }

    /// Returns true if past-end positions should be rejected (strict mode).
    pub fn should_reject_position_past_end(&self) -> bool {
        self.position_past_end_action().should_reject()
    }

    /// Returns true if a reduced-capability (no-genomic-data) degradation
    /// should be rejected — i.e. strict mode. Unlike the registry-backed
    /// spec warnings, `ReducedCapabilityNoGenome` is an *environmental*
    /// limitation rather than an input defect, so it is not user-overridable
    /// per errors-axis; it is simply promoted to an error in strict mode and
    /// surfaced as a warning-plus-best-effort otherwise (#1012 item 2).
    pub fn should_reject_reduced_capability(&self) -> bool {
        self.error_config.mode == ErrorMode::Strict
    }

    /// Returns true if past-end positions should emit warnings (lenient mode).
    pub fn should_warn_position_past_end(&self) -> bool {
        self.position_past_end_action().should_warn()
    }

    /// Get the resolved action for `IntronicOnBareTranscript` (W4007) —
    /// fires when an intronic offset appears on a bare transcript reference
    /// (`NM_` c. / `NR_` n. with `genomic_context: None`). See #486 EINTRONIC.
    pub fn intronic_bare_transcript_action(&self) -> ResolvedAction {
        self.error_config
            .action_for(ErrorType::IntronicOnBareTranscript)
    }

    /// Returns true if an intronic offset on a bare transcript should be
    /// rejected (strict mode default / errors-axis override).
    pub fn should_reject_intronic_bare_transcript(&self) -> bool {
        self.intronic_bare_transcript_action().should_reject()
    }

    /// Returns true if an intronic offset on a bare transcript should emit a
    /// warning (lenient mode default; silent suppresses).
    pub fn should_warn_intronic_bare_transcript(&self) -> bool {
        self.intronic_bare_transcript_action().should_warn()
    }

    /// Get the resolved action for `OverlapConflictingEdits` (W5002).
    pub fn overlap_conflict_action(&self) -> ResolvedAction {
        self.error_config
            .action_for(ErrorType::OverlapConflictingEdits)
    }

    /// Returns true if cis-allele edits with coincident reference
    /// bounds should be rejected (strict mode default). Closes #395
    /// item 6 — previously the `overlap.rs:88` emit site unconditionally
    /// pushed the warning, bypassing the registry's
    /// `always_warn_if_not_rejected` policy table that declared
    /// Strict→Reject.
    pub fn should_reject_overlap_conflict(&self) -> bool {
        self.overlap_conflict_action().should_reject()
    }

    /// Get the resolved action for `UnresolvableCentromere` (W4005).
    pub fn unresolvable_centromere_action(&self) -> ResolvedAction {
        self.error_config
            .action_for(ErrorType::UnresolvableCentromere)
    }

    /// Returns true if an unresolvable `cen` position should be rejected
    /// (strict mode default). A centromere is an assembly-annotated region
    /// with no sequence-derivable base, so it cannot be normalized; strict
    /// mode promotes the `UnresolvableSpecialPosition` warning to an error
    /// rather than silently echoing the input. See #488.
    pub fn should_reject_unresolvable_centromere(&self) -> bool {
        self.unresolvable_centromere_action().should_reject()
    }

    /// Get the resolved action for `TranscriptFlankNotDescribable` (W4006).
    pub fn transcript_flank_action(&self) -> ResolvedAction {
        self.error_config
            .action_for(ErrorType::TranscriptFlankNotDescribable)
    }

    /// Returns true if a telomere marker resolving to a transcript-flank
    /// position on a genomic-reference c. should be rejected (strict default).
    pub fn should_reject_transcript_flank(&self) -> bool {
        self.transcript_flank_action().should_reject()
    }

    /// Get the resolved action for `IncompleteCdsStartReference` (W5004) —
    /// fires when a `c.`/`p.`/`r.` variant is described against a transcript
    /// whose 5' CDS is annotated incomplete (`cds_start_NF`). See #972 Task 5.
    pub fn incomplete_cds_start_action(&self) -> ResolvedAction {
        self.error_config
            .action_for(ErrorType::IncompleteCdsStartReference)
    }

    /// Returns true if a `c.`/`p.`/`r.` variant against a `cds_start_NF`
    /// transcript should be rejected (strict mode default).
    pub fn should_reject_incomplete_cds_start(&self) -> bool {
        self.incomplete_cds_start_action().should_reject()
    }

    /// Returns true if a `c.`/`p.`/`r.` variant against a `cds_start_NF`
    /// transcript should emit a warning (lenient mode default; silent
    /// suppresses).
    pub fn should_warn_incomplete_cds_start(&self) -> bool {
        self.incomplete_cds_start_action().should_warn()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = NormalizeConfig::default();
        assert_eq!(config.shuffle_direction, ShuffleDirection::ThreePrime);
        assert!(!config.cross_boundaries);
        // Default is lenient (warn but don't reject)
        assert!(!config.should_reject_ref_mismatch());
        assert!(config.should_warn_ref_mismatch());
    }

    #[test]
    fn test_strict_config() {
        let config = NormalizeConfig::strict();
        assert!(config.should_reject_ref_mismatch());
        assert!(!config.should_warn_ref_mismatch());
    }

    #[test]
    fn test_lenient_config() {
        let config = NormalizeConfig::lenient();
        assert!(!config.should_reject_ref_mismatch());
        assert!(config.should_warn_ref_mismatch());
    }

    #[test]
    fn test_silent_config() {
        let config = NormalizeConfig::silent();
        assert!(!config.should_reject_ref_mismatch());
        assert!(!config.should_warn_ref_mismatch());
    }

    #[test]
    fn test_error_override() {
        // Start with lenient, override RefSeqMismatch to reject
        let config = NormalizeConfig::lenient()
            .with_error_override(ErrorType::RefSeqMismatch, ErrorOverride::Reject);
        assert!(config.should_reject_ref_mismatch());
    }

    #[test]
    fn test_direction_parsing() {
        assert_eq!(
            "3prime".parse::<ShuffleDirection>().unwrap(),
            ShuffleDirection::ThreePrime
        );
        assert_eq!(
            "5prime".parse::<ShuffleDirection>().unwrap(),
            ShuffleDirection::FivePrime
        );
    }

    #[test]
    #[allow(deprecated)]
    fn test_skip_validation_deprecated() {
        let config = NormalizeConfig::default().skip_validation();
        // skip_validation sets RefSeqMismatch to SilentCorrect
        assert!(!config.should_reject_ref_mismatch());
        assert!(!config.should_warn_ref_mismatch());
    }
}
