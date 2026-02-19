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

    /// Whether to prevent overlaps in compound variant normalization
    ///
    /// When true, normalization will check if variants in an allele would overlap
    /// after shifting and constrain the normalization to prevent collisions.
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

    /// Enable overlap prevention in compound variants
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
