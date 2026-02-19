//! HGVS nomenclature version tracking
//!
//! This module tracks HGVS nomenclature versions and provides information
//! about version-specific features and deprecated syntax.
//!
//! Current supported version: HGVS Nomenclature 20.05 (May 2020)
//! Specification: <https://varnomen.hgvs.org/>

use std::fmt;

/// HGVS nomenclature version
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HgvsVersion {
    /// Year (e.g., 20 for 2020)
    pub year: u16,
    /// Month (e.g., 05 for May)
    pub month: u8,
}

impl HgvsVersion {
    /// Create a new HGVS version
    pub const fn new(year: u16, month: u8) -> Self {
        Self { year, month }
    }

    /// HGVS 15.11 (November 2015) - First versioned recommendation
    pub const V15_11: Self = Self::new(15, 11);

    /// HGVS 19.01 (January 2019)
    pub const V19_01: Self = Self::new(19, 1);

    /// HGVS 20.05 (May 2020) - Current version
    pub const V20_05: Self = Self::new(20, 5);

    /// Current version supported by this library
    pub const CURRENT: Self = Self::V20_05;

    /// Check if this version is at least as new as another
    pub fn at_least(&self, other: Self) -> bool {
        *self >= other
    }
}

impl Default for HgvsVersion {
    fn default() -> Self {
        Self::CURRENT
    }
}

impl fmt::Display for HgvsVersion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:02}.{:02}", self.year, self.month)
    }
}

/// Feature introduced in a specific HGVS version
#[derive(Debug, Clone)]
pub struct VersionedFeature {
    /// Feature name
    pub name: &'static str,
    /// Description
    pub description: &'static str,
    /// Version introduced
    pub introduced: HgvsVersion,
    /// Version deprecated (if any)
    pub deprecated: Option<HgvsVersion>,
}

impl VersionedFeature {
    /// Create a new versioned feature
    pub const fn new(
        name: &'static str,
        description: &'static str,
        introduced: HgvsVersion,
    ) -> Self {
        Self {
            name,
            description,
            introduced,
            deprecated: None,
        }
    }

    /// Create a deprecated feature
    pub const fn deprecated(
        name: &'static str,
        description: &'static str,
        introduced: HgvsVersion,
        deprecated: HgvsVersion,
    ) -> Self {
        Self {
            name,
            description,
            introduced,
            deprecated: Some(deprecated),
        }
    }

    /// Check if this feature is supported in a version
    pub fn is_supported(&self, version: HgvsVersion) -> bool {
        version.at_least(self.introduced)
    }

    /// Check if this feature is deprecated in a version
    pub fn is_deprecated(&self, version: HgvsVersion) -> bool {
        self.deprecated.is_some_and(|dep| version.at_least(dep))
    }
}

/// Known HGVS features with version information
pub mod features {
    use super::*;

    /// Protein extension notation (p.Met1ext-5)
    pub const PROTEIN_EXTENSION: VersionedFeature = VersionedFeature::new(
        "protein_extension",
        "Protein N-terminal and C-terminal extensions",
        HgvsVersion::V15_11,
    );

    /// Repeated sequence notation with brackets
    pub const REPEAT_NOTATION: VersionedFeature = VersionedFeature::new(
        "repeat_notation",
        "Repeated sequence notation with [n] syntax",
        HgvsVersion::V15_11,
    );

    /// Allele notation with brackets and semicolons
    pub const ALLELE_NOTATION: VersionedFeature = VersionedFeature::new(
        "allele_notation",
        "Compound heterozygous notation [a];[b]",
        HgvsVersion::V15_11,
    );

    /// Mosaic notation with /
    pub const MOSAIC_NOTATION: VersionedFeature = VersionedFeature::new(
        "mosaic_notation",
        "Mosaic variant notation with /",
        HgvsVersion::V15_11,
    );

    /// Chimeric notation with //
    pub const CHIMERIC_NOTATION: VersionedFeature = VersionedFeature::new(
        "chimeric_notation",
        "Chimeric variant notation with //",
        HgvsVersion::V15_11,
    );

    /// Uncertain positions with parentheses
    pub const UNCERTAIN_POSITIONS: VersionedFeature = VersionedFeature::new(
        "uncertain_positions",
        "Uncertain positions with (pos) notation",
        HgvsVersion::V15_11,
    );

    /// Reference sequence required
    pub const REFERENCE_REQUIRED: VersionedFeature = VersionedFeature::new(
        "reference_required",
        "Reference sequence mandatory in all descriptions",
        HgvsVersion::V15_11,
    );

    /// Gene conversion notation
    pub const GENE_CONVERSION: VersionedFeature = VersionedFeature::new(
        "gene_conversion",
        "Gene conversion notation with con prefix",
        HgvsVersion::V15_11,
    );
}

/// Deprecated syntax patterns
pub mod deprecated {
    use super::*;

    /// Use of X for stop codon (use Ter instead)
    pub const STOP_AS_X: VersionedFeature = VersionedFeature::deprecated(
        "stop_as_x",
        "Using X for stop codon - use Ter instead",
        HgvsVersion::V15_11,
        HgvsVersion::V15_11,
    );

    /// Use of * for stop codon in protein (use Ter instead)
    pub const STOP_AS_STAR: VersionedFeature = VersionedFeature::deprecated(
        "stop_as_star",
        "Using * for stop codon - use Ter instead",
        HgvsVersion::V15_11,
        HgvsVersion::V15_11,
    );

    /// IVS notation for introns (use c.123+1 instead)
    pub const IVS_NOTATION: VersionedFeature = VersionedFeature::deprecated(
        "ivs_notation",
        "IVS notation for introns - use intronic offsets instead",
        HgvsVersion::V15_11,
        HgvsVersion::V15_11,
    );

    /// Old frameshift notation (fsX instead of fsTer)
    pub const FS_X_NOTATION: VersionedFeature = VersionedFeature::deprecated(
        "fs_x_notation",
        "fsX notation - use fsTer instead",
        HgvsVersion::V15_11,
        HgvsVersion::V15_11,
    );
}

/// Get the current supported HGVS version
pub fn current_version() -> HgvsVersion {
    HgvsVersion::CURRENT
}

/// Check if a feature is supported in the current version
pub fn is_feature_supported(feature: &VersionedFeature) -> bool {
    feature.is_supported(HgvsVersion::CURRENT)
}

/// Get all known features
pub fn known_features() -> Vec<&'static VersionedFeature> {
    vec![
        &features::PROTEIN_EXTENSION,
        &features::REPEAT_NOTATION,
        &features::ALLELE_NOTATION,
        &features::MOSAIC_NOTATION,
        &features::CHIMERIC_NOTATION,
        &features::UNCERTAIN_POSITIONS,
        &features::REFERENCE_REQUIRED,
        &features::GENE_CONVERSION,
    ]
}

/// Get all deprecated features
pub fn deprecated_features() -> Vec<&'static VersionedFeature> {
    vec![
        &deprecated::STOP_AS_X,
        &deprecated::STOP_AS_STAR,
        &deprecated::IVS_NOTATION,
        &deprecated::FS_X_NOTATION,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version_comparison() {
        assert!(HgvsVersion::V20_05 > HgvsVersion::V19_01);
        assert!(HgvsVersion::V19_01 > HgvsVersion::V15_11);
        assert!(HgvsVersion::V20_05.at_least(HgvsVersion::V15_11));
    }

    #[test]
    fn test_version_display() {
        assert_eq!(HgvsVersion::V20_05.to_string(), "20.05");
        assert_eq!(HgvsVersion::V15_11.to_string(), "15.11");
    }

    #[test]
    fn test_current_version() {
        assert_eq!(current_version(), HgvsVersion::V20_05);
    }

    #[test]
    fn test_feature_supported() {
        assert!(features::PROTEIN_EXTENSION.is_supported(HgvsVersion::V20_05));
        assert!(features::PROTEIN_EXTENSION.is_supported(HgvsVersion::V15_11));
        assert!(!features::PROTEIN_EXTENSION.is_supported(HgvsVersion::new(14, 1)));
    }

    #[test]
    fn test_deprecated_feature() {
        assert!(deprecated::STOP_AS_X.is_deprecated(HgvsVersion::V20_05));
        assert!(deprecated::IVS_NOTATION.is_deprecated(HgvsVersion::CURRENT));
    }

    #[test]
    fn test_known_features() {
        let features = known_features();
        assert!(!features.is_empty());
        assert!(features.iter().all(|f| is_feature_supported(f)));
    }

    #[test]
    fn test_deprecated_features() {
        let deps = deprecated_features();
        assert!(!deps.is_empty());
    }
}
