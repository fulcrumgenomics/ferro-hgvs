//! Legacy format support for deprecated HGVS notation.
//!
//! This module provides support for parsing deprecated and legacy HGVS notation
//! formats that are still commonly encountered in clinical databases and literature.
//!
//! # Supported Legacy Formats
//!
//! - **IVS notation**: `c.IVS4+1G>A` → `c.459+1G>A`
//! - **Old substitution syntax**: `c.100_102>ATG` → `c.100_102delinsATG`
//! - **Arrow protein substitution**: `p.V600>E` → `p.Val600Glu`
//! - **Number-first protein**: `p.600V>E` → `p.Val600Glu`
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::legacy::{LegacyParser, LegacyConfig};
//!
//! // Parse a modern variant (passes through)
//! let parser = LegacyParser::new(LegacyConfig::default());
//! let result = parser.parse("NM_000088.3:c.459del").unwrap();
//! assert!(!result.is_legacy());
//!
//! // Parse a legacy protein variant (with accession)
//! let result = parser.parse("NP_000079.2:p.V600>E").unwrap();
//! assert!(result.is_legacy());
//! assert_eq!(result.modern_notation, Some("NP_000079.2:p.Val600Glu".to_string()));
//! ```

mod ivs;
mod protein;
mod substitution;

pub use ivs::{parse_ivs, IvsPosition};
pub use protein::{parse_legacy_protein, LegacyProteinFormat};
pub use substitution::{parse_old_substitution, OldSubstitutionFormat};

use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;

/// Legacy format type detected during parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LegacyFormat {
    /// IVS notation (e.g., `c.IVS4+1G>A`).
    IvsNotation,
    /// Old substitution with arrow (e.g., `c.100_102>ATG`).
    OldSubstitutionArrow,
    /// Single position insertion (e.g., `c.100insA`).
    SinglePositionInsertion,
    /// Number-first protein notation (e.g., `p.600V>E`).
    NumberFirstProtein,
    /// Arrow in protein substitution (e.g., `p.V600>E`).
    ArrowProteinSubstitution,
    /// Chromosome-only format (e.g., `chr1:12345:A:G`).
    ChromosomeOnly,
    /// Unversioned accession (e.g., `NM_000088:c.459del`).
    UnversionedAccession,
}

impl LegacyFormat {
    /// Get the deprecation warning message for this format.
    pub fn deprecation_message(&self) -> &'static str {
        match self {
            LegacyFormat::IvsNotation => {
                "IVS notation is deprecated; use c.N+M format (e.g., c.459+1G>A)"
            }
            LegacyFormat::OldSubstitutionArrow => {
                "> for multi-base substitution is deprecated; use delins syntax"
            }
            LegacyFormat::SinglePositionInsertion => {
                "Single position insertion is deprecated; use two positions (e.g., c.100_101insA)"
            }
            LegacyFormat::NumberFirstProtein => {
                "Number-first protein notation is deprecated; use Aa###Aa format"
            }
            LegacyFormat::ArrowProteinSubstitution => {
                "Arrow in protein substitution is deprecated; use direct notation (e.g., p.Val600Glu)"
            }
            LegacyFormat::ChromosomeOnly => {
                "Chromosome-only format is non-standard; use NC_ accession with g. prefix"
            }
            LegacyFormat::UnversionedAccession => {
                "Unversioned accession is deprecated; specify version (e.g., NM_000088.3)"
            }
        }
    }
}

/// Result of parsing a legacy format.
#[derive(Debug, Clone)]
pub struct LegacyParseResult {
    /// Original input string.
    pub original: String,
    /// Parsed variant (if successful).
    pub parsed: Option<HgvsVariant>,
    /// Legacy formats detected.
    pub legacy_formats: Vec<LegacyFormat>,
    /// Deprecation warnings.
    pub warnings: Vec<String>,
    /// Modern HGVS notation (if converted).
    pub modern_notation: Option<String>,
}

impl LegacyParseResult {
    /// Create a new result for a successfully parsed modern variant.
    pub fn modern(original: String, variant: HgvsVariant) -> Self {
        let modern = format!("{}", variant);
        Self {
            original,
            parsed: Some(variant),
            legacy_formats: vec![],
            warnings: vec![],
            modern_notation: Some(modern),
        }
    }

    /// Create a new result for a legacy format.
    pub fn legacy(
        original: String,
        variant: HgvsVariant,
        formats: Vec<LegacyFormat>,
        warn: bool,
    ) -> Self {
        let modern = format!("{}", variant);
        let warnings = if warn {
            formats
                .iter()
                .map(|f| f.deprecation_message().to_string())
                .collect()
        } else {
            vec![]
        };

        Self {
            original,
            parsed: Some(variant),
            legacy_formats: formats,
            warnings,
            modern_notation: Some(modern),
        }
    }

    /// Check if any legacy formats were detected.
    pub fn is_legacy(&self) -> bool {
        !self.legacy_formats.is_empty()
    }

    /// Check if parsing was successful.
    pub fn is_ok(&self) -> bool {
        self.parsed.is_some()
    }
}

/// Configuration for legacy parsing.
#[derive(Debug, Clone)]
pub struct LegacyConfig {
    /// Enable IVS notation parsing.
    pub parse_ivs: bool,
    /// Enable old substitution syntax.
    pub parse_old_substitution: bool,
    /// Enable legacy protein notation.
    pub parse_legacy_protein: bool,
    /// Enable chromosome-only format.
    pub parse_chromosome_only: bool,
    /// Emit deprecation warnings.
    pub warn_deprecated: bool,
    /// Auto-convert to modern format.
    pub auto_convert: bool,
}

impl Default for LegacyConfig {
    fn default() -> Self {
        Self {
            parse_ivs: true,
            parse_old_substitution: true,
            parse_legacy_protein: true,
            parse_chromosome_only: true,
            warn_deprecated: true,
            auto_convert: true,
        }
    }
}

impl LegacyConfig {
    /// Create a config that only parses modern formats.
    pub fn modern_only() -> Self {
        Self {
            parse_ivs: false,
            parse_old_substitution: false,
            parse_legacy_protein: false,
            parse_chromosome_only: false,
            warn_deprecated: false,
            auto_convert: false,
        }
    }

    /// Create a config that parses all formats without warnings.
    pub fn permissive() -> Self {
        Self {
            parse_ivs: true,
            parse_old_substitution: true,
            parse_legacy_protein: true,
            parse_chromosome_only: true,
            warn_deprecated: false,
            auto_convert: true,
        }
    }
}

/// Exon boundary data for IVS conversion.
#[derive(Debug, Clone)]
pub struct ExonData {
    /// Exon boundaries as (cds_start, cds_end) pairs.
    pub exons: Vec<(i64, i64)>,
}

impl ExonData {
    /// Create new exon data.
    pub fn new(exons: Vec<(i64, i64)>) -> Self {
        Self { exons }
    }

    /// Get the CDS position for an IVS notation.
    ///
    /// IVS N refers to the intron between exon N and exon N+1.
    /// - `IVS4+1` = 1 base after end of exon 4
    /// - `IVS4-2` = 2 bases before start of exon 5
    pub fn ivs_to_cds(&self, ivs: &IvsPosition) -> Result<(i64, Option<i64>), FerroError> {
        if ivs.intron == 0 || ivs.intron as usize > self.exons.len() {
            return Err(FerroError::InvalidCoordinates {
                msg: format!("Invalid intron number: IVS{}", ivs.intron),
            });
        }

        let exon_idx = ivs.intron as usize - 1;

        if ivs.offset > 0 {
            // Donor side: relative to end of exon N
            let exon_end = self.exons[exon_idx].1;
            Ok((exon_end, Some(ivs.offset)))
        } else if ivs.offset < 0 {
            // Acceptor side: relative to start of exon N+1
            if exon_idx + 1 >= self.exons.len() {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!("No exon after IVS{}", ivs.intron),
                });
            }
            let next_exon_start = self.exons[exon_idx + 1].0;
            Ok((next_exon_start, Some(ivs.offset)))
        } else {
            Err(FerroError::InvalidCoordinates {
                msg: "IVS offset cannot be 0".to_string(),
            })
        }
    }
}

/// Parser for legacy HGVS formats.
#[derive(Debug, Clone)]
pub struct LegacyParser {
    config: LegacyConfig,
    exon_data: Option<ExonData>,
}

impl LegacyParser {
    /// Create a new legacy parser with the given configuration.
    pub fn new(config: LegacyConfig) -> Self {
        Self {
            config,
            exon_data: None,
        }
    }

    /// Create a legacy parser with exon data for IVS conversion.
    pub fn with_exon_data(config: LegacyConfig, exon_data: ExonData) -> Self {
        Self {
            config,
            exon_data: Some(exon_data),
        }
    }

    /// Parse a variant string, handling legacy formats.
    pub fn parse(&self, input: &str) -> Result<LegacyParseResult, FerroError> {
        let input = input.trim();

        // Try modern parser first
        if let Ok(variant) = parse_hgvs(input) {
            return Ok(LegacyParseResult::modern(input.to_string(), variant));
        }

        // Try legacy protein formats
        if self.config.parse_legacy_protein {
            if let Some(result) = self.try_legacy_protein(input)? {
                return Ok(result);
            }
        }

        // Try old substitution format
        if self.config.parse_old_substitution {
            if let Some(result) = self.try_old_substitution(input)? {
                return Ok(result);
            }
        }

        // If we have exon data and IVS parsing is enabled, try that
        if self.config.parse_ivs {
            if let Some(result) = self.try_ivs_notation(input)? {
                return Ok(result);
            }
        }

        // Nothing worked
        Err(FerroError::Parse {
            pos: 0,
            msg: format!("Failed to parse variant: {}", input),
            diagnostic: None,
        })
    }

    /// Try to parse as legacy protein format.
    fn try_legacy_protein(&self, input: &str) -> Result<Option<LegacyParseResult>, FerroError> {
        // Check if it looks like a protein variant
        if !input.contains("p.") {
            return Ok(None);
        }

        // Extract the protein part
        let parts: Vec<&str> = input.splitn(2, "p.").collect();
        if parts.len() != 2 {
            return Ok(None);
        }

        let prefix = parts[0];
        let prot_part = parts[1];

        // Try legacy formats
        if let Some((legacy_format, converted)) = protein::convert_legacy_protein(prot_part) {
            // Reconstruct and try to parse
            let modern_str = format!("{}p.{}", prefix, converted);

            if let Ok(variant) = parse_hgvs(&modern_str) {
                return Ok(Some(LegacyParseResult::legacy(
                    input.to_string(),
                    variant,
                    vec![legacy_format],
                    self.config.warn_deprecated,
                )));
            }
        }

        Ok(None)
    }

    /// Try to parse as old substitution format.
    fn try_old_substitution(&self, input: &str) -> Result<Option<LegacyParseResult>, FerroError> {
        // Look for patterns like c.100_102>ATG (should be c.100_102delinsATG)
        if let Some(converted) = substitution::convert_old_substitution(input) {
            if let Ok(variant) = parse_hgvs(&converted) {
                return Ok(Some(LegacyParseResult::legacy(
                    input.to_string(),
                    variant,
                    vec![LegacyFormat::OldSubstitutionArrow],
                    self.config.warn_deprecated,
                )));
            }
        }

        Ok(None)
    }

    /// Try to parse IVS notation.
    fn try_ivs_notation(&self, input: &str) -> Result<Option<LegacyParseResult>, FerroError> {
        // Check if input contains IVS
        if !input.contains("IVS") {
            return Ok(None);
        }

        // We need exon data to convert IVS to CDS positions
        let exon_data = match &self.exon_data {
            Some(data) => data,
            None => {
                return Err(FerroError::InvalidCoordinates {
                    msg: "Exon data required for IVS notation conversion".to_string(),
                });
            }
        };

        // Try to convert IVS notation
        if let Some(converted) = ivs::convert_ivs_notation(input, exon_data)? {
            if let Ok(variant) = parse_hgvs(&converted) {
                return Ok(Some(LegacyParseResult::legacy(
                    input.to_string(),
                    variant,
                    vec![LegacyFormat::IvsNotation],
                    self.config.warn_deprecated,
                )));
            }
        }

        Ok(None)
    }
}

/// Parse a variant with legacy format support.
///
/// Convenience function that uses the default legacy configuration.
pub fn parse_with_legacy(input: &str) -> Result<LegacyParseResult, FerroError> {
    let parser = LegacyParser::new(LegacyConfig::default());
    parser.parse(input)
}

/// Parse a variant with legacy format support and exon data.
///
/// Use this when you need to convert IVS notation to modern format.
pub fn parse_with_legacy_and_exons(
    input: &str,
    exon_data: ExonData,
) -> Result<LegacyParseResult, FerroError> {
    let parser = LegacyParser::with_exon_data(LegacyConfig::default(), exon_data);
    parser.parse(input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modern_variant_passthrough() {
        let parser = LegacyParser::new(LegacyConfig::default());
        let result = parser.parse("NM_000088.3:c.459del").unwrap();

        assert!(!result.is_legacy());
        assert!(result.is_ok());
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn test_legacy_format_deprecation_messages() {
        assert!(!LegacyFormat::IvsNotation.deprecation_message().is_empty());
        assert!(!LegacyFormat::OldSubstitutionArrow
            .deprecation_message()
            .is_empty());
        assert!(!LegacyFormat::ArrowProteinSubstitution
            .deprecation_message()
            .is_empty());
    }

    #[test]
    fn test_legacy_config_default() {
        let config = LegacyConfig::default();
        assert!(config.parse_ivs);
        assert!(config.parse_legacy_protein);
        assert!(config.warn_deprecated);
    }

    #[test]
    fn test_legacy_config_modern_only() {
        let config = LegacyConfig::modern_only();
        assert!(!config.parse_ivs);
        assert!(!config.parse_legacy_protein);
    }

    #[test]
    fn test_exon_data_ivs_to_cds() {
        // Create exon data: exon 1 at 1-100, exon 2 at 201-300
        let exon_data = ExonData::new(vec![(1, 100), (201, 300)]);

        // IVS1+5 should be position 100 with offset +5
        let ivs = IvsPosition {
            intron: 1,
            offset: 5,
        };
        let (base, offset) = exon_data.ivs_to_cds(&ivs).unwrap();
        assert_eq!(base, 100);
        assert_eq!(offset, Some(5));

        // IVS1-10 should be position 201 with offset -10
        let ivs = IvsPosition {
            intron: 1,
            offset: -10,
        };
        let (base, offset) = exon_data.ivs_to_cds(&ivs).unwrap();
        assert_eq!(base, 201);
        assert_eq!(offset, Some(-10));
    }

    #[test]
    fn test_exon_data_ivs_invalid() {
        let exon_data = ExonData::new(vec![(1, 100)]);

        // IVS0 is invalid
        let ivs = IvsPosition {
            intron: 0,
            offset: 5,
        };
        assert!(exon_data.ivs_to_cds(&ivs).is_err());

        // IVS2 is invalid (only 1 exon)
        let ivs = IvsPosition {
            intron: 2,
            offset: 5,
        };
        assert!(exon_data.ivs_to_cds(&ivs).is_err());
    }

    #[test]
    fn test_legacy_protein_arrow_with_accession() {
        // Test the full flow: NP_000079.2:p.V600>E should parse correctly
        let parser = LegacyParser::new(LegacyConfig::default());
        let result = parser.parse("NP_000079.2:p.V600>E");

        // The protein module converts V600>E to Val600Glu
        assert!(result.is_ok(), "Should parse successfully");
        let result = result.unwrap();
        assert!(result.is_ok(), "Should have a parsed variant");
        assert!(result.is_legacy(), "Should be detected as legacy");
        assert_eq!(
            result.modern_notation,
            Some("NP_000079.2:p.Val600Glu".to_string())
        );
        assert!(result
            .legacy_formats
            .contains(&LegacyFormat::ArrowProteinSubstitution));
    }

    #[test]
    fn test_legacy_protein_number_first_with_accession() {
        // Test the number-first format: NP_000001.1:p.600V>E
        let parser = LegacyParser::new(LegacyConfig::default());
        let result = parser.parse("NP_000001.1:p.600V>E");

        assert!(result.is_ok(), "Should parse successfully");
        let result = result.unwrap();
        assert!(result.is_ok(), "Should have a parsed variant");
        assert!(result.is_legacy(), "Should be detected as legacy");
        assert_eq!(
            result.modern_notation,
            Some("NP_000001.1:p.Val600Glu".to_string())
        );
    }
}
