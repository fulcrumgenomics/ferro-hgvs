//! Description extractor for generating HGVS from sequences.
//!
//! This module provides functionality to compare a reference sequence with an
//! observed sequence and generate HGVS variant descriptions.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::extractor::{DescriptionExtractor, ExtractorConfig};
//!
//! let extractor = DescriptionExtractor::new(ExtractorConfig::default());
//! let result = extractor.extract("ATGCATGC", "ATGAATGC").unwrap();
//!
//! assert_eq!(result.variants.len(), 1);
//! assert_eq!(result.hgvs_strings[0], "4C>A");
//! ```

mod align;
mod classify;
mod generate;

pub use align::{align_sequences, EditOp, RawEdit};
pub use classify::{classify_edit, is_duplication, is_inversion, EditClassification};
pub use generate::{generate_hgvs, ExtractedVariant, ExtractionResult};

use crate::error::FerroError;

/// Configuration for the description extractor.
#[derive(Debug, Clone)]
pub struct ExtractorConfig {
    /// Maximum sequence length for alignment (default: 100000).
    pub max_length: usize,
    /// Whether to detect duplications (default: true).
    pub detect_duplications: bool,
    /// Whether to detect inversions (default: true).
    pub detect_inversions: bool,
    /// Whether to detect repeat units (default: false for simplicity).
    pub detect_repeats: bool,
    /// Minimum repeat unit length (default: 1).
    pub min_repeat_unit: usize,
}

impl Default for ExtractorConfig {
    fn default() -> Self {
        Self {
            max_length: 100_000,
            detect_duplications: true,
            detect_inversions: true,
            detect_repeats: false,
            min_repeat_unit: 1,
        }
    }
}

/// Description extractor for generating HGVS from sequence comparison.
#[derive(Debug, Clone)]
pub struct DescriptionExtractor {
    config: ExtractorConfig,
}

impl DescriptionExtractor {
    /// Create a new description extractor with the given configuration.
    pub fn new(config: ExtractorConfig) -> Self {
        Self { config }
    }

    /// Create a description extractor with default configuration.
    pub fn with_defaults() -> Self {
        Self::new(ExtractorConfig::default())
    }

    /// Extract HGVS descriptions from reference and observed sequences.
    ///
    /// # Arguments
    ///
    /// * `reference` - The reference sequence
    /// * `observed` - The observed (variant) sequence
    ///
    /// # Returns
    ///
    /// An `ExtractionResult` containing all detected variants and their HGVS strings.
    pub fn extract(&self, reference: &str, observed: &str) -> Result<ExtractionResult, FerroError> {
        self.extract_internal(reference, observed, None)
    }

    /// Extract HGVS descriptions with an accession prefix.
    ///
    /// # Arguments
    ///
    /// * `accession` - The accession to use (e.g., "NC_000001.11")
    /// * `reference` - The reference sequence
    /// * `observed` - The observed (variant) sequence
    pub fn extract_with_accession(
        &self,
        accession: &str,
        reference: &str,
        observed: &str,
    ) -> Result<ExtractionResult, FerroError> {
        self.extract_internal(reference, observed, Some(accession))
    }

    fn extract_internal(
        &self,
        reference: &str,
        observed: &str,
        accession: Option<&str>,
    ) -> Result<ExtractionResult, FerroError> {
        // Validate input
        if reference.is_empty() {
            return Err(FerroError::InvalidCoordinates {
                msg: "Reference sequence is empty".to_string(),
            });
        }

        if reference.len() > self.config.max_length || observed.len() > self.config.max_length {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Sequence too long: ref={}, obs={}, max={}",
                    reference.len(),
                    observed.len(),
                    self.config.max_length
                ),
            });
        }

        // Align sequences to find edits
        let raw_edits = align_sequences(reference, observed);

        // Classify and generate HGVS for each edit
        let mut variants = Vec::new();
        let mut hgvs_strings = Vec::new();

        for edit in raw_edits {
            let classification = classify_edit(
                &edit,
                reference,
                self.config.detect_duplications,
                self.config.detect_inversions,
            );

            let hgvs = generate_hgvs(&edit, &classification, accession);

            variants.push(ExtractedVariant {
                position: edit.ref_start,
                ref_seq: edit.ref_seq.clone(),
                obs_seq: edit.obs_seq.clone(),
                classification,
                hgvs: hgvs.clone(),
            });

            hgvs_strings.push(hgvs);
        }

        Ok(ExtractionResult {
            reference_length: reference.len() as u64,
            observed_length: observed.len() as u64,
            variants,
            hgvs_strings,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_substitution() {
        let extractor = DescriptionExtractor::with_defaults();
        let result = extractor.extract("ATGC", "ATAC").unwrap();

        assert_eq!(result.variants.len(), 1);
        assert_eq!(result.hgvs_strings[0], "3G>A");
    }

    #[test]
    fn test_deletion() {
        let extractor = DescriptionExtractor::with_defaults();
        let result = extractor.extract("ATGCAT", "ATCAT").unwrap();

        assert_eq!(result.variants.len(), 1);
        assert!(result.hgvs_strings[0].contains("del"));
    }

    #[test]
    fn test_insertion() {
        let extractor = DescriptionExtractor::with_defaults();
        let result = extractor.extract("ATGC", "ATGGC").unwrap();

        assert_eq!(result.variants.len(), 1);
        assert!(result.hgvs_strings[0].contains("ins"));
    }

    #[test]
    fn test_with_accession() {
        let extractor = DescriptionExtractor::with_defaults();
        let result = extractor
            .extract_with_accession("NC_000001.11", "ATGC", "ATAC")
            .unwrap();

        assert!(result.hgvs_strings[0].starts_with("NC_000001.11:g."));
    }

    #[test]
    fn test_no_change() {
        let extractor = DescriptionExtractor::with_defaults();
        let result = extractor.extract("ATGC", "ATGC").unwrap();

        assert!(result.variants.is_empty());
    }

    #[test]
    fn test_empty_reference_error() {
        let extractor = DescriptionExtractor::with_defaults();
        let result = extractor.extract("", "ATGC");

        assert!(result.is_err());
    }
}
