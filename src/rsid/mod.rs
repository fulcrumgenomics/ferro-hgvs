//! rsID/dbSNP lookup functionality.
//!
//! This module provides traits and implementations for looking up
//! rsIDs from dbSNP and converting them to HGVS notation.
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::rsid::{parse_rsid, RsIdResult};
//!
//! // Parse rsID string to numeric value
//! let rsid_num = parse_rsid("rs121913529").unwrap();
//! assert_eq!(rsid_num, 121913529);
//! ```

pub mod vcf;

use crate::error::FerroError;
use crate::reference::transcript::GenomeBuild;

/// Result of rsID lookup.
#[derive(Debug, Clone)]
pub struct RsIdResult {
    /// rsID (e.g., "rs121913529").
    pub rsid: String,
    /// Chromosome/contig (e.g., "chr7" or "7").
    pub contig: String,
    /// 1-based position.
    pub position: u64,
    /// Reference allele.
    pub reference: String,
    /// Alternate allele.
    pub alternate: String,
    /// Genome build.
    pub build: GenomeBuild,
    /// HGVS notation if available.
    pub hgvs: Option<String>,
    /// Allele frequency if available.
    pub allele_frequency: Option<f64>,
    /// Clinical significance if available.
    pub clinical_significance: Option<String>,
}

impl RsIdResult {
    /// Create a new rsID result.
    pub fn new(
        rsid: String,
        contig: String,
        position: u64,
        reference: String,
        alternate: String,
        build: GenomeBuild,
    ) -> Self {
        Self {
            rsid,
            contig,
            position,
            reference,
            alternate,
            build,
            hgvs: None,
            allele_frequency: None,
            clinical_significance: None,
        }
    }

    /// Check if this is a substitution (SNV).
    pub fn is_snv(&self) -> bool {
        self.reference.len() == 1 && self.alternate.len() == 1
    }

    /// Check if this is a deletion.
    pub fn is_deletion(&self) -> bool {
        self.reference.len() > self.alternate.len()
    }

    /// Check if this is an insertion.
    pub fn is_insertion(&self) -> bool {
        self.reference.len() < self.alternate.len()
    }

    /// Generate simple genomic HGVS notation.
    pub fn to_hgvs(&self) -> String {
        let prefix = format!("g.{}", self.position);

        if self.is_snv() {
            format!("{}{}>{}", prefix, self.reference, self.alternate)
        } else if self.is_deletion() {
            let del_len = self.reference.len() - self.alternate.len();
            if del_len == 1 {
                format!("{}del", self.position + self.alternate.len() as u64)
            } else {
                let start = self.position + self.alternate.len() as u64;
                let end = start + del_len as u64 - 1;
                format!("g.{}_{}del", start, end)
            }
        } else if self.is_insertion() {
            let ins_seq = &self.alternate[self.reference.len()..];
            format!(
                "g.{}_{}ins{}",
                self.position + self.reference.len() as u64 - 1,
                self.position + self.reference.len() as u64,
                ins_seq
            )
        } else {
            // Complex delins
            format!("g.{}delins{}", prefix, self.alternate)
        }
    }
}

/// Metadata about the rsID lookup source.
#[derive(Debug, Clone)]
pub struct RsIdMetadata {
    /// Source name (e.g., "dbSNP VCF").
    pub source: String,
    /// Genome build.
    pub build: GenomeBuild,
    /// Version string.
    pub version: String,
    /// Number of variants indexed.
    pub variant_count: u64,
}

/// Trait for rsID lookup implementations.
pub trait RsIdLookup {
    /// Look up rsID and return matching variants.
    fn lookup(&self, rsid: &str) -> Result<Vec<RsIdResult>, FerroError>;

    /// Check if rsID exists.
    fn contains(&self, rsid: &str) -> bool;

    /// Get metadata about the data source.
    fn metadata(&self) -> &RsIdMetadata;
}

/// Parse rsID string to numeric value.
///
/// Accepts formats: "rs121913529" or "121913529".
pub fn parse_rsid(rsid: &str) -> Result<u64, FerroError> {
    let s = rsid.strip_prefix("rs").unwrap_or(rsid);
    s.parse().map_err(|_| FerroError::InvalidCoordinates {
        msg: format!("Invalid rsID: {}", rsid),
    })
}

/// Format numeric rsID to string with "rs" prefix.
pub fn format_rsid(rsid_num: u64) -> String {
    format!("rs{}", rsid_num)
}

/// Simple in-memory rsID lookup for testing.
#[derive(Debug, Clone)]
pub struct InMemoryRsIdLookup {
    entries: std::collections::HashMap<u64, Vec<RsIdResult>>,
    metadata: RsIdMetadata,
}

impl InMemoryRsIdLookup {
    /// Create a new in-memory lookup.
    pub fn new(build: GenomeBuild) -> Self {
        Self {
            entries: std::collections::HashMap::new(),
            metadata: RsIdMetadata {
                source: "In-memory".to_string(),
                build,
                version: "1.0".to_string(),
                variant_count: 0,
            },
        }
    }

    /// Add an rsID entry.
    pub fn add(&mut self, result: RsIdResult) -> Result<(), FerroError> {
        let rsid_num = parse_rsid(&result.rsid)?;
        self.entries.entry(rsid_num).or_default().push(result);
        self.metadata.variant_count = self.entries.values().map(|v| v.len() as u64).sum();
        Ok(())
    }

    /// Create with common test variants.
    pub fn with_test_data() -> Self {
        let mut lookup = Self::new(GenomeBuild::GRCh38);

        // BRAF V600E - rs121913529
        let _ = lookup.add(RsIdResult {
            rsid: "rs121913529".to_string(),
            contig: "chr7".to_string(),
            position: 140753336,
            reference: "A".to_string(),
            alternate: "T".to_string(),
            build: GenomeBuild::GRCh38,
            hgvs: Some("NC_000007.14:g.140753336A>T".to_string()),
            allele_frequency: Some(0.0001),
            clinical_significance: Some("pathogenic".to_string()),
        });

        // BRCA1 185delAG - rs80357906
        let _ = lookup.add(RsIdResult {
            rsid: "rs80357906".to_string(),
            contig: "chr17".to_string(),
            position: 43124027,
            reference: "GAG".to_string(),
            alternate: "G".to_string(),
            build: GenomeBuild::GRCh38,
            hgvs: Some("NC_000017.11:g.43124028_43124029del".to_string()),
            allele_frequency: Some(0.0001),
            clinical_significance: Some("pathogenic".to_string()),
        });

        lookup
    }
}

impl RsIdLookup for InMemoryRsIdLookup {
    fn lookup(&self, rsid: &str) -> Result<Vec<RsIdResult>, FerroError> {
        let rsid_num = parse_rsid(rsid)?;
        self.entries
            .get(&rsid_num)
            .cloned()
            .ok_or_else(|| FerroError::ReferenceNotFound {
                id: rsid.to_string(),
            })
    }

    fn contains(&self, rsid: &str) -> bool {
        parse_rsid(rsid)
            .map(|n| self.entries.contains_key(&n))
            .unwrap_or(false)
    }

    fn metadata(&self) -> &RsIdMetadata {
        &self.metadata
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_rsid() {
        assert_eq!(parse_rsid("rs121913529").unwrap(), 121913529);
        assert_eq!(parse_rsid("121913529").unwrap(), 121913529);
        assert!(parse_rsid("rsABC").is_err());
    }

    #[test]
    fn test_format_rsid() {
        assert_eq!(format_rsid(121913529), "rs121913529");
    }

    #[test]
    fn test_rsid_result_is_snv() {
        let snv = RsIdResult::new(
            "rs1".to_string(),
            "chr1".to_string(),
            100,
            "A".to_string(),
            "G".to_string(),
            GenomeBuild::GRCh38,
        );
        assert!(snv.is_snv());
        assert!(!snv.is_deletion());
        assert!(!snv.is_insertion());
    }

    #[test]
    fn test_rsid_result_is_deletion() {
        let del = RsIdResult::new(
            "rs1".to_string(),
            "chr1".to_string(),
            100,
            "AG".to_string(),
            "A".to_string(),
            GenomeBuild::GRCh38,
        );
        assert!(!del.is_snv());
        assert!(del.is_deletion());
        assert!(!del.is_insertion());
    }

    #[test]
    fn test_rsid_result_is_insertion() {
        let ins = RsIdResult::new(
            "rs1".to_string(),
            "chr1".to_string(),
            100,
            "A".to_string(),
            "AG".to_string(),
            GenomeBuild::GRCh38,
        );
        assert!(!ins.is_snv());
        assert!(!ins.is_deletion());
        assert!(ins.is_insertion());
    }

    #[test]
    fn test_in_memory_lookup() {
        let lookup = InMemoryRsIdLookup::with_test_data();

        // Test BRAF V600E
        let results = lookup.lookup("rs121913529").unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].contig, "chr7");
        assert!(results[0].is_snv());

        // Test BRCA1
        let results = lookup.lookup("rs80357906").unwrap();
        assert_eq!(results.len(), 1);
        assert!(results[0].is_deletion());

        // Test not found
        assert!(lookup.lookup("rs999999999").is_err());
    }

    #[test]
    fn test_contains() {
        let lookup = InMemoryRsIdLookup::with_test_data();

        assert!(lookup.contains("rs121913529"));
        assert!(lookup.contains("121913529")); // Without prefix
        assert!(!lookup.contains("rs999999999"));
    }

    #[test]
    fn test_to_hgvs_snv() {
        let snv = RsIdResult::new(
            "rs1".to_string(),
            "chr1".to_string(),
            100,
            "A".to_string(),
            "G".to_string(),
            GenomeBuild::GRCh38,
        );
        assert_eq!(snv.to_hgvs(), "g.100A>G");
    }

    #[test]
    fn test_to_hgvs_insertion() {
        let ins = RsIdResult::new(
            "rs1".to_string(),
            "chr1".to_string(),
            100,
            "A".to_string(),
            "AG".to_string(),
            GenomeBuild::GRCh38,
        );
        assert_eq!(ins.to_hgvs(), "g.100_101insG");
    }

    #[test]
    fn test_metadata() {
        let lookup = InMemoryRsIdLookup::with_test_data();
        let meta = lookup.metadata();

        assert_eq!(meta.source, "In-memory");
        assert_eq!(meta.build, GenomeBuild::GRCh38);
        assert_eq!(meta.variant_count, 2);
    }
}
