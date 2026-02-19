//! VCF record representation
//!
//! This module provides a VCF record type for representing variants
//! from VCF (Variant Call Format) files.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

use crate::reference::transcript::GenomeBuild;

/// A single VCF record representing one variant
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VcfRecord {
    /// Chromosome name (e.g., "chr1", "1", "X", "chrM")
    pub chrom: String,

    /// 1-based position of the first base in the reference allele
    pub pos: u64,

    /// Variant identifier (e.g., rsID), None if "."
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,

    /// Reference allele (non-empty string of ACGTN)
    pub reference: String,

    /// Alternate allele(s) - at least one for variant records
    pub alternate: Vec<String>,

    /// Phred-scaled quality score, None if "." or not present
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quality: Option<f32>,

    /// Filter status: None means PASS, Some contains filter name(s)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,

    /// INFO field key-value pairs
    #[serde(default)]
    pub info: HashMap<String, InfoValue>,

    /// FORMAT field specification (e.g., "GT:DP:GQ")
    #[serde(skip_serializing_if = "Option::is_none")]
    pub format: Option<String>,

    /// Sample genotype data, one HashMap per sample
    #[serde(default)]
    pub samples: Vec<HashMap<String, String>>,

    /// Genome build/assembly (for reference context)
    #[serde(default)]
    pub genome_build: GenomeBuild,
}

/// INFO field value types
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum InfoValue {
    /// Flag (presence indicates true)
    Flag,
    /// Integer value
    Integer(i64),
    /// Float value
    Float(f64),
    /// String value
    String(String),
    /// Character value
    Character(char),
    /// Multiple integer values
    IntegerArray(Vec<i64>),
    /// Multiple float values
    FloatArray(Vec<f64>),
    /// Multiple string values
    StringArray(Vec<String>),
}

impl fmt::Display for InfoValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InfoValue::Flag => Ok(()),
            InfoValue::Integer(v) => write!(f, "{}", v),
            InfoValue::Float(v) => write!(f, "{}", v),
            InfoValue::String(v) => write!(f, "{}", v),
            InfoValue::Character(v) => write!(f, "{}", v),
            InfoValue::IntegerArray(v) => {
                let s: Vec<_> = v.iter().map(|x| x.to_string()).collect();
                write!(f, "{}", s.join(","))
            }
            InfoValue::FloatArray(v) => {
                let s: Vec<_> = v.iter().map(|x| x.to_string()).collect();
                write!(f, "{}", s.join(","))
            }
            InfoValue::StringArray(v) => write!(f, "{}", v.join(",")),
        }
    }
}

impl VcfRecord {
    /// Create a new VCF record with minimal required fields
    pub fn new(chrom: String, pos: u64, reference: String, alternate: Vec<String>) -> Self {
        Self {
            chrom,
            pos,
            id: None,
            reference,
            alternate,
            quality: None,
            filter: None,
            info: HashMap::new(),
            format: None,
            samples: Vec::new(),
            genome_build: GenomeBuild::default(),
        }
    }

    /// Create a VCF record for a SNV (single nucleotide variant)
    pub fn snv(chrom: &str, pos: u64, reference: char, alternate: char) -> Self {
        Self::new(
            chrom.to_string(),
            pos,
            reference.to_string(),
            vec![alternate.to_string()],
        )
    }

    /// Create a VCF record for a deletion
    pub fn deletion(chrom: &str, pos: u64, reference: &str) -> Self {
        let anchor = reference.chars().next().unwrap_or('N');
        Self::new(
            chrom.to_string(),
            pos,
            reference.to_string(),
            vec![anchor.to_string()],
        )
    }

    /// Create a VCF record for an insertion
    pub fn insertion(chrom: &str, pos: u64, anchor: char, inserted: &str) -> Self {
        Self::new(
            chrom.to_string(),
            pos,
            anchor.to_string(),
            vec![format!("{}{}", anchor, inserted)],
        )
    }

    /// Get the end position (1-based, inclusive) of the reference allele
    pub fn end_pos(&self) -> u64 {
        self.pos + self.reference.len() as u64 - 1
    }

    /// Check if this is a SNV (single nucleotide variant)
    pub fn is_snv(&self) -> bool {
        self.reference.len() == 1
            && self.alternate.len() == 1
            && self.alternate[0].len() == 1
            && self.reference != self.alternate[0]
    }

    /// Check if this is a simple insertion
    pub fn is_insertion(&self) -> bool {
        self.alternate.len() == 1
            && self.reference.len() == 1
            && self.alternate[0].len() > 1
            && self.alternate[0].starts_with(&self.reference)
    }

    /// Check if this is a simple deletion
    pub fn is_deletion(&self) -> bool {
        self.alternate.len() == 1
            && self.reference.len() > 1
            && self.alternate[0].len() == 1
            && self.reference.starts_with(&self.alternate[0])
    }

    /// Check if this is a complex variant (delins, MNV, etc.)
    pub fn is_complex(&self) -> bool {
        !self.is_snv() && !self.is_insertion() && !self.is_deletion()
    }

    /// Check if this is a multi-allelic variant
    pub fn is_multiallelic(&self) -> bool {
        self.alternate.len() > 1
    }

    /// Get the variant type as a string
    pub fn variant_type(&self) -> &'static str {
        if self.is_snv() {
            "SNV"
        } else if self.is_insertion() {
            "INS"
        } else if self.is_deletion() {
            "DEL"
        } else if self.is_multiallelic() {
            "MULTI"
        } else {
            "COMPLEX"
        }
    }

    /// Set the variant ID (e.g., rsID)
    pub fn with_id(mut self, id: &str) -> Self {
        self.id = Some(id.to_string());
        self
    }

    /// Set the quality score
    pub fn with_quality(mut self, quality: f32) -> Self {
        self.quality = Some(quality);
        self
    }

    /// Set the filter field
    pub fn with_filter(mut self, filter: &str) -> Self {
        self.filter = Some(filter.to_string());
        self
    }

    /// Add an INFO field
    pub fn with_info(mut self, key: &str, value: InfoValue) -> Self {
        self.info.insert(key.to_string(), value);
        self
    }

    /// Set the genome build
    pub fn with_genome_build(mut self, build: GenomeBuild) -> Self {
        self.genome_build = build;
        self
    }

    /// Get an INFO field value as a string
    pub fn get_info_str(&self, key: &str) -> Option<String> {
        self.info.get(key).map(|v| v.to_string())
    }

    /// Check if the variant passes all filters
    pub fn passes_filters(&self) -> bool {
        self.filter.is_none() || self.filter.as_deref() == Some("PASS")
    }

    /// Normalize chromosome name to standard format (add "chr" prefix if missing)
    pub fn normalized_chrom(&self) -> String {
        if self.chrom.starts_with("chr") {
            self.chrom.clone()
        } else {
            format!("chr{}", self.chrom)
        }
    }

    /// Get chromosome name without "chr" prefix
    pub fn bare_chrom(&self) -> &str {
        self.chrom.strip_prefix("chr").unwrap_or(&self.chrom)
    }
}

impl fmt::Display for VcfRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Format as VCF line
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.pos,
            self.id.as_deref().unwrap_or("."),
            self.reference,
            self.alternate.join(","),
            self.quality.map_or(".".to_string(), |q| q.to_string()),
            self.filter.as_deref().unwrap_or("PASS"),
        )?;

        // INFO field
        if self.info.is_empty() {
            write!(f, "\t.")?;
        } else {
            let info_str: Vec<String> = self
                .info
                .iter()
                .map(|(k, v)| {
                    if matches!(v, InfoValue::Flag) {
                        k.clone()
                    } else {
                        format!("{}={}", k, v)
                    }
                })
                .collect();
            write!(f, "\t{}", info_str.join(";"))?;
        }

        // FORMAT and samples if present
        if let Some(format) = &self.format {
            write!(f, "\t{}", format)?;
            for sample in &self.samples {
                let fields: Vec<_> = format
                    .split(':')
                    .map(|key| sample.get(key).map(|s| s.as_str()).unwrap_or("."))
                    .collect();
                write!(f, "\t{}", fields.join(":"))?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_record() {
        let record = VcfRecord::new(
            "chr1".to_string(),
            12345,
            "A".to_string(),
            vec!["G".to_string()],
        );

        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.pos, 12345);
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["G"]);
        assert!(record.id.is_none());
        assert!(record.quality.is_none());
        assert!(record.passes_filters());
    }

    #[test]
    fn test_snv() {
        let record = VcfRecord::snv("chr1", 100, 'A', 'G');
        assert!(record.is_snv());
        assert!(!record.is_insertion());
        assert!(!record.is_deletion());
        assert_eq!(record.variant_type(), "SNV");
    }

    #[test]
    fn test_deletion() {
        let record = VcfRecord::deletion("chr1", 100, "ATG");
        assert!(record.is_deletion());
        assert!(!record.is_snv());
        assert_eq!(record.variant_type(), "DEL");
        assert_eq!(record.reference, "ATG");
        assert_eq!(record.alternate, vec!["A"]);
    }

    #[test]
    fn test_insertion() {
        let record = VcfRecord::insertion("chr1", 100, 'A', "TG");
        assert!(record.is_insertion());
        assert!(!record.is_snv());
        assert_eq!(record.variant_type(), "INS");
        assert_eq!(record.reference, "A");
        assert_eq!(record.alternate, vec!["ATG"]);
    }

    #[test]
    fn test_end_pos() {
        let snv = VcfRecord::snv("chr1", 100, 'A', 'G');
        assert_eq!(snv.end_pos(), 100);

        let del = VcfRecord::deletion("chr1", 100, "ATG");
        assert_eq!(del.end_pos(), 102); // 100, 101, 102
    }

    #[test]
    fn test_multiallelic() {
        let record = VcfRecord::new(
            "chr1".to_string(),
            100,
            "A".to_string(),
            vec!["G".to_string(), "T".to_string()],
        );
        assert!(record.is_multiallelic());
        assert_eq!(record.variant_type(), "MULTI");
    }

    #[test]
    fn test_builder_methods() {
        let record = VcfRecord::snv("chr1", 100, 'A', 'G')
            .with_id("rs12345")
            .with_quality(30.0)
            .with_filter("PASS")
            .with_genome_build(GenomeBuild::GRCh38);

        assert_eq!(record.id, Some("rs12345".to_string()));
        assert_eq!(record.quality, Some(30.0));
        assert_eq!(record.filter, Some("PASS".to_string()));
        assert_eq!(record.genome_build, GenomeBuild::GRCh38);
    }

    #[test]
    fn test_info_value_display() {
        assert_eq!(format!("{}", InfoValue::Integer(42)), "42");
        assert_eq!(format!("{}", InfoValue::Float(1.23)), "1.23");
        assert_eq!(format!("{}", InfoValue::String("test".to_string())), "test");
        assert_eq!(
            format!("{}", InfoValue::IntegerArray(vec![1, 2, 3])),
            "1,2,3"
        );
    }

    #[test]
    fn test_chrom_normalization() {
        let record1 = VcfRecord::snv("1", 100, 'A', 'G');
        assert_eq!(record1.normalized_chrom(), "chr1");
        assert_eq!(record1.bare_chrom(), "1");

        let record2 = VcfRecord::snv("chr1", 100, 'A', 'G');
        assert_eq!(record2.normalized_chrom(), "chr1");
        assert_eq!(record2.bare_chrom(), "1");
    }

    #[test]
    fn test_display() {
        let record = VcfRecord::snv("chr1", 12345, 'A', 'G')
            .with_id("rs123")
            .with_quality(30.0);

        let s = format!("{}", record);
        assert!(s.contains("chr1"));
        assert!(s.contains("12345"));
        assert!(s.contains("rs123"));
        assert!(s.contains("A"));
        assert!(s.contains("G"));
    }

    #[test]
    fn test_passes_filters() {
        let pass1 = VcfRecord::snv("chr1", 100, 'A', 'G');
        assert!(pass1.passes_filters());

        let pass2 = VcfRecord::snv("chr1", 100, 'A', 'G').with_filter("PASS");
        assert!(pass2.passes_filters());

        let fail = VcfRecord::snv("chr1", 100, 'A', 'G').with_filter("LowQual");
        assert!(!fail.passes_filters());
    }
}
