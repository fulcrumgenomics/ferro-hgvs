//! VCF-based rsID lookup.
//!
//! Provides rsID lookup from dbSNP or other VCF files containing variant IDs.
//!
//! # Example
//!
//! ```no_run
//! use ferro_hgvs::rsid::vcf::VcfRsIdLookup;
//! use ferro_hgvs::rsid::RsIdLookup;  // Import trait to use lookup method
//! use ferro_hgvs::reference::transcript::GenomeBuild;
//!
//! // Load from a dbSNP VCF file
//! let lookup = VcfRsIdLookup::from_vcf("dbsnp.vcf.gz", GenomeBuild::GRCh38).unwrap();
//!
//! // Look up an rsID
//! let results = lookup.lookup("rs121913529").unwrap();
//! for result in results {
//!     println!("{}: {}:{} {}>{}",
//!         result.rsid, result.contig, result.position,
//!         result.reference, result.alternate
//!     );
//! }
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::FerroError;
use crate::reference::transcript::GenomeBuild;

use super::{parse_rsid, RsIdLookup, RsIdMetadata, RsIdResult};

/// VCF-based rsID lookup.
///
/// Loads variant data from VCF files (like dbSNP) and provides rsID lookup.
/// Supports both plain text and gzip-compressed VCF files.
#[derive(Debug, Clone)]
pub struct VcfRsIdLookup {
    /// Indexed entries: rsID numeric value -> list of results
    entries: HashMap<u64, Vec<RsIdResult>>,
    /// Metadata about the data source
    metadata: RsIdMetadata,
}

impl VcfRsIdLookup {
    /// Create a new empty VCF lookup.
    pub fn new(build: GenomeBuild) -> Self {
        Self {
            entries: HashMap::new(),
            metadata: RsIdMetadata {
                source: "VCF".to_string(),
                build,
                version: "unknown".to_string(),
                variant_count: 0,
            },
        }
    }

    /// Load rsID lookup from a VCF file.
    ///
    /// Automatically detects gzip compression based on `.gz` extension.
    pub fn from_vcf<P: AsRef<Path>>(path: P, build: GenomeBuild) -> Result<Self, FerroError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open VCF file '{}': {}", path.display(), e),
        })?;

        let filename = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        if path.extension().is_some_and(|ext| ext == "gz") {
            let decoder = GzDecoder::new(file);
            let reader = BufReader::new(decoder);
            Self::from_reader(reader, build, filename)
        } else {
            let reader = BufReader::new(file);
            Self::from_reader(reader, build, filename)
        }
    }

    /// Load rsID lookup from a buffered reader.
    pub fn from_reader<R: BufRead>(
        reader: R,
        build: GenomeBuild,
        source_name: &str,
    ) -> Result<Self, FerroError> {
        let mut lookup = Self::new(build);
        lookup.metadata.source = format!("VCF: {}", source_name);

        let mut version = None;

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read VCF line: {}", e),
            })?;

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            // Handle header lines
            if line.starts_with('#') {
                // Extract version from ##fileformat line
                if let Some(ver) = line.strip_prefix("##fileformat=") {
                    version = Some(ver.to_string());
                }
                // Extract dbSNP version if present
                if let Some(build_id) = line.strip_prefix("##dbSNP_BUILD_ID=") {
                    version = Some(format!("dbSNP {}", build_id));
                }
                continue;
            }

            // Parse data line
            if let Some(results) = parse_vcf_line(&line, build) {
                for result in results {
                    if let Ok(rsid_num) = parse_rsid(&result.rsid) {
                        lookup.entries.entry(rsid_num).or_default().push(result);
                    }
                }
            }
        }

        lookup.metadata.version = version.unwrap_or_else(|| "unknown".to_string());
        lookup.metadata.variant_count = lookup.entries.values().map(|v| v.len() as u64).sum();

        Ok(lookup)
    }

    /// Add a single rsID entry.
    pub fn add(&mut self, result: RsIdResult) -> Result<(), FerroError> {
        let rsid_num = parse_rsid(&result.rsid)?;
        self.entries.entry(rsid_num).or_default().push(result);
        self.metadata.variant_count = self.entries.values().map(|v| v.len() as u64).sum();
        Ok(())
    }

    /// Get the number of unique rsIDs indexed.
    pub fn rsid_count(&self) -> usize {
        self.entries.len()
    }

    /// Get the total number of variant entries (may be more than rsID count
    /// for multi-allelic variants).
    pub fn variant_count(&self) -> u64 {
        self.metadata.variant_count
    }
}

impl RsIdLookup for VcfRsIdLookup {
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

/// Parse a VCF data line and extract rsID results.
fn parse_vcf_line(line: &str, build: GenomeBuild) -> Option<Vec<RsIdResult>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return None;
    }

    let chrom = fields[0];
    let pos: u64 = fields[1].parse().ok()?;
    let id_field = fields[2];
    let reference = fields[3];
    let alt_field = fields[4];

    // Skip if no ID or ID is "."
    if id_field.is_empty() || id_field == "." {
        return None;
    }

    // VCF field indices: CHROM(0), POS(1), ID(2), REF(3), ALT(4), QUAL(5), FILTER(6), INFO(7)
    // Parse INFO field for allele frequency if present
    let allele_frequency = if fields.len() > 7 {
        parse_allele_frequency(fields[7]) // INFO field
    } else {
        None
    };

    // Parse clinical significance from INFO if present
    let clinical_significance = if fields.len() > 7 {
        parse_clinical_significance(fields[7]) // INFO field
    } else {
        None
    };

    // Handle multiple IDs (semicolon-separated)
    let ids: Vec<&str> = id_field.split(';').collect();

    // Handle multiple alternate alleles (comma-separated)
    let alts: Vec<&str> = alt_field.split(',').collect();

    let mut results = Vec::new();

    // Create result for each rsID + alt allele combination
    for id in &ids {
        // Only process rsIDs
        if !id.starts_with("rs") {
            continue;
        }

        for alt in &alts {
            // Skip symbolic alleles like <DEL>, <INS>, etc.
            if alt.starts_with('<') && alt.ends_with('>') {
                continue;
            }

            let result = RsIdResult {
                rsid: (*id).to_string(),
                contig: chrom.to_string(),
                position: pos,
                reference: reference.to_string(),
                alternate: (*alt).to_string(),
                build,
                hgvs: None,
                allele_frequency,
                clinical_significance: clinical_significance.clone(),
            };
            results.push(result);
        }
    }

    if results.is_empty() {
        None
    } else {
        Some(results)
    }
}

/// Parse allele frequency from INFO field.
fn parse_allele_frequency(info: &str) -> Option<f64> {
    // Look for AF= or CAF= field
    for field in info.split(';') {
        if let Some(af_str) = field.strip_prefix("AF=") {
            // Take first value if comma-separated
            let first = af_str.split(',').next()?;
            return first.parse().ok();
        }
        if let Some(caf_str) = field.strip_prefix("CAF=") {
            // CAF is ref,alt1,alt2,... - take second value (first alt)
            let mut parts = caf_str.split(',');
            parts.next(); // Skip ref frequency
            let alt_freq = parts.next()?;
            if alt_freq != "." {
                return alt_freq.parse().ok();
            }
        }
    }
    None
}

/// Parse clinical significance from INFO field.
fn parse_clinical_significance(info: &str) -> Option<String> {
    for field in info.split(';') {
        // dbSNP uses CLNSIG
        if let Some(clnsig) = field.strip_prefix("CLNSIG=") {
            return Some(clnsig.to_string());
        }
        // ClinVar uses CLNSIG as well
        if let Some(clinical) = field.strip_prefix("CLINICAL_SIGNIFICANCE=") {
            return Some(clinical.to_string());
        }
    }
    None
}

/// Streaming VCF rsID lookup for large files.
///
/// Instead of loading the entire VCF into memory, this searches the file
/// on each lookup. Suitable for very large files where memory is constrained.
#[derive(Debug)]
pub struct StreamingVcfRsIdLookup {
    path: std::path::PathBuf,
    metadata: RsIdMetadata,
}

impl StreamingVcfRsIdLookup {
    /// Create a new streaming lookup from a VCF file.
    pub fn new<P: AsRef<Path>>(path: P, build: GenomeBuild) -> Result<Self, FerroError> {
        let path = path.as_ref();

        // Verify file exists
        if !path.exists() {
            return Err(FerroError::Io {
                msg: format!("VCF file not found: {}", path.display()),
            });
        }

        let filename = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        Ok(Self {
            path: path.to_path_buf(),
            metadata: RsIdMetadata {
                source: format!("VCF (streaming): {}", filename),
                build,
                version: "unknown".to_string(),
                variant_count: 0, // Unknown for streaming
            },
        })
    }

    /// Search for an rsID in the VCF file.
    fn search_file(&self, rsid_num: u64) -> Result<Vec<RsIdResult>, FerroError> {
        let file = File::open(&self.path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open VCF file: {}", e),
        })?;

        let reader: Box<dyn BufRead> = if self.path.extension().is_some_and(|ext| ext == "gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let rsid_str = format!("rs{}", rsid_num);
        let mut results = Vec::new();

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read VCF line: {}", e),
            })?;

            // Skip headers
            if line.starts_with('#') {
                continue;
            }

            // Quick check if this line might contain our rsID
            if !line.contains(&rsid_str) {
                continue;
            }

            // Parse the line
            if let Some(line_results) = parse_vcf_line(&line, self.metadata.build) {
                for result in line_results {
                    if result.rsid == rsid_str {
                        results.push(result);
                    }
                }
            }
        }

        if results.is_empty() {
            Err(FerroError::ReferenceNotFound { id: rsid_str })
        } else {
            Ok(results)
        }
    }
}

impl RsIdLookup for StreamingVcfRsIdLookup {
    fn lookup(&self, rsid: &str) -> Result<Vec<RsIdResult>, FerroError> {
        let rsid_num = parse_rsid(rsid)?;
        self.search_file(rsid_num)
    }

    fn contains(&self, rsid: &str) -> bool {
        self.lookup(rsid).is_ok()
    }

    fn metadata(&self) -> &RsIdMetadata {
        &self.metadata
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_VCF: &str = r#"##fileformat=VCFv4.3
##dbSNP_BUILD_ID=156
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	140753336	rs121913529	A	T	.	.	AF=0.0001;CLNSIG=pathogenic
chr17	43124027	rs80357906	GAG	G	.	.	AF=0.0001
chr1	12345	rs12345;rs67890	A	G,T	.	.	AF=0.01,0.02
chr1	100	.	A	G	.	.	.
chr1	200	rs999	C	<DEL>	.	.	.
"#;

    #[test]
    fn test_from_reader() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        assert_eq!(lookup.metadata.version, "dbSNP 156");
        assert!(lookup.rsid_count() >= 3); // rs121913529, rs80357906, rs12345, rs67890
    }

    #[test]
    fn test_lookup_snv() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        let results = lookup.lookup("rs121913529").unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].contig, "chr7");
        assert_eq!(results[0].position, 140753336);
        assert_eq!(results[0].reference, "A");
        assert_eq!(results[0].alternate, "T");
        assert!(results[0].is_snv());
        assert_eq!(results[0].allele_frequency, Some(0.0001));
        assert_eq!(
            results[0].clinical_significance,
            Some("pathogenic".to_string())
        );
    }

    #[test]
    fn test_lookup_deletion() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        let results = lookup.lookup("rs80357906").unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].contig, "chr17");
        assert!(results[0].is_deletion());
    }

    #[test]
    fn test_lookup_multiallelic() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        // rs12345 should have two entries (one for each alt allele)
        let results = lookup.lookup("rs12345").unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].alternate, "G");
        assert_eq!(results[1].alternate, "T");
    }

    #[test]
    fn test_lookup_multiple_ids() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        // Both rs12345 and rs67890 should be found
        assert!(lookup.contains("rs12345"));
        assert!(lookup.contains("rs67890"));
    }

    #[test]
    fn test_lookup_not_found() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        let result = lookup.lookup("rs99999999");
        assert!(result.is_err());
    }

    #[test]
    fn test_contains() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        assert!(lookup.contains("rs121913529"));
        assert!(lookup.contains("121913529")); // Without prefix
        assert!(!lookup.contains("rs99999999"));
    }

    #[test]
    fn test_skips_missing_id() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        // The line with "." as ID should not create any entries
        // So we just verify other lookups work
        assert!(lookup.rsid_count() >= 3);
    }

    #[test]
    fn test_skips_symbolic_alleles() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        // rs999 has <DEL> as alt, should not have any results
        assert!(!lookup.contains("rs999"));
    }

    #[test]
    fn test_to_hgvs() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        let results = lookup.lookup("rs121913529").unwrap();
        assert_eq!(results[0].to_hgvs(), "g.140753336A>T");
    }

    #[test]
    fn test_add_entry() {
        let mut lookup = VcfRsIdLookup::new(GenomeBuild::GRCh38);

        let result = RsIdResult::new(
            "rs12345".to_string(),
            "chr1".to_string(),
            100,
            "A".to_string(),
            "G".to_string(),
            GenomeBuild::GRCh38,
        );

        lookup.add(result).unwrap();
        assert!(lookup.contains("rs12345"));
        assert_eq!(lookup.rsid_count(), 1);
    }

    #[test]
    fn test_parse_allele_frequency() {
        assert_eq!(parse_allele_frequency("AF=0.01"), Some(0.01));
        assert_eq!(parse_allele_frequency("DP=100;AF=0.05"), Some(0.05));
        assert_eq!(parse_allele_frequency("AF=0.1,0.2,0.3"), Some(0.1));
        assert_eq!(parse_allele_frequency("CAF=0.99,0.01"), Some(0.01));
        assert_eq!(parse_allele_frequency("DP=100"), None);
    }

    #[test]
    fn test_parse_clinical_significance() {
        assert_eq!(
            parse_clinical_significance("CLNSIG=pathogenic"),
            Some("pathogenic".to_string())
        );
        assert_eq!(
            parse_clinical_significance("DP=100;CLNSIG=benign"),
            Some("benign".to_string())
        );
        assert_eq!(parse_clinical_significance("DP=100"), None);
    }

    #[test]
    fn test_metadata() {
        let reader = BufReader::new(TEST_VCF.as_bytes());
        let lookup = VcfRsIdLookup::from_reader(reader, GenomeBuild::GRCh38, "test.vcf").unwrap();

        let meta = lookup.metadata();
        assert!(meta.source.contains("test.vcf"));
        assert_eq!(meta.build, GenomeBuild::GRCh38);
        assert!(meta.variant_count > 0);
    }
}
