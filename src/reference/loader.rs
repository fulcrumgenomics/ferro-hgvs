//! Transcript database loading
//!
//! This module provides functionality to load transcript annotations from
//! standard file formats (GFF3, GTF) commonly used for genome annotations.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::FerroError;
use crate::reference::transcript::{GenomeBuild, ManeStatus, Transcript};

/// Statistics about MANE transcript coverage in the database
#[derive(Debug, Clone, Default)]
pub struct ManeStats {
    /// Total number of unique genes
    pub total_genes: usize,
    /// Number of genes with MANE Select transcript
    pub genes_with_mane_select: usize,
    /// Number of genes with MANE Plus Clinical transcripts
    pub genes_with_mane_plus_clinical: usize,
    /// Total number of MANE Plus Clinical transcripts
    pub total_mane_plus_clinical_transcripts: usize,
}

impl ManeStats {
    /// Get MANE Select coverage percentage
    pub fn mane_select_coverage(&self) -> f64 {
        if self.total_genes == 0 {
            0.0
        } else {
            (self.genes_with_mane_select as f64 / self.total_genes as f64) * 100.0
        }
    }
}

/// A database of transcripts indexed for efficient lookup
#[derive(Debug, Default)]
pub struct TranscriptDb {
    /// Transcripts indexed by accession ID
    transcripts: HashMap<String, Transcript>,
    /// Index from gene symbol to transcript IDs
    gene_index: HashMap<String, Vec<String>>,
    /// Index from chromosome:region to transcript IDs for overlap queries
    region_index: HashMap<String, Vec<(u64, u64, String)>>,
    /// Index from gene symbol to MANE Select transcript ID
    mane_select_index: HashMap<String, String>,
    /// Index from gene symbol to MANE Plus Clinical transcript IDs
    mane_plus_clinical_index: HashMap<String, Vec<String>>,
    /// Genome build
    pub genome_build: GenomeBuild,
}

impl TranscriptDb {
    /// Create a new empty transcript database
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new database with a specified genome build
    pub fn with_build(build: GenomeBuild) -> Self {
        Self {
            genome_build: build,
            ..Default::default()
        }
    }

    /// Add a transcript to the database
    pub fn add(&mut self, transcript: Transcript) {
        let id = transcript.id.clone();

        // Index by gene symbol
        if let Some(ref gene) = transcript.gene_symbol {
            self.gene_index
                .entry(gene.clone())
                .or_default()
                .push(id.clone());

            // Index MANE transcripts
            match transcript.mane_status {
                ManeStatus::Select => {
                    self.mane_select_index.insert(gene.clone(), id.clone());
                }
                ManeStatus::PlusClinical => {
                    self.mane_plus_clinical_index
                        .entry(gene.clone())
                        .or_default()
                        .push(id.clone());
                }
                ManeStatus::None => {}
            }
        }

        // Index by genomic region
        if let (Some(ref chrom), Some(start), Some(end)) = (
            &transcript.chromosome,
            transcript.genomic_start,
            transcript.genomic_end,
        ) {
            self.region_index
                .entry(chrom.clone())
                .or_default()
                .push((start, end, id.clone()));
        }

        self.transcripts.insert(id, transcript);
    }

    /// Get a transcript by its accession ID
    pub fn get(&self, id: &str) -> Option<&Transcript> {
        self.transcripts.get(id)
    }

    /// Get all transcripts for a gene
    pub fn get_by_gene(&self, gene: &str) -> Vec<&Transcript> {
        self.gene_index
            .get(gene)
            .map(|ids| {
                ids.iter()
                    .filter_map(|id| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all transcripts overlapping a genomic position
    pub fn get_by_position(&self, chrom: &str, pos: u64) -> Vec<&Transcript> {
        self.region_index
            .get(chrom)
            .map(|regions| {
                regions
                    .iter()
                    .filter(|(start, end, _)| pos >= *start && pos <= *end)
                    .filter_map(|(_, _, id)| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all transcripts overlapping a genomic region
    pub fn get_by_region(&self, chrom: &str, start: u64, end: u64) -> Vec<&Transcript> {
        self.region_index
            .get(chrom)
            .map(|regions| {
                regions
                    .iter()
                    .filter(|(tx_start, tx_end, _)| {
                        // Check for overlap
                        *tx_start <= end && *tx_end >= start
                    })
                    .filter_map(|(_, _, id)| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get the MANE Select transcript for a gene
    ///
    /// Returns the single representative transcript designated by NCBI/EBI
    /// for this gene, if available.
    pub fn get_mane_select(&self, gene: &str) -> Option<&Transcript> {
        self.mane_select_index
            .get(gene)
            .and_then(|id| self.transcripts.get(id))
    }

    /// Get all MANE Plus Clinical transcripts for a gene
    ///
    /// Returns additional clinically relevant transcripts beyond MANE Select.
    pub fn get_mane_plus_clinical(&self, gene: &str) -> Vec<&Transcript> {
        self.mane_plus_clinical_index
            .get(gene)
            .map(|ids| {
                ids.iter()
                    .filter_map(|id| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all MANE transcripts for a gene (Select + Plus Clinical)
    pub fn get_mane_transcripts(&self, gene: &str) -> Vec<&Transcript> {
        let mut result = Vec::new();

        if let Some(select) = self.get_mane_select(gene) {
            result.push(select);
        }

        result.extend(self.get_mane_plus_clinical(gene));
        result
    }

    /// Get the best transcript for a gene, preferring MANE Select
    ///
    /// Priority order:
    /// 1. MANE Select
    /// 2. MANE Plus Clinical (first one)
    /// 3. Any coding transcript (first one)
    /// 4. Any transcript (first one)
    pub fn get_preferred_transcript(&self, gene: &str) -> Option<&Transcript> {
        // Try MANE Select first
        if let Some(tx) = self.get_mane_select(gene) {
            return Some(tx);
        }

        // Try MANE Plus Clinical
        let plus_clinical = self.get_mane_plus_clinical(gene);
        if !plus_clinical.is_empty() {
            return Some(plus_clinical[0]);
        }

        // Fall back to gene transcripts
        let transcripts = self.get_by_gene(gene);

        // Prefer coding transcripts
        if let Some(coding) = transcripts.iter().find(|tx| tx.is_coding()) {
            return Some(coding);
        }

        // Any transcript
        transcripts.first().copied()
    }

    /// Check if a gene has a MANE Select transcript
    pub fn has_mane_select(&self, gene: &str) -> bool {
        self.mane_select_index.contains_key(gene)
    }

    /// Get statistics about MANE coverage
    pub fn mane_stats(&self) -> ManeStats {
        let total_genes = self.gene_index.len();
        let genes_with_mane_select = self.mane_select_index.len();
        let genes_with_mane_plus = self.mane_plus_clinical_index.len();
        let total_mane_plus_transcripts: usize = self
            .mane_plus_clinical_index
            .values()
            .map(|v| v.len())
            .sum();

        ManeStats {
            total_genes,
            genes_with_mane_select,
            genes_with_mane_plus_clinical: genes_with_mane_plus,
            total_mane_plus_clinical_transcripts: total_mane_plus_transcripts,
        }
    }

    /// Get the number of transcripts in the database
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Check if the database is empty
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Get all transcript IDs
    pub fn ids(&self) -> impl Iterator<Item = &String> {
        self.transcripts.keys()
    }

    /// Iterate over all transcripts
    pub fn iter(&self) -> impl Iterator<Item = (&String, &Transcript)> {
        self.transcripts.iter()
    }
}

/// Detect genome build from file header or content
pub fn detect_genome_build<P: AsRef<Path>>(path: P) -> Result<GenomeBuild, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open file: {}", e),
    })?;
    let reader = BufReader::new(file);

    for line in reader.lines().take(100) {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        // Check for genome build indicators
        let line_lower = line.to_lowercase();
        if line_lower.contains("grch37") || line_lower.contains("hg19") {
            return Ok(GenomeBuild::GRCh37);
        }
        if line_lower.contains("grch38") || line_lower.contains("hg38") {
            return Ok(GenomeBuild::GRCh38);
        }
    }

    Ok(GenomeBuild::Unknown)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, Strand};
    use std::sync::OnceLock;

    #[test]
    fn test_transcript_db_new() {
        let db = TranscriptDb::new();
        assert!(db.is_empty());
    }

    #[test]
    fn test_transcript_db_add_and_get() {
        let mut db = TranscriptDb::new();

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        assert_eq!(db.len(), 1);
        assert!(db.get("NM_000088.3").is_some());
    }

    #[test]
    fn test_transcript_db_get_by_gene() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let transcripts = db.get_by_gene("COL1A1");
        assert_eq!(transcripts.len(), 2);
    }

    #[test]
    fn test_transcript_db_get_by_position() {
        let mut db = TranscriptDb::new();

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Position within range
        let matches = db.get_by_position("chr1", 1500);
        assert_eq!(matches.len(), 1);

        // Position outside range
        let matches = db.get_by_position("chr1", 500);
        assert_eq!(matches.len(), 0);

        // Different chromosome
        let matches = db.get_by_position("chr2", 1500);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn test_mane_select_lookup() {
        let mut db = TranscriptDb::new();

        let tx_mane = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx_other = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_mane);
        db.add(tx_other);

        // Should find MANE Select
        assert!(db.has_mane_select("COL1A1"));
        let mane = db.get_mane_select("COL1A1");
        assert!(mane.is_some());
        assert_eq!(mane.unwrap().id, "NM_000088.3");

        // Preferred transcript should be MANE Select
        let preferred = db.get_preferred_transcript("COL1A1");
        assert!(preferred.is_some());
        assert_eq!(preferred.unwrap().id, "NM_000088.3");
    }

    #[test]
    fn test_mane_stats() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000099.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let stats = db.mane_stats();
        assert_eq!(stats.total_genes, 2);
        assert_eq!(stats.genes_with_mane_select, 1);
        assert!((stats.mane_select_coverage() - 50.0).abs() < 0.01);
    }

    // ===== ManeStats Tests =====

    #[test]
    fn test_mane_stats_default() {
        let stats = ManeStats::default();
        assert_eq!(stats.total_genes, 0);
        assert_eq!(stats.genes_with_mane_select, 0);
        assert_eq!(stats.genes_with_mane_plus_clinical, 0);
        assert_eq!(stats.total_mane_plus_clinical_transcripts, 0);
    }

    #[test]
    fn test_mane_stats_coverage_empty() {
        let stats = ManeStats::default();
        assert_eq!(stats.mane_select_coverage(), 0.0);
    }

    #[test]
    fn test_mane_stats_coverage_full() {
        let stats = ManeStats {
            total_genes: 100,
            genes_with_mane_select: 100,
            genes_with_mane_plus_clinical: 50,
            total_mane_plus_clinical_transcripts: 60,
        };
        assert!((stats.mane_select_coverage() - 100.0).abs() < 0.01);
    }

    // ===== TranscriptDb Extended Tests =====

    #[test]
    fn test_transcript_db_with_build() {
        let db = TranscriptDb::with_build(GenomeBuild::GRCh37);
        assert!(matches!(db.genome_build, GenomeBuild::GRCh37));
    }

    #[test]
    fn test_transcript_db_get_by_region() {
        let mut db = TranscriptDb::new();

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Overlapping region
        let matches = db.get_by_region("chr1", 1500, 1600);
        assert_eq!(matches.len(), 1);

        // Adjacent region (no overlap)
        let matches = db.get_by_region("chr1", 500, 999);
        assert_eq!(matches.len(), 0);

        // Fully containing region
        let matches = db.get_by_region("chr1", 900, 2100);
        assert_eq!(matches.len(), 1);

        // Different chromosome
        let matches = db.get_by_region("chr2", 1000, 2000);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn test_transcript_db_mane_plus_clinical() {
        let mut db = TranscriptDb::new();

        let tx_plus1 = Transcript {
            id: "NM_000100.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx_plus2 = Transcript {
            id: "NM_000101.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_plus1);
        db.add(tx_plus2);

        let plus_clinical = db.get_mane_plus_clinical("BRCA1");
        assert_eq!(plus_clinical.len(), 2);

        // No MANE Select for this gene
        assert!(db.get_mane_select("BRCA1").is_none());

        // Get all MANE transcripts
        let all_mane = db.get_mane_transcripts("BRCA1");
        assert_eq!(all_mane.len(), 2);
    }

    #[test]
    fn test_transcript_db_preferred_transcript_fallback() {
        let mut db = TranscriptDb::new();

        // Non-MANE, non-coding transcript
        let tx_noncoding = Transcript {
            id: "NR_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        // Coding transcript (preferred over non-coding)
        let tx_coding = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_noncoding);
        db.add(tx_coding);

        // Should prefer coding transcript
        let preferred = db.get_preferred_transcript("TEST");
        assert!(preferred.is_some());
        assert!(preferred.unwrap().is_coding());
    }

    #[test]
    fn test_transcript_db_ids_iter() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("GENE1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("GENE2".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let ids: Vec<_> = db.ids().collect();
        assert_eq!(ids.len(), 2);

        let iter_count = db.iter().count();
        assert_eq!(iter_count, 2);
    }

    // ===== Gene Lookup Edge Cases =====

    #[test]
    fn test_transcript_db_nonexistent_gene() {
        let db = TranscriptDb::new();

        assert!(db.get_by_gene("NONEXISTENT").is_empty());
        assert!(db.get_mane_select("NONEXISTENT").is_none());
        assert!(db.get_mane_plus_clinical("NONEXISTENT").is_empty());
        assert!(db.get_mane_transcripts("NONEXISTENT").is_empty());
        assert!(db.get_preferred_transcript("NONEXISTENT").is_none());
        assert!(!db.has_mane_select("NONEXISTENT"));
    }

    #[test]
    fn test_transcript_db_no_gene_symbol() {
        let mut db = TranscriptDb::new();

        // Transcript without gene symbol
        let tx = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Should still be findable by ID and position
        assert!(db.get("NM_000001.1").is_some());
        assert!(!db.get_by_position("chr1", 1500).is_empty());
    }

    // ===== Additional TranscriptDb Tests =====

    #[test]
    fn test_get_mane_transcripts() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let mane_txs = db.get_mane_transcripts("BRCA1");
        assert_eq!(mane_txs.len(), 2); // Both MANE Select and Plus Clinical
    }

    #[test]
    fn test_get_preferred_transcript_priority() {
        let mut db = TranscriptDb::new();

        // Add transcripts in reverse priority order
        let tx_none = Transcript {
            id: "NM_000003.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_none);

        // Without MANE, should return the coding transcript
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000003.1");

        // Add MANE Plus Clinical
        let tx_plus = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        db.add(tx_plus);

        // Should prefer MANE Plus Clinical
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000002.1");

        // Add MANE Select
        let tx_select = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        db.add(tx_select);

        // Should prefer MANE Select
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000001.1");
    }
}
