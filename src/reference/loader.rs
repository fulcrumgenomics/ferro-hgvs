//! Transcript database loading
//!
//! This module provides functionality to load transcript annotations from
//! standard file formats (GFF3, GTF) commonly used for genome annotations.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::OnceLock;

use crate::error::FerroError;
use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};

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

/// GFF3 attribute parser
struct Gff3Attributes {
    id: Option<String>,
    parent: Option<String>,
    name: Option<String>,
    gene_name: Option<String>,
    transcript_id: Option<String>,
    mane_status: ManeStatus,
    refseq_match: Option<String>,
    ensembl_match: Option<String>,
}

impl Gff3Attributes {
    fn parse(attr_str: &str) -> Self {
        let mut attrs = Self {
            id: None,
            parent: None,
            name: None,
            gene_name: None,
            transcript_id: None,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
        };

        for part in attr_str.split(';') {
            let part = part.trim();
            if let Some((key, value)) = part.split_once('=') {
                let value = url_decode(value);
                match key {
                    "ID" => attrs.id = Some(value),
                    "Parent" => attrs.parent = Some(value),
                    "Name" => attrs.name = Some(value),
                    "gene" | "gene_name" => attrs.gene_name = Some(value),
                    "transcript_id" => attrs.transcript_id = Some(value),
                    // MANE status parsing - various formats used in GFF3/GTF files
                    "tag" => {
                        if value.contains("MANE_Select") || value.contains("MANE Select") {
                            attrs.mane_status = ManeStatus::Select;
                        } else if value.contains("MANE_Plus_Clinical")
                            || value.contains("MANE Plus Clinical")
                        {
                            attrs.mane_status = ManeStatus::PlusClinical;
                        }
                    }
                    "MANE" | "mane_status" => {
                        let lower = value.to_lowercase();
                        if lower.contains("select") {
                            attrs.mane_status = ManeStatus::Select;
                        } else if lower.contains("plus") || lower.contains("clinical") {
                            attrs.mane_status = ManeStatus::PlusClinical;
                        }
                    }
                    // Cross-reference accessions
                    "RefSeq" | "refseq_id" => attrs.refseq_match = Some(value),
                    "Ensembl" | "ensembl_id" | "ENST" => attrs.ensembl_match = Some(value),
                    _ => {}
                }
            }
        }

        attrs
    }
}

/// Simple URL decoding for GFF3 attributes
fn url_decode(s: &str) -> String {
    let mut result = String::with_capacity(s.len());
    let mut chars = s.chars();

    while let Some(c) = chars.next() {
        if c == '%' {
            let hex: String = chars.by_ref().take(2).collect();
            if let Ok(byte) = u8::from_str_radix(&hex, 16) {
                result.push(byte as char);
            } else {
                result.push('%');
                result.push_str(&hex);
            }
        } else {
            result.push(c);
        }
    }

    result
}

/// Parse strand from a single character
fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        _ => Strand::Plus, // Default to plus for unknown
    }
}

/// Load transcripts from a GFF3 file
///
/// GFF3 format has tab-separated columns:
/// 1. seqid (chromosome)
/// 2. source
/// 3. type (gene, mRNA, exon, CDS, etc.)
/// 4. start (1-based)
/// 5. end (1-based, inclusive)
/// 6. score
/// 7. strand (+/-)
/// 8. phase
/// 9. attributes (ID=..;Name=..;Parent=..)
pub fn load_gff3<P: AsRef<Path>>(path: P) -> Result<TranscriptDb, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open GFF3 file: {}", e),
    })?;
    let reader = BufReader::new(file);

    let mut db = TranscriptDb::new();

    // Track transcripts being built
    let mut tx_builders: HashMap<String, TranscriptBuilder> = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        // Skip comments and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            // Check for genome build in header
            if line.contains("genome-build") {
                if line.contains("GRCh37") || line.contains("hg19") {
                    db.genome_build = GenomeBuild::GRCh37;
                } else if line.contains("GRCh38") || line.contains("hg38") {
                    db.genome_build = GenomeBuild::GRCh38;
                }
            }
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let seqid = fields[0];
        let feature_type = fields[2];
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand = parse_strand(fields[6]);
        let attrs = Gff3Attributes::parse(fields[8]);

        match feature_type {
            "mRNA" | "transcript" | "primary_transcript" => {
                if let Some(id) = attrs.id.or(attrs.transcript_id) {
                    let builder = TranscriptBuilder {
                        id: id.clone(),
                        gene_symbol: attrs.gene_name.or(attrs.name),
                        chromosome: Some(seqid.to_string()),
                        strand,
                        genomic_start: start,
                        genomic_end: end,
                        exons: Vec::new(),
                        cds_ranges: Vec::new(),
                        mane_status: attrs.mane_status,
                        refseq_match: attrs.refseq_match,
                        ensembl_match: attrs.ensembl_match,
                    };
                    tx_builders.insert(id, builder);
                }
            }
            "exon" => {
                if let Some(parent) = attrs.parent {
                    // Handle multiple parents (comma-separated)
                    for parent_id in parent.split(',') {
                        if let Some(builder) = tx_builders.get_mut(parent_id) {
                            builder.exons.push((start, end));
                        }
                    }
                }
            }
            "CDS" => {
                if let Some(parent) = attrs.parent {
                    for parent_id in parent.split(',') {
                        if let Some(builder) = tx_builders.get_mut(parent_id) {
                            builder.cds_ranges.push((start, end));
                        }
                    }
                }
            }
            _ => {}
        }
    }

    // Build transcripts from builders
    for (_, builder) in tx_builders {
        if let Some(transcript) = builder.build(db.genome_build) {
            db.add(transcript);
        }
    }

    Ok(db)
}

/// Load transcripts from a GTF file
///
/// GTF format is similar to GFF but uses different attribute format
pub fn load_gtf<P: AsRef<Path>>(path: P) -> Result<TranscriptDb, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open GTF file: {}", e),
    })?;
    let reader = BufReader::new(file);

    let mut db = TranscriptDb::new();
    let mut tx_builders: HashMap<String, TranscriptBuilder> = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let seqid = fields[0];
        let feature_type = fields[2];
        let start: u64 = fields[3].parse().unwrap_or(0);
        let end: u64 = fields[4].parse().unwrap_or(0);
        let strand = parse_strand(fields[6]);
        let attrs = parse_gtf_attributes(fields[8]);

        let transcript_id = attrs.get("transcript_id").cloned();
        let gene_name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("gene_id"))
            .cloned();

        // Parse MANE status from GTF "tag" attribute
        let mane_status = attrs
            .get("tag")
            .map(|tag| {
                if tag.contains("MANE_Select") || tag.contains("MANE Select") {
                    ManeStatus::Select
                } else if tag.contains("MANE_Plus_Clinical") || tag.contains("MANE Plus Clinical") {
                    ManeStatus::PlusClinical
                } else {
                    ManeStatus::None
                }
            })
            .unwrap_or(ManeStatus::None);

        match feature_type {
            "transcript" => {
                if let Some(id) = transcript_id {
                    let builder = TranscriptBuilder {
                        id: id.clone(),
                        gene_symbol: gene_name,
                        chromosome: Some(seqid.to_string()),
                        strand,
                        genomic_start: start,
                        genomic_end: end,
                        exons: Vec::new(),
                        cds_ranges: Vec::new(),
                        mane_status,
                        refseq_match: None,
                        ensembl_match: None,
                    };
                    tx_builders.insert(id, builder);
                }
            }
            "exon" => {
                if let Some(id) = transcript_id {
                    // Create transcript if not seen yet
                    let builder =
                        tx_builders
                            .entry(id.clone())
                            .or_insert_with(|| TranscriptBuilder {
                                id,
                                gene_symbol: gene_name.clone(),
                                chromosome: Some(seqid.to_string()),
                                strand,
                                genomic_start: start,
                                genomic_end: end,
                                exons: Vec::new(),
                                cds_ranges: Vec::new(),
                                mane_status,
                                refseq_match: None,
                                ensembl_match: None,
                            });

                    // Update bounds
                    builder.genomic_start = builder.genomic_start.min(start);
                    builder.genomic_end = builder.genomic_end.max(end);
                    builder.exons.push((start, end));
                }
            }
            "CDS" => {
                if let Some(id) = transcript_id {
                    if let Some(builder) = tx_builders.get_mut(&id) {
                        builder.cds_ranges.push((start, end));
                    }
                }
            }
            _ => {}
        }
    }

    // Build transcripts
    for (_, builder) in tx_builders {
        if let Some(transcript) = builder.build(db.genome_build) {
            db.add(transcript);
        }
    }

    Ok(db)
}

/// Parse GTF attribute string
fn parse_gtf_attributes(attr_str: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();

    for part in attr_str.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }

        // GTF format: key "value"
        let mut iter = part.splitn(2, ' ');
        if let (Some(key), Some(value)) = (iter.next(), iter.next()) {
            // Remove quotes from value
            let value = value.trim_matches('"').to_string();
            attrs.insert(key.to_string(), value);
        }
    }

    attrs
}

/// Builder for constructing transcripts from annotation features
struct TranscriptBuilder {
    id: String,
    gene_symbol: Option<String>,
    chromosome: Option<String>,
    strand: Strand,
    genomic_start: u64,
    genomic_end: u64,
    exons: Vec<(u64, u64)>,
    cds_ranges: Vec<(u64, u64)>,
    mane_status: ManeStatus,
    refseq_match: Option<String>,
    ensembl_match: Option<String>,
}

impl TranscriptBuilder {
    fn build(mut self, genome_build: GenomeBuild) -> Option<Transcript> {
        if self.exons.is_empty() {
            return None;
        }

        // Sort exons by position
        self.exons.sort_by_key(|(start, _)| *start);

        // Calculate CDS bounds
        let (cds_start, cds_end) = if !self.cds_ranges.is_empty() {
            let cds_start = self.cds_ranges.iter().map(|(s, _)| *s).min().unwrap();
            let cds_end = self.cds_ranges.iter().map(|(_, e)| *e).max().unwrap();
            (Some(cds_start), Some(cds_end))
        } else {
            (None, None)
        };

        // Build exon structures
        let mut tx_pos = 1u64;
        let exons: Vec<Exon> = self
            .exons
            .iter()
            .enumerate()
            .map(|(i, (g_start, g_end))| {
                let exon_len = g_end - g_start + 1;
                let exon = Exon::with_genomic(
                    (i + 1) as u32,
                    tx_pos,
                    tx_pos + exon_len - 1,
                    *g_start,
                    *g_end,
                );
                tx_pos += exon_len;
                exon
            })
            .collect();

        // Calculate CDS start/end in transcript coordinates
        let (tx_cds_start, tx_cds_end) =
            if let (Some(g_cds_start), Some(g_cds_end)) = (cds_start, cds_end) {
                let mut cds_tx_start = None;
                let mut cds_tx_end = None;

                for exon in &exons {
                    if let (Some(g_start), Some(g_end)) = (exon.genomic_start, exon.genomic_end) {
                        // Check if CDS start is in this exon
                        if g_cds_start >= g_start && g_cds_start <= g_end {
                            let offset = g_cds_start - g_start;
                            cds_tx_start = Some(exon.start + offset);
                        }
                        // Check if CDS end is in this exon
                        if g_cds_end >= g_start && g_cds_end <= g_end {
                            let offset = g_cds_end - g_start;
                            cds_tx_end = Some(exon.start + offset);
                        }
                    }
                }

                (cds_tx_start, cds_tx_end)
            } else {
                (None, None)
            };

        // Calculate sequence length from exons
        let seq_len: u64 = exons.iter().map(|e| e.end - e.start + 1).sum();

        Some(Transcript {
            id: self.id,
            gene_symbol: self.gene_symbol,
            strand: self.strand,
            sequence: "N".repeat(seq_len as usize),
            cds_start: tx_cds_start,
            cds_end: tx_cds_end,
            exons,
            chromosome: self.chromosome,
            genomic_start: Some(self.genomic_start),
            genomic_end: Some(self.genomic_end),
            genome_build,
            mane_status: self.mane_status,
            refseq_match: self.refseq_match,
            ensembl_match: self.ensembl_match,
            cached_introns: OnceLock::new(),
        })
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
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        let tx_other = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000099.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let stats = db.mane_stats();
        assert_eq!(stats.total_genes, 2);
        assert_eq!(stats.genes_with_mane_select, 1);
        assert!((stats.mane_select_coverage() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_parse_gtf_attributes() {
        let attrs = parse_gtf_attributes(
            r#"gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_name "TEST";"#,
        );

        assert_eq!(attrs.get("gene_id"), Some(&"ENSG00000000001".to_string()));
        assert_eq!(
            attrs.get("transcript_id"),
            Some(&"ENST00000000001".to_string())
        );
        assert_eq!(attrs.get("gene_name"), Some(&"TEST".to_string()));
    }

    #[test]
    fn test_url_decode() {
        assert_eq!(url_decode("hello%20world"), "hello world");
        assert_eq!(url_decode("no%2Fslash"), "no/slash");
        assert_eq!(url_decode("plain"), "plain");
    }

    #[test]
    fn test_parse_strand() {
        assert_eq!(parse_strand("+"), Strand::Plus);
        assert_eq!(parse_strand("-"), Strand::Minus);
        assert_eq!(parse_strand("."), Strand::Plus); // Default
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
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        let tx_plus2 = Transcript {
            id: "NM_000101.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        // Coding transcript (preferred over non-coding)
        let tx_coding = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("GENE2".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let ids: Vec<_> = db.ids().collect();
        assert_eq!(ids.len(), 2);

        let iter_count = db.iter().count();
        assert_eq!(iter_count, 2);
    }

    // ===== Gff3Attributes Tests =====

    #[test]
    fn test_gff3_attributes_parse_basic() {
        let attrs = Gff3Attributes::parse("ID=gene123;Name=TestGene;Parent=parent123");
        assert_eq!(attrs.id, Some("gene123".to_string()));
        assert_eq!(attrs.name, Some("TestGene".to_string()));
        assert_eq!(attrs.parent, Some("parent123".to_string()));
    }

    #[test]
    fn test_gff3_attributes_parse_gene_name() {
        let attrs = Gff3Attributes::parse("ID=tx1;gene=BRCA1;gene_name=BRCA1");
        assert_eq!(attrs.gene_name, Some("BRCA1".to_string()));
    }

    #[test]
    fn test_gff3_attributes_parse_mane_select() {
        let attrs = Gff3Attributes::parse("ID=tx1;tag=MANE_Select");
        assert!(matches!(attrs.mane_status, ManeStatus::Select));
    }

    #[test]
    fn test_gff3_attributes_parse_mane_plus_clinical() {
        let attrs = Gff3Attributes::parse("ID=tx1;tag=MANE_Plus_Clinical");
        assert!(matches!(attrs.mane_status, ManeStatus::PlusClinical));
    }

    #[test]
    fn test_gff3_attributes_parse_mane_attribute() {
        let attrs = Gff3Attributes::parse("ID=tx1;MANE=Select");
        assert!(matches!(attrs.mane_status, ManeStatus::Select));

        let attrs = Gff3Attributes::parse("ID=tx1;mane_status=plus_clinical");
        assert!(matches!(attrs.mane_status, ManeStatus::PlusClinical));
    }

    #[test]
    fn test_gff3_attributes_parse_cross_references() {
        let attrs = Gff3Attributes::parse("ID=tx1;RefSeq=NM_000088.3;Ensembl=ENST00000000001");
        assert_eq!(attrs.refseq_match, Some("NM_000088.3".to_string()));
        assert_eq!(attrs.ensembl_match, Some("ENST00000000001".to_string()));
    }

    #[test]
    fn test_gff3_attributes_parse_url_encoded() {
        let attrs = Gff3Attributes::parse("ID=tx1;Name=Test%20Gene");
        assert_eq!(attrs.name, Some("Test Gene".to_string()));
    }

    // ===== URL Decode Extended Tests =====

    #[test]
    fn test_url_decode_special_chars() {
        assert_eq!(url_decode("%3D"), "=");
        assert_eq!(url_decode("%3E"), ">");
        assert_eq!(url_decode("%3C"), "<");
        assert_eq!(url_decode("%26"), "&");
    }

    #[test]
    fn test_url_decode_invalid() {
        // Invalid hex should pass through
        assert_eq!(url_decode("%ZZ"), "%ZZ");
        assert_eq!(url_decode("%"), "%");
        // Single hex digit: "2" is valid hex (0x02), so it decodes to that character
        assert_eq!(url_decode("%2"), "\u{02}");
    }

    // ===== TranscriptBuilder Tests =====

    #[test]
    fn test_transcript_builder_empty_exons() {
        let builder = TranscriptBuilder {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            chromosome: Some("chr1".to_string()),
            strand: Strand::Plus,
            genomic_start: 1000,
            genomic_end: 2000,
            exons: Vec::new(), // Empty
            cds_ranges: Vec::new(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
        };

        // Should return None for empty exons
        assert!(builder.build(GenomeBuild::GRCh38).is_none());
    }

    #[test]
    fn test_transcript_builder_with_cds() {
        let builder = TranscriptBuilder {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            chromosome: Some("chr1".to_string()),
            strand: Strand::Plus,
            genomic_start: 1000,
            genomic_end: 2000,
            exons: vec![(1000, 1100), (1500, 1600), (1800, 2000)],
            cds_ranges: vec![(1050, 1100), (1500, 1550)],
            mane_status: ManeStatus::Select,
            refseq_match: Some("NM_000001.1".to_string()),
            ensembl_match: Some("ENST00000001".to_string()),
        };

        let transcript = builder.build(GenomeBuild::GRCh38);
        assert!(transcript.is_some());

        let tx = transcript.unwrap();
        assert_eq!(tx.id, "NM_000001.1");
        assert_eq!(tx.gene_symbol, Some("TEST".to_string()));
        assert!(tx.is_coding());
        assert_eq!(tx.exons.len(), 3);
        assert!(matches!(tx.mane_status, ManeStatus::Select));
    }

    #[test]
    fn test_transcript_builder_unsorted_exons() {
        let builder = TranscriptBuilder {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            chromosome: Some("chr1".to_string()),
            strand: Strand::Plus,
            genomic_start: 1000,
            genomic_end: 2000,
            exons: vec![(1500, 1600), (1000, 1100), (1800, 2000)], // Unsorted
            cds_ranges: Vec::new(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
        };

        let transcript = builder.build(GenomeBuild::GRCh38);
        assert!(transcript.is_some());

        let tx = transcript.unwrap();
        // Exons should be sorted by position
        assert!(tx.exons[0].genomic_start.unwrap() < tx.exons[1].genomic_start.unwrap());
    }

    // ===== GTF Attributes Extended Tests =====

    #[test]
    fn test_parse_gtf_attributes_empty() {
        let attrs = parse_gtf_attributes("");
        assert!(attrs.is_empty());
    }

    #[test]
    fn test_parse_gtf_attributes_with_mane_tag() {
        let attrs = parse_gtf_attributes(
            r#"gene_id "ENSG00000000001"; transcript_id "ENST00000000001"; tag "MANE_Select";"#,
        );
        assert_eq!(attrs.get("tag"), Some(&"MANE_Select".to_string()));
    }

    #[test]
    fn test_parse_gtf_attributes_whitespace() {
        let attrs = parse_gtf_attributes(r#"  gene_id "GENE1" ;  transcript_id "TX1"  ;"#);
        assert_eq!(attrs.get("gene_id"), Some(&"GENE1".to_string()));
        assert_eq!(attrs.get("transcript_id"), Some(&"TX1".to_string()));
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Should still be findable by ID and position
        assert!(db.get("NM_000001.1").is_some());
        assert!(!db.get_by_position("chr1", 1500).is_empty());
    }

    // ===== I/O Tests with Temporary Files =====

    #[test]
    fn test_load_gff3_basic() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(gff3_file, "##genome-build GRCh38").unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;Name=TEST;gene=TESTGENE"
        )
        .unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\texon\t1000\t1100\t.\t+\t.\tID=exon1;Parent=tx1"
        )
        .unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\texon\t1500\t1600\t.\t+\t.\tID=exon2;Parent=tx1"
        )
        .unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tCDS\t1050\t1100\t.\t+\t0\tID=cds1;Parent=tx1"
        )
        .unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();

        assert_eq!(db.len(), 1);
        assert_eq!(db.genome_build, GenomeBuild::GRCh38);

        let tx = db.get("tx1").unwrap();
        assert_eq!(tx.gene_symbol, Some("TESTGENE".to_string()));
        assert_eq!(tx.exons.len(), 2);
        assert!(tx.is_coding());
    }

    #[test]
    fn test_load_gff3_with_mane() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;gene=BRCA1;tag=MANE_Select"
        )
        .unwrap();
        writeln!(gff3_file, "chr1\t.\texon\t1000\t2000\t.\t+\t.\tParent=tx1").unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();

        let tx = db.get("tx1").unwrap();
        assert!(matches!(tx.mane_status, ManeStatus::Select));
        assert!(db.has_mane_select("BRCA1"));
    }

    #[test]
    fn test_load_gff3_primary_transcript() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tprimary_transcript\t1000\t2000\t.\t+\t.\tID=tx1;gene=TEST"
        )
        .unwrap();
        writeln!(gff3_file, "chr1\t.\texon\t1000\t2000\t.\t+\t.\tParent=tx1").unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();
        assert_eq!(db.len(), 1);
    }

    #[test]
    fn test_load_gff3_transcript_type() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\ttranscript\t1000\t2000\t.\t+\t.\tID=tx1;gene=TEST"
        )
        .unwrap();
        writeln!(gff3_file, "chr1\t.\texon\t1000\t2000\t.\t+\t.\tParent=tx1").unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();
        assert_eq!(db.len(), 1);
    }

    #[test]
    fn test_load_gff3_grch37() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(gff3_file, "##genome-build GRCh37").unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;gene=TEST"
        )
        .unwrap();
        writeln!(gff3_file, "chr1\t.\texon\t1000\t2000\t.\t+\t.\tParent=tx1").unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();
        assert_eq!(db.genome_build, GenomeBuild::GRCh37);
    }

    #[test]
    fn test_load_gff3_multiple_parents() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;gene=TEST"
        )
        .unwrap();
        writeln!(
            gff3_file,
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx2;gene=TEST"
        )
        .unwrap();
        // Exon belongs to both transcripts
        writeln!(
            gff3_file,
            "chr1\t.\texon\t1000\t2000\t.\t+\t.\tParent=tx1,tx2"
        )
        .unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();

        // Both transcripts should have the exon
        let tx1 = db.get("tx1").unwrap();
        let tx2 = db.get("tx2").unwrap();
        assert_eq!(tx1.exons.len(), 1);
        assert_eq!(tx2.exons.len(), 1);
    }

    #[test]
    fn test_load_gtf_basic() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gtf_file = NamedTempFile::new().unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	transcript	1000	2000	.	+	.	gene_id "TESTGENE"; transcript_id "tx1";"#
        )
        .unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1000	1100	.	+	.	gene_id "TESTGENE"; transcript_id "tx1";"#
        )
        .unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1500	1600	.	+	.	gene_id "TESTGENE"; transcript_id "tx1";"#
        )
        .unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	CDS	1050	1100	.	+	0	gene_id "TESTGENE"; transcript_id "tx1";"#
        )
        .unwrap();
        gtf_file.flush().unwrap();

        let db = load_gtf(gtf_file.path()).unwrap();

        assert_eq!(db.len(), 1);

        let tx = db.get("tx1").unwrap();
        assert_eq!(tx.gene_symbol, Some("TESTGENE".to_string()));
        assert_eq!(tx.exons.len(), 2);
        assert!(tx.is_coding());
    }

    #[test]
    fn test_load_gtf_with_mane() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gtf_file = NamedTempFile::new().unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	transcript	1000	2000	.	+	.	gene_id "BRCA1"; transcript_id "tx1"; tag "MANE_Select";"#
        )
        .unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1000	2000	.	+	.	gene_id "BRCA1"; transcript_id "tx1";"#
        )
        .unwrap();
        gtf_file.flush().unwrap();

        let db = load_gtf(gtf_file.path()).unwrap();

        let tx = db.get("tx1").unwrap();
        assert!(matches!(tx.mane_status, ManeStatus::Select));
    }

    #[test]
    fn test_load_gtf_exon_before_transcript() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // In some GTF files, exons appear before the transcript feature
        let mut gtf_file = NamedTempFile::new().unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1000	1100	.	+	.	gene_id "TESTGENE"; transcript_id "tx1";"#
        )
        .unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1500	1600	.	+	.	gene_id "TESTGENE"; transcript_id "tx1";"#
        )
        .unwrap();
        gtf_file.flush().unwrap();

        let db = load_gtf(gtf_file.path()).unwrap();

        // Transcript should still be created from exon features
        assert_eq!(db.len(), 1);
        let tx = db.get("tx1").unwrap();
        assert_eq!(tx.exons.len(), 2);
    }

    #[test]
    fn test_load_gtf_gene_name_fallback() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        // Use gene_id when gene_name is not present
        let mut gtf_file = NamedTempFile::new().unwrap();
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1000	2000	.	+	.	gene_id "GENEID"; transcript_id "tx1";"#
        )
        .unwrap();
        gtf_file.flush().unwrap();

        let db = load_gtf(gtf_file.path()).unwrap();

        let tx = db.get("tx1").unwrap();
        assert_eq!(tx.gene_symbol, Some("GENEID".to_string()));
    }

    #[test]
    fn test_load_gtf_skip_comments() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gtf_file = NamedTempFile::new().unwrap();
        writeln!(gtf_file, "# Comment line").unwrap();
        writeln!(gtf_file).unwrap(); // Empty line
        writeln!(
            gtf_file,
            r#"chr1	.	exon	1000	2000	.	+	.	gene_id "TEST"; transcript_id "tx1";"#
        )
        .unwrap();
        gtf_file.flush().unwrap();

        let db = load_gtf(gtf_file.path()).unwrap();
        assert_eq!(db.len(), 1);
    }

    #[test]
    fn test_load_gff3_skip_invalid_lines() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut gff3_file = NamedTempFile::new().unwrap();
        writeln!(gff3_file, "##gff-version 3").unwrap();
        writeln!(gff3_file, "short_line").unwrap(); // Invalid - too few fields
        writeln!(
            gff3_file,
            "chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=tx1;gene=TEST"
        )
        .unwrap();
        writeln!(gff3_file, "chr1\t.\texon\t1000\t2000\t.\t+\t.\tParent=tx1").unwrap();
        gff3_file.flush().unwrap();

        let db = load_gff3(gff3_file.path()).unwrap();
        assert_eq!(db.len(), 1); // Invalid line should be skipped
    }

    // ===== Additional TranscriptDb Tests =====

    #[test]
    fn test_get_mane_transcripts() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            sequence: "ATGC".to_string(),
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
            cached_introns: OnceLock::new(),
        };
        db.add(tx_select);

        // Should prefer MANE Select
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000001.1");
    }

    #[test]
    fn test_parse_strand_all_cases() {
        assert!(matches!(parse_strand("+"), Strand::Plus));
        assert!(matches!(parse_strand("-"), Strand::Minus));
        assert!(matches!(parse_strand("."), Strand::Plus)); // Default
        assert!(matches!(parse_strand("?"), Strand::Plus)); // Default
    }

    #[test]
    fn test_transcript_builder_minus_strand() {
        let builder = TranscriptBuilder {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            chromosome: Some("chr1".to_string()),
            strand: Strand::Minus,
            genomic_start: 1000,
            genomic_end: 2000,
            exons: vec![(1800, 2000), (1500, 1600), (1000, 1100)], // Minus strand order
            cds_ranges: Vec::new(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
        };

        let transcript = builder.build(GenomeBuild::GRCh38);
        assert!(transcript.is_some());

        let tx = transcript.unwrap();
        assert!(matches!(tx.strand, Strand::Minus));
    }
}
