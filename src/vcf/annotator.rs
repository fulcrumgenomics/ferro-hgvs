//! Multi-isoform variant annotator
//!
//! This module provides functionality to annotate VCF records with HGVS
//! notations from all overlapping transcripts.

use std::collections::HashMap;

use crate::error::FerroError;
use crate::reference::loader::TranscriptDb;

use super::record::VcfRecord;
use super::to_hgvs::{HgvsAnnotation, VcfToHgvsConverter};

/// Annotation result for a single variant across all transcripts
#[derive(Debug, Clone)]
pub struct MultiIsoformAnnotation {
    /// Original VCF record
    pub vcf: VcfRecord,
    /// All annotations grouped by gene symbol
    pub annotations_by_gene: HashMap<String, Vec<HgvsAnnotation>>,
    /// All annotations (flat list)
    pub annotations: Vec<HgvsAnnotation>,
    /// Warnings during annotation
    pub warnings: Vec<String>,
}

impl MultiIsoformAnnotation {
    /// Get annotations for a specific gene
    pub fn get_gene(&self, gene: &str) -> Option<&Vec<HgvsAnnotation>> {
        self.annotations_by_gene.get(gene)
    }

    /// Get all gene symbols that have annotations
    pub fn genes(&self) -> impl Iterator<Item = &String> {
        self.annotations_by_gene.keys()
    }

    /// Get the number of annotations
    pub fn annotation_count(&self) -> usize {
        self.annotations.len()
    }

    /// Get the primary (MANE/canonical) annotation if identifiable
    ///
    /// Priority order:
    /// 1. MANE Select transcript
    /// 2. MANE Plus Clinical transcript
    /// 3. Any coding transcript (alphabetically first)
    /// 4. Any transcript (first available)
    pub fn primary_annotation(&self) -> Option<&HgvsAnnotation> {
        // Priority 1: MANE Select
        if let Some(mane_select) = self.annotations.iter().find(|a| a.is_mane_select()) {
            return Some(mane_select);
        }

        // Priority 2: MANE Plus Clinical
        if let Some(mane_plus) = self.annotations.iter().find(|a| a.is_mane_plus_clinical()) {
            return Some(mane_plus);
        }

        // Priority 3: Any coding transcript (alphabetically first by accession)
        if let Some(coding) = self
            .annotations
            .iter()
            .filter(|a| a.is_coding)
            .min_by_key(|a| a.transcript_accession.as_ref())
        {
            return Some(coding);
        }

        // Priority 4: Any annotation
        self.annotations.first()
    }

    /// Get all MANE annotations (Select and Plus Clinical)
    pub fn mane_annotations(&self) -> Vec<&HgvsAnnotation> {
        self.annotations.iter().filter(|a| a.is_mane()).collect()
    }

    /// Get the MANE Select annotation if available
    pub fn mane_select_annotation(&self) -> Option<&HgvsAnnotation> {
        self.annotations.iter().find(|a| a.is_mane_select())
    }

    /// Get all HGVS strings as a list
    pub fn hgvs_strings(&self) -> Vec<String> {
        self.annotations.iter().map(|a| a.hgvs_string()).collect()
    }

    /// Get HGVS strings grouped by gene
    pub fn hgvs_by_gene(&self) -> HashMap<String, Vec<String>> {
        self.annotations_by_gene
            .iter()
            .map(|(gene, anns)| (gene.clone(), anns.iter().map(|a| a.hgvs_string()).collect()))
            .collect()
    }
}

/// Multi-isoform variant annotator
///
/// Annotates VCF records with HGVS notation from all overlapping transcripts.
pub struct MultiIsoformAnnotator<'a> {
    /// Transcript database
    db: &'a TranscriptDb,
    /// Whether to include non-coding transcripts
    include_non_coding: bool,
    /// Whether to include intronic variants
    include_intronic: bool,
}

impl<'a> MultiIsoformAnnotator<'a> {
    /// Create a new annotator
    pub fn new(db: &'a TranscriptDb) -> Self {
        Self {
            db,
            include_non_coding: true,
            include_intronic: true,
        }
    }

    /// Configure whether to include non-coding transcripts
    pub fn include_non_coding(mut self, include: bool) -> Self {
        self.include_non_coding = include;
        self
    }

    /// Configure whether to include intronic variants
    pub fn include_intronic(mut self, include: bool) -> Self {
        self.include_intronic = include;
        self
    }

    /// Annotate a VCF record with all overlapping transcripts
    pub fn annotate(&self, vcf: &VcfRecord) -> Result<MultiIsoformAnnotation, FerroError> {
        let mut annotations = Vec::new();
        let mut annotations_by_gene: HashMap<String, Vec<HgvsAnnotation>> = HashMap::new();
        let mut warnings = Vec::new();

        // Normalize chromosome name for lookup
        let chrom = vcf.bare_chrom();

        // Find overlapping transcripts
        let transcripts = self.db.get_by_region(chrom, vcf.pos, vcf.end_pos());

        // Also try with chr prefix if no matches
        let transcripts = if transcripts.is_empty() {
            self.db
                .get_by_region(&vcf.normalized_chrom(), vcf.pos, vcf.end_pos())
        } else {
            transcripts
        };

        if transcripts.is_empty() {
            warnings.push(format!(
                "No overlapping transcripts found for {}:{}-{}",
                chrom,
                vcf.pos,
                vcf.end_pos()
            ));
        }

        // Annotate with each transcript
        for transcript in transcripts {
            // Skip non-coding if configured
            if !self.include_non_coding && !transcript.is_coding() {
                continue;
            }

            let converter = VcfToHgvsConverter::new(transcript);

            match converter.convert(vcf) {
                Ok(mut anns) => {
                    for ann in anns.drain(..) {
                        // Skip intronic if configured
                        if !self.include_intronic && ann.is_intronic {
                            continue;
                        }

                        // Group by gene
                        let gene = ann
                            .gene_symbol
                            .clone()
                            .unwrap_or_else(|| "unknown".to_string());
                        annotations_by_gene
                            .entry(gene)
                            .or_default()
                            .push(ann.clone());

                        annotations.push(ann);
                    }
                }
                Err(e) => {
                    warnings.push(format!("Failed to annotate with {}: {}", transcript.id, e));
                }
            }
        }

        // Sort annotations by gene and then by transcript ID
        annotations.sort_by(|a, b| {
            let gene_a = a.gene_symbol.as_deref().unwrap_or("");
            let gene_b = b.gene_symbol.as_deref().unwrap_or("");
            gene_a.cmp(gene_b).then_with(|| {
                let tx_a = a.transcript_accession.as_deref().unwrap_or("");
                let tx_b = b.transcript_accession.as_deref().unwrap_or("");
                tx_a.cmp(tx_b)
            })
        });

        Ok(MultiIsoformAnnotation {
            vcf: vcf.clone(),
            annotations_by_gene,
            annotations,
            warnings,
        })
    }

    /// Annotate multiple VCF records
    pub fn annotate_batch(
        &self,
        records: &[VcfRecord],
    ) -> Vec<Result<MultiIsoformAnnotation, FerroError>> {
        records.iter().map(|vcf| self.annotate(vcf)).collect()
    }
}

/// Configuration for annotation output
#[derive(Debug, Clone, Default)]
pub struct AnnotationConfig {
    /// Include non-coding transcripts
    pub include_non_coding: bool,
    /// Include intronic variants
    pub include_intronic: bool,
    /// Maximum number of annotations per variant
    pub max_annotations: Option<usize>,
    /// Prefer specific transcript biotypes
    pub preferred_biotypes: Vec<String>,
}

impl AnnotationConfig {
    /// Create a default configuration
    pub fn new() -> Self {
        Self {
            include_non_coding: true,
            include_intronic: true,
            max_annotations: None,
            preferred_biotypes: Vec::new(),
        }
    }

    /// Configure for clinical use (coding transcripts only)
    pub fn clinical() -> Self {
        Self {
            include_non_coding: false,
            include_intronic: false,
            max_annotations: Some(10),
            preferred_biotypes: vec!["protein_coding".to_string()],
        }
    }
}

/// Result of annotating a batch of VCF records
#[derive(Debug)]
pub struct BatchAnnotationResult {
    /// Successfully annotated records
    pub annotated: Vec<MultiIsoformAnnotation>,
    /// Records that failed annotation
    pub failed: Vec<(VcfRecord, FerroError)>,
    /// Total records processed
    pub total: usize,
}

impl BatchAnnotationResult {
    /// Get the success rate
    pub fn success_rate(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            self.annotated.len() as f64 / self.total as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
    use std::sync::OnceLock;

    fn create_test_db() -> TranscriptDb {
        let mut db = TranscriptDb::with_build(GenomeBuild::GRCh38);

        // Add test transcript - MANE Select
        let tx1 = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 2000, 2099),
            ],
            cds_start: Some(50),
            cds_end: Some(150),
            sequence: "ATGC".repeat(50),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        // Add test transcript - not MANE
        let tx2 = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 250, 2000, 2149),
            ],
            cds_start: Some(50),
            cds_end: Some(200),
            sequence: "ATGC".repeat(62),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2149),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);
        db
    }

    #[test]
    fn test_multi_isoform_annotator() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // Should have annotations from both transcripts
        assert_eq!(result.annotation_count(), 2);

        // Both should be for COL1A1
        assert!(result.genes().any(|g| g == "COL1A1"));
    }

    #[test]
    fn test_no_overlapping_transcripts() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        // Position outside any transcript
        let vcf = VcfRecord::snv("chr1", 5000, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        assert_eq!(result.annotation_count(), 0);
        assert!(!result.warnings.is_empty());
    }

    #[test]
    fn test_hgvs_strings() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        let strings = result.hgvs_strings();
        assert_eq!(strings.len(), 2);
    }

    #[test]
    fn test_annotation_config_clinical() {
        let config = AnnotationConfig::clinical();
        assert!(!config.include_non_coding);
        assert!(!config.include_intronic);
        assert!(config.max_annotations.is_some());
    }

    #[test]
    fn test_mane_select_prioritization() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // Should have 2 annotations
        assert_eq!(result.annotation_count(), 2);

        // Primary annotation should be MANE Select (NM_000088.3)
        let primary = result.primary_annotation().unwrap();
        assert!(primary.is_mane_select());
        assert_eq!(
            primary.transcript_accession,
            Some("NM_000088.3".to_string())
        );

        // mane_select_annotation should return the same
        let mane = result.mane_select_annotation().unwrap();
        assert_eq!(mane.transcript_accession, primary.transcript_accession);
    }

    #[test]
    fn test_get_gene() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // Should find COL1A1
        let col1a1_anns = result.get_gene("COL1A1");
        assert!(col1a1_anns.is_some());
        assert_eq!(col1a1_anns.unwrap().len(), 2);

        // Should not find nonexistent gene
        assert!(result.get_gene("NONEXISTENT").is_none());
    }

    #[test]
    fn test_genes_iterator() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        let genes: Vec<_> = result.genes().collect();
        assert_eq!(genes.len(), 1);
        assert!(genes.contains(&&"COL1A1".to_string()));
    }

    #[test]
    fn test_mane_annotations() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        let mane_anns = result.mane_annotations();
        // Only one MANE Select transcript
        assert_eq!(mane_anns.len(), 1);
        assert!(mane_anns[0].is_mane_select());
    }

    #[test]
    fn test_hgvs_by_gene() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        let by_gene = result.hgvs_by_gene();
        assert!(by_gene.contains_key("COL1A1"));
        assert_eq!(by_gene.get("COL1A1").unwrap().len(), 2);
    }

    #[test]
    fn test_annotator_include_non_coding() {
        let db = create_test_db();

        // With non-coding excluded
        let annotator = MultiIsoformAnnotator::new(&db).include_non_coding(false);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // All our test transcripts are coding, so should still have annotations
        assert!(result.annotation_count() > 0);
    }

    #[test]
    fn test_annotator_include_intronic() {
        let db = create_test_db();

        // With intronic excluded
        let annotator = MultiIsoformAnnotator::new(&db).include_intronic(false);

        // Position in intron (between exons)
        let vcf = VcfRecord::snv("chr1", 1500, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // Intronic variants should be filtered out
        for ann in &result.annotations {
            assert!(!ann.is_intronic);
        }
    }

    #[test]
    fn test_annotate_batch() {
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);

        let records = vec![
            VcfRecord::snv("chr1", 1050, 'A', 'G'),
            VcfRecord::snv("chr1", 1060, 'T', 'C'),
            VcfRecord::snv("chr1", 5000, 'A', 'G'), // No overlapping transcripts
        ];

        let results = annotator.annotate_batch(&records);
        assert_eq!(results.len(), 3);

        // First two should have annotations
        assert!(results[0].is_ok());
        assert!(results[0].as_ref().unwrap().annotation_count() > 0);

        assert!(results[1].is_ok());
        assert!(results[1].as_ref().unwrap().annotation_count() > 0);

        // Third should succeed but have no annotations
        assert!(results[2].is_ok());
        assert_eq!(results[2].as_ref().unwrap().annotation_count(), 0);
    }

    #[test]
    fn test_batch_annotation_result_success_rate() {
        let result = BatchAnnotationResult {
            annotated: vec![],
            failed: vec![],
            total: 0,
        };
        assert!((result.success_rate() - 0.0).abs() < 0.001);

        // Create a mock annotation for testing
        let db = create_test_db();
        let annotator = MultiIsoformAnnotator::new(&db);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let ann = annotator.annotate(&vcf).unwrap();

        let result = BatchAnnotationResult {
            annotated: vec![ann],
            failed: vec![],
            total: 2,
        };
        assert!((result.success_rate() - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_annotation_config_new() {
        let config = AnnotationConfig::new();
        assert!(config.include_non_coding);
        assert!(config.include_intronic);
        assert!(config.max_annotations.is_none());
        assert!(config.preferred_biotypes.is_empty());
    }

    #[test]
    fn test_primary_annotation_priority_coding_fallback() {
        // Create a DB with only non-MANE coding transcripts
        let mut db = TranscriptDb::with_build(GenomeBuild::GRCh38);

        let tx1 = Transcript {
            id: "NM_000090.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            exons: vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
            cds_start: Some(10),
            cds_end: Some(90),
            sequence: "ATGC".repeat(25),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000091.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            exons: vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
            cds_start: Some(10),
            cds_end: Some(90),
            sequence: "ATGC".repeat(25),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let annotator = MultiIsoformAnnotator::new(&db);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // Primary should be alphabetically first coding transcript
        let primary = result.primary_annotation().unwrap();
        assert_eq!(
            primary.transcript_accession,
            Some("NM_000090.1".to_string())
        );
    }

    #[test]
    fn test_primary_annotation_empty() {
        let db = TranscriptDb::with_build(GenomeBuild::GRCh38);
        let annotator = MultiIsoformAnnotator::new(&db);

        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // No annotations, primary should be None
        assert!(result.primary_annotation().is_none());
    }

    #[test]
    fn test_mane_select_annotation_none() {
        // Create a DB with only non-MANE transcripts
        let mut db = TranscriptDb::with_build(GenomeBuild::GRCh38);

        let tx = Transcript {
            id: "NM_000090.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            exons: vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
            cds_start: Some(10),
            cds_end: Some(90),
            sequence: "ATGC".repeat(25),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        let annotator = MultiIsoformAnnotator::new(&db);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // No MANE Select
        assert!(result.mane_select_annotation().is_none());
    }

    #[test]
    fn test_primary_annotation_mane_plus_clinical() {
        // Create a DB with MANE Plus Clinical but no MANE Select
        let mut db = TranscriptDb::with_build(GenomeBuild::GRCh38);

        let tx = Transcript {
            id: "NM_000090.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            exons: vec![Exon::with_genomic(1, 1, 100, 1000, 1099)],
            cds_start: Some(10),
            cds_end: Some(90),
            sequence: "ATGC".repeat(25),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        let annotator = MultiIsoformAnnotator::new(&db);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let result = annotator.annotate(&vcf).unwrap();

        // Primary should be MANE Plus Clinical
        let primary = result.primary_annotation().unwrap();
        assert!(primary.is_mane_plus_clinical());
    }
}
