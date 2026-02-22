//! Multi-transcript projection.
//!
//! Project genomic variants onto all overlapping transcripts to generate
//! corresponding c./n. variants.

use crate::data::cdot::CdotMapper;
use crate::data::mapping::{CoordinateMapper, MappingInfo};
use crate::error::FerroError;
use crate::hgvs::location::{CdsPos, GenomePos};
use std::cmp::Ordering;

/// Status of a transcript for clinical prioritization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ManeStatus {
    /// MANE Select - primary clinical transcript
    Select,
    /// MANE Plus Clinical - additional clinical transcript
    PlusClinical,
    /// Not a MANE transcript
    None,
}

impl ManeStatus {
    /// Get priority order (lower is higher priority).
    fn priority(&self) -> u8 {
        match self {
            ManeStatus::Select => 0,
            ManeStatus::PlusClinical => 1,
            ManeStatus::None => 2,
        }
    }
}

/// Projection result for one transcript.
#[derive(Debug, Clone)]
pub struct TranscriptProjection {
    /// Transcript accession (e.g., "NM_000088.3").
    pub transcript_id: String,
    /// Gene symbol if available.
    pub gene_symbol: Option<String>,
    /// CDS position on this transcript.
    pub cds_position: CdsPos,
    /// MANE status.
    pub mane_status: ManeStatus,
    /// Whether this is the canonical transcript.
    pub is_canonical: bool,
    /// Mapping information (exon numbers, region, etc.).
    pub mapping_info: MappingInfo,
    /// CDS length for this transcript.
    pub cds_length: Option<u64>,
}

impl TranscriptProjection {
    /// Format as HGVS c./n. string.
    pub fn to_hgvs_string(&self) -> String {
        let coord_type =
            if self.mapping_info.in_cds || self.mapping_info.in_5utr || self.mapping_info.in_3utr {
                "c"
            } else {
                "n"
            };
        format!(
            "{}:{}.{}",
            self.transcript_id, coord_type, self.cds_position
        )
    }
}

/// All projections for a genomic position.
#[derive(Debug, Clone)]
pub struct ProjectionResult {
    /// Source contig.
    pub source_contig: String,
    /// Source position (1-based).
    pub source_pos: u64,
    /// All transcript projections, sorted by priority.
    pub projections: Vec<TranscriptProjection>,
}

impl ProjectionResult {
    /// Get the highest priority projection (MANE Select if available).
    pub fn best(&self) -> Option<&TranscriptProjection> {
        self.projections.first()
    }

    /// Get all MANE transcripts.
    pub fn mane_only(&self) -> Vec<&TranscriptProjection> {
        self.projections
            .iter()
            .filter(|p| p.mane_status != ManeStatus::None)
            .collect()
    }

    /// Get projections for a specific gene.
    pub fn for_gene(&self, gene: &str) -> Vec<&TranscriptProjection> {
        self.projections
            .iter()
            .filter(|p| p.gene_symbol.as_deref() == Some(gene))
            .collect()
    }

    /// Check if position maps to any coding region.
    pub fn is_coding(&self) -> bool {
        self.projections.iter().any(|p| p.mapping_info.in_cds)
    }
}

/// Projector for mapping genomic positions to transcripts.
#[derive(Debug, Clone)]
pub struct Projector {
    mapper: CoordinateMapper,
    mane_select: std::collections::HashSet<String>,
    mane_plus_clinical: std::collections::HashSet<String>,
    canonical: std::collections::HashSet<String>,
}

impl Projector {
    /// Create a new projector from a cdot mapper.
    pub fn new(cdot: CdotMapper) -> Self {
        Self {
            mapper: CoordinateMapper::new(cdot),
            mane_select: std::collections::HashSet::new(),
            mane_plus_clinical: std::collections::HashSet::new(),
            canonical: std::collections::HashSet::new(),
        }
    }

    /// Create with MANE transcript annotations.
    pub fn with_mane(mut self, select: Vec<String>, plus_clinical: Vec<String>) -> Self {
        self.mane_select = select.into_iter().collect();
        self.mane_plus_clinical = plus_clinical.into_iter().collect();
        self
    }

    /// Create with canonical transcript annotations.
    pub fn with_canonical(mut self, canonical: Vec<String>) -> Self {
        self.canonical = canonical.into_iter().collect();
        self
    }

    /// Get the underlying coordinate mapper.
    pub fn mapper(&self) -> &CoordinateMapper {
        &self.mapper
    }

    /// Project a genomic position to all overlapping transcripts.
    ///
    /// # Arguments
    ///
    /// * `contig` - Contig/chromosome name
    /// * `pos` - 1-based genomic position
    ///
    /// # Returns
    ///
    /// All transcript projections, sorted by clinical priority.
    pub fn project(&self, contig: &str, pos: u64) -> Result<ProjectionResult, FerroError> {
        // Find overlapping transcripts
        let transcripts = self.mapper.find_overlapping_transcripts(contig, pos);

        if transcripts.is_empty() {
            return Ok(ProjectionResult {
                source_contig: contig.to_string(),
                source_pos: pos,
                projections: Vec::new(),
            });
        }

        let genome_pos = GenomePos::new(pos);
        let mut projections = Vec::new();

        for (tx_id, tx) in transcripts {
            match self.mapper.genome_to_cds(tx_id, &genome_pos) {
                Ok(result) => {
                    let projection = TranscriptProjection {
                        transcript_id: tx_id.to_string(),
                        gene_symbol: tx.gene_name.clone(),
                        cds_position: result.variant,
                        mane_status: self.get_mane_status(tx_id),
                        is_canonical: self.canonical.contains(tx_id),
                        mapping_info: result.info,
                        cds_length: tx.cds_length(),
                    };
                    projections.push(projection);
                }
                Err(_) => {
                    // Position might be in intron or outside CDS
                    // Skip this transcript
                    continue;
                }
            }
        }

        // Sort by priority
        self.sort_projections(&mut projections);

        Ok(ProjectionResult {
            source_contig: contig.to_string(),
            source_pos: pos,
            projections,
        })
    }

    /// Project to a specific transcript only.
    pub fn project_to_transcript(
        &self,
        _contig: &str,
        pos: u64,
        transcript_id: &str,
    ) -> Result<TranscriptProjection, FerroError> {
        let tx = self
            .mapper
            .cdot()
            .get_transcript(transcript_id)
            .ok_or_else(|| FerroError::ReferenceNotFound {
                id: transcript_id.to_string(),
            })?;

        let genome_pos = GenomePos::new(pos);
        let result = self.mapper.genome_to_cds(transcript_id, &genome_pos)?;

        Ok(TranscriptProjection {
            transcript_id: transcript_id.to_string(),
            gene_symbol: tx.gene_name.clone(),
            cds_position: result.variant,
            mane_status: self.get_mane_status(transcript_id),
            is_canonical: self.canonical.contains(transcript_id),
            mapping_info: result.info,
            cds_length: tx.cds_length(),
        })
    }

    /// Get MANE status for a transcript.
    fn get_mane_status(&self, transcript_id: &str) -> ManeStatus {
        // Strip version for lookup
        let base_id = transcript_id.split('.').next().unwrap_or(transcript_id);

        if self.mane_select.contains(transcript_id) || self.mane_select.contains(base_id) {
            ManeStatus::Select
        } else if self.mane_plus_clinical.contains(transcript_id)
            || self.mane_plus_clinical.contains(base_id)
        {
            ManeStatus::PlusClinical
        } else {
            ManeStatus::None
        }
    }

    /// Sort projections by clinical priority.
    fn sort_projections(&self, projections: &mut [TranscriptProjection]) {
        projections.sort_by(|a, b| {
            // 1. MANE status (Select > PlusClinical > None)
            let mane_cmp = a.mane_status.priority().cmp(&b.mane_status.priority());
            if mane_cmp != Ordering::Equal {
                return mane_cmp;
            }

            // 2. Canonical status
            let canonical_cmp = b.is_canonical.cmp(&a.is_canonical);
            if canonical_cmp != Ordering::Equal {
                return canonical_cmp;
            }

            // 3. Longer CDS preferred
            let len_cmp = b.cds_length.cmp(&a.cds_length);
            if len_cmp != Ordering::Equal {
                return len_cmp;
            }

            // 4. Alphabetical by transcript ID
            a.transcript_id.cmp(&b.transcript_id)
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::CdotTranscript;
    use crate::reference::Strand;

    fn create_test_mapper() -> CdotMapper {
        let mut cdot = CdotMapper::new();

        // Transcript 1 - "MANE Select"
        let tx1 = CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![
                [1000, 1100, 0, 100],
                [2000, 2200, 100, 300],
                [3000, 3150, 300, 450],
            ],
            cds_start: Some(50),
            cds_end: Some(400),
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        };
        cdot.add_transcript("NM_000001.1".to_string(), tx1);

        // Transcript 2 - different isoform
        let tx2 = CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1150, 0, 150], [2000, 2200, 150, 350]],
            cds_start: Some(50),
            cds_end: Some(300),
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        };
        cdot.add_transcript("NM_000001.2".to_string(), tx2);

        cdot
    }

    #[test]
    fn test_project_single_position() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot);

        let result = projector.project("chr1", 1050).unwrap();

        // Both transcripts should overlap position 1050
        assert!(!result.projections.is_empty());
        assert!(result.projections.len() <= 2);
    }

    #[test]
    fn test_project_with_mane() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot).with_mane(vec!["NM_000001.1".to_string()], vec![]);

        let result = projector.project("chr1", 1050).unwrap();

        // MANE Select should be first
        if let Some(first) = result.best() {
            assert_eq!(first.mane_status, ManeStatus::Select);
        }
    }

    #[test]
    fn test_project_no_overlap() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot);

        let result = projector.project("chr1", 5000).unwrap();

        assert!(result.projections.is_empty());
    }

    #[test]
    fn test_project_to_specific_transcript() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot);

        let result = projector
            .project_to_transcript("chr1", 1050, "NM_000001.1")
            .unwrap();

        assert_eq!(result.transcript_id, "NM_000001.1");
        assert_eq!(result.gene_symbol, Some("TESTGENE".to_string()));
    }

    #[test]
    fn test_mane_status() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot).with_mane(
            vec!["NM_000001.1".to_string()],
            vec!["NM_000001.2".to_string()],
        );

        assert_eq!(projector.get_mane_status("NM_000001.1"), ManeStatus::Select);
        assert_eq!(
            projector.get_mane_status("NM_000001.2"),
            ManeStatus::PlusClinical
        );
        assert_eq!(projector.get_mane_status("NM_000099.1"), ManeStatus::None);
    }

    #[test]
    fn test_projection_result_methods() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot)
            .with_mane(vec!["NM_000001.1".to_string()], vec![])
            .with_canonical(vec!["NM_000001.1".to_string()]);

        let result = projector.project("chr1", 1050).unwrap();

        // Test is_coding
        assert!(result.is_coding());

        // Test for_gene
        let gene_results = result.for_gene("TESTGENE");
        assert!(!gene_results.is_empty());

        // Test mane_only
        let mane_results = result.mane_only();
        assert!(!mane_results.is_empty());
    }

    #[test]
    fn test_to_hgvs_string() {
        let cdot = create_test_mapper();
        let projector = Projector::new(cdot);

        let result = projector
            .project_to_transcript("chr1", 1050, "NM_000001.1")
            .unwrap();

        let hgvs = result.to_hgvs_string();
        assert!(hgvs.starts_with("NM_000001.1:"));
        assert!(hgvs.contains("."));
    }
}
