//! Coordinate mapping between c. and g. coordinates.
//!
//! This module provides functions to convert between CDS (c.) and genomic (g.)
//! coordinates using cdot transcript alignments.

use crate::data::cdot::{CdotMapper, CdotTranscript, CdsPosition};
use crate::error::FerroError;
use crate::hgvs::interval::{CdsInterval, GenomeInterval};
use crate::hgvs::location::{CdsPos, GenomePos};
use crate::reference::Strand;

/// Result of a coordinate mapping operation.
#[derive(Debug, Clone)]
pub struct MappingResult<T> {
    /// The mapped variant.
    pub variant: T,
    /// Information about the mapping.
    pub info: MappingInfo,
}

/// Information about a mapping operation.
#[derive(Debug, Clone, Default)]
pub struct MappingInfo {
    /// The transcript used for mapping.
    pub transcript_id: Option<String>,
    /// Exon number(s) involved.
    pub exon_numbers: Vec<u32>,
    /// Whether the variant is in the CDS.
    pub in_cds: bool,
    /// Whether the variant is in the 5' UTR.
    pub in_5utr: bool,
    /// Whether the variant is in the 3' UTR.
    pub in_3utr: bool,
    /// Whether the variant is intronic.
    pub is_intronic: bool,
    /// Distance to nearest splice site (if intronic).
    pub splice_distance: Option<i64>,
}

/// Coordinate mapper for c. ↔ g. conversions.
#[derive(Debug, Clone)]
pub struct CoordinateMapper {
    cdot: CdotMapper,
}

impl CoordinateMapper {
    /// Create a new coordinate mapper from a cdot mapper.
    pub fn new(cdot: CdotMapper) -> Self {
        Self { cdot }
    }

    /// Get the underlying cdot mapper.
    pub fn cdot(&self) -> &CdotMapper {
        &self.cdot
    }

    /// Convert a CDS position to a genomic position using the cdot mapper's
    /// primary build.
    ///
    /// # Arguments
    ///
    /// * `transcript_id` - The transcript accession (e.g., "NM_000088.3")
    /// * `cds_pos` - The CDS position to convert
    ///
    /// # Returns
    ///
    /// The genomic position and mapping information.
    pub fn cds_to_genome(
        &self,
        transcript_id: &str,
        cds_pos: &CdsPos,
    ) -> Result<MappingResult<GenomePos>, FerroError> {
        self.cds_to_genome_on_build(transcript_id, cds_pos, None)
    }

    /// Build-aware [`cds_to_genome`](Self::cds_to_genome): when `build` is
    /// `Some`, the transcript is resolved via
    /// [`CdotMapper::get_transcript_on_build`] so a multi-build cdot load
    /// returns the alignment for the requested build (rather than silently
    /// using the primary build's view). Used by the projector to honor an
    /// input variant's NG/NC parent (issue #389).
    ///
    /// When `build` is `None`, this is equivalent to
    /// [`cds_to_genome`](Self::cds_to_genome) — exposed so call sites that
    /// computed an `Option<&str>` build hint once can pass it through
    /// without branching on `None`.
    pub fn cds_to_genome_on_build(
        &self,
        transcript_id: &str,
        cds_pos: &CdsPos,
        build: Option<&str>,
    ) -> Result<MappingResult<GenomePos>, FerroError> {
        let tx = match build {
            Some(b) => self.cdot.get_transcript_on_build(transcript_id, b),
            None => self.cdot.get_transcript(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;

        let (genome_pos, info) = self.cds_pos_to_genome_pos(tx, cds_pos, transcript_id)?;

        Ok(MappingResult {
            variant: genome_pos,
            info,
        })
    }

    /// Convert a CDS interval to a genomic interval.
    pub fn cds_interval_to_genome(
        &self,
        transcript_id: &str,
        interval: &CdsInterval,
    ) -> Result<MappingResult<GenomeInterval>, FerroError> {
        let tx = self.cdot.get_transcript(transcript_id).ok_or_else(|| {
            FerroError::ReferenceNotFound {
                id: transcript_id.to_string(),
            }
        })?;

        // Extract positions from UncertainBoundary
        let start_cds = interval
            .start
            .inner()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "CDS interval start position is unknown or a range boundary".to_string(),
            })?;
        let end_cds = interval
            .end
            .inner()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "CDS interval end position is unknown or a range boundary".to_string(),
            })?;

        let (start_genome, start_info) =
            self.cds_pos_to_genome_pos(tx, start_cds, transcript_id)?;
        let (end_genome, end_info) = self.cds_pos_to_genome_pos(tx, end_cds, transcript_id)?;

        // On minus strand, start and end are swapped
        let (start, end) = match tx.strand {
            Strand::Plus => (start_genome, end_genome),
            Strand::Minus => (end_genome, start_genome),
            Strand::Unknown => {
                return Err(FerroError::ConversionError {
                    msg: "strand unknown for transcript".into(),
                })
            }
        };

        let mut info = start_info;
        info.exon_numbers.extend(end_info.exon_numbers);
        info.exon_numbers.sort();
        info.exon_numbers.dedup();

        Ok(MappingResult {
            variant: GenomeInterval::new(start, end),
            info,
        })
    }

    /// Convert a genomic position to a CDS position using the cdot mapper's
    /// primary build.
    ///
    /// # Arguments
    ///
    /// * `transcript_id` - The transcript accession (e.g., "NM_000088.3")
    /// * `genome_pos` - The genomic position to convert
    ///
    /// # Returns
    ///
    /// The CDS position and mapping information.
    pub fn genome_to_cds(
        &self,
        transcript_id: &str,
        genome_pos: &GenomePos,
    ) -> Result<MappingResult<CdsPos>, FerroError> {
        self.genome_to_cds_on_build(transcript_id, genome_pos, None)
    }

    /// Build-aware [`genome_to_cds`](Self::genome_to_cds): when `build` is
    /// `Some`, the transcript is resolved via
    /// [`CdotMapper::get_transcript_on_build`] so a multi-build cdot load
    /// returns the alignment for the requested build. Used by the
    /// projector to honor a g. input's chromosome version (issue #389).
    ///
    /// When `build` is `None`, this is equivalent to
    /// [`genome_to_cds`](Self::genome_to_cds) — exposed so call sites that
    /// computed an `Option<&str>` build hint once can pass it through
    /// without branching on `None`.
    pub fn genome_to_cds_on_build(
        &self,
        transcript_id: &str,
        genome_pos: &GenomePos,
        build: Option<&str>,
    ) -> Result<MappingResult<CdsPos>, FerroError> {
        let tx = match build {
            Some(b) => self.cdot.get_transcript_on_build(transcript_id, b),
            None => self.cdot.get_transcript(transcript_id),
        }
        .ok_or_else(|| FerroError::ReferenceNotFound {
            id: transcript_id.to_string(),
        })?;

        let (cds_pos, info) = self.genome_pos_to_cds_pos(tx, genome_pos, transcript_id)?;

        Ok(MappingResult {
            variant: cds_pos,
            info,
        })
    }

    /// Convert a genomic interval to a CDS interval.
    pub fn genome_interval_to_cds(
        &self,
        transcript_id: &str,
        interval: &GenomeInterval,
    ) -> Result<MappingResult<CdsInterval>, FerroError> {
        let tx = self.cdot.get_transcript(transcript_id).ok_or_else(|| {
            FerroError::ReferenceNotFound {
                id: transcript_id.to_string(),
            }
        })?;

        // Extract positions from UncertainBoundary
        let interval_start =
            interval
                .start
                .inner()
                .ok_or_else(|| FerroError::InvalidCoordinates {
                    msg: "Genomic interval start position is unknown or a range boundary"
                        .to_string(),
                })?;
        let interval_end = interval
            .end
            .inner()
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: "Genomic interval end position is unknown or a range boundary".to_string(),
            })?;

        // On minus strand, start and end are swapped relative to transcript
        let (start, end) = match tx.strand {
            Strand::Plus => (interval_start, interval_end),
            Strand::Minus => (interval_end, interval_start),
            Strand::Unknown => {
                return Err(FerroError::ConversionError {
                    msg: "strand unknown for transcript".into(),
                })
            }
        };

        let (start_cds, start_info) = self.genome_pos_to_cds_pos(tx, start, transcript_id)?;
        let (end_cds, end_info) = self.genome_pos_to_cds_pos(tx, end, transcript_id)?;

        let mut info = start_info;
        info.exon_numbers.extend(end_info.exon_numbers);
        info.exon_numbers.sort();
        info.exon_numbers.dedup();

        Ok(MappingResult {
            variant: CdsInterval::new(start_cds, end_cds),
            info,
        })
    }

    /// Find all transcripts overlapping a genomic position.
    pub fn find_overlapping_transcripts(
        &self,
        contig: &str,
        pos: u64,
    ) -> Vec<(&str, &CdotTranscript)> {
        self.cdot.transcripts_at_position(contig, pos)
    }

    /// Internal: Convert CDS position to genome position.
    fn cds_pos_to_genome_pos(
        &self,
        tx: &CdotTranscript,
        cds_pos: &CdsPos,
        transcript_id: &str,
    ) -> Result<(GenomePos, MappingInfo), FerroError> {
        let mut info = MappingInfo {
            transcript_id: Some(transcript_id.to_string()),
            ..Default::default()
        };

        // Resolve the transcript-coordinate base position, honoring the `utr3`
        // flag on the input. cdot's `cds_to_tx` only accepts an integer and
        // treats every `cds_pos > 0` as a CDS-interior position, which silently
        // maps `c.*N` to `c.N`. For 3'UTR positions we instead compute the tx
        // coord directly from cds_end (cdot's cds_end is 0-based exclusive, so
        // `c.*1` lives at `cds_end`).
        let cds_to_tx_aware = |base: i64| -> Result<u64, FerroError> {
            if cds_pos.utr3 {
                let cds_end = tx.cds_end.ok_or_else(|| FerroError::InvalidCoordinates {
                    msg: format!(
                        "transcript {transcript_id} has no cds_end; cannot resolve c.*{base}"
                    ),
                })?;
                if base < 1 {
                    return Err(FerroError::InvalidCoordinates {
                        msg: format!("Invalid 3' UTR position c.*{base}: must be >= 1"),
                    });
                }
                Ok(cds_end + (base as u64) - 1)
            } else {
                tx.cds_to_tx(base)
                    .ok_or_else(|| FerroError::InvalidCoordinates {
                        msg: format!("Invalid CDS position: {base}"),
                    })
            }
        };

        // Handle intronic positions
        if let Some(offset) = cds_pos.offset {
            // First convert the base position to genomic
            let tx_pos = cds_to_tx_aware(cds_pos.base)?;

            // Find the exon containing this position
            let exon = tx.exon_for_tx_pos(tx_pos);
            let base_genome =
                tx.tx_to_genome(tx_pos)
                    .ok_or_else(|| FerroError::InvalidCoordinates {
                        msg: format!("Cannot map tx position {} to genome", tx_pos),
                    })?;

            // Apply the intronic offset
            let genome_pos = match tx.strand {
                Strand::Plus => {
                    if offset > 0 {
                        // After exon end
                        base_genome + offset as u64
                    } else {
                        // Before exon start (from the perspective of the intron)
                        (base_genome as i64 + offset) as u64
                    }
                }
                Strand::Minus => {
                    if offset > 0 {
                        // After exon end (minus strand = before in genomic coords)
                        base_genome - offset as u64
                    } else {
                        // Before exon start (minus strand = after in genomic coords)
                        (base_genome as i64 - offset) as u64
                    }
                }
                Strand::Unknown => {
                    return Err(FerroError::ConversionError {
                        msg: "strand unknown for transcript".into(),
                    })
                }
            };

            info.is_intronic = true;
            info.splice_distance = Some(offset);
            if let Some(e) = exon {
                info.exon_numbers.push(e.number);
            }

            return Ok((GenomePos::new(genome_pos), info));
        }

        // Non-intronic position
        let tx_pos = cds_to_tx_aware(cds_pos.base)?;

        let genome_pos = tx
            .tx_to_genome(tx_pos)
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: format!("Cannot map tx position {} to genome", tx_pos),
            })?;

        // Determine region
        if let Some(cds_region) = tx.tx_to_cds(tx_pos) {
            match cds_region {
                CdsPosition::FivePrimeUtr(_) => info.in_5utr = true,
                CdsPosition::Cds(_) => info.in_cds = true,
                CdsPosition::ThreePrimeUtr(_) => info.in_3utr = true,
            }
        }

        if let Some(exon) = tx.exon_for_tx_pos(tx_pos) {
            info.exon_numbers.push(exon.number);
        }

        Ok((GenomePos::new(genome_pos), info))
    }

    /// Internal: Convert genome position to CDS position.
    fn genome_pos_to_cds_pos(
        &self,
        tx: &CdotTranscript,
        genome_pos: &GenomePos,
        transcript_id: &str,
    ) -> Result<(CdsPos, MappingInfo), FerroError> {
        // `transcript_id` is part of the parameter list for API stability but
        // no in-tree consumer reads `MappingInfo::transcript_id` — leaving it
        // None avoids a `String::from` allocation per call (~400 calls per
        // SNP project_all on chr17 fan-out).
        let _ = transcript_id;
        let mut info = MappingInfo::default();

        // Check if position is in an exon. `locate_genome_pos` returns both
        // the exon and the converted tx position in a single scan; before
        // c12 we called `genome_to_tx` (one scan) and `exon_for_genome_pos`
        // (a second scan) separately.
        if let Some((tx_pos, exon)) = tx.locate_genome_pos(genome_pos.base) {
            info.exon_numbers.push(exon.number);

            // Convert to CDS position
            let cds_pos = tx
                .tx_to_cds(tx_pos)
                .ok_or_else(|| FerroError::InvalidCoordinates {
                    msg: format!("Cannot convert tx position {} to CDS", tx_pos),
                })?;

            match cds_pos {
                CdsPosition::FivePrimeUtr(offset) => {
                    info.in_5utr = true;
                    Ok((
                        CdsPos {
                            base: -offset,
                            offset: None,
                            utr3: false,
                        },
                        info,
                    ))
                }
                CdsPosition::Cds(pos) => {
                    info.in_cds = true;
                    Ok((
                        CdsPos {
                            base: pos,
                            offset: None,
                            utr3: false,
                        },
                        info,
                    ))
                }
                CdsPosition::ThreePrimeUtr(offset) => {
                    info.in_3utr = true;
                    // Get the last CDS position
                    let cds_end = tx.cds_end.unwrap_or(0);
                    let cds_start = tx.cds_start.unwrap_or(0);
                    let last_cds_pos = (cds_end - cds_start) as i64;
                    Ok((
                        CdsPos {
                            base: last_cds_pos,
                            offset: Some(offset),
                            utr3: true,
                        },
                        info,
                    ))
                }
            }
        } else {
            // Position is intronic - find the nearest exon boundary
            let exons = tx.get_exons();
            let mut best_distance: Option<i64> = None;
            let mut best_exon: Option<usize> = None;
            let mut is_after_exon = false;

            for (i, exon) in exons.iter().enumerate() {
                // Distance to exon end (for positions after exon)
                if genome_pos.base >= exon.genome_end {
                    let dist = (genome_pos.base - exon.genome_end + 1) as i64;
                    if best_distance.is_none() || dist < best_distance.unwrap() {
                        best_distance = Some(dist);
                        best_exon = Some(i);
                        is_after_exon = true;
                    }
                }
                // Distance to exon start (for positions before exon)
                if genome_pos.base < exon.genome_start {
                    let dist = (exon.genome_start - genome_pos.base) as i64;
                    if best_distance.is_none() || dist < best_distance.unwrap() {
                        best_distance = Some(dist);
                        best_exon = Some(i);
                        is_after_exon = false;
                    }
                }
            }

            if let (Some(dist), Some(exon_idx)) = (best_distance, best_exon) {
                let exon = &exons[exon_idx];
                info.is_intronic = true;
                info.exon_numbers.push(exon.number);
                info.splice_distance = Some(dist);

                // Get the CDS position at the exon boundary
                let boundary_tx_pos = if is_after_exon {
                    exon.tx_end - 1 // Last position of exon
                } else {
                    exon.tx_start // First position of exon
                };

                let cds_region = tx.tx_to_cds(boundary_tx_pos);

                // Calculate intronic offset based on strand
                let offset = match tx.strand {
                    Strand::Plus => {
                        if is_after_exon {
                            dist
                        } else {
                            -dist
                        }
                    }
                    Strand::Minus => {
                        if is_after_exon {
                            -dist
                        } else {
                            dist
                        }
                    }
                    Strand::Unknown => {
                        return Err(FerroError::ConversionError {
                            msg: "strand unknown for transcript".into(),
                        })
                    }
                };

                // Convert boundary position to CDS
                if let Some(cds_region) = cds_region {
                    match cds_region {
                        CdsPosition::FivePrimeUtr(utr_offset) => {
                            info.in_5utr = true;
                            Ok((
                                CdsPos {
                                    base: -utr_offset,
                                    offset: Some(offset),
                                    utr3: false,
                                },
                                info,
                            ))
                        }
                        CdsPosition::Cds(pos) => {
                            info.in_cds = false; // In intron, not CDS
                            Ok((
                                CdsPos {
                                    base: pos,
                                    offset: Some(offset),
                                    utr3: false,
                                },
                                info,
                            ))
                        }
                        CdsPosition::ThreePrimeUtr(utr_offset) => {
                            info.in_3utr = true;
                            let cds_end = tx.cds_end.unwrap_or(0);
                            let cds_start = tx.cds_start.unwrap_or(0);
                            let last_cds_pos = (cds_end - cds_start) as i64;
                            Ok((
                                CdsPos {
                                    base: last_cds_pos,
                                    offset: Some(utr_offset + offset),
                                    utr3: true,
                                },
                                info,
                            ))
                        }
                    }
                } else {
                    Err(FerroError::InvalidCoordinates {
                        msg: format!(
                            "Cannot determine CDS position for intronic position {}",
                            genome_pos.base
                        ),
                    })
                }
            } else {
                Err(FerroError::InvalidCoordinates {
                    msg: format!(
                        "Genomic position {} is outside transcript range",
                        genome_pos.base
                    ),
                })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::cdot::CdotTranscript;

    fn create_test_mapper() -> CoordinateMapper {
        let mut cdot = CdotMapper::new();

        // Simple transcript with 3 exons
        let tx = CdotTranscript {
            gene_name: Some("TEST".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![
                [1000, 1100, 0, 100],   // Exon 1: 100bp
                [2000, 2200, 100, 300], // Exon 2: 200bp
                [3000, 3150, 300, 450], // Exon 3: 150bp
            ],
            cds_start: Some(50), // CDS starts at tx pos 50
            cds_end: Some(400),  // CDS ends at tx pos 400
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        };
        cdot.add_transcript("NM_TEST.1".to_string(), tx);

        // Minus strand transcript
        let tx_minus = CdotTranscript {
            gene_name: Some("MINUS".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Minus,
            exons: vec![
                [3000, 3150, 0, 150],   // Exon 1 (3' on genome)
                [2000, 2200, 150, 350], // Exon 2
                [1000, 1100, 350, 450], // Exon 3 (5' on genome)
            ],
            cds_start: Some(50),
            cds_end: Some(400),
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        };
        cdot.add_transcript("NM_MINUS.1".to_string(), tx_minus);

        CoordinateMapper::new(cdot)
    }

    #[test]
    fn test_cds_to_genome_plus_strand() {
        let mapper = create_test_mapper();

        // c.1 (first CDS position) -> genome
        let result = mapper
            .cds_to_genome(
                "NM_TEST.1",
                &CdsPos {
                    base: 1,
                    offset: None,
                    utr3: false,
                },
            )
            .unwrap();

        // c.1 is at tx pos 50, which is at genome pos 1050
        assert_eq!(result.variant.base, 1050);
        assert!(result.info.in_cds);
        assert!(!result.info.is_intronic);
    }

    #[test]
    fn test_cds_to_genome_second_exon() {
        let mapper = create_test_mapper();

        // c.51 (51bp into CDS, which is at tx pos 100, start of exon 2)
        let result = mapper
            .cds_to_genome(
                "NM_TEST.1",
                &CdsPos {
                    base: 51,
                    offset: None,
                    utr3: false,
                },
            )
            .unwrap();

        // c.51 is at tx pos 100, which is at genome pos 2000
        assert_eq!(result.variant.base, 2000);
        assert!(result.info.in_cds);
    }

    #[test]
    fn test_cds_to_genome_5utr() {
        let mapper = create_test_mapper();

        // c.-10 (10bp into 5' UTR)
        let result = mapper
            .cds_to_genome(
                "NM_TEST.1",
                &CdsPos {
                    base: -10,
                    offset: None,
                    utr3: false,
                },
            )
            .unwrap();

        // c.-10 is at tx pos 40, which is at genome pos 1040
        assert_eq!(result.variant.base, 1040);
        assert!(result.info.in_5utr);
    }

    #[test]
    fn test_genome_to_cds_exonic() {
        let mapper = create_test_mapper();

        // Genome pos 1050 -> c.1
        let result = mapper
            .genome_to_cds("NM_TEST.1", &GenomePos::new(1050))
            .unwrap();

        assert_eq!(result.variant.base, 1);
        assert!(result.variant.offset.is_none());
        assert!(result.info.in_cds);
    }

    #[test]
    fn test_genome_to_cds_intronic() {
        let mapper = create_test_mapper();

        // Genome pos 1500 -> intronic (between exon 1 and 2)
        let result = mapper
            .genome_to_cds("NM_TEST.1", &GenomePos::new(1500))
            .unwrap();

        // Should be intronic
        assert!(result.info.is_intronic);
        assert!(result.variant.offset.is_some());
    }

    #[test]
    fn test_genome_to_cds_5utr() {
        let mapper = create_test_mapper();

        // Genome pos 1020 -> 5' UTR (tx pos 20, which is c.-30)
        let result = mapper
            .genome_to_cds("NM_TEST.1", &GenomePos::new(1020))
            .unwrap();

        assert!(result.info.in_5utr);
        assert!(result.variant.base < 0);
    }

    #[test]
    fn test_transcript_not_found() {
        let mapper = create_test_mapper();

        let result = mapper.cds_to_genome(
            "NM_NONEXISTENT.1",
            &CdsPos {
                base: 1,
                offset: None,
                utr3: false,
            },
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_find_overlapping_transcripts() {
        let mapper = create_test_mapper();

        // Position within both transcripts
        let results = mapper.find_overlapping_transcripts("NC_000001.11", 2100);
        assert_eq!(results.len(), 2);

        // Position outside all transcripts
        let results = mapper.find_overlapping_transcripts("NC_000001.11", 5000);
        assert!(results.is_empty());
    }

    /// Boundary semantics of the SuperIntervals index must match the previous
    /// linear-scan implementation: `pos >= min && pos < max` (half-open). Both
    /// fixtures span [1000, 3150) on chr1, so:
    ///   pos = 1000  → in (inclusive lower bound)
    ///   pos = 3149  → in (last position covered by 3150 exclusive)
    ///   pos = 3150  → out (exclusive upper bound)
    ///   pos = 999   → out
    #[test]
    fn test_find_overlapping_transcripts_boundary_semantics() {
        let mapper = create_test_mapper();
        assert_eq!(
            mapper
                .find_overlapping_transcripts("NC_000001.11", 1000)
                .len(),
            2,
            "inclusive lower bound"
        );
        assert_eq!(
            mapper
                .find_overlapping_transcripts("NC_000001.11", 3149)
                .len(),
            2,
            "inclusive upper bound (3150 - 1)"
        );
        assert!(
            mapper
                .find_overlapping_transcripts("NC_000001.11", 3150)
                .is_empty(),
            "exclusive upper bound"
        );
        assert!(
            mapper
                .find_overlapping_transcripts("NC_000001.11", 999)
                .is_empty(),
            "below lower bound"
        );
    }

    /// The query uses each transcript's `[min_exon_start, max_exon_end)` span,
    /// so intronic positions inside that span (between two exons of the same
    /// transcript) still return the transcript. The downstream
    /// `genome_to_cds` is what tells the caller it landed in an intron.
    #[test]
    fn test_find_overlapping_transcripts_returns_intronic_spans() {
        let mapper = create_test_mapper();
        // Position 1500 is in the intron between exon 1 [1000,1100) and exon 2
        // [2000,2200). It's inside the transcript span [1000, 3150).
        let results = mapper.find_overlapping_transcripts("NC_000001.11", 1500);
        assert_eq!(results.len(), 2, "intronic position must return both txs");
    }

    /// After `add_transcript` is called following a prior query, the new
    /// transcript must show up in subsequent queries — i.e. the OnceCell
    /// invalidation works.
    #[test]
    fn test_find_overlapping_transcripts_invalidates_after_add() {
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_FIRST.1".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_000001.11".to_string(),
                strand: Strand::Plus,
                exons: vec![[1000, 1100, 0, 100]],
                cds_start: Some(0),
                cds_end: Some(100),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let mapper = CoordinateMapper::new(cdot);
        assert_eq!(
            mapper
                .find_overlapping_transcripts("NC_000001.11", 1050)
                .len(),
            1,
        );

        // Add a second transcript via a fresh mapper (mimicking the typical
        // build-once flow but exercising the invalidation explicitly).
        let mut cdot2 = mapper.cdot().clone();
        cdot2.add_transcript(
            "NM_SECOND.1".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_000001.11".to_string(),
                strand: Strand::Plus,
                exons: vec![[1000, 1100, 0, 100]],
                cds_start: Some(0),
                cds_end: Some(100),
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        let mapper2 = CoordinateMapper::new(cdot2);
        assert_eq!(
            mapper2
                .find_overlapping_transcripts("NC_000001.11", 1050)
                .len(),
            2,
            "second transcript must be visible after add",
        );
    }
}
