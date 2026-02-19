//! VCF to HGVS conversion
//!
//! This module provides functionality to convert VCF records to HGVS notation.

use std::str::FromStr;

use crate::convert::mapper::CoordinateMapper;
use crate::convert::noncoding::IntronicConsequence;
use crate::error::FerroError;
use crate::hgvs::edit::{Base, InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::{CdsInterval, GenomeInterval, TxInterval};
use crate::hgvs::location::GenomePos;
use crate::hgvs::variant::{Accession, CdsVariant, GenomeVariant, HgvsVariant, LocEdit, TxVariant};
use crate::reference::transcript::{IntronPosition, ManeStatus, Transcript};

use super::record::VcfRecord;

/// Result of converting a VCF record to HGVS
#[derive(Debug, Clone)]
pub struct HgvsAnnotation {
    /// The HGVS variant representation
    pub variant: HgvsVariant,
    /// Gene symbol if known
    pub gene_symbol: Option<String>,
    /// Transcript accession
    pub transcript_accession: Option<String>,
    /// Whether this is from a coding transcript
    pub is_coding: bool,
    /// Whether this variant is intronic
    pub is_intronic: bool,
    /// MANE (Matched Annotation from NCBI and EBI) status
    pub mane_status: ManeStatus,
    /// Intron number if the variant is intronic (1-based)
    pub intron_number: Option<u32>,
    /// Intronic position details if applicable
    pub intron_position: Option<IntronPosition>,
    /// Predicted intronic consequence if applicable
    pub intronic_consequence: Option<IntronicConsequence>,
}

impl HgvsAnnotation {
    /// Get the HGVS string representation
    pub fn hgvs_string(&self) -> String {
        self.variant.to_string()
    }

    /// Check if this annotation is from a MANE Select transcript
    pub fn is_mane_select(&self) -> bool {
        self.mane_status.is_select()
    }

    /// Check if this annotation is from a MANE Plus Clinical transcript
    pub fn is_mane_plus_clinical(&self) -> bool {
        self.mane_status.is_plus_clinical()
    }

    /// Check if this annotation is from any MANE transcript
    pub fn is_mane(&self) -> bool {
        self.mane_status.is_mane()
    }

    /// Check if this is a deep intronic variant (>50bp from exon)
    pub fn is_deep_intronic(&self) -> bool {
        self.intron_position
            .as_ref()
            .map(|p| p.is_deep_intronic())
            .unwrap_or(false)
    }

    /// Check if this variant affects a canonical splice site
    pub fn affects_splice_site(&self) -> bool {
        self.intronic_consequence
            .as_ref()
            .map(|c| c.affects_canonical_splice_site())
            .unwrap_or(false)
    }

    /// Get the impact level of the intronic consequence
    pub fn intronic_impact(&self) -> Option<&'static str> {
        self.intronic_consequence.as_ref().map(|c| c.impact())
    }

    /// Get the SO term for the intronic consequence
    pub fn intronic_so_term(&self) -> Option<&'static str> {
        self.intronic_consequence.as_ref().map(|c| c.so_term())
    }

    /// Get distance from nearest exon if intronic
    pub fn distance_from_exon(&self) -> Option<u64> {
        self.intron_position
            .as_ref()
            .map(|p| p.offset.unsigned_abs())
    }
}

/// Converter for VCF records to HGVS notation
pub struct VcfToHgvsConverter<'a> {
    transcript: &'a Transcript,
    mapper: CoordinateMapper<'a>,
}

impl<'a> VcfToHgvsConverter<'a> {
    /// Create a new converter for a specific transcript
    pub fn new(transcript: &'a Transcript) -> Self {
        let mapper = CoordinateMapper::new(transcript);
        Self { transcript, mapper }
    }

    /// Convert a VCF record to HGVS notation
    ///
    /// Returns an annotation for each alternate allele in the VCF record.
    pub fn convert(&self, vcf: &VcfRecord) -> Result<Vec<HgvsAnnotation>, FerroError> {
        let mut annotations = Vec::new();

        for alt in &vcf.alternate {
            let annotation = self.convert_single_allele(vcf, alt)?;
            annotations.push(annotation);
        }

        Ok(annotations)
    }

    /// Convert a single alternate allele to HGVS
    fn convert_single_allele(
        &self,
        vcf: &VcfRecord,
        alt: &str,
    ) -> Result<HgvsAnnotation, FerroError> {
        // Determine the edit type from REF/ALT
        let (edit, start_pos, end_pos) = self.determine_edit(&vcf.reference, alt)?;

        // Calculate genomic positions
        let genomic_start = vcf.pos + start_pos;
        let genomic_end = vcf.pos + end_pos;

        // Try to map to transcript coordinates
        if self.transcript.has_genomic_coords() {
            self.convert_with_transcript(genomic_start, genomic_end, edit)
        } else {
            // Fall back to genomic-only representation
            self.convert_genomic_only(vcf, genomic_start, genomic_end, edit)
        }
    }

    /// Convert using transcript coordinates (c. or n.)
    fn convert_with_transcript(
        &self,
        genomic_start: u64,
        genomic_end: u64,
        edit: NaEdit,
    ) -> Result<HgvsAnnotation, FerroError> {
        // Map to transcript coordinates (with intronic offset support)
        let tx_start_opt = self.mapper.genomic_to_tx(genomic_start)?;
        let tx_end_opt = if genomic_start == genomic_end {
            tx_start_opt
        } else {
            self.mapper.genomic_to_tx(genomic_end)?
        };

        // Check if the variant is intronic
        let is_intronic = tx_start_opt.is_none() || tx_end_opt.is_none();

        // Get intronic position info if applicable
        let intron_position = if is_intronic {
            self.mapper.get_intron_position(genomic_start)
        } else {
            None
        };

        // If we have a coding transcript, use CDS coordinates
        if self.transcript.is_coding() {
            self.create_cds_annotation(
                genomic_start,
                genomic_end,
                edit,
                is_intronic,
                intron_position,
            )
        } else {
            self.create_tx_annotation_with_intron(
                genomic_start,
                genomic_end,
                edit,
                is_intronic,
                intron_position,
            )
        }
    }

    /// Create a CDS annotation (c. notation)
    fn create_cds_annotation(
        &self,
        genomic_start: u64,
        genomic_end: u64,
        edit: NaEdit,
        is_intronic: bool,
        intron_position: Option<IntronPosition>,
    ) -> Result<HgvsAnnotation, FerroError> {
        // Map to CDS coordinates (with intronic offset support for intronic positions)
        let (cds_start, cds_end) = if is_intronic && intron_position.is_some() {
            // Use intronic mapping to get CDS position with offset
            let start = self.mapper.genomic_to_cds_with_intron(genomic_start)?;
            let end = if genomic_start == genomic_end {
                start
            } else {
                self.mapper.genomic_to_cds_with_intron(genomic_end)?
            };
            (Some(start), Some(end))
        } else {
            // Standard mapping for exonic positions
            let start = self.mapper.genomic_to_cds(genomic_start)?;
            let end = if genomic_start == genomic_end {
                start
            } else {
                self.mapper.genomic_to_cds(genomic_end)?
            };
            (start, end)
        };

        // Create interval based on whether positions are available
        let interval = match (&cds_start, &cds_end) {
            (Some(start), Some(end)) => {
                if start == end {
                    CdsInterval::point(*start)
                } else {
                    CdsInterval::new(*start, *end)
                }
            }
            (Some(start), None) | (None, Some(start)) => CdsInterval::point(*start),
            (None, None) => {
                // Should not happen with intronic support, but handle gracefully
                return Err(FerroError::ConversionError {
                    msg: "Cannot map position to CDS coordinates".to_string(),
                });
            }
        };

        let accession = self.build_accession();
        let variant = HgvsVariant::Cds(CdsVariant {
            accession,
            gene_symbol: self.transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        });

        // Calculate intronic consequence if applicable
        let intronic_consequence = cds_start.as_ref().and_then(|p| {
            if p.is_intronic() {
                IntronicConsequence::from_cds_pos(p)
            } else {
                None
            }
        });

        Ok(HgvsAnnotation {
            variant,
            gene_symbol: self.transcript.gene_symbol.clone(),
            transcript_accession: Some(self.transcript.id.clone()),
            is_coding: true,
            is_intronic,
            mane_status: self.transcript.mane_status,
            intron_number: intron_position.as_ref().map(|p| p.intron_number),
            intron_position,
            intronic_consequence,
        })
    }

    /// Create a non-coding transcript annotation with intronic support (n. notation)
    fn create_tx_annotation_with_intron(
        &self,
        genomic_start: u64,
        genomic_end: u64,
        edit: NaEdit,
        is_intronic: bool,
        intron_position: Option<IntronPosition>,
    ) -> Result<HgvsAnnotation, FerroError> {
        // Map to transcript coordinates (with intronic offset support for intronic positions)
        let (tx_start, tx_end) = if is_intronic && intron_position.is_some() {
            // Use intronic mapping to get transcript position with offset
            let start = self.mapper.genomic_to_tx_with_intron(genomic_start)?;
            let end = if genomic_start == genomic_end {
                start
            } else {
                self.mapper.genomic_to_tx_with_intron(genomic_end)?
            };
            (Some(start), Some(end))
        } else {
            // Standard mapping for exonic positions
            let start = self.mapper.genomic_to_tx(genomic_start)?;
            let end = if genomic_start == genomic_end {
                start
            } else {
                self.mapper.genomic_to_tx(genomic_end)?
            };
            (start, end)
        };

        let interval = match (&tx_start, &tx_end) {
            (Some(start), Some(end)) => {
                if start == end {
                    TxInterval::point(*start)
                } else {
                    TxInterval::new(*start, *end)
                }
            }
            (Some(start), None) | (None, Some(start)) => TxInterval::point(*start),
            (None, None) => {
                return Err(FerroError::ConversionError {
                    msg: "Cannot map position to transcript coordinates".to_string(),
                });
            }
        };

        let accession = self.build_accession();
        let variant = HgvsVariant::Tx(TxVariant {
            accession,
            gene_symbol: self.transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        });

        // Calculate intronic consequence if applicable
        let intronic_consequence = tx_start.as_ref().and_then(|p| {
            if p.is_intronic() {
                IntronicConsequence::from_tx_pos(p)
            } else {
                None
            }
        });

        Ok(HgvsAnnotation {
            variant,
            gene_symbol: self.transcript.gene_symbol.clone(),
            transcript_accession: Some(self.transcript.id.clone()),
            is_coding: false,
            is_intronic,
            mane_status: self.transcript.mane_status,
            intron_number: intron_position.as_ref().map(|p| p.intron_number),
            intron_position,
            intronic_consequence,
        })
    }

    /// Convert to genomic-only representation (g. notation)
    fn convert_genomic_only(
        &self,
        vcf: &VcfRecord,
        genomic_start: u64,
        genomic_end: u64,
        edit: NaEdit,
    ) -> Result<HgvsAnnotation, FerroError> {
        let interval = if genomic_start == genomic_end {
            GenomeInterval::point(GenomePos::new(genomic_start))
        } else {
            GenomeInterval::new(GenomePos::new(genomic_start), GenomePos::new(genomic_end))
        };

        // Build a genomic accession from chromosome
        let accession = Accession::new("NC", vcf.normalized_chrom(), None);

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession,
            gene_symbol: self.transcript.gene_symbol.clone(),
            loc_edit: LocEdit::new(interval, edit),
        });

        Ok(HgvsAnnotation {
            variant,
            gene_symbol: self.transcript.gene_symbol.clone(),
            transcript_accession: Some(self.transcript.id.clone()),
            is_coding: false,
            is_intronic: false,
            mane_status: self.transcript.mane_status,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        })
    }

    /// Determine the edit type from REF and ALT alleles
    ///
    /// Returns the edit type and the relative position offsets for start and end
    fn determine_edit(&self, reference: &str, alt: &str) -> Result<(NaEdit, u64, u64), FerroError> {
        // Trim common prefix and suffix
        let (ref_trimmed, alt_trimmed, prefix_len, _suffix_len) = trim_common_bases(reference, alt);

        let ref_len = ref_trimmed.len();
        let alt_len = alt_trimmed.len();

        let edit = if ref_len == 1 && alt_len == 1 {
            // SNV
            let ref_base =
                Base::from_char(ref_trimmed.chars().next().unwrap()).ok_or_else(|| {
                    FerroError::ConversionError {
                        msg: format!("Invalid reference base: {}", ref_trimmed),
                    }
                })?;
            let alt_base =
                Base::from_char(alt_trimmed.chars().next().unwrap()).ok_or_else(|| {
                    FerroError::ConversionError {
                        msg: format!("Invalid alternate base: {}", alt_trimmed),
                    }
                })?;
            NaEdit::Substitution {
                reference: ref_base,
                alternative: alt_base,
            }
        } else if ref_len > 0 && alt_len == 0 {
            // Deletion
            let seq = if ref_len <= 10 {
                Some(
                    Sequence::from_str(ref_trimmed).map_err(|_| FerroError::ConversionError {
                        msg: format!("Invalid sequence: {}", ref_trimmed),
                    })?,
                )
            } else {
                None
            };
            NaEdit::Deletion {
                sequence: seq,
                length: None,
            }
        } else if ref_len == 0 && alt_len > 0 {
            // Insertion
            let seq = Sequence::from_str(alt_trimmed).map_err(|_| FerroError::ConversionError {
                msg: format!("Invalid sequence: {}", alt_trimmed),
            })?;
            NaEdit::Insertion {
                sequence: InsertedSequence::Literal(seq),
            }
        } else {
            // Delins (deletion-insertion)
            let seq = Sequence::from_str(alt_trimmed).map_err(|_| FerroError::ConversionError {
                msg: format!("Invalid sequence: {}", alt_trimmed),
            })?;
            NaEdit::Delins {
                sequence: InsertedSequence::Literal(seq),
            }
        };

        // Calculate position offsets
        let start_offset = prefix_len as u64;
        let end_offset = if ref_len == 0 {
            // Insertion: position between two bases
            start_offset
        } else {
            start_offset + ref_len as u64 - 1
        };

        Ok((edit, start_offset, end_offset))
    }

    /// Build an accession from the transcript ID
    fn build_accession(&self) -> Accession {
        // Parse transcript ID to extract prefix, number, and version
        // Format: PREFIX_NUMBER.VERSION (e.g., NM_000088.3) or ENST00000000.VERSION
        let id = &self.transcript.id;

        if let Some((base, version)) = id.rsplit_once('.') {
            let version_num: Option<u32> = version.parse().ok();
            if let Some((prefix, number)) = base.split_once('_') {
                Accession::new(prefix, number, version_num)
            } else if base.starts_with("ENST") || base.starts_with("ENSG") {
                // Ensembl format: ENST00000000000
                let prefix = &base[..4];
                let number = &base[4..];
                Accession::with_style(prefix, number, version_num, true)
            } else {
                Accession::new(base, "", version_num)
            }
        } else if let Some((prefix, number)) = id.split_once('_') {
            Accession::new(prefix, number, None)
        } else {
            Accession::new(id.as_str(), "", None)
        }
    }
}

/// Trim common prefix and suffix from REF and ALT
///
/// Returns (trimmed_ref, trimmed_alt, prefix_len, suffix_len)
fn trim_common_bases<'a>(reference: &'a str, alt: &'a str) -> (&'a str, &'a str, usize, usize) {
    let ref_bytes = reference.as_bytes();
    let alt_bytes = alt.as_bytes();

    // Find common prefix
    let mut prefix_len = 0;
    while prefix_len < ref_bytes.len()
        && prefix_len < alt_bytes.len()
        && ref_bytes[prefix_len] == alt_bytes[prefix_len]
    {
        prefix_len += 1;
    }

    // Find common suffix (but leave at least one base if needed)
    let ref_remaining = &ref_bytes[prefix_len..];
    let alt_remaining = &alt_bytes[prefix_len..];

    let mut suffix_len = 0;
    while suffix_len < ref_remaining.len()
        && suffix_len < alt_remaining.len()
        && ref_remaining[ref_remaining.len() - 1 - suffix_len]
            == alt_remaining[alt_remaining.len() - 1 - suffix_len]
    {
        suffix_len += 1;
    }

    let ref_trimmed = &reference[prefix_len..reference.len() - suffix_len];
    let alt_trimmed = &alt[prefix_len..alt.len() - suffix_len];

    (ref_trimmed, alt_trimmed, prefix_len, suffix_len)
}

/// Convert a VCF record to a genomic HGVS variant (g. notation)
///
/// This is a simple conversion that doesn't require transcript information.
pub fn vcf_to_genomic_hgvs(vcf: &VcfRecord, alt_index: usize) -> Result<GenomeVariant, FerroError> {
    let alt = vcf
        .alternate
        .get(alt_index)
        .ok_or_else(|| FerroError::ConversionError {
            msg: format!(
                "Alt index {} out of bounds (have {} alts)",
                alt_index,
                vcf.alternate.len()
            ),
        })?;

    let (ref_trimmed, alt_trimmed, prefix_len, _suffix_len) =
        trim_common_bases(&vcf.reference, alt);

    let ref_len = ref_trimmed.len();
    let alt_len = alt_trimmed.len();

    let edit = if ref_len == 1 && alt_len == 1 {
        let ref_base = Base::from_char(ref_trimmed.chars().next().unwrap()).ok_or_else(|| {
            FerroError::ConversionError {
                msg: format!("Invalid reference base: {}", ref_trimmed),
            }
        })?;
        let alt_base = Base::from_char(alt_trimmed.chars().next().unwrap()).ok_or_else(|| {
            FerroError::ConversionError {
                msg: format!("Invalid alternate base: {}", alt_trimmed),
            }
        })?;
        NaEdit::Substitution {
            reference: ref_base,
            alternative: alt_base,
        }
    } else if ref_len > 0 && alt_len == 0 {
        let seq = if ref_len <= 10 {
            Some(
                Sequence::from_str(ref_trimmed).map_err(|_| FerroError::ConversionError {
                    msg: format!("Invalid sequence: {}", ref_trimmed),
                })?,
            )
        } else {
            None
        };
        NaEdit::Deletion {
            sequence: seq,
            length: None,
        }
    } else if ref_len == 0 && alt_len > 0 {
        let seq = Sequence::from_str(alt_trimmed).map_err(|_| FerroError::ConversionError {
            msg: format!("Invalid sequence: {}", alt_trimmed),
        })?;
        NaEdit::Insertion {
            sequence: InsertedSequence::Literal(seq),
        }
    } else {
        let seq = Sequence::from_str(alt_trimmed).map_err(|_| FerroError::ConversionError {
            msg: format!("Invalid sequence: {}", alt_trimmed),
        })?;
        NaEdit::Delins {
            sequence: InsertedSequence::Literal(seq),
        }
    };

    // Build genomic interval
    let start_pos = vcf.pos + prefix_len as u64;
    let end_pos = if ref_len == 0 {
        start_pos
    } else {
        start_pos + ref_len as u64 - 1
    };

    let interval = if start_pos == end_pos {
        GenomeInterval::point(GenomePos::new(start_pos))
    } else {
        GenomeInterval::new(GenomePos::new(start_pos), GenomePos::new(end_pos))
    };

    // Build accession from chromosome
    // Use bare chrom for now - a more complete implementation would
    // map to RefSeq assembly accessions (NC_000001.11 for chr1, etc.)
    let accession = Accession::new("NC", vcf.bare_chrom(), None);

    Ok(GenomeVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, GenomeBuild, Strand};
    use std::sync::OnceLock;

    fn create_test_transcript() -> Transcript {
        Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 2000, 2099),
                Exon::with_genomic(3, 201, 300, 3000, 3099),
            ],
            cds_start: Some(50),
            cds_end: Some(250),
            sequence: "ATGC".repeat(75),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(3099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_trim_common_bases() {
        // SNV
        let (ref_t, alt_t, pre, suf) = trim_common_bases("A", "G");
        assert_eq!(ref_t, "A");
        assert_eq!(alt_t, "G");
        assert_eq!(pre, 0);
        assert_eq!(suf, 0);

        // Deletion with anchor
        let (ref_t, alt_t, pre, _suf) = trim_common_bases("ATG", "A");
        assert_eq!(ref_t, "TG");
        assert_eq!(alt_t, "");
        assert_eq!(pre, 1);

        // Insertion with anchor
        let (ref_t, alt_t, pre, _suf) = trim_common_bases("A", "ATG");
        assert_eq!(ref_t, "");
        assert_eq!(alt_t, "TG");
        assert_eq!(pre, 1);
    }

    #[test]
    fn test_vcf_to_genomic_snv() {
        let vcf = VcfRecord::snv("chr1", 12345, 'A', 'G');
        let variant = vcf_to_genomic_hgvs(&vcf, 0).unwrap();
        assert_eq!(variant.to_string(), "NC_1:g.12345A>G");
    }

    #[test]
    fn test_vcf_to_genomic_deletion() {
        let vcf = VcfRecord::deletion("chr1", 12345, "ATG");
        let variant = vcf_to_genomic_hgvs(&vcf, 0).unwrap();
        assert_eq!(variant.to_string(), "NC_1:g.12346_12347delTG");
    }

    #[test]
    fn test_vcf_to_genomic_insertion() {
        let vcf = VcfRecord::insertion("chr1", 12345, 'A', "TG");
        let variant = vcf_to_genomic_hgvs(&vcf, 0).unwrap();
        assert_eq!(variant.to_string(), "NC_1:g.12346insTG");
    }

    #[test]
    fn test_converter_snv() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        // SNV at position 1050 (within exon 1 at tx pos 51, CDS pos 2)
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let annotations = converter.convert(&vcf).unwrap();

        assert_eq!(annotations.len(), 1);
        let ann = &annotations[0];
        assert!(ann.is_coding);
        assert!(!ann.is_intronic);
        assert_eq!(ann.gene_symbol, Some("COL1A1".to_string()));
    }

    #[test]
    fn test_determine_edit_snv() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        let (edit, start, end) = converter.determine_edit("A", "G").unwrap();
        assert!(matches!(edit, NaEdit::Substitution { .. }));
        assert_eq!(start, 0);
        assert_eq!(end, 0);
    }

    #[test]
    fn test_determine_edit_deletion() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        let (edit, start, end) = converter.determine_edit("ATG", "A").unwrap();
        assert!(matches!(edit, NaEdit::Deletion { .. }));
        assert_eq!(start, 1);
        assert_eq!(end, 2);
    }

    #[test]
    fn test_determine_edit_insertion() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        let (edit, start, end) = converter.determine_edit("A", "ATG").unwrap();
        assert!(matches!(edit, NaEdit::Insertion { .. }));
        assert_eq!(start, 1);
        assert_eq!(end, 1);
    }

    #[test]
    fn test_determine_edit_delins() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        let (edit, _start, _end) = converter.determine_edit("ATG", "CCC").unwrap();
        assert!(matches!(edit, NaEdit::Delins { .. }));
    }

    #[test]
    fn test_build_accession_refseq() {
        let mut transcript = create_test_transcript();
        transcript.id = "NM_000088.3".to_string();
        let converter = VcfToHgvsConverter::new(&transcript);
        let acc = converter.build_accession();
        assert_eq!(&*acc.prefix, "NM");
        assert_eq!(&*acc.number, "000088");
        assert_eq!(acc.version, Some(3));
    }

    #[test]
    fn test_build_accession_ensembl() {
        let mut transcript = create_test_transcript();
        transcript.id = "ENST00000000123.4".to_string();
        let converter = VcfToHgvsConverter::new(&transcript);
        let acc = converter.build_accession();
        assert_eq!(&*acc.prefix, "ENST");
        assert_eq!(&*acc.number, "00000000123");
        assert_eq!(acc.version, Some(4));
        assert!(acc.ensembl_style);
    }

    #[test]
    fn test_build_accession_no_version() {
        let mut transcript = create_test_transcript();
        transcript.id = "NM_000088".to_string();
        let converter = VcfToHgvsConverter::new(&transcript);
        let acc = converter.build_accession();
        assert_eq!(&*acc.prefix, "NM");
        assert_eq!(&*acc.number, "000088");
        assert_eq!(acc.version, None);
    }

    #[test]
    fn test_build_accession_plain() {
        let mut transcript = create_test_transcript();
        transcript.id = "TRANSCRIPT123".to_string();
        let converter = VcfToHgvsConverter::new(&transcript);
        let acc = converter.build_accession();
        assert_eq!(&*acc.prefix, "TRANSCRIPT123");
        assert_eq!(&*acc.number, "");
    }

    #[test]
    fn test_trim_common_bases_delins() {
        let (ref_t, alt_t, pre, _suf) = trim_common_bases("ATGC", "AGGG");
        assert_eq!(pre, 1); // "A" is common prefix
        assert_eq!(ref_t, "TGC");
        assert_eq!(alt_t, "GGG");
    }

    #[test]
    fn test_trim_common_bases_suffix() {
        let (ref_t, alt_t, pre, suf) = trim_common_bases("ATGC", "AGGC");
        // Common prefix "A", common suffix "GC" (2 chars)
        assert_eq!(pre, 1);
        assert_eq!(suf, 2);
        assert_eq!(ref_t, "T");
        assert_eq!(alt_t, "G");
    }

    #[test]
    fn test_trim_common_bases_identical() {
        let (ref_t, alt_t, pre, suf) = trim_common_bases("ATGC", "ATGC");
        assert_eq!(ref_t, "");
        assert_eq!(alt_t, "");
        assert_eq!(pre, 4);
        assert_eq!(suf, 0); // Suffix is found after prefix is removed
    }

    #[test]
    fn test_trim_common_bases_empty() {
        let (ref_t, alt_t, pre, suf) = trim_common_bases("", "");
        assert_eq!(ref_t, "");
        assert_eq!(alt_t, "");
        assert_eq!(pre, 0);
        assert_eq!(suf, 0);
    }

    #[test]
    fn test_hgvs_annotation_methods() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NM_000088.3:c.459A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::Select,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        assert_eq!(ann.hgvs_string(), "NM_000088.3:c.459A>G");
        assert!(ann.is_mane_select());
        assert!(!ann.is_mane_plus_clinical());
        assert!(ann.is_mane());
        assert!(!ann.is_deep_intronic());
        assert!(!ann.affects_splice_site());
        assert!(ann.intronic_impact().is_none());
        assert!(ann.intronic_so_term().is_none());
        assert!(ann.distance_from_exon().is_none());
    }

    #[test]
    fn test_hgvs_annotation_mane_plus_clinical() {
        use crate::hgvs::parser::parse_hgvs;

        let variant = parse_hgvs("NM_000088.3:c.459A>G").unwrap();
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: false,
            mane_status: ManeStatus::PlusClinical,
            intron_number: None,
            intron_position: None,
            intronic_consequence: None,
        };

        assert!(!ann.is_mane_select());
        assert!(ann.is_mane_plus_clinical());
        assert!(ann.is_mane());
    }

    #[test]
    fn test_hgvs_annotation_intronic() {
        use crate::convert::noncoding::IntronicConsequence;
        use crate::hgvs::parser::parse_hgvs;
        use crate::reference::transcript::IntronBoundary;

        let variant = parse_hgvs("NM_000088.3:c.100+50A>G").unwrap();
        let intron_position = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 50,
            tx_boundary_pos: 100,
            intron_length: 1000,
        };
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: true,
            mane_status: ManeStatus::Select,
            intron_number: Some(1),
            intron_position: Some(intron_position),
            intronic_consequence: Some(IntronicConsequence::NearSpliceSiteVariant),
        };

        assert!(!ann.is_deep_intronic());
        assert!(!ann.affects_splice_site());
        assert_eq!(ann.intronic_impact(), Some("MODIFIER"));
        assert_eq!(ann.intronic_so_term(), Some("intron_variant"));
        assert_eq!(ann.distance_from_exon(), Some(50));
    }

    #[test]
    fn test_hgvs_annotation_deep_intronic() {
        use crate::convert::noncoding::IntronicConsequence;
        use crate::hgvs::parser::parse_hgvs;
        use crate::reference::transcript::IntronBoundary;

        let variant = parse_hgvs("NM_000088.3:c.100+500A>G").unwrap();
        let intron_position = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 500,
            tx_boundary_pos: 100,
            intron_length: 1000,
        };
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: true,
            mane_status: ManeStatus::Select,
            intron_number: Some(1),
            intron_position: Some(intron_position),
            intronic_consequence: Some(IntronicConsequence::IntronVariant),
        };

        assert!(ann.is_deep_intronic());
        assert!(!ann.affects_splice_site());
        assert_eq!(ann.intronic_impact(), Some("MODIFIER"));
    }

    #[test]
    fn test_hgvs_annotation_affects_splice_site() {
        use crate::convert::noncoding::IntronicConsequence;
        use crate::hgvs::parser::parse_hgvs;
        use crate::reference::transcript::IntronBoundary;

        let variant = parse_hgvs("NM_000088.3:c.100+1A>G").unwrap();
        let intron_position = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 1,
            tx_boundary_pos: 100,
            intron_length: 1000,
        };
        let ann = HgvsAnnotation {
            variant,
            gene_symbol: Some("COL1A1".to_string()),
            transcript_accession: Some("NM_000088.3".to_string()),
            is_coding: true,
            is_intronic: true,
            mane_status: ManeStatus::Select,
            intron_number: Some(1),
            intron_position: Some(intron_position),
            intronic_consequence: Some(IntronicConsequence::SpliceDonorVariant),
        };

        assert!(!ann.is_deep_intronic());
        assert!(ann.affects_splice_site());
        assert_eq!(ann.intronic_impact(), Some("HIGH"));
        assert_eq!(ann.intronic_so_term(), Some("splice_donor_variant"));
    }

    #[test]
    fn test_vcf_to_genomic_hgvs_invalid_alt_index() {
        let vcf = VcfRecord::snv("chr1", 12345, 'A', 'G');
        let result = vcf_to_genomic_hgvs(&vcf, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_vcf_to_genomic_hgvs_delins() {
        let mut vcf = VcfRecord::snv("chr1", 12345, 'A', 'G');
        vcf.reference = "ATG".to_string();
        vcf.alternate = vec!["GGG".to_string()];
        let variant = vcf_to_genomic_hgvs(&vcf, 0).unwrap();
        let result_str = variant.to_string();
        assert!(result_str.contains("delins"));
    }

    #[test]
    fn test_converter_deletion() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        // Deletion at position 1050
        let vcf = VcfRecord::deletion("chr1", 1049, "ATG");
        let annotations = converter.convert(&vcf).unwrap();

        assert_eq!(annotations.len(), 1);
        let ann = &annotations[0];
        assert!(ann.is_coding);
    }

    #[test]
    fn test_converter_insertion() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        // Insertion at position 1050
        let vcf = VcfRecord::insertion("chr1", 1049, 'A', "TG");
        let annotations = converter.convert(&vcf).unwrap();

        assert_eq!(annotations.len(), 1);
        let ann = &annotations[0];
        assert!(ann.is_coding);
    }

    #[test]
    fn test_converter_noncoding_transcript() {
        // Create a non-coding transcript
        let transcript = Transcript {
            id: "NR_000001.1".to_string(),
            gene_symbol: Some("NCRNA".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 2000, 2099),
            ],
            cds_start: None, // Non-coding
            cds_end: None,
            sequence: "ATGC".repeat(50),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        let converter = VcfToHgvsConverter::new(&transcript);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let annotations = converter.convert(&vcf).unwrap();

        assert_eq!(annotations.len(), 1);
        let ann = &annotations[0];
        assert!(!ann.is_coding);
        assert!(ann.hgvs_string().contains(":n."));
    }

    #[test]
    fn test_converter_no_genomic_coords() {
        // Create transcript without genomic coordinates
        let transcript = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![Exon::new(1, 1, 100), Exon::new(2, 101, 200)],
            cds_start: Some(50),
            cds_end: Some(150),
            sequence: "ATGC".repeat(50),
            chromosome: None, // No genomic coords
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        let converter = VcfToHgvsConverter::new(&transcript);
        let vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        let annotations = converter.convert(&vcf).unwrap();

        assert_eq!(annotations.len(), 1);
        let ann = &annotations[0];
        // Should fall back to genomic notation
        assert!(ann.hgvs_string().contains(":g."));
    }

    #[test]
    fn test_converter_multi_allele() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        // Multi-allelic VCF record
        let mut vcf = VcfRecord::snv("chr1", 1050, 'A', 'G');
        vcf.alternate = vec!["G".to_string(), "T".to_string(), "C".to_string()];

        let annotations = converter.convert(&vcf).unwrap();
        assert_eq!(annotations.len(), 3);
    }

    #[test]
    fn test_determine_edit_long_deletion() {
        let transcript = create_test_transcript();
        let converter = VcfToHgvsConverter::new(&transcript);

        // Long deletion (>10bp) should not include sequence
        let reference = "AATGCATGCATGCATGC"; // 17 bp
        let alt = "A"; // 16 bp deletion
        let (edit, _start, _end) = converter.determine_edit(reference, alt).unwrap();
        assert!(matches!(edit, NaEdit::Deletion { sequence: None, .. }));
    }
}
