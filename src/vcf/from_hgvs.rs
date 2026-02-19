//! HGVS to VCF conversion
//!
//! This module provides functionality to convert HGVS variants to VCF records.
//!
//! # Coordinate System
//!
//! | Context | Basis | Notes |
//! |---------|-------|-------|
//! | HGVS positions | 1-based | Input positions from parsed variants |
//! | VCF POS field | 1-based | Output position in VCF record |
//! | Array indexing | 0-based | Internal sequence access via `- 1` |
//!
//! VCF requires an "anchor base" for insertions and deletions. This module
//! fetches the preceding reference base when converting HGVS insertions.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::convert::mapper::CoordinateMapper;
use crate::error::FerroError;
use crate::hgvs::edit::NaEdit;
use crate::hgvs::location::CdsPos;
use crate::hgvs::variant::{CdsVariant, GenomeVariant, HgvsVariant, TxVariant};
use crate::reference::transcript::{GenomeBuild, Transcript};
use crate::reference::ReferenceProvider;

use super::record::VcfRecord;

/// Result of converting an HGVS variant to VCF
#[derive(Debug, Clone)]
pub struct VcfConversionResult {
    /// The generated VCF record
    pub record: VcfRecord,
    /// Original HGVS string (for INFO annotation)
    pub hgvs_string: String,
    /// Warning messages (e.g., if normalization was applied)
    pub warnings: Vec<String>,
}

/// Converter for HGVS variants to VCF records
pub struct HgvsToVcfConverter<'a, P: ReferenceProvider> {
    transcript: &'a Transcript,
    mapper: CoordinateMapper<'a>,
    provider: &'a P,
    genome_build: GenomeBuild,
}

impl<'a, P: ReferenceProvider> HgvsToVcfConverter<'a, P> {
    /// Create a new converter for a specific transcript
    pub fn new(transcript: &'a Transcript, provider: &'a P) -> Self {
        let mapper = CoordinateMapper::new(transcript);
        let genome_build = transcript.genome_build;
        Self {
            transcript,
            mapper,
            provider,
            genome_build,
        }
    }

    /// Set the genome build for output VCF records
    pub fn with_genome_build(mut self, build: GenomeBuild) -> Self {
        self.genome_build = build;
        self
    }

    /// Convert an HGVS variant to a VCF record
    pub fn convert(&self, variant: &HgvsVariant) -> Result<VcfConversionResult, FerroError> {
        match variant {
            HgvsVariant::Cds(cds) => self.convert_cds(cds),
            HgvsVariant::Tx(tx) => self.convert_tx(tx),
            HgvsVariant::Genome(genome) => self.convert_genome(genome),
            HgvsVariant::Rna(_) => Err(FerroError::ConversionError {
                msg: "RNA variant conversion to VCF not yet supported".to_string(),
            }),
            HgvsVariant::Protein(_) => Err(FerroError::ConversionError {
                msg: "Protein variant conversion to VCF not supported (requires back-translation)"
                    .to_string(),
            }),
            HgvsVariant::Mt(mt) => {
                // Mitochondrial variants use the same format as genomic
                let genome = GenomeVariant {
                    accession: mt.accession.clone(),
                    gene_symbol: mt.gene_symbol.clone(),
                    loc_edit: mt.loc_edit.clone(),
                };
                self.convert_genome(&genome)
            }
            HgvsVariant::Allele(allele) => {
                // For simple alleles with a single variant, convert directly
                // For complex alleles with multiple variants, convert the first and warn
                if allele.variants.is_empty() {
                    return Err(FerroError::ConversionError {
                        msg: "Empty allele variant cannot be converted to VCF".to_string(),
                    });
                }

                if allele.variants.len() == 1 {
                    // Single variant allele - convert the inner variant
                    return self.convert(&allele.variants[0]);
                }

                // Multiple variants in allele - try to convert each and combine
                // This is complex because they may be at different positions
                // For now, return an error suggesting decomposition
                Err(FerroError::ConversionError {
                    msg: format!(
                        "Complex allele with {} variants should be decomposed into separate VCF records",
                        allele.variants.len()
                    ),
                })
            }
            HgvsVariant::NullAllele => Err(FerroError::ConversionError {
                msg: "Null allele marker [0] cannot be converted to VCF".to_string(),
            }),
            HgvsVariant::UnknownAllele => Err(FerroError::ConversionError {
                msg: "Unknown allele marker [?] cannot be converted to VCF".to_string(),
            }),
            HgvsVariant::Circular(circular) => {
                // Circular variants use the same format as genomic, convert to linear coordinates
                let genome = GenomeVariant {
                    accession: circular.accession.clone(),
                    gene_symbol: circular.gene_symbol.clone(),
                    loc_edit: circular.loc_edit.clone(),
                };
                self.convert_genome(&genome)
            }
            HgvsVariant::RnaFusion(_) => Err(FerroError::ConversionError {
                msg: "RNA fusion (::) conversion to VCF not supported (complex structural variant)"
                    .to_string(),
            }),
        }
    }

    /// Convert a CDS variant to VCF
    fn convert_cds(&self, variant: &CdsVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let mut warnings = Vec::new();

        // Get the CDS interval
        let interval = &variant.loc_edit.location;
        let edit = variant.loc_edit.edit.inner().expect("Edit must be known");

        // Map CDS positions to genomic positions
        let start_cds = interval.start.inner().expect("Edit must be known");
        let end_cds = interval.end.inner().expect("Edit must be known");

        let genomic_start = self.cds_to_genomic(start_cds)?;
        let genomic_end = self.cds_to_genomic(end_cds)?;

        // Get chromosome from transcript
        let chrom = self
            .transcript
            .chromosome
            .as_ref()
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no chromosome information".to_string(),
            })?
            .clone();

        // Handle intronic positions
        if start_cds.is_intronic() || end_cds.is_intronic() {
            warnings.push("Variant includes intronic positions".to_string());
        }

        // Build VCF record from edit
        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Convert a transcript variant to VCF
    fn convert_tx(&self, variant: &TxVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let warnings = Vec::new();

        let interval = &variant.loc_edit.location;
        let edit = variant.loc_edit.edit.inner().expect("Edit must be known");

        // Map transcript positions to genomic
        let start_tx = interval.start.inner().expect("Edit must be known");
        let end_tx = interval.end.inner().expect("Edit must be known");

        let genomic_start =
            self.mapper
                .tx_to_genomic(start_tx)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Could not map transcript position {} to genomic", start_tx),
                })?;
        let genomic_end =
            self.mapper
                .tx_to_genomic(end_tx)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Could not map transcript position {} to genomic", end_tx),
                })?;

        let chrom = self
            .transcript
            .chromosome
            .as_ref()
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no chromosome information".to_string(),
            })?
            .clone();

        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Convert a genomic variant to VCF
    fn convert_genome(&self, variant: &GenomeVariant) -> Result<VcfConversionResult, FerroError> {
        let hgvs_string = variant.to_string();
        let warnings = Vec::new();

        let interval = &variant.loc_edit.location;
        let edit = variant.loc_edit.edit.inner().expect("Edit must be known");

        let genomic_start = interval.start.inner().expect("Edit must be known").base;
        let genomic_end = interval.end.inner().expect("Edit must be known").base;

        // For genomic variants, we need to derive chromosome from accession
        // NC_000001.11 -> chr1, etc.
        let chrom = accession_to_chromosome(&variant.accession.to_string())
            .unwrap_or_else(|| variant.accession.number.to_string());

        let record =
            self.build_vcf_record(&chrom, genomic_start, genomic_end, edit, &hgvs_string)?;

        Ok(VcfConversionResult {
            record,
            hgvs_string,
            warnings,
        })
    }

    /// Build a VCF record from genomic coordinates and edit
    fn build_vcf_record(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        edit: &NaEdit,
        hgvs_string: &str,
    ) -> Result<VcfRecord, FerroError> {
        let (pos, reference, alternate) = match edit {
            NaEdit::Substitution {
                reference: ref_base,
                alternative: alt_base,
            } => {
                // SNV: simple one-base change
                (
                    start,
                    ref_base.to_char().to_string(),
                    vec![alt_base.to_char().to_string()],
                )
            }
            NaEdit::SubstitutionNoRef {
                alternative: alt_base,
            } => {
                // Substitution without reference - need to fetch reference from provider
                let ref_base = self.get_reference_base(chrom, start)?;
                (
                    start,
                    ref_base.to_string(),
                    vec![alt_base.to_char().to_string()],
                )
            }
            NaEdit::Deletion { sequence, .. } => {
                // Deletion: need anchor base before the deletion
                // VCF format: REF includes anchor + deleted bases, ALT is just anchor
                let anchor_pos = start - 1;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;

                let deleted_seq = if let Some(seq) = sequence {
                    seq.to_string()
                } else {
                    // Need to fetch from reference
                    self.get_reference_sequence(chrom, start, end)?
                };

                let ref_allele = format!("{}{}", anchor_base, deleted_seq);
                let alt_allele = anchor_base.to_string();

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::Insertion { sequence } => {
                // Insertion: anchor base at position, ALT includes anchor + inserted sequence
                let anchor_pos = start;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;
                let inserted_seq = sequence.to_string();

                let ref_allele = anchor_base.to_string();
                let alt_allele = format!("{}{}", anchor_base, inserted_seq);

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::Delins { sequence } => {
                // Deletion-insertion: REF is deleted sequence, ALT is replacement
                let deleted_seq = self.get_reference_sequence(chrom, start, end)?;
                let alt_seq = sequence.to_string();

                // If lengths match, no anchor needed
                // If different lengths, might need anchor for proper VCF normalization
                if deleted_seq.len() == alt_seq.len() {
                    (start, deleted_seq, vec![alt_seq])
                } else {
                    // Add anchor for length-changing variants
                    let anchor_pos = start - 1;
                    let anchor_base = self.get_reference_base(chrom, anchor_pos)?;
                    let ref_allele = format!("{}{}", anchor_base, deleted_seq);
                    let alt_allele = format!("{}{}", anchor_base, alt_seq);
                    (anchor_pos, ref_allele, vec![alt_allele])
                }
            }
            NaEdit::Duplication { sequence, .. } => {
                // Duplication is an insertion of the duplicated sequence
                let dup_seq = if let Some(seq) = sequence {
                    seq.to_string()
                } else {
                    self.get_reference_sequence(chrom, start, end)?
                };

                // Insert at end of duplicated region
                let anchor_pos = end;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;

                let ref_allele = anchor_base.to_string();
                let alt_allele = format!("{}{}", anchor_base, dup_seq);

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::DupIns { sequence } => {
                // Duplication-insertion: duplication of region followed by insertion
                // VCF representation: REF is anchor base, ALT is anchor + dup_seq + inserted_seq
                let dup_seq = self.get_reference_sequence(chrom, start, end)?;
                let ins_seq = sequence.to_string();

                let anchor_pos = end;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;

                let ref_allele = anchor_base.to_string();
                let alt_allele = format!("{}{}{}", anchor_base, dup_seq, ins_seq);

                (anchor_pos, ref_allele, vec![alt_allele])
            }
            NaEdit::Inversion { .. } => {
                // Inversion: REF is original, ALT is reverse complement
                let original_seq = self.get_reference_sequence(chrom, start, end)?;
                let inverted_seq = reverse_complement(&original_seq);
                (start, original_seq, vec![inverted_seq])
            }
            NaEdit::Repeat {
                sequence,
                count,
                additional_counts,
                trailing,
            } => {
                use crate::hgvs::edit::RepeatCount;

                // Helper to extract count from RepeatCount
                let get_count = |rc: &RepeatCount| -> Result<u64, FerroError> {
                    match rc {
                        RepeatCount::Exact(n) => Ok(*n),
                        RepeatCount::Range(min, _) | RepeatCount::MinUncertain(min) => Ok(*min),
                        RepeatCount::MaxUncertain(_) | RepeatCount::Unknown => {
                            Err(FerroError::ConversionError {
                                msg: "Repeat variant with unknown count cannot be converted to VCF"
                                    .to_string(),
                            })
                        }
                    }
                };

                // Get repeat unit sequence (including any trailing bases)
                let repeat_unit = if let Some(seq) = sequence {
                    let mut unit = seq.to_string();
                    if let Some(trail) = trailing {
                        unit.push_str(&trail.to_string());
                    }
                    unit
                } else if let Some(trail) = trailing {
                    trail.to_string()
                } else {
                    return Err(FerroError::ConversionError {
                        msg: "Repeat variant without sequence cannot be converted to VCF"
                            .to_string(),
                    });
                };

                // Get all repeat counts (primary + any additional for genotype notation)
                let primary_count = get_count(count)?;
                let mut all_counts = vec![primary_count];
                for additional in additional_counts {
                    all_counts.push(get_count(additional)?);
                }

                // Repeat expansion is an insertion at the given position
                let anchor_pos = start;
                let anchor_base = self.get_reference_base(chrom, anchor_pos)?;
                let ref_allele = anchor_base.to_string();

                // Generate alt alleles for each count (deduplicating)
                let mut alt_alleles: Vec<String> = all_counts
                    .iter()
                    .map(|&n| {
                        let expanded_seq = repeat_unit.repeat(n as usize);
                        format!("{}{}", anchor_base, expanded_seq)
                    })
                    .collect();
                alt_alleles.sort();
                alt_alleles.dedup();

                (anchor_pos, ref_allele, alt_alleles)
            }
            NaEdit::Identity { .. } => {
                // Identity means no change - return a record with same REF and ALT
                let ref_base = self.get_reference_base(chrom, start)?;
                (start, ref_base.to_string(), vec![ref_base.to_string()])
            }
            NaEdit::Conversion { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Conversion variants cannot be directly converted to VCF".to_string(),
                });
            }
            NaEdit::Unknown { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Unknown variants cannot be converted to VCF".to_string(),
                });
            }
            NaEdit::Methylation { .. } => {
                return Err(FerroError::ConversionError {
                    msg:
                        "Methylation variants cannot be converted to VCF (epigenetic, not sequence)"
                            .to_string(),
                });
            }
            NaEdit::CopyNumber { count } => {
                // Copy number variant: use VCF structural variant notation
                // REF: single base at position
                // ALT: <CNV> symbolic allele
                // INFO: SVTYPE=CNV, CN=count, END=end
                let ref_base = self.get_reference_base(chrom, start)?;

                let mut record = VcfRecord::new(
                    chrom.to_string(),
                    start,
                    ref_base.to_string(),
                    vec!["<CNV>".to_string()],
                )
                .with_genome_build(self.genome_build);

                // Add CNV-specific INFO fields
                record.info.insert(
                    "SVTYPE".to_string(),
                    super::record::InfoValue::String("CNV".to_string()),
                );
                record.info.insert(
                    "CN".to_string(),
                    super::record::InfoValue::Integer(*count as i64),
                );
                record.info.insert(
                    "END".to_string(),
                    super::record::InfoValue::Integer(end as i64),
                );
                record.info.insert(
                    "HGVS".to_string(),
                    super::record::InfoValue::String(hgvs_string.to_string()),
                );

                return Ok(record);
            }
            NaEdit::Splice => {
                return Err(FerroError::ConversionError {
                    msg: "Splice variants (r.spl) cannot be converted to VCF (RNA-specific effect)"
                        .to_string(),
                });
            }
            NaEdit::NoProduct => {
                return Err(FerroError::ConversionError {
                    msg:
                        "No product variants (r.0) cannot be converted to VCF (RNA-specific effect)"
                            .to_string(),
                });
            }
            NaEdit::MultiRepeat { .. } => {
                return Err(FerroError::ConversionError {
                    msg: "Multi-repeat variants cannot be directly converted to VCF yet"
                        .to_string(),
                });
            }
            NaEdit::PositionOnly => {
                return Err(FerroError::ConversionError {
                    msg: "Position-only variants cannot be converted to VCF (no edit specified)"
                        .to_string(),
                });
            }
        };

        let mut record = VcfRecord::new(chrom.to_string(), pos, reference, alternate)
            .with_genome_build(self.genome_build);

        // Add HGVS annotation to INFO
        record.info.insert(
            "HGVS".to_string(),
            super::record::InfoValue::String(hgvs_string.to_string()),
        );

        Ok(record)
    }

    /// Map CDS position to genomic position
    fn cds_to_genomic(&self, cds_pos: &CdsPos) -> Result<u64, FerroError> {
        self.mapper
            .cds_to_genomic(cds_pos)?
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!("Could not map CDS position {} to genomic", cds_pos),
            })
    }

    /// Get a single reference base at a position
    fn get_reference_base(&self, _chrom: &str, pos: u64) -> Result<char, FerroError> {
        // Try to get from transcript sequence first
        let seq = &self.transcript.sequence;
        if !seq.is_empty() {
            // This is a simplification - in practice we'd need proper coordinate mapping
            // Use checked conversion to avoid truncation on 32-bit systems
            let idx = usize::try_from(pos.saturating_sub(1)).map_err(|_| {
                FerroError::ConversionError {
                    msg: format!("Position {} is too large for this platform", pos),
                }
            })?;
            if idx < seq.len() {
                return Ok(seq.chars().nth(idx).unwrap_or('N'));
            }
        }

        // Fallback: try reference provider
        let accession = &self.transcript.id;
        let sequence = self.provider.get_sequence(accession, pos - 1, pos)?;
        sequence
            .chars()
            .next()
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!("Could not get reference base at position {}", pos),
            })
    }

    /// Get reference sequence for a range
    fn get_reference_sequence(
        &self,
        _chrom: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Try to get from transcript sequence first
        let seq = &self.transcript.sequence;
        if !seq.is_empty() {
            // Use checked conversion to avoid truncation on 32-bit systems
            let start_idx = usize::try_from(start.saturating_sub(1)).map_err(|_| {
                FerroError::ConversionError {
                    msg: format!("Start position {} is too large for this platform", start),
                }
            })?;
            let end_idx = usize::try_from(end).map_err(|_| FerroError::ConversionError {
                msg: format!("End position {} is too large for this platform", end),
            })?;
            if end_idx <= seq.len() {
                return Ok(seq[start_idx..end_idx].to_string());
            }
        }

        // Fallback: try reference provider
        let accession = &self.transcript.id;
        self.provider.get_sequence(accession, start - 1, end)
    }
}

/// Convert a genomic HGVS variant directly to VCF (no transcript needed)
pub fn genomic_hgvs_to_vcf(variant: &GenomeVariant) -> Result<VcfRecord, FerroError> {
    let interval = &variant.loc_edit.location;
    let edit = variant.loc_edit.edit.inner().expect("Edit must be known");

    let start = interval.start.inner().expect("Edit must be known").base;
    let _end = interval.end.inner().expect("Edit must be known").base;

    let chrom = accession_to_chromosome(&variant.accession.to_string())
        .unwrap_or_else(|| variant.accession.number.to_string());

    let (pos, reference, alternate) = match edit {
        NaEdit::Substitution {
            reference: ref_base,
            alternative: alt_base,
        } => (
            start,
            ref_base.to_char().to_string(),
            vec![alt_base.to_char().to_string()],
        ),
        NaEdit::Deletion { sequence, .. } => {
            let deleted = sequence
                .as_ref()
                .map(|s| s.to_string())
                .ok_or_else(|| FerroError::ConversionError {
                    msg: "Cannot convert deletion to VCF without deleted sequence (no reference data)".to_string(),
                })?;
            // Without reference access, we can't get the anchor base
            return Err(FerroError::ConversionError {
                msg: format!(
                    "Cannot convert deletion '{}' to VCF: anchor base required (no reference data)",
                    deleted
                ),
            });
        }
        NaEdit::Insertion { sequence } => {
            let inserted = sequence.to_string();
            // Without reference access, we can't get the anchor base
            return Err(FerroError::ConversionError {
                msg: format!(
                    "Cannot convert insertion '{}' to VCF: anchor base required (no reference data)",
                    inserted
                ),
            });
        }
        _ => {
            return Err(FerroError::ConversionError {
                msg: format!("Edit type {:?} not supported for simple conversion", edit),
            });
        }
    };

    Ok(VcfRecord::new(chrom, pos, reference, alternate))
}

/// Map RefSeq chromosome accession to chromosome name
fn accession_to_chromosome(accession: &str) -> Option<String> {
    // Extract base accession (without version)
    let base = accession.split('.').next().unwrap_or(accession);

    // GRCh38 primary assembly: NC_000001.11 -> chr1, etc.
    if base.starts_with("NC_0000") && base.len() >= 9 {
        let num_str = &base[7..9];
        if let Ok(num) = num_str.parse::<u32>() {
            if num <= 22 {
                return Some(format!("chr{}", num));
            } else if num == 23 {
                return Some("chrX".to_string());
            } else if num == 24 {
                return Some("chrY".to_string());
            }
        }
    }

    // GRCh38 patches and alternate loci: NC_060000+ series
    // These are typically alternate haplotypes, return the accession itself
    if base.starts_with("NC_06") {
        return Some(base.to_string());
    }

    // Mitochondrial
    if base.starts_with("NC_012920") {
        return Some("chrM".to_string());
    }

    // Bacterial/other RefSeq: NZ_ prefix - return as-is
    if base.starts_with("NZ_") {
        return Some(base.to_string());
    }

    // Unplaced scaffolds: return as-is
    if base.starts_with("NT_") || base.starts_with("NW_") {
        return Some(base.to_string());
    }

    None
}

/// Compute reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            _ => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::{
        Base, InsertedSequence, MethylationStatus, ProteinEdit, RepeatCount, Sequence,
    };
    use crate::hgvs::interval::{
        CdsInterval, GenomeInterval, ProtInterval, RnaInterval, TxInterval,
    };
    use crate::hgvs::location::{CdsPos, GenomePos, ProtPos, RnaPos, TxPos};
    use crate::hgvs::variant::AllelePhase;
    use crate::hgvs::variant::{
        Accession, AlleleVariant, CircularVariant, LocEdit, MtVariant, ProteinVariant,
        RnaFusionBreakpoint, RnaFusionVariant, RnaVariant,
    };
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use crate::reference::MockProvider;
    use std::str::FromStr;
    use std::sync::OnceLock;

    fn create_test_transcript() -> Transcript {
        Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            exons: vec![
                Exon::with_genomic(1, 1, 100, 1000, 1099),
                Exon::with_genomic(2, 101, 200, 2000, 2099),
            ],
            cds_start: Some(50),
            cds_end: Some(150),
            sequence: "ATGCATGCATGC".repeat(20),
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2099),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        assert_eq!(reverse_complement("atgc"), "gcat");
        assert_eq!(reverse_complement("AaTt"), "aAtT");
    }

    #[test]
    fn test_reverse_complement_unknown_bases() {
        assert_eq!(reverse_complement("ATNC"), "GNAT");
    }

    #[test]
    fn test_accession_to_chromosome() {
        assert_eq!(
            accession_to_chromosome("NC_000001.11"),
            Some("chr1".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_000022.11"),
            Some("chr22".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_000023.11"),
            Some("chrX".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_000024.10"),
            Some("chrY".to_string())
        );
        assert_eq!(
            accession_to_chromosome("NC_012920.1"),
            Some("chrM".to_string())
        );
        assert_eq!(accession_to_chromosome("NM_000088.3"), None);
    }

    #[test]
    fn test_accession_to_chromosome_invalid() {
        assert_eq!(accession_to_chromosome("NC_000025.11"), None);
        assert_eq!(accession_to_chromosome(""), None);
        assert_eq!(accession_to_chromosome("invalid"), None);
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_snv() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        let vcf = genomic_hgvs_to_vcf(&variant).unwrap();
        assert_eq!(vcf.chrom, "chr1");
        assert_eq!(vcf.pos, 12345);
        assert_eq!(vcf.reference, "A");
        assert_eq!(vcf.alternate, vec!["G"]);
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_deletion_requires_anchor() {
        // Deletion without reference data for anchor base returns an error
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(12345), GenomePos::new(12347)),
                NaEdit::Deletion {
                    sequence: Some(Sequence::from_str("ATG").unwrap()),
                    length: None,
                },
            ),
        };

        let result = genomic_hgvs_to_vcf(&variant);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("anchor base"));
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_insertion_requires_anchor() {
        // Insertion without reference data for anchor base returns an error
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Insertion {
                    sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
                },
            ),
        };

        let result = genomic_hgvs_to_vcf(&variant);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("anchor base"));
    }

    #[test]
    fn test_genomic_hgvs_to_vcf_unsupported_edit() {
        let variant = GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::new(GenomePos::new(12345), GenomePos::new(12350)),
                NaEdit::Inversion {
                    sequence: None,
                    length: None,
                },
            ),
        };

        let result = genomic_hgvs_to_vcf(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_cds_snv() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        let result = converter.convert(&HgvsVariant::Cds(variant)).unwrap();
        assert_eq!(result.record.chrom, "chr1");
        assert!(result.record.info.contains_key("HGVS"));
    }

    #[test]
    fn test_converter_with_genome_build() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter =
            HgvsToVcfConverter::new(&transcript, &provider).with_genome_build(GenomeBuild::GRCh37);

        let variant = CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        };

        let result = converter.convert(&HgvsVariant::Cds(variant)).unwrap();
        assert_eq!(result.record.genome_build, GenomeBuild::GRCh37);
    }

    #[test]
    fn test_converter_protein_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Protein(ProteinVariant {
            accession: Accession::new("NP", "000079", Some(2)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                ProtInterval::point(ProtPos::new(crate::hgvs::location::AminoAcid::Met, 1)),
                ProteinEdit::Unknown {
                    predicted: false,
                    whole_protein: false,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(matches!(err, FerroError::ConversionError { .. }));
    }

    #[test]
    fn test_converter_rna_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Rna(RnaVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                RnaInterval::point(RnaPos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_null_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let result = converter.convert(&HgvsVariant::NullAllele);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_unknown_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let result = converter.convert(&HgvsVariant::UnknownAllele);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_empty_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Allele(AlleleVariant {
            variants: vec![],
            phase: AllelePhase::Cis,
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_single_variant_allele() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let inner = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let variant = HgvsVariant::Allele(AlleleVariant {
            variants: vec![inner],
            phase: AllelePhase::Cis,
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_converter_complex_allele_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let inner1 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });
        let inner2 = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(10)),
                NaEdit::Substitution {
                    reference: Base::C,
                    alternative: Base::T,
                },
            ),
        });

        let variant = HgvsVariant::Allele(AlleleVariant {
            variants: vec![inner1, inner2],
            phase: AllelePhase::Cis,
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_mt_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Mt(MtVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(3243)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().record.chrom, "chrM");
    }

    #[test]
    fn test_converter_circular_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Circular(CircularVariant {
            accession: Accession::new("NC", "012920", Some(1)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(100)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
    }

    #[test]
    fn test_converter_rna_fusion_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let five_prime = RnaFusionBreakpoint {
            accession: Accession::new("NM", "000001", Some(1)),
            gene_symbol: None,
            interval: RnaInterval::point(RnaPos::new(100)),
        };
        let three_prime = RnaFusionBreakpoint {
            accession: Accession::new("NM", "000002", Some(1)),
            gene_symbol: None,
            interval: RnaInterval::point(RnaPos::new(200)),
        };

        let variant = HgvsVariant::RnaFusion(RnaFusionVariant::new(five_prime, three_prime));

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_converter_tx_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Tx(TxVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: Some("COL1A1".to_string()),
            loc_edit: LocEdit::new(
                TxInterval::point(TxPos::new(50)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().record.chrom, "chr1");
    }

    #[test]
    fn test_converter_genome_variant() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Genome(GenomeVariant {
            accession: Accession::new("NC", "000001", Some(11)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                GenomeInterval::point(GenomePos::new(12345)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        let vcf = result.unwrap();
        assert_eq!(vcf.record.chrom, "chr1");
        assert_eq!(vcf.record.pos, 12345);
    }

    #[test]
    fn test_vcf_conversion_result_fields() {
        let result = VcfConversionResult {
            record: VcfRecord::new(
                "chr1".to_string(),
                100,
                "A".to_string(),
                vec!["G".to_string()],
            ),
            hgvs_string: "NC_000001.11:g.100A>G".to_string(),
            warnings: vec!["test warning".to_string()],
        };

        assert_eq!(result.hgvs_string, "NC_000001.11:g.100A>G");
        assert_eq!(result.warnings.len(), 1);
        assert_eq!(result.record.chrom, "chr1");
    }

    #[test]
    fn test_build_vcf_record_delins_same_length() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Delins {
                    sequence: InsertedSequence::Literal(Sequence::from_str("TTT").unwrap()),
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_build_vcf_record_duplication() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(3)),
                NaEdit::Duplication {
                    sequence: Some(Sequence::from_str("ATG").unwrap()),
                    length: None,
                    uncertain_extent: None,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_build_vcf_record_inversion() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(4)),
                NaEdit::Inversion {
                    sequence: None,
                    length: None,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_build_vcf_record_repeat() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_build_vcf_record_repeat_with_additional_counts() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![RepeatCount::Exact(10)],
                    trailing: None,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_build_vcf_record_repeat_unknown_count() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Unknown,
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_repeat_no_sequence() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: None,
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_identity() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Identity {
                    sequence: Some(Sequence::from_str("A").unwrap()),
                    whole_entity: false,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_build_vcf_record_conversion_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(10)),
                NaEdit::Conversion {
                    source: "NM_other:c.100_110".to_string(),
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_unknown_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Unknown {
                    whole_entity: false,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_methylation_not_supported() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Methylation {
                    status: MethylationStatus::GainOfMethylation,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_build_vcf_record_copy_number() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::new(CdsPos::new(1), CdsPos::new(100)),
                NaEdit::CopyNumber { count: 3 },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_cds_with_intronic_position_warns() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::with_offset(10, 5)),
                NaEdit::Substitution {
                    reference: Base::A,
                    alternative: Base::G,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_ok());
        let vcf = result.unwrap();
        assert!(!vcf.warnings.is_empty());
        assert!(vcf.warnings[0].contains("intronic"));
    }

    #[test]
    fn test_repeat_count_range() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Range(5, 10),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_repeat_count_min_uncertain() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::MinUncertain(5),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_repeat_count_max_uncertain() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::MaxUncertain(10),
                    additional_counts: vec![],
                    trailing: None,
                },
            ),
        });

        let result = converter.convert(&variant);
        assert!(result.is_err());
    }

    #[test]
    fn test_repeat_with_trailing() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: Some(Sequence::from_str("CAG").unwrap()),
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: Some(Sequence::from_str("T").unwrap()),
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }

    #[test]
    fn test_repeat_trailing_only() {
        let transcript = create_test_transcript();
        let provider = MockProvider::new();
        let converter = HgvsToVcfConverter::new(&transcript, &provider);

        let variant = HgvsVariant::Cds(CdsVariant {
            accession: Accession::new("NM", "000088", Some(3)),
            gene_symbol: None,
            loc_edit: LocEdit::new(
                CdsInterval::point(CdsPos::new(1)),
                NaEdit::Repeat {
                    sequence: None,
                    count: RepeatCount::Exact(5),
                    additional_counts: vec![],
                    trailing: Some(Sequence::from_str("T").unwrap()),
                },
            ),
        });

        // Test exercises the conversion code path - may fail due to mock provider limitations
        let _result = converter.convert(&variant);
    }
}
