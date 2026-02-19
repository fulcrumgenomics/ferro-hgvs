//! Coordinate mapper for converting between coordinate systems
//!
//! # Coordinate Systems
//!
//! This module handles conversions between multiple coordinate systems:
//!
//! | System | Basis | Notes |
//! |--------|-------|-------|
//! | Genomic | 0-based | Half-open intervals for exons |
//! | Transcript (tx) | 1-based | `TxPos.base` is 1-based |
//! | CDS (c.) | 1-based | `CdsPos.base` is 1-based, negative for 5'UTR |
//! | Protein (p.) | 1-based | `ProtPos.number` is 1-based |
//!
//! ## Key conversions:
//! - CDS → Tx: `cds_to_tx()` - accounts for CDS start/end and UTR regions
//! - Tx → CDS: `tx_to_cds()` - inverse of above
//! - Genomic → Tx: `genomic_to_tx()` - uses exon boundaries
//! - Tx → Genomic: `tx_to_genomic()` - accounts for strand
//! - CDS → Protein: `cds_to_protein()` - divides by 3 for codon position
//!
//! ## Intronic positions:
//! CDS positions can have intronic offsets (e.g., c.100+5, c.200-10).
//! Use `*_with_intron()` methods to handle these cases.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::hgvs::location::{AminoAcid, CdsPos, ProtPos, TxPos};
use crate::reference::transcript::{IntronPosition, Strand, Transcript};

/// Maps coordinates between different systems for a transcript
pub struct CoordinateMapper<'a> {
    transcript: &'a Transcript,
}

impl<'a> CoordinateMapper<'a> {
    /// Create a new mapper for a transcript
    pub fn new(transcript: &'a Transcript) -> Self {
        Self { transcript }
    }

    /// Convert CDS position to transcript position
    ///
    /// CDS positions are relative to the start of the coding sequence,
    /// while transcript positions are absolute on the transcript.
    ///
    /// This method is exon-aware: it properly maps CDS positions through
    /// exon boundaries, accounting for gaps in cdot's coordinate system.
    pub fn cds_to_tx(&self, pos: &CdsPos) -> Result<TxPos, FerroError> {
        let cds_start = self
            .transcript
            .cds_start
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS".to_string(),
            })?;

        let cds_end = self
            .transcript
            .cds_end
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS end".to_string(),
            })?;

        let tx_base = if pos.utr3 {
            // 3' UTR: *N = CDS_end + N
            // For UTR, use exon-aware mapping from CDS end
            self.cds_to_tx_exon_aware(cds_end as i64, pos.base)?
        } else if pos.base < 1 {
            // 5' UTR: -N = CDS_start - N (where N is positive in HGVS but stored as negative)
            // e.g., c.-3 means 3 bases before CDS start, so offset is -3
            self.cds_to_tx_exon_aware(cds_start as i64, pos.base)?
        } else {
            // Normal CDS position - use exon-aware mapping
            // c.1 = cds_start (offset 0), c.2 = cds_start + 1 (offset 1), etc.
            self.cds_to_tx_exon_aware(cds_start as i64, pos.base - 1)?
        };

        Ok(TxPos {
            base: tx_base,
            offset: pos.offset,
            downstream: false,
        })
    }

    /// Convert CDS position to tx position using exon-aware mapping
    ///
    /// This walks through exons, counting bases in the CDS to properly
    /// handle gaps in cdot's coordinate system.
    ///
    /// # Arguments
    /// * `start_tx` - Starting transcript position (e.g., cds_start)
    /// * `offset` - Number of bases to offset from start (can be negative for 5' UTR)
    fn cds_to_tx_exon_aware(&self, start_tx: i64, offset: i64) -> Result<i64, FerroError> {
        // If no exon data or offset is 0, use simple calculation
        if self.transcript.exons.is_empty() || offset == 0 {
            return Ok(start_tx + offset);
        }

        // Sort exons by transcript position
        let mut sorted_exons: Vec<_> = self.transcript.exons.iter().collect();
        sorted_exons.sort_by_key(|e| e.start);

        // Check if exons are contiguous (no gaps in tx coordinates)
        // Most real transcripts from cdot are contiguous
        let has_gaps = sorted_exons.windows(2).any(|w| w[0].end + 1 != w[1].start);
        if !has_gaps {
            // No gaps - use simple calculation
            return Ok(start_tx + offset);
        }

        // Exons have gaps - use exon-aware mapping
        if offset > 0 {
            // Moving forward through exons
            let mut remaining = offset;
            let mut current_tx = start_tx as u64;

            // Find which exon contains start_tx
            let mut found_start = false;
            for exon in &sorted_exons {
                if current_tx >= exon.start && current_tx <= exon.end {
                    // Start is in this exon
                    found_start = true;
                    let bases_in_exon = exon.end - current_tx;
                    if remaining <= bases_in_exon as i64 {
                        // Target is in this exon
                        return Ok((current_tx + remaining as u64) as i64);
                    }
                    remaining -= bases_in_exon as i64;
                    // Move to next exon's start (skip gap)
                    current_tx = exon.end + 1;
                } else if found_start && current_tx < exon.start {
                    // We're in a gap, jump to this exon's start
                    let bases_in_exon = exon.end - exon.start + 1;
                    if remaining <= bases_in_exon as i64 {
                        return Ok((exon.start + remaining as u64 - 1) as i64);
                    }
                    remaining -= bases_in_exon as i64;
                    current_tx = exon.end + 1;
                }
            }
            // If we get here, position is beyond all exons
            Ok(current_tx as i64 + remaining - 1)
        } else {
            // Moving backward through exons (for 5' UTR)
            let mut remaining = -offset;
            let mut current_tx = start_tx as u64;

            // Find which exon contains start_tx and go backward
            let mut found_start = false;
            for exon in sorted_exons.iter().rev() {
                if current_tx >= exon.start && current_tx <= exon.end {
                    found_start = true;
                    let bases_before = current_tx - exon.start;
                    if remaining <= bases_before as i64 {
                        return Ok((current_tx - remaining as u64) as i64);
                    }
                    remaining -= bases_before as i64;
                    current_tx = exon.start - 1;
                } else if found_start && current_tx > exon.end {
                    // We're in a gap, jump to this exon's end
                    let bases_in_exon = exon.end - exon.start + 1;
                    if remaining <= bases_in_exon as i64 {
                        return Ok((exon.end - remaining as u64 + 1) as i64);
                    }
                    remaining -= bases_in_exon as i64;
                    current_tx = exon.start - 1;
                }
            }
            Ok(current_tx as i64 - remaining + 1)
        }
    }

    /// Convert transcript position to CDS position
    ///
    /// This method is exon-aware: it properly maps tx positions to CDS positions
    /// accounting for gaps in cdot's coordinate system.
    pub fn tx_to_cds(&self, pos: &TxPos) -> Result<CdsPos, FerroError> {
        let cds_start = self
            .transcript
            .cds_start
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS".to_string(),
            })? as i64;
        let cds_end = self
            .transcript
            .cds_end
            .ok_or_else(|| FerroError::ConversionError {
                msg: "Transcript has no CDS end".to_string(),
            })? as i64;

        let base = pos.base;

        if base < cds_start {
            // 5' UTR - count bases backward from CDS start through exons
            let cds_offset = self.tx_to_cds_exon_aware(base, cds_start)?;
            Ok(CdsPos {
                base: cds_offset, // Will be negative for 5' UTR
                offset: pos.offset,
                utr3: false,
            })
        } else if base > cds_end {
            // 3' UTR - count bases forward from CDS end through exons
            let cds_offset = self.tx_to_cds_exon_aware(base, cds_end)?;
            Ok(CdsPos {
                base: cds_offset - 1, // Adjust for 3' UTR notation (*1, *2, etc.)
                offset: pos.offset,
                utr3: true,
            })
        } else {
            // Within CDS - use exon-aware mapping
            let cds_pos = self.tx_to_cds_exon_aware(base, cds_start)?;
            Ok(CdsPos {
                base: cds_pos,
                offset: pos.offset,
                utr3: false,
            })
        }
    }

    /// Convert tx position to CDS position using exon-aware mapping
    ///
    /// This counts actual exonic bases between start_tx and target_tx,
    /// skipping any gaps in cdot's coordinate system.
    ///
    /// Returns the number of CDS bases from start to target (can be negative).
    fn tx_to_cds_exon_aware(&self, target_tx: i64, start_tx: i64) -> Result<i64, FerroError> {
        // If no exon data, use simple calculation
        if self.transcript.exons.is_empty() {
            // For 5' UTR (target < start), formula is target - start
            // For CDS (target >= start), formula is target - start + 1
            if target_tx >= start_tx {
                return Ok(target_tx - start_tx + 1);
            } else {
                return Ok(target_tx - start_tx);
            }
        }

        // Sort exons by transcript position
        let mut sorted_exons: Vec<_> = self.transcript.exons.iter().collect();
        sorted_exons.sort_by_key(|e| e.start);

        // Check if exons are contiguous (no gaps in tx coordinates)
        // Most real transcripts from cdot are contiguous
        let has_gaps = sorted_exons.windows(2).any(|w| w[0].end + 1 != w[1].start);
        if !has_gaps {
            // No gaps - use simple calculation
            // For 5' UTR (target < start), formula is target - start
            // For CDS (target >= start), formula is target - start + 1
            if target_tx >= start_tx {
                return Ok(target_tx - start_tx + 1);
            } else {
                return Ok(target_tx - start_tx);
            }
        }

        // Exons have gaps - use exon-aware counting
        if target_tx >= start_tx {
            // Counting forward
            let mut cds_count: i64 = 0;
            let mut counting = false;

            for exon in &sorted_exons {
                let exon_start = exon.start as i64;
                let exon_end = exon.end as i64;

                if !counting && exon_end >= start_tx && exon_start <= start_tx {
                    // Start counting from here
                    counting = true;
                    let count_start = start_tx.max(exon_start);
                    let count_end = target_tx.min(exon_end);
                    if count_end >= count_start {
                        cds_count += count_end - count_start + 1;
                    }
                    if target_tx <= exon_end {
                        return Ok(cds_count);
                    }
                } else if counting {
                    // Continue counting in subsequent exons
                    if target_tx < exon_start {
                        // Target is in a gap before this exon
                        return Ok(cds_count);
                    }
                    let count_end = target_tx.min(exon_end);
                    cds_count += count_end - exon_start + 1;
                    if target_tx <= exon_end {
                        return Ok(cds_count);
                    }
                }
            }
            Ok(cds_count)
        } else {
            // Counting backward (for 5' UTR)
            let mut cds_count: i64 = 0;
            let mut counting = false;

            for exon in sorted_exons.iter().rev() {
                let exon_start = exon.start as i64;
                let exon_end = exon.end as i64;

                if !counting && exon_start <= start_tx && exon_end >= start_tx {
                    // Start counting backward from here
                    counting = true;
                    let count_end = start_tx.min(exon_end);
                    let count_start = target_tx.max(exon_start);
                    if count_end >= count_start {
                        cds_count -= count_end - count_start + 1;
                    }
                    if target_tx >= exon_start {
                        return Ok(cds_count + 1); // Adjust for 1-based
                    }
                } else if counting {
                    if target_tx > exon_end {
                        return Ok(cds_count);
                    }
                    let count_start = target_tx.max(exon_start);
                    cds_count -= exon_end - count_start + 1;
                    if target_tx >= exon_start {
                        return Ok(cds_count + 1);
                    }
                }
            }
            Ok(cds_count)
        }
    }

    /// Convert CDS position to protein position
    ///
    /// Protein positions are (CDS position + 2) / 3 (rounded up)
    pub fn cds_to_protein(&self, pos: &CdsPos) -> Result<ProtPos, FerroError> {
        if pos.base < 1 || pos.utr3 {
            return Err(FerroError::ConversionError {
                msg: "Cannot convert UTR position to protein".to_string(),
            });
        }

        if pos.is_intronic() {
            return Err(FerroError::ConversionError {
                msg: "Cannot convert intronic position to protein".to_string(),
            });
        }

        // Protein position is 1-indexed, and each codon is 3 bases
        let aa_number = ((pos.base - 1) / 3 + 1) as u64;

        // Get the amino acid at this position if we have the sequence
        // For now, use Xaa (unknown)
        let aa = self.get_amino_acid_at(aa_number).unwrap_or(AminoAcid::Xaa);

        Ok(ProtPos::new(aa, aa_number))
    }

    /// Get the amino acid at a protein position (if sequence is available)
    fn get_amino_acid_at(&self, position: u64) -> Option<AminoAcid> {
        let cds_start = self.transcript.cds_start?;
        let cds_end = self.transcript.cds_end?;

        // Get codon start position (0-based in sequence)
        let codon_start = cds_start as usize - 1 + (position as usize - 1) * 3;
        let codon_end = codon_start + 3;

        if codon_end > cds_end as usize || codon_end > self.transcript.sequence.len() {
            return None;
        }

        let codon = &self.transcript.sequence[codon_start..codon_end];
        Self::translate_codon(codon)
    }

    /// Translate a codon to an amino acid
    fn translate_codon(codon: &str) -> Option<AminoAcid> {
        match codon.to_uppercase().as_str() {
            "TTT" | "TTC" => Some(AminoAcid::Phe),
            "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => Some(AminoAcid::Leu),
            "ATT" | "ATC" | "ATA" => Some(AminoAcid::Ile),
            "ATG" => Some(AminoAcid::Met),
            "GTT" | "GTC" | "GTA" | "GTG" => Some(AminoAcid::Val),
            "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => Some(AminoAcid::Ser),
            "CCT" | "CCC" | "CCA" | "CCG" => Some(AminoAcid::Pro),
            "ACT" | "ACC" | "ACA" | "ACG" => Some(AminoAcid::Thr),
            "GCT" | "GCC" | "GCA" | "GCG" => Some(AminoAcid::Ala),
            "TAT" | "TAC" => Some(AminoAcid::Tyr),
            "TAA" | "TAG" | "TGA" => Some(AminoAcid::Ter),
            "CAT" | "CAC" => Some(AminoAcid::His),
            "CAA" | "CAG" => Some(AminoAcid::Gln),
            "AAT" | "AAC" => Some(AminoAcid::Asn),
            "AAA" | "AAG" => Some(AminoAcid::Lys),
            "GAT" | "GAC" => Some(AminoAcid::Asp),
            "GAA" | "GAG" => Some(AminoAcid::Glu),
            "TGT" | "TGC" => Some(AminoAcid::Cys),
            "TGG" => Some(AminoAcid::Trp),
            "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => Some(AminoAcid::Arg),
            "GGT" | "GGC" | "GGA" | "GGG" => Some(AminoAcid::Gly),
            _ => Some(AminoAcid::Xaa), // Unknown
        }
    }

    /// Convert genomic position to transcript position
    ///
    /// Returns None if the genomic position is in an intron.
    /// The transcript must have genomic coordinates set on its exons.
    pub fn genomic_to_tx(&self, genomic_pos: u64) -> Result<Option<TxPos>, FerroError> {
        // Check if transcript has genomic coordinates
        if !self.transcript.has_genomic_coords() {
            return Err(FerroError::ConversionError {
                msg: "Transcript does not have genomic coordinates".to_string(),
            });
        }

        // Find which exon contains this genomic position
        for exon in &self.transcript.exons {
            let (g_start, g_end) = match (exon.genomic_start, exon.genomic_end) {
                (Some(s), Some(e)) => (s, e),
                _ => continue,
            };

            // Check if position is within this exon's genomic coordinates
            if genomic_pos >= g_start && genomic_pos <= g_end {
                // Calculate offset within exon based on strand
                let tx_pos = match self.transcript.strand {
                    Strand::Plus => {
                        // Plus strand: genomic position increases with transcript position
                        let offset_in_exon = genomic_pos - g_start;
                        exon.start + offset_in_exon
                    }
                    Strand::Minus => {
                        // Minus strand: genomic position decreases with transcript position
                        let offset_in_exon = g_end - genomic_pos;
                        exon.start + offset_in_exon
                    }
                };

                return Ok(Some(TxPos::new(tx_pos as i64)));
            }
        }

        // Position is in an intron (between exons)
        Ok(None)
    }

    /// Convert genomic position to transcript position with intronic offset support
    ///
    /// Unlike `genomic_to_tx`, this method returns a TxPos with an offset for
    /// intronic positions instead of returning None.
    pub fn genomic_to_tx_with_intron(&self, genomic_pos: u64) -> Result<TxPos, FerroError> {
        // First check if it's exonic
        if let Some(tx_pos) = self.genomic_to_tx(genomic_pos)? {
            return Ok(tx_pos);
        }

        // It's intronic - find the intron and calculate offset
        if let Some((_intron, intron_pos)) = self.transcript.find_intron_at_genomic(genomic_pos) {
            return Ok(TxPos::with_offset(
                intron_pos.tx_boundary_pos as i64,
                intron_pos.offset,
            ));
        }

        Err(FerroError::ConversionError {
            msg: format!(
                "Genomic position {} is outside transcript bounds",
                genomic_pos
            ),
        })
    }

    /// Get intronic position information for a genomic position
    ///
    /// Returns None if the position is exonic or outside the transcript.
    pub fn get_intron_position(&self, genomic_pos: u64) -> Option<IntronPosition> {
        self.transcript
            .find_intron_at_genomic(genomic_pos)
            .map(|(_, pos)| pos)
    }

    /// Check if a genomic position is intronic
    pub fn is_intronic_at_genomic(&self, genomic_pos: u64) -> bool {
        self.transcript
            .find_intron_at_genomic(genomic_pos)
            .is_some()
    }

    /// Convert transcript position to genomic position
    ///
    /// Returns None if the transcript position is not covered by exons with genomic coords.
    pub fn tx_to_genomic(&self, tx_pos: &TxPos) -> Result<Option<u64>, FerroError> {
        // Check if transcript has genomic coordinates
        if !self.transcript.has_genomic_coords() {
            return Err(FerroError::ConversionError {
                msg: "Transcript does not have genomic coordinates".to_string(),
            });
        }

        // Find which exon contains this transcript position
        // Only positive transcript positions can be in exons
        if tx_pos.base < 1 {
            return Ok(None);
        }
        let tx_base = tx_pos.base as u64;

        for exon in &self.transcript.exons {
            if tx_base >= exon.start && tx_base <= exon.end {
                let (g_start, g_end) = match (exon.genomic_start, exon.genomic_end) {
                    (Some(s), Some(e)) => (s, e),
                    _ => continue,
                };

                // Calculate offset within exon
                let offset_in_exon = tx_base - exon.start;

                let genomic_pos = match self.transcript.strand {
                    Strand::Plus => {
                        // Plus strand: transcript position increases with genomic position
                        g_start + offset_in_exon
                    }
                    Strand::Minus => {
                        // Minus strand: transcript position increases as genomic position decreases
                        g_end - offset_in_exon
                    }
                };

                return Ok(Some(genomic_pos));
            }
        }

        // Position not found in exons
        Ok(None)
    }

    /// Convert genomic position to CDS position
    ///
    /// Returns None if the position is intronic.
    pub fn genomic_to_cds(&self, genomic_pos: u64) -> Result<Option<CdsPos>, FerroError> {
        // First convert to transcript position
        let tx_pos = self.genomic_to_tx(genomic_pos)?;

        match tx_pos {
            Some(pos) => Ok(Some(self.tx_to_cds(&pos)?)),
            None => Ok(None),
        }
    }

    /// Convert genomic position to CDS position with intronic offset support
    ///
    /// Unlike `genomic_to_cds`, this method returns a CdsPos with an offset for
    /// intronic positions instead of returning None.
    pub fn genomic_to_cds_with_intron(&self, genomic_pos: u64) -> Result<CdsPos, FerroError> {
        // Get tx position with intron support
        let tx_pos = self.genomic_to_tx_with_intron(genomic_pos)?;

        // Convert to CDS, preserving the offset
        self.tx_to_cds(&tx_pos)
    }

    /// Convert CDS position to genomic position
    pub fn cds_to_genomic(&self, cds_pos: &CdsPos) -> Result<Option<u64>, FerroError> {
        // First convert to transcript position
        let tx_pos = self.cds_to_tx(cds_pos)?;
        self.tx_to_genomic(&tx_pos)
    }

    /// Get the chromosome name for this transcript
    pub fn chromosome(&self) -> Option<&str> {
        self.transcript.chromosome.as_deref()
    }

    /// Get the strand for this transcript
    pub fn strand(&self) -> Strand {
        self.transcript.strand
    }

    /// Convert CDS position with intronic offset to genomic position
    ///
    /// For intronic variants like c.100+5 or c.200-10, this calculates the
    /// genomic position by:
    /// 1. Converting the CDS base position to transcript position
    /// 2. Using the intron mapping to find the genomic position
    ///
    /// # Returns
    /// The genomic position, or an error if the transcript lacks genomic coordinates
    pub fn cds_to_genomic_with_intron(&self, cds_pos: &CdsPos) -> Result<u64, FerroError> {
        // First convert to transcript position
        let tx_pos = self.cds_to_tx(cds_pos)?;

        // If not intronic, use the standard method
        if !cds_pos.is_intronic() {
            return self
                .tx_to_genomic(&tx_pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Position {} not found in exons", tx_pos.base),
                });
        }

        // For intronic positions, use the intronic_to_genomic method
        let offset = cds_pos.offset.ok_or_else(|| FerroError::ConversionError {
            msg: "Expected intronic offset".to_string(),
        })?;

        self.transcript
            .intronic_to_genomic(tx_pos.base as u64, offset)
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Could not convert intronic position {}+{} to genomic",
                    tx_pos.base, offset
                ),
            })
    }

    /// Convert transcript position with intronic offset to genomic position
    pub fn tx_to_genomic_with_intron(&self, tx_pos: &TxPos) -> Result<u64, FerroError> {
        // If not intronic, use the standard method
        if tx_pos.offset.is_none() {
            return self
                .tx_to_genomic(tx_pos)?
                .ok_or_else(|| FerroError::ConversionError {
                    msg: format!("Position {} not found in exons", tx_pos.base),
                });
        }

        let offset = tx_pos.offset.unwrap();

        self.transcript
            .intronic_to_genomic(tx_pos.base as u64, offset)
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!(
                    "Could not convert intronic position {}+{} to genomic",
                    tx_pos.base, offset
                ),
            })
    }

    /// Convert genomic position to CDS position with intronic offset support
    pub fn genomic_to_cds_intronic(&self, genomic_pos: u64) -> Result<CdsPos, FerroError> {
        // First check if it's exonic
        if let Some(cds_pos) = self.genomic_to_cds(genomic_pos)? {
            return Ok(cds_pos);
        }

        // It's intronic - get the transcript boundary and offset
        let (tx_boundary, offset) = self
            .transcript
            .genomic_to_intronic(genomic_pos)
            .ok_or_else(|| FerroError::ConversionError {
                msg: format!("Genomic position {} is not within transcript", genomic_pos),
            })?;

        // Convert transcript boundary to CDS
        let cds_boundary = self.tx_to_cds(&TxPos::new(tx_boundary as i64))?;

        Ok(CdsPos {
            base: cds_boundary.base,
            offset: Some(offset),
            utr3: cds_boundary.utr3,
        })
    }
}

/// Result of genomic coordinate lookup
#[derive(Debug, Clone, PartialEq)]
pub struct GenomicLocation {
    /// Chromosome name
    pub chromosome: String,
    /// Genomic position (1-based, HGVS g. format)
    pub position: u64,
    /// Strand orientation
    pub strand: Strand,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            // 5' UTR (5bp) + CDS (30bp) + 3' UTR (5bp) = 40bp
            sequence: "AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string(),
            cds_start: Some(6),
            cds_end: Some(35),
            exons: vec![Exon::new(1, 1, 38)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_cds_to_tx_normal() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.1 should be tx position 6
        let result = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
        assert_eq!(result.base, 6);
    }

    #[test]
    fn test_cds_to_tx_5utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.-3 should be tx position 3
        let result = mapper.cds_to_tx(&CdsPos::new(-3)).unwrap();
        assert_eq!(result.base, 3);
    }

    #[test]
    fn test_cds_to_tx_3utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.*2 should be tx position 37 (cds_end=35 + 2)
        let result = mapper.cds_to_tx(&CdsPos::utr3(2)).unwrap();
        assert_eq!(result.base, 37);
    }

    #[test]
    fn test_tx_to_cds_normal() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // tx position 10 should be c.5 (10 - 6 + 1)
        let result = mapper.tx_to_cds(&TxPos::new(10)).unwrap();
        assert_eq!(result.base, 5);
        assert!(!result.utr3);
    }

    #[test]
    fn test_tx_to_cds_5utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // tx position 3 should be c.-3
        let result = mapper.tx_to_cds(&TxPos::new(3)).unwrap();
        assert_eq!(result.base, -3);
    }

    #[test]
    fn test_tx_to_cds_3utr() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // tx position 37 should be c.*2
        let result = mapper.tx_to_cds(&TxPos::new(37)).unwrap();
        assert_eq!(result.base, 2);
        assert!(result.utr3);
    }

    #[test]
    fn test_cds_to_protein() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // c.1, c.2, c.3 should all map to p.1
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(1)).unwrap().number, 1);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(2)).unwrap().number, 1);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(3)).unwrap().number, 1);

        // c.4, c.5, c.6 should all map to p.2
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(4)).unwrap().number, 2);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(5)).unwrap().number, 2);
        assert_eq!(mapper.cds_to_protein(&CdsPos::new(6)).unwrap().number, 2);
    }

    #[test]
    fn test_translate_codon() {
        assert_eq!(
            CoordinateMapper::translate_codon("ATG"),
            Some(AminoAcid::Met)
        );
        assert_eq!(
            CoordinateMapper::translate_codon("TAA"),
            Some(AminoAcid::Ter)
        );
        assert_eq!(
            CoordinateMapper::translate_codon("GGG"),
            Some(AminoAcid::Gly)
        );
    }

    fn make_genomic_transcript_plus() -> Transcript {
        // Plus strand transcript with 3 exons
        // Exon 1: tx 1-10, genomic 1000-1009
        // Exon 2: tx 11-20, genomic 2000-2009
        // Exon 3: tx 21-30, genomic 3000-3009
        Transcript {
            id: "NM_GENOMIC.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: "A".repeat(30),
            cds_start: Some(5),
            cds_end: Some(25),
            exons: vec![
                Exon::with_genomic(1, 1, 10, 1000, 1009),
                Exon::with_genomic(2, 11, 20, 2000, 2009),
                Exon::with_genomic(3, 21, 30, 3000, 3009),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(3009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    fn make_genomic_transcript_minus() -> Transcript {
        // Minus strand transcript with 3 exons
        // Exon 1: tx 1-10, genomic 3009-3000 (reversed)
        // Exon 2: tx 11-20, genomic 2009-2000 (reversed)
        // Exon 3: tx 21-30, genomic 1009-1000 (reversed)
        Transcript {
            id: "NM_GENOMIC_MINUS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Minus,
            sequence: "A".repeat(30),
            cds_start: Some(5),
            cds_end: Some(25),
            exons: vec![
                Exon::with_genomic(1, 1, 10, 3000, 3009),
                Exon::with_genomic(2, 11, 20, 2000, 2009),
                Exon::with_genomic(3, 21, 30, 1000, 1009),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(3009),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_genomic_to_tx_plus_strand() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // Exon 1: genomic 1000 -> tx 1
        assert_eq!(mapper.genomic_to_tx(1000).unwrap(), Some(TxPos::new(1)));
        // Exon 1: genomic 1005 -> tx 6
        assert_eq!(mapper.genomic_to_tx(1005).unwrap(), Some(TxPos::new(6)));
        // Exon 2: genomic 2000 -> tx 11
        assert_eq!(mapper.genomic_to_tx(2000).unwrap(), Some(TxPos::new(11)));
        // Exon 3: genomic 3009 -> tx 30
        assert_eq!(mapper.genomic_to_tx(3009).unwrap(), Some(TxPos::new(30)));
    }

    #[test]
    fn test_genomic_to_tx_minus_strand() {
        let tx = make_genomic_transcript_minus();
        let mapper = CoordinateMapper::new(&tx);

        // Exon 1: genomic 3009 -> tx 1 (minus strand, higher genomic = earlier tx)
        assert_eq!(mapper.genomic_to_tx(3009).unwrap(), Some(TxPos::new(1)));
        // Exon 1: genomic 3000 -> tx 10
        assert_eq!(mapper.genomic_to_tx(3000).unwrap(), Some(TxPos::new(10)));
        // Exon 2: genomic 2009 -> tx 11
        assert_eq!(mapper.genomic_to_tx(2009).unwrap(), Some(TxPos::new(11)));
    }

    #[test]
    fn test_genomic_to_tx_intronic() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // Position 1500 is between exon 1 (1009) and exon 2 (2000) - intronic
        assert_eq!(mapper.genomic_to_tx(1500).unwrap(), None);
    }

    #[test]
    fn test_tx_to_genomic_plus_strand() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // tx 1 -> genomic 1000
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(1)).unwrap(), Some(1000));
        // tx 15 -> genomic 2004
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(15)).unwrap(), Some(2004));
    }

    #[test]
    fn test_tx_to_genomic_minus_strand() {
        let tx = make_genomic_transcript_minus();
        let mapper = CoordinateMapper::new(&tx);

        // tx 1 -> genomic 3009 (minus strand)
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(1)).unwrap(), Some(3009));
        // tx 10 -> genomic 3000
        assert_eq!(mapper.tx_to_genomic(&TxPos::new(10)).unwrap(), Some(3000));
    }

    #[test]
    fn test_genomic_to_cds() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // genomic 1004 -> tx 5 -> c.1 (cds_start = 5)
        let cds_pos = mapper.genomic_to_cds(1004).unwrap().unwrap();
        assert_eq!(cds_pos.base, 1);
        assert!(!cds_pos.utr3);
    }

    #[test]
    fn test_cds_to_genomic() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        // c.1 -> tx 5 -> genomic 1004
        assert_eq!(mapper.cds_to_genomic(&CdsPos::new(1)).unwrap(), Some(1004));
    }

    #[test]
    fn test_no_genomic_coords_error() {
        let tx = make_test_transcript(); // No genomic coords
        let mapper = CoordinateMapper::new(&tx);

        assert!(mapper.genomic_to_tx(1000).is_err());
        assert!(mapper.tx_to_genomic(&TxPos::new(1)).is_err());
    }

    /// Create a transcript with gaps between exons (like cdot format)
    /// This simulates how cdot encodes transcripts with virtual intron positions
    fn make_transcript_with_gaps() -> Transcript {
        // Transcript with 3 exons and gaps between them
        // Exon 1: tx 2-94 (93 bases)
        // Gap at tx 95
        // Exon 2: tx 96-193 (98 bases), CDS starts at tx 114
        // Gap at tx 194
        // Exon 3: tx 195-247 (53 bases)
        Transcript {
            id: "NM_GAPS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: "A".repeat(244), // Not used for coordinate tests
            cds_start: Some(114),      // 1-based CDS start
            cds_end: Some(247),        // 1-based CDS end
            exons: vec![
                Exon::new(1, 2, 94),    // Exon 1: tx 2-94
                Exon::new(2, 96, 193),  // Exon 2: tx 96-193 (gap at 95)
                Exon::new(3, 195, 247), // Exon 3: tx 195-247 (gap at 194)
            ],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_cds_to_tx_with_gaps_within_exon() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // CDS 1 = tx 114 (cds_start)
        let result = mapper.cds_to_tx(&CdsPos::new(1)).unwrap();
        assert_eq!(result.base, 114);

        // CDS 10 = tx 123 (still within exon 2)
        let result = mapper.cds_to_tx(&CdsPos::new(10)).unwrap();
        assert_eq!(result.base, 123);
    }

    #[test]
    fn test_cds_to_tx_with_gaps_at_exon_boundary() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // CDS 80 = tx 193 (last base of exon 2)
        // 114 + 79 = 193
        let result = mapper.cds_to_tx(&CdsPos::new(80)).unwrap();
        assert_eq!(result.base, 193);

        // CDS 81 = tx 195 (first base of exon 3, skipping gap at tx 194)
        let result = mapper.cds_to_tx(&CdsPos::new(81)).unwrap();
        assert_eq!(result.base, 195);

        // CDS 82 = tx 196 (second base of exon 3)
        let result = mapper.cds_to_tx(&CdsPos::new(82)).unwrap();
        assert_eq!(result.base, 196);
    }

    #[test]
    fn test_cds_to_tx_with_gaps_crossing_multiple_exons() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // CDS 133 = tx 247 (last base of exon 3)
        // CDS 81-133 is in exon 3 (53 bases)
        // tx 195 + 52 = tx 247
        let result = mapper.cds_to_tx(&CdsPos::new(133)).unwrap();
        assert_eq!(result.base, 247);
    }

    #[test]
    fn test_tx_to_cds_with_gaps_within_exon() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // tx 114 = CDS 1
        let result = mapper.tx_to_cds(&TxPos::new(114)).unwrap();
        assert_eq!(result.base, 1);

        // tx 123 = CDS 10
        let result = mapper.tx_to_cds(&TxPos::new(123)).unwrap();
        assert_eq!(result.base, 10);
    }

    #[test]
    fn test_tx_to_cds_with_gaps_at_exon_boundary() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // tx 193 = CDS 80 (last base of exon 2)
        let result = mapper.tx_to_cds(&TxPos::new(193)).unwrap();
        assert_eq!(result.base, 80);

        // tx 195 = CDS 81 (first base of exon 3)
        let result = mapper.tx_to_cds(&TxPos::new(195)).unwrap();
        assert_eq!(result.base, 81);

        // tx 196 = CDS 82
        let result = mapper.tx_to_cds(&TxPos::new(196)).unwrap();
        assert_eq!(result.base, 82);
    }

    #[test]
    fn test_tx_to_cds_with_gaps_end_of_exon3() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // tx 247 = CDS 133 (last base of exon 3)
        let result = mapper.tx_to_cds(&TxPos::new(247)).unwrap();
        assert_eq!(result.base, 133);
    }

    #[test]
    fn test_cds_tx_roundtrip_with_gaps() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

        // Test roundtrip for various CDS positions
        for cds_pos in [1, 10, 80, 81, 82, 100, 133] {
            let tx_pos = mapper.cds_to_tx(&CdsPos::new(cds_pos)).unwrap();
            let back = mapper.tx_to_cds(&tx_pos).unwrap();
            assert_eq!(
                back.base, cds_pos,
                "Roundtrip failed for CDS {}: got tx {}, back to CDS {}",
                cds_pos, tx_pos.base, back.base
            );
        }
    }
}

#[cfg(test)]
mod intronic_debug_tests {
    use super::*;
    use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
    use std::sync::OnceLock;

    /// Create a transcript mimicking NM_003742.4 exon structure (first 14 exons)
    /// Based on annotation file which uses 0-based half-open coordinates:
    /// CDS: start=127, end=4093
    /// Exon 12: [1324, 1435)
    /// Exon 13: [1435, 1561)
    /// Exon 14: [1561, 1765)
    fn create_nm003742_like_transcript() -> Transcript {
        // After multi_fasta.rs conversion (0-based to 1-based):
        // start = original_start + 1, end = original_end (half-open becomes inclusive)
        let exons = vec![
            Exon {
                number: 1,
                start: 1,
                end: 100,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 2,
                start: 101,
                end: 203,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 3,
                start: 204,
                end: 225,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 4,
                start: 226,
                end: 277,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 5,
                start: 278,
                end: 516,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 6,
                start: 517,
                end: 604,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 7,
                start: 605,
                end: 738,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 8,
                start: 739,
                end: 910,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 9,
                start: 911,
                end: 1035,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 10,
                start: 1036,
                end: 1210,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 11,
                start: 1211,
                end: 1324,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 12,
                start: 1325,
                end: 1435,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 13,
                start: 1436,
                end: 1561,
                genomic_start: None,
                genomic_end: None,
            },
            Exon {
                number: 14,
                start: 1562,
                end: 1765,
                genomic_start: None,
                genomic_end: None,
            },
        ];

        Transcript {
            id: "NM_003742.4".to_string(),
            gene_symbol: Some("ABCB11".to_string()),
            strand: Strand::Minus,
            sequence: String::new(), // Not needed for coordinate tests
            cds_start: Some(128),    // 127 + 1 = 128 (1-based)
            cds_end: Some(4093),     // Same (half-open end = inclusive 1-based)
            exons,
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn debug_nm003742_intronic() {
        let transcript = create_nm003742_like_transcript();

        eprintln!("\n=== NM_003742.4-like Transcript ===");
        eprintln!("cds_start: {:?}", transcript.cds_start);
        eprintln!("cds_end: {:?}", transcript.cds_end);
        eprintln!("num exons: {}", transcript.exons.len());

        // Check for gaps
        let mut sorted_exons: Vec<_> = transcript.exons.iter().collect();
        sorted_exons.sort_by_key(|e| e.start);

        eprintln!("\nExons 11-14:");
        for e in &sorted_exons[10..14] {
            eprintln!("  Exon {}: tx {}..{}", e.number, e.start, e.end);
        }

        let has_gaps = sorted_exons.windows(2).any(|w| w[0].end + 1 != w[1].start);
        eprintln!("\nHas gaps in tx coords: {}", has_gaps);

        // Calculate expected tx for c.1435
        let cds_start = transcript.cds_start.unwrap();
        let expected_tx = cds_start as i64 + 1434;
        eprintln!("\nFor c.1435:");
        eprintln!("  cds_start = {}", cds_start);
        eprintln!("  expected tx = {} + 1434 = {}", cds_start, expected_tx);

        // Use the mapper
        let mapper = CoordinateMapper::new(&transcript);
        let cds_pos = CdsPos {
            base: 1435,
            offset: None,
            utr3: false,
        };
        let tx_result = mapper.cds_to_tx(&cds_pos);
        eprintln!("  actual tx from mapper: {:?}", tx_result);

        // Verify the calculation
        assert!(tx_result.is_ok());
        let tx_pos = tx_result.unwrap();

        // c.1435 should map to tx 128 + 1434 = 1562
        assert_eq!(
            tx_pos.base, 1562,
            "c.1435 should map to tx 1562, but got tx {}",
            tx_pos.base
        );
    }

    /// Test that cdot tx coordinates are correctly interpreted as 1-based.
    ///
    /// The cdot format uses 1-based tx_start and tx_end coordinates (the first exon
    /// starts at position 1, not 0). This test verifies that assumption.
    #[test]
    fn test_cdot_tx_coordinates_are_1_based() {
        use crate::data::cdot::CdotMapper;
        use std::path::PathBuf;

        // Load the real cdot data if available
        let cdot_path = PathBuf::from("benchmark-output/cdot/cdot-0.2.32.refseq.GRCh38.json");
        if !cdot_path.exists() {
            // Skip test if cdot file is not available
            return;
        }

        let cdot = CdotMapper::from_json_file(&cdot_path).expect("Failed to load cdot");
        let tx = cdot
            .get_transcript("NM_003742.4")
            .expect("NM_003742.4 not in cdot");

        // First exon should have tx_start = 1 (1-based, not 0-based)
        let first_exon = &tx.exons[0];
        assert_eq!(
            first_exon[2], 1,
            "cdot tx_start for first exon should be 1 (1-based), got {}",
            first_exon[2]
        );

        // CDS start should be 0-based (127 for NM_003742.4, which becomes 128 when converted to 1-based)
        assert_eq!(
            tx.cds_start,
            Some(127),
            "cdot cds_start should be 0-based (127 for NM_003742.4)"
        );
    }
}
