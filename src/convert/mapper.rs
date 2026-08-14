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
//! # The two frames, and which one `c.`/`n.` lives in
//!
//! The conversions above split into two frames that must not borrow each
//! other's arithmetic. Getting that wrong is #1619.
//!
//! **The sequence frame — `cds_to_tx` / `tx_to_cds`.** A `c.` or `n.` position
//! names a base of the transcript reference sequence and nothing else.
//! `docs/background/numbering.md:21` anchors `c.1` on "the **`A`** of the `ATG`
//! translation initiation (start) codon"; `:52` numbers `n.` "from the first to
//! the last nucleotide of the reference sequence"; and `:40` states that a
//! coding DNA reference sequence "**does not contain** intron or 5' and 3' gene
//! flanking sequences". So the axis is a plain count of the transcript's own
//! bases, and `c.N` is `cds_start + N - 1` on the flat sequence — the same
//! bytes `ReferenceProvider::get_sequence` serves and the normalizer indexes.
//! **These two methods therefore do not read the exon table at all.**
//!
//! Two shapes have each tried to make them read it, and both were wrong:
//!
//! - a CIGAR insertion, i.e. transcript bases absent from the genome. `c.`
//!   counts them, so no adjustment (#944 — `NM_015120.4` `c.87` is `n.198`,
//!   not `n.201`).
//! - a transcript-coordinate gap, i.e. transcript bases that align to no exon.
//!   `c.` counts those too (#1619 — `NM_033517.1` carries a real 39-base cdot
//!   hole, and RefSeq annotates its CDS as the contiguous `1..5196` across it,
//!   matching the 1731-aa `NP_277052.1`). Walking the exon list moved every
//!   `c.` position on that accession, so `hgvs_to_spdi` and the normalizer
//!   named different bases of the same transcript.
//!
//! **The genome frame — `genomic_to_tx` / `tx_to_genomic` and everything built
//! on them.** Mapping between a contig and a transcript genuinely is exon- and
//! CIGAR-aware, and stays so. That is where the exon table belongs, and where a
//! gap correctly means "this transcript base aligns to nothing" rather than
//! "renumber the positions after it".
//!
//! The composed conversions (`genomic_to_cds`, `cds_to_genomic`) are one step
//! of each: exon-aware across genome↔tx, flat across tx↔CDS.
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
    /// **This is a sequence-frame conversion and it does not consult the exon
    /// table** — see the module's "two frames" note. `c.`/`n.` numbering counts
    /// the transcript's own bases, so
    /// `c.N` is `cds_start + N - 1` on the flat transcript, with no adjustment
    /// for exon geometry, for a CIGAR insertion (#944), or for a
    /// transcript-coordinate gap in the genome alignment (#1619).
    ///
    /// One input where this is **not** interchangeable with the normalizer's
    /// own sequence-frame copy (`Normalizer::cds_to_tx_pos`): `c.0`, which is
    /// not a valid HGVS position but is reachable (`CDS_BASE_UNKNOWN` is `0`,
    /// and cis-collapse can build a `base == 0` anchor). It falls into the
    /// `base < 1` branch here and answers `cds_start`, where the normalizer
    /// continues the formula and answers `cds_start - 1`. Which is right needs
    /// its own ruling; see #1772. A `c.-N` further upstream than the transcript
    /// is long diverges too: this method returns a `TxPos` with a *negative*
    /// `base` (which [`Self::fold_non_intronic_utr_offset`] relies on), where
    /// the normalizer's copy errors.
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
            // 3' UTR: *N = cds_end + N
            cds_end as i64 + pos.base
        } else if pos.base < 1 {
            // 5' UTR: -N = cds_start - N (N is positive in HGVS, stored negative),
            // e.g. c.-3 is 3 bases before the CDS start.
            cds_start as i64 + pos.base
        } else {
            // Normal CDS position: c.1 = cds_start, c.2 = cds_start + 1, …
            cds_start as i64 + pos.base - 1
        };

        Ok(TxPos {
            base: tx_base,
            offset: pos.offset,
            downstream: false,
        })
    }

    /// Convert transcript position to CDS position
    ///
    /// The inverse of [`Self::cds_to_tx`], and on the same flat sequence axis:
    /// no exon walk, no CIGAR adjustment. The one input on which the round trip
    /// does not close is `c.0` — it converts to `cds_start` and comes back as
    /// `c.1`, because HGVS has no `c.0` for it to come back as.
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
            // 5' UTR: `c.-N`, N bases 5' of the CDS start.
            Ok(CdsPos {
                base: base - cds_start, // negative for 5' UTR
                offset: pos.offset,
                utr3: false,
                special: None,
            })
        } else if base > cds_end {
            // 3' UTR: `c.*N`, N bases 3' of the last CDS base.
            Ok(CdsPos {
                base: base - cds_end,
                offset: pos.offset,
                utr3: true,
                special: None,
            })
        } else {
            // Within the CDS.
            Ok(CdsPos {
                base: base - cds_start + 1,
                offset: pos.offset,
                utr3: false,
                special: None,
            })
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

        // Need cached transcript bases to translate; coordinate-only
        // transcripts return `None`.
        let seq = self.transcript.sequence.as_deref()?;
        if codon_end > cds_end as usize || codon_end > seq.len() {
            return None;
        }

        let codon = &seq[codon_start..codon_end];
        Self::translate_codon(codon)
    }

    /// Translate a codon to an amino acid.
    ///
    /// Total: every one of the 22 match arms (21 named plus the `_`
    /// fall-through) returns `Some`, so this never returns `None` despite
    /// the `Option` signature. The `_` arm is what makes it total: it
    /// coerces *every* unrecognised codon — including `N`-containing ones
    /// from real reference sequences — to `Some(Xaa)`, indistinguishable
    /// from a genuine unknown residue.
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
    ///
    /// **This is a genome-frame conversion and it does consult the exon
    /// table** — see the module's "two frames" note, and contrast
    /// [`Self::cds_to_tx`], which must not. The two answer different questions
    /// and therefore *must* disagree on a transcript whose alignment carries a
    /// transcript-coordinate gap; `tests/it/axis_frame_disagreement.rs` pins
    /// that, so collapsing them onto one implementation fails loudly rather
    /// than reintroducing #1619.
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
                    Strand::Unknown => return Ok(None),
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
    ///
    /// **Genome frame**, like [`Self::genomic_to_tx`] and unlike
    /// [`Self::tx_to_cds`]: the alignment answer of `None` means "this
    /// transcript base aligns to nothing", which is a fact about the alignment
    /// and not about the transcript's own numbering — that base still has a
    /// `c.`/`n.` number.
    ///
    /// `None` is not exclusively that answer, though, so do not read it as
    /// "unaligned" without ruling the other two out: it is also returned for a
    /// non-positive `tx_pos.base`, and for a `Strand::Unknown` transcript,
    /// whose bases are all reported unaligned however complete its exon table
    /// is. A caller using `tx_to_genomic(..).is_some()` as its definition of
    /// "aligned" — `tests/it/axis_frame_disagreement.rs`'s genome-frame walk
    /// does — is relying on the transcript having a known strand.
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
                    Strand::Unknown => return Ok(None),
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

    /// Fold a UTR offset that has *no enclosing intron* into a plain
    /// (non-intronic) CDS coordinate.
    ///
    /// A `+offset`/`-offset` is only meaningful when an intron spans the
    /// boundary. In the UTRs that intron can be absent for two reasons:
    /// - the boundary is the transcript terminus, so the offset points into the
    ///   contiguous genome past the final/first exon (`c.*824+10` 3' of the last
    ///   exon → `c.*834`); or
    /// - the transcript model simply has no intron there (e.g. a single-exon
    ///   supplemental transcript), so the offset is a spurious decoration.
    ///
    /// In both cases the offset is not intronic and folds linearly into the
    /// base coordinate, matching the canonical interpretation and letting the
    /// standard (non-intronic) machinery map and normalize the position instead
    /// of the intron-window shuffle that has no intron to anchor on (#760).
    ///
    /// # Returns
    /// The folded position, or `None` when it should keep its existing
    /// handling: no offset, a position in the CDS proper (where a missing
    /// intron means a genuinely invalid description, not an extended
    /// coordinate), or a genuine enclosing intron.
    pub fn fold_non_intronic_utr_offset(&self, pos: &CdsPos) -> Option<CdsPos> {
        let offset = pos.offset?;
        if offset == 0 {
            return None;
        }

        // Only UTR offsets fold. A mid-CDS offset with no intron is invalid, not
        // an extended coordinate, so it keeps erroring rather than folding.
        let in_utr = pos.utr3 || pos.base < 0;
        if !in_utr {
            return None;
        }

        let tx = self.cds_to_tx(pos).ok()?;

        // A negative tx base sits upstream of transcript base 1 (a deep 5'UTR
        // offset, e.g. `c.-N` with N past the transcript start), so no interior
        // intron can enclose it. Skip the intron probe — which is `u64`-indexed
        // and cannot represent a pre-transcript base — and fold linearly below.
        // A genuine enclosing intron only applies for a non-negative tx base.
        if tx.base >= 0 {
            let tx_base = tx.base as u64;
            // A genuine intron encloses the offset → leave it to the intronic path.
            if self
                .transcript
                .find_intron_at_tx_boundary(tx_base, offset)
                .is_some()
            {
                return None;
            }
        }

        // Transcript coordinates run 5'->3' regardless of strand, so a signed
        // offset shifts the boundary directly; the (contiguous) sequence past
        // the boundary maps back through the standard CDS conversion, which
        // accepts coordinates beyond the transcript termini.
        let folded_tx = tx.base.checked_add(offset)?;
        self.tx_to_cds(&TxPos::new(folded_tx)).ok()
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
            special: None,
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
    use crate::data::CigarOp;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            // 5'UTR (5bp) + CDS (30bp) + 3'UTR (3bp) = 38bp.
            // NOTE: cds_start = 6 is the `T` of the visible `ATG` (whose `A` is at
            // 0-based 4), so the first codon is seq[5..8] == "TGC" and p.1 is Cys,
            // not Met. Use `make_codon_sweep_transcript` for amino-acid assertions.
            sequence: Some("AAAAATGCCCAAAGGGTTTAGGCCCAAAGGGTTATAAA".to_string()),
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
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    /// One codon per `translate_codon` match arm, 21 codons, all 21 amino acids
    /// distinct. Three properties are load-bearing — do not "tidy" them:
    ///
    /// 1. All 21 amino acids differ, so a mutated `codon_start` reads a codon
    ///    that translates differently rather than coinciding.
    /// 2. The 5'UTR's first codon is `CCC` (Pro), never `ATG`. In
    ///    `get_amino_acid_at`'s `codon_start` arithmetic
    ///    (`cds_start as usize - 1 + (position as usize - 1) * 3`), a
    ///    `+ -> *` mutant degenerates `codon_start` to 0 at p.1 and reads
    ///    seq[0..3]; a Met-translating 5'UTR would let it survive.
    /// 3. The 3'UTR translates (`GGG` -> Gly). In `get_amino_acid_at`'s
    ///    past-CDS guard (`codon_end > cds_end as usize || codon_end >
    ///    seq.len()`), a `|| -> &&` mutant stops it firing at p.22 and reads
    ///    seq[68..71]; Gly there is what distinguishes it from the `Xaa` the
    ///    original returns.
    ///
    /// Layout: 5'UTR 5bp (cds_start = 6), CDS 63bp (cds_end = 68), 3'UTR 6bp,
    /// total 74bp.
    fn make_codon_sweep_transcript() -> Transcript {
        Transcript {
            id: "NM_CODONS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            sequence: Some(
                "CCCCC".to_string()
                    + "ATGTTTCTGATTGTGTCACCCACAGCCTAT"
                    + "CATCAAAATAAAGATGAATGTTGGCGCGGA"
                    + "TAA"
                    + "GGGTTT",
            ),
            cds_start: Some(6),
            cds_end: Some(68),
            exons: vec![Exon::new(1, 1, 74)],
            ..Default::default()
        }
    }

    /// A deliberately malformed transcript whose `cds_end` exceeds its sequence
    /// length. This is the ONLY way to reach `get_amino_acid_at`'s second guard
    /// clause: under `||`, `codon_end > seq.len()` is unreachable whenever
    /// `codon_end > cds_end` has already fired.
    ///
    /// Sequence is 12bp, cds_start = 1, cds_end = 100.
    ///
    /// Do not reuse this fixture with `validate_cds_pos`: its
    /// `transcript.sequence_length() - cds_end` is unchecked `u64`
    /// subtraction, and `12 - 100` here underflows and panics in debug.
    fn make_cds_end_past_sequence_transcript() -> Transcript {
        Transcript {
            id: "NM_SHORTSEQ.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            sequence: Some("ATGTTTCTGATT".to_string()),
            cds_start: Some(1),
            cds_end: Some(100),
            exons: vec![Exon::new(1, 1, 12)],
            ..Default::default()
        }
    }

    /// Pins the amino acid at every codon of a 21-arm sweep.
    ///
    /// `test_cds_to_protein` asserts only `.number`, leaving the amino acid
    /// itself unchecked: `get_amino_acid_at` returning `None` is invisible
    /// there because `cds_to_protein` swallows it with `.unwrap_or(Xaa)`,
    /// and every `translate_codon` arm falls through to `_ => Some(Xaa)`
    /// when deleted, so a mutated arm and a correct one are indistinguishable
    /// by number alone. This test checks `.aa` directly so both are caught.
    #[test]
    fn cds_to_protein_pins_every_codon_arm() {
        let tx = make_codon_sweep_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let expected = [
            AminoAcid::Met,
            AminoAcid::Phe,
            AminoAcid::Leu,
            AminoAcid::Ile,
            AminoAcid::Val,
            AminoAcid::Ser,
            AminoAcid::Pro,
            AminoAcid::Thr,
            AminoAcid::Ala,
            AminoAcid::Tyr,
            AminoAcid::His,
            AminoAcid::Gln,
            AminoAcid::Asn,
            AminoAcid::Lys,
            AminoAcid::Asp,
            AminoAcid::Glu,
            AminoAcid::Cys,
            AminoAcid::Trp,
            AminoAcid::Arg,
            AminoAcid::Gly,
            AminoAcid::Ter,
        ];

        for (index, aa) in expected.iter().enumerate() {
            let protein_number = (index + 1) as u64;
            let cds_base = (index as i64) * 3 + 1;
            let got = mapper.cds_to_protein(&CdsPos::new(cds_base)).unwrap();
            assert_eq!(got.number, protein_number, "p.{} number", protein_number);
            assert_eq!(
                got.aa, *aa,
                "p.{} amino acid (c.{})",
                protein_number, cds_base
            );
        }
    }

    /// p.21 is the last codon (`codon_end == cds_end`, 68) and must translate;
    /// p.22 runs past the CDS and must yield `Xaa`.
    ///
    /// The p.22 case is what kills a `|| -> &&` mutant of `get_amino_acid_at`'s
    /// guard (the `codon_end > cds_end as usize || codon_end > seq.len()`
    /// early return in `get_amino_acid_at`): with `&&` the guard does
    /// not fire, the function reads the 3'UTR's `GGG`, and `cds_to_protein`
    /// returns Gly instead of Xaa.
    ///
    /// The two indirect assertions above observe the guard only through
    /// `cds_to_protein`'s `.unwrap_or(AminoAcid::Xaa)`, which — exactly as
    /// `get_amino_acid_at_guards_a_cds_end_past_the_sequence`'s doc comment
    /// diagnoses for its own case — collapses five distinct causes of
    /// `None` into one `Xaa` value: a missing `cds_start`, a missing
    /// `cds_end`, a missing cached `sequence`, `codon_end > cds_end` (the
    /// primary guard clause), and `codon_end > seq.len()` (the secondary
    /// one). `translate_codon` is not a sixth cause — it is total (every
    /// match arm, including `_`, yields `Some`, coercing any unrecognised
    /// codon to `Xaa` rather than `None`), so it can never contribute a
    /// `None` here. That test unmasks the *secondary*,
    /// malformed-transcript-only guard clause (`codon_end > seq.len()`) by
    /// calling `get_amino_acid_at` directly; this test exercises the
    /// *primary* clause (`codon_end > cds_end as usize`), the one that fires
    /// on a well-formed transcript, so it gets the same direct treatment
    /// here.
    #[test]
    fn cds_to_protein_stops_at_the_cds_end() {
        let tx = make_codon_sweep_transcript();
        let mapper = CoordinateMapper::new(&tx);

        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(61)).unwrap().aa,
            AminoAcid::Ter,
            "c.61 is the last codon and must still translate"
        );
        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(64)).unwrap().aa,
            AminoAcid::Xaa,
            "c.64 is past cds_end and must not read into the 3'UTR"
        );

        // Direct, unmasked: `get_amino_acid_at` itself must draw the p.21/p.22
        // line at `codon_end > cds_end`, not just leave `cds_to_protein` to
        // paper over it with `Xaa`.
        assert_eq!(
            mapper.get_amino_acid_at(21),
            Some(AminoAcid::Ter),
            "get_amino_acid_at(21) must translate unmasked, codon_end 68 == cds_end 68"
        );
        assert_eq!(
            mapper.get_amino_acid_at(22),
            None,
            "get_amino_acid_at(22) must return None unmasked, codon_end 71 > cds_end 68"
        );
    }

    /// Exercises `get_amino_acid_at`'s SECOND guard clause
    /// (`codon_end > seq.len()`), which is unreachable on a well-formed
    /// transcript because the first clause fires first.
    ///
    /// cds_end = 100 but the sequence is 12bp, so:
    ///   p.1 codon_end = 3  (< 12)  -> translates
    ///   p.4 codon_end = 12 (== 12) -> translates, the boundary
    ///   p.5 codon_end = 15 (> 12)  -> Xaa
    ///
    /// The assertions below observe the guard only through `cds_to_protein`,
    /// which coerces every `None` from `get_amino_acid_at` — regardless of
    /// which of at least five distinct causes produced it — into `Xaa` via
    /// `cds_to_protein`'s `.unwrap_or(AminoAcid::Xaa)` fallback. That masks
    /// the guard: the assertion cannot fail for the *right* reason. The
    /// direct call below observes `get_amino_acid_at` unmasked, alongside
    /// the existing indirect ones.
    #[test]
    fn get_amino_acid_at_guards_a_cds_end_past_the_sequence() {
        let tx = make_cds_end_past_sequence_transcript();
        let mapper = CoordinateMapper::new(&tx);

        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(1)).unwrap().aa,
            AminoAcid::Met,
            "codon_end 3 is well inside a 12bp sequence"
        );
        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(10)).unwrap().aa,
            AminoAcid::Ile,
            "codon_end 12 == seq.len() is the last readable codon"
        );
        assert_eq!(
            mapper.cds_to_protein(&CdsPos::new(13)).unwrap().aa,
            AminoAcid::Xaa,
            "codon_end 15 > seq.len() must not read out of bounds"
        );

        // Direct, unmasked: `get_amino_acid_at` itself must return `None`
        // past the end (p.5, codon_end 15 > seq.len() 12), not just leave
        // `cds_to_protein` to paper over it with `Xaa`.
        assert_eq!(
            mapper.get_amino_acid_at(5),
            None,
            "get_amino_acid_at(5) must return None unmasked, codon_end 15 > seq.len() 12"
        );
    }

    #[test]
    fn test_cds_to_tx_normal() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        // Pins `make_test_transcript`'s length (5'UTR 5bp + CDS 30bp + 3'UTR
        // 3bp = 38bp) against the hardcoded `38` literal below — that catches
        // the sequence literal changing size, not the comment drifting from
        // it (both copies once said `= 40bp` while the sequence stayed 38bp,
        // and this assertion would not have caught that). The fixture is
        // duplicated with the same data but different code in
        // `tests/it/convert_tests.rs` (a `Transcript { .. }` literal here vs.
        // `Transcript::new(..)` there). It is also not sequence-specific:
        // `sequence_length()` falls back to summing exon lengths when
        // `sequence` is `None`, and this fixture's `Exon::new(1, 1, 38)`
        // sums to 38 too, so dropping the sequence would still pass.
        assert_eq!(tx.sequence_length(), 38, "make_test_transcript is 38bp");

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

    /// `c.*3+5` on a single-exon transcript whose 3'UTR ends at `*3` has no
    /// enclosing intron — it is the extended 3'UTR coordinate `c.*8`. Folding
    /// lets the standard (non-intronic) machinery map it. Regression for #760.
    #[test]
    fn fold_three_prime_utr_offset_at_terminus() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // cds_end=35, single exon ends at tx 38 → *3 is the final transcript base.
        let pos = CdsPos {
            base: 3,
            offset: Some(5),
            utr3: true,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, 8);
        assert!(folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// A 3'UTR offset *interior* to a single-exon transcript (no intron there at
    /// all, the degenerate supplemental-transcript shape behind the #760 corpus
    /// row) still folds: `c.*2+1` → `c.*3`.
    #[test]
    fn fold_three_prime_utr_offset_interior() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // *2 maps to tx 37, interior to the single exon (which ends at tx 38).
        let pos = CdsPos {
            base: 2,
            offset: Some(1),
            utr3: true,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, 3);
        assert!(folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// `c.-5-3` anchored at the first transcript base extends 5' past the
    /// terminus to `c.-8`.
    #[test]
    fn fold_five_prime_utr_offset_at_terminus() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // cds_start=6 → c.-5 maps to tx 1, the first transcript base.
        let pos = CdsPos {
            base: -5,
            offset: Some(-3),
            utr3: false,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, -8);
        assert!(!folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// A deep 5'UTR offset whose anchor base already maps upstream of
    /// transcript base 1 (a negative tx base) must still fold linearly rather
    /// than short-circuit into the intronic path: `c.-7-2` (tx -1) → `c.-9`.
    /// Regression for the negative-tx-base 5'UTR edge of #760.
    #[test]
    fn fold_five_prime_utr_offset_negative_tx_base() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // cds_start=6 → c.-7 maps to tx -1, upstream of transcript base 1.
        let pos = CdsPos {
            base: -7,
            offset: Some(-2),
            utr3: false,
            special: None,
        };
        let folded = mapper
            .fold_non_intronic_utr_offset(&pos)
            .expect("should fold");
        assert_eq!(folded.base, -9);
        assert!(!folded.utr3);
        assert_eq!(folded.offset, None);
    }

    /// A mid-CDS offset with no intron is a genuinely invalid description, not
    /// an extended coordinate — it must not be folded.
    #[test]
    fn fold_offset_declines_cds_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // c.10+5 maps to tx 15, inside the CDS — not a UTR offset.
        let pos = CdsPos {
            base: 10,
            offset: Some(5),
            utr3: false,
            special: None,
        };
        assert!(mapper.fold_non_intronic_utr_offset(&pos).is_none());
    }

    /// A non-intronic position (no offset) is never folded.
    #[test]
    fn fold_offset_declines_non_intronic_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);
        assert!(mapper
            .fold_non_intronic_utr_offset(&CdsPos::utr3(2))
            .is_none());
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
        // An unrecognised codon falls through to the `_ => Some(Xaa)` arm.
        assert_eq!(
            CoordinateMapper::translate_codon("NNN"),
            Some(AminoAcid::Xaa)
        );
    }

    fn make_genomic_transcript_plus() -> Transcript {
        // Plus strand transcript with 3 exons
        // Exon 1: tx 1-10, genomic 1000-1009
        // Exon 2: tx 11-20, genomic 2000-2009
        // Exon 3: tx 21-30, genomic 3000-3009
        Transcript {
            cds_start_incomplete: false,
            id: "NM_GENOMIC.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(30)),
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
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    fn make_genomic_transcript_minus() -> Transcript {
        // Minus strand transcript with 3 exons
        // Exon 1: tx 1-10, genomic 3009-3000 (reversed)
        // Exon 2: tx 11-20, genomic 2009-2000 (reversed)
        // Exon 3: tx 21-30, genomic 1009-1000 (reversed)
        Transcript {
            cds_start_incomplete: false,
            id: "NM_GENOMIC_MINUS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Minus,
            sequence: Some("A".repeat(30)),
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
            protein_id: None,
            exon_cigars: Vec::new(),
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

    /// A transcript whose exons do not tile transcript space: transcript
    /// positions 1, 95 and 194 fall in no exon.
    ///
    /// Real cdot data does carry such holes — 58 of 474,818 GRCh38 multi-exon
    /// builds (0.012%), the smallest of them 23 bases and none of them one base
    /// — so the shape is not itself malformed. The **sizes** here are synthetic:
    /// one base per junction, which is a shape upstream never produces. That is
    /// deliberate and is now the point of the helper. Under the flat sequence
    /// axis the gap size cannot affect a `c.`/`n.` answer at all, so a hole of
    /// the least plausible size is the sharpest available check that nothing
    /// reads the exon table on this path.
    ///
    /// It does **not** describe cdot's *encoding*. This helper previously
    /// claimed to "simulate how cdot encodes transcripts with virtual intron
    /// positions"; cdot has no such encoding. Its `cds_start_i`/`cds_end_i` are
    /// 1-based inclusive coordinates on the *spliced* transcript, and a hole
    /// means transcript bases that align to no genomic exon — not extra
    /// coordinates that `c.` numbering skips.
    ///
    /// The false claim was load-bearing: it was the stated justification for the
    /// one-base-per-junction holes that
    /// `tests/fixtures/sequences/normalization_transcripts.json` carried in every
    /// multi-exon record, which drifted `NM_000492.4:c.1520` ten bases 3' (#1619).
    /// `tests/it/normalization_transcripts_exon_contract.rs` now guards that
    /// fixture; the three tests below use this helper to pin that the walk those
    /// holes used to trigger is gone.
    fn make_transcript_with_gaps() -> Transcript {
        // Transcript with 3 exons and gaps between them
        // Exon 1: tx 2-94 (93 bases)
        // Gap at tx 95
        // Exon 2: tx 96-193 (98 bases), CDS starts at tx 114
        // Gap at tx 194
        // Exon 3: tx 195-247 (53 bases)
        Transcript {
            cds_start_incomplete: false,
            id: "NM_GAPS.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(244)), // Not used for coordinate tests
            cds_start: Some(114),            // 1-based CDS start
            cds_end: Some(247),              // 1-based CDS end
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
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    /// #1619: a `c.` position is `cds_start + N - 1` on the flat transcript.
    /// The exon holes at tx 95 and 194 must change nothing — `c.81` is tx 194,
    /// the very base that lies in a hole, because that base is in the
    /// transcript sequence and `c.` numbering counts it.
    #[test]
    fn cds_to_tx_ignores_a_transcript_coordinate_gap() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);
        // cds_start = 114.
        for (cds, expected_tx) in [
            (1, 114),
            (10, 123),
            (80, 193),
            (81, 194),
            (82, 195),
            (133, 246),
        ] {
            assert_eq!(
                mapper.cds_to_tx(&CdsPos::new(cds)).unwrap().base,
                expected_tx,
                "c.{cds} must be the flat offset; the exon walk answered {}",
                expected_tx + 1,
            );
        }
    }

    /// The inverse, on the same axis. tx 194 is a real transcript base and so
    /// has a `c.` number (`c.81`); it is only the *genome* that has nothing
    /// there.
    #[test]
    fn tx_to_cds_ignores_a_transcript_coordinate_gap() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);
        for (tx_base, expected_cds) in [
            (114, 1),
            (123, 10),
            (193, 80),
            (194, 81),
            (195, 82),
            (247, 134),
        ] {
            assert_eq!(
                mapper.tx_to_cds(&TxPos::new(tx_base)).unwrap().base,
                expected_cds,
                "n.{tx_base} must be the flat CDS offset",
            );
        }
    }

    #[test]
    fn cds_tx_round_trip_across_a_transcript_coordinate_gap() {
        let tx = make_transcript_with_gaps();
        let mapper = CoordinateMapper::new(&tx);

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

    fn make_cigar_insertion_transcript() -> Transcript {
        Transcript {
            cds_start_incomplete: false,
            id: "NM_INS.1".to_string(),
            gene_symbol: Some("INS".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(40)),
            cds_start: Some(1),
            cds_end: Some(40),
            exons: vec![Exon::new(1, 1, 40)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            // M5 I3 M32: 3 inserted transcript bases at tx 6,7,8.
            exon_cigars: vec![Some(vec![
                CigarOp::Match(5),
                CigarOp::Insertion(3),
                CigarOp::Match(32),
            ])],
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn cds_to_tx_is_transcript_native_across_cigar_insertion() {
        // #944: c./n. numbering counts inserted transcript bases, so CDS->tx is
        // naive with NO CIGAR adjustment. cds_start=1 => c.N -> tx N (1-based).
        let tx = make_cigar_insertion_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // Before the insertion — unaffected.
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(3)).unwrap().base, 3);
        // After the 3bp insertion — must NOT shift by +3 (bug produced 13).
        assert_eq!(mapper.cds_to_tx(&CdsPos::new(10)).unwrap().base, 10);
    }

    #[test]
    fn cds_tx_round_trip_across_cigar_insertion() {
        // tx_to_cds must invert cds_to_tx with no CIGAR adjustment (#944).
        let tx = make_cigar_insertion_transcript();
        let mapper = CoordinateMapper::new(&tx);
        // Round-trip a post-insertion position: c.10 -> tx 10 -> c.10.
        let tx_pos = mapper.cds_to_tx(&CdsPos::new(10)).unwrap();
        assert_eq!(tx_pos.base, 10);
        assert_eq!(mapper.tx_to_cds(&tx_pos).unwrap().base, 10);
        // Direct n->c: tx 10 -> c.10 (bug produced c.7).
        assert_eq!(mapper.tx_to_cds(&TxPos::new(10)).unwrap().base, 10);
    }

    /// `chromosome()` had no assertion anywhere, so `-> None`, `-> Some("")`
    /// and `-> Some("xyzzy")` all survived.
    #[test]
    fn chromosome_reports_the_transcript_contig() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);
        assert_eq!(mapper.chromosome(), Some("chr1"));
    }

    /// `get_intron_position` has one production caller (in `src/vcf/to_hgvs.rs`)
    /// and no test reached it, so the whole-function `-> None` survived.
    ///
    /// Exon 1 ends at genomic 1009 (tx 10); genomic 1012 is 3 bases into
    /// intron 1.
    #[test]
    fn get_intron_position_reports_an_intronic_offset() {
        let tx = make_genomic_transcript_plus();
        let mapper = CoordinateMapper::new(&tx);

        let intron = mapper
            .get_intron_position(1012)
            .expect("genomic 1012 lies in intron 1");
        assert_eq!(intron.intron_number, 1, "genomic 1012 lies in intron 1");
        assert_eq!(intron.offset, 3, "3 bases past the exon 1 boundary");
        assert_eq!(intron.tx_boundary_pos, 10, "exon 1 ends at tx 10");

        // An exonic position is not intronic.
        assert!(mapper.get_intron_position(1005).is_none());
    }

    /// `fold_non_intronic_utr_offset` folds only 3'UTR and negative-base
    /// (5'UTR) positions. `base == 0` is neither of those — it is
    /// `CDS_BASE_UNKNOWN` (in `src/hgvs/location.rs`), the `c.?`
    /// unknown-position sentinel, not an arithmetic edge of the UTR check.
    /// It must still decline to fold, and doing so also happens to pin
    /// `pos.base < 0` against a `<= 0` mutant, since 0 is where the two
    /// diverge. Base 0 is the interesting input precisely because the two
    /// functions disagree about it: this gate classifies 5'UTR as
    /// `pos.base < 0`, while `cds_to_tx` (above) uses `pos.base < 1`,
    /// so base 0 is 5'UTR to `cds_to_tx` but not to this gate — that
    /// divergence is the whole reason this test exists.
    #[test]
    fn fold_offset_declines_unknown_position() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let pos = CdsPos::unknown(Some(5));

        // The decline happens at the `in_utr` gate, three lines before
        // `cds_to_tx` is even called: `pos.utr3 || pos.base < 0` is
        // `false || (0 < 0)` = `false` for base 0 (`CDS_BASE_UNKNOWN`), so
        // `fold_non_intronic_utr_offset` returns `None` right there and the
        // `c.?` sentinel never reaches the conversion at all.
        // (`cds_to_tx`'s own handling of this sentinel — it coerces rather
        // than declines — is pinned separately by
        // `cds_to_tx_coerces_the_unknown_sentinel` below.)
        assert!(mapper.fold_non_intronic_utr_offset(&pos).is_none());
    }

    /// `cds_to_tx` does not honour the `c.?` sentinel the way
    /// `merge::simple_cds_pos` (in `src/normalize/merge.rs`) and
    /// `VariantProjector::project_to_genomic` (in `src/project/projector.rs`)
    /// do — both decline it explicitly. Here it silently coerces
    /// `CdsPos::unknown(Some(5))` (base `CDS_BASE_UNKNOWN` = 0) into
    /// `TxPos { base: 6, offset: Some(5) }`, bit-identical to what `c.1+5`
    /// produces (`cds_start` is 6 on `make_test_transcript`).
    ///
    /// The quirk is not contained by those declining callers:
    /// `spdi::convert::resolve_cds_to_tx` (in `src/spdi/convert.rs`) guards
    /// only `is_intronic()`, which is `false` for `c.?` (offset is `None`),
    /// so the sentinel passes straight through and is emitted as a concrete
    /// SPDI position for what is, by definition, an unknown coordinate.
    /// This pins today's behaviour, not a requirement — if `cds_to_tx`
    /// later grows an `is_unknown()` guard and declines this input, update
    /// this assertion to match; do not revert the guard to keep it green.
    #[test]
    fn cds_to_tx_coerces_the_unknown_sentinel() {
        let tx = make_test_transcript();
        let mapper = CoordinateMapper::new(&tx);

        let pos = CdsPos::unknown(Some(5));

        let result = mapper
            .cds_to_tx(&pos)
            .expect("cds_to_tx coerces the c.? sentinel rather than declining it");
        assert_eq!(
            result,
            TxPos {
                base: 6,
                offset: Some(5),
                downstream: false,
            },
            "cds_to_tx coerces the c.? sentinel rather than declining it (today's observed behaviour, not a requirement)"
        );
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
            cds_start_incomplete: false,
            id: "NM_003742.4".to_string(),
            gene_symbol: Some("ABCB11".to_string()),
            strand: Strand::Minus,
            sequence: None,       // Not needed for coordinate tests
            cds_start: Some(128), // 127 + 1 = 128 (1-based)
            cds_end: Some(4093),  // Same (half-open end = inclusive 1-based)
            exons,
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::default(),
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
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
            special: None,
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

    /// A minimal cdot document spelled in cdot's own RAW conventions: genomic
    /// bounds 0-based half-open (`alt_start` inclusive, `alt_end` exclusive)
    /// and cDNA bounds 1-based **inclusive**.
    ///
    /// | exon | raw genomic    | raw tx (1-based incl) | expected internal tx |
    /// |------|----------------|-----------------------|----------------------|
    /// | 1    | `[5000, 5100)` | `1..100`              | `[0, 100)`           |
    /// | 2    | `[6000, 6050)` | `101..150`            | `[100, 150)`         |
    ///
    /// Two exons deliberately. The second is what makes the *half-open* half of
    /// the claim testable at all: under 0-based half-open storage exon 2's
    /// `tx_start` equals exon 1's `tx_end` exactly, where 1-based-inclusive
    /// storage would leave it one higher.
    ///
    /// `start_codon`/`stop_codon` are deliberately ABSENT. `from_genome_build`
    /// copies those through verbatim when they are present, so a fixture
    /// carrying them would bypass the conversion under test and assert nothing
    /// about the coordinate basis. Omitting them forces the genomic-CDS
    /// fallback — the arm that has to reconcile raw genomic bounds against the
    /// already-converted exon table, which is where #742's off-by-one lived.
    const RAW_CDOT_BASIS_FIXTURE: &str = r#"{
        "transcripts": {
            "NM_BASIS_TEST.1": {
                "gene_name": "BASISTEST",
                "genome_builds": {
                    "GRCh38": {
                        "contig": "NC_000001.11",
                        "strand": "+",
                        "cds_start": 5010,
                        "cds_end": 6041,
                        "exons": [
                            [5000, 5100, 1, 1, 100, "M100"],
                            [6000, 6050, 2, 101, 150, "M50"]
                        ]
                    }
                }
            }
        }
    }"#;

    /// cdot tx coordinates are stored 0-based half-open (#742).
    ///
    /// On load, `from_genome_build` converts cdot's raw 1-based-inclusive tx
    /// bounds into the engine's HGVS convention: tx 0-based half-open, the
    /// first exon starting at 0 rather than 1. Violating that convention is a
    /// guaranteed panic on the first exon, so the guard must not be able to
    /// stop running quietly.
    ///
    /// This runs on every machine, against a programmatically-built fixture.
    /// Its predecessor asserted only when an operator-created
    /// `benchmark-output/` symlink was present, and returned silently
    /// otherwise — i.e. it asserted nothing on essentially every run,
    /// including CI. The real-data cross-check lives on in
    /// `cdot_tx_basis_cross_checked_against_real_cdot` below, which is an
    /// addition to this test rather than a substitute for it.
    #[test]
    fn cdot_tx_coordinates_are_0_based_half_open() {
        use crate::data::cdot::CdotMapper;

        let cdot = CdotMapper::from_reader_with_build(RAW_CDOT_BASIS_FIXTURE.as_bytes(), "GRCh38")
            .expect("basis fixture should parse");
        let tx = cdot
            .get_transcript("NM_BASIS_TEST.1")
            .expect("basis fixture defines NM_BASIS_TEST.1");

        assert_eq!(tx.exons.len(), 2, "fixture should yield both exons");

        // 0-based: raw tx 1..100 (1-based inclusive) becomes [0, 100).
        let (first_start, first_end) = (tx.exons[0][2], tx.exons[0][3]);
        assert_eq!(
            (first_start, first_end),
            (0, 100),
            "exon 1 raw tx 1..100 (1-based inclusive) should convert to 0-based \
             half-open [0, 100), got [{first_start}, {first_end})"
        );

        // 0-based: raw tx 101..150 becomes [100, 150).
        let (second_start, second_end) = (tx.exons[1][2], tx.exons[1][3]);
        assert_eq!(
            (second_start, second_end),
            (100, 150),
            "exon 2 raw tx 101..150 (1-based inclusive) should convert to 0-based \
             half-open [100, 150), got [{second_start}, {second_end})"
        );

        // Half-open: adjoining exons meet at one shared value. Under
        // 1-based-inclusive storage exon 2 would start at exon 1's end + 1.
        assert_eq!(
            second_start, first_end,
            "half-open tx bounds should adjoin: exon 2 tx_start ({second_start}) must equal \
             exon 1 tx_end ({first_end}); a difference of 1 means the bounds are stored \
             1-based inclusive"
        );

        // CDS start is 0-based inclusive in transcript space. The fixture's
        // 5'UTR is 10 bases (raw genomic CDS start 5010 against an exon opening
        // at raw genomic 5000), so c.1 sits at tx offset 10 — it would be 11 if
        // the CDS bound were stored 1-based.
        assert_eq!(
            tx.cds_start,
            Some(10),
            "cds_start should be the 0-based tx offset of c.1 (10, after a 10-base 5'UTR)"
        );

        // `cds_end` is deliberately not asserted here. The genomic-CDS fallback
        // applies a further `+ 1` to the mapped end whose correctness depends on
        // whether cdot's raw genomic `cds_end` is inclusive or exclusive — the
        // struct doc says exclusive while the arithmetic implies inclusive. That
        // is a separate question from the coordinate basis this test pins, and
        // pinning today's value here would cement whichever reading is wrong.
    }

    /// Why the real-cdot lookup did or did not produce a path.
    ///
    /// The three non-`Found` arms are kept apart on purpose: they call for
    /// three different messages, and collapsing them into one `None` is how a
    /// skip notice comes to state a condition that does not hold. Telling an
    /// operator whose `benchmark-output/` symlink dangles or whose reference
    /// directory is unreadable to "symlink benchmark-output/ to a prepared
    /// reference" is wrong advice — they already did.
    enum CdotLookup {
        /// A `cdot-*.refseq.GRCh38.json` was found at this path.
        Found(std::path::PathBuf),
        /// The directory does not exist. The default-checkout and CI case, and
        /// also what a dangling symlink reports (both surface as `NotFound`).
        Absent,
        /// The directory exists but could not be read — permissions, or an I/O
        /// error. A wired-up-but-broken reference, not an absent one.
        Unreadable(std::io::Error),
        /// The directory was read and holds no matching cdot JSON.
        NoMatch,
    }

    /// Locate a RefSeq GRCh38 cdot JSON under `dir` without pinning a cdot
    /// release.
    ///
    /// Pattern-matching rather than naming a version is the point: the caller's
    /// predecessor hardcoded `cdot-0.2.32.refseq.GRCh38.json`, so regenerating
    /// the prepared reference at any newer cdot release (0.2.33 is current)
    /// would have returned it to never running, silently.
    ///
    /// When several releases sit side by side the pick is the
    /// **lexicographically last** filename, which is deterministic but is *not*
    /// "the newest release" — `cdot-0.2.9.…` sorts after `cdot-0.2.32.…`.
    /// Any real cdot answers the question this test asks, so the ordering only
    /// has to be stable, not semantic.
    fn locate_refseq_grch38_cdot(dir: &std::path::Path) -> CdotLookup {
        let entries = match std::fs::read_dir(dir) {
            Ok(entries) => entries,
            Err(err) if err.kind() == std::io::ErrorKind::NotFound => return CdotLookup::Absent,
            Err(err) => return CdotLookup::Unreadable(err),
        };
        let mut matches: Vec<std::path::PathBuf> = entries
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| {
                path.file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| {
                        name.starts_with("cdot-") && name.ends_with(".refseq.GRCh38.json")
                    })
            })
            .collect();
        matches.sort();
        match matches.pop() {
            Some(path) => CdotLookup::Found(path),
            None => CdotLookup::NoMatch,
        }
    }

    /// Cross-check the 0-based basis against a real prepared reference.
    ///
    /// Kept as an *addition* to `cdot_tx_coordinates_are_0_based_half_open`,
    /// never as a substitute. It is not a gate — `benchmark-output/` is an
    /// operator-created symlink to a prepared reference, absent in a default
    /// checkout and in CI — so it names the specific condition it hit and
    /// returns.
    ///
    /// **The skip notice does not reach a passing run's console.** Both
    /// `cargo test` and `cargo nextest run` capture a test's stdout/stderr and
    /// replay it only for failures, so the `eprintln!` below is visible only
    /// under `--nocapture` / `--no-capture` (or nextest's
    /// `--success-output immediate`). Do not read it as making the skip
    /// self-announcing: what stops this test from being the silent no-op its
    /// predecessor was is that the basis is now guarded unconditionally by
    /// `cdot_tx_coordinates_are_0_based_half_open`, which needs no reference at
    /// all. Turning the skip itself into a hard failure would need an opt-in
    /// switch (an env var an operator sets when they expect the reference to be
    /// wired), which is deliberately not added here.
    ///
    /// **Read what this does and does not cover.** `CdotMapper::load` prefers a
    /// sibling `.rkyv` cache over the JSON, and a prepared reference normally
    /// ships one, so on a *warm* reference this deserializes coordinates that
    /// were converted at prepare time rather than re-running
    /// `from_genome_build`. It therefore guards the prepared *artifact's*
    /// stored basis — a mis-prepared or stale reference — and NOT the
    /// conversion. Measured, not assumed: inverting the conversion to emit
    /// 1-based leaves this test green while
    /// `cdot_tx_coordinates_are_0_based_half_open` fails.
    ///
    /// That warm-reference reading is a description of the common case, not a
    /// guarantee: `load` falls back to parsing the JSON whenever the sibling
    /// `.rkyv` is missing *or older than the JSON* (its `cache_is_fresh`
    /// check), and the fallback then writes a fresh multi-hundred-MB `.rkyv`
    /// into the operator's reference directory as a side effect. So on a cold
    /// or stale reference this test does pay the ~489 MiB parse, and re-running
    /// it is cheap only from the second run on. That is the reason not to
    /// *force* the JSON arm, and it is also why the fixture above — which
    /// covers the conversion on every machine, in microseconds — is the
    /// primary guard rather than this one.
    ///
    /// The fixture's own faithfulness was checked against this same file rather
    /// than assumed: real cdot spells raw exon tx bounds 1-based **inclusive**
    /// (adjacent exons read `…4521]`, `[4522…`) with the first exon's raw
    /// `tx_start` at 1 — exactly the spelling `RAW_CDOT_BASIS_FIXTURE` uses.
    ///
    /// It asserts a structural invariant over a sample of real transcripts
    /// rather than one accession's `cds_start`. The pinned
    /// `NM_003742.4 cds_start == 127` it replaces was live external data that
    /// upstream can revise, which would have failed this test for a reason
    /// having nothing to do with the coordinate basis.
    ///
    /// Exon-to-exon tx contiguity is deliberately NOT asserted here: real cdot
    /// carries genuine transcript-coordinate holes (58 of 474,818 GRCh38
    /// multi-exon builds, 23–2718 bases), so contiguity is a property of the
    /// fixture's geometry and not of cdot at large. `tx_start == 0` on the
    /// first exon is universal, and is the discriminating one anyway — it reads
    /// 1 under 1-based storage.
    #[test]
    fn cdot_tx_basis_cross_checked_against_real_cdot() {
        use crate::data::cdot::CdotMapper;
        use std::path::PathBuf;

        let cdot_dir = PathBuf::from("benchmark-output/cdot");
        // Each arm states the condition it actually hit. A single "no cdot file
        // found" message would misdiagnose the two broken-setup cases as an
        // absent reference and hand the operator a fix they have already applied.
        const SKIP_PREAMBLE: &str = "cdot_tx_basis_cross_checked_against_real_cdot: skipping — \
             this is a cross-check only; the basis itself is guarded unconditionally by \
             cdot_tx_coordinates_are_0_based_half_open.";
        let cdot_path = match locate_refseq_grch38_cdot(&cdot_dir) {
            CdotLookup::Found(path) => path,
            CdotLookup::Absent => {
                eprintln!(
                    "{SKIP_PREAMBLE} {} does not exist (a dangling benchmark-output symlink \
                     reports the same). To run it, symlink benchmark-output/ to a prepared \
                     reference.",
                    cdot_dir.display()
                );
                return;
            }
            CdotLookup::Unreadable(err) => {
                eprintln!(
                    "{SKIP_PREAMBLE} {} exists but could not be read: {err}. The reference is \
                     wired up and BROKEN, not absent — re-linking benchmark-output/ will not \
                     fix this.",
                    cdot_dir.display()
                );
                return;
            }
            CdotLookup::NoMatch => {
                eprintln!(
                    "{SKIP_PREAMBLE} {} was read but holds no cdot-*.refseq.GRCh38.json. A \
                     prepared reference should carry one; check that prepare completed and \
                     that the RefSeq GRCh38 cdot was not stored gzipped-only.",
                    cdot_dir.display()
                );
                return;
            }
        };

        let cdot = CdotMapper::load(&cdot_path).unwrap_or_else(|err| {
            panic!("failed to load real cdot at {}: {err}", cdot_path.display())
        });

        // A sample is enough: the invariant is per-transcript, and loading the
        // full RefSeq cdot already dominates this test's cost.
        //
        // Sorted before sampling. `transcripts_on_contig` hands back
        // `contig_index` order, which on the JSON arm is `HashMap` iteration
        // order and therefore differs from process to process — an unsorted
        // `take` would compare a different 200 transcripts on every run, so a
        // failure would not reproduce and the covered set would be unknowable.
        const SAMPLE_SIZE: usize = 200;
        let mut accessions = cdot.transcripts_on_contig("NC_000001.11");
        accessions.sort_unstable();
        let attempted = accessions.len().min(SAMPLE_SIZE);
        let mut checked = 0usize;
        let mut unresolved: Vec<&str> = Vec::new();
        let mut exonless: Vec<&str> = Vec::new();
        for accession in accessions.iter().take(SAMPLE_SIZE) {
            let Some(tx) = cdot.get_transcript(accession) else {
                unresolved.push(*accession);
                continue;
            };
            let Some(first_exon) = tx.exons.first() else {
                exonless.push(*accession);
                continue;
            };
            assert_eq!(
                first_exon[2],
                0,
                "{accession}: first exon tx_start should be 0 (0-based half-open) in {}, got {}",
                cdot_path.display(),
                first_exon[2]
            );
            checked += 1;
        }

        // A zero denominator would let this report success having compared
        // nothing — the same failure mode the silent `return` had.
        assert!(
            checked > 0,
            "no transcripts checked from {} — expected transcripts on NC_000001.11, \
             so this run proved nothing",
            cdot_path.display()
        );

        // ...and a near-zero denominator would let it report success having
        // compared almost nothing, which `checked > 0` alone does not exclude:
        // 199 of 200 could fall through the two `continue` arms and the test
        // would still pass. Both arms are structurally unreachable — every
        // accession `transcripts_on_contig` yields is a key of the very
        // `transcripts` map `get_transcript` looks in first, and every cdot
        // genome build deserializes with a non-empty exon list — so a drop here
        // is information about the reference, not noise to absorb into the
        // count. Account for it rather than letting it thin the sample
        // silently.
        assert_eq!(
            checked,
            attempted,
            "compared {checked} of {attempted} sampled transcripts from {}: {} did not resolve \
             via get_transcript ({unresolved:?}) and {} carried an empty exon table \
             ({exonless:?}) — an unaccounted drop makes the comparison count meaningless",
            cdot_path.display(),
            unresolved.len(),
            exonless.len()
        );
    }
}
