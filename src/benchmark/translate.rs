//! Shared CDS translation utilities for benchmark tools.
//!
//! This module provides functions for translating coding DNA sequences (CDS)
//! to protein sequences using the standard genetic code. Used by both
//! mutalyzer cache preparation and biocommons SeqRepo preparation.

use crate::backtranslate::codon::{Codon, CodonTable};

/// Translate a CDS sequence to protein sequence.
///
/// Uses the standard genetic code. Stops at the first stop codon.
/// Returns None if the sequence contains invalid nucleotides.
///
/// # Arguments
/// * `cds_sequence` - The coding DNA sequence (must be complete codons)
///
/// # Returns
/// The translated amino acid sequence as a string (single-letter codes),
/// or None if the sequence contains invalid nucleotides.
pub fn translate_cds_to_protein(cds_sequence: &str) -> Option<String> {
    let table = CodonTable::standard();
    let mut protein = String::new();

    // Process complete codons only (ignore trailing 1-2 nucleotides)
    let codon_count = cds_sequence.len() / 3;

    for i in 0..codon_count {
        let start = i * 3;
        let codon_str = &cds_sequence[start..start + 3];
        if let Some(codon) = Codon::parse(codon_str) {
            if table.is_stop(&codon) {
                break; // Stop at stop codon
            }
            if let Some(aa) = table.amino_acid_for(&codon) {
                use crate::hgvs::location::AminoAcid;
                let _: &AminoAcid = aa; // type assertion
                protein.push(aa.to_one_letter());
            } else {
                return None; // Unknown codon
            }
        } else {
            return None; // Invalid nucleotide (e.g., N)
        }
    }

    Some(protein)
}

/// Derive protein sequence from transcript sequence and CDS coordinates.
///
/// # Arguments
/// * `transcript_seq` - Full transcript sequence
/// * `cds_start` - 0-based CDS start position
/// * `cds_end` - 0-based CDS end position (exclusive)
///
/// # Returns
/// Translated protein sequence, or None if coordinates invalid or translation fails
pub fn derive_protein_from_transcript(
    transcript_seq: &str,
    cds_start: usize,
    cds_end: usize,
) -> Option<String> {
    if cds_end > transcript_seq.len() || cds_start >= cds_end {
        return None;
    }

    let cds = &transcript_seq[cds_start..cds_end];
    translate_cds_to_protein(cds)
}

/// Convert 1-based inclusive CDS coordinates to 0-based half-open.
///
/// GenBank/HGVS uses 1-based inclusive coordinates; Rust uses 0-based half-open.
///
/// # Arguments
/// * `cds_start` - 1-based inclusive start position
/// * `cds_end` - 1-based inclusive end position
///
/// # Returns
/// Tuple of (0-based start, 0-based exclusive end)
pub fn cds_coords_to_indices(cds_start: u64, cds_end: u64) -> (usize, usize) {
    ((cds_start - 1) as usize, cds_end as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate_simple() {
        // ATG = Met, TGG = Trp, TAA = Stop
        assert_eq!(
            translate_cds_to_protein("ATGTGGTAA"),
            Some("MW".to_string())
        );
    }

    #[test]
    fn test_translate_no_stop() {
        // No stop codon - translates all codons
        assert_eq!(translate_cds_to_protein("ATGTGG"), Some("MW".to_string()));
    }

    #[test]
    fn test_translate_invalid() {
        // N is ambiguous, should fail
        assert_eq!(translate_cds_to_protein("ATGNTG"), None);
    }

    #[test]
    fn test_derive_protein() {
        // 5'UTR (3bp) + CDS (9bp) + 3'UTR (3bp)
        let transcript = "AAAATGTGGTAAGGG";
        assert_eq!(
            derive_protein_from_transcript(transcript, 3, 12),
            Some("MW".to_string())
        );
    }

    #[test]
    fn test_derive_protein_invalid_coords() {
        let transcript = "ATGTGGTAA";
        // Start >= end
        assert_eq!(derive_protein_from_transcript(transcript, 5, 3), None);
        // End > length
        assert_eq!(derive_protein_from_transcript(transcript, 0, 100), None);
    }

    #[test]
    fn test_cds_coords_conversion() {
        // 1-based [1, 10] -> 0-based [0, 10)
        assert_eq!(cds_coords_to_indices(1, 10), (0, 10));
        // 1-based [100, 500] -> 0-based [99, 500)
        assert_eq!(cds_coords_to_indices(100, 500), (99, 500));
    }
}
