//! Protein (p.) coordinate handling
//!
//! # Coordinate System
//!
//! | Position Type | Basis | Notes |
//! |---------------|-------|-------|
//! | `ProtPos.number` | 1-based | Amino acid position |
//! | CDS range return | 1-based | `(start, end)` inclusive |
//! | Codon frame | 1-based | Position within codon (1, 2, or 3) |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::hgvs::location::ProtPos;
use crate::reference::transcript::Transcript;

/// Validate a protein position against a transcript
pub fn validate_prot_pos(pos: &ProtPos, transcript: &Transcript) -> Result<(), FerroError> {
    let cds_length = transcript
        .cds_length()
        .ok_or_else(|| FerroError::ConversionError {
            msg: "Transcript has no CDS".to_string(),
        })?;

    // Protein length is CDS length / 3 (including stop codon)
    let protein_length = cds_length / 3;

    if pos.number < 1 {
        return Err(FerroError::InvalidCoordinates {
            msg: "Protein position must be >= 1".to_string(),
        });
    }

    if pos.number > protein_length {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "Protein position {} exceeds protein length {}",
                pos.number, protein_length
            ),
        });
    }

    Ok(())
}

/// Convert protein position to CDS position range
///
/// Returns the CDS positions for the codon encoding this amino acid
pub fn protein_to_cds_range(pos: &ProtPos) -> (i64, i64) {
    let start = ((pos.number - 1) * 3 + 1) as i64;
    let end = (pos.number * 3) as i64;
    (start, end)
}

/// Get the codon frame (1, 2, or 3) for a CDS position
pub fn get_codon_frame(cds_pos: i64) -> u8 {
    match (cds_pos - 1) % 3 {
        0 => 1,
        1 => 2,
        _ => 3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::location::AminoAcid;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: "ATGCCCAAAGGGTTTTAG".to_string(), // 18 bases = 6 codons
            cds_start: Some(1),
            cds_end: Some(18),
            exons: vec![Exon::new(1, 1, 18)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_validate_prot_pos() {
        let tx = make_test_transcript();

        // Valid positions
        assert!(validate_prot_pos(&ProtPos::new(AminoAcid::Met, 1), &tx).is_ok());
        assert!(validate_prot_pos(&ProtPos::new(AminoAcid::Ter, 6), &tx).is_ok());

        // Invalid position
        assert!(validate_prot_pos(&ProtPos::new(AminoAcid::Xaa, 10), &tx).is_err());
    }

    #[test]
    fn test_protein_to_cds_range() {
        // First codon: c.1_3
        let (start, end) = protein_to_cds_range(&ProtPos::new(AminoAcid::Met, 1));
        assert_eq!(start, 1);
        assert_eq!(end, 3);

        // Second codon: c.4_6
        let (start, end) = protein_to_cds_range(&ProtPos::new(AminoAcid::Pro, 2));
        assert_eq!(start, 4);
        assert_eq!(end, 6);
    }

    #[test]
    fn test_get_codon_frame() {
        assert_eq!(get_codon_frame(1), 1);
        assert_eq!(get_codon_frame(2), 2);
        assert_eq!(get_codon_frame(3), 3);
        assert_eq!(get_codon_frame(4), 1);
        assert_eq!(get_codon_frame(5), 2);
        assert_eq!(get_codon_frame(6), 3);
    }
}
