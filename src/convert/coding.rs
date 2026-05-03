//! CDS (c.) coordinate handling
//!
//! # Coordinate System
//!
//! | Position Type | Basis | Notes |
//! |---------------|-------|-------|
//! | `CdsPos.base` | 1-based | Positive for CDS, negative for 5'UTR |
//! | `cds_start`, `cds_end` | 1-based | From transcript metadata |
//! | Return values | 0-based | For array indexing |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::hgvs::location::CdsPos;
use crate::reference::transcript::Transcript;

/// Validate a CDS position against a transcript
pub fn validate_cds_pos(pos: &CdsPos, transcript: &Transcript) -> Result<(), FerroError> {
    let cds_start = transcript
        .cds_start
        .ok_or_else(|| FerroError::ConversionError {
            msg: "Transcript has no CDS".to_string(),
        })?;
    let cds_end = transcript
        .cds_end
        .ok_or_else(|| FerroError::ConversionError {
            msg: "Transcript has no CDS end".to_string(),
        })?;

    let cds_length = (cds_end - cds_start + 1) as i64;

    if pos.utr3 {
        // 3' UTR position
        let utr3_length = transcript.sequence_length() - cds_end;
        if pos.base < 1 || pos.base > utr3_length as i64 {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "3' UTR position *{} is out of range (max *{})",
                    pos.base, utr3_length
                ),
            });
        }
    } else if pos.base < 0 {
        // 5' UTR position
        let utr5_length = (cds_start - 1) as i64;
        if pos.base.abs() > utr5_length {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "5' UTR position {} is out of range (min -{})",
                    pos.base, utr5_length
                ),
            });
        }
    } else if pos.base > cds_length {
        // CDS position out of range
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "CDS position {} is out of range (max {})",
                pos.base, cds_length
            ),
        });
    }

    // Validate intronic offset if present
    if let Some(offset) = pos.offset {
        if offset.abs() > 1_000_000 {
            // Sanity check for very large offsets
            return Err(FerroError::InvalidCoordinates {
                msg: format!("Intronic offset {} is unreasonably large", offset),
            });
        }
    }

    Ok(())
}

/// Convert a CDS position to a 0-based transcript position
pub fn cds_to_transcript_pos(pos: &CdsPos, transcript: &Transcript) -> Result<u64, FerroError> {
    validate_cds_pos(pos, transcript)?;

    let cds_start = transcript.cds_start.unwrap(); // Already validated
    let cds_end = transcript.cds_end.unwrap();

    let tx_pos = if pos.utr3 {
        cds_end + pos.base as u64
    } else if pos.base < 0 {
        // 5'UTR: HGVS numbering skips c.0 (c.-1 is the base immediately
        // upstream of c.1), so c.-N maps to tx position cds_start - N.
        // Issue #97.
        (cds_start as i64 + pos.base) as u64
    } else if pos.base == 0 {
        // c.0 is not a valid HGVS position; preserve the legacy mapping
        // (last 5'UTR base, equivalent to c.-1) rather than failing.
        cds_start.saturating_sub(1)
    } else {
        cds_start + pos.base as u64 - 1
    };

    Ok(tx_pos)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: "AAAAATGCCCAAAGGGTTTTAAAAAA".to_string(), // 26 bases
            cds_start: Some(6),
            cds_end: Some(20),
            exons: vec![Exon::new(1, 1, 26)],
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
    fn test_validate_cds_pos_normal() {
        let tx = make_test_transcript();
        assert!(validate_cds_pos(&CdsPos::new(1), &tx).is_ok());
        assert!(validate_cds_pos(&CdsPos::new(15), &tx).is_ok());
    }

    #[test]
    fn test_validate_cds_pos_out_of_range() {
        let tx = make_test_transcript();
        // CDS is 15bp (positions 1-15)
        assert!(validate_cds_pos(&CdsPos::new(20), &tx).is_err());
    }

    #[test]
    fn test_validate_cds_pos_5utr() {
        let tx = make_test_transcript();
        // 5' UTR is 5bp (positions -1 to -5)
        assert!(validate_cds_pos(&CdsPos::new(-3), &tx).is_ok());
        assert!(validate_cds_pos(&CdsPos::new(-10), &tx).is_err());
    }

    #[test]
    fn test_validate_cds_pos_3utr() {
        let tx = make_test_transcript();
        // 3' UTR is 6bp (positions *1 to *6)
        assert!(validate_cds_pos(&CdsPos::utr3(3), &tx).is_ok());
        assert!(validate_cds_pos(&CdsPos::utr3(10), &tx).is_err());
    }

    #[test]
    fn test_cds_to_transcript_pos_5utr_and_legacy_c0() {
        // Pin the issue #97 fix and the legacy c.0 mapping in the public
        // helper, mirroring Normalizer::cds_to_tx_pos.
        //
        // 5'UTR (issue #97): HGVS skips c.0, so c.-N maps to tx position
        // cds_start - N — i.e. (cds_start as i64 + pos.base) as u64.
        let tx = make_test_transcript();
        let cds_start = tx.cds_start.unwrap(); // 6

        let result = cds_to_transcript_pos(&CdsPos::new(-1), &tx).unwrap();
        assert_eq!(result, (cds_start as i64 + (-1)) as u64);

        // Legacy c.0 mapping: treat as the last 5'UTR base (== c.-1) via
        // cds_start.saturating_sub(1) rather than failing.
        let result = cds_to_transcript_pos(&CdsPos::new(0), &tx).unwrap();
        assert_eq!(result, cds_start.saturating_sub(1));

        // Defensive: cds_start == 0 must not underflow. saturating_sub
        // pins the result at 0 instead of wrapping to u64::MAX.
        let tx_zero_cds = Transcript {
            id: "NM_TEST_ZERO.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: "ATGCCCAAAGGG".to_string(),
            cds_start: Some(0),
            cds_end: Some(5),
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        let result = cds_to_transcript_pos(&CdsPos::new(0), &tx_zero_cds).unwrap();
        assert_eq!(result, tx_zero_cds.cds_start.unwrap().saturating_sub(1));
        assert_eq!(result, 0);
    }
}
