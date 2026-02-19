//! Boundary detection for normalization
//!
//! Determines the boundaries within which a variant can be shuffled.

use crate::error::FerroError;
use crate::normalize::config::NormalizeConfig;
use crate::reference::transcript::Transcript;

/// Boundaries for variant shuffling
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Boundaries {
    /// Leftmost position (cannot shift left of this)
    pub left: u64,
    /// Rightmost position (cannot shift right of this)
    pub right: u64,
}

impl Boundaries {
    pub fn new(left: u64, right: u64) -> Self {
        Self { left, right }
    }

    /// Check if a position is within bounds (inclusive on both ends)
    ///
    /// This uses a closed interval [left, right], meaning both boundary
    /// positions are considered to be within bounds.
    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.left && pos <= self.right
    }
}

/// Get the shuffling boundaries for a genomic position
pub fn get_genomic_boundaries(
    _transcript: &Transcript,
    pos: u64,
    config: &NormalizeConfig,
) -> Result<Boundaries, FerroError> {
    // For genomic variants without transcript context,
    // use a window around the position
    let window = config.window_size;
    let left = pos.saturating_sub(window);
    let right = pos.saturating_add(window);

    Ok(Boundaries::new(left, right))
}

/// Get the shuffling boundaries for a CDS position
///
/// This considers exon boundaries and UTR/CDS boundaries
/// unless cross_boundaries is enabled.
pub fn get_cds_boundaries(
    transcript: &Transcript,
    tx_pos: u64,
    config: &NormalizeConfig,
) -> Result<Boundaries, FerroError> {
    if config.cross_boundaries {
        // If crossing is allowed, use the full transcript
        return Ok(Boundaries::new(1, transcript.sequence_length()));
    }

    // Find which exon contains this position
    for exon in &transcript.exons {
        if tx_pos >= exon.start && tx_pos <= exon.end {
            return Ok(Boundaries::new(exon.start, exon.end));
        }
    }

    // Position is intronic - this shouldn't happen for exonic variants
    Err(FerroError::InvalidCoordinates {
        msg: format!("Position {} is not within an exon", tx_pos),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCATGCATGC".to_string(),
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![Exon::new(1, 1, 4), Exon::new(2, 5, 8), Exon::new(3, 9, 12)],
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
    fn test_boundaries_contains() {
        let bounds = Boundaries::new(10, 20);
        assert!(bounds.contains(10));
        assert!(bounds.contains(15));
        assert!(bounds.contains(20));
        assert!(!bounds.contains(9));
        assert!(!bounds.contains(21));
    }

    #[test]
    fn test_cds_boundaries_within_exon() {
        let transcript = make_test_transcript();
        let config = NormalizeConfig::default();

        let bounds = get_cds_boundaries(&transcript, 2, &config).unwrap();
        assert_eq!(bounds.left, 1);
        assert_eq!(bounds.right, 4);
    }

    #[test]
    fn test_cds_boundaries_cross_allowed() {
        let transcript = make_test_transcript();
        let config = NormalizeConfig::default().allow_crossing_boundaries();

        let bounds = get_cds_boundaries(&transcript, 2, &config).unwrap();
        assert_eq!(bounds.left, 1);
        assert_eq!(bounds.right, 12);
    }
}
