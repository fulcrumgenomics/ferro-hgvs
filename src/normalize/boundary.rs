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

/// Get the shuffling boundaries for a CDS position.
///
/// Two boundary kinds are intersected:
///
///   1. **Exon bound.** When `cross_boundaries=false`, the shuffle is
///      confined to the exon containing `tx_pos`; otherwise the bound is
///      the full transcript `[1, sequence_length]`.
///
///   2. **CDS↔UTR axis bound (always applied; closes #337).** The 5'UTR
///      / CDS / 3'UTR transitions change the HGVS coordinate sub-axis
///      (`c.-N` vs `c.<N>` vs `c.*N`), so a 3'-rule shuffle must not
///      cross them — doing so would silently re-classify the variant
///      onto a different axis. The axis bound depends on which region
///      `tx_pos` lies in:
///        - 5'UTR (`tx_pos < cds_start`): `[1, cds_start - 1]`
///        - CDS   (`cds_start <= tx_pos <= cds_end`): `[cds_start, cds_end]`
///        - 3'UTR (`tx_pos > cds_end`): `[cds_end + 1, sequence_length]`
///
///      For non-coding transcripts (no `cds_start`/`cds_end`) the axis
///      bound is the full transcript.
///
/// The returned `Boundaries` is the intersection. Errors when `tx_pos`
/// falls outside every exon under `cross_boundaries=false` (i.e. is
/// intronic — exonic-shuffle code shouldn't reach here in that case).
pub fn get_cds_boundaries(
    transcript: &Transcript,
    tx_pos: u64,
    config: &NormalizeConfig,
) -> Result<Boundaries, FerroError> {
    let (exon_left, exon_right) = if config.cross_boundaries {
        (1u64, transcript.sequence_length())
    } else {
        let mut found = None;
        for exon in &transcript.exons {
            if tx_pos >= exon.start && tx_pos <= exon.end {
                found = Some((exon.start, exon.end));
                break;
            }
        }
        match found {
            Some(p) => p,
            None => {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!("Position {} is not within an exon", tx_pos),
                })
            }
        }
    };

    let seq_len = transcript.sequence_length();
    let (axis_left, axis_right) = match (transcript.cds_start, transcript.cds_end) {
        (Some(s), Some(e)) if e >= s => {
            if tx_pos < s {
                (1, s.saturating_sub(1))
            } else if tx_pos > e {
                (e + 1, seq_len)
            } else {
                (s, e)
            }
        }
        // Non-coding transcript: no axis sub-region to respect.
        _ => (1, seq_len),
    };

    Ok(Boundaries::new(
        axis_left.max(exon_left),
        axis_right.min(exon_right),
    ))
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
            sequence: Some("ATGCATGCATGC".to_string()),
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
            protein_id: None,
            exon_cigars: Vec::new(),
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
