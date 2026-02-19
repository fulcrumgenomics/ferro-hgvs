//! Non-coding transcript (n.) coordinate handling
//!
//! Provides validation and utilities for transcript coordinates
//! used in non-coding RNA variants (n. notation), as well as
//! intronic variant handling for all transcript types.
//!
//! # Coordinate System
//!
//! | Position Type | Basis | Notes |
//! |---------------|-------|-------|
//! | `TxPos.base` | 1-based | Transcript position (negative for upstream) |
//! | `TxPos.offset` | Signed | Intronic offset (+N for donor, -N for acceptor) |
//! | Exon boundaries | 1-based | From `Exon.start`/`Exon.end` |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::hgvs::location::{CdsPos, TxPos};
use crate::reference::transcript::{IntronPosition, SpliceSiteType, Transcript};

/// Validate a transcript position against a transcript
pub fn validate_tx_pos(pos: &TxPos, transcript: &Transcript) -> Result<(), FerroError> {
    let seq_len = transcript.sequence.len() as i64;

    // Check if position is within transcript bounds
    // Note: negative positions are allowed for upstream positions (e.g., n.-30)
    // but they're not in the actual transcript sequence
    if pos.base < 1 || pos.base > seq_len {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "Position {} is outside transcript bounds (1-{})",
                pos.base, seq_len
            ),
        });
    }

    // Check if offset is within bounds (intronic positions)
    if let Some(offset) = pos.offset {
        // For non-coding transcripts, offsets are typically not used
        // but we still validate they're reasonable
        if offset.abs() > 100000 {
            return Err(FerroError::InvalidCoordinates {
                msg: format!("Offset {} is unreasonably large", offset),
            });
        }
    }

    Ok(())
}

/// Check if a position is intronic (has an offset)
pub fn is_intronic(pos: &TxPos) -> bool {
    pos.offset.is_some() && pos.offset != Some(0)
}

/// Get the exon containing a position, if any
pub fn get_exon_for_position(
    transcript: &Transcript,
    pos: u64,
) -> Option<&crate::reference::transcript::Exon> {
    transcript.exon_at(pos)
}

/// Calculate the distance to the nearest exon boundary
pub fn distance_to_exon_boundary(transcript: &Transcript, pos: u64) -> Option<(i64, bool)> {
    if let Some(exon) = transcript.exon_at(pos) {
        let dist_to_start = pos as i64 - exon.start as i64;
        let dist_to_end = exon.end as i64 - pos as i64;

        if dist_to_start <= dist_to_end {
            Some((dist_to_start, true)) // true = closer to 5' end of exon
        } else {
            Some((dist_to_end, false)) // false = closer to 3' end of exon
        }
    } else {
        None // Position is in intron
    }
}

/// Predicted consequence of an intronic variant
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntronicConsequence {
    /// Variant at canonical splice donor site (+1, +2)
    /// High likelihood of disrupting splicing
    SpliceDonorVariant,
    /// Variant at canonical splice acceptor site (-1, -2)
    /// High likelihood of disrupting splicing
    SpliceAcceptorVariant,
    /// Variant in extended splice donor region (+3 to +6)
    /// May affect splicing efficiency
    SpliceDonorRegionVariant,
    /// Variant in polypyrimidine tract / acceptor region (-3 to -12)
    /// May affect splicing efficiency
    SpliceAcceptorRegionVariant,
    /// Variant in splice region (+7 to +20 or -13 to -20)
    /// Less likely to affect splicing
    SpliceRegionVariant,
    /// Variant near splice site (21-50bp from exon)
    /// Unlikely to affect canonical splicing
    NearSpliceSiteVariant,
    /// Deep intronic variant (>50bp from exon)
    /// May create cryptic splice sites
    IntronVariant,
}

impl IntronicConsequence {
    /// Create from an IntronPosition
    pub fn from_intron_position(pos: &IntronPosition) -> Self {
        match pos.splice_site_type() {
            SpliceSiteType::DonorCanonical => Self::SpliceDonorVariant,
            SpliceSiteType::AcceptorCanonical => Self::SpliceAcceptorVariant,
            SpliceSiteType::DonorExtended => Self::SpliceDonorRegionVariant,
            SpliceSiteType::AcceptorExtended => Self::SpliceAcceptorRegionVariant,
            SpliceSiteType::DonorRegion | SpliceSiteType::AcceptorRegion => {
                Self::SpliceRegionVariant
            }
            SpliceSiteType::NearSplice => Self::NearSpliceSiteVariant,
            SpliceSiteType::DeepIntronic => Self::IntronVariant,
        }
    }

    /// Create from a CDS position with intronic offset
    pub fn from_cds_pos(pos: &CdsPos) -> Option<Self> {
        if !pos.is_intronic() {
            return None;
        }

        let offset = pos.offset?;
        let abs_offset = offset.abs();

        Some(if offset > 0 {
            // 5' end of intron (splice donor side)
            if abs_offset <= 2 {
                Self::SpliceDonorVariant
            } else if abs_offset <= 6 {
                Self::SpliceDonorRegionVariant
            } else if abs_offset <= 20 {
                Self::SpliceRegionVariant
            } else if abs_offset <= 50 {
                Self::NearSpliceSiteVariant
            } else {
                Self::IntronVariant
            }
        } else {
            // 3' end of intron (splice acceptor side)
            if abs_offset <= 2 {
                Self::SpliceAcceptorVariant
            } else if abs_offset <= 12 {
                Self::SpliceAcceptorRegionVariant
            } else if abs_offset <= 20 {
                Self::SpliceRegionVariant
            } else if abs_offset <= 50 {
                Self::NearSpliceSiteVariant
            } else {
                Self::IntronVariant
            }
        })
    }

    /// Create from a transcript position with intronic offset
    pub fn from_tx_pos(pos: &TxPos) -> Option<Self> {
        if !pos.is_intronic() {
            return None;
        }

        let offset = pos.offset?;
        let abs_offset = offset.abs();

        Some(if offset > 0 {
            if abs_offset <= 2 {
                Self::SpliceDonorVariant
            } else if abs_offset <= 6 {
                Self::SpliceDonorRegionVariant
            } else if abs_offset <= 20 {
                Self::SpliceRegionVariant
            } else if abs_offset <= 50 {
                Self::NearSpliceSiteVariant
            } else {
                Self::IntronVariant
            }
        } else if abs_offset <= 2 {
            Self::SpliceAcceptorVariant
        } else if abs_offset <= 12 {
            Self::SpliceAcceptorRegionVariant
        } else if abs_offset <= 20 {
            Self::SpliceRegionVariant
        } else if abs_offset <= 50 {
            Self::NearSpliceSiteVariant
        } else {
            Self::IntronVariant
        })
    }

    /// Get the SO (Sequence Ontology) term for this consequence
    pub fn so_term(&self) -> &'static str {
        match self {
            Self::SpliceDonorVariant => "splice_donor_variant",
            Self::SpliceAcceptorVariant => "splice_acceptor_variant",
            Self::SpliceDonorRegionVariant => "splice_donor_region_variant",
            Self::SpliceAcceptorRegionVariant => "splice_polypyrimidine_tract_variant",
            Self::SpliceRegionVariant => "splice_region_variant",
            Self::NearSpliceSiteVariant => "intron_variant",
            Self::IntronVariant => "intron_variant",
        }
    }

    /// Get the impact level (high, moderate, low, modifier)
    pub fn impact(&self) -> &'static str {
        match self {
            Self::SpliceDonorVariant | Self::SpliceAcceptorVariant => "HIGH",
            Self::SpliceDonorRegionVariant | Self::SpliceAcceptorRegionVariant => "LOW",
            Self::SpliceRegionVariant => "LOW",
            Self::NearSpliceSiteVariant | Self::IntronVariant => "MODIFIER",
        }
    }

    /// Check if this consequence affects canonical splice sites
    pub fn affects_canonical_splice_site(&self) -> bool {
        matches!(self, Self::SpliceDonorVariant | Self::SpliceAcceptorVariant)
    }

    /// Check if this consequence may affect splicing
    pub fn may_affect_splicing(&self) -> bool {
        !matches!(self, Self::NearSpliceSiteVariant | Self::IntronVariant)
    }
}

impl std::fmt::Display for IntronicConsequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.so_term())
    }
}

/// Classify intronic position by distance from exon
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntronicRegion {
    /// At canonical splice site (1-2bp from exon)
    CanonicalSpliceSite,
    /// In extended splice region (3-20bp from exon)
    ExtendedSpliceRegion,
    /// Near splice site (21-50bp from exon)
    NearSpliceSite,
    /// Deep intronic (>50bp from exon)
    DeepIntronic,
}

impl IntronicRegion {
    /// Create from an offset value
    pub fn from_offset(offset: i64) -> Self {
        let abs_offset = offset.abs();
        if abs_offset <= 2 {
            Self::CanonicalSpliceSite
        } else if abs_offset <= 20 {
            Self::ExtendedSpliceRegion
        } else if abs_offset <= 50 {
            Self::NearSpliceSite
        } else {
            Self::DeepIntronic
        }
    }

    /// Create from a CDS position
    pub fn from_cds_pos(pos: &CdsPos) -> Option<Self> {
        pos.offset.map(Self::from_offset)
    }

    /// Create from a transcript position
    pub fn from_tx_pos(pos: &TxPos) -> Option<Self> {
        pos.offset.map(Self::from_offset)
    }
}

/// Get the intron number for a position if it's intronic
pub fn get_intron_number_for_genomic(transcript: &Transcript, genomic_pos: u64) -> Option<u32> {
    transcript
        .find_intron_at_genomic(genomic_pos)
        .map(|(_, pos)| pos.intron_number)
}

/// Check if a CDS position is in a clinically significant splice region
pub fn is_clinically_significant_splice_position(pos: &CdsPos) -> bool {
    if let Some(consequence) = IntronicConsequence::from_cds_pos(pos) {
        consequence.may_affect_splicing()
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, ManeStatus, Strand};
    use std::sync::OnceLock;

    fn make_test_transcript() -> Transcript {
        Transcript {
            id: "NR_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: "A".repeat(100),
            cds_start: None, // Non-coding
            cds_end: None,
            exons: vec![
                Exon::new(1, 1, 30),
                Exon::new(2, 50, 80),
                Exon::new(3, 90, 100),
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
    fn test_validate_tx_pos_valid() {
        let tx = make_test_transcript();

        // Valid positions
        assert!(validate_tx_pos(&TxPos::new(1), &tx).is_ok());
        assert!(validate_tx_pos(&TxPos::new(50), &tx).is_ok());
        assert!(validate_tx_pos(&TxPos::new(100), &tx).is_ok());
    }

    #[test]
    fn test_validate_tx_pos_invalid() {
        let tx = make_test_transcript();

        // Invalid positions
        assert!(validate_tx_pos(&TxPos::new(0), &tx).is_err());
        assert!(validate_tx_pos(&TxPos::new(101), &tx).is_err());
    }

    #[test]
    fn test_is_intronic() {
        let pos_exonic = TxPos::new(50);
        let pos_intronic = TxPos::with_offset(50, 5);
        let pos_zero_offset = TxPos::with_offset(50, 0);

        assert!(!is_intronic(&pos_exonic));
        assert!(is_intronic(&pos_intronic));
        assert!(!is_intronic(&pos_zero_offset));
    }

    #[test]
    fn test_get_exon_for_position() {
        let tx = make_test_transcript();

        // Position in exon 1
        let exon = get_exon_for_position(&tx, 15);
        assert!(exon.is_some());
        assert_eq!(exon.unwrap().number, 1);

        // Position in intron
        let exon = get_exon_for_position(&tx, 40);
        assert!(exon.is_none());

        // Position in exon 2
        let exon = get_exon_for_position(&tx, 65);
        assert!(exon.is_some());
        assert_eq!(exon.unwrap().number, 2);
    }

    #[test]
    fn test_distance_to_exon_boundary() {
        let tx = make_test_transcript();

        // Position near start of exon 1
        let dist = distance_to_exon_boundary(&tx, 5);
        assert!(dist.is_some());
        let (d, near_start) = dist.unwrap();
        assert_eq!(d, 4); // 5 - 1 = 4
        assert!(near_start);

        // Position near end of exon 1
        let dist = distance_to_exon_boundary(&tx, 28);
        assert!(dist.is_some());
        let (d, near_start) = dist.unwrap();
        assert_eq!(d, 2); // 30 - 28 = 2
        assert!(!near_start);

        // Position in intron
        let dist = distance_to_exon_boundary(&tx, 40);
        assert!(dist.is_none());
    }

    #[test]
    fn test_intronic_consequence_from_cds_pos() {
        // Splice donor +1
        let pos = CdsPos::with_offset(100, 1);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceDonorVariant));

        // Splice donor +5
        let pos = CdsPos::with_offset(100, 5);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceDonorRegionVariant)
        );

        // Splice acceptor -2
        let pos = CdsPos::with_offset(100, -2);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceAcceptorVariant)
        );

        // Splice acceptor region -10
        let pos = CdsPos::with_offset(100, -10);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceAcceptorRegionVariant)
        );

        // Splice region +15
        let pos = CdsPos::with_offset(100, 15);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceRegionVariant));

        // Near splice site +30
        let pos = CdsPos::with_offset(100, 30);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::NearSpliceSiteVariant)
        );

        // Deep intronic +100
        let pos = CdsPos::with_offset(100, 100);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::IntronVariant));

        // Non-intronic position
        let pos = CdsPos::new(100);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, None);
    }

    #[test]
    fn test_intronic_consequence_so_terms() {
        assert_eq!(
            IntronicConsequence::SpliceDonorVariant.so_term(),
            "splice_donor_variant"
        );
        assert_eq!(
            IntronicConsequence::SpliceAcceptorVariant.so_term(),
            "splice_acceptor_variant"
        );
        assert_eq!(
            IntronicConsequence::IntronVariant.so_term(),
            "intron_variant"
        );
    }

    #[test]
    fn test_intronic_consequence_impact() {
        // High impact for canonical splice sites
        assert_eq!(IntronicConsequence::SpliceDonorVariant.impact(), "HIGH");
        assert_eq!(IntronicConsequence::SpliceAcceptorVariant.impact(), "HIGH");

        // Low impact for extended regions
        assert_eq!(
            IntronicConsequence::SpliceDonorRegionVariant.impact(),
            "LOW"
        );
        assert_eq!(IntronicConsequence::SpliceRegionVariant.impact(), "LOW");

        // Modifier for deep intronic
        assert_eq!(IntronicConsequence::IntronVariant.impact(), "MODIFIER");
    }

    #[test]
    fn test_intronic_consequence_affects_splice() {
        assert!(IntronicConsequence::SpliceDonorVariant.affects_canonical_splice_site());
        assert!(IntronicConsequence::SpliceAcceptorVariant.affects_canonical_splice_site());
        assert!(!IntronicConsequence::IntronVariant.affects_canonical_splice_site());

        assert!(IntronicConsequence::SpliceDonorVariant.may_affect_splicing());
        assert!(IntronicConsequence::SpliceRegionVariant.may_affect_splicing());
        assert!(!IntronicConsequence::IntronVariant.may_affect_splicing());
    }

    #[test]
    fn test_intronic_region_from_offset() {
        assert_eq!(
            IntronicRegion::from_offset(1),
            IntronicRegion::CanonicalSpliceSite
        );
        assert_eq!(
            IntronicRegion::from_offset(-2),
            IntronicRegion::CanonicalSpliceSite
        );
        assert_eq!(
            IntronicRegion::from_offset(10),
            IntronicRegion::ExtendedSpliceRegion
        );
        assert_eq!(
            IntronicRegion::from_offset(-15),
            IntronicRegion::ExtendedSpliceRegion
        );
        assert_eq!(
            IntronicRegion::from_offset(35),
            IntronicRegion::NearSpliceSite
        );
        assert_eq!(
            IntronicRegion::from_offset(100),
            IntronicRegion::DeepIntronic
        );
    }

    #[test]
    fn test_is_clinically_significant_splice() {
        // Canonical splice sites are significant
        assert!(is_clinically_significant_splice_position(
            &CdsPos::with_offset(100, 1)
        ));
        assert!(is_clinically_significant_splice_position(
            &CdsPos::with_offset(100, -2)
        ));

        // Extended regions are significant
        assert!(is_clinically_significant_splice_position(
            &CdsPos::with_offset(100, 5)
        ));

        // Deep intronic is not clinically significant for splicing
        assert!(!is_clinically_significant_splice_position(
            &CdsPos::with_offset(100, 100)
        ));

        // Non-intronic is not significant
        assert!(!is_clinically_significant_splice_position(&CdsPos::new(
            100
        )));
    }

    #[test]
    fn test_validate_tx_pos_with_large_offset() {
        let tx = make_test_transcript();

        // Large offset should fail
        let pos = TxPos::with_offset(50, 100001);
        assert!(validate_tx_pos(&pos, &tx).is_err());

        // Reasonable offset should pass
        let pos = TxPos::with_offset(50, 1000);
        assert!(validate_tx_pos(&pos, &tx).is_ok());

        // Negative offset should work too
        let pos = TxPos::with_offset(50, -1000);
        assert!(validate_tx_pos(&pos, &tx).is_ok());
    }

    #[test]
    fn test_intronic_consequence_from_tx_pos() {
        // Splice donor +1
        let pos = TxPos::with_offset(50, 1);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceDonorVariant));

        // Splice donor +2
        let pos = TxPos::with_offset(50, 2);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceDonorVariant));

        // Splice donor region +5
        let pos = TxPos::with_offset(50, 5);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceDonorRegionVariant)
        );

        // Splice region +15
        let pos = TxPos::with_offset(50, 15);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceRegionVariant));

        // Near splice site +35
        let pos = TxPos::with_offset(50, 35);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::NearSpliceSiteVariant)
        );

        // Deep intronic +100
        let pos = TxPos::with_offset(50, 100);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::IntronVariant));

        // Splice acceptor -1
        let pos = TxPos::with_offset(50, -1);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceAcceptorVariant)
        );

        // Splice acceptor region -10
        let pos = TxPos::with_offset(50, -10);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceAcceptorRegionVariant)
        );

        // Splice region -15
        let pos = TxPos::with_offset(50, -15);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceRegionVariant));

        // Near splice site -35
        let pos = TxPos::with_offset(50, -35);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::NearSpliceSiteVariant)
        );

        // Non-intronic
        let pos = TxPos::new(50);
        let consequence = IntronicConsequence::from_tx_pos(&pos);
        assert!(consequence.is_none());
    }

    #[test]
    fn test_intronic_consequence_from_intron_position() {
        use crate::reference::transcript::{IntronBoundary, IntronPosition};

        // Splice donor (canonical)
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 1,
            tx_boundary_pos: 100,
            intron_length: 1000,
        };
        let consequence = IntronicConsequence::from_intron_position(&pos);
        assert_eq!(consequence, IntronicConsequence::SpliceDonorVariant);

        // Splice acceptor (canonical)
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::ThreePrime,
            offset: -2,
            tx_boundary_pos: 200,
            intron_length: 1000,
        };
        let consequence = IntronicConsequence::from_intron_position(&pos);
        assert_eq!(consequence, IntronicConsequence::SpliceAcceptorVariant);

        // Deep intronic
        let pos = IntronPosition {
            intron_number: 1,
            boundary: IntronBoundary::FivePrime,
            offset: 100,
            tx_boundary_pos: 100,
            intron_length: 1000,
        };
        let consequence = IntronicConsequence::from_intron_position(&pos);
        assert_eq!(consequence, IntronicConsequence::IntronVariant);
    }

    #[test]
    fn test_intronic_consequence_display() {
        assert_eq!(
            format!("{}", IntronicConsequence::SpliceDonorVariant),
            "splice_donor_variant"
        );
        assert_eq!(
            format!("{}", IntronicConsequence::SpliceAcceptorVariant),
            "splice_acceptor_variant"
        );
        assert_eq!(
            format!("{}", IntronicConsequence::SpliceDonorRegionVariant),
            "splice_donor_region_variant"
        );
        assert_eq!(
            format!("{}", IntronicConsequence::SpliceAcceptorRegionVariant),
            "splice_polypyrimidine_tract_variant"
        );
        assert_eq!(
            format!("{}", IntronicConsequence::SpliceRegionVariant),
            "splice_region_variant"
        );
        assert_eq!(
            format!("{}", IntronicConsequence::NearSpliceSiteVariant),
            "intron_variant"
        );
        assert_eq!(
            format!("{}", IntronicConsequence::IntronVariant),
            "intron_variant"
        );
    }

    #[test]
    fn test_intronic_consequence_may_affect_splicing() {
        // These may affect splicing
        assert!(IntronicConsequence::SpliceDonorVariant.may_affect_splicing());
        assert!(IntronicConsequence::SpliceAcceptorVariant.may_affect_splicing());
        assert!(IntronicConsequence::SpliceDonorRegionVariant.may_affect_splicing());
        assert!(IntronicConsequence::SpliceAcceptorRegionVariant.may_affect_splicing());
        assert!(IntronicConsequence::SpliceRegionVariant.may_affect_splicing());

        // These do not affect splicing
        assert!(!IntronicConsequence::NearSpliceSiteVariant.may_affect_splicing());
        assert!(!IntronicConsequence::IntronVariant.may_affect_splicing());
    }

    #[test]
    fn test_intronic_region_from_cds_pos() {
        let pos = CdsPos::with_offset(100, 1);
        assert_eq!(
            IntronicRegion::from_cds_pos(&pos),
            Some(IntronicRegion::CanonicalSpliceSite)
        );

        let pos = CdsPos::with_offset(100, 10);
        assert_eq!(
            IntronicRegion::from_cds_pos(&pos),
            Some(IntronicRegion::ExtendedSpliceRegion)
        );

        let pos = CdsPos::with_offset(100, 30);
        assert_eq!(
            IntronicRegion::from_cds_pos(&pos),
            Some(IntronicRegion::NearSpliceSite)
        );

        let pos = CdsPos::with_offset(100, 100);
        assert_eq!(
            IntronicRegion::from_cds_pos(&pos),
            Some(IntronicRegion::DeepIntronic)
        );

        // Non-intronic
        let pos = CdsPos::new(100);
        assert_eq!(IntronicRegion::from_cds_pos(&pos), None);
    }

    #[test]
    fn test_intronic_region_from_tx_pos() {
        let pos = TxPos::with_offset(50, 2);
        assert_eq!(
            IntronicRegion::from_tx_pos(&pos),
            Some(IntronicRegion::CanonicalSpliceSite)
        );

        let pos = TxPos::with_offset(50, -15);
        assert_eq!(
            IntronicRegion::from_tx_pos(&pos),
            Some(IntronicRegion::ExtendedSpliceRegion)
        );

        let pos = TxPos::with_offset(50, 40);
        assert_eq!(
            IntronicRegion::from_tx_pos(&pos),
            Some(IntronicRegion::NearSpliceSite)
        );

        let pos = TxPos::with_offset(50, -60);
        assert_eq!(
            IntronicRegion::from_tx_pos(&pos),
            Some(IntronicRegion::DeepIntronic)
        );

        // Non-intronic
        let pos = TxPos::new(50);
        assert_eq!(IntronicRegion::from_tx_pos(&pos), None);
    }

    #[test]
    fn test_get_intron_number_for_genomic() {
        // Create transcript with genomic coords
        let tx = Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: "A".repeat(200),
            cds_start: Some(1),
            cds_end: Some(200),
            exons: vec![
                Exon::with_genomic(1, 1, 50, 1000, 1049),
                Exon::with_genomic(2, 51, 100, 2000, 2049),
                Exon::with_genomic(3, 101, 200, 3000, 3099),
            ],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(3099),
            genome_build: Default::default(),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        // Position in intron 1 (between exon 1 and 2)
        let intron_num = get_intron_number_for_genomic(&tx, 1500);
        assert_eq!(intron_num, Some(1));

        // Position in intron 2 (between exon 2 and 3)
        let intron_num = get_intron_number_for_genomic(&tx, 2500);
        assert_eq!(intron_num, Some(2));

        // Position in exon should return None
        let intron_num = get_intron_number_for_genomic(&tx, 1025);
        assert_eq!(intron_num, None);
    }

    #[test]
    fn test_intronic_consequence_acceptor_extended_boundary() {
        // Test boundary at -12 (acceptor extended)
        let pos = CdsPos::with_offset(100, -12);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceAcceptorRegionVariant)
        );

        // Test -13 (splice region)
        let pos = CdsPos::with_offset(100, -13);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceRegionVariant));
    }

    #[test]
    fn test_intronic_consequence_donor_extended_boundary() {
        // Test +6 (donor extended)
        let pos = CdsPos::with_offset(100, 6);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::SpliceDonorRegionVariant)
        );

        // Test +7 (splice region)
        let pos = CdsPos::with_offset(100, 7);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::SpliceRegionVariant));
    }

    #[test]
    fn test_intronic_consequence_near_boundary() {
        // Test +50 (near splice)
        let pos = CdsPos::with_offset(100, 50);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(
            consequence,
            Some(IntronicConsequence::NearSpliceSiteVariant)
        );

        // Test +51 (deep intronic)
        let pos = CdsPos::with_offset(100, 51);
        let consequence = IntronicConsequence::from_cds_pos(&pos);
        assert_eq!(consequence, Some(IntronicConsequence::IntronVariant));
    }
}
