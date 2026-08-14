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
    // `sequence_length` is exon-derived, so this works for both
    // FASTA-backed and coordinate-only transcripts.
    let seq_len = transcript.sequence_length() as i64;

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

    // An unknown offset (`n.100+?` / `n.100-?`) is not a distance, so it is
    // not a coordinate this can validate — say so rather than reporting it as
    // an out-of-range number (#1767).
    if pos.has_unknown_offset() {
        return Err(FerroError::InvalidCoordinates {
            msg: "Unknown intronic offset (`+?` / `-?`) denotes no transcript coordinate"
                .to_string(),
        });
    }

    // Check if offset is within bounds (intronic positions)
    if let Some(offset) = pos.offset {
        // For non-coding transcripts, offsets are typically not used
        // but we still validate they're reasonable. `unsigned_abs` because
        // `i64::MIN.abs()` overflows.
        if offset.unsigned_abs() > 100_000 {
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
    /// Create from a [`SpliceSiteType`] bin.
    ///
    /// The distance thresholds are **not** restated here — they live in one
    /// place, [`SpliceSiteType::from_distance_on_side`]. This is only the
    /// projection of that bin onto this enum, which is coarser: the two
    /// `*Region` bins collapse into [`Self::SpliceRegionVariant`].
    pub fn from_splice_site_type(site: SpliceSiteType) -> Self {
        match site {
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

    /// Create from an IntronPosition
    pub fn from_intron_position(pos: &IntronPosition) -> Self {
        Self::from_splice_site_type(pos.splice_site_type())
    }

    /// Create from a CDS position with intronic offset.
    ///
    /// Returns `None` for a non-intronic position, and for an unknown offset
    /// (`c.100+?` / `c.100-?`): every variant of this enum states a distance
    /// from the exon boundary, and the sentinels denote a position unbounded in
    /// one direction, so there is no distance to state (#1767). `is_intronic`
    /// is satisfied by a sentinel, so the `has_unknown_offset` guard is what
    /// keeps it out — not the `is_intronic` check beside it.
    ///
    /// Positive offsets are the donor side, zero and negative the acceptor
    /// side — see [`SpliceSiteType::from_signed_offset`], which selects the side
    /// and hands off to [`SpliceSiteType::from_distance_on_side`], where the
    /// thresholds live. (`is_intronic` already excludes `Some(0)`, so the sign
    /// rule's zero case is unreachable from here.)
    pub fn from_cds_pos(pos: &CdsPos) -> Option<Self> {
        if !pos.is_intronic() || pos.has_unknown_offset() {
            return None;
        }
        Some(Self::from_splice_site_type(
            SpliceSiteType::from_signed_offset(pos.offset?),
        ))
    }

    /// Create from a transcript position with intronic offset.
    ///
    /// Returns `None` for a non-intronic position, and declines an unknown
    /// offset (`n.100+?` / `n.100-?`) for the same reason as
    /// [`IntronicConsequence::from_cds_pos`] (#1767). Same ladder and same sign
    /// convention as that method.
    pub fn from_tx_pos(pos: &TxPos) -> Option<Self> {
        if !pos.is_intronic() || pos.has_unknown_offset() {
            return None;
        }
        Some(Self::from_splice_site_type(
            SpliceSiteType::from_signed_offset(pos.offset?),
        ))
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
    /// Classify an offset into its distance band, or decline.
    ///
    /// Returns `None` for the unknown-offset sentinels (`+?` / `-?`, carried as
    /// [`OFFSET_UNKNOWN_POSITIVE`] / [`OFFSET_UNKNOWN_NEGATIVE`]). Every variant
    /// of this enum names a **distance band**, and a sentinel states that the
    /// distance is unknown, so none of the four is a true answer for it.
    ///
    /// #1767 fixed the arithmetic here — the magnitude is taken with
    /// `unsigned_abs`, so no `i64` can overflow it — but left the signature
    /// promising an answer this function cannot always give: a sentinel landed
    /// in `DeepIntronic`, a distance claim its input does not support. #1841
    /// took that break. The precondition used to be documented ("screen with
    /// `has_unknown_offset` before calling"); it is now in the type.
    ///
    /// The bands are the side-agnostic collapse of the same ladder
    /// [`SpliceSiteType::from_distance_on_side`] owns, reached through
    /// [`SpliceSiteType::from_signed_offset`]: the donor/acceptor distinction is
    /// dropped and the extended + region bins merge into
    /// [`Self::ExtendedSpliceRegion`]. It used to restate a partial copy of the
    /// thresholds (2/20/50, with the 6- and 12-rungs missing), which is why it
    /// derives rather than repeats them — the rungs it drops must be dropped by
    /// the mapping, not by a second and independently-editable ladder.
    ///
    /// [`OFFSET_UNKNOWN_POSITIVE`]: crate::hgvs::parser::position::OFFSET_UNKNOWN_POSITIVE
    /// [`OFFSET_UNKNOWN_NEGATIVE`]: crate::hgvs::parser::position::OFFSET_UNKNOWN_NEGATIVE
    pub fn from_offset(offset: i64) -> Option<Self> {
        if crate::hgvs::parser::position::is_unknown_offset(offset) {
            return None;
        }
        Some(match SpliceSiteType::from_signed_offset(offset) {
            SpliceSiteType::DonorCanonical | SpliceSiteType::AcceptorCanonical => {
                Self::CanonicalSpliceSite
            }
            SpliceSiteType::DonorExtended
            | SpliceSiteType::AcceptorExtended
            | SpliceSiteType::DonorRegion
            | SpliceSiteType::AcceptorRegion => Self::ExtendedSpliceRegion,
            SpliceSiteType::NearSplice => Self::NearSpliceSite,
            SpliceSiteType::DeepIntronic => Self::DeepIntronic,
        })
    }

    /// Create from a CDS position.
    ///
    /// `None` when the position is not intronic, and when its offset is unknown.
    /// The unknown-offset rule is [`IntronicRegion::from_offset`]'s own since
    /// #1841, so this no longer screens for it separately — one rule, one site,
    /// which is the discipline #1767 applied to `is_unknown_offset` itself.
    ///
    /// The sibling ladder [`IntronicConsequence`] still screens in *its*
    /// position constructors rather than in a scalar entry point, and that
    /// asymmetry is deliberate rather than a missed cleanup: it exposes no
    /// public `from_offset`, so there is no scalar site to move the rule into.
    pub fn from_cds_pos(pos: &CdsPos) -> Option<Self> {
        pos.offset.and_then(Self::from_offset)
    }

    /// Create from a transcript position.
    ///
    /// Declines an unknown offset, as [`IntronicRegion::from_cds_pos`] does.
    pub fn from_tx_pos(pos: &TxPos) -> Option<Self> {
        pos.offset.and_then(Self::from_offset)
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
            cds_start_incomplete: false,
            id: "NR_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(100)),
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
            protein_id: None,
            exon_cigars: Vec::new(),
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

    /// The #1841 break: the sentinels are out of the enum's domain, and the
    /// signature now says so. `i64::MIN + 1` and `i64::MAX - 1` are the
    /// discriminating neighbours — they are *not* sentinels, so they must still
    /// classify (as `DeepIntronic`), which is what keeps the decline keyed on
    /// the sentinel rather than on "a very large magnitude".
    #[test]
    fn test_intronic_region_from_offset_declines_an_unknown_offset() {
        assert_eq!(IntronicRegion::from_offset(i64::MIN), None);
        assert_eq!(IntronicRegion::from_offset(i64::MAX), None);

        assert_eq!(
            IntronicRegion::from_offset(i64::MIN + 1),
            Some(IntronicRegion::DeepIntronic)
        );
        assert_eq!(
            IntronicRegion::from_offset(i64::MAX - 1),
            Some(IntronicRegion::DeepIntronic)
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
            cds_start_incomplete: false,
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("A".repeat(200)),
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
            protein_id: None,
            exon_cigars: Vec::new(),
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

    /// The offset sanity limit is `offset.unsigned_abs() > 100_000`, so
    /// 100_000 is legal on both the positive and negative side. The negative
    /// side is the routine acceptor-side case, and it needs its own
    /// assertion: the existing positive-only pins leave the magnitude check
    /// free to be dropped without failing the suite, since a large negative
    /// offset would then sail through unchecked.
    ///
    /// The uncertain-offset sentinels `-?` / `+?` (`i64::MIN` / `i64::MAX`,
    /// declared in `src/hgvs/parser/position.rs`) are NOT pinned here or in
    /// `splice_ladder_boundaries`.
    ///
    /// This paragraph previously recorded that as an OPEN production gap —
    /// that `i64::MIN.abs()` overflows, so `validate_tx_pos`,
    /// `IntronicConsequence::from_{cds,tx}_pos` and
    /// `IntronicRegion::from_offset` all panicked in debug on `-?` and filed
    /// `+?` as deep-intronic, none of them calling `has_unknown_offset`
    /// (issue #1087). That is no longer true: `validate_tx_pos` and both
    /// `IntronicConsequence` constructors decline the sentinels via
    /// `pos.has_unknown_offset()`, `IntronicRegion::from_offset` declines via
    /// `is_unknown_offset` and returns `None`, and every magnitude check uses
    /// `unsigned_abs()`. The class is pinned by
    /// `unknown_offset_sentinels_are_declined_not_classified` below.
    #[test]
    fn validate_tx_pos_offset_limit_boundary() {
        let tx = make_test_transcript();

        // Carried over from the deleted `test_validate_tx_pos_with_large_offset`
        // (a routine offset on each side of zero, well inside the limit), so no
        // assertion from that test is lost.
        assert!(
            validate_tx_pos(&TxPos::with_offset(50, 1000), &tx).is_ok(),
            "offset 1000 is well within the limit"
        );
        assert!(
            validate_tx_pos(&TxPos::with_offset(50, -1000), &tx).is_ok(),
            "offset -1000 is well within the limit"
        );

        let at_limit = TxPos::with_offset(50, 100_000);
        assert!(
            validate_tx_pos(&at_limit, &tx).is_ok(),
            "offset 100_000 is at the limit"
        );

        let past_limit = TxPos::with_offset(50, 100_001);
        assert!(
            validate_tx_pos(&past_limit, &tx).is_err(),
            "offset 100_001 is past it"
        );

        let at_negative_limit = TxPos::with_offset(50, -100_000);
        assert!(
            validate_tx_pos(&at_negative_limit, &tx).is_ok(),
            "offset -100_000 is at the limit"
        );

        let past_negative_limit = TxPos::with_offset(50, -100_001);
        assert!(
            validate_tx_pos(&past_negative_limit, &tx).is_err(),
            "offset -100_001 is past it"
        );
    }

    /// Every rung of all three splice ladders, at its exact boundary and
    /// boundary+1, on both the donor (positive offset) and acceptor
    /// (negative offset) side.
    ///
    /// `IntronicConsequence::from_cds_pos`'s donor chain (`abs_offset <= 2`
    /// / `<= 6` / `<= 20` / `<= 50`) and its acceptor chain (`<= 2` / `<= 12`
    /// / `<= 20` / `<= 50`), together with `::from_tx_pos`'s byte-identical
    /// donor and acceptor chains, are copies of one donor/acceptor ladder:
    /// donor breaks at 2 / 6 / 20 / 50, acceptor at 2 / 12 / 20 / 50 — the
    /// donor and acceptor thresholds differ (6 vs 12), so they are not
    /// mirror images of each other. `IntronicRegion::from_offset`'s
    /// `abs_offset <= 2` / `<= 20` / `<= 50` chain is a third, shorter
    /// ladder that is symmetric in the offset's sign (2 / 20 / 50 on both
    /// sides) and has no 6 or 12 break.
    /// A hand-unrolled `#[test]` per rung would be the fourth near-verbatim
    /// copy of this shape in this file; this table drives all three
    /// functions from one set of rows instead.
    #[test]
    fn splice_ladder_boundaries() {
        use IntronicConsequence::{
            IntronVariant, NearSpliceSiteVariant, SpliceAcceptorRegionVariant,
            SpliceAcceptorVariant, SpliceDonorRegionVariant, SpliceDonorVariant,
            SpliceRegionVariant,
        };
        use IntronicRegion::{
            CanonicalSpliceSite, DeepIntronic, ExtendedSpliceRegion, NearSpliceSite,
        };

        // (offset, expected IntronicConsequence, expected IntronicRegion)
        let rows: &[(i64, IntronicConsequence, IntronicRegion)] = &[
            // Donor side (offset > 0): consequence ladder breaks at 2/6/20/50.
            (1, SpliceDonorVariant, CanonicalSpliceSite),
            (2, SpliceDonorVariant, CanonicalSpliceSite),
            (3, SpliceDonorRegionVariant, ExtendedSpliceRegion),
            (6, SpliceDonorRegionVariant, ExtendedSpliceRegion),
            (7, SpliceRegionVariant, ExtendedSpliceRegion),
            (20, SpliceRegionVariant, ExtendedSpliceRegion),
            (21, NearSpliceSiteVariant, NearSpliceSite),
            (50, NearSpliceSiteVariant, NearSpliceSite),
            (51, IntronVariant, DeepIntronic),
            // Acceptor side (offset < 0): consequence ladder breaks at 2/12/20/50.
            (-1, SpliceAcceptorVariant, CanonicalSpliceSite),
            (-2, SpliceAcceptorVariant, CanonicalSpliceSite),
            (-3, SpliceAcceptorRegionVariant, ExtendedSpliceRegion),
            (-12, SpliceAcceptorRegionVariant, ExtendedSpliceRegion),
            (-13, SpliceRegionVariant, ExtendedSpliceRegion),
            (-20, SpliceRegionVariant, ExtendedSpliceRegion),
            (-21, NearSpliceSiteVariant, NearSpliceSite),
            (-50, NearSpliceSiteVariant, NearSpliceSite),
            (-51, IntronVariant, DeepIntronic),
        ];

        for &(offset, expected_consequence, expected_region) in rows {
            let cds_pos = CdsPos::with_offset(100, offset);
            assert_eq!(
                IntronicConsequence::from_cds_pos(&cds_pos),
                Some(expected_consequence),
                "from_cds_pos at offset {offset}"
            );

            let tx_pos = TxPos::with_offset(50, offset);
            assert_eq!(
                IntronicConsequence::from_tx_pos(&tx_pos),
                Some(expected_consequence),
                "from_tx_pos at offset {offset}"
            );

            // #1841 made `from_offset` total by returning `Option`: the unknown-offset
            // sentinels are outside the enum's domain. Every row here is an ordinary
            // offset, so each still classifies.
            assert_eq!(
                IntronicRegion::from_offset(offset),
                Some(expected_region),
                "from_offset at offset {offset}"
            );
        }
    }

    /// `from_cds_pos` and `from_tx_pos` return `None` for a non-intronic
    /// position (no offset); `splice_ladder_boundaries` never exercises this
    /// branch, since every row there carries an offset.
    ///
    /// `offset: Some(0)` is the second, distinct shape reaching that guard —
    /// both `c.100+0` and `n.50+0` parse — and what it pins is the guard
    /// itself, not the ladder behind it. `is_intronic` is
    /// `offset.is_some() && offset != Some(0)` (defined in
    /// `src/hgvs/location.rs`, once for `CdsPos` and once for `TxPos`), so a
    /// zero offset is not intronic
    /// and both functions decline at `!pos.is_intronic()` before any rung
    /// runs; deleting the `!= Some(0)` clause turns these two rows into
    /// `SpliceAcceptorVariant`.
    ///
    /// It deliberately does NOT pin the donor/acceptor split (`offset > 0`)
    /// against an `>= 0` mutant: a zero offset never reaches that branch, so
    /// `> 0` and `>= 0` are observationally identical and no input can
    /// distinguish them.
    #[test]
    fn intronic_consequence_declines_a_non_intronic_position() {
        assert_eq!(
            IntronicConsequence::from_cds_pos(&CdsPos::new(100)),
            None,
            "a CdsPos with no offset is not intronic"
        );
        assert_eq!(
            IntronicConsequence::from_tx_pos(&TxPos::new(50)),
            None,
            "a TxPos with no offset is not intronic"
        );
        assert_eq!(
            IntronicConsequence::from_cds_pos(&CdsPos::with_offset(100, 0)),
            None,
            "a CdsPos with an explicit zero offset is not intronic"
        );
        assert_eq!(
            IntronicConsequence::from_tx_pos(&TxPos::with_offset(50, 0)),
            None,
            "a TxPos with an explicit zero offset is not intronic"
        );
    }

    /// The uncertain-offset sentinels `-?` / `+?` (`OFFSET_UNKNOWN_NEGATIVE` /
    /// `OFFSET_UNKNOWN_POSITIVE` = `i64::MIN` / `i64::MAX`) are out of the
    /// domain of every intronic classifier on this axis, and each must say so
    /// rather than measure them.
    ///
    /// Two failure modes are pinned at once, both of which this code has
    /// carried historically:
    ///
    /// 1. **Panic.** `i64::MIN.abs()` overflows, so any of these reached with
    ///    a plain `.abs()` panics in debug on `-?`. Reaching the assertion at
    ///    all proves the magnitude checks use `unsigned_abs()`.
    /// 2. **Silent misclassification.** `+?` is `i64::MAX`, which compares
    ///    greater than every ladder rung, so an unguarded classifier files it
    ///    as ordinary deep-intronic — a measured distance where there is no
    ///    distance to measure.
    ///
    /// `IntronicRegion::from_offset` is the #1841 half (it declines via
    /// `is_unknown_offset` and returns `None`); the other three decline via
    /// `pos.has_unknown_offset()` (#1767).
    #[test]
    fn unknown_offset_sentinels_are_declined_not_classified() {
        use crate::hgvs::parser::position::{OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE};

        let tx = make_test_transcript();

        for (label, offset) in [
            ("-?", OFFSET_UNKNOWN_NEGATIVE),
            ("+?", OFFSET_UNKNOWN_POSITIVE),
        ] {
            // Checks the MESSAGE, not just `is_err()`: both sentinels exceed
            // the `unsigned_abs() > 100_000` limit on their own, so an
            // `is_err()`-only assertion stays green with the
            // `has_unknown_offset()` decline deleted — it would merely have
            // been reclassified as a distance that is too large.
            let err = validate_tx_pos(&TxPos::with_offset(50, offset), &tx)
                .expect_err("must be declined, not measured");
            let msg = err.to_string();
            assert!(
                msg.contains("denotes no transcript coordinate"),
                "validate_tx_pos must decline the {label} sentinel AS A SENTINEL; got: {msg}"
            );
            assert!(
                !msg.contains("unreasonably large"),
                "validate_tx_pos must not refuse {label} as a mere magnitude; got: {msg}"
            );
            assert_eq!(
                IntronicConsequence::from_cds_pos(&CdsPos::with_offset(100, offset)),
                None,
                "from_cds_pos must decline the {label} sentinel, not classify it"
            );
            assert_eq!(
                IntronicConsequence::from_tx_pos(&TxPos::with_offset(50, offset)),
                None,
                "from_tx_pos must decline the {label} sentinel, not classify it"
            );
            assert_eq!(
                IntronicRegion::from_offset(offset),
                None,
                "from_offset must decline the {label} sentinel, not classify it"
            );
        }
    }
}
