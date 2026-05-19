//! Boundary detection for normalization
//!
//! Determines the boundaries within which a variant can be shuffled.

use crate::error::FerroError;
use crate::normalize::config::NormalizeConfig;
use crate::reference::transcript::Transcript;

/// Boundaries for variant shuffling.
///
/// **Coordinate convention.** Half-open `[left, right)` in 0-based
/// coordinates: `left` is the inclusive minimum and `right` is the
/// exclusive maximum. All production constructors in `src/normalize/`
/// build `Boundaries` with this contract — see e.g.
/// `Boundaries::new(0, ref_seq.len())` in `normalize_genome` and the
/// 0-based bounds returned by [`get_cds_boundaries`]. The CDS↔UTR axis
/// clamp added in #337 explicitly relies on this half-open semantic so
/// that `right == cds_end_1b` excludes the first 3'UTR base from a CDS
/// shuffle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Boundaries {
    /// Inclusive 0-based minimum position.
    pub left: u64,
    /// Exclusive 0-based maximum position (one past the last reachable
    /// position).
    pub right: u64,
}

impl Boundaries {
    pub fn new(left: u64, right: u64) -> Self {
        Self { left, right }
    }

    /// Check if a position is within the half-open bounds `[left, right)`.
    /// `right` itself is **not** contained — it's one past the last
    /// reachable position by the half-open convention documented on
    /// [`Boundaries`].
    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.left && pos < self.right
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
/// **Coordinate system.** Returns 0-based bounds matching the
/// [`crate::normalize::shuffle::shuffle`] convention: `left` is the
/// 0-based inclusive minimum reachable start position and `right` is
/// the 0-based exclusive maximum reachable end position. Tx positions
/// on the `Transcript` struct are 1-based-inclusive, so the conversion
/// to 0-based is performed here (`left = tx_start_1b - 1`, `right =
/// tx_end_1b`). Other long-standing call sites that build `Boundaries`
/// inline (`Boundaries::new(1, seq_len)` etc.) inherit the pre-existing
/// 1-based-input quirk; fixing those is out of scope for #337.
///
/// Two boundary kinds are intersected:
///
///   1. **Exon bound.** When `cross_boundaries=false`, the shuffle is
///      confined to the exon containing `tx_pos`; otherwise the bound is
///      the full transcript.
///
///   2. **CDS↔UTR axis bound (always applied; closes #337).** The 5'UTR
///      / CDS / 3'UTR transitions change the HGVS coordinate sub-axis
///      (`c.-N` vs `c.<N>` vs `c.*N`), so a 3'-rule shuffle must not
///      cross them — doing so would silently re-classify the variant
///      onto a different axis. The axis bound depends on which region
///      `tx_pos` lies in:
///
///        - 5'UTR (`tx_pos < cds_start`): tx 1-based `[1, cds_start - 1]`
///        - CDS   (`cds_start <= tx_pos <= cds_end`): `[cds_start, cds_end]`
///        - 3'UTR (`tx_pos > cds_end`): `[cds_end + 1, sequence_length]`
///
///      For non-coding transcripts (no `cds_start`/`cds_end`) the axis
///      bound is the full transcript.
///
/// The returned `Boundaries` is the 0-based intersection. Errors when
/// `tx_pos` falls outside every exon under `cross_boundaries=false`
/// (i.e. is intronic — exonic-shuffle code shouldn't reach here).
pub fn get_cds_boundaries(
    transcript: &Transcript,
    tx_pos: u64,
    config: &NormalizeConfig,
) -> Result<Boundaries, FerroError> {
    Ok(get_cds_boundaries_with_axis_info(transcript, tx_pos, config)?.clamped)
}

/// 0-based shuffle bounds plus the un-clamped exon bound, for callers
/// that need to detect whether the CDS↔UTR axis clamp tightened the
/// shuffle range (closes-after: #349). `clamped` is the axis ∩ exon
/// intersection (what `get_cds_boundaries` returns); `exon` is the
/// exon-only bound on the same exon-spanning rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdsBoundariesWithAxis {
    /// 0-based shuffle bounds = axis ∩ exon. Use this for the shuffle.
    pub clamped: Boundaries,
    /// 0-based exon-only bounds (without the CDS↔UTR axis clamp). Used
    /// only to detect whether the axis clamp was operative.
    pub exon: Boundaries,
    /// Axis region that `tx_pos` lies in: `"5utr"`, `"cds"`, `"3utr"`,
    /// or `"none"` for non-coding transcripts.
    pub axis_region: AxisRegion,
}

/// Coordinate sub-axis a tx-frame position lies in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisRegion {
    /// `tx_pos < cds_start` (HGVS `c.-N`).
    FiveUtr,
    /// `cds_start <= tx_pos <= cds_end` (HGVS `c.<N>`).
    Cds,
    /// `tx_pos > cds_end` (HGVS `c.*N`).
    ThreeUtr,
    /// Non-coding transcript (no `cds_start`/`cds_end`).
    None,
}

/// Like [`get_cds_boundaries`] but also returns the un-clamped exon
/// bound and the axis region of `tx_pos`. Used by `normalize_cds` to
/// emit `AxisClampApplied` / `CrossAxisVariantNotShuffled` warnings.
pub fn get_cds_boundaries_with_axis_info(
    transcript: &Transcript,
    tx_pos: u64,
    config: &NormalizeConfig,
) -> Result<CdsBoundariesWithAxis, FerroError> {
    let to_0based =
        |left_1b: u64, right_1b: u64| -> (u64, u64) { (left_1b.saturating_sub(1), right_1b) };

    let (exon_left, exon_right) = if config.cross_boundaries {
        to_0based(1, transcript.sequence_length())
    } else {
        let mut found = None;
        for exon in &transcript.exons {
            if tx_pos >= exon.start && tx_pos <= exon.end {
                found = Some(to_0based(exon.start, exon.end));
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
    let (axis_left, axis_right, axis_region) = match (transcript.cds_start, transcript.cds_end) {
        (Some(s), Some(e)) if e >= s => {
            if tx_pos < s {
                let (l, r) = to_0based(1, s.saturating_sub(1));
                (l, r, AxisRegion::FiveUtr)
            } else if tx_pos > e {
                let (l, r) = to_0based(e + 1, seq_len);
                (l, r, AxisRegion::ThreeUtr)
            } else {
                let (l, r) = to_0based(s, e);
                (l, r, AxisRegion::Cds)
            }
        }
        _ => {
            let (l, r) = to_0based(1, seq_len);
            (l, r, AxisRegion::None)
        }
    };

    Ok(CdsBoundariesWithAxis {
        clamped: Boundaries::new(axis_left.max(exon_left), axis_right.min(exon_right)),
        exon: Boundaries::new(exon_left, exon_right),
        axis_region,
    })
}

/// Resolve `tx_pos` to its coordinate sub-axis on `transcript`. Used
/// alongside [`get_cds_boundaries_with_axis_info`] when callers only
/// need the axis of a position (e.g. cross-axis detection in #350).
pub fn axis_region_of(transcript: &Transcript, tx_pos: u64) -> AxisRegion {
    match (transcript.cds_start, transcript.cds_end) {
        (Some(s), Some(e)) if e >= s => {
            if tx_pos < s {
                AxisRegion::FiveUtr
            } else if tx_pos > e {
                AxisRegion::ThreeUtr
            } else {
                AxisRegion::Cds
            }
        }
        _ => AxisRegion::None,
    }
}

impl AxisRegion {
    /// Lowercase short label for warning fields.
    pub fn label(&self) -> &'static str {
        match self {
            AxisRegion::FiveUtr => "5utr",
            AxisRegion::Cds => "cds",
            AxisRegion::ThreeUtr => "3utr",
            AxisRegion::None => "none",
        }
    }
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
        // Half-open `[10, 20)`: 10 is in, 19 is in (last reachable),
        // 20 is the exclusive upper bound (NOT contained).
        let bounds = Boundaries::new(10, 20);
        assert!(bounds.contains(10));
        assert!(bounds.contains(15));
        assert!(bounds.contains(19));
        assert!(!bounds.contains(20));
        assert!(!bounds.contains(9));
        assert!(!bounds.contains(21));
    }

    #[test]
    fn test_cds_boundaries_within_exon() {
        let transcript = make_test_transcript();
        let config = NormalizeConfig::default();

        // Exon 1 is 1-based [1, 4]; converted to 0-based-inclusive /
        // 0-based-exclusive boundaries that's (0, 4).
        let bounds = get_cds_boundaries(&transcript, 2, &config).unwrap();
        assert_eq!(bounds.left, 0);
        assert_eq!(bounds.right, 4);
    }

    #[test]
    fn test_cds_boundaries_cross_allowed() {
        let transcript = make_test_transcript();
        let config = NormalizeConfig::default().allow_crossing_boundaries();

        // make_test_transcript() has cds_start=1, cds_end=12 — the entire
        // transcript IS the CDS. Axis bound = exon bound = full transcript,
        // converted to 0-based that's (0, 12).
        let bounds = get_cds_boundaries(&transcript, 2, &config).unwrap();
        assert_eq!(bounds.left, 0);
        assert_eq!(bounds.right, 12);
    }
}
