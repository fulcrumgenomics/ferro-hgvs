//! Position types for different coordinate systems
//!
//! HGVS supports multiple coordinate systems:
//! - Genomic (g.): 1-based positions on a chromosome/contig
//! - Coding (c.): Positions relative to CDS start, with intron offsets
//! - Transcript (n.): Positions on non-coding transcript
//! - RNA (r.): Positions on RNA
//! - Protein (p.): Amino acid positions

use serde::{Deserialize, Serialize};
use std::fmt;

/// Special genome position markers
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SpecialPosition {
    /// p-arm telomere (short arm end)
    Pter,
    /// q-arm telomere (long arm end)
    Qter,
    /// Centromere
    Cen,
}

impl fmt::Display for SpecialPosition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SpecialPosition::Pter => write!(f, "pter"),
            SpecialPosition::Qter => write!(f, "qter"),
            SpecialPosition::Cen => write!(f, "cen"),
        }
    }
}

/// Genomic position (g. coordinates)
///
/// Simple 1-based position on a genomic reference sequence.
/// Can also represent special positions like pter, qter, or cen.
/// Supports offsets for uncertain position notation (e.g., g.12345-? or g.67890+?),
/// using the same syntax as intronic offsets but applied to genomic coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomePos {
    /// 1-based position (0 if special position is set)
    pub base: u64,
    /// Special position marker (pter, qter, cen)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub special: Option<SpecialPosition>,
    /// Optional offset from the base position (e.g., -? for uncertain upstream, +? for uncertain downstream)
    /// Uses i64::MIN for uncertain negative offset (-?) and i64::MAX for uncertain positive offset (+?)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offset: Option<i64>,
}

impl GenomePos {
    pub fn new(base: u64) -> Self {
        Self {
            base,
            special: None,
            offset: None,
        }
    }

    pub fn with_offset(base: u64, offset: i64) -> Self {
        Self {
            base,
            special: None,
            offset: Some(offset),
        }
    }

    /// Create a pter (p-arm telomere) position
    pub fn pter() -> Self {
        Self {
            base: 0,
            special: Some(SpecialPosition::Pter),
            offset: None,
        }
    }

    /// Create a qter (q-arm telomere) position
    pub fn qter() -> Self {
        Self {
            base: 0,
            special: Some(SpecialPosition::Qter),
            offset: None,
        }
    }

    /// Create a cen (centromere) position
    pub fn cen() -> Self {
        Self {
            base: 0,
            special: Some(SpecialPosition::Cen),
            offset: None,
        }
    }

    /// Check if this is a special position (pter, qter, cen)
    pub fn is_special(&self) -> bool {
        self.special.is_some()
    }
}

/// Sentinel value for unknown positive offset (+?)
pub const GENOME_OFFSET_UNKNOWN_POSITIVE: i64 = i64::MAX;
/// Sentinel value for unknown negative offset (-?)
pub const GENOME_OFFSET_UNKNOWN_NEGATIVE: i64 = i64::MIN;

impl fmt::Display for GenomePos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(special) = &self.special {
            write!(f, "{}", special)
        } else {
            write!(f, "{}", self.base)?;
            if let Some(offset) = self.offset {
                if offset == GENOME_OFFSET_UNKNOWN_POSITIVE {
                    write!(f, "+?")?;
                } else if offset == GENOME_OFFSET_UNKNOWN_NEGATIVE {
                    write!(f, "-?")?;
                } else if offset > 0 {
                    write!(f, "+{}", offset)?;
                } else if offset < 0 {
                    write!(f, "{}", offset)?;
                }
            }
            Ok(())
        }
    }
}

/// CDS position (c. coordinates)
///
/// Position relative to the start of the coding sequence.
/// Can include intronic offsets (e.g., c.100+5 or c.100-10).
/// Negative positions are upstream of CDS start.
/// Positions with * prefix are downstream of stop codon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CdsPos {
    /// Base position in CDS (can be negative for 5' UTR)
    pub base: i64,
    /// Intronic offset (+ for downstream, - for upstream of exon)
    pub offset: Option<i64>,
    /// Whether this is a 3' UTR position (uses * notation)
    pub utr3: bool,
}

/// Sentinel value for unknown CDS base position (?)
pub const CDS_BASE_UNKNOWN: i64 = 0;

impl CdsPos {
    /// Create a simple exonic CDS position
    pub fn new(base: i64) -> Self {
        Self {
            base,
            offset: None,
            utr3: false,
        }
    }

    /// Create an unknown position (represented as ? in HGVS)
    /// Optionally with an offset (e.g., ?-232)
    pub fn unknown(offset: Option<i64>) -> Self {
        Self {
            base: CDS_BASE_UNKNOWN,
            offset,
            utr3: false,
        }
    }

    /// Check if this is an unknown position (?)
    pub fn is_unknown(&self) -> bool {
        self.base == CDS_BASE_UNKNOWN && !self.utr3
    }

    /// Create a CDS position with intronic offset
    pub fn with_offset(base: i64, offset: i64) -> Self {
        Self {
            base,
            offset: Some(offset),
            utr3: false,
        }
    }

    /// Create a 3' UTR position
    pub fn utr3(base: i64) -> Self {
        Self {
            base,
            offset: None,
            utr3: true,
        }
    }

    /// Check if this position is intronic
    pub fn is_intronic(&self) -> bool {
        self.offset.is_some() && self.offset != Some(0)
    }

    /// Check if this position is in 5' UTR
    pub fn is_5utr(&self) -> bool {
        !self.utr3 && self.base < 1 && self.offset.is_none()
    }

    /// Check if this position is in 3' UTR
    pub fn is_3utr(&self) -> bool {
        self.utr3
    }
}

impl fmt::Display for CdsPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.utr3 {
            write!(f, "*{}", self.base)?;
        } else if self.base == CDS_BASE_UNKNOWN {
            // Unknown position
            write!(f, "?")?;
        } else {
            write!(f, "{}", self.base)?;
        }
        if let Some(offset) = self.offset {
            // Handle sentinel values for uncertain offsets
            if offset == i64::MAX {
                write!(f, "+?")?;
            } else if offset == i64::MIN {
                write!(f, "-?")?;
            } else if offset >= 0 {
                write!(f, "+{}", offset)?;
            } else {
                write!(f, "{}", offset)?;
            }
        }
        Ok(())
    }
}

/// Transcript position (n. coordinates)
///
/// Position on a non-coding transcript. Negative positions represent
/// positions upstream of the transcript start (e.g., n.-30 is 30 bases
/// before the transcript start). Downstream positions use * notation
/// (e.g., n.*5 is 5 bases after the transcript end).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct TxPos {
    /// Position on transcript (can be negative for upstream positions)
    pub base: i64,
    /// Intronic offset
    pub offset: Option<i64>,
    /// Whether this is a downstream position (uses * notation)
    pub downstream: bool,
}

impl TxPos {
    pub fn new(base: i64) -> Self {
        Self {
            base,
            offset: None,
            downstream: false,
        }
    }

    pub fn with_offset(base: i64, offset: i64) -> Self {
        Self {
            base,
            offset: Some(offset),
            downstream: false,
        }
    }

    /// Create a downstream position (n.*5 notation)
    pub fn downstream(base: i64) -> Self {
        Self {
            base,
            offset: None,
            downstream: true,
        }
    }

    /// Create a downstream position with offset (n.*5+10 notation)
    pub fn downstream_with_offset(base: i64, offset: i64) -> Self {
        Self {
            base,
            offset: Some(offset),
            downstream: true,
        }
    }

    pub fn is_intronic(&self) -> bool {
        self.offset.is_some() && self.offset != Some(0)
    }

    /// Check if this is an upstream position (negative base)
    pub fn is_upstream(&self) -> bool {
        self.base < 0 && !self.downstream
    }

    /// Check if this is a downstream position (uses * notation)
    pub fn is_downstream(&self) -> bool {
        self.downstream
    }
}

impl fmt::Display for TxPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.downstream {
            write!(f, "*{}", self.base)?;
        } else {
            write!(f, "{}", self.base)?;
        }
        if let Some(offset) = self.offset {
            if offset >= 0 {
                write!(f, "+{}", offset)?;
            } else {
                write!(f, "{}", offset)?;
            }
        }
        Ok(())
    }
}

/// RNA position (r. coordinates)
///
/// Position on an RNA sequence (lowercase nucleotides).
/// Similar to CdsPos, supports 5' UTR (negative base) and 3' UTR (*base) positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RnaPos {
    /// Position relative to start codon (can be negative for 5' UTR)
    pub base: i64,
    /// Intronic offset
    pub offset: Option<i64>,
    /// True if position is in 3' UTR (uses *N notation)
    pub utr3: bool,
}

impl RnaPos {
    /// Create a simple RNA position
    pub fn new(base: i64) -> Self {
        Self {
            base,
            offset: None,
            utr3: false,
        }
    }

    /// Create an RNA position with intronic offset
    pub fn with_offset(base: i64, offset: i64) -> Self {
        Self {
            base,
            offset: Some(offset),
            utr3: false,
        }
    }

    /// Create a 3' UTR position
    pub fn utr3(base: i64) -> Self {
        Self {
            base,
            offset: None,
            utr3: true,
        }
    }

    /// Check if this position is intronic
    pub fn is_intronic(&self) -> bool {
        self.offset.is_some() && self.offset != Some(0)
    }

    /// Check if this position is in 5' UTR
    pub fn is_5utr(&self) -> bool {
        !self.utr3 && self.base < 1 && self.offset.is_none()
    }

    /// Check if this position is in 3' UTR
    pub fn is_3utr(&self) -> bool {
        self.utr3
    }
}

impl fmt::Display for RnaPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.utr3 {
            write!(f, "*{}", self.base)?;
        } else {
            write!(f, "{}", self.base)?;
        }
        if let Some(offset) = self.offset {
            if offset >= 0 {
                write!(f, "+{}", offset)?;
            } else {
                write!(f, "{}", offset)?;
            }
        }
        Ok(())
    }
}

/// Amino acid enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AminoAcid {
    Ala, // A
    Arg, // R
    Asn, // N
    Asp, // D
    Cys, // C
    Gln, // Q
    Glu, // E
    Gly, // G
    His, // H
    Ile, // I
    Leu, // L
    Lys, // K
    Met, // M
    Phe, // F
    Pro, // P
    Pyl, // O (pyrrolysine)
    Sec, // U (selenocysteine)
    Ser, // S
    Thr, // T
    Trp, // W
    Tyr, // Y
    Val, // V
    Ter, // * (stop codon)
    Xaa, // X (unknown)
}

impl AminoAcid {
    /// Parse from 3-letter code
    pub fn from_three_letter(s: &str) -> Option<Self> {
        match s {
            "Ala" => Some(Self::Ala),
            "Arg" => Some(Self::Arg),
            "Asn" => Some(Self::Asn),
            "Asp" => Some(Self::Asp),
            "Cys" => Some(Self::Cys),
            "Gln" => Some(Self::Gln),
            "Glu" => Some(Self::Glu),
            "Gly" => Some(Self::Gly),
            "His" => Some(Self::His),
            "Ile" => Some(Self::Ile),
            "Leu" => Some(Self::Leu),
            "Lys" => Some(Self::Lys),
            "Met" => Some(Self::Met),
            "Phe" => Some(Self::Phe),
            "Pro" => Some(Self::Pro),
            "Pyl" => Some(Self::Pyl),
            "Sec" => Some(Self::Sec),
            "Ser" => Some(Self::Ser),
            "Thr" => Some(Self::Thr),
            "Trp" => Some(Self::Trp),
            "Tyr" => Some(Self::Tyr),
            "Val" => Some(Self::Val),
            "Ter" => Some(Self::Ter),
            "Xaa" => Some(Self::Xaa),
            _ => None,
        }
    }

    /// Get 3-letter code
    pub fn to_three_letter(&self) -> &'static str {
        match self {
            Self::Ala => "Ala",
            Self::Arg => "Arg",
            Self::Asn => "Asn",
            Self::Asp => "Asp",
            Self::Cys => "Cys",
            Self::Gln => "Gln",
            Self::Glu => "Glu",
            Self::Gly => "Gly",
            Self::His => "His",
            Self::Ile => "Ile",
            Self::Leu => "Leu",
            Self::Lys => "Lys",
            Self::Met => "Met",
            Self::Phe => "Phe",
            Self::Pro => "Pro",
            Self::Pyl => "Pyl",
            Self::Sec => "Sec",
            Self::Ser => "Ser",
            Self::Thr => "Thr",
            Self::Trp => "Trp",
            Self::Tyr => "Tyr",
            Self::Val => "Val",
            Self::Ter => "Ter",
            Self::Xaa => "Xaa",
        }
    }

    /// Get 1-letter code
    pub fn to_one_letter(&self) -> char {
        match self {
            Self::Ala => 'A',
            Self::Arg => 'R',
            Self::Asn => 'N',
            Self::Asp => 'D',
            Self::Cys => 'C',
            Self::Gln => 'Q',
            Self::Glu => 'E',
            Self::Gly => 'G',
            Self::His => 'H',
            Self::Ile => 'I',
            Self::Leu => 'L',
            Self::Lys => 'K',
            Self::Met => 'M',
            Self::Phe => 'F',
            Self::Pro => 'P',
            Self::Pyl => 'O',
            Self::Sec => 'U',
            Self::Ser => 'S',
            Self::Thr => 'T',
            Self::Trp => 'W',
            Self::Tyr => 'Y',
            Self::Val => 'V',
            Self::Ter => '*',
            Self::Xaa => 'X',
        }
    }

    /// Parse from 1-letter code (uppercase only)
    ///
    /// HGVS notation uses uppercase for 1-letter amino acid codes.
    /// Lowercase letters are reserved for other notations like "fs" (frameshift)
    /// and "ext" (extension).
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::hgvs::location::AminoAcid;
    ///
    /// assert_eq!(AminoAcid::from_one_letter('V'), Some(AminoAcid::Val));
    /// assert_eq!(AminoAcid::from_one_letter('v'), None); // lowercase not accepted
    /// assert_eq!(AminoAcid::from_one_letter('*'), Some(AminoAcid::Ter));
    /// ```
    pub fn from_one_letter(c: char) -> Option<Self> {
        match c {
            'A' => Some(Self::Ala),
            'R' => Some(Self::Arg),
            'N' => Some(Self::Asn),
            'D' => Some(Self::Asp),
            'C' => Some(Self::Cys),
            'Q' => Some(Self::Gln),
            'E' => Some(Self::Glu),
            'G' => Some(Self::Gly),
            'H' => Some(Self::His),
            'I' => Some(Self::Ile),
            'L' => Some(Self::Leu),
            'K' => Some(Self::Lys),
            'M' => Some(Self::Met),
            'F' => Some(Self::Phe),
            'O' => Some(Self::Pyl),
            'P' => Some(Self::Pro),
            'U' => Some(Self::Sec),
            'S' => Some(Self::Ser),
            'T' => Some(Self::Thr),
            'W' => Some(Self::Trp),
            'Y' => Some(Self::Tyr),
            'V' => Some(Self::Val),
            '*' => Some(Self::Ter),
            'X' => Some(Self::Xaa),
            _ => None,
        }
    }
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_three_letter())
    }
}

/// Protein position (p. coordinates)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ProtPos {
    /// Amino acid at this position
    pub aa: AminoAcid,
    /// 1-based position in protein
    pub number: u64,
}

impl ProtPos {
    pub fn new(aa: AminoAcid, number: u64) -> Self {
        Self { aa, number }
    }
}

impl fmt::Display for ProtPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.aa, self.number)
    }
}

/// IVS (Intervening Sequence) position notation
///
/// Older notation for intronic positions that's still commonly used in
/// clinical settings. e.g., IVS1+5 means position +5 in intron 1.
/// This is equivalent to c.N+5 where N is the last base of the upstream exon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct IvsPos {
    /// Intron number (1-based)
    pub intron: u32,
    /// Offset from the boundary (positive for 5' end, negative for 3' end)
    pub offset: i64,
}

impl IvsPos {
    /// Create a new IVS position
    pub fn new(intron: u32, offset: i64) -> Self {
        Self { intron, offset }
    }

    /// Check if this is a 5' boundary position (positive offset)
    pub fn is_5prime(&self) -> bool {
        self.offset > 0
    }

    /// Check if this is a 3' boundary position (negative offset)
    pub fn is_3prime(&self) -> bool {
        self.offset < 0
    }

    /// Convert to CDS position if the transcript boundary is known
    ///
    /// Requires the transcript position of the exon boundary.
    pub fn to_cds_pos(&self, boundary_cds_pos: i64, is_utr3: bool) -> CdsPos {
        CdsPos {
            base: boundary_cds_pos,
            offset: Some(self.offset),
            utr3: is_utr3,
        }
    }

    /// Check if this is a deep intronic position (>50bp from exon)
    pub fn is_deep_intronic(&self) -> bool {
        self.offset.abs() > 50
    }

    /// Check if this is at a canonical splice site
    pub fn is_canonical_splice_site(&self) -> bool {
        self.offset.abs() <= 2
    }
}

impl fmt::Display for IvsPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.offset >= 0 {
            write!(f, "IVS{}+{}", self.intron, self.offset)
        } else {
            write!(f, "IVS{}{}", self.intron, self.offset)
        }
    }
}

/// Extension trait for CdsPos to convert to/from IVS notation
pub trait IvsNotation {
    /// Convert to IVS notation if this is an intronic position
    fn to_ivs(&self, intron_number: u32) -> Option<IvsPos>;

    /// Check if this position can be represented in IVS notation
    fn has_ivs_notation(&self) -> bool;
}

impl IvsNotation for CdsPos {
    fn to_ivs(&self, intron_number: u32) -> Option<IvsPos> {
        self.offset.map(|offset| IvsPos::new(intron_number, offset))
    }

    fn has_ivs_notation(&self) -> bool {
        self.is_intronic()
    }
}

impl IvsNotation for TxPos {
    fn to_ivs(&self, intron_number: u32) -> Option<IvsPos> {
        self.offset.map(|offset| IvsPos::new(intron_number, offset))
    }

    fn has_ivs_notation(&self) -> bool {
        self.is_intronic()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genome_pos_display() {
        let pos = GenomePos::new(12345);
        assert_eq!(format!("{}", pos), "12345");
    }

    #[test]
    fn test_cds_pos_simple() {
        let pos = CdsPos::new(100);
        assert_eq!(format!("{}", pos), "100");
        assert!(!pos.is_intronic());
    }

    #[test]
    fn test_cds_pos_with_offset() {
        let pos = CdsPos::with_offset(100, 5);
        assert_eq!(format!("{}", pos), "100+5");
        assert!(pos.is_intronic());

        let pos = CdsPos::with_offset(100, -5);
        assert_eq!(format!("{}", pos), "100-5");
        assert!(pos.is_intronic());
    }

    #[test]
    fn test_cds_pos_utr3() {
        let pos = CdsPos::utr3(50);
        assert_eq!(format!("{}", pos), "*50");
        assert!(pos.is_3utr());
    }

    #[test]
    fn test_cds_pos_5utr() {
        let pos = CdsPos::new(-10);
        assert!(pos.is_5utr());
    }

    #[test]
    fn test_prot_pos_display() {
        let pos = ProtPos::new(AminoAcid::Met, 1);
        assert_eq!(format!("{}", pos), "Met1");
    }

    #[test]
    fn test_amino_acid_codes() {
        let aa = AminoAcid::Met;
        assert_eq!(aa.to_three_letter(), "Met");
        assert_eq!(aa.to_one_letter(), 'M');
        assert_eq!(AminoAcid::from_three_letter("Met"), Some(AminoAcid::Met));
    }

    #[test]
    fn test_ivs_pos_new() {
        let pos = IvsPos::new(1, 5);
        assert_eq!(pos.intron, 1);
        assert_eq!(pos.offset, 5);
        assert!(pos.is_5prime());
        assert!(!pos.is_3prime());
    }

    #[test]
    fn test_ivs_pos_display() {
        // 5' boundary notation
        let pos = IvsPos::new(1, 5);
        assert_eq!(format!("{}", pos), "IVS1+5");

        // 3' boundary notation
        let pos = IvsPos::new(2, -10);
        assert_eq!(format!("{}", pos), "IVS2-10");
    }

    #[test]
    fn test_ivs_pos_to_cds_pos() {
        let ivs = IvsPos::new(1, 5);
        let cds = ivs.to_cds_pos(100, false);

        assert_eq!(cds.base, 100);
        assert_eq!(cds.offset, Some(5));
        assert!(!cds.utr3);
    }

    #[test]
    fn test_ivs_pos_deep_intronic() {
        let shallow = IvsPos::new(1, 10);
        assert!(!shallow.is_deep_intronic());
        assert!(!shallow.is_canonical_splice_site());

        let canonical = IvsPos::new(1, 2);
        assert!(canonical.is_canonical_splice_site());

        let deep = IvsPos::new(1, 100);
        assert!(deep.is_deep_intronic());
    }

    #[test]
    fn test_cds_pos_to_ivs() {
        let cds = CdsPos::with_offset(100, 5);
        let ivs = cds.to_ivs(1);

        assert!(ivs.is_some());
        let ivs = ivs.unwrap();
        assert_eq!(ivs.intron, 1);
        assert_eq!(ivs.offset, 5);
    }

    #[test]
    fn test_cds_pos_to_ivs_non_intronic() {
        let cds = CdsPos::new(100); // No offset
        let ivs = cds.to_ivs(1);

        assert!(ivs.is_none());
    }

    #[test]
    fn test_ivs_notation_trait() {
        let intronic_cds = CdsPos::with_offset(100, 5);
        assert!(intronic_cds.has_ivs_notation());

        let exonic_cds = CdsPos::new(100);
        assert!(!exonic_cds.has_ivs_notation());

        let intronic_tx = TxPos::with_offset(100, -10);
        assert!(intronic_tx.has_ivs_notation());
    }
}
