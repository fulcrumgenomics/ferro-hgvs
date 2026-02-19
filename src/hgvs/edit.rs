//! Edit types for nucleic acid and protein variants
//!
//! Edits describe the actual change (substitution, deletion, insertion, etc.)

use super::location::AminoAcid;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// DNA/RNA nucleotide base (including IUPAC ambiguity codes)
///
/// Uses `#[repr(u8)]` with ASCII discriminants for zero-cost `as u8` conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(u8)]
pub enum Base {
    // Standard bases
    A = b'A',
    C = b'C',
    G = b'G',
    T = b'T',
    U = b'U', // RNA

    // IUPAC ambiguity codes
    R = b'R', // A or G (puRine)
    Y = b'Y', // C or T (pYrimidine)
    S = b'S', // G or C (Strong)
    W = b'W', // A or T (Weak)
    K = b'K', // G or T (Keto)
    M = b'M', // A or C (aMino)
    B = b'B', // C, G, or T (not A)
    D = b'D', // A, G, or T (not C)
    H = b'H', // A, C, or T (not G)
    V = b'V', // A, C, or G (not T)
    N = b'N', // Any base
}

impl Base {
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            // Standard bases
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' => Some(Base::T),
            'U' => Some(Base::U),
            // IUPAC ambiguity codes
            'R' => Some(Base::R),
            'Y' => Some(Base::Y),
            'S' => Some(Base::S),
            'W' => Some(Base::W),
            'K' => Some(Base::K),
            'M' => Some(Base::M),
            'B' => Some(Base::B),
            'D' => Some(Base::D),
            'H' => Some(Base::H),
            'V' => Some(Base::V),
            'N' => Some(Base::N),
            _ => None,
        }
    }

    /// Convert to ASCII character. Zero-cost due to `#[repr(u8)]`.
    #[inline]
    pub fn to_char(self) -> char {
        // SAFETY: All discriminants are valid ASCII uppercase letters
        self as u8 as char
    }

    /// Convert to ASCII byte. Zero-cost due to `#[repr(u8)]`.
    #[inline]
    pub fn to_u8(self) -> u8 {
        self as u8
    }

    /// Convert to lowercase ASCII character (for RNA display).
    #[inline]
    pub fn to_lowercase_char(self) -> char {
        (self as u8).to_ascii_lowercase() as char
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// Nucleotide sequence
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Sequence(pub Vec<Base>);

impl Sequence {
    pub fn new(bases: Vec<Base>) -> Self {
        Self(bases)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn bases(&self) -> &[Base] {
        &self.0
    }

    /// Format as lowercase string (for RNA display).
    pub fn to_lowercase_string(&self) -> String {
        self.0.iter().map(|b| b.to_lowercase_char()).collect()
    }
}

impl FromStr for Sequence {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let bases: Option<Vec<Base>> = s.chars().map(Base::from_char).collect();
        bases.map(Self).ok_or("Invalid base character")
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in &self.0 {
            write!(f, "{}", base)?;
        }
        Ok(())
    }
}

/// Repeat count for repeat variants
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RepeatCount {
    /// Exact count (e.g., [12])
    Exact(u64),
    /// Range of counts (e.g., [10_15])
    Range(u64, u64),
    /// Lower bound with uncertain upper (e.g., [10_?])
    MinUncertain(u64),
    /// Upper bound with uncertain lower (e.g., [?_20])
    MaxUncertain(u64),
    /// Fully uncertain (e.g., [?_?] or [?])
    Unknown,
}

/// A single unit in a complex tandem repeat (sequence + count)
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct RepeatUnit {
    /// The repeated sequence
    pub sequence: Sequence,
    /// The repeat count
    pub count: RepeatCount,
}

impl fmt::Display for RepeatUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.sequence, self.count)
    }
}

impl fmt::Display for RepeatCount {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RepeatCount::Exact(n) => write!(f, "[{}]", n),
            RepeatCount::Range(min, max) => write!(f, "[{}_{}]", min, max),
            RepeatCount::MinUncertain(min) => write!(f, "[{}_?]", min),
            RepeatCount::MaxUncertain(max) => write!(f, "[?_{}]", max),
            RepeatCount::Unknown => write!(f, "[?]"),
        }
    }
}

/// Inserted sequence specification
///
/// Represents the different ways an inserted sequence can be specified in HGVS notation.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum InsertedSequence {
    /// Literal sequence (e.g., insATG)
    Literal(Sequence),
    /// Unknown count of unspecified bases (e.g., ins10, ins(10))
    Count(u64),
    /// Unknown range of unspecified bases (e.g., ins(10_20))
    Range(u64, u64),
    /// Repeated base with count (e.g., insA[10], insN[15], insN[15_30])
    Repeat { base: Base, count: RepeatCount },
    /// Repeated sequence with count (e.g., delinsTCGGCAGCGGCACAGCGAGG[13])
    /// Used when a multi-base sequence is repeated a specified number of times
    SequenceRepeat {
        sequence: Sequence,
        count: RepeatCount,
    },
    /// Complex insertion with multiple parts (e.g., ins[A[10];T])
    Complex(Vec<InsertedPart>),
    /// Named sequence element (e.g., insAluYb8, delinsAluYa5)
    Named(String),
    /// Reference to another location (e.g., delinsNC_000022.11:g.100_200)
    /// This is used for conversions specified with delins syntax
    Reference(String),
    /// Position range reference in same sequence (e.g., delins100833642_100833702)
    /// This refers to positions in the same reference sequence
    PositionRange { start: u64, end: u64 },
    /// Position range with inversion (e.g., delins86116_86422inv)
    PositionRangeInv { start: u64, end: u64 },
    /// Uncertain sequence (e.g., delins(?))
    Uncertain,
    /// Empty insertion (e.g., just "ins" without sequence - rare but valid)
    Empty,
}

/// Part of a complex inserted sequence
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum InsertedPart {
    /// Literal sequence part
    Literal(Sequence),
    /// Repeated base part (e.g., A[10], N[15_30])
    Repeat { base: Base, count: RepeatCount },
    /// Position range reference (e.g., 401_419 in c.419_420ins[T;401_419])
    PositionRange { start: u64, end: u64 },
    /// Position range with inversion (e.g., 45043709_45310738inv)
    PositionRangeInv { start: u64, end: u64 },
    /// CDS position range with intronic offsets (e.g., 244-8_249 in c.249_250ins[N[2800];244-8_249])
    /// Stored as string to preserve offset notation
    CdsPositionRange(String),
    /// External sequence reference (e.g., KT192064.1:1_310 or NC_000020.11:g.2823027_2826302)
    ExternalRef(String),
}

impl InsertedPart {
    /// Format with lowercase nucleotides (for RNA display).
    pub fn to_rna_string(&self) -> String {
        match self {
            InsertedPart::Literal(seq) => seq.to_lowercase_string(),
            InsertedPart::Repeat { base, count } => {
                format!("{}{}", base.to_lowercase_char(), count)
            }
            InsertedPart::PositionRange { start, end } => format!("{}_{}", start, end),
            InsertedPart::PositionRangeInv { start, end } => format!("{}_{}inv", start, end),
            InsertedPart::CdsPositionRange(range) => range.clone(),
            InsertedPart::ExternalRef(reference) => reference.clone(),
        }
    }
}

impl fmt::Display for InsertedPart {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InsertedPart::Literal(seq) => write!(f, "{}", seq),
            InsertedPart::Repeat { base, count } => write!(f, "{}{}", base, count),
            InsertedPart::PositionRange { start, end } => write!(f, "{}_{}", start, end),
            InsertedPart::PositionRangeInv { start, end } => write!(f, "{}_{}inv", start, end),
            InsertedPart::CdsPositionRange(range) => write!(f, "{}", range),
            InsertedPart::ExternalRef(reference) => write!(f, "{}", reference),
        }
    }
}

impl fmt::Display for InsertedSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InsertedSequence::Literal(seq) => write!(f, "{}", seq),
            InsertedSequence::Count(n) => write!(f, "{}", n),
            InsertedSequence::Range(min, max) => write!(f, "({}_{})", min, max),
            InsertedSequence::Repeat { base, count } => write!(f, "{}{}", base, count),
            InsertedSequence::SequenceRepeat { sequence, count } => {
                write!(f, "{}{}", sequence, count)
            }
            InsertedSequence::Complex(parts) => {
                write!(f, "[")?;
                for (i, part) in parts.iter().enumerate() {
                    if i > 0 {
                        write!(f, ";")?;
                    }
                    write!(f, "{}", part)?;
                }
                write!(f, "]")
            }
            InsertedSequence::Named(name) => write!(f, "{}", name),
            InsertedSequence::Reference(reference) => write!(f, "[{}]", reference),
            InsertedSequence::PositionRange { start, end } => write!(f, "{}_{}", start, end),
            InsertedSequence::PositionRangeInv { start, end } => {
                write!(f, "{}_{}inv", start, end)
            }
            InsertedSequence::Uncertain => write!(f, "(?)"),
            InsertedSequence::Empty => Ok(()),
        }
    }
}

impl InsertedSequence {
    /// Create a literal sequence insertion
    pub fn literal(seq: Sequence) -> Self {
        InsertedSequence::Literal(seq)
    }

    /// Check if this is a literal sequence
    pub fn is_literal(&self) -> bool {
        matches!(self, InsertedSequence::Literal(_))
    }

    /// Get the literal sequence if this is one
    pub fn as_literal(&self) -> Option<&Sequence> {
        match self {
            InsertedSequence::Literal(seq) => Some(seq),
            _ => None,
        }
    }

    /// Get the bases if this is a literal sequence
    /// Returns None for count, range, repeat, complex, or empty variants
    pub fn bases(&self) -> Option<&[Base]> {
        match self {
            InsertedSequence::Literal(seq) => Some(seq.bases()),
            _ => None,
        }
    }

    /// Get the length of this inserted sequence, if known
    pub fn len(&self) -> Option<usize> {
        match self {
            InsertedSequence::Literal(seq) => Some(seq.len()),
            InsertedSequence::Count(n) => Some(*n as usize),
            InsertedSequence::Repeat { count, .. } => match count {
                RepeatCount::Exact(n) => Some(*n as usize),
                _ => None, // Unknown exact length for ranges/uncertain
            },
            InsertedSequence::SequenceRepeat { sequence, count } => match count {
                RepeatCount::Exact(n) => Some(sequence.len() * (*n as usize)),
                _ => None, // Unknown exact length for ranges/uncertain
            },
            InsertedSequence::Range(_, _) => None, // Unknown exact length
            InsertedSequence::Complex(_) => None,  // Would need to sum parts
            InsertedSequence::Named(_) => None,    // Unknown length for named elements
            InsertedSequence::Reference(_) => None, // Requires lookup
            InsertedSequence::PositionRange { start, end } => Some((end - start + 1) as usize),
            InsertedSequence::PositionRangeInv { start, end } => Some((end - start + 1) as usize),
            InsertedSequence::Uncertain => None, // Unknown
            InsertedSequence::Empty => Some(0),
        }
    }

    /// Check if this is an empty insertion
    pub fn is_empty(&self) -> bool {
        match self {
            InsertedSequence::Literal(seq) => seq.is_empty(),
            InsertedSequence::Named(name) => name.is_empty(),
            InsertedSequence::Empty => true,
            _ => false,
        }
    }

    /// Format with lowercase nucleotides (for RNA display).
    pub fn to_rna_string(&self) -> String {
        match self {
            InsertedSequence::Literal(seq) => seq.to_lowercase_string(),
            InsertedSequence::Count(n) => n.to_string(),
            InsertedSequence::Range(min, max) => format!("({}_{})", min, max),
            InsertedSequence::Repeat { base, count } => {
                format!("{}{}", base.to_lowercase_char(), count)
            }
            InsertedSequence::SequenceRepeat { sequence, count } => {
                format!("{}{}", sequence.to_lowercase_string(), count)
            }
            InsertedSequence::Complex(parts) => {
                let mut s = String::from("[");
                for (i, part) in parts.iter().enumerate() {
                    if i > 0 {
                        s.push(';');
                    }
                    s.push_str(&part.to_rna_string());
                }
                s.push(']');
                s
            }
            InsertedSequence::Named(name) => name.clone(),
            InsertedSequence::Reference(reference) => format!("[{}]", reference),
            InsertedSequence::PositionRange { start, end } => format!("{}_{}", start, end),
            InsertedSequence::PositionRangeInv { start, end } => format!("{}_{}inv", start, end),
            InsertedSequence::Uncertain => String::from("(?)"),
            InsertedSequence::Empty => String::new(),
        }
    }
}

/// Methylation status for epigenetic variants
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MethylationStatus {
    /// Gain of methylation (|gom)
    GainOfMethylation,
    /// Loss of methylation (|lom)
    LossOfMethylation,
    /// Methylation unchanged (|met=)
    Unchanged,
}

impl fmt::Display for MethylationStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MethylationStatus::GainOfMethylation => write!(f, "|gom"),
            MethylationStatus::LossOfMethylation => write!(f, "|lom"),
            MethylationStatus::Unchanged => write!(f, "|met="),
        }
    }
}

/// Uncertain duplication extent
///
/// Used for duplications with uncertain extent specification (e.g., dup(731_741), dup?, dup(?))
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum UncertainDupExtent {
    /// Completely unknown extent (dup? or dup(?))
    Unknown,
    /// Range of positions specifying uncertain extent (dup(731_741))
    Range(u64, u64),
}

impl fmt::Display for UncertainDupExtent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            UncertainDupExtent::Unknown => write!(f, "?"),
            UncertainDupExtent::Range(start, end) => write!(f, "({}_{})", start, end),
        }
    }
}

/// Nucleic acid edit types
///
/// Represents the different types of edits that can occur in DNA/RNA sequences.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NaEdit {
    /// Substitution: single base change (e.g., A>G)
    Substitution { reference: Base, alternative: Base },

    /// Substitution without reference base (e.g., >A, >G)
    /// This non-standard notation specifies only the alternative allele
    SubstitutionNoRef { alternative: Base },

    /// Deletion: removal of bases (e.g., del, delA, del3, del101)
    Deletion {
        sequence: Option<Sequence>,
        /// Explicit length of deleted region (e.g., del101 means 101 bases deleted)
        length: Option<u64>,
    },

    /// Insertion: addition of bases between two positions (e.g., insATG, ins10, insA[10])
    Insertion { sequence: InsertedSequence },

    /// Deletion-insertion: replacement of bases (e.g., delinsATG, delinsA[10])
    Delins { sequence: InsertedSequence },

    /// Duplication: copy of bases (e.g., dup, dupATG, dup101, dup(731_741), dup?)
    Duplication {
        sequence: Option<Sequence>,
        /// Explicit length of duplicated region (e.g., dup101 means 101 bases duplicated)
        length: Option<u64>,
        /// Uncertain extent of duplicated region (e.g., dup(731_741) or dup?)
        uncertain_extent: Option<UncertainDupExtent>,
    },

    /// Duplication-Insertion: duplication followed by insertion (e.g., dupinsCTCA)
    /// This is a non-standard but commonly used notation in databases like ClinVar
    DupIns { sequence: InsertedSequence },

    /// Inversion: reversal of bases (e.g., inv, inv3, invATG)
    Inversion {
        /// Optional sequence that was inverted
        sequence: Option<Sequence>,
        /// Explicit length of inverted region (e.g., inv3)
        length: Option<u64>,
    },

    /// Repeat: short tandem repeat (e.g., [12] or [10_15] or [3]T or A[6][1])
    /// VEP sometimes outputs repeats with trailing nucleotides (e.g., c.212-18[3]T)
    /// Genotype notation uses multiple counts (e.g., A[6][1] for heterozygous)
    Repeat {
        sequence: Option<Sequence>,
        count: RepeatCount,
        /// Additional repeat counts for genotype notation (e.g., [6][1] has second count [1])
        additional_counts: Vec<RepeatCount>,
        /// Trailing sequence after the repeat count (VEP-specific notation)
        trailing: Option<Sequence>,
    },

    /// Multi-repeat: complex tandem repeat with multiple sequence/count pairs
    /// (e.g., GT[2]GC[2]GTGCATGAGTGTGCG[1])
    /// Used when a repeat region contains multiple different repeat units
    MultiRepeat {
        /// The sequence of repeat units (sequence + count pairs)
        units: Vec<RepeatUnit>,
    },

    /// Identity: no change (e.g., =)
    /// When whole_entity is true, this represents "c.=" (no position), similar to p.=
    Identity {
        sequence: Option<Sequence>,
        /// True if this represents whole transcript/CDS identity (c.=), false for position-specific (c.100=)
        whole_entity: bool,
    },

    /// Conversion: gene conversion where sequence is copied from another location
    /// (e.g., conNM_000001.1:c.100_200)
    Conversion {
        /// The source of the converted sequence (e.g., "NM_000001.1:c.100_200")
        source: String,
    },

    /// Unknown effect: ? (e.g., c.?)
    /// When whole_entity is true, represents "c.?" (whole entity unknown effect)
    Unknown {
        /// True if this represents whole transcript/CDS unknown (c.?), false for position-specific
        whole_entity: bool,
    },

    /// Methylation: epigenetic modification (e.g., |gom, |lom, |met=)
    Methylation { status: MethylationStatus },

    /// Copy number: specifies the number of copies of a region (e.g., copy2, copy4)
    CopyNumber { count: u64 },

    /// Aberrant splicing (RNA-specific): r.spl
    /// Indicates that the RNA is aberrantly spliced, but the exact effect is unknown
    Splice,

    /// No RNA product (RNA-specific): r.0
    /// Indicates that no RNA is produced (similar to p.0 for protein)
    NoProduct,

    /// Position-only notation: a location without an edit type
    /// Used in ClinVar to indicate a position of interest (e.g., c.5238_5240)
    /// This is non-standard HGVS but appears in real-world data
    PositionOnly,
}

impl NaEdit {
    /// Create a whole-entity identity (c.= or n.=)
    pub fn whole_entity_identity() -> Self {
        Self::Identity {
            sequence: None,
            whole_entity: true,
        }
    }

    /// Create a position-specific identity (e.g., c.100=)
    pub fn position_identity() -> Self {
        Self::Identity {
            sequence: None,
            whole_entity: false,
        }
    }

    /// Check if this is a whole-entity identity edit
    pub fn is_whole_entity_identity(&self) -> bool {
        matches!(
            self,
            NaEdit::Identity {
                whole_entity: true,
                ..
            }
        )
    }

    /// Create a whole-entity unknown (c.? or n.?)
    pub fn whole_entity_unknown() -> Self {
        Self::Unknown { whole_entity: true }
    }

    /// Check if this is a whole-entity unknown edit
    pub fn is_whole_entity_unknown(&self) -> bool {
        matches!(self, NaEdit::Unknown { whole_entity: true })
    }

    /// Check if this is an RNA splice edit (r.spl)
    pub fn is_splice(&self) -> bool {
        matches!(self, NaEdit::Splice)
    }

    /// Check if this is an RNA no product edit (r.0)
    pub fn is_no_product(&self) -> bool {
        matches!(self, NaEdit::NoProduct)
    }

    /// Check if this is a whole-entity edit (identity, unknown, splice, or no product)
    /// These edits don't require a position to be displayed
    pub fn is_whole_entity(&self) -> bool {
        self.is_whole_entity_identity()
            || self.is_whole_entity_unknown()
            || self.is_splice()
            || self.is_no_product()
    }

    /// Format this edit with lowercase nucleotides (for RNA display).
    pub fn to_rna_string(&self) -> String {
        match self {
            NaEdit::Substitution {
                reference,
                alternative,
            } => {
                format!(
                    "{}>{}",
                    reference.to_lowercase_char(),
                    alternative.to_lowercase_char()
                )
            }
            NaEdit::SubstitutionNoRef { alternative } => {
                format!(">{}", alternative.to_lowercase_char())
            }
            NaEdit::Deletion { sequence, length } => {
                let mut s = String::from("del");
                if let Some(seq) = sequence {
                    s.push_str(&seq.to_lowercase_string());
                } else if let Some(len) = length {
                    s.push_str(&len.to_string());
                }
                s
            }
            NaEdit::Insertion { sequence } => {
                format!("ins{}", sequence.to_rna_string())
            }
            NaEdit::Delins { sequence } => {
                format!("delins{}", sequence.to_rna_string())
            }
            NaEdit::Duplication {
                sequence,
                length,
                uncertain_extent,
            } => {
                let mut s = String::from("dup");
                if let Some(seq) = sequence {
                    s.push_str(&seq.to_lowercase_string());
                } else if let Some(len) = length {
                    s.push_str(&len.to_string());
                }
                if let Some(extent) = uncertain_extent {
                    s.push_str(&extent.to_string());
                }
                s
            }
            NaEdit::DupIns { sequence } => {
                format!("dupins{}", sequence.to_rna_string())
            }
            NaEdit::Inversion { sequence, length } => {
                let mut s = String::from("inv");
                if let Some(seq) = sequence {
                    s.push_str(&seq.to_lowercase_string());
                } else if let Some(len) = length {
                    s.push_str(&len.to_string());
                }
                s
            }
            NaEdit::Repeat {
                sequence,
                count,
                additional_counts,
                trailing,
            } => {
                let mut s = String::new();
                if let Some(seq) = sequence {
                    s.push_str(&seq.to_lowercase_string());
                }
                s.push_str(&count.to_string());
                for additional in additional_counts {
                    s.push_str(&additional.to_string());
                }
                if let Some(trail) = trailing {
                    s.push_str(&trail.to_lowercase_string());
                }
                s
            }
            NaEdit::MultiRepeat { units } => units
                .iter()
                .map(|u| format!("{}{}", u.sequence.to_lowercase_string(), u.count))
                .collect(),
            NaEdit::Conversion { source } => {
                format!("con{}", source)
            }
            NaEdit::Identity { sequence, .. } => {
                let mut s = String::new();
                if let Some(seq) = sequence {
                    s.push_str(&seq.to_lowercase_string());
                }
                s.push('=');
                s
            }
            NaEdit::Unknown { .. } => String::from("?"),
            NaEdit::Methylation { status } => status.to_string(),
            NaEdit::CopyNumber { count } => format!("copy{}", count),
            NaEdit::Splice => String::from("spl"),
            NaEdit::NoProduct => String::from("0"),
            NaEdit::PositionOnly => String::new(),
        }
    }
}

impl fmt::Display for NaEdit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NaEdit::Substitution {
                reference,
                alternative,
            } => {
                write!(f, "{}>{}", reference, alternative)
            }
            NaEdit::SubstitutionNoRef { alternative } => {
                write!(f, ">{}", alternative)
            }
            NaEdit::Deletion { sequence, length } => {
                write!(f, "del")?;
                if let Some(seq) = sequence {
                    write!(f, "{}", seq)?;
                } else if let Some(len) = length {
                    write!(f, "{}", len)?;
                }
                Ok(())
            }
            NaEdit::Insertion { sequence } => {
                write!(f, "ins{}", sequence)
            }
            NaEdit::Delins { sequence } => {
                write!(f, "delins{}", sequence)
            }
            NaEdit::Duplication {
                sequence,
                length,
                uncertain_extent,
            } => {
                write!(f, "dup")?;
                if let Some(seq) = sequence {
                    write!(f, "{}", seq)?;
                } else if let Some(len) = length {
                    write!(f, "{}", len)?;
                }
                if let Some(extent) = uncertain_extent {
                    write!(f, "{}", extent)?;
                }
                Ok(())
            }
            NaEdit::DupIns { sequence } => {
                write!(f, "dupins{}", sequence)
            }
            NaEdit::Inversion { sequence, length } => {
                write!(f, "inv")?;
                if let Some(seq) = sequence {
                    write!(f, "{}", seq)?;
                } else if let Some(len) = length {
                    write!(f, "{}", len)?;
                }
                Ok(())
            }
            NaEdit::Repeat {
                sequence,
                count,
                additional_counts,
                trailing,
            } => {
                if let Some(seq) = sequence {
                    write!(f, "{}", seq)?;
                }
                write!(f, "{}", count)?;
                for additional in additional_counts {
                    write!(f, "{}", additional)?;
                }
                if let Some(trail) = trailing {
                    write!(f, "{}", trail)?;
                }
                Ok(())
            }
            NaEdit::MultiRepeat { units } => {
                for unit in units {
                    write!(f, "{}", unit)?;
                }
                Ok(())
            }
            NaEdit::Identity { sequence, .. } => {
                // HGVS standard format: G= (sequence before equals)
                if let Some(seq) = sequence {
                    write!(f, "{}", seq)?;
                }
                write!(f, "=")
            }
            NaEdit::Conversion { source } => {
                write!(f, "con{}", source)
            }
            NaEdit::Unknown { .. } => {
                write!(f, "?")
            }
            NaEdit::Methylation { status } => {
                write!(f, "{}", status)
            }
            NaEdit::CopyNumber { count } => {
                write!(f, "copy{}", count)
            }
            NaEdit::Splice => {
                write!(f, "spl")
            }
            NaEdit::NoProduct => {
                write!(f, "0")
            }
            NaEdit::PositionOnly => {
                // Position-only: no edit to display (position is handled by the location)
                Ok(())
            }
        }
    }
}

/// Amino acid sequence
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AminoAcidSeq(pub Vec<AminoAcid>);

impl AminoAcidSeq {
    pub fn new(aas: Vec<AminoAcid>) -> Self {
        Self(aas)
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl fmt::Display for AminoAcidSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for aa in &self.0 {
            write!(f, "{}", aa)?;
        }
        Ok(())
    }
}

/// Extension direction for protein extensions
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ExtDirection {
    /// N-terminal extension (upstream of Met1)
    NTerminal,
    /// C-terminal extension (downstream of Ter)
    CTerminal,
}

/// Protein edit types
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ProteinEdit {
    /// Substitution: single amino acid change (e.g., p.Val600Glu)
    Substitution {
        reference: AminoAcid,
        alternative: AminoAcid,
    },

    /// Deletion: removal of amino acids (e.g., p.Lys23del, p.Lys228_Met259del32, p.Gln367_Gly372delQRGEPG)
    Deletion {
        /// Optional amino acid sequence that was deleted (e.g., Lys in p.Lys903delLys)
        sequence: Option<AminoAcidSeq>,
        /// Explicit length of deleted region (e.g., 32 in p.Lys228_Met259del32)
        count: Option<u64>,
    },

    /// Insertion: addition of amino acids (e.g., p.Lys23_Leu24insGln)
    Insertion { sequence: AminoAcidSeq },

    /// Deletion-insertion: replacement of amino acids
    Delins { sequence: AminoAcidSeq },

    /// Duplication: copy of amino acids (e.g., p.Lys23dup)
    Duplication,

    /// Frameshift: reading frame change (e.g., p.Lys23fs, p.Lys23fsTer45, p.Arg97ProfsTer23)
    Frameshift {
        /// The new amino acid at the frameshift position (e.g., Pro in p.Arg97ProfsTer23)
        new_aa: Option<AminoAcid>,
        /// Position of new termination codon (e.g., 23 in p.Arg97ProfsTer23)
        ter_pos: Option<u64>,
    },

    /// Extension: extension beyond normal start/stop (e.g., p.Met1ext-5, p.Ter110Glnext*17)
    Extension {
        /// The new amino acid at the extension position (e.g., Gln in p.Ter110Glnext*17)
        new_aa: Option<AminoAcid>,
        direction: ExtDirection,
        count: Option<i64>,
    },

    /// Identity: no change (e.g., p.= or p.(=) for predicted)
    Identity {
        /// True if this is a predicted no-change (p.(=)), false for certain (p.=)
        predicted: bool,
        /// True if this represents whole protein identity (p.=), false for position-specific (p.Val600=)
        /// When true, the location in loc_edit should be ignored in display
        whole_protein: bool,
    },

    /// Unknown effect (e.g., p.?, p.(?), or p.Met1?)
    Unknown {
        /// True if this is a predicted unknown (p.(?)), false for certain (p.?)
        predicted: bool,
        /// True if this represents whole-protein unknown (p.?), false for position-specific (p.Met1?)
        whole_protein: bool,
    },

    /// No protein produced (e.g., p.0 or p.0? for predicted)
    /// Indicates the variant results in no protein being produced
    NoProtein {
        /// True if this is a predicted no-protein (p.0?), false for certain (p.0)
        predicted: bool,
    },

    /// Repeat: amino acid repeat (e.g., p.223PA[3], p.179_180AP[9])
    /// Used for short tandem repeat expansions in protein notation
    /// The sequence uses single-letter amino acid codes
    Repeat {
        /// The repeated amino acid sequence (single-letter codes)
        sequence: AminoAcidSeq,
        /// The repeat count
        count: RepeatCount,
    },

    /// MultiRepeat: multiple adjacent repeat units (e.g., p.255_258A[4]SAAAA[1])
    /// Used for complex repeat expansions with different repeat units
    MultiRepeat {
        /// The repeated amino acid sequences with their counts
        units: Vec<(AminoAcidSeq, RepeatCount)>,
    },
}

impl ProteinEdit {
    /// Create a whole-protein identity (p.=)
    pub fn whole_protein_identity() -> Self {
        Self::Identity {
            predicted: false,
            whole_protein: true,
        }
    }

    /// Create a predicted whole-protein identity (p.(=))
    pub fn whole_protein_identity_predicted() -> Self {
        Self::Identity {
            predicted: true,
            whole_protein: true,
        }
    }

    /// Create a position-specific identity (e.g., p.Val600=)
    pub fn position_identity() -> Self {
        Self::Identity {
            predicted: false,
            whole_protein: false,
        }
    }

    /// Check if this is a whole-protein identity edit
    pub fn is_whole_protein_identity(&self) -> bool {
        matches!(
            self,
            ProteinEdit::Identity {
                whole_protein: true,
                ..
            }
        )
    }

    /// Create a no-protein edit (p.0)
    pub fn no_protein() -> Self {
        Self::NoProtein { predicted: false }
    }

    /// Create a predicted no-protein edit (p.0?)
    pub fn no_protein_predicted() -> Self {
        Self::NoProtein { predicted: true }
    }

    /// Check if this is a no-protein edit
    pub fn is_no_protein(&self) -> bool {
        matches!(self, ProteinEdit::NoProtein { .. })
    }

    /// Create a whole-protein unknown (p.?)
    pub fn whole_protein_unknown() -> Self {
        Self::Unknown {
            predicted: false,
            whole_protein: true,
        }
    }

    /// Create a predicted whole-protein unknown (p.(?))
    pub fn whole_protein_unknown_predicted() -> Self {
        Self::Unknown {
            predicted: true,
            whole_protein: true,
        }
    }

    /// Create a position-specific unknown (e.g., p.Met1?)
    pub fn position_unknown() -> Self {
        Self::Unknown {
            predicted: false,
            whole_protein: false,
        }
    }

    /// Check if this is a whole-protein unknown edit
    pub fn is_whole_protein_unknown(&self) -> bool {
        matches!(
            self,
            ProteinEdit::Unknown {
                whole_protein: true,
                ..
            }
        )
    }
}

impl fmt::Display for ProteinEdit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ProteinEdit::Substitution {
                reference: _,
                alternative,
            } => {
                write!(f, "{}", alternative)
            }
            ProteinEdit::Deletion { sequence, count } => {
                write!(f, "del")?;
                if let Some(seq) = sequence {
                    write!(f, "{}", seq)?;
                } else if let Some(n) = count {
                    write!(f, "{}", n)?;
                }
                Ok(())
            }
            ProteinEdit::Insertion { sequence } => write!(f, "ins{}", sequence),
            ProteinEdit::Delins { sequence } => write!(f, "delins{}", sequence),
            ProteinEdit::Duplication => write!(f, "dup"),
            ProteinEdit::Frameshift { new_aa, ter_pos } => {
                if let Some(aa) = new_aa {
                    write!(f, "{}", aa)?;
                }
                write!(f, "fs")?;
                if let Some(pos) = ter_pos {
                    write!(f, "Ter{}", pos)?;
                }
                Ok(())
            }
            ProteinEdit::Extension {
                new_aa,
                direction,
                count,
            } => {
                if let Some(aa) = new_aa {
                    write!(f, "{}", aa)?;
                }
                write!(f, "ext")?;
                match (direction, count) {
                    (ExtDirection::NTerminal, Some(n)) => write!(f, "{}", n)?,
                    (ExtDirection::CTerminal, Some(n)) => write!(f, "*{}", n)?,
                    (ExtDirection::CTerminal, None) => write!(f, "*?")?,
                    (ExtDirection::NTerminal, None) => {} // No suffix for unknown N-terminal
                }
                Ok(())
            }
            ProteinEdit::Identity {
                predicted: false, ..
            } => write!(f, "="),
            ProteinEdit::Identity {
                predicted: true, ..
            } => write!(f, "(=)"),
            ProteinEdit::Unknown {
                predicted: false, ..
            } => write!(f, "?"),
            ProteinEdit::Unknown {
                predicted: true, ..
            } => write!(f, "(?)"),
            ProteinEdit::NoProtein { predicted: false } => write!(f, "0"),
            ProteinEdit::NoProtein { predicted: true } => write!(f, "0?"),
            ProteinEdit::Repeat { sequence, count } => {
                write!(f, "{}{}", sequence, count)
            }
            ProteinEdit::MultiRepeat { units } => {
                for (sequence, count) in units {
                    write!(f, "{}{}", sequence, count)?;
                }
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_substitution_display() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        assert_eq!(format!("{}", edit), "A>G");
    }

    #[test]
    fn test_deletion_display() {
        let edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        assert_eq!(format!("{}", edit), "del");

        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        assert_eq!(format!("{}", edit), "delATG");

        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(101),
        };
        assert_eq!(format!("{}", edit), "del101");
    }

    #[test]
    fn test_insertion_display() {
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
        };
        assert_eq!(format!("{}", edit), "insATG");
    }

    #[test]
    fn test_delins_display() {
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
        };
        assert_eq!(format!("{}", edit), "delinsATG");
    }

    #[test]
    fn test_duplication_display() {
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        };
        assert_eq!(format!("{}", edit), "dup");

        let edit = NaEdit::Duplication {
            sequence: None,
            length: Some(101),
            uncertain_extent: None,
        };
        assert_eq!(format!("{}", edit), "dup101");
    }

    #[test]
    fn test_repeat_display() {
        let edit = NaEdit::Repeat {
            sequence: Some(Sequence::from_str("CAG").unwrap()),
            count: RepeatCount::Exact(12),
            additional_counts: vec![],
            trailing: None,
        };
        assert_eq!(format!("{}", edit), "CAG[12]");

        // Test with trailing sequence (VEP notation)
        let edit_with_trailing = NaEdit::Repeat {
            sequence: None,
            count: RepeatCount::Exact(3),
            additional_counts: vec![],
            trailing: Some(Sequence::from_str("T").unwrap()),
        };
        assert_eq!(format!("{}", edit_with_trailing), "[3]T");

        // Test with genotype notation (nested repeats)
        let edit_genotype = NaEdit::Repeat {
            sequence: Some(Sequence::from_str("A").unwrap()),
            count: RepeatCount::Exact(6),
            additional_counts: vec![RepeatCount::Exact(1)],
            trailing: None,
        };
        assert_eq!(format!("{}", edit_genotype), "A[6][1]");
    }

    #[test]
    fn test_sequence_from_str() {
        let seq = Sequence::from_str("ATGC").unwrap();
        assert_eq!(seq.len(), 4);
        assert_eq!(format!("{}", seq), "ATGC");
    }

    #[test]
    fn test_protein_substitution() {
        let edit = ProteinEdit::Substitution {
            reference: AminoAcid::Val,
            alternative: AminoAcid::Glu,
        };
        assert_eq!(format!("{}", edit), "Glu");
    }

    #[test]
    fn test_protein_frameshift() {
        let edit = ProteinEdit::Frameshift {
            new_aa: None,
            ter_pos: Some(45),
        };
        assert_eq!(format!("{}", edit), "fsTer45");

        let edit = ProteinEdit::Frameshift {
            new_aa: None,
            ter_pos: None,
        };
        assert_eq!(format!("{}", edit), "fs");

        // Frameshift with new amino acid
        let edit = ProteinEdit::Frameshift {
            new_aa: Some(crate::hgvs::location::AminoAcid::Pro),
            ter_pos: Some(23),
        };
        assert_eq!(format!("{}", edit), "ProfsTer23");
    }

    // ============== Base Tests ==============

    #[test]
    fn test_base_from_char_standard() {
        assert_eq!(Base::from_char('A'), Some(Base::A));
        assert_eq!(Base::from_char('C'), Some(Base::C));
        assert_eq!(Base::from_char('G'), Some(Base::G));
        assert_eq!(Base::from_char('T'), Some(Base::T));
        assert_eq!(Base::from_char('U'), Some(Base::U));
    }

    #[test]
    fn test_base_from_char_lowercase() {
        assert_eq!(Base::from_char('a'), Some(Base::A));
        assert_eq!(Base::from_char('c'), Some(Base::C));
        assert_eq!(Base::from_char('g'), Some(Base::G));
        assert_eq!(Base::from_char('t'), Some(Base::T));
        assert_eq!(Base::from_char('u'), Some(Base::U));
    }

    #[test]
    fn test_base_from_char_iupac() {
        assert_eq!(Base::from_char('R'), Some(Base::R));
        assert_eq!(Base::from_char('Y'), Some(Base::Y));
        assert_eq!(Base::from_char('S'), Some(Base::S));
        assert_eq!(Base::from_char('W'), Some(Base::W));
        assert_eq!(Base::from_char('K'), Some(Base::K));
        assert_eq!(Base::from_char('M'), Some(Base::M));
        assert_eq!(Base::from_char('B'), Some(Base::B));
        assert_eq!(Base::from_char('D'), Some(Base::D));
        assert_eq!(Base::from_char('H'), Some(Base::H));
        assert_eq!(Base::from_char('V'), Some(Base::V));
        assert_eq!(Base::from_char('N'), Some(Base::N));
    }

    #[test]
    fn test_base_from_char_invalid() {
        assert_eq!(Base::from_char('X'), None);
        assert_eq!(Base::from_char('Z'), None);
        assert_eq!(Base::from_char('1'), None);
        assert_eq!(Base::from_char(' '), None);
    }

    #[test]
    fn test_base_to_char() {
        assert_eq!(Base::A.to_char(), 'A');
        assert_eq!(Base::C.to_char(), 'C');
        assert_eq!(Base::G.to_char(), 'G');
        assert_eq!(Base::T.to_char(), 'T');
        assert_eq!(Base::U.to_char(), 'U');
        assert_eq!(Base::N.to_char(), 'N');
    }

    #[test]
    fn test_base_display() {
        assert_eq!(format!("{}", Base::A), "A");
        assert_eq!(format!("{}", Base::N), "N");
    }

    #[test]
    fn test_base_roundtrip() {
        for base in [
            Base::A,
            Base::C,
            Base::G,
            Base::T,
            Base::U,
            Base::R,
            Base::Y,
            Base::N,
        ] {
            let c = base.to_char();
            assert_eq!(Base::from_char(c), Some(base));
        }
    }

    // ============== Sequence Tests ==============

    #[test]
    fn test_sequence_new() {
        let seq = Sequence::new(vec![Base::A, Base::T, Base::G]);
        assert_eq!(seq.len(), 3);
    }

    #[test]
    fn test_sequence_len() {
        let seq = Sequence::from_str("ATGCATGC").unwrap();
        assert_eq!(seq.len(), 8);
    }

    #[test]
    fn test_sequence_is_empty() {
        let empty = Sequence::new(vec![]);
        assert!(empty.is_empty());

        let nonempty = Sequence::from_str("A").unwrap();
        assert!(!nonempty.is_empty());
    }

    #[test]
    fn test_sequence_bases() {
        let seq = Sequence::from_str("ATG").unwrap();
        let bases = seq.bases();
        assert_eq!(bases.len(), 3);
        assert_eq!(bases[0], Base::A);
        assert_eq!(bases[1], Base::T);
        assert_eq!(bases[2], Base::G);
    }

    #[test]
    fn test_sequence_from_str_invalid() {
        let result = Sequence::from_str("ATXG");
        assert!(result.is_err());
    }

    #[test]
    fn test_sequence_display() {
        let seq = Sequence::from_str("ATGN").unwrap();
        assert_eq!(format!("{}", seq), "ATGN");
    }

    // ============== RepeatCount Tests ==============

    #[test]
    fn test_repeat_count_exact_display() {
        assert_eq!(format!("{}", RepeatCount::Exact(12)), "[12]");
        assert_eq!(format!("{}", RepeatCount::Exact(1)), "[1]");
        assert_eq!(format!("{}", RepeatCount::Exact(100)), "[100]");
    }

    #[test]
    fn test_repeat_count_range_display() {
        assert_eq!(format!("{}", RepeatCount::Range(10, 15)), "[10_15]");
    }

    #[test]
    fn test_repeat_count_min_uncertain_display() {
        assert_eq!(format!("{}", RepeatCount::MinUncertain(10)), "[10_?]");
    }

    #[test]
    fn test_repeat_count_max_uncertain_display() {
        assert_eq!(format!("{}", RepeatCount::MaxUncertain(20)), "[?_20]");
    }

    #[test]
    fn test_repeat_count_unknown_display() {
        assert_eq!(format!("{}", RepeatCount::Unknown), "[?]");
    }

    // ============== InsertedSequence Tests ==============

    #[test]
    fn test_inserted_sequence_literal() {
        let seq = Sequence::from_str("ATG").unwrap();
        let ins = InsertedSequence::literal(seq);
        assert!(ins.is_literal());
    }

    #[test]
    fn test_inserted_sequence_is_literal() {
        let ins = InsertedSequence::Count(10);
        assert!(!ins.is_literal());

        let ins = InsertedSequence::Literal(Sequence::from_str("A").unwrap());
        assert!(ins.is_literal());
    }

    #[test]
    fn test_inserted_sequence_as_literal() {
        let seq = Sequence::from_str("ATG").unwrap();
        let ins = InsertedSequence::Literal(seq);
        assert!(ins.as_literal().is_some());
        assert_eq!(ins.as_literal().unwrap().len(), 3);

        let ins = InsertedSequence::Count(5);
        assert!(ins.as_literal().is_none());
    }

    #[test]
    fn test_inserted_sequence_bases() {
        let ins = InsertedSequence::Literal(Sequence::from_str("ATGC").unwrap());
        assert_eq!(ins.bases().unwrap().len(), 4);

        let ins = InsertedSequence::Count(10);
        assert!(ins.bases().is_none());
    }

    #[test]
    fn test_inserted_sequence_len() {
        assert_eq!(
            InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()).len(),
            Some(3)
        );
        assert_eq!(InsertedSequence::Count(10).len(), Some(10));
        assert_eq!(
            InsertedSequence::Repeat {
                base: Base::A,
                count: RepeatCount::Exact(5)
            }
            .len(),
            Some(5)
        );
        assert_eq!(InsertedSequence::Range(5, 10).len(), None);
        assert_eq!(InsertedSequence::Empty.len(), Some(0));
    }

    #[test]
    fn test_inserted_sequence_is_empty() {
        assert!(InsertedSequence::Empty.is_empty());
        assert!(InsertedSequence::Literal(Sequence::new(vec![])).is_empty());
        assert!(!InsertedSequence::Literal(Sequence::from_str("A").unwrap()).is_empty());
        assert!(!InsertedSequence::Count(5).is_empty());
    }

    #[test]
    fn test_inserted_sequence_count_display() {
        assert_eq!(format!("{}", InsertedSequence::Count(10)), "10");
    }

    #[test]
    fn test_inserted_sequence_range_display() {
        assert_eq!(format!("{}", InsertedSequence::Range(10, 20)), "(10_20)");
    }

    #[test]
    fn test_inserted_sequence_repeat_display() {
        assert_eq!(
            format!(
                "{}",
                InsertedSequence::Repeat {
                    base: Base::A,
                    count: RepeatCount::Exact(10)
                }
            ),
            "A[10]"
        );
    }

    #[test]
    fn test_inserted_sequence_complex_display() {
        let parts = vec![
            InsertedPart::Literal(Sequence::from_str("ATG").unwrap()),
            InsertedPart::Repeat {
                base: Base::N,
                count: RepeatCount::Exact(5),
            },
        ];
        let ins = InsertedSequence::Complex(parts);
        assert_eq!(format!("{}", ins), "[ATG;N[5]]");
    }

    #[test]
    fn test_inserted_sequence_empty_display() {
        assert_eq!(format!("{}", InsertedSequence::Empty), "");
    }

    // ============== InsertedPart Tests ==============

    #[test]
    fn test_inserted_part_literal_display() {
        let part = InsertedPart::Literal(Sequence::from_str("ATG").unwrap());
        assert_eq!(format!("{}", part), "ATG");
    }

    #[test]
    fn test_inserted_part_repeat_display() {
        let part = InsertedPart::Repeat {
            base: Base::A,
            count: RepeatCount::Exact(10),
        };
        assert_eq!(format!("{}", part), "A[10]");
    }

    #[test]
    fn test_inserted_part_position_range_display() {
        let part = InsertedPart::PositionRange {
            start: 401,
            end: 419,
        };
        assert_eq!(format!("{}", part), "401_419");
    }

    // ============== MethylationStatus Tests ==============

    #[test]
    fn test_methylation_status_display() {
        assert_eq!(format!("{}", MethylationStatus::GainOfMethylation), "|gom");
        assert_eq!(format!("{}", MethylationStatus::LossOfMethylation), "|lom");
        assert_eq!(format!("{}", MethylationStatus::Unchanged), "|met=");
    }

    // ============== UncertainDupExtent Tests ==============

    #[test]
    fn test_uncertain_dup_extent_unknown_display() {
        assert_eq!(format!("{}", UncertainDupExtent::Unknown), "?");
    }

    #[test]
    fn test_uncertain_dup_extent_range_display() {
        assert_eq!(
            format!("{}", UncertainDupExtent::Range(731, 741)),
            "(731_741)"
        );
    }

    // ============== Additional NaEdit Tests ==============

    #[test]
    fn test_na_edit_identity_display() {
        let edit = NaEdit::Identity {
            sequence: None,
            whole_entity: false,
        };
        assert_eq!(format!("{}", edit), "=");

        let edit = NaEdit::Identity {
            sequence: Some(Sequence::from_str("A").unwrap()),
            whole_entity: false,
        };
        assert_eq!(format!("{}", edit), "A=");
    }

    #[test]
    fn test_na_edit_whole_entity_identity() {
        let edit = NaEdit::whole_entity_identity();
        assert_eq!(format!("{}", edit), "=");
    }

    #[test]
    fn test_na_edit_unknown_display() {
        let edit = NaEdit::Unknown {
            whole_entity: false,
        };
        assert_eq!(format!("{}", edit), "?");

        let edit = NaEdit::Unknown { whole_entity: true };
        assert_eq!(format!("{}", edit), "?");
    }

    #[test]
    fn test_na_edit_conversion_display() {
        let edit = NaEdit::Conversion {
            source: "NM_000001.1:c.100_200".to_string(),
        };
        assert_eq!(format!("{}", edit), "conNM_000001.1:c.100_200");
    }

    #[test]
    fn test_na_edit_methylation_display() {
        let edit = NaEdit::Methylation {
            status: MethylationStatus::GainOfMethylation,
        };
        assert_eq!(format!("{}", edit), "|gom");
    }

    #[test]
    fn test_na_edit_copy_number_display() {
        let edit = NaEdit::CopyNumber { count: 2 };
        assert_eq!(format!("{}", edit), "copy2");

        let edit = NaEdit::CopyNumber { count: 4 };
        assert_eq!(format!("{}", edit), "copy4");
    }

    #[test]
    fn test_na_edit_splice_display() {
        let edit = NaEdit::Splice;
        assert_eq!(format!("{}", edit), "spl");
    }

    #[test]
    fn test_na_edit_no_product_display() {
        let edit = NaEdit::NoProduct;
        assert_eq!(format!("{}", edit), "0");
    }

    #[test]
    fn test_na_edit_inversion_display() {
        let edit = NaEdit::Inversion {
            sequence: None,
            length: None,
        };
        assert_eq!(format!("{}", edit), "inv");

        let edit = NaEdit::Inversion {
            sequence: None,
            length: Some(3),
        };
        assert_eq!(format!("{}", edit), "inv3");

        let edit = NaEdit::Inversion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        assert_eq!(format!("{}", edit), "invATG");
    }

    #[test]
    fn test_na_edit_duplication_with_extent() {
        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: Some(UncertainDupExtent::Unknown),
        };
        assert_eq!(format!("{}", edit), "dup?");

        let edit = NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: Some(UncertainDupExtent::Range(731, 741)),
        };
        assert_eq!(format!("{}", edit), "dup(731_741)");
    }

    #[test]
    fn test_na_edit_insertion_count() {
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Count(10),
        };
        assert_eq!(format!("{}", edit), "ins10");
    }

    #[test]
    fn test_na_edit_insertion_range() {
        let edit = NaEdit::Insertion {
            sequence: InsertedSequence::Range(5, 10),
        };
        assert_eq!(format!("{}", edit), "ins(5_10)");
    }

    #[test]
    fn test_na_edit_delins_with_repeat() {
        let edit = NaEdit::Delins {
            sequence: InsertedSequence::Repeat {
                base: Base::A,
                count: RepeatCount::Exact(10),
            },
        };
        assert_eq!(format!("{}", edit), "delinsA[10]");
    }

    // ============== ProteinEdit Tests ==============

    #[test]
    fn test_protein_edit_deletion() {
        // Plain deletion
        let edit = ProteinEdit::Deletion {
            sequence: None,
            count: None,
        };
        assert_eq!(format!("{}", edit), "del");

        // Deletion with count
        let edit = ProteinEdit::Deletion {
            sequence: None,
            count: Some(32),
        };
        assert_eq!(format!("{}", edit), "del32");

        // Deletion with sequence
        let edit = ProteinEdit::Deletion {
            sequence: Some(AminoAcidSeq(vec![AminoAcid::Lys])),
            count: None,
        };
        assert_eq!(format!("{}", edit), "delLys");
    }

    #[test]
    fn test_protein_edit_duplication() {
        let edit = ProteinEdit::Duplication;
        assert_eq!(format!("{}", edit), "dup");
    }

    #[test]
    fn test_protein_edit_identity() {
        let edit = ProteinEdit::Identity {
            predicted: false,
            whole_protein: true,
        };
        assert_eq!(format!("{}", edit), "=");

        let edit = ProteinEdit::Identity {
            predicted: true,
            whole_protein: true,
        };
        assert_eq!(format!("{}", edit), "(=)");
    }

    #[test]
    fn test_protein_edit_unknown() {
        let edit = ProteinEdit::Unknown {
            predicted: false,
            whole_protein: true,
        };
        assert_eq!(format!("{}", edit), "?");

        let edit = ProteinEdit::Unknown {
            predicted: true,
            whole_protein: true,
        };
        assert_eq!(format!("{}", edit), "(?)");
    }

    #[test]
    fn test_protein_edit_no_protein() {
        let edit = ProteinEdit::NoProtein { predicted: false };
        assert_eq!(format!("{}", edit), "0");

        let edit = ProteinEdit::NoProtein { predicted: true };
        assert_eq!(format!("{}", edit), "0?");
    }

    #[test]
    fn test_protein_edit_insertion() {
        let edit = ProteinEdit::Insertion {
            sequence: AminoAcidSeq(vec![AminoAcid::Ala, AminoAcid::Gly]),
        };
        assert_eq!(format!("{}", edit), "insAlaGly");
    }

    #[test]
    fn test_protein_edit_extension() {
        let edit = ProteinEdit::Extension {
            direction: ExtDirection::NTerminal,
            new_aa: Some(AminoAcid::Met),
            count: Some(-5),
        };
        // Extension format may vary based on implementation
        let formatted = format!("{}", edit);
        assert!(formatted.contains("ext"));
    }

    #[test]
    fn test_protein_edit_helper_methods() {
        let identity = ProteinEdit::whole_protein_identity();
        assert!(identity.is_whole_protein_identity());

        let predicted = ProteinEdit::whole_protein_identity_predicted();
        assert!(predicted.is_whole_protein_identity());

        let no_protein = ProteinEdit::no_protein();
        assert!(matches!(
            no_protein,
            ProteinEdit::NoProtein { predicted: false }
        ));
    }

    #[test]
    fn test_protein_edit_delins() {
        let edit = ProteinEdit::Delins {
            sequence: AminoAcidSeq(vec![AminoAcid::Val]),
        };
        assert_eq!(format!("{}", edit), "delinsVal");
    }
}
