//! SPDI (Sequence Position Deletion Insertion) variant representation.
//!
//! SPDI is a contextual variant representation format developed by NCBI that
//! provides an unambiguous way to describe sequence variants. The format is:
//!
//! ```text
//! sequence:position:deletion:insertion
//! ```
//!
//! Where:
//! - `sequence`: Reference sequence identifier (e.g., NC_000001.11)
//! - `position`: 0-based position (interbase coordinate)
//! - `deletion`: Deleted sequence (can be empty for insertions)
//! - `insertion`: Inserted sequence (can be empty for deletions)
//!
//! # Examples
//!
//! ```
//! use ferro_hgvs::spdi::SpdiVariant;
//!
//! // Parse an SPDI expression
//! let spdi: SpdiVariant = "NC_000001.11:12345:A:G".parse().unwrap();
//! assert_eq!(spdi.sequence, "NC_000001.11");
//! assert_eq!(spdi.position, 12345);
//! assert_eq!(spdi.deletion, "A");
//! assert_eq!(spdi.insertion, "G");
//!
//! // Format back to string
//! assert_eq!(spdi.to_string(), "NC_000001.11:12345:A:G");
//! ```
//!
//! # References
//!
//! - [NCBI SPDI Service](https://www.ncbi.nlm.nih.gov/variation/notation/)
//! - [SPDI Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7523648/)

pub mod convert;
mod parser;

pub use convert::{hgvs_to_spdi_simple, spdi_to_hgvs, ConversionError};
pub use parser::{parse_spdi, SpdiParseError};

use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// A variant in SPDI (Sequence Position Deletion Insertion) format.
///
/// SPDI provides an unambiguous contextual representation of sequence variants
/// using 0-based interbase coordinates.
///
/// # Coordinate System
///
/// SPDI uses 0-based interbase coordinates, meaning:
/// - Position 0 is before the first base
/// - Position 1 is between the first and second bases
/// - For a substitution at position 12345 (1-based), SPDI position is 12344
///
/// # Examples
///
/// | Variant Type | SPDI | Description |
/// |--------------|------|-------------|
/// | Substitution | `NC:12344:A:G` | A→G at position 12345 (1-based) |
/// | Deletion | `NC:99:ATG:` | Delete ATG starting at position 100 |
/// | Insertion | `NC:100::ATG` | Insert ATG after position 100 |
/// | Delins | `NC:99:ATG:TCA` | Replace ATG with TCA |
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SpdiVariant {
    /// Reference sequence identifier (e.g., "NC_000001.11").
    pub sequence: String,

    /// 0-based interbase position.
    pub position: u64,

    /// Deleted sequence (empty string for pure insertions).
    pub deletion: String,

    /// Inserted sequence (empty string for pure deletions).
    pub insertion: String,
}

impl SpdiVariant {
    /// Create a new SPDI variant.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Reference sequence identifier
    /// * `position` - 0-based interbase position
    /// * `deletion` - Deleted sequence (can be empty)
    /// * `insertion` - Inserted sequence (can be empty)
    ///
    /// # Examples
    ///
    /// ```
    /// use ferro_hgvs::spdi::SpdiVariant;
    ///
    /// let snv = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
    /// let del = SpdiVariant::new("NC_000001.11", 99, "ATG", "");
    /// let ins = SpdiVariant::new("NC_000001.11", 100, "", "ATG");
    /// ```
    pub fn new(
        sequence: impl Into<String>,
        position: u64,
        deletion: impl Into<String>,
        insertion: impl Into<String>,
    ) -> Self {
        Self {
            sequence: sequence.into(),
            position,
            deletion: deletion.into(),
            insertion: insertion.into(),
        }
    }

    /// Create a substitution variant.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Reference sequence identifier
    /// * `position` - 0-based position
    /// * `reference` - Reference base(s)
    /// * `alternate` - Alternate base(s)
    pub fn substitution(
        sequence: impl Into<String>,
        position: u64,
        reference: impl Into<String>,
        alternate: impl Into<String>,
    ) -> Self {
        Self::new(sequence, position, reference, alternate)
    }

    /// Create a deletion variant.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Reference sequence identifier
    /// * `position` - 0-based position
    /// * `deleted` - Deleted sequence
    pub fn deletion(
        sequence: impl Into<String>,
        position: u64,
        deleted: impl Into<String>,
    ) -> Self {
        Self::new(sequence, position, deleted, "")
    }

    /// Create an insertion variant.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Reference sequence identifier
    /// * `position` - 0-based position (insertion point)
    /// * `inserted` - Inserted sequence
    pub fn insertion(
        sequence: impl Into<String>,
        position: u64,
        inserted: impl Into<String>,
    ) -> Self {
        Self::new(sequence, position, "", inserted)
    }

    /// Create a deletion-insertion (delins) variant.
    ///
    /// # Arguments
    ///
    /// * `sequence` - Reference sequence identifier
    /// * `position` - 0-based position
    /// * `deleted` - Deleted sequence
    /// * `inserted` - Inserted sequence
    pub fn delins(
        sequence: impl Into<String>,
        position: u64,
        deleted: impl Into<String>,
        inserted: impl Into<String>,
    ) -> Self {
        Self::new(sequence, position, deleted, inserted)
    }

    /// Returns true if this is a substitution (SNV or MNV).
    pub fn is_substitution(&self) -> bool {
        !self.deletion.is_empty()
            && !self.insertion.is_empty()
            && self.deletion.len() == self.insertion.len()
    }

    /// Returns true if this is a pure deletion.
    pub fn is_deletion(&self) -> bool {
        !self.deletion.is_empty() && self.insertion.is_empty()
    }

    /// Returns true if this is a pure insertion.
    pub fn is_insertion(&self) -> bool {
        self.deletion.is_empty() && !self.insertion.is_empty()
    }

    /// Returns true if this is a deletion-insertion (delins).
    pub fn is_delins(&self) -> bool {
        !self.deletion.is_empty()
            && !self.insertion.is_empty()
            && self.deletion.len() != self.insertion.len()
    }

    /// Returns true if this represents an identity (no change).
    pub fn is_identity(&self) -> bool {
        self.deletion == self.insertion
    }

    /// Returns the variant type as a string.
    pub fn variant_type(&self) -> &'static str {
        if self.is_identity() {
            "identity"
        } else if self.is_substitution() {
            "substitution"
        } else if self.is_deletion() {
            "deletion"
        } else if self.is_insertion() {
            "insertion"
        } else {
            "delins"
        }
    }

    /// Convert 0-based SPDI position to 1-based HGVS position.
    ///
    /// For a single nucleotide at SPDI position N, the 1-based position is N+1.
    pub fn to_one_based_position(&self) -> u64 {
        self.position + 1
    }

    /// Create from 1-based position (convenience for HGVS conversion).
    ///
    /// Converts 1-based HGVS position to 0-based SPDI position internally.
    pub fn with_one_based_position(
        sequence: impl Into<String>,
        one_based_pos: u64,
        deletion: impl Into<String>,
        insertion: impl Into<String>,
    ) -> Self {
        // Convert 1-based HGVS position to 0-based SPDI position
        Self::new(
            sequence,
            one_based_pos.saturating_sub(1), // 1-based → 0-based
            deletion,
            insertion,
        )
    }
}

impl fmt::Display for SpdiVariant {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}",
            self.sequence, self.position, self.deletion, self.insertion
        )
    }
}

impl FromStr for SpdiVariant {
    type Err = SpdiParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_spdi(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spdi_new() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 12344);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_spdi_substitution() {
        let spdi = SpdiVariant::substitution("NC_000001.11", 12344, "A", "G");
        assert!(spdi.is_substitution());
        assert!(!spdi.is_deletion());
        assert!(!spdi.is_insertion());
        assert!(!spdi.is_delins());
        assert_eq!(spdi.variant_type(), "substitution");
    }

    #[test]
    fn test_spdi_deletion() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "ATG");
        assert!(spdi.is_deletion());
        assert!(!spdi.is_substitution());
        assert!(!spdi.is_insertion());
        assert_eq!(spdi.variant_type(), "deletion");
    }

    #[test]
    fn test_spdi_insertion() {
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert!(spdi.is_insertion());
        assert!(!spdi.is_substitution());
        assert!(!spdi.is_deletion());
        assert_eq!(spdi.variant_type(), "insertion");
    }

    #[test]
    fn test_spdi_delins() {
        let spdi = SpdiVariant::delins("NC_000001.11", 99, "ATG", "TTCC");
        assert!(spdi.is_delins());
        assert!(!spdi.is_substitution());
        assert_eq!(spdi.variant_type(), "delins");
    }

    #[test]
    fn test_spdi_identity() {
        let spdi = SpdiVariant::new("NC_000001.11", 100, "A", "A");
        assert!(spdi.is_identity());
        assert_eq!(spdi.variant_type(), "identity");
    }

    #[test]
    fn test_spdi_display() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
    }

    #[test]
    fn test_spdi_display_deletion() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "ATG");
        assert_eq!(spdi.to_string(), "NC_000001.11:99:ATG:");
    }

    #[test]
    fn test_spdi_display_insertion() {
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert_eq!(spdi.to_string(), "NC_000001.11:100::ATG");
    }

    #[test]
    fn test_spdi_from_str() {
        let spdi: SpdiVariant = "NC_000001.11:12344:A:G".parse().unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 12344);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_spdi_to_one_based() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        assert_eq!(spdi.to_one_based_position(), 12345);
    }

    #[test]
    fn test_spdi_with_one_based() {
        let spdi = SpdiVariant::with_one_based_position("NC_000001.11", 12345, "A", "G");
        assert_eq!(spdi.position, 12344);
    }

    #[test]
    fn test_spdi_serialize() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        let json = serde_json::to_string(&spdi).unwrap();
        assert!(json.contains("NC_000001.11"));
        assert!(json.contains("12344"));
    }

    #[test]
    fn test_spdi_deserialize() {
        let json = r#"{"sequence":"NC_000001.11","position":12344,"deletion":"A","insertion":"G"}"#;
        let spdi: SpdiVariant = serde_json::from_str(json).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 12344);
    }

    #[test]
    fn test_spdi_hash() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(SpdiVariant::new("NC_000001.11", 12344, "A", "G"));
        assert!(set.contains(&SpdiVariant::new("NC_000001.11", 12344, "A", "G")));
    }
}
