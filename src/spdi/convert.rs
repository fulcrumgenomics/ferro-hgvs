//! HGVS to SPDI conversion.
//!
//! This module provides conversion functions between HGVS and SPDI formats.
//!
//! # Coordinate Systems
//!
//! - HGVS uses 1-based coordinates (see [`OneBasedPos`])
//! - SPDI uses 0-based interbase coordinates (see [`ZeroBasedPos`])
//!
//! For a variant at genomic position 12345 (1-based):
//! - HGVS: `NC_000001.11:g.12345A>G`
//! - SPDI: `NC_000001.11:12344:A:G` (0-based)
//!
//! # Supported Conversions
//!
//! | HGVS | SPDI |
//! |------|------|
//! | Substitution `g.12345A>G` | `seq:12344:A:G` |
//! | Deletion `g.100_102del` | `seq:99:NNN:` (requires reference) |
//! | Insertion `g.100_101insATG` | `seq:100::ATG` |
//! | Delins `g.100_102delinsATG` | `seq:99:NNN:ATG` (requires reference) |
//! | Duplication `g.100_102dup` | `seq:102::NNN` (requires reference) |
//!
//! Note: Deletions, delins, and duplications require reference sequence data
//! to determine the deleted sequence.
//!
//! [`OneBasedPos`]: crate::coords::OneBasedPos
//! [`ZeroBasedPos`]: crate::coords::ZeroBasedPos

use super::SpdiVariant;
use crate::coords::{OneBasedPos, ZeroBasedPos};
use crate::error::FerroError;
use crate::hgvs::edit::{InsertedSequence, NaEdit, Sequence};
use crate::hgvs::interval::Interval;
use crate::hgvs::location::GenomePos;
use crate::hgvs::parser::accession::parse_accession;
use crate::hgvs::variant::{GenomeVariant, HgvsVariant, LocEdit};

/// Error type for conversion failures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConversionError {
    /// The variant type is not supported for conversion.
    UnsupportedVariantType {
        /// Description of the unsupported type.
        description: String,
    },
    /// Missing reference sequence data needed for conversion.
    MissingReferenceData {
        /// Description of what data is missing.
        description: String,
    },
    /// The edit type is not supported for conversion.
    UnsupportedEditType {
        /// Description of the unsupported edit.
        description: String,
    },
    /// Invalid position or interval.
    InvalidPosition {
        /// Description of the position error.
        description: String,
    },
    /// Invalid accession format.
    InvalidAccession {
        /// Description of the accession error.
        description: String,
    },
}

impl std::fmt::Display for ConversionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConversionError::UnsupportedVariantType { description } => {
                write!(
                    f,
                    "unsupported variant type for conversion: {}",
                    description
                )
            }
            ConversionError::MissingReferenceData { description } => {
                write!(f, "missing reference data: {}", description)
            }
            ConversionError::UnsupportedEditType { description } => {
                write!(f, "unsupported edit type for conversion: {}", description)
            }
            ConversionError::InvalidPosition { description } => {
                write!(f, "invalid position: {}", description)
            }
            ConversionError::InvalidAccession { description } => {
                write!(f, "invalid accession: {}", description)
            }
        }
    }
}

impl std::error::Error for ConversionError {}

impl From<ConversionError> for FerroError {
    fn from(err: ConversionError) -> Self {
        FerroError::ConversionError {
            msg: err.to_string(),
        }
    }
}

/// Helper to convert a Sequence to a String.
fn sequence_to_string(seq: &Sequence) -> String {
    seq.to_string()
}

/// Helper to convert an InsertedSequence to a String (for literal sequences only).
fn inserted_sequence_to_string(seq: &InsertedSequence) -> Option<String> {
    match seq {
        InsertedSequence::Literal(s) => Some(s.to_string()),
        _ => None,
    }
}

/// Helper to get start position from an interval.
fn get_start_pos(interval: &Interval<GenomePos>) -> Option<u64> {
    interval.start.inner().map(|p| p.base)
}

/// Helper to get end position from an interval.
fn get_end_pos(interval: &Interval<GenomePos>) -> Option<u64> {
    interval.end.inner().map(|p| p.base)
}

/// Convert a genomic HGVS variant to SPDI format.
///
/// This is a "simple" conversion that works for substitutions where the
/// reference and alternate sequences are explicitly stated in the HGVS.
/// For deletions, insertions, and duplications, use the version that takes
/// a reference provider.
///
/// # Arguments
///
/// * `variant` - A genomic HGVS variant (g. coordinate system)
///
/// # Returns
///
/// * `Ok(SpdiVariant)` - Successfully converted variant
/// * `Err(ConversionError)` - Conversion failed
///
/// # Examples
///
/// ```
/// use ferro_hgvs::spdi::convert::hgvs_to_spdi_simple;
/// use ferro_hgvs::parse_hgvs;
///
/// let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
/// let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
/// assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
/// ```
pub fn hgvs_to_spdi_simple(variant: &HgvsVariant) -> Result<SpdiVariant, ConversionError> {
    match variant {
        HgvsVariant::Genome(g) => genome_to_spdi_simple(g),
        _ => Err(ConversionError::UnsupportedVariantType {
            description: format!(
                "only genomic (g.) variants can be directly converted to SPDI, got {}",
                variant.variant_type()
            ),
        }),
    }
}

/// Convert a genomic variant to SPDI (simple conversion).
fn genome_to_spdi_simple(variant: &GenomeVariant) -> Result<SpdiVariant, ConversionError> {
    let sequence = variant.accession.to_string();
    let interval = &variant.loc_edit.location;

    // Get the edit (unwrap from Mu)
    let edit = variant
        .loc_edit
        .edit
        .inner()
        .ok_or_else(|| ConversionError::InvalidPosition {
            description: "cannot convert variant with unknown edit".to_string(),
        })?;

    // Get start position
    let start_pos = get_start_pos(interval).ok_or_else(|| ConversionError::InvalidPosition {
        description: "cannot convert variant with unknown start position".to_string(),
    })?;

    // Validate and wrap in type-safe 1-based position
    let hgvs_pos_ob =
        OneBasedPos::try_new(start_pos).ok_or_else(|| ConversionError::InvalidPosition {
            description: "position 0 is not valid in HGVS".to_string(),
        })?;

    // Convert 1-based HGVS position to 0-based SPDI position using type-safe conversion
    let spdi_pos_zb: ZeroBasedPos = hgvs_pos_ob.to_zero_based();
    let spdi_pos = spdi_pos_zb.value();

    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => {
            // Substitution: g.12345A>G -> seq:12344:A:G
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                reference.to_string(),
                alternative.to_string(),
            ))
        }
        NaEdit::Insertion { sequence: inserted } => {
            // Insertion: g.100_101insATG -> seq:100::ATG
            // For insertion, HGVS uses positions flanking the insertion point
            // SPDI position is the 0-based position where insertion happens
            let ins_str = inserted_sequence_to_string(inserted).ok_or_else(|| {
                ConversionError::MissingReferenceData {
                    description: "insertion sequence is not a literal sequence".to_string(),
                }
            })?;
            Ok(SpdiVariant::new(sequence, spdi_pos, "", ins_str))
        }
        NaEdit::Duplication {
            sequence: dup_seq, ..
        } => {
            // Duplication with sequence: g.100_102dupATG -> seq:102::ATG
            // (insertion of the duplicated sequence after the original)
            if let Some(seq) = dup_seq {
                // Position for dup is at the end of the duplicated region
                let end_pos = get_end_pos(interval).unwrap_or(start_pos);
                // Convert 1-based HGVS end position to 0-based SPDI position
                let end_pos_ob = OneBasedPos::new(end_pos);
                let spdi_end_zb = end_pos_ob.to_zero_based();
                Ok(SpdiVariant::new(
                    sequence,
                    spdi_end_zb.value(),
                    "",
                    sequence_to_string(seq),
                ))
            } else {
                Err(ConversionError::MissingReferenceData {
                    description:
                        "duplication sequence not provided; reference data needed to determine duplicated bases"
                            .to_string(),
                })
            }
        }
        NaEdit::Deletion {
            sequence: del_seq, ..
        } => {
            // Deletion: g.100_102del -> seq:99:ATG:
            // Need the deleted sequence
            if let Some(seq) = del_seq {
                Ok(SpdiVariant::new(
                    sequence,
                    spdi_pos,
                    sequence_to_string(seq),
                    "",
                ))
            } else {
                Err(ConversionError::MissingReferenceData {
                    description:
                        "deleted sequence not provided; reference data needed to determine deleted bases"
                            .to_string(),
                })
            }
        }
        NaEdit::Delins { sequence: _ins_seq } => {
            // Delins: g.100_102delinsATG -> seq:99:XYZ:ATG
            // We need to know the deleted sequence to create a valid SPDI
            // Without reference data, we cannot determine what was deleted
            // Calculate deletion length from interval (saturating to avoid overflow)
            let end_pos = get_end_pos(interval).unwrap_or(start_pos);
            let del_len = (end_pos - start_pos).saturating_add(1) as usize;

            // Without reference data, we can't know the deleted sequence
            Err(ConversionError::MissingReferenceData {
                description: format!(
                    "Cannot convert delins to SPDI: deleted sequence of length {} is unknown (no reference data)",
                    del_len
                ),
            })
        }
        NaEdit::Identity {
            sequence: id_seq, ..
        } => {
            // Identity: g.100= or g.100A=
            let ref_base = id_seq.as_ref().map(sequence_to_string).unwrap_or_default();
            Ok(SpdiVariant::new(
                sequence,
                spdi_pos,
                ref_base.clone(),
                ref_base,
            ))
        }
        NaEdit::Inversion { .. } => Err(ConversionError::UnsupportedEditType {
            description: "inversions cannot be represented in SPDI format".to_string(),
        }),
        NaEdit::Repeat { .. } => Err(ConversionError::UnsupportedEditType {
            description: "repeat variants cannot be directly converted to SPDI".to_string(),
        }),
        NaEdit::CopyNumber { .. } => Err(ConversionError::UnsupportedEditType {
            description: "copy number variants cannot be represented in SPDI format".to_string(),
        }),
        NaEdit::Conversion { .. } => Err(ConversionError::UnsupportedEditType {
            description: "conversion variants cannot be represented in SPDI format".to_string(),
        }),
        _ => Err(ConversionError::UnsupportedEditType {
            description: format!("unsupported edit type: {:?}", edit),
        }),
    }
}

/// Convert an SPDI variant to HGVS genomic format.
///
/// # Arguments
///
/// * `spdi` - The SPDI variant to convert
///
/// # Returns
///
/// * `Ok(HgvsVariant)` - Successfully converted variant
/// * `Err(ConversionError)` - Conversion failed
///
/// # Examples
///
/// ```
/// use ferro_hgvs::spdi::{SpdiVariant, convert::spdi_to_hgvs};
///
/// let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
/// let hgvs = spdi_to_hgvs(&spdi).unwrap();
/// assert_eq!(hgvs.to_string(), "NC_000001.11:g.12345A>G");
/// ```
pub fn spdi_to_hgvs(spdi: &SpdiVariant) -> Result<HgvsVariant, ConversionError> {
    // Parse the accession using the HGVS parser
    let accession = parse_accession(&spdi.sequence)
        .map(|(_, acc)| acc)
        .map_err(|_| ConversionError::InvalidAccession {
            description: format!("could not parse accession: {}", spdi.sequence),
        })?;

    // Convert 0-based SPDI position to 1-based HGVS position using type-safe conversion
    let spdi_pos_zb = ZeroBasedPos::new(spdi.position);
    let hgvs_pos_ob = spdi_pos_zb.to_one_based();
    let hgvs_pos = hgvs_pos_ob.value();

    // Determine the edit type based on deletion and insertion
    let (interval, edit) = if spdi.is_identity() {
        // Identity
        let seq = if spdi.deletion.is_empty() {
            None
        } else {
            Some(string_to_sequence(&spdi.deletion)?)
        };
        (
            Interval::point(GenomePos::new(hgvs_pos)),
            NaEdit::Identity {
                sequence: seq,
                whole_entity: false,
            },
        )
    } else if spdi.deletion.len() == 1 && spdi.insertion.len() == 1 {
        // SNV substitution
        let ref_base = char_to_base(spdi.deletion.chars().next().unwrap())?;
        let alt_base = char_to_base(spdi.insertion.chars().next().unwrap())?;
        (
            Interval::point(GenomePos::new(hgvs_pos)),
            NaEdit::Substitution {
                reference: ref_base,
                alternative: alt_base,
            },
        )
    } else if spdi.is_deletion() {
        // Pure deletion
        let del_len = spdi.deletion.len();
        let del_seq = string_to_sequence(&spdi.deletion)?;
        let interval = if del_len > 1 {
            Interval::new(
                GenomePos::new(hgvs_pos),
                GenomePos::new(hgvs_pos + del_len as u64 - 1),
            )
        } else {
            Interval::point(GenomePos::new(hgvs_pos))
        };
        (
            interval,
            NaEdit::Deletion {
                sequence: Some(del_seq),
                length: None,
            },
        )
    } else if spdi.is_insertion() {
        // Pure insertion
        // SPDI position is where insertion happens (interbase)
        // HGVS: g.pos_pos+1insX means insert between pos and pos+1
        let ins_seq = string_to_sequence(&spdi.insertion)?;
        (
            Interval::new(GenomePos::new(hgvs_pos), GenomePos::new(hgvs_pos + 1)),
            NaEdit::Insertion {
                sequence: InsertedSequence::Literal(ins_seq),
            },
        )
    } else {
        // Delins (different lengths or MNV)
        let del_len = spdi.deletion.len();
        let ins_seq = string_to_sequence(&spdi.insertion)?;
        let interval = if del_len > 1 {
            Interval::new(
                GenomePos::new(hgvs_pos),
                GenomePos::new(hgvs_pos + del_len as u64 - 1),
            )
        } else {
            Interval::point(GenomePos::new(hgvs_pos))
        };
        (
            interval,
            NaEdit::Delins {
                sequence: InsertedSequence::Literal(ins_seq),
            },
        )
    };

    Ok(HgvsVariant::Genome(GenomeVariant {
        accession,
        gene_symbol: None,
        loc_edit: LocEdit::new(interval, edit),
    }))
}

/// Helper to convert a string to a Sequence.
fn string_to_sequence(s: &str) -> Result<Sequence, ConversionError> {
    s.parse().map_err(|_| ConversionError::InvalidPosition {
        description: format!("invalid sequence: {}", s),
    })
}

/// Helper to convert a char to a Base.
fn char_to_base(c: char) -> Result<crate::hgvs::edit::Base, ConversionError> {
    crate::hgvs::edit::Base::from_char(c).ok_or_else(|| ConversionError::InvalidPosition {
        description: format!("invalid base character: {}", c),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::parser::parse_hgvs;

    // HGVS to SPDI tests

    #[test]
    fn test_hgvs_to_spdi_substitution() {
        let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 12344);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "G");
    }

    #[test]
    fn test_hgvs_to_spdi_insertion() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_101insATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 99); // 0-based, position before insertion
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_deletion_with_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_hgvs_to_spdi_deletion_without_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_delins_without_ref() {
        // Without reference data, delins with unknown deletion sequence returns error
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("unknown"));
    }

    #[test]
    fn test_hgvs_to_spdi_duplication_with_seq() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dupATG").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        // Dup becomes insertion after the duplicated region
        assert_eq!(spdi.position, 101); // end position, 0-based
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_identity() {
        let hgvs = parse_hgvs("NC_000001.11:g.100A=").unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "A");
        assert_eq!(spdi.insertion, "A");
    }

    #[test]
    fn test_hgvs_to_spdi_unsupported_coding() {
        let hgvs = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedVariantType { .. })
        ));
    }

    #[test]
    fn test_hgvs_to_spdi_unsupported_inversion() {
        let hgvs = parse_hgvs("NC_000001.11:g.100_200inv").unwrap();
        let result = hgvs_to_spdi_simple(&hgvs);
        assert!(matches!(
            result,
            Err(ConversionError::UnsupportedEditType { .. })
        ));
    }

    // SPDI to HGVS tests

    #[test]
    fn test_spdi_to_hgvs_substitution() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.12345A>G");
    }

    #[test]
    fn test_spdi_to_hgvs_deletion() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn test_spdi_to_hgvs_insertion() {
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.101_102insATG");
    }

    #[test]
    fn test_spdi_to_hgvs_delins() {
        let spdi = SpdiVariant::delins("NC_000001.11", 99, "ATG", "TTCC");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100_102delinsTTCC");
    }

    #[test]
    fn test_spdi_to_hgvs_identity() {
        let spdi = SpdiVariant::new("NC_000001.11", 99, "A", "A");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100A=");
    }

    #[test]
    fn test_spdi_to_hgvs_single_del() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 99, "A");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(hgvs.to_string(), "NC_000001.11:g.100delA");
    }

    // Roundtrip tests

    #[test]
    fn test_roundtrip_substitution() {
        let original = "NC_000001.11:g.12345A>G";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), original);
    }

    #[test]
    fn test_roundtrip_insertion() {
        let original = "NC_000001.11:g.100_101insATG";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), original);
    }

    #[test]
    fn test_roundtrip_deletion_with_seq() {
        let original = "NC_000001.11:g.100_102delATG";
        let hgvs = parse_hgvs(original).unwrap();
        let spdi = hgvs_to_spdi_simple(&hgvs).unwrap();
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), original);
    }

    #[test]
    fn test_error_display() {
        let err = ConversionError::UnsupportedVariantType {
            description: "test".to_string(),
        };
        assert!(err.to_string().contains("unsupported variant type"));

        let err = ConversionError::MissingReferenceData {
            description: "test".to_string(),
        };
        assert!(err.to_string().contains("missing reference data"));
    }

    // =========================================================================
    // P1: Reference-dependent conversion tests
    // =========================================================================

    /// Mock genomic reference for testing reference-dependent conversions.
    struct MockGenomicRef {
        /// Map of (accession, start, end) -> sequence
        sequences: std::collections::HashMap<(String, u64, u64), String>,
    }

    impl MockGenomicRef {
        fn new() -> Self {
            let mut sequences = std::collections::HashMap::new();

            // Add test sequences for NC_000001.11
            // Positions are 0-based for internal storage
            // Sequence at g.100_102 (1-based) = indices 99-101 (0-based) = "ATG"
            sequences.insert(("NC_000001.11".to_string(), 99, 102), "ATG".to_string());
            // Sequence at g.200_205 (1-based) = indices 199-204 (0-based) = "GATTACA"
            sequences.insert(
                ("NC_000001.11".to_string(), 199, 206),
                "GATTACA".to_string(),
            );
            // Sequence at g.1000_1009 (1-based) = "AAACCCGGGT"
            sequences.insert(
                ("NC_000001.11".to_string(), 999, 1009),
                "AAACCCGGGT".to_string(),
            );

            Self { sequences }
        }

        fn get_sequence(&self, accession: &str, start: u64, end: u64) -> Option<String> {
            self.sequences
                .get(&(accession.to_string(), start, end))
                .cloned()
        }
    }

    /// Convert deletion HGVS to SPDI using reference data.
    fn hgvs_to_spdi_with_ref(
        variant: &HgvsVariant,
        reference: &MockGenomicRef,
    ) -> Result<SpdiVariant, ConversionError> {
        match variant {
            HgvsVariant::Genome(g) => {
                let sequence = g.accession.to_string();
                let interval = &g.loc_edit.location;

                let edit =
                    g.loc_edit
                        .edit
                        .inner()
                        .ok_or_else(|| ConversionError::InvalidPosition {
                            description: "cannot convert variant with unknown edit".to_string(),
                        })?;

                let start_pos =
                    get_start_pos(interval).ok_or_else(|| ConversionError::InvalidPosition {
                        description: "cannot convert variant with unknown start position"
                            .to_string(),
                    })?;
                let end_pos = get_end_pos(interval).unwrap_or(start_pos);

                // Convert 1-based HGVS position to 0-based SPDI position
                let start_pos_ob = OneBasedPos::new(start_pos);
                let spdi_pos = start_pos_ob.to_zero_based().value();

                match edit {
                    NaEdit::Deletion { sequence: None, .. } => {
                        // Fetch sequence from reference
                        let del_seq = reference
                            .get_sequence(&sequence, spdi_pos, end_pos)
                            .ok_or_else(|| ConversionError::MissingReferenceData {
                                description: format!(
                                    "could not fetch sequence for {}:{}-{}",
                                    sequence, start_pos, end_pos
                                ),
                            })?;
                        Ok(SpdiVariant::new(sequence, spdi_pos, del_seq, ""))
                    }
                    NaEdit::Delins { sequence: ins_seq } => {
                        // Fetch deleted sequence from reference
                        let del_seq = reference
                            .get_sequence(&sequence, spdi_pos, end_pos)
                            .ok_or_else(|| ConversionError::MissingReferenceData {
                                description: format!(
                                    "could not fetch sequence for {}:{}-{}",
                                    sequence, start_pos, end_pos
                                ),
                            })?;
                        let ins_str = inserted_sequence_to_string(ins_seq).ok_or_else(|| {
                            ConversionError::MissingReferenceData {
                                description: "delins inserted sequence is not a literal"
                                    .to_string(),
                            }
                        })?;
                        Ok(SpdiVariant::new(sequence, spdi_pos, del_seq, ins_str))
                    }
                    NaEdit::Duplication { sequence: None, .. } => {
                        // Fetch duplicated sequence from reference
                        let dup_seq = reference
                            .get_sequence(&sequence, spdi_pos, end_pos)
                            .ok_or_else(|| ConversionError::MissingReferenceData {
                                description: format!(
                                    "could not fetch sequence for {}:{}-{}",
                                    sequence, start_pos, end_pos
                                ),
                            })?;
                        // Duplication becomes insertion after the region
                        // Convert 1-based HGVS end position to 0-based SPDI position
                        let end_pos_ob = OneBasedPos::new(end_pos);
                        let spdi_end = end_pos_ob.to_zero_based().value();
                        Ok(SpdiVariant::new(sequence, spdi_end, "", dup_seq))
                    }
                    // Fall back to simple conversion for other types
                    _ => hgvs_to_spdi_simple(variant),
                }
            }
            _ => hgvs_to_spdi_simple(variant),
        }
    }

    #[test]
    fn test_deletion_without_seq_with_reference() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();

        // This would fail with simple conversion
        assert!(hgvs_to_spdi_simple(&hgvs).is_err());

        // But works with reference data
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        assert_eq!(spdi.position, 99);
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_delins_with_reference_provides_actual_deletion() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();

        // Simple conversion without reference data returns an error
        let simple = hgvs_to_spdi_simple(&hgvs);
        assert!(simple.is_err());
        assert!(simple.unwrap_err().to_string().contains("unknown"));

        // Reference-based conversion gets actual sequence
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");
    }

    #[test]
    fn test_delins_roundtrip_with_reference() {
        let reference = MockGenomicRef::new();
        let original_hgvs = parse_hgvs("NC_000001.11:g.100_102delinsTTCC").unwrap();

        // Convert to SPDI with reference
        let spdi = hgvs_to_spdi_with_ref(&original_hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "TTCC");

        // Convert back to HGVS
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delinsTTCC");
    }

    #[test]
    fn test_deletion_roundtrip_with_reference() {
        let reference = MockGenomicRef::new();
        let original_hgvs = parse_hgvs("NC_000001.11:g.100_102del").unwrap();

        // Convert to SPDI with reference
        let spdi = hgvs_to_spdi_with_ref(&original_hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");

        // Convert back to HGVS
        let back = spdi_to_hgvs(&spdi).unwrap();
        // Note: HGVS back includes the sequence
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delATG");
    }

    #[test]
    fn test_duplication_without_seq_with_reference() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102dup").unwrap();

        // This would fail with simple conversion
        assert!(hgvs_to_spdi_simple(&hgvs).is_err());

        // But works with reference data
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.sequence, "NC_000001.11");
        // Dup is insertion after the region
        assert_eq!(spdi.position, 101); // end position 0-based
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_long_deletion_with_reference() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.1000_1009del").unwrap();

        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "AAACCCGGGT");
        assert_eq!(spdi.deletion.len(), 10);
    }

    #[test]
    fn test_reference_missing_region() {
        let reference = MockGenomicRef::new();
        // Position not in mock reference
        let hgvs = parse_hgvs("NC_000001.11:g.50000_50002del").unwrap();

        let result = hgvs_to_spdi_with_ref(&hgvs, &reference);
        assert!(matches!(
            result,
            Err(ConversionError::MissingReferenceData { .. })
        ));
    }

    #[test]
    fn test_substitution_falls_through_to_simple() {
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();

        // Substitution doesn't need reference, falls through
        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");
    }

    #[test]
    fn test_mnv_delins_with_reference() {
        // Multi-nucleotide variant (MNV) as delins
        let reference = MockGenomicRef::new();
        let hgvs = parse_hgvs("NC_000001.11:g.100_102delinsGGG").unwrap();

        let spdi = hgvs_to_spdi_with_ref(&hgvs, &reference).unwrap();
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "GGG");

        // MNV same length should roundtrip
        let back = spdi_to_hgvs(&spdi).unwrap();
        assert_eq!(back.to_string(), "NC_000001.11:g.100_102delinsGGG");
    }

    // =========================================================================
    // P3: SPDI edge case tests
    // =========================================================================

    #[test]
    fn test_spdi_empty_deletion_insertion() {
        // Pure insertion has empty deletion
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert!(spdi.is_insertion());
        assert!(!spdi.is_deletion());
        assert!(!spdi.is_identity());
        assert_eq!(spdi.deletion, "");
        assert_eq!(spdi.insertion, "ATG");
    }

    #[test]
    fn test_spdi_empty_insertion_deletion() {
        // Pure deletion has empty insertion
        let spdi = SpdiVariant::deletion("NC_000001.11", 100, "ATG");
        assert!(spdi.is_deletion());
        assert!(!spdi.is_insertion());
        assert!(!spdi.is_identity());
        assert_eq!(spdi.deletion, "ATG");
        assert_eq!(spdi.insertion, "");
    }

    #[test]
    fn test_spdi_both_empty_is_identity() {
        // Both empty is identity at position (unusual but valid)
        let spdi = SpdiVariant::new("NC_000001.11", 100, "", "");
        assert!(spdi.is_identity());
        assert!(!spdi.is_insertion());
        assert!(!spdi.is_deletion());
    }

    #[test]
    fn test_spdi_single_base_insertion() {
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "A");
        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 1);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("ins"));
    }

    #[test]
    fn test_spdi_single_base_deletion() {
        let spdi = SpdiVariant::deletion("NC_000001.11", 100, "A");
        assert!(spdi.is_deletion());
        assert_eq!(spdi.deletion.len(), 1);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("del"));
    }

    #[test]
    fn test_spdi_long_insertion_100bp() {
        // Test 100bp insertion (common structural variant size)
        let long_seq = "A".repeat(100);
        let spdi = SpdiVariant::insertion("NC_000001.11", 12345, &long_seq);

        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 100);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("ins"));
        // Verify the sequence is preserved
        assert!(hgvs.to_string().ends_with(&format!("ins{}", long_seq)));
    }

    #[test]
    fn test_spdi_long_deletion_100bp() {
        // Test 100bp deletion
        let long_seq = "ACGT".repeat(25); // 100bp
        let spdi = SpdiVariant::deletion("NC_000001.11", 12345, &long_seq);

        assert!(spdi.is_deletion());
        assert_eq!(spdi.deletion.len(), 100);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("del"));
    }

    #[test]
    fn test_spdi_long_indel_asymmetric() {
        // Delete 50bp, insert 150bp (net +100bp)
        let del_seq = "A".repeat(50);
        let ins_seq = "G".repeat(150);
        let spdi = SpdiVariant::delins("NC_000001.11", 12345, &del_seq, &ins_seq);

        assert!(!spdi.is_insertion());
        assert!(!spdi.is_deletion());
        assert_eq!(spdi.deletion.len(), 50);
        assert_eq!(spdi.insertion.len(), 150);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("delins"));
    }

    #[test]
    fn test_spdi_very_long_insertion_1000bp() {
        // Test 1000bp insertion (larger structural variant)
        let long_seq = "ACGT".repeat(250); // 1000bp
        let spdi = SpdiVariant::insertion("NC_000001.11", 50000, &long_seq);

        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 1000);

        // Should still convert to HGVS
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("ins"));
    }

    #[test]
    fn test_spdi_position_zero() {
        // Position 0 is valid in SPDI (0-based)
        let spdi = SpdiVariant::new("NC_000001.11", 0, "A", "G");
        assert_eq!(spdi.position, 0);

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        // HGVS position should be 1 (1-based)
        assert!(hgvs.to_string().contains("g.1A>G"));
    }

    #[test]
    fn test_spdi_position_max() {
        // Very large position (near chromosome end)
        let spdi = SpdiVariant::new("NC_000001.11", 248956421, "A", "G");

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("248956422")); // 1-based
    }

    #[test]
    fn test_spdi_lowercase_sequence() {
        // SPDI should handle lowercase (though uppercase is standard)
        let spdi = SpdiVariant::new("NC_000001.11", 100, "a", "g");

        // The variant should work even with lowercase
        assert_eq!(spdi.deletion, "a");
        assert_eq!(spdi.insertion, "g");
    }

    #[test]
    fn test_spdi_mixed_case_sequence() {
        // Mixed case sequence
        let spdi = SpdiVariant::new("NC_000001.11", 100, "AtGc", "GcTa");

        assert_eq!(spdi.deletion, "AtGc");
        assert_eq!(spdi.insertion, "GcTa");
    }

    #[test]
    fn test_spdi_n_bases_in_sequence() {
        // N (unknown) bases in sequence
        let spdi = SpdiVariant::new("NC_000001.11", 100, "ANG", "TNC");

        assert_eq!(spdi.deletion, "ANG");
        assert_eq!(spdi.insertion, "TNC");
    }

    #[test]
    fn test_spdi_complex_repeat_sequence() {
        // Repeat sequence (e.g., microsatellite)
        let repeat = "CAG".repeat(30); // 90bp CAG repeat
        let spdi = SpdiVariant::insertion("NC_000004.12", 3074876, &repeat);

        assert!(spdi.is_insertion());
        assert_eq!(spdi.insertion.len(), 90);
    }

    #[test]
    fn test_spdi_to_hgvs_delins_single_base_del() {
        // Single base deletion with multi-base insertion
        let spdi = SpdiVariant::delins("NC_000001.11", 100, "A", "TTTT");

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        // Single base del + ins = delins at single position
        assert!(hgvs.to_string().contains("delinsTTTT"));
    }

    #[test]
    fn test_spdi_to_hgvs_delins_single_base_ins() {
        // Multi-base deletion with single base insertion
        let spdi = SpdiVariant::delins("NC_000001.11", 100, "AAAA", "T");

        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        assert!(hgvs.to_string().contains("delinsT"));
    }

    #[test]
    fn test_spdi_different_chromosome_formats() {
        // Various accession formats should work
        let test_cases = vec![
            ("NC_000001.11", "NC_000001.11"), // Standard RefSeq
            ("NC_000023.11", "NC_000023.11"), // X chromosome
            ("NC_000024.10", "NC_000024.10"), // Y chromosome
            ("NC_012920.1", "NC_012920.1"),   // Mitochondrial
        ];

        for (input_acc, expected_acc) in test_cases {
            let spdi = SpdiVariant::new(input_acc, 100, "A", "G");
            assert_eq!(spdi.sequence, expected_acc);

            let hgvs = spdi_to_hgvs(&spdi).unwrap();
            assert!(hgvs.to_string().starts_with(expected_acc));
        }
    }

    #[test]
    fn test_spdi_roundtrip_preserves_case_normalized() {
        // Create SPDI substitution (SNV) which roundtrips cleanly
        let spdi = SpdiVariant::new("NC_000001.11", 100, "A", "G");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();
        let back = hgvs_to_spdi_simple(&hgvs).unwrap();

        // Should preserve the sequence
        assert_eq!(back.deletion, "A");
        assert_eq!(back.insertion, "G");

        // Test multi-base delins - without reference data, the roundtrip
        // cannot reconstruct the deletion sequence and returns an error
        let spdi_delins = SpdiVariant::new("NC_000001.11", 100, "ACGT", "TGCA");
        let hgvs_delins = spdi_to_hgvs(&spdi_delins).unwrap();
        let back_delins = hgvs_to_spdi_simple(&hgvs_delins);

        // Without reference, this should fail since the deletion sequence is unknown
        assert!(back_delins.is_err());
    }

    #[test]
    fn test_spdi_empty_seq_insertion_roundtrip() {
        // Empty deletion (pure insertion) roundtrip
        let spdi = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        let hgvs = spdi_to_hgvs(&spdi).unwrap();

        // HGVS format for insertion is pos_pos+1insX
        assert!(hgvs.to_string().contains("ins"));

        // Roundtrip back
        let back = hgvs_to_spdi_simple(&hgvs).unwrap();
        assert_eq!(back.deletion, "");
        assert_eq!(back.insertion, "ATG");
    }

    #[test]
    fn test_hgvs_to_spdi_various_accession_types() {
        // Test conversion from various accession types
        let test_variants = vec![
            "NC_000001.11:g.12345A>G", // Chromosome
            "NC_000023.11:g.12345A>G", // X chromosome
            "NC_012920.1:g.12345A>G",  // Mitochondrial
        ];

        for variant_str in test_variants {
            let hgvs = parse_hgvs(variant_str).unwrap();
            let spdi = hgvs_to_spdi_simple(&hgvs);
            assert!(spdi.is_ok(), "Failed for: {}", variant_str);
        }
    }

    #[test]
    fn test_spdi_display_format() {
        let spdi = SpdiVariant::new("NC_000001.11", 12344, "A", "G");
        assert_eq!(spdi.to_string(), "NC_000001.11:12344:A:G");

        let spdi_del = SpdiVariant::deletion("NC_000001.11", 100, "ATG");
        assert_eq!(spdi_del.to_string(), "NC_000001.11:100:ATG:");

        let spdi_ins = SpdiVariant::insertion("NC_000001.11", 100, "ATG");
        assert_eq!(spdi_ins.to_string(), "NC_000001.11:100::ATG");
    }

    #[test]
    fn test_spdi_identity_various_lengths() {
        // Single base identity
        let spdi1 = SpdiVariant::new("NC_000001.11", 100, "A", "A");
        assert!(spdi1.is_identity());

        // Multi-base identity (unusual but valid)
        let spdi2 = SpdiVariant::new("NC_000001.11", 100, "ATG", "ATG");
        assert!(spdi2.is_identity());

        // Empty identity
        let spdi3 = SpdiVariant::new("NC_000001.11", 100, "", "");
        assert!(spdi3.is_identity());
    }
}
