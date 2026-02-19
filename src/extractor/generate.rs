//! HGVS string generation from classified edits.
//!
//! # Coordinate System
//!
//! | Field | Basis | Notes |
//! |-------|-------|-------|
//! | `ExtractedVariant.position` | 1-based | HGVS g. position |
//! | `RawEdit.ref_start`, `RawEdit.ref_end` | 1-based | From alignment module |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use super::align::RawEdit;
use super::classify::EditClassification;

/// Result of extracting variants from sequences.
#[derive(Debug, Clone)]
pub struct ExtractionResult {
    /// Length of the reference sequence.
    pub reference_length: u64,
    /// Length of the observed sequence.
    pub observed_length: u64,
    /// All extracted variants.
    pub variants: Vec<ExtractedVariant>,
    /// HGVS strings for all variants.
    pub hgvs_strings: Vec<String>,
}

/// A single extracted variant.
#[derive(Debug, Clone)]
pub struct ExtractedVariant {
    /// Position in reference (1-based).
    pub position: u64,
    /// Reference sequence at this position.
    pub ref_seq: String,
    /// Observed sequence at this position.
    pub obs_seq: String,
    /// Classification of the edit.
    pub classification: EditClassification,
    /// HGVS string representation.
    pub hgvs: String,
}

/// Generate an HGVS string from an edit and its classification.
///
/// # Arguments
///
/// * `edit` - The raw edit with position information
/// * `classification` - The classified edit type
/// * `accession` - Optional accession to prefix (e.g., "NC_000001.11")
pub fn generate_hgvs(
    edit: &RawEdit,
    classification: &EditClassification,
    accession: Option<&str>,
) -> String {
    let position_str = format_position(edit);
    let edit_str = format_edit(classification);

    match accession {
        Some(acc) => format!("{}:g.{}{}", acc, position_str, edit_str),
        None => format!("{}{}", position_str, edit_str),
    }
}

/// Format the position part of the HGVS string.
fn format_position(edit: &RawEdit) -> String {
    let ref_len = edit.ref_seq.len();
    let obs_len = edit.obs_seq.len();

    match (ref_len, obs_len) {
        // Substitution or single deletion
        (1, _) => format!("{}", edit.ref_start),

        // Multi-base deletion or delins
        (r, _) if r > 1 => format!("{}_{}", edit.ref_start, edit.ref_start + r as u64 - 1),

        // Insertion - position is between bases
        (0, _) => {
            if edit.ref_start == 0 {
                // Insertion at very beginning
                "0_1".to_string()
            } else {
                format!("{}_{}", edit.ref_start, edit.ref_start + 1)
            }
        }

        _ => format!("{}", edit.ref_start),
    }
}

/// Format the edit part of the HGVS string.
fn format_edit(classification: &EditClassification) -> String {
    match classification {
        EditClassification::Substitution { ref_base, alt_base } => {
            format!("{}>{}", ref_base, alt_base)
        }

        EditClassification::Deletion { deleted: _ } => "del".to_string(),

        EditClassification::Insertion { inserted } => {
            format!("ins{}", inserted)
        }

        EditClassification::Duplication { sequence: _ } => "dup".to_string(),

        EditClassification::Delins {
            deleted: _,
            inserted,
        } => {
            format!("delins{}", inserted)
        }

        EditClassification::Inversion { sequence: _ } => "inv".to_string(),

        EditClassification::Repeat { unit, count } => {
            format!("{}[{}]", unit, count)
        }
    }
}

/// Format a duplication with proper positioning.
///
/// Duplications in HGVS should reference the original sequence being duplicated,
/// not the insertion point.
#[allow(dead_code)]
pub fn format_duplication_hgvs(ref_start: u64, sequence: &str, accession: Option<&str>) -> String {
    let len = sequence.len() as u64;
    let dup_start = ref_start - len;
    let dup_end = ref_start - 1;

    let position_str = if len == 1 {
        format!("{}", dup_start)
    } else {
        format!("{}_{}", dup_start, dup_end)
    };

    match accession {
        Some(acc) => format!("{}:g.{}dup", acc, position_str),
        None => format!("{}dup", position_str),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_edit(ref_start: u64, ref_seq: &str, obs_seq: &str) -> RawEdit {
        RawEdit {
            ref_start,
            ref_end: ref_start + ref_seq.len() as u64,
            obs_start: ref_start,
            obs_end: ref_start + obs_seq.len() as u64,
            ref_seq: ref_seq.to_string(),
            obs_seq: obs_seq.to_string(),
        }
    }

    #[test]
    fn test_generate_substitution() {
        let edit = make_edit(3, "G", "A");
        let class = EditClassification::Substitution {
            ref_base: 'G',
            alt_base: 'A',
        };

        assert_eq!(generate_hgvs(&edit, &class, None), "3G>A");
    }

    #[test]
    fn test_generate_substitution_with_accession() {
        let edit = make_edit(3, "G", "A");
        let class = EditClassification::Substitution {
            ref_base: 'G',
            alt_base: 'A',
        };

        assert_eq!(
            generate_hgvs(&edit, &class, Some("NC_000001.11")),
            "NC_000001.11:g.3G>A"
        );
    }

    #[test]
    fn test_generate_single_deletion() {
        let edit = make_edit(3, "G", "");
        let class = EditClassification::Deletion {
            deleted: "G".to_string(),
        };

        assert_eq!(generate_hgvs(&edit, &class, None), "3del");
    }

    #[test]
    fn test_generate_multi_deletion() {
        let edit = make_edit(3, "GC", "");
        let class = EditClassification::Deletion {
            deleted: "GC".to_string(),
        };

        assert_eq!(generate_hgvs(&edit, &class, None), "3_4del");
    }

    #[test]
    fn test_generate_insertion() {
        let edit = make_edit(3, "", "T");
        let class = EditClassification::Insertion {
            inserted: "T".to_string(),
        };

        assert_eq!(generate_hgvs(&edit, &class, None), "3_4insT");
    }

    #[test]
    fn test_generate_delins() {
        let edit = make_edit(3, "GC", "TT");
        let class = EditClassification::Delins {
            deleted: "GC".to_string(),
            inserted: "TT".to_string(),
        };

        assert_eq!(generate_hgvs(&edit, &class, None), "3_4delinsTT");
    }

    #[test]
    fn test_generate_inversion() {
        let edit = make_edit(3, "ATG", "CAT");
        let class = EditClassification::Inversion {
            sequence: "ATG".to_string(),
        };

        assert_eq!(generate_hgvs(&edit, &class, None), "3_5inv");
    }

    #[test]
    fn test_format_duplication() {
        // Dup of position 3 (single base)
        assert_eq!(format_duplication_hgvs(4, "G", None), "3dup");

        // Dup of positions 2-3 (two bases)
        assert_eq!(format_duplication_hgvs(4, "TG", None), "2_3dup");

        // With accession
        assert_eq!(
            format_duplication_hgvs(4, "G", Some("NC_000001.11")),
            "NC_000001.11:g.3dup"
        );
    }
}
