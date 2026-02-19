//! Edit classification for HGVS generation.
//!
//! Classifies raw edits into HGVS-compatible types.

use super::align::RawEdit;

/// Classification of an edit for HGVS notation.
#[derive(Debug, Clone, PartialEq)]
pub enum EditClassification {
    /// Single nucleotide substitution (e.g., A>G).
    Substitution {
        /// Reference base.
        ref_base: char,
        /// Alternate base.
        alt_base: char,
    },
    /// Deletion of one or more bases.
    Deletion {
        /// Deleted sequence.
        deleted: String,
    },
    /// Insertion of one or more bases.
    Insertion {
        /// Inserted sequence.
        inserted: String,
    },
    /// Deletion-insertion (replacement).
    Delins {
        /// Deleted sequence.
        deleted: String,
        /// Inserted sequence.
        inserted: String,
    },
    /// Duplication of adjacent sequence.
    Duplication {
        /// Duplicated sequence.
        sequence: String,
    },
    /// Inversion of sequence.
    Inversion {
        /// Inverted sequence.
        sequence: String,
    },
    /// Repeated sequence with count.
    Repeat {
        /// Repeat unit.
        unit: String,
        /// Number of repeats.
        count: u32,
    },
}

/// Classify a raw edit into an HGVS-compatible type.
///
/// # Arguments
///
/// * `edit` - The raw edit to classify
/// * `reference` - The full reference sequence (for duplication detection)
/// * `detect_dup` - Whether to detect duplications
/// * `detect_inv` - Whether to detect inversions
pub fn classify_edit(
    edit: &RawEdit,
    reference: &str,
    detect_dup: bool,
    detect_inv: bool,
) -> EditClassification {
    let ref_len = edit.ref_seq.len();
    let obs_len = edit.obs_seq.len();

    match (ref_len, obs_len) {
        // Single substitution
        (1, 1) => EditClassification::Substitution {
            ref_base: edit.ref_seq.chars().next().unwrap(),
            alt_base: edit.obs_seq.chars().next().unwrap(),
        },

        // Deletion
        (r, 0) if r > 0 => EditClassification::Deletion {
            deleted: edit.ref_seq.clone(),
        },

        // Insertion - check for duplication
        (0, a) if a > 0 => {
            if detect_dup && is_duplication(reference, edit.ref_start, &edit.obs_seq) {
                EditClassification::Duplication {
                    sequence: edit.obs_seq.clone(),
                }
            } else {
                EditClassification::Insertion {
                    inserted: edit.obs_seq.clone(),
                }
            }
        }

        // Delins - check for inversion
        (r, a) if r > 0 && a > 0 => {
            if detect_inv && r == a && is_inversion(&edit.ref_seq, &edit.obs_seq) {
                EditClassification::Inversion {
                    sequence: edit.ref_seq.clone(),
                }
            } else if r == 1 && a == 1 {
                // Single base change
                EditClassification::Substitution {
                    ref_base: edit.ref_seq.chars().next().unwrap(),
                    alt_base: edit.obs_seq.chars().next().unwrap(),
                }
            } else {
                EditClassification::Delins {
                    deleted: edit.ref_seq.clone(),
                    inserted: edit.obs_seq.clone(),
                }
            }
        }

        _ => {
            // Fallback - shouldn't happen
            EditClassification::Delins {
                deleted: edit.ref_seq.clone(),
                inserted: edit.obs_seq.clone(),
            }
        }
    }
}

/// Check if an insertion is a duplication of the preceding sequence.
///
/// According to HGVS, an insertion is a duplication if the inserted sequence
/// is identical to the sequence immediately 5' of the insertion point.
pub fn is_duplication(reference: &str, position: u64, inserted: &str) -> bool {
    if inserted.is_empty() {
        return false;
    }

    let pos = position as usize;
    let ins_len = inserted.len();

    // Check if there's enough sequence before the insertion point
    if pos < ins_len {
        return false;
    }

    // Get the sequence immediately before the insertion point
    let preceding = &reference[pos - ins_len..pos];

    preceding == inserted
}

/// Check if a delins is actually an inversion.
///
/// An inversion occurs when the inserted sequence is the reverse complement
/// of the deleted sequence.
pub fn is_inversion(deleted: &str, inserted: &str) -> bool {
    if deleted.len() != inserted.len() || deleted.is_empty() {
        return false;
    }

    // Check if inserted is reverse complement of deleted
    let rev_comp = reverse_complement(deleted);
    rev_comp == inserted
}

/// Compute the reverse complement of a DNA sequence.
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            other => other,
        })
        .collect()
}

/// Find a repeat unit in a sequence.
///
/// Returns the unit and count if the sequence is composed of repeating units.
#[allow(dead_code)]
pub fn find_repeat_unit(sequence: &str) -> Option<(String, u32)> {
    if sequence.is_empty() {
        return None;
    }

    let len = sequence.len();

    // Try unit lengths from 1 to half the sequence length
    for unit_len in 1..=len / 2 {
        if !len.is_multiple_of(unit_len) {
            continue;
        }

        let unit = &sequence[0..unit_len];
        let count = len / unit_len;

        // Check if entire sequence is this unit repeated
        let is_repeat = (0..count).all(|i| {
            let start = i * unit_len;
            &sequence[start..start + unit_len] == unit
        });

        if is_repeat {
            return Some((unit.to_string(), count as u32));
        }
    }

    None
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
    fn test_classify_substitution() {
        let edit = make_edit(3, "G", "A");
        let class = classify_edit(&edit, "ATGC", true, true);

        match class {
            EditClassification::Substitution { ref_base, alt_base } => {
                assert_eq!(ref_base, 'G');
                assert_eq!(alt_base, 'A');
            }
            _ => panic!("Expected Substitution"),
        }
    }

    #[test]
    fn test_classify_deletion() {
        let edit = make_edit(3, "GC", "");
        let class = classify_edit(&edit, "ATGCAT", true, true);

        match class {
            EditClassification::Deletion { deleted } => {
                assert_eq!(deleted, "GC");
            }
            _ => panic!("Expected Deletion"),
        }
    }

    #[test]
    fn test_classify_insertion() {
        let edit = make_edit(3, "", "T");
        let class = classify_edit(&edit, "ATGC", true, true);

        match class {
            EditClassification::Insertion { inserted } => {
                assert_eq!(inserted, "T");
            }
            _ => panic!("Expected Insertion"),
        }
    }

    #[test]
    fn test_classify_duplication() {
        // Insert "G" at position 3 (1-based) in "ATGC"
        // Position 3 is 'G', so inserting "G" after it is a duplication
        // The preceding character (at index 2, position 3) is 'G'
        let edit = make_edit(3, "", "G");
        let class = classify_edit(&edit, "ATGC", true, true);

        match class {
            EditClassification::Duplication { sequence } => {
                assert_eq!(sequence, "G");
            }
            _ => panic!("Expected Duplication, got {:?}", class),
        }
    }

    #[test]
    fn test_classify_delins() {
        let edit = make_edit(3, "GC", "TT");
        let class = classify_edit(&edit, "ATGCAT", true, true);

        match class {
            EditClassification::Delins { deleted, inserted } => {
                assert_eq!(deleted, "GC");
                assert_eq!(inserted, "TT");
            }
            _ => panic!("Expected Delins"),
        }
    }

    #[test]
    fn test_is_duplication() {
        // "ATGCAT" - inserting "CAT" at position 7 (after the full sequence)
        // Preceding "CAT" (positions 4-6) matches
        assert!(is_duplication("ATGCAT", 6, "CAT"));

        // Not a duplication - different sequence
        assert!(!is_duplication("ATGCAT", 6, "TTT"));

        // Not enough preceding sequence
        assert!(!is_duplication("AT", 2, "ATGC"));
    }

    #[test]
    fn test_is_inversion() {
        // "ATG" inverted is "CAT" (reverse complement)
        assert!(is_inversion("ATG", "CAT"));

        // Not an inversion
        assert!(!is_inversion("ATG", "TTT"));

        // Different lengths
        assert!(!is_inversion("ATG", "CA"));
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATGC"), "GCAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn test_find_repeat_unit() {
        assert_eq!(find_repeat_unit("ATAT"), Some(("AT".to_string(), 2)));
        assert_eq!(find_repeat_unit("AAA"), Some(("A".to_string(), 3)));
        assert_eq!(find_repeat_unit("ATCATCATC"), Some(("ATC".to_string(), 3)));
        assert_eq!(find_repeat_unit("ATGC"), None); // No repeating unit
        assert_eq!(find_repeat_unit(""), None);
    }
}
