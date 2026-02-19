//! Sequence alignment for variant detection.
//!
//! Implements a simple diff algorithm to find differences between sequences.
//!
//! # Coordinate System
//!
//! | Field | Basis | Notes |
//! |-------|-------|-------|
//! | `RawEdit.ref_start`, `RawEdit.ref_end` | 1-based | Converted from 0-based via `index_to_hgvs_pos()` |
//! | `RawEdit.obs_start`, `RawEdit.obs_end` | 1-based | Converted from 0-based via `index_to_hgvs_pos()` |
//! | Internal array indices | 0-based | Used during alignment |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::coords::index_to_hgvs_pos;

/// Raw edit from sequence comparison.
#[derive(Debug, Clone, PartialEq)]
pub struct RawEdit {
    /// Start position in reference (1-based).
    pub ref_start: u64,
    /// End position in reference (1-based, inclusive).
    pub ref_end: u64,
    /// Start position in observed (1-based).
    pub obs_start: u64,
    /// End position in observed (1-based, inclusive).
    pub obs_end: u64,
    /// Reference sequence at this position.
    pub ref_seq: String,
    /// Observed sequence at this position.
    pub obs_seq: String,
}

/// Edit operation type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    /// Match/equal.
    Match,
    /// Substitution.
    Replace,
    /// Insertion in observed.
    Insert,
    /// Deletion from reference.
    Delete,
}

/// Align two sequences and find all differences.
///
/// Uses a simple O(n) algorithm that finds common prefix/suffix and then
/// identifies the differing region. For more complex cases with multiple
/// variants, it iterates through to find all changes.
///
/// # Arguments
///
/// * `reference` - The reference sequence
/// * `observed` - The observed (variant) sequence
///
/// # Returns
///
/// A vector of `RawEdit` describing all differences.
pub fn align_sequences(reference: &str, observed: &str) -> Vec<RawEdit> {
    if reference == observed {
        return Vec::new();
    }

    let ref_bytes = reference.as_bytes();
    let obs_bytes = observed.as_bytes();

    // Use a character-by-character comparison approach
    // This is O(n) and handles most common cases efficiently
    find_edits(ref_bytes, obs_bytes)
}

/// Find all edits between two sequences.
fn find_edits(reference: &[u8], observed: &[u8]) -> Vec<RawEdit> {
    let mut edits = Vec::new();
    let mut ref_pos = 0usize;
    let mut obs_pos = 0usize;

    while ref_pos < reference.len() || obs_pos < observed.len() {
        // Find next mismatch
        while ref_pos < reference.len()
            && obs_pos < observed.len()
            && reference[ref_pos] == observed[obs_pos]
        {
            ref_pos += 1;
            obs_pos += 1;
        }

        if ref_pos >= reference.len() && obs_pos >= observed.len() {
            break;
        }

        // We have a mismatch - determine the edit type
        let edit_start_ref = ref_pos;
        let edit_start_obs = obs_pos;

        // Look ahead to find where sequences re-synchronize
        let (ref_end, obs_end) = find_resync_point(reference, observed, ref_pos, obs_pos);

        // Convert byte slices to strings. DNA/RNA sequences are always ASCII
        // (A, C, G, T, U, N and similar), so byte-to-char conversion is safe.
        let ref_seq: String = reference[edit_start_ref..ref_end]
            .iter()
            .map(|&b| b as char)
            .collect();
        let obs_seq: String = observed[edit_start_obs..obs_end]
            .iter()
            .map(|&b| b as char)
            .collect();

        // Apply 3' normalization (right-shift) for insertions and deletions
        let (norm_ref_start, norm_ref_end, norm_ref_seq, norm_obs_seq) =
            normalize_3prime(reference, edit_start_ref, ref_end, &ref_seq, &obs_seq);

        // COORDINATE CONVERSION: 0-based half-open [start, end) → 1-based inclusive [start, end]
        // - Start: 0-based index + 1 = 1-based position
        // - End: 0-based exclusive = 1-based inclusive (no conversion needed)
        //   Because: exclusive end = last_index + 1 = 1-based position of last element
        edits.push(RawEdit {
            ref_start: index_to_hgvs_pos(norm_ref_start), // 0-based → 1-based
            ref_end: if norm_ref_end > norm_ref_start {
                norm_ref_end as u64 // 0-based exclusive = 1-based inclusive
            } else {
                index_to_hgvs_pos(norm_ref_start) // Single position case
            },
            obs_start: index_to_hgvs_pos(edit_start_obs), // 0-based → 1-based
            obs_end: if obs_end > edit_start_obs {
                obs_end as u64 // 0-based exclusive = 1-based inclusive
            } else {
                index_to_hgvs_pos(edit_start_obs) // Single position case
            },
            ref_seq: norm_ref_seq,
            obs_seq: norm_obs_seq,
        });

        ref_pos = ref_end;
        obs_pos = obs_end;
    }

    edits
}

/// Find the point where sequences re-synchronize after a mismatch.
fn find_resync_point(
    reference: &[u8],
    observed: &[u8],
    ref_start: usize,
    obs_start: usize,
) -> (usize, usize) {
    let ref_remaining = reference.len() - ref_start;
    let obs_remaining = observed.len() - obs_start;

    // Simple case: both sequences have content but mismatch
    if ref_remaining > 0 && obs_remaining > 0 {
        // Try to find a common suffix to re-sync
        // Look for the smallest edit that explains the difference

        // Check for simple substitution
        if ref_remaining >= 1 && obs_remaining >= 1 {
            // Check if next characters match
            let check_len = std::cmp::min(ref_remaining, obs_remaining);
            for i in 1..=check_len {
                if i < ref_remaining
                    && i < obs_remaining
                    && reference[ref_start + i] == observed[obs_start + i]
                {
                    // Found re-sync point - this is a substitution or small indel
                    return (ref_start + i, obs_start + i);
                }
            }
        }

        // Check for deletion (missing bases in observed)
        for del_len in 1..=ref_remaining {
            if ref_start + del_len < reference.len() {
                let ref_after_del = ref_start + del_len;
                if reference[ref_after_del..] == observed[obs_start..] {
                    return (ref_after_del, obs_start);
                }
                // Check partial match
                let remaining_ref = reference.len() - ref_after_del;
                let remaining_obs = observed.len() - obs_start;
                let check_len = std::cmp::min(remaining_ref, remaining_obs);
                if check_len > 0
                    && reference[ref_after_del..ref_after_del + check_len]
                        == observed[obs_start..obs_start + check_len]
                {
                    return (ref_start + del_len, obs_start);
                }
            }
        }

        // Check for insertion (extra bases in observed)
        for ins_len in 1..=obs_remaining {
            if obs_start + ins_len < observed.len() {
                let obs_after_ins = obs_start + ins_len;
                if reference[ref_start..] == observed[obs_after_ins..] {
                    return (ref_start, obs_after_ins);
                }
                // Check partial match
                let remaining_ref = reference.len() - ref_start;
                let remaining_obs = observed.len() - obs_after_ins;
                let check_len = std::cmp::min(remaining_ref, remaining_obs);
                if check_len > 0
                    && reference[ref_start..ref_start + check_len]
                        == observed[obs_after_ins..obs_after_ins + check_len]
                {
                    return (ref_start, obs_start + ins_len);
                }
            }
        }
    }

    // Fall through: consume rest of both sequences
    (reference.len(), observed.len())
}

/// Apply 3' normalization (right-shift) to the edit.
fn normalize_3prime(
    reference: &[u8],
    ref_start: usize,
    ref_end: usize,
    ref_seq: &str,
    obs_seq: &str,
) -> (usize, usize, String, String) {
    // Only normalize pure insertions or deletions
    if !ref_seq.is_empty() && !obs_seq.is_empty() {
        // This is a substitution or delins - no shifting
        return (ref_start, ref_end, ref_seq.to_string(), obs_seq.to_string());
    }

    let seq_to_shift = if ref_seq.is_empty() { obs_seq } else { ref_seq };

    if seq_to_shift.is_empty() {
        return (ref_start, ref_end, ref_seq.to_string(), obs_seq.to_string());
    }

    let shift_bytes = seq_to_shift.as_bytes();
    let mut new_start = ref_start;
    let mut new_end = ref_end;

    // Shift right while the next base matches the first base of the indel
    while new_end < reference.len() && reference[new_end] == shift_bytes[0] {
        // Rotate the sequence
        let first = shift_bytes[0];
        if shift_bytes.iter().all(|&b| b == first) {
            // All bases are the same - just shift position
            new_start += 1;
            new_end += 1;
        } else {
            break;
        }
    }

    if ref_seq.is_empty() {
        (new_start, new_end, String::new(), obs_seq.to_string())
    } else {
        (new_start, new_end, ref_seq.to_string(), String::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identical_sequences() {
        let edits = align_sequences("ATGC", "ATGC");
        assert!(edits.is_empty());
    }

    #[test]
    fn test_single_substitution() {
        let edits = align_sequences("ATGC", "ATAC");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_start, 3);
        assert_eq!(edits[0].ref_seq, "G");
        assert_eq!(edits[0].obs_seq, "A");
    }

    #[test]
    fn test_single_deletion() {
        let edits = align_sequences("ATGC", "ATC");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "G");
        assert_eq!(edits[0].obs_seq, "");
    }

    #[test]
    fn test_single_insertion() {
        let edits = align_sequences("ATGC", "ATGGC");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "");
        assert_eq!(edits[0].obs_seq, "G");
    }

    #[test]
    fn test_multi_base_deletion() {
        let edits = align_sequences("ATGCAT", "ATAT");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "GC");
        assert_eq!(edits[0].obs_seq, "");
    }

    #[test]
    fn test_delins() {
        let edits = align_sequences("ATGC", "ATTC");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "G");
        assert_eq!(edits[0].obs_seq, "T");
    }

    #[test]
    fn test_multiple_substitutions() {
        let edits = align_sequences("ATGCAT", "CTGCGT");
        assert_eq!(edits.len(), 2);
        assert_eq!(edits[0].ref_start, 1);
        assert_eq!(edits[0].ref_seq, "A");
        assert_eq!(edits[0].obs_seq, "C");
        assert_eq!(edits[1].ref_seq, "A");
        assert_eq!(edits[1].obs_seq, "G");
    }

    #[test]
    fn test_complete_replacement() {
        let edits = align_sequences("AAAA", "TTTT");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "AAAA");
        assert_eq!(edits[0].obs_seq, "TTTT");
    }

    #[test]
    fn test_insertion_at_end() {
        let edits = align_sequences("ATG", "ATGC");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "");
        assert_eq!(edits[0].obs_seq, "C");
    }

    #[test]
    fn test_deletion_at_end() {
        let edits = align_sequences("ATGC", "ATG");
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].ref_seq, "C");
        assert_eq!(edits[0].obs_seq, "");
    }
}
