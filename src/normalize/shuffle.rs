//! 3'/5' shuffling algorithm
//!
//! Implementation of the variant shuffling algorithm for normalization.
//!
//! # Coordinate System
//!
//! This module uses **0-based half-open intervals** internally:
//!
//! | Parameter | Basis | Notes |
//! |-----------|-------|-------|
//! | `start` | 0-based | Inclusive start position |
//! | `end` | 0-based | Exclusive end position |
//! | `boundaries.left` | 0-based | Inclusive left limit |
//! | `boundaries.right` | 0-based | Exclusive right limit |
//!
//! Callers must convert from 1-based HGVS positions before calling `shuffle()`.
//! Use `hgvs_pos_to_index()` from [`crate::coords`] for this conversion.

use crate::normalize::boundary::Boundaries;
use crate::normalize::config::ShuffleDirection;

/// Result of a shuffle operation
///
/// All positions are 0-based, matching the input coordinate system.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ShuffleResult {
    /// New start position (0-based, inclusive)
    pub start: u64,
    /// New end position (0-based, exclusive)
    pub end: u64,
    /// Whether the variant was moved
    pub shifted: bool,
}

/// Shuffle a variant towards 3' or 5' end
///
/// This is the core algorithm for HGVS normalization.
/// For insertions and deletions in repetitive regions,
/// the variant position is shifted to the rightmost (3')
/// or leftmost (5') valid position.
///
/// # Arguments
///
/// * `ref_seq` - Reference sequence bytes
/// * `alt_seq` - Alternate sequence bytes (empty for deletions)
/// * `start` - 0-based start position
/// * `end` - 0-based end position (exclusive)
/// * `boundaries` - Shuffling boundaries
/// * `direction` - Direction to shuffle (3' or 5')
///
/// # Returns
///
/// New start and end positions after shuffling
pub fn shuffle(
    ref_seq: &[u8],
    alt_seq: &[u8],
    start: u64,
    end: u64,
    boundaries: &Boundaries,
    direction: ShuffleDirection,
) -> ShuffleResult {
    let mut new_start = start;
    let mut new_end = end;

    match direction {
        ShuffleDirection::ThreePrime => {
            // Shuffle right (3')
            while new_end < boundaries.right {
                let ref_idx = new_end as usize;
                if ref_idx >= ref_seq.len() {
                    break;
                }

                // For deletion: check if we can shift right
                if alt_seq.is_empty() {
                    // The base at the end of the deletion matches the base after
                    let del_start_idx = new_start as usize;
                    if del_start_idx < ref_seq.len() && ref_seq[del_start_idx] == ref_seq[ref_idx] {
                        new_start += 1;
                        new_end += 1;
                    } else {
                        break;
                    }
                } else {
                    // For insertion: check if inserted sequence can shift
                    let alt_idx = ((new_end - start) % alt_seq.len() as u64) as usize;
                    if ref_seq[ref_idx] == alt_seq[alt_idx] {
                        new_start += 1;
                        new_end += 1;
                    } else {
                        break;
                    }
                }
            }
        }
        ShuffleDirection::FivePrime => {
            // Shuffle left (5')
            while new_start > boundaries.left {
                let check_idx = new_start - 1;
                let ref_idx = check_idx as usize;
                if ref_idx >= ref_seq.len() {
                    break;
                }

                // For deletion: check if we can shift left
                if alt_seq.is_empty() {
                    let del_end_idx = (new_end - 1) as usize;
                    if del_end_idx < ref_seq.len() && ref_seq[del_end_idx] == ref_seq[ref_idx] {
                        new_start -= 1;
                        new_end -= 1;
                    } else {
                        break;
                    }
                } else {
                    // For insertion: check if inserted sequence can shift
                    let alt_idx =
                        alt_seq.len() - 1 - ((start - new_start) % alt_seq.len() as u64) as usize;
                    if ref_idx < ref_seq.len() && ref_seq[ref_idx] == alt_seq[alt_idx] {
                        new_start -= 1;
                        new_end -= 1;
                    } else {
                        break;
                    }
                }
            }
        }
    }

    ShuffleResult {
        start: new_start,
        end: new_end,
        shifted: new_start != start,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shuffle_deletion_3prime() {
        // Sequence: ATGGGGGCAT
        // Delete one G at position 3 (0-based): should shift to position 6
        let ref_seq = b"ATGGGGGCAT";
        let boundaries = Boundaries::new(0, 10);

        let result = shuffle(
            ref_seq,
            &[], // deletion
            3,   // start
            4,   // end (exclusive)
            &boundaries,
            ShuffleDirection::ThreePrime,
        );

        assert!(result.shifted);
        assert_eq!(result.start, 6);
        assert_eq!(result.end, 7);
    }

    #[test]
    fn test_shuffle_deletion_5prime() {
        // Sequence: ATGGGGGCAT
        // Delete one G at position 6 (0-based): should shift to position 2
        let ref_seq = b"ATGGGGGCAT";
        let boundaries = Boundaries::new(0, 10);

        let result = shuffle(
            ref_seq,
            &[], // deletion
            6,   // start
            7,   // end
            &boundaries,
            ShuffleDirection::FivePrime,
        );

        assert!(result.shifted);
        assert_eq!(result.start, 2);
        assert_eq!(result.end, 3);
    }

    #[test]
    fn test_no_shuffle_needed() {
        // Sequence: ATGCATGCAT
        // Delete T at position 2: no adjacent T's to shift to
        let ref_seq = b"ATGCATGCAT";
        let boundaries = Boundaries::new(0, 10);

        let result = shuffle(
            ref_seq,
            &[],
            2,
            3,
            &boundaries,
            ShuffleDirection::ThreePrime,
        );

        assert!(!result.shifted);
        assert_eq!(result.start, 2);
        assert_eq!(result.end, 3);
    }

    #[test]
    fn test_shuffle_respects_boundary() {
        // Sequence: ATGGGGGCAT
        // Delete G, but boundary prevents full shift
        let ref_seq = b"ATGGGGGCAT";
        let boundaries = Boundaries::new(0, 5); // Restrict right boundary

        let result = shuffle(
            ref_seq,
            &[],
            3,
            4,
            &boundaries,
            ShuffleDirection::ThreePrime,
        );

        assert!(result.shifted);
        assert_eq!(result.end, 5); // Stopped at boundary
    }
}
