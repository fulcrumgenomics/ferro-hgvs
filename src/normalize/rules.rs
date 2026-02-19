//! HGVS normalization rules
//!
//! Additional rules for HGVS-compliant normalization beyond shuffling.
//!
//! # Coordinate Systems
//!
//! This module uses two coordinate systems:
//!
//! | Context | Basis | Type |
//! |---------|-------|------|
//! | Array indexing | 0-based | `usize` |
//! | HGVS output positions | 1-based | `u64` |
//!
//! Functions that take `pos: u64` parameters use 0-based positions for internal
//! array indexing. Return values marked with "1-indexed" comments use HGVS-style
//! 1-based positions.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::coords::index_to_hgvs_pos;
use crate::hgvs::edit::NaEdit;

/// Check if an edit type needs normalization
pub fn needs_normalization(edit: &NaEdit) -> bool {
    matches!(
        edit,
        NaEdit::Deletion { .. }
            | NaEdit::Insertion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Repeat { .. }
    )
}

/// Check if an insertion should be represented as a duplication
///
/// In HGVS, if an inserted sequence is identical to the sequence
/// immediately 5' (before) or 3' (after) of the insertion point,
/// the variant should be represented as a duplication.
///
/// IMPORTANT: This only applies to INSERTIONS, not deletions.
/// Deletions always stay as deletions - they just shift 3'.
pub fn insertion_is_duplication(ref_seq: &[u8], pos: u64, inserted_seq: &[u8]) -> bool {
    let ins_len = inserted_seq.len();
    let pos_idx = pos as usize;

    // Check if sequence before the insertion point matches the inserted sequence
    // (this would be a duplication of the preceding sequence)
    if pos_idx >= ins_len {
        let before_start = pos_idx - ins_len;
        if ref_seq[before_start..pos_idx] == inserted_seq[..] {
            return true;
        }
    }

    // Also check if sequence after matches (for 5' duplications)
    if pos_idx + ins_len <= ref_seq.len() && ref_seq[pos_idx..pos_idx + ins_len] == inserted_seq[..]
    {
        return true;
    }

    false
}

/// Determine the canonical representation of an indel
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CanonicalForm {
    /// Standard deletion
    Deletion,
    /// Duplication (del of a repeated sequence)
    Duplication,
    /// Deletion-insertion
    Delins,
    /// Standard insertion
    Insertion,
    /// Repeat notation (e.g., A[9] for homopolymer)
    Repeat {
        /// The repeated base
        base: u8,
        /// Total count after the edit
        count: u64,
        /// Start position (1-based HGVS, inclusive) of the repeat region
        start: u64,
        /// End position (1-based HGVS, inclusive) of the repeat region
        end: u64,
    },
}

/// Result of analyzing a potential repeat region
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatAnalysis {
    /// Whether this is a homopolymer (single-base repeat)
    pub is_homopolymer: bool,
    /// The repeated base (for homopolymer)
    pub base: Option<u8>,
    /// Start position (0-based) of the repeat in reference
    pub ref_start: usize,
    /// End position (0-based, exclusive) of the repeat in reference
    pub ref_end: usize,
    /// Count of repeats in reference
    pub ref_count: u64,
}

/// Analyze a position to find the extent of any homopolymer repeat
///
/// Given a position in the reference, finds the full extent of any
/// single-base repeat that includes this position.
///
/// Returns None if the position is not in a repeat (count < 2).
pub fn find_homopolymer_at(ref_seq: &[u8], pos: usize) -> Option<RepeatAnalysis> {
    if pos >= ref_seq.len() {
        return None;
    }

    let base = ref_seq[pos];

    // Find start of the repeat (scan left)
    let mut start = pos;
    while start > 0 && ref_seq[start - 1] == base {
        start -= 1;
    }

    // Find end of the repeat (scan right)
    let mut end = pos + 1;
    while end < ref_seq.len() && ref_seq[end] == base {
        end += 1;
    }

    let count = (end - start) as u64;

    // Only consider it a repeat if count >= 2
    if count < 2 {
        return None;
    }

    Some(RepeatAnalysis {
        is_homopolymer: true,
        base: Some(base),
        ref_start: start,
        ref_end: end,
        ref_count: count,
    })
}

/// Check if an insertion in a homopolymer should become repeat notation
///
/// When inserting bases that are identical to an existing homopolymer repeat,
/// the result should be expressed as repeat notation (e.g., A[9]).
///
/// Returns the total count if this should be repeat notation, None otherwise.
pub fn insertion_to_repeat(
    ref_seq: &[u8],
    pos: u64,
    inserted_seq: &[u8],
) -> Option<(u8, u64, u64, u64)> {
    // Only handle single-base insertions for now
    // (multi-base tandem repeats are more complex)
    if inserted_seq.is_empty() {
        return None;
    }

    // Check if all inserted bases are the same
    let first = inserted_seq[0];
    if !inserted_seq.iter().all(|&b| b == first) {
        return None;
    }

    let pos_idx = pos as usize;

    // Find the homopolymer at or near this insertion position
    // For an insertion between positions X and X+1, we need to find any adjacent
    // homopolymer of the same base and pick the largest one.
    // Check multiple positions and pick the best (largest) homopolymer
    let mut best_analysis: Option<RepeatAnalysis> = None;

    // Check position before insertion (X-1)
    if pos_idx > 0 && ref_seq.get(pos_idx - 1) == Some(&first) {
        if let Some(analysis) = find_homopolymer_at(ref_seq, pos_idx - 1) {
            if analysis.base == Some(first) {
                best_analysis = Some(analysis);
            }
        }
    }

    // Check position at insertion (X)
    if ref_seq.get(pos_idx) == Some(&first) {
        if let Some(analysis) = find_homopolymer_at(ref_seq, pos_idx) {
            if analysis.base == Some(first)
                && (best_analysis.is_none()
                    || analysis.ref_count > best_analysis.as_ref().unwrap().ref_count)
            {
                best_analysis = Some(analysis);
            }
        }
    }

    // Check position after insertion (X+1) - for cases where insertion extends a following tract
    if ref_seq.get(pos_idx + 1) == Some(&first) {
        if let Some(analysis) = find_homopolymer_at(ref_seq, pos_idx + 1) {
            if analysis.base == Some(first)
                && (best_analysis.is_none()
                    || analysis.ref_count > best_analysis.as_ref().unwrap().ref_count)
            {
                best_analysis = Some(analysis);
            }
        }
    }

    let analysis = best_analysis;

    if let Some(analysis) = analysis {
        if analysis.base == Some(first) {
            // Total count = reference count + inserted count
            let total_count = analysis.ref_count + inserted_seq.len() as u64;
            // Return (base, count, start, end) where start/end are 1-based
            // Per HGVS, positions refer to the reference repeat tract, not expanded
            return Some((
                first,
                total_count,
                index_to_hgvs_pos(analysis.ref_start),
                index_to_hgvs_pos(analysis.ref_start + analysis.ref_count as usize - 1),
            ));
        }
    }

    None
}

/// Get the complement of a DNA base
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'G' | b'g' => b'C',
        b'C' | b'c' => b'G',
        b'U' | b'u' => b'A', // RNA
        _ => base,           // N and other IUPAC codes
    }
}

/// Shorten an inversion by removing outer bases that cancel out
///
/// When the first base of an inversion is complementary to the last base,
/// inverting them produces no net change, so they can be excluded.
/// This is repeated until no more cancellation is possible.
///
/// Returns the shortened (start, end) positions (0-indexed), or None if
/// the inversion reduces to identity (0 or 1 base).
pub fn shorten_inversion(ref_seq: &[u8], start: usize, end: usize) -> Option<(usize, usize)> {
    if start >= end || end > ref_seq.len() {
        return None;
    }

    let mut s = start;
    let mut e = end;

    // Keep shrinking while outer bases are complementary
    while s < e {
        let first = ref_seq[s];
        let last = ref_seq[e - 1]; // e is exclusive

        // Check if first base is complement of last base
        if complement(first) == last {
            s += 1;
            e -= 1;
        } else {
            break;
        }
    }

    // If inversion shrinks to 0 or 1 base, it's identity
    if e <= s + 1 {
        return None; // Becomes identity - no inversion needed
    }

    Some((s, e))
}

/// Check if a delins should be represented as a duplication
///
/// In HGVS, if a delins deletes N bases and inserts 2N bases where the
/// inserted sequence is the deleted sequence repeated twice, the net effect
/// is a duplication.
///
/// Example: c.5delinsGG where position 5 has G → del G, ins GG = dup G
///
/// Arguments:
/// - ref_seq: The reference sequence (0-indexed)
/// - start: Start position (0-indexed, inclusive)
/// - end: End position (0-indexed, exclusive)
/// - inserted_seq: The sequence being inserted
///
/// Returns: true if this delins should become a duplication
pub fn delins_is_duplication(
    ref_seq: &[u8],
    start: usize,
    end: usize,
    inserted_seq: &[u8],
) -> bool {
    // Get the deleted region
    if start >= end || end > ref_seq.len() {
        return false;
    }

    let deleted_len = end - start;
    let deleted_seq = &ref_seq[start..end];

    // For delins to be a dup, inserted must be 2x the deleted length
    // and the inserted sequence must be the deleted sequence repeated twice
    if inserted_seq.len() != 2 * deleted_len {
        return false;
    }

    // Check that inserted = deleted + deleted
    let (first_half, second_half) = inserted_seq.split_at(deleted_len);
    first_half == deleted_seq && second_half == deleted_seq
}

/// Check if a duplication in a homopolymer should become repeat notation
///
/// When duplicating bases within a homopolymer repeat,
/// the result should be expressed as repeat notation.
/// Result of duplication to repeat conversion
#[derive(Debug, Clone)]
pub enum DupToRepeatResult {
    /// Single-base homopolymer repeat
    Homopolymer {
        base: u8,
        count: u64,
        start: u64, // 1-based (HGVS), use index_to_hgvs_pos() for conversion
        end: u64,   // 1-based (HGVS), inclusive
    },
    /// Multi-base tandem repeat
    TandemRepeat {
        unit: Vec<u8>,
        count: u64,
        start: u64, // 1-based (HGVS), use index_to_hgvs_pos() for conversion
        end: u64,   // 1-based (HGVS), inclusive
    },
}

pub fn duplication_to_repeat(ref_seq: &[u8], start: u64, end: u64) -> Option<DupToRepeatResult> {
    let start_idx = start as usize;
    let end_idx = end as usize;

    if start_idx >= ref_seq.len() || end_idx > ref_seq.len() || start_idx >= end_idx {
        return None;
    }

    let dup_seq = &ref_seq[start_idx..end_idx];
    if dup_seq.is_empty() {
        return None;
    }

    let dup_len = dup_seq.len();

    // Check if all duplicated bases are the same (homopolymer)
    // IMPORTANT: Only convert to repeat notation when duplicating 2+ bases.
    // Single-base duplications stay as simple dups (e.g., c.5266dup stays as dup, not C[4])
    let first = dup_seq[0];
    if dup_len >= 2 && dup_seq.iter().all(|&b| b == first) {
        // Find the full homopolymer extent
        if let Some(analysis) = find_homopolymer_at(ref_seq, start_idx) {
            if analysis.base == Some(first) {
                let total_count = analysis.ref_count + dup_len as u64;
                return Some(DupToRepeatResult::Homopolymer {
                    base: first,
                    count: total_count,
                    start: index_to_hgvs_pos(analysis.ref_start),
                    end: index_to_hgvs_pos(analysis.ref_start + analysis.ref_count as usize - 1),
                });
            }
        }
    }

    // Check for multi-base tandem repeat
    // The duplicated sequence could be multiple copies of a repeat unit
    // IMPORTANT: Only convert to repeat notation when duplicating MULTIPLE copies
    // of the repeat unit. Single-copy duplications stay as simple dups.
    // Try different unit lengths from 1 to half the dup length
    for unit_len in 1..=dup_len / 2 {
        if !dup_len.is_multiple_of(unit_len) {
            continue;
        }

        let unit = &dup_seq[0..unit_len];
        let copies_in_dup = dup_len / unit_len;

        // Only convert to repeat if duplicating 2+ copies of the unit
        if copies_in_dup < 2 {
            continue;
        }

        // Check if dup_seq is made of repeated copies of unit
        let is_repeat = (0..copies_in_dup).all(|i| {
            let chunk = &dup_seq[i * unit_len..(i + 1) * unit_len];
            chunk == unit
        });

        if !is_repeat {
            continue;
        }

        // Found a repeat unit. Now find the full extent in the reference.
        if let Some((ref_count, rep_start, rep_end)) =
            count_tandem_repeats(ref_seq, start_idx, unit)
        {
            // Total count = reference count + duplicated copies
            let total_count = ref_count + copies_in_dup as u64;
            // rep_end is exclusive (0-based), so last position is rep_end - 1
            return Some(DupToRepeatResult::TandemRepeat {
                unit: unit.to_vec(),
                count: total_count,
                start: index_to_hgvs_pos(rep_start),
                end: index_to_hgvs_pos(rep_end - 1),
            });
        }
    }

    // Note: We do NOT convert single-copy tandem repeat duplications to repeat notation
    // e.g., c.360_362dupGCA stays as dup, not GCA[8]

    None
}

/// Find tandem repeat at a position for a given repeat unit
///
/// Given a position and a repeat unit sequence, finds how many times that
/// sequence is repeated at that position in the reference.
///
/// Returns (count, start, end) where:
/// - count is the number of repeats found
/// - start is the 0-indexed start of the repeat region
/// - end is the 0-indexed exclusive end
pub fn count_tandem_repeats(
    ref_seq: &[u8],
    pos: usize,
    repeat_unit: &[u8],
) -> Option<(u64, usize, usize)> {
    if repeat_unit.is_empty() || pos >= ref_seq.len() {
        return None;
    }

    let unit_len = repeat_unit.len();

    // Check if the repeat unit matches at the position
    if pos + unit_len > ref_seq.len() {
        return None;
    }

    // Try to find repeats starting at or before this position
    // First, scan backwards to find the start of the repeat region
    let mut start = pos;
    while start >= unit_len {
        let candidate = &ref_seq[start - unit_len..start];
        if candidate == repeat_unit {
            start -= unit_len;
        } else {
            break;
        }
    }

    // Now count forward from the start
    let mut end = start;
    let mut count = 0u64;
    while end + unit_len <= ref_seq.len() {
        if &ref_seq[end..end + unit_len] == repeat_unit {
            count += 1;
            end += unit_len;
        } else {
            break;
        }
    }

    if count >= 1 {
        Some((count, start, end))
    } else {
        None
    }
}

/// Result of repeat normalization
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepeatNormResult {
    /// Convert to deletion (specified count < reference count)
    Deletion {
        start: u64, // 1-based (HGVS), inclusive start of region to delete
        end: u64,   // 1-based (HGVS), inclusive end of region to delete
    },
    /// Convert to duplication (specified count = reference count + 1)
    Duplication {
        start: u64, // 1-based (HGVS), inclusive start of duplicated region
        end: u64,   // 1-based (HGVS), inclusive end of duplicated region
        sequence: Vec<u8>,
    },
    /// Keep as repeat notation with canonical position
    Repeat {
        start: u64, // 1-based (HGVS), inclusive start
        end: u64,   // 1-based (HGVS), inclusive end
        sequence: Vec<u8>,
        count: u64,
    },
    /// No change needed
    Unchanged,
}

/// Normalize a repeat variant
///
/// Given a repeat notation like CAT[1], determines the appropriate
/// normalized representation by comparing to the reference.
///
/// Arguments:
/// - ref_seq: The reference sequence (0-indexed)
/// - pos: The position (0-indexed) where the repeat is specified
/// - repeat_unit: The repeated sequence (e.g., b"CAT")
/// - specified_count: The count specified in the variant
///
/// Returns the normalized representation
pub fn normalize_repeat(
    ref_seq: &[u8],
    pos: usize,
    repeat_unit: &[u8],
    specified_count: u64,
) -> RepeatNormResult {
    // Count how many times the repeat unit appears in the reference
    let Some((ref_count, ref_start, ref_end)) = count_tandem_repeats(ref_seq, pos, repeat_unit)
    else {
        return RepeatNormResult::Unchanged;
    };

    let unit_len = repeat_unit.len() as u64;

    if specified_count < ref_count {
        // Convert to deletion - we're removing (ref_count - specified_count) copies
        let del_count = ref_count - specified_count;
        // Delete from 3' end (HGVS convention)
        let del_len = del_count * unit_len;
        // ref_end is exclusive (0-based), so last position is ref_end - 1
        let del_end_idx = ref_end - 1;
        let del_start_idx = ref_end - del_len as usize;
        RepeatNormResult::Deletion {
            start: index_to_hgvs_pos(del_start_idx),
            end: index_to_hgvs_pos(del_end_idx),
        }
    } else if specified_count == ref_count + 1 {
        // Convert to duplication - we're adding exactly one copy
        // The duplicated region is the last copy in the reference
        // ref_end is exclusive, so last position is ref_end - 1
        let dup_end_idx = ref_end - 1;
        let dup_start_idx = ref_end - repeat_unit.len();
        RepeatNormResult::Duplication {
            start: index_to_hgvs_pos(dup_start_idx),
            end: index_to_hgvs_pos(dup_end_idx),
            sequence: repeat_unit.to_vec(),
        }
    } else if specified_count == ref_count {
        // Same as reference - this is identity (no change)
        RepeatNormResult::Unchanged
    } else {
        // Keep as repeat notation with canonical position
        // The repeat region describes the REFERENCE tract (per HGVS spec)
        // ref_end is exclusive, so last position is ref_end - 1
        RepeatNormResult::Repeat {
            start: index_to_hgvs_pos(ref_start),
            end: index_to_hgvs_pos(ref_end - 1),
            sequence: repeat_unit.to_vec(),
            count: specified_count,
        }
    }
}

/// Get the canonical form for an edit
///
/// HGVS rules:
/// - Deletions ALWAYS stay as deletions (just shift 3')
/// - Insertions become duplications if they match adjacent sequence
/// - Duplications stay as duplications
/// - Delins stays as delins
pub fn get_canonical_form(edit: &NaEdit, ref_seq: &[u8], start: u64, _end: u64) -> CanonicalForm {
    use crate::hgvs::edit::InsertedSequence;

    match edit {
        NaEdit::Deletion { .. } => {
            // Deletions always stay as deletions - never convert to dup
            CanonicalForm::Deletion
        }
        NaEdit::Insertion { sequence } => {
            // Insertions may become duplications if they match adjacent sequence
            // Only check Literal sequences (others like Count, Range don't have actual bases)
            if let InsertedSequence::Literal(seq) = sequence {
                // Convert Sequence (Vec<Base>) to bytes
                let seq_bytes: Vec<u8> = seq.bases().iter().map(|b| *b as u8).collect();
                if insertion_is_duplication(ref_seq, start, &seq_bytes) {
                    return CanonicalForm::Duplication;
                }
            }
            CanonicalForm::Insertion
        }
        NaEdit::Delins { .. } => CanonicalForm::Delins,
        NaEdit::Duplication { .. } => CanonicalForm::Duplication,
        _ => CanonicalForm::Deletion, // Default fallback
    }
}

/// Apply minimal notation rules to an edit without requiring reference data.
///
/// This function applies HGVS minimal notation rules:
/// - Deletions: remove explicit sequence and length (del12 → del, delATG → del)
/// - Delins: remove explicit deleted sequence from delins notation
///   (delATGinsGGG → delinsGGG)
/// - Duplications: remove explicit sequence and length (dupATG → dup, dup12 → dup)
///
/// This is useful for canonicalizing variants even when reference sequence
/// is not available for full 3'/5' normalization.
pub fn canonicalize_edit(edit: &NaEdit) -> NaEdit {
    match edit {
        NaEdit::Deletion { .. } => NaEdit::Deletion {
            sequence: None,
            length: None,
        },
        NaEdit::Duplication {
            uncertain_extent, ..
        } => NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: uncertain_extent.clone(),
        },
        NaEdit::Delins { sequence } => {
            // Keep the inserted sequence but remove any explicit deletion info
            // that might be embedded in the sequence description
            NaEdit::Delins {
                sequence: sequence.clone(),
            }
        }
        // Other edits pass through unchanged
        _ => edit.clone(),
    }
}

/// Apply minimal notation to a variant without reference data.
/// Returns true if the edit was modified.
pub fn should_canonicalize(edit: &NaEdit) -> bool {
    match edit {
        NaEdit::Deletion { sequence, length } => sequence.is_some() || length.is_some(),
        NaEdit::Duplication {
            sequence, length, ..
        } => sequence.is_some() || length.is_some(),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_needs_normalization() {
        assert!(needs_normalization(&NaEdit::Deletion {
            sequence: None,
            length: None,
        }));
        assert!(needs_normalization(&NaEdit::Duplication {
            sequence: None,
            length: None,
            uncertain_extent: None,
        }));
        assert!(!needs_normalization(&NaEdit::Inversion {
            sequence: None,
            length: None,
        }));
    }

    #[test]
    fn test_insertion_is_duplication() {
        // Sequence: ATGATGATG
        // Inserting ATG at position 3 (after first ATG) is a dup
        let ref_seq = b"ATGATGATG";
        assert!(insertion_is_duplication(ref_seq, 3, b"ATG"));

        // Inserting TGA at position 3 is not a dup
        assert!(!insertion_is_duplication(ref_seq, 3, b"TGA"));

        // Inserting ATG at position 6 (after second ATG) is a dup
        assert!(insertion_is_duplication(ref_seq, 6, b"ATG"));
    }

    #[test]
    fn test_deletion_stays_deletion() {
        // Deletions should ALWAYS stay as deletions, never become dups
        let ref_seq = b"ATGATGATG";
        let del_edit = NaEdit::Deletion {
            sequence: None,
            length: None,
        };
        // Even if the deleted sequence matches preceding sequence,
        // it should remain a deletion
        assert_eq!(
            get_canonical_form(&del_edit, ref_seq, 3, 6),
            CanonicalForm::Deletion
        );
    }

    #[test]
    fn test_insertion_becomes_dup() {
        use crate::hgvs::edit::{InsertedSequence, Sequence};
        use std::str::FromStr;

        let ref_seq = b"ATGATGATG";
        let ins_edit = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("ATG").unwrap()),
        };
        // Inserting ATG after first ATG (pos 3) should become dup
        assert_eq!(
            get_canonical_form(&ins_edit, ref_seq, 3, 3),
            CanonicalForm::Duplication
        );

        // Inserting TGA should stay as insertion
        let ins_edit2 = NaEdit::Insertion {
            sequence: InsertedSequence::Literal(Sequence::from_str("TGA").unwrap()),
        };
        assert_eq!(
            get_canonical_form(&ins_edit2, ref_seq, 3, 3),
            CanonicalForm::Insertion
        );
    }

    // =========================================================================
    // HOMOPOLYMER REPEAT TESTS
    // =========================================================================

    #[test]
    fn test_find_homopolymer_at() {
        // Sequence: GGGAAAAAGGG (4 A's at positions 3-6, 0-indexed)
        let ref_seq = b"GGGAAAAAGGG";

        // Find homopolymer at position 4 (middle of A's)
        let result = find_homopolymer_at(ref_seq, 4);
        assert!(result.is_some());
        let analysis = result.unwrap();
        assert!(analysis.is_homopolymer);
        assert_eq!(analysis.base, Some(b'A'));
        assert_eq!(analysis.ref_start, 3);
        assert_eq!(analysis.ref_end, 8); // exclusive
        assert_eq!(analysis.ref_count, 5);

        // Position in G's at start
        let result = find_homopolymer_at(ref_seq, 1);
        assert!(result.is_some());
        let analysis = result.unwrap();
        assert_eq!(analysis.base, Some(b'G'));
        assert_eq!(analysis.ref_count, 3);

        // Single base (no repeat) - should return None
        let single_seq = b"ATGC";
        assert!(find_homopolymer_at(single_seq, 0).is_none());
    }

    #[test]
    fn test_insertion_to_repeat() {
        // Sequence: GGGAAAAAGGG (5 A's at positions 3-7, 0-indexed)
        let ref_seq = b"GGGAAAAAGGG";

        // Inserting AA at position 8 (after the A's) should become A[7]
        // Position 8 is 0-indexed, which is after the last A
        let result = insertion_to_repeat(ref_seq, 8, b"AA");
        assert!(result.is_some());
        let (base, count, start, end) = result.unwrap();
        assert_eq!(base, b'A');
        assert_eq!(count, 7); // 5 original + 2 inserted
        assert_eq!(start, 4); // 1-indexed start of A region in reference
        assert_eq!(end, 8); // 1-indexed end of A region in reference (per HGVS, positions refer to reference tract)

        // Inserting T (non-matching) should return None
        let result = insertion_to_repeat(ref_seq, 8, b"T");
        assert!(result.is_none());

        // Inserting mixed bases should return None
        let result = insertion_to_repeat(ref_seq, 8, b"AT");
        assert!(result.is_none());
    }

    #[test]
    fn test_duplication_to_repeat() {
        // Sequence: GGGAAAAAGGG (5 A's at positions 3-7, 0-indexed)
        let ref_seq = b"GGGAAAAAGGG";

        // Duplicating 2 A's (positions 3-5, 0-indexed) should become A[7]
        let result = duplication_to_repeat(ref_seq, 3, 5);
        assert!(result.is_some());
        match result.unwrap() {
            DupToRepeatResult::Homopolymer {
                base, count, start, ..
            } => {
                assert_eq!(base, b'A');
                assert_eq!(count, 7); // 5 original + 2 duplicated
                assert_eq!(start, 4); // 1-indexed start
            }
            _ => panic!("Expected Homopolymer result"),
        }

        // Duplicating non-homopolymer region should return None (if not a tandem repeat)
        // ATGCXYZ has no repeats, so duplicating ATG should return None
        let non_repeat_seq = b"ATGCXYZ";
        let result = duplication_to_repeat(non_repeat_seq, 0, 3);
        assert!(result.is_none());
    }

    #[test]
    fn test_duplication_to_tandem_repeat() {
        // Sequence with GCA repeats
        // String: AAAAAGCAGCAGCAGCAGCAGCAGCAGCAAAAA (33 chars)
        // 5 A's + 8 GCAs (24 chars) + 4 A's = 33 chars
        // 8 GCA repeats at 0-indexed positions 5-28 (exclusive 29)
        let ref_seq = b"AAAAAGCAGCAGCAGCAGCAGCAGCAGCAAAAA";

        // Duplicating single GCA (one copy) should NOT become repeat notation
        // It should stay as a simple dup per HGVS rules
        let result = duplication_to_repeat(ref_seq, 5, 8);
        assert!(result.is_none(), "Single-copy dup should not become repeat");

        // Duplicating GCAGCA (2 copies) SHOULD become repeat notation
        let result = duplication_to_repeat(ref_seq, 5, 11);
        assert!(result.is_some());
        match result.unwrap() {
            DupToRepeatResult::TandemRepeat {
                unit,
                count,
                start,
                end,
            } => {
                assert_eq!(unit, b"GCA");
                assert_eq!(count, 10); // 8 original + 2 duplicated
                assert_eq!(start, 6); // 1-indexed
                assert_eq!(end, 29); // 1-indexed end of 8 GCAs (0-indexed 28 + 1)
            }
            _ => panic!("Expected TandemRepeat result"),
        }
    }

    // =========================================================================
    // TANDEM REPEAT COUNTING TESTS
    // =========================================================================

    #[test]
    fn test_count_tandem_repeats_basic() {
        // Sequence: GGGCATCATCATGGG (3 CAT repeats at positions 3-11)
        let ref_seq = b"GGGCATCATCATGGG";

        let result = count_tandem_repeats(ref_seq, 3, b"CAT");
        assert!(result.is_some());
        let (count, start, end) = result.unwrap();
        assert_eq!(count, 3);
        assert_eq!(start, 3);
        assert_eq!(end, 12); // exclusive

        // Try from middle of the repeat
        let result = count_tandem_repeats(ref_seq, 6, b"CAT");
        assert!(result.is_some());
        let (count, start, end) = result.unwrap();
        assert_eq!(count, 3);
        assert_eq!(start, 3);
        assert_eq!(end, 12);
    }

    #[test]
    fn test_count_tandem_repeats_single_base() {
        // Sequence: GGGAAAAAAGGG (6 A's)
        let ref_seq = b"GGGAAAAAAGGG";

        let result = count_tandem_repeats(ref_seq, 5, b"A");
        assert!(result.is_some());
        let (count, start, end) = result.unwrap();
        assert_eq!(count, 6);
        assert_eq!(start, 3);
        assert_eq!(end, 9);
    }

    #[test]
    fn test_count_tandem_repeats_no_match() {
        let ref_seq = b"GGGAAAAAAGGG";
        // XYZ doesn't appear in the sequence
        let result = count_tandem_repeats(ref_seq, 5, b"XYZ");
        assert!(result.is_none());
    }

    // =========================================================================
    // REPEAT NORMALIZATION TESTS
    // =========================================================================

    #[test]
    fn test_normalize_repeat_to_deletion() {
        // Sequence: GGGCATCATCATCATGGG (4 CAT repeats at positions 3-14)
        // Specifying CAT[1] should become deletion of 3 CATs (9 bases)
        let ref_seq = b"GGGCATCATCATCATGGG";

        let result = normalize_repeat(ref_seq, 3, b"CAT", 1);
        match result {
            RepeatNormResult::Deletion { start, end } => {
                // Should delete positions 7-15 (1-indexed), which is 3 CATs
                assert_eq!(end - start + 1, 9, "Should delete 9 bases (3 CATs)");
            }
            _ => panic!("Expected Deletion, got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_to_duplication() {
        // Sequence: GGGCATCATGGG (2 CAT repeats at positions 3-8)
        // Specifying CAT[3] (ref is 2, so 2+1=3) should become duplication
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(ref_seq, 3, b"CAT", 3);
        match result {
            RepeatNormResult::Duplication {
                start,
                end,
                sequence,
            } => {
                assert_eq!(sequence, b"CAT");
                assert_eq!(end - start + 1, 3, "Should duplicate 3 bases (1 CAT)");
            }
            _ => panic!("Expected Duplication, got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_stays_repeat() {
        // Sequence: GGGCATCATGGG (2 CAT repeats)
        // Specifying CAT[5] (ref is 2, 5 > 2+1) should stay as repeat
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(ref_seq, 3, b"CAT", 5);
        match result {
            RepeatNormResult::Repeat {
                count, sequence, ..
            } => {
                assert_eq!(sequence, b"CAT");
                assert_eq!(count, 5);
            }
            _ => panic!("Expected Repeat, got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_unchanged() {
        // Sequence: GGGCATCATGGG (2 CAT repeats)
        // Specifying CAT[2] (same as ref) should be unchanged
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(ref_seq, 3, b"CAT", 2);
        assert!(matches!(result, RepeatNormResult::Unchanged));
    }

    // =========================================================================
    // INVERSION SHORTENING TESTS
    // =========================================================================

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'N'), b'N'); // N stays N
    }

    #[test]
    fn test_shorten_inversion_basic() {
        // Sequence: ATGCAT (positions 0-5)
        // A(0) is complement of T(5) - cancel
        // T(1) is complement of A(4) - cancel
        // G(2) is complement of C(3) - cancel
        // All bases cancel -> becomes identity
        let seq = b"ATGCAT";
        let result = shorten_inversion(seq, 0, 6);
        assert!(
            result.is_none(),
            "Fully complementary inversion should become identity"
        );
    }

    #[test]
    fn test_shorten_inversion_partial() {
        // Sequence: ATGGAT (positions 0-5)
        // A(0) is complement of T(5) - cancel
        // T(1) is complement of A(4) - cancel
        // G(2) is NOT complement of G(3) - stop
        // Result: positions 2-4
        let seq = b"ATGGAT";
        let result = shorten_inversion(seq, 0, 6);
        assert!(result.is_some());
        let (s, e) = result.unwrap();
        assert_eq!(s, 2);
        assert_eq!(e, 4);
    }

    #[test]
    fn test_shorten_inversion_no_change() {
        // Sequence: GGCC (positions 0-3)
        // G(0) is complement of C(3) - cancel
        // G(1) is complement of C(2) - cancel
        // All cancel -> identity
        let seq = b"GGCC";
        let result = shorten_inversion(seq, 0, 4);
        assert!(result.is_none());

        // Sequence: GATT (positions 0-3)
        // G(0) is NOT complement of T(3) - no cancellation
        let seq2 = b"GATT";
        let result2 = shorten_inversion(seq2, 0, 4);
        assert!(result2.is_some());
        let (s, e) = result2.unwrap();
        assert_eq!(s, 0);
        assert_eq!(e, 4);
    }

    // =========================================================================
    // CANONICALIZATION TESTS
    // =========================================================================

    #[test]
    fn test_canonicalize_deletion_with_length() {
        use crate::hgvs::edit::Sequence;
        use std::str::FromStr;

        // del12 -> del (remove explicit length)
        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(12),
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Deletion {
                sequence: None,
                length: None
            }
        ));

        // delATG -> del (remove explicit sequence)
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Deletion {
                sequence: None,
                length: None
            }
        ));
    }

    #[test]
    fn test_canonicalize_duplication_with_length() {
        use crate::hgvs::edit::Sequence;
        use std::str::FromStr;

        // dup12 -> dup
        let edit = NaEdit::Duplication {
            sequence: None,
            length: Some(12),
            uncertain_extent: None,
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Duplication {
                sequence: None,
                length: None,
                ..
            }
        ));

        // dupATG -> dup
        let edit = NaEdit::Duplication {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
            uncertain_extent: None,
        };
        let canonical = canonicalize_edit(&edit);
        assert!(matches!(
            canonical,
            NaEdit::Duplication {
                sequence: None,
                length: None,
                ..
            }
        ));
    }

    #[test]
    fn test_should_canonicalize() {
        // Deletion with length should be canonicalized
        assert!(should_canonicalize(&NaEdit::Deletion {
            sequence: None,
            length: Some(12)
        }));

        // Deletion without length shouldn't be canonicalized
        assert!(!should_canonicalize(&NaEdit::Deletion {
            sequence: None,
            length: None
        }));

        // Substitution shouldn't be canonicalized
        use crate::hgvs::edit::Base;
        assert!(!should_canonicalize(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G
        }));
    }
}
