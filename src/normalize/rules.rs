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
    if matches!(
        edit,
        NaEdit::Deletion { .. }
            | NaEdit::Insertion { .. }
            | NaEdit::Duplication { .. }
            | NaEdit::Delins { .. }
            | NaEdit::Repeat { .. }
    ) {
        return true;
    }
    // A substitution with ref == alt is degenerate and must be rewritten to
    // identity (`=`); route it through normalization so the rule fires.
    // Real substitutions (ref != alt) keep the fast no-op path.
    matches!(
        edit,
        NaEdit::Substitution { reference, alternative } if reference == alternative
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

    // Guard against positions beyond the reference sequence
    if ref_seq.is_empty() || pos_idx > ref_seq.len() {
        return false;
    }

    // Check if sequence before the insertion point matches the inserted sequence
    // (this would be a duplication of the preceding sequence)
    if pos_idx >= ins_len && pos_idx <= ref_seq.len() {
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

/// Check if an insertion in a tandem repeat should become repeat notation.
///
/// Per HGVS spec: when an insertion adds 2 or more copies of a tandem
/// repeat unit already present adjacent to the insertion point, the
/// variant is expressed as repeat notation `unit[N+k]`. A single-copy
/// addition is left as a duplication (handled by `insertion_is_duplication`).
///
/// `is_coding` enables the spec's codon-frame exception (repeated.md): in
/// c. context, repeat notation requires `unit_len % 3 == 0`. When the gate
/// blocks the rewrite this returns `None` and the caller falls back to a
/// literal `ins` description.
///
/// Returns `Some((first_byte, total_count, ref_start_1based,
/// ref_end_1based, unit_bytes))` when repeat notation applies, `None`
/// otherwise. `first_byte` is preserved for backwards-compat with the
/// homopolymer caller; `unit_bytes` is the full tandem unit.
pub fn insertion_to_repeat(
    ref_seq: &[u8],
    pos: u64,
    inserted_seq: &[u8],
    is_coding: bool,
) -> Option<(u8, u64, u64, u64, Vec<u8>)> {
    if inserted_seq.is_empty() {
        return None;
    }

    let base_unit = smallest_repeat_unit(inserted_seq);
    let added_copies = (inserted_seq.len() / base_unit.len()) as u64;
    if added_copies < 2 {
        // Single-copy addition is a duplication, not repeat notation.
        return None;
    }

    // HGVS spec (repeated.md): in c. context, repeat notation requires
    // unit_len % 3 == 0. Falls back to literal `ins` when the gate blocks.
    if is_coding && !base_unit.len().is_multiple_of(3) {
        return None;
    }

    // The inserted sequence may be written as any cyclic rotation of the
    // reference repeat unit (e.g. `insCAGCAG` against a `GCA` tract). Try
    // every rotation of `base_unit` and pick the one that yields the
    // largest contiguous tandem run anchored at the insertion point.
    let u_len = base_unit.len();
    let mut best: Option<(Vec<u8>, usize, u64)> = None;
    for r in 0..u_len {
        let mut rotated = Vec::with_capacity(u_len);
        rotated.extend_from_slice(&base_unit[r..]);
        rotated.extend_from_slice(&base_unit[..r]);
        if let Some((ref_start, ref_count)) = find_tandem_extent(ref_seq, pos as usize, &rotated) {
            if ref_count > 0 && best.as_ref().is_none_or(|(_, _, bc)| ref_count > *bc) {
                best = Some((rotated, ref_start, ref_count));
            }
        }
    }

    let (unit, ref_start, ref_count) = best?;
    let total_count = ref_count + added_copies;
    let ref_end = ref_start + ref_count as usize * unit.len() - 1;
    Some((
        unit[0],
        total_count,
        index_to_hgvs_pos(ref_start),
        index_to_hgvs_pos(ref_end),
        unit,
    ))
}

/// Smallest unit `U` such that `seq = U * k` for some integer k ≥ 1.
pub(crate) fn smallest_repeat_unit(seq: &[u8]) -> &[u8] {
    let n = seq.len();
    for u in 1..=n {
        if !n.is_multiple_of(u) {
            continue;
        }
        let unit = &seq[..u];
        if seq.chunks_exact(u).all(|c| c == unit) {
            return unit;
        }
    }
    seq // fallback (unreachable for n>0)
}

/// Find the maximal tandem run of `unit` near position `pos` in
/// `ref_seq`. `pos` is the 0-based index of the base immediately 5' of
/// the insertion point — i.e., the insertion is between `pos` and
/// `pos+1`. The tract may not be unit-aligned with `pos`, so we probe
/// anchor offsets in a window around `pos` and pick the run that
/// abuts/contains the insertion point.
///
/// Returns `(ref_start_0based, count)` of the chosen tandem run.
fn find_tandem_extent(ref_seq: &[u8], pos: usize, unit: &[u8]) -> Option<(usize, u64)> {
    let u = unit.len();
    if u == 0 {
        return None;
    }

    // Insertion point in "between-bases" coords: between pos and pos+1.
    let ins_point = pos + 1;

    // Probe candidate anchor offsets within [ins_point - u, ins_point].
    // For each anchor `a`, compute the maximal tandem run that includes
    // anchor `a` (walk left and right in steps of `u`). Accept runs that
    // abut or span the insertion point. Among accepted runs, pick the
    // largest count.
    let lo = ins_point.saturating_sub(u);
    let hi = ins_point.min(ref_seq.len());

    let mut best: Option<(usize, u64)> = None;

    for anchor in lo..=hi {
        if anchor + u > ref_seq.len() {
            continue;
        }
        if &ref_seq[anchor..anchor + u] != unit {
            continue;
        }

        // Anchor identifies one unit-aligned occurrence; extend to the full tract.
        // The `?`-via-let-else is defensive: we just verified the anchor matches,
        // so `extend_tandem_tract` shouldn't return None here.
        let Some(TandemTract {
            start,
            end,
            ref_count: count,
        }) = extend_tandem_tract(ref_seq, anchor..anchor + u, unit)
        else {
            continue;
        };

        // Require the run [start, end) to abut/span ins_point: the
        // insertion point must lie within [start, end].
        if ins_point < start || ins_point > end {
            continue;
        }

        match best {
            None => best = Some((start, count)),
            Some((_, bc)) if count > bc => best = Some((start, count)),
            _ => {}
        }
    }

    best
}

/// Information about a maximal tandem repeat tract in a reference sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TandemTract {
    /// 0-based, inclusive start of the tract.
    pub start: usize,
    /// 0-based, exclusive end of the tract.
    pub end: usize,
    /// Number of full `unit`-length copies in the tract.
    pub ref_count: u64,
}

/// Walk left from `anchor.start` and right from `anchor.end` extending
/// unit-by-unit while `ref_seq` matches `unit`. The anchor span must lie on
/// `unit`-length boundaries — i.e., `anchor.len() % unit.len() == 0` — and
/// `ref_seq[anchor]` must be `unit` repeated some non-negative number of times.
///
/// Returns the maximal tract enclosing the anchor, or `None` if the anchor
/// itself isn't unit-periodic (defensive — `deletion_to_repeat` already checks
/// before calling, but we don't trust future callers).
pub(crate) fn extend_tandem_tract(
    ref_seq: &[u8],
    anchor: std::ops::Range<usize>,
    unit: &[u8],
) -> Option<TandemTract> {
    let u = unit.len();
    if u == 0 || anchor.start > anchor.end || anchor.end > ref_seq.len() {
        return None;
    }
    if !(anchor.end - anchor.start).is_multiple_of(u) {
        return None;
    }
    // Defensive: verify the anchor itself is unit-periodic.
    if !ref_seq[anchor.start..anchor.end]
        .chunks_exact(u)
        .all(|chunk| chunk == unit)
    {
        return None;
    }

    // Walk left in steps of u.
    let mut start = anchor.start;
    while start >= u && &ref_seq[start - u..start] == unit {
        start -= u;
    }
    // Walk right in steps of u.
    let mut end = anchor.end;
    while end + u <= ref_seq.len() && &ref_seq[end..end + u] == unit {
        end += u;
    }

    let ref_count = ((end - start) / u) as u64;
    Some(TandemTract {
        start,
        end,
        ref_count,
    })
}

/// Description of how a deletion can be re-expressed as repeat notation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DelToRepeatResult {
    /// The phase-aligned tandem unit (length 1 for homopolymer, longer for
    /// multi-base tandem).
    pub unit: Vec<u8>,
    /// Post-deletion unit count (`N - k`).
    pub count: u64,
    /// 1-based HGVS start of the canonical reference tract (inclusive).
    pub start: u64,
    /// 1-based HGVS end of the canonical reference tract (inclusive).
    pub end: u64,
}

/// Compute the canonical `unit[N-k]` description for a post-3'-shift deletion,
/// or return `None` if the deletion does not meet B2 trigger conditions.
///
/// # Shuffle phase-alignment lemma (why no rotation iteration)
///
/// The 3' shuffle (`shuffle.rs`) extends `del_end` while
/// `ref[del_start] == ref[del_end]`, terminating exactly when the period
/// breaks. For a deletion sitting inside a tandem tract with period `p`, the
/// post-shift `del_slice = ref_seq[del_start..del_end]` is therefore
/// unit-aligned with the surrounding tract, and `smallest_repeat_unit(del_slice)`
/// equals the tract's natural unit. So `extend_tandem_tract` walking from the
/// post-shift `del_start` finds the full tract on the r=0 phase — no rotation
/// iteration is required here.
///
/// # Trigger conditions
///
/// All five must hold; otherwise returns `None` (caller emits plain `del`):
///
/// 1. `del_len % unit_len == 0` (length is a multiple of the unit's period).
/// 2. The reference tract enclosing the deletion has `ref_count >= 2`.
/// 3. `k = del_len / unit_len >= 2` (single-unit removal stays as `del`).
/// 4. `post_count = ref_count - k >= 1` (full tract removal stays as `del`).
/// 5. The deletion lies entirely within the tract (handled by `extend_tandem_tract`'s
///    defensive periodicity check on the anchor).
///
/// The bounds check rejects 0-length deletions (`del_start == del_end`).
/// Callers passing insertion-point-shaped zero-width ranges should use
/// `extend_tandem_tract` directly instead.
pub(crate) fn deletion_to_repeat(
    ref_seq: &[u8],
    del_start: usize,
    del_end: usize,
    is_coding: bool,
) -> Option<DelToRepeatResult> {
    if del_start >= del_end || del_end > ref_seq.len() {
        return None;
    }
    let del_slice = &ref_seq[del_start..del_end];
    let unit_slice = smallest_repeat_unit(del_slice);
    let p = unit_slice.len();
    if p == 0 || !(del_end - del_start).is_multiple_of(p) {
        return None;
    }

    // HGVS spec (repeated.md): in c. context, repeat notation requires
    // unit_len % 3 == 0. Falls back to plain del when the gate blocks.
    if is_coding && !p.is_multiple_of(3) {
        return None;
    }

    let tract = extend_tandem_tract(ref_seq, del_start..del_end, unit_slice)?;
    if tract.ref_count < 2 {
        return None;
    }

    let k = ((del_end - del_start) / p) as u64;
    if k < 2 {
        return None;
    }

    let post_count = tract.ref_count - k;
    if post_count == 0 {
        return None;
    }

    Some(DelToRepeatResult {
        unit: unit_slice.to_vec(),
        count: post_count,
        start: index_to_hgvs_pos(tract.start),
        end: index_to_hgvs_pos(tract.end - 1),
    })
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

/// Check if a delins should be represented as identity
///
/// Per the HGVS spec, a delins whose inserted sequence equals the deleted
/// reference produces no change and must be expressed using identity
/// notation (`=`). The rule applies for any matching length.
///
/// Examples:
/// - g.1000delinsG where ref[1000] = G → g.1000=
/// - g.1000_1002delinsGCA where ref[1000..1002] = GCA → g.1000_1002=
pub fn delins_is_identity(ref_seq: &[u8], start: usize, end: usize, inserted_seq: &[u8]) -> bool {
    if start >= end || end > ref_seq.len() {
        return false;
    }
    let deleted_len = end - start;
    if inserted_seq.len() != deleted_len {
        return false;
    }
    &ref_seq[start..end] == inserted_seq
}

/// Check if a delins should be represented as a substitution
///
/// Per the HGVS edit-type priority (substitution > deletion > inversion >
/// duplication > insertion), a delins that replaces a single base with a
/// *different* single base must be expressed as a substitution.
///
/// Example: g.1000delinsA where ref[1000] is G → g.1000G>A
///
/// Returns `None` for the same-base case (e.g. `g.1000delinsG` where ref=G),
/// which is identity rather than substitution; that rewrite is handled by
/// `delins_is_identity`.
pub fn delins_is_substitution(
    ref_seq: &[u8],
    start: usize,
    end: usize,
    inserted_seq: &[u8],
) -> Option<(crate::hgvs::edit::Base, crate::hgvs::edit::Base)> {
    use crate::hgvs::edit::Base;

    if end != start + 1 || end > ref_seq.len() || inserted_seq.len() != 1 {
        return None;
    }
    let ref_byte = ref_seq[start];
    let alt_byte = inserted_seq[0];
    if ref_byte == alt_byte {
        return None;
    }
    Some((
        Base::from_char(ref_byte as char)?,
        Base::from_char(alt_byte as char)?,
    ))
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
    /// Codon-frame gate triggered: structural conditions for repeat
    /// notation are met, but `unit_len % 3 != 0` in c. context, so the
    /// canonical form is a literal `ins` of the duplicated sequence at
    /// the 3' flanking position. `start` and `end` are the two flanking
    /// 1-based HGVS positions (last tract base and the next base).
    GatedInsertion {
        start: u64,
        end: u64,
        sequence: Vec<u8>,
    },
}

pub fn duplication_to_repeat(
    ref_seq: &[u8],
    start: u64,
    end: u64,
    is_coding: bool,
) -> Option<DupToRepeatResult> {
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
                // Codon-frame gate (repeated.md): homopolymer unit_len=1 is
                // never codon-aligned, so c. context blocks repeat notation.
                // Emit GatedInsertion so the caller renders as `ins<dup_seq>`
                // at the 3' flanking position of the reference tract.
                if is_coding {
                    let last_tract_idx = analysis.ref_start + analysis.ref_count as usize - 1;
                    return Some(DupToRepeatResult::GatedInsertion {
                        start: index_to_hgvs_pos(last_tract_idx),
                        end: index_to_hgvs_pos(last_tract_idx) + 1,
                        sequence: dup_seq.to_vec(),
                    });
                }
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
            // Codon-frame gate (repeated.md): in c., repeat notation requires
            // unit_len % 3 == 0. Structural conditions are met but the gate
            // forces a literal `ins<dup_seq>` at the 3' tract flanking
            // position instead.
            if is_coding && !unit_len.is_multiple_of(3) {
                let last_tract_idx = rep_end - 1;
                return Some(DupToRepeatResult::GatedInsertion {
                    start: index_to_hgvs_pos(last_tract_idx),
                    end: index_to_hgvs_pos(last_tract_idx) + 1,
                    sequence: dup_seq.to_vec(),
                });
            }
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
    /// Convert to insertion. Used when the codon-frame gate blocks repeat
    /// notation in c. context for an expansion of >=2 unit copies. `start`
    /// and `end` are the two flanking 1-based HGVS positions; `sequence` is
    /// the literal inserted sequence (canonical unit repeated `specified -
    /// ref_count` times).
    Insertion {
        start: u64,
        end: u64,
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
    is_coding: bool,
) -> RepeatNormResult {
    // Match `count_tandem_repeats`'s pre-refactor contract: an empty unit
    // is meaningless and falls through to `Unchanged`. Without this guard,
    // `smallest_repeat_unit(b"")` returns `b""` and the period division below
    // panics.
    if repeat_unit.is_empty() {
        return RepeatNormResult::Unchanged;
    }

    // Canonicalize to the smallest period so callers spelling a non-minimal
    // unit (e.g. `ATAT[1]` over an `AT[4]` tract) reach the right branch.
    // Without this, a contraction misses B2 because k is computed in the
    // caller-spelled unit (1 ATAT removed) instead of the canonical unit
    // (2 ATs removed).
    let canonical_unit = smallest_repeat_unit(repeat_unit);
    let copies_per_input_unit = (repeat_unit.len() / canonical_unit.len()) as u64;
    let specified_count = specified_count * copies_per_input_unit;

    // Count how many times the canonical unit appears in the reference
    let Some((ref_count, ref_start, ref_end)) = count_tandem_repeats(ref_seq, pos, canonical_unit)
    else {
        return RepeatNormResult::Unchanged;
    };

    let unit_len = canonical_unit.len() as u64;
    // HGVS spec (repeated.md): in c. context, repeat notation requires
    // unit_len % 3 == 0. When the gate blocks, contraction-with-survivors
    // routes to Deletion and expansion-of->=2-copies routes to Insertion.
    let codon_blocks_repeat = is_coding && !canonical_unit.len().is_multiple_of(3);

    if specified_count < ref_count {
        let k = ref_count - specified_count;
        if k >= 2 && specified_count >= 1 && !codon_blocks_repeat {
            // B2 (symmetric with A7): >=2 unit reduction with surviving units → repeat
            RepeatNormResult::Repeat {
                start: index_to_hgvs_pos(ref_start),
                end: index_to_hgvs_pos(ref_end - 1),
                sequence: canonical_unit.to_vec(),
                count: specified_count,
            }
        } else {
            // 1-unit reduction, full tract removal, or codon-frame-gated → deletion
            // (HGVS prioritization: deletion outranks unranked repeat[0])
            let del_len = (k as usize) * unit_len as usize;
            let del_end_idx = ref_end - 1;
            let del_start_idx = ref_end - del_len;
            RepeatNormResult::Deletion {
                start: index_to_hgvs_pos(del_start_idx),
                end: index_to_hgvs_pos(del_end_idx),
            }
        }
    } else if specified_count == ref_count + 1 {
        // Convert to duplication - we're adding exactly one copy.
        // dup is always permitted; the spec exception only forbids `[N]`.
        // The duplicated region is the last copy in the reference.
        // ref_end is exclusive, so last position is ref_end - 1
        let dup_end_idx = ref_end - 1;
        let dup_start_idx = ref_end - canonical_unit.len();
        RepeatNormResult::Duplication {
            start: index_to_hgvs_pos(dup_start_idx),
            end: index_to_hgvs_pos(dup_end_idx),
            sequence: canonical_unit.to_vec(),
        }
    } else if specified_count == ref_count {
        // Same as reference - this is identity (no change)
        RepeatNormResult::Unchanged
    } else if codon_blocks_repeat {
        // Expansion of >=2 copies in c. with non-codon-aligned unit: spec
        // mandates `ins<literal>` form (e.g., c.1741_1742insTATATATA), not
        // `[N]`. Insertion point is between the last base of the reference
        // tract (ref_end - 1, 0-based) and the next base (ref_end).
        let added_copies = specified_count - ref_count;
        let mut inserted = Vec::with_capacity((added_copies as usize) * canonical_unit.len());
        for _ in 0..added_copies {
            inserted.extend_from_slice(canonical_unit);
        }
        let flank_left = index_to_hgvs_pos(ref_end - 1);
        let flank_right = flank_left + 1;
        RepeatNormResult::Insertion {
            start: flank_left,
            end: flank_right,
            sequence: inserted,
        }
    } else {
        // Default expansion (>=2 copies, not gated): repeat notation.
        // The repeat region describes the REFERENCE tract (per HGVS spec)
        // ref_end is exclusive, so last position is ref_end - 1
        RepeatNormResult::Repeat {
            start: index_to_hgvs_pos(ref_start),
            end: index_to_hgvs_pos(ref_end - 1),
            sequence: canonical_unit.to_vec(),
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
        // A4: a substitution where ref == alt (e.g. `c.100A>A`) is degenerate;
        // the HGVS v21 spec calls the form "not allowed" (recommendations/DNA/
        // other.md) and gives `c.100=` as the canonical alternative. The rule
        // is purely syntactic on the edit's stated bases, so it lives here in
        // the no-reference canonicalization path alongside the other minimal-
        // notation rewrites — it must fire regardless of provider availability.
        NaEdit::Substitution {
            reference,
            alternative,
        } if reference == alternative => NaEdit::position_identity(),
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
        // Companion to the A4 arm in `canonicalize_edit`: route degenerate
        // substitutions through the no-reference canonicalize path so the
        // rewrite fires even when the provider is empty.
        NaEdit::Substitution {
            reference,
            alternative,
        } => reference == alternative,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_needs_normalization() {
        use crate::hgvs::edit::Base;

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
        // Real substitutions stay on the no-op fast path.
        assert!(!needs_normalization(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        }));
        // Degenerate `A>A` substitutions must enter normalization so they
        // can be rewritten to identity (`=`).
        assert!(needs_normalization(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A,
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
    fn test_insertion_is_duplication_pos_beyond_ref() {
        // Regression: inverted-range insertions (start > end) like
        // NC_000011.10:g.5238138_5153222insTATTT produce a pos far
        // beyond the reference sequence length.  Must return false,
        // not panic.
        let ref_seq = b"ATGATGATG";
        assert!(!insertion_is_duplication(ref_seq, 95, b"TATTT"));
        assert!(!insertion_is_duplication(b"", 95, b"TATTT"));
        assert!(!insertion_is_duplication(b"", 0, b"A"));
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
    fn test_delins_is_identity() {
        let ref_seq = b"ACGTACGT";

        // Single-base same → identity
        assert!(delins_is_identity(ref_seq, 3, 4, b"T"));

        // Multi-base same → identity
        assert!(delins_is_identity(ref_seq, 1, 4, b"CGT"));

        // Whole region, full match → identity
        assert!(delins_is_identity(ref_seq, 0, 8, b"ACGTACGT"));

        // Single-base differ → not identity
        assert!(!delins_is_identity(ref_seq, 3, 4, b"A"));

        // Multi-base differ at one position → not identity
        assert!(!delins_is_identity(ref_seq, 1, 4, b"CGA"));

        // Length mismatch (1→2) → not identity
        assert!(!delins_is_identity(ref_seq, 3, 4, b"TT"));

        // Length mismatch (2→1) → not identity
        assert!(!delins_is_identity(ref_seq, 1, 3, b"C"));

        // Empty insert → not identity
        assert!(!delins_is_identity(ref_seq, 3, 4, b""));

        // OOB end → false (no panic)
        assert!(!delins_is_identity(ref_seq, 7, 9, b"TT"));

        // Inverted range → false
        assert!(!delins_is_identity(ref_seq, 5, 3, b"CG"));
    }

    #[test]
    fn test_delins_is_substitution() {
        use crate::hgvs::edit::Base;

        let ref_seq = b"ACGTACGT";
        // ref[3] = T, insert A → T>A
        assert_eq!(
            delins_is_substitution(ref_seq, 3, 4, b"A"),
            Some((Base::T, Base::A))
        );

        // Single-base delins where ref == alt → not a substitution
        // (would be identity, handled separately)
        assert_eq!(delins_is_substitution(ref_seq, 3, 4, b"T"), None);

        // Multi-base delete → not a single-base substitution
        assert_eq!(delins_is_substitution(ref_seq, 3, 5, b"A"), None);

        // Single-base delete with multi-base insert → not a substitution
        assert_eq!(delins_is_substitution(ref_seq, 3, 4, b"AT"), None);

        // Bounds: end past reference → None (no panic)
        assert_eq!(delins_is_substitution(ref_seq, 8, 9, b"A"), None);

        // Empty insert → None
        assert_eq!(delins_is_substitution(ref_seq, 3, 4, b""), None);
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
        // Sequence: GGGAAAAAGGG (5 A's at positions 3-7, 0-indexed inclusive)
        let ref_seq = b"GGGAAAAAGGG";

        // `pos` is the 0-based index of the base immediately 5' of the
        // insertion point. pos=7 means insert between index 7 and 8 — i.e.,
        // immediately after the last A in the homopolymer. Inserting AA here
        // should become A[7] (5 ref + 2 inserted).
        let result = insertion_to_repeat(ref_seq, 7, b"AA", false);
        assert!(result.is_some());
        let (base, count, start, end, _unit) = result.unwrap();
        assert_eq!(base, b'A');
        assert_eq!(count, 7); // 5 original + 2 inserted
        assert_eq!(start, 4); // 1-indexed start of A region in reference
        assert_eq!(end, 8); // 1-indexed end of A region in reference (per HGVS, positions refer to reference tract)

        // Single-copy inserts (added_copies < 2) are duplications, not repeats.
        let result = insertion_to_repeat(ref_seq, 7, b"A", false);
        assert!(result.is_none());

        // Inserting T (non-matching) should return None
        let result = insertion_to_repeat(ref_seq, 7, b"T", false);
        assert!(result.is_none());

        // Inserting mixed bases should return None
        let result = insertion_to_repeat(ref_seq, 7, b"AT", false);
        assert!(result.is_none());

        // Multi-base tandem unit: insert ACAC into ACAC tract → AC[4].
        // ref_seq has ACAC at indices 0-3. Insert ACAC just before index 4
        // (pos=3, between A at index 3 and base at index 4). Expected:
        // 2 ref AC units + 2 inserted AC units = AC[4] at ref indices 0..3.
        let ref_ac = b"ACACGGG";
        let result = insertion_to_repeat(ref_ac, 3, b"ACAC", false);
        assert!(result.is_some());
        let (base, count, start, end, unit) = result.unwrap();
        assert_eq!(base, b'A');
        assert_eq!(count, 4);
        assert_eq!(start, 1); // 1-indexed
        assert_eq!(end, 4); // 1-indexed
        assert_eq!(unit, b"AC");
    }

    #[test]
    fn test_duplication_to_repeat() {
        // Sequence: GGGAAAAAGGG (5 A's at positions 3-7, 0-indexed)
        let ref_seq = b"GGGAAAAAGGG";

        // Duplicating 2 A's (positions 3-5, 0-indexed) should become A[7]
        let result = duplication_to_repeat(ref_seq, 3, 5, false);
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
        let result = duplication_to_repeat(non_repeat_seq, 0, 3, false);
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
        let result = duplication_to_repeat(ref_seq, 5, 8, false);
        assert!(result.is_none(), "Single-copy dup should not become repeat");

        // Duplicating GCAGCA (2 copies) SHOULD become repeat notation
        let result = duplication_to_repeat(ref_seq, 5, 11, false);
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
        // Sequence: GGGCATCATCATCATGGG (4 CAT repeats at positions 3-14, 0-indexed)
        // Specifying CAT[1]: ref_count=4, specified=1, k=3 >= 2, post=1 >= 1 → B2 → Repeat
        let ref_seq = b"GGGCATCATCATCATGGG";

        let result = normalize_repeat(ref_seq, 3, b"CAT", 1, false);
        match result {
            RepeatNormResult::Repeat {
                sequence, count, ..
            } => {
                assert_eq!(sequence, b"CAT");
                assert_eq!(count, 1, "Should reflect specified count of 1");
            }
            _ => panic!("Expected Repeat (B2), got {:?}", result),
        }
    }

    #[test]
    fn test_normalize_repeat_to_duplication() {
        // Sequence: GGGCATCATGGG (2 CAT repeats at positions 3-8)
        // Specifying CAT[3] (ref is 2, so 2+1=3) should become duplication
        let ref_seq = b"GGGCATCATGGG";

        let result = normalize_repeat(ref_seq, 3, b"CAT", 3, false);
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

        let result = normalize_repeat(ref_seq, 3, b"CAT", 5, false);
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

        let result = normalize_repeat(ref_seq, 3, b"CAT", 2, false);
        assert!(matches!(result, RepeatNormResult::Unchanged));
    }

    #[test]
    fn test_normalize_repeat_empty_unit_is_unchanged() {
        // Pre-refactor, `normalize_repeat` delegated emptiness handling to
        // `count_tandem_repeats`, which returns `None` on an empty unit, yielding
        // `Unchanged`. The B2 canonicalization step now calls
        // `smallest_repeat_unit` first, which returns the empty slice unchanged
        // for an empty input — leaving a divide-by-zero on `repeat_unit.len() /
        // canonical_unit.len()` unless we guard up front. This test pins the
        // pre-refactor contract.
        let ref_seq = b"GGGCATCATGGG";
        let result = normalize_repeat(ref_seq, 3, b"", 1, false);
        assert!(matches!(result, RepeatNormResult::Unchanged));
    }

    #[test]
    fn test_normalize_repeat_canonicalizes_non_minimal_unit() {
        // Caller-spelled `ATAT[1]` over an `AT[4]` tract is a contraction:
        // canonical unit AT, ref_count=4, specified_canonical=2, k=2 → B2
        // emits `AT[2]` at the canonical tract span. Without canonicalization
        // this would fall to a 1-unit (ATAT) reduction → deletion (k<2).
        let ref_seq = b"GGGATATATATGGG"; // AT-tract at indices 3..11 (4 AT)

        let result = normalize_repeat(ref_seq, 3, b"ATAT", 1, false);
        match result {
            RepeatNormResult::Repeat {
                start,
                end,
                sequence,
                count,
            } => {
                assert_eq!(sequence, b"AT", "Should emit canonical (smallest) unit");
                assert_eq!(count, 2, "Specified ATAT[1] = 2 canonical AT copies");
                assert_eq!(start, 4, "Canonical tract starts at HGVS pos 4");
                assert_eq!(end, 11, "Canonical tract ends at HGVS pos 11");
            }
            _ => panic!("Expected canonical AT[2] Repeat, got {:?}", result),
        }
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
    // EXTEND_TANDEM_TRACT TESTS
    // =========================================================================

    #[test]
    fn test_extend_tandem_tract_homopolymer() {
        // ref: "TTAAAATT", anchor [2,4) (the first AA), unit "A"
        let ref_seq = b"TTAAAATT";
        let tract = extend_tandem_tract(ref_seq, 2..4, b"A").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 6);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_multi_base_unit() {
        // ref: "TTGCAGCAGCATT", anchor [5,8) (middle GCA), unit "GCA"
        let ref_seq = b"TTGCAGCAGCATT";
        let tract = extend_tandem_tract(ref_seq, 5..8, b"GCA").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 11);
        assert_eq!(tract.ref_count, 3);
    }

    #[test]
    fn test_extend_tandem_tract_anchor_at_5prime_edge() {
        // ref: "AAAATT", anchor [0,2), unit "A"
        let ref_seq = b"AAAATT";
        let tract = extend_tandem_tract(ref_seq, 0..2, b"A").unwrap();
        assert_eq!(tract.start, 0);
        assert_eq!(tract.end, 4);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_anchor_at_3prime_edge() {
        // ref: "TTAAAA", anchor [4,6), unit "A"
        let ref_seq = b"TTAAAA";
        let tract = extend_tandem_tract(ref_seq, 4..6, b"A").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 6);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_zero_width_anchor() {
        // Zero-width anchor (insertion-point semantics): the anchor itself is
        // trivially unit-periodic; the helper just walks left/right.
        // ref: "TTAAAATT", anchor [4,4), unit "A"
        let ref_seq = b"TTAAAATT";
        let tract = extend_tandem_tract(ref_seq, 4..4, b"A").unwrap();
        assert_eq!(tract.start, 2);
        assert_eq!(tract.end, 6);
        assert_eq!(tract.ref_count, 4);
    }

    #[test]
    fn test_extend_tandem_tract_anchor_not_unit_periodic() {
        // ref: "TTAGAATT", anchor [2,4) is "AG" — not "A"-periodic.
        let ref_seq = b"TTAGAATT";
        assert!(extend_tandem_tract(ref_seq, 2..4, b"A").is_none());
    }

    #[test]
    fn test_extend_tandem_tract_anchor_length_not_multiple_of_unit() {
        // ref: "AAA", anchor [0,3), unit "AA" — len 3 % 2 == 1.
        let ref_seq = b"AAA";
        assert!(extend_tandem_tract(ref_seq, 0..3, b"AA").is_none());
    }

    // =========================================================================
    // DELETION_TO_REPEAT TESTS
    // =========================================================================

    #[test]
    fn test_deletion_to_repeat_homopolymer_two_removed() {
        // ref "TTAAAAATT" (5 A's at indices 2..7). Delete 2 A's at [4..6).
        // After 3' shift, del lands at [5..7). post-shift slice "AA", unit "A", k=2.
        // Tract [2..7), ref_count=5, post_count=3 → A[3] at HGVS [3..7] (1-based).
        let ref_seq = b"TTAAAAATT";
        let r = deletion_to_repeat(ref_seq, 5, 7, false).expect("should fire");
        assert_eq!(r.unit, b"A");
        assert_eq!(r.count, 3);
        assert_eq!(r.start, 3); // 1-based HGVS
        assert_eq!(r.end, 7);
    }

    #[test]
    fn test_deletion_to_repeat_multi_base_tandem_two_removed() {
        // ref "TTGCAGCAGCATT" (3 GCAs at [2..11)). Delete 2 GCAs at [5..11).
        // post-shift slice "GCAGCA", unit "GCA", k=2. Tract [2..11), ref_count=3,
        // post_count=1 → GCA[1] at HGVS [3..11] (1-based).
        let ref_seq = b"TTGCAGCAGCATT";
        let r = deletion_to_repeat(ref_seq, 5, 11, false).expect("should fire");
        assert_eq!(r.unit, b"GCA");
        assert_eq!(r.count, 1);
        assert_eq!(r.start, 3);
        assert_eq!(r.end, 11);
    }

    #[test]
    fn test_deletion_to_repeat_one_unit_removed_returns_none() {
        // k=1 → stays as del.
        // ref "TTAAAAATT", delete 1 A at [6..7).
        let ref_seq = b"TTAAAAATT";
        assert!(deletion_to_repeat(ref_seq, 6, 7, false).is_none());
    }

    #[test]
    fn test_deletion_to_repeat_full_tract_removal_returns_none() {
        // post_count == 0 → stays as del.
        // ref "TTAATT", delete both A's at [2..4). ref_count=2, k=2, post_count=0.
        let ref_seq = b"TTAATT";
        assert!(deletion_to_repeat(ref_seq, 2, 4, false).is_none());
    }

    #[test]
    fn test_deletion_to_repeat_non_tandem_returns_none() {
        // Single isolated unit, ref_count < 2.
        // ref "TTGCATT", delete "GCA" at [2..5). smallest_repeat_unit("GCA")="GCA",
        // ref_count=1, returns None.
        let ref_seq = b"TTGCATT";
        assert!(deletion_to_repeat(ref_seq, 2, 5, false).is_none());
    }

    #[test]
    fn test_deletion_to_repeat_finer_periodicity() {
        // ref "TTATATATATATT" (5 ATs at [2..12)). Delete "ATAT" at [8..12).
        // smallest_repeat_unit("ATAT")="AT", k=2, ref_count=5, post_count=3 → AT[3].
        let ref_seq = b"TTATATATATATT";
        let r = deletion_to_repeat(ref_seq, 8, 12, false).expect("should fire");
        assert_eq!(r.unit, b"AT");
        assert_eq!(r.count, 3);
        assert_eq!(r.start, 3); // 1-based HGVS [3..12]
        assert_eq!(r.end, 12);
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

        // Real substitutions stay on the no-op fast path
        use crate::hgvs::edit::Base;
        assert!(!should_canonicalize(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G
        }));

        // Degenerate `A>A` substitutions enter the no-reference canonicalize
        // path so they can be rewritten to identity (`=`) even when the
        // provider has no transcript / genomic sequence loaded.
        assert!(should_canonicalize(&NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A
        }));
    }

    #[test]
    fn test_canonicalize_edit_degenerate_substitution_to_identity() {
        // A4: c.100A>A is "not allowed" per HGVS v21 spec; canonicalize to `=`.
        use crate::hgvs::edit::Base;

        let degenerate = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::A,
        };
        assert_eq!(canonicalize_edit(&degenerate), NaEdit::position_identity());

        // Real SNVs pass through unchanged.
        let real_sub = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        assert_eq!(canonicalize_edit(&real_sub), real_sub);
    }

    // =========================================================================
    // Codon-frame gate tests (#81 B1)
    //
    // HGVS spec (docs/recommendations/DNA/repeated.md, §Notes):
    // > using a coding DNA reference sequence ("c." description) a Repeated
    // > sequence variant description can be used only for repeat units with
    // > a length which is a multiple of 3, i.e. which can not affect the
    // > reading frame.
    // =========================================================================

    #[test]
    fn test_insertion_to_repeat_codon_frame_gate_blocks_a_in_coding() {
        // Reference: 5-A homopolymer flanked by Cs. Insert AA → would normally
        // emit A[7], but unit_len=1 is not a multiple of 3 in c. context,
        // so the gate returns None.
        let ref_seq = b"CAAAAAC";
        let result = insertion_to_repeat(ref_seq, 5, b"AA", true);
        assert!(
            result.is_none(),
            "is_coding=true + unit_len=1 must return None"
        );
    }

    #[test]
    fn test_insertion_to_repeat_codon_frame_gate_blocks_at_in_coding() {
        // Reference: AT[3] tandem flanked by Cs. Insert ATAT → unit_len=2,
        // gate blocks in coding context.
        let ref_seq = b"CATATATC";
        let result = insertion_to_repeat(ref_seq, 6, b"ATAT", true);
        assert!(
            result.is_none(),
            "is_coding=true + unit_len=2 must return None"
        );
    }

    #[test]
    fn test_insertion_to_repeat_codon_frame_gate_passes_cag_in_coding() {
        // Reference: CAG[3] tandem. Insert CAGCAG → unit_len=3, codon-aligned,
        // gate passes; result is Some(...) carrying CAG[5].
        let ref_seq = b"CCAGCAGCAGT";
        let result = insertion_to_repeat(ref_seq, 9, b"CAGCAG", true);
        assert!(
            result.is_some(),
            "is_coding=true + unit_len=3 must allow rewrite"
        );
        let (_first, count, _start, _end, unit) = result.unwrap();
        assert_eq!(count, 5, "expected CAG[5]");
        assert_eq!(unit, b"CAG");
    }

    #[test]
    fn test_insertion_to_repeat_gate_no_op_in_genomic() {
        // Same A-homopolymer case as the blocking test, but is_coding=false
        // → gate is a no-op, repeat rewrite proceeds.
        let ref_seq = b"CAAAAAC";
        let result = insertion_to_repeat(ref_seq, 5, b"AA", false);
        assert!(result.is_some(), "is_coding=false must not gate");
    }

    #[test]
    fn test_deletion_to_repeat_codon_frame_gate_blocks_a_in_coding() {
        // 5-A tract, delete 2 A's. Span 2..4 covers two of the As so that
        // without the codon-frame gate the function would return Some(A[3]);
        // with the gate (coding) it must return None. Span 2..3 would also
        // return None via the unrelated `k < 2` early exit, so it would not
        // discriminate the gate.
        let ref_seq = b"CAAAAAC";
        let result = deletion_to_repeat(ref_seq, 2, 4, true);
        assert!(
            result.is_none(),
            "is_coding=true + unit_len=1 must return None"
        );
    }

    #[test]
    fn test_deletion_to_repeat_codon_frame_gate_passes_cag_in_coding() {
        // ref "CCAGCAGCAGT": 3-CAG tract at indices 1..10. Delete 2 CAGs at
        // [1..7) (6 bases CAGCAG). With codon-aligned unit, gate passes.
        let ref_seq = b"CCAGCAGCAGT";
        let result = deletion_to_repeat(ref_seq, 1, 7, true);
        assert!(
            result.is_some(),
            "is_coding=true + unit_len=3 must allow rewrite"
        );
    }

    #[test]
    fn test_duplication_to_repeat_codon_frame_gate_routes_a_to_gated_insertion() {
        // 4-A tract at indices 1..5 (0-based). Duplicate 2 A's at positions 1..3.
        // Under the gate, structural conditions for repeat notation are met
        // but unit_len=1 in c. is forbidden, so the result routes to a
        // GatedInsertion that the caller renders as `ins<dup_seq>`.
        let ref_seq = b"CAAAAC";
        let result = duplication_to_repeat(ref_seq, 1, 3, true);
        match result {
            Some(DupToRepeatResult::GatedInsertion { sequence, .. }) => {
                assert_eq!(sequence, b"AA", "sequence is the duplicated literal");
            }
            other => panic!("expected GatedInsertion, got {:?}", other),
        }
    }

    #[test]
    fn test_duplication_to_repeat_codon_frame_gate_passes_cag_in_coding() {
        // 3-CAG tract. Duplicate 2 CAGs.
        let ref_seq = b"CCAGCAGCAGT";
        let result = duplication_to_repeat(ref_seq, 1, 7, true);
        assert!(
            result.is_some(),
            "is_coding=true + unit_len=3 must allow rewrite"
        );
    }

    #[test]
    fn test_normalize_repeat_codon_frame_gate_routes_contraction_to_deletion() {
        // 5-A tract, specified A[3] in coding → must NOT emit Repeat; emits Deletion.
        let ref_seq = b"CAAAAAC";
        let result = normalize_repeat(ref_seq, 1, b"A", 3, true);
        match result {
            RepeatNormResult::Deletion { .. } => {}
            other => panic!("expected Deletion under gate, got {:?}", other),
        }
    }

    #[test]
    fn test_normalize_repeat_codon_frame_gate_routes_expansion_to_insertion() {
        // 5-A tract, specified A[8] in coding → must NOT emit Repeat; emits Insertion.
        let ref_seq = b"CAAAAAC";
        let result = normalize_repeat(ref_seq, 1, b"A", 8, true);
        match result {
            RepeatNormResult::Insertion { sequence, .. } => {
                assert_eq!(sequence, b"AAA", "3 extra A's");
            }
            other => panic!("expected Insertion under gate, got {:?}", other),
        }
    }

    #[test]
    fn test_normalize_repeat_codon_frame_gate_passes_through_dup_branch() {
        // 5-A tract, specified A[6] in coding → +1 copy = dup, gate doesn't change this.
        let ref_seq = b"CAAAAAC";
        let result = normalize_repeat(ref_seq, 1, b"A", 6, true);
        match result {
            RepeatNormResult::Duplication { .. } => {}
            other => panic!("expected Duplication, got {:?}", other),
        }
    }
}
