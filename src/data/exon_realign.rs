//! Sequence-aware correction of transcript↔genome exon coordinate mapping.
//!
//! cdot's per-exon alignment is normally enough to map a genomic position to a
//! transcript position by plain offset arithmetic (`tx_off = genome_off`).
//! That arithmetic is correct only when the exon's transcript and genome
//! sequences are the same length and identical base-for-base. Some cdot
//! snapshots (e.g. cdot 0.2.32 for `NM_000532.5` / PCCB) record an exon as
//! *ungapped* (`exon_cigars[i] == None`) even though the transcript FASTA and
//! the genome FASTA genuinely differ by an indel. Without the gap in the
//! alignment, the plain arithmetic places every position after the indel off by
//! the indel's size, and the projected `c.` coordinate is wrong (issue #644).
//!
//! This module re-derives the transcript↔genome map for one exon directly from
//! the two sequences via a small global alignment, and exposes the corrective
//! delta to apply to the naive transcript offset. It is used only as a fallback
//! when cdot reports the exon ungapped but the sequences disagree, so the common
//! (identical-sequence) case is untouched.
//!
//! # Orientation
//!
//! Both slices are expected in **transcript 5'→3' orientation**: the caller
//! reverse-complements the genome slice on a minus-strand transcript before
//! calling in, mirroring the offset convention used throughout
//! [`crate::data::cdot`]. Offsets are then measured from the exon's
//! transcript-5' end on both axes.

/// Result of correcting a single genome offset against the realigned exon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExonOffsetCorrection {
    /// The corrected transcript offset for the queried genome offset, expressed
    /// as a signed delta to add to the *naive* transcript offset (which equals
    /// the genome offset). `0` means no correction is needed.
    Delta(i64),
    /// The genome offset falls strictly inside a region of genome bases that
    /// have no transcript counterpart (a genome insertion relative to the
    /// transcript). Such a position has no well-defined transcript coordinate,
    /// so the caller should decline rather than emit a wrong one.
    InsideGenomeOnlyGap,
}

/// Maximum exon slice length (in bases) for which we attempt realignment.
///
/// The fallback alignment is `O(G * T)` time and memory. Real exons that
/// trigger this path are short (the PCCB exon is 219 bp), and an exon whose
/// two slices differ by a large amount is far more likely to be a data error
/// (wrong sequence pairing) than a genuine small indel — declining is safer
/// than spending quadratic effort to "correct" garbage. Beyond this bound we
/// return `None` (decline).
const MAX_REALIGN_LEN: usize = 4096;

/// Re-derive the transcript offset for a genome offset within one ungapped-but-
/// divergent exon, returning the correction to apply to the naive offset.
///
/// `genome` and `tx` are the exon's genome and transcript slices, both in
/// transcript 5'→3' orientation (see module docs). `genome_offset` is the
/// 0-based offset of the queried position from the exon's transcript-5' end on
/// the genome axis (i.e. the value the naive arithmetic would also use as the
/// transcript offset).
///
/// Returns:
/// - `Some(Delta(d))` when the position maps to a well-defined transcript
///   offset `genome_offset + d`.
/// - `Some(InsideGenomeOnlyGap)` when the position sits strictly inside a
///   genome-only gap and has no transcript counterpart.
/// - `None` when realignment is not attempted or is not confident enough to
///   correct (sequences identical/equal-length with no indel, a slice empty,
///   slices too long, or the queried offset is out of range) — the caller
///   should fall back to the naive arithmetic unchanged.
///
/// When the two slices are byte-identical this returns `None` (nothing to
/// correct), so the hot path for normal exons is a single length check plus an
/// equality scan.
pub fn correct_genome_offset(
    genome: &[u8],
    tx: &[u8],
    genome_offset: usize,
) -> Option<ExonOffsetCorrection> {
    // Nothing to correct when the sequences already agree exactly. This is the
    // overwhelmingly common case (cdot's alignment is right); short-circuit it.
    if genome == tx {
        return None;
    }
    // Need both sequences to align against; an empty slice is a data problem we
    // decline rather than guess at.
    if genome.is_empty() || tx.is_empty() {
        return None;
    }
    // Guard the quadratic cost / reject implausible pairings.
    if genome.len() > MAX_REALIGN_LEN || tx.len() > MAX_REALIGN_LEN {
        return None;
    }
    if genome_offset >= genome.len() {
        return None;
    }

    let alignment = align(genome, tx);
    correction_at(&alignment, genome_offset)
}

/// Result of correcting a single transcript offset against the realigned exon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TxOffsetCorrection {
    /// The corrected genome offset for the queried transcript offset, expressed
    /// as a signed delta to add to the *naive* genome offset (which equals the
    /// transcript offset). `0` means no correction is needed.
    Delta(i64),
    /// The transcript offset falls strictly inside a region of transcript bases
    /// that have no genome counterpart (a transcript insertion relative to the
    /// genome). Such a position has no well-defined genome coordinate, so the
    /// caller should decline rather than emit a wrong one.
    InsideTxOnlyGap,
}

/// The transcript→genome mirror of [`correct_genome_offset`], used by the
/// c.→g. projection path so a coordinate corrected on the way in round-trips on
/// the way out (issue #644).
///
/// `genome` and `tx` are the exon's slices in transcript 5'→3' orientation (see
/// module docs). `tx_offset` is the 0-based offset of the queried position from
/// the exon's transcript-5' end on the transcript axis. Returns the correction
/// to apply to the naive genome offset (which equals `tx_offset`), or `None`
/// when realignment is not attempted or not confident enough to correct.
pub fn correct_tx_offset(genome: &[u8], tx: &[u8], tx_offset: usize) -> Option<TxOffsetCorrection> {
    if genome == tx {
        return None;
    }
    if genome.is_empty() || tx.is_empty() {
        return None;
    }
    if genome.len() > MAX_REALIGN_LEN || tx.len() > MAX_REALIGN_LEN {
        return None;
    }
    if tx_offset >= tx.len() {
        return None;
    }
    let alignment = align(genome, tx);
    tx_correction_at(&alignment, tx_offset)
}

/// A single column of a genome↔transcript alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Col {
    /// Both axes consume a base (match or mismatch).
    Match,
    /// Genome consumes a base, transcript does not (genome-only base).
    GenomeOnly,
    /// Transcript consumes a base, genome does not (transcript-only base).
    TxOnly,
}

/// Global (Needleman–Wunsch) alignment of `genome` against `tx` with unit
/// mismatch and gap costs. Returns the column path in genome/tx order.
///
/// Linear-gap unit-cost NW is sufficient here: the divergences we correct are
/// small indels in otherwise-identical exon sequence, where the minimum-edit
/// alignment recovers the true indel placement. Traceback prefers diagonal
/// (match) moves on ties so a run of identical bases is consumed on the
/// diagonal rather than as paired gaps.
fn align(genome: &[u8], tx: &[u8]) -> Vec<Col> {
    let n = genome.len();
    let m = tx.len();
    let width = m + 1;
    // cost[i*width + j] = edit distance between genome[..i] and tx[..j].
    let mut cost = vec![0u32; (n + 1) * width];
    for i in 0..=n {
        cost[i * width] = i as u32;
    }
    // First row: edit distance of empty genome vs tx[..j] is j. Indexed
    // assignment mirrors the `cost[i * width]` column init above.
    #[allow(clippy::needless_range_loop)]
    for j in 0..=m {
        cost[j] = j as u32;
    }
    for i in 1..=n {
        for j in 1..=m {
            let sub =
                cost[(i - 1) * width + (j - 1)] + if genome[i - 1] == tx[j - 1] { 0 } else { 1 };
            let del = cost[(i - 1) * width + j] + 1; // genome-only
            let ins = cost[i * width + (j - 1)] + 1; // tx-only
            cost[i * width + j] = sub.min(del).min(ins);
        }
    }

    // Traceback from (n, m) to (0, 0).
    let mut path = Vec::with_capacity(n.max(m));
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let diag = cost[(i - 1) * width + (j - 1)];
            let here = cost[i * width + j];
            let matched = genome[i - 1] == tx[j - 1];
            // Prefer the diagonal move when it is consistent with `here`.
            if (matched && diag == here) || (!matched && diag + 1 == here) {
                path.push(Col::Match);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && cost[(i - 1) * width + j] + 1 == cost[i * width + j] {
            path.push(Col::GenomeOnly);
            i -= 1;
        } else {
            // j > 0 is guaranteed here because the loop condition holds and the
            // genome-only branch did not apply.
            path.push(Col::TxOnly);
            j -= 1;
        }
    }
    path.reverse();
    path
}

/// Walk the alignment columns, tracking genome and transcript offsets, and
/// report the correction for `genome_offset`.
fn correction_at(path: &[Col], genome_offset: usize) -> Option<ExonOffsetCorrection> {
    let mut g = 0usize; // genome bases consumed
    let mut t = 0usize; // transcript bases consumed
    for &col in path {
        match col {
            Col::Match => {
                if g == genome_offset {
                    return Some(ExonOffsetCorrection::Delta(t as i64 - g as i64));
                }
                g += 1;
                t += 1;
            }
            Col::GenomeOnly => {
                if g == genome_offset {
                    // The queried genome base aligns to no transcript base.
                    return Some(ExonOffsetCorrection::InsideGenomeOnlyGap);
                }
                g += 1;
            }
            Col::TxOnly => {
                t += 1;
            }
        }
    }
    None
}

/// Walk the alignment columns, tracking genome and transcript offsets, and
/// report the genome correction for `tx_offset`.
fn tx_correction_at(path: &[Col], tx_offset: usize) -> Option<TxOffsetCorrection> {
    let mut g = 0usize;
    let mut t = 0usize;
    for &col in path {
        match col {
            Col::Match => {
                if t == tx_offset {
                    return Some(TxOffsetCorrection::Delta(g as i64 - t as i64));
                }
                g += 1;
                t += 1;
            }
            Col::TxOnly => {
                if t == tx_offset {
                    // The queried transcript base aligns to no genome base.
                    return Some(TxOffsetCorrection::InsideTxOnlyGap);
                }
                t += 1;
            }
            Col::GenomeOnly => {
                g += 1;
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn delta(genome: &str, tx: &str, off: usize) -> Option<ExonOffsetCorrection> {
        correct_genome_offset(genome.as_bytes(), tx.as_bytes(), off)
    }

    #[test]
    fn identical_sequences_need_no_correction() {
        assert_eq!(delta("ACGTACGT", "ACGTACGT", 3), None);
    }

    #[test]
    fn genome_extra_prefix_shifts_following_positions() {
        // Genome has 2 extra bases ("CA") at the start; the rest matches the
        // transcript. This is the PCCB exon-start shape.
        let common = "GGTGCGCCGGTAGGGG";
        let genome = format!("CA{common}");
        let tx = common.to_string();
        // Positions inside the genome-only prefix have no transcript base.
        assert_eq!(
            delta(&genome, &tx, 0),
            Some(ExonOffsetCorrection::InsideGenomeOnlyGap)
        );
        assert_eq!(
            delta(&genome, &tx, 1),
            Some(ExonOffsetCorrection::InsideGenomeOnlyGap)
        );
        // The first matched genome base (offset 2) maps to tx offset 0 → -2.
        assert_eq!(
            delta(&genome, &tx, 2),
            Some(ExonOffsetCorrection::Delta(-2))
        );
        // A position deeper in the common region keeps the -2 correction.
        assert_eq!(
            delta(&genome, &tx, 9),
            Some(ExonOffsetCorrection::Delta(-2))
        );
    }

    #[test]
    fn tx_extra_suffix_shifts_positions_after_it() {
        // Transcript has an extra base in the middle (tx insertion). A genome
        // position after the insertion must shift UP by the insertion length.
        let genome = "ACGTACGT";
        let tx = "ACGTXACGT"; // one extra 'X' after offset 4
                              // Before the insertion: no correction.
        assert_eq!(delta(genome, tx, 0), Some(ExonOffsetCorrection::Delta(0)));
        assert_eq!(delta(genome, tx, 3), Some(ExonOffsetCorrection::Delta(0)));
        // At/after the insertion the genome offset maps one tx base higher.
        assert_eq!(delta(genome, tx, 4), Some(ExonOffsetCorrection::Delta(1)));
        assert_eq!(delta(genome, tx, 7), Some(ExonOffsetCorrection::Delta(1)));
    }

    #[test]
    fn two_indels_reproducer_shape() {
        // The PCCB exon shape: 2 extra genome bases at the start AND 1 extra
        // transcript base at the end. genome = "CA" + COMMON, tx = COMMON + "A".
        let common = "GGTGCGCCGGTAGGGGACGCGCCGGCACAGCAA";
        let genome = format!("CA{common}");
        let tx = format!("{common}A");
        // A position in the common middle: corrected by -2 (the genome-extra
        // prefix), unaffected by the trailing tx-extra base.
        let off = 2 + 10; // genome offset 12, inside common region
        assert_eq!(
            delta(&genome, &tx, off),
            Some(ExonOffsetCorrection::Delta(-2))
        );
        // The genome-only prefix bases decline.
        assert_eq!(
            delta(&genome, &tx, 0),
            Some(ExonOffsetCorrection::InsideGenomeOnlyGap)
        );
    }

    #[test]
    fn empty_slice_declines() {
        assert_eq!(delta("", "ACGT", 0), None);
        assert_eq!(delta("ACGT", "", 0), None);
    }

    #[test]
    fn out_of_range_offset_declines() {
        assert_eq!(delta("CAACGT", "ACGT", 99), None);
    }

    fn tx_delta(genome: &str, tx: &str, off: usize) -> Option<TxOffsetCorrection> {
        correct_tx_offset(genome.as_bytes(), tx.as_bytes(), off)
    }

    #[test]
    fn tx_correction_mirrors_genome_correction() {
        // Genome-extra prefix: genome = "CA" + COMMON, tx = COMMON.
        let common = "GGTGCGCCGGTAGGGG";
        let genome = format!("CA{common}");
        let tx = common.to_string();
        // tx offset 0 maps to genome offset 2 → +2 (mirror of the -2 forward).
        assert_eq!(
            tx_delta(&genome, &tx, 0),
            Some(TxOffsetCorrection::Delta(2))
        );
        assert_eq!(
            tx_delta(&genome, &tx, 7),
            Some(TxOffsetCorrection::Delta(2))
        );
    }

    #[test]
    fn tx_only_insertion_declines() {
        // Transcript has an extra base with no genome counterpart.
        let genome = "ACGTACGT";
        let tx = "ACGTXACGT"; // 'X' at tx offset 4 is tx-only
        assert_eq!(
            tx_delta(genome, tx, 4),
            Some(TxOffsetCorrection::InsideTxOnlyGap)
        );
        // Bases before/after the insertion map cleanly.
        assert_eq!(tx_delta(genome, tx, 3), Some(TxOffsetCorrection::Delta(0)));
        // tx offset 5 (the 'A' after 'X') maps to genome offset 4 → -1.
        assert_eq!(tx_delta(genome, tx, 5), Some(TxOffsetCorrection::Delta(-1)));
    }

    #[test]
    fn forward_and_inverse_round_trip() {
        // For every matched genome offset, the forward delta and inverse delta
        // must be exact negatives, so a position corrected in one direction
        // round-trips through the other.
        let genome = "CAGGTGCGCCGGTAGGGGA";
        let tx = "GGTGCGCCGGTAGGGGAT";
        for g_off in 0..genome.len() {
            if let Some(ExonOffsetCorrection::Delta(fwd)) =
                correct_genome_offset(genome.as_bytes(), tx.as_bytes(), g_off)
            {
                let t_off = (g_off as i64 + fwd) as usize;
                let inv = correct_tx_offset(genome.as_bytes(), tx.as_bytes(), t_off);
                assert_eq!(
                    inv,
                    Some(TxOffsetCorrection::Delta(-fwd)),
                    "genome offset {g_off}: forward {fwd}, inverse mismatch"
                );
            }
        }
    }

    #[test]
    fn equal_length_single_substitution_needs_no_shift() {
        // Same length, one mismatch: every position still maps 1:1 (delta 0),
        // because a substitution does not move coordinates.
        let genome = "ACGTACGT";
        let tx = "ACGAACGT"; // mismatch at offset 3
        assert_eq!(delta(genome, tx, 0), Some(ExonOffsetCorrection::Delta(0)));
        assert_eq!(delta(genome, tx, 3), Some(ExonOffsetCorrection::Delta(0)));
        assert_eq!(delta(genome, tx, 7), Some(ExonOffsetCorrection::Delta(0)));
    }
}
