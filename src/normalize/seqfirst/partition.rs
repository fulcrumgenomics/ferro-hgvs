//! Unit 3: cut the reference block into disjoint member blocks.
//!
//! Members are runs of *change*, where a reference position counts as unchanged
//! only when the same match **edge** is present in every minimal alignment
//! (node agreement is not enough — see `align.rs`'s `Dominators` doc), and a
//! junction counts as changed when every minimal alignment inserts there.
//! Consecutive runs merge when separated by fewer than [`MIN_SEPARATION`]
//! unchanged bases.
//!
//! Members come out disjoint **by construction** — they are a partition — so
//! the overlap and ordering defects the repair passes exist to fix become
//! unrepresentable rather than repaired.

use super::align::AlignmentDag;
use super::MIN_SEPARATION;

/// One member of the canonicalized allele, as a half-open reference span
/// relative to the start of the reference block.
///
/// A pure insertion has `ref_start == ref_end`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MemberBlock {
    pub(crate) ref_start: u32,
    pub(crate) ref_end: u32,
}

/// Partition the reference block at [`MIN_SEPARATION`].
pub(crate) fn partition_members(dag: &AlignmentDag) -> Vec<MemberBlock> {
    partition_members_with(dag, MIN_SEPARATION)
}

/// Partition the reference block, merging runs separated by fewer than
/// `min_separation` **unchanged reference bases**.
///
/// The unit of separation is unchanged bases, not event indices. Two events at
/// reference offsets `p` and `q` are separated by `q - p - 1` bases. Comparing
/// `q - p` instead splits the spec's own worked example (`LRG_199t1:c.850_901`)
/// into two members where the spec describes one; the other calibration cases
/// score identically under both readings, so that error is invisible without
/// that case.
pub(crate) fn partition_members_with(dag: &AlignmentDag, min_separation: u32) -> Vec<MemberBlock> {
    let dominators = dag.dominators();
    let ref_len = dag.ref_len();

    // An event is a reference position no minimal alignment agrees is matched,
    // or a junction every minimal alignment inserts at. Both live in reference
    // offset space: junction `k` sits between reference bases `k - 1` and `k`.
    let matched_ref = dominators.matched_ref();
    let changed: Vec<u32> = (0..ref_len)
        .filter(|p| matched_ref.binary_search(p).is_err())
        .collect();

    let mut events: Vec<u32> = changed.clone();
    events.extend(dominators.forced_ins_junctions.iter().copied());
    events.sort_unstable();
    events.dedup();

    // Group consecutive events into runs, merging across short gaps.
    let mut runs: Vec<(u32, u32)> = Vec::new();
    for event in events {
        match runs.last_mut() {
            // `event - last_end - 1` is the count of unchanged bases between
            // the previous event and this one.
            Some((_, last_end))
                if event.saturating_sub(*last_end).saturating_sub(1) < min_separation =>
            {
                *last_end = event;
            }
            _ => runs.push((event, event)),
        }
    }

    runs.into_iter()
        .map(|(start, end)| MemberBlock {
            ref_start: start,
            // A run's last event is either a changed reference position, which
            // occupies `end..end + 1`, or a forced-insertion junction, which
            // occupies no reference base at all and so ends at `end`.
            //
            // Widening an insertion-only end to `end + 1` would claim a
            // reference base that every minimal alignment matches, giving the
            // member a non-minimal extent. Both readings denote the same
            // sequence and round-trip, so the round-trip invariant does not
            // catch this. The exact-span assertions in
            // `insertion_junctions_do_not_claim_a_matched_base` and
            // `an_insertion_only_member_has_an_empty_reference_span` do.
            //
            // When a position is both changed and a forced-insertion junction,
            // the changed reading wins — it is the wider of the two and the
            // base really does change.
            ref_end: if changed.binary_search(&end).is_ok() {
                (end + 1).min(ref_len)
            } else {
                end
            },
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The calibration corpus. Seven cases; 7/7 is the bar this design cleared
    /// and no other candidate did (score-margin and match-density both failed
    /// here — see `separations_are_meaningful` in `merge.rs` for their
    /// remains).
    const LRG199_REF: &[u8] = b"CAGGGATATGAGAGAACTTCTTCCCCTAAGCCTCGATTCAAGAGCTATGCCT";

    #[test]
    fn calibration_corpus_partitions_as_measured() {
        for (name, r, a, expected) in [
            // Dominance alone is not enough here: every reference position is
            // matched, and it is the separation rule that merges the two
            // forced insertion junctions.
            (
                "#1260 two insC at adjacent gaps",
                &b"AAAAAAA"[..],
                &b"AACACAAAA"[..],
                1,
            ),
            // Separation is not enough here: there are no dominators at all,
            // so nothing can legitimately be split.
            (
                "#1262 sub + del vs spanning delins",
                &b"AAAAAAA"[..],
                &b"ACAAAA"[..],
                1,
            ),
            // The spec's own worked example. One dominator, which would split
            // the block if the separation rule did not merge across it.
            (
                "delins.md:44 LRG_199t1 c.850_901",
                LRG199_REF,
                &b"TTCCTCGATGCCTG"[..],
                1,
            ),
            // Genuine splits must survive.
            (
                "control two subs 12 apart",
                &b"ACGTACGTACGTACGT"[..],
                &b"ATGTACGTACGTACTT"[..],
                2,
            ),
            (
                "control two subs exactly 2 apart",
                &b"ACGTACGT"[..],
                &b"ATGTTCGT"[..],
                2,
            ),
            (
                "control two insertions 5 apart",
                &b"ACGTACGTAC"[..],
                &b"ACGTTACGTGAC"[..],
                2,
            ),
        ] {
            let dag = AlignmentDag::build(r, a);
            let members = partition_members(&dag);
            assert_eq!(members.len(), expected, "{name}: got {members:?}");
        }
    }

    /// The one cheap case that distinguishes counting unchanged bases from
    /// comparing event indices.
    ///
    /// Two substitutions with exactly **one** unchanged base between them must
    /// merge, because one is fewer than `MIN_SEPARATION`. Comparing event
    /// indices instead gives `4 - 2 = 2`, which is not less than 2, so the
    /// wrong rule splits this into two members.
    ///
    /// Without this test and the `LRG_199` row above, that off-by-one is
    /// invisible — the other five calibration cases score the same either way.
    #[test]
    fn one_unchanged_base_between_runs_merges() {
        let dag = AlignmentDag::build(b"GACTGACTGA", b"GAATAACTGA");
        assert_eq!(
            partition_members(&dag).len(),
            1,
            "one unchanged base is fewer than MIN_SEPARATION, so these merge"
        );
    }

    #[test]
    fn two_unchanged_bases_between_runs_split() {
        // The boundary from the other side: exactly MIN_SEPARATION unchanged
        // bases must NOT merge.
        let dag = AlignmentDag::build(b"GACTGACTGA", b"GAATGTCTGA");
        assert_eq!(partition_members(&dag).len(), 2);
    }

    #[test]
    fn an_unchanged_block_has_no_members() {
        let dag = AlignmentDag::build(b"ACGTACGT", b"ACGTACGT");
        assert!(partition_members(&dag).is_empty());
    }

    #[test]
    fn members_are_disjoint_and_ascending() {
        let dag = AlignmentDag::build(b"ACGTACGTACGTACGT", b"ATGTACGTACGTACTT");
        let members = partition_members(&dag);
        for pair in members.windows(2) {
            assert!(
                pair[0].ref_end <= pair[1].ref_start,
                "members overlap or are out of order: {members:?}"
            );
        }
    }

    #[test]
    fn insertion_junctions_do_not_claim_a_matched_base() {
        // #1260 collapses to one member spanning the two forced junctions and
        // the single matched base between them: reference 2..3, exactly.
        //
        // The failure this pins is widening the end to 4 because the run's last
        // event is at offset 3. Offset 3 is an insertion *junction*, not a
        // changed base, so it occupies no reference base. Getting this wrong
        // still yields one member — the count looks right, and both the wider
        // and narrower spans denote the same sequence, so the round-trip
        // invariant does not catch it either. This exact-span assertion does.
        let dag = AlignmentDag::build(b"AAAAAAA", b"AACACAAAA");
        let members = partition_members(&dag);
        assert_eq!(
            members,
            vec![MemberBlock {
                ref_start: 2,
                ref_end: 3
            }]
        );
    }

    #[test]
    fn an_insertion_only_member_has_an_empty_reference_span() {
        // Two insertions 5 apart: the second is a pure insertion at junction 8
        // with an empty reference span, which must be reported and not dropped.
        let dag = AlignmentDag::build(b"ACGTACGTAC", b"ACGTTACGTGAC");
        let members = partition_members(&dag);
        assert_eq!(
            members,
            vec![
                MemberBlock {
                    ref_start: 3,
                    ref_end: 4
                },
                MemberBlock {
                    ref_start: 8,
                    ref_end: 8
                },
            ]
        );
    }
}
