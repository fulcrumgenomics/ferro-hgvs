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

use super::align::{AlignmentDag, Dominators, Step};
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
///
/// Computes [`AlignmentDag::dominators`] itself, so it is the convenient
/// entry point for a one-off partition. A caller that also needs
/// [`member_alt_spans`] on the same DAG should instead compute the
/// `Dominators` once and call [`partition_members_with`] and
/// [`member_alt_spans`] directly, to avoid sweeping the DAG twice.
// The only item in this module with no production caller: `normalize_allele`'s
// shadow comparison needs `member_alt_spans` on the same DAG, so it takes the
// `partition_members_with` route that shares one `Dominators` sweep. This
// convenience wrapper is used by tests and is kept as the documented one-off
// entry point. Scoped rather than allowed module-wide, so a helper that becomes
// unreachable when the shadow is promoted still warns.
#[allow(dead_code)]
pub(crate) fn partition_members(dag: &AlignmentDag) -> Vec<MemberBlock> {
    let dominators = dag.dominators();
    partition_members_with(dag, &dominators, MIN_SEPARATION)
}

/// Partition the reference block, merging runs separated by fewer than
/// `min_separation` **unchanged reference bases**.
///
/// Takes `dominators` rather than computing it, so a caller that also needs
/// [`member_alt_spans`] can sweep the DAG once and pass the same `Dominators`
/// to both.
///
/// The unit of separation is unchanged bases, not event indices: two events are
/// separated by the count of reference bases lying strictly between the extents
/// they occupy. Comparing event indices instead splits the spec's own worked
/// example (`LRG_199t1:c.850_901`) into two members where the spec describes
/// one; the other calibration cases score identically under both readings, so
/// that error is invisible without that case.
///
/// The left extent is therefore *exclusive* and depends on what the event is: a
/// changed position `p` occupies `p..p + 1`, so an event at `q` is `q - p - 1`
/// bases away, but an insertion junction at `k` occupies no reference base at
/// all, so an event at `q` is `q - k` bases away. Subtracting one
/// unconditionally undercounts after a junction and merges runs that a
/// separation of one unchanged base should keep apart.
pub(crate) fn partition_members_with(
    dag: &AlignmentDag,
    dominators: &Dominators,
    min_separation: u32,
) -> Vec<MemberBlock> {
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

    // Where an event's extent ends, exclusive: a changed reference position
    // occupies `end..end + 1`, while a forced-insertion junction occupies no
    // reference base at all and so ends at `end`.
    //
    // Widening an insertion-only end to `end + 1` would claim a reference base
    // that every minimal alignment matches, giving the member a non-minimal
    // extent. Both readings denote the same sequence and round-trip, so the
    // round-trip invariant does not catch this. The exact-span assertions in
    // `insertion_junctions_do_not_claim_a_matched_base` and
    // `an_insertion_only_member_has_an_empty_reference_span` do.
    //
    // When a position is both changed and a forced-insertion junction, the
    // changed reading wins — it is the wider of the two and the base really
    // does change.
    let exclusive_end = |end: u32| -> u32 {
        if changed.binary_search(&end).is_ok() {
            // `end` is itself a changed position, so `end < ref_len` always (a
            // changed position is by definition within the block); the `+ 1`
            // can never reach past `ref_len`, so no clamp is needed.
            debug_assert!(end < ref_len, "a changed position must be within the block");
            end + 1
        } else {
            end
        }
    };

    // Group consecutive events into runs, merging across short gaps.
    let mut runs: Vec<(u32, u32)> = Vec::new();
    for event in events {
        match runs.last_mut() {
            // The unchanged bases between the run so far and this event are the
            // ones from where the run's extent ends up to this event's start.
            Some((_, last_end))
                if event.saturating_sub(exclusive_end(*last_end)) < min_separation =>
            {
                *last_end = event;
            }
            _ => runs.push((event, event)),
        }
    }

    runs.into_iter()
        .map(|(start, end)| MemberBlock {
            ref_start: start,
            ref_end: exclusive_end(end),
        })
        .collect()
}

/// The alternate-block span each member consumes, parallel to `members`.
///
/// Takes `dominators` rather than computing it, so a caller that also calls
/// [`partition_members_with`] on the same DAG can sweep it once and pass the
/// same `Dominators` to both — `dominators` must be `dag.dominators()` (or
/// equal to it), since the two are read together below.
///
/// Each member's *start* is reached across a run of reference bases that every
/// minimal alignment matches, so the alternate cursor advances one-for-one
/// there. Each member's *end* is read from the dominator that terminates it.
///
/// Returns `None` in any of three cases, each meaning a bug upstream rather
/// than an unusual input:
///
/// 1. a member's end is neither the end of the block nor a dominator match;
/// 2. the members are not an ascending, non-overlapping partition (a member's
///    `ref_start` precedes the previous member's `ref_end`, or a member's own
///    `ref_end` precedes its `ref_start`);
/// 3. the computed span would be inverted (`alt_end < alt_start`) — this can
///    only arise alongside a malformed `members` input, since a well-formed
///    partition's alternate cursor never runs backward.
///
/// # Why the end must come from a dominator
///
/// It is tempting to walk one minimal alignment and record the alternate
/// coordinate at each reference boundary. That is wrong: away from dominators,
/// different minimal alignments assign different alternate coordinates to the
/// same reference boundary, so a walked path yields spans that rebuild the
/// wrong alternate block. `AAAAAAA -> AACACAAAA` and `ACGTACGTAC ->
/// ACGTTACGTGAC` both fail the round-trip that way, and both pass when the end
/// is taken from the dominator's own cell.
///
/// A member's `ref_end` is always either `ref_len` or a dominator-matched
/// position: a run ending at a changed position `p` has `ref_end = p + 1`, and
/// `p + 1` cannot itself be changed or it would be in the same run; a run
/// ending at an insertion junction `k` has `ref_end = k`, and `k` is matched or
/// it would have been a changed position.
pub(crate) fn member_alt_spans(
    dag: &AlignmentDag,
    dominators: &Dominators,
    members: &[MemberBlock],
) -> Option<Vec<(u32, u32)>> {
    let alt_at_boundary = |r: u32| -> Option<u32> {
        if r == dag.ref_len() {
            Some(dag.alt_len())
        } else {
            dominators.alt_at(r)
        }
    };

    let mut spans = Vec::with_capacity(members.len());
    let mut prev_ref_end = 0u32;
    let mut prev_alt_end = 0u32;
    for member in members {
        if member.ref_start < prev_ref_end || member.ref_end < member.ref_start {
            return None;
        }
        // The matched run between the previous member and this one advances
        // both cursors equally.
        let alt_start = prev_alt_end + (member.ref_start - prev_ref_end);
        let alt_end = alt_at_boundary(member.ref_end)?;
        if alt_end < alt_start {
            return None;
        }
        spans.push((alt_start, alt_end));
        prev_ref_end = member.ref_end;
        prev_alt_end = alt_end;
    }
    Some(spans)
}

/// One maximal run of change along a single alignment: the reference and
/// alternate spans it consumes, both half-open and both relative to the start of
/// their block.
///
/// A pure insertion has `ref_start == ref_end`; a pure deletion has
/// `alt_start == alt_end`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CanonicalRun {
    ref_start: u32,
    ref_end: u32,
    alt_start: u32,
    alt_end: u32,
}

/// Whether a run is open at this point of the walk.
///
/// The cost of a change step depends on it — the first change after a match
/// opens a new run and costs one member, every change after that is free — so it
/// is part of the shortest-path state rather than something recoverable from the
/// cell alone.
const RUN_CLOSED: usize = 0;
/// A run is open: the next change step extends it rather than opening another.
const RUN_OPEN: usize = 1;

/// The **canonical** alignment: among the minimal alignments the DAG encodes,
/// one that splits the block into the fewest members.
///
/// This is the second of the two rules `mutalyzer-algebra` distinguishes, and
/// the one ferro did not have. [`partition_members_with`] implements the
/// *dominator* rule — cut at the steps common to **every** minimal alignment,
/// which is algebra's `local_supremal`. This implements algebra's `canonical`:
///
/// > Traverse (BFS) the LCS graph to extract the canonical variant, i.e.
/// > minimize the number of separate variants within an allele.
///
/// The two answer different questions and neither subsumes the other. The
/// dominator rule describes the *union* of what every minimal alignment changes,
/// so it is independent of which alignment you pick; this rule picks one
/// alignment and describes exactly what it changes.
///
/// # Why minimising members is a second key and not a competing one
///
/// Every path through [`AlignmentDag`] is distance-minimal by construction —
/// that is what the DAG encodes — so ranking its paths by member count cannot
/// trade edit distance away. Minimal distance stays the primary key; member
/// count is applied strictly within it.
///
/// # Cost
///
/// One reverse topological sweep plus one forward walk, both `Θ(n·m)` — the same
/// order as [`AlignmentDag::build`] and [`AlignmentDag::dominators`]. The cost of
/// a *path* is the number of maximal change runs on it, not its edge count, so
/// the sweep carries the two-valued run state described above rather than a
/// plain per-cell distance.
///
/// # The tie-break is an implementer's choice
///
/// Member-count-minimal alignments are **not** unique. `AAC -> AC` deletes
/// either of the two `A`s for one member either way, and nothing in the HGVS
/// spec reaches the choice. The **3'-most** alignment is taken, matching
/// `general.md:41` ("the most 3' position possible of the reference sequence is
/// arbitrarily assigned to have been changed") and the tie-break
/// [`AlignmentDag`]'s live counterpart already applies.
///
/// It is realised as a preference order on the forward walk: at every cell, among
/// the out-edges that preserve the minimum, take `Match` first. Preferring the
/// match defers change for as long as the alignment allows, which *is* the 3'
/// placement — on `AAC -> AC` it matches reference offset 0 and deletes offset 1.
///
/// Ties that survive that — a cell whose only optimal steps are two different
/// *kinds* of change — are broken by [`AlignmentDag::out_edges`]'s own order
/// (`Match`, `Sub`, `Del`, `Ins`). That residual order is arbitrary and is
/// deliberately not dressed up as a rule: a run's ref and alt spans are fixed by
/// where it starts and ends, and the order of steps *within* it does not move
/// either. It is stated only so the function is deterministic by construction
/// rather than by accident.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CanonicalAlignment {
    runs: Vec<CanonicalRun>,
}

impl CanonicalAlignment {
    /// Select the canonical alignment of `dag`.
    ///
    /// `dag`'s blocks must have had their common flanks trimmed, for the same
    /// reason [`AlignmentDag::build`] states: an untrimmed flank base is
    /// absorbed into a neighbouring member and changes the answer.
    pub(crate) fn of(dag: &AlignmentDag) -> Self {
        let ref_len = dag.ref_len();
        let alt_len = dag.alt_len();
        let width = alt_len as usize + 1;
        let size = (ref_len as usize + 1) * width;
        let index = |i: u32, j: u32| i as usize * width + j as usize;
        let sink = (ref_len, alt_len);

        // `to_go[state * size + cell]` is the fewest additional members any
        // minimal alignment can reach the sink with, from that cell and that run
        // state. `u32::MAX` marks "not yet computed"; every cell in the DAG
        // reaches the sink, so no cell keeps it.
        let mut to_go = vec![u32::MAX; 2 * size];
        to_go[RUN_CLOSED * size + index(sink.0, sink.1)] = 0;
        to_go[RUN_OPEN * size + index(sink.0, sink.1)] = 0;

        // `cells()` is topological (ascending `i + j`), and every edge strictly
        // increases `i + j`, so walking it backwards visits every successor
        // before its predecessor.
        let cells: Vec<(u32, u32)> = dag.cells().collect();
        for &(i, j) in cells.iter().rev() {
            if (i, j) == sink {
                continue;
            }
            for state in [RUN_CLOSED, RUN_OPEN] {
                let mut best = u32::MAX;
                for (next_i, next_j, step) in dag.out_edges(i, j) {
                    let (cost, next_state) = step_cost(step, state);
                    let onward = to_go[next_state * size + index(next_i, next_j)];
                    if onward != u32::MAX {
                        best = best.min(onward + cost);
                    }
                }
                to_go[state * size + index(i, j)] = best;
            }
        }

        // Walk forward, taking at each cell the first out-edge that preserves
        // the minimum. `out_edges` yields `Match` before any change step, so
        // "first" is the 3'-most preference documented above.
        let mut runs: Vec<CanonicalRun> = Vec::new();
        let (mut i, mut j) = (0u32, 0u32);
        let mut state = RUN_CLOSED;
        while (i, j) != sink {
            let target = to_go[state * size + index(i, j)];
            let chosen = dag.out_edges(i, j).find(|&(next_i, next_j, step)| {
                let (cost, next_state) = step_cost(step, state);
                let onward = to_go[next_state * size + index(next_i, next_j)];
                onward != u32::MAX && onward + cost == target
            });
            // Unreachable: `to_go` was computed from these same edges, so some
            // edge attains it. Declining beats panicking inside a library, and a
            // truncated answer cannot escape: `canonicalize_from_sequence`
            // re-applies the rebuilt edits and compares them against `result`
            // (`merge.rs`, the `if reapplied != result { return None; }` just
            // before it returns `Some(rebuilt)`). That is a *runtime* refusal
            // carried into release builds, not the `debug_assert_eq!` beside it,
            // so a partial partition declines rather than being emitted.
            let Some((next_i, next_j, step)) = chosen else {
                debug_assert!(false, "no out-edge of ({i},{j}) attains its own optimum");
                break;
            };
            if step == Step::Match {
                state = RUN_CLOSED;
            } else {
                if state == RUN_CLOSED {
                    runs.push(CanonicalRun {
                        ref_start: i,
                        ref_end: i,
                        alt_start: j,
                        alt_end: j,
                    });
                }
                if let Some(run) = runs.last_mut() {
                    run.ref_end = next_i;
                    run.alt_end = next_j;
                }
                state = RUN_OPEN;
            }
            (i, j) = (next_i, next_j);
        }

        Self { runs }
    }

    /// The member blocks of the canonical alignment, ascending and disjoint.
    pub(crate) fn members(&self) -> Vec<MemberBlock> {
        self.runs
            .iter()
            .map(|run| MemberBlock {
                ref_start: run.ref_start,
                ref_end: run.ref_end,
            })
            .collect()
    }

    /// The alternate-block span each member consumes, parallel to
    /// [`CanonicalAlignment::members`].
    ///
    /// Unlike [`member_alt_spans`] this needs no dominator lookup and cannot
    /// decline: the spans are read straight off the chosen path, which pins one
    /// alternate coordinate per boundary by construction.
    pub(crate) fn alt_spans(&self) -> Vec<(u32, u32)> {
        self.runs
            .iter()
            .map(|run| (run.alt_start, run.alt_end))
            .collect()
    }
}

/// The member cost of taking `step` with a run in `state`, and the state it
/// leaves behind.
///
/// A match closes any open run and costs nothing. A change step costs one member
/// when it opens a run and nothing when it extends one — which is what makes the
/// path cost the number of maximal change runs rather than the number of edges.
fn step_cost(step: Step, state: usize) -> (u32, usize) {
    match step {
        Step::Match => (0, RUN_CLOSED),
        Step::Sub | Step::Del | Step::Ins => (u32::from(state == RUN_CLOSED), RUN_OPEN),
    }
}

/// Partition the reference block by the **canonical** rule: the member blocks of
/// a member-count-minimal minimal alignment.
///
/// The counterpart of [`partition_members_with`] for the other of the two rules
/// — see [`CanonicalAlignment`] for what distinguishes them, and for the
/// tie-break this makes among equally-few-member alignments.
///
/// It takes no `min_separation`, and that absence is deliberate. Merging runs
/// across a short unchanged gap **widens** a member past what the alignment
/// changes, so a merged partition claims more changed columns than the block's
/// edit distance: #1260's `AAAAAAA -> AACACAAAA` is distance 2, and merging its
/// two insertion runs across the single matched base between them yields one
/// member spanning `2..3` with a 3-base payload — 3 changed columns for a
/// distance-2 block. Distance-minimality is the property this rule exists to
/// have, so the separation rule stays where it belongs, on the dominator side.
// No caller outside this module's own tests. `partition_block_canonical` does
// not route through here — it calls `CanonicalAlignment::of` and then
// `members()`/`alt_spans()` directly, because it needs both halves and this
// returns only the first — and `dump_partitions` reaches the canonical rule
// through `dev_partitioners`, which wraps `partition_block_canonical`. So the
// sole call site is `round_trips_exhaustively_over_a_small_alphabet` below.
// Scoped rather than allowed module-wide, so a helper that becomes unreachable
// still warns.
#[allow(dead_code)]
pub(crate) fn partition_members_canonical(dag: &AlignmentDag) -> Vec<MemberBlock> {
    CanonicalAlignment::of(dag).members()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The calibration corpus. Six cases; 6/6 is the bar this design cleared
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

    /// An insertion junction occupies no reference base, so the run it ends
    /// stops *at* that offset rather than after it, and the base at that offset
    /// is still unchanged and still separates it from whatever follows.
    ///
    /// #1260's block has forced insertion junctions at reference offsets 2 and
    /// 3, with the base at offset 2 unchanged between them: the separation is
    /// `3 - 2 = 1` unchanged base, not `3 - 2 - 1 = 0`.
    ///
    /// At `MIN_SEPARATION` (2) both readings merge, which is why neither the
    /// exhaustive round-trip sweep nor the proptest catches the difference —
    /// the members denote the same sequence either way. Only a threshold of 1
    /// separates them, and there the undercount is the difference between
    /// merging and splitting.
    #[test]
    fn separation_after_an_insertion_junction_counts_the_base_it_does_not_occupy() {
        let dag = AlignmentDag::build(b"AAAAAAA", b"AACACAAAA");
        let dominators = dag.dominators();
        assert_eq!(
            partition_members_with(&dag, &dominators, 1).len(),
            2,
            "one unchanged base is not fewer than a separation threshold of 1"
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

    #[test]
    fn lrg199_partitions_to_the_exact_whole_block_span() {
        // The design's headline case (`delins.md:44`, `LRG_199t1:c.850_901`) is
        // otherwise pinned only by member *count* (the corpus row above). Only
        // reference offset 48 is a dominator match, so the two runs either
        // side of it merge across the single unchanged base between them
        // (fewer than `MIN_SEPARATION`), giving one member spanning the whole
        // 52 nt block.
        let dag = AlignmentDag::build(LRG199_REF, b"TTCCTCGATGCCTG");
        let members = partition_members(&dag);
        assert_eq!(
            members,
            vec![MemberBlock {
                ref_start: 0,
                ref_end: 52
            }]
        );
    }

    #[test]
    fn min_separation_changes_whether_runs_merge() {
        // "control two subs exactly 2 apart": at the default MIN_SEPARATION
        // (2), the gap of exactly 2 unchanged bases does not merge, giving two
        // members. Raising `min_separation` to 3 must merge them into one,
        // pinning that the parameter is actually load-bearing rather than
        // dead flexibility.
        let dag = AlignmentDag::build(b"ACGTACGT", b"ATGTTCGT");
        let dominators = dag.dominators();
        assert_eq!(
            partition_members_with(&dag, &dominators, 2).len(),
            2,
            "gap of exactly 2 must not merge at min_separation 2"
        );
        assert_eq!(
            partition_members_with(&dag, &dominators, 3).len(),
            1,
            "gap of exactly 2 must merge at min_separation 3"
        );
    }

    #[test]
    fn member_alt_spans_matches_measured_values() {
        // Verified exact spans, not just member boundaries: #1260's single
        // member and both members of the "two insertions 5 apart" case.
        for (name, r, a, expected_spans) in [
            (
                "#1260 two insC at adjacent gaps",
                &b"AAAAAAA"[..],
                &b"AACACAAAA"[..],
                vec![(2u32, 5u32)],
            ),
            (
                "two insertions 5 apart",
                &b"ACGTACGTAC"[..],
                &b"ACGTTACGTGAC"[..],
                vec![(3u32, 5u32), (9u32, 10u32)],
            ),
        ] {
            let dag = AlignmentDag::build(r, a);
            let dominators = dag.dominators();
            let members = partition_members_with(&dag, &dominators, MIN_SEPARATION);
            let spans = member_alt_spans(&dag, &dominators, &members)
                .unwrap_or_else(|| panic!("{name}: member_alt_spans returned None"));
            assert_eq!(spans, expected_spans, "{name}: alt spans");
        }
    }

    #[test]
    fn member_alt_spans_rejects_overlapping_members() {
        // `member_alt_spans` takes an arbitrary `&[MemberBlock]`, not only ones
        // `partition_members` could have produced, and its contract promises
        // `None` for a members list that is not an ascending, non-overlapping
        // partition. Exercise that directly rather than trying to coax an
        // impossible members list out of `partition_members` itself.
        let dag = AlignmentDag::build(b"ACGT", b"ACGT");
        let dominators = dag.dominators();
        let overlapping = [
            MemberBlock {
                ref_start: 0,
                ref_end: 2,
            },
            MemberBlock {
                ref_start: 1,
                ref_end: 3,
            },
        ];
        assert_eq!(member_alt_spans(&dag, &dominators, &overlapping), None);
    }

    /// Rebuild the alternate block from the reference block plus the members'
    /// payloads. This is the invariant that makes the partition meaningful: if
    /// it does not hold, the members denote a different sequence.
    fn round_trips(reference: &[u8], alt: &[u8]) -> bool {
        let dag = AlignmentDag::build(reference, alt);
        let dominators = dag.dominators();
        let members = partition_members_with(&dag, &dominators, MIN_SEPARATION);
        let Some(alt_spans) = member_alt_spans(&dag, &dominators, &members) else {
            return false;
        };
        let mut rebuilt = Vec::new();
        let mut ref_cursor = 0usize;
        for (member, (alt_start, alt_end)) in members.iter().zip(&alt_spans) {
            rebuilt.extend_from_slice(&reference[ref_cursor..member.ref_start as usize]);
            rebuilt.extend_from_slice(&alt[*alt_start as usize..*alt_end as usize]);
            ref_cursor = member.ref_end as usize;
        }
        rebuilt.extend_from_slice(&reference[ref_cursor..]);
        rebuilt == alt
    }

    #[test]
    fn calibration_corpus_round_trips() {
        for (r, a) in [
            (&b"AAAAAAA"[..], &b"AACACAAAA"[..]),
            (&b"AAAAAAA"[..], &b"ACAAAA"[..]),
            (LRG199_REF, &b"TTCCTCGATGCCTG"[..]),
            (&b"ACGTACGTACGTACGT"[..], &b"ATGTACGTACGTACTT"[..]),
            (&b"ACGTACGT"[..], &b"ATGTTCGT"[..]),
            (&b"ACGTACGTAC"[..], &b"ACGTTACGTGAC"[..]),
            (&b"GACTGACTGA"[..], &b"GAATAACTGA"[..]),
            (&b"ACGT"[..], &b"ACGT"[..]),
            (&b""[..], &b"ACGT"[..]),
            (&b"ACGT"[..], &b""[..]),
        ] {
            assert!(
                round_trips(r, a),
                "{} -> {} does not round-trip",
                String::from_utf8_lossy(r),
                String::from_utf8_lossy(a)
            );
        }
    }

    #[test]
    fn round_trips_exhaustively_over_a_small_alphabet() {
        let all = small_alphabet_words(5);
        for r in &all {
            for a in &all {
                assert!(
                    round_trips(r, a),
                    "{} -> {} does not round-trip",
                    String::from_utf8_lossy(r),
                    String::from_utf8_lossy(a)
                );
            }
        }
    }

    /// Every word over `{A, C}` up to `max_len`, shortest first.
    ///
    /// A two-letter alphabet and short blocks maximise the number of *distinct*
    /// minimal alignments per block, which is the regime the canonical rule is
    /// about — over four letters most blocks have one minimal alignment and the
    /// selection has nothing to select.
    fn small_alphabet_words(max_len: u32) -> Vec<Vec<u8>> {
        fn words(len: u32) -> Vec<Vec<u8>> {
            if len == 0 {
                return vec![Vec::new()];
            }
            let mut out = Vec::new();
            for shorter in words(len - 1) {
                for base in *b"AC" {
                    let mut w = shorter.clone();
                    w.push(base);
                    out.push(w);
                }
            }
            out
        }
        (0..=max_len).flat_map(words).collect()
    }

    /// Changed columns a member set claims: `max(ref_span, alt_len)` summed,
    /// mirroring `merge::changed_columns_of_pieces`.
    fn changed_columns(members: &[MemberBlock], alt_spans: &[(u32, u32)]) -> u32 {
        members
            .iter()
            .zip(alt_spans)
            .map(|(member, (alt_start, alt_end))| {
                (member.ref_end - member.ref_start).max(alt_end - alt_start)
            })
            .sum()
    }

    /// The canonical rule may not describe more change than the block contains.
    ///
    /// This is the property that separates it from the dominator rule and the
    /// reason it takes no `min_separation`: its members come from **one**
    /// alignment, so their changed columns must total exactly the block's
    /// Levenshtein distance. Anything that merges or widens a member breaks this
    /// — see `partition_members_canonical`'s doc for the worked #1260 case that
    /// a separation threshold would push to 3 columns on a distance-2 block.
    ///
    /// Exhaustive over `{A, C}` up to length 6 both sides (16 129 pairs): a
    /// two-letter alphabet is where alternative minimal alignments are dense, so
    /// a selection rule that silently picks a non-minimal path shows up here and
    /// nowhere cheaper.
    #[test]
    fn canonical_members_claim_exactly_the_blocks_edit_distance() {
        let all = small_alphabet_words(6);
        let mut checked = 0usize;
        for r in &all {
            for a in &all {
                let dag = AlignmentDag::build(r, a);
                let canonical = CanonicalAlignment::of(&dag);
                assert_eq!(
                    changed_columns(&canonical.members(), &canonical.alt_spans()),
                    dag.edit_distance(),
                    "{} -> {} claims a non-minimal number of changed columns",
                    String::from_utf8_lossy(r),
                    String::from_utf8_lossy(a)
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 127 * 127, "expected the full sweep, ran {checked}");
    }

    /// The members must denote the block they were derived from.
    ///
    /// Distance-minimality above says the members are *small enough*; this says
    /// they are *right*. Without it a rule that dropped a run entirely would
    /// score better on changed columns, not worse.
    #[test]
    fn canonical_members_rebuild_the_alternate_block() {
        let all = small_alphabet_words(6);
        for r in &all {
            for a in &all {
                let dag = AlignmentDag::build(r, a);
                let canonical = CanonicalAlignment::of(&dag);
                let members = canonical.members();
                let spans = canonical.alt_spans();
                let mut rebuilt = Vec::new();
                let mut ref_cursor = 0usize;
                for (member, (alt_start, alt_end)) in members.iter().zip(&spans) {
                    assert!(
                        member.ref_start as usize >= ref_cursor,
                        "members overlap or are out of order: {members:?}"
                    );
                    rebuilt.extend_from_slice(&r[ref_cursor..member.ref_start as usize]);
                    rebuilt.extend_from_slice(&a[*alt_start as usize..*alt_end as usize]);
                    ref_cursor = member.ref_end as usize;
                }
                rebuilt.extend_from_slice(&r[ref_cursor..]);
                assert_eq!(
                    rebuilt,
                    *a,
                    "{} -> {} does not round-trip",
                    String::from_utf8_lossy(r),
                    String::from_utf8_lossy(a)
                );
            }
        }
    }

    /// The tie-break, pinned on the smallest block that has one.
    ///
    /// `AAC -> AC` deletes either `A` for one member either way, so member count
    /// cannot choose and `general.md:41`'s 3' rule does: the deleted base is
    /// reference offset **1**, not 0. Pinned because the tie-break is an
    /// implementer's choice — nothing in the spec forces it — so a change to it
    /// is a change of representation and must be a deliberate edit rather than a
    /// side effect of touching the walk order.
    #[test]
    fn canonical_breaks_a_member_count_tie_three_prime_most() {
        let dag = AlignmentDag::build(b"AAC", b"AC");
        assert_eq!(
            partition_members_canonical(&dag),
            vec![MemberBlock {
                ref_start: 1,
                ref_end: 2
            }]
        );
    }

    /// The canonical rule can report **more** members than the dominator rule,
    /// and this pins the smallest block where it does.
    ///
    /// It is tempting to assume that minimising members bounds it below the
    /// dominator partition. It does not, and the two rules are not comparable at
    /// all: the dominator rule merges runs separated by fewer than
    /// [`MIN_SEPARATION`] unchanged bases, which the canonical rule has no
    /// counterpart for and deliberately does not want — merging is what would
    /// cost it distance-minimality.
    ///
    /// `A -> CAC` is #1260's own block, trimmed, and it is the shortest block in
    /// the `{A, C}` sweep where they diverge. It is distance 2 (one `C` inserted
    /// either side of the `A`) and admits exactly one minimal alignment, so this
    /// is not a tie-break disagreement: the dominator rule merges the two forced
    /// insertion junctions across the single matched base between them into one
    /// member claiming **3** changed columns, while the canonical rule reports
    /// them as two pure insertions claiming **2**.
    ///
    /// Recorded rather than asserted away because it is the whole reason both
    /// rules are worth measuring. Neither is a better approximation of the other:
    /// one is minimal in members-after-merging and the other in changed columns,
    /// and the spec does not settle which a description should be.
    #[test]
    fn canonical_can_report_more_members_than_the_dominator_rule() {
        let dag = AlignmentDag::build(b"A", b"CAC");
        let dominators = dag.dominators();
        assert_eq!(
            partition_members_with(&dag, &dominators, MIN_SEPARATION),
            vec![MemberBlock {
                ref_start: 0,
                ref_end: 1
            }],
            "the dominator rule merges both insertions into one spanning member"
        );
        let canonical = CanonicalAlignment::of(&dag);
        assert_eq!(
            canonical.members(),
            vec![
                MemberBlock {
                    ref_start: 0,
                    ref_end: 0
                },
                MemberBlock {
                    ref_start: 1,
                    ref_end: 1
                },
            ],
            "the canonical rule keeps them as two pure insertions"
        );
        assert_eq!(
            changed_columns(&canonical.members(), &canonical.alt_spans()),
            dag.edit_distance(),
            "and only the canonical partition claims exactly the block's distance"
        );
    }

    proptest::proptest! {
        /// The same invariant on random four-letter sequences.
        #[test]
        fn round_trips_on_random_blocks(
            reference in proptest::collection::vec(
                proptest::sample::select(vec![b'A', b'C', b'G', b'T']), 0..40usize),
            alt in proptest::collection::vec(
                proptest::sample::select(vec![b'A', b'C', b'G', b'T']), 0..40usize),
        ) {
            proptest::prop_assert!(round_trips(&reference, &alt));
        }
    }
}
