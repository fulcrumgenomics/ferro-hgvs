//! An insertion or duplication must not carry its junction over a sibling.
//!
//! The sibling-crossing clamp governs edits that consume the bases under their
//! span. An insertion or duplication consumes none — it adds sequence at the
//! junction between two bases — so it is excluded there, deliberately: a member
//! landing flush against one is the adjacency the collapse pass exists to catch
//! (#999, #1135). A duplication is the one partial exception, and only in the
//! 5' direction: the bases under its span are what it copies, so a sibling
//! shifting 5' onto them is bounded (`blocks_sibling_shift`, #1286). Reaching
//! it from the 3' side is still the adjacency above.
//!
//! It can still cross. The junction 3'-shifts through a tract like anything
//! else, and moving it past a base a sibling edits changes what the allele
//! denotes:
//!
//! ```text
//! reference    ACAAAAAAAACGTACGTACG        A-run at 3-10
//! g.[2_3insA;5A>G]  ->  g.[5A>G;10dup]     junction moved 2|3 -> 10|11
//! input applied     ACAAAGAAAAACGTACGTACG
//! output applied    ACAAGAAAAAACGTACGTACG  ← the substituted base moved
//! ```
//!
//! Well-formed, disjoint, warning-free, and wrong — the silent shape.
//!
//! The same defect reaches the repeat path. A `dup` or `ins` inside a tandem
//! tract canonicalises to a copy count over the **whole tract**, and that span
//! can swallow a sibling outright:
//!
//! ```text
//! reference    TTTTTTTTTAATATATTTTA
//! g.[1_2dup;4T>A]  ->  g.[1_9T[11];4T>A]   `4T>A` sits inside `1_9`
//! ```
//!
//! The demotion pass re-spelled only repeats that grew from a *deletion*, so
//! these were left. Both are fixed here: the demotion covers `dup`/`ins`-grown
//! repeats too, and the junction is then bounded at the sibling's 5' edge —
//! still flush against it, so the #999 collapse keeps firing.

use crate::common::cis_apply_oracle::{
    apply, assert_normalizes_preserving, assert_normalizes_preserving_in, normalize, normalize_in,
    sweep_seeds, sweep_sequences,
};
use ferro_hgvs::ShuffleDirection;

/// Nine `T` at positions 1-9 — a tract long enough to canonicalise to a repeat.
const TRACT: &str = "TTTTTTTTTAATATATTTTA";

/// An `A`-run at positions 3-10 — long enough for a junction to travel, short
/// enough that a two-copy insertion does not reach repeat notation.
const RUN: &str = "ACAAAAAAAACGTACGTACG";

/// Sequence-changing cases the exhaustive sweep below finds in the `dup` +
/// `del` + 5'-shuffle shape. **Zero**, and pinned as an exact count so the
/// shape stays measured separately from the rest of the sweep.
///
/// This was 74 — the residual tracked as #1266 — until `blocks_sibling_shift`
/// made a duplication a barrier to a sibling's shift *and* made a duplication's
/// own shift subject to the clamp. A rise means a new sequence-changing defect
/// in that shape specifically.
const FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES: usize = 0;

// The 5'-shuffle insertion-junction shape — the insertion half of #1267 — used
// to be excluded from this file's sequence-preservation assertion and counted
// instead, as `TRACT_FIVE_PRIME_INSERTION_SEQUENCE_CHANGES_BY_POSITION =
// [(15, 100), (36, 16)]`: 116 cases pinned per sibling position so the shape
// stayed measured while it was broken.
//
// It is fixed, so the exclusion is gone and those cases are asserted like every
// other — a regression there now reports the offending descriptions rather than
// a bare count. The two halves of that 116 needed separate bounds, and measuring
// per-position is what showed the second half existed: bounding the junction
// against base-claiming siblings alone took the map to `{15: 50, 36: 8}`,
// exactly half, because the other half's sibling was a `dup` or an `ins` —
// junction-occupying, so invisible to `claims_bases`. Mirroring the
// junction-vs-junction bound (#1290) closed the rest.

#[test]
fn an_insertion_does_not_shift_past_a_substituted_base() {
    // The junction may travel the `A`-run only as far as the base before the
    // substitution; clamped there it is flush against it, and the two coalesce.
    assert_normalizes_preserving(RUN, "TEMPLATE:g.[2_3insA;5A>G]", "TEMPLATE:g.5_6insG");
}

#[test]
fn a_duplication_does_not_shift_past_a_substituted_base() {
    assert_normalizes_preserving(RUN, "TEMPLATE:g.[4dup;7A>G]", "TEMPLATE:g.7_8insG");
}

#[test]
fn a_repeat_grown_from_a_duplication_does_not_span_a_sibling() {
    // `1_2dup` canonicalised to `1_9T[11]`, whose span covers the substituted
    // base at 4. Demoted back to a duplication, clamped, then coalesced.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2dup;4T>A]", "TEMPLATE:g.5_6insAT");
}

#[test]
fn a_repeat_grown_from_a_duplication_does_not_span_a_deletion() {
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1_2dup;4del]", "TEMPLATE:g.9dup");
}

#[test]
fn a_repeat_grown_from_an_insertion_does_not_span_a_sibling() {
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[2_3insTT;5T>A]", "TEMPLATE:g.6_7insAT");
}

#[test]
fn a_lone_duplication_still_reaches_repeat_notation() {
    // Negative control: no sibling, nothing to span, so the tract-wide repeat
    // is exactly right and must survive.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.1_2dup", "TEMPLATE:g.1_9T[11]");
}

#[test]
fn an_insertion_may_still_land_flush_against_its_sibling() {
    // Negative control for #999, the rule this clamp must not break. The
    // inserted `C` 3'-shifts to a dup at 306, landing at the junction
    // `306|307` — the last position the clamp permits, since `307G>T` claims
    // base 307 and the junction has not passed it. The two then coalesce.
    let mut seq = vec![b'A'; 1000];
    for (i, b) in "CATCCTCGCTCCT".bytes().enumerate() {
        seq[299 + i] = b;
    }
    let seq = String::from_utf8(seq).unwrap();
    assert_eq!(
        normalize(&seq, "TEMPLATE:g.[305_306insC;307G>T]"),
        "TEMPLATE:g.307delinsCT"
    );
}

#[test]
fn an_insertion_reaching_the_end_of_its_run_still_collapses() {
    // Negative control: the junction travels the whole `A`-run and lands flush
    // against `11C>G`, which is outside the run — no crossing, so the two
    // coalesce exactly as before.
    assert_normalizes_preserving(RUN, "TEMPLATE:g.[2_3insA;11C>G]", "TEMPLATE:g.11delinsAG");
}

/// A 5'-shuffled insertion merging with a deletion sibling must not move the
/// base it contributes (#1342).
///
/// Surfaced by giving the sweep's first-member `ins` a payload that is not the
/// span's own reference base. In the `T`-run at 1-9, inserting an `A` at the
/// `2|3` junction and deleting a `T` at 4 is net "one `T` becomes an `A`" — and
/// the `A` lands at position 3. It cannot shuffle 5': an `A` in a `T`-run is
/// unique, so moving it is a different sequence, not a different spelling of
/// the same one. The 3' direction already answered `g.3T>A`.
#[test]
fn an_insertion_merging_with_a_deletion_keeps_its_base_in_place() {
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insA;4del]",
        "TEMPLATE:g.3T>A",
        ShuffleDirection::FivePrime,
    );
}

/// Negative controls for the 5' junction barrier: it must bound only the shifts
/// that actually reorder the two edits, never a member that may legitimately
/// travel to the 5' end of its tract.
///
/// This is the failure mode the barrier's own assertion cannot see. A member
/// held back too far still *preserves sequence*, so the exhaustive sweep below
/// would pass while every answer quietly under-shifted. These pin the shifts
/// that must still complete.
#[test]
fn the_five_prime_junction_barrier_does_not_over_clamp() {
    // Payload equals the run base, so the insertion and the deletion cancel:
    // the pair must still reduce all the way to an identity rather than being
    // stranded as two members either side of the junction.
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insT;4del]",
        "TEMPLATE:g.2_3=",
        ShuffleDirection::FivePrime,
    );
    // Net +1 `T` in the run. The 5'-most spelling is a duplication at 1, which
    // is 5' of the `2|3` junction — reachable because the merged member no
    // longer straddles it.
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insTT;5del]",
        "TEMPLATE:g.1dup",
        ShuffleDirection::FivePrime,
    );
    // The deletion is outside the tract entirely, so there is no sweep and
    // nothing for the barrier to bound — both members stay put.
    assert_normalizes_preserving_in(
        TRACT,
        "TEMPLATE:g.[2_3insA;12del]",
        "TEMPLATE:g.[2_3insA;12del]",
        ShuffleDirection::FivePrime,
    );
}

/// The tract `[4_5insC;5_6dup]` shuffles through: an `AT`-alternating run at
/// 3-11, long enough for a two-base duplication to travel and short enough that
/// a third member placed past 11 sits clear of it.
const DUP_RUN: &str = "TTAATATATAATAATATTAT";

/// The pair alone. Both members land inside one run of change, so #1235's
/// sequence-first derivation answers the whole allele from the bases it denotes
/// and emits a single insertion.
///
/// **Re-blessed from `TEMPLATE:g.[4_5insC;4_5dup]`.** The merged form is derived
/// from the sequence rather than assembled per member; `[4_5insC;4_5dup]` was
/// only ever the per-member repair's spelling of the same bases. Verified
/// sequence-preserving: `assert_normalizes_preserving_in` asserts the string
/// *before* it applies both descriptions, so the failing string assertion left
/// the preservation check unrun — with the string re-blessed it executes, and
/// `apply(DUP_RUN, …)` gives `TTAACTATATATAATAATATTAT` for input and output
/// alike.
///
/// **This row no longer bounds anything.** It was #1357's negative control for
/// the over-clamp guard on the 5' junction barrier, and it passes now whether
/// that guard runs or not — the output contains no duplication at all, so the
/// property in its old name (`…_still_reaches_its_five_prime_most_position`) is
/// not exercised here. That coverage moved to
/// [`a_third_member_clear_of_the_tract_keeps_the_duplication_reaching_its_five_prime_most_position`],
/// which is the row that still goes red when the guard is removed.
#[test]
fn a_multi_base_duplication_beside_an_insertion_merges_from_the_sequence() {
    assert_normalizes_preserving_in(
        DUP_RUN,
        "TEMPLATE:g.[4_5insC;5_6dup]",
        "TEMPLATE:g.3_4insACT",
        ShuffleDirection::FivePrime,
    );
}

/// #1357's discriminator, restored on the nucleotide axis.
///
/// The barrier bounds a member by its **span start**, which is the right edge to
/// bound only for a member that claims bases. A `dup` carries its own junction at
/// its *end*, so bounding its start would hold it back by `len - 1` more than the
/// invariant requires — which is why the 5' barrier carries
/// `.filter(|_| a.junction.is_none())`.
///
/// `3_4dup` really is illegal here — it would share interbase 4 with the
/// insertion — so `4_5dup` is the 5'-most legal spelling, and the member must
/// still reach it. Without the guard it is stranded at `5_6dup`, its input
/// position: sequence-preserving, so neither exhaustive sweep can see it, and
/// silently non-canonical. The sweeps are blind to this class twice over — their
/// only oracle is `apply`, and their second member is always a single-base `dup`,
/// exactly the length at which `len - 1` is zero and the over-clamp disappears.
///
/// The third member is what keeps the property visible. `15del` stops the
/// *merge* without touching the shift: the block is no longer one run of change,
/// so the derivation declines to collapse the allele and the per-member pipeline
/// decides the output again, which is where the barrier lives. The `dup` then
/// reaches `4_5dup` and the spelling is observable.
///
/// **Distance measured, not assumed.** Sweeping every third-member position on
/// `DUP_RUN`: at `10del` and `11del` the member joins the block and the allele
/// merges again (to `g.[4_5insC;10_11insT]`) — a test written there would be
/// vacuous in a new way, pinning a merged form while claiming to pin a barrier.
/// `12del` and `13del` are clear but only 1-2 nt past the tract; at `14del` and
/// `18del` the third member itself 5'-shifts, which muddies the row. `15del` is
/// 4 nt clear, stays exactly where it is written, and reaches the tie.
///
/// **Fails both ways.** Deleting `.filter(|_| a.junction.is_none())` from the 5'
/// branch of `clamp_sibling_crossing_shifts` turns this output into
/// `g.[4_5insC;5_6dup;15del]`; restoring it turns it back. (Removing the whole
/// 5' `across_junctions` barrier instead leaves this row green — that half is
/// guarded by `an_insertion_merging_with_a_deletion_keeps_its_base_in_place`.)
#[test]
fn a_third_member_clear_of_the_tract_keeps_the_duplication_reaching_its_five_prime_most_position() {
    assert_normalizes_preserving_in(
        DUP_RUN,
        "TEMPLATE:g.[4_5insC;5_6dup;15del]",
        "TEMPLATE:g.[4_5insC;4_5dup;15del]",
        ShuffleDirection::FivePrime,
    );
}

#[test]
fn no_two_member_allele_normalizes_to_a_different_sequence() {
    // The widened class. The sweep in `repeat_span_sibling_overlap.rs` uses a
    // deletion as the first member; this one adds duplications and insertions,
    // which reach the same tracts through the junction rather than through a
    // consumed span.
    //
    // Two failure modes are counted separately: an *overlapping* output, which
    // a consumer can at least reject, and a *sequence-changing* one, which is
    // silent. Before this fix the same sweep reported 1,318 of the first and
    // 1,033 of the second; it now reports 0 and 0.
    //
    // The `dup` + `del` + 5' shape is still *tallied* on its own, pinned at
    // zero by `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES` and characterised in
    // `a_five_prime_duplication_beside_a_deletion_keeps_its_sequence`. It was
    // the last residual class here — 74 cases, closed by `blocks_sibling_shift`
    // — and keeping the separate tally says which shape regressed if one does.
    let mut checked = 0usize;
    let mut residual = 0usize;
    let mut skipped = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut changed: Vec<String> = Vec::new();

    // 24 seeds, not 48. Widening the second member from 2 shapes to 4 (#1283)
    // doubles the work per sequence, and this is already the suite's slowest
    // test — measured 45s before, 117s at 48 seeds. Halving the seed count
    // holds the runtime and the *case count* at roughly what they were while
    // covering twice the shapes, which is the trade the issue asks for: its
    // five gaps are all shape gaps, not sequence-diversity gaps, and every
    // blocking defect found so far lived in a shape the generator could not
    // emit rather than in a sequence it did not draw.
    //
    // One cost to keep on the record: `FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES`
    // below is an exact count over *this* corpus, so halving the seeds halves
    // the sequence diversity backing that pinned zero. It stays zero here, and
    // a defect would have to live only in a dropped sequence *and* in the
    // dup+del+5' shape to escape — but the guard is weaker than it was, and
    // #1295's seed knob is what restores it without paying the runtime.
    //
    // #1295 delivered that knob, so the full count goes back to 48 and the
    // diversity behind the pinned zero is restored. `sweep_seeds` returns the
    // full 48 when CI asks (`FERRO_SWEEP_SEEDS=full`) and a 4-seed prefix
    // otherwise, so a local run no longer pays for it: this test was 79.6s of an
    // 86.6s local suite at 24 seeds, and the prefix is a strict subset of the
    // corpus, not a different one.
    let seeds = sweep_seeds(48);
    for seq in sweep_sequences(seeds) {
        for first_start in 2..=13usize {
            for first_len in 1..=2usize {
                let first_end = first_start + first_len - 1;
                let span = if first_len == 1 {
                    format!("{first_start}")
                } else {
                    format!("{first_start}_{first_end}")
                };
                // #1342: the payload must not equal either reference base
                // flanking the insertion point, for the same reason the second
                // member's payload excludes both of its neighbours below.
                // This entry used to insert `seq[first_start - 1..first_end]`,
                // the span's own reference bases — so for `first_len == 1` it
                // denoted exactly `{first_start}dup`, the entry directly above
                // it. Half the `first_len` values contributed two distinct
                // first-member shapes while the array read as three, and the
                // single-base span is the case where a junction has the least
                // context to disambiguate it.
                //
                // `first_start <= 13` and `sweep_sequences` yields 20-base
                // sequences, so index `first_start` — the 3' neighbour — is
                // always in bounds.
                let first_base = seq.as_bytes()[first_start - 1] as char;
                let first_next_base = seq.as_bytes()[first_start] as char;
                let first_ins_payload = ['A', 'C', 'G', 'T']
                    .into_iter()
                    .find(|b| *b != first_base && *b != first_next_base)
                    .expect("four bases, at most two excluded");
                let inserted: String = std::iter::repeat_n(first_ins_payload, first_len).collect();
                let firsts = [
                    format!("{span}del"),
                    format!("{span}dup"),
                    format!("{first_start}_{}ins{inserted}", first_start + 1),
                ];
                for first in firsts {
                    for second_start in first_end + 2..=19usize {
                        let base = seq.as_bytes()[second_start - 1] as char;
                        let alt = if base == 'A' { 'G' } else { 'A' };
                        // #1283 gap 1: the second member used to be only
                        // `{del, sub}`, so no case in ~276,000 ever placed a
                        // `dup` or an `ins` downstream of the first member —
                        // a whole class of sibling interactions (#1279's
                        // shape) was unreachable while the file read as
                        // exhaustive.
                        //
                        // The payload is not a reference base adjacent to the
                        // junction, which the first attempt at this list got
                        // wrong twice over. Inserting the base immediately 5' of
                        // the junction denotes the same edit as
                        // `{second_start}dup` — the entry right above it — and
                        // inserting the base immediately 3' of it denotes
                        // `{second_start + 1}dup`. Either way the list reads as
                        // four shapes while carrying three, which is precisely
                        // the gap-1 defect this block exists to remove. So the
                        // payload excludes *both* neighbours rather than only
                        // the 5' one.
                        //
                        // (`firsts` builds its own payload the dup-equivalent
                        // way for `first_len == 1`. That predates this PR and is
                        // tracked as #1342 rather than widened into it.)
                        //
                        // `second_start <= 19` and `sweep_sequences` yields
                        // 20-base sequences (its argument is the *seed* count,
                        // not a length), so index `second_start` — the 3'
                        // neighbour — is at worst the last valid byte.
                        let next_base = seq.as_bytes()[second_start] as char;
                        let second_ins_payload = ['A', 'C', 'G', 'T']
                            .into_iter()
                            .find(|b| *b != base && *b != next_base)
                            .expect("four bases, at most two excluded");
                        for second in [
                            format!("{second_start}del"),
                            format!("{second_start}{base}>{alt}"),
                            format!("{second_start}dup"),
                            format!("{second_start}_{}ins{second_ins_payload}", second_start + 1),
                        ] {
                            for direction in
                                [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime]
                            {
                                let input = format!("TEMPLATE:g.[{first};{second}]");
                                let output = normalize_in(&seq, &input, direction);
                                checked += 1;
                                let Some(want) = apply(&seq, &input) else {
                                    skipped += 1;
                                    continue;
                                };
                                // The dup + del + 5' shape is *counted*, not
                                // excluded. Skipping it would drop roughly
                                // a twelfth of the sweep — some 23,000
                                // cases — out of the sequence-preservation
                                // check, and a new defect anywhere in that
                                // shape would then pass silently. Counting
                                // keeps it measured: the tally is pinned at
                                // zero, so a regression fails just as the
                                // main assertion would, while still naming
                                // the shape it landed in.
                                let residual_shape = first.ends_with("dup")
                                    && second.ends_with("del")
                                    && direction == ShuffleDirection::FivePrime;
                                match apply(&seq, &output) {
                                    None if overlapping.len() < 10 => {
                                        overlapping.push(format!("{seq}: {input} -> {output}"))
                                    }
                                    None => {}
                                    Some(got) if got != want && residual_shape => {
                                        residual += 1;
                                    }
                                    Some(got) if got != want && changed.len() < 10 => {
                                        changed.push(format!(
                                            "{seq}: {input} [{direction:?}] -> {output} (want {want}, got {got})"
                                        ));
                                    }
                                    Some(_) => {}
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Per seed, not absolute (#1295): the seed count is now a knob, so a fixed
    // floor would either fail at the default prefix or be vacuous at the full
    // corpus. Measured at 11,520 cases per seed — the loop bounds below are what
    // fix that, and they do not depend on the sequences drawn — so this floor
    // sits deliberately loose, at roughly a third, and guards against the sweep
    // being hollowed out by a lost loop rather than against exact drift.
    const CASES_PER_SEED_FLOOR: usize = 4_000;
    let floor = CASES_PER_SEED_FLOOR * seeds as usize;
    assert!(
        checked > floor,
        "sweep covered too little: {checked} cases over {seeds} seeds (floor {floor})"
    );
    // `checked` counts enumerated cases, but a case whose *input* will not
    // apply contributes nothing to either property below. Bound that share, so
    // the sweep cannot go quietly hollow while still clearing the floor.
    assert!(
        skipped * 10 < checked,
        "too many cases skipped as unconvertible: {skipped} of {checked}"
    );
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {checked} cases: {overlapping:#?}"
    );
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} cases: {changed:#?}"
    );
    // Pinned exactly at zero. Any rise is a new sequence-changing defect in the
    // dup + del + 5' shape; there is no longer a tolerated residual class here.
    assert_eq!(
        residual, FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES,
        "the dup + del + 5' shape changed"
    );
}

#[test]
fn a_five_prime_duplication_beside_a_deletion_keeps_its_sequence() {
    // #1266. Under 5' shuffle a duplication beside a deletion used to cancel
    // wrongly, emitting a shorter duplication in the wrong place:
    //
    //   reference        TTATTTAAATAAAAATAAAA
    //   g.[6_7dup;9del]  applies to  TTATTTATAATAAAAATAAAA
    //   normalized to    g.4dup      TTATTTTAAATAAAAATAAAA   <- different bases
    //
    // This was the whole 74-case residual of the sweep above. The deletion was
    // free to shift 5' through the duplicated span because a `Duplication`
    // reports claiming no reference bases — yet the bases under its span are
    // exactly what it copies. `blocks_sibling_shift` makes it a barrier.
    let seq = "TTATTTAAATAAAAATAAAA";
    assert_normalizes_preserving_in(
        seq,
        "TEMPLATE:g.[6_7dup;9del]",
        "TEMPLATE:g.7_8insT",
        ShuffleDirection::FivePrime,
    );
}

#[test]
fn a_five_prime_duplication_does_not_cross_an_upstream_sibling() {
    // #1267, the part `blocks_sibling_shift` reaches. Under 5' shuffle a
    // duplication's span travelled toward an **upstream** sibling and crossed
    // it, moving the base the sibling had edited:
    //
    //   reference        ACAAAAAAAACGTACGTACG        A-run at 3-10
    //   g.[4A>G;9dup]    applies to  ACAGAAAAAAACGTACGTACG
    //   normalized to    g.4_5insG   ACAAGAAAAAACGTACGTACG   <- moved base
    //
    // The sweep above cannot reach this: it always places the sibling
    // *downstream*, so under 5' shuffle the member travels away from it. These
    // were found by hand and stand in for the sweep until its generator is
    // widened to upstream siblings (#1283).
    //
    // The expectations are the *merged* forms since #1235: the pair is now
    // derived from the sequence it denotes rather than assembled per member, and
    // a substitution abutting an A-run duplication is one inserted `G`. The
    // property under test is unchanged — the sequence assertion below is what
    // catches a crossing — and `g.[4A>G;5dup]` was only ever the per-member
    // repair's spelling of the same bases.
    let seq = "ACAAAAAAAACGTACGTACG";
    for (input, expected) in [
        ("TEMPLATE:g.[4A>G;9dup]", "TEMPLATE:g.3_4insG"),
        ("TEMPLATE:g.[5A>G;10dup]", "TEMPLATE:g.4_5insG"),
    ] {
        assert_normalizes_preserving_in(seq, input, expected, ShuffleDirection::FivePrime);
    }
}

#[test]
fn an_insertion_junction_does_not_cross_an_upstream_sibling() {
    // The rest of #1267: the member written as an insertion. Under 5' shuffle
    // its payload travelled 5' past an upstream sibling, so the allele denoted
    // different bases while staying well-formed, disjoint and warning-free:
    //
    //   reference             ACAAAAAAAACGTACGTACG        A-run at 3-10
    //   g.[4A>G;9_10insA]     applies to  ACAGAAAAAAACGTACGTACG
    //   emitted  g.4_5insG                ACAAGAAAAAACGTACGTACG   <- moved base
    //
    // `blocks_sibling_shift` does not reach it, and neither does the span clamp:
    // the member is an `Insertion` on the way in but canonicalises to a
    // `Duplication` (`g.9_10insA` alone -> `g.3dup`), so it *grew*, and a member
    // that grew is refused there for reasons #1266/#1279 measured. What bounds
    // it is the junction clamp's 5' half, which governs **every** junction mover
    // regardless of how the input spelled it — restricting it to inputs written
    // as insertions was tried and refuted by measurement (see
    // `clamp_sibling_crossing_junctions`: over 73,376 duplication-mover cases the
    // unrestricted bound leaves zero sequence changes against the restricted
    // form's twelve).
    //
    // Member order is an input the normalizer must be indifferent to, so both
    // orders are asserted.
    let seq = "ACAAAAAAAACGTACGTACG";
    for input in ["TEMPLATE:g.[4A>G;9_10insA]", "TEMPLATE:g.[9_10insA;4A>G]"] {
        let actual = normalize_in(seq, input, ShuffleDirection::FivePrime);
        // Pinned, not merely sequence-preserving. The `apply` comparison below
        // is satisfied by an output identical to the input — which is exactly
        // what a clamp that failed to fire would produce — so on its own it
        // cannot tell "bounded correctly" from "nothing happened". The junction
        // must land at `3|4`, one short of the substituted base at 4, where it
        // coalesces with the sibling.
        assert_eq!(
            actual, "TEMPLATE:g.3_4insG",
            "{input} must bound its junction at the sibling's 5' edge"
        );
        assert_ne!(actual, input, "the clamp must actually move the junction");
        assert_eq!(
            apply(seq, &actual).expect("output applies"),
            apply(seq, input).expect("input applies"),
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

/// `(CAG)x10` at 1-based 11..=40, flanked by `T`. The tract the tandem sweep
/// uses; a junction can travel it without the payload being in phase everywhere,
/// which is what makes it the corpus for the 5' insertion-junction shape.
fn cag_tract() -> String {
    format!("{}{}{}", "T".repeat(10), "CAG".repeat(10), "T".repeat(20))
}

#[test]
fn an_insertion_junction_does_not_cross_an_upstream_junction_sibling() {
    // The other half of #1267's 5' shape. Here the upstream sibling adds
    // sequence at a junction instead of claiming bases, so `claims_bases` — the
    // bound that stops the substitution case above — does not see it at all:
    //
    //   g.[16_17insCAG;15dup]  ->  g.[11_13dup;15dup]
    //
    // The two members stay disjoint (11-13 against 15), so nothing flags it, but
    // the payload has moved from 3' of the sibling to 5' of it. Two junctions'
    // relative order is the only thing fixing the order of their payloads, which
    // is the same reasoning the 3' half already applies (#1290) — mirrored.
    let seq = cag_tract();
    for input in [
        "TEMPLATE:g.[16_17insCAG;15dup]",
        "TEMPLATE:g.[15dup;16_17insCAG]",
        "TEMPLATE:g.[37_38insCAG;36_37insTT]",
        "TEMPLATE:g.[36_37insTT;37_38insCAG]",
    ] {
        let actual = normalize_in(&seq, input, ShuffleDirection::FivePrime);
        assert_eq!(
            apply(&seq, &actual).expect("output applies"),
            apply(&seq, input).expect("input applies"),
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

#[test]
fn a_duplication_junction_does_not_cross_an_upstream_junction_sibling() {
    // The 5' junction bound applies to a **duplication** mover too, not only one
    // the input spelled as an insertion.
    //
    // Restricting it to input insertions looks like the conservative reading of
    // #1259 — a duplication's span and junction move together under a 5' shift,
    // and a blanket mirror of the 3' rule was measured there to turn 80
    // previously-correct outputs silently wrong. But that mirror bounded the
    // junction against a sibling the member was moving *away* from, which is a
    // different rule from this one. Restricting by input spelling instead leaves
    // the shape below broken: a `dup` whose junction sweeps past an upstream
    // sibling's junction, reordering the two payloads.
    //
    //   reference   GGGGGGGG + (A x 24) + GGGGGGGG
    //   g.[3dup;2_3insTT]   applies to  GGTTGGGGGGGAAA…
    //   unbounded           g.[1dup;2_3insTT]   GGGTTGGGGGG AAA…   <- payloads swapped
    //
    // Measured over 73,376 duplication-mover cases with an upstream sibling:
    // bounding every junction mover leaves **zero** sequence changes, while
    // restricting the bound to input insertions leaves **12** — all of this
    // shape. So the unrestricted rule is both simpler and strictly more correct,
    // and it does not decide the output from the input's encoding.
    let seq = format!("{}{}{}", "G".repeat(8), "A".repeat(24), "G".repeat(8));
    for input in [
        "TEMPLATE:g.[3dup;2_3insTT]",
        "TEMPLATE:g.[2_3insTT;3dup]",
        "TEMPLATE:g.[4dup;2_3insTT]",
        "TEMPLATE:g.[2_3insTT;4dup]",
    ] {
        let actual = normalize_in(&seq, input, ShuffleDirection::FivePrime);
        assert_eq!(
            apply(&seq, &actual).expect("output applies"),
            apply(&seq, input).expect("input applies"),
            "{input} -> {actual} changed the denoted sequence"
        );
    }
}

#[test]
fn a_deletion_does_not_shift_across_an_insertion_junction() {
    // The converse of the junction clamp above, and the case it did not cover.
    //
    // An insertion-like member is not a *barrier* for a span-claiming sibling —
    // it consumes no base, so `claims_reference_bases` says it claims nothing —
    // but its junction is still a position the sibling must not slide past.
    // Sliding a deletion across the point where a sibling adds sequence
    // reorders the two edits:
    //
    //   reference (A-run at 258-263)
    //   g.[258del;259_260insC]  ->  g.[259_260insC;263del]
    //   input applied   …T A C A A A A T…
    //   output applied  …T A A C A A A T…   ← the inserted base moved
    //
    // Clamped, the deletion stops at the junction, becomes adjacent to the
    // insertion, and the two coalesce.
    let mut bases = vec![b'T'; 300];
    for position in 258..=263 {
        bases[position - 1] = b'A';
    }
    let seq = String::from_utf8(bases).unwrap();

    // All three spellings of the same variant reach one form.
    for input in [
        "TEMPLATE:g.[258del;259_260insC]",
        "TEMPLATE:g.[259del;259_260insC]",
        "TEMPLATE:g.259delinsC",
    ] {
        assert_normalizes_preserving(&seq, input, "TEMPLATE:g.259A>C");
    }

    // A junction further into the tract clamps to that junction instead.
    assert_normalizes_preserving(&seq, "TEMPLATE:g.[258del;261_262insC]", "TEMPLATE:g.261A>C");

    // Negative control: with no sibling the deletion still shifts fully.
    assert_normalizes_preserving(&seq, "TEMPLATE:g.258del", "TEMPLATE:g.263del");
}

/// A non-homopolymer tandem tract, swept at every rotation and junction
/// (#1283 gap 2).
///
/// The sweep above cannot reach this class at all, and the reason is structural
/// rather than a matter of drawing more sequences:
///
/// * its corpus is 20 bp of `{A,T}` or `{A,C,G,T}` noise, which contains no
///   tandem tract long enough to canonicalise to a **non-homopolymer**
///   `Repeat`; and
/// * its insertion payloads are always a single copy, and a single-copy
///   insertion inside a tract canonicalises to a `dup` — so the multi-copy
///   repeat path is never entered.
///
/// #1280 lives exactly there: an out-of-phase payload in a `(CAG)` tract needs
/// a >=6-unit tract *and* a 2+-copy insertion. This sweeps a 10-unit `(CAG)`
/// tract, inserting 1..=3 copies of the unit at **every rotation** (`CAG`,
/// `AGC`, `GCA`) and at **every junction** across the tract, each against a
/// sibling on either side and in both shuffle directions.
///
/// The oracle is the existing SPDI apply: whatever ferro chooses to call the
/// result, it must denote the same bases. That is why this is a generator
/// change and not a new assertion — correctness is already decided.
#[test]
fn no_tandem_tract_allele_normalizes_to_a_different_sequence() {
    // 10 copies of CAG at 1-based 11..=40, flanked so members can sit outside
    // the tract on either side.
    const UNIT: &str = "CAG";
    const COPIES: usize = 10;
    let tract_start = 11usize;
    let tract: String = UNIT.repeat(COPIES);
    let tract_end = tract_start + tract.len() - 1; // 40
    let seq = format!("{}{tract}{}", "T".repeat(tract_start - 1), "T".repeat(20));

    let mut checked = 0usize;
    let mut skipped = 0usize;
    let mut overlapping: Vec<String> = Vec::new();
    let mut changed: Vec<String> = Vec::new();

    // Every rotation of the unit: an insertion whose payload is out of phase
    // with the adjacent reference is the shape #1280 is about.
    let rotations = [UNIT.to_string(), "AGC".to_string(), "GCA".to_string()];

    // `tract_start - 1 ..= tract_end`, not `tract_start..tract_end`: the two
    // boundary junctions (`10_11` and `40_41`) are where an insertion can shift
    // out of the tract entirely, which is the shape closest to #1267. Excluding
    // them left the sweep claiming "every junction across the tract" while
    // covering only the strictly-interior ones.
    for junction in tract_start - 1..=tract_end {
        for rotation in &rotations {
            for copies in 1..=3usize {
                let payload = rotation.repeat(copies);
                let first = format!("{junction}_{}ins{payload}", junction + 1);
                // A sibling on each side of the tract, and two *inside* it.
                // The inside pair is the point: an insertion in a tract
                // canonicalises to a copy count over the **whole** tract, and
                // that span can swallow a sibling sitting within it — the
                // shape `a_repeat_grown_from_an_insertion_does_not_span_a_sibling`
                // pins by hand for a homopolymer, here swept over a
                // non-homopolymer one at every rotation.
                //
                // All four sibling shapes, not just `{sub, del}`. Restricting
                // the sibling to those two is exactly the gap-1 defect this PR
                // removes from the two-member sweep above, and it would be
                // reintroduced here otherwise: a `dup` or an `ins` sibling
                // never appearing means the tract class reads as exhaustive
                // over sibling *placement* while covering half the sibling
                // *shapes*. `dup` matters most — it is the one edit whose span
                // the 5' clamp treats as claimed bases (`blocks_sibling_shift`,
                // #1286), so it is the shape most likely to interact with a
                // repeat span.
                let inside_5p = tract_start + 4;
                let inside_3p = tract_end - 4;
                let base_at = |pos: usize| seq.as_bytes()[pos - 1] as char;
                // Each sibling names its own position, so the shapes are grouped
                // by comment rather than carried alongside a redundant index.
                let siblings = [
                    format!("{}T>A", tract_start - 2), // 5' of tract
                    format!("{}del", tract_start - 2),
                    format!("{}dup", tract_start - 2),
                    format!("{}_{}insA", tract_start - 2, tract_start - 1),
                    format!("{}T>A", tract_end + 2), // 3' of tract
                    format!("{}del", tract_end + 2),
                    format!("{}dup", tract_end + 2),
                    format!("{}_{}insA", tract_end + 2, tract_end + 3),
                    format!("{inside_5p}{}>T", base_at(inside_5p)), // inside
                    format!("{inside_5p}dup"),
                    format!("{inside_3p}del"),
                    // `T` is absent from a `(CAG)` tract, so this stays an
                    // `Insertion`. Inserting `base_at(inside_3p)` here — the
                    // reference base immediately 5' of the junction — denotes the
                    // same edit as `{inside_3p}dup`, so the array would have
                    // claimed four sibling shapes while carrying three.
                    format!("{inside_3p}_{}insTT", inside_3p + 1),
                ];
                for sibling in siblings {
                    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
                        // Author the sibling first or second: member order is
                        // an input the normalizer must be indifferent to.
                        for input in [
                            format!("TEMPLATE:g.[{first};{sibling}]"),
                            format!("TEMPLATE:g.[{sibling};{first}]"),
                        ] {
                            let output = normalize_in(&seq, &input, direction);
                            checked += 1;
                            let Some(want) = apply(&seq, &input) else {
                                skipped += 1;
                                continue;
                            };
                            // Every shape is asserted, including the 5'
                            // insertion-junction one that used to be excluded and
                            // counted (see the note above the sweep). A sibling
                            // upstream of the junction is the #1267 shape; it is
                            // bounded now, so a sequence change there is a
                            // regression like any other and reports the
                            // description rather than incrementing a count.
                            match apply(&seq, &output) {
                                None if overlapping.len() < 10 => {
                                    overlapping.push(format!("{input} -> {output}"));
                                }
                                None => {}
                                Some(got) if got != want && changed.len() < 10 => {
                                    changed.push(format!(
                                        "{input} [{direction:?}] -> {output} \
                                         (want {want}, got {got})"
                                    ));
                                }
                                Some(_) => {}
                            }
                        }
                    }
                }
            }
        }
    }

    // 13,392 cases as written (31 junctions x 3 rotations x 3 copies x 12
    // siblings x 2 directions x 2 member orders). The floor sits ~20% under
    // that: low enough not to be brittle, high enough that dropping a rotation
    // or a copy count breaches it rather than passing quietly on a thinner
    // sweep.
    assert!(
        checked > 11_000,
        "tract sweep covered too little: {checked}"
    );
    // Same hollowness bound the two-member sweep carries: a case whose *input*
    // does not apply proves nothing, so cap that share rather than letting the
    // sweep pass by skipping.
    assert!(
        skipped * 10 < checked,
        "too many tract cases skipped as unconvertible: {skipped} of {checked}"
    );
    assert!(
        overlapping.is_empty(),
        "overlapping or unconvertible output in {checked} tract cases: {overlapping:#?}"
    );
    // No excluded shape: this covers the 5' insertion-junction cases too, which
    // is where the 116 pinned residuals used to live.
    assert!(
        changed.is_empty(),
        "sequence-changing normalizations in {checked} tract cases: {changed:#?}"
    );
}
