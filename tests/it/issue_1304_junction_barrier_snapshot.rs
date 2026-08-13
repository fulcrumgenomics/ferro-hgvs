//! A junction barrier must read its sibling's payload from the snapshot the
//! span came from.
//!
//! `clamp_sibling_crossing_junctions` compares this member's payload against
//! each sibling's to decide whether crossing is observable — commuting payloads
//! may cross, others may not (#1290). It considers each sibling from *both*
//! snapshots, before and after per-member normalization, so a sibling that moved
//! still counts as a barrier where it started.
//!
//! It read every payload off the **after** variant, including for spans taken
//! from the **before** snapshot. A duplication reads its payload from the
//! reference under its own span, so that pairing reads bases belonging to
//! neither:
//!
//! ```text
//! reference  ("ACGT" x 64) + "GCATGAAAAT" + ("ACGT" x 64)
//! g.[260_261insGA;261_262insA;264del]
//!   per-member ->  g.[261_262dup;265dup;265del]
//!   the `insA` sibling's *before* span 261_262 measured against its *after*
//!   `265dup` reports the payload `GA`, not `A`
//!
//!   intended  G C A T G A G A A A A T
//!   emitted   G C A T G A A G A A A T
//! ```
//!
//! `GA` commutes with this member's own `GA`, so the barrier vanished and the
//! `insGA` junction swept from 260 past the sibling at 261 — reordering the two
//! insertions and denoting different bases, with both members well-formed and
//! disjoint. Nothing else catches it.

use crate::common::synthetic::assert_padded_preserving;

const CORE: &str = "GCATGAAAAT";

#[test]
fn a_moved_sibling_shape_now_merges_from_the_sequence() {
    // #1304's original shape. It no longer reaches the barrier at all: since
    // #1235 — and specifically since `main` removed the input-separator veto
    // (#1345), which had been refusing this derivation for merging across the
    // base the input left at 263 — the whole allele is derived from the sequence
    // it denotes, and that sequence is one `GA` repeated at 261_262.
    //
    // This is a strictly stronger answer than the barrier's, not a weaker one:
    // the barrier stops a member from *shifting* somewhere that would denote
    // different bases, whereas the derivation never assembles members at all, so
    // there is nothing to shift. `assert_padded_preserving` proves the bases
    // (its applier is not the normalizer), and #1304's two spellings now agree —
    // see `cis_spelling_confluence_gap.rs`, where the pair moved to `CONVERGED`.
    //
    // Renamed because the old name asserted the opposite of what it now pins.
    //
    // Still a guard on the barrier, and measured to be: stubbing
    // `clamp_sibling_crossing_junctions` to return immediately turns this output
    // into `g.[261_262dup;265=]`, so the row goes red with the mechanism off and
    // green with it on. It is the only row in this file that still does — the
    // pair-alone case below merges either way now.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA;264del]");
    assert_eq!(output, "NC_TEST.1:g.261_262dup");
}

#[test]
fn the_pair_alone_merges_instead_of_staying_barred() {
    // Without the deletion, so the barrier is the only thing under test.
    //
    // Since #1235 the two-member spelling merges: with nothing else in the
    // allele the pair denotes one three-base insertion, and the merged form is
    // derived from that sequence rather than assembled per member — so the
    // barrier, which decides where a member may *shift* to, no longer decides
    // this allele's output.
    //
    // Which means this row no longer guards the barrier: with
    // `clamp_sibling_crossing_junctions` stubbed out it still passes. It is kept
    // as a value pin on the merged form, not as barrier coverage — that lives in
    // `a_moved_sibling_shape_now_merges_from_the_sequence` above, which is the
    // one row here that still goes red when the mechanism is removed.
    //
    // (An earlier revision of this comment claimed the deletion case "stays
    // three members and the barrier holds". That was true when it was written
    // and stopped being true when `main` removed the input-separator veto
    // (#1345); the deletion case now merges to a single `dup`.)
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA]");
    assert_eq!(output, "NC_TEST.1:g.262_263insGAA");
}

#[test]
fn a_third_member_clear_of_the_tract_keeps_the_pair_barred() {
    // The discriminator `the_pair_alone_...` above used to be, restored on the
    // nucleotide axis. That row stopped bounding the barrier once the pair began
    // merging from the sequence: it now passes whether the barrier runs or not,
    // so on its own it no longer says anything about crossing.
    //
    // A third member restores it by stopping the *merge* without touching the
    // shift. `270del` is 4 nt clear of the `GCATGAAAAT` tract the two payloads
    // shuffle through, so the derivation declines to collapse the allele — the
    // block is no longer one run of change — and the per-member pipeline decides
    // the output again, which is where the barrier lives. The pair then stays
    // exactly where the barrier puts it: `insGA` may not sweep from 260 past the
    // `insA` at 261, because `GA` and `A` do not commute and their order is
    // observable.
    //
    // The distance matters in both directions, and the near end is a trap worth
    // naming: at `266del` and `267del` the third member is close enough to join
    // the block, and the allele merges again (to `g.[261_262dup;266T>A]` and
    // `g.[263A>G;266T>A;267_268insAT]` respectively) — a test written there would
    // be vacuous in a new way, pinning a merged form while claiming to pin a
    // barrier. `270del` and `272del` both sit clear; 270 is used here.
    //
    // Verified to fail both ways, which is the only thing that makes it a guard:
    // stubbing `clamp_sibling_crossing_junctions` to return immediately turns
    // this output into `g.[261_262dup;265dup;270del]` — the `insGA` swept past
    // its sibling — and restoring the mechanism turns it back.
    //
    // # RE-PINNED BY THE PARTITION DEFAULT FLIP (#1835), AND THE BARRIER
    // COVERAGE IS LOST WITH IT
    //
    // Was `NC_TEST.1:g.[260_261insGA;261_262insA;270del]`; now
    // `NC_TEST.1:g.[262_263insGAA;270del]` — the same single insertion
    // `the_pair_alone_merges_instead_of_staying_barred` above already pins
    // WITHOUT the third member. The third member has stopped blocking the
    // derivation, so the device this whole test rests on no longer works.
    //
    // LICENSED BY `contiguous-insertion-split-by-a-blocked-derivation` (decided).
    // Its ruling is that `general.md:34` is stated over "two variants" and a
    // locus like this one carries ONE: aligning the reference against the denoted
    // sequence leaves a single contiguous 3 nt insertion (`GAA`). So `:34` never
    // reached this pair, and the two-member spelling survived "because a
    // derivation was blocked, not because a clause requires it". The record
    // verified that by applying both spellings through `hgvs_to_spdi`,
    // independently of the normalizer — and `assert_padded_preserving` re-checks
    // the same equality on every run here.
    //
    // THIS IS THE THIRD ROW IN THIS FILE TO STOP GUARDING THE BARRIER, and the
    // sequence is worth reading as a whole because it is now complete. The
    // comment on `the_pair_alone_merges_instead_of_staying_barred` records that
    // *it* stopped bounding the barrier when #1235 made the pair merge, and moved
    // the coverage HERE by adding a third member. That escape has now closed too:
    // with the output at one insertion, stubbing
    // `clamp_sibling_crossing_junctions` cannot change it, so the "fails both
    // ways" paragraph above is stale as a claim about the default arm. The one
    // row that still goes red when the mechanism is removed is
    // `a_moved_sibling_shape_now_merges_from_the_sequence`; if that row ever
    // merges too, this file guards nothing and the loss must be replaced rather
    // than re-pinned again.
    //
    // Recorded as a real coverage loss, not papered over. The barrier is still
    // reachable under `FERRO_PARTITION=live`, which is part of why that arm stays
    // selectable by name; re-establishing it on the default arm needs a shape
    // whose members survive the derivation, which this core does not supply.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA;270del]");
    assert_eq!(output, "NC_TEST.1:g.[262_263insGAA;270del]");
}

#[test]
fn the_barred_pair_does_not_depend_on_authored_order() {
    // The barrier reads both snapshots, so its answer must not depend on the
    // order the members were written in — the same independence the #1261/#1301
    // discriminators carry. Same allele as above, authored backwards.
    //
    // #1835: was `NC_TEST.1:g.[260_261insGA;261_262insA;270del]`; now the merged
    // `NC_TEST.1:g.[262_263insGAA;270del]`, for the reason recorded on the row
    // above.
    //
    // THE INDEPENDENCE PROPERTY SURVIVES, and is arguably strengthened. What this
    // row asserts is that the two authored orders reach the SAME output, and they
    // still do — keep this string identical to the row above, since their
    // equality is the whole test. Authored order is now irrelevant for a stronger
    // reason than before: the answer is derived from the resulting sequence, which
    // is order-invariant by construction, rather than from a barrier reading two
    // snapshots. It is no longer evidence about the barrier specifically.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[270del;261_262insA;260_261insGA]");
    assert_eq!(output, "NC_TEST.1:g.[262_263insGAA;270del]");
}

#[test]
fn each_insertion_alone_still_shifts_to_its_own_most_three_prime_form() {
    // The guard: the barrier is a *sibling* effect. On its own each member
    // still reaches its standalone 3'-most spelling, so the clamp is not
    // suppressing the 3' rule.
    for (input, expected) in [
        ("NC_TEST.1:g.260_261insGA", "NC_TEST.1:g.261_262dup"),
        ("NC_TEST.1:g.261_262insA", "NC_TEST.1:g.265dup"),
    ] {
        assert_eq!(assert_padded_preserving(CORE, input), expected);
    }
}
