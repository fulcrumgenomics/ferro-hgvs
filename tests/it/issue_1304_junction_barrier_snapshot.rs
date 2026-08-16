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

use crate::common::synthetic::{assert_padded_preserving, extended_core};

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
    // shift, so the per-member pipeline decides the output again — which is
    // where the barrier lives. The pair then stays exactly where the barrier
    // puts it: `insGA` may not sweep from 260 past the `insA` at 261, because
    // `GA` and `A` do not commute and their order is observable.
    //
    // **Distance is no longer what stops the merge.** This row used `270del`, 4
    // nt clear of the `GCATGAAAAT` tract, on the reasoning that the block was
    // then more than one run of change. The partition default flip (**#1835**)
    // removed that: the derivation now answers the whole allele however far off
    // the third member sits, and the row went vacuous — it read
    // `g.[262_263insGAA;270del]`, the pair merged. The member is therefore
    // placed past `MAX_CANONICAL_WINDOW` (4096) instead, which refuses the
    // window before any alignment runs; see
    // `synthetic::extended_core` for why extending the core cannot disturb the
    // pair. The near-end trap it replaces is worth keeping: at `266del` and
    // `267del` the third member joined the block outright.
    //
    // Verified to fail both ways, which is the only thing that makes it a guard:
    // stubbing `clamp_sibling_crossing_junctions` to return immediately turns
    // this output into `g.[261_262dup;265dup;4500del]` — the `insGA` swept past
    // its sibling — and restoring the mechanism turns it back.
    let output = assert_padded_preserving(
        &extended_core(CORE, 4400),
        "NC_TEST.1:g.[260_261insGA;261_262insA;4500del]",
    );
    assert_eq!(output, "NC_TEST.1:g.[260_261insGA;261_262insA;4500del]");
}

#[test]
fn the_barred_pair_does_not_depend_on_authored_order() {
    // The barrier reads both snapshots, so its answer must not depend on the
    // order the members were written in — the same independence the #1261/#1301
    // discriminators carry. Same allele as above, authored backwards.
    let output = assert_padded_preserving(
        &extended_core(CORE, 4400),
        "NC_TEST.1:g.[4500del;261_262insA;260_261insGA]",
    );
    assert_eq!(output, "NC_TEST.1:g.[260_261insGA;261_262insA;4500del]");
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
