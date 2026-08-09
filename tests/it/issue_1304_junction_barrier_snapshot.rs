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

/// # RED on this stack, and DELIBERATELY not re-blessed
///
/// Ferro emits `g.[261_262dup;265=]` here rather than the pinned
/// `g.261_262dup` — an identity member appearing beside the real one, which is
/// its own known class and not an answer worth pinning.
///
/// **Measured attribution.** The row is green under `FERRO_PARTITION=live`, so
/// it belongs to the partition-preserving-arm defect cluster documented in
/// `tests/it/cis_junction_crossing_shift.rs`'s module doc, and it fails
/// identically on the base branch (#1571) — nothing on this corpus branch moved
/// it. Note also that `g.[261_262dup;265=]` is precisely what the note below
/// records for the barrier **stubbed out**, so re-blessing to it would pin the
/// mechanism-off answer as the expectation and destroy this row's standing as a
/// guard. Left as it is until the arm is fixed.
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
fn the_pair_alone_keeps_both_members_and_each_reaches_its_own_three_prime_form() {
    // Without the deletion, so the barrier is the only thing under test.
    //
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
    // and renamed: the old name asserted the merge the ruling removes.
    //
    //   was  NC_TEST.1:g.262_263insGAA          one member, re-derived
    //   now  NC_TEST.1:g.[261_262dup;265dup]    two members, as authored
    //
    // The pair adds sequence at the junctions `260|261` and `261|262`, one
    // unchanged reference base (261, the `G`) apart, so `general.md:34` says to
    // describe them individually. The single-member answer was the
    // re-derivation of the partition from the resulting sequence, which is the
    // move the ruling removes; each member now 3'-shifts inside its own
    // territory and types as a `dup` — the same two spellings
    // `each_insertion_alone_still_shifts_to_its_own_most_three_prime_form`
    // pins for the members in isolation.
    //
    // **This row DOES guard the barrier again, contrary to the note it
    // replaces.** Measured on this branch by stubbing
    // `clamp_sibling_crossing_junctions` to return immediately: the output
    // becomes `g.[260_261insGA;261_262insA]` — the members do not move at all —
    // so the row is red with the mechanism off and green with it on. The old
    // note ("with the mechanism stubbed out it still passes") was measured
    // before the pair stopped merging and is no longer true.
    //
    // Sequence unchanged: `assert_padded_preserving` applies input and output
    // through `hgvs_to_spdi`, independently of the normalizer, and also requires
    // the members to render disjoint and ascending. Measured: same denoted
    // bases, output a fixed point.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA]");
    assert_eq!(output, "NC_TEST.1:g.[261_262dup;265dup]");
}

#[test]
fn a_third_member_clear_of_the_tract_leaves_the_pair_shifting_independently() {
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
    // and renamed: the old name asserted that the pair stays put, and it does
    // not.
    //
    //   was  NC_TEST.1:g.[260_261insGA;261_262insA;270del]   input, unmoved
    //   now  NC_TEST.1:g.[261_262dup;265dup;270del]          three members, shifted
    //
    // `270del` is 4 nt clear of the `GCATGAAAAT` tract the two payloads shuffle
    // through, so it does not move them; what it does is keep the allele at
    // three members. All three survive — the members are two or more unchanged
    // nucleotides apart, which `general.md:34` describes individually — so the
    // partition is exactly what the input asserted, and each member then reaches
    // its own most-3' spelling inside its own territory. That is why the pair's
    // two spellings here are byte-identical to the ones
    // `the_pair_alone_keeps_both_members_and_each_reaches_its_own_three_prime_form`
    // and `each_insertion_alone_still_shifts_to_its_own_most_three_prime_form`
    // pin, which is the cross-check that the third member is not deciding them.
    //
    // The distance matters in both directions, and the near end is a trap worth
    // naming: at `266del` and `267del` the third member is close enough to join
    // the block. `270del` and `272del` both sit clear; 270 is used here.
    //
    // **Still a guard on the barrier, and the direction of the old measurement
    // is now inverted — re-measured on this branch rather than carried over.**
    // Stubbing `clamp_sibling_crossing_junctions` to return immediately turns
    // this output into `g.[260_261insGA;261_262insA;270del]`, the *old* pinned
    // string: with the mechanism off the members do not move at all. So the row
    // is red with the mechanism removed and green with it on, which is what
    // makes it a guard. The comment this replaces recorded the opposite
    // assignment of the two strings; it was measured before the ruling landed.
    //
    // The relative order of the two payloads is preserved by the shift — `GA`
    // still lands 5' of `A` — and the sequence is unchanged, proved by
    // `assert_padded_preserving`'s `hgvs_to_spdi` applier rather than by the
    // normalizer. Measured: same denoted bases, members disjoint and ascending,
    // output a fixed point.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA;270del]");
    assert_eq!(output, "NC_TEST.1:g.[261_262dup;265dup;270del]");
}

#[test]
fn the_shifted_trio_does_not_depend_on_authored_order() {
    // The barrier reads both snapshots, so its answer must not depend on the
    // order the members were written in — the same independence the #1261/#1301
    // discriminators carry. Same allele as above, authored backwards.
    //
    // RE-BLESSED alongside its forward twin under
    // `partition-is-the-unit-of-normalization` (DECIDED, 2026-08-08,
    // `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`), and
    // renamed with it.
    //
    //   was  NC_TEST.1:g.[260_261insGA;261_262insA;270del]
    //   now  NC_TEST.1:g.[261_262dup;265dup;270del]
    //
    // The *property* this row owns — that the answer is a function of the
    // members and not of the order they were typed in — is untouched by the
    // ruling and still holds: the expectation here is byte-identical to the
    // forward authoring's, which is the whole assertion. Only the string the two
    // authorings agree on moved, for the reason given on the forward test.
    //
    // Sequence unchanged, via `assert_padded_preserving`'s independent applier.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[270del;261_262insA;260_261insGA]");
    assert_eq!(output, "NC_TEST.1:g.[261_262dup;265dup;270del]");
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
