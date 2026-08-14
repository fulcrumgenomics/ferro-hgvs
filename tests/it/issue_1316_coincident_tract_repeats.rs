//! A repeat that grew its tract by more bases than the tract holds is still
//! demotable — as an insertion, not a duplication.
//!
//! [`demote_repeats_spanning_siblings`] undoes the sibling-unaware re-spelling
//! of a member as a copy count, because a repeat spans its whole tract and so
//! hides the member from every junction-based pass. It re-spells back into the
//! edit the member grew out of: a deletion, or a duplication of the tract's
//! 3'-most bases.
//!
//! Neither fits when a member *added* more bases than the reference tract holds:
//!
//! ```text
//! reference  ("ACGT" x 64) + "CAGCCAGTCAGCGCATCAG" + ("ACGT" x 64)
//!            the `A` at 262 is a one-base tract, between `C` at 261 and `G` at 263
//!
//! g.[261_262insAA;262_263insAA]
//!   each member alone grows that tract to three, and is spelled `g.262A[3]`
//!   the demotion wants `g.261_262dup` — two bases of a one-base tract
//! ```
//!
//! So the pass bailed and both members stayed repeats. A repeat carries no
//! junction (`junction_of` returns `None` for `NaEdit::Repeat`), so
//! `coalesce_members_at_one_junction` never saw them either, and the allele
//! rendered as two identical members claiming one tract — an overlap the SPDI
//! apply oracle declines (#1316).
//!
//! The one-base payload escaped only by the spelling it happened to get:
//! `g.261_262insA` normalizes to `g.262dup`, a duplication, which *does* carry a
//! junction, so the pair coalesced and re-spelled as one repeat.
//!
//! The added bases are expressible — as an insertion at the tract's 3' junction.
//! Emitting that gives both members a junction, and the existing merge finishes
//! the job: `g.262_263insAA` twice, coalesced to `insAAAA`, re-spelled `g.262A[5]`.

use crate::common::synthetic::assert_padded_preserving;

/// The core the proptest shrank to, and the reference the seed describes.
const CORE: &str = "CAGCCAGTCAGCGCATCAG";

#[test]
fn two_insertions_growing_one_tract_combine_their_copies() {
    // #1316. Both members grow the one-base `A` tract at 262 by two, so the
    // tract holds five copies — one member, not the same member twice.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[261_262insAA;262_263insAA]");
    assert_eq!(output, "NC_TEST.1:g.262A[5]");
}

#[test]
fn the_seeds_deletion_survives_the_combination() {
    // The proptest's own case, which carries a distant deletion. It is outside
    // the tract and must come through untouched.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[261_262insAA;262_263insAA;269del]");
    // **Re-pinned, and the row now does what this test's own docstring asks.**
    // The sibling `two_insertions_growing_one_tract_combine_their_copies` pins
    // `g.262A[5]` for the same two insertions without the deletion, and the
    // comment above says the distant deletion "is outside the tract and must
    // come through untouched". Both hold here now: the tract renders as one
    // repeat and `269del` is carried through verbatim.
    //
    // An intermediate revision of this branch pinned
    // `g.[263_265delinsAAA;268_269delinsTCAGC]` and explained it as the repeat
    // notation not surviving re-derivation. That was measured against a base
    // where `MAX_SPLIT_BLOCK` was 1024; it dissolved the repeat AND rewrote the
    // deletion, satisfying neither half of the intent stated above. It is
    // withdrawn rather than reworded.
    assert_eq!(output, "NC_TEST.1:g.[262A[5];269del]");
}

#[test]
fn a_one_base_payload_still_reaches_the_same_form() {
    // The shape that already worked, via the `dup` spelling rather than the
    // demotion. Pinned so the two routes stay in agreement.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[261_262insA;262_263insA]");
    assert_eq!(output, "NC_TEST.1:g.262A[3]");
}

#[test]
fn a_multi_base_unit_combines_the_same_way() {
    // The unit is read off the repeat rather than assumed to be one base: the
    // `CA` at 259_260 is a one-copy tract, and each member adds two copies of
    // it. The payload is the unit repeated, not the tract's first base.
    let output = assert_padded_preserving("TTCAGG", "NC_TEST.1:g.[258_259insCACA;260_261insCACA]");
    assert_eq!(output, "NC_TEST.1:g.259_260CA[5]");
}

#[test]
fn a_lone_repeat_growing_past_its_tract_is_left_alone() {
    // The guard: with no sibling to span, the repeat spelling is correct and
    // must survive. Demoting unconditionally would undo B2 for every member.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.261_262insAA");
    assert_eq!(output, "NC_TEST.1:g.262A[3]");
}
