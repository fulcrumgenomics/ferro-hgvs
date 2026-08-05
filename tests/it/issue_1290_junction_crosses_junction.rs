//! A junction must not shift past another junction.
//!
//! `clamp_sibling_crossing_junctions` bounds a member's junction against a
//! sibling's **bases**. Two junction-occupying members have none, so neither was
//! a barrier to the other and one could sweep straight past:
//!
//! ```text
//! reference  ("ACGT" x 64) + "ATACAGAAAATCAGGGCATA" + ("ACGT" x 64)
//! g.[263_264insA;265_266insC]  ->  g.[265_266insC;266dup]
//!   intended  G A A A A C A T
//!   emitted   G A A A C A A T      <- the inserted A crossed the C
//! ```
//!
//! Both members survive and land on *different* junctions, so there is no
//! overlap to detect: the output is well-formed, disjoint and warning-free, and
//! simply denotes different bases. The quietest member of the family that
//! #1276/#1277 (a dup occupying a junction), #1286 (two members on one junction)
//! and #1287 (a junction inside a repeat's tract) also belong to.
//!
//! The barrier stops one position **short** of the sibling's junction. Landing
//! *on* it makes the two share a junction, which is the case with no defined
//! payload order at all.
//!
//! It applies only when the two payloads do **not commute** (`p ++ q != q ++ p`).
//! Commuting payloads have no observable order — `A` and `AA` in a poly-A tract —
//! so crossing is harmless there, and letting them meet is better: they merge
//! into one member (#1286), which is the canonical form the sequence-first
//! derivation also reaches. `coalesce_members_at_one_junction` uses the same
//! predicate, so a pair is either kept apart and ordered, or allowed together
//! and merged — never left overlapping.

use crate::common::synthetic::assert_padded_preserving;

#[test]
fn a_junction_does_not_cross_a_non_commuting_junction() {
    // #1290. The `A` stops one short of the `C`'s junction, keeping the two in
    // their authored order.
    //
    // Since #1235 the pair does not stay apart at all: the merged form is derived
    // from the sequence rather than assembled per member, and these two denote
    // one two-base insertion — the `#1290` row of
    // `cis_spelling_confluence_gap.rs`, now converged. The crossing this test
    // guards against would change the sequence, which
    // `assert_padded_preserving` still checks against an independent applier.
    let output = assert_padded_preserving(
        "ATACAGAAAATCAGGGCATA",
        "NC_TEST.1:g.[263_264insA;265_266insC]",
    );
    assert_eq!(output, "NC_TEST.1:g.266_267insCA");
}

#[test]
fn the_barrier_does_not_depend_on_the_authored_order() {
    let forward = assert_padded_preserving(
        "ATACAGAAAATCAGGGCATA",
        "NC_TEST.1:g.[263_264insA;265_266insC]",
    );
    let reverse = assert_padded_preserving(
        "ATACAGAAAATCAGGGCATA",
        "NC_TEST.1:g.[265_266insC;263_264insA]",
    );
    assert_eq!(forward, reverse);
}

#[test]
fn commuting_payloads_still_meet_and_merge() {
    // The guard on the barrier. Two `A` payloads in a poly-A tract commute, so
    // they are allowed to reach one junction and merge — one member, not two.
    let output = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    assert_eq!(output, "NC_TEST.1:g.257_263A[9]");
}

#[test]
fn a_merged_payload_still_commutes_with_a_remaining_one() {
    // Three insertions merge to one member only because `A` commutes with the
    // merged `AA` — an equality test would stop after the first pair and leave
    // two members.
    let output = assert_padded_preserving(
        "AAAAAA",
        "NC_TEST.1:g.[258_259insA;259_260insA;260_261insA]",
    );
    assert_eq!(output, "NC_TEST.1:g.257_263A[10]");
}
