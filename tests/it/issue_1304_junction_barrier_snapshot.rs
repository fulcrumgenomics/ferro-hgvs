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
fn a_moved_sibling_still_bars_the_crossing_it_started_at() {
    // #1304. The `insGA` junction may not sweep from 260 past the `insA` at
    // 261: the payloads do not commute, so their order is observable.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA;264del]");
    assert_eq!(output, "NC_TEST.1:g.[260_261insGA;261_262insA;265del]");
}

#[test]
fn the_pair_alone_is_barred_too() {
    // Without the deletion, so the barrier is the only thing under test.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260_261insGA;261_262insA]");
    assert_eq!(output, "NC_TEST.1:g.[260_261insGA;261_262insA]");
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
