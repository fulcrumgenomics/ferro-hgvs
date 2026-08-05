//! Two members that arrive at one junction are ordered by where they came
//! from, not by how they happen to be spelled.
//!
//! [`coalesce_members_at_one_junction`] merges members that settle on one
//! interbase (#1286). It required their payloads to **commute**, as an
//! order-independence guard: `p ++ q == q ++ p` means the concatenation cannot
//! depend on the order the members happen to be in. When they do not commute it
//! declined, leaving the pair sharing a junction with no defined order:
//!
//! ```text
//! reference  ("ACGT" x 64) + "CAGGGATCAT" + ("ACGT" x 64)
//!            G at 259-261, A at 262, T at 263
//!
//! g.[260del;261_262insGA;262_263insGA]  ->  g.[261_262dup;262dup]
//!   intended  C A G G G A A G A T C A T
//!   emitted   C A G G G A G A A T C A T
//! ```
//!
//! Both members occupy junction 262 — `262dup` carrying `A`, `261_262dup`
//! carrying `GA` — and the renderer concatenates them in span order, giving
//! `GA ++ A` where the input means `A ++ GA`.
//!
//! Neither member *moved* to reach that junction, so
//! `clamp_sibling_crossing_junctions` never saw a crossing to bound: `261del`
//! and `261_262insGA` collapsed into `261delinsGA` (which canonicalizes to
//! `262dup`, absorbing the deletion — it is not dropped), while
//! `262_263insGA` re-spelled to `261_262dup` in place.
//!
//! Separating them is not available. Moving `261_262dup` back to junction 261
//! is a legal rotation — its last base `A` equals `ref[262]` — but the resulting
//! `[261_262insAG;262dup]` denotes `GAGAAT`, not `GAAGAT`. Both members belong
//! at 262; only their order is wrong.
//!
//! The order is recoverable from the `before` snapshot the sibling passes
//! already take: `261delinsGA` sat at 261 and `262_263insGA` at 262, so the
//! 5'-most original contributes its payload first (#1323).

use crate::common::synthetic::assert_padded_preserving;

/// The core the proptest shrank to, and the reference the seed describes.
const CORE: &str = "CAGGGATCAT";

#[test]
fn two_members_reaching_one_junction_keep_their_original_order() {
    // #1323. The seed's own case. Both members land on junction 262 and the
    // pair must concatenate 5'-most-original first: `A ++ GA`, not `GA ++ A`.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[260del;261_262insGA;262_263insGA]");
    assert_eq!(output, "NC_TEST.1:g.262_263insAGA");
}

#[test]
fn a_commuting_pair_still_merges_as_before() {
    // The guard on the relaxed order-independence test. In a poly-A tract every
    // payload is a power of `A`, so both orders agree and the merge is the one
    // #1286 installed — ordering by origin must not disturb it.
    let output = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    assert_eq!(output, "NC_TEST.1:g.257_263A[9]");
}

#[test]
fn a_non_commuting_pair_kept_apart_by_the_clamp_still_merges() {
    // #1312's case, which the clamp resolves by keeping the two members on
    // distinct junctions. That must keep happening: merging is the repair for a
    // pair that *cannot* be separated, not a licence to fold together members
    // the clamp already ordered.
    //
    // The clamp still keeps the two members on distinct junctions — that is
    // what this test is about and it is unchanged. What moved is the allele's
    // *output*: since #1235 it is derived from the sequence the pair denotes
    // rather than assembled from the clamped members, and that sequence is one
    // four-base insertion (the `#1312` row of `cis_spelling_confluence_gap.rs`,
    // now converged). Merging by derivation is not the "licence to fold
    // together" this comment warns about: the clamp's job is to decide where a
    // member may shift, and it still refuses.
    let output = assert_padded_preserving("TAAAACCA", "NC_TEST.1:g.[260_261insAC;261_262insAC]");
    assert_eq!(output, "NC_TEST.1:g.262_263insAACC");
}
