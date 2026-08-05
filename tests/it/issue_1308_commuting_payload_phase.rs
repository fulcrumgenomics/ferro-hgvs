//! Commuting payloads may meet, but not sweep past a sibling that stays put.
//!
//! `clamp_sibling_crossing_junctions` exempts a sibling whose payload
//! *commutes* with this member's: `p ++ q == q ++ p` means the two have no
//! observable order, so letting them reach one junction is better than keeping
//! them apart — they merge into a single member (#1286).
//!
//! That reasoning is about the order of the two payloads **relative to each
//! other**, once they meet. It says nothing about a member sweeping past a
//! sibling that does not move and landing beyond it, where different reference
//! bases now separate them:
//!
//! ```text
//! reference  ("ACGT" x 64) + "CAGAAGATGAATAA" + ("ACGT" x 64)
//! g.[263_264insTG;264_265insTG]
//!   the second member does not move; the first 3'-shifts 263 -> 265
//!
//!   intended  C A G A A G A T G T T G G A A T A A
//!   emitted   C A G A A G A T T G G T G A A T A A
//! ```
//!
//! Two well-formed members on distinct junctions, disjoint and in order — and
//! denoting different bases. So a commuting sibling bounds this member at the
//! junction it **settles** on, and only its post-normalization position counts,
//! since where it started is exactly the order commuting makes unobservable.
//!
//! Landing on that bound can still be out of phase, and then the member may not
//! go there either. The clamp walks back from the bound to the most 3' junction
//! this member's payload is in phase with — `before_junction` always is, since
//! that is where the shift began.
//!
//! The bound must be tested against the payload this member would carry **at**
//! that junction, not the one it carries now — landing there rotates it, and a
//! rotation can destroy the commuting identity. That half is #1312.

use crate::common::synthetic::assert_padded_preserving;

#[test]
fn a_commuting_payload_pair_merges_without_sweeping_past_its_sibling() {
    // #1308. `264_265insTG` does not move, so the first member may not shift
    // from 263 to 265 across it.
    let output =
        assert_padded_preserving("CAGAAGATGAATAA", "NC_TEST.1:g.[263_264insTG;264_265insTG]");
    // Not the input's own spelling any more. The shift clamp still refuses to
    // sweep the first member across the second — that is what this file is
    // about, and it is unchanged — but since #1235 the allele's output is
    // derived from the sequence the pair denotes rather than assembled from the
    // clamped members, and that sequence is one four-base insertion. It is the
    // `#1308` row of `cis_spelling_confluence_gap.rs`, now converged.
    assert_eq!(output, "NC_TEST.1:g.265_266insTTGG");
}

#[test]
fn commuting_payloads_that_settle_together_still_merge() {
    // The guard on the bound. Both members shift to the same junction in a
    // poly-A tract, so the exemption must still let them meet and merge — one
    // member, not two. Bounding at the sibling's *pre* junction would split it.
    let output = assert_padded_preserving("AAAAAA", "NC_TEST.1:g.[258_259insA;259_260insA]");
    assert_eq!(output, "NC_TEST.1:g.257_263A[9]");
}

#[test]
fn three_commuting_payloads_still_reach_one_member() {
    let output = assert_padded_preserving(
        "AAAAAA",
        "NC_TEST.1:g.[258_259insA;259_260insA;260_261insA]",
    );
    assert_eq!(output, "NC_TEST.1:g.257_263A[10]");
}
