//! A member that cancelled to `=` must not be left overlapping the member that
//! absorbed it.
//!
//! An identity member says a position was examined and found unchanged. It is
//! not authored here — it is derived: a `del` and an `ins` that merge into a
//! `delins` restating the reference cancel to `=`, and the member that absorbed
//! the change then grows over the position it names.
//!
//! ```text
//! reference  ("ACGT" x 64) + "GCATGAAAAT" + ("ACGT" x 64)
//! g.[261_262insAA;263del;264_265insA]  ->  g.[262_265A[6];265=]
//! ```
//!
//! `262_265A[6]` alone denotes the input exactly. The `265=` residue claims a
//! base the repeat covers, so the pair overlaps and the apply oracle declines
//! the output — #1235's second acceptance criterion.
//!
//! Dropping identity members outright would be wrong. `g.[1002=;1005del]` is a
//! real description: the `=` records a position that was checked. Only an
//! identity member a sibling *claims the bases of* is a contradiction, and only
//! that one is dropped.

use crate::common::synthetic::assert_padded_preserving;

const CORE: &str = "GCATGAAAAT";

#[test]
fn a_cancelled_member_does_not_survive_inside_the_repeat_that_absorbed_it() {
    // #1297. Three members reduce to one repeat; the `=` left behind by the
    // cancelling pair must go with them.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[261_262insAA;263del;264_265insA]");
    assert_eq!(output, "NC_TEST.1:g.262_265A[6]");
}

#[test]
fn a_wholly_self_cancelling_allele_still_reports_no_change() {
    // The guard on the other side: when *everything* cancels there is no
    // sibling claiming anything, so the identity is the answer and is kept.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[263del;264_265insA]");
    assert_eq!(output, "NC_TEST.1:g.266_267=");
}

#[test]
fn an_identity_member_beside_an_unrelated_deletion_is_kept() {
    // `=` records a position that was examined. A sibling elsewhere contradicts
    // nothing, so both members stay — dropping identities outright would lose
    // real information.
    let output = assert_padded_preserving("GCATGAAAATCCTTGG", "NC_TEST.1:g.[258=;266del]");
    assert_eq!(output, "NC_TEST.1:g.[258=;266del]");
}
