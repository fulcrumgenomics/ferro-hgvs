//! A duplication's span may not contain a sibling's junction, any more than a
//! repeat's tract may.
//!
//! `demote_repeats_spanning_siblings` undoes the sibling-unaware re-spelling of
//! a member as a copy count, because a repeat spans its whole tract and can
//! swallow a sibling. For a member that *grew* the tract it re-spells to a
//! duplication over the tract's 3'-most bases — and when the growth fills the
//! tract exactly, that duplication spans the very same bases the repeat did:
//!
//! ```text
//! reference  ("ACGT" x 64) + "AACAGTAAAATAT" + ("ACGT" x 64)
//!            the `A` tract at 263_266 is four bases, between `T` at 262 and `T` at 267
//!
//! g.[263_264insAC;265_266insAA;266_267insAA]  ->  g.[263_266dup;264_265insCA]
//!   the dup spans 263-266; the insertion adds its bases at interbase 264,
//!   strictly inside it, so the pair overlaps
//! ```
//!
//! The two `insAA` members merge into the tract-wide repeat correctly — that
//! pair on its own gives `g.263_266A[8]`. The demotion then hands back a
//! duplication that swallows the third member's junction exactly as the repeat
//! had, and nothing narrows it: `clamp_sibling_crossing_shifts` bounds a member
//! against a sibling's *bases* and an insertion has none, while
//! `clamp_sibling_crossing_junctions` governs the insertion's own junction
//! rather than the duplication's span.
//!
//! `respell_colliding_duplications` is the pass that exists for exactly this —
//! a duplication whose span collides becomes the equivalent insertion, which
//! claims nothing and so cannot collide. It only tested siblings that *claim
//! bases*, so an insertion's junction inside the span did not count as a
//! collision. #1287 added that same junction half to the repeat path; this adds
//! it to the duplication path (#1320).
//!
//! Re-spelled, `g.263_266dup` becomes `g.266_267insAAAA` — the copy at the
//! junction 3' of the span it copied — and the two members are then disjoint and
//! in order.

use crate::common::synthetic::assert_padded_preserving;

/// The core the proptest shrank to, and the reference the seed describes.
const CORE: &str = "AACAGTAAAATAT";

#[test]
fn a_dup_does_not_swallow_a_siblings_junction() {
    // #1320. The seed's own case. The duplication the demotion produces spans
    // 263-266 and the sibling adds at interbase 264, inside it.
    //
    // Since #1235 the three members reach one: the merged form is derived from
    // the sequence rather than assembled per member, and this allele denotes a
    // single six-base insertion. It is the `#1320` row of
    // `cis_spelling_confluence_gap.rs`, now converged. The swallowing this file
    // guards against cannot hide inside one member — with nothing to swallow,
    // `assert_padded_preserving` still proves the sequence is the input's.
    let output =
        assert_padded_preserving(CORE, "NC_TEST.1:g.[263_264insAC;265_266insAA;266_267insAA]");
    assert_eq!(output, "NC_TEST.1:g.264_265insCAAAAA");
}

#[test]
fn the_pair_without_the_swallowed_sibling_still_uses_the_repeat() {
    // The guard on the demotion: with no sibling inside the tract there is
    // nothing to demote, and the tract-wide repeat is the right spelling.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[265_266insAA;266_267insAA]");
    assert_eq!(output, "NC_TEST.1:g.263_266A[8]");
}

#[test]
fn a_junction_at_the_dups_five_prime_end_is_left_alone() {
    // The 5' boundary, and #1301's own reproduction over its own core.
    // `264_265dup` sorts to 264, the same position the insertion adds at, so
    // `junction_rank` already orders the pair — there is no sort/act
    // discrepancy to repair, and re-spelling here would undo #1301.
    //
    // Since #1235 the pair merges before either spelling is chosen — the merged
    // form is derived from the sequence, not assembled per member — so the
    // re-spelling this test bars is not reached on this input. #1301's own file
    // records the same movement, and the 5' boundary the predicate must not
    // cross is still pinned by `the_pair_without_the_swallowed_sibling_still_uses_the_repeat`
    // and `a_deletion_colliding_with_a_dups_bases_still_respells`.
    let output = assert_padded_preserving("GCATGAAAAT", "NC_TEST.1:g.[263_264insAC;264_265insAA]");
    assert_eq!(output, "NC_TEST.1:g.264_265insCAAA");
}

#[test]
fn two_duplications_sharing_a_start_keep_both_dup_spellings() {
    // The same 5' boundary between two duplications, which is #1261's
    // reproduction: `258dup`'s junction is 258 and `258_259dup` starts there.
    // A `>=` test would re-spell the wider one and lose the pinned form.
    //
    // Since #1235 both dup spellings are gone, because there is only one member
    // left to spell: the merged form is derived from the sequence the pair
    // denotes — one three-base insertion — rather than assembled per member.
    // #1261's own file records the same movement.
    let output =
        assert_padded_preserving("CAGTATGCAGGCAA", "NC_TEST.1:g.[258_259insA;259_260insAG]");
    assert_eq!(output, "NC_TEST.1:g.258_259insAGA");
}

#[test]
fn a_deletion_colliding_with_a_dups_bases_still_respells() {
    // The pre-existing half of the collision test, which fires on a sibling
    // that claims bases. Widening the predicate must not lose it.
    // The pair collapses to one member, so there is no second member left to
    // claim the dup's bases. Pinned exactly, as the sibling test above is: the
    // disjunction this replaced ("not both `dup` and `del`") also passed for an
    // output that dropped the repair entirely.
    let output = assert_padded_preserving("GCATGAAAAT", "NC_TEST.1:g.[261_262insAA;263del]");
    assert_eq!(output, "NC_TEST.1:g.265dup");
}
