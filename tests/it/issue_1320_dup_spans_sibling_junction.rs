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
fn a_junction_at_the_dups_five_prime_end_merges_instead_of_respelling() {
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
fn a_junction_at_the_dups_five_prime_end_is_left_alone_beside_a_third_member() {
    // The 5' boundary itself, restored: the deletion at 270 sits clear of the
    // A-run the payloads shift through, so the pair still re-spells exactly as
    // it did, but the derivation declines the three-member group instead of
    // merging it. `264_265dup` and the insertion share a junction at 264, which
    // `junction_rank` already orders, so the predicate must *not* re-spell the
    // dup here — doing so would undo #1301.
    //
    // # RE-PINNED BY THE PARTITION DEFAULT FLIP (#1835)
    //
    // Was `NC_TEST.1:g.[264_265insCA;264_265dup;270del]`; now
    // `NC_TEST.1:g.[264_265insCAAA;270del]` — the same single insertion the
    // sibling test above already pins without the third member, so the third
    // member no longer blocks the derivation.
    //
    // LICENSED BY `duplication-must-ranks-the-label-not-the-partition` (decided,
    // operator ruling 2026-08-13), which names this test by its full path. It is
    // one of the two rows that record classes as inside
    // `contiguous-insertion-split-by-a-blocked-derivation`'s STATED reach rather
    // than merely compatible with it: that record's REPRESENTATION EFFECT
    // paragraph predicted this geometry by description in advance — "any stored
    // allele whose members are a junction-adjacent `ins`+`dup` pair sitting
    // beside a third member far enough away to block the derivation" moves "to a
    // single insertion" — and names the pinned
    // `[264_265insCA;264_265dup;270del]` as that geometry with the third member
    // 5 nt clear. So this is a predicted gap closing, not a discovered move.
    //
    // The `dup` is not lost to a MUST being ignored. `duplication.md:18` is
    // applied per piece of the partition re-derived from the resulting sequence,
    // and the single `CAAA` insertion is not a copy of the reference bases
    // immediately 5' of its insertion point, so `:18` never fires on it.
    //
    // WHAT IS NO LONGER GUARDED HERE. The 5'-boundary predicate this test was
    // written for — that `junction_rank` already orders a dup sharing its
    // insertion's junction, so re-spelling would undo #1301 — is not exercised by
    // a one-member output. It stays covered by
    // `the_pair_without_the_swallowed_sibling_still_uses_the_repeat` and
    // `a_deletion_colliding_with_a_dups_bases_still_respells`, which the sibling
    // test above already names for the same reason.
    let output = assert_padded_preserving(
        "GCATGAAAAT",
        "NC_TEST.1:g.[263_264insAC;264_265insAA;270del]",
    );
    assert_eq!(output, "NC_TEST.1:g.[264_265insCAAA;270del]");
}

#[test]
fn two_insertions_sharing_a_start_merge_into_one_member() {
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
fn two_duplications_sharing_a_start_keep_both_dup_spellings_beside_a_third_member() {
    // The same 5' boundary between two duplications — #1261's reproduction —
    // restored by a third member the derivation will not merge them with.
    // `258dup`'s junction is 258 and `258_259dup` starts there; a `>=` test
    // would re-spell the wider one and lose both pinned forms.
    //
    // # RE-PINNED BY THE PARTITION DEFAULT FLIP (#1835)
    //
    // Was `NC_TEST.1:g.[258dup;258_259dup;268del]`; now
    // `NC_TEST.1:g.[258_259insAGA;268del]`. The third member no longer blocks the
    // derivation, so the pair collapses to the single insertion the sibling test
    // `two_insertions_sharing_a_start_merge_into_one_member` already pins without
    // it.
    //
    // LICENSED BY `duplication-must-ranks-the-label-not-the-partition` (decided,
    // 2026-08-13). This is the `dup`+`dup`-sharing-a-start shape that record
    // names as the one "no decided record had reasoned about at all" before it —
    // so unlike its `issue_1301`/`issue_1320` siblings this row is *supplied* by
    // that ruling rather than predicted by an earlier one. Its twin is
    // `issue_1261_cis_member_order::two_duplications_sharing_a_start_render_in_junction_order_beside_a_third_member`,
    // which moves identically over the same core.
    //
    // `duplication.md:18` is not deviated from: applied per piece of the
    // re-derived partition, the single `AGA` insertion is not a copy of the
    // reference bases immediately 5' of its insertion point, so the MUST does not
    // fire on it.
    //
    // The `>=` regression this test names is no longer reachable through this
    // input, for the same reason as its sibling above; the note there records
    // where that coverage still lives.
    let output = assert_padded_preserving(
        "CAGTATGCAGGCAA",
        "NC_TEST.1:g.[258_259insA;259_260insAG;268del]",
    );
    assert_eq!(output, "NC_TEST.1:g.[258_259insAGA;268del]");
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
