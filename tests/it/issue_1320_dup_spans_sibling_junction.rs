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
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`).
    //
    //   was  NC_TEST.1:g.264_265insCAAAAA               one member, re-derived
    //   now  NC_TEST.1:g.[264_265insCA;266_267insAAAA]  two members
    //
    // The old expectation collapsed all three authored blocks into one by
    // re-deriving the partition from the resulting sequence. Under the ruling
    // the `263_264insAC` block survives on its own — the reference bases at 264
    // and 265 sit unchanged between it and the next member, a separation
    // `general.md:34` says to describe individually.
    //
    // The three-to-two collapse that *does* happen is the ruling's other
    // licensed move, MERGE: `265_266insAA` and `266_267insAA` both 3'-shift into
    // the same `A` tract and land on one junction, i.e. separation zero, which
    // is below the genomic floor. Note this is the floor applied to the
    // *canonically placed* members rather than to the written ones — the reading
    // campaign issue #65 records as still open (D2). It is flagged here rather
    // than settled: if that question is answered the other way this row moves to
    // three members, and it should redden rather than pass quietly.
    //
    // The swallowing this file guards against cannot hide in either member:
    // `assert_padded_preserving` applies input and output through `hgvs_to_spdi`
    // — not through the normalizer — and a tract covering the sibling's junction
    // denotes no sequence, which the applier declines. Measured: same denoted
    // bases, members disjoint and ascending, output a fixed point.
    let output =
        assert_padded_preserving(CORE, "NC_TEST.1:g.[263_264insAC;265_266insAA;266_267insAA]");
    assert_eq!(output, "NC_TEST.1:g.[264_265insCA;266_267insAAAA]");
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
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
    // and renamed: the old name asserted the merge, and the merge is gone.
    //
    //   was  NC_TEST.1:g.264_265insCAAA               one member, re-derived
    //   now  NC_TEST.1:g.[264_265insCA;264_265dup]    two members, as authored
    //
    // The pair's members occupy the junctions `263|264` and `264|265`, one
    // unchanged reference base (264) apart, which `general.md:34` says to
    // describe individually. So the previous single-member answer was the
    // re-derivation the ruling removes, and the two members now survive — which
    // restores exactly the shared-span tie this test is about: both span
    // `264_265`, so start and end tie and only `junction_rank` separates them.
    // This row is therefore *stronger* coverage than it was, not weaker: it went
    // from pinning a merged form that never reached the sort to pinning the sort's
    // own output.
    //
    // `assert_padded_preserving` proves the bases through `hgvs_to_spdi`,
    // independently of the normalizer, and additionally requires the members to
    // render disjoint and ascending — the property #1301 is about. Measured:
    // same denoted sequence, output a fixed point.
    let output = assert_padded_preserving("GCATGAAAAT", "NC_TEST.1:g.[263_264insAC;264_265insAA]");
    assert_eq!(output, "NC_TEST.1:g.[264_265insCA;264_265dup]");
}

#[test]
fn a_junction_at_the_dups_five_prime_end_is_left_alone_beside_a_third_member() {
    // The 5' boundary itself, restored: the deletion at 270 sits clear of the
    // A-run the payloads shift through, so the pair still re-spells exactly as
    // it did, but the derivation declines the three-member group instead of
    // merging it. `264_265dup` and the insertion share a junction at 264, which
    // `junction_rank` already orders, so the predicate must *not* re-spell the
    // dup here — doing so would undo #1301.
    let output = assert_padded_preserving(
        "GCATGAAAAT",
        "NC_TEST.1:g.[263_264insAC;264_265insAA;270del]",
    );
    assert_eq!(output, "NC_TEST.1:g.[264_265insCA;264_265dup;270del]");
}

#[test]
fn two_insertions_sharing_a_start_keep_both_dup_spellings() {
    // The same 5' boundary between two duplications, which is #1261's
    // reproduction: `258dup`'s junction is 258 and `258_259dup` starts there.
    // A `>=` test would re-spell the wider one and lose the pinned form.
    //
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
    // and renamed: the old name asserted the merge the ruling removes.
    //
    //   was  NC_TEST.1:g.258_259insAGA          one member, re-derived
    //   now  NC_TEST.1:g.[258dup;258_259dup]    two members, as authored
    //
    // The input's two insertions sit at the junctions `258|259` and `259|260`
    // with the reference base at 259 unchanged between them — separation one,
    // which `general.md:34` describes individually. The single-member answer was
    // the re-derivation from the resulting sequence that the ruling removes.
    //
    // Both `dup` spellings are back, which is what this test was written to
    // guard, so this row now exercises the shared-start tie again rather than
    // pinning a form that never reached it. The identical pair beside a third
    // member (`two_duplications_sharing_a_start_keep_both_dup_spellings_beside_a_third_member`)
    // pins the same two strings, and the two agreeing is itself evidence the
    // answer comes from the key rather than from the third member.
    //
    // `assert_padded_preserving` proves the bases and the disjoint-ascending
    // rendering through `hgvs_to_spdi`, not through the normalizer. Measured:
    // same denoted sequence, output a fixed point.
    let output =
        assert_padded_preserving("CAGTATGCAGGCAA", "NC_TEST.1:g.[258_259insA;259_260insAG]");
    assert_eq!(output, "NC_TEST.1:g.[258dup;258_259dup]");
}

#[test]
fn two_duplications_sharing_a_start_keep_both_dup_spellings_beside_a_third_member() {
    // The same 5' boundary between two duplications — #1261's reproduction —
    // restored by a third member the derivation will not merge them with.
    // `258dup`'s junction is 258 and `258_259dup` starts there; a `>=` test
    // would re-spell the wider one and lose both pinned forms.
    let output = assert_padded_preserving(
        "CAGTATGCAGGCAA",
        "NC_TEST.1:g.[258_259insA;259_260insAG;268del]",
    );
    assert_eq!(output, "NC_TEST.1:g.[258dup;258_259dup;268del]");
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
