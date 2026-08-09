//! A cis member that *reduces* while it shifts must still be clamped at its
//! siblings.
//!
//! `clamp_sibling_crossing_shifts` compares each member's pre- and post-
//! normalization span and pulls back one that rotated over a sibling. It
//! admitted only a **pure translation** — `a.end - b.end == delta` — and skipped
//! anything whose span length changed, on the grounds that a length change is
//! not a rotation this pass can reason about positionally.
//!
//! A member can do both at once. When an adjacent substitution and deletion
//! merge, the resulting `delins` reduces to a plain deletion *and* 5'-shifts
//! through the tract in the same pass:
//!
//! ```text
//! reference        TTTTTTTTTAATATATTTTAATAT      nine-T tract at 1-9
//! g.9_10delinsA -> g.1del                        span 9..10 becomes 1..1
//! ```
//!
//! `delta` is -8 at the start but -9 at the end, so the guard rejected it and
//! the member crossed eight positions unchecked — straight over a `1del`
//! sibling and onto it:
//!
//! ```text
//! g.[1del;9_10delinsA]  ->  g.[1del;1del]        two members claiming base 1
//! ```
//!
//! `parse_hgvs` accepts that output, so `FERRO_ASSERT_REPARSE` does not catch
//! it; the SPDI apply oracle correctly declines it (`None`) while the input
//! applies fine. #1281 reported it at three members, which is simply the
//! shortest route from clean input to such a member — the sub/del pair has to
//! merge before it can reduce.
//!
//! The fix admits a member that **shrank** while its start moved, alongside the
//! pure translation. Both sweep a window this pass can bound, and the existing
//! cap ("never move past where the member started") still applies. A member
//! that *grew* is still refused: an insertion canonicalising to a duplication
//! moves its span and its junction together, and bounding that mis-places the
//! copy (#1266, #1279).
//!
//! Clamped at the sibling, the member lands merely adjacent to it, and the next
//! pass's merge coalesces the two into the correct single deletion —
//! `g.[1del;9_10delinsA]` -> `g.1_2del` under the 5' shuffle these cases use.
//! Fed that adjacent pair directly and shuffled 3' instead, the merge reaches
//! the same two-base deletion, which on this nine-`T` tract the repeat path
//! spells `g.[1del;2del]` -> `g.1_9T[7]`.

use crate::common::cis_apply_oracle::{
    assert_normalizes_preserving, assert_normalizes_preserving_in,
};
use ferro_hgvs::ShuffleDirection;

/// Nine `T` at 1-9 — long enough for a reduced member to shift a clear eight
/// positions before reaching the sibling at position 1.
const TRACT: &str = "TTTTTTTTTAATATATTTTAATAT";

/// Assert 5'-shuffled `input` normalizes to `expected` and denotes the same
/// bases, through the oracle's shared direction-aware assertion.
fn assert_five_prime_preserving(seq: &str, input: &str, expected: &str) {
    assert_normalizes_preserving_in(seq, input, expected, ShuffleDirection::FivePrime);
}

#[test]
fn a_reducing_member_does_not_shift_onto_a_sibling() {
    // The minimal shape: two members, the second reduces 9_10 -> 1 as it shifts.
    assert_five_prime_preserving(TRACT, "TEMPLATE:g.[1del;9_10delinsA]", "TEMPLATE:g.1_2del");
}

#[test]
fn the_three_member_allele_from_the_issue_keeps_its_substitution() {
    // #1281 as filed. The `9T>A`/`11del` pair merges, then reduces and shifts.
    assert_five_prime_preserving(TRACT, "TEMPLATE:g.[1del;9T>A;11del]", "TEMPLATE:g.1_2del");
}

#[test]
fn the_second_three_member_allele_from_the_issue_keeps_its_substitution() {
    // Pin the settled spelling, not merely the absence of the broken one:
    // `g.[4del;4del]` before the fix, the coalesced two-base deletion after.
    assert_five_prime_preserving(
        "TTATTTAAATAAAAATAAAATTTT",
        "TEMPLATE:g.[4del;6T>A;8del]",
        "TEMPLATE:g.4_5del",
    );
}

#[test]
fn the_clamped_member_still_coalesces_with_its_sibling() {
    // The clamp leaves the two merely adjacent; the merge must then close them
    // into one deletion rather than leaving `[1del;2del]` standing.
    assert_normalizes_preserving(TRACT, "TEMPLATE:g.[1del;2del]", "TEMPLATE:g.1_9T[7]");
}

#[test]
fn a_reducing_member_with_no_sibling_in_its_path_is_untouched() {
    // The relaxation must not pull back a member that crossed nothing: with the
    // sibling far 3' of the tract, the reduced member keeps its 5'-most spot.
    assert_five_prime_preserving(
        TRACT,
        "TEMPLATE:g.[9_10delinsA;20A>C]",
        "TEMPLATE:g.[1del;20A>C]",
    );
}

#[test]
fn the_three_prime_direction_keeps_the_substitution_and_the_far_deletion_apart() {
    // The 3' half of #1281's three-member allele.
    //
    // RE-BLESSED under `partition-is-the-unit-of-normalization` (DECIDED,
    // 2026-08-08, `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`),
    // and renamed: "is unchanged" stopped being true of the string, and the
    // property the row actually owns is that the two far-apart blocks are not
    // fused.
    //
    //   was  TEMPLATE:g.1_9T[7]              one member, re-derived
    //   now  TEMPLATE:g.[8_9delinsA;11del]   two members
    //
    // The input is `g.[1del;9T>A;11del]` over a nine-`T` tract. `1del` and
    // `9T>A` are adjacent once the deletion 3'-shifts through the tract, so they
    // merge — separation zero, which is the MERGE move the ruling licenses — and
    // become `8_9delinsA`. `11del` is one unchanged base away from that block —
    // the block ends at 9 and `11del` takes 11, so base 10 alone sits between
    // them — which `general.md:34` says to describe individually, so it stays its
    // own member. Separation one already clears that floor, so the ruling is
    // unchanged by the correction; only the count was wrong. The old `g.1_9T[7]` fused all three into one repeat spelling,
    // which required re-deriving the partition from the resulting sequence.
    //
    // The sequence is unchanged: `assert_normalizes_preserving` applies both
    // descriptions through the shared `cis_apply_oracle` walk — `hgvs_to_spdi`
    // plus a `claimed_from` cursor, not the normalizer — and declines an
    // overlapping output. Measured: both denote `TTTTTTTAATATATTTTAATAT`, and
    // the output is a fixed point.
    //
    // The clamp this file is about is untouched. Its own rows are the 5' ones
    // above, which are unmoved, and `the_clamped_member_still_coalesces_with_its_sibling`
    // still pins `g.[1del;2del]` -> `g.1_9T[7]` — i.e. the repeat spelling is
    // still reached where the input actually asserts one adjacent pair.
    assert_normalizes_preserving(
        TRACT,
        "TEMPLATE:g.[1del;9T>A;11del]",
        "TEMPLATE:g.[8_9delinsA;11del]",
    );
}
