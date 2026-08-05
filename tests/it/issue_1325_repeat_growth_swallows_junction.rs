//! A repeat whose growth exceeds its tract is demotable when it swallows a
//! sibling's junction — as an insertion, which the clamp can then bound.
//!
//! [`demote_repeats_spanning_siblings`] re-spells a tract-wide copy count back
//! into the edit it grew out of, so that the junction-based passes can see the
//! member again. It offers two target forms, a deletion or a duplication of the
//! tract's 3'-most bases, and neither can express growth larger than the tract:
//!
//! ```text
//! reference  ("ACGT" x 64) + "GATCATAAATTCAGC" + ("ACGT" x 64)
//!            the `A` tract at 263_265 is three bases, between `T` at 262 and `T` at 266
//!
//! g.[262_263insAA;263_264insAA;264_265insC]
//! ```
//!
//! The first two members grow that three-base tract by four, so the equivalent
//! duplication would start at `265 - 4 + 1 = 262` — 5' of the tract's own start.
//! The pass declined, and the repeat kept a span that swallowed the third
//! member's junction at 264, which is an overlap (#1325).
//!
//! What makes the decline load-bearing rather than merely incomplete is *when*
//! it happens. The first merge pass reaches the right answer on its own: the
//! clamp bounds both duplications short of the sibling's junction and the
//! coalesce combines them, giving `g.[263_264insAAAA;264_265insC]`. The second
//! pass then re-normalizes each member in isolation, which re-spells
//! `263_264insAAAA` as the tract-wide `263_265A[7]` — and a repeat carries no
//! junction (`junction_of` returns `None` for `NaEdit::Repeat`), so the clamp
//! no longer sees it. The demotion is the only thing that can hand the member
//! back to the clamp, so declining loses a bound that had already been
//! computed correctly.
//!
//! The growth *is* expressible, as an insertion at the tract's 3' junction —
//! `repeat_growth_as_insertion`, the payload #1316 already builds. Emitting it
//! restores the junction, and `clamp_sibling_crossing_junctions` runs
//! immediately afterwards and pulls it back to the last position before the
//! sibling, which is where the first pass had it.
//!
//! Scoped to the junction case on purpose. An insertion claims no bases and so
//! blocks no sibling's shift, which is why demoting a repeat to one *before*
//! the clamps released a sibling deletion that #1296 had deliberately clamped
//! out of the tract. That shape trips `spans_sibling_bases`; this one trips
//! `spans_sibling_junction`, and only the latter takes the new route.

use crate::common::synthetic::assert_padded_preserving;

/// The core the proptest shrank to, and the reference the seed describes.
const CORE: &str = "GATCATAAATTCAGC";

#[test]
fn a_repeat_growing_past_its_tract_releases_a_swallowed_junction() {
    // #1325, the proptest's own case. The four added `A`s sit 5' of the `C`,
    // so the insertion carrying them stays at `263_264` rather than reaching
    // the tract's 3' end.
    let output =
        assert_padded_preserving(CORE, "NC_TEST.1:g.[262_263insAA;263_264insAA;264_265insC]");
    assert_eq!(output, "NC_TEST.1:g.[263_264insAAAA;264_265insC]");
}

#[test]
fn the_same_growth_with_no_sibling_stays_a_repeat() {
    // The guard. With nothing to swallow, the copy count is the correct
    // spelling and B2 must not be undone — the same invariant #1316 pinned.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[262_263insAA;263_264insAA]");
    assert_eq!(output, "NC_TEST.1:g.263_265A[7]");
}

#[test]
fn a_sibling_flush_against_the_tracts_end_still_coalesces() {
    // The 3' end is exclusive: a junction at the tract's own end is flush
    // against it rather than inside it, so the repeat keeps its span and the
    // pair stays disjoint. Already correct before #1325 — pinned so the new
    // route does not widen into it and demote a member that needs no repair.
    let output =
        assert_padded_preserving(CORE, "NC_TEST.1:g.[262_263insAA;263_264insAA;265_266insC]");
    assert_eq!(output, "NC_TEST.1:g.[263_265A[7];265_266insC]");
}
