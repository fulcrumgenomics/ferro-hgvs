//! A repeat that outgrew its tract must release a sibling deletion it swallowed
//! — but only one it has *actually* swallowed.
//!
//! #1325 gave [`demote_repeats_spanning_siblings`] a third target form for a
//! repeat whose growth exceeds its tract: the insertion that growth stands for,
//! which restores a junction the clamps can bound. It was scoped to the case
//! where the repeat swallows a sibling's *junction*, because an insertion
//! claims no bases and so blocks no sibling's shift — demoting before the
//! clamps run is what released the deletion #1296 had clamped out of a tract.
//!
//! That left the base-claiming branch declining, and #1394 is what it left:
//!
//! ```text
//! reference  ("ACGT" x 64) + "CAGGCAAACAGTGAAG" + ("ACGT" x 64)
//!            the `A` tract at 262_264 is three bases, between `C` at 261 and `C` at 265
//!
//! g.[262del;263_264insAA;264_265insAA]  ->  g.[262_264A[7];264del]
//!   the deletion claims 264, inside the tract 262-264
//! ```
//!
//! Same structure as #1325. The first merge pass reaches the right answer —
//! the clamp holds the deletion out of the duplications' span and the coalesce
//! combines them, giving `g.[263del;264_265insAAAA]` — and the next pass
//! re-normalizes each member alone, re-spelling the insertion as the tract-wide
//! `262_264A[7]`, whose span swallows the deletion outright.
//!
//! ## Why the gate is a *current* overlap, not a predicted one
//!
//! `spans_sibling_bases` deliberately reads both the pre- and post-pass
//! snapshots, so that a sibling which moved into the tract still counts. That
//! is right for the deletion and duplication target forms, which preserve the
//! barrier. It is wrong for the insertion form, which destroys it: in #1296 the
//! deletion meets the tract only in the *pre* snapshot, and demoting on that
//! evidence removes the very barrier that is holding it at 272_273 — the tract
//! at 274 is then free to be shifted into, which is the defect #1296 fixed.
//!
//! So the insertion route asks a stricter question: are these two members
//! overlapping *now*, in the post snapshot? #1394 is (264 lies in 262-264);
//! #1296 is not (272_273 against a tract at 274). A repeat that has already
//! swallowed a sibling has no barrier left worth preserving — the pair denotes
//! no sequence at all — whereas one that merely might is still doing its job.

use crate::common::synthetic::assert_padded_preserving;

/// The core the proptest shrank to, and the reference the seed describes.
const CORE: &str = "CAGGCAAACAGTGAAG";

#[test]
fn a_repeat_growing_past_its_tract_releases_a_swallowed_deletion() {
    // #1394, the proptest's own case. Three bases of tract, one deleted and
    // four added, so the tract holds six `A`s in the end.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[262del;263_264insAA;264_265insAA]");
    assert_eq!(output, "NC_TEST.1:g.262_264A[6]");
}

#[test]
fn the_two_insertions_alone_still_reach_the_repeat() {
    // Without the deletion there is nothing to swallow, so the copy count over
    // the whole tract is the correct spelling and must survive.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[263_264insAA;264_265insAA]");
    assert_eq!(output, "NC_TEST.1:g.262_264A[7]");
}
