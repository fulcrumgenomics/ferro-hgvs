//! A repeat's tract must not swallow a sibling's **junction**.
//!
//! `demote_repeats_spanning_siblings` re-spells a repeat back to the plain edit
//! it grew from when the tract-wide span it canonicalised to would cover a
//! sibling. It asked only whether a sibling *claims bases*:
//!
//! ```text
//! reference  ("ACGT" x 64) + "ATACAGAAAATCAGGGCATA" + ("ACGT" x 64)
//! g.[261_262insGA;263_264insAA]  ->  g.[262_263dup;263_266A[6]]
//! ```
//!
//! The duplication claims no bases — it occupies the junction at 263 — so the
//! demotion did not fire, and the repeat kept a span (262-266) containing that
//! junction. The two members overlap, and the pair has no well-defined resulting
//! sequence: `parse_hgvs` accepts the output, the SPDI apply oracle declines it.
//!
//! A copy count over `262_266` cannot express another member adding sequence in
//! the middle of those bases, so a junction strictly inside the tract is
//! swallowed exactly as claimed bases are. The 3' end stays exclusive: a
//! junction at the tract's end is *flush against* it rather than inside it,
//! which is the adjacency the collapse pass exists to catch (#999, #1135) — the
//! same strictness `detect_insertion_overlaps` uses for an interior junction
//! (#1276).

use crate::common::cis_apply_oracle::assert_normalizes_preserving;
use crate::common::synthetic::assert_padded_preserving;

#[test]
fn a_repeat_does_not_swallow_a_siblings_junction() {
    // #1287. Demoted, the repeat becomes a duplication that clears the sibling.
    //
    // Since #1235 the two members never reach the demotion: the merged form is
    // derived from the sequence rather than assembled per member, and this pair
    // denotes one four-base insertion. That is the `#1287` row of
    // `cis_spelling_confluence_gap.rs`, which this convergence moved out of the
    // divergent table. The swallowing this test guards against is still caught —
    // `assert_padded_preserving` applies both spellings independently of the
    // normalizer.
    let output = assert_padded_preserving(
        "ATACAGAAAATCAGGGCATA",
        "NC_TEST.1:g.[261_262insGA;263_264insAA]",
    );
    assert_eq!(output, "NC_TEST.1:g.263_264insGAAA");
}

#[test]
fn the_demotion_does_not_depend_on_the_authored_order() {
    let forward = assert_padded_preserving(
        "ATACAGAAAATCAGGGCATA",
        "NC_TEST.1:g.[261_262insGA;263_264insAA]",
    );
    let reverse = assert_padded_preserving(
        "ATACAGAAAATCAGGGCATA",
        "NC_TEST.1:g.[263_264insAA;261_262insGA]",
    );
    assert_eq!(forward, reverse);
}

#[test]
fn a_lone_insertion_still_reaches_repeat_notation() {
    // The guard that keeps the demotion from firing when there is no sibling to
    // swallow: a lone insertion inside a tract must still canonicalise to a
    // copy count, which is what `repeated.md` B2 asks for.
    assert_normalizes_preserving(
        "TTTTTTTTTAATATATTTTA",
        "TEMPLATE:g.2_3insTT",
        "TEMPLATE:g.1_9T[11]",
    );
}
