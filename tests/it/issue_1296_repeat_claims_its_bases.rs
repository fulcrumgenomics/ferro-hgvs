//! A repeat claims the bases under its tract, and must say so.
//!
//! `claims_reference_bases` is how every sibling-repair pass asks "can this
//! member collide with that one". It excluded `NaEdit::Repeat`, so a repeat was
//! invisible as a barrier — and it has no junction either, so the junction rule
//! did not see it. A sibling deletion could 3'-shift straight into the tract:
//!
//! ```text
//! reference  ("ACGT" x 64) + "AAAAAAATAATCGCAACAGAAG" + ("ACGT" x 64)
//! g.[272_273del;274_275insAA]  ->  g.[273_274del;274A[3]]
//!   both members claim 274 — the apply oracle declines the output
//! ```
//!
//! The exclusion was deliberate and tracked by #1269: a repeat that says it
//! claims its bases becomes a barrier everywhere at once, and that had not been
//! measured. #1296 is the defect it left behind, so the measurement is here.
//!
//! It does **not** make the two passes that stood in for it redundant. Disabling
//! `demote_repeats_spanning_siblings` with the predicate corrected fails 12
//! tests including both committed sweeps, and disabling
//! `respell_colliding_duplications` fails two more. Correcting the predicate
//! removes the reason a *third* workaround would have been needed, not the two
//! that exist.

use crate::common::synthetic::assert_padded_preserving;

const CORE: &str = "AAAAAAATAATCGCAACAGAAG";

#[test]
fn a_deletion_does_not_shift_into_a_repeats_tract() {
    // #1296. The deletion of `AC` at 272_273 may 3'-shift to 273_274 on its
    // own — `AC` and `CA` leave the same bases behind — but 274 is where the
    // insertion's repeat starts, so the shift must stop.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[272_273del;274_275insAA]");
    assert_eq!(output, "NC_TEST.1:g.[272_273del;274A[3]]");
}

#[test]
fn the_three_member_form_reaches_the_same_result() {
    // The shape the confluence property actually generated, before it reduced
    // to the pair above: two single-base deletions that merge first.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[272del;273del;274_275insAA]");
    assert_eq!(output, "NC_TEST.1:g.[272_273del;274A[3]]");
}

#[test]
fn a_deletion_clear_of_the_tract_still_shifts() {
    // The guard on the barrier: a repeat bounds a sibling that would land in
    // its tract, and nothing else. This deletion never reaches 274.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[272del;274_275insAA]");
    assert_eq!(output, "NC_TEST.1:g.[272del;274A[3]]");
}
