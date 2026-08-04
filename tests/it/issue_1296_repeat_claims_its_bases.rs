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
//!
//! # The barrier is no longer what decides these two outputs
//!
//! The sequence-first pass used to be locked out of every one of these inputs:
//! the per-member pipeline *promotes* `274_275insAA` to `274A[3]`, and
//! `collect_canonical_edits` refused any group carrying a repeat, so the
//! derivation could never run on the form the pipeline had just produced. With
//! repeats lowered for derivation, it runs — and the first two rows below reach
//! the single substitution the resulting sequence actually states
//! (`g.273C>A`), which is neither of the two members it came from.
//!
//! That makes those two rows **weaker guards for the barrier than they were**:
//! the derivation reaches `273C>A` from the mis-shifted `[273_274del;274A[3]]`
//! too, so they would no longer redden if `claims_reference_bases` stopped
//! seeing repeats. `a_deletion_clear_of_the_tract_still_shifts` is unaffected
//! by that (its answer keeps two members), and the barrier's own coverage is
//! the 12 + 2 tests named above plus
//! `repeat_span_sibling_overlap::no_two_member_allele_normalizes_to_overlapping_members`.

use crate::common::synthetic::assert_padded_preserving;

const CORE: &str = "AAAAAAATAATCGCAACAGAAG";

#[test]
fn a_deletion_does_not_shift_into_a_repeats_tract() {
    // #1296. The deletion of `AC` at 272_273 may 3'-shift to 273_274 on its
    // own — `AC` and `CA` leave the same bases behind — but 274 is where the
    // insertion's repeat starts, so the shift must stop.
    //
    // Re-blessed from `g.[272_273del;274A[3]]`: the sequence-first pass now
    // reads the promoted repeat instead of refusing it, and what these two
    // members denote between them is one base changing. `AACAG` -> `AAAAG`
    // trims to `C` -> `A` at 273.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[272_273del;274_275insAA]");
    assert_eq!(output, "NC_TEST.1:g.273C>A");
}

#[test]
fn the_three_member_form_reaches_the_same_result() {
    // The shape the confluence property actually generated, before it reduced
    // to the pair above: two single-base deletions that merge first.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[272del;273del;274_275insAA]");
    assert_eq!(output, "NC_TEST.1:g.273C>A");
}

#[test]
fn a_deletion_clear_of_the_tract_still_shifts() {
    // The guard on the barrier: a repeat bounds a sibling that would land in
    // its tract, and nothing else. This deletion never reaches 274.
    //
    // Re-blessed from `g.[272del;274A[3]]` for the same reason as the rows
    // above, but this one does not collapse to a single member: `AC` -> `CAA`
    // over 272_273 is two pieces however it is aligned, and the alignment the
    // derivation picks is the one that changes fewer reference columns
    // (insert `C`, then `273C>A` — one substituted column; the two alternatives
    // each cost two).
    //
    // # This row trades a repeat spelling for an insertion spelling
    //
    // Say so plainly, because `open-issues.md:88` explicitly declines to choose
    // between the two: it states that `g.186_191del` and `g.123_191CAG[21]` are
    // *both* correct and that HGVS "require[s] further specifications to
    // indicate when to use which format". So this is a policy-bearing value,
    // not a mechanical consequence, and it is the only such row on the branch.
    //
    // What defends it is that the trade is a **sibling effect, not a demotion
    // of repeat output**. Measured on this same core:
    //
    // ```text
    // g.274_275insAA           -> g.274A[3]              lone: the repeat spelling is still reached
    // g.274A[3]                -> g.274A[3]              lone: and it is a fixed point
    // g.[272del;274_275insAA]  -> g.[271_272insC;273C>A] this row
    // g.[272del;274A[3]]       -> g.[271_272insC;273C>A] and the promoted spelling agrees
    // ```
    //
    // Repeat notation is still produced, and still stable, for this exact
    // tract. What changed is that when a sibling brings the tract into one
    // contiguous run of change, the allele is answered as a whole — and the
    // last two lines are #1235's acceptance criterion met, since those two
    // spellings did not agree before.
    //
    // The open question a reviewer may still want settled: ferro has no stated
    // policy on whether a submitter-named repeat tract should survive a cis
    // merge. If it should, the fix belongs in the piece renderer and this value
    // is wrong; if the allele's sequence governs, this value is right.
    let output = assert_padded_preserving(CORE, "NC_TEST.1:g.[272del;274_275insAA]");
    assert_eq!(output, "NC_TEST.1:g.[271_272insC;273C>A]");
}

#[test]
fn a_repeat_member_does_not_lock_the_derivation_out() {
    // The lockout itself, stated on the two spellings that produce it. The
    // per-member pipeline turns `274_275insAA` into `274A[3]`, so the *output*
    // of one pass became an input the derivation refused — permanently, for
    // every variant the promotion touches.
    //
    // The value pin is the load-bearing half. Agreement alone is satisfied by
    // the lockout too (both spellings converge on the un-derived
    // `[272_273del;274A[3]]`), so it is `g.273C>A` — the form only the
    // derivation can reach — that tells the two states apart.
    let promoted = assert_padded_preserving(CORE, "NC_TEST.1:g.[272_273del;274A[3]]");
    let plain = assert_padded_preserving(CORE, "NC_TEST.1:g.[272_273del;274_275insAA]");
    assert_eq!(
        promoted, plain,
        "the promoted spelling and the one it was promoted from must agree"
    );
    assert_eq!(promoted, "NC_TEST.1:g.273C>A");
}
