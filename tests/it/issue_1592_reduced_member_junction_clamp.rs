//! A member that *reduces* to an insertion must still be bounded by its siblings.
//!
//! The two sibling clamps in `normalize_allele` have complementary jobs, and
//! #1592 is the shape that fell between them.
//!
//! * `clamp_sibling_crossing_shifts` governs a member that consumes the bases
//!   under its span, and refuses one that **grew** — deliberately, per
//!   #1266/#1279: an insertion canonicalising to a duplication moves its span
//!   and its junction together, so bounding the span mis-places the copy. It
//!   hands that shape to the junction clamp.
//! * `clamp_sibling_crossing_junctions` governs a member that adds sequence at a
//!   junction — but it needs a junction on *both* snapshots, and `junction_of`
//!   reports one only for `ins` / `dup` / `dupins`.
//!
//! So a member that arrives as a `delins` and normalizes into a `dup` has no
//! `before` junction, is skipped by the pass that owns it, and has already been
//! declined by the pass that could see it. Its junction then sweeps as far as
//! the tract allows — straight over a sibling's bases:
//!
//! ```text
//! reference   ("ACGT" x 64) + "TCCAGCAGATAT" + ("ACGT" x 64)   core base 1 at g.257
//! g.[262_263insAG;264G>T;266T>A]  ->  g.265A[3]
//!   intended  C A G A T A A A T
//!   emitted   C A G A A A T A T     <- the payload crossed the 266T>A
//! ```
//!
//! The output is well-formed, its members are disjoint, no warning is raised,
//! and it is a fixed point that re-parses and is in bounds — so all three seam
//! oracles (`FERRO_ASSERT_IDEMPOTENT`, `FERRO_ASSERT_REPARSE`,
//! `FERRO_ASSERT_IN_BOUNDS`) pass on it. None of them compares the sequence the
//! output denotes with the sequence the input denotes, which is the only
//! property this defect violates. That is why every assertion here goes through
//! `cis_apply_oracle::apply_with`, which splices the reference from each
//! member's SPDI triple **without** the normalizer.
//!
//! The junction a reducing member's payload lands at is not stated anywhere, but
//! it is bounded: it lies inside the member's own footprint `[start - 1, end]`.
//! Taking the edge of that footprint on the side the member travelled from is
//! the conservative before-junction the clamp needs, and it is what the fix
//! supplies.
//!
//! # What the pinned strings below are, and are not
//!
//! Each test asserts an **exact** output string, so it is worth being explicit
//! about which half of that assertion is the guard and where the other half's
//! authority comes from — the two are not the same, and reading the string as
//! self-justifying is the failure mode this note exists to prevent.
//!
//! * **The question.** A member that reduces to an insertion can legally be
//!   spelled at more than one position in a tract: every landing site inside the
//!   swept window denotes bases, and several of them denote the *same* bases.
//!   Sequence equivalence therefore establishes **validity, not canonicality** —
//!   `apply_with` can tell you an output is wrong, and can never tell you it is
//!   the one to emit.
//! * **The ruling.** The junction is bounded at the sibling's edge: it may reach
//!   that edge but not pass it. That is the #1267 rule, stated on
//!   `clamp_sibling_crossing_junctions` in `src/normalize/merge.rs` and applied
//!   symmetrically in both shuffle directions — which is why
//!   `the_crossing_is_bounded_under_a_five_prime_shuffle_too` asserts the mirror
//!   rather than assuming it.
//! * **The authority.** House rule, not a spec clause, and it should be cited as
//!   such. No clause in `assets/hgvs-nomenclature` speaks about one cis member
//!   bounding another's shift; the governing record is
//!   `rulings[canonical-form-choice-when-both-legal]` — *re-derivation governs*,
//!   so among legal spellings ferro emits what falls out of the resulting
//!   sequence subject to the explicit tie-breaks, rather than the input's
//!   spelling or the previously-shipped string. The clamp is what keeps that
//!   re-derivation from selecting a spelling that denotes different bases.
//!
//! So: the **sequence identity** in `normalized_preserving_in` is the load-bearing
//! assertion, and it fails if the allele stops denoting its input's bases. The
//! string equality is a *record* of the canonical choice under the rule above —
//! a change to it is a representation change to be declared and re-blessed
//! deliberately, not a test to repair by pasting the new output.

use crate::common::cis_apply_oracle;
use ferro_hgvs::ShuffleDirection;

/// The core every case here is drawn against. Two lone `A`s two bases apart
/// (`g.265` and `g.267`) inside an `AT` tract is what makes a mis-anchored
/// repeat member land on the wrong one while still rendering legally.
const CORE: &str = "TCCAGCAGATAT";

/// [`common::cis_apply_oracle::normalized_preserving_in`] bound to [`CORE`].
///
/// The helper is shared rather than copied: it encodes a deliberate contract —
/// the `apply_with` sequence check runs before the caller compares strings — and
/// a contract kept in two places can be weakened in one of them.
fn normalized_preserving_in(input: &str, direction: ShuffleDirection) -> String {
    cis_apply_oracle::normalized_preserving_in(CORE, input, direction)
}

/// [`normalized_preserving_in`] in the default 3' direction.
fn normalized_preserving(input: &str) -> String {
    cis_apply_oracle::normalized_preserving(CORE, input)
}

#[test]
fn a_member_reducing_to_an_insertion_stops_at_its_sibling() {
    // The minimal reproducer. The `262_263insAG` merges with the `264G>T` into
    // `264delinsGAT`, which reduces to a pure two-base insertion and 3'-shifts
    // through the `ATATA` tract — past the `266T>A`, which is the crossing.
    let output = normalized_preserving("NC_TEST.1:g.[262_263insAG;264G>T;266T>A]");
    assert_eq!(output, "NC_TEST.1:g.267A[3]");
}

#[test]
fn the_spanning_delins_spelling_reaches_the_same_answer() {
    // The same haplotype written as one `delins`. It takes the sequence-first
    // derivation rather than the per-member merge, so it is a second entrance to
    // the same defect: `g.[259del;265A[3]]` on `main`, naming the lone `A` at
    // 265 instead of the one at 267.
    let output = normalized_preserving("NC_TEST.1:g.258_266delinsCAGCAGATAA");
    assert_eq!(output, "NC_TEST.1:g.[259del;267A[3]]");
}

#[test]
fn the_bracketed_spelling_of_that_haplotype_agrees() {
    // The spelling the confluence proptest draws. It already reached the right
    // answer, and is pinned here so the fix is shown to converge the two
    // spellings rather than merely to move one of them.
    //
    // **How it reaches it has changed, and the old reason is gone.** It used to
    // be that the input-relative weight bound in `canonicalize_from_sequence`
    // declined the derivation for this spelling, so the per-member pipeline
    // answered and the input's own partition survived. That bound is deleted
    // (`rulings[derivation-may-not-be-bounded-by-the-inputs-spelling]`), so this
    // spelling is now re-derived like any other and reaches
    // `g.[259del;267A[3]]` by the same route the spanning spelling above does —
    // #1592's own fix, the reduced-member junction clamp this file is named for.
    // The two spellings therefore agree because they are derived together, not
    // because one of them was exempted.
    let output = normalized_preserving("NC_TEST.1:g.[258del;266_267insAA]");
    assert_eq!(output, "NC_TEST.1:g.[259del;267A[3]]");
}

#[test]
fn the_crossing_is_bounded_under_a_five_prime_shuffle_too() {
    // The clamp is direction-symmetric, so the mirror is asserted rather than
    // assumed. Only the sequence identity is load-bearing here; the string is
    // recorded so a change to it is visible.
    let output = normalized_preserving_in(
        "NC_TEST.1:g.[262_263insAG;264G>T;266T>A]",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(output, "NC_TEST.1:g.267A[3]");
}

#[test]
fn the_answer_does_not_depend_on_the_authored_member_order() {
    let forward = normalized_preserving("NC_TEST.1:g.[262_263insAG;264G>T;266T>A]");
    let reverse = normalized_preserving("NC_TEST.1:g.[266T>A;264G>T;262_263insAG]");
    assert_eq!(forward, reverse);
}

#[test]
fn a_member_with_no_sibling_in_its_way_still_reaches_its_three_prime_most_form() {
    // The guard on the fix: the bound must come from a sibling, not from the
    // reducing member's own shape. Drop the `266T>A` and the same merged
    // `264delinsGAT` is free to travel the whole `ATATA` tract, which it must
    // still do.
    let output = normalized_preserving("NC_TEST.1:g.[262_263insAG;264G>T]");
    assert_eq!(output, "NC_TEST.1:g.268_269dup");
}
