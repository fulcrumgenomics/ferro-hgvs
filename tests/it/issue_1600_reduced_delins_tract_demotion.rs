//! A reducing `delins` that lands in a tandem tract must not swallow a sibling.
//!
//! #1600 is #1592's family one step further along, and it is a **wrong-sequence**
//! defect: on `main` a spanning `delins` normalizes to a description that denotes
//! different bases than its input.
//!
//! ```text
//! reference   ("ACGT" x 64) + "GAACAGCAGAAGCGA" + ("ACGT" x 64)   core base 1 at g.257
//! ref 257..272 = GAACAGCAGAAGCGAA
//!
//! g.261_267delinsGCAGAAAC  ->  g.[261del;266_267insCA]
//!   intended  ... G A A C   G C A G A A   A C   G C G A ...
//!   emitted   ... G A A C   G C A G A     C A   A G C G A ...
//! ```
//!
//! The mechanism runs through three passes, each individually right:
//!
//! 1. the sequence-first derivation reaches `g.[261del;263_264insAG;265G>A;267A>C]`,
//!    and its round-trip guard confirms that denotes the input's sequence;
//! 2. `merge_consecutive_edits` merges the insertion with the neighbouring
//!    substitution into `g.265delinsGAA` — a **`delins` that reduces**, which is
//!    #1592's signature;
//! 3. the per-member pipeline reduces it to a two-base insertion at gap 265,
//!    3'-shifts it through the `A` tract to gap 267, and `insertion_to_repeat`
//!    spells the result as a copy count over the whole tract: `g.266_267A[4]`.
//!    That span covers `g.267`, which the sibling `267A>C` claims.
//!
//! Nothing then pulls it back. `clamp_sibling_crossing_shifts` declines a member
//! that grew (`[265,265]` -> `[266,267]`), deliberately, per #1266/#1279;
//! `clamp_sibling_crossing_junctions` needs a junction and `junction_of` reports
//! none for a `Repeat`; and `demote_repeats_spanning_siblings` — the pass that
//! exists to convert such a repeat into a form the clamps *can* see — detected
//! the spanning correctly and then bailed, because its `RepeatSource` match had
//! arms for `Deletion`, `Duplication` and `Insertion` only. A `delins` before-form
//! yielded no source, so the member was skipped by all three.
//!
//! The output is well-formed, its members are disjoint, no warning is raised, and
//! it re-parses, is in bounds and is a fixed point — so all three seam oracles
//! (`FERRO_ASSERT_IDEMPOTENT`, `FERRO_ASSERT_REPARSE`, `FERRO_ASSERT_IN_BOUNDS`)
//! pass on it. None of them compares the denoted sequence, which is the only
//! property it violates, so every assertion here goes through
//! `cis_apply_oracle::apply_with`: it splices the reference from each member's
//! SPDI triple **without** the normalizer. `EquivalenceChecker` cannot serve —
//! `check()` normalizes both operands, so comparing a normalized result against
//! anything returns `Identical` tautologically.
//!
//! # What the pinned strings are, and where their authority comes from
//!
//! Each test below asserts an **exact** output. Two different claims are stacked
//! in that one assertion and they have different warrants, so they are worth
//! separating — reading the string as self-justifying is the trap.
//!
//! * **The question.** Once the reduced insertion lands in the `AA` tract, more
//!   than one rendering denotes the input's bases. Sequence equivalence settles
//!   **validity, not canonicality**: `apply_with` can prove an output wrong, and
//!   can never prove it is the one to emit. That is exactly the defect this
//!   module documents — every wrong output above is *also* well-formed,
//!   in-bounds, re-parsing and idempotent.
//! * **The ruling.** A repeat whose tract spans a sibling's bases is demoted to
//!   the plain edit it grew from, so the sibling stays visible to the clamps;
//!   the member is then bounded rather than allowed to sweep the tract. #1600 is
//!   the `Delins` source of that rule and #1592 is its sibling.
//! * **The authority.** A house rule, **not** a spec clause — no clause in
//!   `assets/hgvs-nomenclature` speaks about one cis member bounding another's
//!   shift, and per `CLAUDE.md` the cross-zone comparison this rests on is
//!   ferro's policy rather than compliance. The governing record is
//!   `rulings[canonical-form-choice-when-both-legal]` (`decided`): among legal
//!   descriptions ferro re-derives from the resulting sequence and emits what
//!   falls out, rather than preserving the input's spelling or the
//!   previously-shipped string.
//!
//! So the **sequence identity** is the load-bearing assertion and fails if the
//! allele stops denoting its input's bases; the string equality is a *record* of
//! the canonical choice under that rule. A change to it is a representation
//! change to declare and re-bless deliberately — not a test to repair by pasting
//! in the new output.

use crate::common::cis_apply_oracle::apply_with;
use crate::common::synthetic::{padded, SyntheticBuilder};
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

/// The core every case here is drawn against.
///
/// Two features make the shape reachable: the `AGCAGAAG` stretch gives the
/// sequence-first derivation an interior alignment match to split on, and the
/// `AA` at `g.266_267` is the two-base tract the reduced insertion shifts into
/// and gets spelled as a copy count over.
const CORE: &str = "GAACAGCAGAAGCGA";

/// Normalize `input` against [`CORE`] in `direction`, assert the output denotes
/// the same bases the input does, and return the rendered output.
///
/// The sequence check runs **before** the caller compares strings, matching
/// `issue_1592_reduced_member_junction_clamp.rs` and for the same reason: a
/// string mismatch asserted first would mean the sequence check never ran, so a
/// re-blessed expectation could silently carry a changed sequence with it. Here
/// the sequence is the assertion and the string is the record.
fn normalized_preserving_in(input: &str, direction: ShuffleDirection) -> String {
    let provider = SyntheticBuilder::genomic(CORE).build();
    let reference = padded(CORE);
    let normalizer = Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    );
    let output = normalizer
        .normalize(&parse_hgvs(input).expect("parse"))
        .expect("normalize")
        .to_string();

    let from_input = apply_with(&provider, &reference, input).expect("input applies");
    let from_output = apply_with(&provider, &reference, &output).unwrap_or_else(|| {
        panic!("`{input}` -> `{output}` has no resulting sequence (overlapping or unconvertible)")
    });
    assert_eq!(
        from_output, from_input,
        "`{input}` -> `{output}` denotes a different sequence"
    );
    output
}

/// [`normalized_preserving_in`] in the default 3' direction.
fn normalized_preserving(input: &str) -> String {
    normalized_preserving_in(input, ShuffleDirection::ThreePrime)
}

#[test]
fn a_reducing_delins_landing_in_a_tract_stops_at_its_sibling() {
    // The reproducer. On `main` this emits `g.[261del;266_267insCA]`, which
    // denotes `...GCAGA CA A GCGA...` where the input denotes
    // `...GCAGAA AC GCGA...` — the `apply_with` comparison in the helper is what
    // fails there, before the string is ever compared.
    let output = normalized_preserving("NC_TEST.1:g.261_267delinsGCAGAAAC");
    assert_eq!(output, "NC_TEST.1:g.[261del;267_268insAC]");
}

#[test]
fn the_bracketed_spelling_of_that_haplotype_agrees() {
    // The same haplotype written as two members. It was already correct — the
    // input-relative weight bound in `canonicalize_from_sequence` declines the
    // re-derivation for this spelling — and is pinned here so the fix is shown to
    // *converge* the two spellings rather than merely to move one of them.
    let output = normalized_preserving("NC_TEST.1:g.[261del;267_268insAC]");
    assert_eq!(output, "NC_TEST.1:g.[261del;267_268insAC]");
}

#[test]
fn a_wider_spanning_delins_over_the_same_block_agrees_too() {
    // Padding the `delins` with unchanged flanks on both sides is a different
    // window and a different input weight, and must not be a different answer.
    let output = normalized_preserving("NC_TEST.1:g.260_268delinsCGCAGAAACG");
    assert_eq!(output, "NC_TEST.1:g.[261del;267_268insAC]");
}

#[test]
fn the_crossing_is_bounded_under_a_five_prime_shuffle_too() {
    // The 5' direction reached a *different* wrong answer on `main`
    // (`g.[261del;265_266insAC]`), so the mirror is asserted rather than assumed.
    let output = normalized_preserving_in(
        "NC_TEST.1:g.261_267delinsGCAGAAAC",
        ShuffleDirection::FivePrime,
    );
    assert_eq!(output, "NC_TEST.1:g.[261del;267_268insAC]");
}

#[test]
fn the_answer_does_not_depend_on_the_authored_member_order() {
    let forward = normalized_preserving("NC_TEST.1:g.[261del;267_268insAC]");
    let reverse = normalized_preserving("NC_TEST.1:g.[267_268insAC;261del]");
    assert_eq!(forward, reverse);
}

#[test]
fn a_lone_reducing_delins_still_reaches_its_repeat_spelling() {
    // The guard against over-clamping, and the reason the fix belongs in
    // `demote_repeats_spanning_siblings` rather than in the tract spelling: with
    // no sibling under the tract, the copy-count form is correct and must
    // survive. `g.265delinsGAA` is the exact member step 2 above builds.
    let output = normalized_preserving("NC_TEST.1:g.265delinsGAA");
    assert_eq!(output, "NC_TEST.1:g.266_267A[4]");
}

#[test]
fn a_sibling_clear_of_the_tract_does_not_bound_it() {
    // Same member, now in an allele — but with the sibling outside the tract, so
    // the demotion must not fire and the repeat spelling must still stand. This
    // is the case that fails if the new `Delins` arm is made unconditional
    // instead of gated on the tract actually spanning a sibling.
    let output = normalized_preserving("NC_TEST.1:g.[261del;265delinsGAA]");
    assert_eq!(output, "NC_TEST.1:g.[261del;266_267A[4]]");

    let distant = normalized_preserving("NC_TEST.1:g.[258A>T;265delinsGAA]");
    assert_eq!(distant, "NC_TEST.1:g.[258A>T;266_267A[4]]");
}
