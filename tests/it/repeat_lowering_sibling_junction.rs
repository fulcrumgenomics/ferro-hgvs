//! Lowering a repeat to an insertion at its tract's 3' junction is only sound
//! while no sibling sits inside the tract.
//!
//! `lower_repeat_edits` rewrites a repeat member into the `ins`/`del` it
//! denotes, and places that difference at the tract's 3' end. Its reasoning is
//! that a tract is a whole number of copies of one unit, so it does not matter
//! which end the copies are added at — the resulting sequence is the same.
//!
//! That holds for a lone repeat and fails, silently, when another member adds
//! sequence *inside* the tract:
//!
//! ```text
//! reference  ("ACGT" x 64) + "GATCATAAATTCAGC" + ("ACGT" x 64)
//!            the `A` tract at 263_265 is three bases
//!
//! g.[263_265A[7];264_265insC]
//!   the `C` lands between the tract's own bases, so the four added `A`s
//!   are 5' of it — but lowering puts them at the tract's 3' end instead
//!
//!   intended   ...TAAAAAACAT...
//!   lowered    ...TAACAAAAAT...
//! ```
//!
//! Nothing downstream can catch it. By then the repeat is an ordinary `Ins` and
//! the tract it came from is gone, so the sequence-first pass canonicalizes the
//! wrong sequence and returns a well-formed description of bases the input never
//! denoted — `g.264_265insCAAAA`. That is worse than the overlap it replaced: an
//! overlapping pair is at least visibly malformed, and the apply oracle declines
//! it, whereas this parses, re-parses, is in bounds and is a fixed point.
//!
//! Such a pair overlaps and has no single resulting sequence, so the derivation
//! now refuses rather than picking an order — the same answer
//! `apply_edits_to_window` already gives an insertion interior to a `del`,
//! `delins` or `inv` span (#486).
//!
//! The refusal hands the allele back to the per-member pipeline, which is where
//! the *separate* defect that produced this member set lives: a repeat whose
//! growth exceeds its tract cannot be demoted, so it keeps a span that swallows
//! its sibling's junction. That is tracked and fixed elsewhere; what this module
//! pins is only that the derivation never answers with the wrong bases.

use crate::common::cis_apply_oracle::{apply, normalize};
use crate::common::synthetic::padded;

/// The core the confluence proptest shrank to.
const CORE: &str = "GATCATAAATTCAGC";

/// The input's own resulting sequence, and the normalized output's — where the
/// output has one at all.
fn denotations(input: &str) -> (String, Option<String>) {
    let reference = padded(CORE);
    let from_input = apply(&reference, input).expect("the input must apply");
    let output = normalize(&reference, input);
    (from_input, apply(&reference, &output))
}

#[test]
fn a_repeat_lowered_past_a_sibling_never_denotes_other_bases() {
    // The invariant, stated over the sequence rather than over a string: the
    // output may decline to be well-formed (that is the sibling defect, not
    // this one), but if it denotes anything at all it must denote the input's
    // own bases.
    let input = "TEMPLATE:g.[262_263insAA;263_264insAA;264_265insC]";
    let (intended, denoted) = denotations(input);
    if let Some(denoted) = denoted {
        assert_eq!(
            denoted, intended,
            "`{input}` normalized to a description of different bases"
        );
    }
}

#[test]
fn the_pre_lowered_member_set_is_refused_rather_than_reordered() {
    // The same shape entered directly, which is what the pass actually sees on
    // its second visit: a tract-wide repeat with the sibling's junction inside
    // it. `g.264_265insCAAAA` is the specific wrong answer this guards.
    let input = "TEMPLATE:g.[263_265A[7];264_265insC]";
    let reference = padded(CORE);
    let output = normalize(&reference, input);
    assert!(
        !output.contains("insCAAAA"),
        "the derivation reordered the payloads instead of refusing: {output}"
    );
    // Absence of the wrong answer is not presence of the right one — a refusal
    // that produced some *other* wrong string would satisfy the assert above.
    // So pin the refusal itself: the members come back exactly as authored.
    //
    // Pinned by string rather than by denotation because this input does not
    // apply: its two members claim overlapping territory, so the apply oracle
    // declines the pair. That is the correct answer for a conflicting input,
    // and it is precisely why there is no denotation to compare against — the
    // sibling defect (#1325) is what produces the member set, and it is fixed
    // elsewhere.
    assert_eq!(
        output, input,
        "the refusal must hand back the authored members unchanged"
    );
}

/// A **duplication** sibling reaches the guard through the other arm of
/// `junction_of`, and that arm had no test.
///
/// A `dup` writes at the junction 3' of what it copies (`duplication.md:5`), so
/// it occupies a junction exactly as an `ins` does. The interior case must be
/// refused; the flush case must still lower, because without it a guard that
/// simply refused every `dup` sibling would pass.
///
/// The interior input is pinned by its output string rather than by its
/// denotation: `g.[263_265A[7];264dup]` does not apply, because the two members
/// claim overlapping territory and the apply oracle declines the pair outright.
/// That is the correct answer for a conflicting input and is exactly why there
/// is nothing to compare denotations against.
#[test]
fn a_duplication_sibling_is_treated_like_an_insertion() {
    let reference = padded(CORE);

    let interior = "TEMPLATE:g.[263_265A[7];264dup]";
    let output = normalize(&reference, interior);
    assert!(
        !output.contains("insCAAAA") && !output.contains("A[8]"),
        "`{interior}` lowered past its `dup` sibling instead of refusing: {output}"
    );

    let flush = "TEMPLATE:g.[262_263insAA;263_264insAA;266dup]";
    let (intended, denoted) = denotations(flush);
    assert_eq!(
        denoted.as_deref(),
        Some(intended.as_str()),
        "`{flush}` sits clear of the tract, so it must still derive and denote \
         its own bases"
    );
}

#[test]
fn a_repeat_with_no_sibling_inside_its_tract_still_lowers() {
    // The guard on the refusal: it must key on a junction *inside* the tract,
    // not on the mere presence of a repeat and a sibling. Here the insertion
    // sits at the tract's 3' end, which is adjacency rather than interference,
    // and the derivation must still run.
    let input = "TEMPLATE:g.[262_263insAA;263_264insAA;265_266insC]";
    let (intended, denoted) = denotations(input);
    assert_eq!(
        denoted.as_deref(),
        Some(intended.as_str()),
        "`{input}` must still derive, and must denote its own bases"
    );
    // Sequence preservation alone cannot tell "derived" from "refused" — a
    // refusal preserves the sequence too, so this test would pass on a guard
    // that refused everything. Pin the derived form itself: a single merged
    // member is what lowering produces and what a refusal never would.
    let output = normalize(&padded(CORE), input);
    assert_eq!(
        output, "TEMPLATE:g.[263_265A[7];265_266insC]",
        "`{input}` must still derive — the two insertions become the repeat — \
         rather than being handed back to the per-member pipeline intact"
    );
}
