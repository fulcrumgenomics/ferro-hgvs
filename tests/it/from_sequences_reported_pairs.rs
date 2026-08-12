//! `from_sequences` against the externally-reported confluence pairs.
//!
//! These nine rows are the only instances of the #1235 defect class reported
//! from outside the project — found by running a real pipeline, not by a
//! generator — and they are the material this whole surface was built for. They
//! are guarded next door as a **ratchet**: `reported_confluence_pairs` pins how
//! many of them `normalize` converges, and that number is currently **0 under
//! both directions**.
//!
//! This module asks the same question of `from_sequences`, and gets a different
//! kind of answer. Convergence here is not a score to be improved on, so it is
//! asserted outright rather than censused — and the reason is worth stating
//! precisely, because it is the mechanism of the whole design rather than a
//! property of these nine rows.
//!
//! # Why convergence here is structural
//!
//! [`crate::common::cis_apply_oracle::apply`] and
//! `reported_confluence_pairs::every_reported_pair_denotes_one_sequence`
//! establish, independently of the normalizer, that the two spellings of a pair
//! denote one sequence. `Normalizer::to_sequences` returns the window that
//! sequence occupies, which is a function of the *denoted bases* and not of how
//! they were written — so it hands both spellings **byte-identical**
//! `(position, reference, alternate)` triples, at every pad including zero.
//! `from_sequences` then sees one input, and one output follows by arithmetic.
//!
//! That is the point rather than a weakness of the test: erasing the spelling is
//! what this surface is *for*. But it means the load-bearing claim is about
//! `to_sequences`, so
//! [`to_sequences_hands_both_spellings_of_a_pair_one_window`] pins it directly
//! rather than leaving it an unexamined premise — were it ever to stop holding,
//! the convergence assertions here would start failing with no indication of
//! which half had moved.
//!
//! # What the two modules together say, and what they do not
//!
//! Do not read this as "`from_sequences` beats `normalize`". They answer
//! different questions: `normalize` is handed a description and owes rules 2 and
//! 3 (`README.md`), which need the reference; this is handed bases and owes
//! rules 1 and 4. The pairs are exactly the case where being handed a
//! description is the problem, because the difference between the two spellings
//! is precisely what should not survive. Discarding the spelling is not a better
//! rule — it is a different input.
//!
//! So nothing here makes a fix for #1235 unnecessary. A caller who *has* a
//! description and wants it normalized still needs `normalize` to converge, and
//! the ratchet next door is still the number that has to move.

use crate::common::cis_apply_oracle::{apply, normalizer_for};
use crate::reported_confluence_pairs::{REPORTED_PAIRS, TEMPLATE};
use ferro_hgvs::{
    from_sequences_detailed, FromSequencesOptions, HgvsVariant, MockProvider, Normalizer,
    SequencePair, ShuffleDirection,
};

/// Flank asked of `to_sequences`, in bases.
///
/// `128` matches `merge::CANONICAL_PAD`, the room the normalizer's own
/// derivation gives its 3' placement, and is what a caller should ask for. `6`
/// is here because a BAM-derived window is whatever the read covered, and a
/// caller has no obligation to hand over flank it does not have.
///
/// **Neither is 0, and that is a finding rather than a convenience.** An
/// inserted payload placed against the window's first base can only be anchored
/// at `position - 1` — HGVS writes an insertion between two positions — which a
/// zero-pad window does not contain, so `from_sequences` refuses. How much flank
/// clears it depends on the direction: `1420-v2` needs one base under
/// `ThreePrime` and more than one under `FivePrime`, where the payload is placed
/// 5'-most and walks to the start of its ambiguous run. `6` clears every row in
/// both directions. The boundary itself is asserted in
/// [`a_zero_pad_window_refuses_an_insertion_on_its_five_prime_edge`] rather than
/// worked around silently here.
///
/// Each pad is compared **against itself** — spelling A at pad 6 against
/// spelling B at pad 6. Comparing across pads would be asserting against
/// read-dependence, which is a documented property of a window-local derivation
/// rather than a defect: with little flank the 3' placement may genuinely lie
/// outside the bases supplied. `window_is_final` on [`SequencePair`] and
/// `placement_bounded_by_window` on the result are what report that to a caller.
const PADS: [u64; 2] = [6, 128];

/// Both directions. A derivation that converges only under the default 3'
/// direction has not shown the property — 5' is a supported option and reaches
/// the same partitioner, which is the same argument
/// `reported_confluence_pairs::DIRECTIONS` makes.
const DIRECTIONS: [ShuffleDirection; 2] =
    [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime];

/// The window `to_sequences` hands `from_sequences` for one spelling.
///
/// Panics rather than returning `None`. Every input here is a committed row on a
/// committed contig, so a failure to serve one is a defect in this pipeline and
/// must not be filed as a skip — a skip that reads as a pass is the failure mode
/// these corpus tests exist to remove.
fn window(normalizer: &Normalizer<MockProvider>, spelling: &str, pad: u64) -> SequencePair {
    let variant = ferro_hgvs::parse_hgvs(spelling).expect("committed row parses");
    normalizer
        .to_sequences(&variant, pad)
        .unwrap_or_else(|e| panic!("to_sequences({spelling}, pad={pad}): {e}"))
}

/// One spelling through the round trip: description -> sequences -> description.
fn derive(
    normalizer: &Normalizer<MockProvider>,
    spelling: &str,
    pad: u64,
    direction: ShuffleDirection,
) -> (HgvsVariant, bool) {
    let pair = window(normalizer, spelling, pad);
    let options = FromSequencesOptions::default().with_direction(direction);
    let derived = from_sequences_detailed(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &options,
    )
    .unwrap_or_else(|e| panic!("from_sequences({spelling}, pad={pad}): {e}"));
    (derived.variant, derived.placement_bounded_by_window)
}

/// **A zero-pad window refuses an insertion resting on its 5' edge**, and says
/// what to do about it.
///
/// `1420-v2` is the row that exposed this: at pad 0 the window starts at 38,
/// `g.[37dup;41del]` derives to `g.[37_38insA;41del]`, and position 37 is
/// outside it. HGVS writes an insertion between two positions, so there is no
/// other spelling available — and had the window instead started at position 1
/// of the sequence, the anchor would not exist at all.
///
/// Pinned as a refusal, not as a `should_panic`: the message is the deliverable
/// here. `from_sequences`'s own round-trip check catches this shape anyway, but
/// reports it as an internal invariant failure, which sends a caller looking for
/// a bug in ferro rather than adding one base of flank.
#[test]
fn a_zero_pad_window_refuses_an_insertion_on_its_five_prime_edge() {
    let normalizer = normalizer_for(TEMPLATE, ShuffleDirection::ThreePrime);
    let pair = window(&normalizer, "TEMPLATE:g.[37dup;41del]", 0);
    assert_eq!(
        (pair.position, pair.reference.len()),
        (38, 4),
        "the row no longer produces the window this test is about"
    );

    let error = ferro_hgvs::from_sequences(
        &pair.accession,
        pair.position,
        &pair.reference,
        &pair.alternate,
        &FromSequencesOptions::default(),
    )
    .expect_err("an insertion on the 5' edge has no anchor inside the window")
    .to_string();

    for expected in ["5'", "position 37", "flank"] {
        assert!(
            error.contains(expected),
            "the refusal does not name {expected:?}, so it does not tell the caller \
             what to do: {error}"
        );
    }

    // The message's remedy — more 5' flank — has to actually work, or it is a
    // dead end dressed up as advice.
    let padded = window(&normalizer, "TEMPLATE:g.[37dup;41del]", 1);
    ferro_hgvs::from_sequences(
        &padded.accession,
        padded.position,
        &padded.reference,
        &padded.alternate,
        &FromSequencesOptions::default(),
    )
    .expect("one base of 5' flank clears it under ThreePrime, as the refusal promises");

    // …and the message's other half: under `FivePrime` the same one base does
    // **not** clear it, because the payload is placed 5'-most and walks to the
    // start of the `AA` run at 36. This is why the refusal says "more 5' flank"
    // rather than naming a number, and pinning it here is what stops the message
    // from being quietly narrowed back to "one base" by someone reading only the
    // assertion above.
    let five_prime = FromSequencesOptions::default().with_direction(ShuffleDirection::FivePrime);
    ferro_hgvs::from_sequences(
        &padded.accession,
        padded.position,
        &padded.reference,
        &padded.alternate,
        &five_prime,
    )
    .expect_err("under FivePrime one base is not enough — see the refusal's own wording");
}

/// **`to_sequences` hands both spellings of a pair one window.**
///
/// The premise every other assertion in this module rests on, pinned rather than
/// assumed — see the module docs. Compared on all three fields, because agreeing
/// on the bases while disagreeing on `position` would still hand
/// `from_sequences` two different problems.
#[test]
fn to_sequences_hands_both_spellings_of_a_pair_one_window() {
    let normalizer = normalizer_for(TEMPLATE, ShuffleDirection::ThreePrime);
    let mut compared = 0usize;

    for pad in PADS {
        for (label, a, b) in REPORTED_PAIRS {
            let left = window(&normalizer, a, pad);
            let right = window(&normalizer, b, pad);
            compared += 1;
            assert_eq!(
                (&left.position, &left.reference, &left.alternate),
                (&right.position, &right.reference, &right.alternate),
                "{label} (pad={pad}): the two spellings were handed different windows, so \
                 convergence below would no longer follow from the input"
            );
        }
    }

    assert_eq!(
        compared,
        REPORTED_PAIRS.len() * PADS.len(),
        "the reported-pair table did not load in full"
    );
}

/// **Both spellings of every reported pair derive to one description.**
///
/// Asserted at every pad and in both directions. The count of rows evaluated is
/// asserted beside it so a table that failed to load cannot pass as a clean
/// result.
#[test]
fn every_reported_pair_derives_to_one_description() {
    let mut evaluated = 0usize;
    let mut divergent: Vec<String> = Vec::new();

    for direction in DIRECTIONS {
        let normalizer = normalizer_for(TEMPLATE, direction);
        for pad in PADS {
            for (label, a, b) in REPORTED_PAIRS {
                let (left, _) = derive(&normalizer, a, pad, direction);
                let (right, _) = derive(&normalizer, b, pad, direction);
                evaluated += 1;
                if left.to_string() != right.to_string() {
                    divergent.push(format!(
                        "{label} (pad={pad}, {direction:?}): {a} -> {left} but {b} -> {right}"
                    ));
                }
            }
        }
    }

    assert_eq!(
        evaluated,
        REPORTED_PAIRS.len() * PADS.len() * DIRECTIONS.len(),
        "the reported-pair table did not load in full"
    );
    assert!(
        divergent.is_empty(),
        "{} of {evaluated} reported pairs reached two descriptions:\n{}",
        divergent.len(),
        divergent.join("\n")
    );
}

/// **Every derived description denotes the sequence its input denoted.**
///
/// Checked through [`crate::common::cis_apply_oracle::apply`], which converts
/// each member to its SPDI triple and splices the reference. That walk is not
/// the normalizer and not the derivation, so it cannot agree with an output
/// merely because `from_sequences` produced it.
///
/// This is the assertion that survives the module's structural convergence: two
/// spellings reaching one description says nothing about whether that
/// description is *right*, and converging on a wrong sequence would be worse
/// than the divergence the ratchet next door tracks. There is no acceptable
/// non-zero value, so it is asserted rather than counted.
#[test]
fn every_reported_pair_derivation_denotes_its_own_sequence() {
    let mut checked = 0usize;

    for direction in DIRECTIONS {
        let normalizer = normalizer_for(TEMPLATE, direction);
        for pad in PADS {
            for (label, a, b) in REPORTED_PAIRS {
                let expected = apply(TEMPLATE, a).expect("the pair applies");
                for spelling in [a, b] {
                    let (derived, _) = derive(&normalizer, spelling, pad, direction);
                    let rendered = derived.to_string();
                    checked += 1;
                    assert_eq!(
                        apply(TEMPLATE, &rendered).as_deref(),
                        Some(expected.as_str()),
                        "{label} (pad={pad}, {direction:?}): {spelling} -> {rendered} denotes a \
                         different sequence"
                    );
                }
            }
        }
    }

    assert_eq!(
        checked,
        REPORTED_PAIRS.len() * 2 * PADS.len() * DIRECTIONS.len(),
        "the reported-pair table did not load in full"
    );
}

/// **Every derived description is a fixed point of its own derivation.**
///
/// Deriving from the derived form's own window must reproduce it. This is the
/// property that makes "derive now, `normalize` later, or never" safe to offer:
/// a caller that stores a derived string and re-derives it from a fresh read of
/// the same locus gets the same string back.
///
/// Unlike the assertions above, this one is **not** structural. The
/// re-derivation starts from the derived form's own span, which is generally
/// narrower than the input's, so `to_sequences` hands the second pass a
/// different window from the first — which is exactly the case where a
/// derivation could move, and is why this is worth asserting separately.
#[test]
fn every_reported_pair_derivation_is_idempotent() {
    for direction in DIRECTIONS {
        let normalizer = normalizer_for(TEMPLATE, direction);
        for pad in PADS {
            for (label, a, b) in REPORTED_PAIRS {
                for spelling in [a, b] {
                    let (derived, _) = derive(&normalizer, spelling, pad, direction);
                    let rendered = derived.to_string();
                    let (again, _) = derive(&normalizer, &rendered, pad, direction);
                    assert_eq!(
                        again.to_string(),
                        rendered,
                        "{label} (pad={pad}, {direction:?}): {rendered} is not a fixed point"
                    );
                }
            }
        }
    }
}
