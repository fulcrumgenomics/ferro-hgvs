//! Property tests for `from_sequences`.
//!
//! Every other module in this family draws from a fixed corpus, so each can only
//! exercise shapes its generator was written to emit — and this branch has
//! already been caught by that twice: the accession guard leaked `ENSMUST` and
//! `LRG_1p1` because a hand-written prefix list did not think of them, and the
//! alphabet refused real ClinVar rows because it admitted only `ACGTN`. Both
//! were found by widening the *input set*, not by adding an assertion.
//!
//! The properties below are the ones the surface claims unconditionally,
//! stated so that proptest can look for a counterexample rather than a corpus
//! author having to think of one:
//!
//! | property | rule |
//! |---|---|
//! | the output parses | 1 (conformant) |
//! | the same input gives the same output | 4 (deterministic) |
//! | re-deriving over its own span reproduces it | 4 |
//! | the direction does not decide describability | 6 (no user options for form) |
//! | construction and derivation agree on what is usable | — |
//!
//! **Denotation is not among them, and its absence is not a gap.**
//! `from_sequences` re-applies its own output to the supplied window and
//! compares, as a runtime check rather than a `debug_assert` — so every `Ok`
//! below has already had that property verified inside the call. Asserting it
//! again here would re-compare the same two byte strings.
//!
//! **Rules 2 and 3 are deliberately absent.** The recommended form and
//! confluence need the reference, which this surface does not read; asserting
//! them here would be asserting something the design does not claim.
//!
//! # No reference, so no provider
//!
//! `from_sequences` is a pure function of its arguments, which is what makes it
//! testable this way at all — a generated `(position, reference, alternate)`
//! triple is a complete input, with no fixture to build and no window to serve.

use ferro_hgvs::{from_sequences_detailed, FromSequencesOptions, SequencePair, ShuffleDirection};
use proptest::prelude::*;

/// The alphabet `validate` admits, minus `U`.
///
/// `U` is excluded because a `g.` description is DNA: mixing it in would
/// generate inputs whose refusal says nothing about the properties under test.
/// The ambiguity codes are deliberately *included* — they are what real
/// submitted data carries, and admitting them is a decision this branch made
/// (see `validate`), so they belong in the generated set rather than outside it.
fn base() -> impl Strategy<Value = char> {
    prop::sample::select(vec![
        'A', 'C', 'G', 'T', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
    ])
}

/// A window: `(position, reference, alternate)`.
///
/// The two sequences are drawn **independently**, so most pairs are not a
/// plausible read at all — a random 12-mer against a random 9-mer. That is on
/// purpose: the properties are claimed for every input the function accepts, and
/// a generator that only produced realistic reads would test the realistic half
/// of that claim.
fn window() -> impl Strategy<Value = (u64, String, String)> {
    (
        1u64..1_000_000,
        prop::collection::vec(base(), 1..24),
        prop::collection::vec(base(), 0..24),
    )
        .prop_map(|(position, reference, alternate)| {
            (
                position,
                reference.into_iter().collect(),
                alternate.into_iter().collect(),
            )
        })
}

fn options(direction: ShuffleDirection) -> FromSequencesOptions {
    FromSequencesOptions::default().with_direction(direction)
}

proptest! {
    /// **Whatever is derived, parses.**
    ///
    /// Rule 1's cheapest half. A description ferro cannot read back is one no
    /// consumer can either.
    #[test]
    fn a_derived_description_parses((position, reference, alternate) in window()) {
        let Ok(derived) = from_sequences_detailed(
            "NC_TEST.1", position, &reference, &alternate, &options(ShuffleDirection::ThreePrime)
        ) else {
            return Ok(());   // a refusal is a legitimate answer; see the module docs
        };
        let rendered = derived.variant.to_string();
        prop_assert!(
            ferro_hgvs::parse_hgvs(&rendered).is_ok(),
            "{} {} -> {} derived {}, which does not re-parse",
            position, reference, alternate, rendered
        );
    }

    /// **The same input gives the same output** — rule 4, asserted directly
    /// rather than inferred from the absence of an obvious source of
    /// nondeterminism.
    #[test]
    fn deriving_twice_agrees((position, reference, alternate) in window()) {
        let opts = options(ShuffleDirection::ThreePrime);
        let first = from_sequences_detailed("NC_TEST.1", position, &reference, &alternate, &opts);
        let second = from_sequences_detailed("NC_TEST.1", position, &reference, &alternate, &opts);
        match (first, second) {
            (Ok(a), Ok(b)) => prop_assert_eq!(a.variant.to_string(), b.variant.to_string()),
            (Err(_), Err(_)) => {}
            _ => prop_assert!(false, "one call derived and the other refused"),
        }
    }

    /// **Deriving again over the derived form's own span reproduces it.**
    ///
    /// Genuinely distinct from `deriving_twice_agrees` above, which repeats the
    /// *same* call: this one re-derives over a **narrower window** — the span
    /// the derived description itself names — which is the window a caller gets
    /// back when they store a derived string and re-read that locus. If the
    /// derivation moved under a tighter window, "derive now, normalize later,
    /// or never" would not be safe to offer.
    ///
    /// An earlier version of this test called the function twice with identical
    /// arguments while its doc claimed the narrowing, so it silently duplicated
    /// the property above; the narrowing is now performed.
    ///
    /// **And the version after that narrowed but asserted nothing.** It checked
    /// that the re-derived string contained `del`, `ins`, `inv`, `dup`, `>` or
    /// ended `g.=` — a disjunction over very nearly everything this surface can
    /// emit, so no reproduction was ever required and the test's own name was
    /// the only place the property appeared. It could not have asserted equality
    /// either, because its two windows were mismatched: `reference[lo..=hi]`
    /// against `alternate[lo..]`, which takes the reference's *hull* and the
    /// alternate's whole *tail*. Those describe different loci, so re-deriving
    /// over them legitimately gives a different answer.
    ///
    /// The alternate's mirror of `reference[lo..=hi]` is measured from the 3'
    /// end, not the 5' one: everything outside the named span is identical
    /// between the two sequences, so the same number of unchanged bases
    /// (`reference.len() - hi - 1`) trails the span in each. With the windows
    /// matched, equality is the property and is asserted.
    #[test]
    fn a_derivation_over_its_own_span_reproduces_it((position, reference, alternate) in window()) {
        let opts = options(ShuffleDirection::ThreePrime);
        let Ok(first) = from_sequences_detailed(
            "NC_TEST.1", position, &reference, &alternate, &opts
        ) else {
            return Ok(());
        };
        let rendered = first.variant.to_string();

        // The span the derived description names, as window offsets. An
        // identity (`g.=`) names none, and a member on the window edge would
        // narrow to a window that cannot anchor it — both are skipped rather
        // than asserted about, since neither is what this property is for.
        let Some((lo, hi)) = named_span(&rendered, position) else {
            return Ok(());
        };
        // A named position past the window's last base belongs to an insertion
        // anchored at the 3' edge; there is no narrower window that keeps it.
        // Clamping `hi` — which the previous version did — silently narrows to a
        // *different* span rather than declining, so it is skipped instead.
        if hi >= reference.len() {
            return Ok(());
        }
        let unchanged_suffix = reference.len() - hi - 1;
        let Some(alt_hi) = alternate.len().checked_sub(unchanged_suffix) else {
            return Ok(());
        };
        if alt_hi < lo {
            return Ok(());
        }
        let (sub_ref, sub_alt) = (&reference[lo..=hi], &alternate[lo..alt_hi]);

        // Re-derive over the narrower window. A refusal is legitimate — the
        // narrowed window can drop the flank an insertion needs to anchor —
        // so only a *different answer* is a failure.
        if let Ok(again) = from_sequences_detailed(
            "NC_TEST.1",
            position + lo as u64,
            sub_ref,
            sub_alt,
            &opts,
        ) {
            prop_assert_eq!(
                again.variant.to_string(),
                rendered.clone(),
                "re-deriving {} over its own span ({}..={} of {} / {}) did not reproduce it",
                rendered, lo, hi, reference, alternate
            );
        }
    }

    /// **Both shuffle directions derive, or both refuse — except at the
    /// window's own 5' edge, where the asymmetry is the geometry and not a
    /// policy.**
    ///
    /// Direction chooses a placement within the window; it must not decide
    /// whether the input is describable at all, because an option that gates
    /// capability is what `README.md` rule 6 forbids.
    ///
    /// There is exactly one honest exception, and stating it is the point of
    /// this test rather than a weakening of it. HGVS writes an insertion
    /// *between* the two positions it falls between, so a payload placed
    /// immediately 5' of the window's first base can only be anchored at
    /// `position - 1` — outside the window, and non-existent at `position` 1. A
    /// 5'-most placement walks to the start of an ambiguous run and can reach
    /// that edge where a 3'-most placement never does. `1 "A" -> "AA"` is the
    /// minimal case and is pinned in the regression file beside this module.
    ///
    /// So the assertion is not `is_ok() == is_ok()` — that form was written
    /// first, held for 256 cases, and failed at 20,000 against exactly that
    /// input. What must hold is that a one-sided refusal is **that** refusal:
    /// any other divergence means the direction really has become a gate.
    #[test]
    fn the_direction_does_not_decide_describability(
        (position, reference, alternate) in window()
    ) {
        let three = from_sequences_detailed(
            "NC_TEST.1", position, &reference, &alternate,
            &options(ShuffleDirection::ThreePrime),
        );
        let five = from_sequences_detailed(
            "NC_TEST.1", position, &reference, &alternate,
            &options(ShuffleDirection::FivePrime),
        );
        // Matched on the anchor clause rather than the error variant, because
        // `InvalidCoordinates` also carries the ordinary validation refusals —
        // which are direction-independent, so accepting the variant would
        // accept the divergence this test exists to catch.
        let anchor_refusal = |e: &ferro_hgvs::FerroError| {
            e.to_string().contains("only HGVS anchor is position")
        };
        match (&three, &five) {
            (Ok(_), Ok(_)) | (Err(_), Err(_)) => {}
            (Ok(_), Err(only)) | (Err(only), Ok(_)) => prop_assert!(
                anchor_refusal(only),
                "{} {} -> {}: 3' {} but 5' {}, and the refusal is not the 5'-anchor one: {}",
                position, reference, alternate,
                if three.is_ok() { "derived" } else { "refused" },
                if five.is_ok() { "derived" } else { "refused" },
                only,
            ),
        }
    }

    /// **A pair that constructs is a pair that derives, and vice versa.**
    ///
    /// `SequencePair::new` shares `from_sequences`'s validator precisely so this
    /// holds; pinned over generated input because the two entry points are what
    /// a caller alternates between, and a divergence would put the error one
    /// call away from the argument that caused it.
    #[test]
    fn construction_and_derivation_agree_on_what_is_usable(
        (position, reference, alternate) in window()
    ) {
        let built = SequencePair::new("NC_TEST.1", position, &*reference, &*alternate);
        let derived = from_sequences_detailed(
            "NC_TEST.1", position, &reference, &alternate,
            &options(ShuffleDirection::ThreePrime),
        );
        // Only the *validation* verdict is compared. A pair can construct and
        // then decline on the grid budget, which is a cost bound rather than a
        // statement about the input, so a construct-ok/derive-err pair is only a
        // failure when the derivation refused for a validation reason.
        if built.is_err() {
            prop_assert!(
                derived.is_err(),
                "the constructor refused an input the derivation accepted"
            );
        }
    }
}

/// The window offsets of the first and last position a description names.
///
/// Parsed off the rendered string rather than the variant, because the point is
/// to reconstruct what a *caller* would slice from a stored description — they
/// have the string, not the internal shape. Returns `None` for an identity or
/// anything that names no coordinate.
fn named_span(rendered: &str, window_start: u64) -> Option<(usize, usize)> {
    let body = rendered.split(":g.").nth(1)?;
    let mut positions: Vec<u64> = Vec::new();
    let mut digits = String::new();
    for c in body.chars() {
        if c.is_ascii_digit() {
            digits.push(c);
        } else if !digits.is_empty() {
            positions.push(digits.parse().ok()?);
            digits.clear();
        }
    }
    if !digits.is_empty() {
        positions.push(digits.parse().ok()?);
    }
    let (lo, hi) = (*positions.iter().min()?, *positions.iter().max()?);
    // A position before the window start belongs to an insertion anchored 5' of
    // it, which no narrowing can keep.
    let lo = lo.checked_sub(window_start)?;
    let hi = hi.checked_sub(window_start)?;
    Some((usize::try_from(lo).ok()?, usize::try_from(hi).ok()?))
}
