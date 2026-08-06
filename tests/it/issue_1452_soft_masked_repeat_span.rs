//! #1452 — a soft-masked repeat tract made the two spellings of one repeat
//! denote different sequences again.
//!
//! This is the #1431 defect reopened by a different door. #1431 taught
//! `hgvs_to_spdi` that a single-position anchor (`g.259CAG[5]`) names the whole
//! tandem run, not one unit of it, by searching the reference for that run. The
//! search compares raw reference bytes against a unit that has already been
//! upper-cased for SPDI output, so on a **soft-masked** (lowercase) tract it
//! matched nothing.
//!
//! What makes it worth a regression test rather than a one-line fold is that it
//! did **not** surface as an error. The no-run branch falls back to the
//! unit-wide span at the anchor — deliberately, so the caller's
//! "does not match repeat unit" diagnostic reports a genuine mismatch — and the
//! caller then upper-cases those bases before checking them. One unit always
//! matches one unit, so the check passed and a **truncated** triple came out:
//!
//! | spelling | uppercase `CAGCAGCAG` | soft-masked `cagcagcag` |
//! |---|---|---|
//! | `g.259CAG[5]` | `258:CAGCAGCAG:CAGCAGCAGCAGCAG` | `258:CAG:CAGCAGCAGCAGCAG` |
//! | `g.259_267CAG[5]` | `258:CAGCAGCAG:CAGCAGCAGCAGCAG` | `258:CAGCAGCAG:CAGCAGCAGCAGCAG` |
//!
//! The masked anchored form applies to seven copies, not five: the two tract
//! units the truncated span left untouched survive the replacement. Both
//! spellings converted, both parse, both are in bounds and both are fixed
//! points — so, exactly as in #1431, none of `FERRO_ASSERT_REPARSE`,
//! `FERRO_ASSERT_IN_BOUNDS` or `FERRO_ASSERT_IDEMPOTENT` can see it.
//!
//! Measured over anchored conversions of `A`/`AT`/`CAG`/`GATC` tracts of 2..=5
//! copies at every anchor position and two target counts, 560 rows in all. 336
//! of those anchor out of phase with the tract and are refused — identically in
//! both case conventions, which the test below pins. Of the 224 that remain,
//! 112 are soft-masked and 112 are not, and the split was total: **every
//! soft-masked row disagreed with its range spelling; no unmasked row did**.
//! After the fix, neither does.
//!
//! **The issue's stated cause is not the one that was firing.** #1452 named two
//! case-sensitive comparisons, the unit-match check
//! (`del_str != unit_str.repeat(pre_count)`) and the tract search. The
//! unit-match check is *already* case-insensitive — both of its operands go
//! through `apply_alphabet`, which upper-cases — so the explicit-range spelling
//! was never affected and no failure here is a conversion *error*. Only the
//! tract search was folding-blind, and its symptom is a wrong answer rather
//! than a refusal. The issue's warning that fixing one of the two comparisons
//! alone would be worse than fixing neither therefore does not apply: there is
//! one site, and folding it is what makes the two spellings agree.
//!
//! Reference FASTAs are routinely soft-masked and masked bytes reach this
//! codebase's normalization from ordinary input
//! (`issue_1318_soft_masked_delins.rs`), so this is not a synthetic condition.

use ferro_hgvs::spdi::hgvs_to_spdi;
use ferro_hgvs::{parse_hgvs, MockProvider};

/// 256 `N`s, so core position 1 is `g.257`. Borrowed from
/// `issue_1431_single_position_repeat_anchor.rs`: `N` matches no unit used
/// here, so the flanks can never look like part of a tract.
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{PAD}{core}{PAD}")
}

fn provider(sequence: &str) -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("TEMPLATE", sequence.to_string());
    provider
}

/// The SPDI triple, rendered.
fn spdi(sequence: &str, input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    hgvs_to_spdi(&variant, &provider(sequence))
        .unwrap_or_else(|e| panic!("`{input}` must convert to SPDI: {e}"))
        .to_string()
}

/// The conversion error, for the cases that must still be refused.
fn spdi_err(sequence: &str, input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    match hgvs_to_spdi(&variant, &provider(sequence)) {
        Ok(triple) => panic!("`{input}` must not convert, got {triple}"),
        Err(e) => e.to_string(),
    }
}

/// A 3-copy `CAG` tract at 259-267, soft-masked; `GG` flanks keep it isolated.
const MASKED_CAG: &str = "GGcagcagcagGG";
/// The byte-identical tract with no masking, as the control.
const UPPER_CAG: &str = "GGCAGCAGCAGGG";

/// The headline: on a soft-masked tract the two spellings of one repeat still
/// denote one sequence.
///
/// Before the fix the anchored form gave `258:CAG:CAGCAGCAGCAGCAG` here — a
/// one-unit deletion span, so the replacement left the tract at seven copies
/// while the description said five.
#[test]
fn the_two_spellings_agree_on_a_soft_masked_tract() {
    let seq = padded(MASKED_CAG);
    let expected = "TEMPLATE:258:CAGCAGCAG:CAGCAGCAGCAGCAG";
    assert_eq!(spdi(&seq, "TEMPLATE:g.259CAG[5]"), expected);
    assert_eq!(spdi(&seq, "TEMPLATE:g.259_267CAG[5]"), expected);
}

/// Masking is not information: the triple must not depend on it.
///
/// Stated as an equality against the uppercase twin rather than as a second
/// literal, because the property is "the reference's case convention does not
/// reach the output", and a literal on each side would pass even if both had
/// moved together.
#[test]
fn a_soft_masked_tract_converts_as_its_uppercase_twin() {
    let masked = padded(MASKED_CAG);
    let upper = padded(UPPER_CAG);
    for input in [
        "TEMPLATE:g.259CAG[5]",
        "TEMPLATE:g.259_267CAG[5]",
        // A count *below* the reference count contracts the tract; the tract
        // still has to be found at its full extent for that to be right.
        "TEMPLATE:g.259CAG[1]",
        "TEMPLATE:g.259_267CAG[1]",
    ] {
        assert_eq!(
            spdi(&masked, input),
            spdi(&upper, input),
            "`{input}` must not depend on the reference's masking"
        );
    }
}

/// The same, for a single-base unit — the shape where the truncation is easiest
/// to mistake for a correct answer, because a one-base span looks like exactly
/// what a one-position anchor should produce.
#[test]
fn a_soft_masked_single_base_tract_spans_the_whole_run() {
    let seq = padded("GGaaaGG");
    let expected = "TEMPLATE:258:AAA:AAAAA";
    assert_eq!(spdi(&seq, "TEMPLATE:g.259A[5]"), expected);
    assert_eq!(spdi(&seq, "TEMPLATE:g.259_261A[5]"), expected);
}

/// The discriminating case: folding the window must not make a span that is
/// genuinely not the declared repeat convert anyway.
///
/// The obvious over-general fix — relax the unit-match check instead of the
/// search — would accept this, because the fallback span is one unit wide and
/// one unit always matches itself. `ATGCAT` is divisible by `AT` and is not
/// `AT[3]`, so it must still be refused, in either case convention.
#[test]
fn a_non_repeat_span_is_still_refused_when_soft_masked() {
    for core in ["GGatgcatGG", "GGATGCATGG"] {
        let seq = padded(core);
        let err = spdi_err(&seq, "TEMPLATE:g.259_264AT[3]");
        assert!(
            err.contains("does not match repeat unit"),
            "`{core}` must still be refused as a non-repeat span, got: {err}"
        );
    }
}

/// The second discriminating case: an anchor that lands out of phase with the
/// tract is refused, and is refused *identically* whether or not the reference
/// is masked.
///
/// This is the pre-existing behaviour the fix must leave alone. `g.260CAG[5]`
/// points at the second base of the run, where no `CAG` copy begins, so the
/// search finds nothing and the unit-wide fallback span (`AGC`) is correctly
/// judged a mismatch. Pinning it here keeps a later widening of the search from
/// quietly inventing a tract at an anchor the spec's start-only format does not
/// name.
#[test]
fn an_out_of_phase_anchor_is_refused_in_both_case_conventions() {
    let masked = spdi_err(&padded(MASKED_CAG), "TEMPLATE:g.260CAG[5]");
    let upper = spdi_err(&padded(UPPER_CAG), "TEMPLATE:g.260CAG[5]");
    assert_eq!(masked, upper, "masking must not change the diagnostic");
    // The text changed in #1492: the decline now names the anchor the caller
    // wrote (260) instead of the unit-wide fallback span (260-262) the tract
    // search invented. The *behaviour* this test exists for — refused, and
    // refused identically in both case conventions — is unchanged.
    assert!(
        masked.contains("no CAG repeat is anchored at TEMPLATE:260"),
        "unexpected diagnostic: {masked}"
    );
    assert!(
        !masked.contains("260-262"),
        "the decline must not quote a span the caller never wrote: {masked}"
    );
}

/// The whole census the module doc quotes, as a standing check: no anchored
/// conversion of a tandem tract may disagree with the range spelling of that
/// same tract, on either case convention.
///
/// Written as a sweep rather than as more literals because the defect was
/// uniform across units and copy counts — it was 112 rows, not one shape — and
/// a per-shape test would have been satisfied by a fix that only handled the
/// unit it was written for.
#[test]
fn no_anchored_conversion_disagrees_with_its_range_spelling() {
    let mut compared = 0usize;
    let mut disagreed = Vec::new();

    for unit in ["A", "AT", "CAG", "GATC"] {
        for copies in 2..=5usize {
            let tract = unit.repeat(copies);
            for masked in [false, true] {
                let core = if masked {
                    format!("GG{}GG", tract.to_ascii_lowercase())
                } else {
                    format!("GG{tract}GG")
                };
                let seq = padded(&core);
                let tract_start = 259u64;
                let tract_end = tract_start + tract.len() as u64 - 1;

                for n_post in [1usize, copies + 2] {
                    let ranged = format!("TEMPLATE:g.{tract_start}_{tract_end}{unit}[{n_post}]");
                    let expected = spdi(&seq, &ranged);
                    // Only in-phase anchors name a tract start; an out-of-phase
                    // one is refused, which the test above pins.
                    for anchor in (tract_start..=tract_end).step_by(unit.len()) {
                        let anchored = format!("TEMPLATE:g.{anchor}{unit}[{n_post}]");
                        compared += 1;
                        let actual = spdi(&seq, &anchored);
                        if actual != expected {
                            disagreed.push(format!(
                                "{core}: {anchored} -> {actual}, {ranged} -> {expected}"
                            ));
                        }
                    }
                }
            }
        }
    }

    assert!(
        disagreed.is_empty(),
        "{} anchored conversions disagree with their range spelling:\n{}",
        disagreed.len(),
        disagreed.join("\n")
    );
    // Second, and only after the substantive assertion: a floor on the corpus,
    // so a generator change that quietly stops emitting cases cannot turn this
    // into a green test that compares nothing.
    assert_eq!(compared, 224, "corpus size changed; re-read the census");
}
