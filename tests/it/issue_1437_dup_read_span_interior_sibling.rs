//! #1437 — `detect_insertion_overlaps` keyed a `dup` on the span it **reads**,
//! so a sibling interior to that span was judged by its own edit *kind* rather
//! than by whether anything actually collided.
//!
//! #1411 established the principle: a conflict keys on what a member **writes**,
//! not on what it reads. A `dup` writes its copy directly 3' of the original
//! (`duplication.md:5`), so `g.4_9dup` writes at the junction `9|10` and reads
//! bases 4-9 without touching them. Fix 1(a) applied that to
//! `detect_overlap_conflicts`; branch (b) of `detect_insertion_overlaps` kept
//! the old reading, so the two detectors disagreed.
//!
//! The measured asymmetry on `main`, every sibling strictly interior to the
//! `dup`'s read span and disjoint from its write:
//!
//! ```text
//! g.[4_9dup;6T>C]     ACCEPTED   g.5_6insCTTTTT
//! g.[4_9dup;8T>C]     ACCEPTED   g.7_8insCTTTTT
//! g.[4_9dup;9T>C]     ACCEPTED   g.8_9insCTTTTT
//! g.[4_9dup;6_7insT]  REJECTED   W5002 "g.4_9 has 2 coincident cis-allele edits (dup, ins)"
//! ```
//!
//! Only the insertion was rejected, because a substitution is not an
//! `NaEdit::Insertion` and so branch (b) never looked at it. Nothing about the
//! two writes differed.
//!
//! Resolved by **accepting**: both writes are disjoint, so the composition is
//! unique. That is also the cheap direction for representation stability — the
//! rejected shapes have no shipped normalized form, while the three accepted
//! ones do, so resolving the other way would move output consumers already hold.

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// The 24 bp fixture shared with `strict_rejection_survives_normalization` and
/// `issue_1307_terminal_dup_respell`.
const SEQUENCE: &str = "TTTTTTTTTAATATATTTTAATAC";

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_TEST.1", SEQUENCE.to_string());
    provider
}

fn strict_accepts(input: &str) -> bool {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::with_config(provider(), NormalizeConfig::strict())
        .normalize(&variant)
        .is_ok()
}

fn lenient(input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::new(provider())
        .normalize(&variant)
        .expect("lenient normalization must not reject")
        .to_string()
}

/// The headline: an insertion interior to a `dup`'s read span is accepted, like
/// the substitutions beside it always were.
///
/// Written as one loop over all four siblings so the property under test is the
/// *agreement*, not four independent verdicts. A regression that re-rejected the
/// insertion would still leave the three substitutions passing if they were
/// asserted separately.
#[test]
fn every_sibling_interior_to_a_dups_read_span_is_accepted() {
    for input in [
        "NC_TEST.1:g.[4_9dup;6T>C]",
        "NC_TEST.1:g.[4_9dup;8T>C]",
        "NC_TEST.1:g.[4_9dup;9T>C]",
        "NC_TEST.1:g.[4_9dup;6_7insT]",
    ] {
        assert!(
            strict_accepts(input),
            "`{input}`: the `dup` writes at the junction 9|10 and the sibling \
             writes inside 4-9, so the two do not collide and strict mode must \
             accept — whatever kind the sibling is"
        );
    }
}

/// The three substitution answers are unchanged, pinned literally.
///
/// This is the stability half. The argument for resolving #1437 by accepting
/// rather than by rejecting is that these three are *shipped* output while the
/// insertion's rejection is not, so a regression here is the expensive one.
#[test]
fn the_substitution_answers_are_unchanged() {
    for (input, expected) in [
        ("NC_TEST.1:g.[4_9dup;6T>C]", "NC_TEST.1:g.5_6insCTTTTT"),
        ("NC_TEST.1:g.[4_9dup;8T>C]", "NC_TEST.1:g.7_8insCTTTTT"),
        ("NC_TEST.1:g.[4_9dup;9T>C]", "NC_TEST.1:g.8_9insCTTTTT"),
    ] {
        assert_eq!(lenient(input), expected, "`{input}` must be unchanged");
    }
}

/// A `repeat` is exempted with the `dup`, not after it.
///
/// The per-member pipeline chooses between the two spellings by axis — the CDS
/// axis rewrites `c.5_6insAA` to `c.5_6dup`, the genomic axis rewrites the same
/// shape to `g.1005_1006A[4]` — so exempting one and not the other would make
/// the detector spelling-sensitive, which is the defect the junction
/// registration exists to prevent.
#[test]
fn a_repeat_gets_the_same_treatment_as_a_dup() {
    // `10_11` is `AA`; the tract spelled as a repeat reads 10-11 and writes at
    // the junction after it, exactly as `10_11dup` does.
    let repeat = "NC_TEST.1:g.[10_11A[3];10_11insT]";
    let dup = "NC_TEST.1:g.[10_11dup;10_11insT]";
    // Both accepted, not merely equal. Parity alone is satisfied by both being
    // *rejected*, which is precisely the state this change moves away from — so
    // an `assert_eq!` on the two booleans would keep passing if the exemption
    // regressed, and would be testing nothing about the fix.
    assert!(
        strict_accepts(dup),
        "a `dup` writes at the junction 3' of the span it reads, so a sibling \
         interior to that read span collides with nothing and strict must accept"
    );
    assert!(
        strict_accepts(repeat),
        "the `repeat` spelling of the same shape must get the same verdict as \
         the `dup` one, or the detector is sensitive to which spelling the \
         per-member pipeline happened to choose"
    );
}

/// The discriminating half: a span edit that really does write across its range
/// still conflicts with an interior insertion.
///
/// An `inv` rewrites the bases it spans in place, so an interior junction has a
/// meaningful position and the pair is genuinely ambiguous. A blanket "spans
/// never conflict with interior insertions" would pass every assertion above and
/// break this.
#[test]
fn an_inversion_still_conflicts_with_an_interior_insertion() {
    assert!(
        !strict_accepts("NC_TEST.1:g.[4_9inv;6_7insT]"),
        "an inversion writes across its own span, so an interior insertion is \
         still a conflict"
    );
}

/// Two `dup`s writing at one junction still conflict — that is branch (a), and
/// removing the read span from branch (b) must not disturb it.
#[test]
fn two_duplications_at_one_junction_still_conflict() {
    assert!(
        !strict_accepts("NC_TEST.1:g.[4_9dup;7_9dup]"),
        "both duplications write at the junction 9|10 with no defined order, \
         which branch (a) reports"
    );
}
