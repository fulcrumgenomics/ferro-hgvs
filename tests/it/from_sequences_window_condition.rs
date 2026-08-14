//! The exact condition under which two windows derive the same description.
//!
//! `from_sequences` reads no reference, so nothing it emits may shift outside
//! the bases the caller supplied. That makes the answer a function of the
//! window as well as of the variant, and the question this module settles is
//! *precisely when* the window stops mattering.
//!
//! # The condition
//!
//! Call the **ambiguity interval** of a variant the closed range over which its
//! change can be placed — operationally, from its 5'-most placement to its
//! 3'-most, which for a tandem tract is the tract itself (plus, for an inserted
//! payload, the base immediately 5' of it, since that is where a 5'-most
//! insertion anchors).
//!
//! > **Two windows that both contain the whole ambiguity interval derive the
//! > same description. A window that cuts it derives a different one, placed at
//! > the window's own edge.**
//!
//! Both halves are pinned below, and the second is what makes the first worth
//! stating: if truncating windows happened to agree anyway, "spanning" would be
//! an unnecessary hypothesis rather than the actual boundary.
//!
//! # Why this is not a defect
//!
//! A truncated window's answer is **still correct**. `g.14del` and `g.15del`
//! over a 4-base `A` run denote the same bases and carry the same canonical
//! SPDI; the truncated one is simply not the 3'-most spelling. So the cost of a
//! short read is a **rule 2** miss (recommended form) and a **rule 3** miss
//! (confluence across inputs) — never rule 1 or rule 4, which
//! `crate::from_sequences` owes unconditionally and which hold here in every
//! row below. Rules 2 and 3 are exactly the two the design assigns to
//! `normalize`, because both need the reference.
//!
//! This was worth measuring rather than asserting: the module doc on
//! `from_sequences` originally claimed the window-edge flag fires "exactly and
//! only when the answer is read-dependent", and the `A4` rows below disprove it
//! — window `9-15` sits flush against the tract's end, is flagged, and is
//! nonetheless the same answer a whole-sequence derivation gives.
//!
//! # No provider
//!
//! Every assertion here calls the free function, which holds no reference. The
//! module is therefore hermetic and needs no fixture, which is the point of the
//! surface it tests.

use ferro_hgvs::{from_sequences_detailed, FromSequencesOptions};
use std::collections::BTreeSet;

/// A synthetic contig carrying four tandem tracts of differing unit length and
/// copy number, separated by unique sequence so no tract's ambiguity interval
/// touches another's.
const SEQUENCE: &str = concat!(
    "GGATTACAGGC", // 1-11    unique
    "AAAA",        // 12-15   homopolymer, unit 1
    "GCC",         // 16-18   unique
    "TGTGTGTG",    // 19-26   dinucleotide, unit 2
    "AGG",         // 27-29   unique
    "CATTAGCCT",   // 30-38   unique
    "AAAAAAAA",    // 39-46   homopolymer, unit 1, longer
    "CGCGCG",      // 47-52   dinucleotide, unit 2
    "TTAGC",       // 53-57   unique
);

/// `(label, first, last, unit)` — 1-based inclusive tract bounds.
///
/// The bounds are asserted against `SEQUENCE` before use by
/// [`the_tract_table_matches_the_sequence`]. A miscounted offset here would
/// silently test a different tract and still pass, which is how the first
/// version of this sweep reported four spurious divergences.
const TRACTS: &[(&str, u64, u64, &str)] = &[
    ("A4", 12, 15, "A"),
    ("TG", 19, 26, "TG"),
    ("A8", 39, 46, "A"),
    ("CG", 47, 52, "CG"),
];

/// Bases `[lo, hi]`, 1-based inclusive.
fn at(lo: u64, hi: u64) -> String {
    SEQUENCE[(lo - 1) as usize..hi as usize].to_string()
}

/// Derive over a window, returning `(description, placement_bounded_by_window)`.
fn derive(lo: u64, reference: &str, alternate: &str) -> (String, bool) {
    let derived = from_sequences_detailed(
        "NC_TEST.1",
        lo,
        reference,
        alternate,
        &FromSequencesOptions::default(),
    )
    .unwrap_or_else(|e| panic!("from_sequences at {lo} over {reference} -> {alternate}: {e}"));
    (
        derived.variant.to_string(),
        derived.placement_bounded_by_window(),
    )
}

/// What a read covering `[lo, hi]` observes when the sample carries one fewer
/// copy of `unit`: the reference window with the first in-tract occurrence of
/// the unit removed.
///
/// `None` when the window holds no whole copy to remove — there is no variant to
/// observe, so there is nothing to ask about.
fn read_with_one_fewer_copy(lo: u64, hi: u64, tract_start: u64, unit: &str) -> Option<String> {
    let reference = at(lo, hi);
    let search_from = (tract_start.saturating_sub(lo)) as usize;
    let found = reference.get(search_from..)?.find(unit)? + search_from;
    let mut alternate = reference.clone();
    alternate.replace_range(found..found + unit.len(), "");
    Some(alternate)
}

/// The tract bounds in [`TRACTS`] are what [`SEQUENCE`] actually contains.
///
/// Asserted separately and first. The rest of this module is meaningless if a
/// tract's coordinates are off by even one base, and an off-by-one produces a
/// *plausible* wrong answer rather than an obvious one.
#[test]
fn the_tract_table_matches_the_sequence() {
    assert_eq!(SEQUENCE.len(), 57);
    for (label, first, last, unit) in TRACTS {
        let tract = at(*first, *last);
        assert_eq!(
            tract,
            unit.repeat(tract.len() / unit.len()),
            "{label}: {first}-{last} is {tract:?}, not a tandem repeat of {unit:?}"
        );
        // A tract that continued past `last` would make every "spanning" window
        // below a truncating one, and the sweep would silently invert.
        assert!(
            !SEQUENCE[*last as usize..].starts_with(unit),
            "{label}: the tract continues past {last}"
        );
    }
}

/// **Windows that span the whole tract all derive the same description.**
///
/// The sufficiency half. Both flanks are varied independently, including the
/// zero-flank window that is exactly the tract, so the claim is about
/// containment rather than about having room to spare.
#[test]
fn windows_that_span_the_tract_agree() {
    let mut checked = 0usize;
    for (label, first, last, unit) in TRACTS {
        let mut answers = BTreeSet::new();
        for lo in (first.saturating_sub(3).max(1))..=*first {
            for hi in *last..=(*last + 3).min(SEQUENCE.len() as u64) {
                let Some(alternate) = read_with_one_fewer_copy(lo, hi, *first, unit) else {
                    continue;
                };
                answers.insert(derive(lo, &at(lo, hi), &alternate).0);
                checked += 1;
            }
        }
        assert_eq!(
            answers.len(),
            1,
            "{label}: spanning windows reached {} descriptions: {answers:?}",
            answers.len()
        );
    }
    assert!(
        checked >= 4 * 16,
        "the sweep evaluated only {checked} windows, so its zero divergences say little"
    );
}

/// **A window that cuts the tract derives a different description, placed at
/// its own right edge.**
///
/// The necessity half, and the sharper claim: the truncated answers are not
/// merely "sometimes different", they track the window edge exactly. That is
/// what identifies the cause as the missing bases rather than as instability.
///
/// Every such answer is still a correct description of the same variant — see
/// the module docs — so this pins a rule 2 boundary, not a defect.
#[test]
fn a_window_that_cuts_the_tract_places_the_change_at_its_own_edge() {
    for (label, first, last, unit) in TRACTS {
        let lo = (first.saturating_sub(3)).max(1);
        let unit_len = unit.len() as u64;
        let mut truncated = BTreeSet::new();

        for hi in (first + unit_len - 1)..*last {
            let Some(alternate) = read_with_one_fewer_copy(lo, hi, *first, unit) else {
                continue;
            };
            let (rendered, edge) = derive(lo, &at(lo, hi), &alternate);
            // The change is placed at the last base the window holds.
            assert!(
                rendered.contains(&hi.to_string()),
                "{label}: window {lo}-{hi} cuts the tract but derived {rendered}, which does \
                 not name its own right edge {hi}"
            );
            assert!(
                edge,
                "{label}: window {lo}-{hi} cuts the tract, so its placement is bounded by the \
                 window and must be flagged; {rendered} was not"
            );
            truncated.insert(rendered);
        }

        let spanning = read_with_one_fewer_copy(lo, *last, *first, unit)
            .map(|alternate| derive(lo, &at(lo, *last), &alternate).0)
            .expect("the spanning window observes the variant");
        assert!(
            !truncated.is_empty(),
            "{label}: no truncating window was tested, so this proves nothing"
        );
        assert!(
            !truncated.contains(&spanning),
            "{label}: a truncating window reached the spanning answer {spanning}, so spanning \
             is not the boundary this module claims it is"
        );
    }
}

/// **The window-edge flag is conservative, in the direction that is safe.**
///
/// It is raised whenever the derivation *could* have been clipped, which
/// includes windows that span the tract exactly and are therefore already
/// correct. Distinguishing those two would require knowing what lies outside
/// the window — i.e. the reference — which this surface deliberately does not
/// read.
///
/// So the flag means **"a wider window could move this"**, never "this is
/// wrong". Pinned because the original doc comment claimed the stronger reading
/// and the `A4` row below refutes it.
#[test]
fn the_window_edge_flag_is_conservative_not_a_correctness_claim() {
    let (first, last, unit) = (12u64, 15u64, "A");

    // Flush against the tract's 3' end: flagged, and already the right answer.
    let flush = read_with_one_fewer_copy(9, last, first, unit).expect("observes the variant");
    let (flush_rendered, flush_edge) = derive(9, &at(9, last), &flush);
    assert!(
        flush_edge,
        "a window ending at the tract's last base is flagged"
    );

    // Room to spare: not flagged, and the same answer.
    let roomy = read_with_one_fewer_copy(9, last + 3, first, unit).expect("observes the variant");
    let (roomy_rendered, roomy_edge) = derive(9, &at(9, last + 3), &roomy);
    assert!(
        !roomy_edge,
        "a window with 3' flank beyond the tract is not flagged"
    );

    assert_eq!(
        flush_rendered, roomy_rendered,
        "the flagged window and the unflagged one must agree — which is precisely why the \
         flag cannot be read as 'this answer is wrong'"
    );
}
