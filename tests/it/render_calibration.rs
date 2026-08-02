//! Rendering rules pinned on single-member inputs.
//!
//! Partitioning decides *where* members are; rendering decides how each is
//! spelled — the 3' rule (`general.md:41`) and prioritisation (`general.md:56`:
//! substitution > deletion > inversion > duplication > insertion). A change can
//! partition correctly and still spell wrongly, and nothing else in the suite
//! tests that axis systematically.
//!
//! Single-member only, deliberately: with one member there is no sibling
//! interaction, so a failure here is unambiguously a rendering fault.
//!
//! # What these rows do *not* cover
//!
//! Stated explicitly so the file is not read as wider than it is. All five
//! prioritisation outcomes now have a row: substitution, deletion, inversion,
//! duplication and insertion each appear as the winning spelling. What is *not*
//! covered is the full pairwise ladder — each row shows one outcome winning,
//! not that it beats every lower-ranked alternative in turn.
//!
//! **Repeat/STR notation is exercised, on both sides of the boundary where
//! ferro starts emitting it.** One unit added or removed stays a `dup`/`del`:
//! over core `CAGCAGCAGG`, a fourth `CAG` copy (whether spelled
//! `g.265_266insCAG` or `g.263_265dup`) renders as `g.263_265dup`, and dropping
//! a copy (`g.263_265del`) renders as `g.263_265del` — never `CAG[4]`/`CAG[2]`,
//! which is consistent with prioritisation (duplication outranks insertion).
//! Two units is where the whole run is described as one repeat instead, and
//! `repeat_notation_describes_the_whole_run` pins that.
//!
//! Which side of that boundary a change falls on is a rendering decision, and
//! it is the one place this axis describes a whole ambiguous run as a single
//! member. It does so only as `dup` or `NNN[k]` — never as a `delins` spanning
//! the run — so those rows are the guard against a minimisation pass that
//! narrows more than it should.

use crate::common::cis_apply_oracle::normalize;
use crate::common::synthetic::padded;

/// `(rule, core, input, expected)`.
const CASES: &[(&str, &str, &str, &str)] = &[
    // 3' rule: a deletion inside a homopolymer takes the most 3' position.
    (
        "3' rule, deletion in a run",
        "CAAAAG",
        "TEMPLATE:g.258del",
        "TEMPLATE:g.261del",
    ),
    // 3' rule: an insertion's junction travels to the tract's 3' end.
    (
        "3' rule, insertion in a run",
        "CAAAAG",
        "TEMPLATE:g.257_258insA",
        "TEMPLATE:g.261dup",
    ),
    // Prioritisation: an insertion duplicating the preceding base is a dup.
    (
        "ins -> dup",
        "CTAGG",
        "TEMPLATE:g.258_259insT",
        "TEMPLATE:g.258dup",
    ),
    // Prioritisation: a one-base delins restating a different base is a sub.
    (
        "delins -> sub",
        "CTAGG",
        "TEMPLATE:g.258delinsG",
        "TEMPLATE:g.258T>G",
    ),
    // A delins whose replacement is the reverse complement is an inversion
    // (`DNA/inversion.md:5`: "more than one nucleotide replacing the original
    // sequence is the reverse complement of the original sequence").
    //
    // The task brief's original row here used core `CTAGCTAG` with
    // `g.258_261delinsCTAG`, expecting `g.258_261inv`. That expectation was
    // wrong on the arithmetic, not on the spec: over this file's padding
    // (`padded`, `PAD_OFFSET` = 256), position 258 is the core's *second*
    // base, so `g.258_261` deletes core bases `TAGC`, not `CTAG`.
    // revcomp(`TAGC`) = `GCTA`, which is not the stated replacement `CTAG` —
    // so that delins never described an inversion in the first place, and
    // ferro leaving it spelled as `delins` was the correct call. (`CTAG`
    // itself is a palindrome — revcomp(`CTAG`) = `CTAG` — so no shift of the
    // same span could have produced a genuine inversion either; the row
    // needed a non-self-complementary 4-mer.) Replaced with `CATCGG` /
    // `delinsCGAT`, where the deleted span `ATCG` and the replacement `CGAT`
    // really are reverse complements of one another, to exercise the same
    // rule with a case that is actually an inversion.
    (
        "delins -> inv",
        "CATCGG",
        "TEMPLATE:g.258_261delinsCGAT",
        "TEMPLATE:g.258_261inv",
    ),
    // Prioritisation: deletion's own slot. Over core `CTAGG` (257..261 =
    // C,T,A,G,G), `g.258_259` is `TA` and the replacement is `A`, so the net
    // edit removes one base and must render as a deletion — not as a delins
    // restating what it kept. Deleting 258 (`T`) yields `CAGG`, which is the
    // sequence `delinsA` denotes.
    (
        "delins -> del",
        "CTAGG",
        "TEMPLATE:g.258_259delinsA",
        "TEMPLATE:g.258del",
    ),
    // Prioritisation: insertion's own slot — the lowest rank, so it is what
    // survives when nothing higher applies. `CC` neither duplicates the
    // flanking bases nor inverts them, so it must stay spelled as an
    // insertion rather than being re-spelled as a dup.
    (
        "ins stays ins",
        "CTAGG",
        "TEMPLATE:g.258_259insCC",
        "TEMPLATE:g.258_259insCC",
    ),
];

#[test]
fn rendering_rules_are_pinned() {
    for (rule, core, input, expected) in CASES {
        let seq = padded(core);
        let got = normalize(&seq, input);
        assert_eq!(
            &got, expected,
            "{rule}: `{input}` over core `{core}` rendered as `{got}`, expected `{expected}`"
        );
    }
}

/// Where a change inside a homopolymer is described as the whole run, and where
/// it is not.
///
/// One unit is localised (`g.263del`, `g.261dup` — the same answers rows 1 and 2
/// of `CASES` pin); two units is spelled as the whole tract in repeat notation,
/// which is *wider* than the change itself. Both directions grow: a 2-base
/// deletion out of a 7-A run and a 2-base insertion into a 5-A run.
///
/// The grown forms are the ones to watch. They are the only members on this axis
/// that describe more reference than they change, and they are legitimate
/// because `NNN[k]` is an edit type of its own — a `delins` over the same span
/// would not be. Any pass that narrows members to their differences must leave
/// these four rows exactly as they are.
#[test]
fn repeat_notation_describes_the_whole_run() {
    for (core, input, expected) in [
        // Core `AAAAAA` at 257..262, and the pad's trailing `A` at 263 extends
        // it, so the run is 257..263 — seven A.
        ("AAAAAA", "TEMPLATE:g.260del", "TEMPLATE:g.263del"),
        ("AAAAAA", "TEMPLATE:g.260_261del", "TEMPLATE:g.257_263A[5]"),
        // Core `CAAAAG`: a four-A run at 258..261.
        ("CAAAAG", "TEMPLATE:g.260_261insA", "TEMPLATE:g.261dup"),
        (
            "CAAAAG",
            "TEMPLATE:g.260_261insAA",
            "TEMPLATE:g.258_261A[6]",
        ),
    ] {
        let seq = padded(core);
        let got = normalize(&seq, input);
        assert_eq!(
            &got, expected,
            "`{input}` over core `{core}` rendered as `{got}`, expected `{expected}`"
        );
    }
}

/// Rendering must not depend on which equivalent spelling was supplied: two
/// spellings of one single-member variant must render identically.
#[test]
fn rendering_does_not_depend_on_the_input_spelling() {
    for (core, a, b, expected) in [
        (
            "CTAGG",
            "TEMPLATE:g.258_259insT",
            "TEMPLATE:g.258dup",
            "TEMPLATE:g.258dup",
        ),
        (
            "CAAAAG",
            "TEMPLATE:g.258del",
            "TEMPLATE:g.261del",
            "TEMPLATE:g.261del",
        ),
        // A tandem-repeat expansion, reached from both spellings: inserting a
        // fourth `CAG` after the tract, versus duplicating the third copy.
        // Both denote `CAGCAGCAGCAGG` and must render alike.
        (
            "CAGCAGCAGG",
            "TEMPLATE:g.265_266insCAG",
            "TEMPLATE:g.263_265dup",
            "TEMPLATE:g.263_265dup",
        ),
    ] {
        let seq = padded(core);
        let (na, nb) = (normalize(&seq, a), normalize(&seq, b));
        assert_eq!(na, nb, "`{a}` -> `{na}` but `{b}` -> `{nb}`");
        // Agreement alone is not enough: two spellings converging on the same
        // *wrong* form satisfies it. Rows 1 and 2 are anchored elsewhere (in
        // `CASES` and in `repeat_notation_describes_the_whole_run`); the third was
        // anchored nowhere, so its expected answer — stated in this module's own
        // doc — is pinned here.
        assert_eq!(
            na, expected,
            "both spellings agree, but on `{na}` rather than the pinned `{expected}`"
        );
    }
}
