//! Per-spelling verdicts for the externally reported partition family
//! (#1419, #1420, #1421).
//!
//! # What these rows are, and why they were not already runnable
//!
//! Three issues report the same defect from a downstream pipeline: one variant,
//! written two ways, normalizes to two different strings. Each issue body states
//! four things about each of its rows, in prose:
//!
//! 1. what ferro prints for **each** of the two spellings,
//! 2. which of the two the HGVS recommendations make canonical, and on which
//!    line,
//! 3. that `EquivalenceChecker` still calls the two one variant
//!    (`SequenceMatch`) — which is what makes the divergence a *representation*
//!    problem rather than a correctness one, and what leaves consumers a
//!    working fallback,
//! 4. the reference bases the argument in (2) rests on.
//!
//! Only (1)'s *pairwise disagreement* was guarded. `reported_confluence_pairs`
//! pins that the two spellings must converge, and deliberately declines to pin
//! *which* string they converge on — pinning a winner there would freeze an open
//! product decision (#1235) by accident. That is the right call for a ratchet on
//! a future fix, and it leaves the four claims above untested: today's answers,
//! the spec argument for them, the equivalence fallback, and the reference.
//!
//! This module tests those. It is the *present-tense* half — what ferro does
//! now and why that is or is not what the recommendations ask for — where
//! `reported_confluence_pairs` is the *future-tense* half, ratcheting the fix.
//! Neither subsumes the other and they assert disjoint things.
//!
//! # Convention: characterization, not `#[ignore]`
//!
//! Nine of the eighteen rows below print something the cited recommendation
//! argues against. They are pinned **as they behave today**, with the wanted
//! answer written into the row, rather than being expressed as an `#[ignore]`d
//! test of the wanted answer.
//!
//! That is deliberate and it is the whole point. An `#[ignore]`d test never
//! runs, so it cannot notice the day the answer moves — and for this family the
//! *move* is the event that matters. A normalized string is the key a downstream
//! consumer stores read counts against, so any change to one of these rows is a
//! migration for somebody, whether it is a fix or a regression. Pinned, it
//! surfaces in a diff and has to be argued for in review. Ignored, it happens
//! silently and is discovered downstream.
//!
//! [`OPEN_GAPS`] counts the disagreeing rows so the count cannot drift away from
//! the table: closing one means editing a row *and* the census, deliberately.
//!
//! # The reference
//!
//! All three issues ran on one synthetic 125 nt contig, reproduced here verbatim
//! so these rows measure the bases they were reported against.
//! [`reference_bases_are_what_the_reports_state`] asserts that premise before
//! anything else; without it every row below could be measuring a different
//! locus and still look green.

use crate::common::cis_apply_oracle::{normalize, provider};
use ferro_hgvs::equivalence::{EquivalenceChecker, EquivalenceLevel};
use ferro_hgvs::parse_hgvs;

/// The 125 nt contig every row runs against. Position 1 is index 0.
///
/// Synthetic — it names no real sequence. Reproduced from the reports so these
/// rows measure the same bases the reports did.
const TEMPLATE: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

/// Whether ferro's answer for one spelling is the form its issue argues for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    /// Ferro prints the form the cited recommendation makes canonical.
    Canonical,
    /// Ferro prints something else, and the row records what was wanted.
    Gap,
}

/// One spelling of one reported variant.
struct Row {
    /// `<issue>-<row>`, matching the label used in `reported_confluence_pairs`
    /// so the two modules' rows can be lined up by eye.
    label: &'static str,
    /// The spelling as authored in the issue.
    input: &'static str,
    /// What ferro prints for it today. Measured, never transcribed.
    output: &'static str,
    verdict: Verdict,
    /// Why that is or is not the canonical form, citing the line of the pinned
    /// spec checkout the issue cites.
    argument: &'static str,
}

/// The eighteen spellings — two per reported pair.
///
/// Both spellings of a pair are listed rather than only the losing one, because
/// the defect is not "this string is wrong": both strings are well-formed and
/// both denote the same sequence. The defect is that the canonical form is
/// reachable from one spelling and not from the other, and you cannot state that
/// without both halves. See
/// [`each_pair_reaches_its_canonical_form_from_exactly_one_spelling`].
const REPORTED_ROWS: &[Row] = &[
    // -- #1419 -- a net deletion spelled as two deletions, versus as one span.
    //
    // Neither spelling is a delins across an unchanged base, so no local rule
    // is violated by either; the issue's own framing is that this is pure
    // partition non-uniqueness. The canonical argument is prioritisation.
    Row {
        label: "1419-r1/cis",
        input: "TEMPLATE:g.[19_23del;27_33del]",
        output: "TEMPLATE:g.[19_23del;27_33del]",
        verdict: Verdict::Gap,
        argument: "Wanted `[19_30del;33T>G]`. general.md:56 prioritises \
                   (1) substitution over (2) deletion, so the change at 33 \
                   should be exposed as a substitution rather than absorbed \
                   into a deletion. The cis spelling hides it; the equivalent \
                   span reaches it (see `1419-r1/span`), so the partition is \
                   reachable and simply is not reached from here.",
    },
    Row {
        label: "1419-r1/span",
        input: "TEMPLATE:g.19_33delinsCGG",
        output: "TEMPLATE:g.[19_30del;33T>G]",
        verdict: Verdict::Canonical,
        argument: "The form #1419 names canonical: re-derived from the \
                   resulting sequence, exposing the substitution at 33 per \
                   general.md:56, with the retained bases placed as 3' as \
                   possible per general.md:41.",
    },
    Row {
        label: "1419-r2/cis",
        input: "TEMPLATE:g.[19_22del;26_36del]",
        output: "TEMPLATE:g.[19_22del;26_36del]",
        verdict: Verdict::Gap,
        argument: "Wanted `[19_33del;36A>G]`, for the same general.md:56 \
                   reason as `1419-r1/cis`, one locus over.",
    },
    Row {
        label: "1419-r2/span",
        input: "TEMPLATE:g.19_36delinsGCG",
        output: "TEMPLATE:g.[19_33del;36A>G]",
        verdict: Verdict::Canonical,
        argument: "The form #1419 names canonical for row 2.",
    },
    Row {
        label: "1419-r3/cis",
        input: "TEMPLATE:g.[19_24del;28_33del]",
        output: "TEMPLATE:g.[19_24del;28_33del]",
        verdict: Verdict::Gap,
        argument: "Wanted `[19T>G;22_33del]`. Same general.md:56 reason, but \
                   note the exposed substitution lands at the 5' end here \
                   rather than the 3' end — which is why #1419 argues no \
                   per-end rule fixes the family.",
    },
    Row {
        label: "1419-r3/span",
        input: "TEMPLATE:g.19_33delinsGGA",
        output: "TEMPLATE:g.[19T>G;22_33del]",
        verdict: Verdict::Canonical,
        argument: "The form #1419 names canonical for row 3.",
    },
    // -- #1420 -- one row per member type. v2/v3 must reduce, v4 must coalesce:
    // the corrections point in opposite directions, which is the issue's
    // argument that no per-member-type rule reconciles them.
    Row {
        label: "1420-v2/cis",
        input: "TEMPLATE:g.[37dup;41del]",
        output: "TEMPLATE:g.[37dup;41del]",
        verdict: Verdict::Gap,
        argument: "Wanted `[38T>A;40_41delinsTG]`. Reference 38-41 is `TTGC` \
                   and the result is `ATTG`, so 38 and 40-41 change and 39 \
                   does not (asserted in `reported_spans_change_the_columns_\
                   the_reports_state`). general.md:56 prioritises \
                   (1) substitution over (4) duplication, so the change at 38 \
                   must not be spelled as a `dup`.",
    },
    Row {
        label: "1420-v2/span",
        input: "TEMPLATE:g.38_41delinsATTG",
        output: "TEMPLATE:g.[38T>A;40_41delinsTG]",
        verdict: Verdict::Canonical,
        argument: "The form #1420 names canonical for v2: substitution exposed \
                   at 38, and the span separated at the unchanged 39 per \
                   delins.md:17.",
    },
    Row {
        label: "1420-v3/cis",
        input: "TEMPLATE:g.[36_37insC;40del]",
        output: "TEMPLATE:g.[36_37insC;40del]",
        verdict: Verdict::Gap,
        argument: "Wanted `[37_38delinsCA;40G>T]`. general.md:56 prioritises \
                   (1) substitution over (5) insertion, so the substitution at \
                   40 must be exposed rather than left inside an `ins`+`del` \
                   pair.",
    },
    Row {
        label: "1420-v3/span",
        input: "TEMPLATE:g.37_40delinsCATT",
        output: "TEMPLATE:g.[37_38delinsCA;40G>T]",
        verdict: Verdict::Canonical,
        argument: "The form #1420 names canonical for v3: reference 37-40 is \
                   `ATTG` against `CATT`, so 37-38 change consecutively \
                   (delins.md:16), 39 does not, and 40 is a substitution \
                   separated from them (delins.md:17).",
    },
    Row {
        label: "1420-v4/cis",
        input: "TEMPLATE:g.[21delinsGC;24del]",
        output: "TEMPLATE:g.[21delinsGC;24del]",
        verdict: Verdict::Gap,
        argument: "Wanted `21_24delinsGCTG` — the opposite direction to v2/v3. \
                   Reference 21-24 is `ATGC` against `GCTG`, so all four \
                   positions change, and delins.md:16 says changes involving \
                   two or more CONSECUTIVE nucleotides are one delins. The \
                   split spelling asserts 22-23 are unchanged; they are not.",
    },
    Row {
        label: "1420-v4/span",
        input: "TEMPLATE:g.21_24delinsGCTG",
        output: "TEMPLATE:g.21_24delinsGCTG",
        verdict: Verdict::Canonical,
        argument: "The form #1420 names canonical for v4: one delins over four \
                   consecutive changed positions, per delins.md:16.",
    },
    // -- #1421 -- net insertions. Here the SPLIT spelling is the canonical one,
    // the reverse of #1419/#1420, so a fix that simply preferred whichever form
    // has fewer members would get this family backwards.
    Row {
        label: "1421-n1/split",
        input: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        output: "TEMPLATE:g.[29C>A;32_33delinsACATACTG]",
        verdict: Verdict::Canonical,
        argument: "The form #1421 names canonical: the substitution at 29 and \
                   the delins at 32-33 are separated by two unchanged \
                   nucleotides, and delins.md:17 says variants separated by one \
                   or more nucleotides are described individually and not as a \
                   delins.",
    },
    Row {
        label: "1421-n1/span",
        input: "TEMPLATE:g.29_33delinsAACACATACTG",
        output: "TEMPLATE:g.29_33delinsAACACATACTG",
        verdict: Verdict::Gap,
        argument: "Wanted `[29C>A;32_33delinsACATACTG]`. This span merges \
                   across the unchanged 30-31, which delins.md:17 forbids. The \
                   delins.md:18 exception (two variants separated by ONE \
                   nucleotide together affecting one amino acid) does not \
                   reach it: the separation is two nucleotides, and a `g.` \
                   description has no amino acid. See \
                   `the_1421_spans_separate_by_two_nucleotides_not_one`.",
    },
    Row {
        label: "1421-n2/split",
        input: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        output: "TEMPLATE:g.[32G>T;35_36delinsGAATCGAC]",
        verdict: Verdict::Canonical,
        argument: "The form #1421 names canonical for row 2; unchanged \
                   interior 33-34.",
    },
    Row {
        label: "1421-n2/span",
        input: "TEMPLATE:g.32_36delinsTTGGAATCGAC",
        output: "TEMPLATE:g.32_36delinsTTGGAATCGAC",
        verdict: Verdict::Gap,
        argument: "Wanted `[32G>T;35_36delinsGAATCGAC]`, per delins.md:17.",
    },
    Row {
        label: "1421-n3/split",
        input: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        output: "TEMPLATE:g.[34G>T;37_38delinsCCTTTACG]",
        verdict: Verdict::Canonical,
        argument: "The form #1421 names canonical for row 3; unchanged \
                   interior 35-36.",
    },
    Row {
        label: "1421-n3/span",
        input: "TEMPLATE:g.34_38delinsTCACCTTTACG",
        output: "TEMPLATE:g.34_38delinsTCACCTTTACG",
        verdict: Verdict::Gap,
        argument: "Wanted `[34G>T;37_38delinsCCTTTACG]`, per delins.md:17.",
    },
];

/// How many of [`REPORTED_ROWS`] print something their issue argues against.
///
/// Nine of eighteen — exactly one per reported pair, which is the family's
/// shape rather than a coincidence; see
/// [`each_pair_reaches_its_canonical_form_from_exactly_one_spelling`].
///
/// **This may only go down.** Up means a change has moved a row *away* from the
/// cited recommendation.
const OPEN_GAPS: usize = 9;

/// The `<issue>-<row>` half of an `<issue>-<row>/<spelling>` label.
fn pair_of(label: &'static str) -> &'static str {
    label
        .split('/')
        .next()
        .expect("labels are `<issue>-<row>/<spelling>`")
}

/// [`REPORTED_ROWS`] read two at a time, checked to really be nine pairs.
///
/// `chunks(2)` on its own trusts the table's row order, and that trust is not
/// free. Reorder [`REPORTED_ROWS`] so a spelling sits beside a *different*
/// pair's spelling and
/// [`each_pair_reaches_its_canonical_form_from_exactly_one_spelling`] still
/// passes — measured, not assumed: every assertion in it compares a row against
/// its own pinned output, so a mispaired neighbour changes nothing it looks at.
/// Only [`every_reported_pair_is_still_one_variant_by_equivalence`] notices, and
/// it reports the mispairing as `NotEquivalent` — that is, as a regression in
/// the equivalence fallback, which is the wrong diagnosis and sends a reader to
/// the wrong file.
///
/// Splitting the label costs nothing and makes the table error itself the
/// failure.
fn reported_pairs() -> impl Iterator<Item = (&'static Row, &'static Row)> {
    REPORTED_ROWS.chunks(2).map(|chunk| {
        let [a, b] = chunk else {
            panic!("REPORTED_ROWS must hold both spellings of every pair")
        };
        assert_eq!(
            pair_of(a.label),
            pair_of(b.label),
            "{} and {} are not two spellings of one reported pair; \
             REPORTED_ROWS has fallen out of pair order",
            a.label,
            b.label
        );
        (a, b)
    })
}

/// Normalize against [`TEMPLATE`] in the shipped default 3' direction.
///
/// 3' only, and on purpose: this module pins the strings the reports observed,
/// and the reports ran through the default. `reported_confluence_pairs` is the
/// module that sweeps both directions, because *convergence* has to hold in
/// both for a fix to count, whereas the canonical string itself is
/// direction-dependent by construction (general.md:41).
fn normalized(input: &str) -> String {
    normalize(TEMPLATE, input)
}

/// The premise: the bases the reports' arguments rest on are the bases here.
///
/// Asserted first and separately. Every canonical argument in
/// [`REPORTED_ROWS`] is a claim about which reference positions change, so if
/// the contig above were transcribed wrong, the rows would still pass — they
/// only compare ferro against itself — while every `argument` field silently
/// became fiction.
#[test]
fn reference_bases_are_what_the_reports_state() {
    assert_eq!(TEMPLATE.len(), 125, "contig length");
    for (first, last, expected, source) in [
        (21usize, 24usize, "ATGC", "#1420 v4: reference 21-24"),
        (29, 33, "CACGT", "#1421: positions 29-33"),
        (37, 40, "ATTG", "#1420 v3: reference 37-40"),
        (38, 41, "TTGC", "#1420 v2: reference 38-41"),
    ] {
        assert_eq!(
            &TEMPLATE[first - 1..last],
            expected,
            "{source} should be {expected}"
        );
    }
}

/// Every spelling still prints what the reports recorded for it.
///
/// The eighteen pins. A failure here is a representation change on a shape a
/// downstream consumer keys read counts on, so it is meant to be argued for in
/// review rather than re-blessed: the message names the verdict and the wanted
/// form so a reviewer can tell a fix from a regression without reopening the
/// issue.
#[test]
fn every_reported_spelling_still_normalizes_as_the_reports_recorded() {
    for row in REPORTED_ROWS {
        assert_eq!(
            normalized(row.input),
            row.output,
            "{} moved.\n  input:   {}\n  pinned:  {}\n  verdict: {:?}\n  {}",
            row.label,
            row.input,
            row.output,
            row.verdict,
            row.argument
        );
    }
}

/// Each reported pair reaches its canonical form from exactly one spelling.
///
/// This is the defect in one assertion, and it is sharper than "the two
/// spellings disagree". Both spellings denote one variant; the canonical form
/// is a property of the variant, so it should be reachable from either. It is
/// reachable from exactly one, which makes reachability a function of how the
/// variant was *written* — the spelling-dependence root-caused on #1419 to a
/// weight bound whose threshold is the input's own spelling.
///
/// While this test passes, the family is unfixed. It is written to fail the day
/// any pair reaches its canonical form from both spellings, so that progress is
/// as loud as regression.
#[test]
fn each_pair_reaches_its_canonical_form_from_exactly_one_spelling() {
    for (a, b) in reported_pairs() {
        assert_ne!(
            a.verdict, b.verdict,
            "{} and {} have the same verdict; a pair is one canonical spelling \
             and one gap",
            a.label, b.label
        );

        let (canonical, gap) = match a.verdict {
            Verdict::Canonical => (a, b),
            Verdict::Gap => (b, a),
        };

        assert_eq!(
            normalized(canonical.input),
            canonical.output,
            "{}: the canonical spelling stopped reaching the canonical form",
            canonical.label
        );
        assert_ne!(
            normalized(gap.input),
            canonical.output,
            "{} now reaches `{}` too, so the pair has converged on the form \
             its issue argues for. That is the fix: flip this row's verdict to \
             Canonical, lower OPEN_GAPS, and say in the PR which stored \
             spelling moved.",
            gap.label,
            canonical.output
        );
    }
}

/// Every gap row is returned exactly as it was authored.
///
/// Deliberately *not* the same claim as `reported_confluence_pairs`'
/// `every_reported_output_is_a_fixed_point`, which normalizes twice and checks
/// the second pass changes nothing. That would also be satisfied by a row whose
/// first pass moved it somewhere new and then settled.
///
/// This asserts the stronger thing the reports actually observed: the input is
/// *already* a fixed point — ferro does not touch it. That distinction is what
/// rules out the benign reading of this family. A non-canonical spelling that
/// moved, even to some third string, would be an under-applied pass. One that is
/// returned untouched is a second canonical form, and a second canonical form is
/// what splits a consumer's counts across two keys.
#[test]
fn every_gap_row_is_returned_exactly_as_authored() {
    for row in REPORTED_ROWS {
        if row.verdict != Verdict::Gap {
            continue;
        }
        let output = normalized(row.input);
        assert_eq!(
            output, row.input,
            "{}: the reported spelling is no longer retained verbatim. If it \
             moved to the canonical form the family is fixed; if it moved \
             somewhere else this is a new representation nobody asked for.",
            row.label
        );
    }
}

/// `EquivalenceChecker` still calls each reported pair one variant.
///
/// All three issues state this, and it is the load-bearing mitigation rather
/// than a detail: it is what makes the divergence a representation problem
/// instead of a correctness one, and it is the fallback the issues tell
/// consumers to compare on while the family is open. If it regressed, the
/// reports' "use `EquivalenceChecker` instead" advice would quietly stop
/// working — and nothing else in the suite covers these shapes.
///
/// `SequenceMatch` exactly, not merely `is_equivalent()`: the two spellings
/// normalize to different strings here, so anything stronger would mean the
/// pair had converged and the rows above are stale.
#[test]
fn every_reported_pair_is_still_one_variant_by_equivalence() {
    let checker = EquivalenceChecker::new(provider(TEMPLATE));
    for (a, b) in reported_pairs() {
        let left = parse_hgvs(a.input).unwrap_or_else(|e| panic!("{}: {e}", a.label));
        let right = parse_hgvs(b.input).unwrap_or_else(|e| panic!("{}: {e}", b.label));
        let result = checker
            .check(&left, &right)
            .unwrap_or_else(|e| panic!("{} vs {}: {e}", a.label, b.label));
        assert_eq!(
            result.level,
            EquivalenceLevel::SequenceMatch,
            "{} vs {}: expected SequenceMatch, got {:?}. If this is now \
             NormalizedMatch the pair has converged and the pinned rows are \
             stale; anything else is a regression in the fallback the reports \
             direct consumers to.",
            a.label,
            b.label,
            result.level
        );
    }
}

/// `(label, first changed position, reference bases, replacement)` for the
/// three #1420 rows whose replacement is the same length as the span it
/// replaces.
///
/// Equal length is what makes the changed/unchanged column pattern well defined
/// without an alignment, which is why only #1420 appears here: #1419's rows are
/// net deletions and #1421's are net insertions, so their columns do not line
/// up one-to-one and the pattern is a choice rather than a fact.
const EQUAL_LENGTH_SPANS: &[(&str, usize, &str, &str)] = &[
    ("1420-v2", 38, "TTGC", "ATTG"),
    ("1420-v3", 37, "ATTG", "CATT"),
    ("1420-v4", 21, "ATGC", "GCTG"),
];

/// The #1420 arguments are checked against the reference, not taken on trust.
///
/// Each row's canonical form follows from *which columns change*: v2 and v3
/// have an unchanged interior column, so delins.md:17 separates them; v4 has
/// none, so delins.md:16 coalesces it. That distinction is the reason v4's
/// correction points the opposite way to v2's and v3's, and it is the one part
/// of #1420's reasoning that is a fact about the sequence rather than a reading
/// of the recommendations — so it is computed here rather than asserted in a
/// comment.
#[test]
fn reported_spans_change_the_columns_the_reports_state() {
    // `X` = changed, `.` = unchanged, one character per position in the span.
    let expected: &[(&str, &str)] = &[
        ("1420-v2", "X.XX"),
        ("1420-v3", "XX.X"),
        ("1420-v4", "XXXX"),
    ];
    // `zip` stops at the shorter side, so a span added to one table and not the
    // other would be skipped rather than reported.
    assert_eq!(
        EQUAL_LENGTH_SPANS.len(),
        expected.len(),
        "every equal-length span needs a column pattern"
    );

    for ((label, first, reference, replacement), (echo, pattern)) in
        EQUAL_LENGTH_SPANS.iter().zip(expected)
    {
        assert_eq!(label, echo, "the two tables fell out of order");
        assert_eq!(
            reference.len(),
            replacement.len(),
            "{label}: this table is only meaningful for equal-length spans"
        );
        assert_eq!(
            &TEMPLATE[first - 1..first - 1 + reference.len()],
            *reference,
            "{label}: reference bases at {first} are not what the report states"
        );

        let observed: String = reference
            .bytes()
            .zip(replacement.bytes())
            .map(|(r, a)| if r == a { '.' } else { 'X' })
            .collect();
        assert_eq!(
            &observed, pattern,
            "{label}: {reference} -> {replacement} changes columns {observed}, \
             not {pattern} as the report states"
        );
    }

    // The consequence, stated so it is not left implicit: two rows have an
    // unchanged interior column and one does not, which is exactly why #1420's
    // three corrections cannot share a per-member-type rule.
    let separable = expected
        .iter()
        .filter(|(_, pattern)| pattern.trim_matches('X').contains('.'))
        .count();
    assert_eq!(
        separable, 2,
        "expected v2 and v3 to be separable under delins.md:17 and v4 to be \
         one run under delins.md:16"
    );
}

/// `(label, last position of the leading substitution, first position of the
/// trailing delins)` for #1421's three split spellings.
const ISSUE_1421_MEMBER_BOUNDS: &[(&str, usize, usize)] = &[
    ("1421-n1", 29, 32),
    ("1421-n2", 32, 35),
    ("1421-n3", 34, 37),
];

/// #1421's rows fall under delins.md:17, not under the delins.md:18 exception.
///
/// Worth its own test because the exception is the only reading under which the
/// retained spans would be correct, and it turns on a number: it covers two
/// variants separated by **one** nucleotide that together affect one amino
/// acid. Every #1421 row is separated by two, and is a `g.` description with no
/// amino acid, so it fails the exception twice over.
///
/// The interior positions are unchanged by construction — they lie strictly
/// between the split spelling's two members and are named by neither, so a cis
/// allele leaves them alone. What needs checking is the width.
///
/// The two member positions are read back out of the corresponding row in
/// [`REPORTED_ROWS`] before the arithmetic runs. Without that the test would be
/// `32 - 29 - 1 == 2` over its own constants — true of the table rather than of
/// the reported variants, and still green if somebody moved a row's coordinates.
#[test]
fn the_1421_spans_separate_by_two_nucleotides_not_one() {
    for (label, substitution, delins_start) in ISSUE_1421_MEMBER_BOUNDS {
        let split_label = format!("{label}/split");
        let row = REPORTED_ROWS
            .iter()
            .find(|row| row.label == split_label)
            .unwrap_or_else(|| panic!("no row labelled `{split_label}`"));

        // The split spelling is `<acc>:g.[<sub-position><ref>><alt>;<delins-start>_...]`,
        // so these two anchors pin the bounds below to the row's own coordinates.
        assert!(
            row.input.contains(&format!("g.[{substitution}")),
            "{label}: the bounds table says the substitution is at \
             {substitution}, but the row's first member is `{}`",
            row.input
        );
        assert!(
            row.input.contains(&format!(";{delins_start}_")),
            "{label}: the bounds table says the delins starts at \
             {delins_start}, but the row's second member is `{}`",
            row.input
        );

        let separation = delins_start - substitution - 1;
        assert_eq!(
            separation, 2,
            "{label}: members at {substitution} and {delins_start} are \
             separated by {separation} nt; the delins.md:18 exception covers \
             one, so a separation of two puts this row squarely under \
             delins.md:17"
        );
    }
}

/// The two reported-family modules describe the same nine pairs, on the same
/// bases.
///
/// Each module holds its own copy of the contig and of all eighteen spellings,
/// and the [`Row::label`] field promises the rows line up with
/// `reported_confluence_pairs`. Nothing enforced that promise, so editing one
/// table would leave the two modules silently measuring different variants while
/// both stayed green — this one pinning verdicts for one set of spellings and
/// `reported_confluence_pairs` ratcheting convergence for another. Since a fix
/// is expected to move rows in *both* modules together, the two falling out of
/// step is a live hazard rather than a theoretical one.
#[test]
fn the_two_reported_modules_describe_the_same_pairs() {
    assert_eq!(
        TEMPLATE,
        crate::reported_confluence_pairs::TEMPLATE,
        "the two modules' contigs have diverged, so their rows are no longer \
         about the same bases"
    );

    let pairs = crate::reported_confluence_pairs::REPORTED_PAIRS;
    assert_eq!(
        pairs.len() * 2,
        REPORTED_ROWS.len(),
        "one module has gained or lost a reported pair"
    );

    for ((label, a, b), (row_a, row_b)) in pairs.iter().zip(reported_pairs()) {
        assert_eq!(
            *label,
            pair_of(row_a.label),
            "the two modules list the reported pairs in different orders"
        );
        assert_eq!(
            [*a, *b],
            [row_a.input, row_b.input],
            "{label}: the two modules disagree about how this pair is spelled"
        );
    }
}

/// The census, so the count and the table cannot drift apart.
#[test]
fn the_open_gap_census_holds() {
    let gaps = REPORTED_ROWS
        .iter()
        .filter(|row| row.verdict == Verdict::Gap)
        .count();
    assert_eq!(
        gaps, OPEN_GAPS,
        "the number of rows disagreeing with their cited recommendation moved \
         (now {gaps}). Down is progress: lower OPEN_GAPS and name the moved \
         representation in the PR. Up means a change pushed a row away from \
         the spec."
    );
    assert_eq!(
        REPORTED_ROWS.len(),
        18,
        "REPORTED_ROWS must stay two spellings per reported pair"
    );
}
