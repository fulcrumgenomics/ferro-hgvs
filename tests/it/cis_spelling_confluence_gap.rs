//! Two spellings of one variant that normalize to two different strings.
//!
//! #1235's criterion 1 requires every encoding of a variant to reach one
//! normalized string. These eight pairs do not, and each pair's *split*
//! spelling is an expectation this repository blessed and shipped.
//!
//! The pairs were found by deriving each variant's minimal-alignment partition
//! from the resulting sequence, rendering the derived single-block form, and
//! normalizing both spellings. Every pair below was verified to denote the same
//! sequence by an applier independent of the normalizer.
//!
//! Each row asserts the two spellings still DIVERGE. This test failing means a
//! pair converged, and that pair's blessed expectation should be re-blessed to
//! the converged form.
//!
//! Both tables are measured in both directions. Every row above was blessed
//! against the 3' rule, but confluence is a property of the normalizer rather
//! than of one shuffle direction, and `--direction 5prime` is a supported public
//! option — so [`DIVERGENT_UNDER_FIVE_PRIME`] records the [`CONVERGED`] rows that
//! converge under 3' and diverge under 5', and
//! [`the_eight_spelling_pairs_still_diverge_under_five_prime`] pins that the
//! [`DIVERGENT`] eight diverge under 5' too. Both use the same assert-then-flip
//! idiom: a red assertion is a result, not a maintenance chore.

use crate::common::cis_apply_oracle::{apply, normalize, normalize_in};
use crate::common::synthetic::padded;
use ferro_hgvs::ShuffleDirection;

/// `(issue, core, split spelling, merged spelling)`.
const DIVERGENT: &[(&str, &str, &str, &str)] = &[
    (
        "#1287",
        "ATACAGAAAATCAGGGCATA",
        "TEMPLATE:g.[261_262insGA;263_264insAA]",
        "TEMPLATE:g.263_264insGAAA",
    ),
    (
        "#1290",
        "ATACAGAAAATCAGGGCATA",
        "TEMPLATE:g.[263_264insA;265_266insC]",
        "TEMPLATE:g.266_267insCA",
    ),
    (
        "#1296",
        "AAAAAAATAATCGCAACAGAAG",
        "TEMPLATE:g.[272_273del;274_275insAA]",
        "TEMPLATE:g.273delinsA",
    ),
    (
        "#1301",
        "GCATGAAAAT",
        "TEMPLATE:g.[263_264insAC;264_265insAA]",
        "TEMPLATE:g.264_265insCAAA",
    ),
    (
        "#1304",
        "GCATGAAAAT",
        "TEMPLATE:g.[260_261insGA;261_262insA;264del]",
        "TEMPLATE:g.262_263insGA",
    ),
    (
        "#1308",
        "CAGAAGATGAATAA",
        "TEMPLATE:g.[263_264insTG;264_265insTG]",
        "TEMPLATE:g.265_266insTTGG",
    ),
    (
        "#1312",
        "TAAAACCA",
        "TEMPLATE:g.[260_261insAC;261_262insAC]",
        "TEMPLATE:g.262_263insAACC",
    ),
    (
        "#1320",
        "AACAGTAAAATAT",
        "TEMPLATE:g.[263_264insAC;265_266insAA;266_267insAA]",
        "TEMPLATE:g.264_265insCAAAAA",
    ),
];

/// `(issue, core, split spelling, merged spelling)` for pairs that agree today.
const CONVERGED: &[(&str, &str, &str, &str)] = &[
    (
        "#1286",
        "AAAAAA",
        "TEMPLATE:g.[258_259insA;259_260insA]",
        "TEMPLATE:g.263_264insAA",
    ),
    (
        "#1297",
        "GCATGAAAAT",
        "TEMPLATE:g.[261_262insAA;263del;264_265insA]",
        "TEMPLATE:g.265_266insAA",
    ),
    (
        "#1316",
        "CAGCCAGTCAGCGCATCAG",
        "TEMPLATE:g.[261_262insAA;262_263insAA]",
        "TEMPLATE:g.262_263insAAAA",
    ),
    (
        "#1321",
        "TCCCAGAAAAT",
        "TEMPLATE:g.[261_262insGA;262_263insA;263del]",
        "TEMPLATE:g.263_264insGA",
    ),
    (
        "#1323",
        "CAGGGATCAT",
        "TEMPLATE:g.[260del;261_262insGA;262_263insGA]",
        "TEMPLATE:g.262_263insAGA",
    ),
];

/// Both spellings in every row must denote the same sequence, or the row is a
/// broken pin rather than evidence about the normalizer.
#[test]
fn every_pinned_pair_denotes_one_variant() {
    for (issue, core, split, merged) in DIVERGENT.iter().chain(CONVERGED) {
        let seq = padded(core);
        let a = apply(&seq, split)
            .unwrap_or_else(|| panic!("{issue}: split spelling `{split}` does not apply"));
        let b = apply(&seq, merged)
            .unwrap_or_else(|| panic!("{issue}: merged spelling `{merged}` does not apply"));
        assert_eq!(
            a, b,
            "{issue}: `{split}` and `{merged}` are not the same variant"
        );
    }
}

#[test]
fn the_eight_spelling_pairs_still_diverge() {
    // The count "eight" is asserted in three places' prose — this test's name,
    // the module doc above, and `tests/it/rewrite_target_corpus.rs` — and in none
    // of them executably. Adding or removing a row would leave all three wrong
    // and silent. `splitter_reproducer_corpus.rs` guards its own table the same
    // way.
    assert_eq!(
        DIVERGENT.len(),
        8,
        "row count changed; update this test's name, the module doc, and \
         tests/it/rewrite_target_corpus.rs's reference to these eight pairs"
    );
    for (issue, core, split, merged) in DIVERGENT {
        let seq = padded(core);
        let (a, b) = (normalize(&seq, split), normalize(&seq, merged));
        assert_ne!(
            a, b,
            "{issue} appears fixed — `{split}` and `{merged}` now both normalize to \
             `{a}`. Re-bless {issue}'s expectation to the converged form and delete \
             this row."
        );
    }
}

#[test]
fn converged_pairs_stay_converged() {
    for (issue, core, split, merged) in CONVERGED {
        let seq = padded(core);
        let (a, b) = (normalize(&seq, split), normalize(&seq, merged));
        assert_eq!(
            a, b,
            "{issue} regressed — `{split}` -> `{a}` but `{merged}` -> `{b}`; these \
             two spellings of one variant must agree."
        );
    }
}

/// Rows from [`CONVERGED`] that converge under the 3' rule but **not** under 5'.
///
/// Confluence is a property of the normalizer, not of one shuffle direction, and
/// `--direction 5prime` is a supported public option on both the CLI
/// (`src/bin/ferro.rs`) and the Python bindings. Every row above was blessed
/// against the 3' direction only, so these gaps were never measured.
///
/// **#1321 has left this set.** Its 5' divergence was the cancelled-member
/// residue: `g.[262_263insA;263del]` merges to `g.263delinsA`, which restates
/// the reference and so renders as `g.263=`, and the split spelling kept that
/// `=` while the merged spelling never grew one. Dropping the residue where the
/// merge creates it converges the pair on `g.261_262dup` — see
/// `issue_1321_identity_inside_a_duplication.rs`.
///
/// Both spellings still denote the input's bases and both are stable fixed
/// points, so this is criterion 1 of #1235 — non-confluence — and nothing else.
/// That is why neither oracle sees it: `FERRO_ASSERT_IDEMPOTENT` re-normalizes a
/// single spelling and finds it stable, and `FERRO_ASSERT_REPARSE` finds it
/// well-formed. Only comparing two spellings of one variant exposes it.
const DIVERGENT_UNDER_FIVE_PRIME: &[&str] = &["#1323"];

#[test]
fn the_five_prime_confluence_gap_is_unchanged() {
    // Assert-then-flip, the same contract the rows above carry: each listed row
    // is asserted to still DIVERGE, and each unlisted row to still CONVERGE. A
    // red assertion here is a *result*, not a maintenance chore — move the row
    // between the two sets rather than widening the list.
    //
    // The cause is not the shuffle direction reaching the derivation, which was
    // measured and refuted: threading `ShuffleDirection` into
    // `canonicalize_from_sequence` changes 298 of 39,600 swept outputs and all of
    // them are sequence-preserving both before and after, so it is a lateral
    // re-spelling with no correctness content. For these two rows the derivation
    // instead *refuses*: `canonicalize_from_sequence` returns `None`, behind
    // three stacked gates — the `collect_canonical_edits` catch-all (reached via
    // a redundant `=` member the 5' pipeline emits), then the
    // `changed_columns_of_pieces` weight bound, then `needs_unsupported_form`.
    // Each masks the next, which is the pattern #1235 and #1345 both describe.
    //
    // Every listed id must still name a `CONVERGED` row. The loop below reaches
    // an entry only through that table, so an id naming no row is inert: it
    // asserts nothing, and the set goes on claiming to pin a divergence that is
    // no longer measured. Deleting a row and leaving its id here is exactly how
    // that happens, and it is silent — verified by adding a bogus id, which left
    // all six tests in this file green. Same class of drift as the row count
    // pinned in `the_eight_spelling_pairs_still_diverge`.
    for issue in DIVERGENT_UNDER_FIVE_PRIME {
        assert!(
            CONVERGED.iter().any(|(row, ..)| row == issue),
            "{issue} is listed in DIVERGENT_UNDER_FIVE_PRIME but names no CONVERGED \
             row, so it pins nothing. Delete the entry, or restore the row it names."
        );
    }
    for (issue, core, split, merged) in CONVERGED {
        let seq = padded(core);
        let (a, b) = (
            normalize_in(&seq, split, ShuffleDirection::FivePrime),
            normalize_in(&seq, merged, ShuffleDirection::FivePrime),
        );
        if DIVERGENT_UNDER_FIVE_PRIME.contains(issue) {
            assert_ne!(
                a, b,
                "{issue} now CONVERGES under 5' shuffle — both spellings reach \
                 `{a}`. Remove it from DIVERGENT_UNDER_FIVE_PRIME; the gap is \
                 closing."
            );
        } else {
            assert_eq!(
                a, b,
                "{issue} has started diverging under 5' shuffle — `{split}` -> \
                 `{a}` but `{merged}` -> `{b}`. Two spellings of one variant must \
                 agree in both directions; add it to DIVERGENT_UNDER_FIVE_PRIME \
                 only with a measured reason."
            );
        }
    }
}

#[test]
fn the_eight_spelling_pairs_still_diverge_under_five_prime() {
    // The 5' half of `the_eight_spelling_pairs_still_diverge`. Those eight rows
    // were blessed against the 3' rule only, so without this the DIVERGENT table
    // had no 5' coverage at all and the module doc's "both directions" claim
    // covered only CONVERGED.
    //
    // All eight diverge under 5' as well, but they reach *different* spellings
    // than under 3' (#1296's merged form settles at `g.273C>A` here against
    // `g.273delinsA` there), so this is measuring the 5' pipeline rather than
    // restating the 3' result. Assert-then-flip: a row that starts converging is
    // the gap closing, and belongs in the message below rather than deleted.
    for (issue, core, split, merged) in DIVERGENT {
        let seq = padded(core);
        let (a, b) = (
            normalize_in(&seq, split, ShuffleDirection::FivePrime),
            normalize_in(&seq, merged, ShuffleDirection::FivePrime),
        );
        assert_ne!(
            a, b,
            "{issue} now CONVERGES under 5' shuffle — `{split}` and `{merged}` both \
             reach `{a}`. The gap is closing; re-bless {issue} and move it to \
             CONVERGED if it converges under 3' too."
        );
    }
}

#[test]
fn the_five_prime_divergences_are_non_confluence_and_nothing_worse() {
    // The severity claim in DIVERGENT_UNDER_FIVE_PRIME's doc, made executable. If
    // one of these ever starts changing the denoted sequence or ceases to be a
    // fixed point, that is a different and much worse defect than the one
    // recorded here, and it must not hide behind a known non-confluence.
    //
    // Every output is now checked directly, and two changes were needed to get
    // there. The one output that could not be put to the applier was `#1321`'s,
    // which carried a cancelled-member `=` residue: #1362 gave an identity's SPDI
    // triple the span it claims, so it no longer collides with a sibling's
    // insertion junction and such an output is expressible at all; #1379 then
    // stopped the residue being emitted, so this corpus no longer produces one.
    // The identity-stripping rewrite and the tolerated-skip count this loop used
    // to carry are both gone, and an inexpressible output is now a failure rather
    // than a known exception.
    for (issue, core, split, merged) in CONVERGED {
        if !DIVERGENT_UNDER_FIVE_PRIME.contains(issue) {
            continue;
        }
        let seq = padded(core);
        let want = apply(&seq, split).expect("input applies");
        for input in [split, merged] {
            let out = normalize_in(&seq, input, ShuffleDirection::FivePrime);
            assert_eq!(
                normalize_in(&seq, &out, ShuffleDirection::FivePrime),
                out,
                "{issue}: `{out}` is no longer a fixed point"
            );
            let got = apply(&seq, &out).unwrap_or_else(|| {
                panic!(
                    "{issue}: `{out}` is inexpressible to the applier. Since #1362 \
                     every output here is expressible, and since #1379 none carries \
                     a cancelled-member `=` residue at all, so this is a new blind \
                     spot rather than a known one — find what shape the applier \
                     declines before treating it as tolerable."
                )
            });
            assert_eq!(
                got, want,
                "{issue}: `{input}` -> `{out}` no longer denotes the input's bases"
            );
        }
    }
}
