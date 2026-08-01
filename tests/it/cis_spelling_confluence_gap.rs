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

use crate::common::cis_apply_oracle::{apply, normalize};
use crate::common::synthetic::padded;

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
