//! An insertion that 3'-shifts into an identity member's span leaves an allele
//! that denotes no sequence (#2281).
//!
//! ```text
//! local 370..390  GGCCCT CCGCGG TGCCCCTGA      376_381 = CCGCGG, its own reverse complement
//!
//! g.[375_376insCA;376_381inv;385_386del]
//!   before  g.[376_377insAC;376_381=;384_387C[2]]
//!   after   g.[376_377insAC;384_387C[2]]
//! ```
//!
//! The palindromic `inv` normalizes to `376_381=`. The insertion at its 5' edge
//! 3'-shifts (`T|CA|CC` == `TC|AC|C`) to junction 376|377, which is inside the
//! identity's span. The applier then has no resulting sequence for the allele,
//! and `FERRO_ASSERT_SEQUENCE` fires on it. The input applies cleanly.
//!
//! `merge::drop_identity_members_covered_by_siblings` already dropped an
//! identity that a shift-blocking sibling overlaps (#1297, #1321, #1416). It
//! now also drops one with an insertion's junction strictly inside it. An
//! insertion at either edge still leaves the identity alone, and the tests
//! below pin both sides of that line.
//!
//! These cases use the GRCh38 window `adjacency_confluence_proptest` edits, so
//! they are the shape that property generated. That property cannot catch this
//! without `FERRO_ASSERT_SEQUENCE`: its soundness half compares sequences only
//! when the output applies, and an output that denotes nothing does not.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer, ShuffleDirection};

use crate::common::cis_apply_oracle::apply_with;
use crate::common::hg38_window::{local_desc, provider, HG38_WINDOW};

/// Normalize `body` against the window in strict mode, assert the output
/// denotes the same sequence as the input under the independent applier, and
/// return the output.
fn normalized_preserving(body: &str, direction: ShuffleDirection) -> String {
    let input = local_desc(body);
    let variant = parse_hgvs(&input).expect("input parses");
    let output = Normalizer::with_config(
        provider(),
        NormalizeConfig::default()
            .with_error_mode(ErrorMode::Strict)
            .with_direction(direction),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("{input} was refused: {e}"))
    .to_string();
    let p = provider();
    let before = apply_with(&p, HG38_WINDOW, &input);
    assert!(before.is_some(), "the input {input} must denote a sequence");
    assert_eq!(
        before,
        apply_with(&p, HG38_WINDOW, &output),
        "{input} normalized to {output}, which does not denote the input's sequence"
    );
    output
}

#[test]
fn the_reported_allele_drops_the_identity_the_insertion_lands_in() {
    for body in [
        "[375_376insCA;376_381inv;385_386del]",
        "[385_386del;376_381inv;375_376insCA]",
    ] {
        assert_eq!(
            normalized_preserving(body, ShuffleDirection::ThreePrime),
            "NC_TESTWIN.1:g.[376_377insAC;384_387C[2]]",
            "{body}"
        );
    }
}

#[test]
fn a_palindromic_inversion_beside_an_insertion() {
    assert_eq!(
        normalized_preserving("[375_376insCA;376_381inv]", ShuffleDirection::ThreePrime),
        "NC_TESTWIN.1:g.376_377insAC"
    );
}

#[test]
fn an_authored_identity_with_an_insertion_shifting_into_it() {
    // The identity need not come from an `inv`: an authored `=` takes the same
    // path, and so does the broken output itself, which used to be a fixed
    // point.
    for body in ["[375_376insCA;376_381=]", "[376_377insAC;376_381=]"] {
        let output = normalize_strict(body);
        assert_eq!(output, "NC_TESTWIN.1:g.376_377insAC", "{body}");
    }
    assert_eq!(
        normalized_preserving("[375_376insCA;376_381=]", ShuffleDirection::ThreePrime),
        "NC_TESTWIN.1:g.376_377insAC"
    );
}

#[test]
fn only_the_identity_with_the_insertion_inside_it_is_dropped() {
    // A sibling outside the identity does not stop the drop, and is kept.
    assert_eq!(
        normalized_preserving(
            "[375_376insCA;376_381=;390_391insT]",
            ShuffleDirection::ThreePrime
        ),
        "NC_TESTWIN.1:g.[376_377insAC;390_391insT]"
    );
}

/// The discriminating side: an insertion at an identity's **edge** is outside
/// the bases it names, so the identity stays. A rule that dropped any identity
/// an insertion touches would pass every test above and fail these.
#[test]
fn an_insertion_at_an_identitys_edge_leaves_it_alone() {
    // 3' edge, the insertion having shifted there from the 5' side.
    assert_eq!(
        normalized_preserving("[375_376insCA;376=]", ShuffleDirection::ThreePrime),
        "NC_TESTWIN.1:g.[376=;376_377insAC]"
    );
    // 3' edge of a multi-base identity, authored there.
    assert_eq!(
        normalized_preserving("[376_381=;381_382insA]", ShuffleDirection::ThreePrime),
        "NC_TESTWIN.1:g.[376_381=;381_382insA]"
    );
    // And an identity with no insertion near it is kept, as before.
    assert_eq!(
        normalized_preserving("[376_381inv;390del]", ShuffleDirection::ThreePrime),
        "NC_TESTWIN.1:g.[376_381=;390del]"
    );
}

/// Plain strict-mode normalization, for inputs that themselves denote no
/// sequence and so cannot go through [`normalized_preserving`].
fn normalize_strict(body: &str) -> String {
    let variant = parse_hgvs(&local_desc(body)).expect("input parses");
    Normalizer::with_config(
        provider(),
        NormalizeConfig::default().with_error_mode(ErrorMode::Strict),
    )
    .normalize(&variant)
    .unwrap_or_else(|e| panic!("{body} was refused: {e}"))
    .to_string()
}
