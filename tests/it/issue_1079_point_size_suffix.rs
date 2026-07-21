//! MUST-level rejection of the size-suffix form on a single-position anchor
//! (`g.123del3`, `g.123dup6`, `r.123del6`, `p.Arg45del6`) — issue #1079.
//!
//! # Spec
//!
//! `recommendations/checklist.md:49`:
//!
//! > Descriptions like `g.123del3` are not allowed, correct is `g.123_125del`.
//!
//! `recommendations/DNA/deletion.md:117`:
//!
//! > No, a deletion of more than one residue should mention the first and last
//! > residue deleted, separated using the range symbol ("_", underscore), e.g.,
//! > `NG_012232.1:g.123_128del` and not `NG_012232.1:g.123del6`.
//!
//! `recommendations/RNA/deletion.md:49` states the same for `r.123del6`, and
//! `recommendations/DNA/duplication.md:140` for `g.123dup6`.
//!
//! "not allowed" is a MUST-level prohibition under the spec's RFC 2119 reading
//! (`recommendations/style.md:9`), so the default (strict) [`parse_hgvs`] entry
//! point rejects it on every axis, including the protein axis (`p.Arg45del6`).
//!
//! # Deletion is repairable; duplication is not
//!
//! For a deletion the spec states the canonical replacement outright — a size
//! `N` at position `p` means `p`..`p+N-1` — so lenient/silent modes rewrite
//! `g.123del3` to `g.123_125del` (W3011). For a duplication the spec
//! explicitly declines to disambiguate:
//!
//! > Note also that from the description `g.123dup6` it is not clear whether
//! > the duplication starts **at** position `g.123` (so `g.123_128dup`) or
//! > **after** position 123 (so `g.124_129dup`).
//!
//! so there is no safe rewrite and every mode rejects (W3023).
//!
//! # Scope
//!
//! Only a **single-position** anchor is affected. A range anchor already names
//! both endpoints (`c.100_102del3`), which is what the rule requires; the
//! redundant size there remains a W3011/W3016 concern, not a hard rejection.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

/// Every axis's spelling of "size suffix on a single-position anchor".
const POINT_DEL_SIZE: &[&str] = &[
    "NC_000023.11:g.123del3",
    "NG_012232.1:g.123del6",
    "NM_004006.2:c.76del3",
    "NR_002196.2:n.123del6",
    "NM_004006.3:r.123del6",
    "NC_012920.1:m.3243del4",
];

const POINT_DUP_SIZE: &[&str] = &[
    "NC_000023.11:g.123dup6",
    "NM_004006.2:c.76dup3",
    "NM_004006.3:r.123dup6",
];

// =====================================================================
// Default (strict) parse rejects the point size-suffix form
// =====================================================================

#[test]
fn default_parse_rejects_point_del_size_suffix() {
    for input in POINT_DEL_SIZE {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "checklist.md:49 forbids {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn default_parse_rejects_point_dup_size_suffix() {
    for input in POINT_DUP_SIZE {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "duplication.md:140 forbids {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn default_parse_rejects_protein_point_del_count() {
    let result = parse_hgvs("NP_003997.1:p.Arg45del6");
    assert!(
        result.is_err(),
        "a single-residue anchor with a size count names only one endpoint; got {:?}",
        result.map(|v| v.to_string())
    );
}

#[test]
fn rejection_names_the_canonical_range_form() {
    let msg = parse_hgvs("NC_000023.11:g.123del3")
        .unwrap_err()
        .to_string();
    assert!(
        msg.contains("123_125del"),
        "the diagnostic should offer the canonical repair; got: {msg}"
    );
}

// =====================================================================
// Lenient / silent modes canonicalise the deletion form
// =====================================================================

#[test]
fn lenient_canonicalises_point_del_size_suffix() {
    for (input, expected) in [
        ("NC_000023.11:g.123del3", "NC_000023.11:g.123_125del"),
        ("NG_012232.1:g.123del6", "NG_012232.1:g.123_128del"),
        ("NM_004006.2:c.76del3", "NM_004006.2:c.76_78del"),
        ("NC_012920.1:m.3243del4", "NC_012920.1:m.3243_3246del"),
    ] {
        let out = parse_hgvs_with_config(input, ErrorConfig::lenient())
            .unwrap_or_else(|e| panic!("lenient must repair {input:?}: {e}"));
        assert_eq!(out.result.to_string(), expected, "input {input:?}");
    }
}

#[test]
fn size_one_deletion_is_redundant_not_forbidden() {
    // `DNA/deletion.md:117` scopes the rule to a deletion "of more than one
    // residue". `g.123del1` names exactly the anchor residue, so it is
    // redundant rather than incomplete and stays a soft W3011 concern.
    assert!(parse_hgvs("NC_000023.11:g.123del1").is_ok());
}

#[test]
fn silent_canonicalises_point_del_size_suffix() {
    let out = parse_hgvs_with_config("NG_012232.1:g.123del6", ErrorConfig::silent())
        .expect("silent must repair the size-suffix form");
    assert_eq!(out.result.to_string(), "NG_012232.1:g.123_128del");
}

#[test]
fn lenient_rejects_point_dup_size_suffix() {
    // duplication.md:140 declares the form ambiguous, so there is no safe
    // repair — lenient rejects rather than guessing an interpretation.
    for input in POINT_DUP_SIZE {
        assert!(
            parse_hgvs_with_config(input, ErrorConfig::lenient()).is_err(),
            "lenient must not guess an interpretation for {input:?}"
        );
    }
}

#[test]
fn compound_allele_members_are_handled_like_standalone() {
    // Each `del<N>` member of a compound allele is treated the same as a
    // standalone description: strict rejects the size-suffix form, and
    // lenient canonicalises every member to the range form (#1079).
    let compound = "NM_004006.2:c.[100del3;200del4]";
    assert!(
        parse_hgvs(compound).is_err(),
        "strict must reject the size-suffix form in every allele member"
    );
    let out = parse_hgvs_with_config(compound, ErrorConfig::lenient())
        .unwrap_or_else(|e| panic!("lenient must repair every member of {compound:?}: {e}"));
    assert_eq!(
        out.result.to_string(),
        "NM_004006.2:c.[100_102del;200_203del]"
    );
}

// =====================================================================
// Canonical and range-anchored forms are untouched
// =====================================================================

#[test]
fn canonical_range_forms_still_parse() {
    for input in [
        "NG_012232.1:g.123_128del",
        "NG_012232.1:g.123_128dup",
        "NG_012232.1:g.123del",
        "NG_012232.1:g.123dup",
        "NP_003997.1:p.Arg45del",
        "NP_003997.1:p.Lys228_Met259del32",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} names both endpoints (or needs none) and must parse"
        );
    }
}

#[test]
fn range_anchored_size_suffix_is_not_hard_rejected() {
    // `c.100_102del3` names both endpoints, satisfying the rule this file
    // enforces. It stays a soft (W3011/W3016) concern.
    assert!(parse_hgvs("NM_004006.2:c.100_102del3").is_ok());
}

#[test]
fn oversized_point_del_size_suffix_does_not_overflow() {
    // The canonical range end is `position + size - 1`. A `u64::MAX` anchor
    // makes that arithmetic overflow, which must not panic the preprocessor's
    // point-del-size rewrite: the token is declined (no safe range to offer)
    // and the anchor is still rejected downstream. Exercised in every mode.
    let oversized = "NC_000023.11:g.18446744073709551615del3";
    for config in [ErrorConfig::lenient(), ErrorConfig::silent()] {
        assert!(
            parse_hgvs_with_config(oversized, config).is_err(),
            "{oversized:?} must be rejected, not panic"
        );
    }
    assert!(
        parse_hgvs(oversized).is_err(),
        "{oversized:?} must be rejected, not panic"
    );
}
