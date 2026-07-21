//! MUST-level rejection of a one-nucleotide inversion — `g.234inv`,
//! `r.234inv` (#1079).
//!
//! # Spec
//!
//! `recommendations/DNA/inversion.md:16`:
//!
//! > by definition, the region inverted (`positions_inverted`) contains
//! > **more than one nucleotide**. The description `g.234inv` is therefore
//! > not allowed; a one-nucleotide inversion should be described as a
//! > substitution.
//!
//! "not allowed" is MUST-level under the spec's RFC 2119 reading
//! (`recommendations/style.md:9`).
//!
//! # Why reject rather than canonicalise
//!
//! The spec's replacement is a substitution, whose reference and alternative
//! bases are the base at that position and its complement. Neither is
//! recoverable from the description alone — it takes the reference sequence,
//! which the parser does not have — so there is no canonical form to rewrite
//! to and every mode rejects.

use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
use ferro_hgvs::parse_hgvs;

const POINT_INVERSIONS: &[&str] = &[
    "NC_000023.11:g.234inv",
    "NM_004006.3:r.234inv",
    "NM_004006.2:c.76inv",
    "NR_002196.2:n.123inv",
    "NC_012920.1:m.3243inv",
];

#[test]
fn default_parse_rejects_single_nucleotide_inversion() {
    for input in POINT_INVERSIONS {
        let result = parse_hgvs(input);
        assert!(
            result.is_err(),
            "inversion.md:16 forbids {input:?}; got {:?}",
            result.map(|v| v.to_string())
        );
    }
}

#[test]
fn every_mode_rejects_single_nucleotide_inversion() {
    for config in [
        ErrorConfig::strict(),
        ErrorConfig::lenient(),
        ErrorConfig::silent(),
    ] {
        assert!(
            parse_hgvs_with_config("NC_000023.11:g.234inv", config).is_err(),
            "no mode may accept a one-nucleotide inversion"
        );
    }
}

#[test]
fn rejection_points_at_the_substitution_form() {
    let msg = parse_hgvs("NC_000023.11:g.234inv").unwrap_err().to_string();
    assert!(
        msg.contains("substitution"),
        "the diagnostic should name the substitution form; got: {msg}"
    );
}

/// A range that collapses to one nucleotide is the same description written
/// differently and is rejected too.
#[test]
fn degenerate_range_inversion_is_rejected() {
    assert!(parse_hgvs("NC_000023.11:g.234_234inv").is_err());
}

// =====================================================================
// Multi-nucleotide inversions — every spec example must keep parsing
// =====================================================================

#[test]
fn multi_nucleotide_inversions_still_parse() {
    for input in [
        "NC_000001.11:g.1234_2345inv",
        "NC_000023.10:g.32361330_32361333inv",
        "NM_004006.2:c.5657_5660inv",
        "NM_004006.3:r.123_345inv",
        "NM_003079.4:c.374_395inv",
    ] {
        assert!(
            parse_hgvs(input).is_ok(),
            "{input:?} inverts more than one nucleotide and must parse"
        );
    }
}

/// An inverted *inserted range* (`ins<a>_<b>inv`, the spec's inverted-
/// duplication form) is an insertion, not an inversion edit, so the rule
/// does not reach it even though its anchor spans two positions.
#[test]
fn inverted_insertion_source_is_untouched() {
    for input in [
        "NM_004006.3:r.234_235ins123_234inv",
        "NM_004006.2:c.849_850ins850_900inv",
    ] {
        assert!(parse_hgvs(input).is_ok(), "{input:?} must still parse");
    }
}
