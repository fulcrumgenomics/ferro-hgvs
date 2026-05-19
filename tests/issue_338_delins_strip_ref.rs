//! Regression test for issue #338: ferro's `normalize()` output for a
//! top-level `delins` retains the explicit deleted bases. The HGVS spec
//! recommends stripping them — §DNA/delins.md: "the recommendation is
//! not to describe the variant as `delTinsGA`, … this description is
//! longer, it contains redundant information, and chances to make an
//! error increase."
//!
//! Biocommons case (from #324 baseline-failures):
//!   `NC_000009.11:g.36233991_36233992delCAinsAC` (3prime+cross)
//!     expected: `…delinsAC`
//!     ferro:    `…delCAinsAC`

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

fn normalize_str(input: &str) -> String {
    // MockProvider has no genomic data — that's fine. The strip-ref
    // canonicalization is purely syntactic and must fire regardless of
    // whether reference data is available. (Other rules like 3'-shift
    // require bases; the strip-ref rule does not.)
    let provider = MockProvider::new();
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs(input).expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    format!("{}", normalized)
}

#[test]
fn genome_delins_strips_explicit_deleted_sequence() {
    // The motivating biocommons case from #324. The two-base delins's
    // alt `AC` is the canonical short form; the explicit `delCA` is
    // redundant and the spec says drop it.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233992delCAinsAC"),
        "NC_000009.11:g.36233991_36233992delinsAC",
    );
}

#[test]
fn genome_delins_strips_explicit_deleted_length() {
    // Numeric form: `del2insAC` should also drop the redundant length.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233992del2insAC"),
        "NC_000009.11:g.36233991_36233992delinsAC",
    );
}

#[test]
fn cds_delins_strips_explicit_deleted_sequence() {
    // Same rule applies to coding variants. The spec example from
    // §DNA/delins.md is `c.6775_6777delinsC` (canonical) vs
    // `c.6775_6777delGAGinsC` (not recommended).
    assert_eq!(
        normalize_str("NM_TEST.1:c.6775_6777delGAGinsC"),
        "NM_TEST.1:c.6775_6777delinsC",
    );
}

#[test]
fn genome_delins_preserves_short_form_round_trip() {
    // Already-canonical input must round-trip unchanged.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233992delinsAC"),
        "NC_000009.11:g.36233991_36233992delinsAC",
    );
}
