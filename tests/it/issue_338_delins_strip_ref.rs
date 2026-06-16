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

// Follow-up #352: `NaEdit::Inversion { sequence, length }` should
// follow the same minimal-notation rule as Deletion/Duplication/Delins
// — the spec recommends `g.100_105inv` over `g.100_105invATGCC` or
// `g.100_105inv5`. The Inversion arm in `canonicalize_edit` /
// `should_canonicalize` strips both fields.
#[test]
fn genome_inversion_strips_explicit_inverted_sequence() {
    // §HGVS v21.0 DNA/inversion.md: "the recommendation is not to
    // describe the inverted nucleotide sequence." Strip the explicit
    // bases regardless of provider availability.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233995invATGCC"),
        "NC_000009.11:g.36233991_36233995inv",
    );
}

#[test]
fn genome_inversion_strips_explicit_inverted_length() {
    // Numeric form `inv5` is the length variant of the same redundancy.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233995inv5"),
        "NC_000009.11:g.36233991_36233995inv",
    );
}

#[test]
fn cds_inversion_strips_explicit_inverted_sequence() {
    // Same rule on coding variants.
    assert_eq!(
        normalize_str("NM_TEST.1:c.100_105invATGCCG"),
        "NM_TEST.1:c.100_105inv",
    );
}

#[test]
fn genome_inversion_preserves_short_form_round_trip() {
    // Already-canonical `inv` must round-trip unchanged.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233995inv"),
        "NC_000009.11:g.36233991_36233995inv",
    );
}

#[test]
fn genome_delins_with_count_insert_strips_deleted_sequence() {
    // Non-literal InsertedSequence (Count) with explicit deleted bases —
    // surfaces the WITH-provider gap caught in code review. With a real
    // provider available, the Literal arm in normalize_na_edit's Delins
    // dispatch trims/canonicalizes, but the non-literal fall-through
    // previously returned the edit unchanged (preserving deleted bases).
    // The fix routes the non-literal pass-through through canonicalize_edit
    // so the same strip-redundant-deleted-bases rule applies.
    //
    // MockProvider has no genomic data, so this still exercises the
    // no-provider canonicalize-only path. The provider-backed counterpart
    // is below in `genome_delins_with_count_insert_strips_deleted_sequence_with_provider`.
    assert_eq!(
        normalize_str("NC_000009.11:g.36233991_36233992delCAins10"),
        "NC_000009.11:g.36233991_36233992delins10",
    );
}

#[test]
fn genome_delins_with_count_insert_strips_deleted_sequence_with_provider() {
    // Provider-backed regression for the non-literal `Delins` branch in
    // `normalize_na_edit` (the line-2067 fall-through). Seeded genomic
    // sequence makes `ref_bytes` non-empty so the WITH-provider path
    // executes, distinct from the no-provider canonicalize-only path
    // exercised above. Closes the CodeRabbit nit on #344.
    let mut provider = MockProvider::new();
    // Seed a synthetic genomic contig. Bases are indexed 1-based by
    // MockProvider::get_sequence. Positions 5 ('C') + 6 ('A') match the
    // explicit `delCA` so the variant is well-formed against the seeded
    // reference and the WITH-provider path runs.
    provider.add_genomic_sequence("NC_TEST.1", "AAAACATTT");
    let normalizer = Normalizer::with_config(
        provider,
        NormalizeConfig::default().with_direction(ShuffleDirection::ThreePrime),
    );
    let variant = parse_hgvs("NC_TEST.1:g.5_6delCAins10").expect("parse");
    let normalized = normalizer.normalize(&variant).expect("normalize");
    assert_eq!(format!("{}", normalized), "NC_TEST.1:g.5_6delins10");
}
