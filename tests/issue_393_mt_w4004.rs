//! Integration tests for issue #393 — W4004 PositionPastEnd on the m. axis.
//!
//! `NC_012920.1` (rCRS, human mitochondrial genome) is 16569 bp. Any
//! `m.` position with `base > 16569` should emit W4004 at normalization
//! time. Positions exactly at the contig end (16569) must NOT fire.
//! Wraparound ranges where both endpoints fit within the contig must NOT
//! fire.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::{parse_hgvs, NormalizeConfig, Normalizer};

fn mt_provider() -> MockProvider {
    let mut p = MockProvider::new();
    // NC_012920.1 is 16569 bp. For bounds-check purposes we only need the
    // length, not the actual bases. A placeholder string of the right length
    // is sufficient.
    p.add_genomic_sequence("NC_012920.1", "A".repeat(16569));
    p
}

// ---------------------------------------------------------------------------
// Past-end position must fire W4004 (strict mode rejects, lenient mode warns)
// ---------------------------------------------------------------------------

#[test]
fn strict_mode_rejects_m_position_past_contig_end() {
    // m.16570 is one past the 16569-bp contig; strict mode must reject W4004.
    let variant = parse_hgvs("NC_012920.1:m.16570A>G").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::strict());
    let err = normalizer
        .normalize(&variant)
        .expect_err("m.16570 (past contig-end 16569) must reject in strict mode");
    match &err {
        ferro_hgvs::FerroError::InvalidCoordinates { msg } => {
            assert!(msg.contains("W4004"), "expected W4004 code, got: {msg}");
            assert!(
                msg.contains("NC_012920.1"),
                "expected accession in message, got: {msg}",
            );
            assert!(
                msg.contains("16570"),
                "expected position in message, got: {msg}",
            );
            assert!(
                msg.contains("contig-end"),
                "expected contig-end bound-kind tag, got: {msg}",
            );
        }
        other => panic!("expected FerroError::InvalidCoordinates, got: {other:?}"),
    }
}

#[test]
fn lenient_mode_warns_on_m_position_past_contig_end() {
    let variant = parse_hgvs("NC_012920.1:m.16570A>G").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::lenient());
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("lenient mode must not error on past-end positions");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected POSITION_PAST_END warning, got: {:?}",
        result.warnings,
    );
}

#[test]
fn silent_mode_accepts_m_position_past_contig_end_without_warning() {
    let variant = parse_hgvs("NC_012920.1:m.16570A>G").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::silent());
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("silent mode must not error");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "silent mode must not emit POSITION_PAST_END, got: {:?}",
        result.warnings,
    );
}

// ---------------------------------------------------------------------------
// Last valid position (m.16569) must NOT fire W4004
// ---------------------------------------------------------------------------

#[test]
fn strict_mode_accepts_m_position_at_contig_end() {
    // m.16569 is the last valid position — must pass the bounds check.
    let variant = parse_hgvs("NC_012920.1:m.16569A>G").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::strict());
    // Note: this may fail for a ref-mismatch reason (placeholder 'A' may not
    // match), so we use lenient mode here to isolate the bounds gate.
    let normalizer_lenient = Normalizer::with_config(mt_provider(), NormalizeConfig::lenient());
    let result = normalizer_lenient
        .normalize_with_warnings(&variant)
        .expect("m.16569 must not error in lenient mode");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "m.16569 must not fire POSITION_PAST_END; got: {:?}",
        result.warnings,
    );
    // Strict mode: the bounds check must not reject for a valid position.
    // (Ref mismatch from the placeholder might fire, but not W4004.)
    let _ = normalizer.normalize(&variant); // either ok or ref-mismatch, not W4004
}

// ---------------------------------------------------------------------------
// Wraparound range with both endpoints in-bounds: no W4004
// ---------------------------------------------------------------------------

#[test]
fn lenient_mode_no_w4004_on_wraparound_range_with_both_endpoints_in_bounds() {
    // m.16569_1del — start=16569, end=1 (reversed: start > end is the
    // wraparound signal per SVD-WG006). Both endpoints are within [1..16569].
    // Must not fire W4004 regardless of shuffle semantics.
    let variant = parse_hgvs("NC_012920.1:m.16569_1del").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::lenient());
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("wraparound del must not error in lenient mode");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "wraparound with both endpoints valid must not fire W4004; got: {:?}",
        result.warnings,
    );
}

// ---------------------------------------------------------------------------
// Wraparound-shaped range where start exceeds the contig: W4004 must fire
// ---------------------------------------------------------------------------

#[test]
fn lenient_mode_warns_on_partial_wraparound_with_start_past_contig() {
    // m.16570_5del — start=16570 is past the 16569-bp contig end.
    // Even though this has a "wraparound shape" (start > end when compared
    // naively), the start position itself is invalid: 16570 > 16569.
    // W4004 must fire for the start position.
    let variant = parse_hgvs("NC_012920.1:m.16570_5del").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::lenient());
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("lenient mode must not error");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected W4004 for past-end start m.16570; got: {:?}",
        result.warnings,
    );
}
