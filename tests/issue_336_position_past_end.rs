//! Regression test for issue #336: `normalize()` must reject `c.`/`n.`
//! variants whose positions lie past the transcript's CDS-end or
//! transcript-end. Biocommons rejects with `HGVSError`; ferro previously
//! accepted and could even shift the position further past the end.
//!
//! Concrete biocommons cases that motivated this issue (from #324's
//! `tests/fixtures/biocommons-normalize/baseline-failures/normalized.txt`):
//!   - NM_001001656.1:c.946G>C       (CDS-end = 945)
//!   - NM_001001656.1:c.946dup
//!   - NM_001001656.1:c.935_946del
//!
//! These tests construct a synthetic `NM_TEST.1` transcript with a short
//! CDS so the past-end cases can be verified without a real manifest.

use ferro_hgvs::error_handling::ErrorMode;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, FerroError, MockProvider, NormalizeConfig, Normalizer};

/// Single-exon synthetic transcript:
///   tx positions 1..20  (20 bases)
///   5'UTR: 1..3   (`AAA`)
///   CDS:   4..12  (CDS length = 9; `ATGAAATAG` — start codon then ATAG-ish ending with stop)
///   3'UTR: 13..20 (`CCCCCCCC`)
///
/// In c. coords:
///   c.-3..c.-1  → tx 1..3   (5'UTR)
///   c.1..c.9    → tx 4..12  (CDS; c.9 is the last valid plain-int position)
///   c.*1..c.*8  → tx 13..20 (3'UTR)
fn provider_with_short_cds_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "AAAATGAAATAGCCCCCCCC".to_string();
    assert_eq!(sequence.len(), 20, "fixture length must be 20");
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        sequence,
        Some(4),
        Some(12),
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    provider
}

/// Helper: assert `err` is a `FerroError::InvalidCoordinates` whose message
/// names the W4004 code and the offending accession + position. Pinning both
/// the variant kind and the message contents catches accidental renames in
/// `FerroError` and silent reshapes of the error-message format.
fn assert_past_end_error(err: &FerroError, expected_position: &str) {
    match err {
        FerroError::InvalidCoordinates { msg } => {
            assert!(
                msg.contains("W4004"),
                "expected W4004 code in message, got: {msg}",
            );
            assert!(
                msg.contains("NM_TEST.1"),
                "expected accession NM_TEST.1 in message, got: {msg}",
            );
            assert!(
                msg.contains(expected_position),
                "expected position {expected_position} in message, got: {msg}",
            );
        }
        other => panic!("expected FerroError::InvalidCoordinates, got: {other:?}"),
    }
}

#[test]
fn strict_mode_rejects_c_position_past_cds_end() {
    // CDS length = 9, so c.10 is past the CDS end. Strict mode must reject
    // with a typed FerroError so callers can flag the input.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.10G>C").expect("parse");

    let err = normalizer
        .normalize(&variant)
        .expect_err("c.10 (past CDS-end 9) must reject in strict mode");
    assert_past_end_error(&err, "c.10");
}

#[test]
fn strict_mode_rejects_c_position_past_cds_end_for_dup() {
    // Mirrors the biocommons `NM_001001656.1:c.946dup` case. Duplications
    // go through `needs_normalization`, but the bounds check runs *before*
    // that short-circuit, so dup at past-end positions is rejected the
    // same way as substitutions and deletions.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.10dup").expect("parse");

    let err = normalizer
        .normalize(&variant)
        .expect_err("c.10dup (past CDS-end 9) must reject in strict mode");
    assert_past_end_error(&err, "c.10");
}

#[test]
fn lenient_mode_warns_on_c_position_past_cds_end() {
    // Lenient mode accepts the input but emits a PositionPastEnd warning.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs("NM_TEST.1:c.10G>C").expect("parse");

    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("lenient mode must not error on past-end positions");

    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected POSITION_PAST_END warning, got: {:?}",
        result.warnings
    );
}

#[test]
fn silent_mode_accepts_c_position_past_cds_end_without_warning() {
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::silent());
    let variant = parse_hgvs("NM_TEST.1:c.10G>C").expect("parse");

    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("silent mode must not error");

    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "silent mode must not emit POSITION_PAST_END warning, got: {:?}",
        result.warnings,
    );
}

#[test]
fn strict_mode_accepts_c_position_at_cds_end() {
    // c.9 is the last valid plain-integer CDS position.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.9G>C").expect("parse");
    normalizer
        .normalize(&variant)
        .expect("c.9 must pass strict-mode bounds check");
}

#[test]
fn strict_mode_rejects_range_end_past_cds_end() {
    // c.5_10del — start in-bounds, end past CDS-end. Strict rejects.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.5_10del").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("range ending past CDS-end must reject in strict mode");
    assert_past_end_error(&err, "c.10");
}

#[test]
fn strict_mode_rejects_range_both_endpoints_past_cds_end() {
    // c.10_11del — both endpoints past CDS-end. Both produce warnings;
    // strict mode reports the first one (the start position).
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.10_11del").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("c.10_11del (both past CDS-end) must reject in strict mode");
    // The first past-end warning (start = c.10) is the one promoted to the error.
    assert_past_end_error(&err, "c.10");
}

#[test]
fn lenient_mode_emits_two_warnings_for_range_both_past_end() {
    // Same input as above, but lenient mode keeps both warnings in the
    // result vec rather than promoting only the first to an error.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs("NM_TEST.1:c.10_11del").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("lenient mode must not error");
    let past_end: Vec<_> = result
        .warnings
        .iter()
        .filter(|w| w.code() == "POSITION_PAST_END")
        .collect();
    assert_eq!(
        past_end.len(),
        2,
        "lenient mode must emit one POSITION_PAST_END warning per offending endpoint; got: {:?}",
        result.warnings,
    );
}

#[test]
fn strict_mode_rejects_utr3_position_past_transcript_end() {
    // 3'UTR length = 8 (tx 13..20). c.*9 is past the transcript end.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.*9G>C").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("c.*9 (past transcript-end) must reject in strict mode");
    assert_past_end_error(&err, "c.*9");
}

#[test]
fn strict_mode_accepts_utr3_position_at_transcript_end() {
    // c.*8 → tx 20 (last base of transcript). In-bounds.
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.*8G>C").expect("parse");
    normalizer
        .normalize(&variant)
        .expect("c.*8 must pass strict-mode bounds check");
}

#[test]
fn strict_mode_skips_check_when_transcript_has_no_cds_bounds() {
    // Non-coding transcript: no cds_start / cds_end. The bounds helper
    // conservatively returns None when CDS metadata is missing, so the
    // bounds check is skipped and normalization proceeds.
    //
    // A `c.<N>` variant against a transcript without CDS bounds would
    // typically be caught upstream (the parser still accepts it, but
    // downstream conversion fails), but the bounds check itself must not
    // panic. This test pins the conservative-skip path.
    let mut provider = MockProvider::new();
    let sequence = "AAAATGAAATAGCCCCCCCC".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NR_TEST.1".to_string(),
        Some("NCRNA".to_string()),
        Strand::Plus,
        sequence,
        None, // cds_start
        None, // cds_end
        vec![Exon::new(1, 1, len)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);

    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    // Use an `n.` variant to keep the input semantically consistent with
    // a non-coding transcript. The bounds check is currently c.-only, so
    // this exercises the "n. is out of scope" path: no warning emitted,
    // no rejection.
    let variant = parse_hgvs("NR_TEST.1:n.10G>C").expect("parse");
    let result = normalizer
        .normalize_with_warnings(&variant)
        .expect("non-coding transcript path must not panic");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "n. variant must not emit POSITION_PAST_END (out of scope for #336); got: {:?}",
        result.warnings,
    );
}

#[test]
fn config_constructors_are_callable() {
    // Sanity: the test fixtures rely on `NormalizeConfig::strict()` /
    // `::lenient()` / `::silent()` existing; this guards against future
    // renames or accidental removal of the convenience constructors.
    let _strict = NormalizeConfig::strict();
    let _lenient = NormalizeConfig::lenient();
    let _silent = NormalizeConfig::silent();
    // Also check that the underlying ErrorMode enum is reachable from tests.
    let _modes = [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent];
}
