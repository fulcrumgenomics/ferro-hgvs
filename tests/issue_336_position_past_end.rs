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
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

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

    let msg = format!("{err:?}");
    assert!(
        msg.contains("past")
            || msg.contains("PositionPastEnd")
            || msg.contains("out of")
            || msg.contains("OutOfBounds"),
        "expected error message to mention past-end / out-of-bounds, got: {msg}",
    );
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
    let msg = format!("{err:?}");
    assert!(
        msg.contains("past")
            || msg.contains("PositionPastEnd")
            || msg.contains("out of")
            || msg.contains("OutOfBounds"),
        "expected past-end / out-of-bounds message, got: {msg}",
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
    let msg = format!("{err:?}");
    assert!(
        msg.contains("past")
            || msg.contains("PositionPastEnd")
            || msg.contains("out of")
            || msg.contains("OutOfBounds"),
        "expected past-end / out-of-bounds message, got: {msg}",
    );
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
fn _strict_mode_compiles_with_default_constructor() {
    // Sanity: the test fixtures rely on `NormalizeConfig::strict()` / `::lenient()`
    // / `::silent()` existing; this guards future renames.
    let _strict = NormalizeConfig::strict();
    let _lenient = NormalizeConfig::lenient();
    let _silent = NormalizeConfig::silent();
    // Also check that the underlying ErrorMode enum is reachable from tests.
    let _modes = [ErrorMode::Strict, ErrorMode::Lenient, ErrorMode::Silent];
}
