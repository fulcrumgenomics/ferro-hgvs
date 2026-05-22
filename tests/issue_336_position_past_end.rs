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
    assert_past_end_error_with_accession(err, "NM_TEST.1", expected_position);
}

/// Variant that takes an explicit accession so the same strict structural
/// checks (variant kind, W-code, accession, full coordinate token) apply to
/// non-`NM_TEST.1` fixtures (`NM_5UTR.1`, `NR_TEST.1`, …).
fn assert_past_end_error_with_accession(
    err: &FerroError,
    accession: &str,
    expected_position: &str,
) {
    match err {
        FerroError::InvalidCoordinates { msg } => {
            assert!(
                msg.contains("W4004"),
                "expected W4004 code in message, got: {msg}",
            );
            assert!(
                msg.contains(accession),
                "expected accession {accession} in message, got: {msg}",
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
        .normalize_with_diagnostics(&variant)
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
        .normalize_with_diagnostics(&variant)
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
        .normalize_with_diagnostics(&variant)
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
fn strict_mode_n_in_bounds_does_not_emit_warning() {
    // Non-coding transcript: no cds_start / cds_end. An in-bounds n.
    // position must pass the new check_tx_pos_past_end check silently
    // — no panic, no warning.
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NR_TEST.1:n.10G>C").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("in-bounds n. variant must not error");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "in-bounds n. variant must not emit POSITION_PAST_END; got: {:?}",
        result.warnings,
    );
}

// ----------------------------------------------------------------------------
// 5'UTR c.-N coverage (folded from #348)
// ----------------------------------------------------------------------------

/// Transcript with a real 5'UTR: tx 1-5 (5'UTR `AAAAA`), CDS 6-14 (`ATGAAATAG`),
/// 3'UTR 15-20 (`CCCCCC`). cds_start=6 → 5'UTR length = 5.
fn provider_with_short_5utr_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "AAAAAATGAAATAGCCCCCC".to_string();
    assert_eq!(sequence.len(), 20);
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NM_5UTR.1".to_string(),
        Some("UTR".to_string()),
        Strand::Plus,
        sequence,
        Some(6),
        Some(14),
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
fn strict_mode_rejects_c_minus_position_past_5utr_start() {
    // 5'UTR length = 5, so `c.-6` is past the 5'UTR start. Strict rejects.
    let provider = provider_with_short_5utr_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_5UTR.1:c.-6G>C").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("c.-6 (past 5'UTR start, length 5) must reject in strict mode");
    assert_past_end_error_with_accession(&err, "NM_5UTR.1", "c.-6");
    // The 5'UTR-start bound kind has a distinct tag in the message —
    // pin it separately so the n. / cds-end / 5utr-start arms stay
    // distinguishable.
    let FerroError::InvalidCoordinates { msg } = &err else {
        unreachable!()
    };
    assert!(
        msg.contains("5utr-start"),
        "expected 5utr-start bound-kind tag, got: {msg}",
    );
}

#[test]
fn strict_mode_accepts_c_minus_position_at_5utr_start() {
    // `c.-5` = the first 5'UTR base. In-bounds.
    let provider = provider_with_short_5utr_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_5UTR.1:c.-5G>C").expect("parse");
    normalizer
        .normalize(&variant)
        .expect("c.-5 must pass strict-mode bounds check");
}

#[test]
fn lenient_mode_warns_on_c_minus_position_past_5utr_start() {
    let provider = provider_with_short_5utr_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs("NM_5UTR.1:c.-6G>C").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must not error");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected POSITION_PAST_END warning, got: {:?}",
        result.warnings,
    );
}

// ----------------------------------------------------------------------------
// n. transcript variant coverage (folded from #347)
// ----------------------------------------------------------------------------

/// Non-coding transcript: 20 bases, no cds_start/cds_end.
fn provider_with_short_noncoding_transcript() -> MockProvider {
    let mut provider = MockProvider::new();
    let sequence = "ACGTACGTACGTACGTACGT".to_string();
    let len = sequence.len() as u64;
    let transcript = Transcript::new(
        "NR_TEST.1".to_string(),
        Some("NCRNA".to_string()),
        Strand::Plus,
        sequence,
        None,
        None,
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
fn strict_mode_rejects_n_position_past_transcript_end() {
    // Transcript length = 20, so n.21 is past the end. Strict rejects.
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NR_TEST.1:n.21G>C").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("n.21 (past transcript-end 20) must reject in strict mode");
    assert_past_end_error_with_accession(&err, "NR_TEST.1", "n.21");
    let FerroError::InvalidCoordinates { msg } = &err else {
        unreachable!()
    };
    assert!(
        msg.contains("transcript-end"),
        "expected transcript-end bound-kind tag, got: {msg}",
    );
}

#[test]
fn strict_mode_accepts_n_position_at_transcript_end() {
    // n.20 = the last base of the transcript. In-bounds.
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NR_TEST.1:n.20G>C").expect("parse");
    normalizer
        .normalize(&variant)
        .expect("n.20 must pass strict-mode bounds check");
}

#[test]
fn lenient_mode_warns_on_n_position_past_transcript_end() {
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs("NR_TEST.1:n.21G>C").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must not error");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected POSITION_PAST_END warning on n., got: {:?}",
        result.warnings,
    );
}

#[test]
fn strict_mode_rejects_n_range_end_past_transcript_end() {
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NR_TEST.1:n.15_25del").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("range ending past transcript-end must reject");
    // End position 25 is the offending one (start=15 is in-bounds).
    assert_past_end_error_with_accession(&err, "NR_TEST.1", "n.25");
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

// ----------------------------------------------------------------------------
// silent-mode coverage for the c.-N and n. branches (CodeRabbit nit)
// ----------------------------------------------------------------------------

#[test]
fn silent_mode_accepts_c_minus_position_past_5utr_start_without_warning() {
    let provider = provider_with_short_5utr_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::silent());
    let variant = parse_hgvs("NM_5UTR.1:c.-6G>C").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("silent mode must not error");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "silent mode must not emit POSITION_PAST_END for past-5UTR-start, got: {:?}",
        result.warnings,
    );
}

#[test]
fn silent_mode_accepts_n_position_past_transcript_end_without_warning() {
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::silent());
    let variant = parse_hgvs("NR_TEST.1:n.21G>C").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("silent mode must not error");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "silent mode must not emit POSITION_PAST_END for past-tx-end n., got: {:?}",
        result.warnings,
    );
}

// ----------------------------------------------------------------------------
// con (SVD-WG009) fast-path coverage: rewritten c.<past>conT / n.<past>conT
// must still hit the W4004 bounds gate (was a silent bypass before the
// recursion fix in normalize_cds / normalize_tx).
// ----------------------------------------------------------------------------

#[test]
fn strict_mode_rejects_c_position_past_cds_end_for_con_rewrite() {
    let provider = provider_with_short_cds_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    // `c.11conT` rewrites to `c.11delinsT`. CDS length is 9 (last valid: c.9);
    // position 11 must surface W4004 even though the input arrives via the
    // `con` fast path.
    let variant = parse_hgvs("NM_TEST.1:c.11conT").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("past-CDS-end con must reject");
    assert_past_end_error(&err, "c.11");
}

#[test]
fn strict_mode_rejects_n_position_past_transcript_end_for_con_rewrite() {
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    // `n.21conT` rewrites to `n.21delinsT`. transcript length 20 → position
    // 21 must surface W4004 via the n. bounds gate.
    let variant = parse_hgvs("NR_TEST.1:n.21conT").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("past-tx-end n. con must reject");
    assert_past_end_error_with_accession(&err, "NR_TEST.1", "n.21");
}

// ----------------------------------------------------------------------------
// `n.*N` downstream positions are out of scope for W4004 — same skip as
// `simple_tx_pos`, otherwise check_tx_pos_past_end would compare the raw
// downstream base against transcript length and emit a false positive.
// ----------------------------------------------------------------------------

#[test]
fn lenient_mode_skips_n_downstream_position() {
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    // `n.*5` is downstream notation; the bounds check must not fire on the
    // raw base value (5) — and equally must not fire on the offset into the
    // post-transcript region.
    let variant = parse_hgvs("NR_TEST.1:n.*5G>C").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("downstream n. must not error in lenient mode");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "POSITION_PAST_END must not fire on downstream n. pos, got: {:?}",
        result.warnings,
    );
}

#[test]
fn lenient_mode_preserves_downstream_n_indel_without_normalization() {
    // `n.*5_*7del` is a downstream del. The bounds-check helper skips
    // downstream positions, but `normalize_tx` must also short-circuit
    // them — otherwise it feeds the raw downstream base (5/7) into the
    // in-transcript shift logic and normalizes against the wrong window.
    // Expectation: variant passes through canonicalize-only, no
    // POSITION_PAST_END warning, and Display preserves the input form.
    let provider = provider_with_short_noncoding_transcript();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs("NR_TEST.1:n.*5_*7del").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("downstream n. indel must not error in lenient mode");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "POSITION_PAST_END must not fire on downstream n. indel, got: {:?}",
        result.warnings,
    );
    // The canonical form must keep the `*` (downstream) prefix on both
    // endpoints — i.e. the indel was not normalized against the in-
    // transcript window.
    let display = format!("{}", result.result);
    assert!(
        display.contains("*5") && display.contains("*7"),
        "downstream endpoints must survive normalization; got: {display}",
    );
}

// ----------------------------------------------------------------------------
// `Transcript::cds_length` / `utr3_length` / `utr5_length` /
// `sequence_length` all return sensible bounds even when no sequence bytes
// are cached (they fall back to `cds_end - cds_start + 1` / the exon-sum
// transcript length). The bounds gate must fire on coordinate-only
// transcripts — a coordinate check must not silently degrade into a
// sequence-availability check.
// ----------------------------------------------------------------------------

fn provider_with_short_cds_transcript_no_sequence() -> MockProvider {
    let mut provider = MockProvider::new();
    // Same structural layout as `provider_with_short_cds_transcript`
    // (CDS 4..12, length 9; exons span tx 1..20), but no cached sequence.
    let transcript = Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TEST".to_string()),
        Strand::Plus,
        None::<String>,
        Some(4),
        Some(12),
        vec![Exon::new(1, 1, 20)],
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

fn provider_with_short_noncoding_transcript_no_sequence() -> MockProvider {
    let mut provider = MockProvider::new();
    let transcript = Transcript::new(
        "NR_TEST.1".to_string(),
        Some("NCRNA".to_string()),
        Strand::Plus,
        None::<String>,
        None,
        None,
        vec![Exon::new(1, 1, 20)],
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
fn strict_mode_rejects_c_past_cds_end_on_coordinate_only_transcript() {
    // No cached sequence — but CDS 4..12 still implies length 9. The bounds
    // gate must fire regardless of whether bases are loaded.
    let provider = provider_with_short_cds_transcript_no_sequence();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.10G>C").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("c.10 past CDS-end 9 must reject even without cached sequence");
    assert_past_end_error(&err, "c.10");
}

#[test]
fn strict_mode_rejects_utr3_past_end_on_coordinate_only_transcript() {
    // No cached sequence — exon-sum transcript length is 20, CDS end 12,
    // so 3'UTR length = 8. c.*9 must reject.
    let provider = provider_with_short_cds_transcript_no_sequence();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NM_TEST.1:c.*9G>C").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("c.*9 past 3'UTR-end 8 must reject even without cached sequence");
    assert_past_end_error(&err, "c.*9");
}

#[test]
fn strict_mode_rejects_n_past_transcript_end_on_coordinate_only_transcript() {
    // No cached sequence — exon-sum transcript length is 20. n.21 must reject.
    let provider = provider_with_short_noncoding_transcript_no_sequence();
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs("NR_TEST.1:n.21G>C").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("n.21 past transcript-end 20 must reject even without cached sequence");
    assert_past_end_error_with_accession(&err, "NR_TEST.1", "n.21");
}
