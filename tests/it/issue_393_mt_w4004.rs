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
            // Axis-aware wording: m. lives on a contig, not a transcript.
            assert!(
                msg.contains("contig") && !msg.contains("transcript"),
                "expected mt-axis wording to say `contig` not `transcript`; got: {msg}",
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
        .normalize_with_diagnostics(&variant)
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
        .normalize_with_diagnostics(&variant)
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
        .normalize_with_diagnostics(&variant)
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
    // A ref-mismatch failure is acceptable (placeholder bases may not match),
    // but the error message must NOT mention W4004 / POSITION_PAST_END.
    match normalizer.normalize(&variant) {
        Ok(_) => {}
        Err(ferro_hgvs::FerroError::InvalidCoordinates { msg }) => {
            assert!(
                !msg.contains("W4004") && !msg.contains("POSITION_PAST_END"),
                "m.16569 must not be treated as past-end in strict mode; got: {msg}",
            );
        }
        Err(_) => {
            // other strict-mode rejections (e.g. ref mismatch) are allowed
        }
    }
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
        .normalize_with_diagnostics(&variant)
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
        .normalize_with_diagnostics(&variant)
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

// ---------------------------------------------------------------------------
// Same-base endpoints with distinct offsets must still check both ends
// ---------------------------------------------------------------------------

#[test]
fn lenient_mode_w4004_when_only_offset_endpoint_distinguishes_positions() {
    // `m.16570+1_16570del` — start has a non-standard `+1` offset (skipped by
    // `check_mt_pos_past_end`'s offset guard), end is plain `m.16570` which
    // exceeds the 16569-bp contig. A naive `if start.base != end.base` dedupe
    // would skip the end position because the bases match, and W4004 would
    // be silently lost. The endpoint distinctness check must compare the
    // full `(base, offset)` tuple — mirrors the c./n. dedupe pattern.
    let variant = parse_hgvs("NC_012920.1:m.16570+1_16570del").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::lenient());
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must not error");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected W4004 for past-end plain end m.16570 (start has +1 offset); got: {:?}",
        result.warnings,
    );
}

// ---------------------------------------------------------------------------
// con (SVD-WG009) fast-path coverage: rewritten m.<past>conT must still hit
// the W4004 bounds gate. Mirrors the c./n. axis tests
// `strict_mode_rejects_c_position_past_cds_end_for_con_rewrite` and
// `strict_mode_rejects_n_position_past_transcript_end_for_con_rewrite` in
// `tests/issue_336_position_past_end.rs`. The bounds check must fire on the
// rewritten `delins` form so `m.16570conT` does not slip past W4004 via the
// `con` fast path.
// ---------------------------------------------------------------------------

#[test]
fn strict_mode_rejects_m_position_past_contig_end_for_con_rewrite() {
    // `m.16570conT` rewrites to `m.16570delinsT`. Position 16570 exceeds the
    // 16569-bp mitochondrial contig and must surface W4004 even though the
    // input arrives via the SVD-WG009 `con` fast path. Requires the
    // `con -> delins` rewrite to recurse BEFORE the past-end bounds gate so
    // the rewritten delins re-enters `normalize_mt` and the gate fires on
    // the canonical form. (Without recursion-first ordering the early-return
    // would emit W4004 but leave the variant in `con` form — non-canonical.)
    let variant = parse_hgvs("NC_012920.1:m.16570conT").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::strict());
    let err = normalizer
        .normalize(&variant)
        .expect_err("m.16570conT (past contig-end 16569) must reject in strict mode");
    match &err {
        ferro_hgvs::FerroError::InvalidCoordinates { msg } => {
            assert!(msg.contains("W4004"), "expected W4004 code, got: {msg}");
            assert!(
                msg.contains("16570"),
                "expected position 16570 in message, got: {msg}",
            );
        }
        other => panic!("expected FerroError::InvalidCoordinates, got: {other:?}"),
    }
}

#[test]
fn lenient_mode_canonicalizes_con_to_delins_even_when_past_contig_end() {
    // Lenient mode: `m.16570conT` must (a) emit W4004 for the past-end
    // position AND (b) canonicalize the `con` edit to its `delins` form.
    // Pins the recursion-first ordering of the `con -> delins` rewrite
    // relative to the bounds gate — if the gate runs first, the returned
    // variant is left in non-canonical `con` form.
    let variant = parse_hgvs("NC_012920.1:m.16570conT").expect("parse must succeed");
    let normalizer = Normalizer::with_config(mt_provider(), NormalizeConfig::lenient());
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must not error on past-end con");
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code() == "POSITION_PAST_END"),
        "expected POSITION_PAST_END warning for m.16570conT, got: {:?}",
        result.warnings,
    );
    let display = format!("{}", result.result);
    assert!(
        display.contains("delins") && !display.contains("con"),
        "expected con -> delins canonicalization in lenient mode; got: {display}",
    );
}
