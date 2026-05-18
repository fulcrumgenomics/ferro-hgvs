//! Tests for #266 — W3016 LengthMismatch soft-validation warning.
//!
//! ferro silently accepts `del/dup/inv` inputs where the explicit
//! reference sequence length does not match the position range length.
//! This PR adds the W3016 SVA code, the `detect_length_mismatch`
//! corrector, and wires it as preprocessor Phase 13a (`warn_accept`).

use ferro_hgvs::error::ErrorCode;
use ferro_hgvs::error_handling::ErrorConfig;
use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_with_config};

fn warning_codes(input: &str) -> Vec<String> {
    let r = parse_hgvs_lenient(input).unwrap_or_else(|e| panic!("lenient parse {input:?}: {e}"));
    r.warnings
        .iter()
        .map(|w| w.error_type.code().to_string())
        .collect()
}

#[track_caller]
fn assert_w3016(input: &str) {
    let codes = warning_codes(input);
    assert!(
        codes.iter().any(|c| c == "W3016"),
        "expected W3016 for {input:?}, got {codes:?}"
    );
}

#[track_caller]
fn assert_no_w3016(input: &str) {
    let codes = warning_codes(input);
    assert!(
        !codes.iter().any(|c| c == "W3016"),
        "unexpected W3016 for {input:?}, got {codes:?}"
    );
}

// =============================================================================
// SECTION 1 — Mismatch cases (lenient emits W3016)
// =============================================================================

mod mismatch_emits_w3016 {
    use super::*;

    /// 11 positions, 10 bases.
    #[test]
    fn del_long_range_short_ref() {
        assert_w3016("NC_000001.11:g.100_110delAAAATTTGCC");
    }

    /// 6 positions, 5 bases.
    #[test]
    fn inv_long_range_short_ref() {
        assert_w3016("NC_000001.11:g.100_105invTAGCA");
    }

    /// 3 positions, 5 bases (ref too long for range).
    #[test]
    fn dup_short_range_long_ref() {
        assert_w3016("NC_000001.11:g.100_102dupACGTA");
    }

    /// CDS coord system also covered.
    #[test]
    fn cds_del_mismatch() {
        assert_w3016("NM_000088.3:c.100_110delAAAATTTGCC");
    }
}

// =============================================================================
// SECTION 2 — Consistent cases (no warning)
// =============================================================================

mod consistent_no_w3016 {
    use super::*;

    /// 10 positions, 10 bases.
    #[test]
    fn del_consistent() {
        assert_no_w3016("NC_000001.11:g.100_109delAAAATTTGCC");
    }

    /// 5 positions, 5 bases.
    #[test]
    fn inv_consistent() {
        assert_no_w3016("NC_000001.11:g.100_104invTAGCA");
    }

    /// 3 positions, 3 bases.
    #[test]
    fn dup_consistent() {
        assert_no_w3016("NC_000001.11:g.100_102dupATG");
    }
}

// =============================================================================
// SECTION 3 — Skipped cases
// =============================================================================

mod skipped {
    use super::*;

    /// del without explicit ref seq — nothing to compare.
    #[test]
    fn del_without_ref_skipped() {
        assert_no_w3016("NC_000001.11:g.100_110del");
    }

    /// dup without explicit ref seq.
    #[test]
    fn dup_without_ref_skipped() {
        assert_no_w3016("NC_000001.11:g.100_102dup");
    }

    /// Insertion — ref seq is the inserted seq, not the spanning ref.
    /// The position range is the *anchor pair*, not a span; len-rule
    /// doesn't apply.
    #[test]
    fn ins_skipped() {
        assert_no_w3016("NC_000001.11:g.100_101insATG");
    }

    /// Substitution — no range, no comparison.
    #[test]
    fn substitution_skipped() {
        assert_no_w3016("NC_000001.11:g.100A>G");
    }

    /// Offset-bearing endpoints — provider needed to compute range
    /// length on intronic offsets; detector bails.
    #[test]
    fn offset_bearing_skipped() {
        assert_no_w3016("NM_000088.3:c.100+5_100+10del");
    }

    /// `delins<ins>` without explicit `del<ref>` — comparing the inserted
    /// length against the range length is the wrong rule; detector bails.
    #[test]
    fn delins_implicit_del_skipped() {
        assert_no_w3016("NC_000001.11:g.100_102delinsTTCC");
    }

    /// Legacy `del<ref>ins<alt>` form (matched as `EditKind::Del` because
    /// `delins` requires the keyword to come first). The ref-seq scan must
    /// stop at the `ins` boundary so that valid inputs aren't mis-flagged
    /// as W3016. Regression for CodeRabbit feedback on PR #272.
    #[test]
    fn del_ref_ins_alt_explicit_form_not_flagged() {
        assert_no_w3016("NC_000001.11:g.100_102delATGinsT");
    }
}

// =============================================================================
// SECTION 4 — Strict mode rejects
// =============================================================================

mod strict_rejects {
    use super::*;

    #[test]
    fn strict_config_del_mismatch_rejected() {
        let err =
            parse_hgvs_with_config("NC_000001.11:g.100_110delAAAATTTGCC", ErrorConfig::strict())
                .expect_err("expected strict mode to reject length mismatch");
        assert_eq!(
            err.code(),
            Some(ErrorCode::InvalidEdit),
            "expected E1004 InvalidEdit, got: {err}"
        );
    }

    #[test]
    fn strict_config_inv_mismatch_rejected() {
        let err = parse_hgvs_with_config("NC_000001.11:g.100_105invTAGCA", ErrorConfig::strict())
            .expect_err("expected strict mode to reject length mismatch");
        assert_eq!(
            err.code(),
            Some(ErrorCode::InvalidEdit),
            "expected E1004 InvalidEdit, got: {err}"
        );
    }
}

// =============================================================================
// SECTION 5 — Compound bracket allele coverage gap
// =============================================================================
//
// Inside `accession:g.[edit1;edit2]`, member edits inherit the coord
// marker from the outer prefix; there's no per-member `g.` for the
// detector to anchor on. The current implementation does not fire on
// bracketed members. Pinned as a documented limitation.

mod compound_bracket_gap {
    use super::*;

    #[test]
    fn cis_allele_member_mismatch_not_detected() {
        // 11 positions, 10 bases on the first member — but the detector
        // doesn't fire because the coord marker is outside the bracket.
        assert_no_w3016("NC_000001.11:g.[100_110delAAAATTTGCC;200A>G]");
    }
}
