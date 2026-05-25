//! Strict-mode rejection + lenient-mode warn/correct for the soft-prohibition
//! forms `dup<N>`, `dup<seq>`, `del<seq>` (#460).
//!
//! Per HGVS v21 spec:
//! - `dup<N>` (e.g. `c.20_21dup2`): non-canonical per `duplication.md:140-143`
//!   Q&A "No, a duplication of more than one nucleotide should give the
//!   position of the first and last nucleotide duplicated". W3023 (DupSizeSuffix)
//!   uses `warn_accept` — strict rejects, lenient warns without rewrite,
//!   silent accepts.
//! - `dup<seq>` (e.g. `c.20_23dupTAGA`): non-canonical per
//!   `duplication.md:35-36` "the recommendation is not to ... chances to
//!   make an error increases". W3024 (DupExplicitSeq) uses
//!   `standard_correctable` — strict rejects, lenient rewrites + warns,
//!   silent rewrites silently.
//! - `del<seq>` (e.g. `g.33344590_33344592delGAT`): same register per
//!   `deletion.md:30-31`. W3025 (DelExplicitSeq), same mapping as W3024.
//!
//! ferro's existing W3011 (DelSizeSuffix) is pinned by
//! `tests/strict_del_size_suffix_mode.rs` (PR #459). This file pins the
//! three new sibling codes. Closes #460.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;

fn strict() -> ErrorConfig {
    ErrorConfig::strict()
}

fn lenient() -> ErrorConfig {
    ErrorConfig::lenient()
}

fn assert_strict_rejects(input: &str, expected_in_msg: &str) {
    let result = parse_hgvs_with_config(input, strict());
    assert!(
        result.is_err(),
        "strict must reject {input:?}; got: {:?}",
        result.map(|r| r.result.to_string())
    );
    let msg = result.unwrap_err().to_string();
    assert!(
        msg.to_lowercase().contains(&expected_in_msg.to_lowercase()),
        "strict rejection of {input:?} should mention {expected_in_msg:?}; got: {msg}"
    );
}

fn assert_lenient_emits(input: &str, expected_display: &str, kind: ErrorType) {
    let result = parse_hgvs_with_config(input, lenient())
        .unwrap_or_else(|e| panic!("lenient must accept {input:?}: {e}"));
    assert_eq!(
        result.result.to_string(),
        expected_display,
        "lenient Display mismatch for {input:?}",
    );
    assert!(
        result.warnings.iter().any(|w| w.error_type == kind),
        "lenient {input:?} should emit {kind:?}; got: {:?}",
        result.warnings,
    );
}

// =====================================================================
// W3023 DupSizeSuffix — warn_accept (no rewrite)
// =====================================================================

#[test]
fn strict_rejects_dup_size_suffix_genomic() {
    assert_strict_rejects("NC_000023.11:g.123dup6", "size-count suffix");
}

#[test]
fn strict_rejects_dup_size_suffix_cds() {
    assert_strict_rejects("NM_004006.2:c.20_21dup2", "size-count suffix");
}

#[test]
fn strict_rejects_dup_size_suffix_rna() {
    assert_strict_rejects("NM_004006.2:r.20_21dup2", "size-count suffix");
}

#[test]
fn lenient_warns_but_preserves_dup_size_suffix() {
    // warn_accept: input is preserved through to Display.
    let result = parse_hgvs_with_config("NM_004006.2:c.20_21dup2", lenient())
        .expect("lenient must accept dup<N>");
    assert!(result
        .warnings
        .iter()
        .any(|w| w.error_type == ErrorType::DupSizeSuffix));
    // Display preserves the legacy form (NaEdit::Duplication.length is populated).
    assert_eq!(result.result.to_string(), "NM_004006.2:c.20_21dup2");
}

// =====================================================================
// W3024 DupExplicitSeq — standard_correctable (rewrite drops seq)
// =====================================================================

#[test]
fn strict_rejects_dup_explicit_seq_cds() {
    assert_strict_rejects("NM_004006.2:c.20_23dupTAGA", "explicit sequence");
}

#[test]
fn strict_rejects_dup_explicit_seq_single_base() {
    assert_strict_rejects("NM_004006.2:c.20dupT", "explicit sequence");
}

#[test]
fn strict_rejects_dup_explicit_seq_rna() {
    assert_strict_rejects("NM_004006.2:r.6_8dupugc", "explicit sequence");
}

#[test]
fn lenient_rewrites_dup_explicit_seq_dna() {
    assert_lenient_emits(
        "NM_004006.2:c.20_23dupTAGA",
        "NM_004006.2:c.20_23dup",
        ErrorType::DupExplicitSeq,
    );
}

#[test]
fn lenient_rewrites_dup_explicit_seq_rna() {
    assert_lenient_emits(
        "NM_004006.2:r.6_8dupugc",
        "NM_004006.2:r.6_8dup",
        ErrorType::DupExplicitSeq,
    );
}

// =====================================================================
// W3025 DelExplicitSeq — standard_correctable (rewrite drops seq)
// =====================================================================

#[test]
fn strict_rejects_del_explicit_seq_genomic() {
    assert_strict_rejects(
        "NC_000023.11:g.33344590_33344592delGAT",
        "explicit sequence",
    );
}

#[test]
fn strict_rejects_del_explicit_seq_single_base() {
    assert_strict_rejects("NC_000023.11:g.33344591delA", "explicit sequence");
}

#[test]
fn strict_rejects_del_explicit_seq_rna() {
    assert_strict_rejects("NM_004006.2:r.6_8deluug", "explicit sequence");
}

#[test]
fn lenient_rewrites_del_explicit_seq_dna() {
    assert_lenient_emits(
        "NC_000023.11:g.33344590_33344592delGAT",
        "NC_000023.11:g.33344590_33344592del",
        ErrorType::DelExplicitSeq,
    );
}

#[test]
fn lenient_rewrites_del_explicit_seq_rna() {
    assert_lenient_emits(
        "NM_004006.2:r.6_8deluug",
        "NM_004006.2:r.6_8del",
        ErrorType::DelExplicitSeq,
    );
}

#[test]
fn strict_does_not_reject_delins() {
    // delins<seq> is canonical for delins; must NOT fire W3025.
    let result = parse_hgvs_with_config("NC_000001.11:g.123delinsATG", strict())
        .expect("strict must accept canonical delins");
    assert_eq!(result.result.to_string(), "NC_000001.11:g.123delinsATG");
}

// =====================================================================
// Canonical inputs pass cleanly in both modes
// =====================================================================

#[test]
fn strict_accepts_canonical_dup() {
    let out = parse_hgvs_with_config("NM_004006.2:c.20_21dup", strict())
        .expect("canonical c.20_21dup must pass strict");
    assert_eq!(out.result.to_string(), "NM_004006.2:c.20_21dup");
}

#[test]
fn strict_accepts_canonical_del() {
    let out = parse_hgvs_with_config("NC_000023.11:g.33344591del", strict())
        .expect("canonical g.33344591del must pass strict");
    assert_eq!(out.result.to_string(), "NC_000023.11:g.33344591del");
}

#[test]
fn lenient_accepts_canonical_silent() {
    let out = parse_hgvs_with_config("NM_004006.2:c.20_23dup", lenient())
        .expect("canonical c.20_23dup must pass lenient silently");
    assert!(
        out.warnings
            .iter()
            .all(|w| w.error_type != ErrorType::DupExplicitSeq
                && w.error_type != ErrorType::DupSizeSuffix),
        "canonical must not emit W3023/W3024",
    );
}

// =====================================================================
// Diagnostic surface check
// =====================================================================

#[test]
fn diagnostic_includes_canonical_hint() {
    let result = parse_hgvs_with_config("NC_000023.11:g.33344591delA", strict());
    let err = result.unwrap_err();
    let formatted = format!("{err:?}");
    assert!(
        formatted.contains("Drop the redundant sequence")
            || formatted.contains("g.<start>_<end>del"),
        "diagnostic should hint at canonical form; got: {formatted}"
    );
}

// =====================================================================
// User can opt into strict-only behavior via per-error override
// =====================================================================

#[test]
fn lenient_with_dup_explicit_seq_override_rejects() {
    let cfg =
        ErrorConfig::lenient().with_override(ErrorType::DupExplicitSeq, ErrorOverride::Reject);
    let result = parse_hgvs_with_config("NM_004006.2:c.20_23dupTAGA", cfg);
    assert!(
        result.is_err(),
        "override should reject {:?}",
        result.map(|r| r.result.to_string())
    );
}
