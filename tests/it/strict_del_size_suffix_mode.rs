//! Strict-mode rejection of the pre-2016 `del<N>` size-count suffix (#447).
//!
//! Per HGVS v21 `checklist.md:49`:
//!
//! > `g.123del3` is invalid: pre-2016 size-suffix form. Must be
//! > `g.123_125del`.
//!
//! ferro's default `parse_hgvs` is **lenient** — it accepts the legacy
//! form, round-trips it (preserving the `<N>` as `NaEdit::Deletion.length`),
//! and emits W3011 (DelSizeSuffix) as a non-rewriting warning. ferro cannot
//! synthesise the end position safely because offset/intronic semantics
//! defeat naive `start + length - 1` arithmetic.
//!
//! For callers that want strict spec conformance, override
//! `ErrorType::DelSizeSuffix` to `Reject` on the config:
//!
//! ```rust,no_run
//! use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
//! use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
//!
//! let strict_del_size = ErrorConfig::lenient()
//!     .with_override(ErrorType::DelSizeSuffix, ErrorOverride::Reject);
//! parse_hgvs_with_config("NG_012232.1:g.123del6", strict_del_size)
//!     .unwrap_err();  // rejected per checklist.md:49
//! ```
//!
//! This file pins both the strict-mode rejection and the default lenient
//! pass-through so neither side regresses silently. Closes #447.
//!
//! See follow-up issue (separately filed) for the related `dup<N>` /
//! `dup<seq>` / `del<seq>` soft-prohibition forms per
//! `duplication.md:35,50,140-143`.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;

fn strict_del_size() -> ErrorConfig {
    ErrorConfig::lenient().with_override(ErrorType::DelSizeSuffix, ErrorOverride::Reject)
}

fn lenient() -> ErrorConfig {
    ErrorConfig::lenient()
}

fn assert_strict_rejects(input: &str, note: &str) {
    let result = parse_hgvs_with_config(input, strict_del_size());
    assert!(
        result.is_err(),
        "strict del-size mode must reject {input:?} ({note}); got: {:?}",
        result.map(|o| o.result.to_string())
    );
    let msg = result.unwrap_err().to_string();
    // The preprocessor's canonical rejection text always contains
    // "size-count suffix" — assert only on that diagnostic phrase rather
    // than test-input-specific literals like "del6"/"del3".
    assert!(
        msg.contains("size-count suffix"),
        "rejection should cite the size-count form; got: {msg}"
    );
}

fn assert_lenient_accepts(input: &str, expected_display: &str) {
    let out = parse_hgvs_with_config(input, lenient())
        .unwrap_or_else(|e| panic!("lenient must accept {input:?}: {e}"));
    assert_eq!(out.result.to_string(), expected_display);
}

/// Like [`assert_lenient_accepts`], but for the legacy `del<N>` size-count
/// form: lenient mode must accept AND surface the W3011/DelSizeSuffix
/// warning (non-rewriting). Acceptance alone would miss a regression
/// where `del<N>` is silently accepted with the flag dropped.
fn assert_lenient_accepts_with_size_warning(input: &str, expected_display: &str) {
    let out = parse_hgvs_with_config(input, lenient())
        .unwrap_or_else(|e| panic!("lenient must accept {input:?}: {e}"));
    assert_eq!(out.result.to_string(), expected_display);
    assert!(
        out.warnings
            .iter()
            .any(|w| w.error_type == ErrorType::DelSizeSuffix),
        "lenient parse of {input:?} must emit a DelSizeSuffix (W3011) warning; got {:?}",
        out.warnings
            .iter()
            .map(|w| w.error_type)
            .collect::<Vec<_>>(),
    );
}

// =====================================================================
// Strict mode REJECTS `del<N>` on every coordinate axis
// =====================================================================

#[test]
fn strict_rejects_genomic_del_size_suffix() {
    // `checklist.md:49` canonical example.
    assert_strict_rejects("NG_012232.1:g.123del6", "g.<pos>del<size>");
}

#[test]
fn strict_rejects_cds_del_size_suffix() {
    assert_strict_rejects("NM_004006.2:c.76del3", "c.<pos>del<size>");
}

#[test]
fn strict_rejects_noncoding_tx_del_size_suffix() {
    assert_strict_rejects("NR_002196.2:n.123del6", "n.<pos>del<size>");
}

#[test]
fn strict_rejects_rna_del_size_suffix() {
    assert_strict_rejects("NR_002196.2:r.123del6", "r.<pos>del<size>");
}

#[test]
fn strict_rejects_mito_del_size_suffix() {
    assert_strict_rejects("NC_012920.1:m.3243del4", "m.<pos>del<size>");
}

#[test]
fn strict_rejects_large_size_count() {
    // Even legitimate-looking size counts are rejected — the form itself
    // is invalid, not the magnitude.
    assert_strict_rejects("NG_012232.1:g.123del1000", "large size count");
}

// =====================================================================
// Default lenient mode accepts `del<N>` (existing intentional behavior —
// preserved for interop with real-world ClinVar / pre-2016 submissions)
// =====================================================================

#[test]
fn lenient_accepts_genomic_del_size_suffix() {
    assert_lenient_accepts_with_size_warning("NG_012232.1:g.123del6", "NG_012232.1:g.123del6");
}

#[test]
fn lenient_accepts_cds_del_size_suffix() {
    assert_lenient_accepts_with_size_warning("NM_004006.2:c.76del3", "NM_004006.2:c.76del3");
}

#[test]
fn lenient_accepts_mito_del_size_suffix() {
    assert_lenient_accepts_with_size_warning("NC_012920.1:m.3243del4", "NC_012920.1:m.3243del4");
}

// =====================================================================
// Canonical range-form deletions pass cleanly in BOTH modes
// =====================================================================

#[test]
fn strict_accepts_canonical_genomic_range_del() {
    let out = parse_hgvs_with_config("NG_012232.1:g.123_128del", strict_del_size())
        .expect("canonical range form must pass strict mode");
    assert_eq!(out.result.to_string(), "NG_012232.1:g.123_128del");
}

#[test]
fn strict_accepts_canonical_cds_range_del() {
    let out = parse_hgvs_with_config("NM_004006.2:c.76_78del", strict_del_size())
        .expect("canonical c. range form must pass strict mode");
    assert_eq!(out.result.to_string(), "NM_004006.2:c.76_78del");
}

#[test]
fn strict_accepts_single_position_del() {
    // Single-residue deletion needs no size suffix and no range — it's
    // canonical as-is.
    let out = parse_hgvs_with_config("NG_012232.1:g.123del", strict_del_size())
        .expect("single-residue del must pass strict mode");
    assert_eq!(out.result.to_string(), "NG_012232.1:g.123del");
}

#[test]
fn lenient_accepts_canonical_range_del() {
    assert_lenient_accepts("NG_012232.1:g.123_128del", "NG_012232.1:g.123_128del");
}

// =====================================================================
// Cross-cutting: the diagnostic carries the proper hint
// =====================================================================

#[test]
fn strict_rejection_diagnostic_includes_canonical_hint() {
    let result = parse_hgvs_with_config("NG_012232.1:g.123del6", strict_del_size());
    let err = result.unwrap_err();
    let formatted = format!("{err:?}");
    assert!(
        formatted.contains("g.<start>_<end>del") || formatted.contains("start>_<end>"),
        "diagnostic should hint at canonical range form; got: {formatted}"
    );
}
