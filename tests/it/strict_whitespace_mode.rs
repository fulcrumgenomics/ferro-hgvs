//! Strict-whitespace parse mode (#449).
//!
//! Per HGVS v21 `general.md:96`: "spaces are not permitted in any HGVS
//! description". ferro's default `parse_hgvs` is lenient — it trims
//! outer whitespace and strips spaces inside allele brackets for interop
//! with real-world submissions (CDS-Var, ClinVar, VCF INFO field
//! manually-curated entries).
//!
//! For callers that want strict spec conformance, override
//! `ErrorType::ExtraWhitespace` to `Reject` on the config:
//!
//! ```rust,no_run
//! use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
//! use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
//!
//! let strict_ws = ErrorConfig::lenient()
//!     .with_override(ErrorType::ExtraWhitespace, ErrorOverride::Reject);
//! parse_hgvs_with_config("  NM_004006.2:c.76A>G", strict_ws)
//!     .unwrap_err();  // rejected per general.md:96
//! ```
//!
//! This file pins both the strict-mode rejection and the default
//! lenient pass-through so neither side regresses silently. Closes #449.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorOverride, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;

fn strict_ws() -> ErrorConfig {
    ErrorConfig::lenient().with_override(ErrorType::ExtraWhitespace, ErrorOverride::Reject)
}

fn lenient() -> ErrorConfig {
    ErrorConfig::lenient()
}

fn assert_strict_rejects(input: &str, note: &str) {
    let result = parse_hgvs_with_config(input, strict_ws());
    assert!(
        result.is_err(),
        "strict whitespace mode must reject {input:?} ({note}); got: {:?}",
        result.map(|o| o.result.to_string())
    );
    let msg = result.unwrap_err().to_string();
    assert!(
        msg.contains("whitespace") || msg.contains("Whitespace"),
        "rejection should cite whitespace; got: {msg}"
    );
}

fn assert_lenient_accepts(input: &str, expected_display: &str) {
    let out = parse_hgvs_with_config(input, lenient())
        .unwrap_or_else(|e| panic!("lenient must accept {input:?}: {e}"));
    assert_eq!(out.result.to_string(), expected_display);
}

// =====================================================================
// Strict mode REJECTS whitespace anywhere
// =====================================================================

#[test]
fn strict_rejects_leading_whitespace_on_identifier() {
    assert_strict_rejects("  NM_004006.2:c.76A>G", "leading whitespace");
}

#[test]
fn strict_rejects_trailing_whitespace_on_identifier() {
    assert_strict_rejects("NM_004006.2:c.76A>G  ", "trailing whitespace");
}

#[test]
fn strict_rejects_whitespace_inside_allele_bracket_dna() {
    assert_strict_rejects(
        "NM_004006.2:c.[76A>C; 80T>G]",
        "space after `;` inside allele bracket",
    );
}

#[test]
fn strict_rejects_whitespace_inside_allele_bracket_rna() {
    assert_strict_rejects("NM_004006.2:r.[76a>c; 80u>g]", "RNA allele-bracket space");
}

#[test]
fn strict_rejects_whitespace_inside_variant_body() {
    assert_strict_rejects("NM_004006.2:c. 76 A > G", "spaces around tokens");
}

#[test]
fn strict_rejects_tab_whitespace() {
    assert_strict_rejects("NM_004006.2:c.76A>G\t", "trailing tab");
}

// =====================================================================
// Default lenient mode trims outer / strips bracket whitespace
// (existing intentional ferro behavior — preserved by this change)
// =====================================================================

#[test]
fn lenient_trims_leading_whitespace() {
    assert_lenient_accepts("  NM_004006.2:c.76A>G", "NM_004006.2:c.76A>G");
}

#[test]
fn lenient_trims_trailing_whitespace() {
    assert_lenient_accepts("NM_004006.2:c.76A>G  ", "NM_004006.2:c.76A>G");
}

#[test]
fn lenient_strips_whitespace_inside_allele_bracket() {
    assert_lenient_accepts(
        "NM_004006.2:c.[76A>C; 80T>G]",
        "NM_004006.2:c.[76A>C;80T>G]",
    );
}

// =====================================================================
// Spec-canonical input round-trips in both modes
// =====================================================================

#[test]
fn strict_accepts_canonical_no_whitespace() {
    let out = parse_hgvs_with_config("NM_004006.2:c.76A>G", strict_ws())
        .expect("canonical input must pass strict mode");
    assert_eq!(out.result.to_string(), "NM_004006.2:c.76A>G");
}

#[test]
fn strict_accepts_canonical_allele() {
    let out = parse_hgvs_with_config("NM_004006.2:c.[76A>C;80T>G]", strict_ws())
        .expect("canonical allele input must pass strict mode");
    assert_eq!(out.result.to_string(), "NM_004006.2:c.[76A>C;80T>G]");
}
