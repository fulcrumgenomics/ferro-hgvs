//! Legacy stop-codon glyphs in protein frameshift descriptions (W3009 /
//! W3010 ModeBehavior pins).
//!
//! HGVS v21 endorses both `Ter` (three-letter) and `*` (one-letter) as
//! valid stop-codon glyphs (`background/standards.md:213`,
//! `protein/frameshift.md` examples). The legacy pre-v15.11 spelling `X`
//! is explicitly **not** valid: `checklist.md:63` says *"`Ter` or `*`
//! should be used to indicate a translation stop codon; the `X` should
//! not be used."*
//!
//! ferro's mode-aware error handling distinguishes these correctly:
//! - `Ter` and `*` (with digits OR `?`) are parser-native and round-trip
//!   silently in every mode.
//! - `X` (with digits OR `?`) triggers **W3010 DeprecatedFrameshiftX**:
//!   strict rejects, lenient rewrites → `Ter` + warning, silent rewrites
//!   silently.
//! - `*` ALSO triggers **W3009 DeprecatedFrameshiftStar** at the
//!   preprocessor level (ferro is stricter than the spec here, treating
//!   `*` as deprecated and canonicalizing to `Ter`). Same mode mapping.
//!
//! This file pins all four lanes — `fsXN`, `fsX?`, `fs*N`, `fs*?` — in
//! both lenient (warn-and-correct) and strict (reject) modes, plus the
//! canonical `fsTerN` / `fsTer?` round-trip. Closes the gap surfaced
//! while investigating #81 D1: previously `fsX?` was rejected outright
//! by the parser (corrector required digits after `X`), inconsistent
//! with the digit-bearing `fsXN` lenient path.

use ferro_hgvs::error_handling::{ErrorConfig, ErrorType};
use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;

const ACC: &str = "NP_000079.2";

fn lenient() -> ErrorConfig {
    ErrorConfig::lenient()
}

fn strict() -> ErrorConfig {
    ErrorConfig::strict()
}

fn lenient_rewrites_to(input: &str, expected_display: &str, expected_warning: ErrorType) {
    let result = parse_hgvs_with_config(input, lenient())
        .unwrap_or_else(|e| panic!("lenient must accept {input:?}: {e}"));
    assert_eq!(
        result.result.to_string(),
        expected_display,
        "lenient Display mismatch for {input:?}",
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.error_type == expected_warning),
        "lenient {input:?} should emit {expected_warning:?}; got: {:?}",
        result.warnings,
    );
}

fn strict_rejects(input: &str, expected_in_msg: &str) {
    let result = parse_hgvs_with_config(input, strict());
    assert!(
        result.is_err(),
        "strict must reject {input:?}; got: {:?}",
        result.map(|r| r.result.to_string()),
    );
    let msg = result.unwrap_err().to_string();
    assert!(
        msg.contains(expected_in_msg),
        "strict rejection of {input:?} should mention {expected_in_msg:?}; got: {msg}"
    );
}

fn lenient_silent_roundtrip(input: &str) {
    let result = parse_hgvs_with_config(input, lenient())
        .unwrap_or_else(|e| panic!("lenient must accept {input:?}: {e}"));
    assert_eq!(result.result.to_string(), input, "round-trip mismatch");
    assert!(
        result
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::DeprecatedFrameshiftStar
                && w.error_type != ErrorType::DeprecatedFrameshiftX),
        "canonical {input:?} should not emit W3009/W3010; got: {:?}",
        result.warnings,
    );
}

// =====================================================================
// W3010 — fsX (legacy 'X' spelling, checklist.md:63 'should not be used')
// =====================================================================

#[test]
fn lenient_rewrites_fs_x_n_with_digits() {
    lenient_rewrites_to(
        &format!("{ACC}:p.Arg97fsX23"),
        &format!("{ACC}:p.Arg97fsTer23"),
        ErrorType::DeprecatedFrameshiftX,
    );
}

#[test]
fn lenient_rewrites_fs_x_question_uncertain() {
    // Corrector handles the `?` (uncertain stop position) lane, not just
    // digits. The deprecated `X?` rewrites to the canonical unknown-stop
    // marker `fsTer?` (preserved, no longer collapsed to bare `fs`).
    lenient_rewrites_to(
        &format!("{ACC}:p.Arg97fsX?"),
        &format!("{ACC}:p.Arg97fsTer?"),
        ErrorType::DeprecatedFrameshiftX,
    );
}

#[test]
fn lenient_rewrites_profs_x_n_with_new_aa() {
    lenient_rewrites_to(
        &format!("{ACC}:p.Arg97ProfsX23"),
        &format!("{ACC}:p.Arg97ProfsTer23"),
        ErrorType::DeprecatedFrameshiftX,
    );
}

#[test]
fn strict_rejects_fs_x_n_per_checklist_63() {
    strict_rejects(&format!("{ACC}:p.Arg97fsX23"), "Deprecated");
}

#[test]
fn strict_rejects_fs_x_question_per_checklist_63() {
    strict_rejects(&format!("{ACC}:p.Arg97fsX?"), "Deprecated");
}

// =====================================================================
// W3009 — fs* (one-letter Star spelling, spec-valid but ferro deprecates)
// =====================================================================

#[test]
fn lenient_rewrites_fs_star_n_with_digits() {
    lenient_rewrites_to(
        &format!("{ACC}:p.Arg97fs*23"),
        &format!("{ACC}:p.Arg97fsTer23"),
        ErrorType::DeprecatedFrameshiftStar,
    );
}

#[test]
fn lenient_rewrites_fs_star_question_uncertain() {
    // `fs*?` goes through the corrector (same mode mapping as `fs*N`); the
    // deprecated `*?` rewrites to the canonical unknown-stop marker `fsTer?`
    // (preserved, no longer collapsed to bare `fs`).
    lenient_rewrites_to(
        &format!("{ACC}:p.Arg97fs*?"),
        &format!("{ACC}:p.Arg97fsTer?"),
        ErrorType::DeprecatedFrameshiftStar,
    );
}

#[test]
fn strict_rejects_fs_star_n() {
    strict_rejects(&format!("{ACC}:p.Arg97fs*23"), "Deprecated");
}

#[test]
fn strict_rejects_fs_star_question() {
    strict_rejects(&format!("{ACC}:p.Arg97fs*?"), "Deprecated");
}

// =====================================================================
// Canonical `Ter` forms: parser-native, silent in BOTH modes
// =====================================================================

#[test]
fn lenient_canonical_fs_ter_n_silent_roundtrip() {
    lenient_silent_roundtrip(&format!("{ACC}:p.Arg97ProfsTer23"));
}

#[test]
fn lenient_canonical_fs_ter_question_silent_roundtrip() {
    // Uncertain-stop with canonical Ter glyph: `fsTer?` round-trips unchanged
    // (the explicit unknown-stop marker is preserved; the short form `fs`
    // remains an accepted equivalent input per frameshift.md `p.Ile327Argfs*?`).
    let input = format!("{ACC}:p.Arg97fsTer?");
    let result = parse_hgvs_with_config(&input, lenient()).expect("parse");
    let display = result.result.to_string();
    assert!(
        display == input || display == format!("{ACC}:p.Arg97fs"),
        "fsTer? should round-trip or canonicalize to short form; got {display}",
    );
    assert!(
        result.warnings.is_empty()
            || result
                .warnings
                .iter()
                .all(|w| w.error_type != ErrorType::DeprecatedFrameshiftStar
                    && w.error_type != ErrorType::DeprecatedFrameshiftX),
        "canonical fsTer? should not emit W3009/W3010",
    );
}

#[test]
fn lenient_short_form_fs_silent_roundtrip() {
    lenient_silent_roundtrip(&format!("{ACC}:p.Arg97fs"));
}

#[test]
fn strict_accepts_canonical_fs_ter_n() {
    let result =
        parse_hgvs_with_config(&format!("{ACC}:p.Arg97ProfsTer23"), strict()).expect("parse");
    assert_eq!(
        result.result.to_string(),
        format!("{ACC}:p.Arg97ProfsTer23"),
    );
}

#[test]
fn strict_accepts_short_form_fs() {
    let result = parse_hgvs_with_config(&format!("{ACC}:p.Arg97fs"), strict()).expect("parse");
    assert_eq!(result.result.to_string(), format!("{ACC}:p.Arg97fs"));
}

#[test]
fn strict_accepts_canonical_fs_ter_question() {
    // `fsTer?` is the uncertain canonical lane this PR's W3009/W3010
    // corrector touches: `Ter` (not `*`/`X`) with an unknown new-stop
    // offset. Strict mode must accept it without emitting the
    // deprecated-glyph warnings — it is already canonical.
    let input = format!("{ACC}:p.Arg97fsTer?");
    let result = parse_hgvs_with_config(&input, strict()).expect("parse");
    let display = result.result.to_string();
    assert!(
        display == input || display == format!("{ACC}:p.Arg97fs"),
        "strict fsTer? should round-trip or canonicalize to the short form; got {display}",
    );
    assert!(
        result
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::DeprecatedFrameshiftStar
                && w.error_type != ErrorType::DeprecatedFrameshiftX),
        "canonical fsTer? must not emit W3009/W3010; got {:?}",
        result
            .warnings
            .iter()
            .map(|w| w.error_type)
            .collect::<Vec<_>>(),
    );
}
