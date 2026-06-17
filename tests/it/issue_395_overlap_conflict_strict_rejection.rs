//! Issue #395 item 6 — `W5002 OverlapConflict` promoted to error in
//! strict mode.
//!
//! `src/error_handling/registry.rs:1190` declared W5002's
//! `mode_behavior` as `ModeBehavior::always_warn_if_not_rejected()`
//! (Strict→Reject, Lenient→WarnAccept, Silent→WarnAccept). But the
//! emit site at `src/normalize/overlap.rs:88` pushes the warning
//! unconditionally, bypassing the policy table. Result: strict mode
//! still emitted the warning and accepted the variant.
//!
//! This file pins:
//!   - Strict mode rejects W5002 with `FerroError::InvalidCoordinates`
//!     citing `OverlapConflictingEdits / W5002`.
//!   - Lenient mode preserves the input and emits the warning.
//!   - Silent mode still emits the warning (the emit site at
//!     `overlap.rs:88` is unconditional); only strict-mode promotion
//!     to error was the gap. Callers using `normalize_with_diagnostics`
//!     get the warning in any mode, while `normalize` only differs
//!     for strict (Err) vs lenient/silent (Ok with warning attached
//!     to the dropped vec).
//!
//! # Spec basis
//!
//! HGVS spec does not define a canonical form for two cis-allele edits
//! sharing identical reference bounds — the variant is ambiguous. Per
//! the registry: "ferro preserves the input verbatim and emits this
//! warning". Strict mode promotes to a typed error so callers can flag
//! the input.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

fn provider() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000001.11", "A".repeat(200));
    p
}

#[test]
fn strict_mode_rejects_coincident_bounds_with_w5002_error() {
    // Two substitutions at the same g.100 position — coincident bounds.
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::strict());
    let variant = parse_hgvs("NC_000001.11:g.[100A>C;100A>G]").expect("parse");
    let err = normalizer
        .normalize(&variant)
        .expect_err("strict mode must reject coincident-bounds cis edits");
    match err {
        FerroError::InvalidCoordinates { msg } => {
            assert!(
                msg.contains("W5002") || msg.contains("OverlapConflictingEdits"),
                "strict-mode error message must reference W5002 / OverlapConflictingEdits; got: {msg}",
            );
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}

#[test]
fn lenient_mode_emits_w5002_warning_and_preserves_input() {
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::lenient());
    let variant = parse_hgvs("NC_000001.11:g.[100A>C;100A>G]").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("lenient mode must accept coincident-bounds cis edits");
    let has_warning = result
        .warnings
        .iter()
        .any(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }));
    assert!(
        has_warning,
        "lenient mode must emit OverlapConflict (W5002) warning; got: {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
    // "Preserves input" contract: the normalized variant equals the
    // parsed input verbatim. Per the registry text and the
    // `OverlapConflictingEdits` docstring, ferro does NOT canonicalize
    // coincident-bounds cis edits — the input is returned unchanged.
    assert_eq!(
        result.result, variant,
        "lenient mode must preserve the parsed input verbatim; got: {}",
        result.result,
    );
}

#[test]
fn strict_mode_accepts_non_coincident_cis_edits() {
    // Negative control: cis edits at DIFFERENT positions must not
    // trigger the overlap rejection.
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::strict());
    let variant = parse_hgvs("NC_000001.11:g.[100A>C;101A>G]").expect("parse");
    normalizer.normalize(&variant).expect(
        "strict mode must accept non-coincident cis edits; the W5002 trigger is coincident-bounds only",
    );
}

#[test]
fn silent_mode_still_emits_w5002_warning() {
    // The emit site at `src/normalize/overlap.rs:88` is unconditional —
    // silent mode does not suppress the warning, only changes
    // strict-mode behavior from Reject to WarnAccept. Callers using
    // `normalize_with_diagnostics` get the warning regardless of mode;
    // `normalize` only differs by whether the warning is promoted
    // (strict → Err) or accepted (lenient/silent → Ok).
    let normalizer = Normalizer::with_config(provider(), NormalizeConfig::silent());
    let variant = parse_hgvs("NC_000001.11:g.[100A>C;100A>G]").expect("parse");
    let result = normalizer
        .normalize_with_diagnostics(&variant)
        .expect("silent mode must accept coincident-bounds cis edits");
    let has_warning = result
        .warnings
        .iter()
        .any(|w| matches!(w, NormalizationWarning::OverlapConflict { .. }));
    assert!(
        has_warning,
        "silent mode still emits the OverlapConflict warning at the detection site; \
         only strict-mode promotion is gated by mode. Got: {:?}",
        result.warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
    // "Preserves input" contract: silent mode, like lenient, returns
    // the parsed input verbatim — ferro does NOT canonicalize
    // coincident-bounds cis edits in any non-strict mode.
    assert_eq!(
        result.result, variant,
        "silent mode must preserve the parsed input verbatim; got: {}",
        result.result,
    );
}
