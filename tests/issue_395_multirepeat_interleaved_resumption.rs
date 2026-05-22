//! Issue #395 item 1 — extend `validate_multirepeat_tract` to validate
//! trailing `Exact` units after a non-`Exact` middle unit, by
//! right-anchoring a suffix walk in addition to the existing
//! left-anchored prefix walk.
//!
//! # Current behavior (pre-PR)
//!
//! `src/normalize/validate.rs::validate_multirepeat_tract` walks units
//! left-to-right and `break`s at the first non-`Exact` count. For
//! `[Exact(3), Range(1,5), Exact(2)]`, the trailing `Exact(2)` is
//! never validated even though its bases at the *right* edge of the
//! reference tract are deterministic.
//!
//! # New behavior
//!
//! - Leading-`Exact` prefix is still validated left-anchored.
//! - **New:** trailing-`Exact` suffix is validated right-anchored
//!   (last `sum(trailing_exact_lens)` bytes of the ref span).
//! - The non-`Exact` middle is still skipped.
//! - When the declared prefix + suffix lengths would together exceed
//!   the reference span, that's a content/length mismatch (the suffix
//!   would overlap the prefix or run off the start of the tract).
//!
//! # Spec basis
//!
//! `assets/hgvs-nomenclature/docs/recommendations/DNA/repeated.md`:
//! every unit's bases are part of the reference tract. Validating
//! *any* deterministic prefix/suffix matches the spec; what's added
//! here is the symmetric right-anchored case.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}
const CORE_START: u64 = PAD.len() as u64 + 1;

fn g_provider(accession: &str, padded_seq: &str) -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(accession, padded_seq);
    p
}

fn normalize_strict(provider: MockProvider, input: &str) -> Result<String, FerroError> {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs(input).expect("parse");
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

fn normalize_lenient(provider: MockProvider, input: &str) -> (String, Vec<NormalizationWarning>) {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs(input).expect("parse");
    let r = normalizer
        .normalize_with_warnings(&variant)
        .expect("lenient normalize");
    (format!("{}", r.result), r.warnings)
}

fn has_refseq_mismatch(warnings: &[NormalizationWarning]) -> bool {
    warnings
        .iter()
        .any(|w| matches!(w, NormalizationWarning::RefSeqMismatch { .. }))
}

fn is_ref_mismatch(err: &FerroError) -> bool {
    matches!(err, FerroError::ReferenceMismatch { .. })
}

// =============================================================================
// Trailing-Exact suffix validation (the new behavior)
// =============================================================================

/// Reference matches both prefix AND suffix: must accept.
///
/// `[CTG[2], TTG[1..3], ACG[2]]` declared over a 15 bp tract whose
/// content is `CTGCTG.TTGTTG.ACGACG` (rendered with the middle's actual
/// count = 2). Prefix CTGCTG (6 bp) matches the start; suffix ACGACG
/// (6 bp) matches the end; the middle 3 bp `TTG` is in the declared
/// range and skipped.
#[test]
fn interleaved_prefix_and_trailing_exact_consistent_accepted_strict() {
    // 2 CTG (6) + 1 TTG (3) + 2 ACG (6) = 15 bp.
    let core = "CTGCTGTTGACGACG";
    assert_eq!(core.len(), 15);
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[2]", CORE_START, end);
    let out =
        normalize_strict(provider, &input).expect("matching prefix + matching suffix must accept");
    // Just sanity that it produced something parseable.
    assert!(out.contains("CTG[2]") && out.contains("ACG[2]"));
}

/// Suffix bases DON'T match the reference: must reject in strict mode
/// (the new validation path catches this).
#[test]
fn trailing_exact_suffix_mismatch_rejected_strict() {
    // Declared: CTG[2] (6) + TTG[1_3] (skip) + ACG[2] (6 expected at right end).
    // Actual end: GGGGGG (not ACGACG) → suffix mismatch.
    let core = "CTGCTGTTGGGGGGG"; // 15 bp: 6 CTG + 3 TTG + 6 G
    assert_eq!(core.len(), 15);
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[2]", CORE_START, end);
    let err = normalize_strict(provider, &input)
        .expect_err("trailing-Exact suffix mismatch must reject in strict mode");
    assert!(
        is_ref_mismatch(&err),
        "expected ReferenceMismatch from trailing-Exact suffix check; got: {err:?}",
    );
}

/// Lenient mode: suffix mismatch emits a `RefSeqMismatch` warning.
#[test]
fn trailing_exact_suffix_mismatch_warns_lenient() {
    let core = "CTGCTGTTGGGGGGG";
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[2]", CORE_START, end);
    let (_, warnings) = normalize_lenient(provider, &input);
    assert!(
        has_refseq_mismatch(&warnings),
        "expected RefSeqMismatch warning in lenient mode; got: {:?}",
        warnings.iter().map(|w| w.code()).collect::<Vec<_>>(),
    );
}

/// Prefix matches but the span is too short to fit both prefix +
/// suffix simultaneously: that's a length mismatch (the suffix would
/// overlap the prefix). Must reject.
#[test]
fn prefix_plus_suffix_overflow_span_rejected_strict() {
    // Declare CTG[2] (6 bp prefix) + TTG[1_3] (skip) + ACG[3] (9 bp
    // suffix) but actual span is only 13 bp — prefix + suffix = 15 bp
    // exceeds the available span (no room for ANY ambiguous middle).
    let core = "CTGCTGTTGGGGGG"; // 14 bp
    assert_eq!(core.len(), 14);
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[3]", CORE_START, end);
    let err = normalize_strict(provider, &input)
        .expect_err("prefix + suffix together exceeding span must reject (overlap)");
    assert!(
        is_ref_mismatch(&err),
        "expected ReferenceMismatch for prefix+suffix overflow; got: {err:?}",
    );
}

// =============================================================================
// Negative controls — preserve existing behavior
// =============================================================================

/// All-Exact regression: full validation still works (prefix walk
/// reaches end, suffix walk has nothing left to check).
#[test]
fn all_exact_consistent_accepted_strict_regression() {
    // 2 CTG (6) + 1 TTG (3) + 11 CTG (33) = 42 bp.
    let core = "CTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[11]", CORE_START, end);
    normalize_strict(provider, &input).expect("all-Exact must still accept");
}

/// First-unit-non-Exact, no trailing Exact: no validation (no
/// deterministic prefix or suffix).
#[test]
fn no_anchored_units_skips_validation() {
    // First unit Range, last unit also Range → no validation anchor on
    // either side. Reference bases are arbitrary; must accept.
    let core = "CTGCTGCTGTTGACGACGACG"; // 21 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!(
        "NC_000001.11:g.{}_{}CTG[1_5]TTG[1]ACG[1_5]",
        CORE_START, end
    );
    normalize_strict(provider, &input)
        .expect("no anchored units → no validation; must accept arbitrary tract");
}

/// Leading-Exact-only (no trailing Exact): the prefix path still
/// works exactly as before — confirms the new suffix path doesn't
/// break the prefix-only case.
#[test]
fn leading_exact_only_prefix_validation_regression() {
    // CTG[2] (6 bp prefix) + TTG[1_3] (skip) — no trailing Exact.
    let core = "CTGCTGTTGTTG"; // 12 bp: 6 CTG + 6 TTG (2 TTGs in the range)
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]", CORE_START, end);
    normalize_strict(provider, &input)
        .expect("leading-Exact-only with matching prefix must accept (regression)");
}

/// Leading prefix mismatch is still detected (regression for the
/// existing prefix-validation path).
#[test]
fn leading_exact_prefix_mismatch_still_rejected_strict() {
    // CTG[2] declared but the actual leading 6 bp is GGGGGG.
    let core = "GGGGGGTTGACGACG"; // 15 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[2]", CORE_START, end);
    let err = normalize_strict(provider, &input)
        .expect_err("leading-Exact prefix mismatch must still reject");
    assert!(is_ref_mismatch(&err), "expected RefMismatch; got: {err:?}");
}
