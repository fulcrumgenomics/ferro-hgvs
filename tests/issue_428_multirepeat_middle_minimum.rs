//! Issue #428 — tighten `validate_multirepeat_tract` so the
//! `prefix_bytes.len() + suffix_bytes.len() == actual_bytes.len()`
//! boundary is no longer accepted unconditionally. PR #421 (#395 item 1)
//! introduced the right-anchored suffix walk + the prefix/suffix length
//! check and acknowledged inline that the `==` case was lax: it lets
//! through declarations like `[Exact(2), Range(1,3), Exact(2)]` over a
//! 12-bp span, where the middle's minimum (1 copy = 3 bp) plus prefix
//! (6) plus suffix (6) = 15 bp > 12 bp span.
//!
//! This issue closes that slack by computing the minimum bp the
//! middle must occupy from each non-`Exact` unit's count semantics:
//!
//!   - `Range(lo, _)`              → minimum `lo` copies
//!   - `UncertainRange(lo, _)`     → minimum `lo` copies
//!   - `MinUncertain(lo)`          → minimum `lo` copies
//!   - `MaxUncertain(_)`, `Unknown` → minimum 0 copies
//!
//! Reject when `prefix + middle_min + suffix > span`.
//!
//! Existing accept-`==` pin `interleaved_prefix_and_trailing_exact_consistent_accepted_strict`
//! (15-bp span = prefix 6 + middle_min 3 + suffix 6) continues to pass:
//! the middle min for `Range(1,3)` is exactly 3 bp, which fills the
//! span cleanly.

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
// New behavior: nonzero-min middle rejects `prefix + suffix == span`
// =============================================================================

/// `[Exact(2), Range(1,3), Exact(2)]` over a 12-bp span. Before #428
/// this was accepted because `prefix(6) + suffix(6) == 12`. The middle's
/// `Range(1,3)` requires ≥ 1 copy (3 bp), so the declared tract demands
/// at least 15 bp of span — must reject in strict mode.
#[test]
fn nonzero_min_middle_rejects_equal_span_strict() {
    // 12 bp: 6 CTG + 6 ACG. No room for the required ≥1 TTG middle copy.
    let core = "CTGCTGACGACG";
    assert_eq!(core.len(), 12);
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[2]", CORE_START, end);
    let err = normalize_strict(provider, &input).expect_err(
        "[Exact(2), Range(1,3), Exact(2)] over 12-bp span must reject in strict — \
         middle minimum (3 bp) cannot fit",
    );
    assert!(
        is_ref_mismatch(&err),
        "expected ReferenceMismatch from middle-minimum overflow; got: {err:?}",
    );
}

/// Lenient mode emits a `RefSeqMismatch` warning for the same input.
#[test]
fn nonzero_min_middle_warns_equal_span_lenient() {
    let core = "CTGCTGACGACG"; // 12 bp
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

// =============================================================================
// Still-accepted shapes: zero-min middle variants
// =============================================================================

/// `[Exact(2), Unknown, Exact(2)]` over a 12-bp span: `Unknown` admits
/// 0 copies, so the declared minimum is `6 + 0 + 6 = 12` — must accept.
#[test]
fn unknown_middle_admits_zero_count_accepted_strict() {
    let core = "CTGCTGACGACG"; // 12 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[?]ACG[2]", CORE_START, end);
    normalize_strict(provider, &input)
        .expect("[Exact(2), Unknown, Exact(2)] over a 12-bp span must accept (0 copies admitted)");
}

/// `[Exact(2), Range(0,3), Exact(2)]` over a 12-bp span: the explicit
/// range starts at 0, so 0 copies are admitted. Declared minimum is
/// `6 + 0 + 6 = 12` — must accept.
#[test]
fn range_starting_at_zero_admits_zero_count_accepted_strict() {
    let core = "CTGCTGACGACG"; // 12 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[0_3]ACG[2]", CORE_START, end);
    normalize_strict(provider, &input).expect(
        "[Exact(2), Range(0,3), Exact(2)] over a 12-bp span must accept \
         (lo=0 admits 0 copies)",
    );
}

/// `[Exact(2), MaxUncertain(3), Exact(2)]` over a 12-bp span: `[?_3]`
/// admits 0 copies. Must accept.
#[test]
fn max_uncertain_middle_admits_zero_count_accepted_strict() {
    let core = "CTGCTGACGACG"; // 12 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[?_3]ACG[2]", CORE_START, end);
    normalize_strict(provider, &input).expect(
        "[Exact(2), MaxUncertain(3), Exact(2)] over 12-bp span must accept \
         (MaxUncertain admits 0 copies)",
    );
}

// =============================================================================
// Boundary regressions: the existing accept-`==` pin still passes
// =============================================================================

/// The lax-acceptance pin from #395 item 1: `[Exact(2), Range(1,3),
/// Exact(2)]` over a **15**-bp span. middle_min (3 bp) exactly fills
/// the gap → still accepted under the tightened check.
#[test]
fn existing_consistent_accepted_pin_still_accepts_under_tightening() {
    let core = "CTGCTGTTGACGACG"; // 15 bp = 6 + 3 + 6
    assert_eq!(core.len(), 15);
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_3]ACG[2]", CORE_START, end);
    normalize_strict(provider, &input).expect(
        "15-bp span exactly accommodates prefix(6) + middle_min(3) + suffix(6) — \
         must still accept",
    );
}

/// `MinUncertain` (`[lo_?]`) carries the same `lo` floor as `Range(lo, ...)`.
/// `[Exact(2), MinUncertain(2), Exact(2)]` over a 12-bp span: middle
/// minimum is 2 copies = 6 bp, so declared min = 6 + 6 + 6 = 18 bp.
/// Must reject the 12-bp span.
#[test]
fn min_uncertain_middle_rejects_short_span_strict() {
    let core = "CTGCTGACGACG"; // 12 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[2_?]ACG[2]", CORE_START, end);
    let err = normalize_strict(provider, &input).expect_err(
        "[Exact(2), MinUncertain(2), Exact(2)] over 12-bp span must reject — \
         middle minimum (6 bp) plus anchors (12 bp) overflows",
    );
    assert!(
        is_ref_mismatch(&err),
        "expected ReferenceMismatch from middle-minimum overflow; got: {err:?}",
    );
}

/// `UncertainRange(lo, _)` — `[(lo_hi)]` form — carries the same
/// minimum `lo`. The spec-recommended uncertainty form
/// (`assets/hgvs-nomenclature/docs/recommendations/DNA/repeated.md` —
/// `c.1210-33_1210-6GT[(9_13)]T[(4_8)]`) accepts `lo > 0`. This
/// rejection complements `nonzero_min_middle_rejects_equal_span_strict`
/// which covers the explicit-`Range` shape; this one covers the
/// uncertain-paren shape.
#[test]
fn uncertain_range_middle_rejects_short_span_strict() {
    let core = "CTGCTGACGACG"; // 12 bp
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!(
        "NC_000001.11:g.{}_{}CTG[2]TTG[(1_3)]ACG[2]",
        CORE_START, end
    );
    let err = normalize_strict(provider, &input).expect_err(
        "[Exact(2), UncertainRange(1,3), Exact(2)] over 12-bp span must reject — \
         the uncertain-paren range still pins a `lo=1` floor",
    );
    assert!(
        is_ref_mismatch(&err),
        "expected ReferenceMismatch for UncertainRange middle; got: {err:?}",
    );
}

/// Interior `Exact(n)` units between two non-`Exact` middle units
/// contribute their full bp to the middle minimum. Without this
/// accounting the interior would be incorrectly skipped.
///
/// `[Exact(2), Range(1,3), Exact(1), Range(1,3), Exact(2)]` declared
/// over an 18-bp span. prefix = CTG[2] (6 bp), suffix = ACG[2] (6 bp),
/// middle min = TTG[1] + GGG[1] + TTG[1] = 3 + 3 + 3 = 9 bp; declared
/// min = 6 + 9 + 6 = 21 bp. 18 < 21, must reject.
#[test]
fn multi_non_exact_middle_with_interior_exact_rejects_strict() {
    // 18 bp ref: 6 CTG + 3 TTG + 3 GGG + 0 (middle TTG must be skipped) + 6 ACG
    // But that's only 18 bp not 21; we declare 21+ via middle min.
    let core = "CTGCTGTTGGGGTTGACGACG"; // 21 bp (just shy of test) — but we want SHORT span
                                        // Actually use a span shorter than the declared minimum.
    let short_core = "CTGCTGTTGGGGACGACG"; // 18 bp
    assert_eq!(short_core.len(), 18);
    let padded_seq = padded(short_core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (short_core.len() as u64) - 1;
    let input = format!(
        "NC_000001.11:g.{}_{}CTG[2]TTG[1_3]GGG[1]TTG[1_3]ACG[2]",
        CORE_START, end
    );
    let err = normalize_strict(provider, &input).expect_err(
        "[Exact(2), Range(1,3), Exact(1), Range(1,3), Exact(2)] over 18-bp span must reject — \
         interior Exact(1) is counted into middle_min (3 bp), pushing declared min to 21 bp",
    );
    assert!(
        is_ref_mismatch(&err),
        "expected ReferenceMismatch from interior-Exact middle accounting; got: {err:?}",
    );
    // Without proper interior accounting, declared_min would be 6+3+3+6 = 18,
    // exactly equal to the span, and would accept incorrectly.
    let _ = core;
}

/// `MinUncertain(lo)` accept case: `[Exact(2), MinUncertain(1), Exact(2)]`
/// over a 15-bp span. middle_min = 1 copy × 3 bp = 3 bp, declared min =
/// 6 + 3 + 6 = 15 bp. Must accept (boundary case for `MinUncertain`).
#[test]
fn min_uncertain_middle_accepts_when_span_matches_floor_strict() {
    let core = "CTGCTGTTGACGACG"; // 15 bp
    assert_eq!(core.len(), 15);
    let padded_seq = padded(core);
    let provider = g_provider("NC_000001.11", &padded_seq);
    let end = CORE_START + (core.len() as u64) - 1;
    let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1_?]ACG[2]", CORE_START, end);
    normalize_strict(provider, &input).expect(
        "[Exact(2), MinUncertain(1), Exact(2)] over 15-bp span must accept — \
         middle_min (3 bp) exactly fills the gap",
    );
}
