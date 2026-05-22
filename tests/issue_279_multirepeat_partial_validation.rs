//! Audit for issue #279 — partial-validation mode for
//! `NaEdit::MultiRepeat` when some units carry non-`Exact` counts.
//!
//! Background (PR #215, closing #214): the MultiRepeat content-vs-span
//! consistency check fires only when *every* unit has
//! `RepeatCount::Exact(_)`. Mixed shapes with `Range`, `MinUncertain`,
//! `MaxUncertain`, or `Unknown` counts were skipped entirely because
//! the total expected length is ambiguous.
//!
//! This audit pins the partial-validation behavior: for a description
//! like `g.START_ENDACG[3]GT[2..5]`, the leading `Exact` units form a
//! known prefix whose bases the reference must match — only the tail
//! that begins at the first non-`Exact` unit is ambiguous. So we
//! validate the prefix and skip the tail.
//!
//! Contract:
//!   1. All-`Exact` MultiRepeat — full validation (regression guard,
//!      identical to PR #215).
//!   2. `Exact` prefix + non-`Exact` tail — validate the prefix bases
//!      against the reference; do NOT validate length or tail bases.
//!      A mismatch in the prefix surfaces as `RefSeqMismatch` (strict
//!      reject / lenient warn) using the same flow as PR #215.
//!   3. First AND last unit non-`Exact` — no validation at all (no
//!      anchored end), identical to today's skip behavior.
//!   4. Mismatch in the unvalidated middle bases does NOT fire — the
//!      middle is genuinely ambiguous and the normalizer cannot tell a
//!      "wrong" middle from one whose count happens to be at the
//!      lower bound of the declared range.
//!
//! Issue #395 item 1 extended this contract: a trailing-`Exact`
//! suffix is now also right-anchored validated. The
//! "first-unit-non-Exact-but-last-is-Exact" case (previously skipped
//! per (3)) is now covered by
//! `tests/issue_395_multirepeat_interleaved_resumption.rs`.

use ferro_hgvs::error::FerroError;
use ferro_hgvs::normalize::NormalizationWarning;
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

// =============================================================================
// Test infrastructure (kept local to this file — mirrors issue_214's layout)
// =============================================================================

/// 256 bp of `N` flanking the variant tract on each side. `N` is inert
/// against the byte-equality scans the normalizer / validator uses for
/// repeat tracts (no IUPAC base equals `N` under
/// `eq_ignore_ascii_case`), so the boundaries never spuriously extend
/// the apparent tract.
const PAD: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\
     NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

fn padded(core: &str) -> String {
    format!("{}{}{}", PAD, core, PAD)
}

/// HGVS position of the first base of `core` in `padded(core)`.
/// Derived from `PAD.len() + 1`: PAD is 256 bp of `N`, so the first
/// core base is the 257th overall (1-based HGVS coordinates).
const CORE_START: u64 = PAD.len() as u64 + 1;

fn g_provider(accession: &str, padded_seq: &str) -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence(accession, padded_seq);
    p
}

fn normalize_strict(provider: MockProvider, input: &str) -> Result<String, FerroError> {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::strict());
    let variant = parse_hgvs(input).expect("parse failed");
    normalizer.normalize(&variant).map(|v| format!("{}", v))
}

fn normalize_lenient_warnings(
    provider: MockProvider,
    input: &str,
) -> (String, Vec<NormalizationWarning>) {
    let normalizer = Normalizer::with_config(provider, NormalizeConfig::lenient());
    let variant = parse_hgvs(input).expect("parse failed");
    let r = normalizer
        .normalize_with_warnings(&variant)
        .expect("lenient normalize failed");
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
// SECTION 1 — Regression guard: all-Exact MultiRepeat still validates
// =============================================================================

mod all_exact_regression {
    use super::*;

    /// All-Exact mixed-repeat over a matching reference — must pass
    /// (same contract as PR #215's `consistent_multirepeat_round_trips_strict`).
    #[test]
    fn all_exact_consistent_accepted_strict() {
        // 2 CTG + 1 TTG + 11 CTG = 14 trimers = 42 bp.
        let core = "CTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
        assert_eq!(core.len(), 42, "test setup: mixed tract must be 42 bp");
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + 42 - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[11]", CORE_START, end);
        normalize_strict(provider, &input)
            .expect("all-Exact consistent multirepeat must accept (regression for #215)");
    }

    /// All-Exact mixed-repeat with span ≠ sum — must still reject
    /// (regression).
    #[test]
    fn all_exact_span_mismatch_rejected_strict() {
        let core = "CTGCTGTTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        // Declare a span 24 bp shorter than the 42 bp sum.
        let short_end = CORE_START + 24 - 1;
        let input = format!(
            "NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[11]",
            CORE_START, short_end
        );
        let err = normalize_strict(provider, &input)
            .expect_err("all-Exact span-vs-sum mismatch must still reject");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 2 — Exact prefix + Range tail: prefix validated, tail skipped
// =============================================================================

mod exact_prefix_range_tail {
    use super::*;

    /// `g.START_END CTG[2]TTG[2..5]` over reference
    /// `CTGCTGTTGTTGTTG` (2 CTG + 3 TTG): the leading Exact prefix
    /// `CTG[2]` matches the reference's first 6 bp, and the tail
    /// `TTG[2..5]` is a legal range that covers the remaining 3 TTGs.
    /// Strict mode must accept (prefix validation passes, tail is
    /// genuinely ambiguous and skipped).
    #[test]
    fn matching_prefix_with_range_tail_accepted_strict() {
        let core = "CTGCTGTTGTTGTTG"; // 6 + 9 = 15 bp
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[(2_5)]", CORE_START, end);
        normalize_strict(provider, &input)
            .expect("matching Exact prefix + Range tail must be accepted (partial validation)");
    }

    /// Same shape, mismatch in the prefix — strict must reject.
    /// Reference's first 6 bp are `CTGCAG`, which is NOT `CTG[2]`, so
    /// the partial-validation prefix check fires.
    #[test]
    fn prefix_mismatch_with_range_tail_rejected_strict() {
        // First trimer is CTG but the second is CAG: prefix doesn't
        // match CTG[2].
        let core = "CTGCAGTTGTTGTTG"; // bad prefix, otherwise legal
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[(2_5)]", CORE_START, end);
        let err = normalize_strict(provider, &input)
            .expect_err("strict mode must reject when the Exact prefix bases mismatch reference");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// Prefix mismatch in lenient mode emits a `RefSeqMismatch`
    /// warning and keeps the input verbatim (same lenient surface as
    /// PR #215).
    #[test]
    fn prefix_mismatch_with_range_tail_warns_lenient() {
        let core = "CTGCAGTTGTTGTTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[(2_5)]", CORE_START, end);
        let (_out, warnings) = normalize_lenient_warnings(provider, &input);
        assert!(
            has_refseq_mismatch(&warnings),
            "expected RefSeqMismatch warning, got {:?}",
            warnings
        );
    }

    /// Mismatch is in the UNvalidated tail — partial validation must
    /// NOT fire because the tail (`TTG[2..5]`) has an ambiguous total
    /// length. The reference here packs `CAG` triplets where the
    /// description expects `TTG`s; the prefix `CTG[2]` is fine, so the
    /// partial check passes and the description is accepted.
    #[test]
    fn tail_mismatch_does_not_fire_strict() {
        // Prefix bases (first 6 bp) are CTGCTG — matches CTG[2].
        // Tail bases are CAGCAGCAG — does NOT match TTG, but the tail
        // is unvalidated.
        let core = "CTGCTGCAGCAGCAG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[(2_5)]", CORE_START, end);
        let _out = normalize_strict(provider, &input).expect(
            "tail mismatch must NOT fire — partial validation only checks the Exact prefix",
        );
    }
}

// =============================================================================
// SECTION 3 — Exact prefix + Unknown / MinUncertain tail
// =============================================================================

mod exact_prefix_unknown_tail {
    use super::*;

    /// `CTG[2]TTG[?]` (Unknown tail count) — prefix `CTG[2]` matches
    /// the reference's first 6 bp; tail is fully ambiguous. Must
    /// accept in strict mode.
    #[test]
    fn matching_prefix_with_unknown_tail_accepted_strict() {
        let core = "CTGCTGTTGTTG"; // anything after the prefix is fine
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[?]", CORE_START, end);
        normalize_strict(provider, &input)
            .expect("matching prefix + Unknown tail must be accepted (partial validation)");
    }

    /// Prefix mismatch with Unknown tail — must reject.
    #[test]
    fn prefix_mismatch_with_unknown_tail_rejected_strict() {
        // First trimer is CAG (not CTG) — prefix CTG[2] doesn't match.
        let core = "CAGCTGTTGTTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[?]", CORE_START, end);
        let err = normalize_strict(provider, &input)
            .expect_err("strict mode must reject prefix mismatch with Unknown tail");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// `CTG[2]TTG[3_?]` (MinUncertain tail count) — prefix
    /// validates, tail skipped. Must accept.
    #[test]
    fn matching_prefix_with_min_uncertain_tail_accepted_strict() {
        let core = "CTGCTGTTGTTGTTG"; // 2 CTG + 3 TTG
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[(3_?)]", CORE_START, end);
        normalize_strict(provider, &input)
            .expect("matching prefix + MinUncertain tail must be accepted");
    }
}

// =============================================================================
// SECTION 4 — First unit non-Exact: no validation possible
// =============================================================================

mod first_unit_non_exact {
    use super::*;

    /// First AND last unit non-Exact — neither end is anchored.
    /// Validation must be skipped entirely. A genuinely WRONG-looking
    /// reference must still be accepted in strict mode.
    ///
    /// (Pre-#395 this also held for `[Range, Exact]` shapes; #395
    /// item 1 added right-anchored validation, so the trailing-Exact
    /// variant is now covered by
    /// `issue_395_multirepeat_interleaved_resumption`.)
    #[test]
    fn first_and_last_unit_range_no_validation_strict() {
        // Reference has CAG triplets — clearly not CTG repeated — but
        // because neither end is anchored, no validation fires.
        let core = "CAGCAGCAGCAGCAG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[(2_5)]TTG[(1_3)]", CORE_START, end);
        let _ = normalize_strict(provider, &input).expect(
            "first-and-last-unit-non-Exact must skip validation (no anchored end to check)",
        );
    }

    /// Both ends `Unknown` count — same skip.
    #[test]
    fn first_and_last_unit_unknown_no_validation_strict() {
        let core = "CAGCAGCAGCAGCAG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[?]TTG[?]", CORE_START, end);
        let _ = normalize_strict(provider, &input)
            .expect("first-and-last-unit Unknown must skip validation (no anchored end to check)");
    }
}

// =============================================================================
// SECTION 5 — Multiple Exact units before the first non-Exact (>1 prefix unit)
// =============================================================================

mod multi_exact_prefix {
    use super::*;

    /// Two `Exact` units in the prefix, then a `Range` tail — the
    /// prefix is `CTG[2]TTG[1]` = 9 bp; reference's first 9 bp must
    /// match `CTGCTGTTG`.
    #[test]
    fn two_exact_units_then_range_accepted_strict() {
        // 2 CTG + 1 TTG = 9 bp prefix; then 3 more CTG = 9 bp tail.
        let core = "CTGCTGTTGCTGCTGCTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!(
            "NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[(2_5)]",
            CORE_START, end
        );
        normalize_strict(provider, &input).expect("two-Exact-prefix + Range tail must accept");
    }

    /// Same shape but the SECOND prefix unit mismatches the reference
    /// (`TTG` declared but ref has `CTG` there). The prefix walker
    /// must catch this even though the first unit matches.
    #[test]
    fn second_prefix_unit_mismatch_rejected_strict() {
        // Reference: 3 CTG + 3 CTG = 18 CTG. Description's second
        // prefix unit `TTG[1]` doesn't match the reference's 3rd
        // trimer (which is CTG, not TTG).
        let core = "CTGCTGCTGCTGCTGCTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!(
            "NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[(2_5)]",
            CORE_START, end
        );
        let err = normalize_strict(provider, &input)
            .expect_err("second prefix unit's bases mismatch reference — strict must reject");
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }
}

// =============================================================================
// SECTION 6 — Boundary tightness: last byte of the validated prefix
// =============================================================================
//
// The partial-validation path slices the reference as
// `&actual_bytes[..prefix_bytes.len()]` and compares byte-by-byte
// against `prefix_bytes`. Flipping ONE base at index
// `prefix_bytes.len() - 1` (the final byte of the validated prefix)
// must surface a mismatch. This locks the slice's exclusive upper bound
// — an off-by-one that shortened the slice would silently miss this
// flip, and an off-by-one that lengthened it would read past the
// prefix into the (ambiguous) tail.

mod prefix_boundary_tightness {
    use super::*;

    /// `CTG[2]TTG[1]CTG[(2_5)]` — Exact prefix is `CTGCTGTTG` (9 bp).
    /// The last byte of the prefix is at index 8 (the final `G` of
    /// `TTG`). Flip just that base to `A`: reference becomes
    /// `CTGCTGTTACTGCTGCTG` — the first 8 prefix bytes are correct, the
    /// 9th is wrong, the tail bytes (indices 9..) are irrelevant to the
    /// check. Strict mode must reject.
    #[test]
    fn flip_last_prefix_byte_rejected_strict() {
        // Prefix declared: CTG[2]TTG[1] = "CTGCTGTTG" (9 bp).
        // Core differs from the matching reference ONLY at index 8:
        // expected G, given A. Tail (indices 9..=17) matches a plausible
        // CTG run so we don't conflate prefix vs tail effects.
        let core = "CTGCTGTTACTGCTGCTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!(
            "NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[(2_5)]",
            CORE_START, end
        );
        let err = normalize_strict(provider, &input).expect_err(
            "flipping the last byte of the validated prefix must surface a ReferenceMismatch",
        );
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// Same single-byte flip in lenient mode — must emit a
    /// `RefSeqMismatch` warning (mirroring the lenient surface for the
    /// other prefix-mismatch cases).
    #[test]
    fn flip_last_prefix_byte_warns_lenient() {
        let core = "CTGCTGTTACTGCTGCTG";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let end = CORE_START + (core.len() as u64) - 1;
        let input = format!(
            "NC_000001.11:g.{}_{}CTG[2]TTG[1]CTG[(2_5)]",
            CORE_START, end
        );
        let (_out, warnings) = normalize_lenient_warnings(provider, &input);
        assert!(
            has_refseq_mismatch(&warnings),
            "expected RefSeqMismatch warning, got {:?}",
            warnings
        );
    }
}

// =============================================================================
// SECTION 7 — Reference span shorter than the declared Exact prefix
// =============================================================================
//
// The partial-validation arm has a defensive branch for
// `actual_bytes.len() < prefix_bytes.len()`: the slice
// `&actual_bytes[..prefix_bytes.len()]` would panic, and silently
// truncating the comparison would mask a genuine length / content
// mismatch. This branch surfaces a `RefSeqMismatch` instead. Pin both
// strict and lenient behavior here.

mod prefix_longer_than_span {
    use super::*;

    /// Declared description is `CTG[2]TTG[?]`. The Exact prefix is
    /// `CTG[2]` (6 bp), but the declared span covers only 4 bp — so the
    /// prefix runs off the end of the reference window. Strict must
    /// reject with `ReferenceMismatch`.
    #[test]
    fn short_span_for_prefix_rejected_strict() {
        // Reference content under the 4 bp window is `CTGC` (the first
        // 4 bp of a `CTGCTG...` tract). The leading bases agree with
        // the prefix, but the window is genuinely too short to contain
        // the full 6 bp prefix.
        let core = "CTGCTGNNNNNNNN";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        // Declare a 4 bp span (positions CORE_START..=CORE_START+3).
        let short_end = CORE_START + 4 - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[?]", CORE_START, short_end);
        let err = normalize_strict(provider, &input).expect_err(
            "reference span shorter than declared Exact prefix must be rejected in strict mode",
        );
        assert!(
            is_ref_mismatch(&err),
            "expected ReferenceMismatch, got {:?}",
            err
        );
    }

    /// Same short-span shape in lenient mode — must emit a
    /// `RefSeqMismatch` warning.
    #[test]
    fn short_span_for_prefix_warns_lenient() {
        let core = "CTGCTGNNNNNNNN";
        let padded_seq = padded(core);
        let provider = g_provider("NC_000001.11", &padded_seq);
        let short_end = CORE_START + 4 - 1;
        let input = format!("NC_000001.11:g.{}_{}CTG[2]TTG[?]", CORE_START, short_end);
        let (_out, warnings) = normalize_lenient_warnings(provider, &input);
        assert!(
            has_refseq_mismatch(&warnings),
            "expected RefSeqMismatch warning, got {:?}",
            warnings
        );
    }
}
