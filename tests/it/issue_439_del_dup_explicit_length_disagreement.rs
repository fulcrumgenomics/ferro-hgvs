//! Issue #439 — extend W3016 length-mismatch detection to numeric-
//! length-suffix forms of `del` / `dup` / `inv`.
//!
//! The pre-#439 `try_detect_length_mismatch_at` walked `take_ref_seq_run`
//! (alphabetic IUPAC bases only) after the edit keyword. For shapes
//! like `c.100_102del3` / `c.100_101dup3` / `c.100_104inv3` —
//! parser-supported numeric-length-suffix forms — the ref-seq run was
//! empty and the detector silently let any range/length disagreement
//! through. This file pins the extended behavior: the detector now
//! also parses a numeric suffix as the declared length and compares
//! against the position-interval span. Disagreement surfaces W3016.
//!
//! `detect_del_size_suffix` (W3011, "deprecated size-suffix syntax")
//! continues to fire independently on these shapes — the two
//! diagnostics are complementary, not duplicative. A `c.100_102del3`
//! input emits W3011 (deprecated) regardless of agreement; W3016 only
//! fires when the size disagrees with the range.
//!
//! Spec basis: `assets/hgvs-nomenclature/docs/recommendations/general.md`
//! range semantics — for explicit-length / explicit-ref-seq forms,
//! `length == end - start + 1` must hold.

use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;

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
// Disagreement on the numeric-length-suffix shape → W3016
// =============================================================================

/// `c.100_102del5`: range 3 bp, declared 5 bp → mismatch.
#[test]
fn del_numeric_length_disagrees_with_range_emits_w3016() {
    assert_w3016("NM_000088.3:c.100_102del5");
}

/// `c.100_101dup3`: range 2 bp, declared 3 bp → mismatch.
#[test]
fn dup_numeric_length_disagrees_with_range_emits_w3016() {
    assert_w3016("NM_000088.3:c.100_101dup3");
}

/// `c.100_104inv2`: range 5 bp, declared 2 bp → mismatch.
#[test]
fn inv_numeric_length_disagrees_with_range_emits_w3016() {
    assert_w3016("NM_000088.3:c.100_104inv2");
}

/// `g.100_120del6`: range 21 bp, declared 6 bp → mismatch. Genomic
/// axis routes through the same detector — pin symmetry.
#[test]
fn genome_del_numeric_length_disagrees_emits_w3016() {
    assert_w3016("NC_000001.11:g.100_120del6");
}

// =============================================================================
// Agreement on the numeric-length-suffix shape → no W3016
// =============================================================================

/// `c.100_102del3`: range 3 bp, declared 3 bp → consistent. Must not
/// emit W3016 (W3011 deprecation still fires; pinned in
/// `del_size_with_matching_range_no_w3016_hit` below).
#[test]
fn del_numeric_length_matches_range_no_w3016() {
    assert_no_w3016("NM_000088.3:c.100_102del3");
}

/// `c.100_101dup2`: range 2 bp, declared 2 bp → consistent.
#[test]
fn dup_numeric_length_matches_range_no_w3016() {
    assert_no_w3016("NM_000088.3:c.100_101dup2");
}

/// `c.100_104inv5`: range 5 bp, declared 5 bp → consistent.
#[test]
fn inv_numeric_length_matches_range_no_w3016() {
    assert_no_w3016("NM_000088.3:c.100_104inv5");
}

// =============================================================================
// Interaction with W3011 (DelSizeSuffix deprecation)
// =============================================================================

/// `c.100del3`: single-position with size suffix. Per the existing
/// W3011 contract this is deprecated syntax (use `c.100_102del`
/// instead). #439 doesn't touch this path: no range to compare
/// against, so W3016 doesn't fire. W3011 still does.
#[test]
fn del_size_single_position_no_w3016_w3011_only() {
    let codes = warning_codes("NM_000088.3:c.100del3");
    assert!(
        codes.iter().any(|c| c == "W3011"),
        "expected W3011 (deprecated size-suffix), got {codes:?}",
    );
    assert!(
        !codes.iter().any(|c| c == "W3016"),
        "must not emit W3016 for single-position del<N>; got {codes:?}",
    );
}

/// `c.100_102del3`: range + size suffix, agreement → W3011 fires
/// (deprecation), W3016 does not. The deprecation aspect is
/// orthogonal to length agreement.
#[test]
fn del_size_with_matching_range_no_w3016_hit() {
    let codes = warning_codes("NM_000088.3:c.100_102del3");
    assert!(
        codes.iter().any(|c| c == "W3011"),
        "expected W3011 for `del<N>` deprecation regardless of range agreement, got {codes:?}",
    );
    assert!(
        !codes.iter().any(|c| c == "W3016"),
        "must not emit W3016 when declared length matches range; got {codes:?}",
    );
}

/// `c.100_102del5`: range + size suffix, DISagreement → both W3011
/// (deprecation) and W3016 (length mismatch) fire. Pin the
/// composition: the two phases are independent and complementary.
#[test]
fn del_size_with_disagreeing_range_emits_both_w3011_and_w3016() {
    let codes = warning_codes("NM_000088.3:c.100_102del5");
    assert!(
        codes.iter().any(|c| c == "W3011"),
        "expected W3011 for `del<N>` deprecation, got {codes:?}",
    );
    assert!(
        codes.iter().any(|c| c == "W3016"),
        "expected W3016 for the length disagreement, got {codes:?}",
    );
}

// =============================================================================
// Negative controls — existing detector behavior preserved
// =============================================================================

/// `c.100_109delAAAATTTGCC`: alphabetic ref seq path (existing pre-
/// #439 detection). Must continue to emit W3016 only when the bases
/// disagree (range 10 != 10 → no warning).
#[test]
fn explicit_alphabetic_ref_consistent_no_w3016_regression() {
    assert_no_w3016("NC_000001.11:g.100_109delAAAATTTGCC");
}

/// `delins<ins>` shape: no `<del>` ref-seq to compare against, so the
/// detector continues to skip (#439 doesn't change this). The
/// inserted-seq length is unrelated to the deleted-range length per
/// the spec.
#[test]
fn delins_short_form_still_skipped() {
    assert_no_w3016("NC_000001.11:g.100_102delinsATGCAT");
}

// =============================================================================
// Same-anchor intronic numeric suffix (cross-branch correctness)
// =============================================================================

/// `c.100+5_100+10del3`: both endpoints intronic, same anchor, same
/// sign → `range_len = 6` (handled by the (b) branch of the
/// detector). Numeric suffix says 3 → mismatch. Pin the cross-branch
/// correctness: the same-anchor intronic shape must compose cleanly
/// with the new numeric-suffix check.
#[test]
fn intronic_same_anchor_numeric_length_disagrees_emits_w3016() {
    assert_w3016("NM_000088.3:c.100+5_100+10del3");
}

/// `c.100+5_100+10del6`: same shape, range and length both 6 → no
/// W3016. Regression for the (b) branch composition.
#[test]
fn intronic_same_anchor_numeric_length_matches_no_w3016() {
    assert_no_w3016("NM_000088.3:c.100+5_100+10del6");
}

// =============================================================================
// dup<N> / inv<N> only emit W3016 (no W3011 — that's del-only)
// =============================================================================

/// `c.100_101dup3` and `c.100_104inv5`: `detect_del_size_suffix` only
/// matches the `del` keyword (see `b"del"`-only check in its impl).
/// `dup<N>` and `inv<N>` therefore do NOT receive a W3011
/// deprecation warning. The W3016 length-mismatch check is the only
/// surface diagnostic for these shapes. Pin the asymmetry so a future
/// extension of W3011 to dup/inv doesn't silently change the W3016
/// composition contract.
#[test]
fn dup_numeric_length_disagrees_emits_w3016_only_no_w3011() {
    let codes = warning_codes("NM_000088.3:c.100_101dup3");
    assert!(
        codes.iter().any(|c| c == "W3016"),
        "expected W3016 for dup<N> disagreement, got {codes:?}",
    );
    assert!(
        !codes.iter().any(|c| c == "W3011"),
        "W3011 is del-only; must not fire for dup<N>, got {codes:?}",
    );
}

#[test]
fn inv_numeric_length_disagrees_emits_w3016_only_no_w3011() {
    let codes = warning_codes("NM_000088.3:c.100_104inv2");
    assert!(
        codes.iter().any(|c| c == "W3016"),
        "expected W3016 for inv<N> disagreement, got {codes:?}",
    );
    assert!(
        !codes.iter().any(|c| c == "W3011"),
        "W3011 is del-only; must not fire for inv<N>, got {codes:?}",
    );
}

// =============================================================================
// Compound-allele inner-member numeric suffix
// =============================================================================

/// `c.[100A>G;200_205del3]`: outer-axis compound bracket with inner
/// numeric-suffix member. The outer-scan walks inner members
/// inheriting the outer accession+axis (per #429's axis-tracking fix),
/// and dispatches `try_detect_length_mismatch_at` for each. The
/// numeric-suffix branch must fire on the inner `200_205del3`
/// disagreement (range 6, declared 3).
#[test]
fn compound_allele_inner_numeric_length_disagrees_emits_w3016() {
    assert_w3016("NM_000088.3:c.[100A>G;200_205del3]");
}

/// `c.[100A>G;200_205del6]`: inner agreement (range 6, declared 6) —
/// must not emit W3016.
#[test]
fn compound_allele_inner_numeric_length_matches_no_w3016() {
    assert_no_w3016("NM_000088.3:c.[100A>G;200_205del6]");
}
