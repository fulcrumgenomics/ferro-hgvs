#![cfg(feature = "web-service")]
//! Issue #427 — extend the #394-item-1 frameshift-classifier fix
//! (which addressed `NaEdit::Delins`) to `NaEdit::Deletion` and
//! `NaEdit::Duplication`.
//!
//! The pre-#427 Deletion / Duplication arms in
//! `src/service/handlers/effect.rs::analyze_na_edit` used
//! `unwrap_or(1)` as the fallback when neither `sequence` nor `length`
//! was set on the parsed edit. Canonical short-form variants like
//! `c.100_102del` (parser sets `sequence: None, length: None`) thus
//! fell into the `unwrap_or(1)` branch and were classified as
//! `is_frameshift = true` because `1 % 3 != 0` — incorrect: the
//! actual span is 3 bp, in-frame.
//!
//! This issue mirrors the #394 item 1 fix: prefer the explicit
//! `sequence` / `length`, fall back to the position-interval
//! `span_len`, and conservative-skip (`is_frameshift = false`,
//! `ref_len = 0`) when none of those is available.
//!
//! # Spec basis
//!
//! `assets/hgvs-nomenclature/docs/recommendations/protein/frameshift.md`:
//! frameshift is a CDS-level concept tied to `net_delta % 3 != 0`,
//! where `net_delta = alt_len - ref_len`. For deletion alt_len = 0;
//! for duplication alt_len = ref_len * 2 → net_delta = ref_len, so
//! the in-frame test reduces to `ref_len % 3 != 0` in both cases.

//! **Test scope note:** the `tx_del_*` and `tx_dup_*` tests below
//! invoke `analyze_na_edit` directly with `Some(span)`. The production
//! `Tx` caller in `analyze_variant_effect` (`effect.rs:177`) passes
//! `None` because the tx-axis effect description only uses
//! `edit_type` and discards the frameshift signal. The tests here
//! exercise the helper directly, pinning its tx-axis behavior in case
//! a future caller starts using the frameshift signal — they do NOT
//! reflect a current production code path.

use ferro_hgvs::service::handlers::effect::{
    analyze_na_edit, span_len_from_cds_interval, span_len_from_genome_interval,
    span_len_from_tx_interval,
};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

// =============================================================================
// Deletion — span-derived classification
// =============================================================================

/// `c.100_102del`: parser sets `sequence: None, length: None`. The
/// span is 3 bp → `is_frameshift = false` via `span_len`. Before #427
/// this returned `true` because `unwrap_or(1)` defaulted to 1 bp.
#[test]
fn cds_del_in_frame_3bp_via_span_len() {
    let v = parse_hgvs("NM_000001.1:c.100_102del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "deletion");
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 0);
    assert!(!is_fs, "3 bp deletion is in-frame");
}

/// `c.100_101del`: 2 bp deletion → frameshift.
#[test]
fn cds_del_frameshift_2bp_via_span_len() {
    let v = parse_hgvs("NM_000001.1:c.100_101del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "deletion");
    assert_eq!(ref_len, 2);
    assert!(is_fs, "2 bp deletion is a frameshift");
}

/// Single-position `c.100del`: span_len = 1 (start == end). 1 bp
/// deletion → frameshift.
#[test]
fn cds_del_single_position_is_frameshift() {
    let v = parse_hgvs("NM_000001.1:c.100del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 1);
    assert!(is_fs, "1 bp deletion is a frameshift");
}

/// Intronic range `c.100+5_120del`: span is undefined from base
/// coordinates alone (`+5` offset can't resolve without intron
/// length). `span_len = None`, `sequence: None, length: None` →
/// conservative skip (`is_frameshift = false`).
#[test]
fn cds_del_intronic_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.100+5_120del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(
        span.is_none(),
        "span_len must be None for intronic-offset endpoints",
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(
        !is_fs,
        "intronic deletion → conservative is_frameshift=false (no panic)",
    );
}

/// Unknown CDS start endpoint (`c.?_123del`): `span_len = None`
/// (the `c.?` placeholder has no integer base), and the parser sets
/// `sequence: None, length: None`, so the fallback chain bottoms out at
/// the conservative skip — `ref_len = 0`, `alt_len = 0`,
/// `is_frameshift = false`. Pins that an unknown-endpoint deletion is
/// never mis-classified as an in-frame (or frameshift) call.
#[test]
fn cds_del_unknown_endpoint_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.?_123del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(
        span.is_none(),
        "span_len must be None for an unknown (c.?) endpoint",
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, span);
    assert_eq!(kind, "deletion");
    assert_eq!(ref_len, 0, "undecidable span reports ref_len = 0");
    assert_eq!(alt_len, 0, "undecidable span reports alt_len = 0");
    assert!(
        !is_fs,
        "unknown-endpoint deletion → conservative is_frameshift=false"
    );
}

/// Unknown CDS end endpoint on a duplication (`c.123_?dup`): same
/// conservative-skip contract as the deletion case above.
#[test]
fn cds_dup_unknown_endpoint_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.123_?dup").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(
        span.is_none(),
        "span_len must be None for an unknown endpoint"
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, span);
    assert_eq!(kind, "duplication");
    assert_eq!(ref_len, 0);
    assert_eq!(alt_len, 0);
    assert!(
        !is_fs,
        "unknown-endpoint duplication → conservative is_frameshift=false"
    );
}

/// Tx (n.) axis: 3 bp deletion → in-frame.
#[test]
fn tx_del_in_frame_via_span_len() {
    let v = parse_hgvs("NR_001234.1:n.100_102del").unwrap();
    let HgvsVariant::Tx(tv) = v else {
        panic!("expected Tx")
    };
    let span = span_len_from_tx_interval(&tv.loc_edit.location).expect("span computable");
    let edit = tv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 3);
    assert!(!is_fs);
}

/// Tx axis: 2 bp deletion → frameshift.
#[test]
fn tx_del_frameshift_2bp_via_span_len() {
    let v = parse_hgvs("NR_001234.1:n.100_101del").unwrap();
    let HgvsVariant::Tx(tv) = v else {
        panic!("expected Tx")
    };
    let span = span_len_from_tx_interval(&tv.loc_edit.location).expect("span computable");
    let edit = tv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 2);
    assert!(is_fs);
}

/// 5'UTR-only deletion (`c.-10_-8del`): both endpoints in the 5'UTR,
/// `span_len` IS computable (3 bp). Should be classified as in-frame.
#[test]
fn cds_del_5utr_only_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.-10_-8del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(!is_fs);
}

/// Mixed-axis deletion (`c.-1_1del`): span_len = None (5'UTR ↔ CDS
/// boundary). Conservative skip.
#[test]
fn cds_del_mixed_axis_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.-1_1del").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(span.is_none(), "5'UTR↔CDS span returns None");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(!is_fs);
}

// =============================================================================
// Duplication — same fallback chain
// =============================================================================

/// `c.100_102dup`: 3 bp duplication, in-frame. `alt_len = 2 * ref_len`.
#[test]
fn cds_dup_in_frame_3bp_via_span_len() {
    let v = parse_hgvs("NM_000001.1:c.100_102dup").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "duplication");
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 6);
    assert!(!is_fs);
}

/// `c.100_101dup`: 2 bp duplication → frameshift.
#[test]
fn cds_dup_frameshift_2bp_via_span_len() {
    let v = parse_hgvs("NM_000001.1:c.100_101dup").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 2);
    assert_eq!(alt_len, 4);
    assert!(is_fs);
}

/// Single-position `c.100dup`: 1 bp duplication → frameshift.
#[test]
fn cds_dup_single_position_is_frameshift() {
    let v = parse_hgvs("NM_000001.1:c.100dup").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 1);
    assert!(is_fs);
}

/// Intronic duplication: conservative skip.
#[test]
fn cds_dup_intronic_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.100+5_120dup").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(span.is_none());
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(!is_fs);
}

/// Genome axis duplication: 3 bp → in-frame.
#[test]
fn genome_dup_in_frame_via_span_len() {
    let v = parse_hgvs("NC_000001.11:g.100_102dup").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 6);
    assert!(!is_fs);
}

// =============================================================================
// Explicit sequence/length take priority over span_len (regression)
// =============================================================================

/// `c.100_102delAAA`: parser sets `sequence = Some("AAA")` →
/// `len = 3`. Span_len would also give 3 — same answer, but pin that
/// the explicit form takes priority via the fallback chain.
#[test]
fn cds_del_explicit_sequence_takes_priority_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_102delAAA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, span);
    assert_eq!(ref_len, 3);
    assert!(!is_fs);
}

/// `c.100_102del3`: parser sets `length = Some(3)`, `sequence = None`.
/// The explicit `length` takes priority over `span_len` via the
/// fallback chain. Pin in-frame classification.
#[test]
fn cds_del_explicit_length_takes_priority_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_102del3").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, _) = analyze_na_edit(edit, span);
    assert_eq!(kind, "deletion");
    assert_eq!(ref_len, 3);
    assert!(!is_fs);
}

/// `c.100_101del2`: 2-bp explicit length, frameshift.
#[test]
fn cds_del_explicit_length_frameshift_2bp() {
    let v = parse_hgvs("NM_000001.1:c.100_101del2").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, _) = analyze_na_edit(edit, span);
    assert_eq!(ref_len, 2);
    assert!(is_fs);
}

/// `c.100_102dupAAA`: explicit-sequence duplication, in-frame.
#[test]
fn cds_dup_explicit_sequence_takes_priority_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_102dupAAA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, span);
    assert_eq!(kind, "duplication");
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 6);
    assert!(!is_fs);
}

/// Genome deletion: 3-bp in-frame via span_len.
#[test]
fn genome_del_in_frame_via_span_len() {
    let v = parse_hgvs("NC_000001.11:g.100_102del").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "deletion");
    assert_eq!(ref_len, 3);
    assert!(!is_fs);
}
