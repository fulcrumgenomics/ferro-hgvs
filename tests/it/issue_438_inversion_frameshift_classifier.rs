#![cfg(feature = "web-service")]
//! Issue #438 — extend the conservative-skip frameshift classifier
//! contract (established for Delins by #394 item 1, and for Deletion
//! and Duplication by #427) to the Inversion arm.
//!
//! The pre-#438 Inversion arm in
//! `src/service/handlers/effect.rs::analyze_na_edit` hardcoded
//! `is_frameshift = true` for every inversion. Per
//! `assets/hgvs-nomenclature/docs/recommendations/DNA/inversion.md`,
//! an inversion reverse-complements its declared range in place — the
//! DNA-level net length change is **zero** by construction. A
//! frameshift requires `net_delta % 3 != 0`, so an inversion is
//! always in-frame at the DNA level regardless of span.
//!
//! Inversion may still disrupt protein function (the inverted bases
//! code for different amino acids), but that's a missense / in-frame
//! change, not a frameshift. The fix matches the
//! Delins/Deletion/Duplication arms' contract: `ref_len = alt_len =
//! span_len` when computable, conservative-skip (`0, 0`) when not.

use ferro_hgvs::service::handlers::effect::{
    analyze_na_edit, span_len_from_cds_interval, span_len_from_genome_interval,
    span_len_from_tx_interval,
};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

// =============================================================================
// CDS axis — in-frame regardless of span
// =============================================================================

/// `c.100_102inv`: 3-bp inversion, in-frame.
#[test]
fn cds_inv_3bp_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_102inv").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "inversion");
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 3, "inversion: alt_len == ref_len by construction");
    assert!(
        !is_fs,
        "inversion is always in-frame at the DNA level (net delta = 0)",
    );
}

/// `c.100_101inv`: 2-bp inversion. Span 2 → frameshift would naively
/// fire on `len % 3 != 0`, but inversion's net delta is 0, so the
/// classifier MUST report in-frame. This is the key bug the previous
/// `is_frameshift = true` hardcode masked.
#[test]
fn cds_inv_2bp_still_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_101inv").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 2);
    assert_eq!(alt_len, 2);
    assert!(
        !is_fs,
        "2-bp inversion is in-frame; previous hardcode mis-classified as frameshift",
    );
}

/// 5-bp inversion: span doesn't divide 3, but inversion is still in-frame.
#[test]
fn cds_inv_5bp_still_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_104inv").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 5);
    assert_eq!(alt_len, 5);
    assert!(!is_fs);
}

// =============================================================================
// Genome axis — same in-frame contract
// =============================================================================

#[test]
fn genome_inv_in_frame_regardless_of_span() {
    let v = parse_hgvs("NC_000001.11:g.100_104inv").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "inversion");
    assert_eq!(ref_len, 5);
    assert_eq!(alt_len, 5);
    assert!(!is_fs);
}

// =============================================================================
// Tx (n.) axis — same in-frame contract
// =============================================================================

#[test]
fn tx_inv_in_frame_regardless_of_span() {
    let v = parse_hgvs("NR_001234.1:n.100_104inv").unwrap();
    let HgvsVariant::Tx(tv) = v else {
        panic!("expected Tx")
    };
    let span = span_len_from_tx_interval(&tv.loc_edit.location).expect("span computable");
    let edit = tv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, Some(span));
    assert_eq!(ref_len, 5);
    assert_eq!(alt_len, 5);
    assert!(!is_fs);
}

// =============================================================================
// Conservative-skip when span is undecidable
// =============================================================================

/// Intronic inversion: span_len = None → `ref_len = 0, alt_len = 0,
/// is_frameshift = false`. Matches the conservative-skip contract of
/// the other arms.
#[test]
fn cds_inv_intronic_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.100+5_120inv").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(span.is_none(), "intronic-offset span must return None",);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, span);
    assert_eq!(ref_len, 0);
    assert_eq!(alt_len, 0);
    assert!(
        !is_fs,
        "intronic inversion: conservative skip → not frameshift",
    );
}

/// Mixed-axis inversion (`c.-1_1inv`): span_len = None (5'UTR ↔ CDS).
/// Conservative-skip.
#[test]
fn cds_inv_mixed_axis_conservative_skip() {
    let v = parse_hgvs("NM_000001.1:c.-1_1inv").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(span.is_none(), "5'UTR-to-CDS span must return None");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(!is_fs);
}

// =============================================================================
// Explicit sequence / length take priority over span_len (regression)
// =============================================================================

/// `c.100_102invATG`: parser sets `sequence = Some("ATG")` → `len = 3`.
/// Span_len also gives 3 — same answer here, but pin that the
/// explicit form takes priority via the fallback chain (mirrors
/// `cds_del_explicit_sequence_takes_priority_in_frame` in
/// `tests/issue_427_del_dup_frameshift_classifier.rs`).
#[test]
fn cds_inv_explicit_sequence_takes_priority_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_102invATG").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, ref_len, alt_len) = analyze_na_edit(edit, span);
    assert_eq!(kind, "inversion");
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 3);
    assert!(!is_fs);
}

/// `c.100_102inv3`: parser sets `length = Some(3)`, `sequence = None`.
/// The explicit `length` takes priority over `span_len`.
#[test]
fn cds_inv_explicit_length_takes_priority_in_frame() {
    let v = parse_hgvs("NM_000001.1:c.100_102inv3").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, ref_len, alt_len) = analyze_na_edit(edit, span);
    assert_eq!(ref_len, 3);
    assert_eq!(alt_len, 3);
    assert!(!is_fs);
}
