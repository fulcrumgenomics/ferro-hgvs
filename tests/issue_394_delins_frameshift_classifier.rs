//! Issue #394 item 1 — effect classifier classifies delins frameshift
//! correctly.
//!
//! `src/service/handlers/effect.rs:216-219` previously hardcoded
//! `is_frameshift = true` for every `NaEdit::Delins`, with the comment
//! "Can't determine ref length easily". The information IS available
//! via the position interval (`end - start + 1`). This file pins the
//! spec-correct behavior.
//!
//! # Spec basis
//!
//! `assets/hgvs-nomenclature/docs/recommendations/protein/frameshift.md`:
//! frameshift is a CDS-level concept tied to `net_delta % 3 != 0`,
//! where `net_delta = alt_len - ref_len`.
//!
//! # Span-length computability
//!
//! The new `span_len_from_*_interval` helpers compute span length only
//! when both endpoints are integer-only (no offset / no special marker
//! / non-mixed axis sides). Anything else returns `None` and the
//! delins frameshift signal falls back to `false` — conservative, safer
//! than the previous always-true default.
//!
//! The `analyze_na_edit` helper is `pub` so this integration test can
//! exercise it directly without invoking the async `predict_effect`
//! HTTP handler.

use ferro_hgvs::service::handlers::effect::{
    analyze_na_edit, span_len_from_cds_interval, span_len_from_genome_interval,
};
use ferro_hgvs::{parse_hgvs, HgvsVariant};

#[test]
fn delins_in_frame_genomic_same_length() {
    // g.10_12delinsATG: ref_len = 3, alt_len = 3 → in-frame.
    let v = parse_hgvs("NC_000001.11:g.10_12delinsATG").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "delins");
    assert!(!is_fs, "3->3 delins is in-frame");
}

#[test]
fn delins_frameshift_genomic_3_to_1() {
    // g.10_12delinsA: ref_len = 3, alt_len = 1 → net_delta = -2 → frameshift.
    let v = parse_hgvs("NC_000001.11:g.10_12delinsA").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (kind, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert_eq!(kind, "delins");
    assert!(is_fs, "3->1 delins (-2) is a frameshift");
}

#[test]
fn delins_in_frame_genomic_3_to_6() {
    // g.10_12delinsATGCAT: ref_len = 3, alt_len = 6 → net_delta = +3 → in-frame.
    let v = parse_hgvs("NC_000001.11:g.10_12delinsATGCAT").unwrap();
    let HgvsVariant::Genome(gv) = v else {
        panic!("expected Genome")
    };
    let span = span_len_from_genome_interval(&gv.loc_edit.location).expect("span computable");
    let edit = gv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(!is_fs, "3->6 delins (+3) is in-frame");
}

#[test]
fn delins_frameshift_cds_2_to_3() {
    // c.10_11delinsATG: ref_len = 2, alt_len = 3 → net_delta = +1 → frameshift.
    let v = parse_hgvs("NM_000001.1:c.10_11delinsATG").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location).expect("span computable");
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, Some(span));
    assert!(is_fs, "2->3 delins (+1) is a frameshift");
}

#[test]
fn delins_intronic_span_unknown_frameshift_false() {
    // c.10+5_10+7delinsA: both endpoints intronic. Span can't be
    // computed from base alone (needs intron length); span_len = None →
    // conservative `is_frameshift = false`. Pins the new behavior
    // against the prior always-true default.
    let v = parse_hgvs("NM_000001.1:c.10+5_10+7delinsA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    assert!(
        span.is_none(),
        "intronic-endpoint interval must produce None span_len",
    );
    let edit = cv.loc_edit.edit.inner().unwrap();
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(
        !is_fs,
        "intronic delins with unknown span_len must conservatively report \
         is_frameshift = false (NOT the prior hardcoded true)",
    );
}

#[test]
fn delins_single_position_alt1_is_substitution_shape_but_classifier_says_delins() {
    // c.10delinsA: ref_len = 1, alt_len = 1. This is the substitution
    // canonical form per spec, but the parser may preserve it as
    // delins. Whichever form survives, the frameshift signal must be
    // `false`.
    let v = parse_hgvs("NM_000001.1:c.10delinsA").unwrap();
    let HgvsVariant::Cds(cv) = v else {
        panic!("expected Cds")
    };
    let edit = cv.loc_edit.edit.inner().unwrap();
    let span = span_len_from_cds_interval(&cv.loc_edit.location);
    let (_, is_fs, _, _) = analyze_na_edit(edit, span);
    assert!(
        !is_fs,
        "1->1 delins is in-frame regardless of canonical form"
    );
}
