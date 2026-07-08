#![cfg(feature = "dev")]
//! Issue #394 follow-up — `predict_protein_consequence` must honor the
//! `is_frameshift` decision from `analyze_na_edit`, NOT recompute it
//! from `(alt_len - ref_len)`.
//!
//! Bug (CodeRabbit, post-merge follow-up review of #419): the previous
//! body of `predict_protein_consequence` overwrote the classifier's
//! decision with `(alt_len as i64 - ref_len as i64).rem_euclid(3) != 0`.
//! That undoes the conservative fallback for shapes where the
//! classifier returned `(_, false, 0, alt_len)` because span_len was
//! unknown. Concretely: an intronic-endpoint or unknown-span delins
//! with a single-nt literal insert (alt_len = 1, ref_len = 0
//! placeholder) used to be re-flagged as a frameshift at the protein
//! consequence step.
//!
//! This file pins the through-pass via a CDS variant that exercises
//! the classifier's `None`-span branch and asserts the resulting
//! `protein_consequence.is_frameshift` is `false`.

use axum::{extract::State, response::Json};
use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::service::config::ServiceConfig;
use ferro_hgvs::service::handlers::effect::predict_effect;
use ferro_hgvs::service::server::{AppState, HealthCache};
use ferro_hgvs::service::tools::ToolManager;
use ferro_hgvs::service::types::EffectRequest;
use std::sync::Arc;

/// Build a minimal `AppState` with a cdot mapper containing `NM_999999.1`
/// so `predict_protein_consequence` does not short-circuit on
/// `state.cdot.as_ref()?`.
fn state_with_test_transcript() -> AppState {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_999999.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );

    AppState {
        tool_manager: Arc::new(ToolManager::empty()),
        config: Arc::new(ServiceConfig::default()),
        cdot: Some(Arc::new(cdot)),
        reference: None,
        liftover: None,
        health_cache: HealthCache::default(),
    }
}

#[tokio::test]
async fn predict_effect_honors_classifier_is_frameshift_for_unknown_span_delins() {
    // `c.4_5+2delinsT` — endpoint `c.5+2` is intronic, so the cds
    // interval helper returns `None` for span_len. `analyze_na_edit`
    // therefore reports `(is_frameshift = false, ref_len = 0,
    // alt_len = 1)`. The previous body of `predict_protein_consequence`
    // would recompute `(1 - 0).rem_euclid(3) != 0` → `true` and emit a
    // `p.(?2fs)` despite the classifier's "unknown, report false"
    // decision. After the fix the protein consequence must NOT be a
    // frameshift.
    let state = state_with_test_transcript();
    let request = EffectRequest {
        hgvs: "NM_999999.1:c.4_5+2delinsT".to_string(),
        include_nmd: false,
    };

    let response = predict_effect(State(state), Json(request))
        .await
        .expect("predict_effect must succeed for valid input")
        .0;

    assert!(
        response.error.is_none(),
        "predict_effect must parse cleanly; got error: {:?}",
        response.error,
    );
    let protein = response
        .protein_consequence
        .expect("CDS delins must produce a protein_consequence");
    assert!(
        !protein.is_frameshift,
        "protein_consequence.is_frameshift must mirror the classifier's `false` decision \
         when span_len is None; was {:?}, hgvs_p = {}",
        protein.is_frameshift, protein.hgvs_p,
    );
    assert!(
        !protein.hgvs_p.contains("fs"),
        "hgvs_p must not render the frameshift suffix when classifier said false; got {}",
        protein.hgvs_p,
    );
}

#[tokio::test]
async fn predict_effect_in_frame_delins_with_known_span_remains_in_frame() {
    // Sanity baseline: `c.4_6delinsAAA` — same length on both sides
    // (ref_len = 3, alt_len = 3), classifier reports
    // `is_frameshift = false`. The previous (broken) recompute would
    // also have arrived at `false` here, so this test pins the
    // post-fix behavior on a known-good shape so the passthrough
    // doesn't silently regress on it.
    let state = state_with_test_transcript();
    let request = EffectRequest {
        hgvs: "NM_999999.1:c.4_6delinsAAA".to_string(),
        include_nmd: false,
    };

    let response = predict_effect(State(state), Json(request))
        .await
        .expect("predict_effect must succeed")
        .0;

    let protein = response
        .protein_consequence
        .expect("CDS delins must produce a protein_consequence");
    assert!(
        !protein.is_frameshift,
        "3->3 delins is in-frame; got is_frameshift = {}",
        protein.is_frameshift,
    );
}

#[tokio::test]
async fn predict_effect_frameshift_delins_with_known_span_is_frameshift() {
    // Sanity baseline: `c.4_6delinsAAAA` — ref_len = 3, alt_len = 4,
    // net_delta = +1, classifier reports `is_frameshift = true`. The
    // passthrough must surface that to the protein consequence layer.
    let state = state_with_test_transcript();
    let request = EffectRequest {
        hgvs: "NM_999999.1:c.4_6delinsAAAA".to_string(),
        include_nmd: false,
    };

    let response = predict_effect(State(state), Json(request))
        .await
        .expect("predict_effect must succeed")
        .0;

    let protein = response
        .protein_consequence
        .expect("CDS delins must produce a protein_consequence");
    assert!(
        protein.is_frameshift,
        "3->4 delins (+1) is a frameshift; got is_frameshift = {}, hgvs_p = {}",
        protein.is_frameshift, protein.hgvs_p,
    );
}
