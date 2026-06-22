#![cfg(feature = "dev")]
//! Issue #806 — real amino-acid resolution and junction-based NMD in the
//! service effect handler, exercised end-to-end against a real
//! `MultiFastaProvider` (sequence-bearing transcripts + cdot CDS).
//!
//! Gated on `FERRO_MANIFEST` (skips in manifest-less CI, like the other
//! axis/real-reference tests). It verifies the production wiring added in
//! #806: `AppState.reference` is consulted by `predict_effect`, so a real
//! missense resolves to concrete residues (never the legacy `?` placeholder),
//! and the NMD prediction is produced for a real transcript.

use axum::{extract::State, response::Json};
use ferro_hgvs::reference::provider::ReferenceProvider;
use ferro_hgvs::service::config::ServiceConfig;
use ferro_hgvs::service::handlers::effect::predict_effect;
use ferro_hgvs::service::server::{AppState, HealthCache};
use ferro_hgvs::service::tools::ToolManager;
use ferro_hgvs::service::types::EffectRequest;
use ferro_hgvs::MultiFastaProvider;
use std::path::PathBuf;
use std::sync::Arc;

/// The manifest path from `FERRO_MANIFEST`, or `None` when the env var is unset.
///
/// Returns `Some` for any non-empty value — even a missing path — so that an
/// explicit opt-in that points at an invalid manifest fails loudly rather than
/// silently skipping. Only an *unset* `FERRO_MANIFEST` is a legitimate skip
/// (manifest-less CI); the caller asserts the path exists.
fn manifest_path() -> Option<PathBuf> {
    std::env::var_os("FERRO_MANIFEST").map(PathBuf::from)
}

/// Build an `AppState` whose `reference` and `cdot` both come from the real
/// manifest-loaded provider — mirroring what `create_app` does from the ferro
/// tool's reference_dir.
fn real_state(provider: Arc<MultiFastaProvider>) -> AppState {
    let cdot = provider
        .cdot_mapper()
        .expect("manifest provider must carry a cdot mapper")
        .clone();
    AppState {
        tool_manager: Arc::new(ToolManager::empty()),
        config: Arc::new(ServiceConfig::default()),
        cdot: Some(Arc::new(cdot)),
        reference: Some(provider as Arc<dyn ReferenceProvider + Send + Sync>),
        liftover: None,
        health_cache: HealthCache::default(),
    }
}

#[tokio::test]
async fn real_missense_resolves_concrete_residues() {
    let Some(path) = manifest_path() else {
        eprintln!("issue-806: skipping — FERRO_MANIFEST unset");
        return;
    };
    // FERRO_MANIFEST is an explicit opt-in: once set, a missing manifest or one
    // lacking the required DMD transcript is a test *failure*, not a skip —
    // otherwise this issue-#806 regression goes green without exercising #806.
    assert!(
        path.is_file(),
        "FERRO_MANIFEST points to a missing manifest: {}",
        path.display()
    );
    let provider = Arc::new(
        MultiFastaProvider::from_manifest(&path)
            .unwrap_or_else(|e| panic!("from_manifest({}) failed: {e}", path.display())),
    );
    provider
        .get_transcript("NM_004006.2")
        .expect("FERRO_MANIFEST must include NM_004006.2 (DMD) for issue #806");
    let state = real_state(provider);

    // DMD NM_004006.2: a substitution early in the CDS. We assert only that the
    // residues are REAL three-letter codes (no `?`), not a specific identity —
    // the point is that the seam reaches real sequence.
    let response = predict_effect(
        State(state),
        Json(EffectRequest {
            hgvs: "NM_004006.2:c.5T>G".to_string(),
            include_nmd: true,
        }),
    )
    .await
    .expect("effect prediction should succeed");

    let pc = response
        .0
        .protein_consequence
        .expect("a coding substitution should yield a protein consequence");

    // DMD NP_003997.1 begins Met-Leu-Trp..., so CDS codon 2 (c.4-6) is Leu
    // (CTT). `c.5T>G` rewrites the middle base CTT -> CGT = Arg — a missense.
    // Pinned against the residue identity, not just the format, so a regression
    // in the resolution seam (wrong codon offset, frame, or strand) is caught.
    // Oracle: the DMD reference protein N-terminus (M-L-W-W-...).
    assert_eq!(pc.ref_aa, "Leu", "reference residue at DMD codon 2");
    assert_eq!(pc.alt_aa, "Arg", "c.5T>G (CTT->CGT) substitutes Leu->Arg");
    assert_eq!(pc.position, 2, "protein position is codon 2");
    assert_eq!(
        pc.hgvs_p, "NP_003997.1:p.(Leu2Arg)",
        "fully-resolved HGVS-p from the engine"
    );
    assert!(!pc.is_frameshift, "a substitution is not a frameshift");

    // A missense introduces no PTC, so junction-based NMD must not predict NMD
    // (and must NOT be the old 0.9-fraction heuristic's output).
    let nmd = response
        .0
        .nmd_prediction
        .expect("NMD prediction should be produced for a real coding variant");
    assert!(!nmd.predicted, "missense introduces no PTC -> no NMD");
    assert_eq!(
        nmd.reason, "Variant does not introduce a premature termination codon",
        "NMD reason should match the no-PTC junction contract"
    );
}
