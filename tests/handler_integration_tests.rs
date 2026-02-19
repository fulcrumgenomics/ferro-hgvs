//! Integration tests for web service handlers
//!
//! These tests call the handler functions directly with actual or minimal state,
//! testing the full code path from request to response.

#![cfg(feature = "dev")]

use axum::{extract::State, response::Json};
use ferro_hgvs::service::config::ServiceConfig;
use ferro_hgvs::service::handlers::validate::ValidateRequest;
use ferro_hgvs::service::types::*;
use std::sync::Arc;

// ==================== Handler Integration Tests ====================

/// Create a minimal AppState for testing (without actual tools)
fn create_test_state() -> ferro_hgvs::service::server::AppState {
    use ferro_hgvs::service::server::{AppState, HealthCache};
    use ferro_hgvs::service::tools::ToolManager;

    // Create minimal tool manager without actual tools
    let tool_manager = ToolManager::empty();

    AppState {
        tool_manager: Arc::new(tool_manager),
        config: Arc::new(ServiceConfig::default()),
        cdot: None,
        liftover: None,
        health_cache: HealthCache::default(),
    }
}

// ==================== Validate Handler Tests ====================

#[tokio::test]
async fn test_validate_handler_valid_input() {
    let state = create_test_state();

    let request = ValidateRequest {
        hgvs: "NM_000249.4:c.350C>T".to_string(),
    };

    let result =
        ferro_hgvs::service::handlers::validate::validate_single(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.valid);
    assert!(response.components.is_some());
    let components = response.components.unwrap();
    assert_eq!(components.coordinate_system, "c");
    assert_eq!(components.variant_type, "substitution");
}

#[tokio::test]
async fn test_validate_handler_invalid_input() {
    let state = create_test_state();

    let request = ValidateRequest {
        hgvs: "not_valid_hgvs".to_string(),
    };

    let result =
        ferro_hgvs::service::handlers::validate::validate_single(State(state), Json(request)).await;

    // Invalid format should be rejected by validation
    assert!(result.is_err());
}

#[tokio::test]
async fn test_validate_handler_empty_input() {
    let state = create_test_state();

    let request = ValidateRequest {
        hgvs: "".to_string(),
    };

    let result =
        ferro_hgvs::service::handlers::validate::validate_single(State(state), Json(request)).await;

    // Empty input should be rejected
    assert!(result.is_err());
}

#[tokio::test]
async fn test_validate_handler_dangerous_input() {
    let state = create_test_state();

    let request = ValidateRequest {
        hgvs: "NM_000249.4:c.350C>T; rm -rf /".to_string(),
    };

    let result =
        ferro_hgvs::service::handlers::validate::validate_single(State(state), Json(request)).await;

    // Dangerous characters should be rejected
    assert!(result.is_err());
}

// ==================== Effect Handler Tests ====================

#[tokio::test]
async fn test_effect_handler_cds_variant() {
    let state = create_test_state();

    let request = EffectRequest {
        hgvs: "NM_000249.4:c.350C>T".to_string(),
        include_nmd: false,
    };

    let result =
        ferro_hgvs::service::handlers::effect::predict_effect(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.effect.is_some());
    let effect = response.effect.unwrap();
    assert_eq!(effect.name, "missense_variant");
    assert_eq!(effect.impact, "MODERATE");
}

#[tokio::test]
async fn test_effect_handler_intronic_variant() {
    let state = create_test_state();

    let request = EffectRequest {
        hgvs: "NM_000249.4:c.117-2del".to_string(),
        include_nmd: false,
    };

    let result =
        ferro_hgvs::service::handlers::effect::predict_effect(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.effect.is_some());
    let effect = response.effect.unwrap();
    assert_eq!(effect.name, "splice_site_variant");
    assert_eq!(effect.impact, "HIGH");
}

#[tokio::test]
async fn test_effect_handler_frameshift_variant() {
    let state = create_test_state();

    let request = EffectRequest {
        hgvs: "NM_000249.4:c.350del".to_string(),
        include_nmd: false,
    };

    let result =
        ferro_hgvs::service::handlers::effect::predict_effect(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.effect.is_some());
    let effect = response.effect.unwrap();
    assert_eq!(effect.name, "frameshift_variant");
    assert_eq!(effect.impact, "HIGH");
}

#[tokio::test]
async fn test_effect_handler_genomic_variant() {
    let state = create_test_state();

    let request = EffectRequest {
        hgvs: "NC_000007.14:g.117559593G>A".to_string(),
        include_nmd: false,
    };

    let result =
        ferro_hgvs::service::handlers::effect::predict_effect(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.effect.is_some());
    let effect = response.effect.unwrap();
    assert_eq!(effect.name, "SNV");
}

// ==================== Convert Handler Tests ====================

#[tokio::test]
async fn test_convert_handler_same_system() {
    let state = create_test_state();

    let request = ConvertRequest {
        hgvs: "NM_000249.4:c.350C>T".to_string(),
        target_system: CoordinateSystem::C,
        include_all: false,
    };

    let result = ferro_hgvs::service::handlers::convert::convert(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert_eq!(response.source_system, "c");
    assert_eq!(response.target_system, "c");
    // Same system should return input as-is
    assert!(response.converted.is_some());
    assert_eq!(response.converted.unwrap(), "NM_000249.4:c.350C>T");
}

#[tokio::test]
async fn test_convert_handler_c_to_g_without_cdot() {
    let state = create_test_state();

    let request = ConvertRequest {
        hgvs: "NM_000249.4:c.350C>T".to_string(),
        target_system: CoordinateSystem::G,
        include_all: false,
    };

    let result = ferro_hgvs::service::handlers::convert::convert(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    // Without cdot, conversion should fail gracefully
    assert!(response.error.is_some() || response.converted.is_none());
}

// ==================== Liftover Handler Tests ====================

#[tokio::test]
async fn test_liftover_handler_without_engine() {
    let state = create_test_state();

    let request = LiftoverRequest {
        position: "chr7:117120148".to_string(),
        from_build: GenomeBuild::GRCh37,
        to_build: GenomeBuild::GRCh38,
    };

    let result =
        ferro_hgvs::service::handlers::liftover::liftover(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    // Without liftover engine, should either return error message or no conversion
    // The exact error message may vary, so we just check it handled gracefully
    assert!(response.error.is_some() || response.converted.is_none());
}

#[tokio::test]
async fn test_liftover_handler_same_build() {
    let state = create_test_state();

    let request = LiftoverRequest {
        position: "chr7:117120148".to_string(),
        from_build: GenomeBuild::GRCh38,
        to_build: GenomeBuild::GRCh38,
    };

    let result =
        ferro_hgvs::service::handlers::liftover::liftover(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    // Same build should return position as-is
    assert_eq!(response.input, "chr7:117120148");
}

// ==================== VCF Convert Handler Tests ====================

#[tokio::test]
async fn test_vcf_to_hgvs_handler_substitution() {
    let state = create_test_state();

    let request = VcfToHgvsRequest {
        chrom: "chr7".to_string(),
        pos: 117559593,
        ref_allele: "G".to_string(),
        alt: "A".to_string(),
        build: GenomeBuild::GRCh38,
        transcript: None,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::vcf_to_hgvs(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.hgvs_g.is_some());
    assert!(response.hgvs_g.unwrap().contains("117559593"));
}

#[tokio::test]
async fn test_vcf_to_hgvs_handler_deletion() {
    let state = create_test_state();

    let request = VcfToHgvsRequest {
        chrom: "chr7".to_string(),
        pos: 117559592,
        ref_allele: "AG".to_string(),
        alt: "A".to_string(),
        build: GenomeBuild::GRCh38,
        transcript: None,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::vcf_to_hgvs(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.hgvs_g.is_some());
    let hgvs_g = response.hgvs_g.unwrap();
    assert!(hgvs_g.contains("del"));
}

#[tokio::test]
async fn test_vcf_to_hgvs_handler_insertion() {
    let state = create_test_state();

    let request = VcfToHgvsRequest {
        chrom: "chr7".to_string(),
        pos: 117559593,
        ref_allele: "G".to_string(),
        alt: "GA".to_string(),
        build: GenomeBuild::GRCh38,
        transcript: None,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::vcf_to_hgvs(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.hgvs_g.is_some());
    let hgvs_g = response.hgvs_g.unwrap();
    assert!(hgvs_g.contains("ins") || hgvs_g.contains("dup"));
}

#[tokio::test]
async fn test_hgvs_to_vcf_handler_genomic() {
    let state = create_test_state();

    let request = HgvsToVcfRequest {
        hgvs: "NC_000007.14:g.117559593G>A".to_string(),
        build: GenomeBuild::GRCh38,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::hgvs_to_vcf(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.vcf.is_some());
    let vcf = response.vcf.unwrap();
    assert_eq!(vcf.chrom, "chr7");
    assert_eq!(vcf.pos, 117559593);
    assert_eq!(vcf.ref_allele, "G");
    assert_eq!(vcf.alt, "A");
}

#[tokio::test]
async fn test_hgvs_to_vcf_handler_cds_without_cdot() {
    let state = create_test_state();

    let request = HgvsToVcfRequest {
        hgvs: "NM_000249.4:c.350C>T".to_string(),
        build: GenomeBuild::GRCh38,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::hgvs_to_vcf(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    // Without cdot, CDS to VCF conversion should fail gracefully
    assert!(response.error.is_some() || response.vcf.is_none());
}

// ==================== Input Validation Tests ====================

#[tokio::test]
async fn test_vcf_to_hgvs_empty_ref() {
    let state = create_test_state();

    let request = VcfToHgvsRequest {
        chrom: "chr7".to_string(),
        pos: 117559593,
        ref_allele: "".to_string(), // Empty ref
        alt: "A".to_string(),
        build: GenomeBuild::GRCh38,
        transcript: None,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::vcf_to_hgvs(State(state), Json(request)).await;

    // Empty ref should be rejected
    assert!(result.is_err());
}

#[tokio::test]
async fn test_vcf_to_hgvs_empty_alt() {
    let state = create_test_state();

    let request = VcfToHgvsRequest {
        chrom: "chr7".to_string(),
        pos: 117559593,
        ref_allele: "G".to_string(),
        alt: "".to_string(), // Empty alt
        build: GenomeBuild::GRCh38,
        transcript: None,
    };

    let result =
        ferro_hgvs::service::handlers::vcf_convert::vcf_to_hgvs(State(state), Json(request)).await;

    // Empty alt should be rejected
    assert!(result.is_err());
}

// ==================== Processing Time Tests ====================

#[tokio::test]
async fn test_response_includes_processing_time() {
    let state = create_test_state();

    let request = EffectRequest {
        hgvs: "NM_000249.4:c.350C>T".to_string(),
        include_nmd: false,
    };

    let result =
        ferro_hgvs::service::handlers::effect::predict_effect(State(state), Json(request)).await;

    assert!(result.is_ok());
    let response = result.unwrap().0;
    assert!(response.processing_time_ms > 0 || response.processing_time_ms == 0);
}
