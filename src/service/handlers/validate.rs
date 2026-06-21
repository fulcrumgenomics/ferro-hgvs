//! HGVS validation endpoint - validates syntax and shows component breakdown

use axum::{extract::State, http::StatusCode, response::Json};
use std::time::Instant;

use crate::hgvs::interval::interval_is_wraparound;
use crate::service::{
    server::AppState,
    types::{extract_variant_details, ErrorResponse, ServiceError, ValidateResponse},
    validation::validate_hgvs as security_validate_hgvs,
};

/// Request for single variant validation
#[derive(Debug, serde::Deserialize)]
pub struct ValidateRequest {
    /// The HGVS variant string to validate
    pub hgvs: String,
}

/// Validate a single HGVS variant and return component breakdown
///
/// This endpoint parses the HGVS string and returns:
/// - Whether the syntax is valid
/// - Any validation errors or warnings
/// - A breakdown of the variant components (reference, position, edit type, etc.)
pub async fn validate_single(
    State(_state): State<AppState>,
    Json(request): Json<ValidateRequest>,
) -> Result<Json<ValidateResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // First, do security validation on the input
    if let Err(validation_error) = security_validate_hgvs(&request.hgvs) {
        let error = ServiceError::InvalidHgvs(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Parse using ferro's lenient parser on a blocking task so it does not run
    // on the async runtime thread.
    let hgvs_str = request.hgvs.clone();
    let parse_result =
        tokio::task::spawn_blocking(move || crate::hgvs::parser::parse_hgvs_lenient(&hgvs_str))
            .await
            .map_err(|e| {
                let error = ServiceError::InternalError(format!("Task error: {}", e));
                (StatusCode::INTERNAL_SERVER_ERROR, Json(error.to_response()))
            })?;

    let elapsed_ms = start.elapsed().as_millis() as u64;

    match parse_result {
        Ok(result) => {
            // Extract component details
            let components = extract_variant_details(&result.result);
            let wraps_origin = variant_wraps_origin(&result.result);

            // Collect any warnings as potential issues
            let errors = if result.warnings.is_empty() {
                None
            } else {
                Some(result.warnings.iter().map(|w| w.message.clone()).collect())
            };

            Ok(Json(ValidateResponse {
                input: request.hgvs,
                valid: true,
                errors,
                components,
                processing_time_ms: elapsed_ms,
                wraps_origin,
            }))
        }
        Err(e) => Ok(Json(ValidateResponse {
            input: request.hgvs,
            valid: false,
            errors: Some(vec![e.to_string()]),
            components: None,
            processing_time_ms: elapsed_ms,
            wraps_origin: false,
        })),
    }
}

/// Return `true` if the variant is a wraparound range on a circular contig
/// (i.e. an `m.` or `o.` variant whose start position exceeds its end position
/// per SVD-WG006). Always returns `false` for linear coordinate systems.
fn variant_wraps_origin(variant: &crate::hgvs::variant::HgvsVariant) -> bool {
    use crate::hgvs::variant::HgvsVariant;
    match variant {
        HgvsVariant::Mt(v) => interval_is_wraparound(&v.loc_edit.location),
        HgvsVariant::Circular(v) => interval_is_wraparound(&v.loc_edit.location),
        _ => false,
    }
}

/// Validate an HGVS string and return a [`ValidateResponse`] with component
/// breakdown and the `wraps_origin` flag.
///
/// Unlike the async `validate_single` axum handler, this function is
/// synchronous and suitable for unit tests and non-HTTP callers.
///
/// `wraps_origin` is derived from the lenient parse: `detect_swapped_positions`
/// is axis-aware and leaves `m.`/`o.` wraparound ranges untouched, so the
/// reversed `<high>_<low>` form that marks a circular-contig wraparound
/// survives preprocessing intact.
pub fn validate_hgvs(hgvs: &str) -> ValidateResponse {
    match crate::hgvs::parser::parse_hgvs_lenient(hgvs) {
        Ok(result) => {
            let components = extract_variant_details(&result.result);
            let wraps_origin = variant_wraps_origin(&result.result);
            let errors = if result.warnings.is_empty() {
                None
            } else {
                Some(result.warnings.iter().map(|w| w.message.clone()).collect())
            };
            ValidateResponse {
                input: hgvs.to_string(),
                valid: true,
                errors,
                components,
                processing_time_ms: 0,
                wraps_origin,
            }
        }
        Err(e) => ValidateResponse {
            input: hgvs.to_string(),
            valid: false,
            errors: Some(vec![e.to_string()]),
            components: None,
            processing_time_ms: 0,
            wraps_origin: false,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Per-arm coverage of the shared `extract_variant_details` /
    // `extract_na_edit_info` helpers lives in `crate::service::types`. Here we
    // cover the validate handler's own surface: that `validate_hgvs` wires the
    // shared component breakdown into the response and that `variant_wraps_origin`
    // is axis-aware.

    #[test]
    fn test_validate_hgvs_populates_components() {
        let response = validate_hgvs("NM_000249.4:c.350C>T");
        assert!(response.valid);
        let components = response.components.expect("components populated");
        assert_eq!(components.coordinate_system, "c");
        assert_eq!(components.variant_type, "substitution");
        assert_eq!(components.reference, "NM_000249.4");
    }

    #[test]
    fn test_validate_hgvs_allele_has_no_components() {
        // Compound alleles have no single coordinate/position to flatten, so the
        // shared extractor returns None and the response carries no breakdown.
        let response = validate_hgvs("NM_000088.3:c.[10A>G;20C>T]");
        assert!(response.valid);
        assert!(response.components.is_none());
    }

    #[test]
    fn test_validate_hgvs_invalid_input() {
        let response = validate_hgvs("not a variant");
        assert!(!response.valid);
        assert!(response.components.is_none());
        assert!(!response.wraps_origin);
    }

    #[test]
    fn test_variant_wraps_origin_linear_axis_is_false() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350C>T").unwrap();
        assert!(!variant_wraps_origin(&result.result));
    }
}
