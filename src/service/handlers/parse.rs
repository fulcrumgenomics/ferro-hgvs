//! HGVS parsing endpoints

use axum::{extract::State, http::StatusCode, response::Json};

use crate::service::{
    server::AppState,
    types::{
        BatchRequest, BatchResponse, ErrorResponse, ServiceError, SingleRequest, SingleResponse,
    },
    validation::{validate_hgvs, validate_hgvs_batch, validate_optional_timeout},
};

/// Parse a single HGVS variant
pub async fn parse_single(
    State(state): State<AppState>,
    Json(request): Json<SingleRequest>,
) -> Result<Json<SingleResponse>, (StatusCode, Json<ErrorResponse>)> {
    // Comprehensive input validation for security
    if let Err(validation_error) = validate_hgvs(&request.hgvs) {
        let error = ServiceError::InvalidHgvs(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Validate timeout if provided
    if let Err(validation_error) = validate_optional_timeout(request.timeout_seconds) {
        let error = ServiceError::BadRequest(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Convert tools to vector of strings for ToolManager (which still expects strings)
    let tools_as_strings = request.tools.as_ref().map(|tools| {
        tools
            .iter()
            .map(|tool| tool.to_string())
            .collect::<Vec<String>>()
    });

    // Parse using tool manager
    match state
        .tool_manager
        .parse_single(
            &request.hgvs,
            tools_as_strings.as_deref(),
            request.timeout_seconds,
        )
        .await
    {
        Ok(response) => Ok(Json(response)),
        Err(e) => {
            let status_code =
                StatusCode::from_u16(e.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR);
            Err((status_code, Json(e.to_response())))
        }
    }
}

/// Parse multiple HGVS variants
pub async fn parse_batch(
    State(state): State<AppState>,
    Json(request): Json<BatchRequest>,
) -> Result<Json<BatchResponse>, (StatusCode, Json<ErrorResponse>)> {
    // Comprehensive batch validation for security
    if let Err(validation_error) = validate_hgvs_batch(&request.variants) {
        let error = ServiceError::BadRequest(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Validate timeout if provided
    if let Err(validation_error) = validate_optional_timeout(request.timeout_seconds) {
        let error = ServiceError::BadRequest(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Convert tools to vector of strings for ToolManager (which still expects strings)
    let tools_as_strings = request.tools.as_ref().map(|tools| {
        tools
            .iter()
            .map(|tool| tool.to_string())
            .collect::<Vec<String>>()
    });

    // Parse using tool manager
    match state
        .tool_manager
        .parse_batch(
            &request.variants,
            tools_as_strings.as_deref(),
            request.timeout_seconds,
        )
        .await
    {
        Ok(response) => Ok(Json(response)),
        Err(e) => {
            let status_code =
                StatusCode::from_u16(e.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR);
            Err((status_code, Json(e.to_response())))
        }
    }
}
