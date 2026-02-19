//! Service information endpoints

use axum::{extract::State, response::Json};
use serde_json::{json, Value};

use crate::service::server::AppState;

/// Provide basic service information
pub async fn service_info(State(state): State<AppState>) -> Json<Value> {
    let enabled_tools = state.config.enabled_tools();

    Json(json!({
        "service": "ferro-hgvs-web",
        "version": env!("CARGO_PKG_VERSION"),
        "description": "Multi-tool HGVS variant normalization web service",
        "available_tools": enabled_tools,
        "endpoints": {
            "parse": {
                "single": "POST /api/v1/parse",
                "batch": "POST /api/v1/batch/parse"
            },
            "normalize": {
                "single": "POST /api/v1/normalize",
                "batch": "POST /api/v1/batch/normalize"
            },
            "health": {
                "service": "GET /api/v1/health",
                "tools": "GET /api/v1/tools/status"
            }
        },
        "documentation": {
            "repository": "https://github.com/fulcrumgenomics/ferro-hgvs",
            "api_docs": "/api/v1/info"
        }
    }))
}
