//! Web server setup using Axum framework

use axum::{
    extract::DefaultBodyLimit,
    http::{header, StatusCode},
    response::{Html, IntoResponse, Json},
    routing::{get, post},
    Router,
};
use std::sync::Arc;
use tokio::sync::RwLock;

use crate::data::cdot::CdotMapper;
use crate::liftover::Liftover;
use crate::service::{
    config::ServiceConfig,
    handlers,
    tools::ToolManager,
    types::{DetailedHealthResponse, ErrorResponse, ServiceError},
};

/// Cached health check results
#[derive(Clone)]
pub struct HealthCache {
    /// Cached detailed health response
    pub detailed: Arc<RwLock<Option<DetailedHealthResponse>>>,
    /// Timestamp of last update
    pub last_updated: Arc<RwLock<Option<chrono::DateTime<chrono::Utc>>>>,
}

impl Default for HealthCache {
    fn default() -> Self {
        Self {
            detailed: Arc::new(RwLock::new(None)),
            last_updated: Arc::new(RwLock::new(None)),
        }
    }
}

/// Application state shared across handlers
#[derive(Clone)]
pub struct AppState {
    /// Tool manager for executing HGVS operations
    pub tool_manager: Arc<ToolManager>,
    /// Service configuration
    pub config: Arc<ServiceConfig>,
    /// Cached health check results
    pub health_cache: HealthCache,
    /// Optional transcript data for coordinate conversion
    pub cdot: Option<Arc<CdotMapper>>,
    /// Optional liftover engine for genome build conversion
    pub liftover: Option<Arc<Liftover>>,
}

/// Create the Axum application with all routes and middleware
pub fn create_app(config: ServiceConfig) -> Result<(Router, AppState), ServiceError> {
    // Initialize tool manager
    let tool_manager = Arc::new(ToolManager::new(&config)?);

    // Load optional transcript data (CdotMapper)
    let cdot = if let Some(cdot_path) = &config.data.cdot_path {
        if cdot_path.exists() {
            tracing::info!("Loading cdot transcript data from {}", cdot_path.display());
            match load_cdot(cdot_path) {
                Ok(mapper) => {
                    tracing::info!("Loaded {} transcripts from cdot", mapper.transcript_count());
                    Some(Arc::new(mapper))
                }
                Err(e) => {
                    tracing::warn!(
                        "Failed to load cdot data: {}. Coordinate conversion will be limited.",
                        e
                    );
                    None
                }
            }
        } else {
            tracing::warn!(
                "Cdot path {} does not exist. Coordinate conversion will be limited.",
                cdot_path.display()
            );
            None
        }
    } else {
        tracing::debug!("No cdot path configured. Coordinate conversion will be limited.");
        None
    };

    // Load optional liftover chain files
    let liftover = if let Some(liftover_config) = &config.data.liftover {
        if liftover_config.grch37_to_38.exists() && liftover_config.grch38_to_37.exists() {
            tracing::info!("Loading liftover chain files");
            match Liftover::from_files(&liftover_config.grch37_to_38, &liftover_config.grch38_to_37)
            {
                Ok(lo) => {
                    tracing::info!("Liftover chain files loaded successfully");
                    Some(Arc::new(lo))
                }
                Err(e) => {
                    tracing::warn!(
                        "Failed to load liftover chain files: {}. Liftover will be unavailable.",
                        e
                    );
                    None
                }
            }
        } else {
            tracing::warn!("Liftover chain files not found. Liftover will be unavailable.");
            None
        }
    } else {
        tracing::debug!("No liftover chain files configured. Liftover will be unavailable.");
        None
    };

    // Create shared application state
    let state = AppState {
        tool_manager,
        config: Arc::new(config.clone()),
        health_cache: HealthCache::default(),
        cdot,
        liftover,
    };

    // Parse max request size
    let max_size = parse_size(&config.server.max_request_size)
        .map_err(|e| ServiceError::ConfigError(format!("Invalid max_request_size: {}", e)))?;

    // Build router
    let mut app = Router::new()
        // Web UI routes
        .route("/", get(index_handler))
        .route("/static/css/styles.css", get(styles_css_handler))
        .route("/static/js/main.js", get(main_js_handler))
        // Health endpoints
        .route("/health", get(handlers::health::health_check))
        .route(
            "/health/detailed",
            get(handlers::health::detailed_health_check),
        )
        .route("/api/v1/health", get(handlers::health::health_check))
        .route(
            "/api/v1/health/detailed",
            get(handlers::health::detailed_health_check),
        )
        .route("/api/v1/tools/status", get(handlers::health::tools_status))
        // Single variant endpoints
        .route("/api/v1/parse", post(handlers::parse::parse_single))
        .route(
            "/api/v1/validate",
            post(handlers::validate::validate_single),
        )
        .route(
            "/api/v1/normalize",
            post(handlers::normalize::normalize_single),
        )
        // Batch endpoints
        .route("/api/v1/batch/parse", post(handlers::parse::parse_batch))
        .route(
            "/api/v1/batch/normalize",
            post(handlers::normalize::normalize_batch),
        )
        // Coordinate conversion endpoint
        .route("/api/v1/convert", post(handlers::convert::convert))
        // Effect prediction endpoint
        .route("/api/v1/effect", post(handlers::effect::predict_effect))
        // Liftover endpoint
        .route("/api/v1/liftover", post(handlers::liftover::liftover))
        // VCF conversion endpoints
        .route(
            "/api/v1/vcf-to-hgvs",
            post(handlers::vcf_convert::vcf_to_hgvs),
        )
        .route(
            "/api/v1/hgvs-to-vcf",
            post(handlers::vcf_convert::hgvs_to_vcf),
        )
        // API info endpoint (separate from web UI)
        .route("/api/v1/info", get(handlers::info::service_info))
        // Handle 404s
        .fallback(handle_404)
        // Add shared state
        .with_state(state.clone());

    // Add middleware layers (simplified for now)
    app = app.layer(DefaultBodyLimit::max(max_size));

    Ok((app, state))
}

/// Handle 404 errors
async fn handle_404() -> (StatusCode, Json<ErrorResponse>) {
    let error = ServiceError::BadRequest("Endpoint not found".to_string());
    (StatusCode::NOT_FOUND, Json(error.to_response()))
}

/// Serve the main web UI page
async fn index_handler() -> Html<&'static str> {
    Html(include_str!("web/templates/index.html"))
}

/// Serve CSS stylesheet
async fn styles_css_handler() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/css; charset=utf-8")],
        include_str!("web/static/css/styles.css"),
    )
}

/// Serve JavaScript
async fn main_js_handler() -> impl IntoResponse {
    (
        [(
            header::CONTENT_TYPE,
            "application/javascript; charset=utf-8",
        )],
        include_str!("web/static/js/main.js"),
    )
}

/// Parse size strings like "10MB", "1GB", etc.
fn parse_size(size_str: &str) -> Result<usize, String> {
    let size_str = size_str.to_uppercase();

    // Check longer suffixes first to avoid partial matches
    if let Some(num_str) = size_str.strip_suffix("GB") {
        let num: usize = num_str
            .parse()
            .map_err(|_| format!("Invalid size format: {}", size_str))?;
        return Ok(num * 1024 * 1024 * 1024);
    }

    if let Some(num_str) = size_str.strip_suffix("MB") {
        let num: usize = num_str
            .parse()
            .map_err(|_| format!("Invalid size format: {}", size_str))?;
        return Ok(num * 1024 * 1024);
    }

    if let Some(num_str) = size_str.strip_suffix("KB") {
        let num: usize = num_str
            .parse()
            .map_err(|_| format!("Invalid size format: {}", size_str))?;
        return Ok(num * 1024);
    }

    if let Some(num_str) = size_str.strip_suffix("B") {
        return num_str
            .parse::<usize>()
            .map_err(|_| format!("Invalid size format: {}", size_str));
    }

    // Try parsing as plain number (bytes)
    size_str
        .parse::<usize>()
        .map_err(|_| format!("Invalid size format: {}", size_str))
}

/// Health check interval (15 minutes)
const HEALTH_CHECK_INTERVAL_SECS: u64 = 15 * 60;

/// Spawn background task to periodically update health cache
pub fn spawn_health_check_task(state: AppState) {
    tokio::spawn(async move {
        // Run initial health check immediately
        update_health_cache(&state).await;

        // Then run periodically
        let mut interval =
            tokio::time::interval(std::time::Duration::from_secs(HEALTH_CHECK_INTERVAL_SECS));
        interval.tick().await; // First tick is immediate, skip it since we just ran

        loop {
            interval.tick().await;
            tracing::info!("Running periodic health check...");
            update_health_cache(&state).await;
        }
    });
}

/// Update the health cache with fresh results
async fn update_health_cache(state: &AppState) {
    match handlers::health::run_detailed_health_check(state).await {
        Ok(response) => {
            let now = chrono::Utc::now();
            *state.health_cache.detailed.write().await = Some(response);
            *state.health_cache.last_updated.write().await = Some(now);
            tracing::info!("Health cache updated at {}", now);
        }
        Err(e) => {
            tracing::error!("Failed to update health cache: {:?}", e);
        }
    }
}

/// Load CdotMapper from a file path (handles both .json and .json.gz)
fn load_cdot(path: &std::path::Path) -> Result<CdotMapper, ServiceError> {
    let path_str = path.to_string_lossy();
    if path_str.ends_with(".gz") {
        CdotMapper::from_json_gz(path)
            .map_err(|e| ServiceError::ConfigError(format!("Failed to load cdot: {}", e)))
    } else {
        CdotMapper::from_json_file(path)
            .map_err(|e| ServiceError::ConfigError(format!("Failed to load cdot: {}", e)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_size() {
        assert_eq!(parse_size("100").unwrap(), 100);
        assert_eq!(parse_size("100B").unwrap(), 100);
        assert_eq!(parse_size("1KB").unwrap(), 1024);
        assert_eq!(parse_size("10MB").unwrap(), 10 * 1024 * 1024);
        assert_eq!(parse_size("1GB").unwrap(), 1024 * 1024 * 1024);

        // Case insensitive
        assert_eq!(parse_size("10mb").unwrap(), 10 * 1024 * 1024);

        // Invalid formats
        assert!(parse_size("invalid").is_err());
        assert!(parse_size("10XB").is_err());
    }
}
