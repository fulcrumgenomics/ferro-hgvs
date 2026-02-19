//! Mutalyzer tool service implementation

use std::sync::Arc;

use crate::benchmark::types::ParseResult;
use crate::service::{
    config::{MutalyzerConfig, MutalyzerMode},
    tools::{http_client::EnhancedHttpClient, HgvsToolService},
    types::{analyze_error_structured, health_check::HealthCheckResult, ServiceError, ToolName},
};

use serde::Deserialize;

/// Response from Mutalyzer normalize endpoint
#[derive(Debug, Deserialize)]
pub struct NormalizeResponse {
    /// The normalized HGVS description
    pub normalized_description: Option<String>,

    /// Original input
    pub input_description: Option<String>,

    /// Errors if any
    #[serde(default)]
    pub errors: Vec<MutalyzerError>,

    /// Warnings if any
    #[serde(default)]
    pub warnings: Vec<MutalyzerWarning>,

    /// Custom field containing detailed error info (returned on 422 errors)
    pub custom: Option<MutalyzerCustom>,

    /// Top-level message (e.g., "Errors encountered. Check the 'custom' field.")
    pub message: Option<String>,
}

/// Custom field in mutalyzer response containing detailed error information
#[derive(Debug, Deserialize)]
pub struct MutalyzerCustom {
    /// Detailed error messages
    #[serde(default)]
    pub errors: Vec<MutalyzerDetailedError>,
}

/// Detailed error from mutalyzer custom field
#[derive(Debug, Deserialize)]
pub struct MutalyzerDetailedError {
    pub code: Option<String>,
    pub details: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct MutalyzerError {
    pub code: Option<String>,
    pub message: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct MutalyzerWarning {
    pub code: Option<String>,
    pub message: Option<String>,
}

/// Mutalyzer tool service
pub struct MutalyzerService {
    /// Enhanced HTTP client with connection pooling and rate limiting (only used in API mode)
    client: Option<Arc<EnhancedHttpClient>>,
    /// Base URL for Mutalyzer API (only used in API mode)
    base_url: String,
    /// Configuration
    config: MutalyzerConfig,
}

impl MutalyzerService {
    /// Create a new MutalyzerService
    pub fn new(config: &MutalyzerConfig) -> Result<Self, ServiceError> {
        // Only create HTTP client for API mode
        let client = match config.mode {
            MutalyzerMode::Api => {
                let http_client = EnhancedHttpClient::new(
                    config.connection_pool.as_ref(),
                    config.circuit_breaker.as_ref(),
                    config.rate_limit_ms,
                )?;
                Some(Arc::new(http_client))
            }
            MutalyzerMode::Local => {
                // Validate that mutalyzer Python package is available
                if !crate::benchmark::mutalyzer::has_mutalyzer_normalizer() {
                    return Err(ServiceError::ConfigError(
                        "Mutalyzer local mode requires Python 'mutalyzer' package. Install with: pip install mutalyzer".to_string()
                    ));
                }
                None
            }
        };

        Ok(Self {
            client,
            base_url: config.api_url.trim_end_matches('/').to_string(),
            config: config.clone(),
        })
    }

    /// Normalize a single HGVS expression using the HTTP API
    async fn normalize_via_api(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        let client = self.client.as_ref().ok_or_else(|| {
            ServiceError::InternalError("HTTP client not available in local mode".to_string())
        })?;

        let encoded = urlencoding::encode(hgvs);
        // mutalyzer.nl uses /api/normalize/, local instances use /normalize/
        let url = if self.base_url.contains("mutalyzer.nl") {
            format!("{}/api/normalize/{}", self.base_url, encoded)
        } else {
            format!("{}/normalize/{}", self.base_url, encoded)
        };

        let response = client.get(&url).await?;
        let status = response.status();

        // Try to parse response body even for non-success status codes
        // Mutalyzer returns JSON with error details in 422 responses
        let body: NormalizeResponse = match response.json().await {
            Ok(body) => body,
            Err(e) => {
                // If we can't parse the body, return the HTTP status
                if !status.is_success() {
                    return Ok(ParseResult {
                        input: hgvs.to_string(),
                        success: false,
                        output: None,
                        error: Some(format!("HTTP {}", status)),
                        error_category: Some("http_error".to_string()),
                        ref_mismatch: None,
                        details: None,
                    });
                }
                return Err(ServiceError::InternalError(format!(
                    "Failed to parse response: {}",
                    e
                )));
            }
        };

        // Check for errors in the custom field (returned on 422 errors)
        if let Some(custom) = &body.custom {
            if !custom.errors.is_empty() {
                let error_msgs: Vec<String> = custom
                    .errors
                    .iter()
                    .filter_map(|e| match (&e.code, &e.details) {
                        (Some(code), Some(details)) => Some(format!("{}: {}", code, details)),
                        (Some(code), None) => Some(code.clone()),
                        (None, Some(details)) => Some(details.clone()),
                        (None, None) => None,
                    })
                    .collect();
                if !error_msgs.is_empty() {
                    return Ok(ParseResult {
                        input: hgvs.to_string(),
                        success: false,
                        output: None,
                        error: Some(error_msgs.join("; ")),
                        error_category: Some("mutalyzer_error".to_string()),
                        ref_mismatch: None,
                        details: None,
                    });
                }
            }
        }

        // Check for HTTP error status that wasn't handled by custom errors
        if !status.is_success() {
            let error_msg = body.message.unwrap_or_else(|| format!("HTTP {}", status));
            return Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: Some(error_msg),
                error_category: Some("http_error".to_string()),
                ref_mismatch: None,
                details: None,
            });
        }

        if !body.errors.is_empty() {
            let error_msg = body
                .errors
                .iter()
                .filter_map(|e| e.message.as_ref())
                .cloned()
                .collect::<Vec<_>>()
                .join("; ");
            return Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: Some(error_msg),
                error_category: Some("mutalyzer_error".to_string()),
                ref_mismatch: None,
                details: None,
            });
        }

        match body.normalized_description {
            Some(normalized) => Ok(ParseResult {
                input: hgvs.to_string(),
                success: true,
                output: Some(normalized),
                error: None,
                error_category: None,
                ref_mismatch: None,
                details: None,
            }),
            None => Ok(ParseResult {
                input: hgvs.to_string(),
                success: false,
                output: None,
                error: Some("No normalized description returned".to_string()),
                error_category: Some("empty_response".to_string()),
                ref_mismatch: None,
                details: None,
            }),
        }
    }

    /// Normalize a single HGVS expression using local Python subprocess
    async fn normalize_via_local(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        let hgvs = hgvs.to_string();
        let settings_file = self.config.settings_file.clone();
        let allow_network = self.config.allow_network;

        tokio::task::spawn_blocking(move || {
            crate::benchmark::mutalyzer::normalize_single(
                &hgvs,
                settings_file
                    .as_ref()
                    .map(|p| p.to_string_lossy())
                    .as_deref(),
                allow_network,
            )
        })
        .await
        .map_err(|e| ServiceError::InternalError(format!("Task join error: {}", e)))?
        .map_err(|e| ServiceError::InternalError(format!("Mutalyzer error: {}", e)))
    }

    /// Normalize a single HGVS expression (dispatches based on mode)
    async fn normalize_variant(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        match self.config.mode {
            MutalyzerMode::Api => self.normalize_via_api(hgvs).await,
            MutalyzerMode::Local => self.normalize_via_local(hgvs).await,
        }
    }

    /// Health check for local mode
    async fn health_check_local(&self) -> HealthCheckResult {
        // For local mode, verify the Python package is available
        if !crate::benchmark::mutalyzer::has_mutalyzer_normalizer() {
            return HealthCheckResult::Unhealthy {
                reason: "Python 'mutalyzer' package not available".to_string(),
            };
        }

        // Try a simple normalization to verify it works
        match self.normalize_via_local("NM_000088.3:c.589G>T").await {
            Ok(result) => {
                if result.success {
                    HealthCheckResult::Healthy
                } else if let Some(error) = &result.error {
                    // Check if it's a cache miss (expected in offline mode without network)
                    if error.contains("Network access") || error.contains("cache") {
                        HealthCheckResult::Degraded {
                            reason:
                                "Local mode working but may need network or cache for some variants"
                                    .to_string(),
                        }
                    } else {
                        HealthCheckResult::Degraded {
                            reason: format!(
                                "Local mode available but test variant failed: {}",
                                error
                            ),
                        }
                    }
                } else {
                    HealthCheckResult::Degraded {
                        reason: "Local mode available but test variant returned no output"
                            .to_string(),
                    }
                }
            }
            Err(e) => HealthCheckResult::Unhealthy {
                reason: format!("Local mode error: {}", e),
            },
        }
    }
}

#[async_trait::async_trait]
impl HgvsToolService for MutalyzerService {
    async fn parse(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        // For mutalyzer, parsing and normalization use the same interface
        self.normalize_variant(hgvs).await
    }

    async fn normalize(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        self.normalize_variant(hgvs).await
    }

    async fn health_check(&self) -> HealthCheckResult {
        // Dispatch to mode-specific health check
        match self.config.mode {
            MutalyzerMode::Local => return self.health_check_local().await,
            MutalyzerMode::Api => {
                // Continue with API health check below
            }
        }

        // API mode health check using a known-valid variant
        let test_variant = "NM_000088.3:c.589G>T";

        match self.normalize(test_variant).await {
            Ok(result) => {
                if result.success {
                    // Tool processed variant successfully
                    HealthCheckResult::Healthy
                } else if let Some(error) = &result.error {
                    // Analyze error using structured error categorization
                    let error_category = analyze_error_structured(ToolName::Mutalyzer, error);

                    match error_category {
                        // These are acceptable failures that indicate the tool is working
                        crate::service::types::StructuredErrorCategory::Reference(_) => {
                            HealthCheckResult::Degraded {
                                reason: "Reference database issues, but tool is responding"
                                    .to_string(),
                            }
                        }
                        crate::service::types::StructuredErrorCategory::Validation(_) => {
                            HealthCheckResult::Degraded {
                                reason: "Tool validation rules active, service operational"
                                    .to_string(),
                            }
                        }
                        // Parse errors or tool errors indicate more serious issues
                        crate::service::types::StructuredErrorCategory::Parse(_)
                        | crate::service::types::StructuredErrorCategory::Tool(_)
                        | crate::service::types::StructuredErrorCategory::Internal => {
                            HealthCheckResult::Unhealthy {
                                reason: format!("Tool error: {}", error),
                            }
                        }
                        crate::service::types::StructuredErrorCategory::Timeout => {
                            HealthCheckResult::Unhealthy {
                                reason: "Tool not responding within timeout".to_string(),
                            }
                        }
                    }
                } else {
                    HealthCheckResult::Unhealthy {
                        reason: "Tool failed without error message".to_string(),
                    }
                }
            }
            Err(e) => HealthCheckResult::Unhealthy {
                reason: format!("Service error: {}", e),
            },
        }
    }

    fn tool_name(&self) -> ToolName {
        ToolName::Mutalyzer
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mutalyzer_service_creation_api_mode() {
        let config = MutalyzerConfig {
            enabled: true,
            mode: MutalyzerMode::Api,
            api_url: "http://localhost:8082".to_string(),
            timeout_seconds: 30,
            rate_limit_ms: Some(50),
            health_check_interval: Some(60),
            connection_pool: None,
            circuit_breaker: None,
            settings_file: None,
            allow_network: false,
        };

        let service = MutalyzerService::new(&config).expect("Failed to create service");
        assert_eq!(service.tool_name(), ToolName::Mutalyzer);
        assert_eq!(service.config.api_url, "http://localhost:8082");
        assert_eq!(service.config.mode, MutalyzerMode::Api);
        assert!(service.client.is_some()); // API mode should have HTTP client
    }

    #[tokio::test]
    async fn test_mutalyzer_health_check_timeout() {
        let config = MutalyzerConfig {
            enabled: true,
            mode: MutalyzerMode::Api,
            api_url: "http://nonexistent.invalid:12345".to_string(), // This should fail
            timeout_seconds: 1,                                      // Short timeout
            rate_limit_ms: Some(50),
            health_check_interval: Some(60),
            connection_pool: None,
            circuit_breaker: None,
            settings_file: None,
            allow_network: false,
        };

        let service = MutalyzerService::new(&config).expect("Failed to create service");
        let result = service.health_check().await;

        // Should return unhealthy status
        match result {
            HealthCheckResult::Unhealthy { .. } => (),
            other => panic!("Expected Unhealthy health check result, got: {:?}", other),
        }
    }

    #[test]
    fn test_mutalyzer_mode_default() {
        let mode = MutalyzerMode::default();
        assert_eq!(mode, MutalyzerMode::Api);
    }
}
