//! Biocommons/HGVS tool service implementation

use url::Url;

use crate::benchmark::{
    biocommons::{normalize_single, BiocommonsLocalConfig},
    types::ParseResult,
};
use crate::service::{
    config::BiocommonsConfig,
    tools::HgvsToolService,
    types::{analyze_error_structured, health_check::HealthCheckResult, ServiceError, ToolName},
};

/// Biocommons/HGVS tool service
pub struct BiocommonsService {
    /// Configuration (stored for future use)
    _config: BiocommonsConfig,
    /// Runtime configuration for biocommons
    biocommons_config: BiocommonsLocalConfig,
    /// UTA schema name
    uta_schema: String,
}

impl BiocommonsService {
    /// Create a new BiocommonsService
    pub fn new(config: &BiocommonsConfig) -> Result<Self, ServiceError> {
        // Convert service config to biocommons local config
        let biocommons_config = BiocommonsLocalConfig {
            uta_container_name: config
                .docker_container
                .clone()
                .unwrap_or_else(|| "ferro-uta".to_string()),
            uta_image_tag: "uta_20210129b".to_string(), // Default UTA image
            uta_port: extract_port_from_url(&config.uta_url)?,
            seqrepo_dir: config.seqrepo_path.clone(),
            seqrepo_instance: "2021-01-29".to_string(), // Default instance
        };

        // Validate paths exist
        if !config.seqrepo_path.exists() {
            return Err(ServiceError::ConfigError(format!(
                "SeqRepo path does not exist: {}",
                config.seqrepo_path.display()
            )));
        }

        Ok(Self {
            _config: config.clone(),
            biocommons_config,
            uta_schema: config.uta_schema.clone(),
        })
    }

    /// Run biocommons normalization on a single variant
    async fn run_biocommons(
        &self,
        hgvs: &str,
        _is_normalize: bool,
    ) -> Result<ParseResult, ServiceError> {
        let hgvs = hgvs.to_string();
        let config = self.biocommons_config.clone();
        let uta_schema = self.uta_schema.clone();

        // Run biocommons in a blocking task since it uses subprocesses
        tokio::task::spawn_blocking(move || {
            // Construct UTA database URL using the configured schema
            let uta_db_url = format!(
                "postgresql://anonymous:anonymous@localhost:{}/uta/{}",
                config.uta_port, uta_schema
            );

            // Use the normalize_single function from the benchmark module
            let result = normalize_single(
                &hgvs,
                Some(&uta_db_url),
                Some(&config.seqrepo_dir.to_string_lossy()),
                None, // LRG mapping file not supported yet
            );

            match result {
                Ok(parse_result) => Ok(parse_result),
                Err(e) => {
                    // Convert FerroError to a ParseResult with error
                    Ok(ParseResult {
                        input: hgvs,
                        success: false,
                        output: None,
                        error: Some(e.to_string()),
                        error_category: Some("biocommons_error".to_string()),
                        ref_mismatch: None,
                        details: None,
                    })
                }
            }
        })
        .await
        .map_err(|e| ServiceError::InternalError(format!("Task join error: {}", e)))?
    }
}

#[async_trait::async_trait]
impl HgvsToolService for BiocommonsService {
    async fn parse(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        // For biocommons, parsing and normalization use the same interface
        // The tool does both parsing and basic validation
        self.run_biocommons(hgvs, false).await
    }

    async fn normalize(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        self.run_biocommons(hgvs, true).await
    }

    async fn health_check(&self) -> HealthCheckResult {
        // Use a well-known transcript that's likely to exist in most UTA databases
        let test_variant = "NM_000088.3:c.589G>T";

        match self.run_biocommons(test_variant, false).await {
            Ok(result) => {
                if result.success {
                    // Tool processed variant successfully
                    HealthCheckResult::Healthy
                } else if let Some(error) = &result.error {
                    // Analyze error using structured error categorization
                    let error_category = analyze_error_structured(ToolName::Biocommons, error);

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
        ToolName::Biocommons
    }
}

/// Extract port number from a PostgreSQL URL with security validation
fn extract_port_from_url(url: &str) -> Result<u16, ServiceError> {
    // Use proper URL parsing instead of naive string splitting
    let parsed = Url::parse(url)
        .map_err(|e| ServiceError::ConfigError(format!("Invalid URL format '{}': {}", url, e)))?;

    // Security: validate scheme - only allow postgresql connections
    if parsed.scheme() != "postgresql" {
        return Err(ServiceError::ConfigError(format!(
            "Only postgresql URLs are allowed, got scheme: {}",
            parsed.scheme()
        )));
    }

    // Security: validate host - only allow localhost for security
    let host = parsed.host_str().unwrap_or("");
    if host != "localhost" && host != "127.0.0.1" {
        return Err(ServiceError::ConfigError(format!(
            "Only localhost connections are allowed for security, got host: {}",
            host
        )));
    }

    // Extract and validate port
    let port = parsed.port().unwrap_or(5432);

    // Security: ensure port is in valid range (not system ports)
    if port < 1024 {
        return Err(ServiceError::ConfigError(format!(
            "Port must be >= 1024 for security, got: {}",
            port
        )));
    }

    Ok(port)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_extract_port_from_url() {
        // Valid localhost URLs should work
        let url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b";
        assert_eq!(extract_port_from_url(url).unwrap(), 5432);

        let url2 = "postgresql://user:pass@127.0.0.1:1234/db";
        assert_eq!(extract_port_from_url(url2).unwrap(), 1234);

        // Default port should be used when not specified
        let url3 = "postgresql://user:pass@localhost/db";
        assert_eq!(extract_port_from_url(url3).unwrap(), 5432);
    }

    #[test]
    fn test_extract_port_from_url_security_validation() {
        // Non-postgresql schemes should be rejected
        let url = "mysql://user:pass@localhost:3306/db";
        assert!(extract_port_from_url(url).is_err());

        // Non-localhost hosts should be rejected for security
        let url2 = "postgresql://user:pass@example.com:5432/db";
        assert!(extract_port_from_url(url2).is_err());

        // System ports should be rejected for security
        let url3 = "postgresql://user:pass@localhost:22/db";
        assert!(extract_port_from_url(url3).is_err());

        // Invalid URLs should be rejected
        let url4 = "not-a-valid-url";
        assert!(extract_port_from_url(url4).is_err());
    }

    #[test]
    fn test_biocommons_service_creation() {
        let config = BiocommonsConfig {
            enabled: true,
            uta_url: "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b"
                .to_string(),
            uta_schema: "uta_20210129b".to_string(),
            seqrepo_path: PathBuf::from("/nonexistent/seqrepo/path"), // Non-existent path for test
            docker_container: Some("test-uta".to_string()),
            parallel_workers: Some(1),
            env_vars: None,
        };

        // This will fail because the seqrepo path does not exist
        let result = BiocommonsService::new(&config);
        assert!(result.is_err());
    }
}
