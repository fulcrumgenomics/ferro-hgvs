//! Configuration for the HGVS web service

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

/// Main service configuration
#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct ServiceConfig {
    /// Server configuration
    pub server: ServerConfig,
    /// Tool configurations
    pub tools: ToolConfigs,
    /// Data source configurations
    #[serde(default)]
    pub data: DataConfig,
}

/// Data source configuration for advanced features
#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct DataConfig {
    /// Path to cdot transcript JSON file (for coordinate conversion)
    pub cdot_path: Option<PathBuf>,
    /// Liftover chain file configuration
    pub liftover: Option<LiftoverConfig>,
}

/// Liftover chain file configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct LiftoverConfig {
    /// Path to hg19ToHg38.over.chain.gz (GRCh37 to GRCh38)
    pub grch37_to_38: PathBuf,
    /// Path to hg38ToHg19.over.chain.gz (GRCh38 to GRCh37)
    pub grch38_to_37: PathBuf,
}

/// Server configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ServerConfig {
    /// Host to bind to (default: "0.0.0.0")
    pub host: String,
    /// Port to listen on (default: 3000)
    pub port: u16,
    /// Maximum request size (default: "10MB")
    pub max_request_size: String,
    /// Request timeout in seconds (default: 60)
    pub request_timeout_seconds: u64,
    /// Enable CORS (default: true)
    pub enable_cors: bool,
    /// Enable request tracing (default: true)
    pub enable_tracing: bool,
    /// Maximum concurrent batch processing tasks (default: 10)
    pub max_concurrent_batches: Option<usize>,
    /// Maximum variants per batch request (default: 1000)
    pub max_batch_size: Option<usize>,
}

/// Configuration for all tools
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ToolConfigs {
    /// Ferro tool configuration
    pub ferro: Option<FerroConfig>,
    /// Mutalyzer tool configuration
    pub mutalyzer: Option<MutalyzerConfig>,
    /// Biocommons tool configuration
    pub biocommons: Option<BiocommonsConfig>,
    /// HGVS-RS tool configuration
    pub hgvs_rs: Option<HgvsRsConfig>,
}

/// Ferro tool configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct FerroConfig {
    /// Whether tool is enabled
    pub enabled: bool,
    /// Path to ferro reference data directory
    pub reference_dir: PathBuf,
    /// Number of parallel workers (default: number of CPU cores)
    pub parallel_workers: Option<usize>,
    /// Shuffle direction (3prime or 5prime, default: 3prime)
    pub shuffle_direction: Option<String>,
    /// Error handling mode (strict, lenient, silent, default: lenient)
    pub error_mode: Option<String>,
}

/// Mutalyzer operation mode
#[derive(Debug, Clone, Deserialize, Serialize, Default, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum MutalyzerMode {
    /// HTTP API mode (calls mutalyzer.nl or local mutalyzer server)
    #[default]
    Api,
    /// Local subprocess mode (uses Python mutalyzer package directly)
    Local,
}

/// Mutalyzer tool configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct MutalyzerConfig {
    /// Whether tool is enabled
    pub enabled: bool,
    /// Operation mode: "api" for HTTP API, "local" for Python subprocess
    #[serde(default)]
    pub mode: MutalyzerMode,
    /// Mutalyzer API URL (used in API mode)
    pub api_url: String,
    /// Request timeout in seconds (default: 30)
    pub timeout_seconds: u32,
    /// Rate limiting delay in milliseconds (default: 50, used in API mode)
    pub rate_limit_ms: Option<u64>,
    /// Health check interval in seconds (default: 60)
    pub health_check_interval: Option<u64>,
    /// Connection pool configuration (used in API mode)
    pub connection_pool: Option<ConnectionPoolConfig>,
    /// Circuit breaker configuration (used in API mode)
    pub circuit_breaker: Option<CircuitBreakerConfig>,
    /// Path to mutalyzer settings file (used in local mode)
    pub settings_file: Option<PathBuf>,
    /// Allow network access in local mode (default: false for offline/cache-only)
    #[serde(default)]
    pub allow_network: bool,
}

/// HTTP connection pool configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ConnectionPoolConfig {
    /// Maximum number of connections in pool (default: 10)
    pub max_connections: Option<usize>,
    /// Connection idle timeout in seconds (default: 30)
    pub idle_timeout_seconds: Option<u64>,
    /// Keep-alive timeout in seconds (default: 90)
    pub keep_alive_seconds: Option<u64>,
    /// Enable HTTP/2 (default: true)
    pub enable_http2: Option<bool>,
}

/// Circuit breaker configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CircuitBreakerConfig {
    /// Failure threshold before opening circuit (default: 5)
    pub failure_threshold: Option<u32>,
    /// Recovery timeout in seconds (default: 60)
    pub recovery_timeout_seconds: Option<u64>,
    /// Success threshold for closing circuit (default: 3)
    pub success_threshold: Option<u32>,
}

/// Default UTA schema name
fn default_uta_schema() -> String {
    "uta_20210129b".to_string()
}

/// Biocommons tool configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct BiocommonsConfig {
    /// Whether tool is enabled
    pub enabled: bool,
    /// UTA database URL
    pub uta_url: String,
    /// UTA database schema (default: "uta_20210129b")
    #[serde(default = "default_uta_schema")]
    pub uta_schema: String,
    /// SeqRepo directory path
    pub seqrepo_path: PathBuf,
    /// Docker container name for UTA (default: "ferro-uta")
    pub docker_container: Option<String>,
    /// Number of parallel workers (default: 1)
    pub parallel_workers: Option<usize>,
    /// Environment variables to pass to subprocess
    pub env_vars: Option<HashMap<String, String>>,
}

/// HGVS-RS tool configuration
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct HgvsRsConfig {
    /// Whether tool is enabled
    pub enabled: bool,
    /// UTA database URL
    pub uta_url: String,
    /// UTA database schema (default: "uta_20210129b")
    pub uta_schema: String,
    /// SeqRepo directory path
    pub seqrepo_path: PathBuf,
    /// LRG mapping file path (optional)
    pub lrg_mapping_file: Option<PathBuf>,
    /// Number of parallel workers (default: number of CPU cores)
    pub parallel_workers: Option<usize>,
}

impl Default for ServerConfig {
    fn default() -> Self {
        Self {
            host: "0.0.0.0".to_string(),
            port: 3000,
            max_request_size: "10MB".to_string(),
            request_timeout_seconds: 60,
            enable_cors: true,
            enable_tracing: true,
            max_concurrent_batches: Some(10),
            max_batch_size: Some(1000),
        }
    }
}

impl Default for ToolConfigs {
    fn default() -> Self {
        Self {
            ferro: Some(FerroConfig::default()),
            mutalyzer: Some(MutalyzerConfig::default()),
            biocommons: None, // Disabled by default (requires setup)
            hgvs_rs: None,    // Disabled by default (requires setup)
        }
    }
}

impl Default for FerroConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            reference_dir: PathBuf::from("./reference"),
            parallel_workers: None,
            shuffle_direction: Some("3prime".to_string()),
            error_mode: Some("lenient".to_string()),
        }
    }
}

impl Default for MutalyzerConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            mode: MutalyzerMode::default(),
            api_url: "http://localhost:8082".to_string(),
            timeout_seconds: 30,
            rate_limit_ms: Some(50),
            health_check_interval: Some(60),
            connection_pool: Some(ConnectionPoolConfig::default()),
            circuit_breaker: Some(CircuitBreakerConfig::default()),
            settings_file: None,
            allow_network: false,
        }
    }
}

impl Default for ConnectionPoolConfig {
    fn default() -> Self {
        Self {
            max_connections: Some(10),
            idle_timeout_seconds: Some(30),
            keep_alive_seconds: Some(90),
            enable_http2: Some(true),
        }
    }
}

impl Default for CircuitBreakerConfig {
    fn default() -> Self {
        Self {
            failure_threshold: Some(5),
            recovery_timeout_seconds: Some(60),
            success_threshold: Some(3),
        }
    }
}

impl ServiceConfig {
    /// Load configuration from TOML file
    pub fn from_file(path: &std::path::Path) -> Result<Self, Box<dyn std::error::Error>> {
        let content = std::fs::read_to_string(path)?;
        let config: ServiceConfig = toml::from_str(&content)?;
        Ok(config)
    }

    /// Save configuration to TOML file
    pub fn to_file(&self, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
        let content = toml::to_string_pretty(self)?;
        std::fs::write(path, content)?;
        Ok(())
    }

    /// Get list of enabled tools
    pub fn enabled_tools(&self) -> Vec<String> {
        let mut tools = Vec::new();

        if let Some(ferro) = &self.tools.ferro {
            if ferro.enabled {
                tools.push("ferro".to_string());
            }
        }

        if let Some(mutalyzer) = &self.tools.mutalyzer {
            if mutalyzer.enabled {
                tools.push("mutalyzer".to_string());
            }
        }

        if let Some(biocommons) = &self.tools.biocommons {
            if biocommons.enabled {
                tools.push("biocommons".to_string());
            }
        }

        if let Some(hgvs_rs) = &self.tools.hgvs_rs {
            if hgvs_rs.enabled {
                tools.push("hgvs-rs".to_string());
            }
        }

        tools
    }

    /// Check if a tool is enabled
    pub fn is_tool_enabled(&self, tool: &str) -> bool {
        match tool {
            "ferro" => self.tools.ferro.as_ref().is_some_and(|c| c.enabled),
            "mutalyzer" => self.tools.mutalyzer.as_ref().is_some_and(|c| c.enabled),
            "biocommons" => self.tools.biocommons.as_ref().is_some_and(|c| c.enabled),
            "hgvs-rs" => self.tools.hgvs_rs.as_ref().is_some_and(|c| c.enabled),
            _ => false,
        }
    }

    /// Validate configuration
    pub fn validate(&self) -> Result<(), String> {
        // Validate server config
        if self.server.port == 0 {
            return Err("Server port must be greater than 0".to_string());
        }

        // Validate at least one tool is enabled
        if self.enabled_tools().is_empty() {
            return Err("At least one tool must be enabled".to_string());
        }

        // Validate ferro config
        if let Some(ferro) = &self.tools.ferro {
            if ferro.enabled && !ferro.reference_dir.exists() {
                return Err(format!(
                    "Ferro reference directory does not exist: {}",
                    ferro.reference_dir.display()
                ));
            }
        }

        // Validate biocommons config
        if let Some(biocommons) = &self.tools.biocommons {
            if biocommons.enabled && !biocommons.seqrepo_path.exists() {
                return Err(format!(
                    "Biocommons seqrepo directory does not exist: {}",
                    biocommons.seqrepo_path.display()
                ));
            }
        }

        // Validate hgvs-rs config
        if let Some(hgvs_rs) = &self.tools.hgvs_rs {
            if hgvs_rs.enabled && !hgvs_rs.seqrepo_path.exists() {
                return Err(format!(
                    "HGVS-RS seqrepo directory does not exist: {}",
                    hgvs_rs.seqrepo_path.display()
                ));
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = ServiceConfig::default();
        assert_eq!(config.server.host, "0.0.0.0");
        assert_eq!(config.server.port, 3000);
        assert!(!config.enabled_tools().is_empty());
    }

    #[test]
    fn test_enabled_tools() {
        let mut config = ServiceConfig::default();

        // Disable all tools
        config.tools.ferro = Some(FerroConfig {
            enabled: false,
            ..FerroConfig::default()
        });
        config.tools.mutalyzer = Some(MutalyzerConfig {
            enabled: false,
            ..MutalyzerConfig::default()
        });
        config.tools.biocommons = None;
        config.tools.hgvs_rs = None;

        assert!(config.enabled_tools().is_empty());

        // Enable ferro
        config.tools.ferro = Some(FerroConfig {
            enabled: true,
            ..FerroConfig::default()
        });
        assert_eq!(config.enabled_tools(), vec!["ferro"]);
    }
}
