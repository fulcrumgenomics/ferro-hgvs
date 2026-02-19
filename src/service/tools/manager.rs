//! Tool manager for coordinating multiple HGVS tools

use std::collections::HashMap;
use std::sync::Arc;
use std::time::Duration;

use crate::service::{
    config::ServiceConfig,
    tools::{
        biocommons::BiocommonsService, ferro::FerroService, hgvs_rs::HgvsRsService,
        mutalyzer::MutalyzerService, run_tools_on_variant, HgvsToolService,
    },
    types::{BatchResponse, ServiceError, SingleResponse, VariantBatchResult},
};

/// Manages multiple HGVS tool services
pub struct ToolManager {
    /// Available tool services by name
    tools: HashMap<String, Arc<dyn HgvsToolService>>,
    /// Configuration
    config: Arc<ServiceConfig>,
}

impl ToolManager {
    /// Create a new ToolManager and initialize all enabled tools
    pub fn new(config: &ServiceConfig) -> Result<Self, ServiceError> {
        let mut tools: HashMap<String, Arc<dyn HgvsToolService>> = HashMap::new();

        // Initialize ferro if enabled
        if let Some(ferro_config) = &config.tools.ferro {
            if ferro_config.enabled {
                tracing::info!("Initializing ferro tool");
                match FerroService::new(ferro_config) {
                    Ok(service) => {
                        tools.insert("ferro".to_string(), Arc::new(service));
                        tracing::info!("Ferro tool initialized successfully");
                    }
                    Err(e) => {
                        tracing::warn!("Failed to initialize ferro tool: {}", e);
                        if tools.is_empty() {
                            return Err(ServiceError::ConfigError(
                                "No tools could be initialized".to_string(),
                            ));
                        }
                    }
                }
            }
        }

        // Initialize mutalyzer if enabled
        if let Some(mutalyzer_config) = &config.tools.mutalyzer {
            if mutalyzer_config.enabled {
                tracing::info!("Initializing mutalyzer tool");
                match MutalyzerService::new(mutalyzer_config) {
                    Ok(service) => {
                        tools.insert("mutalyzer".to_string(), Arc::new(service));
                        tracing::info!("Mutalyzer tool initialized successfully");
                    }
                    Err(e) => {
                        tracing::warn!("Failed to initialize mutalyzer tool: {}", e);
                        // Don't fail completely - just continue without this tool
                    }
                }
            }
        }

        // Initialize biocommons if enabled
        if let Some(biocommons_config) = &config.tools.biocommons {
            if biocommons_config.enabled {
                tracing::info!("Initializing biocommons tool");
                match BiocommonsService::new(biocommons_config) {
                    Ok(service) => {
                        tools.insert("biocommons".to_string(), Arc::new(service));
                        tracing::info!("Biocommons tool initialized successfully");
                    }
                    Err(e) => {
                        tracing::warn!("Failed to initialize biocommons tool: {}", e);
                        // Don't fail completely - just continue without this tool
                    }
                }
            }
        }

        // Initialize hgvs-rs if enabled
        // Note: hgvs-rs uses the postgres crate which calls block_on internally.
        // This conflicts with tokio runtime, so we must create the service in a
        // separate thread to avoid "Cannot start a runtime from within a runtime" panic.
        if let Some(hgvs_rs_config) = &config.tools.hgvs_rs {
            if hgvs_rs_config.enabled {
                tracing::info!("Initializing hgvs-rs tool");
                let config_clone = hgvs_rs_config.clone();

                // Spawn a dedicated thread for hgvs-rs initialization to avoid runtime conflict
                let handle = std::thread::spawn(move || HgvsRsService::new(&config_clone));

                match handle.join() {
                    Ok(Ok(service)) => {
                        tools.insert("hgvs-rs".to_string(), Arc::new(service));
                        tracing::info!("HGVS-RS tool initialized successfully");
                    }
                    Ok(Err(e)) => {
                        tracing::warn!("Failed to initialize hgvs-rs tool: {}", e);
                        // Don't fail completely - just continue without this tool
                    }
                    Err(_) => {
                        tracing::warn!("hgvs-rs initialization thread panicked");
                        // Don't fail completely - just continue without this tool
                    }
                }
            }
        }

        if tools.is_empty() {
            return Err(ServiceError::ConfigError(
                "No tools are enabled or could be initialized".to_string(),
            ));
        }

        tracing::info!(
            "Initialized {} tools: {:?}",
            tools.len(),
            tools.keys().collect::<Vec<_>>()
        );

        Ok(Self {
            tools,
            config: Arc::new(config.clone()),
        })
    }

    /// Create an empty ToolManager for testing purposes
    ///
    /// This creates a manager with no tools initialized, useful for testing
    /// handlers that don't require actual tool calls.
    #[cfg(any(test, feature = "dev"))]
    pub fn empty() -> Self {
        Self {
            tools: HashMap::new(),
            config: Arc::new(ServiceConfig::default()),
        }
    }

    /// Get list of available tool names
    pub fn available_tools(&self) -> Vec<String> {
        self.tools.keys().cloned().collect()
    }

    /// Get all tool services (for comprehensive testing)
    pub fn get_all_tools(&self) -> Vec<(String, Arc<dyn HgvsToolService>)> {
        self.tools
            .iter()
            .map(|(name, service)| (name.clone(), service.clone()))
            .collect()
    }

    /// Check if a tool is available
    pub fn is_tool_available(&self, tool: &str) -> bool {
        self.tools.contains_key(tool)
    }

    /// Get filtered list of tools based on request
    pub fn get_requested_tools(
        &self,
        requested: Option<&[String]>,
    ) -> Vec<(String, Arc<dyn HgvsToolService>)> {
        match requested {
            Some(tool_names) => tool_names
                .iter()
                .filter_map(|name| {
                    self.tools
                        .get(name)
                        .map(|service| (name.clone(), service.clone()))
                })
                .collect(),
            None => {
                // Return all available tools
                self.tools
                    .iter()
                    .map(|(name, service)| (name.clone(), service.clone()))
                    .collect()
            }
        }
    }

    /// Parse a single HGVS variant using specified tools
    pub async fn parse_single(
        &self,
        hgvs: &str,
        tools: Option<&[String]>,
        timeout_seconds: Option<u32>,
    ) -> Result<SingleResponse, ServiceError> {
        let timeout_duration = Duration::from_secs(timeout_seconds.unwrap_or(30) as u64);
        let requested_tools = self.get_requested_tools(tools);

        if requested_tools.is_empty() {
            return Err(ServiceError::BadRequest(
                "No valid tools specified or available".to_string(),
            ));
        }

        let hgvs_arc: Arc<str> = Arc::from(hgvs);
        let result = run_tools_on_variant(
            &requested_tools,
            &hgvs_arc,
            timeout_duration,
            |tool, hgvs_str| async move { tool.parse(&hgvs_str).await },
        )
        .await;

        Ok(SingleResponse {
            input: result.input,
            results: result.results,
            agreement: result.agreement,
            processing_time_ms: result.total_time.as_millis() as u64,
        })
    }

    /// Normalize a single HGVS variant using specified tools
    pub async fn normalize_single(
        &self,
        hgvs: &str,
        tools: Option<&[String]>,
        timeout_seconds: Option<u32>,
    ) -> Result<SingleResponse, ServiceError> {
        let timeout_duration = Duration::from_secs(timeout_seconds.unwrap_or(30) as u64);
        let requested_tools = self.get_requested_tools(tools);

        if requested_tools.is_empty() {
            return Err(ServiceError::BadRequest(
                "No valid tools specified or available".to_string(),
            ));
        }

        let hgvs_arc: Arc<str> = Arc::from(hgvs);
        let result = run_tools_on_variant(
            &requested_tools,
            &hgvs_arc,
            timeout_duration,
            |tool, hgvs_str| async move { tool.normalize(&hgvs_str).await },
        )
        .await;

        Ok(SingleResponse {
            input: result.input,
            results: result.results,
            agreement: result.agreement,
            processing_time_ms: result.total_time.as_millis() as u64,
        })
    }

    /// Parse multiple HGVS variants using specified tools
    pub async fn parse_batch(
        &self,
        variants: &[String],
        tools: Option<&[String]>,
        timeout_seconds: Option<u32>,
    ) -> Result<BatchResponse, ServiceError> {
        let start_time = std::time::Instant::now();
        let timeout_duration = Duration::from_secs(timeout_seconds.unwrap_or(30) as u64);
        let requested_tools = self.get_requested_tools(tools);

        if requested_tools.is_empty() {
            return Err(ServiceError::BadRequest(
                "No valid tools specified or available".to_string(),
            ));
        }

        let mut results = Vec::new();
        let mut successful_variants = 0;

        // Process variants in parallel with configurable concurrency limits
        let concurrent_limit = self.config.server.max_concurrent_batches.unwrap_or(10);
        let semaphore = Arc::new(tokio::sync::Semaphore::new(concurrent_limit));

        // Convert to Arc to avoid cloning expensive data structures
        let tools_arc = Arc::new(requested_tools);

        let mut tasks = Vec::new();
        for hgvs in variants {
            // Use Arc<str> to avoid string cloning
            let hgvs_arc: Arc<str> = Arc::from(hgvs.as_str());
            let tools_ref = Arc::clone(&tools_arc);
            let semaphore = semaphore.clone();

            let task = tokio::spawn(async move {
                let _permit = semaphore.acquire().await.unwrap();

                run_tools_on_variant(
                    &tools_ref,
                    &hgvs_arc,
                    timeout_duration,
                    |tool, hgvs_str| async move { tool.parse(&hgvs_str).await },
                )
                .await
            });

            tasks.push(task);
        }

        // Collect all results
        for task in tasks {
            match task.await {
                Ok(result) => {
                    if result.agreement.successful_tools > 0 {
                        successful_variants += 1;
                    }

                    results.push(VariantBatchResult {
                        input: result.input,
                        results: result.results,
                        agreement: result.agreement,
                    });
                }
                Err(e) => {
                    tracing::error!("Batch processing task failed: {}", e);
                    // Continue processing other variants
                }
            }
        }

        let total_time = start_time.elapsed();

        Ok(BatchResponse {
            total_variants: variants.len(),
            successful_variants,
            results,
            total_processing_time_ms: total_time.as_millis() as u64,
        })
    }

    /// Normalize multiple HGVS variants using specified tools
    pub async fn normalize_batch(
        &self,
        variants: &[String],
        tools: Option<&[String]>,
        timeout_seconds: Option<u32>,
    ) -> Result<BatchResponse, ServiceError> {
        let start_time = std::time::Instant::now();
        let timeout_duration = Duration::from_secs(timeout_seconds.unwrap_or(30) as u64);
        let requested_tools = self.get_requested_tools(tools);

        if requested_tools.is_empty() {
            return Err(ServiceError::BadRequest(
                "No valid tools specified or available".to_string(),
            ));
        }

        let mut results = Vec::new();
        let mut successful_variants = 0;

        // Process variants in parallel with configurable concurrency limits
        let concurrent_limit = self.config.server.max_concurrent_batches.unwrap_or(10);
        let semaphore = Arc::new(tokio::sync::Semaphore::new(concurrent_limit));

        // Convert to Arc to avoid cloning expensive data structures
        let tools_arc = Arc::new(requested_tools);

        let mut tasks = Vec::new();
        for hgvs in variants {
            // Use Arc<str> to avoid string cloning
            let hgvs_arc: Arc<str> = Arc::from(hgvs.as_str());
            let tools_ref = Arc::clone(&tools_arc);
            let semaphore = semaphore.clone();

            let task = tokio::spawn(async move {
                let _permit = semaphore.acquire().await.unwrap();

                run_tools_on_variant(
                    &tools_ref,
                    &hgvs_arc,
                    timeout_duration,
                    |tool, hgvs_str| async move { tool.normalize(&hgvs_str).await },
                )
                .await
            });

            tasks.push(task);
        }

        // Collect all results
        for task in tasks {
            match task.await {
                Ok(result) => {
                    if result.agreement.successful_tools > 0 {
                        successful_variants += 1;
                    }

                    results.push(VariantBatchResult {
                        input: result.input,
                        results: result.results,
                        agreement: result.agreement,
                    });
                }
                Err(e) => {
                    tracing::error!("Batch processing task failed: {}", e);
                    // Continue processing other variants
                }
            }
        }

        let total_time = start_time.elapsed();

        Ok(BatchResponse {
            total_variants: variants.len(),
            successful_variants,
            results,
            total_processing_time_ms: total_time.as_millis() as u64,
        })
    }

    /// Check health of all tools
    pub async fn health_check_all(
        &self,
    ) -> HashMap<String, crate::service::types::health_check::HealthCheckResult> {
        let mut results = HashMap::new();

        for (name, tool) in &self.tools {
            let result = tool.health_check().await;
            results.insert(name.clone(), result);
        }

        results
    }

    /// Get mode information for a tool (e.g., "api" or "local" for mutalyzer)
    pub fn get_tool_mode(&self, tool_name: &str) -> Option<String> {
        match tool_name {
            "mutalyzer" => self.config.tools.mutalyzer.as_ref().map(|c| match c.mode {
                crate::service::config::MutalyzerMode::Api => "api".to_string(),
                crate::service::config::MutalyzerMode::Local => "local".to_string(),
            }),
            // Other tools don't have modes currently
            _ => None,
        }
    }
}
