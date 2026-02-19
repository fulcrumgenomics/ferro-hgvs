//! Tool integration layer for multi-tool HGVS processing
//!
//! This module provides a unified interface for running HGVS parse and normalize
//! operations across multiple tools:
//! - ferro (native Rust)
//! - mutalyzer (HTTP API)
//! - biocommons/hgvs (Python subprocess)
//! - hgvs-rs (native Rust)

use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};
use tokio::time::timeout;

use crate::benchmark::types::ParseResult;
use crate::service::types::{
    categorize_error, health_check::HealthCheckResult, AgreementSummary, ServiceError, ToolName,
    ToolResult,
};

pub mod biocommons;
pub mod ferro;
pub mod hgvs_rs;
pub mod http_client;
pub mod manager;
pub mod mutalyzer;

// Re-export the main ToolManager
pub use manager::ToolManager;

/// Trait for HGVS tool services with standardized health checking
#[async_trait::async_trait]
pub trait HgvsToolService: Send + Sync {
    /// Parse an HGVS variant string
    async fn parse(&self, hgvs: &str) -> Result<ParseResult, ServiceError>;

    /// Normalize an HGVS variant string
    async fn normalize(&self, hgvs: &str) -> Result<ParseResult, ServiceError>;

    /// Check if the tool is available and healthy using standardized system
    async fn health_check(&self) -> HealthCheckResult;

    /// Get the tool name as enum
    fn tool_name(&self) -> ToolName;

    /// Get the tool name as string (for backward compatibility)
    fn tool_name_str(&self) -> &'static str {
        self.tool_name().as_str()
    }
}

/// Result from running multiple tools on a single variant
pub struct MultiToolResult {
    pub input: String,
    pub results: Vec<ToolResult>,
    pub agreement: AgreementSummary,
    pub total_time: Duration,
}

/// Convert ParseResult to ToolResult using secure ToolName enum
pub fn parse_result_to_tool_result(
    tool: ToolName,
    result: Result<ParseResult, ServiceError>,
    elapsed: Duration,
) -> ToolResult {
    match result {
        Ok(parse_result) => ToolResult {
            tool,
            success: parse_result.success,
            output: parse_result.output,
            error: parse_result.error.clone(),
            error_category: parse_result
                .error
                .as_ref()
                .map(|e| categorize_error(tool, e)),
            elapsed_ms: elapsed.as_millis() as u64,
            details: parse_result.details,
        },
        Err(e) => ToolResult {
            tool,
            success: false,
            output: None,
            error: Some(e.to_string()),
            error_category: Some(categorize_error(tool, &e.to_string())),
            elapsed_ms: elapsed.as_millis() as u64,
            details: None,
        },
    }
}

/// Analyze agreement across tool results
pub fn analyze_agreement(results: &[ToolResult]) -> AgreementSummary {
    let successful_results: Vec<_> = results.iter().filter(|r| r.success).collect();
    let failed_results: Vec<_> = results.iter().filter(|r| !r.success).collect();

    // Group successful results by output using secure ToolName enum
    let mut outputs: HashMap<String, Vec<ToolName>> = HashMap::new();
    for result in &successful_results {
        if let Some(output) = &result.output {
            outputs.entry(output.clone()).or_default().push(result.tool);
        }
    }

    // Check if all successful tools agree
    let all_agree = outputs.len() <= 1;

    AgreementSummary {
        all_agree,
        successful_tools: successful_results.len(),
        failed_tools: failed_results.len(),
        outputs,
    }
}

/// Run multiple tools on a single variant with timeout (optimized for reduced cloning)
pub async fn run_tools_on_variant<F, Fut>(
    tools: &[(String, Arc<dyn HgvsToolService>)],
    hgvs: &Arc<str>,
    timeout_duration: Duration,
    operation: F,
) -> MultiToolResult
where
    F: Fn(Arc<dyn HgvsToolService>, Arc<str>) -> Fut + Clone + Send + Sync + 'static,
    Fut: std::future::Future<Output = Result<ParseResult, ServiceError>> + Send + 'static,
{
    let start_time = Instant::now();
    let mut results = Vec::new();

    // Run all tools in parallel
    let mut tasks = Vec::new();

    for (tool_name, tool_service) in tools {
        let tool_name = tool_name.clone();
        let tool_service = tool_service.clone();
        let hgvs_arc = Arc::clone(hgvs); // Use Arc reference instead of string cloning
        let operation = operation.clone();

        let task = tokio::spawn(async move {
            let start = Instant::now();
            let result = timeout(timeout_duration, operation(tool_service, hgvs_arc)).await;
            let elapsed = start.elapsed();

            let final_result = match result {
                Ok(tool_result) => tool_result,
                Err(_) => Err(ServiceError::Timeout),
            };

            // Parse tool name to ToolName enum for secure type handling
            let tool_enum = ToolName::parse(&tool_name).unwrap_or(ToolName::Ferro);
            parse_result_to_tool_result(tool_enum, final_result, elapsed)
        });

        tasks.push(task);
    }

    // Collect all results
    for task in tasks {
        match task.await {
            Ok(tool_result) => results.push(tool_result),
            Err(e) => {
                // Task panicked or was cancelled
                results.push(ToolResult {
                    tool: ToolName::Ferro, // Default fallback for unknown tool
                    success: false,
                    output: None,
                    error: Some(format!("Task error: {}", e)),
                    error_category: Some(categorize_error(
                        ToolName::Ferro,
                        &format!("Task error: {}", e),
                    )),
                    elapsed_ms: 0,
                    details: None,
                });
            }
        }
    }

    let agreement = analyze_agreement(&results);
    let total_time = start_time.elapsed();

    MultiToolResult {
        input: hgvs.to_string(),
        results,
        agreement,
        total_time,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analyze_agreement_all_agree() {
        let results = vec![
            ToolResult {
                tool: ToolName::Ferro,
                success: true,
                output: Some("NM_000088.3:c.589G>T".to_string()),
                error: None,
                error_category: None,
                elapsed_ms: 10,
                details: None,
            },
            ToolResult {
                tool: ToolName::Mutalyzer,
                success: true,
                output: Some("NM_000088.3:c.589G>T".to_string()),
                error: None,
                error_category: None,
                elapsed_ms: 150,
                details: None,
            },
        ];

        let agreement = analyze_agreement(&results);
        assert!(agreement.all_agree);
        assert_eq!(agreement.successful_tools, 2);
        assert_eq!(agreement.failed_tools, 0);
        assert_eq!(agreement.outputs.len(), 1);
    }

    #[test]
    fn test_analyze_agreement_disagree() {
        let results = vec![
            ToolResult {
                tool: ToolName::Ferro,
                success: true,
                output: Some("NM_000088.3:c.589G>T".to_string()),
                error: None,
                error_category: None,
                elapsed_ms: 10,
                details: None,
            },
            ToolResult {
                tool: ToolName::Mutalyzer,
                success: true,
                output: Some("NM_000088.3:c.590G>T".to_string()), // Different output
                error: None,
                error_category: None,
                elapsed_ms: 150,
                details: None,
            },
        ];

        let agreement = analyze_agreement(&results);
        assert!(!agreement.all_agree);
        assert_eq!(agreement.successful_tools, 2);
        assert_eq!(agreement.failed_tools, 0);
        assert_eq!(agreement.outputs.len(), 2);
    }

    #[test]
    fn test_analyze_agreement_with_failures() {
        let results = vec![
            ToolResult {
                tool: ToolName::Ferro,
                success: true,
                output: Some("NM_000088.3:c.589G>T".to_string()),
                error: None,
                error_category: None,
                elapsed_ms: 10,
                details: None,
            },
            ToolResult {
                tool: ToolName::Mutalyzer,
                success: false,
                output: None,
                error: Some("Parse error".to_string()),
                error_category: Some("parse_error".to_string()),
                elapsed_ms: 150,
                details: None,
            },
        ];

        let agreement = analyze_agreement(&results);
        assert!(agreement.all_agree); // Only one successful result
        assert_eq!(agreement.successful_tools, 1);
        assert_eq!(agreement.failed_tools, 1);
        assert_eq!(agreement.outputs.len(), 1);
    }
}
