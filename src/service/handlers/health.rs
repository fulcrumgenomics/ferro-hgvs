//! Health check endpoints

use axum::{extract::State, http::StatusCode, response::Json};
use std::sync::Arc;

use crate::service::{
    server::AppState,
    tools::HgvsToolService,
    types::{
        health_check::HealthCheckResult, DetailedHealthResponse, ErrorResponse, HealthResponse,
        TestCategory, TestResult, TestStatus, ToolName, ToolStatus, ToolTestResults,
    },
};

/// Test variants organized by category for comprehensive health checks
struct HealthTestSuite {
    reference_types: Vec<(&'static str, &'static str)>,
    coordinate_types: Vec<(&'static str, &'static str)>,
    variant_types: Vec<(&'static str, &'static str)>,
}

impl HealthTestSuite {
    fn new() -> Self {
        Self {
            reference_types: vec![
                ("NM_ transcript", "NM_000088.4:c.589G>T"),
                ("NR_ non-coding", "NR_024540.1:n.100T>G"),
                ("NC_ genomic", "NC_000017.11:g.43044295T>A"),
                ("NG_ gene region", "NG_007400.1:g.8638G>T"),
                // NP_ and LRG_ may not work with all tools
            ],
            coordinate_types: vec![
                ("c. coding", "NM_000088.4:c.589G>T"),
                ("g. genomic", "NC_000017.11:g.43044295T>A"),
                ("n. non-coding", "NR_024540.1:n.100T>G"),
                // Intronic position
                ("c.+intronic", "NM_000088.4:c.588+1G>A"),
            ],
            variant_types: vec![
                ("Substitution", "NM_000088.4:c.589G>T"),
                ("Deletion", "NM_000088.4:c.589del"),
                ("Insertion", "NM_000088.4:c.589_590insA"),
                ("Duplication", "NM_000088.4:c.589dup"),
                ("Delins", "NM_000088.4:c.589delinsAT"),
            ],
        }
    }

    fn all_tests(&self) -> Vec<(&'static str, &'static str, &'static str)> {
        let mut tests = Vec::new();
        for (name, variant) in &self.reference_types {
            tests.push(("Reference Types", *name, *variant));
        }
        for (name, variant) in &self.coordinate_types {
            tests.push(("Coordinate Types", *name, *variant));
        }
        for (name, variant) in &self.variant_types {
            tests.push(("Variant Types", *name, *variant));
        }
        tests
    }
}

/// Check overall service health - returns cached basic health from periodic background task
pub async fn health_check(
    State(state): State<AppState>,
) -> Result<Json<HealthResponse>, (StatusCode, Json<ErrorResponse>)> {
    // Return cached basic health from detailed results if available
    let cached = state.health_cache.detailed.read().await;
    if let Some(ref response) = *cached {
        return Ok(Json(response.basic.clone()));
    }
    drop(cached);

    // If no cache yet, return a minimal response indicating startup
    // (background task will populate the cache shortly)
    Ok(Json(HealthResponse {
        status: "starting".to_string(),
        available_tools: Vec::new(),
        unavailable_tools: Vec::new(),
        tools: Vec::new(),
    }))
}

/// Get detailed tool status - returns cached tool statuses from periodic background task
pub async fn tools_status(
    State(state): State<AppState>,
) -> Result<Json<Vec<ToolStatus>>, (StatusCode, Json<ErrorResponse>)> {
    // Return cached tool statuses from detailed results if available
    let cached = state.health_cache.detailed.read().await;
    if let Some(ref response) = *cached {
        return Ok(Json(response.basic.tools.clone()));
    }
    drop(cached);

    // If no cache yet, return empty list
    Ok(Json(Vec::new()))
}

/// Detailed health check - returns cached results from periodic background task
pub async fn detailed_health_check(
    State(state): State<AppState>,
) -> Result<Json<DetailedHealthResponse>, (StatusCode, Json<ErrorResponse>)> {
    // Return cached results if available
    let cached = state.health_cache.detailed.read().await;
    if let Some(ref response) = *cached {
        return Ok(Json(response.clone()));
    }
    drop(cached);

    // If no cache yet, return basic health with empty test results
    // (background task will populate the cache shortly)
    let basic = get_basic_health(&state).await;
    Ok(Json(DetailedHealthResponse {
        basic,
        test_results: Vec::new(),
    }))
}

/// Run detailed health check and return raw result (used by background task)
pub async fn run_detailed_health_check(
    state: &crate::service::server::AppState,
) -> Result<DetailedHealthResponse, String> {
    // Get basic health
    let basic = get_basic_health(state).await;

    // Get test suite
    let test_suite = HealthTestSuite::new();
    let all_tests = test_suite.all_tests();

    // Get all tool services
    let tools = state.tool_manager.get_all_tools();

    // Run tests for each tool
    let mut test_results = Vec::new();

    for (tool_name, tool_service) in tools {
        let tool_enum = ToolName::parse(&tool_name).unwrap_or(ToolName::Ferro);
        let mode = state.tool_manager.get_tool_mode(&tool_name);

        // Group results by category
        let mut category_results: std::collections::HashMap<&str, Vec<TestResult>> =
            std::collections::HashMap::new();
        let mut passed_count = 0;
        let mut applicable_count = 0; // Tests that are not N/A

        for (category, test_name, variant) in &all_tests {
            let mut result = run_single_test(&tool_service, variant).await;
            result.name = test_name.to_string();

            // Count passed and applicable tests (N/A doesn't count toward either)
            if result.status != TestStatus::Na {
                applicable_count += 1;
                if result.status == TestStatus::Pass {
                    passed_count += 1;
                }
            }

            category_results.entry(category).or_default().push(result);
        }

        // Convert to categories
        let categories = vec![
            TestCategory {
                name: "Reference Types".to_string(),
                tests: category_results
                    .remove("Reference Types")
                    .unwrap_or_default(),
            },
            TestCategory {
                name: "Coordinate Types".to_string(),
                tests: category_results
                    .remove("Coordinate Types")
                    .unwrap_or_default(),
            },
            TestCategory {
                name: "Variant Types".to_string(),
                tests: category_results.remove("Variant Types").unwrap_or_default(),
            },
        ];

        test_results.push(ToolTestResults {
            tool: tool_enum,
            passed: passed_count,
            total: applicable_count, // Only count applicable tests (for health rate)
            total_tests: all_tests.len(), // All tests including N/A (for coverage rate)
            mode,
            categories,
        });
    }

    Ok(DetailedHealthResponse {
        basic,
        test_results,
    })
}

/// Check if an error message indicates "not applicable" (unsupported or missing data)
/// vs an actual failure
fn is_na_error(error: &str) -> bool {
    let error_lower = error.to_lowercase();

    // First, check for patterns that indicate REAL failures (not N/A)
    // These take precedence over N/A patterns
    let failure_patterns = [
        "esequencemismatch", // mutalyzer: reference base mismatch
        "sequence mismatch",
        "mismatch",
        "invalid",
    ];

    if failure_patterns
        .iter()
        .any(|pattern| error_lower.contains(pattern))
    {
        return false;
    }

    // Patterns that indicate N/A (not supported or missing data)
    let na_patterns = [
        // Missing reference/transcript (be specific to avoid false positives)
        "reference not found",
        "transcript not found",
        "sequence not found",
        "accession not found",
        "unknown reference",
        "unknown transcript",
        "no transcript",
        "could not find reference",
        "could not find transcript",
        "does not exist",
        "no such file",
        "no such reference",
        // Unsupported features
        "not supported",
        "unsupported",
        "not implemented",
        "unimplemented",
        // Missing data/database issues
        "no data for",
        "missing data",
        "not in database",
        "not in uta",
        // Network/cache issues (data not available locally)
        "network access required",
        "network access disabled",
        "offline mode",
        // Connection/timeout for external services (treat as N/A for health check)
        "connection refused",
        "connection timed out",
        "timed out",
        // Specific tool messages
        "ereferenceerror",           // mutalyzer reference error (not mismatch)
        "eintronic", // mutalyzer: intronic position on non-genomic reference (needs --rewrite-intronic)
        "cannot normalize intronic", // hgvs-rs: intronic variants not supported
        "problem accessing data", // hgvs-rs: reference data not available in SeqRepo/UTA
        "no exon",
        "no cds",
        // Tool/package not available
        "not installed",       // biocommons: pip install hgvs not run
        "feature not enabled", // hgvs-rs: compiled without feature
    ];

    na_patterns
        .iter()
        .any(|pattern| error_lower.contains(pattern))
}

/// Run a single test variant against a tool
async fn run_single_test(tool: &Arc<dyn HgvsToolService>, variant: &str) -> TestResult {
    match tool.normalize(variant).await {
        Ok(result) => {
            if result.success {
                TestResult {
                    name: String::new(),
                    variant: variant.to_string(),
                    status: TestStatus::Pass,
                    passed: true,
                    error: None,
                }
            } else {
                // Check if the error indicates N/A or actual failure
                let error_msg = result.error.as_deref().unwrap_or("");
                let status = if is_na_error(error_msg) {
                    TestStatus::Na
                } else {
                    TestStatus::Fail
                };
                TestResult {
                    name: String::new(),
                    variant: variant.to_string(),
                    status,
                    passed: false,
                    error: result.error,
                }
            }
        }
        Err(e) => {
            let error_msg = e.to_string();
            let status = if is_na_error(&error_msg) {
                TestStatus::Na
            } else {
                TestStatus::Fail
            };
            TestResult {
                name: String::new(),
                variant: variant.to_string(),
                status,
                passed: false,
                error: Some(error_msg),
            }
        }
    }
}

/// Get basic health response (extracted for reuse)
async fn get_basic_health(state: &AppState) -> HealthResponse {
    let health_results = state.tool_manager.health_check_all().await;

    let mut available_tools = Vec::new();
    let mut unavailable_tools = Vec::new();
    let mut tool_statuses = Vec::new();

    for (tool_name_str, result) in health_results {
        let tool_name = ToolName::parse(&tool_name_str).unwrap_or(ToolName::Ferro);

        let (available, status) = match &result {
            HealthCheckResult::Healthy => {
                available_tools.push(tool_name);
                (true, "Available".to_string())
            }
            HealthCheckResult::Degraded { reason } => {
                available_tools.push(tool_name);
                (true, format!("Degraded: {}", reason))
            }
            HealthCheckResult::Unhealthy { reason } => {
                unavailable_tools.push(tool_name);
                (false, format!("Unhealthy: {}", reason))
            }
        };

        tool_statuses.push(ToolStatus {
            tool: tool_name,
            available,
            status,
            last_check: chrono::Utc::now().to_rfc3339(),
            mode: state.tool_manager.get_tool_mode(&tool_name_str),
        });
    }

    let overall_status = if unavailable_tools.is_empty() {
        "healthy"
    } else if available_tools.is_empty() {
        "unhealthy"
    } else {
        "degraded"
    };

    HealthResponse {
        status: overall_status.to_string(),
        available_tools,
        unavailable_tools,
        tools: tool_statuses,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_na_error_reference_not_found() {
        assert!(is_na_error("Reference not found for NM_000001"));
        assert!(is_na_error("Transcript not found in UTA"));
        assert!(is_na_error("Could not find reference"));
        assert!(is_na_error("Unknown reference NC_000001"));
    }

    #[test]
    fn test_is_na_error_unsupported() {
        assert!(is_na_error("Feature not supported"));
        assert!(is_na_error("This is unsupported"));
        assert!(is_na_error("Not implemented yet"));
    }

    #[test]
    fn test_is_na_error_network() {
        assert!(is_na_error("Connection refused"));
        assert!(is_na_error("Connection timed out"));
        assert!(is_na_error("Request timed out after 30 seconds"));
    }

    #[test]
    fn test_is_na_error_tool_specific() {
        assert!(is_na_error("Cannot normalize intronic variants"));
        assert!(is_na_error("Problem accessing data in SeqRepo"));
        assert!(is_na_error("hgvs not installed"));
        assert!(is_na_error("Feature not enabled"));
        assert!(is_na_error("EReferenceError: unknown reference"));
        assert!(is_na_error("EIntronic: cannot process"));
    }

    #[test]
    fn test_is_na_error_mismatch_is_failure() {
        // Mismatch errors should NOT be classified as N/A
        assert!(!is_na_error("ESequenceMismatch: expected C got G"));
        assert!(!is_na_error("Sequence mismatch at position 100"));
        assert!(!is_na_error("Reference mismatch"));
    }

    #[test]
    fn test_is_na_error_invalid_is_failure() {
        // Invalid errors should NOT be classified as N/A
        assert!(!is_na_error("Invalid HGVS format"));
        assert!(!is_na_error("Invalid position"));
    }

    #[test]
    fn test_is_na_error_general_errors() {
        // General errors that don't match patterns should NOT be N/A
        assert!(!is_na_error("Some random error"));
        assert!(!is_na_error("Failed to process"));
    }

    #[test]
    fn test_health_test_suite_creation() {
        let suite = HealthTestSuite::new();

        // Verify we have tests in all categories
        assert!(!suite.reference_types.is_empty());
        assert!(!suite.coordinate_types.is_empty());
        assert!(!suite.variant_types.is_empty());

        // Verify all_tests combines them
        let all = suite.all_tests();
        assert!(all.len() >= 10);
    }
}
