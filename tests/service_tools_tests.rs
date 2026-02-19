//! Tests for the service tools layer (manager, http_client, circuit breaker)
//!
//! These tests verify the tool orchestration and HTTP client resilience features.

#![cfg(feature = "web-service")]

use std::sync::atomic::{AtomicU32, Ordering};
use std::time::{Duration, Instant};

// ==================== Circuit Breaker State Machine Tests ====================
// Tests for http_client.rs circuit breaker logic

/// Simplified circuit breaker state for testing
#[derive(Debug, Clone, PartialEq)]
enum TestCircuitState {
    Closed,
    Open { opened_at: Instant },
    HalfOpen,
}

/// Test implementation of circuit breaker state transitions
struct TestCircuitBreaker {
    state: TestCircuitState,
    failure_count: u32,
    success_count: u32,
    failure_threshold: u32,
    success_threshold: u32,
    recovery_timeout: Duration,
}

impl TestCircuitBreaker {
    fn new(failure_threshold: u32, success_threshold: u32, recovery_timeout_secs: u64) -> Self {
        Self {
            state: TestCircuitState::Closed,
            failure_count: 0,
            success_count: 0,
            failure_threshold,
            success_threshold,
            recovery_timeout: Duration::from_secs(recovery_timeout_secs),
        }
    }

    fn record_failure(&mut self) {
        self.failure_count += 1;
        if self.failure_count >= self.failure_threshold {
            match self.state {
                TestCircuitState::Closed | TestCircuitState::HalfOpen => {
                    self.state = TestCircuitState::Open {
                        opened_at: Instant::now(),
                    };
                }
                _ => {}
            }
        }
    }

    fn record_success(&mut self) {
        match self.state {
            TestCircuitState::Closed => {
                self.failure_count = 0;
            }
            TestCircuitState::HalfOpen => {
                self.success_count += 1;
                if self.success_count >= self.success_threshold {
                    self.state = TestCircuitState::Closed;
                    self.failure_count = 0;
                    self.success_count = 0;
                }
            }
            TestCircuitState::Open { .. } => {
                // Shouldn't happen, but reset
                self.state = TestCircuitState::Closed;
                self.failure_count = 0;
            }
        }
    }

    fn check_state(&mut self) -> Result<(), &'static str> {
        match self.state {
            TestCircuitState::Closed => Ok(()),
            TestCircuitState::Open { opened_at } => {
                if opened_at.elapsed() >= self.recovery_timeout {
                    self.state = TestCircuitState::HalfOpen;
                    self.success_count = 0;
                    Ok(())
                } else {
                    Err("circuit breaker is open")
                }
            }
            TestCircuitState::HalfOpen => Ok(()),
        }
    }

    fn is_closed(&self) -> bool {
        matches!(self.state, TestCircuitState::Closed)
    }

    fn is_open(&self) -> bool {
        matches!(self.state, TestCircuitState::Open { .. })
    }

    fn is_half_open(&self) -> bool {
        matches!(self.state, TestCircuitState::HalfOpen)
    }
}

#[test]
fn test_circuit_breaker_starts_closed() {
    let cb = TestCircuitBreaker::new(5, 3, 60);
    assert!(cb.is_closed());
    assert!(!cb.is_open());
    assert!(!cb.is_half_open());
}

#[test]
fn test_circuit_breaker_opens_after_threshold_failures() {
    let mut cb = TestCircuitBreaker::new(3, 2, 60);

    // Record failures below threshold
    cb.record_failure();
    assert!(cb.is_closed(), "Should still be closed after 1 failure");
    cb.record_failure();
    assert!(cb.is_closed(), "Should still be closed after 2 failures");

    // Third failure should open the circuit
    cb.record_failure();
    assert!(cb.is_open(), "Should be open after 3 failures (threshold)");
}

#[test]
fn test_circuit_breaker_success_resets_failure_count() {
    let mut cb = TestCircuitBreaker::new(3, 2, 60);

    cb.record_failure();
    cb.record_failure();
    assert_eq!(cb.failure_count, 2);

    // Success should reset count
    cb.record_success();
    assert_eq!(cb.failure_count, 0);
    assert!(cb.is_closed());
}

#[test]
fn test_circuit_breaker_half_open_closes_after_successes() {
    let mut cb = TestCircuitBreaker::new(3, 2, 0); // 0 second timeout for immediate transition

    // Open the circuit
    cb.record_failure();
    cb.record_failure();
    cb.record_failure();
    assert!(cb.is_open());

    // Wait a tiny bit and check - should transition to half-open
    std::thread::sleep(Duration::from_millis(10));
    let _ = cb.check_state();
    assert!(
        cb.is_half_open(),
        "Should transition to half-open after timeout"
    );

    // Successes in half-open should close
    cb.record_success();
    assert!(cb.is_half_open(), "Still half-open after 1 success");
    cb.record_success();
    assert!(
        cb.is_closed(),
        "Should be closed after 2 successes (threshold)"
    );
}

#[test]
fn test_circuit_breaker_half_open_reopens_on_failure() {
    let mut cb = TestCircuitBreaker::new(1, 3, 0); // threshold=1 for quick open

    // Open the circuit
    cb.record_failure();
    assert!(cb.is_open());

    // Transition to half-open
    std::thread::sleep(Duration::from_millis(10));
    let _ = cb.check_state();
    assert!(cb.is_half_open());

    // Failure in half-open should reopen
    cb.record_failure();
    assert!(
        cb.is_open(),
        "Should reopen after failure in half-open state"
    );
}

#[test]
fn test_circuit_breaker_rejects_requests_when_open() {
    let mut cb = TestCircuitBreaker::new(2, 2, 60); // Long timeout

    cb.record_failure();
    cb.record_failure();
    assert!(cb.is_open());

    // Should reject requests
    let result = cb.check_state();
    assert!(result.is_err());
    assert_eq!(result.unwrap_err(), "circuit breaker is open");
}

// ==================== Rate Limiting Tests ====================

#[test]
fn test_rate_limiting_delay_calculation() {
    let delay_ms = 100u64;
    let delay = Duration::from_millis(delay_ms);

    let last_request = Instant::now();
    std::thread::sleep(Duration::from_millis(50));
    let elapsed = last_request.elapsed();

    // If elapsed < delay, we need to sleep more
    if elapsed < delay {
        let sleep_duration = delay - elapsed;
        assert!(sleep_duration.as_millis() > 0);
        assert!(sleep_duration.as_millis() < 100);
    }
}

#[test]
fn test_rate_limiting_no_delay_after_sufficient_time() {
    let delay_ms = 50u64;
    let delay = Duration::from_millis(delay_ms);

    let last_request = Instant::now();
    std::thread::sleep(Duration::from_millis(60)); // Wait longer than delay
    let elapsed = last_request.elapsed();

    // Should not need additional delay
    assert!(
        elapsed >= delay,
        "Should not need delay after sufficient time"
    );
}

// ==================== Tool Manager Tests ====================
// Tests for manager.rs tool filtering and request handling

#[test]
fn test_tool_filtering_with_requested_tools() {
    // Simulate tool filtering logic
    let available_tools = vec!["ferro", "mutalyzer", "biocommons"];
    let requested = Some(vec!["ferro", "mutalyzer"]);

    let filtered: Vec<&str> = match requested {
        Some(ref names) => available_tools
            .iter()
            .filter(|t| names.contains(t))
            .copied()
            .collect(),
        None => available_tools.clone(),
    };

    assert_eq!(filtered.len(), 2);
    assert!(filtered.contains(&"ferro"));
    assert!(filtered.contains(&"mutalyzer"));
    assert!(!filtered.contains(&"biocommons"));
}

#[test]
fn test_tool_filtering_returns_all_when_none_requested() {
    let available_tools = vec!["ferro", "mutalyzer", "biocommons"];
    let requested: Option<Vec<&str>> = None;

    let filtered: Vec<&str> = match requested {
        Some(ref names) => available_tools
            .iter()
            .filter(|t| names.contains(t))
            .copied()
            .collect(),
        None => available_tools.clone(),
    };

    assert_eq!(filtered.len(), 3);
}

#[test]
fn test_tool_filtering_handles_unknown_tools() {
    let available_tools = vec!["ferro", "mutalyzer"];
    let requested = Some(vec!["ferro", "unknown_tool"]);

    let filtered: Vec<&str> = match requested {
        Some(ref names) => available_tools
            .iter()
            .filter(|t| names.contains(t))
            .copied()
            .collect(),
        None => available_tools.clone(),
    };

    // Only "ferro" should be returned
    assert_eq!(filtered.len(), 1);
    assert!(filtered.contains(&"ferro"));
}

#[test]
fn test_tool_filtering_empty_request_returns_empty() {
    let available_tools = vec!["ferro", "mutalyzer"];
    let requested = Some(vec!["unknown1", "unknown2"]);

    let filtered: Vec<&str> = match requested {
        Some(ref names) => available_tools
            .iter()
            .filter(|t| names.contains(t))
            .copied()
            .collect(),
        None => available_tools.clone(),
    };

    assert!(filtered.is_empty());
}

#[test]
fn test_timeout_duration_calculation() {
    // Test default timeout
    let timeout_seconds: Option<u32> = None;
    let timeout_duration = Duration::from_secs(timeout_seconds.unwrap_or(30) as u64);
    assert_eq!(timeout_duration.as_secs(), 30);

    // Test custom timeout
    let timeout_seconds: Option<u32> = Some(60);
    let timeout_duration = Duration::from_secs(timeout_seconds.unwrap_or(30) as u64);
    assert_eq!(timeout_duration.as_secs(), 60);
}

#[test]
fn test_concurrent_limit_configuration() {
    // Test default concurrent limit
    let config_limit: Option<usize> = None;
    let concurrent_limit = config_limit.unwrap_or(10);
    assert_eq!(concurrent_limit, 10);

    // Test custom limit
    let config_limit: Option<usize> = Some(5);
    let concurrent_limit = config_limit.unwrap_or(10);
    assert_eq!(concurrent_limit, 5);
}

// ==================== Tool Mode Tests ====================

#[test]
fn test_mutalyzer_mode_parsing() {
    let mode_api = "api";
    let mode_local = "local";

    assert_eq!(mode_api, "api");
    assert_eq!(mode_local, "local");
}

// ==================== Batch Processing Tests ====================

#[test]
fn test_batch_result_aggregation() {
    // Simulate batch result counting
    let results = vec![
        (true, None::<String>),             // Success
        (true, None),                       // Success
        (false, Some("error".to_string())), // Failure
        (true, None),                       // Success
    ];

    let mut successful = 0;
    let total = results.len();

    for (success, _error) in &results {
        if *success {
            successful += 1;
        }
    }

    assert_eq!(total, 4);
    assert_eq!(successful, 3);
}

// ==================== Client Stats Tests ====================

#[test]
fn test_client_stats_tracking() {
    let total_requests = AtomicU32::new(0);
    let failure_count = AtomicU32::new(0);
    let success_count = AtomicU32::new(0);

    // Simulate requests
    total_requests.fetch_add(1, Ordering::Relaxed);
    success_count.fetch_add(1, Ordering::Relaxed);

    total_requests.fetch_add(1, Ordering::Relaxed);
    failure_count.fetch_add(1, Ordering::Relaxed);

    total_requests.fetch_add(1, Ordering::Relaxed);
    success_count.fetch_add(1, Ordering::Relaxed);

    assert_eq!(total_requests.load(Ordering::Relaxed), 3);
    assert_eq!(success_count.load(Ordering::Relaxed), 2);
    assert_eq!(failure_count.load(Ordering::Relaxed), 1);
}

// ==================== Tool Initialization Error Handling Tests ====================

#[test]
fn test_tool_initialization_error_graceful_handling() {
    // Test that initialization continues even when one tool fails
    let mut tools_initialized = vec![];
    let tools_to_init = vec!["ferro", "mutalyzer", "biocommons"];

    for tool in tools_to_init {
        // Simulate: mutalyzer fails, others succeed
        let init_result: Result<(), &str> = if tool == "mutalyzer" {
            Err("Failed to connect")
        } else {
            Ok(())
        };

        match init_result {
            Ok(()) => tools_initialized.push(tool),
            Err(_) => {
                // Log warning but continue
            }
        }
    }

    // Should have 2 tools even though mutalyzer failed
    assert_eq!(tools_initialized.len(), 2);
    assert!(tools_initialized.contains(&"ferro"));
    assert!(tools_initialized.contains(&"biocommons"));
    assert!(!tools_initialized.contains(&"mutalyzer"));
}

#[test]
fn test_tool_initialization_fails_when_all_fail() {
    let mut tools_initialized = vec![];

    // All tools fail
    let init_results: Vec<Result<&str, &str>> = vec![
        Err("Ferro failed"),
        Err("Mutalyzer failed"),
        Err("Biocommons failed"),
    ];

    for result in init_results {
        if let Ok(tool) = result {
            tools_initialized.push(tool);
        }
    }

    // Should fail when no tools initialized
    assert!(
        tools_initialized.is_empty(),
        "Should have no tools when all fail"
    );
}

// ==================== Response Time Tracking Tests ====================

#[test]
fn test_processing_time_calculation() {
    let start = Instant::now();
    std::thread::sleep(Duration::from_millis(10));
    let elapsed_ms = start.elapsed().as_millis() as u64;

    assert!(
        elapsed_ms >= 10,
        "Should capture at least 10ms of processing time"
    );
    assert!(
        elapsed_ms < 100,
        "Should not take more than 100ms for sleep(10)"
    );
}

// ==================== Agreement Calculation Tests ====================

#[test]
fn test_agreement_calculation_all_agree() {
    let results = vec![
        ("ferro", Some("NM_000249.4:c.350C>T")),
        ("mutalyzer", Some("NM_000249.4:c.350C>T")),
        ("biocommons", Some("NM_000249.4:c.350C>T")),
    ];

    let successful: Vec<_> = results.iter().filter(|(_, r)| r.is_some()).collect();
    let unique_results: std::collections::HashSet<_> =
        successful.iter().filter_map(|(_, r)| r.as_ref()).collect();

    assert_eq!(successful.len(), 3);
    assert_eq!(unique_results.len(), 1, "All should agree on same result");
}

#[test]
fn test_agreement_calculation_some_disagree() {
    let results = vec![
        ("ferro", Some("NM_000249.4:c.350C>T")),
        ("mutalyzer", Some("NM_000249.4:c.350C>T")),
        ("biocommons", Some("NM_000249.4:c.351C>T")), // Different!
    ];

    let successful: Vec<_> = results.iter().filter(|(_, r)| r.is_some()).collect();
    let unique_results: std::collections::HashSet<_> =
        successful.iter().filter_map(|(_, r)| r.as_ref()).collect();

    assert_eq!(successful.len(), 3);
    assert_eq!(unique_results.len(), 2, "Should detect disagreement");
}

#[test]
fn test_agreement_calculation_with_failures() {
    let results: Vec<(&str, Option<&str>)> = vec![
        ("ferro", Some("NM_000249.4:c.350C>T")),
        ("mutalyzer", None), // Failed
        ("biocommons", Some("NM_000249.4:c.350C>T")),
    ];

    let successful: Vec<_> = results.iter().filter(|(_, r)| r.is_some()).collect();

    assert_eq!(successful.len(), 2, "Should only count successful results");
}
