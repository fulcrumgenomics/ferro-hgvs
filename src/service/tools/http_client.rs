//! Advanced HTTP client with connection pooling and circuit breaker

use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use reqwest::Client;
use tokio::sync::{Mutex, Semaphore};
use tokio::time::sleep;

use crate::service::{
    config::{CircuitBreakerConfig, ConnectionPoolConfig},
    types::ServiceError,
};

/// Circuit breaker state
#[derive(Debug, Clone, PartialEq)]
pub enum CircuitState {
    Closed,
    Open { opened_at: Instant },
    HalfOpen,
}

/// Advanced HTTP client with connection pooling and circuit breaker
#[derive(Debug)]
pub struct EnhancedHttpClient {
    client: Client,
    /// Rate limiting semaphore
    rate_limiter: Arc<Semaphore>,
    /// Rate limiting delay between requests
    rate_limit_delay: Option<Duration>,
    /// Last request time for rate limiting
    last_request_time: Arc<Mutex<Option<Instant>>>,
    /// Circuit breaker state
    circuit_state: Arc<Mutex<CircuitState>>,
    /// Circuit breaker configuration
    circuit_config: CircuitBreakerConfig,
    /// Failure counter
    failure_count: Arc<AtomicU32>,
    /// Success counter (for half-open state)
    success_count: Arc<AtomicU32>,
    /// Total request counter
    request_count: Arc<AtomicU64>,
}

impl EnhancedHttpClient {
    /// Create a new enhanced HTTP client
    pub fn new(
        pool_config: Option<&ConnectionPoolConfig>,
        circuit_config: Option<&CircuitBreakerConfig>,
        rate_limit_ms: Option<u64>,
    ) -> Result<Self, ServiceError> {
        // Build reqwest client with connection pooling
        let mut client_builder = Client::builder();

        if let Some(pool) = pool_config {
            if let Some(max_conn) = pool.max_connections {
                // Note: reqwest doesn't expose max_connections directly in public API
                // This is a limitation we'll document
                client_builder = client_builder.pool_max_idle_per_host(max_conn);
            }

            if let Some(idle_timeout) = pool.idle_timeout_seconds {
                client_builder =
                    client_builder.pool_idle_timeout(Duration::from_secs(idle_timeout));
            }

            if let Some(keep_alive) = pool.keep_alive_seconds {
                client_builder = client_builder.tcp_keepalive(Duration::from_secs(keep_alive));
            }

            // Note: http2_prior_knowledge() assumes server speaks HTTP/2 without negotiation
            // This breaks HTTPS connections that use ALPN. Only use for local HTTP servers.
            // For HTTPS, let the TLS ALPN negotiation handle protocol selection.
            // if pool.enable_http2 == Some(true) && !url_is_https {
            //     client_builder = client_builder.http2_prior_knowledge();
            // }
        }

        // Enable compression
        client_builder = client_builder.gzip(true).deflate(true);

        let client = client_builder.build().map_err(|e| {
            ServiceError::InternalError(format!("Failed to create HTTP client: {}", e))
        })?;

        // Create rate limiter - allow 1 request at a time for strict rate limiting
        let rate_limiter = Arc::new(Semaphore::new(1));
        let rate_limit_delay = rate_limit_ms.map(Duration::from_millis);

        // Initialize circuit breaker
        let circuit_config = circuit_config.cloned().unwrap_or_default();

        Ok(Self {
            client,
            rate_limiter,
            rate_limit_delay,
            last_request_time: Arc::new(Mutex::new(None)),
            circuit_state: Arc::new(Mutex::new(CircuitState::Closed)),
            circuit_config,
            failure_count: Arc::new(AtomicU32::new(0)),
            success_count: Arc::new(AtomicU32::new(0)),
            request_count: Arc::new(AtomicU64::new(0)),
        })
    }

    /// Make an HTTP GET request with circuit breaker protection
    pub async fn get(&self, url: &str) -> Result<reqwest::Response, ServiceError> {
        self.execute_request(|| self.client.get(url)).await
    }

    /// Execute HTTP request with full protection (rate limiting + circuit breaker)
    async fn execute_request<F>(
        &self,
        request_builder: F,
    ) -> Result<reqwest::Response, ServiceError>
    where
        F: FnOnce() -> reqwest::RequestBuilder,
    {
        // Check circuit breaker state
        self.check_circuit_breaker().await?;

        // Acquire rate limiting permit
        let _permit = self
            .rate_limiter
            .acquire()
            .await
            .map_err(|e| ServiceError::InternalError(format!("Rate limiter error: {}", e)))?;

        // Apply rate limiting delay
        if let Some(delay) = self.rate_limit_delay {
            let mut last_time = self.last_request_time.lock().await;
            if let Some(last) = *last_time {
                let elapsed = last.elapsed();
                if elapsed < delay {
                    let sleep_duration = delay - elapsed;
                    sleep(sleep_duration).await;
                }
            }
            *last_time = Some(Instant::now());
        }

        // Execute the request
        self.request_count.fetch_add(1, Ordering::Relaxed);
        let result = request_builder()
            .send()
            .await
            .map_err(|e| ServiceError::InternalError(format!("HTTP request failed: {}", e)));

        // Update circuit breaker based on result
        match &result {
            Ok(response) => {
                if response.status().is_success() {
                    self.record_success().await;
                } else {
                    self.record_failure().await;
                }
            }
            Err(_) => {
                self.record_failure().await;
            }
        }

        result
    }

    /// Check circuit breaker state and potentially fail fast
    async fn check_circuit_breaker(&self) -> Result<(), ServiceError> {
        let mut state = self.circuit_state.lock().await;

        match *state {
            CircuitState::Closed => Ok(()),
            CircuitState::Open { opened_at } => {
                let recovery_timeout =
                    Duration::from_secs(self.circuit_config.recovery_timeout_seconds.unwrap_or(60));

                if opened_at.elapsed() >= recovery_timeout {
                    // Transition to half-open
                    *state = CircuitState::HalfOpen;
                    self.success_count.store(0, Ordering::Relaxed);
                    tracing::info!("Circuit breaker transitioning to half-open state");
                    Ok(())
                } else {
                    Err(ServiceError::CircuitBreakerOpen)
                }
            }
            CircuitState::HalfOpen => Ok(()),
        }
    }

    /// Record a successful request
    async fn record_success(&self) {
        let mut state = self.circuit_state.lock().await;

        match *state {
            CircuitState::Closed => {
                // Reset failure count on success
                self.failure_count.store(0, Ordering::Relaxed);
            }
            CircuitState::HalfOpen => {
                let success_count = self.success_count.fetch_add(1, Ordering::Relaxed) + 1;
                let success_threshold = self.circuit_config.success_threshold.unwrap_or(3);

                if success_count >= success_threshold {
                    // Transition back to closed
                    *state = CircuitState::Closed;
                    self.failure_count.store(0, Ordering::Relaxed);
                    self.success_count.store(0, Ordering::Relaxed);
                    tracing::info!(
                        "Circuit breaker closed after {} successful requests",
                        success_count
                    );
                }
            }
            CircuitState::Open { .. } => {
                // Should not happen, but reset if it does
                *state = CircuitState::Closed;
                self.failure_count.store(0, Ordering::Relaxed);
            }
        }
    }

    /// Record a failed request
    async fn record_failure(&self) {
        let failure_count = self.failure_count.fetch_add(1, Ordering::Relaxed) + 1;
        let failure_threshold = self.circuit_config.failure_threshold.unwrap_or(5);

        if failure_count >= failure_threshold {
            let mut state = self.circuit_state.lock().await;
            match *state {
                CircuitState::Closed | CircuitState::HalfOpen => {
                    *state = CircuitState::Open {
                        opened_at: Instant::now(),
                    };
                    tracing::warn!("Circuit breaker opened after {} failures", failure_count);
                }
                CircuitState::Open { .. } => {
                    // Already open, do nothing
                }
            }
        }
    }

    /// Get current circuit breaker state for monitoring
    pub async fn circuit_state(&self) -> CircuitState {
        self.circuit_state.lock().await.clone()
    }

    /// Get request statistics
    pub fn get_stats(&self) -> ClientStats {
        ClientStats {
            total_requests: self.request_count.load(Ordering::Relaxed),
            failure_count: self.failure_count.load(Ordering::Relaxed),
            success_count: self.success_count.load(Ordering::Relaxed),
        }
    }
}

/// HTTP client statistics
#[derive(Debug, Clone)]
pub struct ClientStats {
    pub total_requests: u64,
    pub failure_count: u32,
    pub success_count: u32,
}
