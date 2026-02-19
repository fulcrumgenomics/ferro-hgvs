//! Batch processor implementation.

use std::time::{Duration, Instant};

use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::Normalizer;
use crate::reference::ReferenceProvider;

/// Configuration for batch processing.
#[derive(Debug, Clone)]
pub struct BatchConfig {
    /// Whether to continue processing on errors.
    pub continue_on_error: bool,
    /// Callback frequency (call progress callback every N items).
    pub progress_interval: usize,
}

impl Default for BatchConfig {
    fn default() -> Self {
        Self {
            continue_on_error: true,
            progress_interval: 100,
        }
    }
}

impl BatchConfig {
    /// Create a new batch configuration with defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Configure whether to continue on errors.
    pub fn continue_on_error(mut self, continue_on_error: bool) -> Self {
        self.continue_on_error = continue_on_error;
        self
    }

    /// Set the progress callback interval.
    pub fn progress_interval(mut self, interval: usize) -> Self {
        self.progress_interval = interval;
        self
    }
}

/// Progress information for batch operations.
#[derive(Debug, Clone)]
pub struct BatchProgress {
    /// Total items to process.
    pub total: usize,
    /// Items processed so far.
    pub processed: usize,
    /// Successful items so far.
    pub success: usize,
    /// Failed items so far.
    pub errors: usize,
    /// Time elapsed since start.
    pub elapsed: Duration,
}

impl BatchProgress {
    /// Calculate completion percentage.
    pub fn percent(&self) -> f64 {
        if self.total == 0 {
            100.0
        } else {
            (self.processed as f64 / self.total as f64) * 100.0
        }
    }

    /// Calculate processing rate (items per second).
    ///
    /// Returns 0.0 if no time has elapsed yet. This is intentional because:
    /// 1. Division by zero would produce infinity, which isn't useful
    /// 2. A rate of 0.0 correctly indicates we can't estimate throughput yet
    /// 3. Callers should check for 0.0 before using the rate for calculations
    pub fn items_per_second(&self) -> f64 {
        let secs = self.elapsed.as_secs_f64();
        if secs < f64::EPSILON {
            0.0
        } else {
            self.processed as f64 / secs
        }
    }

    /// Estimate remaining time based on current rate.
    pub fn estimated_remaining(&self) -> Option<Duration> {
        let rate = self.items_per_second();
        if rate == 0.0 {
            return None;
        }
        let remaining_items = self.total.saturating_sub(self.processed);
        let remaining_secs = remaining_items as f64 / rate;
        Some(Duration::from_secs_f64(remaining_secs))
    }
}

/// Result of a single item in a batch operation.
#[derive(Debug, Clone)]
pub enum ItemResult<T> {
    /// Successful processing.
    Ok(T),
    /// Failed processing with error.
    Err {
        /// Original input (if available).
        input: Option<String>,
        /// Error that occurred.
        error: FerroError,
    },
}

impl<T> ItemResult<T> {
    /// Check if this is a success.
    pub fn is_ok(&self) -> bool {
        matches!(self, ItemResult::Ok(_))
    }

    /// Check if this is an error.
    pub fn is_err(&self) -> bool {
        matches!(self, ItemResult::Err { .. })
    }

    /// Get the success value if present.
    pub fn ok(self) -> Option<T> {
        match self {
            ItemResult::Ok(v) => Some(v),
            ItemResult::Err { .. } => None,
        }
    }

    /// Get the error if present.
    pub fn err(self) -> Option<FerroError> {
        match self {
            ItemResult::Ok(_) => None,
            ItemResult::Err { error, .. } => Some(error),
        }
    }
}

/// Result of a batch operation.
#[derive(Debug)]
pub struct BatchResult<T> {
    /// Individual results for each item.
    pub results: Vec<ItemResult<T>>,
    /// Total processing time.
    pub duration: Duration,
}

impl<T> BatchResult<T> {
    /// Create a new batch result.
    pub fn new(results: Vec<ItemResult<T>>, duration: Duration) -> Self {
        Self { results, duration }
    }

    /// Get the total number of items processed.
    pub fn total(&self) -> usize {
        self.results.len()
    }

    /// Get the number of successful items.
    pub fn success_count(&self) -> usize {
        self.results.iter().filter(|r| r.is_ok()).count()
    }

    /// Get the number of failed items.
    pub fn error_count(&self) -> usize {
        self.results.iter().filter(|r| r.is_err()).count()
    }

    /// Calculate success rate as a percentage.
    pub fn success_rate(&self) -> f64 {
        if self.results.is_empty() {
            100.0
        } else {
            (self.success_count() as f64 / self.results.len() as f64) * 100.0
        }
    }

    /// Calculate processing rate (items per second).
    ///
    /// Returns 0.0 if the duration is too short to provide a meaningful rate.
    pub fn items_per_second(&self) -> f64 {
        let secs = self.duration.as_secs_f64();
        if secs < f64::EPSILON {
            0.0
        } else {
            self.results.len() as f64 / secs
        }
    }

    /// Get only successful results.
    pub fn successes(self) -> Vec<T> {
        self.results.into_iter().filter_map(|r| r.ok()).collect()
    }

    /// Get only errors.
    pub fn errors(self) -> Vec<FerroError> {
        self.results.into_iter().filter_map(|r| r.err()).collect()
    }

    /// Check if all items were successful.
    pub fn all_ok(&self) -> bool {
        self.results.iter().all(|r| r.is_ok())
    }

    /// Check if any items failed.
    pub fn has_errors(&self) -> bool {
        self.results.iter().any(|r| r.is_err())
    }
}

impl<T: Clone> BatchResult<T> {
    /// Get successful results as references.
    pub fn success_refs(&self) -> Vec<&T> {
        self.results
            .iter()
            .filter_map(|r| match r {
                ItemResult::Ok(v) => Some(v),
                ItemResult::Err { .. } => None,
            })
            .collect()
    }
}

/// Batch processor for HGVS variants.
///
/// Provides high-level APIs for batch parsing and normalization
/// with progress tracking and error aggregation.
pub struct BatchProcessor<P: ReferenceProvider> {
    normalizer: Normalizer<P>,
    config: BatchConfig,
}

impl<P: ReferenceProvider> BatchProcessor<P> {
    /// Create a new batch processor.
    ///
    /// # Arguments
    ///
    /// * `provider` - Reference sequence provider
    pub fn new(provider: P) -> Self {
        Self {
            normalizer: Normalizer::new(provider),
            config: BatchConfig::default(),
        }
    }

    /// Create a new batch processor with configuration.
    ///
    /// # Arguments
    ///
    /// * `provider` - Reference sequence provider
    /// * `config` - Batch processing configuration
    pub fn with_config(provider: P, config: BatchConfig) -> Self {
        Self {
            normalizer: Normalizer::new(provider),
            config,
        }
    }

    /// Get the current configuration.
    pub fn config(&self) -> &BatchConfig {
        &self.config
    }

    /// Parse multiple HGVS strings.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of HGVS strings to parse
    ///
    /// # Returns
    ///
    /// Batch result containing parsed variants or errors.
    pub fn parse<S: AsRef<str>>(&self, variants: &[S]) -> BatchResult<HgvsVariant> {
        self.parse_with_progress(variants, |_| {})
    }

    /// Parse multiple HGVS strings with progress callback.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of HGVS strings to parse
    /// * `progress_fn` - Callback function called with progress updates
    pub fn parse_with_progress<S, F>(
        &self,
        variants: &[S],
        mut progress_fn: F,
    ) -> BatchResult<HgvsVariant>
    where
        S: AsRef<str>,
        F: FnMut(BatchProgress),
    {
        let start = Instant::now();
        let total = variants.len();
        let mut results = Vec::with_capacity(total);
        let mut success = 0;
        let mut errors = 0;

        for (i, input) in variants.iter().enumerate() {
            match parse_hgvs(input.as_ref()) {
                Ok(variant) => {
                    results.push(ItemResult::Ok(variant));
                    success += 1;
                }
                Err(error) => {
                    results.push(ItemResult::Err {
                        input: Some(input.as_ref().to_string()),
                        error,
                    });
                    errors += 1;
                }
            }

            // Progress callback
            if (i + 1) % self.config.progress_interval == 0 || i + 1 == total {
                progress_fn(BatchProgress {
                    total,
                    processed: i + 1,
                    success,
                    errors,
                    elapsed: start.elapsed(),
                });
            }
        }

        BatchResult::new(results, start.elapsed())
    }

    /// Normalize multiple variants.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to normalize
    pub fn normalize(&self, variants: &[HgvsVariant]) -> BatchResult<HgvsVariant> {
        self.normalize_with_progress(variants, |_| {})
    }

    /// Normalize multiple variants with progress callback.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of variants to normalize
    /// * `progress_fn` - Callback function called with progress updates
    pub fn normalize_with_progress<F>(
        &self,
        variants: &[HgvsVariant],
        mut progress_fn: F,
    ) -> BatchResult<HgvsVariant>
    where
        F: FnMut(BatchProgress),
    {
        let start = Instant::now();
        let total = variants.len();
        let mut results = Vec::with_capacity(total);
        let mut success = 0;
        let mut errors = 0;

        for (i, variant) in variants.iter().enumerate() {
            match self.normalizer.normalize(variant) {
                Ok(normalized) => {
                    results.push(ItemResult::Ok(normalized));
                    success += 1;
                }
                Err(error) => {
                    results.push(ItemResult::Err {
                        input: Some(variant.to_string()),
                        error,
                    });
                    errors += 1;
                }
            }

            // Progress callback
            if (i + 1) % self.config.progress_interval == 0 || i + 1 == total {
                progress_fn(BatchProgress {
                    total,
                    processed: i + 1,
                    success,
                    errors,
                    elapsed: start.elapsed(),
                });
            }
        }

        BatchResult::new(results, start.elapsed())
    }

    /// Parse and normalize in one operation.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of HGVS strings to parse and normalize
    pub fn parse_and_normalize<S: AsRef<str>>(&self, variants: &[S]) -> BatchResult<HgvsVariant> {
        self.parse_and_normalize_with_progress(variants, |_| {})
    }

    /// Parse and normalize with progress callback.
    ///
    /// # Arguments
    ///
    /// * `variants` - Slice of HGVS strings to parse and normalize
    /// * `progress_fn` - Callback function called with progress updates
    pub fn parse_and_normalize_with_progress<S, F>(
        &self,
        variants: &[S],
        mut progress_fn: F,
    ) -> BatchResult<HgvsVariant>
    where
        S: AsRef<str>,
        F: FnMut(BatchProgress),
    {
        let start = Instant::now();
        let total = variants.len();
        let mut results = Vec::with_capacity(total);
        let mut success = 0;
        let mut errors = 0;

        for (i, input) in variants.iter().enumerate() {
            let result = parse_hgvs(input.as_ref()).and_then(|v| self.normalizer.normalize(&v));

            match result {
                Ok(normalized) => {
                    results.push(ItemResult::Ok(normalized));
                    success += 1;
                }
                Err(error) => {
                    results.push(ItemResult::Err {
                        input: Some(input.as_ref().to_string()),
                        error,
                    });
                    errors += 1;
                }
            }

            // Progress callback
            if (i + 1) % self.config.progress_interval == 0 || i + 1 == total {
                progress_fn(BatchProgress {
                    total,
                    processed: i + 1,
                    success,
                    errors,
                    elapsed: start.elapsed(),
                });
            }
        }

        BatchResult::new(results, start.elapsed())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MockProvider;

    fn processor() -> BatchProcessor<MockProvider> {
        BatchProcessor::new(MockProvider::with_test_data())
    }

    #[test]
    fn test_batch_config_default() {
        let config = BatchConfig::default();
        assert!(config.continue_on_error);
        assert_eq!(config.progress_interval, 100);
    }

    #[test]
    fn test_batch_config_builder() {
        let config = BatchConfig::new()
            .continue_on_error(false)
            .progress_interval(50);

        assert!(!config.continue_on_error);
        assert_eq!(config.progress_interval, 50);
    }

    #[test]
    fn test_batch_progress_percent() {
        let progress = BatchProgress {
            total: 100,
            processed: 50,
            success: 45,
            errors: 5,
            elapsed: Duration::from_secs(1),
        };
        assert!((progress.percent() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_batch_progress_items_per_second() {
        let progress = BatchProgress {
            total: 100,
            processed: 50,
            success: 45,
            errors: 5,
            elapsed: Duration::from_secs(1),
        };
        assert!((progress.items_per_second() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_batch_progress_estimated_remaining() {
        let progress = BatchProgress {
            total: 100,
            processed: 50,
            success: 45,
            errors: 5,
            elapsed: Duration::from_secs(1),
        };
        let remaining = progress.estimated_remaining().unwrap();
        // 50 items remaining at 50 items/sec = 1 second
        assert!((remaining.as_secs_f64() - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_item_result_ok() {
        let result: ItemResult<HgvsVariant> =
            ItemResult::Ok(parse_hgvs("NM_000088.3:c.10A>G").unwrap());
        assert!(result.is_ok());
        assert!(!result.is_err());
    }

    #[test]
    fn test_item_result_err() {
        let result: ItemResult<HgvsVariant> = ItemResult::Err {
            input: Some("invalid".to_string()),
            error: FerroError::parse(0, "test error"),
        };
        assert!(!result.is_ok());
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_batch() {
        let processor = processor();
        let variants = vec!["NM_000088.3:c.10A>G", "NC_000001.11:g.12345A>G"];

        let result = processor.parse(&variants);
        assert_eq!(result.total(), 2);
        assert_eq!(result.success_count(), 2);
        assert_eq!(result.error_count(), 0);
        assert!(result.all_ok());
    }

    #[test]
    fn test_parse_batch_with_errors() {
        let processor = processor();
        let variants = vec![
            "NM_000088.3:c.10A>G",
            "invalid variant",
            "NC_000001.11:g.12345A>G",
        ];

        let result = processor.parse(&variants);
        assert_eq!(result.total(), 3);
        assert_eq!(result.success_count(), 2);
        assert_eq!(result.error_count(), 1);
        assert!(result.has_errors());
    }

    #[test]
    fn test_parse_batch_with_progress() {
        let processor = processor();
        let variants = vec!["NM_000088.3:c.10A>G", "NC_000001.11:g.12345A>G"];

        let mut progress_count = 0;
        let result = processor.parse_with_progress(&variants, |_| {
            progress_count += 1;
        });

        assert_eq!(result.total(), 2);
        assert!(progress_count > 0);
    }

    #[test]
    fn test_normalize_batch() {
        let processor = processor();
        let variants: Vec<_> = vec!["NM_000088.3:c.10A>G", "NC_000001.11:g.12345A>G"]
            .into_iter()
            .map(|s| parse_hgvs(s).unwrap())
            .collect();

        let result = processor.normalize(&variants);
        assert_eq!(result.total(), 2);
    }

    #[test]
    fn test_parse_and_normalize_batch() {
        let processor = processor();
        let variants = vec!["NM_000088.3:c.10A>G", "NC_000001.11:g.12345A>G"];

        let result = processor.parse_and_normalize(&variants);
        assert_eq!(result.total(), 2);
    }

    #[test]
    fn test_batch_result_successes() {
        let processor = processor();
        let variants = vec!["NM_000088.3:c.10A>G", "invalid", "NC_000001.11:g.12345A>G"];

        let result = processor.parse(&variants);
        let successes = result.successes();
        assert_eq!(successes.len(), 2);
    }

    #[test]
    fn test_batch_result_errors() {
        let processor = processor();
        let variants = vec!["NM_000088.3:c.10A>G", "invalid", "NC_000001.11:g.12345A>G"];

        let result = processor.parse(&variants);
        let errors = result.errors();
        assert_eq!(errors.len(), 1);
    }

    #[test]
    fn test_batch_result_success_rate() {
        let processor = processor();
        let variants = vec!["NM_000088.3:c.10A>G", "invalid"];

        let result = processor.parse(&variants);
        assert!((result.success_rate() - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_empty_batch() {
        let processor = processor();
        let variants: Vec<&str> = vec![];

        let result = processor.parse(&variants);
        assert_eq!(result.total(), 0);
        assert!(result.all_ok());
        assert!((result.success_rate() - 100.0).abs() < 0.01);
    }
}
