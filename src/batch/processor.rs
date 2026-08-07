//! Batch processor implementation.

use std::time::{Duration, Instant};

use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::Normalizer;
use crate::reference::ReferenceProvider;

/// Items processed per chunk by the CLI's batch path.
///
/// The value `run_batch` has used since #687, exported so the CLI and the library
/// share one definition rather than two that can drift. #975 asks for that
/// chunk-and-preserve-order engine to be reused, and it is — see
/// [`STREAM_CHUNK_ITEMS`] for why the streaming API chunks at a different size.
pub const BATCH_CHUNK_ITEMS: usize = 8192;

/// Items read, mapped and drained per chunk by the streaming batch API (#975).
///
/// Larger than [`BATCH_CHUNK_ITEMS`], and sized by measurement rather than
/// inherited. #975 asks for two things at once — peak RSS bounded and flat in
/// input size, *and* throughput within noise of the `Vec`-based path — and 8192
/// delivers only the first. Measured on 1M genomic substitutions, three reps:
///
/// | chunk | peak RSS @ 4M inputs | wall clock vs `Vec` |
/// |---|---|---|
/// | 8 192 | 17.6 MiB | +8 % (0.179 s against 0.166 s) |
/// | 65 536 | 82.9 MiB | parity (0.163 s against 0.164 s) |
///
/// At 8192 the per-chunk overhead — refill, a result `Vec` per chunk, and rayon
/// fanning out over only 8192 items at a time — is a measurable fraction of the
/// per-item work, which for a bare parse is very small. The CLI does not see this
/// because its per-item work (normalize plus formatting plus a write) dwarfs the
/// chunk overhead, which is exactly why the two constants differ instead of one
/// being bent to serve both.
///
/// 65 536 keeps the property that matters: RSS measured 81.3 / 90.6 / 82.9 /
/// 84.1 MiB at 250k / 1M / 4M / 8M inputs — flat across a 32x range, against the
/// `Vec` path's 157 MiB / 610 MiB / 2419 MiB over the first three.
pub const STREAM_CHUNK_ITEMS: usize = 64 * 1024;

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

/// Parallel batch processing. Each input is independent, so the batch is mapped
/// across a rayon pool; `rayon`'s ordered `collect` preserves input order, so
/// results are identical to the serial methods (same order, same per-item
/// success/error classification). Requires the provider to be `Sync`.
///
/// `num_threads` semantics: `0` uses rayon's global pool (all available cores),
/// `1` runs serially (no pool, no fan-out overhead), and `N > 1` runs on a
/// dedicated pool of `N` threads. These have no progress-callback variants —
/// ordered progress under parallelism is ill-defined; use the serial
/// `*_with_progress` methods when a callback is needed.
#[cfg(feature = "parallel")]
impl<P: ReferenceProvider + Sync> BatchProcessor<P> {
    /// Parse multiple HGVS strings in parallel. Order-preserving; equivalent to
    /// [`parse`](Self::parse) but spread across `num_threads` workers.
    pub fn parse_parallel<S: AsRef<str> + Sync>(
        &self,
        variants: &[S],
        num_threads: usize,
    ) -> BatchResult<HgvsVariant> {
        let start = Instant::now();
        let map = |input: &S| match parse_hgvs(input.as_ref()) {
            Ok(v) => ItemResult::Ok(v),
            Err(error) => ItemResult::Err {
                input: Some(input.as_ref().to_string()),
                error,
            },
        };
        let results = self.map_collect(variants, num_threads, map);
        BatchResult::new(results, start.elapsed())
    }

    /// Parse and normalize multiple HGVS strings in parallel. Order-preserving;
    /// equivalent to [`parse_and_normalize`](Self::parse_and_normalize) but
    /// spread across `num_threads` workers.
    pub fn parse_and_normalize_parallel<S: AsRef<str> + Sync>(
        &self,
        variants: &[S],
        num_threads: usize,
    ) -> BatchResult<HgvsVariant> {
        let start = Instant::now();
        let map = |input: &S| match parse_hgvs(input.as_ref())
            .and_then(|v| self.normalizer.normalize(&v))
        {
            Ok(normalized) => ItemResult::Ok(normalized),
            Err(error) => ItemResult::Err {
                input: Some(input.as_ref().to_string()),
                error,
            },
        };
        let results = self.map_collect(variants, num_threads, map);
        BatchResult::new(results, start.elapsed())
    }

    /// Parse a *stream* of HGVS strings, yielding results in input order (#975).
    ///
    /// The streaming counterpart to [`parse_parallel`](Self::parse_parallel), and
    /// the reason it exists: that method takes `&[S]` and returns a `Vec`, so both
    /// the whole input and the whole output are resident at once by contract. No
    /// amount of internal chunking can bound peak memory behind such a signature
    /// — #975 says so explicitly — which is why this takes an `IntoIterator` and
    /// returns an `Iterator` instead.
    ///
    /// Memory is bounded by one chunk of inputs plus one chunk of results, not by
    /// the length of the stream. Order is preserved exactly as the `Vec`-based
    /// path preserves it, because each chunk is mapped with the same
    /// order-preserving `collect` and drained before the next is read.
    ///
    /// The existing `Vec`-based methods are untouched.
    pub fn parse_streaming<'a, I, S>(
        &'a self,
        variants: I,
        num_threads: usize,
    ) -> impl Iterator<Item = ItemResult<HgvsVariant>> + 'a
    where
        I: IntoIterator<Item = S> + 'a,
        S: AsRef<str> + Send + Sync + 'a,
    {
        self.map_streaming(variants, num_threads, |input: &S| {
            match parse_hgvs(input.as_ref()) {
                Ok(variant) => ItemResult::Ok(variant),
                Err(error) => ItemResult::Err {
                    input: Some(input.as_ref().to_string()),
                    error,
                },
            }
        })
    }

    /// Parse and normalize a *stream* of HGVS strings, yielding results in input
    /// order (#975). See [`parse_streaming`](Self::parse_streaming) for why the
    /// signature is shaped this way.
    pub fn parse_and_normalize_streaming<'a, I, S>(
        &'a self,
        variants: I,
        num_threads: usize,
    ) -> impl Iterator<Item = ItemResult<HgvsVariant>> + 'a
    where
        I: IntoIterator<Item = S> + 'a,
        S: AsRef<str> + Send + Sync + 'a,
    {
        self.map_streaming(variants, num_threads, |input: &S| {
            match parse_hgvs(input.as_ref()).and_then(|v| self.normalizer.normalize(&v)) {
                Ok(normalized) => ItemResult::Ok(normalized),
                Err(error) => ItemResult::Err {
                    input: Some(input.as_ref().to_string()),
                    error,
                },
            }
        })
    }

    /// Shared engine behind the streaming methods: read a bounded chunk, map it
    /// order-preservingly, drain it, repeat.
    ///
    /// The engine is [`crate::batch::streaming::StreamingMap`], shared with
    /// [`crate::parallel`]'s streaming functions so there is one chunk-and-drain
    /// loop rather than two. Its shape is the CLI's `run_batch`/`flush_chunk`
    /// pair (#687), which #975 asks to reuse.
    fn map_streaming<'a, I, S, F>(
        &'a self,
        variants: I,
        num_threads: usize,
        f: F,
    ) -> impl Iterator<Item = ItemResult<HgvsVariant>> + 'a
    where
        I: IntoIterator<Item = S> + 'a,
        S: AsRef<str> + Send + Sync + 'a,
        F: Fn(&S) -> ItemResult<HgvsVariant> + Sync + Send + 'a,
    {
        crate::batch::streaming::StreamingMap::new(variants, num_threads, STREAM_CHUNK_ITEMS, f)
    }

    /// Map `f` over `variants`, preserving order. Serial for `num_threads == 1`,
    /// the global rayon pool for `0`, or a dedicated `num_threads`-thread pool
    /// otherwise (falling back to the global pool if one can't be built).
    fn map_collect<S, F>(
        &self,
        variants: &[S],
        num_threads: usize,
        f: F,
    ) -> Vec<ItemResult<HgvsVariant>>
    where
        S: AsRef<str> + Sync,
        F: Fn(&S) -> ItemResult<HgvsVariant> + Sync + Send,
    {
        use rayon::prelude::*;
        if num_threads == 1 {
            return variants.iter().map(f).collect();
        }
        if num_threads == 0 {
            return variants.par_iter().map(f).collect();
        }
        match rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
        {
            Ok(pool) => pool.install(|| variants.par_iter().map(f).collect()),
            Err(err) => {
                log::warn!(
                    "failed to build a dedicated {num_threads}-thread pool ({err}); \
                     falling back to the global rayon pool"
                );
                variants.par_iter().map(f).collect()
            }
        }
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

    /// The streaming batch methods must agree item-for-item with the `Vec`-based
    /// ones (#975), across the thread settings that select different engines
    /// (serial, global pool, dedicated pool).
    #[cfg(feature = "parallel")]
    #[test]
    fn streaming_batch_matches_the_vec_api_exactly() {
        let processor = processor();
        let inputs: Vec<String> = (0..BATCH_CHUNK_ITEMS + 11)
            .map(|i| {
                match i % 4 {
                    0 => "NM_000088.3:c.459A>G",
                    1 => "NC_000001.11:g.12345A>G",
                    2 => "definitely not hgvs",
                    _ => "NP_000079.2:p.Val600Glu",
                }
                .to_string()
            })
            .collect();
        let render = |r: &ItemResult<HgvsVariant>| match r {
            ItemResult::Ok(v) => format!("ok:{v}"),
            ItemResult::Err { input, error } => format!("err:{input:?}:{error}"),
        };
        for threads in [1usize, 0, 2] {
            let expected: Vec<String> = processor
                .parse_parallel(&inputs, threads)
                .results
                .iter()
                .map(render)
                .collect();
            let streamed: Vec<String> = processor
                .parse_streaming(inputs.iter(), threads)
                .map(|r| render(&r))
                .collect();
            assert_eq!(streamed, expected, "parse, num_threads={threads}");

            let expected: Vec<String> = processor
                .parse_and_normalize_parallel(&inputs, threads)
                .results
                .iter()
                .map(render)
                .collect();
            let streamed: Vec<String> = processor
                .parse_and_normalize_streaming(inputs.iter(), threads)
                .map(|r| render(&r))
                .collect();
            assert_eq!(streamed, expected, "parse+normalize, num_threads={threads}");
        }
    }
}
