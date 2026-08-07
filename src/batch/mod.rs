//! Batch processing for HGVS variant operations.
//!
//! This module provides high-level APIs for batch parsing and normalization
//! of multiple HGVS variants, with progress tracking, error aggregation,
//! and detailed statistics.
//!
//! # Examples
//!
//! ## Basic Usage
//!
//! ```
//! use ferro_hgvs::batch::{BatchProcessor, BatchResult};
//! use ferro_hgvs::MockProvider;
//!
//! let provider = MockProvider::with_test_data();
//! let processor = BatchProcessor::new(provider);
//!
//! let variants = vec![
//!     "NM_000088.3:c.10A>G",
//!     "NC_000001.11:g.12345A>G",
//! ];
//!
//! let result = processor.parse(&variants);
//! println!("Parsed {}/{} variants", result.success_count(), result.total());
//! ```
//!
//! ## With Progress Callback
//!
//! ```
//! use ferro_hgvs::batch::{BatchProcessor, BatchConfig};
//! use ferro_hgvs::MockProvider;
//!
//! let provider = MockProvider::with_test_data();
//! let processor = BatchProcessor::new(provider);
//!
//! let variants = vec!["NM_000088.3:c.10A>G"];
//!
//! let result = processor.parse_with_progress(&variants, |progress| {
//!     println!("Progress: {:.1}%", progress.percent());
//! });
//! ```
//!
//! # Features
//!
//! - **Flexible Processing**: Parse-only, normalize-only, or combined operations
//! - **Progress Tracking**: Optional callbacks for monitoring long-running batches
//! - **Error Aggregation**: Collect all errors without stopping the batch
//! - **Statistics**: Detailed success/failure counts and processing rates
//! - **Parallel Support**: multi-threaded processing via the `*_parallel` and
//!   `*_streaming` methods, available in every build

mod processor;

// Ungated, matching its consumers: the streaming impl block in `processor.rs`
// and the streaming functions in `crate::parallel`. Both were once
// `#[cfg(feature = "parallel")]` and this followed them, which made the module
// dead in a featureless build — dead code that `dead_code`'s warn-by-default
// level let through every clippy job, since they all pass `--features dev` or
// `--all-features`. #1507 removed those gates (rayon is unconditional, so they
// gated visibility rather than capability), which removes the dead-code
// condition at its source.
pub(crate) mod streaming;

pub use processor::{
    BatchConfig, BatchProcessor, BatchProgress, BatchResult, ItemResult, BATCH_CHUNK_ITEMS,
    STREAM_CHUNK_ITEMS,
};
