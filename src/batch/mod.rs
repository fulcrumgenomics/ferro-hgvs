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
//! - **Parallel Support**: Enable `parallel` feature for multi-threaded processing

mod processor;

pub use processor::{BatchConfig, BatchProcessor, BatchProgress, BatchResult, ItemResult};
