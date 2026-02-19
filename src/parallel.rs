//! Parallel processing support for ferro-hgvs
//!
//! This module provides parallel variants of parsing and normalization
//! operations using rayon. Enable with the `parallel` feature.
//!
//! # Example
//!
//! ```no_run
//! # #[cfg(feature = "parallel")]
//! # fn main() {
//! use ferro_hgvs::parallel::{parse_hgvs_parallel, normalize_parallel};
//! use ferro_hgvs::{MockProvider, Normalizer};
//!
//! let variants = vec![
//!     "NM_000088.3:c.459A>G",
//!     "NC_000001.11:g.12345A>G",
//!     "NP_000079.2:p.Val600Glu",
//! ];
//!
//! // Parse in parallel
//! let parsed: Vec<_> = parse_hgvs_parallel(&variants)
//!     .into_iter()
//!     .filter_map(|r| r.ok())
//!     .collect();
//!
//! // Normalize in parallel (requires a provider)
//! let provider = MockProvider::with_test_data();
//! let normalizer = Normalizer::new(provider);
//! let _normalized = normalize_parallel(&normalizer, &parsed);
//! # }
//! # #[cfg(not(feature = "parallel"))]
//! # fn main() {}
//! ```

use rayon::prelude::*;

use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::Normalizer;
use crate::reference::ReferenceProvider;

/// Parse multiple HGVS variant strings in parallel
///
/// Returns a vector of results, one for each input string.
/// Order is preserved.
pub fn parse_hgvs_parallel<S: AsRef<str> + Sync>(
    variants: &[S],
) -> Vec<Result<HgvsVariant, FerroError>> {
    variants
        .par_iter()
        .map(|s| parse_hgvs(s.as_ref()))
        .collect()
}

/// Parse multiple HGVS variant strings in parallel, filtering errors
///
/// Returns only successfully parsed variants. Useful when you want to
/// skip invalid variants without error handling.
pub fn parse_hgvs_parallel_ok<S: AsRef<str> + Sync>(variants: &[S]) -> Vec<HgvsVariant> {
    variants
        .par_iter()
        .filter_map(|s| parse_hgvs(s.as_ref()).ok())
        .collect()
}

/// Normalize multiple variants in parallel
///
/// Returns a vector of results, one for each input variant.
/// Order is preserved.
pub fn normalize_parallel<P: ReferenceProvider + Sync>(
    normalizer: &Normalizer<P>,
    variants: &[HgvsVariant],
) -> Vec<Result<HgvsVariant, FerroError>> {
    variants
        .par_iter()
        .map(|v| normalizer.normalize(v))
        .collect()
}

/// Normalize multiple variants in parallel, filtering errors
///
/// Returns only successfully normalized variants.
pub fn normalize_parallel_ok<P: ReferenceProvider + Sync>(
    normalizer: &Normalizer<P>,
    variants: &[HgvsVariant],
) -> Vec<HgvsVariant> {
    variants
        .par_iter()
        .filter_map(|v| normalizer.normalize(v).ok())
        .collect()
}

/// Parse and normalize in a single parallel operation
///
/// Parses each string and normalizes it. Returns results for all operations.
pub fn parse_and_normalize_parallel<P: ReferenceProvider + Sync, S: AsRef<str> + Sync>(
    normalizer: &Normalizer<P>,
    variants: &[S],
) -> Vec<Result<HgvsVariant, FerroError>> {
    variants
        .par_iter()
        .map(|s| {
            let variant = parse_hgvs(s.as_ref())?;
            normalizer.normalize(&variant)
        })
        .collect()
}

/// Configuration for parallel batch processing
#[derive(Debug, Clone)]
pub struct ParallelConfig {
    /// Chunk size for parallel processing
    pub chunk_size: usize,
    /// Number of threads (0 = use rayon default)
    pub num_threads: usize,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        Self {
            chunk_size: 1000,
            num_threads: 0,
        }
    }
}

impl ParallelConfig {
    /// Create a new parallel configuration
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the chunk size for batched processing
    pub fn with_chunk_size(mut self, size: usize) -> Self {
        self.chunk_size = size;
        self
    }

    /// Set the number of threads
    pub fn with_num_threads(mut self, threads: usize) -> Self {
        self.num_threads = threads;
        self
    }
}

/// Statistics from parallel processing
#[derive(Debug, Clone, Default)]
pub struct ParallelStats {
    /// Total items processed
    pub total: usize,
    /// Successfully processed
    pub success: usize,
    /// Failed to process
    pub errors: usize,
}

impl ParallelStats {
    /// Create new stats
    pub fn new() -> Self {
        Self::default()
    }

    /// Calculate success rate as a percentage
    pub fn success_rate(&self) -> f64 {
        if self.total == 0 {
            0.0
        } else {
            (self.success as f64 / self.total as f64) * 100.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MockProvider;

    #[test]
    fn test_parse_parallel() {
        let variants = vec![
            "NM_000088.3:c.459A>G",
            "NC_000001.11:g.12345A>G",
            "NP_000079.2:p.Val600Glu",
        ];

        let results = parse_hgvs_parallel(&variants);
        assert_eq!(results.len(), 3);
        assert!(results.iter().all(|r| r.is_ok()));
    }

    #[test]
    fn test_parse_parallel_ok() {
        let variants = vec![
            "NM_000088.3:c.459A>G",
            "invalid variant",
            "NC_000001.11:g.12345A>G",
        ];

        let results = parse_hgvs_parallel_ok(&variants);
        assert_eq!(results.len(), 2); // Only valid variants
    }

    #[test]
    fn test_normalize_parallel() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variants: Vec<HgvsVariant> = vec![
            parse_hgvs("NM_000088.3:c.10A>G").unwrap(),
            parse_hgvs("NC_000001.11:g.12345A>G").unwrap(),
        ];

        let results = normalize_parallel(&normalizer, &variants);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_parse_and_normalize_parallel() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variants = vec!["NM_000088.3:c.10A>G", "NC_000001.11:g.12345A>G"];

        let results = parse_and_normalize_parallel(&normalizer, &variants);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_parallel_stats() {
        let stats = ParallelStats {
            total: 100,
            success: 95,
            errors: 5,
        };

        assert!((stats.success_rate() - 95.0).abs() < 0.01);
    }

    // Parallel stress tests

    #[test]
    fn test_stress_parse_1000_variants() {
        // Generate 1000 variants
        let variants: Vec<String> = (1..=1000)
            .map(|i| format!("NM_000088.3:c.{}A>G", i))
            .collect();

        let results = parse_hgvs_parallel(&variants);
        assert_eq!(results.len(), 1000);
        assert!(results.iter().all(|r| r.is_ok()));
    }

    #[test]
    fn test_stress_parse_with_mixed_errors() {
        // Generate mix of valid and invalid variants
        let variants: Vec<String> = (1..=500)
            .flat_map(|i| {
                vec![
                    format!("NM_000088.3:c.{}A>G", i), // valid
                    format!("invalid_variant_{}", i),  // invalid
                ]
            })
            .collect();

        let results = parse_hgvs_parallel(&variants);
        assert_eq!(results.len(), 1000);

        let successes = results.iter().filter(|r| r.is_ok()).count();
        let errors = results.iter().filter(|r| r.is_err()).count();
        assert_eq!(successes, 500);
        assert_eq!(errors, 500);
    }

    #[test]
    fn test_stress_parse_order_preserved() {
        // Verify order is maintained in parallel processing
        let variants: Vec<String> = (1..=100)
            .map(|i| format!("NM_000088.3:c.{}A>G", i))
            .collect();

        let results = parse_hgvs_parallel(&variants);

        // Check that position numbers are in order
        for (i, result) in results.iter().enumerate() {
            let variant = result.as_ref().unwrap();
            let expected_pos = i + 1;
            // Verify the variant string contains expected position
            let s = variant.to_string();
            assert!(s.contains(&format!("c.{}A>G", expected_pos)));
        }
    }

    #[test]
    fn test_stress_normalize_500_variants() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variants: Vec<HgvsVariant> = (1..=500)
            .map(|i| parse_hgvs(&format!("NM_000088.3:c.{}A>G", i)).unwrap())
            .collect();

        let results = normalize_parallel(&normalizer, &variants);
        assert_eq!(results.len(), 500);
    }

    #[test]
    fn test_stress_parse_and_normalize_combined() {
        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variants: Vec<String> = (1..=500)
            .map(|i| format!("NM_000088.3:c.{}A>G", i))
            .collect();

        let results = parse_and_normalize_parallel(&normalizer, &variants);
        assert_eq!(results.len(), 500);
    }

    #[test]
    fn test_stress_diverse_variant_types() {
        // Test parallel parsing of diverse variant types
        let variants = vec![
            // Substitutions at various positions
            "NM_000088.3:c.100A>G",
            "NM_000088.3:c.200C>T",
            "NM_000088.3:c.300G>A",
            // Deletions
            "NM_000088.3:c.100del",
            "NM_000088.3:c.100_102del",
            // Insertions
            "NM_000088.3:c.100_101insATG",
            // Duplications
            "NM_000088.3:c.100dup",
            "NM_000088.3:c.100_102dup",
            // Delins
            "NM_000088.3:c.100delinsGGG",
            // Genomic variants
            "NC_000001.11:g.12345A>G",
            "NC_000001.11:g.12345_12350del",
            // Protein variants
            "NP_000079.2:p.Val600Glu",
            "NP_000079.2:p.Ala100Ter",
        ];

        // Repeat each 100 times
        let all_variants: Vec<&str> = variants
            .iter()
            .cycle()
            .take(variants.len() * 100)
            .copied()
            .collect();

        let results = parse_hgvs_parallel(&all_variants);
        assert_eq!(results.len(), variants.len() * 100);
        assert!(results.iter().all(|r| r.is_ok()));
    }

    #[test]
    fn test_stress_empty_input() {
        let variants: Vec<&str> = vec![];
        let results = parse_hgvs_parallel(&variants);
        assert!(results.is_empty());
    }

    #[test]
    fn test_stress_single_item() {
        let variants = vec!["NM_000088.3:c.10A>G"];
        let results = parse_hgvs_parallel(&variants);
        assert_eq!(results.len(), 1);
        assert!(results[0].is_ok());
    }

    #[test]
    fn test_stress_all_errors() {
        let variants: Vec<String> = (1..=100)
            .map(|i| format!("invalid_variant_{}", i))
            .collect();

        let results = parse_hgvs_parallel(&variants);
        assert_eq!(results.len(), 100);
        assert!(results.iter().all(|r| r.is_err()));
    }

    #[test]
    fn test_parallel_config_variations() {
        // Test different configurations
        let config1 = ParallelConfig::new()
            .with_chunk_size(100)
            .with_num_threads(2);
        assert_eq!(config1.chunk_size, 100);
        assert_eq!(config1.num_threads, 2);

        let config2 = ParallelConfig::new()
            .with_chunk_size(10000)
            .with_num_threads(0);
        assert_eq!(config2.chunk_size, 10000);
        assert_eq!(config2.num_threads, 0);
    }

    #[test]
    fn test_parallel_stats_edge_cases() {
        // Empty stats
        let empty = ParallelStats::new();
        assert_eq!(empty.total, 0);
        assert_eq!(empty.success_rate(), 0.0);

        // All success
        let all_success = ParallelStats {
            total: 100,
            success: 100,
            errors: 0,
        };
        assert!((all_success.success_rate() - 100.0).abs() < 0.01);

        // All errors
        let all_errors = ParallelStats {
            total: 100,
            success: 0,
            errors: 100,
        };
        assert!((all_errors.success_rate()).abs() < 0.01);
    }

    #[test]
    fn test_stress_concurrent_throughput() {
        use std::time::Instant;

        // Generate a large number of variants
        let variants: Vec<String> = (1..=2000)
            .map(|i| format!("NM_000088.3:c.{}A>G", i))
            .collect();

        let start = Instant::now();
        let results = parse_hgvs_parallel(&variants);
        let duration = start.elapsed();

        assert_eq!(results.len(), 2000);
        assert!(results.iter().all(|r| r.is_ok()));

        // Just verify it completes in reasonable time (< 5 seconds)
        assert!(duration.as_secs() < 5);
    }
}
