//! Parallel processing support for ferro-hgvs
//!
//! This module provides parallel variants of parsing and normalization
//! operations using rayon. It is available in every build — rayon is an
//! unconditional dependency, so there is no feature to enable.
//!
//! # Example
//!
//! ```no_run
//! use ferro_hgvs::parallel::{parse_hgvs_parallel, normalize_parallel};
//! use ferro_hgvs::{JsonProvider, Normalizer};
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
//! let provider = JsonProvider::with_test_data();
//! let normalizer = Normalizer::new(provider);
//! let _normalized = normalize_parallel(&normalizer, &parsed);
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

/// Parse a *stream* of HGVS strings in parallel, yielding results in input order
/// (#975).
///
/// The streaming counterpart to [`parse_hgvs_parallel`]. That function takes
/// `&[S]` and returns a `Vec`, so the whole input and the whole output are
/// resident at once by contract — no internal chunking can bound peak memory
/// behind such a signature, which is #975's central point. This takes an
/// `IntoIterator` and yields as it goes, so memory is a function of the chunk
/// size rather than of the stream's length.
///
/// Order is preserved, and each item is the same `Result` the `Vec`-based
/// function would have produced at that position.
///
/// ```no_run
/// # fn main() -> std::io::Result<()> {
/// use std::io::{BufRead, BufReader};
/// use ferro_hgvs::parallel::parse_hgvs_streaming;
///
/// let file = BufReader::new(std::fs::File::open("variants.txt")?);
/// // The file's lines are never all in memory at once.
/// for result in parse_hgvs_streaming(file.lines().map_while(Result::ok), 0) {
///     if let Ok(variant) = result {
///         println!("{variant}");
///     }
/// }
/// # Ok(())
/// # }
/// ```
pub fn parse_hgvs_streaming<I, S>(
    variants: I,
    num_threads: usize,
) -> impl Iterator<Item = Result<HgvsVariant, FerroError>>
where
    I: IntoIterator<Item = S>,
    S: AsRef<str> + Send + Sync,
{
    crate::batch::streaming::StreamingMap::new(
        variants,
        num_threads,
        crate::batch::STREAM_CHUNK_ITEMS,
        |s: &S| parse_hgvs(s.as_ref()),
    )
}

/// Parse and normalize a *stream* of HGVS strings in parallel, yielding results
/// in input order (#975). The streaming counterpart to
/// [`parse_and_normalize_parallel`]; see [`parse_hgvs_streaming`] for the
/// rationale.
pub fn parse_and_normalize_streaming<'a, P, I, S>(
    normalizer: &'a Normalizer<P>,
    variants: I,
    num_threads: usize,
) -> impl Iterator<Item = Result<HgvsVariant, FerroError>> + 'a
where
    P: ReferenceProvider + Sync,
    I: IntoIterator<Item = S> + 'a,
    S: AsRef<str> + Send + Sync + 'a,
{
    crate::batch::streaming::StreamingMap::new(
        variants,
        num_threads,
        crate::batch::STREAM_CHUNK_ITEMS,
        move |s: &S| {
            let variant = parse_hgvs(s.as_ref())?;
            normalizer.normalize(&variant)
        },
    )
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

    // ------------------------------------------------------------------
    // Streaming API (#975)
    // ------------------------------------------------------------------

    /// A mixed corpus: valid variants across several axes, plus inputs that must
    /// fail, so the streaming path is checked on errors and not only on successes.
    fn streaming_corpus(n: usize) -> Vec<String> {
        let templates = [
            "NM_000088.3:c.459A>G",
            "NC_000001.11:g.12345A>G",
            "NP_000079.2:p.Val600Glu",
            "not a variant",
            "NM_000088.3:c.",
        ];
        (0..n)
            .map(|i| templates[i % templates.len()].to_string())
            .collect()
    }

    /// The streaming API must yield exactly what the `Vec`-based API returns, in
    /// the same order — that is the whole contract, and #975 states it as an
    /// acceptance criterion ("output is byte-identical and in the same order").
    ///
    /// Sized past one chunk (`STREAM_CHUNK_ITEMS`, which is what this path
    /// chunks at) so chunk refills, not just the first chunk, are exercised.
    /// Compared as rendered strings because
    /// `FerroError` is not `PartialEq`, and rendering is what a caller sees.
    #[test]
    fn streaming_parse_matches_the_vec_api_exactly() {
        let inputs = streaming_corpus(crate::batch::STREAM_CHUNK_ITEMS * 2 + 37);
        let render = |r: &Result<HgvsVariant, FerroError>| match r {
            Ok(v) => format!("ok:{v}"),
            Err(e) => format!("err:{e}"),
        };
        let expected: Vec<String> = parse_hgvs_parallel(&inputs).iter().map(render).collect();
        for threads in [1usize, 0, 2] {
            let streamed: Vec<String> = parse_hgvs_streaming(inputs.iter(), threads)
                .map(|r| render(&r))
                .collect();
            assert_eq!(
                streamed, expected,
                "streaming with num_threads={threads} diverged from the Vec API"
            );
        }
    }

    /// The same equality for the normalize path, which is the one the Python
    /// bindings actually drive.
    #[test]
    fn streaming_parse_and_normalize_matches_the_vec_api_exactly() {
        use crate::reference::JsonProvider;
        let normalizer = Normalizer::new(JsonProvider::with_test_data());
        let inputs = streaming_corpus(crate::batch::STREAM_CHUNK_ITEMS + 5);
        let render = |r: &Result<HgvsVariant, FerroError>| match r {
            Ok(v) => format!("ok:{v}"),
            Err(e) => format!("err:{e}"),
        };
        let expected: Vec<String> = parse_and_normalize_parallel(&normalizer, &inputs)
            .iter()
            .map(render)
            .collect();
        for threads in [1usize, 0, 2] {
            let streamed: Vec<String> =
                parse_and_normalize_streaming(&normalizer, inputs.iter(), threads)
                    .map(|r| render(&r))
                    .collect();
            assert_eq!(streamed, expected, "num_threads={threads}");
        }
    }

    /// The input iterator must be consumed *lazily*, one chunk ahead at most.
    ///
    /// This is the property that makes peak memory bounded, and it is invisible to
    /// an output-equality test: a streaming API that internally collected its
    /// input first would pass every assertion above while bounding nothing. The
    /// counter makes the laziness itself the assertion.
    #[test]
    fn the_input_stream_is_consumed_lazily() {
        use std::cell::Cell;
        let pulled = Cell::new(0usize);
        let total = crate::batch::STREAM_CHUNK_ITEMS * 3;
        let corpus = streaming_corpus(total);
        let lazy = corpus.iter().inspect(|_| pulled.set(pulled.get() + 1));
        let mut stream = parse_hgvs_streaming(lazy, 1);

        // Nothing is read until the first `next()`.
        assert_eq!(pulled.get(), 0, "constructing the stream must read nothing");

        let _ = stream.next().expect("first result");
        assert_eq!(
            pulled.get(),
            crate::batch::STREAM_CHUNK_ITEMS,
            "the first result must cost exactly one chunk of input, not the whole stream"
        );

        // Draining the rest of the first chunk reads nothing further.
        for _ in 1..crate::batch::STREAM_CHUNK_ITEMS {
            let _ = stream.next().expect("first chunk");
        }
        assert_eq!(pulled.get(), crate::batch::STREAM_CHUNK_ITEMS);

        // Crossing into the second chunk reads exactly one more chunk.
        let _ = stream.next().expect("second chunk");
        assert_eq!(pulled.get(), crate::batch::STREAM_CHUNK_ITEMS * 2);

        let remaining = stream.count();
        assert_eq!(remaining, total - crate::batch::STREAM_CHUNK_ITEMS - 1);
        assert_eq!(
            pulled.get(),
            total,
            "every input must eventually be read once"
        );
    }

    /// A short stream (fewer items than one chunk) and an empty one are both
    /// handled — the boundary the chunk loop is most likely to get wrong.
    #[test]
    fn short_and_empty_streams_are_handled() {
        let empty: Vec<String> = Vec::new();
        assert_eq!(parse_hgvs_streaming(empty.iter(), 0).count(), 0);
        let three = streaming_corpus(3);
        assert_eq!(parse_hgvs_streaming(three.iter(), 0).count(), 3);
        // Exactly one chunk: the refill after it must terminate rather than spin.
        let exact = streaming_corpus(crate::batch::STREAM_CHUNK_ITEMS);
        assert_eq!(
            parse_hgvs_streaming(exact.iter(), 2).count(),
            crate::batch::STREAM_CHUNK_ITEMS
        );
    }
}
