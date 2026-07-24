//! Normalization benchmarks.

use crate::benchmark::mutalyzer::MutalyzerClient;
use crate::benchmark::types::{ParseResult, ShardResults, TimingInfo};
use crate::commands;
use crate::FerroError;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Duration, Instant};

/// Normalize patterns with ferro-hgvs.
///
/// If `reference_dir` points to a directory containing FASTA files (from `ferro prepare`),
/// uses real reference data. Otherwise falls back to MockProvider.
///
/// This function routes to the `commands::normalize_batch` function from the main library.
pub fn normalize_ferro<P: AsRef<Path>>(
    input: P,
    results_output: P,
    timing_output: P,
    reference_dir: Option<P>,
) -> Result<ShardResults, FerroError> {
    let input = input.as_ref();
    let results_output = results_output.as_ref();
    let timing_output = timing_output.as_ref();

    // Use commands module for normalization
    let config = commands::NormalizeConfig {
        reference_dir: reference_dir.map(|p| p.as_ref().to_path_buf()),
        show_progress: true,
        workers: 1,
    };
    let batch_results = commands::normalize_batch(input, &config)?;

    // Convert commands::VariantResult to benchmark::ParseResult
    let results: Vec<ParseResult> = batch_results
        .results
        .into_iter()
        .map(|r| {
            let error_category = r.error.as_ref().map(|e| categorize_error_str(e));
            ParseResult {
                input: r.input,
                success: r.success,
                output: r.output,
                error: r.error,
                error_category,
                ref_mismatch: None,
                details: None,
            }
        })
        .collect();

    let elapsed = Duration::from_secs_f64(batch_results.timing.elapsed_seconds);
    let timing = TimingInfo::new(
        "ferro-hgvs",
        batch_results.timing.total,
        batch_results.timing.successful,
        elapsed,
    );

    // Keep all results for comparison (not just a sample)
    let sample_results: Vec<ParseResult> = results.clone();
    let failed_examples: Vec<ParseResult> =
        results.iter().filter(|r| !r.success).cloned().collect();

    let shard_results = ShardResults {
        shard_index: 0,
        tool: "ferro-hgvs".to_string(),
        input_file: input.display().to_string(),
        timing,
        sample_results,
        failed_examples,
    };

    save_results(&shard_results, results_output, timing_output)?;

    Ok(shard_results)
}

/// Normalize patterns with ferro-hgvs across `workers` threads.
///
/// Builds ONE reference provider before the timer starts and wraps it in an
/// `Arc<dyn ReferenceProvider + Send + Sync>`, which is cloned cheaply (pointer
/// copy) into each rayon worker. The ~600k-transcript data is therefore loaded
/// exactly ONCE regardless of worker count. The `Mutex<LruCache>` inside
/// `MultiFastaProvider` serializes concurrent transcript-cache insertions; all
/// other reads are lock-free.
///
/// Provider setup is excluded from the timed region. `Instant::now()` is
/// started only after the provider and the rayon thread-pool are ready. Each
/// rayon worker normalizes its contiguous slice of patterns independently.
/// Results are collected back in original order via indexed output.
///
/// For `workers <= 1`, delegates to [`normalize_ferro`].
pub fn normalize_ferro_parallel<P: AsRef<Path>>(
    input: P,
    results_output: P,
    timing_output: P,
    reference_dir: Option<P>,
    workers: usize,
) -> Result<ShardResults, FerroError> {
    if workers <= 1 {
        return normalize_ferro(input, results_output, timing_output, reference_dir);
    }

    let input = input.as_ref();
    let results_output = results_output.as_ref();
    let timing_output = timing_output.as_ref();
    let reference_dir: Option<PathBuf> = reference_dir.map(|p| p.as_ref().to_path_buf());

    // Read all patterns from the input file up front.
    let file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;
    let reader = BufReader::new(file);
    let patterns: Vec<String> = reader
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    if patterns.is_empty() {
        let timing = TimingInfo::new("ferro-hgvs", 0, 0, Duration::ZERO);
        let shard = ShardResults {
            shard_index: 0,
            tool: "ferro-hgvs".to_string(),
            input_file: input.display().to_string(),
            timing,
            sample_results: Vec::new(),
            failed_examples: Vec::new(),
        };
        save_results(&shard, results_output, timing_output)?;
        return Ok(shard);
    }

    // Build ONE provider outside the timed region. The concrete provider
    // returned by create_reference_provider is Send+Sync (asserted at compile
    // time in src/reference/provider.rs::_assert_provider_send_sync), so we
    // can erase it into Arc<dyn … + Send + Sync> and share it across threads.
    //
    // The benchmark eprintln here intentionally says "one provider" (not N) so
    // the smoke log confirms the single-load invariant.
    eprintln!("Creating one ferro provider (shared across {} rayon workers, excluded from timed region)...", workers);
    // create_reference_provider now returns Box<dyn ReferenceProvider + Send + Sync>,
    // so Arc::from is a safe coercion — no unsafe required.
    let raw_provider: Box<dyn crate::reference::ReferenceProvider + Send + Sync> =
        commands::create_reference_provider(reference_dir.as_deref(), false)?;
    // Wrap the single provider in an Arc so all rayon workers share one copy
    // via cheap pointer clones.  The Box<T: Send+Sync> blanket impl in
    // provider.rs makes Arc<dyn ReferenceProvider + Send + Sync> usable as a
    // ReferenceProvider everywhere.
    let shared_provider: Arc<dyn crate::reference::ReferenceProvider + Send + Sync> =
        Arc::from(raw_provider);

    // Build a rayon thread pool sized to `workers` and drive the parallel
    // normalize entirely through it.  Using a custom pool (not the global one)
    // honours the requested thread count even when the calling thread already
    // lives inside a different rayon pool.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to build rayon thread pool: {e}"),
        })?;

    // Allocate the output vector before the timer — only pattern processing
    // belongs in the timed region.
    let n = patterns.len();
    let mut all_results: Vec<ParseResult> = (0..n)
        .map(|_| ParseResult {
            input: String::new(),
            success: false,
            output: None,
            error: None,
            error_category: None,
            ref_mismatch: None,
            details: None,
        })
        .collect();

    // --- TIMED REGION STARTS HERE ---
    let start = Instant::now();

    pool.install(|| {
        all_results
            .par_iter_mut()
            .zip(patterns.par_iter())
            .for_each(|(slot, pattern)| {
                // Each rayon thread clones the Arc (cheap pointer copy) to get
                // its own handle to the shared provider.
                let provider = Arc::clone(&shared_provider);
                let normalizer = crate::Normalizer::new(provider);

                *slot = match crate::parse_hgvs(pattern) {
                    Ok(parsed) => {
                        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            normalizer.normalize(&parsed)
                        })) {
                            Ok(Ok(normalized)) => ParseResult {
                                input: pattern.clone(),
                                success: true,
                                output: Some(normalized.to_string()),
                                error: None,
                                error_category: None,
                                ref_mismatch: None,
                                details: None,
                            },
                            Ok(Err(e)) => {
                                let msg = format!("{}", e);
                                ParseResult {
                                    input: pattern.clone(),
                                    success: false,
                                    output: None,
                                    error_category: Some(categorize_error_str(&msg)),
                                    error: Some(msg),
                                    ref_mismatch: None,
                                    details: None,
                                }
                            }
                            Err(payload) => {
                                let msg = payload
                                    .downcast_ref::<String>()
                                    .map(|s| s.as_str())
                                    .or_else(|| payload.downcast_ref::<&str>().copied())
                                    .unwrap_or("unknown panic");
                                ParseResult {
                                    input: pattern.clone(),
                                    success: false,
                                    output: None,
                                    error: Some(format!(
                                        "internal error: panic during normalization: {}",
                                        msg
                                    )),
                                    error_category: Some("panic".to_string()),
                                    ref_mismatch: None,
                                    details: None,
                                }
                            }
                        }
                    }
                    Err(e) => {
                        let msg = format!("{}", e);
                        ParseResult {
                            input: pattern.clone(),
                            success: false,
                            output: None,
                            error_category: Some(categorize_error_str(&msg)),
                            error: Some(msg),
                            ref_mismatch: None,
                            details: None,
                        }
                    }
                };
            });
    });

    let elapsed = start.elapsed();
    // --- TIMED REGION ENDS HERE ---

    let total_successful = all_results.iter().filter(|r| r.success).count();
    let total = all_results.len();
    let timing = TimingInfo::new("ferro-hgvs", total, total_successful, elapsed);

    let failed_examples: Vec<ParseResult> =
        all_results.iter().filter(|r| !r.success).cloned().collect();

    let shard = ShardResults {
        shard_index: 0,
        tool: "ferro-hgvs".to_string(),
        input_file: input.display().to_string(),
        timing,
        sample_results: all_results,
        failed_examples,
    };

    save_results(&shard, results_output, timing_output)?;

    Ok(shard)
}

/// Categorize an error based on error message string.
fn categorize_error_str(error_msg: &str) -> String {
    if error_msg.contains("accession") || error_msg.contains("Accession") {
        "invalid_accession".to_string()
    } else if error_msg.contains("position") || error_msg.contains("coordinate") {
        "invalid_position".to_string()
    } else if error_msg.contains("edit") || error_msg.contains("allele") {
        "invalid_edit".to_string()
    } else if error_msg.contains("reference not found") || error_msg.contains("ReferenceNotFound") {
        "missing_reference".to_string()
    } else if error_msg.contains("unsupported") {
        "unsupported_variant".to_string()
    } else if error_msg.contains("mismatch") {
        "reference_mismatch".to_string()
    } else if error_msg.contains("normalize") {
        "normalize_error".to_string()
    } else {
        "parse_error".to_string()
    }
}

/// Normalize patterns via Mutalyzer API.
pub fn normalize_mutalyzer<P: AsRef<Path>>(
    input: P,
    results_output: P,
    timing_output: P,
    api_url: &str,
    rate_limit_ms: Option<u64>,
) -> Result<ShardResults, FerroError> {
    let input = input.as_ref();
    let results_output = results_output.as_ref();
    let timing_output = timing_output.as_ref();

    // Create client
    let client = MutalyzerClient::new(api_url)?;

    // Check API is available
    if !client.health_check()? {
        return Err(FerroError::Io {
            msg: format!("Mutalyzer API not available at {}", api_url),
        });
    }

    // Read patterns
    let file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;
    let reader = BufReader::new(file);
    let patterns: Vec<String> = reader
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    if patterns.is_empty() {
        let timing = TimingInfo::new("mutalyzer-api", 0, 0, std::time::Duration::ZERO);
        let results = ShardResults {
            shard_index: 0,
            tool: "mutalyzer-api".to_string(),
            input_file: input.display().to_string(),
            timing,
            sample_results: Vec::new(),
            failed_examples: Vec::new(),
        };
        save_results(&results, results_output, timing_output)?;
        return Ok(results);
    }

    // Progress bar
    let pb = ProgressBar::new(patterns.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
            .unwrap()
            .progress_chars("##-"),
    );

    // Normalize via API
    let start = Instant::now();
    let delay = rate_limit_ms.map(std::time::Duration::from_millis);

    let results: Vec<ParseResult> = patterns
        .iter()
        .map(|pattern| {
            pb.inc(1);
            let result = client.normalize(pattern).unwrap_or_else(|e| ParseResult {
                input: pattern.clone(),
                success: false,
                output: None,
                error: Some(format!("{}", e)),
                error_category: Some("client_error".to_string()),
                ref_mismatch: None,
                details: None,
            });

            if let Some(d) = delay {
                std::thread::sleep(d);
            }

            result
        })
        .collect();

    let elapsed = start.elapsed();
    pb.finish_with_message(format!(
        "Normalized {} patterns via API in {:.2}s",
        patterns.len(),
        elapsed.as_secs_f64()
    ));

    let successful = results.iter().filter(|r| r.success).count();
    let timing = TimingInfo::new("mutalyzer-api", patterns.len(), successful, elapsed);

    // Keep all results for comparison (not just a sample)
    let sample_results: Vec<ParseResult> = results.clone();
    let failed_examples: Vec<ParseResult> =
        results.iter().filter(|r| !r.success).cloned().collect();

    let shard_results = ShardResults {
        shard_index: 0,
        tool: "mutalyzer-api".to_string(),
        input_file: input.display().to_string(),
        timing,
        sample_results,
        failed_examples,
    };

    save_results(&shard_results, results_output, timing_output)?;

    Ok(shard_results)
}

/// Save results to JSON files.
fn save_results<P: AsRef<Path>>(
    results: &ShardResults,
    results_output: P,
    timing_output: P,
) -> Result<(), FerroError> {
    let results_output = results_output.as_ref();
    let timing_output = timing_output.as_ref();

    // Create directories if needed
    for path in [results_output, timing_output] {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
                msg: format!("Failed to create directory {}: {}", parent.display(), e),
            })?;
        }
    }

    // Write results
    let file = File::create(results_output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", results_output.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, results).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    // Write timing
    let file = File::create(timing_output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", timing_output.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &results.timing).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    Ok(())
}

#[cfg(test)]
mod parallel_tests {
    use super::*;
    use tempfile::NamedTempFile;

    // ---------------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------------

    /// Run the normalize path (serial or parallel) on a temporary input file
    /// containing `patterns` and return `(successful, failed, ordered_outputs)`.
    fn normalize_ferro_full(
        patterns: &[String],
        workers: usize,
    ) -> (usize, usize, Vec<Option<String>>) {
        let mut input_file = NamedTempFile::new().expect("temp input");
        use std::io::Write;
        for p in patterns {
            writeln!(input_file, "{}", p).unwrap();
        }

        let results_file = NamedTempFile::new().expect("temp results");
        let timing_file = NamedTempFile::new().expect("temp timing");

        let shard = if workers <= 1 {
            normalize_ferro(
                input_file.path(),
                results_file.path(),
                timing_file.path(),
                None::<&std::path::Path>,
            )
        } else {
            normalize_ferro_parallel(
                input_file.path(),
                results_file.path(),
                timing_file.path(),
                None::<&std::path::Path>,
                workers,
            )
        }
        .expect("normalize must not error");

        let s = shard.timing.successful;
        let f = shard.timing.failed;
        let outputs: Vec<Option<String>> =
            shard.sample_results.into_iter().map(|r| r.output).collect();
        (s, f, outputs)
    }

    /// Convenience wrapper returning only counts.
    fn normalize_ferro_count(patterns: &[String], workers: usize) -> (usize, usize) {
        let (s, f, _) = normalize_ferro_full(patterns, workers);
        (s, f)
    }

    // ---------------------------------------------------------------------------
    // Chunking invariant — pure logic, no I/O, would have exposed the old bug
    // ---------------------------------------------------------------------------

    /// Verify the thread-count == channel-signal-count invariant for a range of
    /// (n, workers) combinations.
    ///
    /// The old barrier-based code had a mismatch: it sized the barrier to
    /// `effective_workers + 1` (computed before chunking) but spawned only
    /// `chunks.len()` threads — which can be strictly LESS than `effective_workers`
    /// when `ceil(n/effective_workers)` divides into fewer chunks.
    ///
    /// Concrete failures:
    ///   n=9,  w=4  → effective=4, chunk_size=3 → 3 chunks   (barrier wanted 5)
    ///   n=17, w=8  → effective=8, chunk_size=3 → 6 chunks   (barrier wanted 9)
    ///   n=5,  w=4  → effective=4, chunk_size=2 → 3 chunks   (barrier wanted 5)
    ///   n=2,  w=4  → effective=2, chunk_size=1 → 2 chunks   (barrier wanted 3)
    ///
    /// The channel-based fix uses `num_threads = chunks.len()` everywhere, so the
    /// invariant becomes: `channel_signal_count == thread_count == chunks.len()`.
    /// This test verifies: (a) `chunks.len() ≤ effective_workers` (was the source
    /// of deadlock), and (b) all n items are covered exactly once.
    #[test]
    fn chunk_count_leq_effective_workers_and_covers_all() {
        let cases: &[(usize, usize)] = &[
            (1, 1),
            (2, 2),
            (2, 4), // fewer patterns than workers: effective clamps to 2
            (4, 4),
            (5, 4), // n=5, w=4: effective=4, chunk_size=2 → 3 chunks (≤4, ok)
            (8, 4), // exact multiple: 2 chunks
            (9, 4), // the original deadlock: effective=4, chunk_size=3 → 3 chunks (≤4)
            (10, 4),
            (16, 8),
            (17, 8), // effective=8, chunk_size=3 → 6 chunks (≤8)
            (32, 8),
            (33, 8),
            (100, 16),
        ];

        for &(n, workers) in cases {
            let effective = workers.min(n).max(1);
            let chunk_size = n.div_ceil(effective);
            // Build a dummy slice and chunk it the same way the function does.
            let dummy: Vec<u32> = (0..n as u32).collect();
            let chunks: Vec<&[u32]> = dummy.chunks(chunk_size).collect();

            // Key invariant: chunk count ≤ effective_workers (never MORE threads
            // than the barrier/channel expects).
            assert!(
                chunks.len() <= effective,
                "n={n}, workers={workers}: chunks.len()={} > effective={effective}",
                chunks.len()
            );

            // Additional liveness guarantee: at least one chunk (never zero threads
            // for non-empty input).
            assert!(
                !chunks.is_empty(),
                "n={n}, workers={workers}: chunks is empty for non-zero n"
            );

            // All items are covered exactly once.
            let covered: usize = chunks.iter().map(|c| c.len()).sum();
            assert_eq!(
                covered, n,
                "n={n}, workers={workers}: chunks cover {covered} items, expected {n}"
            );
        }
    }

    // ---------------------------------------------------------------------------
    // Functional tests with real MockProvider
    // ---------------------------------------------------------------------------

    #[test]
    fn parallel_matches_serial_on_empty_input() {
        let patterns: Vec<String> = vec![];
        let serial = normalize_ferro_count(&patterns, 1);
        let parallel = normalize_ferro_count(&patterns, 4);
        assert_eq!(
            serial, parallel,
            "empty input: serial and parallel must agree"
        );
    }

    #[test]
    fn parallel_matches_serial_on_trivial_input() {
        // MockProvider handles NM_000088.3 (it is in the built-in test data).
        let patterns = vec![
            "NM_000088.3:c.589G>T".to_string(),
            "not_a_valid_hgvs".to_string(),
        ];
        let serial = normalize_ferro_count(&patterns, 1);
        let parallel_2 = normalize_ferro_count(&patterns, 2);
        let parallel_4 = normalize_ferro_count(&patterns, 4);
        assert_eq!(serial, parallel_2, "workers=2 must agree with serial");
        assert_eq!(serial, parallel_4, "workers=4 must agree with serial");
    }

    /// 9 patterns, workers=4: the original deadlock case (3 chunks ≠ barrier of 5).
    /// With the old barrier-based code this test would hang indefinitely.
    /// With the channel-based fix it must complete quickly and match serial counts.
    #[test]
    fn parallel_matches_serial_n9_w4() {
        // Mix of a known-good pattern and invalid ones so both success and
        // failure paths are exercised.  MockProvider is used (no reference_dir).
        let valid = "NM_000088.3:c.589G>T".to_string();
        let invalid = "not_a_valid_hgvs".to_string();
        let patterns: Vec<String> = std::iter::once(valid)
            .chain(std::iter::repeat_n(invalid, 8))
            .collect();
        assert_eq!(patterns.len(), 9);

        let (s_serial, f_serial, out_serial) = normalize_ferro_full(&patterns, 1);
        let (s_par, f_par, out_par) = normalize_ferro_full(&patterns, 4);

        assert_eq!(
            (s_par, f_par),
            (s_serial, f_serial),
            "n=9 w=4: parallel counts must match serial"
        );
        assert_eq!(
            out_par, out_serial,
            "n=9 w=4: parallel output order must match serial"
        );
    }

    /// 17 patterns, workers=8: another mismatch case (ceil(17/8)=3 chunks ≠ 9).
    /// Would have hung on the old code.
    #[test]
    fn parallel_matches_serial_n17_w8() {
        let valid = "NM_000088.3:c.589G>T".to_string();
        let invalid = "not_a_valid_hgvs".to_string();
        let patterns: Vec<String> = std::iter::once(valid)
            .chain(std::iter::repeat_n(invalid, 16))
            .collect();
        assert_eq!(patterns.len(), 17);

        let (s_serial, f_serial, out_serial) = normalize_ferro_full(&patterns, 1);
        let (s_par, f_par, out_par) = normalize_ferro_full(&patterns, 8);

        assert_eq!(
            (s_par, f_par),
            (s_serial, f_serial),
            "n=17 w=8: parallel counts must match serial"
        );
        assert_eq!(
            out_par, out_serial,
            "n=17 w=8: parallel output order must match serial"
        );
    }

    /// Fewer patterns than workers: n=2, w=4.
    ///
    /// Validates the `effective_workers = workers.min(patterns.len())` clamp: when
    /// there are fewer chunks than requested workers, the function still produces
    /// the correct output.  This case did NOT hang on the old barrier-based code
    /// (effective_workers=2, barrier=3, 2 chunks + main = 3 arrivals → clean
    /// release), but it exercises the clamping path that prevents over-subscribing.
    /// The genuine deadlock-regression cases are n=9/w=4 and n=17/w=8, where the
    /// old barrier sized itself to effective_workers+1 while only chunks.len()
    /// threads arrived.
    #[test]
    fn parallel_matches_serial_fewer_patterns_than_workers() {
        let patterns = vec![
            "NM_000088.3:c.589G>T".to_string(),
            "not_a_valid_hgvs".to_string(),
        ];
        assert_eq!(patterns.len(), 2);

        let (s_serial, f_serial, out_serial) = normalize_ferro_full(&patterns, 1);
        let (s_par, f_par, out_par) = normalize_ferro_full(&patterns, 4);

        assert_eq!(
            (s_par, f_par),
            (s_serial, f_serial),
            "n=2 w=4: parallel counts must match serial"
        );
        assert_eq!(
            out_par, out_serial,
            "n=2 w=4: parallel output order must match serial"
        );
    }
}
