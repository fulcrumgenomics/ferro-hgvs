//! Bulk mutalyzer benchmark for processing large pattern sets.
//!
//! This module provides tools for running mutalyzer normalization on
//! large pattern sets (e.g., all of ClinVar) with:
//! - Streaming JSONL.gz output
//! - Incremental processing (skip successful, retry failed)
//! - Progress reporting
//! - Batch processing for memory efficiency

use crate::benchmark::types::{
    BenchmarkMetadata, BenchmarkSummary, ExistingResults, ExistingStats, MutalyzerBenchmarkConfig,
    ParseResult,
};
use crate::error::FerroError;
use chrono::Utc;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

/// Load existing results from a gzipped JSONL file.
///
/// Parses each line as a `ParseResult` and categorizes patterns as
/// successful (to skip) or failed (to retry).
pub fn load_existing_results<P: AsRef<Path>>(path: P) -> Result<ExistingResults, FerroError> {
    let path = path.as_ref();
    eprintln!("Loading existing results from: {}", path.display());

    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open existing results file: {}", e),
    })?;

    let reader: Box<dyn BufRead> = if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "gzip")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut results = ExistingResults::new();
    let mut line_count = 0;
    let mut parse_errors = 0;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if line.trim().is_empty() {
            continue;
        }

        line_count += 1;

        match serde_json::from_str::<ParseResult>(&line) {
            Ok(result) => {
                let input = result.input.clone();
                if result.success {
                    results.successful.insert(input, result);
                } else {
                    results.failed.insert(input, result);
                }
            }
            Err(e) => {
                parse_errors += 1;
                if parse_errors <= 5 {
                    eprintln!("Warning: Failed to parse line {}: {}", line_count, e);
                }
            }
        }

        if line_count % 1_000_000 == 0 {
            eprintln!(
                "  Loaded {} results ({} successful, {} failed)",
                line_count,
                results.successful_count(),
                results.failed_count()
            );
        }
    }

    eprintln!(
        "Loaded {} total results: {} successful (skip), {} failed (retry)",
        line_count,
        results.successful_count(),
        results.failed_count()
    );

    if parse_errors > 0 {
        eprintln!("Warning: {} lines failed to parse", parse_errors);
    }

    Ok(results)
}

/// JSONL writer - uses gzip if path ends in .gz, otherwise plain text.
/// Plain text is safer for Ctrl-C resume since it doesn't need finalization.
enum JsonlWriter {
    Compressed {
        writer: GzEncoder<BufWriter<File>>,
        count: usize,
    },
    Plain {
        writer: BufWriter<File>,
        count: usize,
    },
}

impl JsonlWriter {
    fn new<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let path = path.as_ref();
        let file = File::create(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create output file: {}", e),
        })?;
        let buf_writer = BufWriter::new(file);

        // Use gzip if path ends in .gz, otherwise plain text
        let is_gzip = path
            .extension()
            .is_some_and(|ext| ext == "gz" || ext == "gzip");

        if is_gzip {
            Ok(Self::Compressed {
                writer: GzEncoder::new(buf_writer, Compression::default()),
                count: 0,
            })
        } else {
            Ok(Self::Plain {
                writer: buf_writer,
                count: 0,
            })
        }
    }

    fn write(&mut self, result: &ParseResult) -> Result<(), FerroError> {
        let json = serde_json::to_string(result).map_err(|e| FerroError::Json {
            msg: format!("Failed to serialize result: {}", e),
        })?;
        match self {
            Self::Compressed { writer, count } => {
                writeln!(writer, "{}", json).map_err(|e| FerroError::Io {
                    msg: format!("Failed to write result: {}", e),
                })?;
                *count += 1;
            }
            Self::Plain { writer, count } => {
                writeln!(writer, "{}", json).map_err(|e| FerroError::Io {
                    msg: format!("Failed to write result: {}", e),
                })?;
                *count += 1;
            }
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), FerroError> {
        match self {
            Self::Compressed { writer, .. } => {
                writer.flush().map_err(|e| FerroError::Io {
                    msg: format!("Failed to flush output: {}", e),
                })?;
            }
            Self::Plain { writer, .. } => {
                writer.flush().map_err(|e| FerroError::Io {
                    msg: format!("Failed to flush output: {}", e),
                })?;
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<usize, FerroError> {
        match self {
            Self::Compressed { writer, count } => {
                writer.finish().map_err(|e| FerroError::Io {
                    msg: format!("Failed to finish gzip stream: {}", e),
                })?;
                Ok(count)
            }
            Self::Plain { mut writer, count } => {
                writer.flush().map_err(|e| FerroError::Io {
                    msg: format!("Failed to flush output: {}", e),
                })?;
                Ok(count)
            }
        }
    }
}

/// Run the mutalyzer benchmark.
///
/// Processes patterns in batches, writes results to gzipped JSONL,
/// and reports progress periodically.
pub fn run_benchmark(config: &MutalyzerBenchmarkConfig) -> Result<BenchmarkSummary, FerroError> {
    let start_time = Utc::now();
    let start_instant = Instant::now();

    // Load input patterns
    eprintln!("Loading patterns from: {}", config.input_path.display());
    let all_patterns = load_patterns(&config.input_path)?;
    eprintln!("Loaded {} patterns", all_patterns.len());

    // Load existing results if specified
    let existing = if let Some(existing_path) = &config.existing_path {
        Some(load_existing_results(existing_path)?)
    } else {
        None
    };

    // Filter patterns to process
    let patterns_to_process: Vec<&String> = if let Some(ref existing) = existing {
        all_patterns
            .iter()
            .filter(|p| {
                // Always skip successful patterns
                if existing.successful.contains_key(*p) {
                    return false;
                }
                // Skip failed patterns if skip_failed is set
                if config.skip_failed && existing.failed.contains_key(*p) {
                    return false;
                }
                true
            })
            .collect()
    } else {
        all_patterns.iter().collect()
    };

    let total_to_process = patterns_to_process.len();
    eprintln!(
        "Will process {} patterns{}",
        total_to_process,
        if let Some(ref e) = existing {
            if config.skip_failed {
                format!(
                    " (skipping {} successful, skipping {} failed)",
                    e.successful_count(),
                    e.failed_count()
                )
            } else {
                format!(
                    " (skipping {} successful, retrying {} failed)",
                    e.successful_count(),
                    e.failed_count()
                )
            }
        } else {
            String::new()
        }
    );

    // Open output file
    let mut writer = JsonlWriter::new(&config.output_path)?;

    // Copy successful results from existing file
    let mut copied_successful = 0;
    let mut copied_failed = 0;
    if let Some(ref existing) = existing {
        eprintln!(
            "Copying {} successful results from existing file...",
            existing.successful_count()
        );
        for result in existing.successful.values() {
            writer.write(result)?;
            copied_successful += 1;
        }

        // Also copy failed results if skip_failed is set
        if config.skip_failed {
            eprintln!(
                "Copying {} failed results from existing file (skip_failed=true)...",
                existing.failed_count()
            );
            for result in existing.failed.values() {
                writer.write(result)?;
                copied_failed += 1;
            }
        }

        writer.flush()?;
        eprintln!(
            "Copied {} successful, {} failed results",
            copied_successful, copied_failed
        );
    }

    // Process in batches
    let mut processed = 0;
    let mut successful = 0;
    let mut failed = 0;
    let mut error_counts: HashMap<String, usize> = HashMap::new();
    let mut last_progress = Instant::now();
    let mut retry_successes = 0;

    // Convert settings file to absolute path for Python subprocesses
    let settings_file = config.settings_file.as_ref().map(|p| {
        if p.is_absolute() {
            p.to_string_lossy().to_string()
        } else {
            std::env::current_dir()
                .map(|cwd| cwd.join(p).to_string_lossy().to_string())
                .unwrap_or_else(|_| p.to_string_lossy().to_string())
        }
    });

    for batch in patterns_to_process.chunks(config.batch_size) {
        let batch_patterns: Vec<String> = batch.iter().map(|s| (*s).clone()).collect();
        let batch_size = batch_patterns.len();

        // Run mutalyzer on this batch
        let (results, _elapsed, batch_error_counts) = run_mutalyzer_normalize_batch(
            &batch_patterns,
            config.workers,
            settings_file.as_deref(),
            config.allow_network,
        )?;

        // Write results and count successes/failures
        for result in &results {
            writer.write(result)?;
            if result.success {
                successful += 1;
                // Track retry successes
                if let Some(ref existing) = existing {
                    if existing.failed.contains_key(&result.input) {
                        retry_successes += 1;
                    }
                }
            } else {
                failed += 1;
            }
        }

        // Aggregate error counts
        for (category, count) in batch_error_counts {
            *error_counts.entry(category).or_insert(0) += count;
        }

        processed += batch_size;

        // Flush after each batch for Ctrl-C resume support
        writer.flush()?;

        // Progress reporting
        if last_progress.elapsed().as_secs() >= config.progress_interval {
            let elapsed = start_instant.elapsed().as_secs_f64();
            let throughput = processed as f64 / elapsed;
            let eta_seconds = if throughput > 0.0 {
                ((total_to_process - processed) as f64 / throughput) as u64
            } else {
                0
            };

            eprintln!(
                "[{}] Batch {}/{} ({:.1}%)",
                Utc::now().format("%Y-%m-%d %H:%M:%S"),
                processed / config.batch_size,
                total_to_process.div_ceil(config.batch_size),
                100.0 * processed as f64 / total_to_process as f64
            );
            eprintln!(
                "  Processed: {} / {}",
                format_number(processed),
                format_number(total_to_process)
            );
            eprintln!(
                "  Successful: {} ({:.1}%)",
                format_number(successful),
                100.0 * successful as f64 / processed as f64
            );
            eprintln!(
                "  Failed: {} ({:.1}%)",
                format_number(failed),
                100.0 * failed as f64 / processed as f64
            );
            eprintln!("  Throughput: {:.1} p/s", throughput);
            eprintln!("  Elapsed: {}", format_duration(elapsed as u64));
            eprintln!("  ETA: {}", format_duration(eta_seconds));

            last_progress = Instant::now();
            writer.flush()?;
        }
    }

    // Finalize output
    let written = writer.finish()?;
    let elapsed = start_instant.elapsed();

    // Build summary
    let summary = BenchmarkSummary {
        metadata: BenchmarkMetadata {
            start_time,
            end_time: Some(Utc::now()),
            total_patterns: all_patterns.len(),
            workers: config.workers,
            batch_size: config.batch_size,
            allow_network: config.allow_network,
            existing_file: config
                .existing_path
                .as_ref()
                .map(|p| p.to_string_lossy().to_string()),
        },
        processed,
        successful: successful + copied_successful,
        failed,
        elapsed_seconds: elapsed.as_secs_f64(),
        throughput: processed as f64 / elapsed.as_secs_f64(),
        error_counts,
        existing_stats: existing.as_ref().map(|e| ExistingStats {
            skipped: e.successful_count(),
            retried: e.failed_count(),
            retry_successes,
        }),
    };

    eprintln!("\n=== Benchmark Complete ===");
    eprintln!("Total patterns: {}", format_number(all_patterns.len()));
    eprintln!("Processed: {}", format_number(processed));
    eprintln!(
        "Successful: {} ({:.1}%)",
        format_number(summary.successful),
        100.0 * summary.successful as f64 / all_patterns.len() as f64
    );
    eprintln!(
        "Failed: {} ({:.1}%)",
        format_number(summary.failed),
        100.0 * summary.failed as f64 / all_patterns.len() as f64
    );
    eprintln!("Elapsed: {}", format_duration(elapsed.as_secs()));
    eprintln!("Throughput: {:.1} patterns/sec", summary.throughput);
    eprintln!(
        "Output: {} ({} results)",
        config.output_path.display(),
        written
    );

    Ok(summary)
}

/// Load patterns from a file (one per line).
fn load_patterns<P: AsRef<Path>>(path: P) -> Result<Vec<String>, FerroError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open patterns file: {}", e),
    })?;

    let reader: Box<dyn BufRead> = if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "gzip")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let patterns: Vec<String> = reader
        .lines()
        .filter_map(|line| {
            line.ok().and_then(|l| {
                let trimmed = l.trim().to_string();
                if trimmed.is_empty() {
                    None
                } else {
                    Some(trimmed)
                }
            })
        })
        .collect();

    Ok(patterns)
}

/// Run mutalyzer normalization on a batch of patterns.
///
/// This is a simplified version that calls the existing sharded subprocess approach.
#[allow(clippy::type_complexity)]
fn run_mutalyzer_normalize_batch(
    patterns: &[String],
    workers: usize,
    settings_file: Option<&str>,
    allow_network: bool,
) -> Result<
    (
        Vec<ParseResult>,
        std::time::Duration,
        HashMap<String, usize>,
    ),
    FerroError,
> {
    use std::fs::File;
    use std::io::BufRead;
    use std::io::BufReader as IoBufReader;
    use std::io::Write as IoWrite;

    // Python script for mutalyzer normalization (same as in compare.rs)
    const MUTALYZER_NORMALIZE_SCRIPT: &str = r#"
import sys
import json
import time
import os
import warnings

# Suppress urllib3 SSL warnings
warnings.filterwarnings('ignore', category=Warning, module='urllib3')

# CRITICAL: Set MUTALYZER_SETTINGS env var BEFORE importing mutalyzer
# mutalyzer-retriever reads config from a file specified by MUTALYZER_SETTINGS
if len(sys.argv) > 3 and sys.argv[3]:
    settings_file = sys.argv[3]
    if os.path.exists(settings_file):
        os.environ['MUTALYZER_SETTINGS'] = settings_file

from mutalyzer_hgvs_parser import to_model
try:
    from mutalyzer.normalizer import normalize as mutalyzer_normalize
    from mutalyzer.description import Description
    HAS_NORMALIZER = True
except ImportError:
    HAS_NORMALIZER = False

# Block network access if requested (after imports to avoid breaking ssl module)
NETWORK_BLOCKED = len(sys.argv) > 4 and sys.argv[4].lower() != 'true'
if NETWORK_BLOCKED:
    import socket
    _original_socket = socket.socket
    def _blocked_socket(*args, **kwargs):
        raise OSError("Network access blocked for fair benchmarking")
    socket.socket = _blocked_socket

def normalize_variant(hgvs_string):
    try:
        if not HAS_NORMALIZER:
            model = to_model(hgvs_string)
            if model is None:
                return {"success": False, "error": "Parse failed", "error_category": "PARSE_ERROR"}
            from mutalyzer_hgvs_parser import to_hgvs
            return {"success": True, "output": to_hgvs(model)}

        # Description expects a string, not a model - it parses internally
        description = Description(hgvs_string)
        description.normalize()

        if description.errors:
            error = description.errors[0]
            code = error.get('code', 'UNKNOWN')
            msg = error.get('details', str(error))
            category = code if code else 'UNKNOWN_ERROR'
            return {"success": False, "error": msg, "error_category": category}

        result = description.output()
        if result and result.get('normalized_description'):
            return {"success": True, "output": result['normalized_description']}
        else:
            return {"success": False, "error": "No normalized output", "error_category": "NO_OUTPUT"}

    except Exception as e:
        error_str = str(e)
        if 'Network access blocked' in error_str:
            category = 'NETWORK_BLOCKED'
        elif 'LRG' in error_str:
            category = 'LRG_NETWORK_REQUIRED'
        elif 'retrieve' in error_str.lower() or 'fetch' in error_str.lower():
            category = 'CACHE_MISS'
        else:
            category = 'EXCEPTION'
        return {"success": False, "error": error_str, "error_category": category}

import gc
import resource

input_file = sys.argv[1]
output_file = sys.argv[2]

error_counts = {}
successful = 0
failed = 0
total = 0

def get_mem_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024

start = time.time()
start_mem = get_mem_mb()

# Stream results to JSONL to avoid memory accumulation
with open(output_file, 'w') as out_f:
    with open(input_file, 'r') as in_f:
        for line in in_f:
            pattern = line.strip()
            if not pattern:
                continue

            result = normalize_variant(pattern)
            result['input'] = pattern
            total += 1

            # Write immediately, don't accumulate
            out_f.write(json.dumps(result) + '\n')
            out_f.flush()  # Force write to disk

            if result['success']:
                successful += 1
            else:
                failed += 1
                cat = result.get('error_category', 'UNKNOWN')
                error_counts[cat] = error_counts.get(cat, 0) + 1

            # GC every 100 patterns, but only log every 5000
            if total % 100 == 0:
                gc.collect()
            if total % 5000 == 0:
                mem = get_mem_mb()
                print(f"[Worker] Processed {total}, Mem: {mem:.1f}MB", file=sys.stderr, flush=True)

elapsed = time.time() - start
end_mem = get_mem_mb()

# Write summary as last line (marked with _summary key)
summary = {
    '_summary': True,
    'tool': 'mutalyzer',
    'total_patterns': total,
    'successful': successful,
    'failed': failed,
    'elapsed_seconds': elapsed,
    'patterns_per_second': total / elapsed if elapsed > 0 else 0,
    'error_counts': error_counts,
    'start_mem_mb': start_mem,
    'end_mem_mb': end_mem
}
with open(output_file, 'a') as out_f:
    out_f.write(json.dumps(summary) + '\n')

# Silent completion - main process logs batch progress
"#;

    /// Summary line from mutalyzer worker subprocess (marked with _summary: true).
    #[derive(Debug, serde::Deserialize)]
    #[allow(dead_code)]
    struct MutalyzerSummary {
        #[serde(rename = "_summary")]
        is_summary: bool,
        #[serde(default)]
        error_counts: HashMap<String, usize>,
    }

    let temp_dir = tempfile::tempdir().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp directory: {}", e),
    })?;

    // Shard patterns
    let shards = shard_patterns(patterns, workers);

    // Write shard input files
    let mut shard_paths = Vec::new();
    for (i, shard) in shards.iter().enumerate() {
        let input_path = temp_dir.path().join(format!("shard_{}_input.txt", i));
        let output_path = temp_dir.path().join(format!("shard_{}_output.json", i));

        let mut file = File::create(&input_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create shard file: {}", e),
        })?;
        for pattern in shard {
            writeln!(file, "{}", pattern).map_err(|e| FerroError::Io {
                msg: format!("Failed to write pattern: {}", e),
            })?;
        }

        shard_paths.push((input_path, output_path));
    }

    // Run workers in parallel
    let start = Instant::now();

    // Spawn all processes
    let mut children: Vec<(usize, std::process::Child)> = Vec::new();
    for (i, (input, output)) in shard_paths.iter().enumerate() {
        let mut cmd = std::process::Command::new("python3");
        cmd.args(["-c", MUTALYZER_NORMALIZE_SCRIPT]);
        cmd.arg(input.display().to_string());
        cmd.arg(output.display().to_string());
        if let Some(settings) = settings_file {
            cmd.arg(settings);
        } else {
            cmd.arg("");
        }
        cmd.arg(if allow_network { "true" } else { "false" });

        let child = cmd.spawn().map_err(|e| FerroError::Io {
            msg: format!("Failed to spawn Python: {}", e),
        })?;
        children.push((i, child));
    }

    // Wait for all processes
    for (_i, mut child) in children {
        let status = child.wait().map_err(|e| FerroError::Io {
            msg: format!("Failed to wait for Python: {}", e),
        })?;
        if !status.success() {
            return Err(FerroError::Io {
                msg: format!(
                    "Mutalyzer normalizer failed with exit code: {:?}",
                    status.code()
                ),
            });
        }
    }

    let elapsed = start.elapsed();

    // Collect results from JSONL files
    let mut all_results = Vec::new();
    let mut aggregated_error_counts: HashMap<String, usize> = HashMap::new();

    for (_, output_path) in &shard_paths {
        let file = File::open(output_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open shard output: {}", e),
        })?;
        let reader = IoBufReader::new(file);

        // Parse JSONL: each line is either a result or the summary (last line)
        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read line: {}", e),
            })?;

            if line.trim().is_empty() {
                continue;
            }

            // Check if this is the summary line
            if line.contains("\"_summary\"") {
                let summary: MutalyzerSummary =
                    serde_json::from_str(&line).map_err(|e| FerroError::Json {
                        msg: format!("Failed to parse summary: {}", e),
                    })?;
                for (category, count) in summary.error_counts {
                    *aggregated_error_counts.entry(category).or_insert(0) += count;
                }
            } else {
                // Regular result line
                let result: ParseResult =
                    serde_json::from_str(&line).map_err(|e| FerroError::Json {
                        msg: format!("Failed to parse result: {}", e),
                    })?;
                all_results.push(result);
            }
        }
    }

    Ok((all_results, elapsed, aggregated_error_counts))
}

/// Shard patterns across workers.
fn shard_patterns(patterns: &[String], workers: usize) -> Vec<Vec<String>> {
    let mut shards: Vec<Vec<String>> = (0..workers).map(|_| Vec::new()).collect();

    for (i, pattern) in patterns.iter().enumerate() {
        shards[i % workers].push(pattern.clone());
    }

    shards
}

/// Format a number with thousands separators.
fn format_number(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

/// Format a duration in human-readable form.
fn format_duration(seconds: u64) -> String {
    if seconds < 60 {
        format!("{}s", seconds)
    } else if seconds < 3600 {
        format!("{}m {}s", seconds / 60, seconds % 60)
    } else {
        format!("{}h {}m", seconds / 3600, (seconds % 3600) / 60)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_existing_results() {
        let mut results = ExistingResults::new();
        results.successful.insert(
            "NM_000001.1:c.1A>G".to_string(),
            ParseResult {
                input: "NM_000001.1:c.1A>G".to_string(),
                success: true,
                output: Some("NM_000001.1:c.1A>G".to_string()),
                error: None,
                error_category: None,
                ref_mismatch: None,
                details: None,
            },
        );
        results.failed.insert(
            "NM_000002.1:c.2A>G".to_string(),
            ParseResult {
                input: "NM_000002.1:c.2A>G".to_string(),
                success: false,
                output: None,
                error: Some("Test error".to_string()),
                error_category: Some("TEST".to_string()),
                ref_mismatch: None,
                details: None,
            },
        );

        assert!(!results.should_process("NM_000001.1:c.1A>G")); // Skip successful
        assert!(results.should_process("NM_000002.1:c.2A>G")); // Retry failed
        assert!(results.should_process("NM_000003.1:c.3A>G")); // Process new
    }

    #[test]
    fn test_format_number() {
        assert_eq!(format_number(0), "0");
        assert_eq!(format_number(999), "999");
        assert_eq!(format_number(1000), "1,000");
        assert_eq!(format_number(1234567), "1,234,567");
    }

    #[test]
    fn test_format_duration() {
        assert_eq!(format_duration(0), "0s");
        assert_eq!(format_duration(59), "59s");
        assert_eq!(format_duration(60), "1m 0s");
        assert_eq!(format_duration(3661), "1h 1m");
    }

    #[test]
    fn test_shard_patterns() {
        let patterns: Vec<String> = (0..10).map(|i| format!("pattern_{}", i)).collect();
        let shards = shard_patterns(&patterns, 3);

        assert_eq!(shards.len(), 3);
        assert_eq!(shards[0].len(), 4); // 0, 3, 6, 9
        assert_eq!(shards[1].len(), 3); // 1, 4, 7
        assert_eq!(shards[2].len(), 3); // 2, 5, 8
    }
}
