//! Normalization benchmarks.

use crate::benchmark::mutalyzer::MutalyzerClient;
use crate::benchmark::types::{ParseResult, ShardResults, TimingInfo};
use crate::commands;
use crate::FerroError;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
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
