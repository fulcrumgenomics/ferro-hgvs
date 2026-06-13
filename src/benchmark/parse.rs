//! Parsing benchmarks for ferro-hgvs.

use crate::benchmark::types::{ParseResult, ShardResults, TimingInfo};
use crate::commands;
use crate::parse_hgvs;
use crate::FerroError;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use std::time::{Duration, Instant};

/// Parse a shard of patterns with ferro-hgvs.
///
/// Returns timing information and sample results for analysis.
pub fn parse_ferro<P: AsRef<Path>>(
    input: P,
    results_output: P,
    timing_output: P,
    shard_index: usize,
) -> Result<ShardResults, FerroError> {
    let input = input.as_ref();
    let results_output = results_output.as_ref();
    let timing_output = timing_output.as_ref();

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
        let timing = TimingInfo::new("ferro-hgvs", 0, 0, std::time::Duration::ZERO);
        let results = ShardResults {
            shard_index,
            tool: "ferro-hgvs".to_string(),
            input_file: input.display().to_string(),
            timing,
            sample_results: Vec::new(),
            failed_examples: Vec::new(),
        };
        save_results(&results, results_output, timing_output)?;
        return Ok(results);
    }

    // Parse all patterns and time it
    let start = Instant::now();

    let results: Vec<ParseResult> = patterns
        .iter()
        .map(|pattern| match parse_hgvs(pattern) {
            Ok(parsed) => ParseResult {
                input: pattern.clone(),
                success: true,
                output: Some(parsed.to_string()),
                error: None,
                error_category: None,
                ref_mismatch: None,
                details: None,
            },
            Err(e) => ParseResult {
                input: pattern.clone(),
                success: false,
                output: None,
                error: Some(format!("{}", e)),
                error_category: Some(categorize_error(&e)),
                ref_mismatch: None,
                details: None,
            },
        })
        .collect();

    let elapsed = start.elapsed();

    let successful = results.iter().filter(|r| r.success).count();
    let timing = TimingInfo::new("ferro-hgvs", patterns.len(), successful, elapsed);

    // Collect samples
    let sample_results: Vec<ParseResult> = results.iter().take(100).cloned().collect();

    let failed_examples: Vec<ParseResult> =
        results.iter().filter(|r| !r.success).cloned().collect();

    let shard_results = ShardResults {
        shard_index,
        tool: "ferro-hgvs".to_string(),
        input_file: input.display().to_string(),
        timing,
        sample_results,
        failed_examples,
    };

    save_results(&shard_results, results_output, timing_output)?;

    Ok(shard_results)
}

/// Parse patterns with ferro and output in unified format.
///
/// Returns ToolParseOutput which is the standard format across all tools.
pub fn parse_ferro_unified<P: AsRef<Path>>(
    input: P,
    output: P,
) -> Result<super::types::ToolParseOutput, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // Use the commands module for parsing
    let batch_results = commands::parse_batch(input, false)?;

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
    let tool_output = super::types::ToolParseOutput::new("ferro", results, elapsed);

    // Write to file
    let file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &tool_output).map_err(|e| FerroError::Json {
        msg: format!("Failed to write results: {}", e),
    })?;

    Ok(tool_output)
}

/// Parse patterns in parallel using Rayon.
pub fn parse_ferro_parallel<P: AsRef<Path>>(
    input: P,
    results_output: P,
    timing_output: P,
    shard_index: usize,
) -> Result<ShardResults, FerroError> {
    let input = input.as_ref();
    let results_output = results_output.as_ref();
    let timing_output = timing_output.as_ref();

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
        let timing = TimingInfo::new("ferro-hgvs", 0, 0, std::time::Duration::ZERO);
        let results = ShardResults {
            shard_index,
            tool: "ferro-hgvs".to_string(),
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

    // Parse all patterns in parallel
    let start = Instant::now();

    let results: Vec<ParseResult> = patterns
        .par_iter()
        .map(|pattern| {
            pb.inc(1);
            match parse_hgvs(pattern) {
                Ok(parsed) => ParseResult {
                    input: pattern.clone(),
                    success: true,
                    output: Some(parsed.to_string()),
                    error: None,
                    error_category: None,
                    ref_mismatch: None,
                    details: None,
                },
                Err(e) => ParseResult {
                    input: pattern.clone(),
                    success: false,
                    output: None,
                    error: Some(format!("{}", e)),
                    error_category: Some(categorize_error(&e)),
                    ref_mismatch: None,
                    details: None,
                },
            }
        })
        .collect();

    let elapsed = start.elapsed();
    pb.finish_with_message(format!(
        "Parsed {} patterns in {:.2}s",
        patterns.len(),
        elapsed.as_secs_f64()
    ));

    let successful = results.iter().filter(|r| r.success).count();
    let timing = TimingInfo::new("ferro-hgvs", patterns.len(), successful, elapsed);

    // Collect samples
    let sample_results: Vec<ParseResult> = results.iter().take(100).cloned().collect();

    let failed_examples: Vec<ParseResult> =
        results.iter().filter(|r| !r.success).cloned().collect();

    let shard_results = ShardResults {
        shard_index,
        tool: "ferro-hgvs".to_string(),
        input_file: input.display().to_string(),
        timing,
        sample_results,
        failed_examples,
    };

    save_results(&shard_results, results_output, timing_output)?;

    Ok(shard_results)
}

/// Parse a slice of HGVS patterns in parallel across `workers` threads using a
/// dedicated Rayon thread pool, returning `(successful, failed, elapsed)`.
///
/// ferro parse is pure (`parse_hgvs` holds no provider state), so rayon works
/// without any `Send`/`Sync` complexity. A sized pool is constructed on each
/// call so that the caller controls the degree of parallelism independently of
/// the global pool.
pub fn parse_ferro_count_parallel(
    patterns: &[String],
    workers: usize,
) -> (usize, usize, std::time::Duration) {
    use rayon::prelude::*;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers.max(1))
        .build()
        .expect("rayon pool");
    let start = std::time::Instant::now();
    let (ok, err): (usize, usize) = pool.install(|| {
        patterns
            .par_iter()
            .map(|p| {
                if parse_hgvs(p).is_ok() {
                    (1usize, 0usize)
                } else {
                    (0usize, 1usize)
                }
            })
            .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1))
    });
    (ok, err, start.elapsed())
}

/// Categorize an error for analysis.
fn categorize_error(error: &FerroError) -> String {
    match error {
        FerroError::Parse { msg, .. } => {
            if msg.contains("accession") || msg.contains("Accession") {
                "invalid_accession".to_string()
            } else if msg.contains("position") || msg.contains("coordinate") {
                "invalid_position".to_string()
            } else if msg.contains("edit") || msg.contains("allele") {
                "invalid_edit".to_string()
            } else {
                "parse_error".to_string()
            }
        }
        FerroError::ReferenceNotFound { .. } => "missing_reference".to_string(),
        FerroError::InvalidCoordinates { .. } => "invalid_position".to_string(),
        FerroError::UnsupportedVariant { .. } => "unsupported_variant".to_string(),
        FerroError::ReferenceMismatch { .. } => "reference_mismatch".to_string(),
        _ => "other_error".to_string(),
    }
}

/// Categorize an error based on error message string.
/// Used when converting from commands module results.
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
    } else {
        "parse_error".to_string()
    }
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
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    /// parse_ferro_count_parallel at W=1 and W=4 must produce the same
    /// (ok, err) counts as iterating with parse_hgvs directly, regardless
    /// of timing. Uses patterns that require no reference data (parse is pure).
    #[test]
    fn parse_parallel_matches_serial_counts() {
        let patterns: Vec<String> = vec![
            "NM_000088.3:c.589G>T".to_string(),
            "NC_000001.11:g.12345A>G".to_string(),
            "NM_000088.3:c.1del".to_string(),
            "NP_000079.2:p.Gly197Arg".to_string(), // protein — valid parse
            "not_a_valid_hgvs".to_string(),        // invalid
            "NM_000088.3:c.badposition".to_string(), // invalid
        ];

        // Compute the ground-truth serial counts.
        let serial_ok: usize = patterns.iter().filter(|p| parse_hgvs(p).is_ok()).count();
        let serial_err: usize = patterns.len() - serial_ok;

        for workers in [1, 4] {
            let (ok, err, elapsed) = parse_ferro_count_parallel(&patterns, workers);
            assert_eq!(
                ok, serial_ok,
                "workers={workers}: ok mismatch (got {ok}, expected {serial_ok})"
            );
            assert_eq!(
                err, serial_err,
                "workers={workers}: err mismatch (got {err}, expected {serial_err})"
            );
            assert!(
                elapsed.as_nanos() > 0,
                "workers={workers}: elapsed should be positive"
            );
        }
    }

    #[test]
    fn test_parse_ferro() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let results = dir.path().join("results.json");
        let timing = dir.path().join("timing.json");

        // Create input
        let mut f = File::create(&input).unwrap();
        writeln!(f, "NC_000001.11:g.12345A>G").unwrap();
        writeln!(f, "NM_000088.3:c.589G>T").unwrap();
        writeln!(f, "invalid pattern").unwrap();

        let result = parse_ferro(&input, &results, &timing, 0).unwrap();

        assert_eq!(result.timing.total_patterns, 3);
        assert_eq!(result.timing.successful, 2);
        assert_eq!(result.timing.failed, 1);
        assert!(result.failed_examples.len() == 1);
    }
}
