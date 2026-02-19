//! Result aggregation and comparison.

use crate::benchmark::types::{
    AggregatedResults, AgreementStats, DisagreementExample, NormalizationComparison, ParseResult,
    ParsingComparison, ShardResults, TimingInfo,
};
use crate::FerroError;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

/// Collate parsing results from multiple shards.
pub fn collate_parsing<P: AsRef<Path>>(
    ferro_dir: P,
    mutalyzer_dir: Option<P>,
    output: P,
    dataset_name: &str,
) -> Result<ParsingComparison, FerroError> {
    let ferro_dir = ferro_dir.as_ref();
    let output = output.as_ref();

    // Load ferro-hgvs results
    let ferro_timings = load_timing_files(ferro_dir)?;
    let ferro_results = AggregatedResults::from_timings(&ferro_timings);

    // Load Mutalyzer results if available
    let (mutalyzer_results, speedup, agreement) = if let Some(mut_dir) = mutalyzer_dir {
        let mut_timings = load_timing_files(mut_dir.as_ref())?;
        let mut_results = AggregatedResults::from_timings(&mut_timings);

        let speedup = if mut_results.throughput > 0.0 {
            Some(ferro_results.throughput / mut_results.throughput)
        } else {
            None
        };

        // Compute agreement if we have detailed results
        let agreement = compute_parsing_agreement(ferro_dir, mut_dir.as_ref())?;

        (Some(mut_results), speedup, agreement)
    } else {
        (None, None, None)
    };

    let comparison = ParsingComparison {
        dataset: dataset_name.to_string(),
        ferro_hgvs: ferro_results,
        mutalyzer: mutalyzer_results,
        speedup,
        agreement,
    };

    // Save result
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &comparison).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    eprintln!(
        "Collated parsing results for {}: ferro={:.4}% pass, throughput={:.0}/s",
        dataset_name,
        comparison.ferro_hgvs.pass_rate * 100.0,
        comparison.ferro_hgvs.throughput
    );

    Ok(comparison)
}

/// Collate normalization results.
pub fn collate_normalization<P: AsRef<Path>>(
    ferro_results: P,
    mutalyzer_results: Option<P>,
    output: P,
    dataset_name: &str,
) -> Result<NormalizationComparison, FerroError> {
    let ferro_results_path = ferro_results.as_ref();
    let output = output.as_ref();

    // Load ferro-hgvs results
    let ferro_data: ShardResults = load_json(ferro_results_path)?;
    let ferro_agg = AggregatedResults::from_timings(std::slice::from_ref(&ferro_data.timing));

    // Load Mutalyzer results if available
    let (mutalyzer_agg, agreement) = if let Some(mut_path) = mutalyzer_results {
        let mut_data: ShardResults = load_json(mut_path.as_ref())?;
        let mut_agg = AggregatedResults::from_timings(std::slice::from_ref(&mut_data.timing));

        // Compute agreement
        let agreement = compute_normalization_agreement(&ferro_data, &mut_data);

        (Some(mut_agg), Some(agreement))
    } else {
        (None, None)
    };

    let comparison = NormalizationComparison {
        dataset: dataset_name.to_string(),
        sample_size: ferro_data.timing.total_patterns,
        ferro_hgvs: ferro_agg,
        mutalyzer: mutalyzer_agg,
        agreement,
    };

    // Save result
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &comparison).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    eprintln!(
        "Collated normalization results for {}: ferro={:.2}% pass",
        dataset_name,
        comparison.ferro_hgvs.pass_rate * 100.0
    );

    Ok(comparison)
}

/// Load timing files from a directory.
fn load_timing_files<P: AsRef<Path>>(dir: P) -> Result<Vec<TimingInfo>, FerroError> {
    let dir = dir.as_ref();
    let mut timings = Vec::new();

    for entry in std::fs::read_dir(dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to read directory {}: {}", dir.display(), e),
    })? {
        let entry = entry.map_err(|e| FerroError::Io {
            msg: format!("Failed to read entry: {}", e),
        })?;
        let path = entry.path();

        if path.extension().is_some_and(|e| e == "json")
            && path
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.contains("timing"))
        {
            let timing: TimingInfo = load_json(&path)?;
            timings.push(timing);
        }
    }

    Ok(timings)
}

/// Load a JSON file.
fn load_json<T: serde::de::DeserializeOwned, P: AsRef<Path>>(path: P) -> Result<T, FerroError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;
    let reader = BufReader::new(file);
    serde_json::from_reader(reader).map_err(|e| FerroError::Json {
        msg: format!("Failed to parse {}: {}", path.display(), e),
    })
}

/// Compute agreement statistics between parsing results.
fn compute_parsing_agreement<P: AsRef<Path>>(
    ferro_dir: P,
    mutalyzer_dir: P,
) -> Result<Option<AgreementStats>, FerroError> {
    let ferro_dir = ferro_dir.as_ref();
    let mutalyzer_dir = mutalyzer_dir.as_ref();

    // Load results files (not timing files)
    let ferro_results = load_results_files(ferro_dir)?;
    let mutalyzer_results = load_results_files(mutalyzer_dir)?;

    if ferro_results.is_empty() || mutalyzer_results.is_empty() {
        return Ok(None);
    }

    // Build lookup maps by input
    let ferro_map: HashMap<&str, &ParseResult> = ferro_results
        .iter()
        .flat_map(|sr| sr.sample_results.iter())
        .chain(
            ferro_results
                .iter()
                .flat_map(|sr| sr.failed_examples.iter()),
        )
        .map(|r| (r.input.as_str(), r))
        .collect();

    let mutalyzer_map: HashMap<&str, &ParseResult> = mutalyzer_results
        .iter()
        .flat_map(|sr| sr.sample_results.iter())
        .chain(
            mutalyzer_results
                .iter()
                .flat_map(|sr| sr.failed_examples.iter()),
        )
        .map(|r| (r.input.as_str(), r))
        .collect();

    // Compute agreement
    let mut both_success = 0usize;
    let mut both_fail = 0usize;
    let mut ferro_only_success = 0usize;
    let mut mutalyzer_only_success = 0usize;
    let mut agreements = 0usize;
    let mut disagreement_examples = Vec::new();

    for (input, ferro_result) in &ferro_map {
        if let Some(mut_result) = mutalyzer_map.get(input) {
            match (ferro_result.success, mut_result.success) {
                (true, true) => {
                    both_success += 1;
                    // Check if outputs match
                    if ferro_result.output == mut_result.output {
                        agreements += 1;
                    } else if disagreement_examples.len() < 50 {
                        disagreement_examples.push(DisagreementExample {
                            input: input.to_string(),
                            ferro_output: ferro_result.output.clone().unwrap_or_default(),
                            mutalyzer_output: mut_result.output.clone().unwrap_or_default(),
                        });
                    }
                }
                (true, false) => ferro_only_success += 1,
                (false, true) => mutalyzer_only_success += 1,
                (false, false) => both_fail += 1,
            }
        }
    }

    let agreement_rate = if both_success > 0 {
        agreements as f64 / both_success as f64
    } else {
        0.0
    };

    Ok(Some(AgreementStats {
        both_success,
        both_fail,
        ferro_only_success,
        mutalyzer_only_success,
        agreements,
        agreement_rate,
        disagreement_examples,
    }))
}

/// Load results files from a directory.
fn load_results_files<P: AsRef<Path>>(dir: P) -> Result<Vec<ShardResults>, FerroError> {
    let dir = dir.as_ref();
    let mut results = Vec::new();

    for entry in std::fs::read_dir(dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to read directory {}: {}", dir.display(), e),
    })? {
        let entry = entry.map_err(|e| FerroError::Io {
            msg: format!("Failed to read entry: {}", e),
        })?;
        let path = entry.path();

        if path.extension().is_some_and(|e| e == "json")
            && !path
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.contains("timing"))
        {
            if let Ok(result) = load_json::<ShardResults, _>(&path) {
                results.push(result);
            }
        }
    }

    Ok(results)
}

/// Compute agreement between normalization results.
fn compute_normalization_agreement(
    ferro_data: &ShardResults,
    mutalyzer_data: &ShardResults,
) -> AgreementStats {
    // Build lookup maps
    let ferro_map: HashMap<&str, &ParseResult> = ferro_data
        .sample_results
        .iter()
        .chain(ferro_data.failed_examples.iter())
        .map(|r| (r.input.as_str(), r))
        .collect();

    let mutalyzer_map: HashMap<&str, &ParseResult> = mutalyzer_data
        .sample_results
        .iter()
        .chain(mutalyzer_data.failed_examples.iter())
        .map(|r| (r.input.as_str(), r))
        .collect();

    let mut both_success = 0usize;
    let mut both_fail = 0usize;
    let mut ferro_only_success = 0usize;
    let mut mutalyzer_only_success = 0usize;
    let mut agreements = 0usize;
    let mut disagreement_examples = Vec::new();

    for (input, ferro_result) in &ferro_map {
        if let Some(mut_result) = mutalyzer_map.get(input) {
            match (ferro_result.success, mut_result.success) {
                (true, true) => {
                    both_success += 1;
                    if ferro_result.output == mut_result.output {
                        agreements += 1;
                    } else if disagreement_examples.len() < 50 {
                        disagreement_examples.push(DisagreementExample {
                            input: input.to_string(),
                            ferro_output: ferro_result.output.clone().unwrap_or_default(),
                            mutalyzer_output: mut_result.output.clone().unwrap_or_default(),
                        });
                    }
                }
                (true, false) => ferro_only_success += 1,
                (false, true) => mutalyzer_only_success += 1,
                (false, false) => both_fail += 1,
            }
        }
    }

    let agreement_rate = if both_success > 0 {
        agreements as f64 / both_success as f64
    } else {
        0.0
    };

    AgreementStats {
        both_success,
        both_fail,
        ferro_only_success,
        mutalyzer_only_success,
        agreements,
        agreement_rate,
        disagreement_examples,
    }
}
