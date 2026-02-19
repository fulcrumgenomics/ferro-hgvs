//! Workflow orchestration for comparison benchmarks.

use crate::benchmark::collate::{collate_normalization, collate_parsing};
use crate::benchmark::extract::{copy_text, extract_clinvar, extract_json};
use crate::benchmark::mutalyzer::{
    has_mutalyzer_parser, run_mutalyzer_parser_subprocess, MutalyzerClient,
};
use crate::benchmark::normalize::{normalize_ferro, normalize_mutalyzer};
use crate::benchmark::parse::parse_ferro;
use crate::benchmark::report::{generate_readme_tables, generate_report, generate_summary};
use crate::benchmark::sample::stratified_sample;
use crate::benchmark::shard::shard_dataset;
use crate::benchmark::types::{ComparisonConfig, DatasetConfig, DatasetFormat, SummaryConfig};
use crate::FerroError;
use rayon::prelude::*;
use std::fs;
use std::path::{Path, PathBuf};

/// Run the complete comparison workflow.
pub fn run_comparison(config: &ComparisonConfig) -> Result<(), FerroError> {
    let output_dir = &config.output_dir;
    fs::create_dir_all(output_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create output directory: {}", e),
    })?;

    eprintln!("Starting benchmark with {} cores", config.cores);
    eprintln!("Output directory: {}", output_dir.display());

    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.cores)
        .build_global()
        .ok(); // Ignore if already initialized

    // Phase 1: Extract patterns from each dataset
    eprintln!("\n=== Phase 1: Extraction ===");
    for dataset in &config.datasets {
        extract_dataset(dataset, output_dir)?;
    }

    // Phase 2: Shard datasets for parallel processing
    eprintln!("\n=== Phase 2: Sharding ===");
    for dataset in &config.datasets {
        shard_dataset_workflow(&dataset.name, output_dir, config.cores)?;
    }

    // Phase 3: Parse with ferro-hgvs
    eprintln!("\n=== Phase 3: Parsing (ferro-hgvs) ===");
    for dataset in &config.datasets {
        parse_ferro_workflow(&dataset.name, output_dir, config.cores)?;
    }

    // Phase 4: Parse with Mutalyzer (if enabled)
    if config.include_mutalyzer && has_mutalyzer_parser() {
        eprintln!("\n=== Phase 4: Parsing (Mutalyzer) ===");
        for dataset in &config.datasets {
            parse_mutalyzer_workflow(&dataset.name, output_dir, config.cores)?;
        }
    } else if config.include_mutalyzer {
        eprintln!("\n=== Phase 4: Skipping Mutalyzer parsing (not available) ===");
    }

    // Phase 5: Collate parsing results
    eprintln!("\n=== Phase 5: Collating parsing results ===");
    for dataset in &config.datasets {
        collate_parsing_workflow(&dataset.name, output_dir, config.include_mutalyzer)?;
    }

    // Phase 6: Sample for normalization
    eprintln!("\n=== Phase 6: Sampling for normalization ===");
    for dataset in &config.datasets {
        sample_for_normalization(&dataset.name, output_dir, config.normalization_sample_size)?;
    }

    // Phase 7: Normalize with ferro-hgvs
    eprintln!("\n=== Phase 7: Normalization (ferro-hgvs) ===");
    for dataset in &config.datasets {
        normalize_ferro_workflow(&dataset.name, output_dir, config.reference_path.as_ref())?;
    }

    // Phase 8: Normalize with Mutalyzer API (if enabled)
    if config.include_mutalyzer {
        if let Some(ref api_url) = config.mutalyzer_api_url {
            eprintln!("\n=== Phase 8: Normalization (Mutalyzer API) ===");
            for dataset in &config.datasets {
                normalize_mutalyzer_workflow(&dataset.name, output_dir, api_url)?;
            }
        } else {
            eprintln!("\n=== Phase 8: Skipping Mutalyzer normalization (no API URL) ===");
        }
    }

    // Phase 9: Collate normalization results
    eprintln!("\n=== Phase 9: Collating normalization results ===");
    for dataset in &config.datasets {
        collate_normalization_workflow(&dataset.name, output_dir, config.include_mutalyzer)?;
    }

    // Phase 10: Generate summary and reports
    eprintln!("\n=== Phase 10: Generating reports ===");
    let summary_config = SummaryConfig {
        cores: config.cores,
        include_mutalyzer: config.include_mutalyzer,
        normalization_sample_size: config.normalization_sample_size,
    };

    let parsing_dir = output_dir.join("comparisons").join("parsing");
    let normalization_dir = output_dir.join("comparisons").join("normalization");
    let summary_path = output_dir.join("summary.json");

    let summary = generate_summary(
        &parsing_dir,
        &normalization_dir,
        &summary_path,
        summary_config,
    )?;

    let report_path = output_dir.join("report.md");
    generate_report(&summary, &report_path)?;

    let readme_tables_path = output_dir.join("readme_tables.md");
    generate_readme_tables(&summary, &readme_tables_path)?;

    eprintln!("\n=== Comparison Complete ===");
    eprintln!("Summary: {}", summary_path.display());
    eprintln!("Report: {}", report_path.display());
    eprintln!("README tables: {}", readme_tables_path.display());

    Ok(())
}

/// Extract patterns from a dataset.
fn extract_dataset(dataset: &DatasetConfig, output_dir: &Path) -> Result<(), FerroError> {
    let output_path = output_dir
        .join("extracted")
        .join(format!("{}.txt", dataset.name));

    // Check if already extracted (timestamp-based caching)
    if is_up_to_date(&dataset.source, &output_path) {
        eprintln!("  {} already extracted, skipping", dataset.name);
        return Ok(());
    }

    eprintln!(
        "  Extracting {} from {}",
        dataset.name,
        dataset.source.display()
    );

    match dataset.format {
        DatasetFormat::ClinvarTsv => {
            extract_clinvar(&dataset.source, &output_path)?;
        }
        DatasetFormat::TestCasesJson => {
            extract_json(&dataset.source, &output_path, dataset.format, "input")?;
        }
        DatasetFormat::JsonArray => {
            extract_json(&dataset.source, &output_path, dataset.format, "hgvs")?;
        }
        DatasetFormat::PlainText => {
            copy_text(&dataset.source, &output_path)?;
        }
    }

    Ok(())
}

/// Shard a dataset for parallel processing.
fn shard_dataset_workflow(
    dataset_name: &str,
    output_dir: &Path,
    num_shards: usize,
) -> Result<(), FerroError> {
    let input = output_dir
        .join("extracted")
        .join(format!("{}.txt", dataset_name));
    let shard_dir = output_dir.join("shards").join(dataset_name);

    // Check if already sharded
    let first_shard = shard_dir.join("shard_0.txt");
    if is_up_to_date(&input, &first_shard) {
        eprintln!("  {} already sharded, skipping", dataset_name);
        return Ok(());
    }

    eprintln!("  Sharding {} into {} shards", dataset_name, num_shards);
    shard_dataset(&input, &shard_dir, num_shards)?;

    Ok(())
}

/// Parse all shards with ferro-hgvs.
fn parse_ferro_workflow(
    dataset_name: &str,
    output_dir: &Path,
    num_shards: usize,
) -> Result<(), FerroError> {
    let shard_dir = output_dir.join("shards").join(dataset_name);
    let results_dir = output_dir
        .join("results")
        .join("parsing")
        .join("ferro")
        .join(dataset_name);

    fs::create_dir_all(&results_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create results directory: {}", e),
    })?;

    // Collect shards to process
    let shards: Vec<(usize, PathBuf)> = (0..num_shards)
        .map(|i| (i, shard_dir.join(format!("shard_{}.txt", i))))
        .filter(|(_, path)| path.exists())
        .collect();

    eprintln!(
        "  Parsing {} with ferro-hgvs ({} shards)",
        dataset_name,
        shards.len()
    );

    // Process shards in parallel
    shards.par_iter().try_for_each(|(shard_idx, shard_path)| {
        let results_path = results_dir.join(format!("shard_{}_results.json", shard_idx));
        let timing_path = results_dir.join(format!("shard_{}_timing.json", shard_idx));

        // Skip if already processed
        if is_up_to_date(shard_path, &timing_path) {
            return Ok(());
        }

        parse_ferro(shard_path, &results_path, &timing_path, *shard_idx)?;
        Ok::<_, FerroError>(())
    })?;

    Ok(())
}

/// Parse all shards with Mutalyzer parser (via Python subprocess).
fn parse_mutalyzer_workflow(
    dataset_name: &str,
    output_dir: &Path,
    num_shards: usize,
) -> Result<(), FerroError> {
    let shard_dir = output_dir.join("shards").join(dataset_name);
    let results_dir = output_dir
        .join("results")
        .join("parsing")
        .join("mutalyzer")
        .join(dataset_name);

    fs::create_dir_all(&results_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create results directory: {}", e),
    })?;

    // Collect shards to process
    let shards: Vec<(usize, PathBuf)> = (0..num_shards)
        .map(|i| (i, shard_dir.join(format!("shard_{}.txt", i))))
        .filter(|(_, path)| path.exists())
        .collect();

    eprintln!(
        "  Parsing {} with Mutalyzer ({} shards)",
        dataset_name,
        shards.len()
    );

    // Process shards sequentially (Python subprocess)
    for (shard_idx, shard_path) in &shards {
        let output_path = results_dir.join(format!("shard_{}_results.json", shard_idx));

        // Skip if already processed
        if is_up_to_date(shard_path, &output_path) {
            continue;
        }

        run_mutalyzer_parser_subprocess(
            &shard_path.display().to_string(),
            &output_path.display().to_string(),
        )?;
    }

    Ok(())
}

/// Collate parsing results for a dataset.
fn collate_parsing_workflow(
    dataset_name: &str,
    output_dir: &Path,
    include_mutalyzer: bool,
) -> Result<(), FerroError> {
    let ferro_dir = output_dir
        .join("results")
        .join("parsing")
        .join("ferro")
        .join(dataset_name);

    let mutalyzer_dir = if include_mutalyzer {
        let dir = output_dir
            .join("results")
            .join("parsing")
            .join("mutalyzer")
            .join(dataset_name);
        if dir.exists() {
            Some(dir)
        } else {
            None
        }
    } else {
        None
    };

    let output_path = output_dir
        .join("comparisons")
        .join("parsing")
        .join(format!("{}.json", dataset_name));

    eprintln!("  Collating parsing results for {}", dataset_name);
    collate_parsing(
        &ferro_dir,
        mutalyzer_dir.as_ref(),
        &output_path,
        dataset_name,
    )?;

    Ok(())
}

/// Sample patterns for normalization testing.
fn sample_for_normalization(
    dataset_name: &str,
    output_dir: &Path,
    sample_size: usize,
) -> Result<(), FerroError> {
    let input = output_dir
        .join("extracted")
        .join(format!("{}.txt", dataset_name));
    let output = output_dir
        .join("samples")
        .join(format!("{}_sample.txt", dataset_name));

    // Skip if already sampled
    if is_up_to_date(&input, &output) {
        eprintln!("  {} already sampled, skipping", dataset_name);
        return Ok(());
    }

    eprintln!("  Sampling {} patterns from {}", sample_size, dataset_name);
    // Exclude protein patterns for normalization samples since protein normalization
    // is handled inconsistently across tools
    stratified_sample(&input, &output, sample_size, 42, true)?;

    Ok(())
}

/// Normalize samples with ferro-hgvs.
fn normalize_ferro_workflow(
    dataset_name: &str,
    output_dir: &Path,
    reference_dir: Option<&PathBuf>,
) -> Result<(), FerroError> {
    let input = output_dir
        .join("samples")
        .join(format!("{}_sample.txt", dataset_name));
    let results_path = output_dir
        .join("results")
        .join("normalization")
        .join("ferro")
        .join(format!("{}_results.json", dataset_name));
    let timing_path = output_dir
        .join("results")
        .join("normalization")
        .join("ferro")
        .join(format!("{}_timing.json", dataset_name));

    if !input.exists() {
        eprintln!("  Skipping {} normalization (no sample)", dataset_name);
        return Ok(());
    }

    // Skip if already processed
    if is_up_to_date(&input, &timing_path) {
        eprintln!("  {} already normalized (ferro), skipping", dataset_name);
        return Ok(());
    }

    eprintln!("  Normalizing {} with ferro-hgvs", dataset_name);
    normalize_ferro(
        input,
        results_path,
        timing_path,
        reference_dir.map(|p| p.to_path_buf()),
    )?;

    Ok(())
}

/// Normalize samples with Mutalyzer API.
fn normalize_mutalyzer_workflow(
    dataset_name: &str,
    output_dir: &Path,
    api_url: &str,
) -> Result<(), FerroError> {
    let input = output_dir
        .join("samples")
        .join(format!("{}_sample.txt", dataset_name));
    let results_path = output_dir
        .join("results")
        .join("normalization")
        .join("mutalyzer")
        .join(format!("{}_results.json", dataset_name));
    let timing_path = output_dir
        .join("results")
        .join("normalization")
        .join("mutalyzer")
        .join(format!("{}_timing.json", dataset_name));

    if !input.exists() {
        eprintln!("  Skipping {} normalization (no sample)", dataset_name);
        return Ok(());
    }

    // Skip if already processed
    if is_up_to_date(&input, &timing_path) {
        eprintln!(
            "  {} already normalized (mutalyzer), skipping",
            dataset_name
        );
        return Ok(());
    }

    // Check API availability
    let client = MutalyzerClient::new(api_url)?;
    if !client.health_check()? {
        eprintln!("  Mutalyzer API not available, skipping {}", dataset_name);
        return Ok(());
    }

    eprintln!("  Normalizing {} with Mutalyzer API", dataset_name);
    normalize_mutalyzer(&input, &results_path, &timing_path, api_url, Some(100))?;

    Ok(())
}

/// Collate normalization results for a dataset.
fn collate_normalization_workflow(
    dataset_name: &str,
    output_dir: &Path,
    include_mutalyzer: bool,
) -> Result<(), FerroError> {
    let ferro_results = output_dir
        .join("results")
        .join("normalization")
        .join("ferro")
        .join(format!("{}_results.json", dataset_name));

    if !ferro_results.exists() {
        eprintln!(
            "  Skipping {} normalization collation (no results)",
            dataset_name
        );
        return Ok(());
    }

    let mutalyzer_results = if include_mutalyzer {
        let path = output_dir
            .join("results")
            .join("normalization")
            .join("mutalyzer")
            .join(format!("{}_results.json", dataset_name));
        if path.exists() {
            Some(path)
        } else {
            None
        }
    } else {
        None
    };

    let output_path = output_dir
        .join("comparisons")
        .join("normalization")
        .join(format!("{}.json", dataset_name));

    eprintln!("  Collating normalization results for {}", dataset_name);
    collate_normalization(
        &ferro_results,
        mutalyzer_results.as_ref(),
        &output_path,
        dataset_name,
    )?;

    Ok(())
}

/// Check if output is up-to-date relative to input (timestamp-based caching).
fn is_up_to_date<P: AsRef<Path>, Q: AsRef<Path>>(input: P, output: Q) -> bool {
    let input = input.as_ref();
    let output = output.as_ref();

    if !output.exists() {
        return false;
    }

    let input_modified = match fs::metadata(input).and_then(|m| m.modified()) {
        Ok(t) => t,
        Err(_) => return false,
    };

    let output_modified = match fs::metadata(output).and_then(|m| m.modified()) {
        Ok(t) => t,
        Err(_) => return false,
    };

    output_modified >= input_modified
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_is_up_to_date() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let output = dir.path().join("output.txt");

        // Neither exists
        assert!(!is_up_to_date(&input, &output));

        // Only input exists
        fs::write(&input, "test").unwrap();
        assert!(!is_up_to_date(&input, &output));

        // Both exist, output newer
        std::thread::sleep(std::time::Duration::from_millis(10));
        fs::write(&output, "test").unwrap();
        assert!(is_up_to_date(&input, &output));

        // Input updated
        std::thread::sleep(std::time::Duration::from_millis(10));
        fs::write(&input, "updated").unwrap();
        assert!(!is_up_to_date(&input, &output));
    }
}
