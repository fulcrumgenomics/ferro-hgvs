//! Report generation for benchmark results.

use crate::benchmark::types::{
    AggregateStats, ComparisonSummary, NormalizationComparison, ParsingComparison, SummaryConfig,
};
use crate::FerroError;
use chrono::Utc;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

/// Generate a full benchmark summary from comparison files.
pub fn generate_summary<P: AsRef<Path>>(
    parsing_dir: P,
    normalization_dir: P,
    output: P,
    config: SummaryConfig,
) -> Result<ComparisonSummary, FerroError> {
    let parsing_dir = parsing_dir.as_ref();
    let normalization_dir = normalization_dir.as_ref();
    let output = output.as_ref();

    // Load parsing comparisons
    let mut parsing = HashMap::new();
    if parsing_dir.exists() {
        for entry in std::fs::read_dir(parsing_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to read {}: {}", parsing_dir.display(), e),
        })? {
            let entry = entry.map_err(|e| FerroError::Io {
                msg: format!("Failed to read entry: {}", e),
            })?;
            let path = entry.path();
            if path.extension().is_some_and(|e| e == "json") {
                if let Ok(comp) = load_json::<ParsingComparison>(&path) {
                    parsing.insert(comp.dataset.clone(), comp);
                }
            }
        }
    }

    // Load normalization comparisons
    let mut normalization = HashMap::new();
    if normalization_dir.exists() {
        for entry in std::fs::read_dir(normalization_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to read {}: {}", normalization_dir.display(), e),
        })? {
            let entry = entry.map_err(|e| FerroError::Io {
                msg: format!("Failed to read entry: {}", e),
            })?;
            let path = entry.path();
            if path.extension().is_some_and(|e| e == "json") {
                if let Ok(comp) = load_json::<NormalizationComparison>(&path) {
                    normalization.insert(comp.dataset.clone(), comp);
                }
            }
        }
    }

    // Compute aggregate stats
    let total_patterns: usize = parsing.values().map(|c| c.ferro_hgvs.total_patterns).sum();
    let total_ferro_time: f64 = parsing
        .values()
        .map(|c| c.ferro_hgvs.total_time_seconds)
        .sum();
    let total_mutalyzer_time: f64 = parsing
        .values()
        .filter_map(|c| c.mutalyzer.as_ref())
        .map(|m| m.total_time_seconds)
        .sum();

    let ferro_throughput = if total_ferro_time > 0.0 {
        total_patterns as f64 / total_ferro_time
    } else {
        0.0
    };

    let mutalyzer_throughput = if total_mutalyzer_time > 0.0 {
        Some(total_patterns as f64 / total_mutalyzer_time)
    } else {
        None
    };

    let speedup = mutalyzer_throughput.map(|mt| ferro_throughput / mt);

    let aggregate = AggregateStats {
        total_patterns,
        ferro_throughput,
        mutalyzer_throughput,
        speedup,
    };

    let summary = ComparisonSummary {
        generated: Utc::now(),
        config,
        parsing,
        normalization,
        aggregate,
    };

    // Save summary
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &summary).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    Ok(summary)
}

/// Generate a Markdown report from a benchmark summary.
pub fn generate_report<P: AsRef<Path>>(
    summary: &ComparisonSummary,
    output: P,
) -> Result<(), FerroError> {
    let output = output.as_ref();

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let mut writer = BufWriter::new(file);

    // Header
    writeln!(writer, "# ferro-hgvs vs Mutalyzer Benchmark Report")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "*Generated: {}*",
        summary.generated.format("%Y-%m-%d %H:%M:%S UTC")
    )?;
    writeln!(writer)?;

    // Executive Summary
    writeln!(writer, "## Executive Summary")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "- **Total patterns tested**: {}",
        format_with_commas(summary.aggregate.total_patterns)
    )?;
    writeln!(
        writer,
        "- **ferro-hgvs throughput**: {:.0} patterns/second",
        summary.aggregate.ferro_throughput
    )?;
    if let Some(mt) = summary.aggregate.mutalyzer_throughput {
        writeln!(
            writer,
            "- **Mutalyzer throughput**: {:.0} patterns/second",
            mt
        )?;
    }
    if let Some(speedup) = summary.aggregate.speedup {
        writeln!(writer, "- **Speedup**: {:.0}x faster", speedup)?;
    }
    writeln!(writer)?;

    // Parsing Results
    if !summary.parsing.is_empty() {
        writeln!(writer, "## Parsing Results by Dataset")?;
        writeln!(writer)?;
        writeln!(
            writer,
            "| Dataset | Total | ferro-hgvs | Mutalyzer | Speedup |"
        )?;
        writeln!(
            writer,
            "|---------|------:|:-------------:|:---------:|--------:|"
        )?;

        for (name, comp) in &summary.parsing {
            let mut_rate = comp
                .mutalyzer
                .as_ref()
                .map(|m| format!("{:.2}%", m.pass_rate * 100.0))
                .unwrap_or_else(|| "-".to_string());
            let speedup = comp
                .speedup
                .map(|s| format!("{:.0}x", s))
                .unwrap_or_else(|| "-".to_string());

            writeln!(
                writer,
                "| {} | {} | {:.4}% | {} | {} |",
                name,
                format_with_commas(comp.ferro_hgvs.total_patterns),
                comp.ferro_hgvs.pass_rate * 100.0,
                mut_rate,
                speedup
            )?;
        }
        writeln!(writer)?;
    }

    // Normalization Results
    if !summary.normalization.is_empty() {
        writeln!(writer, "## Normalization Results")?;
        writeln!(writer)?;
        writeln!(
            writer,
            "| Dataset | Sample Size | ferro-hgvs | Mutalyzer | Agreement |"
        )?;
        writeln!(
            writer,
            "|---------|------------:|:-------------:|:---------:|:---------:|"
        )?;

        for (name, comp) in &summary.normalization {
            let mut_rate = comp
                .mutalyzer
                .as_ref()
                .map(|m| format!("{:.2}%", m.pass_rate * 100.0))
                .unwrap_or_else(|| "-".to_string());
            let agreement = comp
                .agreement
                .as_ref()
                .map(|a| format!("{:.2}%", a.agreement_rate * 100.0))
                .unwrap_or_else(|| "-".to_string());

            writeln!(
                writer,
                "| {} | {} | {:.2}% | {} | {} |",
                name,
                format_with_commas(comp.sample_size),
                comp.ferro_hgvs.pass_rate * 100.0,
                mut_rate,
                agreement
            )?;
        }
        writeln!(writer)?;
    }

    // Methodology
    writeln!(writer, "## Methodology")?;
    writeln!(writer)?;
    writeln!(writer, "### Configuration")?;
    writeln!(writer, "- **Cores**: {}", summary.config.cores)?;
    writeln!(
        writer,
        "- **Include Mutalyzer**: {}",
        summary.config.include_mutalyzer
    )?;
    writeln!(
        writer,
        "- **Normalization sample size**: {}",
        format_with_commas(summary.config.normalization_sample_size)
    )?;
    writeln!(writer)?;

    writeln!(writer, "### Parsing Benchmark")?;
    writeln!(
        writer,
        "- Patterns extracted from source files (ClinVar, gnomAD, etc.)"
    )?;
    writeln!(writer, "- Patterns sharded for parallel processing")?;
    writeln!(writer, "- Timing measured for the parse operation only")?;
    writeln!(writer)?;

    writeln!(writer, "### Normalization Benchmark")?;
    writeln!(
        writer,
        "- Stratified sampling to maintain pattern type distribution"
    )?;
    writeln!(writer, "- ferro-hgvs: Uses mock reference provider")?;
    writeln!(writer, "- Mutalyzer: Uses local Docker API")?;
    writeln!(writer)?;

    eprintln!("Generated report: {}", output.display());
    Ok(())
}

/// Generate README-ready markdown tables.
pub fn generate_readme_tables<P: AsRef<Path>>(
    summary: &ComparisonSummary,
    output: P,
) -> Result<(), FerroError> {
    let output = output.as_ref();

    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "<!-- AUTO-GENERATED BY ferro-benchmark -->")?;
    writeln!(
        writer,
        "<!-- Generated: {} -->",
        summary.generated.format("%Y-%m-%d")
    )?;
    writeln!(writer)?;

    // Comparison table
    writeln!(writer, "## Comparison with Other Parsers")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "| Parser | Language | Throughput | Relative Speed |"
    )?;
    writeln!(
        writer,
        "|--------|----------|------------|----------------|"
    )?;
    writeln!(
        writer,
        "| **ferro-hgvs** | Rust | ~{:.0}/s | **1x (baseline)** |",
        summary.aggregate.ferro_throughput
    )?;
    if let Some(mt) = summary.aggregate.mutalyzer_throughput {
        let speedup = summary.aggregate.speedup.unwrap_or(1.0);
        writeln!(
            writer,
            "| mutalyzer-hgvs-parser | Python | ~{:.0}/s | ~{:.0}x slower |",
            mt, speedup
        )?;
    }
    writeln!(writer)?;

    // Validation table
    if !summary.parsing.is_empty() {
        writeln!(writer, "## Large-Scale Validation")?;
        writeln!(writer)?;
        writeln!(writer, "| Source | Total Patterns | ferro-hgvs Pass Rate |")?;
        writeln!(
            writer,
            "|--------|---------------:|:-----------------------:|"
        )?;

        let mut total = 0usize;
        let mut total_success = 0usize;

        for (name, comp) in &summary.parsing {
            writeln!(
                writer,
                "| {} | {} | {:.4}% |",
                name,
                format_with_commas(comp.ferro_hgvs.total_patterns),
                comp.ferro_hgvs.pass_rate * 100.0
            )?;
            total += comp.ferro_hgvs.total_patterns;
            total_success += comp.ferro_hgvs.successful;
        }

        let total_rate = if total > 0 {
            total_success as f64 / total as f64 * 100.0
        } else {
            0.0
        };

        writeln!(
            writer,
            "| **Total** | **{}** | **{:.4}%** |",
            format_with_commas(total),
            total_rate
        )?;
        writeln!(writer)?;
    }

    eprintln!("Generated README tables: {}", output.display());
    Ok(())
}

fn load_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, FerroError> {
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;
    let reader = BufReader::new(file);
    serde_json::from_reader(reader).map_err(|e| FerroError::Json {
        msg: format!("Failed to parse {}: {}", path.display(), e),
    })
}

/// Format a number with commas as thousands separator.
fn format_with_commas(n: usize) -> String {
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
