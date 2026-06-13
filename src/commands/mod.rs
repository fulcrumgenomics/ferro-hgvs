//! Batch commands for parsing and normalizing HGVS variants.
//!
//! This module provides batch processing functionality for the ferro CLI,
//! including support for timing output and parallel processing.

use crate::reference::ReferenceProvider;
use crate::{parse_hgvs, FerroError, MockProvider, MultiFastaProvider, Normalizer};
use indicatif::{ProgressBar, ProgressStyle};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::{Duration, Instant};

/// Result of processing a single variant.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantResult {
    /// Input HGVS string
    pub input: String,

    /// Whether processing succeeded
    pub success: bool,

    /// Output (parsed or normalized) if successful
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<String>,

    /// Error message if failed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// Timing information for a batch operation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingInfo {
    /// Total variants processed
    pub total: usize,

    /// Successfully processed
    pub successful: usize,

    /// Failed
    pub failed: usize,

    /// Total elapsed time in seconds
    pub elapsed_seconds: f64,

    /// Throughput (variants per second)
    pub variants_per_second: f64,
}

impl TimingInfo {
    /// Create timing info from measurements.
    pub fn new(total: usize, successful: usize, elapsed: Duration) -> Self {
        let elapsed_secs = elapsed.as_secs_f64();
        let throughput = if elapsed_secs > f64::EPSILON {
            total as f64 / elapsed_secs
        } else {
            0.0
        };

        Self {
            total,
            successful,
            failed: total - successful,
            elapsed_seconds: elapsed_secs,
            variants_per_second: throughput,
        }
    }
}

/// Results from a batch operation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchResults {
    /// All results
    pub results: Vec<VariantResult>,

    /// Timing information
    pub timing: TimingInfo,
}

/// Configuration for batch normalization.
#[derive(Debug, Clone)]
pub struct NormalizeConfig {
    /// Reference directory (with manifest.json)
    pub reference_dir: Option<std::path::PathBuf>,

    /// Show progress bar
    pub show_progress: bool,

    /// Number of workers for parallel processing
    pub workers: usize,
}

impl Default for NormalizeConfig {
    fn default() -> Self {
        Self {
            reference_dir: None,
            show_progress: true,
            workers: 1,
        }
    }
}

/// Create a reference provider from a directory.
///
/// Handles:
/// - Directory with manifest.json → MultiFastaProvider::from_manifest
/// - Directory without manifest → MultiFastaProvider::from_directory
/// - Falls back to MockProvider with test data
pub fn create_reference_provider(
    reference_dir: Option<&Path>,
) -> Result<Box<dyn ReferenceProvider + Send + Sync>, FerroError> {
    match reference_dir {
        Some(ref_path) => {
            let manifest_path = ref_path.join("manifest.json");
            if manifest_path.exists() {
                eprintln!("Using reference data from {}", ref_path.display());
                let provider = MultiFastaProvider::from_manifest(&manifest_path)?;
                Ok(Box::new(provider))
            } else if ref_path.is_dir() {
                eprintln!("Using reference data from directory {}", ref_path.display());
                let provider = MultiFastaProvider::from_directory(ref_path)?;
                Ok(Box::new(provider))
            } else {
                eprintln!("Reference path not recognized, using test data");
                Ok(Box::new(MockProvider::with_test_data()))
            }
        }
        None => {
            eprintln!("No reference data provided, using mock provider");
            Ok(Box::new(MockProvider::with_test_data()))
        }
    }
}

/// Report which path the cdot cache loads from for a prepared reference dir.
///
/// Resolves the `cdot_json` entry of `<reference_dir>/manifest.json` (relative
/// entries are resolved against the manifest's directory, mirroring
/// [`MultiFastaProvider::from_manifest`]) and loads it via
/// [`CdotMapper::load_with_source`], returning the [`CdotLoadSource`]. Returns
/// `Ok(None)` when there is nothing to report — no manifest, no `cdot_json`
/// entry, or the referenced cdot file is missing.
///
/// `ferro check` uses this to surface a silent fast-path → JSON-fallback
/// regression (#585), and the nightly perf gate uses it as a timing-free
/// co-assertion that the prepared cache is actually on the fast path.
pub fn cdot_cache_load_source(
    reference_dir: &Path,
) -> Result<Option<crate::data::CdotLoadSource>, FerroError> {
    let manifest_path = reference_dir.join("manifest.json");
    if !manifest_path.exists() {
        return Ok(None);
    }

    let base_dir = manifest_path.parent().unwrap_or(Path::new("."));
    let file = File::open(&manifest_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open manifest {}: {}", manifest_path.display(), e),
    })?;
    let manifest: serde_json::Value =
        serde_json::from_reader(file).map_err(|e| FerroError::Io {
            msg: format!(
                "Failed to parse manifest {}: {}",
                manifest_path.display(),
                e
            ),
        })?;

    let Some(cdot_str) = manifest.get("cdot_json").and_then(|v| v.as_str()) else {
        return Ok(None);
    };
    let cdot_path = {
        let p = std::path::PathBuf::from(cdot_str);
        if p.is_absolute() {
            p
        } else {
            base_dir.join(p)
        }
    };
    if !cdot_path.exists() {
        return Ok(None);
    }

    let (_mapper, source) = crate::data::CdotMapper::load_with_source(&cdot_path)?;
    Ok(Some(source))
}

/// Normalize variants from a file.
///
/// # Arguments
/// * `input` - Path to input file (one variant per line)
/// * `config` - Normalization configuration
///
/// # Returns
/// Batch results with all normalize results and timing info.
pub fn normalize_batch<P: AsRef<Path>>(
    input: P,
    config: &NormalizeConfig,
) -> Result<BatchResults, FerroError> {
    let input = input.as_ref();

    // Create normalizer with appropriate reference provider
    let provider = create_reference_provider(config.reference_dir.as_deref())?;
    let normalizer = Normalizer::new(provider);

    // Read variants
    let variants = read_variants(input)?;

    if variants.is_empty() {
        return Ok(BatchResults {
            results: Vec::new(),
            timing: TimingInfo::new(0, 0, Duration::ZERO),
        });
    }

    // Progress bar
    let pb = if config.show_progress {
        let pb = ProgressBar::new(variants.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("##-"),
        );
        Some(pb)
    } else {
        None
    };

    // Normalize all variants
    let start = Instant::now();

    let results: Vec<VariantResult> = variants
        .iter()
        .map(|variant| {
            if let Some(ref pb) = pb {
                pb.inc(1);
            }

            match parse_hgvs(variant) {
                Ok(parsed) => {
                    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        normalizer.normalize(&parsed)
                    })) {
                        Ok(Ok(normalized)) => VariantResult {
                            input: variant.clone(),
                            success: true,
                            output: Some(normalized.to_string()),
                            error: None,
                        },
                        Ok(Err(e)) => VariantResult {
                            input: variant.clone(),
                            success: false,
                            output: None,
                            error: Some(format!("{}", e)),
                        },
                        Err(payload) => {
                            let msg = payload
                                .downcast_ref::<String>()
                                .map(|s| s.as_str())
                                .or_else(|| payload.downcast_ref::<&str>().copied())
                                .unwrap_or("unknown panic");
                            VariantResult {
                                input: variant.clone(),
                                success: false,
                                output: None,
                                error: Some(format!(
                                    "internal error: panic during normalization: {}",
                                    msg
                                )),
                            }
                        }
                    }
                }
                Err(e) => VariantResult {
                    input: variant.clone(),
                    success: false,
                    output: None,
                    error: Some(format!("{}", e)),
                },
            }
        })
        .collect();

    let elapsed = start.elapsed();

    if let Some(pb) = pb {
        pb.finish_with_message(format!(
            "Normalized {} variants in {:.2}s",
            variants.len(),
            elapsed.as_secs_f64()
        ));
    }

    let successful = results.iter().filter(|r| r.success).count();
    let timing = TimingInfo::new(variants.len(), successful, elapsed);

    Ok(BatchResults { results, timing })
}

/// Parse variants from a file (validation only, no normalization).
///
/// # Arguments
/// * `input` - Path to input file (one variant per line)
/// * `show_progress` - Whether to show a progress bar
///
/// # Returns
/// Batch results with all parse results and timing info.
pub fn parse_batch<P: AsRef<Path>>(
    input: P,
    show_progress: bool,
) -> Result<BatchResults, FerroError> {
    let input = input.as_ref();

    // Read variants
    let variants = read_variants(input)?;

    if variants.is_empty() {
        return Ok(BatchResults {
            results: Vec::new(),
            timing: TimingInfo::new(0, 0, Duration::ZERO),
        });
    }

    // Progress bar
    let pb = if show_progress {
        let pb = ProgressBar::new(variants.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("##-"),
        );
        Some(pb)
    } else {
        None
    };

    // Parse all variants
    let start = Instant::now();

    let results: Vec<VariantResult> = variants
        .iter()
        .map(|variant| {
            if let Some(ref pb) = pb {
                pb.inc(1);
            }

            match parse_hgvs(variant) {
                Ok(parsed) => VariantResult {
                    input: variant.clone(),
                    success: true,
                    output: Some(parsed.to_string()),
                    error: None,
                },
                Err(e) => VariantResult {
                    input: variant.clone(),
                    success: false,
                    output: None,
                    error: Some(format!("{}", e)),
                },
            }
        })
        .collect();

    let elapsed = start.elapsed();

    if let Some(pb) = pb {
        pb.finish_with_message(format!(
            "Parsed {} variants in {:.2}s",
            variants.len(),
            elapsed.as_secs_f64()
        ));
    }

    let successful = results.iter().filter(|r| r.success).count();
    let timing = TimingInfo::new(variants.len(), successful, elapsed);

    Ok(BatchResults { results, timing })
}

/// Read variants from a file (one per line).
fn read_variants<P: AsRef<Path>>(path: P) -> Result<Vec<String>, FerroError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;

    let reader = BufReader::new(file);
    let variants: Vec<String> = reader
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();

    Ok(variants)
}

/// Write batch results to a JSON file.
pub fn write_results<P: AsRef<Path>>(results: &BatchResults, path: P) -> Result<(), FerroError> {
    let path = path.as_ref();

    // Create parent directory if needed
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", path.display(), e),
    })?;

    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, results).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    Ok(())
}

/// Write timing info to a JSON file.
pub fn write_timing<P: AsRef<Path>>(timing: &TimingInfo, path: P) -> Result<(), FerroError> {
    let path = path.as_ref();

    // Create parent directory if needed
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory {}: {}", parent.display(), e),
        })?;
    }

    let file = File::create(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", path.display(), e),
    })?;

    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, timing).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    Ok(())
}

/// Write results in text format (one result per line).
pub fn write_results_text<W: Write>(
    results: &BatchResults,
    mut writer: W,
) -> Result<(), FerroError> {
    for result in &results.results {
        if result.success {
            if let Some(ref output) = result.output {
                writeln!(writer, "{}", output).map_err(|e| FerroError::Io {
                    msg: format!("Failed to write output: {}", e),
                })?;
            }
        } else {
            let error = result.error.as_deref().unwrap_or("unknown error");
            writeln!(writer, "ERROR: {} - {}", result.input, error).map_err(|e| {
                FerroError::Io {
                    msg: format!("Failed to write output: {}", e),
                }
            })?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_parse_batch_empty() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        File::create(&input).unwrap();

        let result = parse_batch(&input, false).unwrap();
        assert!(result.results.is_empty());
        assert_eq!(result.timing.total, 0);
    }

    #[test]
    fn test_parse_batch_valid() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let mut f = File::create(&input).unwrap();
        writeln!(f, "NM_000088.3:c.100A>G").unwrap();
        writeln!(f, "NM_000088.3:c.200del").unwrap();

        let result = parse_batch(&input, false).unwrap();
        assert_eq!(result.results.len(), 2);
        assert_eq!(result.timing.total, 2);
        assert!(result.results[0].success);
        assert!(result.results[1].success);
    }

    #[test]
    fn test_parse_batch_invalid() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("input.txt");
        let mut f = File::create(&input).unwrap();
        writeln!(f, "not a valid variant").unwrap();

        let result = parse_batch(&input, false).unwrap();
        assert_eq!(result.results.len(), 1);
        assert!(!result.results[0].success);
        assert!(result.results[0].error.is_some());
    }

    #[test]
    fn test_timing_info() {
        let timing = TimingInfo::new(100, 90, Duration::from_secs(2));
        assert_eq!(timing.total, 100);
        assert_eq!(timing.successful, 90);
        assert_eq!(timing.failed, 10);
        assert!((timing.elapsed_seconds - 2.0).abs() < 0.001);
        assert!((timing.variants_per_second - 50.0).abs() < 0.001);
    }

    const MINIMAL_CDOT_JSON: &str = r#"
    {
        "transcripts": {
            "NM_000088.3": {
                "gene_name": "COL1A1",
                "contig": "NC_000017.11",
                "strand": "+",
                "exons": [[50184096, 50184169, 0, 73]],
                "cds_start": 10,
                "cds_end": 60
            }
        }
    }
    "#;

    #[test]
    fn cdot_cache_load_source_is_none_without_manifest() {
        let dir = TempDir::new().unwrap();
        assert!(cdot_cache_load_source(dir.path()).unwrap().is_none());
    }

    #[test]
    fn cdot_cache_load_source_is_none_without_cdot_entry() {
        let dir = TempDir::new().unwrap();
        std::fs::write(
            dir.path().join("manifest.json"),
            r#"{"genome_fasta": "x.fa"}"#,
        )
        .unwrap();
        assert!(cdot_cache_load_source(dir.path()).unwrap().is_none());
    }

    #[test]
    fn cdot_cache_load_source_reports_archive_for_fresh_cache() {
        use crate::data::CdotLoadSource;
        let dir = TempDir::new().unwrap();
        let cdot_json = dir.path().join("cdot.json");
        std::fs::write(&cdot_json, MINIMAL_CDOT_JSON).unwrap();
        // Prime the cache: the first load self-heals a fresh sibling .rkyv.
        crate::data::CdotMapper::load(&cdot_json).unwrap();
        std::fs::write(
            dir.path().join("manifest.json"),
            r#"{"cdot_json": "cdot.json"}"#,
        )
        .unwrap();

        assert_eq!(
            cdot_cache_load_source(dir.path()).unwrap(),
            Some(CdotLoadSource::Archive),
            "a prepared cdot with a fresh sibling archive must report archive"
        );
    }

    #[test]
    fn cdot_cache_load_source_reports_json_fallback_without_archive() {
        use crate::data::CdotLoadSource;
        let dir = TempDir::new().unwrap();
        let cdot_json = dir.path().join("cdot.json");
        std::fs::write(&cdot_json, MINIMAL_CDOT_JSON).unwrap();
        // No sibling .rkyv exists yet → the load must fall back to JSON.
        std::fs::write(
            dir.path().join("manifest.json"),
            r#"{"cdot_json": "cdot.json"}"#,
        )
        .unwrap();

        assert_eq!(
            cdot_cache_load_source(dir.path()).unwrap(),
            Some(CdotLoadSource::JsonFallback),
            "a cdot with no usable archive must report json-fallback, not archive"
        );
    }
}
