//! Native Rust hgvs-rs integration for HGVS normalization comparison.
//!
//! This module provides a wrapper around the hgvs-rs crate for normalizing HGVS variants
//! using the same UTA and SeqRepo infrastructure as biocommons/hgvs.
//!
//! # Usage
//!
//! ```bash
//! # Check if hgvs-rs is properly configured
//! ferro-benchmark check-hgvs-rs --uta-db-url postgresql://... --seqrepo-path /path/to/seqrepo
//!
//! # Compare normalization with hgvs-rs
//! ferro-benchmark compare-normalize --validator hgvs-rs --uta-db-url postgresql://...
//! ```

use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;

use hgvs::data::uta_sr::{Config as UtaSrConfig, Provider};
use hgvs::mapper::variant::{Config as MapperConfig, Mapper};
use hgvs::normalizer::{Config as NormalizerConfig, Direction, Normalizer};
use hgvs::parser::{HgvsVariant, NoRef};
use hgvs::validator::IntrinsicValidator;

use crate::FerroError;

/// Result type for hgvs-rs normalization batch operations.
/// Contains (results, elapsed_time, error_counts).
pub type HgvsRsNormalizeResult = (
    Vec<crate::benchmark::ParseResult>,
    std::time::Duration,
    std::collections::HashMap<String, usize>,
);

/// Configuration for hgvs-rs normalization.
#[derive(Debug, Clone)]
pub struct HgvsRsConfig {
    /// PostgreSQL connection URL for UTA database
    pub uta_db_url: String,
    /// UTA database schema (e.g., "uta_20210129")
    pub uta_db_schema: String,
    /// Path to SeqRepo directory (e.g., "/usr/local/share/seqrepo/2021-01-29")
    pub seqrepo_path: String,
    /// Optional LRG-to-RefSeq mapping file for LRG transcript translation
    pub lrg_mapping_file: Option<String>,
}

impl Default for HgvsRsConfig {
    fn default() -> Self {
        Self {
            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta".to_string(),
            uta_db_schema: "uta_20210129b".to_string(),
            seqrepo_path: "/usr/local/share/seqrepo/2021-01-29".to_string(),
            lrg_mapping_file: None,
        }
    }
}

/// LRG-to-RefSeq mapping for pattern translation.
#[derive(Debug, Clone, Default)]
pub struct LrgMapping {
    mapping: std::collections::HashMap<String, String>,
}

impl LrgMapping {
    /// Load LRG-to-RefSeq mapping from a file.
    ///
    /// The file format is tab-separated with columns:
    /// LRG, HGNC_SYMBOL, REFSEQ_GENOMIC, LRG_TRANSCRIPT, REFSEQ_TRANSCRIPT, ...
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        use std::io::BufRead;

        let file = std::fs::File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open LRG mapping file: {}", e),
        })?;
        let reader = std::io::BufReader::new(file);

        let mut mapping = std::collections::HashMap::new();
        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read LRG mapping line: {}", e),
            })?;

            // Skip comments and empty lines
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                continue;
            }

            // Fields: LRG, HGNC_SYMBOL, REFSEQ_GENOMIC, LRG_TRANSCRIPT, REFSEQ_TRANSCRIPT
            let lrg_id = fields[0]; // e.g., "LRG_1"
            let lrg_transcript = fields[3]; // e.g., "t1"
            let refseq_transcript = fields[4]; // e.g., "NM_000088.3"

            // Skip if RefSeq transcript is empty or "-"
            if refseq_transcript.is_empty() || refseq_transcript == "-" {
                continue;
            }

            // Create full LRG transcript ID: "LRG_1t1"
            let lrg_full = format!("{}{}", lrg_id, lrg_transcript);
            mapping.insert(lrg_full, refseq_transcript.to_string());
        }

        Ok(Self { mapping })
    }

    /// Translate an LRG transcript pattern to RefSeq.
    ///
    /// Returns (translated_pattern, original_lrg_accession) if translation occurred,
    /// or (original_pattern, None) if no translation was needed.
    pub fn translate_pattern(&self, pattern: &str) -> (String, Option<String>) {
        // Match LRG transcript patterns like "LRG_199t1:c.79del"
        if let Some(colon_pos) = pattern.find(':') {
            let accession = &pattern[..colon_pos];
            if accession.starts_with("LRG_") && accession.contains('t') {
                if let Some(refseq) = self.mapping.get(accession) {
                    let translated = format!("{}{}", refseq, &pattern[colon_pos..]);
                    return (translated, Some(accession.to_string()));
                }
            }
        }
        (pattern.to_string(), None)
    }
}

/// Result of normalizing a single variant with hgvs-rs.
#[derive(Debug, Clone)]
pub struct HgvsRsResult {
    /// The input HGVS string
    pub input: String,
    /// Whether normalization succeeded
    pub success: bool,
    /// The normalized output (if successful)
    pub output: Option<String>,
    /// Error message (if failed)
    pub error: Option<String>,
}

/// A wrapper around hgvs-rs components for normalization.
pub struct HgvsRsNormalizer {
    mapper: Mapper,
}

impl HgvsRsNormalizer {
    /// Create a new hgvs-rs normalizer with the given configuration.
    pub fn new(config: &HgvsRsConfig) -> Result<Self, FerroError> {
        let uta_sr_config = UtaSrConfig {
            db_url: config.uta_db_url.clone(),
            db_schema: config.uta_db_schema.clone(),
            seqrepo_path: config.seqrepo_path.clone(),
        };

        let provider = Arc::new(Provider::new(uta_sr_config).map_err(|e| FerroError::Io {
            msg: format!("Failed to create hgvs-rs provider: {}", e),
        })?);

        let mapper_config = MapperConfig::default();
        let mapper = Mapper::new(&mapper_config, provider);

        Ok(Self { mapper })
    }

    /// Normalize a single HGVS variant string.
    pub fn normalize(&self, hgvs_str: &str) -> HgvsRsResult {
        use std::panic::{self, AssertUnwindSafe};

        let input = hgvs_str.to_string();

        // Parse the HGVS string
        let parsed = match HgvsVariant::from_str(hgvs_str) {
            Ok(var) => var,
            Err(e) => {
                return HgvsRsResult {
                    input,
                    success: false,
                    output: None,
                    error: Some(format!("Parse error: {}", e)),
                };
            }
        };

        // Create normalizer with 3' shifting (standard HGVS)
        let validator = Arc::new(IntrinsicValidator::new(true));
        let normalizer = Normalizer::new(
            &self.mapper,
            self.mapper.provider(),
            validator,
            NormalizerConfig {
                shuffle_direction: Direction::FiveToThree,
                cross_boundaries: false,
                replace_reference: true,
                ..Default::default()
            },
        );

        // Normalize the variant, catching any panics from upstream hgvs crate
        let normalize_result =
            panic::catch_unwind(AssertUnwindSafe(|| normalizer.normalize(&parsed)));

        match normalize_result {
            Ok(Ok(normalized)) => {
                let output = format!("{}", NoRef(&normalized));
                HgvsRsResult {
                    input,
                    success: true,
                    output: Some(output),
                    error: None,
                }
            }
            Ok(Err(e)) => HgvsRsResult {
                input,
                success: false,
                output: None,
                error: Some(format!("Normalization error: {}", e)),
            },
            Err(panic_info) => {
                // Extract panic message if possible
                let panic_msg = if let Some(s) = panic_info.downcast_ref::<&str>() {
                    s.to_string()
                } else if let Some(s) = panic_info.downcast_ref::<String>() {
                    s.clone()
                } else {
                    "unknown panic".to_string()
                };
                HgvsRsResult {
                    input,
                    success: false,
                    output: None,
                    error: Some(format!("Normalization panic: {}", panic_msg)),
                }
            }
        }
    }

    /// Normalize multiple HGVS variant strings.
    pub fn normalize_batch(&self, hgvs_strs: &[String]) -> Vec<HgvsRsResult> {
        hgvs_strs.iter().map(|s| self.normalize(s)).collect()
    }
}

/// Check if hgvs-rs can connect to UTA and SeqRepo.
pub fn check_hgvs_rs_available(config: &HgvsRsConfig) -> Result<(), FerroError> {
    // Try to create a provider - this will fail if UTA or SeqRepo is not accessible
    let uta_sr_config = UtaSrConfig {
        db_url: config.uta_db_url.clone(),
        db_schema: config.uta_db_schema.clone(),
        seqrepo_path: config.seqrepo_path.clone(),
    };

    Provider::new(uta_sr_config).map_err(|e| FerroError::Io {
        msg: format!("hgvs-rs not available: {}", e),
    })?;

    Ok(())
}

/// Run hgvs-rs normalization on multiple patterns (sequential).
///
/// For parallel normalization, use `run_hgvs_rs_normalize_parallel`.
/// Returns (results, elapsed_time, error_counts).
pub fn run_hgvs_rs_normalize(
    patterns: &[String],
    config: &HgvsRsConfig,
) -> Result<HgvsRsNormalizeResult, crate::FerroError> {
    use std::collections::HashMap;
    use std::time::Instant;

    // Load LRG mapping if provided
    let lrg_mapping = if let Some(ref path) = config.lrg_mapping_file {
        Some(LrgMapping::load(path)?)
    } else {
        None
    };

    let start = Instant::now();

    // Create the normalizer
    let normalizer = HgvsRsNormalizer::new(config)?;

    // Run normalization sequentially
    let results: Vec<crate::benchmark::ParseResult> = patterns
        .iter()
        .map(|pattern| {
            // Translate LRG patterns to RefSeq if mapping available
            let (translated_pattern, original_lrg) = if let Some(ref mapping) = lrg_mapping {
                mapping.translate_pattern(pattern)
            } else {
                (pattern.clone(), None)
            };

            let result = normalizer.normalize(&translated_pattern);

            // If we translated an LRG pattern, report with original accession
            let output = if original_lrg.is_some() && result.success {
                result.output.map(|o| {
                    if let Some(ref lrg) = original_lrg {
                        if let Some(colon_pos) = o.find(':') {
                            format!("{}{}", lrg, &o[colon_pos..])
                        } else {
                            o
                        }
                    } else {
                        o
                    }
                })
            } else {
                result.output
            };

            crate::benchmark::ParseResult {
                input: pattern.clone(), // Always report original input
                success: result.success,
                output,
                error: result.error.clone(),
                error_category: result.error.map(|e| categorize_hgvs_rs_error(&e)),
                ref_mismatch: None,
                details: None,
            }
        })
        .collect();

    let elapsed = start.elapsed();

    // Aggregate error counts
    let mut error_counts: HashMap<String, usize> = HashMap::new();
    for result in &results {
        if let Some(ref category) = result.error_category {
            *error_counts.entry(category.clone()).or_insert(0) += 1;
        }
    }

    Ok((results, elapsed, error_counts))
}

/// Run hgvs-rs normalization on multiple patterns in parallel.
///
/// Uses Rayon to parallelize across workers, with each worker creating its own
/// database connection. This achieves significant speedup for bulk normalization.
///
/// Returns (results, elapsed_time, error_counts).
pub fn run_hgvs_rs_normalize_parallel(
    patterns: &[String],
    config: &HgvsRsConfig,
    workers: usize,
) -> Result<HgvsRsNormalizeResult, crate::FerroError> {
    use rayon::prelude::*;
    use std::collections::HashMap;
    use std::time::Instant;

    if workers <= 1 {
        return run_hgvs_rs_normalize(patterns, config);
    }

    // Load LRG mapping if provided
    let lrg_mapping = if let Some(ref path) = config.lrg_mapping_file {
        Some(LrgMapping::load(path)?)
    } else {
        None
    };

    let start = Instant::now();

    // Calculate chunk size for parallel processing
    let chunk_size = patterns.len().div_ceil(workers).max(1);

    // Process chunks in parallel, each chunk creates its own normalizer
    // This avoids Tokio runtime conflicts with thread-local storage
    let chunk_results: Vec<Vec<(usize, crate::benchmark::ParseResult)>> = patterns
        .chunks(chunk_size)
        .enumerate()
        .collect::<Vec<_>>()
        .into_par_iter()
        .map(|(chunk_idx, chunk)| {
            // Create normalizer for this chunk
            let normalizer = match HgvsRsNormalizer::new(config) {
                Ok(n) => n,
                Err(e) => {
                    // Return error results for all patterns in this chunk
                    return chunk
                        .iter()
                        .enumerate()
                        .map(|(i, pattern)| {
                            let global_idx = chunk_idx * chunk_size + i;
                            (
                                global_idx,
                                crate::benchmark::ParseResult {
                                    input: pattern.clone(),
                                    success: false,
                                    output: None,
                                    error: Some(format!("Failed to create normalizer: {}", e)),
                                    error_category: Some("CONNECTION_ERROR".to_string()),
                                    ref_mismatch: None,
                                    details: None,
                                },
                            )
                        })
                        .collect();
                }
            };

            // Process all patterns in this chunk with the same normalizer
            chunk
                .iter()
                .enumerate()
                .map(|(i, pattern)| {
                    let global_idx = chunk_idx * chunk_size + i;

                    // Translate LRG patterns to RefSeq if mapping available
                    let (translated_pattern, original_lrg) = if let Some(ref mapping) = lrg_mapping
                    {
                        mapping.translate_pattern(pattern)
                    } else {
                        (pattern.clone(), None)
                    };

                    let hgvs_result = normalizer.normalize(&translated_pattern);

                    // If we translated an LRG pattern, report with original accession
                    let output = if original_lrg.is_some() && hgvs_result.success {
                        // Replace RefSeq back to LRG in output for consistency
                        hgvs_result.output.map(|o| {
                            if let Some(ref lrg) = original_lrg {
                                if let Some(colon_pos) = o.find(':') {
                                    format!("{}{}", lrg, &o[colon_pos..])
                                } else {
                                    o
                                }
                            } else {
                                o
                            }
                        })
                    } else {
                        hgvs_result.output
                    };

                    (
                        global_idx,
                        crate::benchmark::ParseResult {
                            input: pattern.clone(), // Always report original input
                            success: hgvs_result.success,
                            output,
                            error: hgvs_result.error.clone(),
                            error_category: hgvs_result.error.map(|e| categorize_hgvs_rs_error(&e)),
                            ref_mismatch: None,
                            details: None,
                        },
                    )
                })
                .collect()
        })
        .collect();

    let elapsed = start.elapsed();

    // Flatten and sort results by original index to preserve ordering
    let mut indexed_results: Vec<(usize, crate::benchmark::ParseResult)> =
        chunk_results.into_iter().flatten().collect();
    indexed_results.sort_by_key(|(idx, _)| *idx);
    let results: Vec<crate::benchmark::ParseResult> =
        indexed_results.into_iter().map(|(_, r)| r).collect();

    // Aggregate error counts
    let mut error_counts: HashMap<String, usize> = HashMap::new();
    for result in &results {
        if let Some(ref category) = result.error_category {
            *error_counts.entry(category.clone()).or_insert(0) += 1;
        }
    }

    Ok((results, elapsed, error_counts))
}

/// Categorize hgvs-rs errors into standard categories.
fn categorize_hgvs_rs_error(error: &str) -> String {
    if error.contains("Parse error") {
        "PARSE_ERROR".to_string()
    } else if error.contains("panic") || error.contains("Panic") {
        "PANIC".to_string()
    } else if error.contains("validation") || error.contains("Validation") {
        "VALIDATION_ERROR".to_string()
    } else if error.contains("connection") || error.contains("Connection") {
        "CONNECTION_ERROR".to_string()
    } else if error.contains("SeqRepo") || error.contains("sequence") {
        "SEQREPO_ERROR".to_string()
    } else if error.contains("transcript") || error.contains("Transcript") {
        "TRANSCRIPT_ERROR".to_string()
    } else if error.contains("intronic") || error.contains("Intronic") {
        "INTRONIC_ERROR".to_string()
    } else if error.contains("protein") || error.contains("Protein") {
        "PROTEIN_ERROR".to_string()
    } else {
        "OTHER".to_string()
    }
}

/// Check if the SeqRepo path exists and appears valid.
pub fn check_seqrepo_path(path: &Path) -> bool {
    if !path.exists() {
        return false;
    }

    // Check for expected SeqRepo structure
    let aliases_db = path.join("aliases.sqlite3");
    let sequences_dir = path.join("sequences");

    aliases_db.exists() && sequences_dir.exists()
}

/// Parse a UTA database URL into components.
pub fn parse_uta_url(url: &str) -> Option<(String, u16, String, String)> {
    // Format: postgresql://user:pass@host:port/dbname
    let url = url.strip_prefix("postgresql://")?;

    // Split at @ to separate credentials from host
    let parts: Vec<&str> = url.splitn(2, '@').collect();
    if parts.len() != 2 {
        return None;
    }

    // Parse host:port/dbname
    let host_db: Vec<&str> = parts[1].splitn(2, '/').collect();
    if host_db.len() != 2 {
        return None;
    }

    let host_port: Vec<&str> = host_db[0].splitn(2, ':').collect();
    let host = host_port[0].to_string();
    let port: u16 = host_port
        .get(1)
        .and_then(|p| p.parse().ok())
        .unwrap_or(5432);
    let dbname = host_db[1].to_string();

    // Parse user:pass
    let creds: Vec<&str> = parts[0].splitn(2, ':').collect();
    let user = creds[0].to_string();

    Some((host, port, dbname, user))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_uta_url() {
        let url = "postgresql://anonymous:anonymous@localhost:5432/uta";
        let result = parse_uta_url(url);
        assert!(result.is_some());

        let (host, port, dbname, user) = result.unwrap();
        assert_eq!(host, "localhost");
        assert_eq!(port, 5432);
        assert_eq!(dbname, "uta");
        assert_eq!(user, "anonymous");
    }

    #[test]
    fn test_parse_uta_url_remote() {
        let url = "postgresql://anonymous:anonymous@uta.biocommons.org:5432/uta";
        let result = parse_uta_url(url);
        assert!(result.is_some());

        let (host, port, dbname, user) = result.unwrap();
        assert_eq!(host, "uta.biocommons.org");
        assert_eq!(port, 5432);
        assert_eq!(dbname, "uta");
        assert_eq!(user, "anonymous");
    }

    #[test]
    fn test_default_config() {
        let config = HgvsRsConfig::default();
        assert!(config.uta_db_url.contains("localhost"));
        assert_eq!(config.uta_db_schema, "uta_20210129b");
    }
}
