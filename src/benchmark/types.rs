//! Types and configuration for comparison benchmarks.

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Duration;

/// Configuration for running a comparison benchmark.
#[derive(Debug, Clone)]
pub struct ComparisonConfig {
    /// Number of parallel workers (shards)
    pub cores: usize,

    /// Output directory for results
    pub output_dir: PathBuf,

    /// Input datasets to process
    pub datasets: Vec<DatasetConfig>,

    /// Whether to run parsing benchmarks
    pub run_parsing: bool,

    /// Whether to run normalization benchmarks
    pub run_normalization: bool,

    /// Whether to include Mutalyzer comparison
    pub include_mutalyzer: bool,

    /// Mutalyzer API URL (for normalization)
    pub mutalyzer_api_url: Option<String>,

    /// Sample size for normalization (full dataset is too slow for Mutalyzer API)
    pub normalization_sample_size: usize,

    /// Path to reference JSON for ferro-hgvs normalization
    pub reference_path: Option<PathBuf>,

    /// Whether to skip steps if outputs already exist
    pub skip_existing: bool,
}

impl Default for ComparisonConfig {
    fn default() -> Self {
        Self {
            cores: 12,
            output_dir: PathBuf::from("benchmark_results"),
            datasets: Vec::new(),
            run_parsing: true,
            run_normalization: true,
            include_mutalyzer: false,
            mutalyzer_api_url: Some("http://localhost:8082".to_string()),
            normalization_sample_size: 10_000,
            reference_path: None,
            skip_existing: true,
        }
    }
}

/// Configuration for a single dataset.
#[derive(Debug, Clone)]
pub struct DatasetConfig {
    /// Name of the dataset (used for file naming)
    pub name: String,

    /// Path to source file
    pub source: PathBuf,

    /// Format of the source file
    pub format: DatasetFormat,

    /// Description for reports
    pub description: String,
}

/// Supported dataset formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatasetFormat {
    /// ClinVar TSV format (hgvs4variation.txt.gz)
    ClinvarTsv,

    /// JSON with test_cases array containing objects with "input" field
    TestCasesJson,

    /// Plain text, one HGVS per line
    PlainText,

    /// JSON array of strings or objects with configurable field
    JsonArray,
}

/// Information about a reference sequence mismatch that was corrected.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RefMismatchInfo {
    /// What the input claimed as reference
    pub stated_ref: String,
    /// What the actual reference sequence has
    pub actual_ref: String,
    /// Position info
    pub position: String,
    /// Whether the mismatch was auto-corrected
    pub corrected: bool,
}

/// Result of parsing a single pattern.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParseResult {
    /// Input HGVS string
    pub input: String,

    /// Whether parsing succeeded
    pub success: bool,

    /// Parsed/normalized output (if successful)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub output: Option<String>,

    /// Error message (if failed)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,

    /// Error category for analysis
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error_category: Option<String>,

    /// Reference mismatch info (if ferro corrected a bad reference)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_mismatch: Option<RefMismatchInfo>,

    /// Parsed variant details (coordinate system, type, position, etc.)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub details: Option<ParsedVariantDetails>,
}

/// Detailed breakdown of a parsed HGVS variant
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParsedVariantDetails {
    /// Reference accession (e.g., "NM_000088.3")
    pub reference: String,
    /// Coordinate system (g, c, n, r, p, m, o)
    pub coordinate_system: String,
    /// Variant type (substitution, deletion, insertion, delins, duplication, inversion, repeat)
    pub variant_type: String,
    /// Position information
    pub position: PositionDetails,
    /// Reference/deleted sequence (if applicable)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub deleted: Option<String>,
    /// Inserted/alternate sequence (if applicable)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inserted: Option<String>,
    /// Whether the variant was shifted during normalization
    #[serde(skip_serializing_if = "Option::is_none")]
    pub was_shifted: Option<bool>,
    /// Original position before shifting (if shifted)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub original_position: Option<String>,
}

/// Position details for a variant
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PositionDetails {
    /// Start position (numeric part)
    pub start: i64,
    /// End position (for ranges/intervals)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end: Option<i64>,
    /// Intronic offset (e.g., +5 in c.123+5A>G)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offset: Option<i32>,
    /// Full position string (e.g., "123+5", "100_102")
    pub display: String,
}

/// Timing information for a benchmark run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingInfo {
    /// Tool name
    pub tool: String,

    /// Total patterns processed
    pub total_patterns: usize,

    /// Patterns successfully processed
    pub successful: usize,

    /// Patterns that failed
    pub failed: usize,

    /// Total elapsed time in seconds
    pub elapsed_seconds: f64,

    /// Throughput (patterns per second)
    pub patterns_per_second: f64,

    /// Average time per pattern in nanoseconds
    pub avg_ns_per_pattern: f64,
}

impl TimingInfo {
    pub fn new(tool: &str, total: usize, successful: usize, elapsed: Duration) -> Self {
        let elapsed_secs = elapsed.as_secs_f64();
        // Use epsilon to avoid division by very small values that would produce
        // misleadingly large or infinite throughput numbers
        let throughput = if elapsed_secs > f64::EPSILON {
            total as f64 / elapsed_secs
        } else {
            0.0
        };
        let avg_ns = if total > 0 {
            elapsed.as_nanos() as f64 / total as f64
        } else {
            0.0
        };

        Self {
            tool: tool.to_string(),
            total_patterns: total,
            successful,
            failed: total - successful,
            elapsed_seconds: elapsed_secs,
            patterns_per_second: throughput,
            avg_ns_per_pattern: avg_ns,
        }
    }
}

/// Results from a shard parsing run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShardResults {
    /// Shard index
    pub shard_index: usize,

    /// Tool name
    pub tool: String,

    /// Input file path
    pub input_file: String,

    /// Timing information
    pub timing: TimingInfo,

    /// Sample of results (limited to save space)
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub sample_results: Vec<ParseResult>,

    /// Sample of failed patterns for analysis
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub failed_examples: Vec<ParseResult>,
}

/// Comparison of parsing results between tools.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParsingComparison {
    /// Dataset name
    pub dataset: String,

    /// ferro-hgvs results
    pub ferro_hgvs: AggregatedResults,

    /// Mutalyzer results (if available)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mutalyzer: Option<AggregatedResults>,

    /// Speedup (ferro throughput / mutalyzer throughput)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub speedup: Option<f64>,

    /// Agreement analysis
    #[serde(skip_serializing_if = "Option::is_none")]
    pub agreement: Option<AgreementStats>,
}

/// Aggregated results across all shards.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregatedResults {
    /// Total patterns processed
    pub total_patterns: usize,

    /// Successfully processed
    pub successful: usize,

    /// Failed
    pub failed: usize,

    /// Pass rate (0.0 - 1.0)
    pub pass_rate: f64,

    /// Total time across all shards
    pub total_time_seconds: f64,

    /// Overall throughput
    pub throughput: f64,
}

impl AggregatedResults {
    pub fn from_timings(timings: &[TimingInfo]) -> Self {
        let total: usize = timings.iter().map(|t| t.total_patterns).sum();
        let successful: usize = timings.iter().map(|t| t.successful).sum();
        let total_time: f64 = timings.iter().map(|t| t.elapsed_seconds).sum();

        Self {
            total_patterns: total,
            successful,
            failed: total - successful,
            pass_rate: if total > 0 {
                successful as f64 / total as f64
            } else {
                0.0
            },
            total_time_seconds: total_time,
            throughput: if total_time > 0.0 {
                total as f64 / total_time
            } else {
                0.0
            },
        }
    }
}

/// Agreement statistics between two tools.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgreementStats {
    /// Both tools succeeded
    pub both_success: usize,

    /// Both tools failed
    pub both_fail: usize,

    /// Only ferro-hgvs succeeded
    pub ferro_only_success: usize,

    /// Only Mutalyzer succeeded
    pub mutalyzer_only_success: usize,

    /// Among patterns where both succeeded, how many produced identical output
    pub agreements: usize,

    /// Agreement rate among both-success cases
    pub agreement_rate: f64,

    /// Sample of disagreements for analysis
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub disagreement_examples: Vec<DisagreementExample>,
}

/// Example of a disagreement between tools.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DisagreementExample {
    pub input: String,
    pub ferro_output: String,
    pub mutalyzer_output: String,
}

/// Statistics about reference sequence mismatches.
///
/// Tracks how many patterns had reference mismatches that ferro corrected,
/// and calculates both strict and lenient agreement rates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceStats {
    /// Number of patterns with reference mismatches that ferro corrected
    pub corrections_count: usize,

    /// Agreement rate if ferro used strict mode (rejected mismatches)
    pub strict_agreement_rate: f64,

    /// Agreement rate with ferro in lenient mode (corrected mismatches)
    pub lenient_agreement_rate: f64,

    /// Patterns that had reference mismatches (limited sample)
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub patterns_with_corrections: Vec<String>,
}

/// Overall benchmark summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonSummary {
    /// Timestamp when the benchmark was generated
    pub generated: DateTime<Utc>,

    /// Configuration used
    pub config: SummaryConfig,

    /// Parsing results by dataset
    pub parsing: HashMap<String, ParsingComparison>,

    /// Normalization results by dataset
    pub normalization: HashMap<String, NormalizationComparison>,

    /// Aggregate statistics
    pub aggregate: AggregateStats,
}

/// Simplified config for inclusion in summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SummaryConfig {
    pub cores: usize,
    pub include_mutalyzer: bool,
    pub normalization_sample_size: usize,
}

/// Aggregate statistics across all datasets.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregateStats {
    pub total_patterns: usize,
    pub ferro_throughput: f64,
    pub mutalyzer_throughput: Option<f64>,
    pub speedup: Option<f64>,
}

/// Normalization comparison results.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NormalizationComparison {
    pub dataset: String,
    pub sample_size: usize,
    pub ferro_hgvs: AggregatedResults,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mutalyzer: Option<AggregatedResults>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub agreement: Option<AgreementStats>,
}

/// Comparison mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CompareMode {
    /// Normalize variants
    Normalize,
    /// Parse only (no normalization)
    Parse,
}

/// External validator to compare against.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum Validator {
    /// Mutalyzer - Independent Python implementation with algebra-based normalizer
    #[default]
    Mutalyzer,
    /// biocommons/hgvs - Canonical Python HGVS implementation
    Biocommons,
    /// hgvs-rs - Rust port of biocommons/hgvs
    HgvsRs,
}

impl std::fmt::Display for Validator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Validator::Mutalyzer => write!(f, "mutalyzer"),
            Validator::Biocommons => write!(f, "biocommons"),
            Validator::HgvsRs => write!(f, "hgvs-rs"),
        }
    }
}

impl std::str::FromStr for Validator {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "mutalyzer" => Ok(Validator::Mutalyzer),
            "biocommons" | "biocommons-hgvs" | "hgvs" => Ok(Validator::Biocommons),
            "hgvs-rs" | "hgvsrs" => Ok(Validator::HgvsRs),
            _ => Err(format!(
                "Unknown validator: {}. Valid options: mutalyzer, biocommons, hgvs-rs",
                s
            )),
        }
    }
}

impl std::fmt::Display for CompareMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CompareMode::Normalize => write!(f, "normalize"),
            CompareMode::Parse => write!(f, "parse"),
        }
    }
}

/// Result of a comparison between ferro-hgvs and mutalyzer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonResult {
    /// Comparison mode (normalize or parse)
    pub mode: CompareMode,

    /// Timestamp when comparison was run
    pub timestamp: DateTime<Utc>,

    /// Number of patterns compared
    pub sample_size: usize,

    /// ferro-hgvs timing and results
    pub ferro: ComparisonToolResult,

    /// Mutalyzer timing and results
    pub mutalyzer: ComparisonToolResult,

    /// Speedup factor (mutalyzer_time / ferro_time)
    pub speedup: f64,

    /// Agreement statistics
    pub agreement: AgreementStats,

    /// Examples of differences (limited sample)
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub differences: Vec<DisagreementExample>,

    /// Reference mismatch statistics
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reference_stats: Option<ReferenceStats>,

    /// Cache coverage statistics
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cache_stats: Option<CacheStats>,
}

/// Result from a single tool in a comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonToolResult {
    /// Total elapsed time in seconds
    pub elapsed_seconds: f64,

    /// Number of successfully processed patterns
    pub successful: usize,

    /// Number of failed patterns
    pub failed: usize,

    /// Throughput (patterns per second)
    pub throughput: f64,

    /// Error counts by category (for failed patterns)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error_counts: Option<HashMap<String, usize>>,
}

/// Statistics about cache coverage for a comparison run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheStats {
    /// Total unique accessions needed
    pub total_accessions: usize,

    /// Accessions missing from ferro cache
    pub ferro_missing: usize,

    /// Accessions missing from mutalyzer cache
    pub mutalyzer_missing: usize,

    /// Sample of missing accessions (limited)
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub ferro_missing_examples: Vec<String>,

    /// Sample of missing accessions (limited)
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub mutalyzer_missing_examples: Vec<String>,
}

impl ComparisonToolResult {
    pub fn new(total: usize, successful: usize, elapsed: std::time::Duration) -> Self {
        let elapsed_secs = elapsed.as_secs_f64();
        let throughput = if elapsed_secs > 0.0 {
            total as f64 / elapsed_secs
        } else {
            0.0
        };
        Self {
            elapsed_seconds: elapsed_secs,
            successful,
            failed: total - successful,
            throughput,
            error_counts: None,
        }
    }

    pub fn with_error_counts(mut self, counts: HashMap<String, usize>) -> Self {
        self.error_counts = Some(counts);
        self
    }
}

/// Pattern category for stratified sampling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PatternCategory {
    GenomicSnv,
    GenomicDel,
    GenomicIns,
    GenomicDup,
    GenomicDelins,
    GenomicInv,
    CodingSnv,
    CodingIntronic,
    CodingDel,
    CodingIns,
    CodingDup,
    ProteinSub,
    ProteinFs,
    ProteinExt,
    NonCoding,
    Mitochondrial,
    Rna,
    Circular,
    Repeat,
    Uncertain,
    Other,
}

impl PatternCategory {
    /// Check if this category is a protein category.
    ///
    /// Protein patterns (NP_*, :p.) are handled inconsistently by normalization tools:
    /// - ferro/hgvs-rs: return unchanged (no-op) as success
    /// - biocommons: explicitly rejects as "unsupported"
    /// - mutalyzer: requires network for NP_→NM_ lookup
    ///
    /// For fair normalization comparisons, protein patterns should be excluded.
    pub fn is_protein(&self) -> bool {
        matches!(
            self,
            PatternCategory::ProteinSub | PatternCategory::ProteinFs | PatternCategory::ProteinExt
        )
    }

    /// Categorize an HGVS pattern string.
    pub fn categorize(pattern: &str) -> Self {
        // Check coordinate type first
        if pattern.contains(":g.") {
            if pattern.ends_with('>') || pattern.chars().filter(|&c| c == '>').count() == 1 {
                return Self::GenomicSnv;
            }
            if pattern.contains("del") && pattern.contains("ins") {
                return Self::GenomicDelins;
            }
            if pattern.contains("del") {
                return Self::GenomicDel;
            }
            if pattern.contains("ins") {
                return Self::GenomicIns;
            }
            if pattern.contains("dup") {
                return Self::GenomicDup;
            }
            if pattern.contains("inv") {
                return Self::GenomicInv;
            }
            return Self::Other;
        }

        if pattern.contains(":c.") {
            // Check for intronic
            if pattern.contains('+') || pattern.contains('-') {
                // Make sure it's actually intronic offset, not a range
                let after_dot = pattern.split(":c.").nth(1).unwrap_or("");
                if after_dot.contains('+') || after_dot.chars().any(|c| c == '-') {
                    return Self::CodingIntronic;
                }
            }
            if pattern.ends_with('>') || pattern.chars().filter(|&c| c == '>').count() == 1 {
                return Self::CodingSnv;
            }
            if pattern.contains("del") {
                return Self::CodingDel;
            }
            if pattern.contains("ins") {
                return Self::CodingIns;
            }
            if pattern.contains("dup") {
                return Self::CodingDup;
            }
            return Self::Other;
        }

        if pattern.contains(":p.") {
            if pattern.contains("fs") {
                return Self::ProteinFs;
            }
            if pattern.contains("ext") {
                return Self::ProteinExt;
            }
            return Self::ProteinSub;
        }

        if pattern.contains(":n.") {
            return Self::NonCoding;
        }

        if pattern.contains(":m.") {
            return Self::Mitochondrial;
        }

        if pattern.contains(":r.") {
            return Self::Rna;
        }

        if pattern.contains(":o.") {
            return Self::Circular;
        }

        if pattern.contains('[') && pattern.contains(']') {
            return Self::Repeat;
        }

        if pattern.contains("(?)") || pattern.contains('?') {
            return Self::Uncertain;
        }

        Self::Other
    }
}

/// Combined result for a single pattern from both ferro and mutalyzer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CombinedPatternResult {
    /// Input HGVS string
    pub input: String,

    /// Whether ferro succeeded
    pub ferro_success: bool,

    /// Ferro output (if successful)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ferro_output: Option<String>,

    /// Ferro error (if failed)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ferro_error: Option<String>,

    /// Whether mutalyzer succeeded
    pub mutalyzer_success: bool,

    /// Mutalyzer output (if successful)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mutalyzer_output: Option<String>,

    /// Mutalyzer error (if failed)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mutalyzer_error: Option<String>,

    /// Whether outputs match (None if either failed)
    #[serde(rename = "match")]
    pub outputs_match: Option<bool>,
}

/// Detailed per-pattern results for both ferro and mutalyzer.
///
/// This is saved separately from the summary comparison result to allow
/// detailed analysis of individual patterns and their outcomes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetailedResults {
    /// Comparison mode (normalize or parse)
    pub mode: CompareMode,

    /// Timestamp when comparison was run
    pub timestamp: DateTime<Utc>,

    /// Number of patterns compared
    pub sample_size: usize,

    /// Combined per-pattern results
    pub results: Vec<CombinedPatternResult>,
}

// =============================================================================
// Mutalyzer Benchmark Types
// =============================================================================

/// Configuration for running a mutalyzer bulk benchmark.
#[derive(Debug, Clone)]
pub struct MutalyzerBenchmarkConfig {
    /// Path to input patterns file
    pub input_path: PathBuf,

    /// Path to output JSONL.gz file
    pub output_path: PathBuf,

    /// Number of parallel workers
    pub workers: usize,

    /// Patterns per batch for progress reporting
    pub batch_size: usize,

    /// Path to mutalyzer settings file
    pub settings_file: Option<PathBuf>,

    /// Whether to allow network access
    pub allow_network: bool,

    /// Path to existing results file (for incremental processing)
    pub existing_path: Option<PathBuf>,

    /// Progress reporting interval in seconds
    pub progress_interval: u64,

    /// Whether to skip failed patterns (don't retry them)
    pub skip_failed: bool,
}

impl Default for MutalyzerBenchmarkConfig {
    fn default() -> Self {
        Self {
            input_path: PathBuf::new(),
            output_path: PathBuf::new(),
            workers: 24,
            batch_size: 10_000,
            settings_file: None,
            allow_network: false,
            existing_path: None,
            progress_interval: 30,
            skip_failed: false,
        }
    }
}

/// Existing results loaded from a previous benchmark run.
///
/// Used to skip patterns that already succeeded and retry failed patterns.
/// Stores full ParseResult for successful patterns to preserve the output.
#[derive(Debug, Clone)]
pub struct ExistingResults {
    /// Patterns that succeeded (full result with output)
    pub successful: HashMap<String, ParseResult>,

    /// Patterns that failed (full result with error info)
    pub failed: HashMap<String, ParseResult>,
}

impl ExistingResults {
    /// Create an empty ExistingResults
    pub fn new() -> Self {
        Self {
            successful: HashMap::new(),
            failed: HashMap::new(),
        }
    }

    /// Check if a pattern should be processed.
    ///
    /// Returns true if the pattern failed previously or is new.
    /// Returns false if the pattern already succeeded.
    pub fn should_process(&self, pattern: &str) -> bool {
        !self.successful.contains_key(pattern)
    }

    /// Number of successful patterns
    pub fn successful_count(&self) -> usize {
        self.successful.len()
    }

    /// Number of failed patterns
    pub fn failed_count(&self) -> usize {
        self.failed.len()
    }
}

impl Default for ExistingResults {
    fn default() -> Self {
        Self::new()
    }
}

/// Metadata for a benchmark run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkMetadata {
    /// Start time of the benchmark
    pub start_time: DateTime<Utc>,

    /// End time of the benchmark (if complete)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub end_time: Option<DateTime<Utc>>,

    /// Total patterns in input
    pub total_patterns: usize,

    /// Number of workers used
    pub workers: usize,

    /// Batch size used
    pub batch_size: usize,

    /// Whether network access was allowed
    pub allow_network: bool,

    /// Path to existing results file (if incremental)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub existing_file: Option<String>,
}

/// Summary of a benchmark run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkSummary {
    /// Metadata about the run
    pub metadata: BenchmarkMetadata,

    /// Number of patterns processed in this run
    pub processed: usize,

    /// Number of successful patterns (new + from existing)
    pub successful: usize,

    /// Number of failed patterns
    pub failed: usize,

    /// Total elapsed time in seconds
    pub elapsed_seconds: f64,

    /// Throughput (patterns per second)
    pub throughput: f64,

    /// Error counts by category
    pub error_counts: HashMap<String, usize>,

    /// Statistics from existing results (if incremental)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub existing_stats: Option<ExistingStats>,
}

/// Statistics about existing results used in incremental processing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExistingStats {
    /// Number of patterns skipped (already successful)
    pub skipped: usize,

    /// Number of patterns retried (previously failed)
    pub retried: usize,

    /// Number of retried patterns that succeeded
    pub retry_successes: usize,
}

/// Unified output format for parse tool results.
///
/// All parse tools (ferro, mutalyzer, biocommons, hgvs-rs) should output
/// this format to enable consistent comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolParseOutput {
    /// Tool name (e.g., "ferro", "mutalyzer", "biocommons", "hgvs-rs")
    pub tool: String,

    /// Total patterns processed
    pub total_patterns: usize,

    /// Patterns successfully parsed
    pub successful: usize,

    /// Patterns that failed to parse
    pub failed: usize,

    /// Total elapsed time in seconds
    pub elapsed_seconds: f64,

    /// Throughput (patterns per second)
    pub throughput: f64,

    /// All parse results
    pub results: Vec<ParseResult>,
}

impl ToolParseOutput {
    /// Create a new ToolParseOutput from results
    pub fn new(tool: &str, results: Vec<ParseResult>, elapsed: std::time::Duration) -> Self {
        let total = results.len();
        let successful = results.iter().filter(|r| r.success).count();
        let elapsed_secs = elapsed.as_secs_f64();
        let throughput = if elapsed_secs > 0.0 {
            total as f64 / elapsed_secs
        } else {
            0.0
        };

        Self {
            tool: tool.to_string(),
            total_patterns: total,
            successful,
            failed: total - successful,
            elapsed_seconds: elapsed_secs,
            throughput,
            results,
        }
    }
}
