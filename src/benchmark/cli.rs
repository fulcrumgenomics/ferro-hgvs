// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! Clap CLI definitions for the `ferro-benchmark` binary.
//!
//! Extracted from `src/bin/benchmark.rs` so integration tests can exercise the
//! parser in-process via [`Cli::try_parse_from`] / [`clap::CommandFactory::command`]
//! without spawning the binary as a subprocess.

use clap::{ArgAction, Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

/// Tool to use for parsing/normalization
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Tool {
    /// ferro-hgvs (native Rust)
    Ferro,
    /// mutalyzer (Python)
    Mutalyzer,
    /// biocommons/hgvs (Python)
    Biocommons,
    /// hgvs-rs (native Rust, requires hgvs-rs feature)
    HgvsRs,
    /// All tools (ferro → mutalyzer → biocommons → hgvs-rs)
    All,
}

impl std::fmt::Display for Tool {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Tool::Ferro => write!(f, "ferro"),
            Tool::Mutalyzer => write!(f, "mutalyzer"),
            Tool::Biocommons => write!(f, "biocommons"),
            Tool::HgvsRs => write!(f, "hgvs-rs"),
            Tool::All => write!(f, "all"),
        }
    }
}

/// Single-tool variant for subcommands (`parse`, `normalize`) where running
/// against every tool isn't supported. Mirrors `Tool` minus the `All`
/// variant so clap rejects `all` at parse time instead of letting the
/// command fall through to a runtime error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ToolNonAll {
    /// ferro-hgvs (native Rust)
    Ferro,
    /// mutalyzer (Python)
    Mutalyzer,
    /// biocommons/hgvs (Python)
    Biocommons,
    /// hgvs-rs (native Rust, requires hgvs-rs feature)
    HgvsRs,
}

impl From<ToolNonAll> for Tool {
    fn from(t: ToolNonAll) -> Self {
        match t {
            ToolNonAll::Ferro => Tool::Ferro,
            ToolNonAll::Mutalyzer => Tool::Mutalyzer,
            ToolNonAll::Biocommons => Tool::Biocommons,
            ToolNonAll::HgvsRs => Tool::HgvsRs,
        }
    }
}

/// Genome assembly option for prepare ferro
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GenomeOption {
    /// Don't download any genome
    None,
    /// Download GRCh37 only
    Grch37,
    /// Download GRCh38 only
    Grch38,
    /// Download both GRCh37 and GRCh38
    All,
}

pub fn parse_genome_option(s: &str) -> Result<GenomeOption, String> {
    match s.to_lowercase().as_str() {
        "none" => Ok(GenomeOption::None),
        "grch37" => Ok(GenomeOption::Grch37),
        "grch38" => Ok(GenomeOption::Grch38),
        "all" => Ok(GenomeOption::All),
        _ => Err(format!(
            "Invalid genome option '{}'. Valid options: none, grch37, grch38, all",
            s
        )),
    }
}

#[derive(Parser)]
#[command(name = "ferro-benchmark")]
#[command(author, version, about = "Benchmark ferro-hgvs against other HGVS tools (mutalyzer, biocommons, hgvs-rs)", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
pub enum Commands {
    // Run command removed - use individual commands instead

    // OLD FLAT COMMANDS REMOVED - now in subcommand groups:
    // - ExtractClinvar, ExtractJson, Shard, Sample → extract subcommand
    // - CollateParsing, CollateNormalization → collate subcommand
    // - GenerateSummary, GenerateReport, GenerateReadmeTables → generate subcommand
    // - CompareNormalize, CompareParse → compare subcommand
    // - StartUta, StopUta → setup subcommand (Uta, Seqrepo integrated into prepare)
    // - GenerateBiocommonsSettings → generate subcommand
    // - BuildTranscriptMapping removed - integrated into prepare mutalyzer
    // PopulateProteinCache removed - use --clinvar flag with prepare mutalyzer

    // =========================================================================
    // NEW UNIFIED CLI COMMANDS
    // =========================================================================
    /// Prepare reference data for a tool
    ///
    /// Downloads and sets up the necessary reference data for the specified tool.
    /// Each tool has different data requirements:
    /// - ferro: RefSeq transcripts + cdot metadata
    /// - mutalyzer: Cache built from ferro data (requires ferro first)
    /// - biocommons: UTA database + SeqRepo (shared with hgvs-rs)
    /// - hgvs-rs: UTA database + SeqRepo (shared with biocommons)
    #[command(name = "prepare")]
    PrepareTool {
        /// Tool to prepare data for
        #[arg(value_enum)]
        tool: Tool,

        /// Output directory for reference data
        #[arg(short, long, default_value = "benchmark-output")]
        output_dir: PathBuf,

        /// Genome assembly to download (ferro only)
        ///
        /// Options: none, grch37, grch38, all (case-insensitive)
        /// Default: all (~4GB total for both assemblies)
        #[arg(long, default_value = "all", value_parser = parse_genome_option)]
        genome: GenomeOption,

        /// Skip RefSeqGene sequences (ferro only, ~600MB)
        ///
        /// By default, RefSeqGene sequences are downloaded.
        #[arg(long, action = ArgAction::SetTrue)]
        no_refseqgene: bool,

        /// Skip LRG sequences (ferro: ~50MB download, mutalyzer: LRG annotation enhancement)
        ///
        /// By default, LRG data is downloaded/enhanced.
        #[arg(long, action = ArgAction::SetTrue)]
        no_lrg: bool,

        /// Force re-download even if files exist
        #[arg(long, action = ArgAction::SetTrue)]
        force: bool,

        /// Ferro reference directory (required for mutalyzer)
        #[arg(long)]
        ferro_reference: Option<PathBuf>,

        /// Pattern file for accession extraction (ferro stores for auto-detection)
        ///
        /// For ferro: Stores patterns file in reference for provenance. Other tools
        /// auto-detect it from {ferro-reference}/patterns/.
        /// For mutalyzer/biocommons/hgvs-rs: Override auto-detected patterns.
        #[arg(long)]
        patterns: Option<PathBuf>,

        /// ClinVar hgvs4variation.txt.gz for protein cache (ferro stores for auto-detection)
        ///
        /// For ferro: Symlinks ClinVar file to reference for provenance. Mutalyzer
        /// auto-detects it from {ferro-reference}/clinvar/.
        /// For mutalyzer: Override auto-detected clinvar file.
        #[arg(long)]
        clinvar: Option<PathBuf>,

        /// Skip protein sequence cache population (mutalyzer only)
        ///
        /// By default, if --clinvar is provided, protein sequences are fetched.
        /// Use this flag to skip protein cache population.
        #[arg(long, action = ArgAction::SetTrue)]
        no_proteins: bool,

        /// Number of cache shards for parallel workers (mutalyzer only)
        ///
        /// Partitions the mutalyzer cache into N shards so workers can load
        /// smaller subsets. Use 32 for up to 32 workers, reducing memory per
        /// worker from ~18GB to ~560MB. Patterns are routed to shards by
        /// hashing their accession.
        #[arg(long, default_value_t = 32)]
        shards: usize,

        /// Port for UTA PostgreSQL (biocommons/hgvs-rs only)
        #[arg(long, default_value_t = 5432)]
        uta_port: u16,

        /// SeqRepo directory path (biocommons/hgvs-rs only)
        #[arg(long)]
        seqrepo_dir: Option<PathBuf>,

        /// UTA dump file path (biocommons/hgvs-rs only)
        ///
        /// If UTA is not running and this path is provided, automatically set up UTA.
        /// Download from: https://dl.biocommons.org/uta/uta_20210129b.pgd.gz
        #[arg(long)]
        uta_dump: Option<PathBuf>,

        /// Skip loading ferro transcripts into SeqRepo (biocommons/hgvs-rs only)
        ///
        /// By default, ALL transcripts from ferro's transcripts/ directory (~270K sequences)
        /// are loaded into SeqRepo, making biocommons/hgvs-rs fully self-contained.
        /// Use this flag to skip transcript loading (only load supplemental sequences).
        #[arg(long, action = ArgAction::SetTrue)]
        no_load_transcripts: bool,

        /// Skip loading cdot transcript alignments into UTA (biocommons/hgvs-rs only)
        ///
        /// By default, transcript alignments from cdot are loaded into UTA for
        /// transcripts that are missing from the UTA snapshot. This enables hgvs-rs
        /// to normalize variants for newer transcripts. Use this flag to skip loading.
        #[arg(long, action = ArgAction::SetTrue)]
        no_load_alignments: bool,
    },

    /// Check if a tool and its dependencies are properly configured
    ///
    /// Verifies that the specified tool is available and its reference data is set up:
    /// - ferro: Check reference data (transcripts, cdot, optional genome)
    /// - mutalyzer: Check parser and normalizer packages, cache directory
    /// - biocommons: Check hgvs package, UTA connection, SeqRepo
    /// - hgvs-rs: Check SeqRepo path, UTA connection, provider initialization
    #[command(name = "check")]
    CheckTool {
        /// Tool to check
        #[arg(value_enum)]
        tool: Tool,

        /// Reference data directory (ferro)
        #[arg(long, default_value = "benchmark-output")]
        reference: PathBuf,

        /// Mutalyzer cache/settings directory (mutalyzer)
        #[arg(long)]
        mutalyzer_settings: Option<PathBuf>,

        /// UTA database URL (biocommons/hgvs-rs)
        #[arg(long)]
        uta_db_url: Option<String>,

        /// SeqRepo directory path (biocommons/hgvs-rs)
        #[arg(long)]
        seqrepo_path: Option<PathBuf>,
    },

    /// Parse HGVS patterns with a tool
    ///
    /// Parses HGVS patterns using the specified tool and outputs results.
    #[command(name = "parse")]
    ParseTool {
        /// Tool to use for parsing
        #[arg(value_enum)]
        tool: ToolNonAll,

        /// Input patterns file
        #[arg(short, long)]
        input: PathBuf,

        /// Output results JSON file
        #[arg(short, long)]
        output: PathBuf,

        /// Timing output JSON file (optional)
        #[arg(short, long)]
        timing: Option<PathBuf>,
    },

    /// Normalize HGVS patterns with a tool
    ///
    /// Normalizes HGVS patterns using the specified tool and outputs results.
    #[command(name = "normalize")]
    NormalizeTool {
        /// Tool to use for normalization
        #[arg(value_enum)]
        tool: ToolNonAll,

        /// Input patterns file
        #[arg(short, long)]
        input: PathBuf,

        /// Output results JSON file
        #[arg(short, long)]
        output: PathBuf,

        /// Timing output JSON file (optional)
        #[arg(short, long)]
        timing: Option<PathBuf>,

        /// Reference data directory (ferro)
        #[arg(long)]
        reference: Option<PathBuf>,

        /// Mutalyzer settings file (mutalyzer)
        #[arg(long)]
        mutalyzer_settings: Option<PathBuf>,

        /// Biocommons settings file (biocommons)
        #[arg(long)]
        biocommons_settings: Option<PathBuf>,

        /// LRG-to-RefSeq mapping file for LRG transcript translation (biocommons, hgvs-rs)
        #[arg(long)]
        lrg_mapping: Option<PathBuf>,

        /// SeqRepo directory path (hgvs-rs)
        #[arg(long)]
        seqrepo_path: Option<PathBuf>,

        /// Number of parallel workers
        #[arg(short = 'j', long, default_value_t = 1)]
        workers: usize,

        /// Skip intronic variant rewriting for mutalyzer (default: rewrite enabled)
        #[arg(long)]
        no_rewrite_intronic: bool,

        /// Allow network access for benchmarking remote vs local performance
        ///
        /// By default, network access is blocked to ensure reproducible local-only
        /// benchmarking. Enable this flag to measure baseline performance without
        /// ferro's reference preparation, demonstrating the speedup from caching.
        ///
        /// Warning: Network mode is significantly slower and requires internet access.
        #[arg(long)]
        allow_network: bool,

        /// Use in-memory provider for hgvs-rs (pre-loads all data, eliminates DB I/O from timing)
        #[arg(long)]
        in_memory: bool,
    },

    // CompareTool → compare results subcommand

    // =========================================================================
    // GROUPED SUBCOMMANDS (Phase 2)
    // =========================================================================
    /// Compare pre-computed results between tools
    #[command(subcommand)]
    Compare(CompareCommands),

    /// Run live benchmarks (sample, parse/normalize, compare)
    #[command(subcommand)]
    Benchmark(BenchmarkCommands),

    /// Extract and sample data
    #[command(subcommand)]
    Extract(ExtractCommands),

    /// Set up external services
    #[command(subcommand)]
    Setup(SetupCommands),

    /// Generate reports and configs
    #[command(subcommand)]
    Generate(GenerateCommands),

    /// Collate sharded results
    #[command(subcommand)]
    Collate(CollateCommands),
}

// =============================================================================
// SUBCOMMAND ENUMS (Phase 2)
// =============================================================================

/// Compare subcommands (for pre-computed results)
#[derive(Subcommand)]
pub enum CompareCommands {
    /// Compare results from multiple tools
    Results {
        /// Comparison mode (parse or normalize)
        #[arg(value_enum)]
        mode: CompareMode,

        /// Result files to compare (at least 2)
        #[arg(required = true, num_args = 2..)]
        results: Vec<PathBuf>,

        /// Output comparison JSON file
        #[arg(short, long, default_value = "comparison.json")]
        output: PathBuf,

        /// Enable equivalence checking for normalization (semantic comparison)
        #[arg(long, action = ArgAction::SetTrue)]
        equivalence: bool,
    },

    /// Show arbitration decisions
    Arbitration {
        /// Path to arbitration database JSON file
        #[arg(short, long, default_value = "tests/fixtures/arbitration.json")]
        database: PathBuf,

        /// Filter by category
        #[arg(long)]
        category: Option<String>,

        /// Show detailed reasoning
        #[arg(long, action = ArgAction::SetTrue)]
        verbose: bool,
    },
}

/// Arguments for the `benchmark matrix` subcommand.
#[derive(clap::Args, Debug)]
#[allow(clippy::struct_excessive_bools)]
pub struct MatrixArgs {
    /// Population file (one HGVS expression per line, plain text).
    #[arg(long)]
    pub population: PathBuf,

    /// Output results JSON file.
    #[arg(long, default_value = "data/benchmark/perf_results.json")]
    pub output: PathBuf,

    /// Manual sample size override for all (tool, op) pairs.
    ///
    /// When set, calibration is skipped and every tool uses this fixed N.
    /// Useful for smoke runs. When omitted, N is calibrated per (tool, op)
    /// using --target-seconds / --min-sample / --max-sample.
    #[arg(long)]
    pub sample_size: Option<usize>,

    /// Target wall-time per measurement used to calibrate per-(tool,op) N.
    ///
    /// Ignored when --sample-size is set.
    #[arg(long, default_value_t = 3.0)]
    pub target_seconds: f64,

    /// Minimum sample size (floor) for calibrated N.
    ///
    /// Ensures statistically meaningful samples even for very slow tools.
    /// Ignored when --sample-size is set.
    #[arg(long, default_value_t = 50)]
    pub min_sample: usize,

    /// Maximum sample size (ceiling) for calibrated N.
    ///
    /// Bounds memory use for very fast tools.
    /// Ignored when --sample-size is set.
    #[arg(long, default_value_t = 2_000_000)]
    pub max_sample: usize,

    /// Cap on the ferro full-population pass (0 = whole population).
    #[arg(long, default_value_t = 0)]
    pub ferro_full_n: usize,

    /// Repetitions per (tool, worker-count) cell.
    #[arg(long, default_value_t = 5)]
    pub reps: u32,

    /// Base random seed; rep r uses seed + r.
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Worker counts for the cross-tool comparison tables (comma-separated).
    #[arg(long, value_delimiter = ',', default_value = "1,8")]
    pub workers: Vec<usize>,

    /// ferro thread-scaling counts (comma-separated).
    #[arg(long, value_delimiter = ',', default_value = "1,2,4,8")]
    pub ferro_threads: Vec<usize>,

    /// Operations to run (comma-separated: parse, normalize).
    #[arg(long, value_delimiter = ',', default_value = "parse,normalize")]
    pub operations: Vec<String>,

    /// ferro reference data directory (required for normalization).
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// mutalyzer settings file (mutalyzer_settings.conf).
    #[arg(long)]
    pub mutalyzer_settings: Option<PathBuf>,

    /// biocommons settings file (biocommons_settings.conf).
    #[arg(long)]
    pub biocommons_settings: Option<PathBuf>,

    /// SeqRepo data directory path (biocommons / hgvs-rs).
    #[arg(long)]
    pub seqrepo_path: Option<PathBuf>,

    /// UTA database URL including schema suffix
    /// (e.g. postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b).
    #[arg(long)]
    pub uta_db_url: Option<String>,

    /// Exclude protein variants (NP_*, :p.) from samples — recommended for normalize.
    #[arg(long, default_value_t = true)]
    pub exclude_protein: bool,

    /// Human-readable machine label recorded in provenance (e.g. "Apple M2 Max").
    #[arg(long, default_value = "unknown")]
    pub machine: String,
}

/// Benchmark subcommands (live sampling and comparison)
#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
pub enum BenchmarkCommands {
    /// Benchmark parsing between ferro-hgvs and mutalyzer
    Parse {
        /// Input patterns file (plain text, ClinVar TSV, or gzipped)
        #[arg(short, long)]
        input: PathBuf,

        /// Output comparison results JSON file
        #[arg(short, long, default_value = "benchmark_parse.json")]
        output: PathBuf,

        /// Number of patterns to sample (0 = all patterns)
        #[arg(short, long, default_value_t = 1000)]
        sample_size: usize,

        /// Number of parallel workers for mutalyzer
        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        workers: usize,

        /// Random seed for sampling
        #[arg(long, default_value_t = 42)]
        seed: u64,

        /// Error handling mode for ferro: strict, lenient, or silent
        #[arg(long, default_value = "strict", value_parser = ["strict", "lenient", "silent"])]
        error_mode: String,
    },

    /// Run the full performance matrix (tools × operations × worker counts × reps)
    /// over a shared seeded stratified sample and write perf_results.json.
    Matrix(MatrixArgs),

    /// Benchmark normalization between ferro-hgvs and an external validator
    Normalize {
        /// Input patterns file (plain text, ClinVar TSV, or gzipped)
        #[arg(short, long)]
        input: PathBuf,

        /// Output comparison results JSON file
        #[arg(short, long, default_value = "benchmark_normalize.json")]
        output: PathBuf,

        /// External validator to compare against
        #[arg(long, default_value = "mutalyzer")]
        validator: String,

        /// Number of patterns to sample (0 = all patterns)
        #[arg(short, long, default_value_t = 1000)]
        sample_size: usize,

        /// Number of parallel workers for external validator
        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        workers: usize,

        /// Random seed for sampling
        #[arg(long, default_value_t = 42)]
        seed: u64,

        /// Reference data directory for ferro-hgvs
        #[arg(short, long)]
        reference: Option<PathBuf>,

        /// Mutalyzer settings file path
        #[arg(long)]
        mutalyzer_settings: Option<String>,

        /// Transcript→chromosome mapping file
        #[arg(long)]
        transcript_mapping: Option<PathBuf>,

        /// Include unnormalizable patterns
        #[arg(long)]
        include_unnormalizable: bool,

        /// Include genomic (NC_) patterns
        #[arg(long)]
        include_genomic: bool,

        /// Include protein (NP_) patterns
        #[arg(long)]
        include_protein: bool,

        /// Allow mutalyzer network requests
        #[arg(long)]
        allow_mutalyzer_network: bool,

        /// Detailed per-pattern output
        #[arg(long)]
        detailed_output: Option<PathBuf>,

        /// Pre-computed mutalyzer results
        #[arg(long)]
        existing_mutalyzer_results: Option<PathBuf>,

        /// UTA database URL for biocommons
        #[arg(long)]
        uta_db_url: Option<String>,

        /// Biocommons settings file
        #[arg(long)]
        biocommons_settings: Option<PathBuf>,

        /// SeqRepo path for hgvs-rs
        #[arg(long)]
        hgvs_rs_seqrepo_path: Option<String>,

        /// UTA schema for hgvs-rs
        #[arg(long, default_value = "uta_20210129b")]
        hgvs_rs_uta_schema: Option<String>,

        /// Arbitration database file
        #[arg(long)]
        arbitration: Option<PathBuf>,

        /// Error handling mode for ferro: strict, lenient, or silent
        #[arg(long, default_value = "strict", value_parser = ["strict", "lenient", "silent"])]
        error_mode: String,
    },
}

/// Extract subcommands
#[derive(Subcommand)]
pub enum ExtractCommands {
    /// Extract patterns from ClinVar TSV
    Clinvar {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Extract patterns from JSON file
    Json {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value = "test_cases_json")]
        format: String,
    },

    /// Extract protein accessions from ClinVar
    Proteins {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Create stratified sample
    Sample {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(short, long, default_value_t = 1000)]
        size: usize,
        #[arg(long)]
        seed: Option<u64>,
        /// Exclude protein patterns (NP_*, :p.) - recommended for normalize samples
        #[arg(long)]
        exclude_protein: bool,
    },

    /// Shard dataset for parallel processing
    Shard {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output_dir: PathBuf,
        #[arg(short, long, default_value_t = num_cpus::get())]
        num_shards: usize,
    },
}

/// Setup subcommands
#[derive(Subcommand)]
pub enum SetupCommands {
    /// Set up local UTA database via Docker
    Uta {
        #[arg(long, default_value = "ferro-uta")]
        container_name: String,
        #[arg(long, default_value = "uta_20210129b")]
        image_tag: String,
        #[arg(long, default_value_t = 5432)]
        port: u16,
        #[arg(long)]
        uta_dump: Option<PathBuf>,
        #[arg(long, action = ArgAction::SetTrue)]
        force: bool,
        /// Seconds to wait for the UTA schema to be ready (default 300).
        #[arg(long, default_value_t = 300)]
        uta_ready_timeout_secs: u64,
    },

    /// Set up local SeqRepo repository
    Seqrepo {
        #[arg(short, long, default_value = "seqrepo")]
        output_dir: PathBuf,
        #[arg(long, default_value = "2021-01-29")]
        instance: String,
    },

    /// Start local UTA Docker container
    #[command(name = "start-uta")]
    StartUta {
        #[arg(long, default_value = "ferro-uta")]
        container_name: String,
    },

    /// Stop local UTA Docker container
    #[command(name = "stop-uta")]
    StopUta {
        #[arg(long, default_value = "ferro-uta")]
        container_name: String,
        #[arg(long, action = ArgAction::SetTrue)]
        force: bool,
    },
}

/// Generate subcommands
#[derive(Subcommand)]
pub enum GenerateCommands {
    /// Generate summary from comparison files
    Summary {
        #[arg(long)]
        normalization: Vec<PathBuf>,
        #[arg(long)]
        parsing: Vec<PathBuf>,
        #[arg(short, long, default_value = "summary.json")]
        output: PathBuf,
    },

    /// Generate markdown report from summary
    Report {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, default_value = "REPORT.md")]
        output: PathBuf,
    },

    /// Generate README tables from summary
    Readme {
        #[arg(short, long)]
        input: PathBuf,
    },
}

/// Collate subcommands
#[derive(Subcommand)]
pub enum CollateCommands {
    /// Collate parsing results from shards
    Parsing {
        #[arg(long)]
        ferro_dir: PathBuf,
        #[arg(long)]
        mutalyzer_dir: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value = "unknown")]
        dataset: String,
    },

    /// Collate normalization results
    Normalization {
        #[arg(long)]
        ferro_results: PathBuf,
        #[arg(long)]
        mutalyzer_results: Option<PathBuf>,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value = "unknown")]
        dataset: String,
    },
}

/// Comparison mode for compare command
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum CompareMode {
    /// Compare parse results
    Parse,
    /// Compare normalize results
    Normalize,
}
