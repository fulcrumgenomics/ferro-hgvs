// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! ferro-benchmark: CLI for benchmarking ferro-hgvs against other HGVS tools.
//!
//! Supports comparison with: mutalyzer, biocommons/hgvs, and hgvs-rs.
//!
//! Build with: cargo build --release --features benchmark
//! Run with: ferro-benchmark --help

use clap::{ArgAction, Parser, Subcommand, ValueEnum};
use ferro_hgvs::benchmark::{
    biocommons::{
        check_docker_available, check_seqrepo, check_uta_connection,
        fetch_and_load_missing_accessions, has_biocommons_normalizer, load_biocommons_settings,
        load_sequences_to_seqrepo, run_biocommons_normalizer_parallel,
        run_biocommons_normalizer_subprocess, setup_seqrepo, setup_uta, start_uta, stop_uta,
        write_biocommons_settings, BiocommonsLocalConfig,
    },
    cache::{
        build_transcript_chromosome_mapping, extract_all_accessions_from_file, fetch_fasta_to_file,
        load_transcript_mapping, populate_mutalyzer_cache,
    },
    collate_normalization, collate_parsing,
    compare::rewrite_with_genomic_context,
    compare::{compare_normalize, compare_parse, run_mutalyzer_normalize_parallel, CompareConfig},
    extract_clinvar, extract_json,
    mutalyzer::{has_mutalyzer_normalizer, has_mutalyzer_parser, run_mutalyzer_parser_subprocess},
    normalize_ferro, prepare_references,
    report::{generate_readme_tables, generate_report, generate_summary},
    sample::stratified_sample,
    shard::shard_dataset,
    types::*,
    PrepareConfig,
};
use ferro_hgvs::check::{check_reference, print_check_summary};
use ferro_hgvs::prepare::{
    detect_clinvar_from_ferro_reference, detect_patterns_from_ferro_reference,
};
use std::path::{Path, PathBuf};

#[cfg(feature = "hgvs-rs")]
use ferro_hgvs::benchmark::{check_hgvs_rs_available, check_seqrepo_path, HgvsRsConfig};

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

fn parse_genome_option(s: &str) -> Result<GenomeOption, String> {
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
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
enum Commands {
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
        tool: Tool,

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
        tool: Tool,

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
enum CompareCommands {
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

/// Benchmark subcommands (live sampling and comparison)
#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
enum BenchmarkCommands {
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
enum ExtractCommands {
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
enum SetupCommands {
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
enum GenerateCommands {
    /// Generate summary from comparison files
    Summary {
        #[arg(long)]
        normalization: Vec<PathBuf>,
        #[arg(long)]
        parsing: Vec<PathBuf>,
        #[arg(short, long, default_value = "summary.json")]
        output: PathBuf,
        #[arg(long, action = ArgAction::SetTrue)]
        detailed: bool,
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
enum CollateCommands {
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

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        // =========================================================================
        // GROUPED SUBCOMMAND HANDLERS (Phase 2)
        // =========================================================================
        Commands::Compare(cmd) => match cmd {
            CompareCommands::Results {
                mode,
                results,
                output,
                equivalence,
            } => run_compare_tool(mode, &results, &output, equivalence),

            CompareCommands::Arbitration {
                database,
                category,
                verbose,
            } => run_show_arbitration(&database, category.as_deref(), verbose),
        },

        Commands::Benchmark(cmd) => match cmd {
            BenchmarkCommands::Parse {
                input,
                output,
                sample_size,
                workers,
                seed,
                error_mode,
            } => {
                let config = CompareConfig {
                    sample_size,
                    workers,
                    seed,
                    mutalyzer_settings: None,
                    uta_db_url: None,
                    biocommons_settings: None,
                    hgvs_rs_seqrepo_path: None,
                    hgvs_rs_uta_schema: None,
                    reference_dir: None,
                    transcript_mapping: None,
                    include_unnormalizable: true,
                    include_genomic: true,
                    include_protein: true,
                    allow_mutalyzer_network: false,
                    detailed_output: None,
                    existing_mutalyzer_results: None,
                    validator: Validator::Mutalyzer,
                    error_mode,
                };
                compare_parse(&input, &output, &config).map(|_| ())
            }

            BenchmarkCommands::Normalize {
                input,
                output,
                validator,
                sample_size,
                workers,
                seed,
                reference,
                mutalyzer_settings,
                transcript_mapping,
                include_unnormalizable,
                include_genomic,
                include_protein,
                allow_mutalyzer_network,
                detailed_output,
                existing_mutalyzer_results,
                uta_db_url,
                biocommons_settings,
                hgvs_rs_seqrepo_path,
                hgvs_rs_uta_schema,
                arbitration,
                error_mode,
            } => {
                let validator_enum: Validator = validator.parse().unwrap_or_else(|e| {
                    eprintln!("Error: {}", e);
                    std::process::exit(1);
                });

                let config = CompareConfig {
                    validator: validator_enum,
                    sample_size,
                    workers,
                    seed,
                    mutalyzer_settings,
                    uta_db_url,
                    biocommons_settings,
                    hgvs_rs_seqrepo_path,
                    hgvs_rs_uta_schema,
                    reference_dir: reference,
                    transcript_mapping,
                    include_unnormalizable,
                    include_genomic,
                    include_protein,
                    allow_mutalyzer_network,
                    detailed_output,
                    existing_mutalyzer_results,
                    error_mode,
                };
                compare_normalize(&input, &output, &config).and_then(|result| {
                    if let Some(arb_path) = arbitration {
                        run_apply_arbitration(&arb_path, &result)
                    } else {
                        Ok(())
                    }
                })
            }
        },

        Commands::Extract(cmd) => match cmd {
            ExtractCommands::Clinvar { input, output } => {
                extract_clinvar(&input, &output).map(|_| ())
            }

            ExtractCommands::Json {
                input,
                output,
                format,
            } => {
                let (fmt, field_name) = match format.as_str() {
                    "test_cases_json" => (DatasetFormat::TestCasesJson, "input"),
                    "json_array" => (DatasetFormat::JsonArray, "hgvs"),
                    _ => {
                        eprintln!(
                            "Unknown format: {}. Use 'test_cases_json' or 'json_array'",
                            format
                        );
                        std::process::exit(1);
                    }
                };
                extract_json(&input, &output, fmt, field_name).map(|_| ())
            }

            ExtractCommands::Proteins { input, output } => {
                run_extract_protein_accessions(&input, &output)
            }

            ExtractCommands::Sample {
                input,
                output,
                size,
                seed,
                exclude_protein,
            } => stratified_sample(&input, &output, size, seed.unwrap_or(42), exclude_protein)
                .map(|_| ()),

            ExtractCommands::Shard {
                input,
                output_dir,
                num_shards,
            } => shard_dataset(&input, &output_dir, num_shards).map(|_| ()),
        },

        Commands::Setup(cmd) => match cmd {
            SetupCommands::Uta {
                container_name,
                image_tag,
                port,
                uta_dump,
                force,
            } => {
                let config = BiocommonsLocalConfig {
                    uta_container_name: container_name,
                    uta_image_tag: image_tag,
                    uta_port: port,
                    ..Default::default()
                };
                match setup_uta(&config, force, uta_dump.as_deref()) {
                    Ok(result) => {
                        println!("\nUTA setup complete!");
                        println!("  Container: {}", result.container_name);
                        println!("  Port: {}", result.port);
                        println!("  UTA URL: {}", result.uta_db_url);
                        println!("\nTo use: ferro-benchmark compare normalize --validator biocommons --uta-db-url '{}'", result.uta_db_url);
                        Ok(())
                    }
                    Err(e) => {
                        eprintln!("Error setting up UTA: {}", e);
                        std::process::exit(1);
                    }
                }
            }

            SetupCommands::Seqrepo {
                output_dir,
                instance,
            } => {
                let config = BiocommonsLocalConfig {
                    seqrepo_dir: output_dir.clone(),
                    seqrepo_instance: instance,
                    ..Default::default()
                };
                match setup_seqrepo(&config, false) {
                    Ok(result) => {
                        println!("\nSeqRepo setup complete!");
                        println!("  Directory: {}", result.seqrepo_dir.display());
                        println!("  Instance: {}", result.instance);
                        Ok(())
                    }
                    Err(e) => {
                        eprintln!("Error setting up SeqRepo: {}", e);
                        std::process::exit(1);
                    }
                }
            }

            SetupCommands::StartUta { container_name } => match start_uta(&container_name) {
                Ok(()) => {
                    println!("Container '{}' started", container_name);
                    Ok(())
                }
                Err(e) => {
                    eprintln!("Error starting container: {}", e);
                    std::process::exit(1);
                }
            },

            SetupCommands::StopUta {
                container_name,
                force: _,
            } => match stop_uta(&container_name) {
                Ok(()) => {
                    println!("Container '{}' stopped", container_name);
                    Ok(())
                }
                Err(e) => {
                    eprintln!("Error stopping container: {}", e);
                    std::process::exit(1);
                }
            },
        },

        Commands::Generate(cmd) => match cmd {
            GenerateCommands::Summary {
                normalization,
                parsing,
                output,
                detailed: _,
            } => {
                // Simplified: just use the first dirs if provided
                let parsing_dir = parsing.first().cloned().unwrap_or_default();
                let norm_dir = normalization.first().cloned().unwrap_or_default();
                let config = SummaryConfig {
                    cores: 1,
                    include_mutalyzer: false,
                    normalization_sample_size: 1000,
                };
                generate_summary(&parsing_dir, &norm_dir, &output, config).map(|_| ())
            }

            GenerateCommands::Report { input, output } => load_summary(&input)
                .and_then(|summary_data| generate_report(&summary_data, &output)),

            GenerateCommands::Readme { input } => {
                let output = PathBuf::from("README_TABLES.md");
                load_summary(&input)
                    .and_then(|summary_data| generate_readme_tables(&summary_data, &output))
            }
        },

        Commands::Collate(cmd) => match cmd {
            CollateCommands::Parsing {
                ferro_dir,
                mutalyzer_dir,
                output,
                dataset,
            } => collate_parsing(&ferro_dir, mutalyzer_dir.as_ref(), &output, &dataset).map(|_| ()),

            CollateCommands::Normalization {
                ferro_results,
                mutalyzer_results,
                output,
                dataset,
            } => collate_normalization(
                &ferro_results,
                mutalyzer_results.as_ref(),
                &output,
                &dataset,
            )
            .map(|_| ()),
        },

        // =========================================================================
        // BuildTranscriptMapping and PopulateProteinCache removed (Phase 3)

        // =========================================================================
        // NEW UNIFIED CLI COMMAND HANDLERS
        // =========================================================================
        Commands::PrepareTool {
            tool,
            output_dir,
            genome,
            no_refseqgene,
            no_lrg,
            force,
            ferro_reference,
            patterns,
            clinvar,
            no_proteins,
            shards,
            uta_port,
            seqrepo_dir,
            uta_dump,
            no_load_transcripts,
            no_load_alignments,
        } => run_prepare_tool(
            tool,
            &output_dir,
            genome,
            no_refseqgene,
            no_lrg,
            force,
            ferro_reference.as_ref(),
            patterns.as_ref(),
            clinvar.as_ref(),
            no_proteins,
            shards,
            uta_port,
            seqrepo_dir.as_ref(),
            uta_dump.as_ref(),
            no_load_transcripts,
            no_load_alignments,
        ),

        Commands::ParseTool {
            tool,
            input,
            output,
            timing,
        } => run_parse_tool(tool, &input, &output, timing.as_ref()),

        Commands::NormalizeTool {
            tool,
            input,
            output,
            timing,
            reference,
            mutalyzer_settings,
            biocommons_settings,
            lrg_mapping,
            seqrepo_path,
            workers,
            no_rewrite_intronic,
            allow_network,
        } => run_normalize_tool(
            tool,
            &input,
            &output,
            timing.as_ref(),
            reference.as_ref(),
            mutalyzer_settings.as_ref(),
            biocommons_settings.as_ref(),
            lrg_mapping.as_ref(),
            seqrepo_path.as_ref(),
            workers,
            no_rewrite_intronic,
            allow_network,
        ),

        // CompareTool now handled by Compare::Results
        Commands::CheckTool {
            tool,
            reference,
            mutalyzer_settings,
            uta_db_url,
            seqrepo_path,
        } => run_check_tool(
            tool,
            &reference,
            mutalyzer_settings.as_ref(),
            uta_db_url.as_deref(),
            seqrepo_path.as_ref(),
        ),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

// run_command and parse_dataset_arg removed - Run command was removed

#[allow(clippy::too_many_arguments)]
fn run_populate_cache(
    reference_dir: &PathBuf,
    cache_dir: &PathBuf,
    filter: Option<&PathBuf>,
    patterns: Option<&PathBuf>,
    clinvar: Option<&PathBuf>,
    lrg: bool,
    dry_run: bool,
    fetch_only: bool,
    enhance_nc: bool,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::cache::{
        enhance_nc_annotations, extract_all_accessions_from_file, extract_clinvar_accessions,
        fetch_lrg_sequences, fetch_missing_accessions,
    };
    use ferro_hgvs::data::cdot::CdotMapper;

    // Ensure cache directory exists
    std::fs::create_dir_all(cache_dir).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to create cache dir: {}", e),
    })?;

    // Collect all accessions from patterns and/or clinvar
    let mut all_accessions = std::collections::HashSet::new();

    if let Some(pattern_file) = patterns {
        eprintln!("Extracting accessions from patterns file...");
        let accessions = extract_all_accessions_from_file(pattern_file)?;
        eprintln!("  Found {} unique accessions", accessions.len());
        all_accessions.extend(accessions);
    }

    if let Some(clinvar_file) = clinvar {
        eprintln!("Extracting accessions from ClinVar file...");
        let accessions = extract_clinvar_accessions(clinvar_file)?;
        eprintln!("  Found {} unique accessions", accessions.len());
        all_accessions.extend(accessions);
    }

    let accessions: Vec<String> = all_accessions.into_iter().collect();

    // In dry-run mode, just show what would be fetched
    if dry_run {
        if accessions.is_empty() && !lrg {
            eprintln!("--dry-run requires --patterns, --clinvar, or --lrg");
            std::process::exit(1);
        }

        // Find missing accessions
        let mut missing = Vec::new();
        let mut missing_lrg = Vec::new();

        for acc in &accessions {
            let seq_path = cache_dir.join(format!("{}.sequence", acc));
            if !seq_path.exists() {
                if acc.starts_with("LRG_") {
                    missing_lrg.push(acc.as_str());
                } else {
                    missing.push(acc.as_str());
                }
            }
        }

        // Count missing LRG if --lrg flag
        if lrg {
            for i in 1..=1400 {
                let lrg_id = format!("LRG_{}", i);
                let seq_path = cache_dir.join(format!("{}.sequence", lrg_id));
                if !seq_path.exists() && !missing_lrg.iter().any(|&x| x == lrg_id) {
                    missing_lrg.push(Box::leak(lrg_id.into_boxed_str()));
                }
            }
        }

        println!("\n=== Dry Run Summary ===");
        println!("Total accessions in patterns: {}", accessions.len());
        println!(
            "Missing NCBI accessions to fetch: {}",
            missing.iter().filter(|a| !a.starts_with("CM_")).count()
        );
        println!("Missing LRG accessions to fetch: {}", missing_lrg.len());

        if !missing.is_empty() {
            println!("\nSample missing NCBI accessions (first 20):");
            for acc in missing.iter().filter(|a| !a.starts_with("CM_")).take(20) {
                println!("  {}", acc);
            }
        }

        println!("\nTo fetch these accessions, run without --dry-run");
        return Ok(());
    }

    // In fetch-only mode, skip local cache population
    if !fetch_only {
        // Populate from local reference data
        let stats = populate_mutalyzer_cache(reference_dir, cache_dir, filter)?;
        println!("\nMutalyzer cache populated:");
        println!("  Transcripts processed: {}", stats.transcripts_processed);
        println!("  Cache entries created: {}", stats.cache_entries_created);
        println!(
            "  Legacy sequences fetched: {}",
            stats.legacy_sequences_fetched
        );
    } else if accessions.is_empty() && !lrg {
        eprintln!("--fetch-only requires --patterns, --clinvar, or --lrg");
        std::process::exit(1);
    }

    // Fetch missing accessions from NCBI
    let mut total_fetched = 0usize;

    if !accessions.is_empty() {
        eprintln!(
            "\nChecking {} accessions against cache...",
            accessions.len()
        );
        let fetched = fetch_missing_accessions(&accessions, cache_dir)?;
        total_fetched += fetched;
        println!("  Fetched {} accessions from NCBI", fetched);
    }

    // Fetch LRG sequences (from ferro reference first, then EBI)
    if lrg {
        eprintln!("\nFetching LRG sequences...");
        let lrg_fetched = fetch_lrg_sequences(cache_dir, Some(reference_dir))?;
        total_fetched += lrg_fetched;
        println!("  Added {} LRG sequences to cache", lrg_fetched);
    }

    // Write settings file
    let settings_path = cache_dir.join("mutalyzer_settings.conf");
    let settings_content = format!(
        "MUTALYZER_CACHE_DIR = {}\nMUTALYZER_FILE_CACHE_ADD = false\n",
        cache_dir.display()
    );
    std::fs::write(&settings_path, settings_content).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;

    // Enhance NC_ annotations with transcript selectors if requested
    if enhance_nc {
        eprintln!("\nEnhancing NC_ annotations with transcript selectors...");

        // Load cdot from manifest
        let manifest_path = reference_dir.join("manifest.json");
        if !manifest_path.exists() {
            eprintln!(
                "Warning: Cannot enhance NC_ annotations - manifest not found at {}",
                manifest_path.display()
            );
            eprintln!("Run with --reference-dir pointing to prepared reference data");
        } else {
            let manifest: serde_json::Value = {
                let file = std::fs::File::open(&manifest_path).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to open manifest: {}", e),
                    }
                })?;
                serde_json::from_reader(file).map_err(|e| ferro_hgvs::FerroError::Json {
                    msg: format!("Failed to parse manifest: {}", e),
                })?
            };

            if let Some(cdot_path) = manifest.get("cdot_json").and_then(|v| v.as_str()) {
                // Join with reference_dir since manifest paths are relative
                let full_cdot_path = reference_dir.join(cdot_path);
                eprintln!("Loading cdot from {}...", full_cdot_path.display());
                let cdot = if cdot_path.ends_with(".gz") {
                    CdotMapper::from_json_gz(&full_cdot_path)?
                } else {
                    CdotMapper::from_json_file(&full_cdot_path)?
                };

                let enhanced = enhance_nc_annotations(cache_dir, &cdot)?;
                println!("Enhanced {} NC_ annotation files", enhanced);
            } else {
                eprintln!("Warning: No cdot_json in manifest, skipping NC_ enhancement");
            }
        }
    }

    println!("\n=== Cache Summary ===");
    println!("Cache directory: {}", cache_dir.display());
    println!("Settings file: {}", settings_path.display());
    if total_fetched > 0 {
        println!("Total accessions fetched: {}", total_fetched);
    }
    println!("\nTo use this cache with compare commands:");
    println!("  --mutalyzer-settings {}", settings_path.display());
    Ok(())
}

fn run_enhance_lrg_annotations(
    cache_dir: &Path,
    ferro_reference: Option<&Path>,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::enhance_lrg_annotations;

    let enhanced = enhance_lrg_annotations(cache_dir, ferro_reference)?;

    println!("\n=== Enhancement Summary ===");
    println!("Cache directory: {}", cache_dir.display());
    println!(
        "Enhanced {} LRG annotation files with exon/CDS structure",
        enhanced
    );
    Ok(())
}

fn load_summary(path: &PathBuf) -> Result<ComparisonSummary, ferro_hgvs::FerroError> {
    let file = std::fs::File::open(path).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;
    let reader = std::io::BufReader::new(file);
    serde_json::from_reader(reader).map_err(|e| ferro_hgvs::FerroError::Json {
        msg: format!("Failed to parse {}: {}", path.display(), e),
    })
}

fn run_extract_protein_accessions(
    input_path: &PathBuf,
    output_path: &PathBuf,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::cache::extract_clinvar_protein_accessions;
    use std::io::Write;

    let accessions = extract_clinvar_protein_accessions(input_path)?;

    // Write to output file
    let mut file = std::fs::File::create(output_path).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to create output file: {}", e),
    })?;

    for acc in &accessions {
        writeln!(file, "{}", acc).map_err(|e| ferro_hgvs::FerroError::Io {
            msg: format!("Failed to write accession: {}", e),
        })?;
    }

    println!("\n=== Extraction Summary ===");
    println!("Protein accessions found: {}", accessions.len());
    println!("Output file: {}", output_path.display());

    // Show breakdown by type
    let np_count = accessions.iter().filter(|a| a.starts_with("NP_")).count();
    let xp_count = accessions.iter().filter(|a| a.starts_with("XP_")).count();
    let yp_count = accessions.iter().filter(|a| a.starts_with("YP_")).count();
    let lrg_count = accessions.iter().filter(|a| a.starts_with("LRG_")).count();
    let uniprot_count = accessions.len() - np_count - xp_count - yp_count - lrg_count;

    println!("\nBreakdown:");
    println!("  NP_ (RefSeq protein): {}", np_count);
    println!("  XP_ (predicted protein): {}", xp_count);
    println!("  YP_ (NCBI protein): {}", yp_count);
    println!("  LRG protein: {}", lrg_count);
    println!("  UniProt: {}", uniprot_count);

    println!("\nTo populate protein cache:");
    println!(
        "  ferro-benchmark populate-protein-cache -i {} -o protein_cache",
        output_path.display()
    );

    Ok(())
}

fn run_populate_protein_cache(
    input_path: &PathBuf,
    cache_dir: &PathBuf,
    is_clinvar: bool,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::cache::{
        extract_clinvar_protein_accessions, populate_protein_cache,
    };
    use std::io::BufRead;

    // Get accessions either from ClinVar or from accession list file
    let accessions: Vec<String> = if is_clinvar {
        extract_clinvar_protein_accessions(input_path)?
    } else {
        // Read accessions from file (one per line)
        let file = std::fs::File::open(input_path).map_err(|e| ferro_hgvs::FerroError::Io {
            msg: format!("Failed to open accession file: {}", e),
        })?;
        let reader = std::io::BufReader::new(file);
        reader
            .lines()
            .map_while(Result::ok)
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect()
    };

    eprintln!("Found {} protein accessions to process", accessions.len());

    // Populate the cache
    let stats = populate_protein_cache(&accessions, cache_dir)?;

    println!("\n=== Protein Cache Summary ===");
    println!("Cache directory: {}", cache_dir.display());
    println!("Total accessions: {}", stats.total);
    println!("Already cached: {}", stats.already_cached);
    println!("Fetched from NCBI: {}", stats.ncbi_fetched);
    println!("Fetched from UniProt: {}", stats.uniprot_fetched);
    if stats.failed > 0 {
        println!("Failed/skipped: {}", stats.failed);
    }

    println!("\nTo use this cache with normalization, set up a provider");
    println!("that implements get_protein_sequence() using these cached FASTAs.");

    Ok(())
}

fn run_show_arbitration(
    database_path: &Path,
    category_filter: Option<&str>,
    verbose: bool,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::{ArbitrationCategory, ArbitrationDatabase};

    let db = ArbitrationDatabase::load(database_path)?;
    let stats = db.stats();

    println!("\n=== Arbitration Database ===");
    println!("File: {}", database_path.display());
    println!("Version: {}", db.version);
    println!("\nStatistics:");
    println!("  Total decisions:   {}", stats.total);
    println!("  Ferro correct:     {}", stats.ferro_correct);
    println!("  Mutalyzer correct: {}", stats.mutalyzer_correct);
    println!("  Equivalent:        {}", stats.equivalent);
    println!("  Both incorrect:    {}", stats.both_incorrect);
    println!("  Unknown:           {}", stats.unknown);

    // Get decisions to display
    let decisions: Vec<_> = match category_filter {
        Some("ferro_correct") => db.by_category(ArbitrationCategory::FerroCorrect),
        Some("mutalyzer_correct") => db.by_category(ArbitrationCategory::MutalyzerCorrect),
        Some("equivalent") => db.by_category(ArbitrationCategory::Equivalent),
        Some("both_incorrect") => db.by_category(ArbitrationCategory::BothIncorrect),
        Some("unknown") => db.by_category(ArbitrationCategory::Unknown),
        Some(other) => {
            eprintln!(
                "Unknown category: {}. Use: ferro_correct, mutalyzer_correct, equivalent, both_incorrect, unknown",
                other
            );
            std::process::exit(1);
        }
        None => db.decisions.values().collect(),
    };

    if !decisions.is_empty() {
        println!("\nDecisions:");
        for (i, decision) in decisions.iter().enumerate() {
            println!("\n{}. {}", i + 1, decision.input);
            println!("   Category: {}", decision.category);
            println!("   Ferro:     {}", decision.ferro_output);
            println!("   Mutalyzer: {}", decision.mutalyzer_output);
            if verbose {
                println!("   Reason: {}", decision.reason);
                if let Some(ref spec) = decision.spec_reference {
                    println!("   Spec: {}", spec);
                }
            }
        }
    }

    Ok(())
}

fn run_apply_arbitration(
    arbitration_path: &Path,
    result: &ComparisonResult,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::{apply_arbitrations, ArbitrationDatabase};

    let db = ArbitrationDatabase::load(arbitration_path)?;

    // Convert differences to the format expected by apply_arbitrations
    let differences: Vec<(String, String, String)> = result
        .differences
        .iter()
        .map(|d| {
            (
                d.input.clone(),
                d.ferro_output.clone(),
                d.mutalyzer_output.clone(),
            )
        })
        .collect();

    let summary = apply_arbitrations(&db, &differences);
    summary.print_report();

    Ok(())
}

fn run_normalize_biocommons(
    input: &Path,
    results_path: &Path,
    timing_path: &Path,
    uta_db_url: Option<&str>,
    seqrepo_dir: Option<&str>,
    lrg_mapping_file: Option<&str>,
    workers: usize,
) -> Result<(), ferro_hgvs::FerroError> {
    use std::time::Instant;

    // Check if biocommons is available
    if !has_biocommons_normalizer() {
        eprintln!("biocommons/hgvs is NOT available");
        eprintln!("Install with: pip install hgvs");
        std::process::exit(1);
    }

    let start = Instant::now();

    // Run normalization (parallel or sequential based on workers)
    if workers > 1 {
        run_biocommons_normalizer_parallel(
            input.to_str().unwrap(),
            results_path.to_str().unwrap(),
            uta_db_url,
            seqrepo_dir,
            lrg_mapping_file,
            workers,
        )?;
    } else {
        run_biocommons_normalizer_subprocess(
            input.to_str().unwrap(),
            results_path.to_str().unwrap(),
            uta_db_url,
            seqrepo_dir,
            lrg_mapping_file,
        )?;
    }

    let elapsed = start.elapsed();

    // Read results to get counts
    let results_content =
        std::fs::read_to_string(results_path).map_err(|e| ferro_hgvs::FerroError::Io {
            msg: format!("Failed to read results: {}", e),
        })?;
    let results_data: serde_json::Value =
        serde_json::from_str(&results_content).map_err(|e| ferro_hgvs::FerroError::Json {
            msg: format!("Failed to parse results JSON: {}", e),
        })?;

    let total = results_data["total_patterns"].as_u64().unwrap_or(0) as usize;
    let successful = results_data["successful"].as_u64().unwrap_or(0) as usize;

    // Write timing info
    let timing_info = TimingInfo::new("biocommons-hgvs", total, successful, elapsed);
    let timing_json =
        serde_json::to_string_pretty(&timing_info).map_err(|e| ferro_hgvs::FerroError::Json {
            msg: format!("Failed to serialize timing: {}", e),
        })?;
    std::fs::write(timing_path, timing_json).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to write timing file: {}", e),
    })?;

    println!(
        "Normalized {} patterns in {:.2}s ({:.0} patterns/s)",
        total,
        elapsed.as_secs_f64(),
        total as f64 / elapsed.as_secs_f64()
    );
    println!(
        "Successful: {} ({:.1}%)",
        successful,
        100.0 * successful as f64 / total as f64
    );

    Ok(())
}

#[cfg(feature = "hgvs-rs")]
#[allow(clippy::too_many_arguments)]
fn run_normalize_hgvs_rs(
    input: &Path,
    results_path: &Path,
    timing_path: &Path,
    uta_db_url: &str,
    uta_db_schema: &str,
    seqrepo_path: &str,
    lrg_mapping: Option<&Path>,
    workers: usize,
) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::{run_hgvs_rs_normalize_parallel, HgvsRsConfig};
    use std::io::BufRead;

    // Read patterns from input file
    let file = std::fs::File::open(input).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to open input file: {}", e),
    })?;
    let reader = std::io::BufReader::new(file);
    let patterns: Vec<String> = reader
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();

    let total = patterns.len();
    eprintln!(
        "Normalizing {} patterns with hgvs-rs ({} workers)...",
        total, workers
    );

    let config = HgvsRsConfig {
        uta_db_url: uta_db_url.to_string(),
        uta_db_schema: uta_db_schema.to_string(),
        seqrepo_path: seqrepo_path.to_string(),
        lrg_mapping_file: lrg_mapping.map(|p| p.to_string_lossy().to_string()),
    };

    // Run normalization (parallel or sequential)
    let (results, elapsed, error_counts) =
        run_hgvs_rs_normalize_parallel(&patterns, &config, workers)?;

    // Count successes
    let successful = results.iter().filter(|r| r.success).count();

    // Build results JSON
    let results_json = serde_json::json!({
        "tool": "hgvs-rs",
        "total_patterns": total,
        "successful": successful,
        "failed": total - successful,
        "elapsed_seconds": elapsed.as_secs_f64(),
        "throughput": total as f64 / elapsed.as_secs_f64(),
        "workers": workers,
        "error_counts": error_counts,
        "results": results.iter().map(|r| {
            serde_json::json!({
                "input": r.input,
                "success": r.success,
                "output": r.output,
                "error": r.error,
                "error_category": r.error_category,
            })
        }).collect::<Vec<_>>(),
    });

    // Write results
    std::fs::write(
        results_path,
        serde_json::to_string_pretty(&results_json).map_err(|e| ferro_hgvs::FerroError::Json {
            msg: format!("Failed to serialize results: {}", e),
        })?,
    )
    .map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to write results file: {}", e),
    })?;

    // Write timing info
    let timing_info = TimingInfo::new("hgvs-rs", total, successful, elapsed);
    let timing_json =
        serde_json::to_string_pretty(&timing_info).map_err(|e| ferro_hgvs::FerroError::Json {
            msg: format!("Failed to serialize timing: {}", e),
        })?;
    std::fs::write(timing_path, timing_json).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to write timing file: {}", e),
    })?;

    println!(
        "Normalized {} patterns in {:.2}s ({:.0} patterns/s)",
        total,
        elapsed.as_secs_f64(),
        total as f64 / elapsed.as_secs_f64()
    );
    println!(
        "Successful: {} ({:.1}%)",
        successful,
        100.0 * successful as f64 / total as f64
    );
    if workers > 1 {
        println!("Workers: {}", workers);
    }

    Ok(())
}

/// Derive protein sequences from transcripts and store in ferro reference.
///
/// This allows `prepare mutalyzer` to use pre-derived proteins without
/// needing to derive them itself or fetch from NCBI.
fn derive_proteins_for_ferro(ferro_dir: &Path) -> Result<(), ferro_hgvs::FerroError> {
    use ferro_hgvs::benchmark::translate::derive_protein_from_transcript;
    use ferro_hgvs::data::cdot::CdotMapper;
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader};

    // Create proteins directory
    let proteins_dir = ferro_dir.join("proteins");
    std::fs::create_dir_all(&proteins_dir).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to create proteins dir: {}", e),
    })?;

    // Load cdot for CDS coordinates
    let cdot_path = ferro_dir.join("cdot/cdot-0.2.32.refseq.GRCh38.json");
    if !cdot_path.exists() {
        eprintln!("  Warning: cdot not found, skipping protein derivation");
        return Ok(());
    }
    println!("  Loading cdot for CDS coordinates...");
    let cdot = CdotMapper::from_json_file(&cdot_path)?;

    // Load supplemental CDS metadata if available
    let mut supplemental_cds: HashMap<String, (Option<u64>, Option<u64>, String)> = HashMap::new();
    let supplemental_dir = ferro_dir.join("supplemental");
    if supplemental_dir.exists() {
        for entry in std::fs::read_dir(&supplemental_dir)
            .into_iter()
            .flatten()
            .flatten()
        {
            let path = entry.path();
            if path.extension().map(|e| e == "cds").unwrap_or(false) {
                if let Ok(file) = std::fs::File::open(&path) {
                    for line in BufReader::new(file).lines().map_while(Result::ok) {
                        let parts: Vec<&str> = line.split('\t').collect();
                        if parts.len() >= 4 {
                            let acc = parts[0].to_string();
                            let start = parts[1].parse().ok();
                            let end = parts[2].parse().ok();
                            let gene = parts[3].to_string();
                            supplemental_cds.insert(acc, (start, end, gene));
                        }
                    }
                }
            }
        }
    }
    if !supplemental_cds.is_empty() {
        println!(
            "  Loaded CDS metadata for {} supplemental transcripts",
            supplemental_cds.len()
        );
    }

    // Load transcript sequences from FASTA files
    println!("  Loading transcript sequences...");
    let mut transcript_seqs: HashMap<String, String> = HashMap::new();

    // Load from transcripts directory
    let transcripts_dir = ferro_dir.join("transcripts");
    if transcripts_dir.exists() {
        for entry in std::fs::read_dir(&transcripts_dir)
            .into_iter()
            .flatten()
            .flatten()
        {
            let path = entry.path();
            if path
                .extension()
                .map(|e| e == "fna" || e == "fa" || e == "fasta")
                .unwrap_or(false)
            {
                load_fasta_sequences(&path, &mut transcript_seqs)?;
            }
        }
    }

    // Load from supplemental directory
    if supplemental_dir.exists() {
        for entry in std::fs::read_dir(&supplemental_dir)
            .into_iter()
            .flatten()
            .flatten()
        {
            let path = entry.path();
            if path
                .extension()
                .map(|e| e == "fna" || e == "fa" || e == "fasta")
                .unwrap_or(false)
            {
                load_fasta_sequences(&path, &mut transcript_seqs)?;
            }
        }
    }

    println!("  Loaded {} transcript sequences", transcript_seqs.len());

    // Derive proteins
    let mut proteins_derived = 0usize;
    let mut proteins_skipped = 0usize;

    for (tx_id, seq) in &transcript_seqs {
        // Skip non-coding transcripts
        if !tx_id.starts_with("NM_") && !tx_id.starts_with("XM_") {
            continue;
        }

        // Get CDS coordinates from cdot or supplemental metadata
        let (cds_start, cds_end, protein_id) = if let Some(tx) = cdot.get_transcript(tx_id) {
            let start = tx.cds_start.unwrap_or(0) as usize;
            let end = tx.cds_end.unwrap_or(seq.len() as u64) as usize;
            let prot_id = tx
                .protein
                .clone()
                .unwrap_or_else(|| tx_id.replace("NM_", "NP_").replace("XM_", "XP_"));
            (start, end, prot_id)
        } else if let Some((start, end, _gene)) = supplemental_cds.get(tx_id) {
            let start = start.map(|s| (s - 1) as usize).unwrap_or(0);
            let end = end.map(|e| e as usize).unwrap_or(seq.len());
            let prot_id = tx_id.replace("NM_", "NP_").replace("XM_", "XP_");
            (start, end, prot_id)
        } else {
            continue;
        };

        // Check if already derived
        let prot_path = proteins_dir.join(format!("{}.faa", protein_id));
        if prot_path.exists() {
            proteins_skipped += 1;
            continue;
        }

        // Derive protein sequence
        if let Some(protein_seq) = derive_protein_from_transcript(seq, cds_start, cds_end) {
            // Write protein to FASTA file
            let content = format!(">{}\n{}\n", protein_id, protein_seq);
            std::fs::write(&prot_path, content).map_err(|e| ferro_hgvs::FerroError::Io {
                msg: format!("Failed to write protein: {}", e),
            })?;
            proteins_derived += 1;
        }
    }

    println!(
        "  Derived {} protein sequences ({} already exist)",
        proteins_derived, proteins_skipped
    );
    println!("  Proteins stored in: {}", proteins_dir.display());

    Ok(())
}

/// Load sequences from a FASTA file into a HashMap
fn load_fasta_sequences(
    path: &Path,
    seqs: &mut std::collections::HashMap<String, String>,
) -> Result<(), ferro_hgvs::FerroError> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to open FASTA: {}", e),
    })?;
    let reader = BufReader::new(file);

    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.map_err(|e| ferro_hgvs::FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if let Some(header) = line.strip_prefix('>') {
            // Save previous sequence
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    seqs.insert(id, std::mem::take(&mut current_seq));
                }
            }
            // Parse new ID (take first whitespace-delimited token after >)
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            current_id = Some(id);
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }

    // Save last sequence
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            seqs.insert(id, current_seq);
        }
    }

    Ok(())
}

// =============================================================================
// NEW UNIFIED CLI HANDLERS
// =============================================================================

/// Print summary of what prepare ferro will do
fn print_ferro_summary(
    genome: &GenomeOption,
    no_refseqgene: bool,
    no_lrg: bool,
    output_dir: &Path,
) {
    println!("=== Preparing ferro reference data ===");
    println!("Output directory: {}", output_dir.display());
    println!();
    println!("Downloads:");
    println!("  [✓] Transcripts (RefSeq + Ensembl)");
    println!("  [✓] cdot transcript data");

    match genome {
        GenomeOption::All => {
            println!("  [✓] GRCh38 genome (~3GB)");
            println!("  [✓] GRCh37 genome (~900MB)");
        }
        GenomeOption::Grch38 => {
            println!("  [✓] GRCh38 genome (~3GB)");
            println!("  [ ] GRCh37 genome (--genome grch38)");
        }
        GenomeOption::Grch37 => {
            println!("  [ ] GRCh38 genome (--genome grch37)");
            println!("  [✓] GRCh37 genome (~900MB)");
        }
        GenomeOption::None => {
            println!("  [ ] GRCh38 genome (--genome none)");
            println!("  [ ] GRCh37 genome (--genome none)");
        }
    }

    if no_refseqgene {
        println!("  [ ] RefSeqGene sequences (--no-refseqgene)");
    } else {
        println!("  [✓] RefSeqGene sequences (~600MB)");
    }

    if no_lrg {
        println!("  [ ] LRG sequences (--no-lrg)");
    } else {
        println!("  [✓] LRG sequences (~50MB)");
    }

    println!();
}

/// Print summary of what prepare mutalyzer will do
fn print_mutalyzer_summary(
    ferro_ref: &Path,
    cache_dir: &Path,
    patterns: Option<&PathBuf>,
    clinvar: Option<&PathBuf>,
    no_lrg: bool,
    no_proteins: bool,
) {
    println!("=== Preparing mutalyzer cache ===");
    println!("Ferro reference: {}", ferro_ref.display());
    println!("Cache directory: {}", cache_dir.display());
    println!();
    println!("Steps:");
    println!("  [✓] Populate cache from ferro reference");
    println!("  [✓] Generate transcript→chromosome mapping (for intronic rewriting)");
    println!("  [✓] Enhance NC_ annotations with transcript mappings");

    if no_lrg {
        println!("  [ ] Enhance LRG annotations (--no-lrg)");
    } else {
        println!("  [✓] Enhance LRG annotations");
    }

    if let Some(p) = patterns {
        println!("  [✓] Pre-cache accessions from: {}", p.display());
    }

    if no_proteins {
        println!("  [ ] Populate protein cache (--no-proteins)");
    } else {
        let ferro_proteins = ferro_ref.join("proteins");
        let protein_count = std::fs::read_dir(&ferro_proteins)
            .into_iter()
            .flatten()
            .flatten()
            .filter(|e| {
                e.path()
                    .extension()
                    .map(|ext| ext == "faa")
                    .unwrap_or(false)
            })
            .count();

        if protein_count > 0 {
            println!("  [✓] Copy {} proteins from ferro reference", protein_count);
        } else if clinvar.is_some() {
            println!(
                "  [✓] Derive proteins from ClinVar (legacy, consider 'prepare ferro --clinvar')"
            );
        } else {
            println!("  [!] No protein source (run 'prepare ferro --clinvar' first)");
        }
    }

    println!();
}

/// Print summary of what prepare biocommons will do
fn print_biocommons_summary(seqrepo_dir: &Path, uta_port: u16) {
    println!("=== Preparing biocommons (UTA + SeqRepo) ===");
    println!("SeqRepo directory: {}", seqrepo_dir.display());
    println!("UTA port: {}", uta_port);
    println!();
    println!("Steps:");
    println!("  [✓] Check/setup UTA database (Docker)");
    println!("  [✓] Check/setup SeqRepo repository");
    println!("  [✓] Generate biocommons settings file");
    println!();
}

/// Print summary of what prepare hgvs-rs will do
fn print_hgvs_rs_summary(seqrepo_dir: &Path, uta_port: u16) {
    println!("=== Preparing hgvs-rs (UTA + SeqRepo) ===");
    println!("SeqRepo directory: {}", seqrepo_dir.display());
    println!("UTA port: {}", uta_port);
    println!();
    println!("Steps:");
    println!("  [✓] Check/setup UTA database (Docker)");
    println!("  [✓] Check/setup SeqRepo repository");
    println!();
}

/// Create a shard manifest for the mutalyzer cache.
///
/// The manifest records the number of shards and shard distribution stats.
/// At normalize time, patterns are grouped by their accession's shard hash,
/// so each worker only loads ~1/N of the accessions into memory (mutalyzer
/// loads sequences on-demand).
///
/// This approach avoids creating millions of symlinks - we just create a
/// manifest file and route patterns at normalize time.
fn shard_mutalyzer_cache(
    cache_dir: &Path,
    num_shards: usize,
) -> Result<(), ferro_hgvs::FerroError> {
    use std::collections::HashSet;
    use std::hash::{Hash, Hasher};

    // Collect all accessions from cache files
    // Files are named: {accession}.sequence, {accession}.annotations
    let mut accessions: HashSet<String> = HashSet::new();

    println!("  Scanning cache directory for accessions...");
    for entry in std::fs::read_dir(cache_dir).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to read cache directory: {}", e),
    })? {
        let entry = entry.map_err(|e| ferro_hgvs::FerroError::Io {
            msg: format!("Failed to read directory entry: {}", e),
        })?;
        let path = entry.path();

        // Skip directories and non-cache files
        if path.is_dir() {
            continue;
        }

        let filename = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n,
            None => continue,
        };

        // Extract accession from filename
        if let Some(accession) = filename
            .strip_suffix(".sequence")
            .or_else(|| filename.strip_suffix(".annotations"))
        {
            accessions.insert(accession.to_string());
        }
    }

    println!("  Found {} unique accessions in cache", accessions.len());

    // Hash function: stable hash of accession string
    fn accession_to_shard(accession: &str, num_shards: usize) -> usize {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        accession.hash(&mut hasher);
        (hasher.finish() as usize) % num_shards
    }

    // Count accessions per shard (don't store the full mapping - too large)
    let mut shard_counts = vec![0usize; num_shards];
    for accession in &accessions {
        let shard_id = accession_to_shard(accession, num_shards);
        shard_counts[shard_id] += 1;
    }

    // Write shard manifest (lightweight - just stats, not full mapping)
    let manifest_path = cache_dir.join("shard_manifest.json");
    let manifest_json = serde_json::json!({
        "num_shards": num_shards,
        "total_accessions": accessions.len(),
        "shard_counts": shard_counts,
    });
    let manifest_file =
        std::fs::File::create(&manifest_path).map_err(|e| ferro_hgvs::FerroError::Io {
            msg: format!("Failed to create manifest: {}", e),
        })?;
    serde_json::to_writer_pretty(manifest_file, &manifest_json).map_err(|e| {
        ferro_hgvs::FerroError::Json {
            msg: format!("Failed to write manifest: {}", e),
        }
    })?;

    // Write single settings file pointing to full cache
    let abs_cache_dir = cache_dir
        .canonicalize()
        .unwrap_or_else(|_| cache_dir.to_path_buf());
    let settings_path = cache_dir.join("mutalyzer_settings.conf");
    let settings_content = format!(
        "# Mutalyzer settings\n\
         MUTALYZER_CACHE_DIR={}\n\
         MUTALYZER_FILE_CACHE_ADD=False\n",
        abs_cache_dir.display()
    );
    std::fs::write(&settings_path, settings_content).map_err(|e| ferro_hgvs::FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;

    // Report distribution
    let min_count = shard_counts.iter().min().unwrap_or(&0);
    let max_count = shard_counts.iter().max().unwrap_or(&0);
    let avg_count = accessions.len() / num_shards;

    println!(
        "  Created manifest for {} shards (~{} accessions/shard, range {}-{})",
        num_shards, avg_count, min_count, max_count
    );
    println!("  Manifest: {}", manifest_path.display());

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_prepare_tool(
    tool: Tool,
    output_dir: &Path,
    genome: GenomeOption,
    no_refseqgene: bool,
    no_lrg: bool,
    force: bool,
    ferro_reference: Option<&PathBuf>,
    patterns: Option<&PathBuf>,
    clinvar: Option<&PathBuf>,
    no_proteins: bool,
    shards: usize,
    uta_port: u16,
    seqrepo_dir: Option<&PathBuf>,
    uta_dump: Option<&PathBuf>,
    no_load_transcripts: bool,
    no_load_alignments: bool,
) -> Result<(), ferro_hgvs::FerroError> {
    match tool {
        Tool::Ferro => {
            // Print summary of what will be done
            print_ferro_summary(&genome, no_refseqgene, no_lrg, output_dir);

            let config = PrepareConfig {
                output_dir: output_dir.to_path_buf(),
                download_transcripts: true,
                download_genome: matches!(genome, GenomeOption::Grch38 | GenomeOption::All),
                download_genome_grch37: matches!(genome, GenomeOption::Grch37 | GenomeOption::All),
                download_refseqgene: !no_refseqgene,
                download_lrg: !no_lrg,
                download_cdot: true,
                skip_existing: !force,
                clinvar_file: None,
                patterns_file: None,
                dry_run: false,
            };
            let manifest = prepare_references(&config)?;
            println!("\n=== Preparation complete ===");
            println!("Manifest: {}/manifest.json", output_dir.display());
            println!("\nReference data prepared:");
            println!("  Transcripts: {}", manifest.transcript_count);
            println!("  Prefixes: {}", manifest.available_prefixes.join(", "));
            if let Some(cdot) = &manifest.cdot_json {
                println!("  cdot: {}", cdot.display());
            }

            // Fetch missing accessions from patterns if provided
            if let Some(patterns_path) = patterns {
                println!("\n=== Fetching missing accessions from patterns ===");

                // Copy patterns file to ferro reference for provenance and auto-detection
                let patterns_dir = output_dir.join("patterns");
                std::fs::create_dir_all(&patterns_dir).map_err(|e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to create patterns dir: {}", e),
                })?;
                let patterns_filename = patterns_path
                    .file_name()
                    .unwrap_or_else(|| std::ffi::OsStr::new("patterns.txt"));
                let stored_patterns = patterns_dir.join(patterns_filename);
                std::fs::copy(patterns_path, &stored_patterns).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to copy patterns file: {}", e),
                    }
                })?;
                println!(
                    "  Stored patterns file for provenance: {}",
                    stored_patterns.display()
                );

                // Extract all accessions from patterns
                let accessions = extract_all_accessions_from_file(patterns_path)?;
                eprintln!("Found {} unique accessions in patterns", accessions.len());

                // Filter to nucleotide transcripts (NG_, NM_, NR_, XM_, XR_)
                // XM_/XR_ are predicted (model) transcripts that may appear in ClinVar
                let nucleotide_accessions: Vec<String> = accessions
                    .into_iter()
                    .filter(|a| {
                        a.starts_with("NG_")
                            || a.starts_with("NM_")
                            || a.starts_with("NR_")
                            || a.starts_with("XM_")
                            || a.starts_with("XR_")
                    })
                    .collect();

                if !nucleotide_accessions.is_empty() {
                    println!(
                        "  {} nucleotide accessions to fetch",
                        nucleotide_accessions.len()
                    );
                    let supplemental_dir = output_dir.join("supplemental");
                    std::fs::create_dir_all(&supplemental_dir).map_err(|e| {
                        ferro_hgvs::FerroError::Io {
                            msg: format!("Failed to create supplemental dir: {}", e),
                        }
                    })?;
                    let supplemental_fasta = supplemental_dir.join("patterns_transcripts.fna");

                    // Fetch accessions from NCBI and save to FASTA
                    match fetch_fasta_to_file(&nucleotide_accessions, &supplemental_fasta, false) {
                        Ok(fetched) => {
                            println!(
                                "  Fetched {} accessions to {}",
                                fetched,
                                supplemental_fasta.display()
                            );
                        }
                        Err(e) => {
                            eprintln!("  Warning: Failed to fetch some accessions: {}", e);
                        }
                    }
                } else {
                    println!("  No nucleotide accessions found in patterns to fetch");
                }
            }

            // Store ClinVar file reference if provided (symlink to avoid copying large file)
            if let Some(clinvar_path) = clinvar {
                println!("\n=== Storing ClinVar reference ===");
                let clinvar_dir = output_dir.join("clinvar");
                std::fs::create_dir_all(&clinvar_dir).map_err(|e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to create clinvar dir: {}", e),
                })?;

                let clinvar_filename = clinvar_path
                    .file_name()
                    .unwrap_or_else(|| std::ffi::OsStr::new("hgvs4variation.txt.gz"));
                let stored_clinvar = clinvar_dir.join(clinvar_filename);

                // Use symlink to avoid copying the large file
                // First, get absolute path of source
                let abs_clinvar_path = if clinvar_path.is_absolute() {
                    clinvar_path.to_path_buf()
                } else {
                    std::env::current_dir()
                        .map_err(|e| ferro_hgvs::FerroError::Io {
                            msg: format!("Failed to get current dir: {}", e),
                        })?
                        .join(clinvar_path)
                };

                // Remove existing symlink if present
                if stored_clinvar.exists() || stored_clinvar.is_symlink() {
                    std::fs::remove_file(&stored_clinvar).ok();
                }

                #[cfg(unix)]
                std::os::unix::fs::symlink(&abs_clinvar_path, &stored_clinvar).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to symlink clinvar file: {}", e),
                    }
                })?;

                #[cfg(not(unix))]
                std::fs::copy(clinvar_path, &stored_clinvar).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to copy clinvar file: {}", e),
                    }
                })?;

                println!(
                    "  Stored ClinVar reference: {} -> {}",
                    stored_clinvar.display(),
                    abs_clinvar_path.display()
                );

                // Derive proteins from transcripts
                println!("\n=== Deriving protein sequences from transcripts ===");
                derive_proteins_for_ferro(output_dir)?;
            }

            Ok(())
        }

        Tool::Mutalyzer => {
            // Prerequisite: ferro reference data must exist
            let ferro_ref = ferro_reference.ok_or_else(|| ferro_hgvs::FerroError::Io {
                msg: "Mutalyzer preparation requires ferro reference data.\n\n\
                      First run: ferro-benchmark prepare ferro --output-dir <dir>\n\
                      Then run:  ferro-benchmark prepare mutalyzer --ferro-reference <dir>"
                    .to_string(),
            })?;

            // Check if ferro reference exists
            let manifest_path = ferro_ref.join("manifest.json");
            if !manifest_path.exists() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: format!(
                        "Ferro reference data not found at {}\n\n\
                         First run: ferro-benchmark prepare ferro --output-dir {}",
                        manifest_path.display(),
                        ferro_ref.display()
                    ),
                });
            }

            let cache_dir = output_dir.to_path_buf();

            // Auto-detect patterns from ferro reference if not explicitly provided
            let detected_patterns = if patterns.is_none() {
                detect_patterns_from_ferro_reference(ferro_ref)
            } else {
                None
            };
            let effective_patterns = patterns.or(detected_patterns.as_ref());
            if let Some(path) = detected_patterns.as_ref() {
                println!(
                    "Auto-detected patterns from ferro reference: {}",
                    path.display()
                );
            }

            // Auto-detect clinvar from ferro reference if not explicitly provided
            let detected_clinvar = if clinvar.is_none() {
                detect_clinvar_from_ferro_reference(ferro_ref)
            } else {
                None
            };
            let effective_clinvar = clinvar.or(detected_clinvar.as_ref());
            if let Some(path) = detected_clinvar.as_ref() {
                println!(
                    "Auto-detected ClinVar from ferro reference: {}",
                    path.display()
                );
            }

            // Check if ferro has pre-derived proteins
            let ferro_proteins = ferro_ref.join("proteins");
            let has_ferro_proteins = ferro_proteins.exists()
                && std::fs::read_dir(&ferro_proteins)
                    .into_iter()
                    .flatten()
                    .flatten()
                    .any(|e| {
                        e.path()
                            .extension()
                            .map(|ext| ext == "faa")
                            .unwrap_or(false)
                    });

            // Protein cache validation: require explicit choice (unless auto-detected or ferro has proteins)
            if !has_ferro_proteins && effective_clinvar.is_none() && !no_proteins {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "Protein cache requires a source.\n\n\
                          Either:\n\
                          - Run 'prepare ferro --clinvar <path>' first (recommended), or\n\
                          - Provide --clinvar <path> to this command, or\n\
                          - Use --no-proteins to explicitly skip protein cache population."
                        .to_string(),
                });
            }

            // Print summary of what will be done
            print_mutalyzer_summary(
                ferro_ref,
                &cache_dir,
                effective_patterns,
                effective_clinvar,
                no_lrg,
                no_proteins,
            );

            // Run populate cache with enhance_nc enabled
            // Pass patterns file to extract and pre-cache all referenced accessions (including NP_ proteins)
            run_populate_cache(
                ferro_ref,
                &cache_dir,
                None, // filter
                effective_patterns,
                None, // clinvar - handled separately below
                !no_lrg,
                false, // dry_run
                false, // fetch_only
                true,  // enhance_nc
            )?;

            // Generate transcript→chromosome mapping for intronic variant rewriting
            println!("\nGenerating transcript→chromosome mapping...");
            let cdot_path = ferro_ref.join("cdot");
            let cdot_file = std::fs::read_dir(&cdot_path).ok().and_then(|entries| {
                entries.filter_map(|e| e.ok()).find_map(|entry| {
                    let name = entry.file_name().to_string_lossy().to_string();
                    if name.contains("grch38") && name.ends_with(".json.gz") {
                        Some(entry.path())
                    } else {
                        None
                    }
                })
            });

            if let Some(cdot_file) = cdot_file {
                let mapping_path = cache_dir.join("transcript_mapping.json");
                match build_transcript_chromosome_mapping(
                    &cdot_file,
                    None::<&PathBuf>, // gene2refseq
                    &mapping_path,
                    "GRCh38",
                ) {
                    Ok(stats) => {
                        println!("  Generated {} transcript mappings", stats.total);
                    }
                    Err(e) => {
                        eprintln!("  Warning: Failed to generate mapping: {}", e);
                    }
                }
            } else {
                eprintln!("  Warning: cdot file not found, skipping mapping generation");
            }

            // Also enhance LRG annotations unless skipped
            if !no_lrg {
                println!("\nEnhancing LRG annotations...");
                run_enhance_lrg_annotations(&cache_dir, Some(ferro_ref))?;
            }

            // Copy proteins from ferro reference if available
            if !no_proteins {
                let ferro_proteins = ferro_ref.join("proteins");
                let mutalyzer_proteins = cache_dir.join("proteins");

                if ferro_proteins.exists() && ferro_proteins.is_dir() {
                    // Count proteins in ferro reference
                    let protein_count: usize = std::fs::read_dir(&ferro_proteins)
                        .into_iter()
                        .flatten()
                        .flatten()
                        .filter(|e| {
                            e.path()
                                .extension()
                                .map(|ext| ext == "faa")
                                .unwrap_or(false)
                        })
                        .count();

                    if protein_count > 0 {
                        println!(
                            "\nCopying {} protein sequences from ferro reference...",
                            protein_count
                        );
                        std::fs::create_dir_all(&mutalyzer_proteins).map_err(|e| {
                            ferro_hgvs::FerroError::Io {
                                msg: format!("Failed to create proteins directory: {}", e),
                            }
                        })?;

                        let mut copied = 0;
                        for entry in std::fs::read_dir(&ferro_proteins)
                            .into_iter()
                            .flatten()
                            .flatten()
                        {
                            let path = entry.path();
                            if path.extension().map(|ext| ext == "faa").unwrap_or(false) {
                                let file_name = match path.file_name() {
                                    Some(name) => name,
                                    None => continue, // Skip paths without a file name
                                };
                                let dest = mutalyzer_proteins.join(file_name);
                                std::fs::copy(&path, &dest).map_err(|e| {
                                    ferro_hgvs::FerroError::Io {
                                        msg: format!("Failed to copy protein file: {}", e),
                                    }
                                })?;
                                copied += 1;
                            }
                        }
                        println!("  Copied {} protein sequences to mutalyzer cache", copied);
                    } else {
                        println!(
                            "\nNo proteins in ferro reference (run 'prepare ferro --clinvar')"
                        );
                    }
                } else if let Some(clinvar_path) = effective_clinvar {
                    // Fall back to deriving proteins from clinvar (legacy behavior)
                    println!("\nDeriving protein sequences from ClinVar (consider using 'prepare ferro --clinvar' first)...");
                    run_populate_protein_cache(clinvar_path, &mutalyzer_proteins, true)?;
                } else {
                    println!("\nNo protein source available. Run 'prepare ferro --clinvar' first.");
                }
            }

            // Shard the cache for parallel workers
            if shards > 1 {
                println!("\n=== Sharding cache for {} workers ===", shards);
                shard_mutalyzer_cache(&cache_dir, shards)?;
            } else {
                // Generate single mutalyzer settings file with absolute path
                let settings_path = cache_dir.join("mutalyzer_settings.conf");
                let abs_cache_dir = cache_dir
                    .canonicalize()
                    .unwrap_or_else(|_| cache_dir.clone());
                let settings_content = format!(
                    "# Mutalyzer settings\n\
                     MUTALYZER_CACHE_DIR={}\n\
                     MUTALYZER_FILE_CACHE_ADD=True\n",
                    abs_cache_dir.display()
                );
                std::fs::write(&settings_path, settings_content).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to write settings: {}", e),
                    }
                })?;
            }

            println!("\n=== Mutalyzer preparation complete ===");
            println!("Cache directory: {}", cache_dir.display());
            if shards > 1 {
                println!("Shards: {} (use -j N where N <= {})", shards, shards);
                println!("Manifest: {}/shard_manifest.json", cache_dir.display());
            } else {
                println!(
                    "Settings file: {}/mutalyzer_settings.conf",
                    cache_dir.display()
                );
            }
            Ok(())
        }

        Tool::Biocommons | Tool::HgvsRs => {
            // Prerequisites: Docker must be available
            if !check_docker_available() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "Docker is required for biocommons/hgvs-rs setup.\n\n\
                          Install Docker: https://docs.docker.com/get-docker/"
                        .to_string(),
                });
            }

            let seqrepo = seqrepo_dir.ok_or_else(|| ferro_hgvs::FerroError::Io {
                msg: "SeqRepo directory is required.\n\n\
                      Use: ferro-benchmark prepare biocommons --seqrepo-dir <path>"
                    .to_string(),
            })?;

            // Print summary of what will be done
            if tool == Tool::Biocommons {
                print_biocommons_summary(seqrepo, uta_port);
            } else {
                print_hgvs_rs_summary(seqrepo, uta_port);
            }

            let tool_name = if tool == Tool::Biocommons {
                "biocommons"
            } else {
                "hgvs-rs"
            };

            // Step 1: Check/setup UTA
            println!("\n=== Step 1: UTA Database ===");
            let uta_url = format!(
                "postgresql://anonymous:anonymous@localhost:{}/uta/uta_20210129b",
                uta_port
            );

            // Check if UTA is already running
            if check_uta_connection(Some(&uta_url)) {
                println!("UTA database is already available at port {}", uta_port);
            } else if let Some(dump_path) = uta_dump {
                // Auto-setup UTA from provided dump file
                println!("UTA database not found. Setting up from dump file...");
                let config = BiocommonsLocalConfig {
                    uta_container_name: "ferro-uta".to_string(),
                    uta_image_tag: "uta_20210129b".to_string(),
                    uta_port,
                    ..Default::default()
                };
                match setup_uta(&config, force, Some(dump_path)) {
                    Ok(result) => {
                        println!("UTA setup complete!");
                        println!("  Container: {}", result.container_name);
                        println!("  Port: {}", result.port);
                    }
                    Err(e) => {
                        return Err(ferro_hgvs::FerroError::Io {
                            msg: format!("Failed to set up UTA: {}", e),
                        });
                    }
                }
            } else {
                println!("UTA database not found. Either:");
                println!("\n  Option A: Provide --uta-dump to auto-setup:");
                println!("    1. Download UTA dump (requires browser due to human verification):");
                println!("       https://dl.biocommons.org/uta/uta_20210129b.pgd.gz");
                println!("    2. Re-run with: ferro-benchmark prepare {} --seqrepo-dir {} --uta-dump <path>",
                         tool_name, seqrepo.display());
                println!("\n  Option B: Use setup command separately:");
                println!(
                    "    ferro-benchmark setup uta --port {} --uta-dump <path>",
                    uta_port
                );
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "UTA setup required - see instructions above".to_string(),
                });
            }

            // Step 2: Setup SeqRepo
            println!("\n=== Step 2: SeqRepo ===");
            let seqrepo_instance = seqrepo.join("2021-01-29");
            if seqrepo_instance.exists() && !force {
                println!("SeqRepo already exists at {}", seqrepo.display());
            } else {
                println!("Setting up SeqRepo at {}...", seqrepo.display());
                let config = BiocommonsLocalConfig {
                    uta_container_name: "ferro-uta".to_string(),
                    uta_image_tag: "uta_20210129b".to_string(),
                    uta_port,
                    seqrepo_dir: seqrepo.clone(),
                    seqrepo_instance: "2021-01-29".to_string(),
                };
                setup_seqrepo(&config, force)?;
            }

            // Step 3: Load additional sequences from ferro reference (if provided)
            if let Some(ferro_ref) = ferro_reference {
                println!("\n=== Step 3: Loading Additional Sequences ===");
                if no_load_transcripts {
                    println!("Loading genomic (NC_), RefSeqGene (NG_), and LRG sequences...");
                } else {
                    println!("Loading ALL sequences (transcripts, genomic, RefSeqGene, LRG)...");
                }
                load_sequences_to_seqrepo(&seqrepo_instance, ferro_ref, !no_load_transcripts)?;
            }

            // Auto-detect patterns from ferro reference if not explicitly provided
            let detected_patterns = if patterns.is_none() {
                ferro_reference.and_then(|fr| detect_patterns_from_ferro_reference(fr))
            } else {
                None
            };
            let effective_patterns: Option<&PathBuf> = patterns.or(detected_patterns.as_ref());
            if let Some(path) = detected_patterns.as_ref() {
                println!(
                    "Auto-detected patterns from ferro reference: {}",
                    path.display()
                );
            }

            // Step 3b: Fetch and load missing accessions from patterns (if provided or auto-detected)
            if let Some(patterns_path) = effective_patterns {
                println!("\n=== Step 3b: Fetching Missing Accessions from Patterns ===");
                match fetch_and_load_missing_accessions(
                    patterns_path,
                    &seqrepo_instance,
                    ferro_reference.as_ref().map(|p| p.as_path()),
                ) {
                    Ok(fetched) if fetched > 0 => {
                        println!("Fetched and loaded {} missing accessions", fetched);
                    }
                    Ok(_) => {
                        println!("All pattern accessions already in SeqRepo");
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to fetch missing accessions: {}", e);
                    }
                }
            }

            // Step 3c: Load cdot alignments into UTA (unless disabled)
            if !no_load_alignments {
                if let Some(ferro_ref) = ferro_reference {
                    let cdot_dir = ferro_ref.join("cdot");
                    // Find the JSON file (not .gz)
                    if let Ok(entries) = std::fs::read_dir(&cdot_dir) {
                        let cdot_file =
                            entries.filter_map(|e| e.ok()).map(|e| e.path()).find(|p| {
                                p.extension().map(|e| e == "json").unwrap_or(false)
                                    && !p.to_string_lossy().ends_with(".gz")
                            });

                        if let Some(cdot_path) = cdot_file {
                            println!(
                                "\n=== Step 3c: Loading Missing Transcript Alignments into UTA ==="
                            );
                            let config = ferro_hgvs::benchmark::UtaLoadConfig::default();

                            match ferro_hgvs::benchmark::load_cdot_alignments_to_uta(
                                &cdot_path, &config,
                            ) {
                                Ok(result) => {
                                    println!(
                                        "Loaded {} transcripts into UTA",
                                        result.transcripts_loaded
                                    );
                                    println!(
                                        "  Skipped (already exist): {}",
                                        result.transcripts_skipped
                                    );
                                    println!("  Exon sets created: {}", result.exon_sets_created);
                                    println!("  Exons created: {}", result.exons_created);
                                    println!("  Alignments created: {}", result.exon_alns_created);
                                    if !result.errors.is_empty() {
                                        eprintln!(
                                            "  Warnings: {} errors (see below)",
                                            result.errors.len()
                                        );
                                        for err in result.errors.iter().take(5) {
                                            eprintln!("    {}", err);
                                        }
                                        if result.errors.len() > 5 {
                                            eprintln!(
                                                "    ... and {} more",
                                                result.errors.len() - 5
                                            );
                                        }
                                    }
                                }
                                Err(e) => {
                                    eprintln!("Warning: Failed to load cdot alignments: {}", e);
                                    eprintln!("  hgvs-rs may fail on newer transcripts not in UTA snapshot");
                                }
                            }
                        }
                    }
                }
            }

            // Step 4: Generate settings file
            // Write settings to ferro reference directory if provided, otherwise output_dir
            println!("\n=== Step 4: Settings File ===");
            let settings_dir = match ferro_reference {
                Some(fr) => fr.clone(),
                None => output_dir.to_path_buf(),
            };
            std::fs::create_dir_all(&settings_dir).map_err(|e| ferro_hgvs::FerroError::Io {
                msg: format!("Failed to create settings directory: {}", e),
            })?;
            let settings_path = settings_dir.join(format!("{}_settings.conf", tool_name));
            write_biocommons_settings(&settings_path, &uta_url, &seqrepo_instance)?;
            println!("Settings written to {}", settings_path.display());

            // Step 5: Verify setup
            println!("\n=== Verification ===");
            if check_seqrepo(&seqrepo_instance) {
                println!("✓ SeqRepo is accessible");
            } else {
                println!("✗ SeqRepo verification failed");
            }

            if check_uta_connection(Some(&uta_url)) {
                println!("✓ UTA database is accessible");
            } else {
                println!("✗ UTA database verification failed");
            }

            println!("\n=== {} preparation complete ===", tool_name);
            println!("Settings file: {}", settings_path.display());
            println!(
                "\nTo normalize: ferro-benchmark normalize {} -i patterns.txt -o results.json \\",
                tool_name
            );
            if tool == Tool::Biocommons {
                println!(
                    "              --biocommons-settings {}",
                    settings_path.display()
                );
            } else {
                println!(
                    "              --seqrepo-path {}",
                    seqrepo_instance.display()
                );
            }

            Ok(())
        }

        Tool::All => {
            // Prepare all tools in sequence: ferro → mutalyzer → biocommons → hgvs-rs
            println!("=== Preparing All Tools ===\n");

            // Determine paths
            let ferro_output = output_dir.join("ferro");
            let mutalyzer_output = output_dir.join("mutalyzer");
            let seqrepo_dir_path = seqrepo_dir
                .cloned()
                .unwrap_or_else(|| output_dir.join("seqrepo"));

            // Step 1: Prepare ferro
            println!("╔════════════════════════════════════════════════════════════════╗");
            println!("║  Step 1/4: Preparing ferro                                     ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");
            run_prepare_tool(
                Tool::Ferro,
                &ferro_output,
                genome.clone(),
                no_refseqgene,
                no_lrg,
                force,
                None,
                patterns,
                None,
                true, // no_proteins not relevant for ferro
                1,    // shards not relevant for ferro
                uta_port,
                None,
                None,
                true, // no_load_transcripts not relevant for ferro
                true, // no_load_alignments not relevant for ferro
            )?;
            println!("\n✓ Ferro preparation complete\n");

            // Step 2: Prepare mutalyzer
            println!("╔════════════════════════════════════════════════════════════════╗");
            println!("║  Step 2/4: Preparing mutalyzer                                 ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");
            run_prepare_tool(
                Tool::Mutalyzer,
                &mutalyzer_output,
                genome.clone(),
                no_refseqgene,
                no_lrg,
                force,
                Some(&ferro_output),
                patterns,
                clinvar,
                no_proteins,
                shards, // use shards parameter for mutalyzer
                uta_port,
                None,
                None,
                true, // no_load_transcripts not relevant for mutalyzer
                true, // no_load_alignments not relevant for mutalyzer
            )?;
            println!("\n✓ Mutalyzer preparation complete\n");

            // Step 3: Prepare biocommons (also sets up UTA and SeqRepo)
            println!("╔════════════════════════════════════════════════════════════════╗");
            println!("║  Step 3/4: Preparing biocommons                                ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");
            run_prepare_tool(
                Tool::Biocommons,
                output_dir,
                genome.clone(),
                no_refseqgene,
                no_lrg,
                force,
                Some(&ferro_output),
                patterns,
                None,
                true,
                1, // shards not relevant for biocommons
                uta_port,
                Some(&seqrepo_dir_path),
                uta_dump,
                no_load_transcripts, // Pass through from parent command
                no_load_alignments,  // Pass through from parent command
            )?;
            println!("\n✓ Biocommons preparation complete\n");

            // Step 4: Prepare hgvs-rs (shares UTA and SeqRepo with biocommons)
            println!("╔════════════════════════════════════════════════════════════════╗");
            println!("║  Step 4/4: Preparing hgvs-rs                                   ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");
            run_prepare_tool(
                Tool::HgvsRs,
                output_dir,
                genome,
                no_refseqgene,
                no_lrg,
                force,
                Some(&ferro_output),
                patterns,
                None,
                true,
                1, // shards not relevant for hgvs-rs
                uta_port,
                Some(&seqrepo_dir_path),
                None, // UTA already set up by biocommons
                true, // transcripts already loaded by biocommons step
                true, // alignments already loaded by biocommons step
            )?;
            println!("\n✓ hgvs-rs preparation complete\n");

            // Summary
            println!("╔════════════════════════════════════════════════════════════════╗");
            println!("║  All Tools Prepared Successfully!                              ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");
            println!("Tool outputs:");
            println!("  ferro:      {}", ferro_output.display());
            println!("  mutalyzer:  {}", mutalyzer_output.display());
            println!(
                "  biocommons: {}/biocommons_settings.conf",
                output_dir.display()
            );
            println!(
                "  hgvs-rs:    {}/hgvs-rs_settings.conf",
                output_dir.display()
            );
            println!(
                "\nTo check all tools: ferro-benchmark check all --reference {} \\",
                output_dir.display()
            );
            println!(
                "                      --mutalyzer-settings {}/mutalyzer_settings.conf \\",
                mutalyzer_output.display()
            );
            println!(
                "                      --seqrepo-path {}/2021-01-29",
                seqrepo_dir_path.display()
            );

            Ok(())
        }
    }
}

/// Parse HGVS patterns with a specific tool
fn run_parse_tool(
    tool: Tool,
    input: &Path,
    output: &Path,
    timing: Option<&PathBuf>,
) -> Result<(), ferro_hgvs::FerroError> {
    // Prerequisite: input file must exist
    if !input.exists() {
        return Err(ferro_hgvs::FerroError::Io {
            msg: format!(
                "Input file not found: {}\n\n\
                 Create a patterns file with one HGVS expression per line, or extract from ClinVar:\n\
                 ferro-benchmark extract-clinvar -i clinvar/hgvs4variation.txt.gz -o patterns.txt",
                input.display()
            ),
        });
    }

    match tool {
        Tool::Ferro => {
            use ferro_hgvs::benchmark::parse_ferro_unified;
            println!("Parsing patterns with ferro...");
            let result = parse_ferro_unified(input, output)?;
            println!(
                "Parsed {} patterns in {:.2}s",
                result.total_patterns, result.elapsed_seconds
            );
            println!(
                "Successful: {} ({:.1}%)",
                result.successful,
                100.0 * result.successful as f64 / result.total_patterns as f64
            );
            Ok(())
        }

        Tool::Mutalyzer => {
            // Check if mutalyzer parser is available
            if !has_mutalyzer_parser() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "Mutalyzer parser is not available.\n\n\
                          Install with: pixi global install mutalyzer-hgvs-parser\n\
                          Or: pip install mutalyzer-hgvs-parser"
                        .to_string(),
                });
            }

            println!("Parsing patterns with mutalyzer...");
            run_mutalyzer_parser_subprocess(input.to_str().unwrap(), output.to_str().unwrap())?;

            // Read results to show summary
            let results: serde_json::Value =
                serde_json::from_str(&std::fs::read_to_string(output).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to read results: {}", e),
                    }
                })?)
                .map_err(|e| ferro_hgvs::FerroError::Json {
                    msg: format!("Failed to parse results: {}", e),
                })?;

            let total = results["total_patterns"].as_u64().unwrap_or(0) as usize;
            let successful = results["successful"].as_u64().unwrap_or(0) as usize;
            let elapsed = results["elapsed_seconds"].as_f64().unwrap_or(0.0);

            println!("Parsed {} patterns in {:.2}s", total, elapsed);
            println!(
                "Successful: {} ({:.1}%)",
                successful,
                100.0 * successful as f64 / total as f64
            );
            Ok(())
        }

        Tool::Biocommons => {
            // biocommons uses its own parser internally, no separate parse command
            // We can still parse by running the normalizer with a null-op
            if !has_biocommons_normalizer() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "biocommons/hgvs is not available.\n\n\
                          Install with: pip install hgvs"
                        .to_string(),
                });
            }

            println!("Parsing patterns with biocommons...");
            println!("Note: biocommons parsing is done via its normalizer (parse-only mode)");

            // Python script - outputs unified ToolParseOutput format
            let script = r#"
import sys
import json
import time
from hgvs.parser import Parser

parser = Parser()
patterns = [line.strip() for line in open(sys.argv[1]) if line.strip()]

results = []
successful = 0
start = time.time()

for pattern in patterns:
    try:
        parsed = parser.parse_hgvs_variant(pattern)
        results.append({
            "input": pattern,
            "success": True,
            "output": str(parsed),
            "error": None
        })
        successful += 1
    except Exception as e:
        results.append({
            "input": pattern,
            "success": False,
            "output": None,
            "error": str(e)[:200]
        })

elapsed = time.time() - start

# Unified output format matching ToolParseOutput
output = {
    "tool": "biocommons",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": len(patterns) - successful,
    "elapsed_seconds": elapsed,
    "throughput": len(patterns) / elapsed if elapsed > 0 else 0,
    "results": results
}

with open(sys.argv[2], 'w') as f:
    json.dump(output, f, indent=2)

print(f"Parsed {len(patterns)} patterns in {elapsed:.2f}s")
print(f"Successful: {successful} ({100*successful/len(patterns):.1f}%)")
"#;

            let script_path = std::env::temp_dir().join("biocommons_parse.py");
            std::fs::write(&script_path, script).map_err(|e| ferro_hgvs::FerroError::Io {
                msg: format!("Failed to write script: {}", e),
            })?;

            let status = std::process::Command::new("python3")
                .arg(&script_path)
                .arg(input)
                .arg(output)
                .status()
                .map_err(|e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to run biocommons parser: {}", e),
                })?;

            if !status.success() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "biocommons parser failed".to_string(),
                });
            }

            // Generate timing file
            if let Some(timing_path) = timing {
                let results: serde_json::Value =
                    serde_json::from_str(&std::fs::read_to_string(output).unwrap_or_default())
                        .unwrap_or_default();

                let total = results["total_patterns"].as_u64().unwrap_or(0) as usize;
                let successful = results["successful"].as_u64().unwrap_or(0) as usize;
                let elapsed = std::time::Duration::from_secs_f64(
                    results["elapsed_seconds"].as_f64().unwrap_or(0.0),
                );

                let timing_info = TimingInfo::new("biocommons", total, successful, elapsed);
                std::fs::write(
                    timing_path,
                    serde_json::to_string_pretty(&timing_info).unwrap(),
                )
                .ok();
            }

            Ok(())
        }

        Tool::HgvsRs => {
            #[cfg(not(feature = "hgvs-rs"))]
            {
                Err(ferro_hgvs::FerroError::Io {
                    msg: "hgvs-rs feature not enabled.\n\n\
                          Rebuild with: cargo build --release --features benchmark,hgvs-rs"
                        .to_string(),
                })?
            }

            #[cfg(feature = "hgvs-rs")]
            {
                println!("Parsing patterns with hgvs-rs...");

                // hgvs-rs uses HgvsVariant::from_str for parsing
                use hgvs::parser::HgvsVariant;
                use std::io::{BufRead, BufReader};
                use std::str::FromStr;
                use std::time::Instant;

                let file = std::fs::File::open(input).map_err(|e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to open input: {}", e),
                })?;
                let reader = BufReader::new(file);
                let patterns: Vec<String> = reader.lines().map_while(Result::ok).collect();

                let mut results = Vec::new();
                let mut successful = 0;

                let start = Instant::now();
                for pattern in &patterns {
                    match HgvsVariant::from_str(pattern) {
                        Ok(parsed) => {
                            results.push(serde_json::json!({
                                "input": pattern,
                                "success": true,
                                "output": format!("{}", parsed),
                                "error": null
                            }));
                            successful += 1;
                        }
                        Err(e) => {
                            results.push(serde_json::json!({
                                "input": pattern,
                                "success": false,
                                "output": null,
                                "error": format!("{}", e)
                            }));
                        }
                    }
                }
                let elapsed = start.elapsed();
                let elapsed_secs = elapsed.as_secs_f64();
                let throughput = if elapsed_secs > 0.0 {
                    patterns.len() as f64 / elapsed_secs
                } else {
                    0.0
                };

                // Unified output format matching ToolParseOutput
                let output_json = serde_json::json!({
                    "tool": "hgvs-rs",
                    "total_patterns": patterns.len(),
                    "successful": successful,
                    "failed": patterns.len() - successful,
                    "elapsed_seconds": elapsed_secs,
                    "throughput": throughput,
                    "results": results
                });

                std::fs::write(output, serde_json::to_string_pretty(&output_json).unwrap())
                    .map_err(|e| ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to write output: {}", e),
                    })?;

                if let Some(timing_path) = timing {
                    let timing_info =
                        TimingInfo::new("hgvs-rs", patterns.len(), successful, elapsed);
                    std::fs::write(
                        timing_path,
                        serde_json::to_string_pretty(&timing_info).unwrap(),
                    )
                    .map_err(|e| ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to write timing: {}", e),
                    })?;
                }

                println!(
                    "Parsed {} patterns in {:.2}s",
                    patterns.len(),
                    elapsed.as_secs_f64()
                );
                println!(
                    "Successful: {} ({:.1}%)",
                    successful,
                    100.0 * successful as f64 / patterns.len() as f64
                );
                Ok(())
            }
        }

        Tool::All => Err(ferro_hgvs::FerroError::Io {
            msg: "'all' is not supported for parsing.\n\n\
                      Use individual tools: ferro-benchmark parse ferro/mutalyzer/biocommons/hgvs-rs"
                .to_string(),
        }),
    }
}

/// Normalize HGVS patterns with a specific tool
#[allow(clippy::too_many_arguments)]
fn run_normalize_tool(
    tool: Tool,
    input: &Path,
    output: &Path,
    timing: Option<&PathBuf>,
    reference: Option<&PathBuf>,
    mutalyzer_settings: Option<&PathBuf>,
    biocommons_settings: Option<&PathBuf>,
    lrg_mapping: Option<&PathBuf>,
    seqrepo_path: Option<&PathBuf>,
    workers: usize,
    no_rewrite_intronic: bool,
    allow_network: bool,
) -> Result<(), ferro_hgvs::FerroError> {
    // Prerequisite: input file must exist
    if !input.exists() {
        return Err(ferro_hgvs::FerroError::Io {
            msg: format!(
                "Input file not found: {}\n\n\
                 Create a patterns file with one HGVS expression per line, or extract from ClinVar:\n\
                 ferro-benchmark extract-clinvar -i clinvar/hgvs4variation.txt.gz -o patterns.txt",
                input.display()
            ),
        });
    }

    let timing_path = timing
        .cloned()
        .unwrap_or_else(|| output.with_extension("timing.json"));

    match tool {
        Tool::Ferro => {
            // Prerequisite: reference data must exist
            let ref_dir = reference.ok_or_else(|| ferro_hgvs::FerroError::Io {
                msg: "Ferro normalization requires reference data.\n\n\
                      First run: ferro prepare --output-dir benchmark-output\n\
                      Or:        ferro-benchmark prepare ferro --output-dir benchmark-output\n\n\
                      Then run:  ferro normalize -i patterns.txt --reference benchmark-output\n\
                      Or:        ferro-benchmark normalize ferro -i patterns.txt -o results.json --reference benchmark-output"
                    .to_string(),
            })?;

            let manifest_path = ref_dir.join("manifest.json");
            if !manifest_path.exists() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: format!(
                        "Ferro reference data not found at {}\n\n\
                         First run: ferro prepare --output-dir {}\n\
                         Or:        ferro-benchmark prepare ferro --output-dir {}",
                        manifest_path.display(),
                        ref_dir.display(),
                        ref_dir.display()
                    ),
                });
            }

            println!("Normalizing patterns with ferro...");
            if workers > 1 {
                println!(
                    "Note: Parallel ferro normalization not yet implemented, using single worker"
                );
            }
            normalize_ferro(input, output, &timing_path, Some(ref_dir)).map(|_| ())
        }

        Tool::Mutalyzer => {
            // Prerequisite: mutalyzer settings must exist
            let settings = mutalyzer_settings.ok_or_else(|| ferro_hgvs::FerroError::Io {
                msg: "Mutalyzer normalization requires settings file.\n\n\
                      First run: ferro-benchmark prepare mutalyzer --ferro-reference <ferro-dir> --output-dir mutalyzer_cache\n\
                      Then run:  ferro-benchmark normalize mutalyzer -i patterns.txt -o results.json \\
                                 --mutalyzer-settings mutalyzer_cache/mutalyzer_settings.conf".to_string(),
            })?;

            if !settings.exists() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: format!(
                        "Mutalyzer settings not found at {}\n\n\
                         First run: ferro-benchmark prepare mutalyzer --ferro-reference <ferro-dir>",
                        settings.display()
                    ),
                });
            }

            // Check if mutalyzer normalizer is available
            if !has_mutalyzer_normalizer() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "Mutalyzer normalizer is not available.\n\n\
                          Install with: pixi global install mutalyzer-hgvs-parser mutalyzer-normalizer\n\
                          Or: pip install mutalyzer-hgvs-parser mutalyzer-normalizer".to_string(),
                });
            }

            // Load patterns from input file
            let patterns: Vec<String> = std::fs::read_to_string(input)
                .map_err(|e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to read input file: {}", e),
                })?
                .lines()
                .map(|l| l.trim().to_string())
                .filter(|l| !l.is_empty())
                .collect();

            // Auto-load transcript mapping for intronic rewriting unless --no-rewrite-intronic
            let patterns_to_normalize = if no_rewrite_intronic {
                patterns.clone()
            } else {
                // Look for mapping file in same directory as settings
                let mapping_path = settings.parent().map(|p| p.join("transcript_mapping.json"));

                if let Some(ref path) = mapping_path {
                    if path.exists() {
                        match load_transcript_mapping(path) {
                            Ok(mapping) => {
                                eprintln!(
                                    "Loaded {} transcript mappings for intronic rewriting",
                                    mapping.len()
                                );
                                let rewritten = rewrite_with_genomic_context(&patterns, &mapping);
                                let rewritten_count = patterns
                                    .iter()
                                    .zip(rewritten.iter())
                                    .filter(|(orig, new)| orig != new)
                                    .count();
                                if rewritten_count > 0 {
                                    eprintln!(
                                        "Rewrote {} intronic variants with genomic context",
                                        rewritten_count
                                    );
                                }
                                rewritten
                            }
                            Err(e) => {
                                eprintln!("Warning: Failed to load transcript mapping: {}", e);
                                patterns.clone()
                            }
                        }
                    } else {
                        patterns.clone()
                    }
                } else {
                    patterns.clone()
                }
            };

            let effective_workers = if workers > 1 { workers } else { 1 };
            println!(
                "Normalizing {} patterns with mutalyzer ({} workers)...",
                patterns.len(),
                effective_workers
            );

            // Use parallel normalization
            let (results, elapsed, _error_counts) = run_mutalyzer_normalize_parallel(
                &patterns_to_normalize,
                effective_workers,
                Some(settings.to_str().unwrap()),
                allow_network,
            )?;

            let total = patterns.len();
            let successful = results.iter().filter(|r| r.success).count();
            let failed = total - successful;

            // Write results in the same format as single-threaded version
            let output_data = serde_json::json!({
                "tool": "mutalyzer",
                "total_patterns": total,
                "successful": successful,
                "failed": failed,
                "elapsed_seconds": elapsed.as_secs_f64(),
                "patterns_per_second": if elapsed.as_secs_f64() > 0.0 {
                    total as f64 / elapsed.as_secs_f64()
                } else {
                    0.0
                },
                "network_calls": 0,
                "network_stats": [],
                "results": results.iter().map(|r| {
                    if r.success {
                        serde_json::json!({
                            "input": r.input,
                            "success": true,
                            "output": r.output.as_deref().unwrap_or("")
                        })
                    } else {
                        serde_json::json!({
                            "input": r.input,
                            "success": false,
                            "error": r.error.as_deref().unwrap_or("Unknown error")
                        })
                    }
                }).collect::<Vec<_>>()
            });

            std::fs::write(output, serde_json::to_string_pretty(&output_data).unwrap()).map_err(
                |e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to write results: {}", e),
                },
            )?;

            let timing_info = TimingInfo::new("mutalyzer", total, successful, elapsed);
            std::fs::write(
                &timing_path,
                serde_json::to_string_pretty(&timing_info).unwrap(),
            )
            .map_err(|e| ferro_hgvs::FerroError::Io {
                msg: format!("Failed to write timing: {}", e),
            })?;

            println!(
                "Completed: {}/{} successful ({:.1}%) in {:.1}s ({:.0} patterns/sec)",
                successful,
                total,
                100.0 * successful as f64 / total as f64,
                elapsed.as_secs_f64(),
                total as f64 / elapsed.as_secs_f64()
            );
            println!("Results: {}", output.display());
            println!("Timing: {}", timing_path.display());
            Ok(())
        }

        Tool::Biocommons => {
            // Check if biocommons is available
            if !has_biocommons_normalizer() {
                return Err(ferro_hgvs::FerroError::Io {
                    msg: "biocommons/hgvs is not available.\n\n\
                          Install with: pip install hgvs"
                        .to_string(),
                });
            }

            // Get UTA URL and SeqRepo dir from settings or default
            let (uta_db_url, seqrepo_dir) = if let Some(settings) = biocommons_settings {
                if !settings.exists() {
                    return Err(ferro_hgvs::FerroError::Io {
                        msg: format!(
                            "Biocommons settings not found at {}\n\n\
                             First run: ferro-benchmark prepare biocommons --seqrepo-dir <path>",
                            settings.display()
                        ),
                    });
                }
                let settings_content = load_biocommons_settings(settings)?;
                (
                    Some(settings_content.uta_db_url),
                    Some(settings_content.seqrepo_dir.to_string_lossy().to_string()),
                )
            } else {
                (None, None)
            };

            // Get LRG mapping file path if provided
            let lrg_mapping_str = lrg_mapping.map(|p| p.to_string_lossy().to_string());

            println!("Normalizing patterns with biocommons...");
            run_normalize_biocommons(
                input,
                output,
                &timing_path,
                uta_db_url.as_deref(),
                seqrepo_dir.as_deref(),
                lrg_mapping_str.as_deref(),
                workers,
            )
        }

        Tool::HgvsRs => {
            #[cfg(not(feature = "hgvs-rs"))]
            {
                let _ = seqrepo_path; // Silence unused variable warning
                Err(ferro_hgvs::FerroError::Io {
                    msg: "hgvs-rs feature not enabled.\n\n\
                          Rebuild with: cargo build --release --features benchmark,hgvs-rs"
                        .to_string(),
                })?
            }

            #[cfg(feature = "hgvs-rs")]
            {
                let seqrepo = seqrepo_path.ok_or_else(|| ferro_hgvs::FerroError::Io {
                    msg: "hgvs-rs normalization requires SeqRepo path.\n\n\
                          First run: ferro-benchmark prepare hgvs-rs --seqrepo-dir <path>\n\
                          Then run:  ferro-benchmark normalize hgvs-rs -i patterns.txt -o results.json \\
                                     --seqrepo-path <path>/2021-01-29".to_string(),
                })?;

                if !seqrepo.exists() {
                    return Err(ferro_hgvs::FerroError::Io {
                        msg: format!(
                            "SeqRepo not found at {}\n\n\
                             First run: ferro-benchmark prepare hgvs-rs --seqrepo-dir <parent-dir>",
                            seqrepo.display()
                        ),
                    });
                }

                // Default UTA configuration
                let uta_db_url = "postgresql://anonymous:anonymous@localhost:5432/uta";
                let uta_db_schema = "uta_20210129b";

                println!("Normalizing patterns with hgvs-rs...");
                run_normalize_hgvs_rs(
                    input,
                    output,
                    &timing_path,
                    uta_db_url,
                    uta_db_schema,
                    seqrepo.to_str().unwrap(),
                    lrg_mapping.as_ref().map(|p| p.as_path()),
                    workers,
                )
            }
        }

        Tool::All => {
            Err(ferro_hgvs::FerroError::Io {
                msg: "'all' is not supported for normalization.\n\n\
                      Use individual tools: ferro-benchmark normalize ferro/mutalyzer/biocommons/hgvs-rs"
                    .to_string(),
            })
        }
    }
}

/// Compute and print overlap tables and upset plot data for multi-tool comparison
fn print_multiway_overlap_tables(
    tool_names: &[String],
    success_sets: &[std::collections::HashSet<String>],
    all_patterns: &std::collections::HashSet<String>,
) {
    use std::collections::HashSet;

    // Build failure sets
    let fail_sets: Vec<HashSet<String>> = success_sets
        .iter()
        .map(|s| all_patterns.difference(s).cloned().collect())
        .collect();

    // Print success overlap table
    println!("\n{}", "=".repeat(70));
    println!("SUCCESS OVERLAP TABLE");
    println!("{}", "=".repeat(70));
    println!();

    // Header
    print!("{:<45} |", "Condition");
    for name in tool_names {
        print!(" {:>10}", name);
    }
    println!();
    println!("{}", "-".repeat(45 + 2 + tool_names.len() * 11));

    // Row 1: first tool succeeded (ferro if present)
    let first_tool = &tool_names[0];
    let base = &success_sets[0];
    print!("{:<45} |", format!("{} succeeded", first_tool));
    for ss in success_sets.iter() {
        print!(" {:>10}", base.intersection(ss).count());
    }
    println!();

    // Additional rows: first tool AND each other tool succeeded
    for (i, name) in tool_names.iter().enumerate().skip(1) {
        let base: HashSet<String> = success_sets[0]
            .intersection(&success_sets[i])
            .cloned()
            .collect();
        print!("{:<45} |", format!("{} AND {} succeeded", first_tool, name));
        for ss in success_sets.iter() {
            print!(" {:>10}", base.intersection(ss).count());
        }
        println!();
    }

    // Print failure overlap table
    println!();
    println!("{}", "=".repeat(70));
    println!("FAILURE OVERLAP TABLE");
    println!("{}", "=".repeat(70));
    println!();

    // Header
    print!("{:<45} |", "Condition");
    for name in tool_names {
        print!(" {:>10}", name);
    }
    println!();
    println!("{}", "-".repeat(45 + 2 + tool_names.len() * 11));

    // Row 1: first tool failed
    let base = &fail_sets[0];
    print!("{:<45} |", format!("{} failed", first_tool));
    for fs in fail_sets.iter() {
        print!(" {:>10}", base.intersection(fs).count());
    }
    println!();

    // Additional rows: first tool AND each other tool failed
    for (i, name) in tool_names.iter().enumerate().skip(1) {
        let base: HashSet<String> = fail_sets[0].intersection(&fail_sets[i]).cloned().collect();
        print!("{:<45} |", format!("{} AND {} failed", first_tool, name));
        for fs in fail_sets.iter() {
            print!(" {:>10}", base.intersection(fs).count());
        }
        println!();
    }

    // Print upset plot data for success
    println!();
    println!("{}", "=".repeat(70));
    println!("UPSET PLOT DATA - SUCCESS (exclusive sets)");
    println!("{}", "=".repeat(70));
    println!();
    println!("{:<50} | Count", "Combination (ONLY these succeeded)");
    println!("{}", "-".repeat(60));

    let mut upset_success: Vec<(Vec<usize>, usize)> = Vec::new();

    // Generate all combinations
    let n = tool_names.len();
    for mask in 1..(1 << n) {
        let combo: Vec<usize> = (0..n).filter(|&i| (mask >> i) & 1 == 1).collect();

        // Patterns that succeeded in exactly these tools
        let mut in_these: HashSet<String> = all_patterns.clone();
        for &i in &combo {
            in_these = in_these.intersection(&success_sets[i]).cloned().collect();
        }
        for (i, ss) in success_sets.iter().enumerate() {
            if !combo.contains(&i) {
                in_these = in_these.difference(ss).cloned().collect();
            }
        }

        if !in_these.is_empty() {
            upset_success.push((combo, in_these.len()));
        }
    }

    // Check for "none succeeded"
    let any_success: HashSet<String> = success_sets
        .iter()
        .fold(HashSet::new(), |acc, s| acc.union(s).cloned().collect());
    let none_succeeded: HashSet<String> = all_patterns.difference(&any_success).cloned().collect();
    if !none_succeeded.is_empty() {
        upset_success.push((vec![], none_succeeded.len()));
    }

    // Sort by count descending
    upset_success.sort_by(|a, b| b.1.cmp(&a.1));

    for (combo, count) in &upset_success {
        let label = if combo.is_empty() {
            "none".to_string()
        } else {
            combo
                .iter()
                .map(|&i| tool_names[i].as_str())
                .collect::<Vec<_>>()
                .join(" & ")
        };
        println!("{:<50} | {:>5}", label, count);
    }

    // Print upset plot data for failure
    println!();
    println!("{}", "=".repeat(70));
    println!("UPSET PLOT DATA - FAILURE (exclusive sets)");
    println!("{}", "=".repeat(70));
    println!();
    println!("{:<50} | Count", "Combination (ONLY these failed)");
    println!("{}", "-".repeat(60));

    let mut upset_fail: Vec<(Vec<usize>, usize)> = Vec::new();

    for mask in 1..(1 << n) {
        let combo: Vec<usize> = (0..n).filter(|&i| (mask >> i) & 1 == 1).collect();

        let mut in_these: HashSet<String> = all_patterns.clone();
        for &i in &combo {
            in_these = in_these.intersection(&fail_sets[i]).cloned().collect();
        }
        for (i, fs) in fail_sets.iter().enumerate() {
            if !combo.contains(&i) {
                in_these = in_these.difference(fs).cloned().collect();
            }
        }

        if !in_these.is_empty() {
            upset_fail.push((combo, in_these.len()));
        }
    }

    // Check for "none failed" (all succeeded)
    let any_fail: HashSet<String> = fail_sets
        .iter()
        .fold(HashSet::new(), |acc, s| acc.union(s).cloned().collect());
    let none_failed: HashSet<String> = all_patterns.difference(&any_fail).cloned().collect();
    if !none_failed.is_empty() {
        upset_fail.push((vec![], none_failed.len()));
    }

    upset_fail.sort_by(|a, b| b.1.cmp(&a.1));

    for (combo, count) in &upset_fail {
        let label = if combo.is_empty() {
            "none (all succeeded)".to_string()
        } else {
            combo
                .iter()
                .map(|&i| tool_names[i].as_str())
                .collect::<Vec<_>>()
                .join(" & ")
        };
        println!("{:<50} | {:>5}", label, count);
    }

    // Print summary counts
    println!();
    println!("{}", "=".repeat(70));
    println!("SUMMARY COUNTS");
    println!("{}", "=".repeat(70));
    println!("Total patterns: {}", all_patterns.len());
    for (i, name) in tool_names.iter().enumerate() {
        println!(
            "{}: success={}, fail={}",
            name,
            success_sets[i].len(),
            fail_sets[i].len()
        );
    }
}

/// Compare results from multiple tools
fn run_compare_tool(
    mode: CompareMode,
    results: &[PathBuf],
    output: &Path,
    equivalence: bool,
) -> Result<(), ferro_hgvs::FerroError> {
    // Prerequisite: at least 2 result files
    if results.len() < 2 {
        return Err(ferro_hgvs::FerroError::Io {
            msg: "At least 2 result files are required for comparison.\n\n\
                  Usage: ferro-benchmark compare <parse|normalize> file1.json file2.json [file3.json ...]".to_string(),
        });
    }

    // Check all files exist
    for result_file in results {
        if !result_file.exists() {
            return Err(ferro_hgvs::FerroError::Io {
                msg: format!("Result file not found: {}", result_file.display()),
            });
        }
    }

    match mode {
        CompareMode::Parse => {
            println!("Comparing parse results from {} files...", results.len());

            // Load all results
            let mut all_results: Vec<(String, serde_json::Value)> = Vec::new();
            for result_file in results {
                let content = std::fs::read_to_string(result_file).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to read {}: {}", result_file.display(), e),
                    }
                })?;
                let data: serde_json::Value =
                    serde_json::from_str(&content).map_err(|e| ferro_hgvs::FerroError::Json {
                        msg: format!("Failed to parse {}: {}", result_file.display(), e),
                    })?;
                let tool_name = data["tool"].as_str().unwrap_or("unknown").to_string();
                all_results.push((tool_name, data));
            }

            let tool_names: Vec<String> =
                all_results.iter().map(|(name, _)| name.clone()).collect();

            // Compare results
            let mut comparison = serde_json::json!({
                "mode": "parse",
                "files": results.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
                "tools": tool_names.clone(),
            });

            // Extract statistics from each tool
            let mut stats = Vec::new();
            for (tool_name, data) in &all_results {
                stats.push(serde_json::json!({
                    "tool": tool_name,
                    "total": data["total_patterns"],
                    "successful": data["successful"],
                    "failed": data["failed"],
                    "elapsed_seconds": data["elapsed_seconds"],
                }));
            }
            comparison["statistics"] = serde_json::Value::Array(stats);

            // Extract result arrays from each tool
            let result_arrays: Vec<Option<&Vec<serde_json::Value>>> = all_results
                .iter()
                .map(|(_, data)| data["results"].as_array())
                .collect();

            // Check if all tools have results arrays
            if result_arrays.iter().all(|r| r.is_some()) {
                let arrays: Vec<&Vec<serde_json::Value>> =
                    result_arrays.into_iter().map(|r| r.unwrap()).collect();

                // Check all arrays have same length
                let lengths: Vec<usize> = arrays.iter().map(|a| a.len()).collect();
                if lengths.iter().all(|&len| len == lengths[0]) {
                    let total_patterns = lengths[0];

                    // Build success sets for overlap tables
                    let mut success_sets: Vec<std::collections::HashSet<String>> = Vec::new();
                    let mut all_patterns: std::collections::HashSet<String> =
                        std::collections::HashSet::new();

                    for arr in &arrays {
                        let mut success_set = std::collections::HashSet::new();
                        for item in arr.iter() {
                            let input = item["input"].as_str().unwrap_or("").to_string();
                            all_patterns.insert(input.clone());
                            if item["success"].as_bool().unwrap_or(false) {
                                success_set.insert(input);
                            }
                        }
                        success_sets.push(success_set);
                    }

                    // Compare outputs across all tools
                    let mut all_succeed = 0;
                    let mut all_agree = 0;
                    let mut disagreements: Vec<serde_json::Value> = Vec::new();

                    for i in 0..total_patterns {
                        // Check if all tools succeeded for this pattern
                        let all_success = arrays
                            .iter()
                            .all(|arr| arr[i]["success"].as_bool().unwrap_or(false));

                        if all_success {
                            all_succeed += 1;

                            // Get outputs from all tools
                            let outputs: Vec<&str> = arrays
                                .iter()
                                .map(|arr| arr[i]["output"].as_str().unwrap_or(""))
                                .collect();

                            // Check if all outputs are the same
                            let first_output = outputs[0];
                            let all_same = outputs.iter().all(|&o| o == first_output);

                            if all_same {
                                all_agree += 1;
                            } else {
                                // Record disagreement
                                let input = arrays[0][i]["input"].as_str().unwrap_or("");
                                let mut disagreement = serde_json::json!({
                                    "input": input,
                                });
                                for (j, tool_name) in tool_names.iter().enumerate() {
                                    disagreement[tool_name] =
                                        serde_json::Value::String(outputs[j].to_string());
                                }
                                disagreements.push(disagreement);
                            }
                        }
                    }

                    // Add output agreement stats to comparison
                    comparison["output_agreement"] = serde_json::json!({
                        "all_tools_succeed": all_succeed,
                        "all_tools_agree": all_agree,
                        "disagreements_count": disagreements.len(),
                        "agreement_rate": if all_succeed > 0 {
                            all_agree as f64 / all_succeed as f64
                        } else {
                            0.0
                        },
                        "disagreements": disagreements,
                    });

                    println!(
                        "Output agreement: {}/{} ({:.1}%)",
                        all_agree,
                        all_succeed,
                        if all_succeed > 0 {
                            100.0 * all_agree as f64 / all_succeed as f64
                        } else {
                            0.0
                        }
                    );

                    // Print overlap tables and upset plot data
                    print_multiway_overlap_tables(&tool_names, &success_sets, &all_patterns);
                }
            }

            // Write comparison
            std::fs::write(output, serde_json::to_string_pretty(&comparison).unwrap()).map_err(
                |e| ferro_hgvs::FerroError::Io {
                    msg: format!("Failed to write comparison: {}", e),
                },
            )?;

            println!("Comparison written to {}", output.display());
            Ok(())
        }

        CompareMode::Normalize => {
            println!(
                "Comparing normalize results from {} files...",
                results.len()
            );

            if equivalence {
                println!("Equivalence checking enabled");
            }

            // For now, do pairwise comparison
            // In the future, we could do multi-way comparison
            if results.len() == 2 {
                // Load both result files and compare
                let result1 = std::fs::read_to_string(&results[0]).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to read {}: {}", results[0].display(), e),
                    }
                })?;
                let result2 = std::fs::read_to_string(&results[1]).map_err(|e| {
                    ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to read {}: {}", results[1].display(), e),
                    }
                })?;

                let data1: serde_json::Value =
                    serde_json::from_str(&result1).map_err(|e| ferro_hgvs::FerroError::Json {
                        msg: format!("Failed to parse {}: {}", results[0].display(), e),
                    })?;
                let data2: serde_json::Value =
                    serde_json::from_str(&result2).map_err(|e| ferro_hgvs::FerroError::Json {
                        msg: format!("Failed to parse {}: {}", results[1].display(), e),
                    })?;

                let tool1 = data1["tool"].as_str().unwrap_or("tool1");
                let tool2 = data2["tool"].as_str().unwrap_or("tool2");

                // Build comparison from results (support both "results" and "sample_results" keys)
                let results1 = data1["results"]
                    .as_array()
                    .or_else(|| data1["sample_results"].as_array());
                let results2 = data2["results"]
                    .as_array()
                    .or_else(|| data2["sample_results"].as_array());

                let mut agreements = 0;
                let mut disagreements = Vec::new();
                let mut both_success = 0;
                let mut tool1_only = 0;
                let mut tool2_only = 0;
                let mut both_fail = 0;

                if let (Some(r1), Some(r2)) = (results1, results2) {
                    for (res1, res2) in r1.iter().zip(r2.iter()) {
                        let s1 = res1["success"].as_bool().unwrap_or(false);
                        let s2 = res2["success"].as_bool().unwrap_or(false);

                        match (s1, s2) {
                            (true, true) => {
                                both_success += 1;
                                let out1 = res1["output"].as_str().unwrap_or("");
                                let out2 = res2["output"].as_str().unwrap_or("");
                                if out1 == out2 {
                                    agreements += 1;
                                } else {
                                    disagreements.push(serde_json::json!({
                                        "input": res1["input"],
                                        tool1: out1,
                                        tool2: out2,
                                    }));
                                }
                            }
                            (true, false) => tool1_only += 1,
                            (false, true) => tool2_only += 1,
                            (false, false) => both_fail += 1,
                        }
                    }
                }

                let total = both_success + tool1_only + tool2_only + both_fail;
                let comparison = serde_json::json!({
                    "mode": "normalize",
                    "files": [results[0].display().to_string(), results[1].display().to_string()],
                    "tools": [tool1, tool2],
                    "total_patterns": total,
                    "both_success": both_success,
                    format!("{}_only_success", tool1): tool1_only,
                    format!("{}_only_success", tool2): tool2_only,
                    "both_fail": both_fail,
                    "agreements": agreements,
                    "disagreements_count": disagreements.len(),
                    "agreement_rate": if both_success > 0 { agreements as f64 / both_success as f64 } else { 0.0 },
                    "disagreements": disagreements,
                });

                std::fs::write(output, serde_json::to_string_pretty(&comparison).unwrap())
                    .map_err(|e| ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to write comparison: {}", e),
                    })?;

                println!("\n=== Comparison Summary ===");
                println!("Tools: {} vs {}", tool1, tool2);
                println!("Total patterns: {}", total);
                println!("Both succeeded: {}", both_success);
                println!(
                    "Agreements: {} ({:.1}%)",
                    agreements,
                    if both_success > 0 {
                        100.0 * agreements as f64 / both_success as f64
                    } else {
                        0.0
                    }
                );
                println!("Disagreements: {}", disagreements.len());
                println!("\nComparison written to {}", output.display());

                // Print overlap tables for 2 tools
                let tool_names = vec![tool1.to_string(), tool2.to_string()];
                let mut all_patterns_set: std::collections::HashSet<String> =
                    std::collections::HashSet::new();
                let mut success_sets: Vec<std::collections::HashSet<String>> = vec![
                    std::collections::HashSet::new(),
                    std::collections::HashSet::new(),
                ];

                if let (Some(r1), Some(r2)) = (results1, results2) {
                    for res in r1.iter() {
                        let input = res["input"].as_str().unwrap_or("").to_string();
                        all_patterns_set.insert(input.clone());
                        if res["success"].as_bool().unwrap_or(false) {
                            success_sets[0].insert(input);
                        }
                    }
                    for res in r2.iter() {
                        let input = res["input"].as_str().unwrap_or("").to_string();
                        all_patterns_set.insert(input.clone());
                        if res["success"].as_bool().unwrap_or(false) {
                            success_sets[1].insert(input);
                        }
                    }
                }

                print_multiway_overlap_tables(&tool_names, &success_sets, &all_patterns_set);
            } else {
                // Multi-way comparison (3+ tools)
                // Load all results
                let mut all_results: Vec<(String, serde_json::Value)> = Vec::new();
                for result_file in results {
                    let content = std::fs::read_to_string(result_file).map_err(|e| {
                        ferro_hgvs::FerroError::Io {
                            msg: format!("Failed to read {}: {}", result_file.display(), e),
                        }
                    })?;
                    let data: serde_json::Value = serde_json::from_str(&content).map_err(|e| {
                        ferro_hgvs::FerroError::Json {
                            msg: format!("Failed to parse {}: {}", result_file.display(), e),
                        }
                    })?;
                    let tool_name = data["tool"].as_str().unwrap_or("unknown").to_string();
                    all_results.push((tool_name, data));
                }

                let tool_names: Vec<String> =
                    all_results.iter().map(|(name, _)| name.clone()).collect();

                // Build success sets
                let mut success_sets: Vec<std::collections::HashSet<String>> = Vec::new();
                let mut all_patterns: std::collections::HashSet<String> =
                    std::collections::HashSet::new();

                for (_, data) in &all_results {
                    let results_arr = data["results"]
                        .as_array()
                        .or_else(|| data["sample_results"].as_array());

                    let mut success_set = std::collections::HashSet::new();
                    if let Some(arr) = results_arr {
                        for item in arr.iter() {
                            let input = item["input"].as_str().unwrap_or("").to_string();
                            all_patterns.insert(input.clone());
                            if item["success"].as_bool().unwrap_or(false) {
                                success_set.insert(input);
                            }
                        }
                    }
                    success_sets.push(success_set);
                }

                // Print overlap tables
                print_multiway_overlap_tables(&tool_names, &success_sets, &all_patterns);

                // Write comparison JSON
                let comparison = serde_json::json!({
                    "mode": "normalize",
                    "files": results.iter().map(|p| p.display().to_string()).collect::<Vec<_>>(),
                    "tools": tool_names,
                    "total_patterns": all_patterns.len(),
                });

                std::fs::write(output, serde_json::to_string_pretty(&comparison).unwrap())
                    .map_err(|e| ferro_hgvs::FerroError::Io {
                        msg: format!("Failed to write comparison: {}", e),
                    })?;

                println!("\nComparison written to {}", output.display());
            }

            Ok(())
        }
    }
}

/// Check if a tool and its dependencies are properly configured
fn run_check_tool(
    tool: Tool,
    reference: &Path,
    mutalyzer_settings: Option<&PathBuf>,
    uta_db_url: Option<&str>,
    seqrepo_path: Option<&PathBuf>,
) -> Result<(), ferro_hgvs::FerroError> {
    match tool {
        Tool::Ferro => {
            // Route to the main library's check module
            let result = check_reference(reference);

            if result.valid {
                print_check_summary(&result, reference);
                println!("\n✓ Ferro is ready!");
                Ok(())
            } else {
                print_check_summary(&result, reference);
                println!(
                    "\nRun: ferro prepare --output-dir {}\nOr:  ferro-benchmark prepare ferro --output-dir {}",
                    reference.display(),
                    reference.display()
                );
                std::process::exit(1);
            }
        }

        Tool::Mutalyzer => {
            println!("Checking mutalyzer configuration...\n");
            let mut all_ok = true;

            // Check parser
            print!("Parser (mutalyzer-hgvs-parser): ");
            if has_mutalyzer_parser() {
                println!("✓ available");
            } else {
                println!("✗ NOT available");
                println!("  Install with: pip install mutalyzer-hgvs-parser");
                all_ok = false;
            }

            // Check normalizer
            print!("Normalizer (mutalyzer-algebra): ");
            if has_mutalyzer_normalizer() {
                println!("✓ available");
            } else {
                println!("✗ NOT available");
                println!("  Install with: pip install mutalyzer");
                all_ok = false;
            }

            // Check cache directory if settings provided
            if let Some(settings_path) = mutalyzer_settings {
                print!("Cache settings: ");
                if settings_path.exists() {
                    println!("✓ found at {}", settings_path.display());

                    // Check if the cache directory exists
                    let settings_content =
                        std::fs::read_to_string(settings_path).unwrap_or_default();
                    if let Some(cache_line) = settings_content
                        .lines()
                        .find(|l| l.starts_with("MUTALYZER_CACHE_DIR"))
                    {
                        let cache_dir = cache_line.split('=').nth(1).map(|s| s.trim());
                        if let Some(dir) = cache_dir {
                            let dir_path = PathBuf::from(dir);
                            if dir_path.exists() {
                                // Count files in cache
                                if let Ok(entries) = std::fs::read_dir(&dir_path) {
                                    let count = entries.count();
                                    println!("  Cache directory: {} ({} files)", dir, count);
                                }
                            } else {
                                println!("  Cache directory: {} (NOT FOUND)", dir);
                                all_ok = false;
                            }
                        }
                    }
                } else {
                    println!("✗ NOT found at {}", settings_path.display());
                    println!("  Run: ferro-benchmark prepare mutalyzer --output-dir <cache-dir> --ferro-reference <reference-dir>");
                    all_ok = false;
                }
            } else {
                println!("Cache settings: not specified (use --mutalyzer-settings)");
            }

            if all_ok {
                println!("\n✓ Mutalyzer is ready!");
                Ok(())
            } else {
                println!("\n✗ Mutalyzer setup is incomplete.");
                std::process::exit(1);
            }
        }

        Tool::Biocommons => {
            println!("Checking biocommons configuration...\n");
            let mut all_ok = true;

            // Check package
            print!("Package (hgvs): ");
            if has_biocommons_normalizer() {
                println!("✓ available");
            } else {
                println!("✗ NOT available");
                println!("  Install with: pip install hgvs");
                all_ok = false;
            }

            // Check Docker
            print!("Docker: ");
            if check_docker_available() {
                println!("✓ available");
            } else {
                println!("✗ NOT available");
                all_ok = false;
            }

            // Check UTA connection
            print!("UTA database: ");
            let uta_url = uta_db_url.unwrap_or("remote uta.biocommons.org");
            if check_uta_connection(uta_db_url) {
                println!("✓ connected ({})", uta_url);
            } else {
                println!("✗ NOT connected ({})", uta_url);
                if uta_db_url.is_some() {
                    println!("  Check that UTA is running: ferro-benchmark start-uta");
                } else {
                    println!("  For local UTA: ferro-benchmark setup-uta");
                }
                all_ok = false;
            }

            // Check SeqRepo
            print!("SeqRepo: ");
            if let Some(path) = seqrepo_path {
                if check_seqrepo(path) {
                    println!("✓ available ({})", path.display());
                } else {
                    println!("✗ NOT available ({})", path.display());
                    all_ok = false;
                }
            } else {
                println!("not specified (use --seqrepo-path)");
                println!("  For local SeqRepo: ferro-benchmark setup-seqrepo --seqrepo-dir <path>");
                all_ok = false;
            }

            if all_ok {
                println!("\n✓ Biocommons is ready!");
                Ok(())
            } else {
                println!("\n✗ Biocommons setup is incomplete.");
                println!(
                    "\nFor full setup: ferro-benchmark prepare biocommons --seqrepo-dir <path>"
                );
                std::process::exit(1);
            }
        }

        Tool::HgvsRs => {
            #[cfg(not(feature = "hgvs-rs"))]
            {
                // Silence unused variable warnings
                let _ = (uta_db_url, seqrepo_path);
                println!("✗ hgvs-rs feature not enabled.\n");
                println!("Rebuild with: cargo build --release --features benchmark,hgvs-rs");
                std::process::exit(1);
            }

            #[cfg(feature = "hgvs-rs")]
            {
                println!("Checking hgvs-rs configuration...\n");
                let mut all_ok = true;

                // Check SeqRepo path
                print!("SeqRepo: ");
                if let Some(path) = seqrepo_path {
                    if check_seqrepo_path(path) {
                        println!("✓ available ({})", path.display());
                    } else {
                        println!("✗ NOT available ({})", path.display());
                        all_ok = false;
                    }
                } else {
                    println!("not specified (use --seqrepo-path)");
                    all_ok = false;
                }

                // Check UTA connection
                print!("UTA database: ");
                // The Python biocommons library requires schema appended to URL
                let uta_url = uta_db_url
                    .unwrap_or("postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b");
                if check_uta_connection(Some(uta_url)) {
                    println!("✓ connected ({})", uta_url);
                } else {
                    println!("✗ NOT connected ({})", uta_url);
                    all_ok = false;
                }

                // Try to create provider
                if all_ok {
                    print!("Provider initialization: ");
                    // hgvs-rs expects URL without schema appended (unlike biocommons)
                    let hgvs_rs_url = "postgresql://anonymous:anonymous@localhost:5432/uta";
                    let config = HgvsRsConfig {
                        uta_db_url: hgvs_rs_url.to_string(),
                        uta_db_schema: "uta_20210129b".to_string(),
                        seqrepo_path: seqrepo_path
                            .map(|p| p.display().to_string())
                            .unwrap_or_default(),
                        lrg_mapping_file: None,
                    };
                    match check_hgvs_rs_available(&config) {
                        Ok(()) => {
                            println!("✓ successful");
                        }
                        Err(e) => {
                            println!("✗ FAILED");
                            println!("  Error: {}", e);
                            all_ok = false;
                        }
                    }
                }

                if all_ok {
                    println!("\n✓ hgvs-rs is ready!");
                    Ok(())
                } else {
                    println!("\n✗ hgvs-rs setup is incomplete.");
                    println!(
                        "\nFor full setup: ferro-benchmark prepare hgvs-rs --seqrepo-dir <path>"
                    );
                    std::process::exit(1);
                }
            }
        }

        Tool::All => {
            // Check all tools and summarize results
            println!("╔════════════════════════════════════════════════════════════════╗");
            println!("║  Checking All Tools                                            ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");

            let mut all_passed = true;
            let mut results: Vec<(&str, bool)> = Vec::new();

            // Check ferro
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            println!("Checking ferro...");
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
            let ferro_ref = reference.join("ferro");
            let ferro_ok = check_reference(&ferro_ref).valid;
            if ferro_ok {
                println!("✓ Ferro is ready ({})", ferro_ref.display());
            } else {
                println!("✗ Ferro not ready at {}", ferro_ref.display());
                all_passed = false;
            }
            results.push(("ferro", ferro_ok));

            // Check mutalyzer
            println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            println!("Checking mutalyzer...");
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
            let mutalyzer_ok = if let Some(settings) = mutalyzer_settings {
                has_mutalyzer_parser() && has_mutalyzer_normalizer() && settings.exists()
            } else {
                let default_settings = reference.join("mutalyzer/mutalyzer_settings.conf");
                has_mutalyzer_parser() && has_mutalyzer_normalizer() && default_settings.exists()
            };
            if mutalyzer_ok {
                println!("✓ Mutalyzer is ready");
            } else {
                println!("✗ Mutalyzer not ready");
                if !has_mutalyzer_parser() {
                    println!("  - Parser not available (pip install mutalyzer-hgvs-parser)");
                }
                if !has_mutalyzer_normalizer() {
                    println!("  - Normalizer not available (pip install mutalyzer)");
                }
                all_passed = false;
            }
            results.push(("mutalyzer", mutalyzer_ok));

            // Check biocommons
            println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            println!("Checking biocommons...");
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
            let biocommons_ok = if has_biocommons_normalizer() && check_docker_available() {
                let seqrepo_check = seqrepo_path.is_some_and(|p| check_seqrepo(p));
                let uta_check = check_uta_connection(uta_db_url);
                seqrepo_check && uta_check
            } else {
                false
            };
            if biocommons_ok {
                println!("✓ Biocommons is ready");
            } else {
                println!("✗ Biocommons not ready");
                if !has_biocommons_normalizer() {
                    println!("  - Package not available (pip install hgvs)");
                }
                if !check_docker_available() {
                    println!("  - Docker not available");
                }
                if seqrepo_path.is_none() {
                    println!("  - SeqRepo path not specified");
                }
                all_passed = false;
            }
            results.push(("biocommons", biocommons_ok));

            // Check hgvs-rs
            println!("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
            println!("Checking hgvs-rs...");
            println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
            #[cfg(not(feature = "hgvs-rs"))]
            let hgvs_rs_ok = {
                println!("✗ hgvs-rs feature not enabled");
                println!("  Rebuild with: cargo build --release --features benchmark,hgvs-rs");
                false
            };
            #[cfg(feature = "hgvs-rs")]
            let hgvs_rs_ok = {
                if let Some(sr_path) = seqrepo_path {
                    if check_seqrepo_path(sr_path) && check_uta_connection(uta_db_url) {
                        let config = HgvsRsConfig {
                            uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta"
                                .to_string(),
                            uta_db_schema: "uta_20210129b".to_string(),
                            seqrepo_path: sr_path.display().to_string(),
                            lrg_mapping_file: None,
                        };
                        check_hgvs_rs_available(&config).is_ok()
                    } else {
                        false
                    }
                } else {
                    false
                }
            };
            if hgvs_rs_ok {
                println!("✓ hgvs-rs is ready");
            } else {
                println!("✗ hgvs-rs not ready");
                all_passed = false;
            }
            results.push(("hgvs-rs", hgvs_rs_ok));

            // Summary
            println!("\n╔════════════════════════════════════════════════════════════════╗");
            println!("║  Summary                                                       ║");
            println!("╚════════════════════════════════════════════════════════════════╝\n");
            for (tool_name, ok) in &results {
                let status = if *ok { "✓ Ready" } else { "✗ Not ready" };
                println!("  {:<12} {}", tool_name, status);
            }

            if all_passed {
                println!("\n✓ All tools are ready!");
                Ok(())
            } else {
                let ready_count = results.iter().filter(|(_, ok)| *ok).count();
                println!("\n✗ {} of 4 tools ready", ready_count);
                println!("\nRun 'ferro-benchmark prepare all' to set up all tools.");
                std::process::exit(1);
            }
        }
    }
}
