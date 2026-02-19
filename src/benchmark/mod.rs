//! Comparison benchmark infrastructure for comparing ferro-hgvs with Mutalyzer.
//!
//! This module provides tools for:
//! - Extracting HGVS patterns from ClinVar and other sources
//! - Sharding datasets for parallel processing
//! - Benchmarking parsing and normalization
//! - Comparing results between ferro-hgvs and Mutalyzer
//! - Generating reports
//!
//! # Usage
//!
//! Build with the benchmark feature:
//! ```bash
//! cargo build --release --features benchmark
//! ```
//!
//! ## Quick: Run 1M Pattern Parsing Benchmark
//!
//! Extract patterns from ClinVar, sample 1M, and benchmark parsing:
//! ```bash
//! # Extract all HGVS patterns from ClinVar
//! ferro-benchmark extract-clinvar -i data/clinvar/hgvs4variation.txt.gz -o /tmp/clinvar_patterns.txt
//!
//! # Create a stratified sample of 1M patterns
//! ferro-benchmark sample -i /tmp/clinvar_patterns.txt -o /tmp/clinvar_1M.txt -s 1000000 --seed 42
//!
//! # Run the parsing benchmark
//! ferro-benchmark parse-ferro -i /tmp/clinvar_1M.txt -r /tmp/results.json -t /tmp/timing.json
//!
//! # View timing results
//! cat /tmp/timing.json
//! ```
//!
//! ## Full Benchmark (with Mutalyzer comparison)
//!
//! Run the complete comparison workflow:
//! ```bash
//! ferro-benchmark run --cores 12 --output-dir benchmark_results
//! ```

pub mod accessions;
pub mod arbitration;
pub mod biocommons;
pub mod cache;
pub mod collate;
pub mod compare;
pub mod extract;
#[cfg(feature = "hgvs-rs")]
pub mod hgvs_rs;
pub mod mutalyzer;
pub mod normalize;
pub mod parse;
pub mod report;
pub mod runner;
pub mod sample;
pub mod shard;
pub mod translate;
pub mod types;
pub mod uta_loader;
pub mod workflow;

pub use accessions::AccessionSources;
pub use arbitration::{
    apply_arbitrations, ArbitrationBuilder, ArbitrationCategory, ArbitrationDatabase,
    ArbitrationDecision, ArbitrationStats, ArbitrationSummary,
};
pub use biocommons::{
    check_container_exists, check_container_running, check_docker_available, check_seqrepo,
    check_uta_connection, check_uta_local, fetch_and_load_missing_accessions,
    has_biocommons_normalizer, load_biocommons_settings, load_sequences_to_seqrepo,
    normalize_single as normalize_biocommons_single, run_biocommons_normalizer_parallel,
    run_biocommons_normalizer_subprocess, setup_seqrepo, setup_uta, start_uta, stop_uta,
    write_biocommons_settings, BiocommonsLocalConfig, BiocommonsSettings, SeqRepoSetupResult,
    UtaSetupResult,
};
pub use cache::{
    build_transcript_chromosome_mapping, download_gene2refseq, enhance_lrg_annotations,
    enhance_nc_annotations, extract_accessions_from_pattern, extract_all_accessions_from_file,
    extract_cdot_mappings, extract_clinvar_accessions, extract_gene2refseq_mappings,
    fetch_fasta_to_file, fetch_lrg_sequences, fetch_missing_accessions, load_transcript_mapping,
    populate_mutalyzer_cache, CacheStats, MappingStats, SupplementalMetadata,
    SupplementalTranscript, TranscriptMapping,
};
pub use collate::{collate_normalization, collate_parsing};
pub use compare::{
    compare_normalize, compare_parse, run_mutalyzer_normalize_parallel, CompareConfig,
};
pub use extract::{extract_clinvar, extract_json};
pub use mutalyzer::MutalyzerClient;
pub use normalize::{normalize_ferro, normalize_mutalyzer};
pub use parse::{parse_ferro, parse_ferro_unified};
pub use runner::{load_existing_results, run_benchmark};
// Re-export from the main prepare module (avoiding code duplication)
pub use crate::prepare::{
    check_references, prepare_references, LegacyMetadata, LegacyTranscript, PrepareConfig,
    ReferenceManifest,
};
pub use report::generate_report;
pub use sample::stratified_sample;
pub use shard::shard_dataset;
pub use types::*;
pub use uta_loader::{load_cdot_alignments_to_uta, UtaLoadConfig, UtaLoadResult};
pub use workflow::run_comparison;

#[cfg(feature = "hgvs-rs")]
pub use hgvs_rs::{
    check_hgvs_rs_available, check_seqrepo_path, run_hgvs_rs_normalize,
    run_hgvs_rs_normalize_parallel, HgvsRsConfig, HgvsRsNormalizer, HgvsRsResult, LrgMapping,
};
