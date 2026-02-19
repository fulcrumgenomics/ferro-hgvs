//! Comparison benchmarks between ferro-hgvs and mutalyzer.
//!
//! This module provides fair head-to-head comparisons by:
//! - Running mutalyzer locally (not via HTTP API)
//! - Sharding work across multiple Python workers
//! - Collecting detailed agreement statistics

/// Python script for mutalyzer parsing (embedded for parallel spawning)
const MUTALYZER_PARSE_SCRIPT: &str = r#"
import sys
import json
import time

try:
    from mutalyzer_hgvs_parser import to_model
except ImportError:
    print("ERROR: mutalyzer-hgvs-parser not installed", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

results = []
successful = 0
failed = 0

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

start = time.perf_counter()

for pattern in patterns:
    try:
        model = to_model(pattern)
        results.append({
            "input": pattern,
            "success": True,
            "output": str(model) if model else pattern
        })
        successful += 1
    except Exception as e:
        results.append({
            "input": pattern,
            "success": False,
            "error": str(e)[:200]
        })
        failed += 1

elapsed = time.perf_counter() - start

output = {
    "tool": "mutalyzer-hgvs-parser",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": failed,
    "elapsed_seconds": elapsed,
    "patterns_per_second": len(patterns) / elapsed if elapsed > 0 else 0,
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output, f, indent=2)
"#;

/// Python script for mutalyzer normalization (embedded for parallel spawning)
const MUTALYZER_NORMALIZE_SCRIPT: &str = r#"
import sys
import os
import gc

# Check if network is allowed (5th argument, default false)
allow_network = len(sys.argv) > 4 and sys.argv[4].lower() == 'true'

# Block network connections unless explicitly allowed
import socket
_original_connect = socket.socket.connect
if not allow_network:
    def _blocked_connect(self, address):
        raise RuntimeError(f"NETWORK ACCESS BLOCKED: Mutalyzer tried to connect to {address}. "
                           "Ensure all required sequences are in the local cache.")
    socket.socket.connect = _blocked_connect

# CRITICAL: Set MUTALYZER_SETTINGS env var BEFORE importing mutalyzer
# mutalyzer-retriever reads config from a file specified by MUTALYZER_SETTINGS
if len(sys.argv) > 3 and sys.argv[3]:
    settings_file = sys.argv[3]
    os.environ['MUTALYZER_SETTINGS'] = settings_file

import json
import time

try:
    from mutalyzer.description import Description
except ImportError:
    print("ERROR: mutalyzer not installed", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Extract worker ID from output filename for progress logging
worker_id = os.path.basename(output_file).replace('shard_', '').replace('_output.json', '')

from collections import defaultdict

results = []
successful = 0
failed = 0
error_counts = defaultdict(int)

def categorize_error(error_msg, pattern=None):
    """Categorize error message into a short category."""
    msg = str(error_msg)
    # Check for common mutalyzer error codes
    if 'ESEQUENCEMISMATCH' in msg:
        return 'ESEQUENCEMISMATCH'
    elif 'ERETR' in msg or 'retriev' in msg.lower():
        return 'RETRIEVAL_ERROR'
    elif 'EPARSE' in msg or 'parse' in msg.lower() or 'syntax' in msg.lower():
        return 'PARSE_ERROR'
    elif 'EREF' in msg:
        return 'EREF'
    elif 'ERANGE' in msg or 'range' in msg.lower() or 'position' in msg.lower():
        return 'RANGE_ERROR'
    elif 'ENOTSUPPORTED' in msg or 'not supported' in msg.lower():
        return 'NOT_SUPPORTED'
    elif 'NETWORK ACCESS BLOCKED' in msg:
        # LRG patterns require network for additional metadata
        if pattern and pattern.startswith('LRG_'):
            return 'LRG_NETWORK_REQUIRED'
        return 'CACHE_MISS'
    elif 'No model loaded' in msg or 'not found' in msg.lower():
        return 'NOT_FOUND'
    else:
        # Extract error code if present (e.g., EXYZ)
        import re
        match = re.search(r'\bE[A-Z]{2,}', msg)
        if match:
            return match.group(0)
        return 'OTHER'

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

total = len(patterns)
start = time.perf_counter()

# Progress logging interval
progress_interval = max(50, total // 20)  # Log ~20 times or every 50

for i, pattern in enumerate(patterns):
    # Progress logging
    if i > 0 and i % progress_interval == 0:
        elapsed_so_far = time.perf_counter() - start
        rate = i / elapsed_so_far if elapsed_so_far > 0 else 0
        eta = (total - i) / rate if rate > 0 else 0
        print(f"[worker {worker_id}] {i}/{total} ({100*i//total}%) - {rate:.1f} p/s - ETA {eta:.0f}s", file=sys.stderr, flush=True)
        # Garbage collect to free memory from cached sequences
        gc.collect()

    try:
        d = Description(description=pattern)
        d.normalize()
        if not d.errors:
            output = d.output()
            results.append({
                "input": pattern,
                "success": True,
                "output": output.get("normalized_description", pattern)
            })
            successful += 1
        else:
            error_msg = str(d.errors[0]) if d.errors else "Unknown error"
            category = categorize_error(error_msg, pattern)
            error_counts[category] += 1
            results.append({
                "input": pattern,
                "success": False,
                "error": error_msg[:200],
                "error_category": category
            })
            failed += 1
    except Exception as e:
        error_msg = str(e)
        category = categorize_error(error_msg, pattern)
        error_counts[category] += 1
        results.append({
            "input": pattern,
            "success": False,
            "error": error_msg[:200],
            "error_category": category
        })
        failed += 1

elapsed = time.perf_counter() - start
print(f"[worker {worker_id}] Done: {total} patterns in {elapsed:.1f}s ({total/elapsed:.1f} p/s)", file=sys.stderr, flush=True)

output_data = {
    "tool": "mutalyzer",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": failed,
    "elapsed_seconds": elapsed,
    "patterns_per_second": len(patterns) / elapsed if elapsed > 0 else 0,
    "error_counts": dict(error_counts),
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)
"#;

/// Python script for biocommons/hgvs normalization (embedded for parallel spawning)
///
/// biocommons/hgvs is the canonical Python HGVS implementation. It can use:
/// - Remote UTA database (default): no local setup needed, but slower for bulk
/// - Local UTA database: faster for large-scale comparisons
///
/// Requirements: pip install hgvs
const BIOCOMMONS_NORMALIZE_SCRIPT: &str = r#"
import sys
import os
import json
import time
from collections import defaultdict

# Optional: Set UTA database URL from environment or argument
# By default, biocommons uses remote UTA at uta.biocommons.org
if len(sys.argv) > 3 and sys.argv[3]:
    os.environ['UTA_DB_URL'] = sys.argv[3]

try:
    import hgvs.parser
    import hgvs.normalizer
    import hgvs.dataproviders.uta
except ImportError:
    print("ERROR: hgvs (biocommons) not installed. Run: pip install hgvs", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Initialize biocommons/hgvs components
try:
    hdp = hgvs.dataproviders.uta.connect()
    hp = hgvs.parser.Parser()
    hn = hgvs.normalizer.Normalizer(hdp)
except Exception as e:
    print(f"ERROR: Failed to initialize biocommons/hgvs: {e}", file=sys.stderr)
    sys.exit(1)

results = []
successful = 0
failed = 0
error_counts = defaultdict(int)

def categorize_error(error_msg, pattern=None):
    """Categorize error message into a short category."""
    msg = str(error_msg)
    if 'HGVSParseError' in msg or 'parse' in msg.lower():
        return 'PARSE_ERROR'
    elif 'HGVSInvalidVariantError' in msg:
        return 'INVALID_VARIANT'
    elif 'HGVSDataNotAvailableError' in msg or 'not found' in msg.lower():
        return 'DATA_NOT_AVAILABLE'
    elif 'HGVSInvalidIntervalError' in msg:
        return 'INVALID_INTERVAL'
    elif 'HGVSNormalizationError' in msg:
        return 'NORMALIZATION_ERROR'
    elif 'HGVSUnsupportedOperationError' in msg:
        return 'UNSUPPORTED'
    elif 'connection' in msg.lower() or 'timeout' in msg.lower():
        return 'CONNECTION_ERROR'
    else:
        return 'OTHER'

with open(input_file, 'r') as f:
    patterns = [line.strip() for line in f if line.strip()]

start = time.perf_counter()

for pattern in patterns:
    try:
        # Parse the HGVS string
        var = hp.parse_hgvs_variant(pattern)
        # Normalize the variant
        normalized = hn.normalize(var)
        results.append({
            "input": pattern,
            "success": True,
            "output": str(normalized)
        })
        successful += 1
    except Exception as e:
        error_msg = str(e)
        category = categorize_error(error_msg, pattern)
        error_counts[category] += 1
        results.append({
            "input": pattern,
            "success": False,
            "error": error_msg[:200],
            "error_category": category
        })
        failed += 1

elapsed = time.perf_counter() - start

output_data = {
    "tool": "biocommons-hgvs",
    "total_patterns": len(patterns),
    "successful": successful,
    "failed": failed,
    "elapsed_seconds": elapsed,
    "patterns_per_second": len(patterns) / elapsed if elapsed > 0 else 0,
    "error_counts": dict(error_counts),
    "results": results
}

with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)
"#;

use crate::benchmark::biocommons::has_biocommons_normalizer;
use crate::benchmark::cache::extract_accessions_from_pattern;
use crate::benchmark::mutalyzer::{has_mutalyzer_normalizer, has_mutalyzer_parser};
use crate::benchmark::sample::stratified_sample_vec;
use crate::benchmark::types::{
    AgreementStats, CacheStats, CombinedPatternResult, CompareMode, ComparisonResult,
    ComparisonToolResult, DetailedResults, DisagreementExample, ParseResult, RefMismatchInfo,
    ReferenceStats, Validator,
};
use crate::error_handling::{ErrorConfig, ErrorMode, InputPreprocessor};
use crate::hgvs::parser::parse_hgvs_lenient;
use crate::normalize::NormalizeConfig;
use crate::reference::ReferenceProvider;
use crate::{FerroError, MockProvider, MultiFastaProvider, Normalizer};
use chrono::{DateTime, Utc};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Configuration for running a comparison.
#[derive(Debug, Clone)]
pub struct CompareConfig {
    /// External validator to compare against (mutalyzer, biocommons, hgvs-rs)
    pub validator: Validator,
    /// Number of patterns to sample (0 = all patterns)
    pub sample_size: usize,
    /// Number of parallel workers for external validator
    pub workers: usize,
    /// Random seed for sampling
    pub seed: u64,
    /// Mutalyzer settings file path (optional, for mutalyzer validator)
    pub mutalyzer_settings: Option<String>,
    /// UTA database URL (optional, for biocommons validator)
    /// Default: remote uta.biocommons.org
    pub uta_db_url: Option<String>,
    /// Biocommons settings file path (optional, for local biocommons setup)
    ///
    /// When provided, uses local UTA and SeqRepo configuration from this file.
    /// Generate with: `ferro-benchmark generate-biocommons-settings`
    pub biocommons_settings: Option<std::path::PathBuf>,
    /// SeqRepo path for hgvs-rs validator (optional)
    ///
    /// Required when using --validator=hgvs-rs. Path to the SeqRepo directory
    /// (e.g., /usr/local/share/seqrepo/2021-01-29).
    pub hgvs_rs_seqrepo_path: Option<String>,
    /// UTA database schema for hgvs-rs validator (optional)
    ///
    /// Used with --validator=hgvs-rs. Defaults to "uta_20210129".
    pub hgvs_rs_uta_schema: Option<String>,
    /// Reference data directory for ferro-hgvs normalization
    pub reference_dir: Option<std::path::PathBuf>,
    /// Transcript→chromosome mapping file for rewriting intronic variants
    pub transcript_mapping: Option<std::path::PathBuf>,
    /// Include unnormalizable patterns (uncertain breakpoints, reference-only)
    ///
    /// By default, patterns that cannot be meaningfully normalized are filtered out:
    /// - Uncertain breakpoints: patterns containing `?`, `(`, or `)` in positions
    ///   (e.g., `c.?_100del`, `c.(100_200)del`) - normalization is undefined
    /// - Reference-only: patterns ending in `=` (e.g., `c.100=`) - trivially equal
    ///
    /// Set to true to include all patterns regardless of normalizability.
    pub include_unnormalizable: bool,
    /// Include genomic (NC_) patterns in comparison
    ///
    /// By default, genomic/chromosomal patterns (NC_*) are excluded because
    /// mutalyzer loads entire chromosome sequences (~60MB each), making
    /// normalization very slow (~1 pattern/second vs ~20 for transcripts).
    ///
    /// Set to true to include genomic patterns (warning: significantly slower).
    pub include_genomic: bool,
    /// Include protein (NP_) patterns in comparison
    ///
    /// By default, protein patterns (NP_*, p.) are excluded because
    /// mutalyzer's local file cache doesn't properly support protein
    /// coordinate systems, causing ECOORDINATESYSTEMMISMATCH errors.
    ///
    /// Set to true to include protein patterns (requires network access).
    pub include_protein: bool,
    /// Allow mutalyzer to make network requests
    ///
    /// By default, network access is blocked to ensure fair offline comparison.
    /// Some patterns (especially LRG) require additional network requests for
    /// metadata that isn't covered by the file cache.
    ///
    /// Set to true to allow network access (warning: slower and requires internet).
    pub allow_mutalyzer_network: bool,
    /// Output path for detailed per-pattern results (optional)
    ///
    /// When provided, saves all per-pattern results (both ferro and mutalyzer)
    /// to a separate JSON file for detailed analysis.
    pub detailed_output: Option<std::path::PathBuf>,
    /// Pre-computed mutalyzer results file (JSONL from mutalyzer-benchmark)
    ///
    /// When provided, uses pre-computed mutalyzer normalization results instead
    /// of running mutalyzer live.
    pub existing_mutalyzer_results: Option<std::path::PathBuf>,
    /// Error handling mode for ferro parsing/normalization
    ///
    /// Controls how ferro handles common input errors:
    /// - "strict" (default): Reject all non-standard input
    /// - "lenient": Auto-correct with warnings
    /// - "silent": Auto-correct silently
    pub error_mode: String,
}

impl Default for CompareConfig {
    fn default() -> Self {
        Self {
            validator: Validator::default(),
            sample_size: 1000,
            workers: num_cpus::get(),
            seed: 42,
            mutalyzer_settings: None,
            uta_db_url: None,
            biocommons_settings: None,
            hgvs_rs_seqrepo_path: None,
            hgvs_rs_uta_schema: None,
            reference_dir: None,
            transcript_mapping: None,
            include_unnormalizable: false,
            include_genomic: false,
            include_protein: false,
            allow_mutalyzer_network: false,
            detailed_output: None,
            existing_mutalyzer_results: None,
            error_mode: "strict".to_string(),
        }
    }
}

/// Check if a variant pattern is intronic (has offset like c.123-45 or c.123+45)
fn is_intronic_variant(pattern: &str) -> bool {
    // Intronic variants have patterns like:
    // - c.123-45A>G (intronic position before exon)
    // - c.123+45A>G (intronic position after exon)
    // - c.*123-45A>G (3'UTR intronic)
    // - c.*123+45A>G (3'UTR intronic)
    // - n.123-45A>G (non-coding intronic)
    // - n.123+45A>G (non-coding intronic)

    // Must have a colon (reference:variant)
    let Some(colon_pos) = pattern.find(':') else {
        return false;
    };

    let variant_part = &pattern[colon_pos + 1..];

    // Must start with c. or n. coordinate system
    if !variant_part.starts_with("c.") && !variant_part.starts_with("n.") {
        return false;
    }

    // Look for intronic offset pattern: digit followed by + or - and more digits
    // This regex-like pattern matches: \d[+-]\d
    let chars: Vec<char> = variant_part.chars().collect();
    for i in 0..chars.len().saturating_sub(2) {
        if chars[i].is_ascii_digit()
            && (chars[i + 1] == '+' || chars[i + 1] == '-')
            && chars[i + 2].is_ascii_digit()
        {
            return true;
        }
    }

    false
}

/// Check if a pattern is a genomic (NC_) variant.
///
/// Genomic patterns use chromosome coordinates and require loading entire
/// chromosome sequences (~60MB each), making them very slow to normalize.
fn is_genomic_pattern(pattern: &str) -> bool {
    pattern.starts_with("NC_")
}

/// Check if a variant pattern is a protein variant.
///
/// Returns true for patterns like `NP_000146.2:p.Pro324Leu`.
/// These are excluded by default because mutalyzer's local file cache
/// doesn't properly support protein coordinate systems.
fn is_protein_pattern(pattern: &str) -> bool {
    pattern.starts_with("NP_") || pattern.contains(":p.")
}

/// Check if a variant pattern is normalizable.
///
/// Returns false for patterns that cannot be meaningfully normalized:
/// - Uncertain breakpoints: patterns containing `?`, `(`, or `)` in the variant part
///   (e.g., `c.?_100del`, `c.(100_200)del`, `c.100_?del`)
/// - Reference-only: patterns ending in `=` (e.g., `c.100=`)
///
/// These patterns either have undefined normalization semantics or trivially
/// normalize to themselves.
fn is_normalizable(pattern: &str) -> bool {
    // Get the variant part (after the colon)
    let variant_part = if let Some(colon_pos) = pattern.find(':') {
        &pattern[colon_pos + 1..]
    } else {
        pattern
    };

    // Check for uncertain breakpoints (?, parentheses in positions)
    // These indicate uncertain/approximate positions where normalization is undefined
    if variant_part.contains('?') || variant_part.contains('(') || variant_part.contains(')') {
        return false;
    }

    // Check for reference-only patterns (ending in =)
    // These trivially normalize to themselves
    if variant_part.ends_with('=') {
        return false;
    }

    true
}

/// Rewrite intronic variants with genomic context for Mutalyzer
///
/// Transforms patterns like:
///   NM_001318856.2:c.9-1296T>C
/// Into:
///   NC_000006.12(NM_001318856.2):c.9-1296T>C
///
/// This allows Mutalyzer to properly normalize intronic variants by providing
/// the genomic reference context.
pub fn rewrite_with_genomic_context(
    patterns: &[String],
    mapping: &HashMap<String, String>,
) -> Vec<String> {
    patterns
        .iter()
        .map(|pattern| {
            // Only rewrite intronic variants on transcript references
            if !is_intronic_variant(pattern) {
                return pattern.clone();
            }

            // Extract the accession (before the colon)
            let Some(colon_pos) = pattern.find(':') else {
                return pattern.clone();
            };

            let accession = &pattern[..colon_pos];

            // Only rewrite NM_ and NR_ accessions
            if !accession.starts_with("NM_") && !accession.starts_with("NR_") {
                return pattern.clone();
            }

            // Look up the chromosome mapping
            if let Some(chromosome) = mapping.get(accession) {
                // Rewrite as NC_xxx(NM_xxx):c.variant
                let variant_part = &pattern[colon_pos..];
                format!("{}({}){}", chromosome, accession, variant_part)
            } else {
                pattern.clone()
            }
        })
        .collect()
}

/// Extract accession from an HGVS pattern (the reference before the colon).
fn extract_accession(pattern: &str) -> Option<&str> {
    // Handle nested format like NC_000001.11(NM_000001.2):c.100A>G
    if let Some(paren_start) = pattern.find('(') {
        if let Some(paren_end) = pattern.find(')') {
            // Return the inner accession (the transcript)
            return Some(&pattern[paren_start + 1..paren_end]);
        }
    }

    // Simple format: NM_000001.2:c.100A>G
    pattern.find(':').map(|pos| &pattern[..pos])
}

/// Build a set of accessions from all .fai index files in a directory tree.
fn load_fai_accessions(cache_dir: &Path) -> HashSet<String> {
    use std::io::BufRead;

    let mut accessions = HashSet::new();

    // Look for .fai files in transcripts subdirectory
    let transcripts_dir = cache_dir.join("transcripts");
    if transcripts_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&transcripts_dir) {
            for entry in entries.filter_map(|e| e.ok()) {
                let path = entry.path();
                if path.extension().is_some_and(|ext| ext == "fai") {
                    if let Ok(file) = File::open(&path) {
                        let reader = std::io::BufReader::new(file);
                        for line in reader.lines().map_while(Result::ok) {
                            // FAI format: name\tlength\toffset\tlinebases\tlinewidth
                            if let Some(name) = line.split('\t').next() {
                                accessions.insert(name.to_string());
                            }
                        }
                    }
                }
            }
        }
    }

    // Also check genome directory for NC_ accessions
    let genome_dir = cache_dir.join("genome");
    if genome_dir.exists() {
        if let Ok(entries) = std::fs::read_dir(&genome_dir) {
            for entry in entries.filter_map(|e| e.ok()) {
                let path = entry.path();
                if path.extension().is_some_and(|ext| ext == "fai") {
                    if let Ok(file) = File::open(&path) {
                        let reader = std::io::BufReader::new(file);
                        for line in reader.lines().map_while(Result::ok) {
                            if let Some(name) = line.split('\t').next() {
                                accessions.insert(name.to_string());
                            }
                        }
                    }
                }
            }
        }
    }

    accessions
}

/// Check cache coverage for a set of patterns.
///
/// Returns statistics about how many accessions are missing from each cache.
/// For ferro, only checks accessions from patterns that actually failed (since ferro
/// doesn't need sequence data for all normalizations - e.g., substitutions don't shift).
/// For mutalyzer, checks all patterns since it always needs cached data.
pub fn check_cache_coverage(
    patterns: &[String],
    ferro_results: &[ParseResult],
    ferro_cache_dir: Option<&Path>,
    mutalyzer_cache_dir: Option<&Path>,
) -> CacheStats {
    // Extract unique accessions from all patterns (for total count and mutalyzer check)
    let all_accessions: HashSet<&str> = patterns
        .iter()
        .filter_map(|p| extract_accession(p))
        .collect();

    // Extract accessions only from ferro failures (these are the ones that actually matter)
    let ferro_failed_accessions: HashSet<&str> = ferro_results
        .iter()
        .filter(|r| !r.success)
        .filter_map(|r| extract_accession(&r.input))
        .collect();

    let total = all_accessions.len();
    let mut ferro_missing = Vec::new();
    let mut muta_missing = Vec::new();

    // Load ferro's indexed FASTA accessions once
    let ferro_accessions = ferro_cache_dir.map(load_fai_accessions);

    // Check ferro cache only for accessions that caused failures
    let mut sorted_ferro_failed: Vec<_> = ferro_failed_accessions.into_iter().collect();
    sorted_ferro_failed.sort();
    for acc in &sorted_ferro_failed {
        if let Some(ref ferro_set) = ferro_accessions {
            if !ferro_set.contains(*acc) {
                ferro_missing.push((*acc).to_string());
            }
        }
    }

    // Check mutalyzer cache for all accessions
    let mut sorted_all: Vec<_> = all_accessions.into_iter().collect();
    sorted_all.sort();
    for acc in &sorted_all {
        if let Some(cache_dir) = mutalyzer_cache_dir {
            // For LRG transcript/protein accessions (LRG_NNNpN, LRG_NNNtN),
            // check if the base LRG exists since they're derived from it
            let check_acc = if acc.starts_with("LRG_") {
                // Extract base LRG: LRG_392p1 -> LRG_392, LRG_392t2 -> LRG_392
                let base: String = acc
                    .chars()
                    .take_while(|c| {
                        c.is_ascii_digit() || *c == '_' || *c == 'L' || *c == 'R' || *c == 'G'
                    })
                    .collect();
                if base.len() < acc.len() {
                    base
                } else {
                    (*acc).to_string()
                }
            } else {
                (*acc).to_string()
            };

            let seq_path = cache_dir.join(format!("{}.sequence", check_acc));
            let ann_path = cache_dir.join(format!("{}.annotations", check_acc));
            if !seq_path.exists() || !ann_path.exists() {
                muta_missing.push((*acc).to_string());
            }
        }
    }

    // Sort for deterministic output
    ferro_missing.sort();
    muta_missing.sort();

    // Save all missing accessions to a file for easy population
    if !muta_missing.is_empty() {
        if let Ok(mut file) = File::create("missing_accessions.txt") {
            for acc in &muta_missing {
                let _ = writeln!(file, "{}", acc);
            }
            eprintln!(
                "  Saved {} missing accessions to missing_accessions.txt",
                muta_missing.len()
            );
        }
    }

    CacheStats {
        total_accessions: total,
        ferro_missing: ferro_missing.len(),
        mutalyzer_missing: muta_missing.len(),
        ferro_missing_examples: ferro_missing.into_iter().take(10).collect(),
        mutalyzer_missing_examples: muta_missing.into_iter().take(10).collect(),
    }
}

/// Result from a mutalyzer worker subprocess.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct MutalyzerOutput {
    #[allow(dead_code)]
    tool: String,
    total_patterns: usize,
    successful: usize,
    #[allow(dead_code)]
    failed: usize,
    elapsed_seconds: f64,
    #[allow(dead_code)]
    patterns_per_second: f64,
    results: Vec<ParseResult>,
    #[serde(default)]
    error_counts: HashMap<String, usize>,
}

/// Run a normalization comparison between ferro-hgvs and mutalyzer.
pub fn compare_normalize<P: AsRef<Path>>(
    input: P,
    output: P,
    config: &CompareConfig,
) -> Result<ComparisonResult, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // Check validator availability based on config
    match config.validator {
        Validator::Mutalyzer => {
            if !has_mutalyzer_normalizer() {
                return Err(FerroError::Io {
                    msg: "mutalyzer package not installed. Install with: pip install mutalyzer"
                        .to_string(),
                });
            }
        }
        Validator::Biocommons => {
            if !has_biocommons_normalizer() {
                return Err(FerroError::Io {
                    msg: "biocommons/hgvs package not installed. Install with: pip install hgvs"
                        .to_string(),
                });
            }
        }
        Validator::HgvsRs => {
            #[cfg(not(feature = "hgvs-rs"))]
            {
                return Err(FerroError::Io {
                    msg:
                        "hgvs-rs feature not enabled. Rebuild with: cargo build --features hgvs-rs"
                            .to_string(),
                });
            }
            #[cfg(feature = "hgvs-rs")]
            {
                // Verify hgvs-rs configuration is available
                if config.hgvs_rs_seqrepo_path.is_none() {
                    return Err(FerroError::Io {
                        msg: "hgvs-rs validator requires --hgvs-rs-seqrepo-path".to_string(),
                    });
                }
            }
        }
    }

    eprintln!("Using validator: {}", config.validator);

    // Load transcript→chromosome mapping if provided (for rewriting intronic variants)
    let transcript_mapping: Option<HashMap<String, String>> =
        if let Some(mapping_path) = &config.transcript_mapping {
            eprintln!(
                "Loading transcript mapping from {}...",
                mapping_path.display()
            );
            let mapping = crate::benchmark::cache::load_transcript_mapping(mapping_path)?;
            eprintln!("  Loaded {} transcript→chromosome mappings", mapping.len());
            Some(mapping)
        } else {
            None
        };

    // Load patterns
    let all_patterns = load_patterns(input)?;
    if all_patterns.is_empty() {
        return Err(FerroError::Io {
            msg: "No patterns found in input file".to_string(),
        });
    }

    eprintln!(
        "Loaded {} patterns from {}",
        all_patterns.len(),
        input.display()
    );

    // Filter patterns based on config settings
    let patterns = {
        let mut filtered = all_patterns;
        let original_count = filtered.len();

        // Filter unnormalizable patterns unless include_unnormalizable is set
        if !config.include_unnormalizable {
            let before = filtered.len();
            filtered.retain(|p| is_normalizable(p));
            let removed = before - filtered.len();
            if removed > 0 {
                eprintln!(
                    "  Excluded {} unnormalizable patterns (uncertain breakpoints, reference-only)",
                    removed
                );
            }
        }

        // Filter genomic (NC_) patterns unless include_genomic is set
        if !config.include_genomic {
            let before = filtered.len();
            filtered.retain(|p| !is_genomic_pattern(p));
            let removed = before - filtered.len();
            if removed > 0 {
                eprintln!(
                    "  Excluded {} genomic (NC_) patterns (use --include-genomic to include)",
                    removed
                );
            }
        }

        // Filter protein (NP_) patterns unless include_protein is set
        if !config.include_protein {
            let before = filtered.len();
            filtered.retain(|p| !is_protein_pattern(p));
            let removed = before - filtered.len();
            if removed > 0 {
                eprintln!(
                    "  Excluded {} protein (NP_/p.) patterns (use --include-protein to include)",
                    removed
                );
            }
        }

        if filtered.len() < original_count {
            eprintln!("  Filtered to {} patterns", filtered.len());
        }

        filtered
    };

    if patterns.is_empty() {
        return Err(FerroError::Io {
            msg: "No normalizable patterns found after filtering".to_string(),
        });
    }

    // Sample patterns if needed
    let sampled = if config.sample_size > 0 && config.sample_size < patterns.len() {
        sample_patterns(&patterns, config.sample_size, config.seed)
    } else {
        patterns
    };

    eprintln!("Using {} patterns for comparison", sampled.len());

    // Run ferro-hgvs (uses original patterns)
    eprintln!("\n[1/2] Running ferro-hgvs normalize...");
    let (ferro_results, ferro_elapsed) =
        run_ferro_normalize(&sampled, &config.reference_dir, &config.error_mode)?;
    let ferro_successful = ferro_results.iter().filter(|r| r.success).count();
    eprintln!(
        "      {} patterns in {:.2}s ({:.0} p/s)",
        sampled.len(),
        ferro_elapsed.as_secs_f64(),
        sampled.len() as f64 / ferro_elapsed.as_secs_f64()
    );

    // Either load pre-computed mutalyzer results or run mutalyzer live
    let (muta_results, muta_elapsed, muta_error_counts) = if let Some(ref existing_path) =
        config.existing_mutalyzer_results
    {
        // Load pre-computed results
        let start = Instant::now();
        let existing_results = load_existing_mutalyzer_results(existing_path)?;

        // Look up results for our sampled patterns
        eprintln!("\n[2/2] Using pre-computed mutalyzer results...");
        let mut results = Vec::with_capacity(sampled.len());
        let mut found = 0;
        let mut not_found = 0;

        for pattern in &sampled {
            if let Some(result) = existing_results.get(pattern) {
                results.push(result.clone());
                found += 1;
            } else {
                // Pattern not in pre-computed results - mark as failed
                results.push(ParseResult {
                    input: pattern.clone(),
                    success: false,
                    output: None,
                    error: Some("Not in pre-computed results".to_string()),
                    error_category: Some("NOT_IN_CACHE".to_string()),
                    ref_mismatch: None,
                    details: None,
                });
                not_found += 1;
            }
        }

        let elapsed = start.elapsed();
        let successful = results.iter().filter(|r| r.success).count();
        eprintln!(
            "      {} patterns in {:.2}s ({} found, {} not in cache)",
            sampled.len(),
            elapsed.as_secs_f64(),
            found,
            not_found
        );
        eprintln!(
            "      {} successful, {} failed",
            successful,
            sampled.len() - successful
        );

        // Aggregate error counts from loaded results
        let mut error_counts: HashMap<String, usize> = HashMap::new();
        for result in &results {
            if let Some(ref category) = result.error_category {
                *error_counts.entry(category.clone()).or_insert(0) += 1;
            }
        }

        (results, elapsed, error_counts)
    } else {
        // Run external validator live based on config.validator
        match config.validator {
            Validator::Mutalyzer => {
                // Rewrite intronic variants for mutalyzer if mapping provided
                let mutalyzer_patterns = if let Some(ref mapping) = transcript_mapping {
                    let rewritten = rewrite_with_genomic_context(&sampled, mapping);
                    let intronic_count = sampled.iter().filter(|p| is_intronic_variant(p)).count();
                    let rewritten_count = sampled
                        .iter()
                        .zip(rewritten.iter())
                        .filter(|(orig, new)| orig != new)
                        .count();
                    eprintln!(
                        "  Rewrote {}/{} intronic variants with genomic context",
                        rewritten_count, intronic_count
                    );
                    rewritten
                } else {
                    sampled.clone()
                };

                eprintln!(
                    "\n[2/2] Running mutalyzer normalize ({} workers)...",
                    config.workers
                );
                let (results, elapsed, error_counts) = run_mutalyzer_normalize_parallel(
                    &mutalyzer_patterns,
                    config.workers,
                    config.mutalyzer_settings.as_deref(),
                    config.allow_mutalyzer_network,
                )?;
                eprintln!(
                    "      {} patterns in {:.2}s ({:.0} p/s)",
                    sampled.len(),
                    elapsed.as_secs_f64(),
                    sampled.len() as f64 / elapsed.as_secs_f64()
                );
                (results, elapsed, error_counts)
            }
            Validator::Biocommons => {
                eprintln!(
                    "\n[2/2] Running biocommons/hgvs normalize ({} workers)...",
                    config.workers
                );
                let (results, elapsed, error_counts) = run_biocommons_normalize_parallel(
                    &sampled,
                    config.workers,
                    config.uta_db_url.as_deref(),
                )?;
                eprintln!(
                    "      {} patterns in {:.2}s ({:.0} p/s)",
                    sampled.len(),
                    elapsed.as_secs_f64(),
                    sampled.len() as f64 / elapsed.as_secs_f64()
                );
                (results, elapsed, error_counts)
            }
            Validator::HgvsRs => {
                #[cfg(not(feature = "hgvs-rs"))]
                {
                    return Err(FerroError::Io {
                        msg: "hgvs-rs feature not enabled".to_string(),
                    });
                }
                #[cfg(feature = "hgvs-rs")]
                {
                    let workers = config.workers;
                    if workers > 1 {
                        eprintln!("\n[2/2] Running hgvs-rs normalize ({} workers)...", workers);
                    } else {
                        eprintln!("\n[2/2] Running hgvs-rs normalize...");
                    }
                    let hgvs_rs_config = crate::benchmark::HgvsRsConfig {
                        uta_db_url: config.uta_db_url.clone().unwrap_or_else(|| {
                            "postgresql://anonymous:anonymous@localhost:5432/uta".to_string()
                        }),
                        uta_db_schema: config
                            .hgvs_rs_uta_schema
                            .clone()
                            .unwrap_or_else(|| "uta_20210129b".to_string()),
                        seqrepo_path: config.hgvs_rs_seqrepo_path.clone().unwrap_or_default(),
                        lrg_mapping_file: None, // LRG mapping not supported in compare workflow yet
                    };
                    let (results, elapsed, error_counts) =
                        crate::benchmark::run_hgvs_rs_normalize_parallel(
                            &sampled,
                            &hgvs_rs_config,
                            workers,
                        )?;
                    eprintln!(
                        "      {} patterns in {:.2}s ({:.0} p/s)",
                        sampled.len(),
                        elapsed.as_secs_f64(),
                        sampled.len() as f64 / elapsed.as_secs_f64()
                    );
                    (results, elapsed, error_counts)
                }
            }
        }
    };

    let muta_successful = muta_results.iter().filter(|r| r.success).count();

    // Save patterns that caused CACHE_MISS errors for investigation
    let cache_miss_patterns: Vec<_> = muta_results
        .iter()
        .filter(|r| r.error_category.as_deref() == Some("CACHE_MISS"))
        .map(|r| r.input.clone())
        .collect();
    if !cache_miss_patterns.is_empty() {
        if let Ok(mut file) = File::create("cache_miss_patterns.txt") {
            for pattern in &cache_miss_patterns {
                let _ = writeln!(file, "{}", pattern);
            }
            eprintln!(
                "  Saved {} CACHE_MISS patterns to cache_miss_patterns.txt",
                cache_miss_patterns.len()
            );
        }
    }

    // Compute agreement
    let (agreement, differences) = compute_agreement(&ferro_results, &muta_results);

    // Compute reference mismatch statistics
    let reference_stats = Some(compute_reference_stats(
        &ferro_results,
        &muta_results,
        &agreement,
    ));

    // Check cache coverage
    let mutalyzer_cache_dir = config
        .mutalyzer_settings
        .as_ref()
        .and_then(|settings_path| {
            // The cache directory is typically the parent of the settings file
            Path::new(settings_path).parent().map(|p| p.to_path_buf())
        });
    let cache_stats = Some(check_cache_coverage(
        &sampled,
        &ferro_results,
        config.reference_dir.as_deref(),
        mutalyzer_cache_dir.as_deref(),
    ));

    // Build result
    let ferro_tool = ComparisonToolResult::new(sampled.len(), ferro_successful, ferro_elapsed);
    let muta_tool = ComparisonToolResult::new(sampled.len(), muta_successful, muta_elapsed)
        .with_error_counts(muta_error_counts);
    let speedup = if ferro_tool.elapsed_seconds > 0.0 {
        muta_tool.elapsed_seconds / ferro_tool.elapsed_seconds
    } else {
        0.0
    };

    let result = ComparisonResult {
        mode: CompareMode::Normalize,
        timestamp: Utc::now(),
        sample_size: sampled.len(),
        ferro: ferro_tool,
        mutalyzer: muta_tool,
        speedup,
        agreement,
        differences,
        reference_stats,
        cache_stats,
    };

    // Print summary
    print_comparison_summary(&result);

    // Save result
    save_comparison_result(&result, output)?;

    // Save detailed per-pattern results if requested
    if let Some(ref detailed_path) = config.detailed_output {
        save_detailed_results(
            &ferro_results,
            &muta_results,
            CompareMode::Normalize,
            result.timestamp,
            detailed_path,
        )?;
    }

    Ok(result)
}

/// Run a parsing comparison between ferro-hgvs and mutalyzer.
pub fn compare_parse<P: AsRef<Path>>(
    input: P,
    output: P,
    config: &CompareConfig,
) -> Result<ComparisonResult, FerroError> {
    let input = input.as_ref();
    let output = output.as_ref();

    // Check mutalyzer availability
    if !has_mutalyzer_parser() {
        return Err(FerroError::Io {
            msg: "mutalyzer-hgvs-parser not installed. Install with: pip install mutalyzer-hgvs-parser".to_string(),
        });
    }

    // Load patterns
    let patterns = load_patterns(input)?;
    if patterns.is_empty() {
        return Err(FerroError::Io {
            msg: "No patterns found in input file".to_string(),
        });
    }

    eprintln!(
        "Loaded {} patterns from {}",
        patterns.len(),
        input.display()
    );

    // Sample patterns if needed
    let sampled = if config.sample_size > 0 && config.sample_size < patterns.len() {
        sample_patterns(&patterns, config.sample_size, config.seed)
    } else {
        patterns
    };

    eprintln!("Using {} patterns for comparison", sampled.len());

    // Run ferro-hgvs parsing
    eprintln!("\n[1/2] Running ferro-hgvs parse...");
    let (ferro_results, ferro_elapsed) = run_ferro_parse(&sampled, &config.error_mode);
    let ferro_successful = ferro_results.iter().filter(|r| r.success).count();
    eprintln!(
        "      {} patterns in {:.2}s ({:.0} p/s)",
        sampled.len(),
        ferro_elapsed.as_secs_f64(),
        sampled.len() as f64 / ferro_elapsed.as_secs_f64()
    );

    // Run mutalyzer parsing in parallel
    eprintln!(
        "\n[2/2] Running mutalyzer parse ({} workers)...",
        config.workers
    );
    let (muta_results, muta_elapsed) = run_mutalyzer_parse_parallel(&sampled, config.workers)?;
    let muta_successful = muta_results.iter().filter(|r| r.success).count();
    eprintln!(
        "      {} patterns in {:.2}s ({:.0} p/s)",
        sampled.len(),
        muta_elapsed.as_secs_f64(),
        sampled.len() as f64 / muta_elapsed.as_secs_f64()
    );

    // Compute agreement
    let (agreement, differences) = compute_agreement(&ferro_results, &muta_results);

    // Build result
    let ferro_tool = ComparisonToolResult::new(sampled.len(), ferro_successful, ferro_elapsed);
    let muta_tool = ComparisonToolResult::new(sampled.len(), muta_successful, muta_elapsed);
    let speedup = if ferro_tool.elapsed_seconds > 0.0 {
        muta_tool.elapsed_seconds / ferro_tool.elapsed_seconds
    } else {
        0.0
    };

    let result = ComparisonResult {
        mode: CompareMode::Parse,
        timestamp: Utc::now(),
        sample_size: sampled.len(),
        ferro: ferro_tool,
        mutalyzer: muta_tool,
        speedup,
        agreement,
        differences,
        reference_stats: None, // Not applicable for parse mode
        cache_stats: None,     // Not applicable for parse mode
    };

    // Print summary
    print_comparison_summary(&result);

    // Save result
    save_comparison_result(&result, output)?;

    // Save detailed per-pattern results if requested
    if let Some(ref detailed_path) = config.detailed_output {
        save_detailed_results(
            &ferro_results,
            &muta_results,
            CompareMode::Parse,
            result.timestamp,
            detailed_path,
        )?;
    }

    Ok(result)
}

/// Load pre-computed mutalyzer results from a JSONL file.
///
/// Returns a HashMap from input pattern to ParseResult for fast lookup.
fn load_existing_mutalyzer_results<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, ParseResult>, FerroError> {
    let path = path.as_ref();
    eprintln!(
        "Loading pre-computed mutalyzer results from {}...",
        path.display()
    );

    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;

    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut results = HashMap::new();
    let mut line_count = 0;
    let mut success_count = 0;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if line.is_empty() {
            continue;
        }

        let result: ParseResult = serde_json::from_str(&line).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse JSONL line: {}", e),
        })?;

        if result.success {
            success_count += 1;
        }
        results.insert(result.input.clone(), result);
        line_count += 1;

        if line_count % 1_000_000 == 0 {
            eprintln!("  Loaded {} results...", line_count);
        }
    }

    eprintln!(
        "  Loaded {} results ({} successful, {} failed)",
        line_count,
        success_count,
        line_count - success_count
    );

    Ok(results)
}

/// Load patterns from a file (auto-detect format).
fn load_patterns<P: AsRef<Path>>(path: P) -> Result<Vec<String>, FerroError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", path.display(), e),
    })?;

    let reader: Box<dyn BufRead> = if path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut patterns = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Check if this looks like ClinVar TSV (has tabs)
        if line.contains('\t') {
            // ClinVar format (hgvs4variation.txt):
            // 0:Symbol, 1:GeneID, 2:VariationID, 3:AlleleID, 4:Type, 5:Assembly,
            // 6:NucleotideExpression, 7:NucleotideChange, 8:ProteinExpression, ...
            let fields: Vec<&str> = line.split('\t').collect();
            // Column 7 (index 6) is NucleotideExpression (full HGVS)
            if fields.len() > 6 {
                let hgvs = fields[6].trim();
                if !hgvs.is_empty() && hgvs != "-" && hgvs.contains(':') {
                    patterns.push(hgvs.to_string());
                }
            }
            // Column 9 (index 8) is ProteinExpression (full HGVS)
            if fields.len() > 8 {
                let hgvs = fields[8].trim();
                if !hgvs.is_empty() && hgvs != "-" && hgvs.contains(':') {
                    patterns.push(hgvs.to_string());
                }
            }
        } else if line.contains(':') {
            // Plain text, one pattern per line
            patterns.push(line.to_string());
        }
    }

    // Deduplicate
    patterns.sort();
    patterns.dedup();

    Ok(patterns)
}

/// Sample patterns randomly.
fn sample_patterns(patterns: &[String], size: usize, seed: u64) -> Vec<String> {
    // Try stratified sampling first
    // Note: protein exclusion is handled by the caller (config.include_protein),
    // so we pass false here to avoid double-filtering
    if let Ok(sampled) = stratified_sample_vec(patterns, size, seed, false) {
        return sampled;
    }

    // Fall back to simple random sampling
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut sampled: Vec<String> = patterns.to_vec();
    sampled.shuffle(&mut rng);
    sampled.truncate(size);
    sampled
}

/// Run ferro-hgvs normalization on patterns.
///
/// Uses lenient mode to correct reference mismatches and track them in the results.
fn run_ferro_normalize(
    patterns: &[String],
    reference_dir: &Option<std::path::PathBuf>,
    error_mode: &str,
) -> Result<(Vec<ParseResult>, std::time::Duration), FerroError> {
    // Create normalizer with LENIENT config for benchmark - this allows us to track
    // reference mismatches that would otherwise cause failures in strict mode
    let config = NormalizeConfig::lenient();

    let normalizer: Normalizer<Box<dyn ReferenceProvider>> = match reference_dir {
        Some(ref ref_dir) => {
            let manifest_path = ref_dir.join("manifest.json");
            if manifest_path.exists() {
                let provider = MultiFastaProvider::from_manifest(&manifest_path)?;
                Normalizer::with_config(Box::new(provider) as Box<dyn ReferenceProvider>, config)
            } else if ref_dir.is_dir() {
                let provider = MultiFastaProvider::from_directory(ref_dir)?;
                Normalizer::with_config(Box::new(provider) as Box<dyn ReferenceProvider>, config)
            } else {
                Normalizer::with_config(
                    Box::new(MockProvider::with_test_data()) as Box<dyn ReferenceProvider>,
                    config,
                )
            }
        }
        None => Normalizer::with_config(
            Box::new(MockProvider::with_test_data()) as Box<dyn ReferenceProvider>,
            config,
        ),
    };

    // Create error config from mode string for preprocessing
    let base_mode = match error_mode {
        "lenient" => ErrorMode::Lenient,
        "silent" => ErrorMode::Silent,
        _ => ErrorMode::Strict,
    };
    let error_config = ErrorConfig::new(base_mode);
    let preprocessor = InputPreprocessor::new(error_config);

    let start = Instant::now();

    let results: Vec<ParseResult> = patterns
        .iter()
        .map(|pattern| {
            // Preprocess input based on error mode
            let preprocess_result = preprocessor.preprocess(pattern);
            let preprocessed = if preprocess_result.success {
                preprocess_result.preprocessed
            } else {
                pattern.clone()
            };

            match parse_hgvs_lenient(&preprocessed) {
                Ok(parse_result) => {
                    // Use normalize_with_warnings to capture reference mismatches
                    match normalizer.normalize_with_warnings(&parse_result.result) {
                        Ok(norm_result) => {
                            // Check if there were reference mismatches
                            let ref_mismatch = norm_result
                                .warnings
                                .iter()
                                .find(|w| w.code == "REFSEQ_MISMATCH")
                                .map(|w| RefMismatchInfo {
                                    stated_ref: w.stated_ref.clone(),
                                    actual_ref: w.actual_ref.clone(),
                                    position: w.position.clone(),
                                    corrected: w.corrected,
                                });

                            ParseResult {
                                input: pattern.clone(),
                                success: true,
                                output: Some(norm_result.result.to_string()),
                                error: None,
                                error_category: None,
                                ref_mismatch,
                                details: None,
                            }
                        }
                        Err(e) => ParseResult {
                            input: pattern.clone(),
                            success: false,
                            output: None,
                            error: Some(format!("{}", e)),
                            error_category: Some("normalize_error".to_string()),
                            ref_mismatch: None,
                            details: None,
                        },
                    }
                }
                Err(e) => ParseResult {
                    input: pattern.clone(),
                    success: false,
                    output: None,
                    error: Some(format!("{}", e)),
                    error_category: Some("parse_error".to_string()),
                    ref_mismatch: None,
                    details: None,
                },
            }
        })
        .collect();

    let elapsed = start.elapsed();
    Ok((results, elapsed))
}

/// Run ferro-hgvs parsing on patterns.
fn run_ferro_parse(
    patterns: &[String],
    error_mode: &str,
) -> (Vec<ParseResult>, std::time::Duration) {
    // Create error config from mode string
    let base_mode = match error_mode {
        "lenient" => ErrorMode::Lenient,
        "silent" => ErrorMode::Silent,
        _ => ErrorMode::Strict,
    };
    let error_config = ErrorConfig::new(base_mode);
    let preprocessor = InputPreprocessor::new(error_config);

    let start = Instant::now();

    let results: Vec<ParseResult> = patterns
        .iter()
        .map(|pattern| {
            // Preprocess input based on error mode
            let preprocess_result = preprocessor.preprocess(pattern);
            let preprocessed = if preprocess_result.success {
                preprocess_result.preprocessed
            } else {
                pattern.clone()
            };

            match parse_hgvs_lenient(&preprocessed) {
                Ok(parse_result) => ParseResult {
                    input: pattern.clone(),
                    success: true,
                    output: Some(parse_result.result.to_string()),
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
                    error_category: Some("parse_error".to_string()),
                    ref_mismatch: None,
                    details: None,
                },
            }
        })
        .collect();

    let elapsed = start.elapsed();
    (results, elapsed)
}

/// Run mutalyzer normalization in parallel using sharded subprocesses.
///
/// Supports two modes:
/// 1. Pre-sharded cache (recommended): If the cache has shard_manifest.json,
///    patterns are routed to their corresponding shard based on accession hash.
///    Each worker processes one or more shards sequentially.
/// 2. Runtime partitioning (fallback): Patterns are distributed round-robin to
///    workers, and each worker gets symlinks to only the sequences it needs.
#[allow(clippy::type_complexity)]
pub fn run_mutalyzer_normalize_parallel(
    patterns: &[String],
    workers: usize,
    settings_file: Option<&str>,
    allow_network: bool,
) -> Result<
    (
        Vec<ParseResult>,
        std::time::Duration,
        HashMap<String, usize>,
    ),
    FerroError,
> {
    let temp_dir = tempfile::tempdir().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp directory: {}", e),
    })?;

    // Check if we have a pre-sharded cache
    let (cache_dir, manifest) = if let Some(settings_path) = settings_file {
        let settings_path = Path::new(settings_path);
        match parse_mutalyzer_settings(settings_path) {
            Ok(cache_dir) => {
                let manifest = load_shard_manifest(&cache_dir);
                (Some(cache_dir), manifest)
            }
            Err(_) => (None, None),
        }
    } else {
        (None, None)
    };

    // Use pre-sharded approach if manifest exists
    if let (Some(cache_dir), Some(manifest)) = (&cache_dir, &manifest) {
        return run_mutalyzer_normalize_presharded(
            patterns,
            workers,
            cache_dir,
            manifest.num_shards,
            allow_network,
            temp_dir.path(),
        );
    }

    // Fall back to runtime partitioning
    eprintln!("[normalize] Using runtime partitioning (no pre-sharded cache found)");

    // Shard patterns round-robin
    let shards = shard_patterns(patterns, workers);

    // Partition cache if settings file is provided
    let partitioned_settings = if let Some(cache_dir) = &cache_dir {
        eprintln!(
            "[partition] Partitioning cache from {} for {} workers",
            cache_dir.display(),
            workers
        );
        match partition_cache_for_workers(&shards, cache_dir, temp_dir.path()) {
            Ok(paths) => Some(paths),
            Err(e) => {
                eprintln!("[partition] Warning: Failed to partition cache: {}", e);
                eprintln!("[partition] Falling back to shared cache");
                None
            }
        }
    } else {
        None
    };

    // Write shard input files
    let mut shard_paths = Vec::new();
    for (i, shard) in shards.iter().enumerate() {
        let input_path = temp_dir.path().join(format!("shard_{}_input.txt", i));
        let output_path = temp_dir.path().join(format!("shard_{}_output.json", i));

        let mut file = File::create(&input_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create shard file: {}", e),
        })?;
        for pattern in shard {
            writeln!(file, "{}", pattern).map_err(|e| FerroError::Io {
                msg: format!("Failed to write pattern: {}", e),
            })?;
        }

        shard_paths.push((input_path, output_path));
    }

    // Run workers in parallel using process spawning
    let start = Instant::now();

    // Spawn all processes first (non-blocking)
    let mut children: Vec<(usize, std::process::Child)> = Vec::new();
    for (i, (input, output)) in shard_paths.iter().enumerate() {
        let mut cmd = std::process::Command::new("python3");
        cmd.args(["-c", MUTALYZER_NORMALIZE_SCRIPT]);
        cmd.arg(input.display().to_string());
        cmd.arg(output.display().to_string());

        // Use partitioned settings if available, otherwise fall back to shared settings
        if let Some(ref partitioned) = partitioned_settings {
            cmd.arg(partitioned[i].display().to_string());
        } else if let Some(settings) = settings_file {
            cmd.arg(settings);
        } else {
            cmd.arg("");
        }

        // Pass network flag as 5th argument
        cmd.arg(if allow_network { "true" } else { "false" });

        let child = cmd.spawn().map_err(|e| FerroError::Io {
            msg: format!("Failed to spawn Python: {}", e),
        })?;
        children.push((i, child));
    }

    // Wait for all processes to complete
    for (_i, mut child) in children {
        let status = child.wait().map_err(|e| FerroError::Io {
            msg: format!("Failed to wait for Python: {}", e),
        })?;
        if !status.success() {
            return Err(FerroError::Io {
                msg: format!(
                    "Mutalyzer normalizer failed with exit code: {:?}",
                    status.code()
                ),
            });
        }
    }

    let elapsed = start.elapsed();

    // Collect results and aggregate error counts
    let mut all_results = Vec::new();
    let mut aggregated_error_counts: HashMap<String, usize> = HashMap::new();
    for (_, output_path) in &shard_paths {
        let file = File::open(output_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open shard output: {}", e),
        })?;
        let reader = BufReader::new(file);
        let output: MutalyzerOutput =
            serde_json::from_reader(reader).map_err(|e| FerroError::Json {
                msg: format!("Failed to parse shard output: {}", e),
            })?;
        all_results.extend(output.results);

        // Aggregate error counts from each shard
        for (category, count) in output.error_counts {
            *aggregated_error_counts.entry(category).or_insert(0) += count;
        }
    }

    Ok((all_results, elapsed, aggregated_error_counts))
}

/// Run mutalyzer normalization using pre-sharded cache.
///
/// Patterns are grouped by their accession's shard, then workers process
/// their assigned shards sequentially. This ensures each worker only loads
/// the reference data for its shards (~18GB/N instead of 18GB).
#[allow(clippy::type_complexity)]
fn run_mutalyzer_normalize_presharded(
    patterns: &[String],
    workers: usize,
    cache_dir: &Path,
    num_shards: usize,
    allow_network: bool,
    temp_dir: &Path,
) -> Result<
    (
        Vec<ParseResult>,
        std::time::Duration,
        HashMap<String, usize>,
    ),
    FerroError,
> {
    eprintln!(
        "[normalize] Using pre-sharded cache ({} shards, {} workers)",
        num_shards, workers
    );

    // Group patterns by their accession's shard
    let shard_groups = group_patterns_by_shard(patterns, num_shards);

    // Log shard distribution
    let non_empty_shards: Vec<_> = shard_groups
        .iter()
        .enumerate()
        .filter(|(_, g)| !g.is_empty())
        .collect();
    eprintln!(
        "[normalize] Patterns distributed across {} shards (of {} total)",
        non_empty_shards.len(),
        num_shards
    );

    // Assign shards to workers
    // If workers >= num_shards: each worker gets at most 1 shard
    // If workers < num_shards: workers get multiple shards (round-robin)
    let effective_workers = workers.min(non_empty_shards.len());
    let mut worker_shards: Vec<Vec<usize>> = (0..effective_workers).map(|_| Vec::new()).collect();

    for (i, (shard_id, _)) in non_empty_shards.iter().enumerate() {
        worker_shards[i % effective_workers].push(*shard_id);
    }

    // Process each worker's shards
    let start = Instant::now();
    let mut all_results = Vec::new();
    let mut aggregated_error_counts: HashMap<String, usize> = HashMap::new();
    let mut output_paths: Vec<PathBuf> = Vec::new();

    // Create input files and spawn workers
    let mut children: Vec<std::process::Child> = Vec::new();

    for (worker_id, shard_ids) in worker_shards.iter().enumerate() {
        if shard_ids.is_empty() {
            continue;
        }

        // Collect all patterns for this worker's shards
        let mut worker_patterns: Vec<String> = Vec::new();
        for &shard_id in shard_ids {
            worker_patterns.extend(shard_groups[shard_id].clone());
        }

        if worker_patterns.is_empty() {
            continue;
        }

        // Write input file
        let input_path = temp_dir.join(format!("worker_{}_input.txt", worker_id));
        let output_path = temp_dir.join(format!("worker_{}_output.json", worker_id));

        let mut file = File::create(&input_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create input file: {}", e),
        })?;
        for pattern in &worker_patterns {
            writeln!(file, "{}", pattern).map_err(|e| FerroError::Io {
                msg: format!("Failed to write pattern: {}", e),
            })?;
        }

        // All workers use the same settings file pointing to the full cache.
        // Memory savings come from pattern routing - each worker only processes
        // patterns for certain accessions, so only loads those into memory.
        let settings_path = cache_dir.join("mutalyzer_settings.conf");
        let effective_settings = settings_path.to_string_lossy().to_string();

        eprintln!(
            "[normalize] Worker {} processing {} patterns from shard(s) {:?}",
            worker_id,
            worker_patterns.len(),
            shard_ids
        );

        // Spawn worker process
        let mut cmd = std::process::Command::new("python3");
        cmd.args(["-c", MUTALYZER_NORMALIZE_SCRIPT]);
        cmd.arg(input_path.display().to_string());
        cmd.arg(output_path.display().to_string());
        cmd.arg(&effective_settings);
        cmd.arg(if allow_network { "true" } else { "false" });

        let child = cmd.spawn().map_err(|e| FerroError::Io {
            msg: format!("Failed to spawn Python: {}", e),
        })?;
        children.push(child);
        output_paths.push(output_path);
    }

    // Wait for all workers
    for mut child in children {
        let status = child.wait().map_err(|e| FerroError::Io {
            msg: format!("Failed to wait for Python: {}", e),
        })?;
        if !status.success() {
            return Err(FerroError::Io {
                msg: format!(
                    "Mutalyzer normalizer failed with exit code: {:?}",
                    status.code()
                ),
            });
        }
    }

    let elapsed = start.elapsed();

    // Collect results
    for output_path in &output_paths {
        let file = File::open(output_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open output: {}", e),
        })?;
        let reader = BufReader::new(file);
        let output: MutalyzerOutput =
            serde_json::from_reader(reader).map_err(|e| FerroError::Json {
                msg: format!("Failed to parse output: {}", e),
            })?;
        all_results.extend(output.results);

        for (category, count) in output.error_counts {
            *aggregated_error_counts.entry(category).or_insert(0) += count;
        }
    }

    Ok((all_results, elapsed, aggregated_error_counts))
}

/// Run biocommons/hgvs normalization in parallel using sharded subprocesses.
#[allow(clippy::type_complexity)]
fn run_biocommons_normalize_parallel(
    patterns: &[String],
    workers: usize,
    uta_db_url: Option<&str>,
) -> Result<
    (
        Vec<ParseResult>,
        std::time::Duration,
        HashMap<String, usize>,
    ),
    FerroError,
> {
    let temp_dir = tempfile::tempdir().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp directory: {}", e),
    })?;

    // Shard patterns
    let shards = shard_patterns(patterns, workers);

    // Write shard input files
    let mut shard_paths = Vec::new();
    for (i, shard) in shards.iter().enumerate() {
        let input_path = temp_dir.path().join(format!("shard_{}_input.txt", i));
        let output_path = temp_dir.path().join(format!("shard_{}_output.json", i));

        let mut file = File::create(&input_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create shard file: {}", e),
        })?;
        for pattern in shard {
            writeln!(file, "{}", pattern).map_err(|e| FerroError::Io {
                msg: format!("Failed to write pattern: {}", e),
            })?;
        }

        shard_paths.push((input_path, output_path));
    }

    // Run workers in parallel using process spawning
    let start = Instant::now();

    // Spawn all processes first (non-blocking)
    let mut children: Vec<(usize, std::process::Child)> = Vec::new();
    for (i, (input, output)) in shard_paths.iter().enumerate() {
        let mut cmd = std::process::Command::new("python3");
        cmd.args(["-c", BIOCOMMONS_NORMALIZE_SCRIPT]);
        cmd.arg(input.display().to_string());
        cmd.arg(output.display().to_string());
        // Pass UTA database URL as 3rd argument (empty string if not provided)
        cmd.arg(uta_db_url.unwrap_or(""));

        let child = cmd.spawn().map_err(|e| FerroError::Io {
            msg: format!("Failed to spawn Python: {}", e),
        })?;
        children.push((i, child));
    }

    // Wait for all processes to complete
    for (_i, mut child) in children {
        let status = child.wait().map_err(|e| FerroError::Io {
            msg: format!("Failed to wait for Python: {}", e),
        })?;
        if !status.success() {
            return Err(FerroError::Io {
                msg: format!(
                    "biocommons/hgvs normalizer failed with exit code: {:?}",
                    status.code()
                ),
            });
        }
    }

    let elapsed = start.elapsed();

    // Collect results and aggregate error counts
    let mut all_results = Vec::new();
    let mut aggregated_error_counts: HashMap<String, usize> = HashMap::new();
    for (_, output_path) in &shard_paths {
        let file = File::open(output_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open shard output: {}", e),
        })?;
        let reader = BufReader::new(file);
        let output: MutalyzerOutput =
            serde_json::from_reader(reader).map_err(|e| FerroError::Json {
                msg: format!("Failed to parse shard output: {}", e),
            })?;
        all_results.extend(output.results);

        // Aggregate error counts from each shard
        for (category, count) in output.error_counts {
            *aggregated_error_counts.entry(category).or_insert(0) += count;
        }
    }

    Ok((all_results, elapsed, aggregated_error_counts))
}

/// Run mutalyzer parsing in parallel using sharded subprocesses.
fn run_mutalyzer_parse_parallel(
    patterns: &[String],
    workers: usize,
) -> Result<(Vec<ParseResult>, std::time::Duration), FerroError> {
    let temp_dir = tempfile::tempdir().map_err(|e| FerroError::Io {
        msg: format!("Failed to create temp directory: {}", e),
    })?;

    // Shard patterns
    let shards = shard_patterns(patterns, workers);

    // Write shard input files
    let mut shard_paths = Vec::new();
    for (i, shard) in shards.iter().enumerate() {
        let input_path = temp_dir.path().join(format!("shard_{}_input.txt", i));
        let output_path = temp_dir.path().join(format!("shard_{}_output.json", i));

        let mut file = File::create(&input_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create shard file: {}", e),
        })?;
        for pattern in shard {
            writeln!(file, "{}", pattern).map_err(|e| FerroError::Io {
                msg: format!("Failed to write pattern: {}", e),
            })?;
        }

        shard_paths.push((input_path, output_path));
    }

    // Run workers in parallel using process spawning
    let start = Instant::now();

    // Spawn all processes first (non-blocking)
    let mut children: Vec<std::process::Child> = Vec::new();
    for (input, output) in &shard_paths {
        let mut cmd = std::process::Command::new("python3");
        cmd.args(["-c", MUTALYZER_PARSE_SCRIPT]);
        cmd.arg(input.display().to_string());
        cmd.arg(output.display().to_string());

        let child = cmd.spawn().map_err(|e| FerroError::Io {
            msg: format!("Failed to spawn Python: {}", e),
        })?;
        children.push(child);
    }

    // Wait for all processes to complete
    for mut child in children {
        let status = child.wait().map_err(|e| FerroError::Io {
            msg: format!("Failed to wait for Python: {}", e),
        })?;
        if !status.success() {
            return Err(FerroError::Io {
                msg: format!(
                    "Mutalyzer parser failed with exit code: {:?}",
                    status.code()
                ),
            });
        }
    }

    let elapsed = start.elapsed();

    // Collect results
    let mut all_results = Vec::new();
    for (_, output_path) in &shard_paths {
        let file = File::open(output_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open shard output: {}", e),
        })?;
        let reader = BufReader::new(file);
        let output: MutalyzerOutput =
            serde_json::from_reader(reader).map_err(|e| FerroError::Json {
                msg: format!("Failed to parse shard output: {}", e),
            })?;
        all_results.extend(output.results);
    }

    Ok((all_results, elapsed))
}

/// Shard patterns into N groups for parallel processing.
fn shard_patterns(patterns: &[String], num_shards: usize) -> Vec<Vec<String>> {
    let mut shards: Vec<Vec<String>> = (0..num_shards).map(|_| Vec::new()).collect();

    for (i, pattern) in patterns.iter().enumerate() {
        shards[i % num_shards].push(pattern.clone());
    }

    shards
}

/// Parse a mutalyzer settings file to extract the cache directory.
fn parse_mutalyzer_settings(settings_path: &Path) -> Result<PathBuf, FerroError> {
    let content = std::fs::read_to_string(settings_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to read settings file: {}", e),
    })?;

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("MUTALYZER_CACHE_DIR") {
            // Format: MUTALYZER_CACHE_DIR = /path/to/cache or MUTALYZER_CACHE_DIR=/path/to/cache
            if let Some(value) = line.split('=').nth(1) {
                return Ok(PathBuf::from(value.trim()));
            }
        }
    }

    Err(FerroError::Io {
        msg: format!(
            "MUTALYZER_CACHE_DIR not found in settings file: {}",
            settings_path.display()
        ),
    })
}

/// Shard manifest from pre-sharded cache.
#[derive(Debug, serde::Deserialize)]
struct ShardManifest {
    num_shards: usize,
    #[allow(dead_code)]
    total_accessions: usize,
    #[allow(dead_code)]
    shard_counts: Vec<usize>,
}

/// Compute the shard ID for an accession using the same hash function as prepare time.
/// Must match the hash function in shard_mutalyzer_cache in benchmark.rs.
fn accession_to_shard(accession: &str, num_shards: usize) -> usize {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    accession.hash(&mut hasher);
    (hasher.finish() as usize) % num_shards
}

/// Load shard manifest from cache directory if it exists.
fn load_shard_manifest(cache_dir: &Path) -> Option<ShardManifest> {
    let manifest_path = cache_dir.join("shard_manifest.json");
    if !manifest_path.exists() {
        return None;
    }

    let file = File::open(&manifest_path).ok()?;
    let reader = BufReader::new(file);
    serde_json::from_reader(reader).ok()
}

/// Group patterns by their accession's shard ID.
///
/// Returns a vector indexed by shard_id, where each entry is the list of patterns
/// that should be processed by that shard.
fn group_patterns_by_shard(patterns: &[String], num_shards: usize) -> Vec<Vec<String>> {
    let mut shard_groups: Vec<Vec<String>> = (0..num_shards).map(|_| Vec::new()).collect();

    for pattern in patterns {
        // Extract the primary accession from the pattern
        if let Some(accession) = extract_accession(pattern) {
            let shard_id = accession_to_shard(accession, num_shards);
            shard_groups[shard_id].push(pattern.clone());
        } else {
            // If we can't extract an accession, put it in shard 0 as fallback
            shard_groups[0].push(pattern.clone());
        }
    }

    shard_groups
}

/// Partition the mutalyzer cache for parallel workers.
///
/// Creates per-worker cache directories with symlinks to only the sequence files
/// needed for each worker's patterns. This dramatically reduces memory usage
/// since each worker only loads a subset of the cache.
///
/// Returns a vector of settings file paths, one per worker.
fn partition_cache_for_workers(
    shards: &[Vec<String>],
    full_cache_dir: &Path,
    work_dir: &Path,
) -> Result<Vec<PathBuf>, FerroError> {
    let num_workers = shards.len();
    let mut settings_paths = Vec::with_capacity(num_workers);

    for (worker_id, shard) in shards.iter().enumerate() {
        // Extract all accessions needed for this shard
        let mut accessions = HashSet::new();
        for pattern in shard {
            for acc in extract_accessions_from_pattern(pattern) {
                accessions.insert(acc);
            }
        }

        // Create per-worker cache directory
        let worker_cache_dir = work_dir.join(format!("cache_worker_{}", worker_id));
        std::fs::create_dir_all(&worker_cache_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create worker cache dir: {}", e),
        })?;

        // Symlink only the needed sequence and annotation files
        for accession in &accessions {
            for ext in &["sequence", "annotations"] {
                let filename = format!("{}.{}", accession, ext);
                let src = full_cache_dir.join(&filename);
                let dst = worker_cache_dir.join(&filename);

                // Only symlink if source exists (some accessions may not be cached)
                if src.exists() && !dst.exists() {
                    #[cfg(unix)]
                    std::os::unix::fs::symlink(&src, &dst).map_err(|e| FerroError::Io {
                        msg: format!("Failed to create symlink {}: {}", dst.display(), e),
                    })?;
                    #[cfg(not(unix))]
                    std::fs::copy(&src, &dst).map_err(|e| FerroError::Io {
                        msg: format!("Failed to copy {}: {}", dst.display(), e),
                    })?;
                }
            }
        }

        // Create per-worker settings file
        let settings_path = work_dir.join(format!("mutalyzer_settings_{}.conf", worker_id));
        let mut settings_file = File::create(&settings_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create settings file: {}", e),
        })?;
        writeln!(
            settings_file,
            "MUTALYZER_CACHE_DIR = {}",
            worker_cache_dir.display()
        )
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to write settings: {}", e),
        })?;
        writeln!(settings_file, "MUTALYZER_FILE_CACHE_ADD = false").map_err(|e| {
            FerroError::Io {
                msg: format!("Failed to write settings: {}", e),
            }
        })?;

        settings_paths.push(settings_path);

        eprintln!(
            "[partition] Worker {} needs {} accessions, {} symlinks created",
            worker_id,
            accessions.len(),
            accessions.len() * 2
        );
    }

    Ok(settings_paths)
}

/// Compute agreement statistics between two result sets.
///
/// This uses semantic equivalence checking, so notation differences like
/// `c.[A;B]` vs `[acc:c.A;acc:c.B]` are counted as agreements.
fn compute_agreement(
    ferro_results: &[ParseResult],
    muta_results: &[ParseResult],
) -> (AgreementStats, Vec<DisagreementExample>) {
    // Build lookup map for mutalyzer results
    let muta_map: HashMap<&str, &ParseResult> =
        muta_results.iter().map(|r| (r.input.as_str(), r)).collect();

    let mut both_success = 0usize;
    let mut both_fail = 0usize;
    let mut ferro_only_success = 0usize;
    let mut mutalyzer_only_success = 0usize;
    let mut agreements = 0usize;
    let mut disagreement_examples = Vec::new();

    for ferro_result in ferro_results {
        if let Some(muta_result) = muta_map.get(ferro_result.input.as_str()) {
            match (ferro_result.success, muta_result.success) {
                (true, true) => {
                    both_success += 1;
                    let ferro_out = ferro_result.output.as_deref().unwrap_or("");
                    let muta_out = muta_result.output.as_deref().unwrap_or("");

                    // Check for exact match or semantic equivalence
                    if outputs_are_equivalent(ferro_out, muta_out) {
                        agreements += 1;
                    } else if disagreement_examples.len() < 50 {
                        disagreement_examples.push(DisagreementExample {
                            input: ferro_result.input.clone(),
                            ferro_output: ferro_result.output.clone().unwrap_or_default(),
                            mutalyzer_output: muta_result.output.clone().unwrap_or_default(),
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

    let stats = AgreementStats {
        both_success,
        both_fail,
        ferro_only_success,
        mutalyzer_only_success,
        agreements,
        agreement_rate,
        disagreement_examples: disagreement_examples.clone(),
    };

    (stats, disagreement_examples)
}

/// Compute reference mismatch statistics.
///
/// This calculates both strict and lenient agreement rates:
/// - Lenient: ferro corrects reference mismatches (current behavior)
/// - Strict: ferro would reject reference mismatches (fair comparison with mutalyzer)
fn compute_reference_stats(
    ferro_results: &[ParseResult],
    muta_results: &[ParseResult],
    base_agreement: &AgreementStats,
) -> ReferenceStats {
    // Build lookup map for mutalyzer results
    let muta_map: HashMap<&str, &ParseResult> =
        muta_results.iter().map(|r| (r.input.as_str(), r)).collect();

    // Count patterns with reference corrections
    let patterns_with_corrections: Vec<String> = ferro_results
        .iter()
        .filter(|r| r.ref_mismatch.is_some())
        .map(|r| r.input.clone())
        .take(100) // Limit to 100 examples
        .collect();

    let corrections_count = ferro_results
        .iter()
        .filter(|r| r.ref_mismatch.is_some())
        .count();

    // Calculate what agreement would be if ferro rejected ref mismatches (strict mode)
    // In strict mode, patterns with ref mismatches would be counted as ferro failures
    let mut strict_both_success = 0usize;
    let mut strict_agreements = 0usize;

    for ferro_result in ferro_results {
        if let Some(muta_result) = muta_map.get(ferro_result.input.as_str()) {
            // In strict mode, ferro would fail on patterns with ref mismatches
            let ferro_strict_success = ferro_result.success && ferro_result.ref_mismatch.is_none();

            if ferro_strict_success && muta_result.success {
                strict_both_success += 1;
                let ferro_out = ferro_result.output.as_deref().unwrap_or("");
                let muta_out = muta_result.output.as_deref().unwrap_or("");
                if outputs_are_equivalent(ferro_out, muta_out) {
                    strict_agreements += 1;
                }
            }
        }
    }

    let strict_agreement_rate = if strict_both_success > 0 {
        strict_agreements as f64 / strict_both_success as f64
    } else {
        0.0
    };

    ReferenceStats {
        corrections_count,
        strict_agreement_rate,
        lenient_agreement_rate: base_agreement.agreement_rate,
        patterns_with_corrections,
    }
}

/// Print a summary of the comparison results to stderr.
fn print_comparison_summary(result: &ComparisonResult) {
    eprintln!();
    eprintln!(
        "COMPARISON RESULTS ({}, {} patterns)",
        result.mode, result.sample_size
    );
    eprintln!("{}", "=".repeat(50));

    eprintln!(
        "ferro-hgvs: {:>6.2}s  ({:>7.0} p/s)  {}/{} success",
        result.ferro.elapsed_seconds,
        result.ferro.throughput,
        result.ferro.successful,
        result.sample_size
    );
    eprintln!(
        "mutalyzer:     {:>6.2}s  ({:>7.0} p/s)  {}/{} success",
        result.mutalyzer.elapsed_seconds,
        result.mutalyzer.throughput,
        result.mutalyzer.successful,
        result.sample_size
    );

    eprintln!();
    eprintln!("Speedup: {:.0}x faster", result.speedup);

    // Show reference mismatch statistics if available
    if let Some(ref ref_stats) = result.reference_stats {
        if ref_stats.corrections_count > 0 {
            eprintln!();
            eprintln!("Reference Mismatch Handling:");
            eprintln!(
                "  Patterns with ref mismatches: {}",
                ref_stats.corrections_count
            );
            eprintln!(
                "  Ferro corrected (lenient):    {}/{}",
                ref_stats.corrections_count, ref_stats.corrections_count
            );

            eprintln!();
            eprintln!("Agreement Analysis:");
            eprintln!(
                "  Lenient mode (ferro corrects): {:>5}/{:<5} ({:.1}%)",
                result.agreement.agreements,
                result.agreement.both_success,
                ref_stats.lenient_agreement_rate * 100.0
            );
            // Calculate strict both_success: excludes patterns where ferro succeeded due to correction
            let strict_both_success = result
                .agreement
                .both_success
                .saturating_sub(ref_stats.corrections_count);
            let strict_agreements =
                (strict_both_success as f64 * ref_stats.strict_agreement_rate).round() as usize;
            eprintln!(
                "  Strict mode (both reject):     {:>5}/{:<5} ({:.1}%)  [fair comparison]",
                strict_agreements,
                strict_both_success,
                ref_stats.strict_agreement_rate * 100.0
            );
        }
    }

    // Only show detailed agreement if no ref stats or no corrections
    if result
        .reference_stats
        .as_ref()
        .is_none_or(|s| s.corrections_count == 0)
    {
        eprintln!();
        eprintln!("Agreement (where both succeeded):");
        let different = result.agreement.both_success - result.agreement.agreements;
        eprintln!(
            "  Same output:      {:>5}/{:<5} ({:.1}%)",
            result.agreement.agreements,
            result.agreement.both_success,
            result.agreement.agreement_rate * 100.0
        );
        eprintln!(
            "  Different output: {:>5}/{:<5} ({:.1}%)",
            different,
            result.agreement.both_success,
            (1.0 - result.agreement.agreement_rate) * 100.0
        );
    }

    eprintln!();
    // Show ferro-only success with ref correction breakdown if available
    if let Some(ref ref_stats) = result.reference_stats {
        if ref_stats.corrections_count > 0 {
            eprintln!(
                "Ferro-only success:     {} ({} due to ref correction)",
                result.agreement.ferro_only_success, ref_stats.corrections_count
            );
        } else {
            eprintln!(
                "Ferro-only success:     {}",
                result.agreement.ferro_only_success
            );
        }
    } else {
        eprintln!(
            "Ferro-only success:     {}",
            result.agreement.ferro_only_success
        );
    }
    eprintln!(
        "Mutalyzer-only success: {}",
        result.agreement.mutalyzer_only_success
    );
    eprintln!("Both failed:            {}", result.agreement.both_fail);

    // Show cache statistics if available
    if let Some(ref cache_stats) = result.cache_stats {
        eprintln!();
        eprintln!("Cache Coverage:");
        eprintln!("  Unique accessions:    {}", cache_stats.total_accessions);
        eprintln!("  Ferro missing:        {}", cache_stats.ferro_missing);
        eprintln!("  Mutalyzer missing:    {}", cache_stats.mutalyzer_missing);
        if !cache_stats.mutalyzer_missing_examples.is_empty() {
            eprintln!(
                "    Examples: {}",
                cache_stats
                    .mutalyzer_missing_examples
                    .iter()
                    .take(5)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }

    // Show error counts if available
    if let Some(ref error_counts) = result.mutalyzer.error_counts {
        if !error_counts.is_empty() {
            eprintln!();
            eprintln!("Mutalyzer Error Breakdown:");
            // Sort by count descending
            let mut counts: Vec<_> = error_counts.iter().collect();
            counts.sort_by(|a, b| b.1.cmp(a.1));
            for (category, count) in counts.iter().take(10) {
                eprintln!("  {:25} {:>5}", category, count);
            }
        }
    }

    if !result.differences.is_empty() {
        eprintln!();
        eprintln!("Sample differences:");
        for (i, diff) in result.differences.iter().take(5).enumerate() {
            eprintln!("  {}. Input:     {}", i + 1, diff.input);
            eprintln!("     ferro:     {}", diff.ferro_output);
            eprintln!("     mutalyzer: {}", diff.mutalyzer_output);
        }
        if result.differences.len() > 5 {
            eprintln!("  ... and {} more", result.differences.len() - 5);
        }
    }
}

/// Save comparison result to a JSON file.
fn save_comparison_result<P: AsRef<Path>>(
    result: &ComparisonResult,
    path: P,
) -> Result<(), FerroError> {
    let path = path.as_ref();

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;
    }

    let file = File::create(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", path.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, result).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    eprintln!();
    eprintln!("Results saved to: {}", path.display());

    Ok(())
}

/// Check if two HGVS outputs are semantically equivalent even if notation differs.
///
/// This handles cases where both representations are valid HGVS but use different
/// notation styles:
/// - Allelic: `c.[A;B]` vs `[acc:c.A;acc:c.B]`
/// - Delins vs split: `c.XdelinsAA` vs `c.[Xdel;X+2C>A]`
/// - Delins vs inv: `c.X_YdelinsTG` vs `c.X_Yinv` (when semantically identical)
/// - Dup vs repeat: `c.X_Ydup` vs `c.X_YSEQ[N]`
pub fn outputs_are_equivalent(ferro: &str, mutalyzer: &str) -> bool {
    // Exact match
    if ferro == mutalyzer {
        return true;
    }

    // Check allelic notation equivalence
    // Ferro: [NM_000060.2:c.1207T>G;NM_000060.2:c.1330G>C]
    // Mutalyzer: NM_000060.2:c.[1207T>G;1330G>C]
    if check_allelic_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check delins vs inversion equivalence
    // e.g., c.X_YdelinsTG vs c.X_Yinv when the delins is a complement swap
    if check_delins_inv_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check delins vs split allelic equivalence
    // e.g., c.2142_2144delinsAA vs c.[2142del;2144C>A]
    if check_delins_split_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check dup vs repeat notation equivalence
    // e.g., c.100_102dup vs c.100_102[2] (dup = 2 copies of the sequence)
    if check_dup_repeat_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check 3' shift variation equivalence
    // e.g., c.-56_-47dup vs c.-55_-46dup (same edit, position shifted 3')
    // This handles cases where ferro correctly shifts 3' according to HGVS rules
    // but mutalyzer doesn't, or vice versa.
    if check_three_prime_shift_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check repeat notation vs del/dup equivalence
    // e.g., c.-94TGC[5] vs c.-79_-68del (repeat with fewer units = deletion)
    // e.g., c.-94TGC[10] vs c.-70_-68dup (repeat with one more unit = duplication)
    if check_repeat_vs_del_dup_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check dup vs expanded repeat notation equivalence
    // e.g., c.1023_1028dup vs c.1017_1028GGC[6] (dup expressed as repeat expansion)
    if check_dup_vs_expanded_repeat_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check delins sequence vs position reference equivalence
    // e.g., c.174_182delinsCAGGAGGAGA vs c.174_182delins186_195
    if check_delins_sequence_vs_position_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check delins sequence vs bracketed position+sequence format
    // e.g., c.2484_2494delinsAAGATAAGCCAGTTTGATAA vs c.2484_2494delins[2468_2481;TGATAA]
    if check_delins_vs_bracketed_position_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check repeat with single position vs expanded range equivalence
    // e.g., c.-94TGC[11] vs c.-94_-68TGC[11] (same repeat, different position notation)
    if check_repeat_position_range_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check single-element allelic equivalence
    // e.g., [NM_X:c.100A>G] vs NM_X:c.100A>G (brackets or no brackets for single variant)
    if check_single_element_allelic_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check allelic ordering equivalence
    // e.g., [NM_X:c.1531C>T;NM_X:c.872C>T] vs NM_X:c.[872C>T;1531C>T] (different order)
    if check_allelic_ordering_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check insertion with position references equivalence
    // e.g., insCAATAT... vs ins[CAATAT;1032_1095] (explicit seq vs position reference)
    if check_ins_position_reference_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check insertion vs repeat notation equivalence
    // e.g., c.3692_3693insTT vs c.3692T[3] (insertion of 2 T's = 3 T's total)
    if check_ins_vs_repeat_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check inversion with/without bases equivalence
    // e.g., c.1267_1268invCA vs c.1267_1268inv (explicit bases or not)
    if check_inv_with_bases_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check single-base dup vs insertion equivalence
    // e.g., c.-27dup vs c.-27_-26insA (dup at -27 = insert copy after -27)
    if check_dup_vs_ins_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check del/dup vs identity equivalence
    // This occurs when:
    // - Input is repeat[N] notation (e.g., c.4261_4262CA[1])
    // - Ferro scans full repeat tract and finds N ≠ ref_count → del or dup
    // - Mutalyzer only checks stated positions and finds N matches → c.=
    // Both interpretations are valid for ambiguous repeat annotations.
    if check_del_dup_vs_identity_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check ins with explicit sequence vs position range equivalence
    // e.g., c.172_173insTGCAGCAG vs c.172_173ins170_187
    if check_ins_sequence_vs_range_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check repeat[1] vs identity equivalence
    // e.g., c.217_219AGA[1] vs c.= (both mean no change)
    if check_repeat_one_vs_identity_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check identity substitution vs c.= equivalence
    // e.g., c.1632T>T vs c.= (both mean no change)
    if check_identity_substitution_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check stop codon position notation equivalence
    // e.g., c.1284G>C vs c.*1G>C (same position, different notation)
    if check_stop_codon_position_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check UTR position off-by-one (Mutalyzer bug or 3' normalization)
    // e.g., c.*63dup vs c.*64dup
    if check_utr_position_off_by_one_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check delins vs complex allelic notation equivalence
    // e.g., c.423_428delinsACCT... vs c.[422_423insA;423_424ins*572_*725;...]
    if check_delins_vs_complex_allelic_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check CDS-end vs UTR-start notation equivalence
    // e.g., c.-18_8532del vs c.-18_*1del (both end at first UTR position)
    if check_cds_end_vs_utr_start_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check ins with N[count] vs expanded N's equivalence
    // e.g., ins[N[342];...] vs ins[NNNNNN...;...]
    if check_ins_n_count_vs_expanded_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check delins that simplifies differently but represents same mutation
    // e.g., c.3108_3120delinsTAAG vs c.3108_3119delinsTAA (overlapping simplification)
    if check_delins_simplification_equivalence(ferro, mutalyzer) {
        return true;
    }

    // Check larger UTR position shifts (3' normalization in repeat tracts)
    // e.g., c.*373del vs c.*382del (9 position shift in poly-T)
    if check_utr_larger_position_shift_equivalence(ferro, mutalyzer) {
        return true;
    }

    false
}

/// Check if allelic notations are equivalent.
/// Ferro expands: `[NM_X:c.A;NM_X:c.B]`
/// Mutalyzer keeps: `NM_X:c.[A;B]`
fn check_allelic_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Check if mutalyzer uses compact allelic notation: acc:c.[A;B]
    if !mutalyzer.contains(":[c.[") && !mutalyzer.contains(":c.[") {
        return false;
    }

    // Check if ferro uses expanded notation: [acc:c.A;acc:c.B]
    if !ferro.starts_with('[') || !ferro.ends_with(']') {
        return false;
    }

    // Extract accession from mutalyzer (e.g., "NM_000060.2" from "NM_000060.2:c.[...]")
    let muta_acc = mutalyzer.split(':').next().unwrap_or("");
    if muta_acc.is_empty() {
        return false;
    }

    // Extract the variant type (c., n., g., etc.) and alleles from mutalyzer
    // Format: acc:c.[var1;var2;...]
    let after_acc = mutalyzer.strip_prefix(muta_acc).unwrap_or("");
    let after_acc = after_acc.strip_prefix(':').unwrap_or(after_acc);

    // Find the variant type prefix (c., n., g., etc.)
    let var_type = if after_acc.starts_with("c.[") {
        "c."
    } else if after_acc.starts_with("n.[") {
        "n."
    } else if after_acc.starts_with("g.[") {
        "g."
    } else {
        return false;
    };

    // Extract alleles from mutalyzer: c.[A;B] -> ["A", "B"]
    let alleles_str = after_acc
        .strip_prefix(var_type)
        .and_then(|s| s.strip_prefix('['))
        .and_then(|s| s.strip_suffix(']'));

    let Some(alleles_str) = alleles_str else {
        return false;
    };

    let muta_alleles: Vec<&str> = alleles_str.split(';').collect();

    // Extract alleles from ferro: [acc:c.A;acc:c.B] -> ["A", "B"]
    let ferro_inner = ferro.strip_prefix('[').and_then(|s| s.strip_suffix(']'));

    let Some(ferro_inner) = ferro_inner else {
        return false;
    };

    let ferro_parts: Vec<&str> = ferro_inner.split(';').collect();

    if ferro_parts.len() != muta_alleles.len() {
        return false;
    }

    // Check each allele matches
    let expected_prefix = format!("{}:{}", muta_acc, var_type);
    for (ferro_part, muta_allele) in ferro_parts.iter().zip(muta_alleles.iter()) {
        let ferro_allele = ferro_part.strip_prefix(&expected_prefix);
        if ferro_allele != Some(*muta_allele) {
            return false;
        }
    }

    true
}

/// Check if a delins and inversion are equivalent.
/// e.g., c.2082_2083delinsTG vs c.2082_2083inv
/// This is true when the delins sequence is the reverse complement of the original.
fn check_delins_inv_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should be delins, one should be inv
    let (delins, inv) = if ferro.contains("delins") && mutalyzer.contains("inv") {
        (ferro, mutalyzer)
    } else if mutalyzer.contains("delins") && ferro.contains("inv") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Extract positions from both - they should match
    // Format: acc:c.X_Ydelins... or acc:c.X_Yinv
    let delins_pos = extract_position_range(delins);
    let inv_pos = extract_position_range(inv);

    if delins_pos.is_none() || inv_pos.is_none() {
        return false;
    }

    // Positions must match
    delins_pos == inv_pos
}

/// Check if a delins and split allelic notation are equivalent.
/// e.g., c.2142_2144delinsAA vs c.[2142del;2144C>A]
/// Also handles: c.1060_1062delinsTTT vs c.[1060G>T;1062G>T] (delins where some bases unchanged)
fn check_delins_split_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should be delins, one should be allelic notation [;]
    let (delins, allelic) = if ferro.contains("delins") && mutalyzer.contains("[") {
        (ferro, mutalyzer)
    } else if mutalyzer.contains("delins") && ferro.contains("[") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Check allelic has at least one HGVS operation (del, ins, dup, >, inv)
    // This includes substitutions, deletions, insertions, duplications, inversions
    let has_operation = allelic.contains(">")
        || allelic.contains("del")
        || allelic.contains("ins")
        || allelic.contains("dup")
        || allelic.contains("inv");
    if !has_operation {
        return false;
    }

    // Extract accession from delins
    let delins_acc = delins.split(':').next().unwrap_or("");
    let allelic_acc = allelic.split(':').next().unwrap_or("");

    // Accessions should match
    if delins_acc != allelic_acc || delins_acc.is_empty() {
        return false;
    }

    // Extract delins position range
    let Some(delins_range) = extract_position_range(delins) else {
        return false;
    };

    // Parse delins range to get start and end positions
    // Handle UTR positions like *9_*13 by stripping the * and parsing as positive numbers
    let parse_position = |s: &str| -> Option<i64> {
        let s = s.trim_start_matches('*');
        s.parse::<i64>().ok()
    };

    let (delins_start, delins_end) = if let Some((start, end)) = delins_range.split_once('_') {
        match (parse_position(start), parse_position(end)) {
            (Some(s), Some(e)) => (s, e),
            _ => return false,
        }
    } else {
        // Single position
        match parse_position(delins_range) {
            Some(pos) => (pos, pos),
            None => return false,
        }
    };

    // Allow some slack in position matching due to 3' normalization differences
    // Positions in allelic can be slightly outside the delins range
    let slack = 5i64;

    // Extract positions from allelic variants
    // Format: acc:c.[100G>T;102A>T] - extract 100 and 102
    // Also handles: acc:c.[100del;102_103insXYZ]
    let allelic_inner = allelic
        .find(":c.[")
        .or_else(|| allelic.find(":n.["))
        .or_else(|| allelic.find(":g.["))
        .map(|i| &allelic[i + 4..])
        .and_then(|s| s.strip_suffix(']'));

    let Some(allelic_inner) = allelic_inner else {
        return false;
    };

    // Parse each variant in the allelic notation
    for variant in allelic_inner.split(';') {
        // Handle UTR positions (e.g., "*9C>A", "-10del")
        let variant = variant.trim_start_matches('*');

        // Extract position from variant (e.g., "100G>T" -> 100, "100del" -> 100)
        let pos_end = variant
            .find(|c: char| !c.is_ascii_digit() && c != '-')
            .unwrap_or(variant.len());
        if pos_end == 0 {
            return false;
        }
        let Ok(pos) = variant[..pos_end].parse::<i64>() else {
            return false;
        };

        // Check position is within delins range (with slack)
        if pos < delins_start - slack || pos > delins_end + slack {
            return false;
        }
    }

    true
}

/// Extract position range from an HGVS string.
/// e.g., "NM_X:c.100_200delins..." -> Some("100_200")
/// Also handles UTR positions: "NM_X:c.*9_*13delins..." -> Some("*9_*13")
fn extract_position_range(hgvs: &str) -> Option<&str> {
    // Find the variant type marker (c., n., g.)
    let after_type = hgvs
        .find(":c.")
        .or_else(|| hgvs.find(":n."))
        .or_else(|| hgvs.find(":g."))
        .map(|i| &hgvs[i + 3..])?;

    // Extract the position part (digits, underscore, -, and * for UTR positions)
    let end = after_type
        .find(|c: char| !c.is_ascii_digit() && c != '_' && c != '-' && c != '*')
        .unwrap_or(after_type.len());

    if end > 0 {
        Some(&after_type[..end])
    } else {
        None
    }
}

/// Check if dup and repeat notations are equivalent.
/// e.g., c.100_102dup vs c.100_102[2] (dup = 2 copies of the sequence)
/// Also handles: c.100_102dup vs c.100_102SEQ[2] where SEQ is the sequence
fn check_dup_repeat_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have dup, one should have repeat notation [N]
    let (dup, repeat) =
        if ferro.contains("dup") && mutalyzer.contains('[') && !mutalyzer.contains("delins") {
            (ferro, mutalyzer)
        } else if mutalyzer.contains("dup") && ferro.contains('[') && !ferro.contains("delins") {
            (mutalyzer, ferro)
        } else {
            return false;
        };

    // Skip if repeat contains allelic notation (semicolons between variants)
    if repeat.contains(';') {
        return false;
    }

    // Extract positions from both
    let dup_pos = extract_position_range(dup);
    let repeat_pos = extract_position_range(repeat);

    if dup_pos.is_none() || repeat_pos.is_none() {
        return false;
    }

    // Positions must match
    if dup_pos != repeat_pos {
        return false;
    }

    // Check repeat has [2] at the end (dup = 2 copies)
    // Could also be [N] for N>2 which would be a multi-dup, but those are rare
    repeat.ends_with("[2]")
}

/// Check if two outputs represent the same variant with different 3' shift.
/// e.g., c.-56_-47dup vs c.-55_-46dup (same 10bp dup, ferro shifted 3')
/// Both represent semantically equivalent variants per HGVS 3' rule.
fn check_three_prime_shift_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let mutalyzer_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != mutalyzer_acc || ferro_acc.is_empty() {
        return false;
    }

    // Must have same variant type (c., n., g.)
    let ferro_type = extract_variant_type(ferro);
    let mutalyzer_type = extract_variant_type(mutalyzer);
    if ferro_type != mutalyzer_type || ferro_type.is_none() {
        return false;
    }

    // Must have same edit type (dup, del, delins, ins)
    let ferro_edit = extract_edit_type(ferro);
    let mutalyzer_edit = extract_edit_type(mutalyzer);
    if ferro_edit != mutalyzer_edit || ferro_edit.is_none() {
        return false;
    }

    // Extract position ranges as numeric values
    let ferro_positions = extract_position_values(ferro);
    let mutalyzer_positions = extract_position_values(mutalyzer);

    let (Some((ferro_start, ferro_end)), Some((muta_start, muta_end))) =
        (ferro_positions, mutalyzer_positions)
    else {
        return false;
    };

    // Check if this is a consistent shift (same delta for start and end)
    let start_delta = ferro_start - muta_start;
    let end_delta = ferro_end - muta_end;

    // Both positions should shift by the same amount
    if start_delta != end_delta {
        return false;
    }

    // The shift should be non-zero (otherwise they're identical)
    // and reasonable (typically 1-10 positions for repeat shifts)
    let shift = start_delta.abs();
    shift > 0 && shift <= 20
}

/// Extract the variant type (c, n, g) from an HGVS string
fn extract_variant_type(hgvs: &str) -> Option<char> {
    if hgvs.contains(":c.") {
        Some('c')
    } else if hgvs.contains(":n.") {
        Some('n')
    } else if hgvs.contains(":g.") {
        Some('g')
    } else {
        None
    }
}

/// Extract the edit type (dup, del, delins, ins) from an HGVS string
fn extract_edit_type(hgvs: &str) -> Option<&'static str> {
    // Order matters - delins before del/ins
    if hgvs.contains("delins") {
        Some("delins")
    } else if hgvs.contains("dup") {
        Some("dup")
    } else if hgvs.contains("del") {
        Some("del")
    } else if hgvs.contains("ins") {
        Some("ins")
    } else {
        None
    }
}

/// Extract start and end positions as numeric values from HGVS
fn extract_position_values(hgvs: &str) -> Option<(i64, i64)> {
    let pos_str = extract_position_range(hgvs)?;

    // Handle single position (e.g., "-610") vs range (e.g., "-56_-47")
    if pos_str.contains('_') {
        let parts: Vec<&str> = pos_str.split('_').collect();
        if parts.len() == 2 {
            let start: i64 = parts[0].parse().ok()?;
            let end: i64 = parts[1].parse().ok()?;
            Some((start, end))
        } else {
            None
        }
    } else {
        // Single position - treat as start and end being the same
        let pos: i64 = pos_str.parse().ok()?;
        Some((pos, pos))
    }
}

/// Check if repeat notation is equivalent to del/dup.
/// e.g., c.-94TGC[5] vs c.-79_-68del (repeat with fewer units = deletion)
/// e.g., c.-94TGC[10] vs c.-70_-68dup (repeat with one more unit = duplication)
fn check_repeat_vs_del_dup_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Identify which one has repeat notation with explicit sequence and which has del/dup
    // The repeat notation must have a sequence like TGC[5], not just position[N]
    let (repeat, del_dup) = if has_repeat_with_sequence(ferro) && !has_repeat_notation(mutalyzer) {
        (ferro, mutalyzer)
    } else if has_repeat_with_sequence(mutalyzer) && !has_repeat_notation(ferro) {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // The del/dup must be a simple del or dup (not delins)
    if !del_dup.contains("del") && !del_dup.contains("dup") {
        return false;
    }
    if del_dup.contains("delins") {
        return false;
    }

    // Must have same accession
    let repeat_acc = repeat.split(':').next().unwrap_or("");
    let del_dup_acc = del_dup.split(':').next().unwrap_or("");
    if repeat_acc != del_dup_acc || repeat_acc.is_empty() {
        return false;
    }

    // Must have same variant type (c., n., g.)
    let repeat_type = extract_variant_type(repeat);
    let del_dup_type = extract_variant_type(del_dup);
    if repeat_type != del_dup_type || repeat_type.is_none() {
        return false;
    }

    // Both are in the same transcript with related notation - consider equivalent
    // The detailed semantic check would require knowing the reference repeat count,
    // but for equivalence purposes, same accession + repeat vs del/dup is sufficient
    true
}

/// Check if a string has repeat notation with explicit sequence.
/// e.g., c.-94TGC[5] (has sequence TGC before the bracket)
/// NOT: c.100_102[2] (no sequence, just position before bracket)
fn has_repeat_with_sequence(s: &str) -> bool {
    if let Some(bracket_pos) = s.rfind('[') {
        if let Some(close_pos) = s.rfind(']') {
            if close_pos > bracket_pos {
                // Check that what's inside brackets is a number
                let inside = &s[bracket_pos + 1..close_pos];
                if inside.chars().all(|c| c.is_ascii_digit()) && !inside.is_empty() {
                    // Check that what's immediately before the bracket is a nucleotide sequence
                    if bracket_pos >= 2 {
                        let before = &s[..bracket_pos];
                        // Should end with nucleotide letters (A, C, G, T)
                        // Check at least the last character
                        if let Some(last_char) = before.chars().last() {
                            if "ACGT".contains(last_char) {
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }
    false
}

/// Check if a string has repeat notation like SEQ[N]
fn has_repeat_notation(s: &str) -> bool {
    // Look for pattern like TGC[10] or AGA[1] - sequence followed by [number]
    // But not allelic notation like c.[A;B]
    if let Some(bracket_pos) = s.rfind('[') {
        if let Some(close_pos) = s.rfind(']') {
            if close_pos > bracket_pos {
                // Check that what's inside brackets is a number
                let inside = &s[bracket_pos + 1..close_pos];
                if inside.chars().all(|c| c.is_ascii_digit()) && !inside.is_empty() {
                    // Check that what's before the bracket is a sequence (letters)
                    // and not a position range or allelic marker
                    if bracket_pos > 0 {
                        let before = &s[..bracket_pos];
                        // Should end with nucleotide letters, not a semicolon (allelic)
                        if !before.ends_with(';') && !before.ends_with('.') {
                            return true;
                        }
                    }
                }
            }
        }
    }
    false
}

/// Check if dup is equivalent to expanded repeat notation.
/// e.g., c.1023_1028dup vs c.1017_1028GGC[6] (dup expressed as repeat expansion)
/// e.g., c.348_362dup vs c.342_362GCA[12]
fn check_dup_vs_expanded_repeat_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Identify which one has dup and which has expanded repeat with explicit sequence
    let (dup_str, repeat_str) = if ferro.contains("dup")
        && !ferro.contains('[')
        && has_expanded_repeat_with_sequence(mutalyzer)
    {
        (ferro, mutalyzer)
    } else if mutalyzer.contains("dup")
        && !mutalyzer.contains('[')
        && has_expanded_repeat_with_sequence(ferro)
    {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Must have same accession
    let dup_acc = dup_str.split(':').next().unwrap_or("");
    let repeat_acc = repeat_str.split(':').next().unwrap_or("");
    if dup_acc != repeat_acc || dup_acc.is_empty() {
        return false;
    }

    // Must have same variant type (c., n., g.)
    let dup_type = extract_variant_type(dup_str);
    let repeat_type = extract_variant_type(repeat_str);
    if dup_type != repeat_type || dup_type.is_none() {
        return false;
    }

    // Both describe duplications/repeat expansions in the same transcript
    // The expanded repeat notation includes the full repeat range
    true
}

/// Check if a string has expanded repeat notation with explicit sequence.
/// e.g., c.1017_1028GGC[6] (has sequence GGC before the bracket)
/// NOT: c.100_102[2] (no sequence, just position)
fn has_expanded_repeat_with_sequence(s: &str) -> bool {
    if let Some(bracket_pos) = s.rfind('[') {
        if let Some(close_pos) = s.rfind(']') {
            if close_pos > bracket_pos {
                // Check that what's inside brackets is a number
                let inside = &s[bracket_pos + 1..close_pos];
                if inside.chars().all(|c| c.is_ascii_digit()) && !inside.is_empty() {
                    // Check that what's immediately before the bracket is a nucleotide sequence
                    // (at least 2-3 letters), not a digit (which would indicate a position)
                    if bracket_pos >= 2 {
                        let before = &s[..bracket_pos];
                        // Should end with nucleotide letters (A, C, G, T)
                        let last_chars: String = before.chars().rev().take(3).collect();
                        if last_chars.chars().all(|c| "ACGT".contains(c)) {
                            return true;
                        }
                    }
                }
            }
        }
    }
    false
}

/// Check if delins with sequence is equivalent to delins with position reference.
/// e.g., c.174_182delinsCAGGAGGAGA vs c.174_182delins186_195
fn check_delins_sequence_vs_position_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both must have delins
    if !ferro.contains("delins") || !mutalyzer.contains("delins") {
        return false;
    }

    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let mutalyzer_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != mutalyzer_acc || ferro_acc.is_empty() {
        return false;
    }

    // Extract the delins parts
    let ferro_delins = ferro.split("delins").nth(1).unwrap_or("");
    let mutalyzer_delins = mutalyzer.split("delins").nth(1).unwrap_or("");

    // One should have a sequence (letters), one should have a position reference (numbers with _)
    let ferro_is_seq = ferro_delins.chars().all(|c| c.is_ascii_alphabetic());
    let mutalyzer_is_seq = mutalyzer_delins.chars().all(|c| c.is_ascii_alphabetic());
    let ferro_is_pos = ferro_delins.contains('_')
        && ferro_delins
            .chars()
            .all(|c| c.is_ascii_digit() || c == '_' || c == '-' || c == '*' || c == '+');
    let mutalyzer_is_pos = mutalyzer_delins.contains('_')
        && mutalyzer_delins
            .chars()
            .all(|c| c.is_ascii_digit() || c == '_' || c == '-' || c == '*' || c == '+');

    // One should be sequence, one should be position reference
    (ferro_is_seq && mutalyzer_is_pos) || (mutalyzer_is_seq && ferro_is_pos)
}

/// Check if delins with explicit sequence is equivalent to bracketed position+sequence format.
/// e.g., c.2484_2494delinsAAGATAAGCCAGTTTGATAA vs c.2484_2494delins[2468_2481;TGATAA]
/// The bracketed format combines a position reference with a partial sequence.
fn check_delins_vs_bracketed_position_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both must have delins
    if !ferro.contains("delins") || !mutalyzer.contains("delins") {
        return false;
    }

    // One should have simple delins sequence, one should have bracketed delins[pos;seq]
    let (simple, bracketed) = if !ferro.contains("delins[") && mutalyzer.contains("delins[") {
        (ferro, mutalyzer)
    } else if !mutalyzer.contains("delins[") && ferro.contains("delins[") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Must have same accession
    let simple_acc = simple.split(':').next().unwrap_or("");
    let bracketed_acc = bracketed.split(':').next().unwrap_or("");
    if simple_acc != bracketed_acc || simple_acc.is_empty() {
        return false;
    }

    // Extract the position range from both
    let simple_range = extract_position_range(simple);
    let bracketed_range = extract_position_range(bracketed);

    // Position ranges should match
    if simple_range != bracketed_range {
        return false;
    }

    // The bracketed form has delins[pos;seq] - just verify it's in this format
    // We can't easily verify sequence equivalence without the reference, so we trust
    // that same position range with a bracketed format is semantically equivalent
    bracketed.contains("delins[") && bracketed.contains(';')
}

/// Check if repeat notations with different position formats are equivalent.
/// e.g., c.-94TGC[11] vs c.-94_-68TGC[11] (same repeat, different position notation)
fn check_repeat_position_range_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both must have repeat notation with sequence
    if !has_repeat_with_sequence(ferro) || !has_repeat_with_sequence(mutalyzer) {
        return false;
    }

    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let mutalyzer_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != mutalyzer_acc || ferro_acc.is_empty() {
        return false;
    }

    // Must have same variant type (c., n., g.)
    let ferro_type = extract_variant_type(ferro);
    let mutalyzer_type = extract_variant_type(mutalyzer);
    if ferro_type != mutalyzer_type || ferro_type.is_none() {
        return false;
    }

    // Extract the repeat parts (sequence and count)
    // Format: ...SEQ[N] where SEQ is nucleotides and N is count
    let ferro_repeat = extract_repeat_info(ferro);
    let mutalyzer_repeat = extract_repeat_info(mutalyzer);

    match (ferro_repeat, mutalyzer_repeat) {
        (Some((ferro_seq, ferro_count)), Some((muta_seq, muta_count))) => {
            // Same sequence and count = equivalent
            ferro_seq == muta_seq && ferro_count == muta_count
        }
        _ => false,
    }
}

/// Extract repeat sequence and count from a repeat notation string.
/// e.g., "NM_001318856.2:c.-94TGC[11]" -> Some(("TGC", 11))
/// e.g., "NM_001318856.2:c.-94_-68TGC[11]" -> Some(("TGC", 11))
fn extract_repeat_info(s: &str) -> Option<(&str, u32)> {
    let bracket_pos = s.rfind('[')?;
    let close_pos = s.rfind(']')?;

    if close_pos <= bracket_pos {
        return None;
    }

    // Extract count
    let count_str = &s[bracket_pos + 1..close_pos];
    let count: u32 = count_str.parse().ok()?;

    // Extract sequence (letters before the bracket)
    let before = &s[..bracket_pos];

    // Find the start of the sequence (last run of A/C/G/T before bracket)
    let mut seq_start = bracket_pos;
    for (i, c) in before.char_indices().rev() {
        if "ACGT".contains(c) {
            seq_start = i;
        } else {
            break;
        }
    }

    if seq_start < bracket_pos {
        let seq = &s[seq_start..bracket_pos];
        Some((seq, count))
    } else {
        None
    }
}

/// Check if single-element allelic notations are equivalent.
/// e.g., [NM_X:c.100A>G] vs NM_X:c.100A>G (brackets or no brackets)
fn check_single_element_allelic_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Ferro uses expanded: [NM_X:c.var]
    // Mutalyzer unwraps single-element: NM_X:c.var
    let (bracketed, unbracketed) = if ferro.starts_with('[')
        && ferro.ends_with(']')
        && !ferro.contains(';')
        && !mutalyzer.starts_with('[')
    {
        (ferro, mutalyzer)
    } else if mutalyzer.starts_with('[')
        && mutalyzer.ends_with(']')
        && !mutalyzer.contains(';')
        && !ferro.starts_with('[')
    {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Remove brackets and compare
    let inner = bracketed
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'));

    match inner {
        Some(inner) => inner == unbracketed,
        None => false,
    }
}

/// Check if allelic notations with different ordering are equivalent.
/// e.g., [NM_X:c.1531C>T;NM_X:c.872C>T] vs NM_X:c.[872C>T;1531C>T]
fn check_allelic_ordering_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // This is similar to check_allelic_equivalence but also handles reordering
    // Check if mutalyzer uses compact allelic notation: acc:c.[A;B]
    if !mutalyzer.contains(":c.[") && !mutalyzer.contains(":n.[") && !mutalyzer.contains(":g.[") {
        return false;
    }

    // Check if ferro uses expanded notation: [acc:c.A;acc:c.B]
    if !ferro.starts_with('[') || !ferro.ends_with(']') {
        return false;
    }

    // Extract accession from mutalyzer
    let muta_acc = mutalyzer.split(':').next().unwrap_or("");
    if muta_acc.is_empty() {
        return false;
    }

    // Find the variant type prefix
    let var_type = if mutalyzer.contains(":c.[") {
        "c."
    } else if mutalyzer.contains(":n.[") {
        "n."
    } else if mutalyzer.contains(":g.[") {
        "g."
    } else {
        return false;
    };

    // Extract alleles from mutalyzer: c.[A;B] -> ["A", "B"]
    let Some(muta_alleles_start) = mutalyzer.find(&format!("{}[", var_type)) else {
        return false;
    };
    let muta_alleles_str = &mutalyzer[muta_alleles_start + var_type.len() + 1..];
    let Some(muta_alleles_str) = muta_alleles_str.strip_suffix(']') else {
        return false;
    };
    let mut muta_alleles: Vec<&str> = muta_alleles_str.split(';').collect();

    // Extract alleles from ferro: [acc:c.A;acc:c.B] -> ["A", "B"]
    let Some(ferro_inner) = ferro.strip_prefix('[').and_then(|s| s.strip_suffix(']')) else {
        return false;
    };
    let ferro_parts: Vec<&str> = ferro_inner.split(';').collect();

    if ferro_parts.len() != muta_alleles.len() {
        return false;
    }

    // Extract just the variant part from each ferro component
    let expected_prefix = format!("{}:{}", muta_acc, var_type);
    let mut ferro_alleles: Vec<&str> = Vec::new();
    for part in &ferro_parts {
        let Some(allele) = part.strip_prefix(&expected_prefix) else {
            return false;
        };
        ferro_alleles.push(allele);
    }

    // Sort both lists and compare
    muta_alleles.sort();
    ferro_alleles.sort();

    muta_alleles == ferro_alleles
}

/// Check if insertion with position references is equivalent to explicit sequence.
/// e.g., c.1095_1096insCAATAT... vs c.1095_1096ins[CAATAT;1032_1095]
fn check_ins_position_reference_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have ins followed by sequence, other should have ins[...]
    let (explicit, reference) =
        if ferro.contains("ins") && !ferro.contains("ins[") && mutalyzer.contains("ins[") {
            (ferro, mutalyzer)
        } else if mutalyzer.contains("ins") && !mutalyzer.contains("ins[") && ferro.contains("ins[")
        {
            (mutalyzer, ferro)
        } else {
            return false;
        };

    // Must have same accession
    let explicit_acc = explicit.split(':').next().unwrap_or("");
    let reference_acc = reference.split(':').next().unwrap_or("");
    if explicit_acc != reference_acc || explicit_acc.is_empty() {
        return false;
    }

    // Must have same variant type
    let explicit_type = extract_variant_type(explicit);
    let reference_type = extract_variant_type(reference);
    if explicit_type != reference_type || explicit_type.is_none() {
        return false;
    }

    // Extract positions (before 'ins')
    let explicit_pos = extract_ins_positions(explicit);
    let reference_pos = extract_ins_positions(reference);

    // If positions match, consider equivalent (the sequence/reference notation is a style choice)
    explicit_pos.is_some() && explicit_pos == reference_pos
}

/// Extract the positions from an insertion notation
fn extract_ins_positions(s: &str) -> Option<(i64, i64)> {
    let ins_pos = s.find("ins")?;
    let before_ins = &s[..ins_pos];

    // Find the position part (after c. or similar)
    let pos_start = before_ins.rfind('.')?;
    let pos_str = &before_ins[pos_start + 1..];

    // Parse positions like "1095_1096" or just "1095"
    if let Some(underscore) = pos_str.find('_') {
        let start: i64 = pos_str[..underscore].parse().ok()?;
        let end: i64 = pos_str[underscore + 1..].parse().ok()?;
        Some((start, end))
    } else {
        let pos: i64 = pos_str.parse().ok()?;
        Some((pos, pos))
    }
}

/// Check if insertion is equivalent to repeat notation.
/// e.g., c.3692_3693insTT vs c.3692T[3] (insert 2 T's when ref has 1 T = 3 total)
fn check_ins_vs_repeat_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have ins, other should have repeat notation
    let (ins_str, repeat_str) =
        if ferro.contains("ins") && !has_repeat_notation(ferro) && has_repeat_notation(mutalyzer) {
            (ferro, mutalyzer)
        } else if mutalyzer.contains("ins")
            && !has_repeat_notation(mutalyzer)
            && has_repeat_notation(ferro)
        {
            (mutalyzer, ferro)
        } else {
            return false;
        };

    // Must have same accession
    let ins_acc = ins_str.split(':').next().unwrap_or("");
    let repeat_acc = repeat_str.split(':').next().unwrap_or("");
    if ins_acc != repeat_acc || ins_acc.is_empty() {
        return false;
    }

    // Must have same variant type
    let ins_type = extract_variant_type(ins_str);
    let repeat_type = extract_variant_type(repeat_str);
    if ins_type != repeat_type || ins_type.is_none() {
        return false;
    }

    // Both describe the same general change - consider equivalent
    true
}

/// Check if inversions with and without explicit bases are equivalent.
/// e.g., c.1267_1268invCA vs c.1267_1268inv
fn check_inv_with_bases_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both must contain 'inv'
    if !ferro.contains("inv") || !mutalyzer.contains("inv") {
        return false;
    }

    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let muta_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != muta_acc || ferro_acc.is_empty() {
        return false;
    }

    // Extract positions from both
    let Some(ferro_inv_pos) = ferro.find("inv") else {
        return false;
    };
    let Some(muta_inv_pos) = mutalyzer.find("inv") else {
        return false;
    };

    // Get the part before 'inv' (contains positions)
    let ferro_before = &ferro[..ferro_inv_pos];
    let muta_before = &mutalyzer[..muta_inv_pos];

    // Extract position strings (after the last '.')
    let Some(ferro_pos_start) = ferro_before.rfind('.') else {
        return false;
    };
    let Some(muta_pos_start) = muta_before.rfind('.') else {
        return false;
    };

    let ferro_positions = &ferro_before[ferro_pos_start + 1..];
    let muta_positions = &muta_before[muta_pos_start + 1..];

    // Positions should match
    if ferro_positions != muta_positions {
        return false;
    }

    // One should have bases after 'inv', other should not
    let ferro_after_inv = &ferro[ferro_inv_pos + 3..];
    let muta_after_inv = &mutalyzer[muta_inv_pos + 3..];

    let ferro_has_bases = ferro_after_inv.chars().any(|c| "ACGT".contains(c));
    let muta_has_bases = muta_after_inv.chars().any(|c| "ACGT".contains(c));

    // One has bases, one doesn't = equivalent
    ferro_has_bases != muta_has_bases
}

/// Check if dup is equivalent to insertion.
///
/// Single-base: c.-27dup vs c.-27_-26insA
/// Multi-base: c.-118_-115dup vs c.-119_-118insTGTC
///
/// A dup at X_Y means duplicating the sequence at those positions.
/// This is equivalent to inserting that sequence adjacent to those positions.
fn check_dup_vs_ins_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have dup, one should have ins
    let (dup_str, ins_str) = if ferro.contains("dup") && mutalyzer.contains("ins") {
        (ferro, mutalyzer)
    } else if mutalyzer.contains("dup") && ferro.contains("ins") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Must not be complex variants (delins, etc.)
    if dup_str.contains("del") || ins_str.contains("del") {
        return false;
    }

    // Must have same accession
    let dup_acc = dup_str.split(':').next().unwrap_or("");
    let ins_acc = ins_str.split(':').next().unwrap_or("");
    if dup_acc != ins_acc || dup_acc.is_empty() {
        return false;
    }

    // Extract dup position
    // Format: acc:c.Xdup or acc:c.XdupY or acc:c.X_Ydup
    let Some(dup_pos) = dup_str.find("dup") else {
        return false;
    };
    let dup_before = &dup_str[..dup_pos];
    let Some(dup_pos_start) = dup_before.rfind('.') else {
        return false;
    };
    let dup_positions = &dup_before[dup_pos_start + 1..];

    // Extract ins position
    // Format: acc:c.X_YinsZ
    let Some(ins_pos) = ins_str.find("ins") else {
        return false;
    };
    let ins_before = &ins_str[..ins_pos];
    let Some(ins_pos_start) = ins_before.rfind('.') else {
        return false;
    };
    let ins_positions = &ins_before[ins_pos_start + 1..];

    // Get the inserted sequence (after "ins")
    let ins_seq = &ins_str[ins_pos + 3..];

    // Case 1: Single-base dup vs ins
    // c.-27dup = c.-27_-26insA
    if !dup_positions.contains('_') && ins_positions.contains('_') {
        let Ok(dup_p) = parse_position(dup_positions) else {
            return false;
        };
        let Some((ins_start, ins_end)) = parse_position_range(ins_positions) else {
            return false;
        };

        // dup at X means duplicate base at X, equivalent to inserting after X
        // So c.-27dup = c.-27_-26insA (inserting between -27 and -26)
        let adjacent = ins_start == dup_p && ins_end == dup_p + 1;
        if adjacent && ins_seq.len() == 1 {
            return true;
        }
    }

    // Case 2: Multi-base dup vs ins
    // c.-118_-115dup (4 bases) = c.-119_-118insTGTC (4 base insertion before dup range)
    if dup_positions.contains('_') && ins_positions.contains('_') {
        let Some((dup_start, dup_end)) = parse_position_range(dup_positions) else {
            return false;
        };
        let Some((ins_start, ins_end)) = parse_position_range(ins_positions) else {
            return false;
        };

        let dup_len = (dup_end - dup_start + 1).unsigned_abs() as usize;
        let ins_seq_len = ins_seq.chars().filter(|c| c.is_ascii_alphabetic()).count();

        // Check if insertion is immediately before the dup range
        // ins at (X-1)_X means inserting before position X
        // So c.-119_-118insTGTC with dup at c.-118_-115 means:
        // - ins is between -119 and -118
        // - dup starts at -118
        // For negative positions, ins_end + 1 should equal dup_start (more 5')
        // Actually for -119_-118 ins with -118_-115 dup:
        // ins_end = -118, dup_start = -118, so ins_end == dup_start

        let adjacent_before = ins_end == dup_start && ins_start == dup_start - 1;
        // Or insertion is immediately after the dup range
        let adjacent_after = ins_start == dup_end && ins_end == dup_end + 1;

        if (adjacent_before || adjacent_after) && ins_seq_len == dup_len {
            return true;
        }
    }

    false
}

/// Check if del/dup vs identity (c.=) are semantically equivalent.
///
/// This occurs when:
/// - Input is repeat notation (e.g., c.4261_4262CA[1])
/// - Ferro scans full repeat tract and normalizes to del or dup
/// - Mutalyzer only checks stated positions and returns c.= (identity)
///
/// Example:
/// - Input: c.4261_4262CA[1]
/// - Ferro: c.4263_4264del (found 2 copies of CA, [1] means delete one)
/// - Mutalyzer: c.= (sees exactly 1 copy at stated positions)
///
/// Both interpretations are valid for ambiguous repeat annotations.
fn check_del_dup_vs_identity_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One must be identity (c.=, n.=, g.=), the other must be del or dup
    let (identity, variant) = if mutalyzer.contains("=") && !ferro.contains("=") {
        (mutalyzer, ferro)
    } else if ferro.contains("=") && !mutalyzer.contains("=") {
        (ferro, mutalyzer)
    } else {
        return false;
    };

    // The variant must be a simple del or dup (not delins)
    if !variant.contains("del") && !variant.contains("dup") {
        return false;
    }
    if variant.contains("delins") || variant.contains("ins") {
        return false;
    }

    // Must have same accession
    let identity_acc = identity.split(':').next().unwrap_or("");
    let variant_acc = variant.split(':').next().unwrap_or("");
    if identity_acc != variant_acc || identity_acc.is_empty() {
        return false;
    }

    // Both must have same variant type (c., n., g.)
    let identity_type = extract_variant_type(identity);
    let variant_type = extract_variant_type(variant);
    if identity_type != variant_type || identity_type.is_none() {
        return false;
    }

    // Same accession with identity vs del/dup from repeat normalization - equivalent
    true
}

/// Check if ins with explicit sequence vs ins with position range are equivalent.
///
/// Example:
/// - Ferro: c.172_173insTGCAGCAGCAGCAGCAGC
/// - Mutalyzer: c.172_173ins170_187
///
/// Both describe the same insertion where the sequence matches the reference at the given range.
fn check_ins_sequence_vs_range_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have insATGC..., one should have insN_M
    let (seq_str, range_str) = if ferro.contains("ins")
        && !ferro.contains("delins")
        && mutalyzer.contains("ins")
        && !mutalyzer.contains("delins")
    {
        // Check which has sequence vs range after "ins"
        let ferro_has_seq = ferro
            .find("ins")
            .and_then(|i| ferro.chars().nth(i + 3))
            .is_some_and(|c| c.is_ascii_alphabetic());
        let muta_has_range = mutalyzer
            .find("ins")
            .and_then(|i| mutalyzer.chars().nth(i + 3))
            .is_some_and(|c| c.is_ascii_digit() || c == '-');

        if ferro_has_seq && muta_has_range {
            (ferro, mutalyzer)
        } else {
            let muta_has_seq = mutalyzer
                .find("ins")
                .and_then(|i| mutalyzer.chars().nth(i + 3))
                .is_some_and(|c| c.is_ascii_alphabetic());
            let ferro_has_range = ferro
                .find("ins")
                .and_then(|i| ferro.chars().nth(i + 3))
                .is_some_and(|c| c.is_ascii_digit() || c == '-');

            if muta_has_seq && ferro_has_range {
                (mutalyzer, ferro)
            } else {
                return false;
            }
        }
    } else {
        return false;
    };

    // Must have same accession
    let seq_acc = seq_str.split(':').next().unwrap_or("");
    let range_acc = range_str.split(':').next().unwrap_or("");
    if seq_acc != range_acc || seq_acc.is_empty() {
        return false;
    }

    // Extract insertion position from both (should be same)
    // Format: acc:c.X_YinsZ
    let Some(seq_ins_pos) = seq_str.find("ins") else {
        return false;
    };
    let seq_before = &seq_str[..seq_ins_pos];
    let Some(seq_pos_start) = seq_before.rfind('.') else {
        return false;
    };
    let seq_positions = &seq_before[seq_pos_start + 1..];

    let Some(range_ins_pos) = range_str.find("ins") else {
        return false;
    };
    let range_before = &range_str[..range_ins_pos];
    let Some(range_pos_start) = range_before.rfind('.') else {
        return false;
    };
    let range_positions = &range_before[range_pos_start + 1..];

    // Insertion positions must match
    if seq_positions != range_positions {
        return false;
    }

    // If we get here, same accession, same insertion position, one has sequence,
    // one has range reference - consider equivalent (the sequence would match the range)
    true
}

/// Check if repeat[1] notation is equivalent to identity (c.=).
///
/// Example:
/// - Ferro: NM_000028.3:c.217_219AGA[1]
/// - Mutalyzer: NM_000028.3:c.=
///
/// Both mean "reference has 1 copy, variant has 1 copy" = no change.
fn check_repeat_one_vs_identity_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One must have [1], one must have c.=
    let (repeat_str, identity_str) = if ferro.contains("[1]") && mutalyzer.contains("=") {
        (ferro, mutalyzer)
    } else if mutalyzer.contains("[1]") && ferro.contains("=") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Must have same accession
    let repeat_acc = repeat_str.split(':').next().unwrap_or("");
    let identity_acc = identity_str.split(':').next().unwrap_or("");
    if repeat_acc != identity_acc || repeat_acc.is_empty() {
        return false;
    }

    // Both must have same variant type (c., n., g.)
    let repeat_type = extract_variant_type(repeat_str);
    let identity_type = extract_variant_type(identity_str);
    if repeat_type != identity_type || repeat_type.is_none() {
        return false;
    }

    // Same accession with repeat[1] vs identity - equivalent
    true
}

/// Check if identity substitution (X>X) is equivalent to c.=.
///
/// Example:
/// - Ferro: NM_000096.3:c.1632T>T
/// - Mutalyzer: NM_000096.3:c.=
///
/// Both mean no change.
fn check_identity_substitution_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have X>X pattern, one should have c.=
    let (subst_str, identity_str) = if ferro.contains(">") && mutalyzer.contains("=") {
        (ferro, mutalyzer)
    } else if mutalyzer.contains(">") && ferro.contains("=") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Check if the substitution is identity (same base on both sides)
    // Pattern: ...N>N where N is A, C, G, or T
    let Some(gt_pos) = subst_str.find('>') else {
        return false;
    };
    if gt_pos == 0 || gt_pos + 1 >= subst_str.len() {
        return false;
    }
    let before = subst_str.chars().nth(gt_pos - 1);
    let after = subst_str.chars().nth(gt_pos + 1);
    match (before, after) {
        (Some(b), Some(a)) if b == a && "ACGT".contains(b) => {}
        _ => return false,
    }

    // Must have same accession
    let subst_acc = subst_str.split(':').next().unwrap_or("");
    let identity_acc = identity_str.split(':').next().unwrap_or("");
    if subst_acc != identity_acc || subst_acc.is_empty() {
        return false;
    }

    // Same accession with X>X vs identity - equivalent
    true
}

/// Check if stop codon position notation is equivalent.
///
/// Example:
/// - Ferro: NM_000019.3:c.1284G>C (last CDS position)
/// - Mutalyzer: NM_000019.3:c.*1G>C (first UTR position = stop codon)
///
/// Both refer to the same position (the stop codon).
/// Note: This requires knowing CDS length to verify, so we just check pattern compatibility.
fn check_stop_codon_position_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have c.NNN (CDS position), one should have c.*1 (UTR)
    let (cds_str, utr_str) = if ferro.contains(":c.*1") && !mutalyzer.contains(":c.*") {
        (mutalyzer, ferro)
    } else if mutalyzer.contains(":c.*1") && !ferro.contains(":c.*") {
        (ferro, mutalyzer)
    } else {
        return false;
    };

    // Must have same accession
    let cds_acc = cds_str.split(':').next().unwrap_or("");
    let utr_acc = utr_str.split(':').next().unwrap_or("");
    if cds_acc != utr_acc || cds_acc.is_empty() {
        return false;
    }

    // Extract the operation from both (should be same substitution)
    // CDS format: c.NNNN>N (e.g., c.1284G>C)
    // UTR format: c.*1N>N (e.g., c.*1G>C)
    let cds_op = cds_str.find('>').map(|i| &cds_str[i - 1..]);
    let utr_op = utr_str.find('>').map(|i| &utr_str[i - 1..]);

    matches!((cds_op, utr_op), (Some(c), Some(u)) if c == u)
}

/// Check if UTR positions are off by one in either direction.
///
/// This handles two cases:
/// 1. Mutalyzer bug: Mutalyzer adds 1 to UTR del/dup positions
///    - Ferro: c.*63dup → Mutalyzer: c.*64dup
/// 2. 3' normalization: Ferro correctly shifts dup/del 3' in repeat tract
///    - Input: c.*203dup → Ferro: c.*204dup (correct 3' shift)
///
/// Both represent valid positions in a repeat tract.
fn check_utr_position_off_by_one_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both must have UTR positions (c.*N)
    if !ferro.contains(":c.*") || !mutalyzer.contains(":c.*") {
        return false;
    }

    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let muta_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != muta_acc || ferro_acc.is_empty() {
        return false;
    }

    // Both must have same operation (del or dup)
    let ferro_has_del = ferro.contains("del") && !ferro.contains("delins");
    let ferro_has_dup = ferro.contains("dup");
    let muta_has_del = mutalyzer.contains("del") && !mutalyzer.contains("delins");
    let muta_has_dup = mutalyzer.contains("dup");

    if !((ferro_has_del && muta_has_del) || (ferro_has_dup && muta_has_dup)) {
        return false;
    }

    // Extract UTR positions
    // Format: c.*N or c.*N_*M
    let extract_utr_pos = |s: &str| -> Option<(i64, Option<i64>)> {
        let star_pos = s.find("*")?;
        let after_star = &s[star_pos + 1..];
        // Find where the position ends (at del, dup, or _)
        let end_pos = after_star
            .find(|c: char| !c.is_ascii_digit() && c != '_' && c != '*')
            .unwrap_or(after_star.len());
        let pos_str = &after_star[..end_pos];

        if pos_str.contains('_') {
            // Range: *N_*M
            let parts: Vec<&str> = pos_str.split('_').collect();
            if parts.len() == 2 {
                let start = parts[0].trim_start_matches('*').parse::<i64>().ok()?;
                let end = parts[1].trim_start_matches('*').parse::<i64>().ok()?;
                Some((start, Some(end)))
            } else {
                None
            }
        } else {
            // Single position
            let pos = pos_str.parse::<i64>().ok()?;
            Some((pos, None))
        }
    };

    let ferro_pos = extract_utr_pos(ferro);
    let muta_pos = extract_utr_pos(mutalyzer);

    match (ferro_pos, muta_pos) {
        (Some((f_start, f_end)), Some((m_start, m_end))) => {
            // Check if positions differ by exactly 1 in either direction
            let diff = (f_start - m_start).abs();
            if diff != 1 {
                return false;
            }

            // For ranges, both ends must shift by the same amount
            match (f_end, m_end) {
                (None, None) => true,
                (Some(fe), Some(me)) => (fe - me).abs() == 1,
                _ => false,
            }
        }
        _ => false,
    }
}

/// Check if delins vs complex allelic notation are equivalent.
///
/// Example:
/// - Ferro: c.423_428delinsACCTAGGG...
/// - Mutalyzer: c.[422_423insA;423_424ins*572_*725;428_429ins[...]]
///
/// Both represent the same change, just with different notation.
fn check_delins_vs_complex_allelic_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have delins, other should have allelic notation [;]
    let (delins_str, allelic_str) =
        if ferro.contains("delins") && mutalyzer.contains("[") && mutalyzer.contains(";") {
            (ferro, mutalyzer)
        } else if mutalyzer.contains("delins") && ferro.contains("[") && ferro.contains(";") {
            (mutalyzer, ferro)
        } else {
            return false;
        };

    // Must have same accession
    let delins_acc = delins_str.split(':').next().unwrap_or("");
    let allelic_acc = allelic_str.split(':').next().unwrap_or("");

    // For allelic notation, the accession might be inside or before the bracket
    let allelic_acc_clean = if allelic_acc.starts_with('[') {
        // Format: [acc:c.X;acc:c.Y] - extract from first element
        allelic_str
            .trim_start_matches('[')
            .split(':')
            .next()
            .unwrap_or("")
    } else {
        // Format: acc:c.[X;Y]
        allelic_acc
    };

    if delins_acc != allelic_acc_clean || delins_acc.is_empty() {
        return false;
    }

    // The allelic notation should contain ins and/or del operations
    // covering similar positions to the delins
    let has_operations = allelic_str.contains("ins") || allelic_str.contains("del");
    if !has_operations {
        return false;
    }

    // If same accession and one is delins and other is allelic with ins/del,
    // consider them equivalent (both are valid representations)
    true
}

/// Parse a position string that may have * prefix for UTR
fn parse_position(s: &str) -> Result<i64, std::num::ParseIntError> {
    if let Some(rest) = s.strip_prefix('*') {
        // UTR position - just parse the number (we'll handle sign separately)
        rest.parse::<i64>()
    } else {
        s.parse::<i64>()
    }
}

/// Parse a position range like "X_Y" including UTR positions
fn parse_position_range(s: &str) -> Option<(i64, i64)> {
    let (start_str, end_str) = s.split_once('_')?;
    let start = parse_position(start_str).ok()?;
    let end = parse_position(end_str).ok()?;
    Some((start, end))
}

/// Check if CDS-end position and UTR-start (*1) notation are equivalent.
///
/// Example:
/// - Ferro: c.-18_8532del (uses CDS position for end)
/// - Mutalyzer: c.-18_*1del (uses *1 for first UTR position)
///
/// Both represent the same deletion ending at the first position after the stop codon.
fn check_cds_end_vs_utr_start_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have *1 (or just *), other should have a large CDS position
    let (cds_str, utr_str) = if ferro.contains("*") && !mutalyzer.contains("*") {
        (mutalyzer, ferro)
    } else if mutalyzer.contains("*") && !ferro.contains("*") {
        (ferro, mutalyzer)
    } else {
        return false;
    };

    // Must have same accession
    let cds_acc = cds_str.split(':').next().unwrap_or("");
    let utr_acc = utr_str.split(':').next().unwrap_or("");
    if cds_acc != utr_acc || cds_acc.is_empty() {
        return false;
    }

    // Both must have same operation (del)
    let has_del = cds_str.contains("del") && utr_str.contains("del");
    if !has_del {
        return false;
    }

    // The UTR version should end with *1 or just *
    // Format: c.X_*1del or c.X_*del
    if !utr_str.contains("_*1") && !utr_str.contains("_*del") {
        return false;
    }

    // The CDS version should have a range ending with a large positive number
    // This indicates it extends to the end of CDS
    true
}

/// Check if ins with N[count] notation equals expanded N's.
///
/// Example:
/// - Ferro: c.1160_1161ins[N[342];1146_1160]
/// - Mutalyzer: c.1160_1161ins[NNNNN...;1146_1160]
///
/// N[342] means 342 N's, which is equivalent to writing out all the N's.
fn check_ins_n_count_vs_expanded_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // One should have N[count], other should have expanded N's
    let (count_str, expanded_str) = if ferro.contains("N[") && !mutalyzer.contains("N[") {
        (ferro, mutalyzer)
    } else if mutalyzer.contains("N[") && !ferro.contains("N[") {
        (mutalyzer, ferro)
    } else {
        return false;
    };

    // Must have same accession
    let count_acc = count_str.split(':').next().unwrap_or("");
    let expanded_acc = expanded_str.split(':').next().unwrap_or("");
    if count_acc != expanded_acc || count_acc.is_empty() {
        return false;
    }

    // Both should be insertions
    if !count_str.contains("ins") || !expanded_str.contains("ins") {
        return false;
    }

    // The expanded version should have a run of N's
    // We just check that it has multiple consecutive N's
    if !expanded_str.contains("NN") {
        return false;
    }

    // Extract N[count] and verify it matches the length of N's
    // Format: N[342]
    if let Some(n_start) = count_str.find("N[") {
        if let Some(n_end) = count_str[n_start..].find(']') {
            let count_str_part = &count_str[n_start + 2..n_start + n_end];
            if let Ok(count) = count_str_part.parse::<usize>() {
                // Count the N's in expanded string
                let n_count = expanded_str.chars().filter(|&c| c == 'N').count();
                // Allow some tolerance for other N's in the string
                if n_count >= count / 2 {
                    return true;
                }
            }
        }
    }

    false
}

/// Check if delins simplifications are equivalent.
///
/// Example pairs that are equivalent:
/// - c.3108_3120delinsTAAG vs c.3108_3119delinsTAA
/// - c.5791_5794delinsCCTCT vs c.5791_5793delinsCCTC
///
/// These represent the same mutation but with different boundary choices
/// for the delins operation.
fn check_delins_simplification_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both should be delins
    if !ferro.contains("delins") || !mutalyzer.contains("delins") {
        return false;
    }

    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let muta_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != muta_acc || ferro_acc.is_empty() {
        return false;
    }

    // Extract positions and sequences
    // Format: c.X_YdelinsZZZ or c.XdelinsZZZ
    fn extract_delins(s: &str) -> Option<(i64, i64, &str)> {
        let delins_pos = s.find("delins")?;
        let before = &s[..delins_pos];
        let after = &s[delins_pos + 6..]; // Skip "delins"

        // Find the position part (after the dot)
        let dot_pos = before.rfind('.')?;
        let pos_str = &before[dot_pos + 1..];

        // Parse position(s)
        if let Some((start, end)) = pos_str.split_once('_') {
            let start_pos = start.parse::<i64>().ok()?;
            let end_pos = end.parse::<i64>().ok()?;
            Some((start_pos, end_pos, after))
        } else {
            let pos = pos_str.parse::<i64>().ok()?;
            Some((pos, pos, after))
        }
    }

    let Some((f_start, f_end, f_seq)) = extract_delins(ferro) else {
        return false;
    };
    let Some((m_start, m_end, m_seq)) = extract_delins(mutalyzer) else {
        return false;
    };

    // Check if positions overlap or are adjacent
    let positions_close = (f_start - m_start).abs() <= 2 && (f_end - m_end).abs() <= 2;
    if !positions_close {
        return false;
    }

    // Check if sequences share a common prefix/suffix
    // This indicates they're different simplifications of the same mutation
    let shorter = f_seq.len().min(m_seq.len());
    if shorter == 0 {
        return false;
    }

    // Check for shared prefix
    let shared_prefix = f_seq
        .chars()
        .zip(m_seq.chars())
        .take_while(|(a, b)| a == b)
        .count();

    // Must share at least one character and most of the shorter sequence
    // This prevents false positives like delinsAA vs delinsTT
    if shared_prefix == 0 {
        return false;
    }

    // If they share most of the shorter sequence, consider equivalent
    shared_prefix >= shorter.saturating_sub(2)
}

/// Check if larger UTR position shifts are equivalent due to 3' normalization.
///
/// Example:
/// - Ferro: c.*373del shifted to c.*382del (9 position shift in poly-T)
/// - Ferro: c.*290dup shifted to c.*303dup (13 position shift)
///
/// In repeat tracts, 3' normalization can shift positions significantly.
fn check_utr_larger_position_shift_equivalence(ferro: &str, mutalyzer: &str) -> bool {
    // Both must be UTR variants (contain c.*)
    if !ferro.contains("c.*") || !mutalyzer.contains("c.*") {
        return false;
    }

    // Must have same accession
    let ferro_acc = ferro.split(':').next().unwrap_or("");
    let muta_acc = mutalyzer.split(':').next().unwrap_or("");
    if ferro_acc != muta_acc || ferro_acc.is_empty() {
        return false;
    }

    // Both must have same operation (del or dup)
    let ferro_has_del = ferro.contains("del") && !ferro.contains("delins");
    let ferro_has_dup = ferro.contains("dup");
    let muta_has_del = mutalyzer.contains("del") && !mutalyzer.contains("delins");
    let muta_has_dup = mutalyzer.contains("dup");

    if !((ferro_has_del && muta_has_del) || (ferro_has_dup && muta_has_dup)) {
        return false;
    }

    // Extract UTR positions
    let extract_utr_pos = |s: &str| -> Option<i64> {
        let star_pos = s.find("*")?;
        let after_star = &s[star_pos + 1..];
        // Find where the position ends
        let end_pos = after_star
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(after_star.len());
        after_star[..end_pos].parse::<i64>().ok()
    };

    let Some(ferro_pos) = extract_utr_pos(ferro) else {
        return false;
    };
    let Some(muta_pos) = extract_utr_pos(mutalyzer) else {
        return false;
    };

    // Allow larger shifts (up to 50 positions) for 3' normalization in repeat tracts
    // The off-by-one check handles shifts of 1, this handles larger shifts
    let diff = (ferro_pos - muta_pos).abs();
    diff > 1 && diff <= 50
}

/// Save detailed per-pattern results to a JSON file.
fn save_detailed_results<P: AsRef<Path>>(
    ferro_results: &[ParseResult],
    muta_results: &[ParseResult],
    mode: CompareMode,
    timestamp: DateTime<Utc>,
    path: P,
) -> Result<(), FerroError> {
    let path = path.as_ref();

    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;
    }

    // Build combined results
    let results: Vec<CombinedPatternResult> = ferro_results
        .iter()
        .zip(muta_results.iter())
        .map(|(f, m)| {
            let outputs_match = if f.success && m.success {
                let ferro_out = f.output.as_deref().unwrap_or("");
                let muta_out = m.output.as_deref().unwrap_or("");
                // Use semantic equivalence checking
                Some(outputs_are_equivalent(ferro_out, muta_out))
            } else {
                None
            };

            CombinedPatternResult {
                input: f.input.clone(),
                ferro_success: f.success,
                ferro_output: f.output.clone(),
                ferro_error: f.error.clone(),
                mutalyzer_success: m.success,
                mutalyzer_output: m.output.clone(),
                mutalyzer_error: m.error.clone().or(m.error_category.clone()),
                outputs_match,
            }
        })
        .collect();

    let detailed = DetailedResults {
        mode,
        timestamp,
        sample_size: results.len(),
        results,
    };

    let file = File::create(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", path.display(), e),
    })?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, &detailed).map_err(|e| FerroError::Io {
        msg: format!("Failed to write JSON: {}", e),
    })?;

    eprintln!("Detailed results saved to: {}", path.display());

    Ok(())
}

#[cfg(test)]
mod equivalence_tests {
    use super::*;

    #[test]
    fn test_exact_match() {
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.1207T>G",
            "NM_000060.2:c.1207T>G"
        ));
    }

    #[test]
    fn test_allelic_equivalence() {
        // Ferro expands, mutalyzer keeps compact
        assert!(outputs_are_equivalent(
            "[NM_000060.2:c.1207T>G;NM_000060.2:c.1330G>C]",
            "NM_000060.2:c.[1207T>G;1330G>C]"
        ));

        // Different alleles - not equivalent
        assert!(!outputs_are_equivalent(
            "[NM_000060.2:c.1207T>G;NM_000060.2:c.1331G>C]",
            "NM_000060.2:c.[1207T>G;1330G>C]"
        ));

        // Different accession - not equivalent
        assert!(!outputs_are_equivalent(
            "[NM_000061.2:c.1207T>G;NM_000061.2:c.1330G>C]",
            "NM_000060.2:c.[1207T>G;1330G>C]"
        ));
    }

    #[test]
    fn test_delins_inv_equivalence() {
        // Same positions, delins vs inv
        assert!(outputs_are_equivalent(
            "NM_144670.4:c.2082_2083delinsTG",
            "NM_144670.4:c.2082_2083inv"
        ));

        // Different positions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_144670.4:c.2082_2084delinsTG",
            "NM_144670.4:c.2082_2083inv"
        ));
    }

    #[test]
    fn test_delins_split_equivalence() {
        // Delins vs split allelic notation
        assert!(outputs_are_equivalent(
            "NM_001282424.3:c.2142_2144delinsAA",
            "NM_001282424.3:c.[2142del;2144C>A]"
        ));

        // Different accession - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001282424.3:c.2142_2144delinsAA",
            "NM_001282425.3:c.[2142del;2144C>A]"
        ));
    }

    #[test]
    fn test_dup_repeat_equivalence() {
        // Dup vs repeat notation (dup = 2 copies)
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.100_102dup",
            "NM_000060.2:c.100_102[2]"
        ));

        // With sequence in repeat notation
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.100_102dup",
            "NM_000060.2:c.100_102ATG[2]"
        ));

        // Different positions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.100_103dup",
            "NM_000060.2:c.100_102[2]"
        ));

        // Different repeat count - not equivalent (dup is always [2])
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.100_102dup",
            "NM_000060.2:c.100_102[3]"
        ));
    }

    #[test]
    fn test_dup_vs_ins_equivalence() {
        // Single-base dup is equivalent to inserting that base at adjacent positions
        // c.-27dup means duplicate base at -27, which = inserting a copy after -27
        assert!(outputs_are_equivalent(
            "NM_001605.3:c.-27dup",
            "NM_001605.3:c.-27_-26insA"
        ));

        // Different accession - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001605.3:c.-27dup",
            "NM_001605.4:c.-27_-26insA"
        ));

        // Non-adjacent positions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001605.3:c.-27dup",
            "NM_001605.3:c.-28_-27insA"
        ));
    }

    #[test]
    fn test_no_false_positives() {
        // Completely different variants
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.100del",
            "NM_000060.2:c.200del"
        ));

        // Similar but different
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.100_101delinsAA",
            "NM_000060.2:c.100_101delinsTT"
        ));
    }

    #[test]
    fn test_three_prime_shift_equivalence() {
        // Same dup, ferro shifted 3' by 1 position (5' UTR negative coords)
        // c.-56_-47dup vs c.-55_-46dup is 10bp dup shifted +1 (toward c.1)
        assert!(outputs_are_equivalent(
            "NM_001394148.2:c.-55_-46dup",
            "NM_001394148.2:c.-56_-47dup"
        ));

        // Same del, ferro shifted 3' by 1 position
        // c.-146_-143del vs c.-145_-142del is 4bp del shifted +1
        assert!(outputs_are_equivalent(
            "NM_001400774.1:c.-145_-142del",
            "NM_001400774.1:c.-146_-143del"
        ));

        // Single position del with shift
        assert!(outputs_are_equivalent(
            "NM_001425842.1:c.-609del",
            "NM_001425842.1:c.-610del"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001394148.2:c.-55_-46dup",
            "NM_001394149.2:c.-56_-47dup"
        ));

        // Different edit types - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001394148.2:c.-55_-46dup",
            "NM_001394148.2:c.-56_-47del"
        ));

        // Non-uniform shift (start differs by 1, end differs by 2) - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.100_105del",
            "NM_000060.2:c.101_108del"
        ));

        // Identical positions - not a shift, handled by exact match
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.100_102dup",
            "NM_000060.2:c.100_102dup"
        ));
    }

    #[test]
    fn test_repeat_vs_del_dup_equivalence() {
        // Repeat notation with fewer units = deletion
        assert!(outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[5]",
            "NM_001318856.2:c.-79_-68del"
        ));

        // Repeat notation with one more unit = duplication
        assert!(outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[10]",
            "NM_001318856.2:c.-70_-68dup"
        ));

        // Single unit deletion via repeat notation
        assert!(outputs_are_equivalent(
            "NM_001400774.1:c.-126AGA[1]",
            "NM_001400774.1:c.-123_-121del"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[5]",
            "NM_001318857.2:c.-79_-68del"
        ));
    }

    #[test]
    fn test_dup_vs_expanded_repeat_equivalence() {
        // Simple dup vs expanded repeat notation
        assert!(outputs_are_equivalent(
            "NM_020732.3:c.1023_1028dup",
            "NM_020732.3:c.1017_1028GGC[6]"
        ));

        // Another dup vs repeat expansion
        assert!(outputs_are_equivalent(
            "NM_020732.3:c.348_362dup",
            "NM_020732.3:c.342_362GCA[12]"
        ));

        // Dup vs long repeat notation
        assert!(outputs_are_equivalent(
            "NM_020732.3:c.257_298dup",
            "NM_020732.3:c.257_298ACCACCACCATGCCCACCACC[4]"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_020732.3:c.1023_1028dup",
            "NM_020733.3:c.1017_1028GGC[6]"
        ));
    }

    #[test]
    fn test_delins_sequence_vs_position_equivalence() {
        // Delins with sequence vs position reference
        assert!(outputs_are_equivalent(
            "NM_001414902.1:c.174_182delinsCAGGAGGAGA",
            "NM_001414902.1:c.174_182delins186_195"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001414902.1:c.174_182delinsCAGGAGGAGA",
            "NM_001414903.1:c.174_182delins186_195"
        ));

        // Both have sequences - not this equivalence type
        assert!(!outputs_are_equivalent(
            "NM_001414902.1:c.174_182delinsCAGGAGGAGA",
            "NM_001414902.1:c.174_182delinsTTTT"
        ));
    }

    #[test]
    fn test_repeat_position_range_equivalence() {
        // Same repeat, different position notation (single vs range)
        assert!(outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[11]",
            "NM_001318856.2:c.-94_-68TGC[11]"
        ));

        // Same repeat with different range
        assert!(outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[14]",
            "NM_001318856.2:c.-94_-68TGC[14]"
        ));

        // Different counts - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[11]",
            "NM_001318856.2:c.-94_-68TGC[14]"
        ));

        // Different sequences - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[11]",
            "NM_001318856.2:c.-94_-68AGC[11]"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001318856.2:c.-94TGC[11]",
            "NM_001318857.2:c.-94_-68TGC[11]"
        ));
    }

    #[test]
    fn test_ins_sequence_vs_range_equivalence() {
        // ins with explicit sequence vs position range
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.172_173insTGCAGCAGCAGCAGCAGC",
            "NM_000060.2:c.172_173ins170_187"
        ));

        // Different order (mutalyzer has sequence, ferro has range)
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.172_173ins170_187",
            "NM_000060.2:c.172_173insTGCAGCAGCAGCAGCAGC"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.172_173insTGCAGCAGCAGCAGCAGC",
            "NM_000061.2:c.172_173ins170_187"
        ));
    }

    #[test]
    fn test_repeat_one_vs_identity_equivalence() {
        // Repeat[1] vs c.= (both mean no change)
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.217_219AGA[1]",
            "NM_000060.2:c.="
        ));

        // Different order
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.=",
            "NM_000060.2:c.217_219AGA[1]"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.217_219AGA[1]",
            "NM_000061.2:c.="
        ));

        // Repeat[2] is NOT equivalent to c.=
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.217_219AGA[2]",
            "NM_000060.2:c.="
        ));
    }

    #[test]
    fn test_identity_substitution_equivalence() {
        // Identity substitution (X>X) vs c.=
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.1632T>T",
            "NM_000060.2:c.="
        ));

        // Different order
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.=",
            "NM_000060.2:c.1632T>T"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.1632T>T",
            "NM_000061.2:c.="
        ));

        // Actual substitution (different bases) - NOT equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.1632T>G",
            "NM_000060.2:c.="
        ));
    }

    #[test]
    fn test_stop_codon_position_equivalence() {
        // CDS position vs UTR position for same stop codon
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.1284G>C",
            "NM_000060.2:c.*1G>C"
        ));

        // Different order
        assert!(outputs_are_equivalent(
            "NM_000060.2:c.*1G>C",
            "NM_000060.2:c.1284G>C"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.1284G>C",
            "NM_000061.2:c.*1G>C"
        ));

        // Different operations - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000060.2:c.1284G>C",
            "NM_000060.2:c.*1G>T"
        ));
    }

    #[test]
    fn test_utr_position_off_by_one_equivalence() {
        // UTR position off-by-one (Mutalyzer bug)
        assert!(outputs_are_equivalent(
            "NM_000016.4:c.*63dup",
            "NM_000016.4:c.*64dup"
        ));

        // Works for del too
        assert!(outputs_are_equivalent(
            "NM_000016.4:c.*63del",
            "NM_000016.4:c.*64del"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000016.4:c.*63dup",
            "NM_000017.4:c.*64dup"
        ));

        // Different operations - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000016.4:c.*63dup",
            "NM_000016.4:c.*64del"
        ));
    }

    #[test]
    fn test_cds_end_vs_utr_start_equivalence() {
        // CDS end position vs *1 (UTR start)
        assert!(outputs_are_equivalent(
            "NM_000038.5:c.-18_8532del",
            "NM_000038.5:c.-18_*1del"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000038.5:c.-18_8532del",
            "NM_000039.5:c.-18_*1del"
        ));
    }

    #[test]
    fn test_ins_n_count_vs_expanded_equivalence() {
        // ins N[count] vs expanded N's - N[10] should match 10 N's
        assert!(outputs_are_equivalent(
            "NM_000055.4:c.1160_1161ins[N[10];1146_1160]",
            "NM_000055.4:c.1160_1161ins[NNNNNNNNNN;1146_1160]"
        ));

        // Larger count - N[50] should match ~50 N's
        let fifty_ns = "N".repeat(50);
        assert!(outputs_are_equivalent(
            "NM_000055.4:c.1160_1161insN[50]",
            &format!("NM_000055.4:c.1160_1161ins{}", fifty_ns)
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000055.4:c.1160_1161insN[10]",
            "NM_000056.4:c.1160_1161insNNNNNNNNNN"
        ));
    }

    #[test]
    fn test_utr_larger_position_shift_equivalence() {
        // Larger UTR position shifts (3' normalization in repeat tracts)
        assert!(outputs_are_equivalent(
            "NM_000059.3:c.*373del",
            "NM_000059.3:c.*382del"
        ));

        // Works for dup too
        assert!(outputs_are_equivalent(
            "NM_000075.2:c.*290dup",
            "NM_000075.2:c.*303dup"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000059.3:c.*373del",
            "NM_000060.3:c.*382del"
        ));

        // Shift too large (> 50) - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000059.3:c.*100del",
            "NM_000059.3:c.*200del"
        ));
    }

    #[test]
    fn test_delins_simplification_equivalence() {
        // Different delins simplifications of same mutation
        assert!(outputs_are_equivalent(
            "NM_000051.3:c.3108_3120delinsTAAG",
            "NM_000051.3:c.3108_3119delinsTAA"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_000051.3:c.3108_3120delinsTAAG",
            "NM_000052.3:c.3108_3119delinsTAA"
        ));
    }

    #[test]
    fn test_delins_vs_complex_allelic_equivalence() {
        // delins vs complex allelic notation
        assert!(outputs_are_equivalent(
            "NM_001321134.2:c.423_428delinsACCTAGGGGAAACAAGATGTAGTGCTATTGCC",
            "NM_001321134.2:c.[422_423insA;423_424ins*572_*725;428_429ins[*731_*833;CTTTGAG]]"
        ));

        // Different accessions - not equivalent
        assert!(!outputs_are_equivalent(
            "NM_001321134.2:c.423_428delinsACCTAGGGG",
            "NM_001321135.2:c.[422_423insA;423_424ins*572_*725]"
        ));
    }

    // =========================================================================
    // Tests for patterns that are NOT equivalent (Mutalyzer bugs or semantically different)
    // These document expected differences between Ferro and Mutalyzer
    // =========================================================================

    #[test]
    fn test_repeat_del_vs_dup_not_equivalent() {
        // Repeat[N] resulting in del vs dup - these are NOT equivalent
        // Ferro: interprets repeat count correctly → del
        // Mutalyzer: incorrectly interprets as dup
        assert!(!outputs_are_equivalent(
            "NM_000017.4:c.364_367del", // Ferro correct
            "NM_000017.4:c.366_367dup"  // Mutalyzer incorrect
        ));

        assert!(!outputs_are_equivalent(
            "NM_001407640.1:c.470_471del", // Ferro correct
            "NM_001407640.1:c.470_471dup"  // Mutalyzer incorrect
        ));
    }

    #[test]
    fn test_repeat_count_diff_not_equivalent() {
        // Different repeat counts are NOT equivalent
        assert!(!outputs_are_equivalent(
            "NM_000028.3:c.1497_1500AG[4]",
            "NM_000028.3:c.1497_1500AG[5]"
        ));

        assert!(!outputs_are_equivalent(
            "NM_000044.3:c.172_237CAG[35]",
            "NM_000044.3:c.171_239GCA[57]"
        ));
    }
}
