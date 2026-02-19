//! Mutalyzer cache population from local reference data.
//!
//! This module creates mutalyzer-compatible cache files from our downloaded
//! FASTA and cdot data, enabling fair local-only comparisons.
//!
//! # Coordinate System
//!
//! This module handles mixed coordinate systems from cdot:
//!
//! | Field | Basis | Notes |
//! |-------|-------|-------|
//! | cdot genomic (genome_start/genome_end) | 0-based | Half-open `[start, end)` |
//! | cdot transcript (tx_start/tx_end) | 1-based | Inclusive |
//! | cdot CDS (cds_start/cds_end) | 0-based | Half-open |
//! | Mutalyzer cache coordinates | Mixed | Follows mutalyzer conventions |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

#![allow(clippy::type_complexity)]

use crate::benchmark::translate::derive_protein_from_transcript;
use crate::coords::hgvs_pos_to_index;
use crate::data::cdot::{CdotMapper, CdotTranscript};
use crate::FerroError;
use rayon::prelude::*;
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Get NCBI API key from environment variable.
///
/// When set, allows 10 requests/second instead of 3 requests/second.
/// Get a key at: https://www.ncbi.nlm.nih.gov/account/settings/
fn get_ncbi_api_key() -> Option<String> {
    std::env::var("NCBI_API_KEY").ok().filter(|k| !k.is_empty())
}

/// Build efetch URL with optional API key.
fn build_efetch_url(db: &str, ids: &str, rettype: &str, retmode: &str) -> String {
    let base = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&id={}&rettype={}&retmode={}",
        db, ids, rettype, retmode
    );
    if let Some(key) = get_ncbi_api_key() {
        format!("{}&api_key={}", base, key)
    } else {
        base
    }
}

/// Build efetch URL using WebEnv history server (for large batch fetches).
#[allow(dead_code)]
fn build_efetch_url_with_history(
    db: &str,
    web_env: &str,
    query_key: &str,
    retstart: usize,
    retmax: usize,
    rettype: &str,
    retmode: &str,
) -> String {
    let base = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={}&WebEnv={}&query_key={}&retstart={}&retmax={}&rettype={}&retmode={}",
        db, web_env, query_key, retstart, retmax, rettype, retmode
    );
    if let Some(key) = get_ncbi_api_key() {
        format!("{}&api_key={}", base, key)
    } else {
        base
    }
}

/// Post IDs to NCBI history server, returns (WebEnv, query_key).
/// This allows efficient batch fetching without resending IDs each time.
#[allow(dead_code)]
fn epost_ids(db: &str, ids: &[&String]) -> Result<(String, String), FerroError> {
    use std::process::Command;

    let url = if let Some(key) = get_ncbi_api_key() {
        format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db={}&api_key={}",
            db, key
        )
    } else {
        format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi?db={}",
            db
        )
    };

    let id_list = ids.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(",");

    let output = Command::new("curl")
        .args(["-s", "-X", "POST", "-d", &format!("id={}", id_list), &url])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to POST to NCBI epost: {}", e),
        })?;

    if !output.status.success() {
        return Err(FerroError::Io {
            msg: "NCBI epost request failed".to_string(),
        });
    }

    let response = String::from_utf8_lossy(&output.stdout);

    // Parse XML response for WebEnv and QueryKey
    let web_env = response
        .split("<WebEnv>")
        .nth(1)
        .and_then(|s| s.split("</WebEnv>").next())
        .ok_or_else(|| FerroError::Io {
            msg: "Failed to parse WebEnv from epost response".to_string(),
        })?
        .to_string();

    let query_key = response
        .split("<QueryKey>")
        .nth(1)
        .and_then(|s| s.split("</QueryKey>").next())
        .ok_or_else(|| FerroError::Io {
            msg: "Failed to parse QueryKey from epost response".to_string(),
        })?
        .to_string();

    Ok((web_env, query_key))
}

/// Populate mutalyzer cache from local reference data.
///
/// Reads FASTA files and cdot transcript metadata to create mutalyzer-compatible
/// cache files (.sequence and .annotations).
///
/// If `filter_file` is provided, only cache entries for accessions listed in that file.
pub fn populate_mutalyzer_cache<P: AsRef<Path>>(
    reference_dir: P,
    cache_dir: P,
    filter_file: Option<P>,
) -> Result<CacheStats, FerroError> {
    let reference_dir = reference_dir.as_ref();
    let cache_dir = cache_dir.as_ref();

    // Create cache directory
    std::fs::create_dir_all(cache_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create cache dir {}: {}", cache_dir.display(), e),
    })?;

    // Load filter list if provided
    let filter_set: Option<HashSet<String>> = if let Some(filter_path) = filter_file {
        let filter_path = filter_path.as_ref();
        eprintln!("Loading filter list from {}...", filter_path.display());
        let file = File::open(filter_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open filter file: {}", e),
        })?;
        let reader = BufReader::new(file);
        let set: HashSet<String> = reader
            .lines()
            .map_while(Result::ok)
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect();
        eprintln!("  Loaded {} accessions to filter", set.len());
        Some(set)
    } else {
        None
    };

    // Load manifest to find cdot file
    let manifest_path = reference_dir.join("manifest.json");
    if !manifest_path.exists() {
        return Err(FerroError::Io {
            msg: format!(
                "Manifest not found at {}. Run 'ferro-benchmark prepare' first.",
                manifest_path.display()
            ),
        });
    }

    let manifest: serde_json::Value = {
        let file = File::open(&manifest_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open manifest: {}", e),
        })?;
        serde_json::from_reader(file).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse manifest: {}", e),
        })?
    };

    // Load cdot
    let cdot_relative = manifest
        .get("cdot_json")
        .and_then(|v| v.as_str())
        .ok_or_else(|| FerroError::Io {
            msg: "No cdot_json in manifest. Run 'ferro-benchmark prepare' with cdot.".to_string(),
        })?;
    let cdot_path = reference_dir.join(cdot_relative);

    eprintln!("Loading cdot from {}...", cdot_path.display());
    let mut cdot = if cdot_relative.ends_with(".gz") {
        CdotMapper::from_json_gz(&cdot_path)?
    } else {
        CdotMapper::from_json_file(&cdot_path)?
    };
    let transcript_count = cdot.transcript_count();
    eprintln!("  Loaded {} transcripts", transcript_count);

    // Load LRG to RefSeq mapping if available
    if let Some(lrg_mapping_relative) = manifest.get("lrg_refseq_mapping").and_then(|v| v.as_str())
    {
        let path = reference_dir.join(lrg_mapping_relative);
        if path.exists() {
            let count = cdot.load_lrg_mapping(&path)?;
            eprintln!("  Loaded {} LRG to RefSeq mappings", count);
        }
    }

    // Load transcript sequences from FASTA files
    eprintln!("Loading transcript sequences from FASTA files...");
    let mut sequences = load_fasta_sequences_parallel(reference_dir)?;
    eprintln!("  Loaded {} transcript sequences", sequences.len());

    // Load supplemental transcripts from ferro reference (patterns_transcripts.fna)
    let supplemental_dir = reference_dir.join("supplemental");
    let mut supplemental_cds: HashMap<String, (Option<u64>, Option<u64>, Option<String>)> =
        HashMap::new();
    if supplemental_dir.exists() {
        let supplemental = load_fasta_sequences_parallel(&supplemental_dir)?;
        let count = supplemental.len();
        for (id, seq) in supplemental {
            sequences.push((id, seq));
        }
        if count > 0 {
            eprintln!("  Loaded {} supplemental transcript sequences", count);
        }

        // Load supplemental metadata for CDS coordinates
        let metadata_path = supplemental_dir.join("patterns_transcripts.metadata.json");
        if metadata_path.exists() {
            let file = File::open(&metadata_path).map_err(|e| FerroError::Io {
                msg: format!("Failed to open {}: {}", metadata_path.display(), e),
            })?;
            let metadata: SupplementalMetadata =
                serde_json::from_reader(file).map_err(|e| FerroError::Json {
                    msg: format!("Failed to parse {}: {}", metadata_path.display(), e),
                })?;
            for (acc, info) in metadata.transcripts {
                supplemental_cds.insert(acc, (info.cds_start, info.cds_end, info.gene_symbol));
            }
            eprintln!(
                "  Loaded CDS metadata for {} supplemental transcripts",
                supplemental_cds.len()
            );
        }
    }

    // Load genomic sequences from genome subdirectory (NC_*)
    let genome_dir = reference_dir.join("genome");
    if genome_dir.exists() {
        eprintln!("Loading genomic sequences from {}...", genome_dir.display());
        let genomic_seqs = load_genomic_sequences(&genome_dir)?;
        eprintln!("  Loaded {} genomic sequences", genomic_seqs.len());
        sequences.extend(genomic_seqs);
    }

    // Load RefSeqGene sequences from refseqgene subdirectory (NG_*)
    let refseqgene_dir = reference_dir.join("refseqgene");
    if refseqgene_dir.exists() {
        eprintln!(
            "Loading RefSeqGene sequences from {}...",
            refseqgene_dir.display()
        );
        let refseqgene_seqs = load_refseqgene_sequences(&refseqgene_dir)?;
        eprintln!("  Loaded {} RefSeqGene sequences", refseqgene_seqs.len());
        sequences.extend(refseqgene_seqs);
    }

    // Fetch legacy GenBank sequences (directly to cache dir, not loaded into memory)
    // Prefers ferro reference over NCBI when legacy_genbank.fna is available
    eprintln!("Fetching legacy GenBank sequences...");
    let legacy_accessions: Vec<String> = LEGACY_GENBANK_ACCESSIONS
        .iter()
        .map(|s| s.to_string())
        .collect();
    let legacy_fetched =
        fetch_legacy_genbank_sequences(&legacy_accessions, cache_dir, Some(reference_dir))?;
    eprintln!("  Fetched {} legacy GenBank sequences", legacy_fetched);

    // Fetch legacy transcript versions (older RefSeq versions not in cdot)
    // Prefers ferro reference over NCBI when legacy transcripts are available
    eprintln!("Fetching legacy transcript versions...");
    let legacy_transcripts_fetched =
        fetch_legacy_transcript_versions(cache_dir, Some(reference_dir))?;
    eprintln!(
        "  Fetched {} legacy transcript versions",
        legacy_transcripts_fetched
    );

    // Filter sequences if filter set provided
    let sequences_to_process: Vec<_> = if let Some(ref filter) = filter_set {
        sequences
            .into_iter()
            .filter(|(id, _)| filter.contains(id))
            .collect()
    } else {
        sequences.into_iter().collect()
    };

    let total_to_process = sequences_to_process.len();
    eprintln!(
        "Creating {} mutalyzer cache entries ({} threads)...",
        total_to_process,
        rayon::current_num_threads()
    );

    // Progress counter
    let created = AtomicUsize::new(0);
    let errors = AtomicUsize::new(0);

    // Create cache entries in parallel
    sequences_to_process.par_iter().for_each(|(tx_id, seq)| {
        let result = create_cache_entry(cache_dir, tx_id, seq, cdot.get_transcript(tx_id));

        if result.is_ok() {
            let count = created.fetch_add(1, Ordering::Relaxed) + 1;
            if count.is_multiple_of(50000) {
                eprintln!("  Created {} cache entries...", count);
            }
        } else {
            errors.fetch_add(1, Ordering::Relaxed);
        }
    });

    let cache_entries_created = created.load(Ordering::Relaxed);
    let error_count = errors.load(Ordering::Relaxed);

    eprintln!(
        "Created {} cache entries in {} ({} errors)",
        cache_entries_created,
        cache_dir.display(),
        error_count
    );

    // Derive protein sequences from coding transcripts
    eprintln!("Deriving protein sequences from transcripts...");
    let mut proteins_derived = 0usize;
    let mut proteins_skipped = 0usize;

    for (tx_id, seq) in &sequences_to_process {
        // Skip non-coding transcripts
        if !tx_id.starts_with("NM_") && !tx_id.starts_with("XM_") {
            continue;
        }

        // Get CDS coordinates from cdot or supplemental metadata
        let (cds_start, cds_end, protein_id) = if let Some(tx) = cdot.get_transcript(tx_id) {
            // cdot CDS coordinates are 0-based, half-open
            let start = tx.cds_start.unwrap_or(0) as usize;
            let end = tx.cds_end.unwrap_or(seq.len() as u64) as usize;
            let prot_id = tx
                .protein
                .clone()
                .unwrap_or_else(|| tx_id.replace("NM_", "NP_").replace("XM_", "XP_"));
            (start, end, prot_id)
        } else if let Some((start, end, _gene)) = supplemental_cds.get(tx_id) {
            // Supplemental metadata CDS coordinates are 1-based, inclusive
            // Convert to 0-based, half-open
            let start = start.map(|s| (s - 1) as usize).unwrap_or(0);
            let end = end.map(|e| e as usize).unwrap_or(seq.len());
            let prot_id = tx_id.replace("NM_", "NP_").replace("XM_", "XP_");
            (start, end, prot_id)
        } else {
            continue; // No CDS info available
        };

        // Check if protein cache entry already exists
        let prot_seq_path = cache_dir.join(format!("{}.sequence", protein_id));
        if prot_seq_path.exists() {
            proteins_skipped += 1;
            continue;
        }

        // Derive protein sequence
        if let Some(protein_seq) = derive_protein_from_transcript(seq, cds_start, cds_end) {
            // Write protein sequence file
            if let Err(e) = std::fs::write(&prot_seq_path, &protein_seq) {
                eprintln!(
                    "  Warning: Failed to write {}: {}",
                    prot_seq_path.display(),
                    e
                );
                continue;
            }

            // Write protein annotations file
            let annotations = build_protein_annotations(&protein_id, protein_seq.len());
            let ann_path = cache_dir.join(format!("{}.annotations", protein_id));
            if let Err(e) = std::fs::write(
                &ann_path,
                serde_json::to_string(&annotations).unwrap_or_default(),
            ) {
                eprintln!("  Warning: Failed to write {}: {}", ann_path.display(), e);
                continue;
            }

            proteins_derived += 1;
        }
    }

    eprintln!(
        "  Derived {} protein sequences ({} already cached)",
        proteins_derived, proteins_skipped
    );

    // Write settings file
    let settings_path = cache_dir.join("mutalyzer_settings.conf");
    let settings_content = format!(
        "MUTALYZER_CACHE_DIR = {}\nMUTALYZER_FILE_CACHE_ADD = false\n",
        cache_dir.display()
    );
    std::fs::write(&settings_path, settings_content).map_err(|e| FerroError::Io {
        msg: format!("Failed to write settings: {}", e),
    })?;
    eprintln!("Settings file: {}", settings_path.display());

    Ok(CacheStats {
        transcripts_processed: total_to_process,
        sequences_matched: total_to_process,
        cache_entries_created,
        legacy_sequences_fetched: legacy_fetched + legacy_transcripts_fetched,
    })
}

/// Create a single cache entry (sequence + annotations files).
fn create_cache_entry(
    cache_dir: &Path,
    tx_id: &str,
    seq: &str,
    tx_data: Option<&CdotTranscript>,
) -> Result<(), FerroError> {
    // Skip if already exists (makes re-runs fast)
    let seq_path = cache_dir.join(format!("{}.sequence", tx_id));
    let ann_path = cache_dir.join(format!("{}.annotations", tx_id));
    if seq_path.exists() && ann_path.exists() {
        return Ok(());
    }

    // Create sequence file
    std::fs::write(&seq_path, seq).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", seq_path.display(), e),
    })?;

    // Create annotations file (compact JSON, not pretty-printed)
    let annotations = if let Some(tx) = tx_data {
        build_annotations_from_cdot(tx_id, seq.len(), tx)
    } else {
        build_minimal_annotations(tx_id, seq.len())
    };

    let ann_file = File::create(&ann_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", ann_path.display(), e),
    })?;
    serde_json::to_writer(ann_file, &annotations).map_err(|e| FerroError::Io {
        msg: format!("Failed to write annotations: {}", e),
    })?;

    Ok(())
}

/// Statistics from cache population.
#[derive(Debug, Clone)]
pub struct CacheStats {
    pub transcripts_processed: usize,
    pub sequences_matched: usize,
    pub cache_entries_created: usize,
    pub legacy_sequences_fetched: usize,
}

/// Load sequences from a single FASTA file into a HashMap.
///
/// Returns a map from sequence ID (first word of header) to sequence string.
fn load_fasta_sequences(fasta_path: &Path) -> Result<HashMap<String, String>, FerroError> {
    let file = File::open(fasta_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", fasta_path.display(), e),
    })?;
    let reader = BufReader::new(file);

    let mut sequences = HashMap::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if let Some(header) = line.strip_prefix('>') {
            // Save previous sequence
            if let Some(id) = current_id.take() {
                if !current_seq.is_empty() {
                    sequences.insert(id, std::mem::take(&mut current_seq));
                }
            }
            // Parse new ID (first word of header)
            current_id = header.split_whitespace().next().map(|s| s.to_string());
        } else if current_id.is_some() {
            current_seq.push_str(line.trim());
        }
    }

    // Save final sequence
    if let Some(id) = current_id {
        if !current_seq.is_empty() {
            sequences.insert(id, current_seq);
        }
    }

    Ok(sequences)
}

/// Load all sequences from FASTA files in a directory (parallel file reading).
fn load_fasta_sequences_parallel<P: AsRef<Path>>(
    dir: P,
) -> Result<Vec<(String, String)>, FerroError> {
    let dir = dir.as_ref();
    let mut sequences = HashMap::new();

    // Find all .fna files
    let fasta_files: Vec<_> = std::fs::read_dir(dir)
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to read {}: {}", dir.display(), e),
        })?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .is_some_and(|e| e == "fna" || e == "fa" || e == "fasta")
        })
        .collect();

    // Also check subdirectories
    let mut all_fasta = fasta_files;
    for entry in std::fs::read_dir(dir).into_iter().flatten().flatten() {
        let path = entry.path();
        if path.is_dir() {
            if let Ok(subdir) = std::fs::read_dir(&path) {
                for subentry in subdir.flatten() {
                    let subpath = subentry.path();
                    if subpath
                        .extension()
                        .is_some_and(|e| e == "fna" || e == "fa" || e == "fasta")
                    {
                        all_fasta.push(subpath);
                    }
                }
            }
        }
    }

    for fasta_path in all_fasta {
        let file = File::open(&fasta_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", fasta_path.display(), e),
        })?;
        let reader = BufReader::new(file);

        let mut current_id: Option<String> = None;
        let mut current_seq = String::new();

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read line: {}", e),
            })?;

            if let Some(header) = line.strip_prefix('>') {
                // Save previous sequence
                if let Some(id) = current_id.take() {
                    if !current_seq.is_empty() {
                        sequences.insert(id, std::mem::take(&mut current_seq));
                    }
                }
                // Parse new ID (format: ">NM_000088.4 ..." or ">NM_000088.4")
                current_id = header.split_whitespace().next().map(|s| s.to_string());
            } else if current_id.is_some() {
                current_seq.push_str(line.trim());
            }
        }

        // Save last sequence
        if let Some(id) = current_id {
            if !current_seq.is_empty() {
                sequences.insert(id, current_seq);
            }
        }
    }

    // Convert to Vec for parallel iteration
    Ok(sequences.into_iter().collect())
}

/// Build mutalyzer-compatible annotations JSON from cdot transcript data.
///
/// Coordinate conversion:
/// - cdot exon coordinates are 1-indexed, inclusive (GenBank-style)
/// - mutalyzer expects 0-indexed, half-open coordinates
/// - Convert: start = cdot_start - 1, end = cdot_end (no change for half-open)
///
/// mol_type is set based on accession prefix:
/// - NM_/XM_ = mRNA (coding transcripts)
/// - NR_/XR_ = ncRNA (non-coding transcripts)
/// - NG_/NC_/LRG_ = genomic DNA
///
/// CDS is only added for coding transcripts (NM_/XM_) that have CDS coordinates.
fn build_annotations_from_cdot(
    tx_id: &str,
    seq_len: usize,
    tx: &CdotTranscript,
) -> serde_json::Value {
    let gene_name = tx.gene_name.clone().unwrap_or_else(|| tx_id.to_string());

    // Determine mol_type, feature_type, and whether this is a coding transcript based on accession prefix
    let (mol_type, feature_type, is_coding) = if tx_id.starts_with("NM_")
        || tx_id.starts_with("XM_")
    {
        ("mRNA", "mRNA", true)
    } else if tx_id.starts_with("NR_") || tx_id.starts_with("XR_") {
        ("ncRNA", "ncRNA", false)
    } else if tx_id.starts_with("NG_") || tx_id.starts_with("NC_") || tx_id.starts_with("LRG_") {
        ("genomic DNA", "gene", false)
    } else {
        ("mRNA", "mRNA", true) // Default to mRNA for unknown types
    };

    // Build exon features from cdot exon data
    // cdot exons (internal format): [genome_start, genome_end, tx_start, tx_end]
    // tx_start and tx_end are 1-indexed inclusive, convert to 0-indexed half-open
    let mut exon_features: Vec<serde_json::Value> = Vec::new();
    for (i, exon) in tx.exons.iter().enumerate() {
        // Convert 1-indexed inclusive to 0-indexed half-open:
        // - start: subtract 1 (1-indexed -> 0-indexed)
        // - end: keep as-is (inclusive -> exclusive)
        let tx_start = hgvs_pos_to_index(exon[2]); // 1-indexed → 0-indexed
        let tx_end = exon[3] as usize; // Don't subtract 1 for half-open end
        exon_features.push(json!({
            "type": "exon",
            "location": {
                "type": "range",
                "start": {"type": "point", "position": tx_start},
                "end": {"type": "point", "position": tx_end},
                "strand": 1
            },
            "id": format!("id-{}-{}", gene_name, i + 1)
        }));
    }

    // Build mRNA/transcript features
    // Only add CDS for coding transcripts that have CDS coordinates defined
    let mut mrna_features: Vec<serde_json::Value> = Vec::new();

    if is_coding && tx.cds_start.is_some() && tx.cds_end.is_some() {
        // cdot CDS positions (start_codon/stop_codon) are 0-indexed
        let cds_start = tx.cds_start.map(|p| p as usize).unwrap_or(0);
        let cds_end = tx.cds_end.map(|p| p as usize).unwrap_or(seq_len);
        let protein_id = tx
            .protein
            .clone()
            .unwrap_or_else(|| format!("{}_protein", tx_id));

        mrna_features.push(json!({
            "type": "CDS",
            "location": {
                "type": "range",
                "start": {"type": "point", "position": cds_start},
                "end": {"type": "point", "position": cds_end},
                "strand": 1
            },
            "id": protein_id,
            "features": []
        }));
    }

    mrna_features.extend(exon_features);

    json!({
        "id": tx_id,
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": seq_len}
        },
        "qualifiers": {
            "name": gene_name,
            "mol_type": mol_type
        },
        "features": [{
            "type": "gene",
            "location": {
                "type": "range",
                "start": {"type": "point", "position": 0},
                "end": {"type": "point", "position": seq_len},
                "strand": 1
            },
            "id": gene_name.clone(),
            "qualifiers": {"name": gene_name},
            "features": [{
                "id": tx_id,
                "type": feature_type,
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 0},
                    "end": {"type": "point", "position": seq_len}
                },
                "features": mrna_features
            }]
        }]
    })
}

/// Load genomic sequences (NC_*, NT_*, NW_*) from genome FASTA files.
/// Includes primary chromosomes (NC_*), contigs (NT_*), and patches/scaffolds (NW_*).
fn load_genomic_sequences<P: AsRef<Path>>(dir: P) -> Result<Vec<(String, String)>, FerroError> {
    let dir = dir.as_ref();
    let mut sequences = Vec::new();

    // Find genome FASTA files (GRCh38.fna, GRCh37.fna, or GCF_* assembly files)
    let fasta_files: Vec<_> = std::fs::read_dir(dir)
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to read {}: {}", dir.display(), e),
        })?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .is_some_and(|e| e == "fna" || e == "fa" || e == "fasta")
                && p.file_name().and_then(|n| n.to_str()).is_some_and(|n| {
                    n.starts_with("GCF_") || n.starts_with("GRCh") || n.contains("genomic")
                })
        })
        .collect();

    for fasta_path in fasta_files {
        eprintln!("  Reading {}...", fasta_path.display());
        let file = File::open(&fasta_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", fasta_path.display(), e),
        })?;
        let reader = BufReader::with_capacity(8 * 1024 * 1024, file); // 8MB buffer for large genome

        let mut current_id: Option<String> = None;
        let mut current_seq = String::new();
        let mut skipped = 0usize;

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read line: {}", e),
            })?;

            if let Some(header) = line.strip_prefix('>') {
                // Save previous sequence if it's a genomic accession (NC_*, NT_*, NW_*)
                if let Some(id) = current_id.take() {
                    if is_genomic_accession(&id) && !current_seq.is_empty() {
                        sequences.push((id, std::mem::take(&mut current_seq)));
                    } else {
                        skipped += 1;
                        current_seq.clear();
                    }
                }
                // Parse new ID - keep NC_*, NT_*, NW_*
                let id = header.split_whitespace().next().map(|s| s.to_string());
                if let Some(ref acc) = id {
                    if is_genomic_accession(acc) {
                        current_id = id;
                    } else {
                        current_id = None;
                    }
                }
            } else if current_id.is_some() {
                current_seq.push_str(line.trim());
            }
        }

        // Save last sequence
        if let Some(id) = current_id {
            if is_genomic_accession(&id) && !current_seq.is_empty() {
                sequences.push((id, current_seq));
            } else {
                skipped += 1;
            }
        }

        if skipped > 0 {
            eprintln!("    Skipped {} non-genomic sequences", skipped);
        }
    }

    Ok(sequences)
}

/// Check if an accession is a genomic sequence type we want to cache.
fn is_genomic_accession(acc: &str) -> bool {
    acc.starts_with("NC_") || acc.starts_with("NT_") || acc.starts_with("NW_")
}

/// Fetch legacy GenBank sequences from ferro reference or NCBI E-utilities.
/// These are old accessions like U31929.1, M75126.1, etc.
///
/// If `ferro_reference` is provided, tries to load from ferro's legacy_genbank.fna first.
pub fn fetch_legacy_genbank_sequences(
    accessions: &[String],
    output_dir: &Path,
    ferro_reference: Option<&Path>,
) -> Result<usize, FerroError> {
    use std::process::Command;

    if accessions.is_empty() {
        return Ok(0);
    }

    std::fs::create_dir_all(output_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create legacy dir: {}", e),
    })?;

    // Load GenBank sequences from ferro's legacy_genbank.fna if available
    let ferro_genbank_sequences: HashMap<String, String> = if let Some(ferro_ref) = ferro_reference
    {
        let genbank_fasta = ferro_ref.join("supplemental/legacy_genbank.fna");
        if genbank_fasta.exists() {
            load_fasta_sequences(&genbank_fasta).unwrap_or_default()
        } else {
            HashMap::new()
        }
    } else {
        HashMap::new()
    };

    if !ferro_genbank_sequences.is_empty() {
        eprintln!(
            "  Found {} GenBank sequences in ferro reference",
            ferro_genbank_sequences.len()
        );
    }

    let mut fetched = 0;
    let mut from_ferro = 0;

    for acc in accessions {
        // Skip if already cached
        let seq_path = output_dir.join(format!("{}.sequence", acc));
        if seq_path.exists() {
            fetched += 1;
            continue;
        }

        // Try to get from ferro's GenBank sequences first
        if let Some(sequence) = ferro_genbank_sequences.get(acc) {
            // Write the sequence to cache
            std::fs::write(&seq_path, sequence).map_err(|e| FerroError::Io {
                msg: format!("Failed to write {}: {}", acc, e),
            })?;

            // Create minimal annotations
            let annotations = build_minimal_annotations(acc, sequence.len());
            let ann_path = output_dir.join(format!("{}.annotations", acc));
            let ann_file = std::fs::File::create(&ann_path).map_err(|e| FerroError::Io {
                msg: format!("Failed to create annotations for {}: {}", acc, e),
            })?;
            serde_json::to_writer(ann_file, &annotations).map_err(|e| FerroError::Io {
                msg: format!("Failed to write annotations for {}: {}", acc, e),
            })?;

            from_ferro += 1;
            fetched += 1;
            continue;
        }

        // Fall back to NCBI fetch
        eprintln!("  Fetching {} from NCBI...", acc);

        // Use efetch to get FASTA
        let url = build_efetch_url("nuccore", acc, "fasta", "text");

        let output = Command::new("curl")
            .args(["-s", "-L", &url])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to fetch {}: {}", acc, e),
            })?;

        if !output.status.success() {
            eprintln!("    Warning: Failed to fetch {}", acc);
            continue;
        }

        let fasta = String::from_utf8_lossy(&output.stdout);

        // Parse FASTA to extract sequence
        let mut seq = String::new();
        for line in fasta.lines() {
            if !line.starts_with('>') && !line.is_empty() {
                seq.push_str(line.trim());
            }
        }

        if seq.is_empty() {
            eprintln!("    Warning: Empty sequence for {}", acc);
            continue;
        }

        // Write sequence file
        std::fs::write(&seq_path, &seq).map_err(|e| FerroError::Io {
            msg: format!("Failed to write {}: {}", seq_path.display(), e),
        })?;

        // Write minimal annotations
        let ann_path = output_dir.join(format!("{}.annotations", acc));
        let annotations = build_minimal_annotations(acc, seq.len());
        let mut ann_file = File::create(&ann_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create {}: {}", ann_path.display(), e),
        })?;
        serde_json::to_writer(&mut ann_file, &annotations).map_err(|e| FerroError::Io {
            msg: format!("Failed to write annotations: {}", e),
        })?;

        fetched += 1;

        // Rate limit to avoid NCBI throttling
        std::thread::sleep(std::time::Duration::from_millis(350));
    }

    if from_ferro > 0 {
        eprintln!(
            "  Loaded {} GenBank sequences from ferro reference",
            from_ferro
        );
    }

    Ok(fetched)
}

/// Known legacy GenBank accessions used in ClinVar.
/// These are fetched on-demand from NCBI.
pub const LEGACY_GENBANK_ACCESSIONS: &[&str] = &[
    "U31929.1",
    "M75126.1",
    "M22590.1",
    "AJ298105.1",
    "AF319569.1",
];

/// Legacy transcript versions not in current cdot (GRCh38).
/// These are older RefSeq versions referenced in ClinVar/test data
/// that have been superseded by newer versions in cdot.
/// Format: (accession_to_fetch, optional_fallback_if_retired)
pub const LEGACY_TRANSCRIPT_VERSIONS: &[(&str, Option<&str>)] = &[
    // Older transcript versions from ClinVar test patterns
    ("NM_000051.3", None),                // ATM - superseded by .4
    ("NM_000174.4", None),                // FSHR - superseded by .5
    ("NM_000246.3", None),                // MEN1 - superseded by .4
    ("NM_000267.3", None),                // NF1 - superseded by .4
    ("NM_000280.5", None),                // PAX6 - superseded by .6
    ("NM_000426.3", None),                // LAMA2 - superseded by .4
    ("NM_000449.3", None),                // RET - superseded by .5
    ("NM_001127221.1", None),             // BBS1
    ("NM_001271223.2", None),             // DNMT3A
    ("NM_001458.4", None),                // FLNC - superseded by .5
    ("NM_002529.3", None),                // NTRK1 - superseded by .4
    ("NM_002878.3", None),                // RAD51D - superseded by .4
    ("NM_004006.2", None),                // DMD - superseded by .3
    ("NM_004082.4", None),                // DCDC2 - superseded by .5
    ("NM_004321.7", None),                // KIF1A - superseded by .8
    ("NM_004364.4", None),                // CEBPA - superseded by .5
    ("NM_015125.4", None),                // PALB2 - superseded by .5
    ("NM_015243.2", None),                // WDR19
    ("NM_015247.2", None),                // NLRP3 - superseded by .3
    ("NM_017890.4", None),                // CPLANE1 - superseded by .5
    ("NM_018400.3", None),                // NFKB2 - superseded by .4
    ("NM_018993.3", None),                // TMEM216 - superseded by .4
    ("NM_019023.4", None),                // ACADVL - superseded by .5
    ("NM_021946.4", None),                // DISC1 - superseded by .5
    ("NM_030962.3", None),                // TRPM4 - superseded by .4
    ("NM_030984.5", None),                // CALN1
    ("NM_032790.3", None),                // PRRC2A - superseded by .4
    ("NM_033517.2", Some("NM_033517.1")), // SHANK3 - .2 retired, use .1 as fallback
    ("NM_138387.3", None),                // G6PC3 - superseded by .4
    ("NM_152383.4", None),                // PRPF31
    ("NM_172201.1", None),                // CHRNA9
    ("NM_177438.2", None),                // ARHGAP31 - superseded by .3
    ("NM_203447.3", None),                // TPCN2 - superseded by .4
];

/// Fetch legacy transcript versions, preferring ferro reference over NCBI.
/// These are older RefSeq versions not in current cdot.
pub fn fetch_legacy_transcript_versions(
    output_dir: &Path,
    ferro_reference: Option<&Path>,
) -> Result<usize, FerroError> {
    use std::process::Command;

    let mut fetched = 0;
    let mut from_ferro = 0;

    // Load legacy transcripts from ferro's legacy_transcripts.fna if available
    let ferro_legacy_sequences: HashMap<String, String> = if let Some(ferro_ref) = ferro_reference {
        let legacy_fasta = ferro_ref.join("supplemental/legacy_transcripts.fna");
        if legacy_fasta.exists() {
            load_fasta_sequences(&legacy_fasta).unwrap_or_default()
        } else {
            HashMap::new()
        }
    } else {
        HashMap::new()
    };

    if !ferro_legacy_sequences.is_empty() {
        eprintln!(
            "  Found {} legacy transcripts in ferro reference",
            ferro_legacy_sequences.len()
        );
    }

    for (accession, fallback) in LEGACY_TRANSCRIPT_VERSIONS {
        // Skip if already cached
        let seq_path = output_dir.join(format!("{}.sequence", accession));
        if seq_path.exists() {
            fetched += 1;
            continue;
        }

        // Try to get from ferro's legacy transcripts first
        if let Some(sequence) = ferro_legacy_sequences.get(*accession) {
            // Write the sequence to cache
            std::fs::write(&seq_path, sequence).map_err(|e| FerroError::Io {
                msg: format!("Failed to write {}: {}", accession, e),
            })?;

            // Create minimal annotations
            let annotations = build_minimal_annotations(accession, sequence.len());
            let ann_path = output_dir.join(format!("{}.annotations", accession));
            let ann_file = std::fs::File::create(&ann_path).map_err(|e| FerroError::Io {
                msg: format!("Failed to create annotations for {}: {}", accession, e),
            })?;
            serde_json::to_writer(ann_file, &annotations).map_err(|e| FerroError::Io {
                msg: format!("Failed to write annotations for {}: {}", accession, e),
            })?;

            from_ferro += 1;
            fetched += 1;
            continue;
        }

        eprintln!("  Fetching {} from NCBI...", accession);

        // Try to fetch the accession
        let url = build_efetch_url("nuccore", accession, "gb", "text");

        let output = Command::new("curl")
            .args(["-s", "-L", &url])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to fetch {}: {}", accession, e),
            })?;

        let genbank_text = String::from_utf8_lossy(&output.stdout);

        // Check if we got an error (e.g., "Error" or empty/invalid response)
        let fetch_failed = !output.status.success()
            || genbank_text.contains("Error")
            || genbank_text.contains("Failed to retrieve")
            || !genbank_text.contains("LOCUS");

        if fetch_failed {
            // If there's a fallback, try to use it
            if let Some(fallback_acc) = fallback {
                eprintln!(
                    "    {} not available, using fallback {}",
                    accession, fallback_acc
                );

                // Check if fallback already exists
                let fallback_seq_path = output_dir.join(format!("{}.sequence", fallback_acc));
                if !fallback_seq_path.exists() {
                    // Fetch the fallback
                    let fallback_url = build_efetch_url("nuccore", fallback_acc, "gb", "text");
                    let fallback_output = Command::new("curl")
                        .args(["-s", "-L", &fallback_url])
                        .output()
                        .map_err(|e| FerroError::Io {
                            msg: format!("Failed to fetch fallback {}: {}", fallback_acc, e),
                        })?;

                    let fallback_genbank = String::from_utf8_lossy(&fallback_output.stdout);
                    if let Some((seq, cds_start, cds_end, gene_name)) =
                        parse_genbank_record(&fallback_genbank)
                    {
                        write_cache_entry(
                            output_dir,
                            fallback_acc,
                            &seq,
                            cds_start,
                            cds_end,
                            &gene_name,
                        )?;
                    }
                    std::thread::sleep(std::time::Duration::from_millis(350));
                }

                // Create symlinks from requested accession to fallback
                let ann_path = output_dir.join(format!("{}.annotations", accession));
                let fallback_seq = format!("{}.sequence", fallback_acc);
                let fallback_ann = format!("{}.annotations", fallback_acc);

                // Use symlinks on Unix
                #[cfg(unix)]
                {
                    use std::os::unix::fs::symlink;
                    let _ = symlink(&fallback_seq, &seq_path);
                    let _ = symlink(&fallback_ann, &ann_path);
                }
                #[cfg(not(unix))]
                {
                    // On non-Unix, copy the files
                    let _ = std::fs::copy(output_dir.join(&fallback_seq), &seq_path);
                    let _ = std::fs::copy(output_dir.join(&fallback_ann), &ann_path);
                }

                fetched += 1;
            } else {
                eprintln!("    Warning: Failed to fetch {}", accession);
            }
            continue;
        }

        // Parse GenBank record
        if let Some((seq, cds_start, cds_end, gene_name)) = parse_genbank_record(&genbank_text) {
            write_cache_entry(output_dir, accession, &seq, cds_start, cds_end, &gene_name)?;
            fetched += 1;
        } else {
            eprintln!("    Warning: Failed to parse GenBank for {}", accession);
        }

        // Rate limit
        std::thread::sleep(std::time::Duration::from_millis(350));
    }

    if from_ferro > 0 {
        eprintln!(
            "  Loaded {} legacy transcripts from ferro reference",
            from_ferro
        );
    }

    Ok(fetched)
}

/// Parse a GenBank record to extract sequence and CDS info.
fn parse_genbank_record(genbank: &str) -> Option<(String, Option<usize>, Option<usize>, String)> {
    let mut sequence = String::new();
    let mut in_origin = false;
    let mut gene_name = String::new();
    let mut cds_start: Option<usize> = None;
    let mut cds_end: Option<usize> = None;

    for line in genbank.lines() {
        if line.starts_with("ORIGIN") {
            in_origin = true;
            continue;
        }

        if in_origin {
            if line.starts_with("//") {
                break;
            }
            // Parse sequence lines (format: "   123 acgtacgt acgtacgt...")
            for part in line.split_whitespace() {
                if part.chars().all(|c| c.is_ascii_alphabetic()) {
                    sequence.push_str(&part.to_uppercase());
                }
            }
        } else {
            // Parse gene name
            if line.contains("/gene=") {
                if let Some(start) = line.find("/gene=\"") {
                    let rest = &line[start + 7..];
                    if let Some(end) = rest.find('"') {
                        gene_name = rest[..end].to_string();
                    }
                }
            }
            // Parse CDS location
            if line.trim_start().starts_with("CDS") {
                // Simple parsing for "CDS             123..456" or "CDS             join(..."
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let loc = parts[1];
                    // Handle simple range "123..456"
                    if let Some(range) = loc.strip_prefix("join(") {
                        // For join(), take first segment
                        if let Some(first_seg) = range.split(',').next() {
                            if let Some((s, e)) = first_seg.split_once("..") {
                                cds_start = s.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                                cds_end = e.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                            }
                        }
                    } else if let Some((s, e)) = loc.split_once("..") {
                        cds_start = s.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                        cds_end = e.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                    }
                }
            }
        }
    }

    if sequence.is_empty() {
        return None;
    }

    // Convert to 0-indexed (1-based → 0-based)
    cds_start = cds_start.map(|p| hgvs_pos_to_index(p as u64));
    cds_end = cds_end.map(|p| hgvs_pos_to_index(p as u64));

    Some((sequence, cds_start, cds_end, gene_name))
}

/// Write a cache entry (sequence + annotations).
fn write_cache_entry(
    output_dir: &Path,
    accession: &str,
    sequence: &str,
    cds_start: Option<usize>,
    cds_end: Option<usize>,
    gene_name: &str,
) -> Result<(), FerroError> {
    let seq_path = output_dir.join(format!("{}.sequence", accession));
    std::fs::write(&seq_path, sequence).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", seq_path.display(), e),
    })?;

    let gene = if gene_name.is_empty() {
        accession.to_string()
    } else {
        gene_name.to_string()
    };

    let seq_len = sequence.len();
    let cds_s = cds_start.unwrap_or(0);
    let cds_e = cds_end.unwrap_or(seq_len);

    let annotations = json!({
        "id": accession,
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": seq_len}
        },
        "qualifiers": {
            "name": gene,
            "mol_type": "mRNA"
        },
        "features": [{
            "type": "gene",
            "location": {
                "type": "range",
                "start": {"type": "point", "position": 0},
                "end": {"type": "point", "position": seq_len},
                "strand": 1
            },
            "id": gene,
            "qualifiers": {"name": gene.clone()},
            "features": [{
                "id": accession,
                "type": "mRNA",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 0},
                    "end": {"type": "point", "position": seq_len}
                },
                "features": [{
                    "type": "CDS",
                    "id": format!("{}_protein", accession),
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": cds_s},
                        "end": {"type": "point", "position": cds_e},
                        "strand": 1
                    },
                    "features": []
                }]
            }]
        }]
    });

    let ann_path = output_dir.join(format!("{}.annotations", accession));
    let ann_file = File::create(&ann_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", ann_path.display(), e),
    })?;
    serde_json::to_writer(ann_file, &annotations).map_err(|e| FerroError::Io {
        msg: format!("Failed to write annotations: {}", e),
    })?;

    Ok(())
}

/// Load RefSeqGene sequences (NG_*) from refseqgene FASTA files.
fn load_refseqgene_sequences<P: AsRef<Path>>(dir: P) -> Result<Vec<(String, String)>, FerroError> {
    let dir = dir.as_ref();
    let mut sequences = Vec::new();

    // Find RefSeqGene FASTA files
    let fasta_files: Vec<_> = std::fs::read_dir(dir)
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to read {}: {}", dir.display(), e),
        })?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .is_some_and(|e| e == "fna" || e == "fa" || e == "fasta")
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.contains("refseqgene") || n.starts_with("NG_"))
        })
        .collect();

    for fasta_path in fasta_files {
        eprintln!(
            "  Reading {}...",
            fasta_path.file_name().unwrap_or_default().to_string_lossy()
        );
        let file = File::open(&fasta_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", fasta_path.display(), e),
        })?;
        let reader = BufReader::with_capacity(4 * 1024 * 1024, file);

        let mut current_id: Option<String> = None;
        let mut current_seq = String::new();

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read line: {}", e),
            })?;

            if let Some(header) = line.strip_prefix('>') {
                // Save previous sequence if it's a RefSeqGene (NG_*)
                if let Some(id) = current_id.take() {
                    if id.starts_with("NG_") && !current_seq.is_empty() {
                        sequences.push((id, std::mem::take(&mut current_seq)));
                    } else {
                        current_seq.clear();
                    }
                }
                // Parse new ID
                let id = header.split_whitespace().next().map(|s| s.to_string());
                if let Some(ref acc) = id {
                    if acc.starts_with("NG_") {
                        current_id = id;
                    } else {
                        current_id = None;
                    }
                }
            } else if current_id.is_some() {
                current_seq.push_str(line.trim());
            }
        }

        // Save last sequence
        if let Some(id) = current_id {
            if id.starts_with("NG_") && !current_seq.is_empty() {
                sequences.push((id, current_seq));
            }
        }
    }

    Ok(sequences)
}

/// Extract all accessions referenced in an HGVS pattern.
///
/// This includes:
/// - Primary reference (e.g., `NM_000059.4` from `NM_000059.4:c.123del`)
/// - Embedded references in `ins[...]`, `delins[...]` (e.g., `PP887427.1` from `ins[PP887427.1:g.1_1518]`)
/// - References in complex alleles
pub fn extract_accessions_from_pattern(pattern: &str) -> Vec<String> {
    let mut accessions = Vec::new();
    let bytes = pattern.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    while i < len {
        // Skip non-ASCII bytes (e.g., en-dash in malformed patterns)
        if !bytes[i].is_ascii() {
            i += 1;
            continue;
        }

        // Check for RefSeq prefixes: NC_, NG_, NM_, NR_, NP_, XM_, XR_, XP_, NT_, NW_, CM_
        // CM_ = GenBank chromosome assemblies (redundant with NC_, but we track them)
        if i + 3 <= len && bytes[i + 1].is_ascii() && bytes[i + 2].is_ascii() {
            let prefix = &bytes[i..i + 3];
            if matches!(
                prefix,
                b"NC_"
                    | b"NG_"
                    | b"NM_"
                    | b"NR_"
                    | b"NP_"
                    | b"XM_"
                    | b"XR_"
                    | b"XP_"
                    | b"NT_"
                    | b"NW_"
                    | b"CM_"
            ) {
                if let Some(acc) = extract_refseq_accession(&pattern[i..]) {
                    accessions.push(acc);
                    i += 3;
                    continue;
                }
            }
        }

        // Check for LRG_ prefix
        if i + 4 <= len
            && bytes[i + 1].is_ascii()
            && bytes[i + 2].is_ascii()
            && bytes[i + 3].is_ascii()
            && &bytes[i..i + 4] == b"LRG_"
        {
            if let Some(acc) = extract_lrg_accession(&pattern[i..]) {
                accessions.push(acc);
                i += 4;
                continue;
            }
        }

        // Check for GenBank accessions: two uppercase letters followed by 6 digits and version
        // e.g., PP887427.1, KY923049.1
        if i + 2 <= len
            && bytes[i].is_ascii_uppercase()
            && bytes[i + 1].is_ascii_uppercase()
            && (i + 2 >= len || !bytes[i + 2].is_ascii_alphabetic())
        {
            if let Some(acc) = extract_genbank_accession(&pattern[i..]) {
                accessions.push(acc);
                i += 2;
                continue;
            }
        }

        i += 1;
    }

    // Deduplicate while preserving order
    let mut seen = std::collections::HashSet::new();
    accessions.retain(|acc| seen.insert(acc.clone()));

    accessions
}

/// Extract a RefSeq accession (e.g., NM_000059.4) starting at the given position.
fn extract_refseq_accession(s: &str) -> Option<String> {
    // Format: XX_NNNNNN.V where XX is prefix, N is digits, V is version
    let bytes = s.as_bytes();
    if bytes.len() < 4 {
        return None;
    }

    // Skip the prefix (already validated)
    let mut i = 3;

    // Read digits
    let digit_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == digit_start {
        return None; // No digits after prefix
    }

    // Expect a dot
    if i >= bytes.len() || bytes[i] != b'.' {
        return None;
    }
    i += 1;

    // Read version digits
    let version_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == version_start {
        return None; // No version digits
    }

    Some(s[..i].to_string())
}

/// Extract an LRG accession (e.g., LRG_292, LRG_292t1, LRG_292p1).
fn extract_lrg_accession(s: &str) -> Option<String> {
    // Format: LRG_N[tN|pN] where N is digits
    let bytes = s.as_bytes();
    if bytes.len() < 5 {
        return None;
    }

    // Skip "LRG_"
    let mut i = 4;

    // Read digits
    let digit_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == digit_start {
        return None; // No digits
    }

    // Optionally read t/p followed by digits
    if i < bytes.len() && (bytes[i] == b't' || bytes[i] == b'p') {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
    }

    Some(s[..i].to_string())
}

/// Extract a GenBank accession (e.g., PP887427.1) starting at the given position.
fn extract_genbank_accession(s: &str) -> Option<String> {
    // Format: LLNNNNNN.V where L is letter, N is digit, V is version
    let bytes = s.as_bytes();
    if bytes.len() < 10 {
        // Minimum: 2 letters + 6 digits + dot + 1 version digit
        return None;
    }

    // First two chars must be uppercase letters (already checked)
    let mut i = 2;

    // Read exactly 6 digits
    let digit_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i - digit_start != 6 {
        return None; // Must be exactly 6 digits for GenBank
    }

    // Expect a dot
    if i >= bytes.len() || bytes[i] != b'.' {
        return None;
    }
    i += 1;

    // Read version digits
    let version_start = i;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    if i == version_start {
        return None; // No version digits
    }

    Some(s[..i].to_string())
}

/// Extract all unique accessions from a file of HGVS patterns.
pub fn extract_all_accessions_from_file<P: AsRef<Path>>(
    pattern_file: P,
) -> Result<Vec<String>, FerroError> {
    let file = File::open(pattern_file.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open pattern file: {}", e),
    })?;
    let reader = BufReader::new(file);

    let mut all_accessions = std::collections::HashSet::new();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        for acc in extract_accessions_from_pattern(line) {
            all_accessions.insert(acc);
        }
    }

    let mut result: Vec<String> = all_accessions.into_iter().collect();
    result.sort();
    Ok(result)
}

/// Fetch missing accessions from NCBI and add them to the cache.
///
/// Checks which accessions are missing from the cache and fetches them from NCBI.
/// Uses batched requests (up to 200 IDs per request) for efficiency.
/// Returns the number of accessions fetched.
pub fn fetch_missing_accessions(
    accessions: &[String],
    cache_dir: &Path,
) -> Result<usize, FerroError> {
    const BATCH_SIZE: usize = 200; // NCBI efetch supports larger batches but 200 is safer

    let mut fetched = 0;
    let mut missing_nuccore: Vec<&String> = Vec::new();
    let mut missing_protein: Vec<&String> = Vec::new();
    let mut skipped_cm = 0usize;
    let mut skipped_lrg = 0usize;
    let mut skipped_nc = 0usize;
    let mut skipped_draft = 0usize;

    // Find which accessions are missing and group by database
    for acc in accessions {
        // Skip CM_ accessions - they are GenBank chromosome assemblies
        // which are redundant with NC_ accessions we already have
        // (e.g., CM000663.1 = NC_000001.10 for GRCh37 chr1)
        if acc.starts_with("CM_") {
            skipped_cm += 1;
            continue;
        }

        // Skip LRG_ accessions - they are downloaded from EBI separately
        if acc.starts_with("LRG_") {
            skipped_lrg += 1;
            continue;
        }

        // Skip NC_ accessions - genomic chromosomes are in genome FASTAs
        if acc.starts_with("NC_") {
            skipped_nc += 1;
            continue;
        }

        // Skip NT_/NW_ accessions - draft genomic scaffolds, rarely needed
        if acc.starts_with("NT_") || acc.starts_with("NW_") {
            skipped_draft += 1;
            continue;
        }

        let seq_path = cache_dir.join(format!("{}.sequence", acc));
        if !seq_path.exists() {
            // Protein accessions go to the protein database
            if acc.starts_with("NP_") || acc.starts_with("XP_") || acc.starts_with("YP_") {
                missing_protein.push(acc);
            } else {
                missing_nuccore.push(acc);
            }
        }
    }

    // Log skipped accessions
    if skipped_cm > 0 {
        eprintln!(
            "  Skipped {} CM_ accessions (redundant with NC_ chromosomes)",
            skipped_cm
        );
    }
    if skipped_lrg > 0 {
        eprintln!(
            "  Skipped {} LRG_ accessions (fetched separately from EBI)",
            skipped_lrg
        );
    }
    if skipped_nc > 0 {
        eprintln!("  Skipped {} NC_ accessions (in genome FASTAs)", skipped_nc);
    }
    if skipped_draft > 0 {
        eprintln!(
            "  Skipped {} NT_/NW_ accessions (draft scaffolds)",
            skipped_draft
        );
    }

    if missing_nuccore.is_empty() && missing_protein.is_empty() {
        eprintln!("  All accessions already in cache");
        return Ok(0);
    }

    // Summary of what we'll fetch
    eprintln!(
        "\n=== Fetch Summary ===\n  Nucleotide: {} missing\n  Protein:    {} missing\n  Total:      {} to fetch\n",
        missing_nuccore.len(),
        missing_protein.len(),
        missing_nuccore.len() + missing_protein.len()
    );

    // Fetch nuccore accessions in batches
    if !missing_nuccore.is_empty() {
        let total_batches = missing_nuccore.len().div_ceil(BATCH_SIZE);
        eprintln!(
            "Fetching {} nucleotide accessions in {} batches...",
            missing_nuccore.len(),
            total_batches
        );
        let mut nuc_fetched = 0usize;
        for (batch_idx, batch) in missing_nuccore.chunks(BATCH_SIZE).enumerate() {
            let batch_fetched = fetch_batch(cache_dir, "nuccore", batch)?;
            nuc_fetched += batch_fetched;
            fetched += batch_fetched;
            eprintln!(
                "  Batch {}/{}: fetched {}/{} (total: {}/{})",
                batch_idx + 1,
                total_batches,
                batch_fetched,
                batch.len(),
                nuc_fetched,
                missing_nuccore.len()
            );
            // Rate limit between batches (NCBI allows 3 req/sec without API key)
            std::thread::sleep(std::time::Duration::from_millis(350));
        }
        if nuc_fetched < missing_nuccore.len() {
            eprintln!(
                "  Warning: {} nucleotide accessions could not be fetched",
                missing_nuccore.len() - nuc_fetched
            );
        }
    }

    // Skip protein fetching - proteins are derived from transcript CDS translation
    if !missing_protein.is_empty() {
        eprintln!(
            "\n  Skipping {} protein accessions (derived from transcript CDS)",
            missing_protein.len()
        );
    }

    eprintln!("\n=== Fetch Complete ===");
    eprintln!("  Total fetched: {}", fetched);
    Ok(fetched)
}

/// Fetch a batch of accessions from NCBI and write to cache.
fn fetch_batch(cache_dir: &Path, db: &str, accessions: &[&String]) -> Result<usize, FerroError> {
    use std::process::Command;

    if accessions.is_empty() {
        return Ok(0);
    }

    // Build comma-separated ID list
    let ids: Vec<&str> = accessions.iter().map(|s| s.as_str()).collect();
    let id_list = ids.join(",");

    // Use efetch to get GenBank format for all IDs at once
    let url = build_efetch_url(db, &id_list, "gb", "text");

    let output = Command::new("curl")
        .args(["-s", "-L", &url])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to fetch batch: {}", e),
        })?;

    if !output.status.success() {
        eprintln!("    Warning: Batch fetch failed");
        return Ok(0);
    }

    let genbank_text = String::from_utf8_lossy(&output.stdout);

    // Track which accessions we successfully fetched vs need FASTA fallback
    let mut fetched = 0;
    let mut need_fasta: Vec<String> = Vec::new();

    // Split response into individual records (each starts with "LOCUS")
    for record in genbank_text.split("\nLOCUS ") {
        let record = if record.starts_with("LOCUS ") {
            record.to_string()
        } else if fetched == 0 && need_fasta.is_empty() && record.contains("LOCUS ") {
            // First record might not have leading newline
            record.to_string()
        } else {
            format!("LOCUS {}", record)
        };

        // Extract accession from VERSION line
        let accession = record
            .lines()
            .find(|l| l.starts_with("VERSION"))
            .and_then(|l| l.split_whitespace().nth(1))
            .map(|s| s.to_string());

        if let Some(acc) = accession {
            // Records with CONTIG (e.g., NG_) don't have inline sequence
            if !record.contains("ORIGIN") {
                need_fasta.push(acc);
                continue;
            }

            if let Some((seq, cds_start, cds_end, gene_name)) = parse_genbank_record(&record) {
                if write_cache_entry(cache_dir, &acc, &seq, cds_start, cds_end, &gene_name).is_ok()
                {
                    fetched += 1;
                }
            }
        }
    }

    // Fetch sequences for CONTIG-based records (NG_, etc.) using FASTA
    if !need_fasta.is_empty() {
        let fasta_fetched = fetch_fasta_fallback(cache_dir, db, &need_fasta)?;
        fetched += fasta_fetched;
    }

    Ok(fetched)
}

/// Fetch sequences via FASTA for records that use CONTIG (no inline sequence).
fn fetch_fasta_fallback(
    cache_dir: &Path,
    db: &str,
    accessions: &[String],
) -> Result<usize, FerroError> {
    use std::process::Command;

    if accessions.is_empty() {
        return Ok(0);
    }

    let id_list = accessions.join(",");
    let url = build_efetch_url(db, &id_list, "fasta", "text");

    let output = Command::new("curl")
        .args(["-s", "-L", &url])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to fetch FASTA: {}", e),
        })?;

    if !output.status.success() {
        return Ok(0);
    }

    let fasta_text = String::from_utf8_lossy(&output.stdout);

    // Parse FASTA records
    let mut fetched = 0;
    let mut current_acc: Option<String> = None;
    let mut current_seq = String::new();

    for line in fasta_text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            // Save previous record if any
            if let Some(acc) = current_acc.take() {
                if !current_seq.is_empty()
                    && write_nucleotide_fasta_entry(cache_dir, &acc, &current_seq).is_ok()
                {
                    fetched += 1;
                }
            }
            // Parse header - format varies but accession is usually first word
            current_acc = header.split_whitespace().next().map(|s| s.to_string());
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.push_str(line.trim());
        }
    }

    // Don't forget the last record
    if let Some(acc) = current_acc {
        if !current_seq.is_empty()
            && write_nucleotide_fasta_entry(cache_dir, &acc, &current_seq).is_ok()
        {
            fetched += 1;
        }
    }

    Ok(fetched)
}

/// Supplemental transcript metadata extracted from GenBank records.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SupplementalTranscript {
    /// Transcript accession (e.g., "NM_033517.1")
    pub id: String,
    /// Gene symbol (e.g., "SHANK3")
    pub gene_symbol: Option<String>,
    /// CDS start position (1-based, inclusive)
    pub cds_start: Option<u64>,
    /// CDS end position (1-based, inclusive)
    pub cds_end: Option<u64>,
    /// Sequence length
    pub sequence_length: usize,
}

/// Supplemental metadata for transcripts not in cdot.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SupplementalMetadata {
    /// When this metadata was generated
    pub generated_at: String,
    /// Transcript metadata keyed by accession
    pub transcripts: std::collections::HashMap<String, SupplementalTranscript>,
}

/// Fetch sequences from NCBI in GenBank format, extract CDS metadata, and write to files.
///
/// This fetches GenBank format (not FASTA) to get CDS coordinates along with sequences.
/// Writes both a FASTA file and a supplemental metadata JSON file.
///
/// Returns the number of sequences successfully fetched.
pub fn fetch_fasta_to_file(
    accessions: &[String],
    output_file: &Path,
    dry_run: bool,
) -> Result<usize, FerroError> {
    use std::io::Write;
    use std::process::Command;

    const BATCH_SIZE: usize = 200; // Larger batches for efficiency (4x fewer requests)

    if accessions.is_empty() {
        return Ok(0);
    }

    // Separate nucleotide and protein accessions
    let mut nuccore: Vec<&String> = Vec::new();
    let mut protein: Vec<&String> = Vec::new();
    let mut skipped_lrg = 0usize;
    let mut skipped_cm = 0usize;

    for acc in accessions {
        if acc.starts_with("LRG_") {
            skipped_lrg += 1;
            continue;
        }
        if acc.starts_with("CM_") {
            skipped_cm += 1;
            continue;
        }
        if acc.starts_with("NP_") || acc.starts_with("XP_") || acc.starts_with("YP_") {
            protein.push(acc);
        } else {
            nuccore.push(acc);
        }
    }

    if skipped_lrg > 0 {
        eprintln!(
            "  Skipped {} LRG_ accessions (not available from NCBI efetch)",
            skipped_lrg
        );
    }
    if skipped_cm > 0 {
        eprintln!(
            "  Skipped {} CM_ accessions (redundant with NC_ chromosomes)",
            skipped_cm
        );
    }

    // Skip protein accessions for now - ferro primarily needs nucleotide
    if !protein.is_empty() {
        eprintln!(
            "  Skipped {} protein accessions (NP_/XP_/YP_) - not needed for ferro",
            protein.len()
        );
    }

    if nuccore.is_empty() {
        eprintln!("  No nucleotide accessions to fetch");
        return Ok(0);
    }

    if dry_run {
        eprintln!(
            "  Would fetch {} nucleotide sequences with CDS metadata",
            nuccore.len()
        );
        return Ok(0);
    }

    // Create/open output file for appending
    let mut fasta_file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(output_file)
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", output_file.display(), e),
        })?;

    // Metadata file path (same directory, different name)
    let metadata_path = output_file.with_extension("metadata.json");

    // Load existing metadata if present
    let mut metadata = if metadata_path.exists() {
        let file = File::open(&metadata_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", metadata_path.display(), e),
        })?;
        serde_json::from_reader(file).unwrap_or_else(|_| SupplementalMetadata {
            generated_at: chrono::Utc::now().to_rfc3339(),
            transcripts: std::collections::HashMap::new(),
        })
    } else {
        SupplementalMetadata {
            generated_at: chrono::Utc::now().to_rfc3339(),
            transcripts: std::collections::HashMap::new(),
        }
    };

    // Filter out already-fetched accessions
    let original_count = nuccore.len();
    let nuccore: Vec<&String> = nuccore
        .into_iter()
        .filter(|acc| !metadata.transcripts.contains_key(*acc))
        .collect();

    if nuccore.is_empty() {
        eprintln!(
            "  All {} accessions already fetched (skipping)",
            original_count
        );
        return Ok(0);
    }

    if original_count != nuccore.len() {
        eprintln!(
            "  Skipping {} already-fetched accessions, {} remaining to fetch",
            original_count - nuccore.len(),
            nuccore.len()
        );
    }

    let mut fetched = 0usize;
    let total_batches = nuccore.len().div_ceil(BATCH_SIZE);

    // Rate limiting: 3 req/sec without API key, 10 req/sec with API key
    let has_api_key = get_ncbi_api_key().is_some();
    let rate_delay_ms = if has_api_key { 100 } else { 350 };

    eprintln!(
        "  Fetching {} nucleotide sequences in {} batches (batch size: {})...",
        nuccore.len(),
        total_batches,
        BATCH_SIZE
    );
    if has_api_key {
        eprintln!("  Using NCBI API key (10 req/sec rate limit)");
    } else {
        eprintln!("  No NCBI_API_KEY set (3 req/sec rate limit - set key for 3x faster)");
    }

    let start_time = std::time::Instant::now();

    for (batch_idx, batch) in nuccore.chunks(BATCH_SIZE).enumerate() {
        let id_list: Vec<&str> = batch.iter().map(|s| s.as_str()).collect();
        let ids = id_list.join(",");

        // Track which accessions we successfully fetch
        let batch_accessions: std::collections::HashSet<&str> = id_list.iter().copied().collect();
        let mut fetched_accessions: std::collections::HashSet<String> =
            std::collections::HashSet::new();

        // Fetch GenBank format to get CDS coordinates
        let url = build_efetch_url("nuccore", &ids, "gb", "text");

        let output = Command::new("curl")
            .args(["-s", "-L", &url])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to fetch batch: {}", e),
            })?;

        if !output.status.success() {
            eprintln!(
                "    Warning: Batch {}/{} failed",
                batch_idx + 1,
                total_batches
            );
            continue;
        }

        let genbank_text = String::from_utf8_lossy(&output.stdout);

        // Parse GenBank records
        let mut batch_count = 0;
        for record in genbank_text.split("\nLOCUS ") {
            // First record keeps original format; subsequent records need "LOCUS " prefix
            let record = if record.starts_with("LOCUS ")
                || (batch_count == 0 && record.contains("LOCUS "))
            {
                record.to_string()
            } else {
                format!("LOCUS {}", record)
            };

            // Extract accession from VERSION line
            let accession = record
                .lines()
                .find(|l| l.starts_with("VERSION"))
                .and_then(|l| l.split_whitespace().nth(1))
                .map(|s| s.to_string());

            if let Some(acc) = accession {
                if let Some((seq, cds_start, cds_end, gene_name)) =
                    parse_genbank_record_full(&record)
                {
                    // Write FASTA entry
                    let header = if gene_name.is_empty() {
                        format!(">{}", acc)
                    } else {
                        format!(">{} {}", acc, gene_name)
                    };
                    writeln!(fasta_file, "{}", header).map_err(|e| FerroError::Io {
                        msg: format!("Failed to write to {}: {}", output_file.display(), e),
                    })?;

                    // Write sequence in 70-char lines
                    for chunk in seq.as_bytes().chunks(70) {
                        writeln!(fasta_file, "{}", std::str::from_utf8(chunk).unwrap()).map_err(
                            |e| FerroError::Io {
                                msg: format!("Failed to write to {}: {}", output_file.display(), e),
                            },
                        )?;
                    }

                    // Add metadata entry
                    metadata.transcripts.insert(
                        acc.clone(),
                        SupplementalTranscript {
                            id: acc.clone(),
                            gene_symbol: if gene_name.is_empty() {
                                None
                            } else {
                                Some(gene_name)
                            },
                            cds_start,
                            cds_end,
                            sequence_length: seq.len(),
                        },
                    );

                    fetched_accessions.insert(acc);
                    batch_count += 1;
                }
            }
        }

        // Fallback: For CON records (no ORIGIN section), fetch FASTA format
        // This handles NG_ RefSeqGene records which are constructed from components
        let missing_from_batch: Vec<&str> = batch_accessions
            .iter()
            .filter(|acc| !fetched_accessions.contains(**acc))
            .copied()
            .collect();

        if !missing_from_batch.is_empty() {
            // Fetch missing accessions using FASTA format
            let missing_ids = missing_from_batch.join(",");
            let fasta_url = build_efetch_url("nuccore", &missing_ids, "fasta", "text");

            if let Ok(fasta_output) = Command::new("curl").args(["-s", "-L", &fasta_url]).output() {
                if fasta_output.status.success() {
                    let fasta_text = String::from_utf8_lossy(&fasta_output.stdout);

                    // Parse FASTA records
                    let mut current_acc = String::new();
                    let mut current_seq = String::new();
                    let mut current_desc = String::new();

                    for line in fasta_text.lines() {
                        if let Some(header) = line.strip_prefix('>') {
                            // Save previous record if any
                            if !current_acc.is_empty() && !current_seq.is_empty() {
                                // Write to FASTA file
                                if current_desc.is_empty() {
                                    writeln!(fasta_file, ">{}", current_acc).ok();
                                } else {
                                    writeln!(fasta_file, ">{} {}", current_acc, current_desc).ok();
                                }
                                for chunk in current_seq.as_bytes().chunks(70) {
                                    writeln!(fasta_file, "{}", std::str::from_utf8(chunk).unwrap())
                                        .ok();
                                }

                                // Add metadata (no CDS info for CON records)
                                metadata.transcripts.insert(
                                    current_acc.clone(),
                                    SupplementalTranscript {
                                        id: current_acc.clone(),
                                        gene_symbol: if current_desc.is_empty() {
                                            None
                                        } else {
                                            // Extract gene name from description if present
                                            current_desc
                                                .split('(')
                                                .next()
                                                .map(|s| s.trim().to_string())
                                        },
                                        cds_start: None,
                                        cds_end: None,
                                        sequence_length: current_seq.len(),
                                    },
                                );

                                batch_count += 1;
                            }

                            // Parse new header: >ACC.VER description
                            let parts: Vec<&str> = header.splitn(2, ' ').collect();
                            current_acc = parts[0].to_string();
                            current_desc = parts.get(1).unwrap_or(&"").to_string();
                            current_seq.clear();
                        } else {
                            current_seq.push_str(line.trim());
                        }
                    }

                    // Save last record
                    if !current_acc.is_empty() && !current_seq.is_empty() {
                        if current_desc.is_empty() {
                            writeln!(fasta_file, ">{}", current_acc).ok();
                        } else {
                            writeln!(fasta_file, ">{} {}", current_acc, current_desc).ok();
                        }
                        for chunk in current_seq.as_bytes().chunks(70) {
                            writeln!(fasta_file, "{}", std::str::from_utf8(chunk).unwrap()).ok();
                        }

                        metadata.transcripts.insert(
                            current_acc.clone(),
                            SupplementalTranscript {
                                id: current_acc.clone(),
                                gene_symbol: if current_desc.is_empty() {
                                    None
                                } else {
                                    current_desc.split('(').next().map(|s| s.trim().to_string())
                                },
                                cds_start: None,
                                cds_end: None,
                                sequence_length: current_seq.len(),
                            },
                        );

                        batch_count += 1;
                    }
                }
            }

            // Small delay between GenBank and FASTA requests
            std::thread::sleep(std::time::Duration::from_millis(rate_delay_ms));
        }

        fetched += batch_count;

        // Progress with ETA every 10 batches or at the end
        if (batch_idx + 1) % 10 == 0 || batch_idx + 1 == total_batches {
            let elapsed = start_time.elapsed().as_secs_f64();
            let batches_done = batch_idx + 1;
            let batches_remaining = total_batches - batches_done;
            let secs_per_batch = elapsed / batches_done as f64;
            let eta_secs = (batches_remaining as f64 * secs_per_batch) as u64;
            let eta_min = eta_secs / 60;
            let eta_sec = eta_secs % 60;

            eprintln!(
                "    Batch {}/{} ({:.1}%): {} sequences fetched, ETA: {}m {}s",
                batches_done,
                total_batches,
                (batches_done as f64 / total_batches as f64) * 100.0,
                fetched,
                eta_min,
                eta_sec
            );
        }

        // Rate limiting - respect NCBI limits
        if batch_idx < total_batches - 1 {
            std::thread::sleep(std::time::Duration::from_millis(rate_delay_ms));
        }
    }

    // Save metadata
    let metadata_file = File::create(&metadata_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", metadata_path.display(), e),
    })?;
    serde_json::to_writer_pretty(metadata_file, &metadata).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", metadata_path.display(), e),
    })?;
    eprintln!("  Wrote metadata to {}", metadata_path.display());

    Ok(fetched)
}

/// Parse a GenBank record to extract sequence and CDS info (1-based coordinates).
///
/// Returns (sequence, cds_start, cds_end, gene_name) with 1-based coordinates.
fn parse_genbank_record_full(genbank: &str) -> Option<(String, Option<u64>, Option<u64>, String)> {
    let mut sequence = String::new();
    let mut in_origin = false;
    let mut gene_name = String::new();
    let mut cds_start: Option<u64> = None;
    let mut cds_end: Option<u64> = None;

    for line in genbank.lines() {
        if line.starts_with("ORIGIN") {
            in_origin = true;
            continue;
        }

        if in_origin {
            if line.starts_with("//") {
                break;
            }
            // Parse sequence lines (format: "   123 acgtacgt acgtacgt...")
            for part in line.split_whitespace() {
                if part.chars().all(|c| c.is_ascii_alphabetic()) {
                    sequence.push_str(&part.to_uppercase());
                }
            }
        } else {
            // Parse gene name
            if line.contains("/gene=") {
                if let Some(start) = line.find("/gene=\"") {
                    let rest = &line[start + 7..];
                    if let Some(end) = rest.find('"') {
                        gene_name = rest[..end].to_string();
                    }
                }
            }
            // Parse CDS location - we want 1-based coordinates
            if line.trim_start().starts_with("CDS") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let loc = parts[1];
                    // Handle complement() wrapper
                    let loc = loc
                        .strip_prefix("complement(")
                        .and_then(|s| s.strip_suffix(')'))
                        .unwrap_or(loc);

                    // Handle join() for multi-exon CDS - get overall range
                    if let Some(range) = loc.strip_prefix("join(") {
                        let range = range.strip_suffix(')').unwrap_or(range);
                        let segments: Vec<&str> = range.split(',').collect();
                        if let Some(first) = segments.first() {
                            if let Some((s, _)) = first.split_once("..") {
                                cds_start = s.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                            }
                        }
                        if let Some(last) = segments.last() {
                            if let Some((_, e)) = last.split_once("..") {
                                cds_end = e.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                            }
                        }
                    } else if let Some((s, e)) = loc.split_once("..") {
                        cds_start = s.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                        cds_end = e.trim_matches(|c: char| !c.is_numeric()).parse().ok();
                    }
                }
            }
        }
    }

    if sequence.is_empty() {
        return None;
    }

    Some((sequence, cds_start, cds_end, gene_name))
}

/// Write a nucleotide sequence from FASTA to cache with minimal annotations.
fn write_nucleotide_fasta_entry(
    cache_dir: &Path,
    accession: &str,
    sequence: &str,
) -> Result<(), FerroError> {
    // Write sequence file
    let seq_path = cache_dir.join(format!("{}.sequence", accession));
    std::fs::write(&seq_path, sequence).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", seq_path.display(), e),
    })?;

    // Write minimal annotations
    let annotations = build_minimal_annotations(accession, sequence.len());
    let ann_path = cache_dir.join(format!("{}.annotations", accession));
    let ann_json = serde_json::to_string_pretty(&annotations).map_err(|e| FerroError::Json {
        msg: format!("Failed to serialize annotations: {}", e),
    })?;
    std::fs::write(&ann_path, ann_json).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", ann_path.display(), e),
    })?;

    Ok(())
}

/// Fetch a batch of protein accessions from NCBI and write to cache.
///
/// Uses FASTA format since protein records don't need CDS coordinates.
#[allow(dead_code)]
fn fetch_protein_batch(cache_dir: &Path, accessions: &[&String]) -> Result<usize, FerroError> {
    use std::process::Command;

    if accessions.is_empty() {
        return Ok(0);
    }

    // Build comma-separated ID list
    let ids: Vec<&str> = accessions.iter().map(|s| s.as_str()).collect();
    let id_list = ids.join(",");

    // Use efetch to get FASTA format for proteins
    let url = build_efetch_url("protein", &id_list, "fasta", "text");

    let output = Command::new("curl")
        .args(["-s", "-L", &url])
        .output()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to fetch protein batch: {}", e),
        })?;

    if !output.status.success() {
        eprintln!("    Warning: Protein batch fetch failed");
        return Ok(0);
    }

    let fasta_text = String::from_utf8_lossy(&output.stdout);

    // Parse FASTA records
    let mut fetched = 0;
    let mut current_acc: Option<String> = None;
    let mut current_seq = String::new();

    for line in fasta_text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            // Save previous record if any
            if let Some(acc) = current_acc.take() {
                if !current_seq.is_empty()
                    && write_protein_cache_entry(cache_dir, &acc, &current_seq).is_ok()
                {
                    fetched += 1;
                }
            }
            // Parse new header - format: >NP_000005.3 alpha-2-macroglobulin...
            // Extract accession (first word after >)
            current_acc = header.split_whitespace().next().map(|s| s.to_string());
            current_seq.clear();
        } else if !line.is_empty() {
            current_seq.push_str(line.trim());
        }
    }

    // Don't forget the last record
    if let Some(acc) = current_acc {
        if !current_seq.is_empty()
            && write_protein_cache_entry(cache_dir, &acc, &current_seq).is_ok()
        {
            fetched += 1;
        }
    }

    Ok(fetched)
}

/// Write a protein sequence to the cache with annotations.
#[allow(dead_code)]
fn write_protein_cache_entry(
    cache_dir: &Path,
    accession: &str,
    sequence: &str,
) -> Result<(), FerroError> {
    // Write sequence file
    let seq_path = cache_dir.join(format!("{}.sequence", accession));
    std::fs::write(&seq_path, sequence).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", seq_path.display(), e),
    })?;

    // Write annotations file
    let annotations = build_protein_annotations(accession, sequence.len());
    let ann_path = cache_dir.join(format!("{}.annotations", accession));
    let ann_json = serde_json::to_string_pretty(&annotations).map_err(|e| FerroError::Json {
        msg: format!("Failed to serialize annotations: {}", e),
    })?;
    std::fs::write(&ann_path, ann_json).map_err(|e| FerroError::Io {
        msg: format!("Failed to write {}: {}", ann_path.display(), e),
    })?;

    Ok(())
}

/// Build minimal annotations when no cdot metadata is available.
///
/// mol_type is set based on accession prefix:
/// - NM_/XM_ = mRNA (coding transcripts)
/// - NR_/XR_ = ncRNA (non-coding transcripts)
/// - NG_/NC_/LRG_/chr* = genomic DNA
fn build_minimal_annotations(tx_id: &str, seq_len: usize) -> serde_json::Value {
    // Determine mol_type and feature type based on accession prefix
    let (mol_type, feature_type) = if tx_id.starts_with("NM_") || tx_id.starts_with("XM_") {
        ("mRNA", "mRNA")
    } else if tx_id.starts_with("NR_") || tx_id.starts_with("XR_") {
        ("ncRNA", "ncRNA")
    } else if tx_id.starts_with("NG_")
        || tx_id.starts_with("NC_")
        || tx_id.starts_with("LRG_")
        || tx_id.starts_with("chr")
    {
        ("genomic DNA", "gene")
    } else {
        ("mRNA", "mRNA") // Default to mRNA for unknown types
    };

    json!({
        "id": tx_id,
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": seq_len}
        },
        "qualifiers": {
            "name": tx_id,
            "mol_type": mol_type
        },
        "features": [{
            "type": "gene",
            "location": {
                "type": "range",
                "start": {"type": "point", "position": 0},
                "end": {"type": "point", "position": seq_len},
                "strand": 1
            },
            "id": tx_id,
            "features": [{
                "id": tx_id,
                "type": feature_type,
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 0},
                    "end": {"type": "point", "position": seq_len}
                },
                "features": []
            }]
        }]
    })
}

/// Build protein annotations for mutalyzer cache.
///
/// Creates annotations in the format expected by mutalyzer for protein sequences.
/// The format is simpler than nucleotide annotations - just a polypeptide feature.
fn build_protein_annotations(protein_id: &str, seq_len: usize) -> serde_json::Value {
    // Note: mol_type must be "CDS" (not "protein") for mutalyzer to recognize
    // this as a protein coordinate system. Mutalyzer's coordinate_system_from_mol_type
    // only maps "CDS" and "unassigned" to "p" coordinate system.
    json!({
        "id": protein_id,
        "type": "record",
        "coordinate_system": "p",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": seq_len}
        },
        "qualifiers": {
            "name": protein_id,
            "mol_type": "CDS"
        },
        "features": [{
            "type": "CDS",
            "location": {
                "type": "range",
                "start": {"type": "point", "position": 0},
                "end": {"type": "point", "position": seq_len},
                "strand": 1
            },
            "id": protein_id,
            "qualifiers": {
                "gbkey": "Prot"
            }
        }]
    })
}

/// Build NC_ chromosome annotations with transcript selectors from cdot.
///
/// This creates annotations that allow mutalyzer to normalize intronic variants
/// with genomic context like `NC_000006.12(NM_001318856.2):c.9-1296T>C`.
///
/// The annotation includes all transcripts that map to the chromosome as selectors,
/// with their exon coordinates in genomic space.
fn build_nc_annotations_with_selectors(
    nc_id: &str,
    seq_len: usize,
    transcripts: &[(&str, &CdotTranscript)],
) -> serde_json::Value {
    // Group transcripts by gene
    let mut genes: HashMap<String, Vec<(&str, &CdotTranscript)>> = HashMap::new();
    for (tx_id, tx) in transcripts {
        let gene_name = tx.gene_name.clone().unwrap_or_else(|| tx_id.to_string());
        genes.entry(gene_name).or_default().push((tx_id, tx));
    }

    // Build gene features
    let mut gene_features: Vec<serde_json::Value> = Vec::new();

    for (gene_name, gene_transcripts) in genes {
        // Calculate gene boundaries from all transcript exons
        let mut gene_start = u64::MAX;
        let mut gene_end = 0u64;
        let mut gene_strand = 1i32;

        for (_, tx) in &gene_transcripts {
            gene_strand = if tx.strand == crate::reference::Strand::Plus {
                1
            } else {
                -1
            };
            for exon in &tx.exons {
                gene_start = gene_start.min(exon[0]);
                gene_end = gene_end.max(exon[1]);
            }
        }

        // Build transcript features (mRNA/ncRNA) for this gene
        let mut transcript_features: Vec<serde_json::Value> = Vec::new();

        for (tx_id, tx) in &gene_transcripts {
            // Determine transcript type
            let (tx_type, is_coding) = if tx_id.starts_with("NM_") || tx_id.starts_with("XM_") {
                ("mRNA", true)
            } else if tx_id.starts_with("NR_") || tx_id.starts_with("XR_") {
                ("ncRNA", false)
            } else {
                ("mRNA", true)
            };

            // Calculate transcript boundaries
            let mut tx_start = u64::MAX;
            let mut tx_end = 0u64;
            for exon in &tx.exons {
                tx_start = tx_start.min(exon[0]);
                tx_end = tx_end.max(exon[1]);
            }

            // Build exon features in genomic coordinates
            let mut exon_features: Vec<serde_json::Value> = Vec::new();
            for (i, exon) in tx.exons.iter().enumerate() {
                let exon_start = exon[0] as usize;
                let exon_end = exon[1] as usize;
                exon_features.push(json!({
                    "type": "exon",
                    "id": format!("id-{}-{}", tx_id, i + 1),
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": exon_start},
                        "end": {"type": "point", "position": exon_end},
                        "strand": gene_strand
                    }
                }));
            }

            // Build CDS if coding transcript with CDS coordinates
            let mut tx_inner_features: Vec<serde_json::Value> = Vec::new();
            if is_coding {
                if let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) {
                    // CDS in cdot is in transcript coordinates, we need genomic
                    // For now, calculate from exons
                    let cds_genomic = calculate_genomic_cds(tx, cds_start, cds_end);
                    if let Some((cds_g_start, cds_g_end)) = cds_genomic {
                        let protein_id = tx
                            .protein
                            .clone()
                            .unwrap_or_else(|| format!("{}_protein", tx_id));
                        tx_inner_features.push(json!({
                            "type": "CDS",
                            "id": protein_id,
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": cds_g_start},
                                "end": {"type": "point", "position": cds_g_end},
                                "strand": gene_strand
                            },
                            "features": []
                        }));
                    }
                }
            }

            tx_inner_features.extend(exon_features);

            transcript_features.push(json!({
                "type": tx_type,
                "id": *tx_id,
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": tx_start as usize},
                    "end": {"type": "point", "position": tx_end as usize},
                    "strand": gene_strand
                },
                "features": tx_inner_features
            }));
        }

        gene_features.push(json!({
            "type": "gene",
            "id": gene_name.clone(),
            "qualifiers": {"name": gene_name},
            "location": {
                "type": "range",
                "start": {"type": "point", "position": gene_start as usize},
                "end": {"type": "point", "position": gene_end as usize},
                "strand": gene_strand
            },
            "features": transcript_features
        }));
    }

    json!({
        "id": nc_id,
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": seq_len}
        },
        "qualifiers": {
            "name": nc_id,
            "mol_type": "genomic DNA"
        },
        "features": gene_features
    })
}

/// Calculate genomic CDS coordinates from transcript CDS and exon mapping.
fn calculate_genomic_cds(
    tx: &CdotTranscript,
    cds_start: u64,
    cds_end: u64,
) -> Option<(usize, usize)> {
    // cds_start and cds_end are in transcript coordinates (0-indexed)
    // We need to map them to genomic coordinates using the exon boundaries

    // Sort exons by transcript position
    let mut exons: Vec<_> = tx.exons.iter().collect();
    exons.sort_by_key(|e| e[2]); // Sort by tx_start

    let mut g_cds_start: Option<usize> = None;
    let mut g_cds_end: Option<usize> = None;

    for exon in &exons {
        let g_start = exon[0] as usize;
        let g_end = exon[1] as usize;
        let t_start = exon[2]; // 1-indexed in cdot
        let t_end = exon[3]; // 1-indexed in cdot

        // Convert to 0-indexed transcript coordinates (1-based → 0-based)
        let t_start_0 = hgvs_pos_to_index(t_start) as u64;
        let t_end_0 = hgvs_pos_to_index(t_end) as u64;

        // Check if CDS start falls in this exon
        if g_cds_start.is_none() && cds_start >= t_start_0 && cds_start <= t_end_0 {
            let offset = (cds_start - t_start_0) as usize;
            if tx.strand == crate::reference::Strand::Plus {
                g_cds_start = Some(g_start + offset);
            } else {
                g_cds_start = Some(g_end - offset);
            }
        }

        // Check if CDS end falls in this exon
        if g_cds_end.is_none() && cds_end >= t_start_0 && cds_end <= t_end_0 {
            let offset = (cds_end - t_start_0) as usize;
            if tx.strand == crate::reference::Strand::Plus {
                g_cds_end = Some(g_start + offset);
            } else {
                g_cds_end = Some(g_end - offset);
            }
        }
    }

    match (g_cds_start, g_cds_end) {
        (Some(s), Some(e)) => {
            if s < e {
                Some((s, e))
            } else {
                Some((e, s))
            }
        }
        _ => None,
    }
}

/// Enhance NC_ chromosome annotations with transcript selectors from cdot.
///
/// This function updates existing NC_ annotation files in the cache to include
/// all transcripts that map to each chromosome, enabling mutalyzer to normalize
/// intronic variants with genomic context.
pub fn enhance_nc_annotations<P: AsRef<Path>>(
    cache_dir: P,
    cdot: &CdotMapper,
) -> Result<usize, FerroError> {
    let cache_dir = cache_dir.as_ref();
    let mut enhanced = 0;

    // Find all NC_ annotation files
    let nc_files: Vec<_> = std::fs::read_dir(cache_dir)
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to read cache dir: {}", e),
        })?
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.file_name();
            let name = name.to_string_lossy();
            name.starts_with("NC_") && name.ends_with(".annotations")
        })
        .collect();

    eprintln!(
        "Enhancing {} NC_ annotation files with transcript selectors...",
        nc_files.len()
    );

    for entry in nc_files {
        let ann_path = entry.path();
        let filename = entry.file_name();
        let nc_id = filename
            .to_string_lossy()
            .trim_end_matches(".annotations")
            .to_string();

        // Get transcripts on this chromosome from cdot
        let tx_ids = cdot.transcripts_on_contig(&nc_id);
        if tx_ids.is_empty() {
            continue;
        }

        // Get sequence length from sequence file
        let seq_path = cache_dir.join(format!("{}.sequence", nc_id));
        let seq_len = if seq_path.exists() {
            std::fs::metadata(&seq_path)
                .map(|m| m.len() as usize)
                .unwrap_or(0)
        } else {
            continue;
        };

        // Collect transcript data
        let transcripts: Vec<_> = tx_ids
            .iter()
            .filter_map(|id| cdot.get_transcript(id).map(|tx| (*id, tx)))
            .collect();

        if transcripts.is_empty() {
            continue;
        }

        // Build enhanced annotations
        let annotations = build_nc_annotations_with_selectors(&nc_id, seq_len, &transcripts);

        // Write annotations
        let ann_file = File::create(&ann_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to create {}: {}", ann_path.display(), e),
        })?;
        serde_json::to_writer(ann_file, &annotations).map_err(|e| FerroError::Io {
            msg: format!("Failed to write annotations: {}", e),
        })?;

        enhanced += 1;
        if enhanced % 10 == 0 {
            eprintln!(
                "  Enhanced {} NC_ files ({} transcripts for {})...",
                enhanced,
                transcripts.len(),
                nc_id
            );
        }
    }

    eprintln!("Enhanced {} NC_ annotation files", enhanced);
    Ok(enhanced)
}

/// Extract all unique accessions from a ClinVar hgvs4variation file.
///
/// Handles gzipped files automatically. Extracts accessions from both
/// NucleotideExpression (column 6) and ProteinExpression (column 8).
pub fn extract_clinvar_accessions<P: AsRef<Path>>(
    clinvar_file: P,
) -> Result<Vec<String>, FerroError> {
    use flate2::read::GzDecoder;

    let clinvar_file = clinvar_file.as_ref();

    let file = File::open(clinvar_file).map_err(|e| FerroError::Io {
        msg: format!(
            "Failed to open ClinVar file {}: {}",
            clinvar_file.display(),
            e
        ),
    })?;

    let reader: Box<dyn BufRead> = if clinvar_file.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut all_accessions = HashSet::new();
    let mut line_count = 0usize;

    // Columns containing HGVS expressions (0-indexed)
    let nucleotide_column = 6;
    let protein_column = 8;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Error reading line: {}", e),
        })?;

        // Skip header lines
        if line.starts_with('#') {
            continue;
        }

        line_count += 1;
        if line_count.is_multiple_of(500_000) {
            eprintln!(
                "  Processed {} lines, found {} accessions...",
                line_count,
                all_accessions.len()
            );
        }

        let parts: Vec<&str> = line.split('\t').collect();

        // Extract from both nucleotide and protein columns
        for &column in &[nucleotide_column, protein_column] {
            if let Some(field) = parts.get(column) {
                for pattern in field.split('|') {
                    let pattern = pattern.trim();
                    if pattern.is_empty() || pattern == "-" {
                        continue;
                    }
                    // Extract all accessions from this pattern
                    for acc in extract_accessions_from_pattern(pattern) {
                        all_accessions.insert(acc);
                    }
                }
            }
        }
    }

    eprintln!("  Processed {} lines total", line_count);

    let mut result: Vec<String> = all_accessions.into_iter().collect();
    result.sort();
    Ok(result)
}

/// Fetch all LRG sequences from EBI that are not already in the cache.
///
/// If `ferro_reference` is provided, first copies LRG sequences from ferro's lrg/ directory.
/// Only fetches from EBI for sequences not in ferro or cache.
/// Returns the number of sequences fetched/copied.
pub fn fetch_lrg_sequences(
    cache_dir: &Path,
    ferro_reference: Option<&Path>,
) -> Result<usize, FerroError> {
    use std::process::Command;

    const LRG_MAX_ID: usize = 1400;
    const LRG_BASE_URL: &str = "https://ftp.ebi.ac.uk/pub/databases/lrgex/fasta/";

    let mut copied = 0usize;
    let mut fetched = 0usize;
    let mut skipped = 0usize;

    // First, try to copy from ferro reference if available
    let ferro_lrg_dir = ferro_reference.map(|r| r.join("lrg"));

    for i in 1..=LRG_MAX_ID {
        let lrg_id = format!("LRG_{}", i);
        let seq_path = cache_dir.join(format!("{}.sequence", lrg_id));

        // Skip if already exists in cache
        if seq_path.exists() {
            skipped += 1;
            continue;
        }

        // Try to copy from ferro reference first
        if let Some(ref lrg_dir) = ferro_lrg_dir {
            let ferro_fasta = lrg_dir.join(format!("{}.fasta", lrg_id));
            if ferro_fasta.exists() {
                // Read and parse ferro's FASTA file
                if let Ok(fasta_content) = std::fs::read_to_string(&ferro_fasta) {
                    let mut sequence = String::new();
                    for line in fasta_content.lines() {
                        if !line.starts_with('>') {
                            sequence.push_str(line.trim());
                        }
                    }

                    if !sequence.is_empty() {
                        // Write sequence file
                        if std::fs::write(&seq_path, &sequence).is_ok() {
                            // Write minimal annotations
                            let annotations = build_minimal_annotations(&lrg_id, sequence.len());
                            let annotations_path =
                                cache_dir.join(format!("{}.annotations", lrg_id));
                            if let Ok(json_str) = serde_json::to_string_pretty(&annotations) {
                                let _ = std::fs::write(&annotations_path, json_str);
                            }

                            copied += 1;
                            if copied.is_multiple_of(100) {
                                eprintln!("  Copied {} LRG sequences from ferro...", copied);
                            }
                            continue;
                        }
                    }
                }
            }
        }

        // Fall back to fetching from EBI
        let url = format!("{}LRG_{}.fasta", LRG_BASE_URL, i);

        // Fetch the FASTA file
        let output = Command::new("curl").args(["-s", "-f", "-L", &url]).output();

        let output = match output {
            Ok(o) if o.status.success() => o,
            _ => continue, // LRG doesn't exist or fetch failed
        };

        let fasta_content = String::from_utf8_lossy(&output.stdout);
        if fasta_content.is_empty() || !fasta_content.starts_with('>') {
            continue;
        }

        // Parse FASTA to extract sequence
        let mut sequence = String::new();
        for line in fasta_content.lines() {
            if !line.starts_with('>') {
                sequence.push_str(line.trim());
            }
        }

        if sequence.is_empty() {
            continue;
        }

        // Write sequence file
        if std::fs::write(&seq_path, &sequence).is_err() {
            continue;
        }

        // Write minimal annotations
        let annotations = build_minimal_annotations(&lrg_id, sequence.len());
        let annotations_path = cache_dir.join(format!("{}.annotations", lrg_id));
        if let Ok(json_str) = serde_json::to_string_pretty(&annotations) {
            let _ = std::fs::write(&annotations_path, json_str);
        }

        fetched += 1;
        if fetched.is_multiple_of(100) {
            eprintln!("  Fetched {} LRG sequences from EBI...", fetched);
        }

        // Rate limit to be nice to EBI servers
        std::thread::sleep(std::time::Duration::from_millis(100));
    }

    if skipped > 0 {
        eprintln!("  Skipped {} existing LRG sequences", skipped);
    }
    if copied > 0 {
        eprintln!("  Copied {} LRG sequences from ferro reference", copied);
    }
    if fetched > 0 {
        eprintln!("  Fetched {} LRG sequences from EBI", fetched);
    }

    Ok(copied + fetched)
}

/// Enhance LRG annotations using XML files with full exon/CDS structure.
///
/// If `ferro_reference` is provided, reads XML from ferro's lrg/ directory.
/// Only fetches from EBI for XML files not found in ferro.
///
/// This replaces the minimal annotations with proper mutalyzer-compatible annotations
/// that include transcript selectors (t1, t2, etc.) with exon and CDS coordinates.
pub fn enhance_lrg_annotations(
    cache_dir: &Path,
    ferro_reference: Option<&Path>,
) -> Result<usize, FerroError> {
    use std::process::Command;

    const LRG_MAX_ID: usize = 1400;
    const LRG_XML_URL: &str = "https://ftp.ebi.ac.uk/pub/databases/lrgex/";

    let mut enhanced = 0usize;
    let mut skipped = 0usize;
    let mut from_ferro = 0usize;
    let mut from_ebi = 0usize;
    let mut failed = 0usize;

    // Check ferro reference for cached XML files
    let ferro_lrg_dir = ferro_reference.map(|r| r.join("lrg"));

    eprintln!("Enhancing LRG annotations...");

    for i in 1..=LRG_MAX_ID {
        let lrg_id = format!("LRG_{}", i);
        let seq_path = cache_dir.join(format!("{}.sequence", &lrg_id));

        // Skip if no sequence file exists (LRG doesn't exist)
        if !seq_path.exists() {
            continue;
        }

        // Check if already has enhanced annotations (has "features" with transcript selectors)
        let ann_path = cache_dir.join(format!("{}.annotations", &lrg_id));
        if ann_path.exists() {
            if let Ok(content) = std::fs::read_to_string(&ann_path) {
                // Check if it has proper exon structure (more than just gene features)
                if content.contains("\"type\":\"exon\"") {
                    skipped += 1;
                    continue;
                }
            }
        }

        // Try to read XML from ferro reference first
        let xml_content = if let Some(ref lrg_dir) = ferro_lrg_dir {
            let ferro_xml_path = lrg_dir.join(format!("{}.xml", lrg_id));
            if ferro_xml_path.exists() {
                match std::fs::read_to_string(&ferro_xml_path) {
                    Ok(content) if content.contains("<lrg") => {
                        from_ferro += 1;
                        Some(content)
                    }
                    _ => None,
                }
            } else {
                None
            }
        } else {
            None
        };

        // Fall back to fetching from EBI if not in ferro
        let xml_content = match xml_content {
            Some(content) => content,
            None => {
                let url = format!("{}{}.xml", LRG_XML_URL, lrg_id);
                let output = Command::new("curl").args(["-s", "-f", "-L", &url]).output();

                match output {
                    Ok(o) if o.status.success() => {
                        let content = String::from_utf8_lossy(&o.stdout).to_string();
                        if content.is_empty() || !content.contains("<lrg") {
                            failed += 1;
                            continue;
                        }
                        from_ebi += 1;
                        // Rate limit only for EBI fetches
                        std::thread::sleep(std::time::Duration::from_millis(100));
                        content
                    }
                    _ => {
                        failed += 1;
                        continue;
                    }
                }
            }
        };

        // Parse XML and build annotations
        match parse_lrg_xml_to_annotations(&lrg_id, &xml_content) {
            Ok((genomic_ann, transcript_anns)) => {
                // Write genomic (LRG_N) annotations
                if let Ok(json_str) = serde_json::to_string_pretty(&genomic_ann) {
                    let _ = std::fs::write(&ann_path, json_str);
                }

                // Write transcript annotations and sequences
                for (tx_id, tx_ann, tx_seq) in transcript_anns {
                    let tx_ann_path = cache_dir.join(format!("{}.annotations", tx_id));
                    let tx_seq_path = cache_dir.join(format!("{}.sequence", tx_id));

                    if let Ok(json_str) = serde_json::to_string_pretty(&tx_ann) {
                        let _ = std::fs::write(&tx_ann_path, json_str);
                    }
                    if let Some(seq) = tx_seq {
                        let _ = std::fs::write(&tx_seq_path, seq);
                    }
                }

                enhanced += 1;
                if enhanced.is_multiple_of(50) {
                    eprintln!("  Enhanced {} LRG annotations...", enhanced);
                }
            }
            Err(_) => {
                failed += 1;
            }
        }
    }

    if skipped > 0 {
        eprintln!("  Skipped {} already-enhanced LRG annotations", skipped);
    }
    if from_ferro > 0 {
        eprintln!(
            "  Enhanced {} LRG annotations from ferro XML cache",
            from_ferro
        );
    }
    if from_ebi > 0 {
        eprintln!("  Enhanced {} LRG annotations from EBI", from_ebi);
    }
    if failed > 0 {
        eprintln!("  Failed to enhance {} LRG annotations", failed);
    }

    Ok(enhanced)
}

/// Parse LRG XML and build mutalyzer-compatible annotations.
///
/// Returns (genomic_annotations, Vec<(transcript_id, transcript_annotations, transcript_sequence)>)
#[allow(clippy::type_complexity)]
fn parse_lrg_xml_to_annotations(
    lrg_id: &str,
    xml_content: &str,
) -> Result<
    (
        serde_json::Value,
        Vec<(String, serde_json::Value, Option<String>)>,
    ),
    FerroError,
> {
    // Extract sequence length from the XML
    let seq_len = extract_xml_value(xml_content, "sequence")
        .map(|s| s.len())
        .unwrap_or(0);

    // Extract mol_type (should be "dna")
    let mol_type = extract_xml_attr(xml_content, "mol_type").unwrap_or_else(|| "dna".to_string());

    // Extract transcripts
    let mut transcript_features = Vec::new();
    let mut transcript_annotations = Vec::new();

    // Find all transcript blocks
    let mut pos = 0;
    while let Some(tx_start) = xml_content[pos..].find("<transcript name=\"") {
        let tx_start = pos + tx_start;
        let tx_name_start = tx_start + 18; // len of '<transcript name="'
        let tx_name_end = match xml_content[tx_name_start..].find('"') {
            Some(e) => tx_name_start + e,
            None => break,
        };
        let tx_name = &xml_content[tx_name_start..tx_name_end];
        let tx_id = format!("{}{}", lrg_id, tx_name);

        // Find the end of this transcript block
        let tx_end = match xml_content[tx_start..].find("</transcript>") {
            Some(e) => tx_start + e + 13,
            None => break,
        };
        let tx_block = &xml_content[tx_start..tx_end];

        // Extract transcript coordinates
        let (tx_start_coord, tx_end_coord, tx_strand) =
            extract_lrg_coordinates(tx_block, lrg_id).unwrap_or((0, seq_len, 1));

        // Extract cdna sequence
        let cdna_seq = extract_xml_value(tx_block, "sequence");

        // Extract CDS coordinates
        let (cds_start, cds_end) = extract_lrg_cds_coordinates(tx_block, lrg_id);

        // Extract exons
        let exons = extract_lrg_exons(tx_block, lrg_id);

        // Determine transcript type based on CDS presence
        let tx_mol_type = if cds_start.is_some() { "mRNA" } else { "ncRNA" };

        // Build transcript selector feature for genomic annotation
        let mut tx_feature = json!({
            "type": tx_mol_type,
            "id": tx_name,
            "location": {
                "type": "range",
                "start": {"type": "point", "position": tx_start_coord},
                "end": {"type": "point", "position": tx_end_coord},
                "strand": tx_strand
            },
            "features": []
        });

        // Add exons to transcript feature
        let exon_features: Vec<serde_json::Value> = exons
            .iter()
            .map(|(label, start, end)| {
                json!({
                    "type": "exon",
                    "id": label,
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": *start},
                        "end": {"type": "point", "position": *end},
                        "strand": tx_strand
                    }
                })
            })
            .collect();
        tx_feature["features"] = json!(exon_features);

        // Add CDS if present
        if let (Some(cds_s), Some(cds_e)) = (cds_start, cds_end) {
            tx_feature["cds"] = json!({
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": cds_s},
                    "end": {"type": "point", "position": cds_e}
                }
            });
        }

        transcript_features.push(tx_feature);

        // Build standalone transcript annotation
        let cdna_len = cdna_seq
            .as_ref()
            .map(|s| s.len())
            .unwrap_or(tx_end_coord - tx_start_coord);

        // Calculate CDS in transcript coordinates
        let tx_cds = if let (Some(cds_s), Some(cds_e)) = (cds_start, cds_end) {
            // Convert genomic CDS to transcript coordinates using exons
            calculate_transcript_cds(&exons, cds_s, cds_e, tx_start_coord)
        } else {
            None
        };

        let mut tx_ann = json!({
            "type": "record",
            "id": tx_id,
            "location": {
                "type": "range",
                "start": {"type": "point", "position": 0},
                "end": {"type": "point", "position": cdna_len}
            },
            "qualifiers": {
                "mol_type": tx_mol_type,
                "name": tx_id
            },
            "features": [{
                "type": "gene",
                "id": tx_id,
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 0},
                    "end": {"type": "point", "position": cdna_len},
                    "strand": 1
                },
                "features": [{
                    "type": tx_mol_type,
                    "id": tx_id,
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": cdna_len}
                    },
                    "features": []
                }]
            }]
        });

        // Add CDS to transcript annotation if present
        if let Some((cds_tx_start, cds_tx_end)) = tx_cds {
            if let Some(gene_features) = tx_ann["features"][0]["features"].as_array_mut() {
                if let Some(mrna) = gene_features.first_mut() {
                    mrna["features"] = json!([{
                        "type": "CDS",
                        "id": tx_id,
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": cds_tx_start},
                            "end": {"type": "point", "position": cds_tx_end}
                        }
                    }]);
                }
            }
        }

        transcript_annotations.push((tx_id, tx_ann, cdna_seq));

        pos = tx_end;
    }

    // Build genomic annotation
    let genomic_ann = json!({
        "type": "record",
        "id": lrg_id,
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": seq_len}
        },
        "qualifiers": {
            "mol_type": mol_type,
            "name": lrg_id
        },
        "features": [{
            "type": "gene",
            "id": lrg_id,
            "location": {
                "type": "range",
                "start": {"type": "point", "position": 0},
                "end": {"type": "point", "position": seq_len},
                "strand": 1
            },
            "features": transcript_features
        }]
    });

    Ok((genomic_ann, transcript_annotations))
}

/// Extract a simple XML element value.
fn extract_xml_value(xml: &str, tag: &str) -> Option<String> {
    let start_tag = format!("<{}>", tag);
    let end_tag = format!("</{}>", tag);

    let start = xml.find(&start_tag)?;
    let content_start = start + start_tag.len();
    let end = xml[content_start..].find(&end_tag)?;

    Some(xml[content_start..content_start + end].to_string())
}

/// Extract a simple XML attribute value from content containing the attribute.
fn extract_xml_attr(xml: &str, tag: &str) -> Option<String> {
    let start_tag = format!("<{}>", tag);
    let end_tag = format!("</{}>", tag);

    let start = xml.find(&start_tag)?;
    let content_start = start + start_tag.len();
    let end = xml[content_start..].find(&end_tag)?;

    Some(xml[content_start..content_start + end].to_string())
}

/// Extract LRG coordinates from a block.
fn extract_lrg_coordinates(block: &str, coord_system: &str) -> Option<(usize, usize, i32)> {
    // Look for coordinates with matching coord_system
    let search = format!("coord_system=\"{}\"", coord_system);
    let pos = block.find(&search)?;
    let line_start = block[..pos].rfind('<')?;
    let line_end = block[pos..].find('>')?;
    let coord_line = &block[line_start..pos + line_end + 1];

    let start = extract_number_attr(coord_line, "start")?;
    let end = extract_number_attr(coord_line, "end")?;
    let strand = extract_number_attr(coord_line, "strand").unwrap_or(1) as i32;

    Some((start, end, strand))
}

/// Extract CDS coordinates from a transcript block.
fn extract_lrg_cds_coordinates(block: &str, coord_system: &str) -> (Option<usize>, Option<usize>) {
    let search = "<coding_region>";
    if let Some(cr_start) = block.find(search) {
        let cr_end = block[cr_start..]
            .find("</coding_region>")
            .unwrap_or(block.len() - cr_start);
        let cr_block = &block[cr_start..cr_start + cr_end];

        if let Some((start, end, _)) = extract_lrg_coordinates(cr_block, coord_system) {
            return (Some(start), Some(end));
        }
    }
    (None, None)
}

/// Extract exons from a transcript block.
fn extract_lrg_exons(block: &str, coord_system: &str) -> Vec<(String, usize, usize)> {
    let mut exons = Vec::new();
    let mut pos = 0;

    while let Some(exon_start) = block[pos..].find("<exon label=\"") {
        let exon_start = pos + exon_start;
        let label_start = exon_start + 13;
        let label_end = match block[label_start..].find('"') {
            Some(e) => label_start + e,
            None => break,
        };
        let label = block[label_start..label_end].to_string();

        // Find the exon block end
        let exon_end = match block[exon_start..].find("</exon>") {
            Some(e) => exon_start + e + 7,
            None => match block[exon_start..].find("/>") {
                Some(e) => exon_start + e + 2,
                None => break,
            },
        };
        let exon_block = &block[exon_start..exon_end];

        // Extract coordinates for the specified coord_system
        if let Some((start, end, _)) = extract_lrg_coordinates(exon_block, coord_system) {
            exons.push((label, start, end));
        }

        pos = exon_end;
    }

    exons
}

/// Extract a numeric attribute value from an XML element.
fn extract_number_attr(element: &str, attr: &str) -> Option<usize> {
    let search = format!("{}=\"", attr);
    let pos = element.find(&search)?;
    let start = pos + search.len();
    let end = element[start..].find('"')?;
    element[start..start + end].parse().ok()
}

/// Calculate CDS coordinates in transcript space from genomic coordinates.
fn calculate_transcript_cds(
    exons: &[(String, usize, usize)],
    cds_start: usize,
    cds_end: usize,
    _tx_start: usize,
) -> Option<(usize, usize)> {
    if exons.is_empty() {
        return None;
    }

    let mut tx_pos = 0usize;
    let mut cds_tx_start = None;
    let mut cds_tx_end = None;

    for (_label, exon_start, exon_end) in exons {
        let exon_len = exon_end - exon_start;

        // Check if CDS start is in this exon
        if cds_tx_start.is_none() && cds_start >= *exon_start && cds_start <= *exon_end {
            cds_tx_start = Some(tx_pos + (cds_start - exon_start));
        }

        // Check if CDS end is in this exon
        if cds_tx_end.is_none() && cds_end >= *exon_start && cds_end <= *exon_end {
            cds_tx_end = Some(tx_pos + (cds_end - exon_start));
        }

        tx_pos += exon_len;

        if cds_tx_start.is_some() && cds_tx_end.is_some() {
            break;
        }
    }

    match (cds_tx_start, cds_tx_end) {
        (Some(s), Some(e)) => Some((s, e)),
        _ => None,
    }
}

/// Statistics from building transcript→chromosome mapping
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MappingStats {
    pub from_cdot: usize,
    pub from_gene2refseq: usize,
    pub total: usize,
}

/// Transcript→chromosome mapping file format
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TranscriptMapping {
    pub assembly: String,
    pub mappings: HashMap<String, String>,
    pub stats: MappingStats,
}

/// Extract transcript→chromosome mappings from cdot JSON file
pub fn extract_cdot_mappings<P: AsRef<Path>>(
    cdot_path: P,
    assembly: &str,
) -> Result<HashMap<String, String>, FerroError> {
    let cdot_path = cdot_path.as_ref();
    eprintln!("Loading cdot mappings from {}...", cdot_path.display());

    // Load cdot JSON (supports .gz)
    let file = File::open(cdot_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open cdot file: {}", e),
    })?;

    let value: serde_json::Value = if cdot_path.extension().is_some_and(|e| e == "gz") {
        let decoder = flate2::read::GzDecoder::new(file);
        let reader = BufReader::new(decoder);
        serde_json::from_reader(reader).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse cdot JSON: {}", e),
        })?
    } else {
        let reader = BufReader::new(file);
        serde_json::from_reader(reader).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse cdot JSON: {}", e),
        })?
    };

    // Extract mappings from transcripts
    let mut mappings = HashMap::new();

    if let Some(transcripts) = value.get("transcripts").and_then(|t| t.as_object()) {
        for (acc, tx_data) in transcripts {
            // Navigate to genome_builds -> assembly -> contig
            if let Some(contig) = tx_data
                .get("genome_builds")
                .and_then(|gb| gb.get(assembly))
                .and_then(|asm| asm.get("contig"))
                .and_then(|c| c.as_str())
            {
                mappings.insert(acc.clone(), contig.to_string());
            }
        }
    }

    eprintln!(
        "  Extracted {} transcript→chromosome mappings",
        mappings.len()
    );
    Ok(mappings)
}

/// Download NCBI gene2refseq.gz file if not present
pub fn download_gene2refseq<P: AsRef<Path>>(output_path: P) -> Result<(), FerroError> {
    let output_path = output_path.as_ref();

    if output_path.exists() {
        eprintln!("gene2refseq already exists at {}", output_path.display());
        return Ok(());
    }

    eprintln!("Downloading gene2refseq.gz from NCBI (this may take a while, ~2GB)...");

    let url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz";

    // Use curl for the download
    let status = std::process::Command::new("curl")
        .args(["-o", &output_path.display().to_string(), "-L", "-f", url])
        .status()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to run curl: {}", e),
        })?;

    if !status.success() {
        return Err(FerroError::Io {
            msg: format!(
                "Failed to download gene2refseq.gz (exit code: {:?})",
                status.code()
            ),
        });
    }

    eprintln!("  Downloaded to {}", output_path.display());
    Ok(())
}

/// Extract transcript→chromosome mappings from NCBI gene2refseq.gz
///
/// Format: tab-separated with columns:
/// 0: tax_id (filter for 9606 = human)
/// 3: RNA_nucleotide_accession.version
/// 7: genomic_nucleotide_accession.version
/// 12: assembly (filter for target assembly)
pub fn extract_gene2refseq_mappings<P: AsRef<Path>>(
    gene2refseq_path: P,
    assembly: &str,
) -> Result<HashMap<String, String>, FerroError> {
    let gene2refseq_path = gene2refseq_path.as_ref();
    eprintln!("Extracting {} mappings from gene2refseq...", assembly);

    let file = File::open(gene2refseq_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open gene2refseq: {}", e),
    })?;

    let reader: Box<dyn BufRead> = if gene2refseq_path.extension().is_some_and(|e| e == "gz") {
        Box::new(BufReader::new(flate2::read::GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut mappings = HashMap::new();
    let mut lines_processed = 0u64;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Error reading gene2refseq: {}", e),
        })?;

        // Skip header
        if line.starts_with('#') {
            continue;
        }

        lines_processed += 1;
        if lines_processed.is_multiple_of(5_000_000) {
            eprintln!(
                "  Processed {} million lines...",
                lines_processed / 1_000_000
            );
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 13 {
            continue;
        }

        // Filter for human (tax_id = 9606)
        if cols[0] != "9606" {
            continue;
        }

        // Get transcript accession (column 3)
        let transcript = cols[3];
        if transcript == "-" || (!transcript.starts_with("NM_") && !transcript.starts_with("NR_")) {
            continue;
        }

        // Get genomic accession (column 7)
        let genomic = cols[7];
        if genomic == "-" || !genomic.starts_with("NC_") {
            continue;
        }

        // Filter for target assembly (column 12)
        let asm = cols[12];
        if !asm.contains(assembly) {
            continue;
        }

        // Add mapping (only if not already present - cdot takes priority)
        mappings
            .entry(transcript.to_string())
            .or_insert_with(|| genomic.to_string());
    }

    eprintln!(
        "  Extracted {} transcript→chromosome mappings from gene2refseq",
        mappings.len()
    );
    Ok(mappings)
}

/// Build combined transcript→chromosome mapping file
///
/// Uses cdot as primary source, gene2refseq as fallback for older transcripts.
pub fn build_transcript_chromosome_mapping<P: AsRef<Path>>(
    cdot_path: P,
    gene2refseq_path: Option<P>,
    output_path: P,
    assembly: &str,
) -> Result<MappingStats, FerroError> {
    let output_path = output_path.as_ref();

    // Extract mappings from cdot (primary source)
    let mut mappings = extract_cdot_mappings(&cdot_path, assembly)?;
    let from_cdot = mappings.len();

    // Add mappings from gene2refseq (fallback for older transcripts)
    let from_gene2refseq = if let Some(g2r_path) = gene2refseq_path {
        let g2r_mappings = extract_gene2refseq_mappings(g2r_path, assembly)?;
        let mut added = 0usize;
        for (transcript, chromosome) in g2r_mappings {
            if let std::collections::hash_map::Entry::Vacant(e) = mappings.entry(transcript) {
                e.insert(chromosome);
                added += 1;
            }
        }
        eprintln!("  Added {} mappings from gene2refseq (not in cdot)", added);
        added
    } else {
        0
    };

    let stats = MappingStats {
        from_cdot,
        from_gene2refseq,
        total: mappings.len(),
    };

    // Build output structure
    let output = TranscriptMapping {
        assembly: assembly.to_string(),
        mappings,
        stats: stats.clone(),
    };

    // Write to file
    let file = File::create(output_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create mapping file: {}", e),
    })?;

    serde_json::to_writer_pretty(file, &output).map_err(|e| FerroError::Json {
        msg: format!("Failed to write mapping file: {}", e),
    })?;

    eprintln!(
        "\nMapping file written to {} ({} mappings)",
        output_path.display(),
        stats.total
    );

    Ok(stats)
}

/// Load transcript→chromosome mapping from JSON file
pub fn load_transcript_mapping<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, String>, FerroError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open mapping file: {}", e),
    })?;

    let mapping: TranscriptMapping =
        serde_json::from_reader(BufReader::new(file)).map_err(|e| FerroError::Json {
            msg: format!("Failed to parse mapping file: {}", e),
        })?;

    Ok(mapping.mappings)
}

// ============================================================================
// Protein sequence caching for protein variant validation
// ============================================================================

/// Extract protein accessions from ClinVar hgvs4variation.txt.gz (column 9)
///
/// Returns a list of unique protein accessions found in the file.
/// Handles NP_, XP_, YP_ (RefSeq), UniProt (e.g., P04181), and LRG protein accessions.
pub fn extract_clinvar_protein_accessions<P: AsRef<Path>>(
    clinvar_path: P,
) -> Result<Vec<String>, FerroError> {
    use flate2::read::GzDecoder;

    let clinvar_path = clinvar_path.as_ref();
    eprintln!(
        "Extracting protein accessions from {}...",
        clinvar_path.display()
    );

    let file = File::open(clinvar_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open ClinVar file: {}", e),
    })?;

    let reader: Box<dyn BufRead> = if clinvar_path.to_string_lossy().ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut accessions = HashSet::new();
    let mut line_count = 0u64;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        // Skip headers
        if line.starts_with('#') {
            continue;
        }

        line_count += 1;
        if line_count.is_multiple_of(1_000_000) {
            eprintln!(
                "  Processed {} lines, found {} accessions...",
                line_count,
                accessions.len()
            );
        }

        // Column 9 (0-indexed: 8) contains protein HGVS
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() > 8 {
            let protein_col = fields[8];
            if protein_col != "-" && !protein_col.is_empty() {
                // Extract accession from HGVS like "NP_000079.2:p.Val600Glu"
                if let Some(acc) = protein_col.split(':').next() {
                    accessions.insert(acc.to_string());
                }
            }
        }
    }

    eprintln!(
        "  Found {} unique protein accessions from {} lines",
        accessions.len(),
        line_count
    );

    let mut result: Vec<String> = accessions.into_iter().collect();
    result.sort();
    Ok(result)
}

/// Statistics for protein cache population
#[derive(Debug, Clone, Default)]
pub struct ProteinCacheStats {
    /// Number of NP_ accessions fetched from NCBI
    pub ncbi_fetched: usize,
    /// Number of UniProt accessions fetched
    pub uniprot_fetched: usize,
    /// Number of accessions already cached
    pub already_cached: usize,
    /// Number of accessions that failed to fetch
    pub failed: usize,
    /// Total accessions processed
    pub total: usize,
}

/// Fetch protein sequences and populate the protein cache
///
/// Fetches sequences from NCBI (for NP_, XP_, YP_ accessions) and UniProt
/// (for UniProt accessions like P04181).
pub fn populate_protein_cache<P: AsRef<Path>>(
    accessions: &[String],
    cache_dir: P,
) -> Result<ProteinCacheStats, FerroError> {
    let cache_dir = cache_dir.as_ref();

    // Create cache directory
    std::fs::create_dir_all(cache_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create protein cache dir: {}", e),
    })?;

    let mut stats = ProteinCacheStats {
        total: accessions.len(),
        ..Default::default()
    };

    // Categorize accessions
    let mut ncbi_accessions: Vec<&str> = Vec::new();
    let mut uniprot_accessions: Vec<&str> = Vec::new();

    // Check for mutalyzer-format cache in parent directory
    let parent_dir = cache_dir.parent();

    for acc in accessions {
        // Check if mutalyzer-compatible .sequence file exists in parent dir
        // (that's where mutalyzer_retriever looks for cached entries)
        let already_cached = parent_dir
            .map(|p| p.join(format!("{}.sequence", acc)).exists())
            .unwrap_or(false);

        if already_cached {
            stats.already_cached += 1;
            continue;
        }

        // Categorize by accession type
        if acc.starts_with("NP_") || acc.starts_with("XP_") || acc.starts_with("YP_") {
            ncbi_accessions.push(acc);
        } else if is_uniprot_accession(acc) {
            uniprot_accessions.push(acc);
        }
        // Skip LRG protein accessions for now - they need special handling
    }

    eprintln!(
        "Protein cache population: {} total, {} already cached",
        stats.total, stats.already_cached
    );
    eprintln!(
        "  To fetch: {} NCBI, {} UniProt",
        ncbi_accessions.len(),
        uniprot_accessions.len()
    );

    // Fetch NCBI protein sequences in batches
    if !ncbi_accessions.is_empty() {
        eprintln!(
            "Fetching {} protein sequences from NCBI...",
            ncbi_accessions.len()
        );
        let fetched = fetch_ncbi_protein_batch(&ncbi_accessions, cache_dir)?;
        stats.ncbi_fetched = fetched;
        eprintln!("  Fetched {} NCBI protein sequences", fetched);
    }

    // Fetch UniProt sequences
    if !uniprot_accessions.is_empty() {
        eprintln!(
            "Fetching {} protein sequences from UniProt...",
            uniprot_accessions.len()
        );
        let fetched = fetch_uniprot_batch(&uniprot_accessions, cache_dir)?;
        stats.uniprot_fetched = fetched;
        eprintln!("  Fetched {} UniProt protein sequences", fetched);
    }

    stats.failed = stats.total - stats.already_cached - stats.ncbi_fetched - stats.uniprot_fetched;

    Ok(stats)
}

/// Check if an accession looks like a UniProt accession
fn is_uniprot_accession(acc: &str) -> bool {
    // UniProt accessions are 6-10 alphanumeric characters
    // Primary format: [OPQ][0-9][A-Z0-9]{3}[0-9] (e.g., P04181)
    // Secondary format: [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9] (e.g., A0A024R1R8)
    if acc.len() < 6 || acc.len() > 10 {
        return false;
    }
    // Must not have a version suffix like .1
    if acc.contains('.') {
        return false;
    }
    // Must be alphanumeric
    acc.chars().all(|c| c.is_ascii_alphanumeric())
}

/// Fetch a batch of protein sequences from NCBI
fn fetch_ncbi_protein_batch(accessions: &[&str], cache_dir: &Path) -> Result<usize, FerroError> {
    use std::process::Command;

    const BATCH_SIZE: usize = 100;
    let mut fetched = 0;

    for (batch_idx, batch) in accessions.chunks(BATCH_SIZE).enumerate() {
        let id_list = batch.join(",");

        // Use efetch to get FASTA format
        let url = build_efetch_url("protein", &id_list, "fasta", "text");

        let output = Command::new("curl")
            .args(["-s", "-L", &url])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to fetch protein batch: {}", e),
            })?;

        if !output.status.success() {
            eprintln!(
                "  Warning: Failed to fetch batch {}: {}",
                batch_idx + 1,
                String::from_utf8_lossy(&output.stderr)
            );
            continue;
        }

        let response = String::from_utf8_lossy(&output.stdout);

        // Parse FASTA response and save individual files
        let mut current_acc: Option<String> = None;
        let mut current_seq = String::new();

        for line in response.lines() {
            if let Some(header) = line.strip_prefix('>') {
                // Save previous sequence if any
                if let Some(ref acc) = current_acc {
                    if !current_seq.is_empty() {
                        save_protein_fasta(cache_dir, acc, &current_seq)?;
                        fetched += 1;
                    }
                }
                // Parse new accession from header like ">NP_000079.2 ..."
                current_acc = header.split_whitespace().next().map(|s| s.to_string());
                current_seq.clear();
            } else if current_acc.is_some() {
                current_seq.push_str(line.trim());
            }
        }

        // Save last sequence
        if let Some(ref acc) = current_acc {
            if !current_seq.is_empty() {
                save_protein_fasta(cache_dir, acc, &current_seq)?;
                fetched += 1;
            }
        }

        eprintln!(
            "  Batch {}/{}: fetched sequences",
            batch_idx + 1,
            accessions.len().div_ceil(BATCH_SIZE)
        );

        // Rate limit between batches
        std::thread::sleep(std::time::Duration::from_millis(350));
    }

    Ok(fetched)
}

/// Fetch a batch of protein sequences from UniProt using batch API
fn fetch_uniprot_batch(accessions: &[&str], cache_dir: &Path) -> Result<usize, FerroError> {
    use std::process::Command;

    const BATCH_SIZE: usize = 100; // UniProt allows up to 100 accessions per request
    let mut fetched = 0;

    for (batch_idx, batch) in accessions.chunks(BATCH_SIZE).enumerate() {
        // Build comma-separated accession list for batch API
        let acc_list = batch.join(",");

        // UniProt batch FASTA URL
        let url = format!(
            "https://rest.uniprot.org/uniprotkb/accessions?accessions={}&format=fasta",
            acc_list
        );

        let output = Command::new("curl")
            .args(["-s", "-L", &url])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to fetch UniProt batch: {}", e),
            })?;

        if !output.status.success() {
            eprintln!(
                "  Warning: Failed to fetch UniProt batch {}: {}",
                batch_idx + 1,
                String::from_utf8_lossy(&output.stderr)
            );
            continue;
        }

        let response = String::from_utf8_lossy(&output.stdout);

        // Parse multi-FASTA response and save individual files
        let mut current_acc: Option<String> = None;
        let mut current_seq = String::new();

        for line in response.lines() {
            if let Some(header) = line.strip_prefix('>') {
                // Save previous sequence if any
                if let Some(ref acc) = current_acc {
                    if !current_seq.is_empty() {
                        save_protein_fasta(cache_dir, acc, &current_seq)?;
                        fetched += 1;
                    }
                }
                // Parse accession from header like ">sp|P04181|OATC_HUMAN ..."
                // or ">tr|A0A024R1R8|A0A024R1R8_HUMAN ..."
                current_acc = header
                    .split('|')
                    .nth(1)
                    .or_else(|| header.split_whitespace().next())
                    .map(|s| s.to_string());
                current_seq.clear();
            } else if current_acc.is_some() {
                current_seq.push_str(line.trim());
            }
        }

        // Save last sequence
        if let Some(ref acc) = current_acc {
            if !current_seq.is_empty() {
                save_protein_fasta(cache_dir, acc, &current_seq)?;
                fetched += 1;
            }
        }

        eprintln!(
            "  Batch {}/{}: fetched sequences",
            batch_idx + 1,
            accessions.len().div_ceil(BATCH_SIZE)
        );

        // Rate limit between batches
        std::thread::sleep(std::time::Duration::from_millis(200));
    }

    Ok(fetched)
}

/// Save a protein sequence to the cache in FASTA format
fn save_protein_fasta(cache_dir: &Path, accession: &str, sequence: &str) -> Result<(), FerroError> {
    use std::io::Write;

    // Sanitize accession for filesystem
    let safe_name = accession.replace(['/', '\\', ':', '*', '?', '"', '<', '>', '|'], "_");
    let path = cache_dir.join(format!("{}.fa", safe_name));

    let mut file = File::create(&path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create protein cache file: {}", e),
    })?;

    // Write FASTA format
    writeln!(file, ">{}", accession).map_err(|e| FerroError::Io {
        msg: format!("Failed to write FASTA header: {}", e),
    })?;

    // Write sequence in 60-char lines
    for chunk in sequence.as_bytes().chunks(60) {
        writeln!(file, "{}", std::str::from_utf8(chunk).unwrap_or("")).map_err(|e| {
            FerroError::Io {
                msg: format!("Failed to write sequence: {}", e),
            }
        })?;
    }

    // Also create mutalyzer-compatible .sequence and .annotations files
    // in the parent directory (where mutalyzer_retriever looks for cached entries)
    if let Some(parent) = cache_dir.parent() {
        let seq_path = parent.join(format!("{}.sequence", accession));
        let ann_path = parent.join(format!("{}.annotations", accession));

        // Only create if not already present (don't overwrite existing cache)
        if !seq_path.exists() {
            std::fs::write(&seq_path, sequence).map_err(|e| FerroError::Io {
                msg: format!("Failed to write protein sequence cache: {}", e),
            })?;

            let annotations = build_protein_annotations(accession, sequence.len());
            std::fs::write(
                &ann_path,
                serde_json::to_string(&annotations).unwrap_or_default(),
            )
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to write protein annotations: {}", e),
            })?;
        }
    }

    Ok(())
}
