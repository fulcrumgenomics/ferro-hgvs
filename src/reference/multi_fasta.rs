//! Multi-FASTA reference sequence provider
//!
//! This module provides a reference provider that can load sequences from
//! multiple FASTA files, useful for large reference datasets split across files.
//!
//! # Coordinate Systems
//!
//! This module handles coordinate conversions when loading sequences:
//!
//! | Source | Basis | Target | Notes |
//! |--------|-------|--------|-------|
//! | FASTA/FAI index | 0-based | Internal | Half-open `[start, end)` |
//! | cdot genomic | 0-based | Transcript | `genome_start` is 0-based, converted to 1-based |
//! | cdot tx | 1-based | Transcript | `tx_start`/`tx_end` used directly (already 1-based) |
//! | cdot CDS | 0-based | Transcript | `cds_start` converted via `+ 1` to 1-based |
//!
//! **Important**: The `get_sequence()` method uses 0-based half-open coordinates
//! consistent with FASTA index conventions.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use log::warn;
use lru::LruCache;
use rustc_hash::FxHashMap;

use crate::data::cdot::CdotMapper;
use crate::error::FerroError;
use crate::reference::authoritative::CanonicalOverrides;
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;

/// Index entry for a sequence in a FASTA file
#[derive(Debug, Clone)]
struct FastaIndexEntry {
    /// Path to the FASTA file containing this sequence
    file_path: PathBuf,
    /// Sequence name/accession
    name: String,
    /// Length of the sequence
    length: u64,
    /// Byte offset to the start of sequence data
    offset: u64,
    /// Number of bases per line
    line_bases: u64,
    /// Number of bytes per line (including newline)
    line_bytes: u64,
}

/// Supplemental transcript info for a single transcript.
#[derive(Debug, Clone)]
pub struct SupplementalTranscriptInfo {
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
    pub gene_symbol: Option<String>,
    pub sequence_length: u64,
}

/// Supplemental transcript metadata for transcripts not in cdot.
#[derive(Debug, Clone, Default)]
pub struct SupplementalCdsInfo {
    /// Mapping from accession to transcript info.
    /// CDS coordinates are 1-based inclusive.
    pub transcripts: FxHashMap<String, SupplementalTranscriptInfo>,
}

/// Parse the integer version suffix of an accession (`"NM_003002.4"` -> `Some(4)`),
/// or `None` when there is no trailing numeric `.<n>`.
fn accession_version(accession: &str) -> Option<u32> {
    accession
        .rsplit_once('.')
        .and_then(|(_, v)| v.parse::<u32>().ok())
}

/// Build the `base -> versioned` fallback map from a fully populated FASTA
/// index, deterministically choosing the **highest** version per base
/// accession.
///
/// Building this by iterating the index in hash order with last-writer-wins
/// made version fallback nondeterministic across runs (the index uses a
/// fixed-seed hasher now, but relying on iteration order for *which* version
/// wins is fragile regardless). Keying on `(version, accession)` gives a stable
/// total order: the newest version always wins, independent of insertion or
/// iteration order.
fn build_base_to_versioned(
    index: &FxHashMap<String, FastaIndexEntry>,
) -> FxHashMap<String, String> {
    let mut map: FxHashMap<String, String> = FxHashMap::default();
    for name in index.keys() {
        let Some(base) = name.split('.').next() else {
            continue;
        };
        let replace = match map.get(base) {
            None => true,
            Some(existing) => {
                (accession_version(name), name.as_str())
                    > (accession_version(existing), existing.as_str())
            }
        };
        if replace {
            map.insert(base.to_string(), name.clone());
        }
    }
    map
}

/// Multi-FASTA reference provider
///
/// Provides efficient random access to sequences across multiple FASTA files.
/// Useful for reference datasets split across files (e.g., RefSeq transcripts).
///
/// # Example
///
/// ```ignore
/// use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
///
/// let provider = MultiFastaProvider::from_directory("reference_data/transcripts")?;
/// let seq = provider.get_sequence("NM_000546.6", 0, 100)?;
/// ```
pub struct MultiFastaProvider {
    /// Index mapping sequence names to their locations.
    ///
    /// These index maps use [`FxHashMap`] (a fast, fixed-seed hasher) rather
    /// than the default SipHash: the keys are internal reference accessions,
    /// not attacker-controlled, so DoS resistance is unnecessary, and building
    /// the ~400k+ entry indexes is a large slice of startup. The fixed seed
    /// also makes iteration deterministic across runs.
    index: FxHashMap<String, FastaIndexEntry>,
    /// Index mapping base accession (without version) to versioned accession
    base_to_versioned: FxHashMap<String, String>,
    /// Chromosome aliases (e.g., NC_000001 -> chr1)
    aliases: FxHashMap<String, String>,
    /// Optional cdot transcript metadata (CDS positions, exon coordinates)
    cdot_mapper: Option<CdotMapper>,
    /// Supplemental CDS info for transcripts not in cdot
    supplemental_cds: SupplementalCdsInfo,
    /// Authoritative canonical overrides (issue #520), loaded from the manifest.
    canonical_overrides: CanonicalOverrides,
    /// Protein FASTA index (NP_/XP_ accessions) for get_protein_sequence (#520).
    protein_index: FxHashMap<String, FastaIndexEntry>,
    /// Bounded LRU memoizing resolved transcripts, keyed on the exact resolution
    /// inputs `(id, build_hint)`. Resolving a transcript materializes its full
    /// sequence from the FASTA and rebuilds cdot exon/CDS metadata — by far the
    /// dominant cost in batch normalization. Most real workloads (transcript- or
    /// genome-position-sorted input) revisit the same transcript many times in a
    /// row, so memoizing turns those repeats into cheap `Arc` clones. The result
    /// is byte-identical to recomputing it, so this is purely a timing change.
    /// Keying on `(id, build_hint)` (not the bare accession) keeps multi-build
    /// resolution correct: the same accession on different builds caches
    /// distinctly. See the resolve-cache perf work.
    transcript_cache: TranscriptCache,
}

/// Resolution inputs that uniquely identify a resolved transcript: the requested
/// accession and the optional genome-build hint. Used as the resolve-cache key.
type TranscriptCacheKey = (String, Option<String>);

/// Bounded LRU memoizing resolved transcripts (see [`MultiFastaProvider`]'s
/// `transcript_cache` field). Wrapped in a [`Mutex`] for `&self` access from the
/// [`ReferenceProvider`] trait, which is shared across threads.
type TranscriptCache = Mutex<LruCache<TranscriptCacheKey, Arc<Transcript>>>;

/// Default capacity (number of distinct transcripts) for the resolve cache.
///
/// Sized well above the working set of typical batch runs (ClinVar-scale inputs
/// touch on the order of tens of thousands of distinct transcripts), so locality
/// is captured without eviction churn, while still bounding memory for
/// long-running services that would otherwise accumulate every transcript ever
/// requested.
const TRANSCRIPT_CACHE_CAPACITY: usize = 65_536;

/// Build an empty resolve cache at the default capacity.
fn new_transcript_cache() -> TranscriptCache {
    let capacity =
        NonZeroUsize::new(TRANSCRIPT_CACHE_CAPACITY).expect("cache capacity is non-zero");
    Mutex::new(LruCache::new(capacity))
}

impl MultiFastaProvider {
    /// Create a new multi-FASTA provider from a directory of FASTA files
    ///
    /// Scans the directory for .fna and .fa files with accompanying .fai indexes.
    pub fn from_directory<P: AsRef<Path>>(dir: P) -> Result<Self, FerroError> {
        let dir = dir.as_ref();

        if !dir.exists() {
            return Err(FerroError::Io {
                msg: format!("Reference directory not found: {}", dir.display()),
            });
        }

        let mut index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();

        // Find all FASTA files with indexes
        let entries = std::fs::read_dir(dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to read directory {}: {}", dir.display(), e),
        })?;

        for entry in entries {
            let entry = entry.map_err(|e| FerroError::Io {
                msg: format!("Failed to read directory entry: {}", e),
            })?;

            let path = entry.path();
            let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");

            // Look for .fna, .fa, or .fasta files
            if ext == "fna" || ext == "fa" || ext == "fasta" {
                let fai_path = PathBuf::from(format!("{}.fai", path.display()));

                if fai_path.exists() {
                    // Load index for this file
                    let file_index = load_fai_index(&path, &fai_path)?;

                    for (name, entry) in file_index {
                        index.insert(name, entry);
                    }
                }
            }
        }

        if index.is_empty() {
            return Err(FerroError::Io {
                msg: format!(
                    "No indexed FASTA files found in {}. Run 'ferro-benchmark prepare' first.",
                    dir.display()
                ),
            });
        }

        eprintln!(
            "Loaded {} sequences from {} FASTA files",
            index.len(),
            count_unique_files(&index)
        );

        Ok(Self {
            base_to_versioned: build_base_to_versioned(&index),
            index,
            aliases: build_chromosome_aliases(),
            cdot_mapper: None,
            supplemental_cds: SupplementalCdsInfo::default(),
            canonical_overrides: CanonicalOverrides::default(),
            protein_index: FxHashMap::default(),
            transcript_cache: new_transcript_cache(),
        })
    }

    /// Create a provider from multiple directories (e.g., transcripts + genome)
    pub fn from_directories<P: AsRef<Path>>(dirs: &[P]) -> Result<Self, FerroError> {
        let mut combined_index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();

        for dir in dirs {
            let dir = dir.as_ref();
            if !dir.exists() {
                continue;
            }

            let partial = Self::from_directory(dir)?;
            combined_index.extend(partial.index);
        }

        if combined_index.is_empty() {
            return Err(FerroError::Io {
                msg: "No indexed FASTA files found in any directory".to_string(),
            });
        }

        Ok(Self {
            base_to_versioned: build_base_to_versioned(&combined_index),
            index: combined_index,
            aliases: build_chromosome_aliases(),
            cdot_mapper: None,
            supplemental_cds: SupplementalCdsInfo::default(),
            canonical_overrides: CanonicalOverrides::default(),
            protein_index: FxHashMap::default(),
            transcript_cache: new_transcript_cache(),
        })
    }

    /// Load the GRCh37 cdot referenced by `cdot_grch37_json` (if present) as a
    /// secondary build on `mapper`. Best-effort: a missing key, a missing file,
    /// or a load error is logged and otherwise ignored, mirroring the
    /// LRG/supplemental loads. Returns the number of GRCh37 transcripts loaded
    /// (0 when nothing was loaded) so callers/tests can assert the outcome.
    ///
    /// Factored out of [`Self::from_manifest`] so it can be exercised without a
    /// full hermetic manifest (genome + transcript FASTAs).
    fn load_grch37_secondary_cdot(
        mapper: &mut CdotMapper,
        manifest: &serde_json::Value,
        resolve_path: &impl Fn(&str) -> PathBuf,
    ) -> usize {
        let Some(grch37_path_str) = manifest.get("cdot_grch37_json").and_then(|v| v.as_str())
        else {
            return 0;
        };
        let grch37_path = resolve_path(grch37_path_str);
        if !grch37_path.exists() {
            eprintln!(
                "Warning: cdot_grch37_json path does not exist: {}",
                grch37_path.display()
            );
            return 0;
        }
        eprintln!(
            "Loading GRCh37 cdot transcript metadata from {}...",
            grch37_path.display()
        );
        match mapper.load_secondary_build(&grch37_path, "GRCh37") {
            Ok(count) => {
                eprintln!("Loaded {} GRCh37 transcripts (secondary build)", count);
                count
            }
            Err(e) => {
                eprintln!("Warning: Failed to load GRCh37 cdot: {}", e);
                0
            }
        }
    }

    /// Create a provider from a manifest file
    pub fn from_manifest<P: AsRef<Path>>(manifest_path: P) -> Result<Self, FerroError> {
        let manifest_path = manifest_path.as_ref();

        // Get the directory containing the manifest - paths are relative to this
        let base_dir = manifest_path
            .parent()
            .unwrap_or(Path::new("."))
            .to_path_buf();

        // Helper to resolve a path relative to the manifest's directory
        let resolve_path = |path_str: &str| -> PathBuf {
            let path = PathBuf::from(path_str);
            if path.is_absolute() {
                path
            } else {
                base_dir.join(path)
            }
        };

        let file = File::open(manifest_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open manifest: {}", e),
        })?;

        let manifest: serde_json::Value =
            serde_json::from_reader(file).map_err(|e| FerroError::Io {
                msg: format!("Failed to parse manifest: {}", e),
            })?;

        let mut dirs = Vec::new();

        // Get transcript directory (paths in manifest are relative to manifest location)
        if let Some(fastas) = manifest.get("transcript_fastas").and_then(|v| v.as_array()) {
            if let Some(first) = fastas.first().and_then(|v| v.as_str()) {
                // Remove .gz extension if present
                let path_str = first.strip_suffix(".gz").unwrap_or(first);
                let path = resolve_path(path_str);
                if let Some(transcript_dir) = path.parent() {
                    dirs.push(transcript_dir.to_path_buf());
                }
            }
        }

        // Get genome directory (paths in manifest are relative to manifest location)
        if let Some(genome) = manifest.get("genome_fasta").and_then(|v| v.as_str()) {
            let path = resolve_path(genome);
            if let Some(genome_dir) = path.parent() {
                dirs.push(genome_dir.to_path_buf());
            }
        }

        // Get RefSeqGene directory (NG_ accessions)
        if let Some(fastas) = manifest.get("refseqgene_fastas").and_then(|v| v.as_array()) {
            if let Some(first) = fastas.first().and_then(|v| v.as_str()) {
                let path = resolve_path(first);
                if let Some(refseqgene_dir) = path.parent() {
                    if !dirs.contains(&refseqgene_dir.to_path_buf()) {
                        dirs.push(refseqgene_dir.to_path_buf());
                    }
                }
            }
        }

        // Get LRG directory (LRG_ accessions)
        if let Some(fastas) = manifest.get("lrg_fastas").and_then(|v| v.as_array()) {
            if let Some(first) = fastas.first().and_then(|v| v.as_str()) {
                let path = resolve_path(first);
                if let Some(lrg_dir) = path.parent() {
                    if !dirs.contains(&lrg_dir.to_path_buf()) {
                        dirs.push(lrg_dir.to_path_buf());
                    }
                }
            }
        }

        // Get supplemental directory (missing ClinVar transcripts)
        if let Some(supplemental) = manifest.get("supplemental_fasta").and_then(|v| v.as_str()) {
            let path = resolve_path(supplemental);
            if let Some(supplemental_dir) = path.parent() {
                if !dirs.contains(&supplemental_dir.to_path_buf()) {
                    dirs.push(supplemental_dir.to_path_buf());
                }
            }
        }

        if dirs.is_empty() {
            return Err(FerroError::Io {
                msg: "No FASTA directories found in manifest".to_string(),
            });
        }

        let mut provider = Self::from_directories(&dirs)?;

        // Load cdot transcript metadata if available
        if let Some(cdot_path_str) = manifest.get("cdot_json").and_then(|v| v.as_str()) {
            let cdot_path = resolve_path(cdot_path_str);
            if cdot_path.exists() {
                eprintln!(
                    "Loading cdot transcript metadata from {}...",
                    cdot_path.display()
                );
                match CdotMapper::load(&cdot_path) {
                    Ok(mut mapper) => {
                        eprintln!(
                            "Loaded {} transcripts with CDS metadata",
                            mapper.transcript_count()
                        );

                        // Load LRG to RefSeq mapping if available
                        if let Some(lrg_mapping_path) =
                            manifest.get("lrg_refseq_mapping").and_then(|v| v.as_str())
                        {
                            let lrg_path = resolve_path(lrg_mapping_path);
                            if lrg_path.exists() {
                                match mapper.load_lrg_mapping(&lrg_path) {
                                    Ok(count) => {
                                        eprintln!("Loaded {} LRG to RefSeq mappings", count);
                                    }
                                    Err(e) => {
                                        eprintln!("Warning: Failed to load LRG mapping: {}", e);
                                    }
                                }
                            }
                        }

                        // Load the GRCh37 cdot as a secondary build if the
                        // manifest provides one. GRCh37 and GRCh38 contig
                        // accessions are distinct, so the secondary build's
                        // transcripts are indexed for overlap discovery without
                        // disturbing the primary (GRCh38) build. Best-effort:
                        // warn but do not fail, like the LRG/supplemental loads.
                        Self::load_grch37_secondary_cdot(&mut mapper, &manifest, &resolve_path);

                        provider.cdot_mapper = Some(mapper);
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to load cdot data: {}", e);
                    }
                }
            }
        }

        // Load supplemental metadata if available (for old/superseded transcripts)
        if let Some(supplemental) = manifest.get("supplemental_fasta").and_then(|v| v.as_str()) {
            let supplemental_path = resolve_path(supplemental);
            // Metadata file is {fasta_basename}.metadata.json (e.g., clinvar_transcripts.metadata.json)
            let metadata_path = supplemental_path.with_extension("metadata.json");
            if metadata_path.exists() {
                eprintln!(
                    "Loading supplemental CDS metadata from {}...",
                    metadata_path.display()
                );
                match std::fs::read_to_string(&metadata_path) {
                    Ok(content) => match serde_json::from_str::<serde_json::Value>(&content) {
                        Ok(metadata) => {
                            if let Some(transcripts) =
                                metadata.get("transcripts").and_then(|v| v.as_object())
                            {
                                for (accession, info) in transcripts {
                                    let cds_start = info.get("cds_start").and_then(|v| v.as_u64());
                                    let cds_end = info.get("cds_end").and_then(|v| v.as_u64());
                                    let gene_symbol = info
                                        .get("gene_symbol")
                                        .and_then(|v| v.as_str())
                                        .map(|s| s.to_string());
                                    let sequence_length = info
                                        .get("sequence_length")
                                        .and_then(|v| v.as_u64())
                                        .unwrap_or(0);
                                    provider.supplemental_cds.transcripts.insert(
                                        accession.clone(),
                                        SupplementalTranscriptInfo {
                                            cds_start,
                                            cds_end,
                                            gene_symbol,
                                            sequence_length,
                                        },
                                    );
                                }
                                eprintln!(
                                    "Loaded {} supplemental transcripts with CDS metadata",
                                    provider.supplemental_cds.transcripts.len()
                                );
                            }
                        }
                        Err(e) => {
                            eprintln!("Warning: Failed to parse supplemental metadata: {}", e);
                        }
                    },
                    Err(e) => {
                        eprintln!("Warning: Failed to read supplemental metadata: {}", e);
                    }
                }
            }
        }

        // Protein FASTAs (issue #520, edge 3): index NP_/XP_ sequences for
        // get_protein_sequence — the translated-CDS-vs-canonical-protein check.
        if let Some(protein_paths) = manifest.get("protein_fastas").and_then(|v| v.as_array()) {
            for p in protein_paths.iter().filter_map(|v| v.as_str()) {
                let path = resolve_path(p);
                let fai = PathBuf::from(format!("{}.fai", path.display()));
                if !path.exists() || !fai.exists() {
                    eprintln!(
                        "Warning: protein FASTA or its .fai is missing: {}",
                        path.display()
                    );
                    continue;
                }
                match load_fai_index(&path, &fai) {
                    Ok(idx) => provider.protein_index.extend(idx),
                    Err(e) => eprintln!(
                        "Warning: Failed to index protein FASTA {}: {}",
                        path.display(),
                        e
                    ),
                }
            }
            if !provider.protein_index.is_empty() {
                eprintln!("Loaded {} protein sequences", provider.protein_index.len());
            }
        }

        // Canonical overrides (issue #520): load the authoritative-overrides
        // file and reconcile it into the cdot metadata so both the served
        // transcript and the projection path see the corrected CDS/protein.
        if let Some(overrides_ref) = manifest.get("canonical_overrides").and_then(|v| v.as_str()) {
            let path = resolve_path(overrides_ref);
            match std::fs::read_to_string(&path) {
                Ok(text) => match CanonicalOverrides::from_json(&text) {
                    Ok(ov) => {
                        eprintln!("Loaded {} canonical override records", ov.len());
                        provider.canonical_overrides = ov;
                        provider.reconcile_cdot_with_overrides();
                    }
                    Err(e) => eprintln!("Warning: Failed to parse canonical overrides: {}", e),
                },
                Err(e) => eprintln!(
                    "Warning: Failed to read canonical overrides {}: {}",
                    path.display(),
                    e
                ),
            }
        }

        Ok(provider)
    }

    /// Create a provider from a single FASTA file with cdot metadata
    ///
    /// This is the recommended method for normalizing intronic variants,
    /// as it combines genomic sequence data (from FASTA) with transcript
    /// metadata (from cdot) needed for coordinate conversions.
    ///
    /// # Arguments
    ///
    /// * `fasta_path` - Path to an indexed FASTA file (must have .fai index)
    /// * `cdot_path` - Path to a cdot JSON file (can be gzipped)
    ///
    /// # Example
    ///
    /// ```ignore
    /// let provider = MultiFastaProvider::with_cdot(
    ///     "reference.fa",
    ///     "cdot-0.2.32.refseq.GRCh38.json.gz",
    /// )?;
    /// ```
    pub fn with_cdot<P: AsRef<Path>, Q: AsRef<Path>>(
        fasta_path: P,
        cdot_path: Q,
    ) -> Result<Self, FerroError> {
        let fasta_path = fasta_path.as_ref();
        let cdot_path = cdot_path.as_ref();

        // Load FASTA index
        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        if !fai_path.exists() {
            return Err(FerroError::Io {
                msg: format!(
                    "FASTA index not found: {}. Run 'samtools faidx {}' to create it.",
                    fai_path.display(),
                    fasta_path.display()
                ),
            });
        }

        let index = load_fai_index(fasta_path, &fai_path)?;

        eprintln!("Loaded {} sequences from FASTA", index.len());

        // Load cdot (prefers bincode if available, falls back to JSON)
        let cdot_mapper = CdotMapper::load(cdot_path)?;
        eprintln!(
            "Loaded {} transcripts from cdot",
            cdot_mapper.transcript_count()
        );

        Ok(Self {
            base_to_versioned: build_base_to_versioned(&index),
            index,
            aliases: build_chromosome_aliases(),
            cdot_mapper: Some(cdot_mapper),
            supplemental_cds: SupplementalCdsInfo::default(),
            canonical_overrides: CanonicalOverrides::default(),
            protein_index: FxHashMap::default(),
            transcript_cache: new_transcript_cache(),
        })
    }

    /// Resolve a sequence name, trying various forms
    fn resolve_name(&self, name: &str) -> Option<String> {
        // Direct lookup
        if self.index.contains_key(name) {
            return Some(name.to_string());
        }

        // Try chromosome alias FIRST (before version fallback)
        // This ensures NC_000001.11 maps to chr1 rather than falling back to NC_000001.10
        if let Some(aliased) = self.aliases.get(name) {
            if self.index.contains_key(aliased) {
                return Some(aliased.clone());
            }
        }

        // Try adding/removing chr prefix
        let alt_name = if name.starts_with("chr") {
            name.strip_prefix("chr").unwrap().to_string()
        } else {
            format!("chr{}", name)
        };

        if self.index.contains_key(&alt_name) {
            return Some(alt_name);
        }

        // Try without version (e.g., NM_000546 -> NM_000546.6)
        // This is for transcript accessions, not chromosome accessions
        if !name.contains('.') {
            if let Some(versioned) = self.base_to_versioned.get(name) {
                return Some(versioned.clone());
            }
        }

        // Try with version stripped (fallback for transcripts with different versions)
        // Skip for chromosome accessions that should have been resolved by aliases above
        if let Some(base) = name.split('.').next() {
            // Only use version fallback for transcript accessions, not chromosomes
            // Chromosome accessions (NC_*) should use aliases, not version fallback
            if !base.starts_with("NC_") {
                if let Some(versioned) = self.base_to_versioned.get(base) {
                    return Some(versioned.clone());
                }
            }
        }

        None
    }

    /// Get a sequence from the appropriate FASTA file
    fn get_sequence_from_index(
        &self,
        entry: &FastaIndexEntry,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Validate coordinates
        if start >= entry.length {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Start position {} exceeds sequence length {} for {}",
                    start, entry.length, entry.name
                ),
            });
        }

        let actual_end = end.min(entry.length);
        if start >= actual_end {
            return Ok(String::new());
        }

        // Calculate file position
        let line_start = start / entry.line_bases;
        let byte_offset = start % entry.line_bases;
        let file_offset = entry.offset + line_start * entry.line_bytes + byte_offset;

        // Calculate how many bytes we need to read
        let seq_len = actual_end - start;
        let num_lines = (seq_len + byte_offset).div_ceil(entry.line_bases);
        let bytes_to_read = seq_len + num_lines; // Extra for newlines

        // Read from file
        let mut file = File::open(&entry.file_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open FASTA file: {}", e),
        })?;

        file.seek(SeekFrom::Start(file_offset))
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to seek in FASTA file: {}", e),
            })?;

        let mut buffer = vec![0u8; bytes_to_read as usize];
        file.read_exact(&mut buffer).map_err(|e| FerroError::Io {
            msg: format!("Failed to read from FASTA file: {}", e),
        })?;

        // Strip newlines and ASCII-uppercase in a single byte pass.
        //
        // The previous implementation built a `String` char-by-char with
        // `filter().map(|&b| b as char).collect()` and then called
        // `to_uppercase()`, which walked the string a second time and
        // re-allocated. Both passes show up in samply (~3% of total CPU on
        // the SNP fixture before this change; more on cold workloads).
        let target_len = seq_len as usize;
        let mut sequence_bytes: Vec<u8> = Vec::with_capacity(target_len);
        for &b in &buffer {
            if b != b'\n' && b != b'\r' {
                sequence_bytes.push(b.to_ascii_uppercase());
                if sequence_bytes.len() == target_len {
                    break;
                }
            }
        }
        // FASTA bases are ASCII; `to_ascii_uppercase()` is a no-op on
        // already-uppercase bytes. The `from_utf8` validity check is O(n)
        // but extremely fast for short-ASCII inputs and doesn't reallocate.
        let sequence = String::from_utf8(sequence_bytes).map_err(|e| FerroError::Io {
            msg: format!("FASTA contains non-UTF8 bytes for {}: {}", entry.name, e),
        })?;
        Ok(sequence)
    }

    /// Check if a sequence exists
    pub fn has_sequence(&self, name: &str) -> bool {
        self.resolve_name(name).is_some()
    }

    /// Get the length of a sequence
    pub fn sequence_length(&self, name: &str) -> Option<u64> {
        self.resolve_name(name)
            .and_then(|n| self.index.get(&n))
            .map(|e| e.length)
    }

    /// Get total number of sequences
    pub fn sequence_count(&self) -> usize {
        self.index.len()
    }

    /// Get the cdot mapper if one was loaded with this provider.
    pub fn cdot_mapper(&self) -> Option<&CdotMapper> {
        self.cdot_mapper.as_ref()
    }

    /// Get CDS info from supplemental metadata for old/superseded transcripts.
    /// Creates a synthetic single exon spanning the entire transcript since
    /// mRNA transcripts are already spliced.
    fn get_supplemental_cds_info(&self, accession: &str) -> TranscriptMetadata {
        use crate::reference::transcript::{Exon, Strand};

        if let Some(info) = self.supplemental_cds.transcripts.get(accession) {
            // Create a synthetic single exon spanning the entire transcript.
            // For mRNA transcripts, the entire sequence is exonic (already spliced).
            let exons = if info.sequence_length > 0 {
                vec![Exon {
                    number: 1,
                    start: 1,
                    end: info.sequence_length,
                    genomic_start: None,
                    genomic_end: None,
                }]
            } else {
                Vec::new()
            };

            TranscriptMetadata {
                gene_symbol: info.gene_symbol.clone(),
                strand: Strand::Plus, // Assume plus strand for supplemental transcripts
                cds_start: info.cds_start, // CDS coordinates are already 1-based from GenBank
                cds_end: info.cds_end,
                chromosome: None,
                exons,
                genomic_start: None,
                genomic_end: None,
                protein_id: None,
                exon_cigars: Vec::new(),
            }
        } else {
            TranscriptMetadata::default()
        }
    }

    /// Build a `Transcript` for `id` whose bases are synthesized from the genome
    /// FASTA by walking the cdot exon alignment. Used as a fallback when the
    /// transcript FASTA does not carry the requested versioned accession but
    /// cdot does (closes #331).
    ///
    /// Strand-aware. cdot exon arrays store genome bases in the genome's
    /// 5'→3' orientation regardless of transcript strand, and exons are sorted
    /// by `tx_start`. For a plus-strand transcript that ordering is also
    /// genome-ascending, so concatenating the per-exon genome windows in
    /// tx-order yields the transcript sequence directly. For a minus-strand
    /// transcript the same ordering is *genome-descending* — so we
    /// reverse-complement each exon's genome bases as we go and concatenate in
    /// tx-order. Concatenating first and reverse-complementing the whole string
    /// would emit the second exon's bases before the first exon's (closes #331
    /// review feedback).
    ///
    /// Returns `FerroError::GenomicReferenceNotAvailable` when the cdot-named
    /// contig is missing from the genome FASTA, and `FerroError::InvalidCoordinates`
    /// when any exon's genomic span is degenerate or the cdot record carries
    /// no exons. The cdot-named contig is consulted directly (not version-
    /// stripped); callers that rely on a different contig version should keep
    /// the genome FASTA in sync with cdot.
    ///
    /// `tx.strand` is honored; metadata fields are projected onto the returned
    /// `Transcript` using the same conventions as the FASTA-backed path so
    /// downstream consumers see one shape regardless of which path served
    /// them. `genome_build` is derived from `build_hint` when the caller
    /// supplied one; otherwise from the attached [`CdotMapper`]'s primary
    /// build (this entry point is only reachable through that mapper, so
    /// the cdot-less branch of
    /// [`genome_build_from_hint`](Self::genome_build_from_hint) cannot
    /// fire here). See [#389].
    fn synthesize_transcript_from_cdot(
        &self,
        id: &str,
        tx: &crate::data::cdot::CdotTranscript,
        build_hint: Option<&str>,
    ) -> Result<Transcript, FerroError> {
        use crate::reference::transcript::{Exon as TxExon, ManeStatus, Strand};
        use std::sync::OnceLock;

        if tx.exons.is_empty() {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "cdot record for {} has no exons; cannot synthesize bases",
                    id
                ),
            });
        }

        // Walk cdot's exons in tx-order. Layout: `[e[0], e[1], e[2], e[3]]` ==
        //   [genome_start (HGVS-value, 1-based incl),
        //    genome_end   (HGVS-value, 1-based excl == Exon.genomic_end + 1),
        //    tx_start     (0-based incl),
        //    tx_end       (0-based excl)]
        // (see src/data/cdot.rs CdotTranscript and from_transcripts docs).
        // `get_genomic_sequence` expects a 0-based half-open window, so we
        // subtract 1 from both bounds. For minus-strand transcripts we
        // revcomp per exon before concatenation — see fn doc.
        let is_minus = matches!(tx.strand, Strand::Minus);
        let mut sequence = String::new();
        for e in &tx.exons {
            let g_start_hgvs = e[0];
            let g_end_excl_hgvs = e[1];
            if g_end_excl_hgvs <= g_start_hgvs || g_start_hgvs == 0 {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!(
                        "cdot exon for {} has degenerate genomic span [{}, {})",
                        id, g_start_hgvs, g_end_excl_hgvs
                    ),
                });
            }
            let bases =
                self.get_genomic_sequence(&tx.contig, g_start_hgvs - 1, g_end_excl_hgvs - 1)?;
            if is_minus {
                sequence.push_str(&crate::sequence::reverse_complement(&bases));
            } else {
                sequence.push_str(&bases);
            }
        }

        // #471: the bases above are pulled straight from the genome without
        // applying the cdot CIGAR, so flag — at a severity matching the
        // available alignment evidence — when the reconstruction is at elevated
        // risk of diverging from the deposited transcript. The call site already
        // emitted a soft "synthesizing from cdot" notice; the gapless case adds
        // nothing further.
        match classify_cdot_synthesis_risk(&tx.exon_cigars) {
            CdotSynthesisRisk::UnappliedIndels {
                insertions,
                deletions,
            } => warn!(
                "{id}: cdot exon alignment encodes {insertions} inserted and {deletions} deleted \
                 base(s) that base synthesis does not apply — the synthesized sequence will \
                 diverge from the deposited transcript in length and content; treat normalize \
                 output for {id} as low-confidence (issue #471)"
            ),
            CdotSynthesisRisk::NoAlignmentGapData => warn!(
                "{id}: cdot record carries no per-exon alignment-gap (CIGAR) data — synthesized \
                 bases assume an exact genome match and cannot represent transcript-specific \
                 indels or substitutions, so normalize output for {id} may diverge from the \
                 deposited transcript (issues #471, #400)"
            ),
            CdotSynthesisRisk::Gapless => {}
        }

        // Project cdot exons + CDS metadata onto the transcript struct using the
        // same conventions as the FASTA-backed `get_transcript` path above so
        // downstream consumers see one shape regardless of which path served them.
        let exons: Vec<TxExon> = tx
            .exons
            .iter()
            .enumerate()
            .map(|(i, e)| TxExon {
                number: (i + 1) as u32,
                start: e[2],
                end: e[3],
                // NOTE: this mirrors the (pre-existing, likely off-by-one)
                // projection in the FASTA path above. cdot's `e[0]` is
                // HGVS-1-based per src/data/cdot.rs, so the `+1` here is
                // strictly speaking wrong, but matching the FASTA path keeps
                // downstream behavior consistent. Fixing the convention is
                // out of scope for #331 — track separately if needed.
                genomic_start: Some(e[0] + 1),
                genomic_end: Some(e[1]),
            })
            .collect();

        let (genomic_start, genomic_end) = {
            let min_start = exons.iter().filter_map(|e| e.genomic_start).min();
            let max_end = exons.iter().filter_map(|e| e.genomic_end).max();
            (min_start, max_end)
        };

        Ok(Transcript {
            id: id.to_string(),
            gene_symbol: tx.gene_name.clone(),
            strand: tx.strand,
            sequence: Some(sequence),
            cds_start: tx.cds_start.map(|s| s + 1),
            cds_end: tx.cds_end,
            exons,
            chromosome: Some(tx.contig.clone()),
            genomic_start,
            genomic_end,
            genome_build: self.genome_build_from_hint(build_hint),
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: tx.protein.clone(),
            exon_cigars: tx.exon_cigars.clone(),
            cached_introns: OnceLock::new(),
        })
    }
}

/// Risk that cdot-driven base synthesis ([`MultiFastaProvider::synthesize_transcript_from_cdot`])
/// produces bases that diverge from the deposited transcript.
///
/// cdot's GFF3 `Gap` CIGAR encodes only `M`/`I`/`D` operations (never
/// substitutions), and `synthesize_transcript_from_cdot` reconstructs bases
/// purely from the genome FASTA over each exon's genomic span — it does **not**
/// apply the CIGAR. So the reconstruction is faithful only when the alignment
/// is gapless; missing CIGARs or indel-bearing CIGARs mean the synthesized
/// bases can silently diverge from the GenBank-deposited transcript. See issue
/// #471 (and #400, the `NM_001166478.1` reconstruction mismatch).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CdotSynthesisRisk {
    /// Every exon carries a gapless (`M`-only) CIGAR: the synthesized length is
    /// exact. Only substitution-level differences — which cdot's GFF3 Gap CIGAR
    /// cannot represent — could still diverge. Lowest risk.
    Gapless,
    /// No per-exon alignment-gap data (the CIGAR vector is empty, or at least
    /// one exon has no CIGAR entry): the synthesis assumes an exact genome match
    /// and cannot represent any transcript-specific indel or substitution.
    NoAlignmentGapData,
    /// At least one exon's CIGAR encodes indels that base synthesis does not
    /// apply, so the synthesized length and content provably diverge.
    UnappliedIndels { insertions: u64, deletions: u64 },
}

/// Classify the divergence risk of cdot base synthesis from the per-exon
/// CIGARs. Pure; see [`CdotSynthesisRisk`]. Indels anywhere dominate (strongest,
/// most actionable signal); otherwise fully-covered gapless CIGARs are low risk
/// and any missing per-exon data is elevated.
pub(crate) fn classify_cdot_synthesis_risk(
    exon_cigars: &[Option<Vec<crate::data::cdot::CigarOp>>],
) -> CdotSynthesisRisk {
    use crate::data::cdot::CigarOp;
    let mut insertions = 0u64;
    let mut deletions = 0u64;
    for op in exon_cigars.iter().flatten().flatten() {
        match op {
            CigarOp::Insertion(n) => insertions += n,
            CigarOp::Deletion(n) => deletions += n,
            CigarOp::Match(_) => {}
        }
    }
    if insertions > 0 || deletions > 0 {
        CdotSynthesisRisk::UnappliedIndels {
            insertions,
            deletions,
        }
    } else if !exon_cigars.is_empty()
        && exon_cigars
            .iter()
            .all(|c| matches!(c, Some(ops) if !ops.is_empty()))
    {
        CdotSynthesisRisk::Gapless
    } else {
        CdotSynthesisRisk::NoAlignmentGapData
    }
}

/// Metadata about a transcript gathered from cdot and/or supplemental sources.
///
/// Groups the fields needed to construct a [`Transcript`] from the reference provider,
/// replacing what was previously a 9-element tuple.
#[derive(Default)]
struct TranscriptMetadata {
    gene_symbol: Option<String>,
    strand: crate::reference::transcript::Strand,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    chromosome: Option<String>,
    exons: Vec<crate::reference::transcript::Exon>,
    genomic_start: Option<u64>,
    genomic_end: Option<u64>,
    protein_id: Option<String>,
    exon_cigars: Vec<Option<Vec<crate::data::cdot::CigarOp>>>,
}

impl MultiFastaProvider {
    /// Core transcript-loading routine. When `build_hint` is `Some`, the cdot
    /// lookup is restricted to that genome build via
    /// [`CdotMapper::get_transcript_on_build`]; otherwise it uses the primary
    /// build via [`CdotMapper::get_transcript`].
    ///
    /// This is the shared body of [`get_transcript`](Self::get_transcript) and
    /// [`get_transcript_for_variant`](Self::get_transcript_for_variant). See
    /// issue #332: when the input carries an NG/NC `genomic_context`, the
    /// caller passes the inferred build so the transcript's `chromosome` is
    /// populated against the correct alignment.
    /// Apply any loaded canonical overrides to a freshly-built transcript,
    /// correcting its CDS/protein metadata in place (issue #520). No-op when no
    /// overrides are loaded. Called at the single [`get_transcript_on_build`]
    /// chokepoint so every served transcript — lenient, strict, or
    /// variant-aware — is corrected before it is cached/translated.
    fn apply_canonical_overrides_to(&self, tx: &mut Transcript) {
        if !self.canonical_overrides.is_empty() {
            crate::reference::validate::apply_canonical_overrides(tx, &self.canonical_overrides);
        }
    }

    /// Reconcile the cdot mapper's parallel CDS/`protein` with the canonical
    /// overrides, so the projection path — which reads cdot's fields
    /// preferentially over the served `Transcript` — also sees the corrected
    /// values. Only reconciles a fully-coding authoritative record whose
    /// served FASTA length matches (so the coordinates apply); authoritative
    /// 1-based bounds are converted to cdot's 0-based start / 0-based-exclusive
    /// end convention.
    ///
    /// Reconciles when the coordinates will index canonical bases. The served
    /// length is the FASTA index length when the transcript is FASTA-backed,
    /// else the cdot exon-alignment length (for a cdot-synthesis-only transcript
    /// with no FASTA entry). It reconciles when that length matches the
    /// authoritative length, OR when the override carries the canonical
    /// `sequence` (the served bases are replaced on read regardless). (cdot
    /// exons themselves are not rewritten — the same exon-staleness caveat as
    /// the served-`Transcript` correction.)
    fn reconcile_cdot_with_overrides(&mut self) {
        use crate::data::cdot::CdotTranscript;
        if self.cdot_mapper.is_none() {
            return;
        }
        // Compute corrections under immutable borrows, then apply them.
        let mut corrections: Vec<(String, CdotTranscript)> = Vec::new();
        {
            let cdot = self.cdot_mapper.as_ref().unwrap();
            for (acc, auth) in &self.canonical_overrides.records {
                let (Some(cds_start), Some(cds_end), Some(protein)) =
                    (auth.cds_start, auth.cds_end, auth.protein_id.as_ref())
                else {
                    continue; // not fully coding → don't half-correct
                };
                // 1-based start → 0-based; guard against a degenerate `0` start
                // (unsigned underflow) from a malformed overrides file.
                let Some(cds_start_0based) = cds_start.checked_sub(1) else {
                    continue;
                };
                // Reconcile when the coordinates will apply to canonical bases.
                // The served length is the FASTA index length if the transcript
                // is FASTA-backed, else — for a cdot-synthesis-only transcript
                // with no FASTA entry — the length the read path actually
                // produces via `synthesize_transcript_from_cdot`: the summed exon
                // *genomic* span (`genomic_end - genomic_start`, i.e. `e[1] -
                // e[0]`). The raw `max(tx_end)` can disagree with that served
                // length when an exon's genomic and transcript spans differ, so
                // `apply_canonical_overrides` would gate on a length the server
                // never produced. Empty exons → `None` (synthesis errors), so
                // reconciliation is skipped. Proceed when the length matches the
                // authoritative length, or when the override carries the
                // canonical `sequence` (the served bases are replaced on read
                // regardless).
                let served_len = self.index.get(acc).map(|e| e.length).or_else(|| {
                    cdot.get_transcript(acc).and_then(|t| {
                        if t.exons.is_empty() {
                            None
                        } else {
                            Some(t.exons.iter().map(|e| e[1].saturating_sub(e[0])).sum())
                        }
                    })
                });
                let length_matches = served_len == Some(auth.tx_length);
                if !length_matches && auth.sequence.is_none() {
                    continue;
                }
                // Require an *exact* cdot record for `acc`: `get_transcript`
                // version-falls-back to a sibling, and cloning that sibling
                // under `acc` would carry over the wrong exon alignment / contig
                // metadata. Only reconcile when cdot holds this exact version.
                if cdot.has_transcript_exact(acc) {
                    let existing = cdot
                        .get_transcript(acc)
                        .expect("exact cdot presence checked above");
                    corrections.push((
                        acc.clone(),
                        CdotTranscript {
                            cds_start: Some(cds_start_0based), // 1-based incl → 0-based incl
                            cds_end: Some(cds_end),            // 1-based incl == 0-based excl
                            protein: Some(protein.clone()),
                            ..existing.clone()
                        },
                    ));
                }
            }
        }
        let cdot = self.cdot_mapper.as_mut().unwrap();
        for (acc, rec) in corrections {
            cdot.add_transcript(acc, rec);
        }
    }

    /// Build a transcript and apply canonical overrides — the single chokepoint
    /// every public entry point funnels through.
    fn get_transcript_on_build(
        &self,
        id: &str,
        build_hint: Option<&str>,
    ) -> Result<Transcript, FerroError> {
        let mut tx = self.get_transcript_on_build_inner(id, build_hint)?;
        self.apply_canonical_overrides_to(&mut tx);
        Ok(tx)
    }

    /// Memoized [`Self::get_transcript_on_build`].
    ///
    /// On a cache hit returns a cheap [`Arc`] clone of the previously resolved
    /// transcript; on a miss it resolves (materializing the full sequence and
    /// rebuilding cdot metadata), stores the result, and returns it. The cached
    /// value is exactly what `get_transcript_on_build` would have produced —
    /// including applied canonical overrides — so callers cannot observe any
    /// behavioral difference, only reduced work on repeat lookups. The cache key
    /// `(id, build_hint)` matches the resolution inputs so distinct builds of the
    /// same accession never collide.
    fn get_transcript_on_build_cached(
        &self,
        id: &str,
        build_hint: Option<&str>,
    ) -> Result<Arc<Transcript>, FerroError> {
        let key = (id.to_string(), build_hint.map(str::to_string));

        if let Some(hit) = self
            .transcript_cache
            .lock()
            .expect("transcript cache mutex poisoned")
            .get(&key)
        {
            return Ok(Arc::clone(hit));
        }

        // Resolve outside the lock so concurrent misses for *different*
        // transcripts don't serialize on the expensive FASTA read.
        // NOTE: concurrent same-key misses both resolve; the second put is a
        // no-op in value terms (LruCache overwrites with a semantically
        // identical Arc, and all callers already hold their own Arc clone).
        let tx = Arc::new(self.get_transcript_on_build(id, build_hint)?);

        self.transcript_cache
            .lock()
            .expect("transcript cache mutex poisoned")
            .put(key, Arc::clone(&tx));
        Ok(tx)
    }

    fn get_transcript_on_build_inner(
        &self,
        id: &str,
        build_hint: Option<&str>,
    ) -> Result<Transcript, FerroError> {
        use crate::reference::transcript::{Exon as TxExon, ManeStatus};
        use std::sync::OnceLock;

        if let Some(resolved) = self.resolve_name(id) {
            if let Some(entry) = self.index.get(&resolved) {
                // Get the full sequence
                let sequence = self.get_sequence_from_index(entry, 0, entry.length)?;

                // Try to get metadata from cdot
                let meta = if let Some(ref cdot) = self.cdot_mapper {
                    let cdot_tx_opt = match build_hint {
                        Some(b) => cdot.get_transcript_on_build(&resolved, b),
                        None => cdot.get_transcript(&resolved),
                    };
                    if let Some(tx) = cdot_tx_opt {
                        // Convert cdot exons to transcript exons
                        // cdot internal format: [genome_start, genome_end, tx_start, tx_end]
                        // COORDINATE SYSTEMS:
                        // - tx_start/tx_end: 1-based (use directly, no conversion needed)
                        // - genome_start: 0-based, convert to 1-based by adding 1
                        // - genome_end: 0-based exclusive = 1-based inclusive (no conversion)
                        let exons: Vec<TxExon> = tx
                            .exons
                            .iter()
                            .enumerate()
                            .map(|(i, e)| TxExon {
                                number: (i + 1) as u32,
                                start: e[2], // tx_start: already 1-based
                                end: e[3],   // tx_end: already 1-based
                                genomic_start: Some(e[0] + 1), // genome_start: 0-based → 1-based
                                genomic_end: Some(e[1]), // genome_end: 0-based excl = 1-based incl
                            })
                            .collect();

                        // Validate exon continuity - for mRNA transcripts, tx coordinates
                        // should be contiguous (no gaps). If there are gaps, the cdot data
                        // may be corrupted, so fall back to supplemental data.
                        let has_gaps = exons.windows(2).any(|w| w[1].start != w[0].end + 1);

                        if has_gaps && self.supplemental_cds.transcripts.contains_key(&resolved) {
                            // cdot has gaps, use supplemental data instead
                            warn!(
                                "cdot exon data for {} has gaps in transcript coordinates, \
                                 using supplemental CDS metadata",
                                resolved
                            );
                            self.get_supplemental_cds_info(&resolved)
                        } else {
                            // Compute transcript genomic bounds from exons
                            let (genomic_start, genomic_end) = if !exons.is_empty() {
                                let min_start = exons.iter().filter_map(|e| e.genomic_start).min();
                                let max_end = exons.iter().filter_map(|e| e.genomic_end).max();
                                (min_start, max_end)
                            } else {
                                (None, None)
                            };

                            TranscriptMetadata {
                                gene_symbol: tx.gene_name.clone(),
                                strand: tx.strand,
                                // CDS coordinates: cdot 0-based → transcript 1-based
                                cds_start: tx.cds_start.map(|s| s + 1), // 0-based → 1-based
                                cds_end: tx.cds_end, // 0-based exclusive = 1-based inclusive
                                chromosome: Some(tx.contig.clone()),
                                exons,
                                genomic_start,
                                genomic_end,
                                protein_id: tx.protein.clone(),
                                exon_cigars: tx.exon_cigars.clone(),
                            }
                        }
                    } else {
                        // Check supplemental CDS for old/superseded transcripts
                        if self.supplemental_cds.transcripts.contains_key(&resolved) {
                            warn!(
                                "{} not in cdot data, using supplemental CDS metadata",
                                resolved
                            );
                        }
                        self.get_supplemental_cds_info(&resolved)
                    }
                } else {
                    // Check supplemental CDS for old/superseded transcripts (no cdot mapper)
                    if self.supplemental_cds.transcripts.contains_key(&resolved) {
                        warn!(
                            "{} not in cdot data (no cdot mapper), using supplemental CDS metadata",
                            resolved
                        );
                    }
                    self.get_supplemental_cds_info(&resolved)
                };

                return Ok(Transcript {
                    id: resolved,
                    gene_symbol: meta.gene_symbol,
                    strand: meta.strand,
                    sequence: Some(sequence),
                    cds_start: meta.cds_start,
                    cds_end: meta.cds_end,
                    exons: meta.exons,
                    chromosome: meta.chromosome,
                    genomic_start: meta.genomic_start,
                    genomic_end: meta.genomic_end,
                    genome_build: self.genome_build_from_hint(build_hint),
                    mane_status: ManeStatus::default(),
                    refseq_match: None,
                    ensembl_match: None,
                    protein_id: meta.protein_id,
                    exon_cigars: meta.exon_cigars,
                    cached_introns: OnceLock::new(),
                });
            }
        }

        // cdot-driven fallback: when the transcript FASTA has no entry for
        // `id` (most often because the requested version is missing) but cdot
        // does carry the exon alignment, synthesize transcript bases by
        // walking the alignment against the genome FASTA. Strand-aware.
        // Closes #331.
        //
        // When `build_hint` is `Some`, restrict the cdot lookup to that
        // genome build so the synthesized exon coordinates and `contig`
        // come from the same alignment the caller asked for; otherwise
        // fall back to the primary-build view. See #389.
        if let Some(ref cdot) = self.cdot_mapper {
            let cdot_tx = match build_hint {
                Some(b) => cdot.get_transcript_on_build(id, b),
                None => cdot.get_transcript(id),
            };
            if let Some(tx) = cdot_tx {
                // cdot.get_transcript does its own version-strip fallback
                // (NM_FOO.1 -> NM_FOO.2 if only .2 is loaded). If that fired,
                // the bases we synthesize will come from a *different* cdot
                // record than the caller asked for — surface this so callers
                // don't get silent cross-version mismatches. The base
                // accession comparison is conservative: it warns whenever cdot
                // does not store the exact-version key.
                if !cdot.has_transcript_exact(id) {
                    warn!(
                        "{} not in transcript FASTA and not in cdot at this exact version; \
                         synthesizing bases via cdot version fallback (the cdot record \
                         may correspond to a sibling version)",
                        id
                    );
                } else {
                    warn!(
                        "{} not in transcript FASTA; synthesizing bases from cdot \
                         exon alignment against the genome FASTA",
                        id
                    );
                }
                return self.synthesize_transcript_from_cdot(id, tx, build_hint);
            }
        }

        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    /// Resolve `build_hint` to a [`GenomeBuild`] for a returned [`Transcript`].
    ///
    /// Precedence:
    ///   1. Explicit `build_hint` is consulted first. Recognized values
    ///      `"GRCh37"` / `"GRCh38"` map to the matching enum variant;
    ///      anything else — including the empty string, non-canonical
    ///      casing (`"grch38"`, `"hg19"`), or an unknown build name — maps
    ///      to [`GenomeBuild::Unknown`]. The cdot primary build is *not*
    ///      consulted as a fallback when a hint is supplied; callers
    ///      wanting that behavior must pass `None`.
    ///   2. With `None`, the attached [`CdotMapper`]'s primary build is
    ///      used when known (matched against the same canonical names).
    ///      An attached cdot whose `primary_build()` is `None` (e.g. a
    ///      bincode snapshot that pre-dates build tracking) surfaces as
    ///      [`GenomeBuild::Unknown`] rather than fabricating a default —
    ///      see #389 follow-up: the prior behavior of silently stamping
    ///      `GRCh38` re-introduced the mis-tag this PR removes.
    ///   3. With `None` and no cdot mapper attached, default to GRCh38
    ///      (preserves the historical FASTA-only behavior).
    ///
    /// Strict canonical-casing matches the format produced by
    /// [`infer_build_from_parent`](Self::infer_build_from_parent), which
    /// is the only in-tree caller that synthesizes hints; foreign hint
    /// sources should normalize before calling.
    ///
    /// Centralizing this here keeps the main FASTA-backed path and the
    /// cdot synthesize fallback in agreement (#389).
    fn genome_build_from_hint(
        &self,
        build_hint: Option<&str>,
    ) -> crate::reference::transcript::GenomeBuild {
        use crate::reference::transcript::GenomeBuild;
        // 1. Explicit hint wins.
        if let Some(h) = build_hint {
            return match h {
                "GRCh37" => GenomeBuild::GRCh37,
                "GRCh38" => GenomeBuild::GRCh38,
                _ => GenomeBuild::Unknown,
            };
        }
        // 2/3. Consult the attached cdot (if any). An unknown primary
        // build is surfaced as `Unknown` so callers don't get a wrong
        // canonical stamp on a snapshot that lacks build provenance.
        match self.cdot_mapper.as_ref().and_then(|c| c.primary_build()) {
            Some("GRCh37") => GenomeBuild::GRCh37,
            Some("GRCh38") => GenomeBuild::GRCh38,
            Some(_) => GenomeBuild::Unknown,
            None if self.cdot_mapper.is_some() => GenomeBuild::Unknown,
            None => GenomeBuild::GRCh38,
        }
    }

    /// Infer the genome build for an `NC_*` parent accession.
    ///
    /// Thin wrapper around the shared
    /// [`crate::liftover::aliases::infer_genome_build_from_accession`] —
    /// kept on the impl so existing call sites and tests stay terse;
    /// see the free function for the canonical contract.
    fn infer_build_from_parent(parent: &crate::hgvs::variant::Accession) -> Option<&'static str> {
        crate::liftover::aliases::infer_genome_build_from_accession(parent)
    }

    /// Whether `id` names a transcript whose bases are directly available at
    /// the *exact* requested versioned accession.
    ///
    /// This is the authoritative-source notion of "exact" for the FASTA /
    /// supplemental index: the FASTA index ([`Self::index`]) is keyed by the
    /// full versioned accession (e.g. `NM_212556.4`), and supplemental
    /// transcript FASTAs are loaded into that same map, so a direct
    /// [`HashMap::contains_key`] is exactly "we ship these bases at this
    /// version". It deliberately does *not* consult
    /// [`resolve_name`](Self::resolve_name), which performs version-strip and
    /// chromosome-alias fallbacks — those are the silent substitutions strict
    /// resolution exists to reject.
    ///
    /// Shared by the strict path ([`get_transcript_strict`](Self::get_transcript_strict))
    /// and exposed for callers that want to spot-check version pinning without
    /// materializing a [`Transcript`]. See design pillar 3 (#471 / #478).
    pub fn has_transcript_exact(&self, id: &str) -> bool {
        self.index.contains_key(id)
    }

    /// Resolve a transcript ONLY if its bases are directly available at the
    /// exact requested versioned accession; otherwise return
    /// [`FerroError::TranscriptVersionNotExact`] WITHOUT falling back to a
    /// sibling version or a cdot-genome reconstruction.
    ///
    /// This is an additive, opt-in capability for callers (e.g. a conformance
    /// harness) that must pin the reference version: the lenient
    /// [`get_transcript`](ReferenceProvider::get_transcript) chain (exact FASTA
    /// → version-strip fallback → cdot exon-alignment synthesis) silently
    /// substitutes a sibling version or a reconstructed transcript and only
    /// emits a `warn!`, which produces "phantom" divergences in conformance
    /// suites (e.g. a manifest shipping `NG_029146.2` serving a request for
    /// `NG_029146.1`). Strict resolution turns that soft signal into a hard
    /// `Err`.
    ///
    /// Exactness is determined by [`has_transcript_exact`](Self::has_transcript_exact):
    /// the FASTA/supplemental index must carry the full versioned accession
    /// directly. When the entry is present, this delegates to the same
    /// transcript-construction body as the lenient path
    /// ([`get_transcript_on_build`](Self::get_transcript_on_build)) so the
    /// returned [`Transcript`] is byte-for-byte identical to what the lenient
    /// path would yield for an exact hit — only the *acceptance criterion*
    /// differs. cdot metadata (CDS, exons) is still attached when available; it
    /// is the cdot-genome base *synthesis* fallback that strict mode refuses,
    /// not cdot metadata enrichment of an exact FASTA entry.
    ///
    /// The lenient `get_transcript` and the [`ReferenceProvider`] trait impl
    /// are unchanged.
    ///
    /// Canonical overrides (#520) still apply: strict resolution is about the
    /// version-exactness of the served *bases*, while a canonical override
    /// corrects CDS/`protein` *metadata* on the exact-version record — the two
    /// are orthogonal, so a strict caller still gets the canonical metadata.
    pub fn get_transcript_strict(&self, id: &str) -> Result<Transcript, FerroError> {
        if !self.has_transcript_exact(id) {
            return Err(FerroError::TranscriptVersionNotExact {
                requested: id.to_string(),
            });
        }
        // Exact FASTA entry present: reuse the shared construction body. The
        // exactness gate above guarantees `resolve_name(id)` returns `id`
        // itself (direct index hit wins in `resolve_name`), so this takes the
        // FASTA-backed branch and never the cdot-synthesis fallback.
        self.get_transcript_on_build(id, None)
    }

    /// Run the translated-CDS-vs-canonical-protein check (#520, the strategy-B
    /// depth tier) for each `(transcript_id, protein_accession)` pair, fetching
    /// the served transcript and the canonical protein from this provider's
    /// ingested protein FASTAs. Pairs whose transcript or protein is
    /// unavailable are skipped — translation validation is best-effort
    /// defense-in-depth, and protein data may not be ingested.
    pub fn validate_translations(
        &self,
        pairs: &[(String, String)],
    ) -> Vec<crate::reference::validate::TranscriptAnomaly> {
        use crate::reference::validate::validate_translation_against_protein;
        let mut anomalies = Vec::new();
        for (tx_id, protein_acc) in pairs {
            // Only validate when we serve this *exact* transcript version.
            // `get_transcript` version-falls-back to a sibling (or synthesizes
            // from cdot), and translating those bases against the canonical
            // protein would raise false `TranslationMismatch` anomalies for a
            // version we don't actually serve (the #505 wrong-attribution guard).
            if !self.has_transcript_version_exact(tx_id) {
                continue;
            }
            let Ok(tx) = self.get_transcript(tx_id) else {
                continue;
            };
            // `0, u64::MAX` fetches the full protein (the index read clamps the
            // end to the sequence length).
            let Ok(protein) = self.get_protein_sequence(protein_acc, 0, u64::MAX) else {
                continue;
            };
            if let Some(anomaly) = validate_translation_against_protein(&tx, &protein) {
                anomalies.push(anomaly);
            }
        }
        anomalies
    }
}

impl ReferenceProvider for MultiFastaProvider {
    fn get_transcript(&self, id: &str) -> Result<Arc<Transcript>, FerroError> {
        self.get_transcript_on_build_cached(id, None)
    }

    fn get_transcript_for_variant(
        &self,
        variant: &crate::hgvs::variant::HgvsVariant,
    ) -> Result<Arc<Transcript>, FerroError> {
        let Some(accession) = variant.accession() else {
            // No accession at all (e.g. NullAllele) — fall through to the
            // default impl which produces a clear ReferenceNotFound error.
            return self.get_transcript_on_build_cached(&format!("{}", variant), None);
        };
        let tx_id = accession.transcript_accession();

        // No parent → exactly the existing behavior.
        let Some(parent) = accession.genomic_context.as_deref() else {
            return self.get_transcript_on_build_cached(&tx_id, None);
        };

        // Decide a probe order based on the parent class.
        //   NC_*  : infer build from the parent version (GRCh37/GRCh38).
        //           Probe the inferred build first, then the other build.
        //   NG_*  : build-agnostic. Probe GRCh38 first (mutalyzer policy),
        //           then GRCh37.
        //   other : fall through to plain lookup (no extra information).
        let probe_order: &[&str] = match Self::infer_build_from_parent(parent) {
            Some("GRCh38") => &["GRCh38", "GRCh37"],
            Some("GRCh37") => &["GRCh37", "GRCh38"],
            _ => {
                // NG_* and unrecognized parents land here.
                if &*parent.prefix == "NG" {
                    &["GRCh38", "GRCh37"]
                } else {
                    return self.get_transcript_on_build_cached(&tx_id, None);
                }
            }
        };

        let mut last_err: Option<FerroError> = None;
        for build in probe_order {
            match self.get_transcript_on_build_cached(&tx_id, Some(build)) {
                Ok(tx) if tx.chromosome.is_some() => return Ok(tx),
                Ok(tx) => last_err = Some(FerroError::ReferenceNotFound { id: tx.id.clone() }),
                Err(e) => last_err = Some(e),
            }
        }

        // Final fallback: try without a build hint (uses cdot's primary build
        // and any supplemental/FASTA-only data). Preserves the historical
        // behavior for transcripts that only live in the primary build.
        match self.get_transcript_on_build_cached(&tx_id, None) {
            Ok(tx) => Ok(tx),
            Err(e) => Err(last_err.unwrap_or(e)),
        }
    }

    fn has_transcript_version_exact(&self, id: &str) -> bool {
        // Version-exact = the bases the read path would actually serve for
        // `id` correspond to the *exact* requested version, with no
        // cross-version substitution. Two ways that holds:
        //
        //  1. The FASTA/supplemental index ships the exact versioned accession.
        //  2. The FASTA index has no entry that `resolve_name` would resolve to
        //     (so no version-strip sibling shadows `id`) AND cdot stores the
        //     exon alignment at the exact version — then `get_transcript_on_build`
        //     reconstructs the bases from the genome at the exact version (#331).
        //
        // Condition 2 must check `resolve_name(id).is_none()`: the read path
        // resolves the FASTA (including a version-strip sibling) BEFORE it ever
        // reaches cdot synthesis, so a cdot-exact record shadowed by a FASTA
        // sibling would NOT be the bases actually served. Omitting that check
        // would green-light translating the sibling's bases as `id` — the exact
        // wrong-protein attribution #505 exists to prevent.
        //
        // Scope: the cdot arm consults the *primary-build* transcript map
        // (`CdotMapper::has_transcript_exact`). The bare-transcript protein path
        // (no NG/NC parent) routes through `get_transcript_on_build(id, None)`,
        // which is exactly the primary-build view this models. An NG/NC-parented
        // variant can probe alt-build maps that carry their own version-strip
        // fallback; covering that path's exactness is left to follow-up work.
        if self.has_transcript_exact(id) {
            return true;
        }
        self.resolve_name(id).is_none()
            && self
                .cdot_mapper
                .as_ref()
                .is_some_and(|cdot| cdot.has_transcript_exact(id))
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let resolved = self
            .resolve_name(id)
            .ok_or_else(|| FerroError::ReferenceNotFound { id: id.to_string() })?;

        let entry = self
            .index
            .get(&resolved)
            .ok_or_else(|| FerroError::ReferenceNotFound { id: id.to_string() })?;

        self.get_sequence_from_index(entry, start, end)
    }

    fn get_genomic_sequence(
        &self,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Try to find the contig in the FASTA index
        // First try exact match, then try with/without version
        let resolved =
            self.resolve_name(contig)
                .ok_or_else(|| FerroError::GenomicReferenceNotAvailable {
                    contig: contig.to_string(),
                    start,
                    end,
                })?;

        let entry =
            self.index
                .get(&resolved)
                .ok_or_else(|| FerroError::GenomicReferenceNotAvailable {
                    contig: contig.to_string(),
                    start,
                    end,
                })?;

        self.get_sequence_from_index(entry, start, end)
    }

    fn has_genomic_data(&self) -> bool {
        // Check if we have chromosome/contig sequences in the index
        // These may be named as NC_* (NCBI RefSeq) or chr* (UCSC style)
        self.index
            .keys()
            .any(|k| k.starts_with("NC_") || k.starts_with("chr"))
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        self.sequence_length(id)
            .ok_or_else(|| FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_protein_sequence(
        &self,
        accession: &str,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        // Exact versioned lookup only (no base-accession fallback) — a protein
        // accession is version-specific and callers supply versioned accessions.
        // (`MockProvider` does an unversioned fallback for testing; this
        // provider deliberately does not.)
        let entry = self.protein_index.get(accession).ok_or_else(|| {
            FerroError::ProteinReferenceNotAvailable {
                accession: accession.to_string(),
                start,
                end,
            }
        })?;
        self.get_sequence_from_index(entry, start, end)
    }

    fn has_protein_data(&self) -> bool {
        !self.protein_index.is_empty()
    }
}

/// Load a FASTA index (.fai) file
///
/// Note: This function canonicalizes the FASTA path to resolve symlinks and
/// provide absolute paths, which helps prevent path traversal issues.
fn load_fai_index<P: AsRef<Path>>(
    fasta_path: P,
    fai_path: P,
) -> Result<FxHashMap<String, FastaIndexEntry>, FerroError> {
    // Canonicalize the FASTA path to resolve symlinks and get absolute path
    // This helps prevent path traversal vulnerabilities
    let fasta_path = fasta_path
        .as_ref()
        .canonicalize()
        .map_err(|e| FerroError::Io {
            msg: format!(
                "Failed to canonicalize FASTA path '{}': {}",
                fasta_path.as_ref().display(),
                e
            ),
        })?;
    let fai_path = fai_path.as_ref();

    let file = File::open(fai_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open FAI file: {}", e),
    })?;
    let reader = BufReader::new(file);

    let mut index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read FAI line: {}", e),
        })?;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        // Parse index values with proper error handling instead of silent defaults
        let name = fields[0].to_string();
        let length: u64 = fields[1].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid length '{}' in FAI for sequence '{}'",
                fields[1], name
            ),
        })?;
        let offset: u64 = fields[2].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid offset '{}' in FAI for sequence '{}'",
                fields[2], name
            ),
        })?;
        let line_bases: u64 = fields[3].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid line_bases '{}' in FAI for sequence '{}'",
                fields[3], name
            ),
        })?;
        let line_bytes: u64 = fields[4].parse().map_err(|_| FerroError::Io {
            msg: format!(
                "Invalid line_bytes '{}' in FAI for sequence '{}'",
                fields[4], name
            ),
        })?;

        // Validate that critical fields are non-zero to prevent divide-by-zero
        if line_bases == 0 || line_bytes == 0 {
            return Err(FerroError::Io {
                msg: format!(
                    "Invalid FAI entry for '{}': line_bases={}, line_bytes={} (must be > 0)",
                    name, line_bases, line_bytes
                ),
            });
        }

        let entry = FastaIndexEntry {
            file_path: fasta_path.clone(),
            name: name.clone(),
            length,
            offset,
            line_bases,
            line_bytes,
        };

        index.insert(name, entry);
    }

    Ok(index)
}

/// Count unique files in the index
fn count_unique_files(index: &FxHashMap<String, FastaIndexEntry>) -> usize {
    let files: std::collections::HashSet<_> = index.values().map(|e| &e.file_path).collect();
    files.len()
}

/// Build chromosome name aliases
fn build_chromosome_aliases() -> FxHashMap<String, String> {
    let mut aliases: FxHashMap<String, String> = FxHashMap::default();

    // RefSeq to UCSC chromosome names (for genome)
    let refseq_chroms = [
        ("NC_000001.11", "chr1"),
        ("NC_000002.12", "chr2"),
        ("NC_000003.12", "chr3"),
        ("NC_000004.12", "chr4"),
        ("NC_000005.10", "chr5"),
        ("NC_000006.12", "chr6"),
        ("NC_000007.14", "chr7"),
        ("NC_000008.11", "chr8"),
        ("NC_000009.12", "chr9"),
        ("NC_000010.11", "chr10"),
        ("NC_000011.10", "chr11"),
        ("NC_000012.12", "chr12"),
        ("NC_000013.11", "chr13"),
        ("NC_000014.9", "chr14"),
        ("NC_000015.10", "chr15"),
        ("NC_000016.10", "chr16"),
        ("NC_000017.11", "chr17"),
        ("NC_000018.10", "chr18"),
        ("NC_000019.10", "chr19"),
        ("NC_000020.11", "chr20"),
        ("NC_000021.9", "chr21"),
        ("NC_000022.11", "chr22"),
        ("NC_000023.11", "chrX"),
        ("NC_000024.10", "chrY"),
        ("NC_012920.1", "chrM"),
    ];

    for (refseq, ucsc) in &refseq_chroms {
        aliases.insert(refseq.to_string(), ucsc.to_string());
        // Also map base accession (without version)
        if let Some(base) = refseq.split('.').next() {
            aliases.insert(base.to_string(), ucsc.to_string());
        }
    }

    // Also map bare chromosome names
    for i in 1..=22 {
        aliases.insert(i.to_string(), format!("chr{}", i));
    }
    aliases.insert("X".to_string(), "chrX".to_string());
    aliases.insert("Y".to_string(), "chrY".to_string());
    aliases.insert("M".to_string(), "chrM".to_string());
    aliases.insert("MT".to_string(), "chrM".to_string());

    aliases
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    // ----------------------------------------------------------------------
    // infer_build_from_parent unit coverage (closes review gap on #332)
    // ----------------------------------------------------------------------

    fn nc_accession(num: &str, version: u32) -> crate::hgvs::variant::Accession {
        crate::hgvs::variant::Accession::new("NC", num, Some(version))
    }

    fn ng_accession(num: &str, version: u32) -> crate::hgvs::variant::Accession {
        crate::hgvs::variant::Accession::new("NG", num, Some(version))
    }

    #[test]
    fn infer_build_chr17_grch38() {
        // NC_000017.11 is chr17 on GRCh38.
        assert_eq!(
            MultiFastaProvider::infer_build_from_parent(&nc_accession("000017", 11)),
            Some("GRCh38")
        );
    }

    #[test]
    fn infer_build_chr17_grch37() {
        // NC_000017.10 is chr17 on GRCh37.
        assert_eq!(
            MultiFastaProvider::infer_build_from_parent(&nc_accession("000017", 10)),
            Some("GRCh37")
        );
    }

    #[test]
    fn infer_build_chr14_uses_alias_table_not_version_heuristic() {
        // chr14 NC accessions are NC_000014.9 (GRCh38) / NC_000014.8 (GRCh37).
        // A naive ".11→GRCh38, .10→GRCh37" rule would mis-classify both.
        assert_eq!(
            MultiFastaProvider::infer_build_from_parent(&nc_accession("000014", 9)),
            Some("GRCh38"),
        );
        assert_eq!(
            MultiFastaProvider::infer_build_from_parent(&nc_accession("000014", 8)),
            Some("GRCh37"),
        );
    }

    #[test]
    fn infer_build_malformed_nc_returns_none() {
        // NC_000999.999 is not in the human alias table on any build.
        assert_eq!(
            MultiFastaProvider::infer_build_from_parent(&nc_accession("000999", 999)),
            None,
        );
    }

    #[test]
    fn infer_build_ng_parent_returns_none() {
        // NG accessions are build-agnostic; caller must probe builds.
        assert_eq!(
            MultiFastaProvider::infer_build_from_parent(&ng_accession("012772", 1)),
            None,
        );
    }

    // ---- #471: cdot base-synthesis divergence-risk classifier ----

    #[test]
    fn cdot_synthesis_risk_empty_is_no_gap_data() {
        // No CIGAR data at all (the common cdot case) → cannot represent any
        // transcript-vs-genome difference → elevated risk.
        assert_eq!(
            classify_cdot_synthesis_risk(&[]),
            CdotSynthesisRisk::NoAlignmentGapData
        );
    }

    #[test]
    fn cdot_synthesis_risk_all_none_is_no_gap_data() {
        assert_eq!(
            classify_cdot_synthesis_risk(&[None, None]),
            CdotSynthesisRisk::NoAlignmentGapData
        );
    }

    #[test]
    fn cdot_synthesis_risk_empty_inner_vec_is_no_gap_data() {
        // `Some(vec![])` is a CIGAR slot that is present but carries no operations:
        // there is no alignment data for that exon, so it must elevate to
        // NoAlignmentGapData rather than masquerade as Gapless (per the doc on
        // `CdotSynthesisRisk::NoAlignmentGapData`).
        assert_eq!(
            classify_cdot_synthesis_risk(&[Some(vec![])]),
            CdotSynthesisRisk::NoAlignmentGapData
        );
    }

    #[test]
    fn cdot_synthesis_risk_mixed_empty_inner_vec_is_no_gap_data() {
        use crate::data::cdot::CigarOp;
        // One exon has real alignment data, another has an empty CIGAR vector:
        // we cannot vouch for the whole transcript, so this is elevated, not gapless.
        let cigars = vec![Some(vec![CigarOp::Match(185)]), Some(vec![])];
        assert_eq!(
            classify_cdot_synthesis_risk(&cigars),
            CdotSynthesisRisk::NoAlignmentGapData
        );
    }

    #[test]
    fn cdot_synthesis_risk_all_match_is_gapless() {
        use crate::data::cdot::CigarOp;
        let cigars = vec![
            Some(vec![CigarOp::Match(185)]),
            Some(vec![CigarOp::Match(250)]),
        ];
        assert_eq!(
            classify_cdot_synthesis_risk(&cigars),
            CdotSynthesisRisk::Gapless
        );
    }

    #[test]
    fn cdot_synthesis_risk_partial_coverage_is_no_gap_data() {
        use crate::data::cdot::CigarOp;
        // One exon carries a CIGAR, one does not — we cannot vouch for the whole
        // transcript, so this is elevated, not gapless.
        let cigars = vec![Some(vec![CigarOp::Match(185)]), None];
        assert_eq!(
            classify_cdot_synthesis_risk(&cigars),
            CdotSynthesisRisk::NoAlignmentGapData
        );
    }

    #[test]
    fn cdot_synthesis_risk_insertion_is_unapplied_indels() {
        use crate::data::cdot::CigarOp;
        let cigars = vec![Some(vec![
            CigarOp::Match(185),
            CigarOp::Insertion(3),
            CigarOp::Match(250),
        ])];
        assert_eq!(
            classify_cdot_synthesis_risk(&cigars),
            CdotSynthesisRisk::UnappliedIndels {
                insertions: 3,
                deletions: 0
            }
        );
    }

    #[test]
    fn cdot_synthesis_risk_deletion_is_unapplied_indels() {
        use crate::data::cdot::CigarOp;
        let cigars = vec![Some(vec![CigarOp::Match(100), CigarOp::Deletion(5)])];
        assert_eq!(
            classify_cdot_synthesis_risk(&cigars),
            CdotSynthesisRisk::UnappliedIndels {
                insertions: 0,
                deletions: 5
            }
        );
    }

    #[test]
    fn cdot_synthesis_risk_indels_dominate_missing_data() {
        use crate::data::cdot::CigarOp;
        // An indel on one exon and a missing CIGAR on another → the indel is the
        // strongest, most actionable signal and must win.
        let cigars = vec![Some(vec![CigarOp::Insertion(2)]), None];
        assert_eq!(
            classify_cdot_synthesis_risk(&cigars),
            CdotSynthesisRisk::UnappliedIndels {
                insertions: 2,
                deletions: 0
            }
        );
    }

    #[test]
    fn test_load_fai_index() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        // Create FAI file
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        let index = load_fai_index(&fasta_path, &fai_path).unwrap();

        assert!(index.contains_key("NM_000001.1"));
        let entry = index.get("NM_000001.1").unwrap();
        assert_eq!(entry.length, 10);
        assert_eq!(entry.offset, 13);
    }

    #[test]
    fn test_multi_fasta_provider_from_directory() {
        let dir = tempdir().unwrap();

        // Create FASTA file
        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        assert!(provider.has_sequence("NM_000001.1"));
        assert!(provider.has_sequence("NM_000001")); // Without version
        assert_eq!(provider.sequence_length("NM_000001.1"), Some(10));
    }

    #[test]
    fn test_get_sequence() {
        let dir = tempdir().unwrap();

        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        let seq = provider.get_sequence("NM_000001.1", 0, 4).unwrap();
        assert_eq!(seq, "ATGC");

        let seq = provider.get_sequence("NM_000001", 0, 4).unwrap();
        assert_eq!(seq, "ATGC");
    }

    #[test]
    fn test_chromosome_aliases() {
        let aliases = build_chromosome_aliases();

        assert_eq!(aliases.get("NC_000001.11"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("NC_000001"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("1"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("X"), Some(&"chrX".to_string()));
    }

    #[test]
    fn test_multi_fasta_provider_get_sequence_length_returns_length_for_known_contig() {
        let dir = tempdir().unwrap();

        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        assert_eq!(provider.get_sequence_length("NM_000001.1").unwrap(), 10);
        // Also verify the unversioned alias resolves
        assert_eq!(provider.get_sequence_length("NM_000001").unwrap(), 10);
    }

    #[test]
    fn test_multi_fasta_provider_get_sequence_length_resolves_versioned_chromosome_alias() {
        let dir = tempdir().unwrap();

        // FASTA keyed as chrM (UCSC style); the versioned RefSeq accession
        // `NC_012920.1` must resolve to it via the `NC_012920 -> chrM` alias
        // after the version suffix is stripped.
        let fasta_path = dir.path().join("chrM.fa");
        let fai_path = dir.path().join("chrM.fa.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chrM").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "chrM\t10\t6\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        assert_eq!(provider.get_sequence_length("chrM").unwrap(), 10);
        assert_eq!(provider.get_sequence_length("NC_012920").unwrap(), 10);
        assert_eq!(provider.get_sequence_length("NC_012920.1").unwrap(), 10);
    }

    #[test]
    fn test_multi_fasta_provider_get_sequence_length_errors_for_unknown_id() {
        let dir = tempdir().unwrap();

        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        let err = provider
            .get_sequence_length("NM_NONEXISTENT.1")
            .unwrap_err();
        assert!(matches!(err, FerroError::ReferenceNotFound { .. }));
    }

    /// Helper: write a tempdir genome FASTA with `>NC_TEST.1` carrying 40
    /// bases of `AAAACCCCGGGGTTTT` repeating. Returns both the provider and
    /// the `TempDir` so the caller's drop scope keeps the FASTA on disk for
    /// the provider's lifetime.
    fn build_provider_with_test_genome() -> (MultiFastaProvider, tempfile::TempDir) {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("genome.fna");
        let fai_path = dir.path().join("genome.fna.fai");
        {
            let mut f = File::create(&fasta_path).unwrap();
            writeln!(f, ">NC_TEST.1").unwrap();
            writeln!(f, "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC").unwrap();
        }
        {
            let mut f = File::create(&fai_path).unwrap();
            // name<TAB>length<TAB>offset<TAB>line_bases<TAB>line_bytes
            // Header ">NC_TEST.1\n" is 11 bytes, so seq offset = 11.
            writeln!(f, "NC_TEST.1\t40\t11\t40\t41").unwrap();
        }
        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();
        (provider, dir)
    }

    /// Build a single-exon synthetic `Transcript` over NC_TEST.1[10..30]
    /// (1-based HGVS 11..=30) at the requested strand. Used by both the
    /// plus-strand and minus-strand single-exon fallback tests.
    fn build_single_exon_synthetic_tx(
        strand: crate::reference::transcript::Strand,
    ) -> crate::reference::transcript::Transcript {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Transcript};
        use std::sync::OnceLock;
        Transcript {
            id: "NM_TEST.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(20),
            exons: vec![Exon::with_genomic(1, 1, 20, 11, 30)],
            chromosome: Some("NC_TEST.1".to_string()),
            genomic_start: Some(11),
            genomic_end: Some(30),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    #[test]
    fn test_get_transcript_falls_back_to_cdot_alignment_plus_strand() {
        use crate::reference::transcript::Strand;
        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));

        // NM_TEST.1 is NOT in the FASTA index (only NC_TEST.1 is). cdot has the
        // exon alignment. Pre-fix this returns ReferenceNotFound; post-fix it
        // reconstructs the bases from the genome FASTA.
        let synthesized = provider
            .get_transcript("NM_TEST.1")
            .expect("cdot exon-alignment fallback should fetch transcript");

        assert_eq!(synthesized.id, "NM_TEST.1");
        assert_eq!(
            synthesized.sequence.as_deref(),
            Some("GGTTTTAAAACCCCGGGGTT")
        );
        assert!(matches!(synthesized.strand, Strand::Plus));
        assert_eq!(synthesized.gene_symbol.as_deref(), Some("TEST"));
        assert_eq!(synthesized.exons.len(), 1);
        // CDS metadata projection from cdot must survive the fallback.
        assert_eq!(synthesized.cds_start, Some(1));
        assert_eq!(synthesized.cds_end, Some(20));
        assert_eq!(synthesized.chromosome.as_deref(), Some("NC_TEST.1"));
    }

    #[test]
    fn test_get_transcript_falls_back_to_cdot_alignment_minus_strand() {
        use crate::reference::transcript::Strand;
        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Minus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));

        let synthesized = provider
            .get_transcript("NM_TEST.1")
            .expect("cdot exon-alignment fallback should fetch transcript");

        // Genome bases [10..30] on NC_TEST.1 are `GGTTTTAAAACCCCGGGGTT`;
        // minus-strand transcript orientation is the reverse complement:
        //   complement: CCAAAATTTTGGGGCCCCAA
        //   reversed:   AACCCCGGGGTTTTAAAACC
        assert_eq!(
            synthesized.sequence.as_deref(),
            Some("AACCCCGGGGTTTTAAAACC")
        );
        assert!(matches!(synthesized.strand, Strand::Minus));
    }

    #[test]
    fn test_get_transcript_falls_back_with_multi_exon_plus_strand() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let (mut provider, _kept) = build_provider_with_test_genome();

        // Two exons on + strand. Genome FASTA (1-based HGVS positions):
        //   pos:   1234 5678 9012 3456 7890 1234 5678 9012 3456 7890
        //   bases: AAAA CCCC GGGG TTTT AAAA CCCC GGGG TTTT AAAA CCCC
        //   exon1: tx 1..4   genome 1..4    bases at 1-based HGVS [1..=4]  = "AAAA"
        //   exon2: tx 5..12  genome 15..22  bases at 1-based HGVS [15..=22] = "TTAAAACC"
        // Expected concatenated transcript sequence: "AAAA" + "TTAAAACC" = "AAAATTAAAACC".
        let tx = Transcript {
            id: "NM_MULTI.1".to_string(),
            gene_symbol: Some("MULTI".to_string()),
            strand: Strand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![
                Exon::with_genomic(1, 1, 4, 1, 4),
                Exon::with_genomic(2, 5, 12, 15, 22),
            ],
            chromosome: Some("NC_TEST.1".to_string()),
            genomic_start: Some(1),
            genomic_end: Some(22),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));

        let synthesized = provider
            .get_transcript("NM_MULTI.1")
            .expect("multi-exon fallback should fetch transcript");
        assert_eq!(synthesized.sequence.as_deref(), Some("AAAATTAAAACC"));
        assert_eq!(synthesized.exons.len(), 2);
    }

    #[test]
    fn test_get_transcript_falls_back_with_multi_exon_minus_strand() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        let (mut provider, _kept) = build_provider_with_test_genome();

        // Two exons on **minus** strand. Genome FASTA (1-based HGVS positions):
        //   pos:   1234 5678 9012 3456 7890 1234 5678 9012 3456 7890
        //   bases: AAAA CCCC GGGG TTTT AAAA CCCC GGGG TTTT AAAA CCCC
        //
        // On minus strand the first tx exon corresponds to the *higher* genome
        // coords. Spec the cdot record with that ordering:
        //   exon1 (first in tx): genome 15..22 → bases "TTAAAACC"  (HGVS 15..=22)
        //   exon2 (second in tx): genome 1..4  → bases "AAAA"      (HGVS 1..=4)
        //
        // Tx-orientation per exon = revcomp of the genome bases:
        //   exon1 tx bases = revcomp("TTAAAACC") = "GGTTTTAA"
        //   exon2 tx bases = revcomp("AAAA")     = "TTTT"
        // Full transcript sequence = exon1_tx ++ exon2_tx = "GGTTTTAATTTT".
        //
        // The pre-fix concatenate-then-revcomp logic emits revcomp of
        // ("TTAAAACC" + "AAAA") = revcomp("TTAAAACCAAAA") = "TTTTGGTTTTAA",
        // which inverts the per-exon order — this test locks the per-exon
        // revcomp behavior in place.
        let tx = Transcript {
            id: "NM_MINUS_MULTI.1".to_string(),
            gene_symbol: Some("MINUS_MULTI".to_string()),
            strand: Strand::Minus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![
                Exon::with_genomic(1, 1, 8, 15, 22),
                Exon::with_genomic(2, 9, 12, 1, 4),
            ],
            chromosome: Some("NC_TEST.1".to_string()),
            genomic_start: Some(1),
            genomic_end: Some(22),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));

        let synthesized = provider
            .get_transcript("NM_MINUS_MULTI.1")
            .expect("multi-exon minus-strand fallback should fetch transcript");
        assert_eq!(synthesized.sequence.as_deref(), Some("GGTTTTAATTTT"));
        assert_eq!(synthesized.exons.len(), 2);
    }

    #[test]
    fn test_get_transcript_falls_back_errors_when_genome_contig_absent() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

        // Genome FASTA has SOME other contig but not the one cdot references.
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("decoy.fna");
        let fai_path = dir.path().join("decoy.fna.fai");
        {
            let mut f = File::create(&fasta_path).unwrap();
            writeln!(f, ">NC_DECOY.1").unwrap();
            writeln!(f, "ACGTACGTACGT").unwrap();
        }
        {
            let mut f = File::create(&fai_path).unwrap();
            writeln!(f, "NC_DECOY.1\t12\t12\t12\t13").unwrap();
        }
        let mut provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        let tx = Transcript {
            id: "NM_ORPHAN.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(10),
            exons: vec![Exon::with_genomic(1, 1, 10, 1, 10)],
            chromosome: Some("NC_MISSING.99".to_string()),
            genomic_start: Some(1),
            genomic_end: Some(10),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));

        let err = provider
            .get_transcript("NM_ORPHAN.1")
            .expect_err("missing genome contig must produce a typed error, not a panic");
        assert!(
            matches!(err, FerroError::GenomicReferenceNotAvailable { .. }),
            "expected GenomicReferenceNotAvailable, got {err:?}"
        );
    }

    #[test]
    fn test_get_transcript_falls_back_errors_on_empty_exon_array() {
        use crate::data::cdot::CdotTranscript;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        // Directly seat a CdotTranscript with no exons. `CdotMapper::from_transcripts`
        // would silently skip such an input, so we bypass it and write into the
        // mapper's internal table by going through `new()` + `add_transcript`.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_EMPTY.1".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_TEST.1".to_string(),
                strand: Strand::Plus,
                exons: Vec::new(),
                cds_start: None,
                cds_end: None,
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);

        let err = provider.get_transcript("NM_EMPTY.1").expect_err(
            "empty exon array must produce a typed error, not a silent empty Transcript",
        );
        assert!(
            matches!(err, FerroError::InvalidCoordinates { .. }),
            "expected InvalidCoordinates, got {err:?}"
        );
    }

    // ----------------------------------------------------------------------
    // genome_build propagation (closes #389 item 1)
    //
    // Pre-fix, `synthesize_transcript_from_cdot` hardcoded
    // `genome_build: GenomeBuild::GRCh38` and the main-path return mapped
    // `None` to GRCh38 as well. For a `CdotMapper` whose primary build is
    // GRCh37, both paths mis-tagged transcripts as GRCh38. These tests pin
    // build propagation through the cdot fallback (synthesize) and the
    // main FASTA-backed path, and exercise the centralized
    // `genome_build_from_hint` helper that both paths share.
    // ----------------------------------------------------------------------

    /// `cdot.primary_build == "GRCh37"`, transcript not in the transcript
    /// FASTA index → synthesize fallback runs. With no caller-supplied
    /// `build_hint`, the synthesized `Transcript` must inherit the cdot
    /// primary build, not the historical hardcoded GRCh38.
    #[test]
    fn test_synthesize_falls_back_to_cdot_primary_build_grch37() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::transcript::{GenomeBuild, Strand};

        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts_with_build(
            std::iter::once(&tx),
            "GRCh37",
        ));

        let synthesized = provider
            .get_transcript("NM_TEST.1")
            .expect("cdot exon-alignment fallback should fetch transcript");
        assert_eq!(
            synthesized.genome_build,
            GenomeBuild::GRCh37,
            "synthesized transcript must inherit cdot primary build, not default to GRCh38"
        );
    }

    /// Explicit `build_hint` is threaded through the synthesize fallback so
    /// the resulting `Transcript.genome_build` reflects the hint. Here the
    /// hint matches the cdot primary build (GRCh37 ↔ GRCh37) — exercises
    /// that the parameter is plumbed through, not merely inferred from
    /// the primary build.
    #[test]
    fn test_synthesize_honors_explicit_build_hint_grch37() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::transcript::{GenomeBuild, Strand};

        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts_with_build(
            std::iter::once(&tx),
            "GRCh37",
        ));

        let synthesized = provider
            .get_transcript_on_build("NM_TEST.1", Some("GRCh37"))
            .expect("cdot exon-alignment fallback should fetch transcript");
        assert_eq!(
            synthesized.genome_build,
            GenomeBuild::GRCh37,
            "explicit build_hint must reach Transcript.genome_build"
        );
    }

    /// When the caller asks for a build that cdot did NOT load, the
    /// synthesize fallback must surface `ReferenceNotFound` rather than
    /// quietly synthesizing from the primary-build alignment and
    /// mis-tagging the result. Pre-#389 the fallback used the
    /// build-agnostic `cdot.get_transcript` which would silently succeed.
    #[test]
    fn test_synthesize_rejects_when_cdot_lacks_requested_build() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts_with_build(
            std::iter::once(&tx),
            "GRCh37",
        ));

        let err = provider
            .get_transcript_on_build("NM_TEST.1", Some("GRCh38"))
            .expect_err("cdot has no GRCh38 alignment; synthesize must fail honestly");
        assert!(
            matches!(err, FerroError::ReferenceNotFound { .. }),
            "expected ReferenceNotFound, got {err:?}"
        );
    }

    /// Main FASTA-backed path (NM_TEST.1 *is* in the transcript FASTA),
    /// cdot primary = GRCh37, no caller build hint → resulting
    /// `Transcript.genome_build` must be GRCh37, not the historical default.
    #[test]
    fn test_main_path_falls_back_to_cdot_primary_build_grch37() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::transcript::{GenomeBuild, Strand};

        let (mut provider, _kept) = build_provider_with_genome_and_named_tx("NM_TEST.1");
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts_with_build(
            std::iter::once(&tx),
            "GRCh37",
        ));

        let loaded = provider
            .get_transcript("NM_TEST.1")
            .expect("main path should fetch transcript from FASTA + cdot metadata");
        assert_eq!(
            loaded.genome_build,
            GenomeBuild::GRCh37,
            "main path must inherit cdot primary build when no caller hint is given"
        );
    }

    /// Main path with no cdot mapper at all → fall back to the historical
    /// default (GRCh38). Documents the unchanged behavior for FASTA-only
    /// providers.
    #[test]
    fn test_main_path_no_cdot_defaults_to_grch38() {
        use crate::reference::transcript::GenomeBuild;

        let (provider, _kept) = build_provider_with_genome_and_named_tx("NM_TEST.1");
        let loaded = provider
            .get_transcript("NM_TEST.1")
            .expect("main path should fetch transcript from FASTA (no cdot needed)");
        assert_eq!(
            loaded.genome_build,
            GenomeBuild::GRCh38,
            "FASTA-only provider with no cdot must preserve historical GRCh38 default"
        );
    }

    /// Main FASTA-backed path with an explicit `build_hint = Some("GRCh38")`
    /// against a cdot whose primary build is GRCh37 → the hint must win
    /// over the cdot primary on the main path. Sibling to
    /// `test_main_path_falls_back_to_cdot_primary_build_grch37` which
    /// covers the `None` arm; together they pin both arms of the helper.
    #[test]
    fn test_main_path_explicit_hint_overrides_cdot_primary() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::transcript::{GenomeBuild, Strand};

        let (mut provider, _kept) = build_provider_with_genome_and_named_tx("NM_TEST.1");
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts_with_build(
            std::iter::once(&tx),
            "GRCh37",
        ));

        let loaded = provider
            .get_transcript_on_build("NM_TEST.1", Some("GRCh38"))
            .expect("main path should fetch transcript from FASTA + cdot metadata");
        assert_eq!(
            loaded.genome_build,
            GenomeBuild::GRCh38,
            "explicit build_hint must override the cdot primary build on the main path"
        );
    }

    /// Direct unit coverage of `genome_build_from_hint` for the
    /// unrecognized-hint, casing, and empty-string arms — pre-#389 these
    /// were inline `match` legs; centralizing made them shared, so test
    /// them in one place.
    #[test]
    fn test_genome_build_from_hint_edge_cases() {
        use crate::reference::transcript::GenomeBuild;

        // No cdot mapper: None falls through to the historical GRCh38 default.
        let (provider, _kept) = build_provider_with_test_genome();
        assert_eq!(
            provider.genome_build_from_hint(None),
            GenomeBuild::GRCh38,
            "None with no cdot must default to GRCh38"
        );
        // Recognized hints map to their enum variant regardless of cdot.
        assert_eq!(
            provider.genome_build_from_hint(Some("GRCh37")),
            GenomeBuild::GRCh37
        );
        assert_eq!(
            provider.genome_build_from_hint(Some("GRCh38")),
            GenomeBuild::GRCh38
        );
        // Non-canonical hints — empty string, lowercase, UCSC-style — all
        // map to Unknown rather than silently coercing to GRCh38.
        for hint in &["", "grch38", "GRCH38", "hg19", "hg38", "GRCh99"] {
            assert_eq!(
                provider.genome_build_from_hint(Some(hint)),
                GenomeBuild::Unknown,
                "non-canonical hint {:?} must map to Unknown",
                hint
            );
        }
    }

    /// Build a provider whose FASTA index contains both `NC_TEST.1`
    /// (genome) and the requested transcript accession (so the main path
    /// is taken, not the cdot synthesize fallback).
    fn build_provider_with_genome_and_named_tx(
        tx_id: &str,
    ) -> (MultiFastaProvider, tempfile::TempDir) {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("genome.fna");
        let fai_path = dir.path().join("genome.fna.fai");
        let tx_fasta_path = dir.path().join("tx.fna");
        let tx_fai_path = dir.path().join("tx.fna.fai");
        {
            let mut f = File::create(&fasta_path).unwrap();
            writeln!(f, ">NC_TEST.1").unwrap();
            writeln!(f, "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC").unwrap();
        }
        {
            let mut f = File::create(&fai_path).unwrap();
            writeln!(f, "NC_TEST.1\t40\t11\t40\t41").unwrap();
        }
        {
            let mut f = File::create(&tx_fasta_path).unwrap();
            writeln!(f, ">{}", tx_id).unwrap();
            writeln!(f, "AAAACCCCGGGGTTTTAAAA").unwrap();
        }
        {
            let mut f = File::create(&tx_fai_path).unwrap();
            // Header ">NM_TEST.1\n" -> header byte length = len(tx_id) + 2.
            let header_len = tx_id.len() as u64 + 2;
            writeln!(f, "{}\t20\t{}\t20\t21", tx_id, header_len).unwrap();
        }
        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();
        (provider, dir)
    }

    #[test]
    fn test_get_transcript_falls_back_errors_on_degenerate_exon_span() {
        use crate::data::cdot::CdotTranscript;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        // e[1] == e[0]: zero-length exon. Bypass from_transcripts (which would
        // reject this on input validation) by writing directly into CdotMapper.
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_DEGEN.1".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_TEST.1".to_string(),
                strand: Strand::Plus,
                exons: vec![[10, 10, 0, 0]],
                cds_start: None,
                cds_end: None,
                gene_id: None,
                protein: None,
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);

        let err = provider
            .get_transcript("NM_DEGEN.1")
            .expect_err("degenerate exon span must produce a typed error, not a panic");
        assert!(
            matches!(err, FerroError::InvalidCoordinates { .. }),
            "expected InvalidCoordinates, got {err:?}"
        );
    }

    // ----------------------------------------------------------------------
    // get_transcript_strict — exact-version pinning (#471 / #478 pillar 3)
    //
    // These tests pin the opt-in strict-resolution capability and, crucially,
    // assert that the lenient `get_transcript` chain is UNCHANGED: where strict
    // refuses a sibling-version substitution, lenient still falls back to it.
    // ----------------------------------------------------------------------

    #[test]
    fn build_base_to_versioned_picks_highest_version_deterministically() {
        let entry = |name: &str| FastaIndexEntry {
            file_path: PathBuf::from("x.fna"),
            name: name.to_string(),
            length: 1,
            offset: 0,
            line_bases: 1,
            line_bytes: 2,
        };
        let mut index: FxHashMap<String, FastaIndexEntry> = FxHashMap::default();
        for v in ["NM_000088.3", "NM_000088.4", "NM_000088.10", "NM_000088.2"] {
            index.insert(v.to_string(), entry(v));
        }
        let map = build_base_to_versioned(&index);
        // Highest *numeric* version wins (10 > 4 > 3 > 2) — not lexical, which
        // would pick ".4" over ".10" — and the result is independent of the
        // (hash) iteration order of the index.
        assert_eq!(
            map.get("NM_000088").map(String::as_str),
            Some("NM_000088.10")
        );
    }

    /// Build a provider from a transcript FASTA carrying exactly the given
    /// `(accession, sequence)` entries (one record each). No genome, no cdot.
    /// Returns the provider plus the kept `TempDir`.
    fn build_provider_with_transcripts(
        entries: &[(&str, &str)],
    ) -> (MultiFastaProvider, tempfile::TempDir) {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("tx.fna");
        let fai_path = dir.path().join("tx.fna.fai");
        let mut fasta = File::create(&fasta_path).unwrap();
        let mut fai = File::create(&fai_path).unwrap();
        let mut offset: u64 = 0;
        for (acc, seq) in entries {
            // Header line ">{acc}\n" then "{seq}\n".
            let header_len = acc.len() as u64 + 2; // '>' + acc + '\n'
            let seq_offset = offset + header_len;
            writeln!(fasta, ">{}", acc).unwrap();
            writeln!(fasta, "{}", seq).unwrap();
            // name<TAB>length<TAB>offset<TAB>line_bases<TAB>line_bytes
            let len = seq.len() as u64;
            writeln!(
                fai,
                "{}\t{}\t{}\t{}\t{}",
                acc,
                len,
                seq_offset,
                len,
                len + 1
            )
            .unwrap();
            // Advance past this record: header + seq + newline.
            offset = seq_offset + len + 1;
        }
        drop(fasta);
        drop(fai);
        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();
        (provider, dir)
    }

    #[test]
    fn resolve_cache_memoizes_and_preserves_content() {
        let (provider, _kept) = build_provider_with_transcripts(&[
            ("NM_CACHE.1", "ATGAAATAA"),
            ("NM_OTHER.1", "ATGCAT"),
        ]);

        // The cached result must be byte-identical to an uncached resolve:
        // caching is purely a timing change, never a behavioral one.
        let uncached = provider
            .get_transcript_on_build("NM_CACHE.1", None)
            .unwrap();
        let first = provider.get_transcript("NM_CACHE.1").unwrap();
        assert_eq!(*first, uncached);

        // A second lookup of the same (id, build) returns the *same* Arc — proof
        // the materialize path ran once and the repeat was served from cache.
        let second = provider.get_transcript("NM_CACHE.1").unwrap();
        assert!(
            Arc::ptr_eq(&first, &second),
            "repeat lookup should hit the cache and return the same Arc"
        );

        // A different accession is a distinct entry, not an accidental alias.
        let other = provider.get_transcript("NM_OTHER.1").unwrap();
        assert!(!Arc::ptr_eq(&first, &other));
        assert_eq!(other.sequence.as_deref(), Some("ATGCAT"));

        // A miss still surfaces as an error (not a stale or empty cache hit).
        assert!(provider.get_transcript("NM_MISSING.1").is_err());
    }

    #[test]
    fn resolve_cache_build_hint_key_isolation() {
        // The docstring for get_transcript_on_build_cached states that the
        // cache key is (id, build_hint), so the same accession with different
        // build hints must produce independent entries that never evict each
        // other.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_BUILD.1", "ATGAAACCC")]);

        let grch37 = provider
            .get_transcript_on_build_cached("NM_BUILD.1", Some("GRCh37"))
            .expect("GRCh37 lookup must succeed");
        let grch38 = provider
            .get_transcript_on_build_cached("NM_BUILD.1", Some("GRCh38"))
            .expect("GRCh38 lookup must succeed");

        // Different build hints produce distinct Arc allocations — neither
        // entry evicted the other from the cache.
        assert!(
            !Arc::ptr_eq(&grch37, &grch38),
            "(id, GRCh37) and (id, GRCh38) must be independent cache entries"
        );

        // A repeat GRCh37 lookup is still served from the cache (same Arc).
        let grch37_again = provider
            .get_transcript_on_build_cached("NM_BUILD.1", Some("GRCh37"))
            .expect("repeat GRCh37 lookup must succeed");
        assert!(
            Arc::ptr_eq(&grch37, &grch37_again),
            "repeat (id, GRCh37) lookup should hit the cache and return the same Arc"
        );
    }

    #[test]
    fn test_get_transcript_strict_returns_exact_version() {
        // FASTA carries the exact requested version → strict returns Ok with
        // the right sequence, identical to the lenient path.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_212556.2", "ACGTACGTAC")]);

        let strict = provider
            .get_transcript_strict("NM_212556.2")
            .expect("exact version present → strict must succeed");
        assert_eq!(strict.id, "NM_212556.2");
        assert_eq!(strict.sequence.as_deref(), Some("ACGTACGTAC"));

        // has_transcript_exact agrees.
        assert!(provider.has_transcript_exact("NM_212556.2"));

        // Strict result matches the lenient result for an exact hit.
        let lenient = provider.get_transcript("NM_212556.2").unwrap();
        assert_eq!(strict.id, lenient.id);
        assert_eq!(strict.sequence, lenient.sequence);
    }

    #[test]
    fn test_get_transcript_strict_rejects_sibling_version_but_lenient_falls_back() {
        // FASTA ships only the .4 sibling; the request is for .2.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_212556.4", "TTTTGGGGCC")]);

        // Strict: structured error, no substitution.
        let err = provider
            .get_transcript_strict("NM_212556.2")
            .expect_err("only sibling .4 present → strict must refuse .2");
        match err {
            FerroError::TranscriptVersionNotExact { requested } => {
                assert_eq!(requested, "NM_212556.2");
            }
            other => panic!("expected TranscriptVersionNotExact, got {other:?}"),
        }
        assert!(!provider.has_transcript_exact("NM_212556.2"));
        assert!(provider.has_transcript_exact("NM_212556.4"));

        // Lenient: UNCHANGED — still falls back to the .4 sibling and serves
        // its bases. This is the behavior strict exists to opt out of.
        let lenient = provider
            .get_transcript("NM_212556.2")
            .expect("lenient must still fall back to the sibling version");
        assert_eq!(lenient.id, "NM_212556.4");
        assert_eq!(lenient.sequence.as_deref(), Some("TTTTGGGGCC"));
    }

    #[test]
    fn test_get_transcript_strict_refuses_cdot_synthesis_but_lenient_synthesizes() {
        use crate::reference::transcript::Strand;

        // Genome FASTA only (NC_TEST.1); NM_TEST.1 is NOT in any transcript
        // FASTA. cdot carries the exon alignment, so the lenient path
        // synthesizes bases from the genome — strict must refuse that
        // reconstruction.
        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));

        // Strict: NM_TEST.1 is absent from the FASTA index, so refuse even
        // though cdot could reconstruct it.
        let err = provider
            .get_transcript_strict("NM_TEST.1")
            .expect_err("cdot-genome reconstruction must be refused under strict resolution");
        assert!(
            matches!(err, FerroError::TranscriptVersionNotExact { .. }),
            "expected TranscriptVersionNotExact, got {err:?}"
        );
        assert!(!provider.has_transcript_exact("NM_TEST.1"));

        // Lenient: UNCHANGED — still synthesizes the transcript from cdot.
        let synthesized = provider
            .get_transcript("NM_TEST.1")
            .expect("lenient must still synthesize from cdot exon alignment");
        assert_eq!(
            synthesized.sequence.as_deref(),
            Some("GGTTTTAAAACCCCGGGGTT")
        );
    }

    #[test]
    fn has_transcript_version_exact_true_for_exact_fasta() {
        // The exact requested version is shipped in the FASTA → version-exact.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_212556.2", "ACGTACGTAC")]);
        assert!(provider.has_transcript_version_exact("NM_212556.2"));
    }

    #[test]
    fn has_transcript_version_exact_false_for_sibling_only() {
        // Only the .4 sibling is shipped; a request for .2 would be served by
        // the lenient version-strip substitution. The protein gate (#505) must
        // see this as NOT version-exact so it declines to translate the wrong
        // bases — even though the lenient `get_transcript` still substitutes.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_212556.4", "TTTTGGGGCC")]);
        assert!(!provider.has_transcript_version_exact("NM_212556.2"));
        // The exact sibling itself is version-exact.
        assert!(provider.has_transcript_version_exact("NM_212556.4"));
    }

    #[test]
    fn has_transcript_version_exact_true_for_exact_cdot_synthesis() {
        use crate::reference::transcript::Strand;
        // NM_TEST.1 is absent from the transcript FASTA but cdot carries it at
        // the exact version, so the lenient path reconstructs its bases from
        // the genome. Those bases DO correspond to the requested version, so it
        // is version-exact and the gate must permit protein prediction — unlike
        // strict resolution, which refuses all synthesis. This pins the #331
        // (exact-version synthesis) regression guard called out in the #505
        // design review.
        let (mut provider, _kept) = build_provider_with_test_genome();
        let tx = build_single_exon_synthetic_tx(Strand::Plus);
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));
        // Not in the FASTA index...
        assert!(!provider.has_transcript_exact("NM_TEST.1"));
        // ...but cdot carries the exact version, so it is version-exact.
        assert!(provider.has_transcript_version_exact("NM_TEST.1"));
    }

    #[test]
    fn has_transcript_version_exact_false_when_fasta_sibling_shadows_cdot_exact() {
        use crate::reference::transcript::Strand;
        // The FASTA ships only the .2 sibling; cdot carries the exact .1. The
        // read path (`get_transcript_on_build`) resolves the FASTA sibling via
        // `resolve_name`'s version-strip BEFORE it would ever reach cdot
        // synthesis, so a request for .1 is actually served the .2 bases. The
        // gate must therefore report .1 as NOT version-exact even though cdot
        // has an exact .1 record — otherwise it green-lights translating .2
        // bases as .1, the exact wrong-protein attribution #505 targets.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_TEST.2", "ACGTACGTAC")]);
        let tx = build_single_exon_synthetic_tx(Strand::Plus); // id NM_TEST.1
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));
        // Precondition: cdot really does carry the exact .1 record...
        assert!(provider
            .cdot_mapper()
            .unwrap()
            .has_transcript_exact("NM_TEST.1"));
        // ...but a FASTA sibling (.2) shadows it on the read path, so .1 is not
        // reachable at its exact version.
        assert!(!provider.has_transcript_version_exact("NM_TEST.1"));
    }

    // --- canonical-override wiring (issue #520, edge 2) ---------------------

    #[test]
    fn get_transcript_applies_canonical_override() {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        // FASTA-only transcript (9 nt); the override supplies the canonical CDS
        // + protein. get_transcript must return the corrected metadata.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_OV.2".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;

        let tx = provider.get_transcript("NM_OV.2").unwrap();
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_OV.2"));
    }

    #[test]
    fn reconcile_cdot_converts_authoritative_to_cdot_coords() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        use crate::reference::Strand;

        // FASTA transcript (9 nt) + a cdot record with the WRONG CDS/protein.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_OV.2".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "chr1".to_string(),
                strand: Strand::Plus,
                exons: vec![[0, 9, 0, 9]],
                cds_start: Some(2), // wrong
                cds_end: Some(6),   // wrong
                gene_id: None,
                protein: Some("NP_WRONG.9".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);

        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_OV.2".to_string(),
            tx_length: 9,
            cds_start: Some(1), // authoritative 1-based
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;

        provider.reconcile_cdot_with_overrides();

        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_OV.2")
            .unwrap();
        // 1-based start 1 → cdot 0-based 0; end 9 unchanged (0-based exclusive).
        assert_eq!(rec.cds_start, Some(0));
        assert_eq!(rec.cds_end, Some(9));
        assert_eq!(rec.protein.as_deref(), Some("NP_OV.2"));
    }

    #[test]
    fn reconcile_cdot_skips_on_length_mismatch() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        use crate::reference::Strand;
        // Served FASTA is 9 nt but the authoritative length is 99 → the
        // coordinates don't apply; cdot must be left untouched.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_OV.2".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "chr1".to_string(),
                strand: Strand::Plus,
                exons: vec![[0, 9, 0, 9]],
                cds_start: Some(2),
                cds_end: Some(6),
                gene_id: None,
                protein: Some("NP_WRONG.9".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_OV.2".to_string(),
            tx_length: 99, // mismatch
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;

        provider.reconcile_cdot_with_overrides();

        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_OV.2")
            .unwrap();
        assert_eq!(rec.cds_start, Some(2), "length mismatch → cdot untouched");
        assert_eq!(rec.protein.as_deref(), Some("NP_WRONG.9"));
    }

    #[test]
    fn reconcile_cdot_requires_exact_version_not_sibling() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        use crate::reference::Strand;
        // FASTA carries NM_OV.2 (9 nt). cdot has only the .4 SIBLING (different
        // contig/exons). An override for .2 must NOT clone the .4 record under
        // .2: `get_transcript` version-falls-back to .4, but reconcile requires
        // an *exact* .2 cdot record, so it leaves cdot untouched.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_OV.4".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "chrSIBLING".to_string(),
                strand: Strand::Plus,
                exons: vec![[0, 9, 0, 9]],
                cds_start: Some(2),
                cds_end: Some(6),
                gene_id: None,
                protein: Some("NP_SIB.4".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_OV.2".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;

        provider.reconcile_cdot_with_overrides();

        // No exact .2 in cdot → nothing written under .2.
        assert!(
            !provider
                .cdot_mapper()
                .unwrap()
                .has_transcript_exact("NM_OV.2"),
            "a sibling-only cdot record must not be cloned under the exact version"
        );
        // The .4 sibling is left untouched.
        let sib = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_OV.4")
            .unwrap();
        assert_eq!(sib.contig, "chrSIBLING");
        assert_eq!(sib.protein.as_deref(), Some("NP_SIB.4"));
    }

    #[test]
    fn get_transcript_strict_also_applies_canonical_override() {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        // The exact version is in the FASTA index, so strict succeeds; the
        // override must still be applied (orthogonal to version-exactness).
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_OV.2".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;
        let tx = provider.get_transcript_strict("NM_OV.2").unwrap();
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_OV.2"));
    }

    #[test]
    fn get_transcript_for_variant_applies_canonical_override() {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        // A bare-NM_ coding variant routes through get_transcript_for_variant's
        // no-parent branch → the wrapper → override applied.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_OV.2".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;
        let variant = crate::parse_hgvs("NM_OV.2:c.4C>A").unwrap();
        let tx = provider.get_transcript_for_variant(&variant).unwrap();
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_OV.2"));
    }

    // --- protein FASTA ingestion (issue #520, edge 3) -----------------------

    #[test]
    fn get_protein_sequence_reads_from_protein_index() {
        // build a FASTA entry, then move it into the protein index (a real file
        // on disk backs it, so get_sequence_from_index can read it).
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NP_TEST.1", "MWAPLE")]);
        let entry = provider.index.get("NP_TEST.1").unwrap().clone();
        provider
            .protein_index
            .insert("NP_TEST.1".to_string(), entry);

        assert!(provider.has_protein_data());
        // `0, u64::MAX` fetches the full protein (read clamps end to length).
        assert_eq!(
            provider
                .get_protein_sequence("NP_TEST.1", 0, u64::MAX)
                .unwrap(),
            "MWAPLE"
        );
        assert!(provider.get_protein_sequence("NP_NOPE.9", 0, 10).is_err());
    }

    #[test]
    fn validate_translations_flags_mismatch_end_to_end() {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        use crate::reference::validate::AnomalyKind;
        // NM_T.1 = ATG|TGG|TAA → "MW"; the canonical protein NP_T.1 = "MV".
        let (mut provider, _kept) =
            build_provider_with_transcripts(&[("NM_T.1", "ATGTGGTAA"), ("NP_T.1", "MV")]);
        let np = provider.index.get("NP_T.1").unwrap().clone();
        provider.protein_index.insert("NP_T.1".to_string(), np);
        // The override gives NM_T.1 its CDS so the served transcript is coding.
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_T.1".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_T.1".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;

        let found = provider.validate_translations(&[("NM_T.1".to_string(), "NP_T.1".to_string())]);
        assert_eq!(found.len(), 1, "translated MW vs canonical MV → mismatch");
        assert_eq!(found[0].kind, AnomalyKind::TranslationMismatch);
    }

    #[test]
    fn validate_translations_skips_non_version_exact_sibling() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::Strand;
        // FASTA ships only the .4 sibling (bases translate "MW"); a request for
        // .2 is served the .4 bases via resolve_name's version-strip, so .2 is
        // NOT version-exact. cdot gives .4 a CDS so the served transcript is
        // coding and *would* mismatch the canonical "MV" — but translation
        // validation must SKIP a non-version-exact accession rather than
        // attribute the sibling's bases to the requested version (#505).
        let (mut provider, _kept) =
            build_provider_with_transcripts(&[("NM_T.4", "ATGTGGTAA"), ("NP_T.1", "MV")]);
        let np = provider.index.get("NP_T.1").unwrap().clone();
        provider.protein_index.insert("NP_T.1".to_string(), np);
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_T.4".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_X.1".to_string(),
                strand: Strand::Plus,
                exons: vec![[0, 9, 1, 9]],
                cds_start: Some(0), // 0-based → served 1-based 1
                cds_end: Some(9),
                gene_id: None,
                protein: Some("NP_T.1".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);

        // Precondition: the .4 sibling shadows .2, so .2 is not version-exact.
        assert!(!provider.has_transcript_version_exact("NM_T.2"));
        // Control: served .2 bases would translate to a mismatch if validated.
        let served = provider.get_transcript("NM_T.2").unwrap();
        assert_eq!(served.cds_start, Some(1), "sibling served as coding");
        // The non-version-exact request must be skipped → no anomaly.
        let found = provider.validate_translations(&[("NM_T.2".to_string(), "NP_T.1".to_string())]);
        assert!(
            found.is_empty(),
            "non-version-exact sibling must be skipped, got {found:?}"
        );
    }

    #[test]
    fn reconcile_cdot_uses_exon_sum_for_synthesis_only_transcript() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        use crate::reference::Strand;
        // NM_SYNTH.1 has NO FASTA entry (only the genome contig is indexed); cdot
        // carries it with exon tx_end=9 → exon-sum length 9, matching the
        // authoritative length. With no carried `sequence`, reconciliation must
        // still fire via the cdot exon-sum gate.
        let (mut provider, _kept) = build_provider_with_test_genome();
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_SYNTH.1".to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_TEST.1".to_string(),
                strand: Strand::Plus,
                exons: vec![[0, 9, 0, 9]],
                cds_start: Some(2), // wrong
                cds_end: Some(6),   // wrong
                gene_id: None,
                protein: Some("NP_WRONG.9".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_SYNTH.1".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None, // no carried sequence → exon-sum gate must fire
        });
        provider.canonical_overrides = ov;

        provider.reconcile_cdot_with_overrides();

        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_SYNTH.1")
            .unwrap();
        assert_eq!(rec.cds_start, Some(0)); // reconciled (1-based 1 → 0-based 0)
        assert_eq!(rec.cds_end, Some(9));
        assert_eq!(rec.protein.as_deref(), Some("NP_OV.2"));
    }

    /// Build a synthesis-only provider (genome contig only, no transcript FASTA)
    /// with one cdot record for `acc` given its exons + (wrong) CDS/protein.
    #[cfg(test)]
    fn synth_provider_with_cdot(
        acc: &str,
        exons: Vec<[u64; 4]>,
    ) -> (MultiFastaProvider, tempfile::TempDir) {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::Strand;
        let (mut provider, kept) = build_provider_with_test_genome();
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            acc.to_string(),
            CdotTranscript {
                gene_name: None,
                contig: "NC_TEST.1".to_string(),
                strand: Strand::Plus,
                exons,
                cds_start: Some(2), // wrong
                cds_end: Some(6),   // wrong
                gene_id: None,
                protein: Some("NP_WRONG.9".to_string()),
                exon_cigars: Vec::new(),
            },
        );
        provider.cdot_mapper = Some(cdot);
        (provider, kept)
    }

    fn coding_override(
        acc: &str,
        tx_length: u64,
    ) -> crate::reference::authoritative::CanonicalOverrides {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: acc.to_string(),
            tx_length,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_OV.2".to_string()),
            sequence: None,
        });
        ov
    }

    #[test]
    fn reconcile_cdot_synthesis_only_multi_exon_uses_genomic_span() {
        // Two exons with genomic spans 5 (0..5) and 4 (10..14) → summed genomic
        // span 9 == authoritative length → reconcile. (Here the genomic span and
        // raw max tx_end coincide; the next test pins the case where they don't.)
        let (mut provider, _kept) =
            synth_provider_with_cdot("NM_MX.1", vec![[0, 5, 0, 5], [10, 14, 5, 9]]);
        provider.canonical_overrides = coding_override("NM_MX.1", 9);
        provider.reconcile_cdot_with_overrides();
        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_MX.1")
            .unwrap();
        assert_eq!(rec.cds_start, Some(0));
        assert_eq!(rec.protein.as_deref(), Some("NP_OV.2"));
    }

    #[test]
    fn reconcile_cdot_synthesis_only_uses_served_length_not_raw_tx_end() {
        // A cdot exon whose genomic span (12 bases, 0..12) exceeds its transcript
        // extent (tx_end 9) — the served bases come from the genomic span, so
        // `synthesize_transcript_from_cdot` produces a 12 nt sequence. The
        // authoritative length 12 therefore matches the *served* length, and
        // reconciliation must fire. Keying off the raw max tx_end (9) instead
        // would wrongly skip, gating on a length the server never produced.
        let (mut provider, _kept) = synth_provider_with_cdot("NM_GSPAN.1", vec![[0, 12, 0, 9]]);
        provider.canonical_overrides = coding_override("NM_GSPAN.1", 12);
        provider.reconcile_cdot_with_overrides();
        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_GSPAN.1")
            .unwrap();
        assert_eq!(
            rec.cds_start,
            Some(0),
            "served (genomic-span) length 12 matches the override → reconcile"
        );
        assert_eq!(rec.protein.as_deref(), Some("NP_OV.2"));
    }

    #[test]
    fn reconcile_cdot_synthesis_only_length_mismatch_skips() {
        // cdot exon max tx_end 9 but authoritative length 100, no sequence → skip.
        let (mut provider, _kept) = synth_provider_with_cdot("NM_MM.1", vec![[0, 9, 0, 9]]);
        provider.canonical_overrides = coding_override("NM_MM.1", 100);
        provider.reconcile_cdot_with_overrides();
        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_MM.1")
            .unwrap();
        assert_eq!(rec.cds_start, Some(2), "length mismatch → cdot untouched");
        assert_eq!(rec.protein.as_deref(), Some("NP_WRONG.9"));
    }

    #[test]
    fn reconcile_cdot_synthesis_only_empty_exons_skips() {
        // No exons → max() is None → no length match → skip (no panic).
        let (mut provider, _kept) = synth_provider_with_cdot("NM_EX.1", vec![]);
        provider.canonical_overrides = coding_override("NM_EX.1", 9);
        provider.reconcile_cdot_with_overrides();
        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_EX.1")
            .unwrap();
        assert_eq!(rec.cds_start, Some(2), "empty exons → cdot untouched");
    }

    // ----------------------------------------------------------------------
    // from_manifest: the `cdot_grch37_json` key loads a GRCh37 secondary build.
    // Exercised via the factored-out `load_grch37_secondary_cdot` helper so the
    // test needs only a tiny GRCh37 cdot file + an in-memory manifest, not a
    // full hermetic manifest with genome/transcript FASTAs.
    // ----------------------------------------------------------------------

    fn grch37_cdot_json() -> &'static str {
        r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.10",
                    "strand": "+",
                    "exons": [[48263025, 48263098, 0, 73, "M73"]],
                    "start_codon": 10,
                    "stop_codon": 60
                }
            }
        }
        "#
    }

    #[test]
    fn test_load_grch37_secondary_cdot_loads_when_key_present() {
        let dir = tempdir().unwrap();
        let grch37_path = dir.path().join("cdot-grch37.json");
        std::fs::write(&grch37_path, grch37_cdot_json()).unwrap();

        // Primary mapper = GRCh38 (only the .11 contig known initially).
        let grch38_json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73, "M73"]],
                    "start_codon": 10,
                    "stop_codon": 60
                }
            }
        }
        "#;
        let mut mapper = CdotMapper::from_reader(grch38_json.as_bytes()).unwrap();
        assert_eq!(
            mapper.available_builds_for("NM_000088.3"),
            vec!["GRCh38".to_string()],
            "only GRCh38 before loading the secondary build"
        );

        // Manifest references the GRCh37 cdot by absolute path.
        let manifest = serde_json::json!({
            "cdot_grch37_json": grch37_path.to_str().unwrap(),
        });
        let resolve_path = |p: &str| -> PathBuf { PathBuf::from(p) };

        let loaded =
            MultiFastaProvider::load_grch37_secondary_cdot(&mut mapper, &manifest, &resolve_path);
        assert_eq!(loaded, 1, "one GRCh37 transcript loaded from the manifest");

        let mut builds = mapper.available_builds_for("NM_000088.3");
        builds.sort();
        assert_eq!(builds, vec!["GRCh37".to_string(), "GRCh38".to_string()]);
        // Overlap discovery now routes the GRCh37 contig to the loaded build.
        let hits = mapper.transcripts_at_position("NC_000017.10", 48263050);
        assert_eq!(
            hits.iter().map(|(acc, _)| *acc).collect::<Vec<_>>(),
            vec!["NM_000088.3"]
        );
    }

    #[test]
    fn test_load_grch37_secondary_cdot_noop_when_key_absent() {
        let grch38_json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73, "M73"]]
                }
            }
        }
        "#;
        let mut mapper = CdotMapper::from_reader(grch38_json.as_bytes()).unwrap();
        // Manifest without the GRCh37 key → helper is a no-op.
        let manifest = serde_json::json!({ "cdot_json": "ignored.json" });
        let resolve_path = |p: &str| -> PathBuf { PathBuf::from(p) };
        let loaded =
            MultiFastaProvider::load_grch37_secondary_cdot(&mut mapper, &manifest, &resolve_path);
        assert_eq!(loaded, 0);
        assert_eq!(
            mapper.available_builds_for("NM_000088.3"),
            vec!["GRCh38".to_string()]
        );
    }
}
