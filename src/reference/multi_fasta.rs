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
use std::io::{BufRead, BufReader};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use log::warn;
use lru::LruCache;
use rustc_hash::FxHashMap;

/// Read exactly `buf.len()` bytes from `file` starting at byte `offset`, using a
/// positioned read (no shared file cursor, no `seek`), so a cached read-only
/// handle can be shared across calls and threads. Implemented per-platform:
/// `pread` on Unix, `seek_read` on Windows. On other targets (e.g. WASM, WASI,
/// Redox) the implementation clones the file descriptor and performs a
/// `seek` + `read_exact` on the clone so the caller's handle cursor is
/// unaffected and the cached handle remains shareable.
#[cfg(unix)]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::os::unix::fs::FileExt;
    file.read_exact_at(buf, offset)
}

#[cfg(windows)]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::os::windows::fs::FileExt;
    let mut filled = 0usize;
    while filled < buf.len() {
        let n = file.seek_read(&mut buf[filled..], offset + filled as u64)?;
        if n == 0 {
            return Err(std::io::Error::from(std::io::ErrorKind::UnexpectedEof));
        }
        filled += n;
    }
    Ok(())
}

/// Fallback for targets that are neither Unix nor Windows (e.g. WASM, WASI, Redox).
/// Clones the file descriptor so the seek does not disturb the cached handle's cursor.
#[cfg(not(any(unix, windows)))]
fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::io::{Read, Seek, SeekFrom};
    let mut cloned = file.try_clone()?;
    cloned.seek(SeekFrom::Start(offset))?;
    cloned.read_exact(buf)
}

use crate::data::cdot::CdotMapper;
use crate::error::FerroError;
use crate::reference::authoritative::CanonicalOverrides;
use crate::reference::provider::{GenomicPlacement, ReferenceProvider};
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

/// Whether `index` holds any genomic (chromosome) sequence, named either
/// NCBI RefSeq style (`NC_*`) or UCSC style (`chr*`). Computed once at
/// construction and cached in [`MultiFastaProvider::has_genomic_data`] — the
/// index can hold hundreds of thousands of keys and this predicate is hit on
/// the per-variant normalization path.
fn index_has_genomic_data(index: &FxHashMap<String, FastaIndexEntry>) -> bool {
    index
        .keys()
        .any(|k| k.starts_with("NC_") || k.starts_with("chr"))
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
    /// Cache of open FASTA file handles keyed by path. Sequence reads use a
    /// positioned `read_exact_at` on a reused handle (one syscall) instead of
    /// reopening the file on every read — `File::open` per read was ~10 us and
    /// dominated `get_sequence`/transcript resolution. There are only a handful
    /// of distinct FASTA files (genome, transcripts, LRG), so this is bounded.
    open_files: Mutex<FxHashMap<PathBuf, Arc<File>>>,
    /// Whether the index contains any genomic (chromosome) sequences, computed
    /// once at construction. `has_genomic_data` is queried on the per-variant
    /// normalization path; computing it by scanning all ~270k index keys on
    /// every call dominated the hot loop (~12% of normalize CPU). The index is
    /// immutable after construction, so the answer is a build-time constant.
    has_genomic_data: bool,
    /// Directory holding the LRG XML records (`LRG_<n>.xml`), if the manifest
    /// listed `lrg_xmls`. Used to resolve an LRG parent's chromosomal placement
    /// on demand for #480 (the XML's `main_assembly` `<mapping>` element gives
    /// LRG↔NC_ coordinates), so transcript coordinates projected onto an LRG
    /// parent can be re-expressed in the LRG's own frame.
    lrg_xml_dir: Option<PathBuf>,
    /// Memoizes parsed LRG placements (and negative results) keyed by LRG
    /// accession, so repeated projections against the same LRG parent parse the
    /// XML once. A stored `None` means "parsed, but no usable single-span
    /// main-assembly placement".
    lrg_placement_cache: Mutex<FxHashMap<String, Option<GenomicPlacement>>>,
    /// NG_ RefSeqGene → chromosome placements, parsed once from the NCBI
    /// `GCF_*_refseqgene_alignments.gff3` named by the manifest's
    /// `refseqgene_alignments`, keyed by NG accession.version (#480). Empty when
    /// the manifest carries no alignments file.
    /// Per-NG-version placements, one entry per genome build (GRCh37/GRCh38)
    /// the parsed RefSeqGene alignments cover (#653). Selection by the input's
    /// resolved build happens in `genomic_placement_on_build`.
    refseqgene_placements: FxHashMap<String, Vec<GenomicPlacement>>,
    /// Legacy gene-model selector resolution: gene Symbol (upper-case) →
    /// reference-standard transcript accession, parsed once from NCBI's
    /// `LRG_RefSeqGene` table named by the manifest's `refseqgene_summary`
    /// (#500/#637). Empty when the manifest carries no summary file.
    legacy_gene_models: FxHashMap<String, String>,
    /// Data-driven accession→build table from the prepared reference's own
    /// assembly report(s) (#716). `None` when the manifest names no report, in
    /// which case build inference uses only the hardcoded fallback.
    contig_aliases: Option<crate::liftover::aliases::ContigAliases>,
    /// Per-`NG_`-version hosted-transcript map (#792); `None` when the manifest
    /// names no artifact (resolution falls back to the global reference map).
    ng_hosted: Option<crate::reference::ng_hosted_transcripts::NgHostedTranscripts>,
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

            // Skip hidden files and macOS AppleDouble sidecars (names beginning with '.', e.g.
            // `._LRG_430.fasta`). These are not real FASTAs — an accompanying `._*.fasta.fai` is a
            // binary extended-attribute blob, and attempting to parse it as a FAI index would abort
            // the entire reference load with a "stream did not contain valid UTF-8" error. Such files
            // are commonly introduced when a reference directory is copied through macOS tooling.
            let is_hidden = path
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with('.'));
            if is_hidden {
                continue;
            }

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
            has_genomic_data: index_has_genomic_data(&index),
            index,
            aliases: build_chromosome_aliases(),
            cdot_mapper: None,
            supplemental_cds: SupplementalCdsInfo::default(),
            canonical_overrides: CanonicalOverrides::default(),
            protein_index: FxHashMap::default(),
            transcript_cache: new_transcript_cache(),
            open_files: Mutex::new(FxHashMap::default()),
            lrg_xml_dir: None,
            lrg_placement_cache: Mutex::new(FxHashMap::default()),
            refseqgene_placements: FxHashMap::default(),
            legacy_gene_models: FxHashMap::default(),
            contig_aliases: None,
            ng_hosted: None,
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
            has_genomic_data: index_has_genomic_data(&combined_index),
            index: combined_index,
            aliases: build_chromosome_aliases(),
            cdot_mapper: None,
            supplemental_cds: SupplementalCdsInfo::default(),
            canonical_overrides: CanonicalOverrides::default(),
            protein_index: FxHashMap::default(),
            transcript_cache: new_transcript_cache(),
            open_files: Mutex::new(FxHashMap::default()),
            lrg_xml_dir: None,
            lrg_placement_cache: Mutex::new(FxHashMap::default()),
            refseqgene_placements: FxHashMap::default(),
            legacy_gene_models: FxHashMap::default(),
            contig_aliases: None,
            ng_hosted: None,
        })
    }

    /// Defer loading the GRCh37 cdot referenced by `cdot_grch37_json` (if present)
    /// as a secondary build on `mapper`. Records the path; the ~198k-transcript
    /// load happens lazily on the first GRCh37 lookup (see
    /// `CdotMapper::defer_secondary_build`). Best-effort: a missing key/file is a
    /// no-op, like the eager loader it replaces.
    fn defer_grch37_secondary_cdot(
        mapper: &mut CdotMapper,
        manifest: &serde_json::Value,
        resolve_path: &impl Fn(&str) -> PathBuf,
    ) {
        let Some(grch37_path_str) = manifest.get("cdot_grch37_json").and_then(|v| v.as_str())
        else {
            return;
        };
        let grch37_path = resolve_path(grch37_path_str);
        if !grch37_path.exists() {
            eprintln!(
                "Warning: cdot_grch37_json path does not exist: {}",
                grch37_path.display()
            );
            return;
        }
        eprintln!(
            "Deferring GRCh37 cdot secondary build (loads on first GRCh37 use): {}",
            grch37_path.display()
        );
        mapper.defer_secondary_build("GRCh37", grch37_path);
    }

    /// Wire Ensembl cdot metadata into an already-loaded RefSeq `mapper`,
    /// **deferring** the GRCh38 primary merge (#964). Ensembl transcripts
    /// (`ENST`/`ENSG`/`ENSP`) share the primary genome build, but merging their
    /// ~198k records eagerly at construction made every startup — including a
    /// pure-RefSeq `ferro normalize` — pay for transcripts it may never touch.
    /// Instead the primary Ensembl cdot is registered as a *deferred* source
    /// ([`CdotMapper::defer_ensembl_primary_merge`]) that materializes lazily on
    /// the first `ENS*`-accession lookup; a RefSeq workload never triggers it
    /// (genomic enumeration is RefSeq-scoped and does not consult it). An
    /// optional Ensembl GRCh37 cdot is still
    /// layered eagerly as a secondary build (present only when the manifest wires
    /// it, which the GRCh38 primary reference does not). Best-effort: warns but
    /// does not fail, like the LRG/GRCh37 loads.
    fn defer_ensembl_cdot(
        mapper: &mut CdotMapper,
        manifest: &serde_json::Value,
        resolve_path: &impl Fn(&str) -> PathBuf,
    ) {
        if let Some(ensembl_path_str) = manifest.get("ensembl_cdot_json").and_then(|v| v.as_str()) {
            let ensembl_path = resolve_path(ensembl_path_str);
            if ensembl_path.exists() {
                eprintln!(
                    "Deferring Ensembl cdot primary merge (loads on first Ensembl use): {}",
                    ensembl_path.display()
                );
                mapper.defer_ensembl_primary_merge(ensembl_path);
            } else {
                eprintln!(
                    "Warning: ensembl_cdot_json path does not exist: {}",
                    ensembl_path.display()
                );
            }
        }

        // Layer the Ensembl GRCh37 cdot as a secondary build, if provided.
        if let Some(grch37_str) = manifest
            .get("ensembl_cdot_grch37_json")
            .and_then(|v| v.as_str())
        {
            let grch37_path = resolve_path(grch37_str);
            if grch37_path.exists() {
                match mapper.load_secondary_build(&grch37_path, "GRCh37") {
                    Ok(count) => {
                        eprintln!(
                            "Loaded {} Ensembl GRCh37 transcripts (secondary build)",
                            count
                        );
                    }
                    Err(e) => {
                        eprintln!("Warning: Failed to load Ensembl GRCh37 cdot: {}", e);
                    }
                }
            }
        }
    }

    /// Create a provider from a manifest file.
    ///
    /// This is the standard (runtime) entry point: it loads every manifest key,
    /// including injecting build-time–derived transcript structures
    /// (`derived_transcript_placements`, #790) into the cdot mapper so projection
    /// resolves cdot-absent old versions with real exons.
    pub fn from_manifest<P: AsRef<Path>>(manifest_path: P) -> Result<Self, FerroError> {
        Self::from_manifest_inner(manifest_path, true)
    }

    /// Like [`Self::from_manifest`], but skips injecting the manifest's
    /// `derived_transcript_placements` artifact into the cdot mapper.
    ///
    /// This exists for the #790 build-time producer
    /// (`examples/derive_tx_placements.rs`), which derives those structures by
    /// reasoning over *real* cdot records. If the producer loaded the standard
    /// way once the manifest's `derived_transcript_placements` key is wired to its
    /// own prior output, that output would be injected into cdot and the target
    /// accession (e.g. `NM_003002.2`) would then satisfy `has_transcript_exact`,
    /// causing the producer to skip (decline) it — i.e. in-place regeneration (or
    /// a CI `--check`) would silently produce an empty artifact (#800). Suppressing
    /// the injection makes the producer's cdot view identical to the pre-wired
    /// state, so re-derivation is reliable.
    ///
    /// Note: because the injection is suppressed, a `canonical_overrides` entry
    /// that targets a *derived* old version would not be reconciled here (the
    /// override needs the injected record present to match). That is intentional
    /// and benign for the producer, whose job is to derive structures, not serve
    /// overridden coordinates.
    pub fn from_manifest_without_derived_tx<P: AsRef<Path>>(
        manifest_path: P,
    ) -> Result<Self, FerroError> {
        Self::from_manifest_inner(manifest_path, false)
    }

    /// Shared body of [`Self::from_manifest`] /
    /// [`Self::from_manifest_without_derived_tx`]. When `inject_derived_tx` is
    /// false, the `derived_transcript_placements` injection block is skipped; all
    /// other manifest keys load identically.
    fn from_manifest_inner<P: AsRef<Path>>(
        manifest_path: P,
        inject_derived_tx: bool,
    ) -> Result<Self, FerroError> {
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

        // Get Ensembl cDNA directory (ENST_ accessions)
        if let Some(fastas) = manifest
            .get("ensembl_transcript_fastas")
            .and_then(|v| v.as_array())
        {
            if let Some(first) = fastas.first().and_then(|v| v.as_str()) {
                // Remove .gz extension if present (the indexed sidecar is the
                // decompressed FASTA, as for RefSeq transcripts).
                let path_str = first.strip_suffix(".gz").unwrap_or(first);
                let path = resolve_path(path_str);
                if let Some(ensembl_dir) = path.parent() {
                    if !dirs.contains(&ensembl_dir.to_path_buf()) {
                        dirs.push(ensembl_dir.to_path_buf());
                    }
                }
            }
        }

        // Directory of LRG XML records, used to resolve an LRG parent's
        // chromosomal placement on demand (#480). Derived from the first
        // `lrg_xmls` entry's parent directory.
        let lrg_xml_dir = manifest
            .get("lrg_xmls")
            .and_then(|v| v.as_array())
            .and_then(|xmls| xmls.first())
            .and_then(|v| v.as_str())
            .map(resolve_path)
            .and_then(|p| p.parent().map(Path::to_path_buf));

        // Build the data-driven accession→build table from the assembly
        // report(s) the manifest names (#716). Each report is paired with its
        // build; a read/parse failure warns and is skipped (build inference
        // then uses the hardcoded fallback). An empty table is treated as
        // absent so it never shadows the fallback.
        //
        // This is built BEFORE placement ingestion so the same layered inferer
        // drives both ingest and lookup: a report-only `NC_` accession (one the
        // hardcoded heuristic cannot classify) is then retained and classified
        // during parse/merge instead of being dropped before
        // `select_placement_for_build` ever sees it.
        use crate::reference::transcript::GenomeBuild;
        let report_fields = [
            ("assembly_report", GenomeBuild::GRCh38),
            ("assembly_report_grch37", GenomeBuild::GRCh37),
        ];
        let mut parsed_reports: Vec<(
            GenomeBuild,
            crate::liftover::assembly_report::AssemblyReport,
        )> = Vec::new();
        for (field, build) in report_fields {
            let Some(rel) = manifest.get(field).and_then(|v| v.as_str()) else {
                continue;
            };
            let path = resolve_path(rel);
            match std::fs::read_to_string(&path) {
                Ok(text) => {
                    parsed_reports.push((
                        build,
                        crate::liftover::assembly_report::parse_assembly_report(&text),
                    ));
                }
                Err(e) => warn!(
                    "Failed to read assembly report {}: {}; \
                     build inference uses the hardcoded fallback for this build",
                    path.display(),
                    e
                ),
            }
        }
        let contig_aliases = if parsed_reports.is_empty() {
            None
        } else {
            let pairs: Vec<(
                GenomeBuild,
                &crate::liftover::assembly_report::AssemblyReport,
            )> = parsed_reports.iter().map(|(b, r)| (*b, r)).collect();
            let table = crate::liftover::aliases::ContigAliases::from_assembly_reports(&pairs);
            // Empty parse (0 assembled-molecule rows) → treat as absent. Warn:
            // a manifest named the report(s), so an empty table points at a
            // corrupt/truncated report rather than an intentional absence.
            if table.is_empty() {
                warn!(
                    "Assembly report(s) named by the manifest parsed to 0 assembled-molecule \
                     rows; build inference falls back to the hardcoded version heuristic (the \
                     report file may be corrupt or empty)"
                );
                None
            } else {
                Some(table)
            }
        };

        // The layered inferer used everywhere below: the data-driven table first,
        // the hardcoded heuristic as fallback (#716). Threading it into placement
        // parse/merge keeps ingest consistent with lookup
        // (`MultiFastaProvider::infer_genome_build`).
        let infer_build = |acc: &crate::hgvs::variant::Accession| {
            crate::liftover::aliases::infer_genome_build_layered(contig_aliases.as_ref(), acc)
        };

        // Parse NG_ RefSeqGene→chromosome placements from the NCBI alignment
        // GFF3(s) the manifest names (#480). Two single-build snapshots are
        // merged so an NG_ resolves to its build-appropriate placement: the
        // GRCh38 file (`refseqgene_alignments`) and the GRCh37 file
        // (`refseqgene_alignments_grch37`, #713). Read/parse failures are warned
        // (not silently swallowed): a named-but-unreadable file would otherwise
        // disable NG_ re-anchoring invisibly. An absent field is fine (that build
        // simply has no placements).
        let mut refseqgene_placements: FxHashMap<String, Vec<GenomicPlacement>> =
            FxHashMap::default();
        for field in ["refseqgene_alignments", "refseqgene_alignments_grch37"] {
            let Some(rel) = manifest.get(field).and_then(|v| v.as_str()) else {
                continue;
            };
            let path = resolve_path(rel);
            match std::fs::read_to_string(&path) {
                Ok(gff3) => {
                    let placements = parse_refseqgene_alignments(&gff3, &infer_build);
                    if placements.is_empty() {
                        warn!(
                            "RefSeqGene alignments {} parsed to 0 placements \
                             (this source contributes none)",
                            path.display()
                        );
                    }
                    merge_refseqgene_placements(
                        &mut refseqgene_placements,
                        placements,
                        &infer_build,
                    );
                }
                Err(e) => {
                    warn!(
                        "Failed to read RefSeqGene alignments {}: {}; \
                         NG_ projection disabled for this source",
                        path.display(),
                        e
                    );
                }
            }
        }

        // Derived NG_/LRG_ placements (#728): fill version gaps the authoritative
        // GFF3 snapshots above do not cover (versions absent from NCBI's archived
        // alignments). Merged AFTER the GFF3 sources with first-source-wins, so a
        // GFF3 entry is never overridden — derived placements only add keys/builds
        // the snapshots lack. A read/parse failure is warned, not fatal.
        if let Some(rel) = manifest
            .get("derived_refseqgene_placements")
            .and_then(|v| v.as_str())
        {
            let path = resolve_path(rel);
            match crate::reference::derived_placement::DerivedPlacements::from_json_path(&path) {
                Ok(derived) => {
                    let mut map: FxHashMap<String, Vec<GenomicPlacement>> = FxHashMap::default();
                    for (key, placement) in derived.to_placements() {
                        map.entry(key).or_default().push(placement);
                    }
                    merge_refseqgene_placements(&mut refseqgene_placements, map, &infer_build);
                }
                Err(e) => {
                    warn!(
                        "Failed to read derived NG_ placements {}: {}; \
                         version-gap NG_ projection disabled",
                        path.display(),
                        e
                    );
                }
            }
        }

        // Parse the NCBI `LRG_RefSeqGene` summary into a gene → reference-standard
        // transcript map for legacy gene-model selector resolution (#500/#637).
        // A read failure is warned (not silently swallowed); an absent field is
        // fine (legacy selectors stay un-resolved).
        let legacy_gene_models = match manifest.get("refseqgene_summary").and_then(|v| v.as_str()) {
            None => FxHashMap::default(),
            Some(rel) => {
                let path = resolve_path(rel);
                match std::fs::read_to_string(&path) {
                    Ok(tsv) => crate::reference::legacy_selector::parse_refseqgene_summary(&tsv),
                    Err(e) => {
                        warn!(
                            "Failed to read RefSeqGene summary {}: {}; \
                             legacy gene-model selector resolution disabled",
                            path.display(),
                            e
                        );
                        FxHashMap::default()
                    }
                }
            }
        };

        // Get supplemental directory (missing ClinVar transcripts)
        if let Some(supplemental) = manifest.get("supplemental_fasta").and_then(|v| v.as_str()) {
            let path = resolve_path(supplemental);
            if let Some(supplemental_dir) = path.parent() {
                if !dirs.contains(&supplemental_dir.to_path_buf()) {
                    dirs.push(supplemental_dir.to_path_buf());
                }
            }
        }

        // Version-aware backfill FASTA (#842): deposited sequences for
        // accession.versions present in cdot but absent from the bulk RNA feed.
        // Scanning its directory ensures the backfilled `accession.version` is
        // served by the primary index path (version-exact cdot metadata), so the
        // lossy genome+CIGAR synthesis fallback is never reached for it.
        if let Some(backfill) = manifest
            .get("backfill_transcripts_fasta")
            .and_then(|v| v.as_str())
        {
            let path = resolve_path(backfill);
            let fai = PathBuf::from(format!("{}.fai", path.display()));
            if !path.exists() || !fai.exists() {
                warn!(
                    "Backfill transcript FASTA or its .fai is missing: {}; skipping backfill index",
                    path.display()
                );
            } else if let Some(backfill_dir) = path.parent() {
                if !dirs.contains(&backfill_dir.to_path_buf()) {
                    dirs.push(backfill_dir.to_path_buf());
                }
            }
        }

        if dirs.is_empty() {
            return Err(FerroError::Io {
                msg: "No FASTA directories found in manifest".to_string(),
            });
        }

        let mut provider = Self::from_directories(&dirs)?;
        // Carry the LRG XML directory + NG_ placements so NG_/LRG_ parents can
        // be re-anchored (#480); `from_directories` cannot know about them.
        provider.lrg_xml_dir = lrg_xml_dir;
        provider.refseqgene_placements = refseqgene_placements;
        provider.legacy_gene_models = legacy_gene_models;
        provider.contig_aliases = contig_aliases;

        // Load NG_-hosted transcript map (#792). A read/parse failure is warned
        // (not silently swallowed); an absent field or empty artifact is fine
        // (NG_-parent selector resolution falls back to the global reference map).
        let ng_hosted = match manifest
            .get("ng_hosted_transcripts")
            .and_then(|v| v.as_str())
        {
            None => None,
            Some(rel) => {
                let path = resolve_path(rel);
                match std::fs::read_to_string(&path) {
                    Ok(s) => {
                        match crate::reference::ng_hosted_transcripts::NgHostedTranscripts::from_json(&s) {
                            Ok(nh) if !nh.is_empty() => Some(nh),
                            Ok(_) => None,
                            Err(e) => {
                                warn!(
                                    "Failed to parse ng_hosted_transcripts {}: {}; \
                                     NG_-parent selector resolution falls back to the global map",
                                    path.display(),
                                    e
                                );
                                None
                            }
                        }
                    }
                    Err(e) => {
                        warn!(
                            "Failed to read ng_hosted_transcripts {}: {}; \
                             falling back to the global map",
                            path.display(),
                            e
                        );
                        None
                    }
                }
            }
        };
        provider.ng_hosted = ng_hosted;

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

                        // Defer loading the GRCh37 cdot secondary build if the
                        // manifest provides one. The ~198k-transcript load
                        // happens lazily on the first GRCh37 lookup;
                        // GRCh38/context-free workloads never pay for it.
                        // Best-effort: a missing key/file is a no-op.
                        Self::defer_grch37_secondary_cdot(&mut mapper, &manifest, &resolve_path);

                        // Defer the Ensembl cdot (ENST/ENSG/ENSP) primary merge,
                        // if the manifest provides it (`ferro prepare --ensembl`):
                        // it materializes lazily on first Ensembl use, so a RefSeq
                        // workload never pays for it (#964). Best-effort, like the
                        // LRG/GRCh37 loads above.
                        Self::defer_ensembl_cdot(&mut mapper, &manifest, &resolve_path);

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

        // Derived transcript structures (#790): inject exon→genome maps for
        // old versions cdot lacks, so projection resolves them with real exons.
        // This runs *before* canonical-override reconciliation so a `#520`
        // override keyed to a now-provisioned old version (e.g. `NM_003002.2`)
        // sees the injected record and applies. Injection only gap-fills records
        // absent from cdot (`has_transcript_exact == false`) and reconcile only
        // corrects records present exactly (`has_transcript_exact == true`), so
        // the two are independent except for this ordering.
        //
        // The #790 build-time producer suppresses this injection
        // (`from_manifest_without_derived_tx`) so it derives structures over real
        // cdot records only — see #800. `inject_derived_tx` is true for the
        // runtime `from_manifest` path, so runtime projection is unchanged.
        if inject_derived_tx {
            if let Some(rel) = manifest
                .get("derived_transcript_placements")
                .and_then(|v| v.as_str())
            {
                let path = resolve_path(rel);
                match crate::reference::derived_tx_structure::DerivedTxStructures::from_json_path(
                    &path,
                ) {
                    Ok(derived) => {
                        eprintln!(
                            "Loaded {} derived transcript structures",
                            derived.structures.len()
                        );
                        provider.inject_derived_tx_structures(&derived);
                    }
                    Err(e) => eprintln!(
                        "Warning: Failed to load derived transcript placements {}: {}",
                        path.display(),
                        e
                    ),
                }
            }
        }

        // Canonical overrides (issue #520): load the authoritative-overrides
        // file and reconcile it into the cdot metadata so both the served
        // transcript and the projection path see the corrected CDS/protein.
        // Runs after derived-structure injection so an override targeting a
        // derived old version is applied to the injected record.
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
            has_genomic_data: index_has_genomic_data(&index),
            index,
            aliases: build_chromosome_aliases(),
            cdot_mapper: Some(cdot_mapper),
            supplemental_cds: SupplementalCdsInfo::default(),
            canonical_overrides: CanonicalOverrides::default(),
            protein_index: FxHashMap::default(),
            transcript_cache: new_transcript_cache(),
            open_files: Mutex::new(FxHashMap::default()),
            lrg_xml_dir: None,
            lrg_placement_cache: Mutex::new(FxHashMap::default()),
            refseqgene_placements: FxHashMap::default(),
            legacy_gene_models: FxHashMap::default(),
            contig_aliases: None,
            ng_hosted: None,
        })
    }

    /// Resolve a sequence name, trying various forms
    fn resolve_name(&self, name: &str) -> Option<String> {
        // Direct lookup
        if self.index.contains_key(name) {
            return Some(name.to_string());
        }

        // LRG genomic accession: HGVS `g.` variants reference the genomic
        // sequence by its bare accession `LRG_<N>`, but the FASTA record is
        // indexed under the LRG suffix convention `LRG_<N>g` (transcripts are
        // `LRG_<N>t<k>`, proteins `LRG_<N>p<k>`, both of which match the index
        // directly above). Map the bare genomic accession to its `g`-suffixed
        // record so `LRG_<N>:g.` normalization can fetch reference bases
        // instead of silently no-op'ing (issue #487).
        if is_bare_lrg_genomic_accession(name) {
            let genomic = format!("{name}g");
            if self.index.contains_key(&genomic) {
                return Some(genomic);
            }
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

        // Reuse a cached open handle and read positionally at `file_offset`
        // (one pread, no open/seek/close, no shared cursor) instead of
        // reopening the file on every read.
        let file = {
            let mut handles = self.open_files.lock().map_err(|_| FerroError::Io {
                msg: "open_files mutex poisoned".to_string(),
            })?;
            match handles.get(&entry.file_path) {
                Some(f) => Arc::clone(f),
                None => {
                    let f = Arc::new(File::open(&entry.file_path).map_err(|e| FerroError::Io {
                        msg: format!("Failed to open FASTA file: {}", e),
                    })?);
                    handles.insert(entry.file_path.clone(), Arc::clone(&f));
                    f
                }
            }
        };

        let mut buffer = vec![0u8; bytes_to_read as usize];
        read_exact_at(&file, &mut buffer, file_offset).map_err(|e| FerroError::Io {
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

    /// Whether the FASTA index contains this accession at its **exact**
    /// version, with no version-flex or `chr`-prefix resolution (unlike
    /// [`has_sequence`](Self::has_sequence), which resolves aliases and falls
    /// back to a different version).
    ///
    /// Used by the #629 CDS start-codon check so a cdot CDS coordinate is only
    /// compared against the **matching-version** sequence. A cdot entry whose
    /// exact version is absent from the FASTA is a version-coverage gap
    /// (#645/#653), not a CDS-frame inconsistency, and checking its CDS against
    /// a *different* version's bytes would be a false positive.
    pub fn contains_exact_sequence(&self, name: &str) -> bool {
        self.index.contains_key(name)
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

    /// Return the single transcript `ng` hosts for `gene` from the loaded
    /// `ng_hosted_transcripts` artifact (#792), or `None` when the artifact is
    /// absent, the `(ng, gene)` pair is unknown, or the gene is hosted by more
    /// than one transcript (ambiguous).
    pub(crate) fn ng_hosted_unique(&self, ng: &str, gene: &str) -> Option<String> {
        self.ng_hosted
            .as_ref()?
            .hosted_unique(ng, gene)
            .map(str::to_string)
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

    /// Synthesize the parent-frame (`NG_`/`LRG_`) reference bases for a parent
    /// that has a genomic placement but no FASTA record of its own — e.g. a
    /// #728 *derived* `NG_` whose exact version is absent from the local
    /// RefSeqGene FASTA. The parent sequence is, by construction of the
    /// placement, exactly the placed chromosome (`NC_`) slice (reverse-
    /// complemented for a minus-strand placement), so we serve it by mapping the
    /// requested parent-frame window through the placement and fetching those
    /// `NC_` bases.
    ///
    /// `start`/`end` are a 0-based half-open window in the parent's own frame,
    /// matching [`ReferenceProvider::get_sequence`]. Returns an error when the
    /// window falls outside the placed span.
    fn synthesize_parent_sequence(
        &self,
        placement: &GenomicPlacement,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        if end < start {
            return Err(FerroError::InvalidCoordinates {
                msg: format!("parent sequence window end {end} < start {start}"),
            });
        }
        if start == end {
            return Ok(String::new());
        }
        // The `[start, end)` window is 0-based into the parent *sequence*, so
        // offset 0 is the parent's first base — parent coordinate
        // `parent_start`. The mapping below assumes `parent_start == 1` (it maps
        // the 0-based offset `start` to 1-based parent coordinate `start + 1`).
        // The only reachable caller today is the FASTA-less #728 derived `NG_`,
        // which always has `parent_start == 1` (`derived_placement.rs`); an
        // LRG/RefSeqGene placement with `parent_start > 1` is never routed here.
        // Guard the assumption in debug builds — a `parent_start > 1` placement
        // would otherwise read the wrong window (or fail loud via `parent_to_nc`
        // returning `None`) rather than corrupt bases.
        debug_assert_eq!(
            placement.parent_start, 1,
            "synthesize_parent_sequence assumes parent_start == 1; offset the window \
             by parent_start before mapping to support placement-only parents"
        );
        // 0-based half-open [start, end) → 1-based inclusive parent positions
        // [start+1, end]. Map both endpoints onto the chromosome; on a minus-
        // strand placement the parent runs antiparallel, so the low parent
        // position maps to the *higher* `NC_` coordinate.
        let (nc_a, nc_b) = match (
            placement.parent_to_nc(start + 1),
            placement.parent_to_nc(end),
        ) {
            (Some(a), Some(b)) => (a, b),
            _ => {
                return Err(FerroError::InvalidCoordinates {
                    msg: format!(
                        "parent window [{start}, {end}) falls outside the placed span on {}",
                        placement.nc.full()
                    ),
                });
            }
        };
        let (nc_lo, nc_hi) = if nc_a <= nc_b {
            (nc_a, nc_b)
        } else {
            (nc_b, nc_a)
        };
        // 1-based inclusive [nc_lo, nc_hi] → 0-based half-open for the fetch.
        let bases = self.get_genomic_sequence(&placement.nc.full(), nc_lo - 1, nc_hi)?;
        Ok(match placement.strand {
            crate::reference::Strand::Minus => crate::sequence::reverse_complement(&bases),
            _ => bases,
        })
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
    /// CIGAR-aware (#807). Each exon's per-exon cdot CIGAR is applied while
    /// synthesizing (via [`apply_exon_cigar`], after the exon window is oriented
    /// to the transcript 5'→3' direction): `Deletion` genome bases are skipped so
    /// the served sequence matches the GenBank-deposited transcript, and a CIGAR
    /// `Insertion` — transcript-only bases cdot does not record — makes synthesis
    /// decline (see the `Returns` list) rather than serve a divergent sequence.
    /// Exons with no/empty CIGAR emit their whole genome window unchanged.
    ///
    /// Returns `FerroError::GenomicReferenceNotAvailable` when the cdot-named
    /// contig is missing from the genome FASTA, `FerroError::InvalidCoordinates`
    /// when any exon's genomic span is degenerate, the cdot record carries no
    /// exons, or an exon CIGAR's consumed genome length disagrees with its span,
    /// and `FerroError::TranscriptSequenceUnreconstructable` when an exon CIGAR
    /// carries an insertion. The cdot-named contig is consulted directly (not
    /// version-stripped); callers that rely on a different contig version should
    /// keep the genome FASTA in sync with cdot.
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
        // #807: apply each exon's cdot CIGAR while synthesizing, so the served
        // bases match the GenBank-deposited transcript instead of silently
        // diverging on indels. The CIGAR is walked from the exon's transcript-5'
        // end (see `cigar_deletion_gap_at_genome_pos` / `cigar_insertion_gap_at_tx_pos`):
        //   - `Match(n)`   — keep `n` bases (present in both transcript and genome).
        //   - `Deletion(n)`— skip `n` genome bases (present in genome, absent from
        //     the transcript): applying these is what makes the length/content match.
        //   - `Insertion(n)`— `n` transcript bases with no genome counterpart. cdot
        //     records only the length, never the inserted bases, and `CdotTranscript`
        //     carries no transcript sequence, so they cannot be reconstructed. We
        //     decline (`TranscriptSequenceUnreconstructable`) rather than serve a
        //     sequence provably missing them.
        // To make the CIGAR's tx-5' offset axis coincide with the byte index of the
        // bases we walk, we first orient each exon's genome window to the transcript
        // 5'→3' direction (reverse-complement on the minus strand) and only then
        // apply the CIGAR. Exons with no/empty CIGAR keep the whole window (the
        // pre-#807 behavior — the overwhelming common, non-indel case).
        let is_minus = matches!(tx.strand, Strand::Minus);
        let mut sequence = String::new();
        for (i, e) in tx.exons.iter().enumerate() {
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
            let genome_bases =
                self.get_genomic_sequence(&tx.contig, g_start_hgvs - 1, g_end_excl_hgvs - 1)?;
            // Transcript 5'→3' orientation of this exon's genome window.
            let tx_oriented = if is_minus {
                crate::sequence::reverse_complement(&genome_bases)
            } else {
                genome_bases
            };

            match tx.exon_cigars.get(i) {
                // Indel-bearing or gapless CIGAR present: walk it.
                Some(Some(ops)) if !ops.is_empty() => {
                    // `apply_exon_cigar` only guarantees the CIGAR is consistent
                    // with the *genome* span (Match + Deletion == window length).
                    // It does not tie the emitted *transcript* length back to the
                    // declared exon transcript span `[e[2], e[3])`, which is what
                    // `TxExon.{start,end}` (and downstream CDS coordinates) are
                    // built from below. A record whose Match count disagrees with
                    // `e[3] - e[2]` would yield a sequence whose byte length is
                    // inconsistent with its own exon/CDS metadata, so refuse it
                    // rather than serve mismatched coordinates (#807).
                    let expected_tx_len =
                        e[3].checked_sub(e[2])
                            .ok_or_else(|| FerroError::InvalidCoordinates {
                                msg: format!(
                                    "cdot exon for {id} has invalid transcript span [{}, {})",
                                    e[2], e[3]
                                ),
                            })?;
                    let applied = apply_exon_cigar(id, &tx_oriented, ops)?;
                    if applied.len() as u64 != expected_tx_len {
                        return Err(FerroError::InvalidCoordinates {
                            msg: format!(
                                "cdot CIGAR for {id} emits {} transcript base(s) but the exon \
                                 transcript span is {expected_tx_len}; refusing to synthesize \
                                 from an inconsistent alignment",
                                applied.len()
                            ),
                        });
                    }
                    sequence.push_str(&applied);
                }
                // No per-exon CIGAR data: keep the whole window (pre-#807 behavior).
                _ => sequence.push_str(&tx_oriented),
            }
        }

        // #807/#471: with deletions now applied and insertions declined above, the
        // only residual divergence base synthesis cannot capture is (a) substitution-
        // level differences (cdot's GFF3 Gap CIGAR encodes M/I/D only, never
        // substitutions) and (b) any exon synthesized under the unverified
        // exact-genome-match assumption because it had no usable per-exon CIGAR. The
        // gapless and deletion-applied exons are faithful and add nothing further.
        //
        // The gate is per-exon coverage: a whole-transcript summary of the CIGARs
        // that are *present* would let an applied deletion (or an otherwise-gapless
        // but too-short CIGAR vector) mask a sibling exon that carries no CIGAR at
        // all, wrongly suppressing the warning for a mixed deletion+missing-CIGAR
        // record (#807/#471/#400).
        if cdot_synthesis_has_unverified_exon(tx.exons.len(), &tx.exon_cigars) {
            warn!(
                "{id}: cdot record is missing per-exon alignment-gap (CIGAR) data for one or \
                 more exons — those bases assume an exact genome match and cannot represent \
                 transcript-specific indels or substitutions, so normalize output for {id} may \
                 diverge from the deposited transcript (issues #471, #400)"
            );
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
                // #742: `CdotTranscript.exons` is now HGVS-convention
                // ([genome 1-based incl/excl, tx 0-based incl/excl]); `TxExon`
                // is 1-based-inclusive on both axes. (This also resolves the
                // former #331 genomic-coordinate off-by-one noted here.)
                start: e[2] + 1,             // tx_start: 0-based incl → 1-based incl
                end: e[3],                   // tx_end: 0-based excl = 1-based incl
                genomic_start: Some(e[0]),   // genome_start: already 1-based incl
                genomic_end: Some(e[1] - 1), // genome_end: 1-based excl → 1-based incl
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
            cds_start_incomplete: tx.cds_start_incomplete,
            exon_cigars: tx.exon_cigars.clone(),
            cached_introns: OnceLock::new(),
        })
    }
}

/// Apply one exon's cdot CIGAR to its transcript-5'→3'-oriented genome bases,
/// returning the bases the transcript actually carries for that exon (#807).
///
/// `tx_oriented_bases` must already be in the transcript's 5'→3' direction
/// (reverse-complemented on the minus strand) so the CIGAR's tx-5' offset axis
/// coincides with the byte index — see
/// [`MultiFastaProvider::synthesize_transcript_from_cdot`]. The CIGAR is walked
/// left-to-right from that 5' end:
///   - `Match(n)`    — copy the next `n` bases.
///   - `Deletion(n)` — skip the next `n` genome bases (absent from the transcript).
///   - `Insertion(n)`— transcript-only bases with no genome counterpart; cdot does
///     not record their identity and the genome cannot supply them, so the exon's
///     bases cannot be reconstructed. Returns
///     [`FerroError::TranscriptSequenceUnreconstructable`].
///
/// Returns [`FerroError::InvalidCoordinates`] when the CIGAR's consumed
/// genome length (`Match` + `Deletion`) does not equal the exon's genome span —
/// an internally inconsistent record we refuse rather than mis-slice.
///
/// `ops` is assumed non-empty (the caller only walks present, non-empty CIGARs).
fn apply_exon_cigar(
    id: &str,
    tx_oriented_bases: &str,
    ops: &[crate::data::cdot::CigarOp],
) -> Result<String, FerroError> {
    use crate::data::cdot::CigarOp;

    // Decline the whole synthesis if any exon carries an insertion op. The
    // presence of *any* `Insertion` op — not just a positive total — is the
    // decline trigger: cdot does not record transcript-only base identities, so
    // even a degenerate `Insertion(0)` marks a record we refuse to reconstruct
    // rather than reach the per-op loop's insertion arm. Use checked addition
    // for the reported total: a malformed CIGAR whose insertion counts overflow
    // `u64` is invalid metadata, not a wrapped sum we should trust.
    if ops.iter().any(|op| matches!(op, CigarOp::Insertion(_))) {
        let inserted: u64 = ops
            .iter()
            .try_fold(0u64, |acc, op| {
                let n = match op {
                    CigarOp::Insertion(n) => *n,
                    _ => 0,
                };
                acc.checked_add(n)
            })
            .ok_or_else(|| FerroError::InvalidCoordinates {
                msg: format!("cdot CIGAR for {id} insertion lengths overflow u64"),
            })?;
        return Err(FerroError::TranscriptSequenceUnreconstructable {
            id: id.to_string(),
            insertions: inserted,
        });
    }

    let genome_len = tx_oriented_bases.len() as u64;
    // Checked addition again: a wrapped `Match`/`Deletion` total could spuriously
    // equal `genome_len` and pass span validation while the cursor math below
    // mis-slices, so refuse overflow as invalid metadata.
    let consumed: u64 = ops
        .iter()
        .try_fold(0u64, |acc, op| {
            let n = match op {
                CigarOp::Match(n) | CigarOp::Deletion(n) => *n,
                CigarOp::Insertion(_) => 0,
            };
            acc.checked_add(n)
        })
        .ok_or_else(|| FerroError::InvalidCoordinates {
            msg: format!("cdot CIGAR for {id} genome-consuming lengths overflow u64"),
        })?;
    if consumed != genome_len {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "cdot CIGAR for {id} consumes {consumed} genome base(s) but the exon's genome \
                 span is {genome_len}; refusing to synthesize from an inconsistent alignment"
            ),
        });
    }

    // Genome bases are single-byte ASCII nucleotides, so byte offsets are char
    // boundaries and the `consumed == genome_len` check above guarantees every
    // slice below is in bounds.
    let mut out = String::with_capacity(tx_oriented_bases.len());
    let mut cursor = 0usize;
    for op in ops {
        match op {
            CigarOp::Match(n) => {
                let end = cursor + *n as usize;
                out.push_str(&tx_oriented_bases[cursor..end]);
                cursor = end;
            }
            CigarOp::Deletion(n) => {
                cursor += *n as usize;
            }
            // Any `Insertion` op (including a degenerate `Insertion(0)`) makes
            // the whole synthesis decline above, so the loop never reaches one.
            // Handle it as a no-op anyway — an insertion contributes no genome
            // bases — so a future code path can never panic here.
            CigarOp::Insertion(_) => {}
        }
    }
    Ok(out)
}

/// True when cdot base synthesis emitted at least one of `exon_count` exons under
/// the unverified exact-genome-match assumption — i.e. that exon had no usable
/// per-exon CIGAR: a missing entry (`None`), an empty operation list, or an
/// `exon_cigars` vector shorter than the exon list (so `get(i)` is `None`).
///
/// This is the gate for the best-effort divergence warning in
/// [`MultiFastaProvider::synthesize_transcript_from_cdot`] and mirrors that
/// function's per-exon fallback branch exactly: an exon is synthesized verbatim
/// from the genome window unless its slot is `Some(Some(ops))` with non-empty
/// `ops`. Checking per exon (rather than summarizing the whole CIGAR vector) is
/// what catches a mixed deletion+missing-CIGAR record, where an applied deletion
/// on one exon would otherwise mask a sibling exon that carries no CIGAR at all
/// (#807/#471/#400).
fn cdot_synthesis_has_unverified_exon(
    exon_count: usize,
    exon_cigars: &[Option<Vec<crate::data::cdot::CigarOp>>],
) -> bool {
    (0..exon_count).any(|i| !matches!(exon_cigars.get(i), Some(Some(ops)) if !ops.is_empty()))
}

/// The transcript length [`MultiFastaProvider::synthesize_transcript_from_cdot`]
/// would actually *serve* for `tx`, applying the same per-exon CIGAR rules — or
/// `None` when synthesis would decline rather than emit bases. This mirrors the
/// synthesis path exactly so length-gated callers (e.g.
/// [`MultiFastaProvider::reconcile_cdot_with_overrides`]) compare against the
/// bases the read path produces, not the raw max `tx_end` (#807):
///   - no exons → `None` (synthesis errors).
///   - exon with a present, non-empty CIGAR: a possibly deletion-shortened
///     sequence whose length is the summed `Match` count — but only once the
///     CIGAR is consistent with both the exon's genome span
///     (`Match + Deletion == e[1] - e[0]`) and its transcript span
///     (`Match == e[3] - e[2]`). An insertion CIGAR (unreconstructable) or any
///     inconsistency makes synthesis decline → `None`.
///   - exon with no/empty CIGAR: the whole genome window, `e[1] - e[0]`
///     (pre-#807 behavior).
///
/// Scope mirrors the per-exon CIGAR handling of synthesis, including declining
/// malformed records with degenerate genome coordinates (`e[0] == 0` or
/// `e[1] <= e[0]`), which [`MultiFastaProvider::synthesize_transcript_from_cdot`]
/// rejects before serving any bases. Returning a length for such a record would
/// let override reconciliation compare against bases the read path can never
/// produce, so this helper declines (`None`) on the same degenerate spans.
fn cdot_synthesized_tx_length(tx: &crate::data::cdot::CdotTranscript) -> Option<u64> {
    use crate::data::cdot::CigarOp;
    if tx.exons.is_empty() {
        return None;
    }
    let mut total: u64 = 0;
    for (i, e) in tx.exons.iter().enumerate() {
        // Decline the same degenerate genome spans synthesis refuses, so this
        // length never disagrees with what the read path would serve.
        if e[0] == 0 || e[1] <= e[0] {
            return None;
        }
        let genome_span = e[1] - e[0];
        let exon_len = match tx.exon_cigars.get(i) {
            Some(Some(ops)) if !ops.is_empty() => {
                let tx_span = e[3].checked_sub(e[2])?;
                // Checked accumulation mirrors `apply_exon_cigar`: a malformed
                // CIGAR whose counts overflow `u64` declines (`None`) here too,
                // matching the length the read path would actually serve.
                let (mut matched, mut consumed) = (0u64, 0u64);
                // Mirror `apply_exon_cigar`: the *presence* of any insertion op
                // (including a degenerate `Insertion(0)`) makes the record
                // unreconstructable, so the read path declines and serves no
                // length. Track presence, not just a positive total.
                let mut has_insertion = false;
                for op in ops {
                    match op {
                        CigarOp::Match(n) => {
                            matched = matched.checked_add(*n)?;
                            consumed = consumed.checked_add(*n)?;
                        }
                        CigarOp::Deletion(n) => consumed = consumed.checked_add(*n)?,
                        CigarOp::Insertion(_) => has_insertion = true,
                    }
                }
                // An insertion CIGAR is unreconstructable, and a CIGAR that
                // disagrees with the exon's genome or transcript span is
                // refused — synthesis declines either way, so no length is
                // served.
                if has_insertion || consumed != genome_span || matched != tx_span {
                    return None;
                }
                matched
            }
            // No usable per-exon CIGAR: the whole genome window is served.
            _ => genome_span,
        };
        total = total.checked_add(exon_len)?;
    }
    Some(total)
}

/// Pick the more informative of two transcript-probe errors, preferring whichever
/// is *not* a bare [`FerroError::ReferenceNotFound`]. A typed decline such as
/// [`FerroError::TranscriptSequenceUnreconstructable`] (E2005 — an insertion-CIGAR
/// transcript that cannot be reconstructed) carries a precise "representation
/// unavailable" signal that a generic not-found would mask; once such an error is
/// captured while probing one build, a later `ReferenceNotFound` from probing the
/// other build (or the no-hint fallback) must not overwrite it (#807). Keeps
/// `existing` when it is already informative, otherwise takes `candidate`.
fn prefer_informative_err(existing: Option<FerroError>, candidate: FerroError) -> FerroError {
    match existing {
        Some(e) if !matches!(e, FerroError::ReferenceNotFound { .. }) => e,
        _ => candidate,
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
                // produces via `synthesize_transcript_from_cdot`. That is *not*
                // simply the summed genomic span: since #807 the synthesis path
                // applies each exon's CIGAR, so a deletion-bearing exon serves
                // its `Match` length (after span validation), not its genomic
                // span, and an insertion/inconsistent CIGAR (or empty exons)
                // makes synthesis decline entirely. `cdot_synthesized_tx_length`
                // mirrors those rules, returning `None` when synthesis would not
                // serve bases so reconciliation is skipped. Proceed when the
                // length matches the authoritative length, or when the override
                // carries the canonical `sequence` (the served bases are
                // replaced on read regardless).
                let served_len = self.index.get(acc).map(|e| e.length).or_else(|| {
                    cdot.get_transcript(acc)
                        .and_then(cdot_synthesized_tx_length)
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

    /// Inject build-time–derived exon→genome structures into the cdot mapper so
    /// the unchanged projection path resolves cdot-absent old versions with real
    /// exons (#790). No-op when the cdot mapper is absent. Uses `add_transcript`,
    /// the same lever `reconcile_cdot_with_overrides` uses.
    fn inject_derived_tx_structures(
        &mut self,
        derived: &crate::reference::derived_tx_structure::DerivedTxStructures,
    ) {
        use crate::data::cdot::CdotTranscript;
        use crate::reference::Strand;
        let Some(cdot) = self.cdot_mapper.as_mut() else {
            return;
        };
        for s in &derived.structures {
            let strand = match s.strand.as_str() {
                "+" => Strand::Plus,
                "-" => Strand::Minus,
                _ => continue, // never mis-place on an unparseable strand
            };
            // Don't overwrite a real cdot record if one exists for this exact
            // version (derived structures fill gaps only).
            if cdot.has_transcript_exact(&s.accession) {
                continue;
            }
            cdot.add_transcript(
                s.accession.clone(),
                CdotTranscript {
                    gene_name: s.gene_name.clone(),
                    contig: s.contig.clone(),
                    strand,
                    exons: s.exons.clone(),
                    cds_start: s.cds_start,
                    cds_end: s.cds_end,
                    exon_cigars: vec![None; s.exons.len()],
                    gene_id: None,
                    protein: s.protein.clone(),
                    cds_start_incomplete: false,
                },
            );
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
                    // Version-EXACT cdot lookup (#714): `resolved` is the
                    // concrete accession.version whose *sequence* we just read
                    // from the index. cdot's fuzzy version-fallback would pair
                    // that sequence with a *sibling* version's CDS/exon frame
                    // (e.g. NM_003002.2's bases with .3's CDS, start_codon 84 vs
                    // 61), producing a spec-invalid `c.` description: HGVS numbers
                    // `c.1` from the named version's start codon and requires the
                    // exact version precisely because versions differ in
                    // length/CDS offset (background/numbering.md, refseq.md). On
                    // an exact miss, fall through to the supplemental-CDS path
                    // below rather than adopt a foreign frame.
                    let cdot_tx_opt = match build_hint {
                        Some(b) => cdot.get_transcript_on_build_exact(&resolved, b),
                        None => cdot.get_transcript_exact(&resolved),
                    };
                    if let Some(tx) = cdot_tx_opt {
                        // Convert cdot in-memory exons to transcript exons.
                        // #742: `CdotTranscript.exons` now uses the HGVS
                        // convention — `[genome_start (1-based incl), genome_end
                        // (1-based excl), tx_start (0-based incl), tx_end (0-based
                        // excl)]`. `TxExon` uses 1-based-inclusive coordinates on
                        // both axes, so convert accordingly.
                        let exons: Vec<TxExon> = tx
                            .exons
                            .iter()
                            .enumerate()
                            .map(|(i, e)| TxExon {
                                number: (i + 1) as u32,
                                start: e[2] + 1, // tx_start: 0-based incl → 1-based incl
                                end: e[3],       // tx_end: 0-based excl = 1-based incl
                                genomic_start: Some(e[0]), // genome_start: already 1-based incl
                                genomic_end: Some(e[1] - 1), // genome_end: 1-based excl → incl
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
                    } else if let Some(tx) = (build_hint.is_none())
                        .then(|| cdot.get_transcript_exact_any_build(&resolved))
                        .flatten()
                    {
                        // #718: the exact version is absent from the active cdot
                        // build but present in another loaded build. CDS in
                        // transcript coordinates is build-independent, so adopt
                        // that version's OWN CDS (the correct `c.` frame) rather
                        // than fall to a no-op — and never a fuzzy sibling
                        // (spec-invalid, #714). Synthetic single exon with NO
                        // genomic coordinates: those are build-specific and must
                        // not be borrowed for projection (so the genomic axis
                        // declines rather than emit wrong coordinates).
                        TranscriptMetadata {
                            gene_symbol: tx.gene_name.clone(),
                            strand: tx.strand,
                            cds_start: tx.cds_start.map(|s| s + 1), // 0-based → 1-based
                            cds_end: tx.cds_end,
                            chromosome: None,
                            exons: vec![TxExon {
                                number: 1,
                                start: 1,
                                end: entry.length,
                                genomic_start: None,
                                genomic_end: None,
                            }],
                            genomic_start: None,
                            genomic_end: None,
                            protein_id: None,
                            exon_cigars: Vec::new(),
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
                    cds_start_incomplete: false,
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
                // These bases are reconstructed from the *genome* FASTA via the
                // cdot exon alignment, not read from the deposited transcript
                // sequence — because that sequence is absent from the prepared
                // reference's transcript FASTA. The cdot GFF3 Gap CIGAR encodes
                // only M/I/D, never substitutions, so even with complete CIGAR
                // coverage an aligned ('M') base can differ from the deposited
                // transcript at a substitution site, undetectably (#807/#471/#400).
                // The only fix is to make the authoritative sequence available:
                // add `id` (at the requested version) to the prepared reference's
                // transcript FASTA and re-prepare. We surface that remediation on
                // every synthesis so the divergence is never silent.
                if !cdot.has_transcript_exact(id) {
                    warn!(
                        "{id} is absent from the prepared reference's transcript FASTA and cdot \
                         has no record at this exact version; reconstructing bases from a \
                         sibling-version cdot exon alignment against the genome FASTA. These \
                         genome-derived bases may diverge from the deposited transcript at \
                         substitutions (and at indels for any exon lacking per-exon CIGAR data). \
                         To serve the authoritative sequence, add {id} at the requested version \
                         to the prepared reference's transcript FASTA and re-run `ferro prepare`."
                    );
                } else {
                    warn!(
                        "{id} is absent from the prepared reference's transcript FASTA; \
                         reconstructing bases from the cdot exon alignment against the genome \
                         FASTA. These genome-derived bases may diverge from the deposited \
                         transcript at substitutions (and at indels for any exon lacking \
                         per-exon CIGAR data). To serve the authoritative sequence, add {id} to \
                         the prepared reference's transcript FASTA and re-run `ferro prepare`."
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
    /// Thin `&self` wrapper around [`ReferenceProvider::infer_genome_build`], so the
    /// cdot build-probe path consults this reference's assembly-report-derived table
    /// (when present) before the hardcoded fallback (#716). See
    /// [`crate::liftover::aliases::infer_genome_build_layered`] for the layering contract.
    fn infer_build_from_parent(
        &self,
        parent: &crate::hgvs::variant::Accession,
    ) -> Option<&'static str> {
        self.infer_genome_build(parent)
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

/// Extract the value of an XML attribute written as `name="value"` from a tag
/// slice. Sufficient for the flat, machine-generated LRG XML; not a general XML
/// parser.
fn xml_attr<'a>(tag: &'a str, name: &str) -> Option<&'a str> {
    let needle = format!("{name}=\"");
    let start = tag.find(&needle)? + needle.len();
    let rest = &tag[start..];
    let end = rest.find('"')?;
    Some(&rest[..end])
}

/// Parse an LRG's chromosomal placement from its XML record (#480).
///
/// Reads the `<mapping type="main_assembly" …>` element (the current assembly,
/// e.g. GRCh38) and its single `<mapping_span …>` to build a
/// [`GenomicPlacement`] mapping LRG coordinates onto the `NC_` chromosome.
/// Returns `None` when there is no main-assembly mapping, its `other_id` is not
/// a parseable accession, or it is not a single contiguous span (a single
/// affine placement cannot represent a gapped mapping — the projector then
/// keeps the chromosome coordinates rather than mis-anchor).
fn parse_lrg_main_assembly_placement(xml: &str) -> Option<GenomicPlacement> {
    let mut search = xml;
    loop {
        let open_rel = search.find("<mapping ")?;
        let after_open = &search[open_rel..];
        let tag_end = after_open.find('>')?;
        let open_tag = &after_open[..tag_end];

        let body_start = open_rel + tag_end + 1;
        let close_rel = search[body_start..].find("</mapping>");
        let body = match close_rel {
            Some(c) => &search[body_start..body_start + c],
            None => &search[body_start..],
        };

        if xml_attr(open_tag, "type") == Some("main_assembly") {
            let nc_str = xml_attr(open_tag, "other_id")?;
            let nc = crate::hgvs::parser::accession::parse_accession(nc_str)
                .map(|(_, a)| a)
                .ok()?;

            // Require exactly one span: a single affine transform.
            let spans: Vec<&str> = body
                .match_indices("<mapping_span ")
                .map(|(i, _)| {
                    let s = &body[i..];
                    let e = s.find('>').unwrap_or(s.len());
                    &s[..e]
                })
                .collect();
            if spans.len() != 1 {
                return None;
            }
            let span = spans[0];
            let parent_start: u64 = xml_attr(span, "lrg_start")?.parse().ok()?;
            let o1: u64 = xml_attr(span, "other_start")?.parse().ok()?;
            let o2: u64 = xml_attr(span, "other_end")?.parse().ok()?;
            let strand = match xml_attr(span, "strand") {
                Some("-1") => crate::reference::transcript::Strand::Minus,
                _ => crate::reference::transcript::Strand::Plus,
            };
            let (nc_start, nc_end) = if o1 <= o2 { (o1, o2) } else { (o2, o1) };
            return Some(GenomicPlacement {
                nc,
                parent_start,
                nc_start,
                nc_end,
                strand,
            });
        }

        // Advance past this non-main mapping. If it lacks a `</mapping>` close
        // tag, skip past its open tag and keep scanning rather than aborting the
        // whole search with `?` — a valid `main_assembly` mapping may still
        // follow. `body_start` is strictly past the current `<mapping …>` open
        // tag, so `search` shrinks every iteration (no infinite loop). Real LRG
        // records always close their tags, so this only matters for
        // malformed/truncated input.
        let advance = match close_rel {
            Some(c) => body_start + c + "</mapping>".len(),
            None => body_start,
        };
        search = &search[advance..];
    }
}

/// Extract a GFF3 attribute value (`key=value`, `;`-separated). The value may
/// contain spaces (e.g. `Target=NG_011717.1 1 55255 -`).
fn gff_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|field| {
        // Match the attribute name exactly: split on the first `=` and compare
        // the whole key, so e.g. `Target=` is not matched by a `key` of `Targe`
        // nor a sibling like `Target_note=`.
        let (name, value) = field.trim().split_once('=')?;
        (name == key).then_some(value)
    })
}

/// Parse NCBI RefSeqGene→genome alignments (`GCF_*_refseqgene_alignments.gff3`)
/// into `NG_ accession.version` → [`GenomicPlacement`] (#480).
///
/// Each `match` line places an `NG_` RefSeqGene on a chromosome, e.g.
/// ```text
/// NC_000020.11  RefSeq  match  63404189  63477640  .  +  .  ID=aln;Target=NG_009004.1 1 73452 -;gap_count=0;…
/// ```
/// col 1 is the chromosome (`NC_`), cols 4/5 the chromosome span, and the
/// **`Target`** attribute carries the NG accession, its own coordinate range,
/// and — in its 4th field — the NG's orientation relative to the chromosome.
/// (The GFF `match` strand in col 7 is always `+`; the real orientation is the
/// Target strand.)
///
/// Only **ungapped** (`gap_count=0`) alignments to a chromosome `infer_build`
/// can classify to a build are kept, stored per build so the input's resolved
/// build selects the right placement at lookup (#653). `infer_build` is the
/// caller's layered inferer (data-driven assembly-report table first, hardcoded
/// `NC_`-version heuristic as fallback, #716), so a report-only accession the
/// heuristic alone cannot place is still retained here rather than dropped before
/// lookup. Alt-loci/patch contigs and gapped alignments — which a single affine
/// span cannot represent — are skipped, so the projector declines rather than
/// emit a wrong `NG_` coordinate. At most one placement is kept per (NG version,
/// build); see [`merge_refseqgene_placements`] for combining multiple alignment
/// files (e.g. the separate GRCh38 and GRCh37 RefSeqGene snapshots, #713).
fn parse_refseqgene_alignments(
    gff3: &str,
    infer_build: &impl Fn(&crate::hgvs::variant::Accession) -> Option<&'static str>,
) -> FxHashMap<String, Vec<GenomicPlacement>> {
    let mut out: FxHashMap<String, Vec<GenomicPlacement>> = FxHashMap::default();
    for line in gff3.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        // Only `match` rows place an NG_ parent on a chromosome (col 3 is the
        // GFF3 feature/record type). Skip any other record so a stray row that
        // happens to carry `Target=NG_*` cannot overwrite a real placement.
        if cols[2] != "match" {
            continue;
        }
        let (Ok(a), Ok(b)) = (cols[3].parse::<u64>(), cols[4].parse::<u64>()) else {
            continue;
        };
        let attrs = cols[8];
        // A single affine placement requires an ungapped alignment.
        if gff_attr(attrs, "gap_count") != Some("0") {
            continue;
        }
        let Some(target) = gff_attr(attrs, "Target") else {
            continue;
        };
        let mut tf = target.split_whitespace();
        let (Some(ng), Some(tstart), Some(_tend), Some(tstrand)) =
            (tf.next(), tf.next(), tf.next(), tf.next())
        else {
            continue;
        };
        if !ng.starts_with("NG_") {
            continue;
        }
        let Ok(parent_start) = tstart.parse::<u64>() else {
            continue;
        };
        let strand = if tstrand == "-" {
            crate::reference::transcript::Strand::Minus
        } else {
            crate::reference::transcript::Strand::Plus
        };
        // Keep placements whose build the layered inferer can classify (the
        // data-driven assembly-report table first, the hardcoded NC_-version
        // heuristic as fallback, #716). Alt-loci / fix-patch contigs the inferer
        // cannot place return `None` and are skipped (a single affine span cannot
        // anchor them). Both builds are retained and stored per-build; selection
        // by the input's resolved build happens at lookup (#653).
        let Ok((_, nc)) = crate::hgvs::parser::accession::parse_accession(cols[0]) else {
            continue;
        };
        let Some(build) = infer_build(&nc) else {
            continue;
        };
        let (nc_start, nc_end) = if a <= b { (a, b) } else { (b, a) };
        let placement = GenomicPlacement {
            nc,
            parent_start,
            nc_start,
            nc_end,
            strand,
        };
        // Keep at most one placement per (NG version, build): first wins if a
        // release lists the same alignment twice.
        let entry = out.entry(ng.to_string()).or_default();
        let has_build = entry.iter().any(|p| infer_build(&p.nc) == Some(build));
        if !has_build {
            entry.push(placement);
        }
    }
    out
}

/// Merge the per-`NG_` placement lists parsed from an additional RefSeqGene
/// alignment file into an accumulator, keeping **at most one placement per
/// (NG version, genome build)** — first source wins (#713).
///
/// `ferro prepare` downloads two RefSeqGene→genome alignment snapshots: the
/// GRCh38 release-109 file and the GRCh37 release-105 file. Each is single-build
/// (GRCh38 rows are on `NC_*.11`, GRCh37 on `NC_*.10`-era accessions), so merging
/// them gives each `NG_` both builds, and [`select_placement_for_build`] picks
/// the one matching the input's resolved build. The build is read back from each
/// placement's `NC_` accession via the same injected `infer_build` (the layered
/// inferer, #716), mirroring `parse_refseqgene_alignments`' per-build
/// de-duplication so a build already present from an earlier file is never
/// overwritten.
fn merge_refseqgene_placements(
    acc: &mut FxHashMap<String, Vec<GenomicPlacement>>,
    addition: FxHashMap<String, Vec<GenomicPlacement>>,
    infer_build: &impl Fn(&crate::hgvs::variant::Accession) -> Option<&'static str>,
) {
    for (ng, placements) in addition {
        let entry = acc.entry(ng).or_default();
        for placement in placements {
            let build = infer_build(&placement.nc);
            let has_build = entry.iter().any(|p| infer_build(&p.nc) == build);
            if !has_build {
                entry.push(placement);
            }
        }
    }
}

impl ReferenceProvider for MultiFastaProvider {
    fn genomic_placement(
        &self,
        parent: &crate::hgvs::variant::Accession,
    ) -> Option<GenomicPlacement> {
        // Build-agnostic entry: prefer GRCh38 (cdot's primary / mutalyzer policy).
        self.genomic_placement_on_build(parent, None)
    }

    fn genomic_placement_on_build(
        &self,
        parent: &crate::hgvs::variant::Accession,
        build: Option<&str>,
    ) -> Option<GenomicPlacement> {
        // NG_ RefSeqGene placement comes from the NCBI RefSeqGene→genome
        // alignments (`GCF_*_refseqgene_alignments.gff3`), parsed once at load
        // and stored per genome build (#653). `select_placement_for_build`
        // picks the entry matching the input's resolved build; with no build
        // hint it prefers GRCh38, and with an explicit build it returns `None`
        // (decline) rather than mis-anchor onto a different build's placement.
        if parent.prefix.as_ref() == "NG" {
            let list = self.refseqgene_placements.get(&parent.full())?;
            return crate::reference::provider::select_placement_for_build(list, build, |a| {
                self.infer_genome_build(a)
            });
        }
        // LRG placement comes from the on-hand LRG XMLs (parsed on demand). LRG
        // carries a single main-assembly placement, so the build hint is not
        // used here — a build mismatch is caught downstream by the re-anchor's
        // endpoint-outside-span guard (#655).
        if !parent.is_lrg() {
            return None;
        }
        let dir = self.lrg_xml_dir.as_ref()?;
        let key = parent.full(); // e.g. "LRG_1" (LRG accessions carry no version)

        if let Some(cached) = self
            .lrg_placement_cache
            .lock()
            .expect("lrg placement cache poisoned")
            .get(&key)
        {
            return cached.clone();
        }

        let placement = std::fs::read_to_string(dir.join(format!("{key}.xml")))
            .ok()
            .and_then(|xml| parse_lrg_main_assembly_placement(&xml));
        self.lrg_placement_cache
            .lock()
            .expect("lrg placement cache poisoned")
            .insert(key, placement.clone());
        placement
    }

    fn resolve_legacy_gene_selector(
        &self,
        selector: &str,
        ng_parent: Option<&crate::hgvs::variant::Accession>,
    ) -> Option<String> {
        let ng = ng_parent.map(|a| a.full());
        crate::reference::legacy_selector::resolve_legacy_selector_with_parent(
            selector,
            ng.as_deref(),
            |ng, g| self.ng_hosted_unique(ng, g),
            |g| self.legacy_gene_models.get(g).cloned(),
        )
    }

    fn sole_hosted_transcript(
        &self,
        ng_parent: &crate::hgvs::variant::Accession,
    ) -> Option<String> {
        self.ng_hosted
            .as_ref()?
            .sole_hosted(&ng_parent.full())
            .map(str::to_string)
    }

    fn infer_genome_build(
        &self,
        accession: &crate::hgvs::variant::Accession,
    ) -> Option<&'static str> {
        crate::liftover::aliases::infer_genome_build_layered(
            self.contig_aliases.as_ref(),
            accession,
        )
    }

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
        self.get_transcript_for_accession(accession)
    }

    fn get_transcript_for_accession(
        &self,
        accession: &crate::hgvs::variant::Accession,
    ) -> Result<Arc<Transcript>, FerroError> {
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
        let probe_order: &[&str] = match self.infer_build_from_parent(parent) {
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

        // Track the most informative probe error across builds. A typed decline
        // (e.g. E2005 `TranscriptSequenceUnreconstructable` from an insertion
        // CIGAR) must survive a later build's bare `ReferenceNotFound` so the
        // decline still surfaces for parented variants instead of degrading to a
        // generic not-found (#807).
        let mut last_err: Option<FerroError> = None;
        for build in probe_order {
            match self.get_transcript_on_build_cached(&tx_id, Some(build)) {
                Ok(tx) if tx.chromosome.is_some() => return Ok(tx),
                Ok(tx) => {
                    last_err = Some(prefer_informative_err(
                        last_err,
                        FerroError::ReferenceNotFound { id: tx.id.clone() },
                    ))
                }
                Err(e) => last_err = Some(prefer_informative_err(last_err, e)),
            }
        }

        // Final fallback: try without a build hint (uses cdot's primary build
        // and any supplemental/FASTA-only data). Preserves the historical
        // behavior for transcripts that only live in the primary build. On
        // failure, keep whichever of the accumulated and final errors is more
        // informative, so a no-hint E2005 is not discarded in favor of an
        // earlier `ReferenceNotFound`.
        match self.get_transcript_on_build_cached(&tx_id, None) {
            Ok(tx) => Ok(tx),
            Err(e) => Err(prefer_informative_err(last_err, e)),
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
        // A parent (`NG_`/`LRG_`) that has a genomic placement but no FASTA
        // record of its own — e.g. a #728 derived `NG_` — is served by
        // synthesizing its bases from the chromosome via the placement. Checked
        // before the version-fuzzy `resolve_name` fallback below so an absent
        // *exact* version is reconstructed rather than silently served a sibling
        // version's bases.
        if !self.index.contains_key(id) {
            if let Ok(("", acc)) = crate::hgvs::parser::accession::parse_accession(id) {
                if let Some(placement) = self.genomic_placement(&acc) {
                    return self.synthesize_parent_sequence(&placement, start, end);
                }
            }
        }

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
        // Precomputed at construction (the index is immutable afterward); see
        // `index_has_genomic_data`. Previously this scanned every index key on
        // each call, which dominated the per-variant normalize hot loop.
        self.has_genomic_data
    }

    fn get_sequence_length(&self, id: &str) -> Result<u64, FerroError> {
        // Mirror `get_sequence`: a placement-only parent (#728 derived `NG_`)
        // has no FASTA record, so its length is the placed span on the
        // chromosome. Checked before the version-fuzzy `sequence_length`.
        if !self.index.contains_key(id) {
            if let Ok(("", acc)) = crate::hgvs::parser::accession::parse_accession(id) {
                if let Some(placement) = self.genomic_placement(&acc) {
                    return Ok(placement.nc_end - placement.nc_start + 1);
                }
            }
        }
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

/// Returns `true` for a bare LRG genomic accession (`LRG_<N>` with `N` a
/// non-empty run of ASCII digits and no trailing `t<k>`/`p<k>` selector). These
/// denote the LRG genomic sequence, whose FASTA record is suffixed `g`
/// (`LRG_<N>g`). Transcript (`LRG_<N>t<k>`) and protein (`LRG_<N>p<k>`)
/// accessions return `false` — they match the index directly under their own
/// names and must not be remapped to the genomic record.
fn is_bare_lrg_genomic_accession(name: &str) -> bool {
    match name.strip_prefix("LRG_") {
        Some(digits) => !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()),
        None => false,
    }
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

    /// The bundled hardcoded build inferer, matching what the placement parsers
    /// used before the layered assembly-report table was threaded through (#716).
    /// Tests pass it so their per-build expectations stay pinned to the heuristic.
    fn heuristic_infer(acc: &crate::hgvs::variant::Accession) -> Option<&'static str> {
        crate::liftover::aliases::infer_genome_build_from_accession(acc)
    }

    // ----------------------------------------------------------------------
    // LRG placement parsing (#480)
    // ----------------------------------------------------------------------

    /// Minimal LRG XML with both assemblies and a (multi-span) transcript
    /// mapping, modeled on real records. The parser must pick the single-span
    /// `main_assembly` (GRCh38) mapping and ignore the others.
    const LRG_XML_SAMPLE: &str = r#"
      <mapping coord_system="GRCh37.p13" other_name="17" other_id="NC_000017.10" other_start="48259457" other_end="48284000" type="other_assembly">
        <mapping_span lrg_start="1" lrg_end="24544" other_start="48259457" other_end="48284000" strand="-1" />
      </mapping>
      <mapping coord_system="GRCh38.p13" other_name="17" other_id="NC_000017.11" other_start="50182096" other_end="50206639" type="main_assembly">
        <mapping_span lrg_start="1" lrg_end="24544" other_start="50182096" other_end="50206639" strand="-1" />
      </mapping>
      <mapping coord_system="NM_000088.3" other_name="NM_000088.3" other_id="NM_000088.3" other_start="1" other_end="5927" type="transcript">
        <mapping_span lrg_start="5001" lrg_end="5229" other_start="1" other_end="229" strand="1" />
        <mapping_span lrg_start="6693" lrg_end="6887" other_start="230" other_end="424" strand="1" />
      </mapping>
    "#;

    #[test]
    fn parses_main_assembly_placement_from_lrg_xml() {
        let p = parse_lrg_main_assembly_placement(LRG_XML_SAMPLE)
            .expect("main_assembly mapping should parse");
        assert_eq!(p.nc.to_string(), "NC_000017.11");
        assert_eq!(p.parent_start, 1);
        assert_eq!(p.nc_start, 50182096);
        assert_eq!(p.nc_end, 50206639);
        assert_eq!(p.strand, crate::reference::transcript::Strand::Minus);
        // Spot-check the affine: the LRG (minus strand) base 1 sits at the high
        // chromosome coordinate.
        assert_eq!(p.nc_to_parent(50206639), Some(1));
        assert_eq!(p.nc_to_parent(50182096), Some(24544));
    }

    #[test]
    fn parses_plus_strand_main_assembly_placement_from_lrg_xml() {
        // The shared `LRG_XML_SAMPLE` only covers a minus-strand placement; the
        // plus-strand branch (`strand="1"`, the `Strand::Plus` default) is
        // otherwise exercised only indirectly via the RefSeqGene GFF3 path.
        let xml = r#"
          <mapping coord_system="GRCh38.p13" other_name="1" other_id="NC_000001.11" other_start="1000" other_end="1099" type="main_assembly">
            <mapping_span lrg_start="1" lrg_end="100" other_start="1000" other_end="1099" strand="1" />
          </mapping>
        "#;
        let p = parse_lrg_main_assembly_placement(xml)
            .expect("plus-strand main_assembly mapping should parse");
        assert_eq!(p.nc.to_string(), "NC_000001.11");
        assert_eq!(p.parent_start, 1);
        assert_eq!(p.nc_start, 1000);
        assert_eq!(p.nc_end, 1099);
        assert_eq!(p.strand, crate::reference::transcript::Strand::Plus);
        // Plus strand: the LRG runs parallel to the chromosome, so `nc_to_parent`
        // ascends — the low chromosome coordinate maps to the low parent base.
        assert_eq!(p.nc_to_parent(1000), Some(1));
        assert_eq!(p.nc_to_parent(1050), Some(51));
        assert_eq!(p.nc_to_parent(1099), Some(100));
    }

    #[test]
    fn no_placement_when_main_assembly_mapping_is_absent() {
        let xml = r#"
          <mapping other_id="NC_000017.10" type="other_assembly">
            <mapping_span lrg_start="1" lrg_end="10" other_start="1" other_end="10" strand="1" />
          </mapping>
        "#;
        assert!(parse_lrg_main_assembly_placement(xml).is_none());
    }

    #[test]
    fn main_assembly_found_after_a_malformed_unclosed_mapping() {
        // A truncated/malformed `other_assembly` mapping with no `</mapping>`
        // close tag precedes the (also unclosed, EOF-terminated) main_assembly
        // mapping. The scan must still find the main_assembly placement rather
        // than aborting at the first mapping that lacks a close tag.
        let xml = r#"
          <mapping coord_system="GRCh37" other_id="NC_000017.10" type="other_assembly">
            <mapping_span lrg_start="1" lrg_end="10" other_start="1" other_end="10" strand="1" />
          <mapping coord_system="GRCh38" other_id="NC_000017.11" other_start="100" other_end="109" type="main_assembly">
            <mapping_span lrg_start="1" lrg_end="10" other_start="100" other_end="109" strand="1" />
        "#;
        let p = parse_lrg_main_assembly_placement(xml)
            .expect("main_assembly mapping should parse despite a malformed earlier mapping");
        assert_eq!(p.nc.to_string(), "NC_000017.11");
        assert_eq!(p.parent_start, 1);
        assert_eq!(p.nc_start, 100);
        assert_eq!(p.nc_end, 109);
        assert_eq!(p.strand, crate::reference::transcript::Strand::Plus);
    }

    #[test]
    fn no_placement_when_main_assembly_mapping_is_gapped() {
        // A multi-span main_assembly mapping cannot be a single affine; decline
        // rather than mis-anchor.
        let xml = r#"
          <mapping other_id="NC_000017.11" type="main_assembly">
            <mapping_span lrg_start="1" lrg_end="10" other_start="100" other_end="109" strand="1" />
            <mapping_span lrg_start="11" lrg_end="20" other_start="200" other_end="209" strand="1" />
          </mapping>
        "#;
        assert!(parse_lrg_main_assembly_placement(xml).is_none());
    }

    /// End-to-end provider path: `genomic_placement` reads the LRG XML from the
    /// configured directory and returns its parsed placement (and caches it).
    #[test]
    fn provider_resolves_lrg_placement_from_xml_dir() {
        let (mut provider, dir) = build_provider_with_test_genome();

        // LRG XML directory with one record.
        let xml_dir = dir.path().join("lrg");
        std::fs::create_dir_all(&xml_dir).unwrap();
        std::fs::write(xml_dir.join("LRG_999.xml"), LRG_XML_SAMPLE).unwrap();
        provider.lrg_xml_dir = Some(xml_dir);

        let lrg = crate::hgvs::variant::Accession::new("LRG", "999", None);
        let placement = provider
            .genomic_placement(&lrg)
            .expect("LRG placement should resolve from the XML directory");
        assert_eq!(placement.nc.to_string(), "NC_000017.11");
        assert_eq!(placement.nc_start, 50182096);
        assert_eq!(
            placement.strand,
            crate::reference::transcript::Strand::Minus
        );

        // An NG_ parent with no loaded RefSeqGene alignments gets no placement.
        let ng = crate::hgvs::variant::Accession::new("NG", "007485", Some(1));
        assert!(provider.genomic_placement(&ng).is_none());
    }

    // ----------------------------------------------------------------------
    // RefSeqGene→genome alignment parsing (#480)
    // ----------------------------------------------------------------------

    /// Realistic `GCF_*_refseqgene_alignments.gff3` rows: plus- and minus-strand
    /// NG_ on a GRCh38 primary chromosome, a GRCh37 placement, and a gapped one.
    /// (Tabs are significant.)
    const RSG_GFF3_SAMPLE: &str = "\
##gff-version 3
NC_000001.11\tRefSeq\tmatch\t1000\t1099\t100\t+\t.\tID=aln0;Target=NG_001000.1 1 100 +;gap_count=0\n\
NC_000001.11\tRefSeq\tmatch\t2000\t2099\t100\t+\t.\tID=aln1;Target=NG_002000.1 1 100 -;gap_count=0\n\
NC_000001.10\tRefSeq\tmatch\t500\t599\t100\t+\t.\tID=aln2;Target=NG_003000.1 1 100 -;gap_count=0\n\
NC_000001.11\tRefSeq\tmatch\t3000\t3120\t100\t+\t.\tID=aln3;Target=NG_004000.1 1 100 +;gap_count=2\n\
NC_000001.11\tRefSeq\tcDNA_match\t4000\t4099\t100\t+\t.\tID=aln4;Target=NG_005000.1 1 100 +;gap_count=0\n";

    #[test]
    fn parses_refseqgene_alignments_plus_and_minus() {
        let m = parse_refseqgene_alignments(RSG_GFF3_SAMPLE, &heuristic_infer);

        // GRCh37 is now KEPT (#653); only the gapped and non-`match` rows drop.
        assert_eq!(
            m.len(),
            3,
            "GRCh37 + both GRCh38 ungapped placements kept: {m:?}"
        );
        assert!(
            m.contains_key("NG_003000.1"),
            "GRCh37 placement must be kept (#653)"
        );
        assert!(
            !m.contains_key("NG_004000.1"),
            "gapped placement must be dropped"
        );
        assert!(
            !m.contains_key("NG_005000.1"),
            "non-`match` record (cDNA_match) must be dropped even when it carries Target=NG_*"
        );

        // Each NG_ maps to a per-build list; these samples have one build each.
        let plus = &m["NG_001000.1"][0];
        assert_eq!(plus.nc.to_string(), "NC_000001.11");
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&plus.nc),
            Some("GRCh38")
        );
        assert_eq!(plus.parent_start, 1);
        assert_eq!(plus.nc_start, 1000);
        assert_eq!(plus.nc_end, 1099);
        assert_eq!(plus.strand, crate::reference::transcript::Strand::Plus);
        assert_eq!(plus.nc_to_parent(1000), Some(1));
        assert_eq!(plus.nc_to_parent(1099), Some(100));

        let minus = &m["NG_002000.1"][0];
        assert_eq!(minus.strand, crate::reference::transcript::Strand::Minus);
        // Minus: NG base 1 sits at the high chromosome coordinate.
        assert_eq!(minus.nc_to_parent(2099), Some(1));
        assert_eq!(minus.nc_to_parent(2000), Some(100));

        // The GRCh37 placement is on the GRCh37 chromosome accession (NC_*.10).
        let grch37 = &m["NG_003000.1"][0];
        assert_eq!(grch37.nc.to_string(), "NC_000001.10");
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&grch37.nc),
            Some("GRCh37")
        );
    }

    /// Regression for #716: a "report-only" `NC_` accession — one whose version
    /// the hardcoded heuristic cannot classify, but a data-driven assembly-report
    /// table can — must be retained during alignment parsing, not dropped before
    /// `select_placement_for_build` could ever see it. The placement parser honors
    /// the injected layered inferer, so threading the table through ingest (not
    /// just lookup) is what keeps the data-driven path reachable for these.
    #[test]
    fn parse_keeps_report_only_accession_via_injected_inferer() {
        // `NC_000001.99` is an unknown version: the bundled heuristic declines it.
        const REPORT_ONLY_GFF3: &str = "\
##gff-version 3
NC_000001.99\tRefSeq\tmatch\t1000\t1099\t100\t+\t.\tID=aln0;Target=NG_009000.1 1 100 +;gap_count=0\n";

        // Heuristic-only ingest drops the row (build is unclassifiable) — this is
        // exactly the gap the layering closes.
        let heuristic_only = parse_refseqgene_alignments(REPORT_ONLY_GFF3, &heuristic_infer);
        assert!(
            heuristic_only.is_empty(),
            "the hardcoded heuristic cannot classify NC_000001.99, so it is dropped: {heuristic_only:?}"
        );

        // A data-driven inferer that classifies NC_000001.99 keeps the placement.
        let layered = |acc: &crate::hgvs::variant::Accession| {
            if acc.full() == "NC_000001.99" {
                Some("GRCh38")
            } else {
                heuristic_infer(acc)
            }
        };
        let kept = parse_refseqgene_alignments(REPORT_ONLY_GFF3, &layered);
        let list = kept
            .get("NG_009000.1")
            .expect("report-only NC_ placement must be retained under the injected inferer");
        assert_eq!(list.len(), 1);
        assert_eq!(list[0].nc.to_string(), "NC_000001.99");
    }

    /// `ferro prepare` writes two single-build RefSeqGene alignment files; the
    /// provider merges them so each NG_ carries both builds and
    /// `genomic_placement_on_build` selects the right one (#713). The GRCh38
    /// snapshot lists `NC_*.11` placements; the GRCh37 snapshot `NC_*.10`.
    #[test]
    fn merges_grch38_and_grch37_alignment_files_per_ng() {
        const GRCH38_FILE: &str = "\
##gff-version 3
NC_000001.11\tRefSeq\tmatch\t1000\t1099\t100\t+\t.\tID=a0;Target=NG_007000.1 1 100 +;gap_count=0\n";
        const GRCH37_FILE: &str = "\
##gff-version 3
NC_000001.10\tRefSeq\tmatch\t5000\t5099\t100\t+\t.\tID=a1;Target=NG_007000.1 1 100 +;gap_count=0\n";

        // Merge in prepare's order: GRCh38 first, then GRCh37.
        let mut merged = parse_refseqgene_alignments(GRCH38_FILE, &heuristic_infer);
        merge_refseqgene_placements(
            &mut merged,
            parse_refseqgene_alignments(GRCH37_FILE, &heuristic_infer),
            &heuristic_infer,
        );

        let list = merged.get("NG_007000.1").expect("NG_ present after merge");
        assert_eq!(list.len(), 2, "both builds retained: {list:?}");

        let infer = crate::liftover::aliases::infer_genome_build_from_accession;
        let g38 =
            crate::reference::provider::select_placement_for_build(list, Some("GRCh38"), infer)
                .expect("GRCh38 placement selectable");
        let g37 =
            crate::reference::provider::select_placement_for_build(list, Some("GRCh37"), infer)
                .expect("GRCh37 placement selectable");
        assert_eq!(g38.nc.to_string(), "NC_000001.11");
        assert_eq!(g38.nc_start, 1000);
        assert_eq!(g37.nc.to_string(), "NC_000001.10");
        assert_eq!(g37.nc_start, 5000);

        // Regression: adding the GRCh37 placement must NOT shift the
        // build-agnostic default away from GRCh38 (cdot primary / mutalyzer).
        let default = crate::reference::provider::select_placement_for_build(list, None, infer)
            .expect("default placement selectable");
        assert_eq!(default.nc.to_string(), "NC_000001.11");
    }

    /// `merge_refseqgene_placements` keeps the first source's placement for a
    /// build that is already present (first-source-wins), mirroring the
    /// per-build de-duplication inside `parse_refseqgene_alignments`.
    #[test]
    fn merge_keeps_first_source_per_build() {
        const FIRST: &str = "\
##gff-version 3
NC_000001.11\tRefSeq\tmatch\t1000\t1099\t100\t+\t.\tID=a0;Target=NG_008000.1 1 100 +;gap_count=0\n";
        // A second GRCh38 placement for the same NG_ at a different span.
        const SECOND: &str = "\
##gff-version 3
NC_000001.11\tRefSeq\tmatch\t9000\t9099\t100\t+\t.\tID=a1;Target=NG_008000.1 1 100 +;gap_count=0\n";

        let mut merged = parse_refseqgene_alignments(FIRST, &heuristic_infer);
        merge_refseqgene_placements(
            &mut merged,
            parse_refseqgene_alignments(SECOND, &heuristic_infer),
            &heuristic_infer,
        );

        let list = &merged["NG_008000.1"];
        assert_eq!(list.len(), 1, "only one placement per build is kept");
        assert_eq!(list[0].nc_start, 1000, "first source wins");
    }

    /// An unknown-orientation placement has no defined parent frame, so
    /// `nc_to_parent` declines (`None`) even for an in-span coordinate, rather
    /// than silently treating it as the plus strand.
    #[test]
    fn nc_to_parent_declines_for_unknown_strand() {
        let nc = crate::hgvs::parser::accession::parse_accession("NC_000001.11")
            .expect("valid accession")
            .1;
        let placement = GenomicPlacement {
            nc,
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1099,
            strand: crate::reference::transcript::Strand::Unknown,
        };
        // In-span coordinate, but unknown strand → decline.
        assert_eq!(placement.nc_to_parent(1000), None);
        assert_eq!(placement.nc_to_parent(1050), None);
        // Out-of-span is still None, as before.
        assert_eq!(placement.nc_to_parent(2000), None);
    }

    /// End-to-end provider path: `genomic_placement` returns the parsed NG_
    /// placement for an NG_ parent, and `None` for one with no alignment.
    #[test]
    fn provider_resolves_ng_placement_from_alignments() {
        let (mut provider, _dir) = build_provider_with_test_genome();
        provider.refseqgene_placements =
            parse_refseqgene_alignments(RSG_GFF3_SAMPLE, &heuristic_infer);

        let ng = crate::hgvs::variant::Accession::new("NG", "001000", Some(1));
        let p = provider
            .genomic_placement(&ng)
            .expect("NG_ placement should resolve from the parsed alignments");
        assert_eq!(p.nc.to_string(), "NC_000001.11");
        assert_eq!(p.nc_start, 1000);

        let missing = crate::hgvs::variant::Accession::new("NG", "999999", Some(1));
        assert!(provider.genomic_placement(&missing).is_none());
    }

    // ----------------------------------------------------------------------
    // Hardcoded-heuristic build inference (closes review gap on #332)
    //
    // These exercise the free function
    // `crate::liftover::aliases::infer_genome_build_from_accession` directly —
    // i.e. the version-aware hardcoded alias table, independent of any
    // provider. The `&self` `infer_build_from_parent` method, which routes
    // through the layered data-driven table, is covered separately by
    // `infer_build_from_parent_uses_assembly_report_table` below.
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
            crate::liftover::aliases::infer_genome_build_from_accession(&nc_accession(
                "000017", 11
            )),
            Some("GRCh38")
        );
    }

    #[test]
    fn infer_build_chr17_grch37() {
        // NC_000017.10 is chr17 on GRCh37.
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&nc_accession(
                "000017", 10
            )),
            Some("GRCh37")
        );
    }

    #[test]
    fn infer_build_chr14_uses_alias_table_not_version_heuristic() {
        // chr14 NC accessions are NC_000014.9 (GRCh38) / NC_000014.8 (GRCh37).
        // A naive ".11→GRCh38, .10→GRCh37" rule would mis-classify both.
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&nc_accession("000014", 9)),
            Some("GRCh38"),
        );
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&nc_accession("000014", 8)),
            Some("GRCh37"),
        );
    }

    #[test]
    fn infer_build_malformed_nc_returns_none() {
        // NC_000999.999 is not in the human alias table on any build.
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&nc_accession(
                "000999", 999
            )),
            None,
        );
    }

    #[test]
    fn infer_build_ng_parent_returns_none() {
        // NG accessions are build-agnostic; caller must probe builds.
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&ng_accession("012772", 1)),
            None,
        );
    }

    // ---- #471/#807: cdot base-synthesis unverified-exon warning gate ----
    //
    // `cdot_synthesis_has_unverified_exon` decides whether the best-effort
    // divergence warning fires: it is true exactly when at least one of the
    // transcript's exons would be synthesized verbatim from the genome window
    // because it has no usable per-exon CIGAR (missing/empty entry, or a CIGAR
    // vector shorter than the exon list).

    #[test]
    fn unverified_exon_empty_cigars_is_true() {
        // No CIGAR data at all (the common cdot case) → every exon is synthesized
        // under the unverified exact-genome-match assumption.
        assert!(cdot_synthesis_has_unverified_exon(2, &[]));
    }

    #[test]
    fn unverified_exon_all_none_is_true() {
        assert!(cdot_synthesis_has_unverified_exon(2, &[None, None]));
    }

    #[test]
    fn unverified_exon_empty_inner_vec_is_true() {
        // `Some(vec![])` is a slot that is present but carries no operations: that
        // exon has no usable alignment data, so it is unverified.
        assert!(cdot_synthesis_has_unverified_exon(1, &[Some(vec![])]));
    }

    #[test]
    fn unverified_exon_mixed_real_and_empty_inner_vec_is_true() {
        use crate::data::cdot::CigarOp;
        let cigars = vec![Some(vec![CigarOp::Match(185)]), Some(vec![])];
        assert!(cdot_synthesis_has_unverified_exon(2, &cigars));
    }

    #[test]
    fn unverified_exon_full_gapless_coverage_is_false() {
        use crate::data::cdot::CigarOp;
        // Every exon carries a non-empty CIGAR → all verified → no warning.
        let cigars = vec![
            Some(vec![CigarOp::Match(185)]),
            Some(vec![CigarOp::Match(250)]),
        ];
        assert!(!cdot_synthesis_has_unverified_exon(2, &cigars));
    }

    #[test]
    fn unverified_exon_full_deletion_coverage_is_false() {
        use crate::data::cdot::CigarOp;
        // Deletion CIGARs are applied faithfully, so a fully-covered transcript
        // (even with deletions) has no unverified exon.
        let cigars = vec![
            Some(vec![CigarOp::Match(100), CigarOp::Deletion(5)]),
            Some(vec![CigarOp::Match(250)]),
        ];
        assert!(!cdot_synthesis_has_unverified_exon(2, &cigars));
    }

    #[test]
    fn unverified_exon_partial_coverage_is_true() {
        use crate::data::cdot::CigarOp;
        // One exon carries a CIGAR, one does not → the missing exon is unverified.
        let cigars = vec![Some(vec![CigarOp::Match(185)]), None];
        assert!(cdot_synthesis_has_unverified_exon(2, &cigars));
    }

    #[test]
    fn unverified_exon_mixed_deletion_and_missing_is_true() {
        use crate::data::cdot::CigarOp;
        // #807 regression: an applied deletion on one exon must NOT mask a sibling
        // exon that carries no CIGAR at all — the missing exon is still unverified
        // and the warning must fire.
        let cigars = vec![Some(vec![CigarOp::Match(100), CigarOp::Deletion(5)]), None];
        assert!(cdot_synthesis_has_unverified_exon(2, &cigars));
    }

    #[test]
    fn unverified_exon_cigar_vector_shorter_than_exons_is_true() {
        use crate::data::cdot::CigarOp;
        // #807 regression: a gapless CIGAR covering only the first of two exons
        // leaves the second exon with no CIGAR (`get(1)` is `None`) → unverified,
        // even though the present CIGAR is itself gapless.
        let cigars = vec![Some(vec![CigarOp::Match(185)])];
        assert!(cdot_synthesis_has_unverified_exon(2, &cigars));
    }

    #[test]
    fn unverified_exon_ignores_out_of_range_trailing_slots() {
        use crate::data::cdot::CigarOp;
        // Only the first `exon_count` slots are synthesized; a trailing `None`
        // beyond the exon list is never emitted, so it must not trigger a warning.
        let cigars = vec![Some(vec![CigarOp::Match(185)]), None];
        assert!(!cdot_synthesis_has_unverified_exon(1, &cigars));
    }

    // `cdot_synthesized_tx_length` must report the length
    // `synthesize_transcript_from_cdot` would actually serve, applying the same
    // per-exon CIGAR rules so length-gated callers don't compare against a length
    // the read path never produces (#807).

    /// Minimal cdot record for the length-helper tests (only exons + CIGARs matter).
    fn cdot_tx(
        exons: Vec<[u64; 4]>,
        exon_cigars: Vec<Option<Vec<crate::data::cdot::CigarOp>>>,
    ) -> crate::data::cdot::CdotTranscript {
        use crate::data::cdot::CdotTranscript;
        use crate::reference::Strand;
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: None,
            contig: "NC_TEST.1".to_string(),
            strand: Strand::Plus,
            exons,
            cds_start: None,
            cds_end: None,
            gene_id: None,
            protein: None,
            exon_cigars,
        }
    }

    #[test]
    fn synthesized_tx_length_no_cigars_sums_genome_spans() {
        // No per-exon CIGAR → each exon serves its whole genome window: 5 + 4 = 9.
        let tx = cdot_tx(vec![[1, 6, 0, 5], [10, 14, 5, 9]], Vec::new());
        assert_eq!(cdot_synthesized_tx_length(&tx), Some(9));
    }

    #[test]
    fn synthesized_tx_length_deletion_cigar_uses_match_length() {
        use crate::data::cdot::CigarOp;
        // Genome span 12 (1..13), CIGAR M9 D3: synthesis emits the 9 Match bases,
        // so the served length is 9 (the transcript span), not the 12-base span.
        let tx = cdot_tx(
            vec![[1, 13, 0, 9]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Deletion(3)])],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), Some(9));
    }

    #[test]
    fn synthesized_tx_length_mixed_cigar_and_no_cigar_exons() {
        use crate::data::cdot::CigarOp;
        // Exon 0: span 12, M9 D3 → 9 served. Exon 1: no CIGAR → whole 4-base window.
        let tx = cdot_tx(
            vec![[1, 13, 0, 9], [20, 24, 9, 13]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Deletion(3)]), None],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), Some(13));
    }

    #[test]
    fn synthesized_tx_length_insertion_cigar_declines() {
        use crate::data::cdot::CigarOp;
        // Insertion CIGAR → unreconstructable (E2005) → synthesis declines → None.
        let tx = cdot_tx(
            vec![[1, 10, 0, 11]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Insertion(2)])],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), None);
    }

    #[test]
    fn synthesized_tx_length_zero_length_insertion_declines() {
        use crate::data::cdot::CigarOp;
        // A degenerate `Insertion(0)` carries no count but still marks the
        // record unreconstructable in `apply_exon_cigar`, so this helper must
        // decline (`None`) too — staying consistent with the bases the read
        // path would actually serve (it serves none).
        let tx = cdot_tx(
            vec![[1, 10, 0, 9]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Insertion(0)])],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), None);
    }

    #[test]
    fn synthesized_tx_length_cigar_inconsistent_with_genome_span_declines() {
        use crate::data::cdot::CigarOp;
        // Match+Deletion (10) != genome span (12) → InvalidCoordinates → None.
        let tx = cdot_tx(
            vec![[1, 13, 0, 9]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Deletion(1)])],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), None);
    }

    #[test]
    fn synthesized_tx_length_cigar_inconsistent_with_tx_span_declines() {
        use crate::data::cdot::CigarOp;
        // Match (9) != transcript span (8) even though genome consumption matches
        // → the emitted length would disagree with the exon metadata → None.
        let tx = cdot_tx(
            vec![[1, 13, 0, 8]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Deletion(3)])],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), None);
    }

    #[test]
    fn synthesized_tx_length_empty_or_underflow_exons_decline() {
        // No exons → synthesis errors → None.
        assert_eq!(
            cdot_synthesized_tx_length(&cdot_tx(vec![], Vec::new())),
            None
        );
        // A genome end below the start (underflow) → None rather than a panic.
        assert_eq!(
            cdot_synthesized_tx_length(&cdot_tx(vec![[10, 5, 0, 5]], Vec::new())),
            None
        );
        // A zero genome start (`e[0] == 0`) is a degenerate span synthesis
        // rejects outright, so the helper must decline too — not report
        // `e[1] - e[0]` for bases the read path will never serve.
        assert_eq!(
            cdot_synthesized_tx_length(&cdot_tx(vec![[0, 5, 0, 5]], Vec::new())),
            None
        );
        // An equal start/end (`e[1] == e[0]`) is likewise degenerate → None.
        assert_eq!(
            cdot_synthesized_tx_length(&cdot_tx(vec![[5, 5, 0, 0]], Vec::new())),
            None
        );
    }

    #[test]
    fn synthesized_tx_length_overflowing_cigar_counts_decline() {
        use crate::data::cdot::CigarOp;
        // A `Match` count that overflows `u64` when accumulated must decline
        // rather than wrap, matching `apply_exon_cigar`'s checked arithmetic.
        let tx = cdot_tx(
            vec![[1, 13, 0, 9]],
            vec![Some(vec![CigarOp::Match(u64::MAX), CigarOp::Match(1)])],
        );
        assert_eq!(cdot_synthesized_tx_length(&tx), None);
    }

    // `prefer_informative_err` keeps a typed decline (E2005) from being masked by
    // a bare `ReferenceNotFound` accumulated while probing another build (#807).

    #[test]
    fn prefer_informative_err_keeps_existing_typed_decline() {
        let existing = FerroError::TranscriptSequenceUnreconstructable {
            id: "NM_TEST.1".to_string(),
            insertions: 2,
        };
        let candidate = FerroError::ReferenceNotFound {
            id: "NM_TEST.1".to_string(),
        };
        assert!(matches!(
            prefer_informative_err(Some(existing), candidate),
            FerroError::TranscriptSequenceUnreconstructable { .. }
        ));
    }

    #[test]
    fn prefer_informative_err_replaces_existing_reference_not_found() {
        // A bare not-found gives way to a later typed decline.
        let existing = FerroError::ReferenceNotFound {
            id: "NM_TEST.1".to_string(),
        };
        let candidate = FerroError::TranscriptSequenceUnreconstructable {
            id: "NM_TEST.1".to_string(),
            insertions: 1,
        };
        assert!(matches!(
            prefer_informative_err(Some(existing), candidate),
            FerroError::TranscriptSequenceUnreconstructable { .. }
        ));
    }

    #[test]
    fn prefer_informative_err_none_takes_candidate() {
        let candidate = FerroError::ReferenceNotFound {
            id: "NM_TEST.1".to_string(),
        };
        assert!(matches!(
            prefer_informative_err(None, candidate),
            FerroError::ReferenceNotFound { .. }
        ));
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
    fn test_from_manifest_resolves_ensembl_sequence_and_cdot() {
        use crate::reference::ReferenceProvider;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // Ensembl cDNA FASTA (ENST_) + its .fai. Header ">ENST00000375549.8\n"
        // is 19 bytes (1 + 17 + 1), so the sequence offset is 19.
        fs::write(
            d.join("ensembl_cdna.fna"),
            ">ENST00000375549.8\nACGTACGTACGT\n",
        )
        .unwrap();
        fs::write(
            d.join("ensembl_cdna.fna.fai"),
            "ENST00000375549.8\t12\t19\t12\t13\n",
        )
        .unwrap();

        // Minimal RefSeq cdot (so the cdot-loading block fires) ...
        fs::write(
            d.join("refseq_cdot.json"),
            r#"{"transcripts":{"NM_000088.3":{"gene_name":"COL1A1","contig":"NC_000017.11","strand":"+","exons":[[50184096,50184108,0,12]],"cds_start":0,"cds_end":12}}}"#,
        )
        .unwrap();
        // ... and the Ensembl cdot it merges in.
        fs::write(
            d.join("ensembl_cdot.json"),
            r#"{"transcripts":{"ENST00000375549.8":{"gene_name":"EDN1","contig":"NC_000006.12","strand":"+","exons":[[12290360,12290372,0,12]],"cds_start":0,"cds_end":12}}}"#,
        )
        .unwrap();

        // Manifest references everything by relative path (resolved against its dir).
        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": [],
                "ensembl_transcript_fastas": ["ensembl_cdna.fna"],
                "cdot_json": "refseq_cdot.json",
                "ensembl_cdot_json": "ensembl_cdot.json",
                "transcript_count": 2,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();

        // Sequence resolution for the Ensembl transcript (versioned + length).
        assert_eq!(
            provider.get_sequence("ENST00000375549.8", 0, 4).unwrap(),
            "ACGT"
        );
        assert_eq!(
            provider.get_sequence_length("ENST00000375549.8").unwrap(),
            12
        );
        // Version-less fallback resolves to the highest indexed version.
        assert_eq!(
            provider.get_sequence("ENST00000375549", 0, 4).unwrap(),
            "ACGT"
        );

        // Coordinate mapping: the Ensembl cdot merged into the primary build
        // alongside RefSeq, so both resolve through the same mapper.
        let cdot = provider.cdot_mapper().expect("cdot mapper loaded");
        assert!(
            cdot.get_transcript("ENST00000375549.8").is_some(),
            "Ensembl transcript resolves in cdot"
        );
        assert!(
            cdot.get_transcript("NM_000088.3").is_some(),
            "RefSeq transcript still resolves in cdot"
        );
    }

    #[test]
    fn from_manifest_serves_backfilled_transcript_from_index() {
        // #842: a backfill FASTA named in the manifest (in its own subdir, NOT
        // the bulk transcript dir) must be ingested so the accession.version is
        // served from the primary index path rather than synthesized.
        use crate::reference::ReferenceProvider;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // Bulk transcript FASTA that does NOT contain the superseded version.
        fs::write(d.join("tx.fna"), ">NM_000088.4\nAAAA\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000088.4\t4\t13\t4\t5\n").unwrap();

        // Backfill FASTA in its own dedicated `backfill/` subdir, carrying the
        // superseded version's deposited bases.
        let backfill_dir = d.join("backfill");
        fs::create_dir_all(&backfill_dir).unwrap();
        fs::write(
            backfill_dir.join("backfill_transcripts.fna"),
            ">NM_000088.3\nACGTACGT\n",
        )
        .unwrap();
        fs::write(
            backfill_dir.join("backfill_transcripts.fna.fai"),
            "NM_000088.3\t8\t13\t8\t9\n",
        )
        .unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "backfill_transcripts_fasta": "backfill/backfill_transcripts.fna",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();
        let tx = provider
            .get_transcript("NM_000088.3")
            .expect("backfilled NM_000088.3 is served from the index");
        assert_eq!(
            tx.sequence.as_deref(),
            Some("ACGTACGT"),
            "the deposited backfill bases are served, not synthesized"
        );
    }

    #[test]
    fn from_manifest_skips_missing_backfill_fasta() {
        // #842: a stale `backfill_transcripts_fasta` (named in the manifest but
        // the file/.fai absent on disk) must not brick provider load. Backfill
        // is optional, so the entry is warned and skipped while the rest of the
        // provider loads normally.
        use crate::reference::ReferenceProvider;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // Bulk transcript FASTA the provider must still serve.
        fs::write(d.join("tx.fna"), ">NM_000088.4\nAAAA\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000088.4\t4\t13\t4\t5\n").unwrap();

        // Manifest names a backfill FASTA whose directory exists but whose
        // FASTA/.fai were never written (stale path). Create the dir so the
        // failure mode being guarded — pushing a FASTA-less dir into
        // from_directory — would otherwise trigger.
        fs::create_dir_all(d.join("backfill")).unwrap();
        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "backfill_transcripts_fasta": "backfill/backfill_transcripts.fna",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json"))
            .expect("provider loads despite the stale backfill path");
        assert!(
            provider.get_transcript("NM_000088.4").is_ok(),
            "the bulk transcript still resolves after the backfill entry is skipped"
        );
    }

    #[test]
    fn test_from_manifest_loads_derived_ng_placement() {
        use crate::hgvs::variant::Accession;
        use crate::reference::Strand;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // A FASTA dir is required by from_manifest; a minimal transcript FASTA.
        fs::write(d.join("tx.fna"), ">NM_000001.1\nACGT\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000001.1\t4\t13\t4\t5\n").unwrap();

        // A derived-placements artifact (no GFF3 alignment present → this is the
        // only source of the NG_ placement).
        fs::write(
            d.join("derived.json"),
            r#"{"description":"test","placements":[{"parent":"NG_012337.1","nc":"NC_000011.10","nc_start":112081847,"nc_end":112097794,"strand":"+","anchored_by":"NM_003002.2","mismatch_fraction":0.0}]}"#,
        )
        .unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "derived_refseqgene_placements": "derived.json",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();
        let ng = Accession::new("NG", "012337", Some(1));
        let placement = provider
            .genomic_placement(&ng)
            .expect("derived NG_012337.1 placement is served from the manifest");
        assert_eq!(placement.nc.full(), "NC_000011.10");
        assert_eq!(placement.parent_start, 1);
        assert_eq!(placement.nc_start, 112081847);
        assert_eq!(placement.nc_end, 112097794);
        assert_eq!(placement.strand, Strand::Plus);
    }

    #[test]
    fn get_sequence_synthesizes_derived_ng_bases_from_placement() {
        // #737: a parent (`NG_`) with a derived placement but no FASTA record of
        // its own is served by synthesizing its bases from the chromosome via the
        // placement — plus strand verbatim, minus strand reverse-complemented —
        // and that *exact-version* synthesis takes precedence over the
        // version-fuzzy `resolve_name` fallback to a sibling version's FASTA.
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // Chromosome FASTA. The placed span (1-based [5, 10]) is "AACGGT".
        fs::write(d.join("genome.fna"), ">NC_000099.9\nNNNNAACGGTNNNN\n").unwrap();
        fs::write(d.join("genome.fna.fai"), "NC_000099.9\t14\t13\t14\t15\n").unwrap();

        // A *sibling* version of the plus-strand NG_ that IS present as a FASTA,
        // with deliberately different bases. Version-fuzzy resolution would serve
        // these for a request on NG_999999.1; exact-version synthesis must win.
        fs::write(d.join("sibling.fna"), ">NG_999999.2\nTTTTTT\n").unwrap();
        fs::write(d.join("sibling.fna.fai"), "NG_999999.2\t6\t13\t6\t7\n").unwrap();

        // Derived placements: a plus- and a minus-strand NG_ over the same span.
        fs::write(
            d.join("derived.json"),
            r#"{"description":"test","placements":[
                {"parent":"NG_999999.1","nc":"NC_000099.9","nc_start":5,"nc_end":10,"strand":"+","anchored_by":"","mismatch_fraction":0.0},
                {"parent":"NG_888888.1","nc":"NC_000099.9","nc_start":5,"nc_end":10,"strand":"-","anchored_by":"","mismatch_fraction":0.0}
            ]}"#,
        )
        .unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["genome.fna", "sibling.fna"],
                "derived_refseqgene_placements": "derived.json",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();

        // Length is the placed span (no own FASTA record).
        assert_eq!(provider.get_sequence_length("NG_999999.1").unwrap(), 6);
        assert_eq!(provider.get_sequence_length("NG_888888.1").unwrap(), 6);

        // Plus strand: verbatim chromosome slice.
        assert_eq!(
            provider.get_sequence("NG_999999.1", 0, 6).unwrap(),
            "AACGGT"
        );
        // Sub-window (0-based half-open) maps through the placement.
        assert_eq!(provider.get_sequence("NG_999999.1", 0, 3).unwrap(), "AAC");
        assert_eq!(provider.get_sequence("NG_999999.1", 3, 6).unwrap(), "GGT");

        // Minus strand: reverse complement of the same chromosome slice.
        assert_eq!(
            provider.get_sequence("NG_888888.1", 0, 6).unwrap(),
            "ACCGTT"
        );

        // Exact-version synthesis wins over the sibling NG_999999.2 FASTA — a
        // version-fuzzy result would (wrongly) return "TTTTTT".
        assert_ne!(
            provider.get_sequence("NG_999999.1", 0, 6).unwrap(),
            "TTTTTT"
        );

        // The sibling version that DOES have a FASTA is still served from it.
        assert_eq!(
            provider.get_sequence("NG_999999.2", 0, 6).unwrap(),
            "TTTTTT"
        );
    }

    #[test]
    fn test_from_manifest_layers_ensembl_grch37_secondary_build() {
        // Covers the `ensembl_cdot_grch37_json` → `load_secondary_build(.., "GRCh37")`
        // branch in `load_ensembl_cdot`: a GRCh37 Ensembl cdot named by the
        // manifest must be absorbed as the GRCh37 alt build (distinct GRCh37
        // contig accession), resolvable via `get_transcript_on_build`.
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // Ensembl cDNA FASTA so the provider has at least one FASTA directory.
        fs::write(
            d.join("ensembl_cdna.fna"),
            ">ENST00000375549.8\nACGTACGTACGT\n",
        )
        .unwrap();
        fs::write(
            d.join("ensembl_cdna.fna.fai"),
            "ENST00000375549.8\t12\t19\t12\t13\n",
        )
        .unwrap();

        // Minimal RefSeq cdot so the cdot-loading block fires and creates the
        // primary mapper that `load_ensembl_cdot` merges the Ensembl cdot into.
        fs::write(
            d.join("refseq_cdot.json"),
            r#"{"transcripts":{"NM_000088.3":{"gene_name":"COL1A1","contig":"NC_000017.11","strand":"+","exons":[[50184096,50184108,0,12]],"cds_start":0,"cds_end":12}}}"#,
        )
        .unwrap();
        // Primary (GRCh38) Ensembl cdot: contig is the GRCh38 accession.
        fs::write(
            d.join("ensembl_cdot.json"),
            r#"{"transcripts":{"ENST00000375549.8":{"gene_name":"EDN1","contig":"NC_000006.12","strand":"+","exons":[[12290360,12290372,0,12]],"cds_start":0,"cds_end":12}}}"#,
        )
        .unwrap();
        // Secondary (GRCh37) Ensembl cdot: same accession, GRCh37 contig.
        fs::write(
            d.join("ensembl_cdot_grch37.json"),
            r#"{"transcripts":{"ENST00000375549.8":{"gene_name":"EDN1","contig":"NC_000006.11","strand":"+","exons":[[12290560,12290572,0,12]],"cds_start":0,"cds_end":12}}}"#,
        )
        .unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": [],
                "ensembl_transcript_fastas": ["ensembl_cdna.fna"],
                "cdot_json": "refseq_cdot.json",
                "ensembl_cdot_json": "ensembl_cdot.json",
                "ensembl_cdot_grch37_json": "ensembl_cdot_grch37.json",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();
        let cdot = provider.cdot_mapper().expect("cdot mapper loaded");

        // Primary build resolves to the GRCh38 contig.
        let tx38 = cdot
            .get_transcript_on_build("ENST00000375549.8", "GRCh38")
            .expect("Ensembl transcript on GRCh38 (primary build)");
        assert_eq!(tx38.contig, "NC_000006.12");

        // The GRCh37 secondary build (the branch under test) resolves to the
        // GRCh37 contig — proving the secondary cdot was layered, not dropped.
        let tx37 = cdot
            .get_transcript_on_build("ENST00000375549.8", "GRCh37")
            .expect("Ensembl transcript on GRCh37 (secondary build)");
        assert_eq!(tx37.contig, "NC_000006.11");
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
    fn test_from_directory_skips_macos_appledouble_sidecars() {
        // A reference directory copied through macOS tooling can gain AppleDouble sidecar files
        // (names beginning with `._`). Their `._*.fasta.fai` is a binary xattr blob; parsing it as
        // a FAI index fails with "stream did not contain valid UTF-8". The scanner must skip these
        // and still load the real FASTA rather than aborting the whole reference load.
        let dir = tempdir().unwrap();

        // Real, valid FASTA + index.
        let mut fasta_file = File::create(dir.path().join("real.fna")).unwrap();
        writeln!(fasta_file, ">NM_000001.1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();
        let mut fai_file = File::create(dir.path().join("real.fna.fai")).unwrap();
        writeln!(fai_file, "NM_000001.1\t10\t13\t10\t11").unwrap();

        // AppleDouble sidecars with a `.fasta` extension and non-UTF-8 index content.
        std::fs::write(
            dir.path().join("._real.fasta"),
            [0x00, 0x05, 0x16, 0x07, 0xff],
        )
        .unwrap();
        std::fs::write(
            dir.path().join("._real.fasta.fai"),
            [0xff, 0xfe, 0x00, 0x80],
        )
        .unwrap();

        // Must succeed (not abort on the sidecar) and expose the real sequence.
        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();
        assert!(provider.has_sequence("NM_000001.1"));
        assert_eq!(provider.sequence_length("NM_000001.1"), Some(10));
    }

    #[test]
    fn test_contains_exact_sequence_does_not_version_flex() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fna");
        let fai_path = dir.path().join("test.fna.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NM_000001.3").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NM_000001.3\t10\t13\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        // Exact version present → true.
        assert!(provider.contains_exact_sequence("NM_000001.3"));
        // `has_sequence` version-flexes to .3, but `contains_exact_sequence`
        // must NOT: a different version and the versionless base are exact
        // misses even though the FASTA holds .3.
        assert!(provider.has_sequence("NM_000001.2"));
        assert!(!provider.contains_exact_sequence("NM_000001.2"));
        assert!(provider.has_sequence("NM_000001"));
        assert!(!provider.contains_exact_sequence("NM_000001"));
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

    /// LRG FASTA records follow the LRG suffix convention: the genomic
    /// sequence is `LRG_<N>g`, transcripts are `LRG_<N>t<k>`, proteins
    /// `LRG_<N>p<k>`. HGVS `g.` variants reference the genomic sequence by its
    /// *bare* accession `LRG_<N>` (no `g` suffix), so the provider must map
    /// `LRG_<N>` → the indexed `LRG_<N>g` record. Without this, every
    /// `LRG_<N>:g.` normalization silently no-ops (the fetch errors and
    /// `normalize_genome` echoes the input unchanged) — issue #487's
    /// `LRG_303:g.6932_6933ins…` → `g.6908_6932dup` example among them.
    #[test]
    fn test_lrg_genomic_accession_resolves_to_g_suffixed_record() {
        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("lrg.fasta");
        let fai_path = dir.path().join("lrg.fasta.fai");

        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">LRG_303g").unwrap();
        writeln!(fasta_file, "ACGTACGTAC").unwrap();
        writeln!(fasta_file, ">LRG_303t1").unwrap();
        writeln!(fasta_file, "TTTTGGGGCC").unwrap();
        writeln!(fasta_file, ">LRG_303p1").unwrap();
        writeln!(fasta_file, "MKMKMKMKMK").unwrap();

        // FAI: name, length, offset, linebases, linewidth.
        // ">LRG_303g\n" = 10 bytes → first seq offset 10.
        // record 1: 10 bases + "\n" → next header at 21; ">LRG_303t1\n" = 11
        // bytes → second seq offset 32.
        // record 2: 10 bases + "\n" → next header at 43; ">LRG_303p1\n" = 11
        // bytes → third seq offset 54.
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "LRG_303g\t10\t10\t10\t11").unwrap();
        writeln!(fai_file, "LRG_303t1\t10\t32\t10\t11").unwrap();
        writeln!(fai_file, "LRG_303p1\t10\t54\t10\t11").unwrap();

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        // Bare genomic accession resolves to the `g`-suffixed record.
        let seq = provider
            .get_sequence("LRG_303", 0, 4)
            .expect("LRG_303 (genomic) must resolve to LRG_303g");
        assert_eq!(seq, "ACGT");
        assert_eq!(provider.sequence_length("LRG_303"), Some(10));

        // The transcript accession still resolves directly (not remapped to g).
        let tx = provider
            .get_sequence("LRG_303t1", 0, 4)
            .expect("LRG_303t1 (transcript) must resolve directly");
        assert_eq!(tx, "TTTT");

        // The protein accession likewise resolves directly (not remapped to g).
        let protein = provider
            .get_sequence("LRG_303p1", 0, 4)
            .expect("LRG_303p1 (protein) must resolve directly");
        assert_eq!(protein, "MKMK");
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

    /// Repeated and offset reads through the cached file handle return the
    /// correct, consistent bytes — i.e. positioned `read_exact_at` on a reused
    /// handle behaves like the old open-seek-read per call.
    #[test]
    fn cached_handle_reads_are_correct_and_repeatable() {
        // NC_TEST.1 = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC"
        let (provider, _dir) = build_provider_with_test_genome();
        // First read opens + caches the handle; the second reuses it.
        let a = provider.get_sequence("NC_TEST.1", 0, 8).unwrap();
        let b = provider.get_sequence("NC_TEST.1", 0, 8).unwrap();
        assert_eq!(a, "AAAACCCC");
        assert_eq!(a, b, "repeated read via cached handle must be identical");
        // A different offset on the same cached handle reads the right window.
        assert_eq!(
            provider.get_sequence("NC_TEST.1", 4, 12).unwrap(),
            "CCCCGGGG"
        );
        assert_eq!(provider.get_sequence("NC_TEST.1", 36, 40).unwrap(), "CCCC");
    }

    /// Deleting the backing FASTA after the provider is built triggers a
    /// `FerroError::Io` on the first read (cache miss → `File::open` fails).
    #[test]
    fn cached_handle_open_failure_returns_io_error() {
        let dir = tempdir().unwrap();

        let fasta_path = dir.path().join("genome.fna");
        let fai_path = dir.path().join("genome.fna.fai");

        // NC_TEST.1: 10 bases, sequence starts at byte offset 12 (">NC_TEST.1\n")
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">NC_TEST.1").unwrap();
        writeln!(fasta_file, "AAACCCTTTG").unwrap();
        drop(fasta_file);

        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "NC_TEST.1\t10\t11\t10\t11").unwrap();
        drop(fai_file);

        let provider = MultiFastaProvider::from_directory(dir.path()).unwrap();

        // Delete the FASTA so the cache-miss File::open path fails.
        std::fs::remove_file(&fasta_path).unwrap();

        let err = provider.get_sequence("NC_TEST.1", 0, 5).unwrap_err();
        assert!(
            matches!(err, FerroError::Io { .. }),
            "expected FerroError::Io after backing file is removed, got {err:?}"
        );
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
            cds_start_incomplete: false,
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
            cds_start_incomplete: false,
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
            cds_start_incomplete: false,
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
            cds_start_incomplete: false,
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
                cds_start_incomplete: false,
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
    // #807: cdot CIGAR-aware base synthesis.
    //
    // When the synthesis fallback fires (transcript absent from the FASTA but
    // present in cdot), the per-exon CIGAR must be APPLIED so the served bases
    // match the deposited transcript instead of silently diverging on indels:
    //   - Deletion(n): skip those genome bases (genome-only → absent from tx).
    //   - Insertion(n): tx-only bases cdot does not record → DECLINE.
    //   - gapless / no CIGAR: emit the whole genome window (unchanged).
    //
    // Genome `NC_TEST.1` = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC".
    // Single-exon fixtures use the genome window HGVS 11..=30
    // (0-based [10..30]) = "GGTTTTAAAACCCCGGGGTT" (20 bases), the same window
    // the no-CIGAR fallback tests use, so the only behavioral difference is the
    // CIGAR application.
    // ----------------------------------------------------------------------

    /// Build a CdotMapper seating a single-exon `NM_CIG.1` over genome window
    /// HGVS 11..=30 on `NC_TEST.1`, with the given strand and CIGAR. CIGARs
    /// cannot be injected via `from_transcripts` (it discards them), so we build
    /// the `CdotTranscript` by hand.
    fn seat_single_exon_cigar(
        strand: crate::reference::transcript::Strand,
        cigar: Vec<crate::data::cdot::CigarOp>,
    ) -> CdotMapper {
        use crate::data::cdot::{CdotTranscript, CigarOp};
        // Derive the transcript span from the CIGAR so the fixture is internally
        // consistent: the exon carries one transcript base per Match/Insertion
        // operation (Deletions are genome-only). Synthesis now validates the
        // emitted length against this span (#807), so a hardcoded span would
        // spuriously fail the deletion fixtures (Match-sum < genome span).
        let tx_len: u64 = cigar
            .iter()
            .map(|op| match op {
                CigarOp::Match(n) | CigarOp::Insertion(n) => *n,
                CigarOp::Deletion(_) => 0,
            })
            .sum();
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_CIG.1".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: Some("CIG".to_string()),
                contig: "NC_TEST.1".to_string(),
                // cdot exon layout: [genome_start (1-based incl),
                // genome_end (1-based excl), tx_start (0-based), tx_end (0-based excl)].
                exons: vec![[11, 31, 0, tx_len]],
                strand,
                cds_start: Some(0),
                cds_end: Some(17),
                gene_id: None,
                protein: None,
                exon_cigars: vec![Some(cigar)],
            },
        );
        cdot
    }

    #[test]
    fn test_synthesis_applies_cigar_deletion_plus_strand() {
        use crate::data::cdot::CigarOp;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        // CIGAR M5 D3 M12 over "GGTTTTAAAACCCCGGGGTT": drop the 3 genome bases at
        // offsets [5,8) ("TAA") → "GGTTT" + "AACCCCGGGGTT" = "GGTTTAACCCCGGGGTT".
        provider.cdot_mapper = Some(seat_single_exon_cigar(
            Strand::Plus,
            vec![CigarOp::Match(5), CigarOp::Deletion(3), CigarOp::Match(12)],
        ));

        let tx = provider
            .get_transcript("NM_CIG.1")
            .expect("deletion CIGAR must synthesize, not decline");
        assert_eq!(tx.sequence.as_deref(), Some("GGTTTAACCCCGGGGTT"));
    }

    #[test]
    fn test_synthesis_applies_cigar_deletion_minus_strand() {
        use crate::data::cdot::CigarOp;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        // Minus strand: orient the window to tx-5'→3' first =
        // revcomp("GGTTTTAAAACCCCGGGGTT") = "AACCCCGGGGTTTTAAAACC", THEN apply
        // M5 D3 M12: drop offsets [5,8) ("CGG") → "AACCC" + "GTTTTAAAACC"
        // = "AACCCGGTTTTAAAACC".
        provider.cdot_mapper = Some(seat_single_exon_cigar(
            Strand::Minus,
            vec![CigarOp::Match(5), CigarOp::Deletion(3), CigarOp::Match(12)],
        ));

        let tx = provider
            .get_transcript("NM_CIG.1")
            .expect("deletion CIGAR must synthesize, not decline");
        assert_eq!(tx.sequence.as_deref(), Some("AACCCGGTTTTAAAACC"));
    }

    #[test]
    fn test_synthesis_declines_on_cigar_insertion() {
        use crate::data::cdot::CigarOp;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        // CIGAR M5 I3 M12: 3 transcript-only bases cdot does not record. The
        // genome cannot supply them, so synthesis must DECLINE with a typed
        // error rather than serve a sequence missing those bases.
        provider.cdot_mapper = Some(seat_single_exon_cigar(
            Strand::Plus,
            vec![CigarOp::Match(5), CigarOp::Insertion(3), CigarOp::Match(15)],
        ));

        let err = provider
            .get_transcript("NM_CIG.1")
            .expect_err("insertion CIGAR must decline, not serve a divergent sequence");
        assert!(
            matches!(
                err,
                FerroError::TranscriptSequenceUnreconstructable { ref id, insertions: 3 }
                    if id == "NM_CIG.1"
            ),
            "expected TranscriptSequenceUnreconstructable{{insertions:3}}, got {err:?}"
        );
    }

    #[test]
    fn test_synthesis_rejects_cigar_inconsistent_with_exon_span() {
        use crate::data::cdot::CigarOp;
        use crate::reference::transcript::Strand;

        let (mut provider, _kept) = build_provider_with_test_genome();
        // Exon genome span is 20 bases, but M5 D3 M2 consumes only 10 → the CIGAR
        // is inconsistent with the exon; refuse rather than mis-slice.
        provider.cdot_mapper = Some(seat_single_exon_cigar(
            Strand::Plus,
            vec![CigarOp::Match(5), CigarOp::Deletion(3), CigarOp::Match(2)],
        ));

        let err = provider
            .get_transcript("NM_CIG.1")
            .expect_err("CIGAR not covering the exon genome span must be rejected");
        assert!(
            matches!(err, FerroError::InvalidCoordinates { .. }),
            "expected InvalidCoordinates, got {err:?}"
        );
    }

    #[test]
    fn test_synthesis_rejects_cigar_emitting_wrong_transcript_length() {
        use crate::data::cdot::{CdotTranscript, CigarOp};
        use crate::reference::transcript::Strand;

        // The CIGAR is consistent with the *genome* span (M17 D3 consumes the
        // full 20-base window) but emits only 17 transcript bases, while the
        // declared exon transcript span is [0, 20) = 20. The genome-consumption
        // check passes, so this is caught only by the transcript-length check
        // (#807): serving a 17-base sequence under 20-base exon/CDS coordinates
        // would leave the metadata inconsistent, so it must be refused.
        let (mut provider, _kept) = build_provider_with_test_genome();
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_CIG.1".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: Some("CIG".to_string()),
                contig: "NC_TEST.1".to_string(),
                // Deliberately inconsistent: tx span 20 but the CIGAR yields 17.
                exons: vec![[11, 31, 0, 20]],
                strand: Strand::Plus,
                cds_start: Some(0),
                cds_end: Some(17),
                gene_id: None,
                protein: None,
                exon_cigars: vec![Some(vec![CigarOp::Match(17), CigarOp::Deletion(3)])],
            },
        );
        provider.cdot_mapper = Some(cdot);

        let err = provider.get_transcript("NM_CIG.1").expect_err(
            "a CIGAR emitting fewer transcript bases than the exon span must be rejected",
        );
        assert!(
            matches!(err, FerroError::InvalidCoordinates { .. }),
            "expected InvalidCoordinates, got {err:?}"
        );
    }

    #[test]
    fn test_synthesis_gapless_cigar_matches_no_cigar_output() {
        use crate::data::cdot::CigarOp;
        use crate::reference::transcript::Strand;

        // An explicit gapless (M-only) CIGAR must synthesize byte-for-byte
        // identically to the no-CIGAR fallback (the common, non-indel case is
        // unregressed): both emit the whole genome window "GGTTTTAAAACCCCGGGGTT".
        let (mut provider, _kept) = build_provider_with_test_genome();
        provider.cdot_mapper = Some(seat_single_exon_cigar(
            Strand::Plus,
            vec![CigarOp::Match(20)],
        ));
        let gapless = provider
            .get_transcript("NM_CIG.1")
            .expect("gapless CIGAR must synthesize");
        assert_eq!(gapless.sequence.as_deref(), Some("GGTTTTAAAACCCCGGGGTT"));

        // Cross-check against the no-CIGAR path (existing helper, same window).
        let (mut provider2, _kept2) = build_provider_with_test_genome();
        let tx_no_cigar = build_single_exon_synthetic_tx(Strand::Plus);
        provider2.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx_no_cigar)));
        let no_cigar = provider2
            .get_transcript("NM_TEST.1")
            .expect("no-CIGAR fallback must synthesize");
        assert_eq!(gapless.sequence, no_cigar.sequence);
    }

    #[test]
    fn test_synthesis_mixed_deletion_and_missing_cigar_exon() {
        use crate::data::cdot::{CdotTranscript, CigarOp};
        use crate::reference::transcript::Strand;

        // #807 regression: a two-exon transcript where exon 0 carries a deletion
        // CIGAR (applied) and exon 1 carries no CIGAR (synthesized verbatim). The
        // applied deletion must not mask exon 1's missing data: synthesis still
        // succeeds and emits each exon under its own rule.
        //
        // Genome NC_TEST.1 = "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC".
        //   Exon 0: genome HGVS [11,21) (0-based [10..20]) = "GGTTTTAAAA";
        //           CIGAR M3 D2 M5 drops offsets [3,5) ("TT") → "GGT" + "TAAAA"
        //           = "GGTTAAAA".
        //   Exon 1: genome HGVS [21,31) (0-based [20..30]) = "CCCCGGGGTT";
        //           no CIGAR → whole window "CCCCGGGGTT".
        //   Synthesized = "GGTTAAAA" + "CCCCGGGGTT" = "GGTTAAAACCCCGGGGTT".
        let (mut provider, _kept) = build_provider_with_test_genome();
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_MIX.1".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: Some("MIX".to_string()),
                contig: "NC_TEST.1".to_string(),
                exons: vec![[11, 21, 0, 8], [21, 31, 8, 18]],
                strand: Strand::Plus,
                cds_start: Some(0),
                cds_end: Some(18),
                gene_id: None,
                protein: None,
                exon_cigars: vec![
                    Some(vec![
                        CigarOp::Match(3),
                        CigarOp::Deletion(2),
                        CigarOp::Match(5),
                    ]),
                    None,
                ],
            },
        );
        provider.cdot_mapper = Some(cdot);

        let tx = provider
            .get_transcript("NM_MIX.1")
            .expect("mixed deletion + missing-CIGAR record must synthesize, not decline");
        assert_eq!(tx.sequence.as_deref(), Some("GGTTAAAACCCCGGGGTT"));
    }

    #[test]
    fn test_apply_exon_cigar_pure_helper() {
        use crate::data::cdot::CigarOp;

        // M-only: identity.
        assert_eq!(
            apply_exon_cigar("X", "ACGTACGT", &[CigarOp::Match(8)]).unwrap(),
            "ACGTACGT"
        );
        // Interior deletion: keep first 2 ("AC"), drop next 3 ("GTA"), keep
        // last 3 ("CGT") → "ACCGT".
        assert_eq!(
            apply_exon_cigar(
                "X",
                "ACGTACGT",
                &[CigarOp::Match(2), CigarOp::Deletion(3), CigarOp::Match(3)]
            )
            .unwrap(),
            "ACCGT"
        );
        // Leading deletion.
        assert_eq!(
            apply_exon_cigar("X", "ACGTAC", &[CigarOp::Deletion(2), CigarOp::Match(4)]).unwrap(),
            "GTAC"
        );
        // Trailing deletion.
        assert_eq!(
            apply_exon_cigar("X", "ACGTAC", &[CigarOp::Match(4), CigarOp::Deletion(2)]).unwrap(),
            "ACGT"
        );
        // Any insertion → decline.
        assert!(matches!(
            apply_exon_cigar(
                "X",
                "ACGT",
                &[CigarOp::Match(2), CigarOp::Insertion(1), CigarOp::Match(2)]
            ),
            Err(FerroError::TranscriptSequenceUnreconstructable { insertions: 1, .. })
        ));
        // Degenerate zero-length insertion → decline cleanly, never panic.
        // Pre-fix the `inserted > 0` guard let an `Insertion(0)` slip past and
        // reach the loop's `unreachable!` insertion arm, panicking the
        // synthesizer on adversarial cdot metadata (DoS). Now the *presence* of
        // any insertion op declines with E2005.
        assert!(matches!(
            apply_exon_cigar(
                "X",
                "ACGT",
                &[CigarOp::Match(2), CigarOp::Insertion(0), CigarOp::Match(2)]
            ),
            Err(FerroError::TranscriptSequenceUnreconstructable { insertions: 0, .. })
        ));
        // An `Insertion(0)` as the sole op (no genome-consuming bases) likewise
        // declines rather than panicking.
        assert!(matches!(
            apply_exon_cigar("X", "", &[CigarOp::Insertion(0)]),
            Err(FerroError::TranscriptSequenceUnreconstructable { insertions: 0, .. })
        ));
        // Consumed length != genome span → reject.
        assert!(matches!(
            apply_exon_cigar("X", "ACGTAC", &[CigarOp::Match(2)]),
            Err(FerroError::InvalidCoordinates { .. })
        ));
        // Overflowing insertion total → invalid metadata, not a wrapped sum that
        // could slip past the E2005 branch.
        assert!(matches!(
            apply_exon_cigar(
                "X",
                "ACGT",
                &[CigarOp::Insertion(u64::MAX), CigarOp::Insertion(1)]
            ),
            Err(FerroError::InvalidCoordinates { .. })
        ));
        // Overflowing genome-consuming total → invalid metadata, not a wrapped
        // sum that could spuriously match the genome span.
        assert!(matches!(
            apply_exon_cigar("X", "ACGT", &[CigarOp::Match(u64::MAX), CigarOp::Match(1)]),
            Err(FerroError::InvalidCoordinates { .. })
        ));
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
                cds_start_incomplete: false,
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

    #[test]
    fn has_genomic_data_reflects_presence_of_chromosome_sequences() {
        // A transcript-only index (NM_*) has no genomic (NC_*/chr*) sequences.
        let (tx_only, _d1) = build_provider_with_transcripts(&[
            ("NM_000088.3", "ACGTACGT"),
            ("NR_000001.1", "ACGTACGT"),
        ]);
        assert!(
            !tx_only.has_genomic_data(),
            "transcript-only provider must report no genomic data"
        );

        // An NC_ chromosome accession makes genomic data available.
        let (with_nc, _d2) = build_provider_with_transcripts(&[("NC_000001.11", "ACGTACGTACGT")]);
        assert!(
            with_nc.has_genomic_data(),
            "an NC_ chromosome sequence must report genomic data available"
        );

        // A UCSC-style chr accession likewise counts as genomic.
        let (with_chr, _d3) = build_provider_with_transcripts(&[("chr1", "ACGTACGTACGT")]);
        assert!(
            with_chr.has_genomic_data(),
            "a chr* sequence must report genomic data available"
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
                cds_start_incomplete: false,
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
                cds_start_incomplete: false,
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
                cds_start_incomplete: false,
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

    // --- canonical-sequence reingestion through the provider (issue #791) -----
    //
    // The metadata-only override class (`sequence: None`) is covered above; these
    // exercise the **wrong-sequence** class end-to-end through the provider: an
    // override that carries the canonical bases must cause `get_transcript` to
    // serve those bases (not the wrong served FASTA sequence), and the cdot side
    // must reconcile even when the served length disagrees with the
    // authoritative length. The base-level behavior of `apply_canonical_overrides`
    // is unit-tested in `validate.rs`; this proves the wiring at the provider
    // chokepoint and the no-regression contract for normal / metadata-only paths.

    #[test]
    fn get_transcript_replaces_served_sequence_when_override_carries_canonical_bases() {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        // Served FASTA carries a WRONG, longer sequence (24 nt); the canonical
        // record is 9 nt. The override carries the canonical bases, so the served
        // sequence must be REPLACED on read with the 9 nt canonical sequence
        // (the NM_000193.2 wrong-sequence class, scaled down).
        let served_wrong = "ATGAAATAACCCGGGTTTAAACCCG"; // 25 nt, not canonical
        let canonical = "ATGAAATAA"; // 9 nt canonical
        let (mut provider, _kept) =
            build_provider_with_transcripts(&[("NM_REING.2", served_wrong)]);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_REING.2".to_string(),
            tx_length: canonical.len() as u64,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_REING.2".to_string()),
            sequence: Some(canonical.to_string()),
        });
        provider.canonical_overrides = ov;

        let tx = provider.get_transcript("NM_REING.2").unwrap();
        assert_eq!(
            tx.sequence.as_deref(),
            Some(canonical),
            "the served sequence must be replaced with the canonical bases"
        );
        assert_eq!(
            tx.sequence.as_deref().map(str::len),
            Some(9),
            "served at the canonical 9 nt, not the wrong 25 nt"
        );
        // The CDS/protein metadata are corrected together with the sequence.
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_REING.2"));
        // A length-changing replacement must collapse the exon structure to a
        // single span so the corrected record does not trip the offline
        // `LengthMismatch` check (served-shorter-than-exon-extent).
        let anomalies = crate::reference::validate::validate_transcript_record(&tx);
        assert!(
            anomalies.is_empty(),
            "sequence-corrected transcript must have no self-inflicted anomalies: {anomalies:?}"
        );
    }

    #[test]
    fn metadata_only_override_does_not_alter_served_sequence() {
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        // A `sequence: None` override corrects only CDS/protein; the served bases
        // must be byte-for-byte unchanged (no-regression for the metadata class).
        let served = "ATGAAATAA";
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_META.2", served)]);
        let mut ov = CanonicalOverrides::default();
        ov.insert(AuthoritativeRecord {
            accession: "NM_META.2".to_string(),
            tx_length: 9,
            cds_start: Some(1),
            cds_end: Some(9),
            protein_id: Some("NP_META.2".to_string()),
            sequence: None,
        });
        provider.canonical_overrides = ov;

        let tx = provider.get_transcript("NM_META.2").unwrap();
        assert_eq!(
            tx.sequence.as_deref(),
            Some(served),
            "a metadata-only override must not touch the served sequence"
        );
        assert_eq!((tx.cds_start, tx.cds_end), (Some(1), Some(9)));
        assert_eq!(tx.protein_id.as_deref(), Some("NP_META.2"));
    }

    #[test]
    fn no_override_serves_transcript_unchanged() {
        // Regression guard: with no override loaded, the served sequence is
        // exactly the FASTA bases.
        let served = "ATGCCCGGGTAA";
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_PLAIN.1", served)]);
        let tx = provider.get_transcript("NM_PLAIN.1").unwrap();
        assert_eq!(tx.sequence.as_deref(), Some(served));
    }

    #[test]
    fn reconcile_cdot_proceeds_with_sequence_despite_length_mismatch() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::authoritative::{AuthoritativeRecord, CanonicalOverrides};
        use crate::reference::Strand;
        // Served FASTA is 9 nt but the authoritative length is 99 — a length
        // mismatch that, for a `sequence: None` override, makes reconcile SKIP
        // (see `reconcile_cdot_skips_on_length_mismatch`). When the override
        // carries the canonical bases, reconcile must PROCEED instead (the
        // wrong-sequence class): the served bases will be replaced on read, so
        // the cdot CDS/protein must be reconciled to the canonical values.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("NM_OV.2", "ATGAAATAA")]);
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_OV.2".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
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
            tx_length: 99,      // length mismatch vs served 9 nt
            cds_start: Some(1), // authoritative 1-based
            cds_end: Some(99),  // within tx_length
            protein_id: Some("NP_OV.2".to_string()),
            sequence: Some("A".repeat(99)), // override carries canonical bases
        });
        provider.canonical_overrides = ov;

        provider.reconcile_cdot_with_overrides();

        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_OV.2")
            .unwrap();
        // 1-based start 1 → cdot 0-based 0; end 99 unchanged (0-based exclusive).
        assert_eq!(
            rec.cds_start,
            Some(0),
            "carrying canonical bases must reconcile cdot despite the length mismatch"
        );
        assert_eq!(rec.cds_end, Some(99));
        assert_eq!(rec.protein.as_deref(), Some("NP_OV.2"));
    }

    // --- get_transcript_for_accession (issue #582) ---------------------------

    #[test]
    fn get_transcript_for_accession_with_ng_context() {
        use crate::hgvs::variant::Accession;
        use crate::reference::ReferenceProvider;
        // FASTA-only provider: no cdot, so chromosome is None for all transcripts.
        // With an NG_* genomic_context the implementation probes GRCh38 then
        // GRCh37 (no cdot → no chromosome hit), then falls back to the
        // no-build-hint path, which succeeds because the transcript is in the
        // FASTA index.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_ACC.1", "ATGAAATAA")]);
        let inner = Accession::new("NM", "ACC", Some(1));
        let context = Accession::new("NG", "012345", Some(1));
        let accession = inner.with_genomic_context(context);
        let tx = provider
            .get_transcript_for_accession(&accession)
            .expect("NG_* context must not prevent transcript lookup");
        assert_eq!(tx.id, "NM_ACC.1");
        assert_eq!(tx.sequence.as_deref(), Some("ATGAAATAA"));
    }

    #[test]
    fn get_transcript_for_accession_with_nc_context() {
        use crate::hgvs::variant::Accession;
        use crate::reference::ReferenceProvider;
        // NC_000023.11 → infer_build_from_parent returns GRCh38, so the probe
        // order is [GRCh38, GRCh37].  In a FASTA-only provider neither probe
        // finds a chromosome-bearing transcript, so we fall back to the plain
        // no-build-hint lookup, which still succeeds via the FASTA index.
        let (provider, _kept) = build_provider_with_transcripts(&[("NM_ACC.2", "ATGTGGTAA")]);
        let inner = Accession::new("NM", "ACC", Some(2));
        let context = Accession::new("NC", "000023", Some(11));
        let accession = inner.with_genomic_context(context);
        let tx = provider
            .get_transcript_for_accession(&accession)
            .expect("NC_* context must not prevent transcript lookup");
        assert_eq!(tx.id, "NM_ACC.2");
        assert_eq!(tx.sequence.as_deref(), Some("ATGTGGTAA"));
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
                cds_start_incomplete: false,
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
                cds_start_incomplete: false,
                gene_name: None,
                contig: "NC_TEST.1".to_string(),
                strand: Strand::Plus,
                // 1-based genome bounds (cdot never emits `e[0] == 0`); span 9.
                exons: vec![[1, 10, 0, 9]],
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
        synth_provider_with_cdot_cigars(acc, exons, Vec::new())
    }

    /// As [`synth_provider_with_cdot`], but with explicit per-exon CIGARs so the
    /// served (CIGAR-applied) length can diverge from the raw genomic span (#807).
    #[cfg(test)]
    fn synth_provider_with_cdot_cigars(
        acc: &str,
        exons: Vec<[u64; 4]>,
        exon_cigars: Vec<Option<Vec<crate::data::cdot::CigarOp>>>,
    ) -> (MultiFastaProvider, tempfile::TempDir) {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::Strand;
        let (mut provider, kept) = build_provider_with_test_genome();
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            acc.to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: None,
                contig: "NC_TEST.1".to_string(),
                strand: Strand::Plus,
                exons,
                cds_start: Some(2), // wrong
                cds_end: Some(6),   // wrong
                gene_id: None,
                protein: Some("NP_WRONG.9".to_string()),
                exon_cigars,
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
        // Two exons with genomic spans 5 (1..6) and 4 (10..14) → summed genomic
        // span 9 == authoritative length → reconcile. (Here the genomic span and
        // raw max tx_end coincide; the next test pins the case where they don't.)
        // 1-based genome starts: cdot never emits `e[0] == 0`, and the helper
        // declines that degenerate span (matching synthesis).
        let (mut provider, _kept) =
            synth_provider_with_cdot("NM_MX.1", vec![[1, 6, 0, 5], [10, 14, 5, 9]]);
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
        // A cdot exon whose genomic span (12 bases, 1..13) exceeds its transcript
        // extent (tx_end 9) — the served bases come from the genomic span, so
        // `synthesize_transcript_from_cdot` produces a 12 nt sequence. The
        // authoritative length 12 therefore matches the *served* length, and
        // reconciliation must fire. Keying off the raw max tx_end (9) instead
        // would wrongly skip, gating on a length the server never produced.
        // 1-based genome start (cdot never emits `e[0] == 0`).
        let (mut provider, _kept) = synth_provider_with_cdot("NM_GSPAN.1", vec![[1, 13, 0, 9]]);
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
        // cdot exon served length 9 (1-based genome span 1..10) but authoritative
        // length 100, no sequence → length mismatch → skip.
        let (mut provider, _kept) = synth_provider_with_cdot("NM_MM.1", vec![[1, 10, 0, 9]]);
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

    #[test]
    fn reconcile_cdot_synthesis_only_uses_cigar_applied_length() {
        use crate::data::cdot::CigarOp;
        // Genome span 12 (1..13) but a deletion CIGAR M9 D3: since #807 the read
        // path applies the CIGAR, so the *served* length is the 9-base Match run,
        // not the 12-base genomic span. An override claiming 9 must therefore
        // reconcile (the old genome-span gate would have computed 12 and skipped).
        let (mut provider, _kept) = synth_provider_with_cdot_cigars(
            "NM_DEL.1",
            vec![[1, 13, 0, 9]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Deletion(3)])],
        );
        provider.canonical_overrides = coding_override("NM_DEL.1", 9);
        provider.reconcile_cdot_with_overrides();
        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_DEL.1")
            .unwrap();
        assert_eq!(
            rec.cds_start,
            Some(0),
            "served (CIGAR-applied) length 9 matches override → reconcile"
        );
        assert_eq!(rec.protein.as_deref(), Some("NP_OV.2"));
    }

    #[test]
    fn reconcile_cdot_synthesis_only_cigar_rejects_raw_genomic_span() {
        use crate::data::cdot::CigarOp;
        // Same record, but the override claims the raw 12-base genomic span. The
        // served length is the CIGAR-applied 9, so this must NOT reconcile — the
        // pre-#807 genome-span gate would have wrongly matched 12.
        let (mut provider, _kept) = synth_provider_with_cdot_cigars(
            "NM_DEL2.1",
            vec![[1, 13, 0, 9]],
            vec![Some(vec![CigarOp::Match(9), CigarOp::Deletion(3)])],
        );
        provider.canonical_overrides = coding_override("NM_DEL2.1", 12);
        provider.reconcile_cdot_with_overrides();
        let rec = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_DEL2.1")
            .unwrap();
        assert_eq!(
            rec.cds_start,
            Some(2),
            "raw genomic-span length 12 no longer matches the served length → skip"
        );
        assert_eq!(rec.protein.as_deref(), Some("NP_WRONG.9"));
    }

    // ----------------------------------------------------------------------
    // from_manifest: the `cdot_grch37_json` key defers a GRCh37 secondary build.
    // Exercised via the factored-out `defer_grch37_secondary_cdot` helper so the
    // test needs only a tiny GRCh37 cdot file + an in-memory manifest, not a
    // full hermetic manifest with genome/transcript FASTAs.
    // ----------------------------------------------------------------------

    fn grch37_cdot_json() -> &'static str {
        // Uses the genome_builds nested format (same as real cdot GRCh37 files)
        // so the transcript is stored under alt_build_transcripts["GRCh37"] when
        // loaded, which is the path the deferred secondary mapper consults.
        r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "+",
                            "exons": [[48263025, 48263098, 0, 73, "M73"]],
                            "cds_start": 48263035,
                            "cds_end": 48263085
                        }
                    }
                }
            }
        }
        "#
    }

    #[test]
    fn test_defer_grch37_secondary_cdot_when_key_present() {
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
            "only GRCh38 before deferring the secondary build"
        );

        // Manifest references the GRCh37 cdot by absolute path.
        let manifest = serde_json::json!({
            "cdot_grch37_json": grch37_path.to_str().unwrap(),
        });
        let resolve_path = |p: &str| -> PathBuf { PathBuf::from(p) };

        MultiFastaProvider::defer_grch37_secondary_cdot(&mut mapper, &manifest, &resolve_path);

        // After deferral, the load has NOT happened yet.
        assert!(
            !mapper.deferred_alt_loaded("GRCh37"),
            "deferred, not yet loaded"
        );

        // Triggering a GRCh37 lookup forces the load on demand.
        assert!(
            mapper
                .get_transcript_on_build("NM_000088.3", "GRCh37")
                .is_some(),
            "NM_000088.3 resolvable on GRCh37 after lazy load"
        );

        // The deferred load has now completed.
        assert!(
            mapper.deferred_alt_loaded("GRCh37"),
            "loaded after first GRCh37 use"
        );
    }

    #[test]
    fn test_defer_grch37_secondary_cdot_noop_when_key_absent() {
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
        MultiFastaProvider::defer_grch37_secondary_cdot(&mut mapper, &manifest, &resolve_path);
        // Nothing deferred — the key was absent.
        assert!(
            !mapper.deferred_alt_loaded("GRCh37"),
            "GRCh37 must not be deferred when the cdot_grch37_json key was absent"
        );
        assert_eq!(
            mapper.available_builds_for("NM_000088.3"),
            vec!["GRCh38".to_string()]
        );
    }

    #[test]
    fn from_manifest_loads_assembly_report_contig_aliases() {
        // Prove that MultiFastaProvider::from_manifest ingests the assembly
        // report and wires it into infer_genome_build via ContigAliases.
        //
        // We use NC_000017.99 — a version the hardcoded table does NOT know
        // (the real table only knows .11 for chr17 GRCh38), so if this test
        // passes it must come from the data-driven path, not the fallback.
        use crate::hgvs::variant::Accession;
        use crate::reference::ReferenceProvider;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // A minimal transcript FASTA (required by from_manifest).
        fs::write(d.join("tx.fna"), ">NM_000001.1\nACGT\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000001.1\t4\t13\t4\t5\n").unwrap();

        // Assembly report with a chr17 row using .99 — unknown to the hardcoded table.
        // Format: 10 tab-separated columns matching NCBI's assembly_report layout.
        let report_text = "# Assembly name:  GRCh38.test\n\
            # Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t\
            GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n\
            17\tassembled-molecule\t17\tChromosome\tCM000679.9\t=\tNC_000017.99\t\
            Primary Assembly\t83257441\tchr17\n";
        fs::write(d.join("assembly_report.txt"), report_text).unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "assembly_report": "assembly_report.txt",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();

        // Sanity: hardcoded heuristic does not know .99.
        let future_chr17 = Accession::new("NC", "000017", Some(99));
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&future_chr17),
            None,
            "hardcoded table must not know NC_000017.99 (test premise)"
        );

        // Data-driven path: the loaded ContigAliases table classifies it.
        assert_eq!(
            provider.infer_genome_build(&future_chr17),
            Some("GRCh38"),
            "from_manifest must load the assembly report and wire it into infer_genome_build"
        );
    }

    #[test]
    fn infer_build_from_parent_uses_assembly_report_table() {
        // Direct unit coverage of the `&self` `infer_build_from_parent` method
        // (the probe-order caller at `get_transcript_for_accession`). It must
        // route through the layered data-driven table, NOT the hardcoded
        // heuristic. We again use NC_000017.99 — a version absent from the
        // hardcoded table — so a `Some("GRCh38")` answer can only come from the
        // assembly-report-derived ContigAliases.
        use crate::hgvs::variant::Accession;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // A minimal transcript FASTA (required by from_manifest).
        fs::write(d.join("tx.fna"), ">NM_000001.1\nACGT\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000001.1\t4\t13\t4\t5\n").unwrap();

        // Assembly report with a chr17 row using .99 — unknown to the hardcoded table.
        let report_text = "# Assembly name:  GRCh38.test\n\
            # Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t\
            GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n\
            17\tassembled-molecule\t17\tChromosome\tCM000679.9\t=\tNC_000017.99\t\
            Primary Assembly\t83257441\tchr17\n";
        fs::write(d.join("assembly_report.txt"), report_text).unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "assembly_report": "assembly_report.txt",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();

        let future_chr17 = Accession::new("NC", "000017", Some(99));

        // Test premise: the hardcoded heuristic does not know .99, so a hit here
        // must come from the data-driven table the method routes through.
        assert_eq!(
            crate::liftover::aliases::infer_genome_build_from_accession(&future_chr17),
            None,
            "hardcoded table must not know NC_000017.99 (test premise)"
        );

        // The `&self` method resolves the accession via the layered table.
        assert_eq!(
            provider.infer_build_from_parent(&future_chr17),
            Some("GRCh38"),
            "infer_build_from_parent must route through the assembly-report-derived table"
        );
    }

    /// End-to-end regression for #716: an `NG_` alignment placed on a
    /// report-only `NC_` accession (one the hardcoded heuristic cannot classify)
    /// must survive ingestion and resolve at lookup. This proves the assembly
    /// report is loaded BEFORE placement ingestion and the layered inferer is
    /// threaded through parse/merge — if the report load were still ordered after
    /// ingestion (or parse used the hardcoded fn), the placement would be dropped
    /// and `genomic_placement` would return `None`.
    #[test]
    fn from_manifest_retains_ng_placement_on_report_only_accession() {
        use crate::hgvs::variant::Accession;
        use crate::reference::ReferenceProvider;
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // A minimal transcript FASTA (required by from_manifest).
        fs::write(d.join("tx.fna"), ">NM_000001.1\nACGT\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000001.1\t4\t13\t4\t5\n").unwrap();

        // Assembly report mapping chr17 to NC_000017.99 — a version the hardcoded
        // heuristic does NOT know, so any retained placement on it must come from
        // the data-driven path threaded through parse.
        let report_text = "# Assembly name:  GRCh38.test\n\
            # Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t\
            GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n\
            17\tassembled-molecule\t17\tChromosome\tCM000679.9\t=\tNC_000017.99\t\
            Primary Assembly\t83257441\tchr17\n";
        fs::write(d.join("assembly_report.txt"), report_text).unwrap();

        // A RefSeqGene alignment placing NG_009000.1 on the report-only NC_000017.99.
        let alignments = "##gff-version 3\n\
            NC_000017.99\tRefSeq\tmatch\t1000\t1099\t100\t+\t.\t\
            ID=aln0;Target=NG_009000.1 1 100 +;gap_count=0\n";
        fs::write(d.join("refseqgene_alignments.gff3"), alignments).unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "assembly_report": "assembly_report.txt",
                "refseqgene_alignments": "refseqgene_alignments.gff3",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        let provider = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();

        // The NG_ placement must resolve, anchored on the report-only accession —
        // proof that ingestion retained it via the layered inferer.
        let ng = Accession::new("NG", "009000", Some(1));
        let placement = provider.genomic_placement(&ng).expect(
            "NG_ placement on a report-only NC_ accession must be retained and resolve (#716)",
        );
        assert_eq!(placement.nc.to_string(), "NC_000017.99");
    }

    #[test]
    fn trait_default_infer_genome_build_matches_hardcoded() {
        use crate::hgvs::variant::Accession;
        use crate::reference::provider::ReferenceProvider;
        let provider = crate::reference::mock::MockProvider::new();
        let grch38_chr17 = Accession::new("NC", "000017", Some(11));
        assert_eq!(
            provider.infer_genome_build(&grch38_chr17),
            crate::liftover::aliases::infer_genome_build_from_accession(&grch38_chr17)
        );
    }

    #[test]
    fn from_manifest_loads_ng_hosted_transcripts() {
        use std::fs;
        let dir = tempdir().unwrap();
        let d = dir.path();
        fs::write(d.join("tx.fna"), ">NM_000001.1\nACGT\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000001.1\t4\t13\t4\t5\n").unwrap();
        fs::write(
            d.join("ngh.json"),
            r#"{"schema_version":1,"records":{"NG_012337.1":{"TIMM8B":["NM_012459.2"]}}}"#,
        )
        .unwrap();
        fs::write(
            d.join("manifest.json"),
            r#"{"prepared_at":"t","transcript_fastas":["tx.fna"],"ng_hosted_transcripts":"ngh.json","transcript_count":1,"available_prefixes":[]}"#,
        )
        .unwrap();
        let p = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();
        assert_eq!(
            p.ng_hosted_unique("NG_012337.1", "TIMM8B"),
            Some("NM_012459.2".to_string())
        );
    }

    // --- inject_derived_tx_structures (issue #790) ---------------------------

    #[test]
    fn injects_derived_tx_structure_into_cdot() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::derived_tx_structure::{DerivedTxStructure, DerivedTxStructures};
        use crate::reference::Strand;
        // Build a provider with an empty cdot mapper; inject a derived structure
        // and confirm it is now resolvable at the exact accession.
        let (mut provider, _kept) = build_provider_with_transcripts(&[("DUMMY.1", "A")]);
        provider.cdot_mapper = Some(CdotMapper::new());
        let derived = DerivedTxStructures {
            description: String::new(),
            structures: vec![DerivedTxStructure {
                accession: "NM_TEST.2".to_string(),
                contig: "NC_000011.10".to_string(),
                strand: "+".to_string(),
                exons: vec![[100, 110, 0, 10], [200, 210, 10, 20]],
                cds_start: Some(2),
                cds_end: Some(18),
                gene_name: Some("T".to_string()),
                protein: Some("NP_TEST.1".to_string()),
                anchored_by: "NM_TEST.3".to_string(),
                mismatch_fraction: 0.0,
            }],
        };
        provider.inject_derived_tx_structures(&derived);
        let cdot = provider.cdot_mapper().expect("cdot mapper present");
        assert!(cdot.has_transcript_exact("NM_TEST.2"));
        let t = cdot.get_transcript("NM_TEST.2").unwrap();
        assert_eq!(t.exons, vec![[100, 110, 0, 10], [200, 210, 10, 20]]);
        assert!(matches!(t.strand, Strand::Plus));
        assert_eq!(t.cds_start, Some(2));
        assert_eq!(t.cds_end, Some(18));
        assert_eq!(t.protein.as_deref(), Some("NP_TEST.1"));
    }

    #[test]
    fn inject_skips_unparseable_strand() {
        use crate::data::cdot::CdotMapper;
        use crate::reference::derived_tx_structure::{DerivedTxStructure, DerivedTxStructures};
        let (mut provider, _kept) = build_provider_with_transcripts(&[("DUMMY.1", "A")]);
        provider.cdot_mapper = Some(CdotMapper::new());
        let derived = DerivedTxStructures {
            description: String::new(),
            structures: vec![DerivedTxStructure {
                accession: "NM_BAD.1".to_string(),
                contig: "NC_000001.1".to_string(),
                strand: "?".to_string(), // unparseable
                exons: vec![[1, 10, 0, 9]],
                cds_start: None,
                cds_end: None,
                gene_name: None,
                protein: None,
                anchored_by: "NM_BAD.2".to_string(),
                mismatch_fraction: 0.0,
            }],
        };
        provider.inject_derived_tx_structures(&derived);
        let cdot = provider.cdot_mapper().expect("cdot mapper present");
        // Must be absent — bad strand must not be injected.
        assert!(!cdot.has_transcript_exact("NM_BAD.1"));
    }

    #[test]
    fn inject_does_not_overwrite_existing_cdot_record() {
        use crate::data::cdot::{CdotMapper, CdotTranscript};
        use crate::reference::derived_tx_structure::{DerivedTxStructure, DerivedTxStructures};
        use crate::reference::Strand;
        let (mut provider, _kept) = build_provider_with_transcripts(&[("DUMMY.1", "A")]);
        let mut cdot = CdotMapper::new();
        cdot.add_transcript(
            "NM_REAL.1".to_string(),
            CdotTranscript {
                cds_start_incomplete: false,
                gene_name: Some("REAL".to_string()),
                contig: "NC_000001.1".to_string(),
                strand: Strand::Minus,
                exons: vec![[50, 60, 0, 10]],
                cds_start: Some(1),
                cds_end: Some(9),
                exon_cigars: Vec::new(),
                gene_id: None,
                protein: Some("NP_REAL.1".to_string()),
            },
        );
        provider.cdot_mapper = Some(cdot);
        let derived = DerivedTxStructures {
            description: String::new(),
            structures: vec![DerivedTxStructure {
                accession: "NM_REAL.1".to_string(), // already in cdot
                contig: "NC_000011.10".to_string(),
                strand: "+".to_string(),
                exons: vec![[100, 110, 0, 10]],
                cds_start: Some(2),
                cds_end: Some(8),
                gene_name: Some("DERIVED".to_string()),
                protein: Some("NP_DERIVED.1".to_string()),
                anchored_by: "NM_REAL.2".to_string(),
                mismatch_fraction: 0.0,
            }],
        };
        provider.inject_derived_tx_structures(&derived);
        // The real cdot record must not be overwritten.
        let t = provider
            .cdot_mapper()
            .unwrap()
            .get_transcript("NM_REAL.1")
            .unwrap();
        assert!(
            matches!(t.strand, Strand::Minus),
            "real cdot record must not be overwritten by derived structure"
        );
        assert_eq!(t.gene_name.as_deref(), Some("REAL"));
    }

    /// #800: `from_manifest_without_derived_tx` must NOT inject the manifest's
    /// `derived_transcript_placements` artifact, while `from_manifest` still does.
    /// This is the regression guard for the build-time producer: if it loaded the
    /// standard way once the key points at its own prior output, the target
    /// accession would be present in cdot and the producer would skip (decline) it,
    /// silently emitting an empty artifact on a re-run / CI `--check`.
    #[test]
    fn from_manifest_without_derived_tx_suppresses_injection() {
        use std::fs;

        let dir = tempdir().unwrap();
        let d = dir.path();

        // A FASTA dir is required by from_manifest; a minimal transcript FASTA.
        fs::write(d.join("tx.fna"), ">NM_000001.1\nACGT\n").unwrap();
        fs::write(d.join("tx.fna.fai"), "NM_000001.1\t4\t13\t4\t5\n").unwrap();

        // A real cdot record for the *sibling* version `NM_REAL.3` — the producer
        // anchors to real cdot records, which must always be present regardless of
        // injection suppression.
        fs::write(
            d.join("cdot.json"),
            r#"{
                "transcripts": {
                    "NM_REAL.3": {
                        "gene_name": "REAL",
                        "genome_builds": {
                            "GRCh38": {
                                "contig": "NC_000011.10",
                                "strand": "+",
                                "exons": [[100, 110, 1, 0, 10, "M10"]]
                            }
                        },
                        "start_codon": 2,
                        "stop_codon": 8
                    }
                }
            }"#,
        )
        .unwrap();

        // A derived-placements artifact for the cdot-absent *old* version
        // `NM_REAL.2` — i.e. the producer's own prior output.
        fs::write(
            d.join("derived.json"),
            r#"{
                "description": "test",
                "structures": [{
                    "accession": "NM_REAL.2",
                    "contig": "NC_000011.10",
                    "strand": "+",
                    "exons": [[100, 110, 0, 10]],
                    "cds_start": 2,
                    "cds_end": 8,
                    "gene_name": "REAL",
                    "protein": "NP_REAL.1",
                    "anchored_by": "NM_REAL.3",
                    "mismatch_fraction": 0.0
                }]
            }"#,
        )
        .unwrap();

        fs::write(
            d.join("manifest.json"),
            r#"{
                "prepared_at": "test",
                "transcript_fastas": ["tx.fna"],
                "cdot_json": "cdot.json",
                "derived_transcript_placements": "derived.json",
                "transcript_count": 1,
                "available_prefixes": []
            }"#,
        )
        .unwrap();

        // Standard (runtime) load injects the derived old version into cdot.
        let injected = MultiFastaProvider::from_manifest(d.join("manifest.json")).unwrap();
        let injected_cdot = injected.cdot_mapper().expect("cdot mapper loaded");
        assert!(
            injected_cdot.has_transcript_exact("NM_REAL.2"),
            "from_manifest must inject the derived old version (runtime path unchanged)"
        );
        assert!(
            injected_cdot.has_transcript_exact("NM_REAL.3"),
            "real sibling record must be present"
        );

        // Producer load suppresses the injection: the derived old version is
        // absent, but the real sibling record remains so derivation can proceed.
        let suppressed =
            MultiFastaProvider::from_manifest_without_derived_tx(d.join("manifest.json")).unwrap();
        let suppressed_cdot = suppressed.cdot_mapper().expect("cdot mapper loaded");
        assert!(
            !suppressed_cdot.has_transcript_exact("NM_REAL.2"),
            "from_manifest_without_derived_tx must NOT inject the derived old version (#800)"
        );
        assert!(
            suppressed_cdot.has_transcript_exact("NM_REAL.3"),
            "real sibling record must still be present after injection is suppressed"
        );
    }
}
