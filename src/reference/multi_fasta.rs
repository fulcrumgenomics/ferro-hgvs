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

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use log::warn;

use crate::data::cdot::CdotMapper;
use crate::error::FerroError;
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
    pub transcripts: HashMap<String, SupplementalTranscriptInfo>,
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
    /// Index mapping sequence names to their locations
    index: HashMap<String, FastaIndexEntry>,
    /// Index mapping base accession (without version) to versioned accession
    base_to_versioned: HashMap<String, String>,
    /// Chromosome aliases (e.g., NC_000001 -> chr1)
    aliases: HashMap<String, String>,
    /// Optional cdot transcript metadata (CDS positions, exon coordinates)
    cdot_mapper: Option<CdotMapper>,
    /// Supplemental CDS info for transcripts not in cdot
    supplemental_cds: SupplementalCdsInfo,
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

        let mut index = HashMap::new();
        let mut base_to_versioned = HashMap::new();

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
                        // Map base accession (without version) to versioned
                        if let Some(base) = name.split('.').next() {
                            base_to_versioned.insert(base.to_string(), name.clone());
                        }
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
            index,
            base_to_versioned,
            aliases: build_chromosome_aliases(),
            cdot_mapper: None,
            supplemental_cds: SupplementalCdsInfo::default(),
        })
    }

    /// Create a provider from multiple directories (e.g., transcripts + genome)
    pub fn from_directories<P: AsRef<Path>>(dirs: &[P]) -> Result<Self, FerroError> {
        let mut combined_index = HashMap::new();
        let mut base_to_versioned = HashMap::new();

        for dir in dirs {
            let dir = dir.as_ref();
            if !dir.exists() {
                continue;
            }

            let partial = Self::from_directory(dir)?;
            combined_index.extend(partial.index);
            base_to_versioned.extend(partial.base_to_versioned);
        }

        if combined_index.is_empty() {
            return Err(FerroError::Io {
                msg: "No indexed FASTA files found in any directory".to_string(),
            });
        }

        Ok(Self {
            index: combined_index,
            base_to_versioned,
            aliases: build_chromosome_aliases(),
            cdot_mapper: None,
            supplemental_cds: SupplementalCdsInfo::default(),
        })
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

        // Build base to versioned mapping
        let mut base_to_versioned = HashMap::new();
        for name in index.keys() {
            if let Some(base) = name.split('.').next() {
                base_to_versioned.insert(base.to_string(), name.clone());
            }
        }

        eprintln!("Loaded {} sequences from FASTA", index.len());

        // Load cdot (prefers bincode if available, falls back to JSON)
        let cdot_mapper = CdotMapper::load(cdot_path)?;
        eprintln!(
            "Loaded {} transcripts from cdot",
            cdot_mapper.transcript_count()
        );

        Ok(Self {
            index,
            base_to_versioned,
            aliases: build_chromosome_aliases(),
            cdot_mapper: Some(cdot_mapper),
            supplemental_cds: SupplementalCdsInfo::default(),
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

        // Filter out newlines and take exact length
        let sequence: String = buffer
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .take(seq_len as usize)
            .map(|&b| b as char)
            .collect();

        Ok(sequence.to_uppercase())
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
    /// Strand-aware: cdot stores genome-orientation bases; minus-strand
    /// transcripts get reverse-complemented to transcript orientation here.
    ///
    /// Returns `FerroError::GenomicReferenceNotAvailable` if the cdot-named
    /// contig is missing from the genome FASTA, and `FerroError::InvalidCoordinates`
    /// if any exon's genomic span is degenerate.
    fn synthesize_transcript_from_cdot(
        &self,
        id: &str,
        tx: &crate::data::cdot::CdotTranscript,
    ) -> Result<Transcript, FerroError> {
        use crate::reference::transcript::{Exon as TxExon, GenomeBuild, ManeStatus, Strand};
        use std::sync::OnceLock;

        // Concatenate genome bases per exon in tx-order. cdot stores exons
        // sorted by tx position; the array layout is
        //   [genome_start (HGVS-1-based incl), genome_end (HGVS-1-based excl),
        //    tx_start (0-based incl), tx_end (0-based excl)]
        // so we subtract 1 from both genomic bounds to get the 0-based
        // half-open window `get_genomic_sequence` expects.
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
            sequence.push_str(&bases);
        }

        if matches!(tx.strand, Strand::Minus) {
            sequence = crate::sequence::reverse_complement(&sequence);
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
                genomic_start: Some(e[0] + 1),
                genomic_end: Some(e[1]),
            })
            .collect();

        let (genomic_start, genomic_end) = if !exons.is_empty() {
            let min_start = exons.iter().filter_map(|e| e.genomic_start).min();
            let max_end = exons.iter().filter_map(|e| e.genomic_end).max();
            (min_start, max_end)
        } else {
            (None, None)
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
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            exon_cigars: tx.exon_cigars.clone(),
            cached_introns: OnceLock::new(),
        })
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

impl ReferenceProvider for MultiFastaProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        use crate::reference::transcript::{Exon as TxExon, GenomeBuild, ManeStatus};
        use std::sync::OnceLock;

        if let Some(resolved) = self.resolve_name(id) {
            if let Some(entry) = self.index.get(&resolved) {
                // Get the full sequence
                let sequence = self.get_sequence_from_index(entry, 0, entry.length)?;

                // Try to get metadata from cdot
                let meta = if let Some(ref cdot) = self.cdot_mapper {
                    if let Some(tx) = cdot.get_transcript(&resolved) {
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
                    genome_build: GenomeBuild::GRCh38,
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
        if let Some(ref cdot) = self.cdot_mapper {
            if let Some(tx) = cdot.get_transcript(id) {
                return self.synthesize_transcript_from_cdot(id, tx);
            }
        }

        Err(FerroError::ReferenceNotFound { id: id.to_string() })
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
}

/// Load a FASTA index (.fai) file
///
/// Note: This function canonicalizes the FASTA path to resolve symlinks and
/// provide absolute paths, which helps prevent path traversal issues.
fn load_fai_index<P: AsRef<Path>>(
    fasta_path: P,
    fai_path: P,
) -> Result<HashMap<String, FastaIndexEntry>, FerroError> {
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

    let mut index = HashMap::new();

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
fn count_unique_files(index: &HashMap<String, FastaIndexEntry>) -> usize {
    let files: std::collections::HashSet<_> = index.values().map(|e| &e.file_path).collect();
    files.len()
}

/// Build chromosome name aliases
fn build_chromosome_aliases() -> HashMap<String, String> {
    let mut aliases = HashMap::new();

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

    /// Helper: build a [`MultiFastaProvider`] from a tempdir whose only FASTA is a
    /// genome contig `NC_TEST.1` with the 40 bases `AAAACCCCGGGGTTTT` × 2.5. No
    /// transcript FASTA is created; the synthetic transcript exists only in cdot.
    fn build_provider_with_synthetic_cdot(
        strand: crate::reference::transcript::Strand,
    ) -> MultiFastaProvider {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Transcript};
        use std::sync::OnceLock;

        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("genome.fna");
        let fai_path = dir.path().join("genome.fna.fai");
        {
            let mut f = File::create(&fasta_path).unwrap();
            writeln!(f, ">NC_TEST.1").unwrap();
            // 40 bases of `AAAACCCCGGGGTTTT` repeating
            writeln!(f, "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC").unwrap();
        }
        {
            let mut f = File::create(&fai_path).unwrap();
            // name<TAB>length<TAB>offset<TAB>line_bases<TAB>line_bytes
            // Header ">NC_TEST.1\n" is 11 bytes, so seq offset = 11.
            writeln!(f, "NC_TEST.1\t40\t11\t40\t41").unwrap();
        }
        // Leak the tempdir so the FASTA stays on disk for the lifetime of the
        // provider; the test process exits before any cleanup matters.
        let _kept = Box::leak(Box::new(dir));

        let mut provider = MultiFastaProvider::from_directory(_kept.path()).unwrap();

        // Synthetic transcript:
        //   single exon, tx positions 1..20 (20 bases)
        //   genome positions 11..30 (1-based inclusive on NC_TEST.1)
        //   On + strand, bases = NC_TEST.1[10..30] = "GGTTTTAAAACCCCGGGGTT"
        //   On - strand, bases = reverse_complement of the above
        //                      = "AACCCCGGGGTTTTAAAACC"
        let tx = Transcript {
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
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        provider.cdot_mapper = Some(CdotMapper::from_transcripts(std::iter::once(&tx)));
        provider
    }

    #[test]
    fn test_get_transcript_falls_back_to_cdot_alignment_plus_strand() {
        use crate::reference::transcript::Strand;
        let provider = build_provider_with_synthetic_cdot(Strand::Plus);

        // NM_TEST.1 is NOT in the FASTA index (only NC_TEST.1 is). cdot has the
        // exon alignment. Pre-fix this returns ReferenceNotFound; post-fix it
        // reconstructs the bases from the genome FASTA.
        let tx = provider
            .get_transcript("NM_TEST.1")
            .expect("cdot exon-alignment fallback should fetch transcript");

        assert_eq!(tx.id, "NM_TEST.1");
        assert_eq!(tx.sequence.as_deref(), Some("GGTTTTAAAACCCCGGGGTT"));
        assert!(matches!(tx.strand, Strand::Plus));
        assert_eq!(tx.gene_symbol.as_deref(), Some("TEST"));
        assert_eq!(tx.exons.len(), 1);
    }

    #[test]
    fn test_get_transcript_falls_back_to_cdot_alignment_minus_strand() {
        use crate::reference::transcript::Strand;
        let provider = build_provider_with_synthetic_cdot(Strand::Minus);

        let tx = provider
            .get_transcript("NM_TEST.1")
            .expect("cdot exon-alignment fallback should fetch transcript");

        // Genome bases [10..30] on NC_TEST.1 are `GGTTTTAAAACCCCGGGGTT`;
        // minus-strand transcript orientation is the reverse complement:
        //   complement: CCAAAATTTTGGGGCCCCAA
        //   reversed:   AACCCCGGGGTTTTAAAACC
        assert_eq!(tx.sequence.as_deref(), Some("AACCCCGGGGTTTTAAAACC"));
        assert!(matches!(tx.strand, Strand::Minus));
    }

    #[test]
    fn test_get_transcript_falls_back_with_multi_exon_in_tx_order() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};
        use std::sync::OnceLock;

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
            writeln!(f, "NC_TEST.1\t40\t11\t40\t41").unwrap();
        }
        let _kept = Box::leak(Box::new(dir));
        let mut provider = MultiFastaProvider::from_directory(_kept.path()).unwrap();

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
        let _kept = Box::leak(Box::new(dir));
        let mut provider = MultiFastaProvider::from_directory(_kept.path()).unwrap();

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
}
