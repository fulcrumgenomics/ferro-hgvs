//! FASTA reference sequence provider
//!
//! This module provides functionality to load and query reference sequences
//! from FASTA files.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::error::FerroError;
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::Transcript;

/// Index entry for a sequence in a FASTA file
#[derive(Debug, Clone)]
struct FastaIndexEntry {
    /// Sequence name
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

/// FASTA-based reference sequence provider
///
/// Provides efficient random access to sequences in a FASTA file
/// using an index (either .fai file or built on the fly).
pub struct FastaProvider {
    /// Path to the FASTA file
    path: PathBuf,
    /// Index of sequences
    index: HashMap<String, FastaIndexEntry>,
    /// Chromosome name aliases (for mapping between different naming conventions)
    aliases: HashMap<String, String>,
    /// Transcript data (if loaded)
    transcripts: HashMap<String, Transcript>,
}

impl FastaProvider {
    /// Create a new FASTA provider from a file path
    ///
    /// Will look for an accompanying .fai index file. If not found,
    /// will build an index by scanning the FASTA file.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The file cannot be opened
    /// - The file is gzip-compressed (not supported)
    /// - The file format is invalid
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let path = path.as_ref().to_path_buf();

        // Check for gzip compression
        if is_gzip_file(&path)? {
            return Err(FerroError::Io {
                msg: format!(
                    "FASTA file appears to be gzip-compressed: {}. \
                     Please decompress the file first (e.g., 'gunzip {}' or 'bgzip -d {}').",
                    path.display(),
                    path.display(),
                    path.display()
                ),
            });
        }

        // Try to load .fai index
        let fai_path = path.with_extension("fa.fai");
        let index = if fai_path.exists() {
            load_fai_index(&fai_path)?
        } else {
            let fai_path = PathBuf::from(format!("{}.fai", path.display()));
            if fai_path.exists() {
                load_fai_index(&fai_path)?
            } else {
                // Build index by scanning file
                build_fasta_index(&path)?
            }
        };

        Ok(Self {
            path,
            index,
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        })
    }

    /// Add a chromosome alias
    pub fn add_alias(&mut self, alias: &str, canonical: &str) {
        self.aliases
            .insert(alias.to_string(), canonical.to_string());
    }

    /// Resolve a sequence name to its canonical form
    fn resolve_name(&self, name: &str) -> String {
        // First try direct lookup
        if self.index.contains_key(name) {
            return name.to_string();
        }

        // Try alias
        if let Some(canonical) = self.aliases.get(name) {
            if self.index.contains_key(canonical) {
                return canonical.clone();
            }
        }

        // Try adding/removing "chr" prefix
        let alt_name = if name.starts_with("chr") {
            name.strip_prefix("chr").unwrap().to_string()
        } else {
            format!("chr{}", name)
        };

        if self.index.contains_key(&alt_name) {
            return alt_name;
        }

        // Return original if no match found
        name.to_string()
    }

    /// Get a sequence region from the FASTA file
    fn get_fasta_sequence(&self, name: &str, start: u64, end: u64) -> Result<String, FerroError> {
        let resolved_name = self.resolve_name(name);
        let entry =
            self.index
                .get(&resolved_name)
                .ok_or_else(|| FerroError::ReferenceNotFound {
                    id: name.to_string(),
                })?;

        // Validate coordinates
        if start >= entry.length {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Start position {} exceeds sequence length {} for {}",
                    start, entry.length, name
                ),
            });
        }

        let actual_end = end.min(entry.length);
        if start >= actual_end {
            return Ok(String::new());
        }

        // Validate index entry to prevent divide by zero
        if entry.line_bases == 0 || entry.line_bytes == 0 {
            return Err(FerroError::Io {
                msg: format!(
                    "Invalid FASTA index entry for '{}': line_bases={}, line_bytes={} (must be > 0)",
                    name, entry.line_bases, entry.line_bytes
                ),
            });
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
        let mut file = File::open(&self.path).map_err(|e| FerroError::Io {
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

    /// Get the length of a sequence
    pub fn sequence_length(&self, name: &str) -> Option<u64> {
        let resolved_name = self.resolve_name(name);
        self.index.get(&resolved_name).map(|e| e.length)
    }

    /// Check if a sequence exists
    pub fn has_sequence(&self, name: &str) -> bool {
        let resolved_name = self.resolve_name(name);
        self.index.contains_key(&resolved_name)
    }

    /// Get all sequence names
    pub fn sequence_names(&self) -> impl Iterator<Item = &String> {
        self.index.keys()
    }

    /// Add a transcript to the provider
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.insert(transcript.id.clone(), transcript);
    }
}

impl ReferenceProvider for FastaProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        self.transcripts
            .get(id)
            .cloned()
            .ok_or_else(|| FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        // First try as transcript ID
        if let Some(transcript) = self.transcripts.get(id) {
            // Get sequence from transcript
            let start_idx = start as usize;
            let end_idx = end as usize;
            if end_idx <= transcript.sequence.len() {
                return Ok(transcript.sequence[start_idx..end_idx].to_string());
            }
        }

        // Try as chromosome/contig name
        self.get_fasta_sequence(id, start, end)
    }
}

/// Load a FASTA index (.fai) file
fn load_fai_index<P: AsRef<Path>>(path: P) -> Result<HashMap<String, FastaIndexEntry>, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
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

/// Build a FASTA index by scanning the file
fn build_fasta_index<P: AsRef<Path>>(
    path: P,
) -> Result<HashMap<String, FastaIndexEntry>, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open FASTA file: {}", e),
    })?;
    let mut reader = BufReader::new(file);

    let mut index = HashMap::new();
    let mut current_entry: Option<FastaIndexEntry> = None;
    let mut byte_position = 0u64;
    let mut first_seq_line = true;

    let mut line = String::new();
    loop {
        let line_start = byte_position;
        line.clear();
        let bytes_read = reader.read_line(&mut line).map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if bytes_read == 0 {
            break;
        }

        byte_position += bytes_read as u64;

        if let Some(header) = line.strip_prefix('>') {
            // Save previous entry
            if let Some(entry) = current_entry.take() {
                index.insert(entry.name.clone(), entry);
            }

            // Start new entry
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            current_entry = Some(FastaIndexEntry {
                name,
                length: 0,
                offset: byte_position,
                line_bases: 0,
                line_bytes: 0,
            });
            first_seq_line = true;
        } else if let Some(ref mut entry) = current_entry {
            let seq_len = line.trim_end().len() as u64;
            entry.length += seq_len;

            if first_seq_line && seq_len > 0 {
                entry.offset = line_start;
                entry.line_bases = seq_len;
                entry.line_bytes = bytes_read as u64;
                first_seq_line = false;
            }
        }
    }

    // Save last entry
    if let Some(entry) = current_entry {
        index.insert(entry.name.clone(), entry);
    }

    Ok(index)
}

/// Check if a file is gzip-compressed by reading its magic bytes
///
/// Gzip files start with the magic bytes 0x1f 0x8b
fn is_gzip_file<P: AsRef<Path>>(path: P) -> Result<bool, FerroError> {
    let mut file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open file: {}", e),
    })?;

    let mut magic = [0u8; 2];
    match file.read_exact(&mut magic) {
        Ok(()) => Ok(magic == [0x1f, 0x8b]),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
            // File is too small to be gzipped
            Ok(false)
        }
        Err(e) => Err(FerroError::Io {
            msg: format!("Failed to read file: {}", e),
        }),
    }
}

/// Build default chromosome aliases
fn build_default_aliases() -> HashMap<String, String> {
    let mut aliases = HashMap::new();

    // RefSeq to UCSC/bare chromosome names
    let refseq_chroms = [
        ("NC_000001", "chr1", "1"),
        ("NC_000002", "chr2", "2"),
        ("NC_000003", "chr3", "3"),
        ("NC_000004", "chr4", "4"),
        ("NC_000005", "chr5", "5"),
        ("NC_000006", "chr6", "6"),
        ("NC_000007", "chr7", "7"),
        ("NC_000008", "chr8", "8"),
        ("NC_000009", "chr9", "9"),
        ("NC_000010", "chr10", "10"),
        ("NC_000011", "chr11", "11"),
        ("NC_000012", "chr12", "12"),
        ("NC_000013", "chr13", "13"),
        ("NC_000014", "chr14", "14"),
        ("NC_000015", "chr15", "15"),
        ("NC_000016", "chr16", "16"),
        ("NC_000017", "chr17", "17"),
        ("NC_000018", "chr18", "18"),
        ("NC_000019", "chr19", "19"),
        ("NC_000020", "chr20", "20"),
        ("NC_000021", "chr21", "21"),
        ("NC_000022", "chr22", "22"),
        ("NC_000023", "chrX", "X"),
        ("NC_000024", "chrY", "Y"),
        ("NC_012920", "chrM", "MT"),
    ];

    for (refseq, ucsc, bare) in &refseq_chroms {
        // Map various formats to each other
        aliases.insert(refseq.to_string(), ucsc.to_string());
        aliases.insert(bare.to_string(), ucsc.to_string());
    }

    aliases
}

/// Cached FASTA provider that wraps FastaProvider with an LRU cache
///
/// This wrapper adds caching for sequence lookups to avoid repeated
/// disk I/O for frequently accessed regions.
///
/// # Example
///
/// ```ignore
/// use ferro_hgvs::reference::fasta::{FastaProvider, CachedFastaProvider};
///
/// let provider = FastaProvider::new("reference.fa")?;
/// let cached = CachedFastaProvider::new(provider, 1000);
/// let seq = cached.get_sequence("chr1", 0, 100)?;
/// ```
pub struct CachedFastaProvider {
    /// Underlying FASTA provider
    inner: FastaProvider,
    /// LRU cache for sequence regions: (name, start, end) -> sequence
    cache: crate::cache::LruCache<(String, u64, u64), String>,
}

impl CachedFastaProvider {
    /// Create a new cached FASTA provider
    ///
    /// # Arguments
    ///
    /// * `inner` - The underlying FastaProvider
    /// * `cache_capacity` - Maximum number of sequence regions to cache
    pub fn new(inner: FastaProvider, cache_capacity: usize) -> Self {
        Self {
            inner,
            cache: crate::cache::LruCache::new(cache_capacity),
        }
    }

    /// Create a new cached FASTA provider from a file path
    ///
    /// Uses a default cache capacity of 1000 regions.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        Self::from_path_with_capacity(path, 1000)
    }

    /// Create a new cached FASTA provider from a file path with custom cache capacity
    pub fn from_path_with_capacity<P: AsRef<Path>>(
        path: P,
        cache_capacity: usize,
    ) -> Result<Self, FerroError> {
        let inner = FastaProvider::new(path)?;
        Ok(Self::new(inner, cache_capacity))
    }

    /// Get cache statistics
    pub fn cache_stats(&self) -> crate::cache::CacheStats {
        self.cache.stats()
    }

    /// Clear the cache
    pub fn clear_cache(&self) {
        self.cache.clear();
    }

    /// Get the underlying provider (for transcript operations)
    pub fn inner(&self) -> &FastaProvider {
        &self.inner
    }

    /// Get mutable access to the underlying provider
    pub fn inner_mut(&mut self) -> &mut FastaProvider {
        &mut self.inner
    }

    /// Add a transcript to the underlying provider
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.inner.add_transcript(transcript);
    }

    /// Check if a sequence exists
    pub fn has_sequence(&self, name: &str) -> bool {
        self.inner.has_sequence(name)
    }

    /// Get the length of a sequence
    pub fn sequence_length(&self, name: &str) -> Option<u64> {
        self.inner.sequence_length(name)
    }

    /// Get all sequence names
    pub fn sequence_names(&self) -> impl Iterator<Item = &String> {
        self.inner.sequence_names()
    }
}

impl ReferenceProvider for CachedFastaProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        self.inner.get_transcript(id)
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        // Create cache key
        let key = (id.to_string(), start, end);

        // Check cache first
        if let Some(cached) = self.cache.get(&key) {
            return Ok(cached);
        }

        // Cache miss - fetch from underlying provider
        let sequence = self.inner.get_sequence(id, start, end)?;

        // Store in cache
        self.cache.insert(key, sequence.clone());

        Ok(sequence)
    }
}

/// Memory-mapped FASTA provider for high-performance sequence access
///
/// Uses memory-mapped I/O for efficient access to large FASTA files.
/// The OS handles caching and paging, which is typically more efficient
/// than manual buffering for large files with many random accesses.
///
/// # Example
///
/// ```ignore
/// use ferro_hgvs::reference::fasta::MmapFastaProvider;
///
/// let provider = MmapFastaProvider::new("reference.fa")?;
/// let seq = provider.get_sequence("chr1", 0, 100)?;
/// ```
///
/// # Feature
///
/// This requires the `mmap` feature to be enabled:
/// ```toml
/// ferro-hgvs = { version = "0.1", features = ["mmap"] }
/// ```
#[cfg(feature = "mmap")]
pub struct MmapFastaProvider {
    /// Memory-mapped file
    mmap: memmap2::Mmap,
    /// Index of sequences
    index: HashMap<String, FastaIndexEntry>,
    /// Chromosome name aliases
    aliases: HashMap<String, String>,
    /// Transcript data (if loaded)
    transcripts: HashMap<String, Transcript>,
}

#[cfg(feature = "mmap")]
impl MmapFastaProvider {
    /// Create a new memory-mapped FASTA provider from a file path
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open FASTA file: {}", e),
        })?;

        // Create memory map
        let mmap = unsafe {
            memmap2::Mmap::map(&file).map_err(|e| FerroError::Io {
                msg: format!("Failed to memory-map FASTA file: {}", e),
            })?
        };

        // Try to load accompanying .fai index
        let fai_path = path.with_extension("fa.fai");
        let index = if fai_path.exists() {
            Self::load_fai_index(&fai_path)?
        } else {
            // Build index by scanning the memory-mapped file
            Self::build_index_from_mmap(&mmap)?
        };

        Ok(Self {
            mmap,
            index,
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        })
    }

    /// Load index from .fai file
    fn load_fai_index(path: &Path) -> Result<HashMap<String, FastaIndexEntry>, FerroError> {
        let file = File::open(path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open FAI index: {}", e),
        })?;
        let reader = BufReader::new(file);
        let mut index = HashMap::new();

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read FAI line: {}", e),
            })?;
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let name = parts[0].to_string();
                let length: u64 = parts[1].parse().map_err(|_| FerroError::Io {
                    msg: "Invalid length in FAI index".to_string(),
                })?;
                let offset: u64 = parts[2].parse().map_err(|_| FerroError::Io {
                    msg: "Invalid offset in FAI index".to_string(),
                })?;
                let line_bases: u64 = parts[3].parse().map_err(|_| FerroError::Io {
                    msg: "Invalid line_bases in FAI index".to_string(),
                })?;
                let line_bytes: u64 = parts[4].parse().map_err(|_| FerroError::Io {
                    msg: "Invalid line_bytes in FAI index".to_string(),
                })?;

                index.insert(
                    name.clone(),
                    FastaIndexEntry {
                        name,
                        length,
                        offset,
                        line_bases,
                        line_bytes,
                    },
                );
            }
        }

        Ok(index)
    }

    /// Build index by scanning memory-mapped FASTA data
    fn build_index_from_mmap(
        mmap: &memmap2::Mmap,
    ) -> Result<HashMap<String, FastaIndexEntry>, FerroError> {
        let mut index = HashMap::new();
        let data = &mmap[..];
        let mut i = 0;

        while i < data.len() {
            // Find sequence header
            if data[i] == b'>' {
                let header_start = i + 1;
                // Find end of header line
                while i < data.len() && data[i] != b'\n' {
                    i += 1;
                }
                let header_end = i;
                let header =
                    std::str::from_utf8(&data[header_start..header_end]).map_err(|_| {
                        FerroError::Io {
                            msg: "Invalid UTF-8 in FASTA header".to_string(),
                        }
                    })?;
                let name = header.split_whitespace().next().unwrap_or("").to_string();

                // Move past newline
                i += 1;
                let offset = i as u64;

                // Calculate line parameters from first line
                let first_line_start = i;
                while i < data.len() && data[i] != b'\n' && data[i] != b'>' {
                    i += 1;
                }
                let line_bases = (i - first_line_start) as u64;
                let line_bytes = if i < data.len() && data[i] == b'\n' {
                    line_bases + 1
                } else {
                    line_bases
                };

                // Count total sequence length
                let seq_start = first_line_start;
                while i < data.len() && data[i] != b'>' {
                    i += 1;
                }
                let seq_end = i;

                // Calculate actual sequence length (excluding newlines)
                let mut length = 0u64;
                for &byte in &data[seq_start..seq_end] {
                    if byte != b'\n' && byte != b'\r' {
                        length += 1;
                    }
                }

                index.insert(
                    name.clone(),
                    FastaIndexEntry {
                        name,
                        length,
                        offset,
                        line_bases,
                        line_bytes,
                    },
                );
            } else {
                i += 1;
            }
        }

        Ok(index)
    }

    /// Get sequence from memory-mapped data
    fn get_sequence_from_mmap(
        &self,
        entry: &FastaIndexEntry,
        start: u64,
        end: u64,
    ) -> Result<String, FerroError> {
        if start >= entry.length || end > entry.length || start >= end {
            return Err(FerroError::InvalidCoordinates {
                msg: format!(
                    "Position {}-{} out of range for {} (length {})",
                    start, end, entry.name, entry.length
                ),
            });
        }

        let length = end - start;
        let mut result = String::with_capacity(length as usize);

        // Calculate byte offset for start position
        let full_lines_before = start / entry.line_bases;
        let pos_in_line = start % entry.line_bases;
        let byte_offset = entry.offset + full_lines_before * entry.line_bytes + pos_in_line;

        let data = &self.mmap[..];
        let mut bytes_read = 0u64;
        let mut pos = byte_offset as usize;

        while bytes_read < length && pos < data.len() {
            let byte = data[pos];
            if byte != b'\n' && byte != b'\r' {
                result.push(byte as char);
                bytes_read += 1;
            }
            pos += 1;
        }

        Ok(result.to_uppercase())
    }

    /// Check if a sequence exists
    pub fn has_sequence(&self, name: &str) -> bool {
        // Check direct name first, then alias
        if self.index.contains_key(name) {
            return true;
        }
        if let Some(alias) = self.aliases.get(name) {
            return self.index.contains_key(alias);
        }
        false
    }

    /// Get the length of a sequence
    pub fn sequence_length(&self, name: &str) -> Option<u64> {
        // Check direct name first, then alias
        if let Some(entry) = self.index.get(name) {
            return Some(entry.length);
        }
        if let Some(alias) = self.aliases.get(name) {
            return self.index.get(alias).map(|e| e.length);
        }
        None
    }

    /// Get all sequence names
    pub fn sequence_names(&self) -> impl Iterator<Item = &String> {
        self.index.keys()
    }

    /// Add a transcript to the provider
    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.insert(transcript.id.clone(), transcript);
    }
}

#[cfg(feature = "mmap")]
impl ReferenceProvider for MmapFastaProvider {
    fn get_transcript(&self, id: &str) -> Result<Transcript, FerroError> {
        // Handle versioned and unversioned lookups
        if let Some(tx) = self.transcripts.get(id) {
            return Ok(tx.clone());
        }

        // Try without version
        let base_id = id.split('.').next().unwrap_or(id);
        for (key, tx) in &self.transcripts {
            if key.starts_with(base_id) {
                return Ok(tx.clone());
            }
        }

        Err(FerroError::ReferenceNotFound { id: id.to_string() })
    }

    fn get_sequence(&self, id: &str, start: u64, end: u64) -> Result<String, FerroError> {
        // Check direct name first, then alias
        let entry = if let Some(e) = self.index.get(id) {
            e
        } else if let Some(alias) = self.aliases.get(id) {
            self.index
                .get(alias)
                .ok_or_else(|| FerroError::ReferenceNotFound { id: id.to_string() })?
        } else {
            return Err(FerroError::ReferenceNotFound { id: id.to_string() });
        };

        self.get_sequence_from_mmap(entry, start, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_default_aliases() {
        let aliases = build_default_aliases();
        assert_eq!(aliases.get("NC_000001"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("1"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("X"), Some(&"chrX".to_string()));
    }

    #[test]
    fn test_resolve_name_with_chr_prefix() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: {
                let mut idx = HashMap::new();
                idx.insert(
                    "chr1".to_string(),
                    FastaIndexEntry {
                        name: "chr1".to_string(),
                        length: 1000,
                        offset: 0,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx
            },
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        };

        // Direct match
        assert_eq!(provider.resolve_name("chr1"), "chr1");
        // Without chr prefix
        assert_eq!(provider.resolve_name("1"), "chr1");
    }

    #[test]
    fn test_fai_format() {
        // FAI format: name, length, offset, line_bases, line_bytes
        let fai_line = "chr1\t248956422\t6\t50\t51";
        let fields: Vec<&str> = fai_line.split('\t').collect();

        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1].parse::<u64>().unwrap(), 248956422);
        assert_eq!(fields[2].parse::<u64>().unwrap(), 6);
        assert_eq!(fields[3].parse::<u64>().unwrap(), 50);
        assert_eq!(fields[4].parse::<u64>().unwrap(), 51);
    }

    // ===== FastaProvider Extended Tests =====

    #[test]
    fn test_fasta_provider_add_alias() {
        let mut provider = FastaProvider {
            path: PathBuf::new(),
            index: HashMap::new(),
            aliases: HashMap::new(),
            transcripts: HashMap::new(),
        };

        provider.add_alias("1", "chr1");
        assert_eq!(provider.aliases.get("1"), Some(&"chr1".to_string()));
    }

    #[test]
    fn test_fasta_provider_has_sequence() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: {
                let mut idx = HashMap::new();
                idx.insert(
                    "chr1".to_string(),
                    FastaIndexEntry {
                        name: "chr1".to_string(),
                        length: 1000,
                        offset: 0,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx
            },
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        };

        assert!(provider.has_sequence("chr1"));
        assert!(provider.has_sequence("1")); // Alias
        assert!(!provider.has_sequence("chr99"));
    }

    #[test]
    fn test_fasta_provider_sequence_length() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: {
                let mut idx = HashMap::new();
                idx.insert(
                    "chr1".to_string(),
                    FastaIndexEntry {
                        name: "chr1".to_string(),
                        length: 248956422,
                        offset: 0,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx
            },
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        };

        assert_eq!(provider.sequence_length("chr1"), Some(248956422));
        assert_eq!(provider.sequence_length("1"), Some(248956422)); // Alias
        assert_eq!(provider.sequence_length("chr99"), None);
    }

    #[test]
    fn test_fasta_provider_sequence_names() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: {
                let mut idx = HashMap::new();
                idx.insert(
                    "chr1".to_string(),
                    FastaIndexEntry {
                        name: "chr1".to_string(),
                        length: 1000,
                        offset: 0,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx.insert(
                    "chr2".to_string(),
                    FastaIndexEntry {
                        name: "chr2".to_string(),
                        length: 2000,
                        offset: 1000,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx
            },
            aliases: HashMap::new(),
            transcripts: HashMap::new(),
        };

        let names: Vec<_> = provider.sequence_names().collect();
        assert_eq!(names.len(), 2);
    }

    #[test]
    fn test_fasta_provider_add_transcript() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = FastaProvider {
            path: PathBuf::new(),
            index: HashMap::new(),
            aliases: HashMap::new(),
            transcripts: HashMap::new(),
        };

        let transcript = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCATGCATGC".to_string(),
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        provider.add_transcript(transcript);
        assert!(provider.transcripts.contains_key("NM_000088.3"));
    }

    #[test]
    fn test_fasta_provider_get_transcript() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = FastaProvider {
            path: PathBuf::new(),
            index: HashMap::new(),
            aliases: HashMap::new(),
            transcripts: HashMap::new(),
        };

        let transcript = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCATGCATGC".to_string(),
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        provider.add_transcript(transcript);

        let result = provider.get_transcript("NM_000088.3");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().gene_symbol, Some("COL1A1".to_string()));

        let result = provider.get_transcript("NM_NONEXISTENT");
        assert!(result.is_err());
    }

    #[test]
    fn test_fasta_provider_get_sequence_from_transcript() {
        use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand};
        use std::sync::OnceLock;

        let mut provider = FastaProvider {
            path: PathBuf::new(),
            index: HashMap::new(),
            aliases: HashMap::new(),
            transcripts: HashMap::new(),
        };

        let transcript = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: "ATGCATGCATGC".to_string(),
            cds_start: Some(1),
            cds_end: Some(12),
            exons: vec![Exon::new(1, 1, 12)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            cached_introns: OnceLock::new(),
        };

        provider.add_transcript(transcript);

        let result = provider.get_sequence("NM_000088.3", 0, 4);
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), "ATGC");
    }

    // ===== resolve_name Extended Tests =====

    #[test]
    fn test_resolve_name_without_chr_prefix() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: {
                let mut idx = HashMap::new();
                idx.insert(
                    "1".to_string(),
                    FastaIndexEntry {
                        name: "1".to_string(),
                        length: 1000,
                        offset: 0,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx
            },
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        };

        // Request with chr prefix, index without
        assert_eq!(provider.resolve_name("chr1"), "1");
        // Direct match
        assert_eq!(provider.resolve_name("1"), "1");
    }

    #[test]
    fn test_resolve_name_no_match() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: HashMap::new(),
            aliases: HashMap::new(),
            transcripts: HashMap::new(),
        };

        // No match should return original
        assert_eq!(provider.resolve_name("chrZ"), "chrZ");
    }

    #[test]
    fn test_resolve_name_refseq_alias() {
        let provider = FastaProvider {
            path: PathBuf::new(),
            index: {
                let mut idx = HashMap::new();
                idx.insert(
                    "chr1".to_string(),
                    FastaIndexEntry {
                        name: "chr1".to_string(),
                        length: 1000,
                        offset: 0,
                        line_bases: 50,
                        line_bytes: 51,
                    },
                );
                idx
            },
            aliases: build_default_aliases(),
            transcripts: HashMap::new(),
        };

        // RefSeq accession should resolve to chr1
        assert_eq!(provider.resolve_name("NC_000001"), "chr1");
    }

    // ===== build_default_aliases Tests =====

    #[test]
    fn test_build_default_aliases_chromosomes() {
        let aliases = build_default_aliases();

        // Test RefSeq mappings
        assert_eq!(aliases.get("NC_000001"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("NC_000023"), Some(&"chrX".to_string()));
        assert_eq!(aliases.get("NC_000024"), Some(&"chrY".to_string()));
        assert_eq!(aliases.get("NC_012920"), Some(&"chrM".to_string()));

        // Test bare chromosome mappings
        assert_eq!(aliases.get("1"), Some(&"chr1".to_string()));
        assert_eq!(aliases.get("X"), Some(&"chrX".to_string()));
        assert_eq!(aliases.get("MT"), Some(&"chrM".to_string()));
    }

    // ===== FastaIndexEntry Tests =====

    #[test]
    fn test_fasta_index_entry_debug() {
        let entry = FastaIndexEntry {
            name: "chr1".to_string(),
            length: 248956422,
            offset: 6,
            line_bases: 50,
            line_bytes: 51,
        };

        let debug_str = format!("{:?}", entry);
        assert!(debug_str.contains("chr1"));
        assert!(debug_str.contains("248956422"));
    }

    // ===== I/O Tests with Temporary Files =====

    #[test]
    fn test_build_fasta_index() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut fasta_file = NamedTempFile::new().unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(
            fasta_file,
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        )
        .unwrap();
        writeln!(
            fasta_file,
            "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        )
        .unwrap();
        writeln!(fasta_file, ">chr2").unwrap();
        writeln!(
            fasta_file,
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        )
        .unwrap();
        fasta_file.flush().unwrap();

        let index = build_fasta_index(fasta_file.path()).unwrap();

        assert!(index.contains_key("chr1"));
        assert!(index.contains_key("chr2"));

        let chr1 = index.get("chr1").unwrap();
        assert_eq!(chr1.length, 100); // 50 + 50 bases
        assert_eq!(chr1.line_bases, 50);
        assert_eq!(chr1.line_bytes, 51); // 50 + newline

        let chr2 = index.get("chr2").unwrap();
        assert_eq!(chr2.length, 50);
    }

    #[test]
    fn test_load_fai_index() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut fai_file = NamedTempFile::new().unwrap();
        writeln!(fai_file, "chr1\t248956422\t6\t50\t51").unwrap();
        writeln!(fai_file, "chr2\t242193529\t254235640\t50\t51").unwrap();
        fai_file.flush().unwrap();

        let index = load_fai_index(fai_file.path()).unwrap();

        assert!(index.contains_key("chr1"));
        assert!(index.contains_key("chr2"));

        let chr1 = index.get("chr1").unwrap();
        assert_eq!(chr1.length, 248956422);
        assert_eq!(chr1.offset, 6);
        assert_eq!(chr1.line_bases, 50);
        assert_eq!(chr1.line_bytes, 51);
    }

    #[test]
    fn test_load_fai_index_invalid_lines() {
        use std::io::Write;
        use tempfile::NamedTempFile;

        let mut fai_file = NamedTempFile::new().unwrap();
        writeln!(fai_file, "chr1\t248956422\t6\t50\t51").unwrap();
        writeln!(fai_file, "invalid_line_too_few_fields").unwrap(); // Should be skipped
        writeln!(fai_file, "chr2\t100\t200\t50\t51").unwrap();
        fai_file.flush().unwrap();

        let index = load_fai_index(fai_file.path()).unwrap();

        assert!(index.contains_key("chr1"));
        assert!(index.contains_key("chr2"));
        assert_eq!(index.len(), 2); // Invalid line skipped
    }

    #[test]
    fn test_fasta_provider_new_with_fai() {
        use std::io::Write;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(
            fasta_file,
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        )
        .unwrap();

        // Create FAI file
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "chr1\t50\t6\t50\t51").unwrap();

        let provider = FastaProvider::new(&fasta_path).unwrap();
        assert!(provider.has_sequence("chr1"));
        assert_eq!(provider.sequence_length("chr1"), Some(50));
    }

    #[test]
    fn test_fasta_provider_new_without_fai() {
        use std::io::Write;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");

        // Create FASTA file only (no FAI)
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(
            fasta_file,
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        )
        .unwrap();

        let provider = FastaProvider::new(&fasta_path).unwrap();
        assert!(provider.has_sequence("chr1"));
        assert_eq!(provider.sequence_length("chr1"), Some(50));
    }

    #[test]
    fn test_get_fasta_sequence() {
        use std::io::Write;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(
            fasta_file,
            "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT"
        )
        .unwrap();
        writeln!(
            fasta_file,
            "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
        )
        .unwrap();

        // Create FAI file
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "chr1\t100\t6\t50\t51").unwrap();

        let provider = FastaProvider::new(&fasta_path).unwrap();

        // Get first 10 bases
        let seq = provider.get_fasta_sequence("chr1", 0, 10).unwrap();
        assert_eq!(seq, "ATGCATGCAT");

        // Get bases spanning a line break
        let seq = provider.get_fasta_sequence("chr1", 45, 55).unwrap();
        assert_eq!(seq.len(), 10);

        // Get with alias
        let seq = provider.get_fasta_sequence("1", 0, 4).unwrap();
        assert_eq!(seq, "ATGC");
    }

    #[test]
    fn test_get_fasta_sequence_errors() {
        use std::io::Write;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        // Create FAI file
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "chr1\t10\t6\t10\t11").unwrap();

        let provider = FastaProvider::new(&fasta_path).unwrap();

        // Unknown sequence
        let result = provider.get_fasta_sequence("chrZ", 0, 10);
        assert!(result.is_err());

        // Start beyond length
        let result = provider.get_fasta_sequence("chr1", 100, 110);
        assert!(result.is_err());

        // Start >= end (empty result)
        let seq = provider.get_fasta_sequence("chr1", 5, 5).unwrap();
        assert!(seq.is_empty());
    }

    #[test]
    fn test_fasta_provider_get_sequence_from_fasta() {
        use std::io::Write;
        use tempfile::tempdir;

        let dir = tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        let fai_path = dir.path().join("test.fa.fai");

        // Create FASTA file
        let mut fasta_file = File::create(&fasta_path).unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(fasta_file, "ATGCATGCAT").unwrap();

        // Create FAI file
        let mut fai_file = File::create(&fai_path).unwrap();
        writeln!(fai_file, "chr1\t10\t6\t10\t11").unwrap();

        let provider = FastaProvider::new(&fasta_path).unwrap();

        // Should fall through to FASTA since no transcript with this ID
        let seq = provider.get_sequence("chr1", 0, 4).unwrap();
        assert_eq!(seq, "ATGC");
    }

    #[test]
    fn test_fasta_index_entry_clone() {
        let entry = FastaIndexEntry {
            name: "chr1".to_string(),
            length: 1000,
            offset: 0,
            line_bases: 50,
            line_bytes: 51,
        };

        let cloned = entry.clone();
        assert_eq!(cloned.name, entry.name);
        assert_eq!(cloned.length, entry.length);
    }
}
