//! cdot JSON parsing for transcript alignments.
//!
//! This module handles loading and parsing cdot JSON files which contain
//! transcript-to-genome alignments used for coordinate mapping.
//!
//! Supports both the real cdot format (with nested genome_builds) and
//! a simplified flat format for testing.
//!
//! # Coordinate Systems
//!
//! **IMPORTANT**: cdot uses MIXED coordinate systems!
//!
//! | Field | Basis | Format |
//! |-------|-------|--------|
//! | Genomic coordinates (`genome_start`, `genome_end`) | 0-based | Half-open `[start, end)` |
//! | Transcript coordinates (`tx_start`, `tx_end`) | 1-based | See note below |
//! | CDS coordinates (`cds_start`, `cds_end`) | 0-based | Half-open `[start, end)` |
//!
//! **Note on transcript coordinates**: The cdot format uses 1-based transcript
//! positions (fixed in commit 944a4e9). This is a common source of bugs when
//! developers assume all coordinates are 0-based.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::reference::Strand;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// cdot transcript alignment data (normalized internal representation).
///
/// # Coordinate Systems
/// - Genomic: 0-based half-open `[start, end)`
/// - Transcript (tx_start/tx_end in exons): 1-based
/// - CDS: 0-based half-open `[start, end)` in transcript space
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdotTranscript {
    /// Gene name/symbol (e.g., "BRCA1").
    #[serde(default)]
    pub gene_name: Option<String>,

    /// Genomic contig/chromosome (e.g., "NC_000017.11").
    #[serde(alias = "contig", alias = "chr")]
    pub contig: String,

    /// Strand (+ or -).
    #[serde(deserialize_with = "deserialize_strand")]
    pub strand: Strand,

    /// Exon alignments: [genome_start(0-based), genome_end(0-based excl), tx_start(1-based), tx_end(1-based)].
    pub exons: Vec<[u64; 4]>,

    /// Start of CDS in transcript coordinates (0-based, inclusive).
    #[serde(default)]
    pub cds_start: Option<u64>,

    /// End of CDS in transcript coordinates (0-based, exclusive).
    #[serde(default)]
    pub cds_end: Option<u64>,

    /// Per-exon CIGAR alignment data from cdot gap info (GFF3 Gap format).
    /// Indexed in the same order as `exons` (sorted by tx position).
    #[serde(skip)]
    pub exon_cigars: Vec<Option<Vec<CigarOp>>>,

    /// Gene ID (e.g., "HGNC:1100").
    #[serde(default)]
    pub gene_id: Option<String>,

    /// Protein accession (e.g., "NP_000079.2").
    #[serde(default)]
    pub protein: Option<String>,
}

/// Raw cdot transcript entry as it appears in the JSON file.
/// The real cdot format nests genome-specific data under genome_builds.
#[derive(Debug, Clone, Deserialize)]
#[allow(dead_code)] // Fields used for serde deserialization
struct RawCdotTranscript {
    #[serde(default)]
    gene_name: Option<String>,
    #[serde(default)]
    gene_version: Option<String>,
    #[serde(default)]
    biotype: Option<Vec<String>>,
    #[serde(default)]
    genome_builds: Option<HashMap<String, RawGenomeBuild>>,
    /// CDS start in transcript coordinates (0-indexed)
    #[serde(default)]
    start_codon: Option<u64>,
    /// CDS end in transcript coordinates (0-indexed, exclusive)
    #[serde(default)]
    stop_codon: Option<u64>,
    // Flat format fields (for backwards compatibility)
    #[serde(default)]
    contig: Option<String>,
    #[serde(default)]
    strand: Option<String>,
    #[serde(default)]
    exons: Option<Vec<Vec<serde_json::Value>>>,
    #[serde(default)]
    cds_start: Option<u64>,
    #[serde(default)]
    cds_end: Option<u64>,
}

/// Genome build specific data in cdot format.
#[derive(Debug, Clone, Deserialize)]
struct RawGenomeBuild {
    contig: String,
    #[serde(default)]
    strand: Option<String>,
    /// Exons in cdot format: [genomic_start, genomic_end, exon_num, tx_start, tx_end, gap_info]
    exons: Vec<Vec<serde_json::Value>>,
    /// CDS start in genomic coordinates
    #[serde(default)]
    cds_start: Option<u64>,
    /// CDS end in genomic coordinates
    #[serde(default)]
    cds_end: Option<u64>,
}

impl RawCdotTranscript {
    /// Convert raw transcript to internal format for a specific genome build.
    fn to_transcript(&self, genome_build: &str) -> Option<CdotTranscript> {
        // Try genome_builds first (real cdot format)
        if let Some(builds) = &self.genome_builds {
            if let Some(build) = builds.get(genome_build) {
                return self.from_genome_build(build);
            }
        }

        // Fall back to flat format
        self.from_flat_format()
    }

    /// Convert from nested genome_builds format.
    #[allow(clippy::wrong_self_convention)]
    fn from_genome_build(&self, build: &RawGenomeBuild) -> Option<CdotTranscript> {
        let strand = parse_strand(build.strand.as_deref().unwrap_or("+"))?;

        // Parse exons: [genomic_start, genomic_end, exon_num, tx_start, tx_end, gap_info]
        let mut exon_pairs: Vec<([u64; 4], Option<Vec<CigarOp>>)> = build
            .exons
            .iter()
            .filter_map(|e| {
                if e.len() >= 5 {
                    let exon = [
                        e[0].as_u64()?,
                        e[1].as_u64()?,
                        e[3].as_u64()?, // tx_start is at index 3
                        e[4].as_u64()?, // tx_end is at index 4
                    ];
                    // Parse CIGAR gap info from index 5 if present
                    let cigar = if e.len() > 5 {
                        e[5].as_str().and_then(|s| match parse_cigar(s) {
                            Ok(ops) => Some(ops),
                            Err(err) => {
                                log::warn!("Malformed CIGAR string '{}': {}", s, err);
                                None
                            }
                        })
                    } else {
                        None
                    };
                    Some((exon, cigar))
                } else {
                    None
                }
            })
            .collect();

        // Sort exons by transcript position, keeping CIGARs in sync
        exon_pairs.sort_by_key(|(e, _)| e[2]);
        let (exons, exon_cigars): (Vec<[u64; 4]>, Vec<Option<Vec<CigarOp>>>) =
            exon_pairs.into_iter().unzip();

        // Use transcript-level CDS coordinates (start_codon/stop_codon) if available
        // These are 0-indexed and more reliable than converting from genomic coords
        let (cds_start, cds_end) = if self.start_codon.is_some() || self.stop_codon.is_some() {
            // Use transcript-level coordinates directly (0-indexed)
            (self.start_codon, self.stop_codon)
        } else {
            // Fall back to converting genomic CDS coordinates to transcript coordinates
            match (build.cds_start, build.cds_end) {
                (Some(g_start), Some(g_end)) => {
                    // Find transcript positions for genomic CDS boundaries
                    let tx_cds_start = genomic_to_tx_pos(&exons, g_start, strand);
                    let tx_cds_end = genomic_to_tx_pos(&exons, g_end, strand);

                    match (tx_cds_start, tx_cds_end) {
                        (Some(s), Some(e)) => {
                            // CDS end is exclusive, so add 1
                            // Use saturating_add to prevent overflow at u64::MAX
                            let (start, end) = if s < e {
                                (s, e.saturating_add(1))
                            } else {
                                (e, s.saturating_add(1))
                            };
                            (Some(start), Some(end))
                        }
                        _ => (None, None),
                    }
                }
                _ => (None, None),
            }
        };

        Some(CdotTranscript {
            gene_name: self.gene_name.clone(),
            contig: build.contig.clone(),
            strand,
            exons,
            exon_cigars,
            cds_start,
            cds_end,
            gene_id: None,
            protein: None,
        })
    }

    /// Convert from flat format (for backwards compatibility).
    #[allow(clippy::wrong_self_convention)]
    fn from_flat_format(&self) -> Option<CdotTranscript> {
        let contig = self.contig.as_ref()?.clone();
        let strand = parse_strand(self.strand.as_deref().unwrap_or("+"))?;

        let exons: Vec<[u64; 4]> = self
            .exons
            .as_ref()?
            .iter()
            .filter_map(|e| {
                if e.len() >= 4 {
                    Some([
                        e[0].as_u64()?,
                        e[1].as_u64()?,
                        e[2].as_u64()?,
                        e[3].as_u64()?,
                    ])
                } else {
                    None
                }
            })
            .collect();

        Some(CdotTranscript {
            gene_name: self.gene_name.clone(),
            contig,
            strand,
            exon_cigars: Vec::new(),
            exons,
            cds_start: self.cds_start,
            cds_end: self.cds_end,
            gene_id: None,
            protein: None,
        })
    }
}

/// Parse strand string to Strand enum.
fn parse_strand(s: &str) -> Option<Strand> {
    match s {
        "+" | "1" | "plus" => Some(Strand::Plus),
        "-" | "-1" | "minus" => Some(Strand::Minus),
        _ => None,
    }
}

/// Convert genomic position to transcript position using exon data.
fn genomic_to_tx_pos(exons: &[[u64; 4]], genomic_pos: u64, strand: Strand) -> Option<u64> {
    for e in exons {
        let (g_start, g_end, tx_start, _tx_end) = (e[0], e[1], e[2], e[3]);
        if genomic_pos >= g_start && genomic_pos < g_end {
            let offset = match strand {
                Strand::Plus => genomic_pos - g_start,
                Strand::Minus => g_end - 1 - genomic_pos,
            };
            return Some(tx_start + offset);
        }
    }
    // Position not in an exon - find closest exon boundary
    // This handles CDS boundaries that might be at intron edges
    for e in exons {
        let (g_start, g_end, tx_start, tx_end) = (e[0], e[1], e[2], e[3]);
        if genomic_pos == g_end {
            return Some(tx_end);
        }
        if genomic_pos == g_start.saturating_sub(1) {
            return Some(tx_start.saturating_sub(1));
        }
    }
    None
}

/// Deserialize strand from string.
fn deserialize_strand<'de, D>(deserializer: D) -> Result<Strand, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    match s.as_str() {
        "+" | "1" | "plus" => Ok(Strand::Plus),
        "-" | "-1" | "minus" => Ok(Strand::Minus),
        _ => Err(serde::de::Error::custom(format!("Invalid strand: {}", s))),
    }
}

/// Parsed exon with named fields for clarity.
///
/// # Coordinate Systems
/// - Genomic coordinates: 0-based half-open `[genome_start, genome_end)`
/// - Transcript coordinates: 1-based (NOT 0-based - common source of bugs!)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Exon {
    /// Exon number (1-based).
    pub number: u32,
    /// Start position in genomic coordinates (0-based, inclusive).
    pub genome_start: u64,
    /// End position in genomic coordinates (0-based, exclusive).
    pub genome_end: u64,
    /// Start position in transcript coordinates (1-based, inclusive).
    pub tx_start: u64,
    /// End position in transcript coordinates (1-based, exclusive for range).
    pub tx_end: u64,
}

impl Exon {
    /// Create an exon from cdot array format.
    pub fn from_array(number: u32, arr: [u64; 4]) -> Self {
        Self {
            number,
            genome_start: arr[0],
            genome_end: arr[1],
            tx_start: arr[2],
            tx_end: arr[3],
        }
    }

    /// Get the length of this exon.
    pub fn len(&self) -> u64 {
        self.genome_end - self.genome_start
    }

    /// Check if the exon is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Check if a genomic position is within this exon.
    pub fn contains_genome_pos(&self, pos: u64) -> bool {
        pos >= self.genome_start && pos < self.genome_end
    }

    /// Check if a transcript position is within this exon.
    pub fn contains_tx_pos(&self, pos: u64) -> bool {
        pos >= self.tx_start && pos < self.tx_end
    }
}

impl CdotTranscript {
    /// Get exons as structured objects.
    pub fn get_exons(&self) -> Vec<Exon> {
        self.exons
            .iter()
            .enumerate()
            .map(|(i, arr)| Exon::from_array((i + 1) as u32, *arr))
            .collect()
    }

    /// Get total transcript length (sum of exon lengths).
    pub fn transcript_length(&self) -> u64 {
        self.exons.iter().map(|e| e[3] - e[2]).sum()
    }

    /// Get CDS length if CDS is defined.
    pub fn cds_length(&self) -> Option<u64> {
        match (self.cds_start, self.cds_end) {
            (Some(start), Some(end)) => Some(end - start),
            _ => None,
        }
    }

    /// Find exon containing a transcript position.
    pub fn exon_for_tx_pos(&self, tx_pos: u64) -> Option<Exon> {
        for (i, arr) in self.exons.iter().enumerate() {
            if tx_pos >= arr[2] && tx_pos < arr[3] {
                return Some(Exon::from_array((i + 1) as u32, *arr));
            }
        }
        None
    }

    /// Find exon containing a genomic position.
    pub fn exon_for_genome_pos(&self, genome_pos: u64) -> Option<Exon> {
        for (i, arr) in self.exons.iter().enumerate() {
            if genome_pos >= arr[0] && genome_pos < arr[1] {
                return Some(Exon::from_array((i + 1) as u32, *arr));
            }
        }
        None
    }

    /// Convert transcript position to genomic position.
    pub fn tx_to_genome(&self, tx_pos: u64) -> Option<u64> {
        let exon = self.exon_for_tx_pos(tx_pos)?;
        let offset = tx_pos - exon.tx_start;

        match self.strand {
            Strand::Plus => Some(exon.genome_start + offset),
            Strand::Minus => Some(exon.genome_end - 1 - offset),
        }
    }

    /// Convert genomic position to transcript position.
    pub fn genome_to_tx(&self, genome_pos: u64) -> Option<u64> {
        let exon = self.exon_for_genome_pos(genome_pos)?;

        match self.strand {
            Strand::Plus => {
                let offset = genome_pos - exon.genome_start;
                Some(exon.tx_start + offset)
            }
            Strand::Minus => {
                let offset = exon.genome_end - 1 - genome_pos;
                Some(exon.tx_start + offset)
            }
        }
    }

    /// Convert CDS position (1-based) to transcript position (0-based).
    pub fn cds_to_tx(&self, cds_pos: i64) -> Option<u64> {
        let cds_start = self.cds_start?;

        if cds_pos > 0 {
            // Normal CDS position
            Some(cds_start + (cds_pos as u64 - 1))
        } else if cds_pos < 0 {
            // 5' UTR position (e.g., c.-100)
            let offset = (-cds_pos) as u64;
            if offset <= cds_start {
                Some(cds_start - offset)
            } else {
                None // Beyond transcript start
            }
        } else {
            None // Position 0 is not valid
        }
    }

    /// Convert transcript position (0-based) to CDS position (1-based).
    pub fn tx_to_cds(&self, tx_pos: u64) -> Option<CdsPosition> {
        let cds_start = self.cds_start?;
        let cds_end = self.cds_end?;

        if tx_pos < cds_start {
            // 5' UTR
            let offset = cds_start - tx_pos;
            Some(CdsPosition::FivePrimeUtr(offset as i64))
        } else if tx_pos < cds_end {
            // CDS
            let cds_pos = tx_pos - cds_start + 1;
            Some(CdsPosition::Cds(cds_pos as i64))
        } else {
            // 3' UTR
            let offset = tx_pos - cds_end + 1;
            Some(CdsPosition::ThreePrimeUtr(offset as i64))
        }
    }

    /// Get the intron number for a genomic position.
    /// Returns None if the position is within an exon.
    pub fn intron_for_genome_pos(&self, genome_pos: u64) -> Option<(u32, i64)> {
        let exons = self.get_exons();

        for i in 0..exons.len().saturating_sub(1) {
            let current = &exons[i];
            let next = &exons[i + 1];

            let (intron_start, intron_end) = match self.strand {
                Strand::Plus => (current.genome_end, next.genome_start),
                Strand::Minus => (next.genome_end, current.genome_start),
            };

            if genome_pos >= intron_start && genome_pos < intron_end {
                let intron_num = (i + 1) as u32;
                let offset = match self.strand {
                    Strand::Plus => (genome_pos as i64) - (current.genome_end as i64),
                    Strand::Minus => (current.genome_start as i64) - (genome_pos as i64),
                };
                // Adjust for 0-based to intronic offset convention
                let intronic_offset = if offset >= 0 { offset + 1 } else { offset };
                return Some((intron_num, intronic_offset));
            }
        }

        None
    }
}

/// CDS position type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CdsPosition {
    /// Position in 5' UTR (e.g., c.-100 would be FivePrimeUtr(100)).
    FivePrimeUtr(i64),
    /// Normal CDS position (1-based).
    Cds(i64),
    /// Position in 3' UTR (e.g., c.*100 would be ThreePrimeUtr(100)).
    ThreePrimeUtr(i64),
}

/// A single CIGAR operation from a GFF3 Gap attribute string.
///
/// cdot uses GFF3 Gap format (letter-first, space-separated), e.g. `M185 I3 M250`,
/// which differs from SAM CIGAR format (`185M3I250M`).
///
/// - `M` (Match): alignment match — bases align between transcript and genome.
/// - `I` (Insertion): bases present in the transcript but not the genome.
/// - `D` (Deletion): bases present in the genome but not the transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// Alignment match of `n` bases.
    Match(u64),
    /// Insertion of `n` bases in the transcript (not in genome).
    Insertion(u64),
    /// Deletion of `n` bases from the transcript (present in genome).
    Deletion(u64),
}

/// Parse a GFF3 Gap attribute string into a sequence of CIGAR operations.
///
/// The format is space-separated tokens, each starting with an operation letter
/// (`M`, `I`, or `D`) followed by a length, e.g. `"M185 I3 M250"`.
///
/// Returns an empty vector for empty or whitespace-only input.
///
/// # Errors
///
/// Returns an error if a token has an unrecognized operation letter or a non-numeric length.
pub fn parse_cigar(cigar_str: &str) -> Result<Vec<CigarOp>, FerroError> {
    let trimmed = cigar_str.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }

    trimmed
        .split_whitespace()
        .map(|token| {
            if token.len() < 2 {
                return Err(FerroError::parse(
                    0,
                    format!("Invalid CIGAR token (too short): \'{token}\'"),
                ));
            }
            let (op_char, len_str) = token.split_at(1);
            let length: u64 = len_str.parse().map_err(|_| {
                FerroError::parse(0, format!("Invalid CIGAR length in token: \'{token}\'"))
            })?;
            match op_char {
                "M" => Ok(CigarOp::Match(length)),
                "I" => Ok(CigarOp::Insertion(length)),
                "D" => Ok(CigarOp::Deletion(length)),
                _ => Err(FerroError::parse(
                    0,
                    format!("Unknown CIGAR operation \'{op_char}\' in token: \'{token}\'"),
                )),
            }
        })
        .collect()
}

/// Compute the cumulative insertion offset at a given 1-based transcript position.
///
/// Walks through the CIGAR operations and counts insertion bases that occur
/// before `tx_pos`. Insertions add bases to the transcript that do not exist
/// in the genome, so CDS numbering (which follows the genome) must account
/// for them.
///
/// Returns the total number of insertion bases encountered before `tx_pos`.
pub fn cumulative_insertion_offset(ops: &[CigarOp], tx_pos: u64) -> u64 {
    let mut current_tx: u64 = 0;
    let mut cumulative: u64 = 0;

    for op in ops {
        match op {
            CigarOp::Match(len) => {
                current_tx += len;
            }
            CigarOp::Insertion(len) => {
                // If the entire insertion is before tx_pos, count it all
                if current_tx + len <= tx_pos {
                    cumulative += len;
                }
                current_tx += len;
            }
            CigarOp::Deletion(_) => {
                // Deletions don't advance the transcript position
            }
        }
        if current_tx >= tx_pos {
            break;
        }
    }

    cumulative
}

/// Raw cdot JSON file structure (as it appears on disk).
#[derive(Debug, Clone, Deserialize)]
#[allow(dead_code)] // Fields used for serde deserialization
struct RawCdotFile {
    /// Transcripts indexed by accession (raw format).
    transcripts: HashMap<String, RawCdotTranscript>,

    /// Genome builds available (list, e.g., ["GRCh38"])
    #[serde(default)]
    genome_builds: Option<Vec<String>>,

    /// cdot version
    #[serde(default)]
    cdot_version: Option<String>,

    /// Gene information (we ignore this for now)
    #[serde(default)]
    genes: Option<serde_json::Value>,

    /// Metadata (we ignore this for now)
    #[serde(default)]
    metadata: Option<serde_json::Value>,
}

/// cdot JSON file structure (normalized for internal use).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CdotFile {
    /// Transcripts indexed by accession.
    pub transcripts: HashMap<String, CdotTranscript>,

    /// Genome builds information (optional).
    #[serde(default)]
    pub genome_builds: Option<HashMap<String, serde_json::Value>>,
}

/// Coordinate mapper using cdot data.
#[derive(Debug, Clone)]
pub struct CdotMapper {
    /// Transcripts indexed by accession.
    transcripts: HashMap<String, CdotTranscript>,
    /// Index from contig to transcript IDs that overlap.
    contig_index: HashMap<String, Vec<String>>,
    /// Index from base accession (without version) to versioned accession.
    base_to_versioned: HashMap<String, String>,
    /// LRG transcript to RefSeq transcript mapping (e.g., "LRG_1t1" -> "NM_000088.3").
    lrg_to_refseq: HashMap<String, String>,
}

impl CdotMapper {
    /// Create a new empty mapper.
    pub fn new() -> Self {
        Self {
            transcripts: HashMap::new(),
            contig_index: HashMap::new(),
            base_to_versioned: HashMap::new(),
            lrg_to_refseq: HashMap::new(),
        }
    }

    /// Load from a cdot JSON file.
    pub fn from_json_file<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open cdot file: {}", e),
        })?;
        let reader = BufReader::new(file);
        Self::from_reader(reader)
    }

    /// Load from a gzip-compressed cdot JSON file.
    pub fn from_json_gz<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open cdot file: {}", e),
        })?;
        let reader = BufReader::new(file);
        let decoder = flate2::read::GzDecoder::new(reader);
        let reader = BufReader::new(decoder);
        Self::from_reader(reader)
    }

    /// Load from any reader, defaulting to GRCh38 genome build.
    pub fn from_reader<R: Read>(reader: R) -> Result<Self, FerroError> {
        Self::from_reader_with_build(reader, "GRCh38")
    }

    /// Load from any reader with a specific genome build.
    pub fn from_reader_with_build<R: Read>(
        reader: R,
        genome_build: &str,
    ) -> Result<Self, FerroError> {
        let raw_file: RawCdotFile = serde_json::from_reader(reader)?;
        Ok(Self::from_raw_cdot_file(raw_file, genome_build))
    }

    /// Create from a raw CdotFile, converting to internal format.
    fn from_raw_cdot_file(raw_file: RawCdotFile, genome_build: &str) -> Self {
        let mut mapper = Self::new();

        for (accession, raw_transcript) in raw_file.transcripts {
            if let Some(transcript) = raw_transcript.to_transcript(genome_build) {
                mapper.add_transcript(accession, transcript);
            }
        }

        mapper
    }

    /// Create from a parsed CdotFile (already normalized).
    pub fn from_cdot_file(cdot_file: CdotFile) -> Self {
        let mut mapper = Self::new();

        for (accession, transcript) in cdot_file.transcripts {
            mapper.add_transcript(accession, transcript);
        }

        mapper
    }

    /// Add a transcript to the mapper.
    pub fn add_transcript(&mut self, accession: String, transcript: CdotTranscript) {
        // Update contig index
        self.contig_index
            .entry(transcript.contig.clone())
            .or_default()
            .push(accession.clone());

        // Update base to versioned index (e.g., "NM_000088" -> "NM_000088.4")
        if let Some(base) = accession.split('.').next() {
            self.base_to_versioned
                .insert(base.to_string(), accession.clone());
        }

        self.transcripts.insert(accession, transcript);
    }

    /// Load LRG to RefSeq transcript mapping from a file.
    ///
    /// The file format is tab-separated with columns:
    /// LRG, HGNC_SYMBOL, REFSEQ_GENOMIC, LRG_TRANSCRIPT, REFSEQ_TRANSCRIPT, ENSEMBL_TRANSCRIPT, CCDS
    ///
    /// This creates mappings like "LRG_1t1" -> "NM_000088.3".
    pub fn load_lrg_mapping<P: AsRef<Path>>(&mut self, path: P) -> Result<usize, FerroError> {
        use std::io::BufRead;

        let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open LRG mapping file: {}", e),
        })?;
        let reader = std::io::BufReader::new(file);

        let mut count = 0;
        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read LRG mapping line: {}", e),
            })?;

            // Skip comments and empty lines
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                continue;
            }

            // Fields: LRG, HGNC_SYMBOL, REFSEQ_GENOMIC, LRG_TRANSCRIPT, REFSEQ_TRANSCRIPT
            let lrg_id = fields[0]; // e.g., "LRG_1"
            let lrg_transcript = fields[3]; // e.g., "t1"
            let refseq_transcript = fields[4]; // e.g., "NM_000088.3"

            // Skip if RefSeq transcript is empty or "-"
            if refseq_transcript.is_empty() || refseq_transcript == "-" {
                continue;
            }

            // Create full LRG transcript ID: "LRG_1t1"
            let lrg_full = format!("{}{}", lrg_id, lrg_transcript);
            self.lrg_to_refseq
                .insert(lrg_full, refseq_transcript.to_string());
            count += 1;
        }

        Ok(count)
    }

    /// Get the RefSeq transcript ID for an LRG transcript.
    pub fn get_refseq_for_lrg(&self, lrg_accession: &str) -> Option<&str> {
        self.lrg_to_refseq.get(lrg_accession).map(|s| s.as_str())
    }

    /// Get a transcript by accession.
    ///
    /// For LRG transcripts (e.g., "LRG_1t1"), this will look up the equivalent
    /// RefSeq transcript if an LRG mapping has been loaded.
    ///
    /// Supports version fallback: if "NM_000088.2" is not found but "NM_000088.3"
    /// exists, the latter will be returned.
    pub fn get_transcript(&self, accession: &str) -> Option<&CdotTranscript> {
        // First try direct lookup
        if let Some(tx) = self.transcripts.get(accession) {
            return Some(tx);
        }

        // Try version fallback (e.g., NM_000088.2 -> NM_000088.3)
        if let Some(base) = accession.split('.').next() {
            if let Some(versioned) = self.base_to_versioned.get(base) {
                if let Some(tx) = self.transcripts.get(versioned) {
                    return Some(tx);
                }
            }
        }

        // For LRG transcripts, try to find via RefSeq mapping
        if accession.starts_with("LRG_") {
            if let Some(refseq) = self.lrg_to_refseq.get(accession) {
                // Try direct lookup of RefSeq
                if let Some(tx) = self.transcripts.get(refseq) {
                    return Some(tx);
                }
                // Try version fallback for RefSeq too
                if let Some(base) = refseq.split('.').next() {
                    if let Some(versioned) = self.base_to_versioned.get(base) {
                        if let Some(tx) = self.transcripts.get(versioned) {
                            return Some(tx);
                        }
                    }
                }
            }
        }

        None
    }

    /// Get all transcripts on a contig.
    pub fn transcripts_on_contig(&self, contig: &str) -> Vec<&str> {
        self.contig_index
            .get(contig)
            .map(|ids| ids.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default()
    }

    /// Find transcripts overlapping a genomic position.
    pub fn transcripts_at_position(&self, contig: &str, pos: u64) -> Vec<(&str, &CdotTranscript)> {
        self.contig_index
            .get(contig)
            .map(|ids| {
                ids.iter()
                    .filter_map(|id| {
                        let tx = self.transcripts.get(id)?;
                        // Check if position is within transcript genomic range
                        let (min, max) = tx.exons.iter().fold((u64::MAX, 0), |(min, max), e| {
                            (min.min(e[0]), max.max(e[1]))
                        });
                        if pos >= min && pos < max {
                            Some((id.as_str(), tx))
                        } else {
                            None
                        }
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get the number of transcripts loaded.
    pub fn transcript_count(&self) -> usize {
        self.transcripts.len()
    }

    /// Get all transcript accessions.
    pub fn transcript_ids(&self) -> impl Iterator<Item = &str> {
        self.transcripts.keys().map(|s| s.as_str())
    }
}

impl Default for CdotMapper {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_transcript() -> CdotTranscript {
        // Simple transcript with 3 exons on + strand
        // Exons: [genome_start, genome_end, tx_start, tx_end]
        CdotTranscript {
            gene_name: Some("TEST".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![
                [1000, 1100, 0, 100],   // Exon 1: 100bp
                [2000, 2200, 100, 300], // Exon 2: 200bp
                [3000, 3150, 300, 450], // Exon 3: 150bp
            ],
            cds_start: Some(50), // CDS starts 50bp into transcript
            cds_end: Some(400),  // CDS ends 400bp into transcript
            exon_cigars: Vec::new(),
            gene_id: None,
            protein: None,
        }
    }

    fn minus_strand_transcript() -> CdotTranscript {
        // Transcript on - strand
        CdotTranscript {
            gene_name: Some("MINUS".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Minus,
            exons: vec![
                [3000, 3150, 0, 150],   // Exon 1 (most 3' on genome)
                [2000, 2200, 150, 350], // Exon 2
                [1000, 1100, 350, 450], // Exon 3 (most 5' on genome)
            ],
            cds_start: Some(50),
            cds_end: Some(400),
            exon_cigars: Vec::new(),
            gene_id: None,
            protein: None,
        }
    }

    // CdotTranscript tests
    #[test]
    fn test_transcript_length() {
        let tx = sample_transcript();
        assert_eq!(tx.transcript_length(), 450);
    }

    #[test]
    fn test_cds_length() {
        let tx = sample_transcript();
        assert_eq!(tx.cds_length(), Some(350));
    }

    #[test]
    fn test_exon_for_tx_pos() {
        let tx = sample_transcript();

        // Position in first exon
        let exon = tx.exon_for_tx_pos(50).unwrap();
        assert_eq!(exon.number, 1);

        // Position in second exon
        let exon = tx.exon_for_tx_pos(150).unwrap();
        assert_eq!(exon.number, 2);

        // Position in third exon
        let exon = tx.exon_for_tx_pos(350).unwrap();
        assert_eq!(exon.number, 3);

        // Position beyond transcript
        assert!(tx.exon_for_tx_pos(500).is_none());
    }

    #[test]
    fn test_exon_for_genome_pos() {
        let tx = sample_transcript();

        // Position in first exon
        let exon = tx.exon_for_genome_pos(1050).unwrap();
        assert_eq!(exon.number, 1);

        // Position in second exon
        let exon = tx.exon_for_genome_pos(2100).unwrap();
        assert_eq!(exon.number, 2);

        // Position in intron
        assert!(tx.exon_for_genome_pos(1500).is_none());
    }

    #[test]
    fn test_tx_to_genome_plus_strand() {
        let tx = sample_transcript();

        // First exon
        assert_eq!(tx.tx_to_genome(0), Some(1000));
        assert_eq!(tx.tx_to_genome(50), Some(1050));
        assert_eq!(tx.tx_to_genome(99), Some(1099));

        // Second exon
        assert_eq!(tx.tx_to_genome(100), Some(2000));
        assert_eq!(tx.tx_to_genome(200), Some(2100));

        // Third exon
        assert_eq!(tx.tx_to_genome(300), Some(3000));
    }

    #[test]
    fn test_tx_to_genome_minus_strand() {
        let tx = minus_strand_transcript();

        // First exon (3' on genome)
        assert_eq!(tx.tx_to_genome(0), Some(3149));
        assert_eq!(tx.tx_to_genome(50), Some(3099));
        assert_eq!(tx.tx_to_genome(149), Some(3000));

        // Second exon
        assert_eq!(tx.tx_to_genome(150), Some(2199));
    }

    #[test]
    fn test_genome_to_tx_plus_strand() {
        let tx = sample_transcript();

        // First exon
        assert_eq!(tx.genome_to_tx(1000), Some(0));
        assert_eq!(tx.genome_to_tx(1050), Some(50));

        // Second exon
        assert_eq!(tx.genome_to_tx(2000), Some(100));

        // Intron
        assert!(tx.genome_to_tx(1500).is_none());
    }

    #[test]
    fn test_genome_to_tx_minus_strand() {
        let tx = minus_strand_transcript();

        // First exon (3' on genome)
        assert_eq!(tx.genome_to_tx(3149), Some(0));
        assert_eq!(tx.genome_to_tx(3000), Some(149));

        // Second exon
        assert_eq!(tx.genome_to_tx(2199), Some(150));
    }

    #[test]
    fn test_cds_to_tx() {
        let tx = sample_transcript();

        // Normal CDS position
        assert_eq!(tx.cds_to_tx(1), Some(50)); // c.1 -> tx pos 50
        assert_eq!(tx.cds_to_tx(100), Some(149)); // c.100 -> tx pos 149

        // 5' UTR
        assert_eq!(tx.cds_to_tx(-1), Some(49)); // c.-1 -> tx pos 49
        assert_eq!(tx.cds_to_tx(-50), Some(0)); // c.-50 -> tx pos 0

        // Position 0 invalid
        assert!(tx.cds_to_tx(0).is_none());
    }

    #[test]
    fn test_tx_to_cds() {
        let tx = sample_transcript();

        // 5' UTR
        assert_eq!(tx.tx_to_cds(0), Some(CdsPosition::FivePrimeUtr(50)));
        assert_eq!(tx.tx_to_cds(49), Some(CdsPosition::FivePrimeUtr(1)));

        // CDS
        assert_eq!(tx.tx_to_cds(50), Some(CdsPosition::Cds(1)));
        assert_eq!(tx.tx_to_cds(149), Some(CdsPosition::Cds(100)));
        assert_eq!(tx.tx_to_cds(399), Some(CdsPosition::Cds(350)));

        // 3' UTR
        assert_eq!(tx.tx_to_cds(400), Some(CdsPosition::ThreePrimeUtr(1)));
        assert_eq!(tx.tx_to_cds(449), Some(CdsPosition::ThreePrimeUtr(50)));
    }

    // Exon tests
    #[test]
    fn test_exon_len() {
        let exon = Exon::from_array(1, [1000, 1100, 0, 100]);
        assert_eq!(exon.len(), 100);
        assert!(!exon.is_empty());
    }

    #[test]
    fn test_exon_contains() {
        let exon = Exon::from_array(1, [1000, 1100, 0, 100]);

        assert!(exon.contains_genome_pos(1000));
        assert!(exon.contains_genome_pos(1099));
        assert!(!exon.contains_genome_pos(1100));

        assert!(exon.contains_tx_pos(0));
        assert!(exon.contains_tx_pos(99));
        assert!(!exon.contains_tx_pos(100));
    }

    // CdotMapper tests
    #[test]
    fn test_mapper_new() {
        let mapper = CdotMapper::new();
        assert_eq!(mapper.transcript_count(), 0);
    }

    #[test]
    fn test_mapper_add_transcript() {
        let mut mapper = CdotMapper::new();
        mapper.add_transcript("NM_000088.3".to_string(), sample_transcript());

        assert_eq!(mapper.transcript_count(), 1);
        assert!(mapper.get_transcript("NM_000088.3").is_some());
        // Version fallback: NM_000088.4 falls back to NM_000088.3
        assert!(mapper.get_transcript("NM_000088.4").is_some());
        // Completely different accession should return None
        assert!(mapper.get_transcript("NM_999999.1").is_none());
    }

    #[test]
    fn test_mapper_contig_index() {
        let mut mapper = CdotMapper::new();
        mapper.add_transcript("NM_000088.3".to_string(), sample_transcript());
        mapper.add_transcript("NM_000088.4".to_string(), sample_transcript());

        let tx_ids = mapper.transcripts_on_contig("NC_000001.11");
        assert_eq!(tx_ids.len(), 2);
        assert!(tx_ids.contains(&"NM_000088.3"));
        assert!(tx_ids.contains(&"NM_000088.4"));

        let tx_ids = mapper.transcripts_on_contig("NC_000002.12");
        assert!(tx_ids.is_empty());
    }

    #[test]
    fn test_mapper_transcripts_at_position() {
        let mut mapper = CdotMapper::new();
        mapper.add_transcript("NM_000088.3".to_string(), sample_transcript());

        // Position within transcript
        let results = mapper.transcripts_at_position("NC_000001.11", 2050);
        assert_eq!(results.len(), 1);

        // Position outside transcript
        let results = mapper.transcripts_at_position("NC_000001.11", 5000);
        assert!(results.is_empty());
    }

    #[test]
    fn test_mapper_from_json() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [
                        [50184096, 50184169, 0, 73],
                        [50185022, 50185148, 73, 199]
                    ],
                    "cds_start": 149,
                    "cds_end": 4544
                }
            }
        }
        "#;

        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();
        assert_eq!(mapper.transcript_count(), 1);

        let tx = mapper.get_transcript("NM_000088.3").unwrap();
        assert_eq!(tx.gene_name.as_deref(), Some("COL1A1"));
        assert_eq!(tx.strand, Strand::Plus);
        assert_eq!(tx.exons.len(), 2);
    }

    // =========================================================================
    // CIGAR parsing tests
    // =========================================================================

    #[test]
    fn test_parse_cigar_simple_match() {
        let ops = parse_cigar("M185").unwrap();
        assert_eq!(ops, vec![CigarOp::Match(185)]);
    }

    #[test]
    fn test_parse_cigar_with_insertion() {
        let ops = parse_cigar("M185 I3 M250").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::Match(185),
                CigarOp::Insertion(3),
                CigarOp::Match(250),
            ]
        );
    }

    #[test]
    fn test_parse_cigar_with_deletion() {
        let ops = parse_cigar("M504 D2 M123").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::Match(504),
                CigarOp::Deletion(2),
                CigarOp::Match(123),
            ]
        );
    }

    #[test]
    fn test_parse_cigar_complex() {
        let ops = parse_cigar("M6 D1 M4 I2 M3").unwrap();
        assert_eq!(
            ops,
            vec![
                CigarOp::Match(6),
                CigarOp::Deletion(1),
                CigarOp::Match(4),
                CigarOp::Insertion(2),
                CigarOp::Match(3),
            ]
        );
    }

    #[test]
    fn test_parse_cigar_empty() {
        let ops = parse_cigar("").unwrap();
        assert!(ops.is_empty());
    }

    #[test]
    fn test_parse_cigar_whitespace_only() {
        let ops = parse_cigar("   ").unwrap();
        assert!(ops.is_empty());
    }

    #[test]
    fn test_parse_cigar_invalid_operation() {
        let result = parse_cigar("X5");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_cigar_invalid_length() {
        let result = parse_cigar("Mabc");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_cigar_too_short_token() {
        let result = parse_cigar("M");
        assert!(result.is_err());
    }

    #[test]
    fn test_cumulative_insertion_offset_before_insertion() {
        // CIGAR: M185 I3 M250 — 3bp insertion at tx position 186-188
        // At tx position 100 (before insertion): offset = 0
        let cigar = vec![
            CigarOp::Match(185),
            CigarOp::Insertion(3),
            CigarOp::Match(250),
        ];
        assert_eq!(cumulative_insertion_offset(&cigar, 100), 0);
    }

    #[test]
    fn test_cumulative_insertion_offset_after_insertion() {
        // At tx position 200 (after insertion): offset = 3
        let cigar = vec![
            CigarOp::Match(185),
            CigarOp::Insertion(3),
            CigarOp::Match(250),
        ];
        assert_eq!(cumulative_insertion_offset(&cigar, 200), 3);
    }

    #[test]
    fn test_cumulative_insertion_offset_at_end() {
        // At tx position 437 (end of exon, 185+3+250-1): offset = 3
        let cigar = vec![
            CigarOp::Match(185),
            CigarOp::Insertion(3),
            CigarOp::Match(250),
        ];
        assert_eq!(cumulative_insertion_offset(&cigar, 437), 3);
    }

    #[test]
    fn test_cumulative_insertion_offset_no_insertions() {
        let cigar = vec![CigarOp::Match(500)];
        assert_eq!(cumulative_insertion_offset(&cigar, 250), 0);
    }

    #[test]
    fn test_cumulative_insertion_offset_multiple_insertions() {
        // M100 I2 M100 I5 M100 — two insertions
        let cigar = vec![
            CigarOp::Match(100),
            CigarOp::Insertion(2),
            CigarOp::Match(100),
            CigarOp::Insertion(5),
            CigarOp::Match(100),
        ];
        // Before first insertion
        assert_eq!(cumulative_insertion_offset(&cigar, 50), 0);
        // After first insertion, before second
        assert_eq!(cumulative_insertion_offset(&cigar, 150), 2);
        // After both insertions
        assert_eq!(cumulative_insertion_offset(&cigar, 300), 7);
    }

    #[test]
    fn test_cumulative_insertion_offset_with_deletion() {
        // M100 D5 M100 I3 M100 — deletion doesn't affect insertion offset
        let cigar = vec![
            CigarOp::Match(100),
            CigarOp::Deletion(5),
            CigarOp::Match(100),
            CigarOp::Insertion(3),
            CigarOp::Match(100),
        ];
        // Before insertion (deletions don't advance tx position)
        assert_eq!(cumulative_insertion_offset(&cigar, 150), 0);
        // After insertion
        assert_eq!(cumulative_insertion_offset(&cigar, 250), 3);
    }
}
