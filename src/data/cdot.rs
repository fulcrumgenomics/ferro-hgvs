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
//! | Genomic coordinates (`genome_start`, `genome_end`) | HGVS-value-based | Half-open `[start, end)` where `start` equals the numeric HGVS g. position of the first exonic base |
//! | Transcript coordinates (`tx_start`, `tx_end`) | 0-based | Half-open `[start, end)` |
//! | CDS coordinates (`cds_start`, `cds_end`) | 0-based | Half-open `[start, end)` |
//!
//! **Note on genomic coordinates**: `genome_start` stores the numeric value of the first HGVS
//! g. position in the exon (e.g., an exon starting at g.1000 has `genome_start = 1000`).
//! `genome_end` is one past the last HGVS g. position (exclusive). This differs from a
//! strict 0-based system: an exon at g.1000–g.1008 has `genome_start = 1000, genome_end = 1009`.
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::error::FerroError;
use crate::liftover::aliases::ContigAliases;
use crate::reference::transcript::GenomeBuild;
use crate::reference::Strand;
use bincode::Options;
use once_cell::sync::OnceCell;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;
use superintervals::IntervalMap;

/// Magic bytes prefixing a cdot bincode cache file, followed by a `u32`
/// schema version. The version MUST be bumped whenever the serialized layout
/// of [`CdotMapperSnapshot`] changes (a field added/removed/reordered, or a
/// nested type's layout changes), so a cache written by an older build is
/// detected as stale and refreshed rather than silently mis-deserialized.
///
/// Without this guard a layout change makes `bincode` read past the end of the
/// file (`failed to fill whole buffer`) only *after* a partial read, and the
/// loader falls back to parsing the multi-hundred-MB JSON on every run — ~10x
/// slower — with no indication anything is wrong.
const CDOT_BINCODE_MAGIC: [u8; 4] = *b"FCDT";
const CDOT_BINCODE_VERSION: u32 = 1;

/// cdot transcript alignment data (normalized internal representation).
///
/// # Coordinate Systems
/// - Genomic: HGVS-value-based half-open `[start, end)` (see module-level docs)
/// - Transcript (tx_start/tx_end in exons): 0-based half-open `[start, end)`
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

    /// Exon alignments: [genome_start(HGVS-val), genome_end(HGVS-val excl), tx_start(0-based), tx_end(0-based excl)].
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
    /// Protein accession for the transcript's CDS (e.g. `NP_000068.1`). cdot
    /// records this per-transcript; plumbing it through lets the projector
    /// emit the correct `p.` accession instead of falling back to the
    /// (frequently wrong) `NM_*` → `NP_*` number-preserving inference.
    #[serde(default)]
    protein: Option<String>,
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
            protein: self.protein.clone(),
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
            protein: self.protein.clone(),
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
                Strand::Unknown => return None,
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
/// - Genomic coordinates: HGVS-value-based half-open `[genome_start, genome_end)` where
///   `genome_start` equals the numeric value of the first HGVS g. position in the exon
///   and `genome_end` is one past the last (exclusive). For example, an exon spanning
///   g.1000–g.1008 has `genome_start = 1000`, `genome_end = 1009`.
/// - Transcript coordinates: 0-based half-open `[tx_start, tx_end)` where `tx_start = 0`
///   for the first base of the transcript. `contains_tx_pos` uses `pos >= tx_start && pos < tx_end`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Exon {
    /// Exon number (1-based).
    pub number: u32,
    /// Start position in genomic coordinates (HGVS-value-based, inclusive).
    /// Equals the numeric value of the first HGVS g. position in the exon.
    pub genome_start: u64,
    /// End position in genomic coordinates (HGVS-value-based, exclusive).
    /// Equals genome_start + exon_length.
    pub genome_end: u64,
    /// Start position in transcript coordinates (0-based, inclusive).
    pub tx_start: u64,
    /// End position in transcript coordinates (0-based, exclusive).
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
    ///
    /// Exons in `self.exons` are sorted by `tx_start` (see
    /// `RawCdotTranscript::from_genome_build`). On a plus-strand transcript
    /// that means `genome_start` is monotonically increasing; on a
    /// minus-strand transcript it's monotonically *decreasing* (the
    /// 5'→3' transcript order runs against genomic coordinates). Use that
    /// to early-terminate the scan once we've passed `genome_pos` — for
    /// intronic / off-transcript positions, the previous unconditional
    /// linear walk consulted every exon even though the answer was already
    /// determined a few iterations in.
    pub fn exon_for_genome_pos(&self, genome_pos: u64) -> Option<Exon> {
        for (i, arr) in self.exons.iter().enumerate() {
            if genome_pos >= arr[0] && genome_pos < arr[1] {
                return Some(Exon::from_array((i + 1) as u32, *arr));
            }
            // Early-terminate based on the sort invariant.
            match self.strand {
                Strand::Plus if genome_pos < arr[0] => return None,
                Strand::Minus if genome_pos >= arr[1] => return None,
                _ => {}
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
            Strand::Unknown => None,
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
            Strand::Unknown => None,
        }
    }

    /// Locate `genome_pos` against this transcript's exons in a single scan,
    /// returning the containing exon and the corresponding transcript-space
    /// position.
    ///
    /// Existing callers chained `genome_to_tx` (which internally scans the
    /// exon list) with a second `exon_for_genome_pos` call (another scan).
    /// This helper does both in one pass — relevant in
    /// [`CoordinateMapper::genome_pos_to_cds_pos`] which is invoked twice per
    /// projection and is the dominant per-variant cost on the SNP fixture
    /// after c11.
    pub fn locate_genome_pos(&self, genome_pos: u64) -> Option<(u64, Exon)> {
        for (i, arr) in self.exons.iter().enumerate() {
            if genome_pos >= arr[0] && genome_pos < arr[1] {
                let exon = Exon::from_array((i + 1) as u32, *arr);
                let tx_pos = match self.strand {
                    Strand::Plus => exon.tx_start + (genome_pos - exon.genome_start),
                    Strand::Minus => exon.tx_start + (exon.genome_end - 1 - genome_pos),
                    Strand::Unknown => return None,
                };
                return Some((tx_pos, exon));
            }
            // Same early-termination as in `exon_for_genome_pos`. Exons are
            // sorted by `tx_start`, which is monotone in `genome_start` on
            // plus-strand (increasing) and in `genome_end` on minus-strand
            // (decreasing).
            match self.strand {
                Strand::Plus if genome_pos < arr[0] => return None,
                Strand::Minus if genome_pos >= arr[1] => return None,
                _ => {}
            }
        }
        None
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
                Strand::Unknown => return None,
            };

            if genome_pos >= intron_start && genome_pos < intron_end {
                let intron_num = (i + 1) as u32;
                let offset = match self.strand {
                    Strand::Plus => (genome_pos as i64) - (current.genome_end as i64),
                    Strand::Minus => (current.genome_start as i64) - (genome_pos as i64),
                    Strand::Unknown => return None,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
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

/// Bincode-friendly snapshot of CdotMapper data, without custom JSON deserializers.
/// Used for deserialization only; serialization uses [`CdotMapperSnapshotRef`].
///
/// **Format note (#389 follow-up).** `primary_build` was added after the
/// initial snapshot layout. bincode (1.x) is positional and does not honor
/// `#[serde(default)]` for trailing fields, so snapshots produced by older
/// builds will fail to deserialize against this struct. That failure is
/// non-fatal: [`CdotMapper::load`] catches it and falls back to JSON, which
/// regenerates the `.bin` with the new layout. Callers that pass a `.bin`
/// path directly (no JSON fallback available) must regenerate via
/// `ferro prepare`.
#[derive(Deserialize)]
struct CdotMapperSnapshot {
    transcripts: HashMap<String, CdotTranscriptSnapshot>,
    contig_index: HashMap<String, Vec<String>>,
    contig_alias_to_canonical: HashMap<String, String>,
    base_to_versioned: HashMap<String, String>,
    lrg_to_refseq: HashMap<String, String>,
    /// Genome build the snapshot was produced under, or `None` if the
    /// snapshot pre-dates the build-tracking format (only reachable via a
    /// hand-rolled deserializer skipping this field — bincode itself will
    /// error on a truncated stream and the loader falls back to JSON).
    primary_build: Option<String>,
    /// Per-build transcript maps for non-primary builds (e.g. GRCh37 when the
    /// primary is GRCh38). Required for NG/NC-parent-aware resolution (#332);
    /// omitting it from the cache silently degraded that path on every fast
    /// load. Part of the v1 layout — a cache missing it fails the version/EOF
    /// guard and is regenerated.
    alt_build_transcripts: HashMap<String, HashMap<String, CdotTranscriptSnapshot>>,
}

/// Borrowed view of CdotMapper for zero-copy serialization to bincode.
#[derive(Serialize)]
struct CdotMapperSnapshotRef<'a> {
    transcripts: HashMap<&'a String, CdotTranscriptSnapshotRef<'a>>,
    contig_index: &'a HashMap<String, Vec<String>>,
    contig_alias_to_canonical: &'a HashMap<String, String>,
    base_to_versioned: &'a HashMap<String, String>,
    lrg_to_refseq: &'a HashMap<String, String>,
    primary_build: Option<&'a str>,
    alt_build_transcripts: HashMap<&'a String, HashMap<&'a String, CdotTranscriptSnapshotRef<'a>>>,
}

/// Bincode-friendly snapshot of CdotTranscript, using standard Strand serde.
/// Used for deserialization only; serialization uses [`CdotTranscriptSnapshotRef`].
#[derive(Deserialize)]
struct CdotTranscriptSnapshot {
    gene_name: Option<String>,
    contig: String,
    strand: Strand,
    exons: Vec<[u64; 4]>,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    exon_cigars: Vec<Option<Vec<CigarOp>>>,
    gene_id: Option<String>,
    protein: Option<String>,
}

/// Borrowed view of CdotTranscript for zero-copy serialization to bincode.
#[derive(Serialize)]
struct CdotTranscriptSnapshotRef<'a> {
    gene_name: &'a Option<String>,
    contig: &'a str,
    strand: Strand,
    exons: &'a Vec<[u64; 4]>,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    exon_cigars: &'a Vec<Option<Vec<CigarOp>>>,
    gene_id: &'a Option<String>,
    protein: &'a Option<String>,
}

impl<'a> From<&'a CdotTranscript> for CdotTranscriptSnapshotRef<'a> {
    fn from(tx: &'a CdotTranscript) -> Self {
        Self {
            gene_name: &tx.gene_name,
            contig: &tx.contig,
            strand: tx.strand,
            exons: &tx.exons,
            cds_start: tx.cds_start,
            cds_end: tx.cds_end,
            exon_cigars: &tx.exon_cigars,
            gene_id: &tx.gene_id,
            protein: &tx.protein,
        }
    }
}

impl From<CdotTranscriptSnapshot> for CdotTranscript {
    fn from(snap: CdotTranscriptSnapshot) -> Self {
        Self {
            gene_name: snap.gene_name,
            contig: snap.contig,
            strand: snap.strand,
            exons: snap.exons,
            cds_start: snap.cds_start,
            cds_end: snap.cds_end,
            exon_cigars: snap.exon_cigars,
            gene_id: snap.gene_id,
            protein: snap.protein,
        }
    }
}

/// Coordinate mapper using cdot data.
#[derive(Debug, Clone)]
pub struct CdotMapper {
    /// Transcripts indexed by accession. Entries are populated from whatever
    /// genome build was passed at load time (the "primary build"). Use
    /// [`get_transcript_on_build`](Self::get_transcript_on_build) to fetch a
    /// non-primary build's view of a transcript.
    transcripts: HashMap<String, CdotTranscript>,
    /// Index from contig to transcript IDs that overlap.
    contig_index: HashMap<String, Vec<String>>,
    /// Alias-to-canonical contig name mapping (e.g., "chr7" -> "NC_000007.14").
    /// Allows lookups by UCSC-style names when `contig_index` keys are RefSeq accessions.
    contig_alias_to_canonical: HashMap<String, String>,
    /// Index from base accession (without version) to versioned accession.
    base_to_versioned: HashMap<String, String>,
    /// LRG transcript to RefSeq transcript mapping (e.g., "LRG_1t1" -> "NM_000088.3").
    lrg_to_refseq: HashMap<String, String>,
    /// Per-transcript data on builds OTHER than the primary load build.
    /// Outer key is the genome-build name (e.g. "GRCh37"); inner key is the
    /// transcript accession. Populated when the source cdot JSON nests
    /// per-build data under `genome_builds`. Empty for mappers built without
    /// a multi-build source (bincode snapshots, `from_transcripts`, manual
    /// `add_transcript`).
    ///
    /// The primary build's data lives in [`Self::transcripts`] only; it is
    /// NOT duplicated here. [`get_transcript_on_build`](Self::get_transcript_on_build)
    /// dispatches between the two maps.
    alt_build_transcripts: HashMap<String, HashMap<String, CdotTranscript>>,
    /// The genome-build name corresponding to entries in [`Self::transcripts`],
    /// or `None` if the build is unknown (e.g. a bincode snapshot produced
    /// by an older build that did not persist this field). Used by
    /// [`get_transcript_on_build`](Self::get_transcript_on_build) to
    /// short-circuit when the caller asks for the primary build, and by
    /// [`primary_build`](Self::primary_build) so downstream callers can
    /// honor uncertainty rather than assume a default (#389 follow-up).
    primary_build: Option<String>,
    /// Lazily-built per-contig SuperIntervals index for stabbing queries.
    ///
    /// Built on the first call to `transcripts_at_position` from the contents of
    /// `contig_index` + `transcripts`. Replaces a linear scan that re-folded
    /// every transcript's exon table on every query. The `u32` payload is the
    /// index into `contig_index[contig]` so the original accession can be
    /// recovered for the caller.
    ///
    /// Built eagerly here because every "real" caller (`from_json_file`,
    /// `load`, `from_transcripts`) populates `contig_index` once and then only
    /// queries; `add_transcript` followed by a query in tests still works
    /// because the cell is initialised on first access. If a caller adds
    /// transcripts AFTER a query has populated the cell the new transcripts
    /// will be missing from the index — `add_transcript` clears the cell to
    /// guard against that.
    contig_query_index: OnceCell<HashMap<String, IntervalMap<u32>>>,
    /// Lazily-built per-transcript `(min_genome_start, max_genome_end)`
    /// cache so `VariantProjector::project_single_inner` doesn't re-fold the
    /// exon table on every projection. Sharing the same `OnceCell` build
    /// trigger as `contig_query_index` keeps the two views in sync.
    transcript_genome_spans: OnceCell<HashMap<String, (u64, u64)>>,
}

impl CdotMapper {
    /// Create a new empty mapper (primary build defaults to GRCh38).
    pub fn new() -> Self {
        Self {
            transcripts: HashMap::new(),
            contig_index: HashMap::new(),
            contig_alias_to_canonical: HashMap::new(),
            base_to_versioned: HashMap::new(),
            lrg_to_refseq: HashMap::new(),
            alt_build_transcripts: HashMap::new(),
            primary_build: Some("GRCh38".to_string()),
            contig_query_index: OnceCell::new(),
            transcript_genome_spans: OnceCell::new(),
        }
    }

    /// Build a `CdotMapper` from in-memory `Transcript` records (typically loaded into a
    /// `MockProvider`). Converts the 1-based `Transcript` coordinates into the form expected
    /// by `CdotTranscript`.
    ///
    /// # Coordinate conversions
    ///
    /// `Transcript.Exon` fields are documented as 1-based inclusive (see
    /// `src/reference/transcript.rs`). `CdotTranscript` exon arrays use:
    /// - `exons[][0]` (genome_start): HGVS-value-based — the numeric value of the first
    ///   HGVS g. position in the exon. Conversion: `genome_start = Exon.genomic_start` (no change).
    /// - `exons[][1]` (genome_end): exclusive, equals `Exon.genomic_end + 1`.
    /// - `exons[][2]` (tx_start): 0-based. Conversion: `tx_start = Exon.start - 1`.
    /// - `exons[][3]` (tx_end): 0-based exclusive. Same numeric value as `Exon.end` (1-based
    ///   inclusive end equals 0-based exclusive end).
    /// - `cds_start`: 0-based. Conversion: `cds_start = Transcript.cds_start - 1`.
    /// - `cds_end`: 0-based exclusive. Same numeric value as `Transcript.cds_end`.
    ///
    /// Defaults to GRCh38 for contig alias resolution (so `chr17` ↔ `NC_000017.11`
    /// works after construction). Use [`from_transcripts_with_build`] to specify a
    /// different genome build.
    pub fn from_transcripts<'a, I>(transcripts: I) -> Self
    where
        I: IntoIterator<Item = &'a crate::reference::transcript::Transcript>,
    {
        Self::from_transcripts_with_build(transcripts, "GRCh38")
    }

    /// Like [`from_transcripts`], but uses `genome_build` (e.g. "GRCh37" or
    /// "GRCh38") for contig alias resolution.
    ///
    /// Transcripts with missing or inconsistent coordinates are skipped with a
    /// `log::warn!`, rather than silently coerced to bogus intervals. A transcript
    /// must have a chromosome and, for every exon, all of: `genomic_start`,
    /// `genomic_end >= genomic_start`, `start >= 1` (1-based), `end >= start`.
    pub fn from_transcripts_with_build<'a, I>(transcripts: I, genome_build: &str) -> Self
    where
        I: IntoIterator<Item = &'a crate::reference::transcript::Transcript>,
    {
        let mut mapper = Self::new();
        mapper.primary_build = Some(genome_build.to_string());
        for tx in transcripts {
            let Some(contig) = tx.chromosome.clone() else {
                log::warn!("Skipping transcript {}: missing chromosome", tx.id);
                continue;
            };
            let mut exons: Vec<[u64; 4]> = Vec::with_capacity(tx.exons.len());
            let mut valid = true;
            for e in &tx.exons {
                // Require both genomic bounds and a non-empty, well-ordered interval
                // on both the transcript and genomic axes; 1-based `start == 0` is
                // invalid input and we refuse to silently turn it into `tx_start = 0`.
                let (Some(g_start), Some(g_end_inclusive)) = (e.genomic_start, e.genomic_end)
                else {
                    valid = false;
                    break;
                };
                let Some(t_start) = e.start.checked_sub(1) else {
                    valid = false;
                    break;
                };
                if e.end < e.start || g_end_inclusive < g_start {
                    valid = false;
                    break;
                }
                // genome_end: 1-based inclusive → exclusive (add 1).
                // tx_end: 1-based inclusive end equals 0-based exclusive end numerically.
                exons.push([g_start, g_end_inclusive + 1, t_start, e.end]);
            }
            if !valid || exons.is_empty() {
                log::warn!(
                    "Skipping transcript {}: incomplete or inconsistent exon coordinates",
                    tx.id
                );
                continue;
            }
            let cdot_tx = CdotTranscript {
                gene_name: tx.gene_symbol.clone(),
                contig,
                strand: tx.strand,
                exons,
                // cds_start: 1-based tx → 0-based. Use `checked_sub` so a stray `Some(0)`
                // becomes `None` rather than wrapping to a giant `u64`.
                cds_start: tx.cds_start.and_then(|c| c.checked_sub(1)),
                // cds_end: 1-based tx inclusive → 0-based exclusive (same value)
                cds_end: tx.cds_end,
                gene_id: None,
                // Propagate the ferro reference's explicit protein accession so
                // c. → p. projection works for non-RefSeq transcript IDs without
                // needing a cdot JSON. See #310.
                protein: tx.protein_id.clone(),
                exon_cigars: Vec::new(),
            };
            mapper.add_transcript(tx.id.clone(), cdot_tx);
        }
        mapper.populate_contig_aliases(genome_build);
        mapper
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

    /// Load from a pre-serialized bincode file (much faster than JSON).
    pub fn from_bincode_file<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open cdot bincode file: {}", e),
        })?;
        let file_size = file.metadata().map(|m| m.len()).unwrap_or(0);
        let mut reader = BufReader::new(file);

        // Verify the magic + schema version before deserializing. A mismatch
        // (including a legacy headerless cache, or one written by a build with a
        // different snapshot layout) is reported cleanly so the caller can fall
        // back to JSON and refresh the cache, instead of mis-deserializing.
        let mut header = [0u8; 8];
        reader.read_exact(&mut header).map_err(|e| FerroError::Io {
            msg: format!("cdot bincode too short to contain a header: {}", e),
        })?;
        if header[..4] != CDOT_BINCODE_MAGIC {
            return Err(FerroError::Io {
                msg: "cdot bincode has no valid magic (legacy or non-cdot file)".to_string(),
            });
        }
        let version = u32::from_le_bytes([header[4], header[5], header[6], header[7]]);
        if version != CDOT_BINCODE_VERSION {
            return Err(FerroError::Io {
                msg: format!(
                    "cdot bincode schema version {} != expected {}",
                    version, CDOT_BINCODE_VERSION
                ),
            });
        }

        // Use a size limit to prevent OOM on corrupt files. Allow 20x the file size
        // (bincode is compact, but deserialized HashMap overhead can be significant).
        // Limit deserialization size to prevent OOM on corrupt files (20x file size or 1MB min).
        let size_limit = file_size.saturating_mul(20).max(1024 * 1024);
        let snapshot: CdotMapperSnapshot = bincode::options()
            .with_fixint_encoding()
            .allow_trailing_bytes()
            .with_limit(size_limit)
            .deserialize_from(reader)
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to deserialize cdot bincode: {}", e),
            })?;
        Ok(Self {
            transcripts: snapshot
                .transcripts
                .into_iter()
                .map(|(k, v)| (k, v.into()))
                .collect(),
            contig_index: snapshot.contig_index,
            contig_alias_to_canonical: snapshot.contig_alias_to_canonical,
            base_to_versioned: snapshot.base_to_versioned,
            lrg_to_refseq: snapshot.lrg_to_refseq,
            // Alt-build maps are now part of the serialized format (v1), so
            // NG/NC-parent-aware resolution (#332) works on the fast path too.
            alt_build_transcripts: snapshot
                .alt_build_transcripts
                .into_iter()
                .map(|(build, txs)| (build, txs.into_iter().map(|(k, v)| (k, v.into())).collect()))
                .collect(),
            // Honor the snapshot's recorded primary build; older snapshots
            // that pre-date this field surface as `None` rather than
            // claiming a fabricated GRCh38 (#389 follow-up).
            primary_build: snapshot.primary_build,
            contig_query_index: OnceCell::new(),
            transcript_genome_spans: OnceCell::new(),
        })
    }

    /// Cheap check: does `path` begin with the current cdot bincode magic and
    /// schema version? Lets callers (e.g. `ferro prepare`) decide whether a
    /// cache needs regenerating without paying for a full deserialize.
    pub fn bincode_is_current<P: AsRef<Path>>(path: P) -> bool {
        let Ok(mut file) = File::open(path) else {
            return false;
        };
        let mut header = [0u8; 8];
        if file.read_exact(&mut header).is_err() {
            return false;
        }
        header[..4] == CDOT_BINCODE_MAGIC
            && u32::from_le_bytes([header[4], header[5], header[6], header[7]])
                == CDOT_BINCODE_VERSION
    }

    /// Write the mapper to a bincode file for fast subsequent loading.
    pub fn to_bincode_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FerroError> {
        let snapshot = CdotMapperSnapshotRef {
            transcripts: self
                .transcripts
                .iter()
                .map(|(k, v)| (k, v.into()))
                .collect(),
            contig_index: &self.contig_index,
            contig_alias_to_canonical: &self.contig_alias_to_canonical,
            base_to_versioned: &self.base_to_versioned,
            lrg_to_refseq: &self.lrg_to_refseq,
            primary_build: self.primary_build.as_deref(),
            alt_build_transcripts: self
                .alt_build_transcripts
                .iter()
                .map(|(build, txs)| (build, txs.iter().map(|(k, v)| (k, v.into())).collect()))
                .collect(),
        };
        // Write to a unique temp sibling, then atomically rename, so a crash
        // mid-write or a concurrent reader never observes a half-written cache.
        // The temp name combines the pid (unique across live processes) with a
        // process-global counter (unique across concurrent writers within this
        // process), so two threads writing the same cache never collide.
        use std::sync::atomic::{AtomicU64, Ordering};
        static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);
        let path = path.as_ref();
        let mut tmp = path.to_path_buf();
        let tmp_name = format!(
            "{}.tmp.{}.{}",
            path.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "cdot.bin".to_string()),
            std::process::id(),
            TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        );
        tmp.set_file_name(tmp_name);

        let file = File::create(&tmp).map_err(|e| FerroError::Io {
            msg: format!("Failed to create cdot bincode temp file: {}", e),
        })?;
        let mut writer = std::io::BufWriter::new(file);
        let write_result = (|| -> std::io::Result<()> {
            writer.write_all(&CDOT_BINCODE_MAGIC)?;
            writer.write_all(&CDOT_BINCODE_VERSION.to_le_bytes())?;
            bincode::options()
                .with_fixint_encoding()
                .serialize_into(&mut writer, &snapshot)
                .map_err(std::io::Error::other)?;
            writer.flush()
        })();
        if let Err(e) = write_result {
            let _ = std::fs::remove_file(&tmp);
            return Err(FerroError::Io {
                msg: format!("Failed to serialize cdot bincode: {}", e),
            });
        }
        std::fs::rename(&tmp, path).map_err(|e| {
            let _ = std::fs::remove_file(&tmp);
            FerroError::Io {
                msg: format!("Failed to finalize cdot bincode file: {}", e),
            }
        })
    }

    /// Load a cdot file, preferring a sibling `.bin` (bincode) file if available.
    ///
    /// Given a path to a `.json` or `.json.gz` file, checks for a `.bin` file in the
    /// same directory. If the bincode file exists, loads from it (much faster). Otherwise
    /// falls back to JSON parsing. Also handles `.bin` paths directly.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let path = path.as_ref();
        let path_str = path.to_string_lossy();

        // Direct bincode path — no JSON fallback available
        if path_str.ends_with(".bin") {
            return Self::from_bincode_file(path);
        }

        // Compute the canonical sibling .bin path:
        // - foo.json.gz -> foo.bin (strip .gz then replace .json with .bin)
        // - foo.json    -> foo.bin
        let bin_path = if path_str.ends_with(".json.gz") {
            path.with_extension("").with_extension("bin")
        } else {
            path.with_extension("bin")
        };
        if bin_path.exists() {
            match Self::from_bincode_file(&bin_path) {
                Ok(mapper) => return Ok(mapper),
                Err(e) => {
                    // Loud (stderr, not just `log`) — the binary may not
                    // initialize a logger, and silently eating the ~10x JSON
                    // parse on every run is exactly the failure mode this guards.
                    eprintln!(
                        "warning: cdot bincode cache {} is unusable ({}); \
                         parsing JSON (much slower) and refreshing the cache",
                        bin_path.display(),
                        e
                    );
                }
            }
        }

        // Fall back to JSON.
        let mapper = if path_str.ends_with(".gz") {
            Self::from_json_gz(path)?
        } else {
            Self::from_json_file(path)?
        };

        // Self-heal: best-effort refresh of the sibling bincode cache so the
        // next load takes the fast path. Failures here (e.g. a read-only
        // reference mount) are non-fatal — we already have a usable mapper.
        if let Err(e) = mapper.to_bincode_file(&bin_path) {
            eprintln!(
                "note: could not refresh cdot bincode cache {}: {} (continuing)",
                bin_path.display(),
                e
            );
        }

        Ok(mapper)
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
    ///
    /// Captures both the primary build's data (into [`Self::transcripts`])
    /// and any non-primary builds present in the source `genome_builds`
    /// (into [`Self::alt_build_transcripts`]), so callers can later request
    /// a non-default build via
    /// [`get_transcript_on_build`](Self::get_transcript_on_build). The
    /// non-primary capture is best-effort: transcripts that omit a given
    /// build are simply absent from that build's map.
    fn from_raw_cdot_file(raw_file: RawCdotFile, genome_build: &str) -> Self {
        let mut mapper = Self::new();
        mapper.primary_build = Some(genome_build.to_string());

        for (accession, raw_transcript) in raw_file.transcripts {
            // Stash any non-primary build views first so a missing primary-
            // build entry does not skip the alt builds.
            if let Some(builds) = &raw_transcript.genome_builds {
                for build_name in builds.keys() {
                    if build_name == genome_build {
                        continue;
                    }
                    if let Some(tx) = raw_transcript.to_transcript(build_name) {
                        mapper
                            .alt_build_transcripts
                            .entry(build_name.clone())
                            .or_default()
                            .insert(accession.clone(), tx);
                    }
                }
            }

            if let Some(transcript) = raw_transcript.to_transcript(genome_build) {
                mapper.add_transcript(accession, transcript);
            }
        }

        mapper.populate_contig_aliases(genome_build);

        mapper
    }

    /// Create from a parsed CdotFile (already normalized), defaulting to GRCh38 aliases.
    pub fn from_cdot_file(cdot_file: CdotFile) -> Self {
        Self::from_cdot_file_with_build(cdot_file, "GRCh38")
    }

    /// Create from a parsed CdotFile with a specific genome build for contig aliases.
    pub fn from_cdot_file_with_build(cdot_file: CdotFile, genome_build: &str) -> Self {
        let mut mapper = Self::new();
        mapper.primary_build = Some(genome_build.to_string());

        for (accession, transcript) in cdot_file.transcripts {
            mapper.add_transcript(accession, transcript);
        }

        mapper.populate_contig_aliases(genome_build);

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

        // Any previously-built query index is stale now. The next query will
        // rebuild it from the updated `contig_index` + `transcripts`.
        self.contig_query_index = OnceCell::new();
        self.transcript_genome_spans = OnceCell::new();
    }

    /// Build contig alias mappings from `ContigAliases` so that UCSC-style names
    /// (e.g., "chr7") resolve to the RefSeq accessions used as `contig_index` keys.
    fn populate_contig_aliases(&mut self, genome_build: &str) {
        let build = match genome_build {
            "GRCh37" => GenomeBuild::GRCh37,
            "GRCh38" => GenomeBuild::GRCh38,
            _ => {
                log::warn!(
                    "Unsupported genome build '{}'; contig alias resolution is disabled",
                    genome_build
                );
                return;
            }
        };

        let aliases = ContigAliases::default_human();
        for refseq_contig in self.contig_index.keys() {
            // Map UCSC alias (e.g., "chr7") -> RefSeq key (e.g., "NC_000007.14")
            if let Some(ucsc) = aliases.refseq_to_ucsc(refseq_contig) {
                self.contig_alias_to_canonical
                    .insert(ucsc.to_string(), refseq_contig.clone());
            }
            // Map Ensembl alias (e.g., "7") -> RefSeq key
            if let Some(ensembl) = aliases.refseq_to_ensembl(refseq_contig) {
                self.contig_alias_to_canonical
                    .insert(ensembl.to_string(), refseq_contig.clone());
            }
            // Also allow lookup by the other build's RefSeq accession via resolve
            // (e.g., if someone passes "NC_000007.13" but data is GRCh38)
            // This is handled implicitly since we only index contigs present in the data.
        }

        // Also map any non-RefSeq contig_index keys back through to_refseq
        // in case the cdot data itself uses UCSC names.
        let ucsc_keys: Vec<String> = self
            .contig_index
            .keys()
            .filter(|k| k.starts_with("chr"))
            .cloned()
            .collect();
        for ucsc_key in ucsc_keys {
            if let Some(refseq) = aliases.resolve_to_refseq(&ucsc_key, build) {
                self.contig_alias_to_canonical
                    .insert(refseq.to_string(), ucsc_key.clone());
            }
        }
    }

    /// Resolve a contig name through aliases, returning the canonical key used in `contig_index`.
    fn resolve_contig<'a>(&'a self, contig: &'a str) -> &'a str {
        if self.contig_index.contains_key(contig) {
            contig
        } else {
            self.contig_alias_to_canonical
                .get(contig)
                .map(|s| s.as_str())
                .unwrap_or(contig)
        }
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
    /// Whether the mapper has a transcript at exactly the requested accession
    /// (no version-strip or LRG fallback). Used by callers that need to
    /// distinguish "we have your exact version" from "we resolved to a
    /// sibling version" — e.g. the cdot-driven base synthesis path in
    /// `MultiFastaProvider::get_transcript` (closes #331).
    pub fn has_transcript_exact(&self, accession: &str) -> bool {
        self.transcripts.contains_key(accession)
    }

    /// Name of the genome build whose data populates [`Self::transcripts`],
    /// or `None` if the build is unknown.
    ///
    /// This is the "primary" build set at load time
    /// (e.g. via [`from_transcripts_with_build`](Self::from_transcripts_with_build)
    /// or [`from_cdot_file_with_build`](Self::from_cdot_file_with_build));
    /// any other builds live in [`Self::alt_build_transcripts`] and are
    /// reachable only through [`get_transcript_on_build`](Self::get_transcript_on_build).
    ///
    /// Returns `None` when the mapper was loaded from a bincode snapshot
    /// produced before the build field was persisted — downstream callers
    /// must surface that uncertainty rather than substituting a default
    /// (#389 follow-up).
    pub fn primary_build(&self) -> Option<&str> {
        self.primary_build.as_deref()
    }

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

    /// Get a transcript by accession, restricted to a specific genome build.
    ///
    /// When `build` matches the mapper's primary load build, delegates to
    /// [`get_transcript`](Self::get_transcript). Otherwise consults the
    /// alt-build cache (populated from the source cdot JSON's `genome_builds`
    /// section). Returns `None` if the transcript is absent from the requested
    /// build (or if the alt-build cache is empty, e.g. for bincode-loaded
    /// mappers or `from_transcripts`-built mappers).
    ///
    /// This is the lookup used by NG/NC-parent-aware transcript resolution
    /// (issue #332): when the input is `NG_*(NM_*):c.…` or
    /// `NC_*(NM_*):c.…`, the caller infers a build hint from the parent and
    /// asks for that build's view of the NM, so the resulting transcript
    /// carries the correct chromosome alignment.
    pub fn get_transcript_on_build(&self, accession: &str, build: &str) -> Option<&CdotTranscript> {
        // Treat unknown primary build (bincode snapshots without the field)
        // as a match: `alt_build_transcripts` is empty for those mappers so
        // the only data available lives in `transcripts`, and probing that
        // map is strictly better than returning `None` (#389 follow-up).
        let primary_matches = match self.primary_build.as_deref() {
            Some(pb) => pb == build,
            None => true,
        };
        if primary_matches {
            return self.get_transcript(accession);
        }
        // Alt-build path: resolve LRG → RefSeq before consulting the alt
        // map so build-specific lookups honor the same LRG fallback as
        // [`get_transcript`] (#389 follow-up). Without this, a non-primary
        // LRG lookup would return `None` even when the mapped RefSeq
        // transcript exists in the alt build.
        let resolved: &str = if accession.starts_with("LRG_") {
            self.lrg_to_refseq
                .get(accession)
                .map(String::as_str)
                .unwrap_or(accession)
        } else {
            accession
        };
        let alt = self.alt_build_transcripts.get(build)?;
        if let Some(tx) = alt.get(resolved) {
            return Some(tx);
        }
        // Version fallback within the alt-build map. Deterministically pick
        // the highest numeric version among siblings sharing the base, rather
        // than the first HashMap iteration hit — alt-build maps occasionally
        // hold more than one version of the same base accession, and
        // returning whichever the hasher landed on first was nondeterministic
        // across runs (#389 follow-up).
        if let Some(base) = resolved.split('.').next() {
            let mut best: Option<(u32, &CdotTranscript)> = None;
            for (k, v) in alt {
                if !(k.starts_with(base) && k[base.len()..].starts_with('.')) {
                    continue;
                }
                let ver = k[base.len() + 1..]
                    .split('.')
                    .next()
                    .and_then(|s| s.parse::<u32>().ok())
                    .unwrap_or(0);
                if best.is_none_or(|(cur, _)| ver > cur) {
                    best = Some((ver, v));
                }
            }
            if let Some((_, tx)) = best {
                return Some(tx);
            }
        }
        None
    }

    /// List every genome build that has data for the given transcript.
    ///
    /// Includes the primary build if the transcript is present there. The
    /// returned list is unordered; callers that want a deterministic probe
    /// order (e.g. GRCh38-then-GRCh37) should sort externally.
    pub fn available_builds_for(&self, accession: &str) -> Vec<String> {
        let mut builds = Vec::new();
        // Only claim the primary build when we actually know what it is —
        // a bincode snapshot loaded without a recorded build returns
        // `None` here and we omit it rather than fabricate one
        // (#389 follow-up).
        if let Some(primary) = self.primary_build.clone() {
            if self.get_transcript(accession).is_some() {
                builds.push(primary);
            }
        }
        for (build, alt) in &self.alt_build_transcripts {
            // Direct hit, or versioned-base fallback.
            let hit = alt.contains_key(accession)
                || accession.split('.').next().is_some_and(|base| {
                    alt.keys()
                        .any(|k| k.starts_with(base) && k[base.len()..].starts_with('.'))
                });
            if hit {
                builds.push(build.clone());
            }
        }
        builds
    }

    /// Get all transcripts on a contig. Resolves aliases (e.g., "chr7" -> "NC_000007.14").
    pub fn transcripts_on_contig(&self, contig: &str) -> Vec<&str> {
        self.contig_index
            .get(self.resolve_contig(contig))
            .map(|ids| ids.iter().map(|s| s.as_str()).collect())
            .unwrap_or_default()
    }

    /// Find transcripts overlapping a genomic position. Resolves contig aliases.
    ///
    /// Backed by a per-contig SuperIntervals stab-query index that is built
    /// lazily on first call and cached for the lifetime of the `CdotMapper`
    /// (invalidated by `add_transcript`). Replaces the previous linear scan
    /// that re-folded every transcript's exon table on every query.
    pub fn transcripts_at_position(&self, contig: &str, pos: u64) -> Vec<(&str, &CdotTranscript)> {
        let canonical = self.resolve_contig(contig);
        let by_contig = self
            .contig_query_index
            .get_or_init(|| self.build_query_index());
        let Some(im) = by_contig.get(canonical) else {
            return Vec::new();
        };
        let Some(accessions) = self.contig_index.get(canonical) else {
            return Vec::new();
        };
        // i32 fits every chromosome in any current build (max ~250M < 2.1B).
        // Positions outside i32 cannot overlap a real transcript span anyway.
        let Ok(p) = i32::try_from(pos) else {
            return Vec::new();
        };
        let mut hits: Vec<u32> = Vec::with_capacity(16);
        im.search_stabbed(p, &mut hits);
        hits.into_iter()
            .filter_map(|idx| {
                let acc = accessions.get(idx as usize)?;
                let tx = self.transcripts.get(acc)?;
                Some((acc.as_str(), tx))
            })
            .collect()
    }

    /// Genomic extent `(min_exon_start, max_exon_end)` of a transcript on
    /// the mapper's primary build, taken from a cached side-table that's
    /// built lazily on first call alongside the contig stab-query index.
    ///
    /// Used by `VariantProjector::project_single_inner` instead of folding
    /// `tx.exons.iter().map(|e| e[0]).min()/max()` per call — the inputs are
    /// stable for the life of the `CdotMapper`, so re-computing every time
    /// (~4.4M folds across the SNP fixture) is pure waste.
    pub fn transcript_genome_span(&self, transcript_id: &str) -> Option<(u64, u64)> {
        let spans = self
            .transcript_genome_spans
            .get_or_init(|| self.build_transcript_genome_spans());
        spans.get(transcript_id).copied()
    }

    /// Build-aware variant of [`transcript_genome_span`]: when `build`
    /// matches the primary build the cached side-table is used; otherwise
    /// the span is computed on demand from the alt-build view returned by
    /// [`get_transcript_on_build`]. Returns `None` if cdot has no
    /// alignment for the requested build (or the alignment has no exons).
    ///
    /// The alt-build branch folds the exon table on every call rather
    /// than caching like the primary side-table. This is an acceptable
    /// trade-off until profiling shows otherwise; for callers that
    /// routinely project alt-build inputs against a multi-build cdot it
    /// becomes a per-call O(exons) cost worth caching parallel to
    /// [`Self::transcript_genome_spans`].
    pub fn transcript_genome_span_on_build(
        &self,
        transcript_id: &str,
        build: &str,
    ) -> Option<(u64, u64)> {
        // Mirror the primary-vs-alt routing in `get_transcript_on_build`:
        // unknown primary build (bincode snapshot without the field) means
        // we only have `transcripts` to work with, so use the cached
        // primary-side span (#389 follow-up).
        let primary_matches = match self.primary_build.as_deref() {
            Some(pb) => pb == build,
            None => true,
        };
        if primary_matches {
            return self.transcript_genome_span(transcript_id);
        }
        let tx = self.get_transcript_on_build(transcript_id, build)?;
        let min = tx.exons.iter().map(|e| e[0]).min()?;
        let max = tx.exons.iter().map(|e| e[1]).max()?;
        Some((min, max))
    }

    /// Build the per-transcript span side-table. Same single-pass shape as
    /// `build_query_index` so the two views can't disagree about what
    /// "transcript span" means.
    fn build_transcript_genome_spans(&self) -> HashMap<String, (u64, u64)> {
        let mut out = HashMap::with_capacity(self.transcripts.len());
        for (acc, tx) in &self.transcripts {
            if tx.exons.is_empty() {
                continue;
            }
            let min = tx.exons.iter().map(|e| e[0]).min().unwrap();
            let max = tx.exons.iter().map(|e| e[1]).max().unwrap();
            out.insert(acc.clone(), (min, max));
        }
        out
    }

    /// Construct the per-contig stab-query index. Called at most once per
    /// `CdotMapper` lifetime (until `add_transcript` invalidates the cache).
    ///
    /// SuperIntervals' intervals are end-inclusive (`[start, end]`). To
    /// preserve the old half-open `pos >= min && pos < max` semantics we
    /// insert `[min, max - 1]` and probe with `search_stabbed(pos)`.
    fn build_query_index(&self) -> HashMap<String, IntervalMap<u32>> {
        let mut by_contig: HashMap<String, IntervalMap<u32>> =
            HashMap::with_capacity(self.contig_index.len());
        for (contig, accessions) in &self.contig_index {
            let mut im: IntervalMap<u32> = IntervalMap::new();
            for (idx, acc) in accessions.iter().enumerate() {
                let Some(tx) = self.transcripts.get(acc) else {
                    continue;
                };
                if tx.exons.is_empty() {
                    continue;
                }
                let min = tx.exons.iter().map(|e| e[0]).min().unwrap();
                let max = tx.exons.iter().map(|e| e[1]).max().unwrap();
                if max <= min {
                    continue;
                }
                let (Ok(s), Ok(e)) = (i32::try_from(min), i32::try_from(max - 1)) else {
                    continue;
                };
                im.add(s, e, idx as u32);
            }
            im.build();
            by_contig.insert(contig.clone(), im);
        }
        by_contig
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

    /// c15 early-termination relies on the cdot sort invariant — exons are
    /// in `tx_start`-ascending order, which means `genome_start` is
    /// monotonically increasing on plus-strand and decreasing on
    /// minus-strand. These tests pin both the off-transcript None return and
    /// the in-transcript Some return so a future change to either the sort
    /// invariant or the early-break condition is caught.
    #[test]
    fn test_exon_for_genome_pos_plus_strand_early_term() {
        let tx = sample_transcript(); // exons at [1000,1100), [2000,2200), [3000,3150)
                                      // Below all exons → None.
        assert!(tx.exon_for_genome_pos(500).is_none());
        // Between exon 0 and 1 (intronic) → None.
        assert!(tx.exon_for_genome_pos(1500).is_none());
        // Above all exons → None.
        assert!(tx.exon_for_genome_pos(5000).is_none());
        // Inside each exon → Some(...).
        assert!(tx.exon_for_genome_pos(1050).is_some());
        assert!(tx.exon_for_genome_pos(2100).is_some());
        assert!(tx.exon_for_genome_pos(3100).is_some());
    }

    #[test]
    fn test_exon_for_genome_pos_minus_strand_early_term() {
        let tx = minus_strand_transcript();
        // The minus-strand fixture has exons (in tx-order, i.e. high→low
        // genome): [3000,3150), [2000,2200), [1000,1100).
        assert!(tx.exon_for_genome_pos(5000).is_none()); // above all
        assert!(tx.exon_for_genome_pos(2500).is_none()); // intronic between 3000 and 2200
        assert!(tx.exon_for_genome_pos(500).is_none()); // below all
        assert!(tx.exon_for_genome_pos(3050).is_some()); // first exon (highest)
        assert!(tx.exon_for_genome_pos(2100).is_some()); // middle exon
        assert!(tx.exon_for_genome_pos(1050).is_some()); // last exon (lowest)
    }

    #[test]
    fn test_locate_genome_pos_matches_exon_lookup() {
        // `locate_genome_pos` must agree with `exon_for_genome_pos` +
        // `genome_to_tx` on which exon (and what tx position) a query maps
        // to. Sweep a representative range of positions and require either
        // both succeed with the same exon number or both return None.
        for tx in [sample_transcript(), minus_strand_transcript()] {
            for pos in (500..5500).step_by(37) {
                let combined = tx.locate_genome_pos(pos);
                let separate = match (tx.genome_to_tx(pos), tx.exon_for_genome_pos(pos)) {
                    (Some(tp), Some(ex)) => Some((tp, ex)),
                    _ => None,
                };
                assert_eq!(
                    combined.map(|(tp, ex)| (tp, ex.number)),
                    separate.map(|(tp, ex)| (tp, ex.number)),
                    "disagreement at pos {pos} for {:?}-strand",
                    tx.strand
                );
            }
        }
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
    fn test_contig_alias_lookup_grch38() {
        // Transcripts indexed by RefSeq accession should also be found by UCSC name.
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();

        // Lookup by RefSeq accession (canonical key)
        assert_eq!(mapper.transcripts_on_contig("NC_000017.11").len(), 1);
        // Lookup by UCSC alias
        assert_eq!(mapper.transcripts_on_contig("chr17").len(), 1);
        // Lookup by Ensembl alias
        assert_eq!(mapper.transcripts_on_contig("17").len(), 1);
        // Non-existent contig
        assert!(mapper.transcripts_on_contig("chr99").is_empty());
    }

    #[test]
    fn test_contig_alias_lookup_grch37() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "+",
                            "exons": [[48263025, 48263098, 1, 0, 73, "M73"]]
                        }
                    },
                    "start_codon": 10,
                    "stop_codon": 60
                }
            }
        }
        "#;
        let mapper = CdotMapper::from_reader_with_build(json.as_bytes(), "GRCh37").unwrap();

        assert_eq!(mapper.transcripts_on_contig("NC_000017.10").len(), 1);
        assert_eq!(mapper.transcripts_on_contig("chr17").len(), 1);
        assert_eq!(mapper.transcripts_on_contig("17").len(), 1);

        // Verify the transcript was fully parsed (not dropped due to bad exon format)
        let tx = mapper.get_transcript("NM_000088.3").unwrap();
        assert_eq!(tx.exons.len(), 1);
        assert_eq!(tx.exons[0][0], 48263025); // genome_start
        assert_eq!(tx.exons[0][1], 48263098); // genome_end
    }

    #[test]
    fn test_from_transcripts_populates_contig_aliases_grch38() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(11),
            cds_end: Some(60),
            exons: vec![Exon::with_genomic(1, 1, 73, 50184096, 50184168)],
            chromosome: Some("NC_000017.11".to_string()),
            genomic_start: Some(50184096),
            genomic_end: Some(50184168),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts(std::iter::once(&tx));

        // Canonical RefSeq lookup works.
        assert_eq!(mapper.transcripts_on_contig("NC_000017.11").len(), 1);
        // UCSC alias resolution (default GRCh38).
        assert_eq!(mapper.transcripts_on_contig("chr17").len(), 1);
        // Ensembl alias resolution.
        assert_eq!(mapper.transcripts_on_contig("17").len(), 1);
    }

    #[test]
    fn test_from_transcripts_with_build_grch37_aliases() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(11),
            cds_end: Some(60),
            exons: vec![Exon::with_genomic(1, 1, 73, 48263025, 48263097)],
            chromosome: Some("NC_000017.10".to_string()),
            genomic_start: Some(48263025),
            genomic_end: Some(48263097),
            genome_build: GenomeBuild::GRCh37,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts_with_build(std::iter::once(&tx), "GRCh37");

        assert_eq!(mapper.transcripts_on_contig("NC_000017.10").len(), 1);
        assert_eq!(mapper.transcripts_on_contig("chr17").len(), 1);
        assert_eq!(mapper.transcripts_on_contig("17").len(), 1);
    }

    #[test]
    fn test_from_transcripts_skips_missing_chromosome() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        let tx = Transcript {
            id: "NM_NOCHR.1".to_string(),
            gene_symbol: Some("X".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::with_genomic(1, 1, 9, 1000, 1008)],
            chromosome: None, // ← missing chromosome
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts(std::iter::once(&tx));
        assert_eq!(
            mapper.transcript_count(),
            0,
            "transcript with no chromosome must be skipped, not silently coerced"
        );
    }

    #[test]
    fn test_from_transcripts_skips_exon_missing_genomic_coords() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        // Exon::new(...) leaves genomic_start/genomic_end as None.
        let tx = Transcript {
            id: "NM_NOGENOMIC.1".to_string(),
            gene_symbol: Some("X".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::new(1, 1, 9)],
            chromosome: Some("chr1".to_string()),
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts(std::iter::once(&tx));
        assert_eq!(
            mapper.transcript_count(),
            0,
            "transcript with exon missing genomic coords must be skipped"
        );
    }

    #[test]
    fn test_from_transcripts_skips_exon_with_zero_start() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        // Exon.start is documented as 1-based; 0 is invalid. We must NOT
        // silently coerce it to tx_start = 0 (it would collide with a real
        // first exon and break downstream coordinate math).
        let tx = Transcript {
            id: "NM_ZEROSTART.1".to_string(),
            gene_symbol: Some("X".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            exons: vec![Exon::with_genomic(1, 0, 9, 1000, 1008)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts(std::iter::once(&tx));
        assert_eq!(
            mapper.transcript_count(),
            0,
            "transcript with 1-based start=0 must be skipped"
        );
    }

    #[test]
    fn test_from_transcripts_skips_inverted_exon() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        // end < start on either axis is invalid (empty or inverted interval).
        let tx = Transcript {
            id: "NM_INVERTED.1".to_string(),
            gene_symbol: Some("X".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(1),
            cds_end: Some(9),
            // Genomic end (1000) < genomic start (1008).
            exons: vec![Exon::with_genomic(1, 1, 9, 1008, 1000)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts(std::iter::once(&tx));
        assert_eq!(
            mapper.transcript_count(),
            0,
            "transcript with inverted exon genomic interval must be skipped"
        );
    }

    #[test]
    fn test_from_transcripts_cds_start_zero_yields_none() {
        use crate::reference::transcript::{Exon, ManeStatus, Transcript};
        use crate::reference::Strand as TxStrand;
        use std::sync::OnceLock;

        // tx.cds_start is documented as 1-based; if upstream feeds Some(0)
        // we should yield None instead of wrapping to u64::MAX via saturating_sub.
        let tx = Transcript {
            id: "NM_CDSZERO.1".to_string(),
            gene_symbol: Some("X".to_string()),
            strand: TxStrand::Plus,
            sequence: None,
            cds_start: Some(0),
            cds_end: Some(9),
            exons: vec![Exon::with_genomic(1, 1, 9, 1000, 1008)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(1008),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let mapper = CdotMapper::from_transcripts(std::iter::once(&tx));
        let stored = mapper
            .get_transcript("NM_CDSZERO.1")
            .expect("transcript should be present (exon coords are valid)");
        assert_eq!(
            stored.cds_start, None,
            "cds_start=Some(0) must become None, not a wrapped u64"
        );
    }

    #[test]
    fn test_contig_alias_transcripts_at_position() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();

        // Position within transcript, looked up by UCSC alias
        let results = mapper.transcripts_at_position("chr17", 50184100);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].0, "NM_000088.3");

        // Position outside transcript
        let results = mapper.transcripts_at_position("chr17", 99999999);
        assert!(results.is_empty());
    }

    #[test]
    fn test_contig_alias_chrm() {
        let json = r#"
        {
            "transcripts": {
                "NM_MITO.1": {
                    "gene_name": "MT_TEST",
                    "contig": "NC_012920.1",
                    "strand": "+",
                    "exons": [[100, 200, 0, 100]],
                    "cds_start": 10,
                    "cds_end": 90
                }
            }
        }
        "#;
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();

        assert_eq!(mapper.transcripts_on_contig("NC_012920.1").len(), 1);
        assert_eq!(mapper.transcripts_on_contig("chrM").len(), 1);
        assert_eq!(mapper.transcripts_on_contig("MT").len(), 1);
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
    // Bincode serialization tests
    // =========================================================================

    #[test]
    fn test_bincode_roundtrip() {
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
                },
                "NM_000059.4": {
                    "gene_name": "BRCA2",
                    "contig": "NC_000013.11",
                    "strand": "+",
                    "exons": [[32315474, 32315667, 0, 193]],
                    "cds_start": 40,
                    "cds_end": 150
                }
            }
        }
        "#;
        let original = CdotMapper::from_reader(json.as_bytes()).unwrap();
        assert_eq!(original.transcript_count(), 2);

        // Roundtrip through bincode file
        let temp_dir = tempfile::TempDir::new().unwrap();
        let bin_path = temp_dir.path().join("roundtrip.bin");
        original.to_bincode_file(&bin_path).unwrap();
        let decoded = CdotMapper::from_bincode_file(&bin_path).unwrap();

        assert_eq!(decoded.transcript_count(), 2);
        let tx = decoded.get_transcript("NM_000088.3").unwrap();
        assert_eq!(tx.gene_name.as_deref(), Some("COL1A1"));
        assert_eq!(tx.contig, "NC_000017.11");
        assert_eq!(tx.strand, Strand::Plus);
        assert_eq!(tx.exons.len(), 2);

        let tx2 = decoded.get_transcript("NM_000059.4").unwrap();
        assert_eq!(tx2.gene_name.as_deref(), Some("BRCA2"));

        // Contig index should be preserved
        assert_eq!(decoded.transcripts_on_contig("NC_000017.11").len(), 1);
        assert_eq!(decoded.transcripts_on_contig("NC_000013.11").len(), 1);

        // Base-to-versioned index should be preserved
        assert!(decoded.get_transcript("NM_000088.4").is_some()); // falls back to .3
    }

    #[test]
    fn test_bincode_file_roundtrip() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let original = CdotMapper::from_reader(json.as_bytes()).unwrap();

        let temp_dir = tempfile::TempDir::new().unwrap();
        let bin_path = temp_dir.path().join("test.bin");

        original.to_bincode_file(&bin_path).unwrap();
        assert!(bin_path.exists());

        let loaded = CdotMapper::from_bincode_file(&bin_path).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
        assert!(loaded.get_transcript("NM_000088.3").is_some());
    }

    #[test]
    fn test_concurrent_to_bincode_file_same_path() {
        // Multiple threads writing the same cache path concurrently must not
        // collide on their temp siblings or leave a half-written file: the
        // final cache must always be loadable and complete.
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let original = std::sync::Arc::new(CdotMapper::from_reader(json.as_bytes()).unwrap());

        let temp_dir = tempfile::TempDir::new().unwrap();
        let bin_path = std::sync::Arc::new(temp_dir.path().join("concurrent.bin"));

        let handles: Vec<_> = (0..8)
            .map(|_| {
                let original = std::sync::Arc::clone(&original);
                let bin_path = std::sync::Arc::clone(&bin_path);
                std::thread::spawn(move || original.to_bincode_file(bin_path.as_path()))
            })
            .collect();
        for handle in handles {
            handle
                .join()
                .unwrap()
                .expect("concurrent write should succeed");
        }

        // The final cache is complete and loadable, and no temp siblings leaked.
        let loaded = CdotMapper::from_bincode_file(bin_path.as_path()).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
        assert!(loaded.get_transcript("NM_000088.3").is_some());

        let leftover_tmp: Vec<_> = std::fs::read_dir(temp_dir.path())
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().contains(".tmp."))
            .collect();
        assert!(
            leftover_tmp.is_empty(),
            "temp siblings leaked: {:?}",
            leftover_tmp
        );
    }

    #[test]
    fn test_load_prefers_bincode() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let original = CdotMapper::from_reader(json.as_bytes()).unwrap();

        // Mutate the mapper so bincode content differs from JSON
        let mut from_bin = original.clone();
        from_bin
            .transcripts
            .get_mut("NM_000088.3")
            .unwrap()
            .gene_name = Some("FROM_BIN".to_string());

        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");
        let bin_path = temp_dir.path().join("cdot.bin");

        // Write JSON file (gene_name = "COL1A1")
        std::fs::write(&json_path, json).unwrap();

        // Write bincode file (gene_name = "FROM_BIN")
        from_bin.to_bincode_file(&bin_path).unwrap();

        // load() given the JSON path should find and use the .bin sibling
        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
        assert_eq!(
            loaded
                .get_transcript("NM_000088.3")
                .unwrap()
                .gene_name
                .as_deref(),
            Some("FROM_BIN"),
            "load() should prefer bincode over JSON"
        );
    }

    #[test]
    fn test_load_falls_back_to_json() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");

        // Write only JSON, no bincode
        std::fs::write(&json_path, json).unwrap();

        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
    }

    #[test]
    fn test_load_falls_back_to_json_on_corrupt_bincode() {
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "contig": "NC_000017.11",
                    "strand": "+",
                    "exons": [[50184096, 50184169, 0, 73]],
                    "cds_start": 10,
                    "cds_end": 60
                }
            }
        }
        "#;
        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");
        let bin_path = temp_dir.path().join("cdot.bin");

        // Write valid JSON and corrupt bincode
        std::fs::write(&json_path, json).unwrap();
        std::fs::write(&bin_path, b"not valid bincode data").unwrap();

        // load() should fall back to JSON despite corrupt .bin
        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
        assert!(loaded.get_transcript("NM_000088.3").is_some());

        // Self-heal: the corrupt cache must have been refreshed to a current
        // one, so a second load takes the fast bincode path.
        assert!(
            CdotMapper::bincode_is_current(&bin_path),
            "corrupt cache should be refreshed after the JSON fallback"
        );
        let reloaded = CdotMapper::from_bincode_file(&bin_path).unwrap();
        assert_eq!(reloaded.transcript_count(), 1);
    }

    #[test]
    fn test_bincode_roundtrip_preserves_alt_build() {
        // A transcript with both GRCh37 and GRCh38 alignments. Primary build is
        // GRCh38, so GRCh37 lives in the alt-build map. The fast (bincode) path
        // must preserve it, or NG/NC-parent build resolution (#332) silently
        // breaks whenever the cache is used.
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "+",
                            "exons": [[48263025, 48263098, 1, 0, 73, "M73"]]
                        },
                        "GRCh38": {
                            "contig": "NC_000017.11",
                            "strand": "+",
                            "exons": [[50184096, 50184169, 1, 0, 73, "M73"]]
                        }
                    },
                    "start_codon": 10,
                    "stop_codon": 60
                }
            }
        }
        "#;
        let from_json = CdotMapper::from_reader(json.as_bytes()).unwrap();
        let temp = tempfile::TempDir::new().unwrap();
        let p = temp.path().join("cdot.bin");
        from_json.to_bincode_file(&p).unwrap();
        let from_bin = CdotMapper::from_bincode_file(&p).unwrap();

        // GRCh37 (alt-build) resolution must be identical across both paths.
        let want = from_json
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .map(|t| (t.contig.clone(), t.cds_start));
        let got = from_bin
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .map(|t| (t.contig.clone(), t.cds_start));
        assert_eq!(
            want, got,
            "alt-build resolution must survive bincode roundtrip"
        );
        assert_eq!(got, Some(("NC_000017.10".to_string(), Some(10))));
    }

    #[test]
    fn test_bincode_header_roundtrip_and_staleness() {
        let temp = tempfile::TempDir::new().unwrap();
        let p = temp.path().join("cdot.bin");

        let mut m = CdotMapper::new();
        m.add_transcript("NM_000001.1".to_string(), sample_transcript());
        m.to_bincode_file(&p).unwrap();

        // A freshly written cache carries the current magic+version and loads.
        assert!(CdotMapper::bincode_is_current(&p));
        assert_eq!(
            CdotMapper::from_bincode_file(&p)
                .unwrap()
                .transcript_count(),
            1
        );

        // A legacy / headerless cache is detected as not-current and rejected
        // cleanly, rather than mis-deserialized.
        std::fs::write(&p, b"legacy-bincode-without-any-header-bytes").unwrap();
        assert!(!CdotMapper::bincode_is_current(&p));
        assert!(CdotMapper::from_bincode_file(&p).is_err());

        // Correct magic but a future schema version is also rejected.
        let mut wrong_version = Vec::new();
        wrong_version.extend_from_slice(&CDOT_BINCODE_MAGIC);
        wrong_version.extend_from_slice(&(CDOT_BINCODE_VERSION + 1).to_le_bytes());
        wrong_version.extend_from_slice(&[0u8; 16]);
        std::fs::write(&p, &wrong_version).unwrap();
        assert!(!CdotMapper::bincode_is_current(&p));
        assert!(CdotMapper::from_bincode_file(&p).is_err());
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

    // -------------------------------------------------------------------------
    // Issue #332: multi-build retention. A cdot JSON with `genome_builds` for
    // both GRCh37 and GRCh38 must yield a mapper that can return the
    // transcript on either build via `get_transcript_on_build`.
    // -------------------------------------------------------------------------

    fn multi_build_cdot_json() -> &'static str {
        r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "+",
                            "exons": [[48263025, 48263098, 1, 0, 73, "M73"]]
                        },
                        "GRCh38": {
                            "contig": "NC_000017.11",
                            "strand": "+",
                            "exons": [[50184096, 50184169, 1, 0, 73, "M73"]]
                        }
                    },
                    "start_codon": 10,
                    "stop_codon": 60
                }
            }
        }
        "#
    }

    #[test]
    fn test_get_transcript_on_build_returns_correct_contig() {
        let json = multi_build_cdot_json();
        // Primary build = GRCh38.
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();

        let tx38 = mapper
            .get_transcript_on_build("NM_000088.3", "GRCh38")
            .expect("GRCh38 view must be present");
        assert_eq!(tx38.contig, "NC_000017.11");

        let tx37 = mapper
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .expect("GRCh37 view must be retained even when GRCh38 is primary");
        assert_eq!(tx37.contig, "NC_000017.10");
    }

    #[test]
    fn test_get_transcript_on_build_primary_build_grch37() {
        let json = multi_build_cdot_json();
        // Primary build = GRCh37.
        let mapper = CdotMapper::from_reader_with_build(json.as_bytes(), "GRCh37").unwrap();

        let tx37 = mapper
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .expect("primary build access via get_transcript_on_build");
        assert_eq!(tx37.contig, "NC_000017.10");

        let tx38 = mapper
            .get_transcript_on_build("NM_000088.3", "GRCh38")
            .expect("non-primary GRCh38 view must be retained");
        assert_eq!(tx38.contig, "NC_000017.11");
    }

    #[test]
    fn test_available_builds_for_multi_build_input() {
        let json = multi_build_cdot_json();
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();
        let mut builds = mapper.available_builds_for("NM_000088.3");
        builds.sort();
        assert_eq!(builds, vec!["GRCh37".to_string(), "GRCh38".to_string()]);
    }

    #[test]
    fn test_available_builds_for_unknown_accession() {
        let json = multi_build_cdot_json();
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();
        assert!(mapper.available_builds_for("NM_999999.1").is_empty());
    }

    #[test]
    fn test_get_transcript_on_build_version_fallback_in_alt() {
        // NM_000088.4 falls back to NM_000088.3 inside the alt-build map too.
        let json = multi_build_cdot_json();
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();
        let tx37 = mapper
            .get_transcript_on_build("NM_000088.4", "GRCh37")
            .expect("version fallback inside alt-build map");
        assert_eq!(tx37.contig, "NC_000017.10");
    }

    #[test]
    fn test_get_transcript_on_build_alt_version_fallback_picks_highest() {
        // When the alt-build map holds multiple versions of the same base
        // (e.g. NM_000088.3 and NM_000088.4), the version fallback must
        // deterministically prefer the highest numeric version. Pre-fix this
        // iterated the HashMap and returned whichever versioned key the hasher
        // landed on first.
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "+",
                            "exons": [[48263025, 48263098, 1, 0, 73, "M73"]]
                        },
                        "GRCh38": {
                            "contig": "NC_000017.11",
                            "strand": "+",
                            "exons": [[50184096, 50184169, 1, 0, 73, "M73"]]
                        }
                    },
                    "start_codon": 10,
                    "stop_codon": 60
                },
                "NM_000088.4": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "+",
                            "exons": [[48263025, 48263198, 1, 0, 173, "M173"]]
                        },
                        "GRCh38": {
                            "contig": "NC_000017.11",
                            "strand": "+",
                            "exons": [[50184096, 50184269, 1, 0, 173, "M173"]]
                        }
                    },
                    "start_codon": 20,
                    "stop_codon": 160
                }
            }
        }
        "#;
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();
        // Primary build = GRCh38, so the GRCh37 alt map holds both .3 and .4.
        // Querying with .2 (absent) forces the version fallback.
        let tx = mapper
            .get_transcript_on_build("NM_000088.2", "GRCh37")
            .expect("version fallback must find a sibling");
        assert_eq!(
            tx.cds_start,
            Some(20),
            "alt-build fallback must deterministically pick the highest version (.4 over .3)"
        );
    }

    #[test]
    fn test_protein_accession_parsed_from_cdot() {
        // cdot records the CDS protein accession per transcript. The loader
        // must surface it on `CdotTranscript.protein` so the projector emits
        // the correct `p.` accession (here NP_000068.1) instead of falling
        // back to the NM_* -> NP_* number-preserving inference (NP_000077.4),
        // which is wrong whenever the NP number differs from the NM number.
        let json = r#"
        {
            "transcripts": {
                "NM_000077.4": {
                    "gene_name": "CDKN2A",
                    "protein": "NP_000068.1",
                    "genome_builds": {
                        "GRCh38": {
                            "contig": "NC_000009.12",
                            "strand": "+",
                            "exons": [[21967750, 21968240, 1, 0, 490, "M490"]]
                        }
                    },
                    "start_codon": 10,
                    "stop_codon": 60
                }
            }
        }
        "#;
        let mapper = CdotMapper::from_reader(json.as_bytes()).unwrap();
        let tx = mapper
            .get_transcript_on_build("NM_000077.4", "GRCh38")
            .expect("transcript present");
        assert_eq!(
            tx.protein.as_deref(),
            Some("NP_000068.1"),
            "cdot `protein` field must be plumbed onto CdotTranscript.protein"
        );
    }

    // -------------------------------------------------------------------------
    // #389 follow-up: primary_build() returns Option, and bincode loads that
    // pre-date the primary_build field must surface as None rather than
    // fabricating "GRCh38" (which mis-stamps GRCh37 snapshots).
    // -------------------------------------------------------------------------

    #[test]
    fn test_primary_build_known_after_explicit_build_load() {
        let mapper =
            CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh37")
                .unwrap();
        assert_eq!(mapper.primary_build(), Some("GRCh37"));
    }

    #[test]
    fn test_bincode_roundtrip_preserves_primary_build_grch37() {
        // Persist a GRCh37-primary mapper to bincode and read it back; the
        // build name must survive the round-trip so downstream callers
        // don't mis-stamp transcripts as GRCh38.
        let mapper =
            CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh37")
                .unwrap();
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("roundtrip.bin");
        mapper.to_bincode_file(&path).expect("write bincode");
        let reloaded = CdotMapper::from_bincode_file(&path).expect("read bincode");
        assert_eq!(
            reloaded.primary_build(),
            Some("GRCh37"),
            "primary_build must survive bincode round-trip; pre-fix this returned \"GRCh38\""
        );
    }

    #[test]
    fn test_bincode_roundtrip_preserves_primary_build_grch38() {
        // Symmetric to the GRCh37 case — verifies the field is persisted
        // (rather than constant-defaulted) by exercising the non-default
        // value above and the default value here.
        let mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join("roundtrip.bin");
        mapper.to_bincode_file(&path).expect("write bincode");
        let reloaded = CdotMapper::from_bincode_file(&path).expect("read bincode");
        assert_eq!(reloaded.primary_build(), Some("GRCh38"));
    }

    // -------------------------------------------------------------------------
    // #389 follow-up: build-specific lookups carry the LRG fallback.
    // Pre-fix, asking for an LRG transcript on a non-primary build returned
    // None even when the mapped RefSeq transcript existed in the alt-build
    // map; `get_transcript_on_build` short-circuited before LRG resolution.
    // -------------------------------------------------------------------------

    #[test]
    fn test_get_transcript_on_build_resolves_lrg_to_refseq_on_alt_build() {
        // Build a multi-build mapper and register an LRG → RefSeq mapping.
        let mut mapper =
            CdotMapper::from_reader_with_build(multi_build_cdot_json().as_bytes(), "GRCh38")
                .unwrap();
        mapper
            .lrg_to_refseq
            .insert("LRG_1t1".to_string(), "NM_000088.3".to_string());

        // Primary build: LRG lookup works via the get_transcript LRG path.
        let tx38 = mapper
            .get_transcript_on_build("LRG_1t1", "GRCh38")
            .expect("LRG → RefSeq on primary build");
        assert_eq!(tx38.contig, "NC_000017.11");

        // Alt build (the regression): LRG lookup must resolve to the
        // RefSeq accession before consulting alt_build_transcripts.
        let tx37 = mapper
            .get_transcript_on_build("LRG_1t1", "GRCh37")
            .expect("LRG → RefSeq must work on alt build too (pre-fix this returned None)");
        assert_eq!(tx37.contig, "NC_000017.10");
    }
}
