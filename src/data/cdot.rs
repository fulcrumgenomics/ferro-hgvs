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
use std::path::{Path, PathBuf};
use superintervals::IntervalMap;

/// rkyv-archivable mirror of the cdot cache. Kept fully separate from the
/// domain types (`CdotTranscript`, `Strand`, `CigarOp`) so the rkyv derives —
/// and the choice of an `i8` strand / flattened cigar encoding — stay contained
/// to the cache layer. `from_rkyv_file`/`to_rkyv_file` convert at the boundary.
mod rkyv_cache {
    use super::{CdotMapper, CdotTranscript, CigarOp};
    use crate::error::FerroError;
    use crate::reference::Strand;
    use rkyv::{Archive, Deserialize, Serialize};
    use std::collections::HashMap;

    /// Bump whenever the layout below changes so a stale archive is rejected by
    /// the version check (in addition to rkyv's own structural validation) and
    /// regenerated rather than mis-read.
    pub(super) const RKYV_FORMAT_VERSION: u32 = 1;

    #[derive(Archive, Serialize, Deserialize)]
    pub(super) struct RkyvCigar {
        /// 0 = Match, 1 = Insertion, 2 = Deletion.
        op: u8,
        len: u64,
    }

    #[derive(Archive, Serialize, Deserialize)]
    pub(super) struct RkyvTx {
        gene_name: Option<String>,
        contig: String,
        /// 1 = Plus, -1 = Minus, 0 = Unknown.
        strand: i8,
        exons: Vec<[u64; 4]>,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
        exon_cigars: Vec<Option<Vec<RkyvCigar>>>,
        gene_id: Option<String>,
        protein: Option<String>,
    }

    #[derive(Archive, Serialize, Deserialize)]
    pub(super) struct RkyvSnapshot {
        /// Format version — first field so a mismatch is a clean signal.
        pub(super) format_version: u32,
        pub(super) transcripts: HashMap<String, RkyvTx>,
        pub(super) contig_index: HashMap<String, Vec<String>>,
        pub(super) contig_alias_to_canonical: HashMap<String, String>,
        pub(super) base_to_versioned: HashMap<String, String>,
        pub(super) lrg_to_refseq: HashMap<String, String>,
        pub(super) primary_build: Option<String>,
        pub(super) alt_build_transcripts: HashMap<String, HashMap<String, RkyvTx>>,
    }

    fn strand_to_i8(s: Strand) -> i8 {
        match s {
            Strand::Plus => 1,
            Strand::Minus => -1,
            Strand::Unknown => 0,
        }
    }
    /// Reject unknown tags rather than coercing them: a corrupt `.rkyv` whose
    /// strand byte is out of range must be regenerated, not silently read back
    /// as `Strand::Unknown`.
    fn strand_from_i8(v: i8) -> Result<Strand, FerroError> {
        match v {
            1 => Ok(Strand::Plus),
            -1 => Ok(Strand::Minus),
            0 => Ok(Strand::Unknown),
            other => Err(FerroError::Io {
                msg: format!("cdot rkyv archive has invalid strand tag {other}"),
            }),
        }
    }
    fn cigar_to_rkyv(c: &CigarOp) -> RkyvCigar {
        match *c {
            CigarOp::Match(n) => RkyvCigar { op: 0, len: n },
            CigarOp::Insertion(n) => RkyvCigar { op: 1, len: n },
            CigarOp::Deletion(n) => RkyvCigar { op: 2, len: n },
        }
    }

    pub(super) fn tx_to_rkyv(t: &CdotTranscript) -> RkyvTx {
        RkyvTx {
            gene_name: t.gene_name.clone(),
            contig: t.contig.clone(),
            strand: strand_to_i8(t.strand),
            exons: t.exons.clone(),
            cds_start: t.cds_start,
            cds_end: t.cds_end,
            exon_cigars: t
                .exon_cigars
                .iter()
                .map(|opt| opt.as_ref().map(|v| v.iter().map(cigar_to_rkyv).collect()))
                .collect(),
            gene_id: t.gene_id.clone(),
            protein: t.protein.clone(),
        }
    }

    /// Reject unknown op tags rather than coercing them to `Match`, for the same
    /// reason as [`strand_from_i8`].
    fn cigar_from_archived(a: &ArchivedRkyvCigar) -> Result<CigarOp, FerroError> {
        let len = a.len.to_native();
        match a.op {
            0 => Ok(CigarOp::Match(len)),
            1 => Ok(CigarOp::Insertion(len)),
            2 => Ok(CigarOp::Deletion(len)),
            other => Err(FerroError::Io {
                msg: format!("cdot rkyv archive has invalid cigar op tag {other}"),
            }),
        }
    }

    /// Build a `CdotTranscript` directly from the archived form — a single pass,
    /// avoiding an intermediate owned `RkyvTx`. Returns an error if any archived
    /// enum tag (strand, cigar op) is out of range, so a corrupt archive is
    /// rejected and regenerated rather than materialized with wrong data.
    fn tx_from_archived(a: &ArchivedRkyvTx) -> Result<CdotTranscript, FerroError> {
        Ok(CdotTranscript {
            gene_name: a.gene_name.as_ref().map(|s| s.as_str().to_string()),
            contig: a.contig.as_str().to_string(),
            strand: strand_from_i8(a.strand)?,
            exons: a
                .exons
                .iter()
                .map(|e| {
                    [
                        e[0].to_native(),
                        e[1].to_native(),
                        e[2].to_native(),
                        e[3].to_native(),
                    ]
                })
                .collect(),
            cds_start: a.cds_start.as_ref().map(|v| v.to_native()),
            cds_end: a.cds_end.as_ref().map(|v| v.to_native()),
            exon_cigars: a
                .exon_cigars
                .iter()
                .map(|opt| {
                    opt.as_ref()
                        .map(|v| {
                            v.iter()
                                .map(cigar_from_archived)
                                .collect::<Result<Vec<_>, _>>()
                        })
                        .transpose()
                })
                .collect::<Result<Vec<_>, _>>()?,
            gene_id: a.gene_id.as_ref().map(|s| s.as_str().to_string()),
            protein: a.protein.as_ref().map(|s| s.as_str().to_string()),
        })
    }

    /// Owned cdot maps materialized directly from an archive (one pass).
    pub(super) struct CdotMaps {
        pub transcripts: HashMap<String, CdotTranscript>,
        pub contig_index: HashMap<String, Vec<String>>,
        pub contig_alias_to_canonical: HashMap<String, String>,
        pub base_to_versioned: HashMap<String, String>,
        pub lrg_to_refseq: HashMap<String, String>,
        pub primary_build: Option<String>,
        pub alt_build_transcripts: HashMap<String, HashMap<String, CdotTranscript>>,
    }

    pub(super) fn version_of(a: &ArchivedRkyvSnapshot) -> u32 {
        a.format_version.to_native()
    }

    /// Materialize owned cdot maps from a validated archive in a single pass.
    /// Returns an error if any transcript carries an out-of-range enum tag.
    pub(super) fn maps_from_archived(a: &ArchivedRkyvSnapshot) -> Result<CdotMaps, FerroError> {
        let txs = |m: &rkyv::collections::swiss_table::ArchivedHashMap<
            rkyv::string::ArchivedString,
            ArchivedRkyvTx,
        >|
         -> Result<HashMap<String, CdotTranscript>, FerroError> {
            m.iter()
                .map(|(k, v)| Ok((k.as_str().to_string(), tx_from_archived(v)?)))
                .collect()
        };
        let str_map = |m: &rkyv::collections::swiss_table::ArchivedHashMap<
            rkyv::string::ArchivedString,
            rkyv::string::ArchivedString,
        >|
         -> HashMap<String, String> {
            m.iter()
                .map(|(k, v)| (k.as_str().to_string(), v.as_str().to_string()))
                .collect()
        };
        Ok(CdotMaps {
            transcripts: txs(&a.transcripts)?,
            contig_index: a
                .contig_index
                .iter()
                .map(|(k, v)| {
                    (
                        k.as_str().to_string(),
                        v.iter().map(|s| s.as_str().to_string()).collect(),
                    )
                })
                .collect(),
            contig_alias_to_canonical: str_map(&a.contig_alias_to_canonical),
            base_to_versioned: str_map(&a.base_to_versioned),
            lrg_to_refseq: str_map(&a.lrg_to_refseq),
            primary_build: a.primary_build.as_ref().map(|s| s.as_str().to_string()),
            alt_build_transcripts: a
                .alt_build_transcripts
                .iter()
                .map(|(b, m)| Ok((b.as_str().to_string(), txs(m)?)))
                .collect::<Result<HashMap<_, _>, FerroError>>()?,
        })
    }

    /// Build the rkyv snapshot mirror from a populated mapper (serialize side).
    pub(super) fn snapshot_from_mapper(m: &CdotMapper) -> RkyvSnapshot {
        let map = |src: &HashMap<String, CdotTranscript>| -> HashMap<String, RkyvTx> {
            src.iter()
                .map(|(k, v)| (k.clone(), tx_to_rkyv(v)))
                .collect()
        };
        RkyvSnapshot {
            format_version: RKYV_FORMAT_VERSION,
            transcripts: map(&m.transcripts),
            contig_index: m.contig_index.clone(),
            contig_alias_to_canonical: m.contig_alias_to_canonical.clone(),
            base_to_versioned: m.base_to_versioned.clone(),
            lrg_to_refseq: m.lrg_to_refseq.clone(),
            primary_build: m.primary_build.clone(),
            alt_build_transcripts: m
                .alt_build_transcripts
                .iter()
                .map(|(b, txs)| (b.clone(), map(txs)))
                .collect(),
        }
    }

    /// Serialize a current-version archive carrying a single transcript whose
    /// `strand` / `cigar op` bytes are forced to the given (possibly invalid)
    /// raw tags. Lets tests exercise the archive-tag rejection path without
    /// going through `strand_to_i8` / `cigar_to_rkyv`, whose construction can
    /// only ever produce valid tags. The archive passes structural and version
    /// validation, so any rejection comes purely from tag checking.
    #[cfg(test)]
    pub(super) fn archive_bytes_with_tags(strand: i8, cigar_op: u8) -> Vec<u8> {
        let tx = RkyvTx {
            gene_name: None,
            contig: "NC_000017.11".to_string(),
            strand,
            exons: vec![[50184096, 50184169, 0, 73]],
            cds_start: Some(10),
            cds_end: Some(60),
            exon_cigars: vec![Some(vec![RkyvCigar {
                op: cigar_op,
                len: 73,
            }])],
            gene_id: None,
            protein: None,
        };
        let mut transcripts = HashMap::new();
        transcripts.insert("NM_000088.3".to_string(), tx);
        let snapshot = RkyvSnapshot {
            format_version: RKYV_FORMAT_VERSION,
            transcripts,
            contig_index: HashMap::new(),
            contig_alias_to_canonical: HashMap::new(),
            base_to_versioned: HashMap::new(),
            lrg_to_refseq: HashMap::new(),
            primary_build: None,
            alt_build_transcripts: HashMap::new(),
        };
        rkyv::to_bytes::<rkyv::rancor::Error>(&snapshot)
            .expect("serialize test archive")
            .to_vec()
    }
}

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

/// Parse the integer version suffix of an accession (`"NM_003002.4"` -> `Some(4)`),
/// or `None` when there is no trailing numeric `.<n>`. Used to pick the highest
/// version deterministically when resolving a base accession to a specific one.
fn accession_version(accession: &str) -> Option<u32> {
    accession
        .split('.')
        .nth(1)
        .and_then(|v| v.parse::<u32>().ok())
}

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

    /// Determine whether `genome_pos` falls strictly inside a transcript-genome
    /// CIGAR `Deletion` gap (a genome base with no transcript counterpart) for
    /// the exon at index `exon_idx` in [`Self::exons`] / [`Self::exon_cigars`].
    ///
    /// The CIGAR ops describe the alignment within the exon, walked from the
    /// exon's transcript-5' end. A `Deletion(n)` consumes `n` genome bases and
    /// `0` transcript bases — so any genome position landing inside the
    /// deletion's genome span has no well-defined transcript coordinate. Plain
    /// offset arithmetic (as in [`Self::locate_genome_pos`]) would silently map
    /// such a position to a wrong coordinate; this lets the mapping site refuse
    /// it instead.
    ///
    /// Returns `Some(CigarGap)` describing the deletion when `genome_pos` is
    /// strictly inside one, and `None` when the position is matched, when the
    /// exon has no CIGAR data, or when the position is outside the exon.
    ///
    /// `genome_pos` must lie within the exon's genome span; callers pass the
    /// exon already located by [`Self::locate_genome_pos`].
    pub(crate) fn cigar_deletion_gap_at_genome_pos(
        &self,
        exon_idx: usize,
        genome_pos: u64,
    ) -> Option<CigarGap> {
        let exon = self.exons.get(exon_idx)?;
        let (genome_start, genome_end) = (exon[0], exon[1]);
        if genome_pos < genome_start || genome_pos >= genome_end {
            return None;
        }
        let ops = match self.exon_cigars.get(exon_idx) {
            Some(Some(ops)) if !ops.is_empty() => ops,
            _ => return None,
        };

        // Genome offset along the CIGAR, measured from the exon's tx-5' end
        // (the same orientation the CIGAR is written in). On the plus strand
        // this counts up from `genome_start`; on the minus strand it counts
        // down from `genome_end - 1`. This mirrors the offset
        // `locate_genome_pos` uses to derive the (CIGAR-unaware) tx position.
        let target_genome_offset = match self.strand {
            Strand::Plus => genome_pos - genome_start,
            Strand::Minus => genome_end - 1 - genome_pos,
            Strand::Unknown => return None,
        };

        let mut genome_consumed: u64 = 0;
        for op in ops {
            match op {
                CigarOp::Match(len) => {
                    genome_consumed += len;
                }
                CigarOp::Deletion(len) => {
                    // The deletion covers genome offsets
                    // [genome_consumed, genome_consumed + len). A position
                    // landing there is strictly inside the gap.
                    if target_genome_offset >= genome_consumed
                        && target_genome_offset < genome_consumed + len
                    {
                        return Some(CigarGap {
                            kind: CigarGapKind::Deletion,
                            length: *len,
                            offset_in_gap: target_genome_offset - genome_consumed,
                        });
                    }
                    genome_consumed += len;
                }
                CigarOp::Insertion(_) => {
                    // Insertions consume transcript bases only, not genome.
                }
            }
            if genome_consumed > target_genome_offset {
                break;
            }
        }
        None
    }

    /// Determine whether transcript position `tx_pos` falls strictly inside a
    /// transcript-genome CIGAR `Insertion` gap (a transcript base with no
    /// genome counterpart) for the exon at index `exon_idx`.
    ///
    /// The mirror of [`Self::cigar_deletion_gap_at_genome_pos`] on the
    /// transcript axis: an `Insertion(n)` consumes `n` transcript bases and `0`
    /// genome bases, so a transcript position inside the insertion's span has
    /// no well-defined genome coordinate. Returns `Some(CigarGap)` with
    /// `kind = Insertion` when `tx_pos` is strictly inside such a gap.
    ///
    /// `tx_pos` is the 0-based transcript position; callers pass an exon whose
    /// `[tx_start, tx_end)` span contains it.
    pub(crate) fn cigar_insertion_gap_at_tx_pos(
        &self,
        exon_idx: usize,
        tx_pos: u64,
    ) -> Option<CigarGap> {
        let exon = self.exons.get(exon_idx)?;
        let (tx_start, tx_end) = (exon[2], exon[3]);
        if tx_pos < tx_start || tx_pos >= tx_end {
            return None;
        }
        let ops = match self.exon_cigars.get(exon_idx) {
            Some(Some(ops)) if !ops.is_empty() => ops,
            _ => return None,
        };

        // Transcript offset along the CIGAR, measured from the exon's tx-5'
        // end (the order the CIGAR is written in, on both strands).
        let target_tx_offset = tx_pos - tx_start;

        let mut tx_consumed: u64 = 0;
        for op in ops {
            match op {
                CigarOp::Match(len) => {
                    tx_consumed += len;
                }
                CigarOp::Insertion(len) => {
                    if target_tx_offset >= tx_consumed && target_tx_offset < tx_consumed + len {
                        return Some(CigarGap {
                            kind: CigarGapKind::Insertion,
                            length: *len,
                            offset_in_gap: target_tx_offset - tx_consumed,
                        });
                    }
                    tx_consumed += len;
                }
                CigarOp::Deletion(_) => {
                    // Deletions consume genome bases only, not transcript.
                }
            }
            if tx_consumed > target_tx_offset {
                break;
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

/// Which side of a transcript-genome CIGAR alignment a gap lives on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarGapKind {
    /// A genome base with no transcript counterpart (CIGAR `Deletion`).
    Deletion,
    /// A transcript base with no genome counterpart (CIGAR `Insertion`).
    Insertion,
}

/// Describes a transcript-genome CIGAR indel gap that a queried position falls
/// strictly inside. Produced by the gap-detection helpers on
/// [`CdotTranscript`] so the mapping site can build an informative error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarGap {
    /// Whether the gap is a `Deletion` (genome-only) or `Insertion` (tx-only).
    pub kind: CigarGapKind,
    /// Length of the gap op in bases.
    pub length: u64,
    /// Offset of the queried position from the start of the gap (0-based).
    pub offset_in_gap: u64,
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

/// Which path [`CdotMapper::load_with_source`] took to produce the mapper.
///
/// Makes a silent fast-path → slow-path fallback impossible to ignore at the
/// API level (the root cause of the #585 regression, where a broken binary
/// cache silently re-parsed 512 MB of JSON on every run). A binary cache
/// (`.rkyv` or legacy `.bin`) is the fast `Archive` path; an actual JSON parse
/// is `JsonFallback`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CdotLoadSource {
    /// Loaded from a fast binary archive (`.rkyv`, or legacy `.bin`).
    Archive,
    /// Fell back to parsing the source JSON (slow path).
    JsonFallback,
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
    /// Lazily-built per-contig stab-query index for the NON-primary builds.
    ///
    /// Mirrors [`Self::contig_query_index`] but covers the contigs that only
    /// appear in [`Self::alt_build_transcripts`] (e.g. the GRCh37 `NC_*.11`
    /// accessions when the primary build is GRCh38). GRCh37 and GRCh38 contig
    /// accessions are distinct, so a contig is owned by exactly one build and
    /// `transcripts_at_position` can route to the owning build purely from the
    /// contig name.
    ///
    /// Derived entirely from `alt_build_transcripts`, so it is rebuilt lazily
    /// rather than persisted — this keeps the rkyv/bincode snapshot schema
    /// unchanged (snapshots already carry `alt_build_transcripts`, so a
    /// reloaded mapper rebuilds this index on first query). Cleared alongside
    /// `contig_query_index` whenever the transcript maps change.
    ///
    /// The inner `Vec<String>` is the accession list for the contig in the
    /// owning build's map; the `IntervalMap` payload is an index into it.
    alt_build_query_index: OnceCell<HashMap<String, AltBuildContigStab>>,
    /// Lazily-built per-transcript `(min_genome_start, max_genome_end)`
    /// cache so `VariantProjector::project_single_inner` doesn't re-fold the
    /// exon table on every projection. Sharing the same `OnceCell` build
    /// trigger as `contig_query_index` keeps the two views in sync.
    transcript_genome_spans: OnceCell<HashMap<String, (u64, u64)>>,
    /// Deferred secondary builds: build name -> source cdot path. Set by
    /// `MultiFastaProvider::from_manifest` so a declared secondary (e.g. the
    /// GRCh37 cdot) is not loaded until a lookup for that build needs it. NOT
    /// part of any rkyv/bincode snapshot — purely a runtime handle.
    deferred_alt_sources: HashMap<String, PathBuf>,
    /// Lazily-loaded secondary-build mappers, keyed by build name. `get_or_init`
    /// loads `deferred_alt_sources[build]` on first use; the inner `Option` is
    /// `None` when the source is missing/unreadable (best-effort, never retried).
    lazy_alt_mappers: HashMap<String, OnceCell<Option<Box<CdotMapper>>>>,
}

/// Per-contig stab-query index entry for a non-primary build, used by
/// [`CdotMapper::transcripts_at_position`] to route alt-build contig queries.
///
/// Holds the owning build name, the accession list for the contig in that
/// build's [`CdotMapper::alt_build_transcripts`] map, and a SuperIntervals
/// index whose `u32` payload indexes into `accessions`.
#[derive(Debug, Clone)]
struct AltBuildContigStab {
    /// Genome build that owns this contig (key into `alt_build_transcripts`).
    build: String,
    /// Accessions on this contig, in the order the `IntervalMap` payloads index.
    accessions: Vec<String>,
    /// Stab-query index; payload is an index into `accessions`.
    index: IntervalMap<u32>,
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
            alt_build_query_index: OnceCell::new(),
            transcript_genome_spans: OnceCell::new(),
            deferred_alt_sources: HashMap::new(),
            lazy_alt_mappers: HashMap::new(),
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
            alt_build_query_index: OnceCell::new(),
            transcript_genome_spans: OnceCell::new(),
            deferred_alt_sources: HashMap::new(),
            lazy_alt_mappers: HashMap::new(),
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

    /// `mmap`-backed load from a pre-serialized rkyv archive.
    ///
    /// Maps the archive file into memory (one `unsafe` `memmap2::Mmap::map`
    /// call) then validates and deserializes it in a single pass, rejects stale
    /// or corrupt archives so [`load`](Self::load) can regenerate them, and is
    /// ~3-4x faster than the bincode path because rkyv's layout skips the
    /// byte-stream decode.
    pub fn from_rkyv_file<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        use rkyv::rancor::Error as RkyvError;

        let file = std::fs::File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open cdot rkyv file: {}", e),
        })?;
        // mmap (page-aligned, satisfies rkyv's alignment) so the 260MB+ archive
        // isn't read+copied up front — pages fault in lazily as `deserialize`
        // touches them. The mapping is dropped once we own the data below.
        // SAFETY: concurrent writers use `rename` (new inode) so an active
        // mapping is never invalidated by a parallel cache refresh; external
        // truncation of the file is not guarded against.
        let mmap = unsafe {
            memmap2::Mmap::map(&file).map_err(|e| FerroError::Io {
                msg: format!("Failed to mmap cdot rkyv file: {}", e),
            })?
        };

        let archived = rkyv::access::<rkyv_cache::ArchivedRkyvSnapshot, RkyvError>(&mmap[..])
            .map_err(|e| FerroError::Io {
                msg: format!("cdot rkyv archive failed validation: {}", e),
            })?;
        let version = rkyv_cache::version_of(archived);
        if version != rkyv_cache::RKYV_FORMAT_VERSION {
            return Err(FerroError::Io {
                msg: format!(
                    "cdot rkyv schema version {} != expected {}",
                    version,
                    rkyv_cache::RKYV_FORMAT_VERSION
                ),
            });
        }

        // Materialize owned cdot maps directly from the archive (single pass).
        // Rejects (rather than coerces) any out-of-range enum tag, so a corrupt
        // archive surfaces as an error and is regenerated by the caller.
        let maps = rkyv_cache::maps_from_archived(archived)?;
        Ok(Self {
            transcripts: maps.transcripts,
            contig_index: maps.contig_index,
            contig_alias_to_canonical: maps.contig_alias_to_canonical,
            base_to_versioned: maps.base_to_versioned,
            lrg_to_refseq: maps.lrg_to_refseq,
            alt_build_transcripts: maps.alt_build_transcripts,
            primary_build: maps.primary_build,
            contig_query_index: OnceCell::new(),
            alt_build_query_index: OnceCell::new(),
            transcript_genome_spans: OnceCell::new(),
            deferred_alt_sources: HashMap::new(),
            lazy_alt_mappers: HashMap::new(),
        })
    }

    /// Cheap check: does `path` begin with a current, structurally valid rkyv
    /// archive? Reads and validates (the archive must be parsed to reach the
    /// version field), so it is not free — but far cheaper than a full load,
    /// and lets `prepare` decide whether to regenerate.
    pub fn rkyv_is_current<P: AsRef<Path>>(path: P) -> bool {
        let Ok(file) = std::fs::File::open(path.as_ref()) else {
            return false;
        };
        // SAFETY: concurrent writers use `rename` (new inode) so an active
        // mapping is never invalidated by a parallel cache refresh; external
        // truncation of the file is not guarded against.
        let Ok(mmap) = (unsafe { memmap2::Mmap::map(&file) }) else {
            return false;
        };
        match rkyv::access::<rkyv_cache::ArchivedRkyvSnapshot, rkyv::rancor::Error>(&mmap[..]) {
            Ok(a) => a.format_version.to_native() == rkyv_cache::RKYV_FORMAT_VERSION,
            Err(_) => false,
        }
    }

    /// Write the mapper to an rkyv archive (atomically), for fast loading.
    pub fn to_rkyv_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FerroError> {
        let snapshot = rkyv_cache::snapshot_from_mapper(self);
        let bytes =
            rkyv::to_bytes::<rkyv::rancor::Error>(&snapshot).map_err(|e| FerroError::Io {
                msg: format!("Failed to serialize cdot rkyv: {}", e),
            })?;

        // Unique temp sibling (pid + process-global counter) so concurrent
        // writers in the same process never collide — see `to_bincode_file`.
        use std::sync::atomic::{AtomicU64, Ordering};
        static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);
        let path = path.as_ref();
        let mut tmp = path.to_path_buf();
        let tmp_name = format!(
            "{}.tmp.{}.{}",
            path.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "cdot.rkyv".to_string()),
            std::process::id(),
            TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        );
        tmp.set_file_name(tmp_name);

        if let Err(e) = std::fs::write(&tmp, &bytes) {
            let _ = std::fs::remove_file(&tmp);
            return Err(FerroError::Io {
                msg: format!("Failed to write cdot rkyv temp file: {}", e),
            });
        }
        std::fs::rename(&tmp, path).map_err(|e| {
            let _ = std::fs::remove_file(&tmp);
            FerroError::Io {
                msg: format!("Failed to finalize cdot rkyv file: {}", e),
            }
        })
    }

    /// Load a cdot file using the fastest available cache format.
    ///
    /// Given a path to a `.json` or `.json.gz` file, checks for sibling cache
    /// files in order: `.rkyv` (fastest, preferred) then `.bin` (legacy
    /// bincode). If only `.bin` exists, loads from it and promotes it to
    /// `.rkyv` for future runs. Falls back to JSON if no cache is present or
    /// usable, then writes a fresh `.rkyv` cache. Also accepts direct `.rkyv`
    /// or `.bin` paths.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        Self::load_with_source(path).map(|(mapper, _)| mapper)
    }

    /// Like [`load`](Self::load), but also reports which path produced the
    /// mapper. This closes the API hole that let the #585 regression go silent:
    /// before this, a caller could not tell a fast archive load from a slow
    /// JSON re-parse, so a silently-broken cache (re-parsing 512 MB of JSON on
    /// every run) was invisible. The returned [`CdotLoadSource`] makes that
    /// distinction assertable in tests and observable via `ferro check`.
    ///
    /// A binary cache (`.rkyv` or legacy `.bin`) counts as
    /// [`CdotLoadSource::Archive`] — the "fast path". Only an actual JSON parse
    /// returns [`CdotLoadSource::JsonFallback`]. Note that the JSON-fallback arm
    /// still self-heals the sibling `.rkyv` as a side effect, so a *subsequent*
    /// `load_with_source` on the same path reports `Archive`.
    pub fn load_with_source<P: AsRef<Path>>(path: P) -> Result<(Self, CdotLoadSource), FerroError> {
        let path = path.as_ref();
        let path_str = path.to_string_lossy();

        // Direct archive paths — no JSON fallback available.
        if path_str.ends_with(".rkyv") {
            return Self::from_rkyv_file(path).map(|m| (m, CdotLoadSource::Archive));
        }
        if path_str.ends_with(".bin") {
            return Self::from_bincode_file(path).map(|m| (m, CdotLoadSource::Archive));
        }

        // Sibling cache paths next to the JSON. Replace the *single* trailing
        // extension (`.json`, or `.json` after stripping `.gz`); the accession
        // stem itself contains dots (e.g. `cdot-0.2.32.refseq.GRCh38`) so we
        // must not over-strip.
        let (rkyv_path, bin_path) = if path_str.ends_with(".json.gz") {
            let stem = path.with_extension(""); // foo.json.gz -> foo.json
            (stem.with_extension("rkyv"), stem.with_extension("bin"))
        } else {
            (path.with_extension("rkyv"), path.with_extension("bin"))
        };

        // A sibling cache is only usable if it is at least as new as the source
        // JSON. cdot files are normally version-stamped in the filename, but a
        // re-`prepare` that rewrites the JSON under the same name would
        // otherwise be masked by a stale archive. Mirrors the freshness check
        // the GFF/GTF loader applies via `TranscriptDb::cache_is_fresh`. An
        // unreadable mtime (cache or source) is treated as "not fresh" so we
        // fall through to JSON rather than serve a stale cache.
        let cache_is_fresh = |cache: &Path| -> bool {
            let Ok(cache_mtime) = std::fs::metadata(cache).and_then(|m| m.modified()) else {
                return false;
            };
            let Ok(source_mtime) = std::fs::metadata(path).and_then(|m| m.modified()) else {
                return false;
            };
            cache_mtime >= source_mtime
        };

        // 1. Prefer the rkyv archive (fastest, validated).
        if rkyv_path.exists() && cache_is_fresh(&rkyv_path) {
            match Self::from_rkyv_file(&rkyv_path) {
                Ok(mapper) => return Ok((mapper, CdotLoadSource::Archive)),
                Err(e) => eprintln!(
                    "warning: cdot rkyv cache {} is unusable ({}); trying other formats",
                    rkyv_path.display(),
                    e
                ),
            }
        }

        // 2. Legacy bincode cache — usable, but upgrade to rkyv for next time.
        //    Only promote a *fresh* bincode cache to rkyv; a stale one is
        //    skipped here and the JSON below refreshes both.
        if bin_path.exists() && cache_is_fresh(&bin_path) {
            match Self::from_bincode_file(&bin_path) {
                Ok(mapper) => {
                    if let Err(e) = mapper.to_rkyv_file(&rkyv_path) {
                        eprintln!(
                            "note: could not write cdot rkyv cache {}: {} (continuing)",
                            rkyv_path.display(),
                            e
                        );
                    }
                    return Ok((mapper, CdotLoadSource::Archive));
                }
                Err(e) => eprintln!(
                    "warning: cdot bincode cache {} is unusable ({}); \
                     parsing JSON (much slower) and refreshing the cache",
                    bin_path.display(),
                    e
                ),
            }
        }

        // 3. Fall back to JSON.
        let mapper = if path_str.ends_with(".gz") {
            Self::from_json_gz(path)?
        } else {
            Self::from_json_file(path)?
        };

        // Self-heal: best-effort refresh of the sibling rkyv cache so the next
        // load takes the fast path. Failures (e.g. a read-only mount) are
        // non-fatal — we already have a usable mapper.
        if let Err(e) = mapper.to_rkyv_file(&rkyv_path) {
            eprintln!(
                "note: could not refresh cdot rkyv cache {}: {} (continuing)",
                rkyv_path.display(),
                e
            );
        }

        Ok((mapper, CdotLoadSource::JsonFallback))
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

    /// Load a second cdot file as a NON-primary build, layering its
    /// transcripts onto this mapper without disturbing the primary build.
    ///
    /// The real reference manifest ships GRCh37 and GRCh38 cdot files
    /// separately (single-build each), so after loading the primary (GRCh38)
    /// cdot we layer the GRCh37 file here. The loaded transcripts go into
    /// [`Self::alt_build_transcripts`]`[build_name]` — they are NOT merged into
    /// the primary [`Self::transcripts`] map, because the same accession exists
    /// on both builds with different coordinates and a merge would corrupt the
    /// primary build (collision safety).
    ///
    /// Honors the same archive/JSON path handling as [`Self::load`] (`.rkyv`,
    /// `.bin`, `.json`, `.json.gz`). Returns the number of transcripts loaded
    /// into the secondary build.
    pub fn load_secondary_build<P: AsRef<Path>>(
        &mut self,
        path: P,
        build_name: &str,
    ) -> Result<usize, FerroError> {
        // Reuse the full `load` path so caches and `.gz`/`.rkyv`/`.bin` are
        // handled identically to the primary load. The loaded mapper's
        // `transcripts` map holds the secondary build's data (its primary
        // build); we then absorb it as our alt build.
        let other = Self::load(path.as_ref())?;
        Ok(self.absorb_secondary_build(other, build_name))
    }

    /// Reader-based variant of [`Self::load_secondary_build`] for tests and
    /// in-memory callers. Parses a single-build cdot JSON from `reader` and
    /// layers it as `build_name`.
    pub fn load_secondary_build_from_reader<R: Read>(
        &mut self,
        reader: R,
        build_name: &str,
    ) -> Result<usize, FerroError> {
        let other = Self::from_reader_with_build(reader, build_name)?;
        Ok(self.absorb_secondary_build(other, build_name))
    }

    /// Absorb the transcripts of an independently-loaded `other` mapper as the
    /// non-primary build `build_name`. The `other` mapper's primary
    /// `transcripts` become this mapper's `alt_build_transcripts[build_name]`;
    /// any builds `other` itself nested under `genome_builds` are also merged
    /// in (so a multi-build secondary file still contributes its nested views).
    /// Clears the lazy overlap indexes so they rebuild with the new contigs.
    fn absorb_secondary_build(&mut self, other: Self, build_name: &str) -> usize {
        // Snapshot the build's transcript count so we can report how many this
        // secondary file *newly* contributes (the primary cdot may already have
        // populated some `build_name` views from its own `genome_builds`).
        let before = self
            .alt_build_transcripts
            .get(build_name)
            .map_or(0, |m| m.len());
        let CdotMapper {
            transcripts: other_transcripts,
            alt_build_transcripts: other_alt,
            lrg_to_refseq: other_lrg,
            contig_alias_to_canonical: other_aliases,
            ..
        } = other;

        // The requested build's transcripts can arrive two ways: as `other`'s
        // primary view (when the secondary file's primary build *is*
        // `build_name`), or nested under `other.alt_build_transcripts[build_name]`
        // (when the file nests its data under `genome_builds`, which the real
        // RefSeq cdot files do — there the primary view selected by the default
        // load build is empty). Absorb both into our `build_name` slot.
        {
            let target = self
                .alt_build_transcripts
                .entry(build_name.to_string())
                .or_default();
            for (acc, tx) in other_transcripts {
                target.insert(acc, tx);
            }
        }
        // Fold in any builds the secondary file nested under `genome_builds`,
        // EXCEPT the primary build of this mapper (that data is authoritative
        // here and must not be overwritten by a secondary file's view).
        for (nested_build, txs) in other_alt {
            if Some(nested_build.as_str()) == self.primary_build.as_deref() {
                continue;
            }
            let dest = self.alt_build_transcripts.entry(nested_build).or_default();
            for (acc, tx) in txs {
                dest.insert(acc, tx);
            }
        }

        // Carry over the secondary build's contig aliases (e.g. chr17 ↔
        // NC_000017.10) so alt-build contig lookups can resolve UCSC/Ensembl
        // names too. Do not overwrite an existing primary-build alias.
        for (alias, canonical) in other_aliases {
            self.contig_alias_to_canonical
                .entry(alias)
                .or_insert(canonical);
        }
        // Merge LRG mappings the secondary file may carry; keep existing ones.
        for (lrg, refseq) in other_lrg {
            self.lrg_to_refseq.entry(lrg).or_insert(refseq);
        }

        // The new alt-build contigs invalidate the lazy overlap indexes.
        self.alt_build_query_index = OnceCell::new();
        // Report how many transcripts this secondary file newly contributed to
        // `build_name` (the net new entries). For real (`genome_builds`-nested)
        // cdot files the default-build primary view is empty, so the old count —
        // which only tallied `other.transcripts` — logged a misleading "0" even
        // though the build's transcripts were absorbed via the fold above.
        self.alt_build_transcripts
            .get(build_name)
            .map_or(0, |m| m.len())
            .saturating_sub(before)
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

        // Update base to versioned index (e.g., "NM_000088" -> "NM_000088.4").
        //
        // The cdot file deserializes its transcripts into a `HashMap`, so the
        // order `add_transcript` is called in is non-deterministic across runs.
        // A plain last-writer-wins insert would therefore map a base accession
        // to a *random* one of its versions, which surfaces downstream as a
        // run-to-run flip in version-fallback resolution (e.g. requesting
        // `NM_003002.2` — absent from cdot — would resolve to `.3` or `.4`
        // arbitrarily, shifting CDS coordinates by the versions' length delta).
        // Pick the highest version deterministically so the fallback is stable.
        if let Some(base) = accession.split('.').next() {
            let replace = match self.base_to_versioned.get(base) {
                None => true,
                Some(existing) => {
                    (accession_version(&accession), accession.as_str())
                        > (accession_version(existing), existing.as_str())
                }
            };
            if replace {
                self.base_to_versioned
                    .insert(base.to_string(), accession.clone());
            }
        }

        self.transcripts.insert(accession, transcript);

        // Any previously-built query index is stale now. The next query will
        // rebuild it from the updated `contig_index` + `transcripts`. The
        // alt-build index is derived from `alt_build_transcripts`, which
        // `add_transcript` does not touch, but clear it too so any consumer
        // that mixes `add_transcript` with alt-build mutation stays consistent.
        self.contig_query_index = OnceCell::new();
        self.alt_build_query_index = OnceCell::new();
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
    /// section). If the transcript is not found in the eagerly-captured alt
    /// views, falls back to a lazily-loaded deferred secondary build registered
    /// via [`defer_secondary_build`](Self::defer_secondary_build), if any.
    /// Returns `None` if the transcript is absent from the requested build (or
    /// if neither the alt-build cache nor any deferred secondary covers it,
    /// e.g. for bincode-loaded mappers or `from_transcripts`-built mappers).
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
        // Eagerly-captured alt views for this build, if any. Absent when the
        // build is deferred (consulted below).
        if let Some(alt) = self.alt_build_transcripts.get(build) {
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
        }
        // Not in the eagerly-captured alt views — consult a deferred secondary
        // build (e.g. a manifest's GRCh37 cdot) loaded lazily on first use.
        self.deferred_alt_mapper(build)?
            .get_transcript_on_build(accession, build)
    }

    /// Record a secondary build to load lazily on first use, instead of folding
    /// it in eagerly. `path` is any cdot source `load()` accepts (`.json`,
    /// `.json.gz`, `.rkyv`, `.bin`). The build's transcripts are not read until
    /// `deferred_alt_mapper(build)` is first called.
    ///
    /// Calling this again after the `OnceCell` for `build` has already been
    /// initialized updates the recorded path but has no effect on the
    /// already-loaded mapper — `OnceCell` is init-once and the existing value
    /// is returned on all subsequent calls.
    pub fn defer_secondary_build(&mut self, build: &str, path: PathBuf) {
        self.deferred_alt_sources.insert(build.to_string(), path);
        self.lazy_alt_mappers.entry(build.to_string()).or_default();
        // A newly-deferred build adds alt-build contigs; invalidate the lazy
        // overlap index so a previously force-built index doesn't omit them
        // (mirrors the eager `load_secondary_build` reset).
        self.alt_build_query_index = OnceCell::new();
    }

    /// Return the lazily-loaded secondary mapper for `build`, loading it from its
    /// deferred source on first call. `None` if `build` was never deferred or its
    /// source could not be loaded (best-effort; resolved once and not retried).
    fn deferred_alt_mapper(&self, build: &str) -> Option<&CdotMapper> {
        let cell = self.lazy_alt_mappers.get(build)?;
        cell.get_or_init(|| {
            let path = self.deferred_alt_sources.get(build)?;
            match Self::load(path) {
                Ok(m) => {
                    let n = m.alt_build_transcripts.get(build).map_or(0, |t| t.len());
                    eprintln!("Loaded {n} {build} transcripts (deferred secondary build)");
                    Some(Box::new(m))
                }
                Err(e) => {
                    eprintln!(
                        "warning: failed to load deferred {build} cdot {}: {e}",
                        path.display()
                    );
                    None
                }
            }
        })
        .as_deref()
    }

    /// Test probe: has the deferred-load attempt for `build` completed (whether it loaded or failed)?
    #[cfg(test)]
    pub(crate) fn deferred_alt_loaded(&self, build: &str) -> bool {
        self.lazy_alt_mappers
            .get(build)
            .is_some_and(|c| c.get().is_some())
    }

    /// List every genome build that has data for the given transcript.
    ///
    /// Includes the primary build if the transcript is present there. The
    /// returned list is unordered; callers that want a deterministic probe
    /// order (e.g. GRCh38-then-GRCh37) should sort externally.
    ///
    /// Side effect: this can force-load a deferred secondary build on first
    /// use (via `deferred_alt_mapper`), which may read a large cdot source
    /// (~198k transcripts) from disk. Callers on a hot path should be aware
    /// this is not a cheap metadata-only query.
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
        // Also consider deferred secondary builds (loaded on demand). Only report
        // a deferred build when it actually has the transcript, mirroring the
        // eager check above.
        let deferred: Vec<String> = self.lazy_alt_mappers.keys().cloned().collect();
        for build in deferred {
            if builds.iter().any(|b| b == &build) {
                continue;
            }
            // `get_transcript_on_build` applies its own version-fallback, so we
            // don't repeat the eager loop's inline fallback here.
            if let Some(sub) = self.deferred_alt_mapper(&build) {
                if sub.get_transcript_on_build(accession, &build).is_some() {
                    builds.push(build);
                }
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
        // i32 fits every chromosome in any current build (max ~250M < 2.1B).
        // Positions outside i32 cannot overlap a real transcript span anyway.
        let Ok(p) = i32::try_from(pos) else {
            return Vec::new();
        };

        // Primary build: the contig lives in `contig_index` / `transcripts`.
        let by_contig = self
            .contig_query_index
            .get_or_init(|| self.build_query_index());
        if let (Some(im), Some(accessions)) =
            (by_contig.get(canonical), self.contig_index.get(canonical))
        {
            let mut hits: Vec<u32> = Vec::with_capacity(16);
            im.search_stabbed(p, &mut hits);
            return hits
                .into_iter()
                .filter_map(|idx| {
                    let acc = accessions.get(idx as usize)?;
                    let tx = self.transcripts.get(acc)?;
                    Some((acc.as_str(), tx))
                })
                .collect();
        }

        // Non-primary build: GRCh37/GRCh38 contig accessions are distinct, so
        // the contig itself names the owning build. Route the query into that
        // build's `alt_build_transcripts` map.
        //
        // The index is built lazily and includes contigs from both eagerly-loaded
        // secondary builds (`alt_build_transcripts`) and deferred secondary builds
        // (force-loaded from `lazy_alt_mappers` at index-build time). Transcripts
        // for a given build may be split across both stores (the primary cdot's
        // eager alt-build view and the deferred secondary file), so each hit
        // accession is resolved from both stores before being dropped as missing.
        let alt_by_contig = self
            .alt_build_query_index
            .get_or_init(|| self.build_alt_build_query_index());
        let Some(stab) = alt_by_contig.get(canonical) else {
            return Vec::new();
        };
        let eager_alt = self.alt_build_transcripts.get(&stab.build);
        // Resolve the deferred sub-mapper now (force-loading if needed) so we
        // can check both stores per accession below. If neither store has data
        // for this build the whole query returns empty.
        let deferred_alt = self
            .deferred_alt_mapper(&stab.build)
            .and_then(|sub| sub.alt_build_transcripts.get(&stab.build));
        if eager_alt.is_none() && deferred_alt.is_none() {
            return Vec::new();
        }
        let mut hits: Vec<u32> = Vec::with_capacity(16);
        stab.index.search_stabbed(p, &mut hits);
        hits.into_iter()
            .filter_map(|idx| {
                let acc = stab.accessions.get(idx as usize)?;
                // Check the eager store first; fall back to the deferred store.
                // A transcript can only live in one of the two stores.
                let tx = eager_alt
                    .and_then(|m| m.get(acc))
                    .or_else(|| deferred_alt.and_then(|m| m.get(acc)))?;
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

    /// Construct the per-contig stab-query index for the non-primary builds,
    /// derived from [`Self::alt_build_transcripts`] and any deferred secondary
    /// builds registered via [`Self::defer_secondary_build`] (which are
    /// force-loaded here). Mirrors [`Self::build_query_index`] but groups by
    /// `(build, contig)` so a query on an alt-build contig can recover both the
    /// owning build and the accession.
    ///
    /// Because GRCh37 and GRCh38 contig accessions are distinct, no contig is
    /// claimed by two builds; the per-contig map is therefore unambiguous.
    ///
    /// **Merge semantics**: a single build (e.g. "GRCh37") can have transcripts
    /// in *both* the eager store (`self.alt_build_transcripts`) and in a deferred
    /// secondary mapper, potentially on the *same* contig.  The old two-pass
    /// approach (call `index_one_build` once for eager, once for deferred) used
    /// `HashMap::insert` which silently replaced the eager contig stab with the
    /// deferred one, discarding all eagerly-indexed accessions.  The fix below
    /// gathers all `(build, accession)` pairs across **both** stores first, then
    /// builds exactly one `AltBuildContigStab` per contig from the combined,
    /// sorted accession list.
    fn build_alt_build_query_index(&self) -> HashMap<String, AltBuildContigStab> {
        // Step 1 — collect (build, Vec<accession>) per contig from every source.
        // A contig maps to exactly one build — NC_*.10 = GRCh37, NC_*.11 =
        // GRCh38 — so first-seen wins; the debug_assert guards the invariant
        // in test builds.
        let mut contig_build: HashMap<String, String> = HashMap::new();
        let mut contig_accs: HashMap<String, Vec<String>> = HashMap::new();

        let mut add_transcripts = |build: &str, txs: &HashMap<String, CdotTranscript>| {
            for (acc, tx) in txs {
                let contig = &tx.contig;
                let existing = contig_build
                    .entry(contig.clone())
                    .or_insert_with(|| build.to_string());
                debug_assert_eq!(
                    existing.as_str(),
                    build,
                    "contig {contig} claimed by two builds"
                );
                contig_accs
                    .entry(contig.clone())
                    .or_default()
                    .push(acc.clone());
            }
        };

        // Eager store.
        for (build, transcripts) in &self.alt_build_transcripts {
            add_transcripts(build, transcripts);
        }

        // Deferred secondary builds. Collect the keys first to avoid holding a
        // borrow on `self.lazy_alt_mappers` while calling `deferred_alt_mapper`
        // (which also borrows `self`).
        let deferred: Vec<String> = self.lazy_alt_mappers.keys().cloned().collect();
        for build in &deferred {
            if let Some(sub) = self.deferred_alt_mapper(build) {
                if let Some(transcripts) = sub.alt_build_transcripts.get(build.as_str()) {
                    add_transcripts(build, transcripts);
                }
            }
        }

        // Step 2 — for every contig, sort the combined accession list and build
        // one IntervalMap.  To resolve a transcript's exons we need to look the
        // accession up in a transcript map.  A given accession may live in either
        // the eager store or the deferred sub-mapper; try both.
        let mut by_contig: HashMap<String, AltBuildContigStab> =
            HashMap::with_capacity(contig_accs.len());
        for (contig, mut accessions) in contig_accs {
            // Deterministic payload→accession mapping (HashMap iteration order
            // is random across runs; sort to stabilise the u32 index).
            accessions.sort();
            // An accession present in BOTH the eager alt-build view and the deferred sub-mapper
            // would appear twice in the combined list; dedup keeps the IntervalMap
            // payload->accession mapping 1:1.
            accessions.dedup();

            let build = contig_build
                .get(&contig)
                .expect("contig_build populated in same loop")
                .as_str();

            let eager_map = self.alt_build_transcripts.get(build);
            // Resolve the deferred sub-mapper (force-loading if needed).
            let deferred_map = self
                .deferred_alt_mapper(build)
                .and_then(|sub| sub.alt_build_transcripts.get(build));

            let mut im: IntervalMap<u32> = IntervalMap::new();
            for (idx, acc) in accessions.iter().enumerate() {
                // Look up the transcript in the eager store first, then deferred.
                let tx = eager_map
                    .and_then(|m| m.get(acc))
                    .or_else(|| deferred_map.and_then(|m| m.get(acc)));
                let Some(tx) = tx else {
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
            by_contig.insert(
                contig,
                AltBuildContigStab {
                    build: build.to_string(),
                    accessions,
                    index: im,
                },
            );
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
    fn test_concurrent_to_rkyv_file_same_path() {
        // Multiple threads writing the same rkyv cache path concurrently must
        // not collide on their temp siblings or leave a half-written archive:
        // the final cache must always be loadable and complete.  Mirrors the
        // bincode variant above; `to_rkyv_file` uses the same pid-counter-rename
        // pattern.
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
        let rkyv_path = std::sync::Arc::new(temp_dir.path().join("concurrent.rkyv"));

        let handles: Vec<_> = (0..8)
            .map(|_| {
                let original = std::sync::Arc::clone(&original);
                let rkyv_path = std::sync::Arc::clone(&rkyv_path);
                std::thread::spawn(move || original.to_rkyv_file(rkyv_path.as_path()))
            })
            .collect();
        for handle in handles {
            handle
                .join()
                .unwrap()
                .expect("concurrent rkyv write should succeed");
        }

        // The final archive is complete and loadable, and no temp siblings leaked.
        assert!(
            CdotMapper::rkyv_is_current(rkyv_path.as_path()),
            "rkyv cache must be current after concurrent writes"
        );
        let loaded = CdotMapper::from_rkyv_file(rkyv_path.as_path()).unwrap();
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
    fn test_load_prefers_rkyv_over_bincode() {
        // When both a .rkyv sibling and a .bin sibling exist, load() must pick
        // the .rkyv one. Each on-disk representation carries a unique sentinel in
        // gene_name so the returned value identifies which path was taken.
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

        // Build distinct mappers for each on-disk format.
        let mut from_bin = original.clone();
        from_bin
            .transcripts
            .get_mut("NM_000088.3")
            .unwrap()
            .gene_name = Some("FROM_BIN".to_string());

        let mut from_rkyv = original.clone();
        from_rkyv
            .transcripts
            .get_mut("NM_000088.3")
            .unwrap()
            .gene_name = Some("FROM_RKYV".to_string());

        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");
        let bin_path = temp_dir.path().join("cdot.bin");
        let rkyv_path = temp_dir.path().join("cdot.rkyv");

        std::fs::write(&json_path, json).unwrap();
        from_bin.to_bincode_file(&bin_path).unwrap();
        from_rkyv.to_rkyv_file(&rkyv_path).unwrap();

        // With all three siblings present, load() must prefer .rkyv.
        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(
            loaded
                .get_transcript("NM_000088.3")
                .unwrap()
                .gene_name
                .as_deref(),
            Some("FROM_RKYV"),
            "load() should prefer .rkyv over .bin and JSON"
        );
    }

    #[test]
    fn test_load_skips_stale_rkyv_and_reparses_json() {
        // A sibling .rkyv OLDER than the source JSON must not be served: a
        // re-`prepare` that rewrites the JSON under the same name would
        // otherwise be masked by the stale archive. load() must skip it, reparse
        // the JSON, and refresh the cache.
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

        // A stale rkyv carrying a sentinel gene_name distinct from the JSON.
        let mut stale = original.clone();
        stale.transcripts.get_mut("NM_000088.3").unwrap().gene_name =
            Some("FROM_STALE_RKYV".to_string());

        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");
        let rkyv_path = temp_dir.path().join("cdot.rkyv");

        std::fs::write(&json_path, json).unwrap();
        stale.to_rkyv_file(&rkyv_path).unwrap();

        // Backdate the rkyv so it is strictly older than the JSON source.
        let backdated = std::time::SystemTime::now() - std::time::Duration::from_secs(60);
        std::fs::File::options()
            .write(true)
            .open(&rkyv_path)
            .unwrap()
            .set_modified(backdated)
            .unwrap();

        // load() must ignore the stale rkyv and reparse the JSON (COL1A1).
        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(
            loaded
                .get_transcript("NM_000088.3")
                .unwrap()
                .gene_name
                .as_deref(),
            Some("COL1A1"),
            "load() must skip an rkyv cache older than its source JSON"
        );

        // Self-heal: the stale archive should have been refreshed from JSON, so
        // reading the rkyv directly now yields the JSON's value, not the stale
        // sentinel.
        let refreshed = CdotMapper::from_rkyv_file(&rkyv_path).unwrap();
        assert_eq!(
            refreshed
                .get_transcript("NM_000088.3")
                .unwrap()
                .gene_name
                .as_deref(),
            Some("COL1A1"),
            "load() should refresh the stale rkyv after reparsing the JSON"
        );
    }

    // -------------------------------------------------------------------------
    // Load-source observability (the #585 perf-regression gate). These assert
    // the `CdotLoadSource` tag, not just the data — so a silent fast-path →
    // JSON-fallback regression fails the build deterministically, with no
    // reference data and no timing thresholds.
    // -------------------------------------------------------------------------

    /// A fresh sibling archive is taken AND reported as `Archive`, and the data
    /// round-trips through it. Guards "layout drifted but load still claims the
    /// fast path" (#585): a load that quietly fell back to JSON would report
    /// `JsonFallback` and fail this.
    #[test]
    fn load_with_source_reports_archive_for_fresh_rkyv() {
        let json = multi_build_cdot_json();
        let expected = CdotMapper::from_reader(json.as_bytes()).unwrap();

        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");
        let rkyv_path = temp_dir.path().join("cdot.rkyv");

        // Write the JSON first, then the sibling archive — so the archive is at
        // least as new as the source and `load` prefers it.
        std::fs::write(&json_path, json).unwrap();
        expected.to_rkyv_file(&rkyv_path).unwrap();

        let (loaded, source) = CdotMapper::load_with_source(&json_path).unwrap();
        assert_eq!(
            source,
            CdotLoadSource::Archive,
            "a fresh sibling .rkyv must load via the fast archive path"
        );
        // Data round-trips through the archive.
        assert_eq!(loaded.transcript_count(), expected.transcript_count());
        assert_eq!(loaded.primary_build, expected.primary_build);
    }

    /// A stale archive is observably *not* served: the load falls back to JSON
    /// and reports `JsonFallback` (rather than silently re-parsing forever as in
    /// #585), and the side-effecting self-heal makes the *next* load report
    /// `Archive`.
    #[test]
    fn load_with_source_reports_json_fallback_then_self_heals() {
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

        // A stale rkyv carrying a sentinel distinct from the JSON.
        let mut stale = original.clone();
        stale.transcripts.get_mut("NM_000088.3").unwrap().gene_name =
            Some("FROM_STALE_RKYV".to_string());

        let temp_dir = tempfile::TempDir::new().unwrap();
        let json_path = temp_dir.path().join("cdot.json");
        let rkyv_path = temp_dir.path().join("cdot.rkyv");

        std::fs::write(&json_path, json).unwrap();
        stale.to_rkyv_file(&rkyv_path).unwrap();

        // Backdate the rkyv so it is strictly older than the JSON source.
        let backdated = std::time::SystemTime::now() - std::time::Duration::from_secs(60);
        std::fs::File::options()
            .write(true)
            .open(&rkyv_path)
            .unwrap()
            .set_modified(backdated)
            .unwrap();

        // First load: stale archive skipped → slow JSON path, observably tagged.
        let (first, first_source) = CdotMapper::load_with_source(&json_path).unwrap();
        assert_eq!(
            first_source,
            CdotLoadSource::JsonFallback,
            "a stale archive must surface as a JSON fallback, not be silently served"
        );
        assert_eq!(
            first
                .get_transcript("NM_000088.3")
                .unwrap()
                .gene_name
                .as_deref(),
            Some("COL1A1"),
            "the fallback must serve the real JSON data, not the stale sentinel"
        );

        // Second load: the self-heal refreshed the archive → fast path restored.
        let (_second, second_source) = CdotMapper::load_with_source(&json_path).unwrap();
        assert_eq!(
            second_source,
            CdotLoadSource::Archive,
            "after self-heal the next load must take the fast archive path \
             (not re-parse JSON forever, the #585 failure)"
        );
    }

    #[test]
    fn test_load_bin_upgrades_to_rkyv() {
        // When only a .bin sibling is present (no .rkyv), load() must load from
        // .bin and then write a .rkyv sibling for future runs.
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
        let json_path = temp_dir.path().join("cdot.json");
        let bin_path = temp_dir.path().join("cdot.bin");
        let rkyv_path = temp_dir.path().join("cdot.rkyv");

        std::fs::write(&json_path, json).unwrap();
        original.to_bincode_file(&bin_path).unwrap();

        // No .rkyv exists yet.
        assert!(!rkyv_path.exists(), ".rkyv must not exist before load()");

        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
        assert!(loaded.get_transcript("NM_000088.3").is_some());

        // load() must have written the .rkyv promotion.
        assert!(
            CdotMapper::rkyv_is_current(&rkyv_path),
            "load() should promote .bin to a valid .rkyv sibling"
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
        let rkyv_path = temp_dir.path().join("cdot.rkyv");

        // Write valid JSON and a corrupt bincode cache; no rkyv cache yet.
        std::fs::write(&json_path, json).unwrap();
        std::fs::write(&bin_path, b"not valid bincode data").unwrap();

        // load() should fall back to JSON despite the corrupt .bin.
        let loaded = CdotMapper::load(&json_path).unwrap();
        assert_eq!(loaded.transcript_count(), 1);
        assert!(loaded.get_transcript("NM_000088.3").is_some());

        // Self-heal now writes the fast rkyv cache, so a second load takes the
        // fast path.
        assert!(
            CdotMapper::rkyv_is_current(&rkyv_path),
            "an rkyv cache should be written after the JSON fallback"
        );
        let reloaded = CdotMapper::from_rkyv_file(&rkyv_path).unwrap();
        assert_eq!(reloaded.transcript_count(), 1);
        assert!(reloaded.get_transcript("NM_000088.3").is_some());
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
    fn test_rkyv_roundtrip_preserves_data() {
        // Same transcript as the bincode test, plus a per-exon CIGAR ("M73"),
        // so the rkyv path is verified to preserve exons, strand, CDS, the
        // per-exon cigar, AND the alt-build (GRCh37) map.
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000017.10",
                            "strand": "-",
                            "exons": [[48263025, 48263098, 1, 0, 73, "M73"]]
                        },
                        "GRCh38": {
                            "contig": "NC_000017.11",
                            "strand": "-",
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
        let p = temp.path().join("cdot.rkyv");
        from_json.to_rkyv_file(&p).unwrap();
        assert!(CdotMapper::rkyv_is_current(&p));
        let from_rkyv = CdotMapper::from_rkyv_file(&p).unwrap();

        // Primary-build transcript must be byte-identical across the roundtrip,
        // including strand, exons, CDS, and the per-exon cigars.
        let a = from_json.get_transcript("NM_000088.3").unwrap();
        let b = from_rkyv.get_transcript("NM_000088.3").unwrap();
        assert_eq!(a.contig, b.contig);
        assert_eq!(a.strand, b.strand);
        assert_eq!(a.exons, b.exons);
        assert_eq!((a.cds_start, a.cds_end), (b.cds_start, b.cds_end));
        assert_eq!(a.exon_cigars, b.exon_cigars);
        assert_eq!(a.gene_name, b.gene_name);

        // Alt-build (GRCh37) resolution must survive too.
        let want = from_json
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .map(|t| (t.contig.clone(), t.strand, t.cds_start));
        let got = from_rkyv
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .map(|t| (t.contig.clone(), t.strand, t.cds_start));
        assert_eq!(want, got, "alt-build must survive the rkyv roundtrip");
    }

    #[test]
    fn test_rkyv_rejects_invalid_strand_tag() {
        // A structurally valid, current-version archive whose strand byte is out
        // of range must be rejected (not coerced to Strand::Unknown), so a
        // corrupt cache is regenerated rather than materialized with wrong data.
        let temp = tempfile::TempDir::new().unwrap();
        let p = temp.path().join("cdot.rkyv");
        let bytes = rkyv_cache::archive_bytes_with_tags(7, 0);
        std::fs::write(&p, &bytes).unwrap();

        // The archive itself is current — rejection comes from the tag check.
        assert!(CdotMapper::rkyv_is_current(&p));
        let err = CdotMapper::from_rkyv_file(&p).unwrap_err();
        assert!(
            err.to_string().contains("invalid strand tag"),
            "expected invalid-strand rejection, got: {err}"
        );
    }

    #[test]
    fn test_rkyv_rejects_invalid_cigar_op_tag() {
        // Same as above, but for an out-of-range per-exon cigar op tag.
        let temp = tempfile::TempDir::new().unwrap();
        let p = temp.path().join("cdot.rkyv");
        let bytes = rkyv_cache::archive_bytes_with_tags(1, 9);
        std::fs::write(&p, &bytes).unwrap();

        assert!(CdotMapper::rkyv_is_current(&p));
        let err = CdotMapper::from_rkyv_file(&p).unwrap_err();
        assert!(
            err.to_string().contains("invalid cigar op tag"),
            "expected invalid-cigar-op rejection, got: {err}"
        );
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
    fn load_secondary_build_counts_genome_builds_nested_transcripts() {
        // Regression for the misleading "Loaded 0 GRCh37 transcripts" log. Real
        // cdot files nest their alignments under `genome_builds`, and the
        // path-based `load_secondary_build` loads via `load()` using the default
        // build — so a GRCh37-nested file's *primary* (default-build) view is
        // empty and the transcripts arrive through the alt-build fold. The
        // reported count must reflect what is actually queryable for the build,
        // not the empty primary view (which previously logged "0").
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000001.10",
                            "strand": "+",
                            "exons": [[100, 200, 1, 0, 100, "M100"]]
                        }
                    }
                },
                "NM_000002.1": {
                    "gene_name": "GENE_B",
                    "genome_builds": {
                        "GRCh37": {
                            "contig": "NC_000001.10",
                            "strand": "+",
                            "exons": [[300, 400, 1, 0, 100, "M100"]]
                        }
                    }
                }
            }
        }
        "#;
        // Primary mapper on GRCh38; write the GRCh37 secondary to a temp file
        // (no sibling cache) so `load_secondary_build` parses the JSON via the
        // default-build `load()` path that exhibited the miscount.
        let mut mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        let temp_dir = tempfile::TempDir::new().unwrap();
        let path = temp_dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        let count = mapper.load_secondary_build(&path, "GRCh37").unwrap();
        assert_eq!(
            count, 2,
            "secondary load must count the GRCh37-nested transcripts, not the empty primary view"
        );
        // And the data is genuinely queryable on GRCh37 (not just counted).
        assert!(
            mapper
                .get_transcript_on_build("NM_000002.1", "GRCh37")
                .is_some(),
            "GRCh37-nested transcript must be resolvable after secondary load"
        );
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
    fn test_available_builds_for_includes_deferred_build() {
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000001.10", "strand": "+", "exons": [[100, 200, 1, 0, 100, "M100"]] }
                    }
                }
            }
        }
        "#;
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        let mut lazy = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        lazy.defer_secondary_build("GRCh37", path);

        let mut builds = lazy.available_builds_for("NM_000001.1");
        builds.sort();
        assert_eq!(
            builds,
            vec!["GRCh37".to_string()],
            "deferred GRCh37 must be the sole available build"
        );
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
    fn test_get_transcript_primary_version_fallback_picks_highest() {
        // Same determinism guarantee as the alt-build path, but for the primary
        // build's base->versioned fallback. The cdot file deserializes its
        // transcripts into a HashMap, so `add_transcript` runs in a
        // nondeterministic order; `base_to_versioned` must still map a base
        // accession to its highest version rather than whichever the hasher
        // inserted last. Pre-fix, requesting an absent version (.2) resolved to
        // .3 or .4 arbitrarily across runs, shifting CDS coordinates downstream.
        let json = r#"
        {
            "transcripts": {
                "NM_000088.3": {
                    "gene_name": "COL1A1",
                    "genome_builds": {
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
        // .2 is absent, forcing the base->versioned fallback on the primary build.
        let tx = mapper
            .get_transcript("NM_000088.2")
            .expect("version fallback must find a sibling");
        assert_eq!(
            tx.cds_start,
            Some(20),
            "primary-build base->version fallback must deterministically pick \
             the highest version (.4 over .3)"
        );
    }

    #[test]
    fn test_accession_version_parsing() {
        assert_eq!(accession_version("NM_003002.4"), Some(4));
        assert_eq!(accession_version("NM_003002.10"), Some(10));
        assert_eq!(accession_version("NM_003002"), None);
        assert_eq!(accession_version("LRG_1t1"), None);
        // Versionless accessions both yield `None`; the `>` tie-break comparison
        // in `base_to_versioned` insertion evaluates to `false`, so the existing
        // entry is kept — the fallback is stable for unversioned ids.
        assert_eq!(accession_version("LRG_1t2"), accession_version("LRG_1t1"));
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

    /// Single-exon transcript with a CIGAR `M3 D4 M3`: exon genome `[1000, 1010)`
    /// (10 genome bases), tx `[0, 6)` (6 tx bases). Genome offsets 3..=6 are the
    /// deletion gap; offsets 0..=2 and 7..=9 are matched.
    fn cigar_deletion_transcript(strand: Strand) -> CdotTranscript {
        CdotTranscript {
            gene_name: Some("GAP".to_string()),
            contig: "NC_000001.11".to_string(),
            strand,
            exons: vec![[1000, 1010, 0, 6]],
            cds_start: Some(0),
            cds_end: Some(6),
            exon_cigars: vec![Some(vec![
                CigarOp::Match(3),
                CigarOp::Deletion(4),
                CigarOp::Match(3),
            ])],
            gene_id: None,
            protein: None,
        }
    }

    #[test]
    fn cigar_deletion_gap_detects_inside_only_plus() {
        let tx = cigar_deletion_transcript(Strand::Plus);
        // Inside the deletion: genome positions 1003..=1006.
        for g in 1003..=1006 {
            let gap = tx
                .cigar_deletion_gap_at_genome_pos(0, g)
                .unwrap_or_else(|| panic!("genome {g} should be inside the deletion gap"));
            assert_eq!(gap.kind, CigarGapKind::Deletion);
            assert_eq!(gap.length, 4);
        }
        // Matched (outside the gap): 1000..=1002 and 1007..=1009.
        for g in [1000, 1001, 1002, 1007, 1008, 1009] {
            assert!(
                tx.cigar_deletion_gap_at_genome_pos(0, g).is_none(),
                "genome {g} is matched and must not report a gap"
            );
        }
    }

    #[test]
    fn cigar_deletion_gap_detects_inside_only_minus() {
        let tx = cigar_deletion_transcript(Strand::Minus);
        // On minus strand the CIGAR walks the genome from genome_end-1 down.
        // Deletion genome offsets 3..=6 → positions 1009-6..=1009-3 = 1003..=1006.
        for g in 1003..=1006 {
            let gap = tx
                .cigar_deletion_gap_at_genome_pos(0, g)
                .unwrap_or_else(|| panic!("minus genome {g} should be inside the gap"));
            assert_eq!(gap.kind, CigarGapKind::Deletion);
        }
        for g in [1000, 1001, 1002, 1007, 1008, 1009] {
            assert!(tx.cigar_deletion_gap_at_genome_pos(0, g).is_none());
        }
    }

    #[test]
    fn cigar_deletion_gap_none_without_cigar() {
        let tx = sample_transcript(); // no exon_cigars
        for g in [1000, 1050, 1099] {
            assert!(tx.cigar_deletion_gap_at_genome_pos(0, g).is_none());
        }
    }

    #[test]
    fn cigar_insertion_gap_detects_inside_only() {
        // CIGAR `M3 I2 M15`: exon genome `[1000, 1018)` (18 genome bases),
        // tx `[0, 20)` (20 tx bases). Tx offsets 3..=4 are the insertion gap
        // (transcript bases with no genome counterpart).
        let tx = CdotTranscript {
            gene_name: Some("INS".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1018, 0, 20]],
            cds_start: Some(0),
            cds_end: Some(20),
            exon_cigars: vec![Some(vec![
                CigarOp::Match(3),
                CigarOp::Insertion(2),
                CigarOp::Match(15),
            ])],
            gene_id: None,
            protein: None,
        };
        // Inside the insertion: tx positions 3 and 4.
        for t in [3, 4] {
            let gap = tx
                .cigar_insertion_gap_at_tx_pos(0, t)
                .unwrap_or_else(|| panic!("tx {t} should be inside the insertion gap"));
            assert_eq!(gap.kind, CigarGapKind::Insertion);
            assert_eq!(gap.length, 2);
        }
        // Matched positions before (0..=2) and after (5..) the insertion.
        for t in [0, 1, 2, 5, 6, 19] {
            assert!(
                tx.cigar_insertion_gap_at_tx_pos(0, t).is_none(),
                "tx {t} is matched and must not report an insertion gap"
            );
        }
    }

    // -------------------------------------------------------------------------
    // Multi-build overlap discovery: `transcripts_at_position` must route by
    // contig to the build that owns it. GRCh37 and GRCh38 have distinct contig
    // accessions (NC_000017.10 vs .11), so the contig disambiguates the build.
    // Pre-fix, the stab index was built from the primary build's `contig_index`
    // only, so a GRCh37 contig query returned empty.
    // -------------------------------------------------------------------------

    #[test]
    fn test_transcripts_at_position_routes_to_owning_build_by_contig() {
        // Primary build = GRCh38; the same NM lives on GRCh37 in the alt map.
        let mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();

        // Primary (GRCh38) contig: a position inside the GRCh38 exon must hit.
        let hits38 = mapper.transcripts_at_position("NC_000017.11", 50184100);
        let ids38: Vec<&str> = hits38.iter().map(|(acc, _)| *acc).collect();
        assert_eq!(
            ids38,
            vec!["NM_000088.3"],
            "primary-build overlap must still resolve"
        );
        // The returned transcript carries the GRCh38 alignment.
        assert_eq!(hits38[0].1.contig, "NC_000017.11");

        // Alt (GRCh37) contig: a position inside the GRCh37 exon must hit, and
        // the returned transcript must carry the GRCh37 alignment.
        let hits37 = mapper.transcripts_at_position("NC_000017.10", 48263050);
        let ids37: Vec<&str> = hits37.iter().map(|(acc, _)| *acc).collect();
        assert_eq!(
            ids37,
            vec!["NM_000088.3"],
            "alt-build overlap must resolve via the contig (pre-fix returned empty)"
        );
        assert_eq!(hits37[0].1.contig, "NC_000017.10");

        // A contig owned by neither build resolves to empty.
        assert!(
            mapper
                .transcripts_at_position("NC_000099.9", 1000)
                .is_empty(),
            "unknown contig must return empty"
        );

        // A position outside the alt-build exon on a known alt contig is empty.
        assert!(
            mapper
                .transcripts_at_position("NC_000017.10", 1_000)
                .is_empty(),
            "position outside the alt-build exon must return empty"
        );
    }

    // -------------------------------------------------------------------------
    // Loading a second, single-build cdot file as a non-primary build. The
    // real manifest ships GRCh37 and GRCh38 cdot files separately (single-build
    // each, not nested `genome_builds`), so the loader must be able to layer a
    // second file onto an already-loaded primary mapper.
    // -------------------------------------------------------------------------

    fn single_build_grch38_json() -> &'static str {
        r#"
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
        "#
    }

    fn single_build_grch37_json() -> &'static str {
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
    fn test_load_secondary_build_layers_second_cdot_file() {
        // Primary build = GRCh38, loaded from a single-build file.
        let mut mapper = CdotMapper::from_reader(single_build_grch38_json().as_bytes()).unwrap();
        // Sanity: only GRCh38 known before layering.
        assert_eq!(
            mapper.available_builds_for("NM_000088.3"),
            vec!["GRCh38".to_string()]
        );

        // Layer GRCh37 on top.
        let count = mapper
            .load_secondary_build_from_reader(single_build_grch37_json().as_bytes(), "GRCh37")
            .unwrap();
        assert_eq!(count, 1, "one GRCh37 transcript loaded");

        // Both builds now available, and each returns its own contig.
        let mut builds = mapper.available_builds_for("NM_000088.3");
        builds.sort();
        assert_eq!(builds, vec!["GRCh37".to_string(), "GRCh38".to_string()]);

        assert_eq!(
            mapper
                .get_transcript_on_build("NM_000088.3", "GRCh38")
                .unwrap()
                .contig,
            "NC_000017.11"
        );
        assert_eq!(
            mapper
                .get_transcript_on_build("NM_000088.3", "GRCh37")
                .unwrap()
                .contig,
            "NC_000017.10"
        );

        // The primary `transcripts` map must NOT be clobbered: the GRCh38
        // record's contig is unchanged (collision safety).
        assert_eq!(
            mapper.get_transcript("NM_000088.3").unwrap().contig,
            "NC_000017.11"
        );

        // Overlap discovery routes by contig to the freshly-loaded GRCh37 build.
        let hits = mapper.transcripts_at_position("NC_000017.10", 48263050);
        let ids: Vec<&str> = hits.iter().map(|(acc, _)| *acc).collect();
        assert_eq!(ids, vec!["NM_000088.3"]);
        assert_eq!(hits[0].1.contig, "NC_000017.10");
    }

    // -------------------------------------------------------------------------
    // Deferred secondary-build: lazy-load tests
    // -------------------------------------------------------------------------

    #[test]
    fn deferred_alt_mapper_loads_once_on_first_access() {
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000001.10", "strand": "+", "exons": [[100, 200, 1, 0, 100, "M100"]] }
                    }
                }
            }
        }
        "#;
        let mut mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        mapper.defer_secondary_build("GRCh37", path);
        assert!(
            !mapper.deferred_alt_loaded("GRCh37"),
            "deferral must not load eagerly"
        );

        let tx = mapper
            .deferred_alt_mapper("GRCh37")
            .and_then(|m| m.get_transcript_on_build("NM_000001.1", "GRCh37"));
        assert!(
            tx.is_some(),
            "deferred mapper must expose its GRCh37 transcript"
        );
        assert!(
            mapper.deferred_alt_loaded("GRCh37"),
            "first access must load it"
        );

        assert!(mapper.deferred_alt_mapper("GRCh37").is_some());
        assert!(mapper.deferred_alt_loaded("GRCh37"));
    }

    #[test]
    fn deferred_alt_mapper_missing_source_resolves_to_none() {
        let mut mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        mapper.defer_secondary_build("GRCh37", PathBuf::from("/no/such/cdot.json"));
        assert!(
            mapper.deferred_alt_mapper("GRCh37").is_none(),
            "missing source -> None, no panic"
        );
        assert!(mapper.deferred_alt_loaded("GRCh37"));
        assert!(mapper.deferred_alt_mapper("GRCh99").is_none());
    }

    #[test]
    fn get_transcript_on_build_deferred_matches_eager() {
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000001.10", "strand": "+", "exons": [[100, 200, 1, 0, 100, "M100"]] }
                    }
                }
            }
        }
        "#;
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        let mut eager = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        eager.load_secondary_build(&path, "GRCh37").unwrap();

        let mut lazy = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        lazy.defer_secondary_build("GRCh37", path.clone());

        let e = eager
            .get_transcript_on_build("NM_000001.1", "GRCh37")
            .unwrap();
        let l = lazy
            .get_transcript_on_build("NM_000001.1", "GRCh37")
            .unwrap();
        assert_eq!(e.contig, l.contig);
        assert_eq!(e.contig, "NC_000001.10");
        assert_eq!(e.strand, l.strand);
        assert_eq!(e.exons, l.exons);

        let mut lazy2 = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        lazy2.defer_secondary_build("GRCh37", path);
        // A GRCh37 lookup satisfied by the PRIMARY cdot's own eager GRCh37 alt view
        // must NOT spill into (load) the deferred secondary.
        let _ = lazy2
            .get_transcript_on_build("NM_000088.3", "GRCh37")
            .unwrap();
        assert!(
            !lazy2.deferred_alt_loaded("GRCh37"),
            "an eagerly-satisfiable GRCh37 lookup must not load the deferred secondary"
        );
    }

    #[test]
    fn test_overlap_index_includes_deferred_grch37_contigs() {
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000001.10", "strand": "+", "exons": [[100, 200, 1, 0, 100, "M100"]] }
                    }
                }
            }
        }
        "#;
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        let mut eager = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        eager.load_secondary_build(&path, "GRCh37").unwrap();
        let eager_hits = eager.transcripts_at_position("NC_000001.10", 150);

        let mut lazy = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        lazy.defer_secondary_build("GRCh37", path);
        let lazy_hits = lazy.transcripts_at_position("NC_000001.10", 150);

        let mut e: Vec<&str> = eager_hits.iter().map(|(a, _)| *a).collect();
        let mut l: Vec<&str> = lazy_hits.iter().map(|(a, _)| *a).collect();
        e.sort();
        l.sort();
        assert_eq!(e, l, "deferred overlap must match eager overlap");
        assert!(
            lazy.deferred_alt_loaded("GRCh37"),
            "an overlap query must force-load the deferred secondary"
        );
        assert!(
            l.contains(&"NM_000001.1"),
            "overlap result must include NM_000001.1 at position 150"
        );
    }

    #[test]
    fn test_overlap_merges_eager_and_deferred_same_contig() {
        // The primary (multi_build) cdot puts NM_000088.3 on GRCh37 contig
        // NC_000017.10 (eager alt view). A deferred GRCh37 file adds another
        // transcript on the SAME contig. Both must appear in an overlap query —
        // the deferred index must MERGE with, not replace, the eager contig stab.
        //
        // Regression: the old `index_one_build`-based approach called
        // `HashMap::insert` per contig, so the deferred pass silently overwrote
        // the eager stab for NC_000017.10, discarding NM_000088.3.
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_999999.1": {
                    "gene_name": "GENE_X",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000017.10", "strand": "+", "exons": [[48263025, 48263098, 1, 0, 73, "M73"]] }
                    }
                }
            }
        }
        "#;
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        let mut lazy = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        lazy.defer_secondary_build("GRCh37", path);

        // A position inside NM_000088.3's GRCh37 exon (NC_000017.10:48263025-48263098).
        // NM_999999.1 has an identical exon, so both transcripts overlap pos 48263050.
        let hits = lazy.transcripts_at_position("NC_000017.10", 48263050);
        let mut accs: Vec<&str> = hits.iter().map(|(a, _)| *a).collect();
        accs.sort();
        assert!(
            accs.contains(&"NM_000088.3"),
            "eager GRCh37 transcript must not be dropped by the deferred merge; got: {accs:?}"
        );
        assert!(
            accs.contains(&"NM_999999.1"),
            "deferred GRCh37 transcript must be present; got: {accs:?}"
        );
    }

    #[test]
    fn deferred_grch37_not_loaded_for_grch38_only_sequence() {
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000001.10", "strand": "+", "exons": [[100, 200, 1, 0, 100, "M100"]] }
                    }
                }
            }
        }
        "#;
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();

        // A sequence of GRCh38 / primary-build-only lookups must NOT load the deferred GRCh37 mapper.
        let mut m = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        m.defer_secondary_build("GRCh37", path);
        let _ = m.get_transcript("NM_000088.3");
        let _ = m.get_transcript_on_build("NM_000088.3", "GRCh38");
        assert!(
            !m.deferred_alt_loaded("GRCh37"),
            "GRCh38/primary lookups must not load the GRCh37 secondary"
        );
    }

    #[test]
    fn deferred_grch37_loads_once_under_concurrency() {
        use std::sync::Arc;
        let secondary_grch37_json = r#"
        {
            "transcripts": {
                "NM_000001.1": {
                    "gene_name": "GENE_A",
                    "genome_builds": {
                        "GRCh37": { "contig": "NC_000001.10", "strand": "+", "exons": [[100, 200, 1, 0, 100, "M100"]] }
                    }
                }
            }
        }
        "#;
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("grch37.json");
        std::fs::write(&path, secondary_grch37_json).unwrap();
        let mut mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        mapper.defer_secondary_build("GRCh37", path);
        let shared = Arc::new(mapper);

        let mut handles = Vec::new();
        for _ in 0..8 {
            let m = Arc::clone(&shared);
            handles.push(std::thread::spawn(move || {
                m.get_transcript_on_build("NM_000001.1", "GRCh37")
                    .map(|t| t.contig.clone())
            }));
        }
        for h in handles {
            assert_eq!(h.join().unwrap().as_deref(), Some("NC_000001.10"));
        }
        assert!(shared.deferred_alt_loaded("GRCh37"));
    }

    #[test]
    fn deferred_state_is_not_serialized() {
        // The deferred source is runtime-only: serializing and reloading drops it
        // (from_manifest re-establishes it). Round-trip must not panic and must yield
        // a mapper with no deferred state, while primary data survives.
        let mut mapper = CdotMapper::from_reader(multi_build_cdot_json().as_bytes()).unwrap();
        mapper.defer_secondary_build("GRCh37", PathBuf::from("/tmp/whatever.json"));

        let dir = tempfile::TempDir::new().unwrap();
        let rkyv_path = dir.path().join("m.rkyv");
        mapper.to_rkyv_file(&rkyv_path).unwrap();
        let reloaded = CdotMapper::from_rkyv_file(&rkyv_path).unwrap();

        // Deferred state did not survive serialization.
        assert!(
            !reloaded.deferred_alt_loaded("GRCh37"),
            "deferred state must not survive serialization"
        );
        assert!(
            reloaded.deferred_alt_mapper("GRCh37").is_none(),
            "no deferred source should be reconstructed from a snapshot"
        );
        // Primary data survives normally.
        assert!(reloaded.get_transcript("NM_000088.3").is_some());
    }
}
