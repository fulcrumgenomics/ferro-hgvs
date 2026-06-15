//! Transcript database loading
//!
//! This module provides functionality to load transcript annotations from
//! standard file formats (GFF3, GTF) commonly used for genome annotations.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use bincode::Options;
use serde::{Deserialize, Serialize};
use superintervals::IntervalMap;

use crate::error::FerroError;
use crate::reference::annotation::AnnotationFormat;
use crate::reference::transcript::{GenomeBuild, ManeStatus, Strand, Transcript};

/// Magic + schema version prefixing a TranscriptDb binary cache. Bump the
/// version whenever the serialized layout of [`TranscriptDb`] (or `Transcript`)
/// changes, so a stale cache is detected and regenerated rather than mis-read —
/// the same guard the cdot cache uses.
const TDB_CACHE_MAGIC: [u8; 4] = *b"FTDB";
const TDB_CACHE_VERSION: u32 = 1;

/// Statistics about MANE transcript coverage in the database
#[derive(Debug, Clone, Default)]
pub struct ManeStats {
    /// Total number of unique genes
    pub total_genes: usize,
    /// Number of genes with MANE Select transcript
    pub genes_with_mane_select: usize,
    /// Number of genes with MANE Plus Clinical transcripts
    pub genes_with_mane_plus_clinical: usize,
    /// Total number of MANE Plus Clinical transcripts
    pub total_mane_plus_clinical_transcripts: usize,
}

impl ManeStats {
    /// Get MANE Select coverage percentage
    pub fn mane_select_coverage(&self) -> f64 {
        if self.total_genes == 0 {
            0.0
        } else {
            (self.genes_with_mane_select as f64 / self.total_genes as f64) * 100.0
        }
    }
}

/// Bincode-friendly snapshot of a `Transcript` — the domain type carries many
/// `#[serde(skip_serializing_if)]` attributes, which bincode (a positional
/// format) cannot round-trip, so the cache serializes this no-skip mirror
/// instead. `exon_cigars` and `cached_introns` are intentionally omitted: the
/// former is empty for GFF/GTF-loaded transcripts and the latter is rebuilt
/// lazily.
#[derive(Serialize, Deserialize)]
struct TranscriptSnapshot {
    id: String,
    gene_symbol: Option<String>,
    strand: Strand,
    sequence: Option<String>,
    cds_start: Option<u64>,
    cds_end: Option<u64>,
    exons: Vec<ExonSnapshot>,
    chromosome: Option<String>,
    genomic_start: Option<u64>,
    genomic_end: Option<u64>,
    genome_build: GenomeBuild,
    mane_status: ManeStatus,
    refseq_match: Option<String>,
    ensembl_match: Option<String>,
    protein_id: Option<String>,
}

#[derive(Serialize, Deserialize)]
struct ExonSnapshot {
    number: u32,
    start: u64,
    end: u64,
    genomic_start: Option<u64>,
    genomic_end: Option<u64>,
}

/// Bincode-friendly snapshot of [`TranscriptDb`]. The index maps carry no skip
/// attributes, so only `transcripts` needs the mirror above.
#[derive(Serialize, Deserialize)]
struct TranscriptDbSnapshot {
    transcripts: HashMap<String, TranscriptSnapshot>,
    gene_index: HashMap<String, Vec<String>>,
    region_index: HashMap<String, Vec<(u64, u64, String)>>,
    mane_select_index: HashMap<String, String>,
    mane_plus_clinical_index: HashMap<String, Vec<String>>,
    genome_build: GenomeBuild,
}

impl From<&Transcript> for TranscriptSnapshot {
    fn from(t: &Transcript) -> Self {
        Self {
            id: t.id.clone(),
            gene_symbol: t.gene_symbol.clone(),
            strand: t.strand,
            sequence: t.sequence.clone(),
            cds_start: t.cds_start,
            cds_end: t.cds_end,
            exons: t
                .exons
                .iter()
                .map(|e| ExonSnapshot {
                    number: e.number,
                    start: e.start,
                    end: e.end,
                    genomic_start: e.genomic_start,
                    genomic_end: e.genomic_end,
                })
                .collect(),
            chromosome: t.chromosome.clone(),
            genomic_start: t.genomic_start,
            genomic_end: t.genomic_end,
            genome_build: t.genome_build,
            mane_status: t.mane_status,
            refseq_match: t.refseq_match.clone(),
            ensembl_match: t.ensembl_match.clone(),
            protein_id: t.protein_id.clone(),
        }
    }
}

impl From<TranscriptSnapshot> for Transcript {
    fn from(s: TranscriptSnapshot) -> Self {
        Transcript {
            id: s.id,
            gene_symbol: s.gene_symbol,
            strand: s.strand,
            sequence: s.sequence,
            cds_start: s.cds_start,
            cds_end: s.cds_end,
            exons: s
                .exons
                .into_iter()
                .map(|e| crate::reference::transcript::Exon {
                    number: e.number,
                    start: e.start,
                    end: e.end,
                    genomic_start: e.genomic_start,
                    genomic_end: e.genomic_end,
                })
                .collect(),
            chromosome: s.chromosome,
            genomic_start: s.genomic_start,
            genomic_end: s.genomic_end,
            genome_build: s.genome_build,
            mane_status: s.mane_status,
            refseq_match: s.refseq_match,
            ensembl_match: s.ensembl_match,
            protein_id: s.protein_id,
            ..Default::default()
        }
    }
}

/// A database of transcripts indexed for efficient lookup
#[derive(Debug, Default)]
pub struct TranscriptDb {
    /// Transcripts indexed by accession ID
    transcripts: HashMap<String, Transcript>,
    /// Index from gene symbol to transcript IDs
    gene_index: HashMap<String, Vec<String>>,
    /// Index from chromosome:region to transcript IDs for overlap queries
    region_index: HashMap<String, Vec<(u64, u64, String)>>,
    /// Index from gene symbol to MANE Select transcript ID
    mane_select_index: HashMap<String, String>,
    /// Index from gene symbol to MANE Plus Clinical transcript IDs
    mane_plus_clinical_index: HashMap<String, Vec<String>>,
    /// Genome build
    pub genome_build: GenomeBuild,
    /// Lazily-built per-contig interval index for `get_by_region`, derived from
    /// `region_index`. The flat `region_index` is a linear scan per query; this
    /// turns each overlap query into O(log n + k). Runtime-only — not part of
    /// the serialized snapshot/cache (rebuilt on first query after a load).
    /// `take()`n whenever `region_index` changes (see `add_transcript`).
    region_query_index: OnceLock<HashMap<String, RegionIntervalIndex>>,
}

/// One contig's interval index plus the transcript ids its `u32` payloads point
/// into. Built from `TranscriptDb::region_index` by `build_region_query_index`.
#[derive(Debug)]
struct RegionIntervalIndex {
    /// Interval index over transcript genomic spans; payload indexes `ids`.
    index: IntervalMap<u32>,
    /// Transcript ids in the order the `index` payloads address them.
    ids: Vec<String>,
}

/// Build the per-contig interval index from the flat `region_index`. Each
/// transcript span is added with a `u32` payload indexing the contig's `ids`
/// vec; intervals whose coordinates don't fit `i32` (well beyond any real
/// genomic position) are skipped rather than wrapped.
fn build_region_query_index(
    region_index: &HashMap<String, Vec<(u64, u64, String)>>,
) -> HashMap<String, RegionIntervalIndex> {
    let mut out = HashMap::with_capacity(region_index.len());
    for (contig, regions) in region_index {
        let mut index: IntervalMap<u32> = IntervalMap::new();
        let mut ids: Vec<String> = Vec::with_capacity(regions.len());
        for (start, end, id) in regions {
            let (Ok(s), Ok(e)) = (i32::try_from(*start), i32::try_from(*end)) else {
                continue;
            };
            let payload = ids.len() as u32;
            ids.push(id.clone());
            index.add(s, e, payload);
        }
        index.build();
        out.insert(contig.clone(), RegionIntervalIndex { index, ids });
    }
    out
}

impl TranscriptDb {
    /// Create a new empty transcript database
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new database with a specified genome build
    pub fn with_build(build: GenomeBuild) -> Self {
        Self {
            genome_build: build,
            ..Default::default()
        }
    }

    /// Sibling binary-cache path for a source annotation file at `build` parsed
    /// as `format` (`<source>.<build>.<format>.tdb`). Building the cache key off
    /// the full path (not `with_extension`) avoids over-stripping multi-dot
    /// filenames. The `format` component keeps a cache built for one parser
    /// format from being served for a different requested format when the same
    /// source path is loaded with an explicit `with_format` override.
    pub fn cache_path_for<P: AsRef<Path>>(
        source: P,
        build: GenomeBuild,
        format: AnnotationFormat,
    ) -> PathBuf {
        let format_tag = match format {
            AnnotationFormat::Gff3 => "gff3",
            AnnotationFormat::Gtf => "gtf",
        };
        PathBuf::from(format!(
            "{}.{}.{}.tdb",
            source.as_ref().display(),
            build,
            format_tag
        ))
    }

    /// Whether `cache_path` is a usable, up-to-date cache for `source`: it must
    /// exist, carry the current magic+version, and be at least as new as the
    /// source file (so editing the GFF/GTF invalidates the cache).
    pub fn cache_is_fresh<P: AsRef<Path>, Q: AsRef<Path>>(cache_path: P, source: Q) -> bool {
        let cache_path = cache_path.as_ref();
        let Ok(cache_meta) = std::fs::metadata(cache_path) else {
            return false;
        };
        // Treat an unreadable mtime (cache or source) as "not fresh": if the
        // source file is missing or unreadable we must not serve a stale cache,
        // and `load_annotations` returns early on a fresh cache, so falling
        // through here would mask the source-file error.
        let Ok(cm) = cache_meta.modified() else {
            return false;
        };
        let Ok(sm) = std::fs::metadata(source.as_ref()).and_then(|m| m.modified()) else {
            return false;
        };
        if cm < sm {
            return false; // source edited after the cache was written
        }
        // Verify the header.
        let Ok(mut file) = File::open(cache_path) else {
            return false;
        };
        let mut header = [0u8; 8];
        if file.read_exact(&mut header).is_err() {
            return false;
        }
        header[..4] == TDB_CACHE_MAGIC
            && u32::from_le_bytes([header[4], header[5], header[6], header[7]]) == TDB_CACHE_VERSION
    }

    /// Load a `TranscriptDb` from a binary cache, verifying the magic+version
    /// header first so a stale/foreign file is rejected cleanly.
    pub fn from_cache_file<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
            msg: format!("Failed to open TranscriptDb cache: {}", e),
        })?;
        let file_size = file.metadata().map(|m| m.len()).unwrap_or(0);
        let mut reader = BufReader::new(file);

        let mut header = [0u8; 8];
        reader.read_exact(&mut header).map_err(|e| FerroError::Io {
            msg: format!("TranscriptDb cache too short for a header: {}", e),
        })?;
        if header[..4] != TDB_CACHE_MAGIC {
            return Err(FerroError::Io {
                msg: "TranscriptDb cache has no valid magic".to_string(),
            });
        }
        let version = u32::from_le_bytes([header[4], header[5], header[6], header[7]]);
        if version != TDB_CACHE_VERSION {
            return Err(FerroError::Io {
                msg: format!(
                    "TranscriptDb cache version {} != expected {}",
                    version, TDB_CACHE_VERSION
                ),
            });
        }

        // 40x: TranscriptDb carries five HashMap indexes (transcripts, by_name,
        // by_chrom, canonical_ids, gene_to_transcripts) vs cdot's one, so the
        // in-memory size relative to the serialized file is proportionally larger.
        let size_limit = file_size.saturating_mul(40).max(1024 * 1024);
        let snap: TranscriptDbSnapshot = bincode::options()
            .with_fixint_encoding()
            .allow_trailing_bytes()
            .with_limit(size_limit)
            .deserialize_from(reader)
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to deserialize TranscriptDb cache: {}", e),
            })?;
        Ok(TranscriptDb {
            transcripts: snap
                .transcripts
                .into_iter()
                .map(|(k, v)| (k, v.into()))
                .collect(),
            gene_index: snap.gene_index,
            region_index: snap.region_index,
            mane_select_index: snap.mane_select_index,
            mane_plus_clinical_index: snap.mane_plus_clinical_index,
            genome_build: snap.genome_build,
            region_query_index: OnceLock::new(),
        })
    }

    /// Write this database to a binary cache (atomically: temp + rename).
    pub fn to_cache_file<P: AsRef<Path>>(&self, path: P) -> Result<(), FerroError> {
        // Unique temp sibling (pid + process-global counter) so concurrent
        // writers in the same process never collide on the temp name.
        use std::sync::atomic::{AtomicU64, Ordering};
        static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);
        let path = path.as_ref();
        let mut tmp = path.to_path_buf();
        let tmp_name = format!(
            "{}.tmp.{}.{}",
            path.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "transcriptdb.tdb".to_string()),
            std::process::id(),
            TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        );
        tmp.set_file_name(tmp_name);

        let file = File::create(&tmp).map_err(|e| FerroError::Io {
            msg: format!("Failed to create TranscriptDb cache temp file: {}", e),
        })?;
        let snap = TranscriptDbSnapshot {
            transcripts: self
                .transcripts
                .iter()
                .map(|(k, v)| (k.clone(), TranscriptSnapshot::from(v)))
                .collect(),
            gene_index: self.gene_index.clone(),
            region_index: self.region_index.clone(),
            mane_select_index: self.mane_select_index.clone(),
            mane_plus_clinical_index: self.mane_plus_clinical_index.clone(),
            genome_build: self.genome_build,
        };
        let mut writer = std::io::BufWriter::new(file);
        let result = (|| -> std::io::Result<()> {
            writer.write_all(&TDB_CACHE_MAGIC)?;
            writer.write_all(&TDB_CACHE_VERSION.to_le_bytes())?;
            bincode::options()
                .with_fixint_encoding()
                .serialize_into(&mut writer, &snap)
                .map_err(std::io::Error::other)?;
            writer.flush()
        })();
        if let Err(e) = result {
            let _ = std::fs::remove_file(&tmp);
            return Err(FerroError::Io {
                msg: format!("Failed to serialize TranscriptDb cache: {}", e),
            });
        }
        std::fs::rename(&tmp, path).map_err(|e| {
            let _ = std::fs::remove_file(&tmp);
            FerroError::Io {
                msg: format!("Failed to finalize TranscriptDb cache: {}", e),
            }
        })
    }

    /// Add a transcript to the database
    pub fn add(&mut self, transcript: Transcript) {
        // Invalidate the lazily-built region interval index: `region_index` is
        // about to change, so any cached index is stale and must be rebuilt on
        // the next query.
        self.region_query_index.take();

        let id = transcript.id.clone();

        // Index by gene symbol
        if let Some(ref gene) = transcript.gene_symbol {
            self.gene_index
                .entry(gene.clone())
                .or_default()
                .push(id.clone());

            // Index MANE transcripts
            match transcript.mane_status {
                ManeStatus::Select => {
                    self.mane_select_index.insert(gene.clone(), id.clone());
                }
                ManeStatus::PlusClinical => {
                    self.mane_plus_clinical_index
                        .entry(gene.clone())
                        .or_default()
                        .push(id.clone());
                }
                ManeStatus::None => {}
            }
        }

        // Index by genomic region
        if let (Some(ref chrom), Some(start), Some(end)) = (
            &transcript.chromosome,
            transcript.genomic_start,
            transcript.genomic_end,
        ) {
            self.region_index
                .entry(chrom.clone())
                .or_default()
                .push((start, end, id.clone()));
        }

        self.transcripts.insert(id, transcript);
    }

    /// Get a transcript by its accession ID
    pub fn get(&self, id: &str) -> Option<&Transcript> {
        self.transcripts.get(id)
    }

    /// Get all transcripts for a gene
    pub fn get_by_gene(&self, gene: &str) -> Vec<&Transcript> {
        self.gene_index
            .get(gene)
            .map(|ids| {
                ids.iter()
                    .filter_map(|id| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all transcripts overlapping a genomic position
    pub fn get_by_position(&self, chrom: &str, pos: u64) -> Vec<&Transcript> {
        self.region_index
            .get(chrom)
            .map(|regions| {
                regions
                    .iter()
                    .filter(|(start, end, _)| pos >= *start && pos <= *end)
                    .filter_map(|(_, _, id)| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all transcripts overlapping a genomic region
    pub fn get_by_region(&self, chrom: &str, start: u64, end: u64) -> Vec<&Transcript> {
        // Resolve the per-contig interval index (built once, lazily, from
        // `region_index`), then range-query it. This replaces a linear scan over
        // every transcript on the contig with an O(log n + k) overlap query —
        // worth it because annotation fires this per variant.
        let Some(contig) = self
            .region_query_index
            .get_or_init(|| build_region_query_index(&self.region_index))
            .get(chrom)
        else {
            return Vec::new();
        };
        // superintervals uses i32 coordinates; clamp defensively (human genomic
        // positions fit comfortably, but a malformed huge coordinate must not
        // wrap). An out-of-range query simply returns no hits.
        let (Ok(qs), Ok(qe)) = (i32::try_from(start), i32::try_from(end)) else {
            return Vec::new();
        };
        let mut hits: Vec<u32> = Vec::new();
        contig.index.search_values(qs, qe, &mut hits);
        hits.iter()
            .filter_map(|&i| contig.ids.get(i as usize))
            .filter_map(|id| self.transcripts.get(id))
            .collect()
    }

    /// Get the MANE Select transcript for a gene
    ///
    /// Returns the single representative transcript designated by NCBI/EBI
    /// for this gene, if available.
    pub fn get_mane_select(&self, gene: &str) -> Option<&Transcript> {
        self.mane_select_index
            .get(gene)
            .and_then(|id| self.transcripts.get(id))
    }

    /// Get all MANE Plus Clinical transcripts for a gene
    ///
    /// Returns additional clinically relevant transcripts beyond MANE Select.
    pub fn get_mane_plus_clinical(&self, gene: &str) -> Vec<&Transcript> {
        self.mane_plus_clinical_index
            .get(gene)
            .map(|ids| {
                ids.iter()
                    .filter_map(|id| self.transcripts.get(id))
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Get all MANE transcripts for a gene (Select + Plus Clinical)
    pub fn get_mane_transcripts(&self, gene: &str) -> Vec<&Transcript> {
        let mut result = Vec::new();

        if let Some(select) = self.get_mane_select(gene) {
            result.push(select);
        }

        result.extend(self.get_mane_plus_clinical(gene));
        result
    }

    /// Get the best transcript for a gene, preferring MANE Select
    ///
    /// Priority order:
    /// 1. MANE Select
    /// 2. MANE Plus Clinical (first one)
    /// 3. Any coding transcript (first one)
    /// 4. Any transcript (first one)
    pub fn get_preferred_transcript(&self, gene: &str) -> Option<&Transcript> {
        // Try MANE Select first
        if let Some(tx) = self.get_mane_select(gene) {
            return Some(tx);
        }

        // Try MANE Plus Clinical
        let plus_clinical = self.get_mane_plus_clinical(gene);
        if !plus_clinical.is_empty() {
            return Some(plus_clinical[0]);
        }

        // Fall back to gene transcripts
        let transcripts = self.get_by_gene(gene);

        // Prefer coding transcripts
        if let Some(coding) = transcripts.iter().find(|tx| tx.is_coding()) {
            return Some(coding);
        }

        // Any transcript
        transcripts.first().copied()
    }

    /// Check if a gene has a MANE Select transcript
    pub fn has_mane_select(&self, gene: &str) -> bool {
        self.mane_select_index.contains_key(gene)
    }

    /// Get statistics about MANE coverage
    pub fn mane_stats(&self) -> ManeStats {
        let total_genes = self.gene_index.len();
        let genes_with_mane_select = self.mane_select_index.len();
        let genes_with_mane_plus = self.mane_plus_clinical_index.len();
        let total_mane_plus_transcripts: usize = self
            .mane_plus_clinical_index
            .values()
            .map(|v| v.len())
            .sum();

        ManeStats {
            total_genes,
            genes_with_mane_select,
            genes_with_mane_plus_clinical: genes_with_mane_plus,
            total_mane_plus_clinical_transcripts: total_mane_plus_transcripts,
        }
    }

    /// Get the number of transcripts in the database
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Check if the database is empty
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Get all transcript IDs
    pub fn ids(&self) -> impl Iterator<Item = &String> {
        self.transcripts.keys()
    }

    /// Iterate over all transcripts
    pub fn iter(&self) -> impl Iterator<Item = (&String, &Transcript)> {
        self.transcripts.iter()
    }
}

/// Detect genome build from file header or content
pub fn detect_genome_build<P: AsRef<Path>>(path: P) -> Result<GenomeBuild, FerroError> {
    let file = File::open(path.as_ref()).map_err(|e| FerroError::Io {
        msg: format!("Failed to open file: {}", e),
    })?;
    let reader = BufReader::new(file);

    for line in reader.lines().take(100) {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        // Check for genome build indicators
        let line_lower = line.to_lowercase();
        if line_lower.contains("grch37") || line_lower.contains("hg19") {
            return Ok(GenomeBuild::GRCh37);
        }
        if line_lower.contains("grch38") || line_lower.contains("hg38") {
            return Ok(GenomeBuild::GRCh38);
        }
    }

    Ok(GenomeBuild::Unknown)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::transcript::{Exon, Strand};
    use std::sync::OnceLock;

    #[test]
    fn test_transcript_db_new() {
        let db = TranscriptDb::new();
        assert!(db.is_empty());
    }

    #[test]
    fn test_transcript_db_add_and_get() {
        let mut db = TranscriptDb::new();

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        assert_eq!(db.len(), 1);
        assert!(db.get("NM_000088.3").is_some());
    }

    /// Build a minimal transcript spanning `[start, end]` on `chrom`.
    fn region_tx(id: &str, chrom: &str, start: u64, end: u64) -> Transcript {
        Transcript {
            id: id.to_string(),
            gene_symbol: Some(id.to_string()),
            strand: Strand::Plus,
            sequence: None,
            cds_start: None,
            cds_end: None,
            exons: Vec::new(),
            chromosome: Some(chrom.to_string()),
            genomic_start: Some(start),
            genomic_end: Some(end),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        }
    }

    /// Sorted ids returned by `get_by_region`, for order-independent assertions.
    fn region_ids(db: &TranscriptDb, chrom: &str, start: u64, end: u64) -> Vec<String> {
        let mut ids: Vec<String> = db
            .get_by_region(chrom, start, end)
            .iter()
            .map(|t| t.id.clone())
            .collect();
        ids.sort();
        ids
    }

    #[test]
    fn get_by_region_interval_index_overlap_semantics() {
        let mut db = TranscriptDb::new();
        db.add(region_tx("A", "chr1", 1000, 2000));
        db.add(region_tx("B", "chr1", 1500, 2500));
        db.add(region_tx("C", "chr1", 5000, 6000));
        db.add(region_tx("D", "chr2", 1000, 2000));

        // Overlaps both A and B.
        assert_eq!(region_ids(&db, "chr1", 1800, 1900), vec!["A", "B"]);
        // Touches only A.
        assert_eq!(region_ids(&db, "chr1", 1000, 1100), vec!["A"]);
        // Boundary: query end == A.start and query start == B-region gap — A only
        // (inclusive overlap, matching the previous linear scan).
        assert_eq!(region_ids(&db, "chr1", 500, 1000), vec!["A"]);
        // Between B and C: no overlap.
        assert!(region_ids(&db, "chr1", 3000, 4000).is_empty());
        // Spanning everything on chr1.
        assert_eq!(region_ids(&db, "chr1", 0, 100_000), vec!["A", "B", "C"]);
        // Wrong contig / unknown contig.
        assert_eq!(region_ids(&db, "chr2", 1500, 1600), vec!["D"]);
        assert!(region_ids(&db, "chrX", 1500, 1600).is_empty());
    }

    #[test]
    fn get_by_region_index_rebuilds_after_add() {
        let mut db = TranscriptDb::new();
        db.add(region_tx("A", "chr1", 1000, 2000));
        // Build the lazy index via a query.
        assert_eq!(region_ids(&db, "chr1", 1500, 1600), vec!["A"]);
        // Adding must invalidate the cached index so the new transcript is seen.
        db.add(region_tx("B", "chr1", 1400, 1700));
        assert_eq!(region_ids(&db, "chr1", 1500, 1600), vec!["A", "B"]);
    }

    #[test]
    fn test_transcript_db_cache_roundtrip() {
        let mut db = TranscriptDb::with_build(GenomeBuild::GRCh37);
        db.add(Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Minus,
            sequence: None,
            cds_start: Some(5),
            cds_end: Some(40),
            exons: vec![Exon::new(1, 1, 50), Exon::new(2, 51, 100)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh37,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: Some("NP_000079.2".to_string()),
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        });

        let temp = tempfile::TempDir::new().unwrap();
        let src = temp.path().join("anno.gff");
        std::fs::write(&src, b"# source").unwrap();
        let cache = TranscriptDb::cache_path_for(&src, GenomeBuild::GRCh37, AnnotationFormat::Gff3);
        // cache_path_for must use the stable Display impl ("GRCh37"), not Debug,
        // and append the format tag so caches don't collide across formats.
        assert!(
            cache.to_str().unwrap().ends_with(".GRCh37.gff3.tdb"),
            "cache path must end with .GRCh37.gff3.tdb, got {}",
            cache.display()
        );
        db.to_cache_file(&cache).unwrap();
        assert!(TranscriptDb::cache_is_fresh(&cache, &src));

        let loaded = TranscriptDb::from_cache_file(&cache).unwrap();
        assert_eq!(loaded.len(), db.len());
        assert_eq!(loaded.genome_build, GenomeBuild::GRCh37);
        let a = db.get("NM_000088.3").unwrap();
        let b = loaded.get("NM_000088.3").unwrap();
        assert_eq!(a.strand, b.strand);
        assert_eq!((a.cds_start, a.cds_end), (b.cds_start, b.cds_end));
        assert_eq!(
            (a.genomic_start, a.genomic_end),
            (b.genomic_start, b.genomic_end)
        );
        assert_eq!(a.exons.len(), b.exons.len());
        assert_eq!(a.protein_id, b.protein_id);
        // The region index is rebuilt correctly from the cache.
        assert!(!loaded.get_by_region("chr1", 1500, 1500).is_empty());

        // A magic-less / corrupt cache is rejected cleanly (not mis-read).
        std::fs::write(&cache, b"garbage-without-magic").unwrap();
        assert!(!TranscriptDb::cache_is_fresh(&cache, &src));
        assert!(TranscriptDb::from_cache_file(&cache).is_err());
    }

    #[test]
    fn cache_is_not_fresh_when_source_mtime_unreadable() {
        // A valid, current cache whose source file has gone missing must NOT be
        // reported fresh: load_annotations returns early on a fresh cache, so a
        // false positive here would silently serve stale annotations instead of
        // surfacing the missing-source error.
        let db = TranscriptDb::with_build(GenomeBuild::GRCh38);
        let temp = tempfile::TempDir::new().unwrap();
        let src = temp.path().join("anno.gff");
        std::fs::write(&src, b"# source").unwrap();
        let cache = TranscriptDb::cache_path_for(&src, GenomeBuild::GRCh38, AnnotationFormat::Gff3);
        db.to_cache_file(&cache).unwrap();

        // Cache is fresh while the source exists...
        assert!(TranscriptDb::cache_is_fresh(&cache, &src));

        // ...but not once the source is unreadable/missing.
        std::fs::remove_file(&src).unwrap();
        assert!(
            !TranscriptDb::cache_is_fresh(&cache, &src),
            "an unreadable source mtime must make the cache stale, not fresh"
        );
    }

    #[test]
    fn test_transcript_db_get_by_gene() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let transcripts = db.get_by_gene("COL1A1");
        assert_eq!(transcripts.len(), 2);
    }

    #[test]
    fn test_transcript_db_get_by_position() {
        let mut db = TranscriptDb::new();

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Position within range
        let matches = db.get_by_position("chr1", 1500);
        assert_eq!(matches.len(), 1);

        // Position outside range
        let matches = db.get_by_position("chr1", 500);
        assert_eq!(matches.len(), 0);

        // Different chromosome
        let matches = db.get_by_position("chr2", 1500);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn test_mane_select_lookup() {
        let mut db = TranscriptDb::new();

        let tx_mane = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx_other = Transcript {
            id: "NM_000089.4".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_mane);
        db.add(tx_other);

        // Should find MANE Select
        assert!(db.has_mane_select("COL1A1"));
        let mane = db.get_mane_select("COL1A1");
        assert!(mane.is_some());
        assert_eq!(mane.unwrap().id, "NM_000088.3");

        // Preferred transcript should be MANE Select
        let preferred = db.get_preferred_transcript("COL1A1");
        assert!(preferred.is_some());
        assert_eq!(preferred.unwrap().id, "NM_000088.3");
    }

    #[test]
    fn test_mane_stats() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000099.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let stats = db.mane_stats();
        assert_eq!(stats.total_genes, 2);
        assert_eq!(stats.genes_with_mane_select, 1);
        assert!((stats.mane_select_coverage() - 50.0).abs() < 0.01);
    }

    // ===== ManeStats Tests =====

    #[test]
    fn test_mane_stats_default() {
        let stats = ManeStats::default();
        assert_eq!(stats.total_genes, 0);
        assert_eq!(stats.genes_with_mane_select, 0);
        assert_eq!(stats.genes_with_mane_plus_clinical, 0);
        assert_eq!(stats.total_mane_plus_clinical_transcripts, 0);
    }

    #[test]
    fn test_mane_stats_coverage_empty() {
        let stats = ManeStats::default();
        assert_eq!(stats.mane_select_coverage(), 0.0);
    }

    #[test]
    fn test_mane_stats_coverage_full() {
        let stats = ManeStats {
            total_genes: 100,
            genes_with_mane_select: 100,
            genes_with_mane_plus_clinical: 50,
            total_mane_plus_clinical_transcripts: 60,
        };
        assert!((stats.mane_select_coverage() - 100.0).abs() < 0.01);
    }

    // ===== TranscriptDb Extended Tests =====

    #[test]
    fn test_transcript_db_with_build() {
        let db = TranscriptDb::with_build(GenomeBuild::GRCh37);
        assert!(matches!(db.genome_build, GenomeBuild::GRCh37));
    }

    #[test]
    fn test_transcript_db_get_by_region() {
        let mut db = TranscriptDb::new();

        let tx = Transcript {
            id: "NM_000088.3".to_string(),
            gene_symbol: Some("COL1A1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::default(),
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Overlapping region
        let matches = db.get_by_region("chr1", 1500, 1600);
        assert_eq!(matches.len(), 1);

        // Adjacent region (no overlap)
        let matches = db.get_by_region("chr1", 500, 999);
        assert_eq!(matches.len(), 0);

        // Fully containing region
        let matches = db.get_by_region("chr1", 900, 2100);
        assert_eq!(matches.len(), 1);

        // Different chromosome
        let matches = db.get_by_region("chr2", 1000, 2000);
        assert_eq!(matches.len(), 0);
    }

    #[test]
    fn test_transcript_db_mane_plus_clinical() {
        let mut db = TranscriptDb::new();

        let tx_plus1 = Transcript {
            id: "NM_000100.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx_plus2 = Transcript {
            id: "NM_000101.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_plus1);
        db.add(tx_plus2);

        let plus_clinical = db.get_mane_plus_clinical("BRCA1");
        assert_eq!(plus_clinical.len(), 2);

        // No MANE Select for this gene
        assert!(db.get_mane_select("BRCA1").is_none());

        // Get all MANE transcripts
        let all_mane = db.get_mane_transcripts("BRCA1");
        assert_eq!(all_mane.len(), 2);
    }

    #[test]
    fn test_transcript_db_preferred_transcript_fallback() {
        let mut db = TranscriptDb::new();

        // Non-MANE, non-coding transcript
        let tx_noncoding = Transcript {
            id: "NR_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: None,
            cds_end: None,
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        // Coding transcript (preferred over non-coding)
        let tx_coding = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_noncoding);
        db.add(tx_coding);

        // Should prefer coding transcript
        let preferred = db.get_preferred_transcript("TEST");
        assert!(preferred.is_some());
        assert!(preferred.unwrap().is_coding());
    }

    #[test]
    fn test_transcript_db_ids_iter() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("GENE1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("GENE2".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let ids: Vec<_> = db.ids().collect();
        assert_eq!(ids.len(), 2);

        let iter_count = db.iter().count();
        assert_eq!(iter_count, 2);
    }

    // ===== Gene Lookup Edge Cases =====

    #[test]
    fn test_transcript_db_nonexistent_gene() {
        let db = TranscriptDb::new();

        assert!(db.get_by_gene("NONEXISTENT").is_empty());
        assert!(db.get_mane_select("NONEXISTENT").is_none());
        assert!(db.get_mane_plus_clinical("NONEXISTENT").is_empty());
        assert!(db.get_mane_transcripts("NONEXISTENT").is_empty());
        assert!(db.get_preferred_transcript("NONEXISTENT").is_none());
        assert!(!db.has_mane_select("NONEXISTENT"));
    }

    #[test]
    fn test_transcript_db_no_gene_symbol() {
        let mut db = TranscriptDb::new();

        // Transcript without gene symbol
        let tx = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: None,
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: Some("chr1".to_string()),
            genomic_start: Some(1000),
            genomic_end: Some(2000),
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx);

        // Should still be findable by ID and position
        assert!(db.get("NM_000001.1").is_some());
        assert!(!db.get_by_position("chr1", 1500).is_empty());
    }

    // ===== Additional TranscriptDb Tests =====

    #[test]
    fn test_get_mane_transcripts() {
        let mut db = TranscriptDb::new();

        let tx1 = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        let tx2 = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("BRCA1".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx1);
        db.add(tx2);

        let mane_txs = db.get_mane_transcripts("BRCA1");
        assert_eq!(mane_txs.len(), 2); // Both MANE Select and Plus Clinical
    }

    #[test]
    fn test_get_preferred_transcript_priority() {
        let mut db = TranscriptDb::new();

        // Add transcripts in reverse priority order
        let tx_none = Transcript {
            id: "NM_000003.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::None,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };

        db.add(tx_none);

        // Without MANE, should return the coding transcript
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000003.1");

        // Add MANE Plus Clinical
        let tx_plus = Transcript {
            id: "NM_000002.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::PlusClinical,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        db.add(tx_plus);

        // Should prefer MANE Plus Clinical
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000002.1");

        // Add MANE Select
        let tx_select = Transcript {
            id: "NM_000001.1".to_string(),
            gene_symbol: Some("TEST".to_string()),
            strand: Strand::Plus,
            sequence: Some("ATGC".to_string()),
            cds_start: Some(1),
            cds_end: Some(4),
            exons: vec![Exon::new(1, 1, 4)],
            chromosome: None,
            genomic_start: None,
            genomic_end: None,
            genome_build: GenomeBuild::GRCh38,
            mane_status: ManeStatus::Select,
            refseq_match: None,
            ensembl_match: None,
            protein_id: None,
            exon_cigars: Vec::new(),
            cached_introns: OnceLock::new(),
        };
        db.add(tx_select);

        // Should prefer MANE Select
        let preferred = db.get_preferred_transcript("TEST").unwrap();
        assert_eq!(preferred.id, "NM_000001.1");
    }
}
