//! Reference data preparation for ferro-hgvs normalization.
//!
//! This module provides functionality to download and prepare reference data
//! needed for HGVS variant normalization.

use crate::FerroError;
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Legacy GenBank accessions referenced in ClinVar but not in RefSeq.
///
/// These are historical sequence records (non-RefSeq) that are still
/// referenced in ClinVar variation data.
pub const LEGACY_GENBANK_ACCESSIONS: &[&str] = &[
    "U31929.1",   // BRCA1 historical
    "M75126.1",   // Factor V Leiden
    "M22590.1",   // Beta-globin
    "AJ298105.1", // CFTR
    "AF319569.1", // MLH1
];

/// Configuration for reference data preparation.
#[derive(Debug, Clone)]
pub struct PrepareConfig {
    /// Output directory for reference data
    pub output_dir: PathBuf,
    /// Download RefSeq transcripts
    pub download_transcripts: bool,
    /// Download GRCh38 genome (large ~3GB)
    pub download_genome: bool,
    /// Download GRCh37 genome (large ~3GB) - for older ClinVar patterns
    pub download_genome_grch37: bool,
    /// Download/refresh RefSeqGene sequences (NG_* accessions, ~600MB)
    pub download_refseqgene: bool,
    /// Download LRG sequences from EBI (~1325 files)
    pub download_lrg: bool,
    /// Download cdot transcript metadata (CDS positions, exon coords, ~200MB)
    pub download_cdot: bool,
    /// Skip download if files exist
    pub skip_existing: bool,
    /// ClinVar file to extract missing accessions from
    pub clinvar_file: Option<PathBuf>,
    /// Plain text pattern file to extract missing accessions from
    pub patterns_file: Option<PathBuf>,
    /// Dry run - show what would be fetched without fetching
    pub dry_run: bool,
}

impl Default for PrepareConfig {
    fn default() -> Self {
        Self {
            output_dir: PathBuf::from("ferro-reference"),
            download_transcripts: true,
            download_genome: true,
            download_genome_grch37: false,
            download_refseqgene: false,
            download_lrg: false,
            download_cdot: true,
            skip_existing: true,
            clinvar_file: None,
            patterns_file: None,
            dry_run: false,
        }
    }
}

/// URLs for reference data downloads.
pub mod urls {
    /// RefSeq human transcript RNA sequences (multiple files)
    pub const REFSEQ_RNA_BASE: &str =
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.";
    pub const REFSEQ_RNA_SUFFIX: &str = ".rna.fna.gz";

    /// GRCh38 reference genome (with RefSeq accessions NC_000001.11, etc.)
    pub const GRCH38_GENOME: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz";

    /// GRCh37 reference genome (with RefSeq accessions NC_000001.10, etc.)
    pub const GRCH37_GENOME: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz";

    /// RefSeqGene sequences (NG_* accessions) - 9 files
    pub const REFSEQGENE_BASE: &str = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/";
    pub const REFSEQGENE_COUNT: usize = 9;

    /// cdot transcript database (contains CDS metadata, exon coordinates, etc.)
    /// From https://github.com/SACGF/cdot/releases
    pub const CDOT_REFSEQ_GRCH38: &str = "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/cdot-0.2.32.refseq.GRCh38.json.gz";
    pub const CDOT_REFSEQ_GRCH37: &str = "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/cdot-0.2.32.refseq.GRCh37.json.gz";

    /// LRG (Locus Reference Genomic) sequences from EBI
    /// ~1325 FASTA files
    pub const LRG_FASTA_BASE: &str = "https://ftp.ebi.ac.uk/pub/databases/lrgex/fasta/";
    /// LRG XML files with full annotation structure (exons, CDS, transcripts)
    pub const LRG_XML_BASE: &str = "https://ftp.ebi.ac.uk/pub/databases/lrgex/";
    pub const LRG_MAX_ID: usize = 1400; // Upper bound for LRG IDs to try

    /// LRG to RefSeq transcript mapping file
    /// Maps LRG transcript IDs (e.g., LRG_1t1) to RefSeq transcripts (e.g., NM_000088.3)
    pub const LRG_REFSEQ_MAPPING: &str =
        "https://ftp.ebi.ac.uk/pub/databases/lrgex/list_LRGs_transcripts_xrefs.txt";
}

/// Manifest of prepared reference data.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ReferenceManifest {
    /// When the data was prepared
    pub prepared_at: String,
    /// Transcript FASTA files
    pub transcript_fastas: Vec<PathBuf>,
    /// GRCh38 genome FASTA file (if downloaded)
    pub genome_fasta: Option<PathBuf>,
    /// GRCh37 genome FASTA file (if downloaded)
    #[serde(default)]
    pub genome_grch37_fasta: Option<PathBuf>,
    /// RefSeqGene FASTA files (NG_* accessions)
    #[serde(default)]
    pub refseqgene_fastas: Vec<PathBuf>,
    /// LRG FASTA files (LRG_* accessions)
    #[serde(default)]
    pub lrg_fastas: Vec<PathBuf>,
    /// LRG XML files with full annotation structure
    #[serde(default)]
    pub lrg_xmls: Vec<PathBuf>,
    /// LRG to RefSeq transcript mapping file
    #[serde(default)]
    pub lrg_refseq_mapping: Option<PathBuf>,
    /// cdot transcript metadata JSON (if downloaded)
    pub cdot_json: Option<PathBuf>,
    /// Supplemental FASTA file (missing ClinVar transcripts fetched from NCBI)
    #[serde(default)]
    pub supplemental_fasta: Option<PathBuf>,
    /// Legacy transcript versions FASTA (older versions not in current RefSeq)
    #[serde(default)]
    pub legacy_transcripts_fasta: Option<PathBuf>,
    /// Legacy transcript metadata JSON (CDS coordinates, gene names)
    #[serde(default)]
    pub legacy_transcripts_metadata: Option<PathBuf>,
    /// Legacy GenBank sequences FASTA (non-RefSeq sequences like U31929.1)
    #[serde(default)]
    pub legacy_genbank_fasta: Option<PathBuf>,
    /// Legacy GenBank metadata JSON (CDS coordinates, gene names)
    #[serde(default)]
    pub legacy_genbank_metadata: Option<PathBuf>,
    /// Total number of transcripts
    pub transcript_count: usize,
    /// List of available accession prefixes
    pub available_prefixes: Vec<String>,
}

impl ReferenceManifest {
    /// Convert all paths in the manifest to be relative to the given base directory.
    ///
    /// This ensures the manifest is portable - paths work when running from the
    /// directory containing the manifest, regardless of where `prepare` was run from.
    pub fn make_paths_relative(&mut self, base: &Path) {
        fn strip_prefix(path: &Path, base: &Path) -> PathBuf {
            path.strip_prefix(base)
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|_| path.to_path_buf())
        }

        self.transcript_fastas = self
            .transcript_fastas
            .iter()
            .map(|p| strip_prefix(p, base))
            .collect();

        if let Some(ref p) = self.genome_fasta {
            self.genome_fasta = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.genome_grch37_fasta {
            self.genome_grch37_fasta = Some(strip_prefix(p, base));
        }

        self.refseqgene_fastas = self
            .refseqgene_fastas
            .iter()
            .map(|p| strip_prefix(p, base))
            .collect();

        self.lrg_fastas = self
            .lrg_fastas
            .iter()
            .map(|p| strip_prefix(p, base))
            .collect();

        self.lrg_xmls = self
            .lrg_xmls
            .iter()
            .map(|p| strip_prefix(p, base))
            .collect();

        if let Some(ref p) = self.lrg_refseq_mapping {
            self.lrg_refseq_mapping = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.cdot_json {
            self.cdot_json = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.supplemental_fasta {
            self.supplemental_fasta = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.legacy_transcripts_fasta {
            self.legacy_transcripts_fasta = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.legacy_transcripts_metadata {
            self.legacy_transcripts_metadata = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.legacy_genbank_fasta {
            self.legacy_genbank_fasta = Some(strip_prefix(p, base));
        }

        if let Some(ref p) = self.legacy_genbank_metadata {
            self.legacy_genbank_metadata = Some(strip_prefix(p, base));
        }
    }

    /// Deduplicate paths in all path lists.
    pub fn deduplicate_paths(&mut self) {
        fn dedup_vec(paths: &mut Vec<PathBuf>) {
            let mut seen = HashSet::new();
            paths.retain(|p| seen.insert(p.clone()));
        }

        dedup_vec(&mut self.transcript_fastas);
        dedup_vec(&mut self.refseqgene_fastas);
        dedup_vec(&mut self.lrg_fastas);
        dedup_vec(&mut self.lrg_xmls);
    }
}

/// Metadata for a legacy transcript fetched from GenBank.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct LegacyTranscript {
    /// Accession with version (e.g., "NM_000051.3")
    pub id: String,
    /// Gene symbol if available
    pub gene_symbol: Option<String>,
    /// CDS start position (1-based, inclusive)
    pub cds_start: Option<u64>,
    /// CDS end position (1-based, inclusive)
    pub cds_end: Option<u64>,
    /// Total sequence length
    pub sequence_length: usize,
}

/// Metadata file for legacy transcripts.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct LegacyMetadata {
    /// When the metadata was generated
    pub generated_at: String,
    /// Map from accession to transcript metadata
    pub transcripts: HashMap<String, LegacyTranscript>,
}

/// Prepare reference data for normalization.
pub fn prepare_references(config: &PrepareConfig) -> Result<ReferenceManifest, FerroError> {
    eprintln!(
        "Preparing reference data in {}",
        config.output_dir.display()
    );

    // Create output directory
    fs::create_dir_all(&config.output_dir).map_err(|e| FerroError::Io {
        msg: format!(
            "Failed to create directory {}: {}",
            config.output_dir.display(),
            e
        ),
    })?;

    // Load existing manifest if present, otherwise start fresh
    let manifest_path = config.output_dir.join("manifest.json");
    let mut manifest = if manifest_path.exists() {
        let file = File::open(&manifest_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open existing manifest: {}", e),
        })?;
        let mut existing: ReferenceManifest =
            serde_json::from_reader(file).map_err(|e| FerroError::Io {
                msg: format!("Failed to parse existing manifest: {}", e),
            })?;
        existing.prepared_at = chrono::Utc::now().to_rfc3339();

        // Resolve relative paths in manifest against output directory
        existing.transcript_fastas = existing
            .transcript_fastas
            .into_iter()
            .map(|p| config.output_dir.join(p))
            .collect();
        existing.refseqgene_fastas = existing
            .refseqgene_fastas
            .into_iter()
            .map(|p| config.output_dir.join(p))
            .collect();
        existing.lrg_fastas = existing
            .lrg_fastas
            .into_iter()
            .map(|p| config.output_dir.join(p))
            .collect();
        existing.lrg_xmls = existing
            .lrg_xmls
            .into_iter()
            .map(|p| config.output_dir.join(p))
            .collect();
        existing.genome_fasta = existing.genome_fasta.map(|p| config.output_dir.join(p));
        existing.genome_grch37_fasta = existing
            .genome_grch37_fasta
            .map(|p| config.output_dir.join(p));
        existing.cdot_json = existing.cdot_json.map(|p| config.output_dir.join(p));
        existing.lrg_refseq_mapping = existing
            .lrg_refseq_mapping
            .map(|p| config.output_dir.join(p));
        existing.supplemental_fasta = existing
            .supplemental_fasta
            .map(|p| config.output_dir.join(p));
        existing.legacy_transcripts_fasta = existing
            .legacy_transcripts_fasta
            .map(|p| config.output_dir.join(p));

        eprintln!("  Loaded existing manifest, will merge updates");
        existing
    } else {
        ReferenceManifest {
            prepared_at: chrono::Utc::now().to_rfc3339(),
            transcript_fastas: Vec::new(),
            genome_fasta: None,
            genome_grch37_fasta: None,
            refseqgene_fastas: Vec::new(),
            lrg_fastas: Vec::new(),
            lrg_xmls: Vec::new(),
            lrg_refseq_mapping: None,
            cdot_json: None,
            supplemental_fasta: None,
            legacy_transcripts_fasta: None,
            legacy_transcripts_metadata: None,
            legacy_genbank_fasta: None,
            legacy_genbank_metadata: None,
            transcript_count: 0,
            available_prefixes: Vec::new(),
        }
    };

    // Download transcripts
    if config.download_transcripts {
        eprintln!("\n=== Downloading RefSeq transcripts ===");
        let transcript_dir = config.output_dir.join("transcripts");
        fs::create_dir_all(&transcript_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        // Download RefSeq RNA files (numbered 1-N)
        for i in 1..=20 {
            let filename = format!("human.{}.rna.fna.gz", i);
            let url = format!("{}{}{}", urls::REFSEQ_RNA_BASE, i, urls::REFSEQ_RNA_SUFFIX);
            let output_path = transcript_dir.join(&filename);

            if config.skip_existing && output_path.exists() {
                eprintln!("  Skipping {} (exists)", filename);
                manifest.transcript_fastas.push(output_path);
                continue;
            }

            match download_file(&url, &output_path) {
                Ok(_) => {
                    eprintln!("  Downloaded {}", filename);
                    manifest.transcript_fastas.push(output_path);
                }
                Err(e) => {
                    // File might not exist (numbering ends at some point)
                    if i > 1 {
                        eprintln!("  No more files after human.{}.rna.fna.gz", i - 1);
                        break;
                    }
                    return Err(e);
                }
            }
        }

        // Decompress and index
        eprintln!("\n=== Processing transcript files ===");
        for gz_path in &manifest.transcript_fastas.clone() {
            let fasta_path = gz_path.with_extension("").with_extension("fna");

            if config.skip_existing && fasta_path.exists() {
                eprintln!("  Skipping decompress {} (exists)", fasta_path.display());
            } else {
                decompress_gzip(gz_path, &fasta_path)?;
                eprintln!("  Decompressed {}", fasta_path.display());
            }

            // Index with samtools if available, otherwise use our indexer
            let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
            if config.skip_existing && fai_path.exists() {
                eprintln!("  Skipping index {} (exists)", fai_path.display());
            } else {
                index_fasta(&fasta_path)?;
                eprintln!("  Indexed {}", fasta_path.display());
            }
        }

        // Count transcripts and collect prefixes
        let (count, prefixes) = count_transcripts(&manifest.transcript_fastas)?;
        manifest.transcript_count = count;
        manifest.available_prefixes = prefixes;
        eprintln!("\n  Total transcripts: {}", count);
    }

    // Download genome
    if config.download_genome {
        eprintln!("\n=== Downloading GRCh38 genome (this may take a while) ===");
        let genome_dir = config.output_dir.join("genome");
        fs::create_dir_all(&genome_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        let gz_path = genome_dir.join("GRCh38.fna.gz");
        let fasta_path = genome_dir.join("GRCh38.fna");

        if config.skip_existing && fasta_path.exists() {
            eprintln!("  Skipping genome download (exists)");
            // Ensure index exists even for skipped files
            let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
            if !fai_path.exists() {
                eprintln!("  Indexing existing genome...");
                index_fasta(&fasta_path)?;
            }
        } else {
            download_file(urls::GRCH38_GENOME, &gz_path)?;
            decompress_gzip(&gz_path, &fasta_path)?;
            index_fasta(&fasta_path)?;
        }

        manifest.genome_fasta = Some(fasta_path);
    }

    // Download GRCh37 genome
    if config.download_genome_grch37 {
        eprintln!("\n=== Downloading GRCh37 genome (this may take a while, ~900MB) ===");
        let genome_dir = config.output_dir.join("genome");
        fs::create_dir_all(&genome_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        let gz_path = genome_dir.join("GRCh37.fna.gz");
        let fasta_path = genome_dir.join("GRCh37.fna");

        if config.skip_existing && fasta_path.exists() {
            eprintln!("  Skipping GRCh37 genome download (exists)");
            // Ensure index exists even for skipped files
            let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
            if !fai_path.exists() {
                eprintln!("  Indexing existing GRCh37 genome...");
                index_fasta(&fasta_path)?;
            }
        } else {
            eprintln!("  Downloading GRCh37 genome...");
            download_file(urls::GRCH37_GENOME, &gz_path)?;
            eprintln!("  Decompressing (~3GB)...");
            decompress_gzip(&gz_path, &fasta_path)?;
            eprintln!("  Indexing...");
            index_fasta(&fasta_path)?;
            eprintln!("  Done.");
        }

        manifest.genome_grch37_fasta = Some(fasta_path);
    }

    // Download RefSeqGene sequences
    if config.download_refseqgene {
        eprintln!("\n=== Downloading RefSeqGene sequences (NG_* accessions, ~600MB) ===");
        let refseqgene_dir = config.output_dir.join("refseqgene");
        fs::create_dir_all(&refseqgene_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        for i in 1..=urls::REFSEQGENE_COUNT {
            let filename = format!("refseqgene.{}.genomic.fna.gz", i);
            let url = format!("{}refseqgene.{}.genomic.fna.gz", urls::REFSEQGENE_BASE, i);
            let gz_path = refseqgene_dir.join(&filename);
            let fasta_path = refseqgene_dir.join(format!("refseqgene.{}.genomic.fna", i));

            if config.skip_existing && fasta_path.exists() {
                eprintln!("  Skipping {} (exists)", filename);
                // Ensure index exists even for skipped files
                let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
                if !fai_path.exists() {
                    index_fasta(&fasta_path)?;
                }
                manifest.refseqgene_fastas.push(fasta_path);
                continue;
            }

            eprintln!("  Downloading {}...", filename);
            match download_file(&url, &gz_path) {
                Ok(_) => {
                    decompress_gzip(&gz_path, &fasta_path)?;
                    index_fasta(&fasta_path)?;
                    manifest.refseqgene_fastas.push(fasta_path);
                }
                Err(e) => {
                    eprintln!("  Warning: Failed to download {}: {}", filename, e);
                }
            }
        }
        eprintln!(
            "  Downloaded {} RefSeqGene files",
            manifest.refseqgene_fastas.len()
        );
    }

    // Download LRG sequences and XML annotations from EBI
    if config.download_lrg {
        eprintln!("\n=== Preparing LRG sequences from EBI (~1325 files) ===");
        let lrg_dir = config.output_dir.join("lrg");
        fs::create_dir_all(&lrg_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        let mut fasta_downloaded = 0;
        let mut fasta_skipped = 0;
        let mut xml_downloaded = 0;
        let mut xml_skipped = 0;

        for i in 1..=urls::LRG_MAX_ID {
            let lrg_id = format!("LRG_{}", i);
            let fasta_filename = format!("{}.fasta", lrg_id);
            let xml_filename = format!("{}.xml", lrg_id);
            let fasta_url = format!("{}{}.fasta", urls::LRG_FASTA_BASE, lrg_id);
            let xml_url = format!("{}{}.xml", urls::LRG_XML_BASE, lrg_id);
            let fasta_path = lrg_dir.join(&fasta_filename);
            let xml_path = lrg_dir.join(&xml_filename);

            // Track if this LRG exists (either already downloaded or successfully downloaded)
            let mut lrg_exists = false;

            // Download/skip FASTA
            if config.skip_existing && fasta_path.exists() {
                // Ensure index exists even for skipped files
                let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
                if !fai_path.exists() {
                    index_fasta(&fasta_path)?;
                }
                fasta_skipped += 1;
                manifest.lrg_fastas.push(fasta_path.clone());
                lrg_exists = true;
            } else {
                match download_file(&fasta_url, &fasta_path) {
                    Ok(_) => {
                        index_fasta(&fasta_path)?;
                        manifest.lrg_fastas.push(fasta_path.clone());
                        fasta_downloaded += 1;
                        lrg_exists = true;
                        if fasta_downloaded % 100 == 0 {
                            eprintln!("  Downloaded {} LRG FASTA files...", fasta_downloaded);
                        }
                    }
                    Err(_) => {
                        // LRG doesn't exist (not all LRG IDs are used)
                        let _ = std::fs::remove_file(&fasta_path);
                    }
                }
            }

            // Download/skip XML only if FASTA exists
            if lrg_exists {
                if config.skip_existing && xml_path.exists() {
                    xml_skipped += 1;
                    manifest.lrg_xmls.push(xml_path);
                } else {
                    match download_file(&xml_url, &xml_path) {
                        Ok(_) => {
                            manifest.lrg_xmls.push(xml_path);
                            xml_downloaded += 1;
                            if xml_downloaded % 100 == 0 {
                                eprintln!("  Downloaded {} LRG XML files...", xml_downloaded);
                            }
                        }
                        Err(_) => {
                            // XML download failed - this shouldn't happen if FASTA exists
                            let _ = std::fs::remove_file(&xml_path);
                        }
                    }
                }
            }
        }

        if fasta_skipped > 0 {
            eprintln!("  Skipped {} existing LRG FASTA files", fasta_skipped);
        }
        if fasta_downloaded > 0 {
            eprintln!("  Downloaded {} new LRG FASTA files", fasta_downloaded);
        }
        if xml_skipped > 0 {
            eprintln!("  Skipped {} existing LRG XML files", xml_skipped);
        }
        if xml_downloaded > 0 {
            eprintln!("  Downloaded {} new LRG XML files", xml_downloaded);
        }
        eprintln!("  Indexed all LRG files");

        // Download LRG to RefSeq mapping file
        let mapping_path = lrg_dir.join("lrg_refseq_mapping.txt");
        if config.skip_existing && mapping_path.exists() {
            eprintln!("  Skipping LRG-RefSeq mapping (exists)");
        } else {
            eprintln!("  Downloading LRG-RefSeq mapping...");
            download_file(urls::LRG_REFSEQ_MAPPING, &mapping_path)?;
        }
        manifest.lrg_refseq_mapping = Some(mapping_path);
    }

    // Download cdot transcript metadata
    if config.download_cdot {
        eprintln!("\n=== Downloading cdot transcript metadata (~200MB) ===");
        let cdot_dir = config.output_dir.join("cdot");
        fs::create_dir_all(&cdot_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        let gz_path = cdot_dir.join("cdot-0.2.32.refseq.GRCh38.json.gz");
        let json_path = cdot_dir.join("cdot-0.2.32.refseq.GRCh38.json");

        if config.skip_existing && json_path.exists() {
            eprintln!("  Skipping cdot download (exists)");
        } else {
            eprintln!("  Downloading cdot-0.2.32.refseq.GRCh38.json.gz...");
            download_file(urls::CDOT_REFSEQ_GRCH38, &gz_path)?;
            eprintln!("  Decompressing...");
            decompress_gzip(&gz_path, &json_path)?;
            eprintln!("  Done.");
        }

        manifest.cdot_json = Some(json_path);
    }

    // Fetch missing transcripts from ClinVar or pattern files (requires benchmark feature)
    #[cfg(feature = "benchmark")]
    if config.clinvar_file.is_some() || config.patterns_file.is_some() {
        fetch_supplemental_data(config, &mut manifest)?;
    }

    #[cfg(not(feature = "benchmark"))]
    if config.clinvar_file.is_some() || config.patterns_file.is_some() {
        eprintln!("\n=== Warning: --clinvar and --patterns require benchmark feature ===");
        eprintln!("  Build with: cargo build --features benchmark");
        eprintln!("  Skipping supplemental data fetching.");
    }

    // Fetch legacy GenBank sequences (non-RefSeq accessions referenced in ClinVar)
    eprintln!("\n=== Fetching legacy GenBank sequences ===");
    let supplemental_dir = config.output_dir.join("supplemental");
    fs::create_dir_all(&supplemental_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create supplemental directory: {}", e),
    })?;

    let genbank_fasta = supplemental_dir.join("legacy_genbank.fna");
    if config.skip_existing && genbank_fasta.exists() {
        eprintln!("  Skipping (exists): {}", genbank_fasta.display());
        manifest.legacy_genbank_fasta = Some(genbank_fasta.clone());
        let metadata_path = genbank_fasta.with_extension("metadata.json");
        if metadata_path.exists() {
            manifest.legacy_genbank_metadata = Some(metadata_path);
        }
    } else {
        let genbank_accessions: Vec<String> = LEGACY_GENBANK_ACCESSIONS
            .iter()
            .map(|s| s.to_string())
            .collect();
        let fetched = fetch_legacy_versions(&genbank_accessions, &genbank_fasta, config.dry_run)?;

        if !config.dry_run && fetched > 0 {
            eprintln!("  Indexing legacy GenBank FASTA...");
            index_fasta(&genbank_fasta)?;
            manifest.legacy_genbank_fasta = Some(genbank_fasta.clone());
            let metadata_path = genbank_fasta.with_extension("metadata.json");
            if metadata_path.exists() {
                manifest.legacy_genbank_metadata = Some(metadata_path);
            }
        }
    }

    // Clean up and save manifest
    // Convert paths to be relative to the output directory for portability
    manifest.deduplicate_paths();
    manifest.make_paths_relative(&config.output_dir);

    let manifest_path = config.output_dir.join("manifest.json");
    let file = File::create(&manifest_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create manifest: {}", e),
    })?;
    serde_json::to_writer_pretty(file, &manifest).map_err(|e| FerroError::Io {
        msg: format!("Failed to write manifest: {}", e),
    })?;

    eprintln!("\n=== Preparation complete ===");
    eprintln!("Manifest: {}", manifest_path.display());

    Ok(manifest)
}

/// Fetch supplemental data from ClinVar/patterns (benchmark feature only)
#[cfg(feature = "benchmark")]
fn fetch_supplemental_data(
    config: &PrepareConfig,
    manifest: &mut ReferenceManifest,
) -> Result<(), FerroError> {
    eprintln!("\n=== Fetching missing transcripts ===");

    // Load existing accessions from FAI indexes
    let existing = load_existing_accessions(&config.output_dir)?;
    eprintln!(
        "  Found {} existing accessions in reference data",
        existing.len()
    );

    // Collect accessions from all sources
    let mut all_accessions = HashSet::new();

    // Extract from ClinVar TSV file
    if let Some(ref clinvar_file) = config.clinvar_file {
        let accessions = crate::benchmark::cache::extract_clinvar_accessions(clinvar_file)?;
        eprintln!("  Found {} accessions in ClinVar file", accessions.len());
        all_accessions.extend(accessions);
    }

    // Extract from plain pattern file
    if let Some(ref patterns_file) = config.patterns_file {
        let accessions = crate::benchmark::cache::extract_all_accessions_from_file(patterns_file)?;
        eprintln!("  Found {} accessions in patterns file", accessions.len());
        all_accessions.extend(accessions);
    }

    eprintln!("  {} unique accessions total", all_accessions.len());

    // Find missing accessions
    let missing: Vec<String> = all_accessions
        .into_iter()
        .filter(|acc| !existing.contains(acc))
        .collect();
    eprintln!(
        "  {} accessions are missing from reference data",
        missing.len()
    );

    if !missing.is_empty() {
        let supplemental_dir = config.output_dir.join("supplemental");
        fs::create_dir_all(&supplemental_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create supplemental directory: {}", e),
        })?;

        let supplemental_fasta = supplemental_dir.join("clinvar_transcripts.fna");

        // Fetch missing transcripts
        let fetched = crate::benchmark::cache::fetch_fasta_to_file(
            &missing,
            &supplemental_fasta,
            config.dry_run,
        )?;

        if !config.dry_run && fetched > 0 {
            // Index the supplemental FASTA
            eprintln!("  Indexing supplemental FASTA...");
            index_fasta(&supplemental_fasta)?;
            manifest.supplemental_fasta = Some(supplemental_fasta);
            eprintln!("  Fetched {} supplemental transcripts", fetched);
        }
    }

    // Detect and fetch legacy transcript versions from patterns
    if let Some(ref patterns_path) = config.patterns_file {
        eprintln!("\n=== Detecting legacy transcript versions ===");

        let transcript_dir = config.output_dir.join("transcripts");
        let supplemental_dir = config.output_dir.join("supplemental");
        let supp_dir_opt = if supplemental_dir.exists() {
            Some(supplemental_dir.as_path())
        } else {
            None
        };

        let legacy_versions = detect_legacy_versions(patterns_path, &transcript_dir, supp_dir_opt)?;

        if !legacy_versions.is_empty() {
            eprintln!("  Found {} legacy versions to fetch", legacy_versions.len());

            let supplemental_dir = config.output_dir.join("supplemental");
            fs::create_dir_all(&supplemental_dir).map_err(|e| FerroError::Io {
                msg: format!("Failed to create supplemental directory: {}", e),
            })?;

            let legacy_fasta = supplemental_dir.join("legacy_transcripts.fna");
            let fetched = fetch_legacy_versions(&legacy_versions, &legacy_fasta, config.dry_run)?;

            if !config.dry_run && fetched > 0 {
                // Index the legacy FASTA
                eprintln!("  Indexing legacy transcripts FASTA...");
                index_fasta(&legacy_fasta)?;
                manifest.legacy_transcripts_fasta = Some(legacy_fasta.clone());
                // Metadata file is created by fetch_legacy_versions
                let metadata_path = legacy_fasta.with_extension("metadata.json");
                if metadata_path.exists() {
                    manifest.legacy_transcripts_metadata = Some(metadata_path);
                }
            }
        } else {
            eprintln!("  No legacy versions needed");
        }
    }

    Ok(())
}

/// Check what reference data is available.
pub fn check_references(reference_dir: &Path) -> Result<ReferenceManifest, FerroError> {
    let manifest_path = reference_dir.join("manifest.json");

    if !manifest_path.exists() {
        return Err(FerroError::Io {
            msg: format!(
                "No reference data found at {}. Run 'ferro prepare' first.",
                reference_dir.display()
            ),
        });
    }

    let file = File::open(&manifest_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open manifest: {}", e),
    })?;

    serde_json::from_reader(file).map_err(|e| FerroError::Io {
        msg: format!("Failed to parse manifest: {}", e),
    })
}

/// Print a summary of reference data.
pub fn print_reference_summary(manifest: &ReferenceManifest, reference_dir: &Path) {
    eprintln!("=== Reference Data Summary ===");
    eprintln!("  Directory: {}", reference_dir.display());
    eprintln!("  Prepared at: {}", manifest.prepared_at);
    eprintln!("  Transcripts: {}", manifest.transcript_count);
    eprintln!(
        "  Available prefixes: {}",
        manifest.available_prefixes.join(", ")
    );

    if let Some(ref genome) = manifest.genome_fasta {
        eprintln!("  GRCh38 genome: {}", genome.display());
    }
    if let Some(ref genome) = manifest.genome_grch37_fasta {
        eprintln!("  GRCh37 genome: {}", genome.display());
    }
    if !manifest.refseqgene_fastas.is_empty() {
        eprintln!("  RefSeqGene files: {}", manifest.refseqgene_fastas.len());
    }
    if !manifest.lrg_fastas.is_empty() {
        eprintln!("  LRG files: {}", manifest.lrg_fastas.len());
    }
    if let Some(ref cdot) = manifest.cdot_json {
        eprintln!("  cdot metadata: {}", cdot.display());
    }
    if let Some(ref supp) = manifest.supplemental_fasta {
        eprintln!("  Supplemental transcripts: {}", supp.display());
    }
    if let Some(ref legacy) = manifest.legacy_transcripts_fasta {
        eprintln!("  Legacy transcripts: {}", legacy.display());
    }
    if let Some(ref genbank) = manifest.legacy_genbank_fasta {
        eprintln!("  Legacy GenBank: {}", genbank.display());
    }
}

// ============================================================================
// Helper functions
// ============================================================================

/// Download a file from a URL.
fn download_file(url: &str, output: &Path) -> Result<(), FerroError> {
    // Convert path to string, handling non-UTF-8 paths gracefully
    let output_str = output.to_str().ok_or_else(|| FerroError::Io {
        msg: format!("Path contains invalid UTF-8: {:?}", output),
    })?;

    // Use curl or wget if available (more reliable for large files)
    let curl_result = Command::new("curl")
        .args(["-fSL", "-o", output_str, url])
        .output();

    match curl_result {
        Ok(output_result) if output_result.status.success() => Ok(()),
        _ => {
            // Try wget as fallback
            let wget_result = Command::new("wget")
                .args(["-q", "-O", output_str, url])
                .output();

            match wget_result {
                Ok(output_result) if output_result.status.success() => Ok(()),
                _ => {
                    // Fall back to reqwest
                    download_with_reqwest(url, output)
                }
            }
        }
    }
}

/// Download using reqwest (for when curl/wget aren't available).
fn download_with_reqwest(url: &str, output: &Path) -> Result<(), FerroError> {
    let client = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(3600))
        .build()
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to create HTTP client: {}", e),
        })?;

    let response = client.get(url).send().map_err(|e| FerroError::Io {
        msg: format!("Failed to download {}: {}", url, e),
    })?;

    if !response.status().is_success() {
        return Err(FerroError::Io {
            msg: format!("HTTP {} for {}", response.status(), url),
        });
    }

    let mut file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;

    let content = response.bytes().map_err(|e| FerroError::Io {
        msg: format!("Failed to read response: {}", e),
    })?;

    file.write_all(&content).map_err(|e| FerroError::Io {
        msg: format!("Failed to write file: {}", e),
    })?;

    Ok(())
}

/// Decompress a gzip file.
fn decompress_gzip(input: &Path, output: &Path) -> Result<(), FerroError> {
    use flate2::read::GzDecoder;

    let input_file = File::open(input).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", input.display(), e),
    })?;

    let output_file = File::create(output).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output.display(), e),
    })?;

    let mut decoder = GzDecoder::new(BufReader::new(input_file));
    let mut writer = BufWriter::new(output_file);

    std::io::copy(&mut decoder, &mut writer).map_err(|e| FerroError::Io {
        msg: format!("Failed to decompress: {}", e),
    })?;

    Ok(())
}

/// Index a FASTA file (create .fai).
pub fn index_fasta(fasta_path: &Path) -> Result<(), FerroError> {
    // Convert path to string, handling non-UTF-8 paths gracefully
    let fasta_str = fasta_path.to_str().ok_or_else(|| FerroError::Io {
        msg: format!("FASTA path contains invalid UTF-8: {:?}", fasta_path),
    })?;

    // Try samtools first
    let samtools_result = Command::new("samtools").args(["faidx", fasta_str]).output();

    if let Ok(output) = samtools_result {
        if output.status.success() {
            return Ok(());
        }
    }

    // Fall back to our own indexer
    build_fasta_index(fasta_path)
}

/// Build a FASTA index (.fai file).
fn build_fasta_index(fasta_path: &Path) -> Result<(), FerroError> {
    let file = File::open(fasta_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open {}: {}", fasta_path.display(), e),
    })?;
    let reader = BufReader::new(file);

    let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
    let fai_file = File::create(&fai_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", fai_path.display(), e),
    })?;
    let mut writer = BufWriter::new(fai_file);

    let mut current_name: Option<String> = None;
    let mut seq_length: u64 = 0;
    let mut seq_offset: u64 = 0;
    let mut line_bases: u64 = 0;
    let mut line_bytes: u64 = 0;
    let mut byte_offset: u64 = 0;
    let mut first_seq_line = true;

    for line in reader.lines() {
        let line = line.map_err(|e| FerroError::Io {
            msg: format!("Failed to read line: {}", e),
        })?;

        if let Some(header) = line.strip_prefix('>') {
            // Write previous entry
            if let Some(name) = current_name.take() {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}",
                    name, seq_length, seq_offset, line_bases, line_bytes
                )
                .map_err(|e| FerroError::Io {
                    msg: format!("Failed to write index: {}", e),
                })?;
            }

            // Start new entry
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            current_name = Some(name);
            seq_length = 0;
            seq_offset = byte_offset + line.len() as u64 + 1; // +1 for newline
            first_seq_line = true;
        } else if current_name.is_some() {
            let bases = line.trim_end().len() as u64;
            seq_length += bases;

            if first_seq_line {
                line_bases = bases;
                line_bytes = line.len() as u64 + 1; // +1 for newline
                first_seq_line = false;
            }
        }

        byte_offset += line.len() as u64 + 1; // +1 for newline
    }

    // Write last entry
    if let Some(name) = current_name {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            name, seq_length, seq_offset, line_bases, line_bytes
        )
        .map_err(|e| FerroError::Io {
            msg: format!("Failed to write index: {}", e),
        })?;
    }

    Ok(())
}

/// Count transcripts and collect accession prefixes.
fn count_transcripts(fasta_files: &[PathBuf]) -> Result<(usize, Vec<String>), FerroError> {
    let mut count = 0usize;
    let mut prefixes = HashSet::new();

    let pb = ProgressBar::new(fasta_files.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos}/{len} {msg}")
            .unwrap(),
    );

    for fasta_path in fasta_files {
        // Use the decompressed file
        let fasta_path = fasta_path.with_extension("").with_extension("fna");
        if !fasta_path.exists() {
            continue;
        }

        let file = File::open(&fasta_path).map_err(|e| FerroError::Io {
            msg: format!("Failed to open {}: {}", fasta_path.display(), e),
        })?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(|e| FerroError::Io {
                msg: format!("Failed to read line: {}", e),
            })?;

            if let Some(header) = line.strip_prefix('>') {
                count += 1;

                // Extract accession prefix (e.g., "NM_", "NR_", "XM_")
                if let Some(acc) = header.split('|').nth(3) {
                    if let Some(prefix) = acc.split('_').next() {
                        prefixes.insert(format!("{}_", prefix));
                    }
                }
            }
        }

        pb.inc(1);
    }

    pb.finish_with_message(format!("{} transcripts", count));

    let mut prefix_vec: Vec<String> = prefixes.into_iter().collect();
    prefix_vec.sort();

    Ok((count, prefix_vec))
}

/// Load all accessions from existing FAI index files in the reference directory.
///
/// Scans all .fai files in the transcripts, genome, refseqgene, lrg, and supplemental
/// directories and returns a set of all accession IDs.
#[cfg(feature = "benchmark")]
fn load_existing_accessions(reference_dir: &Path) -> Result<HashSet<String>, FerroError> {
    let mut accessions = HashSet::new();

    // Directories to scan for FAI files
    let dirs = [
        reference_dir.join("transcripts"),
        reference_dir.join("genome"),
        reference_dir.join("refseqgene"),
        reference_dir.join("lrg"),
        reference_dir.join("supplemental"),
    ];

    for dir in &dirs {
        if !dir.exists() {
            continue;
        }

        // Find all .fai files
        let entries = fs::read_dir(dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to read directory {}: {}", dir.display(), e),
        })?;

        for entry in entries.flatten() {
            let path = entry.path();
            if path.extension().and_then(|e| e.to_str()) == Some("fai") {
                // Read FAI file and extract accession names
                let file = File::open(&path).map_err(|e| FerroError::Io {
                    msg: format!("Failed to open {}: {}", path.display(), e),
                })?;
                let reader = BufReader::new(file);

                for line in reader.lines().map_while(Result::ok) {
                    // FAI format: name\tlength\toffset\tline_bases\tline_bytes
                    if let Some(name) = line.split('\t').next() {
                        accessions.insert(name.to_string());
                    }
                }
            }
        }
    }

    Ok(accessions)
}

/// Load all versioned accessions from FAI index files.
///
/// Returns a set of accession.version strings (e.g., "NM_000051.4").
fn load_available_versions(
    transcript_dir: &Path,
    supplemental_dir: Option<&Path>,
) -> HashSet<String> {
    let mut versions = HashSet::new();

    // Load from transcript FAI files
    if let Ok(entries) = fs::read_dir(transcript_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.extension().is_some_and(|e| e == "fai") {
                if let Ok(file) = File::open(&path) {
                    for line in BufReader::new(file).lines().map_while(Result::ok) {
                        if let Some(acc) = line.split('\t').next() {
                            versions.insert(acc.to_string());
                        }
                    }
                }
            }
        }
    }

    // Load from supplemental FAI files if present
    if let Some(supp_dir) = supplemental_dir {
        if let Ok(entries) = fs::read_dir(supp_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.extension().is_some_and(|e| e == "fai") {
                    if let Ok(file) = File::open(&path) {
                        for line in BufReader::new(file).lines().map_while(Result::ok) {
                            if let Some(acc) = line.split('\t').next() {
                                versions.insert(acc.to_string());
                            }
                        }
                    }
                }
            }
        }
    }

    versions
}

/// Extract versioned accessions from patterns file.
///
/// Filters to only transcript accessions (NM_, NR_, XM_, XR_) with version numbers.
fn extract_versioned_accessions_from_patterns(
    patterns_path: &Path,
) -> Result<HashSet<String>, FerroError> {
    let file = File::open(patterns_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to open patterns file: {}", e),
    })?;

    let mut accessions = HashSet::new();

    for line in BufReader::new(file).lines().map_while(Result::ok) {
        // Extract accession from HGVS pattern (e.g., "NM_000051.3:c.123A>G" -> "NM_000051.3")
        if let Some(colon_pos) = line.find(':') {
            let accession = &line[..colon_pos];
            // Check if it's a versioned transcript accession
            let is_transcript = accession.starts_with("NM_")
                || accession.starts_with("NR_")
                || accession.starts_with("XM_")
                || accession.starts_with("XR_");
            let has_version = accession.contains('.');

            if is_transcript && has_version {
                accessions.insert(accession.to_string());
            }
        }
    }

    Ok(accessions)
}

/// Detect legacy transcript versions from patterns that aren't in current RefSeq.
///
/// A "legacy version" is an older version of a transcript (e.g., NM_000051.3)
/// when the current RefSeq only has a newer version (e.g., NM_000051.4).
pub fn detect_legacy_versions(
    patterns_path: &Path,
    transcript_dir: &Path,
    supplemental_dir: Option<&Path>,
) -> Result<Vec<String>, FerroError> {
    // Load all available versions from FAI files
    let available = load_available_versions(transcript_dir, supplemental_dir);
    eprintln!("  Loaded {} available transcript versions", available.len());

    // Extract versioned accessions from patterns
    let pattern_versions = extract_versioned_accessions_from_patterns(patterns_path)?;
    eprintln!(
        "  Found {} versioned accessions in patterns",
        pattern_versions.len()
    );

    // Find versions in patterns that aren't available
    let mut legacy: Vec<String> = pattern_versions
        .into_iter()
        .filter(|acc| {
            // Only consider transcript accessions
            (acc.starts_with("NM_")
                || acc.starts_with("NR_")
                || acc.starts_with("XM_")
                || acc.starts_with("XR_"))
                &&
            // Not in available versions
            !available.contains(acc)
        })
        .collect();

    legacy.sort();
    Ok(legacy)
}

/// Parse a GenBank record to extract sequence and CDS info.
///
/// Returns (sequence, cds_start, cds_end, gene_name) where coordinates are 1-based inclusive.
fn parse_genbank_record(genbank: &str) -> Option<(String, Option<u64>, Option<u64>, String)> {
    let mut sequence = String::new();
    let mut in_origin = false;
    let mut gene_name = String::new();
    let mut cds_start: Option<u64> = None;
    let mut cds_end: Option<u64> = None;

    for line in genbank.lines() {
        // Extract sequence from ORIGIN section
        if line.starts_with("ORIGIN") {
            in_origin = true;
            continue;
        }
        if in_origin {
            if line.starts_with("//") {
                break;
            }
            // Parse sequence lines: "   123 acgtacgt acgtacgt..."
            let seq_part: String = line
                .chars()
                .filter(|c| c.is_alphabetic())
                .collect::<String>()
                .to_uppercase();
            sequence.push_str(&seq_part);
        }

        // Extract gene name from /gene qualifier
        if line.contains("/gene=\"") {
            if let Some(start) = line.find("/gene=\"") {
                if let Some(end) = line[start + 7..].find('"') {
                    gene_name = line[start + 7..start + 7 + end].to_string();
                }
            }
        }

        // Extract CDS coordinates
        // Format: "     CDS             123..456" or "     CDS             join(123..200,300..456)"
        let trimmed = line.trim();
        if trimmed.starts_with("CDS") && !trimmed.contains('/') {
            let coords = trimmed.trim_start_matches("CDS").trim();
            // Handle complement() wrapper
            let coords = coords
                .trim_start_matches("complement(")
                .trim_end_matches(')');
            // Handle join() for multi-exon
            if coords.starts_with("join(") {
                let inner = coords.trim_start_matches("join(").trim_end_matches(')');
                // Get first and last coordinates
                let parts: Vec<&str> = inner.split(',').collect();
                if let Some(first) = parts.first() {
                    if let Some((start, _)) = first.split_once("..") {
                        cds_start = start.trim().trim_start_matches('<').parse().ok();
                    }
                }
                if let Some(last) = parts.last() {
                    if let Some((_, end)) = last.split_once("..") {
                        cds_end = end.trim().trim_end_matches('>').parse().ok();
                    }
                }
            } else if let Some((start, end)) = coords.split_once("..") {
                cds_start = start.trim().trim_start_matches('<').parse().ok();
                cds_end = end.trim().trim_end_matches('>').parse().ok();
            }
        }
    }

    if sequence.is_empty() {
        return None;
    }

    Some((sequence, cds_start, cds_end, gene_name))
}

/// Fetch legacy transcript versions from NCBI and save to FASTA with metadata.
///
/// These are older RefSeq versions that are referenced in patterns
/// but not in the current RefSeq release. Fetches GenBank format to
/// extract CDS coordinates and gene names.
pub fn fetch_legacy_versions(
    legacy_accessions: &[String],
    output_fasta: &Path,
    dry_run: bool,
) -> Result<usize, FerroError> {
    if legacy_accessions.is_empty() {
        return Ok(0);
    }

    if dry_run {
        eprintln!(
            "  [dry run] Would fetch {} legacy versions",
            legacy_accessions.len()
        );
        for acc in legacy_accessions.iter().take(10) {
            eprintln!("    {}", acc);
        }
        if legacy_accessions.len() > 10 {
            eprintln!("    ... and {} more", legacy_accessions.len() - 10);
        }
        return Ok(0);
    }

    eprintln!(
        "  Fetching {} legacy transcript versions (GenBank format for CDS metadata)...",
        legacy_accessions.len()
    );

    let mut fasta_file = File::create(output_fasta).map_err(|e| FerroError::Io {
        msg: format!("Failed to create {}: {}", output_fasta.display(), e),
    })?;

    // Create metadata file path (same as FASTA but with .metadata.json extension)
    let metadata_path = output_fasta.with_extension("metadata.json");
    let mut metadata = LegacyMetadata {
        generated_at: chrono::Utc::now().to_rfc3339(),
        transcripts: HashMap::new(),
    };

    let mut fetched = 0;
    let mut failed = Vec::new();

    // Create progress bar
    let pb = ProgressBar::new(legacy_accessions.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("    [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("##-"),
    );

    for acc in legacy_accessions {
        pb.set_message(acc.clone());

        // Fetch GenBank format (not FASTA) to get CDS coordinates
        let url = format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=gb&retmode=text",
            acc
        );

        let output = Command::new("curl")
            .args(["-s", "-L", &url])
            .output()
            .map_err(|e| FerroError::Io {
                msg: format!("Failed to fetch {}: {}", acc, e),
            })?;

        let genbank_text = String::from_utf8_lossy(&output.stdout);

        // Parse GenBank record
        if output.status.success() && genbank_text.contains("LOCUS") {
            if let Some((seq, cds_start, cds_end, gene_name)) = parse_genbank_record(&genbank_text)
            {
                // Write FASTA format
                let header = if gene_name.is_empty() {
                    acc.clone()
                } else {
                    format!("{} {}", acc, gene_name)
                };
                writeln!(fasta_file, ">{}", header).map_err(|e| FerroError::Io {
                    msg: format!("Failed to write {}: {}", acc, e),
                })?;
                // Write sequence in 70-character lines
                for chunk in seq.as_bytes().chunks(70) {
                    writeln!(fasta_file, "{}", std::str::from_utf8(chunk).unwrap_or("")).map_err(
                        |e| FerroError::Io {
                            msg: format!("Failed to write sequence: {}", e),
                        },
                    )?;
                }

                // Add to metadata
                metadata.transcripts.insert(
                    acc.clone(),
                    LegacyTranscript {
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

                fetched += 1;
            } else {
                failed.push(acc.clone());
            }
        } else {
            failed.push(acc.clone());
        }

        pb.inc(1);

        // Rate limit to avoid NCBI throttling
        std::thread::sleep(std::time::Duration::from_millis(100));
    }

    pb.finish_and_clear();

    // Write metadata JSON file
    let metadata_file = File::create(&metadata_path).map_err(|e| FerroError::Io {
        msg: format!("Failed to create metadata file: {}", e),
    })?;
    serde_json::to_writer_pretty(metadata_file, &metadata).map_err(|e| FerroError::Io {
        msg: format!("Failed to write metadata: {}", e),
    })?;

    if !failed.is_empty() {
        eprintln!("  Warning: Failed to fetch {} accessions:", failed.len());
        for acc in failed.iter().take(5) {
            eprintln!("    {}", acc);
        }
        if failed.len() > 5 {
            eprintln!("    ... and {} more", failed.len() - 5);
        }
    }

    eprintln!("  Fetched {} legacy transcripts with CDS metadata", fetched);
    eprintln!("  Metadata written to: {}", metadata_path.display());

    Ok(fetched)
}

/// Detect patterns file from ferro reference directory.
/// Returns the first .txt file found in {ferro_ref}/patterns/
pub fn detect_patterns_from_ferro_reference(ferro_ref: &Path) -> Option<PathBuf> {
    let patterns_dir = ferro_ref.join("patterns");
    if !patterns_dir.exists() {
        return None;
    }

    std::fs::read_dir(&patterns_dir)
        .ok()?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| p.extension().map(|ext| ext == "txt").unwrap_or(false))
}

/// Detect ClinVar file from ferro reference directory.
/// Returns the first .gz file found in {ferro_ref}/clinvar/
pub fn detect_clinvar_from_ferro_reference(ferro_ref: &Path) -> Option<PathBuf> {
    let clinvar_dir = ferro_ref.join("clinvar");
    if !clinvar_dir.exists() {
        return None;
    }

    std::fs::read_dir(&clinvar_dir)
        .ok()?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| p.extension().map(|ext| ext == "gz").unwrap_or(false))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_build_fasta_index() {
        let dir = TempDir::new().unwrap();
        let fasta_path = dir.path().join("test.fna");

        // Create a simple FASTA file
        let mut f = File::create(&fasta_path).unwrap();
        writeln!(f, ">seq1 description").unwrap();
        writeln!(f, "ATGCATGC").unwrap();
        writeln!(f, "ATGCATGC").unwrap();
        writeln!(f, ">seq2").unwrap();
        writeln!(f, "GGGGCCCC").unwrap();

        build_fasta_index(&fasta_path).unwrap();

        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        assert!(fai_path.exists());

        let content = std::fs::read_to_string(&fai_path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("seq1\t16\t")); // 16 bases total
        assert!(lines[1].starts_with("seq2\t8\t")); // 8 bases
    }

    #[test]
    fn test_prepare_config_default() {
        let config = PrepareConfig::default();
        assert_eq!(config.output_dir, PathBuf::from("ferro-reference"));
        assert!(config.download_transcripts);
        assert!(config.download_genome);
        assert!(!config.download_genome_grch37);
        assert!(!config.download_refseqgene);
        assert!(!config.download_lrg);
        assert!(config.download_cdot);
        assert!(config.skip_existing);
    }

    #[test]
    fn test_detect_patterns_from_ferro_reference_with_patterns() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create patterns directory with a .txt file
        let patterns_dir = ferro_ref.join("patterns");
        std::fs::create_dir_all(&patterns_dir).unwrap();
        let patterns_file = patterns_dir.join("clinvar_patterns.txt");
        std::fs::write(
            &patterns_file,
            "NM_000001.1:c.100A>G\nNM_000002.2:c.200del\n",
        )
        .unwrap();

        // Should detect the patterns file
        let detected = detect_patterns_from_ferro_reference(ferro_ref);
        assert!(detected.is_some());
        assert_eq!(detected.unwrap(), patterns_file);
    }

    #[test]
    fn test_detect_patterns_from_ferro_reference_no_patterns_dir() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // No patterns directory exists
        let detected = detect_patterns_from_ferro_reference(ferro_ref);
        assert!(detected.is_none());
    }

    #[test]
    fn test_detect_patterns_from_ferro_reference_empty_patterns_dir() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create empty patterns directory
        let patterns_dir = ferro_ref.join("patterns");
        std::fs::create_dir_all(&patterns_dir).unwrap();

        // Should not detect any patterns file
        let detected = detect_patterns_from_ferro_reference(ferro_ref);
        assert!(detected.is_none());
    }

    #[test]
    fn test_detect_patterns_from_ferro_reference_ignores_non_txt() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create patterns directory with non-.txt files
        let patterns_dir = ferro_ref.join("patterns");
        std::fs::create_dir_all(&patterns_dir).unwrap();
        std::fs::write(patterns_dir.join("patterns.json"), "{}").unwrap();
        std::fs::write(patterns_dir.join("readme.md"), "# Readme").unwrap();

        // Should not detect non-.txt files
        let detected = detect_patterns_from_ferro_reference(ferro_ref);
        assert!(detected.is_none());
    }

    #[test]
    fn test_detect_patterns_finds_first_txt_file() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create patterns directory with multiple .txt files
        let patterns_dir = ferro_ref.join("patterns");
        std::fs::create_dir_all(&patterns_dir).unwrap();
        std::fs::write(patterns_dir.join("a_patterns.txt"), "pattern1").unwrap();
        std::fs::write(patterns_dir.join("b_patterns.txt"), "pattern2").unwrap();

        // Should detect a .txt file (order may vary due to filesystem)
        let detected = detect_patterns_from_ferro_reference(ferro_ref);
        assert!(detected.is_some());
        let path = detected.unwrap();
        assert!(path.extension().map(|e| e == "txt").unwrap_or(false));
    }

    #[test]
    fn test_detect_clinvar_from_ferro_reference_with_clinvar() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create clinvar directory with a .gz file
        let clinvar_dir = ferro_ref.join("clinvar");
        std::fs::create_dir_all(&clinvar_dir).unwrap();
        let clinvar_file = clinvar_dir.join("hgvs4variation.txt.gz");
        std::fs::write(&clinvar_file, "mock clinvar data").unwrap();

        // Should detect the clinvar file
        let detected = detect_clinvar_from_ferro_reference(ferro_ref);
        assert!(detected.is_some());
        assert_eq!(detected.unwrap(), clinvar_file);
    }

    #[test]
    fn test_detect_clinvar_from_ferro_reference_no_clinvar_dir() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // No clinvar directory exists
        let detected = detect_clinvar_from_ferro_reference(ferro_ref);
        assert!(detected.is_none());
    }

    #[test]
    fn test_detect_clinvar_from_ferro_reference_empty_clinvar_dir() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create empty clinvar directory
        let clinvar_dir = ferro_ref.join("clinvar");
        std::fs::create_dir_all(&clinvar_dir).unwrap();

        // Should not detect any clinvar file
        let detected = detect_clinvar_from_ferro_reference(ferro_ref);
        assert!(detected.is_none());
    }

    #[test]
    fn test_detect_clinvar_from_ferro_reference_ignores_non_gz() {
        let temp_dir = TempDir::new().unwrap();
        let ferro_ref = temp_dir.path();

        // Create clinvar directory with non-.gz files
        let clinvar_dir = ferro_ref.join("clinvar");
        std::fs::create_dir_all(&clinvar_dir).unwrap();
        std::fs::write(clinvar_dir.join("data.txt"), "plain text").unwrap();
        std::fs::write(clinvar_dir.join("readme.md"), "# Readme").unwrap();

        // Should not detect non-.gz files
        let detected = detect_clinvar_from_ferro_reference(ferro_ref);
        assert!(detected.is_none());
    }
}
