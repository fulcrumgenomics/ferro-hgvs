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

pub mod backfill;
pub mod manifest;
pub use manifest::{check_references, ReferenceManifest};

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
    /// Download companion RefSeq protein FASTAs (`*.protein.faa.gz`, ~same size
    /// as the transcript set) for the translated-CDS-vs-canonical-protein check
    /// (issue #520). Off by default — opt in only when you need the check.
    pub download_proteins: bool,
    /// Download GRCh38 genome (large ~3GB)
    pub download_genome: bool,
    /// Download GRCh37 genome (large ~3GB) - for older ClinVar patterns
    pub download_genome_grch37: bool,
    /// Download/refresh RefSeqGene sequences (NG_* accessions, ~600MB)
    pub download_refseqgene: bool,
    /// Download LRG sequences from EBI (~1325 files)
    pub download_lrg: bool,
    /// Download GRCh38 cdot transcript metadata (CDS positions, exon coords, ~200MB)
    pub download_cdot: bool,
    /// Download GRCh37 cdot transcript metadata (CDS positions, exon coords, ~200MB)
    pub download_cdot_grch37: bool,
    /// Download Ensembl cDNA transcript sequences (ENST_ accessions, ~75MB) and
    /// the Ensembl cdot metadata. Off by default — opt in with `--ensembl`.
    pub download_ensembl: bool,
    /// Skip download if files exist
    pub skip_existing: bool,
    /// ClinVar file to extract missing accessions from
    pub clinvar_file: Option<PathBuf>,
    /// Plain text pattern file to extract missing accessions from
    pub patterns_file: Option<PathBuf>,
    /// Newline-delimited file of exact versioned accessions to validate against
    /// their authoritative GenBank records, writing a canonical-overrides file
    /// (issue #520). `None` skips the canonical-validation stage.
    pub validate_canonical_accessions: Option<PathBuf>,
    /// Optional newline-delimited NG_ accession file; when set, prepare derives
    /// `derived_refseqgene_placements.json` and wires the manifest field.
    pub derive_ng_placements: Option<PathBuf>,
    /// Optional newline-delimited `accession.version` file (#842); when set,
    /// prepare backfills each listed transcript that is present in cdot but
    /// absent from the bulk RefSeq RNA FASTA — fetching its deposited sequence
    /// from NCBI and appending it to `backfill/backfill_transcripts.fna`, then
    /// wiring the manifest field. Requires cdot in the same prepare run.
    /// Networked; per-accession failures warn and continue.
    pub backfill_transcripts: Option<PathBuf>,
    /// Genome build selection string ("grch38", "grch37", "all", "none").
    /// Stored so `builds_for_genome` can map it to build names at derivation time.
    pub genome: String,
    /// Dry run - show what would be fetched without fetching
    pub dry_run: bool,
}

impl Default for PrepareConfig {
    fn default() -> Self {
        Self {
            output_dir: PathBuf::from("ferro-reference"),
            download_transcripts: true,
            download_proteins: false,
            download_genome: true,
            download_genome_grch37: false,
            download_refseqgene: false,
            download_lrg: false,
            download_cdot: true,
            download_cdot_grch37: false,
            download_ensembl: false,
            skip_existing: true,
            clinvar_file: None,
            patterns_file: None,
            validate_canonical_accessions: None,
            derive_ng_placements: None,
            backfill_transcripts: None,
            genome: "grch38".to_string(),
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
    pub const REFSEQ_RNA_COUNT: usize = 20;

    /// GRCh38 reference genome (with RefSeq accessions NC_000001.11, etc.)
    pub const GRCH38_GENOME: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz";

    /// GRCh37 reference genome (with RefSeq accessions NC_000001.10, etc.)
    pub const GRCH37_GENOME: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz";

    /// GRCh38 assembly report (RefSeq-accession → assembly map, #716).
    pub const GRCH38_ASSEMBLY_REPORT: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt";
    /// GRCh37 assembly report (#716).
    pub const GRCH37_ASSEMBLY_REPORT: &str = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt";

    /// RefSeqGene sequences (NG_* accessions) - 9 files
    pub const REFSEQGENE_BASE: &str = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/";
    pub const REFSEQGENE_COUNT: usize = 9;

    /// RefSeqGene→genome alignment GFF3 (NG_ chromosomal placements, #480).
    ///
    /// NCBI stopped *updating* the RefSeqGene→genome alignments in 2024 (the
    /// live `H_sapiens/alignments` dir now carries only known/model RefSeq
    /// transcript alignments, and the `RefSeqGene/` dir's GFF3 symlinks are
    /// dangling). The latest **GRCh38** RefSeqGene alignment is the archived
    /// RefSeq-release-109 (GRCh38.p13) snapshot under `alignments/ARCHIVE/all/`.
    /// Primary-chromosome `NC_` accessions are unchanged across GRCh38 patches,
    /// so it remains valid against the `.40` (p14) genome/cdot; `NG_` records
    /// curated after 2021 are simply absent (they decline — never a wrong
    /// coordinate). Revisit if NCBI republishes a maintained feed.
    pub const REFSEQGENE_ALIGNMENTS_FILE: &str =
        "GCF_000001405.39_109.20211119_refseqgene_alignments.gff3";
    pub const REFSEQGENE_ALIGNMENTS_URL: &str =
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/alignments/ARCHIVE/all/GCF_000001405.39_109.20211119_refseqgene_alignments.gff3";

    /// GRCh37 RefSeqGene→genome alignment GFF3 (#653/#713).
    ///
    /// The final **GRCh37** RefSeqGene alignment is the archived RefSeq
    /// annotation-release-105 snapshot (assembly `GCF_000001405.25` = GRCh37.p13,
    /// dated 2020-10-22) under `alignments/ARCHIVE/all/`. NCBI froze the GRCh37
    /// feed, so this is the latest; `NG_` records curated later are simply absent
    /// (they decline — never a wrong coordinate). Merged with the GRCh38 file
    /// (above) so an `NG_` resolves to its build-appropriate placement.
    pub const REFSEQGENE_ALIGNMENTS_GRCH37_FILE: &str =
        "GCF_000001405.25_105.20201022_refseqgene_alignments.gff3";
    pub const REFSEQGENE_ALIGNMENTS_GRCH37_URL: &str =
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/alignments/ARCHIVE/all/GCF_000001405.25_105.20201022_refseqgene_alignments.gff3";

    /// NCBI `LRG_RefSeqGene` association table — maps each gene to its
    /// reference-standard transcript, used to resolve legacy gene-model
    /// selectors (`NG_(GENE_v001):c.…`) to the transcript accession (#500/#637).
    /// A small plain (un-gzipped) TSV alongside the RefSeqGene FASTAs.
    pub const REFSEQGENE_SUMMARY_FILE: &str = "LRG_RefSeqGene";
    pub const REFSEQGENE_SUMMARY_URL: &str =
        "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene";

    /// cdot transcript database (contains CDS metadata, exon coordinates, etc.)
    /// From https://github.com/SACGF/cdot/releases
    pub const CDOT_REFSEQ_GRCH38: &str = "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/cdot-0.2.32.refseq.GRCh38.json.gz";
    pub const CDOT_REFSEQ_GRCH37: &str = "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/cdot-0.2.32.refseq.GRCh37.json.gz";

    /// Ensembl cdot transcript database (ENST/ENSG/ENSP), same SACGF/cdot
    /// release as the RefSeq files above and the identical JSON schema.
    pub const CDOT_ENSEMBL_GRCH38: &str = "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/cdot-0.2.32.ensembl.GRCh38.json.gz";
    pub const CDOT_ENSEMBL_GRCH37: &str = "https://github.com/SACGF/cdot/releases/download/data_v0.2.32/cdot-0.2.32.ensembl.GRCh37.json.gz";

    /// Ensembl cDNA (transcript) sequences (ENST_ accessions). Release 110 is
    /// contemporaneous with cdot `data_v0.2.32`; both GRCh38 and GRCh37 mirrors
    /// share the same file name under their respective release trees.
    pub const ENSEMBL_CDNA_GRCH38: &str = "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz";
    /// GRCh37 cDNA mirror. Intentionally unused by the prepare flow: a
    /// transcript's cDNA sequence is build-independent (only its genomic mapping
    /// differs between GRCh37 and GRCh38), so the GRCh38 download above already
    /// covers both builds. Kept here as the documented source URL should a
    /// build-specific cDNA set ever be needed.
    pub const ENSEMBL_CDNA_GRCH37: &str = "https://ftp.ensembl.org/pub/grch37/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz";

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

    // Pre-flight: deriving NG_ placements needs cdot + a genome build in the
    // same run (the provider reconstructs NC_ exon structure from cdot). Reject
    // a cdot-less selection (e.g. `--genome none` or `--no-cdot`, where neither
    // `download_cdot*` flag is set) up front, before any downloads, instead of
    // running the whole download budget and failing in the late guard in
    // `derive_placements_for_accessions`.
    if config.derive_ng_placements.is_some()
        && !(config.download_cdot || config.download_cdot_grch37)
    {
        return Err(FerroError::Io {
            msg: "--derive-ng-placements requires cdot + a genome build in the same \
                  prepare run; the current selection (e.g. --genome none or --no-cdot) \
                  downloads no cdot metadata, so NG_ placements cannot be derived"
                .to_string(),
        });
    }

    // Create output directory
    fs::create_dir_all(&config.output_dir).map_err(|e| FerroError::Io {
        msg: format!(
            "Failed to create directory {}: {}",
            config.output_dir.display(),
            e
        ),
    })?;

    // Load existing manifest if present, else default
    let mut manifest = ReferenceManifest::load_or_default(&config.output_dir)?;

    // Pre-flight: the version-aware backfill (#842) targets the cdot-vs-FASTA
    // gap, so it needs cdot to know which versions are worth fetching — either
    // downloaded in this run, or already on disk from a prior prepare (so a
    // user can backfill a few accessions against an existing reference without
    // re-downloading cdot).
    //
    // Reject stale referenced cdot paths first. `load_cdot_versions` unions the
    // versions from *every* cdot build the manifest references, so the backfill's
    // "is this version known to cdot?" universe is only complete if all referenced
    // builds load. A manifest that references both GRCh38 and GRCh37 cdot but whose
    // GRCh37 file was deleted/moved would still pass an any-one-exists guard
    // (GRCh38 is present), yet `load_cdot_versions` would only warn and silently
    // compute an incomplete universe — leaving GRCh37-only versions like
    // NM_002001.2 (#802) wrongly classified `not_in_cdot`. So require every
    // referenced cdot path that is not being redownloaded this run to exist on
    // disk, failing up front before any downloads instead of burning the budget
    // and failing later in `backfill_transcripts_step` (before `manifest.save()`).
    if config.backfill_transcripts.is_some() {
        let stale_cdot_paths: Vec<&PathBuf> = [
            (manifest.cdot_json.as_ref(), config.download_cdot),
            (
                manifest.cdot_grch37_json.as_ref(),
                config.download_cdot_grch37,
            ),
        ]
        .into_iter()
        .filter_map(|(path, refreshed_this_run)| match path {
            Some(path) if !refreshed_this_run && !path.exists() => Some(path),
            _ => None,
        })
        .collect();
        if !stale_cdot_paths.is_empty() {
            return Err(FerroError::Io {
                msg: format!(
                    "--backfill-transcripts requires every referenced cdot file to exist \
                     (or be redownloaded this run); missing: {}",
                    stale_cdot_paths
                        .iter()
                        .map(|p| p.display().to_string())
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            });
        }

        // Then reject a truly cdot-less selection (nothing downloaded this run, no
        // cdot referenced on disk) — the stale-path check above already cleared any
        // referenced-but-missing path, so an on-disk reference here is loadable.
        let existing_cdot_on_disk = manifest.cdot_json.as_ref().is_some_and(|p| p.exists())
            || manifest
                .cdot_grch37_json
                .as_ref()
                .is_some_and(|p| p.exists());
        if !(config.download_cdot || config.download_cdot_grch37 || existing_cdot_on_disk) {
            return Err(FerroError::Io {
                msg: "--backfill-transcripts requires cdot, either downloaded in this run or \
                      already present in the reference; the current selection (e.g. \
                      --genome none or --no-cdot against a reference with no cdot) has none, \
                      so the cdot-vs-FASTA version gap cannot be computed"
                    .to_string(),
            });
        }
    }

    // Download transcripts
    if config.download_transcripts {
        eprintln!("\n=== Downloading RefSeq transcripts ===");
        let transcript_dir = config.output_dir.join("transcripts");
        fs::create_dir_all(&transcript_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        // Download RefSeq RNA files (numbered 1-N)
        for i in 1..=urls::REFSEQ_RNA_COUNT {
            let filename = format!("human.{}.rna.fna.gz", i);
            let url = format!("{}{}{}", urls::REFSEQ_RNA_BASE, i, urls::REFSEQ_RNA_SUFFIX);
            let output_path = transcript_dir.join(&filename);

            if config.skip_existing && output_path.exists() {
                eprintln!("  Skipping {} (exists)", filename);
                // Treat skipped pre-existing files as discovered so downstream
                // steps (count_transcripts, manifest persistence) see them.
                if !manifest.transcript_fastas.contains(&output_path) {
                    manifest.transcript_fastas.push(output_path);
                }
                continue;
            }

            match download_file(&url, &output_path) {
                Ok(_) => {
                    eprintln!("  Downloaded {}", filename);
                    if !manifest.transcript_fastas.contains(&output_path) {
                        manifest.transcript_fastas.push(output_path);
                    }
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
        for gz_path in &manifest.transcript_fastas {
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

    // Companion RefSeq protein FASTAs (`*.protein.faa.gz`) for the
    // translated-CDS-vs-canonical-protein check (issue #520, edge 3). Opt-in:
    // the protein set is roughly the size of the RNA set, so it is NOT pulled
    // by a default `ferro prepare`. Same `mRNA_Prot/human.N.` base as the RNA
    // files. Discover the available `.gz` files first (break on the first 404
    // past #1), then decompress + index + record each by its `.faa` path.
    if config.download_proteins {
        eprintln!("\n=== Downloading RefSeq proteins ===");
        let protein_dir = config.output_dir.join("transcripts");
        fs::create_dir_all(&protein_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;

        let mut gz_paths: Vec<PathBuf> = Vec::new();
        for i in 1..=urls::REFSEQ_RNA_COUNT {
            let gz_name = format!("human.{}.protein.faa.gz", i);
            let gz_path = protein_dir.join(&gz_name);
            if config.skip_existing && gz_path.exists() {
                eprintln!("  Skipping {} (exists)", gz_name);
                gz_paths.push(gz_path);
                continue;
            }
            let url = format!("{}{}.protein.faa.gz", urls::REFSEQ_RNA_BASE, i);
            match download_file(&url, &gz_path) {
                Ok(_) => {
                    eprintln!("  Downloaded {}", gz_name);
                    gz_paths.push(gz_path);
                }
                Err(e) => {
                    if i > 1 {
                        eprintln!(
                            "  No more protein files after human.{}.protein.faa.gz",
                            i - 1
                        );
                        break;
                    }
                    return Err(e);
                }
            }
        }

        for gz_path in &gz_paths {
            let faa_path = gz_path.with_extension("").with_extension("faa");
            if !(config.skip_existing && faa_path.exists()) {
                decompress_gzip(gz_path, &faa_path)?;
            }
            let fai_path = PathBuf::from(format!("{}.fai", faa_path.display()));
            if !(config.skip_existing && fai_path.exists()) {
                index_fasta(&faa_path)?;
            }
            if !manifest.protein_fastas.contains(&faa_path) {
                manifest.protein_fastas.push(faa_path);
            }
        }
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

        // Only overwrite the recorded path when a new one is produced: a
        // transient download soft-fail returns `None`, and assigning that would
        // erase a previously valid path and persist a degraded manifest.
        if let Some(path) = record_assembly_report(
            &genome_dir,
            "GRCh38.assembly_report.txt",
            urls::GRCH38_ASSEMBLY_REPORT,
            config.skip_existing,
        )? {
            manifest.assembly_report = Some(path);
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

        // Preserve a previously valid path on a transient soft-fail (`None`);
        // see the GRCh38 call above.
        if let Some(path) = record_assembly_report(
            &genome_dir,
            "GRCh37.assembly_report.txt",
            urls::GRCH37_ASSEMBLY_REPORT,
            config.skip_existing,
        )? {
            manifest.assembly_report_grch37 = Some(path);
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

        // RefSeqGene→genome alignment GFF3s — each NG_ record's chromosomal
        // placement, used to project transcript coordinates into an NG_ parent's
        // own frame (#480). Plain (un-gzipped) GFF3. Two single-build snapshots
        // are fetched and merged downstream so an NG_ resolves to its
        // build-appropriate placement: GRCh38 (release 109) and GRCh37 (release
        // 105, #653/#713). A failed/absent file is non-fatal — that build's
        // placements are simply unavailable (NG_ inputs on it decline).
        for (name, url, is_grch37) in [
            (
                urls::REFSEQGENE_ALIGNMENTS_FILE,
                urls::REFSEQGENE_ALIGNMENTS_URL,
                false,
            ),
            (
                urls::REFSEQGENE_ALIGNMENTS_GRCH37_FILE,
                urls::REFSEQGENE_ALIGNMENTS_GRCH37_URL,
                true,
            ),
        ] {
            let aln_path = refseqgene_dir.join(name);
            let record = |manifest: &mut ReferenceManifest, path: PathBuf| {
                if is_grch37 {
                    manifest.refseqgene_alignments_grch37 = Some(path);
                } else {
                    manifest.refseqgene_alignments = Some(path);
                }
            };
            if config.skip_existing && aln_path.exists() {
                eprintln!("  Skipping {} (exists)", name);
                record(&mut manifest, aln_path);
            } else {
                eprintln!("  Downloading {}...", name);
                match download_file(url, &aln_path) {
                    Ok(_) => record(&mut manifest, aln_path),
                    Err(e) => eprintln!("  Warning: Failed to download {}: {}", name, e),
                }
            }
        }

        // `LRG_RefSeqGene` association table — gene → reference-standard
        // transcript, used to resolve legacy gene-model selectors (#500/#637).
        let summary_name = urls::REFSEQGENE_SUMMARY_FILE;
        let summary_path = refseqgene_dir.join(summary_name);
        if config.skip_existing && summary_path.exists() {
            eprintln!("  Skipping {} (exists)", summary_name);
            manifest.refseqgene_summary = Some(summary_path);
        } else {
            eprintln!("  Downloading {}...", summary_name);
            match download_file(urls::REFSEQGENE_SUMMARY_URL, &summary_path) {
                Ok(_) => manifest.refseqgene_summary = Some(summary_path),
                Err(e) => eprintln!("  Warning: Failed to download {}: {}", summary_name, e),
            }
        }
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
        eprintln!("\n=== Downloading GRCh38 cdot transcript metadata (~200MB) ===");
        let cdot_dir = config.output_dir.join("cdot");
        manifest.cdot_json = Some(download_cdot(
            urls::CDOT_REFSEQ_GRCH38,
            &cdot_dir,
            config.skip_existing,
        )?);
    }

    // Download GRCh37 cdot transcript metadata
    if config.download_cdot_grch37 {
        eprintln!("\n=== Downloading GRCh37 cdot transcript metadata (~200MB) ===");
        let cdot_dir = config.output_dir.join("cdot");
        manifest.cdot_grch37_json = Some(download_cdot(
            urls::CDOT_REFSEQ_GRCH37,
            &cdot_dir,
            config.skip_existing,
        )?);
    }

    // Download Ensembl reference (cDNA sequences + cdot metadata). Opt-in via
    // `--ensembl`: lets variants on Ensembl transcripts/genes (ENST/ENSG/ENSP)
    // resolve instead of no-op'ing.
    if config.download_ensembl {
        let ensembl_dir = config.output_dir.join("ensembl");
        fs::create_dir_all(&ensembl_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create directory: {}", e),
        })?;
        let cdot_dir = config.output_dir.join("cdot");

        // 1. Ensembl cDNA FASTA (ENST_ sequences). Decompress + index so the
        //    provider's directory scan picks it up like the RefSeq transcripts.
        eprintln!("\n=== Downloading Ensembl cDNA transcripts (GRCh38, ~75MB) ===");
        let gz_path = ensembl_dir.join("Homo_sapiens.GRCh38.cdna.all.fa.gz");
        if config.skip_existing && gz_path.exists() {
            eprintln!("  Skipping {} (exists)", gz_path.display());
        } else {
            download_file(urls::ENSEMBL_CDNA_GRCH38, &gz_path)?;
            eprintln!("  Downloaded {}", gz_path.display());
        }
        // `with_extension("")` strips only the `.gz`, leaving the `.fa`.
        let fasta_path = gz_path.with_extension("");
        if config.skip_existing && fasta_path.exists() {
            eprintln!("  Skipping decompress {} (exists)", fasta_path.display());
        } else {
            decompress_gzip(&gz_path, &fasta_path)?;
            eprintln!("  Decompressed {}", fasta_path.display());
        }
        let fai_path = PathBuf::from(format!("{}.fai", fasta_path.display()));
        if config.skip_existing && fai_path.exists() {
            eprintln!("  Skipping index {} (exists)", fai_path.display());
        } else {
            index_fasta(&fasta_path)?;
            eprintln!("  Indexed {}", fasta_path.display());
        }
        if !manifest.ensembl_transcript_fastas.contains(&gz_path) {
            manifest.ensembl_transcript_fastas.push(gz_path);
        }

        // 2. Ensembl cdot metadata (GRCh38 primary; GRCh37 secondary if asked).
        eprintln!("\n=== Downloading Ensembl cdot transcript metadata (GRCh38) ===");
        manifest.ensembl_cdot_json = Some(download_cdot(
            urls::CDOT_ENSEMBL_GRCH38,
            &cdot_dir,
            config.skip_existing,
        )?);
        if config.download_cdot_grch37 {
            eprintln!("\n=== Downloading Ensembl cdot transcript metadata (GRCh37) ===");
            manifest.ensembl_cdot_grch37_json = Some(download_cdot(
                urls::CDOT_ENSEMBL_GRCH37,
                &cdot_dir,
                config.skip_existing,
            )?);
        }
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

    // Canonical-record validation: fetch authoritative GenBank records for the
    // requested accessions and write the canonical-overrides file (issue #520).
    if let Some(ref acc_file) = config.validate_canonical_accessions {
        eprintln!("\n=== Validating against authoritative GenBank records ===");
        let accessions = read_accession_list(acc_file)?;
        let out = config.output_dir.join("canonical_overrides.json");
        let written = fetch_canonical_overrides(&accessions, &out, config.dry_run)?;
        if !config.dry_run && written > 0 {
            manifest.canonical_overrides = Some(out.clone());
            eprintln!(
                "  Wrote {} canonical override records to {}",
                written,
                out.display()
            );
        }
    }

    // Version-aware transcript backfill (#842): fetch deposited sequences for
    // requested accession.versions that are present in cdot but absent from the
    // bulk transcript FASTA, append them to the backfill FASTA, and wire the
    // manifest field so the provider serves them from the primary index path
    // (synthesis is never invoked for those accessions). Networked; opt-in.
    if let Some(acc_file) = config.backfill_transcripts.clone() {
        backfill_transcripts_step(config, &mut manifest, &acc_file)?;
    }

    manifest.save()?;

    // Derive version-independent NG_ placements (#728/#740) when requested.
    // The manifest is already saved above, so the provider can load cdot +
    // genome from disk. Honors skip_existing: an existing output file is left
    // in place unless --force.
    if let Some(acc_file) = &config.derive_ng_placements {
        if config.dry_run {
            eprintln!(
                "  [dry-run] would derive NG_ placements from {}",
                acc_file.display()
            );
        } else {
            let out_path = config.output_dir.join("derived_refseqgene_placements.json");
            let hosted_out_path = config.output_dir.join("ng_hosted_transcripts.json");
            if config.skip_existing && out_path.exists() {
                // Caveat: if this is a pre-#792 reference (derived_refseqgene_placements.json
                // exists but ng_hosted_transcripts.json does not), manifest.ng_hosted_transcripts
                // will remain unset. Re-run with --force to regenerate both artifacts.
                eprintln!("  Skipping NG_ placement derivation (exists)");
                manifest.derived_refseqgene_placements = Some(out_path);
                if hosted_out_path.exists() {
                    manifest.ng_hosted_transcripts = Some(hosted_out_path);
                } else {
                    eprintln!(
                        "  Note: ng_hosted_transcripts.json is absent (pre-#792 reference); \
                         re-run with --force to derive it."
                    );
                }
                manifest.save()?;
            } else {
                eprintln!("  Deriving NG_ placements from {}...", acc_file.display());
                let accessions = read_accession_list(acc_file)?;
                let builds = builds_for_genome(&config.genome);
                let manifest_path = config.output_dir.join("manifest.json");
                let provider = crate::reference::multi_fasta::MultiFastaProvider::from_manifest(
                    &manifest_path,
                )?;
                let description = format!(
                    "Derived NG_/LRG_ placements (#728) produced by `ferro prepare \
                     --derive-ng-placements {}`.",
                    acc_file.display()
                );
                let derive_out =
                    crate::reference::ng_placement_builder::derive_placements_for_accessions(
                        &provider,
                        &accessions,
                        &builds,
                        description,
                    )?;
                std::fs::write(&out_path, derive_out.artifact.to_json()?).map_err(|e| {
                    FerroError::Io {
                        msg: format!("write {}: {e}", out_path.display()),
                    }
                })?;
                eprintln!(
                    "  Derived {} placement(s), {} declined.",
                    derive_out.artifact.placements.len(),
                    derive_out.declined_count
                );
                manifest.derived_refseqgene_placements = Some(out_path);
                // Co-emit ng_hosted_transcripts.json from the same GenBank fetch (#792).
                let hosted_path =
                    write_ng_hosted_transcripts(&derive_out.hosted, &hosted_out_path)?;
                manifest.ng_hosted_transcripts = Some(hosted_path);
                manifest.save()?;
            }
        }
    }

    eprintln!("\n=== Preparation complete ===");
    eprintln!(
        "Manifest: {}",
        manifest.reference_dir.join("manifest.json").display()
    );
    Ok(manifest)
}

/// Genome-build names to derive placements for, from the `--genome` selection.
fn builds_for_genome(genome: &str) -> Vec<String> {
    match genome {
        "grch38" => vec!["GRCh38".to_string()],
        "grch37" => vec!["GRCh37".to_string()],
        "all" => vec!["GRCh38".to_string(), "GRCh37".to_string()],
        // `none` (and any unknown selection) provisions no genome build, so it
        // has no builds to derive placements for. The up-front guard in
        // `prepare_references` already rejects `--derive-ng-placements` with such
        // a selection; returning empty keeps this consistent if reached otherwise.
        _ => Vec::new(),
    }
}

/// Build [`NgHostedTranscripts`](crate::reference::ng_hosted_transcripts::NgHostedTranscripts)
/// from `(NG_acc.version, [(gene, transcript_id)])` records and write it as
/// JSON to `out`, returning the written path (#792).
fn write_ng_hosted_transcripts(
    records: &[(String, Vec<(String, String)>)],
    out: &Path,
) -> Result<PathBuf, FerroError> {
    use crate::reference::ng_hosted_transcripts::NgHostedTranscripts;
    let nh = NgHostedTranscripts::from_records(records.iter().cloned());
    let json = nh.to_json().map_err(|e| FerroError::Io {
        msg: format!("serialize ng_hosted_transcripts: {e}"),
    })?;
    std::fs::write(out, json).map_err(|e| FerroError::Io {
        msg: format!("write {}: {e}", out.display()),
    })?;
    Ok(out.to_path_buf())
}

/// Load the union of `accession.version` keys present in the manifest's cdot
/// JSON file(s) (GRCh38 primary + GRCh37 secondary, whichever are present).
///
/// Used by the version-aware backfill (#842) to know which requested versions
/// cdot actually carries an alignment for. A load failure for one file warns and
/// is skipped (the other file's keys still count).
fn load_cdot_versions(manifest: &ReferenceManifest) -> HashSet<String> {
    let mut versions = HashSet::new();
    for path in [&manifest.cdot_json, &manifest.cdot_grch37_json]
        .into_iter()
        .flatten()
    {
        // Use `load` (not `from_json_file`) so the rkyv cache built by
        // `download_cdot` is reused instead of re-parsing the ~200MB JSON.
        match crate::data::cdot::CdotMapper::load(path) {
            // `all_build_transcript_ids` (not `transcript_ids`) so a version
            // present only on a non-primary build is still seen as in-cdot.
            // `CdotMapper::load` always uses GRCh38 as the primary build, so
            // when it loads the GRCh37 cdot file standalone every transcript —
            // nested under `genome_builds["GRCh37"]` — lands in the alt map and
            // the primary `transcripts` map (all `transcript_ids` reads) is
            // empty. An old version GRCh38 dropped but GRCh37 retains (e.g.
            // NM_002001.2, #802) therefore lives only in the alt map, so
            // `transcript_ids` misses it and it is wrongly classified
            // `not_in_cdot`; `all_build_transcript_ids` surfaces it.
            Ok(mapper) => versions.extend(mapper.all_build_transcript_ids().map(str::to_string)),
            Err(e) => eprintln!(
                "  Warning: failed to load cdot {} for backfill diff: {}",
                path.display(),
                e
            ),
        }
    }
    versions
}

/// Version-aware transcript backfill (#842).
///
/// Reads the requested `accession.version` list, diffs it against the cdot keys
/// and the transcript-FASTA keys, fetches the genuine gap (in cdot, absent from
/// the FASTA) from NCBI, appends it to `backfill/backfill_transcripts.fna`,
/// reindexes, and wires the manifest field. Per-accession fetch failures warn and
/// continue; the prepare run is never aborted by them.
fn backfill_transcripts_step(
    config: &PrepareConfig,
    manifest: &mut ReferenceManifest,
    acc_file: &Path,
) -> Result<(), FerroError> {
    use backfill::{execute_backfill, plan_backfill};

    eprintln!("\n=== Version-aware transcript backfill (#842) ===");
    let requested = read_accession_list(acc_file)?;
    eprintln!("  {} requested accession.version(s)", requested.len());

    let backfill_dir = config.output_dir.join("backfill");
    let backfill_fasta = backfill_dir.join("backfill_transcripts.fna");

    if config.dry_run {
        eprintln!(
            "  [dry-run] would diff {} requested accessions against cdot + the \
             transcript FASTA and backfill the gap into {}",
            requested.len(),
            backfill_fasta.display()
        );
        return Ok(());
    }

    // cdot keys (the alignment-bearing universe worth backfilling).
    let cdot_versions = load_cdot_versions(manifest);
    if cdot_versions.is_empty() {
        return Err(FerroError::Io {
            msg: "no cdot transcripts could be loaded for --backfill-transcripts; \
                  refresh or re-download cdot in this reference before retrying"
                .to_string(),
        });
    }

    // `--force` (skip_existing == false): discard any prior backfill output so
    // every requested accession is re-fetched, and so the append path never
    // produces a duplicate `>accession` record (which would make `.fai` keys
    // collide). In the default incremental mode the prior output is kept and
    // folded into the present set below as the cache.
    let force = !config.skip_existing;
    if force {
        let fai = PathBuf::from(format!("{}.fai", backfill_fasta.display()));
        let _ = fs::remove_file(&backfill_fasta);
        let _ = fs::remove_file(&fai);
    }

    // Present keys: bulk transcript FASTA + supplemental + (incremental only)
    // any prior backfill output, so a re-run skips already-fetched accessions —
    // the incremental cache. The backfill present set is read from the FASTA
    // headers themselves (the file `execute_backfill` appends to), not its
    // derived `.fai`: a prior run interrupted after the append but before
    // reindex — or a deleted `.fai` — would otherwise look empty and re-fetch
    // and re-append the same accession, duplicating its `>accession` record
    // (and colliding `.fai` keys on the subsequent reindex).
    let transcripts_dir = config.output_dir.join("transcripts");
    let supplemental_dir = config.output_dir.join("supplemental");
    let mut present = load_available_versions(&transcripts_dir, Some(&supplemental_dir));
    if !force {
        present.extend(load_fasta_header_accessions(&backfill_fasta));
    }

    let plan = plan_backfill(&requested, &cdot_versions, &present);
    eprintln!(
        "  {} to fetch, {} already present, {} not in cdot",
        plan.to_fetch.len(),
        plan.already_present.len(),
        plan.not_in_cdot.len()
    );
    if !plan.not_in_cdot.is_empty() {
        eprintln!(
            "  Warning: {} requested accession(s) are absent from cdot and will not be \
             backfilled (no exon alignment): {}",
            plan.not_in_cdot.len(),
            plan.not_in_cdot.join(", ")
        );
    }

    if !plan.to_fetch.is_empty() {
        fs::create_dir_all(&backfill_dir).map_err(|e| FerroError::Io {
            msg: format!("Failed to create backfill directory: {}", e),
        })?;
        eprintln!(
            "  Fetching {} deposited transcript sequence(s) from NCBI...",
            plan.to_fetch.len()
        );
        let outcome = execute_backfill(&plan.to_fetch, &backfill_fasta, fetch_transcript_fasta)?;
        eprintln!(
            "  Backfilled {} sequence(s), {} failed",
            outcome.fetched.len(),
            outcome.failed.len()
        );
        if !outcome.failed.is_empty() {
            eprintln!(
                "  Warning: {} accession(s) could not be backfilled (fetch/parse/version \
                 mismatch): {}",
                outcome.failed.len(),
                outcome.failed.join(", ")
            );
        }
    }

    // Wire the manifest field iff a non-empty backfill FASTA exists (newly
    // written this run, or carried over from a prior incremental run); otherwise
    // clear it so a deleted/empty file never leaves a dangling manifest
    // reference. Mirrors the no-orphan rule in fetch_canonical_overrides.
    let has_content = backfill_fasta
        .metadata()
        .map(|m| m.len() > 0)
        .unwrap_or(false);
    if has_content {
        eprintln!("  Indexing backfill FASTA...");
        index_fasta(&backfill_fasta)?;
        manifest.backfill_transcripts_fasta = Some(backfill_fasta);
    } else {
        manifest.backfill_transcripts_fasta = None;
    }

    Ok(())
}

/// Fetch a transcript's deposited sequence from NCBI EFetch as FASTA text.
///
/// Returns the raw FASTA record (`>{accession.version} desc\n<bases>`) on success,
/// or `None` on a fetch failure. `rettype=fasta` (not `gb`): the version-exact CDS
/// frame comes from cdot, so only the bases are needed here. Same 100 ms NCBI
/// rate-limit as the other efetch callers. The injected boundary for
/// [`backfill::execute_backfill`]; tests supply fixtures instead.
fn fetch_transcript_fasta(accession: &str) -> Option<String> {
    let url = format!(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text",
        accession
    );
    // `--fail` turns an HTTP 4xx/5xx into a non-zero exit (treated as a failed
    // fetch below); the timeouts bound a single slow/stalled request so one
    // accession cannot hang the whole prepare run.
    let result = Command::new("curl")
        .args([
            "-sS",
            "-L",
            "--fail",
            "--connect-timeout",
            "10",
            "--max-time",
            "60",
            &url,
        ])
        .output();
    std::thread::sleep(std::time::Duration::from_millis(100)); // NCBI rate limit
    let out = result.ok()?;
    let text = String::from_utf8_lossy(&out.stdout).into_owned();
    (out.status.success() && text.starts_with('>')).then_some(text)
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

// ============================================================================
// ============================================================================
// cdot helpers
// ============================================================================

/// Download and decompress a cdot JSON.gz file from `url` into `cdot_dir`,
/// then serialize to bincode for fast subsequent loading.
/// Returns the path to the decompressed JSON file.
fn download_cdot(url: &str, cdot_dir: &Path, skip_existing: bool) -> Result<PathBuf, FerroError> {
    fs::create_dir_all(cdot_dir).map_err(|e| FerroError::Io {
        msg: format!("Failed to create directory: {}", e),
    })?;

    let gz_filename = Path::new(url).file_name().ok_or_else(|| FerroError::Io {
        msg: format!("cdot URL has no filename component: {}", url),
    })?;
    let gz_path = cdot_dir.join(gz_filename);
    let json_path = gz_path.with_extension("");

    let json_is_fresh = if skip_existing && json_path.exists() {
        eprintln!("  Skipping cdot download (exists)");
        false
    } else {
        if skip_existing && gz_path.exists() {
            eprintln!(
                "  Reusing existing {}...",
                gz_path.file_name().unwrap_or_default().to_string_lossy()
            );
        } else {
            eprintln!(
                "  Downloading {}...",
                gz_path.file_name().unwrap_or_default().to_string_lossy()
            );
            download_file(url, &gz_path)?;
        }
        eprintln!("  Decompressing...");
        decompress_gzip(&gz_path, &json_path)?;
        eprintln!("  Done.");
        true
    };

    // Parse JSON and serialize to an rkyv archive for fast subsequent loading.
    // Only skip when an existing cache is actually *current* — a stale cache
    // (written by a build with a different layout) fails to load and silently
    // forces the slow JSON parse on every run, so it must be regenerated.
    let rkyv_path = json_path.with_extension("rkyv");
    let rkyv_is_current =
        rkyv_path.exists() && crate::data::cdot::CdotMapper::rkyv_is_current(&rkyv_path);
    if !json_is_fresh && skip_existing && rkyv_is_current {
        eprintln!("  Skipping cdot cache conversion (up to date)");
    } else {
        if rkyv_path.exists() && !rkyv_is_current {
            eprintln!("  Existing cdot cache is stale; regenerating...");
        }
        eprintln!("  Converting cdot JSON to rkyv archive for fast loading...");
        let mapper = crate::data::cdot::CdotMapper::from_json_file(&json_path).map_err(|e| {
            FerroError::Io {
                msg: format!("Failed to parse cdot JSON: {}", e),
            }
        })?;
        mapper.to_rkyv_file(&rkyv_path)?;
        eprintln!("  Done.");
    }

    Ok(json_path)
}

// Helper functions
// ============================================================================

/// Download an NCBI assembly report into `genome_dir/<filename>` (unless it
/// already exists and `skip_existing`), returning the path to record in the
/// manifest. A download failure is a soft warning that yields `Ok(None)` — the
/// report is strictly additive (build inference falls back to the hardcoded
/// table), so it must never abort the genome download (#716).
fn record_assembly_report(
    genome_dir: &Path,
    filename: &str,
    url: &str,
    skip_existing: bool,
) -> Result<Option<PathBuf>, FerroError> {
    let dest = genome_dir.join(filename);
    if skip_existing && dest.exists() {
        return Ok(Some(dest));
    }
    match download_file(url, &dest) {
        Ok(()) => Ok(Some(dest)),
        Err(e) => {
            eprintln!(
                "  Warning: failed to download assembly report {url}: {e}; \
                 build inference will use the hardcoded fallback"
            );
            Ok(None)
        }
    }
}

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

/// Read the set of record accessions directly from a FASTA file's headers
/// (`>{accession} description` → `accession`, the first whitespace-delimited
/// token after `>`).
///
/// Unlike [`load_available_versions`], which scans the derived `.fai`, this
/// reads the `.fna` itself — the source of truth for the incremental backfill
/// cache, since [`backfill::execute_backfill`] appends to the `.fna`. A prior
/// run interrupted after the append but before reindex, or a deleted `.fai`,
/// would leave the index empty (or stale) and cause the same accession to be
/// re-fetched and appended a second time. Returns an empty set if the file is
/// absent or unreadable.
fn load_fasta_header_accessions(fasta: &Path) -> HashSet<String> {
    let mut versions = HashSet::new();
    if let Ok(file) = File::open(fasta) {
        for line in BufReader::new(file).lines().map_while(Result::ok) {
            if let Some(rest) = line.strip_prefix('>') {
                if let Some(acc) = rest.split_whitespace().next() {
                    versions.insert(acc.to_string());
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

/// Read a newline-delimited accession list, trimming whitespace and skipping
/// blank lines and `#` comments.
fn read_accession_list(path: &Path) -> Result<Vec<String>, FerroError> {
    let text = std::fs::read_to_string(path).map_err(|e| FerroError::Io {
        msg: format!("Failed to read accession list {}: {}", path.display(), e),
    })?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        // Upper-case so a lowercase entry (`nm_012459.2`) doesn't false-reject:
        // efetch is case-insensitive and returns the canonical upper-case
        // `VERSION`, which the version-match guard compares against exactly.
        .map(str::to_uppercase)
        .collect())
}

/// Fetch authoritative GenBank records for `accessions` and write a
/// canonical-overrides JSON to `output`, returning the number of records
/// written.
///
/// Reuses the same NCBI efetch (`rettype=gb`) + 100 ms rate-limit pattern as
/// [`fetch_legacy_versions`]; record parsing and assembly (including the
/// requested-vs-fetched version check) are delegated to
/// [`build_canonical_overrides`](crate::reference::authoritative::build_canonical_overrides).
pub fn fetch_canonical_overrides(
    accessions: &[String],
    output: &Path,
    dry_run: bool,
) -> Result<usize, FerroError> {
    use crate::reference::authoritative::build_canonical_overrides;

    if accessions.is_empty() {
        return Ok(0);
    }
    if dry_run {
        eprintln!(
            "  [dry run] Would validate {} accessions against authoritative GenBank records",
            accessions.len()
        );
        return Ok(0);
    }

    eprintln!(
        "  Fetching {} authoritative GenBank records for canonical validation...",
        accessions.len()
    );
    let pb = ProgressBar::new(accessions.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("    [{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} {msg}")
            .unwrap()
            .progress_chars("##-"),
    );

    let (overrides, failed) = build_canonical_overrides(accessions, |acc| {
        pb.set_message(acc.to_string());
        let url = format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=gb&retmode=text",
            acc
        );
        let result = Command::new("curl").args(["-s", "-L", &url]).output();
        pb.inc(1);
        std::thread::sleep(std::time::Duration::from_millis(100)); // NCBI rate limit
        let out = result.ok()?;
        let text = String::from_utf8_lossy(&out.stdout).into_owned();
        (out.status.success() && text.contains("LOCUS")).then_some(text)
    });
    pb.finish_and_clear();

    // Don't leave an orphan empty file (and an unreferenced manifest entry) when
    // every accession failed; that mirrors the empty-input path which writes
    // nothing. Still report the failures below.
    if !overrides.is_empty() {
        let json = overrides.to_json().map_err(|e| FerroError::Io {
            msg: format!("Failed to serialize canonical overrides: {e}"),
        })?;
        std::fs::write(output, json).map_err(|e| FerroError::Io {
            msg: format!("Failed to write {}: {}", output.display(), e),
        })?;
    }

    if !failed.is_empty() {
        eprintln!(
            "  Warning: {} accessions could not be validated (fetch/parse/version mismatch):",
            failed.len()
        );
        for acc in failed.iter().take(5) {
            eprintln!("    {}", acc);
        }
        if failed.len() > 5 {
            eprintln!("    ... and {} more", failed.len() - 5);
        }
    }
    Ok(overrides.len())
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
        assert!(!config.download_cdot_grch37);
        assert!(config.skip_existing);
    }

    #[test]
    fn test_backfill_transcripts_step_errors_when_no_cdot_loaded() {
        // #842: an explicit `--backfill-transcripts` run against a reference
        // whose cdot cannot be loaded (empty version universe) must fail loudly
        // rather than succeed as a warning-only no-op — otherwise a stale or
        // corrupt manifest silently produces no backfill output.
        let temp_dir = TempDir::new().unwrap();
        let acc_file = temp_dir.path().join("accessions.txt");
        std::fs::write(&acc_file, "NM_000088.3\n").unwrap();

        let config = PrepareConfig {
            output_dir: temp_dir.path().to_path_buf(),
            dry_run: false,
            ..PrepareConfig::default()
        };
        // Default manifest carries no cdot paths, so `load_cdot_versions`
        // returns an empty set.
        let mut manifest = ReferenceManifest::default();

        let err = backfill_transcripts_step(&config, &mut manifest, &acc_file)
            .expect_err("empty cdot must abort the backfill step");
        match err {
            FerroError::Io { msg } => {
                assert!(
                    msg.contains("no cdot transcripts could be loaded"),
                    "unexpected error message: {msg}"
                );
            }
            other => panic!("expected FerroError::Io, got {other:?}"),
        }
    }

    #[test]
    fn test_load_fasta_header_accessions_reads_unindexed_fna() {
        // #842 regression: the incremental backfill cache must reflect what is
        // actually in backfill_transcripts.fna even when its `.fai` is missing
        // (e.g. a prior run was interrupted after the append but before reindex,
        // or the index was deleted). Reading only the `.fai` would treat the
        // file as empty and re-append the same accession, duplicating the
        // `>accession` record.
        let dir = TempDir::new().unwrap();
        let fasta = dir.path().join("backfill_transcripts.fna");
        let mut f = File::create(&fasta).unwrap();
        writeln!(f, ">NM_000088.3 Homo sapiens collagen").unwrap();
        writeln!(f, "ACGTACGT").unwrap();
        writeln!(f, ">NM_000001.1 some description").unwrap();
        writeln!(f, "GGGGCCCC").unwrap();
        drop(f);
        // Deliberately no `.fai` written alongside the `.fna`.

        let present = load_fasta_header_accessions(&fasta);
        assert_eq!(present.len(), 2);
        assert!(present.contains("NM_000088.3"));
        assert!(present.contains("NM_000001.1"));
    }

    #[test]
    fn test_load_fasta_header_accessions_absent_file_is_empty() {
        let dir = TempDir::new().unwrap();
        let present = load_fasta_header_accessions(&dir.path().join("missing.fna"));
        assert!(present.is_empty());
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

    #[test]
    fn prepare_config_default_has_no_derive_ng() {
        let cfg = PrepareConfig::default();
        assert!(cfg.derive_ng_placements.is_none());
    }

    #[test]
    fn prepare_config_default_has_no_backfill() {
        let cfg = PrepareConfig::default();
        assert!(cfg.backfill_transcripts.is_none());
    }

    #[test]
    fn backfill_requires_cdot_in_same_run() {
        // A cdot-less selection (no GRCh38/GRCh37 cdot) must be rejected up front,
        // before any downloads, mirroring the --derive-ng-placements guard.
        let tmp = tempfile::tempdir().unwrap();
        let acc_file = tmp.path().join("targets.txt");
        std::fs::write(&acc_file, "NM_002001.2\n").unwrap();
        let config = PrepareConfig {
            output_dir: tmp.path().join("ref"),
            download_transcripts: false,
            download_genome: false,
            download_cdot: false,
            download_cdot_grch37: false,
            backfill_transcripts: Some(acc_file),
            ..Default::default()
        };
        let err = prepare_references(&config).expect_err("cdot-less backfill must be rejected");
        assert!(
            format!("{err}").contains("--backfill-transcripts requires cdot"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn backfill_rejects_stale_manifest_cdot_path() {
        // A manifest that *references* cdot but whose cdot file is no longer on
        // disk (deleted/moved) must be rejected up front, just like a cdot-less
        // selection — `is_some()` would pass the guard and let the run burn the
        // download budget only to fail later in backfill_transcripts_step before
        // the manifest is persisted.
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path().join("ref");
        std::fs::create_dir_all(&ref_dir).unwrap();
        // Persist a manifest whose cdot_json points at a file that does not exist.
        let mut manifest = ReferenceManifest::load_or_default(&ref_dir).unwrap();
        manifest.cdot_json = Some(ref_dir.join("cdot.json"));
        manifest.save().unwrap();
        assert!(!ref_dir.join("cdot.json").exists());

        let acc_file = tmp.path().join("targets.txt");
        std::fs::write(&acc_file, "NM_002001.2\n").unwrap();
        let config = PrepareConfig {
            output_dir: ref_dir,
            download_transcripts: false,
            download_genome: false,
            download_cdot: false,
            download_cdot_grch37: false,
            backfill_transcripts: Some(acc_file),
            ..Default::default()
        };
        let err =
            prepare_references(&config).expect_err("stale-cdot-path backfill must be rejected");
        assert!(
            format!("{err}").contains("requires every referenced cdot file to exist"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn backfill_rejects_partial_stale_cdot_when_other_build_present() {
        // A manifest that references *both* cdot builds but whose GRCh37 file was
        // deleted/moved must be rejected even though the GRCh38 file is present.
        // An any-one-exists guard would pass on GRCh38, but `load_cdot_versions`
        // unions every referenced build, so a missing GRCh37 would only warn and
        // compute an incomplete universe — wrongly classifying GRCh37-only
        // versions (e.g. NM_002001.2, #802) as `not_in_cdot`.
        let tmp = tempfile::tempdir().unwrap();
        let ref_dir = tmp.path().join("ref");
        std::fs::create_dir_all(&ref_dir).unwrap();
        // GRCh38 cdot present on disk; GRCh37 cdot referenced but missing.
        let cdot_grch38 = ref_dir.join("cdot_grch38.json");
        std::fs::write(&cdot_grch38, b"").unwrap();
        let mut manifest = ReferenceManifest::load_or_default(&ref_dir).unwrap();
        manifest.cdot_json = Some(cdot_grch38.clone());
        manifest.cdot_grch37_json = Some(ref_dir.join("cdot_grch37.json"));
        manifest.save().unwrap();
        assert!(cdot_grch38.exists());
        assert!(!ref_dir.join("cdot_grch37.json").exists());

        let acc_file = tmp.path().join("targets.txt");
        std::fs::write(&acc_file, "NM_002001.2\n").unwrap();
        let config = PrepareConfig {
            output_dir: ref_dir.clone(),
            download_transcripts: false,
            download_genome: false,
            download_cdot: false,
            download_cdot_grch37: false,
            backfill_transcripts: Some(acc_file),
            ..Default::default()
        };
        let err =
            prepare_references(&config).expect_err("partial-stale-cdot backfill must be rejected");
        let msg = format!("{err}");
        assert!(
            msg.contains("requires every referenced cdot file to exist"),
            "unexpected error: {err}"
        );
        assert!(
            msg.contains("cdot_grch37.json"),
            "error should name the missing GRCh37 cdot path: {err}"
        );
    }

    #[test]
    fn builds_for_genome_maps_selection() {
        assert_eq!(builds_for_genome("grch38"), vec!["GRCh38".to_string()]);
        assert_eq!(builds_for_genome("grch37"), vec!["GRCh37".to_string()]);
        assert_eq!(
            builds_for_genome("all"),
            vec!["GRCh38".to_string(), "GRCh37".to_string()]
        );
        // `none` provisions no genome build, so it has no builds to derive
        // placements for (the up-front guard rejects --derive-ng-placements
        // with this selection).
        assert!(builds_for_genome("none").is_empty());
    }

    #[test]
    fn skip_existing_still_records_assembly_report() {
        let tmp = tempfile::tempdir().unwrap();
        let genome_dir = tmp.path().join("genome");
        std::fs::create_dir_all(&genome_dir).unwrap();
        // Pre-place a genome FASTA + its .fai and the assembly report so the
        // skip_existing branch is taken (no network).
        let fasta = genome_dir.join("GRCh38.fna");
        std::fs::write(&fasta, ">NC_000001.11\nACGT\n").unwrap();
        std::fs::write(
            format!("{}.fai", fasta.display()),
            "NC_000001.11\t4\t13\t4\t5\n",
        )
        .unwrap();
        let report = genome_dir.join("GRCh38.assembly_report.txt");
        std::fs::write(&report, "# Assembly name:  GRCh38.p14\n").unwrap();

        let manifest = record_assembly_report(
            genome_dir.as_path(),
            "GRCh38.assembly_report.txt",
            urls::GRCH38_ASSEMBLY_REPORT,
            /* skip_existing */ true,
        )
        .unwrap();

        assert_eq!(manifest, Some(report));
    }

    #[test]
    fn build_and_write_ng_hosted_writes_artifact_and_returns_path() {
        let tmp = tempfile::tempdir().unwrap();
        let out = tmp.path().join("ng_hosted_transcripts.json");
        let recs = vec![(
            "NG_012337.1".to_string(),
            vec![("TIMM8B".to_string(), "NM_012459.2".to_string())],
        )];
        let p = write_ng_hosted_transcripts(&recs, &out).unwrap();
        assert_eq!(p, out);
        let loaded = crate::reference::ng_hosted_transcripts::NgHostedTranscripts::from_json(
            &std::fs::read_to_string(&out).unwrap(),
        )
        .unwrap();
        assert_eq!(
            loaded.hosted_unique("NG_012337.1", "TIMM8B"),
            Some("NM_012459.2")
        );
    }
}
