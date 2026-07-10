// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! ferro CLI
//!
//! Command-line interface for HGVS variant normalization and VCF annotation.

use clap::{Parser, Subcommand};
use ferro_hgvs::cli::{
    normalize_tsv_row, normalize_tsv_summary, output_error as cli_output_error,
    output_error_with_context as cli_output_error_with_context, parse_genome_build,
    parse_shuffle_direction, parse_vcf_line, process_input_line, OutputFormat,
    NORMALIZE_TSV_HEADER,
};
use ferro_hgvs::config::{apply_override, FerroConfig, OverrideSource};
use ferro_hgvs::error_handling::{
    get_code_info, list_all_codes, list_error_codes, list_warning_codes, ErrorConfig, ErrorMode,
};
use ferro_hgvs::reference::TranscriptDb;
use ferro_hgvs::vcf::{generate_info_header_lines, open_vcf, VcfAnnotator, VcfRecord};
use ferro_hgvs::{parse_hgvs, FerroError, NormalizeConfig, Normalizer, ReferenceProvider};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Parser)]
#[command(name = "ferro")]
#[command(author, version, about = "HGVS variant normalizer and VCF annotator")]
#[command(
    long_about = "Normalize HGVS variant descriptions and annotate VCF files.

Examples:
  ferro normalize 'NM_000088.3:c.459del'
  ferro normalize -i variants.txt
  echo 'NC_000001.11:g.12345A>G' | ferro normalize
  ferro parse 'NM_000088.3:c.459del'
  ferro annotate-vcf -i input.vcf -o output.vcf --gff transcripts.gff3"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Annotate a VCF file with HGVS notation
    AnnotateVcf {
        /// Input VCF file (use - for stdin)
        #[arg(short, long)]
        input: PathBuf,

        /// Output VCF file (use - for stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Transcript annotation file (GFF3 or GTF)
        #[arg(long)]
        gff: Option<PathBuf>,

        /// Reference FASTA file
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Genome build (GRCh37 or GRCh38)
        #[arg(long, default_value = "GRCh38")]
        build: String,

        /// Use VEP-style ANN field instead of separate INFO fields
        #[arg(long)]
        ann_format: bool,

        /// Include all transcripts (default: true)
        #[arg(long, default_value = "true")]
        all_transcripts: bool,

        /// Number of parallel workers for record annotation (default: 1).
        /// Only applies to file input; stdin is always processed serially.
        #[arg(short = 'j', long, default_value = "1")]
        workers: usize,
    },

    /// Liftover genomic coordinates between genome builds
    Liftover {
        /// HGVS variant or position (e.g., NC_000001.10:g.12345A>G or chr1:12345)
        position: Option<String>,

        /// Input file (one position per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Chain file for liftover (e.g., hg19ToHg38.over.chain.gz)
        #[arg(long, required = true)]
        chain: PathBuf,

        /// Reverse chain file (for bidirectional liftover)
        #[arg(long)]
        reverse_chain: Option<PathBuf>,

        /// Source genome build (GRCh37 or GRCh38)
        #[arg(long, default_value = "GRCh37")]
        from: String,

        /// Target genome build (GRCh37 or GRCh38)
        #[arg(long, default_value = "GRCh38")]
        to: String,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,
    },

    /// Generate HGVS description from reference and observed sequences
    Describe {
        /// Reference sequence
        #[arg(long, required_unless_present = "input")]
        reference: Option<String>,

        /// Observed sequence
        #[arg(long, required_unless_present = "input")]
        observed: Option<String>,

        /// Accession to use in HGVS output (e.g., NC_000001.11)
        #[arg(long)]
        accession: Option<String>,

        /// Input file with tab-separated reference and observed sequences
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,

        /// Detect duplications
        #[arg(long, default_value = "true")]
        detect_duplications: bool,

        /// Detect inversions
        #[arg(long, default_value = "true")]
        detect_inversions: bool,
    },

    /// Predict protein effect from variant
    Effect {
        /// Protein variant (e.g., p.Val600Glu) or coding variant (e.g., c.1799T>A)
        variant: Option<String>,

        /// Input file (one variant per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,
    },

    /// Backtranslate protein variant to possible DNA variants
    Backtranslate {
        /// Protein change (e.g., p.Val600Glu or V600E)
        variant: Option<String>,

        /// Reference amino acid (single letter, e.g., V)
        #[arg(long)]
        ref_aa: Option<String>,

        /// Alternate amino acid (single letter, e.g., E)
        #[arg(long)]
        alt_aa: Option<String>,

        /// Input file (one protein change per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,
    },

    /// Convert VCF to HGVS notation
    VcfToHgvs {
        /// Input VCF file (use - for stdin)
        #[arg(short, long)]
        input: PathBuf,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,

        /// Transcript annotation file (GFF3 or GTF)
        #[arg(long)]
        gff: Option<PathBuf>,

        /// Genome build (GRCh37 or GRCh38)
        #[arg(long, default_value = "GRCh38")]
        build: String,
    },

    /// Convert HGVS to VCF format
    HgvsToVcf {
        /// Input HGVS variant or file
        variant: Option<String>,

        /// Input file (one variant per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "vcf", value_parser = ["vcf", "text", "json"])]
        format: String,

        /// Reference data file (JSON)
        #[arg(long)]
        reference: Option<PathBuf>,
    },

    /// Normalize HGVS variants
    Normalize {
        /// HGVS variant to normalize
        variant: Option<String>,

        /// Input file (one variant per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Output format. `tsv` writes a header-bearing table (one row per
        /// input: line, input, normalized, changed, status, detail) plus a
        /// summary line on stderr, for answering "which of my variants
        /// changed?".
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json", "tsv"])]
        format: String,

        /// Shuffle direction (3prime or 5prime)
        #[arg(long, default_value = "3prime")]
        direction: String,

        /// Reference directory (with manifest.json from 'ferro prepare')
        #[arg(long)]
        reference: Option<PathBuf>,

        /// Error handling mode: strict, lenient, or silent
        #[arg(long, default_value = "strict", value_parser = ["strict", "lenient", "silent"])]
        error_mode: String,

        /// Warning codes to ignore (comma-separated, e.g., W1001,W2001)
        #[arg(long, value_delimiter = ',')]
        ignore: Vec<String>,

        /// Warning codes to always reject (comma-separated, e.g., W3003)
        #[arg(long, value_delimiter = ',')]
        reject: Vec<String>,

        /// Output timing information to JSON file
        #[arg(short = 't', long)]
        timing: Option<PathBuf>,

        /// Number of parallel workers (default: 1)
        #[arg(short = 'j', long, default_value = "1")]
        workers: usize,

        /// Hard-fail if the reference's content does not match its recorded
        /// identity (drifted in place). Default: warn and proceed (#1001).
        #[arg(long)]
        strict_reference: bool,
    },

    /// Project an HGVS variant onto a chosen output axis (g/c/n/p/r)
    Project {
        /// HGVS variant to project
        variant: Option<String>,

        /// Output axis: g (genomic), c (coding), n (non-coding), p (protein), r (RNA)
        #[arg(long, required = true, value_parser = ["g", "c", "n", "p", "r"])]
        axis: String,

        /// Transcript accession to project onto (required for a bare g. input)
        #[arg(long)]
        transcript: Option<String>,

        /// Input file (one variant per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,

        /// Reference directory (with manifest.json from 'ferro prepare')
        #[arg(long, required = true)]
        reference: PathBuf,

        /// Hard-fail if the reference's content does not match its recorded
        /// identity (drifted in place). Default: warn and proceed (#1001).
        #[arg(long)]
        strict_reference: bool,

        /// Error handling mode: strict, lenient, or silent
        #[arg(long, default_value = "strict", value_parser = ["strict", "lenient", "silent"])]
        error_mode: String,

        /// Warning codes to ignore (comma-separated, e.g., W1001,W2001)
        #[arg(long, value_delimiter = ',')]
        ignore: Vec<String>,

        /// Warning codes to always reject (comma-separated, e.g., W3003)
        #[arg(long, value_delimiter = ',')]
        reject: Vec<String>,
    },

    /// Arbitrate an HGVS parse/normalize/projection disagreement between
    /// ferro and another tool (default: Mutalyzer). Reports who is
    /// spec-compliant.
    Arbitrate {
        /// HGVS variant to arbitrate.
        variant: Option<String>,

        /// Prepared reference directory (required for normalize/projection).
        #[arg(long, required = true)]
        reference: PathBuf,

        /// The other tool's output to compare against (skip auto-fetch).
        #[arg(long)]
        other_output: Option<String>,

        /// Label for the other tool when --other-output is given.
        #[arg(long, default_value = "provided")]
        other_tool: String,

        /// Mutalyzer base URL for auto-fetch.
        #[arg(long, default_value = "https://mutalyzer.nl")]
        mutalyzer_url: String,

        /// Output format.
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,
    },

    /// File a GitHub bug report from a completed `ferro arbitrate --format
    /// json` bundle where ferro was judged wrong. Shows the full issue body
    /// and requires confirmation before opening a browser.
    BugReport {
        /// Path to a `ferro arbitrate --format json` bundle, or `-` to read
        /// it from stdin. Required in v1 (no standalone fallback).
        #[arg(long)]
        from_arbitration: Option<PathBuf>,

        /// Informational category hint (e.g. from the driving skill). The
        /// bundle's own `category` field is authoritative for the rendered
        /// report; a mismatch only produces a warning.
        #[arg(long)]
        category: Option<String>,

        /// Free-text notes from the reporter, included verbatim in the body.
        #[arg(long)]
        notes: Option<String>,

        /// Never open a browser — always print the prefilled URL instead.
        /// This is also the automatic behavior when `$BROWSER` is unset or
        /// stdin is not an interactive terminal.
        #[arg(long)]
        no_open: bool,

        /// Include ferro version and OS/arch details in the issue body.
        /// Omitted by default to avoid leaking machine details.
        #[arg(long)]
        include_environment: bool,
    },

    /// Parse HGVS variants (validation only)
    Parse {
        /// HGVS variant to parse
        variant: Option<String>,

        /// Input file (one variant per line)
        #[arg(short, long)]
        input: Option<PathBuf>,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,

        /// Error handling mode: strict, lenient, or silent
        #[arg(long, default_value = "strict", value_parser = ["strict", "lenient", "silent"])]
        error_mode: String,

        /// Warning codes to ignore (comma-separated, e.g., W1001,W2001)
        #[arg(long, value_delimiter = ',')]
        ignore: Vec<String>,

        /// Warning codes to always reject (comma-separated, e.g., W3003)
        #[arg(long, value_delimiter = ',')]
        reject: Vec<String>,

        /// Output timing information to JSON file
        #[arg(short = 't', long)]
        timing: Option<PathBuf>,

        /// Number of parallel workers (default: 1)
        #[arg(short = 'j', long, default_value = "1")]
        workers: usize,
    },

    /// Explain an error or warning code
    Explain {
        /// Error or warning code to explain (e.g., E1001, W2001)
        code: Option<String>,

        /// List all codes
        #[arg(long)]
        list: bool,

        /// Show only error codes (with --list)
        #[arg(long)]
        errors_only: bool,

        /// Show only warning codes (with --list)
        #[arg(long)]
        warnings_only: bool,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json", "markdown"])]
        format: String,
    },

    /// Convert GFF3/GTF to transcripts.json format
    ConvertGff {
        /// Input GFF3 or GTF file
        #[arg(short, long)]
        gff: PathBuf,

        /// Reference FASTA file (for extracting sequences)
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Genome build (GRCh37 or GRCh38)
        #[arg(long, default_value = "GRCh38")]
        build: String,

        /// Only include MANE transcripts
        #[arg(long)]
        mane_only: bool,

        /// Filter to specific transcript IDs (comma-separated)
        #[arg(long)]
        transcripts: Option<String>,

        /// Filter to specific gene symbols (comma-separated)
        #[arg(long)]
        genes: Option<String>,

        /// Promote loader warnings to log-level errors and fail if any Severity::Error diagnostic is recorded.
        #[arg(long, conflicts_with = "silent")]
        strict: bool,

        /// Suppress loader diagnostic warnings; still records them in the report.
        #[arg(long, conflicts_with = "strict")]
        silent: bool,

        /// Skip CDS-length and start-codon FASTA-aware validation even when --fasta is supplied.
        #[arg(long)]
        no_validate_fasta: bool,

        /// Also emit a top-level `genomic_sequences` map (the referenced contig bytes
        /// from --fasta, keyed by chromosome) so the resulting transcripts.json is
        /// genome-capable: `Normalizer(reference_json=...)` can then run the
        /// genome-aware normalization rules that need genomic sequence (#1026).
        /// Requires --fasta. Emits full contig sequences, so it can be large for a
        /// whole-genome FASTA — intended for small/synthetic references.
        #[arg(long, requires = "fasta")]
        emit_genomic_sequences: bool,

        /// Write the bounded sample of loader diagnostics as JSON to this path.
        #[arg(long, value_name = "PATH")]
        diagnostics_json: Option<PathBuf>,
    },

    /// Generate HGVS descriptions from components (Name Generator)
    Generate {
        /// Reference sequence accession (e.g., NC_000007.14, NM_000088.3)
        #[arg(long, required = true)]
        accession: String,

        /// Coordinate system: g (genomic), c (coding), n (non-coding), p (protein), r (RNA), m (mitochondrial)
        #[arg(long, required = true, value_parser = ["g", "c", "n", "p", "r", "m"])]
        coord: String,

        /// Variant type: sub, del, ins, dup, delins, inv, repeat
        #[arg(long, required = true, value_parser = ["sub", "del", "ins", "dup", "delins", "inv", "repeat"])]
        variant_type: String,

        /// Start position (1-based)
        #[arg(long, required = true)]
        pos: i64,

        /// End position for intervals (optional, for del/dup/inv/delins)
        #[arg(long)]
        end: Option<i64>,

        /// Reference base(s) for substitution or delins
        #[arg(long)]
        ref_base: Option<String>,

        /// Alternate base(s) for substitution, insertion, or delins
        #[arg(long)]
        alt_base: Option<String>,

        /// Intronic offset (e.g., +5, -10) for c. coordinates
        #[arg(long)]
        offset: Option<i64>,

        /// Repeat count for repeat variants
        #[arg(long)]
        repeat_count: Option<u32>,

        /// Reference amino acid (3-letter code for protein variants)
        #[arg(long)]
        ref_aa: Option<String>,

        /// Alternate amino acid (3-letter code for protein variants)
        #[arg(long)]
        alt_aa: Option<String>,

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
        format: String,

        /// Validate the generated HGVS by parsing it
        #[arg(long, default_value = "true")]
        validate: bool,
    },

    /// Extract HGVS patterns from VEP-annotated VCF files
    ExtractHgvs {
        /// Input VCF file (gzipped or uncompressed)
        input: PathBuf,

        /// Output file for extracted HGVS patterns
        output: PathBuf,

        /// VEP annotation field name
        #[arg(long, default_value = "CSQ")]
        field: String,

        /// Index of HGVSc in VEP annotation (0-based)
        #[arg(long, default_value = "29")]
        hgvsc_idx: usize,

        /// Index of HGVSp in VEP annotation (0-based)
        #[arg(long, default_value = "30")]
        hgvsp_idx: usize,

        /// Only extract patterns starting with this prefix
        #[arg(long)]
        search: Option<String>,
    },

    /// Prepare reference data for normalization
    Prepare {
        /// Output directory for reference data
        #[arg(short, long, default_value = "ferro-reference")]
        output_dir: PathBuf,

        /// Genome build to download: grch38, grch37, all, or none
        #[arg(long, default_value = "grch38", value_parser = ["grch38", "grch37", "all", "none"])]
        genome: String,

        /// Skip transcript FASTA sequences (~2GB)
        #[arg(long)]
        no_transcripts: bool,

        /// Also download companion protein FASTAs (~2GB) for the canonical
        /// translation check (issue #520); off by default
        #[arg(long)]
        proteins: bool,

        /// Skip cdot transcript metadata
        #[arg(long)]
        no_cdot: bool,

        /// Skip RefSeqGene sequences (~600MB)
        #[arg(long)]
        no_refseqgene: bool,

        /// Skip LRG sequences (~50MB)
        #[arg(long)]
        no_lrg: bool,

        /// Force re-download even if files exist
        #[arg(long)]
        force: bool,

        /// Pattern file for deriving supplemental accessions
        #[arg(long)]
        patterns: Option<PathBuf>,

        /// ClinVar file for deriving supplemental accessions
        #[arg(long)]
        clinvar: Option<PathBuf>,

        /// Newline-delimited accession list to validate against authoritative
        /// GenBank records, writing canonical_overrides.json (issue #520)
        #[arg(long)]
        validate_canonical: Option<PathBuf>,

        /// Also download Ensembl reference (cDNA sequences + cdot metadata) so
        /// variants on Ensembl transcripts/genes (ENST/ENSG/ENSP) resolve
        #[arg(long)]
        ensembl: bool,

        /// Newline-delimited file of exact NG_ versions (e.g. `NG_012337.3`) to
        /// derive version-independent genomic placements for at prepare time,
        /// writing `derived_refseqgene_placements.json` and wiring the manifest
        /// field. Requires cdot + genome in the same prepare run. Networked
        /// (NCBI EFetch per accession); per-accession failures warn and continue.
        #[arg(long)]
        derive_ng_placements: Option<PathBuf>,

        /// Newline-delimited file of exact `accession.version`s (e.g.
        /// `NM_002001.2`) to backfill at prepare time (#842): for each that is
        /// present in cdot but absent from the bulk RefSeq RNA FASTA, fetch its
        /// deposited sequence and append it to `backfill/backfill_transcripts.fna`,
        /// wiring the manifest field so the primary path serves it instead of
        /// synthesizing lossy bases from the genome. Requires cdot (downloaded
        /// in this run or already in the reference). Networked (NCBI EFetch per
        /// accession); per-accession failures warn and continue.
        ///
        /// Example: `ferro prepare --backfill-transcripts targets.txt`
        #[arg(long)]
        backfill_transcripts: Option<PathBuf>,

        /// Dry run - show what would be downloaded without downloading
        #[arg(long)]
        dry_run: bool,
    },

    /// Check reference data setup
    Check {
        /// Reference to check: a prepared reference directory (containing
        /// `manifest.json`), or a standalone `transcripts.json` file produced by
        /// `convert-gff`/`build-transcript`.
        #[arg(long, default_value = "ferro-reference")]
        reference: PathBuf,

        /// After checking, load the reference once to build/refresh the on-disk
        /// cdot bincode cache. Pay the one-time cache build here as a setup step
        /// so it doesn't land inside (and slow down the start of) a real run.
        #[arg(long)]
        build_cache: bool,

        /// Validate that every coding transcript's CDS start codon is a
        /// recognized translation initiator — `ATG`, or the alternative
        /// initiators `CTG`/`GTG`/`TTG` — in the served FASTA (issue #629).
        /// Catches a cdot-release / mRNA-FASTA revision mismatch where
        /// `cds_start` has drifted off the true start, which silently breaks
        /// reference-protein translation. Loads the full reference, so it is
        /// opt-in.
        #[arg(long)]
        validate_cds: bool,

        /// Accessions (one per line, `#` comments allowed; versioned or
        /// versionless) exempt from the `--validate-cds` start-codon check —
        /// transcripts with a legitimate non-`ATG` start codon.
        #[arg(long, value_name = "FILE")]
        cds_allowlist: Option<PathBuf>,

        /// Compute and write the content identity into the manifest without a
        /// full re-prepare (#1001). Refuses to overwrite a differing existing
        /// stamp unless `--force`.
        #[arg(long)]
        write_identity: bool,

        /// Allow `--write-identity` to overwrite an existing, differing stamp.
        /// Only meaningful with `--write-identity`, so clap requires it — alone
        /// it would be silently ignored.
        #[arg(long, requires = "write_identity")]
        force: bool,
    },

    /// Build a transcripts.json from a FASTA + CDS coordinates (single-exon).
    ///
    /// Useful for synthetic constructs (plasmids, reporter genes, custom amplicons)
    /// where the user already knows the CDS bounds and does not need a GFF3.
    BuildTranscript {
        /// Path to the FASTA file (indexed or plain; index will be built on the fly if absent).
        #[arg(long)]
        fasta: PathBuf,

        /// CDS start position (1-based inclusive, in transcript coordinates).
        ///
        /// For minus-strand transcripts the sequence is reverse-complemented before
        /// emission, so the position is relative to the reverse-complemented sequence.
        #[arg(long)]
        cds_start: u64,

        /// CDS end position (1-based inclusive, in transcript coordinates).
        ///
        /// For minus-strand transcripts the sequence is reverse-complemented before
        /// emission, so the position is relative to the reverse-complemented sequence.
        #[arg(long)]
        cds_end: u64,

        /// Output path for transcripts.json.
        #[arg(long, short = 'o')]
        output: PathBuf,

        /// Transcript ID (default: FASTA contig name).
        #[arg(long)]
        id: Option<String>,

        /// Strand: + or -.
        #[arg(long, default_value = "+")]
        strand: String,

        /// Contig name to use when the FASTA has multiple contigs.
        #[arg(long)]
        contig: Option<String>,

        /// Optional gene symbol to embed in the transcript record.
        #[arg(long)]
        gene: Option<String>,

        /// Genome build name embedded in the output (default: GRCh38).
        #[arg(long, default_value = "GRCh38")]
        genome_build: String,

        /// Also emit a top-level `genomic_sequences` map (the contig's forward bytes,
        /// keyed by chromosome) so the resulting transcripts.json is genome-capable:
        /// `Normalizer(reference_json=...)` can then run the genome-aware
        /// normalization rules that need genomic sequence (#1026).
        #[arg(long)]
        emit_genomic_sequences: bool,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::AnnotateVcf {
            input,
            output,
            gff,
            fasta: _,
            build,
            ann_format,
            all_transcripts,
            workers,
        } => run_annotate_vcf(
            &input,
            output.as_ref(),
            gff.as_ref(),
            &build,
            ann_format,
            all_transcripts,
            workers,
        ),
        Commands::Liftover {
            position,
            input,
            chain,
            reverse_chain,
            from,
            to,
            format,
        } => run_liftover(
            position.as_deref(),
            input.as_ref(),
            &chain,
            reverse_chain.as_ref(),
            &from,
            &to,
            &format,
        ),
        Commands::Describe {
            reference,
            observed,
            accession,
            input,
            format,
            detect_duplications,
            detect_inversions,
        } => run_describe(
            reference.as_deref(),
            observed.as_deref(),
            accession.as_deref(),
            input.as_ref(),
            &format,
            detect_duplications,
            detect_inversions,
        ),
        Commands::Effect {
            variant,
            input,
            format,
        } => run_effect(variant.as_deref(), input.as_ref(), &format),
        Commands::Backtranslate {
            variant,
            ref_aa,
            alt_aa,
            input,
            format,
        } => run_backtranslate(
            variant.as_deref(),
            ref_aa.as_deref(),
            alt_aa.as_deref(),
            input.as_ref(),
            &format,
        ),
        Commands::VcfToHgvs {
            input,
            format,
            gff,
            build,
        } => run_vcf_to_hgvs(&input, &format, gff.as_ref(), &build),
        Commands::HgvsToVcf {
            variant,
            input,
            format,
            reference,
        } => run_hgvs_to_vcf(
            variant.as_deref(),
            input.as_ref(),
            &format,
            reference.as_ref(),
        ),
        Commands::Normalize {
            variant,
            input,
            output,
            format,
            direction,
            reference,
            error_mode,
            ignore,
            reject,
            timing,
            workers,
            strict_reference,
        } => {
            let config = build_error_config(&error_mode, &ignore, &reject);
            run_normalize(
                variant.as_deref(),
                input.as_ref(),
                output.as_ref(),
                &format,
                &direction,
                reference.as_ref(),
                timing.as_ref(),
                workers,
                &config,
                strict_reference,
            )
        }
        Commands::Project {
            variant,
            axis,
            transcript,
            input,
            output,
            format,
            reference,
            strict_reference,
            error_mode,
            ignore,
            reject,
        } => {
            let config = build_error_config(&error_mode, &ignore, &reject);
            run_project(
                variant.as_deref(),
                &axis,
                transcript.as_deref(),
                input.as_ref(),
                output.as_ref(),
                &format,
                &reference,
                strict_reference,
                &config,
            )
        }
        Commands::Arbitrate {
            variant,
            reference,
            other_output,
            other_tool,
            mutalyzer_url,
            format,
        } => run_arbitrate(
            variant.as_deref(),
            &reference,
            other_output.as_deref(),
            &other_tool,
            &mutalyzer_url,
            &format,
        ),
        Commands::BugReport {
            from_arbitration,
            category,
            notes,
            no_open,
            include_environment,
        } => run_bug_report(
            from_arbitration.as_deref(),
            category.as_deref(),
            notes.as_deref(),
            no_open,
            include_environment,
        ),
        Commands::Parse {
            variant,
            input,
            output,
            format,
            error_mode,
            ignore,
            reject,
            timing,
            workers,
        } => {
            let config = build_error_config(&error_mode, &ignore, &reject);
            run_parse(
                variant.as_deref(),
                input.as_ref(),
                output.as_ref(),
                &format,
                timing.as_ref(),
                workers,
                &config,
            )
        }
        Commands::Explain {
            code,
            list,
            errors_only,
            warnings_only,
            format,
        } => run_explain(code.as_deref(), list, errors_only, warnings_only, &format),
        Commands::ConvertGff {
            gff,
            fasta,
            output,
            build,
            mane_only,
            transcripts,
            genes,
            strict,
            silent,
            no_validate_fasta,
            emit_genomic_sequences,
            diagnostics_json,
        } => run_convert_gff(
            &gff,
            fasta.as_ref(),
            output.as_ref(),
            &build,
            mane_only,
            transcripts.as_deref(),
            genes.as_deref(),
            strict,
            silent,
            no_validate_fasta,
            emit_genomic_sequences,
            diagnostics_json.as_ref(),
        ),
        Commands::Generate {
            accession,
            coord,
            variant_type,
            pos,
            end,
            ref_base,
            alt_base,
            offset,
            repeat_count,
            ref_aa,
            alt_aa,
            format,
            validate,
        } => run_generate(
            &accession,
            &coord,
            &variant_type,
            pos,
            end,
            ref_base.as_deref(),
            alt_base.as_deref(),
            offset,
            repeat_count,
            ref_aa.as_deref(),
            alt_aa.as_deref(),
            &format,
            validate,
        ),
        Commands::ExtractHgvs {
            input,
            output,
            field,
            hgvsc_idx,
            hgvsp_idx,
            search,
        } => run_extract_hgvs(
            &input,
            &output,
            &field,
            hgvsc_idx,
            hgvsp_idx,
            search.as_deref(),
        ),
        Commands::Prepare {
            output_dir,
            genome,
            no_transcripts,
            proteins,
            no_cdot,
            no_refseqgene,
            no_lrg,
            force,
            patterns,
            clinvar,
            validate_canonical,
            ensembl,
            derive_ng_placements,
            backfill_transcripts,
            dry_run,
        } => run_prepare(
            &output_dir,
            &genome,
            no_transcripts,
            proteins,
            no_cdot,
            no_refseqgene,
            no_lrg,
            force,
            patterns.as_deref(),
            clinvar.as_deref(),
            validate_canonical.as_deref(),
            ensembl,
            derive_ng_placements,
            backfill_transcripts,
            dry_run,
        ),
        Commands::Check {
            reference,
            build_cache,
            validate_cds,
            cds_allowlist,
            write_identity,
            force,
        } => run_check(
            &reference,
            build_cache,
            validate_cds,
            cds_allowlist.as_deref(),
            write_identity,
            force,
        ),
        Commands::BuildTranscript {
            fasta,
            cds_start,
            cds_end,
            output,
            id,
            strand,
            contig,
            gene,
            genome_build,
            emit_genomic_sequences,
        } => run_build_transcript(
            &fasta,
            cds_start,
            cds_end,
            &output,
            id.as_deref(),
            &strand,
            contig.as_deref(),
            gene.as_deref(),
            &genome_build,
            emit_genomic_sequences,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_annotate_vcf(
    input: &PathBuf,
    output: Option<&PathBuf>,
    gff: Option<&PathBuf>,
    build: &str,
    ann_format: bool,
    all_transcripts: bool,
    workers: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    // Parse genome build
    let genome_build = parse_genome_build(build);

    // Load transcript database
    let db = if let Some(gff_path) = gff {
        let cfg =
            ferro_hgvs::reference::annotation::LoaderConfig::new().with_genome_build(genome_build);
        let (db, report) = ferro_hgvs::reference::annotation::load_annotations(gff_path, &cfg)?;
        if !report.diagnostics_by_code.is_empty() {
            eprintln!("{}", report.summary_line());
        }
        db
    } else {
        // Create empty database
        TranscriptDb::with_build(genome_build)
    };

    // Create annotator
    let annotator = VcfAnnotator::new(&db)
        .use_ann_field(ann_format)
        .include_all_transcripts(all_transcripts);

    // Set up output. Wrap in a `BufWriter` so the per-record `writeln!`s don't
    // each cost a `write` syscall: annotating a 1M-record VCF previously spent
    // ~35% of its time in the kernel doing one unbuffered write per output line
    // (unlike `normalize`/`parse`, whose writers were already buffered).
    let mut writer: Box<dyn Write> = if let Some(out_path) = output {
        if out_path.to_string_lossy() == "-" {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            Box::new(BufWriter::new(std::fs::File::create(out_path)?))
        }
    } else {
        Box::new(BufWriter::new(io::stdout()))
    };

    // Read VCF and annotate
    if input.to_string_lossy() == "-" {
        // Read from stdin
        let stdin = io::stdin();
        let mut in_header = true;

        // Write header additions
        for line in generate_info_header_lines() {
            writeln!(writer, "{}", line)?;
        }

        for line in stdin.lock().lines() {
            let line = line?;
            if line.starts_with('#') {
                if in_header {
                    // Pass through header lines
                    writeln!(writer, "{}", line)?;
                }
            } else {
                in_header = false;
                // Parse and annotate in place (we own the parsed record).
                if let Ok(mut record) = parse_vcf_line(&line) {
                    annotator.annotate_into(&mut record)?;
                    writeln!(writer, "{}", record)?;
                } else {
                    writeln!(writer, "{}", line)?;
                }
            }
        }
    } else {
        // Read from file using noodles
        let reader = open_vcf(input)?;

        // Write new header with INFO fields
        writeln!(writer, "##fileformat=VCFv4.2")?;
        for line in generate_info_header_lines() {
            writeln!(writer, "{}", line)?;
        }
        writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

        // Process records. Annotation is independent per record, so with
        // `workers > 1` we annotate fixed-size chunks across a rayon pool and
        // write the results in input order — output is byte-identical to the
        // serial path. The first per-record annotation error is propagated
        // (matching the serial `?`).
        if workers > 1 {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()?;
            const ANNOTATE_CHUNK: usize = 4096;
            let mut records = reader.records();
            let mut chunk: Vec<VcfRecord> = Vec::with_capacity(ANNOTATE_CHUNK);
            loop {
                chunk.clear();
                for record in records.by_ref() {
                    chunk.push(record?);
                    if chunk.len() >= ANNOTATE_CHUNK {
                        break;
                    }
                }
                if chunk.is_empty() {
                    break;
                }
                let annotated: Vec<Result<String, FerroError>> = pool.install(|| {
                    chunk
                        .par_iter_mut()
                        .map(|r| annotator.annotate_into(r).map(|()| r.to_string()))
                        .collect()
                });
                for line in annotated {
                    writeln!(writer, "{}", line?)?;
                }
            }
        } else {
            for record in reader.records() {
                let mut record = record?;
                annotator.annotate_into(&mut record)?;
                writeln!(writer, "{}", record)?;
            }
        }
    }

    // Flush the BufWriter explicitly so a write error surfaces here rather than
    // being silently swallowed by the drop-time flush.
    writer.flush()?;
    Ok(())
}

fn run_vcf_to_hgvs(
    input: &PathBuf,
    format: &str,
    gff: Option<&PathBuf>,
    build: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Parse genome build
    let genome_build = parse_genome_build(build);

    // Load transcript database
    let db = if let Some(gff_path) = gff {
        let cfg =
            ferro_hgvs::reference::annotation::LoaderConfig::new().with_genome_build(genome_build);
        let (db, report) = ferro_hgvs::reference::annotation::load_annotations(gff_path, &cfg)?;
        if !report.diagnostics_by_code.is_empty() {
            eprintln!("{}", report.summary_line());
        }
        db
    } else {
        TranscriptDb::with_build(genome_build)
    };

    // Read VCF
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    if input.to_string_lossy() == "-" {
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let line = line?;
            if !line.starts_with('#') && !line.is_empty() {
                if let Ok(record) = parse_vcf_line(&line) {
                    output_vcf_hgvs(&record, &db, format, &mut handle)?;
                }
            }
        }
    } else {
        let reader = open_vcf(input)?;
        for record in reader.records() {
            let record = record?;
            output_vcf_hgvs(&record, &db, format, &mut handle)?;
        }
    }

    Ok(())
}

fn output_vcf_hgvs(
    record: &VcfRecord,
    db: &TranscriptDb,
    format: &str,
    handle: &mut impl Write,
) -> io::Result<()> {
    use ferro_hgvs::vcf::vcf_to_genomic_hgvs;

    // Convert to genomic HGVS
    for i in 0..record.alternate.len() {
        match vcf_to_genomic_hgvs(record, i) {
            Ok(variant) => {
                let hgvs = variant.to_string();
                match format {
                    "json" => {
                        writeln!(
                            handle,
                            r#"{{"chrom": "{}", "pos": {}, "ref": "{}", "alt": "{}", "hgvs": "{}"}}"#,
                            record.chrom,
                            record.pos,
                            record.reference,
                            record.alternate.get(i).unwrap_or(&String::new()),
                            hgvs
                        )?;
                    }
                    _ => {
                        writeln!(
                            handle,
                            "{}:{} {}/{} -> {}",
                            record.chrom,
                            record.pos,
                            record.reference,
                            record.alternate.get(i).unwrap_or(&String::new()),
                            hgvs
                        )?;
                    }
                }
            }
            Err(e) => {
                writeln!(handle, "ERROR: {}:{} - {}", record.chrom, record.pos, e)?;
            }
        }
    }

    // Also try transcript-specific conversion if DB is available
    if !db.is_empty() {
        use ferro_hgvs::vcf::MultiIsoformAnnotator;
        let annotator = MultiIsoformAnnotator::new(db);
        if let Ok(result) = annotator.annotate(record) {
            for ann in &result.annotations {
                let hgvs = ann.hgvs_string();
                match format {
                    "json" => {
                        writeln!(
                            handle,
                            r#"{{"transcript": "{}", "gene": "{}", "hgvs": "{}"}}"#,
                            ann.transcript_accession.as_deref().unwrap_or(""),
                            ann.gene_symbol.as_deref().unwrap_or(""),
                            hgvs
                        )?;
                    }
                    _ => {
                        writeln!(
                            handle,
                            "  {} ({}) -> {}",
                            ann.transcript_accession.as_deref().unwrap_or("unknown"),
                            ann.gene_symbol.as_deref().unwrap_or(""),
                            hgvs
                        )?;
                    }
                }
            }
        }
    }

    Ok(())
}

fn run_hgvs_to_vcf(
    variant: Option<&str>,
    input: Option<&PathBuf>,
    format: &str,
    _reference: Option<&PathBuf>,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::hgvs::variant::HgvsVariant;
    use ferro_hgvs::vcf::genomic_hgvs_to_vcf;

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // Print VCF header if outputting VCF format
    if format == "vcf" {
        writeln!(handle, "##fileformat=VCFv4.2")?;
        writeln!(handle, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    }

    let process_variant = |v: &str,
                           handle: &mut io::StdoutLock|
     -> Result<(), Box<dyn std::error::Error>> {
        let parsed = parse_hgvs(v)?;

        // Convert genomic variants directly
        if let HgvsVariant::Genome(ref gv) = parsed {
            let vcf = genomic_hgvs_to_vcf(gv)?;
            match format {
                "vcf" => {
                    writeln!(handle, "{}", vcf)?;
                }
                "json" => {
                    writeln!(
                        handle,
                        r#"{{"hgvs": "{}", "chrom": "{}", "pos": {}, "ref": "{}", "alt": "{}"}}"#,
                        v,
                        vcf.chrom,
                        vcf.pos,
                        vcf.reference,
                        vcf.alternate.first().unwrap_or(&String::new())
                    )?;
                }
                _ => {
                    writeln!(
                        handle,
                        "{} -> {}:{}:{}/{}",
                        v,
                        vcf.chrom,
                        vcf.pos,
                        vcf.reference,
                        vcf.alternate.first().unwrap_or(&String::new())
                    )?;
                }
            }
        } else {
            writeln!(handle, "ERROR: {} is not a genomic variant (g.)", v)?;
        }
        Ok(())
    };

    if let Some(v) = variant {
        process_variant(v, &mut handle)?;
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                process_variant(trimmed, &mut handle)?;
            }
        }
    } else {
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                process_variant(trimmed, &mut handle)?;
            }
        }
    }

    Ok(())
}

/// Outcome of processing one input line in the batch normalize/parse driver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LineStatus {
    /// Line was blank/comment and produced no output (not counted).
    Skipped,
    /// Variant processed successfully.
    Ok,
    /// Variant failed (error already formatted into `stderr`).
    Err,
}

/// One input line's fully-formatted output, captured into per-line buffers so
/// the batch driver can write results in input order regardless of whether they
/// were produced serially or in parallel.
struct LineResult {
    /// Bytes destined for the order-critical stdout/file stream.
    stdout: Vec<u8>,
    /// Bytes destined for stderr (warnings and/or error diagnostics).
    stderr: Vec<u8>,
    status: LineStatus,
}

impl LineResult {
    fn skipped() -> Self {
        LineResult {
            stdout: Vec::new(),
            stderr: Vec::new(),
            status: LineStatus::Skipped,
        }
    }
}

/// Write the `normalize --format tsv` header row, when that is the format.
///
/// Deliberately called only once the run's input is known to be readable: writing
/// it while opening the output stream made `--input missing.txt` leave behind a
/// header-only artifact byte-identical to the one a *successful* empty-input run
/// produces, with nothing on the output stream to tell the two apart.
fn write_tsv_header<W: Write>(writer: &mut W, tsv: bool) -> io::Result<()> {
    if tsv {
        writeln!(writer, "{}", NORMALIZE_TSV_HEADER)?;
    }
    Ok(())
}

/// Number of input lines processed per batch before their results are written.
/// Bounds memory (only one chunk is in flight) while giving rayon enough work
/// per fan-out to amortize scheduling overhead.
///
/// Re-exported from the library rather than defined twice: #975 built the
/// streaming batch API on this same engine, and two copies of the number would
/// have been two things to drift apart.
const BATCH_CHUNK_LINES: usize = ferro_hgvs::batch::BATCH_CHUNK_ITEMS;

/// Drive a batch of input lines through `item`, writing each line's output in
/// input order. With `workers <= 1` the work runs serially; with `workers > 1`
/// each chunk is mapped across a dedicated rayon pool (order preserved by
/// `par_iter().collect()`), then drained serially. Returns
/// `(total, success, error)` counts (skipped lines are not counted).
fn run_batch<I, F>(
    lines: I,
    writer: &mut dyn Write,
    workers: usize,
    item: &F,
) -> io::Result<(usize, usize, usize)>
where
    I: Iterator<Item = io::Result<String>>,
    F: Fn(usize, &str) -> LineResult + Sync,
{
    let pool = if workers > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .map_err(|e| io::Error::other(e.to_string()))?,
        )
    } else {
        None
    };
    let stderr = io::stderr();

    let mut counts = (0usize, 0usize, 0usize);
    let mut chunk: Vec<(usize, String)> = Vec::with_capacity(BATCH_CHUNK_LINES);
    for (line_num, line) in lines.enumerate() {
        chunk.push((line_num, line?));
        if chunk.len() >= BATCH_CHUNK_LINES {
            flush_chunk(&chunk, pool.as_ref(), item, writer, &stderr, &mut counts)?;
            chunk.clear();
        }
    }
    flush_chunk(&chunk, pool.as_ref(), item, writer, &stderr, &mut counts)?;
    Ok(counts)
}

/// Process one chunk (in parallel when `pool` is `Some`, else serially) and
/// write its results to `writer`/`stderr` in input order, updating
/// `(total, success, error)`.
fn flush_chunk<F>(
    chunk: &[(usize, String)],
    pool: Option<&rayon::ThreadPool>,
    item: &F,
    writer: &mut dyn Write,
    stderr: &io::Stderr,
    counts: &mut (usize, usize, usize),
) -> io::Result<()>
where
    F: Fn(usize, &str) -> LineResult + Sync,
{
    if chunk.is_empty() {
        return Ok(());
    }
    let results: Vec<LineResult> = match pool {
        Some(p) => p.install(|| chunk.par_iter().map(|(ln, l)| item(*ln, l)).collect()),
        None => chunk.iter().map(|(ln, l)| item(*ln, l)).collect(),
    };
    let mut err_handle = stderr.lock();
    for r in &results {
        match r.status {
            LineStatus::Skipped => continue,
            LineStatus::Ok => {
                counts.0 += 1;
                counts.1 += 1;
            }
            LineStatus::Err => {
                counts.0 += 1;
                counts.2 += 1;
            }
        }
        if !r.stdout.is_empty() {
            writer.write_all(&r.stdout)?;
        }
        if !r.stderr.is_empty() {
            err_handle.write_all(&r.stderr)?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn run_normalize(
    variant: Option<&str>,
    input: Option<&PathBuf>,
    output: Option<&PathBuf>,
    format: &str,
    direction: &str,
    reference: Option<&PathBuf>,
    timing: Option<&PathBuf>,
    workers: usize,
    error_config: &ErrorConfig,
    strict_reference: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::commands::create_reference_provider;
    use std::time::Instant;

    let preprocessor = error_config.preprocessor();
    // #1181: the error configuration must reach the *normalizer*, not just the
    // preprocessor, or every `--error-mode` (and every `--ignore`/`--reject`
    // override, which `build_error_config` folds into the same `ErrorConfig`)
    // is silently discarded and all modes behave as lenient.
    // `for_entry_point` takes the error configuration as a required argument,
    // so the omission cannot recur at this call site (#1197). It is not a
    // whole-crate guarantee — see the constructor's docs and the entry-point
    // scan in `tests/it/issue_1197_required_error_config.rs`.
    let config =
        NormalizeConfig::for_entry_point(parse_shuffle_direction(direction), error_config.clone());

    // Create reference provider from directory
    let provider = create_reference_provider(reference.map(|p| p.as_path()), strict_reference)?;
    let normalizer = Normalizer::with_config(provider, config);

    // Print capability summary
    print_normalize_capabilities_dir(reference);

    // Create output writer - either file or stdout
    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(io::stdout()),
    };
    let mut error_count = 0usize;
    let mut success_count = 0usize;
    let mut total_count = 0usize;
    let start = Instant::now();

    // TSV is a *reporting* format rather than a per-variant rendering: every
    // input becomes exactly one row (failures included), so it needs its own
    // per-variant path and a header written before any row. Keeping it entirely
    // out of `process` below is what guarantees `text` and `json` output is
    // byte-for-byte unchanged.
    let tsv = format == "tsv";
    // Counted separately from success/error because only TSV reports it, and the
    // batch driver's per-line closure must be `Sync` (it may run on a rayon
    // pool), so the counter cannot be a captured `&mut usize`.
    let changed_count = AtomicUsize::new(0);

    // Render one input variant as a TSV row, classifying a failure by the phase
    // it happened in. `line` is the 1-based input line number, or `None` for a
    // single positional variant, which has no line context. Warnings go to `err`
    // in the same shape as the text format *and* into the row's `detail`, so a
    // lenient-mode correction is both greppable where it always was and
    // correlatable to the row it explains.
    let tsv_row = |line: Option<usize>, v: &str, err: &mut dyn Write| -> (String, LineStatus) {
        use ferro_hgvs::cli::{NormalizeTsvFailure, NormalizeTsvOutcome};

        let failed = |phase: NormalizeTsvFailure, detail: &str| -> (String, LineStatus) {
            (
                normalize_tsv_row(line, v, NormalizeTsvOutcome::Failed(phase, detail)).row,
                LineStatus::Err,
            )
        };

        let mut preprocess_result = preprocessor.preprocess(v);
        if let Some(rejection) = preprocess_result.take_rejection_error() {
            return failed(NormalizeTsvFailure::Preprocess, &rejection.to_string());
        }
        let corrections: Vec<(&str, &str)> = preprocess_result
            .warnings
            .iter()
            .map(|w| (w.error_type.code(), w.message.as_str()))
            .collect();
        for warning in &preprocess_result.warnings {
            let _ = writeln!(
                err,
                "warning[{}]: {}",
                warning.error_type.code(),
                warning.message
            );
        }
        let parsed = match parse_hgvs(&preprocess_result.preprocessed) {
            Ok(parsed) => parsed,
            Err(e) => return failed(NormalizeTsvFailure::Parse, &e.to_string()),
        };
        match normalizer.normalize(&parsed) {
            Ok(normalized) => {
                let normalized_str = normalized.to_string();
                let row = normalize_tsv_row(
                    line,
                    v,
                    NormalizeTsvOutcome::Normalized {
                        normalized: normalized_str.as_str(),
                        corrections: &corrections,
                    },
                );
                // The verdict comes back from the renderer rather than being
                // recomputed here, so the summary's `changed` count and the
                // column cannot drift apart.
                if row.changed {
                    changed_count.fetch_add(1, Ordering::Relaxed);
                }
                (row.row, LineStatus::Ok)
            }
            Err(e) => failed(NormalizeTsvFailure::Normalize, &e.to_string()),
        }
    };

    // Per-line work, factored so it can run serially or in parallel: `out`
    // receives the success output (the order-critical stdout/file stream) and
    // `err` receives text-mode warnings. Both are written into per-line buffers
    // by the batch driver, so parallel output is byte-identical to serial.
    let process = |v: &str, out: &mut dyn Write, err: &mut dyn Write| -> Result<(), FerroError> {
        // Preprocess input according to error mode
        let mut preprocess_result = preprocessor.preprocess(v);
        if let Some(rejection) = preprocess_result.take_rejection_error() {
            // Preprocessing rejected the input: surface the phase's own
            // diagnosis (message, span, and structured code), not a summary
            // built from the — deliberately empty — warning list (#1162).
            return Err(rejection);
        }

        let parsed = parse_hgvs(&preprocess_result.preprocessed)?;
        let normalized = normalizer.normalize(&parsed)?;
        match format {
            "json" => {
                let corrections: Vec<String> = preprocess_result
                    .warnings
                    .iter()
                    .map(|w| {
                        format!(
                            r#"{{"code":"{}","message":"{}"}}"#,
                            w.error_type.code(),
                            w.message.replace('"', "\\\"")
                        )
                    })
                    .collect();
                writeln!(
                    out,
                    r#"{{"input": "{}", "output": "{}", "status": "ok", "corrections": [{}]}}"#,
                    v,
                    normalized,
                    corrections.join(",")
                )
                .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }
            _ => {
                // Emit warnings to the stderr buffer if any
                for warning in &preprocess_result.warnings {
                    writeln!(
                        err,
                        "warning[{}]: {}",
                        warning.error_type.code(),
                        warning.message
                    )
                    .map_err(|e| FerroError::Io { msg: e.to_string() })?;
                }
                // Format the normalized variant once and reuse it for both the
                // unchanged-vs-changed comparison and the written output. The
                // previous code called `normalized.to_string()` for the compare
                // and then formatted `normalized` again in `writeln!`, paying
                // the `Display` cost twice per line (~8% of a batch run).
                let normalized_str = normalized.to_string();
                if v == normalized_str {
                    writeln!(out, "{}", normalized_str)
                        .map_err(|e| FerroError::Io { msg: e.to_string() })?;
                } else {
                    writeln!(out, "{} -> {}", v, normalized_str)
                        .map_err(|e| FerroError::Io { msg: e.to_string() })?;
                }
            }
        }
        Ok(())
    };

    if let Some(v) = variant {
        // Single-variant path: unchanged (no line context, errors via
        // `output_error`). Warnings go straight to the live stderr.
        total_count += 1;
        let mut stderr = io::stderr();
        if tsv {
            // A failure is a row here too, so `--format tsv` on one positional
            // variant reports the same way as a one-line `--input` file — except
            // for the `line` column, which is empty because there is no file to
            // number lines in.
            write_tsv_header(&mut writer, tsv)?;
            let (row, status) = tsv_row(None, v, &mut stderr);
            writeln!(writer, "{}", row)?;
            if status == LineStatus::Ok {
                success_count += 1;
            } else {
                error_count += 1;
            }
        } else if let Err(e) = process(v, &mut writer, &mut stderr) {
            output_error(v, &e, format)?;
            error_count += 1;
        } else {
            success_count += 1;
        }
    } else {
        // Batch path (file or stdin): each input line is turned into a
        // `LineResult` (success output + warnings/error pre-formatted into
        // per-line buffers) by `item`, then `run_batch` writes them in input
        // order — serially or, when `workers > 1`, across a rayon pool. Because
        // every line is formatted into its own buffers with the exact same
        // code, parallel output is byte-identical to serial.
        let format_owned: OutputFormat = format.parse().unwrap_or_default();
        let item = |line_num: usize, line: &str| -> LineResult {
            let is_first = line_num == 0;
            let Some(variant_str) = process_input_line(line, is_first) else {
                return LineResult::skipped();
            };
            if tsv {
                // One row per line, failures included: the row *is* the error
                // report, so nothing goes to the per-line stderr buffer except
                // preprocessing warnings.
                let mut err = Vec::new();
                let (row, status) = tsv_row(Some(line_num + 1), variant_str, &mut err);
                let mut out = Vec::new();
                let _ = writeln!(out, "{}", row);
                return LineResult {
                    stdout: out,
                    stderr: err,
                    status,
                };
            }
            let mut out = Vec::new();
            let mut err = Vec::new();
            let status = match process(variant_str, &mut out, &mut err) {
                Ok(()) => LineStatus::Ok,
                Err(e) => {
                    // Mirror `output_error_with_line`: format the error (with
                    // 1-based line number) into the per-line stderr buffer.
                    let _ = cli_output_error_with_context(
                        &mut err,
                        variant_str,
                        &e,
                        format_owned,
                        Some(line_num + 1),
                    );
                    LineStatus::Err
                }
            };
            LineResult {
                stdout: out,
                stderr: err,
                status,
            }
        };

        let (t, s, er) = if let Some(input_path) = input {
            // Open before the header is written: a missing input file must not
            // leave a header-only artifact behind (see `write_tsv_header`).
            let file = std::fs::File::open(input_path)?;
            write_tsv_header(&mut writer, tsv)?;
            run_batch(
                io::BufReader::new(file).lines(),
                &mut writer,
                workers,
                &item,
            )?
        } else {
            let stdin = io::stdin();
            let reader = stdin.lock();
            write_tsv_header(&mut writer, tsv)?;
            run_batch(reader.lines(), &mut writer, workers, &item)?
        };
        total_count += t;
        success_count += s;
        error_count += er;
    }

    // The reference provider (cdot maps + FASTA index, ~1.4 GB) is owned by
    // `normalizer`. Running its destructor here only frees memory the OS
    // reclaims at process exit anyway — on a batch run ~10% of wall time is
    // spent walking the millions of map entries to drop them. Both `process` and
    // `tsv_row` borrowed `normalizer`, but the last use of either is above, so
    // those borrows have already ended; leak `normalizer` so the process can exit
    // without the teardown. `writer` is independent (output is passed in, not
    // captured) and still flushes on its own drop, so this does not affect output.
    std::mem::forget(normalizer);

    // The run summary answers "how many of my variants changed?" without the
    // caller having to re-read the table. It goes to stderr so the output stream
    // stays a table whose every line is a row.
    if tsv {
        // Flush before claiming counts: with `--output` the rows are sitting in a
        // `BufWriter` whose `Drop` would swallow a write error, so flushing here
        // is what turns a failed write into a reported failure rather than a
        // summary that asserts rows were emitted. (Bare stdout is a `LineWriter`,
        // already flushed per row, so there this is a no-op.)
        writer.flush()?;
        eprintln!(
            "{}",
            normalize_tsv_summary(
                total_count,
                changed_count.load(Ordering::Relaxed),
                error_count
            )
        );
    }

    let elapsed = start.elapsed();

    // Write timing info if requested
    if let Some(timing_path) = timing {
        let timing_info =
            ferro_hgvs::commands::TimingInfo::new(total_count, success_count, elapsed);
        ferro_hgvs::commands::write_timing(&timing_info, timing_path)?;
    }

    if error_count > 0 {
        Err(format!("{} variant(s) failed to normalize", error_count).into())
    } else {
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn run_project(
    variant: Option<&str>,
    axis: &str,
    transcript: Option<&str>,
    input: Option<&PathBuf>,
    output: Option<&PathBuf>,
    format: &str,
    reference: &Path,
    strict_reference: bool,
    error_config: &ErrorConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::cli::format::{output_project_error, output_projection};
    use ferro_hgvs::cli::project::AxisOutcome;
    use ferro_hgvs::cli::{project_axis, Axis, OutputFormat};
    use ferro_hgvs::data::cdot::CdotMapper;
    use ferro_hgvs::data::projection::Projector;
    use ferro_hgvs::project::VariantProjector;
    use ferro_hgvs::reference::multi_fasta::{LoadOptions, MultiFastaProvider};
    use std::sync::Arc;

    let axis = Axis::parse(axis).expect("clap value_parser guarantees a valid axis");
    // OutputFormat: FromStr (Infallible) — use the inherent str::parse so no
    // `use std::str::FromStr` is needed; Default = Text covers the never-taken err arm.
    let out_format: OutputFormat = format.parse().unwrap_or_default();

    // Build the projector. The Box<dyn ReferenceProvider> from
    // create_reference_provider is not Clone, and VariantProjector requires
    // P: Clone. MultiFastaProvider is itself NOT Clone, but `Arc<T>` has a
    // blanket `impl ReferenceProvider` (reference/provider.rs) and is Clone, so
    // build the concrete provider, extract the owned cdot, then project through
    // an Arc.
    let manifest_path = reference.join("manifest.json");
    let provider = MultiFastaProvider::from_manifest_with_options(
        &manifest_path,
        LoadOptions {
            strict_identity: strict_reference,
            ..Default::default()
        },
    )?;
    // Extract owned cdot BEFORE moving the provider into the Arc.
    let cdot = provider
        .cdot_mapper()
        .cloned()
        .unwrap_or_else(CdotMapper::new);
    let provider = Arc::new(provider);
    // #1182: `project` normalizes its input before projecting, so the error
    // configuration has to reach that normalizer or every normalizer-level
    // diagnostic (W4004 among them) is unreachable — the same defect #1181
    // fixed on `normalize`. This is what makes `--error-mode strict` reach the
    // normalizer here.
    //
    // It does *not* make the whole `ErrorConfig` effective on `project`:
    // unlike `run_normalize` and `run_parse`, this path never builds
    // `error_config.preprocessor()`, so the preprocessing half of `--ignore` /
    // `--reject` is still discarded. That is a pre-existing gap of the same
    // family as #1181, out of scope here; noted so the line above is not read
    // as a stronger claim than it makes.
    //
    // `project` has no direction flag; naming the 3' default explicitly is the
    // price of a constructor that cannot silently default the error
    // configuration (#1197), and it documents the choice at the call site.
    let projector = VariantProjector::new(Projector::new(cdot), provider).with_normalize_config(
        ferro_hgvs::normalize::NormalizeConfig::for_entry_point(
            ferro_hgvs::normalize::ShuffleDirection::ThreePrime,
            error_config.clone(),
        ),
    );

    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(io::stdout()),
    };
    let mut error_count = 0usize;

    let process = |v: &str, writer: &mut dyn Write, line: Option<usize>| -> bool {
        let parsed = match parse_hgvs(v) {
            Ok(p) => p,
            Err(e) => {
                let _ = cli_output_error_with_context(&mut io::stderr(), v, &e, out_format, line);
                return false;
            }
        };
        match project_axis(&projector, &parsed, axis, transcript) {
            Ok(outcome) => {
                // #1182: surface the normalization warnings the projection
                // carried. Text mode sends them to stderr so the result stream
                // stays one line per variant and pipelines are unaffected —
                // matching how `run_normalize` reports its warnings. JSON mode
                // carries them inside the object instead (see
                // `output_projection`), so emitting them here too would
                // duplicate them.
                if out_format != OutputFormat::Json {
                    let warnings = match &outcome {
                        AxisOutcome::Rendered { warnings, .. } => warnings,
                        AxisOutcome::Unavailable { warnings, .. } => warnings,
                    };
                    for w in warnings {
                        let _ = writeln!(io::stderr(), "warning[{}]: {}", w.code, w.message);
                    }
                }
                // A failed write to the result stream (e.g. a full or unwritable
                // --output file) must not be reported as success; count it so the
                // exit code reflects the truncated output. (A broken pipe from a
                // downstream `| head` lands here too and exits nonzero, which is
                // acceptable — the output was genuinely truncated.)
                if let Err(e) = output_projection(writer, v, axis, &outcome, out_format, line) {
                    let _ = writeln!(
                        io::stderr(),
                        "ERROR: failed writing output for {}: {}",
                        v,
                        e
                    );
                    return false;
                }
                true
            }
            Err(e) => {
                // Format-aware so `--format json` stays parseable (parity with the
                // parse-error path above, which uses cli_output_error_with_context).
                let _ = output_project_error(&mut io::stderr(), v, &e.0, out_format, line);
                false
            }
        }
    };

    if let Some(v) = variant {
        if !process(v, &mut writer, None) {
            error_count += 1;
        }
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            if let Some(vs) = process_input_line(&line, line_num == 0) {
                if !process(vs, &mut writer, Some(line_num + 1)) {
                    error_count += 1;
                }
            }
        }
    } else {
        let stdin = io::stdin();
        for (line_num, line) in stdin.lock().lines().enumerate() {
            let line = line?;
            if let Some(vs) = process_input_line(&line, line_num == 0) {
                if !process(vs, &mut writer, Some(line_num + 1)) {
                    error_count += 1;
                }
            }
        }
    }

    // Flush the BufWriter explicitly so a write error surfaces here rather than
    // being silently swallowed by the drop-time flush (parity with run_parse).
    writer.flush()?;

    if error_count > 0 {
        Err(format!("{} variant(s) failed to project", error_count).into())
    } else {
        Ok(())
    }
}

/// Provider type shared by the arbitrate command's normalizer and projector:
/// `Arc<MultiFastaProvider>` so both can own a cheap clone of the same
/// underlying reference (mirrors `run_project`'s provider idiom).
type ArbitrateProvider = std::sync::Arc<ferro_hgvs::reference::multi_fasta::MultiFastaProvider>;

/// Thin wrapper around a `VariantProjector` so `ferro arbitrate` can implement
/// `ferro_hgvs::arbitrate::FrameResolver` for it. The orphan rule blocks
/// implementing that trait directly for `VariantProjector<ArbitrateProvider>`
/// here: both the trait and the type are defined in the `ferro_hgvs` library
/// crate, and `src/bin/ferro.rs` is a separate (binary) crate — neither is
/// local to it. Wrapping the projector in a local newtype makes the impl
/// legal.
struct CliFrameResolver(ferro_hgvs::project::VariantProjector<ArbitrateProvider>);

impl ferro_hgvs::arbitrate::FrameResolver for CliFrameResolver {
    fn resolve(&self, v: &ferro_hgvs::HgvsVariant) -> Result<ferro_hgvs::HgvsVariant, FerroError> {
        ferro_hgvs::arbitrate::oracle::resolve_frame(v, &self.0)
    }
}

/// Obtain the "other tool's" result for `variant`: either the caller-supplied
/// `--other-output` verbatim, or an auto-fetch against a live Mutalyzer
/// instance. See the `MutalyzerResult` -> `OtherResult` mapping in the design
/// spec: a clean normalization maps straight through; an `EINTRONIC` genomic
/// redirect picks the option whose assembly matches the reference build (or
/// the first option when the build can't be determined), since that is the
/// option arbitration should actually compare against; a redirect with no
/// options, or a parse/network failure, reports `status: "unavailable"` with
/// no output so the caller falls back to `--other-output`.
fn fetch_other_result(
    variant: &str,
    other_output: Option<&str>,
    other_tool: &str,
    mutalyzer_url: &str,
    parsed: &ferro_hgvs::HgvsVariant,
    provider: &ArbitrateProvider,
) -> Result<ferro_hgvs::arbitrate::OtherResult, Box<dyn std::error::Error>> {
    use ferro_hgvs::arbitrate::OtherResult;
    use ferro_hgvs::mutalyzer::{MutalyzerClient, MutalyzerStatus};

    if let Some(output) = other_output {
        return Ok(OtherResult {
            tool: other_tool.to_string(),
            status: "ok".to_string(),
            output: Some(output.to_string()),
        });
    }

    let client = MutalyzerClient::new(mutalyzer_url)?;
    let result = client.normalize(variant)?;
    Ok(match result.status {
        MutalyzerStatus::Ok => OtherResult {
            tool: "mutalyzer".to_string(),
            status: "ok".to_string(),
            output: result.normalized,
        },
        MutalyzerStatus::GenomicRedirect => {
            // Prefer the option whose assembly matches the reference build the
            // input's own accession resolves to; fall back to the first option
            // when the build can't be determined (e.g. an accession the
            // reference has no alias data for).
            let build = parsed
                .accession()
                .and_then(|acc| provider.infer_genome_build(acc));
            let chosen = build
                .and_then(|b| {
                    result
                        .genomic_options
                        .iter()
                        .find(|o| o.assembly_id.eq_ignore_ascii_case(b))
                })
                .or_else(|| result.genomic_options.first());
            match chosen {
                Some(opt) => OtherResult {
                    tool: "mutalyzer".to_string(),
                    status: "ok".to_string(),
                    output: Some(opt.description.clone()),
                },
                None => OtherResult {
                    tool: "mutalyzer".to_string(),
                    status: "unavailable".to_string(),
                    output: None,
                },
            }
        }
        MutalyzerStatus::Unavailable | MutalyzerStatus::ParseError => OtherResult {
            tool: "mutalyzer".to_string(),
            status: "unavailable".to_string(),
            output: None,
        },
    })
}

/// Split the leading accession off an SPDI string (`accession:position:del:ins`).
fn spdi_accession(spdi: &str) -> Option<&str> {
    spdi.split(':').next()
}

/// The unversioned root of an accession (`NM_003002.4` -> `NM_003002`).
fn accession_root(accession: &str) -> &str {
    accession.split('.').next().unwrap_or(accession)
}

/// Render the verdict line: the first line of `--format text` output, and the
/// only part that branches on `(verdict, compliance, category)` rather than
/// `verdict` alone. Compliance is meaningful even when `verdict` is
/// `Equivalent` (see module-level design note in the task brief): a dup-vs-ins
/// disagreement where the other tool used the spec-forbidden spelling still
/// reduces to the same edit (`Equivalent`), but ferro's spelling is the one
/// the spec actually mandates and the other tool's is not (or vice versa).
fn render_verdict_line(a: &ferro_hgvs::arbitrate::Arbitration) -> String {
    use ferro_hgvs::arbitrate::{ArbitrationCategory, Compliance, Verdict};

    // Read the tool label from the arbitration result itself (`a.other.tool`)
    // rather than the raw `--other-tool` CLI flag: on an auto-fetch path the
    // flag still holds its "provided" default (it only applies to
    // `--other-output`), while `a.other.tool` is always the tool that
    // actually produced `a.other.output` (e.g. "mutalyzer").
    let other_tool = a.other.tool.as_str();

    match (a.verdict, a.compliance, a.category) {
        (Verdict::Equivalent, Compliance::NotApplicable, _) => {
            "VERDICT: Equivalent — same variant, both forms are valid spellings.".to_string()
        }
        (Verdict::Equivalent, Compliance::Ferro, ArbitrationCategory::FerroCorrect) => format!(
            "VERDICT: Equivalent — same variant; ferro's form is spec-compliant, {other_tool}'s is not."
        ),
        (Verdict::Equivalent, Compliance::Other, ArbitrationCategory::MutalyzerCorrect) => format!(
            "VERDICT: Equivalent — same variant, but ferro's spelling is NOT spec-compliant — \
             the spec mandates {other_tool}'s form (ferro notation bug)."
        ),
        (Verdict::Equivalent, ..) => {
            "VERDICT: Equivalent — same variant.".to_string()
        }
        (Verdict::Different, ..) => {
            "VERDICT: Different — these are different edits to the reference; the governing \
             rule below needs interpretation to say which (if either) is correct."
                .to_string()
        }
        (Verdict::BasisMismatch, ..) => {
            let ferro_acc = a.ferro_spdi.as_deref().and_then(spdi_accession);
            let other_acc = a.other_spdi.as_deref().and_then(spdi_accession);
            match (ferro_acc, other_acc) {
                (Some(fa), Some(oa))
                    if fa != oa && accession_root(fa) == accession_root(oa) =>
                {
                    format!(
                        "VERDICT: BasisMismatch — {other_tool} normalized on {oa}, your input \
                         is {fa} — align versions."
                    )
                }
                (Some(fa), Some(oa)) => format!(
                    "VERDICT: BasisMismatch — ferro reduced to {fa}, {other_tool} reduced to \
                     {oa}: no shared reference basis to compare on."
                ),
                _ => format!(
                    "VERDICT: BasisMismatch — ferro and {other_tool} normalized on different \
                     reference bases."
                ),
            }
        }
        (Verdict::OtherUnparseable, ..) => format!(
            "VERDICT: OtherUnparseable — {other_tool}'s output could not be parsed as HGVS."
        ),
        (Verdict::FerroParseError, ..) => {
            "VERDICT: FerroParseError — ferro's own output could not be parsed as HGVS."
                .to_string()
        }
        (Verdict::Inconclusive, ..) => {
            let reason = a.reason.as_deref().unwrap_or("unknown error");
            let mut line = format!("INCONCLUSIVE: could not adjudicate — {reason}");
            if reason.contains("NG/NC parent") {
                line.push_str(
                    "\n  (a bare NM_ transcript variant may need an explicit NC_(NM_) parent form)",
                );
            }
            line
        }
    }
}

/// Render the full `--format text` body: verdict line, then the governing
/// spec excerpt (if any), then the supporting facts.
fn render_arbitration_text(a: &ferro_hgvs::arbitrate::Arbitration) -> String {
    let mut out = render_verdict_line(a);
    out.push('\n');

    if !a.spec_citations.is_empty() {
        out.push_str("\nGoverning spec passage:\n");
        for c in &a.spec_citations {
            out.push_str(&format!(
                "  [{} — {}] {}\n    \"{}\"\n",
                c.file,
                c.heading,
                c.spec_version,
                c.excerpt.trim()
            ));
        }
    }

    out.push_str("\nFacts:\n");
    out.push_str(&format!("  input:            {}\n", a.input));
    if let Some(ferro_output) = &a.ferro_output {
        out.push_str(&format!("  ferro:            {ferro_output}\n"));
    }
    out.push_str(&format!(
        "  {} [{}]: {}\n",
        a.other.tool,
        a.other.status,
        a.other.output.as_deref().unwrap_or("<none>")
    ));
    if a.other.output.is_none() {
        out.push_str(&format!(
            "  (no {} output to compare — re-run with \
             `--other-output '<paste the tool's normalized description>'`)\n",
            a.other.tool
        ));
    }
    out.push_str(&format!("  category:         {}\n", a.category));
    out.push_str(&format!("  compliance:       {}\n", a.compliance));
    if let Some(spdi) = &a.ferro_spdi {
        out.push_str(&format!("  ferro spdi:       {spdi}\n"));
    }
    if let Some(spdi) = &a.other_spdi {
        out.push_str(&format!("  other spdi:       {spdi}\n"));
    }
    out
}

/// `ferro arbitrate`: adjudicate a parse/normalize/projection disagreement
/// between ferro and another tool (Mutalyzer by default) for a single HGVS
/// input. Builds the same provider+projector as `run_project`, computes
/// ferro's own normalized output, obtains the other tool's output
/// (`--other-output` or a live Mutalyzer fetch), and calls
/// `ferro_hgvs::arbitrate::arbitrate` to produce the verdict.
fn run_arbitrate(
    variant: Option<&str>,
    reference: &Path,
    other_output: Option<&str>,
    other_tool: &str,
    mutalyzer_url: &str,
    format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::arbitrate::arbitrate;
    use ferro_hgvs::data::cdot::CdotMapper;
    use ferro_hgvs::data::projection::Projector;
    use ferro_hgvs::project::VariantProjector;
    use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
    use std::sync::Arc;

    let variant = match variant {
        Some(v) => v.to_string(),
        None => {
            let mut line = String::new();
            io::stdin().read_line(&mut line)?;
            let line = line.trim().to_string();
            if line.is_empty() {
                return Err(
                    "no variant provided: pass it as an argument or pipe one line on stdin".into(),
                );
            }
            line
        }
    };

    // Build the provider + projector (mirrors run_project's idiom): a
    // `Box<dyn ReferenceProvider>` isn't `Clone`, and `VariantProjector`
    // requires `P: Clone`, so build the concrete `MultiFastaProvider`,
    // extract its owned cdot map, then share it through an `Arc`.
    let manifest_path = reference.join("manifest.json");
    let provider = MultiFastaProvider::from_manifest(&manifest_path)?;
    let cdot = provider
        .cdot_mapper()
        .cloned()
        .unwrap_or_else(CdotMapper::new);
    let provider: ArbitrateProvider = Arc::new(provider);
    let projector = VariantProjector::new(Projector::new(cdot), provider.clone());
    let resolver = CliFrameResolver(projector);

    let parsed = parse_hgvs(&variant)?;

    // Ferro's own normalized output for the input variant.
    //
    // `arbitrate` carries no `--error-mode`/`--ignore`/`--reject` flags, so
    // there is no user-supplied `ErrorConfig` to thread through here. Name the
    // lenient one explicitly rather than let a defaulting constructor supply it
    // silently (#1197): arbitration compares ferro's best-effort rendering
    // against another tool's, so refusing a repairable input would report a
    // disagreement that neither tool actually has.
    let config = NormalizeConfig::for_entry_point(
        parse_shuffle_direction("3prime"),
        ferro_hgvs::error_handling::ErrorConfig::lenient(),
    );
    let normalizer = Normalizer::with_config(provider.clone(), config);
    let ferro_output = normalizer.normalize(&parsed).map_err(|e| {
        format!("ferro failed to normalize {variant}: {e} (arbitration needs ferro's own output)")
    })?;
    let ferro_output = ferro_output.to_string();

    let other = fetch_other_result(
        &variant,
        other_output,
        other_tool,
        mutalyzer_url,
        &parsed,
        &provider,
    )?;

    let arbitration = arbitrate(&variant, &ferro_output, other, &provider, Some(&resolver))?;

    match format {
        "json" => {
            println!("{}", serde_json::to_string_pretty(&arbitration)?);
        }
        _ => {
            print!("{}", render_arbitration_text(&arbitration));
        }
    }

    Ok(())
}

/// Drive [`ferro_hgvs::arbitrate::bug_report::BugReport`] from a completed
/// `ferro arbitrate --format json` bundle to file a GitHub bug report.
///
/// Reads the `Arbitration` from `from_arbitration` (`-` = stdin), prints the
/// full rendered issue body so the user sees exactly what will be
/// submitted, then either opens a browser (after an explicit `y`/`yes`
/// confirmation) or prints the prefilled URL. Opening is attempted only when
/// `--no-open` was not given, `$BROWSER` is set, and stdin is an
/// interactive terminal — a headless/piped invocation always degrades to
/// print-only, per the design's "no `$BROWSER` -> print, never error" rule.
/// An over-length rendered body (`BugReportUrl::OverLength`) is written to a
/// `bug-report-<slug>.md` file in the system temp dir (never the cwd, to
/// avoid leaking cwd context) instead of being prefilled.
fn run_bug_report(
    from_arbitration: Option<&Path>,
    category: Option<&str>,
    notes: Option<&str>,
    no_open: bool,
    include_environment: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::arbitrate::bug_report::BugReport;
    use ferro_hgvs::arbitrate::{Arbitration, Verdict};

    let from_arbitration = from_arbitration.ok_or(
        "ferro bug-report requires --from-arbitration <path|-> in this version \
         (pipe in `ferro arbitrate --format json` output; there is no standalone fallback yet)",
    )?;

    let json = if from_arbitration == Path::new("-") {
        let mut buf = String::new();
        io::stdin().read_to_string(&mut buf)?;
        buf
    } else {
        std::fs::read_to_string(from_arbitration).map_err(|e| {
            format!(
                "failed to read arbitration bundle {}: {e}",
                from_arbitration.display()
            )
        })?
    };
    let arbitration: Arbitration = serde_json::from_str(&json).map_err(|e| {
        format!(
            "failed to parse arbitration bundle {} as JSON: {e}",
            from_arbitration.display()
        )
    })?;

    // An inconclusive arbitration never reached a compliance judgment (the
    // oracle/projection step itself errored, e.g. an unplaceable variant) —
    // there is nothing to accuse ferro or the other tool of, so refuse to
    // file a bug report rather than offer one for a non-verdict (#886).
    if arbitration.verdict == Verdict::Inconclusive {
        return Err(format!(
            "ferro bug-report: the arbitration bundle is Inconclusive (reason: {}); \
             nothing to report — a bug report requires a completed verdict",
            arbitration.reason.as_deref().unwrap_or("unknown")
        )
        .into());
    }

    // `category` is an informational hint (e.g. from the driving skill); the
    // bundle's own category is authoritative for the report, so a mismatch
    // is only a warning, never a hard error.
    if let Some(cat) = category {
        let bundle_cat = arbitration.category.to_string();
        if cat != bundle_cat {
            eprintln!(
                "warning: --category {cat} does not match the arbitration bundle's category \
                 `{bundle_cat}`; the bundle's category is what gets reported"
            );
        }
    }

    let report = BugReport {
        arbitration: &arbitration,
        notes,
        include_environment,
        ferro_version: env!("CARGO_PKG_VERSION"),
    };

    println!("{}", report.issue_title());
    println!();
    println!("{}", report.issue_body());
    println!();

    if !should_offer_to_open_browser(no_open) {
        return emit_bug_report_url(&report, false);
    }

    print!("Open a browser to file this issue? [y/N] ");
    io::stdout().flush()?;
    let mut answer = String::new();
    io::stdin().read_line(&mut answer)?;
    let answer = answer.trim().to_ascii_lowercase();
    let confirmed = answer == "y" || answer == "yes";
    if !confirmed {
        println!("Not opening a browser.");
    }
    emit_bug_report_url(&report, confirmed)
}

/// Resolve the report's prefilled new-issue URL — writing the over-length
/// fallback body to a temp-dir file first when the report doesn't fit in a
/// URL — then either open it in a browser (`try_open`) or just print it.
/// The URL (and, for the over-length case, the file path) is always
/// printed, whether or not opening was attempted or succeeded, so the user
/// has it either way.
fn emit_bug_report_url(
    report: &ferro_hgvs::arbitrate::bug_report::BugReport,
    try_open: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::arbitrate::bug_report::BugReportUrl;

    let url = match report.new_issue_url() {
        BugReportUrl::Url(url) => url,
        BugReportUrl::OverLength { url, body } => {
            let path = write_overlength_bug_report(report.arbitration, &body)?;
            println!(
                "The full report is too long to prefill in a URL; the full body was written to: {}",
                path.display()
            );
            url
        }
    };

    if try_open && try_open_browser(&url) {
        return Ok(());
    }
    if try_open {
        println!("Could not open a browser automatically; open this URL yourself:");
    }
    println!("{url}");
    Ok(())
}

/// Whether `run_bug_report` should attempt to open a browser at all: only
/// when the caller didn't pass `--no-open`, `$BROWSER` is set, and stdin is
/// an interactive terminal (so a `y`/`yes` confirmation can actually be
/// answered — a script piping `--from-arbitration -` has already redirected
/// stdin, so this naturally also skips confirmation in that case).
fn should_offer_to_open_browser(no_open: bool) -> bool {
    !no_open && std::env::var_os("BROWSER").is_some() && atty::is(atty::Stream::Stdin)
}

/// Best-effort browser launch via the platform opener (`open` on macOS,
/// `xdg-open` elsewhere). Any failure — missing binary, spawn error,
/// nonzero exit — degrades to `false` rather than propagating an error,
/// since opening a browser is a convenience, not the primary deliverable
/// (the caller always also prints the URL).
fn try_open_browser(url: &str) -> bool {
    let opener = if cfg!(target_os = "macos") {
        "open"
    } else {
        "xdg-open"
    };
    matches!(
        std::process::Command::new(opener).arg(url).status(),
        Ok(status) if status.success()
    )
}

/// Write an over-length bug-report body to `bug-report-<slug>.md` in the
/// system temp dir (never the cwd, to avoid leaking cwd context), returning
/// the path written.
fn write_overlength_bug_report(
    arbitration: &ferro_hgvs::arbitrate::Arbitration,
    body: &str,
) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let slug = slugify(&arbitration.input);
    let path = std::env::temp_dir().join(format!("bug-report-{slug}.md"));
    std::fs::write(&path, body)?;
    Ok(path)
}

/// Lowercase, filesystem-safe slug: runs of non-alphanumeric characters
/// collapse to a single `-`, with no leading/trailing `-`.
fn slugify(s: &str) -> String {
    let mut slug = String::with_capacity(s.len());
    let mut last_was_sep = true; // avoid a leading '-'
    for c in s.chars() {
        if c.is_ascii_alphanumeric() {
            slug.push(c.to_ascii_lowercase());
            last_was_sep = false;
        } else if !last_was_sep {
            slug.push('-');
            last_was_sep = true;
        }
    }
    while slug.ends_with('-') {
        slug.pop();
    }
    slug
}

#[allow(clippy::too_many_arguments)]
fn run_parse(
    variant: Option<&str>,
    input: Option<&PathBuf>,
    output: Option<&PathBuf>,
    format: &str,
    timing: Option<&PathBuf>,
    workers: usize,
    error_config: &ErrorConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::time::Instant;

    let preprocessor = error_config.preprocessor();

    // Create output writer - either file or stdout
    let mut writer: Box<dyn Write> = match output {
        Some(path) => Box::new(BufWriter::new(File::create(path)?)),
        None => Box::new(io::stdout()),
    };
    let mut error_count = 0usize;
    let mut success_count = 0usize;
    let mut total_count = 0usize;
    let start = Instant::now();

    // Per-line work, factored so it can run serially or in parallel: `out`
    // receives the success output (the order-critical stdout/file stream) and
    // `err` receives text-mode warnings. Mirrors `run_normalize`'s `process`.
    let process = |v: &str, out: &mut dyn Write, err: &mut dyn Write| -> Result<(), FerroError> {
        // Preprocess input according to error mode
        let mut preprocess_result = preprocessor.preprocess(v);
        if let Some(rejection) = preprocess_result.take_rejection_error() {
            // Preprocessing rejected the input: surface the phase's own
            // diagnosis (message, span, and structured code), not a summary
            // built from the — deliberately empty — warning list (#1162).
            return Err(rejection);
        }

        let parsed = parse_hgvs(&preprocess_result.preprocessed)?;
        match format {
            "json" => {
                let corrections: Vec<String> = preprocess_result
                    .warnings
                    .iter()
                    .map(|w| {
                        format!(
                            r#"{{"code":"{}","message":"{}"}}"#,
                            w.error_type.code(),
                            w.message.replace('"', "\\\"")
                        )
                    })
                    .collect();
                writeln!(
                    out,
                    r#"{{"input": "{}", "parsed": "{}", "status": "ok", "corrections": [{}]}}"#,
                    v,
                    parsed,
                    corrections.join(",")
                )
                .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }
            _ => {
                // Emit warnings to the stderr buffer if any
                for warning in &preprocess_result.warnings {
                    writeln!(
                        err,
                        "warning[{}]: {}",
                        warning.error_type.code(),
                        warning.message
                    )
                    .map_err(|e| FerroError::Io { msg: e.to_string() })?;
                }
                writeln!(out, "{} -> {}", v, parsed)
                    .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }
        }
        Ok(())
    };

    if let Some(v) = variant {
        // Single-variant path: unchanged (no line context, errors via
        // `output_error`). Warnings go straight to the live stderr.
        total_count += 1;
        let mut stderr = io::stderr();
        if let Err(e) = process(v, &mut writer, &mut stderr) {
            output_error(v, &e, format)?;
            error_count += 1;
        } else {
            success_count += 1;
        }
    } else {
        // Batch path (file or stdin): same order-preserving, optionally parallel
        // driver as `run_normalize`, producing byte-identical output to serial.
        let format_owned: OutputFormat = format.parse().unwrap_or_default();
        let item = |line_num: usize, line: &str| -> LineResult {
            let is_first = line_num == 0;
            let Some(variant_str) = process_input_line(line, is_first) else {
                return LineResult::skipped();
            };
            let mut out = Vec::new();
            let mut err = Vec::new();
            let status = match process(variant_str, &mut out, &mut err) {
                Ok(()) => LineStatus::Ok,
                Err(e) => {
                    let _ = cli_output_error_with_context(
                        &mut err,
                        variant_str,
                        &e,
                        format_owned,
                        Some(line_num + 1),
                    );
                    LineStatus::Err
                }
            };
            LineResult {
                stdout: out,
                stderr: err,
                status,
            }
        };

        let (t, s, er) = if let Some(input_path) = input {
            let file = std::fs::File::open(input_path)?;
            run_batch(
                io::BufReader::new(file).lines(),
                &mut writer,
                workers,
                &item,
            )?
        } else {
            let stdin = io::stdin();
            let reader = stdin.lock();
            run_batch(reader.lines(), &mut writer, workers, &item)?
        };
        total_count += t;
        success_count += s;
        error_count += er;
    }

    let elapsed = start.elapsed();

    // Write timing info if requested
    if let Some(timing_path) = timing {
        let timing_info =
            ferro_hgvs::commands::TimingInfo::new(total_count, success_count, elapsed);
        ferro_hgvs::commands::write_timing(&timing_info, timing_path)?;
    }

    if error_count > 0 {
        Err(format!("{} variant(s) failed to parse", error_count).into())
    } else {
        Ok(())
    }
}

/// Byte threshold above which embedding genomic sequence into a `transcripts.json`
/// is flagged as heavy — the file becomes large to store and parse, and a
/// genome-scale reference is better served by `ferro prepare` + `from_manifest`.
/// Warn on stderr when `--emit-genomic-sequences` embeds a large amount of
/// genomic sequence, so a user does not silently produce a multi-hundred-MB
/// `transcripts.json` from whole-chromosome (or whole-genome) coordinates.
fn warn_if_genomic_sequences_large(total_bytes: u64) {
    use ferro_hgvs::reference::annotation::convert::LARGE_GENOMIC_SEQUENCES_WARN_BYTES;
    if total_bytes > LARGE_GENOMIC_SEQUENCES_WARN_BYTES {
        eprintln!(
            "warning: --emit-genomic-sequences embedded {:.1} MB of genomic sequence; the \
             transcripts.json will be large to store and slow to parse. For a genome-scale \
             reference prefer `ferro prepare` + `Normalizer.from_manifest(...)`.",
            total_bytes as f64 / 1_000_000.0
        );
    }
}

/// Split a comma-separated CLI filter value (`--transcripts`, `--genes`) into
/// trimmed items. Empty items are preserved so the behavior matches the
/// previous inline split; [`convert_gff`](ferro_hgvs::reference::annotation::convert_gff)
/// trims again defensively.
fn split_comma_list(s: &str) -> Vec<String> {
    s.split(',').map(|t| t.trim().to_string()).collect()
}

#[allow(clippy::too_many_arguments)]
fn run_convert_gff(
    gff: &Path,
    fasta: Option<&PathBuf>,
    output: Option<&PathBuf>,
    build: &str,
    mane_only: bool,
    transcripts: Option<&str>,
    genes: Option<&str>,
    strict: bool,
    silent: bool,
    no_validate_fasta: bool,
    emit_genomic_sequences: bool,
    diagnostics_json: Option<&PathBuf>,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::reference::annotation::{convert_gff, ConvertGffConfig};

    // Build the conversion config, mirroring the CLI flags one-to-one. The
    // serialization itself lives in `ferro_hgvs::reference::annotation::convert_gff`
    // so this command and the Python binding produce byte-identical output.
    let config = ConvertGffConfig {
        genome_build: parse_genome_build(build),
        mane_only,
        transcripts: transcripts.map(split_comma_list),
        genes: genes.map(split_comma_list),
        strict,
        silent,
        no_validate_fasta,
        emit_genomic_sequences,
    };

    let outcome = convert_gff(gff, fasta.map(PathBuf::as_path), &config)?;

    // Surface the loader summary line, as before.
    eprintln!("{}", outcome.report.summary_line());

    // Write diagnostics JSON if requested.
    if let Some(path) = diagnostics_json {
        let json = serde_json::to_string_pretty(&outcome.report.sample_diagnostics)
            .map_err(|e| format!("serialize diagnostics: {}", e))?;
        std::fs::write(path, json).map_err(|e| format!("write {}: {}", path.display(), e))?;
    }

    // Reproduce the --emit-genomic-sequences stderr warnings from the outcome.
    if outcome.emit_requested_no_placement {
        eprintln!(
            "warning: --emit-genomic-sequences was set but no emitted transcript has a \
             genomic placement; genomic_sequences not written (reference stays \
             transcript-only)"
        );
    }
    warn_if_genomic_sequences_large(outcome.emitted_genomic_bytes);

    // Write output.
    let mut writer: Box<dyn Write> = if let Some(out_path) = output {
        Box::new(std::fs::File::create(out_path)?)
    } else {
        Box::new(io::stdout())
    };

    writeln!(writer, "{}", serde_json::to_string_pretty(&outcome.json)?)?;

    Ok(())
}

fn run_liftover(
    position: Option<&str>,
    input: Option<&PathBuf>,
    chain: &PathBuf,
    reverse_chain: Option<&PathBuf>,
    from: &str,
    to: &str,
    format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::liftover::{ChainFile, Liftover};

    let source_build = parse_genome_build(from);
    let target_build = parse_genome_build(to);

    // Load chain files
    let forward_chain = ChainFile::from_file(chain)?;
    let reverse_chain_file = if let Some(rev_path) = reverse_chain {
        ChainFile::from_file(rev_path)?
    } else {
        // Create empty chain for one-way liftover
        ChainFile::new()
    };

    let liftover = Liftover::new(forward_chain, reverse_chain_file);

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let process_position =
        |pos_str: &str, handle: &mut io::StdoutLock| -> Result<(), Box<dyn std::error::Error>> {
            // Parse position: chr1:12345 or just chr1 12345 in file
            let (contig, pos) = if pos_str.contains(':') {
                // Try HGVS format first: NC_000001.10:g.12345A>G
                if pos_str.contains(":g.") || pos_str.contains(":m.") {
                    // Parse as HGVS
                    let parts: Vec<&str> = pos_str.split(':').collect();
                    if parts.len() >= 2 {
                        let accession = parts[0];
                        let coord_part = parts[1];
                        // Extract position from g.12345...
                        if let Some(pos_start) = coord_part.find(|c: char| c.is_ascii_digit()) {
                            let pos_end = coord_part[pos_start..]
                                .find(|c: char| !c.is_ascii_digit())
                                .map(|i| pos_start + i)
                                .unwrap_or(coord_part.len());
                            let pos: u64 = coord_part[pos_start..pos_end].parse()?;
                            (accession.to_string(), pos)
                        } else {
                            return Err(format!("Invalid HGVS position: {}", pos_str).into());
                        }
                    } else {
                        return Err(format!("Invalid position format: {}", pos_str).into());
                    }
                } else {
                    // Simple chr:pos format
                    let parts: Vec<&str> = pos_str.split(':').collect();
                    if parts.len() == 2 {
                        (parts[0].to_string(), parts[1].parse()?)
                    } else {
                        return Err(format!("Invalid position format: {}", pos_str).into());
                    }
                }
            } else {
                return Err(format!("Invalid position format: {}", pos_str).into());
            };

            // Perform liftover
            match liftover.lift(source_build, target_build, &contig, pos) {
                Ok(result) => match format {
                    "json" => {
                        writeln!(
                            handle,
                            r#"{{"source": "{}:{}", "target": "{}:{}", "chain_id": {}}}"#,
                            contig, pos, result.target_contig, result.target_pos, result.chain_id
                        )?;
                    }
                    _ => {
                        writeln!(
                            handle,
                            "{}:{} -> {}:{}",
                            contig, pos, result.target_contig, result.target_pos
                        )?;
                    }
                },
                Err(e) => match format {
                    "json" => {
                        writeln!(
                            handle,
                            r#"{{"source": "{}:{}", "error": "{}"}}"#,
                            contig, pos, e
                        )?;
                    }
                    _ => {
                        writeln!(handle, "{}:{} -> ERROR: {}", contig, pos, e)?;
                    }
                },
            }
            Ok(())
        };

    if let Some(pos) = position {
        process_position(pos, &mut handle)?;
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                process_position(trimmed, &mut handle)?;
            }
        }
    } else {
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                process_position(trimmed, &mut handle)?;
            }
        }
    }

    Ok(())
}

fn run_describe(
    reference: Option<&str>,
    observed: Option<&str>,
    accession: Option<&str>,
    input: Option<&PathBuf>,
    format: &str,
    detect_duplications: bool,
    detect_inversions: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::extractor::{DescriptionExtractor, ExtractorConfig};

    let config = ExtractorConfig {
        detect_duplications,
        detect_inversions,
        ..Default::default()
    };
    let extractor = DescriptionExtractor::new(config);

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let process_seqs = |ref_seq: &str,
                        obs_seq: &str,
                        accession: Option<&str>,
                        handle: &mut io::StdoutLock|
     -> Result<(), Box<dyn std::error::Error>> {
        let result = if let Some(acc) = accession {
            extractor.extract_with_accession(acc, ref_seq, obs_seq)?
        } else {
            extractor.extract(ref_seq, obs_seq)?
        };

        if result.variants.is_empty() {
            match format {
                "json" => {
                    writeln!(
                        handle,
                        r#"{{"reference": "{}", "observed": "{}", "hgvs": [], "message": "No variants detected"}}"#,
                        ref_seq, obs_seq
                    )?;
                }
                _ => {
                    writeln!(handle, "No variants detected (sequences are identical)")?;
                }
            }
        } else {
            match format {
                "json" => {
                    let hgvs_array: Vec<String> = result
                        .hgvs_strings
                        .iter()
                        .map(|s| format!("\"{}\"", s))
                        .collect();
                    writeln!(
                        handle,
                        r#"{{"reference_length": {}, "observed_length": {}, "variants": {}, "hgvs": [{}]}}"#,
                        result.reference_length,
                        result.observed_length,
                        result.variants.len(),
                        hgvs_array.join(", ")
                    )?;
                }
                _ => {
                    for hgvs in &result.hgvs_strings {
                        writeln!(handle, "{}", hgvs)?;
                    }
                }
            }
        }
        Ok(())
    };

    if let (Some(ref_seq), Some(obs_seq)) = (reference, observed) {
        process_seqs(ref_seq, obs_seq, accession, &mut handle)?;
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            // Tab-separated: accession\treference\tobserved or just reference\tobserved
            let parts: Vec<&str> = trimmed.split('\t').collect();
            match parts.len() {
                2 => process_seqs(parts[0], parts[1], None, &mut handle)?,
                3 => process_seqs(parts[1], parts[2], Some(parts[0]), &mut handle)?,
                _ => {
                    writeln!(
                        handle,
                        "ERROR: Invalid line format (expected 2 or 3 tab-separated fields): {}",
                        trimmed
                    )?;
                }
            }
        }
    } else {
        return Err("Either --reference and --observed or --input must be provided".into());
    }

    Ok(())
}

fn run_effect(
    variant: Option<&str>,
    input: Option<&PathBuf>,
    format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::effect::{Consequence, EffectPredictor, Impact, ProteinEffect};

    let predictor = EffectPredictor::new();
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let process_variant = |v: &str,
                           handle: &mut io::StdoutLock|
     -> Result<(), Box<dyn std::error::Error>> {
        // Try parsing as simple protein shorthand first: V600E or Val600Glu
        let effect = if let Some((ref_aa, pos, alt_aa)) = parse_simple_protein_variant(v) {
            predictor.classify_amino_acid_change(&ref_aa, &alt_aa, pos)
        } else if let Ok((ref_aa, pos, alt_aa)) = parse_protein_change(v) {
            predictor.classify_amino_acid_change(&ref_aa, &alt_aa, pos)
        } else if let Ok(parsed) = parse_hgvs(v) {
            // Try HGVS parsing
            match &parsed {
                ferro_hgvs::HgvsVariant::Protein(pv) => {
                    if let Some(edit) = pv.loc_edit.edit.inner() {
                        match edit {
                            ferro_hgvs::hgvs::edit::ProteinEdit::Substitution {
                                reference,
                                alternative,
                            } => {
                                let pos = pv
                                    .loc_edit
                                    .location
                                    .start
                                    .inner()
                                    .map(|s| s.number)
                                    .unwrap_or(1);
                                predictor.classify_amino_acid_change(reference, alternative, pos)
                            }
                            ferro_hgvs::hgvs::edit::ProteinEdit::Frameshift { .. } => {
                                predictor.classify_indel(1, 0)
                            }
                            ferro_hgvs::hgvs::edit::ProteinEdit::Deletion { .. } => {
                                predictor.classify_indel(3, 0)
                            }
                            ferro_hgvs::hgvs::edit::ProteinEdit::Insertion { sequence }
                                if sequence.introduces_stop() =>
                            {
                                // `insTer<n>` / `ins*<n>` (or a literal ending in
                                // `Ter`) introduces a premature stop.
                                ProteinEffect {
                                    consequences: vec![Consequence::StopGained],
                                    impact: Consequence::StopGained.impact(),
                                    amino_acid_change: None,
                                    intronic_offset: None,
                                }
                            }
                            ferro_hgvs::hgvs::edit::ProteinEdit::Insertion { .. } => {
                                predictor.classify_indel(0, 3)
                            }
                            ferro_hgvs::hgvs::edit::ProteinEdit::Delins { sequence }
                                if sequence.ends_with_stop() =>
                            {
                                // A delins whose inserted sequence ends in `Ter`
                                // (`delinsLeuTer`) truncates the protein — a
                                // premature stop, mirroring the insertion arm.
                                ProteinEffect {
                                    consequences: vec![Consequence::StopGained],
                                    impact: Consequence::StopGained.impact(),
                                    amino_acid_change: None,
                                    intronic_offset: None,
                                }
                            }
                            _ => ProteinEffect {
                                consequences: vec![Consequence::ProteinAlteringVariant],
                                impact: Impact::Moderate,
                                amino_acid_change: None,
                                intronic_offset: None,
                            },
                        }
                    } else {
                        ProteinEffect {
                            consequences: vec![Consequence::ProteinAlteringVariant],
                            impact: Impact::Moderate,
                            amino_acid_change: None,
                            intronic_offset: None,
                        }
                    }
                }
                ferro_hgvs::HgvsVariant::Cds(cv) => {
                    // Check for intronic offsets
                    if let Some(start_pos) = cv.loc_edit.location.start.inner() {
                        if let Some(offset) = start_pos.offset {
                            predictor.classify_splice_variant(offset)
                        } else {
                            ProteinEffect {
                                consequences: vec![Consequence::CodingSequenceVariant],
                                impact: Impact::Modifier,
                                amino_acid_change: None,
                                intronic_offset: None,
                            }
                        }
                    } else {
                        ProteinEffect {
                            consequences: vec![Consequence::CodingSequenceVariant],
                            impact: Impact::Modifier,
                            amino_acid_change: None,
                            intronic_offset: None,
                        }
                    }
                }
                _ => {
                    return Err(
                        format!("Variant type not supported for effect prediction: {}", v).into(),
                    );
                }
            }
        } else {
            return Err(format!("Could not parse variant: {}", v).into());
        };

        match format {
            "json" => {
                let consequences: Vec<String> = effect
                    .consequences
                    .iter()
                    .map(|c| format!("\"{}\"", c.so_term()))
                    .collect();
                writeln!(
                    handle,
                    r#"{{"variant": "{}", "consequences": [{}], "impact": "{}", "is_high_impact": {}}}"#,
                    v,
                    consequences.join(", "),
                    effect.impact,
                    effect.is_high_impact()
                )?;
            }
            _ => {
                let consequences: Vec<String> = effect
                    .consequences
                    .iter()
                    .map(|c| c.so_term().to_string())
                    .collect();
                writeln!(
                    handle,
                    "{}: {} ({})",
                    v,
                    consequences.join(", "),
                    effect.impact
                )?;
            }
        }
        Ok(())
    };

    if let Some(v) = variant {
        process_variant(v, &mut handle)?;
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                if let Err(e) = process_variant(trimmed, &mut handle) {
                    eprintln!("ERROR: {} - {}", trimmed, e);
                }
            }
        }
    } else {
        let stdin = io::stdin();
        for line in stdin.lock().lines() {
            let line = line?;
            let trimmed = line.trim();
            if !trimmed.is_empty() && !trimmed.starts_with('#') {
                if let Err(e) = process_variant(trimmed, &mut handle) {
                    eprintln!("ERROR: {} - {}", trimmed, e);
                }
            }
        }
    }

    Ok(())
}

/// Parse simple protein variant like V600E or Val600Glu
fn parse_simple_protein_variant(
    v: &str,
) -> Option<(
    ferro_hgvs::hgvs::location::AminoAcid,
    u64,
    ferro_hgvs::hgvs::location::AminoAcid,
)> {
    use ferro_hgvs::hgvs::location::AminoAcid;

    // Try single-letter: V600E
    if v.len() >= 3 {
        let first_char = v.chars().next()?;
        if let Some(ref_aa) = AminoAcid::from_one_letter(first_char) {
            // Find number
            let num_start = 1;
            let num_end = v[1..].find(|c: char| !c.is_ascii_digit())? + 1;
            let pos: u64 = v[num_start..num_end].parse().ok()?;

            let last_char = v.chars().last()?;
            if let Some(alt_aa) = AminoAcid::from_one_letter(last_char) {
                return Some((ref_aa, pos, alt_aa));
            }
        }
    }

    None
}

fn run_backtranslate(
    variant: Option<&str>,
    ref_aa: Option<&str>,
    alt_aa: Option<&str>,
    input: Option<&PathBuf>,
    format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::backtranslate::Backtranslator;
    use ferro_hgvs::hgvs::location::AminoAcid;

    let bt = Backtranslator::standard();
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let process_change = |ref_aa: &AminoAcid,
                          alt_aa: &AminoAcid,
                          label: &str,
                          handle: &mut io::StdoutLock|
     -> Result<(), Box<dyn std::error::Error>> {
        let changes = if *alt_aa == AminoAcid::Ter {
            bt.backtranslate_to_stop(ref_aa)
        } else if *ref_aa == AminoAcid::Ter {
            bt.backtranslate_stop_loss(alt_aa)
        } else {
            bt.backtranslate_substitution(ref_aa, alt_aa)
        };

        match format {
            "json" => {
                let codon_changes: Vec<String> = changes
                    .iter()
                    .map(|c| {
                        format!(
                            r#"{{"ref_codon": "{}", "alt_codon": "{}", "position": {}}}"#,
                            c.ref_codon,
                            c.alt_codon,
                            c.changed_positions.first().unwrap_or(&0)
                        )
                    })
                    .collect();
                writeln!(
                    handle,
                    r#"{{"variant": "{}", "ref_aa": "{}", "alt_aa": "{}", "codon_changes": [{}]}}"#,
                    label,
                    ref_aa,
                    alt_aa,
                    codon_changes.join(", ")
                )?;
            }
            _ => {
                writeln!(handle, "{} ({} -> {}):", label, ref_aa, alt_aa)?;
                for change in &changes {
                    writeln!(
                        handle,
                        "  {} -> {} (position {})",
                        change.ref_codon,
                        change.alt_codon,
                        change
                            .changed_positions
                            .iter()
                            .map(|p| p.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    )?;
                }
                if changes.is_empty() {
                    writeln!(handle, "  (no single-nucleotide changes possible)")?;
                }
            }
        }
        Ok(())
    };

    // If explicit amino acids provided
    if let (Some(ref_str), Some(alt_str)) = (ref_aa, alt_aa) {
        let ref_aa = AminoAcid::from_one_letter(ref_str.chars().next().ok_or("Empty ref_aa")?)
            .or_else(|| parse_three_letter_aa(ref_str))
            .ok_or_else(|| format!("Invalid amino acid: {}", ref_str))?;
        let alt_aa = AminoAcid::from_one_letter(alt_str.chars().next().ok_or("Empty alt_aa")?)
            .or_else(|| parse_three_letter_aa(alt_str))
            .ok_or_else(|| format!("Invalid amino acid: {}", alt_str))?;
        process_change(
            &ref_aa,
            &alt_aa,
            &format!("{}->{}", ref_str, alt_str),
            &mut handle,
        )?;
    } else if let Some(v) = variant {
        // Parse variant like p.Val600Glu or V600E
        let (ref_aa, _pos, alt_aa) = parse_protein_change(v)?;
        process_change(&ref_aa, &alt_aa, v, &mut handle)?;
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            match parse_protein_change(trimmed) {
                Ok((ref_aa, _pos, alt_aa)) => {
                    process_change(&ref_aa, &alt_aa, trimmed, &mut handle)?;
                }
                Err(e) => {
                    eprintln!("ERROR: {} - {}", trimmed, e);
                }
            }
        }
    } else {
        return Err(
            "Either --variant, (--ref-aa and --alt-aa), or --input must be provided".into(),
        );
    }

    Ok(())
}

/// Parse 3-letter amino acid code
fn parse_three_letter_aa(s: &str) -> Option<ferro_hgvs::hgvs::location::AminoAcid> {
    use ferro_hgvs::hgvs::location::AminoAcid;

    let normalized: String = s
        .chars()
        .enumerate()
        .map(|(i, c)| {
            if i == 0 {
                c.to_ascii_uppercase()
            } else {
                c.to_ascii_lowercase()
            }
        })
        .collect();

    match normalized.as_str() {
        "Ala" => Some(AminoAcid::Ala),
        "Arg" => Some(AminoAcid::Arg),
        "Asn" => Some(AminoAcid::Asn),
        "Asp" => Some(AminoAcid::Asp),
        "Cys" => Some(AminoAcid::Cys),
        "Gln" => Some(AminoAcid::Gln),
        "Glu" => Some(AminoAcid::Glu),
        "Gly" => Some(AminoAcid::Gly),
        "His" => Some(AminoAcid::His),
        "Ile" => Some(AminoAcid::Ile),
        "Leu" => Some(AminoAcid::Leu),
        "Lys" => Some(AminoAcid::Lys),
        "Met" => Some(AminoAcid::Met),
        "Phe" => Some(AminoAcid::Phe),
        "Pro" => Some(AminoAcid::Pro),
        "Ser" => Some(AminoAcid::Ser),
        "Thr" => Some(AminoAcid::Thr),
        "Trp" => Some(AminoAcid::Trp),
        "Tyr" => Some(AminoAcid::Tyr),
        "Val" => Some(AminoAcid::Val),
        "Ter" | "*" => Some(AminoAcid::Ter),
        "Sec" => Some(AminoAcid::Sec),
        "Xaa" | "X" => Some(AminoAcid::Xaa),
        _ => None,
    }
}

/// Parse protein change like p.Val600Glu or V600E
fn parse_protein_change(
    v: &str,
) -> Result<
    (
        ferro_hgvs::hgvs::location::AminoAcid,
        u64,
        ferro_hgvs::hgvs::location::AminoAcid,
    ),
    Box<dyn std::error::Error>,
> {
    use ferro_hgvs::hgvs::location::AminoAcid;

    // Remove p. prefix if present
    let v = v.strip_prefix("p.").unwrap_or(v);

    // Try 3-letter code first: Val600Glu
    if v.len() >= 7 {
        if let Some(ref_aa) = parse_three_letter_aa(&v[..3]) {
            let rest = &v[3..];
            let num_end = rest
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(rest.len());
            if num_end > 0 {
                let pos: u64 = rest[..num_end].parse()?;
                let alt_part = &rest[num_end..];
                if alt_part.len() >= 3 {
                    if let Some(alt_aa) = parse_three_letter_aa(&alt_part[..3]) {
                        return Ok((ref_aa, pos, alt_aa));
                    }
                } else if alt_part == "*" {
                    return Ok((ref_aa, pos, AminoAcid::Ter));
                }
            }
        }
    }

    // Try single-letter: V600E
    if v.len() >= 3 {
        let first = v.chars().next().ok_or("Empty variant")?;
        if let Some(ref_aa) = AminoAcid::from_one_letter(first) {
            let rest = &v[1..];
            let num_end = rest
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(rest.len());
            if num_end > 0 {
                let pos: u64 = rest[..num_end].parse()?;
                let alt_char = rest[num_end..].chars().next().ok_or("No alt amino acid")?;
                if let Some(alt_aa) = AminoAcid::from_one_letter(alt_char) {
                    return Ok((ref_aa, pos, alt_aa));
                }
            }
        }
    }

    Err(format!("Could not parse protein change: {}", v).into())
}

#[allow(clippy::too_many_arguments)]
fn run_generate(
    accession: &str,
    coord: &str,
    variant_type: &str,
    pos: i64,
    end: Option<i64>,
    ref_base: Option<&str>,
    alt_base: Option<&str>,
    offset: Option<i64>,
    repeat_count: Option<u32>,
    ref_aa: Option<&str>,
    alt_aa: Option<&str>,
    format: &str,
    validate: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    // Build position string
    let position_str = if let Some(off) = offset {
        if off >= 0 {
            format!("{}+{}", pos, off)
        } else {
            format!("{}{}", pos, off)
        }
    } else {
        pos.to_string()
    };

    // Build interval if end position provided
    let interval_str = if let Some(e) = end {
        let end_str = if let Some(off) = offset {
            if off >= 0 {
                format!("{}+{}", e, off)
            } else {
                format!("{}{}", e, off)
            }
        } else {
            e.to_string()
        };
        format!("{}_{}", position_str, end_str)
    } else {
        position_str
    };

    // Build the HGVS description based on variant type and coordinate system
    let hgvs = if coord == "p" {
        // Protein variant
        let ref_aa_str = ref_aa.unwrap_or("?");
        let alt_aa_str = alt_aa.unwrap_or("?");

        match variant_type {
            "sub" => format!("{}:p.{}{}{}", accession, ref_aa_str, pos, alt_aa_str),
            "del" => {
                if let Some(e) = end {
                    format!(
                        "{}:p.{}{}_{}{}del",
                        accession, ref_aa_str, pos, alt_aa_str, e
                    )
                } else {
                    format!("{}:p.{}{}del", accession, ref_aa_str, pos)
                }
            }
            "ins" => {
                let end_pos = end.unwrap_or(pos + 1);
                format!(
                    "{}:p.{}{}_{}{}ins{}",
                    accession, ref_aa_str, pos, ref_aa_str, end_pos, alt_aa_str
                )
            }
            "dup" => {
                if let Some(e) = end {
                    format!(
                        "{}:p.{}{}_{}{}dup",
                        accession, ref_aa_str, pos, alt_aa_str, e
                    )
                } else {
                    format!("{}:p.{}{}dup", accession, ref_aa_str, pos)
                }
            }
            "delins" => {
                if let Some(e) = end {
                    format!(
                        "{}:p.{}{}_{}{}delins{}",
                        accession, ref_aa_str, pos, ref_aa_str, e, alt_aa_str
                    )
                } else {
                    format!("{}:p.{}{}delins{}", accession, ref_aa_str, pos, alt_aa_str)
                }
            }
            _ => return Err(format!("Unsupported protein variant type: {}", variant_type).into()),
        }
    } else {
        // Nucleotide variant (g, c, n, r, m)
        let ref_str = ref_base.unwrap_or("");
        let alt_str = alt_base.unwrap_or("");

        match variant_type {
            "sub" => {
                if ref_str.is_empty() || alt_str.is_empty() {
                    return Err("Substitution requires --ref-base and --alt-base".into());
                }
                format!(
                    "{}:{}.{}{}>{}",
                    accession, coord, interval_str, ref_str, alt_str
                )
            }
            "del" => {
                if ref_str.is_empty() {
                    format!("{}:{}.{}del", accession, coord, interval_str)
                } else {
                    format!("{}:{}.{}del{}", accession, coord, interval_str, ref_str)
                }
            }
            "ins" => {
                if alt_str.is_empty() {
                    return Err("Insertion requires --alt-base".into());
                }
                // Insertion needs two flanking positions
                let end_pos = end.unwrap_or(pos + 1);
                format!("{}:{}.{}_{}ins{}", accession, coord, pos, end_pos, alt_str)
            }
            "dup" => {
                if ref_str.is_empty() {
                    format!("{}:{}.{}dup", accession, coord, interval_str)
                } else {
                    format!("{}:{}.{}dup{}", accession, coord, interval_str, ref_str)
                }
            }
            "delins" => {
                if alt_str.is_empty() {
                    return Err("Delins requires --alt-base".into());
                }
                format!("{}:{}.{}delins{}", accession, coord, interval_str, alt_str)
            }
            "inv" => {
                format!("{}:{}.{}inv", accession, coord, interval_str)
            }
            "repeat" => {
                let count = repeat_count.unwrap_or(2);
                if ref_str.is_empty() {
                    format!("{}:{}.{}[{}]", accession, coord, interval_str, count)
                } else {
                    format!(
                        "{}:{}.{}{}[{}]",
                        accession, coord, interval_str, ref_str, count
                    )
                }
            }
            _ => return Err(format!("Unsupported variant type: {}", variant_type).into()),
        }
    };

    // Optionally validate by parsing
    let is_valid = if validate {
        match parse_hgvs(&hgvs) {
            Ok(_) => true,
            Err(e) => {
                eprintln!("Warning: Generated HGVS failed validation: {}", e);
                false
            }
        }
    } else {
        true
    };

    // Output the result
    match format {
        "json" => {
            writeln!(
                handle,
                r#"{{"hgvs": "{}", "accession": "{}", "coord": "{}", "type": "{}", "position": {}, "valid": {}}}"#,
                hgvs, accession, coord, variant_type, pos, is_valid
            )?;
        }
        _ => {
            if is_valid {
                writeln!(handle, "{}", hgvs)?;
            } else {
                writeln!(handle, "{} (validation failed)", hgvs)?;
            }
        }
    }

    Ok(())
}

fn run_extract_hgvs(
    input_path: &PathBuf,
    output_path: &PathBuf,
    field_name: &str,
    hgvsc_idx: usize,
    hgvsp_idx: usize,
    search_prefix: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    eprintln!("Reading {}...", input_path.display());
    eprintln!(
        "Looking for {} field, HGVSc at index {}, HGVSp at index {}",
        field_name, hgvsc_idx, hgvsp_idx
    );
    if let Some(prefix) = search_prefix {
        eprintln!("Filtering to patterns starting with: {}", prefix);
    }

    let mut seen: HashSet<String> = HashSet::with_capacity(50_000_000);
    let field_prefix = format!("{}=", field_name);

    // Open input
    let file = File::open(input_path)?;
    let input_str = input_path.to_string_lossy();
    let reader: Box<dyn BufRead> = if input_str.ends_with(".gz") || input_str.ends_with(".bgz") {
        Box::new(BufReader::with_capacity(
            1024 * 1024,
            MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(1024 * 1024, file))
    };

    // Open output
    let out_file = File::create(output_path)?;
    let mut writer = BufWriter::with_capacity(1024 * 1024, out_file);

    let mut lines_processed: u64 = 0;
    let mut patterns_written: usize = 0;

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };

        // Skip headers
        if line.starts_with('#') {
            if line.contains(&format!("ID={}", field_name)) {
                eprintln!(
                    "Found {} header: {}...",
                    field_name,
                    &line[..line.len().min(150)]
                );
            }
            continue;
        }

        lines_processed += 1;
        #[allow(clippy::manual_is_multiple_of)]
        if lines_processed % 500_000 == 0 {
            eprintln!(
                "  Processed {} lines, written {} unique patterns...",
                lines_processed, patterns_written
            );
        }

        // Parse variant line - find INFO field (column 8, 0-indexed 7)
        let mut field_count = 0;
        let mut info_start = 0;
        let mut info_end = line.len();

        for (i, c) in line.char_indices() {
            if c == '\t' {
                field_count += 1;
                if field_count == 7 {
                    info_start = i + 1;
                } else if field_count == 8 {
                    info_end = i;
                    break;
                }
            }
        }

        if field_count < 7 {
            continue;
        }

        let info = &line[info_start..info_end];

        // Find the VEP/CSQ field
        for item in info.split(';') {
            if item.starts_with(&field_prefix) {
                let data = &item[field_prefix.len()..];

                // Process each transcript annotation
                for transcript in data.split(',') {
                    let values: Vec<&str> = transcript.split('|').collect();

                    // Extract HGVSc
                    if values.len() > hgvsc_idx {
                        let hgvsc = values[hgvsc_idx];
                        if !hgvsc.is_empty() && hgvsc.contains(':') {
                            // Strip VEP's appended protein consequence (p.=), (p.?)
                            let clean = if let Some(pos) = hgvsc.find("(p.") {
                                &hgvsc[..pos]
                            } else {
                                hgvsc
                            };

                            // Check search prefix
                            let matches =
                                search_prefix.map(|p| clean.starts_with(p)).unwrap_or(true);

                            if matches && !seen.contains(clean) {
                                seen.insert(clean.to_string());
                                writeln!(writer, "{}", clean)?;
                                patterns_written += 1;
                            }
                        }
                    }

                    // Extract HGVSp
                    if values.len() > hgvsp_idx {
                        let hgvsp_raw = values[hgvsp_idx];
                        if !hgvsp_raw.is_empty() && hgvsp_raw.contains(':') {
                            // URL decode %3D -> = and %3E -> >
                            let hgvsp_decoded = hgvsp_raw.replace("%3D", "=").replace("%3E", ">");

                            // Strip VEP's appended protein consequence (p.=), (p.?)
                            let hgvsp = if let Some(pos) = hgvsp_decoded.find("(p.") {
                                &hgvsp_decoded[..pos]
                            } else {
                                &hgvsp_decoded
                            };

                            // Check search prefix
                            let matches =
                                search_prefix.map(|p| hgvsp.starts_with(p)).unwrap_or(true);

                            if matches && !seen.contains(hgvsp) {
                                seen.insert(hgvsp.to_string());
                                writeln!(writer, "{}", hgvsp)?;
                                patterns_written += 1;
                            }
                        }
                    }
                }
                break;
            }
        }
    }

    writer.flush()?;
    eprintln!("Done! Extracted {} unique HGVS patterns", patterns_written);

    Ok(())
}

/// Build an ErrorConfig from CLI flags, with optional config file loading.
///
/// Config file priority:
/// 1. CLI flags (highest priority)
/// 2. `.ferro.toml` in current directory
/// 3. `~/.config/ferro/config.toml`
/// 4. Built-in defaults (Strict mode)
fn build_error_config(mode: &str, ignore: &[String], reject: &[String]) -> ErrorConfig {
    use ferro_hgvs::error_handling::ErrorOverride;

    // Try to load config file
    if let Some(file_config) = FerroConfig::load() {
        // Merge file config with CLI args (CLI takes precedence)
        // Only use CLI mode if it's not the default "strict"
        let cli_mode = if mode != "strict" { Some(mode) } else { None };
        return file_config.merge_with_cli(cli_mode, ignore, reject);
    }

    // No config file found, use CLI flags only
    let base_mode = match mode {
        "lenient" => ErrorMode::Lenient,
        "silent" => ErrorMode::Silent,
        _ => ErrorMode::Strict,
    };

    let mut config = ErrorConfig::new(base_mode);

    // Apply ignore overrides (silent correct)
    for code in ignore {
        apply_override(
            &mut config,
            code,
            ErrorOverride::SilentCorrect,
            "ignore",
            OverrideSource::Cli,
        );
    }

    // Apply reject overrides
    for code in reject {
        apply_override(
            &mut config,
            code,
            ErrorOverride::Reject,
            "reject",
            OverrideSource::Cli,
        );
    }

    config
}

// `code_to_error_type` was replaced by `ErrorType::from_code`
// (`src/error_handling/types.rs`) so the CLI's `--ignore` / `--reject`
// parser cannot drift away from the registered W-code set the way an
// independent hard-coded table did.

/// Run the explain command.
fn run_explain(
    code: Option<&str>,
    list: bool,
    errors_only: bool,
    warnings_only: bool,
    format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Determine if terminal supports color
    let use_color = atty::is(atty::Stream::Stdout);

    if list {
        // List codes
        let codes = if errors_only {
            list_error_codes()
        } else if warnings_only {
            list_warning_codes()
        } else {
            list_all_codes()
        };

        match format {
            "json" => {
                print!("[");
                for (i, code_info) in codes.iter().enumerate() {
                    if i > 0 {
                        print!(",");
                    }
                    print!("{}", code_info.format_json());
                }
                println!("]");
            }
            "markdown" => {
                println!("# ferro-hgvs Error and Warning Codes\n");
                for code_info in codes {
                    println!("{}", code_info.format_markdown());
                    println!("---\n");
                }
            }
            _ => {
                // Text format - table listing
                println!("{:<8} {:<30} Summary", "Code", "Name");
                println!("{}", "-".repeat(78));
                for code_info in codes {
                    println!(
                        "{:<8} {:<30} {}",
                        code_info.code, code_info.name, code_info.summary
                    );
                }
                println!();
                println!("Run 'ferro explain <CODE>' for detailed documentation and URLs.");
            }
        }
    } else if let Some(code_str) = code {
        // Explain a specific code
        if let Some(code_info) = get_code_info(code_str) {
            match format {
                "json" => println!("{}", code_info.format_json()),
                "markdown" => println!("{}", code_info.format_markdown()),
                _ => print!("{}", code_info.format_terminal(use_color)),
            }
        } else {
            return Err(format!(
                "Unknown code: {}. Use 'ferro explain --list' to see all codes.",
                code_str
            )
            .into());
        }
    } else {
        return Err("Please provide a code to explain or use --list to see all codes.".into());
    }

    Ok(())
}

fn output_error(input: &str, error: &FerroError, format: &str) -> io::Result<()> {
    let stderr = io::stderr();
    let mut handle = stderr.lock();
    cli_output_error(
        &mut handle,
        input,
        error,
        format.parse().unwrap_or_default(),
    )
}

/// Print a simplified capability summary for directory-based reference.
fn print_normalize_capabilities_dir(reference: Option<&PathBuf>) {
    if reference.is_some() {
        return; // Assume directory has all needed reference data
    }

    eprintln!("No reference directory provided, using mock test data.");
    eprintln!("Run 'ferro prepare' to download reference data, then use '--reference <dir>'");
    eprintln!();
}

/// Prepare reference data for normalization.
#[allow(clippy::too_many_arguments)]
fn run_prepare(
    output_dir: &Path,
    genome: &str,
    no_transcripts: bool,
    proteins: bool,
    no_cdot: bool,
    no_refseqgene: bool,
    no_lrg: bool,
    force: bool,
    patterns: Option<&Path>,
    clinvar: Option<&Path>,
    validate_canonical: Option<&Path>,
    ensembl: bool,
    derive_ng_placements: Option<PathBuf>,
    backfill_transcripts: Option<PathBuf>,
    dry_run: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::prepare::{prepare_references, PrepareConfig};

    let config = PrepareConfig {
        output_dir: output_dir.to_path_buf(),
        download_transcripts: !no_transcripts,
        download_proteins: proteins,
        download_genome: genome == "grch38" || genome == "all",
        download_genome_grch37: genome == "grch37" || genome == "all",
        download_refseqgene: !no_refseqgene,
        download_lrg: !no_lrg,
        download_cdot: !no_cdot && (genome == "grch38" || genome == "all"),
        download_cdot_grch37: !no_cdot && (genome == "grch37" || genome == "all"),
        download_ensembl: ensembl,
        skip_existing: !force,
        clinvar_file: clinvar.map(|p| p.to_path_buf()),
        patterns_file: patterns.map(|p| p.to_path_buf()),
        validate_canonical_accessions: validate_canonical.map(|p| p.to_path_buf()),
        derive_ng_placements,
        backfill_transcripts,
        genome: genome.to_string(),
        dry_run,
    };

    prepare_references(&config)?;
    Ok(())
}

/// Check reference data setup.
fn run_check(
    reference: &Path,
    build_cache: bool,
    validate_cds: bool,
    cds_allowlist: Option<&Path>,
    write_identity: bool,
    force: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::check::{
        check_reference, check_transcripts_json, print_check_summary,
        print_transcripts_json_failure, print_transcripts_json_summary,
    };
    use std::ffi::OsStr;

    // A path that doesn't exist is a user error (typically a typo) — say so, rather
    // than telling the user to run `ferro prepare` for a missing directory.
    if !reference.exists() {
        return Err(format!(
            "reference path does not exist: {} (pass a prepared reference directory or a \
             standalone transcripts.json file)",
            reference.display()
        )
        .into());
    }

    // A standalone `transcripts.json` file (the convert-gff/build-transcript output)
    // has no manifest; check it directly instead of reporting it as missing
    // reference data (#1012 comment 2). A `manifest.json` file is really the
    // prepared-reference directory case, so route it to its parent — users
    // naturally point at either. The manifest-only options (`--validate-cds`,
    // `--build-cache`, cdot cache source) do not apply to a standalone file.
    let is_manifest_file =
        reference.is_file() && reference.file_name() == Some(OsStr::new("manifest.json"));
    if reference.is_file() && !is_manifest_file {
        let ignored: Vec<&str> = [
            ("--validate-cds", validate_cds),
            ("--build-cache", build_cache),
            ("--cds-allowlist", cds_allowlist.is_some()),
            ("--write-identity", write_identity),
            ("--force", force),
        ]
        .iter()
        .filter(|(_, set)| *set)
        .map(|(name, _)| *name)
        .collect();
        if !ignored.is_empty() {
            eprintln!(
                "warning: {} do{} not apply to a standalone transcripts.json file and {} ignored",
                ignored.join(" / "),
                if ignored.len() == 1 { "es" } else { "" },
                if ignored.len() == 1 { "is" } else { "are" },
            );
        }
        return match check_transcripts_json(reference) {
            Ok(summary) => {
                print_transcripts_json_summary(&summary);
                Ok(())
            }
            Err(e) => {
                print_transcripts_json_failure(reference, &e);
                Err("Reference data check failed".into())
            }
        };
    }

    // Directory path (or a `manifest.json` file → its parent directory).
    let reference_dir: &Path = if is_manifest_file {
        // `Path::new("manifest.json").parent()` is `Some("")` (an empty path),
        // NOT `None` — so a bare relative `manifest.json` would otherwise resolve
        // the directory to "" and print an empty dir name. Map empty → ".".
        reference
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."))
    } else {
        reference
    };

    let result = check_reference(reference_dir);
    print_check_summary(&result, reference_dir);

    if !result.valid {
        return Err("Reference data check failed".into());
    }

    // Content-identity status (#1001): report drift/unstamped/verified, or
    // (with `--write-identity`) stamp the reference in place without a full
    // re-prepare.
    let manifest_path = reference_dir.join("manifest.json");
    if manifest_path.exists() {
        let v: serde_json::Value = serde_json::from_slice(&std::fs::read(&manifest_path)?)?;
        use ferro_hgvs::reference::multi_fasta::{verify_reference_identity, IdentityStatus};
        if write_identity {
            let computed = ferro_hgvs::prepare::identity::reference_identity(reference_dir, &v);
            if let Some(existing) = v.get("reference_identity").and_then(|x| x.as_str()) {
                if existing != computed && !force {
                    return Err(format!(
                        "reference already stamped {existing} but content computes to \
                         {computed} (drift). Re-run with --force to overwrite."
                    )
                    .into());
                }
            }
            // Load typed, set the directory, and save — `save()` recomputes
            // and writes the stamp from the artifacts under `reference_dir`.
            let mut m =
                ferro_hgvs::prepare::manifest::ReferenceManifest::load_or_default(reference_dir)?;
            m.reference_dir = reference_dir.to_path_buf();
            m.save()?;
            // Report exactly what `save()` persisted (it recomputes and writes its
            // own stamp), not the separately-computed drift-guard value above.
            println!(
                "Wrote reference identity: {}",
                m.reference_identity.as_deref().unwrap_or("<none>")
            );
        } else {
            match verify_reference_identity(&v, reference_dir)? {
                IdentityStatus::Verified => println!("Reference identity: verified"),
                IdentityStatus::Unstamped => println!(
                    "Reference identity: unstamped (run `--write-identity` to enable drift detection)"
                ),
                IdentityStatus::Mismatch { expected, actual } => {
                    return Err(format!(
                        "Reference identity: MISMATCH (recorded {expected}, computed {actual}) — \
                         reference drifted; re-prepare or `--write-identity` if intentional"
                    )
                    .into());
                }
            }
        }

        // Recorded cdot data-release provenance (#1001): surface which pinned
        // cdot release this reference's cdot artifacts came from. When it is
        // absent, distinguish the two very different reasons — a reference
        // legitimately prepared with no cdot at all (`--no-cdot`) is not an old
        // reference, and saying so would be a falsehood.
        for line in cdot_data_version_lines(&v) {
            println!("{line}");
        }
    }

    if validate_cds {
        run_cds_consistency_check(reference_dir, cds_allowlist)?;
    }

    if build_cache {
        // Load the reference once so the cdot bincode cache is built/refreshed
        // now (a stale cache self-heals on this load), keeping subsequent runs
        // on the fast path. The constructed provider is intentionally dropped.
        use ferro_hgvs::commands::create_reference_provider;
        eprintln!("Building/refreshing reference cache (one-time)...");
        let started = std::time::Instant::now();
        let _provider = create_reference_provider(Some(reference_dir), false)?;
        eprintln!("Reference cache ready in {:.1?}.", started.elapsed());
    }

    // Report which path the cdot cache loads from, so a silent fast-path →
    // JSON-fallback regression (the #585 class) is visible rather than only
    // showing up as a slow startup. The nightly perf gate parses this line as a
    // timing-free co-assertion that the prepared cache is on the fast path.
    // Printed after `--build-cache` so it reflects the freshly built cache.
    match ferro_hgvs::commands::cdot_cache_load_source(reference_dir) {
        Ok(Some(ferro_hgvs::data::CdotLoadSource::Archive)) => {
            println!("cdot cache: archive");
        }
        Ok(Some(ferro_hgvs::data::CdotLoadSource::JsonFallback)) => {
            println!("cdot cache: json-fallback");
        }
        Ok(None) => {} // No manifest or cdot entry — nothing to report.
        Err(e) => eprintln!("warning: could not determine cdot cache source: {}", e),
    }

    Ok(())
}

/// The `ferro check` line(s) reporting a reference's recorded cdot data-release
/// provenance (#1001). Split out of `run_check` so every branch is directly
/// testable.
///
/// Reports **per artifact**, because a reference can legitimately mix releases:
/// `ferro prepare --ensembl --no-cdot` on an existing reference refreshes only
/// the Ensembl cdot, leaving the RefSeq cdot at an older pin. A single line
/// could only name one release for both — or, if dropped, claim the reference
/// "predates tracking" when it does not.
///
/// An artifact wired but absent from the map is reported as `not recorded`,
/// without guessing why: it may have been inherited from a run whose build did
/// not track releases, or had its release discarded by a build that cleared the
/// old single-value field on a partial refresh. That is distinct from a
/// reference with no cdot at all, which has nothing to version.
fn cdot_data_version_lines(manifest: &serde_json::Value) -> Vec<String> {
    // The per-artifact map is authoritative whenever it has entries. `ferro
    // prepare` migrates the superseded scalar into it and retires the scalar,
    // so the two normally never coexist — but a hand-edited manifest could
    // carry both, and the scalar must not shadow the more precise map.
    let versions = manifest
        .get("cdot_data_versions")
        .filter(|v| v.as_object().is_some_and(|m| !m.is_empty()));

    // Legacy scalar only (a manifest written before per-artifact tracking and
    // not yet re-prepared): report exactly what it says, attributed to the
    // whole reference, rather than inventing per-artifact detail it lacks.
    if versions.is_none() {
        if let Some(ver) = manifest.get("cdot_data_version").and_then(|x| x.as_str()) {
            return vec![format!("cdot data version: {ver} (all cdot artifacts)")];
        }
    }

    let wired: Vec<&str> = ferro_hgvs::prepare::cdot_artifact_fields()
        .filter(|key| manifest.get(*key).is_some_and(|x| !x.is_null()))
        .collect();
    if wired.is_empty() {
        return vec!["cdot data version: n/a (no cdot artifacts in this reference)".to_string()];
    }

    wired
        .into_iter()
        .map(|key| {
            match versions
                .and_then(|v| v.get(key))
                .and_then(serde_json::Value::as_str)
            {
                Some(ver) => format!("cdot data version [{key}]: {ver}"),
                None => format!("cdot data version [{key}]: not recorded"),
            }
        })
        .collect()
}

/// Scan every coding transcript in the prepared reference set for CDS
/// start-codon consistency (issue #629) and print a summary.
///
/// Loads the full reference (FASTA index + cdot metadata), enumerates the
/// coding transcripts from the cdot mapper, and checks that each one's
/// `cds_start` lands on an `ATG` in the served FASTA. A non-`ATG` start is the
/// data-side signal that the cdot annotation release and the mRNA FASTA release
/// disagree — the root cause behind #625's runtime protein-prediction declines.
/// Inconsistencies are reported as warnings; this does not (yet) fail the check.
fn run_cds_consistency_check(
    reference: &Path,
    cds_allowlist: Option<&Path>,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::reference::multi_fasta::MultiFastaProvider;
    use ferro_hgvs::reference::validate::{load_cds_allowlist, scan_cds_start_codons};
    use indicatif::ProgressBar;

    let allowlist = match cds_allowlist {
        Some(path) => load_cds_allowlist(path)?,
        None => std::collections::BTreeSet::new(),
    };

    let manifest_path = reference.join("manifest.json");
    eprintln!();
    eprintln!("=== CDS start-codon consistency (#629) ===");
    let provider = MultiFastaProvider::from_manifest(&manifest_path)?;

    // Enumerate the coding transcripts from the cdot mapper (CDS bounds
    // populated). Without cdot metadata there are no CDS coordinates to check.
    let Some(cdot) = provider.cdot_mapper() else {
        eprintln!("  cdot metadata not available; nothing to validate.");
        return Ok(());
    };
    // Only check transcripts present in the FASTA at their EXACT version: a
    // cdot CDS coordinate is only meaningful against the matching-version
    // sequence. A cdot version absent from the FASTA is a version-coverage gap
    // (#645/#653), and checking it against a different version's bytes (which
    // the provider's normal version-flex resolution would do) is a false
    // positive, not a #629 CDS-frame inconsistency.
    let coding_ids: Vec<String> = cdot
        .transcript_ids()
        .filter(|id| {
            cdot.get_transcript(id)
                .is_some_and(|t| t.cds_start.is_some() && t.cds_end.is_some())
                && provider.contains_exact_sequence(id)
        })
        .map(str::to_string)
        .collect();

    eprintln!("  Scanning {} coding transcripts...", coding_ids.len());
    let pb = ProgressBar::new(coding_ids.len() as u64);
    let report = scan_cds_start_codons(
        &provider,
        coding_ids.iter().inspect(|_| pb.inc(1)),
        &allowlist,
    );
    pb.finish_and_clear();

    eprintln!("  Examined:      {}", report.coding_examined);
    eprintln!("  Inconsistent:  {}", report.inconsistent.len());
    if report.alternative_start > 0 {
        eprintln!(
            "  Alt start:     {} (recognized CTG/GTG/TTG initiators — legitimate)",
            report.alternative_start
        );
    }
    if report.allowlisted > 0 {
        eprintln!("  Allowlisted:   {}", report.allowlisted);
    }
    if report.ambiguous_skipped > 0 {
        eprintln!(
            "  Ambiguous:     {} (first codon contained N/IUPAC)",
            report.ambiguous_skipped
        );
    }
    if report.unresolved > 0 {
        eprintln!(
            "  Unresolved:    {} (no served sequence or CDS start out of range)",
            report.unresolved
        );
    }
    if report.load_errors > 0 {
        eprintln!(
            "  Load errors:   {} (transcript failed to load from the provider)",
            report.load_errors
        );
    }

    if report.has_inconsistencies() {
        // Cap the listed accessions so a badly-mismatched reference set does not
        // flood the terminal; the count above is the authoritative total.
        const MAX_LISTED: usize = 50;
        eprintln!();
        eprintln!(
            "  WARNING: {} coding transcript(s) have a non-initiator CDS start codon \
             (neither ATG nor CTG/GTG/TTG — cdot CDS coordinates inconsistent with \
             the transcript FASTA):",
            report.inconsistent.len()
        );
        for inc in report.inconsistent.iter().take(MAX_LISTED) {
            eprintln!(
                "    {} (starts with {})",
                inc.transcript_id, inc.observed_start_codon
            );
        }
        if report.inconsistent.len() > MAX_LISTED {
            eprintln!(
                "    ... and {} more",
                report.inconsistent.len() - MAX_LISTED
            );
        }
    } else {
        eprintln!("  All examined coding transcripts start with a recognized initiator.");
    }

    Ok(())
}

/// Build a transcripts.json from a FASTA file and CDS coordinates.
///
/// Creates a single-exon transcript record directly from a FASTA sequence and
/// user-supplied CDS start/end positions (in transcript coordinates), without
/// requiring a GFF3 intermediary. The exon spans the full contig. For plus-strand
/// transcripts, transcript coordinates equal genomic coordinates; for minus-strand
/// the sequence is reverse-complemented and the CDS positions are interpreted
/// against the reverse-complemented (transcript) sequence.
#[allow(clippy::too_many_arguments)]
fn run_build_transcript(
    fasta_path: &Path,
    cds_start: u64,
    cds_end: u64,
    output: &PathBuf,
    id: Option<&str>,
    strand: &str,
    contig: Option<&str>,
    gene: Option<&str>,
    genome_build: &str,
    emit_genomic_sequences: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::reference::annotation::{build_transcript, BuildTranscriptConfig};

    // Build the config, mirroring the CLI flags. The construction and
    // serialization live in `ferro_hgvs::reference::annotation::build_transcript`
    // so this command and the Python binding produce byte-identical output.
    let config = BuildTranscriptConfig {
        cds_start,
        cds_end,
        id: id.map(str::to_string),
        strand: strand.to_string(),
        contig: contig.map(str::to_string),
        gene: gene.map(str::to_string),
        genome_build: genome_build.to_string(),
        emit_genomic_sequences,
    };

    let outcome = build_transcript(fasta_path, &config)?;

    // Warn if a large genome was embedded (a no-op when nothing was emitted).
    warn_if_genomic_sequences_large(outcome.emitted_genomic_bytes);

    std::fs::write(output, serde_json::to_string_pretty(&outcome.json)?)?;
    eprintln!(
        "Wrote 1 transcript ({}) to {}",
        outcome.transcript_id,
        output.display()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_no_cdot_reference_as_nothing_to_version() {
        // A reference prepared with `--no-cdot` is not an *old* reference, so
        // `ferro check` must not claim it predates version tracking (#1001).
        let manifest = serde_json::json!({ "prepared_at": "2026-07-23" });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec!["cdot data version: n/a (no cdot artifacts in this reference)"]
        );
    }

    #[test]
    fn explicit_null_cdot_fields_do_not_count_as_wired() {
        let manifest = serde_json::json!({
            "cdot_json": serde_json::Value::Null,
            "ensembl_cdot_grch37_json": serde_json::Value::Null,
        });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec!["cdot data version: n/a (no cdot artifacts in this reference)"]
        );
    }

    #[test]
    fn reports_a_release_per_wired_artifact() {
        let manifest = serde_json::json!({
            "cdot_json": "cdot/cdot-0.2.32.refseq.GRCh38.json.gz",
            "ensembl_cdot_json": "cdot/cdot-0.2.32.ensembl.GRCh38.json.gz",
            "cdot_data_versions": {
                "cdot_json": "data_v0.2.32",
                "ensembl_cdot_json": "data_v0.2.32",
            },
        });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec![
                "cdot data version [cdot_json]: data_v0.2.32",
                "cdot data version [ensembl_cdot_json]: data_v0.2.32",
            ]
        );
    }

    #[test]
    fn reports_a_mixed_reference_honestly() {
        // The case that motivated per-artifact provenance: `--ensembl
        // --no-cdot` on an existing reference refreshes only the Ensembl cdot
        // (plain `--ensembl` also refreshes GRCh38 RefSeq). Each artifact is
        // reported with the release it actually came from — neither a single
        // wrong release, nor a bogus "predates tracking".
        let manifest = serde_json::json!({
            "cdot_json": "cdot/cdot-0.2.28.refseq.GRCh38.json.gz",
            "ensembl_cdot_json": "cdot/cdot-0.2.32.ensembl.GRCh38.json.gz",
            "cdot_data_versions": {
                "cdot_json": "data_v0.2.28",
                "ensembl_cdot_json": "data_v0.2.32",
            },
        });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec![
                "cdot data version [cdot_json]: data_v0.2.28",
                "cdot data version [ensembl_cdot_json]: data_v0.2.32",
            ]
        );
    }

    #[test]
    fn reports_a_wired_but_unrecorded_artifact_as_not_recorded() {
        let manifest = serde_json::json!({ "cdot_json": "cdot/cdot-0.2.28.refseq.GRCh38.json.gz" });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec!["cdot data version [cdot_json]: not recorded"]
        );
    }

    #[test]
    fn a_non_empty_map_wins_over_a_stale_legacy_scalar() {
        // The exact state the bug produced: a schema-3 reference re-prepared by
        // a build that stamped the map but left the scalar in place. Reporting
        // the scalar first suppressed every per-artifact line and named a
        // release older than the reference actually had. `ferro prepare` now
        // migrates and retires the scalar, but a hand-edited manifest can still
        // carry both, and the map must win.
        let manifest = serde_json::json!({
            "cdot_json": "cdot/cdot-0.2.28.refseq.GRCh38.json.gz",
            "ensembl_cdot_json": "cdot/cdot-0.2.32.ensembl.GRCh38.json.gz",
            "cdot_data_version": "data_v0.2.28",
            "cdot_data_versions": {
                "cdot_json": "data_v0.2.28",
                "ensembl_cdot_json": "data_v0.2.32",
            },
        });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec![
                "cdot data version [cdot_json]: data_v0.2.28",
                "cdot data version [ensembl_cdot_json]: data_v0.2.32",
            ]
        );
    }

    #[test]
    fn an_empty_map_does_not_suppress_the_legacy_scalar() {
        // `skip_serializing_if` keeps an empty map out of the JSON, but a
        // hand-written `{}` must not shadow the scalar either.
        let manifest = serde_json::json!({
            "cdot_json": "cdot/cdot-0.2.32.refseq.GRCh38.json.gz",
            "cdot_data_version": "data_v0.2.32",
            "cdot_data_versions": {},
        });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec!["cdot data version: data_v0.2.32 (all cdot artifacts)"]
        );
    }

    #[test]
    fn legacy_scalar_version_is_still_reported() {
        // A manifest written before per-artifact tracking carries the scalar.
        // Report what it says, scoped to the whole reference — it carries no
        // per-artifact detail to invent.
        let manifest = serde_json::json!({
            "cdot_json": "cdot/cdot-0.2.32.refseq.GRCh38.json.gz",
            "cdot_data_version": "data_v0.2.32",
        });
        assert_eq!(
            cdot_data_version_lines(&manifest),
            vec!["cdot data version: data_v0.2.32 (all cdot artifacts)"]
        );
    }

    #[test]
    fn every_cdot_field_is_reported_when_wired() {
        for key in ferro_hgvs::prepare::cdot_artifact_fields() {
            let manifest = serde_json::json!({ key: "cdot/some-cdot.json.gz" });
            let lines = cdot_data_version_lines(&manifest);
            assert_eq!(lines.len(), 1, "{key} should be reported");
            assert!(
                lines[0].contains(key),
                "{key} should name itself in its line, got: {}",
                lines[0]
            );
        }
    }
}
