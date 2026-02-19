// Copyright (c) 2024-2025 Fulcrum Genomics LLC
// SPDX-License-Identifier: MIT

//! ferro CLI
//!
//! Command-line interface for HGVS variant normalization and VCF annotation.

use clap::{Parser, Subcommand};
use ferro_hgvs::cli::{
    output_error as cli_output_error, output_error_with_context as cli_output_error_with_context,
    parse_genome_build, parse_shuffle_direction, parse_vcf_line, process_input_line,
    reverse_complement,
};
use ferro_hgvs::config::FerroConfig;
use ferro_hgvs::error_handling::{
    get_code_info, list_all_codes, list_error_codes, list_warning_codes, ErrorConfig, ErrorMode,
    ErrorType,
};
use ferro_hgvs::reference::transcript::GenomeBuild;
use ferro_hgvs::reference::FastaProvider;
use ferro_hgvs::reference::TranscriptDb;
use ferro_hgvs::vcf::{generate_info_header_lines, open_vcf, VcfAnnotator, VcfRecord};
use ferro_hgvs::{parse_hgvs, FerroError, NormalizeConfig, Normalizer, ReferenceProvider};
use flate2::read::MultiGzDecoder;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

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

        /// Output format
        #[arg(short = 'f', long, default_value = "text", value_parser = ["text", "json"])]
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

        /// Warning codes to always reject (comma-separated, e.g., W4002)
        #[arg(long, value_delimiter = ',')]
        reject: Vec<String>,

        /// Output timing information to JSON file
        #[arg(short = 't', long)]
        timing: Option<PathBuf>,

        /// Number of parallel workers (default: 1)
        #[arg(short = 'j', long, default_value = "1")]
        workers: usize,
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

        /// Warning codes to always reject (comma-separated, e.g., W4002)
        #[arg(long, value_delimiter = ',')]
        reject: Vec<String>,

        /// Output timing information to JSON file
        #[arg(short = 't', long)]
        timing: Option<PathBuf>,
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

        /// Dry run - show what would be downloaded without downloading
        #[arg(long)]
        dry_run: bool,
    },

    /// Check reference data setup
    Check {
        /// Reference directory to check
        #[arg(long, default_value = "ferro-reference")]
        reference: PathBuf,
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
        } => run_annotate_vcf(
            &input,
            output.as_ref(),
            gff.as_ref(),
            &build,
            ann_format,
            all_transcripts,
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
            )
        }
        Commands::Parse {
            variant,
            input,
            output,
            format,
            error_mode,
            ignore,
            reject,
            timing,
        } => {
            let config = build_error_config(&error_mode, &ignore, &reject);
            run_parse(
                variant.as_deref(),
                input.as_ref(),
                output.as_ref(),
                &format,
                timing.as_ref(),
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
        } => run_convert_gff(
            &gff,
            fasta.as_ref(),
            output.as_ref(),
            &build,
            mane_only,
            transcripts.as_deref(),
            genes.as_deref(),
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
            no_refseqgene,
            no_lrg,
            force,
            patterns,
            clinvar,
            dry_run,
        } => run_prepare(
            &output_dir,
            &genome,
            no_refseqgene,
            no_lrg,
            force,
            patterns.as_deref(),
            clinvar.as_deref(),
            dry_run,
        ),
        Commands::Check { reference } => run_check(&reference),
    }
}

fn run_annotate_vcf(
    input: &PathBuf,
    output: Option<&PathBuf>,
    gff: Option<&PathBuf>,
    build: &str,
    ann_format: bool,
    all_transcripts: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Parse genome build
    let genome_build = parse_genome_build(build);

    // Load transcript database
    let db = if let Some(gff_path) = gff {
        let ext = gff_path.extension().and_then(|e| e.to_str()).unwrap_or("");
        if ext == "gtf" || ext == "gtf.gz" {
            ferro_hgvs::reference::loader::load_gtf(gff_path)?
        } else {
            ferro_hgvs::reference::loader::load_gff3(gff_path)?
        }
    } else {
        // Create empty database
        TranscriptDb::with_build(genome_build)
    };

    // Create annotator
    let annotator = VcfAnnotator::new(&db)
        .use_ann_field(ann_format)
        .include_all_transcripts(all_transcripts);

    // Set up output
    let mut writer: Box<dyn Write> = if let Some(out_path) = output {
        if out_path.to_string_lossy() == "-" {
            Box::new(io::stdout())
        } else {
            Box::new(std::fs::File::create(out_path)?)
        }
    } else {
        Box::new(io::stdout())
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
                // Parse and annotate
                if let Ok(record) = parse_vcf_line(&line) {
                    let annotated = annotator.annotate(&record)?;
                    writeln!(writer, "{}", annotated)?;
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

        // Process records
        for record in reader.records() {
            let record = record?;
            let annotated = annotator.annotate(&record)?;
            writeln!(writer, "{}", annotated)?;
        }
    }

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
        let ext = gff_path.extension().and_then(|e| e.to_str()).unwrap_or("");
        if ext == "gtf" || ext == "gtf.gz" {
            ferro_hgvs::reference::loader::load_gtf(gff_path)?
        } else {
            ferro_hgvs::reference::loader::load_gff3(gff_path)?
        }
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

#[allow(clippy::too_many_arguments)]
fn run_normalize(
    variant: Option<&str>,
    input: Option<&PathBuf>,
    output: Option<&PathBuf>,
    format: &str,
    direction: &str,
    reference: Option<&PathBuf>,
    timing: Option<&PathBuf>,
    _workers: usize, // Reserved for future parallel implementation
    error_config: &ErrorConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::commands::create_reference_provider;
    use std::time::Instant;

    let preprocessor = error_config.preprocessor();
    let config = NormalizeConfig::default().with_direction(parse_shuffle_direction(direction));

    // Create reference provider from directory
    let provider = create_reference_provider(reference.map(|p| p.as_path()))?;
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

    let process = |v: &str, writer: &mut dyn Write| -> Result<(), FerroError> {
        // Preprocess input according to error mode
        let preprocess_result = preprocessor.preprocess(v);
        if !preprocess_result.success {
            // Preprocessing failed (strict mode rejection)
            return Err(FerroError::parse(
                0,
                format!(
                    "Input contains errors that are rejected in strict mode: {}",
                    preprocess_result
                        .warnings
                        .iter()
                        .map(|w| w.message.as_str())
                        .collect::<Vec<_>>()
                        .join("; ")
                ),
            ));
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
                    writer,
                    r#"{{"input": "{}", "output": "{}", "status": "ok", "corrections": [{}]}}"#,
                    v,
                    normalized,
                    corrections.join(",")
                )
                .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }
            _ => {
                // Print warnings to stderr if any
                for warning in &preprocess_result.warnings {
                    eprintln!(
                        "warning[{}]: {}",
                        warning.error_type.code(),
                        warning.message
                    );
                }
                if v == normalized.to_string() {
                    writeln!(writer, "{}", normalized)
                        .map_err(|e| FerroError::Io { msg: e.to_string() })?;
                } else {
                    writeln!(writer, "{} -> {}", v, normalized)
                        .map_err(|e| FerroError::Io { msg: e.to_string() })?;
                }
            }
        }
        Ok(())
    };

    if let Some(v) = variant {
        total_count += 1;
        if let Err(e) = process(v, &mut writer) {
            output_error(v, &e, format)?;
            error_count += 1;
        } else {
            success_count += 1;
        }
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            let is_first = line_num == 0;
            if let Some(variant_str) = process_input_line(&line, is_first) {
                total_count += 1;
                if let Err(e) = process(variant_str, &mut writer) {
                    output_error_with_line(variant_str, &e, format, Some(line_num + 1))?;
                    error_count += 1;
                } else {
                    success_count += 1;
                }
            }
        }
    } else {
        let stdin = io::stdin();
        for (line_num, line) in stdin.lock().lines().enumerate() {
            let line = line?;
            let is_first = line_num == 0;
            if let Some(variant_str) = process_input_line(&line, is_first) {
                total_count += 1;
                if let Err(e) = process(variant_str, &mut writer) {
                    output_error_with_line(variant_str, &e, format, Some(line_num + 1))?;
                    error_count += 1;
                } else {
                    success_count += 1;
                }
            }
        }
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

fn run_parse(
    variant: Option<&str>,
    input: Option<&PathBuf>,
    output: Option<&PathBuf>,
    format: &str,
    timing: Option<&PathBuf>,
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

    let process = |v: &str, writer: &mut dyn Write| -> Result<(), FerroError> {
        // Preprocess input according to error mode
        let preprocess_result = preprocessor.preprocess(v);
        if !preprocess_result.success {
            // Preprocessing failed (strict mode rejection)
            return Err(FerroError::parse(
                0,
                format!(
                    "Input contains errors that are rejected in strict mode: {}",
                    preprocess_result
                        .warnings
                        .iter()
                        .map(|w| w.message.as_str())
                        .collect::<Vec<_>>()
                        .join("; ")
                ),
            ));
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
                    writer,
                    r#"{{"input": "{}", "parsed": "{}", "status": "ok", "corrections": [{}]}}"#,
                    v,
                    parsed,
                    corrections.join(",")
                )
                .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }
            _ => {
                // Print warnings to stderr if any
                for warning in &preprocess_result.warnings {
                    eprintln!(
                        "warning[{}]: {}",
                        warning.error_type.code(),
                        warning.message
                    );
                }
                writeln!(writer, "{} -> {}", v, parsed)
                    .map_err(|e| FerroError::Io { msg: e.to_string() })?;
            }
        }
        Ok(())
    };

    if let Some(v) = variant {
        total_count += 1;
        if let Err(e) = process(v, &mut writer) {
            output_error(v, &e, format)?;
            error_count += 1;
        } else {
            success_count += 1;
        }
    } else if let Some(input_path) = input {
        let file = std::fs::File::open(input_path)?;
        let reader = io::BufReader::new(file);
        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            let is_first = line_num == 0;
            if let Some(variant_str) = process_input_line(&line, is_first) {
                total_count += 1;
                if let Err(e) = process(variant_str, &mut writer) {
                    output_error_with_line(variant_str, &e, format, Some(line_num + 1))?;
                    error_count += 1;
                } else {
                    success_count += 1;
                }
            }
        }
    } else {
        let stdin = io::stdin();
        for (line_num, line) in stdin.lock().lines().enumerate() {
            let line = line?;
            let is_first = line_num == 0;
            if let Some(variant_str) = process_input_line(&line, is_first) {
                total_count += 1;
                if let Err(e) = process(variant_str, &mut writer) {
                    output_error_with_line(variant_str, &e, format, Some(line_num + 1))?;
                    error_count += 1;
                } else {
                    success_count += 1;
                }
            }
        }
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

fn run_convert_gff(
    gff: &PathBuf,
    fasta: Option<&PathBuf>,
    output: Option<&PathBuf>,
    build: &str,
    mane_only: bool,
    transcripts: Option<&str>,
    genes: Option<&str>,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::reference::transcript::ManeStatus;
    use std::collections::HashSet;

    // Parse genome build
    let genome_build = parse_genome_build(build);

    // Load GFF/GTF file
    let ext = gff.extension().and_then(|e| e.to_str()).unwrap_or("");
    let mut db = if ext == "gtf" || gff.to_string_lossy().contains(".gtf") {
        ferro_hgvs::reference::loader::load_gtf(gff)?
    } else {
        ferro_hgvs::reference::loader::load_gff3(gff)?
    };
    db.genome_build = genome_build;

    // Load FASTA if provided and extract sequences
    let fasta_provider = if let Some(fasta_path) = fasta {
        Some(FastaProvider::new(fasta_path)?)
    } else {
        None
    };

    // Parse transcript and gene filters
    let transcript_filter: Option<HashSet<&str>> =
        transcripts.map(|s| s.split(',').map(|t| t.trim()).collect());
    let gene_filter: Option<HashSet<&str>> =
        genes.map(|s| s.split(',').map(|g| g.trim()).collect());

    // Collect transcripts to output
    let mut output_transcripts: Vec<serde_json::Value> = Vec::new();

    for (id, tx) in db.iter() {
        // Apply filters
        if mane_only && tx.mane_status == ManeStatus::None {
            continue;
        }

        if let Some(ref filter) = transcript_filter {
            if !filter.contains(id.as_str()) {
                continue;
            }
        }

        if let Some(ref filter) = gene_filter {
            if let Some(ref gene) = tx.gene_symbol {
                if !filter.contains(gene.as_str()) {
                    continue;
                }
            } else {
                continue;
            }
        }

        // Extract sequence from FASTA if available
        let sequence = if let (Some(ref fasta), Some(ref chrom), Some(_start), Some(_end)) = (
            &fasta_provider,
            &tx.chromosome,
            tx.genomic_start,
            tx.genomic_end,
        ) {
            // Extract exonic sequences
            let mut exon_seqs = Vec::new();
            for exon in &tx.exons {
                if let (Some(g_start), Some(g_end)) = (exon.genomic_start, exon.genomic_end) {
                    // FASTA uses 0-based coordinates, GFF uses 1-based
                    match fasta.get_sequence(chrom, g_start - 1, g_end) {
                        Ok(seq) => exon_seqs.push(seq),
                        Err(_) => exon_seqs.push("N".repeat((g_end - g_start + 1) as usize)),
                    }
                }
            }

            // Join exons and handle strand
            let mut full_seq = exon_seqs.join("");
            if tx.strand == ferro_hgvs::reference::transcript::Strand::Minus {
                full_seq = reverse_complement(&full_seq);
            }
            full_seq
        } else if !tx.sequence.is_empty() && !tx.sequence.chars().all(|c| c == 'N') {
            tx.sequence.clone()
        } else {
            // No FASTA, use placeholder
            "N".repeat(tx.exons.iter().map(|e| e.end - e.start + 1).sum::<u64>() as usize)
        };

        // Build exon objects
        let exons: Vec<serde_json::Value> = tx
            .exons
            .iter()
            .map(|e| {
                let mut exon_obj = serde_json::json!({
                    "number": e.number,
                    "start": e.start,
                    "end": e.end,
                });
                if let Some(g_start) = e.genomic_start {
                    exon_obj["genomic_start"] = serde_json::json!(g_start);
                }
                if let Some(g_end) = e.genomic_end {
                    exon_obj["genomic_end"] = serde_json::json!(g_end);
                }
                exon_obj
            })
            .collect();

        // Build transcript object
        let mut tx_obj = serde_json::json!({
            "id": id,
            "sequence": sequence,
            "exons": exons,
        });

        // Add optional fields
        if let Some(ref gene) = tx.gene_symbol {
            tx_obj["gene_symbol"] = serde_json::json!(gene);
        }
        if let Some(ref chrom) = tx.chromosome {
            tx_obj["chromosome"] = serde_json::json!(chrom);
        }
        if let Some(start) = tx.genomic_start {
            tx_obj["genomic_start"] = serde_json::json!(start);
        }
        if let Some(end) = tx.genomic_end {
            tx_obj["genomic_end"] = serde_json::json!(end);
        }
        if let Some(cds_start) = tx.cds_start {
            tx_obj["cds_start"] = serde_json::json!(cds_start);
        }
        if let Some(cds_end) = tx.cds_end {
            tx_obj["cds_end"] = serde_json::json!(cds_end);
        }

        tx_obj["strand"] = serde_json::json!(match tx.strand {
            ferro_hgvs::reference::transcript::Strand::Plus => "+",
            ferro_hgvs::reference::transcript::Strand::Minus => "-",
        });

        tx_obj["genome_build"] = serde_json::json!(match genome_build {
            GenomeBuild::GRCh37 => "GRCh37",
            GenomeBuild::GRCh38 => "GRCh38",
            GenomeBuild::Unknown => "unknown",
        });

        if tx.mane_status != ManeStatus::None {
            tx_obj["mane_status"] = serde_json::json!(match tx.mane_status {
                ManeStatus::Select => "MANE_Select",
                ManeStatus::PlusClinical => "MANE_Plus_Clinical",
                ManeStatus::None => "none",
            });
        }

        output_transcripts.push(tx_obj);
    }

    // Create output JSON
    let output_json = serde_json::json!({
        "version": "1.0",
        "genome_build": match genome_build {
            GenomeBuild::GRCh37 => "GRCh37",
            GenomeBuild::GRCh38 => "GRCh38",
            GenomeBuild::Unknown => "unknown",
        },
        "transcripts": output_transcripts,
    });

    // Write output
    let mut writer: Box<dyn Write> = if let Some(out_path) = output {
        Box::new(std::fs::File::create(out_path)?)
    } else {
        Box::new(io::stdout())
    };

    writeln!(writer, "{}", serde_json::to_string_pretty(&output_json)?)?;

    eprintln!(
        "Converted {} transcripts from {}",
        output_transcripts.len(),
        gff.display()
    );

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
                            ferro_hgvs::hgvs::edit::ProteinEdit::Insertion { .. } => {
                                predictor.classify_indel(0, 3)
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
        if let Some(error_type) = code_to_error_type(code) {
            config.set_override(error_type, ErrorOverride::SilentCorrect);
        }
    }

    // Apply reject overrides
    for code in reject {
        if let Some(error_type) = code_to_error_type(code) {
            config.set_override(error_type, ErrorOverride::Reject);
        }
    }

    config
}

/// Map a warning code string to an ErrorType.
fn code_to_error_type(code: &str) -> Option<ErrorType> {
    match code.to_uppercase().as_str() {
        "W1001" => Some(ErrorType::LowercaseAminoAcid),
        "W1002" => Some(ErrorType::SingleLetterAminoAcid),
        "W1003" => Some(ErrorType::LowercaseAccessionPrefix),
        "W1004" => Some(ErrorType::MixedCaseEditType),
        "W2001" => Some(ErrorType::WrongDashCharacter),
        "W2002" => Some(ErrorType::WrongQuoteCharacter),
        "W2003" => Some(ErrorType::ExtraWhitespace),
        "W2004" => Some(ErrorType::InvalidUnicodeCharacter),
        "W3001" => Some(ErrorType::MissingVersion),
        "W3002" => Some(ErrorType::ProteinSubstitutionArrow),
        "W3003" => Some(ErrorType::OldSubstitutionSyntax),
        "W3004" => Some(ErrorType::OldAlleleFormat),
        "W3005" => Some(ErrorType::TrailingAnnotation),
        "W3006" => Some(ErrorType::MissingCoordinatePrefix),
        "W4001" => Some(ErrorType::SwappedPositions),
        "W4002" => Some(ErrorType::PositionZero),
        "W5001" => Some(ErrorType::RefSeqMismatch),
        _ => None,
    }
}

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

fn output_error_with_line(
    input: &str,
    error: &FerroError,
    format: &str,
    line_number: Option<usize>,
) -> io::Result<()> {
    let stderr = io::stderr();
    let mut handle = stderr.lock();
    cli_output_error_with_context(
        &mut handle,
        input,
        error,
        format.parse().unwrap_or_default(),
        line_number,
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
    no_refseqgene: bool,
    no_lrg: bool,
    force: bool,
    patterns: Option<&Path>,
    clinvar: Option<&Path>,
    dry_run: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::prepare::{prepare_references, PrepareConfig};

    let config = PrepareConfig {
        output_dir: output_dir.to_path_buf(),
        download_transcripts: true,
        download_genome: genome == "grch38" || genome == "all",
        download_genome_grch37: genome == "grch37" || genome == "all",
        download_refseqgene: !no_refseqgene,
        download_lrg: !no_lrg,
        download_cdot: true,
        skip_existing: !force,
        clinvar_file: clinvar.map(|p| p.to_path_buf()),
        patterns_file: patterns.map(|p| p.to_path_buf()),
        dry_run,
    };

    prepare_references(&config)?;
    Ok(())
}

/// Check reference data setup.
fn run_check(reference: &Path) -> Result<(), Box<dyn std::error::Error>> {
    use ferro_hgvs::check::{check_reference, print_check_summary};

    let result = check_reference(reference);
    print_check_summary(&result, reference);

    if result.valid {
        Ok(())
    } else {
        Err("Reference data check failed".into())
    }
}
