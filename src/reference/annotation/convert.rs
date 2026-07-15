//! GFF3/GTF → `transcripts.json` conversion.
//!
//! [`convert_gff`] is the library entry point behind the `ferro convert-gff`
//! CLI command and the `ferro_hgvs.convert_gff` Python binding. It loads an
//! annotation file via [`load_annotations`], applies the transcript/gene/MANE
//! filters, extracts (or synthesizes) each transcript's spliced sequence, and
//! serializes the result into the `transcripts.json` object shape consumed by
//! [`JsonProvider`] and `Normalizer(reference_json=...)`.
//!
//! Both the CLI (`run_convert_gff`) and the Python binding call this one
//! function so their output is byte-identical for the same inputs/flags — the
//! only difference is I/O framing (where the JSON is written) and the stderr
//! warnings the CLI emits from the returned [`ConvertGffOutcome`] fields.
//!
//! [`JsonProvider`]: crate::reference::mock::JsonProvider

use std::collections::{BTreeMap, HashSet};
use std::path::Path;

use serde_json::{json, Map, Value};

use crate::error::FerroError;
use crate::reference::annotation::{load_annotations, LoaderConfig, LoaderReport};
use crate::reference::fasta::FastaProvider;
use crate::reference::provider::ReferenceProvider;
use crate::reference::transcript::{GenomeBuild, ManeStatus, Strand};
use crate::sequence::reverse_complement;

/// Threshold above which embedding genomic sequence into a `transcripts.json`
/// is worth warning about: a `genomic_sequences` map this large makes the file
/// slow to store and parse. Shared by the CLI's stderr warning and the Python
/// binding so both use the same cutoff.
pub const LARGE_GENOMIC_SEQUENCES_WARN_BYTES: u64 = 50_000_000; // 50 MB

/// Configuration for [`convert_gff`], mirroring the `ferro convert-gff` CLI
/// flags one-to-one.
///
/// Construct with [`ConvertGffConfig::new`] and set fields as needed, or build
/// one directly with a struct literal; all fields are public so callers (the
/// CLI, the Python binding) can populate it however is convenient.
#[derive(Debug, Clone)]
pub struct ConvertGffConfig {
    /// Genome build recorded in the output and used when constructing transcripts.
    pub genome_build: GenomeBuild,
    /// Only emit MANE Select / MANE Plus Clinical transcripts.
    pub mane_only: bool,
    /// Restrict output to these transcript IDs (already split; `None` = no filter).
    pub transcripts: Option<Vec<String>>,
    /// Restrict output to these gene symbols (already split; `None` = no filter).
    pub genes: Option<Vec<String>>,
    /// Promote loader warnings to errors and fail if any error diagnostic is recorded.
    pub strict: bool,
    /// Suppress loader diagnostic warnings (still counted in the report).
    pub silent: bool,
    /// Skip CDS-length / start-codon FASTA-aware validation even when a FASTA is supplied.
    pub no_validate_fasta: bool,
    /// Embed the referenced contig sequences under `genomic_sequences` so the
    /// output is genome-capable. Requires a FASTA (`fasta` argument of
    /// [`convert_gff`]); see #1026.
    pub emit_genomic_sequences: bool,
}

impl ConvertGffConfig {
    /// A default configuration: GRCh38, no filters, lenient errors, FASTA
    /// validation on, no genomic-sequence embedding.
    pub fn new() -> Self {
        Self {
            genome_build: GenomeBuild::GRCh38,
            mane_only: false,
            transcripts: None,
            genes: None,
            strict: false,
            silent: false,
            no_validate_fasta: false,
            emit_genomic_sequences: false,
        }
    }
}

impl Default for ConvertGffConfig {
    fn default() -> Self {
        Self::new()
    }
}

/// Result of a [`convert_gff`] run.
///
/// The serialized `transcripts.json` lives in [`json`](Self::json); the other
/// fields carry the information the CLI needs to reproduce its stderr side
/// effects (a loader summary line and the `--emit-genomic-sequences` warnings)
/// without the library itself writing to stderr.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ConvertGffOutcome {
    /// The `transcripts.json` document (a JSON object with `version`,
    /// `genome_build`, `transcripts`, and — when emitted — `genomic_sequences`).
    pub json: Value,
    /// Diagnostics aggregated while loading the annotation file.
    pub report: LoaderReport,
    /// Total bytes of genomic sequence embedded under `genomic_sequences`
    /// (`0` when `emit_genomic_sequences` was not requested or nothing was
    /// emitted). Callers may warn when this is large.
    pub emitted_genomic_bytes: u64,
    /// `true` when `emit_genomic_sequences` was requested but no emitted
    /// transcript carried a genomic placement, so `genomic_sequences` was not
    /// written and the reference stays transcript-only. Callers may warn.
    pub emit_requested_no_placement: bool,
}

/// Convert a GFF3/GTF annotation file into the `transcripts.json` object shape.
///
/// # Arguments
///
/// * `gff` - Path to the GFF3 or GTF annotation file.
/// * `fasta` - Optional reference FASTA used to extract exonic sequences and,
///   when [`ConvertGffConfig::emit_genomic_sequences`] is set, the embedded
///   contig sequences. Also drives CDS-length / start-codon validation unless
///   [`ConvertGffConfig::no_validate_fasta`] is set.
/// * `config` - Filters, error mode, and genome build; see [`ConvertGffConfig`].
///
/// # Errors
///
/// Returns [`FerroError`] if the annotation file cannot be read/parsed (or, in
/// strict mode, records an error diagnostic), if a supplied FASTA cannot be
/// opened, or if `emit_genomic_sequences` is set but a required contig is
/// absent from — or too short in — the FASTA to cover a transcript placement.
pub fn convert_gff(
    gff: &Path,
    fasta: Option<&Path>,
    config: &ConvertGffConfig,
) -> Result<ConvertGffOutcome, FerroError> {
    let genome_build = config.genome_build;

    // Build loader config with error mode, genome build, and optional FASTA validation.
    let mut cfg = LoaderConfig::new().with_genome_build(genome_build);
    if config.strict {
        cfg = cfg.with_strict();
    } else if config.silent {
        cfg = cfg.with_silent();
    }
    // Pass a FASTA provider to load_annotations for CDS-length / start-codon validation.
    // A second FastaProvider is constructed below for sequence extraction; FastaProvider::new
    // is cheap (reads only the index), so the double construction is acceptable.
    if let Some(fasta_path) = fasta {
        cfg = cfg.with_fasta(FastaProvider::new(fasta_path)?);
    }
    if config.no_validate_fasta {
        cfg = cfg.with_no_validate_fasta();
    }

    // Load GFF/GTF file via load_annotations.
    let (mut db, report) = load_annotations(gff, &cfg)?;
    db.genome_build = genome_build;

    // Load FASTA if provided and extract sequences.
    let fasta_provider = if let Some(fasta_path) = fasta {
        Some(FastaProvider::new(fasta_path)?)
    } else {
        None
    };

    // Parse transcript and gene filters.
    let transcript_filter: Option<HashSet<&str>> = config
        .transcripts
        .as_ref()
        .map(|ids| ids.iter().map(|t| t.trim()).collect());
    let gene_filter: Option<HashSet<&str>> = config
        .genes
        .as_ref()
        .map(|gs| gs.iter().map(|g| g.trim()).collect());

    // Collect transcripts to output.
    let mut output_transcripts: Vec<Value> = Vec::new();

    // For each contig an emitted transcript is placed on, the highest 1-based
    // genomic coordinate its placement references (max over exon `genomic_end` and
    // the transcript's own `genomic_end`). Used to emit `genomic_sequences` when
    // `emit_genomic_sequences` is set, and — crucially — to fail fast at emit time
    // if the FASTA contig cannot cover the placement, using the SAME invariant
    // `JsonProvider::from_json` enforces at load, so convert-gff can never emit a
    // file its own loader would reject (#1026). BTreeMap keeps the output
    // deterministic.
    let mut contig_required_len: BTreeMap<String, u64> = BTreeMap::new();

    for (id, tx) in db.iter() {
        // Apply filters.
        if config.mane_only && tx.mane_status == ManeStatus::None {
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

        // Extract sequence from FASTA if available.
        let sequence = if let (Some(ref fasta), Some(ref chrom), Some(_start), Some(_end)) = (
            &fasta_provider,
            &tx.chromosome,
            tx.genomic_start,
            tx.genomic_end,
        ) {
            // Extract exonic sequences.
            let mut exon_seqs = Vec::new();
            for exon in &tx.exons {
                if let (Some(g_start), Some(g_end)) = (exon.genomic_start, exon.genomic_end) {
                    // FASTA uses 0-based coordinates, GFF uses 1-based.
                    match fasta.get_sequence(chrom, g_start - 1, g_end) {
                        Ok(seq) => exon_seqs.push(seq),
                        Err(_) => exon_seqs.push("N".repeat((g_end - g_start + 1) as usize)),
                    }
                }
            }

            // Join exons and handle strand.
            let mut full_seq = exon_seqs.join("");
            if tx.strand == Strand::Minus {
                full_seq = reverse_complement(&full_seq);
            }
            full_seq
        } else if let Some(seq) = tx
            .sequence
            .as_deref()
            .filter(|s| !s.is_empty() && !s.chars().all(|c| c == 'N'))
        {
            seq.to_string()
        } else {
            // No FASTA, use placeholder.
            "N".repeat(tx.exons.iter().map(|e| e.end - e.start + 1).sum::<u64>() as usize)
        };

        // Build exon objects.
        let exons: Vec<Value> = tx
            .exons
            .iter()
            .map(|e| {
                let mut exon_obj = json!({
                    "number": e.number,
                    "start": e.start,
                    "end": e.end,
                });
                if let Some(g_start) = e.genomic_start {
                    exon_obj["genomic_start"] = json!(g_start);
                }
                if let Some(g_end) = e.genomic_end {
                    exon_obj["genomic_end"] = json!(g_end);
                }
                exon_obj
            })
            .collect();

        // Build transcript object.
        let mut tx_obj = json!({
            "id": id,
            "sequence": sequence,
            "exons": exons,
        });

        // Add optional fields.
        if let Some(ref gene) = tx.gene_symbol {
            tx_obj["gene_symbol"] = json!(gene);
        }
        if let Some(ref chrom) = tx.chromosome {
            tx_obj["chromosome"] = json!(chrom);
            // Track the highest genomic coordinate this transcript's placement
            // needs on its contig, mirroring the loader's `required_len`.
            let required = tx
                .exons
                .iter()
                .filter_map(|e| e.genomic_end)
                .chain(tx.genomic_end)
                .max();
            if let Some(required) = required {
                let entry = contig_required_len.entry(chrom.clone()).or_insert(0);
                *entry = (*entry).max(required);
            }
        }
        if let Some(start) = tx.genomic_start {
            tx_obj["genomic_start"] = json!(start);
        }
        if let Some(end) = tx.genomic_end {
            tx_obj["genomic_end"] = json!(end);
        }
        if let Some(cds_start) = tx.cds_start {
            tx_obj["cds_start"] = json!(cds_start);
        }
        if let Some(cds_end) = tx.cds_end {
            tx_obj["cds_end"] = json!(cds_end);
        }

        tx_obj["strand"] = json!(match tx.strand {
            Strand::Plus => "+",
            Strand::Minus => "-",
            // `Strand` is `#[non_exhaustive]`, but all variants are matched
            // within the defining crate; `Unknown` and any future variant map
            // to the GFF3 unspecified-strand marker.
            _ => ".",
        });

        tx_obj["genome_build"] = json!(match genome_build {
            GenomeBuild::GRCh37 => "GRCh37",
            GenomeBuild::GRCh38 => "GRCh38",
            GenomeBuild::Unknown => "unknown",
        });

        if tx.mane_status != ManeStatus::None {
            tx_obj["mane_status"] = json!(match tx.mane_status {
                ManeStatus::Select => "MANE_Select",
                ManeStatus::PlusClinical => "MANE_Plus_Clinical",
                ManeStatus::None => "none",
            });
        }

        output_transcripts.push(tx_obj);
    }

    // Create output JSON.
    let mut output_json = json!({
        "version": "1.0",
        "genome_build": match genome_build {
            GenomeBuild::GRCh37 => "GRCh37",
            GenomeBuild::GRCh38 => "GRCh38",
            GenomeBuild::Unknown => "unknown",
        },
        "transcripts": output_transcripts,
    });

    // Optionally emit the referenced contig sequences so the transcripts.json is
    // genome-capable. The exon `genomic_start`/`genomic_end` coordinates are
    // absolute contig positions, so we emit each referenced contig's full forward
    // sequence keyed by chromosome name — this is exactly what `JsonProvider`
    // slices by genomic coordinate when running the genome-aware rules (#1026).
    let mut emitted_genomic_bytes: u64 = 0;
    let mut emit_requested_no_placement = false;
    if config.emit_genomic_sequences {
        let fasta = fasta_provider.as_ref().ok_or_else(|| FerroError::Io {
            msg: "--emit-genomic-sequences requires a FASTA to read the contig sequences"
                .to_string(),
        })?;
        if contig_required_len.is_empty() {
            // The flag was requested but no emitted transcript carries a genomic
            // placement, so there is nothing to back a genome with. Signal the
            // caller to warn rather than write an empty `genomic_sequences` map
            // that would look genome-capable but report `has_genomic_data() == false`.
            emit_requested_no_placement = true;
        } else {
            let mut genomic_sequences = Map::new();
            let mut total_bytes: u64 = 0;
            for (contig, required_len) in &contig_required_len {
                let length = fasta
                    .sequence_length(contig)
                    .ok_or_else(|| FerroError::Io {
                        msg: format!(
                            "--emit-genomic-sequences: a transcript is placed on contig '{}', \
                         but it is not present in the FASTA",
                            contig
                        ),
                    })?;
                // Fail fast, with the SAME invariant the loader enforces, so the
                // emitted file always loads. A shorter contig means the FASTA does
                // not cover the annotation's coordinates — typically a sub-region
                // FASTA paired with whole-genome coordinates, or a mismatched
                // assembly.
                if length < *required_len {
                    return Err(FerroError::Io {
                        msg: format!(
                            "--emit-genomic-sequences: contig '{}' is {} bases in the FASTA, but a \
                             transcript placement needs genomic coordinate {}. The FASTA does not \
                             cover the annotation's coordinates (a sub-region FASTA with whole-genome \
                             coordinates, or a different assembly?).",
                            contig, length, required_len
                        ),
                    });
                }
                let sequence = fasta.get_sequence(contig, 0, length)?;
                total_bytes += sequence.len() as u64;
                genomic_sequences.insert(contig.clone(), json!(sequence));
            }
            emitted_genomic_bytes = total_bytes;
            output_json["genomic_sequences"] = Value::Object(genomic_sequences);
        }
    }

    Ok(ConvertGffOutcome {
        json: output_json,
        report,
        emitted_genomic_bytes,
        emit_requested_no_placement,
    })
}
