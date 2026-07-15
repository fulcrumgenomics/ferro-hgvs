//! Single-construct `transcripts.json` builder from a FASTA + CDS bounds.
//!
//! [`build_transcript`] is the library entry point behind the
//! `ferro build-transcript` CLI command and the `ferro_hgvs.build_transcript`
//! Python binding. It reads one contig from a FASTA, wraps it as a single-exon
//! transcript with the given CDS bounds, and serializes the `transcripts.json`
//! object shape consumed by [`JsonProvider`] and `Normalizer(reference_json=...)`.
//!
//! Both the CLI (`run_build_transcript`) and the Python binding call this one
//! function so their output is byte-identical for the same inputs/flags — the
//! only difference is I/O framing (where the JSON is written) and the stderr
//! messages the CLI emits.
//!
//! [`JsonProvider`]: crate::reference::mock::JsonProvider

use std::path::Path;

use serde_json::{json, Map, Value};

use crate::error::FerroError;
use crate::reference::fasta::FastaProvider;
use crate::reference::provider::ReferenceProvider;
use crate::sequence::reverse_complement;

/// Configuration for [`build_transcript`], mirroring the `ferro build-transcript`
/// CLI flags.
///
/// Construct with [`BuildTranscriptConfig::new`] and set fields as needed, or
/// build one directly with a struct literal; all fields are public.
#[derive(Debug, Clone)]
pub struct BuildTranscriptConfig {
    /// CDS start position (1-based inclusive, in transcript coordinates).
    pub cds_start: u64,
    /// CDS end position (1-based inclusive, in transcript coordinates).
    pub cds_end: u64,
    /// Transcript ID; defaults to the FASTA contig name when `None`.
    pub id: Option<String>,
    /// Strand: `"+"` or `"-"`. Minus-strand sequences are reverse-complemented,
    /// so `cds_start`/`cds_end` are relative to the reverse-complemented sequence.
    pub strand: String,
    /// Contig to use when the FASTA has multiple contigs; `None` requires a
    /// single-contig FASTA.
    pub contig: Option<String>,
    /// Optional gene symbol to embed in the transcript record.
    pub gene: Option<String>,
    /// Genome build name embedded verbatim in the output.
    pub genome_build: String,
    /// Also emit a top-level `genomic_sequences` map (the contig's forward bytes)
    /// so the output is genome-capable (#1026).
    pub emit_genomic_sequences: bool,
}

impl BuildTranscriptConfig {
    /// A default configuration for the given CDS bounds: plus strand, GRCh38,
    /// contig/id/gene auto-selected, no genomic-sequence embedding.
    pub fn new(cds_start: u64, cds_end: u64) -> Self {
        Self {
            cds_start,
            cds_end,
            id: None,
            strand: "+".to_string(),
            contig: None,
            gene: None,
            genome_build: "GRCh38".to_string(),
            emit_genomic_sequences: false,
        }
    }
}

/// Result of a [`build_transcript`] run.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct BuildTranscriptOutcome {
    /// The `transcripts.json` document (a JSON object with `version`,
    /// `genome_build`, `transcripts`, and — when emitted — `genomic_sequences`).
    pub json: Value,
    /// The transcript ID that was emitted (the configured `id` or the contig name).
    pub transcript_id: String,
    /// Total bytes of genomic sequence embedded under `genomic_sequences`
    /// (`0` when `emit_genomic_sequences` was not requested). Callers may warn
    /// when this is large.
    pub emitted_genomic_bytes: u64,
}

/// Build a single-construct `transcripts.json` from a FASTA contig + CDS bounds.
///
/// # Arguments
///
/// * `fasta_path` - Path to the FASTA (indexed or plain; the index is built on
///   the fly if absent).
/// * `config` - CDS bounds, strand, contig/id/gene selection, and genome build;
///   see [`BuildTranscriptConfig`].
///
/// # Errors
///
/// Returns [`FerroError`] if the FASTA cannot be opened, the requested/implied
/// contig cannot be resolved, the CDS bounds are out of range, or the strand is
/// not `"+"`/`"-"`.
pub fn build_transcript(
    fasta_path: &Path,
    config: &BuildTranscriptConfig,
) -> Result<BuildTranscriptOutcome, FerroError> {
    let fasta = FastaProvider::new(fasta_path)?;

    // Collect contig names; HashMap iteration order is non-deterministic but we
    // only need the names to detect single vs. multi-contig.
    let contig_names: Vec<String> = fasta.sequence_names().cloned().collect();

    let contig_name: String = match (config.contig.as_deref(), contig_names.as_slice()) {
        (Some(c), _) => {
            if !fasta.has_sequence(c) {
                return Err(FerroError::Io {
                    msg: format!(
                        "Contig '{}' not found in FASTA; available: {}",
                        c,
                        contig_names.join(", ")
                    ),
                });
            }
            c.to_string()
        }
        (None, [single]) => single.clone(),
        (None, []) => {
            return Err(FerroError::Io {
                msg: "FASTA contains no contigs".to_string(),
            })
        }
        (None, _) => {
            return Err(FerroError::Io {
                msg: format!(
                    "FASTA has {} contigs; set `contig` to select one (available: {})",
                    contig_names.len(),
                    contig_names.join(", ")
                ),
            })
        }
    };

    let contig_length = fasta
        .sequence_length(&contig_name)
        .ok_or_else(|| FerroError::Io {
            msg: format!("Could not determine length of contig '{}'", contig_name),
        })?;

    // Validate CDS bounds (1-based inclusive).
    if config.cds_start < 1 || config.cds_end < config.cds_start || config.cds_end > contig_length {
        return Err(FerroError::InvalidCoordinates {
            msg: format!(
                "Invalid CDS bounds: cds_start={}, cds_end={}, contig length={}",
                config.cds_start, config.cds_end, contig_length
            ),
        });
    }

    // Parse and validate strand.
    if config.strand != "+" && config.strand != "-" {
        return Err(FerroError::InvalidCoordinates {
            msg: format!("Invalid strand '{}'; expected + or -", config.strand),
        });
    }

    // Read the full contig sequence (0-based half-open for get_sequence). This is
    // the forward genomic sequence; the exon `genomic_start`/`genomic_end`
    // coordinates index into it, so it is what we emit under `genomic_sequences`.
    let mut genomic_forward = fasta.get_sequence(&contig_name, 0, contig_length)?;

    // Reverse-complement for minus-strand transcripts (transcript-space `sequence`).
    // On the plus strand the transcript sequence equals the forward contig: clone it
    // only when we still need `genomic_forward` for `genomic_sequences`, otherwise
    // move it out to avoid duplicating the whole contig.
    let final_sequence = if config.strand == "-" {
        reverse_complement(&genomic_forward)
    } else if config.emit_genomic_sequences {
        genomic_forward.clone()
    } else {
        std::mem::take(&mut genomic_forward)
    };

    let tx_id = config.id.clone().unwrap_or_else(|| contig_name.clone());

    // Build the exon object: single exon spanning the full contig. For a
    // single-exon synthetic construct, plus-strand transcript coordinates equal
    // genomic coordinates; for minus-strand the emitted `sequence` is
    // reverse-complemented and CDS positions are in transcript space.
    let exon = json!({
        "number": 1,
        "start": 1_u64,
        "end": contig_length,
        "genomic_start": 1_u64,
        "genomic_end": contig_length,
    });

    let mut tx_obj = json!({
        "id": tx_id,
        "strand": config.strand,
        "sequence": final_sequence,
        "cds_start": config.cds_start,
        "cds_end": config.cds_end,
        "exons": [exon],
        "chromosome": contig_name,
        "genomic_start": 1_u64,
        "genomic_end": contig_length,
        "genome_build": config.genome_build,
    });

    if let Some(ref g) = config.gene {
        tx_obj["gene_symbol"] = json!(g);
    }

    let mut output_json = json!({
        "version": "1.0",
        "genome_build": config.genome_build,
        "transcripts": [tx_obj],
    });

    // Optionally emit the contig's forward sequence so the transcripts.json is
    // genome-capable (#1026): the single exon's genomic coordinates index into it.
    // Note this repeats the contig bytes already carried in the transcript
    // `sequence` (identical on plus strand; reverse-complemented there on minus),
    // which is the price of a self-contained genome-capable reference.
    let mut emitted_genomic_bytes: u64 = 0;
    if config.emit_genomic_sequences {
        emitted_genomic_bytes = genomic_forward.len() as u64;
        let mut genomic_sequences = Map::new();
        genomic_sequences.insert(contig_name.clone(), json!(genomic_forward));
        output_json["genomic_sequences"] = Value::Object(genomic_sequences);
    }

    Ok(BuildTranscriptOutcome {
        json: output_json,
        transcript_id: tx_id,
        emitted_genomic_bytes,
    })
}
