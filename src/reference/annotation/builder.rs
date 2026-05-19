//! Transcript construction from feature graph. See spec §6 Stage 4.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use super::diagnostics::{DiagnosticPayload, LoaderDiagnostic, LoaderReport, SourceLocation};
use super::feature::FeatureType;
use super::format_detect::AnnotationFormat;
use super::graph::{FeatureGraph, FeatureId};
use crate::reference::transcript::{Exon, GenomeBuild, ManeStatus, Strand, Transcript};

pub fn build_transcripts(
    graph: &FeatureGraph,
    format: AnnotationFormat,
    genome_build: GenomeBuild,
    source_path: PathBuf,
    report: &mut LoaderReport,
) -> Vec<Transcript> {
    let mut out = Vec::new();

    // Identify transcript-shaped features.
    let mut tx_fids: Vec<FeatureId> = (0..graph.features.len() as u32)
        .map(FeatureId)
        .filter(|fid| graph.feature(*fid).feature_type.is_transcript_like())
        .collect();

    // Also accept gene → CDS prokaryotic form: gene with no transcript-like children
    // but with at least one CDS child.
    for i in 0..graph.features.len() {
        let fid = FeatureId(i as u32);
        let f = graph.feature(fid);
        if !f.feature_type.is_gene_like() {
            continue;
        }
        let children = graph.children_of(fid);
        if children.is_empty() {
            continue;
        }
        let has_any_tx_child = children
            .iter()
            .any(|c| graph.feature(*c).feature_type.is_transcript_like());
        if has_any_tx_child {
            continue;
        }
        let has_cds = children
            .iter()
            .any(|c| matches!(graph.feature(*c).feature_type, FeatureType::Cds));
        if has_cds {
            tx_fids.push(fid);
            report.record(LoaderDiagnostic::warning(
                "W-LOAD-101",
                format!(
                    "GeneAsTranscript: synthesizing transcript from gene {}",
                    f.id.as_deref().unwrap_or("<no-id>")
                ),
                SourceLocation {
                    path: source_path.clone(),
                    line: f.source_line,
                },
                f.id.as_ref().map(|s| s.to_string()),
                DiagnosticPayload::GeneAsTranscript {
                    gene_id: f.id.as_deref().unwrap_or("").to_string(),
                },
            ));
        }
    }

    for tx_fid in tx_fids {
        match build_single(graph, tx_fid, format, genome_build, &source_path, report) {
            Some(t) => {
                report.transcripts_loaded += 1;
                out.push(t);
            }
            None => report.records_dropped += 1,
        }
    }
    out
}

fn build_single(
    graph: &FeatureGraph,
    tx_fid: FeatureId,
    format: AnnotationFormat,
    genome_build: GenomeBuild,
    source_path: &Path,
    report: &mut LoaderReport,
) -> Option<Transcript> {
    let tx_feat = graph.feature(tx_fid);
    let tx_id_str = tx_feat.id.as_deref().unwrap_or("<no-id>").to_string();

    if tx_feat.strand == Strand::Unknown {
        report.record(LoaderDiagnostic::error(
            "E-LOAD-103",
            format!("StrandRequired: transcript {} has strand '.'", tx_id_str),
            SourceLocation {
                path: source_path.to_path_buf(),
                line: tx_feat.source_line,
            },
            Some(tx_id_str.clone()),
            DiagnosticPayload::StrandRequired {
                feature_id: Some(tx_id_str.clone()),
            },
        ));
        return None;
    }

    // Exon-derivation ladder.
    let exon_fids = graph.descendants_of_type(tx_fid, &FeatureType::Exon);
    let cds_fids = graph.descendants_of_type(tx_fid, &FeatureType::Cds);

    let raw_exon_ranges: Vec<(u64, u64)> = if !exon_fids.is_empty() {
        // Step 1: explicit exons.
        exon_fids
            .iter()
            .map(|fid| graph.feature(*fid).range)
            .collect()
    } else {
        // Step 2 candidates: UTR records that, together with CDS, span the transcript.
        let utr5 = graph.descendants_of_type(tx_fid, &FeatureType::FivePrimeUtr);
        let utr3 = graph.descendants_of_type(tx_fid, &FeatureType::ThreePrimeUtr);
        let utr_generic = graph.descendants_of_type(tx_fid, &FeatureType::Utr);
        let has_utr = !utr5.is_empty() || !utr3.is_empty() || !utr_generic.is_empty();

        if has_utr && !cds_fids.is_empty() {
            // Step 2: merge UTR + CDS into synthetic exons.
            let mut all: Vec<(u64, u64)> = utr5
                .iter()
                .chain(utr3.iter())
                .chain(utr_generic.iter())
                .chain(cds_fids.iter())
                .map(|f| graph.feature(*f).range)
                .collect();
            all.sort_by_key(|(s, _)| *s);
            merge_adjacent(&all)
        } else if !cds_fids.is_empty() {
            // Step 3: CDS-as-exon (issue #183).
            //
            // A single CDS child is the only signal we have for UTR sequence:
            // the rest of the mRNA span is flanking UTR. Use the transcript's
            // own range as the single exon so 5'/3' UTRs survive the load.
            // Multiple CDS fragments instead encode the spliced exon structure
            // of a coding transcript, so we keep the existing fallback for
            // that case.
            if cds_fids.len() == 1 {
                vec![tx_feat.range]
            } else {
                cds_fids
                    .iter()
                    .map(|fid| graph.feature(*fid).range)
                    .collect()
            }
        } else {
            // Step 4: drop.
            report.record(LoaderDiagnostic::warning(
                "W-LOAD-100",
                format!(
                    "TranscriptWithoutExons: {} has no exon/UTR/CDS children",
                    tx_id_str
                ),
                SourceLocation {
                    path: source_path.to_path_buf(),
                    line: tx_feat.source_line,
                },
                Some(tx_id_str.clone()),
                DiagnosticPayload::NoExonsDerivable {
                    transcript_id: tx_id_str.clone(),
                },
            ));
            return None;
        }
    };

    // Order exons in transcript-space (5'→3'). For plus strand that's
    // genomic-ascending; for minus strand the 5'-most exon has the highest
    // genomic coordinates so we sort descending. The convention matches
    // cdot (see `genomic_to_tx_pos` in `src/data/cdot.rs`).
    let mut sorted_exons = raw_exon_ranges;
    if tx_feat.strand == Strand::Minus {
        sorted_exons.sort_by(|(a, _), (b, _)| b.cmp(a));
    } else {
        sorted_exons.sort_by_key(|(s, _)| *s);
    }

    // `derive_cds_bounds` needs ranges in genomic-ascending order to find the
    // exon containing a given genomic position; pass a genomic-sorted copy.
    let mut genomic_sorted_exons = sorted_exons.clone();
    genomic_sorted_exons.sort_by_key(|(s, _)| *s);

    let (cds_g_start, cds_g_end) = derive_cds_bounds(
        graph,
        tx_fid,
        &cds_fids,
        &genomic_sorted_exons,
        tx_feat.strand,
        &tx_id_str,
        format,
        source_path,
        report,
    );

    let mut tx_pos = 1u64;
    let exons: Vec<Exon> = sorted_exons
        .iter()
        .enumerate()
        .map(|(i, (g_start, g_end))| {
            let len = g_end - g_start + 1;
            let ex = Exon::with_genomic((i + 1) as u32, tx_pos, tx_pos + len - 1, *g_start, *g_end);
            tx_pos += len;
            ex
        })
        .collect();

    let (tx_cds_start, tx_cds_end) = match (cds_g_start, cds_g_end) {
        (Some(gs), Some(ge)) => {
            // `gs`/`ge` are min/max genomic CDS bounds; on minus strand `ge`
            // is the 5'-most CDS base (smaller tx position) and `gs` is the
            // 3'-most (larger tx position). Mirror the per-exon offset math
            // accordingly and order the result so cds_end >= cds_start.
            let minus = tx_feat.strand == Strand::Minus;
            let mut a = None;
            let mut b = None;
            for ex in &exons {
                let (egs, ege) = (ex.genomic_start.unwrap(), ex.genomic_end.unwrap());
                if gs >= egs && gs <= ege {
                    a = Some(if minus {
                        ex.start + (ege - gs)
                    } else {
                        ex.start + (gs - egs)
                    });
                }
                if ge >= egs && ge <= ege {
                    b = Some(if minus {
                        ex.start + (ege - ge)
                    } else {
                        ex.start + (ge - egs)
                    });
                }
            }
            match (a, b) {
                (Some(x), Some(y)) if x > y => (Some(y), Some(x)),
                _ => (a, b),
            }
        }
        _ => (None, None),
    };

    let chromosome = Some(tx_feat.seqid.to_string());
    let g_start = exons.iter().filter_map(|e| e.genomic_start).min();
    let g_end = exons.iter().filter_map(|e| e.genomic_end).max();

    let attrs = &tx_feat.attrs;
    // gene_symbol resolution differs by format. In GTF, the only spec-defined
    // symbol attributes are `gene_name` (preferred) and `gene_id` (last resort:
    // the ID itself, but better than nothing). In GFF3, `Name` is a generic
    // display attribute that some sources use for the gene symbol, and `gene=`
    // is a GFF3 convention used by RefSeq.
    let gene_symbol = match format {
        AnnotationFormat::Gtf => {
            crate::reference::annotation::feature::attr_get(attrs, "gene_name")
                .or_else(|| crate::reference::annotation::feature::attr_get(attrs, "gene_id"))
                .map(String::from)
        }
        AnnotationFormat::Gff3 => {
            crate::reference::annotation::feature::attr_get(attrs, "gene_name")
                .or_else(|| crate::reference::annotation::feature::attr_get(attrs, "gene"))
                .or_else(|| crate::reference::annotation::feature::attr_get(attrs, "Name"))
                .map(String::from)
        }
    };

    let mane_status = {
        let tag = crate::reference::annotation::feature::attr_get(attrs, "tag").unwrap_or("");
        let mane = crate::reference::annotation::feature::attr_get(attrs, "MANE")
            .or_else(|| crate::reference::annotation::feature::attr_get(attrs, "mane_status"))
            .unwrap_or("");
        if tag.contains("MANE_Select")
            || tag.contains("MANE Select")
            || mane.to_lowercase().contains("select")
        {
            ManeStatus::Select
        } else if tag.contains("MANE_Plus_Clinical")
            || tag.contains("MANE Plus Clinical")
            || mane.to_lowercase().contains("plus")
            || mane.to_lowercase().contains("clinical")
        {
            ManeStatus::PlusClinical
        } else {
            ManeStatus::None
        }
    };

    let refseq_match = crate::reference::annotation::feature::attr_get(attrs, "RefSeq")
        .or_else(|| crate::reference::annotation::feature::attr_get(attrs, "refseq_id"))
        .map(String::from);

    let ensembl_match = crate::reference::annotation::feature::attr_get(attrs, "Ensembl")
        .or_else(|| crate::reference::annotation::feature::attr_get(attrs, "ensembl_id"))
        .map(String::from);

    // Protein accession resolution. GTF carries `protein_id "NP_..."` on each
    // CDS record (and sometimes on the transcript record). GFF3 carries
    // `protein_id=NP_...` on each CDS feature (per GenBank/RefSeq convention).
    // Prefer the transcript-level attribute when present (authoritative);
    // otherwise fall back to the CDS children. All CDS fragments of a single
    // transcript should share the same protein_id (they're all one
    // polypeptide); if they disagree, we emit W-LOAD-120 and leave protein_id
    // unset rather than silently pick an input-order-dependent first value.
    // The projector then falls back to its transcript-id heuristic (see #310).
    let protein_id = if let Some(tx_pid) =
        crate::reference::annotation::feature::attr_get(attrs, "protein_id")
    {
        Some(tx_pid.to_string())
    } else {
        let mut distinct: Vec<String> = Vec::new();
        for fid in &cds_fids {
            if let Some(pid) = crate::reference::annotation::feature::attr_get(
                &graph.feature(*fid).attrs,
                "protein_id",
            ) {
                if !distinct.iter().any(|v| v == pid) {
                    distinct.push(pid.to_string());
                }
            }
        }
        match distinct.len() {
            0 => None,
            1 => Some(distinct.remove(0)),
            _ => {
                report.record(LoaderDiagnostic::warning(
                    "W-LOAD-120",
                    format!(
                        "ConflictingProteinId: transcript {} CDS fragments report \
                         {} distinct protein_id values ({}); protein_id left unset",
                        tx_id_str,
                        distinct.len(),
                        distinct.join(", "),
                    ),
                    SourceLocation {
                        path: source_path.to_path_buf(),
                        line: tx_feat.source_line,
                    },
                    Some(tx_id_str.clone()),
                    DiagnosticPayload::ConflictingProteinId {
                        transcript_id: tx_id_str.clone(),
                        values: distinct,
                    },
                ));
                None
            }
        }
    };

    Some(Transcript {
        id: tx_id_str,
        gene_symbol,
        strand: tx_feat.strand,
        // Coordinate-only load: actual bases are supplied later by a FASTA
        // provider. `sequence_length()` derives from exons, so callers that
        // need transcript length still work.
        sequence: None,
        cds_start: tx_cds_start,
        cds_end: tx_cds_end,
        exons,
        chromosome,
        genomic_start: g_start,
        genomic_end: g_end,
        genome_build,
        mane_status,
        refseq_match,
        ensembl_match,
        protein_id,
        exon_cigars: Vec::new(),
        cached_introns: OnceLock::new(),
    })
}

/// Merge sorted, possibly overlapping, possibly adjacent (touching) intervals into
/// the minimum set of disjoint ranges. Input must be sorted by start.
fn merge_adjacent(ranges: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if ranges.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(ranges.len());
    let mut cur = ranges[0];
    for &(s, e) in &ranges[1..] {
        // Adjacent or overlapping: extend.
        if s <= cur.1 + 1 {
            cur.1 = cur.1.max(e);
        } else {
            out.push(cur);
            cur = (s, e);
        }
    }
    out.push(cur);
    out
}

/// Derive `(cds_g_start_genomic_lo, cds_g_end_genomic_hi)` as a `(lo, hi)` pair in
/// genomic (not biological 5'→3') coordinate order.
///
/// - `cds_start_genomic`: the 5'-most coding base. Uses `start_codon` feature when
///   present; otherwise applies CDS phase to shift the leading CDS edge. Emits
///   W-LOAD-110 (phase applied) or W-LOAD-111 (phase unavailable) when phase is used.
/// - `cds_end_genomic`: the 3'-most coding base **including** the stop codon, matching
///   ferro's downstream convention. Uses `stop_codon` feature when present. For GTF
///   the stop is already included in CDS. For GFF3 the CDS 3' edge is extended by 3 bp
///   (clipped at the exon boundary). Emits W-LOAD-112 (StopCodonAssumed) when
///   extending.
///
/// `sorted_exon_ranges` must be sorted by ascending genomic start; the clipping
/// step uses it to look up the exon containing the 3'-most CDS base.
///
/// Returns `(None, None)` when there are no CDS features.
#[allow(clippy::too_many_arguments)]
fn derive_cds_bounds(
    graph: &FeatureGraph,
    tx_fid: FeatureId,
    cds_fids: &[FeatureId],
    sorted_exon_ranges: &[(u64, u64)],
    strand: Strand,
    tx_id: &str,
    format: AnnotationFormat,
    source_path: &Path,
    report: &mut LoaderReport,
) -> (Option<u64>, Option<u64>) {
    if cds_fids.is_empty() {
        return (None, None);
    }

    let start_codons = graph.descendants_of_type(tx_fid, &FeatureType::StartCodon);
    let stop_codons = graph.descendants_of_type(tx_fid, &FeatureType::StopCodon);

    let plus = matches!(strand, Strand::Plus);

    // 5'-most CDS edge in genomic coords.
    let cds_5prime_genomic: u64 = if plus {
        cds_fids
            .iter()
            .map(|f| graph.feature(*f).range.0)
            .min()
            .unwrap()
    } else {
        cds_fids
            .iter()
            .map(|f| graph.feature(*f).range.1)
            .max()
            .unwrap()
    };

    // Find the "leading" CDS feature (the one whose 5'-most edge matches cds_5prime_genomic).
    let leading_cds = cds_fids
        .iter()
        .copied()
        .find(|f| {
            let r = graph.feature(*f).range;
            if plus {
                r.0 == cds_5prime_genomic
            } else {
                r.1 == cds_5prime_genomic
            }
        })
        .unwrap();

    // Step 1: resolve cds_start_genomic.
    let resolved_start: u64 = if let Some(sc) = start_codons.first() {
        let r = graph.feature(*sc).range;
        if plus {
            r.0
        } else {
            r.1
        }
    } else {
        let phase = graph.feature(leading_cds).phase;
        // Phase per GFF3/GTF spec is exactly one of {0, 1, 2}. The record
        // parser stores phase as `u8`, so out-of-range values (e.g. a
        // malformed `3`) survive parsing. Treat them as if phase were
        // missing rather than blindly shifting the CDS edge by the
        // attacker-controlled offset.
        match phase {
            Some(0) => cds_5prime_genomic,
            Some(p @ (1 | 2)) => {
                report.record(LoaderDiagnostic::warning(
                    "W-LOAD-110",
                    format!("PhaseAppliedToCdsStart: tx {} phase={}", tx_id, p),
                    SourceLocation {
                        path: source_path.to_path_buf(),
                        line: graph.feature(leading_cds).source_line,
                    },
                    Some(tx_id.into()),
                    DiagnosticPayload::PhaseAppliedToCdsStart {
                        transcript_id: tx_id.into(),
                        phase: p,
                    },
                ));
                if plus {
                    cds_5prime_genomic + p as u64
                } else {
                    cds_5prime_genomic.saturating_sub(p as u64)
                }
            }
            Some(_) | None => {
                report.record(LoaderDiagnostic::warning(
                    "W-LOAD-111",
                    format!("PhaseUnavailable: tx {}", tx_id),
                    SourceLocation {
                        path: source_path.to_path_buf(),
                        line: graph.feature(leading_cds).source_line,
                    },
                    Some(tx_id.into()),
                    DiagnosticPayload::PhaseUnavailable {
                        transcript_id: tx_id.into(),
                    },
                ));
                cds_5prime_genomic
            }
        }
    };

    // Step 2: resolve cds_end_genomic (3'-most coding base, stop INCLUDED).
    let cds_3prime_genomic: u64 = if plus {
        cds_fids
            .iter()
            .map(|f| graph.feature(*f).range.1)
            .max()
            .unwrap()
    } else {
        cds_fids
            .iter()
            .map(|f| graph.feature(*f).range.0)
            .min()
            .unwrap()
    };

    let resolved_end: u64 = if let Some(stop) = stop_codons.first() {
        let r = graph.feature(*stop).range;
        if plus {
            r.1
        } else {
            r.0
        }
    } else if matches!(format, AnnotationFormat::Gtf) {
        // GTF already includes stop codon in CDS.
        cds_3prime_genomic
    } else {
        // GFF3: extend by 3 bp on the 3' side, clipping at exon boundary.
        report.record(LoaderDiagnostic::warning(
            "W-LOAD-112",
            format!("StopCodonAssumed: extending CDS for tx {}", tx_id),
            SourceLocation {
                path: source_path.to_path_buf(),
                line: graph.feature(leading_cds).source_line,
            },
            Some(tx_id.into()),
            DiagnosticPayload::StopCodonAssumed {
                transcript_id: tx_id.into(),
            },
        ));
        if plus {
            // Find the exon containing cds_3prime_genomic and clip the extension to its end.
            let extended = cds_3prime_genomic + 3;
            sorted_exon_ranges
                .iter()
                .find(|(es, ee)| cds_3prime_genomic >= *es && cds_3prime_genomic <= *ee)
                .map(|(_, ee)| extended.min(*ee))
                .unwrap_or(cds_3prime_genomic)
        } else {
            let extended = cds_3prime_genomic.saturating_sub(3);
            sorted_exon_ranges
                .iter()
                .find(|(es, ee)| cds_3prime_genomic >= *es && cds_3prime_genomic <= *ee)
                .map(|(es, _)| extended.max(*es))
                .unwrap_or(cds_3prime_genomic)
        }
    };

    // Return as (lo, hi) by GENOMIC coordinate (not biological 5'→3').
    let (lo, hi) = if plus {
        (resolved_start, resolved_end)
    } else {
        (resolved_end, resolved_start)
    };
    (Some(lo), Some(hi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::annotation::format_detect::AnnotationFormat;
    use crate::reference::annotation::graph::FeatureGraph;
    use crate::reference::annotation::record::Gff3Record;
    use crate::reference::transcript::GenomeBuild;

    fn build_db(lines: &[&str], fmt: AnnotationFormat) -> (Vec<Transcript>, LoaderReport) {
        let mut graph = FeatureGraph::new();
        for (i, l) in lines.iter().enumerate() {
            let r = Gff3Record::parse(l, (i + 1) as u64).unwrap().unwrap();
            graph.ingest(r);
        }
        graph.resolve();
        let mut report = LoaderReport::default();
        let txs = build_transcripts(
            &graph,
            fmt,
            GenomeBuild::GRCh38,
            "test.gff".into(),
            &mut report,
        );
        (txs, report)
    }

    #[test]
    fn cds_start_uses_start_codon_when_present() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1",
            "chr1\t.\tstart_codon\t150\t152\t.\t+\t0\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        let tx = &txs[0];
        // exon starts at tx pos 1 = genomic 100, so genomic 150 = tx pos 51.
        assert_eq!(tx.cds_start, Some(51));
    }

    #[test]
    fn cds_start_applies_phase_when_no_start_codon() {
        // Leading CDS at 150 with phase 2 means the first in-frame base is at 152.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t2\tParent=tx1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        let tx = &txs[0];
        assert_eq!(tx.cds_start, Some(53)); // 150 + 2 = 152 → tx pos 53
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-110"), Some(&1));
    }

    /// Phase is spec-defined as {0, 1, 2} but the record parser stores it as
    /// `u8`, so a malformed value (here: 5) survives parsing. The CDS-bounds
    /// resolver must not blindly shift the CDS edge by an out-of-range phase;
    /// it must treat it as "phase unavailable" (W-LOAD-111) and emit no
    /// W-LOAD-110.
    #[test]
    fn cds_start_rejects_out_of_range_phase_as_unavailable() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t5\tParent=tx1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        let tx = &txs[0];
        // Phase ignored: cds_start_genomic stays at 150 → tx pos 51.
        assert_eq!(tx.cds_start, Some(51));
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-110"), None);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-111"), Some(&1));
    }

    #[test]
    fn gff3_cds_end_extends_to_include_stop_when_within_exon() {
        // GFF3 input, CDS 150..450 inside exon 100..500, no stop_codon record.
        // Stop extension by 3 bp → cds_end_genomic = 453, which is inside the exon.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        let tx = &txs[0];
        // genomic 453 → tx pos = 453 - 100 + 1 = 354
        assert_eq!(tx.cds_end, Some(354));
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-112"), Some(&1));
    }

    #[test]
    fn gff3_cds_end_clipped_when_no_room_for_stop_extension() {
        // Issue-#183 case: CDS = exon = mRNA = 100..1200, no stop_codon, no room.
        // Extension should clip to the exon end (1200), with W-LOAD-112.
        let lines = &[
            "seq1\t.\tgene\t100\t1200\t.\t+\t.\tID=gene01;Name=gene01",
            "seq1\t.\tmRNA\t100\t1200\t.\t+\t.\tID=gene01.1;Parent=gene01",
            "seq1\t.\tCDS\t100\t1200\t.\t+\t0\tParent=gene01.1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        let tx = &txs[0];
        // No room: cds_end_genomic clips to 1200 → tx pos 1200 - 100 + 1 = 1101.
        assert_eq!(tx.cds_end, Some(1101));
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-112"), Some(&1));
    }

    #[test]
    fn gff3_cds_end_uses_stop_codon_record_when_present() {
        // When stop_codon record exists, use it directly — no W-LOAD-112.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1",
            "chr1\t.\tstop_codon\t451\t453\t.\t+\t0\tParent=tx1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        let tx = &txs[0];
        assert_eq!(tx.cds_end, Some(354)); // genomic 453 → tx pos 354
                                           // No W-LOAD-112 because stop_codon was explicit.
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-112"), None);
    }

    #[test]
    fn gtf_cds_end_uses_max_cds_verbatim() {
        // GTF input: stop already included in CDS by convention.
        let mut graph = crate::reference::annotation::graph::FeatureGraph::new();
        let lines = &[
            "chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";",
            "chr1\tHAVANA\texon\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";",
            "chr1\tHAVANA\tCDS\t150\t450\t.\t+\t0\tgene_id \"g1\"; transcript_id \"tx1\";",
        ];
        for (i, l) in lines.iter().enumerate() {
            let r = crate::reference::annotation::record::GtfRecord::parse(l, (i + 1) as u64)
                .unwrap()
                .unwrap();
            graph.ingest(r);
        }
        graph.resolve();
        let mut report = LoaderReport::default();
        let txs = build_transcripts(
            &graph,
            AnnotationFormat::Gtf,
            crate::reference::transcript::GenomeBuild::GRCh38,
            "t.gtf".into(),
            &mut report,
        );
        let tx = &txs[0];
        // GTF: cds_end_genomic = 450 → tx pos 351. No W-LOAD-112.
        assert_eq!(tx.cds_end, Some(351));
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-112"), None);
    }

    #[test]
    fn ladder_step1_uses_explicit_exons() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;gene=G1",
            "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1",
            "chr1\t.\texon\t300\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].exons.len(), 2);
        assert_eq!(txs[0].exons[0].genomic_start, Some(100));
        assert_eq!(txs[0].exons[1].genomic_end, Some(500));
    }

    #[test]
    fn ladder_step3_cds_as_exon_for_issue_183() {
        let lines = &[
            "seq1\t.\tgene\t100\t1200\t.\t+\t.\tID=gene01;Name=gene01",
            "seq1\t.\tmRNA\t100\t1200\t.\t+\t.\tID=gene01.1;Parent=gene01",
            "seq1\t.\tCDS\t100\t1200\t.\t+\t0\tParent=gene01.1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(
            txs.len(),
            1,
            "issue #183: single-exon transcript must not be dropped"
        );
        assert_eq!(txs[0].exons.len(), 1);
        assert_eq!(txs[0].exons[0].genomic_start, Some(100));
        assert_eq!(txs[0].exons[0].genomic_end, Some(1200));
    }

    /// Issue #183, second leg: when a single-exon transcript has the
    /// `gene → mRNA → CDS` pattern and the CDS is smaller than the mRNA,
    /// the flanking UTR sequence must survive the load. Mapping CDS spans
    /// directly into exons truncates the UTRs and shifts cds_start/cds_end
    /// to the wrong transcript positions.
    #[test]
    fn ladder_step3_preserves_utr_when_cds_is_subset_of_mrna() {
        let lines = &[
            "seq1\t.\tgene\t100\t1200\t.\t+\t.\tID=g1;Name=g1",
            "seq1\t.\tmRNA\t100\t1200\t.\t+\t.\tID=tx1;Parent=g1",
            "seq1\t.\tCDS\t200\t1000\t.\t+\t0\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        assert_eq!(tx.exons.len(), 1);
        // The exon spans the entire mRNA, not just the CDS.
        assert_eq!(tx.exons[0].genomic_start, Some(100));
        assert_eq!(tx.exons[0].genomic_end, Some(1200));
        // Exon length = 1200 - 100 + 1.
        assert_eq!(tx.exons[0].end, 1101);
        // CDS in transcript space: genomic 200 → tx 101; genomic 1003 (stop-extended) → tx 904.
        assert_eq!(tx.cds_start, Some(101));
        assert_eq!(tx.cds_end, Some(904));
    }

    #[test]
    fn ladder_step4_drops_with_diagnostic() {
        let lines = &["chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1"];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        assert!(txs.is_empty());
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-100"), Some(&1));
    }

    #[test]
    fn unknown_strand_dropped_with_diagnostic() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t.\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t.\t.\tParent=tx1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        assert!(txs.is_empty());
        assert_eq!(report.diagnostics_by_code.get("E-LOAD-103"), Some(&1));
    }

    #[test]
    fn multi_parent_exon_is_attached_to_both_transcripts() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx2",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1,tx2",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs.len(), 2);
        assert!(txs.iter().all(|tx| tx.exons.len() == 1));
    }

    /// Minus-strand transcripts must be assembled in transcript order
    /// (5'→3'), which means the highest-genomic exon is exon 1, and CDS
    /// coordinate math mirrors `(ege - gpos)` instead of `(gpos - egs)`.
    /// See cdot's `genomic_to_tx_pos` for the canonical convention.
    #[test]
    fn minus_strand_exons_and_cds_use_transcript_order() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t-\t.\tID=tx1;gene=G1",
            "chr1\t.\texon\t100\t200\t.\t-\t.\tParent=tx1",
            "chr1\t.\texon\t300\t500\t.\t-\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t-\t0\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        // Exon 1 is the 5'-most exon, which on minus strand has the highest
        // genomic coordinates.
        assert_eq!(tx.exons[0].number, 1);
        assert_eq!(tx.exons[0].genomic_start, Some(300));
        assert_eq!(tx.exons[0].genomic_end, Some(500));
        assert_eq!(tx.exons[0].start, 1);
        assert_eq!(tx.exons[0].end, 201);
        assert_eq!(tx.exons[1].number, 2);
        assert_eq!(tx.exons[1].genomic_start, Some(100));
        assert_eq!(tx.exons[1].genomic_end, Some(200));
        assert_eq!(tx.exons[1].start, 202);
        assert_eq!(tx.exons[1].end, 302);
        // CDS spans genomic 150..=450. On minus strand, the 5'-most CDS base
        // is at genomic 450 (in exon 1, tx position 51); the 3'-most CDS base
        // (after GFF3 stop-codon extension shifts 150 → 147) is at genomic 147
        // (in exon 2, tx position 202 + (200 - 147) = 255).
        assert_eq!(tx.cds_start, Some(51));
        assert_eq!(tx.cds_end, Some(255));
    }

    #[test]
    fn ladder_step2_merges_utrs_and_cds_into_one_exon_when_adjacent() {
        // mRNA with no `exon` records but explicit UTR and CDS children.
        // The three intervals are adjacent so they merge to a single exon.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\tfive_prime_UTR\t100\t149\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1",
            "chr1\t.\tthree_prime_UTR\t451\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        assert_eq!(
            tx.exons.len(),
            1,
            "adjacent UTR+CDS+UTR should merge to one exon"
        );
        assert_eq!(tx.exons[0].genomic_start, Some(100));
        assert_eq!(tx.exons[0].genomic_end, Some(500));
    }

    #[test]
    fn ladder_step2_keeps_separate_exons_when_intervals_have_gaps() {
        // mRNA with UTR/CDS/UTR but with a genomic gap → multiple synthetic exons.
        let lines = &[
            "chr1\t.\tmRNA\t100\t800\t.\t+\t.\tID=tx1",
            "chr1\t.\tfive_prime_UTR\t100\t149\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t300\t.\t+\t0\tParent=tx1",
            "chr1\t.\tCDS\t500\t650\t.\t+\t0\tParent=tx1",
            "chr1\t.\tthree_prime_UTR\t651\t800\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs.len(), 1);
        let tx = &txs[0];
        // The 5'UTR + first CDS are adjacent (100..300), then a gap, then second CDS + 3'UTR adjacent (500..800).
        assert_eq!(tx.exons.len(), 2);
        assert_eq!(tx.exons[0].genomic_start, Some(100));
        assert_eq!(tx.exons[0].genomic_end, Some(300));
        assert_eq!(tx.exons[1].genomic_start, Some(500));
        assert_eq!(tx.exons[1].genomic_end, Some(800));
    }

    #[test]
    fn ladder_step2_with_only_utrs_no_cds_falls_through_to_drop() {
        // No CDS → step 2 doesn't fire; step 3 needs CDS too; step 4 drops.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\tfive_prime_UTR\t100\t250\t.\t+\t.\tParent=tx1",
            "chr1\t.\tthree_prime_UTR\t251\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        assert!(txs.is_empty());
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-100"), Some(&1));
    }

    #[test]
    fn gff3_gene_symbol_extracted_from_gene_attribute() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;gene=GENE1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].gene_symbol.as_deref(), Some("GENE1"));
    }

    #[test]
    fn gff3_gene_name_attribute_takes_precedence_over_name() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Name=fallback;gene_name=PRIMARY",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].gene_symbol.as_deref(), Some("PRIMARY"));
    }

    /// GTF transcripts without `gene_name` must still get a symbol via the
    /// `gene_id` fallback. The previous chain (`gene_name` → `gene` → `Name`)
    /// skipped `gene_id`, so common GTF inputs (e.g. GENCODE basic, Ensembl
    /// stripped down to gene_id/transcript_id) lost the gene symbol entirely.
    #[test]
    fn gtf_gene_symbol_falls_back_to_gene_id_when_no_gene_name() {
        let mut graph = crate::reference::annotation::graph::FeatureGraph::new();
        let lines = &[
            "chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tgene_id \"ENSG_NOSYMBOL\"; transcript_id \"tx1\";",
            "chr1\tHAVANA\texon\t100\t500\t.\t+\t.\tgene_id \"ENSG_NOSYMBOL\"; transcript_id \"tx1\";",
        ];
        for (i, l) in lines.iter().enumerate() {
            let r = crate::reference::annotation::record::GtfRecord::parse(l, (i + 1) as u64)
                .unwrap()
                .unwrap();
            graph.ingest(r);
        }
        graph.resolve();
        let mut report = LoaderReport::default();
        let txs = build_transcripts(
            &graph,
            AnnotationFormat::Gtf,
            crate::reference::transcript::GenomeBuild::GRCh38,
            "t.gtf".into(),
            &mut report,
        );
        assert_eq!(txs[0].gene_symbol.as_deref(), Some("ENSG_NOSYMBOL"));
    }

    #[test]
    fn mane_select_extracted_from_tag() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;tag=MANE_Select",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert!(matches!(
            txs[0].mane_status,
            crate::reference::transcript::ManeStatus::Select
        ));
    }

    #[test]
    fn mane_plus_clinical_extracted_from_tag() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;tag=MANE_Plus_Clinical",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert!(matches!(
            txs[0].mane_status,
            crate::reference::transcript::ManeStatus::PlusClinical
        ));
    }

    #[test]
    fn refseq_and_ensembl_cross_refs_extracted() {
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;RefSeq=NM_000088.3;Ensembl=ENST00000000001",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].refseq_match.as_deref(), Some("NM_000088.3"));
        assert_eq!(txs[0].ensembl_match.as_deref(), Some("ENST00000000001"));
    }

    #[test]
    fn protein_id_extracted_from_cds_feature_gff3() {
        // GFF3: RefSeq/GenBank put `protein_id=NP_...` on every CDS feature
        // of a coding transcript. The builder picks it up from the first CDS
        // child.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1;protein_id=NP_000088.3",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].protein_id.as_deref(), Some("NP_000088.3"));
    }

    #[test]
    fn protein_id_extracted_from_transcript_attr_when_present() {
        // Some GFF3 producers put `protein_id` directly on the mRNA record.
        // The transcript-level attribute wins over CDS-level when both are
        // present.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;protein_id=NP_TXLEVEL.1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1;protein_id=NP_CDSLEVEL.1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].protein_id.as_deref(), Some("NP_TXLEVEL.1"));
    }

    #[test]
    fn protein_id_extracted_from_first_cds_when_multiple_fragments() {
        // Multi-exon coding transcripts have one CDS feature per coding
        // exon; all CDS fragments of a single transcript share the same
        // `protein_id`. The first match suffices.
        let lines = &[
            "chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t300\t.\t+\t.\tParent=tx1",
            "chr1\t.\texon\t500\t900\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t300\t.\t+\t0\tParent=tx1;protein_id=NP_SPLICED.2",
            "chr1\t.\tCDS\t500\t850\t.\t+\t0\tParent=tx1;protein_id=NP_SPLICED.2",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].protein_id.as_deref(), Some("NP_SPLICED.2"));
    }

    #[test]
    fn protein_id_absent_when_no_attribute() {
        // No `protein_id` attribute anywhere → `None`, and the projector's
        // fallback chain (NM_ → NP_ inference, then transcript-id) takes
        // over downstream. See #310.
        let lines = &[
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1",
        ];
        let (txs, _) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].protein_id, None);
    }

    #[test]
    fn protein_id_warns_and_drops_when_cds_fragments_disagree() {
        // GFF/GTF convention: every CDS fragment of one transcript carries
        // the same `protein_id` (it's all one polypeptide). If the upstream
        // file disagrees we must not silently pick the first one and let
        // output depend on iteration order — emit W-LOAD-120 and leave
        // `protein_id` unset so the projector's transcript-id fallback
        // (see #310) takes over downstream.
        let lines = &[
            "chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=tx1",
            "chr1\t.\texon\t100\t300\t.\t+\t.\tParent=tx1",
            "chr1\t.\texon\t500\t900\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t300\t.\t+\t0\tParent=tx1;protein_id=NP_FIRST.1",
            "chr1\t.\tCDS\t500\t850\t.\t+\t0\tParent=tx1;protein_id=NP_SECOND.1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].protein_id, None);
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-120"), Some(&1));
    }

    #[test]
    fn protein_id_transcript_attr_wins_even_when_cds_fragments_disagree() {
        // A transcript-level `protein_id` is authoritative; when present we
        // do not look at CDS records (and so do not emit W-LOAD-120 for any
        // CDS-only disagreement).
        let lines = &[
            "chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=tx1;protein_id=NP_TXLEVEL.1",
            "chr1\t.\texon\t100\t300\t.\t+\t.\tParent=tx1",
            "chr1\t.\texon\t500\t900\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t150\t300\t.\t+\t0\tParent=tx1;protein_id=NP_A.1",
            "chr1\t.\tCDS\t500\t850\t.\t+\t0\tParent=tx1;protein_id=NP_B.1",
        ];
        let (txs, report) = build_db(lines, AnnotationFormat::Gff3);
        assert_eq!(txs[0].protein_id.as_deref(), Some("NP_TXLEVEL.1"));
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-120"), None);
    }

    #[test]
    fn protein_id_extracted_from_cds_feature_gtf() {
        // GTF: `protein_id "NP_..."` is the standard attribute on CDS
        // records produced by both GENCODE and RefSeq.
        let lines = &[
            "chr1\thavana\ttranscript\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";",
            "chr1\thavana\texon\t100\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";",
            "chr1\thavana\tCDS\t150\t450\t.\t+\t0\tgene_id \"g1\"; transcript_id \"tx1\"; protein_id \"ENSP00000000001.1\";",
        ];
        let mut graph = crate::reference::annotation::graph::FeatureGraph::new();
        for (i, l) in lines.iter().enumerate() {
            let r = crate::reference::annotation::record::GtfRecord::parse(l, (i + 1) as u64)
                .unwrap()
                .unwrap();
            graph.ingest(r);
        }
        graph.resolve();
        let mut report = LoaderReport::default();
        let txs = build_transcripts(
            &graph,
            AnnotationFormat::Gtf,
            GenomeBuild::GRCh38,
            "t.gtf".into(),
            &mut report,
        );
        assert_eq!(txs[0].protein_id.as_deref(), Some("ENSP00000000001.1"));
    }
}
