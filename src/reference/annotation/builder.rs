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
    _format: AnnotationFormat,
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

    let cds_ranges: Vec<(u64, u64)> = cds_fids.iter().map(|f| graph.feature(*f).range).collect();
    let (cds_g_start, cds_g_end) = if !cds_ranges.is_empty() {
        let s = cds_ranges.iter().map(|(s, _)| *s).min().unwrap();
        let e = cds_ranges.iter().map(|(_, e)| *e).max().unwrap();
        (Some(s), Some(e))
    } else {
        (None, None)
    };

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

    Some(Transcript {
        id: tx_id_str,
        gene_symbol: None, // Phase 2 wires this
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
        mane_status: ManeStatus::None,
        refseq_match: None,
        ensembl_match: None,
        exon_cigars: Vec::new(),
        cached_introns: OnceLock::new(),
    })
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
            graph.ingest(&r);
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
        // CDS in transcript space: genomic 200 → tx 101; genomic 1000 → tx 901.
        assert_eq!(tx.cds_start, Some(101));
        assert_eq!(tx.cds_end, Some(901));
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
        // is at genomic 150 (in exon 2, tx position 252).
        assert_eq!(tx.cds_start, Some(51));
        assert_eq!(tx.cds_end, Some(252));
    }
}
