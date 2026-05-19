//! End-to-end tests for the variant-projection module.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{VariantProjection, VariantProjector};

/// Build a test (Projector, MockProvider) pair for NM_TEST.1.
///
/// The transcript encodes "ATGCGCTAA" (Met-Arg-Stop, 3 codons) on chr1 plus
/// strand at genomic positions 1000-1008 (1-based inclusive).
///
/// Cdot coordinate layout (0-based):
///   exon: genome [1000, 1009) — tx [0, 9)
///   cds_start = 0 (no 5' UTR), cds_end = 9
///
/// Coordinate mapping (cdot 0-based genome → HGVS c. 1-based):
///   g.1000 → tx_pos 0 → c.1  (A)
///   g.1001 → tx_pos 1 → c.2  (T)
///   g.1002 → tx_pos 2 → c.3  (G)
///   g.1003 → tx_pos 3 → c.4  (C ← ref base for the substitution test)
///   ...
fn fixture() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [genome_start(0-based), genome_end(0-based excl), tx_start, tx_end]
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Genomic sequence: 1000 N's + "ATGCGCTAA" + 100 N's.
    // 0-based index 1000 = 'A', 1001 = 'T', 1002 = 'G', 1003 = 'C', ...
    let prefix = "N".repeat(1000);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    (projector, provider)
}

#[test]
fn end_to_end_missense() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    // g.1003C>A: the 4th base (c.4) of codon 2 "CGC" (Arg) → "AGC" (Ser).
    let result: VariantProjection = vp
        .project("NC_000001.11:g.1003C>A", "NM_TEST.1")
        .expect("projection should succeed");
    assert_eq!(result.transcript_id, "NM_TEST.1");
    let c = result.coding.as_ref().unwrap().to_string();
    assert!(c.contains(":c.4C>A"), "got c. = {}", c);
    let p = result.protein.as_ref().unwrap().to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
    assert!(!result.is_frameshift);
    assert!(!result.is_intronic);
    assert!(!result.is_utr);
}

#[test]
fn end_to_end_no_overlap() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    // chr1:5000 is far outside the transcript at genome [1000, 1009).
    let err = vp.project("NC_000001.11:g.5000A>G", "NM_TEST.1");
    assert!(err.is_err(), "expected error for off-transcript position");
}
