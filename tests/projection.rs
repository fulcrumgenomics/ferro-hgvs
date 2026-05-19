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

// =============================================================================
// Issue #332: `VariantProjector` must also route transcript lookups through
// the variant-aware `ReferenceProvider::get_transcript_for_variant` so the
// projector benefits from the same NG/NC-parent build resolution as
// `Normalizer::normalize`. Confirms the projector resolves an intronic c.
// input that carries an NG_* parent without erroring on the codon-fetch step.
// =============================================================================

/// Intronic-aware fixture: a 3-exon transcript on chr1 plus strand, intron
/// homopolymers, so intronic 3' shifting is well-defined. The projector path
/// only requires that the substitution-codon-fetch and indel-codon-fetch
/// sites use the variant-aware lookup; this fixture exercises an exonic g.
/// substitution against an NG-parented HGVS form to confirm the projector's
/// `provider.get_transcript_for_variant(variant)` plumbing is in place.
fn intronic_fixture() -> (Projector, MockProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // 3 exons, 30bp total, plus 2 introns of 50bp each.
            // Exon 1: g.[1001,1010] tx [0,10)
            // Exon 2: g.[1061,1070] tx [10,20)
            // Exon 3: g.[1121,1130] tx [20,30)
            exons: vec![
                [1000, 1010, 0, 10],
                [1060, 1070, 10, 20],
                [1120, 1130, 20, 30],
            ],
            cds_start: Some(0),
            cds_end: Some(30),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    let tx_seq = "ATGCGCTAAATGCGCTAAATGCGCTAACGT"; // 30 bases
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        TxStrand::Plus,
        tx_seq.to_string(),
        Some(1),
        Some(30),
        vec![
            Exon::with_genomic(1, 1, 10, 1001, 1010),
            Exon::with_genomic(2, 11, 20, 1061, 1070),
            Exon::with_genomic(3, 21, 30, 1121, 1130),
        ],
        Some("chr1".to_string()),
        Some(1001),
        Some(1130),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Plant a deterministic genomic sequence on chr1.
    let mut bytes = vec![b'N'; 2000];
    let plant = |bytes: &mut Vec<u8>, start_1: usize, seq: &str| {
        for (i, b) in seq.bytes().enumerate() {
            bytes[start_1 - 1 + i] = b;
        }
    };
    plant(&mut bytes, 1001, "ATGCGCTAAA");
    plant(&mut bytes, 1011, &"A".repeat(50));
    plant(&mut bytes, 1061, "TGCGCTAAAT");
    plant(&mut bytes, 1071, &"T".repeat(50));
    plant(&mut bytes, 1121, "GCGCTAACGT");
    provider.add_genomic_sequence("chr1", String::from_utf8(bytes).unwrap());
    (projector, provider)
}

#[test]
fn project_intronic_c_dot_succeeds() {
    // Pins that intronic g.→c. projection runs through `Normalizer::normalize`
    // (intronic path, which uses the variant-aware lookup after #332) and the
    // projector's per-codon transcript fetches without erroring.
    //
    // The input here is `g.`-form (the projector's primary supported entry
    // point); the per-codon transcript lookups inside the projector consult
    // `get_transcript_for_variant` on internally-constructed `CdsVariant`s
    // (which do NOT carry `genomic_context` today — that's a separate gap
    // tracked by #328). Direct NG-parented projection coverage will land
    // once #328 ships.
    let (projector, provider) = intronic_fixture();
    let vp = VariantProjector::new(projector, provider);
    let result = vp
        .project("NC_000001.11:g.1012A>T", "NM_TEST.1")
        .expect("projection must succeed for intronic input");
    assert_eq!(result.transcript_id, "NM_TEST.1");
    assert!(result.is_intronic);
    let c = result.coding.as_ref().unwrap().to_string();
    // Pin the position number — a regression that mis-identifies the
    // intronic offset would still match a loose `c.10+` substring.
    assert!(
        c.contains("c.10+3"),
        "expected intronic offset c.10+3; got {}",
        c
    );
    // Intronic substitutions skip protein prediction.
    assert!(result.protein.is_none(), "no p. for intronic substitutions");
}
