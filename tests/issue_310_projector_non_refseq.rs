//! Issue #310 — protein prediction for non-RefSeq transcript IDs and
//! spec-compliant `p.` selector emission.
//!
//! Pins:
//! 1. With `protein_id` set, the projector emits a `p.` prediction using
//!    that accession even when the transcript ID is not RefSeq-style.
//! 2. Without `protein_id`, the projector falls back to using the
//!    transcript ID as the `p.` accession — never silently drops the
//!    prediction just because of an ID format.
//! 3. RefSeq `NM_*` transcripts still infer `NP_*` and now omit the
//!    gene-symbol selector from Display.
//! 4. The parser still accepts the legacy `NP_*(GENE):p.` form and
//!    preserves the selector text in `gene_symbol` for callers that
//!    read it programmatically.

use ferro_hgvs::data::cdot::CdotMapper;
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::{VariantProjection, VariantProjector};

/// Build a mock reference where the transcript ID does NOT match the
/// RefSeq `NM_*` convention. Mirrors the reproduction in the issue body
/// but lets the caller supply or omit `protein_id`.
fn make_provider(transcript_id: &str, protein_id: Option<&str>) -> MockProvider {
    let mut provider = MockProvider::new();
    let tx = Transcript::new(
        transcript_id.to_string(),
        Some("MYGENE".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::with_genomic(1, 1, 9, 1000, 1008)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    )
    .with_protein_id(protein_id.map(str::to_string));
    provider.add_transcript(tx);
    // Genomic sequence: 999 N's + "ATGCGCTAA" + 100 N's so g.1003 is the
    // 4th base of the CDS (= c.4, the first base of codon 2 "CGC" → Arg).
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    provider
}

fn project(provider: MockProvider, variant: &str, transcript_id: &str) -> VariantProjection {
    let cdot = CdotMapper::from_transcripts(provider.all_transcripts());
    let projector = Projector::new(cdot);
    let vp = VariantProjector::new(projector, provider);
    vp.project(variant, transcript_id)
        .expect("projection should succeed")
}

#[test]
fn non_refseq_id_with_explicit_protein_id_uses_protein_id() {
    let provider = make_provider("MYGENE-gene.1", Some("MYGENE-protein.1"));
    let result = project(provider, "chr1:g.1003C>A", "MYGENE-gene.1");
    let p = result
        .protein
        .as_ref()
        .expect("p. should be set via protein_id")
        .to_string();
    assert_eq!(p, "MYGENE-protein.1:p.(Arg2Ser)");
    assert!(!p.contains("(MYGENE)"), "no (GENE) selector on p.: {p}");
}

#[test]
fn non_refseq_id_without_protein_id_falls_back_to_transcript_id() {
    let provider = make_provider("MYGENE-gene.1", None);
    let result = project(provider, "chr1:g.1003C>A", "MYGENE-gene.1");
    let p = result
        .protein
        .as_ref()
        .expect("p. should be set via transcript-id fallback")
        .to_string();
    assert_eq!(p, "MYGENE-gene.1:p.(Arg2Ser)");
}

#[test]
fn refseq_nm_prefix_still_infers_np_and_omits_gene_selector() {
    let provider = make_provider("NM_TEST.1", None);
    let result = project(provider, "chr1:g.1003C>A", "NM_TEST.1");
    let p = result
        .protein
        .as_ref()
        .expect("p. should be set via NM_ → NP_ inference")
        .to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
    assert!(!p.contains("(MYGENE)"), "no (GENE) selector on p.: {p}");
}

#[test]
fn refseq_xm_prefix_still_infers_xp_and_omits_gene_selector() {
    let provider = make_provider("XM_TEST.1", None);
    let result = project(provider, "chr1:g.1003C>A", "XM_TEST.1");
    let p = result
        .protein
        .as_ref()
        .expect("p. should be set via XM_ → XP_ inference")
        .to_string();
    assert_eq!(p, "XP_TEST.1:p.(Arg2Ser)");
    assert!(!p.contains("(MYGENE)"), "no (GENE) selector on p.: {p}");
}

#[test]
fn parser_still_accepts_legacy_np_with_selector() {
    use ferro_hgvs::hgvs::variant::HgvsVariant;
    use ferro_hgvs::parse_hgvs;
    let v = parse_hgvs("NP_TEST.1(MYGENE):p.(Arg2Ser)").expect("legacy form should still parse");
    let p = match v {
        HgvsVariant::Protein(p) => p,
        other => panic!("expected ProteinVariant, got {other:?}"),
    };
    // The parser preserves the selector text in `gene_symbol` for callers
    // that read it programmatically, even though Display now drops it.
    assert_eq!(p.gene_symbol.as_deref(), Some("MYGENE"));
    assert_eq!(p.to_string(), "NP_TEST.1:p.(Arg2Ser)");
}
