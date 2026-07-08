//! #857 — C-terminal extension / stop-loss for a deletion spanning the
//! CDS→3'UTR boundary. CI-runnable (MockProvider, no manifest).
//!
//! Spec: `assets/hgvs-nomenclature/docs/recommendations/protein/extension.md`
//! (a variant disrupting the stop codon yields `p.(Ter<pos><aa>extTer<k>)`,
//! or `extTer?` with no downstream stop). ferro previously dropped a
//! CDS→3'UTR-spanning deletion (e.g. `c.9_*1del`) as an empty projection.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, VariantProjector};

/// `VariantProjector` for `NM_TEST.1` with CDS "ATGAAATAA" (Met-Lys-Ter, stop
/// at c.7_9) and the given 3'UTR appended after the CDS.
fn fixture(utr3: &str) -> VariantProjector<MockProvider> {
    let full = format!("ATGAAATAA{utr3}");
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // CdotTranscript.exons is Vec<[u64;4]>: [g_start, g_end, tx_start, tx_end].
            exons: vec![[1000, 1000 + full.len() as u64, 0, full.len() as u64]],
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
        full.clone(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, full.len() as u64)],
        Some("chr1".to_string()),
        Some(1000),
        Some(1000 + full.len() as u64 - 1),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    VariantProjector::new(projector, provider)
}

fn protein_of(vp: &VariantProjector<MockProvider>, input: &str) -> Option<String> {
    let v = parse_hgvs(input).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    let r = vp
        .project_variant(&v, "NM_TEST.1")
        .unwrap_or_else(|e| panic!("project {input}: {e}"));
    r.protein.as_ref().map(|p| p.to_string())
}

/// Boundary-spanning deletion with a downstream in-frame stop → `extTer1`.
#[test]
fn boundary_del_emits_cterminal_extension() {
    let vp = fixture("TCTAA");
    assert_eq!(
        protein_of(&vp, "NM_TEST.1:c.9_*1del").as_deref(),
        Some("NP_TEST.1:p.(Ter3TyrextTer1)")
    );
}

/// Boundary-spanning deletion with no downstream stop → `extTer?`.
#[test]
fn boundary_del_no_downstream_stop_is_ext_ter_unknown() {
    let vp = fixture("TCGGG");
    assert_eq!(
        protein_of(&vp, "NM_TEST.1:c.9_*1del").as_deref(),
        Some("NP_TEST.1:p.(Ter3TyrextTer?)")
    );
}

/// Regression: an in-CDS stop-disrupting deletion (no `*N` end) must still route
/// through the existing in-CDS path and produce an extension — guards against
/// breaking that path or double-appending the 3'UTR.
#[test]
fn in_cds_stop_disruption_still_uses_existing_path() {
    let vp = fixture("TCTAA");
    // c.7del drops the first base of the stop codon (TAA at c.7_9); the reading
    // frame reads through the 3'UTR with no downstream in-frame stop in this
    // fixture, so the consequence is a fully-specified extension to an unknown
    // terminator. Pin the complete `p.` string so a wrong extension body (e.g.
    // a double-appended 3'UTR) fails the test.
    assert_eq!(
        protein_of(&vp, "NM_TEST.1:c.7del").as_deref(),
        Some("NP_TEST.1:p.(Ter3AsnextTer?)")
    );
}

/// A purely-3'UTR deletion (both ends `*N`, base > 3 to avoid the utr3-blind
/// init-codon path) has no protein consequence.
#[test]
fn pure_3utr_deletion_has_no_protein_consequence() {
    let vp = fixture("TCTAATTTGGG");
    assert_eq!(protein_of(&vp, "NM_TEST.1:c.*4_*5del"), None);
}
