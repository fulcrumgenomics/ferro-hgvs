//! #867: `VariantProjector::project_to_genomic_normalized` returns the
//! spec-canonical 3'-shifted genomic form, vs the raw (un-normalized) pivot
//! `project_to_genomic`. Hermetic MockProvider, no manifest.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::Transcript;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::VariantProjector;

/// Plus-strand single-exon transcript NM_POLY.1 on NC_000001.11 with a poly-A
/// run in the CDS, so a 5'-anchored deletion is non-3'-canonical.
/// CDS = "ATGAAAAACGCTAA" (14 nt): poly-A run c.4..c.8.
fn poly_a_fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_POLY.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("POLYGENE".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1014, 0, 14]],
            cds_start: Some(0),
            cds_end: Some(14),
            gene_id: None,
            protein: Some("NP_POLY.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_POLY.1".to_string(),
        Some("POLYGENE".to_string()),
        TxStrand::Plus,
        "ATGAAAAACGCTAA".to_string(),
        Some(1),
        Some(14),
        vec![Exon::new(1, 1, 14)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1013),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence(
        "NC_000001.11",
        format!("{}ATGAAAAACGCTAA{}", prefix, suffix),
    );

    VariantProjector::new(projector, provider)
}

#[test]
fn project_to_genomic_normalized_3prime_shifts_noncanonical_genome_input() {
    let vp = poly_a_fixture();
    // A non-3'-most g. deletion inside the poly-A run g.1003..1007. The raw pivot
    // is idempotent on Genome input (#785) — it passes the input through
    // un-normalized — whereas the normalized form 3'-shifts to the spec-canonical
    // position (HGVS 3'rule). This is the case #867 is about.
    let v = parse_hgvs("NC_000001.11:g.1003del").unwrap();

    let raw = format!("{}", vp.project_to_genomic(&v).expect("raw pivot"));
    assert_eq!(
        raw, "NC_000001.11:g.1003del",
        "raw project_to_genomic must NOT normalize a Genome input"
    );

    let norm = format!(
        "{}",
        vp.project_to_genomic_normalized(&v).expect("normalized")
    );
    assert_eq!(
        norm, "NC_000001.11:g.1007del",
        "project_to_genomic_normalized must return the 3'-most canonical form"
    );
}

#[test]
fn project_to_genomic_splits_raw_vs_normalized_for_transcript_input() {
    let vp = poly_a_fixture();
    // The transcript-side contract split (the case #867 advertises). The pivot
    // (`project_to_genomic_nc`) does NOT 3'-shift, so a non-3'-most c. deletion in
    // the poly-A run projects to a non-canonical genomic image: c.4del deletes the
    // first A of the run (c.4..c.8 → g.1003..1007) and the raw pivot yields the
    // 5'-anchored g.1003del. `project_to_genomic_normalized` 3'-shifts it to the
    // spec-canonical g.1007del. If `project_to_genomic` re-normalized internally,
    // the raw form would wrongly equal g.1007del — so this catches the contract
    // break that a plus-strand-only "raw == normalized" assertion misses.
    let v = parse_hgvs("NC_000001.11(NM_POLY.1):c.4del").unwrap();
    let raw = format!("{}", vp.project_to_genomic(&v).expect("raw"));
    let norm = format!(
        "{}",
        vp.project_to_genomic_normalized(&v).expect("normalized")
    );
    assert_eq!(
        raw, "NC_000001.11:g.1003del",
        "raw project_to_genomic must NOT 3'-shift the transcript-coordinate projection"
    );
    assert_eq!(
        norm, "NC_000001.11:g.1007del",
        "project_to_genomic_normalized must return the 3'-most canonical form"
    );
    assert_ne!(
        raw, norm,
        "raw and normalized must differ for a non-3'-most transcript input"
    );
}
