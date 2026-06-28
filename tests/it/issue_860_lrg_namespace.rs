//! #860: preserve the input's LRG_t/LRG_p namespace on projection output.
//!
//! Case A — an `LRG_<n>t<k>:c.` transcript-coordinate input: the coding axis
//! already passes through as `LRG_<n>t<k>`, but the protein axis must echo the
//! input's LRG namespace (`LRG_<n>p<k>`) rather than the resolved RefSeq `NP_`.
//! Both forms are spec-valid public references (general.md L16); preserving the
//! input's namespace converges with mutalyzer.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/860>.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, VariantProjector};

fn parse_acc(s: &str) -> ferro_hgvs::hgvs::variant::Accession {
    match parse_hgvs(&format!("{s}:g.1=")).expect("parse accession") {
        ferro_hgvs::HgvsVariant::Genome(g) => g.accession,
        other => panic!("expected genome variant, got {other:?}"),
    }
}

/// Build a projector + provider where `LRG_1t1` resolves to `NM_TEST.1`
/// (protein `NP_TEST.1`), with an `LRG_1` genomic placement. Transcript
/// `ATGCGCTAA` = Met-Arg-Stop; CDS is the whole 9 bp.
fn lrg_fixture() -> (Projector, MockProvider) {
    use std::io::Write;

    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    // Map LRG_1t1 -> NM_TEST.1 (NCBI/LRG column layout).
    let mut mapping = tempfile::NamedTempFile::new().unwrap();
    writeln!(mapping, "# LRG\tHGNC\tREFSEQ_GENOMIC\tLRG_TX\tREFSEQ_TX").unwrap();
    writeln!(
        mapping,
        "LRG_1\tTESTGENE\tNG_000001.1\tt1\tNM_TEST.1\tENST0\tCCDS0"
    )
    .unwrap();
    cdot.load_lrg_mapping(mapping.path())
        .expect("load lrg mapping");
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    // The real MultiFastaProvider resolves `LRG_1t1` to its RefSeq NM_ FASTA via
    // the LRG mapping; MockProvider has no such aliasing, so register the
    // transcript under both ids (same sequence/CDS) — coordinate mapping resolves
    // LRG_1t1 -> NM_TEST.1 through cdot, while protein prediction reads the FASTA
    // the provider returns for the requested id.
    for tx_id in ["NM_TEST.1", "LRG_1t1"] {
        provider.add_transcript(Transcript::new(
            tx_id.to_string(),
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
    }
    let prefix = "N".repeat(999);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    provider.add_genomic_placement(
        "LRG_1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1009,
            strand: Strand::Plus,
        },
    );
    (projector, provider)
}

/// Case A: `LRG_1t1:c.4C>A` (codon 2 CGC->AGC, Arg->Ser) projects with the
/// protein axis named in the input's LRG namespace (`LRG_1p1`), while the
/// coding axis stays `LRG_1t1`.
#[test]
fn lrg_transcript_input_emits_lrg_protein_accession() {
    let (projector, provider) = lrg_fixture();
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("LRG_1t1:c.4C>A").expect("parse");

    // The conformance protein axis uses `project_variant(&v, tx_id)` with the
    // structurally-inferred transcript id (the bare LRG transcript accession).
    let result = vp.project_variant(&variant, "LRG_1t1").expect("projects");
    let protein = result.protein.as_ref().expect("a protein axis");
    assert!(
        protein.to_string().starts_with("LRG_1p1:"),
        "protein accession should echo the input LRG namespace (LRG_1p1), got: {}",
        protein
    );
    let coding = result.coding.as_ref().expect("a coding axis");
    assert!(
        coding.to_string().starts_with("LRG_1t1:"),
        "coding accession should stay LRG_1t1, got: {}",
        coding
    );
}
