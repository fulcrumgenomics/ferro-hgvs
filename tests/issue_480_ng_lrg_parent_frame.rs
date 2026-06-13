//! Issue #480 — `project_to_genomic` must express transcript coordinates in an
//! `NG_`/`LRG_` parent's own frame, not stamp chromosome coordinates under the
//! parent accession.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/480>.
//!
//! cdot aligns transcripts only to chromosomes (`NM_`→`NC_`). When a c./n./r.
//! input is parented on an `NG_` RefSeqGene or an `LRG_`, the resolved
//! coordinate is on the chromosome but the emitted accession is the
//! `NG_`/`LRG_` one — the two frames disagree. Given the parent's chromosomal
//! placement (`GenomicPlacement`), ferro composes NM→NC (cdot) with the affine
//! NC→parent transform and emits the coordinate in the parent's own frame.
//! Without a placement it keeps the chromosome coordinate (prior behavior),
//! never failing the projection.
//!
//! Fixture: NM_TEST.1 on chr1 plus strand, "ATGCGCTAA" at genomic [1000, 1009);
//! `c.4C>A` resolves to chromosome position 1003.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, VariantProjector};

fn base_fixture() -> (Projector, MockProvider) {
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
    let prefix = "N".repeat(1000);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ATGCGCTAA{}", prefix, suffix));
    (projector, provider)
}

// ---------------------------------------------------------------------------
// Pure affine transform
// ---------------------------------------------------------------------------

#[test]
fn placement_maps_chromosome_to_parent_frame_plus_strand() {
    // Parent occupies chr1 [1000, 1009], its base 1 at chr 1000, same orientation.
    let p = GenomicPlacement {
        nc: parse_acc("NC_000001.11"),
        parent_start: 1,
        nc_start: 1000,
        nc_end: 1009,
        strand: Strand::Plus,
    };
    assert_eq!(p.nc_to_parent(1000), Some(1));
    assert_eq!(p.nc_to_parent(1003), Some(4));
    assert_eq!(p.nc_to_parent(1009), Some(10));
    // Outside the placed span declines rather than emitting an out-of-frame pos.
    assert_eq!(p.nc_to_parent(999), None);
    assert_eq!(p.nc_to_parent(1010), None);
}

#[test]
fn placement_maps_chromosome_to_parent_frame_minus_strand() {
    // Parent antiparallel to the chromosome: base 1 at the high chr coordinate.
    let p = GenomicPlacement {
        nc: parse_acc("NC_000001.11"),
        parent_start: 1,
        nc_start: 1000,
        nc_end: 1009,
        strand: Strand::Minus,
    };
    assert_eq!(p.nc_to_parent(1009), Some(1));
    assert_eq!(p.nc_to_parent(1003), Some(7));
    assert_eq!(p.nc_to_parent(1000), Some(10));
}

// ---------------------------------------------------------------------------
// Projection
// ---------------------------------------------------------------------------

#[test]
fn ng_parent_emits_parent_frame_coords_plus_strand() {
    let (projector, mut provider) = base_fixture();
    provider.add_genomic_placement(
        "NG_TEST.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1009,
            strand: Strand::Plus,
        },
    );
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("NG_TEST.1(NM_TEST.1):c.4C>A").expect("parse");

    let g = vp.project_to_genomic(&variant).expect("should project");
    // chr 1003 → parent frame 4 (1003 - 1000 + 1), edit unchanged on plus strand.
    assert_eq!(g.to_string(), "NG_TEST.1:g.4C>A");
}

#[test]
fn ng_parent_emits_parent_frame_coords_minus_strand() {
    let (projector, mut provider) = base_fixture();
    provider.add_genomic_placement(
        "NG_REV.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 1,
            nc_start: 1000,
            nc_end: 1009,
            strand: Strand::Minus,
        },
    );
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("NG_REV.1(NM_TEST.1):c.4C>A").expect("parse");

    let g = vp.project_to_genomic(&variant).expect("should project");
    // chr 1003 → parent frame 7 (1009 - 1003 + 1); the C>A edit reverse-
    // complements to G>T because the parent runs antiparallel to the chromosome.
    assert_eq!(g.to_string(), "NG_REV.1:g.7G>T");
}

#[test]
fn ng_parent_without_placement_keeps_prior_behavior() {
    // No placement registered: the projection keeps the chromosome coordinate
    // under the parent accession (the pre-#480 behavior). This is the
    // no-regression fallback for NG_ records whose placement is not ingested.
    let (projector, provider) = base_fixture();
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("NG_NOPLACE.1(NM_TEST.1):c.4C>A").expect("parse");

    let g = vp.project_to_genomic(&variant).expect("should project");
    assert_eq!(g.to_string(), "NG_NOPLACE.1:g.1003C>A");
}

fn parse_acc(s: &str) -> ferro_hgvs::hgvs::variant::Accession {
    match parse_hgvs(&format!("{s}:g.1=")).expect("parse accession") {
        ferro_hgvs::HgvsVariant::Genome(g) => g.accession,
        other => panic!("expected genome variant, got {other:?}"),
    }
}
