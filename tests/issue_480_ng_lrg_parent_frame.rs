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
//! Without a placement the chromosome coordinate cannot be re-anchored into the
//! parent frame, so `project_to_genomic` declines (`UnsupportedProjection`)
//! rather than stamp a chromosome coordinate under the parent accession, which
//! would be invalid HGVS (#655).
//!
//! Fixture: NM_TEST.1 on chr1 plus strand, "ATGCGCTAA" at genomic [1000, 1009);
//! `c.4C>A` resolves to chromosome position 1003.
//!
//! #646 extends this `NG_`/`LRG_` framing to the *cross-isoform* fan-out: a
//! c./n./r. input with an `NG_` genomic context is de-anchored into the
//! chromosome frame and enumerated against all overlapping transcripts, each
//! re-framed under the parent. That enumeration path is exercised by
//! `project_variant_all_frames_coding_input_with_ng_context` in
//! `src/project/projector.rs` (not in this file).
//!
//! Source (#646): <https://github.com/fulcrumgenomics/ferro-hgvs/issues/646>.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, FerroError, VariantProjector};

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
    // cdot exon genome coords are 1-based HGVS values (exon spans HGVS
    // [1000, 1009)); sequence fetches convert to 0-based by subtracting 1, so
    // HGVS g.1000 is 0-based index 999. Pad with 999 Ns so "ATGCGCTAA" lands at
    // 0-based [999, 1008) and the exon's genome bases match the transcript —
    // otherwise the #644 sequence-aware projection sees a phantom 1-bp indel.
    let prefix = "N".repeat(999);
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
fn ng_parent_without_placement_declines() {
    // No placement registered: cdot carries only the transcript's chromosome
    // (NC_) alignment, so the chromosome coordinate cannot be re-anchored into
    // the NG_ parent's own frame. Emitting the chromosome coordinate under the
    // NG_ accession (e.g. `NG_NOPLACE.1:g.1003C>A`) would be invalid HGVS — the
    // coordinate is a chromosome position, not a position in the NG_ frame — so
    // `project_to_genomic` declines with `UnsupportedProjection` (#655). This
    // mirrors the unit test `project_to_genomic_ng_parent_without_placement_declines`.
    let (projector, provider) = base_fixture();
    let vp = VariantProjector::new(projector, provider);
    let variant = parse_hgvs("NG_NOPLACE.1(NM_TEST.1):c.4C>A").expect("parse");

    let err = vp
        .project_to_genomic(&variant)
        .expect_err("NG_ parent without a placement must decline, not emit invalid HGVS");
    match err {
        FerroError::UnsupportedProjection { reason } => assert!(
            reason.contains("NG_NOPLACE.1") && reason.contains("no chromosomal placement"),
            "expected a no-placement decline naming NG_NOPLACE.1, got: {reason}"
        ),
        other => panic!("expected UnsupportedProjection, got: {other:?}"),
    }
}

fn parse_acc(s: &str) -> ferro_hgvs::hgvs::variant::Accession {
    match parse_hgvs(&format!("{s}:g.1=")).expect("parse accession") {
        ferro_hgvs::HgvsVariant::Genome(g) => g.accession,
        other => panic!("expected genome variant, got {other:?}"),
    }
}

/// A bare LRG transcript input (`LRG_<n>t<m>:c.…`, no genomic-context
/// parenthetical — the form the conformance corpus uses) must project to its
/// genomic LRG (`LRG_<n>`) in the LRG's own frame: the genomic parent is derived
/// structurally from the transcript accession, cdot resolves the LRG transcript
/// to its RefSeq NM_ for the chromosome mapping, and the LRG placement
/// re-anchors the coordinate.
#[test]
fn bare_lrg_transcript_input_reanchors_to_lrg_frame() {
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
    // Map LRG_1t1 -> NM_TEST.1 via a temp mapping file (NCBI/LRG column layout:
    // LRG, HGNC_SYMBOL, REFSEQ_GENOMIC, LRG_TRANSCRIPT, REFSEQ_TRANSCRIPT, ...).
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
    // 999-N pad so HGVS exon coord 1000 lands at 0-based index 999 (see the
    // matching note in `base_fixture`): the exon's genome bases then match the
    // transcript and the #644 sequence-aware projection sees no phantom indel.
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
    let vp = VariantProjector::new(projector, provider);

    // Bare LRG transcript input — no genomic-context parenthetical.
    let variant = parse_hgvs("LRG_1t1:c.4C>A").expect("parse");
    let g = vp
        .project_to_genomic(&variant)
        .expect("bare LRG transcript should project to its genomic LRG");
    assert_eq!(g.to_string(), "LRG_1:g.4C>A");
}
