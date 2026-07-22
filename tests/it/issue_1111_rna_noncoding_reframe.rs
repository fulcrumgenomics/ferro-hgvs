//! Issue #1111: PR-CI coverage for the cross-axis `r.` → `n.` reframe (#1086
//! D2 / #1089).
//!
//! `VariantProjection::noncoding` is contracted to *always* mean the `n.`
//! form. For an `r.` input the projector must therefore reshape it into a
//! genuine `n.` description — an `n.`-prefixed `Tx` variant with the RNA `U`
//! bases mapped to DNA `T` (`projector.rs`'s `noncoding_axis`). The only test
//! that exercised this reframe (`issue_1086_cross_axis_projection.rs::
//! rna_input_noncoding_axis_is_an_n_description`) `return`s early whenever
//! `FERRO_MANIFEST` is unset, so it runs only in the nightly reference-backed
//! job — a regression to the U→T / prefix logic would pass PR CI and surface
//! only nightly.
//!
//! These tests drive the same reframe with an in-memory [`JsonProvider`] +
//! `CdotMapper` fixture (no manifest, no reference download), so PR CI catches
//! the regression. They mirror the `MockProvider` fixtures in
//! `issue_328_projector_accepts_cnr.rs`.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::JsonProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{HgvsVariant, VariantProjector};

/// A single-exon **non-coding** transcript `NR_TEST.1` on chr1 plus strand.
/// The transcript sequence is `ACGTACGTAC`, so `r.4` is a `u` (DNA `T`) — the
/// base the U→T reframe must translate. Non-coding means the projector takes
/// the `!is_coding` early return, whose sole cross-axis output is the reframed
/// `.noncoding`, isolating the reframe under test from protein prediction.
fn noncoding_fixture() -> (Projector, JsonProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NR_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("NCRNA".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            exons: vec![[2000, 2010, 0, 10]],
            cds_start: None,
            cds_end: None,
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = JsonProvider::new();
    provider.add_transcript(Transcript::new(
        "NR_TEST.1".to_string(),
        Some("NCRNA".to_string()),
        TxStrand::Plus,
        "ACGTACGTAC".to_string(),
        None,
        None,
        vec![Exon::new(1, 1, 10)],
        Some("chr1".to_string()),
        Some(2000),
        Some(2009),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    (projector, provider)
}

/// A single-exon **coding** transcript `NM_TEST.1` on chr1 plus strand,
/// sequence `ATGCGCTAA`, CDS = whole transcript. `r.7` is a `u` (DNA `T`),
/// inside the CDS, so the reframe runs on the *coding* path (which also emits
/// `c.`/`p.`) — the same shape as the manifest-gated `LRG_199t1:r.11u>g` case.
fn coding_fixture() -> (Projector, JsonProvider) {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
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

    let mut provider = JsonProvider::new();
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
    (projector, provider)
}

/// The headline reframe: a bare `r.` input on a non-coding transcript must
/// yield a genuine `n.` description on `.noncoding` — an `n.`-prefixed `Tx`
/// variant (not the `r.` input echoed back), with the RNA `u` reference mapped
/// to DNA `T` (uppercase, `T` not `u`). A regression to the U→T / prefix
/// logic (echoing `NR_TEST.1:r.4u>a`) fails these assertions in PR CI.
#[test]
fn rna_input_reframed_to_n_description_noncoding_transcript() {
    let (projector, provider) = noncoding_fixture();
    let vp = VariantProjector::new(projector, provider);
    let proj = vp
        .project("NR_TEST.1:r.4u>a", "NR_TEST.1")
        .expect("bare r. input on a non-coding transcript should project");

    let noncoding = proj.noncoding.expect("noncoding axis must be present");
    assert!(
        matches!(noncoding, HgvsVariant::Tx(_)),
        "noncoding axis must hold an n. (Tx) variant, got {noncoding}",
    );
    // Genuine n. description: `n.` prefix, DNA `T` (uppercase), not `r.`/`u`.
    assert_eq!(noncoding.to_string(), "NR_TEST.1:n.4T>A");

    // Non-coding transcript: no CDS, so no c./p. forms are synthesized.
    assert!(proj.coding.is_none(), "got coding {:?}", proj.coding);
    assert!(proj.protein.is_none(), "got protein {:?}", proj.protein);
}

/// The same reframe on a *coding* transcript (matches the manifest-gated
/// `LRG_199t1:r.11u>g` shape): `.noncoding` is still a genuine `n.` `Tx`
/// description with U→T applied, and the `c.` axis is produced alongside it.
#[test]
fn rna_input_reframed_to_n_description_coding_transcript() {
    let (projector, provider) = coding_fixture();
    let vp = VariantProjector::new(projector, provider);
    let proj = vp
        .project("NM_TEST.1:r.7u>a", "NM_TEST.1")
        .expect("bare r. input on a coding transcript should project");

    let noncoding = proj.noncoding.expect("noncoding axis must be present");
    assert!(
        matches!(noncoding, HgvsVariant::Tx(_)),
        "noncoding axis must hold an n. (Tx) variant, got {noncoding}",
    );
    assert_eq!(noncoding.to_string(), "NM_TEST.1:n.7T>A");

    // The coding axis is produced alongside the reframed n. form.
    let coding = proj
        .coding
        .expect("coding axis must be present for a coding tx");
    assert_eq!(coding.to_string(), "NM_TEST.1:c.7T>A");
}

/// Regression guard on the reframe: an `n.` input already *is* the `n.` form,
/// so it must pass through `.noncoding` unchanged (the reframe only reshapes
/// the `r.` case). Pins that the reframe doesn't corrupt genuine `n.` inputs.
#[test]
fn noncoding_input_noncoding_axis_unchanged() {
    let (projector, provider) = noncoding_fixture();
    let vp = VariantProjector::new(projector, provider);
    let proj = vp
        .project("NR_TEST.1:n.4T>A", "NR_TEST.1")
        .expect("bare n. input on a non-coding transcript should project");

    assert_eq!(
        proj.noncoding.expect("noncoding axis").to_string(),
        "NR_TEST.1:n.4T>A",
    );
}
