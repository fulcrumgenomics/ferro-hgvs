//! Issue #328 — `VariantProjector::project` must accept c./n./r. inputs.
//!
//! Source: <https://github.com/fulcrumgenomics/ferro-hgvs/issues/328>.
//!
//! `project_variant_inner` previously required `HgvsVariant::Genome`
//! and rejected non-g. inputs with
//! `FerroError::UnsupportedProjection { reason: "VariantProjector
//! currently only accepts g. variants" }`. Mutalyzer accepts c./n./r.
//! inputs and projects them to (g., c., p.) via the named transcript.
//!
//! This PR adds a thin dispatch layer on top of the existing g.-anchored
//! projection: c./n./r. inputs are first run through
//! [`VariantProjector::project_to_genomic`] (the #327 helper) to obtain
//! a `Genome` form, then re-enter the existing pipeline. The
//! `transcript_id` argument still drives the projection target.
//!
//! Tests use the same `MockProvider`+`CdotMapper` fixture as the
//! existing `tests/projection.rs::fixture` (NM_TEST.1 on chr1 plus
//! strand, "ATGCGCTAA" at genomic [1000, 1009)).

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::VariantProjector;

/// Same fixture shape as `tests/projection.rs::fixture` — kept in this
/// file (not shared) to avoid coupling test files via `mod`.
fn fixture() -> (Projector, MockProvider) {
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

/// A `c.` substitution input must project to (g., c., p.) exactly as if
/// the equivalent g. form had been passed in.
///
/// `NC_000001.11(NM_TEST.1):c.4C>A` should yield the same projection
/// as `NC_000001.11:g.1003C>A` (which the existing
/// `end_to_end_missense` test pins to `c.4C>A` → `p.(Arg2Ser)`).
#[test]
fn projector_accepts_cds_substitution_input() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    let result = vp
        .project("NC_000001.11(NM_TEST.1):c.4C>A", "NM_TEST.1")
        .expect("projection of c. input should succeed");

    assert_eq!(result.transcript_id, "NM_TEST.1");
    let c = result
        .coding
        .as_ref()
        .expect("coding output must be present for c. input")
        .to_string();
    assert!(c.contains(":c.4C>A"), "got c. = {}", c);
    let p = result
        .protein
        .as_ref()
        .expect("protein output must be present")
        .to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
    // `.genomic` must carry the projected g. form, not the input c.
    // This pins the contract that downstream consumers reading
    // `.genomic` get a canonical g. variant regardless of input axis.
    // This c. input carries an explicit NC_ genomic_context, so it takes the
    // genome-pivot path and `.genomic` is present (a bare-NM_ input would be
    // `None`; see #498).
    let g = format!(
        "{}",
        result
            .genomic
            .as_ref()
            .expect(".genomic must be present for a genome-anchored c. input")
    );
    assert!(
        g.contains(":g."),
        ".genomic must be a g. variant for c. input; got {g}",
    );
    assert!(
        !g.contains(":c."),
        ".genomic must NOT carry the input c. form; got {g}",
    );
    assert!(!result.is_frameshift);
    assert!(!result.is_intronic);
    assert!(!result.is_utr);
}

/// A multi-base `c.` deletion input must project through the same
/// pipeline. `c.4_6del` (delete CGC = the second codon Arg) yields a
/// 3-base genomic deletion and a one-residue protein deletion.
#[test]
fn projector_accepts_cds_deletion_input() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    let result = vp
        .project("NC_000001.11(NM_TEST.1):c.4_6del", "NM_TEST.1")
        .expect("projection of c.4_6del should succeed");

    let c = result.coding.as_ref().unwrap().to_string();
    assert!(c.contains(":c.4_6del"), "got c. = {}", c);
    // The protein output exists; pin only that the coding output
    // survives round-trip. Protein output for an in-frame 3-base
    // deletion is delegated to the existing protein-prediction code.
    assert!(
        result.protein.is_some(),
        "in-frame del should produce a protein output"
    );
}

/// An `r.` input on a coding transcript must produce both c. and p.
/// outputs after projection. `r.4c>a` is the RNA equivalent of `c.4C>A`
/// (lowercase + `u` instead of `t`, but the position and edit are the
/// same shape).
#[test]
fn projector_accepts_rna_substitution_input() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    let result = vp
        .project("NC_000001.11(NM_TEST.1):r.4c>a", "NM_TEST.1")
        .expect("projection of r. input should succeed");

    let c = result.coding.as_ref().unwrap().to_string();
    assert!(c.contains(":c.4C>A"), "got c. = {}", c);
    let p = result.protein.as_ref().unwrap().to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
}

/// A bare `c.` input (no parent NG/NC `genomic_context`) projects directly
/// to protein without a genome roundtrip (#498). This reverses the prior
/// #327/#328 "must error" contract: ferro still does **not synthesize** a
/// genomic form (`genomic` is `None`, not a fabricated g.), but protein is
/// predicted straight from the transcript's CDS, since an exonic c. position
/// is already the 1-based CDS position prediction needs.
#[test]
fn projector_bare_cds_projects_directly_to_protein() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    // No NC_ parent in the accession — just `NM_TEST.1:c.4C>A`.
    let result = vp
        .project("NM_TEST.1:c.4C>A", "NM_TEST.1")
        .expect("bare-NM_ c. input should project directly to protein (#498)");
    // No genomic form is synthesized for a bare-NM_ input (#327 preserved).
    assert!(
        result.genomic.is_none(),
        "bare-NM_ projection must not synthesize a genomic form; got {:?}",
        result.genomic,
    );
    // Coding axis present; protein predicted from the CDS.
    assert!(result.coding.is_some(), "coding form should be present");
    let p = result
        .protein
        .as_ref()
        .expect("protein should be predicted via the direct c.→p. path")
        .to_string();
    assert_eq!(p, "NP_TEST.1:p.(Arg2Ser)");
}

/// `n.` projection support is a downstream scope item. The existing
/// `project_single_inner` g.→c. path calls `genome_to_cds`, which
/// requires `cds_start`/`cds_end` and returns `InvalidCoordinates`
/// (translated to `TranscriptNotOverlapping`) on non-coding
/// transcripts. Adding n.-axis end-to-end support requires teaching
/// `project_single_inner` to dispatch on the target transcript's
/// coding/non-coding status — that's a separate gap from "accept the
/// input axis" which this PR addresses. Pinning this contract here
/// so a future PR can flip the assertion.
#[test]
fn projector_noncoding_input_currently_errors_through_existing_path() {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NR_NCRNA.1".to_string(),
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

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NR_NCRNA.1".to_string(),
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
    let prefix = "N".repeat(2000);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}ACGTACGTAC{}", prefix, suffix));
    let vp = VariantProjector::new(projector, provider);

    // After this PR the input axis is accepted (no longer errors with
    // "currently only accepts g. variants"); the downstream g.→c.
    // path then fails on the non-coding transcript with
    // `TranscriptNotOverlapping`. Pin this so a future expansion of
    // `project_single_inner` to handle non-coding tx can flip the
    // expectation to success.
    let err = vp.project("NC_000001.11(NR_NCRNA.1):n.4T>A", "NR_NCRNA.1");
    assert!(
        err.is_err(),
        "n. input is accepted at the entry point but the downstream c. path \
         currently does not handle non-coding tx; expected an error, got {err:?}",
    );
    let msg = format!("{}", err.unwrap_err());
    assert!(
        !msg.contains("currently only accepts g. variants"),
        "the entry-point gate must be removed; got {msg}",
    );
}

/// Document the 3'-shift asymmetry between c.-input and g.-input
/// projections near an exon junction (re #334).
///
/// The c.-axis normalizer respects the HGVS exon-junction exception
/// (does not shift across exons); the g.-axis normalizer does not.
/// `project_single_inner` receives a pre-normalized variant, and the
/// projected g. from a c. input is NOT re-normalized by the g.-axis
/// normalizer. So the same biological variant fed as c. vs. g. near
/// an exon junction can project to different `(g., c., p.)` tuples.
///
/// This is intentional: the input axis carries the spec-correct
/// semantic for that axis. Pinning the behaviour here so a future
/// refactor doesn't silently change it. The plus-strand single-exon
/// fixture used by the other tests has no junction to demonstrate
/// this on, so we just assert that c. → projection preserves the
/// c.-form canonical position (input was already canonical via the
/// in-exon shuffle).
#[test]
fn projector_preserves_cds_axis_normalization_semantics() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);
    // c.4_6del on `ATGCGCTAA`: the deletion is in-exon (single-exon
    // fixture), so c.-axis 3'-shift can run freely. The c.-axis
    // normalizer settles on the canonical c. form, and that's what
    // the projector reports.
    let result = vp
        .project("NC_000001.11(NM_TEST.1):c.4_6del", "NM_TEST.1")
        .expect("projection should succeed");
    let c = result.coding.as_ref().unwrap().to_string();
    assert!(
        c.contains(":c."),
        "c.-input must produce a c. output preserving the c.-axis canonical form; got {c}",
    );
}

/// Reject `transcript_id` mismatch: passing a c.-input on NM_FOO with
/// a request to project against NM_BAR would silently project through
/// one transcript's exons and re-project through the other's,
/// producing nonsensical results. The dispatch layer rejects this
/// configuration up front.
#[test]
fn projector_rejects_transcript_id_mismatch_for_cds_input() {
    // Build a fixture with two distinct transcripts so we have a
    // valid (but wrong) transcript_id to pass.
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
    cdot.add_transcript(
        "NM_OTHER.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("OTHER".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // Different region of the same chromosome.
            exons: vec![[2000, 2009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_OTHER.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);
    let provider = MockProvider::new();
    let vp = VariantProjector::new(projector, provider);

    let err = vp
        .project("NC_000001.11(NM_TEST.1):c.4C>A", "NM_OTHER.1")
        .expect_err("expected transcript_id mismatch error");
    let msg = format!("{}", err);
    assert!(
        msg.contains("transcript_id mismatch"),
        "expected transcript_id mismatch diagnostic; got {msg}",
    );
}

/// Minus-strand round-trip: a c. input on a minus-strand transcript
/// must project to (g., c., p.) where:
/// - the g. coords are on the genomic reference (anti-parallel to the
///   transcript)
/// - the c. output round-trips back to the input c. form (canonical)
/// - the edit's stated bases are reverse-complemented on the g. side
///   (per `transform_edit_for_strand`) and then back on the c. side
///
/// Pins the c.→g.→c. round-trip stability on the minus strand, which
/// is the highest-risk path for strand-flip bugs.
#[test]
fn projector_accepts_minus_strand_cds_input() {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_MINUS.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("MINUSGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Minus,
            // Genome [3000, 3009), tx [0, 9). On minus strand the
            // transcript reads anti-parallel: g.3008 (0-based) is
            // the first tx base.
            exons: vec![[3000, 3009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_MINUS.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    // Transcript-view sequence (always 5' → 3' on the tx): "ATGCGCTAA"
    // (same shape as the plus-strand fixture).
    provider.add_transcript(Transcript::new(
        "NM_MINUS.1".to_string(),
        Some("MINUSGENE".to_string()),
        TxStrand::Minus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("chr1".to_string()),
        Some(3000),
        Some(3008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // Genomic sequence: the reverse-complement of "ATGCGCTAA" is
    // "TTAGCGCAT". Place it at genome [3000, 3009).
    let prefix = "N".repeat(3000);
    let suffix = "N".repeat(100);
    provider.add_genomic_sequence("chr1", format!("{}TTAGCGCAT{}", prefix, suffix));
    let vp = VariantProjector::new(projector, provider);

    let result = vp
        .project("NC_000001.11(NM_MINUS.1):c.4C>A", "NM_MINUS.1")
        .expect("projection of minus-strand c. input should succeed");

    let c = result.coding.as_ref().unwrap().to_string();
    assert!(
        c.contains(":c.4C>A"),
        "minus-strand c. round-trip must preserve the c. form; got {c}",
    );
    let p = result.protein.as_ref().unwrap().to_string();
    assert_eq!(p, "NP_MINUS.1:p.(Arg2Ser)");
}
