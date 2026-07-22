//! Issue #1086 (defect D3, genomic-input half) — projecting a genomic (`g.`)
//! input onto the `c.`/`n.` axis must retain the genomic accession as the
//! reference context (`<contig>(NM_…):c.…`), not emit a bare `NM_…:c.…`.
//!
//! `refseq.md:47-49` forbids a bare *intronic* `c.` (it is the spec's own
//! literal "not correct" example), and `refseq.md:138-140` endorses the
//! parenthesised form whenever a `c.` is based on a genomic reference sequence.
//! PR #1090 fixed the transcript-input half (retain an *explicit* context); a
//! bare genomic input has no explicit context, so it used to strip to the bare
//! form. It now stamps its own accession instead.
//!
//! PR-CI-runnable: MockProvider single-exon fixture, no reference manifest.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::VariantProjector;

/// NM_TEST.1 on `chr1` (+ strand), CDS = the whole 9-base exon `ATGCGCTAA` at
/// genomic `[1000, 1009)`. Same shape as `issue_328`'s fixture.
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

/// Round-trip: derive the bare genomic form of an exonic variant, then project
/// that genomic input back onto the `c.` axis. The genomic contig must be
/// retained as the reference context, not stripped to a bare `NM_…:c.…`.
#[test]
fn genomic_input_retains_context_on_c_axis() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);

    let g = vp
        .project("NC_000001.11(NM_TEST.1):c.4C>A", "NM_TEST.1")
        .expect("project to genomic")
        .genomic
        .expect("genomic axis")
        .to_string();
    // The genomic axis is bare (no transcript context of its own).
    assert!(!g.contains('('), "genomic form should be bare: {g}");

    let c = vp
        .project(&g, "NM_TEST.1")
        .expect("project genomic input to c.")
        .coding
        .expect("coding axis")
        .to_string();
    assert_eq!(
        c, "NC_000001.11(NM_TEST.1):c.4C>A",
        "a genomic input must retain its accession as the c.-axis context"
    );
}

/// The non-coding axis retains the context too.
#[test]
fn genomic_input_retains_context_on_n_axis() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);

    let g = vp
        .project("NC_000001.11(NM_TEST.1):c.4C>A", "NM_TEST.1")
        .expect("project to genomic")
        .genomic
        .expect("genomic axis")
        .to_string();

    let n = vp
        .project(&g, "NM_TEST.1")
        .expect("project genomic input to n.")
        .noncoding
        .expect("noncoding axis")
        .to_string();
    assert_eq!(
        n, "NC_000001.11(NM_TEST.1):n.4C>A",
        "a genomic input must retain its accession as the n.-axis context"
    );
}

/// A transcript-context allele input keeps its EXPLICIT parent on every member
/// through the per-member projection recursion. This exercises the #1090
/// (transcript-input) half — the input already carries `NC_(NM_)` framing — not
/// this PR's genomic arm; `genomic_allele_input_retains_context_on_members`
/// below covers the genomic case.
#[test]
fn allele_input_retains_context_on_members() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);

    let c = vp
        .project("NC_000001.11(NM_TEST.1):c.[4C>A];[6C>A]", "NM_TEST.1")
        .expect("project allele")
        .coding
        .expect("coding axis")
        .to_string();
    assert_eq!(
        c, "NC_000001.11(NM_TEST.1):c.[4C>A];[6C>A]",
        "every allele member must keep the reference context"
    );
}

/// A bare **genomic** allele retains context on every member: because
/// `project_allele_inner` projects each member through the single-variant path,
/// the #1086 genomic arm fires per-member, so each stamps its own accession.
/// This pins this PR's behavior for the allele case (the transcript-context
/// test above passes even without the fix; this one does not).
#[test]
fn genomic_allele_input_retains_context_on_members() {
    let (projector, provider) = fixture();
    let vp = VariantProjector::new(projector, provider);

    // Derive the bare genomic form of two exonic members, then assemble a bare
    // genomic allele (no explicit parent) from them.
    let member_g = |input: &str| -> String {
        let g = vp
            .project(input, "NM_TEST.1")
            .expect("project to genomic")
            .genomic
            .expect("genomic axis")
            .to_string();
        // strip the accession prefix, keeping the `g.<edit>` payload
        g.split("g.").nth(1).expect("g. payload").to_string()
    };
    let g1 = member_g("NC_000001.11(NM_TEST.1):c.4C>A");
    let g2 = member_g("NC_000001.11(NM_TEST.1):c.6C>A");
    let allele_input = format!("NC_000001.11:g.[{g1}];[{g2}]");
    assert!(
        !allele_input.contains("(NM"),
        "the genomic allele input must be bare (no explicit parent): {allele_input}"
    );

    let c = vp
        .project(&allele_input, "NM_TEST.1")
        .expect("project genomic allele to c.")
        .coding
        .expect("coding axis")
        .to_string();
    assert_eq!(
        c, "NC_000001.11(NM_TEST.1):c.[4C>A];[6C>A]",
        "each genomic allele member must retain its accession as the c.-axis context"
    );
}
