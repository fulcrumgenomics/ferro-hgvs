//! #2280: a predicted cis allele `[(a;b)]` stays predicted on every projected
//! axis. The `(…)` marks a change as predicted rather than observed
//! (`recommendations/general.md:12`, and `:85` for the parentheses), and the
//! wrapper belongs to the whole allele.
//!
//! Before this change `project_allele_inner` rebuilt the `c.`, `n.` and `g.`
//! alleles with `AlleleVariant::new`, which is always certain. A predicted
//! `g.[(1003C>A;1007T>C)]` projected to `c.[4C>A;8T>C]`, `n.[4C>A;8T>C]` and
//! `g.[1003C>A;1007T>C]` beside an `r.[(4c>a;8u>c)]` that kept the wrapper.
//! Where the codon-frame exception then merged the coding allele, the merge
//! came out certain too: `g.(1003_1005delinsAGA)` gave `c.4_6delinsAGA`.
//!
//! #2276 (single edits on `c.`/`n.`) uses the same fixture: one 9-base
//! transcript, `ATGCGCTAA`, on `NC_000001.11:1000-1008`, so `c.N` is `n.N` is
//! `g.(999+N)`. Codon 2 is `CGC` at `c.4_6`.

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::VariantProjector;

fn fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TX1.1".to_string(),
        CdotTranscript {
            cds_start_incomplete: false,
            gene_name: Some("GENE1".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: None,
            exon_cigars: Vec::new(),
        },
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_TX1.1".to_string(),
        Some("GENE1".to_string()),
        TxStrand::Plus,
        "ATGCGCTAA".to_string(),
        Some(1),
        Some(9),
        vec![Exon::new(1, 1, 9)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1008),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    provider.add_genomic_sequence(
        "NC_000001.11",
        format!("{}ATGCGCTAA{}", "N".repeat(999), "N".repeat(100)),
    );
    VariantProjector::new(Projector::new(cdot), provider)
}

/// The `(coding, noncoding, genomic)` axes of projecting `input` onto `NM_TX1.1`.
fn axes(input: &str) -> (Option<String>, Option<String>, Option<String>) {
    let v = ferro_hgvs::parse_hgvs(input).unwrap_or_else(|e| panic!("{input} parses: {e}"));
    let p = fixture()
        .project_variant(&v, "NM_TX1.1")
        .unwrap_or_else(|e| panic!("{input} projects: {e}"));
    (
        p.coding.map(|v| v.to_string()),
        p.noncoding.map(|v| v.to_string()),
        p.genomic.map(|v| v.to_string()),
    )
}

fn some(s: &str) -> Option<String> {
    Some(s.to_string())
}

/// A predicted genomic allele the codon exception does not touch: the wrapper
/// rides onto both transcript axes and stays on the genomic one.
#[test]
fn a_predicted_genomic_allele_stays_predicted_on_every_axis() {
    assert_eq!(
        axes("NC_000001.11:g.[(1003C>A;1007T>C)]"),
        (
            some("NC_000001.11(NM_TX1.1):c.[(4C>A;8T>C)]"),
            some("NC_000001.11(NM_TX1.1):n.[(4C>A;8T>C)]"),
            some("NC_000001.11:g.[(1003C>A;1007T>C)]"),
        ),
    );
}

/// A parent-framed predicted `c.` allele goes through the genome pivot member
/// by member; all three reassembled alleles keep the wrapper.
#[test]
fn a_predicted_parent_framed_coding_allele_stays_predicted_on_every_axis() {
    assert_eq!(
        axes("NC_000001.11(NM_TX1.1):c.[(4C>A;8A>C)]"),
        (
            some("NC_000001.11(NM_TX1.1):c.[(4C>A;8A>C)]"),
            some("NC_000001.11(NM_TX1.1):n.[(4C>A;8A>C)]"),
            some("NC_000001.11:g.[(1003C>A;1007A>C)]"),
        ),
    );
}

/// A bare `NM_` predicted allele, with no genomic axis.
#[test]
fn a_predicted_bare_transcript_allele_stays_predicted() {
    assert_eq!(
        axes("NM_TX1.1:c.[(4C>A;8A>C)]"),
        (
            some("NM_TX1.1:c.[(4C>A;8A>C)]"),
            some("NM_TX1.1:n.[(4C>A;8A>C)]"),
            None,
        ),
    );
}

/// The reported case. The genomic input normalizes to the predicted allele
/// `g.[(1003C>A;1005C>A)]`; the codon-frame exception merges the coding axis
/// to one `delins`, which takes the wrapper on its edit. The `g.` axis is the
/// input's own normalization and stays split (#79).
#[test]
fn a_predicted_codon_change_merged_on_the_coding_axis_stays_predicted() {
    let expected = (
        some("NC_000001.11(NM_TX1.1):c.(4_6delinsAGA)"),
        some("NC_000001.11(NM_TX1.1):n.(4_6delinsAGA)"),
        some("NC_000001.11:g.[(1003C>A;1005C>A)]"),
    );
    assert_eq!(axes("NC_000001.11:g.(1003_1005delinsAGA)"), expected);
    assert_eq!(axes("NC_000001.11:g.[(1003C>A;1005C>A)]"), expected);
}

/// A merge that leaves more than one member: the remaining allele carries the
/// wrapper, not the merged member.
#[test]
fn a_predicted_allele_partly_merged_on_the_coding_axis_stays_predicted() {
    assert_eq!(
        axes("NC_000001.11:g.[(1003C>A;1005C>A;1007A>C)]"),
        (
            some("NC_000001.11(NM_TX1.1):c.[(4_6delinsAGA;8A>C)]"),
            some("NC_000001.11(NM_TX1.1):n.[(4_6delinsAGA;8A>C)]"),
            some("NC_000001.11:g.[(1003C>A;1005C>A;1007A>C)]"),
        ),
    );
}

/// The observed controls: the same changes without the wrapper gain none on
/// any axis, and take the same partition as their predicted forms above.
#[test]
fn observed_alleles_gain_no_wrapper() {
    assert_eq!(
        axes("NC_000001.11:g.[1003C>A;1007T>C]"),
        (
            some("NC_000001.11(NM_TX1.1):c.[4C>A;8T>C]"),
            some("NC_000001.11(NM_TX1.1):n.[4C>A;8T>C]"),
            some("NC_000001.11:g.[1003C>A;1007T>C]"),
        ),
    );
    assert_eq!(
        axes("NC_000001.11:g.[1003C>A;1005C>A]"),
        (
            some("NC_000001.11(NM_TX1.1):c.4_6delinsAGA"),
            some("NC_000001.11(NM_TX1.1):n.4_6delinsAGA"),
            some("NC_000001.11:g.[1003C>A;1005C>A]"),
        ),
    );
    assert_eq!(
        axes("NC_000001.11:g.[1003C>A;1005C>A;1007A>C]"),
        (
            some("NC_000001.11(NM_TX1.1):c.[4_6delinsAGA;8A>C]"),
            some("NC_000001.11(NM_TX1.1):n.[4_6delinsAGA;8A>C]"),
            some("NC_000001.11:g.[1003C>A;1005C>A;1007A>C]"),
        ),
    );
}

/// The `r.` axis of projecting `input` onto `NM_TX1.1`.
fn rna(input: &str) -> Option<String> {
    let v = ferro_hgvs::parse_hgvs(input).unwrap_or_else(|e| panic!("{input} parses: {e}"));
    fixture()
        .project_variant(&v, "NM_TX1.1")
        .unwrap_or_else(|e| panic!("{input} projects: {e}"))
        .rna
        .map(|v| v.to_string())
}

/// A one-member allele whose member is predicted keeps the wrapper on `r.` as
/// well. `predict_rna` stripped each member's `( )` for the allele's `[( )]` to
/// replace, and a one-member allele renders as its bare member, so nothing put
/// it back: `r.4c>a` beside `c.(4C>A)`.
#[test]
fn a_one_member_predicted_allele_stays_predicted_on_the_rna_axis() {
    assert_eq!(
        axes("NC_000001.11:g.[(1003C>A)]"),
        (
            some("NC_000001.11(NM_TX1.1):c.(4C>A)"),
            some("NC_000001.11(NM_TX1.1):n.(4C>A)"),
            some("NC_000001.11:g.(1003C>A)"),
        ),
    );
    assert_eq!(rna("NC_000001.11:g.[(1003C>A)]"), some("NM_TX1.1:r.(4c>a)"));
    assert_eq!(
        rna("NC_000001.11(NM_TX1.1):c.[(4C>A)]"),
        some("NC_000001.11(NM_TX1.1):r.(4c>a)")
    );
}

/// The observed control. An observed one-member allele gains no wrapper on
/// `c.`, `n.` or `g.`, and its `r.` is the same predicted `r.(4c>a)` as the
/// change written as one variant: the projector's `r.` is a prediction from
/// the `c.` axis whatever the input states, and the brackets do not change
/// that. Multi-member alleles already rendered `r.[(…)]` and are unchanged.
#[test]
fn an_observed_one_member_allele_gains_no_wrapper_and_agrees_on_the_rna_axis() {
    assert_eq!(
        axes("NC_000001.11:g.[1003C>A]"),
        (
            some("NC_000001.11(NM_TX1.1):c.4C>A"),
            some("NC_000001.11(NM_TX1.1):n.4C>A"),
            some("NC_000001.11:g.1003C>A"),
        ),
    );
    assert_eq!(
        rna("NC_000001.11:g.[1003C>A]"),
        rna("NC_000001.11:g.1003C>A")
    );
    assert_eq!(rna("NC_000001.11:g.1003C>A"), some("NM_TX1.1:r.(4c>a)"));
    for input in [
        "NC_000001.11:g.[1003C>A;1007T>C]",
        "NC_000001.11:g.[(1003C>A;1007T>C)]",
    ] {
        assert_eq!(rna(input), some("NM_TX1.1:r.[(4c>a;8u>c)]"), "{input}");
    }
}
