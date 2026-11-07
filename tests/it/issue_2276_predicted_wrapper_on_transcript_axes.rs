//! #2276: a predicted (uncertain) variant stays predicted on every projected
//! axis. The `(…)` wrapper says the change was "theoretically deduced" rather
//! than "experimentally determined" (`recommendations/general.md:12`, and `:85`
//! for the parentheses), and a single projection result must not say both.
//!
//! Before this change the `n.` axis was built from the bare edit on every path,
//! and the `c.` axis was too whenever projection went through the genome pivot,
//! so `NC_000001.11(NM_TX1.1):c.(5G>T)` projected to `c.5G>T` and `n.5G>T`
//! beside a `g.(1004G>T)` that kept its wrapper.
//!
//! One 9-base transcript, `ATGCGCTAA`, on `NC_000001.11:1000-1008`, so `c.5`
//! is `n.5` is `g.1004`.

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

/// Bare `NM_` input, the direct `c.` path: `c.` is the input itself, so it
/// always kept the wrapper; `n.` is reframed from it and used to drop it.
#[test]
fn a_predicted_coding_input_stays_predicted_on_the_noncoding_axis() {
    assert_eq!(
        axes("NM_TX1.1:c.(5G>T)"),
        (some("NM_TX1.1:c.(5G>T)"), some("NM_TX1.1:n.(5G>T)"), None),
    );
}

/// Parent-framed input, the genome pivot: both transcript axes used to be
/// rebuilt certain from the pivot while `g.` kept the wrapper.
#[test]
fn a_predicted_parent_framed_input_stays_predicted_on_every_axis() {
    assert_eq!(
        axes("NC_000001.11(NM_TX1.1):c.(5del)"),
        (
            some("NC_000001.11(NM_TX1.1):c.(5del)"),
            some("NC_000001.11(NM_TX1.1):n.(5del)"),
            some("NC_000001.11:g.(1004del)"),
        ),
    );
}

/// A predicted `g.` input carries its wrapper onto both transcript axes.
#[test]
fn a_predicted_genomic_input_stays_predicted_on_the_transcript_axes() {
    assert_eq!(
        axes("NC_000001.11:g.(1004G>T)"),
        (
            some("NC_000001.11(NM_TX1.1):c.(5G>T)"),
            some("NC_000001.11(NM_TX1.1):n.(5G>T)"),
            some("NC_000001.11:g.(1004G>T)"),
        ),
    );
}

/// In a cis allele the wrapper belongs to its member alone: the predicted
/// member keeps it on every axis and the observed member does not gain one.
#[test]
fn a_predicted_allele_member_stays_predicted_and_its_sibling_stays_observed() {
    assert_eq!(
        axes("NC_000001.11(NM_TX1.1):c.[(5G>T);7T>C]"),
        (
            some("NC_000001.11(NM_TX1.1):c.[(5G>T);7T>C]"),
            some("NC_000001.11(NM_TX1.1):n.[(5G>T);7T>C]"),
            some("NC_000001.11:g.[(1004G>T);1006T>C]"),
        ),
    );
}

/// The discriminating case: an observed input gains no wrapper anywhere. An
/// over-eager fix (e.g. always wrapping the transcript axes, by analogy with
/// `r.`, which is always predicted) passes every test above and fails this one.
#[test]
fn an_observed_input_gains_no_wrapper() {
    assert_eq!(
        axes("NC_000001.11(NM_TX1.1):c.5G>T"),
        (
            some("NC_000001.11(NM_TX1.1):c.5G>T"),
            some("NC_000001.11(NM_TX1.1):n.5G>T"),
            some("NC_000001.11:g.1004G>T"),
        ),
    );
}
