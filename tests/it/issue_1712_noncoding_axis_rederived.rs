//! #1712: a projected `n.` axis is normalized on its own non-coding reference
//! (`noncoding-axis-is-re-derived-on-its-own-reference`), not inherited from the
//! form the `c.` axis chose under its reading frame.
//!
//! The `c.` frame merges two changes inside one codon into a `delins`
//! (`DNA/delins.md:18`). The `n.` reference declares no reading frame, so on that
//! axis the same bases re-derive to the two substitutions they are
//! (`general.md:22`). Before this change the `n.` axis carried the `c.` frame's
//! `delins`, which is not a fixed point of normalizing that `n.` description.
//!
//! One 9-base transcript, `ATGCGCTAA`, on `NC_000001.11:1000-1008`.
//! `c.4_6delinsAGA` turns `CGC` into `AGA`: on the `n.` axis that is `4C>A` and
//! `6C>A` with the unchanged `G` between them.

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

fn project(input: &str) -> ferro_hgvs::VariantProjection {
    let v = ferro_hgvs::parse_hgvs(input).unwrap_or_else(|e| panic!("{input} parses: {e}"));
    fixture()
        .project_variant(&v, "NM_TX1.1")
        .unwrap_or_else(|e| panic!("{input} projects: {e}"))
}

fn noncoding(input: &str) -> String {
    project(input)
        .noncoding
        .unwrap_or_else(|| panic!("{input} has an n. axis"))
        .to_string()
}

/// A lone codon `delins` stays a `delins` on `c.` but re-derives to its two
/// substitutions on `n.`, which declares no frame to merge them in.
#[test]
fn a_codon_delins_re_derives_to_its_substitutions_on_the_noncoding_axis() {
    let input = "NC_000001.11(NM_TX1.1):c.4_6delinsAGA";
    assert_eq!(
        project(input).coding.map(|c| c.to_string()).as_deref(),
        Some(input),
        "the c. axis keeps the codon delins its reading frame merges",
    );
    assert_eq!(
        noncoding("NC_000001.11(NM_TX1.1):c.4_6delinsAGA"),
        "NC_000001.11(NM_TX1.1):n.[4C>A;6C>A]",
    );
}

/// In a cis allele the re-derived member folds into the allele as one flat,
/// ordered member list rather than a bracket inside a bracket.
#[test]
fn a_re_derived_member_folds_into_its_cis_allele_on_the_noncoding_axis() {
    assert_eq!(
        noncoding("NC_000001.11(NM_TX1.1):c.[4_6delinsAGA;8A>G]"),
        "NC_000001.11(NM_TX1.1):n.[4C>A;6C>A;8A>G]",
    );
}
