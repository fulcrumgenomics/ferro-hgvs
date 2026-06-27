//! Issue #851 — NG_/LRG in-place c.→g. projection for compound alleles + UTR
//! coordinates. CI-runnable via MockProvider (no manifest).
//! <https://github.com/fulcrumgenomics/ferro-hgvs/issues/851>

use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::hgvs::variant::HgvsVariant;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::provider::GenomicPlacement;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand as TxStrand, Transcript};
use ferro_hgvs::reference::Strand;
use ferro_hgvs::{parse_hgvs, VariantProjector};

fn parse_acc(s: &str) -> ferro_hgvs::hgvs::variant::Accession {
    match parse_hgvs(&format!("{s}:g.1=")).unwrap() {
        HgvsVariant::Genome(g) => g.accession,
        _ => unreachable!(),
    }
}

/// Plus-strand NM_P.1 on NC_000001.11, 11-base transcript "TTATGCGCTAA",
/// 5'UTR "TT" (cds_start=2, 0-based tx), CDS "ATGCGCTAA" at genome [1000,1011).
/// NG_P.1 placement is identity (parent_start==nc_start) over [900,1100], so the
/// parent coordinate equals the chromosome coordinate; the 5' flank g.999 lives
/// inside the placed span. Genomic sequence: 989 Ns, then "GGGGGGGGGG"
/// (g.990–999), then the transcript bases at g.1000.
fn plus_fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_P.1".to_string(),
        CdotTranscript {
            gene_name: Some("PGENE".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Plus,
            exons: vec![[1000, 1011, 0, 11]],
            cds_start: Some(2),
            cds_end: Some(11),
            gene_id: None,
            protein: Some("NP_P.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_P.1".to_string(),
        Some("PGENE".to_string()),
        TxStrand::Plus,
        "TTATGCGCTAA".to_string(),
        Some(3), // 1-based CDS start (tx index 2)
        Some(11),
        vec![Exon::new(1, 1, 11)],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(1010),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // g.1000 == 0-based index 999. Flank "GGGGGGGGGG" at g.990–999.
    let seq = format!(
        "{}{}{}{}",
        "N".repeat(989),
        "GGGGGGGGGG",
        "TTATGCGCTAA",
        "N".repeat(89)
    );
    provider.add_genomic_sequence("NC_000001.11", seq);
    provider.add_genomic_placement(
        "NG_P.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 900,
            nc_start: 900,
            nc_end: 1100,
            strand: Strand::Plus,
        },
    );
    VariantProjector::new(projector, provider)
}

#[test]
fn plus_strand_cis_allele_routes_to_genomic() {
    let vp = plus_fixture();
    // c.1 = CDS base 1 ('A' at g.1002); c.4 = CDS base 4 ('C' at g.1005). Plus
    // strand: no reverse-complement. Identity placement: parent == chromosome.
    let v = parse_hgvs("NG_P.1(NM_P.1):c.[1A>T;4C>A]").expect("parse");
    let g = vp.project_to_genomic(&v).expect("should project");
    assert_eq!(g.to_string(), "NG_P.1:g.[1002A>T;1005C>A]");
}

#[test]
fn equals_member_allele_projects_to_divergent_form() {
    // Pins the behavior rows 4/5 are annotated for: ferro retains the explicit
    // `=` member and brackets; it does NOT drop the `=` or de-bracket.
    let vp = plus_fixture();
    let v = parse_hgvs("NG_P.1(NM_P.1):c.[1=;4del]").expect("parse");
    let g = vp.project_to_genomic(&v).expect("should project");
    assert_eq!(g.to_string(), "NG_P.1:g.[1002=;1005del]");
}

#[test]
fn genome_allele_projects_idempotently() {
    // A genome allele's members are already Genome → each passes through
    // unchanged; the allele reassembles as itself (routing must not error).
    let vp = plus_fixture();
    let v = parse_hgvs("NC_000001.11:g.[1002A>T;1005C>A]").expect("parse");
    let g = vp.project_to_genomic(&v).expect("idempotent");
    assert_eq!(g.to_string(), "NC_000001.11:g.[1002A>T;1005C>A]");
}

#[test]
fn allele_with_unparented_member_declines_whole() {
    // All-or-nothing: a bare-NM_ member (no genomic parent) makes the whole
    // allele decline rather than silently dropping the member. Assert the
    // concrete missing-parent decline (not bare `is_err()`, which would also
    // pass for an unrelated parse/provider failure and not prove the contract).
    let vp = plus_fixture();
    let v = parse_hgvs("NM_P.1:c.[1A>T;4C>A]").expect("parse");
    let err = vp
        .project_to_genomic(&v)
        .expect_err("an unparented member must decline the whole allele");
    let msg = err.to_string();
    assert!(
        msg.contains("no parent reference"),
        "expected an all-or-nothing decline citing the missing genomic parent, got: {msg}"
    );
}

/// Minus-strand NM_M.1 on NC_000001.11, 20-base transcript, cds_start=2 (0-based
/// tx). tx_to_genome(tx) = 2019 - tx (exon genome [2000,2020)). NG_M.1 identity
/// placement over [1900,2100]. c.1→tx2→g.2017; c.5→tx6→g.2013; c.9→tx10→g.2009
/// (descending in input order → C2 must re-sort ascending).
fn minus_fixture() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_M.1".to_string(),
        CdotTranscript {
            gene_name: Some("MGENE".to_string()),
            contig: "NC_000001.11".to_string(),
            strand: Strand::Minus,
            exons: vec![[2000, 2020, 0, 20]],
            cds_start: Some(2),
            cds_end: Some(20),
            gene_id: None,
            protein: Some("NP_M.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let projector = Projector::new(cdot);

    let mut provider = MockProvider::new();
    // Transcript sequence is the reverse-complement of the genome (minus strand);
    // the sort/flank tests assert positions, not bases.
    provider.add_transcript(Transcript::new(
        "NM_M.1".to_string(),
        Some("MGENE".to_string()),
        TxStrand::Minus,
        "ACGTACGTACGTACGTACGT".to_string(),
        Some(3),
        Some(20),
        vec![Exon::new(1, 1, 20)],
        Some("NC_000001.11".to_string()),
        Some(2000),
        Some(2019),
        Default::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    let seq = format!(
        "{}{}{}",
        "N".repeat(1999),
        "ACGTACGTACGTACGTACGT",
        "N".repeat(80)
    );
    provider.add_genomic_sequence("NC_000001.11", seq);
    provider.add_genomic_placement(
        "NG_M.1",
        GenomicPlacement {
            nc: parse_acc("NC_000001.11"),
            parent_start: 1900,
            nc_start: 1900,
            nc_end: 2100,
            strand: Strand::Plus,
        },
    );
    VariantProjector::new(projector, provider)
}

/// Extract each Genome allele member's start position, in order.
fn member_starts(v: &HgvsVariant) -> Vec<u64> {
    match v {
        HgvsVariant::Allele(a) => a
            .variants
            .iter()
            .map(|m| match m {
                HgvsVariant::Genome(g) => {
                    g.loc_edit.location.start.inner().map(|p| p.base).unwrap()
                }
                other => panic!("expected Genome member, got {other:?}"),
            })
            .collect(),
        other => panic!("expected Allele, got {other:?}"),
    }
}

#[test]
fn minus_strand_allele_members_sorted_genomic_order() {
    let vp = minus_fixture();
    // Input c. order is 1,5,9 → genomic 2017,2013,2009 (descending). After C2 the
    // projected members must be ascending: 2009, 2013, 2017.
    let v = parse_hgvs("NG_M.1(NM_M.1):c.[1A>C;5A>C;9A>C]").expect("parse");
    let g = vp.project_to_genomic(&v).expect("should project");
    let starts = member_starts(&g);
    assert!(
        starts.windows(2).all(|w| w[0] < w[1]),
        "members must be ascending by genomic start, got {starts:?}"
    );
    assert_eq!(starts, vec![2009, 2013, 2017]);
}

#[test]
fn trans_allele_members_not_reordered() {
    // A trans allele `[m1];[m2]` encodes haplotype assignment by member order;
    // the genomic-order sort must NOT touch it (only cis is sorted). c.1→g.2017
    // (haplotype A), c.9→g.2009 (haplotype B) must stay in input order.
    let vp = minus_fixture();
    let v = parse_hgvs("NG_M.1(NM_M.1):c.[1A>C];[9A>C]").expect("parse");
    let g = vp.project_to_genomic(&v).expect("should project");
    assert_eq!(
        member_starts(&g),
        vec![2017, 2009],
        "trans order must be preserved"
    );
}
