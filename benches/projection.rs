//! Projection benchmarks (v0.6.0+ subsystem). Establishes a HEAD baseline.
//! Run with: cargo bench --bench projection

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ferro_hgvs::data::cdot::{CdotMapper, CdotTranscript};
use ferro_hgvs::data::projection::Projector;
use ferro_hgvs::project::VariantProjector;
use ferro_hgvs::reference::transcript::{GenomeBuild, ManeStatus};
use ferro_hgvs::reference::{Exon, MockProvider, Strand, Transcript};

/// Plus-strand 9bp CDS "ATGCGCTAA" on chr1 [1000,1009), with NP_TEST.1 protein.
///
/// Layout (0-based, half-open):
///   genome chr1[999,1008) = ATGCGCTAA (genome g.1000..g.1008, HGVS 1-based inclusive)
///   single exon: tx[1,9], genome[1000,1008]
///   CDS: tx.cds_start=1 (Transcript, 1-based) / cdot.cds_start=0 (CdotTranscript, 0-based); cds_end=9
fn plus_projector() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_TEST.1".to_string(),
        CdotTranscript {
            gene_name: Some("TESTGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [genome_start(incl), genome_end(excl), tx_start(0-based), tx_end(0-based excl)]
            exons: vec![[1000, 1009, 0, 9]],
            cds_start: Some(0),
            cds_end: Some(9),
            gene_id: None,
            protein: Some("NP_TEST.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(
        Transcript::new(
            "NM_TEST.1".to_string(),
            Some("TESTGENE".to_string()),
            Strand::Plus,
            Some("ATGCGCTAA".to_string()),
            Some(1),
            Some(9),
            // Exon tx[1,9] maps to genome[1000,1008] (1-based inclusive)
            vec![Exon::with_genomic(1, 1, 9, 1000, 1008)],
            Some("chr1".to_string()),
            Some(1000),
            Some(1008),
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        )
        .with_protein_id(Some("NP_TEST.1".to_string())),
    );
    provider.add_protein("NP_TEST.1", "MR*");
    // chr1 sequence: 999 N's + ATGCGCTAA + 100 N's
    // HGVS positions g.1000..g.1008 = ATGCGCTAA
    provider.add_genomic_sequence(
        "chr1",
        format!("{}{}{}", "N".repeat(999), "ATGCGCTAA", "N".repeat(100)),
    );
    VariantProjector::new(Projector::new(cdot), provider)
}

/// Two-exon transcript with a 990bp intron on chr1, for intronic g→c projection.
///
/// Layout (1-based, HGVS inclusive):
///   exon1: tx[1,10], genome g.1000..g.1009
///   intron: genome g.1010..g.1999
///   exon2: tx[11,20], genome g.2000..g.2009
///   CDS: tx.cds_start=1 (Transcript, 1-based) / cdot.cds_start=0 (CdotTranscript, 0-based); cds_end=18 (exon1[1..10] + exon2[11..18])
///
/// Intronic input g.1015A>G maps to c.10+6A>G (6bp into the intron after exon1).
/// Note: the intron is N-filled, so the stated `A` reference base isn't real. The
/// projector does coordinate-only g→c mapping and does not validate the intronic
/// reference allele here, so the `A` is unchecked by design (the input string works as-is).
fn intronic_projector() -> VariantProjector<MockProvider> {
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_INTR.1".to_string(),
        CdotTranscript {
            gene_name: Some("INTRGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [genome_start(incl), genome_end(excl), tx_start(0-based), tx_end(0-based excl)]
            exons: vec![[1000, 1010, 0, 10], [2000, 2010, 10, 20]],
            cds_start: Some(0),
            cds_end: Some(18),
            gene_id: None,
            protein: Some("NP_INTR.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(Transcript::new(
        "NM_INTR.1".to_string(),
        Some("INTRGENE".to_string()),
        Strand::Plus,
        Some("ATGCGCAAAGGGTAACCC".to_string()),
        Some(1),
        Some(18),
        vec![
            Exon::with_genomic(1, 1, 10, 1000, 1009),
            Exon::with_genomic(2, 11, 20, 2000, 2009),
        ],
        Some("chr1".to_string()),
        Some(1000),
        Some(2009),
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    ));
    // chr1 sequence: 999 N's + exon1(10bp) + intron(990bp N's) + exon2(10bp) + 100 N's
    let exon1 = "ATGCGCAAAG"; // genome g.1000..g.1009
    let intron = "N".repeat(990); // genome g.1010..g.1999
    let exon2 = "GGTAACCCNN"; // genome g.2000..g.2009
    provider.add_genomic_sequence(
        "chr1",
        format!(
            "{}{}{}{}{}",
            "N".repeat(999),
            exon1,
            intron,
            exon2,
            "N".repeat(100)
        ),
    );
    VariantProjector::new(Projector::new(cdot), provider)
}

fn bench_projection(c: &mut Criterion) {
    let plus = plus_projector();
    let intr = intronic_projector();

    // Guard: each input must project on the real path.
    // plus:    g.1003C>A is position 4 of ATGCGCTAA (0-indexed 3 = 'C') → c.4C>A
    // intronic: g.1015A>G is 6bp into the intron after exon1 (c.10+6A>G)
    let cases: Vec<(&str, &VariantProjector<MockProvider>, &str, &str)> = vec![
        ("g_to_c_plus", &plus, "chr1:g.1003C>A", "NM_TEST.1"),
        ("g_to_c_intronic", &intr, "chr1:g.1015A>G", "NM_INTR.1"),
    ];
    for (name, vp, input, tx) in &cases {
        vp.project(input, tx)
            .unwrap_or_else(|e| panic!("projection bench {name:?} ({input:?}) errored: {e}"));
    }

    let mut group = c.benchmark_group("projection");
    for (name, vp, input, tx) in &cases {
        group.bench_with_input(BenchmarkId::new("project", name), input, |b, i| {
            b.iter(|| vp.project(black_box(i), black_box(tx)))
        });
    }
    group.finish();
}

criterion_group!(projection, bench_projection);
criterion_main!(projection);
