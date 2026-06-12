//! Projection benchmarks (v0.6.0+ subsystem). Establishes a HEAD baseline.
//! Run with: cargo bench --bench projection
//!
//! VERSION FLOOR: v0.6.0 — HEAD BASELINE target (NOT swept). The project module is new in v0.6.0. Run: cargo bench --bench projection

use criterion::{
    black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput,
};
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
    provider.add_transcript(
        Transcript::new(
            "NM_INTR.1".to_string(),
            Some("INTRGENE".to_string()),
            Strand::Plus,
            // 20 bases to match the exon tx-coordinate spans (tx 1..=20) and the
            // genomic reconstruction below (exon1 `ATGCGCAAAG` + exon2
            // `GGTAACCCNN`); CDS still ends at 18, leaving the trailing `NN` as 3'UTR.
            Some("ATGCGCAAAGGGTAACCCNN".to_string()),
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
        )
        .with_protein_id(Some("NP_INTR.1".to_string())),
    );
    // Register the protein the CdotTranscript references (parallels
    // `plus_projector`), so the fixture is self-consistent and a future
    // protein-projection bench over this transcript would resolve.
    // CDS tx 1..=18 = ATG CGC AAA GGG TAA CCC -> M R K G * P (stop at codon 5).
    provider.add_protein("NP_INTR.1", "MRKG*P");
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

/// Bare-`NM_` coding transcript with a real multi-codon CDS **plus a 3'UTR**,
/// for direct `c→p` / `c↔n` projection (no genomic context required).
///
/// Sequence/coords cribbed verbatim from the #587 frameshift-into-3'UTR
/// regression test `frameshift_new_stop_in_3utr_emits_fster_count`
/// (`src/project/protein/indel.rs`), so the codon biology is known-valid
/// rather than hand-derived.
///
/// Transcript "ATGAAGAAGAAGTGACTAACCC" (22bp), cds_start/cds_end frame it as:
///   CDS    "ATGAAGAAGAAGTGA" (tx 1..=15) = Met Lys Lys Lys Ter → protein "MKKK*"
///   3'UTR  "CTAACCC"         (tx 16..=22)
///
/// Codon layout (c. = tx pos here; cds_start at tx pos 1, no 5'UTR):
///   c.1-3   ATG  Met (codon 1)
///   c.4-6   AAG  Lys (codon 2)   ← all three bench inputs hit codon 2
///   c.7-9   AAG  Lys (codon 3)
///   c.10-12 AAG  Lys (codon 4)
///   c.13-15 TGA  Ter (codon 5)
///
/// Bench inputs (all guarded to project to a protein):
///   missense SNV          c.4A>C → codon 2 AAG(Lys)→CAG(Gln) ⇒ p.(Lys2Gln)
///   in-frame 3-base del   c.4_6del → delete codon 2 (Lys)    ⇒ p.(Lys2del)
///   frameshift-into-3'UTR c.4del → new stop in the 3'UTR     ⇒ p.(Lys2ArgfsTer5)
///
/// A genome alignment (chr1) is included too so the same fixture could serve a
/// `c↔g` reuse, but the direct `c→p`/`c↔n` paths exercised here do not need it.
fn coding_projector() -> VariantProjector<MockProvider> {
    const CDS_3UTR: &str = "ATGAAGAAGAAGTGACTAACCC"; // 22bp: CDS(15) + 3'UTR(7)
    let mut cdot = CdotMapper::new();
    cdot.add_transcript(
        "NM_CODE.1".to_string(),
        CdotTranscript {
            gene_name: Some("CODEGENE".to_string()),
            contig: "chr1".to_string(),
            strand: Strand::Plus,
            // [genome_start(incl), genome_end(excl), tx_start(0-based), tx_end(0-based excl)]
            exons: vec![[1000, 1022, 0, 22]],
            cds_start: Some(0), // 0-based: CDS starts at tx pos 0 (no 5'UTR)
            cds_end: Some(15),  // 0-based exclusive: CDS spans tx [0,15)
            gene_id: None,
            protein: Some("NP_CODE.1".to_string()),
            exon_cigars: Vec::new(),
        },
    );
    let mut provider = MockProvider::new();
    provider.add_transcript(
        Transcript::new(
            "NM_CODE.1".to_string(),
            Some("CODEGENE".to_string()),
            Strand::Plus,
            Some(CDS_3UTR.to_string()),
            Some(1),  // 1-based inclusive per Transcript convention
            Some(15), // CDS ends at tx pos 15 (inclusive of the stop codon)
            // Exon tx[1,22] maps to genome[1000,1021] (1-based inclusive)
            vec![Exon::with_genomic(1, 1, 22, 1000, 1021)],
            Some("chr1".to_string()),
            Some(1000),
            Some(1021),
            GenomeBuild::default(),
            ManeStatus::default(),
            None,
            None,
        )
        .with_protein_id(Some("NP_CODE.1".to_string())),
    );
    // CDS "ATGAAGAAGAAGTGA" = Met Lys Lys Lys Ter.
    provider.add_protein("NP_CODE.1", "MKKK*");
    // chr1 sequence: 999 N's + transcript(22bp) + 100 N's (g.1000..g.1021 = transcript).
    provider.add_genomic_sequence(
        "chr1",
        format!("{}{}{}", "N".repeat(999), CDS_3UTR, "N".repeat(100)),
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

/// The three bench inputs on `coding_projector()`, each hitting codon 2.
/// (name, c. input string) — all must project to a protein consequence.
const C_TO_P_INPUTS: [(&str, &str); 3] = [
    ("missense_snv", "NM_CODE.1:c.4A>C"),    // p.(Lys2Gln)
    ("inframe_del", "NM_CODE.1:c.4_6del"),   // p.(Lys2del)
    ("frameshift_3utr", "NM_CODE.1:c.4del"), // p.(Lys2ArgfsTer5)
];

const C_TO_P_TX: &str = "NM_CODE.1";

/// `c→p` direct-path benchmarks: a missense SNV, an in-frame 3-base deletion,
/// and the #587 frameshift-into-3'UTR single-base deletion. Each guarded to
/// actually predict a protein (panic otherwise) before benching.
fn bench_c_to_p(c: &mut Criterion) {
    let vp = coding_projector();

    // Guard: every input must project AND yield a protein consequence.
    for (name, input) in &C_TO_P_INPUTS {
        let proj = vp
            .project(input, C_TO_P_TX)
            .unwrap_or_else(|e| panic!("c_to_p bench {name:?} ({input:?}) errored: {e}"));
        assert!(
            proj.protein.is_some(),
            "c_to_p bench {name:?} ({input:?}) did not predict a protein"
        );
    }

    let mut group = c.benchmark_group("c_to_p");
    for (name, input) in &C_TO_P_INPUTS {
        group.bench_with_input(BenchmarkId::new("c_to_p", name), input, |b, i| {
            b.iter(|| vp.project(black_box(i), black_box(C_TO_P_TX)))
        });
    }
    group.finish();
}

/// `c↔n` benchmark: projecting a bare-`c.` SNV populates the n. (transcript-
/// relative) axis even with no genome alignment. Guarded on `.noncoding`.
fn bench_c_to_n(c: &mut Criterion) {
    let vp = coding_projector();
    let input = "NM_CODE.1:c.4A>C";

    // Guard: the noncoding (n.) axis must be populated.
    let proj = vp
        .project(input, C_TO_P_TX)
        .unwrap_or_else(|e| panic!("c_to_n bench ({input:?}) errored: {e}"));
    assert!(
        proj.noncoding.is_some(),
        "c_to_n bench ({input:?}) did not populate the noncoding (n.) axis"
    );

    let mut group = c.benchmark_group("c_to_n");
    group.bench_with_input(BenchmarkId::new("c_to_n", "snv"), &input, |b, i| {
        b.iter(|| vp.project(black_box(i), black_box(C_TO_P_TX)))
    });
    group.finish();
}

/// Cold-vs-warm ref-AA cache benchmark over the frameshift `c→p` input.
///
/// cold: each iteration builds a FRESH `coding_projector()` (cold ref-AA /
///       `RefProteinBundle` translation), then projects once.
/// warm: build the projector ONCE, prime the cache, then project repeatedly.
fn bench_cache(c: &mut Criterion) {
    let input = "NM_CODE.1:c.4del"; // frameshift-into-3'UTR (p.(Lys2ArgfsTer5))

    // Guard via a one-off warm projector.
    let primed = coding_projector();
    primed
        .project(input, C_TO_P_TX)
        .unwrap_or_else(|e| panic!("cache bench ({input:?}) errored: {e}"))
        .protein
        .expect("cache bench input must predict a protein");

    let mut group = c.benchmark_group("cache");

    // Cold: fixture build is the per-iter setup; only the project() is timed.
    group.bench_function("cold", |b| {
        b.iter_batched(
            coding_projector,
            |vp| vp.project(black_box(input), black_box(C_TO_P_TX)),
            BatchSize::SmallInput,
        )
    });

    // Warm: one projector, cache primed once, then projected in a tight loop.
    let warm = coding_projector();
    let _ = warm.project(input, C_TO_P_TX); // prime the ref-AA cache
    group.bench_function("warm", |b| {
        b.iter(|| warm.project(black_box(input), black_box(C_TO_P_TX)))
    });

    group.finish();
}

/// Batch-throughput benchmark: project 75 distinct `c.` inputs (substitutions,
/// deletions, and duplications, varying the CDS position so each is a real,
/// distinct projection) and report variants/sec via `Throughput::Elements`.
///
/// NOTE: `VariantProjector` exposes no batch/parallel projection API (checked
/// `src/project/projector.rs`: `project` / `project_variant` /
/// `project_variant_all` are all single-variant; the `parallel`/rayon feature
/// is not wired into the projector). Batch projection is therefore SERIAL here;
/// a rayon-backed batch API is a future enhancement.
fn bench_batch(c: &mut Criterion) {
    // 75 DISTINCT inputs over the CDS. The CDS is only 15bp, so to get many
    // distinct projections we sweep several edit kinds across every CDS position
    // (c.1..=15) rather than substitutions alone:
    //   - substitution to each of the 3 non-ref bases  (3 per position)
    //   - single-base deletion  (c.{pos}del)           (1 per position)
    //   - single-base duplication  (c.{pos}dup)        (1 per position)
    // Over the 15 CDS positions that is 15*(3+1+1) = 75 distinct input strings.
    // We restrict to the CDS (positions map directly to c.N); the 3'UTR would
    // need c.*N coordinates. Subs reference the real ref base; del/dup are
    // coordinate-only — so each parses and projects (the guard below enforces it).
    const CDS: &[u8] = b"ATGAAGAAGAAGTGA"; // 15bp CDS (tx/c. 1..=15)
    let bases = ['A', 'C', 'G', 'T'];
    let mut inputs: Vec<String> = Vec::new();
    for (idx, &refb) in CDS.iter().enumerate() {
        let c_pos = idx + 1; // 1-based c. position
        let refc = refb as char;
        for &alt in &bases {
            if alt != refc {
                inputs.push(format!("NM_CODE.1:c.{c_pos}{refc}>{alt}"));
            }
        }
        inputs.push(format!("NM_CODE.1:c.{c_pos}del"));
        inputs.push(format!("NM_CODE.1:c.{c_pos}dup"));
    }

    let vp = coding_projector();

    // Guard: every batch input must project (allele-ref mismatches would error).
    for input in &inputs {
        vp.project(input, C_TO_P_TX)
            .unwrap_or_else(|e| panic!("batch bench ({input:?}) errored: {e}"));
    }

    let mut group = c.benchmark_group("batch");
    group.throughput(Throughput::Elements(inputs.len() as u64));
    group.bench_function("serial", |b| {
        b.iter(|| {
            for input in &inputs {
                let _ = black_box(vp.project(black_box(input), black_box(C_TO_P_TX)));
            }
        })
    });
    group.finish();
}

criterion_group!(
    projection,
    bench_projection,
    bench_c_to_p,
    bench_c_to_n,
    bench_cache,
    bench_batch
);
criterion_main!(projection);
