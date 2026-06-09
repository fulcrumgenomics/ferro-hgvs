//! Performance benchmarks for ferro-hgvs features added after v0.4.0.
//!
//! VERSION FLOOR: post-v0.4.0 — HEAD BASELINE target (NOT swept). These exercise features added after v0.4.0: enriched genomic/intronic normalization (Transcript::new signature drift makes it uncompilable at v0.4.0), W30xx correctors (the swapped-position corrector errors at v0.4.0), and trans-allele parsing (errors at v0.4.0). No v0.4.0 baseline exists, so this establishes a forward floor at HEAD only. Run: cargo bench --bench baseline_head

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ferro_hgvs::reference::{Exon, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

// =============================================================================
// Enriched normalization benchmarks (genomic + intronic)
// =============================================================================

/// `with_test_data()` plus a genome-mapped contig and two intron-bearing
/// transcripts (one per strand), so `g.` and intronic normalization run on
/// the real path. Uses the post-v0.4.0 `Transcript::new` signature
/// (`Option<String>` gene name), so this only compiles at HEAD.
fn enriched_provider() -> MockProvider {
    use ferro_hgvs::reference::transcript::{GenomeBuild, ManeStatus};
    let mut p = MockProvider::with_test_data();

    // Plus-strand, two-exon coding transcript with a real intron, mapped to
    // NC_000001.11. Exon1 genome[1000,1009] tx[1,10]; intron genome[1010,1999];
    // exon2 genome[2000,2009] tx[11,20]. CDS tx[1,18]; tx[19,20] is a 2-nt 3'UTR tail.
    p.add_transcript(Transcript::new(
        "NM_INTR.1".to_string(),
        Some("INTRGENE".to_string()),
        Strand::Plus,
        Some("ATGCGCAAAGGGTAACCCNN".to_string()), // 20 bp (18 bp CDS + 2 nt 3'UTR tail)
        Some(1),
        Some(18),
        vec![
            Exon::with_genomic(1, 1, 10, 1000, 1009),
            Exon::with_genomic(2, 11, 20, 2000, 2009),
        ],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(2009),
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    ));

    // Minus-strand counterpart for the #497 minus-strand intronic shuffle path.
    // On minus strand exons are listed 5'→3' in transcript coordinates, but
    // their genomic coords are reversed: exon1 (tx 1..10) maps to genome 2000..2009,
    // exon2 (tx 11..20) maps to genome 1000..1009.
    p.add_transcript(Transcript::new(
        "NM_INTRM.1".to_string(),
        Some("INTRGENM".to_string()),
        Strand::Minus,
        Some("ATGCGCAAAGGGTAACCCNN".to_string()), // 20 bp (18 bp CDS + 2 nt 3'UTR tail)
        Some(1),
        Some(18),
        vec![
            Exon::with_genomic(1, 1, 10, 2000, 2009),
            Exon::with_genomic(2, 11, 20, 1000, 1009),
        ],
        Some("NC_000001.11".to_string()),
        Some(1000),
        Some(2009),
        GenomeBuild::default(),
        ManeStatus::default(),
        None,
        None,
    ));

    // Genomic sequence for NC_000001.11: 999 N's + exon1 + 990 N intron +
    // exon2 + 100 N's, so positions 1000..2009 carry the modeled bases.
    let prefix = "N".repeat(999);
    let exon1 = "ATGCGCAAAG"; // genome 1000..1009
    let intron = "N".repeat(990); // genome 1010..1999
    let exon2 = "GGTAACCCNN"; // genome 2000..2009
    let suffix = "N".repeat(100);
    p.add_genomic_sequence(
        "NC_000001.11",
        format!("{prefix}{exon1}{intron}{exon2}{suffix}"),
    );
    p
}

/// Benchmark genomic + intronic normalization on the enriched provider.
fn bench_normalization_enriched(c: &mut Criterion) {
    let provider = enriched_provider();
    let normalizer = Normalizer::new(provider);

    let test_cases = vec![
        // c. intronic (enriched fixture; #497 shuffle path)
        ("c.intronic_plus_del", "NM_INTR.1:c.10+2del"),
        ("c.intronic_minus_del", "NM_INTRM.1:c.10+2del"),
        // g. (enriched fixture, NC_000001.11)
        ("g.sub", "NC_000001.11:g.1000A>G"),
        ("g.del", "NC_000001.11:g.1001del"),
        ("g.ins", "NC_000001.11:g.1001_1002insTT"),
        ("g.delins", "NC_000001.11:g.1001_1003delinsTT"),
        ("g.inv", "NC_000001.11:g.1001_1003inv"),
    ];

    // Guard: every case must parse AND normalize on the real path. A case
    // that errors here (bad coordinate, missing accession) would otherwise be
    // timed on the error path. Fail loud.
    for (name, s) in &test_cases {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("normalize bench {name:?} parse: {e}"));
        normalizer
            .normalize(&v)
            .unwrap_or_else(|e| panic!("normalize bench {name:?} ({s:?}) errored: {e}"));
    }

    let mut group = c.benchmark_group("normalization_enriched");
    for (name, variant_str) in &test_cases {
        let variant = parse_hgvs(variant_str).unwrap();
        group.bench_function(*name, |b| {
            b.iter(|| normalizer.normalize(black_box(&variant)))
        });
    }
    group.finish();
}

// =============================================================================
// Corrector benchmarks (W30xx machinery)
// =============================================================================

/// Parse-path timing for inputs that exercise the W30xx correction/diagnostic
/// machinery (swapped positions, dup/del soft-prohibition suffix forms).
fn bench_correctors(c: &mut Criterion) {
    // These inputs exercise the W30xx correction/diagnostic machinery:
    //   swapped_positions (W3026): Ok-with-warning in WarnCorrect mode, Err in Reject mode.
    //   dup_count_suffix / del_seq_suffix / dup_seq_suffix: Ok-with-warning under the lenient default.
    // We time the full parse path regardless of Ok/Err — that is intentional, so there is no is_ok guard.
    let cases = vec![
        ("swapped_positions", "NC_000001.11:g.200_100del"),
        ("dup_count_suffix", "NM_000088.3:c.100_102dup3"),
        ("del_seq_suffix", "NM_000088.3:c.100delA"),
        ("dup_seq_suffix", "NM_000088.3:c.100_102dupATG"),
    ];
    let mut group = c.benchmark_group("correctors");
    for (name, input) in &cases {
        group.bench_with_input(BenchmarkId::new("parse", name), input, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
    }
    group.finish();
}

// =============================================================================
// Parsing benchmarks for post-v0.4.0 grammar
// =============================================================================

/// Parse-path timing for grammar constructs that error at v0.4.0.
fn bench_parsing_baseline(c: &mut Criterion) {
    let variants = vec![("allele_trans", "NC_000001.11:g.[100A>G];[200del]")];

    // Guard: every benched string must actually parse. A typo or grammar
    // drift that turns a case into a parse error would otherwise be timed
    // on the error path and silently misreport. Fail loud instead.
    for (name, v) in &variants {
        assert!(
            parse_hgvs(v).is_ok(),
            "parsing bench case {name:?} ({v:?}) did not parse"
        );
    }

    let mut group = c.benchmark_group("parsing_baseline");
    for (name, variant) in &variants {
        group.bench_with_input(BenchmarkId::new("type", name), variant, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
    }
    group.finish();
}

criterion_group!(
    baseline_head,
    bench_normalization_enriched,
    bench_correctors,
    bench_parsing_baseline,
);

criterion_main!(baseline_head);
