//! Performance benchmarks for ferro-hgvs
//!
//! Run with: cargo bench
//! Run specific benchmark: cargo bench -- parsing

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ferro_hgvs::hgvs::parser::parse_hgvs_fast;
use ferro_hgvs::reference::{Exon, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

// =============================================================================
// Parsing benchmarks
// =============================================================================

/// Benchmark HGVS parsing for different variant types
fn bench_parsing(c: &mut Criterion) {
    let variants = vec![
        // Genomic variants
        ("g.sub", "NC_000001.11:g.12345A>G"),
        ("g.del", "NC_000001.11:g.100del"),
        ("g.del_range", "NC_000001.11:g.100_200del"),
        ("g.dup", "NC_000001.11:g.100dup"),
        ("g.dup_range", "NC_000001.11:g.100_200dup"),
        ("g.ins", "NC_000001.11:g.100_101insATG"),
        ("g.delins", "NC_000001.11:g.100_200delinsATG"),
        ("g.inv", "NC_000001.11:g.100_200inv"),
        ("g.repeat", "NC_000001.11:g.100[12]"),
        // CDS variants
        ("c.sub", "NM_000088.3:c.459A>G"),
        ("c.del", "NM_000088.3:c.459del"),
        ("c.ins", "NM_000088.3:c.459_460insATG"),
        ("c.intronic+", "NM_000088.3:c.100+5G>A"),
        ("c.intronic-", "NM_000088.3:c.100-10A>G"),
        ("c.dup", "NM_000088.3:c.100_102dup"),
        // Protein variants
        ("p.sub", "NP_000079.2:p.Val600Glu"),
        ("p.del", "NP_000079.2:p.Val600del"),
        ("p.fs", "NP_000079.2:p.Val600fs"),
        ("p.fsTer", "NP_000079.2:p.Val600fsTer15"),
        ("p.ext", "NP_000079.2:p.Met1ext-5"),
        ("p.identity", "NP_000079.2:p.="),
        ("p.no_protein", "NP_000079.2:p.0"),
        // Non-coding variants
        ("n.sub", "NR_000001.1:n.100A>G"),
        ("n.del", "NR_000001.1:n.100del"),
        // RNA variants
        ("r.sub", "NM_000088.3:r.100a>g"),
        // Ensembl
        ("enst.sub", "ENST00000357033:c.100A>G"),
        // --- coverage gap-fills (churned grammar) ---
        ("m.sub", "NC_012920.1:m.8993T>G"),
        ("allele_cis", "NC_000001.11:g.[100A>G;200del]"),
        ("allele_trans", "NC_000001.11:g.[100A>G];[200del]"),
        ("c.compound_repeat", "NM_000088.3:c.5GA[10]"),
        ("c.uncertain", "NM_000088.3:c.(100_200)del"),
        ("p.predicted", "NP_000079.2:p.(Val600Glu)"),
        ("p.three_letter", "NP_000079.2:p.Val600Glufs*15"),
        ("c.utr3", "NM_000088.3:c.*5A>G"),
        ("c.utr5", "NM_000088.3:c.-5A>G"),
        ("c.deep_intronic", "NM_000088.3:c.100+2000A>G"),
        (
            "g.long_ins",
            "NC_000001.11:g.100_101insATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT",
        ),
        (
            "g.delins_long",
            "NC_000001.11:g.100_120delinsACGTACGTACGTACGTACGT",
        ),
    ];

    // Guard: every benched string must actually parse. A typo or grammar
    // drift that turns a case into a parse error would otherwise be timed
    // on the error path and silently misreport. Fail loud instead.
    for (name, v) in &variants {
        assert!(
            parse_hgvs(v).is_ok(),
            "parsing bench case {name:?} ({v:?}) did not parse"
        );
    }

    let mut group = c.benchmark_group("parsing");

    for (name, variant) in &variants {
        group.bench_with_input(BenchmarkId::new("type", name), variant, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
    }

    group.finish();
}

/// Benchmark parsing by variant string length
fn bench_parsing_by_length(c: &mut Criterion) {
    let variants = vec![
        ("short", "NC_000001.11:g.1A>G"),
        ("medium", "NC_000001.11:g.12345678A>G"),
        ("long_pos", "NC_000001.11:g.123456789_123456799del"),
        ("long_seq", "NC_000001.11:g.100_101insATGCATGCATGC"),
    ];

    let mut group = c.benchmark_group("parsing_length");

    for (name, variant) in &variants {
        let len = variant.len();
        group.throughput(Throughput::Bytes(len as u64));
        group.bench_with_input(BenchmarkId::new("len", name), variant, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
    }

    group.finish();
}

// =============================================================================
// Throughput benchmarks
// =============================================================================

/// Benchmark parsing throughput (variants per second)
fn bench_parsing_throughput(c: &mut Criterion) {
    let variants: Vec<&str> = vec![
        "NC_000001.11:g.12345A>G",
        "NM_000088.3:c.459A>G",
        "NM_000088.3:c.100+5G>A",
        "NP_000079.2:p.Val600Glu",
        "NC_000001.11:g.100del",
        "NM_000088.3:c.100_102dup",
        "NC_000001.11:g.100_101insATG",
        "NP_000079.2:p.Val600fs",
    ];

    let mut group = c.benchmark_group("throughput");

    // Batch of 100 variants
    group.throughput(Throughput::Elements(100));
    group.bench_function("parse_100", |b| {
        b.iter(|| {
            for _ in 0..100 / variants.len() + 1 {
                for variant in &variants {
                    let _ = parse_hgvs(black_box(variant));
                }
            }
        })
    });

    // Batch of 1000 variants
    group.throughput(Throughput::Elements(1000));
    group.bench_function("parse_1000", |b| {
        b.iter(|| {
            for _ in 0..1000 / variants.len() + 1 {
                for variant in &variants {
                    let _ = parse_hgvs(black_box(variant));
                }
            }
        })
    });

    group.finish();
}

// =============================================================================
// Normalization benchmarks
// =============================================================================

/// `with_test_data()` plus a genome-mapped contig and two intron-bearing
/// transcripts (one per strand), so `g.` and intronic normalization run on
/// the real path. Uses only reference APIs present since v0.4.0, so this
/// file overlays onto every release tag during the sweep.
fn enriched_provider() -> MockProvider {
    use ferro_hgvs::reference::transcript::{GenomeBuild, ManeStatus};
    let mut p = MockProvider::with_test_data();

    // Plus-strand, two-exon coding transcript with a real intron, mapped to
    // NC_000001.11. Exon1 genome[1000,1009] tx[1,10]; intron genome[1010,1999];
    // exon2 genome[2000,2009] tx[11,20]. CDS tx[1,18].
    p.add_transcript(Transcript::new(
        "NM_INTR.1".to_string(),
        Some("INTRGENE".to_string()),
        Strand::Plus,
        Some("ATGCGCAAAGGGTAACCC".to_string()), // 18 bp
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
        Some("ATGCGCAAAGGGTAACCC".to_string()),
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

/// Benchmark normalization performance for different variant types
fn bench_normalization(c: &mut Criterion) {
    let provider = enriched_provider();
    let normalizer = Normalizer::new(provider);

    let test_cases = vec![
        // c. plus-strand (NM_000088.3, NM_888888.1 for 3'-shift dup)
        ("c.sub", "NM_000088.3:c.10A>G"),
        ("c.del", "NM_000088.3:c.10del"),
        ("c.del_range", "NM_000088.3:c.10_15del"),
        ("c.ins", "NM_000088.3:c.10_11insATG"),
        ("c.delins", "NM_000088.3:c.10_12delinsAT"),
        ("c.inv", "NM_000088.3:c.10_15inv"),
        ("c.dup", "NM_000088.3:c.10dup"),
        ("c.dup_3shift", "NM_888888.1:c.8dup"),
        // c. minus-strand (NM_999999.1)
        ("c.sub_minus", "NM_999999.1:c.5A>G"),
        ("c.del_minus", "NM_999999.1:c.5del"),
        ("c.ins_minus", "NM_999999.1:c.5_6insAT"),
        // c. intronic (enriched fixture; #497 shuffle path)
        ("c.intronic_plus_del", "NM_INTR.1:c.10+2del"),
        ("c.intronic_minus_del", "NM_INTRM.1:c.10+2del"),
        // protein (NP_000079.2: M1, R97, V600)
        ("p.sub", "NP_000079.2:p.Val600Glu"),
        ("p.del", "NP_000079.2:p.Val600del"),
        ("p.fs", "NP_000079.2:p.Arg97fs"),
        ("p.fsTer", "NP_000079.2:p.Arg97fsTer15"),
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

    let mut group = c.benchmark_group("normalization");
    for (name, variant_str) in &test_cases {
        let variant = parse_hgvs(variant_str).unwrap();
        group.bench_function(*name, |b| {
            b.iter(|| normalizer.normalize(black_box(&variant)))
        });
    }
    group.finish();
}

// =============================================================================
// Display/formatting benchmarks
// =============================================================================

/// Benchmark display/formatting performance
fn bench_display(c: &mut Criterion) {
    let variants: Vec<_> = vec![
        "NC_000001.11:g.12345A>G",
        "NM_000088.3:c.459A>G",
        "NM_000088.3:c.100_102dup",
        "NP_000079.2:p.Val600Glu",
        "NM_000088.3:c.100_101insATGCATGC",
        "NP_000079.2:p.Val600fsTer15",
    ]
    .into_iter()
    .map(|s| parse_hgvs(s).unwrap())
    .collect();

    let mut group = c.benchmark_group("display");

    group.bench_function("single", |b| {
        let variant = &variants[0];
        b.iter(|| format!("{}", black_box(variant)))
    });

    group.bench_function("batch", |b| {
        b.iter(|| {
            for variant in &variants {
                let _ = format!("{}", black_box(variant));
            }
        })
    });

    group.finish();
}

// =============================================================================
// Full pipeline benchmarks
// =============================================================================

/// Benchmark parse + normalize full pipeline
fn bench_full_pipeline(c: &mut Criterion) {
    let provider = MockProvider::with_test_data();
    let normalizer = Normalizer::new(provider);

    let mut group = c.benchmark_group("pipeline");

    // Single variant pipeline
    let variant_str = "NM_000088.3:c.10A>G";
    group.bench_function("single", |b| {
        b.iter(|| {
            let variant = parse_hgvs(black_box(variant_str)).unwrap();
            let _ = normalizer.normalize(&variant);
        })
    });

    // Full roundtrip: parse -> normalize -> display -> parse
    group.bench_function("roundtrip", |b| {
        b.iter(|| {
            let variant = parse_hgvs(black_box(variant_str)).unwrap();
            let normalized = normalizer.normalize(&variant).unwrap();
            let displayed = format!("{}", normalized);
            let _ = parse_hgvs(&displayed);
        })
    });

    group.finish();
}

// =============================================================================
// Comparison benchmarks (for comparing with other libraries)
// =============================================================================

/// Benchmark parsing of clinically relevant variants
fn bench_clinical_variants(c: &mut Criterion) {
    // These are common clinical variant patterns
    let clinical_variants = vec![
        // BRAF V600E (melanoma)
        ("BRAF_V600E", "NP_004324.2:p.Val600Glu"),
        // EGFR exon 19 deletion (lung cancer)
        ("EGFR_del", "NM_005228.5:c.2235_2249del"),
        // BRCA1 frameshift
        ("BRCA1_fs", "NM_007294.4:c.68_69del"),
        // KRAS G12D
        ("KRAS_G12D", "NP_004976.2:p.Gly12Asp"),
        // TP53 splice variant
        ("TP53_splice", "NM_000546.6:c.375+1G>A"),
        // Complex insertion
        ("complex_ins", "NM_000088.3:c.100_101insATGCATGCATGC"),
    ];

    let mut group = c.benchmark_group("clinical");

    for (name, variant) in &clinical_variants {
        group.bench_with_input(BenchmarkId::new("parse", name), variant, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
    }

    group.finish();
}

/// Benchmark fast-path parser vs standard parser
fn bench_fast_path(c: &mut Criterion) {
    // Variants that should use fast path
    let fast_path_variants = vec![
        // RefSeq substitutions (fast path)
        ("refseq_g.sub", "NC_000001.11:g.12345A>G"),
        ("refseq_c.sub", "NM_000088.3:c.459A>G"),
        // Ensembl substitutions (fast path)
        ("ensembl_g.sub", "ENSG00000141510.5:g.12345A>G"),
        ("ensembl_c.sub", "ENST00000357033.1:c.100A>G"),
        // LRG (fast path)
        ("lrg_g.sub", "LRG_1:g.12345A>G"),
        // Assembly (fast path)
        ("grch38_g.sub", "GRCh38(chr1):g.12345A>G"),
    ];

    // Variants that fall back to standard parser
    let fallback_variants = vec![
        ("c.intronic", "NM_000088.3:c.100+5G>A"),
        ("c.utr", "NM_000088.3:c.*100A>G"),
        ("g.del", "NC_000001.11:g.100del"),
        ("g.ins", "NC_000001.11:g.100_101insATG"),
        ("p.sub", "NP_000079.2:p.Val600Glu"),
    ];

    let mut group = c.benchmark_group("fast_path");

    // Benchmark fast-path variants with both parsers
    for (name, variant) in &fast_path_variants {
        group.bench_with_input(BenchmarkId::new("standard", name), variant, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
        group.bench_with_input(BenchmarkId::new("fast", name), variant, |b, v| {
            b.iter(|| parse_hgvs_fast(black_box(v)))
        });
    }

    // Benchmark fallback variants (should be same speed)
    for (name, variant) in &fallback_variants {
        group.bench_with_input(BenchmarkId::new("fallback_std", name), variant, |b, v| {
            b.iter(|| parse_hgvs(black_box(v)))
        });
        group.bench_with_input(BenchmarkId::new("fallback_fast", name), variant, |b, v| {
            b.iter(|| parse_hgvs_fast(black_box(v)))
        });
    }

    group.finish();
}

/// Benchmark error handling (invalid variants)
fn bench_error_handling(c: &mut Criterion) {
    let invalid_variants = vec![
        ("empty", ""),
        ("no_accession", "g.12345A>G"),
        ("invalid_prefix", "XX_000001.1:g.12345A>G"),
        ("invalid_pos", "NC_000001.11:g.0A>G"),
        ("invalid_base", "NC_000001.11:g.12345X>G"),
        ("malformed", "this is not HGVS"),
    ];

    let mut group = c.benchmark_group("errors");

    for (name, variant) in &invalid_variants {
        group.bench_with_input(BenchmarkId::new("parse", name), variant, |b, v| {
            b.iter(|| {
                let _ = parse_hgvs(black_box(v));
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_parsing,
    bench_parsing_by_length,
    bench_parsing_throughput,
    bench_normalization,
    bench_display,
    bench_full_pipeline,
    bench_clinical_variants,
    bench_fast_path,
    bench_error_handling,
);

criterion_main!(benches);
