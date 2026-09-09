//! Ruled-partitioner hot-path benchmark (HEAD baseline; NOT swept).
//!
//! VERSION FLOOR: post-flip HEAD only. This exercises the shipped-default `Ruled`
//! partitioner on the indel-rich shapes where the flip's measured overhead lives
//! (the +14% partitioning-heavy regime), so it is the regression vehicle for the
//! performance work in
//! `reports/2026-09-09-ruled-partitioner-performance-opportunities.md`.
//!
//! Unlike `benchmarks.rs::bench_allele_scaling` (which uses `MockProvider::new()`
//! with no sequences and therefore does not drive the sequence-first ruled
//! partitioner), every case here normalizes a `g.` variant against a real contig,
//! so it reaches `collapse_overlapping_cis_edits` -> the ruled consult ->
//! `canonicalize_from_sequence_with_rule` (seed DAG, the fixpoint engine, the
//! shipping chain, and the coalesce passes). Self-contained: a `MockProvider` with
//! one genomic contig, no external prepared reference.
//!
//! Workflow:
//!   cargo bench --bench ruled_partitioner -- --save-baseline pre-opt
//!   # ...make a change...
//!   cargo bench --bench ruled_partitioner -- --baseline pre-opt

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ferro_hgvs::{parse_hgvs, MockProvider, Normalizer};

/// A 500 bp period-4 `ACGT` contig. The period-4 tandem structure means a dup of
/// a 4-mer extends a reference tandem (so `peel_tandem_dup_beside_change` and
/// `dup_extends_reference_tandem` fire), while any range is invertible, so `RunInv`
/// is reachable. Base at 1-based position `p` is `ACGT[(p - 1) % 4]`.
const CONTIG_LEN: usize = 500;
fn contig() -> String {
    const UNIT: &[u8; 4] = b"ACGT";
    (0..CONTIG_LEN).map(|i| UNIT[i % 4] as char).collect()
}

/// Base (as a char) at 1-based genomic position `p` on the period-4 contig.
fn base_at(p: usize) -> char {
    ['A', 'C', 'G', 'T'][(p - 1) % 4]
}

/// A substitution string `<pos><ref>><alt>` whose reference base matches the contig
/// (so it is not treated as a mismatch/correction), altered to the next base.
fn sub(pos: usize) -> String {
    let refb = base_at(pos);
    let alt = match refb {
        'A' => 'C',
        'C' => 'G',
        'G' => 'T',
        _ => 'A',
    };
    format!("{pos}{refb}>{alt}")
}

fn provider() -> MockProvider {
    let mut p = MockProvider::new();
    p.add_genomic_sequence("NC_000001.11", contig());
    p
}

/// A cis allele of `n` net-deletion delins members spaced a few bases apart — the
/// multi-member collapse-consult + payload-coincidence coalesce path, at scale.
fn cis_delins_allele(n: usize) -> String {
    let mut s = String::from("NC_000001.11:g.[");
    for i in 0..n {
        if i > 0 {
            s.push(';');
        }
        let start = 40 + i * 8;
        // 5-base span -> 3-base payload: net deletion, gap-bearing.
        s.push_str(&format!("{}_{}delinsACG", start, start + 4));
    }
    s.push(']');
    s
}

fn bench_single_shapes(c: &mut Criterion) {
    let normalizer = Normalizer::new(provider());

    // Each drives a different arm of the ruled partitioner.
    let cases: Vec<(&str, String)> = vec![
        // Long net-deletion delins: DAG seed + payload-coincidence coalesce.
        (
            "delins_long",
            "NC_000001.11:g.100_140delinsACGTACGTACGTAC".to_string(),
        ),
        // Tandem dup abutting a change (the peel + solid-run collapse + consult).
        (
            "dup_beside_change",
            format!("NC_000001.11:g.[200_203dup;{}]", sub(204)),
        ),
        // Whole-range inversion (RunInv).
        ("inv_range", "NC_000001.11:g.300_320inv".to_string()),
        // A tight two-member cis pair one base apart (SepZero / codon-frame shape).
        (
            "cis_pair_close",
            format!("NC_000001.11:g.[{};{}]", sub(250), sub(252)),
        ),
        // Dup range that extends the reference tandem (dup_extends_reference_tandem).
        ("dup_tandem", "NC_000001.11:g.400_403dup".to_string()),
    ];

    // Guard: every case must parse AND normalize on the real path, or it would be
    // timed on the error path.
    for (name, s) in &cases {
        let v = parse_hgvs(s).unwrap_or_else(|e| panic!("ruled bench {name:?} parse: {e}"));
        normalizer
            .normalize(&v)
            .unwrap_or_else(|e| panic!("ruled bench {name:?} ({s:?}) errored: {e}"));
    }

    let mut group = c.benchmark_group("ruled_single");
    for (name, s) in &cases {
        let v = parse_hgvs(s).unwrap();
        group.bench_function(*name, |b| b.iter(|| normalizer.normalize(black_box(&v))));
    }
    group.finish();
}

fn bench_cis_scaling(c: &mut Criterion) {
    let normalizer = Normalizer::new(provider());
    let mut group = c.benchmark_group("ruled_cis_scaling");
    // Modest n: the collapse loop + consult is superlinear, so large n dominates
    // the whole bench. These are enough to expose the loop-multiplier work.
    for n in [4usize, 12, 32] {
        let input = cis_delins_allele(n);
        let v = parse_hgvs(&input).unwrap_or_else(|e| panic!("ruled cis bench n={n} parse: {e}"));
        normalizer
            .normalize(&v)
            .unwrap_or_else(|e| panic!("ruled cis bench n={n} normalize: {e}"));
        group.bench_with_input(BenchmarkId::new("cis_delins", n), &v, |b, v| {
            b.iter(|| normalizer.normalize(black_box(v)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_single_shapes, bench_cis_scaling);
criterion_main!(benches);
