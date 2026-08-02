//! Sizing the sequence-first DP.
//!
//! `AlignmentDag::build` is `O(n*m)` in time **and space**, but Criterion has
//! no memory axis — this measures wall-clock time only. It shows where that
//! grows too slow to be acceptable, so the block-size cap in the migration is
//! chosen on time evidence; the memory dimension is not measured here.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

/// Deterministic pseudo-random bases; no RNG dependency, and reproducible
/// across runs so successive measurements are comparable.
fn bases(len: usize, seed: u64) -> Vec<u8> {
    let mut state = seed | 1;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("seqfirst_align_build");
    for &len in &[16usize, 64, 256, 1024, 4096] {
        let reference = bases(len, 1);
        let mut alt = reference.clone();
        // Perturb ~5% of positions so the DAG is non-trivial.
        for k in (0..len).step_by(20) {
            alt[k] = if alt[k] == b'A' { b'C' } else { b'A' };
        }
        group.throughput(Throughput::Elements((len * len) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(len), &len, |b, _| {
            b.iter(|| {
                let dag = ferro_hgvs::normalize::seqfirst::align::AlignmentDag::build(
                    black_box(&reference),
                    black_box(&alt),
                );
                black_box(dag.edit_distance())
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_build);
criterion_main!(benches);
