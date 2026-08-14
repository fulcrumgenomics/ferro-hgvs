//! Sizing the sequence-first DP.
//!
//! `AlignmentDag::build` is `O(n*m)` in time **and space**, but Criterion has
//! no memory axis — this measures wall-clock time only. It shows where that
//! grows too slow to be acceptable, so the block-size cap in the migration is
//! chosen on time evidence; the memory dimension is not measured here (the
//! figures are on #1928: ~18 bytes per grid cell, ~275 MB at the widest
//! accepted block).
//!
//! # Divergence is an AXIS, not a constant (#1928)
//!
//! This file used to fix the alternate at "~5% of positions perturbed" and vary
//! only length. That is one point on a curve whose shape is the whole question,
//! because every cell on a minimal alignment satisfies `|i - j| <= k` for edit
//! distance `k` — so the work that is *reachable* is `Θ(n · k)`, and a benchmark
//! that holds `k / n` fixed makes the cost look inherently quadratic when it is
//! quadratic only in the regime it happens to sample.
//!
//! The three regimes are far apart, and real inputs are at the cheap end:
//!
//! | shape | `k` relative to `n` | where it comes from |
//! |---|---|---|
//! | a real variant | single digits, independent of `n` | a reference and an alternate for one variant differ in a few places |
//! | this file's old fixture | `n / 20` | the 5% perturbation |
//! | uniform-random pair | `~0.6 n` | `dump_normalized_corpus`'s random cores, measured for #107 |
//!
//! So `k` is parameterised directly, and **held constant as `n` grows** in the
//! `by_length_fixed_k` group — that is the shape a real corpus presents, and the
//! one under which a banded implementation should approach linear.
//!
//! # Indels are a second axis
//!
//! The old fixture perturbed by substitution only, so no optimal path ever left
//! the main diagonal and no `OUT_DEL`/`OUT_INS` edge was ever exercised at any
//! length. Substitution-only and indel-bearing divergence are generated
//! separately below: they cost the same to *compute* but produce different DAGs,
//! and only the second has more than one minimal alignment to enumerate.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ferro_hgvs::normalize::seqfirst::align::AlignmentDag;

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

/// Next pseudo-random value, so the perturbation sites are spread rather than
/// periodic. A `step_by` stride puts every difference on a lattice, which is
/// not a property real divergence has.
fn next(state: &mut u64) -> u64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *state >> 33
}

/// An alternate differing from `reference` at exactly `k` positions, by
/// **substitution only** — so the edit distance is `k` and every minimal
/// alignment stays on the main diagonal.
fn diverge_subs(reference: &[u8], k: usize, seed: u64) -> Vec<u8> {
    let mut alt = reference.to_vec();
    if reference.is_empty() || k == 0 {
        return alt;
    }
    let mut state = seed | 1;
    let mut placed = 0;
    // Sample without replacement by rejection; `k` is far below `len` in every
    // configuration here, so the expected number of retries is small.
    let mut seen = vec![false; reference.len()];
    while placed < k {
        let at = (next(&mut state) as usize) % reference.len();
        if seen[at] {
            continue;
        }
        seen[at] = true;
        // Any base other than the one already there, so the difference is real.
        alt[at] = match alt[at] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        placed += 1;
    }
    alt
}

/// An alternate differing from `reference` by `k` **indels** — `k/2` single-base
/// deletions and `k/2` single-base insertions, spread through the block.
///
/// Equal length is deliberate: it keeps the grid square across the two
/// generators so their timings are comparable, and it is the shape the
/// `unchanged-is-read-over-every-minimal-alignment` ruling governs.
fn diverge_indels(reference: &[u8], k: usize, seed: u64) -> Vec<u8> {
    if reference.is_empty() || k < 2 {
        return reference.to_vec();
    }
    let mut state = seed | 1;
    let pairs = k / 2;
    // Delete at `pairs` sites and insert at `pairs` others, walking left to
    // right so the result stays the same length.
    let mut out = Vec::with_capacity(reference.len());
    let stride = reference.len() / (pairs + 1);
    for (i, &b) in reference.iter().enumerate() {
        let slot = i / stride.max(1);
        if slot < pairs && i % stride.max(1) == 0 {
            // Insert a base that differs from the one it precedes, then skip a
            // base — one insertion and one deletion, `stride` apart.
            out.push(b"ACGT"[(next(&mut state) as usize) % 4]);
            continue;
        }
        out.push(b);
    }
    out.truncate(reference.len());
    while out.len() < reference.len() {
        out.push(b'A');
    }
    out
}

/// Cost against **block length, with `k` held constant** — the shape a real
/// corpus presents, and the one under which a banded implementation should
/// approach linear while the current `Θ(n·m)` grid stays quadratic.
fn bench_by_length_fixed_k(c: &mut Criterion) {
    let mut group = c.benchmark_group("seqfirst_align/by_length_fixed_k");
    for &len in &[64usize, 256, 1024, 4096] {
        let reference = bases(len, 1);
        for &k in &[2usize, 8] {
            if k >= len {
                continue;
            }
            let alt = diverge_subs(&reference, k, 7);
            group.throughput(Throughput::Elements((len * len) as u64));
            group.bench_with_input(
                BenchmarkId::from_parameter(format!("len{len}_k{k}")),
                &len,
                |b, _| {
                    b.iter(|| {
                        let dag = AlignmentDag::build(black_box(&reference), black_box(&alt));
                        black_box(dag.edit_distance())
                    })
                },
            );
        }
    }
    group.finish();
}

/// Cost against **divergence, at fixed length** — the axis the old fixture held
/// constant. `k` runs from a real-variant scale up to the uniform-random regime
/// where a band is the whole grid and banding cannot help.
fn bench_by_divergence(c: &mut Criterion) {
    let len = 1024usize;
    let reference = bases(len, 1);
    let mut group = c.benchmark_group("seqfirst_align/by_divergence");
    // 2 and 8 are real-variant scale; 51 is the old fixture's 5%; 614 is the
    // ~0.6n uniform-random regime measured for #107.
    for &k in &[2usize, 8, 32, 51, 128, 512, 614] {
        let alt = diverge_subs(&reference, k, 7);
        group.throughput(Throughput::Elements((len * len) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(format!("k{k}")), &k, |b, _| {
            b.iter(|| {
                let dag = AlignmentDag::build(black_box(&reference), black_box(&alt));
                black_box(dag.edit_distance())
            })
        });
    }
    group.finish();
}

/// The same divergence budget spent on **indels** rather than substitutions.
///
/// Substitution-only divergence keeps every minimal alignment on the diagonal;
/// indels are what make the DAG wide and what `dominators()` exists to
/// enumerate, so the two are not interchangeable as a cost model.
fn bench_indel_vs_substitution(c: &mut Criterion) {
    let len = 1024usize;
    let reference = bases(len, 1);
    let mut group = c.benchmark_group("seqfirst_align/indel_vs_sub");
    for &k in &[8usize, 32, 128] {
        for (kind, alt) in [
            ("sub", diverge_subs(&reference, k, 7)),
            ("indel", diverge_indels(&reference, k, 7)),
        ] {
            group.throughput(Throughput::Elements((len * len) as u64));
            group.bench_with_input(
                BenchmarkId::from_parameter(format!("k{k}_{kind}")),
                &k,
                |b, _| {
                    b.iter(|| {
                        let dag = AlignmentDag::build(black_box(&reference), black_box(&alt));
                        black_box(dag.edit_distance())
                    })
                },
            );
        }
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_by_length_fixed_k,
    bench_by_divergence,
    bench_indel_vs_substitution
);
criterion_main!(benches);
