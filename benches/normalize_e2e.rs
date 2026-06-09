//! End-to-end normalize benchmark over a real `MultiFastaProvider`.
//!
//! This is the surface that actually shows the A2 win (transcript resolved
//! once and shared as an `Arc` instead of rebuilt on every `get_transcript`).
//! The `MockProvider` microbench in `benchmarks.rs` cannot show it — its
//! `get_transcript` is trivially cheap and never rebuilds.
//!
//! Guarded: it skips unless BOTH env vars are set, so it is inert in normal
//! `cargo bench` runs and only measures where real reference data exists
//! (e.g. the nightly, which builds the manifest via `ferro prepare`):
//!   FERRO_BENCH_MANIFEST  = path to a `ferro prepare` manifest.json
//!   FERRO_NORM_CORPUS_TXT = newline-delimited HGVS inputs (one per line)
//!
//! Run:
//!   FERRO_BENCH_MANIFEST=/path/manifest.json \
//!   FERRO_NORM_CORPUS_TXT=/path/inputs.txt \
//!     cargo bench --bench normalize_e2e --features dev

use criterion::{criterion_group, criterion_main, Criterion, Throughput};
use ferro_hgvs::{parse_hgvs, HgvsVariant, MultiFastaProvider, Normalizer};
use std::hint::black_box;

/// Parse up to `limit` HGVS inputs from `FERRO_NORM_CORPUS_TXT`. Unparseable
/// lines are skipped — the benchmark measures `normalize`, not parsing.
fn load_variants(limit: usize) -> Vec<HgvsVariant> {
    let Ok(path) = std::env::var("FERRO_NORM_CORPUS_TXT") else {
        return Vec::new();
    };
    let Ok(text) = std::fs::read_to_string(&path) else {
        return Vec::new();
    };
    text.lines()
        .filter_map(|line| parse_hgvs(line.trim()).ok())
        .take(limit)
        .collect()
}

fn bench_normalize_e2e(c: &mut Criterion) {
    let Ok(manifest) = std::env::var("FERRO_BENCH_MANIFEST") else {
        eprintln!(
            "SKIP normalize_e2e: set FERRO_BENCH_MANIFEST (manifest.json) and \
             FERRO_NORM_CORPUS_TXT (newline-delimited HGVS) to measure"
        );
        return;
    };
    let provider = match MultiFastaProvider::from_manifest(&manifest) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("SKIP normalize_e2e: could not load manifest {manifest}: {e}");
            return;
        }
    };
    let normalizer = Normalizer::new(provider);

    let variants = load_variants(50_000);
    if variants.is_empty() {
        eprintln!("SKIP normalize_e2e: no parseable inputs (set FERRO_NORM_CORPUS_TXT)");
        return;
    }

    let mut group = c.benchmark_group("normalize_e2e");
    group.throughput(Throughput::Elements(variants.len() as u64));
    group.bench_function("normalize_corpus", |b| {
        b.iter(|| {
            let mut ok = 0usize;
            for variant in &variants {
                ok += black_box(normalizer.normalize(black_box(variant))).is_ok() as usize;
            }
            black_box(ok);
        });
    });
    group.finish();
}

criterion_group!(benches, bench_normalize_e2e);
criterion_main!(benches);
