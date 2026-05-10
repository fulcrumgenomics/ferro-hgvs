//! Paraphase Exhaustive Tests
//!
//! Tests against ALL ClinVar variants for Paraphase-supported genes.
//! This is the exhaustive version with 435K variants from 35 segdup genes.
//!
//! Run with: cargo test --release --test paraphase_exhaustive_tests -- --nocapture

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

// Slim deserialization shape: drop fields the test never reads.
// Mostly a readability win; serde's `IgnoredAny` still walks unread
// JSON tokens, so the perf gain is small.
#[derive(Deserialize)]
struct ParaphaseFixture {
    total_target_genes: usize,
    genes_with_variants: usize,
    gene_descriptions: HashMap<String, String>,
    test_cases: Vec<ParaphaseCase>,
}

#[derive(Deserialize)]
struct ParaphaseCase {
    input: String,
    #[serde(rename = "type")]
    coord_type: String,
    gene: String,
}

fn load_fixture() -> Option<ParaphaseFixture> {
    let path = "tests/fixtures/validation/paraphase_genes_exhaustive.json.gz";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    // See cmrg_exhaustive_tests::load_fixture for why we decompress to
    // a Vec and use `from_slice` rather than streaming via `from_reader`.
    let file = File::open(path).expect("Failed to open paraphase_genes_exhaustive.json.gz");
    let mut buf = Vec::new();
    GzDecoder::new(file)
        .read_to_end(&mut buf)
        .expect("Failed to decompress paraphase_genes_exhaustive.json.gz");
    Some(serde_json::from_slice(&buf).expect("Failed to parse paraphase_genes_exhaustive.json.gz"))
}

#[test]
fn test_paraphase_exhaustive_benchmark() {
    let fixture = match load_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping: paraphase_genes_exhaustive.json not found");
            return;
        }
    };

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("Paraphase Exhaustive Benchmark");
    eprintln!("========================================");
    eprintln!("Total target genes: {}", fixture.total_target_genes);
    eprintln!("Genes with variants: {}", fixture.genes_with_variants);
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    // Parse all cases (parallel with rayon when available). Aggregate the
    // per-coord-type and per-gene counts serially over the resulting
    // Vec<bool>; the aggregation pass is microseconds.
    #[cfg(feature = "parallel")]
    let oks: Vec<bool> = fixture
        .test_cases
        .par_iter()
        .map(|case| parse_hgvs(&case.input).is_ok())
        .collect();
    #[cfg(not(feature = "parallel"))]
    let oks: Vec<bool> = fixture
        .test_cases
        .iter()
        .map(|case| parse_hgvs(&case.input).is_ok())
        .collect();

    let mut passed = 0usize;
    let mut by_type: HashMap<String, (usize, usize)> = HashMap::new();
    let mut by_gene: HashMap<String, (usize, usize)> = HashMap::new();
    for (case, ok) in fixture.test_cases.iter().zip(oks.iter()) {
        let type_entry = by_type.entry(case.coord_type.clone()).or_insert((0, 0));
        let gene_entry = by_gene.entry(case.gene.clone()).or_insert((0, 0));
        if *ok {
            passed += 1;
            type_entry.0 += 1;
            gene_entry.0 += 1;
        } else {
            type_entry.1 += 1;
            gene_entry.1 += 1;
        }
    }

    let elapsed = start.elapsed();
    let rate = total as f64 / elapsed.as_secs_f64();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    eprintln!("\nPerformance:");
    eprintln!("  Time: {:.2}s", elapsed.as_secs_f64());
    eprintln!("  Rate: {:.0} variants/sec", rate);

    eprintln!("\nResults:");
    eprintln!("  Passed: {}/{} ({:.2}%)", passed, total, pass_rate);

    eprintln!("\nBy coordinate type:");
    let mut types: Vec<_> = by_type.iter().collect();
    types.sort_by_key(|b| std::cmp::Reverse(b.1 .0 + b.1 .1));
    for (t, (p, f)) in types.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.2}%)", t, p, tot, r);
    }

    eprintln!("\nTop 10 genes by variant count:");
    let mut genes: Vec<_> = by_gene.iter().collect();
    genes.sort_by_key(|b| std::cmp::Reverse(b.1 .0 + b.1 .1));
    for (gene, (p, f)) in genes.iter().take(10) {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        let desc = fixture
            .gene_descriptions
            .get(*gene)
            .map(|s| s.as_str())
            .unwrap_or("");
        eprintln!("  {}: {}/{} ({:.2}%) - {}", gene, p, tot, r, desc);
    }

    eprintln!("\n========================================\n");

    // Benchmark assertion
    assert!(
        pass_rate > 99.0,
        "Paraphase exhaustive pass rate should be >99%, got {:.2}%",
        pass_rate
    );
}
