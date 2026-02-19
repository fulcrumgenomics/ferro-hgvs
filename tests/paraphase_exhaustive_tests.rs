//! Paraphase Exhaustive Tests
//!
//! Tests against ALL ClinVar variants for Paraphase-supported genes.
//! This is the exhaustive version with 435K variants from 35 segdup genes.
//!
//! Run with: cargo test --release --test paraphase_exhaustive_tests -- --nocapture

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParaphaseFixture {
    description: String,
    source: String,
    source_url: String,
    version: String,
    generated: String,
    total_target_genes: usize,
    genes_with_variants: usize,
    genes_missing: usize,
    total_test_cases: usize,
    gene_descriptions: HashMap<String, String>,
    summary: serde_json::Value,
    test_cases: Vec<ParaphaseCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParaphaseCase {
    input: String,
    #[serde(rename = "type")]
    coord_type: String,
    hgvs_type: String,
    variation_id: String,
    gene: String,
    challenge: String,
    #[serde(default)]
    assembly: Option<String>,
    valid: bool,
}

fn load_fixture() -> Option<ParaphaseFixture> {
    let path = "tests/fixtures/validation/paraphase_genes_exhaustive.json.gz";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    let file = File::open(path).expect("Failed to open paraphase_genes_exhaustive.json.gz");
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    Some(
        serde_json::from_reader(reader)
            .expect("Failed to parse paraphase_genes_exhaustive.json.gz"),
    )
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
    let mut passed = 0;
    let mut by_type: HashMap<String, (usize, usize)> = HashMap::new();
    let mut by_gene: HashMap<String, (usize, usize)> = HashMap::new();

    for (i, case) in fixture.test_cases.iter().enumerate() {
        if i % 100000 == 0 && i > 0 {
            let elapsed = start.elapsed();
            let rate = i as f64 / elapsed.as_secs_f64();
            eprintln!("  Progress: {}/{} ({:.0}/sec)", i, total, rate);
        }

        let result = parse_hgvs(&case.input);
        let type_entry = by_type.entry(case.coord_type.clone()).or_insert((0, 0));
        let gene_entry = by_gene.entry(case.gene.clone()).or_insert((0, 0));

        if result.is_ok() {
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
    types.sort_by(|a, b| (b.1 .0 + b.1 .1).cmp(&(a.1 .0 + a.1 .1)));
    for (t, (p, f)) in types.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.2}%)", t, p, tot, r);
    }

    eprintln!("\nTop 10 genes by variant count:");
    let mut genes: Vec<_> = by_gene.iter().collect();
    genes.sort_by(|a, b| (b.1 .0 + b.1 .1).cmp(&(a.1 .0 + a.1 .1)));
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
