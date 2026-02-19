//! GIAB CMRG Exhaustive Tests
//!
//! Tests against ALL ClinVar variants for the 273 CMRG genes.
//! This is the exhaustive version with 4.8M variants.
//!
//! Run with: cargo test --release --test cmrg_exhaustive_tests -- --nocapture

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct CmrgFixture {
    description: String,
    source: String,
    source_url: String,
    reference: String,
    version: String,
    generated: String,
    total_cmrg_genes: usize,
    genes_with_variants: usize,
    genes_missing: usize,
    total_test_cases: usize,
    summary: serde_json::Value,
    test_cases: Vec<CmrgCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct CmrgCase {
    input: String,
    #[serde(rename = "type")]
    coord_type: String,
    hgvs_type: String,
    variation_id: String,
    gene: String,
    #[serde(default)]
    assembly: Option<String>,
    valid: bool,
}

fn load_fixture() -> Option<CmrgFixture> {
    let path = "tests/fixtures/validation/cmrg_genes_exhaustive.json.gz";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    let file = File::open(path).expect("Failed to open cmrg_genes_exhaustive.json.gz");
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    Some(serde_json::from_reader(reader).expect("Failed to parse cmrg_genes_exhaustive.json.gz"))
}

#[test]
fn test_cmrg_exhaustive_benchmark() {
    let fixture = match load_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping: cmrg_genes_exhaustive.json not found");
            return;
        }
    };

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("GIAB CMRG Exhaustive Benchmark");
    eprintln!("========================================");
    eprintln!("Total CMRG genes: {}", fixture.total_cmrg_genes);
    eprintln!("Genes with variants: {}", fixture.genes_with_variants);
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    let mut passed = 0;
    let mut by_type: HashMap<String, (usize, usize)> = HashMap::new();

    for (i, case) in fixture.test_cases.iter().enumerate() {
        if i % 500000 == 0 && i > 0 {
            let elapsed = start.elapsed();
            let rate = i as f64 / elapsed.as_secs_f64();
            eprintln!("  Progress: {}/{} ({:.0}/sec)", i, total, rate);
        }

        let result = parse_hgvs(&case.input);
        let entry = by_type.entry(case.coord_type.clone()).or_insert((0, 0));

        if result.is_ok() {
            passed += 1;
            entry.0 += 1;
        } else {
            entry.1 += 1;
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
    for (t, (p, f)) in by_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.2}%)", t, p, tot, r);
    }

    eprintln!("\n========================================\n");

    // Benchmark assertion
    assert!(
        pass_rate > 99.0,
        "CMRG exhaustive pass rate should be >99%, got {:.2}%",
        pass_rate
    );
}
