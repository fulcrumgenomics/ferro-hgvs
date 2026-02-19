//! ClinVar HGVS Bulk Tests
//!
//! Tests against the comprehensive ClinVar HGVS fixtures:
//! - clinvar_hgvs_500k.json.gz: 500K stratified sample for CI
//! - clinvar_hgvs_unique.json.gz: 4.2M unique variants for extended testing

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinvarHgvsFixture {
    description: String,
    source: String,
    source_url: String,
    version: String,
    generated: String,
    total_test_cases: usize,
    summary: ClinvarHgvsSummary,
    test_cases: Vec<ClinvarHgvsCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinvarHgvsSummary {
    by_coordinate_type: HashMap<String, usize>,
    by_hgvs_type: HashMap<String, usize>,
    #[serde(default)]
    source_stats: HashMap<String, serde_json::Value>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinvarHgvsCase {
    input: String,
    #[serde(rename = "type")]
    coord_type: String,
    hgvs_type: String,
    variation_id: String,
    #[serde(default)]
    gene: Option<String>,
    #[serde(default)]
    assembly: Option<String>,
    #[serde(default)]
    protein_expr: Option<String>,
    valid: bool,
}

fn load_fixture(filename: &str) -> Option<ClinvarHgvsFixture> {
    let path = format!("tests/fixtures/bulk/{}", filename);
    if !std::path::Path::new(&path).exists() {
        return None;
    }
    let file = File::open(&path).unwrap_or_else(|e| panic!("Failed to open {}: {}", filename, e));
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    Some(serde_json::from_reader(reader).unwrap_or_else(|e| {
        panic!(
            "Failed to parse {} (is Git LFS installed?): {}",
            filename, e
        )
    }))
}

// ============================================================================
// 500K Stratified Sample Tests
// ============================================================================

#[test]
fn test_clinvar_hgvs_500k_benchmark() {
    let fixture = match load_fixture("clinvar_hgvs_500k.json.gz") {
        Some(f) => f,
        None => {
            eprintln!("Skipping: clinvar_hgvs_500k.json.gz not found");
            return;
        }
    };

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS 500K Benchmark");
    eprintln!("========================================");
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    let mut passed = 0;
    let mut by_type: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
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
    eprintln!("  Passed: {}/{} ({:.1}%)", passed, total, pass_rate);

    eprintln!("\nBy coordinate type:");
    for (t, (p, f)) in by_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", t, p, tot, r);
    }

    eprintln!("\n========================================\n");

    // Benchmark assertions
    assert!(
        elapsed.as_secs() < 60,
        "500K should parse in under 60 seconds, took {:.1}s",
        elapsed.as_secs_f64()
    );
}

#[test]
fn test_clinvar_hgvs_500k_by_hgvs_type() {
    let fixture = match load_fixture("clinvar_hgvs_500k.json.gz") {
        Some(f) => f,
        None => return,
    };

    let mut by_hgvs_type: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        let entry = by_hgvs_type.entry(case.hgvs_type.clone()).or_insert((0, 0));

        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nClinVar HGVS 500K - By HGVS Type:");
    let mut types: Vec<_> = by_hgvs_type.iter().collect();
    types.sort_by(|a, b| (b.1 .0 + b.1 .1).cmp(&(a.1 .0 + a.1 .1)));

    for (hgvs_type, (passed, failed)) in types.iter() {
        let total = passed + failed;
        let rate = (*passed as f64 / total as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", hgvs_type, passed, total, rate);
    }
}

#[test]
fn test_clinvar_hgvs_500k_sample_failures() {
    let fixture = match load_fixture("clinvar_hgvs_500k.json.gz") {
        Some(f) => f,
        None => return,
    };

    let mut failures: Vec<&ClinvarHgvsCase> = Vec::new();

    for case in &fixture.test_cases {
        if parse_hgvs(&case.input).is_err() && failures.len() < 20 {
            failures.push(case);
        }
    }

    eprintln!("\nSample of parse failures (first 20):");
    for case in failures {
        eprintln!("  {} [{}]", case.input, case.coord_type);
    }
}

// ============================================================================
// 4.2M Unique Variants Tests (Extended)
// ============================================================================

#[test]
fn test_clinvar_hgvs_unique_benchmark() {
    let fixture = match load_fixture("clinvar_hgvs_unique.json.gz") {
        Some(f) => f,
        None => {
            eprintln!("Skipping: clinvar_hgvs_unique.json.gz not found");
            return;
        }
    };

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS Unique Variants Benchmark");
    eprintln!("========================================");
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
    eprintln!("  Passed: {}/{} ({:.1}%)", passed, total, pass_rate);

    eprintln!("\nBy coordinate type:");
    for (t, (p, f)) in by_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", t, p, tot, r);
    }

    eprintln!("\n========================================\n");
}
