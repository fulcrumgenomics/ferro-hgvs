//! ClinVar HGVS Bulk Tests
//!
//! Tests against the comprehensive ClinVar HGVS fixtures:
//! - clinvar_hgvs_500k.json.gz: 500K stratified sample for CI
//! - clinvar_hgvs_unique.json.gz: 4.2M unique variants (broadest coverage)

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

mod common;
use common::failure_expectations::{enforce, FixtureCheck};

const CLINVAR_500K_FAILURE_EXPECTATIONS_PATH: &str =
    "tests/fixtures/bulk/clinvar_hgvs_500k_failure_expectations.json";
const CLINVAR_UNIQUE_FAILURE_EXPECTATIONS_PATH: &str =
    "tests/fixtures/bulk/clinvar_hgvs_unique_failure_expectations.json";

// Slim deserialization shape: borrowed `&'a str` against the
// decompressed JSON buffer. See cmrg_exhaustive_tests for rationale.
#[derive(Deserialize)]
struct ClinvarHgvsFixture<'a> {
    #[serde(borrow)]
    test_cases: Vec<ClinvarHgvsCase<'a>>,
}

#[derive(Deserialize)]
struct ClinvarHgvsCase<'a> {
    #[serde(borrow)]
    input: &'a str,
    #[serde(rename = "type", borrow)]
    coord_type: &'a str,
    #[serde(borrow)]
    hgvs_type: &'a str,
}

fn load_fixture_bytes(filename: &str) -> Option<Vec<u8>> {
    let path = format!("tests/fixtures/bulk/{}", filename);
    if !std::path::Path::new(&path).exists() {
        return None;
    }
    // See cmrg_exhaustive_tests::load_fixture_bytes for why we
    // decompress to a Vec and use `from_slice`.
    let file = File::open(&path).unwrap_or_else(|e| panic!("Failed to open {}: {}", filename, e));
    let mut buf = Vec::new();
    GzDecoder::new(file)
        .read_to_end(&mut buf)
        .unwrap_or_else(|e| panic!("Failed to decompress {}: {}", filename, e));
    Some(buf)
}

// ============================================================================
// 500K Stratified Sample
// ============================================================================

/// Single pass over the 500K stratified ClinVar fixture. Produces timing,
/// per-coord-type and per-hgvs-type breakdowns, and a sample of failures
/// for diagnostics; enforces a per-input failure-expectations snapshot
/// (see `tests/common/failure_expectations.rs`) plus a 60s throughput
/// floor.
///
/// Regenerate the failure-expectations snapshot when the parser or
/// fixture changes:
///
///   UPDATE_FAILURE_EXPECTATIONS=1 \
///     cargo nextest run --features dev test_clinvar_hgvs_500k_benchmark
#[test]
fn test_clinvar_hgvs_500k_benchmark() {
    let buf = match load_fixture_bytes("clinvar_hgvs_500k.json.gz") {
        Some(b) => b,
        None => {
            eprintln!("Skipping: clinvar_hgvs_500k.json.gz not found");
            return;
        }
    };
    let fixture: ClinvarHgvsFixture<'_> = serde_json::from_slice(&buf)
        .expect("Failed to parse clinvar_hgvs_500k.json.gz (is Git LFS installed?)");

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS 500K Benchmark");
    eprintln!("========================================");
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    #[cfg(feature = "parallel")]
    let case_failures: Vec<(&str, String)> = fixture
        .test_cases
        .par_iter()
        .filter_map(|case| {
            parse_hgvs(case.input)
                .err()
                .map(|e| (case.input, e.to_string()))
        })
        .collect();
    #[cfg(not(feature = "parallel"))]
    let case_failures: Vec<(&str, String)> = fixture
        .test_cases
        .iter()
        .filter_map(|case| {
            parse_hgvs(case.input)
                .err()
                .map(|e| (case.input, e.to_string()))
        })
        .collect();

    let failures: BTreeMap<&str, String> = case_failures.into_iter().collect();
    let passed = total - failures.len();
    let elapsed = start.elapsed();
    let rate = total as f64 / elapsed.as_secs_f64();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    let mut by_coord_type: HashMap<&str, (usize, usize)> = HashMap::new();
    let mut by_hgvs_type: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    let mut sample_failures: Vec<&ClinvarHgvsCase<'_>> = Vec::new();

    for case in &fixture.test_cases {
        let coord_entry = by_coord_type.entry(case.coord_type).or_insert((0, 0));
        let hgvs_entry = by_hgvs_type.entry(case.hgvs_type).or_insert((0, 0));
        if failures.contains_key(case.input) {
            coord_entry.1 += 1;
            hgvs_entry.1 += 1;
            if sample_failures.len() < 20 {
                sample_failures.push(case);
            }
        } else {
            coord_entry.0 += 1;
            hgvs_entry.0 += 1;
        }
    }

    eprintln!("\nPerformance:");
    eprintln!("  Time: {:.2}s", elapsed.as_secs_f64());
    eprintln!("  Rate: {:.0} variants/sec", rate);

    eprintln!("\nResults:");
    eprintln!("  Passed: {}/{} ({:.1}%)", passed, total, pass_rate);

    eprintln!("\nBy coordinate type:");
    for (t, (p, f)) in by_coord_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", t, p, tot, r);
    }

    eprintln!("\nBy HGVS type:");
    let mut types: Vec<_> = by_hgvs_type.iter().collect();
    types.sort_by_key(|b| std::cmp::Reverse(b.1 .0 + b.1 .1));
    for (hgvs_type, (p, f)) in types.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", hgvs_type, p, tot, r);
    }

    eprintln!("\nSample of parse failures (first 20):");
    for case in &sample_failures {
        eprintln!("  {} [{}]", case.input, case.coord_type);
    }

    eprintln!("\n========================================\n");

    // Throughput floor.
    assert!(
        elapsed.as_secs() < 60,
        "500K should parse in under 60 seconds, took {:.1}s",
        elapsed.as_secs_f64()
    );

    enforce(
        Path::new(CLINVAR_500K_FAILURE_EXPECTATIONS_PATH),
        "UPDATE_FAILURE_EXPECTATIONS",
        FixtureCheck {
            total_inputs: total,
            failures,
        },
    );
}

// ============================================================================
// 4.2M Unique Variants
// ============================================================================

/// Exhaustive parse over the 4.2M unique ClinVar HGVS strings — the broadest
/// single fixture, including the long tail that the 500K stratified sample
/// drops and inputs outside the CMRG/paraphase gene scopes. Enforces a
/// per-input failure-expectations snapshot (see
/// `tests/common/failure_expectations.rs`).
#[test]
fn test_clinvar_hgvs_unique_benchmark() {
    let buf = match load_fixture_bytes("clinvar_hgvs_unique.json.gz") {
        Some(b) => b,
        None => {
            eprintln!("Skipping: clinvar_hgvs_unique.json.gz not found");
            return;
        }
    };
    let fixture: ClinvarHgvsFixture<'_> = serde_json::from_slice(&buf)
        .expect("Failed to parse clinvar_hgvs_unique.json.gz (is Git LFS installed?)");

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS Unique Variants Benchmark");
    eprintln!("========================================");
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    #[cfg(feature = "parallel")]
    let case_failures: Vec<(&str, String)> = fixture
        .test_cases
        .par_iter()
        .filter_map(|case| {
            parse_hgvs(case.input)
                .err()
                .map(|e| (case.input, e.to_string()))
        })
        .collect();
    #[cfg(not(feature = "parallel"))]
    let case_failures: Vec<(&str, String)> = fixture
        .test_cases
        .iter()
        .filter_map(|case| {
            parse_hgvs(case.input)
                .err()
                .map(|e| (case.input, e.to_string()))
        })
        .collect();

    let failures: BTreeMap<&str, String> = case_failures.into_iter().collect();
    let passed = total - failures.len();
    let elapsed = start.elapsed();
    let rate = total as f64 / elapsed.as_secs_f64();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    let mut by_type: HashMap<&str, (usize, usize)> = HashMap::new();
    for case in &fixture.test_cases {
        let entry = by_type.entry(case.coord_type).or_insert((0, 0));
        if failures.contains_key(case.input) {
            entry.1 += 1;
        } else {
            entry.0 += 1;
        }
    }

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

    enforce(
        Path::new(CLINVAR_UNIQUE_FAILURE_EXPECTATIONS_PATH),
        "UPDATE_FAILURE_EXPECTATIONS",
        FixtureCheck {
            total_inputs: total,
            failures,
        },
    );
}
