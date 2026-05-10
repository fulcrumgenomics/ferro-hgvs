//! GIAB CMRG Exhaustive Tests
//!
//! Tests against ALL ClinVar variants for the 273 CMRG genes.
//! This is the exhaustive version with 4.8M variants.
//!
//! Run with: cargo test --release --test cmrg_exhaustive_tests -- --nocapture

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

// Slim deserialization shape: only the fields the test reads are
// captured. Serde uses `IgnoredAny` for all other keys, which skips
// them without materializing Strings or sub-Values. Mostly a
// readability win (drops `#[allow(dead_code)]` and giant `Debug`
// derives); the perf gain over the previous shape is small (~1-2s on
// the 4.8M-case fixture) because `IgnoredAny` still walks the JSON
// token stream.
#[derive(Deserialize)]
struct CmrgFixture {
    total_cmrg_genes: usize,
    genes_with_variants: usize,
    test_cases: Vec<CmrgCase>,
}

#[derive(Deserialize)]
struct CmrgCase {
    input: String,
    #[serde(rename = "type")]
    coord_type: String,
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
    // Parse all cases (parallel with rayon when available — yields a ~4x
    // wall-clock speedup on a 4 vCPU CI runner). Aggregate the per-coord-
    // type counts serially after the parallel section; the aggregation
    // pass over a Vec<bool> is microseconds.
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
    for (case, ok) in fixture.test_cases.iter().zip(oks.iter()) {
        let entry = by_type.entry(case.coord_type.clone()).or_insert((0, 0));
        if *ok {
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
