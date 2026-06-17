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
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

use crate::common::failure_expectations::{enforce, FixtureCheck};

const FAILURE_EXPECTATIONS_PATH: &str =
    "tests/fixtures/validation/cmrg_genes_failure_expectations.json";

// Slim deserialization shape: only the fields the test reads are
// captured. `&'a str` borrows directly into the decompressed JSON
// buffer (no per-string allocation; serde_json hands back slices
// pointing into the buffer). Combined with `IgnoredAny` for unread
// keys this avoids ~10M `String` allocations on the cmrg fixture.
#[derive(Deserialize)]
struct CmrgFixture<'a> {
    total_cmrg_genes: usize,
    genes_with_variants: usize,
    #[serde(borrow)]
    test_cases: Vec<CmrgCase<'a>>,
}

#[derive(Deserialize)]
struct CmrgCase<'a> {
    #[serde(borrow)]
    input: &'a str,
    #[serde(rename = "type", borrow)]
    coord_type: &'a str,
}

fn load_fixture_bytes() -> Option<Vec<u8>> {
    let path = "tests/fixtures/validation/cmrg_genes_exhaustive.json.gz";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    // Decompress to an in-memory Vec, then deser via `from_slice` rather
    // than `from_reader` over a streaming gzip pipeline. Empirically ~5x
    // faster on this fixture (34s -> 7s in debug mode) — `from_slice`
    // works against a contiguous byte slice and supports zero-copy
    // borrowed `&str` deserialization. The transient ~1 GB buffer is
    // fine on the CI runner (16 GB RAM).
    let file = File::open(path).expect("Failed to open cmrg_genes_exhaustive.json.gz");
    let mut buf = Vec::new();
    GzDecoder::new(file)
        .read_to_end(&mut buf)
        .expect("Failed to decompress cmrg_genes_exhaustive.json.gz");
    Some(buf)
}

#[test]
fn test_cmrg_exhaustive_benchmark() {
    let buf = match load_fixture_bytes() {
        Some(b) => b,
        None => {
            eprintln!("Skipping: cmrg_genes_exhaustive.json not found");
            return;
        }
    };
    let fixture: CmrgFixture<'_> =
        serde_json::from_slice(&buf).expect("Failed to parse cmrg_genes_exhaustive.json.gz");

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("GIAB CMRG Exhaustive Benchmark");
    eprintln!("========================================");
    eprintln!("Total CMRG genes: {}", fixture.total_cmrg_genes);
    eprintln!("Genes with variants: {}", fixture.genes_with_variants);
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    // Parse all cases (parallel with rayon when available). For each
    // failure capture (input, error_string) so the per-input expectations
    // framework can categorize it.
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

    // Per-coord-type tally for diagnostic eprintln only.
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
    eprintln!("  Passed: {}/{} ({:.2}%)", passed, total, pass_rate);

    eprintln!("\nBy coordinate type:");
    for (t, (p, f)) in by_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.2}%)", t, p, tot, r);
    }

    eprintln!("\n========================================\n");

    enforce(
        Path::new(FAILURE_EXPECTATIONS_PATH),
        "UPDATE_FAILURE_EXPECTATIONS",
        FixtureCheck {
            total_inputs: total,
            failures,
        },
    );
}
