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
use std::io::Read;
use std::time::Instant;

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
    // Parse all cases (parallel with rayon when available — yields a ~4x
    // wall-clock speedup on a 4 vCPU CI runner). Aggregate the per-coord-
    // type counts serially after the parallel section.
    #[cfg(feature = "parallel")]
    let oks: Vec<bool> = fixture
        .test_cases
        .par_iter()
        .map(|case| parse_hgvs(case.input).is_ok())
        .collect();
    #[cfg(not(feature = "parallel"))]
    let oks: Vec<bool> = fixture
        .test_cases
        .iter()
        .map(|case| parse_hgvs(case.input).is_ok())
        .collect();

    // Per-coord-type tally with `&str` keys borrowed from the buffer:
    // there are ~3 distinct coord_types, so this stays a tiny map.
    let mut passed = 0usize;
    let mut by_type: HashMap<&str, (usize, usize)> = HashMap::new();
    for (case, ok) in fixture.test_cases.iter().zip(oks.iter()) {
        let entry = by_type.entry(case.coord_type).or_insert((0, 0));
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
