//! ClinVar HGVS Bulk Tests
//!
//! Tests against the comprehensive ClinVar HGVS fixtures:
//! - clinvar_hgvs_500k.json.gz: 500K stratified sample for CI
//! - clinvar_hgvs_unique.json.gz: 4.2M unique variants for extended testing

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

// Slim deserialization shape: drop fields the test never reads.
// Mostly a readability win; serde's `IgnoredAny` still walks unread
// JSON tokens, so the perf gain is small.
#[derive(Deserialize)]
struct ClinvarHgvsFixture {
    test_cases: Vec<ClinvarHgvsCase>,
}

#[derive(Deserialize)]
struct ClinvarHgvsCase {
    input: String,
    #[serde(rename = "type")]
    coord_type: String,
    hgvs_type: String,
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

/// Path to the per-HGVS-type passing-count snapshot. Acts as a strict floor:
/// the test asserts current passed[type] >= snapshot[type] for every type
/// recorded in the snapshot, so per-type regressions hidden by aggregate
/// pass-rate can't slip through.
const PER_TYPE_SNAPSHOT_PATH: &str =
    "tests/fixtures/bulk/clinvar_hgvs_500k_passing_by_hgvs_type.json";

/// Single pass over the 500K stratified ClinVar fixture that produces all
/// three diagnostics the original three tests produced (timing, per-type
/// breakdown, sample failures) plus a per-HGVS-type strict-floor regression
/// guard.
///
/// Regenerate the snapshot when ClinVar fixture changes or when a parser
/// improvement intentionally lifts the per-type counts:
///
///   UPDATE_CLINVAR_500K_SNAPSHOT=1 \
///     cargo nextest run --features dev test_clinvar_hgvs_500k_benchmark
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
    // Parse all cases (parallel with rayon when available). Aggregate
    // per-type counts and the first-20 failure sample serially over the
    // resulting Vec<bool>; the aggregation pass is microseconds.
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
    let mut by_coord_type: HashMap<String, (usize, usize)> = HashMap::new();
    let mut by_hgvs_type: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    let mut sample_failures: Vec<&ClinvarHgvsCase> = Vec::new();

    for (case, ok) in fixture.test_cases.iter().zip(oks.iter()) {
        let coord_entry = by_coord_type
            .entry(case.coord_type.clone())
            .or_insert((0, 0));
        let hgvs_entry = by_hgvs_type.entry(case.hgvs_type.clone()).or_insert((0, 0));

        if *ok {
            passed += 1;
            coord_entry.0 += 1;
            hgvs_entry.0 += 1;
        } else {
            coord_entry.1 += 1;
            hgvs_entry.1 += 1;
            if sample_failures.len() < 20 {
                sample_failures.push(case);
            }
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

    // Existing throughput floor.
    assert!(
        elapsed.as_secs() < 60,
        "500K should parse in under 60 seconds, took {:.1}s",
        elapsed.as_secs_f64()
    );

    // Per-HGVS-type strict-floor regression guard.
    let observed_passing: BTreeMap<String, usize> = by_hgvs_type
        .iter()
        .map(|(t, (p, _))| (t.clone(), *p))
        .collect();
    enforce_per_type_floor(Path::new(PER_TYPE_SNAPSHOT_PATH), &observed_passing);
}

/// Strict-floor check: every HGVS type recorded in the snapshot must have
/// `observed_passing[type] >= snapshot[type]`. New types in the observed map
/// are ignored (they'll be picked up the next time the snapshot is
/// regenerated). Set `UPDATE_CLINVAR_500K_SNAPSHOT=1` to overwrite the
/// snapshot from the current run.
fn enforce_per_type_floor(snapshot_path: &Path, observed: &BTreeMap<String, usize>) {
    if std::env::var_os("UPDATE_CLINVAR_500K_SNAPSHOT").is_some() {
        let json =
            serde_json::to_string_pretty(observed).expect("serialize per-type passing snapshot");
        std::fs::write(snapshot_path, json + "\n")
            .unwrap_or_else(|e| panic!("write {}: {}", snapshot_path.display(), e));
        eprintln!("Updated snapshot: {}", snapshot_path.display());
        return;
    }

    let raw = match std::fs::read_to_string(snapshot_path) {
        Ok(s) => s,
        Err(e) => panic!(
            "{}: {}\nGenerate it with: UPDATE_CLINVAR_500K_SNAPSHOT=1 cargo nextest run \
             --features dev test_clinvar_hgvs_500k_benchmark",
            snapshot_path.display(),
            e
        ),
    };
    let snapshot: BTreeMap<String, usize> = serde_json::from_str(&raw)
        .unwrap_or_else(|e| panic!("parse {}: {}", snapshot_path.display(), e));

    let mut regressions: Vec<(String, usize, usize)> = Vec::new();
    for (hgvs_type, &expected) in &snapshot {
        let actual = observed.get(hgvs_type).copied().unwrap_or(0);
        if actual < expected {
            regressions.push((hgvs_type.clone(), expected, actual));
        }
    }

    if !regressions.is_empty() {
        let mut msg = String::from(
            "ClinVar 500K per-HGVS-type pass count regressed below snapshot floor.\n\
             If this drop is intentional (parser change OR fixture refresh that lowers \
             counts), regenerate the snapshot:\n\
             \n  UPDATE_CLINVAR_500K_SNAPSHOT=1 cargo nextest run --features dev \
             test_clinvar_hgvs_500k_benchmark\n\nRegressions:\n",
        );
        for (t, expected, actual) in &regressions {
            msg.push_str(&format!(
                "  {}: expected >= {}, got {} (delta {})\n",
                t,
                expected,
                actual,
                *actual as i64 - *expected as i64
            ));
        }
        panic!("{}", msg);
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
