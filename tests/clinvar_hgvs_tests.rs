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
use std::io::Read;
use std::path::Path;
use std::time::Instant;

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

    let mut passed = 0usize;
    let mut by_coord_type: HashMap<&str, (usize, usize)> = HashMap::new();
    let mut by_hgvs_type: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    let mut sample_failures: Vec<&ClinvarHgvsCase<'_>> = Vec::new();

    for (case, ok) in fixture.test_cases.iter().zip(oks.iter()) {
        let coord_entry = by_coord_type.entry(case.coord_type).or_insert((0, 0));
        let hgvs_entry = by_hgvs_type.entry(case.hgvs_type).or_insert((0, 0));

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
    let observed_passing: BTreeMap<&str, usize> =
        by_hgvs_type.iter().map(|(t, (p, _))| (*t, *p)).collect();
    enforce_per_type_floor(Path::new(PER_TYPE_SNAPSHOT_PATH), &observed_passing);
}

/// Strict-floor check: every HGVS type recorded in the snapshot must have
/// `observed_passing[type] >= snapshot[type]`. Any HGVS type present in the
/// observed map but missing from the snapshot is also a hard failure — a new
/// bucket can otherwise hide a low pass-rate behind the aggregate floor until
/// someone remembers to regenerate the snapshot. Set
/// `UPDATE_CLINVAR_500K_SNAPSHOT=1` to overwrite the snapshot from the current
/// run.
fn enforce_per_type_floor(snapshot_path: &Path, observed: &BTreeMap<&str, usize>) {
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

    let mut missing_from_snapshot: Vec<&str> = observed
        .keys()
        .copied()
        .filter(|hgvs_type| !snapshot.contains_key(*hgvs_type))
        .collect();
    missing_from_snapshot.sort_unstable();
    if !missing_from_snapshot.is_empty() {
        panic!(
            "ClinVar 500K snapshot is missing HGVS types: {:?}\n\
             Regenerate it with: UPDATE_CLINVAR_500K_SNAPSHOT=1 cargo nextest run \
             --features dev test_clinvar_hgvs_500k_benchmark",
            missing_from_snapshot
        );
    }

    let mut regressions: Vec<(String, usize, usize)> = Vec::new();
    for (hgvs_type, &expected) in &snapshot {
        let actual = observed.get(hgvs_type.as_str()).copied().unwrap_or(0);
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
// Per-HGVS-type floor unit tests
// ============================================================================

#[cfg(test)]
mod per_type_floor_tests {
    use super::enforce_per_type_floor;
    use std::collections::BTreeMap;

    fn write_snapshot(
        dir: &std::path::Path,
        snapshot: &BTreeMap<&str, usize>,
    ) -> std::path::PathBuf {
        let path = dir.join("snapshot.json");
        let json = serde_json::to_string_pretty(snapshot).unwrap();
        std::fs::write(&path, json).unwrap();
        path
    }

    fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
        if let Some(s) = payload.downcast_ref::<String>() {
            s.clone()
        } else if let Some(s) = payload.downcast_ref::<&'static str>() {
            (*s).to_string()
        } else {
            String::new()
        }
    }

    /// Observed map containing an HGVS type missing from the snapshot must
    /// panic so the floor check can't silently skip new buckets.
    #[test]
    fn rejects_observed_keys_missing_from_snapshot() {
        let dir = tempfile::tempdir().unwrap();
        let snapshot: BTreeMap<&str, usize> = BTreeMap::from([("sub", 10), ("dup", 5)]);
        let snapshot_path = write_snapshot(dir.path(), &snapshot);

        let observed: BTreeMap<&str, usize> =
            BTreeMap::from([("sub", 10), ("dup", 5), ("inv", 3), ("delins", 2)]);

        let result = std::panic::catch_unwind(|| {
            enforce_per_type_floor(&snapshot_path, &observed);
        });

        let payload = result.expect_err("should panic when observed has unknown HGVS types");
        let msg = panic_message(payload);
        assert!(msg.contains("inv"), "panic should mention 'inv': {msg}");
        assert!(
            msg.contains("delins"),
            "panic should mention 'delins': {msg}"
        );
        assert!(
            msg.contains("UPDATE_CLINVAR_500K_SNAPSHOT"),
            "panic should mention regen env var: {msg}"
        );
    }

    /// Observed values >= snapshot values for every snapshot key, with no
    /// extra keys in observed, must pass.
    #[test]
    fn accepts_observed_at_or_above_snapshot() {
        let dir = tempfile::tempdir().unwrap();
        let snapshot: BTreeMap<&str, usize> = BTreeMap::from([("sub", 10), ("dup", 5)]);
        let snapshot_path = write_snapshot(dir.path(), &snapshot);

        let observed: BTreeMap<&str, usize> = BTreeMap::from([("sub", 12), ("dup", 5)]);
        enforce_per_type_floor(&snapshot_path, &observed);
    }
}

// ============================================================================
// 4.2M Unique Variants Tests (Extended)
// ============================================================================

/// Exhaustive parse over the 4.2M unique ClinVar HGVS strings — the broadest
/// single fixture, including the long tail that the 500K stratified sample
/// drops and inputs outside the CMRG/paraphase gene scopes. Asserts a >99%
/// pass-rate floor; matches the shape of `cmrg_exhaustive_tests` and
/// `paraphase_exhaustive_tests`.
///
/// A per-input expectations framework (golden snapshot of which inputs are
/// expected to fail and why) is tracked as follow-up work that will replace
/// the coarse >99% threshold with per-input regression detection.
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
    eprintln!("  Passed: {}/{} ({:.1}%)", passed, total, pass_rate);

    eprintln!("\nBy coordinate type:");
    for (t, (p, f)) in by_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", t, p, tot, r);
    }

    eprintln!("\n========================================\n");

    assert!(
        pass_rate > 99.0,
        "ClinVar unique-variants exhaustive pass rate should be >99%, got {:.2}%",
        pass_rate
    );
}
