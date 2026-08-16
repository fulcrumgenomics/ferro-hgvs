//! ClinVar HGVS Bulk Tests
//!
//! Tests against the comprehensive ClinVar HGVS fixtures:
//! - clinvar_hgvs_500k.json.gz: 500K stratified sample for CI
//! - clinvar_hgvs_unique.json.gz: 4.2M unique variants (broadest coverage)

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Read;
use std::time::Instant;

use crate::common::failure_expectations::{enforce, FixtureCheck};
use crate::common::fixture_gen;

const CLINVAR_500K_FAILURE_EXPECTATIONS_PATH: &str =
    "tests/fixtures/bulk/clinvar_hgvs_500k_failure_expectations.json";
const CLINVAR_UNIQUE_FAILURE_EXPECTATIONS_PATH: &str =
    "tests/fixtures/bulk/clinvar_hgvs_unique_failure_expectations.json";
const CLINVAR_UNIQUE_ROUNDTRIP_EXPECTATIONS_PATH: &str =
    "tests/fixtures/bulk/clinvar_hgvs_unique_roundtrip_expectations.json";

/// A generous floor on how many inputs the round-trip guard must actually
/// exercise. It is a **non-vacuity guard**, not the property under test: the
/// property is that `parse → format → parse → format` agrees, and the snapshot
/// below is what pins it. This number only ensures a future breakage that made
/// most of the corpus unparseable (so nothing reached the round trip) cannot
/// masquerade as a clean pass. The fixture currently offers 4,188,224 parseable
/// inputs; the floor sits well below that so ordinary corpus/parser drift does
/// not trip it, while a collapse to a handful would.
const CLINVAR_UNIQUE_ROUNDTRIP_MIN_INPUTS: usize = 4_000_000;

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
    // Absent means "skip" locally and "fail" under `FERRO_REQUIRE_BULK_FIXTURES`,
    // which CI sets — see `common::bulk_fixtures`.
    let path = crate::common::bulk_fixtures::present_or_skip(&path)?;
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
        // `load_fixture_bytes` has already reported the skip, or panicked
        // under `FERRO_REQUIRE_BULK_FIXTURES`.
        None => return,
    };
    let fixture: ClinvarHgvsFixture<'_> =
        serde_json::from_slice(&buf).expect("Failed to parse clinvar_hgvs_500k.json.gz");

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS 500K Benchmark");
    eprintln!("========================================");
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    let case_failures: Vec<(&str, String)> = fixture
        .test_cases
        .par_iter()
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
        // Resolved against the workspace root, not the cwd: nextest sets the cwd
        // to the PACKAGE root, and this module is compiled into
        // `ferro-hgvs-soak-tests` as well as into `it`, where that root is
        // `tests-soak/`. A bare relative path resolves there and the snapshot is
        // reported missing — loudly, but in the one binary that cannot be run
        // from the workspace root. Same reason `common::bulk_fixtures` and
        // `common::fixture_gen` resolve their own paths this way.
        &fixture_gen::fixture_path(CLINVAR_500K_FAILURE_EXPECTATIONS_PATH),
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
        // `load_fixture_bytes` has already reported the skip, or panicked
        // under `FERRO_REQUIRE_BULK_FIXTURES`.
        None => return,
    };
    let fixture: ClinvarHgvsFixture<'_> =
        serde_json::from_slice(&buf).expect("Failed to parse clinvar_hgvs_unique.json.gz");

    let total = fixture.test_cases.len();
    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS Unique Variants Benchmark");
    eprintln!("========================================");
    eprintln!("Total test cases: {}", total);

    let start = Instant::now();
    let case_failures: Vec<(&str, String)> = fixture
        .test_cases
        .par_iter()
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
        // See the note on the 500k call above: workspace-root-relative, because
        // this module now also compiles into the soak driver.
        &fixture_gen::fixture_path(CLINVAR_UNIQUE_FAILURE_EXPECTATIONS_PATH),
        "UPDATE_FAILURE_EXPECTATIONS",
        FixtureCheck {
            total_inputs: total,
            failures,
        },
    );
}

// ============================================================================
// Round trip over the 4.2M unique variants
// ============================================================================

/// Round-trip guard over the 4.2M unique ClinVar HGVS strings (issue #1859):
/// for every input that parses, `format!` it, re-parse the rendering, and
/// require the two renderings to agree — i.e. `parse → format → parse → format`
/// is stable.
///
/// **Why the sibling benchmarks do not cover this.**
/// `test_clinvar_hgvs_unique_benchmark` (and the 500K one) read the same corpora
/// but only *parse*: they assert a parse-failure snapshot and, for the 500K, a
/// throughput floor. Neither ever calls `format!` on a parsed variant, so a
/// formatter that emits something the parser rejects — or accepts as a
/// *different* variant — is invisible to them. The idempotency suite does the
/// round trip, but only over the curated spec/edge-case corpora, not the long
/// tail of what ClinVar actually stores. This test is the intersection that was
/// nominally covered by a test (`test_clinvar_sample_parsing_idempotency`) which
/// read a corpus that never existed and so only ever took its skip branch.
///
/// **The unique corpus, not the 500K stratified sample, on purpose.** Measured
/// while filing #1859: the 500K sample round-trips with zero disagreements,
/// while the broadest fixture surfaces exactly one — a malformed ClinVar input
/// carrying a doubled `ins` token, whose rendering re-wraps as `ins[ins…]` and
/// no longer parses. That single residue is pinned below; it is the long tail
/// the stratified sample drops.
///
/// The residue is tracked per-input via the same failure-expectations snapshot
/// framework the parse-only benchmarks use (`tests/common/failure_expectations.rs`),
/// rather than a percentage floor: a new disagreement is a loud regression, and
/// a previously-disagreeing input that starts round-tripping is flagged as an
/// improvement to re-bless. Regenerate the snapshot after an intended change:
///
///   UPDATE_FAILURE_EXPECTATIONS=1 \
///     cargo nextest run --features dev test_clinvar_hgvs_unique_roundtrip
///
/// Absence of the fixture skips locally and fails under
/// `FERRO_REQUIRE_BULK_FIXTURES` — see `common::bulk_fixtures`.
#[test]
fn test_clinvar_hgvs_unique_roundtrip() {
    let buf = match load_fixture_bytes("clinvar_hgvs_unique.json.gz") {
        Some(b) => b,
        // `load_fixture_bytes` has already reported the skip, or panicked
        // under `FERRO_REQUIRE_BULK_FIXTURES`.
        None => return,
    };
    let fixture: ClinvarHgvsFixture<'_> =
        serde_json::from_slice(&buf).expect("Failed to parse clinvar_hgvs_unique.json.gz");

    let start = Instant::now();

    // One parallel pass: count the inputs that parse (the round-trip domain —
    // parse failures are the sibling benchmark's job, not this one) and collect
    // every disagreement. `fold`/`reduce` keeps only the tiny disagreement list
    // in memory rather than one entry per input.
    let report = fixture
        .test_cases
        .par_iter()
        .fold(RoundTripReport::default, |mut acc, case| {
            let Ok(v1) = parse_hgvs(case.input) else {
                return acc; // does not parse: not a round-trip case
            };
            acc.roundtripped += 1;
            let rendered = format!("{}", v1);
            match parse_hgvs(&rendered) {
                Err(e) => acc
                    .disagreements
                    .push((case.input, format!("re-parse of rendering failed: {e}"))),
                Ok(v2) => {
                    let reparsed = format!("{}", v2);
                    if rendered != reparsed {
                        acc.disagreements.push((
                            case.input,
                            format!("rendering not stable: '{rendered}' vs '{reparsed}'"),
                        ));
                    }
                }
            }
            acc
        })
        .reduce(RoundTripReport::default, RoundTripReport::merge);

    let elapsed = start.elapsed();
    let failures: BTreeMap<&str, String> = report.disagreements.into_iter().collect();

    eprintln!("\n========================================");
    eprintln!("ClinVar HGVS Unique Round-Trip");
    eprintln!("========================================");
    eprintln!("Total cases:        {}", fixture.test_cases.len());
    eprintln!("Round-tripped:      {}", report.roundtripped);
    eprintln!("Disagreements:      {}", failures.len());
    eprintln!("Time:               {:.2}s", elapsed.as_secs_f64());
    eprintln!("========================================\n");

    // Non-vacuity: the guard is only meaningful if it actually ran the round
    // trip over the corpus. Checked before `enforce` so a corpus that collapsed
    // to nothing fails here loudly rather than passing with an empty snapshot.
    assert!(
        report.roundtripped >= CLINVAR_UNIQUE_ROUNDTRIP_MIN_INPUTS,
        "round-trip exercised only {} inputs, expected at least {} — did parsing collapse?",
        report.roundtripped,
        CLINVAR_UNIQUE_ROUNDTRIP_MIN_INPUTS
    );

    enforce(
        // Workspace-root-relative, for the same reason as the benchmarks above:
        // this module also compiles into the soak driver.
        &fixture_gen::fixture_path(CLINVAR_UNIQUE_ROUNDTRIP_EXPECTATIONS_PATH),
        "UPDATE_FAILURE_EXPECTATIONS",
        FixtureCheck {
            // The round-trip population is the inputs that parsed, so that is the
            // denominator the zero-inputs bless guard keys on.
            total_inputs: report.roundtripped,
            failures,
        },
    );
}

/// Accumulator for the round-trip pass: how many inputs were round-tripped and
/// which of them disagreed (`input -> reason`).
#[derive(Default)]
struct RoundTripReport<'a> {
    roundtripped: usize,
    disagreements: Vec<(&'a str, String)>,
}

impl<'a> RoundTripReport<'a> {
    fn merge(mut self, mut other: Self) -> Self {
        self.roundtripped += other.roundtripped;
        self.disagreements.append(&mut other.disagreements);
        self
    }
}
