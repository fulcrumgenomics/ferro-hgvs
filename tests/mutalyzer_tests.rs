//! Tests using mutalyzer/hgvs-parser test data
//!
//! These test cases are extracted from the mutalyzer project
//! to ensure compatibility with their HGVS parsing.

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct TestFixture {
    description: String,
    source: String,
    parsing_tests: Vec<ParsingTestCase>,
    normalization_tests: Vec<NormalizationTestCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParsingTestCase {
    input: String,
    valid: bool,
    #[serde(rename = "type")]
    variant_type: Option<String>,
    supported: bool,
    source: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct NormalizationTestCase {
    input: String,
    normalized: String,
    description: String,
}

fn load_mutalyzer_fixtures() -> TestFixture {
    let content = fs::read_to_string("tests/fixtures/normalization/mutalyzer.json")
        .expect("Failed to read mutalyzer.json");
    serde_json::from_str(&content).expect("Failed to parse mutalyzer.json")
}

#[test]
fn test_mutalyzer_valid_variants() {
    let fixtures = load_mutalyzer_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for test in fixtures.parsing_tests.iter().filter(|t| t.valid) {
        if !test.supported {
            skipped += 1;
            continue;
        }

        let result = parse_hgvs(&test.input);

        match result {
            Ok(variant) => {
                if let Some(ref expected_type) = test.variant_type {
                    let actual_type = match &variant {
                        HgvsVariant::Genome(_) => "Genome",
                        HgvsVariant::Cds(_) => "Cds",
                        HgvsVariant::Tx(_) => "Tx",
                        HgvsVariant::Rna(_) => "Rna",
                        HgvsVariant::Protein(_) => "Protein",
                        HgvsVariant::Mt(_) => "Mt",
                        HgvsVariant::Circular(_) => "Circular",
                        HgvsVariant::RnaFusion(_) => "RnaFusion",
                        HgvsVariant::Allele(_) => "Allele",
                        HgvsVariant::NullAllele => "NullAllele",
                        HgvsVariant::UnknownAllele => "UnknownAllele",
                    };

                    if actual_type == expected_type {
                        passed += 1;
                    } else {
                        failed += 1;
                        failures.push(format!(
                            "{}: expected type {}, got {}",
                            test.input, expected_type, actual_type
                        ));
                    }
                } else {
                    passed += 1;
                }
            }
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: parse error: {}", test.input, e));
            }
        }
    }

    if !failures.is_empty() {
        eprintln!("\nMutalyzer valid variant failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    println!(
        "Mutalyzer valid: {} passed, {} failed, {} skipped (unsupported)",
        passed, failed, skipped
    );

    assert!(
        failed == 0,
        "Mutalyzer valid: {} passed, {} failed",
        passed,
        failed
    );
}

#[test]
fn test_mutalyzer_invalid_variants() {
    let fixtures = load_mutalyzer_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures.parsing_tests.iter().filter(|t| !t.valid) {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(v) => {
                failed += 1;
                failures.push(format!(
                    "{}: should have failed to parse, got {:?}",
                    test.input,
                    format!("{}", v)
                ));
            }
            Err(_) => {
                passed += 1;
            }
        }
    }

    if !failures.is_empty() {
        eprintln!("\nMutalyzer invalid variant failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        failed == 0,
        "Mutalyzer invalid: {} passed, {} failed",
        passed,
        failed
    );
}

#[test]
fn test_mutalyzer_roundtrip() {
    let fixtures = load_mutalyzer_fixtures();

    for test in fixtures
        .parsing_tests
        .iter()
        .filter(|t| t.valid && t.supported)
    {
        let parsed = parse_hgvs(&test.input);
        assert!(
            parsed.is_ok(),
            "Failed to parse: {} - {:?}",
            test.input,
            parsed.err()
        );

        let variant = parsed.unwrap();
        let displayed = format!("{}", variant);

        // Verify it can be reparsed
        let reparsed = parse_hgvs(&displayed);
        assert!(
            reparsed.is_ok(),
            "Roundtrip failed: {} -> {} - {:?}",
            test.input,
            displayed,
            reparsed.err()
        );
    }
}

#[test]
fn test_mutalyzer_coverage_summary() {
    let fixtures = load_mutalyzer_fixtures();

    let total = fixtures.parsing_tests.len();
    let valid_count = fixtures.parsing_tests.iter().filter(|t| t.valid).count();
    let invalid_count = total - valid_count;
    let supported_count = fixtures
        .parsing_tests
        .iter()
        .filter(|t| t.valid && t.supported)
        .count();
    let unsupported_count = fixtures
        .parsing_tests
        .iter()
        .filter(|t| t.valid && !t.supported)
        .count();

    println!("\nMutalyzer test coverage:");
    println!("  Total parsing tests: {}", total);
    println!("  Valid variants: {}", valid_count);
    println!("    - Supported: {}", supported_count);
    println!("    - Unsupported: {}", unsupported_count);
    println!("  Invalid variants: {}", invalid_count);
    println!(
        "  Normalization tests: {}",
        fixtures.normalization_tests.len()
    );
}
