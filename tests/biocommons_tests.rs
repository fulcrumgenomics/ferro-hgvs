//! Tests using biocommons/hgvs test data
//!
//! These test cases are extracted from the biocommons/hgvs repository
//! to ensure compatibility with the reference implementation.

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct TestFixture {
    description: String,
    source: String,
    tests: Vec<TestCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct TestCase {
    input: String,
    valid: bool,
    #[serde(rename = "type")]
    variant_type: Option<String>,
    source: String,
    #[serde(default = "default_supported")]
    supported: bool,
}

fn default_supported() -> bool {
    true
}

fn load_biocommons_fixtures() -> TestFixture {
    let content = fs::read_to_string("tests/fixtures/grammar/biocommons.json")
        .expect("Failed to read biocommons.json");
    serde_json::from_str(&content).expect("Failed to parse biocommons.json")
}

#[test]
fn test_biocommons_valid_variants() {
    let fixtures = load_biocommons_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for test in fixtures.tests.iter().filter(|t| t.valid) {
        // Skip unsupported variants
        if !test.supported {
            skipped += 1;
            continue;
        }

        let result = parse_hgvs(&test.input);

        match result {
            Ok(variant) => {
                // Verify type if specified
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
        eprintln!("\nBiocommons valid variant failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    println!(
        "Biocommons valid: {} passed, {} failed, {} skipped (unsupported)",
        passed, failed, skipped
    );

    assert!(
        failed == 0,
        "Biocommons valid: {} passed, {} failed",
        passed,
        failed
    );
}

#[test]
fn test_biocommons_invalid_variants() {
    let fixtures = load_biocommons_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures.tests.iter().filter(|t| !t.valid) {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => {
                failed += 1;
                failures.push(format!("{}: should have failed to parse", test.input));
            }
            Err(_) => {
                passed += 1;
            }
        }
    }

    if !failures.is_empty() {
        eprintln!("\nBiocommons invalid variant failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        failed == 0,
        "Biocommons invalid: {} passed, {} failed",
        passed,
        failed
    );
}

#[test]
fn test_biocommons_gauntlet_valid() {
    let fixtures = load_biocommons_fixtures();
    let gauntlet_tests: Vec<_> = fixtures
        .tests
        .iter()
        .filter(|t| t.source == "biocommons-gauntlet" && t.valid && t.supported)
        .collect();

    assert!(
        !gauntlet_tests.is_empty(),
        "No supported gauntlet tests found in fixtures"
    );

    let mut passed = 0;
    for test in gauntlet_tests {
        let result = parse_hgvs(&test.input);
        assert!(
            result.is_ok(),
            "Gauntlet test failed to parse: {} - {:?}",
            test.input,
            result.err()
        );

        // Verify roundtrip
        let variant = result.unwrap();
        let displayed = format!("{}", variant);

        // Some variants may normalize during display, so just verify it re-parses
        let reparsed = parse_hgvs(&displayed);
        assert!(
            reparsed.is_ok(),
            "Gauntlet roundtrip failed: {} -> {} - {:?}",
            test.input,
            displayed,
            reparsed.err()
        );
        passed += 1;
    }

    println!("Gauntlet tests: {} passed", passed);
}

#[test]
fn test_biocommons_coverage_summary() {
    let fixtures = load_biocommons_fixtures();

    let total = fixtures.tests.len();
    let valid_count = fixtures.tests.iter().filter(|t| t.valid).count();
    let invalid_count = total - valid_count;

    let by_source: std::collections::HashMap<&str, usize> =
        fixtures
            .tests
            .iter()
            .fold(std::collections::HashMap::new(), |mut acc, t| {
                *acc.entry(t.source.as_str()).or_insert(0) += 1;
                acc
            });

    println!("\nBiocommons test coverage:");
    println!("  Total tests: {}", total);
    println!("  Valid variants: {}", valid_count);
    println!("  Invalid variants: {}", invalid_count);
    println!("  By source:");
    for (source, count) in by_source {
        println!("    {}: {}", source, count);
    }
}
