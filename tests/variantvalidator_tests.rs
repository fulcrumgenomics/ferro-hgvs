//! VariantValidator Oracle Tests
//!
//! These tests validate parsing against known clinically-relevant variants
//! that have been validated by VariantValidator (https://variantvalidator.org/).
//!
//! The test cases represent common pathogenic variants used in clinical testing.

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::fs;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct TestFixture {
    description: String,
    source: String,
    note: String,
    test_cases: Vec<TestCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct TestCase {
    input: String,
    #[serde(rename = "type")]
    variant_type: String,
    valid: bool,
    description: String,
}

fn load_variantvalidator_fixtures() -> TestFixture {
    let content = fs::read_to_string("tests/fixtures/validation/variantvalidator.json")
        .expect("Failed to read variantvalidator.json");
    serde_json::from_str(&content).expect("Failed to parse variantvalidator.json")
}

#[test]
fn test_variantvalidator_clinical_variants() {
    let fixtures = load_variantvalidator_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in &fixtures.test_cases {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(variant) => {
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

                if actual_type == test.variant_type {
                    passed += 1;
                } else {
                    failed += 1;
                    failures.push(format!(
                        "{}: expected type {}, got {} ({})",
                        test.input, test.variant_type, actual_type, test.description
                    ));
                }
            }
            Err(e) => {
                if test.valid {
                    failed += 1;
                    failures.push(format!(
                        "{}: parse error: {} ({})",
                        test.input, e, test.description
                    ));
                } else {
                    passed += 1;
                }
            }
        }
    }

    if !failures.is_empty() {
        eprintln!("\nVariantValidator clinical variant failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    println!(
        "\nVariantValidator clinical variants: {} passed, {} failed",
        passed, failed
    );

    assert!(
        failed == 0,
        "VariantValidator: {} passed, {} failed",
        passed,
        failed
    );
}

#[test]
fn test_variantvalidator_roundtrip() {
    let fixtures = load_variantvalidator_fixtures();

    for test in fixtures.test_cases.iter().filter(|t| t.valid) {
        let parsed = parse_hgvs(&test.input);
        assert!(
            parsed.is_ok(),
            "Failed to parse {}: {:?} ({})",
            test.input,
            parsed.err(),
            test.description
        );

        let variant = parsed.unwrap();
        let displayed = format!("{}", variant);

        // Verify it can be reparsed
        let reparsed = parse_hgvs(&displayed);
        assert!(
            reparsed.is_ok(),
            "Roundtrip failed for {}: {} -> {} ({:?})",
            test.input,
            test.input,
            displayed,
            reparsed.err()
        );
    }
}

/// Test specific clinically important genes
mod clinical_genes {
    use super::*;

    #[test]
    fn test_brca_variants() {
        // BRCA1 and BRCA2 variants commonly tested in oncology
        let variants = vec![
            "NM_007294.4:c.5266dup",  // BRCA1
            "NM_000059.4:c.68_69del", // BRCA2
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse BRCA variant: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Cds(_)));
        }
    }

    #[test]
    fn test_tp53_variants() {
        // TP53 variants - most commonly mutated gene in cancer
        let variants = vec!["NM_000546.6:c.215C>G", "NM_000546.6:c.743G>A"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse TP53 variant: {}", v);
        }
    }

    #[test]
    fn test_braf_variants() {
        // BRAF V600E - common in melanoma
        let variants = vec!["NP_000079.2:p.Val600Glu"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse BRAF variant: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Protein(_)));
        }
    }

    #[test]
    fn test_cftr_variants() {
        // CFTR F508del - most common CF mutation
        let variants = vec!["NM_000492.4:c.1521_1523del"];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse CFTR variant: {}", v);
        }
    }

    #[test]
    fn test_mitochondrial_variants() {
        // Common mitochondrial disease mutations
        let variants = vec![
            "NC_012920.1:m.3243A>G", // MELAS
            "NC_012920.1:m.8344A>G", // MERRF
        ];

        for v in variants {
            let result = parse_hgvs(v);
            assert!(result.is_ok(), "Failed to parse MT variant: {}", v);
            assert!(matches!(result.unwrap(), HgvsVariant::Mt(_)));
        }
    }
}

#[test]
fn test_variantvalidator_coverage_summary() {
    let fixtures = load_variantvalidator_fixtures();

    let total = fixtures.test_cases.len();
    let by_type: std::collections::HashMap<&str, usize> =
        fixtures
            .test_cases
            .iter()
            .fold(std::collections::HashMap::new(), |mut acc, t| {
                *acc.entry(t.variant_type.as_str()).or_insert(0) += 1;
                acc
            });

    println!("\nVariantValidator test coverage:");
    println!("  Total test cases: {}", total);
    println!("  By variant type:");
    for (vtype, count) in by_type {
        println!("    {}: {}", vtype, count);
    }
}
