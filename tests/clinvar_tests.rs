//! ClinVar Parsing Validation Tests
//!
//! These tests validate the parser against real-world HGVS expressions
//! found in ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/).
//!
//! ClinVar is the primary source of clinical variant interpretations and
//! contains a wide variety of HGVS notation patterns, including edge cases
//! and non-standard formats commonly used in clinical practice.

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinVarFixture {
    description: String,
    source: String,
    version: String,
    note: String,
    categories: HashMap<String, String>,
    test_cases: Vec<TestCase>,
    edge_cases: Vec<TestCase>,
    quirks: Quirks,
    statistics: Statistics,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct TestCase {
    input: String,
    #[serde(rename = "type")]
    variant_type: String,
    valid: bool,
    category: String,
    #[serde(default)]
    gene: String,
    description: String,
    #[serde(default)]
    clinvar_id: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Quirks {
    description: String,
    items: Vec<QuirkItem>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct QuirkItem {
    pattern: String,
    examples: Vec<String>,
    handling: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Statistics {
    total_test_cases: usize,
    by_type: HashMap<String, usize>,
    by_category: HashMap<String, usize>,
    genes_covered: usize,
}

fn load_clinvar_fixtures() -> ClinVarFixture {
    let content = fs::read_to_string("tests/fixtures/validation/clinvar.json")
        .expect("Failed to read clinvar.json");
    serde_json::from_str(&content).expect("Failed to parse clinvar.json")
}

/// Main test that validates all ClinVar test cases
#[test]
fn test_clinvar_parsing_validation() {
    let fixtures = load_clinvar_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in &fixtures.test_cases {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(variant) => {
                let actual_type = get_variant_type(&variant);

                if actual_type == test.variant_type {
                    passed += 1;
                } else {
                    failed += 1;
                    failures.push(format!(
                        "{}: expected type {}, got {} [{}] ({})",
                        test.input, test.variant_type, actual_type, test.gene, test.description
                    ));
                }
            }
            Err(e) => {
                if test.valid {
                    failed += 1;
                    failures.push(format!(
                        "{}: parse error: {} [{}] ({})",
                        test.input, e, test.gene, test.description
                    ));
                } else {
                    passed += 1;
                }
            }
        }
    }

    if !failures.is_empty() {
        eprintln!("\nClinVar parsing validation failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    println!(
        "\nClinVar test cases: {} passed, {} failed out of {}",
        passed,
        failed,
        fixtures.test_cases.len()
    );

    assert!(
        failed == 0,
        "ClinVar validation: {} passed, {} failed",
        passed,
        failed
    );
}

/// Test edge cases specifically
#[test]
fn test_clinvar_edge_cases() {
    let fixtures = load_clinvar_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in &fixtures.edge_cases {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(variant) => {
                let actual_type = get_variant_type(&variant);

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
        eprintln!("\nClinVar edge case failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    println!(
        "\nClinVar edge cases: {} passed, {} failed out of {}",
        passed,
        failed,
        fixtures.edge_cases.len()
    );

    assert!(
        failed == 0,
        "ClinVar edge cases: {} passed, {} failed",
        passed,
        failed
    );
}

/// Roundtrip test: parse -> display -> parse again
#[test]
fn test_clinvar_roundtrip() {
    let fixtures = load_clinvar_fixtures();
    let mut roundtrip_failures = Vec::new();

    for test in fixtures.test_cases.iter().filter(|t| t.valid) {
        let parsed = parse_hgvs(&test.input);

        if let Ok(variant) = parsed {
            let displayed = format!("{}", variant);

            // Try to reparse the displayed string
            let reparsed = parse_hgvs(&displayed);
            if let Err(e) = reparsed {
                roundtrip_failures.push(format!(
                    "{} -> {} (reparse failed: {})",
                    test.input, displayed, e
                ));
            }
        }
    }

    if !roundtrip_failures.is_empty() {
        eprintln!("\nClinVar roundtrip failures:");
        for f in &roundtrip_failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        roundtrip_failures.is_empty(),
        "ClinVar roundtrip: {} failures",
        roundtrip_failures.len()
    );
}

/// Test variants by category
mod categories {
    use super::*;

    #[test]
    fn test_substitutions() {
        test_category("substitutions");
    }

    #[test]
    fn test_deletions() {
        test_category("deletions");
    }

    #[test]
    fn test_duplications() {
        test_category("duplications");
    }

    #[test]
    fn test_insertions() {
        test_category("insertions");
    }

    #[test]
    fn test_delins() {
        test_category("delins");
    }

    #[test]
    fn test_intronic() {
        test_category("intronic");
    }

    #[test]
    fn test_splice() {
        test_category("splice");
    }

    #[test]
    fn test_utr() {
        test_category("utr");
    }

    #[test]
    fn test_protein() {
        test_category("protein");
    }

    #[test]
    fn test_mitochondrial() {
        test_category("mitochondrial");
    }

    #[test]
    fn test_uncertain() {
        test_category("uncertain");
    }

    #[test]
    fn test_nonsense() {
        test_category("nonsense");
    }

    fn test_category(category: &str) {
        let fixtures = load_clinvar_fixtures();
        let category_tests: Vec<_> = fixtures
            .test_cases
            .iter()
            .filter(|t| t.category == category)
            .collect();

        let mut passed = 0;
        let mut failed = 0;

        for test in &category_tests {
            let result = parse_hgvs(&test.input);
            if result.is_ok() == test.valid {
                passed += 1;
            } else {
                failed += 1;
                eprintln!(
                    "  {} category failure: {} - {}",
                    category, test.input, test.description
                );
            }
        }

        println!(
            "\nClinVar {} category: {} passed, {} failed",
            category, passed, failed
        );
        assert!(
            failed == 0,
            "ClinVar {} category: {} failures",
            category,
            failed
        );
    }
}

/// Test specific clinical genes commonly found in ClinVar
mod clinical_genes {
    use super::*;

    #[test]
    fn test_brca_variants() {
        test_gene_variants(&["BRCA1", "BRCA2"]);
    }

    #[test]
    fn test_tp53_variants() {
        test_gene_variants(&["TP53"]);
    }

    #[test]
    fn test_lynch_syndrome_genes() {
        test_gene_variants(&["MLH1", "MSH2", "MSH6", "PMS2"]);
    }

    #[test]
    fn test_cftr_variants() {
        test_gene_variants(&["CFTR"]);
    }

    #[test]
    fn test_oncogenes() {
        test_gene_variants(&["BRAF", "KRAS", "EGFR", "PIK3CA"]);
    }

    #[test]
    fn test_myeloid_genes() {
        test_gene_variants(&["JAK2", "DNMT3A", "NPM1", "FLT3", "IDH1", "IDH2"]);
    }

    fn test_gene_variants(genes: &[&str]) {
        let fixtures = load_clinvar_fixtures();
        let gene_tests: Vec<_> = fixtures
            .test_cases
            .iter()
            .filter(|t| genes.iter().any(|g| t.gene == *g))
            .collect();

        if gene_tests.is_empty() {
            println!("No test cases found for genes: {:?}", genes);
            return;
        }

        let mut passed = 0;
        let mut failed = 0;

        for test in &gene_tests {
            let result = parse_hgvs(&test.input);
            if result.is_ok() == test.valid {
                passed += 1;
            } else {
                failed += 1;
                eprintln!(
                    "  {} gene failure: {} - {}",
                    test.gene, test.input, test.description
                );
            }
        }

        println!(
            "\nClinVar genes {:?}: {} passed, {} failed",
            genes, passed, failed
        );
        assert!(
            failed == 0,
            "ClinVar genes {:?}: {} failures",
            genes,
            failed
        );
    }
}

/// Test mitochondrial variants
#[test]
fn test_clinvar_mitochondrial() {
    let fixtures = load_clinvar_fixtures();
    let mt_variants: Vec<_> = fixtures
        .test_cases
        .iter()
        .filter(|t| t.variant_type == "Mt")
        .collect();

    for test in mt_variants {
        let result = parse_hgvs(&test.input);
        assert!(
            result.is_ok(),
            "Failed to parse MT variant {}: {:?} ({})",
            test.input,
            result.err(),
            test.description
        );

        let variant = result.unwrap();
        assert!(
            matches!(variant, HgvsVariant::Mt(_)),
            "Expected Mt variant type for {}",
            test.input
        );
    }
}

/// Test protein stop codon notations (Ter vs *)
#[test]
fn test_clinvar_stop_codon_notations() {
    // Both Ter and * should work for stop codons
    let ter_notation = parse_hgvs("NP_000537.3:p.Arg342Ter");
    let asterisk_notation = parse_hgvs("NP_000537.3:p.Arg342*");

    assert!(ter_notation.is_ok(), "Failed to parse Ter notation");
    assert!(asterisk_notation.is_ok(), "Failed to parse * notation");

    // Both should be protein variants
    assert!(matches!(ter_notation.unwrap(), HgvsVariant::Protein(_)));
    assert!(matches!(
        asterisk_notation.unwrap(),
        HgvsVariant::Protein(_)
    ));
}

/// Test deep intronic variants with large offsets
#[test]
fn test_clinvar_deep_intronic() {
    let deep_intronic = vec![
        ("NM_000492.4:c.3718-2477C>T", "CFTR deep intronic"),
        ("NM_206933.4:c.7595-2144A>G", "USH2A deep intronic"),
        ("NM_000169.3:c.639+919G>A", "GLA deep intronic"),
        ("NM_000152.5:c.-32-13T>G", "GAA deep intronic"),
    ];

    for (input, desc) in deep_intronic {
        let result = parse_hgvs(input);
        assert!(
            result.is_ok(),
            "Failed to parse deep intronic variant {}: {:?} ({})",
            input,
            result.err(),
            desc
        );
    }
}

/// Test uncertainty notation - currently supported patterns
#[test]
fn test_clinvar_uncertainty_notation() {
    // Currently supported uncertainty patterns
    let supported_uncertain = vec![
        "NP_000537.3:p.Met1?",       // Start codon uncertainty - supported
        "NM_000546.6:c.1_?del",      // Unknown end position
        "NM_000546.6:c.?_1182del",   // Unknown start position
        "NP_000537.3:p.?",           // Whole-protein unknown
        "NP_000537.3:p.(Arg248Gln)", // Predicted protein change
    ];

    for input in supported_uncertain {
        let result = parse_hgvs(input);
        assert!(
            result.is_ok(),
            "Failed to parse supported uncertain variant {}: {:?}",
            input,
            result.err()
        );
    }

    // Verify roundtrip for predicted protein change
    let variant = parse_hgvs("NP_000537.3:p.(Arg248Gln)").unwrap();
    assert_eq!(format!("{}", variant), "NP_000537.3:p.(Arg248Gln)");
}

/// Coverage summary for ClinVar tests
#[test]
fn test_clinvar_coverage_summary() {
    let fixtures = load_clinvar_fixtures();

    let total = fixtures.test_cases.len() + fixtures.edge_cases.len();
    let mut by_type: HashMap<&str, usize> = HashMap::new();
    let mut by_category: HashMap<&str, usize> = HashMap::new();
    let mut genes: std::collections::HashSet<&str> = std::collections::HashSet::new();

    for test in &fixtures.test_cases {
        *by_type.entry(test.variant_type.as_str()).or_insert(0) += 1;
        *by_category.entry(test.category.as_str()).or_insert(0) += 1;
        if !test.gene.is_empty() {
            genes.insert(&test.gene);
        }
    }

    for test in &fixtures.edge_cases {
        *by_type.entry(test.variant_type.as_str()).or_insert(0) += 1;
        *by_category.entry(test.category.as_str()).or_insert(0) += 1;
    }

    println!("\n=== ClinVar Test Coverage ===");
    println!("Source: {}", fixtures.source);
    println!("Version: {}", fixtures.version);
    println!("\nTotal test cases: {}", total);
    println!("  - Main test cases: {}", fixtures.test_cases.len());
    println!("  - Edge cases: {}", fixtures.edge_cases.len());
    println!("\nBy variant type:");
    for (vtype, count) in &by_type {
        println!("  {}: {}", vtype, count);
    }
    println!("\nBy category:");
    for (cat, count) in &by_category {
        println!("  {}: {}", cat, count);
    }
    println!("\nGenes covered: {}", genes.len());
}

/// Helper to get variant type string
fn get_variant_type(variant: &HgvsVariant) -> &'static str {
    match variant {
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
    }
}
