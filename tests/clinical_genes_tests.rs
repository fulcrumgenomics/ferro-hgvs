//! Clinical Genes Tests
//!
//! These tests validate parsing against real clinical variants from ClinVar,
//! focusing on pathogenic/likely pathogenic variants from clinically important genes.

use ferro_hgvs::parse_hgvs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinicalGenesFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    genes_queried: Vec<String>,
    total_test_cases: usize,
    summary: ClinicalSummary,
    test_cases: Vec<ClinicalGeneCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinicalSummary {
    by_gene: HashMap<String, usize>,
    by_type: HashMap<String, usize>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinicalGeneCase {
    input: String,
    #[serde(rename = "type")]
    variant_type: String,
    valid: bool,
    category: String,
    gene: String,
    description: String,
    clinvar_id: String,
}

fn load_clinical_genes_fixtures() -> ClinicalGenesFixture {
    let content = fs::read_to_string("tests/fixtures/validation/clinical_genes.json")
        .expect("Failed to read clinical_genes.json");
    serde_json::from_str(&content).expect("Failed to parse clinical_genes.json")
}

/// Test all clinical gene variants (coding variants)
#[test]
fn test_clinical_genes_coding_variants() {
    let fixture = load_clinical_genes_fixtures();

    let coding_cases: Vec<_> = fixture
        .test_cases
        .iter()
        .filter(|tc| tc.variant_type == "coding")
        .collect();

    assert!(
        !coding_cases.is_empty(),
        "Should have coding variants to test"
    );

    let mut passed = 0;

    for case in &coding_cases {
        let result = parse_hgvs(&case.input);
        if result.is_ok() {
            passed += 1;
        }
    }

    let total = coding_cases.len();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    eprintln!(
        "\nClinical Genes - Coding Variants: {}/{} passed ({:.1}%)",
        passed, total, pass_rate
    );

    // Note: Many ClinVar variants use non-standard notation formats
    // This test tracks progress - low pass rate expected initially
    assert!(
        pass_rate >= 15.0,
        "Expected at least 15% pass rate for clinical coding variants, got {:.1}%",
        pass_rate
    );
}

/// Test all clinical gene variants (all types)
#[test]
fn test_clinical_genes_all_variants() {
    let fixture = load_clinical_genes_fixtures();

    let mut passed = 0;
    let mut failed = 0;
    let mut by_type: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        let entry = by_type.entry(case.variant_type.clone()).or_insert((0, 0));

        if result.is_ok() {
            passed += 1;
            entry.0 += 1;
        } else {
            failed += 1;
            entry.1 += 1;
        }
    }

    let total = fixture.test_cases.len();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    eprintln!("\nClinical Genes - All Variants Summary:");
    eprintln!("  Total: {}", total);
    eprintln!("  Passed: {} ({:.1}%)", passed, pass_rate);
    eprintln!("  Failed: {}", failed);

    eprintln!("\nBy variant type:");
    for (vtype, (p, f)) in by_type.iter() {
        let type_total = p + f;
        let type_rate = (*p as f64 / type_total as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", vtype, p, type_total, type_rate);
    }

    // Note: Many ClinVar variants use non-standard notation formats
    // This test tracks progress - low pass rate expected initially
    assert!(
        pass_rate >= 15.0,
        "Expected at least 15% pass rate for clinical variants, got {:.1}%",
        pass_rate
    );
}

/// Test clinical variants by gene
#[test]
fn test_clinical_genes_by_gene() {
    let fixture = load_clinical_genes_fixtures();

    let mut by_gene: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        let entry = by_gene.entry(case.gene.clone()).or_insert((0, 0));

        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nClinical Genes - By Gene:");
    let mut genes: Vec<_> = by_gene.iter().collect();
    genes.sort_by(|a, b| (b.1 .0 + b.1 .1).cmp(&(a.1 .0 + a.1 .1)));

    for (gene, (passed, failed)) in genes.iter().take(15) {
        let total = passed + failed;
        let rate = (*passed as f64 / total as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", gene, passed, total, rate);
    }
}

/// Informational test - show ClinVar IDs for failed variants
#[test]
fn test_clinical_genes_show_failures() {
    let fixture = load_clinical_genes_fixtures();

    let mut failures: Vec<(&str, &str, &str)> = Vec::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        if result.is_err() {
            failures.push((&case.input, &case.gene, &case.clinvar_id));
        }
    }

    if !failures.is_empty() {
        eprintln!("\nClinical Genes - Failed Variants (first 10):");
        for (input, gene, clinvar_id) in failures.iter().take(10) {
            eprintln!("  ClinVar {}: {} ({})", clinvar_id, input, gene);
        }
        eprintln!("  ... and {} more", failures.len().saturating_sub(10));
    }
}
