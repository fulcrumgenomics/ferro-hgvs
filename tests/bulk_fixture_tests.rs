//! Bulk Fixture Tests
//!
//! These tests validate parsing against large-scale fixtures from:
//! - ClinVar bulk variants
//! - gnomAD constrained gene variants
//! - Parse gap analysis results

use ferro_hgvs::parse_hgvs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

// ClinVar Bulk Fixture

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinvarBulkFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    total_test_cases: usize,
    summary: ClinvarBulkSummary,
    test_cases: Vec<ClinvarBulkCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinvarBulkSummary {
    by_coordinate_type: HashMap<String, usize>,
    by_query_source: HashMap<String, usize>,
    by_clinical_significance: HashMap<String, usize>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ClinvarBulkCase {
    input: String,
    #[serde(rename = "type")]
    variant_type: String,
    query_source: String,
    gene: String,
    clinical_significance: String,
    #[serde(default)]
    molecular_consequence: String,
    clinvar_id: String,
    valid: bool,
}

// gnomAD Constrained Fixture

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct GnomadConstrainedFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    #[serde(default)]
    notes: String,
    total_test_cases: usize,
    summary: GnomadSummary,
    test_cases: Vec<GnomadCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct GnomadSummary {
    by_coordinate_type: HashMap<String, usize>,
    by_gene: HashMap<String, usize>,
    by_consequence: HashMap<String, usize>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct GnomadCase {
    input: String,
    #[serde(rename = "type")]
    variant_type: String,
    gene: String,
    #[serde(default)]
    transcript: String,
    consequence: String,
    source: String,
    valid: bool,
}

// Parse Gaps Fixture

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParseGapsFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    summary: ParseGapsSummary,
    categories: HashMap<String, String>,
    test_cases: Vec<ParseGapCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParseGapsSummary {
    total_failures: usize,
    by_category: HashMap<String, usize>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ParseGapCase {
    input: String,
    category: String,
    source: String,
    #[serde(default)]
    description: String,
    #[serde(default)]
    error_sample: String,
    valid: bool,
    supported: bool,
}

fn load_clinvar_bulk_fixture() -> Option<ClinvarBulkFixture> {
    let path = "tests/fixtures/bulk/clinvar_bulk.json";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    let content = fs::read_to_string(path).expect("Failed to read clinvar_bulk.json");
    Some(serde_json::from_str(&content).expect("Failed to parse clinvar_bulk.json"))
}

fn load_gnomad_fixture() -> Option<GnomadConstrainedFixture> {
    let path = "tests/fixtures/bulk/gnomad_constrained.json";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    let content = fs::read_to_string(path).expect("Failed to read gnomad_constrained.json");
    Some(serde_json::from_str(&content).expect("Failed to parse gnomad_constrained.json"))
}

fn load_parse_gaps_fixture() -> Option<ParseGapsFixture> {
    let path = "tests/fixtures/validation/parse_gaps.json";
    if !std::path::Path::new(path).exists() {
        return None;
    }
    let content = fs::read_to_string(path).expect("Failed to read parse_gaps.json");
    Some(serde_json::from_str(&content).expect("Failed to parse parse_gaps.json"))
}

// ============================================================================
// ClinVar Bulk Tests
// ============================================================================

#[test]
fn test_clinvar_bulk_coding_variants() {
    let fixture = match load_clinvar_bulk_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping: clinvar_bulk.json not found");
            return;
        }
    };

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
    let mut by_source: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &coding_cases {
        let result = parse_hgvs(&case.input);
        let entry = by_source.entry(case.query_source.clone()).or_insert((0, 0));

        if result.is_ok() {
            passed += 1;
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    let total = coding_cases.len();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    eprintln!(
        "\nClinVar Bulk - Coding Variants: {}/{} passed ({:.1}%)",
        passed, total, pass_rate
    );
    eprintln!("\nBy query source:");
    for (source, (p, f)) in by_source.iter() {
        let t = p + f;
        let r = (*p as f64 / t as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", source, p, t, r);
    }

    // ClinVar bulk should have reasonable parse rate
    assert!(
        pass_rate >= 10.0,
        "Expected at least 10% pass rate for ClinVar bulk, got {:.1}%",
        pass_rate
    );
}

#[test]
fn test_clinvar_bulk_by_clinical_significance() {
    let fixture = match load_clinvar_bulk_fixture() {
        Some(f) => f,
        None => return,
    };

    let mut by_clin_sig: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        let sig = if case.clinical_significance.is_empty() {
            "unknown".to_string()
        } else {
            case.clinical_significance.to_lowercase()
        };
        let entry = by_clin_sig.entry(sig).or_insert((0, 0));

        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nClinVar Bulk - By Clinical Significance:");
    for (sig, (passed, failed)) in by_clin_sig.iter() {
        let total = passed + failed;
        let rate = (*passed as f64 / total as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", sig, passed, total, rate);
    }
}

// ============================================================================
// gnomAD Constrained Tests
// ============================================================================

#[test]
fn test_gnomad_constrained_all() {
    let fixture = match load_gnomad_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping: gnomad_constrained.json not found");
            return;
        }
    };

    let mut passed = 0;
    let mut by_type: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        let entry = by_type.entry(case.variant_type.clone()).or_insert((0, 0));

        if result.is_ok() {
            passed += 1;
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    let total = fixture.test_cases.len();
    let pass_rate = (passed as f64 / total as f64) * 100.0;

    eprintln!(
        "\ngnomAD Constrained - All Variants: {}/{} passed ({:.1}%)",
        passed, total, pass_rate
    );
    eprintln!("\nBy coordinate type:");
    for (t, (p, f)) in by_type.iter() {
        let tot = p + f;
        let r = (*p as f64 / tot as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", t, p, tot, r);
    }

    // gnomAD constrained uses standard RefSeq accessions - should have high pass rate
    assert!(
        pass_rate >= 80.0,
        "Expected at least 80% pass rate for gnomAD constrained, got {:.1}%",
        pass_rate
    );
}

#[test]
fn test_gnomad_constrained_by_gene() {
    let fixture = match load_gnomad_fixture() {
        Some(f) => f,
        None => return,
    };

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

    eprintln!("\ngnomAD Constrained - By Gene (top 10):");
    let mut genes: Vec<_> = by_gene.iter().collect();
    genes.sort_by(|a, b| (b.1 .0 + b.1 .1).cmp(&(a.1 .0 + a.1 .1)));

    for (gene, (passed, failed)) in genes.iter().take(10) {
        let total = passed + failed;
        let rate = (*passed as f64 / total as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", gene, passed, total, rate);
    }
}

#[test]
fn test_gnomad_constrained_by_consequence() {
    let fixture = match load_gnomad_fixture() {
        Some(f) => f,
        None => return,
    };

    let mut by_consequence: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        let result = parse_hgvs(&case.input);
        let entry = by_consequence
            .entry(case.consequence.clone())
            .or_insert((0, 0));

        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\ngnomAD Constrained - By Consequence:");
    for (consequence, (passed, failed)) in by_consequence.iter() {
        let total = passed + failed;
        let rate = (*passed as f64 / total as f64) * 100.0;
        eprintln!("  {}: {}/{} ({:.1}%)", consequence, passed, total, rate);
    }
}

// ============================================================================
// Parse Gaps Tests (Informational)
// ============================================================================

#[test]
fn test_parse_gaps_informational() {
    let fixture = match load_parse_gaps_fixture() {
        Some(f) => f,
        None => {
            eprintln!("Skipping: parse_gaps.json not found");
            return;
        }
    };

    let mut by_category: HashMap<String, (usize, usize)> = HashMap::new();

    for case in &fixture.test_cases {
        // Parse gaps are expected to fail, but let's see if any now pass
        let result = parse_hgvs(&case.input);
        let entry = by_category.entry(case.category.clone()).or_insert((0, 0));

        if result.is_ok() {
            entry.0 += 1; // Now passes
        } else {
            entry.1 += 1; // Still fails
        }
    }

    let total = fixture.test_cases.len();
    let now_passing: usize = by_category.values().map(|(p, _)| p).sum();
    let still_failing: usize = by_category.values().map(|(_, f)| f).sum();

    eprintln!("\nParse Gaps Analysis:");
    eprintln!("  Total gap patterns: {}", total);
    eprintln!("  Now passing: {} (implementation improved!)", now_passing);
    eprintln!("  Still failing: {}", still_failing);

    eprintln!("\nBy gap category:");
    for (category, (passed, failed)) in by_category.iter() {
        let t = passed + failed;
        eprintln!(
            "  {}: {} passing, {} still failing (of {})",
            category, passed, failed, t
        );
    }

    // This is informational - we expect most to fail
    // But if any are now passing, that's progress!
}

// ============================================================================
// Summary Test
// ============================================================================

#[test]
fn test_bulk_fixtures_summary() {
    eprintln!("\n========================================");
    eprintln!("Bulk Fixtures Summary");
    eprintln!("========================================\n");

    // ClinVar Bulk
    if let Some(fixture) = load_clinvar_bulk_fixture() {
        let mut passed = 0;
        for case in &fixture.test_cases {
            if parse_hgvs(&case.input).is_ok() {
                passed += 1;
            }
        }
        let total = fixture.test_cases.len();
        let rate = (passed as f64 / total as f64) * 100.0;
        eprintln!("ClinVar Bulk: {}/{} ({:.1}%)", passed, total, rate);
    }

    // gnomAD Constrained
    if let Some(fixture) = load_gnomad_fixture() {
        let mut passed = 0;
        for case in &fixture.test_cases {
            if parse_hgvs(&case.input).is_ok() {
                passed += 1;
            }
        }
        let total = fixture.test_cases.len();
        let rate = (passed as f64 / total as f64) * 100.0;
        eprintln!("gnomAD Constrained: {}/{} ({:.1}%)", passed, total, rate);
    }

    // Parse Gaps
    if let Some(fixture) = load_parse_gaps_fixture() {
        let mut passed = 0;
        for case in &fixture.test_cases {
            if parse_hgvs(&case.input).is_ok() {
                passed += 1;
            }
        }
        let total = fixture.test_cases.len();
        let rate = (passed as f64 / total as f64) * 100.0;
        eprintln!(
            "Parse Gaps (now passing): {}/{} ({:.1}%)",
            passed, total, rate
        );
    }

    eprintln!("\n========================================\n");
}
