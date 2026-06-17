//! Edge Case Notation Tests
//!
//! These tests validate parsing against unusual but valid HGVS patterns
//! including legacy notations, complex alleles, uncertain positions,
//! comprehensive protein notation coverage, and reference sequence formats.

use ferro_hgvs::parse_hgvs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

// Unusual Notation Fixtures

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct UnusualNotationFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    categories: HashMap<String, String>,
    test_cases: Vec<UnusualNotationCase>,
    statistics: UnusualStats,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct UnusualNotationCase {
    category: String,
    input: String,
    description: String,
    valid: bool,
    supported: bool,
    #[serde(default)]
    notes: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct UnusualStats {
    total_cases: usize,
    by_category: HashMap<String, usize>,
    supported_count: usize,
    unsupported_count: usize,
}

// Protein Notation Fixtures

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ProteinNotationFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    categories: HashMap<String, String>,
    test_cases: Vec<ProteinNotationCase>,
    amino_acid_codes: AminoAcidCodes,
    statistics: ProteinStats,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ProteinNotationCase {
    category: String,
    input: String,
    description: String,
    valid: bool,
    supported: bool,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct AminoAcidCodes {
    three_letter: Vec<String>,
    one_letter: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ProteinStats {
    total_cases: usize,
    by_category: HashMap<String, usize>,
    supported_count: usize,
    unsupported_count: usize,
}

fn load_unusual_notation_fixtures() -> UnusualNotationFixture {
    let content = fs::read_to_string("tests/fixtures/edge_cases/unusual_notation.json")
        .expect("Failed to read unusual_notation.json");
    serde_json::from_str(&content).expect("Failed to parse unusual_notation.json")
}

fn load_protein_notation_fixtures() -> ProteinNotationFixture {
    let content = fs::read_to_string("tests/fixtures/edge_cases/protein_notation.json")
        .expect("Failed to read protein_notation.json");
    serde_json::from_str(&content).expect("Failed to parse protein_notation.json")
}

/// Test unusual notation patterns (supported variants)
#[test]
fn test_unusual_notation_supported() {
    let fixtures = load_unusual_notation_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures.test_cases.iter().filter(|t| t.supported) {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{} ({}): {}", test.input, test.category, e));
            }
        }
    }

    eprintln!("\nUnusual notation (supported) test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    let total = passed + failed;
    let success_rate = if total > 0 {
        (passed as f64 / total as f64) * 100.0
    } else {
        100.0
    };

    eprintln!("  Success rate: {:.1}%", success_rate);

    // Require at least 80% success rate for supported patterns
    assert!(
        success_rate >= 80.0,
        "Unusual notation success rate {:.1}% below 80% threshold",
        success_rate
    );
}

/// Test unusual notation patterns by category
#[test]
fn test_unusual_notation_by_category() {
    let fixtures = load_unusual_notation_fixtures();

    let mut category_counts: HashMap<&str, (usize, usize)> = HashMap::new();

    for test in &fixtures.test_cases {
        if !test.supported {
            continue;
        }

        let result = parse_hgvs(&test.input);
        let entry = category_counts.entry(&test.category).or_insert((0, 0));
        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nUnusual notation by category (supported only):");
    for (category, (passed, failed)) in &category_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", category, passed, total, rate);
    }
}

/// Test protein notation patterns (supported variants)
#[test]
fn test_protein_notation_supported() {
    let fixtures = load_protein_notation_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures.test_cases.iter().filter(|t| t.supported) {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{} ({}): {}", test.input, test.category, e));
            }
        }
    }

    eprintln!("\nProtein notation (supported) test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    let total = passed + failed;
    let success_rate = if total > 0 {
        (passed as f64 / total as f64) * 100.0
    } else {
        100.0
    };

    eprintln!("  Success rate: {:.1}%", success_rate);

    // Require at least 90% success rate for supported protein patterns
    assert!(
        success_rate >= 90.0,
        "Protein notation success rate {:.1}% below 90% threshold",
        success_rate
    );
}

/// Test protein notation patterns by category
#[test]
fn test_protein_notation_by_category() {
    let fixtures = load_protein_notation_fixtures();

    let mut category_counts: HashMap<&str, (usize, usize)> = HashMap::new();

    for test in &fixtures.test_cases {
        if !test.supported {
            continue;
        }

        let result = parse_hgvs(&test.input);
        let entry = category_counts.entry(&test.category).or_insert((0, 0));
        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nProtein notation by category (supported only):");
    for (category, (passed, failed)) in &category_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", category, passed, total, rate);
    }
}

/// Test all protein notation patterns including unsupported
#[test]
fn test_protein_notation_all() {
    let fixtures = load_protein_notation_fixtures();
    let mut supported_pass = 0;
    let mut supported_fail = 0;
    let mut unsupported_pass = 0;
    let mut unsupported_fail = 0;

    for test in &fixtures.test_cases {
        let result = parse_hgvs(&test.input);

        if test.supported {
            if result.is_ok() {
                supported_pass += 1;
            } else {
                supported_fail += 1;
            }
        } else if result.is_ok() {
            unsupported_pass += 1;
        } else {
            unsupported_fail += 1;
        }
    }

    eprintln!("\nProtein notation overall results:");
    eprintln!(
        "  Supported: {}/{} passed",
        supported_pass,
        supported_pass + supported_fail
    );
    eprintln!(
        "  Unsupported: {}/{} passed (bonus!)",
        unsupported_pass,
        unsupported_pass + unsupported_fail
    );
    eprintln!(
        "  Total patterns: {}",
        supported_pass + supported_fail + unsupported_pass + unsupported_fail
    );
}

/// Test one-letter amino acid codes specifically
#[test]
fn test_one_letter_amino_acid_codes() {
    let fixtures = load_protein_notation_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures
        .test_cases
        .iter()
        .filter(|t| t.category == "one_letter" && t.supported)
    {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", test.input, e));
            }
        }
    }

    eprintln!("\nOne-letter amino acid code test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // All one-letter codes should parse
    assert_eq!(failed, 0, "All one-letter code variants should parse");
}

/// Test frameshift notation specifically
#[test]
fn test_frameshift_notation() {
    let fixtures = load_protein_notation_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures
        .test_cases
        .iter()
        .filter(|t| t.category == "frameshift" && t.supported)
    {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", test.input, e));
            }
        }
    }

    eprintln!("\nFrameshift notation test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // All frameshift notations should parse
    assert_eq!(failed, 0, "All frameshift variants should parse");
}

/// Test extension notation specifically
#[test]
fn test_extension_notation() {
    let fixtures = load_protein_notation_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures
        .test_cases
        .iter()
        .filter(|t| t.category == "extension" && t.supported)
    {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", test.input, e));
            }
        }
    }

    eprintln!("\nExtension notation test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }
}

// Reference Sequence Fixtures

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ReferenceSeqFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    categories: HashMap<String, String>,
    test_cases: Vec<ReferenceSeqCase>,
    statistics: ReferenceSeqStats,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ReferenceSeqCase {
    category: String,
    input: String,
    description: String,
    valid: bool,
    supported: bool,
    #[serde(default)]
    notes: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ReferenceSeqStats {
    total_cases: usize,
    by_category: HashMap<String, usize>,
    supported_count: usize,
    unsupported_count: usize,
}

fn load_reference_seq_fixtures() -> ReferenceSeqFixture {
    let content = fs::read_to_string("tests/fixtures/edge_cases/reference_sequences.json")
        .expect("Failed to read reference_sequences.json");
    serde_json::from_str(&content).expect("Failed to parse reference_sequences.json")
}

/// Test reference sequence patterns (supported variants)
#[test]
fn test_reference_sequences_supported() {
    let fixtures = load_reference_seq_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures.test_cases.iter().filter(|t| t.supported) {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{} ({}): {}", test.input, test.category, e));
            }
        }
    }

    eprintln!("\nReference sequence (supported) test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    let total = passed + failed;
    let success_rate = if total > 0 {
        (passed as f64 / total as f64) * 100.0
    } else {
        100.0
    };

    eprintln!("  Success rate: {:.1}%", success_rate);

    // Require at least 90% success rate for supported reference patterns
    assert!(
        success_rate >= 90.0,
        "Reference sequence success rate {:.1}% below 90% threshold",
        success_rate
    );
}

/// Test reference sequence patterns by category
#[test]
fn test_reference_sequences_by_category() {
    let fixtures = load_reference_seq_fixtures();

    let mut category_counts: HashMap<&str, (usize, usize)> = HashMap::new();

    for test in &fixtures.test_cases {
        if !test.supported {
            continue;
        }

        let result = parse_hgvs(&test.input);
        let entry = category_counts.entry(&test.category).or_insert((0, 0));
        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nReference sequences by category (supported only):");
    for (category, (passed, failed)) in &category_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", category, passed, total, rate);
    }
}

/// Test LRG references specifically
#[test]
fn test_lrg_references() {
    let fixtures = load_reference_seq_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures
        .test_cases
        .iter()
        .filter(|t| t.category == "lrg" && t.supported)
    {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", test.input, e));
            }
        }
    }

    eprintln!("\nLRG reference test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // All LRG references should parse
    assert_eq!(failed, 0, "All LRG references should parse");
}

/// Test mitochondrial references specifically
#[test]
fn test_mitochondrial_references() {
    let fixtures = load_reference_seq_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures
        .test_cases
        .iter()
        .filter(|t| t.category == "mitochondrial" && t.supported)
    {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", test.input, e));
            }
        }
    }

    eprintln!("\nMitochondrial reference test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // All mtDNA references should parse
    assert_eq!(failed, 0, "All mitochondrial references should parse");
}

/// Test nested references specifically
#[test]
fn test_nested_references() {
    let fixtures = load_reference_seq_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for test in fixtures
        .test_cases
        .iter()
        .filter(|t| t.category == "nested" && t.supported)
    {
        let result = parse_hgvs(&test.input);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", test.input, e));
            }
        }
    }

    eprintln!("\nNested reference test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // All nested references should parse
    assert_eq!(failed, 0, "All nested references should parse");
}
