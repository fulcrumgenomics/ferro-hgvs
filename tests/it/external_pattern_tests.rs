//! Tests for external HGVS pattern sources
//!
//! This module tests patterns collected from biocommons, mutalyzer, NCBI dbSNP,
//! and other external sources to verify coverage and identify unsupported patterns.

use ferro_hgvs::parse_hgvs;
use std::collections::HashMap;
use std::fs;

/// Extract HGVS patterns from biocommons grammar_test.tsv
fn extract_grammar_test_patterns() -> Vec<(String, bool)> {
    let content = fs::read_to_string("tests/fixtures/external/biocommons_grammar_test.tsv")
        .expect("Failed to read biocommons_grammar_test.tsv");

    let mut patterns = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }

        let func = parts[0];
        let test = parts[1];
        let valid = parts[2].to_lowercase() == "true";

        // Only extract variant tests (hgvs_variant, c_variant, g_variant, etc.)
        if func.ends_with("_variant") || func == "hgvs_variant" {
            // Handle pipe-separated test cases
            for pattern in test.split('|') {
                let pattern = pattern.trim();
                if !pattern.is_empty() {
                    patterns.push((pattern.to_string(), valid));
                }
            }
        }
    }

    patterns
}

/// Extract HGVS patterns from biocommons gauntlet file
fn extract_gauntlet_patterns() -> Vec<(String, bool)> {
    let content = fs::read_to_string("tests/fixtures/external/biocommons_gauntlet.txt")
        .expect("Failed to read biocommons_gauntlet.txt");

    let mut patterns = Vec::new();

    for line in content.lines() {
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // Skip section headers
        if line.contains("####") {
            continue;
        }

        // Check if it's marked as unsupported
        if line.starts_with("#!unsupported:") {
            let pattern = line.trim_start_matches("#!unsupported:").trim();
            patterns.push((pattern.to_string(), false));
        } else if line.starts_with("#-") {
            // Negative examples
            let pattern = line.trim_start_matches("#-").trim();
            patterns.push((pattern.to_string(), false));
        } else {
            // Regular pattern
            patterns.push((line.to_string(), true));
        }
    }

    patterns
}

/// Extract HGVS patterns from biocommons real.tsv
fn extract_real_tsv_patterns() -> Vec<(String, bool)> {
    let content = fs::read_to_string("tests/fixtures/external/biocommons_real.tsv")
        .expect("Failed to read biocommons_real.tsv");

    let mut patterns = Vec::new();

    for line in content.lines().skip(1) {
        // Skip header
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 4 {
            // Extract HGVSg, HGVSc, HGVSp
            for i in 1..=3 {
                if let Some(pattern) = parts.get(i) {
                    let pattern = pattern.trim();
                    if !pattern.is_empty() {
                        patterns.push((pattern.to_string(), true));
                    }
                }
            }
        }
    }

    patterns
}

/// Extract HGVS patterns from mutalyzer test file
fn extract_mutalyzer_patterns() -> Vec<(String, bool)> {
    let content = fs::read_to_string("tests/fixtures/external/mutalyzer_patterns.txt")
        .expect("Failed to read mutalyzer_patterns.txt");

    let mut patterns = Vec::new();

    for line in content.lines() {
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // All mutalyzer patterns are valid (invalid ones in comments)
        patterns.push((line.to_string(), true));
    }

    patterns
}

/// Extract HGVS patterns from NCBI dbSNP file
fn extract_ncbi_patterns() -> Vec<(String, bool)> {
    let content = fs::read_to_string("tests/fixtures/external/ncbi_dbsnp.txt")
        .expect("Failed to read ncbi_dbsnp.txt");

    let mut patterns = Vec::new();

    for line in content.lines() {
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        patterns.push((line.to_string(), true));
    }

    patterns
}

/// Categorize parsing results
#[derive(Debug, Default)]
struct CategoryStats {
    supported_valid: Vec<String>,   // Expected valid, parsed successfully
    supported_invalid: Vec<String>, // Expected invalid, failed to parse (correct)
    unsupported: Vec<String>,       // Expected valid, failed to parse
    false_positive: Vec<String>,    // Expected invalid, parsed successfully (bug)
}

impl CategoryStats {
    fn summary(&self) -> String {
        format!(
            "Supported valid: {}, Correctly rejected: {}, Unsupported: {}, False positives: {}",
            self.supported_valid.len(),
            self.supported_invalid.len(),
            self.unsupported.len(),
            self.false_positive.len()
        )
    }
}

fn categorize_patterns(patterns: Vec<(String, bool)>) -> CategoryStats {
    let mut stats = CategoryStats::default();

    for (pattern, expected_valid) in patterns {
        let result = parse_hgvs(&pattern);

        match (result.is_ok(), expected_valid) {
            (true, true) => stats.supported_valid.push(pattern),
            (false, false) => stats.supported_invalid.push(pattern),
            (false, true) => stats.unsupported.push(pattern),
            (true, false) => stats.false_positive.push(pattern),
        }
    }

    stats
}

#[test]
fn test_biocommons_grammar_patterns() {
    let patterns = extract_grammar_test_patterns();
    let stats = categorize_patterns(patterns);

    println!("\n=== Biocommons Grammar Test Patterns ===");
    println!("{}", stats.summary());

    if !stats.unsupported.is_empty() {
        println!("\nUnsupported patterns (first 20):");
        for p in stats.unsupported.iter().take(20) {
            println!("  - {}", p);
        }
    }

    if !stats.false_positive.is_empty() {
        println!("\nFalse positives (should have rejected):");
        for p in &stats.false_positive {
            println!("  - {}", p);
        }
    }

    // Don't fail the test, just report
    println!(
        "Total patterns tested: {}",
        stats.supported_valid.len()
            + stats.supported_invalid.len()
            + stats.unsupported.len()
            + stats.false_positive.len()
    );
}

#[test]
fn test_biocommons_gauntlet_patterns() {
    let patterns = extract_gauntlet_patterns();
    let stats = categorize_patterns(patterns);

    println!("\n=== Biocommons Gauntlet Patterns ===");
    println!("{}", stats.summary());

    if !stats.unsupported.is_empty() {
        println!("\nUnsupported patterns:");
        for p in &stats.unsupported {
            println!("  - {}", p);
        }
    }
}

#[test]
fn test_biocommons_real_patterns() {
    let patterns = extract_real_tsv_patterns();
    let stats = categorize_patterns(patterns);

    println!("\n=== Biocommons Real TSV Patterns ===");
    println!("{}", stats.summary());

    if !stats.unsupported.is_empty() {
        println!("\nUnsupported patterns (first 20):");
        for p in stats.unsupported.iter().take(20) {
            println!("  - {}", p);
        }
    }

    // These are real ClinVar variants, so we should support most of them
    let total = stats.supported_valid.len() + stats.unsupported.len();
    let support_rate = stats.supported_valid.len() as f64 / total as f64 * 100.0;
    println!("Support rate: {:.1}%", support_rate);
}

#[test]
fn test_mutalyzer_patterns() {
    let patterns = extract_mutalyzer_patterns();
    let stats = categorize_patterns(patterns);

    println!("\n=== Mutalyzer Patterns ===");
    println!("{}", stats.summary());

    if !stats.unsupported.is_empty() {
        println!("\nUnsupported patterns:");
        for p in &stats.unsupported {
            println!("  - {}", p);
        }
    }
}

#[test]
fn test_ncbi_dbsnp_patterns() {
    let patterns = extract_ncbi_patterns();
    let stats = categorize_patterns(patterns);

    println!("\n=== NCBI dbSNP Patterns ===");
    println!("{}", stats.summary());

    // All NCBI dbSNP patterns should be supported
    assert!(
        stats.unsupported.is_empty(),
        "NCBI dbSNP patterns should all be supported: {:?}",
        stats.unsupported
    );
}

#[test]
fn test_comprehensive_pattern_report() {
    println!("\n============================================================");
    println!("       COMPREHENSIVE HGVS PATTERN SUPPORT REPORT");
    println!("============================================================\n");

    let mut all_unsupported: HashMap<String, Vec<String>> = HashMap::new();
    let mut total_supported = 0;
    let mut total_unsupported = 0;

    // Grammar test
    let patterns = extract_grammar_test_patterns();
    let stats = categorize_patterns(patterns);
    println!("Biocommons Grammar:");
    println!("  {}", stats.summary());
    total_supported += stats.supported_valid.len();
    total_unsupported += stats.unsupported.len();
    all_unsupported.insert("grammar".to_string(), stats.unsupported);

    // Gauntlet
    let patterns = extract_gauntlet_patterns();
    let stats = categorize_patterns(patterns);
    println!("Biocommons Gauntlet:");
    println!("  {}", stats.summary());
    total_supported += stats.supported_valid.len();
    total_unsupported += stats.unsupported.len();
    all_unsupported.insert("gauntlet".to_string(), stats.unsupported);

    // Real TSV
    let patterns = extract_real_tsv_patterns();
    let stats = categorize_patterns(patterns);
    println!("Biocommons Real:");
    println!("  {}", stats.summary());
    total_supported += stats.supported_valid.len();
    total_unsupported += stats.unsupported.len();
    all_unsupported.insert("real".to_string(), stats.unsupported);

    // Mutalyzer
    let patterns = extract_mutalyzer_patterns();
    let stats = categorize_patterns(patterns);
    println!("Mutalyzer:");
    println!("  {}", stats.summary());
    total_supported += stats.supported_valid.len();
    total_unsupported += stats.unsupported.len();
    all_unsupported.insert("mutalyzer".to_string(), stats.unsupported);

    // NCBI
    let patterns = extract_ncbi_patterns();
    let stats = categorize_patterns(patterns);
    println!("NCBI dbSNP:");
    println!("  {}", stats.summary());
    total_supported += stats.supported_valid.len();
    total_unsupported += stats.unsupported.len();
    all_unsupported.insert("ncbi".to_string(), stats.unsupported);

    println!("\n------------------------------------------------------------");
    println!("TOTALS:");
    println!("  Supported patterns: {}", total_supported);
    println!("  Unsupported patterns: {}", total_unsupported);
    let overall_rate =
        total_supported as f64 / (total_supported + total_unsupported) as f64 * 100.0;
    println!("  Overall support rate: {:.1}%", overall_rate);

    println!("\n------------------------------------------------------------");
    println!("ALL UNSUPPORTED PATTERNS:");

    for (source, patterns) in &all_unsupported {
        if !patterns.is_empty() {
            println!("\n[{}]:", source);
            for p in patterns {
                println!("  {}", p);
            }
        }
    }

    println!("\n============================================================");
}

#[test]
fn test_specific_patterns_to_verify() {
    // Test patterns that were marked for verification in the research plan
    let patterns_to_check = vec![
        // Selenocysteine
        ("NP_000001.1:p.Sec25Cys", "Selenocysteine substitution"),
        // Stop codon asterisk vs Ter
        ("NP_000001.1:p.Arg342*", "Stop codon with asterisk"),
        ("NP_000001.1:p.Arg342Ter", "Stop codon with Ter"),
        // Gene symbol in accession
        ("NM_004006.1(DMD):c.22+1A>T", "Gene symbol in accession"),
        // RNA lowercase
        ("NM_000001.1:r.76a>c", "RNA lowercase nucleotides"),
        // Repeat genotypes
        ("NM_000001.1:c.100CAG[14]", "Repeat genotype single"),
    ];

    println!("\n=== Specific Pattern Verification ===\n");

    for (pattern, description) in patterns_to_check {
        let result = parse_hgvs(pattern);
        let status = if result.is_ok() {
            "✓ SUPPORTED"
        } else {
            "✗ UNSUPPORTED"
        };
        println!("{}: {} - {}", status, description, pattern);
        if let Err(e) = result {
            println!("   Error: {}", e);
        }
    }
}
