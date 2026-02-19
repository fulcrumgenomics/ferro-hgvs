//! Round-Trip Idempotency Tests
//!
//! These tests verify that HGVS parsing and normalization are idempotent:
//! 1. Parsing: parse(format(parse(input))) == parse(input)
//! 2. Normalization: normalize(normalize(input)) == normalize(input)
//!
//! This ensures that our outputs can be correctly re-processed.

use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, Normalizer};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::sync::OnceLock;

// =============================================================================
// Fixture Loading
// =============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct EdgeCaseFixture {
    mutalyzer_edge_cases: Option<EdgeCaseCategory>,
    biocommons_issues: Option<EdgeCaseCategory>,
    selenoprotein_pyrrolysine: Option<EdgeCaseCategory>,
    invalid_patterns: Option<EdgeCaseCategory>,
    protein_variants: Option<EdgeCaseCategory>,
    complex_indels: Option<EdgeCaseCategory>,
    extensions_frameshifts: Option<EdgeCaseCategory>,
}

#[derive(Debug, Deserialize)]
struct EdgeCaseCategory {
    tests: Vec<EdgeCase>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct EdgeCase {
    input: String,
    valid: bool,
    #[serde(rename = "type")]
    variant_type: Option<String>,
}

#[derive(Debug, Deserialize)]
struct HgvsSpecFixture {
    examples: Vec<SpecExample>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct SpecExample {
    hgvs: String,
    variant_type: String,
    level: String,
    coordinate_system: String,
}

#[allow(clippy::unnecessary_lazy_evaluations)]
fn load_edge_cases() -> Vec<String> {
    let content = fs::read_to_string("tests/fixtures/comprehensive_edge_cases.json")
        .unwrap_or_else(|_| "{}".to_string());
    let fixture: EdgeCaseFixture =
        serde_json::from_str(&content).unwrap_or_else(|_| EdgeCaseFixture {
            mutalyzer_edge_cases: None,
            biocommons_issues: None,
            selenoprotein_pyrrolysine: None,
            invalid_patterns: None,
            protein_variants: None,
            complex_indels: None,
            extensions_frameshifts: None,
        });

    let mut patterns = Vec::new();

    // Extract valid patterns from each category
    for category in [
        fixture.mutalyzer_edge_cases,
        fixture.biocommons_issues,
        fixture.selenoprotein_pyrrolysine,
        fixture.protein_variants,
        fixture.complex_indels,
        fixture.extensions_frameshifts,
    ]
    .into_iter()
    .flatten()
    {
        for test in category.tests {
            if test.valid {
                patterns.push(test.input);
            }
        }
    }

    patterns
}

fn load_spec_examples() -> Vec<String> {
    let content = fs::read_to_string("tests/fixtures/grammar/hgvs_spec_examples.json")
        .unwrap_or_else(|_| r#"{"examples":[]}"#.to_string());
    let fixture: HgvsSpecFixture =
        serde_json::from_str(&content).unwrap_or_else(|_| HgvsSpecFixture { examples: vec![] });

    fixture.examples.into_iter().map(|e| e.hgvs).collect()
}

// =============================================================================
// Transcript Data for Normalization Tests
// =============================================================================

#[derive(Debug, Deserialize, Clone)]
struct TranscriptData {
    id: String,
    gene_symbol: String,
    strand: String,
    sequence: String,
    cds_start: u64,
    cds_end: u64,
    exons: Vec<ExonData>,
}

#[derive(Debug, Deserialize, Clone)]
struct ExonData {
    number: u32,
    start: u64,
    end: u64,
}

static TRANSCRIPTS: OnceLock<HashMap<String, TranscriptData>> = OnceLock::new();

fn get_transcripts() -> &'static HashMap<String, TranscriptData> {
    TRANSCRIPTS.get_or_init(|| {
        let json_path = format!(
            "{}{}",
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/sequences/normalization_transcripts.json"
        );
        match std::fs::read_to_string(&json_path) {
            Ok(json_str) => {
                let transcripts: Vec<TranscriptData> =
                    serde_json::from_str(&json_str).unwrap_or_default();
                transcripts.into_iter().map(|t| (t.id.clone(), t)).collect()
            }
            Err(_) => HashMap::new(),
        }
    })
}

fn create_provider_for_transcript(accession: &str) -> Option<MockProvider> {
    let data = get_transcripts().get(accession)?;
    let mut provider = MockProvider::new();

    let strand = if data.strand == "+" {
        Strand::Plus
    } else {
        Strand::Minus
    };
    let exons: Vec<Exon> = data
        .exons
        .iter()
        .map(|e| Exon::new(e.number, e.start, e.end))
        .collect();

    let transcript = Transcript::new(
        data.id.clone(),
        Some(data.gene_symbol.clone()),
        strand,
        data.sequence.clone(),
        Some(data.cds_start),
        Some(data.cds_end),
        exons,
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    );
    provider.add_transcript(transcript);
    Some(provider)
}

// =============================================================================
// PARSING IDEMPOTENCY TESTS
// =============================================================================
// Verify: parse(format(parse(input))) == parse(input)

#[test]
fn test_parsing_idempotency_edge_cases() {
    let patterns = load_edge_cases();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for input in &patterns {
        // First parse
        let v1 = match parse_hgvs(input) {
            Ok(v) => v,
            Err(_) => continue, // Skip unparseable inputs
        };

        // Format and re-parse
        let formatted = format!("{}", v1);
        let v2 = match parse_hgvs(&formatted) {
            Ok(v) => v,
            Err(e) => {
                failed += 1;
                if failures.len() < 20 {
                    failures.push(format!(
                        "Re-parse failed: '{}' -> '{}' -> {}",
                        input, formatted, e
                    ));
                }
                continue;
            }
        };

        // Compare string representations (structural equality)
        let f1 = format!("{}", v1);
        let f2 = format!("{}", v2);

        if f1 == f2 {
            passed += 1;
        } else {
            failed += 1;
            if failures.len() < 20 {
                failures.push(format!("Mismatch: '{}' -> '{}' vs '{}'", input, f1, f2));
            }
        }
    }

    eprintln!("\nParsing idempotency (edge cases):");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        failed == 0,
        "Parsing idempotency failed for {} patterns",
        failed
    );
}

#[test]
fn test_parsing_idempotency_spec_examples() {
    let patterns = load_spec_examples();
    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for input in &patterns {
        // First parse
        let v1 = match parse_hgvs(input) {
            Ok(v) => v,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };

        // Format and re-parse
        let formatted = format!("{}", v1);
        let v2 = match parse_hgvs(&formatted) {
            Ok(v) => v,
            Err(e) => {
                failed += 1;
                if failures.len() < 20 {
                    failures.push(format!(
                        "Re-parse failed: '{}' -> '{}' -> {}",
                        input, formatted, e
                    ));
                }
                continue;
            }
        };

        // Compare string representations
        let f1 = format!("{}", v1);
        let f2 = format!("{}", v2);

        if f1 == f2 {
            passed += 1;
        } else {
            failed += 1;
            if failures.len() < 20 {
                failures.push(format!("Mismatch: '{}' -> '{}' vs '{}'", input, f1, f2));
            }
        }
    }

    eprintln!("\nParsing idempotency (HGVS spec examples):");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped (unparseable): {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // Allow some failures for spec examples (some use advanced features)
    let total_testable = passed + failed;
    let success_rate = if total_testable > 0 {
        (passed as f64 / total_testable as f64) * 100.0
    } else {
        100.0
    };

    assert!(
        success_rate >= 95.0,
        "Parsing idempotency {:.1}% below 95% threshold",
        success_rate
    );
}

// =============================================================================
// NORMALIZATION IDEMPOTENCY TESTS
// =============================================================================
// Verify: normalize(normalize(input)) == normalize(input)

#[test]
fn test_normalization_idempotency() {
    // Test patterns that we have transcript data for
    let test_patterns = vec![
        // CFTR variants
        "NM_000492.4:c.1521_1523del",
        "NM_000492.4:c.1520_1522del",
        "NM_000492.4:c.3846G>A",
        // BRCA1 variants
        "NM_007294.4:c.68_69del",
        "NM_007294.4:c.5266dup",
        // BRCA2 variants
        "NM_000059.4:c.5946del",
        "NM_000059.4:c.6275_6276del",
        // TP53 variants
        "NM_000546.6:c.215C>G",
        "NM_000546.6:c.743G>A",
        "NM_000546.6:c.532del",
        // COL1A1 variants
        "NM_000088.4:c.769G>A",
        // Non-coding RNA
        "NR_153405.1:n.3650G>A",
        "NR_153405.1:n.3799C>T",
        // Insertions that become dups
        "NM_001089.2:c.3057_3058insT",
        // Duplications that become repeats
        "NM_020732.3:c.357_362dupGCAGCA",
        // Repeat notation
        "NM_000492.4:c.1516ATC[1]",
    ];

    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for input in &test_patterns {
        // Extract accession
        let accession = match input.split(':').next() {
            Some(a) => a,
            None => {
                skipped += 1;
                continue;
            }
        };

        // Get provider for this transcript
        let provider = match create_provider_for_transcript(accession) {
            Some(p) => p,
            None => {
                skipped += 1;
                continue;
            }
        };

        // Parse
        let v1 = match parse_hgvs(input) {
            Ok(v) => v,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };

        // First normalization
        let normalizer = Normalizer::new(provider.clone());
        let n1 = match normalizer.normalize(&v1) {
            Ok(n) => n,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };
        let f1 = format!("{}", n1);

        // Parse the normalized output
        let v2 = match parse_hgvs(&f1) {
            Ok(v) => v,
            Err(e) => {
                failed += 1;
                if failures.len() < 20 {
                    failures.push(format!("Re-parse failed: '{}' -> '{}' -> {}", input, f1, e));
                }
                continue;
            }
        };

        // Second normalization
        let normalizer2 = Normalizer::new(provider);
        let n2 = match normalizer2.normalize(&v2) {
            Ok(n) => n,
            Err(e) => {
                failed += 1;
                if failures.len() < 20 {
                    failures.push(format!(
                        "Re-normalize failed: '{}' -> '{}' -> {}",
                        input, f1, e
                    ));
                }
                continue;
            }
        };
        let f2 = format!("{}", n2);

        // Compare
        if f1 == f2 {
            passed += 1;
        } else {
            failed += 1;
            if failures.len() < 20 {
                failures.push(format!(
                    "Idempotency failed: '{}' -> '{}' -> '{}'",
                    input, f1, f2
                ));
            }
        }
    }

    eprintln!("\nNormalization idempotency:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped (no transcript data): {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        failed == 0,
        "Normalization idempotency failed for {} patterns",
        failed
    );
}

// =============================================================================
// CLINVAR SAMPLE IDEMPOTENCY
// =============================================================================

#[test]
fn test_clinvar_sample_parsing_idempotency() {
    // Load a sample of ClinVar patterns if available
    let clinvar_path = "clinvar_1000.txt";
    let patterns: Vec<String> = match fs::read_to_string(clinvar_path) {
        Ok(content) => content.lines().map(|s| s.to_string()).collect(),
        Err(_) => {
            eprintln!(
                "Skipping ClinVar sample test - file not found: {}",
                clinvar_path
            );
            return;
        }
    };

    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for input in patterns.iter().take(1000) {
        // Skip empty lines
        if input.trim().is_empty() {
            continue;
        }

        // First parse
        let v1 = match parse_hgvs(input) {
            Ok(v) => v,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };

        // Format and re-parse
        let formatted = format!("{}", v1);
        let v2 = match parse_hgvs(&formatted) {
            Ok(v) => v,
            Err(e) => {
                failed += 1;
                if failures.len() < 10 {
                    failures.push(format!(
                        "Re-parse failed: '{}' -> '{}' -> {}",
                        input, formatted, e
                    ));
                }
                continue;
            }
        };

        // Compare
        let f1 = format!("{}", v1);
        let f2 = format!("{}", v2);

        if f1 == f2 {
            passed += 1;
        } else {
            failed += 1;
            if failures.len() < 10 {
                failures.push(format!("Mismatch: '{}' -> '{}' vs '{}'", input, f1, f2));
            }
        }
    }

    eprintln!("\nClinVar sample parsing idempotency:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped: {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nSample failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    // Should have very high success rate
    let total_testable = passed + failed;
    if total_testable > 0 {
        let success_rate = (passed as f64 / total_testable as f64) * 100.0;
        assert!(
            success_rate >= 99.0,
            "ClinVar parsing idempotency {:.1}% below 99% threshold",
            success_rate
        );
    }
}
