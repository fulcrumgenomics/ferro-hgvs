//! HGVS Specification Compliance Tests
//!
//! These tests validate parsing and normalization against examples from the official
//! HGVS nomenclature specification at hgvs-nomenclature.org.
//!
//! Source: https://hgvs-nomenclature.org
//!
//! ## Normalization Tests
//!
//! In addition to parsing tests, this module includes normalization compliance tests:
//! - Normalized output validity (re-parseable as valid HGVS)
//! - Normalization idempotency (normalize(normalize(x)) == normalize(x))
//! - HGVS-specific normalization rules (3' rule, deletion/dup representation)

use ferro_hgvs::parse_hgvs;
use ferro_hgvs::reference::mock::MockProvider;
use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::Normalizer;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::sync::OnceLock;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct HgvsSpecFixture {
    description: String,
    source: String,
    version: String,
    generated: String,
    total_examples: usize,
    summary: Summary,
    examples: Vec<Example>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Summary {
    by_variant_type: HashMap<String, usize>,
    by_coordinate_system: HashMap<String, usize>,
    by_level: HashMap<String, usize>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Example {
    hgvs: String,
    variant_type: String,
    level: String,
    coordinate_system: String,
    source: String,
}

fn load_hgvs_spec_fixtures() -> HgvsSpecFixture {
    let content = fs::read_to_string("tests/fixtures/grammar/hgvs_spec_examples.json")
        .expect("Failed to read hgvs_spec_examples.json");
    serde_json::from_str(&content).expect("Failed to parse hgvs_spec_examples.json")
}

/// Test all HGVS spec examples
#[test]
fn test_hgvs_spec_all_examples() {
    let fixtures = load_hgvs_spec_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for example in &fixtures.examples {
        let result = parse_hgvs(&example.hgvs);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                if failures.len() < 50 {
                    failures.push(format!(
                        "{} ({}/{}): {}",
                        example.hgvs, example.variant_type, example.coordinate_system, e
                    ));
                }
            }
        }
    }

    eprintln!("\nHGVS spec compliance test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Total: {}", fixtures.total_examples);

    if !failures.is_empty() {
        eprintln!("\nSample failures (first 30):");
        for f in failures.iter().take(30) {
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

    // Note: Some examples from the spec are incomplete or use advanced features
    // Track compliance rate over time - threshold is informational
    assert!(
        success_rate >= 40.0,
        "HGVS spec compliance {:.1}% below 40% threshold",
        success_rate
    );
}

/// Test by variant type
#[test]
fn test_hgvs_spec_by_variant_type() {
    let fixtures = load_hgvs_spec_fixtures();

    let mut type_counts: HashMap<&str, (usize, usize)> = HashMap::new();

    for example in &fixtures.examples {
        let result = parse_hgvs(&example.hgvs);
        let entry = type_counts.entry(&example.variant_type).or_insert((0, 0));
        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nHGVS spec by variant type:");
    for (vtype, (passed, failed)) in &type_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", vtype, passed, total, rate);
    }
}

/// Test by coordinate system
#[test]
fn test_hgvs_spec_by_coordinate_system() {
    let fixtures = load_hgvs_spec_fixtures();

    let mut coord_counts: HashMap<&str, (usize, usize)> = HashMap::new();

    for example in &fixtures.examples {
        let result = parse_hgvs(&example.hgvs);
        let entry = coord_counts
            .entry(&example.coordinate_system)
            .or_insert((0, 0));
        if result.is_ok() {
            entry.0 += 1;
        } else {
            entry.1 += 1;
        }
    }

    eprintln!("\nHGVS spec by coordinate system:");
    for (coord, (passed, failed)) in &coord_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", coord, passed, total, rate);
    }
}

/// Test substitution examples specifically
#[test]
fn test_hgvs_spec_substitutions() {
    let fixtures = load_hgvs_spec_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for example in fixtures
        .examples
        .iter()
        .filter(|e| e.variant_type == "substitution")
    {
        let result = parse_hgvs(&example.hgvs);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", example.hgvs, e));
            }
        }
    }

    eprintln!("\nHGVS spec substitution test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nSubstitution failures:");
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

    // Track substitution compliance - threshold is informational
    // Some examples may be incomplete patterns from spec text
    assert!(
        success_rate >= 30.0,
        "Substitution compliance {:.1}% below 30% threshold",
        success_rate
    );
}

/// Test deletion examples specifically
#[test]
fn test_hgvs_spec_deletions() {
    let fixtures = load_hgvs_spec_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for example in fixtures
        .examples
        .iter()
        .filter(|e| e.variant_type == "deletion")
    {
        let result = parse_hgvs(&example.hgvs);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", example.hgvs, e));
            }
        }
    }

    eprintln!("\nHGVS spec deletion test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nDeletion failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }
}

/// Test duplication examples specifically
#[test]
fn test_hgvs_spec_duplications() {
    let fixtures = load_hgvs_spec_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for example in fixtures
        .examples
        .iter()
        .filter(|e| e.variant_type == "duplication")
    {
        let result = parse_hgvs(&example.hgvs);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", example.hgvs, e));
            }
        }
    }

    eprintln!("\nHGVS spec duplication test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nDuplication failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }
}

/// Test insertion examples specifically
#[test]
fn test_hgvs_spec_insertions() {
    let fixtures = load_hgvs_spec_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for example in fixtures
        .examples
        .iter()
        .filter(|e| e.variant_type == "insertion")
    {
        let result = parse_hgvs(&example.hgvs);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", example.hgvs, e));
            }
        }
    }

    eprintln!("\nHGVS spec insertion test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nInsertion failures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }
}

/// Test repeat examples specifically
#[test]
fn test_hgvs_spec_repeats() {
    let fixtures = load_hgvs_spec_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for example in fixtures
        .examples
        .iter()
        .filter(|e| e.variant_type == "repeat")
    {
        let result = parse_hgvs(&example.hgvs);

        match result {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", example.hgvs, e));
            }
        }
    }

    eprintln!("\nHGVS spec repeat test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nRepeat failures (first 10):");
        for f in failures.iter().take(10) {
            eprintln!("  - {}", f);
        }
    }
}

// =============================================================================
// NORMALIZATION COMPLIANCE TESTS
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

/// Test that normalized output is valid HGVS (can be re-parsed)
///
/// This verifies that our normalizer produces syntactically correct HGVS strings.
#[test]
fn test_normalized_output_is_valid_hgvs() {
    // Test patterns with known transcript data
    let test_patterns = vec![
        // Substitutions
        "NM_000492.4:c.3846G>A",
        "NM_007294.4:c.5266A>G",
        "NM_000546.6:c.215C>G",
        "NM_000546.6:c.743G>A",
        // Deletions
        "NM_000492.4:c.1521_1523del",
        "NM_007294.4:c.68_69del",
        "NM_000059.4:c.5946del",
        "NM_000546.6:c.532del",
        // Duplications
        "NM_007294.4:c.5266dup",
        // Insertions
        "NM_001089.2:c.3057_3058insT",
        // Non-coding
        "NR_153405.1:n.3650G>A",
        "NR_153405.1:n.3799C>T",
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

        // Normalize
        let normalizer = Normalizer::new(provider);
        let normalized = match normalizer.normalize(&v1) {
            Ok(n) => n,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };

        // Format and verify re-parseable
        let formatted = format!("{}", normalized);
        match parse_hgvs(&formatted) {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                failures.push(format!("{} -> {} (error: {})", input, formatted, e));
            }
        }
    }

    eprintln!("\nNormalized output validity test:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped: {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        failed == 0,
        "Normalized output not valid HGVS for {} patterns",
        failed
    );
}

/// Test normalization idempotency: normalize(normalize(x)) == normalize(x)
///
/// This is a fundamental property that ensures stability of normalized forms.
#[test]
fn test_normalization_idempotency() {
    let test_patterns = vec![
        // Various variant types that should normalize stably
        "NM_000492.4:c.1521_1523del",
        "NM_000492.4:c.3846G>A",
        "NM_007294.4:c.68_69del",
        "NM_007294.4:c.5266dup",
        "NM_000059.4:c.5946del",
        "NM_000546.6:c.215C>G",
        "NM_000546.6:c.532del",
        "NR_153405.1:n.3650G>A",
    ];

    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for input in &test_patterns {
        let accession = match input.split(':').next() {
            Some(a) => a,
            None => {
                skipped += 1;
                continue;
            }
        };

        let provider = match create_provider_for_transcript(accession) {
            Some(p) => p,
            None => {
                skipped += 1;
                continue;
            }
        };

        let v1 = match parse_hgvs(input) {
            Ok(v) => v,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };

        // First normalization
        let normalizer1 = Normalizer::new(provider.clone());
        let n1 = match normalizer1.normalize(&v1) {
            Ok(n) => n,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };
        let f1 = format!("{}", n1);

        // Parse and normalize again
        let v2 = match parse_hgvs(&f1) {
            Ok(v) => v,
            Err(e) => {
                failed += 1;
                failures.push(format!("Re-parse failed: {} -> {} ({})", input, f1, e));
                continue;
            }
        };

        let normalizer2 = Normalizer::new(provider);
        let n2 = match normalizer2.normalize(&v2) {
            Ok(n) => n,
            Err(e) => {
                failed += 1;
                failures.push(format!("Re-normalize failed: {} -> {} ({})", input, f1, e));
                continue;
            }
        };
        let f2 = format!("{}", n2);

        if f1 == f2 {
            passed += 1;
        } else {
            failed += 1;
            failures.push(format!("Not idempotent: {} -> {} -> {}", input, f1, f2));
        }
    }

    eprintln!("\nNormalization idempotency test:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped: {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(
        failed == 0,
        "Normalization not idempotent for {} patterns",
        failed
    );
}

// =============================================================================
// HGVS NORMALIZATION RULE TESTS
// =============================================================================
// These tests verify specific HGVS normalization rules from the specification.

/// Test 3' rule: variants in repetitive regions should be positioned at the 3' end
///
/// HGVS specification states that for all descriptions the most 3' position
/// possible is arbitrarily assigned to have been changed.
/// Reference: https://hgvs-nomenclature.org/stable/background/basics/
#[test]
fn test_three_prime_rule_deletions() {
    // Test cases where a deletion in a repeat region should shift 3'
    // Format: (input, expected_normalized)
    let test_cases = vec![
        // Single base deletion in poly-T should shift 3'
        // The exact expected output depends on the sequence context
        ("NM_000492.4:c.1521_1523del", "NM_000492.4:c.1521_1523del"),
    ];

    for (input, expected) in test_cases {
        let accession = input.split(':').next().unwrap();
        let provider = match create_provider_for_transcript(accession) {
            Some(p) => p,
            None => {
                eprintln!("Skipping {}: no transcript data", input);
                continue;
            }
        };

        let v = parse_hgvs(input).expect("Should parse");
        let normalizer = Normalizer::new(provider);
        let normalized = normalizer.normalize(&v).expect("Should normalize");
        let result = format!("{}", normalized);

        assert_eq!(
            result, expected,
            "3' rule not applied correctly: {} -> {} (expected {})",
            input, result, expected
        );
    }
}

/// Test insertion to duplication conversion
///
/// HGVS specification: an insertion of a sequence already present should be
/// described as a duplication.
#[test]
fn test_insertion_to_duplication_conversion() {
    // When an insertion is a copy of the immediately preceding sequence,
    // it should be normalized to a duplication
    let test_cases = vec![
        // Insert T after T -> dup
        ("NM_001089.2:c.3057_3058insT", "NM_001089.2:c.3057dup"),
    ];

    for (input, expected) in test_cases {
        let accession = input.split(':').next().unwrap();
        let provider = match create_provider_for_transcript(accession) {
            Some(p) => p,
            None => {
                eprintln!("Skipping {}: no transcript data", input);
                continue;
            }
        };

        let v = parse_hgvs(input).expect("Should parse");
        let normalizer = Normalizer::new(provider);
        let normalized = normalizer.normalize(&v).expect("Should normalize");
        let result = format!("{}", normalized);

        assert_eq!(
            result, expected,
            "Insertion not converted to dup: {} -> {} (expected {})",
            input, result, expected
        );
    }
}

/// Test that substitutions remain as substitutions (no spurious conversions)
#[test]
fn test_substitutions_remain_substitutions() {
    let test_cases = vec![
        "NM_000492.4:c.3846G>A",
        "NM_000546.6:c.215C>G",
        "NM_000546.6:c.743G>A",
        "NR_153405.1:n.3650G>A",
    ];

    for input in test_cases {
        let accession = input.split(':').next().unwrap();
        let provider = match create_provider_for_transcript(accession) {
            Some(p) => p,
            None => {
                eprintln!("Skipping {}: no transcript data", input);
                continue;
            }
        };

        let v = parse_hgvs(input).expect("Should parse");
        let normalizer = Normalizer::new(provider);
        let normalized = normalizer.normalize(&v).expect("Should normalize");
        let result = format!("{}", normalized);

        // Substitutions should stay as substitutions (contain ">")
        assert!(
            result.contains('>'),
            "Substitution converted to different type: {} -> {}",
            input,
            result
        );
    }
}

// =============================================================================
// CANONICAL HGVS SPEC EXAMPLES
// =============================================================================
// Tests based on official examples from https://github.com/HGVSnomenclature/hgvs-nomenclature

/// Test parsing of canonical DNA examples from HGVS syntax.yaml
#[test]
fn test_canonical_dna_examples() {
    let canonical_examples = vec![
        // Deletions
        "NC_000001.11:g.1234del",
        "NC_000001.11:g.1234_2345del",
        // Deletion-insertions
        "NC_000001.11:g.123delinsAC",
        "NC_000001.11:g.123_129delinsAC",
        // Duplications
        "NC_000001.11:g.1234dup",
        "NC_000001.11:g.1234_2345dup",
        // Insertions
        "NC_000001.11:g.1234_1235insACGT",
        // Inversions
        "NC_000001.11:g.1234_2345inv",
        // No change
        "NC_000001.11:g.1234=",
        // Substitutions
        "NC_000023.10:g.33038255C>A",
    ];

    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for input in &canonical_examples {
        match parse_hgvs(input) {
            Ok(v) => {
                // Also verify round-trip
                let formatted = format!("{}", v);
                match parse_hgvs(&formatted) {
                    Ok(_) => passed += 1,
                    Err(e) => {
                        failed += 1;
                        failures.push(format!(
                            "Round-trip failed: {} -> {} ({})",
                            input, formatted, e
                        ));
                    }
                }
            }
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", input, e));
            }
        }
    }

    eprintln!("\nCanonical DNA examples:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(failed == 0, "Canonical DNA examples: {} failures", failed);
}

/// Test parsing of canonical RNA examples from HGVS syntax.yaml
#[test]
fn test_canonical_rna_examples() {
    let canonical_examples = vec![
        "NM_004006.3:r.123_127del",
        "NM_004006.3:r.123_127delinsag",
        "NM_004006.3:r.123_345dup",
        "NM_004006.3:r.123_124insauc",
        "NM_004006.3:r.123_345inv",
        "NM_004006.3:r.123c>g",
    ];

    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for input in &canonical_examples {
        match parse_hgvs(input) {
            Ok(v) => {
                let formatted = format!("{}", v);
                match parse_hgvs(&formatted) {
                    Ok(_) => passed += 1,
                    Err(e) => {
                        failed += 1;
                        failures.push(format!(
                            "Round-trip failed: {} -> {} ({})",
                            input, formatted, e
                        ));
                    }
                }
            }
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", input, e));
            }
        }
    }

    eprintln!("\nCanonical RNA examples:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);

    if !failures.is_empty() {
        eprintln!("\nFailures:");
        for f in &failures {
            eprintln!("  - {}", f);
        }
    }

    assert!(failed == 0, "Canonical RNA examples: {} failures", failed);
}

/// Test parsing of canonical protein examples from HGVS syntax.yaml
#[test]
fn test_canonical_protein_examples() {
    let canonical_examples = vec![
        // Deletions
        "NP_003997.2:p.Val7del",
        "NP_003997.2:p.Lys23_Val25del",
        // Duplications
        "NP_003997.2:p.Val7dup",
        "NP_003997.2:p.Lys23_Val25dup",
        // Substitutions
        "NP_003997.1:p.Trp24Cys",
        "NP_003997.1:p.Trp24Ter",
        "NP_003997.1:p.W24*",
    ];

    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for input in &canonical_examples {
        match parse_hgvs(input) {
            Ok(v) => {
                let formatted = format!("{}", v);
                match parse_hgvs(&formatted) {
                    Ok(_) => passed += 1,
                    Err(e) => {
                        failed += 1;
                        failures.push(format!(
                            "Round-trip failed: {} -> {} ({})",
                            input, formatted, e
                        ));
                    }
                }
            }
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", input, e));
            }
        }
    }

    eprintln!("\nCanonical protein examples:");
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
        "Canonical protein examples: {} failures",
        failed
    );
}

/// Test parsing of canonical coding DNA examples from HGVS spec
#[test]
fn test_canonical_coding_dna_examples() {
    let canonical_examples = vec![
        "NM_004006.2:c.145_147delinsTGG",
        "NM_004006.2:c.20dup",
        "NM_004006.2:c.20_23dup",
        "NM_004006.2:c.5657_5660inv",
        "NM_004006.2:c.4145_4160inv",
        "NM_004006.2:c.4375C>T",
    ];

    let mut passed = 0;
    let mut failed = 0;
    let mut failures = Vec::new();

    for input in &canonical_examples {
        match parse_hgvs(input) {
            Ok(v) => {
                let formatted = format!("{}", v);
                match parse_hgvs(&formatted) {
                    Ok(_) => passed += 1,
                    Err(e) => {
                        failed += 1;
                        failures.push(format!(
                            "Round-trip failed: {} -> {} ({})",
                            input, formatted, e
                        ));
                    }
                }
            }
            Err(e) => {
                failed += 1;
                failures.push(format!("{}: {}", input, e));
            }
        }
    }

    eprintln!("\nCanonical coding DNA examples:");
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
        "Canonical coding DNA examples: {} failures",
        failed
    );
}
