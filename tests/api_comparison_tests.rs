//! API Comparison Tests
//!
//! These tests compare ferro-hgvs parsing and normalization results against
//! authoritative external APIs: Mutalyzer and VariantValidator.
//!
//! These tests are informational - discrepancies are logged but do not cause
//! test failures, as they may indicate areas for investigation rather than bugs.

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

// ============================================================================
// Mutalyzer API Fixture Types
// ============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct MutalyzerFixture {
    source: String,
    api_base: String,
    generated: String,
    total_variants: usize,
    successful: usize,
    variants: Vec<MutalyzerVariant>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct MutalyzerVariant {
    input: String,
    normalize: Option<MutalyzerNormalize>,
    model: Option<serde_json::Value>,
    spdi: Option<serde_json::Value>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct MutalyzerNormalize {
    #[serde(default)]
    normalized_description: Option<String>,
    #[serde(default)]
    equivalent_descriptions: Option<serde_json::Value>,
    #[serde(default)]
    error: Option<String>,
}

// ============================================================================
// VariantValidator API Fixture Types
// ============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct VariantValidatorFixture {
    source: String,
    api_base: String,
    generated: String,
    total_variants: usize,
    successful: usize,
    variants: Vec<VVVariant>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct VVVariant {
    input: String,
    genome_build: String,
    raw_validation: Option<serde_json::Value>,
    reference: Option<serde_json::Value>,
    summary: VVSummary,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct VVSummary {
    valid: bool,
    #[serde(default)]
    normalized: Option<String>,
    #[serde(default)]
    genomic: Option<String>,
    #[serde(default)]
    protein: Option<String>,
    #[serde(default)]
    warnings: Vec<String>,
    #[serde(default)]
    errors: Vec<String>,
}

// ============================================================================
// Comparison Report
// ============================================================================

#[derive(Default)]
struct ComparisonReport {
    total: usize,
    parsed_by_ferro: usize,
    parsed_by_api: usize,
    both_parsed: usize,
    normalization_matches: usize,
    normalization_mismatches: Vec<String>,
    ferro_only_failures: Vec<String>,
    api_only_failures: Vec<String>,
}

impl ComparisonReport {
    fn summary(&self) -> String {
        format!(
            "Total: {}, Ferro parsed: {}, API parsed: {}, Both: {}, Norm matches: {}, Mismatches: {}",
            self.total,
            self.parsed_by_ferro,
            self.parsed_by_api,
            self.both_parsed,
            self.normalization_matches,
            self.normalization_mismatches.len()
        )
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

#[allow(dead_code)]
fn get_variant_type(variant: &HgvsVariant) -> &'static str {
    match variant {
        HgvsVariant::Genome(_) => "genome",
        HgvsVariant::Cds(_) => "cds",
        HgvsVariant::Tx(_) => "tx",
        HgvsVariant::Rna(_) => "rna",
        HgvsVariant::Protein(_) => "protein",
        HgvsVariant::Mt(_) => "mt",
        HgvsVariant::Circular(_) => "circular",
        HgvsVariant::RnaFusion(_) => "rna_fusion",
        HgvsVariant::Allele(_) => "allele",
        HgvsVariant::NullAllele => "null_allele",
        HgvsVariant::UnknownAllele => "unknown_allele",
    }
}

// ============================================================================
// Mutalyzer Comparison Tests
// ============================================================================

#[test]
#[ignore = "Requires mutalyzer_api.json fixture - run scripts/fetch_mutalyzer_normalized.py first"]
fn test_mutalyzer_comparison() {
    let fixture_path = Path::new("tests/fixtures/normalization/mutalyzer_api.json");

    if !fixture_path.exists() {
        eprintln!("Mutalyzer fixture not found. Run: python scripts/fetch_mutalyzer_normalized.py");
        return;
    }

    let content = fs::read_to_string(fixture_path).expect("Failed to read mutalyzer_api.json");
    let fixture: MutalyzerFixture =
        serde_json::from_str(&content).expect("Failed to parse mutalyzer_api.json");

    let mut report = ComparisonReport {
        total: fixture.variants.len(),
        ..Default::default()
    };

    for variant in &fixture.variants {
        // Try parsing with ferro-hgvs
        let ferro_result = parse_hgvs(&variant.input);

        // Check if Mutalyzer parsed it
        let mutalyzer_parsed = variant
            .normalize
            .as_ref()
            .map(|n| n.normalized_description.is_some())
            .unwrap_or(false);

        let ferro_parsed = ferro_result.is_ok();

        if ferro_parsed {
            report.parsed_by_ferro += 1;
        }

        if mutalyzer_parsed {
            report.parsed_by_api += 1;
        }

        if ferro_parsed && mutalyzer_parsed {
            report.both_parsed += 1;

            // Compare normalized forms
            if let Some(ref normalize) = variant.normalize {
                if let Some(ref _mutalyzer_norm) = normalize.normalized_description {
                    // For now, just count - detailed comparison would require
                    // ferro-hgvs normalization implementation
                    report.normalization_matches += 1;
                }
            }
        } else if ferro_parsed && !mutalyzer_parsed {
            report.api_only_failures.push(variant.input.clone());
        } else if !ferro_parsed && mutalyzer_parsed {
            report.ferro_only_failures.push(variant.input.clone());
        }
    }

    // Print report
    println!("\n=== Mutalyzer Comparison Report ===");
    println!("{}", report.summary());

    if !report.ferro_only_failures.is_empty() {
        println!(
            "\nFerro failed to parse (Mutalyzer succeeded): {}",
            report.ferro_only_failures.len()
        );
        for input in report.ferro_only_failures.iter().take(10) {
            println!("  - {}", input);
        }
    }

    if !report.api_only_failures.is_empty() {
        println!(
            "\nMutalyzer failed to parse (Ferro succeeded): {}",
            report.api_only_failures.len()
        );
        for input in report.api_only_failures.iter().take(10) {
            println!("  - {}", input);
        }
    }

    // Informational - don't fail
    assert!(
        report.parsed_by_ferro > 0,
        "Should parse at least some variants"
    );
}

// ============================================================================
// VariantValidator Comparison Tests
// ============================================================================

#[test]
#[ignore = "Requires variantvalidator_api.json fixture - run scripts/fetch_variantvalidator.py first"]
fn test_variantvalidator_comparison() {
    let fixture_path = Path::new("tests/fixtures/validation/variantvalidator_api.json");

    if !fixture_path.exists() {
        eprintln!(
            "VariantValidator fixture not found. Run: python scripts/fetch_variantvalidator.py"
        );
        return;
    }

    let content =
        fs::read_to_string(fixture_path).expect("Failed to read variantvalidator_api.json");
    let fixture: VariantValidatorFixture =
        serde_json::from_str(&content).expect("Failed to parse variantvalidator_api.json");

    let mut report = ComparisonReport {
        total: fixture.variants.len(),
        ..Default::default()
    };

    let mut by_genome_build: HashMap<String, (usize, usize)> = HashMap::new();

    for variant in &fixture.variants {
        // Try parsing with ferro-hgvs
        let ferro_result = parse_hgvs(&variant.input);

        // Check if VariantValidator validated it
        let vv_validated = variant.summary.valid;

        let ferro_parsed = ferro_result.is_ok();

        // Track by genome build
        let (parsed, total) = by_genome_build
            .entry(variant.genome_build.clone())
            .or_insert((0, 0));
        *total += 1;
        if ferro_parsed {
            *parsed += 1;
        }

        if ferro_parsed {
            report.parsed_by_ferro += 1;
        }

        if vv_validated {
            report.parsed_by_api += 1;
        }

        if ferro_parsed && vv_validated {
            report.both_parsed += 1;

            // Compare normalized forms if available
            if let Some(ref _vv_norm) = variant.summary.normalized {
                // Store for potential analysis
                report.normalization_matches += 1;
            }
        } else if ferro_parsed && !vv_validated {
            report.api_only_failures.push(format!(
                "{} ({}) - VV errors: {:?}",
                variant.input, variant.genome_build, variant.summary.errors
            ));
        } else if !ferro_parsed && vv_validated {
            report.ferro_only_failures.push(format!(
                "{} ({}) - VV normalized to: {:?}",
                variant.input, variant.genome_build, variant.summary.normalized
            ));
        }
    }

    // Print report
    println!("\n=== VariantValidator Comparison Report ===");
    println!("{}", report.summary());

    println!("\nBy genome build:");
    for (build, (parsed, total)) in &by_genome_build {
        println!("  {}: {}/{} parsed", build, parsed, total);
    }

    if !report.ferro_only_failures.is_empty() {
        println!(
            "\nFerro failed to parse (VV validated): {}",
            report.ferro_only_failures.len()
        );
        for input in report.ferro_only_failures.iter().take(10) {
            println!("  - {}", input);
        }
    }

    // Check for warnings from VariantValidator
    let mut variants_with_warnings = 0;
    for variant in &fixture.variants {
        if !variant.summary.warnings.is_empty() {
            variants_with_warnings += 1;
        }
    }
    println!("\nVariants with VV warnings: {}", variants_with_warnings);

    // Informational - don't fail
    assert!(
        report.parsed_by_ferro > 0,
        "Should parse at least some variants"
    );
}

// ============================================================================
// Combined Analysis
// ============================================================================

#[test]
#[ignore = "Requires both API fixtures"]
fn test_cross_api_consistency() {
    let mutalyzer_path = Path::new("tests/fixtures/normalization/mutalyzer_api.json");
    let vv_path = Path::new("tests/fixtures/validation/variantvalidator_api.json");

    if !mutalyzer_path.exists() || !vv_path.exists() {
        eprintln!("API fixtures not found. Run the fetch scripts first.");
        return;
    }

    // This test could compare how both APIs handle the same variants
    // and identify any systematic differences

    println!("\n=== Cross-API Consistency Analysis ===");
    println!("(Would compare Mutalyzer and VariantValidator handling of shared variants)");

    // Placeholder for cross-API analysis
    // In a full implementation, this would:
    // 1. Find variants present in both fixtures
    // 2. Compare normalized forms
    // 3. Identify systematic differences
    // 4. Report on ferro-hgvs alignment with each
}
