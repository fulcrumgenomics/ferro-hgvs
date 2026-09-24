//! API Comparison Tests
//!
//! This test compares ferro-hgvs parsing results against VariantValidator,
//! using the live fetch the weekly External API Validation workflow writes to
//! `tests/fixtures/validation/variantvalidator_api.json`.
//!
//! It is informational: discrepancies are logged but do not fail the test, as
//! they may indicate areas for investigation rather than bugs. A missing fetch,
//! or one in which VariantValidator validated nothing, does fail it: the test
//! runs only in that workflow, right after the fetch, and would otherwise pass
//! having checked nothing.

use ferro_hgvs::parse_hgvs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

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
    api_normalized: usize,
    ferro_only_failures: Vec<String>,
    api_only_failures: Vec<String>,
}

impl ComparisonReport {
    fn summary(&self) -> String {
        format!(
            "Total: {}, Ferro parsed: {}, API parsed: {}, Both: {}, API normalized (not compared): {}",
            self.total,
            self.parsed_by_ferro,
            self.parsed_by_api,
            self.both_parsed,
            self.api_normalized
        )
    }
}

// ============================================================================
// VariantValidator Comparison Tests
// ============================================================================

#[test]
#[ignore = "Requires variantvalidator_api.json fixture - run scripts/fetch_variantvalidator.py first"]
fn test_variantvalidator_comparison() {
    let fixture_path = Path::new("tests/fixtures/validation/variantvalidator_api.json");

    assert!(
        fixture_path.exists(),
        "{} is missing: it is written by scripts/fetch_variantvalidator.py, which the \
         weekly External API Validation workflow runs before this test",
        fixture_path.display()
    );

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

            // Counted only: ferro's normalization is not compared with
            // VariantValidator's here.
            if variant.summary.normalized.is_some() {
                report.api_normalized += 1;
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

    if !report.api_only_failures.is_empty() {
        println!(
            "\nVariantValidator rejected (Ferro parsed): {}",
            report.api_only_failures.len()
        );
        for input in report.api_only_failures.iter().take(10) {
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

    // Discrepancies are informational, but a comparison must have happened: at
    // least one variant parsed by ferro AND validated by VariantValidator. A
    // fetch in which every API call failed (an upstream outage) fails here
    // rather than passing having compared nothing.
    assert!(
        report.both_parsed > 0,
        "no variant was both parsed by ferro and validated by VariantValidator \
         (ferro parsed {}, VariantValidator validated {} of {}), so nothing was compared; \
         an upstream outage, or a broken request in scripts/fetch_variantvalidator.py",
        report.parsed_by_ferro,
        report.parsed_by_api,
        report.total
    );
}
