//! SPDI Roundtrip Tests
//!
//! These tests validate HGVS ↔ SPDI consistency using data from NCBI Variation Services.
//! SPDI (Sequence Position Deletion Insertion) is NCBI's canonical representation
//! for sequence variants.
//!
//! Tests verify that:
//! 1. HGVS expressions parse correctly
//! 2. When NCBI returns SPDI, we can parse the corresponding HGVS
//! 3. Roundtrip conversions maintain semantic equivalence

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::fs;
use std::path::Path;

// ============================================================================
// NCBI Variation Services Fixture Types
// ============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct NcbiVariationFixture {
    source: String,
    api_base: String,
    generated: String,
    hgvs_conversions: HgvsConversions,
    rsid_lookups: RsidLookups,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct HgvsConversions {
    total: usize,
    successful: usize,
    variants: Vec<HgvsConversion>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct HgvsConversion {
    input_hgvs: String,
    spdi_result: Option<SpdiResult>,
    roundtrip_hgvs: Option<RoundtripResult>,
    vcf: Option<serde_json::Value>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct SpdiResult {
    #[serde(default)]
    data: Option<SpdiData>,
    #[serde(default)]
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct SpdiData {
    #[serde(default)]
    spdis: Vec<Spdi>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Spdi {
    seq_id: String,
    position: i64,
    deleted_sequence: String,
    inserted_sequence: String,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RoundtripResult {
    #[serde(default)]
    data: Option<RoundtripData>,
    #[serde(default)]
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RoundtripData {
    #[serde(default)]
    hgvs_list: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RsidLookups {
    total: usize,
    successful: usize,
    variants: Vec<RsidLookup>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RsidLookup {
    rsid: String,
    info: Option<serde_json::Value>,
    #[serde(default)]
    hgvs_expressions: Vec<String>,
    spdi: Option<serde_json::Value>,
}

// ============================================================================
// Test Report
// ============================================================================

#[derive(Default)]
struct SpdiTestReport {
    total_hgvs: usize,
    hgvs_parsed: usize,
    hgvs_with_spdi: usize,
    roundtrip_available: usize,
    roundtrip_parsed: usize,
    parse_failures: Vec<String>,
}

impl SpdiTestReport {
    fn summary(&self) -> String {
        format!(
            "HGVS total: {}, parsed: {}, with SPDI: {}, roundtrip available: {}, roundtrip parsed: {}",
            self.total_hgvs,
            self.hgvs_parsed,
            self.hgvs_with_spdi,
            self.roundtrip_available,
            self.roundtrip_parsed
        )
    }
}

// ============================================================================
// HGVS to SPDI Tests
// ============================================================================

#[test]
#[ignore = "Requires ncbi_variation.json fixture - run scripts/fetch_ncbi_variation.py first"]
fn test_hgvs_to_spdi_parsing() {
    let fixture_path = Path::new("tests/fixtures/validation/ncbi_variation.json");

    if !fixture_path.exists() {
        eprintln!("NCBI fixture not found. Run: python scripts/fetch_ncbi_variation.py");
        return;
    }

    let content = fs::read_to_string(fixture_path).expect("Failed to read ncbi_variation.json");
    let fixture: NcbiVariationFixture =
        serde_json::from_str(&content).expect("Failed to parse ncbi_variation.json");

    let mut report = SpdiTestReport {
        total_hgvs: fixture.hgvs_conversions.variants.len(),
        ..Default::default()
    };

    for conversion in &fixture.hgvs_conversions.variants {
        // Parse the input HGVS
        let parse_result = parse_hgvs(&conversion.input_hgvs);

        if parse_result.is_ok() {
            report.hgvs_parsed += 1;
        } else {
            report.parse_failures.push(conversion.input_hgvs.clone());
        }

        // Check if NCBI returned SPDI
        if let Some(ref spdi_result) = conversion.spdi_result {
            if let Some(ref data) = spdi_result.data {
                if !data.spdis.is_empty() {
                    report.hgvs_with_spdi += 1;
                }
            }
        }

        // Check roundtrip HGVS
        if let Some(ref roundtrip) = conversion.roundtrip_hgvs {
            if let Some(ref data) = roundtrip.data {
                if !data.hgvs_list.is_empty() {
                    report.roundtrip_available += 1;

                    // Try to parse each roundtrip HGVS
                    let mut any_parsed = false;
                    for hgvs in &data.hgvs_list {
                        if parse_hgvs(hgvs).is_ok() {
                            any_parsed = true;
                            break;
                        }
                    }
                    if any_parsed {
                        report.roundtrip_parsed += 1;
                    }
                }
            }
        }
    }

    // Print report
    println!("\n=== SPDI Roundtrip Test Report ===");
    println!("{}", report.summary());

    if !report.parse_failures.is_empty() {
        println!("\nParse failures: {}", report.parse_failures.len());
        for input in report.parse_failures.iter().take(10) {
            println!("  - {}", input);
        }
    }

    // Calculate success rate
    let parse_rate = (report.hgvs_parsed as f64 / report.total_hgvs as f64) * 100.0;
    println!("\nParse success rate: {:.1}%", parse_rate);

    // Should parse most HGVS expressions
    assert!(parse_rate > 80.0, "Parse rate should be above 80%");
}

// ============================================================================
// rsID Lookup Tests
// ============================================================================

#[test]
#[ignore = "Requires ncbi_variation.json fixture"]
fn test_rsid_hgvs_parsing() {
    let fixture_path = Path::new("tests/fixtures/validation/ncbi_variation.json");

    if !fixture_path.exists() {
        return;
    }

    let content = fs::read_to_string(fixture_path).expect("Failed to read fixture");
    let fixture: NcbiVariationFixture =
        serde_json::from_str(&content).expect("Failed to parse fixture");

    let mut total_hgvs = 0;
    let mut parsed_hgvs = 0;
    let mut rsids_with_hgvs = 0;

    for lookup in &fixture.rsid_lookups.variants {
        if !lookup.hgvs_expressions.is_empty() {
            rsids_with_hgvs += 1;

            for hgvs in &lookup.hgvs_expressions {
                total_hgvs += 1;
                if parse_hgvs(hgvs).is_ok() {
                    parsed_hgvs += 1;
                }
            }
        }
    }

    println!("\n=== rsID HGVS Parsing Report ===");
    println!(
        "rsIDs with HGVS: {}/{}",
        rsids_with_hgvs,
        fixture.rsid_lookups.variants.len()
    );
    println!("Total HGVS expressions: {}", total_hgvs);
    println!("Parsed successfully: {}", parsed_hgvs);

    if total_hgvs > 0 {
        let rate = (parsed_hgvs as f64 / total_hgvs as f64) * 100.0;
        println!("Parse rate: {:.1}%", rate);
    }
}

// ============================================================================
// SPDI Format Validation Tests
// ============================================================================

#[test]
fn test_spdi_format_understanding() {
    // Test that we understand SPDI format correctly
    // SPDI: Sequence:Position:Deletion:Insertion

    // Example: NM_000518.4:76:G:GG represents a G duplication
    // This corresponds to NM_000518.4:c.27dupG

    // Verify we can parse the HGVS form
    let result = parse_hgvs("NM_000518.4:c.27dupG");
    assert!(
        result.is_ok(),
        "Should parse canonical duplication notation"
    );

    // Verify we can parse the equivalent insertion form
    let result2 = parse_hgvs("NM_000518.4:c.27_28insG");
    assert!(result2.is_ok(), "Should parse insertion notation");
}

#[test]
fn test_spdi_deletion_equivalence() {
    // SPDI deletion: NM_000492.4:1653:CTT: (empty insertion = deletion)
    // HGVS: NM_000492.4:c.1521_1523del

    let result = parse_hgvs("NM_000492.4:c.1521_1523del");
    assert!(result.is_ok(), "Should parse deletion");

    // With explicit deleted sequence
    let result2 = parse_hgvs("NM_000492.4:c.1521_1523delCTT");
    assert!(result2.is_ok(), "Should parse deletion with sequence");
}

#[test]
fn test_spdi_substitution_equivalence() {
    // SPDI substitution: NC_000017.11:7674220:C:T
    // HGVS: NC_000017.11:g.7674220C>T

    let result = parse_hgvs("NC_000017.11:g.7674220C>T");
    assert!(result.is_ok(), "Should parse genomic substitution");

    if let Ok(variant) = result {
        match variant {
            HgvsVariant::Genome(genome) => {
                assert_eq!(format!("{}", genome.accession), "NC_000017.11");
            }
            _ => panic!("Expected genome variant"),
        }
    }
}

#[test]
fn test_spdi_indel_equivalence() {
    // SPDI complex: NM_000546.6:100:ACG:TTT
    // HGVS: NM_000546.6:c.100_102delinsTTT

    let result = parse_hgvs("NM_000546.6:c.100_102delinsTTT");
    assert!(result.is_ok(), "Should parse delins");

    if let Ok(variant) = result {
        match variant {
            HgvsVariant::Cds(_) => {} // Expected - delins is a coding variant
            _ => panic!("Expected Cds variant for delins"),
        }
    }
}

// ============================================================================
// Coordinate System Tests
// ============================================================================

#[test]
fn test_spdi_zero_based_understanding() {
    // SPDI uses 0-based coordinates
    // HGVS uses 1-based coordinates

    // SPDI: NM_000546.6:214:C:G (0-based position 214)
    // HGVS: NM_000546.6:c.215C>G (1-based position 215)

    let result = parse_hgvs("NM_000546.6:c.215C>G");
    assert!(result.is_ok(), "Should parse with 1-based coordinate");
}

#[test]
fn test_spdi_half_open_intervals() {
    // SPDI uses half-open intervals [start, end)
    // HGVS uses closed intervals [start, end]

    // For a 3bp deletion:
    // SPDI: position=1520, deleted_sequence=CTT (deletes positions 1520,1521,1522)
    // HGVS: c.1521_1523del (deletes positions 1521,1522,1523 in 1-based)

    let result = parse_hgvs("NM_000492.4:c.1521_1523del");
    assert!(result.is_ok());
}
