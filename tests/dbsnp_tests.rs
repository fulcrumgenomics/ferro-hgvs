//! dbSNP rsID Validation Tests
//!
//! These tests validate that ferro-hgvs can correctly parse HGVS expressions
//! associated with clinically important dbSNP rsIDs.
//!
//! The test fixture contains curated rsID → HGVS mappings for variants that are:
//! - Clinically significant (pathogenic, drug response, etc.)
//! - Well-characterized in literature
//! - Commonly encountered in clinical genomics

use ferro_hgvs::{parse_hgvs, HgvsVariant};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

// ============================================================================
// Fixture Types
// ============================================================================

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct DbsnpFixture {
    source: String,
    description: String,
    generated: String,
    total_rsids: usize,
    rsids: Vec<RsidEntry>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct RsidEntry {
    rsid: String,
    gene: String,
    hgvs_c: Option<String>,
    hgvs_p: Option<String>,
    clinical_significance: String,
    condition: String,
    notes: String,
}

// ============================================================================
// Test Report
// ============================================================================

#[derive(Default)]
struct DbsnpTestReport {
    total_rsids: usize,
    hgvs_c_total: usize,
    hgvs_c_parsed: usize,
    hgvs_p_total: usize,
    hgvs_p_parsed: usize,
    by_gene: HashMap<String, (usize, usize)>, // (parsed, total)
    by_significance: HashMap<String, (usize, usize)>,
    failures: Vec<(String, String, String)>, // (rsid, hgvs, error)
}

impl DbsnpTestReport {
    fn add_hgvs_c(&mut self, entry: &RsidEntry, success: bool, error: Option<String>) {
        self.hgvs_c_total += 1;
        if success {
            self.hgvs_c_parsed += 1;
        } else if let Some(hgvs) = &entry.hgvs_c {
            self.failures
                .push((entry.rsid.clone(), hgvs.clone(), error.unwrap_or_default()));
        }

        let (parsed, total) = self.by_gene.entry(entry.gene.clone()).or_insert((0, 0));
        *total += 1;
        if success {
            *parsed += 1;
        }

        let (sig_parsed, sig_total) = self
            .by_significance
            .entry(entry.clinical_significance.clone())
            .or_insert((0, 0));
        *sig_total += 1;
        if success {
            *sig_parsed += 1;
        }
    }

    fn add_hgvs_p(&mut self, entry: &RsidEntry, success: bool, error: Option<String>) {
        self.hgvs_p_total += 1;
        if success {
            self.hgvs_p_parsed += 1;
        } else if let Some(hgvs) = &entry.hgvs_p {
            self.failures
                .push((entry.rsid.clone(), hgvs.clone(), error.unwrap_or_default()));
        }
    }

    fn summary(&self) -> String {
        format!(
            "rsIDs: {}, HGVSc: {}/{}, HGVSp: {}/{}",
            self.total_rsids,
            self.hgvs_c_parsed,
            self.hgvs_c_total,
            self.hgvs_p_parsed,
            self.hgvs_p_total
        )
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

fn load_dbsnp_fixtures() -> DbsnpFixture {
    let content = fs::read_to_string("tests/fixtures/validation/dbsnp_rsid.json")
        .expect("Failed to read dbsnp_rsid.json");
    serde_json::from_str(&content).expect("Failed to parse dbsnp_rsid.json")
}

#[allow(dead_code)]
fn get_variant_type_name(variant: &HgvsVariant) -> &'static str {
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
// Main Validation Test
// ============================================================================

#[test]
fn test_dbsnp_rsid_parsing() {
    let fixture = load_dbsnp_fixtures();
    let mut report = DbsnpTestReport {
        total_rsids: fixture.rsids.len(),
        ..Default::default()
    };

    for entry in &fixture.rsids {
        // Test coding HGVS
        if let Some(ref hgvs_c) = entry.hgvs_c {
            let result = parse_hgvs(hgvs_c);
            match result {
                Ok(_) => report.add_hgvs_c(entry, true, None),
                Err(e) => report.add_hgvs_c(entry, false, Some(e.to_string())),
            }
        }

        // Test protein HGVS
        if let Some(ref hgvs_p) = entry.hgvs_p {
            let result = parse_hgvs(hgvs_p);
            match result {
                Ok(_) => report.add_hgvs_p(entry, true, None),
                Err(e) => report.add_hgvs_p(entry, false, Some(e.to_string())),
            }
        }
    }

    // Print summary
    println!("\n=== dbSNP rsID Parsing Report ===");
    println!("{}", report.summary());

    // By gene
    println!("\nBy gene (top 10):");
    let mut genes: Vec<_> = report.by_gene.iter().collect();
    genes.sort_by(|a, b| b.1 .1.cmp(&a.1 .1));
    for (gene, (parsed, total)) in genes.iter().take(10) {
        println!("  {}: {}/{}", gene, parsed, total);
    }

    // By clinical significance
    println!("\nBy clinical significance:");
    for (sig, (parsed, total)) in &report.by_significance {
        println!("  {}: {}/{}", sig, parsed, total);
    }

    // Failures
    if !report.failures.is_empty() {
        println!("\nParse failures ({}):", report.failures.len());
        for (rsid, hgvs, error) in report.failures.iter().take(20) {
            println!("  {} ({}): {}", rsid, hgvs, error);
        }
    }

    // Calculate success rates
    let c_rate = if report.hgvs_c_total > 0 {
        (report.hgvs_c_parsed as f64 / report.hgvs_c_total as f64) * 100.0
    } else {
        100.0
    };
    let p_rate = if report.hgvs_p_total > 0 {
        (report.hgvs_p_parsed as f64 / report.hgvs_p_total as f64) * 100.0
    } else {
        100.0
    };

    println!("\nHGVSc parse rate: {:.1}%", c_rate);
    println!("HGVSp parse rate: {:.1}%", p_rate);

    // Should have high parse rates for curated clinical variants
    assert!(c_rate > 90.0, "HGVSc parse rate should be above 90%");
    assert!(p_rate > 85.0, "HGVSp parse rate should be above 85%");
}

// ============================================================================
// Gene-Specific Tests
// ============================================================================

#[test]
fn test_brca_variants() {
    let fixture = load_dbsnp_fixtures();

    let brca_entries: Vec<_> = fixture
        .rsids
        .iter()
        .filter(|e| e.gene == "BRCA1" || e.gene == "BRCA2")
        .collect();

    assert!(!brca_entries.is_empty(), "Should have BRCA variants");

    let mut passed = 0;
    let mut total = 0;

    for entry in brca_entries {
        if let Some(ref hgvs) = entry.hgvs_c {
            total += 1;
            if parse_hgvs(hgvs).is_ok() {
                passed += 1;
            }
        }
    }

    println!("BRCA variants: {}/{}", passed, total);
    assert_eq!(passed, total, "All BRCA variants should parse");
}

#[test]
fn test_tp53_variants() {
    let fixture = load_dbsnp_fixtures();

    let tp53_entries: Vec<_> = fixture.rsids.iter().filter(|e| e.gene == "TP53").collect();

    assert!(!tp53_entries.is_empty(), "Should have TP53 variants");

    for entry in tp53_entries {
        if let Some(ref hgvs) = entry.hgvs_c {
            let result = parse_hgvs(hgvs);
            assert!(
                result.is_ok(),
                "TP53 variant {} ({}) should parse: {:?}",
                entry.rsid,
                hgvs,
                result
            );
        }
    }
}

#[test]
fn test_pharmacogenomic_variants() {
    let fixture = load_dbsnp_fixtures();

    let pgx_entries: Vec<_> = fixture
        .rsids
        .iter()
        .filter(|e| e.clinical_significance.contains("Drug response"))
        .collect();

    if pgx_entries.is_empty() {
        return;
    }

    let mut passed = 0;

    for entry in &pgx_entries {
        if let Some(ref hgvs) = entry.hgvs_c {
            if parse_hgvs(hgvs).is_ok() {
                passed += 1;
            }
        }
    }

    println!("Pharmacogenomic variants: {}/{}", passed, pgx_entries.len());
}

// ============================================================================
// Variant Type Tests
// ============================================================================

#[test]
fn test_substitution_variants() {
    let fixture = load_dbsnp_fixtures();

    let substitutions: Vec<_> = fixture
        .rsids
        .iter()
        .filter(|e| e.hgvs_c.as_ref().map(|h| h.contains('>')).unwrap_or(false))
        .collect();

    for entry in substitutions {
        if let Some(ref hgvs) = entry.hgvs_c {
            let result = parse_hgvs(hgvs);
            if let Ok(variant) = result {
                // Substitutions should parse as Cds (c.) or Mt (m.) variants
                assert!(
                    matches!(variant, HgvsVariant::Cds(_) | HgvsVariant::Mt(_)),
                    "Expected CdsVariant or MtVariant for {}",
                    hgvs
                );
            }
        }
    }
}

#[test]
fn test_deletion_variants() {
    let fixture = load_dbsnp_fixtures();

    let deletions: Vec<_> = fixture
        .rsids
        .iter()
        .filter(|e| {
            e.hgvs_c
                .as_ref()
                .map(|h| h.contains("del") && !h.contains("delins"))
                .unwrap_or(false)
        })
        .collect();

    for entry in deletions {
        if let Some(ref hgvs) = entry.hgvs_c {
            let result = parse_hgvs(hgvs);
            if let Ok(variant) = result {
                // Coding deletions should parse as CdsVariant
                assert!(
                    matches!(variant, HgvsVariant::Cds(_)),
                    "Expected CdsVariant for {}",
                    hgvs
                );
            }
        }
    }
}

#[test]
fn test_frameshift_protein_variants() {
    let fixture = load_dbsnp_fixtures();

    let frameshifts: Vec<_> = fixture
        .rsids
        .iter()
        .filter(|e| e.hgvs_p.as_ref().map(|h| h.contains("fs")).unwrap_or(false))
        .collect();

    for entry in frameshifts {
        if let Some(ref hgvs) = entry.hgvs_p {
            let result = parse_hgvs(hgvs);
            assert!(
                result.is_ok(),
                "Frameshift {} ({}) should parse",
                entry.rsid,
                hgvs
            );
        }
    }
}

// ============================================================================
// Specific rsID Tests
// ============================================================================

#[test]
fn test_braf_v600e() {
    // rs121913529 - BRAF V600E
    let result = parse_hgvs("NM_004333.4:c.1799T>A");
    assert!(result.is_ok(), "BRAF V600E coding should parse");

    let result_p = parse_hgvs("NP_004324.2:p.Val600Glu");
    assert!(result_p.is_ok(), "BRAF V600E protein should parse");
}

#[test]
fn test_cftr_f508del() {
    // rs113993960 - CFTR F508del
    let result = parse_hgvs("NM_000492.4:c.1521_1523del");
    assert!(result.is_ok(), "CFTR F508del should parse");

    if let Ok(variant) = result {
        // Coding variant should parse as CdsVariant
        assert!(matches!(variant, HgvsVariant::Cds(_)));
    }
}

#[test]
fn test_sickle_cell() {
    // rs334 - HBB sickle cell
    let result = parse_hgvs("NM_000518.5:c.20A>T");
    assert!(result.is_ok(), "Sickle cell variant should parse");
}

#[test]
fn test_factor_v_leiden() {
    // rs6025 - Factor V Leiden
    let result = parse_hgvs("NM_000130.5:c.1601G>A");
    assert!(result.is_ok(), "Factor V Leiden should parse");
}

#[test]
fn test_mitochondrial_variants() {
    // Mitochondrial variants use m. coordinate system
    let fixture = load_dbsnp_fixtures();

    let mito_entries: Vec<_> = fixture
        .rsids
        .iter()
        .filter(|e| {
            e.hgvs_c
                .as_ref()
                .map(|h| h.contains(":m."))
                .unwrap_or(false)
        })
        .collect();

    for entry in mito_entries {
        if let Some(ref hgvs) = entry.hgvs_c {
            let result = parse_hgvs(hgvs);
            assert!(
                result.is_ok(),
                "Mitochondrial variant {} ({}) should parse",
                entry.rsid,
                hgvs
            );
        }
    }
}

// ============================================================================
// Coverage Statistics
// ============================================================================

#[test]
fn test_coverage_statistics() {
    let fixture = load_dbsnp_fixtures();

    let mut genes: HashMap<String, usize> = HashMap::new();
    let mut conditions: HashMap<String, usize> = HashMap::new();

    for entry in &fixture.rsids {
        *genes.entry(entry.gene.clone()).or_insert(0) += 1;
        *conditions.entry(entry.condition.clone()).or_insert(0) += 1;
    }

    println!("\n=== dbSNP Fixture Coverage ===");
    println!("Total rsIDs: {}", fixture.rsids.len());
    println!("Unique genes: {}", genes.len());
    println!("Unique conditions: {}", conditions.len());

    println!("\nTop genes:");
    let mut gene_list: Vec<_> = genes.iter().collect();
    gene_list.sort_by(|a, b| b.1.cmp(a.1));
    for (gene, count) in gene_list.iter().take(10) {
        println!("  {}: {}", gene, count);
    }
}
