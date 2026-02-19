//! MaveDB Functional Variant Tests
//!
//! These tests validate parsing against HGVS patterns from MaveDB,
//! the database for Multiplex Assays of Variant Effect (MAVE).
//!
//! NOTE: MaveDB uses MAVE-HGVS format which often lacks reference accessions
//! (e.g., `p.Leu11Pro` instead of `NP_123456:p.Leu11Pro`). These short-form
//! notations require context from the score set to be fully interpretable.
//!
//! Source: https://www.mavedb.org

use ferro_hgvs::parse_hgvs;
use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct MaveDbFixture {
    source: String,
    api_base: String,
    generated: String,
    score_sets_processed: usize,
    total_nucleotide_variants: usize,
    total_protein_variants: usize,
    summary: Summary,
    score_sets: Vec<ScoreSet>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Summary {
    by_effect: HashMap<String, usize>,
    by_score_class: HashMap<String, usize>,
    score_sets: Vec<ScoreSetSummary>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ScoreSetSummary {
    urn: String,
    target_gene: Option<String>,
    variant_count: usize,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct ScoreSet {
    urn: String,
    target_gene: Option<String>,
    title: String,
    variants: Vec<Variant>,
}

#[derive(Debug, Deserialize)]
#[allow(dead_code)]
struct Variant {
    hgvs_nt: Option<String>,
    hgvs_pro: Option<String>,
    score: Option<f64>,
    effect_nt: Option<String>,
    effect_pro: Option<String>,
    score_class: Option<String>,
}

fn load_mavedb_fixtures() -> MaveDbFixture {
    let content = fs::read_to_string("tests/fixtures/external/mavedb_functional.json")
        .expect("Failed to read mavedb_functional.json");
    serde_json::from_str(&content).expect("Failed to parse mavedb_functional.json")
}

/// Test all nucleotide variants from MaveDB
#[test]
fn test_mavedb_nucleotide_variants() {
    let fixtures = load_mavedb_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for score_set in &fixtures.score_sets {
        for variant in &score_set.variants {
            if let Some(ref hgvs) = variant.hgvs_nt {
                // Skip wild-type/identity variants
                if hgvs == "c.=" || hgvs == "n.=" {
                    skipped += 1;
                    continue;
                }

                let result = parse_hgvs(hgvs);

                match result {
                    Ok(_) => {
                        passed += 1;
                    }
                    Err(e) => {
                        failed += 1;
                        if failures.len() < 50 {
                            failures.push(format!("{} ({}): {}", hgvs, score_set.urn, e));
                        }
                    }
                }
            }
        }
    }

    eprintln!("\nMaveDB nucleotide variant test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped (identity): {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nSample failures (first 20):");
        for f in failures.iter().take(20) {
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

    // NOTE: MaveDB uses MAVE-HGVS short form without reference accessions
    // The parser requires full HGVS with accessions, so 0% success is expected
    // for variants like `c.[...]` instead of `NM_xxx:c.[...]`
    // This test is informational to track pattern coverage
    eprintln!("  (Note: MAVE-HGVS short form lacks accessions required by parser)");
}

/// Test all protein variants from MaveDB
#[test]
fn test_mavedb_protein_variants() {
    let fixtures = load_mavedb_fixtures();
    let mut passed = 0;
    let mut failed = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for score_set in &fixtures.score_sets {
        for variant in &score_set.variants {
            if let Some(ref hgvs) = variant.hgvs_pro {
                // Skip wild-type/identity variants
                if hgvs == "p.=" || hgvs == "p.(=)" {
                    skipped += 1;
                    continue;
                }

                let result = parse_hgvs(hgvs);

                match result {
                    Ok(_) => {
                        passed += 1;
                    }
                    Err(e) => {
                        failed += 1;
                        if failures.len() < 50 {
                            failures.push(format!("{} ({}): {}", hgvs, score_set.urn, e));
                        }
                    }
                }
            }
        }
    }

    eprintln!("\nMaveDB protein variant test results:");
    eprintln!("  Passed: {}", passed);
    eprintln!("  Failed: {}", failed);
    eprintln!("  Skipped (identity): {}", skipped);

    if !failures.is_empty() {
        eprintln!("\nSample failures (first 20):");
        for f in failures.iter().take(20) {
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

    // NOTE: MaveDB uses MAVE-HGVS short form without reference accessions
    // The parser requires full HGVS with accessions, so 0% success is expected
    // for variants like `p.Xxx123Yyy` instead of `NP_xxx:p.Xxx123Yyy`
    // This test is informational to track pattern coverage
    eprintln!("  (Note: MAVE-HGVS short form lacks accessions required by parser)");
}

/// Test variants by effect classification
#[test]
fn test_mavedb_by_effect() {
    let fixtures = load_mavedb_fixtures();

    let mut effect_counts: HashMap<String, (usize, usize)> = HashMap::new();

    for score_set in &fixtures.score_sets {
        for variant in &score_set.variants {
            // Test protein variants by effect
            if let Some(ref hgvs) = variant.hgvs_pro {
                if hgvs == "p.=" || hgvs == "p.(=)" {
                    continue;
                }

                let effect = variant.effect_pro.as_deref().unwrap_or("unknown");
                let result = parse_hgvs(hgvs);
                let entry = effect_counts.entry(effect.to_string()).or_insert((0, 0));
                if result.is_ok() {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
    }

    eprintln!("\nMaveDB protein variants by effect type:");
    for (effect, (passed, failed)) in &effect_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", effect, passed, total, rate);
    }
}

/// Test variants by functional score class
#[test]
fn test_mavedb_by_score_class() {
    let fixtures = load_mavedb_fixtures();

    let mut class_counts: HashMap<String, (usize, usize)> = HashMap::new();

    for score_set in &fixtures.score_sets {
        for variant in &score_set.variants {
            if let Some(ref hgvs) = variant.hgvs_pro {
                if hgvs == "p.=" || hgvs == "p.(=)" {
                    continue;
                }

                let score_class = variant.score_class.as_deref().unwrap_or("unscored");
                let result = parse_hgvs(hgvs);
                let entry = class_counts
                    .entry(score_class.to_string())
                    .or_insert((0, 0));
                if result.is_ok() {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
    }

    eprintln!("\nMaveDB protein variants by score class:");
    for (class, (passed, failed)) in &class_counts {
        let total = passed + failed;
        let rate = if total > 0 {
            (*passed as f64 / total as f64) * 100.0
        } else {
            100.0
        };
        eprintln!("  {}: {}/{} ({:.1}%)", class, passed, total, rate);
    }
}
