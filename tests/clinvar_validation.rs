//! ClinVar HGVS validation test - ALL 42M expressions
//!
//! This test validates against ALL ClinVar HGVS expressions (~42M).
//! Takes ~23 seconds in release mode, ~7 minutes in debug mode.
//!
//! Run with: cargo test --release --features slow-tests test_clinvar_validation -- --nocapture

#![cfg(feature = "slow-tests")]

use ferro_hgvs::parse_hgvs;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[test]
fn test_clinvar_validation() {
    let file = File::open("tests/fixtures/bulk/hgvs4variation.txt.gz")
        .expect("ClinVar data not found. Run: curl -o tests/fixtures/bulk/hgvs4variation.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/hgvs4variation.txt.gz");
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);

    let mut total = 0;
    let mut passed = 0;
    let mut failures: HashMap<String, Vec<String>> = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        if i < 16 {
            continue;
        } // Skip 16 header lines

        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() >= 9 {
            // Column 6: NucleotideExpression
            if let Some(nuc) = fields.get(6) {
                if !nuc.is_empty() && *nuc != "-" {
                    total += 1;
                    match parse_hgvs(nuc) {
                        Ok(_) => passed += 1,
                        Err(e) => {
                            let key = format!("{:?}", e).chars().take(80).collect::<String>();
                            let entry = failures.entry(key).or_default();
                            if entry.len() < 5 {
                                entry.push(nuc.to_string());
                            }
                        }
                    }
                }
            }

            // Column 8: ProteinExpression
            if let Some(prot) = fields.get(8) {
                if !prot.is_empty() && *prot != "-" {
                    total += 1;
                    match parse_hgvs(prot) {
                        Ok(_) => passed += 1,
                        Err(e) => {
                            let key = format!("{:?}", e).chars().take(80).collect::<String>();
                            let entry = failures.entry(key).or_default();
                            if entry.len() < 5 {
                                entry.push(prot.to_string());
                            }
                        }
                    }
                }
            }
        }
    }

    let rate = (passed as f64 / total as f64) * 100.0;

    println!("\n=== ClinVar Validation Results ===");
    println!("Total HGVS expressions: {}", total);
    println!("Passed: {}", passed);
    println!("Failed: {}", total - passed);
    println!("Coverage: {:.4}%", rate);

    println!("\nTop failure categories:");
    let mut sorted_failures: Vec<_> = failures.iter().collect();
    sorted_failures.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    for (error, examples) in sorted_failures.iter().take(20) {
        println!("  [{}x] {}", examples.len(), error);
        for ex in examples.iter().take(3) {
            println!("       - {}", ex);
        }
    }

    assert!(rate > 99.0, "ClinVar coverage below 99%: {:.4}%", rate);
}
