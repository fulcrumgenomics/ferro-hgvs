//! gnomAD HGVS validation test
//!
//! Run with: cargo test --release test_gnomad_validation -- --nocapture --ignored

use ferro_hgvs::parse_hgvs;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[test]
#[ignore] // Run manually - requires extracted data
fn test_gnomad_validation() {
    let file = File::open("data/gnomad/gnomad_chr21_hgvs.txt")
        .expect("Extract gnomAD HGVS first: python scripts/extract_gnomad_hgvs.py ...");
    let reader = BufReader::new(file);

    let mut total = 0;
    let mut passed = 0;
    let mut failures: HashMap<String, Vec<String>> = HashMap::new();

    for line in reader.lines() {
        let pattern = line.unwrap();
        if pattern.is_empty() {
            continue;
        }

        total += 1;
        match parse_hgvs(&pattern) {
            Ok(_) => passed += 1,
            Err(e) => {
                let key = format!("{:?}", e).chars().take(80).collect::<String>();
                let entry = failures.entry(key).or_default();
                if entry.len() < 5 {
                    entry.push(pattern);
                }
            }
        }
    }

    let rate = (passed as f64 / total as f64) * 100.0;

    println!("\n=== gnomAD chr21 Validation Results ===");
    println!("Total HGVS expressions: {}", total);
    println!("Passed: {}", passed);
    println!("Failed: {}", total - passed);
    println!("Coverage: {:.4}%", rate);

    println!("\nTop failure categories:");
    let mut sorted_failures: Vec<_> = failures.iter().collect();
    sorted_failures.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    for (error, examples) in sorted_failures.iter().take(15) {
        println!("  [{}x] {}", examples.len(), error);
        for ex in examples.iter().take(2) {
            println!("       - {}", ex);
        }
    }

    assert!(rate > 99.0, "gnomAD coverage below 99%: {:.4}%", rate);
}
