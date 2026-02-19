//! CIViC (Clinical Interpretation of Variants in Cancer) validation
//!
//! Tests parsing of curated HGVS expressions from the CIViC database.

use ferro_hgvs::parse_hgvs;
use std::collections::HashMap;
use std::fs;

/// Test parsing of CIViC HGVS patterns
#[test]
#[ignore = "requires data/civic/civic_hgvs.txt"]
fn test_civic_validation() {
    let path = "data/civic/civic_hgvs.txt";

    let content = match fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Could not read {}: {}", path, e);
            eprintln!("Download CIViC data first:");
            eprintln!(
                "  curl -sL 'https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv' | \\"
            );
            eprintln!(
                "  python3 -c \"import csv,sys; [print(p.strip()) for row in csv.DictReader(sys.stdin, delimiter='\\t') for p in row.get('hgvs_descriptions','').split(',') if ':' in p]\" > {}", path
            );
            return;
        }
    };

    let mut passed = 0;
    let mut failed = 0;
    let mut failure_categories: HashMap<String, Vec<String>> = HashMap::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        match parse_hgvs(line) {
            Ok(_) => passed += 1,
            Err(e) => {
                failed += 1;
                let category = format!("{:?}", e).chars().take(80).collect::<String>();
                failure_categories
                    .entry(category)
                    .or_default()
                    .push(line.to_string());
            }
        }
    }

    let total = passed + failed;
    let coverage = if total > 0 {
        (passed as f64 / total as f64) * 100.0
    } else {
        0.0
    };

    println!("\n=== CIViC Validation Results ===");
    println!("Total HGVS expressions: {}", total);
    println!("Passed: {}", passed);
    println!("Failed: {}", failed);
    println!("Coverage: {:.4}%", coverage);

    if !failure_categories.is_empty() {
        println!("\nTop failure categories:");
        let mut sorted: Vec<_> = failure_categories.iter().collect();
        sorted.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

        for (category, examples) in sorted.iter().take(15) {
            println!("  [{}x] {}", examples.len(), category);
            for ex in examples.iter().take(3) {
                println!("       - {}", ex);
            }
        }
    }

    // CIViC has curated clinical data - expect very high coverage
    assert!(
        coverage >= 95.0,
        "CIViC coverage {:.2}% below 95% threshold",
        coverage
    );
}
