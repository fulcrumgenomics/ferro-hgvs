//! ExAC (Exome Aggregation Consortium) validation
//!
//! Tests parsing of VEP-generated HGVS expressions from ExAC GRCh38.

use ferro_hgvs::parse_hgvs;
use std::collections::HashMap;
use std::fs;

/// Test parsing of ExAC HGVS patterns
#[test]
#[ignore = "requires data/exac/exac_grch38_hgvs.txt"]
fn test_exac_validation() {
    let path = "data/exac/exac_grch38_hgvs.txt";

    let content = match fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Could not read {}: {}", path, e);
            eprintln!("Download and extract ExAC data first:");
            eprintln!("  curl -L http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/ExAC.0.3.GRCh38.vcf.gz -o data/exac/ExAC.0.3.GRCh38.vcf.gz");
            eprintln!("  python3 scripts/extract_exac_hgvs.py data/exac/ExAC.0.3.GRCh38.vcf.gz data/exac/exac_grch38_hgvs.txt");
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

    println!("\n=== ExAC GRCh38 Validation Results ===");
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

    // ExAC has VEP annotations - expect very high coverage
    assert!(
        coverage >= 99.0,
        "ExAC coverage {:.2}% below 99% threshold",
        coverage
    );
}
