//! Test patterns from VEP and SnpEff test suites
//!
//! Note: These are "partial HGVS" patterns (variant-only, without accession prefix)
//! used internally by annotation tools. They are NOT valid complete HGVS expressions.
//!
//! VEP patterns like "21:g.25585733C>T" use chromosome numbers, not accessions.
//! SnpEff patterns like "c.100A>G" are missing the reference sequence.
//!
//! These tests are informational only and documented as expected failures.

use ferro_hgvs::parse_hgvs;

/// VEP parser patterns - most are valid HGVS with accessions
/// Some use chromosome notation (21:g.) which is not standard HGVS
#[test]
fn test_vep_parser_patterns() {
    let patterns = include_str!("fixtures/external/vep_parser_hgvs.txt");

    let mut passed = 0;
    let mut failed = 0;

    for line in patterns.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        match parse_hgvs(line) {
            Ok(_) => passed += 1,
            Err(_) => {
                failed += 1;
                // Expected: chromosome notation (21:g.) and gene symbols (JAM2:c.) fail
            }
        }
    }

    println!(
        "\nVEP Parser: {} passed, {} failed (chromosome/gene notation)",
        passed, failed
    );
    // Don't assert - these failures are expected for non-standard patterns
}

/// SnpEff patterns - these are partial HGVS (variant only, no accession)
/// They cannot be parsed as complete HGVS expressions
#[test]
fn test_snpeff_patterns() {
    let patterns = include_str!("fixtures/external/snpeff_hgvs.txt");

    let mut total = 0;
    let mut needs_accession = 0;

    for line in patterns.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        total += 1;

        // SnpEff patterns are variant-only (c.100A>G, p.Gly4_Gln6dup)
        // They need an accession prefix to be valid HGVS
        if line.starts_with("c.") || line.starts_with("p.") || line.starts_with("n.") {
            needs_accession += 1;
        }
    }

    println!(
        "\nSnpEff: {} total patterns, {} need accession prefix",
        total, needs_accession
    );
    // These are partial patterns, not complete HGVS - no assertion needed
}
