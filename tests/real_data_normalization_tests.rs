//! Tests using real reference data from benchmark-output
//!
//! These tests require the benchmark-output directory to be present.

use ferro_hgvs::{parse_hgvs, MultiFastaProvider, Normalizer, ReferenceProvider};
use std::path::Path;

fn get_provider() -> Option<MultiFastaProvider> {
    let manifest_path = Path::new("benchmark-output/manifest.json");
    if manifest_path.exists() {
        MultiFastaProvider::from_manifest(manifest_path).ok()
    } else {
        None
    }
}

#[test]
fn test_nm_001408491_c517dela_should_shift() {
    // This test reproduces the issue found in mutalyzer comparison:
    // NM_001408491.1:c.517delA should normalize to c.518del (3' rule)
    // because c.517 and c.518 are both 'A'
    //
    // Mutalyzer output: NM_001408491.1:c.518del
    // Ferro output (before fix): NM_001408491.1:c.517del

    let Some(provider) = get_provider() else {
        eprintln!("Skipping test: benchmark-output not available");
        return;
    };

    // First, let's get the transcript and examine the sequence
    let transcript = match provider.get_transcript("NM_001408491.1") {
        Ok(tx) => tx,
        Err(e) => {
            eprintln!("Skipping test: transcript not found: {}", e);
            return;
        }
    };

    println!("=== Transcript NM_001408491.1 ===");
    println!("Sequence length: {}", transcript.sequence.len());
    println!("CDS start (1-based): {:?}", transcript.cds_start);
    println!("CDS end: {:?}", transcript.cds_end);
    println!("Exon count: {}", transcript.exons.len());
    for (i, exon) in transcript.exons.iter().enumerate() {
        println!(
            "  Exon {}: {} - {} (number: {})",
            i + 1,
            exon.start,
            exon.end,
            exon.number
        );
    }

    // Calculate transcript positions for c.517 and c.518
    let cds_start = transcript.cds_start.unwrap_or(1);
    let tx_pos_517 = cds_start + 517 - 1; // c.517 -> tx position
    let tx_pos_518 = cds_start + 518 - 1; // c.518 -> tx position

    println!("\nCoordinate mapping:");
    println!("  CDS start: {}", cds_start);
    println!(
        "  c.517 -> tx position {} (0-based: {})",
        tx_pos_517,
        tx_pos_517 - 1
    );
    println!(
        "  c.518 -> tx position {} (0-based: {})",
        tx_pos_518,
        tx_pos_518 - 1
    );

    // Get the bases at these positions
    let seq = transcript.sequence.as_bytes();
    let idx_517 = (tx_pos_517 - 1) as usize; // 0-based
    let idx_518 = (tx_pos_518 - 1) as usize; // 0-based

    if idx_517 < seq.len() && idx_518 < seq.len() {
        println!("\nSequence around c.517-518:");
        // Show context
        let start = idx_517.saturating_sub(3);
        let end = (idx_518 + 4).min(seq.len());
        let context: String = seq[start..end].iter().map(|&b| b as char).collect();
        println!("  ... {} ...", context);
        println!(
            "  c.517 (tx pos {}, 0-based {}): {}",
            tx_pos_517, idx_517, seq[idx_517] as char
        );
        println!(
            "  c.518 (tx pos {}, 0-based {}): {}",
            tx_pos_518, idx_518, seq[idx_518] as char
        );

        // The key assertion: both should be 'A' for the 3' shift to apply
        if seq[idx_517] == b'A' && seq[idx_518] == b'A' {
            println!("\n  CONFIRMED: Both c.517 and c.518 are 'A'");
            println!("  According to 3' rule, c.517delA should shift to c.518del");
        }
    } else {
        eprintln!("Position out of range!");
    }

    // Find which exon contains c.517
    println!("\nExon containing c.517 (tx pos {}):", tx_pos_517);
    for exon in &transcript.exons {
        if tx_pos_517 >= exon.start && tx_pos_517 <= exon.end {
            println!("  Exon {} (tx {}-{})", exon.number, exon.start, exon.end);
            // Check if c.518 is also in this exon
            if tx_pos_518 <= exon.end {
                println!("  c.518 is also in this exon - no boundary issues");
            } else {
                println!("  WARNING: c.518 is NOT in this exon - boundary may prevent shift!");
            }
        }
    }

    // Now test normalization
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_001408491.1:c.517delA").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    let output = format!("{}", result);

    println!("\n=== Normalization Result ===");
    println!("Input:  NM_001408491.1:c.517delA");
    println!("Output: {}", output);

    // The test expectation
    assert!(
        output.contains("c.518del"),
        "Expected c.517delA to normalize to c.518del (3' rule), got: {}",
        output
    );
}

#[test]
fn test_potential_bug_deletion_shift_nm033517() {
    // Investigate: ferro outputs c.1324del, mutalyzer outputs c.1326del
    // One of these is wrong about 3' shifting
    let Some(provider) = get_provider() else {
        eprintln!("Skipping test: benchmark-output not available");
        return;
    };

    let transcript = match provider.get_transcript("NM_033517.1") {
        Ok(tx) => tx,
        Err(e) => {
            eprintln!("Skipping test: transcript not found: {}", e);
            return;
        }
    };

    let cds_start = transcript.cds_start.unwrap_or(1);
    let seq = transcript.sequence.as_bytes();

    // Check c.1324, c.1325, c.1326
    for pos in 1324..=1328 {
        let tx_pos = cds_start + pos - 1;
        let idx = (tx_pos - 1) as usize;
        if idx < seq.len() {
            println!(
                "c.{} (tx {}, idx {}): {}",
                pos, tx_pos, idx, seq[idx] as char
            );
        }
    }

    // The 3' rule: if there are consecutive identical bases, shift to rightmost
    // If c.1324-1326 are all the same base, the deletion should be at c.1326
    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_033517.1:c.1324del").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    println!("\nInput:  NM_033517.1:c.1324del");
    println!("Output: {}", result);
}

#[test]
fn test_potential_bug_delins_shift() {
    // Investigate: ferro shifts delins positions
    // Input:     c.2139_2140delinsTATGCA
    // Ferro:     c.2140_2141delinsTATGCA
    // Mutalyzer: c.2139_2140delinsTATGCA
    //
    // HGVS spec: delins should NOT be 3' shifted like del/dup
    let Some(provider) = get_provider() else {
        eprintln!("Skipping test: benchmark-output not available");
        return;
    };

    let transcript = match provider.get_transcript("NM_001282424.3") {
        Ok(tx) => tx,
        Err(e) => {
            eprintln!("Skipping test: transcript not found: {}", e);
            return;
        }
    };

    let cds_start = transcript.cds_start.unwrap_or(1);
    let seq = transcript.sequence.as_bytes();

    // Check sequence around c.2139-2141
    println!("Sequence around c.2139-2141:");
    for pos in 2137..=2145 {
        let tx_pos = cds_start + pos - 1;
        let idx = (tx_pos - 1) as usize;
        if idx < seq.len() {
            println!(
                "c.{} (tx {}, idx {}): {}",
                pos, tx_pos, idx, seq[idx] as char
            );
        }
    }

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_001282424.3:c.2139_2140delinsTATGCA").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    println!("\nInput:  NM_001282424.3:c.2139_2140delinsTATGCA");
    println!("Output: {}", result);

    // Check if ferro is incorrectly shifting delins
    let output = format!("{}", result);
    if output.contains("2140_2141") {
        println!("\n⚠️  POTENTIAL BUG: Ferro shifted delins from 2139_2140 to 2140_2141");
        println!("    HGVS spec says delins should use original position, not shift");
    }
}

#[test]
fn test_5utr_duplication_shifting() {
    // Investigate: ferro shifts 5'UTR duplications, mutalyzer doesn't
    // Input:     c.-56_-47dup
    // Ferro:     c.-55_-46dup
    // Mutalyzer: c.-56_-47dup
    //
    // Question: Does 3' rule apply to 5'UTR (negative coordinates)?
    // In 5'UTR, the 3' direction is towards the start codon (+1)
    // So -55 is more 3' than -56
    let Some(provider) = get_provider() else {
        eprintln!("Skipping test: benchmark-output not available");
        return;
    };

    let transcript = match provider.get_transcript("NM_001394148.2") {
        Ok(tx) => tx,
        Err(e) => {
            eprintln!("Skipping test: transcript not found: {}", e);
            return;
        }
    };

    let cds_start = transcript.cds_start.unwrap_or(1);
    let seq = transcript.sequence.as_bytes();

    // Check sequence around c.-56 to c.-45
    println!("CDS start: {}", cds_start);
    println!("Sequence around c.-56 to c.-45:");
    for offset in -60i64..=-40 {
        // c.-N means position cds_start - N in 1-based
        let tx_pos = (cds_start as i64 + offset) as u64;
        if tx_pos >= 1 {
            let idx = (tx_pos - 1) as usize;
            if idx < seq.len() {
                println!(
                    "c.{} (tx {}, idx {}): {}",
                    offset, tx_pos, idx, seq[idx] as char
                );
            }
        }
    }

    let normalizer = Normalizer::new(provider);
    let variant = parse_hgvs("NM_001394148.2:c.-56_-47dup").unwrap();
    let result = normalizer.normalize(&variant).unwrap();
    println!("\nInput:  NM_001394148.2:c.-56_-47dup");
    println!("Output: {}", result);

    // Note: The correct behavior depends on HGVS spec interpretation
    // If 3' rule applies to 5'UTR, ferro is correct (shift towards +1)
    // If 3' rule doesn't apply to 5'UTR, mutalyzer is correct
}

#[test]
fn test_compare_with_mutalyzer_sequence() {
    // Compare what we have in cdot vs what mutalyzer reports
    let Some(provider) = get_provider() else {
        eprintln!("Skipping test: benchmark-output not available");
        return;
    };

    let transcript = match provider.get_transcript("NM_001408491.1") {
        Ok(tx) => tx,
        Err(e) => {
            eprintln!("Skipping test: transcript not found: {}", e);
            return;
        }
    };

    // From mutalyzer normalizer output for NM_001408491.1:c.517delA:
    // - Output: NM_001408491.1:c.518del
    // - This means the deleted base shifts from c.517 to c.518
    // - This is only correct if c.517 == c.518 (both are the same base)

    let cds_start = transcript.cds_start.unwrap_or(1);
    let tx_pos_517 = cds_start + 517 - 1;
    let tx_pos_518 = cds_start + 518 - 1;

    let seq = transcript.sequence.as_bytes();
    let idx_517 = (tx_pos_517 - 1) as usize;
    let idx_518 = (tx_pos_518 - 1) as usize;

    if idx_517 < seq.len() && idx_518 < seq.len() {
        let base_517 = seq[idx_517] as char;
        let base_518 = seq[idx_518] as char;

        println!("cdot sequence at c.517: {}", base_517);
        println!("cdot sequence at c.518: {}", base_518);

        // For the mutalyzer result to be correct, both must be 'A'
        assert_eq!(
            base_517, 'A',
            "Expected c.517 to be 'A' per mutalyzer result"
        );
        assert_eq!(
            base_518, 'A',
            "Expected c.518 to be 'A' for 3' shift to be valid"
        );
    }
}
