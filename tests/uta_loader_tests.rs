//! Integration tests for UTA loader.
//!
//! These tests require:
//! - Docker running
//! - ferro-uta container available
//! - Prepared ferro reference with cdot
//!
//! Run with: `cargo test --features benchmark --test uta_loader_tests -- --ignored`

#![cfg(feature = "benchmark")]

use std::path::PathBuf;
use std::process::Command;

use ferro_hgvs::benchmark::{load_cdot_alignments_to_uta, UtaLoadConfig, UtaLoadResult};

/// Test that UTA load config defaults are correct.
#[test]
fn test_uta_load_config_default() {
    let config = UtaLoadConfig::default();
    assert_eq!(config.uta_schema, "uta_20210129b");
    assert_eq!(config.container_name, "ferro-uta");
    assert_eq!(config.origin_name, "cdot-ferro");
    assert_eq!(config.batch_size, 10000);
}

/// Test that UtaLoadResult default is correct.
#[test]
fn test_uta_load_result_default() {
    let result = UtaLoadResult::default();
    assert_eq!(result.transcripts_loaded, 0);
    assert_eq!(result.transcripts_skipped, 0);
    assert_eq!(result.exon_sets_created, 0);
    assert_eq!(result.exons_created, 0);
    assert_eq!(result.exon_alns_created, 0);
    assert!(result.errors.is_empty());
}

/// Test loading cdot alignments into UTA.
#[test]
#[ignore = "Requires Docker UTA container and ferro reference"]
fn test_load_cdot_alignments() {
    let cdot_path = PathBuf::from("data/ferro/cdot/cdot-0.2.32.refseq.GRCh38.json");
    if !cdot_path.exists() {
        eprintln!("Skipping test: cdot file not found at {:?}", cdot_path);
        return;
    }

    let config = UtaLoadConfig::default();
    let result = load_cdot_alignments_to_uta(&cdot_path, &config);

    assert!(result.is_ok(), "Loading should succeed: {:?}", result.err());
    let result = result.unwrap();

    // Should have loaded or skipped transcripts
    assert!(
        result.transcripts_loaded + result.transcripts_skipped > 0,
        "Should have processed at least some transcripts"
    );
    println!(
        "Loaded: {}, Skipped: {}",
        result.transcripts_loaded, result.transcripts_skipped
    );
}

/// Test CLI --no-load-alignments flag is recognized.
#[test]
fn test_no_load_alignments_flag_help() {
    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--features",
            "benchmark",
            "--bin",
            "ferro-benchmark",
            "--",
            "prepare",
            "--help",
        ])
        .output()
        .expect("Failed to run help");

    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains("--no-load-alignments"),
        "Help should show --no-load-alignments flag"
    );
}

/// Verify a loaded transcript exists in UTA.
#[test]
#[ignore = "Requires Docker UTA container with loaded alignments"]
fn test_verify_loaded_transcript() {
    // Check that a newer transcript exists after loading
    let output = Command::new("docker")
        .args([
            "exec",
            "ferro-uta",
            "psql",
            "-U",
            "anonymous",
            "-d",
            "uta",
            "-t",
            "-c",
            "SELECT COUNT(*) FROM uta_20210129b.exon_set WHERE tx_ac = 'NM_001407958.1' AND alt_aln_method = 'cdot'",
        ])
        .output()
        .expect("Failed to query UTA");

    if !output.status.success() {
        eprintln!(
            "Docker query failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let count: i32 = String::from_utf8_lossy(&output.stdout)
        .trim()
        .parse()
        .unwrap_or(0);

    assert!(
        count > 0,
        "Transcript NM_001407958.1 should have exon_set entries after loading"
    );
}

/// Test that the cdot-ferro origin exists in UTA after loading.
#[test]
#[ignore = "Requires Docker UTA container with loaded alignments"]
fn test_cdot_ferro_origin_exists() {
    let output = Command::new("docker")
        .args([
            "exec",
            "ferro-uta",
            "psql",
            "-U",
            "anonymous",
            "-d",
            "uta",
            "-t",
            "-c",
            "SELECT origin_id FROM uta_20210129b.origin WHERE name = 'cdot-ferro'",
        ])
        .output()
        .expect("Failed to query UTA");

    if !output.status.success() {
        eprintln!(
            "Docker query failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return;
    }

    let output_str = String::from_utf8_lossy(&output.stdout);
    let origin_id: Option<i32> = output_str.trim().parse().ok();

    assert!(origin_id.is_some(), "cdot-ferro origin should exist in UTA");
}
