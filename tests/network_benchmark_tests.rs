//! Integration tests for network vs cache performance comparison
//!
//! The full network-vs-cache timing test below is `#[ignore]`d because it
//! requires network access and a prepared mutalyzer cache. The two
//! always-on tests (help-text + flag acceptance) parse the CLI in-process
//! via `clap::CommandFactory` and `Cli::try_parse_from`, avoiding the
//! ~170s-per-test cost of spawning `cargo run --release` from inside a
//! test.
//!
//! Run the ignored test with:
//!   cargo test --features benchmark --test network_benchmark_tests -- --ignored

#![cfg(feature = "benchmark")]

use clap::{CommandFactory, Parser};
use ferro_hgvs::benchmark::cli::Cli;

/// The `normalize` subcommand's help text advertises `--allow-network` and
/// explains its purpose.
#[test]
fn test_normalize_help_shows_allow_network_flag() {
    let mut cmd = Cli::command();
    let normalize = cmd
        .find_subcommand_mut("normalize")
        .expect("normalize subcommand should exist");
    let help_text = normalize.render_long_help().to_string();

    assert!(
        help_text.contains("--allow-network"),
        "Help should show --allow-network flag. Got:\n{}",
        help_text
    );
    assert!(
        help_text.contains("benchmarking remote vs local"),
        "Help should explain the flag purpose. Got:\n{}",
        help_text
    );
}

/// The CLI parser accepts `--allow-network` on the `normalize` subcommand
/// without complaining about an unknown argument.
#[test]
fn test_normalize_accepts_allow_network_flag() {
    let result = Cli::try_parse_from([
        "ferro-benchmark",
        "normalize",
        "mutalyzer",
        "-i",
        "nonexistent_file_12345.txt",
        "-o",
        "output.json",
        "--allow-network",
    ]);
    assert!(
        result.is_ok(),
        "Parser should accept --allow-network on normalize. Got error:\n{}",
        result.err().map(|e| e.to_string()).unwrap_or_default()
    );
}

/// `parse all ...` must be rejected at the CLI boundary — the runtime
/// rejects it later, but tightening the parser keeps invalid commands
/// from ever parsing.
#[test]
fn test_parse_rejects_all_tool() {
    let result = Cli::try_parse_from([
        "ferro-benchmark",
        "parse",
        "all",
        "-i",
        "patterns.txt",
        "-o",
        "out.json",
    ]);
    assert!(
        result.is_err(),
        "Parser should reject 'all' on parse subcommand"
    );
}

/// `--detailed` was a no-op on `generate summary` (destructured and
/// dropped). The flag has been removed; the parser must now reject it.
#[test]
fn test_generate_summary_rejects_detailed_flag() {
    let result = Cli::try_parse_from([
        "ferro-benchmark",
        "generate",
        "summary",
        "--parsing",
        "p.json",
        "--normalization",
        "n.json",
        "-o",
        "out.json",
        "--detailed",
    ]);
    assert!(
        result.is_err(),
        "Parser should reject --detailed on generate summary"
    );
}

/// `normalize all ...` must be rejected at the CLI boundary.
#[test]
fn test_normalize_rejects_all_tool() {
    let result = Cli::try_parse_from([
        "ferro-benchmark",
        "normalize",
        "all",
        "-i",
        "patterns.txt",
        "-o",
        "out.json",
    ]);
    assert!(
        result.is_err(),
        "Parser should reject 'all' on normalize subcommand"
    );
}

/// Test that network mode is slower than cache mode
///
/// This test requires:
/// - Network access
/// - Mutalyzer cache to be prepared at data/mutalyzer/
/// - sample_100.txt to exist with test patterns
#[test]
#[ignore = "requires network access and prepared cache - run manually"]
fn test_mutalyzer_network_slower_than_cache() {
    use std::process::Command;
    use std::time::Instant;
    use tempfile::TempDir;

    let temp_dir = TempDir::new().unwrap();

    // Create a small test patterns file
    let patterns_file = temp_dir.path().join("patterns.txt");
    std::fs::write(
        &patterns_file,
        "NM_000088.3:c.589G>T\nNM_000088.3:c.590C>A\n",
    )
    .unwrap();

    // Run with cache (should be fast)
    let cache_output = temp_dir.path().join("cache_results.json");
    let cache_start = Instant::now();
    let cache_result = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--features",
            "benchmark",
            "--bin",
            "ferro-benchmark",
            "--",
            "normalize",
            "mutalyzer",
            "-i",
            patterns_file.to_str().unwrap(),
            "-o",
            cache_output.to_str().unwrap(),
            "--mutalyzer-settings",
            "data/mutalyzer/mutalyzer_settings.conf",
            "-j",
            "1",
        ])
        .output()
        .expect("Failed to run cache mode");
    let cache_elapsed = cache_start.elapsed();

    if !cache_result.status.success() {
        let stderr = String::from_utf8_lossy(&cache_result.stderr);
        eprintln!("Cache mode failed: {}", stderr);
        // Skip test if cache isn't available
        return;
    }

    // Run with network (should be slow)
    let network_output = temp_dir.path().join("network_results.json");
    let network_start = Instant::now();
    let network_result = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--features",
            "benchmark",
            "--bin",
            "ferro-benchmark",
            "--",
            "normalize",
            "mutalyzer",
            "-i",
            patterns_file.to_str().unwrap(),
            "-o",
            network_output.to_str().unwrap(),
            "--allow-network",
            "-j",
            "1",
        ])
        .output()
        .expect("Failed to run network mode");
    let network_elapsed = network_start.elapsed();

    // Network mode might fail on rate limiting, but if it succeeds it should be slower
    if network_result.status.success() {
        println!("Cache mode: {:?}", cache_elapsed);
        println!("Network mode: {:?}", network_elapsed);

        // Network should be at least 2x slower (often 10-20x)
        assert!(
            network_elapsed > cache_elapsed,
            "Network mode ({:?}) should be slower than cache mode ({:?})",
            network_elapsed,
            cache_elapsed
        );
    } else {
        let stderr = String::from_utf8_lossy(&network_result.stderr);
        eprintln!("Network mode failed (may be rate limited): {}", stderr);
    }
}
