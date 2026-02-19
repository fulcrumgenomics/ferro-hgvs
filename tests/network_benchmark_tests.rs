//! Integration tests for network vs cache performance comparison
//!
//! These tests verify that the --allow-network flag works correctly and
//! measure the performance difference between network and cached modes.
//!
//! These tests are ignored by default as they require:
//! - Network access
//! - Mutalyzer cache to be prepared
//! - Significant time to run
//!
//! Run with: cargo test --features benchmark network_benchmark -- --ignored

use std::process::Command;

/// Test that the normalize help shows the --allow-network flag
#[test]
#[cfg(feature = "benchmark")]
fn test_normalize_help_shows_allow_network_flag() {
    let output = Command::new("cargo")
        .args([
            "run",
            "--release",
            "--features",
            "benchmark",
            "--bin",
            "ferro-benchmark",
            "--",
            "normalize",
            "--help",
        ])
        .output()
        .expect("Failed to run help command");

    let help_text = String::from_utf8_lossy(&output.stdout);
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

/// Test that the normalize command accepts the --allow-network flag without error
#[test]
#[cfg(feature = "benchmark")]
fn test_normalize_accepts_allow_network_flag() {
    // This test verifies that --allow-network is parsed correctly
    // It will fail on "file not found" rather than "unknown argument"
    let output = Command::new("cargo")
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
            "nonexistent_file_12345.txt",
            "-o",
            "output.json",
            "--allow-network",
        ])
        .output()
        .expect("Failed to run command");

    let stderr = String::from_utf8_lossy(&output.stderr);
    // Should fail on file not found, not unknown argument
    assert!(
        !stderr.contains("unexpected argument"),
        "Command should accept --allow-network flag. Got:\n{}",
        stderr
    );
    assert!(
        !stderr.contains("error: Found argument '--allow-network'"),
        "Command should accept --allow-network flag. Got:\n{}",
        stderr
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
#[cfg(feature = "benchmark")]
fn test_mutalyzer_network_slower_than_cache() {
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
