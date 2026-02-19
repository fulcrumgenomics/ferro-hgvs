//! Integration tests for local biocommons/hgvs setup.
//!
//! These tests require local UTA and SeqRepo infrastructure and are ignored by default.
//! Run with: `cargo test --features benchmark --test biocommons_local_tests -- --ignored`

#![cfg(feature = "benchmark")]

use std::path::PathBuf;

use ferro_hgvs::benchmark::{
    check_docker_available, check_seqrepo, check_uta_connection, check_uta_local,
    has_biocommons_normalizer, load_biocommons_settings, normalize_biocommons_single,
    write_biocommons_settings, BiocommonsLocalConfig, BiocommonsSettings,
};

/// Test that Docker is available on the system.
#[test]
#[ignore = "Requires Docker installed"]
fn test_docker_available() {
    let available = check_docker_available();
    assert!(available, "Docker should be available");
}

/// Test local UTA database connection.
#[test]
#[ignore = "Requires local UTA database running"]
fn test_local_uta_connection() {
    // Default UTA port
    let connected = check_uta_local(5432);
    assert!(connected, "Should connect to local UTA at port 5432");
}

/// Test UTA connection with explicit URL.
#[test]
#[ignore = "Requires local UTA database running"]
fn test_uta_connection_with_url() {
    let url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b";
    let connected = check_uta_connection(Some(url));
    assert!(connected, "Should connect to local UTA with explicit URL");
}

/// Test SeqRepo directory validation.
#[test]
#[ignore = "Requires local SeqRepo setup"]
fn test_local_seqrepo_access() {
    // This test requires SEQREPO_DIR environment variable or a known path
    let seqrepo_dir = std::env::var("SEQREPO_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/usr/local/share/seqrepo/2021-01-29"));

    let valid = check_seqrepo(&seqrepo_dir);
    assert!(
        valid,
        "SeqRepo directory should be valid: {:?}",
        seqrepo_dir
    );
}

/// Test settings file write and read round-trip.
#[test]
fn test_settings_round_trip() {
    let temp_dir = tempfile::tempdir().unwrap();
    let settings_path = temp_dir.path().join("biocommons_settings.conf");

    let uta_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b";
    let seqrepo_dir = PathBuf::from("/path/to/seqrepo/2021-01-29");

    // Write settings
    write_biocommons_settings(&settings_path, uta_url, &seqrepo_dir).unwrap();

    // Read settings back
    let settings = load_biocommons_settings(&settings_path).unwrap();

    assert_eq!(settings.uta_db_url, uta_url);
    assert_eq!(settings.seqrepo_dir, seqrepo_dir);
}

/// Test biocommons normalizer availability check.
#[test]
fn test_biocommons_normalizer_check() {
    let available = has_biocommons_normalizer();
    // This test passes whether or not biocommons is installed
    // It just verifies the check function works
    println!(
        "biocommons/hgvs normalizer available: {}",
        if available { "yes" } else { "no" }
    );
}

/// Test local biocommons normalization with a simple variant.
#[test]
#[ignore = "Requires local UTA + SeqRepo setup"]
fn test_local_biocommons_normalization() {
    // Get settings from environment or use defaults
    let uta_url = std::env::var("UTA_DB_URL").unwrap_or_else(|_| {
        "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b".to_string()
    });

    let seqrepo_dir = std::env::var("HGVS_SEQREPO_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/usr/local/share/seqrepo/2021-01-29"));

    // Verify setup is available
    assert!(check_uta_local(5432), "UTA must be running locally");
    assert!(check_seqrepo(&seqrepo_dir), "SeqRepo must be available");

    // Test normalization of a simple variant
    let input = "NM_000546.6:c.215C>G";
    let result = normalize_biocommons_single(input, Some(&uta_url), None, None);

    assert!(result.is_ok(), "Normalization should succeed: {:?}", result);
    let parse_result = result.unwrap();
    assert!(parse_result.success, "Normalization should succeed");
    assert_eq!(
        parse_result.output.as_deref(),
        Some("NM_000546.6:c.215C>G"),
        "Simple substitution should normalize to itself"
    );
}

/// Test local vs remote biocommons results are equivalent.
#[test]
#[ignore = "Requires local UTA + SeqRepo setup"]
fn test_local_vs_remote_equivalence() {
    let local_url = "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b";

    // Variants that should normalize identically
    let test_variants = vec![
        "NM_000546.6:c.215C>G",
        "NM_000059.4:c.68_69del",
        "NM_000546.6:c.742C>T",
    ];

    for variant in test_variants {
        // Normalize with local UTA
        let local_result = normalize_biocommons_single(variant, Some(local_url), None, None);

        // Normalize with remote UTA
        let remote_result = normalize_biocommons_single(variant, None, None, None);

        match (&local_result, &remote_result) {
            (Ok(local), Ok(remote)) => {
                // Compare outputs if both succeeded
                if local.success && remote.success {
                    assert_eq!(
                        local.output, remote.output,
                        "Local and remote should produce same result for {}",
                        variant
                    );
                } else if !local.success && !remote.success {
                    // Both failed to normalize - also equivalent
                } else {
                    panic!(
                        "Success mismatch for {}: local.success={}, remote.success={}",
                        variant, local.success, remote.success
                    );
                }
            }
            _ => {
                panic!(
                    "Result mismatch for {}: local={:?}, remote={:?}",
                    variant, local_result, remote_result
                );
            }
        }
    }
}

/// Test BiocommonsLocalConfig default values.
#[test]
fn test_biocommons_local_config_defaults() {
    let config = BiocommonsLocalConfig {
        uta_container_name: "ferro-uta".to_string(),
        uta_image_tag: "uta_20210129b".to_string(),
        uta_port: 5432,
        seqrepo_dir: PathBuf::from("/path/to/seqrepo"),
        seqrepo_instance: "2021-01-29".to_string(),
    };

    assert_eq!(config.uta_container_name, "ferro-uta");
    assert_eq!(config.uta_image_tag, "uta_20210129b");
    assert_eq!(config.uta_port, 5432);
    assert_eq!(config.seqrepo_instance, "2021-01-29");
}

/// Test BiocommonsSettings struct.
#[test]
fn test_biocommons_settings_struct() {
    let settings = BiocommonsSettings {
        uta_db_url: "postgresql://anonymous:anonymous@localhost:5432/uta/uta_20210129b".to_string(),
        seqrepo_dir: PathBuf::from("/usr/local/share/seqrepo/2021-01-29"),
    };

    assert!(settings.uta_db_url.contains("localhost"));
    assert!(settings.seqrepo_dir.to_string_lossy().contains("seqrepo"));
}

/// Test that settings file with invalid format returns an error.
#[test]
fn test_invalid_settings_format() {
    let temp_dir = tempfile::tempdir().unwrap();
    let settings_path = temp_dir.path().join("invalid_settings.conf");

    // Write invalid content
    std::fs::write(&settings_path, "this is not valid settings format").unwrap();

    let result = load_biocommons_settings(&settings_path);
    assert!(result.is_err(), "Should fail to parse invalid settings");
}

/// Test settings file with missing required fields.
#[test]
fn test_settings_missing_fields() {
    let temp_dir = tempfile::tempdir().unwrap();
    let settings_path = temp_dir.path().join("incomplete_settings.conf");

    // Write settings with only UTA_DB_URL
    std::fs::write(&settings_path, "UTA_DB_URL = postgresql://localhost/uta").unwrap();

    let result = load_biocommons_settings(&settings_path);
    assert!(
        result.is_err(),
        "Should fail when HGVS_SEQREPO_DIR is missing"
    );
}
