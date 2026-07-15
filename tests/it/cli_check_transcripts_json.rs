//! Integration tests for `ferro check` on a standalone transcripts.json file
//! (#1012 comment 2): it must check the file directly instead of reporting
//! "Run 'ferro prepare' first".

use std::io::Write;
use std::process::Command;

fn write_json(name: &str, body: &str) -> tempfile::NamedTempFile {
    let mut tf = tempfile::Builder::new().suffix(name).tempfile().unwrap();
    tf.write_all(body.as_bytes()).unwrap();
    tf.flush().unwrap();
    tf
}

#[test]
fn check_accepts_a_standalone_transcripts_json() {
    let json = r#"{
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [
            {"id": "NM_1.1", "strand": "+", "sequence": "ACGTACGT",
             "exons": [{"number": 1, "start": 1, "end": 8}]}
        ]
    }"#;
    let file = write_json(".transcripts.json", json);

    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args(["check", "--reference", file.path().to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "check on a transcripts.json should succeed: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("transcripts.json Check") && stderr.contains("Transcripts: 1"),
        "expected a transcripts.json summary, got: {stderr}"
    );
    assert!(
        !stderr.contains("Run 'ferro prepare'"),
        "must not tell the user to run ferro prepare for a valid transcripts.json"
    );
}

#[test]
fn check_fails_on_incompatible_transcripts_json_version() {
    let file = write_json(
        ".transcripts.json",
        r#"{"version": "2.0", "transcripts": []}"#,
    );

    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args(["check", "--reference", file.path().to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "check must fail on an incompatible schema version"
    );
    // A nonzero exit alone also passes for unrelated failures; confirm it's the
    // schema-version diagnostic (not a parse/arg error) that reached the user.
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("newer than this build supports"),
        "expected the incompatible-version diagnostic, got: {stderr}"
    );
}

#[test]
fn check_nonexistent_path_reports_not_found_not_prepare() {
    // A typo'd path must say "does not exist", not tell the user to run prepare.
    // Build a missing path under a real tempdir (no committed absolute paths).
    let dir = tempfile::tempdir().unwrap();
    let missing = dir.path().join("no-such-transcripts.json");
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args(["check", "--reference"])
        .arg(&missing)
        .output()
        .unwrap();
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("does not exist"),
        "expected a 'does not exist' error, got: {stderr}"
    );
    assert!(
        !stderr.contains("Run 'ferro prepare'"),
        "must not suggest ferro prepare for a nonexistent path: {stderr}"
    );
}

#[test]
fn check_empty_transcripts_json_fails() {
    let file = write_json(".transcripts.json", r#"{"transcripts": []}"#);
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args(["check", "--reference", file.path().to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "an empty reference must fail the check, not be green-lit"
    );
    // #1014's `from_json` rejects a wholly-empty reference at load with "no usable
    // reference data", which preempts (and is equivalent to) `check`'s own
    // "no usable data" message — either way an empty reference is rejected.
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("no usable reference data"),
        "expected a 'no usable reference data' failure, got: {stderr}"
    );
}

#[test]
fn check_manifest_json_file_routes_to_directory() {
    // Pointing at a `manifest.json` file (not the directory) must be treated as the
    // prepared-reference directory case, not shoved through the transcripts.json
    // loader. A malformed manifest.json therefore yields the manifest-path error,
    // not a transcripts.json "unknown field" serde error.
    let dir = tempfile::tempdir().unwrap();
    let manifest = dir.path().join("manifest.json");
    std::fs::write(&manifest, b"not valid json").unwrap();

    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args(["check", "--reference", manifest.to_str().unwrap()])
        .output()
        .unwrap();
    // The malformed manifest must FAIL (a bare non-zero-status check would also
    // pass on an accidental success), and it must fail through the manifest/
    // reference path — not the transcripts.json check.
    assert!(
        !out.status.success(),
        "a malformed manifest.json must make the check fail"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("Reference data check failed"),
        "a manifest.json path must fail via the reference-data check: {stderr}"
    );
    assert!(
        !stderr.contains("transcripts.json Check"),
        "a manifest.json path must route to the directory path, not the transcripts.json check: {stderr}"
    );
}
