//! Integration test for `ferro convert-gff` flags added in Phase 3.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

fn write_gff3(content: &str) -> NamedTempFile {
    let mut tf = tempfile::Builder::new().suffix(".gff3").tempfile().unwrap();
    tf.write_all(content.as_bytes()).unwrap();
    tf.flush().unwrap();
    tf
}

#[test]
fn strict_mode_fails_on_malformed_coordinate() {
    let f = write_gff3(
        "##gff-version 3\n\
         chr1\t.\tgene\tnot_a_number\t500\t.\t+\t.\tID=g1\n",
    );
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            f.path().to_str().unwrap(),
            "--strict",
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "strict mode should fail on malformed record: stdout={}, stderr={}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn diagnostics_json_writes_file() {
    let f = write_gff3(
        "##gff-version 3\n\
         chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1\n",
    );
    let json_path = tempfile::Builder::new()
        .suffix(".json")
        .tempfile()
        .unwrap()
        .into_temp_path();
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "convert-gff",
            "--gff",
            f.path().to_str().unwrap(),
            "--diagnostics-json",
            json_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "convert-gff failed: stderr={}",
        String::from_utf8_lossy(&out.stderr)
    );
    let body = std::fs::read_to_string(&json_path).expect("diagnostics file");
    let diags: Vec<serde_json::Value> =
        serde_json::from_str(&body).expect("diagnostics JSON parses as array");
    assert!(
        diags
            .iter()
            .any(|d| d.get("code").and_then(|v| v.as_str()) == Some("W-LOAD-100")),
        "expected diagnostic with code=W-LOAD-100, got: {}",
        body
    );
}
