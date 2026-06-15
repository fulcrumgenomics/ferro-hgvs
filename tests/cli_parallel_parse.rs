//! Integration tests for parallel batch `ferro parse` (`-j/--workers`).
//!
//! Like the parallel `normalize` tests, the core property is **invariance**:
//! output must be byte-identical to the serial path regardless of worker count,
//! including ordering across chunk boundaries. `parse` needs no reference data,
//! so these tests are pure-CPU and fully self-contained.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Build an input with `n` distinct variants interleaved with comments and blank
/// lines (which the driver must skip identically in serial and parallel). `n` is
/// larger than the driver's chunk size so ordering spans multiple chunks.
fn write_input(n: usize) -> (NamedTempFile, std::path::PathBuf) {
    let mut tf = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    for i in 1..=n {
        if i % 100 == 0 {
            writeln!(tf, "# comment at {i}").unwrap();
        }
        if i % 50 == 0 {
            writeln!(tf).unwrap();
        }
        writeln!(tf, "NM_000088.3:c.{i}A>G").unwrap();
    }
    tf.flush().unwrap();
    let path = tf.path().to_path_buf();
    (tf, path)
}

/// Run `ferro parse` over `input` with the given worker count, returning stdout.
fn run_parse(input: &std::path::Path, workers: usize) -> Vec<u8> {
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "parse",
            "--error-mode",
            "silent",
            "-j",
            &workers.to_string(),
            "-i",
            input.to_str().unwrap(),
        ])
        .output()
        .expect("run ferro parse");
    assert!(
        out.status.success(),
        "ferro parse -j{workers} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

#[test]
fn parallel_parse_output_is_byte_identical_to_serial() {
    let (_keep, path) = write_input(20_000);
    let serial = run_parse(&path, 1);
    assert!(!serial.is_empty(), "serial run produced no output");
    for workers in [2usize, 4, 8] {
        let parallel = run_parse(&path, workers);
        assert_eq!(
            serial, parallel,
            "parse stdout for -j{workers} differs from serial (-j1)"
        );
    }
}

#[test]
fn parallel_parse_is_deterministic() {
    let (_keep, path) = write_input(12_000);
    let a = run_parse(&path, 8);
    let b = run_parse(&path, 8);
    assert_eq!(a, b, "two -j8 parse runs differ (nondeterministic)");
    assert!(!a.is_empty());
}
