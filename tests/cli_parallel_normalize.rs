//! Integration tests for parallel batch `ferro normalize` (`-j/--workers`).
//!
//! The core correctness property of the parallel batch driver is **invariance**:
//! output must be byte-identical to the serial path regardless of worker count,
//! including input ordering across chunk boundaries. These tests run under the
//! mock provider (no `--reference`), so they need no external data and run in CI.
//! Under the mock provider `normalize` echoes each input variant verbatim, so a
//! reordering bug would surface as a differing stdout byte stream.

use std::io::Write;
use std::process::Command;
use tempfile::NamedTempFile;

/// Mirror of `BATCH_CHUNK_LINES` in `src/bin/ferro.rs` (the batch driver's chunk
/// size). Integration tests cannot import items from a `bin` target, so this is
/// duplicated here; the runtime assertions in each test below verify that the
/// chosen variant count still spans multiple chunks, so a drift in the real
/// constant surfaces as a test failure rather than silently degrading coverage.
const BATCH_CHUNK_LINES: usize = 8192;

/// Build an input with `n` distinct variants interleaved with comments and blank
/// lines (which the driver must skip identically in serial and parallel). `n` is
/// chosen larger than the driver's chunk size so ordering is exercised across
/// multiple chunks.
fn write_input(n: usize) -> (NamedTempFile, std::path::PathBuf) {
    let mut tf = tempfile::Builder::new().suffix(".txt").tempfile().unwrap();
    for i in 1..=n {
        if i % 100 == 0 {
            writeln!(tf, "# comment at {i}").unwrap();
        }
        if i % 50 == 0 {
            writeln!(tf).unwrap(); // blank line
        }
        // Distinct per line so any reordering changes the output byte stream.
        writeln!(tf, "NM_000088.3:c.{i}A>G").unwrap();
    }
    tf.flush().unwrap();
    let path = tf.path().to_path_buf();
    (tf, path)
}

/// Run `ferro normalize` over `input` with the given worker count, returning
/// stdout bytes. Mock provider (no `--reference`).
fn run_normalize(input: &std::path::Path, workers: usize) -> Vec<u8> {
    let bin = env!("CARGO_BIN_EXE_ferro");
    let out = Command::new(bin)
        .args([
            "normalize",
            "--error-mode",
            "silent",
            "-j",
            &workers.to_string(),
            "-i",
            input.to_str().unwrap(),
        ])
        .output()
        .expect("run ferro normalize");
    assert!(
        out.status.success(),
        "ferro normalize -j{workers} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

#[test]
fn parallel_output_is_byte_identical_to_serial() {
    // Choose a variant count that spans multiple chunks so ordering is exercised
    // across chunk boundaries. Asserting it (rather than relying on a comment)
    // means coverage degrades loudly if BATCH_CHUNK_LINES grows past this.
    let n = 20_000;
    assert!(
        n > 2 * BATCH_CHUNK_LINES,
        "n={n} no longer spans multiple chunks (BATCH_CHUNK_LINES={BATCH_CHUNK_LINES}); \
         raise n so ordering is exercised across chunk boundaries"
    );
    let (_keep, path) = write_input(n);

    let serial = run_normalize(&path, 1);
    assert!(!serial.is_empty(), "serial run produced no output");

    for workers in [2usize, 4, 8] {
        let parallel = run_normalize(&path, workers);
        assert_eq!(
            serial, parallel,
            "stdout for -j{workers} differs from serial (-j1): ordering or formatting diverged"
        );
    }
}

#[test]
fn parallel_runs_are_deterministic() {
    let n = 12_000;
    assert!(
        n > BATCH_CHUNK_LINES,
        "n={n} no longer spans multiple chunks (BATCH_CHUNK_LINES={BATCH_CHUNK_LINES})"
    );
    let (_keep, path) = write_input(n);
    let a = run_normalize(&path, 8);
    let b = run_normalize(&path, 8);
    assert_eq!(
        a, b,
        "two -j8 runs produced different output (nondeterministic)"
    );
    assert!(!a.is_empty());
}
