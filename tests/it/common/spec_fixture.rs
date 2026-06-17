//! On-demand regeneration of the HGVS spec-normalization fixture.
//!
//! `tests/fixtures/grammar/hgvs_spec_normalization.json` is a generated build
//! artifact, not a committed file (see `CLAUDE.md`): it is produced by the
//! `generate_spec_fixture` example from the spec submodule, the parser's
//! behavior, and the curated `hgvs_spec_normalization_overrides.json`. Tracking
//! it made every parser PR a merge-conflict magnet, so it is regenerated
//! instead.
//!
//! Tests that read the fixture call [`ensure_spec_fixture`] first. CI
//! regenerates it explicitly before the test run, so in CI the subprocess path
//! below is never taken; it exists as the local-dev safety net for a fresh
//! checkout.

use std::path::PathBuf;
use std::process::Command;
use std::sync::Mutex;

/// Serializes regeneration attempts within a single test binary. Cross-binary
/// races are made safe by writing to a per-process temp file and atomically
/// renaming it into place (see [`ensure_spec_fixture`]).
static GEN_LOCK: Mutex<()> = Mutex::new(());

/// Absolute path to the generated spec-normalization fixture.
pub fn spec_fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/fixtures/grammar/hgvs_spec_normalization.json");
    p
}

/// Ensure the generated spec fixture exists, regenerating it via the
/// `generate_spec_fixture` example when it is missing. Idempotent and safe to
/// call from many tests concurrently: regeneration is serialized in-process and
/// writes atomically, so a partially written file is never observed.
pub fn ensure_spec_fixture() {
    let path = spec_fixture_path();
    if path.exists() {
        return;
    }

    // Recover from a poisoned lock — a prior panic mid-regeneration must not
    // wedge every other test that needs the fixture.
    let _guard = GEN_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    // Another caller may have generated it while we waited for the lock.
    if path.exists() {
        return;
    }

    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let dir = path.parent().expect("fixture path has a parent directory");
    let tmp = dir.join(format!(
        ".hgvs_spec_normalization.{}.tmp",
        std::process::id()
    ));

    let status = Command::new(env!("CARGO"))
        .current_dir(manifest_dir)
        .args([
            "run",
            "--quiet",
            "--features",
            "dev",
            "--example",
            "generate_spec_fixture",
            "--",
            "--output",
        ])
        .arg(&tmp)
        .status()
        .expect("failed to run `generate_spec_fixture` example");
    assert!(
        status.success(),
        "`generate_spec_fixture` exited with failure"
    );

    // Atomically install the result. If a sibling test binary won the race and
    // already created the final file, drop our redundant temp copy.
    if path.exists() {
        let _ = std::fs::remove_file(&tmp);
    } else if let Err(err) = std::fs::rename(&tmp, &path) {
        let _ = std::fs::remove_file(&tmp);
        assert!(
            path.exists(),
            "failed to install generated spec fixture: {err}"
        );
    }
}
