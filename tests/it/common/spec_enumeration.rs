//! On-demand regeneration of the HGVS spec **test enumeration**.
//!
//! `tests/fixtures/grammar/hgvs_spec_enumeration.json` is a generated build
//! artifact, not a committed file — same contract as
//! `hgvs_spec_normalization.json` (see `CLAUDE.md`). It is produced by the
//! `generate_spec_enumeration` example from the spec submodule, ferro's
//! behaviour, and the curated `hgvs_spec_enumeration_overrides.json`. Tracking
//! it would make every parser PR a merge-conflict magnet, so it is regenerated
//! instead. Only the overrides, the generator and the driver are committed.
//!
//! The generator reads `hgvs_spec_normalization.json` for deduplication, so
//! [`ensure_spec_enumeration`] makes sure that fixture exists first.

use std::path::PathBuf;
use std::process::Command;
use std::sync::Mutex;

/// Serializes regeneration attempts within a single test binary. Cross-binary
/// races are made safe by writing to a per-process temp file and atomically
/// renaming it into place.
static GEN_LOCK: Mutex<()> = Mutex::new(());

/// Absolute path to the generated enumeration fixture.
pub fn spec_enumeration_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests/fixtures/grammar/hgvs_spec_enumeration.json");
    p
}

/// Ensure the generated enumeration exists, regenerating it when missing.
/// Idempotent and safe to call concurrently.
pub fn ensure_spec_enumeration() {
    let path = spec_enumeration_path();
    if path.exists() {
        return;
    }

    // The generator dedups against the normalization fixture, so that one has
    // to exist first.
    crate::common::spec_fixture::ensure_spec_fixture();

    // Recover from a poisoned lock — a prior panic mid-regeneration must not
    // wedge every other test that needs the enumeration.
    let _guard = GEN_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    if path.exists() {
        return;
    }

    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let dir = path.parent().expect("fixture path has a parent directory");
    let tmp = dir.join(format!(".hgvs_spec_enumeration.{}.tmp", std::process::id()));

    let status = Command::new(env!("CARGO"))
        .current_dir(manifest_dir)
        .args([
            "run",
            "--quiet",
            "--features",
            "dev",
            "--example",
            "generate_spec_enumeration",
            "--",
            "--output",
        ])
        .arg(&tmp)
        .status()
        .expect("failed to run `generate_spec_enumeration` example");
    assert!(
        status.success(),
        "`generate_spec_enumeration` exited with failure"
    );

    if path.exists() {
        let _ = std::fs::remove_file(&tmp);
    } else if let Err(err) = std::fs::rename(&tmp, &path) {
        let _ = std::fs::remove_file(&tmp);
        assert!(
            path.exists(),
            "failed to install generated spec enumeration: {err}"
        );
    }
}
