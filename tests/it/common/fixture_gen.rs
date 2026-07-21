//! Shared on-demand regeneration for the gitignored generated spec fixtures.
//!
//! Both `hgvs_spec_normalization.json` (see [`super::spec_fixture`]) and
//! `hgvs_spec_enumeration.json` (see [`super::spec_enumeration`]) are generated
//! build artifacts, not committed files (see `CLAUDE.md`): tracking them made
//! every parser PR a merge-conflict magnet. Each is produced by running a
//! `--features dev` example with `--output <tmp>` and atomically renaming the
//! result into place. This module holds the one regeneration flow they share —
//! locking, subprocess execution, temp-file cleanup, atomic rename — so the two
//! callers are thin wrappers that differ only in path, example name, temp stem,
//! and an optional prerequisite fixture to satisfy first.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Mutex;

/// Serializes regeneration attempts within a single test binary, across both
/// fixtures. Cross-binary races are made safe by writing to a per-process temp
/// file and atomically renaming it into place. The lock is never held across
/// `dependency` (which is run to completion first), so a single non-reentrant
/// mutex is safe.
static GEN_LOCK: Mutex<()> = Mutex::new(());

/// Absolute path to a generated fixture, rooted at `CARGO_MANIFEST_DIR` so it is
/// independent of the test's working directory.
pub fn fixture_path(relative: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push(relative);
    p
}

/// Ensure `path` exists, regenerating it via
/// `cargo run --features dev --example <example> -- --output <tmp>` when it is
/// missing. `dependency` runs first (before the lock is taken) to satisfy any
/// fixture the generator itself reads; pass `|| {}` when there is none. `label`
/// names the fixture in panic messages. Idempotent and safe to call from many
/// tests concurrently: regeneration is serialized in-process and writes
/// atomically, so a partially written file is never observed.
pub fn ensure_generated_fixture(
    path: &Path,
    example: &str,
    tmp_stem: &str,
    label: &str,
    dependency: impl FnOnce(),
) {
    if path.exists() {
        return;
    }

    // Satisfy any prerequisite fixture before taking the lock (the generator
    // may read it), so the lock is never held across the dependency call.
    dependency();

    // Recover from a poisoned lock — a prior panic mid-regeneration must not
    // wedge every other test that needs the fixture.
    let _guard = GEN_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    // Another caller may have generated it while we waited for the lock.
    if path.exists() {
        return;
    }

    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let dir = path.parent().expect("fixture path has a parent directory");
    let tmp = dir.join(format!("{tmp_stem}.{}.tmp", std::process::id()));

    let status = Command::new(env!("CARGO"))
        .current_dir(manifest_dir)
        .args([
            "run",
            "--quiet",
            "--features",
            "dev",
            "--example",
            example,
            "--",
            "--output",
        ])
        .arg(&tmp)
        .status()
        .unwrap_or_else(|e| panic!("failed to run `{example}` example: {e}"));
    assert!(status.success(), "`{example}` exited with failure");

    // Atomically install the result. If a sibling test binary won the race and
    // already created the final file, drop our redundant temp copy.
    if path.exists() {
        let _ = std::fs::remove_file(&tmp);
    } else if let Err(err) = std::fs::rename(&tmp, path) {
        let _ = std::fs::remove_file(&tmp);
        assert!(path.exists(), "failed to install generated {label}: {err}");
    }
}
