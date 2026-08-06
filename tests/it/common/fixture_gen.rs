//! Shared on-demand regeneration for the gitignored generated spec fixtures.
//!
//! Both `hgvs_spec_normalization.json` (see [`super::spec_fixture`]) and
//! `hgvs_spec_enumeration.json` (see [`super::spec_enumeration`]) are generated
//! build artifacts, not committed files (see `CLAUDE.md`): tracking them made
//! every parser PR a merge-conflict magnet. Each is produced by running a
//! `--features dev` generator binary with `--output <tmp>` and atomically renaming the
//! result into place. This module holds the one regeneration flow they share —
//! locking, subprocess execution, temp-file cleanup, atomic rename — so the two
//! callers are thin wrappers that differ only in path, generator name, temp stem,
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
/// `cargo run --features dev --bin <generator> -- --output <tmp>` when it is
/// missing. `dependency` runs first (before the lock is taken) to satisfy any
/// fixture the generator itself reads; pass `|| {}` when there is none. `label`
/// names the fixture in panic messages. Idempotent and safe to call from many
/// tests concurrently: regeneration is serialized in-process and writes
/// atomically, so a partially written file is never observed.
pub fn ensure_generated_fixture(
    path: &Path,
    generator: &str,
    tmp_stem: &str,
    label: &str,
    dependency: impl FnOnce(),
) {
    ensure_generated(path, "--bin", generator, &[], tmp_stem, label, dependency);
}

/// [`ensure_generated_fixture`], for a generator that is an `[[example]]`
/// rather than a `[[bin]]`.
///
/// The two spec generators are binaries for a reason their `Cargo.toml` block
/// spells out — tests must locate them with `CARGO_BIN_EXE_*` — and that reason
/// does not apply to a generator no test resolves by path. So the corpus
/// generators stay examples, and the regeneration flow takes the target kind
/// instead of assuming one. `--example` and `--bin` are the only difference;
/// everything that makes this safe to call concurrently is shared.
///
/// `generator_args` are passed through ahead of `--output`, so one generator can
/// back several fixtures that differ only in its flags — which is how the `g,c`
/// and the `n,r` confluence corpora are both produced by
/// `generate_cis_confluence_corpus`. Each such fixture needs its **own** path and
/// `tmp_stem`; reusing either across two argument sets would let whichever test
/// ran first install its corpus under the other's name.
pub fn ensure_generated_example_fixture(
    path: &Path,
    generator: &str,
    generator_args: &[&str],
    tmp_stem: &str,
    label: &str,
    dependency: impl FnOnce(),
) {
    ensure_generated(
        path,
        "--example",
        generator,
        generator_args,
        tmp_stem,
        label,
        dependency,
    );
}

fn ensure_generated(
    path: &Path,
    target_kind: &str,
    generator: &str,
    generator_args: &[&str],
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

    // Deliberately a nested `cargo run`, NOT `CARGO_BIN_EXE_<generator>`.
    //
    // The generator-executing test modules (`spec_generator_preconditions`,
    // `issue_1046_generator_gitdir_leak`) do use `CARGO_BIN_EXE_*`, and should:
    // they are `#[cfg(feature = "dev")]`. This helper is not. It is reached from
    // `hgvs_spec_normalization_tests`, `idempotency_tests` and
    // `spec_enumeration_tests`, which are ungated so that a plain `cargo test`
    // still runs them. `env!` resolves at *compile* time and cargo defines
    // `CARGO_BIN_EXE_*` only for bins it builds — and both generators are
    // `required-features = ["dev"]` — so expanding it here would fail to compile
    // the whole `it` target without `--features dev`. `generator` is also a
    // runtime string, which `env!` cannot take. Resolving the generator at run
    // time is what keeps this callable from ungated modules.
    let status = Command::new(env!("CARGO"))
        .current_dir(manifest_dir)
        .args([
            "run",
            "--quiet",
            "--features",
            "dev",
            target_kind,
            generator,
            "--",
        ])
        .args(generator_args)
        .arg("--output")
        .arg(&tmp)
        .status()
        .unwrap_or_else(|e| panic!("failed to run `{generator}` binary: {e}"));
    assert!(status.success(), "`{generator}` exited with failure");

    // Atomically install the result. If a sibling test binary won the race and
    // already created the final file, drop our redundant temp copy.
    if path.exists() {
        let _ = std::fs::remove_file(&tmp);
    } else if let Err(err) = std::fs::rename(&tmp, path) {
        let _ = std::fs::remove_file(&tmp);
        assert!(path.exists(), "failed to install generated {label}: {err}");
    }
}
