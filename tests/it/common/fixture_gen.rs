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

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::{Mutex, OnceLock};

/// Serializes regeneration attempts within a single test binary, across both
/// fixtures. Cross-binary races are made safe by writing to a per-process temp
/// file and atomically renaming it into place. The lock is never held across
/// `dependency` (which is run to completion first), so a single non-reentrant
/// mutex is safe.
static GEN_LOCK: Mutex<()> = Mutex::new(());

/// The workspace root, which every path in this module is relative to.
///
/// **Not `CARGO_MANIFEST_DIR` directly, and the difference is the whole reason
/// this function exists.** That macro names the package this file is compiled
/// INTO, and it is compiled into two of them: the root `ferro-hgvs` package's
/// `it` target, and the `ferro-hgvs-soak-tests` member, which `#[path]`-includes
/// this file (see `tests-soak/tests/soak/main.rs`). For the first it already
/// names the workspace root; for the second it names `<root>/tests-soak`, and
/// every fixture path built from it would be wrong by one directory — silently,
/// because a missing generated fixture is *regenerated* rather than reported,
/// so the failure would be a nested `cargo run` writing into a directory nobody
/// reads.
///
/// The root is found by ascending to the first ancestor holding `Cargo.lock`,
/// which cargo writes at the workspace root and never in a member. Resolved once
/// and cached: it is read on every fixture lookup and the answer cannot change
/// within a process.
fn workspace_root() -> &'static Path {
    static ROOT: OnceLock<PathBuf> = OnceLock::new();
    ROOT.get_or_init(|| {
        let mut dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        loop {
            if dir.join("Cargo.lock").is_file() {
                return dir;
            }
            assert!(
                dir.pop(),
                "no Cargo.lock in any ancestor of {} — the workspace root cannot \
                 be located, so no generated-fixture path can be built",
                env!("CARGO_MANIFEST_DIR")
            );
        }
    })
}

/// Absolute path to a generated fixture, rooted at the workspace root so it is
/// independent of the test's working directory *and* of which package the
/// caller was compiled into.
pub fn fixture_path(relative: &str) -> PathBuf {
    workspace_root().join(relative)
}

/// Ensure `path` exists **and was produced by the current revision of
/// `generator`**, regenerating it via
/// `cargo run --features dev --bin <generator> -- --output <tmp>` when it is
/// missing or stale. `dependency` runs first (before the lock is taken) to
/// satisfy any fixture the generator itself reads; pass `|| {}` when there is
/// none. `label` names the fixture in panic messages. Idempotent and safe to
/// call from many tests concurrently: regeneration is serialized in-process and
/// writes atomically, so a partially written file is never observed.
///
/// See [`ensure_generated`] for what "the current revision" means and why a
/// fixture carrying no revision stamp is trusted rather than rebuilt.
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

/// The filename suffix of the sidecar that records which generator revision
/// produced a fixture. `<fixture>` -> `<fixture>.genrev`. Gitignored alongside
/// the fixtures themselves.
const REVISION_STAMP_SUFFIX: &str = ".genrev";

/// The revision stamp beside `fixture` — the sidecar path, not its contents.
fn revision_stamp_path(fixture: &Path) -> PathBuf {
    let mut name = fixture
        .file_name()
        .expect("fixture path has a file name")
        .to_os_string();
    name.push(REVISION_STAMP_SUFFIX);
    fixture.with_file_name(name)
}

/// A revision token for `generator`, run with `generator_args`.
///
/// It is a hash of the generator's own source file (`examples/<generator>.rs`,
/// which is where every generator this module drives lives, `[[bin]]` and
/// `[[example]]` alike) together with the arguments this fixture is produced
/// with — the same generator backs several corpora that differ only in its
/// flags. Any edit to that source, or a change to the flags a caller passes,
/// moves the token and so invalidates a cached fixture produced by the old one.
///
/// This is deliberately **not** a hash of every transitive input (the spec
/// submodule, the curated overrides, the shared `examples/common/` modules a
/// generator `#[path]`-includes): those are covered by the committed
/// per-consumer guards the fixtures feed, whereas the failure this token closes
/// is the one nothing else observes — reusing a corpus built by an *older
/// revision of its own generator* after a checkout that changed it. `CLAUDE.md`
/// files this under the "stale local artifact" class.
///
/// `DefaultHasher`, not a crypto hash: this needs a value that changes when the
/// source changes, not collision resistance, and staying in `std` keeps the
/// helper free of a dependency the `ferro-hgvs-soak-tests` member (which
/// `#[path]`-includes this file) would otherwise have to carry. `DefaultHasher`
/// uses fixed keys, so the token is stable across runs; were it ever not, a
/// changed token merely forces one redundant regeneration, never a wrong reuse.
fn generator_revision(generator: &str, generator_args: &[&str]) -> String {
    let source = workspace_root()
        .join("examples")
        .join(format!("{generator}.rs"));
    let bytes = std::fs::read(&source).unwrap_or_else(|e| {
        panic!(
            "failed to read generator source {} to stamp its revision: {e}",
            source.display()
        )
    });
    let mut hasher = DefaultHasher::new();
    bytes.hash(&mut hasher);
    for arg in generator_args {
        arg.hash(&mut hasher);
    }
    format!("{generator} {:016x}", hasher.finish())
}

/// Whether the on-disk `fixture` may be reused without regenerating, given the
/// current generator `revision`.
///
/// The cases, and why the "no stamp" one is deliberate:
/// - **absent** -> not fresh; it must be generated.
/// - **present, no stamp** (the sidecar does not exist) -> fresh (trusted). A
///   fixture with no `.genrev` sidecar was produced outside this helper — by
///   CI's `spec-fixtures` job or a manual `cargo run` — and is trusted as-is.
///   This is what keeps CI on the download-and-reuse fast path: the
///   archive-running jobs have no warm `target/`, so treating a stamp-less
///   fixture as stale would rebuild the whole crate inside every test (the exact
///   regression `generated_fixture_ci_wiring.rs` exists to prevent).
/// - **present, stamp unreadable** (a read error other than "not found" — a
///   corrupt or non-UTF-8 sidecar) -> stale. Only a *missing* sidecar is the
///   trusted case above; a sidecar that exists but cannot be compared is
///   ambiguous, and the safe direction for an ambiguous stamp is to regenerate.
/// - **present, stamped** -> fresh only if the stamp still names `revision`. A
///   fixture this helper generated carries a stamp, so a generator change is
///   caught here and forces regeneration.
fn fixture_is_fresh(fixture: &Path, revision: &str) -> bool {
    if !fixture.exists() {
        return false;
    }
    match std::fs::read_to_string(revision_stamp_path(fixture)) {
        Ok(stamped) => stamped.trim() == revision,
        // A missing sidecar is the trusted "no stamp" case; any other error
        // (corrupt bytes, non-UTF-8, permissions) is ambiguous, so regenerate.
        Err(err) => err.kind() == std::io::ErrorKind::NotFound,
    }
}

/// Record which generator revision produced `fixture`, beside it.
///
/// Best-effort: a missing or partially written stamp only forces a future
/// regeneration (a stamp that fails to equal the current revision is treated as
/// stale), never an incorrect reuse — so a write failure here is safe to ignore.
fn write_revision_stamp(fixture: &Path, revision: &str) {
    let _ = std::fs::write(revision_stamp_path(fixture), revision);
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
    let revision = generator_revision(generator, generator_args);

    if fixture_is_fresh(path, &revision) {
        return;
    }

    // Satisfy any prerequisite fixture before taking the lock (the generator
    // may read it), so the lock is never held across the dependency call.
    dependency();

    // Recover from a poisoned lock — a prior panic mid-regeneration must not
    // wedge every other test that needs the fixture.
    let _guard = GEN_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    // Another caller may have regenerated it while we waited for the lock.
    if fixture_is_fresh(path, &revision) {
        return;
    }

    // The workspace root, not this package's directory: the generators are
    // `[[bin]]`/`[[example]]` targets of `ferro-hgvs`, so a nested `cargo run`
    // has to start where that package is the default one. From
    // `ferro-hgvs-soak-tests`'s own directory it would fail with "no bin target
    // named `generate_spec_fixture`".
    let manifest_dir = workspace_root();
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

    // Atomically install the result, replacing any stale copy a previous
    // revision left behind. `rename` within one directory is atomic on POSIX,
    // so a partially written file is never observed; a sibling test binary that
    // won the race is simply overwritten with our equivalent content.
    match std::fs::rename(&tmp, path) {
        // We installed these exact bytes, so stamp them with the revision that
        // produced them — that is what lets a later run tell a fixture built by
        // a different revision of this generator from a current one.
        Ok(()) => write_revision_stamp(path, &revision),
        // The rename failed but the fixture is present: a sibling test binary
        // won the race and installed it. Drop our redundant temp copy, and
        // deliberately do NOT stamp — we did not write those bytes and cannot
        // vouch for their revision. Stamping them `revision` here would be the
        // one way to mark a stale (or foreign) fixture current. Leaving them
        // unstamped keeps freshness self-correcting: the next lookup trusts a
        // stamp-less file (the sibling stamps its own successful install) or
        // regenerates one carrying a stale stamp.
        Err(err) => {
            let _ = std::fs::remove_file(&tmp);
            assert!(path.exists(), "failed to install generated {label}: {err}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The root every fixture path hangs off must be the WORKSPACE root, not
    /// the including package's directory.
    ///
    /// This compiles into both binaries that include this file, so it is
    /// asserted once per package — which is the point: the `it` copy would pass
    /// under the old `CARGO_MANIFEST_DIR` reading too, and the
    /// `ferro-hgvs-soak-tests` copy is the one that would not. A regression
    /// here does not go red on its own, it makes `ensure_spec_fixture`
    /// regenerate into a directory nobody reads.
    #[test]
    fn the_fixture_root_is_the_workspace_root_from_either_package() {
        let root = workspace_root();
        assert!(
            root.join("Cargo.lock").is_file(),
            "{} holds no Cargo.lock",
            root.display()
        );
        let manifest = std::fs::read_to_string(root.join("Cargo.toml"))
            .unwrap_or_else(|e| panic!("read {}/Cargo.toml: {e}", root.display()));
        assert!(
            manifest.lines().any(|line| line.trim() == "[workspace]"),
            "{} holds a Cargo.toml with no [workspace] table, so it is a member \
             directory rather than the workspace root",
            root.display()
        );
        assert!(
            fixture_path("tests/fixtures/grammar").is_dir(),
            "the generated-fixture directory does not resolve from {}",
            root.display()
        );
    }

    /// A unique, non-existent fixture path in the OS temp dir, plus its sidecar.
    fn scratch_fixture(tag: &str) -> (PathBuf, PathBuf) {
        let base = std::env::temp_dir().join(format!(
            "ferro-fixture-gen-{}-{}-{tag}.json",
            std::process::id(),
            // A monotonic-enough disambiguator so two cases in one process do
            // not collide on the same path.
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        let stamp = revision_stamp_path(&base);
        let _ = std::fs::remove_file(&base);
        let _ = std::fs::remove_file(&stamp);
        (base, stamp)
    }

    /// The revision token is a pure function of the generator source and args:
    /// stable across calls, and it moves when the args move.
    #[test]
    fn a_generator_revision_is_deterministic_and_arg_sensitive() {
        // `generate_cis_confluence_corpus` is a real generator this module
        // drives, so its source resolves.
        let a = generator_revision("generate_cis_confluence_corpus", &["--axes", "g,c"]);
        let b = generator_revision("generate_cis_confluence_corpus", &["--axes", "g,c"]);
        assert_eq!(
            a, b,
            "the same generator + args must hash to the same token"
        );

        let c = generator_revision("generate_cis_confluence_corpus", &["--axes", "n,r"]);
        assert_ne!(
            a, c,
            "a change in the generator's arguments must move its revision token"
        );

        assert!(
            a.starts_with("generate_cis_confluence_corpus "),
            "the token should name its generator for readability, got {a:?}"
        );
    }

    /// The freshness decision — the whole point of the revision check.
    ///
    /// A cached fixture whose stamp names a *different* generator revision is
    /// stale and must be regenerated; one whose stamp matches is reused. This is
    /// the property #1654 asked for, exercised without paying for the (expensive,
    /// `--features dev`) subprocess by testing the decision directly.
    #[test]
    fn a_stamp_mismatch_forces_regeneration_and_a_match_reuses_the_cache() {
        let (fixture, stamp) = scratch_fixture("freshness");

        // Absent fixture: never fresh.
        assert!(
            !fixture_is_fresh(&fixture, "rev-A"),
            "a missing fixture must not be considered fresh"
        );

        std::fs::write(&fixture, b"corpus produced at revision A").expect("write scratch fixture");

        // Present but stamped with an older revision: stale -> regenerate.
        std::fs::write(&stamp, "rev-A").expect("write scratch stamp");
        assert!(
            !fixture_is_fresh(&fixture, "rev-B"),
            "a fixture stamped with a different revision must be regenerated"
        );

        // Present and stamped with the current revision: fresh -> reuse.
        assert!(
            fixture_is_fresh(&fixture, "rev-A"),
            "a fixture stamped with the current revision must be reused"
        );

        std::fs::remove_file(&fixture).ok();
        std::fs::remove_file(&stamp).ok();
    }

    /// A stamp-less fixture is trusted, so CI's downloaded artifacts (and manual
    /// `cargo run` output) stay on the reuse fast path rather than triggering a
    /// crate rebuild inside every test.
    #[test]
    fn a_fixture_with_no_stamp_is_trusted() {
        let (fixture, stamp) = scratch_fixture("no-stamp");
        std::fs::write(&fixture, b"externally produced corpus").expect("write scratch fixture");
        assert!(
            !stamp.exists(),
            "precondition: this fixture has no revision stamp"
        );
        assert!(
            fixture_is_fresh(&fixture, "any-revision"),
            "a fixture with no stamp must be reused, not rebuilt (CI fast path)"
        );
        std::fs::remove_file(&fixture).ok();
    }

    /// An unreadable stamp is treated as stale, not trusted. Only a *missing*
    /// sidecar is the trusted "no stamp" fast path; a sidecar that exists but
    /// cannot be compared (corrupt / non-UTF-8) must force regeneration, or a
    /// stale fixture could bypass the revision check on a garbled stamp.
    #[test]
    fn an_unreadable_stamp_forces_regeneration() {
        let (fixture, stamp) = scratch_fixture("corrupt-stamp");
        std::fs::write(&fixture, b"corpus").expect("write scratch fixture");

        // A non-UTF-8 sidecar: present, but `read_to_string` cannot decode it.
        std::fs::write(&stamp, [0xff, 0xfe, 0x00, 0x9f]).expect("write corrupt stamp");
        assert!(
            stamp.exists(),
            "precondition: the stamp exists but is not valid UTF-8"
        );
        assert!(
            !fixture_is_fresh(&fixture, "any-revision"),
            "a present-but-unreadable stamp must be treated as stale"
        );

        std::fs::remove_file(&fixture).ok();
        std::fs::remove_file(&stamp).ok();
    }

    /// `write_revision_stamp` then `fixture_is_fresh` round-trips: the value a
    /// generation writes is exactly what a later run compares against.
    #[test]
    fn writing_a_stamp_then_reading_it_agrees() {
        let (fixture, stamp) = scratch_fixture("roundtrip");
        std::fs::write(&fixture, b"corpus").expect("write scratch fixture");

        let revision = generator_revision("generate_spec_fixture", &[]);
        write_revision_stamp(&fixture, &revision);

        assert!(stamp.exists(), "the stamp sidecar should have been written");
        assert!(
            fixture_is_fresh(&fixture, &revision),
            "the just-written stamp must read back as fresh"
        );
        assert!(
            !fixture_is_fresh(&fixture, "some-other-revision"),
            "a different revision must still read as stale after a write"
        );

        std::fs::remove_file(&fixture).ok();
        std::fs::remove_file(&stamp).ok();
    }
}
