//! One place that decides what an absent bulk fixture means.
//!
//! The four large corpora — `tests/fixtures/bulk/*.gz` and
//! `tests/fixtures/validation/*_exhaustive.json.gz` — are not in the git tree.
//! They are hosted as assets on the `test-fixtures-v1` release and fetched by
//! `scripts/fetch-test-fixtures.sh` (and by CI's `fixtures` job, which
//! publishes them to the sharded jobs as a workflow artifact).
//!
//! Eight tests across four modules read them, and every one of them **skips
//! green** when its fixture is missing:
//!
//! ```ignore
//! let buf = match load_fixture_bytes(..) {
//!     Some(b) => b,
//!     None => { eprintln!("Skipping: .. not found"); return; }
//! };
//! ```
//!
//! That is the right behaviour for a developer who has not fetched 156 MB, and
//! exactly the wrong behaviour in CI: a fetch that failed quietly would delete
//! the ClinVar, CMRG and Paraphase coverage while making the job look *faster*,
//! which is the worst outcome available here. A speedup that is really coverage
//! loss is indistinguishable from a real speedup unless something asserts the
//! difference.
//!
//! So the skip is made conditional in one place rather than eight.
//! `FERRO_REQUIRE_BULK_FIXTURES=1` turns every such skip into a panic naming the
//! missing path and the command that fixes it. CI sets it in every job that
//! fetches the fixtures, so absence can no longer masquerade as a pass; it stays
//! unset locally, where skipping is the useful default.

use std::path::PathBuf;

/// The environment variable that promotes a fixture-absence skip to a failure.
pub const REQUIRE_ENV: &str = "FERRO_REQUIRE_BULK_FIXTURES";

/// Whether an absent bulk fixture must fail the test rather than skip it.
///
/// True when `FERRO_REQUIRE_BULK_FIXTURES` is set to anything other than the
/// empty string or `0`. Read on every call rather than cached in a `OnceLock`:
/// these tests are few and slow, the read is trivially cheap beside them, and a
/// cache would make the variable unsettable from within a test — which is what
/// this module's own coverage needs.
pub fn fixtures_are_required() -> bool {
    value_requires_fixtures(std::env::var(REQUIRE_ENV).ok().as_deref())
}

/// The predicate itself, over an explicit value so it is testable without
/// touching the process environment.
fn value_requires_fixtures(value: Option<&str>) -> bool {
    matches!(value, Some(v) if !v.is_empty() && v != "0")
}

/// Called when a bulk fixture is not on disk: panics under
/// `FERRO_REQUIRE_BULK_FIXTURES`, otherwise reports the skip on stderr.
///
/// Returning `()` rather than `!` is deliberate — the call sites keep their
/// `return`, so what a reader sees at each one is still an ordinary early exit
/// rather than control flow that depends on an environment variable.
pub fn absent(path: &str) {
    assert!(
        !fixtures_are_required(),
        "{}",
        missing_fixture_message(path)
    );
    eprintln!("Skipping: {path} not found (run scripts/fetch-test-fixtures.sh)");
}

/// The failure text, factored out so it can be asserted on rather than
/// re-typed in a test.
fn missing_fixture_message(path: &str) -> String {
    format!(
        "{REQUIRE_ENV} is set, but the bulk fixture `{path}` is missing.\n\
         \n\
         This test reads a corpus hosted as a release asset, not committed to git.\n\
         Without it the test would SKIP GREEN, silently dropping its coverage — so\n\
         under {REQUIRE_ENV} an absent fixture is a failure instead.\n\
         \n\
         Fetch and verify the fixtures with:\n\
         \n    scripts/fetch-test-fixtures.sh\n\
         \n\
         In CI this means the `fixtures` job's artifact did not reach this job."
    )
}

/// Read a gzip-compressed bulk fixture's presence, applying the policy above.
///
/// `path` is relative to the WORKSPACE ROOT, and the resolved absolute path is
/// what comes back — so the caller opens the file this function tested for,
/// rather than re-resolving the same string against its own working directory.
///
/// **That is the whole reason it returns a path.** It used to answer `Some(())`
/// and leave each caller to `File::open(path)`, which agreed only while every
/// caller's process ran with the workspace root as its cwd. `#[cfg(test)]`
/// callers no longer all do: nextest sets a test binary's cwd to its own
/// PACKAGE root, and `normalize_axis_preserving` is compiled into
/// `ferro-hgvs-soak-tests` as well as into `it` (see
/// `tests-soak/tests/soak/main.rs`), where the same relative string would name
/// `tests-soak/tests/fixtures/…`. Under `FERRO_REQUIRE_BULK_FIXTURES` that is a
/// loud failure; without it, it is the silent one this module exists to
/// prevent — three corpora reported absent and skipped green.
///
/// `None` after recording the skip (or after panicking, under
/// `FERRO_REQUIRE_BULK_FIXTURES`).
pub fn present_or_skip(path: &str) -> Option<PathBuf> {
    let resolved = crate::common::fixture_gen::fixture_path(path);
    if resolved.exists() {
        return Some(resolved);
    }
    absent(path);
    None
}

// --------------------------------------------------------------------------
// Tests.
//
// Plain `#[test]`, deliberately **not** inside a `#[cfg(test)] mod tests`: this
// tree is an integration-test binary, which compiles without `cfg(test)`, so a
// gated module would never run and would read as coverage it does not provide
// (see `rulings.rs`, which records the same, and the repository `CONTRIBUTING.md` on
// committed tests that have never executed).
// --------------------------------------------------------------------------

/// The message must name the missing path and the fix. A guard whose failure
/// text does not say what to run gets worked around rather than fixed.
#[test]
fn the_missing_fixture_message_names_the_path_and_the_remedy() {
    let message = missing_fixture_message("tests/fixtures/bulk/clinvar_hgvs_500k.json.gz");
    assert!(message.contains("tests/fixtures/bulk/clinvar_hgvs_500k.json.gz"));
    assert!(message.contains("scripts/fetch-test-fixtures.sh"));
    assert!(message.contains(REQUIRE_ENV));
}

/// A fixture that exists is never a skip, whatever the variable says.
///
/// It also pins the resolution: the answer must be the workspace-root copy,
/// which is what makes this test meaningful in `ferro-hgvs-soak-tests`, whose
/// working directory is `tests-soak/` rather than the workspace root. Before
/// `present_or_skip` resolved the path itself, this assertion failed there —
/// which is how the silent-skip hazard on the bulk corpora was found.
#[test]
fn an_existing_bulk_fixture_path_is_present() {
    let resolved = present_or_skip("tests/fixtures/README.md")
        .expect("tests/fixtures/README.md is committed and must resolve");
    assert!(resolved.is_absolute(), "{resolved:?} is not absolute");
    assert!(resolved.is_file(), "{resolved:?} is not a file");
}

/// The predicate's contract, pinned against the spellings CI and a developer
/// shell actually produce. `0` and the empty string are the two ways a variable
/// is "set but off", and reading either as *on* would make a local run fail on
/// a missing 156 MB corpus.
///
/// Uses a private helper over an explicit value rather than the process
/// environment: mutating the environment is a process-global side effect and
/// nextest runs tests in threads, so a test that set the real variable could
/// change what a concurrently running sibling sees.
#[test]
fn only_a_non_empty_non_zero_value_requires_the_fixtures() {
    assert!(!value_requires_fixtures(None));
    assert!(!value_requires_fixtures(Some("")));
    assert!(!value_requires_fixtures(Some("0")));
    assert!(value_requires_fixtures(Some("1")));
    assert!(value_requires_fixtures(Some("true")));
}
