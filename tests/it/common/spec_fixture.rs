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
//! (in [`super::fixture_gen`]) is never taken; it exists as the local-dev
//! safety net for a fresh checkout.

use std::path::PathBuf;

use crate::common::fixture_gen;

const FIXTURE_RELATIVE_PATH: &str = "tests/fixtures/grammar/hgvs_spec_normalization.json";

/// Absolute path to the generated spec-normalization fixture.
pub fn spec_fixture_path() -> PathBuf {
    fixture_gen::fixture_path(FIXTURE_RELATIVE_PATH)
}

/// Ensure the generated spec fixture exists, regenerating it via the
/// `generate_spec_fixture` example when it is missing. Idempotent and safe to
/// call from many tests concurrently.
pub fn ensure_spec_fixture() {
    fixture_gen::ensure_generated_fixture(
        &spec_fixture_path(),
        "generate_spec_fixture",
        ".hgvs_spec_normalization",
        "spec fixture",
        || {},
    );
}
