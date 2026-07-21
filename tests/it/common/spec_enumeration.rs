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
//! [`ensure_spec_enumeration`] makes sure that fixture exists first (passed as
//! the prerequisite to [`super::fixture_gen::ensure_generated_fixture`]).

use std::path::PathBuf;

use crate::common::fixture_gen;

const FIXTURE_RELATIVE_PATH: &str = "tests/fixtures/grammar/hgvs_spec_enumeration.json";

/// Absolute path to the generated enumeration fixture.
pub fn spec_enumeration_path() -> PathBuf {
    fixture_gen::fixture_path(FIXTURE_RELATIVE_PATH)
}

/// Ensure the generated enumeration exists, regenerating it when missing.
/// Idempotent and safe to call concurrently. The generator dedups against the
/// normalization fixture, so that one is ensured first.
pub fn ensure_spec_enumeration() {
    fixture_gen::ensure_generated_fixture(
        &spec_enumeration_path(),
        "generate_spec_enumeration",
        ".hgvs_spec_enumeration",
        "spec enumeration",
        crate::common::spec_fixture::ensure_spec_fixture,
    );
}
