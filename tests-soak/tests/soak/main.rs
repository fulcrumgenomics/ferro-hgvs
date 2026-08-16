//! The test driver for the three CI jobs that run off the optimized archive:
//! `soak`, `sweeps` and `censuses`.
//!
//! # Why this binary exists
//!
//! `test-build-soak` builds one archive at `opt-level = 3` because the work
//! those three jobs do is compute-bound — eight shards of 125,000 proptest
//! cases, three exhaustive cross-products, three censuses — and a constant
//! factor on the code under test is multiplied a million times over. It used to
//! build that archive as `--test it`, which compiles the WHOLE 139k-line `it`
//! crate (6,217 tests) at `opt-level = 3` so that about 45 of them can run.
//!
//! Measured on this machine in `user` CPU seconds, since wall clock here is
//! not reproducible: rebuilding the `it` binary alone costs 251.9s while
//! rebuilding `libferro_hgvs` and then `it` costs 368.6s — so the library is
//! ~32% of that job and the test driver ~68%. The runtime `opt-level = 3` buys
//! comes from optimizing the LIBRARY; the driver is loops and `format!`.
//!
//! Cargo profile overrides are keyed on the package and not on the target, so
//! `[profile.soak.package.ferro-hgvs]` reaches `ferro-hgvs`'s lib, its bins and
//! its `it` test target alike. A separate package is therefore the only unit
//! that can be addressed on its own — hence this one, and hence
//! `[profile.soak.package.ferro-hgvs-soak-tests]` in the root manifest.
//!
//! # Why the modules are `#[path]`-included rather than moved here
//!
//! Every module below is compiled from its file under `tests/it/`, unmoved, and
//! is STILL declared in `tests/it/main.rs`. That is deliberate on three counts:
//!
//! * **The guards that scan by directory keep working.** `sweep_filter_invariant`,
//!   `oracle_exclude_invariant`, `generated_fixture_ci_wiring`,
//!   `ruling_citation_currency` and `ruling_guard_field` all enumerate `.rs`
//!   files under `tests/it/` (or under `tests/`) and each fails in the
//!   *flattering* direction — a module that leaves the scan simply stops being
//!   policed. Moving `spec_conformance_axis.rs`, which carries 54 spec
//!   citations, out of `tests/` would have silently emptied two of them.
//! * **The ruling ledger cites test files by path.** Two `decided` records name
//!   guards in `tests/it/cis_junction_crossing_shift.rs` and
//!   `tests/it/spec_conformance_axis.rs`; moving the files would have meant
//!   editing the adjudication ledger to chase a build-layout change.
//! * **`issue_1615_denoted_sequence_oracle` is run twice on purpose** — un-armed
//!   in `test`, armed in `sweeps` (see `ci.yml`). It stays in the `it` binary
//!   for the first and is compiled here for the second, which is exactly what
//!   an in-place include gives for free.
//!
//! A fourth consequence is worth stating because it is the one that would have
//! failed silently: **the pinned proptest seeds keep replaying, unmoved.**
//! `tests/proptest-regressions/*.txt` are RNG seeds resolved from `file!()`, and
//! proptest walks up from the source file to the directory holding `main.rs`,
//! then writes its `proptest-regressions` sibling. Because `#[path]` leaves the
//! `..` segments in `file!()` and the filesystem resolves them, this binary
//! resolves
//! `tests-soak/tests/soak/../../../tests/proptest-regressions/<module>.txt` —
//! the same committed file `it` reads. Verified by appending an unparsable line
//! and reading back the path proptest named in its warning. Had they moved,
//! nothing would have gone red: those seeds are checked in NO other job, since
//! `test` and `test-oracle` both negate `test(proptest)`.
//!
//! The cost of the choice is that membership here is a second list that can
//! drift out of step with `ci.yml`'s filters. That is closed by
//! `tests/it/soak_package_membership.rs`, which fails when a module CI selects
//! from this archive is not compiled into it — the failure that would otherwise
//! be silent, because such a module is negated by `test`/`test-oracle`, absent
//! from this binary, and therefore run nowhere while CI gets faster.
//!
//! # Running it locally
//!
//! It is not in `default-members`, so a plain `cargo nextest run --features dev`
//! does not build it and the local suite is unchanged — those same modules still
//! run from the `it` binary. To exercise this binary itself:
//!
//! ```text
//! cargo nextest run -p ferro-hgvs-soak-tests
//! cargo nextest archive --cargo-profile soak -p ferro-hgvs-soak-tests \
//!     --archive-file nextest-soak.tar.zst
//! ```

// The shared helpers the modules below reach through `crate::common::…`,
// `#[path]`-included from `tests/it/common/` for the same reason the test
// modules are. `#[path]` on a top-level `mod` resolves against the directory of
// the file declaring it, i.e. `tests-soak/tests/soak/`.
//
// `allow(dead_code)` for the reason `tests/it/common/mod.rs` gives: each binary
// that includes the tree compiles all of it and uses a different subset, so the
// unused-item warnings are per-binary noise. This binary uses a much smaller
// subset than `it` does.
#[allow(dead_code)]
#[path = "../../../tests/it/common/bulk_fixtures.rs"]
mod bulk_fixtures;
#[allow(dead_code)]
#[path = "../../../tests/it/common/cis_apply_oracle.rs"]
mod cis_apply_oracle;
#[allow(dead_code)]
#[path = "../../../tests/it/common/failure_expectations.rs"]
mod failure_expectations;
#[allow(dead_code)]
#[path = "../../../tests/it/common/fixture_gen.rs"]
mod fixture_gen;
#[allow(dead_code)]
#[path = "../../../tests/it/common/hg38_window.rs"]
mod hg38_window;
#[allow(dead_code)]
#[path = "../../../tests/it/common/manifest.rs"]
mod manifest;
#[allow(dead_code)]
#[path = "../../../tests/it/common/minimal_alignment.rs"]
mod minimal_alignment;
#[allow(dead_code)]
#[path = "../../../tests/it/common/spec_fixture.rs"]
mod spec_fixture;
#[allow(dead_code)]
#[path = "../../../tests/it/common/synthetic.rs"]
mod synthetic;

/// The `crate::common::<helper>` façade the included files are written against.
///
/// They are compiled here byte-for-byte as `tests/it` compiles them, so the
/// path they use has to resolve here too. Re-exported rather than declared
/// inside this module so that every `#[path]` above is relative to one
/// directory instead of two.
mod common {
    pub(crate) use super::{
        bulk_fixtures, cis_apply_oracle, failure_expectations, fixture_gen, hg38_window, manifest,
        minimal_alignment, spec_fixture, synthetic,
    };
}

// ---------------------------------------------------------------------------
// The modules the three jobs select. Each list is kept in step with `ci.yml` by
// `tests/it/soak_package_membership.rs`.
// ---------------------------------------------------------------------------

// `soak` — `-E 'test(proptest)'`, eight shards of 125,000 cases each.
#[path = "../../../tests/it/adjacency_confluence_proptest.rs"]
mod adjacency_confluence_proptest;
#[path = "../../../tests/it/cis_allele_confluence_proptest.rs"]
mod cis_allele_confluence_proptest;
#[path = "../../../tests/it/from_sequences_proptest.rs"]
mod from_sequences_proptest;
#[path = "../../../tests/it/minimal_alignment_enumeration_proptest.rs"]
mod minimal_alignment_enumeration_proptest;
#[path = "../../../tests/it/normalize_idempotency_proptest.rs"]
mod normalize_idempotency_proptest;

// `sweeps` — `-E "$SWEEP_FILTER + test(issue_1615_denoted_sequence_oracle)"`.
#[path = "../../../tests/it/cis_junction_crossing_shift.rs"]
mod cis_junction_crossing_shift;
#[path = "../../../tests/it/issue_1254_sibling_crossing_shift.rs"]
mod issue_1254_sibling_crossing_shift;
#[path = "../../../tests/it/issue_1615_denoted_sequence_oracle.rs"]
mod issue_1615_denoted_sequence_oracle;
#[path = "../../../tests/it/repeat_span_sibling_overlap.rs"]
mod repeat_span_sibling_overlap;

// `censuses` — `-E "$CENSUS_FILTER"`, split over the job's two steps by whether
// the seam oracles may be armed.
#[path = "../../../tests/it/issue_1542_direction_symmetry.rs"]
mod issue_1542_direction_symmetry;
#[path = "../../../tests/it/normalize_axis_preserving.rs"]
mod normalize_axis_preserving;
#[path = "../../../tests/it/spec_conformance_axis.rs"]
mod spec_conformance_axis;
// The run-backed census instrument tests (#2094), carved out of
// `conformance_census_instrument` so the `censuses` job owns them instead of the
// `Test` shards. Selected by the job's UN-armed step, beside `spec_conformance_axis`
// and for the same reason: `run_census` refuses with a seam oracle armed.
#[path = "../../../tests/it/conformance_census_runs.rs"]
mod conformance_census_runs;
// A parser benchmark rather than a normalization census — it never reaches
// `canonicalize_from_sequence` — which is exactly why it belongs here: its cost
// is a constant factor on gzip decode, JSON deserialization and `parse_hgvs`,
// and the soak profile lowers all three. It was the second-heaviest test in
// `Test oracle (2/4)` at 54.9s, on the shard that ended that job last.
#[path = "../../../tests/it/clinvar_hgvs_tests.rs"]
mod clinvar_hgvs_tests;

// ---------------------------------------------------------------------------
// SUPPORT: compiled because a module above reaches into it, NOT because any job
// selects it. Its own tests run from `it`, where they always have.
//
// `tests/it/soak_package_membership.rs` distinguishes the two categories by
// derivation rather than by a marker here: anything this file compiles that no
// job selects must be named as `crate::<module>::` by something else this file
// compiles. So the exemption cannot be claimed, only earned — and it lapses the
// moment the reference goes away, which is what keeps this from becoming a
// place to park dead compile cost on the critical path.
// ---------------------------------------------------------------------------

// `minimal_alignment_enumeration_proptest::the_closed_form_and_the_enumerator_agree`
// cross-checks the enumerator against
// `issue_1539_split_member_separation::forced_unchanged_columns`, the Θ(n·m)
// closed-form detector.
//
// THE TWO IMPLEMENTATIONS ARE KEPT APART ON PURPOSE, so the obvious tidy-up —
// move the function next to its sibling in `common/minimal_alignment.rs` and
// drop this include — is the one thing not to do. `common/minimal_alignment.rs`
// says so in its own doc, and the `unchanged-is-read-over-every-minimal-alignment`
// ruling names this symbol BY PATH as the detector the ruling came from. That
// path appears in the ledger's rationale, in the blessed
// `docs/NORMALIZATION_CONTRACT.md` rendered from it, and in the helper's doc —
// so relocating the function to satisfy a build-layout constraint would falsify
// a `decided` record's prose in three places and force a contract re-bless.
// Compiling 1117 lines of test module is the cheaper side of that trade by a
// wide margin.
#[path = "../../../tests/it/issue_1539_split_member_separation.rs"]
mod issue_1539_split_member_separation;
