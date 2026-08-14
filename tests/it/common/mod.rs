//! Shared test infrastructure for ferro-hgvs integration tests.
//!
//! Exposed sub-modules:
//! - `bulk_fixtures`: the one place that decides what an ABSENT bulk corpus
//!   means. Those four fixtures are release assets rather than git objects, and
//!   every suite reading them skips green without them, so
//!   `FERRO_REQUIRE_BULK_FIXTURES` promotes the skip to a failure in CI.
//! - `cis_apply_oracle`: the SPDI-based apply oracle the sibling-crossing
//!   cis-allele tests use to check that normalization did not change the
//!   sequence a description denotes, plus `sweep_sequences`, the deterministic
//!   corpus their exhaustive sweeps enumerate over and pin counts against.
//! - `synthetic`: builders for synthetic `MockProvider` fixtures across
//!   coordinate systems, used by the ins/del/dup shift coverage matrices.
//! - `failure_expectations`: per-input expectations framework for the
//!   bulk parser fixtures (cmrg, paraphase, clinvar 500K + unique). See
//!   `failure_expectations.rs` for the snapshot shape and contract;
//!   tracking issue: #174.
//! - `manifest`: the sibling of `bulk_fixtures` for the other input that can
//!   go missing without anything going red — a ferro-prepared reference. Every
//!   reference-aware guard skips green without one, so
//!   `FERRO_REQUIRE_MANIFEST` promotes that skip to a failure in the one job
//!   that prepares a manifest.
//! - `fixture_gen`: the shared on-demand regeneration flow (locking,
//!   subprocess, atomic rename) that `spec_fixture` and `spec_enumeration`
//!   both wrap.
//! - `rulings`: reads the adjudication ledger's `rulings` section into typed
//!   records, so the citation-currency scan and the clause index share one
//!   definition of what a record is.
//! - `spec_fixture`: regenerates the gitignored HGVS spec-normalization
//!   fixture on demand so tests that read it work on a fresh checkout.
//! - `spec_enumeration`: same, for the gitignored exhaustive spec test
//!   enumeration (`hgvs_spec_enumeration.json`).
//!
//! `#![allow(dead_code)]`: each integration test binary that does
//! `mod common;` compiles the full common tree into its binary and the
//! compiler reports unused items as dead. Different binaries use
//! different subsets, so the warnings are per-binary noise rather than
//! a real signal.
#![allow(dead_code)]

pub mod bulk_fixtures;
pub mod cis_apply_oracle;
pub mod failure_expectations;
pub mod fixture_gen;
pub mod manifest;
pub mod rulings;
pub mod spec_enumeration;
pub mod spec_fixture;
pub mod synthetic;
