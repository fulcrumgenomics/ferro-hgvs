//! Shared test infrastructure for ferro-hgvs integration tests.
//!
//! Exposed sub-modules:
//! - `synthetic`: builders for synthetic `MockProvider` fixtures across
//!   coordinate systems, used by the ins/del/dup shift coverage matrices.
//! - `failure_expectations`: per-input expectations framework for the
//!   bulk parser fixtures (cmrg, paraphase, clinvar 500K + unique). See
//!   `failure_expectations.rs` for the snapshot shape and contract;
//!   tracking issue: #174.
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

pub mod failure_expectations;
pub mod spec_enumeration;
pub mod spec_fixture;
pub mod synthetic;
