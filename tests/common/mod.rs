//! Shared test infrastructure for ferro-hgvs integration tests.
//!
//! Exposed sub-modules:
//! - `synthetic`: builders for synthetic `MockProvider` fixtures across
//!   coordinate systems, used by the ins/del/dup shift coverage matrices.
//! - `failure_expectations`: per-input expectations framework for the
//!   bulk parser fixtures (cmrg, paraphase, clinvar 500K + unique). See
//!   `failure_expectations.rs` for the snapshot shape and contract;
//!   tracking issue: #174.
//!
//! `#![allow(dead_code)]`: each integration test binary that does
//! `mod common;` compiles the full common tree into its binary and the
//! compiler reports unused items as dead. Different binaries use
//! different subsets, so the warnings are per-binary noise rather than
//! a real signal.
#![allow(dead_code)]

pub mod failure_expectations;
pub mod synthetic;
