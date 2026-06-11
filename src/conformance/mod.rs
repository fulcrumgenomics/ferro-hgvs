//! Conformance-corpus schema and summary generation.
//!
//! The `mutalyzer` and `biocommons` submodules hold the deserialization schema
//! for `tests/fixtures/{mutalyzer,biocommons}-normalize/cases.json`. They are
//! the single source of truth for those schemas: the integration harnesses and
//! the `examples/generate_conformance_summary.rs` generator both deserialize
//! through them, so the generated `failure-patterns.md` views cannot drift from
//! the schema the harnesses enforce (issue #509).
//!
//! These types live in the library — rather than the test crates — only so the
//! `examples/` generator can share them; they are not part of the runtime API.

pub mod biocommons;
pub mod hgvs_rs_projection;
pub mod mutalyzer;
pub mod schema;
pub mod summary;

pub use schema::{validate_cluster_refs, Cluster};
