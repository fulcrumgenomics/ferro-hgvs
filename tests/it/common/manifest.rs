//! The integration tree's handle on the one place that decides what an ABSENT
//! reference manifest means.
//!
//! The policy, the `FERRO_REQUIRE_MANIFEST` predicate and the failure text all
//! live in [`ferro_hgvs::conformance::manifest_gate`], which is where their
//! reasoning is written down. They are in the library rather than here because
//! two of the manifest-gated guards are `#[cfg(test)]` units inside `src/`
//! itself and cannot see this tree — and a second copy for those two is the
//! drift the whole mechanism exists to prevent.
//!
//! This module is the re-export, so integration modules can keep writing
//! `crate::common::manifest::absent(..)` alongside
//! `crate::common::bulk_fixtures::absent(..)` and read as one convention.

// `allow(unused_imports)`: the same reason `mod.rs` carries `allow(dead_code)` —
// this tree is compiled into one binary and only `absent` has call sites today.
// The other three are the module's stated surface, and a re-export that
// disappears whenever its last caller does is one nobody can find again.
#[allow(unused_imports)]
pub use ferro_hgvs::conformance::manifest_gate::{
    absent, manifest_is_required, unserved, MANIFEST_ENV, REQUIRE_ENV,
};
