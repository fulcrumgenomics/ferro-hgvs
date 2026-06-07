//! Deserialization schema for the biocommons-normalize conformance corpus
//! (`tests/fixtures/biocommons-normalize/cases.json`).
//!
//! Single source of truth for the corpus schema, shared by the integration
//! harness (`tests/biocommons_normalize_tests.rs`) and the summary generator
//! (`examples/generate_conformance_summary.rs`); see the module docs in
//! [`super::mutalyzer`] for the rationale.

use serde::Deserialize;

/// Top-level corpus document: metadata header plus the case list.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub struct Fixture {
    pub description: String,
    pub source: String,
    pub source_commit: String,
    pub biocommons_upstream: String,
    pub license: String,
    pub refreshed_at: String,
    pub cases: Vec<Case>,
}

/// One corpus case: an input, the expected normalized form (or error), the
/// normalization configuration, and any divergence disposition.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub struct Case {
    pub input: String,
    /// `None` means biocommons expected `normalize()` to raise an error.
    /// Mutually informed by `expects_error: true`.
    pub normalized: Option<String>,
    #[serde(default)]
    pub expects_error: bool,
    pub shuffle_direction: String,
    pub cross_boundaries: bool,
    pub source_function: String,
    #[serde(default = "default_true")]
    pub to_test: bool,
    /// Optional disposition for a row where ferro diverges from biocommons.
    /// Absent means "ferro must match biocommons". See [`Disposition`] and
    /// `classify_outcome` (issue #478, conformance-harness redesign).
    #[serde(default)]
    pub disposition: Option<Disposition>,
}

fn default_true() -> bool {
    true
}

/// Why a row diverges from biocommons, and what ferro produces instead.
///
/// Drives the conformance harness's pass/fail decision (issue #478 pillars 1-2).
/// An annotated divergence whose recorded `ferro_output` still holds is
/// accounted for (not a failure). An annotation that no longer holds fails the
/// suite loudly so it gets corrected: a `known_bug` that started matching
/// biocommons (XPASS — the bug appears fixed), or any drift away from the
/// recorded `ferro_output`. This is what makes the burn-down self-correcting —
/// a fixed bug cannot linger silently in the fixture.
#[derive(Debug, Clone, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
#[allow(dead_code)] // `spec_citation` is documentation for humans reading cases.json.
pub enum Disposition {
    /// ferro intentionally diverges from biocommons and is spec-correct.
    AcceptedDivergence {
        /// Human rationale: the HGVS-spec basis or data-availability constraint.
        reason: String,
        /// Optional citation into the HGVS recommendations.
        #[serde(default)]
        spec_citation: Option<String>,
        /// The exact output ferro is expected to produce.
        ferro_output: String,
    },
    /// ferro is wrong; xfail until fixed. If ferro starts matching biocommons
    /// the harness fails (XPASS) so the annotation and the fixed row are
    /// cleaned up.
    KnownBug {
        /// Issue tracking the fix.
        tracking_issue: u64,
        /// The (incorrect) output ferro currently produces.
        ferro_output: String,
    },
}
