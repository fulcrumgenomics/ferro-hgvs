//! Deserialization schema for the biocommons-normalize conformance corpus
//! (`tests/fixtures/biocommons-normalize/cases.json`).
//!
//! Single source of truth for the corpus schema, shared by the integration
//! harness (`tests/biocommons_normalize_tests.rs`) and the summary generator
//! (`examples/generate_conformance_summary.rs`); see the module docs in
//! [`super::mutalyzer`] for the rationale.

use serde::Deserialize;

use super::schema::{validate_cluster_refs, Cluster};

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
    /// Root-cause cluster taxonomy referenced by the per-disposition `cluster`
    /// field. Optional so a corpus with no taxonomy still parses.
    #[serde(default)]
    pub clusters: Vec<Cluster>,
    pub cases: Vec<Case>,
}

impl Fixture {
    /// Every `(input, cluster_id)` pair referenced by a disposition `cluster`.
    pub fn cluster_refs(&self) -> Vec<(&str, &str)> {
        self.cases
            .iter()
            .filter_map(|case| {
                let cluster = case.disposition.as_ref()?.cluster()?;
                Some((case.input.as_str(), cluster))
            })
            .collect()
    }

    /// Validate that every disposition `cluster` ref resolves to a registry
    /// entry (orphan clusters are allowed). See [`validate_cluster_refs`].
    pub fn validate_clusters(&self) -> Result<(), String> {
        validate_cluster_refs(&self.clusters, self.cluster_refs())
    }
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
        /// Root-cause cluster id (see the corpus `clusters` registry).
        #[serde(default)]
        cluster: Option<String>,
    },
    /// ferro is wrong; xfail until fixed. If ferro starts matching biocommons
    /// the harness fails (XPASS) so the annotation and the fixed row are
    /// cleaned up.
    KnownBug {
        /// Issue tracking the fix.
        tracking_issue: u64,
        /// The (incorrect) output ferro currently produces.
        ferro_output: String,
        /// Root-cause cluster id (see the corpus `clusters` registry).
        #[serde(default)]
        cluster: Option<String>,
    },
}

impl Disposition {
    /// The root-cause cluster id this disposition references, if any.
    pub fn cluster(&self) -> Option<&str> {
        match self {
            Disposition::AcceptedDivergence { cluster, .. }
            | Disposition::KnownBug { cluster, .. } => cluster.as_deref(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(clusters: &str, disposition: &str) -> Fixture {
        let json = format!(
            r#"{{"description":"t","source":"t","source_commit":"t",
                 "biocommons_upstream":"t","license":"t","refreshed_at":"t",
                 "clusters":[{clusters}],
                 "cases":[{{"input":"X","normalized":"Y","shuffle_direction":"3prime",
                   "cross_boundaries":false,"source_function":"t","disposition":{disposition}}}]}}"#
        );
        serde_json::from_str(&json).expect("fixture should deserialize")
    }

    #[test]
    fn cluster_ref_resolves_to_registry() {
        let clusters = r#"{"id":"panic","title":"overflow panic","spec_section":"x"}"#;
        let disp =
            r#"{"kind":"known_bug","tracking_issue":472,"ferro_output":"Z","cluster":"panic"}"#;
        let fixture = parse(clusters, disp);
        assert_eq!(fixture.cluster_refs(), vec![("X", "panic")]);
        assert!(fixture.validate_clusters().is_ok());
    }

    #[test]
    fn dangling_cluster_ref_is_rejected() {
        let disp =
            r#"{"kind":"known_bug","tracking_issue":472,"ferro_output":"Z","cluster":"missing"}"#;
        let err = parse("", disp)
            .validate_clusters()
            .expect_err("a dangling cluster ref must be rejected");
        assert!(err.contains("missing"), "{err}");
    }

    #[test]
    fn disposition_without_cluster_is_ok() {
        let disp = r#"{"kind":"known_bug","tracking_issue":472,"ferro_output":"Z"}"#;
        let fixture = parse("", disp);
        assert!(fixture.cluster_refs().is_empty());
        assert!(fixture.validate_clusters().is_ok());
    }
}
