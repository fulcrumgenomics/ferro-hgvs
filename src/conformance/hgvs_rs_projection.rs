//! Deserialization schema for the hgvs-rs projection conformance corpus
//! (`tests/fixtures/hgvs-rs-projection/cases.json`).
//!
//! This is the **third consumer** of the corpus harness, mirroring
//! [`super::mutalyzer`] (the multi-axis, projection-focused template) rather
//! than the single-axis biocommons normalize corpus. The case `input` is a
//! genome-anchored HGVSg (or a bare HGVSc for the direct `c→p` rows); the
//! expected HGVSc/HGVSp/`n.` values drive the `coding` / `protein_description`
//! / `coding_protein_descriptions` / `noncoding` projection axes.
//!
//! These types are the single source of truth for the corpus schema: both the
//! integration harness (`tests/hgvs_rs_projection_tests.rs`) and the summary
//! generator (`examples/generate_conformance_summary.rs`) deserialize
//! `cases.json` through them, so the generated `failure-patterns.md` cannot
//! drift from the schema the harness enforces. They live in the library only so
//! the `examples/` binary can share them.
//!
//! The disposition machinery (`Axis` / `Policy` / `SpecSection` enums and the
//! `AcceptedDivergence` / `KnownBug` / `Improvement` / `SpecCitation` structs)
//! is reused verbatim from [`super::mutalyzer`] — a divergence on a projection
//! axis is triaged exactly the same way as a divergence on a normalize axis.

use serde::Deserialize;

// Reuse the mutalyzer disposition machinery wholesale — these are corpus-
// agnostic and already public.
pub use super::mutalyzer::{
    AcceptedDivergence, Axis, Improvement, KnownBug, Policy, SpecCitation, SpecSection,
};
use super::schema::{validate_cluster_refs, Cluster};
use super::summary::{DispositionKind, MemberRow, SummaryModel};

/// Corpus title used as the generated summary's heading and model title.
pub const CORPUS_TITLE: &str = "hgvs-rs-projection";

/// Top-level corpus document: metadata header plus the case list.
#[derive(Debug, Deserialize)]
#[allow(dead_code)]
pub struct Fixture {
    pub description: String,
    pub source: String,
    pub source_commit: String,
    pub license: String,
    pub refreshed_at: String,
    /// Root-cause cluster taxonomy referenced by the per-disposition `cluster`
    /// field. Optional so a corpus with no taxonomy still parses (this corpus
    /// seeds it empty until an initial manifest run reveals divergences).
    #[serde(default)]
    pub clusters: Vec<Cluster>,
    pub cases: Vec<Case>,
}

impl Fixture {
    /// Every `(input, cluster_id)` pair referenced by a disposition `cluster`.
    pub fn cluster_refs(&self) -> Vec<(&str, &str)> {
        let mut refs = Vec::new();
        for case in &self.cases {
            let input = case.input.as_str();
            for cluster in [
                case.accepted_divergence
                    .as_ref()
                    .and_then(|d| d.cluster.as_deref()),
                case.known_bug.as_ref().and_then(|d| d.cluster.as_deref()),
                case.improvement.as_ref().and_then(|d| d.cluster.as_deref()),
                case.spec_citation
                    .as_ref()
                    .and_then(|d| d.cluster.as_deref()),
            ]
            .into_iter()
            .flatten()
            {
                refs.push((input, cluster));
            }
        }
        refs
    }

    /// Validate that every disposition `cluster` ref resolves to a registry
    /// entry (orphan clusters are allowed). See [`validate_cluster_refs`].
    pub fn validate_clusters(&self) -> Result<(), String> {
        validate_cluster_refs(&self.clusters, self.cluster_refs())
    }

    /// Flatten this corpus into the corpus-agnostic summary model. Dispositions
    /// carry no `ferro_output`, so that column is always absent.
    pub fn to_summary(&self) -> SummaryModel {
        let mut rows = Vec::new();
        for case in &self.cases {
            let input = case.input.as_str();
            if let Some(d) = &case.accepted_divergence {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::AcceptedDivergence,
                    ferro_output: None,
                    tracking_issue: None,
                });
            }
            if let Some(d) = &case.known_bug {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::KnownBug,
                    ferro_output: None,
                    tracking_issue: Some(d.tracking_issue),
                });
            }
            if let Some(d) = &case.improvement {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::Improvement,
                    ferro_output: None,
                    tracking_issue: Some(d.tracking_issue),
                });
            }
            if let Some(d) = &case.spec_citation {
                rows.push(MemberRow {
                    cluster: d.cluster.clone(),
                    input: input.to_string(),
                    axis: d.axis.as_str().to_string(),
                    kind: DispositionKind::SpecCitation,
                    ferro_output: None,
                    tracking_issue: None,
                });
            }
        }
        SummaryModel {
            title: CORPUS_TITLE.to_string(),
            clusters: self.clusters.clone(),
            rows,
        }
    }
}

/// One corpus case: an input plus the expected output on each projection axis
/// and any per-axis divergence disposition.
///
/// Unlike the mutalyzer `Case`, this carries only the fields this corpus uses
/// (`coding`, `protein_description`, `coding_protein_descriptions`,
/// `noncoding`) plus the **new** `coding` axis (`g→c`-only UTR rows). The
/// disposition annotations are the same closed-enum-backed structs the
/// mutalyzer corpus uses.
#[derive(Debug, Deserialize, Clone)]
#[allow(dead_code)]
pub struct Case {
    #[serde(default)]
    pub keywords: Vec<String>,
    pub input: String,
    /// Expected coding (`c.`) projection for `g→c`-only UTR rows (a row with an
    /// HGVSg + HGVSc but no HGVSp). Drives the `coding` axis. **New** field
    /// relative to the mutalyzer schema.
    #[serde(default)]
    pub coding: Option<String>,
    /// Expected protein (`p.`) projection for direct `c→p` rows (a row with an
    /// HGVSc + HGVSp but no HGVSg). Drives the `protein_description` axis.
    #[serde(default)]
    pub protein_description: Option<String>,
    /// Expected `(c., p.)` pairs for the rich `g→c→p` rows. Drives the
    /// `coding_protein_descriptions` axis.
    #[serde(default)]
    pub coding_protein_descriptions: Option<Vec<Vec<String>>>,
    /// Expected non-coding (`n.`) projections for NR_ transcript rows. Drives
    /// the `noncoding` axis.
    #[serde(default)]
    pub noncoding: Option<Vec<String>>,
    #[serde(default = "default_true")]
    pub to_test: bool,
    /// Marks a mismatch on a specific axis as an accepted (non-bug) divergence —
    /// e.g. ferro's spec-preferred protein glyph form vs hgvs-rs's `ext*`. See
    /// [`AcceptedDivergence`].
    #[serde(default)]
    pub accepted_divergence: Option<AcceptedDivergence>,
    /// Marks a mismatch on a specific axis as a known ferro bug (xfail). See
    /// [`KnownBug`].
    #[serde(default)]
    pub known_bug: Option<KnownBug>,
    /// Marks a mismatch on a specific axis as a tracked improvement (ferro's
    /// output is valid but not spec-preferred). See [`Improvement`].
    #[serde(default)]
    pub improvement: Option<Improvement>,
    /// Documentary spec citation when `cases.json`'s expected value has been
    /// corrected to ferro's spec-correct value. See [`SpecCitation`].
    #[serde(default)]
    pub spec_citation: Option<SpecCitation>,
}

fn default_true() -> bool {
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    const BASE: &str =
        r#""description":"t","source":"t","source_commit":"t","license":"t","refreshed_at":"t""#;

    fn parse(clusters: &str, cases: &str) -> Fixture {
        let json = format!("{{{BASE},\"clusters\":[{clusters}],\"cases\":[{cases}]}}");
        serde_json::from_str(&json).expect("fixture should deserialize")
    }

    #[test]
    fn case_parses_each_projection_shape() {
        let cases = r#"
            {"input":"NC_000002.11:g.96780553C>T",
              "coding_protein_descriptions":[["NM_000682.5:c.1345G>A","NP_000673.2:p.Ala449Thr"]]},
            {"input":"NC_000002.11:g.96778665_96778666insA",
              "coding":"NM_000682.5:c.*1879_*1880insT"},
            {"input":"NM_000051.3:c.9170_9171delGA",
              "protein_description":"NP_000042.3:p.Ter3057Pheext*4"},
            {"input":"NC_000001.10:g.12776161G>A",
              "noncoding":["NR_111984.1:n.44G>A"]},
            {"input":"x","to_test":false}
        "#;
        let fixture = parse("", cases);
        assert_eq!(fixture.cases.len(), 5);
        assert!(fixture.cases[0].coding_protein_descriptions.is_some());
        assert_eq!(
            fixture.cases[1].coding.as_deref(),
            Some("NM_000682.5:c.*1879_*1880insT")
        );
        assert_eq!(
            fixture.cases[2].protein_description.as_deref(),
            Some("NP_000042.3:p.Ter3057Pheext*4")
        );
        assert!(fixture.cases[3].noncoding.is_some());
        assert!(fixture.cases[0].to_test);
        assert!(!fixture.cases[4].to_test);
    }

    #[test]
    fn coding_axis_disposition_roundtrips() {
        // The new `Axis::Coding` variant deserializes and renders.
        let cases = r#"
            {"input":"g","coding":"c",
              "accepted_divergence":{"axis":"coding",
                "policy":"ferro-policy-121-gene-symbol-selector","cluster":"sel"}}
        "#;
        let clusters = r#"{"id":"sel","title":"selector","spec_section":"refseq.md"}"#;
        let fixture = parse(clusters, cases);
        let ad = fixture.cases[0]
            .accepted_divergence
            .as_ref()
            .expect("accepted_divergence present");
        assert_eq!(ad.axis, Axis::Coding);
        assert_eq!(ad.axis.as_str(), "coding");
        assert!(fixture.validate_clusters().is_ok());

        let summary = fixture.to_summary();
        assert_eq!(summary.title, "hgvs-rs-projection");
        assert_eq!(summary.rows.len(), 1);
        assert_eq!(summary.rows[0].axis, "coding");
    }

    #[test]
    fn dangling_disposition_cluster_ref_is_rejected() {
        let cases = r#"
            {"input":"g","coding":"c",
              "spec_citation":{"axis":"coding",
                "section":"HGVS protein reference (bare NP)","cluster":"missing"}}
        "#;
        let err = parse("", cases)
            .validate_clusters()
            .expect_err("a dangling disposition cluster ref must be rejected");
        assert!(err.contains("missing"), "{err}");
    }

    /// The committed corpus file parses and has at least one `to_test` case.
    #[test]
    fn loads_cases_json() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/hgvs-rs-projection/cases.json"
        );
        let content = fs::read_to_string(path).expect("read cases.json");
        let fixture: Fixture = serde_json::from_str(&content).expect("parse cases.json");
        assert!(!fixture.cases.is_empty(), "cases.json should have cases");
        assert!(
            fixture.cases.iter().any(|c| c.to_test),
            "cases.json should have at least one to_test case"
        );
        assert_eq!(
            fixture.source_commit, "cf0e5a8d1b9a2ae3e883a58826f397aa6cb52e23",
            "committed cases.json must be pinned to the vendored hgvs-rs SHA"
        );
    }
}
