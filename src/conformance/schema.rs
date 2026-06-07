//! Cross-corpus conformance schema: the cluster taxonomy.
//!
//! A *cluster* is a named cross-case root-cause grouping (the "Pattern A/B/I"
//! taxonomy that the hand-maintained `failure-patterns.md` docs used to carry in
//! prose). Promoting it to data in `cases.json` lets the generated summary
//! derive the cross-case narrative instead of anyone hand-maintaining it (issue
//! #509).

use serde::Deserialize;

/// One root-cause cluster: a stable id plus the editorial title and spec
/// citation that the generated summary renders. Defined once in each corpus's
/// `cases.json` `clusters` registry and referenced by id from each divergence
/// disposition.
#[derive(Debug, Clone, Deserialize)]
#[allow(dead_code)]
pub struct Cluster {
    /// Stable slug referenced by dispositions (e.g. `"refseqgene-selector"`).
    pub id: String,
    /// Human title / root-cause description rendered in the summary.
    pub title: String,
    /// HGVS spec citation for the cluster (a recommendations path or section).
    pub spec_section: String,
    /// Optional editorial note.
    #[serde(default)]
    pub note: Option<String>,
}

/// Validate that every `cluster` reference resolves to a registry entry.
///
/// `refs` yields `(input, cluster_id)` pairs — the `cases.json` input is carried
/// for error context. A reference to an id absent from `clusters` is a hard
/// error (a typo or a deleted cluster). Unreferenced ("orphan") clusters are
/// **allowed**: they hold durable taxonomy whose member rows may not be seeded
/// yet (e.g. the biocommons corpus before its manifest-backed seeding), and the
/// generator renders them as taxonomy-only.
pub fn validate_cluster_refs<'a>(
    clusters: &[Cluster],
    refs: impl IntoIterator<Item = (&'a str, &'a str)>,
) -> Result<(), String> {
    let known: std::collections::HashSet<&str> = clusters.iter().map(|c| c.id.as_str()).collect();
    for (input, cluster_id) in refs {
        if !known.contains(cluster_id) {
            return Err(format!(
                "input {input:?} references unknown cluster id {cluster_id:?}; \
                 every `cluster` ref must resolve to a `clusters` registry entry"
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cluster(id: &str) -> Cluster {
        Cluster {
            id: id.to_string(),
            title: format!("title for {id}"),
            spec_section: "recommendations/DNA/example.md".to_string(),
            note: None,
        }
    }

    #[test]
    fn dangling_cluster_ref_is_rejected() {
        let clusters = vec![cluster("refseqgene-selector")];
        let refs = vec![
            ("NG_1(NM_1):c.1del", "refseqgene-selector"),
            ("NG_2(NM_2):c.2del", "typo-cluster"),
        ];
        let err = validate_cluster_refs(&clusters, refs)
            .expect_err("a reference to an unknown cluster id must be rejected");
        assert!(
            err.contains("typo-cluster"),
            "error should name the dangling id: {err}"
        );
        assert!(
            err.contains("NG_2(NM_2):c.2del"),
            "error should name the offending input for context: {err}"
        );
    }

    #[test]
    fn orphan_cluster_is_allowed() {
        let clusters = vec![cluster("used"), cluster("never-referenced")];
        let refs = vec![("NM_1:c.1del", "used")];
        assert!(
            validate_cluster_refs(&clusters, refs).is_ok(),
            "an unreferenced cluster is durable taxonomy, not an error"
        );
    }

    #[test]
    fn no_references_is_ok() {
        let clusters: Vec<Cluster> = vec![];
        let refs: Vec<(&str, &str)> = vec![];
        assert!(validate_cluster_refs(&clusters, refs).is_ok());
    }
}
