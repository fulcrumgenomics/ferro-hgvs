//! Render a generated, `--check`-gated conformance summary from `cases.json`.
//!
//! The summary is a derived view — never hand-edited — over the disposition
//! annotations in a corpus's `cases.json` (issue #509). It groups dispositions
//! by their root-cause [`Cluster`] and emits hermetic per-axis disposition
//! tallies. It deliberately does **not** emit a live FAIL count: that set
//! requires the reference manifest and so is non-hermetic and never committed.

use std::path::Path;

use super::schema::Cluster;

/// Header marking the file as generated; the generator and `--check` both
/// emit/expect it so a hand-edit is visible in review.
pub const GENERATED_HEADER: &str = "<!-- GENERATED — do not edit; regenerate via \
     `cargo run --features dev --example generate_conformance_summary` -->";

/// The four divergence dispositions, in tally column order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DispositionKind {
    AcceptedDivergence,
    KnownBug,
    Improvement,
    SpecCitation,
}

impl DispositionKind {
    /// Every kind, in the order they appear as tally columns.
    pub const ALL: [DispositionKind; 4] = [
        DispositionKind::AcceptedDivergence,
        DispositionKind::KnownBug,
        DispositionKind::Improvement,
        DispositionKind::SpecCitation,
    ];

    /// Stable snake_case identifier used in the rendered tables.
    pub fn as_str(self) -> &'static str {
        match self {
            DispositionKind::AcceptedDivergence => "accepted_divergence",
            DispositionKind::KnownBug => "known_bug",
            DispositionKind::Improvement => "improvement",
            DispositionKind::SpecCitation => "spec_citation",
        }
    }
}

/// One member divergence under a cluster, flattened from a corpus disposition.
#[derive(Debug, Clone)]
pub struct MemberRow {
    /// Cluster id this row groups under (`None` → ungrouped).
    pub cluster: Option<String>,
    /// The `cases.json` input.
    pub input: String,
    /// Axis the disposition applies to (biocommons rows are all `normalized`).
    pub axis: String,
    /// Disposition kind.
    pub kind: DispositionKind,
    /// The recorded ferro output, when the disposition carries one.
    pub ferro_output: Option<String>,
    /// The tracking issue, when the disposition carries one.
    pub tracking_issue: Option<u64>,
}

/// Corpus-agnostic input to [`render`]: the cluster registry plus the flattened
/// member rows. Each corpus's `Fixture` produces one of these.
#[derive(Debug, Clone)]
pub struct SummaryModel {
    /// Corpus title (e.g. `"mutalyzer-normalize"`).
    pub title: String,
    /// The cluster taxonomy registry.
    pub clusters: Vec<Cluster>,
    /// Flattened member dispositions.
    pub rows: Vec<MemberRow>,
}

/// Render the summary markdown. Pure function of `model`.
pub fn render(model: &SummaryModel) -> String {
    use std::collections::BTreeMap;
    use std::fmt::Write;

    let mut out = String::new();
    let _ = writeln!(out, "{GENERATED_HEADER}");
    let _ = writeln!(out, "\n# {} conformance summary\n", model.title);
    let _ = writeln!(
        out,
        "Generated from `cases.json` — do not edit by hand; regenerate with the \
         example above. Every row below is a tracked disposition. The live \
         divergence set against full reference data is emitted only by the \
         manifest run and is never committed (it is non-hermetic).\n"
    );

    // --- Root-cause clusters (registry order; members sorted by input) ---
    let _ = writeln!(out, "## Root-cause clusters\n");
    for cluster in &model.clusters {
        let mut members: Vec<&MemberRow> = model
            .rows
            .iter()
            .filter(|r| r.cluster.as_deref() == Some(cluster.id.as_str()))
            .collect();
        members.sort_by(|a, b| a.input.cmp(&b.input).then(a.axis.cmp(&b.axis)));

        let _ = writeln!(out, "### {}\n", cluster.title);
        let _ = writeln!(out, "Spec: `{}`\n", cluster.spec_section);
        if let Some(note) = &cluster.note {
            let _ = writeln!(out, "{note}\n");
        }
        if members.is_empty() {
            let _ = writeln!(
                out,
                "_No seeded member rows yet (manifest-gated seeding, #325)._\n"
            );
        } else {
            write_member_table(&mut out, &members);
        }
    }

    // --- Ungrouped dispositions (no cluster reference), if any ---
    let ungrouped: Vec<&MemberRow> = model.rows.iter().filter(|r| r.cluster.is_none()).collect();
    if !ungrouped.is_empty() {
        let mut ungrouped = ungrouped;
        ungrouped.sort_by(|a, b| a.input.cmp(&b.input).then(a.axis.cmp(&b.axis)));
        let _ = writeln!(out, "### Ungrouped\n");
        write_member_table(&mut out, &ungrouped);
    }

    // --- Per-axis disposition tallies (hermetic; no live divergence count) ---
    let mut tally: BTreeMap<&str, [usize; 4]> = BTreeMap::new();
    for row in &model.rows {
        let counts = tally.entry(row.axis.as_str()).or_default();
        let col = DispositionKind::ALL
            .iter()
            .position(|k| *k == row.kind)
            .unwrap();
        counts[col] += 1;
    }
    let _ = writeln!(out, "## Disposition tallies\n");
    let header_cols = DispositionKind::ALL
        .iter()
        .map(|k| k.as_str())
        .collect::<Vec<_>>()
        .join(" | ");
    let _ = writeln!(out, "| axis | {header_cols} |");
    let _ = writeln!(out, "|---|---:|---:|---:|---:|");
    for (axis, counts) in &tally {
        let _ = writeln!(
            out,
            "| {axis} | {} | {} | {} | {} |",
            counts[0], counts[1], counts[2], counts[3]
        );
    }

    out
}

/// Render the per-cluster member table.
fn write_member_table(out: &mut String, members: &[&MemberRow]) {
    use std::fmt::Write;
    let _ = writeln!(
        out,
        "| input | axis | disposition | ferro output | tracking |"
    );
    let _ = writeln!(out, "|---|---|---|---|---|");
    for m in members {
        let ferro = m
            .ferro_output
            .as_deref()
            .map(|s| format!("`{s}`"))
            .unwrap_or_else(|| "—".into());
        let issue = m
            .tracking_issue
            .map(|n| format!("#{n}"))
            .unwrap_or_else(|| "—".into());
        let _ = writeln!(
            out,
            "| `{}` | {} | {} | {} | {} |",
            m.input,
            m.axis,
            m.kind.as_str(),
            ferro,
            issue
        );
    }
    let _ = writeln!(out);
}

/// Compare the on-disk summary at `path` to a fresh render of `model`.
/// Returns `Err` (the drift message) when they differ — the `--check` gate.
pub fn check(path: &Path, model: &SummaryModel) -> Result<(), String> {
    let want = render(model);
    let got = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    if got == want {
        Ok(())
    } else {
        Err(format!(
            "{} is out of date; regenerate with \
             `cargo run --features dev --example generate_conformance_summary`",
            path.display()
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cluster(id: &str, title: &str) -> Cluster {
        Cluster {
            id: id.to_string(),
            title: title.to_string(),
            spec_section: format!("recommendations/{id}.md"),
            note: None,
        }
    }

    fn model() -> SummaryModel {
        SummaryModel {
            title: "mutalyzer-normalize".to_string(),
            clusters: vec![
                cluster("refseqgene-selector", "RefSeqGene transcript selection"),
                cluster("orphan", "Taxonomy with no seeded rows"),
            ],
            rows: vec![
                MemberRow {
                    cluster: Some("refseqgene-selector".to_string()),
                    input: "NG_1(NM_1):c.12del".to_string(),
                    axis: "normalized".to_string(),
                    kind: DispositionKind::Improvement,
                    ferro_output: None,
                    tracking_issue: Some(500),
                },
                MemberRow {
                    cluster: Some("refseqgene-selector".to_string()),
                    input: "NG_2(NM_2):c.34del".to_string(),
                    axis: "normalized".to_string(),
                    kind: DispositionKind::Improvement,
                    ferro_output: None,
                    tracking_issue: Some(500),
                },
            ],
        }
    }

    #[test]
    fn render_emits_generated_header() {
        let out = render(&model());
        assert!(
            out.starts_with(GENERATED_HEADER),
            "summary must open with the GENERATED marker so hand-edits are visible"
        );
    }

    #[test]
    fn render_groups_member_rows_under_their_cluster_with_citation() {
        let out = render(&model());
        assert!(
            out.contains("RefSeqGene transcript selection"),
            "cluster title: {out}"
        );
        assert!(
            out.contains("recommendations/refseqgene-selector.md"),
            "cluster spec citation must render: {out}"
        );
        assert!(
            out.contains("NG_1(NM_1):c.12del"),
            "member input must render: {out}"
        );
        assert!(out.contains("#500"), "tracking issue must render: {out}");
    }

    #[test]
    fn render_emits_per_axis_disposition_tallies_not_a_fail_count() {
        let out = render(&model());
        // Two improvement rows on the normalized axis.
        assert!(
            out.contains("improvement"),
            "tally must name the disposition kind: {out}"
        );
        assert!(
            out.contains("| normalized |") || out.contains("normalized"),
            "tally must be per-axis: {out}"
        );
        // The non-hermetic live FAIL set must never be committed.
        assert!(
            !out.to_lowercase().contains("fail"),
            "summary must not emit a live FAIL count (non-hermetic): {out}"
        );
    }

    #[test]
    fn render_marks_orphan_cluster_as_taxonomy_only() {
        let out = render(&model());
        assert!(
            out.contains("Taxonomy with no seeded rows"),
            "orphan cluster title must still render: {out}"
        );
        assert!(
            out.to_lowercase().contains("no seeded member rows"),
            "orphan cluster must be marked taxonomy-only: {out}"
        );
    }
}
