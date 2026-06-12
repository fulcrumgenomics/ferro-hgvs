//! Tool-support matrix: types + table projection.
//!
//! Source of truth is `docs/tool_support_matrix.json`. This module loads it
//! and renders deterministic projections (markdown tables, website JSON). See
//! `docs/superpowers/specs/2026-06-09-tool-support-matrix-design.md`.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::json;

/// Whether a tool parses/validates a pattern.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Validate {
    /// Parses and performs semantic validation.
    Yes,
    /// Parses/round-trips but does no semantic validation (permissive stub).
    Permissive,
    /// Does not accept the pattern.
    No,
}

/// Whether/how a tool normalizes a pattern.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Normalize {
    /// Genuinely canonicalizes / 3'-shifts.
    Full,
    /// Works but requires network access (cannot be cached).
    NetworkOnly,
    /// Deliberately unsupported (explicit design decision).
    UnsupportedByDesign,
    /// Parses but normalization errors out.
    Errors,
    /// Parses but normalization is a no-op by design.
    No,
}

/// One (tool × pattern) support claim.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cell {
    pub validate: Validate,
    pub normalize: Normalize,
    /// Key into [`Matrix::footnotes`].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub footnote: Option<String>,
    /// Free-text rationale / provenance note (not rendered into tables).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

impl Cell {
    /// True when the tool accepts the pattern (renders the `V` component).
    pub fn validate_active(&self) -> bool {
        !matches!(self.validate, Validate::No)
    }

    /// True when the tool produces a normalized result (renders the `N`
    /// component). `network_only` counts as active.
    pub fn normalize_active(&self) -> bool {
        matches!(self.normalize, Normalize::Full | Normalize::NetworkOnly)
    }

    /// Website glyph: `V/N`, `V`, `N`, or `-`.
    pub fn website_glyph(&self) -> &'static str {
        match (self.validate_active(), self.normalize_active()) {
            (true, true) => "V/N",
            (true, false) => "V",
            (false, true) => "N",
            (false, false) => "-",
        }
    }

    /// Markdown glyph for the normalization-capabilities view: `✓`, `Net`, `✗`.
    pub fn markdown_glyph(&self) -> &'static str {
        match self.normalize {
            Normalize::Full => "✓",
            Normalize::NetworkOnly => "Net",
            Normalize::UnsupportedByDesign | Normalize::Errors | Normalize::No => "✗",
        }
    }
}

/// A pattern row within a category.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Row {
    pub key: String,
    pub label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub example: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spec_url: Option<String>,
    /// Keyed by tool id; must contain exactly the ids in [`Matrix::tools`].
    pub support: BTreeMap<String, Cell>,
    pub last_reviewed: String,
}

/// A group of rows (reference types, coordinate types, variant types).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Category {
    pub id: String,
    pub title: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub spec_url: Option<String>,
    pub rows: Vec<Row>,
}

/// One compared tool.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tool {
    pub id: String,
    pub name: String,
    pub version: String,
    pub url: String,
}

/// One selected row within a view (with an optional relabel).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewSelect {
    pub category: String,
    pub key: String,
    pub label: String,
}

/// A markdown projection (e.g. README "Normalization Capabilities").
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct View {
    pub id: String,
    pub title: String,
    pub row_label_header: String,
    pub select: Vec<ViewSelect>,
}

/// The whole matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Matrix {
    pub schema_version: u32,
    pub tools: Vec<Tool>,
    pub footnotes: BTreeMap<String, String>,
    pub categories: Vec<Category>,
    #[serde(default)]
    pub views: Vec<View>,
}

impl Matrix {
    /// Parse from a JSON string.
    pub fn from_json_str(s: &str) -> Result<Self, String> {
        serde_json::from_str(s).map_err(|e| format!("parse tool_support_matrix: {e}"))
    }

    /// Load from a file path.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let p = path.as_ref();
        let s = std::fs::read_to_string(p).map_err(|e| format!("read {}: {e}", p.display()))?;
        Self::from_json_str(&s)
    }

    /// Tool ids in declared order.
    pub fn tool_ids(&self) -> Vec<&str> {
        self.tools.iter().map(|t| t.id.as_str()).collect()
    }

    /// Fail if any row is missing a tool, references an unknown tool, or cites
    /// an undefined footnote; or a view selects a missing (category, key).
    pub fn validate_invariants(&self) -> Result<(), String> {
        let ids: Vec<&str> = self.tool_ids();
        for cat in &self.categories {
            for row in &cat.rows {
                for id in &ids {
                    if !row.support.contains_key(*id) {
                        return Err(format!(
                            "category {} row {} missing tool {id}",
                            cat.id, row.key
                        ));
                    }
                }
                for (id, cell) in &row.support {
                    if !ids.contains(&id.as_str()) {
                        return Err(format!(
                            "category {} row {} has unknown tool {id}",
                            cat.id, row.key
                        ));
                    }
                    if let Some(f) = &cell.footnote {
                        if !self.footnotes.contains_key(f) {
                            return Err(format!(
                                "category {} row {} tool {id} cites unknown footnote {f}",
                                cat.id, row.key
                            ));
                        }
                    }
                }
            }
        }
        for view in &self.views {
            for sel in &view.select {
                if self.find_row(&sel.category, &sel.key).is_none() {
                    return Err(format!(
                        "view {} selects missing row {}/{}",
                        view.id, sel.category, sel.key
                    ));
                }
            }
        }
        Ok(())
    }

    /// Find a row by (category id, row key).
    pub fn find_row(&self, category: &str, key: &str) -> Option<&Row> {
        self.categories
            .iter()
            .find(|c| c.id == category)?
            .rows
            .iter()
            .find(|r| r.key == key)
    }

    /// Render a markdown view (selected rows, relabeled) into a GitHub-flavored
    /// table plus a footnote legend. Column order follows [`Matrix::tools`].
    ///
    /// Returns `Err` if the view is not found, a row is missing a tool entry,
    /// or a cell cites an unknown footnote. Call
    /// [`Matrix::validate_invariants`] first to catch those conditions early.
    pub fn render_markdown_view(&self, view_id: &str) -> Result<String, String> {
        let view = self
            .views
            .iter()
            .find(|v| v.id == view_id)
            .ok_or_else(|| format!("unknown view {view_id}"))?;

        let mut fa = FootnoteAssigner::new();
        let mut out = String::new();

        // Header row.
        out.push_str(&format!("| {} |", view.row_label_header));
        for t in &self.tools {
            out.push_str(&format!(" {} |", t.name));
        }
        out.push('\n');

        // Alignment row: left-align the label column, center the tool columns,
        // widths matching the header cell text for stable diffs.
        out.push_str(&format!(
            "|{}|",
            "-".repeat(view.row_label_header.len() + 2)
        ));
        for t in &self.tools {
            out.push_str(&format!(":{}:|", "-".repeat(t.name.len())));
        }
        out.push('\n');

        // Body rows.
        for sel in &view.select {
            let row = self.find_row(&sel.category, &sel.key).ok_or_else(|| {
                format!("view {view_id} row {}/{} missing", sel.category, sel.key)
            })?;
            out.push_str(&format!("| {} |", sel.label));
            for t in &self.tools {
                let cell = row.support.get(&t.id).ok_or_else(|| {
                    format!("category/view row {} missing tool {}", row.key, t.id)
                })?;
                let mut g = cell.markdown_glyph().to_string();
                if let Some(f) = &cell.footnote {
                    g.push_str(&fa.marker(f));
                }
                out.push_str(&format!(" {g} |"));
            }
            out.push('\n');
        }

        // Footnote legend (escaped `*` so markdown doesn't italicize).
        let used = fa.used();
        if !used.is_empty() {
            out.push('\n');
            for (marker, key) in used {
                let text = self
                    .footnotes
                    .get(key)
                    .ok_or_else(|| format!("view {view_id} cites unknown footnote {key}"))?;
                let escaped = marker.replace('*', "\\*");
                out.push_str(&format!("{escaped} {text}\n"));
            }
        }

        Ok(out)
    }

    /// Build the render-ready website projection consumed by `main.js`.
    /// Glyphs and footnote markers are pre-computed; JS only builds the DOM.
    ///
    /// # Panics
    ///
    /// Panics if the matrix has not passed [`Matrix::validate_invariants`]
    /// (a row missing a tool, or a cell citing an unknown footnote).
    pub fn render_website_json(&self) -> serde_json::Value {
        let tools: Vec<serde_json::Value> = self
            .tools
            .iter()
            .map(|t| json!({"id": t.id, "name": t.name, "version": t.version, "url": t.url}))
            .collect();

        let categories: Vec<serde_json::Value> = self
            .categories
            .iter()
            .map(|cat| {
                let mut fa = FootnoteAssigner::new();
                let rows: Vec<serde_json::Value> = cat
                    .rows
                    .iter()
                    .map(|row| {
                        let cells: Vec<serde_json::Value> = self
                            .tools
                            .iter()
                            .map(|t| {
                                let cell = row.support.get(&t.id).expect(
                                    "row missing a tool; call Matrix::validate_invariants() before rendering",
                                );
                                let marker = cell
                                    .footnote
                                    .as_ref()
                                    .map(|f| fa.marker(f))
                                    .unwrap_or_default();
                                json!({"glyph": cell.website_glyph(), "footnote": marker})
                            })
                            .collect();
                        json!({
                            "key": row.key,
                            "label": row.label,
                            "example": row.example,
                            "spec_url": row.spec_url,
                            "cells": cells,
                        })
                    })
                    .collect();
                let legend: Vec<serde_json::Value> = fa
                    .used()
                    .into_iter()
                    .map(|(marker, key)| {
                        json!({"marker": marker, "text": self.footnotes.get(key).cloned().expect("cell cites unknown footnote; call Matrix::validate_invariants() before rendering")})
                    })
                    .collect();
                json!({"id": cat.id, "title": cat.title, "spec_url": cat.spec_url, "rows": rows, "footnotes": legend})
            })
            .collect();

        json!({
            "schema_version": self.schema_version,
            "tools": tools,
            "categories": categories,
        })
    }

    /// Return human-readable identifiers of cells whose row `last_reviewed` is
    /// older than `months` before `today` (both `YYYY-MM-DD`). Non-failing.
    pub fn stale_cells(&self, today: &str, months: i64) -> Vec<String> {
        let cutoff = months_before(today, months);
        let mut out = Vec::new();
        for cat in &self.categories {
            for row in &cat.rows {
                if row.last_reviewed.as_str() < cutoff.as_str() {
                    out.push(format!(
                        "{}/{} (last_reviewed {})",
                        cat.id, row.key, row.last_reviewed
                    ));
                }
            }
        }
        out
    }
}

/// Assigns `*`, `**`, `***`, … to footnote keys in first-appearance order
/// within one rendered table.
struct FootnoteAssigner {
    order: Vec<String>,
}

impl FootnoteAssigner {
    fn new() -> Self {
        Self { order: Vec::new() }
    }

    /// Return the marker for a key, registering it on first use.
    fn marker(&mut self, key: &str) -> String {
        let idx = match self.order.iter().position(|k| k == key) {
            Some(i) => i,
            None => {
                self.order.push(key.to_string());
                self.order.len() - 1
            }
        };
        "*".repeat(idx + 1)
    }

    /// Markers used so far, in assignment order, as (marker, key).
    fn used(&self) -> Vec<(String, &str)> {
        self.order
            .iter()
            .enumerate()
            .map(|(i, k)| ("*".repeat(i + 1), k.as_str()))
            .collect()
    }
}

/// Subtract `months` from a `YYYY-MM-DD` date lexically (string comparison is
/// valid for ISO dates). Clamps the day implicitly by keeping it; only the
/// year/month are shifted, which is sufficient for a coarse staleness check.
fn months_before(today: &str, months: i64) -> String {
    let mut parts = today.splitn(3, '-');
    let (y, m, d) = (
        parts
            .next()
            .and_then(|s| s.parse::<i64>().ok())
            .unwrap_or(0),
        parts
            .next()
            .and_then(|s| s.parse::<i64>().ok())
            .unwrap_or(1),
        parts.next().unwrap_or("01"),
    );
    let total = y * 12 + (m - 1) - months;
    let ny = total.div_euclid(12);
    let nm = total.rem_euclid(12) + 1;
    format!("{ny:04}-{nm:02}-{d}")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(v: Validate, n: Normalize) -> Cell {
        Cell {
            validate: v,
            normalize: n,
            footnote: None,
            note: None,
        }
    }

    #[test]
    fn website_glyphs() {
        assert_eq!(cell(Validate::Yes, Normalize::Full).website_glyph(), "V/N");
        assert_eq!(
            cell(Validate::Yes, Normalize::NetworkOnly).website_glyph(),
            "V/N"
        );
        assert_eq!(
            cell(Validate::Permissive, Normalize::Errors).website_glyph(),
            "V"
        );
        assert_eq!(cell(Validate::Yes, Normalize::No).website_glyph(), "V");
        assert_eq!(
            cell(Validate::Yes, Normalize::UnsupportedByDesign).website_glyph(),
            "V"
        );
        assert_eq!(cell(Validate::No, Normalize::No).website_glyph(), "-");
        assert_eq!(cell(Validate::No, Normalize::Full).website_glyph(), "N");
    }

    #[test]
    fn markdown_glyphs() {
        assert_eq!(cell(Validate::Yes, Normalize::Full).markdown_glyph(), "✓");
        assert_eq!(
            cell(Validate::Yes, Normalize::NetworkOnly).markdown_glyph(),
            "Net"
        );
        assert_eq!(
            cell(Validate::Permissive, Normalize::Errors).markdown_glyph(),
            "✗"
        );
        assert_eq!(
            cell(Validate::Yes, Normalize::UnsupportedByDesign).markdown_glyph(),
            "✗"
        );
        assert_eq!(cell(Validate::Yes, Normalize::No).markdown_glyph(), "✗");
    }

    const SAMPLE: &str = r#"{
      "schema_version": 1,
      "tools": [
        {"id":"ferro","name":"ferro","version":"0.6.0","url":"https://x"},
        {"id":"hgvs-rs","name":"hgvs-rs","version":"0.20.2","url":"https://y"}
      ],
      "footnotes": {"net":"needs network"},
      "categories": [
        {"id":"coordinate_types","title":"Coordinate Types","rows":[
          {"key":"p","label":"Protein","example":"p.Gly12Val","last_reviewed":"2026-06-09",
           "support":{"ferro":{"validate":"yes","normalize":"full"},"hgvs-rs":{"validate":"permissive","normalize":"errors"}}}
        ]}
      ],
      "views": []
    }"#;

    #[test]
    fn loads_and_validates() {
        let m = Matrix::from_json_str(SAMPLE).expect("parse");
        m.validate_invariants().expect("invariants");
        assert_eq!(m.tools.len(), 2);
    }

    #[test]
    fn rejects_missing_tool_in_row() {
        // Drop the hgvs-rs cell -> invariant must fail.
        let bad = SAMPLE.replace(
            r#","hgvs-rs":{"validate":"permissive","normalize":"errors"}"#,
            "",
        );
        let m = Matrix::from_json_str(&bad).expect("parse");
        assert!(m.validate_invariants().is_err());
    }

    #[test]
    fn rejects_unknown_footnote() {
        let bad = SAMPLE.replace(
            r#""normalize":"full"}"#,
            r#""normalize":"full","footnote":"nope"}"#,
        );
        let m = Matrix::from_json_str(&bad).expect("parse");
        assert!(m.validate_invariants().is_err());
    }

    fn view_matrix() -> Matrix {
        let json = r#"{
          "schema_version":1,
          "tools":[
            {"id":"ferro","name":"ferro","version":"0","url":"u"},
            {"id":"mutalyzer","name":"mutalyzer","version":"0","url":"u"},
            {"id":"biocommons","name":"biocommons","version":"0","url":"u"},
            {"id":"hgvs-rs","name":"hgvs-rs","version":"0","url":"u"}
          ],
          "footnotes":{
            "rewrite":"enabled by default via genomic-context rewriting; disable with --no-rewrite-intronic.",
            "net":"requires network access for NP_->NM_ lookups."
          },
          "categories":[{"id":"coordinate_types","title":"Coordinate Types","rows":[
            {"key":"intronic","label":"Intronic","last_reviewed":"2026-06-09","support":{
              "ferro":{"validate":"yes","normalize":"full"},
              "mutalyzer":{"validate":"yes","normalize":"full","footnote":"rewrite"},
              "biocommons":{"validate":"yes","normalize":"no"},
              "hgvs-rs":{"validate":"permissive","normalize":"unsupported_by_design"}
            }},
            {"key":"p","label":"Protein","last_reviewed":"2026-06-09","support":{
              "ferro":{"validate":"yes","normalize":"full"},
              "mutalyzer":{"validate":"yes","normalize":"network_only","footnote":"net"},
              "biocommons":{"validate":"yes","normalize":"unsupported_by_design"},
              "hgvs-rs":{"validate":"permissive","normalize":"errors"}
            }}
          ]}],
          "views":[{"id":"normalization_capabilities","title":"Normalization Capabilities",
            "row_label_header":"Pattern Type","select":[
              {"category":"coordinate_types","key":"intronic","label":"Coding (c.) intronic"},
              {"category":"coordinate_types","key":"p","label":"Protein (p.)"}
            ]}]
        }"#;
        Matrix::from_json_str(json).unwrap()
    }

    #[test]
    fn renders_markdown_view() {
        let m = view_matrix();
        let out = m
            .render_markdown_view("normalization_capabilities")
            .unwrap();
        // Header + alignment.
        assert!(out.contains("| Pattern Type | ferro | mutalyzer | biocommons | hgvs-rs |"));
        assert!(out.contains("|--------------|:-----:|:---------:|:----------:|:-------:|"));
        // Intronic row: ferro ✓, mutalyzer ✓ with first footnote (*), biocommons ✗, hgvs-rs ✗.
        assert!(out.contains("| Coding (c.) intronic | ✓ | ✓* | ✗ | ✗ |"));
        // Protein row: mutalyzer Net with second footnote (**).
        assert!(out.contains("| Protein (p.) | ✓ | Net** | ✗ | ✗ |"));
        // Footnote legend, in assignment order.
        assert!(out.contains("\\* enabled by default via genomic-context rewriting"));
        assert!(out.contains("\\*\\* requires network access"));
    }

    #[test]
    fn renders_website_json() {
        let m = view_matrix();
        let v = m.render_website_json();
        assert_eq!(v["tools"][0]["name"], "ferro");
        let cats = v["categories"].as_array().unwrap();
        let coord = &cats[0];
        assert_eq!(coord["title"], "Coordinate Types");
        let rows = coord["rows"].as_array().unwrap();
        // Intronic row glyphs in tool order.
        assert_eq!(rows[0]["label"], "Intronic");
        assert_eq!(rows[0]["cells"][0]["glyph"], "V/N"); // ferro
        assert_eq!(rows[0]["cells"][2]["glyph"], "V"); // biocommons (no normalize)
        assert_eq!(rows[0]["cells"][3]["glyph"], "V"); // hgvs-rs (unsupported)
                                                       // Protein: hgvs-rs errors -> V, biocommons unsupported -> V.
        assert_eq!(rows[1]["cells"][3]["glyph"], "V");
    }

    #[test]
    fn staleness_scan_flags_old() {
        let mut m = view_matrix();
        m.categories[0].rows[0].last_reviewed = "2000-01-01".into();
        let stale = m.stale_cells("2026-06-09", 6);
        assert_eq!(stale.len(), 1);
        assert!(stale[0].contains("intronic"));
    }

    #[test]
    fn months_before_wraps_year() {
        assert_eq!(months_before("2026-01-15", 1), "2025-12-15");
        assert_eq!(months_before("2026-06-09", 6), "2025-12-09");
    }

    /// `stale_cells` is documented "Non-failing": malformed `today` strings must
    /// not panic — they should yield a degenerate cutoff and return normally.
    #[test]
    fn stale_cells_non_failing_on_malformed_today() {
        let m = view_matrix();
        // Empty string, year-only, and year+month all must return without panicking.
        let _ = m.stale_cells("", 6);
        let _ = m.stale_cells("2026", 6);
        let _ = m.stale_cells("2026-06", 6);
    }
}
