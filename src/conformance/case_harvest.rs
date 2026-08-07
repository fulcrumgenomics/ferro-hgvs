//! Deserialization schema for the harvested-case corpus
//! (`tests/fixtures/case-harvest/cases.json`).
//!
//! # What this corpus is
//!
//! An audit of this repository's issues, campaign notes and PR descriptions
//! turned up roughly 150 concrete HGVS reproducers that carry a *behavioural
//! claim* — "this input produces that output, and here is why that is right or
//! wrong" — and that **nothing running in CI can fail on**. Each had been
//! measured once, written into prose, and then left with no executable home:
//! the manifest-backed conformance axes skip without `FERRO_MANIFEST` (which in
//! CI is always), the harvested spec fixture normalizes through an empty
//! `MockProvider` so no sequence-dependent rule can fire, and the synthetic
//! corpora in `examples/dump_normalized_corpus.rs` are structurally unable to
//! build several of the shapes (#1456, #1460, #1478).
//!
//! This corpus is the verified subset, re-measured against the prepared
//! reference and pinned against a **committed** slice of it, so the gate blocks
//! a PR instead of skipping.
//!
//! # Why a curated corpus rather than a generated one
//!
//! Because for the flagship row the generated corpus is structurally blind.
//! #1535 fixed `NM_004006.2:c.76_83inv` and could still report
//! `Representation-Change: 0 rows move over 5,761,302 real expressions` — not
//! because the fix was inert but because no generator in the tree emits a block
//! whose reverse complement coincides with its reference at the interior
//! columns. A hand-written row is the only instrument that can catch that
//! regression, which is exactly the argument for keeping one.
//!
//! # The shape of a row
//!
//! Every [`Case`] carries a [`Citation`] into the vendored spec, because a
//! pinned output with no authority behind it is a change detector rather than a
//! record: it tells a future reader that something moved, not whether the old
//! answer or the new one is right.
//!
//! Three things a row can be:
//!
//! * **pinned** — `expected` is `Some`, and the driver asserts ferro still
//!   produces it.
//! * **confluent** — `confluence_class` is `Some`, and every row in the class
//!   must produce the *same* string. Where the governing ruling is open, the
//!   class deliberately carries no `expected`, so the gate asserts agreement
//!   without freezing which of two legal forms wins.
//! * **red** — `defect` is `Some`. The row is a live reproducer: its `expected`
//!   records today's **wrong** answer so movement is still detected, and the
//!   `#[ignore]`d test that names the defect fails until it is fixed.
//!
//! Single source of truth for the schema, shared by the fixture extractor
//! (`examples/extract_case_harvest_windows.rs`) and the integration harness
//! (`tests/it/case_harvest.rs`) — the same arrangement, and for the same reason,
//! as [`super::spec_enumeration`] and [`super::biocommons`].

use std::path::Path;

use serde::Deserialize;

use crate::FerroError;

/// Repo-relative path to the corpus.
pub const CASES_PATH: &str = "tests/fixtures/case-harvest/cases.json";
/// Repo-relative path to the committed hermetic reference slice.
pub const WINDOWS_PATH: &str = "tests/fixtures/case-harvest/reference-windows.json";
/// Repo-relative path to the vendored HGVS spec checkout the citations resolve
/// against.
pub const SPEC_DIR: &str = "assets/hgvs-nomenclature";

/// Top-level corpus document.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fixture {
    /// What this corpus is, in prose.
    pub description: String,
    /// Vendored spec checkout the citations were verified against.
    pub spec_commit: String,
    pub cases: Vec<Case>,
}

impl Fixture {
    /// Load the corpus from a JSON file.
    pub fn from_json_path(path: &Path) -> Result<Self, FerroError> {
        let content = std::fs::read_to_string(path)?;
        Ok(serde_json::from_str(&content)?)
    }

    /// Every distinct confluence class in the corpus, in first-seen order.
    pub fn confluence_classes(&self) -> Vec<&str> {
        let mut seen: Vec<&str> = Vec::new();
        for case in &self.cases {
            if let Some(class) = case.confluence_class.as_deref() {
                if !seen.contains(&class) {
                    seen.push(class);
                }
            }
        }
        seen
    }

    /// The rows belonging to `class`.
    pub fn members_of(&self, class: &str) -> Vec<&Case> {
        self.cases
            .iter()
            .filter(|case| case.confluence_class.as_deref() == Some(class))
            .collect()
    }
}

/// One harvested reproducer.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// The HGVS description handed to the normalizer.
    pub input: String,
    /// The output this row pins, or `None` for a row whose only assertion is
    /// its confluence class.
    ///
    /// `None` is not an omission — it is the honest record for a pair whose
    /// governing ruling is open. Asserting the two spellings agree is a
    /// correctness claim the spec supports; asserting *which* of two legal forms
    /// they agree on would freeze a representation nobody has adjudicated.
    pub expected: Option<String>,
    /// Rows sharing a class denote one variant and must normalize to one
    /// string.
    #[serde(default)]
    pub confluence_class: Option<String>,
    /// Where the spec says it. Verified verbatim against the checkout.
    pub citation: Citation,
    /// Why this row is here: what is adjudicated, what is merely observed, and
    /// what a prior audit got wrong about it.
    pub note: String,
    /// Set on a row that reproduces a live defect. `expected` then records the
    /// current *wrong* answer, and the `#[ignore]`d guard naming
    /// [`Defect::issue`] fails until the defect is fixed.
    #[serde(default)]
    pub defect: Option<Defect>,
}

impl Case {
    /// Whether this row reproduces a live defect.
    pub fn is_red(&self) -> bool {
        self.defect.is_some()
    }

    /// Whether the corpus-wide separation-0 check must exempt this row.
    ///
    /// Exactly one row is exempt, and only because it *is* the reproducer for
    /// the violation the check looks for. Exempting it here is what lets the
    /// check run over the whole corpus rather than being switched off wholesale.
    pub fn allows_adjacent_members(&self) -> bool {
        matches!(
            self.defect.as_ref().map(|d| d.kind),
            Some(DefectKind::AdjacentMembers)
        )
    }
}

/// The live defect a red row reproduces.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Defect {
    /// Issue reference, e.g. `#1539`. Reproduced verbatim into the `#[ignore]`
    /// reason so a red test names its own ticket.
    pub issue: String,
    /// Which guard is expected to fail.
    pub kind: DefectKind,
    /// One line on what is wrong.
    pub summary: String,
}

/// Which property a red row violates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum DefectKind {
    /// Two spellings of one variant normalize to two different strings.
    Confluence,
    /// The output puts two members on consecutive nucleotides
    /// (`DNA/delins.md:16`).
    AdjacentMembers,
}

/// A pointer into the vendored spec checkout, carrying the text it points at.
///
/// Same shape and same contract as the spec fixture's citations: `clause` is
/// `<path-relative-to-spec-dir>:<line>` or `…:<start>-<end>`, and `quote` must
/// appear verbatim within those lines. The quote is what makes the citation
/// self-verifying — a bare line number resolves against any file long enough to
/// have that line, so bumping the spec submodule could silently leave a
/// citation pointing at unrelated prose.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Citation {
    pub clause: String,
    pub quote: String,
}

impl Citation {
    /// Resolve the cited lines out of `spec_dir` and confirm the quote appears
    /// in them. `Err` carries a message naming what went wrong.
    pub fn verify(&self, spec_dir: &Path) -> Result<(), String> {
        let (path, first, last) = self.parse_clause()?;
        let full = spec_dir.join(path);
        let text = std::fs::read_to_string(&full)
            .map_err(|e| format!("citation {:?}: read {}: {e}", self.clause, full.display()))?;
        let lines: Vec<&str> = text.lines().collect();
        if first == 0 || last > lines.len() || first > last {
            return Err(format!(
                "citation {:?}: {} has {} lines",
                self.clause,
                full.display(),
                lines.len()
            ));
        }
        let cited = lines[first - 1..last].join("\n");
        if cited.contains(&self.quote) {
            return Ok(());
        }
        Err(format!(
            "citation {:?} does not quote the spec. Cited text:\n{cited}\nQuote:\n{}",
            self.clause, self.quote
        ))
    }

    /// Split a `path:line` or `path:start-end` clause into its parts.
    fn parse_clause(&self) -> Result<(&str, usize, usize), String> {
        let (path, lines) = self
            .clause
            .rsplit_once(':')
            .ok_or_else(|| format!("citation {:?} is not `<path>:<line>`", self.clause))?;
        let (first, last) = match lines.split_once('-') {
            Some((a, b)) => (a, b),
            None => (lines, lines),
        };
        let parse = |s: &str| {
            s.parse::<usize>()
                .map_err(|_| format!("citation {:?} has a non-numeric line", self.clause))
        };
        Ok((path, parse(first)?, parse(last)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn citation(clause: &str, quote: &str) -> Citation {
        Citation {
            clause: clause.to_string(),
            quote: quote.to_string(),
        }
    }

    #[test]
    fn a_clause_without_a_line_is_rejected() {
        let err = citation("docs/recommendations/general.md", "anything")
            .parse_clause()
            .unwrap_err();
        assert!(err.contains("is not `<path>:<line>`"), "{err}");
    }

    #[test]
    fn a_single_line_clause_parses_as_a_one_line_range() {
        let cited = citation("a/b.md:41", "q");
        let parsed = cited.parse_clause().unwrap();
        assert_eq!(parsed, ("a/b.md", 41, 41));
    }

    #[test]
    fn a_range_clause_parses_both_endpoints() {
        let cited = citation("a/b.md:21-22", "q");
        let parsed = cited.parse_clause().unwrap();
        assert_eq!(parsed, ("a/b.md", 21, 22));
    }

    #[test]
    fn a_non_numeric_line_is_rejected() {
        let err = citation("a/b.md:top", "q").parse_clause().unwrap_err();
        assert!(err.contains("non-numeric line"), "{err}");
    }

    #[test]
    fn only_an_adjacent_members_defect_exempts_a_row_from_the_separation_check() {
        let mut case = Case {
            input: "x".to_string(),
            expected: None,
            confluence_class: None,
            citation: citation("a/b.md:1", "q"),
            note: "n".to_string(),
            defect: None,
        };
        assert!(!case.is_red());
        assert!(!case.allows_adjacent_members());

        case.defect = Some(Defect {
            issue: "#1".to_string(),
            kind: DefectKind::Confluence,
            summary: "s".to_string(),
        });
        assert!(case.is_red());
        assert!(
            !case.allows_adjacent_members(),
            "a confluence defect must not buy an exemption from delins.md:16"
        );

        case.defect = Some(Defect {
            issue: "#2".to_string(),
            kind: DefectKind::AdjacentMembers,
            summary: "s".to_string(),
        });
        assert!(case.allows_adjacent_members());
    }
}
