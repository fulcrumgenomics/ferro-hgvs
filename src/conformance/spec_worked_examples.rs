//! Deserialization schema for the spec worked-examples corpus
//! (`tests/fixtures/spec-worked-examples/cases.json`).
//!
//! These are the HGVS spec's **own published answers**: rows where the
//! recommendation names an input and states, in the same sentence, what the
//! description must be instead ("use `X` and **not** `Y`"). They are the only
//! normalization rules in the spec whose correct output is not a matter of
//! reading — the document prints it.
//!
//! They are kept apart from the harvested spec fixture
//! (`tests/fixtures/grammar/hgvs_spec_normalization.json`) because that corpus
//! is normalized through `MockProvider::new()` — no reference bases — so no
//! sequence-dependent rule can fire there at all. The 3'rule and the repeat
//! typing restriction are both sequence-dependent, which is why they had never
//! actually been executed. This corpus runs them against real bases, committed
//! as a hermetic slice so the gate blocks a PR rather than skipping.
//!
//! Single source of truth for the schema, shared by the fixture extractor
//! (`examples/extract_spec_worked_example_windows.rs`) and the integration
//! harness (`tests/it/spec_worked_examples.rs`) — the same arrangement, and for
//! the same reason, as [`super::biocommons`].

use std::path::Path;

use serde::Deserialize;

use crate::FerroError;

/// Repo-relative path to the corpus.
pub const CASES_PATH: &str = "tests/fixtures/spec-worked-examples/cases.json";
/// Repo-relative path to the committed hermetic reference slice.
pub const WINDOWS_PATH: &str = "tests/fixtures/spec-worked-examples/reference-windows.json";
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
}

/// One worked example, pinned under each of ferro's three shipped
/// input-hygiene modes.
///
/// Three outputs rather than one because the modes do not agree here, and which
/// mode produces which answer *is* the finding: the single divergence from the
/// spec belongs to the lenient preprocessor's W3013 repair, not to the
/// normalizer, and strict mode rejects two of the spec's own inputs outright.
/// One pinned output would have hidden both facts.
///
/// Every case also carries a [`Citation`], because a pinned output with no
/// authority behind it is a change detector rather than a record: it tells a
/// future reader that something moved, not whether the old answer or the new
/// one is right.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Case {
    /// The HGVS description handed to `normalize`.
    pub input: String,
    /// The spec's own published answer for this input. `None` only where the
    /// spec publishes none — a state no row is currently in, and one that is
    /// harder to be in than it looks: the `LRG_199t1:c.3922dup` row sat here
    /// for exactly that reason until the
    /// `exon-junction-dup-converge-from-the-far-side` ruling found the spec
    /// answering it three times over (`DNA/duplication.md:26`, `:60`, `:148`).
    /// Setting this to `None` is an adjudication, not a fixture edit.
    pub spec_expected: Option<String>,
    /// What ferro produces through the plain library API (`parse_hgvs` +
    /// `Normalizer::new`) — no input-hygiene preprocessing, lenient
    /// normalization. `None` means it errors.
    pub default_output: Option<String>,
    /// What ferro produces through the shipped `--error-mode lenient` path
    /// (preprocessor on). `None` means it errors.
    pub lenient_output: Option<String>,
    /// What ferro produces through the shipped `--error-mode strict` path.
    /// `None` means strict rejects the input.
    pub strict_output: Option<String>,
    /// Where the spec says it. Verified verbatim against the checkout.
    pub citation: Citation,
    /// Why this row is here, and — where the modes disagree, or where a sibling
    /// row is deliberately not asserted against it — why that is the record.
    pub note: String,
}

impl Case {
    /// Whether the shipped lenient path reproduces the spec's published answer.
    /// `false` for a row the spec answers and ferro answers differently, and
    /// for a row the spec does not answer.
    pub fn lenient_matches_spec(&self) -> bool {
        match (&self.spec_expected, &self.lenient_output) {
            (Some(spec), Some(actual)) => spec == actual,
            _ => false,
        }
    }
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
