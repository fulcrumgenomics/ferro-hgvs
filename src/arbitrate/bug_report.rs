//! Bug-report builder: renders a GitHub issue title, a Markdown body, and a
//! prefilled *new-issue* URL from a completed [`Arbitration`] where ferro was
//! judged wrong. Task 12 (CLI) drives this to let a user file a bug with one
//! command.
//!
//! The prefilled URL carries **only** `title` + `body` query parameters.
//! GitHub's `labels=` query param 404s for read-only external users (the
//! entire audience for this feature, verified against GitHub's own
//! behavior) — so suggested labels are embedded in the body itself as an
//! HTML comment (`<!-- suggested labels: ... -->`) instead of a URL param.

use url::Url;

use crate::arbitrate::Arbitration;

/// Base URL for filing a new issue against the ferro-hgvs repo.
const NEW_ISSUE_BASE: &str = "https://github.com/fulcrumgenomics/ferro-hgvs/issues/new";

/// Prefilled-URL length above which browsers/GitHub may truncate or reject
/// the query string; past this, [`BugReport::new_issue_url`] falls back to a
/// short pointer body instead of the full report.
const MAX_URL_LEN: usize = 7500;

/// Builds a GitHub bug-report issue title, Markdown body, and prefilled
/// new-issue URL from a completed arbitration.
///
/// Callers are expected to construct this only for arbitrations where ferro
/// was judged wrong (e.g. `category == MutalyzerCorrect`); the builder
/// itself renders whatever `arbitration` it is given without checking that.
pub struct BugReport<'a> {
    /// The arbitration verdict the report is filed about.
    pub arbitration: &'a Arbitration,
    /// Free-text notes from the reporter, included verbatim if present.
    pub notes: Option<&'a str>,
    /// Whether to include environment details (ferro version, OS/arch) in
    /// the body. Defaults to omitted: this data can leak local machine
    /// information, so callers must opt in explicitly rather than have it
    /// leak by default. Does NOT include a prepared-reference identity —
    /// `BugReport` carries no field for one (no hardcoded placeholder is
    /// emitted in its place).
    pub include_environment: bool,
    /// Ferro's version string, included only when `include_environment`.
    pub ferro_version: &'a str,
}

/// The prefilled new-issue URL, or an over-length fallback when the fully
/// percent-encoded URL would exceed [`MAX_URL_LEN`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BugReportUrl {
    /// The rendered report was short enough to prefill both `title` and
    /// `body` directly.
    Url(String),
    /// The rendered body was too long to fit in a URL. `url` is a prefilled
    /// new-issue link whose body is a short pointer asking the reporter to
    /// paste the full report themselves; `body` is that full Markdown
    /// report for the caller to display or copy separately.
    OverLength { url: String, body: String },
}

impl<'a> BugReport<'a> {
    /// Render the issue title: `arbitration: <input> — ferro vs <tool>`.
    pub fn issue_title(&self) -> String {
        format!(
            "arbitration: {} — ferro vs {}",
            self.arbitration.input, self.arbitration.other.tool
        )
    }

    /// Render the full Markdown issue body: input, ferro output, the other
    /// tool's output, the verdict/compliance/category, governing spec
    /// citations, SPDI forms of both outputs, reporter notes (if any), a
    /// repro command, and — only when `include_environment` is set —
    /// environment details. Always ends with the suggested-labels HTML
    /// comment consumed by [`Self::new_issue_url`]'s label-free URL.
    pub fn issue_body(&self) -> String {
        self.render_body()
    }

    /// Build the prefilled new-issue URL carrying only `title` and `body`
    /// query parameters (never `labels`/`assignees`/etc. — see the module
    /// docs). Falls back to [`BugReportUrl::OverLength`] when the full
    /// report would make the URL exceed [`MAX_URL_LEN`] characters.
    pub fn new_issue_url(&self) -> BugReportUrl {
        let title = self.issue_title();
        let body = self.issue_body();
        let full_url = build_url(&title, &body);
        if full_url.len() <= MAX_URL_LEN {
            return BugReportUrl::Url(full_url);
        }
        let pointer_url = build_url(&title, &Self::over_length_pointer_body());
        BugReportUrl::OverLength {
            url: pointer_url,
            body,
        }
    }

    /// Short body used in place of the full report when the full report
    /// would push the prefilled URL past [`MAX_URL_LEN`].
    fn over_length_pointer_body() -> String {
        "This bug report is too long to prefill in the URL. Please replace this text with the \
full report (from `ferro arbitrate ... --bug-report`) before submitting.\n\n\
<!-- suggested labels: bug, arbitration -->"
            .to_string()
    }

    /// Assemble the body's Markdown sections in order, then join them with
    /// blank lines. Broken out of [`Self::issue_body`] only for naming
    /// clarity — the two are equivalent.
    fn render_body(&self) -> String {
        let a = self.arbitration;
        let mut sections: Vec<String> = Vec::new();

        sections.push(format!("## Input\n\n`{}`", a.input));

        sections.push(format!(
            "## Ferro output\n\n```\n{}\n```",
            a.ferro_output
                .as_deref()
                .unwrap_or("(none — ferro failed to parse)")
        ));

        sections.push(format!(
            "## {} output\n\n```\n{}\n```",
            a.other.tool,
            a.other.output.as_deref().unwrap_or("(none)")
        ));

        sections.push(format!(
            "## Verdict\n\n- Verdict: `{:?}`\n- Compliance: `{:?}`\n- Category: `{}`",
            a.verdict, a.compliance, a.category
        ));

        if !a.spec_citations.is_empty() {
            let citations: Vec<String> = a
                .spec_citations
                .iter()
                .map(|c| {
                    format!(
                        "> {}\n>\n> — *{}*, {} ({})",
                        c.excerpt.trim(),
                        c.heading,
                        c.file,
                        c.spec_version
                    )
                })
                .collect();
            sections.push(format!("## Spec citation\n\n{}", citations.join("\n\n")));
        }

        sections.push(format!(
            "## SPDI\n\n- ferro: `{}`\n- {}: `{}`",
            a.ferro_spdi.as_deref().unwrap_or("(unavailable)"),
            a.other.tool,
            a.other_spdi.as_deref().unwrap_or("(unavailable)"),
        ));

        if let Some(notes) = self.notes {
            sections.push(format!("## Notes\n\n{notes}"));
        }

        if self.include_environment {
            // No `Reference identity` line here: `BugReport` carries no field
            // for a real prepared-reference identity hash, so a hardcoded
            // placeholder string previously stood in for one — misleading,
            // since it read like real data. Removed rather than faked; add
            // a genuine identity field (and this line) if/when one exists.
            sections.push(format!(
                "## Environment\n\n- ferro version: `{}`\n- OS: `{} ({})`",
                self.ferro_version,
                std::env::consts::OS,
                std::env::consts::ARCH,
            ));
        }

        sections.push(format!(
            "## Reproduce\n\n```\nferro arbitrate {:?} --reference <reference-dir> \
--other-tool {:?} --other-output {:?}\n```",
            a.input,
            a.other.tool,
            a.other.output.as_deref().unwrap_or(""),
        ));

        sections.push("<!-- suggested labels: bug, arbitration -->".to_string());

        sections.join("\n\n")
    }
}

/// Build the prefilled new-issue URL for `title`/`body`, percent-encoding
/// both via [`url::Url::query_pairs_mut`] so neither can inject additional
/// query parameters (e.g. a `body` containing `&labels=` stays a literal,
/// encoded value, not a live param).
fn build_url(title: &str, body: &str) -> String {
    let mut url = Url::parse(NEW_ISSUE_BASE).expect("NEW_ISSUE_BASE is a valid, static URL");
    url.query_pairs_mut()
        .append_pair("title", title)
        .append_pair("body", body);
    url.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arbitrate::spec_citations::SpecCitation;
    use crate::arbitrate::{ArbitrationCategory, Compliance, OtherResult, Verdict};

    /// A minimal ferro-is-wrong `Arbitration`: a genuine disagreement where
    /// the other tool's spelling is spec-compliant and ferro's is not.
    fn sample() -> crate::arbitrate::Arbitration {
        Arbitration {
            input: "NC_000001.11:g.51dup".to_string(),
            verdict: Verdict::Different,
            compliance: Compliance::Other,
            category: ArbitrationCategory::MutalyzerCorrect,
            ferro_output: Some("NC_000001.11:g.51dup".to_string()),
            other: OtherResult {
                tool: "mutalyzer".to_string(),
                status: "ok".to_string(),
                output: Some("NC_000001.11:g.52dup".to_string()),
            },
            spec_citations: vec![SpecCitation {
                spec_version: "21.0.4".to_string(),
                file: "duplication.md".to_string(),
                heading: "Duplication".to_string(),
                excerpt: "A duplication must be described using the shortest possible 3' \
                          representation."
                    .to_string(),
            }],
            ferro_spdi: Some("NC_000001.11:50:A:AA".to_string()),
            other_spdi: Some("NC_000001.11:51:A:AA".to_string()),
            reason: None,
        }
    }

    #[test]
    fn url_has_only_title_and_body_params() {
        let a = sample();
        let br = BugReport {
            arbitration: &a,
            notes: None,
            include_environment: false,
            ferro_version: "0.7.1",
        };
        let url = match br.new_issue_url() {
            BugReportUrl::Url(u) => u,
            BugReportUrl::OverLength { url, .. } => url,
        };
        assert!(url.starts_with("https://github.com/fulcrumgenomics/ferro-hgvs/issues/new?"));
        assert!(url.contains("title=") && url.contains("body="));
        assert!(!url.contains("labels=")); // labels 404 for read-only users
        assert!(!url.contains("assignees="));
    }

    #[test]
    fn environment_omitted_by_default() {
        let a = sample();
        let br = BugReport {
            arbitration: &a,
            notes: None,
            include_environment: false,
            ferro_version: "0.7.1",
        };
        let body = br.issue_body();
        assert!(!body.to_lowercase().contains("os:"));
        assert!(!body.contains("reference identity"));
    }

    #[test]
    fn environment_included_when_opted_in() {
        let a = sample();
        let br = BugReport {
            arbitration: &a,
            notes: None,
            include_environment: true,
            ferro_version: "0.7.1",
        };
        let body = br.issue_body();
        assert!(body.to_lowercase().contains("os:"));
        // No fake "reference identity" line: `BugReport` has no field
        // carrying a real identity hash, so the environment section must
        // not claim to report one.
        assert!(!body.to_lowercase().contains("reference identity"));
        assert!(body.contains("0.7.1"));
    }

    #[test]
    fn notes_are_included_when_present() {
        let a = sample();
        let br = BugReport {
            arbitration: &a,
            notes: Some("Observed while normalizing a batch of ClinVar variants."),
            include_environment: false,
            ferro_version: "0.7.1",
        };
        let body = br.issue_body();
        assert!(body.contains("Observed while normalizing a batch of ClinVar variants."));
    }

    #[test]
    fn body_always_carries_suggested_labels_comment() {
        let a = sample();
        let br = BugReport {
            arbitration: &a,
            notes: None,
            include_environment: false,
            ferro_version: "0.7.1",
        };
        let body = br.issue_body();
        assert!(body.contains("<!-- suggested labels: bug, arbitration -->"));
    }

    #[test]
    fn over_length_body_falls_back() {
        let mut a = sample();
        a.ferro_output = Some("X".repeat(9000));
        let br = BugReport {
            arbitration: &a,
            notes: None,
            include_environment: false,
            ferro_version: "0.7.1",
        };
        match br.new_issue_url() {
            BugReportUrl::OverLength { url, body } => {
                assert!(url.len() <= MAX_URL_LEN);
                assert!(!url.contains("labels="));
                assert!(body.contains(&"X".repeat(9000)));
            }
            BugReportUrl::Url(_) => panic!("expected OverLength fallback for a 9000-char output"),
        }
    }

    #[test]
    fn short_report_does_not_fall_back() {
        let a = sample();
        let br = BugReport {
            arbitration: &a,
            notes: None,
            include_environment: false,
            ferro_version: "0.7.1",
        };
        assert!(matches!(br.new_issue_url(), BugReportUrl::Url(_)));
    }
}
