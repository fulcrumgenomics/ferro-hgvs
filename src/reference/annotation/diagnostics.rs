//! Loader diagnostics. See spec §7.

use std::collections::HashMap;
use std::path::PathBuf;

use serde::Serialize;

const SAMPLE_LIMIT: usize = 100;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Severity {
    Error,
    Warning,
    Info,
}

#[derive(Debug, Clone, Serialize)]
pub struct SourceLocation {
    pub path: PathBuf,
    pub line: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub enum DiagnosticPayload {
    MalformedRecord {
        field: &'static str,
        value: String,
    },
    OrphanFeature {
        parent_id: String,
        child_id: Option<String>,
    },
    NoExonsDerivable {
        transcript_id: String,
    },
    PhaseAppliedToCdsStart {
        transcript_id: String,
        phase: u8,
    },
    PhaseUnavailable {
        transcript_id: String,
    },
    StopCodonAssumed {
        transcript_id: String,
    },
    CdsLengthNotMod3 {
        transcript_id: String,
        length: u64,
    },
    NonCanonicalStartCodon {
        transcript_id: String,
        codon: String,
    },
    UnknownFormat,
    InlineFastaIgnored,
    StrandRequired {
        feature_id: Option<String>,
    },
    GeneAsTranscript {
        gene_id: String,
    },
    ConflictingProteinId {
        transcript_id: String,
        values: Vec<String>,
    },
    Other,
}

#[derive(Debug, Clone, Serialize)]
pub struct LoaderDiagnostic {
    pub code: &'static str,
    pub severity: Severity,
    pub message: String,
    pub source: SourceLocation,
    pub feature_id: Option<String>,
    pub payload: DiagnosticPayload,
}

impl LoaderDiagnostic {
    pub fn warning(
        code: &'static str,
        message: impl Into<String>,
        source: SourceLocation,
        feature_id: Option<String>,
        payload: DiagnosticPayload,
    ) -> Self {
        Self {
            code,
            severity: Severity::Warning,
            message: message.into(),
            source,
            feature_id,
            payload,
        }
    }

    pub fn error(
        code: &'static str,
        message: impl Into<String>,
        source: SourceLocation,
        feature_id: Option<String>,
        payload: DiagnosticPayload,
    ) -> Self {
        Self {
            code,
            severity: Severity::Error,
            message: message.into(),
            source,
            feature_id,
            payload,
        }
    }
}

#[derive(Debug, Default, Clone)]
#[non_exhaustive]
pub struct LoaderReport {
    pub transcripts_loaded: usize,
    pub records_dropped: usize,
    pub diagnostics_by_code: HashMap<&'static str, usize>,
    pub sample_diagnostics: Vec<LoaderDiagnostic>,
    /// True if any `Severity::Error` diagnostic has been recorded, including
    /// errors beyond the `sample_diagnostics` cap. Used by strict-mode callers.
    pub has_error: bool,
    /// When `true`, [`Self::record`] is a no-op. Set by [`super::load_annotations`]
    /// when [`super::ErrorMode::Silent`] is configured so that all recording sites
    /// (record parsing, orphan resolution, transcript builder) observe the same
    /// gate without each having to thread `ErrorMode` through.
    pub silent: bool,
}

impl LoaderReport {
    pub fn record(&mut self, diag: LoaderDiagnostic) {
        if self.silent {
            return;
        }
        *self.diagnostics_by_code.entry(diag.code).or_insert(0) += 1;
        if diag.severity == Severity::Error {
            self.has_error = true;
        }
        if self.sample_diagnostics.len() < SAMPLE_LIMIT {
            self.sample_diagnostics.push(diag.clone());
        }
        log::warn!(
            "{} ({}:{}): {}",
            diag.code,
            diag.source.path.display(),
            diag.source.line,
            diag.message
        );
    }

    pub fn summary_line(&self) -> String {
        if self.diagnostics_by_code.is_empty() {
            format!("Loaded {} transcripts", self.transcripts_loaded)
        } else {
            let mut parts: Vec<(&&'static str, &usize)> = self.diagnostics_by_code.iter().collect();
            parts.sort_by(|a, b| b.1.cmp(a.1));
            let summary: Vec<String> = parts.iter().map(|(c, n)| format!("{} {}", n, c)).collect();
            format!(
                "Loaded {} transcripts (dropped {}: {})",
                self.transcripts_loaded,
                self.records_dropped,
                summary.join(", ")
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagnostic_constructs_and_records_severity() {
        let diag = LoaderDiagnostic::warning(
            "W-LOAD-100",
            "TranscriptWithoutExons: tx1 has no children",
            SourceLocation {
                path: "test.gff".into(),
                line: 5,
            },
            Some("tx1".into()),
            DiagnosticPayload::NoExonsDerivable {
                transcript_id: "tx1".into(),
            },
        );
        assert_eq!(diag.code, "W-LOAD-100");
        assert_eq!(diag.severity, Severity::Warning);
        assert_eq!(diag.source.line, 5);
    }

    #[test]
    fn loader_report_aggregates_by_code() {
        let mut report = LoaderReport::default();
        report.record(LoaderDiagnostic::warning(
            "W-LOAD-100",
            "...",
            SourceLocation {
                path: "t.gff".into(),
                line: 1,
            },
            None,
            DiagnosticPayload::Other,
        ));
        report.record(LoaderDiagnostic::warning(
            "W-LOAD-100",
            "...",
            SourceLocation {
                path: "t.gff".into(),
                line: 2,
            },
            None,
            DiagnosticPayload::Other,
        ));
        assert_eq!(report.diagnostics_by_code.get("W-LOAD-100"), Some(&2));
        assert_eq!(report.sample_diagnostics.len(), 2);
    }
}
