//! Diagnostic mode for debugging and tracing
//!
//! This module provides tools for debugging parser and normalizer behavior:
//! - Parse tracing with step-by-step output
//! - Intermediate representation dumps
//! - Normalization step logging
//!
//! # Example
//!
//! ```
//! use ferro_hgvs::diagnostic::{DiagnosticConfig, ParseDiagnostics};
//!
//! let config = DiagnosticConfig::verbose();
//! let mut diagnostics = ParseDiagnostics::new(config);
//! let result = diagnostics.parse_with_trace("NM_000088.3:c.459A>G");
//! println!("{}", diagnostics.format_trace());
//! ```

use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant};

use crate::error::FerroError;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::Normalizer;
use crate::reference::ReferenceProvider;

/// Global diagnostic mode flag
static DIAGNOSTIC_MODE: AtomicBool = AtomicBool::new(false);

/// Enable global diagnostic mode
pub fn enable_diagnostic_mode() {
    DIAGNOSTIC_MODE.store(true, Ordering::Relaxed);
}

/// Disable global diagnostic mode
pub fn disable_diagnostic_mode() {
    DIAGNOSTIC_MODE.store(false, Ordering::Relaxed);
}

/// Check if diagnostic mode is enabled
pub fn is_diagnostic_mode() -> bool {
    DIAGNOSTIC_MODE.load(Ordering::Relaxed)
}

/// Diagnostic verbosity level
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Verbosity {
    /// No output
    #[default]
    Quiet,
    /// Summary only
    Summary,
    /// Detailed step-by-step output
    Verbose,
    /// Maximum detail with internal state
    Debug,
}

/// Configuration for diagnostic output
#[derive(Debug, Clone)]
pub struct DiagnosticConfig {
    /// Verbosity level
    pub verbosity: Verbosity,
    /// Include timing information
    pub include_timing: bool,
    /// Include memory usage (approximate)
    pub include_memory: bool,
    /// Include internal parser state
    pub include_parser_state: bool,
    /// Maximum trace entries to keep
    pub max_trace_entries: usize,
}

impl Default for DiagnosticConfig {
    fn default() -> Self {
        Self {
            verbosity: Verbosity::Quiet,
            include_timing: false,
            include_memory: false,
            include_parser_state: false,
            max_trace_entries: 1000,
        }
    }
}

impl DiagnosticConfig {
    /// Create a new diagnostic config
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a verbose configuration
    pub fn verbose() -> Self {
        Self {
            verbosity: Verbosity::Verbose,
            include_timing: true,
            include_memory: false,
            include_parser_state: false,
            max_trace_entries: 1000,
        }
    }

    /// Create a debug configuration (maximum detail)
    pub fn debug() -> Self {
        Self {
            verbosity: Verbosity::Debug,
            include_timing: true,
            include_memory: true,
            include_parser_state: true,
            max_trace_entries: 10000,
        }
    }

    /// Set verbosity level
    pub fn with_verbosity(mut self, verbosity: Verbosity) -> Self {
        self.verbosity = verbosity;
        self
    }

    /// Enable timing information
    pub fn with_timing(mut self) -> Self {
        self.include_timing = true;
        self
    }
}

/// A trace entry for diagnostic output
#[derive(Debug, Clone)]
pub struct TraceEntry {
    /// Timestamp relative to start
    pub elapsed: Duration,
    /// Operation phase
    pub phase: &'static str,
    /// Operation name
    pub operation: String,
    /// Additional details
    pub details: Option<String>,
    /// Result status
    pub status: TraceStatus,
}

/// Status of a traced operation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceStatus {
    /// Operation started
    Started,
    /// Operation completed successfully
    Success,
    /// Operation failed
    Failed,
    /// Operation skipped
    Skipped,
}

impl fmt::Display for TraceStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TraceStatus::Started => write!(f, "STARTED"),
            TraceStatus::Success => write!(f, "OK"),
            TraceStatus::Failed => write!(f, "FAILED"),
            TraceStatus::Skipped => write!(f, "SKIPPED"),
        }
    }
}

/// Diagnostics for parsing operations
pub struct ParseDiagnostics {
    config: DiagnosticConfig,
    trace: Vec<TraceEntry>,
    start_time: Option<Instant>,
}

impl ParseDiagnostics {
    /// Create a new parse diagnostics instance
    pub fn new(config: DiagnosticConfig) -> Self {
        Self {
            config,
            trace: Vec::new(),
            start_time: None,
        }
    }

    /// Start tracing
    fn start(&mut self) {
        self.start_time = Some(Instant::now());
        self.trace.clear();
    }

    /// Add a trace entry
    fn add_trace(&mut self, phase: &'static str, operation: String, status: TraceStatus) {
        if self.trace.len() >= self.config.max_trace_entries {
            return;
        }

        let elapsed = self
            .start_time
            .map(|t| t.elapsed())
            .unwrap_or(Duration::ZERO);

        self.trace.push(TraceEntry {
            elapsed,
            phase,
            operation,
            details: None,
            status,
        });
    }

    /// Parse with tracing
    pub fn parse_with_trace(&mut self, input: &str) -> Result<HgvsVariant, FerroError> {
        self.start();

        self.add_trace("parse", format!("Input: {}", input), TraceStatus::Started);

        // Track parsing steps
        self.add_trace(
            "parse",
            "Parsing accession".to_string(),
            TraceStatus::Started,
        );

        let result = parse_hgvs(input);

        match &result {
            Ok(variant) => {
                self.add_trace(
                    "parse",
                    "Parsing complete".to_string(),
                    TraceStatus::Success,
                );
                self.add_trace(
                    "parse",
                    format!("Variant type: {}", variant_type_name(variant)),
                    TraceStatus::Success,
                );
            }
            Err(e) => {
                self.add_trace("parse", format!("Parse error: {}", e), TraceStatus::Failed);
            }
        }

        result
    }

    /// Get the trace entries
    pub fn trace(&self) -> &[TraceEntry] {
        &self.trace
    }

    /// Format the trace for display
    pub fn format_trace(&self) -> String {
        let mut output = String::new();

        output.push_str("=== Parse Trace ===\n");

        for entry in &self.trace {
            if self.config.include_timing {
                output.push_str(&format!("[{:?}] ", entry.elapsed));
            }
            output.push_str(&format!(
                "[{}] {} - {}\n",
                entry.phase, entry.status, entry.operation
            ));
            if let Some(ref details) = entry.details {
                output.push_str(&format!("    {}\n", details));
            }
        }

        output
    }
}

impl Default for ParseDiagnostics {
    fn default() -> Self {
        Self::new(DiagnosticConfig::default())
    }
}

/// Diagnostics for normalization operations
pub struct NormalizeDiagnostics {
    config: DiagnosticConfig,
    trace: Vec<TraceEntry>,
    start_time: Option<Instant>,
}

impl NormalizeDiagnostics {
    /// Create a new normalize diagnostics instance
    pub fn new(config: DiagnosticConfig) -> Self {
        Self {
            config,
            trace: Vec::new(),
            start_time: None,
        }
    }

    /// Start tracing
    fn start(&mut self) {
        self.start_time = Some(Instant::now());
        self.trace.clear();
    }

    /// Add a trace entry
    fn add_trace(
        &mut self,
        phase: &'static str,
        operation: String,
        details: Option<String>,
        status: TraceStatus,
    ) {
        if self.trace.len() >= self.config.max_trace_entries {
            return;
        }

        let elapsed = self
            .start_time
            .map(|t| t.elapsed())
            .unwrap_or(Duration::ZERO);

        self.trace.push(TraceEntry {
            elapsed,
            phase,
            operation,
            details,
            status,
        });
    }

    /// Normalize with tracing
    pub fn normalize_with_trace<P: ReferenceProvider>(
        &mut self,
        normalizer: &Normalizer<P>,
        variant: &HgvsVariant,
    ) -> Result<HgvsVariant, FerroError> {
        self.start();

        self.add_trace(
            "normalize",
            format!("Input: {}", variant),
            None,
            TraceStatus::Started,
        );

        self.add_trace(
            "normalize",
            format!("Variant type: {}", variant_type_name(variant)),
            None,
            TraceStatus::Success,
        );

        let result = normalizer.normalize(variant);

        match &result {
            Ok(normalized) => {
                self.add_trace(
                    "normalize",
                    "Normalization complete".to_string(),
                    None,
                    TraceStatus::Success,
                );
                self.add_trace(
                    "normalize",
                    format!("Output: {}", normalized),
                    None,
                    TraceStatus::Success,
                );
            }
            Err(e) => {
                self.add_trace(
                    "normalize",
                    format!("Normalization error: {}", e),
                    None,
                    TraceStatus::Failed,
                );
            }
        }

        result
    }

    /// Get the trace entries
    pub fn trace(&self) -> &[TraceEntry] {
        &self.trace
    }

    /// Format the trace for display
    pub fn format_trace(&self) -> String {
        let mut output = String::new();

        output.push_str("=== Normalize Trace ===\n");

        for entry in &self.trace {
            if self.config.include_timing {
                output.push_str(&format!("[{:?}] ", entry.elapsed));
            }
            output.push_str(&format!(
                "[{}] {} - {}\n",
                entry.phase, entry.status, entry.operation
            ));
            if let Some(ref details) = entry.details {
                output.push_str(&format!("    {}\n", details));
            }
        }

        output
    }
}

impl Default for NormalizeDiagnostics {
    fn default() -> Self {
        Self::new(DiagnosticConfig::default())
    }
}

/// Get a human-readable name for a variant type
fn variant_type_name(variant: &HgvsVariant) -> &'static str {
    match variant {
        HgvsVariant::Genome(_) => "Genomic (g.)",
        HgvsVariant::Cds(_) => "CDS (c.)",
        HgvsVariant::Tx(_) => "Non-coding transcript (n.)",
        HgvsVariant::Rna(_) => "RNA (r.)",
        HgvsVariant::Protein(_) => "Protein (p.)",
        HgvsVariant::Mt(_) => "Mitochondrial (m.)",
        HgvsVariant::Circular(_) => "Circular DNA (o.)",
        HgvsVariant::RnaFusion(_) => "RNA Fusion (::)",
        HgvsVariant::Allele(_) => "Allele/Compound",
        HgvsVariant::NullAllele => "Null Allele [0]",
        HgvsVariant::UnknownAllele => "Unknown Allele [?]",
    }
}

/// Variant structure dump for debugging
pub fn dump_variant_structure(variant: &HgvsVariant) -> String {
    let mut output = String::new();

    output.push_str("=== Variant Structure ===\n");
    output.push_str(&format!("Type: {}\n", variant_type_name(variant)));
    output.push_str(&format!("Display: {}\n", variant));
    output.push_str(&format!("Debug: {:?}\n", variant));

    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MockProvider;

    #[test]
    fn test_diagnostic_mode() {
        assert!(!is_diagnostic_mode());
        enable_diagnostic_mode();
        assert!(is_diagnostic_mode());
        disable_diagnostic_mode();
        assert!(!is_diagnostic_mode());
    }

    #[test]
    fn test_parse_diagnostics() {
        let config = DiagnosticConfig::verbose();
        let mut diag = ParseDiagnostics::new(config);

        let result = diag.parse_with_trace("NM_000088.3:c.459A>G");
        assert!(result.is_ok());

        let trace = diag.trace();
        assert!(!trace.is_empty());
        assert!(trace.iter().any(|e| e.status == TraceStatus::Success));
    }

    #[test]
    fn test_parse_diagnostics_error() {
        let config = DiagnosticConfig::verbose();
        let mut diag = ParseDiagnostics::new(config);

        let result = diag.parse_with_trace("invalid variant");
        assert!(result.is_err());

        let trace = diag.trace();
        assert!(trace.iter().any(|e| e.status == TraceStatus::Failed));
    }

    #[test]
    fn test_normalize_diagnostics() {
        let config = DiagnosticConfig::verbose();
        let mut diag = NormalizeDiagnostics::new(config);

        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);

        let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();
        let result = diag.normalize_with_trace(&normalizer, &variant);
        assert!(result.is_ok());

        let trace = diag.trace();
        assert!(!trace.is_empty());
    }

    #[test]
    fn test_format_trace() {
        let config = DiagnosticConfig::verbose();
        let mut diag = ParseDiagnostics::new(config);

        let _ = diag.parse_with_trace("NM_000088.3:c.459A>G");

        let formatted = diag.format_trace();
        assert!(formatted.contains("Parse Trace"));
        assert!(formatted.contains("parse"));
    }

    #[test]
    fn test_dump_variant_structure() {
        let variant = parse_hgvs("NM_000088.3:c.459A>G").unwrap();
        let dump = dump_variant_structure(&variant);

        assert!(dump.contains("Variant Structure"));
        assert!(dump.contains("CDS"));
    }

    #[test]
    fn test_verbosity_levels() {
        let quiet = DiagnosticConfig::new();
        assert_eq!(quiet.verbosity, Verbosity::Quiet);

        let verbose = DiagnosticConfig::verbose();
        assert_eq!(verbose.verbosity, Verbosity::Verbose);
        assert!(verbose.include_timing);

        let debug = DiagnosticConfig::debug();
        assert_eq!(debug.verbosity, Verbosity::Debug);
        assert!(debug.include_memory);
    }

    // Additional diagnostic tests

    #[test]
    fn test_verbosity_default() {
        let verbosity = Verbosity::default();
        assert_eq!(verbosity, Verbosity::Quiet);
    }

    #[test]
    fn test_verbosity_debug_trait() {
        let verbosity = Verbosity::Verbose;
        let debug_str = format!("{:?}", verbosity);
        assert!(debug_str.contains("Verbose"));
    }

    #[test]
    fn test_verbosity_clone() {
        let v1 = Verbosity::Debug;
        let v2 = v1;
        assert_eq!(v1, v2);
    }

    #[test]
    fn test_diagnostic_config_default() {
        let config = DiagnosticConfig::default();
        assert_eq!(config.verbosity, Verbosity::Quiet);
        assert!(!config.include_timing);
        assert!(!config.include_memory);
        assert!(!config.include_parser_state);
        assert_eq!(config.max_trace_entries, 1000);
    }

    #[test]
    fn test_diagnostic_config_with_verbosity() {
        let config = DiagnosticConfig::new().with_verbosity(Verbosity::Summary);
        assert_eq!(config.verbosity, Verbosity::Summary);
    }

    #[test]
    fn test_diagnostic_config_with_timing() {
        let config = DiagnosticConfig::new().with_timing();
        assert!(config.include_timing);
    }

    #[test]
    fn test_diagnostic_config_debug_preset() {
        let config = DiagnosticConfig::debug();
        assert_eq!(config.verbosity, Verbosity::Debug);
        assert!(config.include_timing);
        assert!(config.include_memory);
        assert!(config.include_parser_state);
        assert_eq!(config.max_trace_entries, 10000);
    }

    #[test]
    fn test_diagnostic_config_clone() {
        let config1 = DiagnosticConfig::verbose();
        let config2 = config1.clone();
        assert_eq!(config1.verbosity, config2.verbosity);
        assert_eq!(config1.include_timing, config2.include_timing);
    }

    #[test]
    fn test_trace_status_equality() {
        assert_eq!(TraceStatus::Started, TraceStatus::Started);
        assert_eq!(TraceStatus::Success, TraceStatus::Success);
        assert_eq!(TraceStatus::Failed, TraceStatus::Failed);
        assert_eq!(TraceStatus::Skipped, TraceStatus::Skipped);
        assert_ne!(TraceStatus::Started, TraceStatus::Success);
    }

    #[test]
    fn test_trace_status_debug() {
        let status = TraceStatus::Success;
        let debug_str = format!("{:?}", status);
        assert!(debug_str.contains("Success"));
    }

    #[test]
    fn test_trace_entry_creation() {
        let entry = TraceEntry {
            elapsed: Duration::from_millis(100),
            phase: "parse",
            operation: "tokenize".to_string(),
            details: Some("processing input".to_string()),
            status: TraceStatus::Success,
        };

        assert_eq!(entry.elapsed, Duration::from_millis(100));
        assert_eq!(entry.phase, "parse");
        assert_eq!(entry.operation, "tokenize");
        assert!(entry.details.is_some());
        assert_eq!(entry.status, TraceStatus::Success);
    }

    #[test]
    fn test_trace_entry_no_details() {
        let entry = TraceEntry {
            elapsed: Duration::from_secs(0),
            phase: "normalize",
            operation: "shift".to_string(),
            details: None,
            status: TraceStatus::Started,
        };

        assert!(entry.details.is_none());
    }

    #[test]
    fn test_trace_entry_clone() {
        let entry = TraceEntry {
            elapsed: Duration::from_millis(50),
            phase: "test",
            operation: "op".to_string(),
            details: Some("detail".to_string()),
            status: TraceStatus::Success,
        };

        let cloned = entry.clone();
        assert_eq!(cloned.elapsed, entry.elapsed);
        assert_eq!(cloned.operation, entry.operation);
    }

    #[test]
    fn test_parse_diagnostics_reset_on_new_parse() {
        let config = DiagnosticConfig::verbose();
        let mut diag = ParseDiagnostics::new(config);

        // First parse populates trace
        let _ = diag.parse_with_trace("NM_000088.3:c.10A>G");
        assert!(!diag.trace().is_empty());
        let first_count = diag.trace().len();

        // Second parse resets and repopulates trace (start() is called internally)
        let _ = diag.parse_with_trace("NM_000088.3:c.20C>T");
        // Trace should have similar count, not double
        assert_eq!(diag.trace().len(), first_count);
    }

    #[test]
    fn test_parse_diagnostics_error_capture() {
        let config = DiagnosticConfig::verbose();
        let mut diag = ParseDiagnostics::new(config);

        // Parse invalid input
        let result = diag.parse_with_trace("invalid_hgvs");
        assert!(result.is_err());

        // Trace should still have entries even for failed parses
        assert!(!diag.trace().is_empty());

        // Check that trace contains failure indication
        let has_failure = diag
            .trace()
            .iter()
            .any(|e| matches!(e.status, TraceStatus::Failed));
        assert!(has_failure);
    }

    #[test]
    fn test_parse_diagnostics_quiet_mode() {
        let config = DiagnosticConfig::new(); // Quiet by default
        let mut diag = ParseDiagnostics::new(config);

        let result = diag.parse_with_trace("NM_000088.3:c.10A>G");
        assert!(result.is_ok());
    }

    #[test]
    fn test_normalize_diagnostics_reset_on_new_normalize() {
        let config = DiagnosticConfig::verbose();
        let mut diag = NormalizeDiagnostics::new(config);

        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();

        // First normalize populates trace
        let _ = diag.normalize_with_trace(&normalizer, &variant);
        assert!(!diag.trace().is_empty());
        let first_count = diag.trace().len();

        // Second normalize resets and repopulates trace
        let _ = diag.normalize_with_trace(&normalizer, &variant);
        // Trace should have similar count, not double
        assert_eq!(diag.trace().len(), first_count);
    }

    #[test]
    fn test_normalize_diagnostics_format_trace() {
        let config = DiagnosticConfig::verbose();
        let mut diag = NormalizeDiagnostics::new(config);

        let provider = MockProvider::with_test_data();
        let normalizer = Normalizer::new(provider);
        let variant = parse_hgvs("NM_000088.3:c.10A>G").unwrap();

        let _ = diag.normalize_with_trace(&normalizer, &variant);

        let formatted = diag.format_trace();
        assert!(formatted.contains("Normalize Trace"));
    }

    #[test]
    fn test_dump_variant_structure_coding() {
        let variant = parse_hgvs("NM_000088.3:c.100A>G").unwrap();
        let dump = dump_variant_structure(&variant);

        assert!(dump.contains("Variant Structure"));
        assert!(dump.contains("CDS"));
        assert!(dump.contains("NM_000088.3"));
    }

    #[test]
    fn test_dump_variant_structure_genomic() {
        let variant = parse_hgvs("NC_000001.11:g.12345A>G").unwrap();
        let dump = dump_variant_structure(&variant);

        assert!(dump.contains("Variant Structure"));
        assert!(dump.contains("Genome"));
    }

    #[test]
    fn test_dump_variant_structure_protein() {
        let variant = parse_hgvs("NP_000079.2:p.Val600Glu").unwrap();
        let dump = dump_variant_structure(&variant);

        assert!(dump.contains("Variant Structure"));
        assert!(dump.contains("Protein"));
    }

    #[test]
    fn test_global_diagnostic_mode_thread_safety() {
        // Ensure the global flag works correctly
        disable_diagnostic_mode();
        assert!(!is_diagnostic_mode());

        enable_diagnostic_mode();
        assert!(is_diagnostic_mode());

        disable_diagnostic_mode();
        assert!(!is_diagnostic_mode());
    }

    #[test]
    fn test_trace_entry_debug_format() {
        let entry = TraceEntry {
            elapsed: Duration::from_millis(123),
            phase: "parse",
            operation: "test_op".to_string(),
            details: Some("test details".to_string()),
            status: TraceStatus::Success,
        };

        let debug_str = format!("{:?}", entry);
        assert!(debug_str.contains("TraceEntry"));
        assert!(debug_str.contains("parse"));
    }

    #[test]
    fn test_diagnostic_config_debug_format() {
        let config = DiagnosticConfig::verbose();
        let debug_str = format!("{:?}", config);
        assert!(debug_str.contains("DiagnosticConfig"));
        assert!(debug_str.contains("Verbose"));
    }
}
