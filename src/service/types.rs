//! Request and response types for the HGVS web service

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Re-export parsed variant types from benchmark module
pub use crate::benchmark::types::{ParsedVariantDetails, PositionDetails};

use crate::hgvs::edit::NaEdit;
use crate::hgvs::variant::HgvsVariant;

/// Extract parsed variant details from an [`HgvsVariant`].
///
/// Returns `Some` with a flattened [`ParsedVariantDetails`] for the seven
/// variant classes that map onto a single coordinate/position (`c`/`g`/`n`/`p`/
/// `r`/`m`/`o`). Returns `None` for the complex classes (allele lists, fusions,
/// genome rings, supernumerary contigs, null/unknown alleles) which have no
/// single coordinate/position to flatten into the [`ParsedVariantDetails`]
/// shape; this `None` is a "no structured breakdown" signal, not an error.
///
/// The complex classes are listed explicitly (no `_` arm) so that adding a new
/// [`HgvsVariant`] arm fails to compile here rather than being silently dropped.
///
/// This is the single shared implementation used by both the validate HTTP
/// handler and the ferro MCP tool so the two surfaces stay in lockstep.
pub(crate) fn extract_variant_details(variant: &HgvsVariant) -> Option<ParsedVariantDetails> {
    match variant {
        HgvsVariant::Cds(v) => {
            let (variant_type, deleted, inserted) = if let Some(edit) = v.loc_edit.edit.inner() {
                extract_na_edit_info(edit)
            } else {
                ("unknown".to_string(), None, None)
            };
            Some(ParsedVariantDetails {
                reference: v.accession.to_string(),
                coordinate_system: "c".to_string(),
                variant_type,
                position: PositionDetails {
                    start: 0, // Position encoded in display string
                    end: None,
                    offset: None,
                    display: v.loc_edit.location.to_string(),
                },
                deleted,
                inserted,
                was_shifted: None,
                original_position: None,
            })
        }
        HgvsVariant::Genome(v) => {
            let (variant_type, deleted, inserted) = if let Some(edit) = v.loc_edit.edit.inner() {
                extract_na_edit_info(edit)
            } else {
                ("unknown".to_string(), None, None)
            };
            Some(ParsedVariantDetails {
                reference: v.accession.to_string(),
                coordinate_system: "g".to_string(),
                variant_type,
                position: PositionDetails {
                    start: 0, // Position encoded in display string
                    end: None,
                    offset: None,
                    display: v.loc_edit.location.to_string(),
                },
                deleted,
                inserted,
                was_shifted: None,
                original_position: None,
            })
        }
        HgvsVariant::Tx(v) => {
            let (variant_type, deleted, inserted) = if let Some(edit) = v.loc_edit.edit.inner() {
                extract_na_edit_info(edit)
            } else {
                ("unknown".to_string(), None, None)
            };
            Some(ParsedVariantDetails {
                reference: v.accession.to_string(),
                coordinate_system: "n".to_string(),
                variant_type,
                position: PositionDetails {
                    start: 0, // Position encoded in display string
                    end: None,
                    offset: None,
                    display: v.loc_edit.location.to_string(),
                },
                deleted,
                inserted,
                was_shifted: None,
                original_position: None,
            })
        }
        HgvsVariant::Protein(v) => Some(ParsedVariantDetails {
            reference: v.accession.to_string(),
            coordinate_system: "p".to_string(),
            variant_type: "protein_change".to_string(),
            position: PositionDetails {
                start: 0,
                end: None,
                offset: None,
                display: v.loc_edit.location.to_string(),
            },
            deleted: None,
            inserted: None,
            was_shifted: None,
            original_position: None,
        }),
        HgvsVariant::Rna(v) => {
            let (variant_type, deleted, inserted) = if let Some(edit) = v.loc_edit.edit.inner() {
                extract_na_edit_info(edit)
            } else {
                ("unknown".to_string(), None, None)
            };
            Some(ParsedVariantDetails {
                reference: v.accession.to_string(),
                coordinate_system: "r".to_string(),
                variant_type,
                position: PositionDetails {
                    start: 0, // Position encoded in display string
                    end: None,
                    offset: None,
                    display: v.loc_edit.location.to_string(),
                },
                deleted,
                inserted,
                was_shifted: None,
                original_position: None,
            })
        }
        HgvsVariant::Mt(v) => {
            let (variant_type, deleted, inserted) = if let Some(edit) = v.loc_edit.edit.inner() {
                extract_na_edit_info(edit)
            } else {
                ("unknown".to_string(), None, None)
            };
            Some(ParsedVariantDetails {
                reference: v.accession.to_string(),
                coordinate_system: "m".to_string(),
                variant_type,
                position: PositionDetails {
                    start: 0, // Position encoded in display string
                    end: None,
                    offset: None,
                    display: v.loc_edit.location.to_string(),
                },
                deleted,
                inserted,
                was_shifted: None,
                original_position: None,
            })
        }
        HgvsVariant::Circular(v) => {
            let (variant_type, deleted, inserted) = if let Some(edit) = v.loc_edit.edit.inner() {
                extract_na_edit_info(edit)
            } else {
                ("unknown".to_string(), None, None)
            };
            Some(ParsedVariantDetails {
                reference: v.accession.to_string(),
                coordinate_system: "o".to_string(),
                variant_type,
                position: PositionDetails {
                    start: 0, // Position encoded in display string
                    end: None,
                    offset: None,
                    display: v.loc_edit.location.to_string(),
                },
                deleted,
                inserted,
                was_shifted: None,
                original_position: None,
            })
        }
        // These classes have no single coordinate/position to flatten into the
        // ParsedVariantDetails shape, so they intentionally have no structured
        // breakdown. Listed explicitly (no `_`) so a future HgvsVariant arm fails to
        // compile here rather than being silently dropped.
        HgvsVariant::RnaFusion(_)
        | HgvsVariant::GenomeRing(_)
        | HgvsVariant::Supernumerary(_)
        | HgvsVariant::Allele(_)
        | HgvsVariant::NullAllele
        | HgvsVariant::UnknownAllele => None,
    }
}

/// Extract the edit type label and the deleted/inserted sequences from a
/// [`NaEdit`], returning `(variant_type, deleted, inserted)`.
pub(crate) fn extract_na_edit_info(edit: &NaEdit) -> (String, Option<String>, Option<String>) {
    match edit {
        NaEdit::Substitution {
            reference,
            alternative,
        } => (
            "substitution".to_string(),
            Some(reference.to_string()),
            Some(alternative.to_string()),
        ),
        NaEdit::SubstitutionNoRef { alternative } => (
            "substitution".to_string(),
            None,
            Some(alternative.to_string()),
        ),
        NaEdit::Deletion { sequence, length } => {
            let deleted = sequence
                .as_ref()
                .map(|s| s.to_string())
                .or_else(|| length.map(|l| format!("{} bp", l)));
            ("deletion".to_string(), deleted, None)
        }
        NaEdit::Insertion { sequence } => {
            ("insertion".to_string(), None, Some(sequence.to_string()))
        }
        NaEdit::Delins {
            sequence,
            deleted,
            deleted_length,
            ..
        } => {
            let deleted = deleted
                .as_ref()
                .map(|s| s.to_string())
                .or_else(|| deleted_length.map(|l| format!("{} bp", l)));
            ("delins".to_string(), deleted, Some(sequence.to_string()))
        }
        NaEdit::Duplication {
            sequence, length, ..
        } => {
            let deleted = sequence
                .as_ref()
                .map(|s| s.to_string())
                .or_else(|| length.map(|l| format!("{} bp", l)));
            ("duplication".to_string(), deleted, None)
        }
        NaEdit::Inversion { sequence, length } => {
            let deleted = sequence
                .as_ref()
                .map(|s| s.to_string())
                .or_else(|| length.map(|l| format!("{} bp", l)));
            ("inversion".to_string(), deleted, None)
        }
        NaEdit::Repeat {
            sequence, count, ..
        } => {
            let seq = sequence.as_ref().map(|s| s.to_string());
            ("repeat".to_string(), seq, Some(format!("{}", count)))
        }
        NaEdit::Identity { .. } => ("identity".to_string(), None, None),
        NaEdit::Unknown { .. } => ("unknown".to_string(), None, None),
        // The remaining edit classes either carry no flat sequence to surface or
        // are non-standard/structural forms without a single deleted/inserted
        // pair. They are listed explicitly (no `_`) so a future `NaEdit` arm
        // fails to compile here rather than being silently mapped to "other".
        NaEdit::NPaddedDeletion { count } => {
            ("deletion".to_string(), Some(format!("N[{}]", count)), None)
        }
        NaEdit::BreakpointInsertion { sequence } => {
            ("insertion".to_string(), None, Some(sequence.to_string()))
        }
        NaEdit::DupIns { sequence } => ("dupins".to_string(), None, Some(sequence.to_string())),
        NaEdit::MultiRepeat { .. } => ("repeat".to_string(), None, None),
        NaEdit::Conversion { .. } => ("conversion".to_string(), None, None),
        NaEdit::Methylation { .. } => ("methylation".to_string(), None, None),
        NaEdit::CopyNumber { .. } => ("copy_number".to_string(), None, None),
        NaEdit::Splice { .. } => ("splice".to_string(), None, None),
        NaEdit::NoProduct => ("no_product".to_string(), None, None),
        NaEdit::PositionOnly => ("position_only".to_string(), None, None),
    }
}

/// Strongly-typed tool names to prevent injection attacks
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ToolName {
    #[serde(rename = "ferro")]
    Ferro,
    #[serde(rename = "mutalyzer")]
    Mutalyzer,
    #[serde(rename = "biocommons")]
    Biocommons,
    #[serde(rename = "hgvs-rs")]
    HgvsRs,
}

impl ToolName {
    /// Get the string representation of the tool name
    pub fn as_str(&self) -> &'static str {
        match self {
            ToolName::Ferro => "ferro",
            ToolName::Mutalyzer => "mutalyzer",
            ToolName::Biocommons => "biocommons",
            ToolName::HgvsRs => "hgvs-rs",
        }
    }

    /// Parse a tool name from a string (case-insensitive)
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "ferro" => Some(ToolName::Ferro),
            "mutalyzer" => Some(ToolName::Mutalyzer),
            "biocommons" => Some(ToolName::Biocommons),
            "hgvs-rs" | "hgvsrs" => Some(ToolName::HgvsRs),
            _ => None,
        }
    }

    /// Get all available tool names
    pub fn all() -> &'static [ToolName] {
        &[
            ToolName::Ferro,
            ToolName::Mutalyzer,
            ToolName::Biocommons,
            ToolName::HgvsRs,
        ]
    }
}

impl std::fmt::Display for ToolName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Structured error categorization to replace string-based matching
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum StructuredErrorCategory {
    /// Parse errors (invalid syntax, format issues)
    Parse(ParseErrorKind),
    /// Reference sequence errors (missing sequences, transcripts)
    Reference(ReferenceErrorKind),
    /// Validation errors (position out of bounds, invalid ranges)
    Validation(ValidationErrorKind),
    /// Tool-specific errors (unavailable, configuration issues)
    Tool(ToolErrorKind),
    /// Timeout errors
    Timeout,
    /// Internal system errors
    Internal,
}

/// Specific types of parsing errors
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParseErrorKind {
    /// Invalid HGVS syntax
    InvalidSyntax,
    /// Invalid or unknown accession number
    InvalidAccession,
    /// Invalid position specification
    InvalidPosition,
    /// Invalid edit/change specification
    InvalidEdit,
    /// Unknown or unsupported variant type
    UnknownVariantType,
}

/// Specific types of reference sequence errors
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReferenceErrorKind {
    /// Sequence not found in database
    SequenceNotFound,
    /// Transcript not found
    TranscriptNotFound,
    /// Reference sequence mismatch
    SequenceMismatch,
    /// Chromosome not found
    ChromosomeNotFound,
    /// Database connection issues
    DatabaseError,
}

/// Specific types of validation errors
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ValidationErrorKind {
    /// Position is out of sequence bounds
    PositionOutOfBounds,
    /// Invalid coordinate range
    InvalidRange,
    /// Variant type not supported by tool
    UnsupportedVariant,
    /// Intronic variants not supported
    IntronicNotSupported,
    /// Protein variants not supported
    ProteinNotSupported,
}

/// Tool-specific error types
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ToolErrorKind {
    /// Tool is unavailable or not responding
    Unavailable,
    /// Tool configuration error
    ConfigurationError,
    /// Tool execution failed
    ExecutionFailed,
    /// Tool version incompatible
    IncompatibleVersion,
}

impl std::fmt::Display for StructuredErrorCategory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            StructuredErrorCategory::Parse(kind) => write!(f, "parse_error_{:?}", kind),
            StructuredErrorCategory::Reference(kind) => write!(f, "reference_error_{:?}", kind),
            StructuredErrorCategory::Validation(kind) => write!(f, "validation_error_{:?}", kind),
            StructuredErrorCategory::Tool(kind) => write!(f, "tool_error_{:?}", kind),
            StructuredErrorCategory::Timeout => write!(f, "timeout"),
            StructuredErrorCategory::Internal => write!(f, "internal_error"),
        }
    }
}

/// Error mode for handling errors during normalization
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ErrorMode {
    /// Suppress all errors, return empty/null on failure
    Silent,
    /// Log warnings but continue processing (default)
    #[default]
    Lenient,
    /// Fail immediately on any error
    Strict,
}

impl ErrorMode {
    /// Get the string representation
    pub fn as_str(&self) -> &'static str {
        match self {
            ErrorMode::Silent => "silent",
            ErrorMode::Lenient => "lenient",
            ErrorMode::Strict => "strict",
        }
    }
}

impl std::fmt::Display for ErrorMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Request for single variant processing
#[derive(Debug, Clone, Deserialize)]
pub struct SingleRequest {
    /// The HGVS variant string to process
    pub hgvs: String,
    /// Optional list of tools to use (default: all available)
    pub tools: Option<Vec<ToolName>>,
    /// Optional timeout in seconds (default: 30)
    pub timeout_seconds: Option<u32>,
    /// Error handling mode (silent, lenient, strict)
    #[serde(default)]
    pub error_mode: ErrorMode,
}

/// Request for batch variant processing
#[derive(Debug, Clone, Deserialize)]
pub struct BatchRequest {
    /// List of HGVS variant strings to process
    pub variants: Vec<String>,
    /// Optional list of tools to use (default: all available)
    pub tools: Option<Vec<ToolName>>,
    /// Optional timeout in seconds per variant (default: 30)
    pub timeout_seconds: Option<u32>,
    /// Error handling mode (silent, lenient, strict)
    #[serde(default)]
    pub error_mode: ErrorMode,
}

/// Response for single variant processing
#[derive(Debug, Serialize)]
pub struct SingleResponse {
    /// Original input HGVS string
    pub input: String,
    /// Results from each tool
    pub results: Vec<ToolResult>,
    /// Agreement analysis
    pub agreement: AgreementSummary,
    /// Total processing time in milliseconds
    pub processing_time_ms: u64,
}

/// Response for batch variant processing
#[derive(Debug, Serialize)]
pub struct BatchResponse {
    /// Number of variants processed
    pub total_variants: usize,
    /// Number of variants successfully processed by at least one tool
    pub successful_variants: usize,
    /// Results for each variant
    pub results: Vec<VariantBatchResult>,
    /// Total processing time in milliseconds
    pub total_processing_time_ms: u64,
}

/// Result for a single variant in batch processing
#[derive(Debug, Serialize)]
pub struct VariantBatchResult {
    /// Original input HGVS string
    pub input: String,
    /// Results from each tool
    pub results: Vec<ToolResult>,
    /// Agreement analysis for this variant
    pub agreement: AgreementSummary,
}

/// Result from a single tool
#[derive(Debug, Serialize)]
pub struct ToolResult {
    /// Tool name (ferro, mutalyzer, biocommons, hgvs-rs)
    pub tool: ToolName,
    /// Whether the tool succeeded
    pub success: bool,
    /// Normalized/parsed output if successful
    pub output: Option<String>,
    /// Error message if failed
    pub error: Option<String>,
    /// Standardized error category
    pub error_category: Option<String>,
    /// Processing time for this tool in milliseconds
    pub elapsed_ms: u64,
    /// Parsed variant details (for ferro tool)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub details: Option<ParsedVariantDetails>,
}

/// Response for validate endpoint
#[derive(Debug, Serialize)]
pub struct ValidateResponse {
    /// Original input HGVS string
    pub input: String,
    /// Whether the variant is syntactically valid
    pub valid: bool,
    /// Validation errors/warnings if any
    #[serde(skip_serializing_if = "Option::is_none")]
    pub errors: Option<Vec<String>>,
    /// Parsed component breakdown (if valid)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub components: Option<ParsedVariantDetails>,
    /// Processing time in milliseconds
    pub processing_time_ms: u64,
    /// True if the input is a wraparound range on a circular reference
    /// (`m.` or `o.`), where the start position is greater than the
    /// end position per SVD-WG006. Omitted from the JSON payload for
    /// linear axes.
    #[serde(skip_serializing_if = "std::ops::Not::not", default)]
    pub wraps_origin: bool,
}

/// Agreement analysis across tools
#[derive(Debug, Serialize)]
pub struct AgreementSummary {
    /// Whether all successful tools produced the same output
    pub all_agree: bool,
    /// Number of tools that succeeded
    pub successful_tools: usize,
    /// Number of tools that failed
    pub failed_tools: usize,
    /// Unique outputs and which tools produced them
    pub outputs: HashMap<String, Vec<ToolName>>,
}

/// Tool health status
#[derive(Debug, Clone, Serialize)]
pub struct ToolStatus {
    /// Tool name
    pub tool: ToolName,
    /// Whether the tool is available
    pub available: bool,
    /// Status message
    pub status: String,
    /// Last health check time (ISO 8601)
    pub last_check: String,
    /// Mode (e.g., "api" or "local" for mutalyzer)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<String>,
}

/// Overall service health response
#[derive(Debug, Clone, Serialize)]
pub struct HealthResponse {
    /// Service status (healthy, degraded, unhealthy)
    pub status: String,
    /// Available tools
    pub available_tools: Vec<ToolName>,
    /// Unavailable tools
    pub unavailable_tools: Vec<ToolName>,
    /// Detailed tool status
    pub tools: Vec<ToolStatus>,
}

/// Detailed health check response with comprehensive test results
#[derive(Debug, Clone, Serialize)]
pub struct DetailedHealthResponse {
    /// Basic health response
    #[serde(flatten)]
    pub basic: HealthResponse,
    /// Detailed test results per tool
    pub test_results: Vec<ToolTestResults>,
}

/// Test results for a single tool
#[derive(Debug, Clone, Serialize)]
pub struct ToolTestResults {
    /// Tool name
    pub tool: ToolName,
    /// Number of passed tests (of applicable tests)
    pub passed: usize,
    /// Number of applicable tests (not N/A) - for health rate
    pub total: usize,
    /// Total number of all tests including N/A - for coverage rate
    pub total_tests: usize,
    /// Mode (e.g., "api" or "local" for mutalyzer)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mode: Option<String>,
    /// Detailed results by category
    pub categories: Vec<TestCategory>,
}

/// Test results for a category (reference types, coordinate types, variant types)
#[derive(Debug, Clone, Serialize)]
pub struct TestCategory {
    /// Category name (e.g., "Reference Types", "Variant Types")
    pub name: String,
    /// Individual test results
    pub tests: Vec<TestResult>,
}

/// Test result status
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum TestStatus {
    /// Test passed
    Pass,
    /// Test failed (tool supports this but returned error)
    Fail,
    /// Not applicable (tool doesn't support this variant type or missing data)
    Na,
}

/// Result of a single test variant
#[derive(Debug, Clone, Serialize)]
pub struct TestResult {
    /// Test name/label (e.g., "NM_ transcript", "Substitution")
    pub name: String,
    /// Example variant tested
    pub variant: String,
    /// Test status: pass, fail, or na
    pub status: TestStatus,
    /// Whether the test passed (for backwards compatibility)
    pub passed: bool,
    /// Error message if failed or N/A
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// Standard error response
#[derive(Debug, Serialize)]
pub struct ErrorResponse {
    /// Error type
    pub error: String,
    /// Human-readable error message
    pub message: String,
    /// Optional additional details
    pub details: Option<serde_json::Value>,
}

/// Service error types
#[derive(Debug, thiserror::Error)]
pub enum ServiceError {
    #[error("Tool unavailable: {0}")]
    ToolUnavailable(String),

    #[error("Invalid HGVS: {0}")]
    InvalidHgvs(String),

    #[error("Request timeout")]
    Timeout,

    #[error("Internal error: {0}")]
    InternalError(String),

    #[error("Configuration error: {0}")]
    ConfigError(String),

    #[error("Bad request: {0}")]
    BadRequest(String),

    #[error("Circuit breaker open - service temporarily unavailable")]
    CircuitBreakerOpen,
}

impl ServiceError {
    /// Convert to HTTP status code
    pub fn status_code(&self) -> u16 {
        match self {
            ServiceError::BadRequest(_) => 400,
            ServiceError::InvalidHgvs(_) => 400,
            ServiceError::Timeout => 408,
            ServiceError::ConfigError(_) => 500,
            ServiceError::InternalError(_) => 500,
            ServiceError::ToolUnavailable(_) => 503,
            ServiceError::CircuitBreakerOpen => 503,
        }
    }

    /// Convert to error response
    pub fn to_response(&self) -> ErrorResponse {
        ErrorResponse {
            error: match self {
                ServiceError::BadRequest(_) => "bad_request".to_string(),
                ServiceError::InvalidHgvs(_) => "invalid_hgvs".to_string(),
                ServiceError::Timeout => "timeout".to_string(),
                ServiceError::ConfigError(_) => "config_error".to_string(),
                ServiceError::InternalError(_) => "internal_error".to_string(),
                ServiceError::ToolUnavailable(_) => "tool_unavailable".to_string(),
                ServiceError::CircuitBreakerOpen => "circuit_breaker_open".to_string(),
            },
            message: self.to_string(),
            details: None,
        }
    }
}

/// Error analysis module with tool-specific analyzers
pub mod error_analysis {
    use super::*;

    /// Trait for analyzing tool-specific errors
    pub trait ErrorAnalyzer {
        fn analyze_error(&self, error: &str) -> StructuredErrorCategory;
    }

    /// Ferro error analyzer using structured error codes
    pub struct FerroErrorAnalyzer;

    impl ErrorAnalyzer for FerroErrorAnalyzer {
        fn analyze_error(&self, error: &str) -> StructuredErrorCategory {
            // Ferro uses structured error codes E1xxx, E2xxx, E3xxx, etc.
            if let Some(error_code) = extract_ferro_error_code(error) {
                match error_code / 1000 {
                    1 => {
                        // E1xxx = Parse errors
                        match error_code {
                            1001 => {
                                StructuredErrorCategory::Parse(ParseErrorKind::InvalidAccession)
                            }
                            1002 => StructuredErrorCategory::Parse(ParseErrorKind::InvalidSyntax),
                            1003 => StructuredErrorCategory::Parse(ParseErrorKind::InvalidPosition),
                            1004 => StructuredErrorCategory::Parse(ParseErrorKind::InvalidEdit),
                            _ => StructuredErrorCategory::Parse(ParseErrorKind::InvalidSyntax),
                        }
                    }
                    2 => {
                        // E2xxx = Reference errors
                        match error_code {
                            2001 => StructuredErrorCategory::Reference(
                                ReferenceErrorKind::SequenceNotFound,
                            ),
                            2002 => StructuredErrorCategory::Reference(
                                ReferenceErrorKind::TranscriptNotFound,
                            ),
                            2003 => StructuredErrorCategory::Reference(
                                ReferenceErrorKind::SequenceMismatch,
                            ),
                            _ => StructuredErrorCategory::Reference(
                                ReferenceErrorKind::SequenceNotFound,
                            ),
                        }
                    }
                    3 => {
                        // E3xxx = Validation errors
                        match error_code {
                            3001 => StructuredErrorCategory::Validation(
                                ValidationErrorKind::PositionOutOfBounds,
                            ),
                            3002 => StructuredErrorCategory::Validation(
                                ValidationErrorKind::InvalidRange,
                            ),
                            3003 => StructuredErrorCategory::Validation(
                                ValidationErrorKind::UnsupportedVariant,
                            ),
                            _ => StructuredErrorCategory::Validation(
                                ValidationErrorKind::PositionOutOfBounds,
                            ),
                        }
                    }
                    _ => StructuredErrorCategory::Internal,
                }
            } else {
                // Fallback to string analysis for non-coded errors
                analyze_generic_error(error)
            }
        }
    }

    /// Mutalyzer error analyzer
    pub struct MutalyzerErrorAnalyzer;

    impl ErrorAnalyzer for MutalyzerErrorAnalyzer {
        fn analyze_error(&self, error: &str) -> StructuredErrorCategory {
            let lower_error = error.to_lowercase();

            if lower_error.contains("parse") || lower_error.contains("syntax") {
                StructuredErrorCategory::Parse(ParseErrorKind::InvalidSyntax)
            } else if lower_error.contains("sequence") && lower_error.contains("not") {
                StructuredErrorCategory::Reference(ReferenceErrorKind::SequenceNotFound)
            } else if lower_error.contains("transcript") && lower_error.contains("not") {
                StructuredErrorCategory::Reference(ReferenceErrorKind::TranscriptNotFound)
            } else if lower_error.contains("esequencemismatch") || lower_error.contains("mismatch")
            {
                StructuredErrorCategory::Reference(ReferenceErrorKind::SequenceMismatch)
            } else if lower_error.contains("position") && lower_error.contains("out") {
                StructuredErrorCategory::Validation(ValidationErrorKind::PositionOutOfBounds)
            } else if lower_error.contains("unavailable") || lower_error.contains("connection") {
                StructuredErrorCategory::Tool(ToolErrorKind::Unavailable)
            } else {
                analyze_generic_error(error)
            }
        }
    }

    /// Biocommons error analyzer
    pub struct BiocommonsErrorAnalyzer;

    impl ErrorAnalyzer for BiocommonsErrorAnalyzer {
        fn analyze_error(&self, error: &str) -> StructuredErrorCategory {
            let lower_error = error.to_lowercase();

            if lower_error.contains("parse") || lower_error.contains("invalid hgvs") {
                StructuredErrorCategory::Parse(ParseErrorKind::InvalidSyntax)
            } else if lower_error.contains("retrieval") || lower_error.contains("not found") {
                StructuredErrorCategory::Reference(ReferenceErrorKind::SequenceNotFound)
            } else if lower_error.contains("range") || lower_error.contains("bounds") {
                StructuredErrorCategory::Validation(ValidationErrorKind::PositionOutOfBounds)
            } else if lower_error.contains("not_supported") || lower_error.contains("unsupported") {
                StructuredErrorCategory::Validation(ValidationErrorKind::UnsupportedVariant)
            } else if lower_error.contains("esequencemismatch") || lower_error.contains("mismatch")
            {
                StructuredErrorCategory::Reference(ReferenceErrorKind::SequenceMismatch)
            } else if lower_error.contains("subprocess") || lower_error.contains("python") {
                StructuredErrorCategory::Tool(ToolErrorKind::ExecutionFailed)
            } else {
                analyze_generic_error(error)
            }
        }
    }

    /// HGVS-RS error analyzer
    pub struct HgvsRsErrorAnalyzer;

    impl ErrorAnalyzer for HgvsRsErrorAnalyzer {
        fn analyze_error(&self, error: &str) -> StructuredErrorCategory {
            let lower_error = error.to_lowercase();

            if lower_error.contains("parse") || lower_error.contains("syntax") {
                StructuredErrorCategory::Parse(ParseErrorKind::InvalidSyntax)
            } else if lower_error.contains("validation") {
                StructuredErrorCategory::Validation(ValidationErrorKind::UnsupportedVariant)
            } else if lower_error.contains("connection") || lower_error.contains("database") {
                StructuredErrorCategory::Reference(ReferenceErrorKind::DatabaseError)
            } else if lower_error.contains("transcript") {
                StructuredErrorCategory::Reference(ReferenceErrorKind::TranscriptNotFound)
            } else if lower_error.contains("intronic") {
                StructuredErrorCategory::Validation(ValidationErrorKind::IntronicNotSupported)
            } else if lower_error.contains("protein") {
                StructuredErrorCategory::Validation(ValidationErrorKind::ProteinNotSupported)
            } else if lower_error.contains("panic") || lower_error.contains("thread") {
                StructuredErrorCategory::Internal
            } else {
                analyze_generic_error(error)
            }
        }
    }

    /// Extract ferro error code from error message (E1234 -> 1234)
    fn extract_ferro_error_code(error: &str) -> Option<u32> {
        // Look for pattern like E1234 or E12345
        for word in error.split_whitespace() {
            if let Some(code_part) = word.strip_prefix('E') {
                if let Ok(code) = code_part.parse::<u32>() {
                    if (1000..=9999).contains(&code) {
                        return Some(code);
                    }
                }
            }
        }
        None
    }

    /// Generic error analysis for common patterns
    fn analyze_generic_error(error: &str) -> StructuredErrorCategory {
        let lower_error = error.to_lowercase();

        if lower_error.contains("timeout") {
            StructuredErrorCategory::Timeout
        } else if lower_error.contains("unavailable") || lower_error.contains("unreachable") {
            StructuredErrorCategory::Tool(ToolErrorKind::Unavailable)
        } else if lower_error.contains("configuration") || lower_error.contains("config") {
            StructuredErrorCategory::Tool(ToolErrorKind::ConfigurationError)
        } else {
            StructuredErrorCategory::Internal
        }
    }

    /// Get the appropriate error analyzer for a tool
    pub fn get_error_analyzer(tool: ToolName) -> Box<dyn ErrorAnalyzer> {
        match tool {
            ToolName::Ferro => Box::new(FerroErrorAnalyzer),
            ToolName::Mutalyzer => Box::new(MutalyzerErrorAnalyzer),
            ToolName::Biocommons => Box::new(BiocommonsErrorAnalyzer),
            ToolName::HgvsRs => Box::new(HgvsRsErrorAnalyzer),
        }
    }
}

/// Analyze error using structured error categories (replaces categorize_error)
pub fn analyze_error_structured(tool: ToolName, error: &str) -> StructuredErrorCategory {
    let analyzer = error_analysis::get_error_analyzer(tool);
    analyzer.analyze_error(error)
}

/// Convert structured error category to string for backward compatibility
pub fn categorize_error(tool: ToolName, error: &str) -> String {
    analyze_error_structured(tool, error).to_string()
}

// ==================== Coordinate Conversion Types ====================

/// Target coordinate system for conversion
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CoordinateSystem {
    /// Coding DNA (c.)
    C,
    /// Genomic (g.)
    G,
    /// Protein (p.)
    P,
    /// Non-coding (n.)
    N,
}

impl CoordinateSystem {
    pub fn as_str(&self) -> &'static str {
        match self {
            CoordinateSystem::C => "c",
            CoordinateSystem::G => "g",
            CoordinateSystem::P => "p",
            CoordinateSystem::N => "n",
        }
    }
}

impl std::fmt::Display for CoordinateSystem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Request for coordinate conversion
#[derive(Debug, Clone, Deserialize)]
pub struct ConvertRequest {
    /// The HGVS variant string to convert
    pub hgvs: String,
    /// Target coordinate system (c, g, p, n)
    pub target_system: CoordinateSystem,
    /// Include all possible conversions (default: false)
    #[serde(default)]
    pub include_all: bool,
}

/// Response for coordinate conversion
#[derive(Debug, Serialize)]
pub struct ConvertResponse {
    /// Original input HGVS string
    pub input: String,
    /// Source coordinate system detected
    pub source_system: String,
    /// Target coordinate system requested
    pub target_system: String,
    /// Primary converted result
    pub converted: Option<String>,
    /// All possible conversions (if include_all was true)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all_conversions: Option<Vec<ConversionResult>>,
    /// Error message if conversion failed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Processing time in milliseconds
    pub processing_time_ms: u64,
}

/// A single conversion result
#[derive(Debug, Serialize)]
pub struct ConversionResult {
    /// Target coordinate system
    pub system: String,
    /// Converted HGVS expression
    pub hgvs: String,
    /// Reference sequence used
    pub reference: String,
}

// ==================== Effect Prediction Types ====================

/// Request for effect prediction
#[derive(Debug, Clone, Deserialize)]
pub struct EffectRequest {
    /// The HGVS variant string to analyze
    pub hgvs: String,
    /// Include NMD (nonsense-mediated decay) prediction (default: false)
    #[serde(default)]
    pub include_nmd: bool,
}

/// Response for effect prediction
#[derive(Debug, Serialize)]
pub struct EffectResponse {
    /// Original input HGVS string
    pub input: String,
    /// Predicted effect using Sequence Ontology terms
    pub effect: Option<SequenceEffect>,
    /// Protein consequence if applicable
    #[serde(skip_serializing_if = "Option::is_none")]
    pub protein_consequence: Option<ProteinConsequence>,
    /// NMD prediction if requested
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nmd_prediction: Option<NmdPrediction>,
    /// Error message if prediction failed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Processing time in milliseconds
    pub processing_time_ms: u64,
}

/// Sequence effect using Sequence Ontology terms
#[derive(Debug, Serialize)]
pub struct SequenceEffect {
    /// Sequence Ontology term ID (e.g., "SO:0001583")
    pub so_term: String,
    /// Human-readable effect name (e.g., "missense_variant")
    pub name: String,
    /// Effect description
    pub description: String,
    /// Impact level (HIGH, MODERATE, LOW, MODIFIER)
    pub impact: String,
}

/// Protein consequence details
#[derive(Debug, Serialize)]
pub struct ProteinConsequence {
    /// Protein HGVS notation (e.g., "p.Arg117His")
    pub hgvs_p: String,
    /// Reference amino acid(s)
    pub ref_aa: String,
    /// Alternate amino acid(s)
    pub alt_aa: String,
    /// Position in protein
    pub position: u64,
    /// Whether this is a frameshift
    pub is_frameshift: bool,
}

/// NMD (Nonsense-Mediated Decay) prediction
#[derive(Debug, Serialize)]
pub struct NmdPrediction {
    /// Whether NMD is predicted
    pub predicted: bool,
    /// Confidence score (0.0-1.0)
    pub confidence: f64,
    /// Reasoning for prediction
    pub reason: String,
}

// ==================== Liftover Types ====================

/// Genome build for liftover
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GenomeBuild {
    #[serde(rename = "GRCh37", alias = "hg19")]
    GRCh37,
    #[serde(rename = "GRCh38", alias = "hg38")]
    GRCh38,
}

impl GenomeBuild {
    pub fn as_str(&self) -> &'static str {
        match self {
            GenomeBuild::GRCh37 => "GRCh37",
            GenomeBuild::GRCh38 => "GRCh38",
        }
    }
}

impl std::fmt::Display for GenomeBuild {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// Request for liftover
#[derive(Debug, Clone, Deserialize)]
pub struct LiftoverRequest {
    /// Genomic position (e.g., "chr7:117120148" or "NC_000007.13:g.117120148")
    pub position: String,
    /// Source genome build
    pub from_build: GenomeBuild,
    /// Target genome build
    pub to_build: GenomeBuild,
}

/// Response for liftover
#[derive(Debug, Serialize)]
pub struct LiftoverResponse {
    /// Original input position
    pub input: String,
    /// Source genome build
    pub from_build: String,
    /// Target genome build
    pub to_build: String,
    /// Converted position
    pub converted: Option<String>,
    /// Converted HGVS genomic notation
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hgvs_g: Option<String>,
    /// Chain file region used
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain_region: Option<String>,
    /// Error message if liftover failed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Processing time in milliseconds
    pub processing_time_ms: u64,
}

// ==================== VCF Conversion Types ====================

/// Request to convert VCF to HGVS
#[derive(Debug, Clone, Deserialize)]
pub struct VcfToHgvsRequest {
    /// Chromosome (e.g., "chr7" or "7")
    pub chrom: String,
    /// Position (1-based)
    pub pos: u64,
    /// Reference allele
    #[serde(rename = "ref")]
    pub ref_allele: String,
    /// Alternate allele
    pub alt: String,
    /// Genome build (default: GRCh38)
    #[serde(default = "default_grch38")]
    pub build: GenomeBuild,
    /// Transcript to use for c. notation (optional)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transcript: Option<String>,
}

fn default_grch38() -> GenomeBuild {
    GenomeBuild::GRCh38
}

/// Request to convert HGVS to VCF
#[derive(Debug, Clone, Deserialize)]
pub struct HgvsToVcfRequest {
    /// HGVS variant string
    pub hgvs: String,
    /// Genome build for output (default: GRCh38)
    #[serde(default = "default_grch38")]
    pub build: GenomeBuild,
}

/// Response for VCF to HGVS conversion
#[derive(Debug, Serialize)]
pub struct VcfToHgvsResponse {
    /// Original VCF representation
    pub vcf: VcfRecord,
    /// Genomic HGVS notation
    pub hgvs_g: Option<String>,
    /// Coding HGVS notation (if transcript provided)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hgvs_c: Option<String>,
    /// Protein HGVS notation (if applicable)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hgvs_p: Option<String>,
    /// Error message if conversion failed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Processing time in milliseconds
    pub processing_time_ms: u64,
}

/// Response for HGVS to VCF conversion
#[derive(Debug, Serialize)]
pub struct HgvsToVcfResponse {
    /// Original HGVS input
    pub input: String,
    /// VCF representation
    pub vcf: Option<VcfRecord>,
    /// Error message if conversion failed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Processing time in milliseconds
    pub processing_time_ms: u64,
}

/// VCF record representation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VcfRecord {
    /// Chromosome
    pub chrom: String,
    /// Position (1-based)
    pub pos: u64,
    /// Reference allele
    #[serde(rename = "ref")]
    pub ref_allele: String,
    /// Alternate allele
    pub alt: String,
    /// Genome build
    pub build: String,
}

/// Standardized health check system
pub mod health_check {
    use super::*;

    /// Result of a health check operation
    #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
    pub enum HealthCheckResult {
        /// Tool is fully operational
        Healthy,
        /// Tool is operational but with limitations or warnings
        Degraded { reason: String },
        /// Tool is not operational
        Unhealthy { reason: String },
    }

    impl HealthCheckResult {
        /// Check if the tool is available for use (healthy or degraded)
        pub fn is_available(&self) -> bool {
            matches!(
                self,
                HealthCheckResult::Healthy | HealthCheckResult::Degraded { .. }
            )
        }

        /// Get status string for API responses
        pub fn status_string(&self) -> &'static str {
            match self {
                HealthCheckResult::Healthy => "healthy",
                HealthCheckResult::Degraded { .. } => "degraded",
                HealthCheckResult::Unhealthy { .. } => "unhealthy",
            }
        }
    }

    /// Trait for standardized health checking
    #[async_trait::async_trait]
    pub trait HealthChecker {
        /// Perform health check for the tool
        async fn check_availability(&self) -> HealthCheckResult;

        /// Get the tool name
        fn tool_name(&self) -> ToolName;
    }

    /// Health check configuration for different tools
    #[derive(Debug, Clone)]
    pub struct HealthCheckConfig {
        /// Test variant to use for health checks
        pub test_variant: String,
        /// Timeout for health check operation
        pub timeout_seconds: u64,
        /// Expected behaviors that indicate healthy operation
        pub expected_behaviors: Vec<ExpectedBehavior>,
    }

    /// Expected behaviors for health check validation
    #[derive(Debug, Clone)]
    pub enum ExpectedBehavior {
        /// Tool should successfully process the test variant
        Success,
        /// Tool may fail but with specific error types that indicate it's working
        AcceptableFailure(Vec<StructuredErrorCategory>),
        /// Tool should respond within timeout (even if it fails)
        Responsive,
    }

    impl Default for HealthCheckConfig {
        fn default() -> Self {
            Self {
                test_variant: "NM_000001.2:c.1A>G".to_string(),
                timeout_seconds: 10,
                expected_behaviors: vec![
                    ExpectedBehavior::Success,
                    ExpectedBehavior::AcceptableFailure(vec![
                        StructuredErrorCategory::Reference(ReferenceErrorKind::SequenceNotFound),
                        StructuredErrorCategory::Reference(ReferenceErrorKind::TranscriptNotFound),
                        StructuredErrorCategory::Validation(
                            ValidationErrorKind::UnsupportedVariant,
                        ),
                    ]),
                    ExpectedBehavior::Responsive,
                ],
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn extract(hgvs: &str) -> Option<ParsedVariantDetails> {
        let result = crate::hgvs::parser::parse_hgvs_lenient(hgvs).unwrap();
        extract_variant_details(&result.result)
    }

    #[test]
    fn test_extract_variant_details_cds() {
        let d = extract("NM_000249.4:c.350C>T").expect("cds details");
        assert_eq!(d.coordinate_system, "c");
        assert_eq!(d.variant_type, "substitution");
        assert_eq!(d.reference, "NM_000249.4");
    }

    #[test]
    fn test_extract_variant_details_genomic() {
        let d = extract("NC_000007.14:g.117559593G>A").expect("genomic details");
        assert_eq!(d.coordinate_system, "g");
        assert_eq!(d.variant_type, "substitution");
    }

    #[test]
    fn test_extract_variant_details_noncoding() {
        let d = extract("NR_000001.1:n.100A>G").expect("noncoding details");
        assert_eq!(d.coordinate_system, "n");
    }

    #[test]
    fn test_extract_variant_details_protein() {
        let d = extract("NP_000240.1:p.Val600Glu").expect("protein details");
        assert_eq!(d.coordinate_system, "p");
        assert_eq!(d.variant_type, "protein_change");
    }

    #[test]
    fn test_extract_variant_details_rna() {
        let d = extract("NR_000001.1:r.100a>g").expect("rna details should be populated");
        assert_eq!(d.coordinate_system, "r");
    }

    #[test]
    fn test_extract_variant_details_mitochondrial() {
        let d =
            extract("NC_012920.1:m.8993T>G").expect("mitochondrial details should be populated");
        assert_eq!(d.coordinate_system, "m");
    }

    #[test]
    fn test_extract_variant_details_circular() {
        let d = extract("J01749.1:o.100A>G").expect("circular details should be populated");
        assert_eq!(d.coordinate_system, "o");
    }

    #[test]
    fn test_extract_variant_details_allele_returns_none() {
        // Compound alleles have no single coordinate/position to flatten into the
        // ParsedVariantDetails shape, so they intentionally return None.
        let result =
            crate::hgvs::parser::parse_hgvs_lenient("NM_000088.3:c.[10A>G;20C>T]").unwrap();
        assert!(matches!(result.result, HgvsVariant::Allele(_)));
        assert!(extract_variant_details(&result.result).is_none());
    }

    fn edit_info(hgvs: &str) -> (String, Option<String>, Option<String>) {
        let result = crate::hgvs::parser::parse_hgvs_lenient(hgvs).unwrap();
        match &result.result {
            HgvsVariant::Cds(v) => {
                let edit = v.loc_edit.edit.inner().expect("inner edit");
                extract_na_edit_info(edit)
            }
            other => panic!("expected a CDS variant, got {other:?}"),
        }
    }

    #[test]
    fn test_extract_na_edit_info_substitution() {
        let (vtype, deleted, inserted) = edit_info("NM_000249.4:c.350C>T");
        assert_eq!(vtype, "substitution");
        assert_eq!(deleted, Some("C".to_string()));
        assert_eq!(inserted, Some("T".to_string()));
    }

    #[test]
    fn test_extract_na_edit_info_deletion() {
        // W3025 (DelExplicitSeq) strips the explicit sequence in lenient mode,
        // so the parsed NaEdit::Deletion has sequence=None after preprocessing.
        let (vtype, deleted, inserted) = edit_info("NM_000249.4:c.350delC");
        assert_eq!(vtype, "deletion");
        assert_eq!(deleted, None);
        assert!(inserted.is_none());
    }

    #[test]
    fn test_extract_na_edit_info_insertion() {
        let (vtype, deleted, inserted) = edit_info("NM_000249.4:c.350_351insATG");
        assert_eq!(vtype, "insertion");
        assert!(deleted.is_none());
        assert_eq!(inserted, Some("ATG".to_string()));
    }

    #[test]
    fn test_extract_na_edit_info_delins() {
        let (vtype, deleted, inserted) = edit_info("NM_000249.4:c.350delinsATG");
        assert_eq!(vtype, "delins");
        assert_eq!(deleted, None, "short form has no explicit deleted");
        assert_eq!(inserted, Some("ATG".to_string()));
    }

    #[test]
    fn test_extract_na_edit_info_delins_with_explicit_deleted_seq() {
        let (vtype, deleted, inserted) = edit_info("NM_000249.4:c.350_352delATGinsTTCC");
        assert_eq!(vtype, "delins");
        assert_eq!(deleted, Some("ATG".to_string()));
        assert_eq!(inserted, Some("TTCC".to_string()));
    }

    #[test]
    fn test_extract_na_edit_info_delins_with_explicit_deleted_length() {
        let (vtype, deleted, inserted) = edit_info("NM_000249.4:c.350_352del3insTA");
        assert_eq!(vtype, "delins");
        assert_eq!(deleted, Some("3 bp".to_string()));
        assert_eq!(inserted, Some("TA".to_string()));
    }

    #[test]
    fn test_extract_na_edit_info_duplication() {
        // W3024 (DupExplicitSeq) strips the explicit sequence in lenient mode,
        // so the parsed NaEdit::Duplication has sequence=None after preprocessing.
        let (vtype, deleted, _inserted) = edit_info("NM_000249.4:c.350dupC");
        assert_eq!(vtype, "duplication");
        assert_eq!(deleted, None);
    }
}
