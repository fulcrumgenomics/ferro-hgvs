//! Error types for ferro-hgvs
//!
//! This module provides comprehensive error handling with:
//! - Error codes for categorization
//! - Source span tracking for error location
//! - Helpful diagnostic messages
//! - "Did you mean?" suggestions where applicable

use std::fmt;
use thiserror::Error;

/// Error codes for categorizing errors
///
/// These codes can be used for programmatic error handling
/// and for documentation lookup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u16)]
pub enum ErrorCode {
    // Parse errors (E1xxx)
    /// Invalid accession format
    InvalidAccession = 1001,
    /// Unknown variant type prefix
    UnknownVariantType = 1002,
    /// Invalid position format
    InvalidPosition = 1003,
    /// Invalid edit format
    InvalidEdit = 1004,
    /// Unexpected end of input
    UnexpectedEnd = 1005,
    /// Unexpected characters
    UnexpectedChar = 1006,
    /// Invalid base/nucleotide
    InvalidBase = 1007,
    /// Invalid amino acid
    InvalidAminoAcid = 1008,

    // Reference errors (E2xxx)
    /// Reference/transcript not found
    ReferenceNotFound = 2001,
    /// Sequence not available
    SequenceNotFound = 2002,
    /// Chromosome/contig not found
    ChromosomeNotFound = 2003,

    // Validation errors (E3xxx)
    /// Position out of bounds
    PositionOutOfBounds = 3001,
    /// Reference sequence mismatch
    ReferenceMismatch = 3002,
    /// Invalid coordinate range
    InvalidRange = 3003,
    /// Exon-intron boundary crossing
    ExonIntronBoundary = 3004,
    /// UTR-CDS boundary crossing
    UtrCdsBoundary = 3005,

    // Normalization errors (E4xxx)
    /// Intronic variant (cannot normalize)
    IntronicVariant = 4001,
    /// Unsupported variant type
    UnsupportedVariant = 4002,

    // Conversion errors (E5xxx)
    /// Coordinate conversion failed
    ConversionFailed = 5001,
    /// No overlapping transcript
    NoOverlappingTranscript = 5002,

    // IO errors (E9xxx)
    /// File IO error
    IoError = 9001,
    /// JSON parsing error
    JsonError = 9002,
}

impl ErrorCode {
    /// Get the error code as a string (e.g., "E1001")
    pub fn as_str(&self) -> String {
        format!("E{:04}", *self as u16)
    }

    /// Get a brief description of this error code
    pub fn description(&self) -> &'static str {
        match self {
            ErrorCode::InvalidAccession => "invalid accession format",
            ErrorCode::UnknownVariantType => "unknown variant type prefix",
            ErrorCode::InvalidPosition => "invalid position format",
            ErrorCode::InvalidEdit => "invalid edit format",
            ErrorCode::UnexpectedEnd => "unexpected end of input",
            ErrorCode::UnexpectedChar => "unexpected character",
            ErrorCode::InvalidBase => "invalid nucleotide base",
            ErrorCode::InvalidAminoAcid => "invalid amino acid",
            ErrorCode::ReferenceNotFound => "reference not found",
            ErrorCode::SequenceNotFound => "sequence not available",
            ErrorCode::ChromosomeNotFound => "chromosome not found",
            ErrorCode::PositionOutOfBounds => "position out of bounds",
            ErrorCode::ReferenceMismatch => "reference sequence mismatch",
            ErrorCode::InvalidRange => "invalid coordinate range",
            ErrorCode::ExonIntronBoundary => "variant crosses exon-intron boundary",
            ErrorCode::UtrCdsBoundary => "variant crosses UTR-CDS boundary",
            ErrorCode::IntronicVariant => "intronic variant not supported",
            ErrorCode::UnsupportedVariant => "unsupported variant type",
            ErrorCode::ConversionFailed => "coordinate conversion failed",
            ErrorCode::NoOverlappingTranscript => "no overlapping transcript",
            ErrorCode::IoError => "file I/O error",
            ErrorCode::JsonError => "JSON parsing error",
        }
    }
}

impl fmt::Display for ErrorCode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

/// A span in the source input indicating error location
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SourceSpan {
    /// Starting byte offset (0-indexed)
    pub start: usize,
    /// Ending byte offset (exclusive)
    pub end: usize,
}

impl SourceSpan {
    /// Create a new source span
    pub fn new(start: usize, end: usize) -> Self {
        Self { start, end }
    }

    /// Create a span for a single position
    pub fn point(pos: usize) -> Self {
        Self {
            start: pos,
            end: pos + 1,
        }
    }

    /// Format the source with the error highlighted
    ///
    /// Returns a string like:
    /// ```text
    /// NM_000088.3:c.459A>G
    ///                ^~~~
    /// ```
    pub fn highlight(&self, source: &str) -> String {
        if source.is_empty() {
            return String::new();
        }

        let safe_start = self.start.min(source.len());
        let safe_end = self.end.min(source.len()).max(safe_start);

        // Build the pointer line
        let mut pointer = String::with_capacity(source.len() + 4);
        for _ in 0..safe_start {
            pointer.push(' ');
        }
        if safe_start < safe_end {
            pointer.push('^');
            for _ in (safe_start + 1)..safe_end {
                pointer.push('~');
            }
        } else {
            pointer.push('^');
        }

        format!("{}\n{}", source, pointer)
    }
}

/// Diagnostic information for an error
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Diagnostic {
    /// Error code
    pub code: Option<ErrorCode>,
    /// Source span for highlighting
    pub span: Option<SourceSpan>,
    /// The original input (for error display)
    pub source: Option<String>,
    /// Helpful hint or suggestion
    pub hint: Option<String>,
    /// "Did you mean?" suggestion
    pub suggestion: Option<String>,
}

impl Diagnostic {
    /// Create a new empty diagnostic
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an error code
    pub fn with_code(mut self, code: ErrorCode) -> Self {
        self.code = Some(code);
        self
    }

    /// Add a source span
    pub fn with_span(mut self, span: SourceSpan) -> Self {
        self.span = Some(span);
        self
    }

    /// Add the original source
    pub fn with_source(mut self, source: impl Into<String>) -> Self {
        self.source = Some(source.into());
        self
    }

    /// Add a hint
    pub fn with_hint(mut self, hint: impl Into<String>) -> Self {
        self.hint = Some(hint.into());
        self
    }

    /// Add a suggestion
    pub fn with_suggestion(mut self, suggestion: impl Into<String>) -> Self {
        self.suggestion = Some(suggestion.into());
        self
    }

    /// Format the diagnostic as a detailed error message
    pub fn format(&self, primary_message: &str) -> String {
        let mut result = String::new();

        // Error code prefix
        if let Some(code) = &self.code {
            result.push_str(&format!("[{}] ", code));
        }

        // Primary message
        result.push_str(primary_message);

        // Source span highlight
        if let (Some(span), Some(source)) = (&self.span, &self.source) {
            result.push_str("\n\n");
            result.push_str(&span.highlight(source));
        }

        // Hint
        if let Some(hint) = &self.hint {
            result.push_str("\n\nHint: ");
            result.push_str(hint);
        }

        // Suggestion
        if let Some(suggestion) = &self.suggestion {
            result.push_str("\n\nDid you mean: ");
            result.push_str(suggestion);
        }

        result
    }
}

/// Main error type for ferro-hgvs operations
#[derive(Error, Debug, Clone, PartialEq)]
pub enum FerroError {
    /// Parse error with position and message
    #[error("Parse error at position {pos}: {msg}")]
    Parse {
        pos: usize,
        msg: String,
        /// Optional diagnostic with additional context
        diagnostic: Option<Box<Diagnostic>>,
    },

    /// Reference sequence or transcript not found
    #[error("Reference not found: {id}")]
    ReferenceNotFound { id: String },

    /// Variant spans an exon-intron boundary
    #[error("Variant spans exon-intron boundary at exon {exon}: {variant}")]
    ExonIntronBoundary { exon: u32, variant: String },

    /// Variant spans a UTR-CDS boundary
    #[error("Variant spans UTR-CDS boundary: {variant}")]
    UtrCdsBoundary { variant: String },

    /// Invalid coordinates provided
    #[error("Invalid coordinates: {msg}")]
    InvalidCoordinates { msg: String },

    /// Unsupported variant type
    #[error("Unsupported variant type: {variant_type}")]
    UnsupportedVariant { variant_type: String },

    /// Intronic variant cannot be normalized (no genomic data)
    #[error("Intronic variant normalization not supported: {variant}")]
    IntronicVariant { variant: String },

    /// Genomic reference data is not available
    #[error("Genomic reference not available for {contig}:{start}-{end}")]
    GenomicReferenceNotAvailable {
        contig: String,
        start: u64,
        end: u64,
    },

    /// Protein reference data is not available
    #[error("Protein reference not available for {accession}:{start}-{end}")]
    ProteinReferenceNotAvailable {
        accession: String,
        start: u64,
        end: u64,
    },

    /// Amino acid mismatch with reference
    #[error("Amino acid mismatch at position {position}: expected {expected}, found {found} in {accession}")]
    AminoAcidMismatch {
        accession: String,
        position: u64,
        expected: String,
        found: String,
    },

    /// Reference sequence mismatch
    #[error("Reference mismatch at {location}: expected {expected}, found {found}")]
    ReferenceMismatch {
        location: String,
        expected: String,
        found: String,
    },

    /// Coordinate conversion error
    #[error("Coordinate conversion error: {msg}")]
    ConversionError { msg: String },

    /// IO error (for file operations)
    #[error("IO error: {msg}")]
    Io { msg: String },

    /// JSON parsing error
    #[error("JSON error: {msg}")]
    Json { msg: String },
}

impl FerroError {
    /// Create a parse error with diagnostic information
    pub fn parse_with_diagnostic(
        pos: usize,
        msg: impl Into<String>,
        diagnostic: Diagnostic,
    ) -> Self {
        FerroError::Parse {
            pos,
            msg: msg.into(),
            diagnostic: Some(Box::new(diagnostic)),
        }
    }

    /// Create a simple parse error without diagnostic
    pub fn parse(pos: usize, msg: impl Into<String>) -> Self {
        FerroError::Parse {
            pos,
            msg: msg.into(),
            diagnostic: None,
        }
    }

    /// Get the error code if available
    pub fn code(&self) -> Option<ErrorCode> {
        match self {
            FerroError::Parse {
                diagnostic: Some(d),
                ..
            } => d.code,
            FerroError::ReferenceNotFound { .. } => Some(ErrorCode::ReferenceNotFound),
            FerroError::ExonIntronBoundary { .. } => Some(ErrorCode::ExonIntronBoundary),
            FerroError::UtrCdsBoundary { .. } => Some(ErrorCode::UtrCdsBoundary),
            FerroError::InvalidCoordinates { .. } => Some(ErrorCode::InvalidRange),
            FerroError::UnsupportedVariant { .. } => Some(ErrorCode::UnsupportedVariant),
            FerroError::IntronicVariant { .. } => Some(ErrorCode::IntronicVariant),
            FerroError::ReferenceMismatch { .. } => Some(ErrorCode::ReferenceMismatch),
            FerroError::ConversionError { .. } => Some(ErrorCode::ConversionFailed),
            FerroError::Io { .. } => Some(ErrorCode::IoError),
            FerroError::Json { .. } => Some(ErrorCode::JsonError),
            _ => None,
        }
    }

    /// Get a formatted error with full diagnostic output
    pub fn detailed_message(&self) -> String {
        match self {
            FerroError::Parse {
                pos,
                msg,
                diagnostic: Some(d),
            } => d.format(&format!("Parse error at position {}: {}", pos, msg)),
            _ => self.to_string(),
        }
    }
}

/// Helper to suggest similar variant type prefixes
pub fn suggest_variant_type(found: &str) -> Option<&'static str> {
    let found_lower = found.to_lowercase();
    let prefixes = [
        ("c", "c."),
        ("g", "g."),
        ("p", "p."),
        ("n", "n."),
        ("r", "r."),
        ("m", "m."),
        ("c:", "c."),
        ("g:", "g."),
        ("p:", "p."),
    ];

    for (pattern, suggestion) in prefixes {
        if found_lower.starts_with(pattern) {
            return Some(suggestion);
        }
    }
    None
}

/// Helper to suggest similar amino acids
pub fn suggest_amino_acid(found: &str) -> Option<&'static str> {
    let suggestions = [
        ("ale", "Ala"),
        ("arg", "Arg"),
        ("asg", "Asn"),
        ("asp", "Asp"),
        ("cis", "Cys"),
        ("cys", "Cys"),
        ("gln", "Gln"),
        ("glu", "Glu"),
        ("gly", "Gly"),
        ("his", "His"),
        ("ile", "Ile"),
        ("leu", "Leu"),
        ("lys", "Lys"),
        ("met", "Met"),
        ("phe", "Phe"),
        ("pro", "Pro"),
        ("sec", "Sec"),
        ("sel", "Sec"),
        ("ser", "Ser"),
        ("thr", "Thr"),
        ("trp", "Trp"),
        ("tyr", "Tyr"),
        ("val", "Val"),
        ("ter", "Ter"),
        ("stop", "Ter"),
        ("end", "Ter"),
        ("xaa", "Xaa"),
        ("unk", "Xaa"),
    ];

    let found_lower = found.to_lowercase();
    for (pattern, suggestion) in suggestions {
        if found_lower.starts_with(pattern) || pattern.starts_with(&found_lower) {
            return Some(suggestion);
        }
    }
    None
}

impl From<std::io::Error> for FerroError {
    fn from(err: std::io::Error) -> Self {
        FerroError::Io {
            msg: err.to_string(),
        }
    }
}

impl From<serde_json::Error> for FerroError {
    fn from(err: serde_json::Error) -> Self {
        FerroError::Json {
            msg: err.to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ErrorCode tests
    #[test]
    fn test_error_code_as_str() {
        assert_eq!(ErrorCode::InvalidAccession.as_str(), "E1001");
        assert_eq!(ErrorCode::UnknownVariantType.as_str(), "E1002");
        assert_eq!(ErrorCode::ReferenceNotFound.as_str(), "E2001");
        assert_eq!(ErrorCode::PositionOutOfBounds.as_str(), "E3001");
        assert_eq!(ErrorCode::IntronicVariant.as_str(), "E4001");
        assert_eq!(ErrorCode::ConversionFailed.as_str(), "E5001");
        assert_eq!(ErrorCode::IoError.as_str(), "E9001");
    }

    #[test]
    fn test_error_code_description() {
        assert_eq!(
            ErrorCode::InvalidAccession.description(),
            "invalid accession format"
        );
        assert_eq!(
            ErrorCode::InvalidPosition.description(),
            "invalid position format"
        );
        assert_eq!(ErrorCode::InvalidEdit.description(), "invalid edit format");
        assert_eq!(
            ErrorCode::UnexpectedEnd.description(),
            "unexpected end of input"
        );
        assert_eq!(
            ErrorCode::UnexpectedChar.description(),
            "unexpected character"
        );
        assert_eq!(
            ErrorCode::InvalidBase.description(),
            "invalid nucleotide base"
        );
        assert_eq!(
            ErrorCode::InvalidAminoAcid.description(),
            "invalid amino acid"
        );
        assert_eq!(
            ErrorCode::SequenceNotFound.description(),
            "sequence not available"
        );
        assert_eq!(
            ErrorCode::ChromosomeNotFound.description(),
            "chromosome not found"
        );
        assert_eq!(
            ErrorCode::ReferenceMismatch.description(),
            "reference sequence mismatch"
        );
        assert_eq!(
            ErrorCode::InvalidRange.description(),
            "invalid coordinate range"
        );
        assert_eq!(
            ErrorCode::ExonIntronBoundary.description(),
            "variant crosses exon-intron boundary"
        );
        assert_eq!(
            ErrorCode::UtrCdsBoundary.description(),
            "variant crosses UTR-CDS boundary"
        );
        assert_eq!(
            ErrorCode::UnsupportedVariant.description(),
            "unsupported variant type"
        );
        assert_eq!(
            ErrorCode::NoOverlappingTranscript.description(),
            "no overlapping transcript"
        );
        assert_eq!(ErrorCode::JsonError.description(), "JSON parsing error");
    }

    #[test]
    fn test_error_code_display() {
        assert_eq!(format!("{}", ErrorCode::InvalidAccession), "E1001");
        assert_eq!(format!("{}", ErrorCode::IoError), "E9001");
    }

    // SourceSpan tests
    #[test]
    fn test_source_span_new() {
        let span = SourceSpan::new(5, 10);
        assert_eq!(span.start, 5);
        assert_eq!(span.end, 10);
    }

    #[test]
    fn test_source_span_point() {
        let span = SourceSpan::point(7);
        assert_eq!(span.start, 7);
        assert_eq!(span.end, 8);
    }

    #[test]
    fn test_source_span_highlight() {
        let span = SourceSpan::new(12, 15);
        let result = span.highlight("NM_000088.3:c.459A>G");
        assert!(result.contains("NM_000088.3:c.459A>G"));
        assert!(result.contains("^"));
    }

    #[test]
    fn test_source_span_highlight_empty_source() {
        let span = SourceSpan::new(0, 5);
        let result = span.highlight("");
        assert_eq!(result, "");
    }

    #[test]
    fn test_source_span_highlight_single_position() {
        let span = SourceSpan::new(0, 1);
        let result = span.highlight("ABC");
        assert!(result.contains("ABC"));
        assert!(result.contains("^"));
    }

    #[test]
    fn test_source_span_highlight_out_of_bounds() {
        let span = SourceSpan::new(100, 200);
        let result = span.highlight("short");
        assert!(result.contains("short"));
    }

    // Diagnostic tests
    #[test]
    fn test_diagnostic_new() {
        let diag = Diagnostic::new();
        assert!(diag.code.is_none());
        assert!(diag.span.is_none());
        assert!(diag.source.is_none());
        assert!(diag.hint.is_none());
        assert!(diag.suggestion.is_none());
    }

    #[test]
    fn test_diagnostic_builder() {
        let diag = Diagnostic::new()
            .with_code(ErrorCode::InvalidAccession)
            .with_span(SourceSpan::new(0, 5))
            .with_source("NM_invalid:c.100A>G")
            .with_hint("Check accession format")
            .with_suggestion("NM_000088.3:c.100A>G");

        assert_eq!(diag.code, Some(ErrorCode::InvalidAccession));
        assert!(diag.span.is_some());
        assert_eq!(diag.source, Some("NM_invalid:c.100A>G".to_string()));
        assert_eq!(diag.hint, Some("Check accession format".to_string()));
        assert_eq!(diag.suggestion, Some("NM_000088.3:c.100A>G".to_string()));
    }

    #[test]
    fn test_diagnostic_format_simple() {
        let diag = Diagnostic::new();
        let result = diag.format("Simple error");
        assert_eq!(result, "Simple error");
    }

    #[test]
    fn test_diagnostic_format_with_code() {
        let diag = Diagnostic::new().with_code(ErrorCode::InvalidAccession);
        let result = diag.format("Invalid format");
        assert!(result.starts_with("[E1001]"));
    }

    #[test]
    fn test_diagnostic_format_full() {
        let diag = Diagnostic::new()
            .with_code(ErrorCode::InvalidAccession)
            .with_span(SourceSpan::new(0, 5))
            .with_source("error_input")
            .with_hint("Try this")
            .with_suggestion("correct_input");

        let result = diag.format("Test message");
        assert!(result.contains("[E1001]"));
        assert!(result.contains("Test message"));
        assert!(result.contains("Hint: Try this"));
        assert!(result.contains("Did you mean: correct_input"));
    }

    // FerroError tests
    #[test]
    fn test_ferro_error_parse() {
        let err = FerroError::parse(10, "unexpected token");
        assert!(matches!(err, FerroError::Parse { pos: 10, .. }));
    }

    #[test]
    fn test_ferro_error_parse_with_diagnostic() {
        let diag = Diagnostic::new().with_code(ErrorCode::InvalidEdit);
        let err = FerroError::parse_with_diagnostic(5, "bad edit", diag);
        assert!(matches!(err, FerroError::Parse { pos: 5, .. }));
    }

    #[test]
    fn test_ferro_error_code() {
        let err = FerroError::ReferenceNotFound {
            id: "NM_123".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::ReferenceNotFound));

        let err = FerroError::ExonIntronBoundary {
            exon: 3,
            variant: "c.100del".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::ExonIntronBoundary));

        let err = FerroError::UtrCdsBoundary {
            variant: "c.1-100del".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::UtrCdsBoundary));

        let err = FerroError::InvalidCoordinates {
            msg: "bad coords".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::InvalidRange));

        let err = FerroError::UnsupportedVariant {
            variant_type: "xyz".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::UnsupportedVariant));

        let err = FerroError::IntronicVariant {
            variant: "c.100+5del".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::IntronicVariant));

        let err = FerroError::ReferenceMismatch {
            location: "c.100".to_string(),
            expected: "A".to_string(),
            found: "G".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::ReferenceMismatch));

        let err = FerroError::ConversionError {
            msg: "conversion failed".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::ConversionFailed));

        let err = FerroError::Io {
            msg: "file error".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::IoError));

        let err = FerroError::Json {
            msg: "json error".to_string(),
        };
        assert_eq!(err.code(), Some(ErrorCode::JsonError));
    }

    #[test]
    fn test_ferro_error_code_from_diagnostic() {
        let diag = Diagnostic::new().with_code(ErrorCode::InvalidPosition);
        let err = FerroError::parse_with_diagnostic(0, "test", diag);
        assert_eq!(err.code(), Some(ErrorCode::InvalidPosition));
    }

    #[test]
    fn test_ferro_error_code_none() {
        let err = FerroError::parse(0, "simple error");
        assert_eq!(err.code(), None);
    }

    #[test]
    fn test_ferro_error_detailed_message() {
        let diag = Diagnostic::new()
            .with_code(ErrorCode::InvalidAccession)
            .with_source("bad_input")
            .with_span(SourceSpan::new(0, 3));
        let err = FerroError::parse_with_diagnostic(0, "invalid", diag);
        let msg = err.detailed_message();
        assert!(msg.contains("[E1001]"));
        assert!(msg.contains("bad_input"));
    }

    #[test]
    fn test_ferro_error_detailed_message_simple() {
        let err = FerroError::ReferenceNotFound {
            id: "NM_test".to_string(),
        };
        let msg = err.detailed_message();
        assert!(msg.contains("NM_test"));
    }

    #[test]
    fn test_ferro_error_display() {
        let err = FerroError::parse(10, "unexpected");
        let display = format!("{}", err);
        assert!(display.contains("10"));
        assert!(display.contains("unexpected"));

        let err = FerroError::ReferenceNotFound {
            id: "test_id".to_string(),
        };
        let display = format!("{}", err);
        assert!(display.contains("test_id"));
    }

    // Suggestion helper tests
    #[test]
    fn test_suggest_variant_type() {
        assert_eq!(suggest_variant_type("c:100A>G"), Some("c."));
        assert_eq!(suggest_variant_type("g:123A>G"), Some("g."));
        assert_eq!(suggest_variant_type("p:Val600Glu"), Some("p."));
        assert_eq!(suggest_variant_type("c100A>G"), Some("c."));
        assert_eq!(suggest_variant_type("C.100A>G"), Some("c."));
        assert_eq!(suggest_variant_type("xyz"), None);
    }

    #[test]
    fn test_suggest_amino_acid() {
        assert_eq!(suggest_amino_acid("Ale"), Some("Ala"));
        assert_eq!(suggest_amino_acid("arg"), Some("Arg"));
        assert_eq!(suggest_amino_acid("asg"), Some("Asn"));
        assert_eq!(suggest_amino_acid("cis"), Some("Cys"));
        assert_eq!(suggest_amino_acid("stop"), Some("Ter"));
        assert_eq!(suggest_amino_acid("end"), Some("Ter"));
        assert_eq!(suggest_amino_acid("unk"), Some("Xaa"));
        assert_eq!(suggest_amino_acid("sel"), Some("Sec"));
        assert_eq!(suggest_amino_acid("xyz"), None);
    }

    // From impl tests
    #[test]
    fn test_from_io_error() {
        let io_err = std::io::Error::new(std::io::ErrorKind::NotFound, "file not found");
        let ferro_err: FerroError = io_err.into();
        assert!(matches!(ferro_err, FerroError::Io { .. }));
        assert!(ferro_err.to_string().contains("not found"));
    }

    #[test]
    fn test_ferro_error_equality() {
        let err1 = FerroError::parse(10, "test");
        let err2 = FerroError::parse(10, "test");
        assert_eq!(err1, err2);

        let err3 = FerroError::parse(11, "test");
        assert_ne!(err1, err3);
    }

    #[test]
    fn test_source_span_equality() {
        let span1 = SourceSpan::new(5, 10);
        let span2 = SourceSpan::new(5, 10);
        assert_eq!(span1, span2);

        let span3 = SourceSpan::new(5, 11);
        assert_ne!(span1, span3);
    }

    #[test]
    fn test_diagnostic_equality() {
        let diag1 = Diagnostic::new().with_code(ErrorCode::InvalidAccession);
        let diag2 = Diagnostic::new().with_code(ErrorCode::InvalidAccession);
        assert_eq!(diag1, diag2);

        let diag3 = Diagnostic::new().with_code(ErrorCode::InvalidEdit);
        assert_ne!(diag1, diag3);
    }

    #[test]
    fn test_error_code_hash() {
        use std::collections::HashSet;
        let mut set = HashSet::new();
        set.insert(ErrorCode::InvalidAccession);
        set.insert(ErrorCode::InvalidEdit);
        assert!(set.contains(&ErrorCode::InvalidAccession));
        assert!(set.contains(&ErrorCode::InvalidEdit));
        assert!(!set.contains(&ErrorCode::InvalidPosition));
    }
}
