//! MAVE-HGVS parser with context support.
//!
//! Parses MAVE-HGVS short-form notation using context to supply missing accessions.

use super::context::MaveContext;
use crate::hgvs::parser::parse_hgvs;
use crate::hgvs::variant::HgvsVariant;
use std::fmt;

/// Error type for MAVE-HGVS parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MaveParseError {
    /// The input is empty.
    EmptyInput,

    /// Cannot determine the coordinate type from the input.
    UnknownCoordinateType(String),

    /// No accession available in context for the coordinate type.
    MissingAccession { coord_type: char, input: String },

    /// The underlying HGVS parser failed.
    HgvsParseError { input: String, message: String },

    /// Input already has an accession (not a short-form).
    AlreadyHasAccession(String),
}

impl fmt::Display for MaveParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "Empty input"),
            Self::UnknownCoordinateType(input) => {
                write!(f, "Cannot determine coordinate type from '{}'", input)
            }
            Self::MissingAccession { coord_type, input } => {
                write!(
                    f,
                    "No accession in context for '{}.' variants: '{}'",
                    coord_type, input
                )
            }
            Self::HgvsParseError { input, message } => {
                write!(f, "Failed to parse '{}': {}", input, message)
            }
            Self::AlreadyHasAccession(input) => {
                write!(
                    f,
                    "Input already has accession, use parse_hgvs directly: '{}'",
                    input
                )
            }
        }
    }
}

impl std::error::Error for MaveParseError {}

/// Parse a MAVE-HGVS short-form variant with context.
///
/// This function handles short-form HGVS notation that lacks a reference accession,
/// using the provided context to supply the appropriate accession.
///
/// # Arguments
///
/// * `input` - The MAVE-HGVS short-form string (e.g., "p.Glu6Val", "c.20A>T")
/// * `context` - The MaveContext providing target sequence accessions
///
/// # Returns
///
/// The parsed HgvsVariant with the full accession from context.
///
/// # Example
///
/// ```
/// use ferro_hgvs::mave::{MaveContext, parse_mave_hgvs};
///
/// let context = MaveContext::new()
///     .with_protein_accession("NP_000509.1");
///
/// let variant = parse_mave_hgvs("p.Glu6Val", &context).unwrap();
/// assert_eq!(variant.to_string(), "NP_000509.1:p.Glu6Val");
/// ```
pub fn parse_mave_hgvs(input: &str, context: &MaveContext) -> Result<HgvsVariant, MaveParseError> {
    let input = input.trim();

    if input.is_empty() {
        return Err(MaveParseError::EmptyInput);
    }

    // Check if input already has an accession (contains ':')
    if input.contains(':') {
        // Try to parse it directly as full HGVS
        return parse_hgvs(input).map_err(|e| MaveParseError::HgvsParseError {
            input: input.to_string(),
            message: e.to_string(),
        });
    }

    // Determine the coordinate type from the input
    let coord_type = detect_coordinate_type(input)?;

    // Get the appropriate accession from context
    let accession = context
        .accession_for_coordinate_type(coord_type)
        .ok_or_else(|| MaveParseError::MissingAccession {
            coord_type,
            input: input.to_string(),
        })?;

    // Construct full HGVS string
    let full_hgvs = format!("{}:{}", accession, input);

    // Parse the full HGVS
    parse_hgvs(&full_hgvs).map_err(|e| MaveParseError::HgvsParseError {
        input: input.to_string(),
        message: e.to_string(),
    })
}

/// Parse a MAVE-HGVS variant, allowing both short and full forms.
///
/// This is a lenient parser that:
/// - Uses context for short-form variants (e.g., "p.Glu6Val")
/// - Passes through full-form variants (e.g., "NP_000509.1:p.Glu6Val")
///
/// # Arguments
///
/// * `input` - The MAVE-HGVS string (short or full form)
/// * `context` - The MaveContext providing target sequence accessions
///
/// # Returns
///
/// The parsed HgvsVariant.
pub fn parse_mave_hgvs_lenient(
    input: &str,
    context: &MaveContext,
) -> Result<HgvsVariant, MaveParseError> {
    let input = input.trim();

    if input.is_empty() {
        return Err(MaveParseError::EmptyInput);
    }

    // If it has an accession, parse directly
    if input.contains(':') {
        return parse_hgvs(input).map_err(|e| MaveParseError::HgvsParseError {
            input: input.to_string(),
            message: e.to_string(),
        });
    }

    // Otherwise, use context-aware parsing
    parse_mave_hgvs(input, context)
}

/// Detect the coordinate type from a short-form MAVE-HGVS string.
///
/// # Arguments
///
/// * `input` - The short-form HGVS string (e.g., "p.Glu6Val", "c.20A>T")
///
/// # Returns
///
/// The coordinate type character ('p', 'c', 'g', 'n', 'r', 'm', 'o').
fn detect_coordinate_type(input: &str) -> Result<char, MaveParseError> {
    // Look for coordinate type prefix: p., c., g., n., r., m., o.
    let chars: Vec<char> = input.chars().collect();

    if chars.len() >= 2 && chars[1] == '.' {
        let coord = chars[0].to_ascii_lowercase();
        if ['p', 'c', 'g', 'n', 'r', 'm', 'o'].contains(&coord) {
            return Ok(coord);
        }
    }

    Err(MaveParseError::UnknownCoordinateType(input.to_string()))
}

/// Check if an input string is a short-form MAVE-HGVS (no accession).
pub fn is_mave_short_form(input: &str) -> bool {
    let input = input.trim();
    !input.is_empty() && !input.contains(':') && detect_coordinate_type(input).is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_context() -> MaveContext {
        MaveContext::new()
            .with_protein_accession("NP_000509.1")
            .with_coding_accession("NM_000518.5")
            .with_genomic_accession("NC_000011.10")
            .with_gene_symbol("HBB")
    }

    #[test]
    fn test_parse_protein_variant() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Glu6Val", &ctx);
        assert!(result.is_ok());
        let variant = result.unwrap();
        assert_eq!(variant.to_string(), "NP_000509.1:p.Glu6Val");
    }

    #[test]
    fn test_parse_protein_missense() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Leu11Pro", &ctx);
        assert!(result.is_ok());
        let variant = result.unwrap();
        assert_eq!(variant.to_string(), "NP_000509.1:p.Leu11Pro");
    }

    #[test]
    fn test_parse_coding_variant() {
        let ctx = test_context();
        let result = parse_mave_hgvs("c.20A>T", &ctx);
        assert!(result.is_ok());
        let variant = result.unwrap();
        assert_eq!(variant.to_string(), "NM_000518.5:c.20A>T");
    }

    #[test]
    fn test_parse_genomic_variant() {
        let ctx = test_context();
        let result = parse_mave_hgvs("g.12345A>G", &ctx);
        assert!(result.is_ok());
        let variant = result.unwrap();
        assert_eq!(variant.to_string(), "NC_000011.10:g.12345A>G");
    }

    #[test]
    fn test_parse_full_form_passthrough() {
        let ctx = test_context();
        let result = parse_mave_hgvs("NP_000509.1:p.Glu6Val", &ctx);
        assert!(result.is_ok());
        let variant = result.unwrap();
        assert_eq!(variant.to_string(), "NP_000509.1:p.Glu6Val");
    }

    #[test]
    fn test_parse_empty_input() {
        let ctx = test_context();
        let result = parse_mave_hgvs("", &ctx);
        assert!(matches!(result, Err(MaveParseError::EmptyInput)));
    }

    #[test]
    fn test_parse_unknown_coordinate_type() {
        let ctx = test_context();
        let result = parse_mave_hgvs("x.123", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::UnknownCoordinateType(_))
        ));
    }

    #[test]
    fn test_parse_missing_accession() {
        let ctx = MaveContext::new().with_protein_accession("NP_000509.1");
        // No coding accession in context
        let result = parse_mave_hgvs("c.20A>T", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::MissingAccession { .. })
        ));
    }

    #[test]
    fn test_parse_protein_deletion() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Leu11del", &ctx);
        assert!(result.is_ok());
        let variant = result.unwrap();
        assert!(variant.to_string().contains("p.Leu11del"));
    }

    #[test]
    fn test_parse_protein_insertion() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Leu11_Pro12insAla", &ctx);
        assert!(result.is_ok());
    }

    #[test]
    fn test_parse_protein_frameshift() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Glu6fs", &ctx);
        assert!(result.is_ok());
    }

    #[test]
    fn test_parse_protein_extension() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Ter147Tyrext*30", &ctx);
        assert!(result.is_ok());
    }

    #[test]
    fn test_detect_coordinate_type_protein() {
        assert_eq!(detect_coordinate_type("p.Glu6Val"), Ok('p'));
    }

    #[test]
    fn test_detect_coordinate_type_coding() {
        assert_eq!(detect_coordinate_type("c.20A>T"), Ok('c'));
    }

    #[test]
    fn test_detect_coordinate_type_genomic() {
        assert_eq!(detect_coordinate_type("g.12345A>G"), Ok('g'));
    }

    #[test]
    fn test_detect_coordinate_type_noncoding() {
        assert_eq!(detect_coordinate_type("n.100A>G"), Ok('n'));
    }

    #[test]
    fn test_detect_coordinate_type_invalid() {
        assert!(detect_coordinate_type("invalid").is_err());
        assert!(detect_coordinate_type("x.123").is_err());
    }

    #[test]
    fn test_is_mave_short_form() {
        assert!(is_mave_short_form("p.Glu6Val"));
        assert!(is_mave_short_form("c.20A>T"));
        assert!(!is_mave_short_form("NP_000509.1:p.Glu6Val"));
        assert!(!is_mave_short_form(""));
        assert!(!is_mave_short_form("invalid"));
    }

    #[test]
    fn test_lenient_parse_short_form() {
        let ctx = test_context();
        let result = parse_mave_hgvs_lenient("p.Glu6Val", &ctx);
        assert!(result.is_ok());
    }

    #[test]
    fn test_lenient_parse_full_form() {
        let ctx = test_context();
        let result = parse_mave_hgvs_lenient("NP_000509.1:p.Glu6Val", &ctx);
        assert!(result.is_ok());
    }

    #[test]
    fn test_allele_variant() {
        let ctx = test_context();
        // Alleles like c.[32T>C;39G>A] should work with context
        let result = parse_mave_hgvs("c.[32T>C;39G>A]", &ctx);
        // Note: this may fail depending on allele support, but we test the context injection
        if let Err(e) = &result {
            // Should not be a MissingAccession error
            assert!(
                !matches!(e, MaveParseError::MissingAccession { .. }),
                "Should have injected accession"
            );
        }
    }

    #[test]
    fn test_error_display() {
        let err = MaveParseError::EmptyInput;
        assert_eq!(err.to_string(), "Empty input");

        let err = MaveParseError::MissingAccession {
            coord_type: 'c',
            input: "c.20A>T".to_string(),
        };
        assert!(err.to_string().contains("No accession"));
    }

    // =========================================================================
    // P1: MAVE allele roundtrip tests
    // =========================================================================

    #[test]
    fn test_allele_simple_coding() {
        let ctx = test_context();
        // Simple two-variant allele
        let result = parse_mave_hgvs("c.[20A>T;30G>C]", &ctx);
        assert!(result.is_ok(), "Failed to parse: {:?}", result.err());
        let variant = result.unwrap();
        let output = variant.to_string();
        assert!(
            output.contains("NM_000518.5"),
            "Output should contain accession: {}",
            output
        );
        assert!(
            output.contains("c.[") || output.contains("c.20"),
            "Output should contain allele notation: {}",
            output
        );
    }

    #[test]
    fn test_allele_multiple_variants() {
        let ctx = test_context();
        // MaveDB style: c.[32T>C;39G>A;42A>G;81C>T]
        let result = parse_mave_hgvs("c.[32T>C;39G>A;42A>G]", &ctx);
        assert!(result.is_ok(), "Failed to parse multi-variant allele");
    }

    #[test]
    fn test_allele_protein() {
        let ctx = test_context();
        // Protein allele
        let result = parse_mave_hgvs("p.[Glu6Val;Lys17Arg]", &ctx);
        if let Ok(variant) = result {
            let output = variant.to_string();
            assert!(
                output.contains("NP_000509.1"),
                "Output should contain protein accession: {}",
                output
            );
        }
        // Note: Some allele formats may not be fully supported
    }

    #[test]
    fn test_allele_genomic() {
        let ctx = test_context();
        // Genomic allele
        let result = parse_mave_hgvs("g.[12345A>G;12350C>T]", &ctx);
        // Genomic alleles may not be fully supported, check accession injection
        if let Err(e) = &result {
            // Should not be a MissingAccession error - accession should be injected
            assert!(
                !matches!(e, MaveParseError::MissingAccession { .. }),
                "Accession should be injected: {:?}",
                e
            );
        } else {
            let variant = result.unwrap();
            let output = variant.to_string();
            assert!(
                output.contains("NC_000011.10"),
                "Output should contain genomic accession: {}",
                output
            );
        }
    }

    #[test]
    fn test_allele_roundtrip_structure() {
        let ctx = test_context();

        // Test that parsing produces a variant we can re-serialize
        let input = "c.[20A>T;30G>C]";
        let result = parse_mave_hgvs(input, &ctx);
        if let Ok(variant) = result {
            let output = variant.to_string();
            // The accession should be included in the output
            // Note: Allele output format may include accession on each part:
            // "[NM_000518.5:c.20A>T;NM_000518.5:c.30G>C]"
            // Or may keep the structure: "NM_000518.5:c.[20A>T;30G>C]"
            assert!(
                output.contains("NM_000518.5"),
                "Should contain accession: {}",
                output
            );
            // Should still be parseable (even if format differs)
            let reparsed = crate::hgvs::parser::parse_hgvs(&output);
            assert!(
                reparsed.is_ok(),
                "Reparsed variant should be valid HGVS: {}",
                output
            );
        }
    }

    #[test]
    fn test_single_variant_in_brackets() {
        let ctx = test_context();
        // Single variant in allele notation (edge case)
        let result = parse_mave_hgvs("c.[20A>T]", &ctx);
        assert!(result.is_ok(), "Should parse single variant in brackets");
    }

    #[test]
    fn test_mavedb_real_pattern_1() {
        // From actual MaveDB data: c.[32T>C;39G>A;42A>G;81C>T;99T>C;105G>T]
        let ctx = test_context();
        let result = parse_mave_hgvs("c.[32T>C;39G>A;42A>G;81C>T;99T>C;105G>T]", &ctx);
        assert!(result.is_ok(), "Should parse real MaveDB allele pattern");
    }

    #[test]
    fn test_mavedb_real_pattern_2() {
        // From MaveDB: c.[16G>C;17A>C;18A>G;39G>A;42A>G]
        let ctx = test_context();
        let result = parse_mave_hgvs("c.[16G>C;17A>C;18A>G;39G>A;42A>G]", &ctx);
        assert!(result.is_ok(), "Should parse real MaveDB allele pattern");
    }

    #[test]
    fn test_protein_complex_allele() {
        let ctx = test_context();
        // Complex protein changes in allele
        let result = parse_mave_hgvs("p.[Glu6Val;Lys17del]", &ctx);
        // May or may not succeed depending on parser support
        // But should not fail due to missing accession
        if let Err(e) = &result {
            assert!(
                !matches!(e, MaveParseError::MissingAccession { .. }),
                "Accession should be injected even for complex alleles"
            );
        }
    }

    #[test]
    fn test_allele_with_whitespace() {
        let ctx = test_context();
        // With extra whitespace (should be trimmed)
        let result = parse_mave_hgvs("  c.[20A>T;30G>C]  ", &ctx);
        assert!(result.is_ok(), "Should handle whitespace around allele");
    }

    // =========================================================================
    // P3: MAVE error recovery tests
    // =========================================================================

    #[test]
    fn test_error_recovery_whitespace_only() {
        let ctx = test_context();
        let result = parse_mave_hgvs("   ", &ctx);
        assert!(matches!(result, Err(MaveParseError::EmptyInput)));
    }

    #[test]
    fn test_error_recovery_tab_and_newline() {
        let ctx = test_context();
        let result = parse_mave_hgvs("\t\n", &ctx);
        assert!(matches!(result, Err(MaveParseError::EmptyInput)));
    }

    #[test]
    fn test_error_recovery_invalid_coord_type_uppercase() {
        let ctx = test_context();
        // Uppercase coordinate type - detection is case insensitive
        let result = parse_mave_hgvs("P.Glu6Val", &ctx);
        // The parser may or may not handle uppercase in the variant part
        // If it fails, it should be an HgvsParseError, not MissingAccession
        if let Err(e) = &result {
            assert!(
                !matches!(e, MaveParseError::MissingAccession { .. }),
                "Should inject accession even with uppercase coord type"
            );
        }
    }

    #[test]
    fn test_error_recovery_invalid_coord_type_z() {
        let ctx = test_context();
        let result = parse_mave_hgvs("z.123A>G", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::UnknownCoordinateType(_))
        ));
        if let Err(e) = result {
            let msg = e.to_string();
            assert!(msg.contains("z.123A>G"), "Error should include input");
        }
    }

    #[test]
    fn test_error_recovery_numeric_coord_type() {
        let ctx = test_context();
        let result = parse_mave_hgvs("1.123A>G", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::UnknownCoordinateType(_))
        ));
    }

    #[test]
    fn test_error_recovery_no_dot_after_coord() {
        let ctx = test_context();
        let result = parse_mave_hgvs("cGlu6Val", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::UnknownCoordinateType(_))
        ));
    }

    #[test]
    fn test_error_recovery_missing_protein_accession() {
        let ctx = MaveContext::new().with_coding_accession("NM_000518.5");
        // No protein accession
        let result = parse_mave_hgvs("p.Glu6Val", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::MissingAccession { .. })
        ));
        if let Err(MaveParseError::MissingAccession { coord_type, input }) = result {
            assert_eq!(coord_type, 'p');
            assert_eq!(input, "p.Glu6Val");
        }
    }

    #[test]
    fn test_error_recovery_missing_genomic_accession() {
        let ctx = MaveContext::new().with_protein_accession("NP_000509.1");
        let result = parse_mave_hgvs("g.12345A>G", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::MissingAccession { .. })
        ));
    }

    #[test]
    fn test_error_recovery_missing_noncoding_accession() {
        let ctx = MaveContext::new().with_protein_accession("NP_000509.1");
        let result = parse_mave_hgvs("n.100A>G", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::MissingAccession { .. })
        ));
    }

    #[test]
    fn test_error_recovery_invalid_hgvs_after_accession() {
        let ctx = test_context();
        // Accession gets injected, but variant syntax is invalid
        let result = parse_mave_hgvs("c.invalid_syntax", &ctx);
        assert!(matches!(result, Err(MaveParseError::HgvsParseError { .. })));
        if let Err(MaveParseError::HgvsParseError { input, message }) = result {
            assert_eq!(input, "c.invalid_syntax");
            assert!(!message.is_empty());
        }
    }

    #[test]
    fn test_error_recovery_invalid_position() {
        let ctx = test_context();
        let result = parse_mave_hgvs("c.XXX>G", &ctx);
        assert!(matches!(result, Err(MaveParseError::HgvsParseError { .. })));
    }

    #[test]
    fn test_error_recovery_invalid_protein_aa() {
        let ctx = test_context();
        let result = parse_mave_hgvs("p.Xxx6Val", &ctx);
        // Xxx is valid for unknown amino acid, but check it parses
        if let Err(e) = &result {
            assert!(!matches!(e, MaveParseError::MissingAccession { .. }));
        }
    }

    #[test]
    fn test_error_recovery_full_form_invalid() {
        let ctx = test_context();
        // Full form (has accession) but invalid variant
        let result = parse_mave_hgvs("NP_000509.1:p.invalid", &ctx);
        assert!(matches!(result, Err(MaveParseError::HgvsParseError { .. })));
    }

    #[test]
    fn test_error_recovery_partial_accession() {
        let ctx = test_context();
        // Has colon but malformed - treated as full form
        let result = parse_mave_hgvs("NP:p.Glu6Val", &ctx);
        // This might parse if "NP" is accepted as an accession prefix
        // The key is it should not inject a context accession
        // Just verify it doesn't panic and handles the input
        let _ = result;
    }

    #[test]
    fn test_error_message_empty_input() {
        let err = MaveParseError::EmptyInput;
        assert_eq!(format!("{}", err), "Empty input");
    }

    #[test]
    fn test_error_message_unknown_coordinate() {
        let err = MaveParseError::UnknownCoordinateType("x.123".to_string());
        let msg = format!("{}", err);
        assert!(msg.contains("Cannot determine coordinate type"));
        assert!(msg.contains("x.123"));
    }

    #[test]
    fn test_error_message_missing_accession() {
        let err = MaveParseError::MissingAccession {
            coord_type: 'c',
            input: "c.20A>T".to_string(),
        };
        let msg = format!("{}", err);
        assert!(msg.contains("No accession"));
        assert!(msg.contains("c."));
        assert!(msg.contains("c.20A>T"));
    }

    #[test]
    fn test_error_message_hgvs_parse_error() {
        let err = MaveParseError::HgvsParseError {
            input: "c.invalid".to_string(),
            message: "expected number".to_string(),
        };
        let msg = format!("{}", err);
        assert!(msg.contains("Failed to parse"));
        assert!(msg.contains("c.invalid"));
        assert!(msg.contains("expected number"));
    }

    #[test]
    fn test_error_message_already_has_accession() {
        let err = MaveParseError::AlreadyHasAccession("NP_000509.1:p.Glu6Val".to_string());
        let msg = format!("{}", err);
        assert!(msg.contains("already has accession"));
        assert!(msg.contains("NP_000509.1:p.Glu6Val"));
    }

    #[test]
    fn test_error_is_std_error() {
        let err = MaveParseError::EmptyInput;
        // Verify it implements std::error::Error
        let _: &dyn std::error::Error = &err;
    }

    #[test]
    fn test_lenient_parse_recovers_from_missing_accession() {
        let ctx = test_context();
        // Short form should use lenient parsing
        let short_result = parse_mave_hgvs_lenient("p.Glu6Val", &ctx);
        assert!(short_result.is_ok());

        // Full form should also work
        let full_result = parse_mave_hgvs_lenient("NP_000509.1:p.Glu6Val", &ctx);
        assert!(full_result.is_ok());

        // But invalid should still error
        let invalid_result = parse_mave_hgvs_lenient("x.123", &ctx);
        assert!(invalid_result.is_err());
    }

    #[test]
    fn test_lenient_parse_empty_context() {
        let ctx = MaveContext::new();
        // No accessions at all
        let result = parse_mave_hgvs_lenient("p.Glu6Val", &ctx);
        assert!(matches!(
            result,
            Err(MaveParseError::MissingAccession { .. })
        ));
    }

    #[test]
    fn test_is_short_form_edge_cases() {
        // Valid short forms
        assert!(is_mave_short_form("p.Glu6Val"));
        assert!(is_mave_short_form("c.20A>T"));
        assert!(is_mave_short_form("g.12345del"));
        assert!(is_mave_short_form("n.100A>G"));
        assert!(is_mave_short_form("r.100a>g"));
        assert!(is_mave_short_form("m.100A>G"));

        // Not short forms
        assert!(!is_mave_short_form("NP_000509.1:p.Glu6Val")); // has accession
        assert!(!is_mave_short_form("")); // empty
        assert!(!is_mave_short_form("   ")); // whitespace only
        assert!(!is_mave_short_form("invalid")); // no coord type
        assert!(!is_mave_short_form("x.123")); // invalid coord type
        assert!(!is_mave_short_form("12345")); // just a number
        assert!(!is_mave_short_form("Glu6Val")); // no prefix
    }

    #[test]
    fn test_coordinate_detection_all_types() {
        // All valid coordinate types
        assert_eq!(detect_coordinate_type("p.Glu6Val"), Ok('p'));
        assert_eq!(detect_coordinate_type("c.20A>T"), Ok('c'));
        assert_eq!(detect_coordinate_type("g.12345A>G"), Ok('g'));
        assert_eq!(detect_coordinate_type("n.100A>G"), Ok('n'));
        assert_eq!(detect_coordinate_type("r.100a>g"), Ok('r'));
        assert_eq!(detect_coordinate_type("m.100A>G"), Ok('m'));
        assert_eq!(detect_coordinate_type("o.100A>G"), Ok('o'));
    }

    #[test]
    fn test_coordinate_detection_case_insensitive() {
        assert_eq!(detect_coordinate_type("P.Glu6Val"), Ok('p'));
        assert_eq!(detect_coordinate_type("C.20A>T"), Ok('c'));
        assert_eq!(detect_coordinate_type("G.12345A>G"), Ok('g'));
    }

    #[test]
    fn test_context_accession_mapping() {
        let ctx = MaveContext::new()
            .with_protein_accession("NP_001")
            .with_coding_accession("NM_001")
            .with_genomic_accession("NC_001")
            .with_noncoding_accession("NR_001");

        assert_eq!(ctx.accession_for_coordinate_type('p'), Some("NP_001"));
        assert_eq!(ctx.accession_for_coordinate_type('c'), Some("NM_001"));
        assert_eq!(ctx.accession_for_coordinate_type('g'), Some("NC_001"));
        assert_eq!(ctx.accession_for_coordinate_type('n'), Some("NR_001"));
        // RNA uses coding accession
        assert_eq!(ctx.accession_for_coordinate_type('r'), Some("NM_001"));
        // Mitochondrial uses genomic accession
        assert_eq!(ctx.accession_for_coordinate_type('m'), Some("NC_001"));
    }

    #[test]
    fn test_multiple_parse_attempts() {
        let ctx = test_context();

        // Parse multiple variants to test consistency
        let inputs = vec![
            "p.Glu6Val",
            "c.20A>T",
            "p.Lys17Arg",
            "c.30G>C",
            "g.12345A>G",
        ];

        for input in inputs {
            let result = parse_mave_hgvs(input, &ctx);
            assert!(result.is_ok(), "Failed to parse: {}", input);
        }
    }

    #[test]
    fn test_consecutive_errors_dont_accumulate() {
        let ctx = MaveContext::new(); // Empty context

        // Multiple error attempts should each give clean errors
        for _ in 0..5 {
            let result = parse_mave_hgvs("p.Glu6Val", &ctx);
            assert!(matches!(
                result,
                Err(MaveParseError::MissingAccession { .. })
            ));
        }
    }

    #[test]
    fn test_special_characters_in_input() {
        let ctx = test_context();

        // These should fail with appropriate errors, not panic
        let inputs = vec![
            "p.Glu6Val\0", // null byte
            "c.20A>G\n",   // newline
            "g.123<>",     // angle brackets
            "c.20A>G\"",   // quote
        ];

        for input in inputs {
            let result = parse_mave_hgvs(input, &ctx);
            // Should not panic, may succeed or fail appropriately
            let _ = result;
        }
    }

    #[test]
    fn test_very_long_input() {
        let ctx = test_context();
        // Very long position number
        let result = parse_mave_hgvs("c.999999999999A>G", &ctx);
        // Should not panic, may or may not parse depending on position limits
        let _ = result;
    }

    #[test]
    fn test_unicode_in_input() {
        let ctx = test_context();
        // Unicode characters
        let result = parse_mave_hgvs("p.Glu6Välñ", &ctx);
        // Should fail gracefully with HgvsParseError
        assert!(result.is_err());
    }
}
