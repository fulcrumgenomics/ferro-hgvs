//! Reference sequence validation for normalization.
//!
//! Validates that stated reference bases in HGVS expressions match
//! the actual reference sequence.
//!
//! # Coordinate System
//!
//! | Parameter | Basis | Notes |
//! |-----------|-------|-------|
//! | `start`, `end` | 1-based | HGVS positions (inclusive) |
//! | Array indexing | 0-based | Converted via `hgvs_pos_to_index()` |
//!
//! For type-safe coordinate handling, see [`crate::coords`].

use crate::coords::hgvs_pos_to_index;
use crate::error::FerroError;
use crate::error_handling::ResolvedAction;
use crate::hgvs::edit::{Base, NaEdit};
use crate::normalize::config::NormalizeConfig;

/// Result of reference validation
#[derive(Debug, Clone)]
pub struct ValidationResult {
    /// Whether validation passed
    pub valid: bool,
    /// Warning message if any (for lenient mode)
    pub warning: Option<String>,
    /// The actual reference sequence (for correction)
    pub actual_ref: Option<String>,
    /// The stated reference sequence
    pub stated_ref: Option<String>,
}

impl ValidationResult {
    /// Create a passing validation result
    pub fn ok() -> Self {
        Self {
            valid: true,
            warning: None,
            actual_ref: None,
            stated_ref: None,
        }
    }

    /// Create a mismatch result
    pub fn mismatch(stated: String, actual: String) -> Self {
        Self {
            valid: false,
            warning: Some(format!(
                "Reference mismatch: stated '{}' but actual is '{}'",
                stated, actual
            )),
            actual_ref: Some(actual),
            stated_ref: Some(stated),
        }
    }
}

/// Validate the reference bases in an edit against the actual sequence.
///
/// Returns a ValidationResult indicating whether the stated reference matches
/// the actual sequence, and provides correction information if not.
///
/// # Arguments
/// * `edit` - The edit to validate
/// * `ref_seq` - The reference sequence (1-indexed, so position 1 is at index 0)
/// * `start` - Start position (1-indexed)
/// * `end` - End position (1-indexed, inclusive)
pub fn validate_reference(edit: &NaEdit, ref_seq: &[u8], start: u64, end: u64) -> ValidationResult {
    match edit {
        NaEdit::Substitution { reference, .. } => validate_single_base(reference, ref_seq, start),
        NaEdit::Deletion { sequence, .. } => {
            if let Some(seq) = sequence {
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else {
                // No sequence stated, nothing to validate
                ValidationResult::ok()
            }
        }
        NaEdit::Delins { .. } => {
            // Delins doesn't state the reference in standard HGVS
            // (though some parsers extract it from formats like "delACinsGT")
            ValidationResult::ok()
        }
        NaEdit::Duplication { sequence, .. } => {
            if let Some(seq) = sequence {
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else {
                // No sequence stated, nothing to validate
                ValidationResult::ok()
            }
        }
        NaEdit::Inversion { sequence, .. } => {
            if let Some(seq) = sequence {
                validate_sequence(seq.bases(), ref_seq, start, end)
            } else {
                ValidationResult::ok()
            }
        }
        // Other edit types don't have stated reference bases
        _ => ValidationResult::ok(),
    }
}

/// Validate a single base against the reference
fn validate_single_base(stated: &Base, ref_seq: &[u8], pos: u64) -> ValidationResult {
    // Convert 1-based HGVS position to 0-based array index
    let idx = hgvs_pos_to_index(pos);

    if idx >= ref_seq.len() {
        return ValidationResult::mismatch(
            stated.to_char().to_string(),
            format!("(position {} out of range)", pos),
        );
    }

    let actual_byte = ref_seq[idx];
    let stated_byte = stated.to_u8();

    // Handle case-insensitive comparison
    if actual_byte.eq_ignore_ascii_case(&stated_byte) {
        ValidationResult::ok()
    } else {
        let actual_char = (actual_byte as char).to_ascii_uppercase();
        ValidationResult::mismatch(stated.to_char().to_string(), actual_char.to_string())
    }
}

/// Validate a sequence against the reference
fn validate_sequence(stated: &[Base], ref_seq: &[u8], start: u64, end: u64) -> ValidationResult {
    // Convert 1-based HGVS positions to 0-based array indices
    // start/end are inclusive in HGVS, so end_idx = end for half-open interval
    let start_idx = hgvs_pos_to_index(start);
    let end_idx = end as usize; // 1-based inclusive end = 0-based exclusive end

    if end_idx > ref_seq.len() {
        let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
        return ValidationResult::mismatch(stated_str, format!("(position {} out of range)", end));
    }

    let actual_bytes = &ref_seq[start_idx..end_idx];

    // Compare lengths first
    if stated.len() != actual_bytes.len() {
        let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
        let actual_str: String = actual_bytes
            .iter()
            .map(|&b| (b as char).to_ascii_uppercase())
            .collect();
        return ValidationResult::mismatch(stated_str, actual_str);
    }

    // Compare each base
    for (stated_base, &actual_byte) in stated.iter().zip(actual_bytes.iter()) {
        let stated_byte = stated_base.to_u8();
        if !actual_byte.eq_ignore_ascii_case(&stated_byte) {
            let stated_str: String = stated.iter().map(|b| b.to_char()).collect();
            let actual_str: String = actual_bytes
                .iter()
                .map(|&b| (b as char).to_ascii_uppercase())
                .collect();
            return ValidationResult::mismatch(stated_str, actual_str);
        }
    }

    ValidationResult::ok()
}

/// Apply reference validation policy based on configuration.
///
/// Returns Ok(()) if validation passes or correction is allowed,
/// Returns Err if validation fails and strict mode rejects it.
///
/// Warnings are printed to stderr in lenient mode.
pub fn apply_validation_policy(
    result: &ValidationResult,
    config: &NormalizeConfig,
    variant_str: &str,
) -> Result<(), FerroError> {
    if result.valid {
        return Ok(());
    }

    let action = config.ref_mismatch_action();

    match action {
        ResolvedAction::Reject => Err(FerroError::ReferenceMismatch {
            location: variant_str.to_string(),
            expected: result.stated_ref.clone().unwrap_or_else(|| "?".to_string()),
            found: result.actual_ref.clone().unwrap_or_else(|| "?".to_string()),
        }),
        ResolvedAction::WarnCorrect => {
            // Print warning to stderr
            if let Some(ref warning) = result.warning {
                eprintln!("Warning: {} in '{}'", warning, variant_str);
            }
            Ok(())
        }
        ResolvedAction::SilentCorrect | ResolvedAction::Accept => {
            // Silently proceed
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hgvs::edit::Sequence;
    use std::str::FromStr;

    #[test]
    fn test_validate_substitution_match() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 1);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_substitution_mismatch() {
        let edit = NaEdit::Substitution {
            reference: Base::G, // stated G
            alternative: Base::A,
        };
        let ref_seq = b"ATGC"; // actual is A at position 1
        let result = validate_reference(&edit, ref_seq, 1, 1);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("G".to_string()));
        assert_eq!(result.actual_ref, Some("A".to_string()));
    }

    #[test]
    fn test_validate_deletion_match() {
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("ATG").unwrap()),
            length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid);
    }

    #[test]
    fn test_validate_deletion_mismatch() {
        let edit = NaEdit::Deletion {
            sequence: Some(Sequence::from_str("GGG").unwrap()),
            length: None,
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(!result.valid);
        assert_eq!(result.stated_ref, Some("GGG".to_string()));
        assert_eq!(result.actual_ref, Some("ATG".to_string()));
    }

    #[test]
    fn test_validate_no_sequence() {
        let edit = NaEdit::Deletion {
            sequence: None,
            length: Some(3),
        };
        let ref_seq = b"ATGC";
        let result = validate_reference(&edit, ref_seq, 1, 3);
        assert!(result.valid); // No stated sequence to validate
    }

    #[test]
    fn test_validate_case_insensitive() {
        let edit = NaEdit::Substitution {
            reference: Base::A,
            alternative: Base::G,
        };
        let ref_seq = b"atgc"; // lowercase
        let result = validate_reference(&edit, ref_seq, 1, 1);
        assert!(result.valid);
    }

    #[test]
    fn test_apply_policy_strict() {
        let result = ValidationResult::mismatch("G".to_string(), "A".to_string());
        let config = NormalizeConfig::strict();
        let err = apply_validation_policy(&result, &config, "c.1G>T");
        assert!(err.is_err());
    }

    #[test]
    fn test_apply_policy_lenient() {
        let result = ValidationResult::mismatch("G".to_string(), "A".to_string());
        let config = NormalizeConfig::lenient();
        let ok = apply_validation_policy(&result, &config, "c.1G>T");
        assert!(ok.is_ok()); // Should pass but emit warning
    }

    #[test]
    fn test_apply_policy_silent() {
        let result = ValidationResult::mismatch("G".to_string(), "A".to_string());
        let config = NormalizeConfig::silent();
        let ok = apply_validation_policy(&result, &config, "c.1G>T");
        assert!(ok.is_ok());
    }
}
