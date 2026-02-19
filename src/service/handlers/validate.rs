//! HGVS validation endpoint - validates syntax and shows component breakdown

use axum::{extract::State, http::StatusCode, response::Json};
use std::time::Instant;

use crate::service::{
    server::AppState,
    types::{ErrorResponse, ParsedVariantDetails, PositionDetails, ServiceError, ValidateResponse},
    validation::validate_hgvs,
};

/// Request for single variant validation
#[derive(Debug, serde::Deserialize)]
pub struct ValidateRequest {
    /// The HGVS variant string to validate
    pub hgvs: String,
}

/// Validate a single HGVS variant and return component breakdown
///
/// This endpoint parses the HGVS string and returns:
/// - Whether the syntax is valid
/// - Any validation errors or warnings
/// - A breakdown of the variant components (reference, position, edit type, etc.)
pub async fn validate_single(
    State(_state): State<AppState>,
    Json(request): Json<ValidateRequest>,
) -> Result<Json<ValidateResponse>, (StatusCode, Json<ErrorResponse>)> {
    let start = Instant::now();

    // First, do security validation on the input
    if let Err(validation_error) = validate_hgvs(&request.hgvs) {
        let error = ServiceError::InvalidHgvs(validation_error.to_string());
        return Err((
            StatusCode::from_u16(error.status_code()).unwrap_or(StatusCode::INTERNAL_SERVER_ERROR),
            Json(error.to_response()),
        ));
    }

    // Parse using ferro's lenient parser
    let hgvs_str = request.hgvs.clone();
    let parse_result =
        tokio::task::spawn_blocking(move || crate::hgvs::parser::parse_hgvs_lenient(&hgvs_str))
            .await
            .map_err(|e| {
                let error = ServiceError::InternalError(format!("Task error: {}", e));
                (StatusCode::INTERNAL_SERVER_ERROR, Json(error.to_response()))
            })?;

    let elapsed_ms = start.elapsed().as_millis() as u64;

    match parse_result {
        Ok(result) => {
            // Extract component details
            let components = extract_variant_details(&result.result);

            // Collect any warnings as potential issues
            let errors = if result.warnings.is_empty() {
                None
            } else {
                Some(result.warnings.iter().map(|w| w.message.clone()).collect())
            };

            Ok(Json(ValidateResponse {
                input: request.hgvs,
                valid: true,
                errors,
                components,
                processing_time_ms: elapsed_ms,
            }))
        }
        Err(e) => Ok(Json(ValidateResponse {
            input: request.hgvs,
            valid: false,
            errors: Some(vec![e.to_string()]),
            components: None,
            processing_time_ms: elapsed_ms,
        })),
    }
}

/// Extract parsed variant details from an HgvsVariant
fn extract_variant_details(
    variant: &crate::hgvs::variant::HgvsVariant,
) -> Option<ParsedVariantDetails> {
    use crate::hgvs::variant::HgvsVariant;

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
                    start: 0,
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
                    start: 0,
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
                    start: 0,
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
                    start: 0,
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
        _ => None,
    }
}

/// Extract edit type and sequences from NaEdit
fn extract_na_edit_info(
    edit: &crate::hgvs::edit::NaEdit,
) -> (String, Option<String>, Option<String>) {
    use crate::hgvs::edit::NaEdit;

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
        NaEdit::Delins { sequence } => ("delins".to_string(), None, Some(sequence.to_string())),
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
        _ => ("other".to_string(), None, None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_variant_details_cds() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350C>T").unwrap();
        let details = extract_variant_details(&result.result);

        assert!(details.is_some());
        let d = details.unwrap();
        assert_eq!(d.coordinate_system, "c");
        assert_eq!(d.variant_type, "substitution");
        assert_eq!(d.reference, "NM_000249.4");
    }

    #[test]
    fn test_extract_variant_details_genomic() {
        let result =
            crate::hgvs::parser::parse_hgvs_lenient("NC_000007.14:g.117559593G>A").unwrap();
        let details = extract_variant_details(&result.result);

        assert!(details.is_some());
        let d = details.unwrap();
        assert_eq!(d.coordinate_system, "g");
        assert_eq!(d.variant_type, "substitution");
    }

    #[test]
    fn test_extract_variant_details_protein() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NP_000240.1:p.Val600Glu").unwrap();
        let details = extract_variant_details(&result.result);

        assert!(details.is_some());
        let d = details.unwrap();
        assert_eq!(d.coordinate_system, "p");
        assert_eq!(d.variant_type, "protein_change");
    }

    #[test]
    fn test_extract_variant_details_mitochondrial() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NC_012920.1:m.8993T>G").unwrap();
        let details = extract_variant_details(&result.result);

        assert!(details.is_some());
        let d = details.unwrap();
        assert_eq!(d.coordinate_system, "m");
    }

    #[test]
    fn test_extract_variant_details_noncoding() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NR_000001.1:n.100A>G").unwrap();
        let details = extract_variant_details(&result.result);

        assert!(details.is_some());
        let d = details.unwrap();
        assert_eq!(d.coordinate_system, "n");
    }

    #[test]
    fn test_extract_na_edit_info_substitution() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350C>T").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (vtype, deleted, inserted) = extract_na_edit_info(edit);
                assert_eq!(vtype, "substitution");
                assert_eq!(deleted, Some("C".to_string()));
                assert_eq!(inserted, Some("T".to_string()));
            }
        }
    }

    #[test]
    fn test_extract_na_edit_info_deletion() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350delC").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (vtype, deleted, inserted) = extract_na_edit_info(edit);
                assert_eq!(vtype, "deletion");
                assert_eq!(deleted, Some("C".to_string()));
                assert!(inserted.is_none());
            }
        }
    }

    #[test]
    fn test_extract_na_edit_info_insertion() {
        let result =
            crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350_351insATG").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (vtype, deleted, inserted) = extract_na_edit_info(edit);
                assert_eq!(vtype, "insertion");
                assert!(deleted.is_none());
                assert_eq!(inserted, Some("ATG".to_string()));
            }
        }
    }

    #[test]
    fn test_extract_na_edit_info_delins() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350delinsATG").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (vtype, _deleted, inserted) = extract_na_edit_info(edit);
                assert_eq!(vtype, "delins");
                assert_eq!(inserted, Some("ATG".to_string()));
            }
        }
    }

    #[test]
    fn test_extract_na_edit_info_duplication() {
        let result = crate::hgvs::parser::parse_hgvs_lenient("NM_000249.4:c.350dupC").unwrap();
        if let crate::hgvs::variant::HgvsVariant::Cds(v) = &result.result {
            if let Some(edit) = v.loc_edit.edit.inner() {
                let (vtype, deleted, _inserted) = extract_na_edit_info(edit);
                assert_eq!(vtype, "duplication");
                assert_eq!(deleted, Some("C".to_string()));
            }
        }
    }
}
