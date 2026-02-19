//! Ferro tool service implementation

use std::sync::Arc;

use crate::benchmark::types::{ParseResult, ParsedVariantDetails, PositionDetails};
use crate::error_handling::{ErrorConfig, ErrorMode};
use crate::hgvs::edit::NaEdit;
use crate::hgvs::variant::HgvsVariant;
use crate::normalize::{NormalizeConfig, Normalizer, ShuffleDirection};
use crate::reference::MultiFastaProvider;
use crate::service::{
    config::FerroConfig,
    tools::HgvsToolService,
    types::{health_check::HealthCheckResult, ServiceError, ToolName},
};

/// Extract parsed variant details from an HgvsVariant
fn extract_variant_details(variant: &HgvsVariant) -> Option<ParsedVariantDetails> {
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
        _ => None, // For alleles, RNA fusions, etc. - skip for now
    }
}

/// Extract edit type and sequences from NaEdit
fn extract_na_edit_info(edit: &NaEdit) -> (String, Option<String>, Option<String>) {
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

/// Ferro tool service
pub struct FerroService {
    /// Normalizer instance with reference provider
    normalizer: Arc<Normalizer<MultiFastaProvider>>,
    /// Configuration
    #[allow(dead_code)]
    config: FerroConfig,
}

impl FerroService {
    /// Create a new FerroService
    pub fn new(config: &FerroConfig) -> Result<Self, ServiceError> {
        // Load reference provider - prefer manifest if available (includes cdot CDS info)
        let manifest_path = std::path::Path::new(&config.reference_dir).join("manifest.json");
        let provider = if manifest_path.exists() {
            MultiFastaProvider::from_manifest(&manifest_path).map_err(|e| {
                ServiceError::ConfigError(format!(
                    "Failed to load ferro reference from manifest: {}",
                    e
                ))
            })?
        } else {
            MultiFastaProvider::from_directory(&config.reference_dir).map_err(|e| {
                ServiceError::ConfigError(format!("Failed to load ferro reference data: {}", e))
            })?
        };

        // Create normalization config
        let shuffle_direction = match config.shuffle_direction.as_deref() {
            Some("5prime") => ShuffleDirection::FivePrime,
            Some("3prime") | None => ShuffleDirection::ThreePrime,
            Some(other) => {
                return Err(ServiceError::ConfigError(format!(
                    "Invalid shuffle_direction '{}', must be '3prime' or '5prime'",
                    other
                )));
            }
        };

        let error_mode = match config.error_mode.as_deref() {
            Some("strict") => ErrorMode::Strict,
            Some("lenient") | None => ErrorMode::Lenient,
            Some("silent") => ErrorMode::Silent,
            Some(other) => {
                return Err(ServiceError::ConfigError(format!(
                    "Invalid error_mode '{}', must be 'strict', 'lenient', or 'silent'",
                    other
                )));
            }
        };

        let normalize_config = NormalizeConfig {
            shuffle_direction,
            cross_boundaries: false, // Keep default
            error_config: ErrorConfig::new(error_mode),
            window_size: 100,       // Keep default
            prevent_overlap: false, // Keep default
        };

        // Create normalizer
        let normalizer = Arc::new(Normalizer::with_config(provider, normalize_config));

        Ok(Self {
            normalizer,
            config: config.clone(),
        })
    }

    /// Parse HGVS using ferro parser
    async fn parse_hgvs(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        // Use ferro's parse functionality
        // We'll run this in a blocking task since ferro is sync
        let hgvs = hgvs.to_string();

        tokio::task::spawn_blocking(move || {
            // Parse the HGVS string
            match crate::hgvs::parser::parse_hgvs_lenient(&hgvs) {
                Ok(parse_result) => {
                    let success = true;
                    let output = Some(parse_result.result.to_string());
                    let error = if parse_result.warnings.is_empty() {
                        None
                    } else {
                        Some(format!("Warnings: {}", parse_result.warnings.len()))
                    };

                    // Extract parsed details
                    let details = extract_variant_details(&parse_result.result);

                    ParseResult {
                        input: hgvs,
                        success,
                        output,
                        error,
                        error_category: None,
                        ref_mismatch: None,
                        details,
                    }
                }
                Err(e) => ParseResult {
                    input: hgvs,
                    success: false,
                    output: None,
                    error: Some(e.to_string()),
                    error_category: None,
                    ref_mismatch: None,
                    details: None,
                },
            }
        })
        .await
        .map_err(|e| ServiceError::InternalError(format!("Task join error: {}", e)))
    }

    /// Normalize HGVS using ferro normalizer
    async fn normalize_hgvs(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        let hgvs = hgvs.to_string();
        let normalizer = self.normalizer.clone();

        tokio::task::spawn_blocking(move || {
            // First parse the HGVS string
            match crate::hgvs::parser::parse_hgvs_lenient(&hgvs) {
                Ok(parse_result) => {
                    let original_str = parse_result.result.to_string();

                    // Then normalize it
                    match normalizer.normalize_with_warnings(&parse_result.result) {
                        Ok(normalize_result) => {
                            let success = true;
                            let normalized_str = normalize_result.result.to_string();
                            let output = Some(normalized_str.clone());
                            let error = if normalize_result.warnings.is_empty() {
                                None
                            } else {
                                Some(format!("Warnings: {}", normalize_result.warnings.len()))
                            };

                            // Extract details from normalized variant and detect shifting
                            let mut details = extract_variant_details(&normalize_result.result);
                            if let Some(ref mut d) = details {
                                // Detect if the variant was shifted
                                let was_shifted = normalized_str != original_str;
                                d.was_shifted = Some(was_shifted);
                                if was_shifted {
                                    // Extract original position for display
                                    if let Some(orig_details) =
                                        extract_variant_details(&parse_result.result)
                                    {
                                        d.original_position = Some(orig_details.position.display);
                                    }
                                }
                            }

                            ParseResult {
                                input: hgvs,
                                success,
                                output,
                                error,
                                error_category: None,
                                ref_mismatch: None,
                                details,
                            }
                        }
                        Err(e) => ParseResult {
                            input: hgvs,
                            success: false,
                            output: None,
                            error: Some(e.to_string()),
                            error_category: None,
                            ref_mismatch: None,
                            details: None,
                        },
                    }
                }
                Err(e) => ParseResult {
                    input: hgvs,
                    success: false,
                    output: None,
                    error: Some(e.to_string()),
                    error_category: None,
                    ref_mismatch: None,
                    details: None,
                },
            }
        })
        .await
        .map_err(|e| ServiceError::InternalError(format!("Task join error: {}", e)))
    }
}

#[async_trait::async_trait]
impl HgvsToolService for FerroService {
    async fn parse(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        self.parse_hgvs(hgvs).await
    }

    async fn normalize(&self, hgvs: &str) -> Result<ParseResult, ServiceError> {
        self.normalize_hgvs(hgvs).await
    }

    async fn health_check(&self) -> HealthCheckResult {
        // Simple health check - try to parse a basic HGVS variant
        // This tests that the ferro parser and normalizer are functioning
        let result = tokio::task::spawn_blocking(move || {
            // Try parsing a simple variant to verify ferro is working
            crate::hgvs::parser::parse_hgvs("NM_000001.2:c.1A>G")
        })
        .await;

        match result {
            Ok(Ok(_)) => HealthCheckResult::Healthy,
            Ok(Err(e)) => {
                // Parse errors are expected if the variant is unknown, but parser is working
                let error_msg = e.to_string();
                if error_msg.contains("transcript") || error_msg.contains("not found") {
                    HealthCheckResult::Degraded {
                        reason: "Parser working but reference data may be incomplete".to_string(),
                    }
                } else {
                    HealthCheckResult::Unhealthy {
                        reason: format!("Parser health check failed: {}", e),
                    }
                }
            }
            Err(e) => HealthCheckResult::Unhealthy {
                reason: format!("Health check task error: {}", e),
            },
        }
    }

    fn tool_name(&self) -> ToolName {
        ToolName::Ferro
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    fn create_test_config() -> Result<(FerroConfig, TempDir), Box<dyn std::error::Error>> {
        let temp_dir = TempDir::new()?;

        // Create a minimal reference structure
        fs::create_dir_all(temp_dir.path().join("transcripts"))?;
        fs::write(
            temp_dir.path().join("manifest.json"),
            r#"{"version": "test", "files": []}"#,
        )?;

        let config = FerroConfig {
            enabled: true,
            reference_dir: temp_dir.path().to_path_buf(),
            parallel_workers: Some(1),
            shuffle_direction: Some("3prime".to_string()),
            error_mode: Some("lenient".to_string()),
        };

        Ok((config, temp_dir))
    }

    #[tokio::test]
    async fn test_ferro_service_creation() {
        let (config, _temp_dir) = create_test_config().expect("Failed to create test config");

        // This might fail if reference data is not available, but the service should be creatable
        let result = FerroService::new(&config);

        // We expect this to fail in test environment without proper reference data
        // but the error should be a config error, not a panic
        if let Err(e) = result {
            assert!(matches!(e, ServiceError::ConfigError(_)));
        }
    }

    #[test]
    fn test_invalid_config() {
        let config = FerroConfig {
            enabled: true,
            reference_dir: std::path::PathBuf::from("/nonexistent/path"),
            parallel_workers: Some(1),
            shuffle_direction: Some("invalid".to_string()),
            error_mode: Some("lenient".to_string()),
        };

        let result = FerroService::new(&config);
        assert!(result.is_err());
    }
}
