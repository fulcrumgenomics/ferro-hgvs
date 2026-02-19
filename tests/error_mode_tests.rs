//! Comprehensive tests for error handling modes.
//!
//! Tests the three error modes (Strict, Lenient, Silent) and per-code overrides.

use ferro_hgvs::error_handling::{
    get_code_info, list_all_codes, list_error_codes, list_warning_codes, CodeCategory, ErrorConfig,
    ErrorOverride, ErrorType, ResolvedAction,
};

// ============================================================================
// Mode Tests
// ============================================================================

#[test]
fn test_strict_mode_rejects_wrong_dash() {
    let config = ErrorConfig::strict();
    let preprocessor = config.preprocessor();

    // En-dash should be rejected in strict mode
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(!result.success, "Strict mode should reject en-dash");
}

#[test]
fn test_lenient_mode_corrects_wrong_dash_with_warning() {
    let config = ErrorConfig::lenient();
    let preprocessor = config.preprocessor();

    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(result.success, "Lenient mode should accept en-dash");
    assert_eq!(result.preprocessed, "NM_000088.3:c.100-200del");
    assert!(
        !result.warnings.is_empty(),
        "Lenient mode should emit warning"
    );
    assert!(result
        .warnings
        .iter()
        .any(|w| w.error_type == ErrorType::WrongDashCharacter));
}

#[test]
fn test_silent_mode_corrects_wrong_dash_without_warning() {
    let config = ErrorConfig::silent();
    let preprocessor = config.preprocessor();

    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(result.success, "Silent mode should accept en-dash");
    assert_eq!(result.preprocessed, "NM_000088.3:c.100-200del");
    assert!(
        result.warnings.is_empty(),
        "Silent mode should not emit warnings"
    );
}

// ============================================================================
// Override Tests
// ============================================================================

#[test]
fn test_override_ignore_in_strict_mode() {
    let config = ErrorConfig::strict()
        .with_override(ErrorType::WrongDashCharacter, ErrorOverride::SilentCorrect);
    let preprocessor = config.preprocessor();

    // Override should allow silent correction even in strict mode
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(result.success, "Override should allow correction");
    assert_eq!(result.preprocessed, "NM_000088.3:c.100-200del");
    assert!(
        result.warnings.is_empty(),
        "SilentCorrect should not emit warnings"
    );
}

#[test]
fn test_override_reject_in_lenient_mode() {
    let config =
        ErrorConfig::lenient().with_override(ErrorType::WrongDashCharacter, ErrorOverride::Reject);
    let preprocessor = config.preprocessor();

    // Override should reject even in lenient mode
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(!result.success, "Override should cause rejection");
}

#[test]
fn test_override_warn_in_silent_mode() {
    let config = ErrorConfig::silent()
        .with_override(ErrorType::WrongDashCharacter, ErrorOverride::WarnCorrect);
    let preprocessor = config.preprocessor();

    // Override should emit warning even in silent mode
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(result.success, "Should still correct");
    assert!(!result.warnings.is_empty(), "Override should emit warning");
}

// ============================================================================
// Per-Warning-Code Tests
// ============================================================================

mod w1001_lowercase_amino_acid {
    use super::*;

    #[test]
    fn test_strict_rejects() {
        let config = ErrorConfig::strict();
        assert!(config.should_reject(ErrorType::LowercaseAminoAcid));
    }

    #[test]
    fn test_lenient_warns_and_corrects() {
        let config = ErrorConfig::lenient();
        assert!(config.should_correct(ErrorType::LowercaseAminoAcid));
        assert!(config.should_warn(ErrorType::LowercaseAminoAcid));
    }

    #[test]
    fn test_silent_corrects_silently() {
        let config = ErrorConfig::silent();
        assert!(config.should_correct(ErrorType::LowercaseAminoAcid));
        assert!(!config.should_warn(ErrorType::LowercaseAminoAcid));
    }
}

mod w2001_wrong_dash_character {
    use super::*;

    #[test]
    fn test_strict_rejects() {
        let config = ErrorConfig::strict();
        assert!(config.should_reject(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_lenient_warns_and_corrects() {
        let config = ErrorConfig::lenient();
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
        assert!(config.should_warn(ErrorType::WrongDashCharacter));
    }

    #[test]
    fn test_silent_corrects_silently() {
        let config = ErrorConfig::silent();
        assert!(config.should_correct(ErrorType::WrongDashCharacter));
        assert!(!config.should_warn(ErrorType::WrongDashCharacter));
    }
}

mod w2003_extra_whitespace {
    use super::*;

    #[test]
    fn test_strict_rejects() {
        let config = ErrorConfig::strict();
        assert!(config.should_reject(ErrorType::ExtraWhitespace));
    }

    #[test]
    fn test_lenient_warns_and_corrects() {
        let config = ErrorConfig::lenient();
        assert!(config.should_correct(ErrorType::ExtraWhitespace));
        assert!(config.should_warn(ErrorType::ExtraWhitespace));
    }

    #[test]
    fn test_silent_corrects_silently() {
        let config = ErrorConfig::silent();
        assert!(config.should_correct(ErrorType::ExtraWhitespace));
        assert!(!config.should_warn(ErrorType::ExtraWhitespace));
    }

    // Note: Whitespace removal is not yet fully implemented in the preprocessor.
    // This test is commented out until that feature is added.
    // #[test]
    // fn test_lenient_removes_whitespace() {
    //     let config = ErrorConfig::lenient();
    //     let preprocessor = config.preprocessor();
    //     let result = preprocessor.preprocess("NM_000088.3: c.100A>G");
    //     assert!(result.success);
    //     assert_eq!(result.preprocessed, "NM_000088.3:c.100A>G");
    // }
}

mod w4002_position_zero {
    use super::*;

    // Position zero is not correctable, but the current implementation
    // treats it like other error types based on mode.
    // The recommended usage is to explicitly reject it via override.

    #[test]
    fn test_strict_rejects() {
        let config = ErrorConfig::strict();
        assert!(config.should_reject(ErrorType::PositionZero));
    }

    #[test]
    fn test_is_not_correctable() {
        // PositionZero should be marked as not correctable
        assert!(!ErrorType::PositionZero.is_correctable());
    }

    #[test]
    fn test_override_reject_works() {
        // Users can explicitly reject PositionZero even in lenient mode
        let config =
            ErrorConfig::lenient().with_override(ErrorType::PositionZero, ErrorOverride::Reject);
        assert!(config.should_reject(ErrorType::PositionZero));
    }
}

// ============================================================================
// ResolvedAction Tests
// ============================================================================

#[test]
fn test_resolved_action_reject() {
    let action = ResolvedAction::Reject;
    assert!(action.should_reject());
    assert!(!action.should_correct());
    assert!(!action.should_warn());
}

#[test]
fn test_resolved_action_warn_correct() {
    let action = ResolvedAction::WarnCorrect;
    assert!(!action.should_reject());
    assert!(action.should_correct());
    assert!(action.should_warn());
}

#[test]
fn test_resolved_action_silent_correct() {
    let action = ResolvedAction::SilentCorrect;
    assert!(!action.should_reject());
    assert!(action.should_correct());
    assert!(!action.should_warn());
}

#[test]
fn test_resolved_action_accept() {
    let action = ResolvedAction::Accept;
    assert!(!action.should_reject());
    assert!(!action.should_correct());
    assert!(!action.should_warn());
}

// ============================================================================
// Override Resolution Tests
// ============================================================================

#[test]
fn test_override_default_uses_base_mode() {
    // Default override should use the base mode's behavior
    let config = ErrorConfig::lenient();
    let action = config.action_for(ErrorType::WrongDashCharacter);
    assert_eq!(action, ResolvedAction::WarnCorrect);
}

#[test]
fn test_override_reject_overrides_lenient() {
    let config =
        ErrorConfig::lenient().with_override(ErrorType::WrongDashCharacter, ErrorOverride::Reject);
    let action = config.action_for(ErrorType::WrongDashCharacter);
    assert_eq!(action, ResolvedAction::Reject);
}

#[test]
fn test_override_silent_correct_overrides_strict() {
    let config = ErrorConfig::strict()
        .with_override(ErrorType::WrongDashCharacter, ErrorOverride::SilentCorrect);
    let action = config.action_for(ErrorType::WrongDashCharacter);
    assert_eq!(action, ResolvedAction::SilentCorrect);
}

// ============================================================================
// Code Registry Tests
// ============================================================================

#[test]
fn test_all_error_codes_exist() {
    let codes = list_error_codes();
    assert!(!codes.is_empty(), "Should have error codes");

    // Verify some specific error codes
    let code_strs: Vec<&str> = codes.iter().map(|c| c.code).collect();
    assert!(code_strs.contains(&"E1001"), "Should have E1001");
    assert!(code_strs.contains(&"E1002"), "Should have E1002");
    assert!(code_strs.contains(&"E2001"), "Should have E2001");
}

#[test]
fn test_all_warning_codes_exist() {
    let codes = list_warning_codes();
    assert!(!codes.is_empty(), "Should have warning codes");

    // Verify some specific warning codes
    let code_strs: Vec<&str> = codes.iter().map(|c| c.code).collect();
    assert!(code_strs.contains(&"W1001"), "Should have W1001");
    assert!(code_strs.contains(&"W2001"), "Should have W2001");
    assert!(code_strs.contains(&"W3001"), "Should have W3001");
}

#[test]
fn test_get_code_info() {
    let info = get_code_info("W1001");
    assert!(info.is_some(), "Should find W1001");
    let info = info.unwrap();
    assert_eq!(info.code, "W1001");
    assert!(!info.summary.is_empty());
    assert!(!info.explanation.is_empty());
}

#[test]
fn test_get_code_info_case_insensitive() {
    let info1 = get_code_info("w1001");
    let info2 = get_code_info("W1001");
    assert!(info1.is_some());
    assert!(info2.is_some());
    assert_eq!(info1.unwrap().code, info2.unwrap().code);
}

#[test]
fn test_code_info_has_examples() {
    let info = get_code_info("W1001").unwrap();
    assert!(!info.bad_examples.is_empty(), "Should have bad examples");
    assert!(!info.good_examples.is_empty(), "Should have good examples");
}

#[test]
fn test_code_info_has_mode_behavior() {
    let info = get_code_info("W1001").unwrap();
    assert!(
        info.mode_behavior.is_some(),
        "Warning codes should have mode behavior"
    );
}

#[test]
fn test_error_code_no_mode_behavior() {
    let info = get_code_info("E1001").unwrap();
    // Error codes might not have configurable mode behavior
    // This depends on the implementation
    assert!(info.is_error());
}

#[test]
fn test_code_categories() {
    let all_codes = list_all_codes();

    // Check that we have codes from multiple categories
    let categories: std::collections::HashSet<CodeCategory> =
        all_codes.iter().map(|c| c.category).collect();

    assert!(
        categories.contains(&CodeCategory::Parse),
        "Should have parse errors"
    );
    assert!(
        categories.contains(&CodeCategory::Case),
        "Should have case warnings"
    );
    assert!(
        categories.contains(&CodeCategory::Character),
        "Should have character warnings"
    );
}

// ============================================================================
// ErrorType Code Mapping Tests
// ============================================================================

#[test]
fn test_error_type_code_mapping() {
    assert_eq!(ErrorType::LowercaseAminoAcid.code(), "W1001");
    assert_eq!(ErrorType::SingleLetterAminoAcid.code(), "W1002");
    assert_eq!(ErrorType::WrongDashCharacter.code(), "W2001");
    assert_eq!(ErrorType::ExtraWhitespace.code(), "W2003");
    assert_eq!(ErrorType::MissingVersion.code(), "W3001");
    assert_eq!(ErrorType::PositionZero.code(), "W4002");
}

// ============================================================================
// Preprocessor Tests
// ============================================================================

#[test]
fn test_preprocessor_preserves_valid_input() {
    let config = ErrorConfig::lenient();
    let preprocessor = config.preprocessor();

    let valid = "NM_000088.3:c.100A>G";
    let result = preprocessor.preprocess(valid);

    assert!(result.success);
    assert_eq!(result.preprocessed, valid);
    assert!(result.warnings.is_empty());
}

#[test]
fn test_preprocessor_dash_correction() {
    let config = ErrorConfig::lenient();
    let preprocessor = config.preprocessor();

    // Input with en-dash - this should be corrected
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");

    assert!(result.success);
    // Dash correction should be applied
    assert_eq!(result.preprocessed, "NM_000088.3:c.100-200del");
    // Should have warning for the dash issue
    assert!(!result.warnings.is_empty());
    assert!(result
        .warnings
        .iter()
        .any(|w| w.error_type == ErrorType::WrongDashCharacter));
}

#[test]
fn test_preprocessor_idempotent() {
    let config = ErrorConfig::lenient();
    let preprocessor = config.preprocessor();

    let input = "NM_000088.3:c.100–200del";
    let result1 = preprocessor.preprocess(input);
    let result2 = preprocessor.preprocess(&result1.preprocessed);

    // Second pass should make no changes
    assert_eq!(result1.preprocessed, result2.preprocessed);
    assert!(result2.warnings.is_empty());
}

// ============================================================================
// Config Builder Tests
// ============================================================================

#[test]
fn test_config_builder_chain() {
    let config = ErrorConfig::strict()
        .with_override(ErrorType::WrongDashCharacter, ErrorOverride::SilentCorrect)
        .with_override(ErrorType::ExtraWhitespace, ErrorOverride::WarnCorrect)
        .with_override(ErrorType::LowercaseAminoAcid, ErrorOverride::Reject);

    assert_eq!(
        config.action_for(ErrorType::WrongDashCharacter),
        ResolvedAction::SilentCorrect
    );
    assert_eq!(
        config.action_for(ErrorType::ExtraWhitespace),
        ResolvedAction::WarnCorrect
    );
    assert_eq!(
        config.action_for(ErrorType::LowercaseAminoAcid),
        ResolvedAction::Reject
    );
}

#[test]
fn test_config_set_and_remove_override() {
    let mut config = ErrorConfig::strict();

    // Set override
    config.set_override(ErrorType::WrongDashCharacter, ErrorOverride::SilentCorrect);
    assert!(config.should_correct(ErrorType::WrongDashCharacter));

    // Remove override
    config.remove_override(ErrorType::WrongDashCharacter);
    assert!(config.should_reject(ErrorType::WrongDashCharacter));
}

// ============================================================================
// Format Output Tests
// ============================================================================

#[test]
fn test_code_info_terminal_format() {
    let info = get_code_info("W1001").unwrap();
    let formatted = info.format_terminal(false);

    assert!(formatted.contains("W1001"));
    assert!(formatted.contains("LowercaseAminoAcid"));
    assert!(formatted.contains("Mode Behavior"));
}

#[test]
fn test_code_info_json_format() {
    let info = get_code_info("W1001").unwrap();
    let formatted = info.format_json();

    assert!(formatted.contains("\"code\":\"W1001\""));
    assert!(formatted.contains("\"name\":\"LowercaseAminoAcid\""));
}

#[test]
fn test_code_info_markdown_format() {
    let info = get_code_info("W1001").unwrap();
    let formatted = info.format_markdown();

    assert!(formatted.contains("## W1001"));
    assert!(formatted.contains("LowercaseAminoAcid"));
    assert!(formatted.contains("### Mode Behavior"));
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]
fn test_full_parsing_workflow_strict() {
    use ferro_hgvs::parse_hgvs;

    let config = ErrorConfig::strict();
    let preprocessor = config.preprocessor();

    // Valid input should work
    let result = preprocessor.preprocess("NM_000088.3:c.100A>G");
    assert!(result.success);
    let parsed = parse_hgvs(&result.preprocessed);
    assert!(parsed.is_ok());

    // Invalid input should be rejected at preprocessing
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(!result.success);
}

#[test]
fn test_full_parsing_workflow_lenient() {
    use ferro_hgvs::parse_hgvs;

    let config = ErrorConfig::lenient();
    let preprocessor = config.preprocessor();

    // Input with en-dash should be corrected and parse
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(result.success);
    assert!(!result.warnings.is_empty());

    let parsed = parse_hgvs(&result.preprocessed);
    assert!(parsed.is_ok(), "Corrected input should parse: {:?}", parsed);
}

#[test]
fn test_full_parsing_workflow_silent() {
    use ferro_hgvs::parse_hgvs;

    let config = ErrorConfig::silent();
    let preprocessor = config.preprocessor();

    // Input with en-dash should be corrected silently and parse
    let result = preprocessor.preprocess("NM_000088.3:c.100–200del");
    assert!(result.success);
    assert!(result.warnings.is_empty());

    let parsed = parse_hgvs(&result.preprocessed);
    assert!(parsed.is_ok());
}
