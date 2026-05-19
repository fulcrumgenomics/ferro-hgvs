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

    #[test]
    fn test_lenient_removes_whitespace_around_colon() {
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NM_000088.3: c.100A>G");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.100A>G");
        assert!(result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::ExtraWhitespace));
    }

    // Mid-expression whitespace (between digits and a base letter, between
    // underscore and the second coordinate, around `;` in an allele list, etc.)
    // is non-canonical per the HGVS spec but lenient mode should accept-and-warn
    // rather than hard-reject. See #128.

    #[test]
    fn test_lenient_strips_whitespace_between_position_and_base() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c.100 A>G");
        assert!(
            result.is_ok(),
            "lenient parse should accept embedded whitespace: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert!(parsed
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::ExtraWhitespace));
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.100A>G");
    }

    #[test]
    fn test_lenient_strips_whitespace_between_underscore_and_position() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c.100_ 200del");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace after `_`: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.100_200del");
    }

    #[test]
    fn test_lenient_strips_whitespace_around_allele_separator() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c.[100A>G ; 200T>C]");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace around `;`: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.[100A>G;200T>C]");
    }

    #[test]
    fn test_lenient_strips_whitespace_inside_brackets() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c.[ 100A>G;200T>C ]");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace inside `[ ]`: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.[100A>G;200T>C]");
    }

    #[test]
    fn test_silent_strips_embedded_whitespace_without_warning() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_silent;
        let result = parse_hgvs_silent("NM_000088.3:c.100 A>G");
        assert!(result.is_ok());
        let parsed = result.unwrap();
        assert!(!parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.100A>G");
    }

    #[test]
    fn test_strict_still_rejects_embedded_whitespace() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
        let config = ErrorConfig::strict();
        let result = parse_hgvs_with_config("NM_000088.3:c.100 A>G", config);
        assert!(
            result.is_err(),
            "strict mode must continue to reject embedded whitespace"
        );
    }

    // Coverage for additional HGVS contexts brainstormed alongside #128:
    // protein, accession-internal, repeat, uncertain-position, compact allele,
    // zero-width invisible characters, multi-run warning counts, idempotency.

    #[test]
    fn test_lenient_strips_whitespace_in_protein() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NP_000079.2:p.Arg 97 Trp");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace in protein: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NP_000079.2:p.Arg97Trp");
    }

    #[test]
    fn test_lenient_strips_whitespace_around_accession_dot() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088 .3:c.459A>G");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace around accession dot: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.459A>G");
    }

    #[test]
    fn test_lenient_strips_whitespace_in_repeat_block() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c.123_125 CAG[10]");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace in repeat: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.123_125CAG[10]");
    }

    #[test]
    fn test_lenient_strips_whitespace_inside_uncertain_position() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c.( 123 _ 127 )delA");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace inside uncertain position: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.(123_127)delA");
    }

    #[test]
    fn test_lenient_strips_whitespace_in_compact_allele() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        let result = parse_hgvs_lenient("NM_000088.3:c. [ 100A>G ; 200T>C ]");
        assert!(
            result.is_ok(),
            "lenient parse should accept whitespace in compact allele: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.[100A>G;200T>C]");
    }

    #[test]
    fn test_lenient_strips_zero_width_space() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        // Zero-width space pasted between accession and colon.
        let result = parse_hgvs_lenient("NM_000088.3\u{200B}:c.100A>G");
        assert!(
            result.is_ok(),
            "lenient parse should accept zero-width space: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.100A>G");
    }

    #[test]
    fn test_lenient_strips_byte_order_mark() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        // BOM/ZWNBSP at start of input — common when copy-pasting from a
        // UTF-8-with-BOM file.
        let result = parse_hgvs_lenient("\u{FEFF}NM_000088.3:c.100A>G");
        assert!(
            result.is_ok(),
            "lenient parse should accept leading BOM: {:?}",
            result
        );
        let parsed = result.unwrap();
        assert!(parsed.has_warnings());
        assert_eq!(parsed.preprocessed_input, "NM_000088.3:c.100A>G");
    }

    #[test]
    fn test_lenient_one_warning_per_run_not_per_char() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        // Three runs of whitespace -> exactly three W2003 warnings.
        let result = parse_hgvs_lenient(" NM_000088.3 :c.100 A>G").unwrap();
        let ws_warnings: Vec<_> = result
            .warnings
            .iter()
            .filter(|w| w.error_type == ErrorType::ExtraWhitespace)
            .collect();
        assert_eq!(
            ws_warnings.len(),
            3,
            "expected one W2003 per run, got: {:?}",
            result.warnings
        );
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.100A>G");
    }

    #[test]
    fn test_lenient_idempotent_on_canonical_form() {
        use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
        // Canonical input -> no whitespace warnings emitted.
        let result = parse_hgvs_lenient("NM_000088.3:c.100A>G").unwrap();
        let ws_warnings: Vec<_> = result
            .warnings
            .iter()
            .filter(|w| w.error_type == ErrorType::ExtraWhitespace)
            .collect();
        assert!(
            ws_warnings.is_empty(),
            "canonical input must not generate W2003"
        );
    }
}

// The `w4002_position_zero` mod (always_reject / not_correctable / override
// behaviour for `ErrorType::PositionZero`) was retired in issue #269. The
// `c.0…` input is now rejected at parse time as `E1003 InvalidPosition`
// — see `test_position_zero_input_still_rejected_with_e1003` below.

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
}

// ============================================================================
// Issue #269: W4002 PositionZero deprecation
// ============================================================================
//
// W4002 was registered but never emitted; the parser's `detect_position_zero`
// already raises `E1003 InvalidPosition` for `c.0…` inputs (a hard error, not
// a soft validation warning). Per the issue's principled long-term resolution
// (Option 3), the W4002 identity is removed from the registry, ErrorType
// enum, Python bindings, and config parser so there is one canonical code
// for the symptom.

#[test]
fn test_position_zero_input_still_rejected_with_e1003() {
    // Behavior preserved across the deprecation: `c.0A>G` is rejected with
    // E1003, the canonical code documented in the audit.
    use ferro_hgvs::hgvs::parser::parse_hgvs_with_config;
    let cfg = ErrorConfig::strict();
    let err = parse_hgvs_with_config("NM_000088.3:c.0A>G", cfg).expect_err("c.0 must be rejected");
    assert_eq!(
        err.code(),
        Some(ferro_hgvs::error::ErrorCode::InvalidPosition),
        "c.0 rejection must carry E1003 InvalidPosition (canonical)"
    );
}

#[test]
fn test_w4002_no_longer_in_registry() {
    assert!(
        get_code_info("W4002").is_none(),
        "W4002 must be removed from the registry per issue #269 deprecation"
    );
    assert!(
        !list_warning_codes().iter().any(|c| c.code == "W4002"),
        "W4002 must not appear in list_warning_codes()"
    );
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

// ============================================================================
// W1001 LowercaseAminoAcid emission tests (#124)
// ============================================================================

mod w1001_lowercase_amino_acid_emission {
    use super::*;

    #[test]
    fn test_lenient_emits_w1001_on_lowercase_aa() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.val600glu");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Val600Glu");
        let codes: Vec<&'static str> = result
            .warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect();
        assert_eq!(codes, vec!["W1001", "W1001"]);
    }

    #[test]
    fn test_strict_rejects_lowercase_aa() {
        let preprocessor = ErrorConfig::strict().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.val600glu");
        assert!(!result.success);
    }

    #[test]
    fn test_silent_corrects_lowercase_aa_no_warning() {
        let preprocessor = ErrorConfig::silent().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.val600glu");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Val600Glu");
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn test_canonical_input_no_w1001() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Val600Glu");
        assert!(result.success);
        assert!(result
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::LowercaseAminoAcid));
    }

    #[test]
    fn test_w1001_idempotent() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let first = preprocessor.preprocess("NP_000079.2:p.val600glu");
        assert!(first.success);
        let second = preprocessor.preprocess(&first.preprocessed);
        assert!(second.success);
        assert!(second
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::LowercaseAminoAcid));
    }
}

// ============================================================================
// W1002 SingleLetterAminoAcid emission tests (#124)
// ============================================================================

mod w1002_single_letter_amino_acid_emission {
    use super::*;

    #[test]
    fn test_lenient_emits_w1002_on_single_letter_aa() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.V600E");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Val600Glu");
        let codes: Vec<&'static str> = result
            .warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect();
        assert_eq!(codes, vec!["W1002", "W1002"]);
    }

    #[test]
    fn test_strict_rejects_single_letter_aa() {
        let preprocessor = ErrorConfig::strict().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.V600E");
        assert!(!result.success);
    }

    #[test]
    fn test_silent_corrects_single_letter_aa_no_warning() {
        let preprocessor = ErrorConfig::silent().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.V600E");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Val600Glu");
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn test_canonical_input_no_w1002() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Val600Glu");
        assert!(result.success);
        assert!(result
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::SingleLetterAminoAcid));
    }

    #[test]
    fn test_w1002_idempotent() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let first = preprocessor.preprocess("NP_000079.2:p.V600E");
        assert!(first.success);
        let second = preprocessor.preprocess(&first.preprocessed);
        assert!(second.success);
        assert!(second
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::SingleLetterAminoAcid));
    }

    #[test]
    fn test_w1001_w1002_combined() {
        // Mixed: lowercase three-letter (`glu`) and a single-letter (`V`).
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.V600glu");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Val600Glu");
        let codes: Vec<&'static str> = result
            .warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect();
        assert!(codes.contains(&"W1001"));
        assert!(codes.contains(&"W1002"));
    }
}

// ============================================================================
// W3001 MissingVersion emission tests (#124)
// ============================================================================

mod w3001_missing_version_emission {
    use super::*;

    #[test]
    fn test_lenient_emits_w3001_on_unversioned_accession() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NM_000088:c.100A>G");
        assert!(
            result.success,
            "lenient should accept unversioned accession"
        );
        // Input is NOT mutated: there is no version to inject.
        assert_eq!(result.preprocessed, "NM_000088:c.100A>G");
        let codes: Vec<&'static str> = result
            .warnings
            .iter()
            .map(|w| w.error_type.code())
            .collect();
        assert_eq!(codes, vec!["W3001"]);
    }

    #[test]
    fn test_strict_rejects_unversioned_accession() {
        let preprocessor = ErrorConfig::strict().preprocessor();
        let result = preprocessor.preprocess("NM_000088:c.100A>G");
        assert!(!result.success);
    }

    #[test]
    fn test_silent_accepts_unversioned_no_warning() {
        let preprocessor = ErrorConfig::silent().preprocessor();
        let result = preprocessor.preprocess("NM_000088:c.100A>G");
        assert!(result.success);
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn test_canonical_versioned_no_w3001() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NM_000088.3:c.100A>G");
        assert!(result.success);
        assert!(result
            .warnings
            .iter()
            .all(|w| w.error_type != ErrorType::MissingVersion));
    }

    #[test]
    fn test_w3001_idempotent_on_unversioned() {
        // Re-preprocessing an unversioned accession STILL warns — the input
        // is not auto-corrected, so the warning is persistent. But it warns
        // at most once per occurrence on each pass.
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let first = preprocessor.preprocess("NM_000088:c.100A>G");
        let second = preprocessor.preprocess(&first.preprocessed);
        assert_eq!(
            first
                .warnings
                .iter()
                .filter(|w| w.error_type == ErrorType::MissingVersion)
                .count(),
            1
        );
        assert_eq!(
            second
                .warnings
                .iter()
                .filter(|w| w.error_type == ErrorType::MissingVersion)
                .count(),
            1
        );
    }

    #[test]
    fn test_w3001_inner_outer_accession() {
        let preprocessor = ErrorConfig::lenient().preprocessor();
        let result = preprocessor.preprocess("NG_012232(NM_004006):c.93+1G>T");
        assert!(result.success);
        let count = result
            .warnings
            .iter()
            .filter(|w| w.error_type == ErrorType::MissingVersion)
            .count();
        assert_eq!(count, 2);
    }
}

// ============================================================================
// Deprecated stop-codon and frameshift forms (issue #125, #81 L2 SVA-003..006)
//
// HGVS spec rationale:
//   - recommendations/checklist.md: "'Ter' or '*' should be used to indicate a
//     translation stop codon; the X should not be used"
//   - recommendations/protein/substitution.md: 'X' is reserved for the "any
//     amino acid" symbol Xaa, so 'X' as stop is ambiguous.
//   - recommendations/protein/frameshift.md: canonical form is 'fsTerN'
//     (e.g. p.Arg123LysfsTer34); 'fs*N' is permitted but 'fsTerN' preferred.
// ============================================================================

mod deprecated_stop_codon_and_frameshift_forms {
    use super::*;
    use ferro_hgvs::parse_hgvs;

    /// Helper: assert that a deprecated input parses to the expected canonical
    /// HGVS string in lenient mode, with exactly one warning of the given type.
    fn assert_corrects_with_warning(input: &str, expected: &str, expected_error_type: ErrorType) {
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();

        let result = preprocessor.preprocess(input);
        assert!(
            result.success,
            "lenient should accept {}: {:?}",
            input, result.error
        );
        assert_eq!(
            result.preprocessed, expected,
            "preprocessed mismatch for {}",
            input
        );
        assert_eq!(result.warnings.len(), 1, "expected 1 warning for {}", input);
        assert_eq!(result.warnings[0].error_type, expected_error_type);

        // The corrected form must round-trip through the full parser.
        let parsed = parse_hgvs(&result.preprocessed)
            .unwrap_or_else(|e| panic!("corrected {} should parse: {:?}", expected, e));
        assert_eq!(parsed.to_string(), expected);
    }

    // --- W3007: deprecated `*` for stop in protein substitution (SVA-003) ---

    #[test]
    fn lenient_corrects_w3007_star_for_stop() {
        assert_corrects_with_warning(
            "NP_000079.2:p.Arg97*",
            "NP_000079.2:p.Arg97Ter",
            ErrorType::DeprecatedStopCodonStar,
        );
    }

    #[test]
    fn strict_rejects_w3007_star_for_stop() {
        let config = ErrorConfig::strict();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97*");
        assert!(!result.success);
    }

    #[test]
    fn silent_corrects_w3007_no_warning() {
        let config = ErrorConfig::silent();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97*");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97Ter");
        assert!(result.warnings.is_empty());
    }

    // --- W3008: deprecated `X` for stop in protein substitution (SVA-004) ---

    #[test]
    fn lenient_corrects_w3008_x_for_stop() {
        assert_corrects_with_warning(
            "NP_000079.2:p.Arg97X",
            "NP_000079.2:p.Arg97Ter",
            ErrorType::DeprecatedStopCodonX,
        );
    }

    #[test]
    fn strict_rejects_w3008_x_for_stop() {
        let config = ErrorConfig::strict();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97X");
        assert!(!result.success);
    }

    #[test]
    fn silent_corrects_w3008_no_warning() {
        let config = ErrorConfig::silent();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97X");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg97Ter");
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn lenient_does_not_flag_xaa_any_amino_acid() {
        // 'Xaa' (the 'any amino acid' symbol) is canonical and must not warn.
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg782Xaa");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.Arg782Xaa");
        assert!(result.warnings.is_empty());
    }

    // --- W3009: deprecated `fs*N` frameshift termination (SVA-006) ---

    #[test]
    fn lenient_corrects_w3009_fs_star_n() {
        assert_corrects_with_warning(
            "NP_000079.2:p.Arg97fs*23",
            "NP_000079.2:p.Arg97fsTer23",
            ErrorType::DeprecatedFrameshiftStar,
        );
    }

    #[test]
    fn lenient_corrects_w3009_with_new_aa() {
        assert_corrects_with_warning(
            "NP_000079.2:p.Arg97Profs*23",
            "NP_000079.2:p.Arg97ProfsTer23",
            ErrorType::DeprecatedFrameshiftStar,
        );
    }

    #[test]
    fn strict_rejects_w3009_fs_star_n() {
        let config = ErrorConfig::strict();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fs*23");
        assert!(!result.success);
    }

    // --- W3010: deprecated `fsXN` frameshift termination (SVA-005) ---

    #[test]
    fn lenient_corrects_w3010_fs_x_n() {
        // Without this PR, ferro hard-rejects this with a generic
        // "unexpected trailing characters" parse error. With the soft-warn
        // wired in, lenient mode rewrites to canonical fsTer23 + parses.
        assert_corrects_with_warning(
            "NP_000079.2:p.Arg97fsX23",
            "NP_000079.2:p.Arg97fsTer23",
            ErrorType::DeprecatedFrameshiftX,
        );
    }

    #[test]
    fn strict_rejects_w3010_fs_x_n() {
        let config = ErrorConfig::strict();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.Arg97fsX23");
        assert!(!result.success);
    }

    // --- Edge cases ---

    #[test]
    fn lenient_canonical_protein_no_warnings() {
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();
        for input in [
            "NP_000079.2:p.Arg97Ter",
            "NP_000079.2:p.Arg97ProfsTer23",
            "NP_000079.2:p.Tyr180fs",
            "NP_000079.2:p.Val600Glu",
            "NP_000079.2:p.Arg782Xaa",
        ] {
            let result = preprocessor.preprocess(input);
            assert!(result.success, "expected success for {}", input);
            assert_eq!(result.preprocessed, input);
            assert!(
                result.warnings.is_empty(),
                "expected no warnings for canonical {}",
                input
            );
        }
    }

    #[test]
    fn lenient_compound_protein_allele_emits_two_warnings() {
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NP_000079.2:p.[Arg97*;Arg100X]");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NP_000079.2:p.[Arg97Ter;Arg100Ter]");
        assert_eq!(result.warnings.len(), 2);
    }

    #[test]
    fn lenient_idempotent_on_corrected_output() {
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();
        let first = preprocessor.preprocess("NP_000079.2:p.Arg97fsX23");
        assert!(first.success);
        let second = preprocessor.preprocess(&first.preprocessed);
        assert!(second.success);
        assert_eq!(second.preprocessed, first.preprocessed);
        assert!(
            second.warnings.is_empty(),
            "second-pass warnings: {:?}",
            second.warnings
        );
    }

    #[test]
    fn lenient_does_not_affect_cds_utr_position() {
        // c.*5A>G is a 3'UTR position offset, not a deprecated stop codon.
        let config = ErrorConfig::lenient();
        let preprocessor = config.preprocessor();
        let result = preprocessor.preprocess("NM_000088.3:c.*5A>G");
        assert!(result.success);
        assert_eq!(result.preprocessed, "NM_000088.3:c.*5A>G");
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn registry_has_all_four_codes() {
        for code in ["W3007", "W3008", "W3009", "W3010"] {
            let info = get_code_info(code).unwrap_or_else(|| panic!("missing {}", code));
            assert_eq!(info.code, code);
            assert!(info.category == CodeCategory::Format);
            assert!(info.mode_behavior.is_some());
        }
    }
}

// W3011 / W3012 / W3013 / W4003 — Issue #127 non-canonical input forms
// ============================================================================

mod w3011_del_size_suffix_emission {
    use super::*;
    use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_silent, parse_hgvs_with_config};

    #[test]
    fn lenient_warns_without_rewrite() {
        // SVA-007: warn_accept — lenient warns but does NOT rewrite the input
        // (we cannot synthesise the end position safely).
        let result = parse_hgvs_lenient("NG_012232.1:g.123del6").unwrap();
        assert!(result.has_warnings());
        assert_eq!(
            result
                .warnings
                .iter()
                .filter(|w| w.error_type == ErrorType::DelSizeSuffix)
                .count(),
            1
        );
        // Input is unchanged (no auto-correction).
        assert_eq!(result.preprocessed_input, "NG_012232.1:g.123del6");
    }

    #[test]
    fn strict_rejects() {
        let result = parse_hgvs_with_config("NG_012232.1:g.123del6", ErrorConfig::strict());
        assert!(result.is_err());
    }

    #[test]
    fn silent_accepts_without_warning() {
        let result = parse_hgvs_silent("NG_012232.1:g.123del6").unwrap();
        assert!(!result.has_warnings());
        assert_eq!(result.preprocessed_input, "NG_012232.1:g.123del6");
    }

    #[test]
    fn canonical_input_does_not_warn() {
        let result = parse_hgvs_lenient("NG_012232.1:g.123_128del").unwrap();
        assert!(!result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::DelSizeSuffix));
    }

    #[test]
    fn idempotent_re_pass() {
        let r1 = parse_hgvs_lenient("NG_012232.1:g.123del6").unwrap();
        let r2 = parse_hgvs_lenient(&r1.preprocessed_input).unwrap();
        assert_eq!(r1.preprocessed_input, r2.preprocessed_input);
        assert!(r2.has_warnings());
    }
}

mod w4003_single_position_range_emission {
    use super::*;
    use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_silent, parse_hgvs_with_config};

    #[test]
    fn lenient_warns_and_collapses_del() {
        let result = parse_hgvs_lenient("NM_000088.3:c.123_123del").unwrap();
        assert!(result.has_warnings());
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.123del");
        assert_eq!(
            result
                .warnings
                .iter()
                .filter(|w| w.error_type == ErrorType::SinglePositionRange)
                .count(),
            1
        );
    }

    #[test]
    fn lenient_warns_and_collapses_dup() {
        let result = parse_hgvs_lenient("NM_000088.3:c.123_123dup").unwrap();
        assert!(result.has_warnings());
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.123dup");
    }

    #[test]
    fn lenient_warns_and_collapses_inv() {
        let result = parse_hgvs_lenient("NM_000088.3:c.100_100inv").unwrap();
        assert!(result.has_warnings());
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.100inv");
    }

    #[test]
    fn strict_rejects() {
        for input in [
            "NM_000088.3:c.123_123del",
            "NM_000088.3:c.123_123dup",
            "NM_000088.3:c.100_100inv",
        ] {
            let result = parse_hgvs_with_config(input, ErrorConfig::strict());
            assert!(result.is_err(), "strict should reject {}", input);
        }
    }

    #[test]
    fn silent_corrects_without_warning() {
        let result = parse_hgvs_silent("NM_000088.3:c.123_123del").unwrap();
        assert!(!result.has_warnings());
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.123del");
    }

    #[test]
    fn canonical_input_does_not_warn() {
        for input in [
            "NM_000088.3:c.123_126del",
            "NM_000088.3:c.123_126dup",
            "NM_000088.3:c.100_102inv",
        ] {
            let result = parse_hgvs_lenient(input).unwrap();
            assert!(
                !result
                    .warnings
                    .iter()
                    .any(|w| w.error_type == ErrorType::SinglePositionRange),
                "canonical input {} should not warn",
                input
            );
        }
    }

    #[test]
    fn compound_allele_two_warnings() {
        let result = parse_hgvs_lenient("NM_000088.3:c.[100_100del;200_200dup]").unwrap();
        let n = result
            .warnings
            .iter()
            .filter(|w| w.error_type == ErrorType::SinglePositionRange)
            .count();
        assert_eq!(n, 2);
    }

    #[test]
    fn idempotent_re_pass() {
        let r1 = parse_hgvs_lenient("NM_000088.3:c.123_123del").unwrap();
        let r2 = parse_hgvs_lenient(&r1.preprocessed_input).unwrap();
        assert_eq!(r2.preprocessed_input, "NM_000088.3:c.123del");
        assert!(!r2
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::SinglePositionRange));
    }
}

mod w3012_empty_delins_emission {
    use super::*;
    use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_silent, parse_hgvs_with_config};

    #[test]
    fn lenient_warns_and_rewrites_to_del() {
        let result = parse_hgvs_lenient("NC_000001.11:g.100_102delins").unwrap();
        assert!(result.has_warnings());
        assert_eq!(result.preprocessed_input, "NC_000001.11:g.100_102del");
        assert_eq!(
            result
                .warnings
                .iter()
                .filter(|w| w.error_type == ErrorType::EmptyDelinsInsert)
                .count(),
            1
        );
    }

    #[test]
    fn strict_rejects() {
        let result = parse_hgvs_with_config("NC_000001.11:g.100_102delins", ErrorConfig::strict());
        assert!(result.is_err());
    }

    #[test]
    fn silent_corrects_without_warning() {
        let result = parse_hgvs_silent("NC_000001.11:g.100_102delins").unwrap();
        assert!(!result.has_warnings());
        assert_eq!(result.preprocessed_input, "NC_000001.11:g.100_102del");
    }

    #[test]
    fn delins_with_payload_does_not_warn() {
        for input in [
            "NC_000001.11:g.100_102delinsATG",
            "NC_000001.11:g.100_102delinsN[12]",
            "NC_000001.11:g.100_102delins[A;G]",
        ] {
            let result = parse_hgvs_lenient(input).unwrap();
            assert!(
                !result
                    .warnings
                    .iter()
                    .any(|w| w.error_type == ErrorType::EmptyDelinsInsert),
                "{} should not warn",
                input
            );
        }
    }

    #[test]
    fn idempotent_re_pass() {
        let r1 = parse_hgvs_lenient("NC_000001.11:g.100_102delins").unwrap();
        let r2 = parse_hgvs_lenient(&r1.preprocessed_input).unwrap();
        assert_eq!(r2.preprocessed_input, "NC_000001.11:g.100_102del");
        assert!(!r2
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::EmptyDelinsInsert));
    }
}

mod w3013_redundant_repeat_label_emission {
    use super::*;
    use ferro_hgvs::hgvs::parser::{parse_hgvs_lenient, parse_hgvs_silent, parse_hgvs_with_config};

    #[test]
    fn lenient_warns_and_strips_label() {
        let result = parse_hgvs_lenient("NM_000088.3:r.100_102cug[4]").unwrap();
        assert!(result.has_warnings());
        assert_eq!(result.preprocessed_input, "NM_000088.3:r.100_102[4]");
        assert_eq!(
            result
                .warnings
                .iter()
                .filter(|w| w.error_type == ErrorType::RedundantRepeatLabel)
                .count(),
            1
        );
    }

    #[test]
    fn strict_rejects() {
        let result = parse_hgvs_with_config("NM_000088.3:r.100_102cug[4]", ErrorConfig::strict());
        assert!(result.is_err());
    }

    #[test]
    fn silent_strips_label_without_warning() {
        let result = parse_hgvs_silent("NM_000088.3:r.100_102cug[4]").unwrap();
        assert!(!result.has_warnings());
        assert_eq!(result.preprocessed_input, "NM_000088.3:r.100_102[4]");
    }

    #[test]
    fn canonical_input_does_not_warn() {
        let result = parse_hgvs_lenient("NM_000088.3:r.100_102[4]").unwrap();
        assert!(!result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::RedundantRepeatLabel));
    }

    #[test]
    fn idempotent_re_pass() {
        let r1 = parse_hgvs_lenient("NM_000088.3:r.100_102cug[4]").unwrap();
        let r2 = parse_hgvs_lenient(&r1.preprocessed_input).unwrap();
        assert_eq!(r2.preprocessed_input, "NM_000088.3:r.100_102[4]");
        assert!(!r2
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::RedundantRepeatLabel));
    }
}

// ============================================================================
// Issue #115 — spec-mandated input rejections
// ============================================================================

mod spec_mandated_rejections {
    use super::*;
    use ferro_hgvs::error::ErrorCode;
    use ferro_hgvs::hgvs::parser::parse_hgvs_lenient;
    use ferro_hgvs::parse_hgvs;

    // ----- W3003 OldSubstitutionSyntax (multi-base subs) -----

    #[test]
    fn w3003_lenient_rewrites_with_ref_bases_and_warns() {
        let result = parse_hgvs_lenient("NM_000088.3:c.79_80GC>TT").unwrap();
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.79_80delinsTT");
        assert!(result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::OldSubstitutionSyntax));
    }

    #[test]
    fn w3003_lenient_rewrites_no_ref_bases() {
        let result = parse_hgvs_lenient("NM_000088.3:c.100_102>ATG").unwrap();
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.100_102delinsATG");
    }

    #[test]
    fn w3003_strict_rejects_multi_base_substitution() {
        let config = ErrorConfig::strict();
        let pp = config.preprocessor();
        let result = pp.preprocess("NM_000088.3:c.79_80GC>TT");
        assert!(!result.success);
        let err = result.error.expect("strict mode must populate error");
        assert_eq!(err.code(), Some(ErrorCode::InvalidEdit));
    }

    #[test]
    fn w3003_canonical_substitution_unaffected() {
        let result = parse_hgvs_lenient("NM_000088.3:c.100A>G").unwrap();
        assert!(result.warnings.is_empty());
        assert_eq!(result.preprocessed_input, "NM_000088.3:c.100A>G");
    }

    // ----- W3014 DeprecatedIvsNotation -----

    #[test]
    fn w3014_strict_rejects_ivs() {
        let config = ErrorConfig::strict();
        let pp = config.preprocessor();
        let result = pp.preprocess("NM_000088.3:c.IVS2+2T>G");
        assert!(!result.success);
        let err = result.error.expect("strict mode must populate error");
        assert_eq!(err.code(), Some(ErrorCode::UnexpectedChar));
    }

    #[test]
    fn w3014_lenient_rejects_ivs_no_autocorrect() {
        let err = parse_hgvs_lenient("NM_000088.3:c.IVS2+2T>G")
            .expect_err("IVS cannot be auto-corrected");
        assert_eq!(err.code(), Some(ErrorCode::UnexpectedChar));
    }

    #[test]
    fn w3014_canonical_intronic_unaffected() {
        let result = parse_hgvs_lenient("NM_000088.3:c.88+2T>G");
        assert!(result.is_ok());
    }

    // ----- W3015 DeprecatedConSyntax -----

    #[test]
    fn w3015_lenient_rewrites_con_to_delins() {
        let result = parse_hgvs_lenient("NM_004006.2:c.100_200conNM_001.1:c.5_105").unwrap();
        assert_eq!(
            result.preprocessed_input,
            "NM_004006.2:c.100_200delinsNM_001.1:c.5_105"
        );
        assert!(result
            .warnings
            .iter()
            .any(|w| w.error_type == ErrorType::DeprecatedConSyntax));
    }

    #[test]
    fn w3015_strict_rejects_con() {
        let config = ErrorConfig::strict();
        let pp = config.preprocessor();
        let result = pp.preprocess("NM_004006.2:c.100_200conNM_001.1:c.5_105");
        assert!(!result.success);
        let err = result.error.expect("strict mode must populate error");
        assert_eq!(err.code(), Some(ErrorCode::InvalidEdit));
    }

    // ----- E3006 SelfCancellingAllele -----

    #[test]
    fn e3006_rejects_spec_example_in_strict() {
        let err = parse_hgvs("NM_004006.2:c.[762_768del;767_774dup]")
            .expect_err("self-cancelling allele must be rejected in strict mode");
        assert_eq!(err.code(), Some(ErrorCode::SelfCancellingAllele));
    }

    #[test]
    fn e3006_rejects_spec_example_in_lenient() {
        // E-codes are unconditional rejections; lenient mode is irrelevant.
        let err = parse_hgvs_lenient("NM_004006.2:c.[762_768del;767_774dup]")
            .expect_err("self-cancelling allele must be rejected in lenient mode");
        assert_eq!(err.code(), Some(ErrorCode::SelfCancellingAllele));
    }

    #[test]
    fn e3006_non_overlapping_allele_accepted() {
        let result = parse_hgvs("NM_004006.2:c.[100_110del;200_210dup]");
        assert!(result.is_ok());
    }

    #[test]
    fn e3006_trans_phase_overlap_accepted() {
        // Trans-phase alleles describe variants on different chromosomes;
        // a del on one allele cannot self-cancel a dup on the other.
        let result = parse_hgvs("NM_004006.2:c.[762_768del];[767_774dup]");
        assert!(
            result.is_ok(),
            "trans-phase overlap should not trigger E3006: {:?}",
            result
        );
    }

    // ----- Registry sanity -----

    #[test]
    fn registry_exposes_all_new_codes() {
        assert!(get_code_info("E3006").is_some());
        assert!(get_code_info("W3014").is_some());
        assert!(get_code_info("W3015").is_some());
    }
}
