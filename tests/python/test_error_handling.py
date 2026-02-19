"""Tests for error handling functionality."""

import ferro_hgvs


class TestErrorConfig:
    """Tests for ErrorConfig class."""

    def test_strict_config(self) -> None:
        config = ferro_hgvs.ErrorConfig.strict()
        assert config.mode == ferro_hgvs.ErrorMode.Strict

    def test_lenient_config(self) -> None:
        config = ferro_hgvs.ErrorConfig.lenient()
        assert config.mode == ferro_hgvs.ErrorMode.Lenient

    def test_silent_config(self) -> None:
        config = ferro_hgvs.ErrorConfig.silent()
        assert config.mode == ferro_hgvs.ErrorMode.Silent

    def test_with_override(self) -> None:
        config = ferro_hgvs.ErrorConfig.lenient()
        new_config = config.with_override(
            ferro_hgvs.ErrorType.LowercaseAminoAcid,
            ferro_hgvs.ErrorOverride.Reject,
        )
        assert new_config.should_reject(ferro_hgvs.ErrorType.LowercaseAminoAcid)


class TestParseLenient:
    """Tests for lenient parsing."""

    def test_parse_lenient_valid(self) -> None:
        result = ferro_hgvs.parse_lenient("NM_000088.3:c.100A>G")
        assert result.variant is not None
        assert not result.has_warnings()
        assert not result.had_corrections()

    def test_parse_lenient_with_corrections(self) -> None:
        # Lowercase amino acid should be auto-corrected in lenient mode
        result = ferro_hgvs.parse_lenient("NP_000079.2:p.glu6val")
        assert result.variant is not None
        # Check if corrections were made
        if result.had_corrections():
            assert result.preprocessed_input != result.original_input

    def test_parse_lenient_with_config(self) -> None:
        # Test that different configs can be passed to parse_lenient
        strict_config = ferro_hgvs.ErrorConfig.strict()
        lenient_config = ferro_hgvs.ErrorConfig.lenient()

        # Both should parse valid input
        result_strict = ferro_hgvs.parse_lenient("NM_000088.3:c.100A>G", strict_config)
        result_lenient = ferro_hgvs.parse_lenient("NM_000088.3:c.100A>G", lenient_config)
        assert result_strict.variant is not None
        assert result_lenient.variant is not None


class TestCorrectionWarning:
    """Tests for CorrectionWarning class."""

    def test_warning_attributes(self) -> None:
        # Parse something that generates warnings
        result = ferro_hgvs.parse_lenient("NP_000079.2:p.glu6val")
        if result.warnings:
            warning = result.warnings[0]
            assert hasattr(warning, "error_type")
            assert hasattr(warning, "message")
            assert hasattr(warning, "original")
            assert hasattr(warning, "corrected")


class TestErrorTypes:
    """Tests for error type enums."""

    def test_error_mode_values(self) -> None:
        assert ferro_hgvs.ErrorMode.Strict is not None
        assert ferro_hgvs.ErrorMode.Lenient is not None
        assert ferro_hgvs.ErrorMode.Silent is not None

    def test_error_override_values(self) -> None:
        assert ferro_hgvs.ErrorOverride.Default is not None
        assert ferro_hgvs.ErrorOverride.Reject is not None
        assert ferro_hgvs.ErrorOverride.WarnCorrect is not None
        assert ferro_hgvs.ErrorOverride.SilentCorrect is not None
        assert ferro_hgvs.ErrorOverride.Accept is not None
