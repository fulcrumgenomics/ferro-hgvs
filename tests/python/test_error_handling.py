"""Tests for error handling functionality."""

import re
from pathlib import Path

import ferro_hgvs

# Location of the committed type stub, relative to this test file
# (<repo>/tests/python/test_error_handling.py -> <repo>/python/ferro_hgvs/__init__.pyi).
_STUB_PATH = Path(__file__).resolve().parents[2] / "python" / "ferro_hgvs" / "__init__.pyi"


def _stub_enum_members(class_name: str) -> dict[str, int]:
    """Parse ``name = value`` members of an ``IntEnum`` class from the type stub.

    Reads the block following ``class <class_name>(IntEnum):`` up to the next
    top-level ``class`` declaration, returning the declared name->value mapping.
    Comment and docstring lines are ignored, so retired discriminants documented
    only in comments are not treated as members.
    """
    text = _STUB_PATH.read_text()
    lines = text.splitlines()
    header = re.compile(rf"^class {re.escape(class_name)}\(IntEnum\):")
    member = re.compile(r"^    (\w+) = (\d+)$")
    members: dict[str, int] = {}
    in_block = False
    for line in lines:
        if header.match(line):
            in_block = True
            continue
        if in_block:
            if line.startswith("class "):
                break
            match = member.match(line)
            if match:
                members[match.group(1)] = int(match.group(2))
    return members


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

    def test_error_type_stub_matches_runtime(self) -> None:
        """The ``ErrorType`` type stub must list exactly the runtime members.

        Guards against the stub drifting out of sync with the Rust
        ``PyErrorType`` enum (issue #630: discriminants 41/42 were missing from
        the stub while present at runtime). CI builds the native extension, so
        the runtime enum is authoritative here.

        ``ErrorType`` is a PyO3 enum, not a Python ``IntEnum``, so it is not
        iterable; its variants are exposed as ``int``-valued class attributes.
        """
        error_type = ferro_hgvs.ErrorType
        runtime = {
            name: int(getattr(error_type, name))
            for name in dir(error_type)
            if not name.startswith("_") and isinstance(getattr(error_type, name), error_type)
        }
        stub = _stub_enum_members("ErrorType")
        assert stub == runtime

    def test_error_type_centromere_and_flank_members(self) -> None:
        """Discriminants 41/42 are reachable as named members (issue #630)."""
        assert int(ferro_hgvs.ErrorType.UnresolvableCentromere) == 41
        assert int(ferro_hgvs.ErrorType.TranscriptFlankNotDescribable) == 42
