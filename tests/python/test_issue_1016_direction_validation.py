"""Tests for issue #1016: the ``direction`` argument is validated instead of
silently defaulting to 3'.

Previously any unrecognized shuffle-direction string (``"5prim"``, ``"five"``,
``"3prine"``, ...) was silently mapped to 3', quietly changing normalization
output. The only documented values are ``"3prime"`` and ``"5prime"`` (plus the
``3``/``5``/``3'``/``5'`` short forms); anything else must now raise
``ValueError`` at the Python boundary, mirroring the ``assembly`` argument's
validate-and-raise behavior.
"""

from pathlib import Path

import pytest

import ferro_hgvs

# A variant that parses and normalizes cleanly against the built-in test data
# (bundled ``NM_000088.3`` transcript), so a successful call reflects an
# accepted direction rather than a coincidental normalization failure.
_VALID_HGVS = "NM_000088.3:c.4C>G"

# Every documented/accepted spelling, case-insensitively.
_ACCEPTED_DIRECTIONS = [
    "3prime",
    "5prime",
    "3",
    "5",
    "3'",
    "5'",
    "3PRIME",
    "5Prime",
]

# Representative unrecognized spellings (typos and unsupported words).
_UNRECOGNIZED_DIRECTIONS = ["5prim", "five", "3prine", "", "left", "threeprime"]


class TestNormalizerConstructor:
    """The ``Normalizer`` constructor validates ``direction``."""

    @pytest.mark.parametrize("direction", _UNRECOGNIZED_DIRECTIONS)
    def test_unrecognized_direction_raises(self, direction: str) -> None:
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.Normalizer(direction=direction)

    @pytest.mark.parametrize("direction", _ACCEPTED_DIRECTIONS)
    def test_accepted_direction_ok(self, direction: str) -> None:
        # Construction must succeed for every documented alias.
        ferro_hgvs.Normalizer(direction=direction)


class TestModuleLevelNormalize:
    """The module-level ``normalize()`` free function validates ``direction``."""

    @pytest.mark.parametrize("direction", _UNRECOGNIZED_DIRECTIONS)
    def test_unrecognized_direction_raises(self, direction: str) -> None:
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.normalize(_VALID_HGVS, direction=direction)

    @pytest.mark.parametrize("direction", _ACCEPTED_DIRECTIONS)
    def test_accepted_direction_ok(self, direction: str) -> None:
        result = ferro_hgvs.normalize(_VALID_HGVS, direction=direction)
        assert isinstance(result, str)


class TestHgvsVariantNormalize:
    """``HgvsVariant.normalize()`` validates ``direction``."""

    @pytest.mark.parametrize("direction", _UNRECOGNIZED_DIRECTIONS)
    def test_unrecognized_direction_raises(self, direction: str) -> None:
        variant = ferro_hgvs.parse(_VALID_HGVS)
        with pytest.raises(ValueError, match="unrecognized direction"):
            variant.normalize(direction=direction)

    @pytest.mark.parametrize("direction", _ACCEPTED_DIRECTIONS)
    def test_accepted_direction_ok(self, direction: str) -> None:
        variant = ferro_hgvs.parse(_VALID_HGVS)
        result = variant.normalize(direction=direction)
        # A recognized direction normalizes to a non-empty rendered variant
        # (`isinstance(str(x), str)` is vacuously true — assert real content).
        assert str(result)


class TestVariantProjectorConstructor:
    """The ``VariantProjector`` constructor validates ``direction`` too — it is
    the fifth Python entry point that routes through the same boundary helper."""

    @pytest.mark.parametrize("direction", _UNRECOGNIZED_DIRECTIONS)
    def test_unrecognized_direction_raises(self, direction: str) -> None:
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.VariantProjector(direction=direction)

    @pytest.mark.parametrize("direction", _ACCEPTED_DIRECTIONS)
    def test_accepted_direction_ok(self, direction: str) -> None:
        # Constructs against built-in test data; a recognized direction must not
        # raise (construction alone exercises the boundary validation).
        ferro_hgvs.VariantProjector(direction=direction)


# The two `from_manifest` entry points route through the same boundary helper, so
# an invalid direction must raise there too — and, crucially, BEFORE the manifest
# is loaded, so the reject path needs no real reference data.
MANIFEST_TINY = (
    Path(__file__).parent.parent / "fixtures" / "python" / "manifest_tiny" / "manifest.json"
)


class TestFromManifestDirection:
    """`Normalizer.from_manifest` / `VariantProjector.from_manifest` validate
    ``direction`` — the two entry points that were otherwise uncovered."""

    @pytest.mark.parametrize("direction", _UNRECOGNIZED_DIRECTIONS)
    def test_normalizer_from_manifest_unrecognized_raises(self, direction: str) -> None:
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY), direction=direction)

    @pytest.mark.parametrize("direction", _UNRECOGNIZED_DIRECTIONS)
    def test_projector_from_manifest_unrecognized_raises(self, direction: str) -> None:
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.VariantProjector.from_manifest(str(MANIFEST_TINY), direction=direction)

    def test_normalizer_from_manifest_accepted_ok(self) -> None:
        # A recognized direction loads the (tiny) manifest without raising.
        ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY), direction="5prime")


class TestDirectionValidatedBeforeReferenceLoad:
    """``direction`` (and ``assembly``) are validated BEFORE any reference is
    loaded, so a bad/missing reference path cannot mask an invalid argument: the
    boundary ``ValueError`` must win over the load's reference error. Uses a
    nonexistent path so the load *would* fail if it were reached first."""

    def test_normalizer_new_bad_json_and_bad_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.Normalizer(reference_json=str(missing), direction="sideways")

    def test_normalizer_from_manifest_bad_path_and_bad_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.Normalizer.from_manifest(str(missing), direction="sideways")

    def test_projector_new_bad_json_and_bad_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.VariantProjector(reference_json=str(missing), direction="sideways")

    def test_projector_from_manifest_bad_path_and_bad_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(ValueError, match="unrecognized direction"):
            ferro_hgvs.VariantProjector.from_manifest(str(missing), direction="sideways")

    def test_projector_bad_assembly_validated_before_bad_path(self, tmp_path: Path) -> None:
        # The `assembly` argument carries the same before-load guarantee.
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(ValueError, match="unrecognized assembly"):
            ferro_hgvs.VariantProjector.from_manifest(str(missing), assembly="bogus")
