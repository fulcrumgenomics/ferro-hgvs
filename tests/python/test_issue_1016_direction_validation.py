"""Tests for issue #1016, as its property survives the removal of the public 5'
surface: no ``direction`` argument may be accepted-and-ignored.

#1016 filed the original shape — any unrecognized shuffle-direction string
(``"5prim"``, ``"five"``, ``"3prine"``, ...) was silently mapped to 3', quietly
changing normalization output — and #1017 fixed it by raising ``ValueError`` at
the Python boundary.

The ``direction=`` keyword has since been **removed** from every Python entry
point. ``README.md`` rule 6 says there are no user options for normalization
form, and a shuffle direction is not orthogonal to the form: it selects the
frame every other rule is evaluated in. ferro normalizes 3', the only direction
the HGVS recommendations describe.

The removal has to preserve #1016's property, and this file is what pins that.
The failure mode a removal invites is the one #1016 was about, only worse:
a keyword that is *quietly dropped* would make ``direction="5prime"`` a
**silent 3' run** at a caller who explicitly asked for 5' — a wrong answer
returned under an argument that appears to have been honoured. PyO3 raises
``TypeError`` for a keyword its signature does not name, so the argument cannot
be silently dropped; every entry point below is asserted to reject it.

Migration for a caller who was passing ``direction="3prime"``: drop the
argument, the behaviour is unchanged. For ``direction="5prime"``: there is no
replacement — that form is no longer produced by any ferro entry point.
"""

from pathlib import Path

import pytest

import ferro_hgvs

# A variant that parses and normalizes cleanly against the built-in test data
# (bundled ``NM_000088.3`` transcript), so a successful call reflects a working
# entry point rather than a coincidental normalization failure.
_VALID_HGVS = "NM_000088.3:c.4C>G"

# Every spelling the removed keyword used to accept, plus the typos #1016 was
# filed about. After removal all of them must fail identically — the keyword is
# rejected by name, so its *value* is never even examined. That uniformity is
# the point: a caller cannot stumble onto a spelling that is quietly tolerated.
_ANY_DIRECTION_VALUE = [
    "3prime",
    "5prime",
    "3",
    "5",
    "3'",
    "5'",
    "3PRIME",
    "5Prime",
    "5prim",
    "five",
    "3prine",
    "",
    "left",
    "threeprime",
    "sideways",
]

# PyO3's messages for a keyword a signature does not name. Matching them keeps
# the assertion honest: a bare ``TypeError`` could come from an arity mistake in
# the test itself, whereas these can only come from the argument being unknown.
#
# There are two forms, and both are measured rather than guessed. A signature
# that still names other keywords reports the offending one
# ("...got an unexpected keyword argument 'direction'"); a signature that names
# *no* keywords at all — ``HgvsVariant.normalize()``, which lost its only one —
# reports "takes no keyword arguments" instead, without naming it. The second is
# a stronger statement, not a weaker one: the method accepts no keyword
# whatsoever, so ``direction=`` cannot be quietly absorbed.
_UNEXPECTED_KWARG = "unexpected keyword argument 'direction'|takes no keyword arguments"


class TestNormalizerConstructor:
    """The ``Normalizer`` constructor rejects ``direction``."""

    @pytest.mark.parametrize("direction", _ANY_DIRECTION_VALUE)
    def test_direction_keyword_rejected(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.Normalizer(direction=direction)

    def test_constructs_without_the_keyword(self) -> None:
        ferro_hgvs.Normalizer()


class TestModuleLevelNormalize:
    """The module-level ``normalize()`` free function rejects ``direction``."""

    @pytest.mark.parametrize("direction", _ANY_DIRECTION_VALUE)
    def test_direction_keyword_rejected(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.normalize(_VALID_HGVS, direction=direction)

    def test_normalizes_without_the_keyword(self) -> None:
        result = ferro_hgvs.normalize(_VALID_HGVS)
        assert isinstance(result, str)
        assert result


class TestHgvsVariantNormalize:
    """``HgvsVariant.normalize()`` rejects ``direction``."""

    @pytest.mark.parametrize("direction", _ANY_DIRECTION_VALUE)
    def test_direction_keyword_rejected(self, direction: str) -> None:
        variant = ferro_hgvs.parse(_VALID_HGVS)
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            variant.normalize(direction=direction)

    def test_normalizes_without_the_keyword(self) -> None:
        variant = ferro_hgvs.parse(_VALID_HGVS)
        result = variant.normalize()
        # Assert real content — `isinstance(str(x), str)` is vacuously true.
        assert str(result)


class TestVariantProjectorConstructor:
    """The ``VariantProjector`` constructor rejects ``direction`` too — it was
    the fifth Python entry point routing through the removed boundary helper."""

    @pytest.mark.parametrize("direction", _ANY_DIRECTION_VALUE)
    def test_direction_keyword_rejected(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.VariantProjector(direction=direction)

    def test_constructs_without_the_keyword(self) -> None:
        # Constructs against built-in test data.
        ferro_hgvs.VariantProjector()


# The two `from_manifest` entry points routed through the same boundary helper,
# so they must reject the keyword too — and, as before, without needing real
# reference data, since an unknown keyword is refused before any body runs.
MANIFEST_TINY = (
    Path(__file__).parent.parent / "fixtures" / "python" / "manifest_tiny" / "manifest.json"
)


class TestFromManifestDirection:
    """`Normalizer.from_manifest` / `VariantProjector.from_manifest` reject
    ``direction`` — the two entry points that were otherwise uncovered."""

    @pytest.mark.parametrize("direction", _ANY_DIRECTION_VALUE)
    def test_normalizer_from_manifest_rejects(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY), direction=direction)

    @pytest.mark.parametrize("direction", _ANY_DIRECTION_VALUE)
    def test_projector_from_manifest_rejects(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.VariantProjector.from_manifest(str(MANIFEST_TINY), direction=direction)

    def test_normalizer_from_manifest_loads_without_the_keyword(self) -> None:
        ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY))


class TestFromSequencesEntryPoints:
    """The four sequence-first entry points carried the keyword as well."""

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_free_from_sequences_rejects(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.from_sequences("NC_000001.11", 100, "A", "AT", direction=direction)

    @pytest.mark.parametrize("direction", ["3prime", "5prime"])
    def test_free_from_sequences_detailed_rejects(self, direction: str) -> None:
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.from_sequences_detailed("NC_000001.11", 100, "A", "AT", direction=direction)

    def test_free_from_sequences_works_without_the_keyword(self) -> None:
        assert str(ferro_hgvs.from_sequences("NC_000001.11", 100, "A", "AT"))


class TestDirectionRejectedBeforeReferenceLoad:
    """The keyword is refused BEFORE any reference is loaded, so a bad/missing
    reference path cannot mask it: the ``TypeError`` must win over the load's
    reference error. Uses a nonexistent path so the load *would* fail if it were
    reached first. This is the guarantee #1016's fix established for the
    argument's ``ValueError``, carried over to its removal — a caller must never
    be told their path is wrong when what is actually wrong is that they asked
    for a direction that no longer exists."""

    def test_normalizer_new_bad_json_and_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.Normalizer(reference_json=str(missing), direction="5prime")

    def test_normalizer_from_manifest_bad_path_and_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.Normalizer.from_manifest(str(missing), direction="5prime")

    def test_projector_new_bad_json_and_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.VariantProjector(reference_json=str(missing), direction="5prime")

    def test_projector_from_manifest_bad_path_and_direction(self, tmp_path: Path) -> None:
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(TypeError, match=_UNEXPECTED_KWARG):
            ferro_hgvs.VariantProjector.from_manifest(str(missing), direction="5prime")

    def test_projector_bad_assembly_validated_before_bad_path(self, tmp_path: Path) -> None:
        # The `assembly` argument still carries the same before-load guarantee,
        # and still reports it as a `ValueError` — it is a live keyword whose
        # *value* was wrong, which is a different fault from a keyword that no
        # longer exists. Keeping this row here pins that the two stay distinct.
        missing = tmp_path / "does_not_exist.json"
        with pytest.raises(ValueError, match="unrecognized assembly"):
            ferro_hgvs.VariantProjector.from_manifest(str(missing), assembly="bogus")
