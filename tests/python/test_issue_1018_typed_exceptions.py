"""Tests for issue #1018: a typed exception hierarchy at the Python boundary.

Previously every ferro failure surfaced as a bare ``ValueError`` or
``RuntimeError`` with a stringified message, discarding the structured
``E####`` / ``W####`` error code and the equivalent mutalyzer diagnostic codes
that ``FerroError`` carries internally.

Now every variant-processing failure is a ``ferro_hgvs.FerroError`` (or a
subclass: ``ParseError``, ``NormalizationError``, ``ReferenceDataError``,
``ProjectionError``), each exposing ``.code`` and ``.mutalyzer_codes``. To stay
non-breaking, each subclass also inherits the built-in exception its call sites
historically raised, so existing ``except ValueError`` / ``except RuntimeError``
callers keep working. Pure argument-validation failures stay plain
``ValueError`` and are intentionally outside the hierarchy.
"""

import copy
import pickle

import pytest

import ferro_hgvs


class TestHierarchyShape:
    """The class relationships are what callers catch against."""

    def test_subclasses_inherit_ferro_error(self) -> None:
        for exc in (
            ferro_hgvs.ParseError,
            ferro_hgvs.NormalizationError,
            ferro_hgvs.ReferenceDataError,
            ferro_hgvs.ProjectionError,
            ferro_hgvs.EquivalenceError,
        ):
            assert issubclass(exc, ferro_hgvs.FerroError)

    def test_builtin_backcompat_bases(self) -> None:
        # Each subclass also inherits the built-in it historically raised, so
        # existing except clauses keep working.
        assert issubclass(ferro_hgvs.ParseError, ValueError)
        assert issubclass(ferro_hgvs.NormalizationError, RuntimeError)
        # Reference and projection failures spanned both ValueError and
        # RuntimeError across their call sites, so they inherit both.
        assert issubclass(ferro_hgvs.ReferenceDataError, ValueError)
        assert issubclass(ferro_hgvs.ReferenceDataError, RuntimeError)
        assert issubclass(ferro_hgvs.ProjectionError, ValueError)
        assert issubclass(ferro_hgvs.ProjectionError, RuntimeError)
        # Equivalence checks historically raised RuntimeError; the typed
        # EquivalenceError must keep catching under ``except RuntimeError``.
        assert issubclass(ferro_hgvs.EquivalenceError, RuntimeError)

    def test_base_is_plain_exception(self) -> None:
        # The base is not a ValueError/RuntimeError itself; only the subclasses
        # opt into that back-compat.
        assert issubclass(ferro_hgvs.FerroError, Exception)
        assert not issubclass(ferro_hgvs.FerroError, ValueError)
        assert not issubclass(ferro_hgvs.FerroError, RuntimeError)


class TestParseError:
    """``parse`` raises a typed ``ParseError``."""

    def test_parse_raises_parse_error(self) -> None:
        with pytest.raises(ferro_hgvs.ParseError):
            ferro_hgvs.parse("this is not a variant")

    def test_parse_error_is_ferro_error(self) -> None:
        with pytest.raises(ferro_hgvs.FerroError):
            ferro_hgvs.parse("this is not a variant")

    def test_parse_error_backcompat_value_error(self) -> None:
        # Existing callers catching ValueError still work.
        with pytest.raises(ValueError):
            ferro_hgvs.parse("this is not a variant")

    def test_parse_error_carries_mutalyzer_codes(self) -> None:
        with pytest.raises(ferro_hgvs.ParseError) as info:
            ferro_hgvs.parse("this is not a variant")
        exc = info.value
        # The underlying FerroError::Parse maps to a non-empty set of
        # mutalyzer codes.
        assert isinstance(exc.mutalyzer_codes, tuple)
        assert len(exc.mutalyzer_codes) > 0
        assert all(isinstance(code, str) for code in exc.mutalyzer_codes)

    def test_str_is_message_only(self) -> None:
        # The structured fields live on attributes; str() is just the message,
        # with no tuple/None leakage.
        with pytest.raises(ferro_hgvs.ParseError) as info:
            ferro_hgvs.parse("this is not a variant")
        assert str(info.value).startswith("Parse error:")


class TestReferenceDataError:
    """Reference lookups raise a typed ``ReferenceDataError`` with a code."""

    def test_unknown_transcript_raises_reference_error(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        with pytest.raises(ferro_hgvs.ReferenceDataError):
            mapper.c_to_n("NONEXISTENT.1", 100)

    def test_reference_error_carries_structured_code(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        with pytest.raises(ferro_hgvs.ReferenceDataError) as info:
            mapper.c_to_n("NONEXISTENT.1", 100)
        exc = info.value
        # "reference not found" is ErrorCode 2001.
        assert exc.code == "E2001"
        assert isinstance(exc.mutalyzer_codes, tuple)

    def test_reference_error_backcompat(self) -> None:
        # Historically "Transcript not found" was a ValueError; the load/prepare
        # sites were RuntimeError. ReferenceDataError satisfies both.
        mapper = ferro_hgvs.CoordinateMapper()
        with pytest.raises(ValueError):
            mapper.c_to_n("NONEXISTENT.1", 100)


class TestNormalizationError:
    """Normalization failures raise a typed ``NormalizationError`` end-to-end.

    Guards the ``ferro_typed("NormalizationError", …)`` raise sites: the class
    name is a bare string literal, so without an end-to-end test a typo (or a
    future rename) would silently fall through to the plain-``RuntimeError``
    fallback in ``ferro_exception`` and go unnoticed.
    """

    # A protein deletion that cannot be normalized against the built-in data.
    _FAILS = "NP_000079.2:p.Gly2del"

    def test_raises_normalization_error(self) -> None:
        with pytest.raises(ferro_hgvs.NormalizationError):
            ferro_hgvs.normalize(self._FAILS)

    def test_backcompat_runtime_error(self) -> None:
        # Historically a normalization failure surfaced as ``RuntimeError``.
        with pytest.raises(RuntimeError):
            ferro_hgvs.normalize(self._FAILS)

    def test_is_ferro_error_but_not_value_error(self) -> None:
        with pytest.raises(ferro_hgvs.FerroError) as info:
            ferro_hgvs.normalize(self._FAILS)
        # NormalizationError deliberately does NOT inherit ValueError.
        assert not isinstance(info.value, ValueError)


class TestProjectionError:
    """Coordinate conversions raise a typed ``ProjectionError`` end-to-end."""

    def test_raises_projection_error(self) -> None:
        mapper = ferro_hgvs.CoordinateMapper()
        # The transcript resolves, but CDS position 0 is not a valid coordinate,
        # so the conversion (not the reference lookup) fails.
        with pytest.raises(ferro_hgvs.ProjectionError):
            mapper.c_to_p("NM_000088.3", 0)

    def test_backcompat_both_bases(self) -> None:
        # Conversion sites historically spanned ValueError and RuntimeError, so
        # ProjectionError must satisfy either ``except`` clause.
        mapper = ferro_hgvs.CoordinateMapper()
        with pytest.raises(ferro_hgvs.ProjectionError) as info:
            mapper.c_to_p("NM_000088.3", 0)
        assert isinstance(info.value, ValueError)
        assert isinstance(info.value, RuntimeError)


class TestEquivalenceError:
    """Equivalence checks raise a typed ``EquivalenceError`` end-to-end."""

    def _incompatible_pair(self) -> tuple[object, object]:
        # A protein variant and a coding variant cannot be compared for
        # equivalence, so the check itself errors.
        return (
            ferro_hgvs.parse("NP_000079.2:p.Gly2Ala"),
            ferro_hgvs.parse("NM_000088.3:c.5G>C"),
        )

    def test_raises_equivalence_error(self) -> None:
        checker = ferro_hgvs.EquivalenceChecker()
        v1, v2 = self._incompatible_pair()
        with pytest.raises(ferro_hgvs.EquivalenceError):
            checker.check(v1, v2)

    def test_backcompat_runtime_error(self) -> None:
        checker = ferro_hgvs.EquivalenceChecker()
        v1, v2 = self._incompatible_pair()
        with pytest.raises(RuntimeError) as info:
            checker.check(v1, v2)
        # EquivalenceError inherits RuntimeError but not ValueError.
        assert not isinstance(info.value, ValueError)


class TestArgumentValidationStaysPlain:
    """Pure argument-validation errors are plain ValueError, not FerroError."""

    def test_unrecognized_direction_is_plain_value_error(self) -> None:
        with pytest.raises(ValueError) as info:
            ferro_hgvs.normalize("NM_000088.3:c.100A>G", direction="sideways")
        assert not isinstance(info.value, ferro_hgvs.FerroError)

    def test_empty_ref_base_is_plain_value_error(self) -> None:
        with pytest.raises(ValueError) as info:
            ferro_hgvs.VcfRecord.snv("chr1", 12345, "", "G")
        assert not isinstance(info.value, ferro_hgvs.FerroError)


class TestManualConstruction:
    """Constructing an exception directly yields sensible defaults."""

    def test_defaults(self) -> None:
        exc = ferro_hgvs.ParseError("boom")
        assert str(exc) == "boom"
        assert exc.code is None
        assert exc.mutalyzer_codes == ()

    def test_explicit_fields(self) -> None:
        exc = ferro_hgvs.NormalizationError("bad", code="E4001", mutalyzer_codes=["EINTRONIC"])
        assert exc.code == "E4001"
        assert exc.mutalyzer_codes == ("EINTRONIC",)


class TestPreservesStructuredFieldsAcrossPickleAndCopy:
    """``code`` / ``mutalyzer_codes`` survive pickling and deep-copying.

    ``BaseException.__init__`` only seeds ``args``, so without the custom
    ``__reduce__`` on ``FerroError`` the default pickler drops the structured
    fields — silently, on any cross-process (e.g. ``multiprocessing``) hop.
    """

    def _sample(self) -> "ferro_hgvs.ReferenceDataError":
        return ferro_hgvs.ReferenceDataError(
            "boom", code="E2001", mutalyzer_codes=["EINTRONIC", "EOUTOFBOUNDARY"]
        )

    def test_pickle_round_trip_preserves_fields(self) -> None:
        original = self._sample()
        restored = pickle.loads(pickle.dumps(original))
        assert type(restored) is ferro_hgvs.ReferenceDataError
        assert str(restored) == "boom"
        assert restored.code == "E2001"
        assert restored.mutalyzer_codes == ("EINTRONIC", "EOUTOFBOUNDARY")
        # Back-compat bases survive too.
        assert isinstance(restored, ValueError)
        assert isinstance(restored, RuntimeError)

    def test_deepcopy_preserves_fields(self) -> None:
        original = self._sample()
        restored = copy.deepcopy(original)
        assert type(restored) is ferro_hgvs.ReferenceDataError
        assert restored.code == "E2001"
        assert restored.mutalyzer_codes == ("EINTRONIC", "EOUTOFBOUNDARY")

    def test_default_fields_round_trip(self) -> None:
        # A bare exception (no code / codes) must also reconstruct cleanly.
        restored = pickle.loads(pickle.dumps(ferro_hgvs.ParseError("x")))
        assert type(restored) is ferro_hgvs.ParseError
        assert restored.code is None
        assert restored.mutalyzer_codes == ()
