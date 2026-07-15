"""Tests for issue #1018: batch error-isolation on the projector.

``VariantProjector.project_many`` / ``project_normalized_many`` historically
aborted the whole batch on the first bad input, raising ``ProjectionError`` and
discarding every other result. The new ``return_exceptions=True`` mode (mirroring
``asyncio.gather(return_exceptions=True)``) instead returns, per input and in
order, either that input's list of projections or the ``ProjectionError``
describing its failure — so one malformed variant no longer loses the batch.
The default (``return_exceptions=False``) keeps the abort-and-raise behavior.
"""

import warnings

import pytest

import ferro_hgvs

# A syntactically valid g. variant (projects to a possibly-empty list — a
# success, not an error) and an unparseable string (fails inside project_all).
GOOD = "NC_000017.11:g.50198002C>T"
BAD = "totally not a variant"

# A parseable coding variant: it is NOT a genomic variant, so the
# already-normalized projection path (`project_normalized_all`) fails on it —
# used to exercise the isolation/fail-fast branches of `project_normalized_many`.
BAD_NORMALIZED = "NM_000088.3:c.4C>G"


def _projector() -> ferro_hgvs.VariantProjector:
    # Built-in test data has no genome; that emits a reduced-capability warning
    # (deduplicated by Python's warning registry across the session, so it is
    # not asserted here) and does not affect batch error-isolation semantics.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return ferro_hgvs.VariantProjector()


class TestProjectManyDefault:
    """The default keeps the historical abort-on-first-error behavior."""

    def test_default_aborts_and_raises_on_bad_input(self) -> None:
        proj = _projector()
        with pytest.raises(ferro_hgvs.ProjectionError):
            proj.project_many([GOOD, BAD, GOOD])

    def test_default_returns_list_of_lists_when_all_good(self) -> None:
        proj = _projector()
        out = proj.project_many([GOOD, GOOD])
        assert isinstance(out, list)
        assert all(isinstance(inner, list) for inner in out)


class TestProjectManyReturnExceptions:
    """``return_exceptions=True`` isolates failures per input."""

    def test_isolates_failures_in_order(self) -> None:
        proj = _projector()
        out = proj.project_many([GOOD, BAD, GOOD], return_exceptions=True)
        assert len(out) == 3
        assert isinstance(out[0], list)
        assert isinstance(out[1], ferro_hgvs.ProjectionError)
        assert isinstance(out[2], list)

    def test_returned_exception_is_typed_and_backcompat(self) -> None:
        proj = _projector()
        out = proj.project_many([BAD], return_exceptions=True)
        exc = out[0]
        # Returned in place, not raised — but still a real ProjectionError, and
        # still a RuntimeError so existing except-clauses recognize it.
        assert isinstance(exc, ferro_hgvs.ProjectionError)
        assert isinstance(exc, RuntimeError)
        assert isinstance(exc, ferro_hgvs.FerroError)
        assert "input[0]" in str(exc)

    def test_returned_exception_preserves_structured_payload(self) -> None:
        # Returning the typed object (rather than a stringified error) exists so
        # the structured payload survives capture: mutalyzer_codes is a tuple and
        # code is either a populated E#### string or None — never lost.
        proj = _projector()
        exc = proj.project_many([BAD], return_exceptions=True)[0]
        assert isinstance(exc, ferro_hgvs.ProjectionError)
        assert isinstance(exc.mutalyzer_codes, tuple)
        assert exc.code is None or (isinstance(exc.code, str) and exc.code.startswith(("E", "W")))

    def test_all_good_returns_no_exceptions(self) -> None:
        proj = _projector()
        out = proj.project_many([GOOD, GOOD], return_exceptions=True)
        assert all(isinstance(inner, list) for inner in out)

    def test_all_fail_batch_isolates_every_input(self) -> None:
        proj = _projector()
        out = proj.project_many([BAD, BAD, BAD], return_exceptions=True)
        assert len(out) == 3
        assert all(isinstance(item, ferro_hgvs.ProjectionError) for item in out)
        # The index in each message tracks its position.
        assert [f"input[{i}]" in str(item) for i, item in enumerate(out)] == [True, True, True]

    def test_empty_batch_returns_empty_list(self) -> None:
        proj = _projector()
        assert proj.project_many([], return_exceptions=True) == []
        assert proj.project_many([]) == []


class TestProjectNormalizedManyReturnExceptions:
    """The already-normalized variant path isolates failures the same way."""

    def test_default_returns_list_of_lists(self) -> None:
        proj = _projector()
        variant = ferro_hgvs.parse(GOOD)
        out = proj.project_normalized_many([variant])
        assert isinstance(out, list)
        assert all(isinstance(inner, list) for inner in out)

    def test_return_exceptions_isolates(self) -> None:
        proj = _projector()
        variant = ferro_hgvs.parse(GOOD)
        out = proj.project_normalized_many([variant], return_exceptions=True)
        assert len(out) == 1
        assert isinstance(out[0], list)

    def test_return_exceptions_isolates_a_real_failure_in_order(self) -> None:
        # Exercises the ERROR-capture branch of project_normalized_many (not just
        # project_many): a coding variant cannot be projected as a genomic one, so
        # its slot must hold a ProjectionError while the good slots stay lists.
        proj = _projector()
        good = ferro_hgvs.parse(GOOD)
        bad = ferro_hgvs.parse(BAD_NORMALIZED)
        out = proj.project_normalized_many([good, bad, good], return_exceptions=True)
        assert len(out) == 3
        assert isinstance(out[0], list)
        assert isinstance(out[1], ferro_hgvs.ProjectionError)
        assert "input[1]" in str(out[1])
        assert isinstance(out[2], list)

    def test_default_aborts_and_raises_on_bad_input(self) -> None:
        # The fail-fast default must still raise on the normalized path.
        proj = _projector()
        good = ferro_hgvs.parse(GOOD)
        bad = ferro_hgvs.parse(BAD_NORMALIZED)
        with pytest.raises(ferro_hgvs.ProjectionError):
            proj.project_normalized_many([good, bad])
