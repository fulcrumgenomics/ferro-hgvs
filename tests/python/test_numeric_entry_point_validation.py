"""``CoordinateMapper``'s numeric entry points must validate what they are handed.

Five methods take raw integers from Python and build position structs directly,
so none of the checks ``parse_hgvs`` applies to a written description reach
them. Three escapes follow, and all three are only observable from here — the
Rust API is reached through the parser, which refuses these values before a
position struct exists.

**Escape 1 — a base of zero.** ``parse_hgvs`` refuses ``c.0`` / ``n.0`` as
``E1003 InvalidPosition`` ("Position 0 is not valid in HGVS notation"), on the
spec's authority: ``background/numbering.md:31`` states there is no nucleotide
``c.0``. That is a claim about what the numbering axis contains, so it holds
however the coordinate arrives. Measured on the unfixed build, over the
two-exon transcript below::

    c_to_g("NM_TEST.1", 0) -> ('chr1', 1011)
    c_to_g("NM_TEST.1", 1) -> ('chr1', 1011)

— one entry point refusing and the other answering ``c.0`` *identically to*
``c.1``.

**Escape 2 — an unknown-offset sentinel used as a distance.** ``c.N+?`` / ``c.N-?``
parse to ``i64::MAX`` / ``i64::MIN``, which are markers for "unknown, unbounded
in this direction", not measured distances. Passed as ``offset`` they reach
``Transcript::intronic_to_genomic``'s unchecked arithmetic — the ``#1087`` /
PR ``#1088`` guard sits in ``src/data/mapping.rs`` and is not on this path.
Measured on the unfixed build::

    c_to_g("NM_TEST.1", 90, 2**63 - 1) -> ('chr1', 9223372036854776907)
    c_to_g("NM_TEST.1", 91, -2**63)    !! pyo3_runtime.PanicException:
                                          attempt to negate with overflow

The negative sentinel overflows in ``-offset`` *before* the cast, so it panics
under ``debug_assertions`` and wraps silently in the ``--release`` wheel CI
ships. The positive sentinel needs no overflow to do damage: it returns a
number that looks exactly like a coordinate.

**Escape 2b — the same class, one step away from the sentinel (#1765).** The
guard for escape 2 keys on sentinel *identity*, which is the right key for a
marker and the wrong one for the arithmetic. ``i64::MIN + 1`` and
``i64::MAX - 1`` are plain numbers, so no widening of the sentinel predicate
reaches them, and both reproduced escape 2 in full. The magnitude condition is
enforced where the magnitude can be judged — ``intronic_to_genomic`` is now
total (``unsigned_abs`` plus checked add/sub) and refuses a result outside the
intron the offset was measured into. See ``TestOversizedOffsetIsRefused``.

**Escape 3 — a value handed back that the caller cannot read.** ``c_to_n``
returns ``tx_pos.offset`` verbatim, so a sentinel comes back as
``9223372036854775807``, indistinguishable from a measured distance; and
``n_to_c`` forwards ``downstream``, which ``cds_to_tx``/``tx_to_cds`` discard,
so ``downstream=True`` and ``downstream=False`` returned the same triple.

The fix is at the boundary in every case: refuse the argument, or refuse to
return the value. Repairing the conversions themselves is separate work on
other branches, and a boundary that accepts a value it cannot mean is a defect
whichever way those land.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

import ferro_hgvs

# The parser's own sentinels, spelled as Python ints. ``c.N+?`` -> ``i64::MAX``,
# ``c.N-?`` -> ``i64::MIN`` (``src/hgvs/parser/position.rs``).
OFFSET_UNKNOWN_POSITIVE = 2**63 - 1
OFFSET_UNKNOWN_NEGATIVE = -(2**63)

# A two-exon plus-strand transcript with genomic coordinates, so ``c_to_g``
# resolves rather than declining for want of a genome. Built programmatically;
# nothing here is committed as data.
#
#   exon 1: tx   1..100, genomic 1001..1100
#   intron:               genomic 1101..2000
#   exon 2: tx 101..200, genomic 2001..2100
#   CDS:    tx  11..190  ->  c.1 = tx 11 = g.1011, c.90 = tx 100 = g.1100
#
# So ``c.90`` is the last exonic base before the intron and ``c.91`` the first
# after it: ``c.90+n`` and ``c.91-n`` are the two intronic arms whose arithmetic
# escape 2 is about.
TRANSCRIPT_ID = "NM_TEST.1"
REFERENCE = {
    "transcripts": [
        {
            "id": TRANSCRIPT_ID,
            "gene_symbol": "TEST",
            "strand": "+",
            "chromosome": "chr1",
            "genomic_start": 1001,
            "genomic_end": 2100,
            "cds_start": 11,
            "cds_end": 190,
            "exons": [
                {
                    "number": 1,
                    "start": 1,
                    "end": 100,
                    "genomic_start": 1001,
                    "genomic_end": 1100,
                },
                {
                    "number": 2,
                    "start": 101,
                    "end": 200,
                    "genomic_start": 2001,
                    "genomic_end": 2100,
                },
            ],
        }
    ]
}


@pytest.fixture
def mapper(tmp_path: Path) -> ferro_hgvs.CoordinateMapper:
    """A mapper over the two-exon transcript above.

    Construction may emit the reduced-capability ``UserWarning`` (the reference
    carries transcript geometry but no genome FASTA), which is expected and
    unrelated to anything under test. It is suppressed rather than asserted:
    the warn-once flag is process-global and shared across every surface, so
    only the first mapper built in a given interpreter emits it and a
    ``pytest.warns`` here would fail on every test after the first.
    """
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(REFERENCE))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return ferro_hgvs.CoordinateMapper(reference_json=str(path))


# ---------------------------------------------------------------------------
# Escape 1: a base of zero
# ---------------------------------------------------------------------------


class TestBaseZeroIsRefused:
    """``c.0`` / ``n.0`` get the parser's answer, not a coordinate."""

    @pytest.mark.parametrize(
        ("method", "args"),
        [
            ("c_to_g", ()),
            ("c_to_p", ()),
            ("c_to_n", ()),
            ("n_to_c", ()),
        ],
    )
    def test_zero_base_raises_e1003(
        self, mapper: ferro_hgvs.CoordinateMapper, method: str, args: tuple[object, ...]
    ) -> None:
        """Every numeric entry point refuses a base of zero, with the parser's code.

        The code is asserted, not just the exception type: ``c_to_p(0)`` already
        failed before this change, but with "Cannot convert UTR position to
        protein" — a wrong diagnosis that a bare ``pytest.raises`` would accept.
        """
        with pytest.raises(ferro_hgvs.ParseError) as excinfo:
            getattr(mapper, method)(TRANSCRIPT_ID, 0, *args)

        assert excinfo.value.code == "E1003"
        assert "Position 0 is not valid in HGVS notation" in str(excinfo.value)

    def test_the_written_and_the_numeric_entry_points_agree(
        self, mapper: ferro_hgvs.CoordinateMapper
    ) -> None:
        """One coordinate, one answer, whichever door it arrives through.

        This is the defect stated as a property: ``parse`` refuses ``c.0`` and
        the numeric entry point answered it, so ferro had two entry points and
        two answers. Only the *refusal* is compared, because the two doors
        genuinely differ in how much they know: ``parse`` reports
        ``ParseError`` with ``code is None`` (the nom parser bails before the
        preprocessor's position-zero phase, which is what attaches ``E1003``
        and what ``ferro parse`` prints on the CLI), while the numeric entry
        point can name the code exactly because it has the integer in hand.
        Pinning ``parse``'s missing code here would be pinning an unrelated
        inconsistency.
        """
        with pytest.raises(ferro_hgvs.ParseError):
            ferro_hgvs.parse(f"{TRANSCRIPT_ID}:c.0A>G")

        with pytest.raises(ferro_hgvs.ParseError):
            mapper.c_to_g(TRANSCRIPT_ID, 0)

    def test_nonzero_neighbours_still_convert(self, mapper: ferro_hgvs.CoordinateMapper) -> None:
        """The guard fires at zero and nowhere else.

        Without this, refusing every base would satisfy the assertions above.
        ``c.-1`` is the 5'UTR base immediately before ``c.1`` and must survive:
        it is a legal coordinate whose sign a naive "must be positive" check
        would reject.
        """
        assert mapper.c_to_g(TRANSCRIPT_ID, 1) == ("chr1", 1011)
        assert mapper.c_to_g(TRANSCRIPT_ID, -1) == ("chr1", 1010)
        assert mapper.c_to_p(TRANSCRIPT_ID, 1) == 1
        assert mapper.c_to_n(TRANSCRIPT_ID, 1) == (11, None)
        assert mapper.n_to_c(TRANSCRIPT_ID, 1) == (-10, None, False)


# ---------------------------------------------------------------------------
# Escape 2: an unknown-offset sentinel used as a distance
# ---------------------------------------------------------------------------


class TestSentinelOffsetIsRefused:
    """A sentinel is not a distance, so no coordinate may be derived from it."""

    @pytest.mark.parametrize("offset", [OFFSET_UNKNOWN_POSITIVE, OFFSET_UNKNOWN_NEGATIVE])
    @pytest.mark.parametrize(("method", "base"), [("c_to_g", 90), ("c_to_g", 91)])
    def test_c_to_g_declines_a_sentinel_offset(
        self,
        mapper: ferro_hgvs.CoordinateMapper,
        method: str,
        base: int,
        offset: int,
    ) -> None:
        """Refused before the arithmetic, in both directions and both arms.

        Reaching the assertion at all is half the test on the negative
        sentinel: a ``PanicException`` subclasses ``BaseException``, so it
        aborts the call rather than raising something ``pytest.raises`` sees.
        """
        with pytest.raises(ferro_hgvs.ProjectionError) as excinfo:
            getattr(mapper, method)(TRANSCRIPT_ID, base, offset)

        assert excinfo.value.code == "E4003"
        assert "unknown intronic offset" in str(excinfo.value)

    @pytest.mark.parametrize("offset", [OFFSET_UNKNOWN_POSITIVE, OFFSET_UNKNOWN_NEGATIVE])
    def test_c_to_n_declines_a_sentinel_offset(
        self, mapper: ferro_hgvs.CoordinateMapper, offset: int
    ) -> None:
        """The value must not come back either — see escape 3."""
        with pytest.raises(ferro_hgvs.ProjectionError) as excinfo:
            mapper.c_to_n(TRANSCRIPT_ID, 90, offset)

        assert excinfo.value.code == "E4003"

    @pytest.mark.parametrize("offset", [OFFSET_UNKNOWN_POSITIVE, OFFSET_UNKNOWN_NEGATIVE])
    def test_n_to_c_declines_a_sentinel_offset(
        self, mapper: ferro_hgvs.CoordinateMapper, offset: int
    ) -> None:
        with pytest.raises(ferro_hgvs.ProjectionError) as excinfo:
            mapper.n_to_c(TRANSCRIPT_ID, 100, offset)

        assert excinfo.value.code == "E4003"

    def test_measured_offsets_still_convert(self, mapper: ferro_hgvs.CoordinateMapper) -> None:
        """A real intronic distance is untouched on both arms.

        The intron is genomic 1101..2000, so ``c.90+5`` is 1105 and ``c.91-5``
        is 1996. A guard that rejected every non-``None`` offset would pass
        every assertion above and fail these.
        """
        assert mapper.c_to_g(TRANSCRIPT_ID, 90, 5) == ("chr1", 1105)
        assert mapper.c_to_g(TRANSCRIPT_ID, 91, -5) == ("chr1", 1996)
        assert mapper.c_to_n(TRANSCRIPT_ID, 90, 5) == (100, 5)
        assert mapper.n_to_c(TRANSCRIPT_ID, 100, 5) == (90, 5, False)


class TestOversizedOffsetIsRefused:
    """The hazard is the offset's magnitude, not its sentinel identity (#1765).

    The guard above keys on the two exact sentinel values, which is the right
    key for a *marker* and the wrong one for the arithmetic downstream. Stepping
    one away from each sentinel reproduced the whole class, measured on this
    fixture on the unfixed build::

        c_to_g(90, 2**63 - 2)   -> ('chr1', 9223372036854776906)
        c_to_g(90, -(2**63) + 1) !! pyo3_runtime.PanicException:
                                    attempt to subtract with overflow
                                    (src/reference/transcript.rs:814)

    Neither value is a sentinel, so no widening of the sentinel predicate could
    reach them. ``Transcript::intronic_to_genomic`` is total and intron-bounded
    instead, so both now decline.
    """

    #: One step in from each sentinel — a plain number, refused on magnitude.
    NEIGHBOURS = [OFFSET_UNKNOWN_POSITIVE - 1, OFFSET_UNKNOWN_NEGATIVE + 1]

    @pytest.mark.parametrize("offset", NEIGHBOURS)
    @pytest.mark.parametrize("base", [90, 91])
    def test_a_sentinel_neighbour_declines_rather_than_panicking(
        self, mapper: ferro_hgvs.CoordinateMapper, base: int, offset: int
    ) -> None:
        """No panic, and no coordinate-shaped number either.

        Reaching the assertion at all is half the test: a ``PanicException``
        subclasses ``BaseException``, so ``pytest.raises(ProjectionError)``
        does not catch it and the pre-fix run errors out here rather than
        failing an assertion.
        """
        with pytest.raises(ferro_hgvs.ProjectionError):
            mapper.c_to_g(TRANSCRIPT_ID, base, offset)

    @pytest.mark.parametrize("offset", [901, -901, 10**6, -(10**6)])
    def test_an_offset_larger_than_its_intron_declines(
        self, mapper: ferro_hgvs.CoordinateMapper, offset: int
    ) -> None:
        """The intron spans genomic 1101..2000, so 900 is the largest honourable
        distance on either arm. An offset past it names no intronic base."""
        with pytest.raises(ferro_hgvs.ProjectionError):
            mapper.c_to_g(TRANSCRIPT_ID, 90 if offset > 0 else 91, offset)

    def test_the_intron_is_still_reachable_end_to_end(
        self, mapper: ferro_hgvs.CoordinateMapper
    ) -> None:
        """The negative control for the bound above: a guard that refused every
        large offset would pass both tests above and fail this one."""
        assert mapper.c_to_g(TRANSCRIPT_ID, 90, 900) == ("chr1", 2000)
        assert mapper.c_to_g(TRANSCRIPT_ID, 91, -900) == ("chr1", 1101)


# ---------------------------------------------------------------------------
# Escape 3: a value handed back that the caller cannot read
# ---------------------------------------------------------------------------


class TestNoUnreadableValueIsReturned:
    """What comes back must mean what a caller would take it to mean."""

    @pytest.mark.parametrize("offset", [OFFSET_UNKNOWN_POSITIVE, OFFSET_UNKNOWN_NEGATIVE])
    def test_c_to_n_never_returns_a_sentinel(
        self, mapper: ferro_hgvs.CoordinateMapper, offset: int
    ) -> None:
        """A sentinel returned as an ``int`` is indistinguishable from a distance.

        Stated as a property of the *return value* rather than of the argument,
        so it still holds if some future conversion path produces a sentinel it
        was not handed.
        """
        try:
            _, returned = mapper.c_to_n(TRANSCRIPT_ID, 90, offset)
        except ferro_hgvs.ProjectionError:
            return  # refused at entry, which is the stronger outcome

        assert returned not in (OFFSET_UNKNOWN_POSITIVE, OFFSET_UNKNOWN_NEGATIVE)

    def test_n_to_c_refuses_the_downstream_flag_it_cannot_honour(
        self, mapper: ferro_hgvs.CoordinateMapper
    ) -> None:
        """``downstream=True`` must not be answered as if it were ``False``.

        ``tx_to_cds`` discards the flag, so before this change the two calls
        returned the same triple and a caller asking about ``n.*195`` was told
        about ``n.195``. Repairing the conversion is separate work; until it
        lands, the honest answer at the boundary is a refusal rather than a
        confidently wrong coordinate.
        """
        with pytest.raises(ferro_hgvs.ProjectionError) as excinfo:
            mapper.n_to_c(TRANSCRIPT_ID, 195, None, True)

        assert excinfo.value.code == "E4003"
        assert "downstream" in str(excinfo.value)

    def test_downstream_false_is_unaffected(self, mapper: ferro_hgvs.CoordinateMapper) -> None:
        """The default path keeps working, positionally and by keyword."""
        assert mapper.n_to_c(TRANSCRIPT_ID, 195) == (5, None, True)
        assert mapper.n_to_c(TRANSCRIPT_ID, 195, None, False) == (5, None, True)
        assert mapper.n_to_c(TRANSCRIPT_ID, 195, downstream=False) == (5, None, True)


# ---------------------------------------------------------------------------
# The docstring defect riding along
# ---------------------------------------------------------------------------


class TestDocstringsStateBothBases:
    """Every numeric entry point names the basis of both sides.

    All five previously left at least one side undocumented: ``c_to_g`` and
    ``c_to_n`` documented neither, and ``c_to_p`` / ``g_to_c`` / ``n_to_c``
    documented the input's basis and not the output's. A caller cannot use an
    integer coordinate whose basis is not stated, and the ambiguity is exactly
    what let these arguments go unvalidated.
    """

    @pytest.mark.parametrize("method", ["c_to_g", "g_to_c", "c_to_p", "c_to_n", "n_to_c"])
    def test_docstring_states_one_based(self, method: str) -> None:
        doc = getattr(ferro_hgvs.CoordinateMapper, method).__doc__
        assert doc is not None, f"{method} has no docstring"
        # Both the Args and the Returns section must say which basis they are
        # on; "1-based" is the vocabulary the rest of this class already uses.
        # Assert the delimiters exist before slicing on them. Without this a
        # missing ``Returns:`` raises a bare ``IndexError: list index out of
        # range`` naming neither the method nor the heading, and — the quieter
        # half — a missing ``Raises:`` does not fail at all: the second split
        # returns a one-element list, so ``returns_section`` silently widens to
        # the whole docstring tail and can be satisfied by text from a section
        # this test is not looking at.
        for heading in ("Returns:", "Raises:"):
            assert heading in doc, f"{method}: docstring has no {heading!r} section"
        args_section = doc.split("Returns:")[0]
        returns_section = doc.split("Returns:")[1].split("Raises:")[0]
        assert "1-based" in args_section, f"{method}: input basis undocumented"
        assert "1-based" in returns_section, f"{method}: output basis undocumented"
