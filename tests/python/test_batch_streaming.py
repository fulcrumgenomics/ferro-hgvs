"""Tests for the streaming batch API (#975).

``BatchProcessor.parse`` / ``parse_and_normalize`` take a fully materialized
``list`` and return a fully materialized ``BatchResult``, so both the input and
the output are resident at once *by contract* — no internal chunking can bound
peak memory behind that signature. ``parse_streaming`` /
``parse_and_normalize_streaming`` consume an iterable and yield, so peak memory
is a function of the chunk size rather than of the input's length.

Two properties are tested, and the second is the one an equality check cannot
see:

1. **Equivalence** — streaming yields exactly what the list API returns, in the
   same order, with the same per-item success/error classification.
2. **Laziness** — the input iterable is pulled from only as results are demanded.
   A streaming API that collected its input first would satisfy (1) completely
   while bounding nothing.

These use the mock provider (``BatchProcessor()`` with no manifest), so they need
no reference data.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs


def _variants(n: int) -> list[str]:
    """`n` inputs with deliberate parse errors at known positions."""
    assert n >= 3
    vs = [f"NM_000088.3:c.{i}A>G" for i in range(1, n + 1)]
    vs[1] = "definitely not a variant"
    vs[n - 1] = "also::not:valid"
    return vs


def _from_list(result) -> list[tuple[bool, str]]:
    """Order-sensitive fingerprint of a ``BatchResult``, per item."""
    out: list[tuple[bool, str]] = []
    errors = dict(result.errors())
    successes = iter(result.successes())
    for i in range(result.total()):
        if i in errors:
            out.append((False, ""))
        else:
            out.append((True, str(next(successes))))
    return out


def _from_stream(stream) -> list[tuple[bool, str]]:
    return [(item.ok, str(item.variant) if item.ok else "") for item in stream]


@pytest.fixture
def processor():
    return ferro_hgvs.BatchProcessor()


@pytest.mark.parametrize("workers", [0, 1, 2])
def test_streaming_parse_matches_the_list_api(processor, workers):
    """Streaming must be item-for-item identical to the list API."""
    variants = _variants(50)
    expected = _from_list(processor.parse(variants, workers=workers))
    streamed = _from_stream(processor.parse_streaming(variants, workers=workers))
    assert streamed == expected


@pytest.mark.parametrize("workers", [0, 1, 2])
def test_streaming_normalize_matches_the_list_api(processor, workers):
    variants = _variants(50)
    expected = _from_list(processor.parse_and_normalize(variants, workers=workers))
    streamed = _from_stream(processor.parse_and_normalize_streaming(variants, workers=workers))
    assert streamed == expected


def _genome_bp(tmp_path: Path) -> "ferro_hgvs.BatchProcessor":
    """A genome-capable BatchProcessor: ``rederive`` reads the bases a variant
    denotes, which the default ``BatchProcessor()`` test data does not carry."""
    seq = ["A"] * 5000
    for i in range(1000, 4000, 7):
        seq[i] = "G"
    ref = {
        "transcripts": [],
        "proteins": {},
        "genomic_sequences": {"NC_000001.11": "".join(seq)},
    }
    path = tmp_path / "genome_ref.json"
    path.write_text(json.dumps(ref))
    return ferro_hgvs.BatchProcessor(reference_json=str(path))


def _genome_variants(n: int) -> list[str]:
    """`n` in-bounds genomic deletions that rederive, with two parse errors."""
    assert n >= 3
    vs = [f"NC_000001.11:g.{100 + i}del" for i in range(n)]
    vs[1] = "definitely not a variant"
    vs[n - 1] = "also::not:valid"
    return vs


@pytest.mark.parametrize("workers", [0, 1, 2])
def test_streaming_rederive_matches_the_list_api(tmp_path: Path, workers: int) -> None:
    bp = _genome_bp(tmp_path)
    variants = _genome_variants(50)
    expected = _from_list(bp.parse_and_rederive(variants, workers=workers))
    # A genome-capable provider means rederive actually succeeds, so this
    # compares streaming vs list over real (Ok) results, not only errors.
    assert any(ok for ok, _ in expected), "expected successful rederivations"
    streamed = _from_stream(bp.parse_and_rederive_streaming(variants, workers=workers))
    assert streamed == expected


@pytest.mark.parametrize("workers", [0, 1, 2])
def test_streaming_rederive_recommended_form_matches_the_list_api(
    tmp_path: Path, workers: int
) -> None:
    bp = _genome_bp(tmp_path)
    variants = _genome_variants(40)
    expected = _from_list(bp.parse_and_rederive(variants, recommended_form=True, workers=workers))
    # As in the sibling above: prove the comparison runs over real (Ok) results,
    # not a vacuously-passing all-error corpus.
    assert any(ok for ok, _ in expected), "expected successful rederivations"
    streamed = _from_stream(
        bp.parse_and_rederive_streaming(variants, recommended_form=True, workers=workers)
    )
    assert streamed == expected


def test_streaming_rederive_rejects_zero_grid(processor: "ferro_hgvs.BatchProcessor") -> None:
    # The 0-grid ValueError is raised at the call, before any item is streamed.
    with pytest.raises(ValueError):
        processor.parse_and_rederive_streaming(["NM_000088.3:c.1A>G"], max_grid_cells=0)


def test_failed_items_carry_their_input(processor):
    """A failure must be reportable without the caller keeping the input list —
    otherwise a streaming consumer would have to retain what it streamed."""
    items = list(processor.parse_streaming(["NM_000088.3:c.1A>G", "nope"]))
    assert items[0].ok and items[0].input is None and items[0].error is None
    assert not items[1].ok
    assert items[1].input == "nope"
    assert items[1].error and items[1].variant is None


def test_the_input_is_consumed_lazily(processor):
    """The property that makes peak memory bounded.

    Nothing may be read from the iterable before the first ``next()``. This is
    what an output-equality test cannot detect.
    """
    pulled: list[str] = []

    def source():
        for v in _variants(10):
            pulled.append(v)
            yield v

    stream = processor.parse_streaming(source())
    assert pulled == [], "constructing the stream must read nothing"
    first = next(stream)
    assert first.ok
    assert pulled, "the first result must have pulled from the source"
    rest = list(stream)
    assert len(rest) == 9
    assert len(pulled) == 10, "every input read exactly once"


def test_the_first_result_costs_a_bounded_prefix_not_the_whole_stream(processor):
    """The half of laziness ``test_the_input_is_consumed_lazily`` cannot see.

    That test uses ten inputs, which is far below the streaming chunk size, so
    the first ``next()`` legitimately consumes all of them. Every one of its
    assertions therefore also holds for an implementation that drains the entire
    iterable on the first ``next()`` — it pins "reads nothing at construction",
    not "reads a bounded prefix", and only the second bounds peak memory.

    Asserting the *contract* (bounded) rather than the chunk constant, so this
    keeps testing the right thing if the chunk size is retuned — which it has
    been, and the Python and Rust paths deliberately use different values.
    """
    total = 30_000  # comfortably several chunks on any plausible chunk size
    pulled: list[str] = []
    exhausted = False

    def source():
        nonlocal exhausted
        for v in _variants(total):
            pulled.append(v)
            yield v
        exhausted = True

    stream = processor.parse_streaming(source())
    assert pulled == [], "constructing the stream must read nothing"

    next(stream)
    after_first = len(pulled)
    assert 0 < after_first < total, (
        f"the first result must cost a bounded prefix, not the whole stream; "
        f"pulled {after_first} of {total}"
    )
    assert not exhausted, "the source generator must not have run to completion"

    # Draining the rest reads every remaining input exactly once, so bounding the
    # prefix has not cost any input.
    consumed = 1 + sum(1 for _ in stream)
    assert consumed == total
    assert len(pulled) == total, "every input read exactly once"
    assert exhausted


def test_the_stream_is_its_own_iterator(processor):
    stream = processor.parse_streaming(["NM_000088.3:c.1A>G"])
    assert iter(stream) is stream
    assert len(list(stream)) == 1
    # Exhausted streams stay exhausted rather than re-polling a spent iterator.
    assert list(stream) == []


def test_accepts_any_iterable(processor):
    """A list, a generator and a plain iterator must all work — the point is to
    accept whatever the caller already has, including a file object."""
    variants = _variants(5)
    for source in (variants, iter(variants), (v for v in variants)):
        assert len(list(processor.parse_streaming(source))) == 5


def test_empty_input_yields_nothing(processor):
    assert list(processor.parse_streaming([])) == []
    assert list(processor.parse_and_normalize_streaming([])) == []
    assert list(processor.parse_and_rederive_streaming([])) == []


def test_non_string_element_raises(processor):
    """Surfaced, not skipped: a silently dropped item would desynchronize the
    caller's own bookkeeping against the stream."""
    with pytest.raises(TypeError):
        list(processor.parse_streaming(["NM_000088.3:c.1A>G", 42]))


def test_non_iterable_argument_raises_at_the_call(processor):
    """Raised where the mistake is, rather than on the first ``next()``."""
    with pytest.raises(TypeError):
        processor.parse_streaming(42)


def test_the_streaming_types_are_exported():
    """``BatchStream`` and ``BatchItem`` are part of the public surface the stubs
    declare, so they belong in ``__all__`` alongside the other batch classes.

    ``from .ferro_hgvs import *`` makes them reachable as attributes either way,
    but ``__all__`` is what governs ``from ferro_hgvs import *`` and what tools
    read as the package's public API — and ``__init__.pyi`` already promises both
    names.
    """
    assert "BatchStream" in ferro_hgvs.__all__
    assert "BatchItem" in ferro_hgvs.__all__
    # The siblings they were added beside, so the set stays coherent.
    for sibling in ("BatchProgress", "BatchResult", "BatchProcessor"):
        assert sibling in ferro_hgvs.__all__
