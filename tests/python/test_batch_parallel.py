"""Tests for the parallel batch API (``BatchProcessor.parse`` / ``parse_and_normalize`` /
``parse_and_rederive``).

The core property is invariance: parallel results must equal serial results
(same order, same per-item success/error classification) for any worker count.
The ``parse`` / ``parse_and_normalize`` tests use the mock provider
(``BatchProcessor()`` with no manifest), so they need no reference data. The
``parse_and_rederive`` tests need a genome-capable provider — ``rederive`` reads
the bases a variant denotes, which the mock test data does not carry — so they
build a small single-contig ``reference_json`` under ``tmp_path``.
"""

import json
from pathlib import Path

import ferro_hgvs


def _variants(n: int) -> list[str]:
    """`n` distinct variants with two deliberate parse errors at known positions."""
    assert n >= 3, "need n >= 3 for two distinct, in-bounds error indices"
    vs = [f"NM_000088.3:c.{i}A>G" for i in range(1, n + 1)]
    vs[1] = "definitely not a variant"
    vs[n - 1] = "also::not:valid"
    return vs


def _key(result) -> tuple:
    """Order-sensitive fingerprint of a BatchResult: counts, the ordered success
    descriptions, and the ordered ``(index, message)`` errors — so a run that
    misplaces or mis-tags a failure fails the comparison even when the success
    list and the counts still match."""
    return (
        result.total(),
        result.success_count(),
        result.error_count(),
        [str(v) for v in result.successes()],
        list(result.errors()),
    )


def test_parse_parallel_matches_serial() -> None:
    bp = ferro_hgvs.BatchProcessor()
    variants = _variants(5000)

    serial = _key(bp.parse(variants, workers=1))
    assert serial[1] > 0, "expected some successful parses"
    assert serial[2] == 2, "expected exactly the two injected errors"

    for workers in (0, 2, 4, 8):
        assert _key(bp.parse(variants, workers=workers)) == serial, (
            f"parse(workers={workers}) differs from serial"
        )


def test_parse_and_normalize_parallel_matches_serial() -> None:
    bp = ferro_hgvs.BatchProcessor()
    variants = _variants(5000)

    serial = _key(bp.parse_and_normalize(variants, workers=1))
    assert serial[1] > 0, "expected some successful parses"
    assert serial[2] == 2, "expected exactly the two injected errors"

    for workers in (0, 2, 4, 8):
        assert _key(bp.parse_and_normalize(variants, workers=workers)) == serial, (
            f"parse_and_normalize(workers={workers}) differs from serial"
        )


def test_parallel_is_deterministic() -> None:
    bp = ferro_hgvs.BatchProcessor()
    variants = _variants(4000)
    a = _key(bp.parse_and_normalize(variants, workers=8))
    b = _key(bp.parse_and_normalize(variants, workers=8))
    assert a == b


def test_default_workers_is_parallel_and_correct() -> None:
    # The default (workers=0, all cores) must produce the same result as serial.
    bp = ferro_hgvs.BatchProcessor()
    variants = _variants(2000)
    assert _key(bp.parse(variants)) == _key(bp.parse(variants, workers=1))


def _genome_bp(tmp_path: Path) -> "ferro_hgvs.BatchProcessor":
    """A genome-capable BatchProcessor over one synthetic contig.

    ``rederive`` reads the bases a variant denotes, so it needs genomic sequence
    data — the default ``BatchProcessor()`` test data has none, and rederive would
    error on every input there.
    """
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
    """`n` in-bounds genomic deletions that rederive successfully, with two parse
    errors at known indices."""
    assert n >= 3, "need n >= 3 for two distinct, in-bounds error indices"
    vs = [f"NC_000001.11:g.{100 + i}del" for i in range(n)]
    vs[1] = "definitely not a variant"
    vs[n - 1] = "also::not:valid"
    return vs


def test_parse_and_rederive_parallel_matches_serial(tmp_path: Path) -> None:
    bp = _genome_bp(tmp_path)
    variants = _genome_variants(2000)

    serial = _key(bp.parse_and_rederive(variants, workers=1))
    assert serial[1] > 0, "expected successful rederivations, not an all-error corpus"
    assert serial[2] == 2, "expected exactly the two injected parse errors"

    for workers in (0, 2, 4, 8):
        assert _key(bp.parse_and_rederive(variants, workers=workers)) == serial, (
            f"parse_and_rederive(workers={workers}) differs from serial"
        )


def _insertion_repeat_bp(tmp_path: Path) -> "ferro_hgvs.BatchProcessor":
    """A genome-capable BatchProcessor over a contig with a lone 'A' at position 7.

    On this contig ``g.7_8insAAA`` derives, unnormalized, as the literal
    ``g.7_8insAAA``; ``recommended_form=True`` routes it through ``normalize``,
    which collapses it to the repeat ``g.7A[4]``. The two settings therefore
    produce different descriptions — which is what lets a test prove the flag
    reaches the parallel path rather than only that the call succeeds. Mirrors
    ``tests/it/rederive.rs::recommended_form_true_routes_through_normalize``.
    """
    ref = {
        "transcripts": [],
        "proteins": {},
        "genomic_sequences": {"NC_TEST.1": "GCTAGCATCGATC"},
    }
    path = tmp_path / "insertion_repeat_ref.json"
    path.write_text(json.dumps(ref))
    return ferro_hgvs.BatchProcessor(reference_json=str(path))


def test_parse_and_rederive_recommended_form_runs_parallel(tmp_path: Path) -> None:
    # recommended_form must reach the parallel rederive path AND change the
    # output: the literal insertion `g.7_8insAAA` collapses to the repeat
    # `g.7A[4]` only when the flag is honoured. Asserting the two settings
    # differ — and stay invariant across worker counts — fails if the binding
    # silently drops the flag.
    bp = _insertion_repeat_bp(tmp_path)
    variants = ["NC_TEST.1:g.7_8insAAA"]
    for workers in (0, 1, 2, 4, 8):
        plain = bp.parse_and_rederive(variants, recommended_form=False, workers=workers)
        recommended = bp.parse_and_rederive(variants, recommended_form=True, workers=workers)
        assert [str(v) for v in plain.successes()] == ["NC_TEST.1:g.7_8insAAA"], (
            f"recommended_form=False (workers={workers}) must keep the literal insertion"
        )
        assert [str(v) for v in recommended.successes()] == ["NC_TEST.1:g.7A[4]"], (
            f"recommended_form=True (workers={workers}) must collapse to repeat notation"
        )


def test_parse_and_rederive_is_deterministic(tmp_path: Path) -> None:
    bp = _genome_bp(tmp_path)
    variants = _genome_variants(1000)
    a = _key(bp.parse_and_rederive(variants, workers=8))
    b = _key(bp.parse_and_rederive(variants, workers=8))
    assert a == b


def test_parse_and_rederive_rejects_zero_grid() -> None:
    import pytest

    # No reference needed: the 0-grid ValueError is raised before any work.
    bp = ferro_hgvs.BatchProcessor()
    with pytest.raises(ValueError):
        bp.parse_and_rederive(["NM_000088.3:c.1A>G"], max_grid_cells=0)
