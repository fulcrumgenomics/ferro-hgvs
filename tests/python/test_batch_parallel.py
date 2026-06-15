"""Tests for the parallel batch API (``BatchProcessor.parse`` / ``parse_and_normalize``).

The core property is invariance: parallel results must equal serial results
(same order, same per-item success/error classification) for any worker count.
These use the mock provider (``BatchProcessor()`` with no manifest), so they need
no reference data.
"""

import ferro_hgvs


def _variants(n: int) -> list[str]:
    """`n` distinct variants with two deliberate parse errors at known positions."""
    assert n >= 3, "need n >= 3 for two distinct, in-bounds error indices"
    vs = [f"NM_000088.3:c.{i}A>G" for i in range(1, n + 1)]
    vs[1] = "definitely not a variant"
    vs[n - 1] = "also::not:valid"
    return vs


def _key(result) -> tuple:
    """Order-sensitive fingerprint of a BatchResult."""
    return (
        result.total(),
        result.success_count(),
        result.error_count(),
        [str(v) for v in result.successes()],
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
