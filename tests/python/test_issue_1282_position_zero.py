"""Issue #1282: a payload that lands before the first base, from Python.

A 5'-shuffled cis allele can denote "the reference with bases added at the very
front". Before this fix the derived member reached ``hgvs_pos_to_index`` with
``start == 0`` and the subtraction underflowed. The Rust suite
(``tests/it/issue_1282_position_zero.rs``) covers the fix at the library level;
this file exists because ``Normalizer.normalize`` is the surface where the
underflow did the most damage.

A Rust panic crosses PyO3 as ``pyo3_runtime.PanicException``, which subclasses
``BaseException`` rather than ``Exception`` — so it slips straight through the
``except Exception`` a caller would reasonably wrap ``normalize`` in. And the
wheel is built ``--release`` (``ci.yml``'s Python Wheel Test job, and every
published artifact), where ``[profile.release]`` sets no ``overflow-checks`` —
so there the subtraction wrapped silently to an index near ``usize::MAX``
rather than panicking at all. Both doors are shut by the same rewrite, and
both are only observable from here.

``tests/python/test_coordinates.py`` covers the other half — the standalone
``hgvs_pos_to_index(0)`` conversion — but that is the helper, not the path that
reached it.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs

CONTIG = "TEMPLATE"

# Leading ``T`` run, so a 5' shuffle drives the payload to the contig start.
# Identical to the Rust fixture's ``SEQ``, so the two suites pin one behaviour.
SEQ = "TTTTTTTTTAATATATTTTA"

# The reported reproducer and two siblings, all of which underflowed. Each
# denotes the reference with an ``A`` added before base 1: the substitution
# rewrites base 1 to ``A`` and the insertion lengthens the ``T`` run. The
# insertion's own start position is irrelevant to the outcome, which is why all
# three collapse to one answer.
REPRODUCERS = [
    "TEMPLATE:g.[3_4insT;1T>A]",
    "TEMPLATE:g.[2_3insT;1T>A]",
    "TEMPLATE:g.[5_6insT;1T>A]",
]

# "Insert ``A`` immediately 5' of base 1" spelled as the boundary delins that
# every other bound in ``src/normalize/mod.rs`` already uses, because
# ``g.0_1ins`` is not a position any HGVS axis has.
BOUNDARY_DELINS = "TEMPLATE:g.1delinsAT"


def _normalizer(tmp_path: Path) -> "ferro_hgvs.Normalizer":
    """Build a 5'-shuffling normalizer over the synthetic reference."""
    reference = {"transcripts": [], "genomic_sequences": {CONTIG: SEQ}}
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(reference))
    return ferro_hgvs.Normalizer(reference_json=str(path), direction="5prime")


@pytest.mark.parametrize("variant", REPRODUCERS)
def test_payload_before_the_first_base_normalizes(tmp_path: Path, variant: str) -> None:
    """The reported crash: a description comes back instead of a panic.

    Reaching the assertion at all is half the test — a ``PanicException`` would
    abort the call rather than return.
    """
    assert _normalizer(tmp_path).normalize(variant) == BOUNDARY_DELINS


def test_the_boundary_delins_is_a_fixed_point(tmp_path: Path) -> None:
    """Re-normalizing the output must not move it, or the clamp is a one-pass patch."""
    normalizer = _normalizer(tmp_path)

    assert normalizer.normalize(BOUNDARY_DELINS) == BOUNDARY_DELINS


def test_an_interior_payload_keeps_its_insertion_spelling(tmp_path: Path) -> None:
    """The clamp fires only at interbase 0.

    Without this, a clamp that fired too eagerly would rewrite valid insertions
    into ``delins`` and still pass the assertions above.
    """
    normalizer = _normalizer(tmp_path)

    assert normalizer.normalize("TEMPLATE:g.3_4insT") == "TEMPLATE:g.1dup"
    assert normalizer.normalize("TEMPLATE:g.1T>A") == "TEMPLATE:g.1T>A"
