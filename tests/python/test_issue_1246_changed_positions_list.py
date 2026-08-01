"""Issue #1246: ``CodonChange.changed_positions`` must be a real ``list[int]``.

The getter returned a Rust ``Vec<u8>``, and pyo3 special-cases that one type to
Python ``bytes`` rather than a list of integers. The published stub declares
``list[int]``, so type-checkers accepted code that then failed at runtime:
``isinstance(x, list)`` was ``False``, ``x == [2]`` was ``False`` (a ``bytes``
never equals a ``list``), and ``json.dumps(x)`` raised ``TypeError``.

Only ``list(x)`` happened to work, because a ``bytes`` iterates as integers —
which is why the mismatch went unnoticed.

These tests assert the *contract the stub advertises*, not the representation,
so they stay meaningful if the underlying Rust integer width ever changes.
"""

import json

import ferro_hgvs


def _a_codon_change() -> ferro_hgvs.CodonChange:
    """A single-nucleotide codon change: Val -> Ala differs at position 2.

    Val is GTx and Ala is GCx across the whole codon table, so *every* codon
    pair this returns changes exactly the second base. Taking ``changes[0]``
    therefore does not depend on the order the backtranslator emits them.
    """
    changes = ferro_hgvs.Backtranslator.standard().backtranslate_substitution("Val", "Ala")
    assert changes, "expected at least one codon change for Val -> Ala"
    return changes[0]


def test_changed_positions_is_a_list() -> None:
    positions = _a_codon_change().changed_positions
    assert isinstance(positions, list)
    assert not isinstance(positions, bytes)


def test_changed_positions_holds_ints() -> None:
    positions = _a_codon_change().changed_positions
    assert positions, "expected at least one changed position"
    assert all(isinstance(position, int) for position in positions)
    # Every Val codon is GTx and every Ala codon GCx, so the substitution
    # changes the second base and only the second base. Asserting the exact
    # value, not just the 1..=3 range, is what catches a regression that keeps
    # the `list[int]` type but reports the wrong index.
    assert positions == [2]


def test_changed_positions_compares_equal_to_a_list() -> None:
    # The failure that motivated the issue: a bytes never equals a list, so
    # stub-conformant comparisons silently returned False.
    positions = _a_codon_change().changed_positions
    assert positions == list(positions)


def test_changed_positions_is_json_serializable() -> None:
    # `json.dumps` raised "Object of type bytes is not JSON serializable".
    positions = _a_codon_change().changed_positions
    assert json.loads(json.dumps(positions)) == positions


def test_changed_positions_is_mutable_like_a_list() -> None:
    # `bytes` has no `.append`; a `list[int]` does, and the stub promises one.
    positions = _a_codon_change().changed_positions
    before = len(positions)
    positions.append(1)
    assert len(positions) == before + 1


def test_changed_positions_matches_num_changes() -> None:
    change = _a_codon_change()
    assert len(change.changed_positions) == change.num_changes()
