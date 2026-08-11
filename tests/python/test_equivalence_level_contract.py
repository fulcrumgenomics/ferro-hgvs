"""``EquivalenceLevel``'s two new predicates, exercised through PyO3.

``is_decided`` and ``is_at_least`` are pinned in Rust by
``tests/it/equivalence_cross_axis_rung.rs``. A Rust test cannot cover what is
specific to the binding, and there are three such things:

* the two new discriminants (``CrossAxisSequenceMatch``, ``Indeterminate``) are
  bound as attributes at all — the enum is hand-mapped in ``src/python.rs``, so
  a member can exist in Rust and be missing from Python;
* ``is_at_least`` takes an ``EquivalenceLevel`` *argument*, which is the only
  place in this enum's surface where a level is converted Python -> Rust;
* the predicates return real ``bool``, not a truthy proxy.

The order itself is asserted here as well as in Rust deliberately. A binding
that silently mapped every level onto one discriminant would satisfy "the
attribute exists" while answering every question wrong.
"""

import ferro_hgvs

LEVEL = ferro_hgvs.EquivalenceLevel


def test_the_new_discriminants_are_bound() -> None:
    """Both levels this rung added are reachable from Python."""
    for name in ("CrossAxisSequenceMatch", "Indeterminate"):
        member = getattr(LEVEL, name)
        assert isinstance(member, LEVEL), f"{name} is not an EquivalenceLevel member"


def test_indeterminate_is_the_only_undecided_level() -> None:
    """``is_decided`` separates "no" from "cannot tell", and only there."""
    assert LEVEL.Indeterminate.is_decided() is False
    for name in (
        "Identical",
        "NormalizedMatch",
        "AccessionVersionDifference",
        "NotEquivalent",
        "SequenceMatch",
        "CrossAxisSequenceMatch",
    ):
        assert getattr(LEVEL, name).is_decided() is True, f"{name} should be decided"


def test_indeterminate_is_not_a_positive_verdict() -> None:
    """Undecidable answers ``False`` to ``is_equivalent`` as well."""
    assert LEVEL.Indeterminate.is_equivalent() is False


def test_the_denotational_order_holds_across_the_binding() -> None:
    """``is_at_least`` converts its argument and ranks the denotational rungs."""
    assert LEVEL.Identical.is_at_least(LEVEL.CrossAxisSequenceMatch) is True
    assert LEVEL.Identical.is_at_least(LEVEL.SequenceMatch) is True
    assert LEVEL.CrossAxisSequenceMatch.is_at_least(LEVEL.SequenceMatch) is True
    # Reflexive on the order, and antisymmetric off it.
    assert LEVEL.SequenceMatch.is_at_least(LEVEL.SequenceMatch) is True
    assert LEVEL.SequenceMatch.is_at_least(LEVEL.CrossAxisSequenceMatch) is False
    assert LEVEL.CrossAxisSequenceMatch.is_at_least(LEVEL.Identical) is False


def test_normalized_match_can_neither_satisfy_nor_be_a_gate() -> None:
    """The circular rung is off the order in both directions.

    A caller writing a confluence gate in Python gets the same protection the
    Rust API has: ``NormalizedMatch`` is defined by the normalizer the gate is
    about, so it can neither clear a floor nor serve as one.
    """
    assert LEVEL.NormalizedMatch.is_at_least(LEVEL.SequenceMatch) is False
    assert LEVEL.NormalizedMatch.is_at_least(LEVEL.NormalizedMatch) is False
    assert LEVEL.Identical.is_at_least(LEVEL.NormalizedMatch) is False
    assert LEVEL.CrossAxisSequenceMatch.is_at_least(LEVEL.NormalizedMatch) is False


def test_the_accession_version_rung_is_off_the_order_but_decided() -> None:
    """Orthogonal, not weaker: apply-equality is undefined across two versions."""
    assert LEVEL.AccessionVersionDifference.is_decided() is True
    assert LEVEL.AccessionVersionDifference.is_at_least(LEVEL.SequenceMatch) is False
    assert LEVEL.Identical.is_at_least(LEVEL.AccessionVersionDifference) is False
