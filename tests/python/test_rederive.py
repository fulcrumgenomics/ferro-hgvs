"""``Normalizer.rederive`` and ``DerivedDescription``'s per-side flags.

``rederive`` expresses a variant as a padded reference/alternate window
(``to_sequences``), derives a description from those bases alone
(``from_sequences``), and — while a member still rests on a window edge that can
still move — doubles the pad and retries. Two spellings of one variant therefore
reach one description, decided by the observed bases.

The properties themselves are judged in Rust, over
``tests/it/rederive.rs`` and the cis corpus. What is checked here is
what only the binding can get wrong: the keyword-only signature, that
``recommended_form`` is threaded through rather than accepted and dropped, which
exception class comes out of which failure, and that the two new
``DerivedDescription`` getters are wired to the fields they name rather than
both to ``placement_bounded_by_window``.
"""

import json
import tempfile
from pathlib import Path

import pytest

import ferro_hgvs

# 1-based:                    1..=10        11..=310        311..=322
#
# A 300-base `A` run flanked by unique sequence on both sides. The run is far
# longer than the 128-base start pad, so a single fetch cannot contain the whole
# roll and the widening loop has to run.
LONG_RUN = "GCTAGCTAGC" + "A" * 300 + "GCTAGCTAGCTA"

# 1-based:            1..=6      7..=14     15..=20
#
# `ACGTACGT` at 7-14, so the same variant can be spelled as a duplication of
# 7_14 or as an insertion of those eight bases after 14.
DUP_CONTIG = "GCTAGC" + "ACGTACGT" + "TCGATC"


def _normalizer(sequence: str) -> ferro_hgvs.Normalizer:
    """A genome-capable provider over one synthetic contig named ``NC_TEST.1``.

    The name is genomic on purpose: ``from_sequences``, which
    ``rederive`` calls, emits ``g.`` and refuses a transcript or
    protein accession, so a non-genomic contig name would be rejected before any
    derivation.

    There is no ``direction=`` here because the bindings no longer take one —
    3' is the only direction ferro shifts. The 5'-shuffle half of the per-side
    split is therefore only reachable from Rust, and is pinned there by
    ``tests/it/rederive.rs``.
    """
    payload = {
        "transcripts": [],
        "proteins": {},
        "genomic_sequences": {"NC_TEST.1": sequence},
    }
    path = Path(tempfile.mkdtemp()) / "reference.json"
    path.write_text(json.dumps(payload))
    return ferro_hgvs.Normalizer(reference_json=str(path))


# ---------------------------------------------------------------------------
# Confluence
# ---------------------------------------------------------------------------


def test_three_spellings_of_one_deletion_converge():
    """The headline contract, at the binding: one variant, three spellings, one
    description — and reaching it requires the loop to widen past its start pad."""
    normalizer = _normalizer(LONG_RUN)
    results = {
        normalizer.rederive(spelling)
        for spelling in ("NC_TEST.1:g.11del", "NC_TEST.1:g.160del", "NC_TEST.1:g.310del")
    }
    assert results == {"NC_TEST.1:g.310del"}


def test_a_duplication_and_its_insertion_spelling_converge():
    """The shape a narrow first window gets silently wrong.

    ``dup`` typing reads the reference bases immediately 5' of the insertion
    point, so a window that does not reach them derives an ``ins`` instead — and
    that mis-typed ``ins`` rests on neither window edge, so neither per-side flag
    fires and the loop would return it without widening.
    """
    normalizer = _normalizer(DUP_CONTIG)
    assert normalizer.rederive("NC_TEST.1:g.7_14dup") == "NC_TEST.1:g.7_14dup"
    assert normalizer.rederive("NC_TEST.1:g.14_15insACGTACGT") == "NC_TEST.1:g.7_14dup"


def test_a_run_reaching_the_contig_end_settles_by_the_sequence():
    """A stall, not an overrun. The run reaches the contig's last base, so the 3'
    edge the placement rests on is the sequence's own — the loop must return
    rather than widen past it or decline it as unbounded."""
    normalizer = _normalizer("GCTAGCTAGC" + "A" * 40)
    assert normalizer.rederive("NC_TEST.1:g.15del") == "NC_TEST.1:g.50del"


# ---------------------------------------------------------------------------
# The signature
# ---------------------------------------------------------------------------


def test_max_grid_cells_and_recommended_form_are_keyword_only():
    normalizer = _normalizer(LONG_RUN)
    with pytest.raises(TypeError):
        normalizer.rederive("NC_TEST.1:g.11del", 4096)  # type: ignore[misc]
    with pytest.raises(TypeError):
        normalizer.rederive("NC_TEST.1:g.11del", normalize_flag=True)  # type: ignore[call-arg]


def test_recommended_form_true_reaches_the_normalizer():
    """``recommended_form`` must be threaded through, not accepted and dropped.

    A lone reference `A` immediately upstream of an insertion of three more
    `A`s is a case where the two flag values give genuinely different strings:
    repeat-notation collapse is a `normalize` (rule 2) pass that
    `from_sequences`'s own derivation never performs (see
    `from_sequences.rs`'s `retype_inversions` docs, which name this exact
    shape as still un-typed there). So `recommended_form=False` must derive
    the literal insertion and `recommended_form=True` must collapse it to
    repeat notation — if the keyword were silently dropped, the second
    assertion would see the first string instead and fail.
    """
    core = "GCTAGC" + "A" + "TCGATC"  # a lone 'A' at position 7, 'T' at 8
    normalizer = _normalizer(core)
    assert normalizer.rederive("NC_TEST.1:g.7_8insAAA") == "NC_TEST.1:g.7_8insAAA"
    assert (
        normalizer.rederive("NC_TEST.1:g.7_8insAAA", recommended_form=True) == "NC_TEST.1:g.7A[4]"
    )


# ---------------------------------------------------------------------------
# Which exception comes out of which failure
# ---------------------------------------------------------------------------


def test_an_unparseable_description_raises_parse_error():
    normalizer = _normalizer(LONG_RUN)
    with pytest.raises(ferro_hgvs.ParseError):
        normalizer.rederive("not an hgvs description")


def test_a_zero_grid_budget_raises_value_error():
    """Checked at the binding and outside the ``FerroError`` hierarchy, as for
    ``from_sequences``."""
    normalizer = _normalizer(LONG_RUN)
    with pytest.raises(ValueError):
        normalizer.rederive("NC_TEST.1:g.11del", max_grid_cells=0)


def test_an_unknown_accession_raises_normalization_error():
    normalizer = _normalizer(LONG_RUN)
    with pytest.raises(ferro_hgvs.NormalizationError):
        normalizer.rederive("NC_ABSENT.1:g.11del")


# ---------------------------------------------------------------------------
# The two new DerivedDescription getters
# ---------------------------------------------------------------------------


def test_the_per_side_flags_are_wired_to_their_own_fields():
    """Both getters returning ``placement_bounded_by_window`` would pass any test
    that only checks the OR, so each side is asserted alone — on two pairs whose
    answers differ, and without reaching for a shuffle direction the bindings no
    longer expose.

    A substitution flush against one end of the window rests on that edge and no
    other, which separates the two flags under the only direction there is.
    """
    # Changed base at the window's last position: 3' edge only.
    derived = ferro_hgvs.from_sequences_detailed("NC_TEST.1", 5, "ACGT", "ACGA")
    assert derived.bounded_at_end is True
    assert derived.bounded_at_start is False
    assert derived.placement_bounded_by_window is True

    # Changed base at the window's first position: 5' edge only.
    derived = ferro_hgvs.from_sequences_detailed("NC_TEST.1", 5, "ACGT", "TCGT")
    assert derived.bounded_at_start is True
    assert derived.bounded_at_end is False
    assert derived.placement_bounded_by_window is True


def test_an_interior_change_flags_neither_side():
    derived = ferro_hgvs.from_sequences_detailed("NC_TEST.1", 5, "TAAAAG", "TAAAG")
    assert derived.bounded_at_start is False
    assert derived.bounded_at_end is False
    assert derived.placement_bounded_by_window is False


def test_the_repr_names_all_three_flags():
    derived = ferro_hgvs.from_sequences_detailed("NC_TEST.1", 5, "ACGT", "ACGA")
    text = repr(derived)
    assert "placement_bounded_by_window=True" in text
    assert "bounded_at_start=False" in text
    assert "bounded_at_end=True" in text
