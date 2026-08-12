"""Deriving an HGVS description from a reference/alternate sequence pair.

The scenario this surface exists for: a pipeline aligns reads, post-processes a
BAM, and wants one canonical, spec-conformant description per variant determined
by the reference and observed bases alone. Today that pipeline has to choose an
HGVS spelling itself and hand it to ``normalize`` — and that choice leaks into
the output.

``from_sequences`` takes the bases instead. It delivers the two normalization
rules that are always achievable, **conformant** and **deterministic**, and
deliberately not the two that need the reference, **recommended form** and
**confluent**. ``Normalizer.to_sequences`` is the inverse, so a caller holding
descriptions rather than bases needs no new plumbing to reach it.

The tests below are the Python half of that contract. The properties themselves
— confluence over the cis corpus's 5 636 **genomic** classes (its other 5 636 are
``c.`` against ``NM_TEST.1``, which the ``g.``-only gate refuses, so they never
enter the comparison; the exclusion is structural and is asserted as such in
``tests/it/from_sequences_corpus.rs``), the nine externally-reported #1419/#1420/
#1421 pairs, denotation preservation over 23 121 derivations — are judged in
Rust, where the corpora live. What is checked here is what only the binding can
get wrong: the keyword-only signature, the 1-based position, which exception
comes out of which failure, and that the two directions actually reach the
partitioner rather than being accepted and dropped.
"""

import json
import tempfile
from pathlib import Path

import pytest

import ferro_hgvs

# 1-based:      1234567890
#
# A 4-base `A` homopolymer run at 12-15, and `GGATTACAGG` repeated at 1-10 and
# 22-31 — so a derivation has something to shuffle and the two directions are
# distinguishable.
#
# The run's coordinates read 11-14 here until 2026-08-12, which is where it sits
# inside `test_the_direction_reaches_the_partitioner`'s *window* (that pair
# starts at 10, not at 1) rather than in this sequence. #1694's Rust counterpart,
# `tests/it/from_sequences_reanchor.rs`, says 12-15 and is right.
SEQUENCE = "GGATTACAGGCAAAAGCCTGAGGATTACAGGCATTAGCCT"


@pytest.fixture(scope="module")
def normalizer():
    """A provider whose bases are known, so expected strings can be written by
    hand rather than read back out of the code under test."""
    reference = {
        "transcripts": [
            {
                "id": "TEMPLATE",
                "gene_symbol": "T",
                "strand": "+",
                "sequence": SEQUENCE,
                "exons": [{"number": 1, "start": 1, "end": len(SEQUENCE)}],
            }
        ]
    }
    path = Path(tempfile.mkdtemp()) / "reference.json"
    path.write_text(json.dumps(reference))
    return ferro_hgvs.Normalizer(reference_json=str(path))


# ---------------------------------------------------------------------------
# The free function: a pure function of its arguments
# ---------------------------------------------------------------------------


def test_derives_a_description_from_bases_alone():
    """No reference, no provider, no Normalizer — just the four arguments."""
    variant = ferro_hgvs.from_sequences("NC_000001.11", 1000, "AGCG", "AG")
    assert str(variant) == "NC_000001.11:g.1002_1003del"


def test_the_position_is_one_based():
    """The single most likely thing for a caller to get wrong, and the reason
    ``SequencePair`` exists rather than ``AppliedVariant`` being handed straight
    across: ``AppliedVariant.start`` is 0-based and this is not.

    Asserted against a coordinate computed by hand: the window starts at 1000,
    the last two bases of ``AGCG`` are dropped, so the deletion is 1002_1003.
    """
    assert str(ferro_hgvs.from_sequences("NC_1", 1, "AGCG", "AG")) == "NC_1:g.3_4del"
    assert str(ferro_hgvs.from_sequences("NC_1", 2, "AGCG", "AG")) == "NC_1:g.4_5del"


def test_the_same_bases_give_the_same_description():
    """Rule 4, at the binding. Not a deep property — the Rust side proves it over
    the corpora — but a binding that leaked a hash seed or an iteration order
    into the output would fail here and nowhere else in this file."""
    calls = [ferro_hgvs.from_sequences("NC_1", 5, "GGCATTAGC", "GGCTTAGC") for _ in range(8)]
    assert len({str(v) for v in calls}) == 1


def test_an_unchanged_pair_is_the_identity_description():
    assert str(ferro_hgvs.from_sequences("NC_1", 5, "GGCA", "GGCA")) == "NC_1:g.="


# ---------------------------------------------------------------------------
# The signature
# ---------------------------------------------------------------------------


def test_the_options_are_keyword_only():
    """``max_grid_cells`` and ``direction`` are keyword-only on purpose: they are
    tuning knobs, and a positional fifth argument would read as part of the
    variant rather than as configuration."""
    with pytest.raises(TypeError):
        ferro_hgvs.from_sequences("NC_1", 1, "AGCG", "AG", 1024)  # type: ignore[misc]


def test_the_direction_reaches_the_partitioner():
    """Both directions must be accepted *and act*. A binding that parsed the
    keyword and then dropped it would pass a "does not raise" test, so this pins
    a pair where the two answers differ — a single-base deletion inside the
    ``AAAA`` run, which 3' places at 14 and 5' at 11. Those are coordinates in
    *this pair's* frame (it starts at 10 and is written out below, not sliced
    from ``SEQUENCE``); the run sits at 12-15 of ``SEQUENCE`` itself."""
    three = ferro_hgvs.from_sequences("NC_1", 10, "CAAAAG", "CAAAG", direction="3prime")
    five = ferro_hgvs.from_sequences("NC_1", 10, "CAAAAG", "CAAAG", direction="5prime")
    assert str(three) == "NC_1:g.14del"
    assert str(five) == "NC_1:g.11del"


def test_an_unrecognized_direction_is_a_value_error():
    with pytest.raises(ValueError, match="unrecognized direction"):
        ferro_hgvs.from_sequences("NC_1", 1, "AGCG", "AG", direction="sideways")


# ---------------------------------------------------------------------------
# Refusals: which exception, and does it say what to do
# ---------------------------------------------------------------------------


def test_a_zero_grid_budget_is_a_value_error():
    """Checked at the binding, before any derivation — so it is a ``ValueError``
    rather than a ``NormalizationError``, and the docs say which is which."""
    with pytest.raises(ValueError, match="max_grid_cells must be positive"):
        ferro_hgvs.from_sequences("NC_1", 1, "AGCG", "AG", max_grid_cells=0)


def test_a_grid_over_budget_refuses_and_names_the_knob():
    """The operator's stated policy is to refuse above the window rather than
    degrade to a weaker rule, so a refusal has to be actionable: it names the
    knob and its per-cell cost."""
    with pytest.raises(ferro_hgvs.NormalizationError) as excinfo:
        ferro_hgvs.from_sequences("NC_1", 1, "ACGTACGTAC", "TGCATGCATG", max_grid_cells=4)
    assert "max_grid_cells" in str(excinfo.value)


def test_a_transcript_accession_is_refused():
    """``NM_TEST.1:g.9_10del`` is well-formed and denotes nothing — a ``g.``
    description on a transcript is not one the recommendations admit
    (``checklist.md:20``). Found by running the corpus, not by a unit test."""
    with pytest.raises(ferro_hgvs.NormalizationError):
        ferro_hgvs.from_sequences("NM_000088.3", 10, "AGCG", "AG")


@pytest.mark.parametrize(
    ("position", "reference", "alternate"),
    [
        (0, "AGCG", "AG"),  # 1-based; 0 names no base
        (1, "", "AG"),  # no interval to describe
        (1, "AGXG", "AG"),  # standards.md:39 — X is alignment-only
    ],
)
def test_unusable_input_refuses(position, reference, alternate):
    with pytest.raises(ferro_hgvs.NormalizationError):
        ferro_hgvs.from_sequences("NC_1", position, reference, alternate)


def test_the_two_refusal_classes_are_the_established_split():
    """Argument-shape errors are a plain ``ValueError``; everything the
    derivation itself refuses is a ``NormalizationError``, which is a
    ``FerroError``.

    That split is the binding's existing convention, not a choice made here —
    ``Normalizer(direction=...)`` has always raised a plain ``ValueError`` for an
    unrecognized direction. Pinned because the obvious wrong guess is that
    ``FerroError`` catches everything from one call, and it does not: a caller
    wrapping ``from_sequences`` needs ``(ValueError, ferro_hgvs.FerroError)``, or
    simply ``Exception``."""
    binding_checks = (
        lambda: ferro_hgvs.from_sequences("NC_1", 1, "AGCG", "AG", max_grid_cells=0),
        lambda: ferro_hgvs.from_sequences("NC_1", 1, "AGCG", "AG", direction="sideways"),
    )
    for call in binding_checks:
        with pytest.raises(ValueError) as excinfo:
            call()
        assert not isinstance(excinfo.value, ferro_hgvs.FerroError)

    with pytest.raises(ferro_hgvs.FerroError):
        ferro_hgvs.from_sequences("NC_1", 0, "AGCG", "AG")


# ---------------------------------------------------------------------------
# from_sequences_detailed
# ---------------------------------------------------------------------------


def test_a_derivation_reports_whether_it_reached_a_window_edge():
    """The one caveat a window-local derivation owes its caller: on an edge, the
    3' placement may lie outside the bases supplied, so a wider window could move
    the answer. Reported rather than silently absorbed.

    The pair is the same deletion twice, once with flank and once without."""
    interior = ferro_hgvs.from_sequences_detailed("NC_1", 10, "CAAAAGCC", "CAAAGCC")
    assert interior.placement_bounded_by_window is False

    flush = ferro_hgvs.from_sequences_detailed("NC_1", 10, "AAAA", "AAA")
    assert flush.placement_bounded_by_window is True
    assert str(flush.variant) == "NC_1:g.13del"


def test_the_detailed_variant_matches_the_plain_call():
    plain = ferro_hgvs.from_sequences("NC_1", 10, "CAAAAGCC", "CAAAGCC")
    detailed = ferro_hgvs.from_sequences_detailed("NC_1", 10, "CAAAAGCC", "CAAAGCC")
    assert str(detailed.variant) == str(plain)


def test_derived_description_has_a_useful_repr():
    rendered = repr(ferro_hgvs.from_sequences_detailed("NC_1", 10, "CAAAAGCC", "CAAAGCC"))
    assert rendered.startswith("DerivedDescription(")
    assert "placement_bounded_by_window" in rendered


# ---------------------------------------------------------------------------
# to_sequences: the inverse
# ---------------------------------------------------------------------------


def test_to_sequences_returns_a_one_based_window(normalizer):
    variant = normalizer.parse("TEMPLATE:n.13del")
    pair = normalizer.to_sequences(variant, pad=4)

    assert pair.accession == "TEMPLATE"
    # The window is the member's span padded on both sides, and `position` is
    # 1-based, so `SEQUENCE[position - 1]` is its first base.
    assert SEQUENCE[pair.position - 1 : pair.position - 1 + len(pair.reference)] == pair.reference
    assert len(pair.alternate) == len(pair.reference) - 1
    assert pair.window_is_final is True


def test_to_sequences_pads_both_sides(normalizer):
    """The 5' half of the pad is load-bearing and was missing at first: ``dup``
    typing reads the reference bases immediately 5' of an insertion point
    (``duplication.md:18``), so a member flush with the window's 5' edge comes
    back as an ``ins`` instead of a ``dup``. Ten corpus classes reached both
    spellings before this was fixed."""
    variant = normalizer.parse("TEMPLATE:n.13del")
    pair = normalizer.to_sequences(variant, pad=4)
    assert pair.position < 13
    assert pair.position - 1 + len(pair.reference) > 13


def test_to_sequences_clamps_the_pad_to_the_sequence(normalizer):
    """Near a sequence start there is less flank than asked for, so the window
    is clamped rather than the request refused.

    ``window_is_final`` stays True through that, and the distinction is the
    point: a window clipped by the *sequence end* has nothing further to read,
    so the roll is settled and a 3' placement against that edge is trustworthy.
    Only a provider stopping short of both makes it False. This test asserted
    the pad-shaped reading first and failed, which is what renamed the field
    from ``full_pad_served``."""
    variant = normalizer.parse("TEMPLATE:n.2del")
    pair = normalizer.to_sequences(variant, pad=100)
    assert pair.position == 1
    assert pair.reference == SEQUENCE
    assert pair.window_is_final is True


def test_to_sequences_refuses_what_apply_to_reference_refuses(normalizer):
    """``to_sequences`` promises ``ProjectionError`` and had no test at all, so
    the promise was unchecked in the one direction it could be wrong — the
    binding maps it through a *literal* class name, and nothing but a test on
    that name notices when the name and the docstring drift apart."""
    with pytest.raises(ferro_hgvs.ProjectionError):
        normalizer.to_sequences(normalizer.parse("NOT_A_TEMPLATE:n.13del"), pad=4)


def test_sequence_pair_has_a_useful_repr(normalizer):
    rendered = repr(normalizer.to_sequences(normalizer.parse("TEMPLATE:n.13del"), pad=4))
    assert rendered.startswith("SequencePair(")
    assert "position=" in rendered


def test_sequence_pairs_compare_by_value(normalizer):
    """The Rust type derives ``PartialEq``; a ``#[pyclass]`` does not forward it,
    so without an explicit ``__eq__`` two pairs holding the same window compared
    unequal by identity — the kind of default an ``assert a != b`` passes over in
    silence."""
    variant = normalizer.parse("TEMPLATE:n.13del")
    one = normalizer.to_sequences(variant, pad=4)
    same = normalizer.to_sequences(variant, pad=4)
    wider = normalizer.to_sequences(variant, pad=6)

    assert one == same
    assert one is not same
    assert one != wider
    assert len({one, same, wider}) == 2


# ---------------------------------------------------------------------------
# The round trip, and the Normalizer method
# ---------------------------------------------------------------------------


def test_a_description_survives_the_round_trip(normalizer):
    """description -> sequences -> description. The two functions check each
    other, which is exactly how the Rust corpus tests turn every committed HGVS
    corpus into a ``from_sequences`` corpus with no new fixtures."""
    for spelling in ["TEMPLATE:n.13del", "TEMPLATE:n.3_5del", "TEMPLATE:n.13dup"]:
        pair = normalizer.to_sequences(normalizer.parse(spelling), pad=16)
        derived = normalizer.from_sequences(
            pair.accession, pair.position, pair.reference, pair.alternate
        )
        # Compared on the bases denoted, not on the string: a derivation may
        # legitimately narrow the span or choose another legal spelling, and
        # `canonical_spdi` is independent of both.
        assert normalizer.canonical_spdi(derived) == normalizer.canonical_spdi(
            normalizer.parse(spelling)
        )


def test_the_method_refuses_an_interval_past_the_sequence_end(normalizer):
    """The free function cannot range-check ``position`` — that needs the
    reference, which would make the provider a hidden input and cost exactly the
    determinism it exists to provide. The method holds a provider, so it does.

    Asserted on the exact class, not on ``FerroError``. The binding promises one
    per method, and a ``FerroError`` assertion passes for all four of them — so
    it cannot catch a docstring that promises a class the code never raises,
    which is exactly what this method's did (``ReferenceDataError`` for both of
    the refusals it adds over the free function)."""
    with pytest.raises(ferro_hgvs.NormalizationError):
        normalizer.from_sequences("TEMPLATE", len(SEQUENCE), "AGCG", "AG")


def test_the_method_refuses_an_unknown_accession(normalizer):
    """The other refusal the method adds over the free function, and the other
    one its docstring promised as ``ReferenceDataError``. That class is this
    binding's for failures *loading* reference data; a per-call lookup is the
    method's own class, as ``normalize_variant`` has always been."""
    with pytest.raises(ferro_hgvs.NormalizationError):
        normalizer.from_sequences("NOT_A_TEMPLATE", 10, "CAAAAG", "CAAAG")


def test_the_method_uses_the_normalizers_own_direction():
    """One direction per normalizer, set at construction — the contract every
    other method here follows, rather than a second place to set it."""
    reference = {
        "transcripts": [
            {
                "id": "TEMPLATE",
                "gene_symbol": "T",
                "strand": "+",
                "sequence": SEQUENCE,
                "exons": [{"number": 1, "start": 1, "end": len(SEQUENCE)}],
            }
        ]
    }
    path = Path(tempfile.mkdtemp()) / "reference.json"
    path.write_text(json.dumps(reference))

    args = ("TEMPLATE", 10, "CAAAAG", "CAAAG")
    three = ferro_hgvs.Normalizer(reference_json=str(path), direction="3prime")
    five = ferro_hgvs.Normalizer(reference_json=str(path), direction="5prime")
    assert str(three.from_sequences(*args)) != str(five.from_sequences(*args))


def test_post_normalizing_is_available_and_is_a_no_op_here(normalizer):
    """``normalize=False`` is the ordinary posture. On measured material the flag
    is a no-op **on these rows** — a 6,000-shape sweep found it moves 8.6% of
    derived descriptions overall, so this is a property of the rows below and
    never of the surface. On them a derived description is a fixed point of ``normalize``
    on the genomic axis — so its value is the diagnostics ``normalize`` raises
    and inheriting its improvements as they land. Pinned as a measurement on
    these rows, never as a law."""
    pair = normalizer.to_sequences(normalizer.parse("TEMPLATE:n.13del"), pad=16)
    plain = normalizer.from_sequences(pair.accession, pair.position, pair.reference, pair.alternate)
    normalized = normalizer.from_sequences(
        pair.accession, pair.position, pair.reference, pair.alternate, normalize=True
    )
    assert str(normalized) == str(plain)


def test_the_method_rejects_a_zero_grid_budget(normalizer):
    with pytest.raises(ValueError, match="max_grid_cells must be positive"):
        normalizer.from_sequences("TEMPLATE", 10, "CAAAAG", "CAAAG", max_grid_cells=0)


# ---------------------------------------------------------------------------
# The exported surface
# ---------------------------------------------------------------------------


def test_the_surface_is_exported():
    """Both classes and both functions must be reachable as ``ferro_hgvs.X`` and
    named in ``__all__`` — a binding registered in the pymodule but missing from
    ``__all__`` is importable and invisible to ``from ferro_hgvs import *``."""
    for name in (
        "from_sequences",
        "from_sequences_detailed",
        "SequencePair",
        "DerivedDescription",
    ):
        assert hasattr(ferro_hgvs, name), f"{name} is not exported"
        assert name in ferro_hgvs.__all__, f"{name} is missing from __all__"


# ---------------------------------------------------------------------------
# Re-anchoring: bounding a derivation to a window it must not leave
# ---------------------------------------------------------------------------


def test_a_pair_can_be_built_from_bases_alone():
    """`SequencePair` was previously only ever *returned*, which put both
    re-anchoring entry points out of reach of the caller they are for: one who
    has bases out of a BAM and no description at all."""
    pair = ferro_hgvs.SequencePair("T", 10, "GCAAAAG", "GCAAAG")
    assert (pair.position, pair.end) == (10, 16)
    assert pair.window_is_final is False


def test_a_pair_that_constructs_is_a_pair_that_derives():
    """The constructor shares `from_sequences`'s check rather than copying it,
    so the refusal lands on the argument that caused it."""
    for position, reference, alternate in [(0, "AGCG", "AG"), (1, "", "AG"), (1, "AGXG", "AG")]:
        with pytest.raises(ferro_hgvs.NormalizationError):
            ferro_hgvs.SequencePair("T", position, reference, alternate)


def test_the_constructor_admits_the_iupac_ambiguity_codes():
    """The alphabet is ``Base::from_char``'s — the IUPAC-IUBMB set of
    ``general.md:48`` — not ``ACGTN``, which is what it was before #1693 widened
    it and what both Python docstrings still claimed. Real ClinVar rows carry the
    ambiguity codes (``NM_000518.4:c.[20A>T;249G>Y]`` among them), so a narrow
    alphabet refuses submitted data.

    ``U`` stays refused, and that exclusion is this surface's own rather than
    ``from_char``'s: this surface emits ``g.``/``m.`` descriptions, which are
    DNA."""
    assert ferro_hgvs.SequencePair("T", 1, "AGYG", "AG").reference == "AGYG"

    for bad in ("AGUG", "AGXG", "AG-G"):
        with pytest.raises(ferro_hgvs.NormalizationError):
            ferro_hgvs.SequencePair("T", 1, bad, "AG")


def test_a_bound_stops_the_roll_at_the_bound():
    """The headline behaviour. Unbounded the deletion rolls 3' to 15; bounded at
    14 it stays at 14."""
    pair = ferro_hgvs.SequencePair("NC_1", 10, "GCAAAAG", "GCAAAG")
    assert str(
        ferro_hgvs.from_sequences("NC_1", pair.position, pair.reference, pair.alternate)
    ) == ("NC_1:g.15del")

    bounded = pair.trim_to(end=14)
    assert (bounded.position, bounded.end) == (10, 14)
    derived = ferro_hgvs.from_sequences_detailed(
        "NC_1", bounded.position, bounded.reference, bounded.alternate
    )
    assert str(derived.variant) == "NC_1:g.14del"
    assert derived.placement_bounded_by_window is True


def test_trim_to_refuses_rather_than_clamping():
    pair = ferro_hgvs.SequencePair("NC_1", 10, "GCAAAAG", "GCAAAAT")
    # Cuts the substituted base at 16.
    with pytest.raises(ferro_hgvs.NormalizationError, match="differ"):
        pair.trim_to(end=15)
    # Widening is not something a reference-free object can do; the message says
    # which method can.
    with pytest.raises(ferro_hgvs.NormalizationError, match="reanchor"):
        pair.trim_to(start=5)


def test_trim_to_treats_case_as_agreement():
    """A soft-masked reference against an upper-case alternate is an ordinary
    pileup. The comparison folds — ``Base.from_char`` does, and the constructor's
    own validation says case is not a reason to refuse — so a byte comparison
    here refused a legal trim claiming the two sequences "first differ" at a
    coordinate where they do not.

    ``trim_to`` fetches nothing, so it cannot manufacture a mixed-case pair and
    leaves the bases as passed. ``Normalizer.reanchor``, which does fetch, folds
    the whole window instead."""
    masked = ferro_hgvs.SequencePair("NC_1", 10, "gcaaaag", "GCAAAG")
    bounded = masked.trim_to(end=14)
    assert (bounded.position, bounded.end) == (10, 14)
    assert bounded.reference == "gcaaa"
    assert str(bounded.derive().variant) == "NC_1:g.14del"

    # A real disagreement is still refused, so the fold is not a blanket accept.
    with pytest.raises(ferro_hgvs.NormalizationError, match="differ"):
        ferro_hgvs.SequencePair("NC_1", 10, "gcaaaag", "GCAAAAT").trim_to(end=15)


def test_reanchor_pads_and_trims(normalizer):
    """One call that widens one edge and narrows the other."""
    pair = normalizer.to_sequences(normalizer.parse("TEMPLATE:n.13del"), pad=2)
    widened = normalizer.reanchor(pair, start=1, end=20)
    assert (widened.position, widened.end) == (1, 20)
    assert widened.reference == SEQUENCE[:20]
    assert len(widened.alternate) == len(widened.reference) - 1


def test_reanchor_keeps_window_is_final_when_the_three_prime_edge_never_moves(normalizer):
    """``window_is_final`` was recomputed as ``end == length`` unconditionally, so
    an *identity* re-anchor downgraded a settled window to False — and no test in
    either language fed a ``to_sequences`` pair to ``reanchor`` and then read the
    flag, which is the only way a True ever reaches it.

    The rule is ``trim_to``'s, stated over coordinates because this method can
    widen as well as narrow: reaching the sequence's own end settles the edge,
    and so does leaving an already-settled edge alone."""
    settled = normalizer.to_sequences(normalizer.parse("TEMPLATE:n.13del"), pad=4)
    assert settled.window_is_final is True
    assert settled.end < len(SEQUENCE), "the fixture must stop short of the sequence end"

    assert normalizer.reanchor(settled).window_is_final is True
    assert normalizer.reanchor(settled, start=1).window_is_final is True
    assert normalizer.reanchor(settled, end=settled.end + 1).window_is_final is False
    assert normalizer.reanchor(settled, end=len(SEQUENCE)).window_is_final is True

    # A caller-built pair carries False by construction, and moving nothing does
    # not settle it.
    raw = ferro_hgvs.SequencePair("TEMPLATE", 10, "GCAAAAG", "GCAAAG")
    assert normalizer.reanchor(raw).window_is_final is False
    assert normalizer.reanchor(raw, end=len(SEQUENCE)).window_is_final is True


def test_reanchor_refuses_to_leave_the_sequence(normalizer):
    """The operator's call: a caller who asked for bases that do not exist has a
    bug upstream, and a silently clamped window would hide it.

    On the exact class, not on ``FerroError``: the docstring promised
    ``ReferenceDataError`` for precisely these two bounds and the binding maps
    every error here through ``NormalizationError``, so a ``FerroError``
    assertion is what let the promise sit green."""
    pair = normalizer.to_sequences(normalizer.parse("TEMPLATE:n.13del"), pad=2)
    whole = normalizer.reanchor(pair, start=1, end=len(SEQUENCE))
    assert (whole.position, whole.end) == (1, len(SEQUENCE))
    assert whole.reference == SEQUENCE
    with pytest.raises(ferro_hgvs.NormalizationError):
        normalizer.reanchor(pair, start=1, end=len(SEQUENCE) + 1)
    with pytest.raises(ferro_hgvs.NormalizationError):
        normalizer.reanchor(pair, start=0)


def test_reanchor_moves_a_windows_edges_and_does_not_relocate_it(normalizer):
    """The constraint the README taught backwards. ``reanchor`` pads or trims
    each edge, in any combination — but the window asked for must overlap the
    pair's own, because the changed bases exist only in the pair. A disjoint
    request is refused, not fetched, and the message says so rather than naming
    ``trim_to`` and two coordinates the caller never passed."""
    pair = normalizer.to_sequences(normalizer.parse("TEMPLATE:n.13del"), pad=2)

    with pytest.raises(ferro_hgvs.NormalizationError, match="does not overlap") as excinfo:
        normalizer.reanchor(pair, start=30, end=len(SEQUENCE))
    assert "trim_to" not in str(excinfo.value)

    # A partial overlap is legal — the boundary is overlap, not containment.
    # [11, 15] narrowed to 13 and widened to 30; the deletion still rolls to the
    # 3' end of the ``A`` run at 12-15, which the new window still contains.
    partial = normalizer.reanchor(pair, start=13, end=30)
    assert (partial.position, partial.end) == (13, 30)
    assert str(partial.derive().variant) == "TEMPLATE:g.15del"


def test_reanchor_upper_cases_the_window_it_returns(normalizer):
    """``reanchor`` reads flank from the provider and keeps the caller's bases in
    the middle, so a soft-masked window would otherwise come back mixed-case — a
    pair no caller wrote. It folds the whole window, exactly as ``to_sequences``
    does. ``trim_to``, which fetches nothing, does not."""
    masked = ferro_hgvs.SequencePair("TEMPLATE", 10, "gcaaaag", "GCAAAG")
    widened = normalizer.reanchor(masked, start=5, end=25)
    assert widened.reference == SEQUENCE[4:25]
    assert widened.alternate == widened.alternate.upper()


def test_reads_anchored_to_one_region_agree():
    """The property the feature exists for, over windows that disagree without
    it."""
    seq = "GGATTACAGGCAAAAGCCTGAGGATTACAGGCATTAGCCT"

    def read(lo, hi):
        reference = seq[lo - 1 : hi]
        at = reference.index("A", 12 - lo)
        return ferro_hgvs.SequencePair("NC_1", lo, reference, reference[:at] + reference[at + 1 :])

    reads = [read(9, 16), read(10, 17), read(11, 18)]
    anchored = {
        str(ferro_hgvs.from_sequences("NC_1", p.position, p.reference, p.alternate))
        for p in (r.trim_to(start=11, end=14) for r in reads)
    }
    assert anchored == {"NC_1:g.14del"}


def test_a_pair_derives_from_itself():
    """`derive` over the four values the pair already carries.

    The point is that a pair returned by `trim_to` carries its **own**
    position, so a caller never pairs a pre-trim position with post-trim bases
    — the mistake the README example used to teach.
    """
    pair = ferro_hgvs.SequencePair("NC_1", 10, "GCAAAAG", "GCAAAG")
    assert str(pair.derive().variant) == "NC_1:g.15del"

    bounded = pair.trim_to(end=14)
    derived = bounded.derive()
    assert str(derived.variant) == "NC_1:g.14del"
    assert derived.placement_bounded_by_window is True

    # Same knobs as the free function, and the same refusals.
    assert str(pair.derive(direction="5prime").variant) == "NC_1:g.12del"
    with pytest.raises(ValueError, match="max_grid_cells must be positive"):
        pair.derive(max_grid_cells=0)
    with pytest.raises(ValueError, match="unrecognized direction"):
        pair.derive(direction="sideways")
