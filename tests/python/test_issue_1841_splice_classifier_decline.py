"""Issue #1841: ``EffectPredictor.classify_splice_variant`` returns ``None`` for an unknown offset.

``c.<base>-?`` / ``c.<base>+?`` are legal HGVS: they say the position is
intronic and decline to say how far from the boundary it sits. The parser
carries that as a sentinel — ``-2**63`` and ``2**63 - 1`` — and the classifier's
whole input is a distance, so there is nothing for it to classify.

Before #1806 this was a live defect: the sentinel's ``.abs()`` wrapped negative
in release, passed the ``<= 2`` rung, and the position was reported as a
canonical splice site with **HIGH** impact. #1806 fixed the answer to
``intron_variant``/MODIFIER, which is true, and recorded that the honest answer
was out of reach because the return type could not decline. #1841 took that
break: the Rust function returns ``Option<ProteinEffect>``, so the Python method
returns ``ProteinEffect | None``.

**This is a breaking change to the Python surface.** A caller that reads
``.consequences`` off the result unconditionally now raises ``AttributeError``
on these two inputs, where it previously read ``intron_variant``.

There were no Python tests for this method at all before this file, in either
direction, which is why the sentinel neighbours below are pinned as well as the
sentinels themselves: a decline keyed on "a very large magnitude" rather than on
the sentinel would satisfy the first two assertions and be wrong.
"""

import pytest

import ferro_hgvs

#: How the parser spells an unknown intronic offset, per
#: ``ferro_hgvs::hgvs::parser::position``.
OFFSET_UNKNOWN_NEGATIVE = -(2**63)
OFFSET_UNKNOWN_POSITIVE = 2**63 - 1

SENTINELS = [OFFSET_UNKNOWN_NEGATIVE, OFFSET_UNKNOWN_POSITIVE]

#: Ordinary magnitudes that merely sit next to the sentinels. They are not
#: sentinels, so they must still classify.
SENTINEL_NEIGHBOURS = [OFFSET_UNKNOWN_NEGATIVE + 1, OFFSET_UNKNOWN_POSITIVE - 1]


@pytest.fixture
def predictor() -> ferro_hgvs.EffectPredictor:
    """A predictor; it holds no reference and no state."""
    return ferro_hgvs.EffectPredictor()


@pytest.mark.parametrize("offset", SENTINELS)
def test_an_unknown_offset_returns_none(predictor: ferro_hgvs.EffectPredictor, offset: int) -> None:
    """The break itself: no answer is expressed as ``None``, not as a weak answer."""
    assert predictor.classify_splice_variant(offset) is None


@pytest.mark.parametrize("offset", SENTINEL_NEIGHBOURS)
def test_a_sentinel_neighbour_still_classifies(
    predictor: ferro_hgvs.EffectPredictor, offset: int
) -> None:
    """The discriminating case: the decline is keyed on the sentinel, not on magnitude."""
    effect = predictor.classify_splice_variant(offset)
    assert effect is not None
    assert effect.consequences == [ferro_hgvs.Consequence.IntronVariant]
    assert effect.impact == ferro_hgvs.Impact.Modifier


@pytest.mark.parametrize(
    ("offset", "consequence", "impact"),
    [
        (1, ferro_hgvs.Consequence.SpliceDonorVariant, ferro_hgvs.Impact.High),
        (-2, ferro_hgvs.Consequence.SpliceAcceptorVariant, ferro_hgvs.Impact.High),
        (5, ferro_hgvs.Consequence.SpliceRegionVariant, ferro_hgvs.Impact.Low),
        (50, ferro_hgvs.Consequence.IntronVariant, ferro_hgvs.Impact.Modifier),
    ],
)
def test_a_measured_offset_still_classifies(
    predictor: ferro_hgvs.EffectPredictor,
    offset: int,
    consequence: "ferro_hgvs.Consequence",
    impact: "ferro_hgvs.Impact",
) -> None:
    """Every ordinary rung is unchanged — only the sentinel class moved."""
    effect = predictor.classify_splice_variant(offset)
    assert effect is not None
    assert effect.consequences == [consequence]
    assert effect.impact == impact
