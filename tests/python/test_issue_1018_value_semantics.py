"""Tests for issue #1018 (P3 polish): parse idempotency + value equality.

Two things:

* A round-trip/idempotency property mirroring the Rust idempotency suite:
  ``parse(str(parse(x))) == parse(x)`` for a corpus of variants across axes.
* ``__eq__`` / ``__hash__`` on the value types that lacked them
  (``VariantProjection``, ``ProteinEffect``), so projection/effect results are
  usable as ``dict`` keys and ``set`` members — matching ``HgvsVariant`` and
  ``SpdiVariant`` which already compared by value.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs

CORPUS = [
    "NM_000088.3:c.459A>G",
    "NC_000001.11:g.1000A>T",
    "NM_000088.3:c.100_102del",
    "NM_000088.3:c.100dup",
    "NP_000079.2:p.Gly100Ala",
    "NM_000088.3:c.100+5G>A",
    "NR_003051.3:n.100A>G",
]


class TestParseIdempotency:
    """``parse(str(parse(x))) == parse(x)`` — re-parsing a rendered variant is a
    fixed point, by value and by hash."""

    @pytest.mark.parametrize("variant", CORPUS)
    def test_roundtrip_is_a_fixed_point(self, variant: str) -> None:
        once = ferro_hgvs.parse(variant)
        twice = ferro_hgvs.parse(str(once))
        assert once == twice
        assert hash(once) == hash(twice)

    def test_hgvs_variant_is_set_usable(self) -> None:
        variants = [ferro_hgvs.parse(v) for v in CORPUS]
        # Each parsed once and again; the set collapses the duplicates.
        again = [ferro_hgvs.parse(str(v)) for v in variants]
        assert len(set(variants + again)) == len(CORPUS)


class TestCrossTypeEqualityDoesNotRaise:
    """`__eq__` against `None` or an unrelated type must return `False`, never
    raise. This rests on PyO3 returning `NotImplemented` when the argument
    cannot be extracted to `&Self`; pin it so a future switch to an explicit
    `__richcmp__` that forgets the mismatch case is caught."""

    def test_hgvs_variant_cross_type(self) -> None:
        v = ferro_hgvs.parse("NM_000088.3:c.459A>G")
        assert (v == None) is False  # noqa: E711 — testing == against None, not `is`
        assert v != "NM_000088.3:c.459A>G"
        assert (v == 5) is False

    def test_protein_effect_cross_type(self) -> None:
        effect = ferro_hgvs.EffectPredictor().classify_indel(3, 0)
        assert (effect == None) is False  # noqa: E711
        assert effect != "frameshift"
        assert (effect == 5) is False

    def test_projection_cross_type(self, projector: ferro_hgvs.VariantProjector) -> None:
        proj = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        assert (proj == None) is False  # noqa: E711
        assert proj != "NM_TEST.1"
        assert (proj == 5) is False
        # A projection and an effect are different types → not equal, no raise.
        effect = ferro_hgvs.EffectPredictor().classify_indel(3, 0)
        assert proj != effect


class TestProteinEffectValueSemantics:
    """``ProteinEffect`` compares and hashes by value."""

    def test_equal_effects_are_equal_and_hash_equal(self) -> None:
        predictor = ferro_hgvs.EffectPredictor()
        a = predictor.classify_indel(3, 0)
        b = predictor.classify_indel(3, 0)
        assert a == b
        assert hash(a) == hash(b)

    def test_distinct_effects_differ(self) -> None:
        predictor = ferro_hgvs.EffectPredictor()
        inframe = predictor.classify_indel(3, 0)
        frameshift = predictor.classify_indel(1, 0)
        assert inframe != frameshift

    def test_effect_is_set_usable(self) -> None:
        predictor = ferro_hgvs.EffectPredictor()
        effects = [
            predictor.classify_indel(3, 0),
            predictor.classify_indel(3, 0),
            predictor.classify_indel(1, 0),
        ]
        assert len(set(effects)) == 2


@pytest.fixture
def projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """Projector over NM_TEST.1 (chr1 g.1000-1008, CDS 'ATGCGCTAA', + strand)."""
    fixture = {
        "transcripts": [
            {
                "id": "NM_TEST.1",
                "gene_symbol": "TESTGENE",
                "strand": "+",
                "sequence": "ATGCGCTAA",
                "cds_start": 1,
                "cds_end": 9,
                "exons": [
                    {
                        "number": 1,
                        "start": 1,
                        "end": 9,
                        "genomic_start": 1000,
                        "genomic_end": 1008,
                    }
                ],
                "chromosome": "chr1",
                "genomic_start": 1000,
                "genomic_end": 1008,
            }
        ],
        "genomic_sequences": {"chr1": "N" * 999 + "ATGCGCTAA" + "N" * 100},
    }
    path = tmp_path / "transcripts.json"
    path.write_text(json.dumps(fixture))
    return ferro_hgvs.VariantProjector(reference_json=str(path))


class TestVariantProjectionValueSemantics:
    """``VariantProjection`` compares and hashes by value."""

    def test_equal_projections_are_equal_and_hash_equal(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        a = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        b = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        assert a == b
        assert hash(a) == hash(b)

    def test_distinct_projections_differ(self, projector: ferro_hgvs.VariantProjector) -> None:
        a = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        c = projector.project("NC_000001.11:g.1004G>T", transcript="NM_TEST.1")
        assert a != c

    def test_projection_is_set_usable(self, projector: ferro_hgvs.VariantProjector) -> None:
        a = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        b = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        c = projector.project("NC_000001.11:g.1004G>T", transcript="NM_TEST.1")
        assert len({a, b, c}) == 2

    def test_warning_carrying_projection_is_hashable(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        # The eq/hash key is the inner value's Debug string, so a projection that
        # carries a `NormalizationWarning` (which is itself unhashable) is still
        # hashable — the Debug-key design's whole point. A coincident cis-allele
        # emits OVERLAP_CONFLICTING_EDITS, giving a warning-bearing projection.
        results = projector.project_all("NC_000001.11:g.[1003C>A;1003C>G]")
        warned = [p for p in results if p.has_warnings()]
        assert warned, "expected a warning-carrying projection from the cis-allele"
        proj = warned[0]
        # Must not raise TypeError, and must be usable in a set.
        assert isinstance(hash(proj), int)
        assert len({proj}) == 1
        assert proj == proj
