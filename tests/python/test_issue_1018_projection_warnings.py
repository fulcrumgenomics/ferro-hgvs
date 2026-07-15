"""Tests for issue #1018: normalization warnings surfaced on a projection.

`VariantProjector` normalizes its input internally and used to discard the
`NormalizationWarning`s that step emits. Every projection now carries them on
`VariantProjection.warnings` (with a `has_warnings()` helper), mirroring
`Normalizer.normalize_with_warnings` — so the same diagnostics (e.g. an
overlap-conflict or an auto-corrected reference mismatch, and, once it lands on
the normalize path, the reduced-capability degrade) are reachable off the
projector and its batch/fan-out methods, not just the free normalizer.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs


@pytest.fixture
def projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """A projector over transcript NM_TEST.1 mapping chr1 g.1000-1008.

    Mirrors the fixture in test_variant_projection.py: CDS "ATGCGCTAA" on the
    plus strand, g.1000 == c.1. g.1003 is the reference 'C'.
    """
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


# A clean, reference-matching single substitution (no warnings) and a
# coincident cis-allele that normalizes with OVERLAP_CONFLICTING_EDITS.
CLEAN = "NC_000001.11:g.1003C>A"
WARNING = "NC_000001.11:g.[1003C>A;1003C>G]"


class TestCleanProjection:
    def test_warnings_empty_when_input_normalizes_cleanly(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        result = projector.project(CLEAN, transcript="NM_TEST.1")
        assert result.warnings == []
        assert result.has_warnings() is False


class TestWarningSurfaced:
    def test_project_carries_the_normalization_warning(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        result = projector.project(WARNING, transcript="NM_TEST.1")
        assert result.has_warnings() is True
        codes = [w.code for w in result.warnings]
        assert "OVERLAP_CONFLICTING_EDITS" in codes

    def test_warning_objects_are_typed_with_code_and_message(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        result = projector.project(WARNING, transcript="NM_TEST.1")
        assert len(result.warnings) > 0
        w = result.warnings[0]
        assert isinstance(w, ferro_hgvs.NormalizationWarning)
        assert isinstance(w.code, str) and w.code
        assert isinstance(w.message, str) and w.message

    def test_project_all_elements_carry_warnings(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        results = projector.project_all(WARNING)
        assert len(results) >= 1
        assert all(
            "OVERLAP_CONFLICTING_EDITS" in [w.code for w in proj.warnings] for proj in results
        )

    def test_project_many_elements_carry_warnings(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        # project_many returns one inner list per input; each projection in each
        # inner list carries its own warnings.
        batch = projector.project_many([WARNING, CLEAN])
        assert [proj.has_warnings() for proj in batch[0]] == [True]
        assert [proj.has_warnings() for proj in batch[1]] == [False]


class TestNormalizedPathHasNoWarnings:
    """The already-normalized entry points do NOT re-normalize, so they emit no
    normalization warnings — every projection reports an empty list. Pins that
    contract so a caller isn't surprised to find warnings absent there."""

    def test_project_normalized_many_reports_empty_warnings(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        # Pre-normalize the warning-bearing input, then project it as
        # already-normalized: the warning arose during the earlier normalize, so
        # the skip-normalize projection surfaces no warnings of its own.
        variant = ferro_hgvs.parse(WARNING)
        batch = projector.project_normalized_many([variant])
        assert len(batch) == 1
        assert all(not proj.has_warnings() for proj in batch[0])
        assert all(proj.warnings == [] for proj in batch[0])
