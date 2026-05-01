"""Tests for Normalizer/EquivalenceChecker/BatchProcessor/CoordinateMapper loader paths."""

from __future__ import annotations

import json
from pathlib import Path

import ferro_hgvs

FIXTURES = Path(__file__).parent.parent / "fixtures"
MANIFEST_TINY = FIXTURES / "python" / "manifest_tiny" / "manifest.json"


class TestObjectFormJson:
    """Object-form reference_json: {transcripts, proteins, genomic_sequences}."""

    def test_normalizer_loads_object_form(self, tmp_path: Path) -> None:
        ref = {
            "transcripts": [
                {
                    "id": "NM_TEST.1",
                    "gene_symbol": "TEST",
                    "strand": "+",
                    "sequence": "ATGCCCAAGGTGCTGCCC",
                    "cds_start": 1,
                    "cds_end": 18,
                    "exons": [{"number": 1, "start": 1, "end": 18}],
                }
            ],
            "proteins": {"NP_TEST.1": "MPKVLP"},
        }
        ref_file = tmp_path / "ref.json"
        ref_file.write_text(json.dumps(ref))

        norm = ferro_hgvs.Normalizer(reference_json=str(ref_file))
        out = norm.normalize("NM_TEST.1:c.1A>G")
        assert "NM_TEST.1" in out

    def test_bare_array_form_still_works(self, tmp_path: Path) -> None:
        ref = [
            {
                "id": "NM_TEST.1",
                "gene_symbol": "TEST",
                "strand": "+",
                "sequence": "ATGCCCAAGGTGCTGCCC",
                "cds_start": 1,
                "cds_end": 18,
                "exons": [{"number": 1, "start": 1, "end": 18}],
            }
        ]
        ref_file = tmp_path / "ref.json"
        ref_file.write_text(json.dumps(ref))

        norm = ferro_hgvs.Normalizer(reference_json=str(ref_file))
        out = norm.normalize("NM_TEST.1:c.1A>G")
        assert "NM_TEST.1" in out


class TestFromManifest:
    """from_manifest constructors backed by MultiFastaProvider."""

    def test_normalizer_from_manifest(self) -> None:
        norm = ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY))
        out = norm.normalize("NM_TEST.1:c.1A>G")
        assert "NM_TEST.1" in out

    def test_equivalence_checker_from_manifest(self) -> None:
        chk = ferro_hgvs.EquivalenceChecker.from_manifest(str(MANIFEST_TINY))
        v1 = ferro_hgvs.parse("NM_TEST.1:c.1A>G")
        v2 = ferro_hgvs.parse("NM_TEST.1:c.1A>G")
        result = chk.check(v1, v2)
        assert result.is_equivalent()

    def test_batch_processor_from_manifest(self) -> None:
        bp = ferro_hgvs.BatchProcessor.from_manifest(str(MANIFEST_TINY))
        result = bp.parse(["NM_TEST.1:c.1A>G", "INVALID"])
        assert result.success_count() == 1
        assert result.error_count() == 1

    def test_coordinate_mapper_from_manifest(self) -> None:
        cm = ferro_hgvs.CoordinateMapper.from_manifest(str(MANIFEST_TINY))
        assert cm.has_transcript("NM_TEST.1")
