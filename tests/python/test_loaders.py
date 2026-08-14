"""Tests for Normalizer/EquivalenceChecker/BatchProcessor/CoordinateMapper loader paths."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

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
            "genomic_sequences": {"NC_TEST.1": "ACGTACGT"},
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
        # The tiny manifest is a FASTA and no cdot, so it serves NM_TEST.1's
        # BASES and resolves no CDS. Normalize on the `n.` axis, which needs no
        # CDS — `background/numbering.md:52` numbers it over the reference
        # sequence's own nucleotides — so this stays a test of the loader
        # plumbing rather than of CDS resolution.
        norm = ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY))
        out = norm.normalize("NM_TEST.1:n.1A>G")
        # The exact value, not just containment: the fixture is deterministic
        # (`transcripts/test.fna` is a single `NM_TEST.1` / `ATGCATGCAT`), so
        # `n.1` is the leading `A` and a substitution has no run to shift along.
        # `assert "NM_TEST.1" in out` would also be satisfied by the input
        # echoed back, which is the state #1870 exists to make unrepresentable.
        assert out == "NM_TEST.1:n.1A>G"

    def test_manifest_without_cds_refuses_the_coding_axis(self) -> None:
        # The sibling half, and the reason the row above moved off `c.`: the
        # decided `rulings[c-description-against-an-unresolvable-cds-is-refused]`
        # refuses a `c.` description against a transcript whose CDS the
        # reference cannot resolve, in every mode, rather than returning the
        # input labelled ok. This asserts the bindings surface that refusal as a
        # typed error rather than swallowing it.
        norm = ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY))
        with pytest.raises(ferro_hgvs.NormalizationError) as excinfo:
            norm.normalize("NM_TEST.1:c.1A>G")
        assert "has no CDS start" in str(excinfo.value)

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
