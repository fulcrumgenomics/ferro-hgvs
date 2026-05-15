"""End-to-end tests for ferro_hgvs.VariantProjector."""

import json
from pathlib import Path

import pytest

import ferro_hgvs


@pytest.fixture
def projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """Build a VariantProjector backed by a minimal hand-crafted reference JSON.

    The transcript NM_TEST.1 maps chr1 genomic positions 1001-1009 (1-based,
    0-based half-open [1000, 1009)) with CDS sequence "ATGCGCTAA"
    (Met-Arg-Stop, 3 codons) on the plus strand.

    Exon coordinate notes
    ---------------------
    The MockProvider uses ``Exon.start`` directly as the cdot ``tx_start``
    field (0-based), so we set ``start=0`` and ``end=8`` to match the
    0-based cdot convention used internally by ``CdotMapper.from_transcripts``.
    ``genomic_start`` and ``genomic_end`` are 1-based inclusive; they are
    converted to 0-based half-open by subtracting 1 from the start only.
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
                        "start": 0,
                        "end": 8,
                        "genomic_start": 1001,
                        "genomic_end": 1009,
                    }
                ],
                "chromosome": "chr1",
                "genomic_start": 1000,
                "genomic_end": 1008,
            }
        ],
        "genomic_sequences": {
            "chr1": "N" * 1000 + "ATGCGCTAA" + "N" * 100,
        },
    }
    path = tmp_path / "transcripts.json"
    path.write_text(json.dumps(fixture))
    return ferro_hgvs.VariantProjector(reference_json=str(path))


class TestVariantProjector:
    def test_missense_substitution_plus_strand(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        # g.1003C>A: codon 2 of "ATGCGCTAA" is "CGC" (Arg); the 4th base (c.4)
        # is C.  Changing C→A gives "AGC" (Ser) — missense.
        result = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        assert result.transcript_id == "NM_TEST.1"
        assert result.gene_symbol == "TESTGENE"
        assert result.c_name is not None
        assert ":c.4C>A" in result.c_name
        assert result.p_name == "NP_TEST.1(TESTGENE):p.(Arg2Ser)"
        assert result.is_frameshift is False
        assert result.is_intronic is False
        assert result.is_utr is False

    def test_deletion_no_protein_is_frameshift(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        # 1-base deletion → frameshift; p. is None (indel p. nomenclature deferred).
        result = projector.project("NC_000001.11:g.1003del", transcript="NM_TEST.1")
        assert result.c_name is not None
        assert "del" in result.c_name
        assert result.p_name is None
        assert result.is_frameshift is True

    def test_unknown_transcript_raises(self, projector: ferro_hgvs.VariantProjector) -> None:
        with pytest.raises(RuntimeError):
            projector.project("NC_000001.11:g.1003C>A", transcript="NM_NOPE.99")

    def test_no_overlap_raises(self, projector: ferro_hgvs.VariantProjector) -> None:
        # chr1:5000 is far outside the transcript at positions 1001-1009.
        with pytest.raises(RuntimeError, match="overlap"):
            projector.project("NC_000001.11:g.5000A>G", transcript="NM_TEST.1")

    def test_repr_includes_g_name(self, projector: ferro_hgvs.VariantProjector) -> None:
        result = projector.project("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        rep = repr(result)
        assert "VariantProjection" in rep
        # repr should include either the g_name or c_name
        assert ":g." in rep or ":c." in rep
