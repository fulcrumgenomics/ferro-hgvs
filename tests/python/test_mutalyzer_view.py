"""End-to-end tests for ferro_hgvs.VariantProjector.mutalyzer.

A Mutalyzer-style, multi-axis view: the normalized form plus the g./c./n./r./p.
axes bundled into one MutalyzerResult. Interface-compatible, not
output-compatible — the values are ferro's own, so these assertions pin ferro's
behaviour, never agreement with real Mutalyzer.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs


def _write_reference(tmp_path: Path) -> Path:
    """Write the minimal NM_TEST.1 reference and return its path.

    NM_TEST.1 maps chr1 g.1000-1008 (1-based inclusive) with CDS "ATGCGCTAA"
    (Met-Arg-Stop) on the plus strand; g.1000 is c.1. Identical to the fixture
    in test_variant_projection.py so the two suites pin the same reference.
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
        "genomic_sequences": {
            "chr1": "N" * 999 + "ATGCGCTAA" + "N" * 100,
        },
    }
    path = tmp_path / "transcripts.json"
    path.write_text(json.dumps(fixture))
    return path


@pytest.fixture
def projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """A default-style VariantProjector over the minimal NM_TEST.1 reference."""
    return ferro_hgvs.VariantProjector(reference_json=str(_write_reference(tmp_path)))


class TestMutalyzerView:
    def test_bundles_every_derivable_axis(self, projector: ferro_hgvs.VariantProjector) -> None:
        # g.1003C>A: codon 2 "CGC" (Arg) -> "AGC" (Ser); missense c.4C>A.
        result = projector.mutalyzer("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")

        assert isinstance(result, ferro_hgvs.MutalyzerResult)
        assert result.input == "NC_000001.11:g.1003C>A"
        assert result.transcript_id == "NM_TEST.1"
        assert result.gene_symbol == "TESTGENE"

        # A genomic input has a genomic axis, and the coding/protein axes are
        # derived onto the transcript.
        assert result.genomic is not None
        assert result.coding is not None
        assert ":c.4C>A" in result.coding
        assert result.protein == "NM_TEST.1:p.(Arg2Ser)"

        # `normalized` is the input's own axis (genomic here), normalized.
        assert result.normalized is not None
        assert result.normalized.startswith("NC_000001.11:g.")

        assert result.has_warnings() is False
        assert result.warnings == []

    def test_str_renders_an_aligned_block(self, projector: ferro_hgvs.VariantProjector) -> None:
        result = projector.mutalyzer("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        text = str(result)
        assert "Input:" in text
        assert "Genomic:" in text
        assert "Protein:" in text
        # The input string appears in the rendered block.
        assert "NC_000001.11:g.1003C>A" in text

    def test_repr_names_the_type(self, projector: ferro_hgvs.VariantProjector) -> None:
        result = projector.mutalyzer("NC_000001.11:g.1003C>A", transcript="NM_TEST.1")
        assert "MutalyzerResult" in repr(result)

    def test_unknown_transcript_raises(self, projector: ferro_hgvs.VariantProjector) -> None:
        # The bindings promise ProjectionError specifically (a subclass of
        # ValueError and RuntimeError); assert that, not a broad superclass.
        with pytest.raises(ferro_hgvs.ProjectionError):
            projector.mutalyzer("NC_000001.11:g.1003C>A", transcript="NM_NOPE.99")

    def test_no_overlap_raises(self, projector: ferro_hgvs.VariantProjector) -> None:
        with pytest.raises(ferro_hgvs.ProjectionError, match="overlap"):
            projector.mutalyzer("NC_000001.11:g.5000A>G", transcript="NM_TEST.1")

    def test_protein_axis_respects_render_style(self, tmp_path: Path) -> None:
        # g.1003C>A is a missense (Arg2Ser) with a p. axis, so its spelling
        # exposes the projector's configured protein style. The mutalyzer view
        # must honour that style, not fall back to the default three-letter form.
        reference = str(_write_reference(tmp_path))
        variant = "NC_000001.11:g.1003C>A"

        default = ferro_hgvs.VariantProjector(reference_json=reference)
        one_letter = ferro_hgvs.VariantProjector(reference_json=reference, amino_acid_code="one")

        default_result = default.mutalyzer(variant, transcript="NM_TEST.1")
        styled_result = one_letter.mutalyzer(variant, transcript="NM_TEST.1")

        assert default_result.protein == "NM_TEST.1:p.(Arg2Ser)"
        assert styled_result.protein == "NM_TEST.1:p.(R2S)"
