"""End-to-end tests for ferro_hgvs.VariantProjector."""

import json
from pathlib import Path

import pytest

import ferro_hgvs


@pytest.fixture
def projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """Build a VariantProjector backed by a minimal hand-crafted reference JSON.

    The transcript NM_TEST.1 maps chr1 genomic positions 1000-1008 (1-based
    inclusive) with CDS sequence "ATGCGCTAA" (Met-Arg-Stop, 3 codons) on the
    plus strand.  g.1000 is the 'A' of the start codon (c.1).

    All coordinate fields follow the documented 1-based inclusive convention:
    - ``Exon.start`` / ``Exon.end``: 1-based transcript positions.
    - ``Exon.genomic_start`` / ``Exon.genomic_end``: 1-based genomic positions.
    - ``cds_start`` / ``cds_end``: 1-based transcript positions.

    ``CdotMapper.from_transcripts`` converts these to the internal cdot
    representation:
    - tx_start  = Exon.start - 1      (1-based → 0-based)
    - tx_end    = Exon.end            (1-based inclusive == 0-based exclusive)
    - g_start   = Exon.genomic_start  (HGVS-value-based, no conversion)
    - g_end     = Exon.genomic_end + 1  (1-based inclusive → exclusive)
    - cds_start = cds_start - 1       (1-based → 0-based)

    This produces the canonical cdot exon [1000, 1009, 0, 9] matching the
    Rust unit tests in ``src/project/projector.rs``.
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
                        "start": 1,  # 1-based, inclusive tx start
                        "end": 9,  # 1-based, inclusive tx end
                        "genomic_start": 1000,  # 1-based, inclusive genomic start
                        "genomic_end": 1008,  # 1-based, inclusive genomic end
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
        assert result.p_name == "NP_TEST.1:p.(Arg2Ser)"
        assert result.is_frameshift is False
        assert result.is_intronic is False
        assert result.is_utr is False

    def test_deletion_frameshift_protein(self, projector: ferro_hgvs.VariantProjector) -> None:
        # g.1003del deletes c.4 (C, first base of Arg codon CGC) — net=-1 → frameshift.
        # Mutated CDS: "ATGGCTAA" → Met-Ala-Stop → p.(Arg2Alafsx).
        result = projector.project("NC_000001.11:g.1003del", transcript="NM_TEST.1")
        assert result.c_name is not None
        assert "del" in result.c_name
        assert result.p_name is not None
        assert "fs" in result.p_name
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


@pytest.fixture
def long_projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """Build a VariantProjector backed by a longer transcript for indel p. tests.

    Transcript NM_LONG.1 on chr2, plus strand, genomic positions 2000-2023 (1-based).
    CDS sequence: "ATGCGCAAAGGGTTTTAA" (18 bp = 6 codons: Met-Arg-Lys-Gly-Phe-Stop).
    g.2000 = c.1 (A of start codon).
    g.2005 = c.6 (A of Lys codon AAA).
    g.2008 = c.9 (last base of Gly codon GGG).

    The full transcript is 24 bp: 3 bp 5'UTR "AAA" + 18 bp CDS + 3 bp 3'UTR "CCC".
    So tx positions 1-3 are UTR, 4-21 are CDS (cds_start=4, cds_end=21),
    and 22-24 are 3'UTR. Genomic: g.2000 = c.1 (tx pos 4 = first CDS base).
    """
    # CDS "ATGCGCAAAGGGTTTTAA":
    #   ATG = Met  (c.1-3,  g.2000-2002)
    #   CGC = Arg  (c.4-6,  g.2003-2005)
    #   AAA = Lys  (c.7-9,  g.2006-2008)
    #   GGG = Gly  (c.10-12, g.2009-2011)
    #   TTT = Phe  (c.13-15, g.2012-2014)
    #   TAA = Stop (c.16-18, g.2015-2017)
    #
    # Full tx seq (24 bp): "AAA" + "ATGCGCAAAGGGTTTTAA" + "CCC"
    # cds_start=4 (1-based), cds_end=21 (1-based).
    cds_seq = "ATGCGCAAAGGGTTTTAA"  # 18 bp
    full_seq = "AAA" + cds_seq + "CCC"  # 24 bp
    # Genomic layout: 3 bp UTR at g.2000-2002, CDS at g.2003-2020, UTR at g.2021-2023.
    # Actually: tx pos 1-3 = UTR (genomic 2000-2002), tx pos 4-21 = CDS (genomic 2003-2020).
    # Genomic start of transcript = 2000, end = 2023 (1-based inclusive).
    fixture = {
        "transcripts": [
            {
                "id": "NM_LONG.1",
                "gene_symbol": "LONGGENE",
                "strand": "+",
                "sequence": full_seq,
                "cds_start": 4,  # 1-based tx pos of first CDS base
                "cds_end": 21,  # 1-based tx pos of last CDS base (stop codon included)
                "exons": [
                    {
                        "number": 1,
                        "start": 1,  # tx pos start (1-based)
                        "end": 24,  # tx pos end (1-based)
                        "genomic_start": 2000,  # g.2000 = 'A' of 5'UTR
                        "genomic_end": 2023,  # g.2023 = last 'C' of 3'UTR
                    }
                ],
                "chromosome": "chr2",
                "genomic_start": 2000,
                "genomic_end": 2023,
            }
        ],
        "genomic_sequences": {
            "chr2": "N" * 2000 + full_seq + "N" * 100,
        },
    }
    path = tmp_path / "long_transcripts.json"
    path.write_text(json.dumps(fixture))
    return ferro_hgvs.VariantProjector(reference_json=str(path))


class TestIndelProteinNomenclature:
    """End-to-end p. nomenclature tests for indel projections on NM_LONG.1.

    CDS "ATGCGCAAAGGGTTTTAA": Met1-Arg2-Lys3-Gly4-Phe5-Stop.
    Genomic positions (1-based):
      g.2003 = c.1 = Met codon start
      g.2006 = c.4 = Arg codon start (CGC)
      g.2009 = c.7 = Lys codon start (AAA)
      g.2012 = c.10 = Gly codon start (GGG)
      g.2015 = c.13 = Phe codon start (TTT)
      g.2018 = c.16 = Stop codon start (TAA)
    """

    def test_del_whole_codon_single_aa(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Delete Arg codon (c.4_6 = g.2006_2008del): p.(Arg2del).
        result = long_projector.project("NC_000002.12:g.2006_2008del", transcript="NM_LONG.1")
        assert result.c_name is not None
        assert "del" in result.c_name
        assert result.p_name == "NP_LONG.1:p.(Arg2del)"
        assert result.is_frameshift is False

    def test_del_two_whole_codons(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Delete Arg+Lys codons (c.4_9 = g.2006_2014del): p.(Arg2_Lys3del).
        result = long_projector.project("NC_000002.12:g.2006_2011del", transcript="NM_LONG.1")
        assert result.c_name is not None
        assert result.p_name == "NP_LONG.1:p.(Arg2_Lys3del)"
        assert result.is_frameshift is False

    def test_del_single_base_frameshift(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Delete one base from Arg codon (g.2006del = c.4del): frameshift.
        result = long_projector.project("NC_000002.12:g.2006del", transcript="NM_LONG.1")
        assert result.p_name is not None
        assert "fs" in result.p_name
        assert "Arg2" in result.p_name
        assert result.is_frameshift is True

    def test_dup_whole_codon(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Duplicate Arg codon (g.2006_2008dup): p.(Arg2dup).
        result = long_projector.project("NC_000002.12:g.2006_2008dup", transcript="NM_LONG.1")
        assert result.c_name is not None
        assert "dup" in result.c_name
        assert result.p_name == "NP_LONG.1:p.(Arg2dup)"
        assert result.is_frameshift is False

    def test_dup_single_base_frameshift(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Duplicate one base (g.2006dup = c.4dup): frameshift.
        result = long_projector.project("NC_000002.12:g.2006dup", transcript="NM_LONG.1")
        assert result.p_name is not None
        assert "fs" in result.p_name
        assert result.is_frameshift is True

    def test_ins_three_bases_codon_boundary(
        self, long_projector: ferro_hgvs.VariantProjector
    ) -> None:
        # Insert GGG at codon boundary c.3_4 (g.2005_2006insGGG):
        # insertion between Met (c.1-3) and Arg (c.4-6) → p.(Met1_Arg2insGly).
        result = long_projector.project("NC_000002.12:g.2005_2006insGGG", transcript="NM_LONG.1")
        assert result.c_name is not None
        assert "ins" in result.c_name
        assert result.p_name is not None
        assert "ins" in result.p_name
        assert result.is_frameshift is False

    def test_ins_single_base_frameshift(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Insert one base (g.2005_2006insA): net=+1 → frameshift.
        result = long_projector.project("NC_000002.12:g.2005_2006insA", transcript="NM_LONG.1")
        assert result.p_name is not None
        assert "fs" in result.p_name
        assert result.is_frameshift is True

    def test_inversion_whole_codon(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Invert Arg codon (g.2006_2008inv = c.4_6inv = CGC → GCG = Ala):
        # p.(Arg2delinsAla).
        result = long_projector.project("NC_000002.12:g.2006_2008inv", transcript="NM_LONG.1")
        assert result.c_name is not None
        assert "inv" in result.c_name
        assert result.p_name is not None
        assert "Arg2" in result.p_name
        assert "delins" in result.p_name
        assert "Ala" in result.p_name
        assert result.is_frameshift is False

    def test_stop_codon_readthrough(self, long_projector: ferro_hgvs.VariantProjector) -> None:
        # Replace stop codon TAA (g.2018_2020) with TGG (Trp):
        # p.(Ter6Trpext*?).
        result = long_projector.project("NC_000002.12:g.2018_2020delinsTGG", transcript="NM_LONG.1")
        assert result.p_name is not None
        assert "Ter6" in result.p_name
        assert "ext" in result.p_name


@pytest.fixture
def two_tx_projector(tmp_path: Path) -> ferro_hgvs.VariantProjector:
    """VariantProjector with two transcripts overlapping chr1:1000-1008.

    NM_TX1.1 and NM_TX2.1 share the same genomic region and CDS sequence
    "ATGCGCTAA" on chr1. NM_TX2.1 is registered as MANE Select via the
    reference JSON (the projector's MANE annotation comes from the Projector
    built at VariantProjector construction time and uses the underlying
    cdot mapper).

    Because the Python VariantProjector does not yet expose MANE annotation,
    tests using this fixture only verify that *both* transcripts are returned;
    they do not assert ordering by MANE status.
    """
    fixture = {
        "transcripts": [
            {
                "id": "NM_TX1.1",
                "gene_symbol": "GENE1",
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
            },
            {
                "id": "NM_TX2.1",
                "gene_symbol": "GENE1",
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
            },
        ],
        "genomic_sequences": {
            "chr1": "N" * 1000 + "ATGCGCTAA" + "N" * 100,
        },
    }
    path = tmp_path / "two_tx.json"
    path.write_text(json.dumps(fixture))
    return ferro_hgvs.VariantProjector(reference_json=str(path))


class TestProjectAll:
    """Smoke tests for project_all and project_normalized_all."""

    def test_project_all_returns_both_transcripts(
        self, two_tx_projector: ferro_hgvs.VariantProjector
    ) -> None:
        results = two_tx_projector.project_all("chr1:g.1003C>A")
        assert len(results) == 2, f"expected 2 projections, got {len(results)}"

    def test_project_all_no_overlap_returns_empty(
        self, two_tx_projector: ferro_hgvs.VariantProjector
    ) -> None:
        # g.5000 is far outside all transcripts.
        results = two_tx_projector.project_all("NC_000001.11:g.5000A>G")
        assert results == [], "expected empty list for non-overlapping position"

    def test_project_all_projections_have_c_name(
        self, two_tx_projector: ferro_hgvs.VariantProjector
    ) -> None:
        results = two_tx_projector.project_all("chr1:g.1003C>A")
        for r in results:
            assert r.c_name is not None
            assert ":c." in r.c_name

    def test_project_normalized_all_same_result_as_project_all(
        self, two_tx_projector: ferro_hgvs.VariantProjector
    ) -> None:
        hgvs_string = "chr1:g.1003C>A"
        via_all = two_tx_projector.project_all(hgvs_string)
        # Parse only; g.1003C>A is already canonical (a CDS substitution can't
        # 3'-shift), so normalize() would be a no-op and skipping it exercises
        # the fast-path contract: project_normalized_all on a parsed-but-not-
        # normalized canonical variant equals project_all on the source string.
        variant = ferro_hgvs.parse(hgvs_string)
        via_normalized_all = two_tx_projector.project_normalized_all(variant)
        assert len(via_all) == len(via_normalized_all)
        for a, b in zip(via_all, via_normalized_all, strict=True):
            assert a.transcript_id == b.transcript_id
            assert a.c_name == b.c_name


class TestProjectNormalized:
    """Smoke tests for project_normalized."""

    def test_project_normalized_matches_project(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        hgvs_string = "NC_000001.11:g.1003C>A"
        via_project = projector.project(hgvs_string, "NM_TEST.1")
        # project_normalized skips re-normalization; for an already-canonical
        # string the results must be identical.
        variant = ferro_hgvs.parse(hgvs_string)
        via_normalized = projector.project_normalized(variant, "NM_TEST.1")
        assert via_project.c_name == via_normalized.c_name
        assert via_project.p_name == via_normalized.p_name
        assert via_project.transcript_id == via_normalized.transcript_id

    def test_project_normalized_returns_variant_projection(
        self, projector: ferro_hgvs.VariantProjector
    ) -> None:
        variant = ferro_hgvs.parse("NC_000001.11:g.1003C>A")
        result = projector.project_normalized(variant, "NM_TEST.1")
        assert isinstance(result, ferro_hgvs.VariantProjection)
        assert result.c_name is not None
