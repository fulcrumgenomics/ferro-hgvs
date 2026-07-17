"""Tests for #1050: protein render-style bindings.

Covers the constructor default (Option A: ``VariantProjector(protein_stop=...,
amino_acid_code=...)``) and the per-call override (Option B:
``VariantProjection.p_name_styled(protein_stop=..., amino_acid_code=...)``).

The stop-codon axis (``protein_stop``: "ter" vs "star") is only
distinguishable in three-letter mode, and only on a variant whose p. name
actually carries a stop token — a plain missense substitution never does. So
``G_FS_INPUT`` (a frameshift, ``NC_000001.11:g.1003del`` → p. name containing
``...fsTer?)``/``...fs*?)``) is used everywhere the ``protein_stop`` axis needs
to be exercised; ``G_INPUT`` (missense, ``NC_000001.11:g.1003C>A`` →
``p.(Arg2Ser)``) is kept only as a plain default-rendering sanity check. Both
reuse ``NM_TEST.1`` and its minimal single-exon reference JSON, taken verbatim
from ``tests/python/test_variant_projection.py``'s ``projector`` fixture (whose
``test_missense_substitution_plus_strand`` / ``test_deletion_frameshift_protein``
assert the same two p. names off the same fixture).

``G_UTR_INPUT`` / ``TRANSCRIPT_LONG`` (for the genuine-``None`` case) reuse
``NM_LONG.1`` and its reference JSON from the same file's ``long_projector``
fixture, which — unlike ``NM_TEST.1`` (no UTR: cds spans the whole transcript)
— has a 5'UTR (tx positions 1-3, cds_start=4), so a genomic substitution there
projects to a real UTR variant with a genuinely ``None`` p_name (verified
empirically: ``g.2001A>G`` gives ``is_utr=True``, ``p_name=None``).

The built-in test-data projector (``VariantProjector()`` with no
``reference_json``) is *not* usable for any of this: every transcript in
``JsonProvider::with_test_data`` (src/reference/mock.rs) carries no
``chromosome``/``genomic_start``/``genomic_end``, so ``.project()`` against it
always raises "Reference not found" regardless of input — there is no
(genomic_string, transcript) pair in the built-in test data that yields any
p_name at all.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs

G_INPUT = "NC_000001.11:g.1003C>A"
G_FS_INPUT = "NC_000001.11:g.1003del"
TRANSCRIPT = "NM_TEST.1"

G_UTR_INPUT = "NC_000002.12:g.2001A>G"
TRANSCRIPT_LONG = "NM_LONG.1"


def _reference_json(tmp_path: Path) -> str:
    """Write the minimal single-exon reference used by test_variant_projection.py's
    ``projector`` fixture (NM_TEST.1, CDS "ATGCGCTAA", g.1000-1008 on chr1,
    codon 2 "CGC"=Arg). ``g.1003C>A`` mutates it to "AGC"=Ser (missense,
    p.(Arg2Ser), no stop token); ``g.1003del`` deletes the codon's first base,
    shifting the frame into "GCT AA" (frameshift, p.(Arg2AlafsTer?), carries a
    stop token). Returns the fixture's path.
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
    return str(path)


def _long_reference_json(tmp_path: Path) -> str:
    """Write the NM_LONG.1 reference used by test_variant_projection.py's
    ``long_projector`` fixture: 3bp 5'UTR + 18bp CDS + 3bp 3'UTR on chr2,
    g.2000-2023. Unlike NM_TEST.1 (no UTR), this carries a real 5'UTR (tx
    positions 1-3, cds_start=4), so a genomic variant there (``g.2001A>G``)
    projects to a genuine UTR variant with p_name == None. Returns the
    fixture's path.
    """
    cds_seq = "ATGCGCAAAGGGTTTTAA"  # 18 bp
    full_seq = "AAA" + cds_seq + "CCC"  # 24 bp
    fixture = {
        "transcripts": [
            {
                "id": "NM_LONG.1",
                "gene_symbol": "LONGGENE",
                "strand": "+",
                "sequence": full_seq,
                "cds_start": 4,
                "cds_end": 21,
                "exons": [
                    {
                        "number": 1,
                        "start": 1,
                        "end": 24,
                        "genomic_start": 2000,
                        "genomic_end": 2023,
                    }
                ],
                "chromosome": "chr2",
                "genomic_start": 2000,
                "genomic_end": 2023,
            }
        ],
        "genomic_sequences": {
            "chr2": "N" * 1999 + full_seq + "N" * 100,
        },
    }
    path = tmp_path / "long_transcripts.json"
    path.write_text(json.dumps(fixture))
    return str(path)


def _proj(
    tmp_path: Path, hgvs_string: str = G_INPUT, **kwargs: object
) -> ferro_hgvs.VariantProjection:
    projector = ferro_hgvs.VariantProjector(reference_json=_reference_json(tmp_path), **kwargs)
    return projector.project(hgvs_string, TRANSCRIPT)


def test_default_missense_unchanged(tmp_path: Path) -> None:
    # Sanity/backward-compat: a plain missense (no stop token) default-renders
    # exactly as it did before #1050 — matches
    # test_variant_projection.py::test_missense_substitution_plus_strand.
    p = _proj(tmp_path).p_name
    assert p == "NM_TEST.1:p.(Arg2Ser)"


def test_default_is_three_letter_ter(tmp_path: Path) -> None:
    p = _proj(tmp_path, G_FS_INPUT).p_name
    assert p == "NM_TEST.1:p.(Arg2AlafsTer?)"
    assert "Ter" in p
    assert "*" not in p


def test_constructor_protein_stop_star_switches_stop_token(tmp_path: Path) -> None:
    # Same frameshift variant; only the stop spelling changes (three-letter
    # residue codes are kept) — proves the protein_stop axis actually does
    # something, not just amino_acid_code.
    p = _proj(tmp_path, G_FS_INPUT, protein_stop="star").p_name
    assert p == "NM_TEST.1:p.(Arg2Alafs*?)"
    assert "*" in p
    assert "Ter" not in p
    assert "Arg" in p  # residues still three-letter


def test_constructor_option_a_one_letter_star(tmp_path: Path) -> None:
    p = _proj(tmp_path, G_FS_INPUT, protein_stop="star", amino_acid_code="one").p_name
    assert p == "NM_TEST.1:p.(R2Afs*?)"
    assert "*" in p
    assert "Ter" not in p
    assert "Arg" not in p and "Ala" not in p  # three-letter codes gone


def test_per_call_option_b_overrides_default(tmp_path: Path) -> None:
    proj = _proj(tmp_path, G_FS_INPUT)  # default projector
    default_p = proj.p_name
    assert default_p == "NM_TEST.1:p.(Arg2AlafsTer?)"
    styled = proj.p_name_styled(protein_stop="star", amino_acid_code="one")
    assert styled == "NM_TEST.1:p.(R2Afs*?)"
    assert styled != default_p


def test_option_b_none_falls_back_to_constructor(tmp_path: Path) -> None:
    proj = _proj(tmp_path, G_FS_INPUT, protein_stop="star", amino_acid_code="one")
    # p_name_styled() with no args returns the constructor style
    assert proj.p_name_styled() == proj.p_name == "NM_TEST.1:p.(R2Afs*?)"


def test_option_b_partial_override(tmp_path: Path) -> None:
    proj = _proj(tmp_path, G_FS_INPUT, protein_stop="star", amino_acid_code="one")
    assert proj.p_name == "NM_TEST.1:p.(R2Afs*?)"
    # override only aa_code back to three-letter, keep constructor's star stop
    styled = proj.p_name_styled(amino_acid_code="three")
    assert styled == "NM_TEST.1:p.(Arg2Alafs*?)"
    assert "Arg" in styled  # three-letter residues restored
    assert "*" in styled  # constructor's star stop retained
    assert "Ter" not in styled


@pytest.mark.parametrize("value", ["nope", "Ter", "TER", "*", "", "THREE"])
def test_invalid_protein_stop_raises(value: str) -> None:
    # Validated before any reference is loaded, so no reference_json is needed.
    # Includes formerly-accepted aliases ("Ter"/"TER"/"*") that #1050 review
    # dropped so the accepted set matches the documented "ter"|"star" contract
    # exactly.
    with pytest.raises(ValueError):
        ferro_hgvs.VariantProjector(protein_stop=value)


@pytest.mark.parametrize("value", ["nope", "3", "1", "", "ONE"])
def test_invalid_amino_acid_code_raises(value: str) -> None:
    # Includes formerly-accepted aliases ("3"/"1") dropped for the same reason.
    with pytest.raises(ValueError):
        ferro_hgvs.VariantProjector(amino_acid_code=value)


def test_p_name_styled_none_passthrough(tmp_path: Path) -> None:
    # A genuine None-p_name projection: NM_LONG.1 has a real 5'UTR (unlike
    # NM_TEST.1, which has none), so a genomic substitution there projects to
    # a UTR variant with p_name == None. p_name_styled must stay None too.
    projector = ferro_hgvs.VariantProjector(reference_json=_long_reference_json(tmp_path))
    proj = projector.project(G_UTR_INPUT, TRANSCRIPT_LONG)
    assert proj.is_utr is True
    assert proj.p_name is None
    assert proj.p_name_styled(protein_stop="star") is None
