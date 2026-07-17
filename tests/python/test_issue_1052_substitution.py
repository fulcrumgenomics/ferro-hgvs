"""Issue #1052: end-to-end repro — substitution reference-base validation is
surfaced via ``normalize_with_warnings`` in lenient mode, and rejected by
strict mode alongside the already-validated deletion half (#1042)."""

import tempfile

import pytest

import ferro_hgvs


def _genomic_reference() -> str:
    """Build a tiny genomic transcripts.json via convert_gff (mirrors the issue)."""
    seq = "GGGG" + "ATGTCTACGCTGGGATAA" + "CCCC"  # 1-based g.8 == 'T'
    ctg, gs, ge = "TEST", 5, 22
    d = tempfile.mkdtemp()
    fasta, gff, out = f"{d}/r.fa", f"{d}/r.gff3", f"{d}/transcripts.json"
    with open(fasta, "w") as fh:
        fh.write(f">{ctg}\n{seq}\n")
    with open(gff, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {ctg} 1 {len(seq)}\n")
        fh.write(
            f"{ctg}\tx\tregion\t1\t{len(seq)}\t.\t+\t.\tID={ctg}:1..{len(seq)};Name={ctg};chromosome={ctg};gbkey=Src;mol_type=genomic DNA\n"
        )
        fh.write(f"{ctg}\tx\tgene\t{gs}\t{ge}\t.\t+\t0\tID={ctg}-gene;Name={ctg}-gene\n")
        fh.write(
            f"{ctg}\tx\tmRNA\t{gs}\t{ge}\t.\t+\t0\tID={ctg}-gene.1;Parent={ctg}-gene;transcript_id=1\n"
        )
        fh.write(
            f"{ctg}\tx\tCDS\t{gs}\t{ge}\t.\t+\t0\tID={ctg}-gene.1.cds;Parent={ctg}-gene.1;protein_id=1\n"
        )
    ferro_hgvs.convert_gff(
        ferro_hgvs.ConvertGffConfig(gff=gff, fasta=fasta, output=out, emit_genomic_sequences=True)
    )
    return out


def test_wrong_ref_substitution_surfaces_refseq_mismatch():
    out = _genomic_reference()
    norm = ferro_hgvs.Normalizer(reference_json=out)
    assert norm.has_genomic_data() is True
    res = norm.normalize_with_warnings("TEST:g.8A>C")  # g.8 is 'T'
    assert norm.normalize("TEST:g.8A>C") == "TEST:g.8A>C"  # unchanged in lenient
    assert any("REFSEQ_MISMATCH" in str(w) for w in res.warnings), list(res.warnings)


def test_correct_ref_substitution_no_warning():
    out = _genomic_reference()
    norm = ferro_hgvs.Normalizer(reference_json=out)
    res = norm.normalize_with_warnings("TEST:g.8T>C")
    assert list(res.warnings) == []


def test_strict_rejects_wrong_ref_sub_and_del():
    out = _genomic_reference()
    strict = ferro_hgvs.ErrorConfig.strict()
    norm = ferro_hgvs.Normalizer(reference_json=out, error_config=strict)
    with pytest.raises(ferro_hgvs.NormalizationError):
        norm.normalize("TEST:g.8A>C")
    with pytest.raises(ferro_hgvs.NormalizationError):
        norm.normalize("TEST:g.8delA")  # del half already validated (#1042)
