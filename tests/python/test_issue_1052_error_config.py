"""Issue #1052: ErrorConfig can be passed to Normalizer / VariantProjector, and
strict mode rejects a wrong-reference-base variant."""

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


def test_normalizer_accepts_error_config_and_rejects_wrong_ref_sub():
    out = _genomic_reference()
    strict = ferro_hgvs.ErrorConfig.strict()
    norm = ferro_hgvs.Normalizer(reference_json=out, error_config=strict)
    with pytest.raises(ferro_hgvs.NormalizationError):
        norm.normalize("TEST:g.8A>C")  # g.8 is 'T'


def test_normalizer_lenient_default_unchanged():
    out = _genomic_reference()
    norm = ferro_hgvs.Normalizer(reference_json=out)  # no error_config
    assert norm.normalize("TEST:g.8A>C") == "TEST:g.8A>C"
    warnings = list(norm.normalize_with_warnings("TEST:g.8A>C").warnings)
    assert any("REFSEQ_MISMATCH" in str(w) for w in warnings), warnings


def test_projector_error_config_flows_through_to_normalization():
    # Constructor must accept the keyword AND actually apply it: a strict
    # projector must reject a wrong-reference g. variant during its
    # normalize-then-project step, while the default (lenient) projector must
    # not. Asserting both directions means the test fails if VariantProjector
    # silently discards error_config. g.8 is 'T'; assert 'A>C'.
    out = _genomic_reference()

    strict_projector = ferro_hgvs.VariantProjector(
        reference_json=out, error_config=ferro_hgvs.ErrorConfig.strict()
    )
    with pytest.raises(ferro_hgvs.ProjectionError):
        strict_projector.project_all("TEST:g.8A>C")

    # Default (lenient) projector: the same wrong-ref variant normalizes with a
    # warning and projects without raising — proving the strict rejection above
    # came from error_config, not from an unrelated failure.
    lenient_projector = ferro_hgvs.VariantProjector(reference_json=out)
    lenient_projector.project_all("TEST:g.8A>C")
