"""Issue #1035 (sibling command): `ferro build-transcript` is exposed in the
Python bindings as `ferro_hgvs.build_transcript` / `ferro_hgvs.BuildTranscriptConfig`.

A synthetic single-construct reference (own FASTA + CDS bounds) can be built and
normalized entirely in-process via the PyPI wheel — no shelling out to the CLI.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs


def _write_fasta(tmp_path: Path, name: str = "construct1", length: int = 60, seq: str = "") -> Path:
    fasta = tmp_path / "construct.fa"
    if not seq:
        seq = "".join("ACGT"[i % 4] for i in range(length))
    fasta.write_text(f">{name}\n{seq}\n")
    return fasta


def test_symbols_are_exported() -> None:
    assert "build_transcript" in ferro_hgvs.__all__
    assert "BuildTranscriptConfig" in ferro_hgvs.__all__
    assert "BuildTranscriptReport" in ferro_hgvs.__all__
    assert hasattr(ferro_hgvs, "build_transcript")


def test_build_writes_transcripts_json_to_output_path(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)
    out = tmp_path / "transcripts.json"

    cfg = ferro_hgvs.BuildTranscriptConfig(
        fasta=str(fasta), cds_start=1, cds_end=60, output=str(out)
    )
    report = ferro_hgvs.build_transcript(cfg)

    assert report.output_path == str(out)
    assert report.transcripts_json is None
    assert report.transcript_id == "construct1"

    doc = json.loads(out.read_text())
    assert doc["version"] == "1.0"
    assert doc["genome_build"] == "GRCh38"
    tx = doc["transcripts"][0]
    assert tx["id"] == "construct1"
    assert tx["cds_start"] == 1
    assert tx["cds_end"] == 60
    assert tx["strand"] == "+"


def test_build_returns_json_in_memory_when_output_is_none(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)

    report = ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=60)
    )
    assert report.output_path is None
    assert report.transcripts_json is not None
    assert json.loads(report.transcripts_json)["transcripts"][0]["id"] == "construct1"


def test_written_file_equals_in_memory_json(tmp_path: Path) -> None:
    """build-transcript writes the pretty JSON with NO trailing newline (the CLI
    uses `std::fs::write`), so the file equals the in-memory JSON exactly."""
    fasta = _write_fasta(tmp_path)

    in_memory = ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=60)
    )
    out = tmp_path / "transcripts.json"
    ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=60, output=str(out))
    )
    assert out.read_text() == in_memory.transcripts_json


def test_custom_id_and_gene(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)
    report = ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(
            fasta=str(fasta), cds_start=1, cds_end=60, id="MYTX.1", gene="GENE1"
        )
    )
    tx = json.loads(report.transcripts_json)["transcripts"][0]
    assert report.transcript_id == "MYTX.1"
    assert tx["id"] == "MYTX.1"
    assert tx["gene_symbol"] == "GENE1"


def test_minus_strand_reverse_complements(tmp_path: Path) -> None:
    # Use an ASYMMETRIC forward sequence. An "ACGT" repeat is its own reverse
    # complement, so it cannot distinguish "revcomp applied" from "revcomp
    # skipped" — this block sequence's revcomp differs from its forward strand.
    forward = "A" * 20 + "C" * 20 + "G" * 20
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    expected = "".join(complement[b] for b in reversed(forward))
    # Guard: the fixture actually exercises reverse-complementation.
    assert expected != forward

    fasta = _write_fasta(tmp_path, seq=forward)
    report = ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=60, strand="-")
    )
    tx = json.loads(report.transcripts_json)["transcripts"][0]
    assert tx["strand"] == "-"
    assert tx["sequence"] == expected
    assert tx["sequence"] != forward


def test_emit_genomic_sequences_makes_reference_genome_capable(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)
    out = tmp_path / "genome_capable.json"

    report = ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(
            fasta=str(fasta),
            cds_start=1,
            cds_end=60,
            output=str(out),
            emit_genomic_sequences=True,
        )
    )
    assert report.emitted_genomic_bytes == 60

    doc = json.loads(out.read_text())
    assert doc["genomic_sequences"]["construct1"] == "".join("ACGT"[i % 4] for i in range(60))

    normalizer = ferro_hgvs.Normalizer(reference_json=str(out))
    assert normalizer.has_genomic_data() is True


def test_built_reference_normalizes_in_process(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)
    out = tmp_path / "transcripts.json"
    ferro_hgvs.build_transcript(
        ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=60, output=str(out))
    )

    tx = json.loads(out.read_text())["transcripts"][0]
    ref_base = tx["sequence"][tx["cds_start"] - 1]  # base at c.1

    normalizer = ferro_hgvs.Normalizer(reference_json=str(out))
    normalized = normalizer.normalize(f"construct1:c.1{ref_base}>N")
    assert normalized.startswith("construct1:c.1")


def test_invalid_strand_raises_value_error(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)
    with pytest.raises(ValueError, match="strand"):
        ferro_hgvs.build_transcript(
            ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=60, strand="*")
        )


def test_invalid_cds_bounds_raises_value_error(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)  # 60 bases
    with pytest.raises(ValueError, match="CDS bounds"):
        ferro_hgvs.build_transcript(
            ferro_hgvs.BuildTranscriptConfig(fasta=str(fasta), cds_start=1, cds_end=999)
        )


def test_missing_contig_raises_reference_data_error(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path)
    with pytest.raises(ferro_hgvs.ReferenceDataError):
        ferro_hgvs.build_transcript(
            ferro_hgvs.BuildTranscriptConfig(
                fasta=str(fasta), cds_start=1, cds_end=60, contig="nope"
            )
        )
