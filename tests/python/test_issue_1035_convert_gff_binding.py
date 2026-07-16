"""Issue #1035: `ferro convert-gff` is exposed in the Python bindings as
`ferro_hgvs.convert_gff` / `ferro_hgvs.ConvertGffConfig`, mirroring
`prepare_reference_data` / `PrepareConfig`.

A synthetic/local reference (own GFF3 + FASTA) can be built and normalized
entirely in-process via the PyPI wheel — no shelling out to the `ferro` CLI.
This exercises the exact Python surface a saturation-mutagenesis pipeline uses.
"""

import json
import warnings
from pathlib import Path

import pytest

import ferro_hgvs

# A single-exon transcript on chr1 with a CDS, mirroring the Rust integration
# fixtures. Built programmatically (no committed test-data files).
GFF3 = (
    "##gff-version 3\n"
    "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=g1;Name=GENE1\n"
    "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=g1;gene=GENE1\n"
    "chr1\t.\texon\t100\t500\t.\t+\t.\tParent=tx1\n"
    "chr1\t.\tCDS\t150\t450\t.\t+\t0\tParent=tx1\n"
    "chr1\t.\tstop_codon\t448\t450\t.\t+\t0\tParent=tx1\n"
)


def _write_gff(tmp_path: Path) -> Path:
    gff = tmp_path / "constructs.gff3"
    gff.write_text(GFF3)
    return gff


def _write_fasta(tmp_path: Path, name: str = "chr1", length: int = 520) -> Path:
    fasta = tmp_path / "constructs.fa"
    seq = "".join("ACGT"[i % 4] for i in range(length))
    fasta.write_text(f">{name}\n{seq}\n")
    return fasta


def test_symbols_are_exported() -> None:
    assert "convert_gff" in ferro_hgvs.__all__
    assert "ConvertGffConfig" in ferro_hgvs.__all__
    assert "ConvertGffReport" in ferro_hgvs.__all__
    assert hasattr(ferro_hgvs, "convert_gff")
    assert hasattr(ferro_hgvs, "ConvertGffConfig")
    assert hasattr(ferro_hgvs, "ConvertGffReport")


def test_convert_writes_transcripts_json_to_output_path(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)
    out = tmp_path / "transcripts.json"

    cfg = ferro_hgvs.ConvertGffConfig(gff=str(gff), output=str(out), build="GRCh38")
    report = ferro_hgvs.convert_gff(cfg)

    assert report.output_path == str(out)
    # When written to a file, the JSON is not also returned in-memory.
    assert report.transcripts_json is None
    assert report.transcript_count == 1

    doc = json.loads(out.read_text())
    assert doc["version"] == "1.0"
    assert doc["genome_build"] == "GRCh38"
    assert [t["id"] for t in doc["transcripts"]] == ["tx1"]


def test_convert_returns_json_in_memory_when_output_is_none(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)

    cfg = ferro_hgvs.ConvertGffConfig(gff=str(gff))  # output defaults to None
    report = ferro_hgvs.convert_gff(cfg)

    assert report.output_path is None
    assert report.transcripts_json is not None
    doc = json.loads(report.transcripts_json)
    assert doc["transcripts"][0]["id"] == "tx1"


def test_written_file_equals_in_memory_json_plus_newline(tmp_path: Path) -> None:
    """The file written for output=<path> is the in-memory JSON followed by a
    trailing newline — same bytes the CLI writes via `writeln!`."""
    gff = _write_gff(tmp_path)

    in_memory = ferro_hgvs.convert_gff(ferro_hgvs.ConvertGffConfig(gff=str(gff)))
    out = tmp_path / "transcripts.json"
    ferro_hgvs.convert_gff(ferro_hgvs.ConvertGffConfig(gff=str(gff), output=str(out)))

    assert out.read_text() == in_memory.transcripts_json + "\n"


def test_built_reference_normalizes_in_process(tmp_path: Path) -> None:
    """The whole point of #1035: build a reference with convert_gff (no
    subprocess) and normalize a real variant against it."""
    gff = _write_gff(tmp_path)
    fasta = _write_fasta(tmp_path)  # FASTA-backed so the transcript sequence is real
    out = tmp_path / "transcripts.json"
    ferro_hgvs.convert_gff(
        ferro_hgvs.ConvertGffConfig(
            gff=str(gff), fasta=str(fasta), output=str(out), validate_fasta=False
        )
    )

    tx = json.loads(out.read_text())["transcripts"][0]
    ref_base = tx["sequence"][tx["cds_start"] - 1]  # base at c.1

    normalizer = ferro_hgvs.Normalizer(reference_json=str(out))
    normalized = normalizer.normalize(f"tx1:c.1{ref_base}>G")
    assert normalized == "tx1:c.1{}>G".format(ref_base)


def test_emit_genomic_sequences_makes_reference_genome_capable(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)
    fasta = _write_fasta(tmp_path)
    out = tmp_path / "genome_capable.json"

    cfg = ferro_hgvs.ConvertGffConfig(
        gff=str(gff),
        fasta=str(fasta),
        output=str(out),
        validate_fasta=False,
        emit_genomic_sequences=True,
    )
    report = ferro_hgvs.convert_gff(cfg)
    assert report.emitted_genomic_bytes == 520

    doc = json.loads(out.read_text())
    assert "genomic_sequences" in doc
    assert doc["genomic_sequences"]["chr1"] == "".join("ACGT"[i % 4] for i in range(520))

    normalizer = ferro_hgvs.Normalizer(reference_json=str(out))
    assert normalizer.has_genomic_data() is True
    assert normalizer.reference_summary()["has_genomic_data"] is True


def test_transcript_only_reference_reports_reduced_capability(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)
    out = tmp_path / "transcript_only.json"
    ferro_hgvs.convert_gff(ferro_hgvs.ConvertGffConfig(gff=str(gff), output=str(out)))

    normalizer = ferro_hgvs.Normalizer(reference_json=str(out))
    assert normalizer.has_genomic_data() is False


def test_gene_filter_accepts_string_and_list_equivalently(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)

    as_string = ferro_hgvs.convert_gff(
        ferro_hgvs.ConvertGffConfig(gff=str(gff), genes="GENE1")
    ).transcripts_json
    as_list = ferro_hgvs.convert_gff(
        ferro_hgvs.ConvertGffConfig(gff=str(gff), genes=["GENE1"])
    ).transcripts_json
    assert as_string == as_list
    assert len(json.loads(as_string)["transcripts"]) == 1

    # A non-matching gene filters everything out.
    none_match = ferro_hgvs.convert_gff(
        ferro_hgvs.ConvertGffConfig(gff=str(gff), genes="NOSUCHGENE")
    )
    assert none_match.transcript_count == 0


def test_transcript_filter_selects_by_id(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)
    report = ferro_hgvs.convert_gff(ferro_hgvs.ConvertGffConfig(gff=str(gff), transcripts=["tx1"]))
    assert report.transcript_count == 1
    report_miss = ferro_hgvs.convert_gff(
        ferro_hgvs.ConvertGffConfig(gff=str(gff), transcripts=["tx_absent"])
    )
    assert report_miss.transcript_count == 0


def test_invalid_error_mode_raises_value_error(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)
    with pytest.raises(ValueError, match="error_mode"):
        ferro_hgvs.convert_gff(ferro_hgvs.ConvertGffConfig(gff=str(gff), error_mode="loud"))


def test_emit_genomic_sequences_without_fasta_raises_value_error(tmp_path: Path) -> None:
    gff = _write_gff(tmp_path)
    with pytest.raises(ValueError, match="fasta"):
        ferro_hgvs.convert_gff(
            ferro_hgvs.ConvertGffConfig(gff=str(gff), emit_genomic_sequences=True)
        )


def test_strict_mode_raises_on_malformed_record(tmp_path: Path) -> None:
    bad = tmp_path / "bad.gff3"
    bad.write_text("##gff-version 3\nchr1\t.\tgene\tnot_a_number\t500\t.\t+\t.\tID=g1\n")
    with pytest.raises(ferro_hgvs.ReferenceDataError):
        ferro_hgvs.convert_gff(ferro_hgvs.ConvertGffConfig(gff=str(bad), error_mode="strict"))


def test_emit_with_nothing_placed_warns_and_stays_transcript_only(tmp_path: Path) -> None:
    """emit_genomic_sequences with nothing left to place (every transcript
    filtered out) warns and does not write genomic_sequences — mirroring the
    CLI's stderr warning."""
    gff = _write_gff(tmp_path)
    fasta = _write_fasta(tmp_path)
    out = tmp_path / "out.json"

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        report = ferro_hgvs.convert_gff(
            ferro_hgvs.ConvertGffConfig(
                gff=str(gff),
                fasta=str(fasta),
                output=str(out),
                validate_fasta=False,
                emit_genomic_sequences=True,
                genes="NOSUCHGENE",  # filter everything out -> nothing placed
            )
        )

    assert report.transcript_count == 0
    assert any("no emitted transcript" in w for w in report.warnings)
    assert any(issubclass(c.category, UserWarning) for c in caught)
    doc = json.loads(out.read_text())
    assert "genomic_sequences" not in doc
