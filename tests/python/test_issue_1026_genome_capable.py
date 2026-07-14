"""Issue #1026: a genome-capable transcripts.json (the shape emitted by
`ferro convert-gff --emit-genomic-sequences`) drives genome-aware normalization
through `Normalizer(reference_json=...)` and reports full capability with no
reduced-capability warning.

This exercises the exact Python surface an integrator uses.
"""

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import ferro_hgvs


def _genome_capable_reference() -> dict:
    """A transcript carrying a placement (chromosome + exon genomic coords) plus
    the backing contig bytes — the object form `--emit-genomic-sequences` writes."""
    return {
        "version": "1.0",
        "genome_build": "GRCh38",
        "transcripts": [
            {
                "id": "NM_1026.1",
                "strand": "+",
                "sequence": "ACGTACGT",
                "chromosome": "chr1",
                "genomic_start": 1,
                "genomic_end": 8,
                "exons": [
                    {
                        "number": 1,
                        "start": 1,
                        "end": 8,
                        "genomic_start": 1,
                        "genomic_end": 8,
                    }
                ],
            }
        ],
        "genomic_sequences": {"chr1": "ACGTACGTACGTACGT"},
    }


def test_genome_capable_reference_reports_full_capability(tmp_path: Path) -> None:
    ref_path = tmp_path / "transcripts.json"
    ref_path.write_text(json.dumps(_genome_capable_reference()))

    n = ferro_hgvs.Normalizer(reference_json=str(ref_path))
    assert n.has_genomic_data() is True
    summary = n.reference_summary()
    assert summary["provider_kind"] == "json"
    assert summary["has_genomic_data"] is True


def test_genome_capable_reference_emits_no_reduced_capability_warning(tmp_path: Path) -> None:
    # The once-per-process warning must be checked in a fresh interpreter.
    ref_path = tmp_path / "transcripts.json"
    ref_path.write_text(json.dumps(_genome_capable_reference()))

    code = """
    import sys
    import warnings
    import ferro_hgvs
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        n = ferro_hgvs.Normalizer(reference_json=sys.argv[1])
    assert n.has_genomic_data() is True
    user = [w for w in caught if issubclass(w.category, UserWarning)]
    assert len(user) == 0, f"expected no UserWarning for a genome-capable reference, got {len(user)}"
    """
    result = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code), str(ref_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"


def test_transcript_only_reference_reports_reduced_capability(tmp_path: Path) -> None:
    ref = _genome_capable_reference()
    del ref["genomic_sequences"]  # placement without bytes → transcript-only
    ref_path = tmp_path / "transcripts_only.json"
    ref_path.write_text(json.dumps(ref))

    n = ferro_hgvs.Normalizer(reference_json=str(ref_path))
    assert n.has_genomic_data() is False
    assert n.reference_summary()["has_genomic_data"] is False


def test_transcript_only_reference_emits_reduced_capability_warning(tmp_path: Path) -> None:
    # The reduced-capability UserWarning must actually FIRE for a transcript-only
    # reference. Checked in a fresh interpreter because the warning is once-per-
    # process and shared across surfaces (an in-process check would be order-flaky).
    ref = _genome_capable_reference()
    del ref["genomic_sequences"]
    ref_path = tmp_path / "transcripts_only.json"
    ref_path.write_text(json.dumps(ref))

    code = """
    import sys
    import warnings
    import ferro_hgvs
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        n = ferro_hgvs.Normalizer(reference_json=sys.argv[1])
    assert n.has_genomic_data() is False
    user = [w for w in caught if issubclass(w.category, UserWarning)]
    assert len(user) == 1, f"expected exactly one reduced-capability UserWarning, got {len(user)}"
    """
    result = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code), str(ref_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"
