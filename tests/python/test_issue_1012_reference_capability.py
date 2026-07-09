"""Tests for #1012: reference-capability introspection and the reduced-capability warning.

Covers:
  * Item 3 — ``Normalizer.has_genomic_data`` / ``has_protein_data`` /
    ``reference_summary`` introspection.
  * Item 1 — the one-time ``UserWarning`` emitted when a reduced-capability
    (no-genomic) reference is loaded.
  * Item 5 — a wholly-empty ``reference_json`` is rejected at load.
"""

from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

import ferro_hgvs

FIXTURES = Path(__file__).parent.parent / "fixtures"
MANIFEST_TINY = FIXTURES / "python" / "manifest_tiny" / "manifest.json"


def _run_child(code: str) -> subprocess.CompletedProcess[str]:
    """Run ``code`` in a fresh interpreter.

    The reduced-capability ``UserWarning`` fires at most once per process, so it
    cannot be asserted reliably in-process alongside other tests that also build
    reduced-capability normalizers. A fresh interpreter starts with the
    warn-once flag unset.
    """
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
    )


def _transcript_ref(tmp_path: Path) -> str:
    """Write a minimal transcript-only reference and return its path."""
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
        ]
    }
    path = tmp_path / "ref.json"
    path.write_text(json.dumps(ref))
    return str(path)


class TestReferenceSummary:
    """Item 3: capability introspection getters."""

    def test_test_data_reports_no_genomic(self) -> None:
        norm = ferro_hgvs.Normalizer()
        assert norm.has_genomic_data() is False
        summary = norm.reference_summary()
        assert summary["provider_kind"] == "test_data"
        assert summary["has_genomic_data"] is False
        assert isinstance(summary["has_protein_data"], bool)

    def test_transcripts_json_reports_kind(self, tmp_path: Path) -> None:
        norm = ferro_hgvs.Normalizer(reference_json=_transcript_ref(tmp_path))
        summary = norm.reference_summary()
        # A reference_json build is labeled "json": such a file may carry
        # genomic-only data, so has_genomic_data (not the label) is the source
        # of truth for capability.
        assert summary["provider_kind"] == "json"
        assert summary["has_genomic_data"] is False

    def test_manifest_reports_kind(self) -> None:
        norm = ferro_hgvs.Normalizer.from_manifest(str(MANIFEST_TINY))
        summary = norm.reference_summary()
        assert summary["provider_kind"] == "manifest"
        assert isinstance(norm.has_genomic_data(), bool)
        assert isinstance(norm.has_protein_data(), bool)

    def test_protein_bearing_reference_reports_has_protein_data(self, tmp_path: Path) -> None:
        # A reference that carries protein sequences reports has_protein_data() True
        # (the positive case — the other tests only cover the False / bool shape).
        ref = {
            "transcripts": [],
            "proteins": {"NP_000001.1": "MAPLE"},
        }
        ref_path = tmp_path / "protein_ref.json"
        ref_path.write_text(json.dumps(ref))
        norm = ferro_hgvs.Normalizer(reference_json=str(ref_path))
        assert norm.has_protein_data() is True
        assert norm.reference_summary()["has_protein_data"] is True


class TestReducedCapabilityWarning:
    """Item 1: one-time UserWarning on a no-genomic build."""

    def test_warns_on_test_data_build(self) -> None:
        code = """
        import pytest
        import ferro_hgvs
        with pytest.warns(UserWarning):
            ferro_hgvs.Normalizer()
        """
        result = _run_child(code)
        assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"

    def test_warning_message_points_at_from_manifest(self) -> None:
        code = """
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.Normalizer()
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected one UserWarning, got {len(user)}"
        msg = str(user[0].message)
        assert "genomic data" in msg, msg
        assert "from_manifest" in msg, msg
        """
        result = _run_child(code)
        assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"

    def test_warns_at_most_once_per_process(self) -> None:
        code = """
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.Normalizer()
            ferro_hgvs.Normalizer()  # a second reduced-capability build
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected exactly one UserWarning, got {len(user)}"
        """
        result = _run_child(code)
        assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"

    def test_manifest_warning_not_masked_by_earlier_generic_warning(self) -> None:
        # The generic (test-data / reference_json) and manifest reduced-capability
        # warnings carry different guidance and are keyed on SEPARATE once-flags:
        # a throwaway test-data Normalizer warning first must NOT suppress a
        # genuinely genome-less manifest warning later in the same process — a
        # real production misconfiguration that would otherwise go unwarned.
        code = f"""
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.Normalizer()                                     # generic variant warns
            ferro_hgvs.Normalizer.from_manifest({str(MANIFEST_TINY)!r})  # manifest variant must ALSO warn
        msgs = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        generic = [m for m in msgs if "reference lacks genomic data" in m]
        manifest = [m for m in msgs if "manifest reference provides no genomic data" in m]
        assert len(generic) == 1, f"expected one generic warning, got {{msgs}}"
        assert len(manifest) == 1, f"expected one manifest warning, got {{msgs}}"
        """
        result = _run_child(code)
        assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"

    def test_failed_warning_emission_does_not_consume_the_once_flag(self) -> None:
        # Under warnings.simplefilter("error") the reduced-capability warning is
        # raised as an exception. That must NOT consume the once-per-process flag:
        # a later construction under a normal filter must still warn. (Regression:
        # the flag was set before emitting, so a raised warning permanently
        # suppressed all future ones.)
        code = """
        import warnings
        import ferro_hgvs
        # First build with warnings-as-errors: construction raises, flag must reset.
        raised = False
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            try:
                ferro_hgvs.Normalizer()
            except UserWarning:
                raised = True
        assert raised, "expected the reduced-capability warning to raise under simplefilter('error')"
        # Now a normal build must still emit the warning (flag was not consumed).
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.Normalizer()
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected the warning to still fire, got {len(user)}"
        """
        result = _run_child(code)
        assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"

    def test_genomic_reference_does_not_warn(self, tmp_path: Path) -> None:
        # A reference that DOES carry genomic data has full capability, so it
        # must not trigger the reduced-capability warning.
        ref = {
            "transcripts": [],
            "genomic_sequences": {"NC_000001.11": "ACGTACGTACGTACGT"},
        }
        ref_path = tmp_path / "gref.json"
        ref_path.write_text(json.dumps(ref))
        # Pass the reference path via argv so the child body needs no
        # interpolation (avoids nesting f-strings).
        code = """
        import sys
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            norm = ferro_hgvs.Normalizer(reference_json=sys.argv[1])
        assert norm.has_genomic_data() is True
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 0, f"expected no UserWarning, got {len(user)}"
        """
        result = subprocess.run(
            [sys.executable, "-c", textwrap.dedent(code), str(ref_path)],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"


class TestEmptyReferenceRejected:
    """Item 5: a wholly-empty reference_json is rejected at load."""

    def test_wholly_empty_reference_json_raises(self, tmp_path: Path) -> None:
        empty = tmp_path / "empty.json"
        empty.write_text("{}")
        with pytest.raises(RuntimeError):
            ferro_hgvs.Normalizer(reference_json=str(empty))
