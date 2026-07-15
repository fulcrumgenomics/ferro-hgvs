"""Tests for #1012 (sibling scope): reference-capability surfacing on the
classes and free functions that mirror ``Normalizer``.

The identical ``reference_json`` / no-args footgun (silent no-genomic fallback,
no introspection) lives on ``VariantProjector``, ``EquivalenceChecker``,
``BatchProcessor``, ``CoordinateMapper``, the free ``normalize()`` function, and
``HgvsVariant.normalize()``. These tests assert:

  * introspection getters (``has_genomic_data`` / ``has_protein_data`` /
    ``reference_summary``) on each of the four classes;
  * the one-time reduced-capability ``UserWarning`` fires at construction / call
    on each surface (asserted in fresh subprocesses, since the warn-once flag is
    process-global and shared across every surface);
  * the ``provider_kind`` label is ``"json"`` (not ``"transcripts_json"``) for a
    ``reference_json`` build.
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


def _run_child(code: str, *args: str) -> subprocess.CompletedProcess[str]:
    """Run ``code`` in a fresh interpreter.

    The reduced-capability ``UserWarning`` fires at most once per process and the
    flag is shared across every surface, so it must be asserted in a fresh
    interpreter that starts with the warn-once flag unset.
    """
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code), *args],
        capture_output=True,
        text=True,
    )


def _assert_child_ok(code: str, *args: str) -> None:
    """Run ``code`` in a fresh interpreter and assert it exited cleanly."""
    result = _run_child(code, *args)
    assert result.returncode == 0, f"stdout={result.stdout!r} stderr={result.stderr!r}"


def _transcript_ref(tmp_path: Path) -> str:
    """Write a minimal transcript-only (no-genomic) reference; return its path."""
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


# Each entry: (class-name, no-arg constructor callable). Every one falls back to
# built-in test data with no genomic sequences.
_CLASSES = [
    ("VariantProjector", lambda: ferro_hgvs.VariantProjector()),
    ("EquivalenceChecker", lambda: ferro_hgvs.EquivalenceChecker()),
    ("BatchProcessor", lambda: ferro_hgvs.BatchProcessor()),
    ("CoordinateMapper", lambda: ferro_hgvs.CoordinateMapper()),
]


class TestSiblingIntrospection:
    """Capability getters on the sibling classes."""

    @pytest.mark.parametrize("name,ctor", _CLASSES, ids=[c[0] for c in _CLASSES])
    def test_test_data_reports_no_genomic(self, name: str, ctor) -> None:  # noqa: ANN001
        obj = ctor()
        assert obj.has_genomic_data() is False
        assert isinstance(obj.has_protein_data(), bool)
        summary = obj.reference_summary()
        assert summary["provider_kind"] == "test_data"
        assert summary["has_genomic_data"] is False
        assert isinstance(summary["has_protein_data"], bool)

    def test_projector_reference_json_labeled_json(self, tmp_path: Path) -> None:
        proj = ferro_hgvs.VariantProjector(reference_json=_transcript_ref(tmp_path))
        summary = proj.reference_summary()
        assert summary["provider_kind"] == "json"
        assert summary["has_genomic_data"] is False

    def test_checker_reference_json_labeled_json(self, tmp_path: Path) -> None:
        checker = ferro_hgvs.EquivalenceChecker(reference_json=_transcript_ref(tmp_path))
        assert checker.reference_summary()["provider_kind"] == "json"

    def test_batch_reference_json_labeled_json(self, tmp_path: Path) -> None:
        batch = ferro_hgvs.BatchProcessor(reference_json=_transcript_ref(tmp_path))
        assert batch.reference_summary()["provider_kind"] == "json"

    def test_mapper_reference_json_labeled_json(self, tmp_path: Path) -> None:
        mapper = ferro_hgvs.CoordinateMapper(reference_json=_transcript_ref(tmp_path))
        assert mapper.reference_summary()["provider_kind"] == "json"

    def test_projector_manifest_reports_kind(self) -> None:
        proj = ferro_hgvs.VariantProjector.from_manifest(str(MANIFEST_TINY))
        summary = proj.reference_summary()
        assert summary["provider_kind"] == "manifest"
        assert isinstance(summary["has_genomic_data"], bool)


class TestSiblingConstructionWarnings:
    """The reduced-capability UserWarning fires once, per surface."""

    @pytest.mark.parametrize(
        "ctor_expr",
        [
            "ferro_hgvs.VariantProjector()",
            "ferro_hgvs.EquivalenceChecker()",
            "ferro_hgvs.BatchProcessor()",
            "ferro_hgvs.CoordinateMapper()",
        ],
    )
    def test_warns_on_test_data_build(self, ctor_expr: str) -> None:
        code = f"""
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            {ctor_expr}
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected one UserWarning, got {{len(user)}}"
        msg = str(user[0].message)
        assert "genomic data" in msg, msg
        assert "from_manifest" in msg, msg
        """
        _assert_child_ok(code)

    def test_shared_global_warns_at_most_once_across_surfaces(self) -> None:
        # The warn-once flag is shared across every surface, so building a
        # projector and then a normalizer must yield exactly one warning.
        code = """
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.VariantProjector()
            ferro_hgvs.EquivalenceChecker()
            ferro_hgvs.Normalizer()
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected exactly one UserWarning, got {len(user)}"
        """
        _assert_child_ok(code)

    def test_projector_reference_json_warns(self, tmp_path: Path) -> None:
        ref = _transcript_ref(tmp_path)
        code = """
        import sys
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.VariantProjector(reference_json=sys.argv[1])
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected one UserWarning, got {len(user)}"
        """
        _assert_child_ok(code, ref)


class TestFreeFunctionWarnings:
    """The free normalize() and HgvsVariant.normalize() also warn."""

    def test_free_normalize_warns(self) -> None:
        code = """
        import warnings
        import ferro_hgvs
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ferro_hgvs.normalize("NM_000088.3:c.100delA")
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected one UserWarning, got {len(user)}"
        msg = str(user[0].message)
        assert "genomic data" in msg, msg
        assert "from_manifest" in msg, msg
        """
        _assert_child_ok(code)

    def test_hgvsvariant_normalize_warns(self) -> None:
        code = """
        import warnings
        import ferro_hgvs
        variant = ferro_hgvs.parse("NM_000088.3:c.100delA")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            variant.normalize()
        user = [w for w in caught if issubclass(w.category, UserWarning)]
        assert len(user) == 1, f"expected one UserWarning, got {len(user)}"
        """
        _assert_child_ok(code)
