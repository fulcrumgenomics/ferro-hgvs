"""Issue #395 item 5 — `NormalizationWarning` exposed to Python.

`src/python.rs` previously had zero references to
`crate::normalize::NormalizationWarning`. Python callers calling
`Normalizer.normalize(...)` got back only the warning-stripped String
representation. The new `Normalizer.normalize_with_warnings(...)` method
surfaces the warnings on a `NormalizeResultWithWarnings` object.

These tests run under `pytest tests/python/` after `maturin develop
--features python`. They cover:

  - The new `NormalizeResultWithWarnings` shape (result + warnings).
  - `NormalizationWarning.code` and `.message` accessors.
  - The empty-warnings case (no-corrections normalize call).
"""

import json

import ferro_hgvs


def _make_reference_json(tmp_path, contig: str, start_1based: int, bases: str) -> str:
    """Write a JSON reference with `bases` at `start_1based..` and return its path."""
    seq = ["A"] * max(2000, start_1based + len(bases) + 200)
    for i, b in enumerate(bases):
        seq[start_1based - 1 + i] = b
    payload = {
        "transcripts": [],
        "proteins": {},
        "genomic_sequences": {contig: "".join(seq)},
    }
    path = tmp_path / "ref.json"
    path.write_text(json.dumps(payload))
    return str(path)


class TestNormalizeWithWarnings:
    """Smoke tests for the new `normalize_with_warnings` Python method."""

    def test_normalize_with_warnings_returns_result_and_warnings(self) -> None:
        """A normalize call that emits no warnings still returns the
        result wrapped in `NormalizeResultWithWarnings` with an empty
        warnings list.

        Uses an in-bounds, reference-matching position on the bundled
        `NM_000088.3` test transcript (CDS sequence
        `ATGCCCAAGGTGCTGCCCCAGATGCTGCCAGTGCTGCTGCTGCTGCTGCTGCTGCTGCTG`,
        CDS length 60). `c.4C>G` matches the reference base at position 4
        and stays within the CDS, so neither `RefSeqMismatch` nor
        `PositionPastEnd` is emitted."""
        normalizer = ferro_hgvs.Normalizer()
        result = normalizer.normalize_with_warnings("NM_000088.3:c.4C>G")
        assert result.has_warnings() is False
        assert result.warnings == []
        # The wrapped variant should be accessible via `.result`.
        assert str(result.result) == "NM_000088.3:c.4C>G"

    def test_normalize_with_warnings_surfaces_overlap_conflict(self, tmp_path) -> None:
        """Two cis-allele substitutions at coincident bounds emit
        `OVERLAP_CONFLICTING_EDITS` (W5002) in lenient mode."""
        ref = _make_reference_json(tmp_path, "NC_000001.11", 100, "A")
        normalizer = ferro_hgvs.Normalizer(reference_json=ref)
        result = normalizer.normalize_with_warnings("NC_000001.11:g.[100A>C;100A>G]")
        codes = [w.code for w in result.warnings]
        assert "OVERLAP_CONFLICTING_EDITS" in codes, (
            f"expected OVERLAP_CONFLICTING_EDITS in warnings; got {codes}"
        )
        assert result.has_warnings() is True

    def test_normalization_warning_has_code_and_message(self, tmp_path) -> None:
        """Each `NormalizationWarning` carries `.code` and `.message`."""
        ref = _make_reference_json(tmp_path, "NC_000001.11", 100, "A")
        normalizer = ferro_hgvs.Normalizer(reference_json=ref)
        result = normalizer.normalize_with_warnings("NC_000001.11:g.[100A>C;100A>G]")
        assert len(result.warnings) > 0
        w = result.warnings[0]
        # Both accessors must be non-empty.
        assert isinstance(w.code, str) and len(w.code) > 0
        assert isinstance(w.message, str) and len(w.message) > 0
        # __repr__ surfaces the code + message.
        repr_str = repr(w)
        assert w.code in repr_str
