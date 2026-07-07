"""Issue #968: the Python provider must forward get_sequence_length (and
infer_genome_build) to its inner provider.

`PyProvider` delegated most `ReferenceProvider` methods to its inner
`Mock`/`MultiFasta` provider but omitted `get_sequence_length`, so it fell
through to the trait default — which returns `Err(ReferenceNotFound)`
unconditionally. Every length-dependent path on the Python side therefore
silently degraded even when the wrapped provider knew the length.

This test exercises the mitochondrial position-past-end bounds check (#393,
W4004 `POSITION_PAST_END`): in the default (lenient) mode, an `m.` position past
the contig length emits a warning — but only if `get_sequence_length` returns the
real contig length. Without the forwarder the length lookup errors and the
warning is silently dropped, so the test fails, which makes it a non-vacuous
guard on the forwarding.

(`infer_genome_build` is forwarded in the same change; it is not covered here
because it is only observable through a `MultiFastaProvider` manifest carrying an
assembly-report contig-alias table plus a projection whose build inference
depends on that alias — see the PR discussion.)
"""

from __future__ import annotations

import json
from pathlib import Path

import ferro_hgvs

MT_ACCESSION = "NC_012920.1"
CONTIG_LENGTH = 100  # short synthetic mtDNA; the real 16569 bp is not needed.


def _reference_json(tmp_path: Path) -> str:
    """Write a MockProvider reference JSON with a short mt contig; return its path."""
    payload = {
        "transcripts": [],
        "proteins": {},
        "genomic_sequences": {MT_ACCESSION: "A" * CONTIG_LENGTH},
    }
    path = tmp_path / "ref.json"
    path.write_text(json.dumps(payload))
    return str(path)


def test_get_sequence_length_forwarded_mt_past_end_warns(tmp_path: Path) -> None:
    """An `m.` position past the contig end warns (W4004) — which requires the
    forwarded `get_sequence_length`. Without the forwarder the length lookup
    errors and the past-end check is silently skipped."""
    norm = ferro_hgvs.Normalizer(reference_json=_reference_json(tmp_path))
    # Position 150 > contig length 100 → past the mt contig 3' end.
    result = norm.normalize_with_warnings(f"{MT_ACCESSION}:m.150A>G")
    codes = [w.code for w in result.warnings]
    assert "POSITION_PAST_END" in codes, (
        f"expected a POSITION_PAST_END (W4004) warning once the contig length is "
        f"available via the forwarded get_sequence_length; got {codes}"
    )


def test_in_bounds_mt_position_does_not_warn(tmp_path: Path) -> None:
    """A control: an in-bounds `m.` position emits no past-end warning, so the
    warning above is driven by the position exceeding the (forwarded) length, not
    emitted unconditionally."""
    norm = ferro_hgvs.Normalizer(reference_json=_reference_json(tmp_path))
    result = norm.normalize_with_warnings(f"{MT_ACCESSION}:m.50A>G")
    codes = [w.code for w in result.warnings]
    assert "POSITION_PAST_END" not in codes, (
        f"an in-bounds m. position must not emit POSITION_PAST_END; got {codes}"
    )
