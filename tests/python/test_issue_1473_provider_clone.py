"""Per-call cost must not scale with reference size (#1473).

Every reference-backed binding builds a fresh normalizer per call --
``Normalizer::with_config(self.provider.clone(), ...)`` -- so whatever the
provider handle holds is copied once per Python call. The ``MultiFasta`` variant
was already an ``Arc`` for exactly that reason; the ``JsonProvider`` variant was a
``Box``, chosen to keep the enum small, and ``JsonProvider`` derives ``Clone``
with ``String``-owning maps. So a ``Normalizer(reference_json=...)`` deep-copied
its entire reference on every call.

Measured before the fix, on a variant whose own work is held constant -- a 6-base
deletion on a 600-base transcript at the start of the contig -- while only the
*unread* remainder of the contig grew:

============  ===================  =========
contig bases   ``normalize_variant``  ``to_spdi``
============  ===================  =========
      10 000               1.2 us     1.0 us
     100 000               6.8 us     6.1 us
   1 000 000              21.8 us    17.5 us
   4 000 000             443.2 us    74.0 us
============  ===================  =========

~9 GB/s, i.e. memcpy speed, for bases no call reads. After the fix the same table
is flat at 0.8 us / 0.5 us.

The assertion below is on the **shape** rather than on a microsecond figure: cost
must not grow with reference size. That is what the defect was, it survives a
change of machine or build profile, and it does not turn a slow CI runner red.
The signal is enormous -- a ~180x ratio before, ~1x after -- so the bound is set
far from both rather than between them.
"""

import json
import tempfile
import time
from pathlib import Path

import pytest

import ferro_hgvs

# The transcript is identical in both references and sits at the front of the
# contig, so every call does the same work; only the tail the call never reads
# differs. 2 MB is enough to make a deep copy unmistakable while keeping the JSON
# fixture cheap to write and parse.
SMALL_CONTIG = 10_000
LARGE_CONTIG = 2_000_000

TRANSCRIPT_LENGTH = 600
DESCRIPTOR = "NM_TEST.1:n.100_105del"

# Pre-fix this ratio was ~180 and post-fix it is ~1. Anything under 5 means cost
# is not tracking reference size, which is the property; anything at 180 means
# the clone is back.
MAX_COST_RATIO = 5.0

# Enough iterations that the per-call figure is not dominated by timer
# resolution, and few enough that the pre-fix cost (~220 us a call at 2 MB) does
# not make a failing run take minutes.
CALLS = 300
WARMUP = 30


def _build_normalizer(contig_length):
    contig = "ACGT" * (contig_length // 4)
    reference = {
        "transcripts": [
            {
                "id": "NM_TEST.1",
                "gene_symbol": "T",
                "strand": "+",
                "sequence": contig[2 : 2 + TRANSCRIPT_LENGTH],
                "chromosome": "chr1",
                "genomic_start": 3,
                "genomic_end": 2 + TRANSCRIPT_LENGTH,
                "exons": [
                    {
                        "number": 1,
                        "start": 1,
                        "end": TRANSCRIPT_LENGTH,
                        "genomic_start": 3,
                        "genomic_end": 2 + TRANSCRIPT_LENGTH,
                    }
                ],
            }
        ],
        "genomic_sequences": {"chr1": contig},
    }
    path = Path(tempfile.mkdtemp()) / "reference.json"
    path.write_text(json.dumps(reference))
    return ferro_hgvs.Normalizer(reference_json=str(path))


@pytest.fixture(scope="module")
def small_normalizer():
    return _build_normalizer(SMALL_CONTIG)


@pytest.fixture(scope="module")
def large_normalizer():
    return _build_normalizer(LARGE_CONTIG)


@pytest.fixture(scope="module")
def variant():
    return ferro_hgvs.parse(DESCRIPTOR)


def _seconds_per_call(call):
    for _ in range(WARMUP):
        call()
    # Minimum of three: contention only ever adds time, so the smallest sample is
    # the least contaminated one.
    return min(_timed(call) for _ in range(3))


def _timed(call):
    started = time.perf_counter()
    for _ in range(CALLS):
        call()
    return (time.perf_counter() - started) / CALLS


@pytest.mark.parametrize(
    "entry_point",
    ["normalize_variant", "to_spdi", "canonical_spdi", "apply_to_reference"],
)
def test_cost_does_not_scale_with_unread_reference(
    small_normalizer, large_normalizer, variant, entry_point
):
    small = _seconds_per_call(lambda: getattr(small_normalizer, entry_point)(variant))
    large = _seconds_per_call(lambda: getattr(large_normalizer, entry_point)(variant))

    ratio = large / small
    assert ratio < MAX_COST_RATIO, (
        f"Normalizer.{entry_point} cost {large * 1e6:.1f} us against a "
        f"{LARGE_CONTIG:,}-base reference and {small * 1e6:.1f} us against a "
        f"{SMALL_CONTIG:,}-base one — a {ratio:.1f}x difference for identical "
        "work on an identical transcript. The binding is copying the whole "
        "reference per call (#1473)."
    )


def test_both_references_give_the_same_answer(small_normalizer, large_normalizer, variant):
    """The two fixtures must differ only in the tail no call reads — otherwise the
    comparison above would be timing two different pieces of work."""
    assert str(small_normalizer.to_spdi(variant)) == str(large_normalizer.to_spdi(variant))
