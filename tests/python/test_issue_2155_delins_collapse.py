"""Issue #2155: a genomic payload-coincidence change collapses to one spanning
``delins`` on both entry points — ``Normalizer.normalize`` and
``ferro_hgvs.from_sequences`` — not only on the coding (``c.``) axis.

The Rust half of this widening is pinned in
``tests/it/issue_2155_payload_coincidence_all_dna.rs`` (the ``normalize``
surface, every DNA axis) and ``tests/it/issue_2155_from_sequences_collapse.rs``
(the ``from_sequences`` surface, ``g.``/``m.``, plus their byte-for-byte
convergence). Neither exercises the Python bindings, and the two surfaces reach
the fix through different code paths — ``normalize`` through
``canonicalize_from_sequence``'s coalesce passes, ``from_sequences`` through the
same passes newly wired into ``derive_block_members`` — so a binding-layer
regression in either forwarding path would not show up in the Rust suite at
all. This file is the case from the issue itself, exercised through the
bindings on the genomic (``g.``) axis, which is what both Rust tests use for
their headline assertion.

Block at 1-based 10_17 of the reference below: ``CTTAGTTA -> AAACAAAC`` (equal
length, a payload coincidence per ``DNA/delins.md:44-47``). Reported
individually, that change fragments into
``g.[10_12delinsAA;14_16delinsCAA;17_18insC]``; the carve-out collapses it to
the single spanning ``g.10_17delinsAAACAAAC``.
"""

import json

import ferro_hgvs

#: Same reference and block as the two Rust test files above, so the expected
#: strings here can be checked against those files rather than re-derived.
REFERENCE = "AGAACCCCCCTTAGTTAAGAACAAAAGCAACAATCTTCGTGGTCCTGG"
#: ``REFERENCE`` with the 1-based 10_17 block replaced by ``AAACAAAC``.
ALT = "AGAACCCCCAAACAAACAGAACAAAAGCAACAATCTTCGTGGTCCTGG"
FRAGMENTED = "TEMPLATE:g.[10_12delinsAA;14_16delinsCAA;17_18insC]"
COLLAPSED = "TEMPLATE:g.10_17delinsAAACAAAC"


def _normalizer(tmp_path) -> "ferro_hgvs.Normalizer":
    """A ``Normalizer`` over a genomic-only reference (no transcripts, just
    ``genomic_sequences`` — a ``g.`` normalization needs nothing else). Uses
    pytest's ``tmp_path`` so the reference file is cleaned up automatically."""
    path = tmp_path / "reference.json"
    path.write_text(json.dumps({"genomic_sequences": {"TEMPLATE": REFERENCE}}))
    return ferro_hgvs.Normalizer(reference_json=str(path))


def test_normalize_collapses_a_genomic_payload_coincidence_change(tmp_path) -> None:
    """The ``Normalizer.normalize`` entry point collapses the fragmented input
    to the single spanning ``delins``."""
    assert _normalizer(tmp_path).normalize(FRAGMENTED) == COLLAPSED


def test_from_sequences_collapses_the_same_change() -> None:
    """The ``from_sequences`` entry point: a pure function of the reference and
    alternate bases, no provider at all."""
    derived = ferro_hgvs.from_sequences("TEMPLATE", 1, REFERENCE, ALT)
    assert str(derived) == COLLAPSED


def test_both_entry_points_agree(tmp_path) -> None:
    """The two surfaces reach the fix through different code paths — see the
    module docstring — so their agreement is itself worth pinning, not only
    each one's own answer."""
    from_normalize = _normalizer(tmp_path).normalize(FRAGMENTED)
    from_derivation = str(ferro_hgvs.from_sequences("TEMPLATE", 1, REFERENCE, ALT))

    assert from_normalize == from_derivation == COLLAPSED
