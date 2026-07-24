"""Issue #1158: end-to-end repro via the Python API.

``EquivalenceChecker.check`` returned ``NotEquivalent`` for two variants with
the same resulting sequence when they used different (but equivalent) HGVS
encodings of a complex indel. A sequence-level rung now reconstructs each
variant's edited sequence against the reference and reports ``SequenceMatch``.

The reference here is fully synthetic (no real accession, coordinates, or
reference bases): the core ``AGTCAGT`` sits at 1-based g.21..=27, padded on both
sides so the normalizer's lookahead window stays in bounds.
"""

import json
from pathlib import Path

import ferro_hgvs

PAD = "GGGGCCCCTTTTAAAACGCG"  # 20 bp of padding on each side
CORE = "AGTCAGT"  # g.21..=27 : 21=A 22=G 23=T 24=C 25=A 26=G 27=T
CONTIG = "TEMPLATE"


def _reference(tmp_path: Path) -> str:
    """Write the synthetic reference under pytest's per-test ``tmp_path``.

    ``tmp_path`` (rather than ``tempfile.mktemp``, which is deprecated as
    insecure and leaves a file behind on every run) so pytest owns creation and
    cleanup.
    """
    seq = PAD + CORE + PAD
    ref = {"transcripts": [], "genomic_sequences": {CONTIG: seq}}
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(ref))
    return str(path)


def _checker(tmp_path: Path) -> "ferro_hgvs.EquivalenceChecker":
    return ferro_hgvs.EquivalenceChecker(reference_json=_reference(tmp_path))


def _level(checker, a: str, b: str) -> "ferro_hgvs.EquivalenceLevel":
    return checker.check(ferro_hgvs.parse(a), ferro_hgvs.parse(b)).level


# A single length-changing delins over g.21..=27 replacing AGTCAGT with GATTA,
# and its decomposition as a cis allele of the same edit. Same resulting
# sequence, different normalized strings.
DELINS = "TEMPLATE:g.21_27delinsGATTA"
ALLELE = "TEMPLATE:g.[21A>G;22G>A;24C>T;26_27del]"


def test_delins_vs_decomposed_allele_is_sequence_match(tmp_path):
    checker = _checker(tmp_path)
    result = checker.check(ferro_hgvs.parse(DELINS), ferro_hgvs.parse(ALLELE))
    assert result.level == ferro_hgvs.EquivalenceLevel.SequenceMatch
    assert result.level.is_equivalent()


def test_sanity_identical_and_normalized_still_work(tmp_path):
    checker = _checker(tmp_path)
    assert _level(checker, "TEMPLATE:g.24C>T", "TEMPLATE:g.24C>T") == (
        ferro_hgvs.EquivalenceLevel.Identical
    )
    # substitution written as a 1-base delins normalizes to the same string
    assert _level(checker, "TEMPLATE:g.21A>G", "TEMPLATE:g.21delinsG") == (
        ferro_hgvs.EquivalenceLevel.NormalizedMatch
    )


def test_different_resulting_sequence_stays_not_equivalent(tmp_path):
    checker = _checker(tmp_path)
    assert _level(checker, DELINS, "TEMPLATE:g.21_27delinsGATTC") == (
        ferro_hgvs.EquivalenceLevel.NotEquivalent
    )
