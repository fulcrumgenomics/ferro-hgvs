"""Issue #1244: the overlapping-allele crash, seen from the Python bindings.

``EquivalenceChecker.check`` panicked (slice index out of bounds) for a cis
allele whose members overlap. The Rust suite
(``tests/it/issue_1244_equivalence_overlap_panic.rs``) covers the fix at the
library level; this file exists because the Python surface is where the panic
did the most damage.

A Rust panic crosses PyO3 as ``pyo3_runtime.PanicException``, which subclasses
``BaseException`` rather than ``Exception`` — so it slips straight through the
``except Exception`` that a caller would reasonably wrap ``check`` in, and takes
down the interpreter's call stack instead of yielding a verdict. Returning
``NotEquivalent`` is not merely a nicer answer here; it is the difference
between a catchable outcome and an uncatchable one.

The reference mirrors the Rust fixture exactly (``SyntheticBuilder::genomic``):
256 bases of ``ACGT`` padding on each side of the core, so core position ``n``
sits at 1-based HGVS ``256 + n`` and the normalizer's lookahead window stays in
bounds.
"""

import json
from pathlib import Path

import ferro_hgvs

PAD = "ACGT" * 64  # 256 bp, matching SyntheticBuilder's PAD_OFFSET
CORE = "CAGCTAGCTGATATATATATGCGCGCGCGC"  # core 11..=20 is an (AT) tandem repeat
CONTIG = "NC_TEST.1"

# A repeat expansion over core 11..=14 followed by a deletion over core 14..=20.
# Both members claim core position 14, so the allele's members overlap and it has
# no single resulting sequence. This is the input from the issue.
OVERLAPPING_ALLELE = "NC_TEST.1:g.[267_270AT[4];270_276del]"

# A substitution at core position 1 — a different locus entirely.
ELSEWHERE = "NC_TEST.1:g.257C>A"


def _checker(tmp_path: Path) -> "ferro_hgvs.EquivalenceChecker":
    """Build a checker over the synthetic reference, written under ``tmp_path``."""
    reference = {"transcripts": [], "genomic_sequences": {CONTIG: PAD + CORE + PAD}}
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(reference))
    return ferro_hgvs.EquivalenceChecker(reference_json=str(path))


def test_overlapping_cis_allele_returns_a_verdict(tmp_path):
    """The reported crash: a verdict comes back instead of a ``PanicException``.

    Reaching the assertions at all is half the test — a panic would abort the
    call rather than return.
    """
    checker = _checker(tmp_path)
    result = checker.check(ferro_hgvs.parse(OVERLAPPING_ALLELE), ferro_hgvs.parse(ELSEWHERE))

    # The two describe different loci, so NotEquivalent is the honest answer.
    assert result.level == ferro_hgvs.EquivalenceLevel.NotEquivalent
    assert not result.level.is_equivalent()


def test_overlapping_cis_allele_verdict_is_symmetric(tmp_path):
    """Argument order must not decide whether the splice is reached."""
    checker = _checker(tmp_path)
    result = checker.check(ferro_hgvs.parse(ELSEWHERE), ferro_hgvs.parse(OVERLAPPING_ALLELE))

    assert result.level == ferro_hgvs.EquivalenceLevel.NotEquivalent


def test_all_equivalent_also_returns_a_verdict(tmp_path):
    """The sibling surface: ``all_equivalent`` is the other door into the crash.

    It is a thin loop over the same comparison, so it raised the same
    ``PanicException`` before the fix — and unlike ``check`` it had no Python
    coverage at all, so nothing would have caught a regression that reopened
    only this door.
    """
    checker = _checker(tmp_path)
    variants = [ferro_hgvs.parse(OVERLAPPING_ALLELE), ferro_hgvs.parse(ELSEWHERE)]

    # Different loci, so the honest answer is False — the point is that one
    # comes back rather than the call aborting.
    assert checker.all_equivalent(variants) is False
