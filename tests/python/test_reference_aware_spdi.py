"""Tests for reference-aware SPDI and apply-to-reference (#1159).

#1159 reports that there is no way to get an encoding-invariant variant key from
Python. The module-level ``hgvs_to_spdi`` does no reference lookup, so it fails on
every edit whose bases only the reference knows — ``del``, ``delins``, ``inv``,
``dup``. And the two fallbacks are unavailable: normalized HGVS strings are not
encoding-invariant for complex indels (#1157), and ``EquivalenceChecker`` was
affected too (#1158).

Three methods answer it, and the tests below separate what each one promises:

* ``Normalizer.to_spdi`` — reference-aware, but a *transliteration*: it preserves
  the caller's partitioning, so it is explicitly **not** a key.
* ``Normalizer.canonical_spdi`` — derived from the resulting bases, so equal iff
  two descriptions denote the same edit. This is the key #1159 asks for.
* ``Normalizer.apply_to_reference`` — the resulting bases themselves, the ground
  truth both of the above rest on.
"""

import json
import tempfile
from pathlib import Path

import pytest

import ferro_hgvs

# 1-based:      1234567890
SEQUENCE = "GGATTACAGGCATTAGCCTGAGGATTACAGGCATTAGCCT"


@pytest.fixture(scope="module")
def normalizer():
    """A transcript-only provider whose bases are known, so expected triples can
    be written by hand rather than read back out of the code under test."""
    reference = {
        "transcripts": [
            {
                "id": "TEMPLATE",
                "gene_symbol": "T",
                "strand": "+",
                "sequence": SEQUENCE,
                "exons": [{"number": 1, "start": 1, "end": len(SEQUENCE)}],
            }
        ]
    }
    path = Path(tempfile.mkdtemp()) / "reference.json"
    path.write_text(json.dumps(reference))
    return ferro_hgvs.Normalizer(reference_json=str(path))


# The four edit types #1159 names as failing without a reference.
NEEDS_REFERENCE = [
    "TEMPLATE:n.3_7delinsGGCTA",
    "TEMPLATE:n.3_5del",
    "TEMPLATE:n.3_5inv",
    "TEMPLATE:n.3dup",
]


@pytest.mark.parametrize("descriptor", NEEDS_REFERENCE)
def test_the_module_level_conversion_still_declines(descriptor):
    """The behaviour #1159 reports, pinned so the new methods are measured against
    a real gap rather than a supposed one."""
    with pytest.raises(ferro_hgvs.ProjectionError) as excinfo:
        ferro_hgvs.hgvs_to_spdi(ferro_hgvs.parse(descriptor))
    assert "reference" in str(excinfo.value).lower()


@pytest.mark.parametrize("descriptor", NEEDS_REFERENCE)
def test_reference_aware_conversion_resolves_them(normalizer, descriptor):
    spdi = normalizer.to_spdi(ferro_hgvs.parse(descriptor))
    assert spdi.sequence == "TEMPLATE"


def test_reference_aware_conversion_resolves_the_expected_bases(normalizer):
    """Pinned, not merely non-raising: a conversion that resolved to the wrong
    bases would satisfy the test above."""
    p = ferro_hgvs.parse
    assert str(normalizer.to_spdi(p("TEMPLATE:n.3_7delinsGGCTA"))) == "TEMPLATE:2:ATTAC:GGCTA"
    assert str(normalizer.to_spdi(p("TEMPLATE:n.3_5del"))) == "TEMPLATE:2:ATT:"
    assert str(normalizer.to_spdi(p("TEMPLATE:n.3_5inv"))) == "TEMPLATE:2:ATT:AAT"


def test_canonical_spdi_is_the_same_for_two_encodings(normalizer):
    """#1159's central requirement: a key equal iff two descriptions denote the
    same edit.

    Reference 3..=7 is ``ATTAC``; replacing it with ``GGCTA`` changes all five
    bases, so the decomposition is five members. ``to_spdi`` gives one triple for
    the spanning form and cannot express the allele at all; ``canonical_spdi``
    gives one answer for both.
    """
    p = ferro_hgvs.parse
    spanning = normalizer.canonical_spdi(p("TEMPLATE:n.3_7delinsGGCTA"))
    decomposed = normalizer.canonical_spdi(p("TEMPLATE:n.[3A>G;4T>G;5T>C;6A>T;7C>A]"))
    assert str(spanning) == str(decomposed)
    assert str(spanning) == "TEMPLATE:2:ATTAC:GGCTA"


def test_canonical_spdi_ignores_member_order(normalizer):
    """A cis allele is a set of edits on one molecule, so order carries no
    meaning and must not change the key."""
    p = ferro_hgvs.parse
    assert str(normalizer.canonical_spdi(p("TEMPLATE:n.[3A>G;7C>A]"))) == str(
        normalizer.canonical_spdi(p("TEMPLATE:n.[7C>A;3A>G]"))
    )


def test_different_edits_get_different_keys(normalizer):
    """A key that collapsed everything would satisfy every test above."""
    p = ferro_hgvs.parse
    keys = {
        str(normalizer.canonical_spdi(p(d)))
        for d in ["TEMPLATE:n.3A>G", "TEMPLATE:n.3A>C", "TEMPLATE:n.4T>G", "TEMPLATE:n.3_5del"]
    }
    assert len(keys) == 4


def test_apply_to_reference_returns_both_windows(normalizer):
    applied = normalizer.apply_to_reference(ferro_hgvs.parse("TEMPLATE:n.3_7delinsGGCTA"))
    assert applied.accession == "TEMPLATE"
    assert applied.start == 2  # 0-based interbase start of n.3
    assert applied.reference == "ATTAC"
    assert applied.resulting == "GGCTA"
    assert "AppliedVariant" in repr(applied)


def test_apply_to_reference_agrees_with_the_sequence_itself(normalizer):
    """The window really is a slice of the reference — otherwise `resulting` would
    be a claim about bases nobody checked."""
    applied = normalizer.apply_to_reference(ferro_hgvs.parse("TEMPLATE:n.3_5del"))
    start, end = applied.start, applied.start + len(applied.reference)
    assert applied.reference == SEQUENCE[start:end]
    assert applied.resulting == ""


@pytest.mark.parametrize(
    "descriptor",
    [
        # Overlapping members: applying them depends on order, so there is no one
        # resulting sequence.
        "TEMPLATE:n.[3_5del;4T>G]",
        # Two molecules, not one sequence.
        "TEMPLATE:n.[3A>G(;)7C>A]",
    ],
)
def test_shapes_without_one_resulting_sequence_are_refused(normalizer, descriptor):
    """Refused, not answered from one arbitrary reading of an ambiguous input."""
    variant = ferro_hgvs.parse(descriptor)
    with pytest.raises(ferro_hgvs.ProjectionError):
        normalizer.canonical_spdi(variant)
    with pytest.raises(ferro_hgvs.ProjectionError):
        normalizer.apply_to_reference(variant)


def test_an_unknown_accession_is_refused(normalizer):
    with pytest.raises(ferro_hgvs.ProjectionError):
        normalizer.canonical_spdi(ferro_hgvs.parse("ABSENT:n.3A>G"))


def test_applied_variant_is_exported():
    """``AppliedVariant`` is a new public class the stubs declare and
    ``apply_to_reference`` returns, so it belongs in ``__all__`` beside
    ``SpdiVariant``.

    ``from .ferro_hgvs import *`` makes it reachable as an attribute either way,
    but ``__all__`` is what governs ``from ferro_hgvs import *`` and what tools
    read as the package's public API.
    """
    assert "AppliedVariant" in ferro_hgvs.__all__
    assert "SpdiVariant" in ferro_hgvs.__all__
