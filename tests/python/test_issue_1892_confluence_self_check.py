"""Issue #1892: the opt-in confluence self-check, via the Python API.

``EquivalenceChecker.check_confluence`` groups a corpus of variants into
equivalence classes under a chosen relation and reports every class whose
members normalize to more than one distinct output — a non-confluence witness.
It is a diagnostic: it reports offending groups and never emits a pass/fail
release verdict.

The reference is fully synthetic: the 21-base ``NC_KEY.1`` contig that
``spdi::apply``'s own tests use, so the grouping can be checked by hand.
"""

import json
from pathlib import Path

import pytest

import ferro_hgvs

CONTIG = "NC_KEY.1"
SEQUENCE = "GGATTACAGGCATTAGCCTGA"  # 1-based: GGCATTAGCCT sits at 9..=19


def _checker(tmp_path: Path) -> "ferro_hgvs.EquivalenceChecker":
    ref = {"transcripts": [], "genomic_sequences": {CONTIG: SEQUENCE}}
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(ref))
    return ferro_hgvs.EquivalenceChecker(reference_json=str(path))


def _variants(descriptors):
    return [ferro_hgvs.parse(d) for d in descriptors]


def test_confluent_corpus_reports_no_violations(tmp_path: Path) -> None:
    checker = _checker(tmp_path)
    # A literal insertion and the same bases named by reference range: apply-equal
    # and converged by the normalizer.
    corpus = _variants([f"{CONTIG}:g.5_6insGGCATTAGCCT", f"{CONTIG}:g.5_6ins9_19"])

    report = checker.check_confluence(corpus)

    assert report["relation"] == "cross_axis"
    assert report["is_confluent"] is True
    assert report["is_complete"] is True
    assert report["violations"] == []
    assert report["skipped"] == []
    assert report["classes_checked"] == 1
    assert report["undecided_pairs"] == 0


def test_distinct_variants_form_separate_classes(tmp_path: Path) -> None:
    checker = _checker(tmp_path)
    # Three inputs, two variants: `g.3delinsG` is `g.3A>G` respelled and must
    # join it. Two inputs and two classes would also be what a grouping that
    # never merges produces, so the count would pin nothing.
    corpus = _variants([f"{CONTIG}:g.3A>G", f"{CONTIG}:g.3A>C", f"{CONTIG}:g.3delinsG"])

    report = checker.check_confluence(corpus, relation="cross_axis")

    assert report["is_confluent"] is True
    assert report["is_complete"] is True
    assert report["classes_checked"] == 2


def test_spdi_relation_is_selectable(tmp_path: Path) -> None:
    checker = _checker(tmp_path)
    corpus = _variants([f"{CONTIG}:g.13_14insT", f"{CONTIG}:g.14dup"])

    report = checker.check_confluence(corpus, relation="spdi")

    assert report["relation"] == "spdi"
    assert report["is_confluent"] is True
    assert report["classes_checked"] == 1


def test_inputs_the_relation_cannot_place_are_skipped(tmp_path: Path) -> None:
    checker = _checker(tmp_path)
    # An unspecified inserted length denotes no resolvable sequence, so it has no
    # SPDI key and cannot be placed.
    corpus = _variants([f"{CONTIG}:g.10_11ins(10)", f"{CONTIG}:g.5_6ins(6)"])

    report = checker.check_confluence(corpus, relation="spdi")

    assert report["violations"] == []
    assert len(report["skipped"]) == 2
    assert all(entry["reason"] for entry in report["skipped"])
    assert all("input" in entry for entry in report["skipped"])
    assert all(entry["kind"] == "unplaceable" for entry in report["skipped"])
    assert report["is_complete"] is False


def test_unknown_relation_raises_value_error(tmp_path: Path) -> None:
    checker = _checker(tmp_path)
    with pytest.raises(ValueError):
        checker.check_confluence(_variants([f"{CONTIG}:g.3A>G"]), relation="bogus")


def test_empty_corpus_is_confluent(tmp_path: Path) -> None:
    checker = _checker(tmp_path)
    report = checker.check_confluence([])
    assert report["is_confluent"] is True
    assert report["is_complete"] is True
    assert report["classes_checked"] == 0
    assert report["undecided_pairs"] == 0


def test_a_class_that_normalizes_apart_is_reported(tmp_path: Path) -> None:
    """The headline behaviour: a violation is reported, with its inputs and outputs.

    `g.[3A>G;10=]` states an unchanged base alongside the substitution and
    `g.3A>G` does not; applying either produces the same bases, so they are one
    class, and the normalizer keeps the `=` member, so the class reaches two
    outputs. Every other test in this file asserts `violations == []`, which
    cannot observe the reporting path at all.
    """
    checker = _checker(tmp_path)
    corpus = _variants([f"{CONTIG}:g.[3A>G;10=]", f"{CONTIG}:g.3A>G"])

    report = checker.check_confluence(corpus)

    assert report["classes_checked"] == 1
    assert report["is_complete"] is True, "the finding is about the corpus, not the reference"
    assert report["is_confluent"] is False
    assert len(report["violations"]) == 1
    assert report["violations"][0]["inputs"] == [
        f"{CONTIG}:g.[3A>G;10=]",
        f"{CONTIG}:g.3A>G",
    ]
    assert report["violations"][0]["outputs"] == [
        f"{CONTIG}:g.3A>G",
        f"{CONTIG}:g.[3A>G;10=]",
    ]


def test_an_unreadable_reference_is_named_rather_than_read_as_clean(
    tmp_path: Path,
) -> None:
    """A reference that cannot serve the bases must not read as a clean corpus.

    Both descriptions are well formed and denote one variant; the contig is
    simply not in the reference. Before this was reported, the run answered
    `is_confluent: True`, `skipped: []` — byte-identical to a corpus that was
    read and found confluent.
    """
    checker = _checker(tmp_path)
    corpus = _variants(["NC_ABSENT.1:g.2del", "NC_ABSENT.1:g.2delC"])

    for relation in ("cross_axis", "spdi"):
        report = checker.check_confluence(corpus, relation=relation)

        assert report["is_complete"] is False, relation
        assert report["skipped"], relation
        assert all(e["kind"] == "unplaceable" for e in report["skipped"]), relation
        # The refusal must be classified against the reference rather than
        # blamed on the description. Keyed on `decline_class`, not on the
        # sentence: sniffing `reason` for "reference" passes whenever that word
        # happens to appear and breaks silently on any rewording, so it tests
        # the prose instead of the classification.
        assert any(e["decline_class"] == "reference_unavailable" for e in report["skipped"]), (
            relation
        )
        # The site is machine-readable too, so a consumer never parses English.
        assert all(e["decline_site"] in (None, "member_conversion") for e in report["skipped"]), (
            relation
        )
