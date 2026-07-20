"""Guard tests for scripts/refresh-mutalyzer-fixtures.py's curation merge.

`build_refresh_payload` is the only thing standing between an upstream refresh
and the loss of the hand-derived rows, the hand-corrected expected values, and
the per-case dispositions recorded on upstream rows. These tests cover the
expected / error / boundary triad with small inline dicts — no network access,
no reference data, no upstream checkout.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "refresh-mutalyzer-fixtures.py"


def _load_refresh_module() -> ModuleType:
    """Import the hyphenated refresh script as a module."""
    spec = importlib.util.spec_from_file_location("refresh_mutalyzer_fixtures", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


refresh = _load_refresh_module()


def _upstream_cases() -> list[dict]:
    """Two upstream rows, one of which upstream has since changed and extended."""
    return [
        {
            "input": "NM_000000.1:c.1A>T",
            "to_test": True,
            "keywords": ["upstream"],
            # Upstream asserts a value we hand-corrected in the existing fixture.
            "normalized": "NM_000000.1:c.1A>T",
            "protein_description": "NP_000000.1:p.(Met1?)",
            # A genuinely new upstream field that must survive the merge.
            "input_spdi": "NM_000000.1:0:A:T",
        },
        {
            "input": "NM_000000.1:c.9del",
            "to_test": True,
            "keywords": ["upstream"],
            "normalized": "NM_000000.1:c.9del",
        },
        {
            # Upstream asserts a normalized form we deliberately removed because
            # the pinned comparator rejects the variant outright.
            "input": "NM_000000.1:g.5del",
            "to_test": True,
            "keywords": ["upstream"],
            "normalized": "NM_000000.1:n.5del",
        },
    ]


def _existing_fixture() -> dict:
    """A curated fixture: corpus curation, dispositions, corrections, hand rows."""
    return {
        "comparator_provenance": {"note": "pinned comparator"},
        "clusters": [{"id": "protein", "title": "Protein consequence"}],
        "cases": [
            {
                "input": "NM_000000.1:c.1A>T",
                "to_test": True,
                "keywords": ["upstream", "ferro: start-codon policy, see #857"],
                "normalized": "NM_000000.1:c.1A>T",
                # Hand-corrected away from upstream.
                "protein_description": "NP_000000.1:p.(Met1Leu)",
                # Per-case disposition recorded on an UPSTREAM row.
                "spec_citation": "substitution.md:12",
                "accepted_divergence": "upstream under-reports the start codon",
            },
            {
                "input": "NM_000000.1:c.9del",
                "to_test": True,
                "keywords": ["upstream"],
                "normalized": "NM_000000.1:c.9del",
            },
            {
                # `normalized` deliberately dropped: the pinned comparator
                # rejects this variant, so there is no normalized form.
                "input": "NM_000000.1:g.5del",
                "to_test": True,
                "keywords": ["upstream"],
                "errors": ["ECOORDINATESYSTEMMISMATCH"],
            },
            {
                "input": "NM_111111.1:c.100_102dup",
                "to_test": True,
                "keywords": [
                    f"{refresh.HAND_ADDED_PREFIX}protein",
                    "Derived from the pinned comparator.",
                ],
                "protein_description": "NP_111111.1:p.(Ala34dup)",
            },
        ],
    }


def test_merge_preserves_hand_rows_dispositions_and_corrections() -> None:
    """Expected case: curation survives and new upstream fields still arrive."""
    payload = refresh.build_refresh_payload(_upstream_cases(), _existing_fixture(), force=True)

    by_input = {c["input"]: c for c in payload["cases"]}
    assert len(payload["cases"]) == 4

    # Corpus-level curation carried forward verbatim.
    assert payload["clusters"] == [{"id": "protein", "title": "Protein consequence"}]
    assert payload["comparator_provenance"] == {"note": "pinned comparator"}

    # The hand-added row survives.
    hand = by_input["NM_111111.1:c.100_102dup"]
    assert refresh._is_hand_added(hand)
    assert hand["protein_description"] == "NP_111111.1:p.(Ala34dup)"

    # Dispositions and the hand-corrected expected value survive on the
    # upstream row, while the new upstream field arrives.
    merged = by_input["NM_000000.1:c.1A>T"]
    assert merged["spec_citation"] == "substitution.md:12"
    assert merged["accepted_divergence"] == "upstream under-reports the start codon"
    assert merged["protein_description"] == "NP_000000.1:p.(Met1Leu)"
    assert merged["input_spdi"] == "NM_000000.1:0:A:T"
    # Curated keyword annotations are appended to upstream's keywords.
    assert merged["keywords"] == ["upstream", "ferro: start-codon policy, see #857"]
    # An uncorrected upstream row is taken as-is.
    assert by_input["NM_000000.1:c.9del"]["normalized"] == "NM_000000.1:c.9del"
    # A deliberately removed expected value stays removed, not resurrected.
    removed = by_input["NM_000000.1:g.5del"]
    assert "normalized" not in removed
    assert removed["errors"] == ["ECOORDINATESYSTEMMISMATCH"]


def test_duplicate_inputs_are_paired_by_occurrence_order() -> None:
    """Boundary: upstream ships some variants twice; curation must not cross over.

    Upstream's occurrences are in the REVERSED order of the curated rows, so this
    only passes if pairing is by occurrence order — a pairing that matched on the
    `normalized` value instead would cross the curation over and fail.
    """
    upstream = [
        {"input": "NG_000000.1:g.26_31del", "keywords": [], "normalized": "b"},
        {"input": "NG_000000.1:g.26_31del", "keywords": [], "normalized": "a"},
    ]
    existing = {
        "cases": [
            {
                "input": "NG_000000.1:g.26_31del",
                "keywords": [],
                "normalized": "a",
                "spec_citation": "first",
            },
            {
                "input": "NG_000000.1:g.26_31del",
                "keywords": [],
                "normalized": "b",
                "spec_citation": "second",
            },
        ]
    }
    cases = refresh.build_refresh_payload(upstream, existing, force=True)["cases"]
    assert [(c["normalized"], c["spec_citation"]) for c in cases] == [
        ("a", "first"),
        ("b", "second"),
    ]


def test_hand_added_marker_is_found_in_any_keyword_position() -> None:
    """Boundary: prepending a keyword must not reclassify a hand row."""
    case = {"keywords": ["cluster:protein", f"{refresh.HAND_ADDED_PREFIX}protein"]}
    assert refresh._is_hand_added(case)
    assert not refresh._is_hand_added({"keywords": ["upstream"]})
    assert not refresh._is_hand_added({})


def test_refuses_when_curation_present_and_not_forced() -> None:
    """Error case: the guard refuses and enumerates what is at stake."""
    with pytest.raises(SystemExit) as excinfo:
        refresh.build_refresh_payload(_upstream_cases(), _existing_fixture(), force=False)
    message = str(excinfo.value)
    assert "refusing to refresh" in message
    assert "1 hand-added case(s)" in message
    assert "1 upstream row(s) carrying per-case dispositions" in message
    assert "2 upstream row(s) whose expected values were hand-corrected" in message


def test_fresh_import_against_empty_fixture_proceeds() -> None:
    """Boundary: a first import has no curation, so no force is needed."""
    payload = refresh.build_refresh_payload(_upstream_cases(), {}, force=False)
    assert [c["input"] for c in payload["cases"]] == [
        "NM_000000.1:c.1A>T",
        "NM_000000.1:c.9del",
        "NM_000000.1:g.5del",
    ]
    assert "clusters" not in payload
    assert "comparator_provenance" not in payload
    assert payload["source_commit"] == refresh.MUTALYZER_SHA


def test_corrupt_existing_fixture_raises_with_context(tmp_path: Path) -> None:
    """Error case: a truncated fixture must not read as 'just use --force'."""
    corrupt = tmp_path / "cases.json"
    corrupt.write_text('{"cases": [')
    with pytest.raises(SystemExit) as excinfo:
        refresh._load_existing(corrupt)
    message = str(excinfo.value)
    assert str(corrupt) in message
    assert "--force" in message
