"""
Tests for `scripts/check_representation_change.py` (#1522).

The check exists because an *absent* trailer and a *considered-and-declined* trailer
looked identical for a whole release cycle. So the assertions that matter here are the
ones separating those two states, and the one pinning that the watched set still matches
`release-plz.toml`'s stated scope — a directory silently dropping out of that tuple would
make this check pass on exactly the changes it exists to catch.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

_MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts" / "check_representation_change.py"
_SPEC = importlib.util.spec_from_file_location("check_representation_change", _MODULE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
check_representation_change = importlib.util.module_from_spec(_SPEC)
sys.modules["check_representation_change"] = check_representation_change
_SPEC.loader.exec_module(check_representation_change)

check = check_representation_change.check
find_declaration = check_representation_change.find_declaration
watched_files = check_representation_change.watched_files
WATCHED_PREFIXES = check_representation_change.WATCHED_PREFIXES


# ---------------------------------------------------------------------------
# The two states the whole check exists to separate
# ---------------------------------------------------------------------------


def test_watched_change_without_a_trailer_fails() -> None:
    """Silence is what this rejects — the v0.13.0 cycle's actual failure mode."""
    ok, message = check(["src/normalize/merge.rs"], "Fixes a thing.\n\nCloses #1.")
    assert not ok
    assert "src/normalize/merge.rs" in message
    assert "Representation-Change: none" in message, "the message must name the way to decline"


def test_watched_change_declaring_none_passes() -> None:
    """Declining is a declaration. #1521 is the real example: 46 changed lines in
    `src/normalize/merge.rs`, every one a comment."""
    ok, _ = check(
        ["src/normalize/merge.rs"],
        "Comments only; zero non-comment lines change.\n\nRepresentation-Change: none",
    )
    assert ok


def test_watched_change_declaring_a_move_passes() -> None:
    ok, message = check(
        ["src/spdi/mod.rs"],
        "Representation-Change: 577 rows move, 360 merge / 205 split / 12 respell",
    )
    assert ok
    assert "577 rows move" in message


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------


def test_unwatched_change_needs_no_declaration() -> None:
    ok, _ = check(["README.md", "tests/it/foo.rs", "src/cli/mod.rs"], "no trailer here")
    assert ok


def test_a_watched_file_among_unwatched_ones_still_requires_a_declaration() -> None:
    """The trigger is *any* watched file, not a majority of them."""
    ok, _ = check(["README.md", "docs/x.md", "src/project/projector.rs"], "no trailer")
    assert not ok


@pytest.mark.parametrize(
    "path",
    [
        "src/normalize/merge.rs",
        "src/hgvs/variant.rs",
        "src/spdi/mod.rs",
        "src/project/projector.rs",
    ],
)
def test_every_watched_prefix_triggers(path: str) -> None:
    assert watched_files([path]) == [path]
    ok, _ = check([path], "no trailer")
    assert not ok, f"{path} must require a declaration"


def test_watched_prefixes_match_the_release_config() -> None:
    """`release-plz.toml`'s reviewer checklist names the directories this must watch.

    Pinned as a set comparison rather than a count: two compensating edits could keep the
    count while swapping a directory out, which is exactly the drift that would make this
    check pass on the changes it exists to catch.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    for prefix in WATCHED_PREFIXES:
        directory = prefix.rstrip("/")
        assert f"`{directory}/`" in config, (
            f"{directory} is watched by the check but not named in release-plz.toml's "
            "checklist; the two must describe the same scope"
        )


# ---------------------------------------------------------------------------
# Trailer parsing
# ---------------------------------------------------------------------------


def test_trailer_is_found_case_insensitively() -> None:
    assert find_declaration("representation-change: none") == "none"
    assert find_declaration("REPRESENTATION-CHANGE: none") == "none"


def test_trailer_is_found_at_the_end_of_a_body() -> None:
    body = "Some description.\n\nCloses #1522.\n\nRepresentation-Change: none\n"
    assert find_declaration(body) == "none"


def test_absent_trailer_reads_as_none_not_as_a_declaration() -> None:
    assert find_declaration("Representation change: none") is None, (
        "the hyphenated trailer is what git-cliff matches; prose must not satisfy it"
    )


def test_an_empty_trailer_value_is_not_a_declaration() -> None:
    """`Representation-Change:` with nothing after it declares nothing."""
    assert find_declaration("Representation-Change:   \n") is None


@pytest.mark.parametrize("value", ["none", "None", "NONE", "no", "n/a", "na"])
def test_decline_spellings_all_pass(value: str) -> None:
    ok, _ = check(["src/normalize/merge.rs"], f"Representation-Change: {value}")
    assert ok


def test_an_empty_trailer_does_not_scavenge_the_next_line() -> None:
    """`\\s` spans newlines, so a permissive class would read `none` off the line below and
    call an empty trailer a declaration."""
    assert find_declaration("Representation-Change:\nnone\n") is None


def test_empty_trailer_on_a_watched_change_fails() -> None:
    ok, _ = check(["src/normalize/merge.rs"], "Representation-Change:   \n")
    assert not ok


def test_both_streams_on_stdin_is_rejected() -> None:
    """Reading both inputs from stdin would silently give the second an empty string,
    which for the declaration reads as "no trailer" — a wrong failure, not an error."""
    with pytest.raises(SystemExit) as excinfo:
        check_representation_change.main(["--changed-files", "-", "--declaration-file", "-"])
    assert excinfo.value.code == 2


# ---------------------------------------------------------------------------
# The changelog consumer must agree with this checker about what "declining" means
# ---------------------------------------------------------------------------


def test_decline_vocabulary_matches_the_changelog_config() -> None:
    """`release-plz.toml`'s exclusion rule and `NONE_VALUES` must accept the same words.

    They are two halves of one decision and they drift silently. If the checker accepts a
    decline the changelog rule does not, that PR passes CI and then renders under
    **Representation changes** as if it moved output — which is #1526, the bug this test
    was written for. If the changelog excludes a word the checker rejects, a contributor
    is told to declare something the changelog then hides.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    match = re.search(r'footer = "\(\?i\)\^Representation-Change:[^"]*?\(([^)]+)\)', config)
    assert match is not None, (
        "no case-insensitive Representation-Change exclusion rule found in release-plz.toml; "
        "the decline vocabulary cannot be checked"
    )
    configured = frozenset(match.group(1).split("|"))
    assert configured == check_representation_change.NONE_VALUES, (
        f"release-plz.toml excludes {sorted(configured)} but the checker treats "
        f"{sorted(check_representation_change.NONE_VALUES)} as declines; a word in one and "
        "not the other either leaks a non-change into the changelog or hides a real one"
    )


def test_both_representation_change_rules_are_case_insensitive() -> None:
    """A case-sensitive *inclusion* rule drops a lowercase `representation-change:` — a real
    disclosure — silently into `Other`, since the checker accepts any casing."""
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    rules = re.findall(r'\{ footer = "([^"]*Representation-Change[^"]*)"', config)
    assert len(rules) == 2, f"expected an exclusion and an inclusion rule, found {rules}"
    for rule in rules:
        assert rule.startswith("(?i)"), (
            f"rule {rule!r} is case-sensitive; the checker's own trailer regex is not, so the "
            "two disagree about whether `REPRESENTATION-CHANGE:` is a declaration"
        )


def test_contributing_documents_the_same_decline_vocabulary() -> None:
    """CONTRIBUTING.md is the third place these words appear, and prose drifts fastest.

    A contributor reads the doc, not the config or this script, so a word documented here
    but not accepted by both is guidance that fails CI, and one accepted but not documented
    is a decline nobody knows they can make.
    """
    doc = (Path(__file__).resolve().parents[2] / "CONTRIBUTING.md").read_text(encoding="utf-8")
    documented = {
        word
        for word in re.findall(r"`([a-z/]+)`(?=[,\s]|$)", doc)
        if word in check_representation_change.NONE_VALUES
    }
    assert documented == check_representation_change.NONE_VALUES, (
        f"CONTRIBUTING.md documents {sorted(documented)} as declines but the checker accepts "
        f"{sorted(check_representation_change.NONE_VALUES)}"
    )


def test_the_decline_exclusion_precedes_the_inclusion() -> None:
    """git-cliff takes the FIRST matching parser, so swapping these two lines silently
    restores #1526 in full.

    Verified against git-cliff 2.13.1 rather than assumed: with the exclusion second, all
    four trailers -- `none`, `NONE`, and both real disclosures -- group under
    **Representation changes**. Every other guard in this file still passes in that state,
    which is why the ordering needs one of its own.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    rules = re.findall(r'\{ footer = "([^"]*Representation-Change[^"]*)"', config)
    assert len(rules) == 2, f"expected an exclusion and an inclusion rule, found {rules}"
    exclusion, inclusion = rules
    assert any(word in exclusion for word in check_representation_change.NONE_VALUES), (
        f"the first Representation-Change rule is {exclusion!r}, which does not name a decline "
        "value; the exclusion must come first or every decline is grouped as a real change"
    )
    assert not any(word in inclusion for word in check_representation_change.NONE_VALUES), (
        f"the second rule is {inclusion!r}, which looks like the exclusion -- the two are reversed"
    )
