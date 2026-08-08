"""
Tests for `scripts/check_changelog_grouping.py`.

The script renders the changelog with the real configuration and audits the result, so its
end-to-end behaviour needs git-cliff and a commit range. What is unit-tested here is the
part that decides *verdicts* — the independent decline rule, and the parsing that maps a
rendered section back to commits — because that is where a mistake would make the audit
silently agree with whatever it was handed.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts" / "check_changelog_grouping.py"
_SPEC = importlib.util.spec_from_file_location("check_changelog_grouping", _MODULE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
check_changelog_grouping = importlib.util.module_from_spec(_SPEC)
sys.modules["check_changelog_grouping"] = check_changelog_grouping
_SPEC.loader.exec_module(check_changelog_grouping)

opens_with_a_decline = check_changelog_grouping.opens_with_a_decline
parse_rendered_groups = check_changelog_grouping.parse_rendered_groups
rendered_headings = check_changelog_grouping.rendered_headings
strip_ordering_prefix = check_changelog_grouping.strip_ordering_prefix
trailer_value = check_changelog_grouping.trailer_value


# ---------------------------------------------------------------------------
# The independent decline rule
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        "none",
        "NONE",
        "none.",
        "none. Tests only; nothing under a watched directory.",
        "no — nothing reaches a normalizer",
        "n/a: docs",
        "na. generated artifact only",
        # Not a decline by the *checker's* rule (no terminator), and that is the point:
        # this rule is coarser on purpose, so that a shared mistake cannot hide.
        "none, except two rows that merge",
    ],
)
def test_a_leading_decline_word_is_a_decline(value: str) -> None:
    assert opens_with_a_decline(value)


@pytest.mark.parametrize(
    "value",
    [
        "577 rows move, 360 merge / 205 split",
        "0 rows move over 5,761,302 real expressions",
        "3 rows of 500,004 move",
        "nothing moves under src/normalize/",
        "not measured yet",
    ],
)
def test_a_description_of_a_move_is_not_a_decline(value: str) -> None:
    assert not opens_with_a_decline(value)


def test_the_rule_shares_no_vocabulary_with_the_checker() -> None:
    """This is the whole design of the audit, so it is pinned rather than left to comments.

    #1555 survived a test written to catch exactly it, because that test compared the two
    halves against *each other* and they were wrong together. A second opinion has to be
    derived differently to be worth anything — here, a plain first-word test with no
    terminator logic.
    """
    source = _MODULE_PATH.read_text(encoding="utf-8")
    for borrowed in ("import check_representation_change", "NONE_VALUES", "DECLINE_RE"):
        assert borrowed not in source, (
            f"the audit reuses {borrowed!r} from the checker; a shared rule cannot catch a "
            "shared mistake, which is exactly how #1555 passed a test written to catch it"
        )


# ---------------------------------------------------------------------------
# Mapping a rendered section back to commits
# ---------------------------------------------------------------------------


_IDS = [prefix.ljust(40, "0") for prefix in ("5616cdb9", "d362a47d", "ea5bd8ce")]

_RENDERED = f"""
## [unreleased] - 2026-08-08

### <!-- 0 -->Representation changes
{_IDS[0]}
{_IDS[1]}

### <!-- 7 -->Other
{_IDS[2]}
"""


def test_groups_are_keyed_without_the_ordering_prefix() -> None:
    """A broken postprocessor must not also break the lookup.

    Keeping them coupled made one failure hide the other: with the prefix present, the
    section lookup found nothing and the audit reported the two real disclosures as
    *unfiled* rather than reporting the declines that were filed. Each failure has to name
    itself, so the key is normalised and the prefix is reported separately.
    """
    groups = parse_rendered_groups(_RENDERED)
    assert set(groups) == {"Representation changes", "Other"}
    assert groups["Representation changes"] == _IDS[:2]
    assert groups["Other"] == _IDS[2:]


def test_headings_are_reported_verbatim() -> None:
    """The prefix check reads the raw headings, which is what release-plz writes."""
    headings = rendered_headings(_RENDERED)
    assert "<!-- 0 -->Representation changes" in headings


@pytest.mark.parametrize(
    ("heading", "expected"),
    [
        ("<!-- 0 -->Representation changes", "Representation changes"),
        ("<!--0-->Other", "Other"),
        ("Representation changes", "Representation changes"),
    ],
)
def test_strip_ordering_prefix(heading: str, expected: str) -> None:
    assert strip_ordering_prefix(heading) == expected


# ---------------------------------------------------------------------------
# Trailer extraction
# ---------------------------------------------------------------------------


def test_trailer_is_read_from_a_full_commit_message() -> None:
    message = "fix(normalize): a thing\n\nBody.\n\nRepresentation-Change: 3 rows move\n"
    assert trailer_value(message) == "3 rows move"


def test_a_message_without_a_trailer_reads_as_none() -> None:
    assert trailer_value("chore: something\n\nNo declaration here.\n") is None
