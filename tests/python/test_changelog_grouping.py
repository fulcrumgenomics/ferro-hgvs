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
import re
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


# ---------------------------------------------------------------------------
# The squash-equivalent preview
# ---------------------------------------------------------------------------
#
# On a `pull_request` the commit this audit must judge does not exist yet: GitHub builds
# the squash subject from the PR title and the body from the PR description, and the
# trailer lives in the description. `ci.yml` synthesizes that commit before running the
# audit. These two pin the halves of that — the message shape the workflow builds, and the
# fact that it builds one at all — because a PR whose merged form misfiles its declaration
# would otherwise pass here and be caught only after it lands.


def test_a_trailer_that_arrives_only_in_the_pr_description_is_read() -> None:
    """The synthesized message is `<title> (#N)` + the description, verbatim.

    The trailer ends up in the footer position the real squash gives it, which is what
    `release-plz.toml`'s parsers match — so reading it here is the same question the
    changelog will ask after the merge.
    """
    title = "fix(normalize): converge the two spellings"
    body = "Why this changes.\n\nRepresentation-Change: 3 rows move, 2 merge / 1 respell"
    squash = f"{title} (#1560)\n\n{body}\n"
    assert trailer_value(squash) == "3 rows move, 2 merge / 1 respell"
    assert not opens_with_a_decline(trailer_value(squash) or "")


def test_a_declining_description_reads_as_a_decline_through_the_same_path() -> None:
    squash = "test(changelog): audit the section (#1560)\n\nBody.\n\nRepresentation-Change: none. Tests only.\n"
    value = trailer_value(squash)
    assert value == "none. Tests only."
    assert opens_with_a_decline(value)


def test_ci_synthesizes_the_squash_commit_before_auditing() -> None:
    """The audit is only as good as the history it is handed.

    Pinned against the workflow text rather than left to review: dropping the synthesis
    step turns this job into one that passes on every PR regardless of its declaration,
    and nothing else in the suite would notice.
    """
    workflow = (Path(__file__).resolve().parents[2] / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    audit = workflow.index("Audit the Representation changes section")
    synthesis = workflow.index("Synthesize the squash commit this PR will become")
    assert synthesis < audit, "the synthesis step must precede the audit step"

    # Scoped to the synthesis step itself, not to everything above the audit. Asserting
    # against the whole preamble would let a match come from an unrelated step — and it is
    # the *dataflow inside this step* that carries the property.
    step = workflow[synthesis:audit]

    # The subject and body must both come from the PR, and through `env` — a PR title is
    # attacker-controlled text and `${{ }}` inside `run` would splice it into the shell.
    assert "PR_TITLE: ${{ github.event.pull_request.title }}" in step
    assert "PR_BODY: ${{ github.event.pull_request.body }}" in step

    # Declaring the variables is not using them. Without this, replacing the message with
    # a literal — `commit --allow-empty -m 'ci: preview'` — leaves every assertion above
    # true while the audit runs against a commit carrying no trailer, so it would pass
    # every PR regardless of its declaration. That is the exact failure this test exists to
    # prevent, so pin the whole chain: env var -> message file -> commit.
    message_file = "squash-message.txt"
    written_from = re.search(
        r"printf\s+\S+\s+(?P<args>[^>]*)>\s*" + re.escape(message_file),
        step,
    )
    assert written_from is not None, (
        f"the synthesis step must build {message_file} with printf; without it the audit "
        "judges a commit that never carried the PR's trailer"
    )
    for variable in ('"$PR_TITLE"', '"$PR_NUMBER"', '"$PR_BODY"'):
        assert variable in written_from.group("args"), (
            f"{variable} is defined in `env` but never reaches {message_file}; the synthesized "
            "commit would not be the squash commit GitHub will write. Note the quotes are "
            "part of the assertion: unquoted, a multi-word PR title word-splits."
        )
    assert re.search(
        r"commit\s+--allow-empty(?:\s+\S+)*\s+--file\s+" + re.escape(message_file), step
    ), (
        f"the empty commit must be created from {message_file} (`--file`), or the message "
        "assembled above is discarded"
    )

    # And the history must be the preview too, not just the message. `actions/checkout`
    # hands us the MERGE ref on a pull_request, so without the soft reset the audited range
    # carries the branch's intermediate commits and the merge commit — none of which
    # survive a squash. Soft, so the PR's tree (and its `release-plz.toml`) is what gets
    # committed and read.
    assert "PR_BASE_SHA: ${{ github.event.pull_request.base.sha }}" in step, (
        "the base sha must reach the step through `env` for the soft reset"
    )
    reset = re.search(r"git\s+reset\s+(?P<mode>--\S+)\s+\"\$PR_BASE_SHA\"", step)
    assert reset is not None, (
        'the step must `git reset` to "$PR_BASE_SHA" before committing, or the audit '
        "judges commits the squash discards"
    )
    assert reset.group("mode") == "--soft", (
        f"the reset must be --soft, not {reset.group('mode')}: a mixed or hard reset would "
        "discard the PR's tree, so the audit would read the base's release-plz.toml and "
        "the synthesized commit would not carry the PR's changes"
    )
    assert step.index(reset.group(0)) < step.index("commit --allow-empty"), (
        "the reset must precede the commit, or the commit lands on the merge ref"
    )
