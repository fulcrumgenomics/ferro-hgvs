"""Unit tests for `scripts/inject_representation_disclosure.py` (#1556).

The script's job is to put the `Representation-Change:` trailer's own text under its
changelog bullet, which the changelog template cannot do. What is worth testing is not the
happy path -- that is visible by running it -- but the three rules that decide *how much* of
the message is the value, each of which was wrong at some point while the script was written
and each of which fails silently:

- where the value **ends** (a CodeRabbit block appended after the trailer is not the value);
- that a bullet with **no** matching commit or **no** trailer is reported rather than skipped;
- that a second run is a **no-op**, since release-plz regenerates the section on every push.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "inject_representation_disclosure.py"


def _load():
    """Import the script by path; `scripts/` is not a package."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("inject_representation_disclosure", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


module = _load()


SECTION = """## [unreleased]

### Representation changes

- *(normalize)* Move a thing ([#1234](https://example.invalid/1234))

### Other

- *(ci)* Something else ([#1235](https://example.invalid/1235))
"""


def test_the_value_stops_at_an_appended_coderabbit_block():
    """The trailer is documented as the last thing in the description and is not.

    GitHub appends CodeRabbit's auto-generated summary to the PR body, so it lands in the
    squash message *after* the trailer. Read to end-of-message, #1484's disclosure came out
    with three paragraphs of "Summary by CodeRabbit" quoted under it as the author's words.
    """
    message = (
        "fix(normalize): move a thing (#1234)\n\n"
        "Body prose.\n\n"
        "Representation-Change: 3 rows move of 500,004.\n"
        "Previously accepted, so a real migration.\n\n"
        "<!-- This is an auto-generated comment: release notes by coderabbit.ai -->\n"
        "## Summary by CodeRabbit\n"
        "- **Bug Fixes**\n"
    )
    value = module.disclosure_value(message)
    assert value == "3 rows move of 500,004.\nPreviously accepted, so a real migration."


def test_the_value_stops_at_the_next_trailer():
    message = (
        "fix: x (#1)\n\nRepresentation-Change: 2 rows move.\ncontinued.\nCloses: #99\ntrailing\n"
    )
    assert module.disclosure_value(message) == "2 rows move.\ncontinued."


def test_a_single_line_value_is_unchanged():
    assert module.disclosure_value("fix: x (#1)\n\nRepresentation-Change: none.\n") == "none."


def test_no_trailer_is_none():
    assert module.disclosure_value("fix: x (#1)\n\nJust a body.\n") is None


def test_an_indented_trailer_is_not_a_trailer():
    """Column 0, matching both existing checkers: an indented line is a continuation of the
    value above it, and reading it as a trailer is what #1573 fixed in the sibling script."""
    assert module.disclosure_value("fix: x (#1)\n\n  Representation-Change: none.\n") is None


def test_the_disclosure_is_attached_under_its_bullet():
    by_pr = {"1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"}
    rewritten, changed, problems = module.inject(SECTION, by_pr)
    assert problems == []
    assert changed == ["#1234"]
    assert "  > 3 rows move." in rewritten
    # And nothing outside the section is touched.
    assert "- *(ci)* Something else ([#1235](https://example.invalid/1235))\n" in rewritten
    assert rewritten.count("  > ") == 1


def test_a_second_run_changes_nothing():
    by_pr = {"1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"}
    once, _, _ = module.inject(SECTION, by_pr)
    twice, changed, problems = module.inject(once, by_pr)
    assert twice == once
    assert changed == []
    assert problems == []


@pytest.mark.parametrize(
    ("by_pr", "expected"),
    [
        ({}, "no commit in this history cites it"),
        (
            {"1234": "fix: move a thing (#1234)\n\nNo trailer here.\n"},
            "no\n`Representation-Change:` trailer",
        ),
    ],
)
def test_a_bullet_that_cannot_be_resolved_is_reported_not_skipped(by_pr, expected):
    """A silently-skipped bullet would reintroduce exactly the gap this script closes: a
    Representation changes entry with no disclosure, indistinguishable from one that needed
    none."""
    _, changed, problems = module.inject(SECTION, by_pr)
    assert changed == []
    assert len(problems) == 1
    assert expected.replace("\n", " ") in problems[0].replace("\n", " ")


def test_an_edited_disclosure_is_reported_rather_than_trusted():
    """The bullet is supposed to carry the trailer's own words, so a block that no longer
    matches its trailer is drift.

    This was silently green: the loop consuming the existing `  > ` lines then `continue`d
    past both the commit lookup and the trailer comparison, so a stale, truncated or
    hand-rewritten disclosure passed `--check` forever. The script could only ever be right
    on the run that first wrote it, which is the failure mode it exists to remove.
    """
    by_pr = {"1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"}
    injected, _, _ = module.inject(SECTION, by_pr)
    tampered = injected.replace("  > 3 rows move.", "  > 9000 rows move, actually.")

    _, changed, problems = module.inject(tampered, by_pr)

    assert changed == []
    assert len(problems) == 1
    assert "no longer matches its `Representation-Change:` trailer" in problems[0]
    # The message must show both sides; "they differ" alone is not actionable.
    assert "9000 rows move" in problems[0]
    assert "3 rows move" in problems[0]


def test_a_trailer_edited_upstream_is_caught_on_the_next_run():
    """The same guard from the other direction: the changelog is untouched and the *trailer*
    moved. Both are drift and both must be reported, because which side someone edited is
    not knowable from the artifact."""
    original = {"1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"}
    injected, _, _ = module.inject(SECTION, original)

    corrected = {"1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 4 rows move.\n"}
    _, changed, problems = module.inject(injected, corrected)

    assert changed == []
    assert len(problems) == 1
    assert "no longer matches" in problems[0]


def test_stdout_and_check_are_refused_together():
    """`--stdout` used to be tested first, so `--stdout --check` printed the rewritten file
    and returned 0 on a stale changelog -- a silent false pass in the one mode whose entire
    job is to fail."""
    with pytest.raises(SystemExit) as excinfo:
        module.main(["--stdout", "--check"])
    assert excinfo.value.code != 0


def test_write_mode_refuses_before_writing_when_a_bullet_is_unresolvable(tmp_path):
    """A partial artifact is not a better outcome than no artifact.

    The write used to happen before the exit code was computed, so a run with one resolvable
    and one unresolvable bullet left a PARTIALLY injected changelog on disk *and* exited 2 --
    the file looks updated and the failure looks transient.
    """
    two_bullets = SECTION.replace(
        "- *(normalize)* Move a thing ([#1234](https://example.invalid/1234))",
        "- *(normalize)* Move a thing ([#1234](https://example.invalid/1234))\n"
        "- *(normalize)* Move another ([#4321](https://example.invalid/4321))",
    )
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(two_bullets, encoding="utf-8")
    before = changelog.read_text(encoding="utf-8")

    # #1234 resolves; #4321 does not.
    module.squash_commits = lambda _repo: {
        "1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"
    }
    code = module.main(["--changelog", str(changelog)])

    assert code == 2
    assert changelog.read_text(encoding="utf-8") == before, "wrote a partially injected file"


def test_a_multi_paragraph_disclosure_is_not_reported_as_drift_against_itself():
    """The blank line inside a disclosure renders as a bare `  >`, not `  > `.

    `wrap_disclosure` rstrips and `DISCLOSURE_PREFIX` carries a trailing space, so a
    collector keyed only on the prefix stops at the first paragraph break. The block then
    reads back one line long and is reported as drift against the very trailer that
    produced it -- a red CI on a file the script itself wrote, which is the same
    unsatisfiable-check failure this module exists to prevent.
    """
    message = (
        "fix: move a thing (#1234)\n\n"
        "Representation-Change: 3 rows move.\n\n"
        "Previously accepted, so a real migration.\n"
    )
    by_pr = {"1234": message}

    once, changed, problems = module.inject(SECTION, by_pr)
    assert changed == ["#1234"]
    assert problems == []
    # The blank line really is emitted bare -- if this stops holding, the guard above is
    # solving a problem that no longer exists and should be re-derived, not kept.
    assert "  >\n" in once

    twice, changed_again, problems_again = module.inject(once, by_pr)
    assert twice == once
    assert changed_again == []
    assert problems_again == [], f"re-run reported drift against its own output: {problems_again}"
