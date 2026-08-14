"""Unit tests for `scripts/inject_representation_disclosure.py` (#1556).

The script's job is to put the `Representation-Change:` trailer's own text under its
changelog bullet, which the changelog template cannot do. What is worth testing is not the
happy path -- that is visible by running it -- but the three rules that decide *how much* of
the message is the value, each of which was wrong at some point while the script was written
and each of which fails silently:

- where the value **ends** (a CodeRabbit block appended after the trailer is not the value);
- that a bullet with **no** matching commit or **no** trailer is reported rather than skipped;
- that a second run is a **no-op**, since release-plz regenerates the section on every push.

`EDITORIAL_CORRECTIONS` (#1886) is tested to the same standard and for the same reason. It
exists because a merged trailer cannot be corrected, so the correction has to be re-rendered
on every run rather than written into `CHANGELOG.md` once -- which makes "a second run is a
no-op" and "a hand-edit is drift" load-bearing for it too, not merely nice.
"""

from __future__ import annotations

import importlib.util
import re
import subprocess
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


def test_write_mode_refuses_before_writing_when_a_bullet_is_unresolvable(tmp_path, monkeypatch):
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

    # #1234 resolves; #4321 does not. Stubbed through `monkeypatch` rather than by rebinding
    # `module.squash_commits` directly: a bare rebind is never undone, so every test after
    # this one in the module silently sees a two-entry history instead of the repository's.
    monkeypatch.setattr(
        module,
        "squash_commits",
        lambda _repo: {
            "1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"
        },
    )
    code = module.main(["--changelog", str(changelog)])

    assert code == 2
    assert changelog.read_text(encoding="utf-8") == before, "wrote a partially injected file"


TRAILERED = {"1234": "fix: move a thing (#1234)\n\nRepresentation-Change: 3 rows move.\n"}


@pytest.fixture
def one_correction(monkeypatch):
    """Register a correction for #1234 for the duration of a test."""
    monkeypatch.setattr(
        module,
        "EDITORIAL_CORRECTIONS",
        {"1234": ["**Editorial correction (#4321).** It was 4 rows, not 3."]},
    )


def test_nothing_is_attached_for_a_pr_with_no_registered_correction():
    """The default is silence. A mechanism that decorated every disclosure would dilute the
    label that makes a correction recognisable as one."""
    rewritten, _, problems = module.inject(SECTION, TRAILERED)
    assert problems == []
    assert "Editorial correction" not in rewritten


def test_a_registered_correction_is_attached_below_the_trailers_own_words(one_correction):
    """Below, and separated -- the trailer's words must stay readable as the trailer's.

    The separator is not cosmetic: without a blank quote line the correction would be a lazy
    continuation of the trailer's last paragraph, i.e. it would render as the author's own
    closing sentence, which is the one thing a correction must never look like.
    """
    rewritten, changed, problems = module.inject(SECTION, TRAILERED)

    assert problems == []
    assert changed == ["#1234"]
    body = rewritten.splitlines()
    start = body.index("  > 3 rows move.")
    assert body[start + 1] == "  >"
    assert body[start + 2] == "  > **Editorial correction (#4321).** It was 4 rows, not 3."


def test_a_correction_is_idempotent(one_correction):
    """release-plz regenerates the section on every push to `main`, so the correction is
    re-attached rather than written once -- which only works if re-running is a no-op."""
    once, _, _ = module.inject(SECTION, TRAILERED)
    twice, changed, problems = module.inject(once, TRAILERED)

    assert twice == once
    assert changed == []
    assert problems == []


@pytest.fixture
def two_corrections(monkeypatch):
    """Register a *two*-paragraph correction, which is the shape actually shipped.

    `EDITORIAL_CORRECTIONS["1835"]` holds two paragraphs, so the second iteration of
    `wrap_correction`'s loop -- the one that emits a bare `  >` *between* corrections, and
    that `_is_disclosure_line`'s second arm has to collect back for idempotence -- is the
    only path the real register uses and was the one path no test ran.
    """
    monkeypatch.setattr(
        module,
        "EDITORIAL_CORRECTIONS",
        {
            "1234": [
                "**Editorial correction (#4321).** It was 4 rows, not 3.",
                "**Editorial correction (#5678).** The count excludes two declined rows.",
            ]
        },
    )


def test_two_corrections_are_each_separated_and_each_labelled(two_corrections):
    """Appending rather than rewriting is the register's stated contract, so two corrections
    to one disclosure must both survive, in order, each behind its own bare `  >`."""
    rewritten, changed, problems = module.inject(SECTION, TRAILERED)

    assert problems == []
    assert changed == ["#1234"]
    body = rewritten.splitlines()
    start = body.index("  > 3 rows move.")
    assert body[start + 1 : start + 6] == [
        "  >",
        "  > **Editorial correction (#4321).** It was 4 rows, not 3.",
        "  >",
        "  > **Editorial correction (#5678).** The count excludes two declined rows.",
        "",
    ]


def test_two_corrections_are_idempotent(two_corrections):
    """The separator between two corrections is a bare `  >`, which only `_is_disclosure_line`'s
    second arm matches. Miss it and the re-read stops at the first paragraph break, so the
    block reads back short and is reported as drift against a changelog nothing touched."""
    once, _, _ = module.inject(SECTION, TRAILERED)
    twice, changed, problems = module.inject(once, TRAILERED)

    assert twice == once
    assert changed == []
    assert problems == []


def test_deleting_a_registered_correction_from_the_changelog_is_reported_as_drift(
    one_correction,
):
    """The register is the source of truth for a correction exactly as the commit message is
    for the trailer, so an edited changelog copy has to fail `--check` the same way. Without
    this, a correction could be silently removed from the released section it corrects."""
    injected, _, _ = module.inject(SECTION, TRAILERED)
    tampered = "\n".join(
        line for line in injected.splitlines() if "Editorial correction" not in line
    )

    _, changed, problems = module.inject(tampered + "\n", TRAILERED)

    assert changed == []
    assert len(problems) == 1
    assert "no longer matches" in problems[0]


def test_every_registered_correction_names_itself_and_cites_an_issue():
    """A guard on the real register, not on a fixture.

    Inside the trailer's blockquote the label is the *only* thing separating editorial text
    from the author's own words, so a paragraph that omits it silently converts a correction
    into an apparent part of the disclosure -- the precise failure this mechanism exists to
    avoid. Checked here rather than left to review, because it is checkable.

    The cite is matched as `(#<digits>)` rather than by prefix: a bare `**Editorial
    correction (#` names no issue at all and would otherwise pass, which is the half of the
    label a reader actually follows.
    """
    for number, paragraphs in module.EDITORIAL_CORRECTIONS.items():
        assert paragraphs, f"#{number} registers an empty correction"
        for paragraph in paragraphs:
            assert re.match(r"\*\*Editorial correction \(#\d+\)", paragraph), (
                f"#{number}: a correction must open by naming itself editorial and citing "
                f"the issue that raised it; got {paragraph[:60]!r}"
            )


def _repository_is_shallow() -> bool:
    """Is this checkout missing history? `git log` then answers about one commit, not the repo."""
    probe = subprocess.run(
        ["git", "rev-parse", "--is-shallow-repository"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return probe.stdout.strip() == "true"


def test_every_registered_correction_key_is_a_pr_number():
    """Hermetic half of the key check, so it runs everywhere the suite does.

    `inject` reads the register through `EDITORIAL_CORRECTIONS.get(number, [])` and `number`
    comes from `bullet_pr_number`, which returns the digits out of `(#N)` / `[#N]`. A key that
    is not those digits -- `"#1835"`, `"1835 "`, `"PR-1835"` -- can never match any bullet, and
    the failure is total silence: nothing renders, nothing is reported, `--check` stays green.
    """
    for number in module.EDITORIAL_CORRECTIONS:
        assert number.isdigit(), (
            f"{number!r} is registered in EDITORIAL_CORRECTIONS but is not a bare PR number, "
            "so `bullet_pr_number` can never produce it and the correction is unreachable"
        )


def test_every_registered_correction_names_a_pr_that_carries_a_trailer():
    """The half that needs history: the key must name a real, trailered PR.

    A transposed digit passes the hermetic check above and is still silent -- and the
    changelog cannot answer it either, because a correction is routinely registered while its
    bullet is still unrendered, in the window before release-plz cuts the section. What *is*
    checkable is the commit: the key must name a PR in this history whose squash message
    carries a `Representation-Change:` trailer, which is the only kind of bullet a correction
    can ever attach to.

    Skipped, loudly and by name, on a shallow checkout. `actions/checkout` defaults to
    `fetch-depth: 1` and only the two jobs that need `<latest tag>..HEAD` override it
    (`ci.yml:335`, `:395`); `Python Wheel Test`, which is where this file runs in CI, does
    not -- so `squash_commits` would see one commit there and this would fail for the
    environment rather than for the register. It arms in every full clone, which is where a
    key is typed and where the pre-push hooks run.
    """
    if _repository_is_shallow():
        pytest.skip("shallow checkout: `git log` cannot see the commit a key names")

    by_pr = module.squash_commits(REPO_ROOT)
    for number in module.EDITORIAL_CORRECTIONS:
        message = by_pr.get(number)
        assert message is not None, (
            f"#{number} has a registered editorial correction but no commit in this history "
            "cites it -- check the PR number for a typo"
        )
        assert module.disclosure_value(message) is not None, (
            f"#{number} has a registered editorial correction but its commit carries no "
            "`Representation-Change:` trailer, so it can never appear under a "
            "Representation changes bullet"
        )


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
