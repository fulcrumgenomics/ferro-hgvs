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
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_MODULE_PATH = _REPO_ROOT / "scripts" / "check_changelog_grouping.py"
_SPEC = importlib.util.spec_from_file_location("check_changelog_grouping", _MODULE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
check_changelog_grouping = importlib.util.module_from_spec(_SPEC)
sys.modules["check_changelog_grouping"] = check_changelog_grouping
_SPEC.loader.exec_module(check_changelog_grouping)

load_changelog_config = check_changelog_grouping.load_changelog_config
opens_with_a_decline = check_changelog_grouping.opens_with_a_decline
parse_rendered_groups = check_changelog_grouping.parse_rendered_groups
render = check_changelog_grouping.render
rendered_headings = check_changelog_grouping.rendered_headings
strip_ordering_prefix = check_changelog_grouping.strip_ordering_prefix
trailer_value = check_changelog_grouping.trailer_value
write_cliff_config = check_changelog_grouping.write_cliff_config


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
        # A missing space after the terminator is an ordinary typo, and it used to flip the
        # verdict: `split()[0]` yielded `none.Tests`, so the checker called this a decline
        # and the audit called it a move — a red build on a correct trailer.
        "none.Tests only.",
        "none—no watched file is touched",
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
        # `CONTRIBUTING.md:103-105` documents both of these as filed as real changes: a comma
        # is not a terminator, and `no rows move` has no terminator at all. "Both err toward
        # listing a change rather than hiding one." The audit used to call them declines,
        # which made them unsatisfiable — release-plz filed them one way and the audit failed
        # them the other, with no trailer text able to satisfy both.
        "none, except two rows that merge",
        "no rows move",
    ],
)
def test_a_description_of_a_move_is_not_a_decline(value: str) -> None:
    assert not opens_with_a_decline(value)


def test_the_audit_and_the_checker_agree_on_every_documented_form() -> None:
    """The two halves must reach the same verdict, or CI is unsatisfiable.

    This is the guard the audit shipped without. Its rule was deliberately coarser than the
    checker's — first word, punctuation stripped, no terminator logic — on the reasoning that
    a second opinion derived the same way is worth nothing. But a *disagreement* is not a
    second opinion, it is a red build: `no rows move` and `none, except two rows that merge`
    are filed as real changes by design (`CONTRIBUTING.md:103-105`) and the audit called them
    declines, so no trailer text could satisfy both halves. Once such a commit is on `main`
    the job is red for every open PR until the next release tag.

    What is still independent is the derivation — this check renders the changelog through
    real git-cliff and reads the result; the checker regexes the PR body. That is pinned
    separately by `test_the_rule_shares_no_vocabulary_with_the_checker`. Agreement on the
    verdict and independence of derivation are different properties, and #1555 is the reason
    to want both: a test comparing only vocabulary found the two halves agreeing while both
    were wrong.
    """
    checker_path = (
        Path(__file__).resolve().parents[2] / "scripts" / "check_representation_change.py"
    )
    spec = importlib.util.spec_from_file_location("_checker_for_agreement", checker_path)
    assert spec is not None and spec.loader is not None
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)

    for value in [
        # declines, including every terminator CONTRIBUTING.md names
        "none",
        "NONE",
        "none.",
        "none. Tests only; nothing under a watched directory.",
        "none: comments only",
        "none; one fixture JSON",
        "none — 0 of 950 rows move",
        "none.Tests only.",
        "no. Nothing reaches a normalizer.",
        "n/a: docs",
        "na. generated artifact only",
        # real changes, including the two the coarse rule got wrong
        "577 rows move, 360 merge / 205 split",
        "0 rows move over 5,761,302 real expressions",
        "no rows move",
        "none, except two rows that merge",
        "nothing moves",
        "not measured yet",
    ]:
        assert opens_with_a_decline(value) == checker.declines(value), (
            f"{value!r}: the audit calls it "
            f"{'a decline' if opens_with_a_decline(value) else 'a move'} and the checker calls "
            f"it {'a decline' if checker.declines(value) else 'a move'}. A disagreement is a "
            "red build on a correct trailer, not a useful second opinion — the independence "
            "that matters is how each derives its answer, not that they answer differently."
        )


def test_the_rule_shares_no_vocabulary_with_the_checker() -> None:
    """The audit must not *import* the checker's rule, so the two can be wrong separately.

    #1555 survived a test written to catch exactly it, because that test compared the two
    halves against *each other* and they were wrong together.

    Note what this does and does not buy, because the original framing over-claimed and that
    over-claim is what shipped the defect above. It keeps the implementations separate. It
    does **not** make the audit's *verdict* independent — the verdict has to agree, or CI is
    unsatisfiable (`test_the_audit_and_the_checker_agree_on_every_documented_form`). The
    independence that catches a shared mistake is in the *derivation*: this check renders the
    changelog through real git-cliff and reads the result, which is the question no
    vocabulary-comparison test was asking.
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
    raw_step = workflow[synthesis:audit]

    # And with `#`-comment lines removed, because YAML comments are not code. Adversarial
    # review found three mutants that passed both this test and `actionlint` purely by
    # matching comment text: commenting the reset out (`# git reset --soft "$PR_BASE_SHA"`)
    # restored the pre-fix behaviour silently, and `# was: git reset --soft "$PR_BASE_SHA"`
    # above a `reset --soft HEAD~0` did the same. Commenting a line out to debug is an
    # ordinary edit, so the assertions must not be satisfiable by the prose that explains
    # them. `step` is what every assertion below reads.
    step = "\n".join(line for line in raw_step.splitlines() if not line.lstrip().startswith("#"))
    assert "#" in raw_step and len(step) < len(raw_step), (
        "comment stripping matched nothing, so it is not protecting anything — the step's "
        "explanatory comments were expected to be present and removed"
    )

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

    # Written exactly once. `[^>]*` below is positional, so a second redirect appended after
    # the good `printf` — `printf 'ci: preview\n' > squash-message.txt` — clobbers the file
    # while the first match still succeeds. Adversarial review found that mutant surviving
    # both this test and `actionlint`.
    writes = len(re.findall(r">\s*" + re.escape(message_file), step))
    assert writes == 1, (
        f"{message_file} is redirected to {writes} times in this step; a second write "
        "clobbers the message assembled from the PR and the audit then judges a commit that "
        "never carried the trailer"
    )

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

    # The tag must be resolved on the MERGE REF, i.e. before the reset moves HEAD's
    # ancestry, and handed to the audit explicitly.
    #
    # `base.sha` is the base tip as of the last PR event while the merge ref uses the
    # *current* base, so once a release lands, an un-resynchronized PR resets to a commit
    # from which the new tag is unreachable. `describe` then answers with the PREVIOUS tag
    # and the audit silently widens back over an entire released cycle, re-judging merged
    # work the PR author cannot touch. Measured: 3 commits audited that way against 1.
    describe = re.search(r"(?P<var>[A-Z_]+)=\$\(git describe --tags --abbrev=0", step)
    assert describe is not None, (
        "the step must resolve the latest tag itself; leaving it to the audit resolves it "
        "after the reset, against the wrong ancestry"
    )
    assert step.index(describe.group(0)) < step.index(reset.group(0)), (
        "the tag must be resolved BEFORE the soft reset — after it, a release that landed "
        "since the PR's base.sha is unreachable and `describe` silently returns the previous "
        "tag, rewinding the audited range over an already-released cycle"
    )
    assert "AUDIT_RANGE=" in step and "GITHUB_ENV" in step, (
        "the resolved range must be exported (AUDIT_RANGE -> $GITHUB_ENV) so the audit step "
        "uses it rather than re-deriving it from the reset HEAD"
    )


# ---------------------------------------------------------------------------
# The grouping of a declining `fix:` (#1557), rendered through real git-cliff
# ---------------------------------------------------------------------------
#
# The bug: a commit whose `Representation-Change:` trailer DECLINES (`none`, ...) was filed
# under **Other** whatever its conventional type, so a `fix:` that correctly declined never
# reached **Fixed**. git-cliff evaluates a parser's fields as OR (not AND), so the naive
# `{ footer = <decline>, message = "^fix" }` cannot express "a declining fix"; the fix is a
# `commit_preprocessor` in `release-plz.toml` that renames a declining trailer's key so the
# commit falls through to its conventional-type parser. These render the REAL configuration
# through real git-cliff and assert where each commit lands — the only test that can, because
# the grouping is git-cliff's, not this script's.


_GIT_CLIFF = shutil.which("git-cliff")
requires_git_cliff = pytest.mark.skipif(
    _GIT_CLIFF is None,
    reason="git-cliff is not installed; the grouping is rendered by git-cliff itself",
)


# A hermetic git environment: the throwaway repo must NOT inherit the developer's global
# config. That config commonly enables commit signing (this repo signs via 1Password), which
# would prompt or hang on `git commit` in a headless test, and could pull in global hooks. The
# identity is supplied through env vars so no `git config` step is needed, and signing is
# forced off. `os.devnull` for the config files isolates global and system config entirely.
_GIT_ENV = {
    **os.environ,
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_SYSTEM": os.devnull,
    "GIT_AUTHOR_NAME": "test",
    "GIT_AUTHOR_EMAIL": "test@example.com",
    "GIT_COMMITTER_NAME": "test",
    "GIT_COMMITTER_EMAIL": "test@example.com",
}


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=repo, capture_output=True, text=True, check=True, env=_GIT_ENV
    ).stdout


def _render_real_config(commits: list[str]) -> dict[str, list[str]]:
    """Build a throwaway repo of `commits`, render it with the real config, return groups.

    Returns `{group name: [commit subject, ...]}` — subjects rather than ids, because a test
    reads better asserting on `"declining fix"` than on a 40-hex sha. The configuration is
    `release-plz.toml`'s own `[changelog]`, written out by the audit's `write_cliff_config`,
    so this exercises exactly what the release will render.
    """
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp)
        _git(repo, "init", "-q")
        # A base commit so the audited range `<base>..HEAD` is a well-formed OID range.
        _git(repo, "commit", "-q", "--no-gpg-sign", "--allow-empty", "-m", "chore: base")
        base = _git(repo, "rev-parse", "HEAD").strip()
        for message in commits:
            _git(repo, "commit", "-q", "--no-gpg-sign", "--allow-empty", "-m", message)

        config = repo / "cliff.toml"
        write_cliff_config(load_changelog_config(_REPO_ROOT / "release-plz.toml"), config)
        rendered = render(f"{base}..HEAD", config, repo)

        subject_of = {
            line.split(" ", 1)[0]: line.split(" ", 1)[1]
            for line in _git(repo, "log", "--format=%H %s").strip().splitlines()
        }
        groups = parse_rendered_groups(rendered)
        return {group: [subject_of[i] for i in ids] for group, ids in groups.items()}


@requires_git_cliff
def test_a_declining_fix_is_filed_under_fixed_not_other() -> None:
    """The #1557 regression: a `fix:` that declined must land under **Fixed**.

    Rendered through the real `release-plz.toml`. Before the fix this subject rendered under
    **Other**; the assertion is on the group git-cliff actually assigns.
    """
    subject = "fix(changelog): a declining fix"
    groups = _render_real_config(
        [f"{subject}\n\nBody.\n\nRepresentation-Change: none. Tests only."]
    )
    assert subject in groups.get("Fixed", []), groups
    assert subject not in groups.get("Other", []), groups
    # And it must NOT be mistaken for a real representation change.
    assert subject not in groups.get("Representation changes", []), groups


@requires_git_cliff
def test_declines_group_by_type_while_real_changes_group_together() -> None:
    """The whole routing, in one render: declines by type, real changes into their section.

    Covers the three cases `release-plz.toml` documents, plus the two documented edge forms
    that read as real changes (`none, except ...` and `no rows move`) so a decline regex that
    swallowed them — the failure this config's regex is built to avoid — would show up here.
    """
    groups = _render_real_config(
        [
            "fix(a): declining fix\n\nRepresentation-Change: none. Tests only.",
            "feat(b): declining feature\n\nRepresentation-Change: none",
            "fix(c): plain fix, no trailer",
            "fix(d): a real move\n\nRepresentation-Change: 3 rows move",
            "fix(e): none-comma is a real change\n\nRepresentation-Change: none, except two rows merge",
            "fix(f): no-rows-move is real\n\nRepresentation-Change: no rows move",
            "chore(g): housekeeping",
        ]
    )
    assert set(groups.get("Fixed", [])) == {
        "fix(a): declining fix",
        "fix(c): plain fix, no trailer",
    }, groups
    assert groups.get("Added", []) == ["feat(b): declining feature"], groups
    assert set(groups.get("Representation changes", [])) == {
        "fix(d): a real move",
        "fix(e): none-comma is a real change",
        "fix(f): no-rows-move is real",
    }, groups
    assert groups.get("Other", []) == ["chore(g): housekeeping"], groups


@requires_git_cliff
def test_declining_commits_are_absent_from_representation_changes() -> None:
    """A declining trailer must never be filed as a real change (the sibling of #1557).

    The preprocessor and the audit share one decline regex, so this and the audit agree by
    construction — but this asserts it against the rendered output, not the regex.
    """
    groups = _render_real_config(
        [
            "fix(a): bare none\n\nRepresentation-Change: none",
            "fix(b): none dot\n\nRepresentation-Change: none.",
            "feat(c): NONE upper\n\nRepresentation-Change: NONE",
            "fix(d): n/a\n\nRepresentation-Change: n/a: docs",
            "fix(e): dash reason\n\nRepresentation-Change: no — nothing reaches a normalizer",
        ]
    )
    section = groups.get("Representation changes", [])
    assert section == [], f"declines leaked into Representation changes: {section}"


def test_config_carries_the_decline_neutralizing_preprocessor() -> None:
    """A structural guard: the fix is a `commit_preprocessors` entry, not a parser.

    Rendering (above) skips where git-cliff is absent, e.g. the Python Wheel Test job. This
    runs everywhere and fails loudly if the preprocessor is dropped or the old
    route-declines-to-a-fixed-group parser is reintroduced.
    """
    changelog = load_changelog_config(_REPO_ROOT / "release-plz.toml")
    preprocessors = changelog.get("commit_preprocessors", [])
    assert any(
        "Representation-Change" in p.get("pattern", "")
        and "none" in p.get("pattern", "")
        and "Declined" in (p.get("replace") or "")
        for p in preprocessors
    ), f"the decline-neutralizing preprocessor is missing: {preprocessors}"
    # The default PR-link preprocessor must be reproduced, or the changelog loses every link
    # (setting commit_preprocessors replaces release-plz's defaults wholesale).
    assert any(
        r"(#" in p.get("pattern", "") and "pull/" in (p.get("replace") or "") for p in preprocessors
    ), f"the PR-link preprocessor is missing, so links would vanish: {preprocessors}"
    # No parser may route a `Representation-Change` footer to a fixed non-section group —
    # that is the pre-#1557 shape that buried a declining `fix:` under Other.
    for parser in changelog.get("commit_parsers", []):
        footer = parser.get("footer", "")
        if "Representation-Change" in footer:
            assert parser.get("group", "").endswith("Representation changes"), (
                f"a footer parser routes a Representation-Change trailer to {parser.get('group')!r}; "
                "declines are handled by the preprocessor now, so the only footer parser must be "
                "the inclusion rule (#1557)"
            )
