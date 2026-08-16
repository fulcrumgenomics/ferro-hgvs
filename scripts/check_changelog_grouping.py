#!/usr/bin/env python3
"""
Check that the changelog's **Representation changes** section contains what it claims.

`release-plz.toml` routes commits carrying a `Representation-Change:` trailer into their own
changelog section, and `scripts/check_representation_change.py` makes PRs carry one. Both
halves have now been wrong in production, in opposite directions:

- **#1526** listed a decline as a representation change (bare `none`).
- **#1555** listed eight of them, because a decline that gave a reason -- `none. <why>` --
  did not match a rule anchored on the bare word. The v0.13.1 section listed 10 entries of
  which 2 were real.

Neither was caught by CI, and the reason is worth stating because it shapes this check.
#1527 added `test_decline_vocabulary_matches_the_changelog_config`, whose whole purpose is
to stop those two halves drifting apart. It passed throughout #1555, because the two halves
had **not** drifted: they named the same four words and were wrong together. An
agreement-based test cannot see a shared mistake.

So this check does not ask whether the two halves agree. It renders the changelog with the
real configuration and asks a deliberately *simpler*, independent question of the result:

    does anything filed under "Representation changes" begin with a word that means no?

That rule shares no code with the checker, and it would have fired on all eight #1555 rows
and on #1526 regardless of what either half believed at the time. It is coarse on purpose --
a second opinion is only useful when it is derived differently.

Usage:
    python scripts/check_changelog_grouping.py                  # since the latest tag
    python scripts/check_changelog_grouping.py --range v0.13.0..HEAD

Exits 0 when the section is consistent, 1 when it is not.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

# `tomllib` is standard from 3.11. The workflow pins 3.12; ruff groups it as third-party
# here only because the project's `target-version` is py310 for the published package.
import tomllib

REPO_ROOT = Path(__file__).resolve().parents[1]

#: The group `release-plz.toml` routes trailered commits into. Matched after the ordering
#: prefix is stripped, so it survives a change to the `<!-- N -->` numbering.
REPRESENTATION_GROUP = "Representation changes"

#: Words that mean "this moves nothing", as the *verdict*.
DECLINE_WORDS = frozenset({"none", "no", "n/a", "na"})

#: Punctuation that may introduce the reason for a decline. A comma is deliberately absent,
#: matching `CONTRIBUTING.md`: it usually introduces a qualification that changes the answer
#: ("none, except two rows"), which is filed as a real change.
DECLINE_TERMINATORS = ".;:—–"

#: A decline: a verdict word, optionally followed by a terminator and a reason, and nothing
#: else. Built from the two constants above.
#:
#: **This rule agrees with `scripts/check_representation_change.py` on the verdict, and it
#: has to.** An earlier revision was deliberately coarser — first word, punctuation stripped,
#: no terminator logic — on the reasoning that "a second opinion derived the same way is not
#: a second opinion". That reasoning confused *deriving the answer differently* with *giving
#: a different answer*, and the coarse rule failed CI on trailer forms `CONTRIBUTING.md`
#: documents as correct, in both directions:
#:
#: - `no rows move` and `none, except two rows that merge` are filed as real changes **by
#:   design** (`CONTRIBUTING.md`: "Both err toward listing a change rather than hiding one"),
#:   and the coarse rule called them declines — so the audit failed, and no trailer text
#:   could satisfy both halves at once.
#: - `none.Tests only.` — a missing space after the terminator — went the other way:
#:   `split()[0]` yields `none.Tests`, which strips to neither a decline nor a description,
#:   so the checker called it a decline and the audit called it a move.
#:
#: Either direction turns `Changelog grouping audit` red for **every** open PR once such a
#: commit is on `main`, until the next release tag, and no PR author can fix it.
#:
#: What stays independent is the *derivation*, which is where the value always was: this
#: check renders the changelog through real git-cliff with the real config and reads the
#: result, where the checker regexes the PR body. #1555 slipped past a test that compared
#: the two halves' **vocabulary** and found them agreeing; rendering is what that test could
#: not do. Sharing no code is still pinned by
#: `test_the_rule_shares_no_vocabulary_with_the_checker`.
DECLINE_RULE = re.compile(
    r"^(?:"
    + "|".join(re.escape(word) for word in sorted(DECLINE_WORDS, key=len, reverse=True))
    + r")[ \t]*(?:["
    + re.escape(DECLINE_TERMINATORS)
    + r"].*)?$",
    re.IGNORECASE | re.DOTALL,
)

#: A body template that renders one commit id per line under its group heading. Rendering
#: ids rather than messages makes the mapping back to commits exact, and keeps this check
#: independent of message formatting, link parsers and preprocessors.
#:
#: `### {{ group | upper_first }}` reproduces how release-plz renders a heading today --
#: with no `striptags`, which is why the ordering prefix reaches the file at all (#1558).
BODY_TEMPLATE = """
{% for group, commits in commits | group_by(attribute="group") %}
### {{ group | upper_first }}
{% for commit in commits %}{{ commit.id }}
{% endfor %}{% endfor %}
"""


def load_changelog_config(path: Path) -> dict[str, Any]:
    """Return the `[changelog]` table of a release-plz config."""
    with path.open("rb") as handle:
        config: dict[str, Any] = tomllib.load(handle)["changelog"]
    return config


def _toml_value(value: object) -> str:
    """Render `value` as a TOML literal.

    Written out rather than taken from a library because the standard library ships a TOML
    reader and no writer, and the alternative is adding a dependency to a CI-only script.
    """
    if isinstance(value, str):
        escaped = (
            value.replace("\\", "\\\\")
            .replace('"', '\\"')
            .replace("\n", "\\n")
            .replace("\r", "\\r")
        )
        return f'"{escaped}"'
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, dict):
        return "{ " + ", ".join(f"{k} = {_toml_value(v)}" for k, v in value.items()) + " }"
    if isinstance(value, list):
        return "[\n  " + ",\n  ".join(_toml_value(item) for item in value) + ",\n]"
    raise TypeError(f"cannot render {type(value).__name__} as TOML")


def write_cliff_config(changelog: dict[str, Any], destination: Path) -> None:
    """Write a git-cliff config carrying the real parsers, preprocessors and postprocessors.

    Only the fields that decide *grouping* are copied. The body is this module's own, so
    that the rendered output is a list of commit ids rather than prose.

    `commit_preprocessors` is copied because grouping now depends on it: a declining
    `Representation-Change:` trailer is neutralized there (its key is renamed) so the commit
    falls through to its conventional-type group instead of the "Representation changes"
    section (#1557). Rendering without it would file every decline as a real change and the
    audit would report a decline in the section it is meant to keep clean.
    """
    lines = ["[changelog]", 'header = ""', "trim = true", f"body = {_toml_value(BODY_TEMPLATE)}"]
    if "postprocessors" in changelog:
        lines.append(f"postprocessors = {_toml_value(changelog['postprocessors'])}")
    lines += [
        "",
        "[git]",
        "conventional_commits = true",
        "filter_unconventional = false",
    ]
    if "commit_preprocessors" in changelog:
        lines.append(f"commit_preprocessors = {_toml_value(changelog['commit_preprocessors'])}")
    lines.append(f"commit_parsers = {_toml_value(changelog['commit_parsers'])}")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render(commit_range: str, config: Path, repo: Path) -> str:
    """Render the changelog for `commit_range` and return it verbatim."""
    result = subprocess.run(
        ["git-cliff", "--config", str(config), "--tag", "unreleased", commit_range],
        cwd=repo,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"git-cliff failed:\n{result.stderr.strip()}")
    return result.stdout


def parse_rendered_groups(rendered: str) -> dict[str, list[str]]:
    """Return `{group name: [commit id, ...]}` from a rendered changelog.

    The key is the heading with any `<!-- N -->` ordering prefix removed, so that a broken
    postprocessor does not also break the lookup. Keeping them coupled made one failure
    hide the other: with #1558 reverted, every heading kept its prefix, the section lookup
    found nothing, and the check reported the two real disclosures as *unfiled* rather
    than reporting the eight declines that were filed. Each failure now names itself.

    Group headings are `###`; the `##` above them is the version, which is why the pattern
    requires three hashes rather than one.
    """
    groups: dict[str, list[str]] = {}
    current: str | None = None
    for line in rendered.splitlines():
        heading = re.match(r"^#{3,}\s+(.*\S)\s*$", line)
        if heading:
            current = strip_ordering_prefix(heading.group(1))
            groups.setdefault(current, [])
        elif re.fullmatch(r"[0-9a-f]{40}", line.strip()) and current is not None:
            groups[current].append(line.strip())
    return groups


def rendered_headings(rendered: str) -> list[str]:
    """Return the heading lines verbatim, prefix included."""
    return [m.group(1) for m in re.finditer(r"^#{3,}\s+(.*\S)\s*$", rendered, re.MULTILINE)]


def strip_ordering_prefix(heading: str) -> str:
    """Remove a `<!-- N -->` group-ordering prefix from `heading`."""
    return re.sub(r"^<!--\s*\d+\s*-->\s*", "", heading).strip()


def commit_bodies(commit_range: str, repo: Path) -> dict[str, str]:
    """Return `{commit id: full message}` for every commit in `commit_range`.

    Framed on NUL, by git itself. An earlier revision appended `\\x1e` after `%B` and split
    on it, which a commit message is free to contain — and a message that did would
    re-frame the records and silently mis-attribute a trailer to the wrong commit, in a
    guard whose whole job is attributing trailers. NUL is the one byte a commit message
    cannot hold, so `-z` (record terminator) plus `%x00` (field separator) makes the
    framing git's problem rather than ours: the output is a flat NUL-separated
    `id, message, id, message, …`.
    """
    result = subprocess.run(
        ["git", "log", "-z", "--format=%H%x00%B", commit_range],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )
    # Trailing terminator yields one empty tail field; every real record is a pair.
    fields = result.stdout.split("\x00")
    if fields and fields[-1] == "":
        fields.pop()
    if len(fields) % 2 != 0:
        raise RuntimeError(
            f"git log -z returned {len(fields)} NUL-separated fields for {commit_range!r}, "
            "which is not an id/message pairing; refusing to guess the framing"
        )
    return {fields[i]: fields[i + 1].strip("\n") for i in range(0, len(fields), 2)}


def trailer_value(message: str) -> str | None:
    """Return the `Representation-Change:` value in `message`, or None.

    **Column 0, matching `scripts/check_representation_change.py` (#1573).** This pattern
    allowed `^[ \\t]*` before the token, which git and git-cliff both treat as a
    *continuation* of the value above it — so the two halves counted different trailers. A PR
    description that quotes the declining example above its own real disclosure was the live
    divergence: the checker read the column-0 disclosure and this read the indented `none`,
    then reported the commit as a decline filed under **Representation changes** and failed
    the build.
    """
    match = re.search(
        r"^Representation-Change:[ \t]*(?P<value>\S.*?)[ \t]*$",
        message,
        re.IGNORECASE | re.MULTILINE,
    )
    return match.group("value") if match else None


def opens_with_a_decline(value: str) -> bool:
    """Return whether `value`'s verdict means "nothing moved".

    The verdict is the first word; a terminator from `DECLINE_TERMINATORS` may introduce a
    reason after it. Anything else following the word — another word, or a comma — makes it a
    description of a move. See `DECLINE_RULE` for why this agrees with the checker rather
    than being deliberately coarser.
    """
    return DECLINE_RULE.match(value.strip()) is not None


def check(commit_range: str, repo: Path) -> tuple[bool, list[str]]:
    """Render the changelog and audit its Representation changes section.

    Returns `(ok, messages)`.
    """
    changelog = load_changelog_config(repo / "release-plz.toml")
    with tempfile.TemporaryDirectory() as tmp:
        config = Path(tmp) / "cliff.toml"
        write_cliff_config(changelog, config)
        rendered = render(commit_range, config, repo)

    groups = parse_rendered_groups(rendered)
    headings = rendered_headings(rendered)
    bodies = commit_bodies(commit_range, repo)
    problems: list[str] = []
    notes: list[str] = []

    listed = groups.get(REPRESENTATION_GROUP, [])
    notes.append(
        f"{len(bodies)} commits in {commit_range}; {len(listed)} listed as representation changes"
    )

    # 1. Nothing that says "no" may be filed as a change. This is #1526 and #1555.
    for commit_id in listed:
        value = trailer_value(bodies.get(commit_id, ""))
        if value is None:
            problems.append(
                f"{commit_id[:8]} is filed under {REPRESENTATION_GROUP} with no trailer at all"
            )
        elif opens_with_a_decline(value):
            problems.append(
                f"{commit_id[:8]} is filed under {REPRESENTATION_GROUP} but its trailer opens "
                f"with a decline: {value[:80]!r}"
            )

    # 2. Nothing that declares a move may be filed anywhere else. The converse failure --
    #    a real disclosure hidden in `Other` -- is the silent one, and is what an
    #    over-eager decline rule would cause.
    for commit_id, message in bodies.items():
        value = trailer_value(message)
        if value is None or opens_with_a_decline(value) or commit_id in listed:
            continue
        if commit_id in {c for ids in groups.values() for c in ids}:
            problems.append(
                f"{commit_id[:8]} declares a move but is not filed under "
                f"{REPRESENTATION_GROUP}: {value[:80]!r}"
            )
        else:
            # release-plz drops commits that touch no packaged file (`.github/` and the
            # rest of Cargo.toml's `exclude`), so absence is not necessarily an error.
            notes.append(
                f"{commit_id[:8]} declares a move and is absent from the changelog "
                f"(expected only if it touches no packaged file): {value[:60]!r}"
            )

    # 3. The ordering prefix must not survive into a heading (#1558).
    for heading in headings:
        if re.search(r"<!--\s*\d+\s*-->", heading):
            problems.append(
                f"heading {heading!r} still carries its ordering prefix; the postprocessor "
                "in release-plz.toml no longer matches what release-plz renders"
            )

    return not problems, notes + problems


def latest_tag(repo: Path) -> str:
    """Return the most recent tag reachable from HEAD.

    Raises `RuntimeError` naming the fix rather than surfacing a `CalledProcessError`
    traceback: on a shallow clone or a history with no tags, `describe` exits 128, and the
    workflow already guards its own `describe` against exactly that.
    """
    result = subprocess.run(
        ["git", "describe", "--tags", "--abbrev=0"],
        cwd=repo,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "no tag is reachable from HEAD, so the audited range cannot be derived: "
            f"`git describe --tags --abbrev=0` exited {result.returncode} "
            f"({result.stderr.strip()}). Pass an explicit `--range`, or fetch tags "
            "(`fetch-depth: 0`) if this is a shallow clone."
        )
    return result.stdout.strip()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--range",
        dest="commit_range",
        help="Commit range to audit. Defaults to `<latest tag>..HEAD`.",
    )
    parser.add_argument(
        "--repo",
        type=Path,
        default=REPO_ROOT,
        help="Repository root holding release-plz.toml.",
    )
    args = parser.parse_args(argv)

    commit_range = args.commit_range or f"{latest_tag(args.repo)}..HEAD"
    try:
        ok, messages = check(commit_range, args.repo)
    except RuntimeError as error:
        print(str(error), file=sys.stderr)
        return 1

    for message in messages:
        print(message, file=sys.stdout if ok else sys.stderr)
    if ok:
        print("Representation changes section is consistent with the trailers.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
