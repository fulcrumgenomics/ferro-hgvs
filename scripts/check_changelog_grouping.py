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

#: Words that mean "this moves nothing", as a *first word*. Deliberately a plain prefix
#: test with no terminator logic: `scripts/check_representation_change.py` owns the precise
#: rule, and a second opinion derived the same way is not a second opinion.
DECLINE_WORDS = frozenset({"none", "no", "n/a", "na"})

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
    """Write a git-cliff config carrying the real parsers and postprocessors.

    Only the fields that decide *grouping* are copied. The body is this module's own, so
    that the rendered output is a list of commit ids rather than prose.
    """
    lines = ["[changelog]", 'header = ""', "trim = true", f"body = {_toml_value(BODY_TEMPLATE)}"]
    if "postprocessors" in changelog:
        lines.append(f"postprocessors = {_toml_value(changelog['postprocessors'])}")
    lines += [
        "",
        "[git]",
        "conventional_commits = true",
        "filter_unconventional = false",
        f"commit_parsers = {_toml_value(changelog['commit_parsers'])}",
    ]
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
    """Return `{commit id: full message}` for every commit in `commit_range`."""
    separator = "\x1e"
    result = subprocess.run(
        ["git", "log", f"--format=%H%x00%B{separator}", commit_range],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )
    bodies = {}
    for record in result.stdout.split(separator):
        if record.strip():
            commit_id, message = record.strip("\n").split("\x00", 1)
            bodies[commit_id] = message
    return bodies


def trailer_value(message: str) -> str | None:
    """Return the `Representation-Change:` value in `message`, or None."""
    match = re.search(
        r"^[ \t]*Representation-Change:[ \t]*(?P<value>\S.*?)[ \t]*$",
        message,
        re.IGNORECASE | re.MULTILINE,
    )
    return match.group("value") if match else None


def opens_with_a_decline(value: str) -> bool:
    """Return whether `value`'s first word means "nothing moved".

    Punctuation is stripped from the first word, so `none.`, `none:` and `none —` all read
    as declines. This is the whole rule: no terminator logic, no vocabulary sharing.
    """
    first_word = value.strip().split()[0] if value.strip() else ""
    return first_word.strip(".,;:!—–").lower() in DECLINE_WORDS


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
    """Return the most recent tag reachable from HEAD."""
    result = subprocess.run(
        ["git", "describe", "--tags", "--abbrev=0"],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
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
