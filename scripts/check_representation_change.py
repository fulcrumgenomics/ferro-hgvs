#!/usr/bin/env python3
"""
Require a `Representation-Change:` trailer on changes that can move normalized output.

Normalized output is a shipped key: a downstream consumer stores read counts against
the normalized HGVS string, so a canonical form that moves between releases re-buckets
what they have already stored. `release-plz.toml` renders a **Representation changes**
changelog section from a `Representation-Change:` trailer on the squash commit, and the
release PR body carries a reviewer checklist for it (#1433).

That machinery works. Nothing was feeding it: across the v0.13.0 cycle, 87 commits
landed and **zero** carried the trailer, while 46 touched the directories below. The
section rendered empty, and an empty section looks identical whether it is empty because
nothing moved or because nobody declared anything — so it went unnoticed for ~90 commits
and the disclosure had to be reconstructed by hand before release (#1522).

This closes that loop: a change under a watched directory must *say* whether it moves
output. Saying "no" is a first-class answer — `Representation-Change: none` — because the
point is to force the judgement, not to force a disclosure. What is rejected is silence.

Usage:
    python scripts/check_representation_change.py \\
        --changed-files changed.txt --declaration-file body.txt

    # or with the declaration on stdin
    git diff --name-only origin/main...HEAD > changed.txt
    python scripts/check_representation_change.py --changed-files changed.txt -

Exits 0 when the change is clear, 1 when a declaration is required and absent.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

#: Directories whose output a change can move. A variant's canonical form is decided in
#: `normalize`; `hgvs` owns how it is spelled back out; `spdi` and `project` re-express it
#: on another axis, and a consumer keys on those renderings too.
WATCHED_PREFIXES: tuple[str, ...] = (
    "src/normalize/",
    "src/hgvs/",
    "src/spdi/",
    "src/project/",
)

#: The trailer git-cliff groups on. Matched at the start of a line, case-insensitively,
#: because a trailer is prose written by a human and `representation-change:` is the same
#: declaration as `Representation-Change:`.
#:
#: The value must open with a non-space and the surrounding whitespace classes are
#: horizontal-only. Both matter: `\s*(?P<value>.+?)\s*$` under `MULTILINE` matches a lone
#: space as the value, so a bare `Representation-Change:` with nothing after it counted as
#: a declaration — the precise hole this check exists to close. And `\s` spans newlines, so
#: a permissive class would let the value be scavenged from the *following* line.
#:
#: **The token must sit at column 0** (#1573). The pattern used to allow `[ \t]*` before it,
#: which made this checker recognise as a *trailer* what git and git-cliff both treat as a
#: *continuation* of the value above it — and `CONTRIBUTING.md` documents indentation as
#: exactly that continuation mechanism. The disagreement is only reachable when a
#: continuation line itself opens with the token, i.e. when a PR description quotes a
#: declining example under its own real disclosure, but there it is the whole defect: the
#: quoted example became a second declaration here and stayed prose for git-cliff. Requiring
#: column 0 makes both halves count the same trailers, which is what lets `check` refuse a
#: genuine duplicate and pass an indented quotation. Measured over every open PR at the time
#: of the change: no description carried an indented trailer, so nothing was reclassified.
TRAILER_RE = re.compile(
    r"^Representation-Change:[ \t]*(?P<value>\S.*?)[ \t]*$",
    re.IGNORECASE | re.MULTILINE,
)

#: Values that declare "this moves nothing". Anything else is read as a description of
#: what moved, which is what lands in the changelog.
NONE_VALUES = frozenset({"none", "no", "n/a", "na"})

#: Punctuation that may separate a decline from the reason for it, so that
#: `none. Tests only; no src/ file is touched` still reads as a decline (#1555).
#:
#: Requiring a bare verdict punished the contributors who explained themselves: in the
#: v0.13.1 cycle 8 commits declined *with a reason* and every one of them was grouped as a
#: real representation change. Explaining a decline is the habit to encourage.
#:
#: `,` is deliberately absent. A comma usually introduces a qualification that changes the
#: verdict -- "none, except two rows" -- where a full stop, semicolon, colon or dash closes
#: it. This set is mirrored in `release-plz.toml`'s exclusion rule and pinned by
#: `test_the_two_decline_rules_agree_value_by_value`.
DECLINE_TERMINATORS = ".;:—–"

#: A declination: one of `NONE_VALUES`, alone or followed by a reason. Built from the two
#: constants above so the vocabulary has one definition rather than a regex restating it.
#: `DOTALL` so a reason spanning lines does not defeat the trailing `.*`.
DECLINE_RE = re.compile(
    r"^(?:"
    + "|".join(re.escape(value) for value in sorted(NONE_VALUES))
    + r")[ \t]*(?:["
    + re.escape(DECLINE_TERMINATORS)
    + r"].*)?$",
    re.IGNORECASE | re.DOTALL,
)


def declines(declaration: str) -> bool:
    """Return whether `declaration` is a decline rather than a description of a move.

    The verdict is the first word; anything after a terminator is the reason for it.
    """
    return DECLINE_RE.match(declaration.strip()) is not None


#: A quantified movement claim -- "3 rows move", "577 rows merge", "2 rows respell",
#: "3 rows of 500,004 move". Reading the verdict from the first word means a trailer can
#: contradict itself in its own reason (`no. 3 rows move`), and that failure is silent and
#: in the dangerous direction: the disclosure disappears from the changelog and nobody is
#: told. This is the tripwire for it.
#:
#: Two deliberate limits, because a looser rule fails honest PRs. **Zero is excluded**: a
#: good decline often quantifies its zero, and `0 rows move over 5,761,302 real
#: expressions` (#1535) must not trip. Excluded by requiring a nonzero digit *somewhere in
#: the count* rather than by rejecting a bare `0`, so `0,000` and `000` are zero here too —
#: an anchored `(?!0\b)` read them as nonzero and failed the decline. **The count must sit
#: immediately before `rows`**,
#: which is what keeps `0 of 950 real cis-allele rows move their normalized string`
#: (#1547) from tripping on the 950 -- a looser pattern fires on that exemplary decline,
#: measured.
#:
#: It catches the phrasing this repository actually uses -- every quantified disclosure in
#: the corpus takes this shape -- and it is not a general contradiction detector. Spelled
#: out numbers and unquantified claims pass. It is a tripwire, not a proof.
MOVEMENT_CLAIM_RE = re.compile(
    r"\b(?=[\d,]*[1-9])\d[\d,]*\s+rows?(?:\s+of\s+[\d,]+)?\s+"
    r"(?:move|moves|moved|merge|merges|merged|split|splits|"
    r"respell|respells|respelled)\b",
    re.IGNORECASE,
)


def contradicted_decline(declaration: str) -> str | None:
    """Return the movement claim a declining `declaration` makes, if it makes one."""
    if not declines(declaration):
        return None
    match = MOVEMENT_CLAIM_RE.search(declaration)
    return match.group(0) if match else None


def watched_files(changed: list[str]) -> list[str]:
    """Return the changed paths that sit under a watched directory."""
    return [p for p in changed if any(p.startswith(prefix) for prefix in WATCHED_PREFIXES)]


def find_declarations(text: str) -> list[str]:
    """Return every trailer value in `text`, in the order they appear.

    Plural because *how many* there are is itself a verdict — see `check`. One trailer is
    the only shape whose meaning both this checker and `release-plz.toml` agree on.
    """
    return [match.group("value") for match in TRAILER_RE.finditer(text)]


def find_declaration(text: str) -> str | None:
    """Return the first trailer's value, or `None` when the text carries no trailer."""
    declarations = find_declarations(text)
    return declarations[0] if declarations else None


def check(changed: list[str], declaration_text: str) -> tuple[bool, str]:
    """
    Decide whether `changed` is adequately declared by `declaration_text`.

    Returns `(ok, message)`. `ok` is False in three cases:

    1. the message carries **more than one** trailer (#1573) — checked first, and ahead of
       the watched-file test, because two trailers are a disagreement with how the
       *changelog* files the commit, which applies to every commit;
    2. a watched file changed and **no** trailer is present at all;
    3. the trailer **declines while describing a move** (`no. 3 rows move`).

    A trailer whose value is `none` otherwise passes, because declining is a declaration.

    Note case 3 is scoped to watched changes, unlike case 1: it sits after the
    watched-file early return, so a docs-only commit whose voluntary trailer contradicts
    itself is not refused. That is the narrower choice — the contradiction is the same
    either way — and it keeps a gratuitous trailer on an unwatched change from failing CI.
    """
    # Refused BEFORE the watched-file test, and deliberately: the harm is in changelog
    # grouping, which applies to every commit, not only to the ones touching a watched
    # directory. This job already runs on every PR.
    declarations = find_declarations(declaration_text)
    if len(declarations) > 1:
        listed = "\n".join(f"  Representation-Change: {value}" for value in declarations)
        return False, (
            f"{len(declarations)} `Representation-Change:` trailers found:\n{listed}\n\n"
            "Keep exactly one. Two trailers do not mean anything consistent, because this\n"
            "check and the changelog resolve them differently: this check reads the FIRST\n"
            "trailer, while git-cliff matches its footer rule against EVERY footer and takes\n"
            "the first rule that matches any of them -- so a decline anywhere wins there,\n"
            "whatever its position. A message disclosing a move and then declining one\n"
            "therefore passes here as a disclosure and is filed under `Other` in the\n"
            "changelog, and the disclosure silently disappears.\n\n"
            "Delete the superseded trailer, or indent it if it is quoted as an example --\n"
            "an indented line is a continuation of the value above it, not a new trailer."
        )

    watched = watched_files(changed)
    if not watched:
        return True, "No file under a watched directory changed; no declaration needed."

    declaration = declarations[0] if declarations else None
    if declaration is None:
        listed = "\n".join(f"  {p}" for p in sorted(watched))
        return False, (
            f"{len(watched)} changed file(s) can move normalized output:\n{listed}\n\n"
            "Add a `Representation-Change:` trailer to the PR description, which is what\n"
            "becomes the squash commit body and what release-plz reads.\n\n"
            "  Representation-Change: <what moved, in which direction, over how many rows,\n"
            "                          and whether the affected inputs were previously\n"
            "                          rejected (free) or accepted (a migration)>\n\n"
            "If the change cannot move output — a comment, a doc, a test, a rename — say so\n"
            "explicitly:\n\n"
            "  Representation-Change: none\n\n"
            "Declining is a declaration and passes this check. Silence does not, because an\n"
            "absent trailer is indistinguishable from an unconsidered one.\n\n"
            "Measure with: cargo run --release --features dev --example dump_normalized_corpus"
        )

    if declines(declaration):
        claim = contradicted_decline(declaration)
        if claim is not None:
            return False, (
                f"This trailer declines and then describes a move:\n\n  {declaration}\n\n"
                f"The verdict is the first word, so this is filed as `none` and the "
                f"disclosure -- {claim!r} -- never reaches the changelog. Nobody is told, "
                "which is the failure this whole mechanism exists to prevent.\n\n"
                "Say which it is. If the change moves output, lead with what moved:\n\n"
                "  Representation-Change: 3 rows of 500,004 move, 2 respell / 1 merge,\n"
                "    all away from the already-shipped form. Previously-accepted inputs,\n"
                "    so a real migration.\n\n"
                "If it moves nothing and the number describes something else -- a corpus "
                "that grew, a measurement quoted from another change -- say so without the "
                "`N rows move` phrasing, which reads as this change's own disclosure."
            )
        return (
            True,
            f"Declared `Representation-Change: {declaration}` over {len(watched)} watched file(s).",
        )
    return (
        True,
        f"Declared a representation change over {len(watched)} watched file(s): {declaration}",
    )


def read_text(path: str) -> str:
    """Read `path`, or stdin when it is `-`."""
    return sys.stdin.read() if path == "-" else Path(path).read_text(encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--changed-files",
        required=True,
        help="File holding one changed path per line, or `-` for stdin.",
    )
    parser.add_argument(
        "--declaration-file",
        required=True,
        help="File holding the PR description / squash commit body, or `-` for stdin.",
    )
    args = parser.parse_args(argv)
    if args.changed_files == "-" and args.declaration_file == "-":
        parser.error("--changed-files and --declaration-file cannot both read stdin")

    changed = [line.strip() for line in read_text(args.changed_files).splitlines() if line.strip()]
    ok, message = check(changed, read_text(args.declaration_file))
    print(message, file=sys.stdout if ok else sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
