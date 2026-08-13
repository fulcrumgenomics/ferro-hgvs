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
#: an anchored `(?!0\b)` read them as nonzero and failed the decline.
#:
#: **The number that matters is the COUNT, never the denominator** (#1647). Both orders
#: occur in this repository's own trailers and both must be read the same way:
#:
#:     3 rows of 500,004 move          <- count first, denominator trailing
#:     3 of 500,004 rows move          <- count first, denominator between count and `rows`
#:     3 of 500,004 corpus rows move   <- and a noun phrase before `rows`
#:
#: This used to be written as "the count must sit immediately before `rows`", which reads
#: the second form's *denominator* as the count. That made `none. 0 of 950 rows move` --
#: the form `CONTRIBUTING.md` and `CLAUDE.md` both offer as the way to quantify a zero, and
#: both state passes -- fire on the 950 and fail the build. The documented form was a merge
#: blocker.
#:
#: Two lookbehinds carry the distinction, and the second is not redundant:
#:
#: * `(?<!of\s)` refuses to read a denominator as a count, which is what rescues
#:   `0 of 950 rows move`.
#: * `(?<![\d,])` stops the match restarting in the MIDDLE of a grouped number. Without it
#:   `0 of 78,298 rows move` still fires -- the comma opens a word boundary, so the scan
#:   resumes at `298`, which is nonzero, sits immediately before `rows`, and is not preceded
#:   by `of `. Measured; the first lookbehind alone leaves that decline failing.
#:
#: The third form is here because a version of this comment claiming "every quantified
#: disclosure in the corpus takes one of the two shapes above" was refuted by a trailer that
#: was already on `main`: #1651's `1826 of 78298 corpus rows move`. One word between the
#: denominator and `rows` and the tripwire went silent -- with the verdict flipped to a
#: decline that is a real move hidden behind `none`, the one direction this mechanism must
#: never fail in. Measured, not reasoned: the two-lookbehind pattern passed that string.
#: So up to three intervening words are allowed. Punctuation still stops the run, which is
#: what keeps `1 issue filed, 0 rows move` from reading `1` as the moving count.
#:
#: It catches the phrasing this repository actually uses -- and it is not a general
#: contradiction detector. Spelled out numbers and unquantified claims pass, as does a
#: nonzero count separated from `rows` by more than three words. It is a tripwire, not a
#: proof, and the way to find the next gap is to run a real trailer through it rather than
#: to read the pattern.
MOVEMENT_CLAIM_RE = re.compile(
    r"\b(?<!of\s)(?<![\d,])(?=[\d,]*[1-9])\d[\d,]*\s+(?:of\s+[\d,]+\s+)?"
    r"(?:[A-Za-z][\w-]*\s+){0,3}rows?"
    r"(?:\s+of\s+[\d,]+)?\s+"
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


#: A line that *looks* like a declaration but that `TRAILER_RE` will not match, and that
#: git-cliff's footer rule and `inject_representation_disclosure.py` will not match either.
#:
#: #1838 is the live instance: its description declares on line 3 as
#: `` `Representation-Change: none` ``, inside a code span. One backtick, and all three
#: consumers see silence while the author believes they declared. Nothing was lost there —
#: the PR touches no watched directory and the value was a decline, which the changelog
#: excludes anyway — but a gated change carrying a real disclosure in that shape is exactly
#: the silent representation change this whole mechanism exists to prevent.
#:
#: The pattern is deliberately looser than `TRAILER_RE` in three directions, each a real
#: Markdown habit rather than an invented one:
#:
#: * **decoration before the token** — a code span, `>`, `#`, `-`/`*`/`+`, `*`/`_`, or plain
#:   indentation. Every one of these renders the line as prose while leaving it looking like
#:   a trailer in the editor.
#: * **decoration around the colon and the value** — `**Representation-Change:** none`.
#: * **the separator** — the hyphen is what git-cliff keys on, so `Representation Change:`
#:   and `Representation_Change:` are as invisible as a backtick.
#:
#: **What keeps it from firing on prose is the anchor, and that is the load-bearing half.**
#: A declaration is a *line*, not a phrase, so the token must open the line once decoration
#: is stripped. Measured over all 39 open PRs, that alone excludes every mention that is not
#: a declaration: #1753's two reviewer-checklist mentions (the release PR) and #1752's
#: "the `Representation-Change:` disclosure above is unaffected" all carry the token
#: mid-sentence and none of them match.
#:
#: **The value requirement is defensive rather than measured-necessary.** Dropping it changes
#: nothing on today's corpus — the same one PR is reported either way, measured — but
#: `Representation-Change:` with nothing after it names the field rather than filling it in,
#: so there is no disclosure to lose and reporting it would be a false accusation. Keep it,
#: and do not cite the release PR as evidence for it; that was an earlier draft's claim and
#: it does not survive re-derivation.
NEAR_MISS_RE = re.compile(
    r"^(?P<lead>[ \t>*_+#`-]*)"
    r"Representation(?P<separator>[-_ ]?)Change"
    # Captured, not merely skipped. `Representation-Change : none` is unreadable to all
    # three consumers and matches here, but its lead is empty and its separator is a
    # hyphen — so with the gap discarded no reason could be named for it and
    # `find_near_misses` dropped it silently, which is the very shape this scan exists to
    # stop being silent about. Same for a tab, and for `Representation-Change*: none`.
    r"(?P<gap>[`*_ \t]*):"
    r"[`*_ \t]*"
    # The value's first character must be a real one. `\S` is not enough: on
    # `` `Representation-Change:` `` the regex backtracks and matches the *closing backtick*
    # as the value, so a bare field mention in a code span was reported as a failed
    # declaration. Excluding the decoration characters here is what makes "has a value" mean
    # "has something other than the punctuation that wrapped it".
    r"[^\s`*_].*?[ \t`*_]*$",
    re.IGNORECASE,
)

#: A list marker is `-`, `*` or `+` followed by whitespace; the same characters *not*
#: followed by whitespace are emphasis. That is Markdown's own rule, and keeping it lets the
#: message say `list marker` or `emphasis` rather than naming a character and leaving the
#: author to work out which construct they wrote.
_LIST_MARKER_RE = re.compile(r"^[-*+][ \t]")


def _near_miss_reason(line: str, lead: str, separator: str, gap: str) -> str:
    """Name what makes `line` unreadable to the strict parser, in the author's terms.

    Several can apply at once — an indented code span inside a list item is one line — so
    every one that applies is reported. The point is that the author can see which character
    to delete; "no declaration found" is what this replaces, and it sends them looking for a
    trailer they have already written.

    **Never returns the empty string for a line `NEAR_MISS_RE` matched**, and that is a
    contract rather than an observation: the caller reports what this names, so a match with
    nothing to say is a near miss dropped in silence. It is enforced by
    `test_every_near_miss_shape_names_a_reason`, which drives every shape through the real
    pattern rather than through a list of the ones anybody thought of.
    """
    reasons: list[str] = []
    if "`" in line:
        reasons.append("wrapped in a code span (a backtick)")
    marker = lead.lstrip(" \t")
    if marker.startswith(">"):
        reasons.append("prefixed by a block quote marker `>`")
    elif _LIST_MARKER_RE.match(marker):
        reasons.append(f"prefixed by a list marker `{marker[0]}`")
    elif marker.startswith("#"):
        reasons.append("prefixed by a heading marker `#`")
    elif marker.startswith(("*", "_")):
        reasons.append(f"wrapped in Markdown emphasis (`{marker[0]}`)")
    elif not marker and lead:
        # Only when the lead is whitespace and nothing else. Testing `lead` alone reported a
        # backtick-led line as "indented", which names a character the author did not type
        # and sends them to delete whitespace that is not the problem.
        reasons.append("indented, which git reads as a continuation of the line above")
    if separator != "-":
        spelt = "a space" if separator == " " else ("an underscore" if separator else "nothing")
        reasons.append(f"spelled with {spelt} instead of a hyphen (the spelling git-cliff keys on)")
    decoration = gap.strip(" \t")
    if decoration and "`" not in decoration:
        # A backtick anywhere on the line is already reported above, so naming it again here
        # would read as two faults where the author typed one.
        reasons.append(f"separated from its colon by `{decoration[0]}`")
    elif gap and not decoration:
        reasons.append("separated from its colon by whitespace")
    if not reasons:
        # Unreached today, and deliberately not an assertion: a matched line the specific
        # rules cannot explain is still unreadable to every consumer, so the honest answer is
        # to report it in general terms rather than to drop it or to crash on it.
        reasons.append("not the exact `Representation-Change:` token at the start of the line")
    return "; ".join(reasons)


def find_near_misses(text: str) -> list[tuple[int, str, str]]:
    """Return `(line number, line, reason)` for every near-miss declaration in `text`.

    Lines `TRAILER_RE` already accepts are skipped — those are real trailers, not near
    misses. Nothing else is: `_near_miss_reason` is contracted to name a reason for every
    line the pattern matched, so a match can no longer be discarded for having nothing to
    say. That discard used to be here, described as unreachable, and it was not —
    `Representation-Change : none`, the tab form, and `Representation-Change*: none` all
    matched, named no reason, and were dropped, which reported them as silence in a scan
    whose whole subject is a declaration nobody reads.

    **A single valid trailer anywhere in `text` silences the whole scan**, and that is a
    property of this function rather than a rule its callers have to remember. An author who
    declared at column 0 has demonstrated they know the form, and the same author routinely
    *discusses* that declaration in prose: #1742 and #1837 each quote their own trailer in a
    code span, on a line beside a real one. `CONTRIBUTING.md` also documents indenting a
    quoted example as the way to keep it out of the trailer count, and the two-trailer
    refusal's own message recommends exactly that — so reporting here would contradict the
    escape hatch this checker already offers.

    Re-measured over all **46** open PRs later on 2026-08-13: the pattern matches a line in
    **five** of them, **two** of which carry a real trailer as well (#1742 and #1837, each
    explaining its own declaration in prose — failing those would be the #1555 mistake in a
    new place), so the scan reports **three**: #1838, the PR that motivated it, plus #1850
    and #1855, which opened afterwards in the same shape.

    Quote the command rather than the figure — the denominator is a property of the day, and
    it moved from 39 to 46 while this PR was in review:

    ```
    gh pr list -R fulcrumgenomics/ferro-hgvs --state open --limit 100 --json number,body
    ```

    Three live instances in one afternoon is the argument for the check. One would not have
    been.
    """
    if find_declarations(text):
        return []

    misses: list[tuple[int, str, str]] = []
    for number, line in enumerate(text.splitlines(), 1):
        if TRAILER_RE.match(line):
            continue
        match = NEAR_MISS_RE.match(line)
        if match is None:
            continue
        reason = _near_miss_reason(
            line, match.group("lead"), match.group("separator"), match.group("gap")
        )
        misses.append((number, line.strip(), reason))
    return misses


def find_declaration(text: str) -> str | None:
    """Return the first trailer's value, or `None` when the text carries no trailer."""
    declarations = find_declarations(text)
    return declarations[0] if declarations else None


def check(changed: list[str], declaration_text: str) -> tuple[bool, str]:
    """
    Decide whether `changed` is adequately declared by `declaration_text`.

    Returns `(ok, message)`. `ok` is False in four cases:

    1. the message carries **more than one** trailer (#1573) — checked first, and ahead of
       the watched-file test, because two trailers are a disagreement with how the
       *changelog* files the commit, which applies to every commit;
    2. the message carries a **near miss** and no readable trailer — a line that looks like
       a declaration but that none of the three consumers can see. Also ahead of the
       watched-file test, and for the same reason as case 1;
    3. a watched file changed and **no** trailer is present at all;
    4. the trailer **declines while describing a move** (`no. 3 rows move`).

    A trailer whose value is `none` otherwise passes, because declining is a declaration.

    Note case 4 is scoped to watched changes, unlike cases 1 and 2: it sits after the
    watched-file early return, so a docs-only commit whose voluntary trailer contradicts
    itself is not refused. That is the narrower choice — the contradiction is the same
    either way — and it keeps a gratuitous trailer on an unwatched change from failing CI.

    Cases 2 and 4 point in opposite directions on that question and both are deliberate. A
    contradicted decline is a trailer the machinery *does* read, so the author's declaration
    is published and merely disagrees with itself; a near miss is a declaration nobody reads
    at all, which is indistinguishable from silence to every consumer and to the author is
    indistinguishable from having declared. Only the second can ship a representation change
    with nothing in the changelog, so only the second is refused everywhere.
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

    # Also refused BEFORE the watched-file test, and for the same reason as the case above:
    # the harm is in what the *changelog* renders, which applies to every commit. Putting it
    # behind the early return would rebuild the exact hole it closes — #1838, the PR that
    # motivated this, touches no watched directory, so it reaches that return before the
    # unreadable trailer would have mattered and passes green today.
    #
    # `find_near_misses` returns nothing when a readable trailer exists, so there is no
    # condition to repeat here — see its docstring for why a valid column-0 trailer silences
    # the scan rather than merely outranking it.
    near_misses = find_near_misses(declaration_text)
    if near_misses:
        listed = "\n".join(
            f"  line {number}: {line}\n    -> {reason}" for number, line, reason in near_misses
        )
        return False, (
            f"{len(near_misses)} line(s) look like a `Representation-Change:` declaration but\n"
            "no consumer can read them:\n\n"
            f"{listed}\n\n"
            "This is not a style note. The trailer must begin at column 0 with the exact\n"
            "token, because three separate consumers match it that way and all three see\n"
            "silence here: this check; git-cliff's footer rule in `release-plz.toml`, so it\n"
            "never reaches the changelog; and `scripts/inject_representation_disclosure.py`,\n"
            "so its text is never attached to the changelog bullet.\n\n"
            "You have declared something and nobody will be told. Write the line with no\n"
            "decoration of any kind:\n\n"
            "  Representation-Change: none\n\n"
            "or, for a real disclosure, what moved and whether the affected inputs were\n"
            "previously rejected (free) or accepted (a migration):\n\n"
            "  Representation-Change: 3 rows of 500,004 move, 2 respell / 1 merge.\n"
            "    Previously-accepted inputs, so a real migration.\n\n"
            "To quote the form as an *example* rather than declare it, indent the line — an\n"
            "indented line is a continuation, which is what keeps it out of the count. Note\n"
            "that indenting is a quoting device only *beside* a real declaration: this scan\n"
            "stands down as soon as one readable trailer exists, so a body whose only\n"
            "`Representation-Change:` line is the indented one is reported here rather than\n"
            "read as a quotation."
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
