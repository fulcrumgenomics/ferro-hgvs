#!/usr/bin/env python3
"""
Require a `Representation-Change:` trailer on changes that can move normalized output.

Normalized output is a shipped key: a downstream consumer stores read counts against
the normalized HGVS string, so a canonical form that moves between releases re-buckets
what they have already stored. `release-plz.toml` renders a **Representation changes**
changelog section from a `Representation-Change:` trailer on the squash commit, and the
release PR body carries a reviewer checklist for it (#1433).

That machinery works. Nothing was feeding it: across the v0.13.0 cycle, 87 commits
landed and **zero** carried the trailer, while 46 touched the four directories watched at
the time (51 touch the six watched now). The section rendered empty, and an empty section
looks identical whether it is empty because nothing moved or because nobody declared
anything — so it went unnoticed for ~90 commits and the disclosure had to be reconstructed
by hand before release (#1522).

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
#: on another axis, and a consumer keys on those renderings too. `reference` decides which
#: bases any of that resolves against, and `error_handling` decides the accept/reject
#: boundary — an input that starts producing a string, or stops, is as visible to a
#: consumer as one whose string changes.
#:
#: **The last two were added on measurement, not on reasoning** (#1853, #1522). The one
#: cycle in which movement was measured without relying on anyone declaring it is v0.13.0:
#: 87 commits landed carrying **zero** trailers (`git rev-list v0.12.0..v0.13.0` counts 88,
#: the extra being release-plz's own release commit), and the **Representation changes** section
#: was reconstructed by hand by normalizing 5,761,302 expressions through both releases and
#: diffing row by row. Its 326,404 newly-normalizing rows attribute to exactly two PRs —
#: #1490, which touches `src/reference/` and nothing else, and #1501, which touches
#: `src/error_handling/` and nothing else. Neither was watched. Merged history since the
#: gate landed agrees: the only real out-of-gate disclosures are #1594 and #1661
#: (`src/error_handling/`, the second a previously-*accepted* string, so a genuine
#: migration) and #1724 (`src/reference/`).
#:
#: **Two directories are deliberately absent, and both get re-proposed.** `src/conformance/`
#: is measurement and adjudication code that by construction cannot move ferro's output, so
#: gating it would demand a declaration about something the directory cannot do — both of
#: its real disclosures (#1623, #1839) open `0 rows move in this PR`. `src/data/` was
#: proposed alongside the two above and declined on **zero marginal catch in both measured
#: populations**: its one real disclosure (#1737) also touches `src/reference/`, so the
#: entry above already gates it; both its merged gate-era commits (#1786, #1770) declined;
#: and it was touched by 0 of the commits in the v0.13.0 blind cycle. The single fact
#: that would earn it a place is a real, non-declining disclosure touching `src/data/` and
#: **not** `src/reference/` — one PR, no argument needed.
#:
#: Widening this tuple adds no required *context* — `Representation change declared` already
#: runs on every PR — so nothing is stranded on "Expected — waiting for status to be
#: reported"; only the verdict changes.
WATCHED_PREFIXES: tuple[str, ...] = (
    "src/normalize/",
    "src/hgvs/",
    "src/spdi/",
    "src/project/",
    "src/reference/",
    "src/error_handling/",
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

#: Opens or closes a Markdown code fence: three or more backticks or tildes, indented by at
#: most three spaces (CommonMark's limit before it becomes an indented code block).
FENCE_RE = re.compile(r"^ {0,3}(?P<fence>`{3,}|~{3,})(?P<info>.*)$")


def fence_match(line: str) -> re.Match[str] | None:
    """Match a fence line, honouring CommonMark's backtick-info restriction.

    A backtick fence's info string may not itself contain a backtick — otherwise an inline
    code span could open a block. Tilde fences carry no such restriction. `FENCE_RE` alone
    accepts both, so ```` ```a`b ```` reads as an opener here while CommonMark (and GitHub,
    and git-cliff's eventual reader) treat it as ordinary text.

    Getting this wrong fails closed rather than open — a spurious opener blanks the rest of
    the body, so a real declaration goes unseen and the author is refused rather than
    credited with a disclosure they did not write. That is the safe direction, and it is
    still the wrong answer.
    """
    match = FENCE_RE.match(line)
    if match is None:
        return None
    if match.group("fence")[0] == "`" and "`" in match.group("info"):
        return None
    return match


def fenced_line_numbers(text: str) -> set[int]:
    """Return the 1-based line numbers inside a fenced code block, fences included.

    `TRAILER_RE` anchors at column 0 and knows nothing about Markdown, so a trailer quoted
    inside a fence reads as a real declaration (#1929). That is the one decoration that
    reaches column 0 with its text intact — an inline code span, `**bold**`, `> blockquote`
    and a four-space indent are all caught already, and correctly *fail* with a near-miss
    hint. This one **passed**, which is worse: the author is recorded as disclosing a
    migration they explicitly said they had not worked out.

    `CONTRIBUTING.md` publishes three such blocks, so the text that triggers it is
    copy-pasteable straight out of the contributor documentation — which is exactly how it
    would be met in practice, by someone quoting the docs to ask what to write.

    The closing fence must use the same character and be at least as long as the opener, so
    a ```` ```` ```` block may contain ```` ``` ```` lines — the shape the issue reproduced
    with. An unclosed fence runs to the end of the text, which is CommonMark's rule and the
    safe direction here: it withholds a declaration rather than inventing one.
    """
    inside: set[int] = set()
    open_fence: str | None = None
    for number, line in enumerate(text.splitlines(), 1):
        match = fence_match(line)
        if open_fence is None:
            if match is not None:
                open_fence = match.group("fence")
                inside.add(number)
            continue
        inside.add(number)
        if (
            match is not None
            and match.group("fence")[0] == open_fence[0]
            and len(match.group("fence")) >= len(open_fence)
            and not match.group("info").strip()
        ):
            open_fence = None
    return inside


def strip_fenced(text: str) -> str:
    """Blank every fenced line, preserving line count so numbering still lines up.

    Blanking rather than deleting is what lets `find_near_misses` report a fenced trailer at
    its real line number in the same pass.
    """
    fenced = fenced_line_numbers(text)
    return "\n".join(
        "" if number in fenced else line for number, line in enumerate(text.splitlines(), 1)
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
#:
#: Matched against the trailer's FIRST LINE only (see `declines`), so it carries no `DOTALL`
#: and no `MULTILINE`: the reason it captures is the remainder of that one line, and any
#: continuation below is the contradiction tripwire's business, not the decline verdict's.
#: This used to carry `DOTALL` and be matched against the whole value, which made the verdict
#: depend on terminating punctuation — a bare `none` with a continuation had nothing to feed
#: the trailing `.*`, so `$` never reached end-of-string and the value was not read as a
#: decline at all (#2031). git-cliff reads the verdict off the trailer's own line, so this
#: now agrees with it by construction rather than only when a terminator happens to be
#: present -- pinned value-by-value by `test_the_two_decline_rules_agree_value_by_value`.
DECLINE_RE = re.compile(
    r"^(?:"
    + "|".join(re.escape(value) for value in sorted(NONE_VALUES))
    + r")[ \t]*(?:["
    + re.escape(DECLINE_TERMINATORS)
    + r"].*)?$",
    re.IGNORECASE,
)


def declines(declaration: str) -> bool:
    """Return whether `declaration` is a decline rather than a description of a move.

    The verdict is read from the FIRST LINE of the value, which is what git-cliff's footer
    rule reads (it is `MULTILINE`, so its `$` anchors on the trailer's own line). Reading the
    whole value instead depended on terminating punctuation to span a continuation, so a bare
    `none` followed by a move on the next line was not read as a decline at all and the
    contradiction tripwire — gated on this verdict — never fired (#2031). The tripwire still
    scans the whole value; only the *verdict* is a first-line question, and separating the two
    is what closes the gap #2027 left.

    The verdict is the first word; anything after a terminator on that line is its reason.
    """
    first_line = declaration.strip().partition("\n")[0]
    return DECLINE_RE.match(first_line) is not None


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
#: the form `CONTRIBUTING.md` offers as the way to quantify a zero, and
#: states a pass -- fire on the 950 and fail the build. The documented form was a merge
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

    Fenced regions are blanked first (#1929), so quoting `CONTRIBUTING.md`'s own example is
    not read as a declaration. The fenced line is not merely ignored — `find_near_misses`
    reports it, because someone who pasted the documented form is asking a question rather
    than staying silent.

    **This is not the set git-cliff sees** — use `find_trailers_as_git_cliff_would` for
    that. The two must stay separate; see its docstring for why collapsing them re-opens
    the hole this fix closed.
    """
    return [match.group("value") for match in TRAILER_RE.finditer(strip_fenced(text))]


def find_trailers_as_git_cliff_would(text: str) -> list[str]:
    """Return every column-0 trailer, **fenced ones included**.

    `git_conventional` parses footers with no Markdown awareness at all — any `word:` line
    is a footer token, which is the same property that makes a prose body shred its footer
    parsing (see `inject_representation_disclosure.py`). So a trailer inside a fence is
    invisible to a reader and fully visible to the changelog.

    **The sibling scripts are deliberately NOT made fence-aware**, and that is the same
    reasoning rather than an oversight. `check_changelog_grouping.py` keeps its own copy of
    the column-0 pattern and `inject_representation_disclosure.py` imports it from there;
    both read *commit messages* and exist to check or reproduce what git-cliff actually did.
    Teaching them about fences would make them disagree with the changelog they audit, which
    is the opposite of their job. Only the PR-body checker chooses a *value*, and only that
    choice needed narrowing.

    That is why the duplicate refusal in `check` counts *this* set rather than
    `find_declarations`. Making both fence-aware looks tidier and is a safety regression:
    a body carrying a quoted example **beside** a real trailer would then pass here while
    git-cliff still saw two footers, and the quoted 577-row example could land in the
    changelog under an author who declared `none`. The refusal is what forces the author to
    indent or remove the quotation, which is the escape hatch `CONTRIBUTING.md` documents.

    So the split is deliberate and the asymmetry is the point: **fenced text may not BE the
    declaration, but it still COUNTS toward how many the changelog will see.**
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

    fenced = fenced_line_numbers(text)
    misses: list[tuple[int, str, str]] = []
    for number, line in enumerate(text.splitlines(), 1):
        if TRAILER_RE.match(line):
            # Column 0 and well formed, so the only thing keeping it from being a
            # declaration is the fence around it (#1929). Report it rather than passing over
            # it: this shape is almost always someone quoting `CONTRIBUTING.md` to ask what
            # to write, and silence would be indistinguishable from having declared nothing.
            if number in fenced:
                misses.append(
                    (
                        number,
                        line.strip(),
                        "inside a fenced code block, so it is documentation rather than a "
                        "declaration; move it to column 0 outside any ``` or ~~~ block",
                    )
                )
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


#: A line that ends a trailer's value because it opens something the value cannot own:
#: another column-0 trailer token, a Markdown heading, or an HTML comment (GitHub appends
#: CodeRabbit's summary after the trailer, opening `<!-- ... -->`).
#:
#: A deliberately SEPARATE copy of the terminator in
#: `inject_representation_disclosure.disclosure_value`. This checker shares no code with
#: `check_changelog_grouping.py` or the injector that imports from it, so the two decline
#: rules can be wrong independently (#1555). The column-0 rule this extends is this module's
#: own `TRAILER_RE`, so no fourth copy of *that* is minted.
CONTINUATION_END_RE = re.compile(r"^([A-Za-z][A-Za-z0-9-]*:[ \t]|#{1,6}\s|<!--)")


def full_declaration_value(text: str) -> str | None:
    """Return the first trailer's value with its continuation lines, or `None`.

    `find_declaration` returns only the first line, because `TRAILER_RE` stops at the
    newline. This repo's trailers are routinely multi-line with UNINDENTED continuation lines
    -- `git interpret-trailers --parse` reports no trailers at all for #1537/#1535/#1547 --
    so a decline on the first line can hide a contradicting `<n> rows move` on a later one
    (#1854). The contradicted-decline tripwire must read the whole value, not the first
    clause, or the disclosure is filed as `none` and never reaches the changelog.

    The value runs from the trailer to the next line it cannot own -- a column-0 trailer
    token, a Markdown heading, or the HTML comment GitHub appends -- the same bound
    `inject_representation_disclosure.disclosure_value` uses, arrived at independently rather
    than imported. Fenced regions are blanked first, matching `find_declarations`, so a
    continuation quoted inside a fence is not scavenged into the value.
    """
    first = find_declaration(text)
    if first is None:
        return None
    lines = strip_fenced(text).splitlines()
    start = next(index for index, line in enumerate(lines) if TRAILER_RE.match(line))
    collected = [first]
    for line in lines[start + 1 :]:
        if CONTINUATION_END_RE.match(line):
            break
        collected.append(line.rstrip())
    while collected and not collected[-1].strip():
        collected.pop()
    return "\n".join(collected)


def check(changed: list[str], declaration_text: str) -> tuple[bool, str]:
    """
    Decide whether `changed` is adequately declared by `declaration_text`.

    Returns `(ok, message)`. `ok` is False in four cases:

    1. the message carries **more than one** trailer (#1573) — checked first, and ahead of
       the watched-file test, because two trailers are a disagreement with how the
       *changelog* files the commit, which applies to every commit. Cases 1 and 2 are not
       exclusive: when every one of those trailers is fenced, case 1's refusal carries
       case 2's near-miss lines as well, because "delete the superseded trailer" names
       nothing the author wrote;
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
    # Counted as GIT-CLIFF would count them, fenced examples included: the refusal exists
    # because this checker reads the first trailer while git-cliff matches every footer, and
    # git-cliff cannot see a fence. Counting the fence-aware set here would let a quoted
    # example sit beside a real `none` and still reach the changelog (#1929).
    as_git_cliff_sees = find_trailers_as_git_cliff_would(declaration_text)
    if len(as_git_cliff_sees) > 1:
        listed = "\n".join(f"  Representation-Change: {value}" for value in as_git_cliff_sees)
        message = (
            f"{len(as_git_cliff_sees)} `Representation-Change:` trailers found:\n{listed}\n\n"
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
        # An EMPTY `declarations` here means every trailer git-cliff counted is fenced —
        # `find_declarations` reads the same pattern over `strip_fenced`, so the two sets
        # differ by exactly the fenced lines. Without this arm that author is told to
        # "keep exactly one" and to "delete the superseded trailer" about text they wrote
        # as *documentation*, and which a Markdown reader renders as code. The advice
        # happens to work — indenting inside the fence drops the line from `TRAILER_RE` —
        # but naming the wrong construct is what sends someone hunting for a declaration
        # they never made. `CONTRIBUTING.md` publishes three such blocks, so quoting two of
        # them to ask which form applies is the ordinary way into this path.
        #
        # `find_near_misses` is the half that holds the right diagnosis and it is
        # unreachable from here otherwise: the duplicate refusal returns first, and the
        # near-miss scan below stands down as soon as one readable trailer exists. Append
        # rather than fall through, because BOTH facts are true and only together do they
        # explain the refusal — the count is what git-cliff will see, and the fence is why
        # the author cannot see it.
        if not declarations:
            fenced = find_near_misses(declaration_text)
            if fenced:
                fenced_lines = "\n".join(
                    f"  line {number}: {line}\n    -> {reason}" for number, line, reason in fenced
                )
                message += (
                    "\n\nEvery one of them is inside a fenced code block, so none of them is a\n"
                    "declaration you can see -- but git-cliff has no Markdown awareness and\n"
                    f"counts all {len(as_git_cliff_sees)}:\n\n"
                    f"{fenced_lines}\n\n"
                    "So there is nothing to delete here. Indent the quoted lines inside the\n"
                    "fence to keep them out of the count, and write your real declaration at\n"
                    "column 0 outside any block."
                )
        return False, message

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

    # Two spans, deliberately different (#2031). The decline VERDICT is a first-line question
    # -- `declines` reads it off the trailer's own line, which is what git-cliff reads -- but
    # the contradiction TRIPWIRE must scan the WHOLE value, because a decline on line 1 can
    # hide a `<n> rows move` on an unindented continuation line (#1854). Reading the verdict
    # off the whole value instead (as #2027 did) made it depend on terminating punctuation, so
    # a BARE `none` with a move continuation was read as neither a decline nor a contradiction.
    # `declines(full_value)` still asks the first-line question -- `full_value`'s first line is
    # this same trailer's -- while `contradicted_decline` searches the whole value below.
    # `declarations` is non-empty here, so a full value exists.
    full_value = full_declaration_value(declaration_text)
    assert full_value is not None
    if declines(full_value):
        claim = contradicted_decline(full_value)
        if claim is not None:
            return False, (
                f"This trailer declines and then describes a move:\n\n  {full_value}\n\n"
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
            f"Declared `Representation-Change: {full_value}` over {len(watched)} watched file(s).",
        )
    return (
        True,
        f"Declared a representation change over {len(watched)} watched file(s): {full_value}",
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
