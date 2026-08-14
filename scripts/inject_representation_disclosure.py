#!/usr/bin/env python3
"""
Attach each `Representation-Change:` trailer's own text to its changelog bullet.

## The problem this closes

`release-plz.toml` routes trailered commits into a **Representation changes** section, but
a changelog bullet is the commit **subject** and nothing else:

    - *(normalize)* never split a delins into members on consecutive nucleotides ([#1537](...))

The trailer behind it -- `3 rows of 500,004 move (0.0006%) -- 2 respell, 1 merge toward the
submitted form. Previously accepted, so a real migration.` -- reaches no reader. The release
checklist asks whether each entry's inputs were previously **rejected** (free) or previously
**accepted** (a migration), and no entry carries the answer, so a reviewer either ticks the
box without checking or opens every linked PR (#1556).

## Why this is not a git-cliff template

Three routes through the template are measured dead in #1556, and none of them is what this
script uses:

1. **`commit.footers` is mis-parsed.** git_conventional treats any `word:` line as a footer
   token, so a prose body shreds it -- for #1523 it produced seven "footers", the first with
   a value containing the entire rest of the body.
2. **`commit.body` does not contain the trailer.** git_conventional splits body from footers,
   so splitting `commit.body` on `\\nRepresentation-Change:` finds **0** occurrences for both
   real v0.13.1 disclosures.
3. **Splitting on the unanchored literal finds only prose mentions**, which is worse than
   nothing: it produces plausible-looking output for some commits and quoted prose for
   others.

The fourth route is not through the template at all. **Read the trailer from the raw commit
message with git, and attach it to the rendered bullet afterwards.** That is exactly what
`scripts/check_changelog_grouping.py` already does in CI to attribute trailers to commits,
using the same column-anchored `trailer_value` this script imports rather than re-implements
-- three copies of one rule drift silently, which is the lesson `CONTRIBUTING.md` records
about the decline vocabulary.

The bullet is matched to its commit by the **PR number** GitHub puts in the squash subject,
which is the one token both sides are guaranteed to share.

## Correcting a disclosure that is already merged

A trailer lives in a commit message, so once it is on `main` it cannot be corrected without
rewriting released history. The changelog copy, on the other hand, is **rendered by this
script on every run** -- which makes the script, not `CHANGELOG.md`, the place a correction
has to live. `EDITORIAL_CORRECTIONS` below is that register: an entry is appended under its
bullet after the trailer's own words, inside the same blockquote, opening with a label that
says whose words it is. Hand-editing `CHANGELOG.md` instead fails twice over -- release-plz
discards it when it regenerates the pending section, and `--check` reports it as drift
against the trailer.

The first of those two reaches only the **pending** section. A released section is never
regenerated, so registering a correction against a block that is already on disk is
*reported* as drift rather than repaired -- `inject` judges an existing block, it does not
rewrite one. Remove that block from `CHANGELOG.md` first and the re-run attaches the
correction with it.

## Usage

    python scripts/inject_representation_disclosure.py --check      # CI: is CHANGELOG.md current?
    python scripts/inject_representation_disclosure.py              # rewrite CHANGELOG.md in place
    python scripts/inject_representation_disclosure.py --stdout     # print, change nothing

Idempotent: a bullet that already carries its disclosure is left alone, so re-running after
release-plz regenerates the pending section is safe and is the intended workflow.

Exits 0 when nothing needed changing, 1 under `--check` when something did, and 2 on a
bullet whose commit or trailer cannot be found -- a silently-skipped bullet would reintroduce
the exact failure this script exists to remove.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from check_changelog_grouping import (  # noqa: E402
    REPRESENTATION_GROUP,
    strip_ordering_prefix,
    trailer_value,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

#: The marker opening an injected disclosure. Both the idempotence check and `--check` key
#: on it, so it must be something no hand-written section would begin a line with.
DISCLOSURE_PREFIX = "  > "

#: A changelog bullet: `- ...` at column 0, possibly wrapped over later lines.
BULLET_RE = re.compile(r"^- (?P<text>.*\S)\s*$")

#: A section heading. `##` is the version, `###` the group.
HEADING_RE = re.compile(r"^(?P<hashes>#{2,})\s+(?P<title>.*\S)\s*$")

#: The PR number GitHub writes into a squash subject, in the bullet's link or bare.
PR_RE = re.compile(r"\(#(?P<number>\d+)\)|\[#(?P<linked>\d+)\]")

#: Editorial corrections to a merged disclosure, keyed by PR number, one string per
#: paragraph.
#:
#: **Why this is a register and not a `CHANGELOG.md` edit.** The trailer is immutable once
#: its commit is on `main`, and correcting it would mean rewriting released history. Editing
#: the changelog copy instead fails twice: release-plz regenerates the pending section on
#: every push to `main`, and `--check` reports the edited block as drift against the trailer
#: it no longer matches. What *is* durable is this script -- it re-renders the block on every
#: run, so a correction registered here survives every regeneration and is what `--check`
#: then demands.
#:
#: **The trailer's own words stay verbatim.** A correction is rendered as a further
#: paragraph in the same blockquote and must open by naming itself as editorial and citing
#: the issue that raised it, because whose words these are is the distinction the release
#: checklist turns on -- not which block they sit in.
#:
#: **Add an entry only to correct a statement of fact**: something the trailer got wrong, or
#: a cost it did not disclose. Not to summarise it, re-word it, or improve its prose. Append
#: a paragraph rather than rewriting one, so two corrections to one disclosure conflict
#: visibly instead of overwriting each other.
EDITORIAL_CORRECTIONS: dict[str, list[str]] = {
    "1835": [
        "**Editorial correction (#1886) — not the trailer's words.** The trailer names\n"
        "`pre-v0.15.0` as the output `FERRO_PARTITION=live` reproduces. The release carrying\n"
        "this change is the one heading this section, so read that as *output from before\n"
        "this change* — every release up to and including v0.13.1. `release-plz.toml` sets\n"
        "`features_always_increment_minor`, so which release a `feat:` lands in is not\n"
        "knowable when its trailer is written; naming one is now against house style, and\n"
        "`CONTRIBUTING.md` asks for the version-free phrasing instead.",
        "**Editorial correction (#1886) — not the trailer's words.** The trailer discloses\n"
        "confluence only. The flip also moves **2 rows of 932** on the HGVS spec corpus,\n"
        "measured per arm at `45820926` against a prepared GRCh38 reference, 3' direction,\n"
        "with the two rows #1846 cannot terminate on excluded: `live` 646, `shadow` 646,\n"
        "`canonical` 642, **`canonical-coalesced` — the new default — 644**. Two rows stop\n"
        "matching the spec's stated form and none starts, and both are named here because\n"
        "the bare count reads as a regression and neither row is one. `DNA/delins.md`\n"
        "publishes one variant as two corpus rows — the spanning\n"
        "`c.850_901delinsTTCCTCGATGCCTG` that `:47` recommends, and the split `:46` calls\n"
        'its "alternative description" — and both echoed themselves before; the new\n'
        "default converges the split onto the spanning form, which is what any change that\n"
        "converges a published pair must cost. The other row,\n"
        "`LRG_199t1:c.992_1004delinsAC`, is harvested from `consultation/SVD-WG010.md` —\n"
        "the **rejected** proposal, whose answer `tests/it/spec_worked_example_rules.rs`\n"
        "already pins as one ferro must not produce — and the new default no longer echoes\n"
        "it either. Read this as neither more nor less conformant: the measurement supports\n"
        "neither claim. 3' only, because #1879 removes the 5' direction from the public\n"
        "surface in this same release.",
    ],
}


def squash_commits(repo: Path) -> dict[str, str]:
    """Return `{PR number: full commit message}` for every commit with a PR number.

    Framed on NUL for the same reason `check_changelog_grouping.commit_bodies` is: a commit
    message may contain any other byte, and mis-framing here would attribute a trailer to
    the wrong entry in a script whose whole job is attribution.

    A PR number seen twice keeps the **newest** commit, which is the one a reader following
    the link lands on.
    """
    result = subprocess.run(
        ["git", "log", "-z", "--format=%H%x00%B"],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    )
    fields = result.stdout.split("\x00")
    if fields and fields[-1] == "":
        fields.pop()
    if len(fields) % 2 != 0:
        raise RuntimeError(
            f"git log -z returned {len(fields)} NUL-separated fields, which is not an "
            "id/message pairing; refusing to guess the framing"
        )

    by_pr: dict[str, str] = {}
    for index in range(0, len(fields), 2):
        message = fields[index + 1]
        subject = message.splitlines()[0] if message else ""
        match = re.search(r"\(#(\d+)\)\s*$", subject)
        if match and match.group(1) not in by_pr:
            by_pr[match.group(1)] = message
    return by_pr


#: A line that ends the disclosure because it opens a block the trailer cannot own: another
#: column-0 trailer token, a Markdown heading, or an HTML comment.
#:
#: **The HTML-comment arm is not defensive, it is load-bearing.** `CONTRIBUTING.md` requires
#: the trailer to be the last thing in the PR description, and in the squash bodies on this
#: repo it is *not*: GitHub appends CodeRabbit's auto-generated summary after it, opening
#: with `<!-- This is an auto-generated comment: release notes by coderabbit.ai -->`. Taken
#: literally to end-of-message, #1484's disclosure came out with three paragraphs of
#: "Summary by CodeRabbit" quoted under it as though the author had written them.
DISCLOSURE_END_RE = re.compile(r"^([A-Za-z][A-Za-z0-9-]*:[ \t]|#{1,6}\s|<!--)")


def disclosure_value(message: str) -> str | None:
    """Return the whole `Representation-Change:` value, continuations included.

    **`trailer_value` alone is not enough here, and using it alone truncated every
    multi-line disclosure to its first clause** -- `371 rows move of 5,629,002 swept; all
    371 were` cut mid-sentence. It is still what decides *whether* there is a trailer and
    where it starts, so the column-0 rule that `scripts/check_representation_change.py` and
    `scripts/check_changelog_grouping.py` share stays in one place; this only extends the
    value past its first line.

    **The continuation rule has to be "to the end of the message", because this repo's
    trailers are not well-formed** -- `git interpret-trailers --parse` reports *no* trailers
    at all for #1537, #1535 or #1547, since their continuation lines are unindented (#1556).
    So indentation cannot delimit the value.

    **Do not rest this on "the trailer is last".** It is not, and `CONTRIBUTING.md` says so
    as of #1556: GitHub appends CodeRabbit's auto-generated summary after the description,
    so end-of-message would quote three paragraphs of "Summary by CodeRabbit" as the
    author's disclosure -- which is exactly what happened to #1484 before the terminator
    existed. Nothing enforces position either; the trailer checker is column-0-anchored and
    `MULTILINE`, and polices *count*, not placement. What actually bounds the value is the
    terminator: the next column-0 `Token: ` line, Markdown heading, or HTML comment.
    """
    first = trailer_value(message)
    if first is None:
        return None
    lines = message.splitlines()
    start = next(
        index
        for index, line in enumerate(lines)
        if line.lower().startswith("representation-change:")
    )
    collected = [first]
    for line in lines[start + 1 :]:
        if DISCLOSURE_END_RE.match(line):
            break
        collected.append(line.rstrip())
    while collected and not collected[-1].strip():
        collected.pop()
    return "\n".join(collected)


def bullet_pr_number(text: str) -> str | None:
    """Return the PR number a bullet cites, or None."""
    numbers = [m.group("number") or m.group("linked") for m in PR_RE.finditer(text)]
    return numbers[-1] if numbers else None


def _is_disclosure_line(line: str) -> bool:
    """Is this line part of an attached disclosure block?

    Both forms count, and missing the second one is a real defect rather than a nicety: a
    disclosure containing a blank line renders that line as a bare `  >`, because
    `wrap_disclosure` rstrips and `DISCLOSURE_PREFIX` carries a trailing space. Matching
    only the prefix therefore stops collecting at the first paragraph break, so a
    multi-paragraph disclosure the script itself wrote reads back as a one-line block and
    is reported as drift against its own trailer -- exit 2, red CI, on a file nothing had
    touched. Measured on a two-paragraph value before this guard existed.
    """
    return line.startswith(DISCLOSURE_PREFIX) or line == DISCLOSURE_PREFIX.rstrip()


def _flatten(block: list[str]) -> str:
    """Render a `  > ` block as one line, for a side-by-side drift message.

    The two sides differ in *wrapping* as often as in content, so a raw diff of the lists
    reports a difference the reader cannot see. Stripping the quote prefix and joining makes
    the actual disagreement visible in the first 160 characters.
    """
    return " ".join(line.strip("> ").strip() for line in block)


def wrap_disclosure(value: str) -> list[str]:
    """Render `value` as blockquote continuation lines under a bullet.

    A blockquote rather than an indented paragraph because the value is frequently several
    sentences and sometimes several lines: quoting it keeps it visibly the *trailer's own
    words* rather than editorial summary, which is the distinction the checklist turns on.
    """
    return [f"{DISCLOSURE_PREFIX}{line}".rstrip() for line in value.splitlines()]


def wrap_correction(paragraphs: list[str]) -> list[str]:
    """Render `EDITORIAL_CORRECTIONS` paragraphs as further lines of the same blockquote.

    Inside the trailer's own blockquote rather than beside it, deliberately. A correction
    that a reader can meet without the statement it corrects is worse than none, and every
    other placement can be separated from it: a paragraph after the block is dropped the next
    time release-plz regenerates the section, and one before it re-orders on the next
    trailer.

    What keeps the trailer's words distinguishable is therefore the **label**, which each
    paragraph carries, not the block boundary -- so each paragraph is separated by a bare
    `  >` and must open by naming itself. That is also why this returns the separator: a
    correction lazily continuing the trailer's last line would read as the author's own
    sentence, which is the one thing this must never do.
    """
    lines: list[str] = []
    for paragraph in paragraphs:
        lines.append(DISCLOSURE_PREFIX.rstrip())
        lines.extend(f"{DISCLOSURE_PREFIX}{line}".rstrip() for line in paragraph.splitlines())
    return lines


def inject(changelog: str, by_pr: dict[str, str]) -> tuple[str, list[str], list[str]]:
    """Return `(rewritten changelog, changed bullets, problems)`.

    Only bullets inside a **Representation changes** group are touched. A group ends at the
    next heading of any depth, so a hand-written released section -- which has prose and no
    bullets -- passes through untouched.
    """
    lines = changelog.splitlines()
    out: list[str] = []
    changed: list[str] = []
    problems: list[str] = []
    in_section = False

    index = 0
    while index < len(lines):
        line = lines[index]
        heading = HEADING_RE.match(line)
        if heading:
            in_section = strip_ordering_prefix(heading.group("title")) == REPRESENTATION_GROUP
            out.append(line)
            index += 1
            continue

        bullet = BULLET_RE.match(line) if in_section else None
        if bullet is None:
            out.append(line)
            index += 1
            continue

        out.append(line)
        index += 1
        # Consume any continuation the bullet already has, so a re-run does not duplicate a
        # disclosure. What is consumed is KEPT for comparison rather than trusted: the
        # bullet is supposed to carry the trailer's own words, so an existing block that no
        # longer matches its trailer is drift, and drift is what this script exists to
        # catch. An earlier revision `continue`d here the moment it saw a `  > ` line, which
        # meant a stale, truncated or hand-rewritten block passed `--check` green forever —
        # the script could only ever be right on the run that first wrote it.
        existing: list[str] = []
        while index < len(lines) and _is_disclosure_line(lines[index]):
            existing.append(lines[index])
            out.append(lines[index])
            index += 1

        number = bullet_pr_number(bullet.group("text"))
        if number is None:
            problems.append(
                f"bullet cites no PR number, so its trailer cannot be found: "
                f"{bullet.group('text')[:70]!r}"
            )
            continue

        message = by_pr.get(number)
        if message is None:
            problems.append(f"#{number} is in the changelog but no commit in this history cites it")
            continue
        value = disclosure_value(message)
        if value is None:
            problems.append(
                f"#{number} is filed under {REPRESENTATION_GROUP} with no "
                "`Representation-Change:` trailer in its commit message"
            )
            continue

        wanted = wrap_disclosure(value) + wrap_correction(EDITORIAL_CORRECTIONS.get(number, []))
        if existing:
            # Already carried a block: judge it rather than trust it. The trailer in the
            # commit message is the single source of truth for its own words — editing the
            # changelog copy instead makes the two disagree with nothing to notice, which is
            # the failure this whole script exists to remove. Fix by correcting the trailer
            # and re-running, not by editing `CHANGELOG.md`; and where the trailer is already
            # merged and so cannot be corrected, by registering the correction in
            # `EDITORIAL_CORRECTIONS`, which is rendered here and so cannot be regenerated
            # away.
            if existing != wanted:
                problems.append(
                    f"#{number}'s disclosure in the changelog no longer matches its "
                    "`Representation-Change:` trailer (plus any registered editorial "
                    "correction). Those two are the source of truth, and a re-run judges an "
                    "existing block rather than rewriting it: revert the changelog edit, or "
                    "delete the block and let release-plz or a re-run of "
                    "`python scripts/inject_representation_disclosure.py` re-attach it.\n"
                    f"  changelog says: {_flatten(existing)[:160]!r}\n"
                    f"  should say:     {_flatten(wanted)[:160]!r}"
                )
            continue

        out.extend(wanted)
        changed.append(f"#{number}")

    return "\n".join(out) + ("\n" if changelog.endswith("\n") else ""), changed, problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=REPO_ROOT, help="Repository root.")
    parser.add_argument(
        "--changelog",
        type=Path,
        help="Changelog to rewrite. Defaults to `<repo>/CHANGELOG.md`.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Report what would change and exit 1 if anything would, writing nothing.",
    )
    parser.add_argument(
        "--stdout",
        action="store_true",
        help=(
            "Write the result to stdout instead of the file. Mutually exclusive with "
            "`--check`, which reports drift through its exit code."
        ),
    )
    args = parser.parse_args(argv)

    # `--stdout` used to be tested first, so `--stdout --check` printed the rewritten file
    # and returned 0 on a stale changelog — a silent false pass in the one mode whose entire
    # job is to fail. They answer different questions and cannot be combined coherently.
    if args.stdout and args.check:
        parser.error("--stdout and --check are mutually exclusive")

    path = args.changelog or (args.repo / "CHANGELOG.md")
    original = path.read_text(encoding="utf-8")
    try:
        rewritten, changed, problems = inject(original, squash_commits(args.repo))
    except RuntimeError as error:
        print(str(error), file=sys.stderr)
        return 2
    except subprocess.CalledProcessError as error:
        # `squash_commits` runs `git log` with `check=True`. Unguarded, a git failure exits
        # 1 with a traceback — colliding with the documented meaning of exit 1 ("--check
        # found drift"), so a broken environment reads as a stale changelog.
        print(f"git failed while reading this history: {error}", file=sys.stderr)
        return 2

    for problem in problems:
        print(problem, file=sys.stderr)

    # Refuse before writing anything. Previously the write happened first and the exit code
    # came after, so a run with one resolvable and one unresolvable bullet left a PARTIALLY
    # injected changelog on disk *and* exited 2 — the worst pairing available, because the
    # file looks updated and the failure looks like a transient. A partial artifact is not
    # a better outcome than no artifact.
    if problems:
        return 2

    if args.stdout:
        sys.stdout.write(rewritten)
        return 0

    if args.check:
        if changed:
            print(
                f"{len(changed)} changelog entr{'y' if len(changed) == 1 else 'ies'} would gain "
                f"a disclosure: {', '.join(changed)}. Run "
                "`python scripts/inject_representation_disclosure.py`.",
                file=sys.stderr,
            )
        # `problems` is already handled above by an early `return 2`, so drift is the only
        # thing left to report here.
        return 1 if changed else 0

    if changed:
        path.write_text(rewritten, encoding="utf-8")
        print(f"attached {len(changed)} disclosure(s): {', '.join(changed)}")
    else:
        print("every Representation changes entry already carries its disclosure")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
