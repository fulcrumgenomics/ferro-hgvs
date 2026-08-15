"""
Tests for `scripts/check_representation_change.py` (#1522).

The check exists because an *absent* trailer and a *considered-and-declined* trailer
looked identical for a whole release cycle. So the assertions that matter here are the
ones separating those two states, and the one pinning that the watched set still matches
`release-plz.toml`'s stated scope — a directory silently dropping out of that tuple would
make this check pass on exactly the changes it exists to catch.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

_MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts" / "check_representation_change.py"
_SPEC = importlib.util.spec_from_file_location("check_representation_change", _MODULE_PATH)
assert _SPEC is not None and _SPEC.loader is not None
check_representation_change = importlib.util.module_from_spec(_SPEC)
sys.modules["check_representation_change"] = check_representation_change
_SPEC.loader.exec_module(check_representation_change)

check = check_representation_change.check
find_declaration = check_representation_change.find_declaration
find_declarations = check_representation_change.find_declarations
find_near_misses = check_representation_change.find_near_misses
watched_files = check_representation_change.watched_files
WATCHED_PREFIXES = check_representation_change.WATCHED_PREFIXES
fenced_line_numbers = check_representation_change.fenced_line_numbers
find_trailers_as_git_cliff_would = check_representation_change.find_trailers_as_git_cliff_would
_REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# The two states the whole check exists to separate
# ---------------------------------------------------------------------------


def test_watched_change_without_a_trailer_fails() -> None:
    """Silence is what this rejects — the v0.13.0 cycle's actual failure mode."""
    ok, message = check(["src/normalize/merge.rs"], "Fixes a thing.\n\nCloses #1.")
    assert not ok
    assert "src/normalize/merge.rs" in message
    assert "Representation-Change: none" in message, "the message must name the way to decline"


def test_watched_change_declaring_none_passes() -> None:
    """Declining is a declaration. #1521 is the real example: 46 changed lines in
    `src/normalize/merge.rs`, every one a comment."""
    ok, _ = check(
        ["src/normalize/merge.rs"],
        "Comments only; zero non-comment lines change.\n\nRepresentation-Change: none",
    )
    assert ok


def test_watched_change_declaring_a_move_passes() -> None:
    ok, message = check(
        ["src/spdi/mod.rs"],
        "Representation-Change: 577 rows move, 360 merge / 205 split / 12 respell",
    )
    assert ok
    assert "577 rows move" in message


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------


def test_unwatched_change_needs_no_declaration() -> None:
    ok, _ = check(["README.md", "tests/it/foo.rs", "src/cli/mod.rs"], "no trailer here")
    assert ok


def test_a_watched_file_among_unwatched_ones_still_requires_a_declaration() -> None:
    """The trigger is *any* watched file, not a majority of them."""
    ok, _ = check(["README.md", "docs/x.md", "src/project/projector.rs"], "no trailer")
    assert not ok


#: One representative path per watched directory. Kept as data because the two tests below
#: read it from opposite sides — every entry must be gated, and the tuple must contain
#: nothing this list does not name — which together pin the watched set exactly, through
#: `check()` rather than by restating the constant.
_A_FILE_IN_EVERY_WATCHED_DIRECTORY = (
    "src/normalize/merge.rs",
    "src/hgvs/variant.rs",
    "src/spdi/mod.rs",
    "src/project/projector.rs",
    # Added by #1853's backtest. `src/reference/` decides which bases a description
    # resolves against and `src/error_handling/` decides the accept/reject boundary; both
    # have moved consumer-visible output while the gate stayed silent -- see
    # `test_the_widened_prefixes_are_the_ones_the_backtest_measured`.
    "src/reference/multi_fasta.rs",
    "src/error_handling/mod.rs",
)

#: Directories that changed output-adjacent code and are deliberately NOT watched. Each
#: carries its reason, because "not watched" is a decision that gets re-proposed.
_A_FILE_IN_EVERY_DELIBERATELY_UNWATCHED_DIRECTORY = (
    # Measurement and adjudication code. It cannot move ferro's output, so requiring a
    # declaration would demand a statement about something the directory cannot do -- both
    # of its real disclosures (#1623 merged, #1839 open) declare "0 rows move in this PR".
    "src/conformance/spec_corpus.rs",
    # Declined on evidence in #1853's backtest: zero marginal catch in both populations.
    # Its one real disclosure (#1737) also touches `src/reference/`, so the widening above
    # already gates it; both its merged gate-era commits (#1786, #1770) declined.
    "src/data/cdot.rs",
    # Zero real disclosures in either population.
    "src/convert/mapper.rs",
)


@pytest.mark.parametrize("path", _A_FILE_IN_EVERY_WATCHED_DIRECTORY)
def test_every_watched_prefix_triggers(path: str) -> None:
    assert watched_files([path]) == [path]
    ok, _ = check([path], "no trailer")
    assert not ok, f"{path} must require a declaration"


@pytest.mark.parametrize("path", _A_FILE_IN_EVERY_DELIBERATELY_UNWATCHED_DIRECTORY)
def test_a_deliberately_unwatched_directory_needs_no_declaration(path: str) -> None:
    """The exclusions are decisions, so pin them the same way the inclusions are pinned.

    Without this, widening the tuple by one more directory is a silent change: nothing
    fails, and the next reader cannot tell an excluded directory from an overlooked one.
    """
    assert watched_files([path]) == []
    ok, _ = check([path], "no trailer")
    assert ok, f"{path} must not require a declaration"


def test_the_widened_prefixes_are_the_ones_the_backtest_measured() -> None:
    """`src/reference/` and `src/error_handling/` are gated; `src/data/` is not (#1853).

    The three were proposed together and the backtest supports two of them. Both additions
    are carried by the v0.13.0 cycle, the one place movement was measured without relying
    on anyone declaring it: its 326,404 newly-normalizing rows attribute to #1490
    (`src/reference/` and nothing else) and #1501 (`src/error_handling/` and nothing else),
    and the cycle shipped with zero trailers. `src/data/` was touched by **0** commits in
    that cycle and has zero marginal catch in either population measured since.

    Asserted through `check()` on synthetic file lists, so it survives any rewrite of how
    the tuple is spelled -- and paired with the coverage assertion below, which is what
    makes it a pin rather than three examples: a seventh directory cannot be added in
    silence.
    """
    for path in ("src/reference/multi_fasta.rs", "src/error_handling/mod.rs"):
        ok, message = check([path], "Fixes a thing.\n\nCloses #1.")
        assert not ok, f"{path} must require a declaration"
        assert path in message, "the failure must name the file that demanded the trailer"
        ok, _ = check([path], "Representation-Change: none")
        assert ok, "declining must still pass on a newly watched directory"

    ok, _ = check(["src/data/cdot.rs"], "no trailer")
    assert ok, (
        "src/data/ was proposed in #1853 and declined on measurement; adding it needs a "
        "real, non-declining disclosure that touches src/data/ and NOT src/reference/, "
        "which did not exist in either measured population"
    )

    # Coverage, not a count. A count is defeated by two compensating edits -- add a
    # seventh prefix, and add a seventh representative under a directory that is already
    # covered -- which leaves the lengths equal, every representative gated by
    # `test_every_watched_prefix_triggers`, and the new directory with no representative
    # at all. That is the same weakness `test_watched_prefixes_match_the_release_config`
    # rejects in its own docstring, so do not reintroduce it here: the map from
    # representatives to prefixes is a bijection only if no two representatives share a
    # prefix, and nothing asserts that.
    unrepresented = [
        prefix
        for prefix in WATCHED_PREFIXES
        if not any(path.startswith(prefix) for path in _A_FILE_IN_EVERY_WATCHED_DIRECTORY)
    ]
    assert not unrepresented, (
        f"{sorted(unrepresented)} is watched but has no representative path in "
        "_A_FILE_IN_EVERY_WATCHED_DIRECTORY; every entry needs one so the gate is pinned "
        "by behaviour rather than by restating the constant"
    )


def test_watched_prefixes_match_the_release_config() -> None:
    """`release-plz.toml`'s reviewer checklist names the directories this must watch.

    A real set comparison, in **both** directions. This test's docstring used to claim it
    was one while the assertion only ran checker -> config, which left the dangerous
    direction open: a prefix silently dropped from `WATCHED_PREFIXES` still satisfied
    "every watched directory is named in the config", and a dropped prefix is exactly the
    drift that makes this check pass on the changes it exists to catch. A count would not
    close it either, since two compensating edits keep the count while swapping a
    directory out.

    The config side is read as every backticked ``src/<dir>/`` token in the file, so the
    comparison keeps working when the prose is re-wrapped -- which it was here, the list
    having grown past one line.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    documented = {f"{directory}/" for directory in re.findall(r"`(src/[a-z_]+)/`", config)}
    assert documented == set(WATCHED_PREFIXES), (
        f"release-plz.toml names {sorted(documented)} but the check watches "
        f"{sorted(WATCHED_PREFIXES)}; the two must describe the same scope. A directory in "
        "the config and not the tuple is a gate that does not exist; one in the tuple and "
        "not the config is a scope the release reviewer is never told to look for."
    )


@pytest.mark.parametrize("document", ["CONTRIBUTING.md", "CLAUDE.md", ".github/workflows/ci.yml"])
def test_the_prose_restatements_name_every_watched_directory(document: str) -> None:
    """The list is restated in prose in three more places; none of them was pinned.

    `release-plz.toml` is pinned above, but the same six directories are also written out
    in `CONTRIBUTING.md` (contributor guidance), `CLAUDE.md` (agent guidance) and the
    `representation-change` job's comment in `ci.yml` -- which says of itself that it
    "must be kept in step with" the constant, an obligation nothing enforced. One rule
    written in several places and then drifting apart is this repository's named recurring
    failure mode, and the sibling `test_contributing_documents_the_same_decline_vocabulary`
    already pins the decline vocabulary for exactly that reason.

    **Containment, not set equality**, and the asymmetry is deliberate. All three documents
    legitimately name directories that are *not* watched -- `src/conformance/` and
    `src/data/`, whose exclusion is itself a decision worth stating -- so equality would
    fail on correct prose. The direction that matters is the other one: a directory added
    to `WATCHED_PREFIXES` and not to the docs is a required check that fails a contributor
    who was never told the scope had grown.
    """
    text = (Path(__file__).resolve().parents[2] / document).read_text(encoding="utf-8")
    undocumented = [prefix for prefix in WATCHED_PREFIXES if f"`{prefix.rstrip('/')}/`" not in text]
    assert not undocumented, (
        f"{document} does not name {sorted(undocumented)}, which "
        f"scripts/check_representation_change.py watches; the check watches "
        f"{sorted(WATCHED_PREFIXES)}. A watched directory missing from the prose is a "
        "required check nobody was told about."
    )


# ---------------------------------------------------------------------------
# Trailer parsing
# ---------------------------------------------------------------------------


def test_trailer_is_found_case_insensitively() -> None:
    assert find_declaration("representation-change: none") == "none"
    assert find_declaration("REPRESENTATION-CHANGE: none") == "none"


def test_trailer_is_found_at_the_end_of_a_body() -> None:
    body = "Some description.\n\nCloses #1522.\n\nRepresentation-Change: none\n"
    assert find_declaration(body) == "none"


def test_absent_trailer_reads_as_none_not_as_a_declaration() -> None:
    assert find_declaration("Representation change: none") is None, (
        "the hyphenated trailer is what git-cliff matches; prose must not satisfy it"
    )


def test_an_empty_trailer_value_is_not_a_declaration() -> None:
    """`Representation-Change:` with nothing after it declares nothing."""
    assert find_declaration("Representation-Change:   \n") is None


@pytest.mark.parametrize("indent", [" ", "  ", "\t", "    "])
def test_an_indented_trailer_is_not_a_trailer(indent: str) -> None:
    """The token must sit at column 0, because that is where git and git-cliff require it.

    `CONTRIBUTING.md` documents indentation as the way to continue a trailer's *value* onto
    another line, so a checker that read an indented line as a new trailer contradicted the
    convention it enforces. A sole indented trailer now fails as absent — which is the
    honest answer, since git-cliff would not group it either and the disclosure would be
    lost rather than merely unread.
    """
    assert find_declaration(f"{indent}Representation-Change: none") is None
    ok, message = check(["src/normalize/merge.rs"], f"{indent}Representation-Change: none")
    assert not ok, "an indented trailer must not satisfy the check"
    assert "Representation-Change: none" in message, "the message must show the column-0 form"


@pytest.mark.parametrize("value", ["none", "None", "NONE", "no", "n/a", "na"])
def test_decline_spellings_all_pass(value: str) -> None:
    ok, _ = check(["src/normalize/merge.rs"], f"Representation-Change: {value}")
    assert ok


# ---------------------------------------------------------------------------
# Exactly one trailer, because two do not mean anything consistent
# ---------------------------------------------------------------------------

#: A real disclosure followed by a *quoted* declining example at column 0. The shape is
#: invited by `CONTRIBUTING.md`, which documents the declining form, and by a corrected
#: trailer added without deleting the superseded one.
_TWO_TRAILERS = (
    "Representation-Change: 577 rows move, 360 merge.\n"
    "\n"
    "For contrast, a declining trailer looks like:\n"
    "Representation-Change: none\n"
)


def test_two_trailers_are_refused() -> None:
    """Neither half of the machinery can resolve two trailers the same way, so refuse.

    Measured against git-cliff 2.13.1 with the real config: this message groups under
    **Other**, because git-cliff matches its footer rule against *every* footer and takes
    the first rule that matches any of them, so a decline wins wherever it sits. This
    checker reads the *first* trailer and calls the same message a real disclosure. The
    commit then passes CI as a disclosure and is filed as a decline, and the disclosure
    disappears from the changelog — the under-reporting direction of #1522, one commit at
    a time.

    It is not an anchoring bug and #1573's `m` does not reach it: the shape reproduces
    identically under `(?i)`, `(?im)` and `(?im)\\A`, because the second line is parsed as
    its own footer and the rule matches at *that* footer's start.

    Refusing rather than picking a winner keeps the rule that a decline must be *stated*,
    never inferred from a message that also says the opposite.
    """
    ok, message = check(["src/normalize/merge.rs"], _TWO_TRAILERS)
    assert not ok, "two trailers must not pass; the changelog and this check disagree on them"
    assert "2 `Representation-Change:` trailers found" in message
    assert "577 rows move, 360 merge." in message, "the message must show what it found"
    assert "none" in message


def test_two_trailers_are_refused_even_with_no_watched_file() -> None:
    """The harm is in changelog grouping, which applies to every commit.

    So this is refused ahead of the watched-file test — a docs-only commit reaches the
    changelog exactly like a `src/normalize/` one.
    """
    ok, message = check(["README.md"], _TWO_TRAILERS)
    assert not ok
    assert "trailers found" in message


def test_two_declining_trailers_are_still_refused() -> None:
    """Even when both agree, one of them is unmaintained. Ambiguity is the defect."""
    ok, _ = check(["src/hgvs/variant.rs"], "Representation-Change: none\nRepresentation-Change: no")
    assert not ok


def test_an_indented_second_trailer_is_a_continuation_not_a_second_trailer() -> None:
    """The escape hatch the refusal message names must actually work.

    Indenting makes the line a continuation of the value above it, which is also how
    git-cliff sees it — measured, the indented form groups under **Representation
    changes** where the column-0 form groups under **Other**. So the fix the message
    recommends is the one that makes both halves agree.
    """
    body = (
        "Representation-Change: 577 rows move, 360 merge.\n"
        "\n"
        "For contrast, a declining trailer looks like:\n"
        "    Representation-Change: none\n"
    )
    assert find_declarations(body) == ["577 rows move, 360 merge."]
    ok, _ = check(["src/normalize/merge.rs"], body)
    assert ok


def test_one_trailer_still_passes() -> None:
    """The ordinary case, pinned so the refusal above cannot widen onto it."""
    for body in (
        "Representation-Change: none",
        "Representation-Change: none. Tests only.",
        "Some prose.\n\nCloses #1.\n\nRepresentation-Change: 577 rows move, 360 merge.\n",
    ):
        ok, _ = check(["src/normalize/merge.rs"], body)
        assert ok, f"{body!r} carries exactly one trailer and must pass"


def test_find_declarations_returns_every_value_in_order() -> None:
    assert find_declarations("no trailer here") == []
    assert find_declarations("Representation-Change: none") == ["none"]
    assert find_declarations(_TWO_TRAILERS) == ["577 rows move, 360 merge.", "none"]


# ---------------------------------------------------------------------------
# A decline may give its reason (#1555)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        "none. Tests only; no file under a watched directory is touched.",
        "none: comments only",
        "none; the diff is one fixture JSON and one new test module",
        "none — 0 of 950 real cis-allele rows move their normalized string",
        "None. Ruling records and their status pins only.",
        "no. Nothing here can reach a normalizer.",
        "n/a: docs",
        "na. generated artifact only",
        "none.\n  A reason spanning two lines.",
    ],
)
def test_a_decline_may_give_its_reason(value: str) -> None:
    """The verdict is the first word; what follows a terminator explains it.

    Requiring a bare `none` punished exactly the contributors who explained themselves:
    8 commits in the v0.13.1 cycle declined with a reason and all 8 were then listed as
    representation changes.
    """
    assert check_representation_change.declines(value), f"{value!r} declines a move"


@pytest.mark.parametrize(
    "value",
    [
        "577 rows move, 360 merge / 205 split / 12 respell",
        "3 rows of 500,004 move (0.0006%) — 2 respell, 1 merge",
        "0 rows move over 5,761,302 real expressions",
        "no rows move",
        "nothing moves under src/normalize/",
        "not measured yet",
        "none, except two rows that merge",
    ],
)
def test_a_description_of_a_move_is_not_a_decline(value: str) -> None:
    """The failure that matters here is the silent one: a real disclosure read as a decline
    disappears from the changelog, and nobody is told.

    `no rows move` and `none, except …` are the near misses. Both begin with a decline word,
    and neither is one — the first has no terminator, and `,` is not a terminator precisely
    because it usually introduces a qualification that changes the verdict.
    """
    assert not check_representation_change.declines(value), f"{value!r} describes a move"


def test_an_empty_trailer_does_not_scavenge_the_next_line() -> None:
    """`\\s` spans newlines, so a permissive class would read `none` off the line below and
    call an empty trailer a declaration."""
    assert find_declaration("Representation-Change:\nnone\n") is None


def test_empty_trailer_on_a_watched_change_fails() -> None:
    ok, _ = check(["src/normalize/merge.rs"], "Representation-Change:   \n")
    assert not ok


def test_both_streams_on_stdin_is_rejected() -> None:
    """Reading both inputs from stdin would silently give the second an empty string,
    which for the declaration reads as "no trailer" — a wrong failure, not an error."""
    with pytest.raises(SystemExit) as excinfo:
        check_representation_change.main(["--changed-files", "-", "--declaration-file", "-"])
    assert excinfo.value.code == 2


# ---------------------------------------------------------------------------
# The changelog consumer must agree with this checker about what "declining" means
# ---------------------------------------------------------------------------


def test_decline_vocabulary_matches_the_changelog_config() -> None:
    """`release-plz.toml`'s exclusion rule and `NONE_VALUES` must accept the same words.

    They are two halves of one decision and they drift silently. If the checker accepts a
    decline the changelog rule does not, that PR passes CI and then renders under
    **Representation changes** as if it moved output — which is #1526, the bug this test
    was written for. If the changelog excludes a word the checker rejects, a contributor
    is told to declare something the changelog then hides.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    # `[a-z]*i[a-z]*` rather than `[a-z]+`: the flag group gained `m` in #1573, so matching
    # `(?i)` literally would read as a rename -- but dropping the `i` requirement entirely
    # would let this find a case-*sensitive* rule while asserting it found the opposite.
    match = re.search(
        r'footer = "\(\?[a-z]*i[a-z]*\)\^Representation-Change:[^"]*?\(([^)]+)\)', config
    )
    assert match is not None, (
        "no case-insensitive Representation-Change exclusion rule found in release-plz.toml; "
        "the decline vocabulary cannot be checked"
    )
    configured = frozenset(match.group(1).split("|"))
    assert configured == check_representation_change.NONE_VALUES, (
        f"release-plz.toml excludes {sorted(configured)} but the checker treats "
        f"{sorted(check_representation_change.NONE_VALUES)} as declines; a word in one and "
        "not the other either leaks a non-change into the changelog or hides a real one"
    )


# ---------------------------------------------------------------------------
# A decline that describes a move contradicts itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        "no. 3 rows move",
        "none. 577 rows move, 360 merge / 205 split",
        "none. 3 rows of 500,004 move",
        "na — but 2 rows respell",
        "none. 12,530 rows move in the synthetic corpus",
        # The past tense of every verb, not only of `move`. `moved` was covered and
        # `merged`/`respelled` were not, so a trailer could contradict itself in the
        # tense a retrospective disclosure is most naturally written in.
        "none. 3 rows merged",
        "no. 12 rows respelled",
        "none. 205 rows moved",
    ],
)
def test_a_decline_that_describes_a_move_fails(value: str) -> None:
    """Reading the verdict from the first word lets a trailer contradict its own reason.

    The failure is silent and in the dangerous direction — filed as `none`, so the
    disclosure never reaches the changelog and nobody is told.
    """
    ok, message = check(["src/normalize/merge.rs"], f"Representation-Change: {value}")
    assert not ok, f"{value!r} declines and then describes a move"
    assert "declines and then describes a move" in message


@pytest.mark.parametrize(
    "value",
    [
        # #1535, verbatim: a decline is entitled to quantify its zero.
        "0 rows move over 5,761,302 real expressions",
        # #1547, verbatim: the count that sits next to `rows` is 950, and it is not the
        # moving set. A looser pattern fires here — measured — and this is an exemplary
        # decline, so failing it would punish the phrasing worth encouraging.
        "none. 0 of 950 real cis-allele rows move their normalized string; 1 of 950 "
        "changes strict verdict from accept to reject, which is the defect being fixed.",
        # #1546, verbatim in shape: the corpus grew; no shipped string moved.
        "none. This change does grow the corpus, 78,028 -> 78,298 rows, which changes the "
        "denominator of every future compare run against it.",
        # #1538, verbatim in shape: a count of adjudicated rows, not of moved ones.
        "none. The ruling ratifies shipped behaviour — v0.12.0 already emits the unsplit "
        "form on 208 of 208 adjudicated rows.",
        # A zero is a zero however it is spelled. The rule excludes zero by requiring a
        # nonzero digit in the count; an anchored `(?!0\\b)` excluded only the bare `0`,
        # so these read as nonzero and failed a decline that moves nothing.
        "none. 0,000 rows move",
        "none. 000 rows merged",
    ],
)
def test_a_numerate_decline_is_not_a_contradiction(value: str) -> None:
    """The guard must not fire on the declines this repository actually writes.

    Every one of these is real, and they are numerate on purpose. Verified against the 15
    declines of the v0.13.1 cycle: this rule fires on none of them.
    """
    assert check_representation_change.contradicted_decline(value) is None, (
        f"{value!r} declines legitimately; failing it would punish a good disclosure"
    )


def test_a_real_disclosure_is_not_a_contradiction() -> None:
    """The guard applies only to declines — a trailer that leads with the move is fine."""
    assert (
        check_representation_change.contradicted_decline(
            "3 rows of 500,004 move (0.0006%) — 2 respell, 1 merge"
        )
        is None
    )


# ---------------------------------------------------------------------------
# #1647: the count is the count, never the denominator
#
# The guard used to read "the number immediately before `rows`" as the moving
# count. In `0 of 950 rows move` that number is the DENOMINATOR, so the form
# both `CONTRIBUTING.md` and `CLAUDE.md` publish as the way to quantify a zero
# — and both state passes — failed the build instead. The documented form was a
# merge blocker.
#
# These are separated from the parametrized declines above because they pin a
# specific past defect rather than the general "declines are numerate" rule.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value",
    [
        # `CONTRIBUTING.md` and `CLAUDE.md`, verbatim. The whole of #1647.
        "none. 0 of 950 rows move",
        "none. 0 of 950 rows move.",
        # A grouped denominator. `(?<!of\\s)` alone does NOT rescue this one: the
        # comma opens a word boundary, so the scan restarts at `298` — nonzero,
        # immediately before `rows`, not preceded by `of `. Measured; it is why
        # the second `(?<![\\d,])` lookbehind exists.
        "none. 0 of 78,298 rows move",
        "none. 0 of 5,761,302 rows respell",
        # Both orders of the same honest zero.
        "none. 0 rows of 78,298 move",
        # A noun phrase between the denominator and `rows`, which the repository
        # uses freely -- #1547's own exemplary decline is this shape.
        "none. 0 of 950 real cis-allele rows move",
        "none. 0 of 78,298 corpus rows move",
    ],
)
def test_a_quantified_zero_in_either_order_passes(value: str) -> None:
    """The documented way to quantify a zero must not be a merge blocker."""
    assert check_representation_change.contradicted_decline(value) is None, (
        f"{value!r} is the documented quantified-zero form and must pass"
    )


@pytest.mark.parametrize(
    "value",
    [
        # Count first, denominator between the count and `rows`. The tripwire has
        # to reach this order too, or a real disclosure vanishes from the
        # changelog behind a decline — the one direction this mechanism must
        # never fail in.
        "no. 3 of 500,004 rows move",
        "none. 1,826 of 78,298 rows move",
        # A noun phrase before `rows`. This is #1651's ACTUAL trailer wording with
        # the verdict flipped to a decline, and the two-lookbehind pattern let it
        # through -- a real 1826-row move hidden behind `none`, which is the one
        # direction this mechanism must never fail in. Measured, then fixed.
        "none. 1826 of 78298 corpus rows move",
        "none. 3 of 500,004 real cis-allele rows move",
        # The order the guard already caught, kept here as the control.
        "no. 3 rows move",
        "none. 2 rows of 500,004 respell",
    ],
)
def test_a_nonzero_count_in_either_order_is_still_a_contradiction(value: str) -> None:
    """Rescuing the zero must not blunt the tripwire for a real move."""
    assert check_representation_change.contradicted_decline(value) is not None, (
        f"{value!r} declines while disclosing a move and must be refused"
    )


@pytest.mark.parametrize(
    "value",
    [
        # Punctuation must stop the intervening-word run, or the leading count of
        # an unrelated clause gets read as the moving count. Without that, `1` here
        # would pair with the `rows move` three words later and this honest decline
        # would fail the build.
        "none. 1 issue filed, 0 rows move",
        "none. 12 clauses cited; 0 rows move",
        "none. 4 tests added, no rows move",
        # No number at all, and the ordinary shape of most declines.
        "none. Tests only; no watched file is touched.",
    ],
)
def test_allowing_words_before_rows_did_not_blunt_the_decline_path(value: str) -> None:
    """The control on the widening: three words, not any distance.

    Reaching a noun phrase before `rows` means a nonzero count elsewhere in the
    sentence could pair with a *different* clause's `rows move`. Punctuation is
    what prevents it, so these pin the boundary rather than the capability.
    """
    assert check_representation_change.contradicted_decline(value) is None, (
        f"{value!r} is an honest decline and must not be read as a disclosure"
    )


def test_the_documented_zero_passes_the_whole_check_not_only_the_predicate() -> None:
    """End-to-end, because #1647 was reported against the CLI's exit code.

    `contradicted_decline` returning `None` is necessary but not sufficient —
    the reported symptom was `check_representation_change.py` exiting 1 on a
    watched change, which is what actually blocks the merge.
    """
    ok, _message = check(
        ["src/normalize/merge.rs"],
        "Representation-Change: none. 0 of 950 rows move.",
    )
    assert ok, "the documented quantified-zero trailer must pass the full check"


def test_the_two_decline_rules_agree_value_by_value() -> None:
    """The vocabulary test above compares word lists; this one compares *verdicts*.

    #1555 is the reason both are needed. The two rules named the same four words and still
    disagreed, because the changelog rule anchored the value (`none$`) while the checker
    lowercased and compared it — so `none. Tests only.` passed CI as a decline and then
    rendered as a representation change. A word-list comparison cannot see that; only
    running both rules over the same values can.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    # Anchored on the trailer name, not on `{ footer = "`: the first such literal in the file
    # is git-cliff's `^changelog: ?ignore` example, quoted in a comment.
    rules = re.findall(r'\{ footer = "([^"]*Representation-Change[^"]*)"', config)
    assert rules, "no Representation-Change exclusion rule found in release-plz.toml"
    # The rule is a TOML basic string, so its backslashes are escaped in the file.
    #
    # Compiled with NO flags of our own: whatever the rule needs it must declare inline,
    # because git-cliff adds none either. Passing `re.MULTILINE` here -- as this test did
    # until #1573 -- makes the harness kinder than production and hides exactly the class
    # of anchoring bug it exists to catch.
    changelog_rule = re.compile(rules[0].replace("\\\\", "\\"))

    values = [
        "none",
        "NONE",
        "none.",
        "none. Tests only; nothing under a watched directory.",
        "none: comments only",
        "none; one fixture JSON",
        "none — 0 of 950 rows move",
        "no. Nothing reaches a normalizer.",
        "n/a: docs",
        "na. generated artifact only",
        "577 rows move, 360 merge / 205 split",
        "0 rows move over 5,761,302 real expressions",
        "no rows move",
        "none, except two rows that merge",
        "nothing moves",
    ]
    for value in values:
        by_changelog = changelog_rule.search(f"Representation-Change: {value}") is not None
        by_checker = check_representation_change.declines(value)
        assert by_changelog == by_checker, (
            f"{value!r}: release-plz.toml calls it "
            f"{'a decline' if by_changelog else 'a real change'} and the checker calls it "
            f"{'a decline' if by_checker else 'a real change'}; the changelog and CI must "
            "agree, or a PR passes the check and is then filed as its opposite"
        )


#: What a bot appends after the trailer. CodeRabbit posts a block like this on essentially
#: every PR, and GitHub builds the squash commit body from the PR description, so it reaches
#: the commit message verbatim and lands *after* the trailer.
_APPENDED_BLOCK = (
    "\n\n<!-- This is an auto-generated comment: release notes by coderabbit.ai -->\n"
    "## Summary by CodeRabbit\n"
    "- Tests: added coverage.\n"
)


@pytest.mark.parametrize(
    ("value", "is_a_decline"),
    [
        ("none", True),
        ("NONE", True),
        ("none. Comment prose only.", True),
        ("none; tests only", True),
        ("n/a: docs", True),
        ("577 rows move, 360 merge / 205 split.", False),
        ("none, except two rows that merge.", False),
    ],
)
def test_a_trailer_keeps_its_verdict_when_text_is_appended_after_it(
    value: str, is_a_decline: bool
) -> None:
    """Trailing text must not invert the verdict — #1573, and #1526 in a new shape.

    git-cliff matches the exclusion rule against the **whole footer value**, which spans
    lines. Anchored with `$` and no `m` flag, the rule only fired when the decline was the
    last thing in the message, so a bare `none` followed by an appended block was filed
    under **Representation changes**. `3ebcdddd` (#1551) is the real instance: it declined,
    CodeRabbit appended release notes, and the v0.14.0 changelog listed it as a change.

    The reason-terminator set from #1555 does not cover this: after a bare decline the next
    character is a newline, not one of `.;:—–`.

    Every other guard in this file inspects the config as a *string* — vocabulary, casing,
    ordering — and all of them passed throughout the bug, because anchoring is a behaviour.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    rules = re.findall(r'\{ footer = "([^"]*Representation-Change[^"]*)"', config)
    assert rules, "no Representation-Change exclusion rule found in release-plz.toml"
    # No flags of our own: git-cliff supplies none, so the rule must declare what it needs.
    changelog_rule = re.compile(rules[0].replace("\\\\", "\\"))

    footer = f"Representation-Change: {value}{_APPENDED_BLOCK}"
    assert (changelog_rule.search(footer) is not None) == is_a_decline, (
        f"{value!r} followed by an appended block is filed as "
        f"{'a real change' if is_a_decline else 'a decline'} by release-plz.toml; text "
        "appended after a trailer must not change what the trailer says"
    )

    # And the checker must reach the same verdict on the same message, since disagreeing is
    # what lets a PR pass CI and then render as its opposite.
    declaration = find_declaration(footer)
    assert declaration is not None, "the checker no longer finds a trailer it used to find"
    assert check_representation_change.declines(declaration) == is_a_decline, (
        f"{value!r}: the checker and release-plz.toml disagree once a block is appended"
    )


def test_the_ordering_prefix_is_stripped_from_rendered_headings() -> None:
    """`<!-- N -->` orders the groups; it must not survive into the changelog text.

    Without the postprocessor every released section carries `### <!-- 0 -->Representation
    changes` in the file. It is an HTML comment, so GitHub renders it invisibly — which is
    how it survived two releases unnoticed.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    match = re.search(
        r'postprocessors = \[\{ pattern = "((?:[^"\\]|\\.)*)", replace = "([^"]*)"', config
    )
    assert match is not None, "no postprocessor found; the ordering prefix would be rendered"
    pattern = re.compile(match.group(1).replace("\\\\", "\\"), re.MULTILINE)
    replacement = match.group(2).replace("${1}", r"\1")

    rendered = "### <!-- 0 -->Representation changes\n### <!-- 7 -->Other\n"
    assert pattern.sub(replacement, rendered) == "### Representation changes\n### Other\n"

    subject = "- *(docs)* explain the <!-- 0 --> marker\n"
    assert pattern.sub(replacement, subject) == subject, (
        "the rule must be anchored to a heading; a marker quoted in a commit subject is text"
    )


def test_both_representation_change_rules_are_case_insensitive() -> None:
    """A case-sensitive *inclusion* rule drops a lowercase `representation-change:` — a real
    disclosure — silently into `Other`, since the checker accepts any casing."""
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    rules = re.findall(r'\{ footer = "([^"]*Representation-Change[^"]*)"', config)
    assert len(rules) == 2, f"expected an exclusion and an inclusion rule, found {rules}"
    for rule in rules:
        # Read the inline flag group rather than matching `(?i)` literally, so that adding
        # another flag -- `m`, in #1573 -- does not read as "the rule became case-sensitive".
        flags = re.match(r"\(\?([a-z]+)\)", rule)
        assert flags is not None and "i" in flags.group(1), (
            f"rule {rule!r} is case-sensitive; the checker's own trailer regex is not, so the "
            "two disagree about whether `REPRESENTATION-CHANGE:` is a declaration"
        )


def test_contributing_documents_the_same_decline_vocabulary() -> None:
    """CONTRIBUTING.md is the third place these words appear, and prose drifts fastest.

    A contributor reads the doc, not the config or this script, so a word documented here
    but not accepted by both is guidance that fails CI, and one accepted but not documented
    is a decline nobody knows they can make.
    """
    doc = (Path(__file__).resolve().parents[2] / "CONTRIBUTING.md").read_text(encoding="utf-8")
    documented = {
        word
        for word in re.findall(r"`([a-z/]+)`(?=[,\s]|$)", doc)
        if word in check_representation_change.NONE_VALUES
    }
    assert documented == check_representation_change.NONE_VALUES, (
        f"CONTRIBUTING.md documents {sorted(documented)} as declines but the checker accepts "
        f"{sorted(check_representation_change.NONE_VALUES)}"
    )


# ---------------------------------------------------------------------------
# A near miss must fail loudly, never read as silence
#
# #1838 is the live instance: its description declares on line 3 as
# `` `Representation-Change: none` `` — inside a code span. The strict parser
# finds zero trailers, git-cliff's footer rule will not match it, and
# `inject_representation_disclosure.py` will not attach it. Nothing was lost
# there, because the PR touches no watched directory and the value was a
# decline. The *shape* is the hazard: an author can believe they declared while
# every consumer of the trailer sees silence, and one backtick is the whole
# difference.
# ---------------------------------------------------------------------------


#: Every near-miss spelling, paired with the word the message must use to name it.
#: Each is a real Markdown habit rather than an invented one: a code span is how a
#: token gets quoted, `>`/`*`/`-`/`#` open the four block constructs that indent
#: nothing visibly, and the un-hyphenated spelling is what an author writes when
#: they are describing the field rather than copying it.
_NEAR_MISSES = [
    # #1838, verbatim.
    ("`Representation-Change: none`", "code span"),
    ("`Representation-Change: 3 rows of 500,004 move`", "code span"),
    ("**Representation-Change:** none", "emphasis"),
    ("*Representation-Change: none*", "emphasis"),
    ("_Representation-Change: none_", "emphasis"),
    ("> Representation-Change: none", "block quote"),
    ("- Representation-Change: none", "list marker"),
    ("* Representation-Change: none", "list marker"),
    ("+ Representation-Change: none", "list marker"),
    ("## Representation-Change: none", "heading"),
    ("    Representation-Change: none", "indent"),
    ("\tRepresentation-Change: none", "indent"),
    # The hyphen is what git-cliff's footer rule keys on, so a space or an
    # underscore is exactly as invisible as a backtick.
    ("Representation Change: none", "spelling"),
    ("Representation_Change: none", "spelling"),
    # Nothing between the token and its colon but a keystroke, and all three
    # consumers anchor on `Representation-Change:` with the colon flush — so these
    # are as unreadable as the rest. They matched the pattern from the start and
    # named no reason, which meant they were dropped without a word.
    ("Representation-Change : none", "separated from its colon"),
    ("Representation-Change\t: none", "separated from its colon"),
    ("Representation-Change*: none", "separated from its colon"),
]


@pytest.mark.parametrize(("line", "reason"), _NEAR_MISSES)
def test_a_near_miss_is_refused_and_the_message_names_it(line: str, reason: str) -> None:
    """A declaration nobody can read must fail, and say *why* it cannot be read.

    Reporting "no declaration" here is the failure this test exists to prevent: the
    author has written one, so the message they need is not "add a trailer" but "the
    trailer you added is invisible, and here is the character responsible".
    """
    ok, message = check(["src/normalize/merge.rs"], f"Some prose.\n\n{line}\n")
    assert not ok, f"{line!r} looks like a declaration and must not pass as silence"
    assert reason in message, f"the message must name the problem ({reason!r}); got:\n{message}"
    assert line.strip() in message, "the message must quote the offending line verbatim"


@pytest.mark.parametrize(("line", "reason"), _NEAR_MISSES)
def test_a_near_miss_is_refused_even_with_no_watched_file(line: str, reason: str) -> None:
    """The early return is the hole, so the check must sit in front of it.

    This is the same reasoning the two-trailer refusal already carries and is placed
    ahead of the watched-file test for: the harm is in what the *changelog* renders,
    which applies to every commit. #1838 is green today only because it touches no
    watched directory — it reaches the early return before the trailer would have
    mattered — so a near-miss check behind that return would rebuild the same hole and
    pass the one PR that motivated it.

    It does make the check stricter on changes it currently ignores entirely. The
    narrower alternative was measured and rejected: over all 39 open PRs it would
    report nothing at all.
    """
    ok, message = check(["README.md", "docs/x.md"], f"Some prose.\n\n{line}\n")
    assert not ok, f"{line!r} must be refused even where no declaration is required"
    assert reason in message


def test_find_near_misses_reports_the_line_and_its_reason() -> None:
    body = "Prose.\n\n`Representation-Change: none`\n"
    misses = find_near_misses(body)
    assert len(misses) == 1
    (line_number, text, reason) = misses[0]
    assert line_number == 3
    assert text == "`Representation-Change: none`"
    assert "code span" in reason


def test_every_near_miss_shape_names_a_reason() -> None:
    """No line the pattern matches may be reported as silence.

    `_NEAR_MISSES` above is a list of the shapes somebody thought of, so it cannot say
    anything about the ones nobody did. This drives the *cross product* of the decorations
    the pattern admits and asserts the invariant directly: a line that
    `NEAR_MISS_RE` matches, and that `TRAILER_RE` does not, is reported with a non-empty
    reason.

    It is not hypothetical. `Representation-Change : none` — one space, hyphen intact, no
    lead — matched the pattern and produced no reason, so it was dropped, which reported an
    unreadable declaration as silence inside the scan written to stop exactly that. Three
    shapes did. A list-based test could only have caught them by happening to list them.
    """
    leads = ["", "  ", "\t", "> ", "- ", "* ", "+ ", "# ", "`", "**", "_"]
    separators = ["-", "_", " ", ""]
    gaps = ["", " ", "\t", "*", "`", "_"]
    tails = ["none", "3 rows of 500,004 move", "none`", "none**"]

    checked = 0
    for lead in leads:
        for separator in separators:
            for gap in gaps:
                for tail in tails:
                    line = f"{lead}Representation{separator}Change{gap}: {tail}"
                    if check_representation_change.TRAILER_RE.match(line):
                        continue  # a real trailer, not a near miss
                    if not check_representation_change.NEAR_MISS_RE.match(line):
                        continue  # the pattern's own business, and not this contract's
                    misses = find_near_misses(f"Prose.\n\n{line}\n")
                    assert misses, f"{line!r} matches the pattern and was reported as silence"
                    assert misses[0][2].strip(), f"{line!r} was reported with an empty reason"
                    checked += 1

    # Anti-vacuity: a pattern change that stopped matching everything would otherwise leave
    # this test green while asserting nothing at all.
    assert checked > 100, f"only {checked} shapes reached the assertion; the battery is hollow"


# --- and now the half that keeps it from becoming a nuisance ---------------
#
# Every case below is drawn from a real open PR body. Measured over all 39 open
# PRs on 2026-08-13, this rule fires on exactly one of them (#1838); without
# either discriminator below it fires on five. The three extra are not defects
# and failing them would make the check something contributors route around.


def test_a_valid_trailer_silences_the_near_miss_scan() -> None:
    """An author who declared at column 0 has demonstrated they know the form.

    Real instances: #1742 line 5 (`` `Representation-Change: none` is unaffected: no
    normalized description moves.``) and #1837 line 84, both *discussing* their own
    declaration in prose while carrying it properly above. Firing on those would
    punish a PR for explaining itself, which is the #1555 mistake in a new place.

    It also preserves the documented escape hatch: `CONTRIBUTING.md` says to indent a
    quoted example, and the two-trailer refusal's own message recommends exactly that.
    """
    body = (
        "Representation-Change: 577 rows move, 360 merge.\n"
        "\n"
        "`Representation-Change: none` is what the declining form looks like.\n"
    )
    assert find_near_misses(body) == [], "a column-0 trailer must silence the scan"
    ok, _ = check(["src/normalize/merge.rs"], body)
    assert ok


@pytest.mark.parametrize(
    "body",
    [
        # #1753, the release PR, verbatim — both of its reviewer-checklist mentions.
        "Any PR that moves a normalized output should carry a `Representation-Change:`\ntrailer.",
        "- [ ] Nothing in it merely *declined*. A decline reads `Representation-Change:\n  none`.",
        # #1752 line 136, verbatim in shape.
        "Both closures are tests only, so the `Representation-Change:` disclosure above\nis unaffected.",
        # The checker's own failure message, pasted into a description.
        "The check said to add a `Representation-Change:` trailer.",
    ],
)
def test_a_mid_sentence_mention_is_not_a_near_miss(body: str) -> None:
    """A declaration is a *line*, not a phrase, so the token must open the line.

    This is the discriminator that does the real work: measured over all 39 open PRs, the
    anchor alone excludes every mention that is not a declaration. All four bodies here are
    real, and firing on them would fail a PR for naming the mechanism it is complying with.
    """
    assert find_near_misses(body) == [], f"{body!r} mentions the field, it does not declare"
    ok, _ = check(["README.md"], body)
    assert ok


def test_a_line_anchored_mention_with_no_value_is_not_a_near_miss() -> None:
    """The value requirement, isolated — the anchor cannot save this one.

    Kept deliberately, and deliberately not claimed as measured-necessary: dropping it
    reports the same single PR over today's 39. What it buys is that naming the field at
    the start of a line — which is how a checklist or a doc edit refers to it — is not
    reported as a declaration that failed, because there is no value to lose.
    """
    assert find_near_misses("`Representation-Change:`\n") == []
    assert find_near_misses("> Representation-Change:\n") == []


def test_an_empty_trailer_is_still_reported_as_absent_not_as_a_near_miss() -> None:
    """`Representation-Change:` at column 0 with no value is the #1522 hole itself, and
    it already has its own verdict. Re-filing it as a near miss would change a settled
    message for no gain."""
    assert find_near_misses("Representation-Change:   \n") == []
    ok, message = check(["src/normalize/merge.rs"], "Representation-Change:   \n")
    assert not ok
    assert "can move normalized output" in message, "the absent-trailer message must survive"


def test_the_1838_body_is_refused_by_the_check_it_slipped_past() -> None:
    """End-to-end on the shape that motivated this, with #1838's real file list.

    It touches one fixture JSON and no watched directory, so it reaches the early
    return today and passes. The point of the fix is that it must not.
    """
    body = (
        "Two rationale-only corrections to decided records, in one PR because they "
        "touch the same file.\n"
        "\n"
        "`Representation-Change: none`\n"
        "\n"
        "Measured, not assumed: the ledger is read only at generation time.\n"
    )
    assert find_declarations(body) == [], "the strict parser sees nothing — that is the defect"
    ok, message = check(["tests/fixtures/grammar/hgvs_spec_normalization_overrides.json"], body)
    assert not ok, "a trailer no consumer can read must not pass as a considered declaration"
    assert "code span" in message
    assert "`Representation-Change: none`" in message


def test_the_decline_exclusion_precedes_the_inclusion() -> None:
    """git-cliff takes the FIRST matching parser, so swapping these two lines silently
    restores #1526 in full.

    Verified against git-cliff 2.13.1 rather than assumed: with the exclusion second, all
    four trailers -- `none`, `NONE`, and both real disclosures -- group under
    **Representation changes**. Every other guard in this file still passes in that state,
    which is why the ordering needs one of its own.
    """
    config = (Path(__file__).resolve().parents[2] / "release-plz.toml").read_text(encoding="utf-8")
    rules = re.findall(r'\{ footer = "([^"]*Representation-Change[^"]*)"', config)
    assert len(rules) == 2, f"expected an exclusion and an inclusion rule, found {rules}"
    exclusion, inclusion = rules
    assert any(word in exclusion for word in check_representation_change.NONE_VALUES), (
        f"the first Representation-Change rule is {exclusion!r}, which does not name a decline "
        "value; the exclusion must come first or every decline is grouped as a real change"
    )
    assert not any(word in inclusion for word in check_representation_change.NONE_VALUES), (
        f"the second rule is {inclusion!r}, which looks like the exclusion -- the two are reversed"
    )


# ---------------------------------------------------------------------------
# A fenced trailer is documentation, not a declaration (#1929)
# ---------------------------------------------------------------------------


def test_a_fenced_trailer_is_not_a_declaration() -> None:
    """The #1929 reproducer: quoting the docs must not file a disclosure.

    This is the one decoration that reached column 0 with its text intact, so unlike an
    inline code span or a blockquote it did not merely go unread — it was read as REAL, and
    the author was recorded as disclosing a 577-row migration they said in the same body
    they had not worked out. A false disclosure, not a missed one.
    """
    body = (
        "This PR refactors a helper.\n\n"
        "For reference, CONTRIBUTING.md says the trailer looks like:\n\n"
        "```\n"
        "Representation-Change: 577 rows move, 360 merge / 205 split / 12 respell.\n"
        "  Previously-accepted inputs, so a real migration for the consumer.\n"
        "```\n\n"
        "I have not worked out the disclosure for this change yet.\n"
    )
    assert find_declaration(body) is None
    passed, message = check(["src/normalize/mod.rs"], body)
    assert not passed
    assert "fenced code block" in message


def test_a_fenced_trailer_is_reported_as_a_near_miss_not_as_silence() -> None:
    """Someone who pasted the documented form is asking a question, not staying silent."""
    body = "x\n\n```\nRepresentation-Change: none\n```\n"
    misses = find_near_misses(body)
    assert len(misses) == 1
    line_number, line, reason = misses[0]
    assert line_number == 4
    assert line.startswith("Representation-Change:")
    assert "fenced" in reason and "column 0" in reason


@pytest.mark.parametrize(
    "opener,closer",
    [("```", "```"), ("~~~", "~~~"), ("````", "````"), ("```rust", "```")],
)
def test_every_fence_spelling_hides_a_trailer(opener: str, closer: str) -> None:
    body = f"x\n\n{opener}\nRepresentation-Change: 9 rows move\n{closer}\n"
    assert find_declaration(body) is None


def test_a_longer_fence_may_contain_a_shorter_one() -> None:
    """The shape the issue reproduced with: a ```` block quoting a ``` block."""
    body = "x\n\n````markdown\n```\nRepresentation-Change: 577 rows move\n```\n````\n"
    assert find_declaration(body) is None


def test_a_backtick_opener_whose_info_carries_a_backtick_is_not_a_fence() -> None:
    """CommonMark forbids a backtick in a BACKTICK fence's info string; tildes allow it.

    Without that restriction an inline code span opens a block, so a line like
    ```` ```see `none` below ```` reads as an opener and blanks everything after it — the
    real declaration below included. It fails closed (the author is refused rather than
    credited with a disclosure they never wrote), which is the safe direction and still
    the wrong answer.

    Both halves are asserted because only together do they show the rule is about the
    fence CHARACTER and not about backticks in info strings generally.
    """
    body = "```see `none` below\n\nRepresentation-Change: none. Tests only.\n"
    assert find_declaration(body) == "none. Tests only."
    passed, _ = check(["src/normalize/mod.rs"], body)
    assert passed

    # A tilde fence has no such restriction, so this one genuinely opens and swallows it.
    tilde = "~~~see `none` below\n\nRepresentation-Change: none. Tests only.\n"
    assert find_declaration(tilde) is None


def test_an_unclosed_fence_runs_to_the_end() -> None:
    """CommonMark's rule, and the safe direction: withhold rather than invent."""
    body = "x\n\n```\nRepresentation-Change: 577 rows move\n"
    assert find_declaration(body) is None


def test_a_real_trailer_after_a_closed_fence_still_counts() -> None:
    """The fence must not swallow the rest of the body."""
    body = "```\nsome code\n```\n\nRepresentation-Change: none. Tests only.\n"
    assert find_declaration(body) == "none. Tests only."
    passed, _ = check(["src/normalize/mod.rs"], body)
    assert passed


def test_a_fenced_example_beside_a_real_trailer_is_still_refused() -> None:
    """The safety property that makes the fence-aware set the WRONG one to count.

    `git_conventional` has no Markdown awareness — any `word:` line is a footer token — so
    git-cliff sees the fenced example too. Were the duplicate refusal made fence-aware, this
    body would pass here while the changelog still had two footers to choose from, and the
    quoted 577-row example could land under an author who declared `none`.

    So a fenced trailer may not BE the declaration and still COUNTS toward how many the
    changelog will see. The refusal is what forces the author to indent or drop the quote.
    """
    body = (
        "```\n"
        "Representation-Change: 577 rows move\n"
        "```\n\n"
        "Representation-Change: none. Tests only.\n"
    )
    assert find_declaration(body) == "none. Tests only."
    assert len(find_trailers_as_git_cliff_would(body)) == 2
    passed, message = check(["src/normalize/mod.rs"], body)
    assert not passed
    assert "trailers found" in message


def test_two_fenced_examples_and_no_real_trailer_are_refused_naming_the_fence() -> None:
    """The one path where the fenced case does not get the fenced message.

    `check` counts trailers as git-cliff would — fences included, deliberately — and
    refuses on more than one BEFORE the near-miss scan runs. When *every* counted trailer
    is fenced, that refusal is the one the author meets, and its advice ("Keep exactly
    one", "Delete the superseded trailer") names a construct they never wrote: a Markdown
    reader sees two code blocks. The advice happens to work, since indenting inside the
    fence drops the line from `TRAILER_RE`, but it sends the author looking for a
    declaration they have not made.

    The body is not contrived. `CONTRIBUTING.md` publishes three such blocks and the whole
    premise of #1929 is someone quoting the docs to ask what to write; quoting two of them
    is how you ask WHICH form applies.

    Both halves of the refusal are asserted, because both are true and neither alone
    explains it: the count is what the changelog will see, and the fence is why the author
    cannot see it.
    """
    body = (
        "Which of these should I write?\n\n"
        "```\n"
        "Representation-Change: none\n"
        "```\n\n"
        "or\n\n"
        "```\n"
        "Representation-Change: 577 rows move, 360 merge / 205 split / 12 respell.\n"
        "```\n"
    )
    assert find_declaration(body) is None, "neither example is readable, so nothing is declared"
    assert len(find_trailers_as_git_cliff_would(body)) == 2, "git-cliff counts both"

    passed, message = check(["src/normalize/mod.rs"], body)
    assert not passed
    assert "trailers found" in message, "the duplicate count is still reported"
    assert "fenced code block" in message, (
        "the duplicate refusal fired without naming the fence, so the author is told to "
        f"delete a trailer that is not one; message was:\n{message}"
    )
    assert "nothing to delete here" in message


def test_the_fenced_note_is_absent_when_a_real_trailer_is_among_the_duplicates() -> None:
    """The negative control for the arm above: a readable trailer means real duplicates.

    Here `Delete the superseded trailer` is exactly right, so the fenced note must not be
    appended — it would tell an author with a genuine duplicate that there is nothing to
    delete.
    """
    body = (
        "```\n"
        "Representation-Change: 577 rows move\n"
        "```\n\n"
        "Representation-Change: none. Tests only.\n"
    )
    passed, message = check(["src/normalize/mod.rs"], body)
    assert not passed
    assert "trailers found" in message
    assert "nothing to delete here" not in message


def test_contributing_md_s_own_examples_are_not_read_as_declarations() -> None:
    """Run the checker over the contributor documentation itself.

    The issue's suggested guard, and it is stronger than a hand-copied fixture: it keeps
    working when the docs are rewritten. `CONTRIBUTING.md` publishes the trailer's form in
    fenced blocks, so before #1929 its own text — pasted by someone asking what to write —
    filed a disclosure.
    """
    contributing = (_REPO_ROOT / "CONTRIBUTING.md").read_text(encoding="utf-8")
    assert "Representation-Change:" in contributing, "the docs stopped documenting the trailer"
    fenced = fenced_line_numbers(contributing)
    quoted_in_a_fence = [
        number
        for number, line in enumerate(contributing.splitlines(), 1)
        if line.startswith("Representation-Change:") and number in fenced
    ]
    assert quoted_in_a_fence, "expected CONTRIBUTING.md to publish the form inside a fence"
    assert find_declaration(contributing) is None, (
        "CONTRIBUTING.md's own examples are read as a declaration; "
        f"fenced trailer lines were {quoted_in_a_fence}"
    )
