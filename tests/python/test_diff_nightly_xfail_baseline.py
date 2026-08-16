"""Unit tests for `scripts/diff_nightly_xfail_baseline.py` (#1998).

The script is the gating signal for `nightly-mutalyzer.yml`: it diffs the
reference-aware run's actual FAIL set (parsed from nextest JUnit) against a
committed baseline of known-failing test ids, and turns the nightly RED when the
two disagree in *either* direction. What is worth pinning is exactly that gate
logic -- the three outcomes and their exit codes -- because each is silent if
wrong:

- an **exact match** must stay green (exit 0), or the ~38 known-xfail cases turn
  the job red every night and the signal is worthless;
- a **new failure** (failing now, not baselined) must go red (exit 1), which is
  the whole point of the baseline;
- an **unexpected pass** (baselined, now passing or no longer observed) must also
  go red (exit 1), so the baseline cannot silently rot away from reality.

The parser is tested against synthetic JUnit rather than a captured nextest
report so the cases stay hermetic and readable.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "diff_nightly_xfail_baseline.py"


def _load():
    """Import the script by path; `scripts/` is not a package."""
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("diff_nightly_xfail_baseline", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    # Register before exec: the script's frozen dataclasses resolve their
    # (stringized) field annotations against `sys.modules[__module__]`, which
    # must exist by the time `@dataclass` runs.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


mod = _load()


def _junit(cases: list[tuple[str, str, str]]) -> str:
    """Build a JUnit document from ``(classname, name, outcome)`` tuples.

    ``outcome`` is one of ``"pass"``, ``"fail"``, ``"error"``, ``"skip"``.
    """
    body = []
    for classname, name, outcome in cases:
        child = {
            "fail": "<failure message='boom'/>",
            "error": "<error message='boom'/>",
            "skip": "<skipped/>",
            "pass": "",
        }[outcome]
        body.append(f'<testcase classname="{classname}" name="{name}">{child}</testcase>')
    inner = "".join(body)
    return f'<testsuites><testsuite name="nextest" tests="{len(cases)}">{inner}</testsuite></testsuites>'


# --------------------------------------------------------------------------- #
# parsing                                                                      #
# --------------------------------------------------------------------------- #


def test_test_id_joins_classname_and_name() -> None:
    assert mod.test_id("ferro-hgvs::it", "mod::a") == "ferro-hgvs::it::mod::a"


def test_test_id_falls_back_to_bare_name() -> None:
    assert mod.test_id("", "mod::a") == "mod::a"


def test_parse_junit_partitions_outcomes() -> None:
    xml = _junit(
        [
            ("bin", "mod::pass_one", "pass"),
            ("bin", "mod::fail_one", "fail"),
            ("bin", "mod::error_one", "error"),
            ("bin", "mod::skip_one", "skip"),
        ]
    )
    results = mod.parse_junit(xml)
    assert results.failed == frozenset({"bin::mod::fail_one", "bin::mod::error_one"})
    assert results.passed == frozenset({"bin::mod::pass_one"})
    assert results.skipped == frozenset({"bin::mod::skip_one"})
    assert "bin::mod::skip_one" in results.observed


def test_parse_junit_raises_on_garbage() -> None:
    with pytest.raises(ValueError):
        mod.parse_junit("not xml <<<")


def test_read_baseline_ignores_comments_and_blanks() -> None:
    text = "# a comment\n\nbin::mod::a\n  bin::mod::b  \n"
    assert mod.read_baseline(text) == frozenset({"bin::mod::a", "bin::mod::b"})


# --------------------------------------------------------------------------- #
# baseline strictness                                                          #
#                                                                              #
# Every shape below used to shrink the baseline in silence, and a silently      #
# smaller baseline is a silently weaker gate: an id that drops out stops being  #
# expected-to-fail, so when it stops failing nothing reports it. Each is now a  #
# `ValueError`, which `main()` turns into exit 2.                               #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "candidate",
    [
        "ferro-hgvs::it::mod::a_test",
        "ferro-hgvs::it::a_test",
        "ferro-hgvs::tests::a_test",
    ],
)
def test_is_test_id_accepts_real_nextest_ids(candidate: str) -> None:
    assert mod.is_test_id(candidate)


@pytest.mark.parametrize(
    "candidate",
    [
        "",
        "bare",
        "two::segments",
        "bin::mod::a  # flaky",
        "bin::mod::a\tb",
        "bin::mod::#a",
        "bin::::a",
        "::mod::a",
        "bin::mod::",
    ],
)
def test_is_test_id_rejects_malformed(candidate: str) -> None:
    assert not mod.is_test_id(candidate)


def test_read_baseline_rejects_an_inline_comment() -> None:
    # The whole line used to be taken as the id, comment and all, so it matched
    # no test in any run -- reported forever as a baselined test that stopped
    # failing, or, if someone re-blessed past it, quietly dropped.
    with pytest.raises(ValueError, match="line 2"):
        mod.read_baseline("bin::mod::a\nbin::mod::b  # known flaky\n")


def test_read_baseline_rejects_a_commented_out_id() -> None:
    # Indistinguishable from an ordinary comment to a parser that skips every
    # `#` line; the committed file carries no comments, so the ambiguity is
    # resolved by refusing rather than by guessing.
    with pytest.raises(ValueError, match="commented-out"):
        mod.read_baseline("bin::mod::a\n# bin::mod::b\n")


def test_read_baseline_still_ignores_a_prose_comment() -> None:
    # Only a comment whose text is a well-formed id is ambiguous. Prose is not.
    assert mod.read_baseline("# regenerate with --update-baseline\nbin::mod::a\n") == frozenset(
        {"bin::mod::a"}
    )


def test_read_baseline_rejects_a_duplicate() -> None:
    with pytest.raises(ValueError, match="repeats"):
        mod.read_baseline("bin::mod::a\nbin::mod::a\n")


def test_main_malformed_baseline_exits_two(tmp_path: Path) -> None:
    # The exit code matters as much as the raise: the workflow step now runs
    # under `set -euo pipefail`, so a 2 fails the job instead of being swallowed
    # by the `tee` the way the 1 was (#2082).
    junit = tmp_path / "junit.xml"
    junit.write_text(_junit([("bin", "mod::a", "fail")]), encoding="utf-8")
    baseline = tmp_path / "baseline.txt"
    baseline.write_text("bin::mod::a  # oops\n", encoding="utf-8")
    assert mod.main(["--junit", str(junit), "--baseline", str(baseline)]) == 2


# --------------------------------------------------------------------------- #
# the COMMITTED baseline                                                       #
#                                                                              #
# Nothing else in PR CI reads this file -- the nightly is the only consumer,    #
# and it runs on a schedule against a reference no pull request has. So a       #
# baseline broken in a PR is invisible until the next night, on a job that      #
# opens an issue when it fails. These checks are cheap and hermetic.            #
# --------------------------------------------------------------------------- #

COMMITTED_BASELINE = (
    REPO_ROOT / "tests" / "fixtures" / "nightly-reference-aware" / "xfail-baseline.txt"
)


def test_the_committed_baseline_parses_and_is_not_empty() -> None:
    ids = mod.read_baseline(COMMITTED_BASELINE.read_text(encoding="utf-8"))
    assert ids, (
        f"{COMMITTED_BASELINE} parsed to an empty set. An empty baseline makes every one of the "
        "nightly's known failures a NEW FAILURE, so the job is red every night and the signal is "
        "worthless -- and it is what a parser that has stopped matching also produces."
    )
    for tid in sorted(ids):
        assert mod.is_test_id(tid), f"{tid!r} is not a well-formed nextest test id"


def test_the_committed_baseline_is_exactly_what_update_baseline_writes() -> None:
    # One property covering sorted, unique, no blank lines, no trailing
    # whitespace and no comments: the file must be byte-identical to what
    # `--update-baseline` would render from the set it parses to. That makes a
    # hand-edit -- which the README warns against, and which is how the two
    # silent-shrink shapes above get introduced -- visible in PR CI rather than
    # on the next night's run.
    text = COMMITTED_BASELINE.read_text(encoding="utf-8")
    assert text == mod.render_baseline(mod.read_baseline(text)), (
        f"{COMMITTED_BASELINE} is not in the canonical form `--update-baseline` writes "
        "(sorted, unique, one id per line, no comments or blanks). Annotations belong in that "
        "directory's README.md -- a regeneration would delete them from here anyway."
    )


def test_update_baseline_round_trips_what_it_writes(tmp_path: Path) -> None:
    junit, baseline = _write(
        tmp_path,
        [("ferro-hgvs::it", "mod::b", "fail"), ("ferro-hgvs::it", "mod::a", "fail")],
        [],
    )
    assert mod.run(junit, baseline, update=True) == 0
    written = baseline.read_text(encoding="utf-8")
    assert mod.read_baseline(written) == frozenset(
        {"ferro-hgvs::it::mod::a", "ferro-hgvs::it::mod::b"}
    )


def test_render_baseline_is_sorted_with_trailing_newline() -> None:
    rendered = mod.render_baseline(frozenset({"bin::c", "bin::a", "bin::b"}))
    assert rendered == "bin::a\nbin::b\nbin::c\n"


# --------------------------------------------------------------------------- #
# diff logic                                                                   #
# --------------------------------------------------------------------------- #


def _results(failed: set[str], passed: set[str] = frozenset(), skipped: set[str] = frozenset()):
    return mod.JUnitResults(
        failed=frozenset(failed),
        passed=frozenset(passed),
        skipped=frozenset(skipped),
    )


def test_exact_match_is_not_drift() -> None:
    baseline = frozenset({"bin::a", "bin::b"})
    results = _results(failed={"bin::a", "bin::b"}, passed={"bin::c"})
    diff = mod.diff_baseline(results, baseline)
    assert not diff.is_drift
    assert diff.new_failures == frozenset()
    assert diff.unexpected_passes == frozenset()


def test_new_failure_is_drift() -> None:
    baseline = frozenset({"bin::a"})
    results = _results(failed={"bin::a", "bin::NEW"})
    diff = mod.diff_baseline(results, baseline)
    assert diff.is_drift
    assert diff.new_failures == frozenset({"bin::NEW"})
    assert diff.unexpected_passes == frozenset()


def test_unexpected_pass_when_baselined_test_passes() -> None:
    baseline = frozenset({"bin::a", "bin::b"})
    results = _results(failed={"bin::a"}, passed={"bin::b"})
    diff = mod.diff_baseline(results, baseline)
    assert diff.is_drift
    assert diff.now_passing == frozenset({"bin::b"})
    assert diff.absent == frozenset()
    assert diff.unexpected_passes == frozenset({"bin::b"})


def test_baselined_but_not_observed_is_drift() -> None:
    baseline = frozenset({"bin::a", "bin::gone"})
    results = _results(failed={"bin::a"})
    diff = mod.diff_baseline(results, baseline)
    assert diff.is_drift
    assert diff.absent == frozenset({"bin::gone"})
    assert diff.unexpected_passes == frozenset({"bin::gone"})


def test_baselined_but_now_skipped_is_drift() -> None:
    # A known-failing test that becomes #[ignore]d stops failing without passing;
    # it must still count as rot, or the baseline can silently drift.
    baseline = frozenset({"bin::a", "bin::ign"})
    results = _results(failed={"bin::a"}, skipped={"bin::ign"})
    diff = mod.diff_baseline(results, baseline)
    assert diff.is_drift
    assert diff.now_skipped == frozenset({"bin::ign"})
    assert diff.now_passing == frozenset()
    assert diff.absent == frozenset()
    assert diff.unexpected_passes == frozenset({"bin::ign"})


# --------------------------------------------------------------------------- #
# end-to-end run() exit codes                                                  #
# --------------------------------------------------------------------------- #


def _write(tmp_path: Path, junit_cases, baseline_ids: list[str]):
    junit = tmp_path / "junit.xml"
    junit.write_text(_junit(junit_cases), encoding="utf-8")
    baseline = tmp_path / "baseline.txt"
    baseline.write_text("".join(f"{i}\n" for i in baseline_ids), encoding="utf-8")
    return junit, baseline


def test_run_exact_match_exits_zero(tmp_path: Path) -> None:
    junit, baseline = _write(
        tmp_path,
        [("bin", "mod::a", "fail"), ("bin", "mod::b", "pass")],
        ["bin::mod::a"],
    )
    assert mod.run(junit, baseline, update=False) == 0


def test_run_new_failure_exits_one(tmp_path: Path) -> None:
    junit, baseline = _write(
        tmp_path,
        [("bin", "mod::a", "fail"), ("bin", "mod::b", "fail")],
        ["bin::mod::a"],
    )
    assert mod.run(junit, baseline, update=False) == 1


def test_run_unexpected_pass_exits_one(tmp_path: Path) -> None:
    junit, baseline = _write(
        tmp_path,
        [("bin", "mod::a", "pass")],
        ["bin::mod::a"],
    )
    assert mod.run(junit, baseline, update=False) == 1


def test_run_baselined_now_skipped_exits_one(tmp_path: Path) -> None:
    junit, baseline = _write(
        tmp_path,
        [("bin", "mod::a", "skip")],
        ["bin::mod::a"],
    )
    assert mod.run(junit, baseline, update=False) == 1


def test_update_baseline_rewrites_from_run(tmp_path: Path) -> None:
    junit, baseline = _write(
        tmp_path,
        [("bin", "mod::a", "fail"), ("bin", "mod::b", "fail"), ("bin", "mod::c", "pass")],
        ["stale::entry"],
    )
    assert mod.run(junit, baseline, update=True) == 0
    assert baseline.read_text(encoding="utf-8") == "bin::mod::a\nbin::mod::b\n"


def test_main_missing_junit_exits_two(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline.txt"
    baseline.write_text("bin::mod::a\n", encoding="utf-8")
    rc = mod.main(["--junit", str(tmp_path / "nope.xml"), "--baseline", str(baseline)])
    assert rc == 2
