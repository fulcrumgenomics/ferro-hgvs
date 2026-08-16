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
