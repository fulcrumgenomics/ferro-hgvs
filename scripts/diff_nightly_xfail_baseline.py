#!/usr/bin/env python3
"""Diff the nightly reference-aware run's FAIL set against a committed baseline.

`nightly-mutalyzer.yml` runs every reference-aware guard against a real
ferro-prepared manifest. ~38 of those tests are known-failing (live divergences
from mutalyzer / biocommons / hgvs-rs, or tracked ferro bugs) and the job is
deliberately *non-gating* on the failures themselves — see #1998 and the
`continue-on-error` comment in that workflow.

The gap #1998 records is that nothing compared that FAIL set against anything, so
a **new** failure was indistinguishable from the ~38 beside it and the job went
green either way. This script closes that gap. It reads the JUnit report nextest
emits under the `nightly-xfail` profile, extracts the exact set of failing test
ids, and diffs it against the committed baseline:

- **new failure** — a test failing now that the baseline does not list. This is
  drift (a regression, or an upstream reference change) and turns the job RED.
- **unexpected pass** — a baselined test that is no longer failing (it passed, or
  is no longer present in the run). The baseline has rotted and must be updated;
  this also turns the job RED, so the baseline cannot silently drift away from
  reality.
- **exact match** — the FAIL set equals the baseline. Green.

Regenerating the baseline after a deliberate change (a burn-down PR that fixes
tests, or an upstream reference refresh) is `--update-baseline`, which rewrites
the baseline file from the current JUnit report.

A test id is `"{classname}::{name}"` as nextest records them in JUnit —
`classname` is the binary id (e.g. `ferro-hgvs::it`) and `name` is the full
module path plus test function (e.g. `mutalyzer_normalize_tests::axis_normalized`).
Including the binary id keeps lib-unit tests (`normalize_tests::…::tests`) and
integration tests distinct.
"""

from __future__ import annotations

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


def test_id(classname: str, name: str) -> str:
    """Canonical id for a JUnit testcase: ``"{classname}::{name}"``.

    ``classname`` may be empty for a harness that does not set it; fall back to
    the bare name so an id is always non-empty.
    """
    classname = classname.strip()
    name = name.strip()
    if classname:
        return f"{classname}::{name}"
    return name


@dataclass(frozen=True)
class JUnitResults:
    """The three disjoint outcome sets parsed from a JUnit report."""

    failed: frozenset[str]
    passed: frozenset[str]
    skipped: frozenset[str]

    @property
    def observed(self) -> frozenset[str]:
        """Every test id the report mentions, regardless of outcome."""
        return self.failed | self.passed | self.skipped


def parse_junit(xml_text: str) -> JUnitResults:
    """Parse a nextest JUnit report into failed / passed / skipped id sets.

    A ``<testcase>`` is *failed* if it carries a ``<failure>`` or ``<error>``
    child, *skipped* if it carries a ``<skipped>`` child, and *passed*
    otherwise. Raises ``ValueError`` if the XML does not parse.
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:  # noqa: TRY003
        raise ValueError(f"could not parse JUnit XML: {exc}") from exc

    failed: set[str] = set()
    passed: set[str] = set()
    skipped: set[str] = set()
    for case in root.iter("testcase"):
        tid = test_id(case.get("classname", ""), case.get("name", ""))
        if case.find("failure") is not None or case.find("error") is not None:
            failed.add(tid)
        elif case.find("skipped") is not None:
            skipped.add(tid)
        else:
            passed.add(tid)
    return JUnitResults(
        failed=frozenset(failed),
        passed=frozenset(passed),
        skipped=frozenset(skipped),
    )


#: The fewest ``::``-separated segments a nextest id can have:
#: ``{binary-id}::{module-path}::{test-fn}``. Three is the real floor, not a
#: round number — a test function directly in an integration target's root is
#: ``ferro-hgvs::it::a_test``, and a unit test in a crate-root ``mod tests`` is
#: ``ferro-hgvs::tests::a_test``.
MINIMUM_ID_SEGMENTS = 3


def is_test_id(candidate: str) -> bool:
    """Whether ``candidate`` is a well-formed nextest test id.

    ``{binary-id}::{module-path}::{test-fn}``: at least
    ``MINIMUM_ID_SEGMENTS`` non-empty ``::``-separated segments, with no
    whitespace and no ``#`` anywhere. Those two characters are what make a
    baseline line ambiguous between an id, a comment, and an id with a comment
    stuck to it — which is the whole of what this predicate exists to refuse.
    """
    if not candidate or "#" in candidate or any(c.isspace() for c in candidate):
        return False
    segments = candidate.split("::")
    if len(segments) < MINIMUM_ID_SEGMENTS:
        return False
    return all(segment and ":" not in segment for segment in segments)


def read_baseline(text: str) -> frozenset[str]:
    """Parse a baseline file: one test id per line, ``#`` comments and blanks ignored.

    **Strict, and deliberately so.** Every way this parser could once shrink the
    baseline it was silent about, and a silently smaller baseline is a silently
    weaker gate: an id that drops out stops being expected-to-fail, so when it
    stops failing nothing notices. The three shapes, all of which now raise
    ``ValueError`` naming the line:

    - an **inline comment** (``…::a_test  # flaky``) — the whole line used to be
      taken as the id, comment and all, so it matched no test ever;
    - a **commented-out id** (``# ferro-hgvs::it::m::t``) — indistinguishable
      from an ordinary comment to a parser that skips every ``#`` line, and the
      committed file carries no comments at all, so the ambiguity is resolved by
      refusing rather than by guessing;
    - a **duplicate**, which the set collapsed without comment.

    Blank lines and comment lines whose text is *not* a test id are still
    ignored; the file's format has not changed, only what it will accept in
    silence.

    Raises:
        ValueError: on any line that is neither blank, nor a plain comment, nor
            a well-formed and not-already-seen test id.
    """
    ids: list[str] = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            commented = line.lstrip("#").strip()
            if is_test_id(commented):
                raise ValueError(
                    f"baseline line {number} is a commented-out test id: {line!r}. "
                    "A baseline is a set of expected failures; commenting one out removes it "
                    "silently, so the test is no longer expected to fail and nothing reports "
                    "it when it stops. Delete the line, or keep it and explain why in "
                    "tests/fixtures/nightly-reference-aware/README.md."
                )
            continue
        if not is_test_id(line):
            raise ValueError(
                f"baseline line {number} is not a well-formed test id: {line!r}. "
                "Expected `{binary-id}::{module-path}::{test-fn}` with no whitespace and no "
                "`#` — an inline comment makes the line match no test at all, which reads as "
                "a baselined test that has stopped failing."
            )
        if line in ids:
            raise ValueError(
                f"baseline line {number} repeats {line!r}. The baseline is a set; a duplicate "
                "is collapsed on read, so the file and the gate disagree about its size."
            )
        ids.append(line)
    return frozenset(ids)


def render_baseline(fail_ids: frozenset[str]) -> str:
    """Render a baseline file body: sorted ids, one per line, trailing newline."""
    return "".join(f"{tid}\n" for tid in sorted(fail_ids))


@dataclass(frozen=True)
class Diff:
    """The outcome of comparing a run's FAIL set against the baseline."""

    new_failures: frozenset[str]
    now_passing: frozenset[str]
    now_skipped: frozenset[str]
    absent: frozenset[str]

    @property
    def is_drift(self) -> bool:
        """True when the run diverges from the baseline in either direction."""
        return bool(self.new_failures or self.unexpected_passes)

    @property
    def unexpected_passes(self) -> frozenset[str]:
        """Every baselined test that is no longer failing — the baseline has rotted.

        This is exactly ``baseline - failed``: a baselined test that passed, was
        skipped (e.g. newly ``#[ignore]``d), or is no longer present. All three
        mean the committed set no longer describes reality, so all three are
        drift — a baselined test silently ceasing to fail is precisely the rot
        the baseline exists to prevent.
        """
        return self.now_passing | self.now_skipped | self.absent


def diff_baseline(results: JUnitResults, baseline: frozenset[str]) -> Diff:
    """Diff a run's outcomes against the baseline of known-failing test ids.

    - ``new_failures`` — failing now, not in the baseline.
    - ``now_passing``  — baselined, but observed passing in this run.
    - ``now_skipped``  — baselined, but skipped in this run (no longer failing).
    - ``absent``       — baselined, but not observed at all (removed / renamed).

    The three latter sets partition ``baseline - failed`` (the outcome sets are
    disjoint), and together they are ``unexpected_passes``: any baselined test
    not failing now is drift, whichever way it stopped.
    """
    new_failures = results.failed - baseline
    now_passing = baseline & results.passed
    now_skipped = baseline & results.skipped
    absent = baseline - results.observed
    return Diff(
        new_failures=frozenset(new_failures),
        now_passing=frozenset(now_passing),
        now_skipped=frozenset(now_skipped),
        absent=frozenset(absent),
    )


def _annotate(level: str, message: str) -> None:
    """Emit a GitHub Actions workflow annotation when running under Actions."""
    if os.environ.get("GITHUB_ACTIONS") == "true":
        print(f"::{level}::{message}")


def format_report(diff: Diff, results: JUnitResults, baseline: frozenset[str]) -> str:
    """Render a human-readable Markdown report of the diff."""
    lines: list[str] = []
    lines.append("## Nightly reference-aware xfail diff")
    lines.append("")
    lines.append(
        f"baseline: {len(baseline)} known-failing · "
        f"run: {len(results.failed)} failing, "
        f"{len(results.passed)} passing, {len(results.skipped)} skipped"
    )
    lines.append("")
    if not diff.is_drift:
        lines.append("Exact match — the FAIL set equals the baseline. No drift.")
        return "\n".join(lines) + "\n"
    if diff.new_failures:
        lines.append(f"### NEW FAILURE ({len(diff.new_failures)})")
        lines.append("A test is failing that the baseline does not list — a regression or drift.")
        lines.extend(f"- `{tid}`" for tid in sorted(diff.new_failures))
        lines.append("")
    if diff.now_passing:
        lines.append(f"### UNEXPECTED PASS ({len(diff.now_passing)})")
        lines.append("A baselined test now passes — update the baseline so it cannot rot.")
        lines.extend(f"- `{tid}`" for tid in sorted(diff.now_passing))
        lines.append("")
    if diff.now_skipped:
        lines.append(f"### BASELINED BUT NOW SKIPPED ({len(diff.now_skipped)})")
        lines.append(
            "A baselined test is skipped (e.g. newly `#[ignore]`d), so it no longer "
            "fails — update the baseline."
        )
        lines.extend(f"- `{tid}`" for tid in sorted(diff.now_skipped))
        lines.append("")
    if diff.absent:
        lines.append(f"### BASELINED BUT NOT OBSERVED ({len(diff.absent)})")
        lines.append(
            "A baselined test was not observed at all (removed or renamed) — update the baseline."
        )
        lines.extend(f"- `{tid}`" for tid in sorted(diff.absent))
        lines.append("")
    lines.append(
        "Regenerate the baseline after a deliberate change with:\n"
        "```\n"
        "python scripts/diff_nightly_xfail_baseline.py \\\n"
        "  --junit target/nextest/nightly-xfail/junit.xml \\\n"
        "  --baseline tests/fixtures/nightly-reference-aware/xfail-baseline.txt \\\n"
        "  --update-baseline\n"
        "```"
    )
    return "\n".join(lines) + "\n"


def run(junit_path: Path, baseline_path: Path, update: bool) -> int:
    """Execute the diff (or the baseline update). Returns the process exit code."""
    junit_text = junit_path.read_text(encoding="utf-8")
    results = parse_junit(junit_text)

    if update:
        rendered = render_baseline(results.failed)
        # Round-trip before writing: whatever is written must be readable by the
        # strict parser above, so `--update-baseline` cannot mint a file that the
        # next run refuses (or, worse, silently shrinks).
        if read_baseline(rendered) != results.failed:
            raise ValueError(
                "the rendered baseline does not read back as the FAIL set it was rendered from"
            )
        baseline_path.parent.mkdir(parents=True, exist_ok=True)
        baseline_path.write_text(rendered, encoding="utf-8")
        print(f"wrote {len(results.failed)} failing test ids to {baseline_path}")
        return 0

    baseline = read_baseline(baseline_path.read_text(encoding="utf-8"))
    diff = diff_baseline(results, baseline)
    print(format_report(diff, results, baseline))

    if diff.new_failures:
        _annotate(
            "error",
            f"{len(diff.new_failures)} new reference-aware failure(s) not in the xfail baseline",
        )
    if diff.unexpected_passes:
        _annotate(
            "error",
            f"{len(diff.unexpected_passes)} baselined test(s) no longer failing — "
            "update the xfail baseline",
        )
    return 1 if diff.is_drift else 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--junit",
        type=Path,
        required=True,
        help="path to the nextest JUnit report (target/nextest/nightly-xfail/junit.xml)",
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        required=True,
        help="path to the committed xfail baseline file",
    )
    parser.add_argument(
        "--update-baseline",
        action="store_true",
        help="rewrite the baseline from the JUnit report's FAIL set instead of diffing",
    )
    args = parser.parse_args(argv)
    try:
        return run(args.junit, args.baseline, args.update_baseline)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
