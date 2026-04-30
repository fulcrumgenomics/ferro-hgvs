"""Build the HGVS spec normalization regression fixture.

Reads the curated extract at
``tests/fixtures/grammar/hgvs_spec_examples_v21.md``, runs every variant
string through ``ferro_hgvs.normalize()``, and writes a structured JSON
fixture at ``tests/fixtures/grammar/hgvs_spec_normalization.json``.

The output captures, for each spec example:

- ``input`` — the raw token from the spec
- ``input_prefixed`` — input with a default reference prefix prepended
  (only when needed, e.g. ``c.78T>C`` → ``NM_004006.2:c.78T>C``)
- ``spec_expected`` — the spec-correct form
- ``current`` — what ferro emits today (``null`` for groups skipped at
  runtime)
- ``status`` — short tag (e.g. ``preserved``, ``transformed``,
  ``diverges``, ``rejected``, ``false-acceptance``)
- ``todo`` — optional pointer at the master tracking issue when current
  output diverges from spec

Companion to the regression test in
``tests/hgvs_spec_normalization_tests.rs`` (issue #84) and the master
spec-compliance tracker (issue #83).

Re-run after any normalization change to refresh the ``current`` field
of affected rows; the fixture is byte-deterministic for unchanged
behavior.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import ferro_hgvs

ROOT = Path(__file__).resolve().parents[1]
SOURCE_MD = ROOT / "tests" / "fixtures" / "grammar" / "hgvs_spec_examples_v21.md"
FIXTURE = ROOT / "tests" / "fixtures" / "grammar" / "hgvs_spec_normalization.json"

SPEC_SOURCE = "https://hgvs-nomenclature.org/stable/"
SPEC_REPO = "https://github.com/HGVSnomenclature/hgvs-nomenclature"
SPEC_VERSION = "21.0"
TRACKING_ISSUE = "https://github.com/fulcrumgenomics/ferro-hgvs/issues/83"
TODO_TAG = f"see {TRACKING_ISSUE}"

# Default reference prefixes per coordinate system, used to upgrade bare
# illustrative fragments (e.g. "c.78T>C") into parseable HGVS strings.
# Picked to match the most-referenced accessions in the spec doc.
DEFAULT_PREFIX = {
    "c": "NM_004006.2",
    "n": "NR_002196.1",
    "r": "NM_004006.3",
    "g": "NC_000023.11",
    "p": "NP_003997.1",
    "m": "NC_012920.1",
    "o": "NC_000023.11",
}

# Sections that need reference data ferro can't access via the JSON
# loader today — see issue #83 §B for context. Cases here are emitted
# with current=null and skipped by the test until #82 lands.
NEEDS_REFERENCE_SECTIONS = {"2.1 3' rule shifting (true normalization)"}

TICK = re.compile(r"`([^`]+)`")
SUBHEAD = re.compile(r"^###\s+(\S+)\s+(.+)$")
TOPHEAD = re.compile(r"^##\s+(\d+)\.\s+(.+)$")
BARE_RE = re.compile(r"^\[?([cnrgpmo])\.")
REF_RE = re.compile(r"^[A-Z]+_?\d+(\.\d+)?(t\d+|p\d+)?[:\(]")


@dataclass
class Case:
    section: str
    kind: str  # compliant | pair | reject
    input: str
    spec_expected: str | None  # None → "<same as input>" for compliant; "<parse-error>" for reject
    input_prefixed: str | None = None
    spec_expected_prefixed: str | None = None
    current: str | None = None
    status: str = ""
    todo: str | None = None


def is_variant_like(s: str) -> bool:
    """Filter out backticked tokens that are clearly not HGVS strings."""
    if len(s) < 3:
        return False
    if re.search(r"[gcrnpmo]\.", s) or s.startswith("[") or "ext" in s:
        return True
    if re.match(r"^[A-Z]+_?\d+(\.\d+)?$", s):
        return True
    if "ins" in s or "del" in s or "dup" in s or "inv" in s:
        return True
    return False


def add_prefix(s: str) -> str | None:
    """Return s with a default reference prefix prepended, or None if none needed."""
    if REF_RE.match(s):
        return None
    m = BARE_RE.match(s)
    if not m:
        return None
    prefix = DEFAULT_PREFIX.get(m.group(1))
    if not prefix:
        return None
    if s.startswith("["):
        return f"[{prefix}:{s[1:]}"
    return f"{prefix}:{s}"


def parse_spec(text: str) -> list[Case]:
    """Walk the curated markdown and emit one Case per backticked variant."""
    cases: list[Case] = []
    section = "(top)"
    in_section = "1"

    for line in text.splitlines():
        if (m := TOPHEAD.match(line)) is not None:
            in_section = m.group(1)
            section = f"{m.group(1)}. {m.group(2)}"
            continue
        if (m := SUBHEAD.match(line)) is not None:
            section = f"{m.group(1)} {m.group(2)}"
            continue
        if not line.lstrip().startswith("- "):
            continue
        body = line.lstrip()[2:]
        ticks = TICK.findall(body)
        if not ticks:
            continue

        if in_section == "1":
            for t in ticks:
                if not is_variant_like(t):
                    continue
                cases.append(Case(section=section, kind="compliant", input=t, spec_expected=None))
        elif in_section == "2":
            if "→" not in body:
                continue
            left, _, right = body.partition("→")
            lt = TICK.findall(left)
            rt = TICK.findall(right)
            if not lt or not rt:
                continue
            inp, exp = lt[0], rt[0]
            if not is_variant_like(inp) or not is_variant_like(exp):
                continue
            cases.append(Case(section=section, kind="pair", input=inp, spec_expected=exp))
        elif in_section == "3":
            for t in ticks:
                if not is_variant_like(t):
                    continue
                cases.append(Case(section=section, kind="reject", input=t, spec_expected=None))
    return cases


def evaluate(cases: list[Case]) -> None:
    """Apply default-prefix logic and run each case through normalize()."""
    for c in cases:
        c.input_prefixed = add_prefix(c.input)
        if c.spec_expected is not None:
            c.spec_expected_prefixed = add_prefix(c.spec_expected)

        if c.section in NEEDS_REFERENCE_SECTIONS:
            c.status = "needs-reference"
            c.todo = TODO_TAG
            continue

        # Reject cases are tested without prefixing — the parser is meant
        # to refuse them as written.
        target = c.input if c.kind == "reject" else (c.input_prefixed or c.input)

        try:
            c.current = ferro_hgvs.normalize(target)
        except Exception:
            c.current = None

        c.status = classify(c)
        if c.status not in ("preserved", "transformed", "rejected"):
            c.todo = TODO_TAG


def classify(c: Case) -> str:
    eff_in = c.input_prefixed or c.input
    eff_exp = c.spec_expected_prefixed or c.spec_expected

    if c.kind == "compliant":
        if c.current is None:
            return "parse-error"
        return "preserved" if c.current == eff_in else "diverges"
    if c.kind == "pair":
        if c.current is None:
            return "parse-error"
        if c.current == eff_exp:
            return "transformed"
        if c.current == eff_in:
            return "unchanged"
        return "diverges"
    if c.kind == "reject":
        return "rejected" if c.current is None else "false-acceptance"
    return "?"


def to_groups(cases: list[Case]) -> list[dict[str, Any]]:
    """Bucket cases into ordered groups, one per spec section."""
    by_section: dict[str, list[Case]] = {}
    order: list[str] = []
    for c in cases:
        if c.section not in by_section:
            by_section[c.section] = []
            order.append(c.section)
        by_section[c.section].append(c)

    groups: list[dict[str, Any]] = []
    for sec in order:
        section_cases = by_section[sec]
        kind = (
            "needs-reference"
            if sec in NEEDS_REFERENCE_SECTIONS
            else section_cases[0].kind
        )
        rows: list[dict[str, Any]] = []
        for c in section_cases:
            row: dict[str, Any] = {
                "input": c.input,
                "spec_expected": expected_for_output(c),
                "current": c.current,
                "status": c.status,
            }
            if c.input_prefixed:
                row["input_prefixed"] = c.input_prefixed
            if c.spec_expected_prefixed:
                row["spec_expected_prefixed"] = c.spec_expected_prefixed
            if c.todo:
                row["todo"] = c.todo
            rows.append(row)
        groups.append({"name": sec, "kind": kind, "cases": rows})
    return groups


def expected_for_output(c: Case) -> str:
    if c.kind == "compliant":
        return c.spec_expected or "<same as input>"
    if c.kind == "reject":
        return "<parse-error>"
    return c.spec_expected or ""


def main() -> None:
    text = SOURCE_MD.read_text()
    cases = parse_spec(text)
    evaluate(cases)

    fixture = {
        "description": (
            "HGVS spec v21.0 examples paired with ferro's current "
            "normalize() output. Used by tests/hgvs_spec_normalization_tests.rs "
            "(issue #84) and tracked under issue #83."
        ),
        "spec_source": SPEC_SOURCE,
        "spec_repo": SPEC_REPO,
        "spec_version": SPEC_VERSION,
        "tracking_issue": TRACKING_ISSUE,
        "ferro_version": ferro_hgvs.__version__,
        "default_prefixes": DEFAULT_PREFIX,
        "groups": to_groups(cases),
    }

    FIXTURE.write_text(json.dumps(fixture, indent=2, sort_keys=False) + "\n")

    summary: dict[tuple[str, str], int] = {}
    for c in cases:
        summary.setdefault((c.kind, c.status), 0)
        summary[(c.kind, c.status)] += 1
    print(f"Wrote {FIXTURE.relative_to(ROOT)} ({len(cases)} cases)")
    for (kind, status), n in sorted(summary.items()):
        print(f"  {kind:10s} {status:20s} {n:>4}")


if __name__ == "__main__":
    main()
