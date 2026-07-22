#!/usr/bin/env python3
"""Refresh tests/fixtures/mutalyzer-normalize/cases.json from upstream mutalyzer.

Pins the upstream commit at MUTALYZER_SHA. Run by hand when we choose to
sync with upstream:

    pixi run python scripts/refresh-mutalyzer-fixtures.py refresh

The script reads tests/variants_set.py from the pinned mutalyzer checkout
WITHOUT modifying that checkout's working tree (uses `git show <sha>:<path>`).
It then exec()s that file in an isolated namespace to extract TESTS_ALL and
serializes the cases to cases.json. Python tuples become JSON arrays; None
becomes null.

The script also exposes a one-shot `annotate-xfails` subcommand used during
initial import to record which axes are expected to fail for which inputs.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

# --- Pin to a specific upstream commit ---------------------------------------
MUTALYZER_SHA = "00045d63feba04ccff18a4d95f4245ad5d8d3092"  # "Bump version to 3.1.1"
MUTALYZER_REPO = Path.home() / "work" / "git" / "mutalyzer"
VARIANTS_SET_REL = "tests/variants_set.py"

# --- Output paths -------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "mutalyzer-normalize"
CASES_JSON = FIXTURE_DIR / "cases.json"
NOTICE = FIXTURE_DIR / "NOTICE"


def read_variants_set_at_sha() -> str:
    """Read variants_set.py from the pinned commit without mutating the
    upstream checkout. Returns the file's source as a string."""
    if not MUTALYZER_REPO.exists():
        sys.exit(f"error: mutalyzer repo not found at {MUTALYZER_REPO}")
    out = subprocess.check_output(
        ["git", "-C", str(MUTALYZER_REPO), "show", f"{MUTALYZER_SHA}:{VARIANTS_SET_REL}"],
        text=True,
    )
    return out


def load_tests_all() -> list[dict]:
    """Exec variants_set.py source in an isolated namespace and return TESTS_ALL."""
    source = read_variants_set_at_sha()
    ns: dict = {}
    exec(compile(source, f"{MUTALYZER_SHA}:{VARIANTS_SET_REL}", "exec"), ns)
    tests_all = ns.get("TESTS_ALL")
    if tests_all is None:
        sys.exit("error: TESTS_ALL not defined in variants_set.py at the pinned SHA")
    return list(tests_all)


def to_jsonable(obj):
    """Recursively convert Python tuples/sets to lists for JSON serialization."""
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(x) for x in obj]
    if isinstance(obj, (set, frozenset)):
        # Preserve order via sorted str repr for determinism across runs.
        return sorted((to_jsonable(x) for x in obj), key=lambda x: json.dumps(x, sort_keys=True))
    if isinstance(obj, dict):
        return {k: to_jsonable(v) for k, v in obj.items()}
    return obj


def normalize_case(case: dict) -> dict:
    """Project an upstream case dict into the cases.json shape.

    - Converts tuples/sets-of-tuples (e.g. coding_protein_descriptions) into
      lists of lists.
    - Ensures every case has an 'xfail' key (empty list by default).
    """
    out = to_jsonable(case)
    out.setdefault("xfail", [])
    return out


HAND_ADDED_PREFIX = "ferro-hand-derived:"


# Per-case curation keys that record *why* a row behaves as it does. These live
# on UPSTREAM rows (they are our disposition of an upstream case), so a refresh
# that replaced upstream rows wholesale would silently destroy them.
DISPOSITION_KEYS = (
    "spec_citation",
    "accepted_divergence",
    "known_bug",
    "improvement",
    "reference_unavailable",
    "accepted_rejection",
)

# Expected-value fields that we sometimes hand-correct away from what upstream
# asserts (upstream is not always spec-correct). A corrected value must survive
# a refresh, or the affected rows flip PASS -> FAIL.
EXPECTED_VALUE_KEYS = (
    "normalized",
    "genomic",
    "protein_description",
    "coding_protein_descriptions",
    "rna_description",
    "noncoding",
    "errors",
    "infos",
)


def _load_existing(path: Path) -> dict:
    """Return the current fixture dict, or {} when absent."""
    if not path.exists():
        return {}
    try:
        with path.open() as fh:
            return json.load(fh)
    except json.JSONDecodeError as exc:
        raise SystemExit(
            f"error: {path} is not valid JSON ({exc}). Fix or restore the file "
            "(e.g. `git checkout -- <path>`) before refreshing. Do NOT rerun with "
            "--force: an unreadable fixture cannot be merged, and --force would "
            "overwrite all hand curation with bare upstream rows."
        ) from exc


def _is_hand_added(case: dict) -> bool:
    """A case is hand-added when ANY of its keywords carries the marker prefix.

    Checking every keyword (not just the first) keeps a hand row identifiable
    if a future edit prepends a keyword — otherwise it would be reclassified as
    upstream and deleted by the next refresh.
    """
    kws = case.get("keywords") or []
    return any(str(k).startswith(HAND_ADDED_PREFIX) for k in kws)


def _merge_case(existing_case: dict, upstream_case: dict) -> tuple[dict, bool, bool, bool]:
    """Overlay one existing row's curation onto its fresh upstream row.

    Returns `(merged_row, carried_disposition, carried_correction,
    carried_removal)`. `carried_correction` flags an expected-value field whose
    curated value overrides upstream's; `carried_removal` flags a row carrying a
    non-empty `removed_keys` marker (a deliberately-dropped expected value) —
    tracked separately because a removal is curation the refresh must reconcile
    even when upstream no longer ships the removed field, in which case there is
    no value difference to flag as a correction. The
    upstream row is the base so genuinely-new upstream fields still arrive;
    onto it we overlay

    - every disposition key present on the existing row, and
    - every expected-value field whose existing value differs from upstream's;
    - every expected-value field the curated row *deliberately removed*, recorded
      explicitly in its `removed_keys` list. Only a key named there is dropped —
      an expected-value field absent from the curated row but NOT named in
      `removed_keys` is treated as genuinely new upstream data and kept, so a
      field upstream added after the row was curated survives the merge. The
      `removed_keys` marker itself is carried forward so the removal keeps
      holding across future refreshes;
    - any keyword the curated row added and upstream does not have. Keyword
      annotations are how a disposition's rationale is recorded, and appending
      (rather than replacing) keeps new upstream keywords as well.
    """
    merged = dict(upstream_case)

    dispositions = {k: existing_case[k] for k in DISPOSITION_KEYS if k in existing_case}
    merged.update(dispositions)

    removed_keys = existing_case.get("removed_keys") or []
    unknown_removed = [k for k in removed_keys if k not in EXPECTED_VALUE_KEYS]
    if unknown_removed:
        raise ValueError(
            f"removed_keys for input {existing_case.get('input')!r} names non-expected-value "
            f"key(s) {unknown_removed}; only {list(EXPECTED_VALUE_KEYS)} may be removed"
        )
    corrected = False
    for key in EXPECTED_VALUE_KEYS:
        in_existing, in_upstream = key in existing_case, key in upstream_case
        if in_existing and existing_case[key] != upstream_case.get(key):
            merged[key] = existing_case[key]
            corrected = True
        elif not in_existing and in_upstream and key in removed_keys:
            del merged[key]
    if removed_keys:
        merged["removed_keys"] = list(removed_keys)

    upstream_keywords = list(upstream_case.get("keywords") or [])
    added_keywords = [
        k for k in (existing_case.get("keywords") or []) if k not in upstream_keywords
    ]
    if added_keywords:
        merged["keywords"] = upstream_keywords + added_keywords

    return merged, bool(dispositions), corrected, bool(removed_keys)


def _pair_with_existing(
    upstream_cases: list[dict], existing_cases: list[dict]
) -> tuple[list[tuple[dict, dict | None]], list[dict]]:
    """Pair each upstream row with the curated row it supersedes.

    Rows are matched on `input`, and on *occurrence order* within a duplicated
    `input` — upstream occasionally ships the same variant twice with different
    expected values, and collapsing those onto one curated row would copy the
    wrong curation onto both.

    Returns `(pairs, unmatched)`. `unmatched` holds the curated (non-hand-added)
    rows that no upstream occurrence claimed — e.g. upstream dropped one of a
    duplicated `input` that carried curation. Those rows must not be dropped
    silently (that is the very curation loss the merge exists to prevent), so
    the caller carries them forward and counts them.
    """
    by_input: dict[str | None, list[dict]] = {}
    for case in existing_cases:
        if not _is_hand_added(case):
            by_input.setdefault(case.get("input"), []).append(case)
    seen: dict[str | None, int] = {}
    pairs: list[tuple[dict, dict | None]] = []
    for upstream in upstream_cases:
        key = upstream.get("input")
        ordinal = seen.get(key, 0)
        seen[key] = ordinal + 1
        candidates = by_input.get(key, [])
        pairs.append((upstream, candidates[ordinal] if ordinal < len(candidates) else None))
    unmatched: list[dict] = []
    for key, candidates in by_input.items():
        claimed = seen.get(key, 0)
        unmatched.extend(candidates[claimed:])
    return pairs, unmatched


def build_refresh_payload(upstream_cases: list[dict], existing: dict, force: bool) -> dict:
    """Merge upstream cases with preserved hand curation.

    Each fresh upstream row is taken as the base and matched to the curated row
    it supersedes (see `_pair_with_existing`); that row's per-case dispositions,
    hand-corrected expected values and added keyword annotations are overlaid
    onto it (see `_merge_case`), so genuinely-new upstream fields still arrive
    while curation survives. Curated rows that no upstream occurrence claims
    (e.g. upstream dropped one of a duplicated `input`) are carried forward
    rather than silently dropped, so their curation is never lost. Hand-added
    rows (marked via HAND_ADDED_PREFIX) and corpus-level curation (clusters,
    comparator_provenance) are carried forward verbatim. Raises when the existing
    file holds curation and force is not set, so the operator can reconcile the
    merge first.
    """
    existing_cases = existing.get("cases", [])
    preserved = [c for c in existing_cases if _is_hand_added(c)]
    pairs, unmatched = _pair_with_existing(upstream_cases, existing_cases)

    merged_cases: list[dict] = []
    n_disp = n_corrected = n_removed = 0
    for upstream, prior in pairs:
        if prior is None:
            merged_cases.append(dict(upstream))
            continue
        merged, carried_disposition, carried_correction, carried_removal = _merge_case(
            prior, upstream
        )
        n_disp += carried_disposition
        n_corrected += carried_correction
        n_removed += carried_removal
        merged_cases.append(merged)
    n_hand = len(preserved)
    n_unmatched = len(unmatched)

    has_curation = bool(
        existing.get("clusters")
        or existing.get("comparator_provenance")
        or preserved
        or n_disp
        or n_corrected
        or n_removed
        or n_unmatched
    )
    if has_curation and not force:
        raise SystemExit(
            "refusing to refresh: cases.json carries hand curation that a refresh "
            "must reconcile. At stake:\n"
            f"  - {n_hand} hand-added case(s) (keyword prefix {HAND_ADDED_PREFIX!r})\n"
            f"  - {n_disp} upstream row(s) carrying per-case dispositions "
            f"({', '.join(DISPOSITION_KEYS)})\n"
            f"  - {n_corrected} upstream row(s) whose expected values were hand-corrected "
            "away from upstream (reverting these flips them PASS -> FAIL)\n"
            f"  - {n_removed} upstream row(s) with deliberately-removed expected values "
            "(removed_keys; resurrecting the removed key flips them PASS -> FAIL)\n"
            f"  - {n_unmatched} curated row(s) upstream no longer ships (paired to no "
            "upstream occurrence; carried forward so their curation is not lost)\n"
            f"  - {len(existing.get('clusters', []))} cluster(s), "
            f"comparator_provenance {'present' if existing.get('comparator_provenance') else 'absent'}\n"
            "Rerun with --force once you have reviewed the merge; --force keeps the "
            "curation above but accepts every other upstream change."
        )

    payload = {
        "description": "Mutalyzer normalizer test cases",
        "source": "https://github.com/mutalyzer/mutalyzer",
        "source_commit": MUTALYZER_SHA,
        "license": "MIT",
        "refreshed_at": date.today().isoformat(),
    }
    if existing.get("comparator_provenance"):
        payload["comparator_provenance"] = existing["comparator_provenance"]
    if existing.get("clusters"):
        payload["clusters"] = existing["clusters"]
    payload["cases"] = merged_cases + unmatched + preserved
    print(
        f"merge: preserved dispositions on {n_disp} row(s), "
        f"corrected expected values on {n_corrected} row(s), "
        f"{n_hand} hand-added row(s), "
        f"{n_unmatched} unmatched curated row(s) carried forward"
    )
    return payload


def cmd_refresh(force: bool) -> None:
    tests_all = load_tests_all()
    cases = [normalize_case(c) for c in tests_all]

    existing = _load_existing(CASES_JSON)
    payload = build_refresh_payload(cases, existing, force=force)

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    with CASES_JSON.open("w") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")

    NOTICE.write_text(
        f"""Mutalyzer normalizer test fixtures
==================================

Source:  https://github.com/mutalyzer/mutalyzer
Commit:  {MUTALYZER_SHA}
License: MIT (Leiden University Medical Center)

cases.json is generated by scripts/refresh-mutalyzer-fixtures.py from
tests/variants_set.py in the upstream repo at the pinned commit above,
AND hand-curated thereafter. Hand curation now includes: the clusters
taxonomy, comparator_provenance, all per-case dispositions, and
hand-derived cases whose first keyword begins "ferro-hand-derived:".

Refresh MERGES rather than replaces. build_refresh_payload() takes each
fresh upstream row as the base, matches it to the existing row with the
same `input`, and overlays that row's per-case dispositions
(spec_citation, accepted_divergence, known_bug, improvement,
reference_unavailable, accepted_rejection) plus any expected-value field
we had hand-corrected away from upstream (normalized, genomic,
protein_description, coding_protein_descriptions, rna_description,
noncoding, errors, infos). A field we deliberately REMOVED is recorded
explicitly in the curated row's "removed_keys" list; only a key named
there is dropped on refresh. An expected-value field absent from the
curated row but NOT in "removed_keys" is treated as genuinely-new
upstream data and kept — so a field upstream adds after a row is curated
still arrives. Keyword annotations the curated row added are appended to
upstream's keywords. New upstream fields therefore still arrive. Hand-added
rows and corpus-level curation are carried forward verbatim. Curated rows
that no upstream occurrence claims (e.g. upstream dropped one of a
duplicated `input`) are carried forward too, never silently dropped. Rows
are paired on `input` and on occurrence order within a duplicated `input`.

The script still REFUSES to run while curation is present, listing what
is at stake, so the merge is reviewed rather than assumed. --force is the
override for when you have reconciled the merge — NOT the routine path.

Expected values for hand-derived cases are DERIVED FROM THE PINNED
COMPARATOR (mutalyzer 3.1.1 / mutalyzer-hgvs-parser 0.4.0), never copied
from upstream or invented — see comparator_provenance.note.

To refresh against a newer upstream commit:
    1. Update MUTALYZER_SHA in scripts/refresh-mutalyzer-fixtures.py
    2. pixi run python scripts/refresh-mutalyzer-fixtures.py refresh
       It will refuse and enumerate the curation at stake. Review that
       list; once you are satisfied the merge is correct, rerun with
       --force to apply it.
    3. Re-run the mutalyzer-normalize tests; annotate any new xfails via:
         pixi run python scripts/refresh-mutalyzer-fixtures.py annotate-xfails \\
           --axis <axis> --from /tmp/ferro-xfail/<axis>.txt
    4. Commit cases.json + NOTICE in a single commit.

The upstream MIT license requires attribution; this file satisfies that
requirement. See the upstream LICENSE for the full text.
"""
    )

    n_to_test = sum(1 for c in payload["cases"] if c.get("to_test"))
    print(
        f"wrote {len(payload['cases'])} cases ({n_to_test} to_test) "
        f"to {CASES_JSON.relative_to(REPO_ROOT)}"
    )


def cmd_annotate_xfails(axis: str, failing_inputs_path: Path) -> None:
    """One-shot: add `axis` to xfail[] for every case whose `input` is listed
    in failing_inputs_path (one input per line).

    When upstream contains multiple cases with the same `input` (it occasionally
    does — same variant, different expected `normalized`/`genomic` etc.), we
    only mark the case if its actual axis output diverges from the expected
    on inspection. The conservative thing is to mark every duplicate;
    over-marked cases will surface as XPASS on the next run and can be
    manually demoted.
    """
    failing = {
        line.strip() for line in failing_inputs_path.read_text().splitlines() if line.strip()
    }
    if not failing:
        print(f"no inputs to annotate (file {failing_inputs_path} was empty)")
        return
    data = json.loads(CASES_JSON.read_text())
    n_updated = 0
    for case in data["cases"]:
        if case.get("input") in failing:
            xfails = case.setdefault("xfail", [])
            if axis not in xfails:
                xfails.append(axis)
                n_updated += 1
    with CASES_JSON.open("w") as fh:
        json.dump(data, fh, indent=2, sort_keys=False)
        fh.write("\n")
    print(f"annotated {n_updated} cases with xfail={axis!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_refresh = sub.add_parser("refresh", help="regenerate cases.json from pinned upstream SHA")
    p_refresh.add_argument(
        "--force",
        action="store_true",
        help="proceed even though cases.json carries hand curation",
    )
    ann = sub.add_parser("annotate-xfails", help="record xfailed inputs for an axis")
    ann.add_argument("--axis", required=True)
    ann.add_argument("--from", dest="from_file", type=Path, required=True)
    args = parser.parse_args()
    if args.cmd == "refresh":
        cmd_refresh(force=args.force)
    elif args.cmd == "annotate-xfails":
        cmd_annotate_xfails(args.axis, args.from_file)


if __name__ == "__main__":
    main()
