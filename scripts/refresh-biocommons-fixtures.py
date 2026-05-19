#!/usr/bin/env python3
"""Refresh tests/fixtures/biocommons-normalize/cases.json.

The biocommons/hgvs normalizer corpus comes from two upstream files:

  1. `tests/test_hgvs_normalizer.py` in biocommons/hgvs — primary corpus
     across `test_c_normalizer`, `test_g_normalizer`, `test_norm_var_near_bound`.
     This file was transcribed verbatim into Rust in the hgvs-rs port at
     `src/normalizer.rs::mod test`. We extract from the Rust port because
     its case-tuple format is easier to regex than the Python AST and the
     port is feature-complete against the Python tests (one BRCA2 case is
     re-classified in hgvs-rs — we restore biocommons' assertRaises
     expectation here).
  2. `tests/issues/test_02xx.py` + `tests/issues/test_03xx.py` — six
     normalize-bearing issue-regression cases that the hgvs-rs port did
     not transcribe. We include them by hand-curated entries below.

The hgvs-rs source is read via `git show <sha>:<path>` so this script does
not mutate the upstream checkout's working tree.

Usage:
    pixi run python scripts/refresh-biocommons-fixtures.py refresh
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

# --- Pin to a specific upstream commit ---------------------------------------
HGVS_RS_SHA = "b6513b3d793615fd770b1e8ae2a87d85d10396b5"  # v0.20.2-1-gb6513b3
HGVS_RS_REPO = Path.home() / "work" / "git" / "hgvs-rs"
NORMALIZER_RS_REL = "src/normalizer.rs"

# Biocommons project doesn't supply a tagged version we pin against — the
# Rust port is the proximate upstream. Recorded for attribution only.
BIOCOMMONS_REPO = "https://github.com/biocommons/hgvs"

# --- Output paths ------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "biocommons-normalize"
CASES_JSON = FIXTURE_DIR / "cases.json"
NOTICE = FIXTURE_DIR / "NOTICE"


# Map hgvs-rs Normalizer variable names to (shuffle_direction, cross_boundaries).
# Derived from the `normalizers()` helper in hgvs-rs src/normalizer.rs at the
# pinned commit:
#   norm    -> Direction::FiveToThree, cross=true   (== ThreePrime + cross=true)
#   norm5   -> Direction::ThreeToFive, cross=true   (== FivePrime  + cross=true)
#   normc   -> Direction::FiveToThree, cross=false  (== ThreePrime + cross=false)
#   norm5c  -> Direction::ThreeToFive, cross=false  (== FivePrime  + cross=false)
VAR_CONFIG = {
    "norm": ("3prime", True),
    "norm5": ("5prime", True),
    "normc": ("3prime", False),
    "norm5c": ("5prime", False),
}
# `cases3` is iterated with `norm`; `cases5` is iterated with `norm5`.
VEC_VAR = {"cases3": "norm", "cases5": "norm5"}


def read_normalizer_rs_at_sha() -> str:
    if not HGVS_RS_REPO.exists():
        sys.exit(f"error: hgvs-rs repo not found at {HGVS_RS_REPO}")
    return subprocess.check_output(
        ["git", "-C", str(HGVS_RS_REPO), "show", f"{HGVS_RS_SHA}:{NORMALIZER_RS_REL}"],
        text=True,
    )


def find_test_functions(text: str) -> list[tuple[str, int]]:
    pattern = re.compile(r"#\[test\]\s*\n\s*fn\s+(normalize_[a-z0-9_]+)\s*\(\s*\)", re.MULTILINE)
    return [(m.group(1), m.end()) for m in pattern.finditer(text)]


def find_fn_body(text: str, after: int) -> str:
    """Return the contents of the {...} block whose `{` appears after `after`."""
    i = text.index("{", after)
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1 : j]
    raise RuntimeError("unbalanced braces")


def split_braced_blocks(body: str):
    """Yield each top-level `{...}` block (depth-1) within `body`."""
    depth = 0
    start = -1
    for i, c in enumerate(body):
        if c == "{":
            if depth == 0:
                start = i
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                yield body[start + 1 : i]


def extract_cases_from_hgvs_rs(text: str) -> list[dict]:
    cases: list[dict] = []
    seen: set[tuple] = set()  # dedupe key: (fn, input, normalized, dir, cross, expects_err)

    fns = find_test_functions(text)
    for fn_name, end in fns:
        body = find_fn_body(text, end)

        # 1. `let casesN = vec![ ("input", "expected"), ... ];` arrays
        for vm in re.finditer(r"let\s+(cases\d)\s*=\s*vec!\[(.+?)\];", body, flags=re.DOTALL):
            var_name = VEC_VAR.get(vm.group(1))
            if not var_name:
                continue
            direction, cross = VAR_CONFIG[var_name]
            for tm in re.finditer(r'\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,?\s*\)', vm.group(2)):
                key = (fn_name, tm.group(1), tm.group(2), direction, cross, False)
                if key in seen:
                    continue
                seen.add(key)
                cases.append(
                    {
                        "input": tm.group(1),
                        "normalized": tm.group(2),
                        "shuffle_direction": direction,
                        "cross_boundaries": cross,
                        "source_function": fn_name,
                        "to_test": True,
                    }
                )

        # 2. Inline `{ ... }` blocks: HgvsVariant::from_str + (is_err | exp + normalize)
        for blk in split_braced_blocks(body):
            im = re.search(r'HgvsVariant::from_str\(\s*"([^"]+)"', blk) or re.search(
                r'HgvsVariant::from_str\(\s*\n\s*"([^"]+)"', blk
            )
            if not im:
                continue
            input_str = im.group(1)

            # 2a. is_err form: assert!(<var>.normalize(&raw).is_err())
            errm = re.search(r"(norm5c|normc|norm5|norm)\.normalize\(\s*&raw\s*\)\.is_err\(\)", blk)
            if errm:
                var = errm.group(1)
                direction, cross = VAR_CONFIG[var]
                key = (fn_name, input_str, None, direction, cross, True)
                if key in seen:
                    continue
                seen.add(key)
                cases.append(
                    {
                        "input": input_str,
                        "normalized": None,
                        "expects_error": True,
                        "shuffle_direction": direction,
                        "cross_boundaries": cross,
                        "source_function": fn_name,
                        "to_test": True,
                    }
                )
                continue

            # 2b. let exp + let res = <var>.normalize(&raw)
            expm = re.search(r'let\s+exp\s*=\s*"([^"]+)"', blk)
            varm = re.search(
                r"let\s+res\s*=\s*(norm5c|normc|norm5|norm)\.normalize\s*\(\s*&raw\s*\)", blk
            )
            if expm and varm:
                var = varm.group(1)
                direction, cross = VAR_CONFIG[var]
                key = (fn_name, input_str, expm.group(1), direction, cross, False)
                if key in seen:
                    continue
                seen.add(key)
                cases.append(
                    {
                        "input": input_str,
                        "normalized": expm.group(1),
                        "shuffle_direction": direction,
                        "cross_boundaries": cross,
                        "source_function": fn_name,
                        "to_test": True,
                    }
                )
                continue

            # 2c. assert_eq!("expected", format!("{}", <var>.normalize(&raw)?))
            altm = re.search(
                r'assert_eq!\(\s*"([^"]+)"\s*,\s*format!\s*\(\s*"\{\}"\s*,\s*(norm5c|normc|norm5|norm)\.normalize',
                blk,
            )
            if altm:
                var = altm.group(2)
                direction, cross = VAR_CONFIG[var]
                key = (fn_name, input_str, altm.group(1), direction, cross, False)
                if key in seen:
                    continue
                seen.add(key)
                cases.append(
                    {
                        "input": input_str,
                        "normalized": altm.group(1),
                        "shuffle_direction": direction,
                        "cross_boundaries": cross,
                        "source_function": fn_name,
                        "to_test": True,
                    }
                )

    return cases


# Cases that biocommons python's `tests/issues/test_02xx.py` and `test_03xx.py`
# carry but the hgvs-rs port did not transcribe. All use the test-file's
# `self.vn = Normalizer(hdp, shuffle_direction=3, cross_boundaries=True)`
# config, except issue_376 which uses `self.hn = Normalizer(hdp)` with
# biocommons' defaults (`cross_boundaries=False` per src/hgvs/_data/defaults.ini).
ISSUE_CASES: list[dict] = [
    {
        "input": "NG_029146.1:g.6494delG",
        "normalized": None,
        "expects_error": True,
        "shuffle_direction": "3prime",
        "cross_boundaries": True,
        "source_function": "issue_293_parser_attribute_assignment_error",
        "to_test": True,
    },
    {
        "input": "NM_000535.5:c.1673_1674inv",
        "normalized": "NM_000535.5:c.1673_1674inv",
        "shuffle_direction": "3prime",
        "cross_boundaries": True,
        "source_function": "issue_324_error_normalizing_simple_inversion",
        "to_test": True,
    },
    {
        "input": "NC_000009.11:g.36233991_36233992delCAinsAC",
        "normalized": "NC_000009.11:g.36233991_36233992delinsAC",
        "shuffle_direction": "3prime",
        "cross_boundaries": True,
        "source_function": "issue_335_delins_canonicalization",
        "to_test": True,
    },
    {
        "input": "NM_000535.5:c.1673_1674delCCinsGG",
        "normalized": "NM_000535.5:c.1673_1674inv",
        "shuffle_direction": "3prime",
        "cross_boundaries": True,
        "source_function": "issue_335_delins_canonicalization",
        "to_test": True,
    },
    {
        "input": "NM_000535.5:c.1_3delATGinsCAT",
        "normalized": "NM_000535.5:c.1_3inv",
        "shuffle_direction": "3prime",
        "cross_boundaries": True,
        "source_function": "issue_335_delins_canonicalization",
        "to_test": True,
    },
    {
        "input": "NM_005877.4:c.1104_1105insGG",
        "normalized": "NM_005877.4:c.1104_1105insGG",
        "shuffle_direction": "3prime",
        "cross_boundaries": False,
        "source_function": "issue_376_ins_no_shift",
        "to_test": True,
    },
]

# hgvs-rs reinterpretation overrides: when the Rust port intentionally diverges
# from biocommons' expectation, we restore the biocommons-faithful entry.
# Currently one case: NM_000059.3:c.7790delAAG — biocommons python's
# test_c_normalizer wraps this in `with self.assertRaises(HGVSError)`, but
# hgvs-rs transcribed it as `assert_eq!("c.7791delA", format!(...))` with a
# comment that they intentionally accept where Python rejects.
HGVS_RS_REINTERPRETATIONS: dict[tuple[str, str, bool], dict] = {
    ("NM_000059.3:c.7790delAAG", "3prime", True): {
        "input": "NM_000059.3:c.7790delAAG",
        "normalized": None,
        "expects_error": True,
        "shuffle_direction": "3prime",
        "cross_boundaries": True,
        "source_function": "normalize_cds_utr_variant",
        "to_test": True,
    },
}


def apply_reinterpretations(cases: list[dict]) -> list[dict]:
    out = []
    for c in cases:
        key = (c["input"], c["shuffle_direction"], c["cross_boundaries"])
        if key in HGVS_RS_REINTERPRETATIONS:
            out.append(HGVS_RS_REINTERPRETATIONS[key])
        else:
            out.append(c)
    return out


def cmd_refresh() -> None:
    src = read_normalizer_rs_at_sha()
    cases = extract_cases_from_hgvs_rs(src)
    cases = apply_reinterpretations(cases)
    cases.extend(ISSUE_CASES)

    # Dedupe by (source_function, input, shuffle_direction, cross_boundaries).
    # If upstream later transcribes one of the hand-curated ISSUE_CASES, the
    # duplicate would silently double-count under the same identity tuple; this
    # guard keeps the first occurrence and fails loudly on a conflicting
    # duplicate so the merge is reviewed rather than masked.
    seen: dict[tuple[str, str, str, bool], dict] = {}
    for c in cases:
        key = (
            c["source_function"],
            c["input"],
            c["shuffle_direction"],
            c["cross_boundaries"],
        )
        existing = seen.get(key)
        if existing is None:
            seen[key] = c
        elif existing != c:
            raise RuntimeError(f"conflicting duplicate case for key={key}")
    cases = list(seen.values())

    # Deterministic ordering: source_function, then input, then direction,
    # then cross_boundaries. Stable across runs.
    cases.sort(
        key=lambda c: (
            c["source_function"],
            c["input"],
            c["shuffle_direction"],
            c["cross_boundaries"],
        )
    )

    payload = {
        "description": "Biocommons normalizer test cases",
        "source": "https://github.com/biocommons/hgvs-rs",
        "source_commit": HGVS_RS_SHA,
        "biocommons_upstream": BIOCOMMONS_REPO,
        "license": "Apache-2.0",
        "refreshed_at": date.today().isoformat(),
        "cases": cases,
    }

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    with CASES_JSON.open("w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False)
        fh.write("\n")

    NOTICE.write_text(
        f"""Biocommons normalizer test fixtures
===================================

Source (proximate):  https://github.com/biocommons/hgvs-rs
Commit:              {HGVS_RS_SHA}
License:             Apache-2.0
Upstream project:    {BIOCOMMONS_REPO}

cases.json is generated by scripts/refresh-biocommons-fixtures.py from
src/normalizer.rs in the hgvs-rs repo at the pinned commit above. The Rust
port transcribed biocommons/hgvs's `tests/test_hgvs_normalizer.py` cases
verbatim, except for one BRCA2 1bp/3bp-deletion mismatch case that hgvs-rs
intentionally reinterpreted as a success; we restore biocommons' assertRaises
expectation here.

Six normalize-bearing cases from biocommons/hgvs `tests/issues/test_02xx.py`
and `test_03xx.py` (issue regressions for #293, #324, #335, #376) are hand-
curated in ISSUE_CASES at the top of the refresh script since the hgvs-rs
port did not transcribe them.

The script reads src/normalizer.rs via `git show <sha>:<path>` and does not
modify the upstream checkout's working tree.

To refresh against a newer upstream commit:
    1. Update HGVS_RS_SHA in scripts/refresh-biocommons-fixtures.py
    2. pixi run python scripts/refresh-biocommons-fixtures.py refresh
    3. Re-bless the mock pin:
         BLESS_MOCK_PIN=1 cargo nextest run --features dev \\
           -E 'test(regression_under_mock_normalized) & binary(biocommons_normalize_tests)'
    4. Commit cases.json + NOTICE + mock-pin together.

Biocommons/hgvs upstream is Apache-2.0; hgvs-rs is also Apache-2.0. Both
licenses require attribution; this file satisfies that requirement.
"""
    )

    n_to_test = sum(1 for c in cases if c.get("to_test"))
    print(f"wrote {len(cases)} cases ({n_to_test} to_test) to {CASES_JSON.relative_to(REPO_ROOT)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("refresh", help="regenerate cases.json from pinned upstream SHA")
    args = parser.parse_args()
    if args.cmd == "refresh":
        cmd_refresh()


if __name__ == "__main__":
    main()
