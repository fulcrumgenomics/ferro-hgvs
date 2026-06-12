#!/usr/bin/env python3
"""Refresh the curated hgvs-rs projection cases.json from a local hgvs-rs checkout.

The biocommons hgvs-rs corpus ships ~22.5k `(HGVSg, HGVSc, HGVSp)` triples under
`tests/data/mapper/*.tsv`, but 92% of them are one near-exhaustive dbSNP file for
a single gene (`gcp/DNAH11-dbSNP.tsv`, ~20.8k rows). Committing the full corpus
would be both a license/size problem and almost pure redundancy. Instead this
script curates a stratified, reproducible ~700-900-case subset:

  * every diverse per-gene / `real_cp` / `noncoding` file taken near-wholesale,
  * the four DNAH11 files reduced to a per-edit-class quota.

The full corpus is NOT committed; it stays in a local, untracked `vendor/hgvs-rs`
checkout and is read by this script only. See tests/fixtures/CORPUS_LAYOUT.md and
docs/superpowers/specs/2026-06-11-projection-cycle1c-hgvs-rs-corpus-design.md.

Usage:

    pixi run python scripts/refresh-hgvs-rs-projection-fixtures.py refresh \\
        --date 2026-06-11

    pixi run python scripts/refresh-hgvs-rs-projection-fixtures.py --check \\
        --date 2026-06-11

`--date` is required (no `datetime.now()`) so the emitted `cases.json` is
byte-reproducible and diffable. `--hgvs-rs-dir` overrides the default checkout
location.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# --- Pin to the vendored hgvs-rs checkout HEAD -------------------------------
# Asserted against `<hgvs-rs-dir>/../../../.git` HEAD when that checkout is a git
# repo, so a vendor bump cannot silently desync the curated cases.json.
HGVS_RS_SHA = "cf0e5a8d1b9a2ae3e883a58826f397aa6cb52e23"

# Default location of the mapper TSVs inside a local hgvs-rs checkout. The
# contributor convention is to clone hgvs-rs under `vendor/hgvs-rs`; this
# directory is .gitignore'd and never committed.
DEFAULT_DIR = "vendor/hgvs-rs/tests/data/mapper"

# --- Output paths -------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "hgvs-rs-projection"
CASES_JSON = FIXTURE_DIR / "cases.json"
NOTICE = FIXTURE_DIR / "NOTICE"


# Byte-deterministic I/O. The committed fixtures must be reproducible
# byte-for-byte regardless of platform locale or newline conventions, so the
# `--check` drift guard cannot pass on one machine and fail on another. Always
# read/write UTF-8 with explicit `\n` line endings (never the platform default).
def read_utf8(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_utf8_lf(path: Path, content: str) -> None:
    # newline="" disables translation so the embedded "\n"s are written verbatim.
    path.write_text(content, encoding="utf-8", newline="")


# --- Stratification quota -----------------------------------------------------
# DNAH11 is the redundant bulk: cap each (file, edit-class) bucket at this many
# rows, chosen deterministically by stable id/HGVSc sort. The other files are
# taken wholesale.
DNAH11_PER_CLASS = 20

# Files whose rows are taken near-wholesale (every usable row). Paths are
# relative to the mapper dir. `real_cp.tsv` sits at the mapper root; the rest
# live under `gcp/`.
WHOLESALE_FILES = [
    "real_cp.tsv",
    "gcp/ADRA2B-dbSNP.tsv",
    "gcp/BAHCC1-dbSNP.tsv",
    "gcp/FOLR3-dbSNP.tsv",
    "gcp/JRK-dbSNP.tsv",
    "gcp/NEFL-dbSNP.tsv",
]

# Non-coding (NR_) transcript file: rows drive the `noncoding` axis.
NONCODING_FILE = "gcp/noncoding.tsv"

# DNAH11 files: stratified by edit/position class.
DNAH11_FILES = [
    "gcp/DNAH11-dbSNP.tsv",
    "gcp/DNAH11-dbSNP-NM_001277115.tsv",
    "gcp/DNAH11-dbSNP-NM_003777.tsv",
    "gcp/DNAH11-HGMD.tsv",
]

# Alignment-discrepancy file: a different schema entirely
# (`class  position  HGVS-input  HGVS-expected`). Imported tagged and
# `to_test: false`; reserved for the 1c-iii AlignmentGap axis.
PROJ_NEAR_DISC_FILE = "proj-near-disc.tsv"


# -----------------------------------------------------------------------------
# TSV reader
# -----------------------------------------------------------------------------
def read_rows(path: Path) -> list[dict[str, str]]:
    """Parse a `(id, HGVSg, HGVSc, HGVSp, ...)` TSV.

    The first non-comment line is the header. Comment lines (starting with `#`)
    and blank lines are skipped. Empty cells become "". Extra columns beyond the
    four canonical ones (e.g. noncoding's `description`/`alternatives`) are kept
    under their header name. Each returned dict also carries `id` (possibly "").
    """
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    for raw in read_utf8(path).splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split("\t")
        if header is None:
            header = fields
            continue
        row = {header[i]: (fields[i] if i < len(fields) else "") for i in range(len(header))}
        rows.append(row)
    return rows


def edit_class(hgvsc: str) -> str:
    """Classify an HGVSc/HGVSn by edit/position for DNAH11 stratification.

    Order matters: position prefixes (5'UTR `c.-`, 3'UTR `c.*`) and intronic
    offsets (`±N`) are tested before edit tokens so a UTR/intronic substitution
    is bucketed by location, not by `>`.
    """
    body = hgvsc.split(":c.", 1)[-1].split(":n.", 1)[-1]
    if body.startswith("-"):
        return "5utr"
    if body.startswith("*"):
        return "3utr"
    # intronic: an offset like `123+4` or `123-7` (a digit-then-sign-then-digit).
    import re

    if re.search(r"\d[+-]\d", body):
        return "intronic"
    if "delins" in body:
        return "delins"
    if "dup" in body:
        return "dup"
    if "inv" in body:
        return "inv"
    if "del" in body:
        return "del"
    if "ins" in body:
        return "ins"
    if ">" in body:
        return "sub"
    return "other"


# -----------------------------------------------------------------------------
# Row -> case mapping
# -----------------------------------------------------------------------------
def base_keywords(row: dict[str, str], source_file: str) -> list[str]:
    kw = []
    rid = (row.get("id") or "").strip()
    if rid:
        kw.append(f"hgvs-rs:{rid}")
    kw.append(f"file:{source_file}")
    return kw


def row_to_case(row: dict[str, str], source_file: str) -> dict | None:
    """Map a `(HGVSg, HGVSc, HGVSp)` row to a cases.json case.

    The *input* axis is the genome-anchored HGVSg whenever present, because ferro
    projects best from a genome-pivot input (a bare HGVSc projects to
    genomic=None). The mapping:

      * g & c & p   -> input=g, coding_protein_descriptions=[[c, p]]  (g->c->p)
      * g & c & !p  -> input=g, coding=c                              (g->c only)
      * !g & c & p  -> input=c, protein_description=p                 (direct c->p)

    Rows with no projectable expectation return None.
    """
    g = (row.get("HGVSg") or "").strip() or None
    c = (row.get("HGVSc") or "").strip() or None
    p = (row.get("HGVSp") or "").strip() or None
    base = {"keywords": base_keywords(row, source_file), "to_test": True}
    if g and c and p:
        return {**base, "input": g, "coding_protein_descriptions": [[c, p]]}
    if g and c and not p:
        return {**base, "input": g, "coding": c}
    if not g and c and p:
        return {**base, "input": c, "protein_description": p}
    return None


def noncoding_row_to_case(row: dict[str, str], source_file: str) -> dict | None:
    """Map a noncoding (NR_) row: input=g, noncoding=[n.] (the g->n axis)."""
    g = (row.get("HGVSg") or "").strip() or None
    c = (row.get("HGVSc") or "").strip() or None  # NR_ accession with n. coordinate
    if not (g and c):
        return None
    return {
        "keywords": base_keywords(row, source_file),
        "to_test": True,
        "input": g,
        "noncoding": [c],
    }


def proj_near_disc_to_case(fields: list[str], source_file: str) -> dict | None:
    """Map a proj-near-disc row (`class  position  input  expected`).

    Imported `to_test: false` and tagged for the 1c-iii AlignmentGap axis; the
    expected genomic projection is carried as a keyword so 1c-iii can wire it.
    """
    if len(fields) < 4:
        return None
    klass, position, src, expected = fields[0], fields[1], fields[2], fields[3]
    if not src.strip():
        return None
    return {
        "keywords": [
            "proj-near-disc",
            f"file:{source_file}",
            f"class:{klass.strip()}",
            f"position:{position.strip()}",
            f"expected:{expected.strip()}",
        ],
        "to_test": False,
        "input": src.strip(),
    }


# -----------------------------------------------------------------------------
# Curation
# -----------------------------------------------------------------------------
def curate(mapper_dir: Path) -> tuple[list[dict], dict[str, int], dict[str, int]]:
    """Build the curated case list plus per-file and per-DNAH11-class counts."""
    cases: list[dict] = []
    per_file: dict[str, int] = {}
    per_class: dict[str, int] = {}

    def emit(case: dict | None, source_file: str) -> None:
        if case is None:
            return
        cases.append((source_file, case))
        per_file[source_file] = per_file.get(source_file, 0) + 1

    # Wholesale gene/real_cp files.
    for rel in WHOLESALE_FILES:
        for row in read_rows(mapper_dir / rel):
            emit(row_to_case(row, rel), rel)

    # Noncoding file (n. axis).
    for row in read_rows(mapper_dir / NONCODING_FILE):
        emit(noncoding_row_to_case(row, NONCODING_FILE), NONCODING_FILE)

    # DNAH11 stratified: cap each (file, edit-class) bucket at DNAH11_PER_CLASS,
    # selecting deterministically by stable (HGVSc, HGVSg, id) sort.
    for rel in DNAH11_FILES:
        rows = read_rows(mapper_dir / rel)
        rows.sort(key=lambda r: (r.get("HGVSc", ""), r.get("HGVSg", ""), r.get("id", "")))
        taken: dict[str, int] = {}
        for row in rows:
            c = (row.get("HGVSc") or "").strip()
            if not c:
                continue
            klass = edit_class(c)
            if taken.get(klass, 0) >= DNAH11_PER_CLASS:
                continue
            case = row_to_case(row, rel)
            if case is None:
                continue
            taken[klass] = taken.get(klass, 0) + 1
            per_class[klass] = per_class.get(klass, 0) + 1
            emit(case, rel)

    # proj-near-disc: different schema; reserved for 1c-iii.
    for raw in read_utf8(mapper_dir / PROJ_NEAR_DISC_FILE).splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        emit(proj_near_disc_to_case(raw.split("\t"), PROJ_NEAR_DISC_FILE), PROJ_NEAR_DISC_FILE)

    # Stable sort by (source_file, case input) for byte-stable diffs, then drop
    # the carried source_file tag (it already lives in keywords).
    cases.sort(key=lambda sc: (sc[0], sc[1]["input"]))
    return [c for _, c in cases], per_file, per_class


# -----------------------------------------------------------------------------
# Emit
# -----------------------------------------------------------------------------
def build_payload(cases: list[dict], date: str) -> dict:
    return {
        "description": "biocommons hgvs-rs projection test corpus (curated subset)",
        "source": "https://github.com/biocommons/hgvs-rs",
        "source_commit": HGVS_RS_SHA,
        "license": "Apache-2.0",
        "refreshed_at": date,
        "clusters": [],
        "cases": cases,
    }


def render_json(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=False) + "\n"


NOTICE_TEXT = f"""hgvs-rs projection test fixtures
================================

Source (proximate):  https://github.com/biocommons/hgvs-rs
Commit:              {HGVS_RS_SHA}
License:             Apache-2.0
Upstream project:    https://github.com/biocommons/hgvs

cases.json is a curated, stratified subset of the hgvs-rs mapper corpus
(`tests/data/mapper/*.tsv`, ~22.5k `(HGVSg, HGVSc, HGVSp)` triples derived from
the biocommons UTA). It is generated by
scripts/refresh-hgvs-rs-projection-fixtures.py from a local, untracked
hgvs-rs checkout at the pinned commit above. The full corpus is NOT committed:
92% of it is one near-exhaustive single-gene dbSNP file, so the script takes the
diverse per-gene / real_cp / noncoding files near-wholesale and reduces the four
DNAH11 files to a deterministic per-edit-class quota.

The case `input` is the genome-anchored HGVSg when a row has one (ferro projects
best from a genome-pivot input); the expected HGVSc/HGVSp drive the coding /
protein / coding_protein_descriptions / noncoding axes.

To refresh against a newer hgvs-rs commit:
    1. Update HGVS_RS_SHA in scripts/refresh-hgvs-rs-projection-fixtures.py
    2. pixi run python scripts/refresh-hgvs-rs-projection-fixtures.py refresh \\
         --hgvs-rs-dir <path-to-hgvs-rs>/tests/data/mapper --date <ISO-date>
    3. Re-bless the mock pins:
         BLESS_MOCK_PIN=1 cargo nextest run --features dev \\
           -E 'test(regression_under_mock) & binary(hgvs_rs_projection_tests)'
    4. Commit cases.json + NOTICE + mock-pin together.

hgvs-rs (biocommons/hgvs-rs) and its upstream biocommons/hgvs are both
Apache-2.0. The vendored checkout ships only LICENSE.txt (no separate upstream
NOTICE to carry forward per Apache-2.0 §4(d)). The Apache-2.0 license requires
attribution; this file satisfies that requirement. See the upstream LICENSE.txt
for the full text.
"""


def assert_sha(mapper_dir: Path) -> None:
    """If the checkout is a git repo, assert its HEAD matches HGVS_RS_SHA."""
    import subprocess

    # mapper_dir is <checkout>/tests/data/mapper.
    checkout = mapper_dir.parent.parent.parent
    if not (checkout / ".git").exists():
        return
    # Fail closed: if the checkout is a git repo but its HEAD cannot be read, we
    # cannot prove it matches the pin, so refuse rather than silently proceeding
    # against an unknown checkout state.
    try:
        head = subprocess.check_output(
            ["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True
        ).strip()
    except FileNotFoundError:
        sys.exit("error: git is not available; cannot verify the hgvs-rs pinned HEAD")
    except subprocess.CalledProcessError as e:
        sys.exit(f"error: failed to resolve the hgvs-rs checkout HEAD: {e}")
    if head != HGVS_RS_SHA:
        sys.exit(
            f"error: hgvs-rs checkout HEAD {head} != pinned HGVS_RS_SHA {HGVS_RS_SHA}; "
            "update the pin or check out the pinned commit before refreshing"
        )


def log_counts(cases: list[dict], per_file: dict[str, int], per_class: dict[str, int]) -> None:
    to_test = sum(1 for c in cases if c.get("to_test"))
    print(f"curated {len(cases)} cases ({to_test} to_test)")
    print("per-file:")
    for f in sorted(per_file):
        print(f"  {f:40s} {per_file[f]}")
    print("DNAH11 per-class (across the four DNAH11 files):")
    for k in sorted(per_class):
        print(f"  {k:12s} {per_class[k]}")


def cmd_refresh(mapper_dir: Path, date: str) -> None:
    assert_sha(mapper_dir)
    cases, per_file, per_class = curate(mapper_dir)
    payload = build_payload(cases, date)
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    write_utf8_lf(CASES_JSON, render_json(payload))
    write_utf8_lf(NOTICE, NOTICE_TEXT)
    log_counts(cases, per_file, per_class)
    print(f"wrote {CASES_JSON.relative_to(REPO_ROOT)} and {NOTICE.relative_to(REPO_ROOT)}")


def cmd_check(mapper_dir: Path, date: str) -> None:
    assert_sha(mapper_dir)
    cases, _, _ = curate(mapper_dir)
    want = render_json(build_payload(cases, date))
    if not CASES_JSON.exists():
        sys.exit(f"error: {CASES_JSON} does not exist; run `refresh` first")
    # Compare bytes (not decoded text) so the drift guard catches encoding/newline
    # differences, not just logical-content differences.
    got = CASES_JSON.read_bytes()
    if got != want.encode("utf-8"):
        sys.exit(
            f"error: {CASES_JSON.relative_to(REPO_ROOT)} is out of date; regenerate with `refresh`"
        )
    print(f"{CASES_JSON.relative_to(REPO_ROOT)} is up to date")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        nargs="?",
        default="refresh",
        choices=["refresh", "check"],
        help="`refresh` (regenerate) or `check` (drift-guard); default refresh",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="alias for the `check` command (drift-guard, writes nothing)",
    )
    parser.add_argument(
        "--hgvs-rs-dir",
        default=DEFAULT_DIR,
        help=f"path to the hgvs-rs mapper TSV dir (default {DEFAULT_DIR})",
    )
    parser.add_argument(
        "--date",
        required=True,
        help="ISO date written to refreshed_at (no datetime.now() for reproducibility)",
    )
    args = parser.parse_args()

    mapper_dir = Path(args.hgvs_rs_dir)
    if not mapper_dir.is_absolute():
        mapper_dir = (REPO_ROOT / mapper_dir).resolve()
    if not mapper_dir.exists():
        sys.exit(f"error: hgvs-rs mapper dir not found at {mapper_dir}")

    if args.check or args.command == "check":
        cmd_check(mapper_dir, args.date)
    else:
        cmd_refresh(mapper_dir, args.date)


if __name__ == "__main__":
    main()
