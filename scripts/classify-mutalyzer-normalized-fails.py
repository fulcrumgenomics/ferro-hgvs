#!/usr/bin/env python3
"""Classify rows of /tmp/ferro-xfail/normalized.tsv into the two patterns
that issue #335 resolves via comparator annotations:

A. ``gene-symbol-selector`` (ferro policy #121): expected and got
   differ EXCLUSIVELY in the substring between the first ``(`` and the
   matching ``)`` at the same position. Both sides must share the same
   prefix up to ``(``, the same ``)`` byte offset, and the same suffix.

B. ``dup-vs-ins`` at coding boundary (HGVS §Prioritization): expected
   ends with ``ins<bases>`` and got ends with ``dup`` (with the same
   accession + coordinate-system prefix). Plus a leniency for the
   reverse shape where expected says ``dup`` (rare).

All other rows are classified ``other`` and printed last so a reviewer
can decide whether they need additional comparator annotations or
upstream fixes.

Usage:
    python3 scripts/classify-mutalyzer-normalized-fails.py \
        [path/to/normalized.tsv]

If no path is given, defaults to ``/tmp/ferro-xfail/normalized.tsv``.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

DEFAULT_TSV = Path("/tmp/ferro-xfail/normalized.tsv")

# Capture the want/got pair embedded in the diagnostic. Diagnostics
# emitted by tests/mutalyzer_normalize_tests.rs have the form
#     expected="..." got="..."
# or
#     expected="..." err=<message>
DIAG_OK = re.compile(r'expected="(?P<want>[^"]*)" got="(?P<got>[^"]*)"')
DIAG_ERR = re.compile(r'expected="(?P<want>[^"]*)" err=')

# Match an HGVS-like description: <accession>(<selector>):<rest>.
# Both want and got must use the SAME accession + coord rest for the
# selector-only divergence to count.
SEL = re.compile(r"^(?P<acc>[^(]+)\((?P<sel>[^)]*)\)(?P<tail>:.*)$")


@dataclass
class Row:
    input: str
    diag: str
    want: str
    got: str


def parse_rows(tsv_path: Path) -> list[Row]:
    rows: list[Row] = []
    for line in tsv_path.read_text().splitlines():
        if "\t" not in line:
            continue
        inp, diag = line.split("\t", 1)
        m_ok = DIAG_OK.search(diag)
        if m_ok is not None:
            rows.append(Row(inp, diag, m_ok["want"], m_ok["got"]))
        else:
            m_err = DIAG_ERR.search(diag)
            if m_err is not None:
                # ferro errored out — neither pattern applies. Capture for
                # the "other" bucket below.
                rows.append(Row(inp, diag, m_err["want"], "<err>"))
    return rows


def is_gene_symbol_selector(want: str, got: str) -> bool:
    """Return True iff want and got differ ONLY in the parenthesized
    selector token. Both must parse as ``acc(sel):tail`` with identical
    ``acc`` and ``tail``. Either side may use a longer/shorter selector
    string — the divergence is in spelling, not structure.
    """
    mw = SEL.match(want)
    mg = SEL.match(got)
    if mw is None or mg is None:
        return False
    return mw["acc"] == mg["acc"] and mw["tail"] == mg["tail"] and mw["sel"] != mg["sel"]


# Match a description ending in ``ins<DNA-bases>`` on the expected side
# and a description ending in ``dup`` on the got side. Captured prefix
# is everything before the edit suffix (accession + position interval)
# so we can require both sides to be talking about the same position.
ENDS_INS = re.compile(r"ins[ACGTN]+$")
ENDS_DUP = re.compile(r"dup$")
INS_WITH_PREFIX = re.compile(r"^(?P<prefix>.+)ins[ACGTN]+$")
DUP_WITH_PREFIX = re.compile(r"^(?P<prefix>.+)dup$")


def is_dup_vs_ins(want: str, got: str) -> bool:
    """Return True iff the two descriptions differ only in ins-vs-dup
    edit form at the same accession+interval prefix. The HGVS-preferred
    canonical form is ``dup`` when the inserted bases duplicate the
    immediately-preceding span, so this bucket covers both the common
    case (mutalyzer keeps ``ins<bases>``, ferro picks ``dup``) and the
    documented reverse case (mutalyzer keeps ``dup``, ferro emits the
    redundant ``ins<bases>``). Prefix equality keeps us from matching
    unrelated variants that only happen to share a suffix shape.
    """
    w_ins = INS_WITH_PREFIX.match(want)
    g_dup = DUP_WITH_PREFIX.match(got)
    if w_ins is not None and g_dup is not None and w_ins["prefix"] == g_dup["prefix"]:
        return True

    # Rare reverse shape: expected uses dup, got uses ins<bases>.
    w_dup = DUP_WITH_PREFIX.match(want)
    g_ins = INS_WITH_PREFIX.match(got)
    return w_dup is not None and g_ins is not None and w_dup["prefix"] == g_ins["prefix"]


def main(argv: list[str]) -> int:
    tsv_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_TSV
    if not tsv_path.exists():
        print(f"error: {tsv_path} does not exist", file=sys.stderr)
        return 2

    rows = parse_rows(tsv_path)
    gss: list[Row] = []
    dup_ins: list[Row] = []
    other: list[Row] = []

    for r in rows:
        if r.got == "<err>":
            other.append(r)
            continue
        if is_dup_vs_ins(r.want, r.got):
            dup_ins.append(r)
        elif is_gene_symbol_selector(r.want, r.got):
            gss.append(r)
        else:
            other.append(r)

    print(f"# gene-symbol-selector (#121 ferro policy): {len(gss)}")
    for r in gss:
        print(r.input)
    print()
    print(f"# dup-vs-ins (HGVS §Prioritization): {len(dup_ins)}")
    for r in dup_ins:
        # Print "<input>\t<got>" so we can fold the corrected `normalized`
        # value into cases.json without re-running the corpus.
        print(f"{r.input}\t{r.got}")
    print()
    print(f"# other: {len(other)}")
    for r in other:
        print(r.input)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
