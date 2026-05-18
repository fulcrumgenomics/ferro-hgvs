# Mutalyzer-normalize failure patterns

This document groups the **495 FAILs** that ferro-hgvs surfaces today against
the imported mutalyzer normalizer fixtures (320 active cases × 8 output axes,
many cases hitting multiple axes). Each section is one root-cause pattern, not
one input — burn-down PRs fix the root cause and demote rows from the
per-axis `baseline-failures/<axis>.txt` files.

Raw per-axis FAIL data is rebuilt by running

    rm -rf /tmp/ferro-xfail
    cargo nextest run --features dev --test mutalyzer_normalize_tests \
      -E 'test(/^axis_/)' --no-capture --no-fail-fast

against a `cases.json` whose `xfail` arrays are empty. The test runner writes
`/tmp/ferro-xfail/<axis>.tsv` (input \t diagnostic) and `<axis>.txt` (input
per line). `scripts/group-mutalyzer-failures.py` (next PR) will re-generate
this Markdown from those TSVs.

## Per-pattern summary

| # | Pattern | Count | Axes | Existing GH refs | Arbitration verdict | Burn-down disposition |
|---|---|---:|---|---|---|---|
| 1 | `c→g API not exposed` | 83 | genomic | (none) | **ferro gap** — c./n. → g. projection is not exposed; HGVS spec requires reverse projection support to back the `equivalent_descriptions["g"]` shape | **new issue**: "expose c./n. → g. projection in `VariantProjector`"; once landed, re-run axis and the 83 fall to 0 |
| 2 | `VariantProjector only accepts g. variants` (`project: error` / `project_all: error`) | 83 (65+18) | protein_description, coding_protein_descriptions | **#310** ([PR #313](https://github.com/fulcrumgenomics/ferro-hgvs/pull/313) in flight) covers non-RefSeq protein prediction; the g.-only restriction is a separate adjacent gap | **ferro gap** — projector should accept c./n. inputs by first projecting back to g., then forward to p. | **new issue**: "VariantProjector::project should accept c./n./r. inputs"; possibly subsumed if c→g API (#1) ships first |
| 3 | `output divergence: other` | 101 | normalized, genomic | many — needs per-row spec arbitration | mixed | **bulk arbitration** below (split into sub-patterns 3a–3e) |
| 4 | `error-code mapping missing` | 50 | errors | (none) | ferro has its own `FerroError` taxonomy; mutalyzer has 30+ specific codes (`ESYNTAXUEOF`, `ENODNA`, `ELENGTHMISMATCH`, `EOUTOFBOUNDARY`, …) | **new issue**: "build mutalyzer↔ferro error-code mapping table for cross-tool diagnostics"; until then, errors axis stays in baseline-failures |
| 5 | `output divergence: coding_protein pair missing` | 33 | coding_protein_descriptions | **#310** / #313 — the gene-symbol-in-selector issue accounts for most of these once c→g exists | **ferro policy** — #121 closed deliberately *preserving* gene-symbol selector; mutalyzer uses `(NM_legacy_alias)`. The pair shape is the same; only the selector spelling differs | Once #313 merges and c→g API ships: adjust comparator to normalize selectors before matching, OR document as accepted divergence |
| 6 | `parse error` | 23 | normalized, genomic, protein_description | check #87 / #83 trackers — likely uncovered spec gaps | **needs per-row spec arbitration** — these are inputs ferro cannot parse at all | Triage each: if parseable per spec → **new issue** under #87; if intentionally rejected → move expected to `errors` axis with spec citation |
| 7 | `info-code surface not wired` | 22 | infos | (none) | ferro doesn't yet emit structured `info` codes the way mutalyzer's `IINFO_*` does | **new issue**: "structured info-code surface in normalize output (mirror W##### warnings)" |
| 8 | `gene-symbol vs alias selector (#121)` | 21 | normalized | **#121** (closed, ferro policy explicitly preserves gene-symbol selector); see also #310/#313 | **ferro policy** — these are *not bugs*. mutalyzer emits `(BRCA2_v001)`, ferro emits `(BRCA2)` | Document in `cases.json` per-case with `accepted_divergence: "ferro-policy-121-gene-symbol-selector"`; comparator will treat selector-only differences as PASS. **No issue filed; no ferro change.** |
| 9 | `3' shift differs by N bases` | 20 | normalized | overlaps **#161**, **#11**, and adjacent intronic/repeat-shift work | **needs per-row spec arbitration** — some are real ferro bugs in 3' rule, some are mutalyzer over-shifts | Triage each against HGVS §SVD-WG009; cluster into 2-4 sub-issues per shape (intronic-boundary, repeat-region, etc.) |
| 10 | `r. prediction not wired` | 16 | rna_description | **#283** (open, [PR #301](https://github.com/fulcrumgenomics/ferro-hgvs/pull/301)) audits c.→r.; **#291** (open, [PR #304](https://github.com/fulcrumgenomics/ferro-hgvs/pull/304)) related r. branch fix | **ferro gap** — runner uses a stub. After #301/#304 merge, wire the runner to ferro's r. prediction surface | Wire runner: replace the `Err("not wired")` stub with real call once those PRs merge |
| 11 | `no chromosome mapping for intronic normalization` | 12 | normalized, genomic | adjacent to **#172** (E3006 intronic ranges) and the recently-closed #98/#100 family | **ferro bug** — when input is `NG_xxx(NM_yyy):c.123-5_…`, ferro requires a chromosome mapping but the cdot mapper has the NG-based projection | **new issue**: "intronic boundary normalization with NG-prefixed transcript inputs" |
| 12 | `panic: arithmetic overflow in normalize` | 9 | normalized, genomic, protein_description | adjacent to **#87** "insertion form" rules and shuffle.rs:86 (called out in #87 body) | **ferro bug** — `c.123insG` (single-pos insertion) panics in `Normalizer::normalize` at `normalize/shuffle.rs:86` per #87's call-out | This is already documented in #87. **No new issue needed.** Add to #87 checklist if not there; demote when fixed. |
| 13 | `n. projection not wired` | 8 | noncoding | same family as #10 (r. prediction) | **ferro gap** — runner stub | Wire runner once n. projection API is settled |
| 14 | `error expected but ferro succeeded` | 7 | errors | needs per-row spec arbitration | **needs per-row spec arbitration** — mutalyzer rejects these as errors; should ferro? | Triage 7 inputs; per-row spec citation; either tighten ferro validation (new issue) or accept divergence |
| 15 | `no chromosome for boundary normalization` | 2 | normalized, genomic | same family as #11 | **ferro bug** | Subsumed into #11's new issue |
| 16 | `normalize: other error` | 2 | normalized, genomic | per-row triage | **needs per-row spec arbitration** | Triage |
| 17 | `could not infer transcript_id` | 2 | protein_description | runner-side issue, not ferro | **runner bug** — for `NG_…:g.…` style inputs the runner can't infer which transcript to project against | Fix in the runner: use `project_all` for protein axis when input has no embedded transcript |
| 18 | `output divergence: expected empty string` | 1 | normalized | upstream fixture quirk | **upstream fixture quirk** — one mutalyzer case has `"normalized": ""` which mutalyzer treats as "no normalized form" | Document in `cases.json` as accepted divergence |

**Total burn-down work:** 13 new ferro issues (rows 1, 2, 4, 6, 7, 9, 10, 11, 13, 14, 17, plus #3 splits + child issues), minus rows that map to existing in-flight PRs (#310/#313, #283/#301, #291/#304, #87, #121).

## Sub-patterns inside `output divergence: other` (101 FAILs)

This pattern is too generic to act on as-is. Eyeballing the diagnostics:

### 3a. `ins[…]` → `ins<seq>` collapse (~30)
**Examples:**
- `NG_007485.1(NM_000077.4):c.161_162ins[ATC]` → expected `c.161_162insATC`
- `NG_008939.1:g.5207_5208ins[GTCCTGTGCTC;ATTATCTGGC]` → expected `c.…insGTCCTGTGCTCATTATCTGGC`
- `NG_008939.1:g.5207_5208ins[4300_4320]` → expected `c.…insGTCCTGTGCTCATTATCTGGC` (range→sequence lookup)

**Verdict:** **ferro gap** — bracketed and range-form inserts should be canonicalized to a single concatenated sequence. The range-form `ins[start_end]` needs reference-sequence lookup, which is a separate feature.

**Disposition:** new issue: "canonicalize bracketed `ins[…]` and range-form `ins[start_end]` to flat sequence form".

### 3b. dup-vs-ins canonicalization at coding-region boundary (~5)
**Examples:**
- `NM_000143.3:c.1_2insCAT` → expected `c.1_2insCAT`, got `c.-1_2dup`
- `NM_000143.3:c.-1_1insCAT` → expected `c.-1_1insCAT`, got `c.-1_2dup`

**Verdict:** **needs spec arbitration** — HGVS spec §insertion vs §duplication is sensitive to whether the inserted seq matches flanking ref bases; ferro is canonicalizing aggressively where mutalyzer keeps insertion form.

**Disposition:** triage 2 representative cases against HGVS §SVD-WG009 + §insertion, decide policy, file 1 issue.

### 3c. Gene-symbol-vs-alias selector in `normalized` axis (~20)
Spillover from row 8 — gene-symbol selector preservation by ferro vs alias in mutalyzer.

**Disposition:** same as row 8 — `accepted_divergence` field in `cases.json`.

### 3d. 5' UTR boundary delins (~15)
**Examples:** several `NM_…:c.<positive>_<positive>delins…` where ferro produces a slightly different position pair than mutalyzer.

**Disposition:** **needs spec arbitration** + likely 1 new issue.

### 3e. Other (~30)
Tail of single-occurrence divergences that need 1:1 spec arbitration.

**Disposition:** spec-arbitrate batch; expected outcome is small (~1-3) new issues plus some `accepted_divergence` entries.

## Cross-reference matrix: existing GH refs ↔ patterns

| GH ref | State | Patterns it touches |
|---|---|---|
| **#310** / [PR #313](https://github.com/fulcrumgenomics/ferro-hgvs/pull/313) | issue open / PR open | rows 2, 5 (protein prediction, gene-symbol selector) |
| **#311** / [PR #312](https://github.com/fulcrumgenomics/ferro-hgvs/pull/312) | open / open | possibly handful of row 3e cases |
| **#121** | closed | row 8 (ferro policy precedent) |
| **#283** / [PR #301](https://github.com/fulcrumgenomics/ferro-hgvs/pull/301) | open / open | row 10 |
| **#291** / [PR #304](https://github.com/fulcrumgenomics/ferro-hgvs/pull/304) | open / open | row 10 |
| **#275** / [PR #292](https://github.com/fulcrumgenomics/ferro-hgvs/pull/292) | open / open | possibly row 9 |
| **#280** / [PR #297](https://github.com/fulcrumgenomics/ferro-hgvs/pull/297) | open / open | row 14 |
| **#83** | open, umbrella | rows 6, 9 (spec-driven gaps) |
| **#87** | open, umbrella | row 12 (explicitly documents the shuffle.rs:86 overflow) |
| **#172** | open | row 11 (intronic ranges adjacent) |
| **#129** | open | possibly m. axis cases (none surfaced so far) |
| **#10** | closed | reference: "Conflicting Mutalyzer outputs" — first ferro/mutalyzer divergence ticket; useful template |
| **#11** | closed | adjacent to row 9 (intron-spanning insertion numbering) |

## After this PR lands

1. File **one** umbrella tracking issue ("tracking: ferro-hgvs ↔ mutalyzer normalize parity") in `fulcrumgenomics/ferro-hgvs`. Body lists rows 1–18 with the disposition column above as the burn-down checklist.
2. File child issues for rows requiring new ferro work (1, 2, 4, 7, 11, 13 — and the 3a/3b/3d sub-issues). **Each child issue gets title + body approval before filing per the project's GitHub guardrail.**
3. Burn-down PRs: one per child issue. Each:
   - Drives a TDD case in `tests/normalize_tests.rs` (CI-runnable, `MockProvider`-backed)
   - Removes the corresponding inputs from `baseline-failures/<axis>.txt`
   - Cross-links the umbrella issue
4. Re-run the runner; demoted inputs become passes; XPASS guard ensures any forgotten demotions are loud.
