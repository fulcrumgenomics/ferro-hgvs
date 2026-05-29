# Mutalyzer-normalize failure patterns

Per-axis root-cause analysis of the divergences ferro-hgvs surfaces against the
imported mutalyzer normalizer fixtures (320 active cases × 8 output axes; many
cases hit multiple axes). The authoritative, continuously-updated tracker is
[#326](https://github.com/fulcrumgenomics/ferro-hgvs/issues/326); this file is
the in-repo snapshot.

> **Regenerated 2026-05-29 from a live `FERRO_MANIFEST` run.** This supersedes
> the earlier version, whose "479 FAILs" header and per-pattern dispositions
> predated the spec arbitration and were optimistic (e.g. it assumed exposing
> the c→g projection API would clear the `genomic` axis; in fact that API emits
> invalid coordinates — see #480).

## Important: "FAIL" means "differs from mutalyzer", not "ferro is wrong"

Mutalyzer is **not** HGVS ground truth. Every divergence must be arbitrated
against the spec (`assets/hgvs-nomenclature`) into one of:

- **(a) ferro bug** — fix ferro.
- **(b) ferro spec-correct / mutalyzer diverges** — annotate `accepted_divergence`
  in `cases.json` (the #335 mechanism); no ferro change.
- **(c) missing feature** — ferro doesn't implement this surface yet.
- **(d) both spec-allowed** — policy call.

## Live run (origin/main, 2026-05-29)

| Axis | FAIL | Verdict | Tracking |
|---|---:|---|---|
| `normalized` | 105 | mostly (a): `ins`→`dup`, tandem-repeat & allele-collapse detection, mito `m.`, 3′ off-by-one; plus panics. ~8–10 (b) `(GENE_v001)` selector (#121). | **#487** (edit-form), **#488** (panics) |
| `genomic` | 100 | (a)/scope: `project_to_genomic` emits **chromosome coords under the `NG_`/`LRG_` accession** — invalid HGVS regardless of mutalyzer. Axis is 100% NG/LRG-relative; cdot has no NM→NG alignment. | **#480** |
| `protein_description` | 71 → **64** | (a): routes c→g→CDS (`projector.rs`), so it sits on top of #480 and needs a `genomic_context` parent. NP accession was wrong (NM→NP inference). 7 wrapper-only cases now spec-correct (bare `NP_`); ~10 genuine prediction divergences + ~54 bare-`NM_` remain. | **PR #483** (NP fix, merged) + **PR #495** (7 bare-NP) |
| `coding_protein_descriptions` | 51 | (a): same projection path; many `got []`. | as above |
| `errors` | 54 | mostly (a): ferro **over-permissively accepts** invalid HGVS (U-in-DNA, single-pos `ins`, coord-system mismatch on `NR_`, seq/length/repeat mismatch, out-of-bounds, intronic-on-bare-`NM_`). ~7 (b) Ensembl/complex-repeat/selector. | **#486** |
| `infos` | 22 → **11** | (b): the `I*` codes are mutalyzer **internal** diagnostics, absent from the HGVS spec. Accepted; 11 genuine divergences remain (`SHUFFLE_APPLIED` over/under-emit, parse failures). | **PR #481** (shipped 22→11) |
| `rna_description` | 16 | (c): no public c.→r. consequence-prediction surface exists. | **#485** |
| `noncoding` | 8 | (c): no public c.→n. surface. | **#485** |
| **Total** | **427** | | |

Merged demotions: **PR #479** removed 49 already-passing `ins[...]` rows
(`normalized` 133→105, `genomic` 120→99); **PR #481** reclassified the 11
non-spec `infos` rows (22→11); **PR #483** fixed the protein NP accession.
In flight: **PR #495** (7 spec-correct bare-`NP_` protein cases). Committed
baseline total is now **415** (→ 408 once #495 merges).

## Regenerating this snapshot

The per-axis FAIL lists under `/tmp/ferro-xfail/<axis>.{txt,tsv}` are produced
by the single-process `regenerate_baselines` test (the per-axis `axis_*` gates
no longer write files):

```sh
FERRO_MANIFEST=/path/to/manifest.json \
  cargo nextest run --features dev -E 'test(regenerate_baselines)' \
  --run-ignored only --no-capture
```

Then refresh the committed `baseline-failures/<axis>.txt` from those outputs.
A future `scripts/group-mutalyzer-failures.py` may auto-generate the table
above from the `.tsv` files; the **verdict** column is human spec arbitration
and is maintained by hand (mirrored in #326).

## Background

- Corpus: 320 active cases × 8 output axes, ported from `mutalyzer/mutalyzer`
  `tests/variants_set.py` at the pinned upstream commit.
- Refresh the corpus: `pixi run python scripts/refresh-mutalyzer-fixtures.py refresh`.
  Bumping the pinned SHA is a deliberate per-PR action; never blind.
- Layout convention: `tests/fixtures/CORPUS_LAYOUT.md`.
