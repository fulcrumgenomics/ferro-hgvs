# Normalize corpus layout

Canonical structure for importing a third-party normalizer test corpus
(mutalyzer, biocommons, HGVS spec, etc.) into ferro-hgvs. Every imported
corpus follows the same layout so contributors learn one shape and apply
it to every fixture.

## Directory structure

```text
tests/fixtures/<upstream>-normalize/
├── cases.json                 # source of truth: upstream expected outputs
├── baseline-failures/         # current FAIL inputs per axis (manifest runs)
│   ├── normalized.txt
│   ├── genomic.txt
│   └── ...
├── mock-pin/                  # ferro's current behavior under MockProvider
│   └── normalized.txt         # one line per case: `<input>\t<ferro_output_or_error>`
├── failure-patterns.md        # spec-arbitrated disposition + linked issues
└── NOTICE                     # upstream attribution + pinned SHA + refresh

scripts/refresh-<upstream>-fixtures.py
                               # reproducible regen from pinned upstream SHA
                               # (uses `git show <sha>:<path>` against a
                               # local checkout; does NOT mutate upstream)

tests/<upstream>_normalize_tests.rs
                               # two test layers (see below)
```

## The two-layer test pattern

Each corpus's test binary runs **two** logical layers against the same
`cases.json`:

### Layer 1 — `regression_under_mock_<axis>` (CI-always)

- Runs every applicable case through ferro under `MockProvider::new()`.
- Diffs ferro's output against `mock-pin/<axis>.txt`.
- Tests refactor-side regressions in mock-mode behavior.
- **Does not assert correctness** — only stability. A `mock-pin` line like
  `INPUT\tINPUT` means "ferro returns the input verbatim under Mock"
  (no shift possible without reference bases). That's the pinned baseline.
- CI runs this every PR.
- Regenerate after intentional behavior changes:
  `BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock_<axis>)'`

### Layer 2 — `axis_<axis>` (manifest-or-skip)

- Runs every applicable case through ferro with a real
  `MultiFastaProvider` loaded from the manifest at `FERRO_MANIFEST` or
  one of the well-known paths.
- Strict-asserts ferro's output equals `cases.json` expected.
- Tests correctness against upstream truth.
- When the manifest is absent (e.g. GitHub Actions CI), the test
  reports `skipping — no manifest` and exits 0.
- Divergences fail the test on dev boxes; FAIL inputs are written to
  `/tmp/ferro-xfail/<axis>.{txt,tsv}` for burn-down tracking.

The committed `baseline-failures/<axis>.txt` is the static ledger of
known FAILs at the time of the most recent burn-down PR. It is
informational — the runner does not consult it — and is meant for
contributors to diff `/tmp/ferro-xfail/<axis>.txt` against, identifying
new regressions vs. known gaps.

## Why two layers?

| | Layer 1 (Mock) | Layer 2 (Manifest) |
|---|---|---|
| Catches refactor regressions | ✓ | partial (only if manifest run) |
| Catches correctness gaps | ✗ | ✓ |
| Runs in CI without external data | ✓ | ✗ (skips) |
| Pin shape | ferro behavior (ephemeral) | upstream truth (durable) |
| Source of truth file | `mock-pin/<axis>.txt` | `cases.json` |

Each layer answers a different question. Both are necessary.

## File schemas

### `cases.json`

Source of truth: upstream expected outputs. **Never contains ferro's
current behavior** — that pollutes the durable record and goes stale.

```json
{
  "description": "string",
  "source": "string",                  // upstream repo URL
  "source_commit": "string",           // pinned SHA
  "license": "string",                 // upstream license (MIT, Apache, etc.)
  "refreshed_at": "ISO8601 string",
  "cases": [
    {
      "keywords": ["string", ...],      // upstream-specific tags
      "input": "string",                 // HGVS variant
      "normalized": "string|null",       // upstream-expected normalize() output
      "genomic": "string|null",          // upstream-expected g. projection
      "protein_description": "string|null",
      "coding_protein_descriptions": [["c.", "p."], ...] | null,
      "rna_description": "string|null",
      "noncoding": ["string", ...] | null,
      "errors": ["CODE", ...] | null,    // upstream-expected error codes
      "infos": ["CODE", ...] | null,
      "to_test": boolean                 // default true; false to disable a row
    }
  ]
}
```

### `mock-pin/<axis>.txt`

Ephemeral. One line per `to_test`-and-axis-applicable case:

```text
<input>\t<ferro_output_or_error>
```

`ferro_output_or_error` is either:
- ferro's normalize output (`{normalized_variant}`)
- `"parse error: <message>"`
- `"normalize error: <message>"`

Inputs appear in `cases.json` order. File ends with a trailing newline.

### `baseline-failures/<axis>.txt`

Informational (not read by the runner). One line per input that
currently FAILs under `axis_<axis>` with a real manifest. Inputs in
arbitrary order. Updated by burn-down PRs that demote inputs as ferro
is fixed.

### `failure-patterns.md`

Free-form markdown grouping current FAILs into root-cause patterns,
each tagged with disposition:
- `ferro-bug` — file a follow-up issue, fix in a burn-down PR
- `upstream-bug` — edit the expected output in `cases.json` with a
  `spec_citation` note
- `accepted-divergence` — ferro policy intentionally differs; record
  reason
- `spec-ambiguous` — HGVS spec doesn't arbitrate; document and accept

Each disposition cross-references the HGVS spec section that arbitrates
the case + any existing or to-be-filed ferro issues.

### `NOTICE`

Upstream attribution (required by most upstream licenses) + pinned SHA
+ refresh procedure.

## Refresh workflow

```sh
# 1. Update <UPSTREAM>_SHA at the top of scripts/refresh-<upstream>-fixtures.py
# 2. Regenerate cases.json from the new SHA:
pixi run python scripts/refresh-<upstream>-fixtures.py refresh
# 3. Rebless the mock-pin against the (possibly larger/smaller) case set:
BLESS_MOCK_PIN=1 cargo nextest run --features dev -E 'test(regression_under_mock_normalized)'
# 4. Run with a manifest if you have one; diff /tmp/ferro-xfail/<axis>.txt
#    against baseline-failures/<axis>.txt to identify new regressions or
#    newly-fixed inputs.
# 5. Commit cases.json + mock-pin/ + NOTICE together; commit baseline-failures
#    + failure-patterns.md updates separately if applicable.
```

## Per-corpus extensions

Some corpora's source tests run the normalizer under multiple configs
within the same fixture (e.g. biocommons exercises every input under
both 3'-shifting and 5'-shifting × cross-boundaries on/off). In that
case, the corpus's case schema adds `shuffle_direction` and
`cross_boundaries` fields per row, and its `mock-pin/<axis>.txt`
format extends to include the config in the pin line:

```
<input>\t<shuffle_direction>\t<cross_boundaries>\t<ferro_output_or_error>
```

This is a per-corpus convention extension, not a corpus-layout breaking
change. The harness for each corpus knows its own pin format. Mutalyzer's
corpus assumes a single fixed config (mutalyzer's own defaults) and
therefore uses the simpler `<input>\t<ferro_output_or_error>` shape.

## Consumers

| Corpus | Path | Per-case config? | Tracking issue |
|---|---|---|---|
| mutalyzer | `tests/fixtures/mutalyzer-normalize/` | no (fixed config) | — |
| biocommons | `tests/fixtures/biocommons-normalize/` | yes (`shuffle_direction`, `cross_boundaries`) | — |

Add a row when importing a new corpus.
