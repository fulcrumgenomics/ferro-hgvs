# Normalize corpus layout

Canonical structure for importing a third-party normalizer test corpus
(mutalyzer, biocommons, HGVS spec, etc.) into ferro-hgvs. Every imported
corpus follows the same layout so contributors learn one shape and apply
it to every fixture.

## Directory structure

```text
tests/fixtures/<upstream>-normalize/
├── cases.json                 # source of truth: upstream expected outputs + dispositions
├── reference-windows.json     # GENERATED hermetic fixture (biocommons): the exact
│                              # reference bases the manifest pass touches; never hand-edit
├── baseline-failures/         # current FAIL inputs per axis (manifest runs)
│   ├── normalized.txt
│   ├── genomic.txt
│   └── ...
├── mock-pin/                  # ferro's current behavior under MockProvider
│   └── normalized.txt         # one line per case: `<input>\t<ferro_output_or_error>`
├── empty-projection/          # per-axis empty/degenerate-projection count baseline (manifest runs)
│   └── <axis>.count           # line 1 = count, optional line 2 = the reference identity it
│                              # was blessed against (#764); a rise fails the axis test (#651)
│                              # only when the live reference matches, else it skips (see below)
├── failure-patterns.md        # GENERATED summary (cluster taxonomy + tallies); never hand-edit
└── NOTICE                     # upstream attribution + pinned SHA + refresh

scripts/refresh-<upstream>-fixtures.py
                               # reproducible regen from pinned upstream SHA
                               # (uses `git show <sha>:<path>` against a
                               # local checkout; does NOT mutate upstream)

examples/extract_<upstream>_windows.rs
                               # GENERATOR for reference-windows.json: records the
                               # bases a manifest-backed normalize pass reads (#478)

tests/<upstream>_normalize_tests.rs
                               # test layers (see below)
```

## The three-layer test pattern

Each corpus's test binary runs **three** logical layers against the same
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
- **Output-quality gate (#651):** the burn-down compares which inputs fail,
  not the *quality* of each row's output, so a populated projection silently
  degrading to empty (e.g. a framed `NG_(NM_):c…/NG_(NP_):p…` pair collapsing
  to `[]`) stays in the same FAIL/annotated bucket and the membership diff
  reports "0 regressions". To catch this, each axis also counts rows where
  ferro produced an **empty/degenerate projection** (no protein predicted, or
  `project_variant_all` yielded zero pairs) and gates a rise in that count
  against the committed `empty-projection/<axis>.count` baseline. Because the
  count is **reference-dependent** (a differently-prepared reference legitimately
  yields a different number), each baseline is *pinned* to the reference it was
  blessed against (#764): line 1 is the count, optional line 2 is the reference
  identity. The gate enforces a rise (panics on a regression) only when the live
  reference matches the pinned identity; on a mismatch — or for a legacy
  single-line file with no identity — it **skips with a notice** instead of
  firing on reference drift. An **absent** file is the reference-independent
  "0 empty projections expected" claim and is always enforced. Regenerate (and
  re-pin) after an intentional change with
  `BLESS_EMPTY_PROJECTION=1 cargo nextest run --features dev -E 'test(axis_)'`.

The committed `baseline-failures/<axis>.txt` is the static ledger of
known FAILs at the time of the most recent burn-down PR. It is
informational — the runner does not consult it — and is meant for
contributors to diff `/tmp/ferro-xfail/<axis>.txt` against, identifying
new regressions vs. known gaps. Once a corpus is fully annotated with
per-case dispositions (`accepted_divergence` / `known_bug` / `improvement`,
see `cases.json` below) the ledger holds **no** unannotated FAILs and is
empty — every divergence lives as a disposition the harness enforces.

### Layer 3 — `axis_<axis>_hermetic` (CI-always merge gate; biocommons)

The same correctness assertion as Layer 2, but against a
`WindowProvider` built from the committed, generated `reference-windows.json`
instead of the out-of-band manifest — so it runs **on every PR in CI** with
zero external data and is the actual merge gate (#478 pillar 4). The fixture
captures *exactly* the reference bases the manifest pass reads (transcripts
whole; padded, clamped genomic windows), so this layer agrees with Layer 2
row-for-row. Layer 2 is demoted to a nightly/dev tier that catches anything the
committed windows miss and is the source `extract_<upstream>_windows` regenerates
from. Pioneered by the biocommons corpus; the other corpora can adopt the same
shape.

## Why these layers?

| | Layer 1 (Mock) | Layer 2 (Manifest) | Layer 3 (Hermetic) |
|---|---|---|---|
| Catches refactor regressions | ✓ | partial (only if manifest run) | ✓ |
| Catches correctness gaps | ✗ | ✓ | ✓ |
| Runs in CI without external data | ✓ | ✗ (skips) | ✓ |
| Is the merge gate | ✗ | ✗ | ✓ |
| Pin shape | ferro behavior (ephemeral) | upstream truth (durable) | upstream truth (durable) |
| Source of truth file | `mock-pin/<axis>.txt` | `cases.json` | `cases.json` + `reference-windows.json` |

Each layer answers a different question.

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

### `empty-projection/<axis>.count`

The reference-pinned output-quality budget (#651/#764) consumed by Layer 2
manifest runs. Up to two lines:

```text
<count>            # line 1: empty/degenerate-projection count baseline
<reference-id>     # line 2 (optional): the reference identity it was blessed against
```

- **Line 1** is the blessed count. An existing file whose first line is not a
  `usize` is a hard error (corruption surfaces immediately).
- **Line 2** is the reference identity (#764): a deterministic, path-independent
  FNV-1a digest of a content signature of the live prepared reference
  (`prepared_at`, `transcript_count`, and the basenames of the cdot / genome /
  transcript artifacts). Blessed by `BLESS_EMPTY_PROJECTION=1`.

Gate semantics:

- **Absent file** → the axis expects zero empty projections (a
  reference-independent claim) → enforced against `0`.
- **Pinned baseline whose identity matches the live reference** → enforce the
  committed count; a rise panics as a populated→empty regression.
- **Pinned baseline whose identity differs** (reference drift), or a **legacy
  single-line file with no identity** → the gate **skips with a notice** rather
  than firing, since the absolute count is not comparable across references.
  The skip notice distinguishes the two: drift → "re-bless on this reference";
  legacy → "migrate to the pinned format". Re-bless with
  `BLESS_EMPTY_PROJECTION=1 cargo nextest run --features dev -E 'test(axis_)'`.

### `failure-patterns.md`

**Generated** — never hand-edit. Produced by `cargo run --features dev
--example generate_conformance_summary` as a derived view over the `cases.json`
cluster taxonomy + per-case dispositions, and verified in CI with `-- --check`
so it cannot drift (#509). It groups each tracked divergence under its
root-cause cluster (title + HGVS spec citation) and emits per-axis disposition
tallies. It deliberately carries **no** live FAIL count: that set is
non-hermetic (it needs the reference manifest) and is emitted only by the
nightly run. To change what it shows, edit the `clusters` registry or the
disposition annotations in `cases.json` and regenerate.

### `reference-windows.json`

**Generated** — never hand-edit. The hermetic fixture behind Layer 3, produced
by `cargo run --features dev --example extract_<upstream>_windows` from the
manifest. A `RecordingProvider` wraps the real `MultiFastaProvider`, runs the
same per-case normalize loop the test does, and records every transcript
resolved and every genomic range read; transcripts are stored whole, genomic
accesses are clustered into disjoint windows padded by a safety margin, and each
contig's true length is captured so `WindowProvider` reproduces the manifest
provider's end-clamping (the short read that drives `CanonicalSplitSkipped` /
W5003) while still erroring on a read inside the contig but outside a captured
window — a genuine too-narrow extraction.

```json
{
  "description": "string",                 // provenance + regen command
  "captured_from": "string",               // manifest prepared_at
  "contig_lengths": { "<contig>": 0 },     // true length per windowed contig
  "transcripts": [ /* whole Transcript records (serde) */ ],
  "genomic": [ { "contig": "string", "start": 0, "bases": "ACGT…" } ]
}
```

Unlike `failure-patterns.md`, its `--check` needs the manifest (it derives from
reference bases, not `cases.json`), so it is a local/nightly guard, not a per-PR
CI check; the per-PR gate is Layer 3 consuming the committed file. Regenerate
whenever `cases.json` or normalize behavior changes the reference access set.

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
# 5. Commit cases.json + mock-pin/ + NOTICE together; regenerate the GENERATED
#    failure-patterns.md (example generate_conformance_summary) in the same PR
#    as any cases.json disposition/cluster change; commit baseline-failures
#    updates separately if applicable.
# 6. If the corpus has a hermetic Layer 3, regenerate reference-windows.json
#    from the manifest (example extract_<upstream>_windows) whenever the case
#    set or normalize behavior changes the reference access set, and confirm
#    axis_<axis>_hermetic agrees with the manifest Layer 2 run.
```

## Per-corpus extensions

Some corpora's source tests run the normalizer under multiple configs
within the same fixture (e.g. biocommons exercises every input under
both 3'-shifting and 5'-shifting × cross-boundaries on/off). In that
case, the corpus's case schema adds `shuffle_direction` and
`cross_boundaries` fields per row, and its `mock-pin/<axis>.txt`
format extends to include the config in the pin line:

```text
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
| hgvs-rs-projection | `tests/fixtures/hgvs-rs-projection/` | no (fixed config) | — |

Add a row when importing a new corpus.

The hgvs-rs-projection corpus additionally carries a data-source provenance
doc, `tests/fixtures/hgvs-rs-projection/ALIGNMENT_SOURCE_DIVERGENCE.md`. It
backs the two `accepted_divergence` clusters (`alignment-source-skew`,
`transcript-selection-vs-uta`): ferro projects on NCBI-canonical RefSeq-GFF
alignments while the corpus expectations come from the 2021 biocommons UTA
snapshot, and the doc records the per-transcript UTA splign CIGARs and the
gate (ungapped → 0% divergence, gapped/boundary → 39–100%) that the structural
quarantine in `tests/hgvs_rs_projection_tests.rs` keys off.
