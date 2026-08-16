# `nightly-reference-aware/xfail-baseline.txt` — the nightly's known-failing set

`nightly-mutalyzer.yml` runs every reference-aware guard in the repository
(`REFERENCE_AWARE_SELECTION`) against a real ferro-prepared manifest. A large,
stable set of those tests is known-failing — live divergences from mutalyzer,
biocommons and hgvs-rs, plus tracked ferro bugs — and the job is deliberately
**non-gating on the failures themselves** (see the `continue-on-error` comment in
that workflow, and #1998 for why).

`xfail-baseline.txt` is the committed snapshot of that known-failing set: one
nextest test id per line (`{binary-id}::{module-path}::{test-fn}`), sorted and
unique. `#`-prefixed lines and blanks are ignored.

## What it gates

The nightly captures its actual FAIL set as JUnit
(`target/nextest/nightly-xfail/junit.xml`, the `nightly-xfail` profile in
`.config/nextest.toml`) and diffs it against this file with
`scripts/diff_nightly_xfail_baseline.py`. That diff — **not** the
`continue-on-error` step conclusion — is the nightly's gating signal:

- a **new failure** (failing now, not in this file) turns the job RED — the drift
  #1998 says a bare list cannot see;
- an **unexpected pass** (a test in this file that now passes, or is no longer
  present) also turns the job RED, so the baseline cannot silently rot;
- an **exact match** is green.

When the nightly goes red on drift, `report-failure` opens a tracking issue.

## Regenerating

This is a **non-hermetic snapshot** — it can only be reproduced on a box with the
reference manifest, and its exact membership depends on the upstream (NCBI /
Ensembl) reference the run prepared, so an upstream refresh can legitimately move
it. After a deliberate change (a burn-down PR that fixes tests, or an upstream
reference refresh), regenerate it from a manifest run:

```bash
FERRO_MANIFEST=<ref>/manifest.json \
  FERRO_REQUIRE_MANIFEST=1 \
  FERRO_ASSERT_IDEMPOTENT=1 FERRO_ASSERT_REPARSE=1 \
  FERRO_ASSERT_IN_BOUNDS=1 FERRO_ASSERT_SEQUENCE=1 \
  cargo nextest run --features dev --profile nightly-xfail -E "<REFERENCE_AWARE_SELECTION>"

python scripts/diff_nightly_xfail_baseline.py \
  --junit target/nextest/nightly-xfail/junit.xml \
  --baseline tests/fixtures/nightly-reference-aware/xfail-baseline.txt \
  --update-baseline
```

Keep the `-E` selection identical to `REFERENCE_AWARE_SELECTION` in
`nightly-mutalyzer.yml`, and keep the four oracles armed — the baseline is the
FAIL set *with those oracles on*, since the nightly runs them.

Do **not** hand-edit this file to silence a red nightly without understanding the
drift first: a new failure is a regression to investigate, not a line to delete.
