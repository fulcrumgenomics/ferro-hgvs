# `nightly-reference-aware/xfail-baseline.txt` — the nightly's known-failing set

`nightly-mutalyzer.yml` runs every reference-aware guard in the repository
(`REFERENCE_AWARE_SELECTION`) against a real ferro-prepared manifest. A large,
stable set of those tests is known-failing — live divergences from mutalyzer,
biocommons and hgvs-rs, plus tracked ferro bugs — and the job is deliberately
**non-gating on the failures themselves** (see the `continue-on-error` comment in
that workflow, and #1998 for why).

`xfail-baseline.txt` is the committed snapshot of that known-failing set: one
nextest test id per line (`{binary-id}::{module-path}::{test-fn}`), sorted and
unique.

The file is **machine-written**. `scripts/diff_nightly_xfail_baseline.py
--update-baseline` renders it, and a test in PR CI
(`tests/python/test_diff_nightly_xfail_baseline.py`) asserts the committed bytes
are exactly what that render produces — sorted, unique, one id per line, no
comments and no blanks. Annotations about *why* a test is xfail therefore belong
in this README, not in the file: a regeneration would delete them from there.

The parser is strict for the same reason. A line it cannot read as an id, a
commented-out id, and a duplicate are all errors rather than skips, because each
of them used to shrink the baseline in silence — and a silently smaller baseline
is a silently weaker gate. An id that drops out stops being expected-to-fail, so
when it stops failing nothing reports it.

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

A missing or unparseable JUnit report is also red: the diff step only runs when
the test step actually executed, so an absent report means nextest died before
writing one.

When the nightly goes red on drift, `report-failure` opens a tracking issue.

## What it does NOT gate: drift *inside* a baselined test

The baseline is a **set of test ids**, so what it detects is drift in the *set of
failing tests* — a 39th failure, or a baselined test that stops failing. It
detects nothing that happens *within* the 38.

That is not a small remainder, because the baselined tests are the highest-fan-in
ones in the selection. Each comparator axis runs a whole corpus and reports one
verdict. Measured on nightly run `31926234266`, one baselined axis carries
**30 pass / 57 FAIL** over its own inputs; if it moved to **0 pass / 87 FAIL**
the id set would be byte-identical and this gate would print
`Exact match — no drift.` So 38 of the 292 selected tests are exempt from
per-input drift, and they are the 38 with the most inputs behind them.

The per-input data is not lost — it is simply not diffed. Each comparator axis's
`AxisTally::finish` writes its FAIL inputs to `/tmp/ferro-xfail/<axis>.txt` and
`.tsv`, which the nightly uploads as the `ferro-xfail-*` artifact. So a movement
of 30/57 to 0/87 is *in the run*, on the axis's summary line and in that
artifact, and a human reading them sees it. Nothing compares them against
anything — which is the same shape as the gap #1998 recorded one level up.

Closing that would mean a second baseline per axis, keyed on an input string
rather than a test id, and a policy for what a corpus refresh does to it. That is
a separate change and is deliberately not this file's job. What this file buys is
the coarse signal that did not exist at all before #1998: a *new* test failing,
or a baselined one going green, is now visible without a human remembering
yesterday's artifact. Do not read it as covering the corpora inside.

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
