# `baseline-failures/<axis>.txt` — live-divergence input lists

Each file lists the `case.input` strings that, as of the parent PR, surface a
divergence between ferro-hgvs and the hgvs-rs oracle on that projection axis
(sorted + unique, one input per line).

These files are **informational** — not read by the test runner or CI. The
enforced gates are the per-case dispositions in `cases.json` (XPASS-guarded) and
the `mock-pin/*.txt` regression pins. The committed live-FAIL set here is a
**non-hermetic snapshot** (it can only be reproduced from a reference manifest).

The files are **seeded empty** at corpus import: no divergence has been triaged
into a disposition yet, and the live FAIL set on a manifest-equipped box is
large and un-categorized. Burn-down PRs that fix ferro behaviour — or that bless
divergences as `accepted_divergence` / `improvement` / `known_bug` in
`cases.json` — populate these snapshots from a real manifest run.

Do **not** hand-maintain counts in this directory. Per-axis disposition tallies
are derived into the generated `../failure-patterns.md`
(`cargo run --features dev --example generate_conformance_summary`); the live
FAIL set is emitted only by a manifest run (under `/tmp/ferro-xfail/`).

## Regenerating after a manifest run

On a dev box with the reference manifest:

```bash
cargo nextest run --features dev -E 'binary(hgvs_rs_projection_tests)'
sort -u /tmp/ferro-xfail/<axis>.txt \
  > tests/fixtures/hgvs-rs-projection/baseline-failures/<axis>.txt
```

Regenerate the whole snapshot only after an upstream refresh
(`scripts/refresh-hgvs-rs-projection-fixtures.py refresh`), because new upstream
cases would otherwise show up as unsynchronized new FAILs.
