# `baseline-failures/<axis>.txt` — live-divergence input lists

Each file lists the `case.input` strings that, as of the parent PR, surface a
divergence between ferro-hgvs and mutalyzer on that axis (sorted + unique, one
input per line).

These files are **informational** — not read by the test runner or CI. The
enforced gates are the per-case dispositions in `cases.json` (XPASS-guarded) and
the `mock-pin/*.txt` regression pins. The committed live-FAIL set here is a
**non-hermetic snapshot** (it can only be reproduced from a reference manifest);
its principled retirement — seeding the untriaged rows into `cases.json` and
deleting these snapshots — is tracked by #325 / #326.

Do **not** hand-maintain counts in this directory. Per-axis disposition tallies
are derived into the generated `../failure-patterns.md`
(`cargo run --features dev --example generate_conformance_summary`); the live
FAIL set is emitted only by the nightly manifest run (under `/tmp/ferro-xfail/`).

## Regenerating after a manifest run

A burn-down PR fixes ferro behaviour, demoting inputs from "broken" to
"passing". On a dev box with the reference manifest:

```bash
cargo nextest run --features dev --test mutalyzer_normalize_tests
sort -u /tmp/ferro-xfail/<axis>.txt \
  > tests/fixtures/mutalyzer-normalize/baseline-failures/<axis>.txt
```

Regenerate the whole snapshot only after an upstream refresh
(`scripts/refresh-mutalyzer-fixtures.py refresh`), because new upstream cases
would otherwise show up as unsynchronized new FAILs.
