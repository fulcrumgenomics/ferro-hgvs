# `baseline-failures/<axis>.txt` — burn-down ledger

Each file lists the `case.input` strings that, **as of the parent PR**, surface
a divergence between ferro-hgvs and mutalyzer on that axis. The list is
sorted + unique, one input per line.

These files are **informational**, not used by the test runner. The runner
(`tests/mutalyzer_normalize_tests.rs`) asserts strictly: any divergence fails
the corresponding `axis_*` test. The baseline is the snapshot we hand to
reviewers and the umbrella tracking issue so they can see precisely what is
broken today and burn it down PR by PR.

## When to update

A burn-down PR fixes some ferro-hgvs behaviour, demoting one or more inputs
from "broken" to "passing". In the same PR:

1. Run `cargo nextest run --features dev --test mutalyzer_normalize_tests`
   on a dev box with the reference manifest.
2. For each axis that now has fewer FAILs than the committed
   `<axis>.txt`, regenerate it:

   ```bash
   sort -u /tmp/ferro-xfail/<axis>.txt \
     > tests/fixtures/mutalyzer-normalize/baseline-failures/<axis>.txt
   ```

3. Commit the smaller list alongside the `src/` fix.

## When to regenerate the whole snapshot

Only after an upstream refresh (`scripts/refresh-mutalyzer-fixtures.py
refresh`), because new upstream cases would otherwise show up as
unsynchronized new FAILs.

## Per-axis counts (current snapshot)

| Axis | FAILs |
|---|---:|
| `normalized` | 149 |
| `genomic` | 120 |
| `protein_description` | 71 |
| `coding_protein_descriptions` | 51 |
| `errors` | 57 |
| `rna_description` | 16 |
| `infos` | 22 |
| `noncoding` | 8 |
| **Total** | **494** |

See `../failure-patterns.md` for root-cause grouping and the burn-down
disposition per pattern.
