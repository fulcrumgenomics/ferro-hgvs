# Benchmarking

ferro ships a benchmark harness for measuring its own throughput and comparing it against other HGVS
tools (mutalyzer, biocommons/hgvs, hgvs-rs). This guide covers running it; for the published figures
and full method, see
[`docs/BENCHMARK_RUNBOOK.md`](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/BENCHMARK_RUNBOOK.md).

## Timing ferro alone

The main `ferro` binary needs no special build. Prepare a reference (see
[Reference data](reference-data.md)), then time a normalize run:

```bash
ferro prepare --output-dir data/ferro
ferro check --reference data/ferro
ferro normalize -i patterns.txt --reference data/ferro -j 8 -t timing.json
```

`-t/--timing` writes timing information to JSON. Two things matter for a fair number:

- **Warm the cache first.** `ferro check --reference data/ferro --build-cache` builds the on-disk
  cdot cache as a setup step, so the one-time build does not land inside the timed region.
- **Sort the input** by transcript accession (or genomic position). ferro caches each resolved
  transcript, so sorted input keeps the working set resident and is markedly faster on large
  batches.

## Comparing against other tools

Tool comparison uses the separate `ferro-benchmark` binary, built with the `benchmark` feature:

```bash
cargo build --release --features benchmark
```

`ferro-benchmark` prepares each tool's reference data (reusing ferro's for transcripts), runs
parse/normalize, and compares results:

```bash
# Prepare another tool against the same transcript data
ferro-benchmark prepare mutalyzer --ferro-reference data/ferro --output-dir data/mutalyzer

# Normalize a shared pattern set with each tool
ferro-benchmark normalize mutalyzer -i patterns.txt -o mutalyzer.json \
  --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf

# Compare two tools' outputs
ferro-benchmark compare results normalize ferro.json mutalyzer.json -o comparison.json
```

Supported tool values: `ferro`, `mutalyzer`, `biocommons`, `hgvs-rs`, and `all`.

The Python-based tools (mutalyzer, biocommons/hgvs, seqrepo) run in a [pixi](https://pixi.sh)
environment defined by the repository's `pixi.toml`; run `pixi shell` to activate it before
preparing them.

## Reading the numbers honestly

All tools are benchmarked **offline** — reference data preloaded locally, network disabled — so the
figures measure parse/normalize compute, not per-variant I/O. Out of the box the other tools resolve
each variant against a remote service (a network round-trip per variant), which is far slower than
any local figure. That local, offline setup is exactly what `ferro prepare` builds. When you quote a
benchmark, say which configuration it was measured in, and confirm the run **succeeded and produced
correct output** before trusting a timing — a crashed run still prints a wall-clock line.
