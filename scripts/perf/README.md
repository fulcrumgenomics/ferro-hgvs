# Performance sweep

`sweep_tags.sh` overlays the swept `benches/benchmarks.rs` onto each release
tag (v0.4.0 → HEAD), runs `cargo bench --save-baseline <tag>` into a shared
target dir, and prints a `critcmp` table localizing any regression to a release
window. Corpus is pinned via `FERRO_BENCH_CORPUS` so all tags parse identical
input.

## Version-floor contract — which bench runs against which tags

Each bench target declares a version floor in its header doc-comment. The rule:
**the sweep only measures operations that existed at v0.4.0**, so a slowdown is a
real regression. Anything added after v0.4.0 has no v0.4.0 baseline (it cannot
"regress"), so it lives in a HEAD-only baseline target that just establishes a
forward floor.

| Bench target | Floor | Role | Run |
|---|---|---|---|
| `benches/benchmarks.rs` | **v0.4.0** | **SWEPT** — regression detection (parsing, c./p. normalize on `with_test_data`, allele scaling, real-corpus throughput) | `sweep_tags.sh` (all tags) |
| `benches/baseline_head.rs` | post-v0.4.0 | HEAD baseline — enriched genomic/intronic normalize (`Transcript::new` signature drift), W30xx correctors, trans-allele parsing | `cargo bench --bench baseline_head` |
| `benches/projection.rs` | v0.6.0 | HEAD baseline — `project` module is new in v0.6.0 | `cargo bench --bench projection` |

When adding a case: if it does not both **compile and resolve at v0.4.0** (verify
with the overlay check in `benchmarks.rs`'s header), it belongs in a baseline
target, not the swept one.

## Run
```bash
cargo install critcmp        # once
git lfs pull                 # ensure the 500k corpus is real, not a pointer
scripts/perf/sweep_tags.sh   # ~5 release builds (LTO); run on a quiet machine
```

## HEAD baseline targets (not swept)
These cover post-v0.4.0 features; run them at HEAD to establish a forward floor:
```bash
cargo bench --bench baseline_head   # enriched g./intronic normalize, correctors, trans-allele
cargo bench --bench projection      # g→c projection (v0.6.0+)
```

## Reading output
`critcmp` columns are tags; rows are bench ids. A row whose time jumps between
two adjacent tag columns points at the PRs merged in that window. Investigate
those, not the whole range.
