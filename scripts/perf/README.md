# Performance sweep

`sweep_tags.sh` overlays the expanded `benches/benchmarks.rs` onto each release
tag (v0.4.0 → HEAD), runs `cargo bench --save-baseline <tag>` into a shared
target dir, and prints a `critcmp` table localizing any regression to a release
window. Corpus is pinned via `FERRO_BENCH_CORPUS` so all tags parse identical
input.

## Run
```
cargo install critcmp        # once
git lfs pull                 # ensure the 500k corpus is real, not a pointer
scripts/perf/sweep_tags.sh   # ~5 release builds (LTO); run on a quiet machine
```

## Projection (HEAD only)
Projection is v0.6.0-new and not part of the sweep:
```
cargo bench --bench projection
```

## Reading output
`critcmp` columns are tags; rows are bench ids. A row whose time jumps between
two adjacent tag columns points at the PRs merged in that window. Investigate
those, not the whole range.
