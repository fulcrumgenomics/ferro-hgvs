#!/usr/bin/env bash
# Five-tag A/B sweep of the portable microbench suite (benches/benchmarks.rs).
# Overlays the current branch's benchmarks.rs onto each tag, runs cargo bench
# with a per-tag criterion baseline into a shared target dir, then critcmp.
#
# Projection (benches/projection.rs) is NOT swept — run it at HEAD only:
#   cargo bench --bench projection
#
# Prereqs: critcmp (cargo install critcmp), git lfs (corpus), a quiet machine.
set -euo pipefail

# Fail fast on a missing dependency before spending minutes on benchmark
# builds: a missing `critcmp` would otherwise only surface after every tag's
# `cargo bench` run has already completed.
for cmd in git cargo critcmp; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "missing required command: $cmd"; exit 1; }
done

TAGS=(v0.4.0 v0.4.1 v0.5.0 v0.6.0 HEAD)
ROOT="$(git rev-parse --show-toplevel)"
WORK="$(mktemp -d)"
SHARED_TARGET="${ROOT}/target/perf-sweep"          # shared so critcmp sees all baselines
CORPUS="${ROOT}/tests/fixtures/bulk/clinvar_hgvs_500k.json.gz"  # pinned corpus
BENCH_SRC="${ROOT}/benches/benchmarks.rs"          # the expanded portable suite

# Track the currently-live worktree so an abnormal exit (Ctrl-C, SIGTERM)
# between `git worktree add` and `git worktree remove` doesn't leave a stale
# registration under .git/worktrees that `git worktree list` would then show.
CURRENT_WT=""

# Clean up WORK on exit (normal or error) so temp dirs never accumulate.
cleanup() {
    [[ -n "${CURRENT_WT:-}" ]] && git worktree remove --force "$CURRENT_WT" 2>/dev/null || true
    git worktree prune 2>/dev/null || true
    rm -rf "$WORK"
}
trap cleanup EXIT

[ -f "$CORPUS" ] || { echo "missing corpus: $CORPUS (git lfs pull?)"; exit 1; }
cp "$BENCH_SRC" "${WORK}/benchmarks.rs"
export FERRO_BENCH_CORPUS="$CORPUS"
export CARGO_TARGET_DIR="$SHARED_TARGET"

for tag in "${TAGS[@]}"; do
    echo "=== $tag ==="
    wt="${WORK}/wt-${tag/\//_}"
    git worktree add --detach "$wt" "$tag" >/dev/null
    CURRENT_WT="$wt"
    # Overlay the portable suite (projection target is not present pre-branch; skip it).
    cp "${WORK}/benchmarks.rs" "${wt}/benches/benchmarks.rs"
    if ( cd "$wt" && cargo bench --bench benchmarks -- --save-baseline "$tag" ); then
        git worktree remove --force "$wt"
        CURRENT_WT=""
    else
        echo "bench failed at $tag"
        git worktree remove --force "$wt"
        CURRENT_WT=""
        exit 1
    fi
done

echo "=== critcmp (all tags) ==="
critcmp "${TAGS[@]}"
echo "Baselines in: ${SHARED_TARGET}/criterion"
