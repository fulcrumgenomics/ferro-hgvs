#!/usr/bin/env bash
#
# Pre-generate the gitignored spec fixtures ONCE, before a `cargo nextest run`,
# so no test process shells out to a nested `cargo run` while the suite is live.
#
# Why this exists (#1608). The two generated build artifacts
#   tests/fixtures/grammar/hgvs_spec_normalization.json
#   tests/fixtures/grammar/hgvs_spec_enumeration.json
# are `.gitignore`d (see CLAUDE.md, "Generated spec fixture"). Their consuming
# tests regenerate them on demand from inside the test body, via
# `tests/it/common/fixture_gen.rs`, which fires a nested
# `cargo run --features dev --bin generate_spec_*`. That nested build is atomic
# and in-process-locked (#1085), so a *partial* read can no longer happen — but
# under a full `cargo nextest run` many test processes start at once, and when
# the artifact is absent several of them race to fire that nested `cargo run`
# concurrently. Those siblings contend the cargo build lock (against each other
# and against other worktrees compiling), and the generator intermittently
# "exits with failure", reddening a spec-fixture consumer on the first full run
# that then passes on an isolated re-run. See the issue for the observed symptom
# (a changed total test count and three spurious `clause_ruling_index` failures).
#
# This script is wired as a nextest `[scripts.setup]` in `.config/nextest.toml`,
# so it runs once, serially, after nextest's build phase and before any test
# process starts. That removes the concurrency, not the fallback:
# `fixture_gen.rs` still regenerates on demand, so a plain `cargo test` (no
# nextest) and any run that still finds a fixture missing keep working.
#
# The binding is SCOPED to the modules that read a fixture, not `all()`, and the
# reason is this script's own failure mode: it needs `assets/hgvs-nomenclature`
# and a crate it can build. `soak` and `sweeps` have neither — they check out
# with no submodule and run pre-built binaries from an archive — so an `all()`
# binding fired there, compiled for 68.7-88.7s, and then died on an empty spec
# checkout, taking 8 soak shards and 3 sweeps to `0 tests run` (run 31931487583).
# `.config/nextest.toml` carries the full account, and
# `tests/it/spec_fixture_setup_filter.rs` keeps the binding equal to the derived
# consumer set so neither direction can drift.
#
# It complements the CI / pre-push pregeneration (#1775), which regenerates the
# same artifacts unconditionally before the archived-binary test jobs; this is
# the local-`nextest` analogue.
#
# Semantics match the on-demand fallback deliberately: generate a fixture only
# when it is ABSENT. That keeps the common case (fixtures already present from a
# prior run) near-free, and — like the fallback — leaves staleness to CI and the
# pre-push hook, which regenerate explicitly. The normalization fixture is
# generated first because the enumeration generator reads it for deduplication.

set -euo pipefail

# Run from the workspace root regardless of nextest's working directory, so the
# relative fixture paths and the default-package `cargo run` both resolve. This
# script lives in `<workspace-root>/scripts/`.
cd "$(dirname "$0")/.."

# The temp file currently being written, and the trap that removes it however we
# leave.
#
# This is an EXIT trap and NOT the `RETURN` trap it replaces, because that one
# never fired. Under `set -e` a failing command inside a function makes the shell
# EXIT rather than return, so the function's `RETURN` trap is skipped entirely —
# in exactly the case it was written for. Measured, not reasoned: the failing
# setup script on run 31931487583's shape left a
# `tests/fixtures/grammar/hgvs_spec_normalization.json.470DPv` behind in the
# checkout, and an isolated repro (`false` in place of the `cargo run`) leaves
# its probe file the same way. `INT`/`TERM` are listed so a Ctrl-C or a job-level
# timeout is covered too.
#
# A scalar rather than an array: at most one temp is live at a time, and
# expanding an EMPTY array under `set -u` is an "unbound variable" error on bash
# 3.2 — which is `/bin/bash` on macOS, i.e. the local runs this script exists to
# serve.
TEMP_FILE=""
cleanup_temp_file() {
  [[ -n "$TEMP_FILE" ]] && rm -f "$TEMP_FILE"
  return 0
}
trap cleanup_temp_file EXIT INT TERM

# Generate a fixture only if it is absent, ATOMICALLY: the generator writes its
# default output in place (`std::fs::write`, non-atomic), so writing straight to
# the destination would leave a partially-written file behind if the run were
# killed mid-write — and the `-f` guard would then treat that stump as valid and
# skip it forever. Instead write to a sibling temp and `mv` into place (atomic
# within one filesystem), mirroring `fixture_gen.rs`'s own temp-then-rename.
generate_if_absent() {
  local destination="$1"
  local generator="$2"

  [[ -f "$destination" ]] && return 0

  echo "pregenerate-spec-fixtures: generating $destination" >&2

  TEMP_FILE="$(mktemp "${destination}.XXXXXX")"
  cargo run --quiet --features dev --bin "$generator" -- --output "$TEMP_FILE"
  mv -f "$TEMP_FILE" "$destination"
  # Renamed, so there is nothing left to clean up; clearing it also keeps the
  # trap from removing the NEXT call's temp should this one have been the last.
  TEMP_FILE=""
}

# The normalization fixture first: the enumeration generator reads it.
generate_if_absent \
  "tests/fixtures/grammar/hgvs_spec_normalization.json" generate_spec_fixture
generate_if_absent \
  "tests/fixtures/grammar/hgvs_spec_enumeration.json" generate_spec_enumeration
