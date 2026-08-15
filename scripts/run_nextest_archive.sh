#!/usr/bin/env bash
#
# Run `cargo nextest run` under a wall-clock bound, and — when that bound fires —
# name the phase the run stopped in.
#
# WHY (#2000). `Test (2/6)` stalled twice in three hours. Both job logs stop at
# the same line, and a healthy run of the same shard shows what should follow:
#
#     healthy   Extracted 66 files ... in 1.09s
#               Starting 1794 tests across 24 binaries      <- 0.101s later
#     stalled   Extracted 66 files ... in 0.82s
#               (nothing, for 72 minutes)
#
# `Starting N tests across M binaries` is printed once nextest has enumerated
# every archived binary — by executing each one with `--list` — and before any
# test body runs. Neither stall reached it, so no test ever started: both were in
# the binary-list phase. The runner's orphan reaping agrees and names the culprit,
# `cargo-nextest` plus exactly one test binary, a DIFFERENT one each time
# (`it-…`, then `ferro_hgvs-…`).
#
# That rules out the three things the issue proposed. A per-test `slow-timeout`
# with `terminate-after` governs test execution, which never began. `--status-level`
# reports on tests, of which there were none. And "which tests land in shard 2" is
# the wrong question: listing enumerates every binary identically on every shard,
# with `--partition` and `-E` applied to the resulting list, so no shard lists less
# than any other.
#
# The root cause — why one child of `cargo-nextest` stops returning, transiently,
# on a tree that runs in 91s on the next attempt — is NOT established. What this
# script buys is that the next occurrence is bounded in minutes rather than by
# GitHub's 6-hour default, and that it says which phase it died in, so the answer
# above does not have to be re-derived from two cancelled logs a third time.
#
# It does NOT cover the whole job. `actions/download-artifact` fetches the
# archive in its own step, outside this wrapper, so a hang there is still bounded
# only by the job's `timeout-minutes` and still arrives as a bare cancellation
# with no annotation. Both recorded #2000 stalls were past that point, so it is
# out of scope here — but do not read a green wrapper as "the job is covered".
#
# Usage:
#   scripts/run_nextest_archive.sh <nextest args...>
#   scripts/run_nextest_archive.sh --classify <logfile>   # used by the guard test
#
# Knob: NEXTEST_STALL_TIMEOUT (default 25m). Sized against a measured 70-94s
# baseline for the shard that stalled, and against stalls of 21 and 73 minutes.
# It is a stall detector, not a performance budget — keep it far above the slowest
# legitimate job so it can never fire on a merely slow run.

set -uo pipefail

# The two lines nextest prints on the way in, and the phases they delimit. A
# healthy archive run prints them in this order, 0.101s apart:
#
#     Extracted 66 files ... in 1.09s          <- archive read + unpack finished
#     Starting 1794 tests across 24 binaries   <- binary enumeration finished
#
# So the run has THREE phases, not two, and the log says which one it died in:
# before `Extracted` is the archive read; between the two is the binary-list
# phase; after `Starting` is test execution.
#
# `tests?` / `files?` are load-bearing: nextest pluralises, so a one-test
# selection prints `Starting 1 test across 1 binary`. A pattern requiring
# `tests` classifies that run as a list-phase stall — i.e. reports a false
# #2000 — and the singular form is exactly what a bisected or `-E`-narrowed
# re-run produces, which is when someone is most likely to be reading this
# output.
STARTED_RUNNING='Starting [0-9]+ tests? across'
EXTRACTION_FINISHED='Extracted [0-9]+ files? to'

# Strip SGR colour codes before matching. `ci.yml` sets `CARGO_TERM_COLOR: always`,
# so the line arrives interleaved —
#     \e[32;1m    Starting\e[0m \e[1m1794\e[0m tests across ...
# — and a grep over the raw bytes returns 0 on a HEALTHY run, which would report
# every genuine in-test hang as a list-phase stall. Pinned against the recorded
# bytes of the real runs by `tests/it/issue_2000_nextest_stall_phase.rs`.
strip_sgr() {
    perl -pe 's/\e\[[0-9;]*m//g' "$1" 2>/dev/null
}

# Four verdicts, because the annotation is the entire deliverable and it must
# never be more confident than its evidence.
#
# `unknown` is not defensive padding: without it an unreadable or empty log
# grep-misses everything and so reads as the earliest phase — i.e. a `mktemp` or
# `tee` failure would be reported as a real stall with no run behind it.
#
# `pre-list` is the same argument one step in, and is the reason this function
# is not a two-way test on `Starting` alone. Extraction happens INSIDE the
# wrapped command (`--extract-to` is a `cargo nextest run` flag), so a wedged
# archive read, a disk stall or a `zstd` hang produces a log with neither line —
# which a two-way test calls `list-phase` and announces as "the #2000 shape",
# telling the reader a per-test timeout cannot see it and asking them to file
# the URL on #2000. None of that is true of a stall in extraction, and it would
# pollute the one issue this script exists to make self-serving. Both recorded
# #2000 stalls reached `Extracted` and stopped there, so requiring that line is
# what makes the `list-phase` verdict match the incidents it names.
classify() {
    if [ ! -s "$1" ]; then
        echo "unknown"
    elif strip_sgr "$1" | grep -qE "$STARTED_RUNNING"; then
        echo "execution-phase"
    elif strip_sgr "$1" | grep -qE "$EXTRACTION_FINISHED"; then
        echo "list-phase"
    else
        echo "pre-list"
    fi
}

if [ "${1:-}" = "--classify" ]; then
    if [ -z "${2:-}" ]; then
        echo "usage: run_nextest_archive.sh --classify <logfile>" >&2
        exit 2
    fi
    classify "$2"
    exit 0
fi

LIMIT="${NEXTEST_STALL_TIMEOUT:-25m}"
# The `X`s must be the LAST characters of the template. With a `.log` suffix after
# them GNU/BusyBox `mktemp` fails outright (`Invalid argument`) — so every CI step
# would have died here — while BSD `mktemp` quietly returns the template verbatim,
# giving concurrent runs one shared fixed filename. Neither is visible on a green
# run, which is why this carries a comment rather than a suffix.
LOG="$(mktemp "${RUNNER_TEMP:-/tmp}/nextest-run-XXXXXX")"
# `rm` on the normal path alone leaks the file whenever the job-level timeout
# SIGTERMs us — which is precisely the scenario this script is written for.
trap 'rm -f "$LOG"' EXIT

# `timeout(1)` is coreutils. CI is Linux so it is always present there, but a
# stock macOS has only `gtimeout`, and a bare `command not found` exits 127 —
# which the caller reads as an ordinary test failure. Degrade to an unbounded
# run with a stated reason rather than failing a local reproduction.
#
# Deliberately a function rather than an argv ARRAY. The obvious spelling builds
# `RUNNER=(timeout …)` / `RUNNER=()` and expands `"${RUNNER[@]}"`, and expanding
# an EMPTY array under `set -u` is an "unbound variable" error on bash 3.2 —
# which is `/bin/bash` on macOS, i.e. the empty-array branch would die on exactly
# the platform the fallback exists to serve. Verified, not assumed.
run_nextest() {
    if command -v timeout > /dev/null 2>&1; then
        timeout --signal=TERM --kill-after=60s "$LIMIT" cargo nextest run "$@"
    else
        cargo nextest run "$@"
    fi
}

if ! command -v timeout > /dev/null 2>&1; then
    echo "run_nextest_archive.sh: no timeout(1) in PATH (install coreutils); running unbounded" >&2
fi

# `tee` keeps the live log identical to an unwrapped run — the output a reviewer
# reads is unchanged — while giving the classifier something to read afterwards.
run_nextest "$@" 2>&1 | tee "$LOG"
status=${PIPESTATUS[0]}

# 124 is GNU timeout's own report; 137 (128+SIGKILL) is `--kill-after` finishing
# the job when SIGTERM was not enough. Reading only the first files a hard-killed
# stall as an ordinary test failure.
if [ "$status" -eq 124 ] || [ "$status" -eq 137 ]; then
    case "$(classify "$LOG")" in
        pre-list)
            echo "::error title=nextest stalled BEFORE extracting the archive::The run exceeded ${LIMIT} without printing its \`Extracted N files\` line, so it stopped while reading or unpacking the archive — before nextest began enumerating any binary. This is NOT the #2000 shape: both recorded #2000 stalls reached \`Extracted\` and stopped after it, so do not file this run on #2000. Check the archive artifact and the runner's disk first."
            ;;
        list-phase)
            echo "::error title=nextest stalled BEFORE running any test::The run printed its \`Extracted N files\` line and then never printed \`Starting N tests across M binaries\`, so extraction finished and it stalled while enumerating the archived test binaries — no test body had started. This is the #2000 shape: a per-test timeout cannot see it, and the stalled binary is not necessarily the same one twice. Re-running usually succeeds; please add the run URL to #2000 rather than only re-running."
            ;;
        execution-phase)
            echo "::error title=nextest stalled DURING test execution::The run announced its test count and then exceeded ${LIMIT}, so a test body is stuck. This is NOT the #2000 shape. The last \`PASS\`/\`SLOW\` lines above bound which test it is; nextest's \`slow-timeout\` in \`.config/nextest.toml\` is the right instrument for this half."
            ;;
        unknown)
            echo "::error title=nextest stalled, phase undetermined::The run exceeded ${LIMIT}, but its captured log is missing or empty, so none of the pre-list/list-phase/execution-phase discriminators could be applied. Treat the phase as UNKNOWN — do not record this as a #2000 occurrence. The capture itself is the thing to check first (\`mktemp\`/\`tee\` under \${RUNNER_TEMP})."
            ;;
    esac
fi

exit "$status"
