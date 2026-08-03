#!/usr/bin/env bash
# Runs the manifest-backed HGVS conformance axis (the `axis_*` test family
# in tests/it/{mutalyzer_normalize,biocommons_normalize,hgvs_rs_projection,...}
# _tests.rs) against a real ferro-prepared reference.
#
# Why this script exists: the manifest-backed conformance tests are gated at
# runtime on `FERRO_MANIFEST` -- when it's unset, or points at a path that
# doesn't exist, the test's own manifest-lookup helper returns `None` and the
# test body returns early without asserting anything. `cargo nextest` reports
# that early return as PASSED, not skipped. So `cargo nextest run --features
# dev -E 'test(axis_)'` run bare is a vacuous, all-green no-op for exactly the
# tests you wanted -- it can look identical to a real conformance pass.
# Confirmed directly: with FERRO_MANIFEST unset the filter reports
# "179 tests run: 179 passed".
#
# This script closes that hole by validating the manifest *before* invoking
# nextest at all, so a missing/unset manifest is a hard failure here -- never
# a silent pass buried inside the test suite. That pre-flight check is what
# does the real work; the test count reported at the end is a weaker guard
# that only proves the filter matched something (see the note below).
#
# Note on the filter: `test(axis_)` is a substring match, so it selects far
# more than the conformance corpora -- ~179 tests, of which only the ~26 in
# mutalyzer_normalize_tests / biocommons_normalize_tests /
# hgvs_rs_projection_tests are manifest-gated. The rest (cli::project::tests,
# hgvs::variant::tests, issue_337_cds_start_clamp, issue_1086_*, ...) merely
# have "axis_" in their names and run identically with or without a manifest.
# They are harmless to run, but they mean a non-zero test count on its own is
# NOT evidence that conformance ran. Do not weaken the pre-flight check on the
# strength of the count.
#
# Usage:
#   FERRO_MANIFEST=/path/to/manifest.json scripts/run_conformance_axis.sh
#   scripts/run_conformance_axis.sh /path/to/manifest.json
#
# The manifest is always operator-supplied -- this script never hardcodes a
# machine-local path. Build one with `ferro prepare` (see `ferro prepare
# --help`, or README.md's Reference Data section).
set -euo pipefail

MANIFEST="${1:-${FERRO_MANIFEST:-}}"

if [[ -z "$MANIFEST" ]]; then
    echo "error: no conformance manifest supplied." >&2
    echo "  Set FERRO_MANIFEST=/path/to/manifest.json, or pass the path as \$1 to this script." >&2
    echo "  Build one with 'ferro prepare' (see 'ferro prepare --help')." >&2
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "error: manifest not found at '$MANIFEST'" >&2
    echo "  FERRO_MANIFEST (or \$1) must point at an existing manifest.json produced by 'ferro prepare'." >&2
    exit 1
fi

export FERRO_MANIFEST="$MANIFEST"
echo "Running conformance axis (axis_ tests) with FERRO_MANIFEST=${FERRO_MANIFEST}"

OUTPUT="$(mktemp)"
trap 'rm -f "$OUTPUT"' EXIT

# Tee so the operator sees live nextest output, then re-derive the summary
# from the captured copy to confirm a non-zero number of tests actually ran
# -- the whole point of this script is to make a vacuous run impossible to
# miss.
set +e
cargo nextest run --features dev -E 'test(axis_)' 2>&1 | tee "$OUTPUT"
STATUS="${PIPESTATUS[0]}"
set -e

# nextest writes "1 test run" (singular) and "N tests run" (plural), so the
# optional `s` is load-bearing: without it a legitimate single-test run parses
# as no run at all and is misreported below as vacuous.
RAN="$(grep -Eo '[0-9]+ tests? run' "$OUTPUT" | tail -1 | grep -Eo '^[0-9]+' || true)"

# A non-zero nextest status is reported as itself. Only claim "vacuous" when
# the run actually succeeded and still matched nothing -- otherwise a build
# break (no summary line at all, so RAN is empty) gets misattributed.
if [[ "$STATUS" -ne 0 ]]; then
    echo "error: conformance axis FAILED (nextest exit ${STATUS}) against ${FERRO_MANIFEST}." >&2
    exit "$STATUS"
fi

if [[ -z "$RAN" || "$RAN" -eq 0 ]]; then
    echo "error: no axis_ tests ran -- this would be a vacuous conformance run." >&2
    exit 1
fi

echo "Conformance axis: ${RAN} axis_ test(s) ran against ${FERRO_MANIFEST}."
