#!/usr/bin/env bash
# Runs the seam-oracle suite locally the way `ci.yml`'s `test-oracle` job runs
# it: the oracle flags that job sets, over that job's selection.
#
# WHY THIS EXISTS. The obvious local command --
#
#     FERRO_ASSERT_IDEMPOTENT=1 cargo nextest run --features dev
#
# -- CANNOT PASS on `main`, and has not been able to for as long as the spec
# corpus has existed. It fails 7 tests, and all 7 are in modules `ci.yml` names
# in `ORACLE_EXCLUDE` precisely so that the armed job never runs them:
#
#   * defect_non_idempotent_outputs (4 tests) and
#     spec_corpus_regressions::an_insertion_at_the_cds_end_is_not_a_fixed_point
#     ASSERT the non-idempotency the oracle PANICS on. `c.*1delinsCTT` ->
#     `c.72_*1insCT` -> `c.72delinsCCT` is a pinned defect; a test that pins it
#     and an oracle that aborts on it cannot both run.
#   * spec_conformance_axis's two censuses COUNT it. The corpus wraps
#     normalization in `catch_unwind`, so a panicking row is filed `declined`
#     and its output never reaches the family's output set -- which does not
#     redden the census, it FLATTERS it. Measured on `main` @ 1dd8148d, both
#     directions, armed against the committed pins:
#
#       3': declined 0 -> 4, non_idempotent_outputs 4 -> 0,
#           converged 9141 -> 9145, split_two 2440 -> 2436
#       5': declined 0 -> 4, non_idempotent_outputs 4 -> 0,
#           converged 8944 -> 8948, split_two 2706 -> 2702
#
#     Four families read as converged in each direction that do not converge.
#
# So a bare armed run is not a coverage gap to be closed by running it anyway;
# it is a measurement taken with two instruments that destroy each other. This
# script applies the same exclusion CI does, so a local armed run is a signal
# rather than a known-red wall.
#
# THE SELECTION AND THE FLAGS ARE READ FROM `ci.yml`, NEVER COPIED. That means
# the WHOLE `-E` expression, not just `ORACLE_EXCLUDE`: `test-oracle` also
# negates `test(proptest)` and `SWEEP_FILTER`, so a run that dropped either
# would execute tests that job does not while claiming to mirror it. A second
# copy of any of it would drift, and a drifted exclusion here fails in the
# flattering direction -- it would exclude a module CI runs, so a local run
# would go green on a defect CI is red on.
# `tests/it/oracle_exclude_invariant.rs` invokes `--print-selection` below and
# compares this extraction against one derived independently in Rust, so a
# `ci.yml` restructure that breaks the awk fails loudly instead of yielding an
# empty filter.
#
# Usage:
#   scripts/run_oracle_suite.sh                    # run it
#   scripts/run_oracle_suite.sh --print-selection  # print what it would run, and stop
#   scripts/run_oracle_suite.sh -E 'test(foo)'     # extra args go through to nextest
#
# The denoted-sequence oracle (`FERRO_ASSERT_SEQUENCE`) is deliberately NOT
# among the flags this mirrors, because `test-oracle` does not set it. Nothing
# has to be done here to keep that true: `oracle_flags` below READS the flag set
# out of that job, so this script arms it on the day the job does and not
# before.
#
# The figure this header used to quote -- "5 further failures, all in
# `spec_corpus_regressions`" -- was about `ORACLE_EXCLUDE`'s modules, which this
# selection does not run, so it never described what it claimed to. Measured
# over THIS selection instead, on `origin/main` @ 674e9c8b: 5 failures, none of
# them in `spec_corpus_regressions` -- 3 in `issue_1487_canonical_window_overflow`
# (an `i64` overflow at `src/convert/mapper.rs:114`, issue #1690) and 2 in
# `stranded_identity_member` (a real fire on a module that PINS a defect). The
# blocking rows are no longer #1618/#1619, both of which are closed and green;
# see `test-oracle`'s own comment in `ci.yml` for the full triage.
#
# RE-MEASURED on `c9207d7e` once #1690 closed (#1990), same selection and flags:
# `10904 tests run: 10902 passed, 2 failed, 306 skipped`. The 3
# `issue_1487_canonical_window_overflow` rows are GONE and nothing new fired, so
# the remaining blocker is `stranded_identity_member` alone -- #1690 is closed
# and is no longer one.
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CI_YML="$REPO_ROOT/.github/workflows/ci.yml"

if [[ ! -f "$CI_YML" ]]; then
    echo "error: cannot find .github/workflows/ci.yml at '$CI_YML'" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

# A `KEY: >-` folded scalar from `ci.yml`'s top-level `env:`, folded back onto
# one line.
#
# Reading only the first line would silently exempt whichever module happens to
# sit on a continuation line, which is the same blind spot
# `oracle_exclude_invariant.rs` exists to close -- and for `ORACLE_EXCLUDE` it
# fails flatteringly, since an unexcluded module there is one this script would
# run armed.
#
# Parameterised by key rather than written twice: `test-oracle`'s selection
# negates `SWEEP_FILTER` as well as `ORACLE_EXCLUDE`, and two copies of one
# folding rule is the drift this script exists to avoid.
#
# `$1` is interpolated into an awk regex, so it must be a bare YAML key with no
# regex metacharacters. Both callers pass a literal; a key containing `.` or `[`
# would need escaping first.
folded_scalar() {
    awk -v key="$1" '
        $0 ~ "^[[:space:]]*" key ":[[:space:]]*[>|]" { found = 1; next }
        found {
            if ($0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/) { exit }
            # A line at or left of the key indent ends the block scalar.
            if ($0 !~ /^[[:space:]][[:space:]][[:space:]]/) { exit }
            gsub(/^[[:space:]]+|[[:space:]]+$/, "")
            printf "%s%s", (n++ ? " " : ""), $0
        }
    ' "$CI_YML"
}

# The `-E` expression `test-oracle`'s own step hands to nextest, verbatim and
# still holding its `$SWEEP_FILTER` / `$ORACLE_EXCLUDE` references.
#
# Read whole rather than reassembled from its parts. An earlier revision of this
# script negated only `ORACLE_EXCLUDE`, so a local run also executed the proptest
# modules and the three exhaustive sweeps -- tests `test-oracle` does not run --
# while claiming to mirror that job. Rebuilding the expression here from a
# hardcoded `not test(proptest) and not (...) and not (...)` would fix today's
# drift by introducing tomorrow's: the shape of CI's filter would then live in
# two files. So the shape is read too, and only the variable references are
# expanded.
oracle_selection_template() {
    awk '
        /^  test-oracle:/ { in_job = 1; next }
        # A non-indented, non-comment line at job-key indent ends the job.
        in_job && /^  [^ #]/ { exit }
        in_job {
            if (match($0, /-E "[^"]*"/)) {
                print substr($0, RSTART + 4, RLENGTH - 5)
                exit
            }
        }
    ' "$CI_YML"
}

# The `FERRO_ASSERT_*` keys the `test-oracle` job's own step sets, in file
# order.
#
# Scoped to the window between that step's `name:` and its `run:`, because the
# comment block inside that window MENTIONS `FERRO_ASSERT_SEQUENCE` in prose to
# explain why it is absent. Comment lines are skipped, so the prose mention is
# not mistaken for a setting -- if it were, this script would arm an oracle CI
# does not and go red on the two rows that comment is about.
oracle_flags() {
    awk '
        /^[[:space:]]*-[[:space:]]*name:[[:space:]]*Run Rust tests with the normalization self-checks/ {
            found = 1; next
        }
        found && /^[[:space:]]*run:/ { exit }
        found && /^[[:space:]]*#/ { next }
        found && /^[[:space:]]*FERRO_ASSERT_[A-Z_]+:/ {
            sub(/:.*/, "")
            gsub(/^[[:space:]]+|[[:space:]]+$/, "")
            print
        }
    ' "$CI_YML"
}

EXCLUDE="$(folded_scalar ORACLE_EXCLUDE)"
SWEEPS="$(folded_scalar SWEEP_FILTER)"
CENSUSES="$(folded_scalar CENSUS_FILTER)"
TEMPLATE="$(oracle_selection_template)"
# A read loop rather than `mapfile`, which is bash 4+: stock macOS still ships
# bash 3.2 as /bin/bash, and a script that dies on `mapfile: command not found`
# would send the reader looking at their shell instead of at the oracle.
FLAGS=()
while IFS= read -r flag; do
    [[ -n "$flag" ]] && FLAGS+=("$flag")
done < <(oracle_flags)

# Refuse a vacuous extraction rather than running something weaker than CI.
# An empty exclusion would run the 7 known-red tests; an empty flag list would
# run the whole suite with no oracle armed at all and report it as an oracle
# pass, which is the worse of the two.
if [[ -z "$EXCLUDE" ]]; then
    echo "error: could not read ORACLE_EXCLUDE from $CI_YML." >&2
    echo "  Its formatting changed; fix the awk in this script rather than inlining a copy." >&2
    exit 1
fi
if [[ "${#FLAGS[@]}" -eq 0 ]]; then
    echo "error: could not read any FERRO_ASSERT_* flag from test-oracle's step in $CI_YML." >&2
    echo "  Its formatting changed; fix the awk in this script rather than inlining a copy." >&2
    exit 1
fi
if [[ -z "$SWEEPS" ]]; then
    echo "error: could not read SWEEP_FILTER from $CI_YML." >&2
    echo "  Its formatting changed; fix the awk in this script rather than inlining a copy." >&2
    exit 1
fi
# `CENSUS_FILTER`'s modules moved to the `censuses` job, which runs them on the
# optimized archive; `test-oracle` negates them. An empty read here would put
# them back into a local armed run -- not a known-red wall like an empty
# `ORACLE_EXCLUDE`, but ~9 minutes of debug-profile census this script is not
# meant to run, which reads as a hang rather than as a misconfiguration.
if [[ -z "$CENSUSES" ]]; then
    echo "error: could not read CENSUS_FILTER from $CI_YML." >&2
    echo "  Its formatting changed; fix the awk in this script rather than inlining a copy." >&2
    exit 1
fi
if [[ -z "$TEMPLATE" ]]; then
    echo "error: could not read test-oracle's -E selection from $CI_YML." >&2
    echo "  Its formatting changed; fix the awk in this script rather than inlining a copy." >&2
    exit 1
fi

# Expand the two references CI's expression carries. Substitution rather than
# `eval`, so nothing in `ci.yml` can execute here.
SELECTION="${TEMPLATE//\$SWEEP_FILTER/$SWEEPS}"
SELECTION="${SELECTION//\$ORACLE_EXCLUDE/$EXCLUDE}"
SELECTION="${SELECTION//\$CENSUS_FILTER/$CENSUSES}"

# Refuse a selection still naming a variable. A new `$FOO` in that expression
# would otherwise reach nextest literally -- which either errors obscurely or,
# worse, parses as a test-name substring and quietly narrows the run.
if [[ "$SELECTION" == *'$'* ]]; then
    echo "error: test-oracle's -E selection references a variable this script does not expand:" >&2
    echo "  $SELECTION" >&2
    echo "  Teach this script to read it from $CI_YML rather than inlining its value." >&2
    exit 1
fi

if [[ "${1:-}" == "--print-selection" ]]; then
    printf 'ORACLE_EXCLUDE=%s\n' "$EXCLUDE"
    printf 'SWEEP_FILTER=%s\n' "$SWEEPS"
    printf 'CENSUS_FILTER=%s\n' "$CENSUSES"
    printf 'SELECTION=%s\n' "$SELECTION"
    for flag in "${FLAGS[@]}"; do printf 'FLAG=%s\n' "$flag"; done
    exit 0
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

for flag in "${FLAGS[@]}"; do export "$flag=1"; done

echo "Running the seam-oracle suite with: ${FLAGS[*]}"
echo "Selecting (as ci.yml's test-oracle job does): $SELECTION"

OUTPUT="$(mktemp)"
trap 'rm -f "$OUTPUT"' EXIT

set +e
cargo nextest run --features dev -E "$SELECTION" "$@" 2>&1 | tee "$OUTPUT"
STATUS="${PIPESTATUS[0]}"
set -e

# nextest writes "1 test run" (singular) and "N tests run" (plural), so the
# optional `s` is load-bearing -- without it a legitimate single-test run parses
# as no run at all. Same rule as `run_conformance_axis.sh`.
RAN="$(grep -Eo '[0-9]+ tests? run' "$OUTPUT" | tail -1 | grep -Eo '^[0-9]+' || true)"

if [[ "$STATUS" -ne 0 ]]; then
    echo "error: the seam-oracle suite FAILED (nextest exit ${STATUS})." >&2
    exit "$STATUS"
fi

# A green run over zero tests is the failure this repo keeps meeting, so the
# denominator is asserted rather than assumed.
if [[ -z "$RAN" || "$RAN" -eq 0 ]]; then
    echo "error: no tests ran -- this would be a vacuous oracle run." >&2
    exit 1
fi

echo "Seam-oracle suite: ${RAN} test(s) ran armed with ${FLAGS[*]}."
