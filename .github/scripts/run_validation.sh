#!/usr/bin/env bash
# Run one external-validation suite so that every way it can go wrong is loud.
#
# Sourced by `.github/workflows/external-validation.yml`. Three failure modes
# are in scope, and the workflow had been silently surviving all three (#1404):
#
#   1. the tests fail            — propagated, because `pipefail` is set and the
#                                  step no longer carries `continue-on-error`
#   2. the target does not exist — `cargo test` exits non-zero, which #725's
#                                  binary consolidation made the live case
#   3. the filter selects zero   — exits 0 with `running 0 tests`, so it is
#      tests                       checked explicitly here
#
# The third is the one worth the script. A suite that runs nothing looks exactly
# like a suite that passes, and it is *faster*, so it reads as an improvement.
#
# Usage: run_validation <module> <output-file> [extra cargo-test args...]
run_validation() {
    local module="$1"
    local output="$2"
    shift 2

    set -euo pipefail

    # `--test it`: every integration test lives in the single `it` binary since
    # #725. `$module` is a filter on the test *name*, which is
    # `<module>::<test>`, so it selects that module's tests and nothing else.
    cargo test --test it "$module" -- --nocapture "$@" 2>&1 | tee "$output"

    # `cargo test` prints one `running N tests` line per binary. Sum them, so a
    # filter that matched nothing is an error rather than a fast pass.
    local ran
    ran=$(awk '/^running [0-9]+ test/ { total += $2 } END { print total + 0 }' "$output")

    if [[ "$ran" -eq 0 ]]; then
        echo "::error::${module}: the filter selected no tests, so this suite validated nothing." >&2
        echo "Check that the module still exists under tests/it/ and that any" >&2
        echo "--ignored flag still matches how its tests are marked." >&2
        return 1
    fi

    echo "${module}: ran ${ran} test(s)."
}
