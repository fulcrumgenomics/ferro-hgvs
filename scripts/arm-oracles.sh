#!/usr/bin/env bash
# Nextest setup script. Arms the oracles for the profile that bound it by writing the
# FERRO_ASSERT_* flags to $NEXTEST_ENV. `oracle-rerun` runs the rows the sequence oracle
# skips, so it gets the other three.
set -euo pipefail
flags=(IDEMPOTENT REPARSE IN_BOUNDS)
[[ "$NEXTEST_PROFILE" == oracle ]] && flags+=(SEQUENCE)
printf 'FERRO_ASSERT_%s=1\n' "${flags[@]}" >> "$NEXTEST_ENV"
