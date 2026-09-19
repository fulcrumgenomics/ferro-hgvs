#!/usr/bin/env bash
# Runs the oracle suite locally with the `oracle` nextest profile, which carries the
# exclusions and arms the flags. CI runs the same profile over its own shard and partition.
# A local run is a superset of CI's slice: it does not subtract the sweeps or the proptests.
# Add -E 'not test(proptest)' if that is too slow. See docs/ORACLES.md, section
# "Running the oracles locally".
set -euo pipefail
case "${1:-}" in
  --print-selection) shift; exec cargo nextest list --features dev --profile oracle "$@" ;;
  *)                 exec cargo nextest run  --features dev --profile oracle "$@" ;;
esac
