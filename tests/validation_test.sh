#!/bin/bash
# Quick validation test for ferro reference optimization.
#
# This script validates that the prepare/check commands work correctly
# by testing with a pre-existing ferro reference (if available) or
# by checking that the commands parse correctly.
#
# For a full end-to-end test that downloads data, use:
#   ./tests/validation_test.sh --full
#
# Usage: ./tests/validation_test.sh [--cleanup] [--full]
#
# Options:
#   --cleanup  Remove test directory after successful run
#   --full     Run full test with downloads (takes several minutes)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
TEST_DIR=$(mktemp -d)
CLEANUP=false
FULL_TEST=false

# Parse arguments
for arg in "$@"; do
    case $arg in
        --cleanup)
            CLEANUP=true
            ;;
        --full)
            FULL_TEST=true
            ;;
    esac
done

# Ensure ferro-benchmark is built
if [[ ! -f "$PROJECT_DIR/target/release/ferro-benchmark" ]]; then
    echo "Building ferro-benchmark..."
    (cd "$PROJECT_DIR" && cargo build --release --features benchmark)
fi

FERRO_BENCHMARK="$PROJECT_DIR/target/release/ferro-benchmark"

echo "============================================================"
echo "  Ferro Reference Validation Test"
echo "============================================================"
echo ""
echo "Test directory: $TEST_DIR"
echo "Mode: $(if $FULL_TEST; then echo 'FULL (with downloads)'; else echo 'QUICK (CLI validation only)'; fi)"
echo ""

# Quick mode: just validate CLI commands work
if [[ "$FULL_TEST" != "true" ]]; then
    echo "=== Quick Validation: CLI Commands ==="
    echo ""

    ERRORS=0

    # Test that prepare commands show help
    if "$FERRO_BENCHMARK" prepare ferro --help > /dev/null 2>&1; then
        echo "  [PASS] prepare ferro --help works"
    else
        echo "  [FAIL] prepare ferro --help failed"
        ERRORS=$((ERRORS + 1))
    fi

    if "$FERRO_BENCHMARK" prepare all --help > /dev/null 2>&1; then
        echo "  [PASS] prepare all --help works"
    else
        echo "  [FAIL] prepare all --help failed"
        ERRORS=$((ERRORS + 1))
    fi

    if "$FERRO_BENCHMARK" check ferro --help > /dev/null 2>&1; then
        echo "  [PASS] check ferro --help works"
    else
        echo "  [FAIL] check ferro --help failed"
        ERRORS=$((ERRORS + 1))
    fi

    if "$FERRO_BENCHMARK" check all --help > /dev/null 2>&1; then
        echo "  [PASS] check all --help works"
    else
        echo "  [FAIL] check all --help failed"
        ERRORS=$((ERRORS + 1))
    fi

    # Check if existing ferro reference is available
    if [[ -d "$PROJECT_DIR/data/ferro" && -f "$PROJECT_DIR/data/ferro/manifest.json" ]]; then
        echo ""
        echo "=== Found existing ferro reference, running check ==="
        if "$FERRO_BENCHMARK" check ferro --reference "$PROJECT_DIR/data/ferro" 2>&1 | grep -q "Ferro is ready"; then
            echo "  [PASS] check ferro passed on existing reference"
        else
            echo "  [WARN] check ferro did not pass (reference may be incomplete)"
        fi
    fi

    echo ""
    echo "============================================================"
    echo "  Summary (Quick Mode)"
    echo "============================================================"
    echo ""

    if [[ $ERRORS -eq 0 ]]; then
        echo "  All CLI validations PASSED"
        echo ""
        echo "  For full end-to-end test with downloads, run:"
        echo "    ./tests/validation_test.sh --full"
        if [[ "$CLEANUP" == "true" ]]; then
            rm -rf "$TEST_DIR"
        fi
        exit 0
    else
        echo "  $ERRORS validation(s) FAILED"
        exit 1
    fi
fi

# Full mode: download and prepare everything
echo "=== Step 1: Prepare Ferro Reference ==="
echo "(This downloads ~1GB of RefSeq data and takes several minutes)"
echo ""

"$FERRO_BENCHMARK" prepare ferro \
    --output-dir "$TEST_DIR/ferro" \
    --genome none \
    --no-refseqgene \
    --no-lrg \
    2>&1 | tee "$TEST_DIR/ferro.log"

# Verify ferro created expected files
echo ""
echo "=== Verifying Ferro Output ==="
ERRORS=0

if [[ -f "$TEST_DIR/ferro/manifest.json" ]]; then
    echo "  [PASS] manifest.json exists"
else
    echo "  [FAIL] manifest.json missing"
    ERRORS=$((ERRORS + 1))
fi

if [[ -d "$TEST_DIR/ferro/transcripts" ]]; then
    echo "  [PASS] transcripts directory exists"
else
    echo "  [FAIL] transcripts directory missing"
    ERRORS=$((ERRORS + 1))
fi

if [[ -d "$TEST_DIR/ferro/cdot" ]]; then
    echo "  [PASS] cdot directory exists"
else
    echo "  [FAIL] cdot directory missing"
    ERRORS=$((ERRORS + 1))
fi

echo ""
echo "=== Step 2: Prepare Mutalyzer Cache ==="
echo ""

# Note: Using --no-proteins since we don't have ClinVar
"$FERRO_BENCHMARK" prepare mutalyzer \
    --ferro-reference "$TEST_DIR/ferro" \
    --output-dir "$TEST_DIR/mutalyzer" \
    --no-proteins \
    2>&1 | tee "$TEST_DIR/mutalyzer.log"

# Verify mutalyzer output
echo ""
echo "=== Verifying Mutalyzer Output ==="

if [[ -f "$TEST_DIR/mutalyzer/mutalyzer_settings.conf" ]]; then
    echo "  [PASS] mutalyzer_settings.conf exists"
else
    echo "  [FAIL] mutalyzer_settings.conf missing"
    ERRORS=$((ERRORS + 1))
fi

# Count cache entries
CACHE_COUNT=$(find "$TEST_DIR/mutalyzer" -name "*.sequence" 2>/dev/null | wc -l | tr -d ' ')
if [[ "$CACHE_COUNT" -gt 0 ]]; then
    echo "  [PASS] Cache contains $CACHE_COUNT sequence entries"
else
    echo "  [WARN] Cache contains no sequence entries (may be expected for minimal test)"
fi

echo ""
echo "=== Step 3: Check Ferro Reference ==="
echo ""

"$FERRO_BENCHMARK" check ferro --reference "$TEST_DIR/ferro" 2>&1 | tee "$TEST_DIR/check_ferro.log"

if grep -q "Ferro is ready" "$TEST_DIR/check_ferro.log"; then
    echo "  [PASS] Ferro check passed"
else
    echo "  [FAIL] Ferro check failed"
    ERRORS=$((ERRORS + 1))
fi

echo ""
echo "============================================================"
echo "  Summary"
echo "============================================================"
echo ""

if [[ $ERRORS -eq 0 ]]; then
    echo "  All validations PASSED"
    echo ""
    if [[ "$CLEANUP" == "true" ]]; then
        echo "  Cleaning up test directory..."
        rm -rf "$TEST_DIR"
    else
        echo "  Test directory preserved at: $TEST_DIR"
    fi
    exit 0
else
    echo "  $ERRORS validation(s) FAILED"
    echo ""
    echo "  Test directory preserved at: $TEST_DIR"
    echo "  Check logs for details:"
    echo "    - $TEST_DIR/ferro.log"
    echo "    - $TEST_DIR/mutalyzer.log"
    echo "    - $TEST_DIR/check_ferro.log"
    exit 1
fi
