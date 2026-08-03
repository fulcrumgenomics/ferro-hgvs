# Scripts

Utility scripts for generating and updating test fixtures used by ferro-hgvs.

## Scripts

| Script | Description | Output Fixture |
|--------|-------------|----------------|
| `extract_hgvs_spec_examples.py` | Scrapes HGVS examples from [hgvs-nomenclature.org](https://hgvs-nomenclature.org) | `tests/fixtures/grammar/hgvs_spec_examples.json` |
| `extract_mutalyzer_github.py` | Extracts test patterns from the [mutalyzer-hgvs-parser](https://github.com/mutalyzer/mutalyzer-hgvs-parser) repo | `tests/fixtures/grammar/mutalyzer_github.json` |
| `fetch_ncbi_variation.py` | Ad-hoc live fetch of HGVS/SPDI conversions from the [NCBI Variation Services API](https://api.ncbi.nlm.nih.gov/variation/v0/). **Does not regenerate the committed fixture** (see note below). | `tests/fixtures/validation/ncbi_variation.live.json` (live; not committed) |
| `fetch_variantvalidator.py` | Fetches validated HGVS variants from the [VariantValidator API](https://rest.variantvalidator.org) | `tests/fixtures/validation/variantvalidator_api.json` |
| `benchmark_python.py` | Benchmarks ferro-hgvs Python bindings against biocommons/hgvs | N/A (prints results) |
| `run_conformance_axis.sh` | Runs the manifest-backed conformance axis (`axis_*` tests) with a validated `FERRO_MANIFEST` | N/A (runs tests) |

## Requirements

The fetch/extract scripts require `requests` and `beautifulsoup4`:

```bash
pip install requests beautifulsoup4
```

The benchmark script requires the ferro-hgvs Python package:

```bash
maturin develop --features python
```

## Usage

All scripts support `--help` for full usage information.  Typical usage:

```bash
# Regenerate HGVS spec examples fixture
python scripts/extract_hgvs_spec_examples.py

# Regenerate mutalyzer fixture (requires cloned repo)
git clone https://github.com/mutalyzer/mutalyzer-hgvs-parser external-repos/mutalyzer-hgvs-parser
python scripts/extract_mutalyzer_github.py

# Ad-hoc live NCBI fetch (rate-limited, takes several minutes).
# NOTE: tests/fixtures/validation/ncbi_variation.json is a curated, hand-authored
# offline fixture that tests/it/spdi_tests.rs parses and asserts semantic invariants
# against (not a byte-level file comparison). This script
# does NOT regenerate it — by default it writes the live response to
# tests/fixtures/validation/ncbi_variation.live.json instead. To inspect live data,
# point --output at a scratch path; do not overwrite the curated fixture.
python scripts/fetch_ncbi_variation.py --output /tmp/ncbi_variation.live.json

# Regenerate VariantValidator fixture (rate-limited, takes several minutes)
python scripts/fetch_variantvalidator.py

# Run Python binding benchmarks
python scripts/benchmark_python.py --iterations 10000

# Run the manifest-backed conformance axis (axis_ tests). Fails loudly
# instead of silently passing if FERRO_MANIFEST is unset or missing --
# see the script's header comment for why that matters.
FERRO_MANIFEST=/path/to/manifest.json scripts/run_conformance_axis.sh
# or:
scripts/run_conformance_axis.sh /path/to/manifest.json
```

### Conformance axis (`axis_*` tests)

The conformance corpora (`tests/it/mutalyzer_normalize_tests.rs`,
`tests/it/biocommons_normalize_tests.rs`,
`tests/it/hgvs_rs_projection_tests.rs`) assert against a real prepared
reference and are gated on the `FERRO_MANIFEST` environment variable: when it
is unset, or points at a path that does not exist, each test's manifest-lookup
helper returns early *without asserting anything*, and `cargo nextest` reports
that as **PASSED**, not skipped. Running
`cargo nextest run --features dev -E 'test(axis_)'` directly, without
`FERRO_MANIFEST` set, is therefore a vacuous all-green run for exactly those
tests, and looks identical to a real conformance pass.

`scripts/run_conformance_axis.sh` closes that hole by validating the manifest
path *before* invoking nextest, hard-failing if it is unset or missing. Build
a manifest with `ferro prepare` (see `ferro prepare --help`); the manifest
itself is never committed and this script never hardcodes a machine-local
path.

**What the printed test count does and does not prove.** `test(axis_)` is a
substring match, so it selects ~179 tests — but only ~26 of them, the three
corpora above, are manifest-gated. The rest (`cli::project::tests`,
`hgvs::variant::tests`, `issue_337_cds_start_clamp`, `issue_1086_*`, …) merely
contain `axis_` in their names and run identically with or without a manifest.
They are harmless to run, but a non-zero count on its own is **not** evidence
that conformance ran; it only shows the filter matched something. The
pre-flight manifest check is the load-bearing guard.
