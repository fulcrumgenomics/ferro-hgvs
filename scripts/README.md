# Scripts

Utility scripts for generating and updating test fixtures used by ferro-hgvs.

## Scripts

| Script | Description | Output Fixture |
|--------|-------------|----------------|
| `extract_hgvs_spec_examples.py` | Scrapes HGVS examples from [hgvs-nomenclature.org](https://hgvs-nomenclature.org) | `tests/fixtures/grammar/hgvs_spec_examples.json` |
| `extract_mutalyzer_github.py` | Extracts test patterns from the [mutalyzer-hgvs-parser](https://github.com/mutalyzer/mutalyzer-hgvs-parser) repo | `tests/fixtures/grammar/mutalyzer_github.json` |
| `fetch_ncbi_variation.py` | Fetches HGVS/SPDI conversions from the [NCBI Variation Services API](https://api.ncbi.nlm.nih.gov/variation/v0/) | `tests/fixtures/validation/ncbi_variation.json` |
| `fetch_variantvalidator.py` | Fetches validated HGVS variants from the [VariantValidator API](https://rest.variantvalidator.org) | `tests/fixtures/validation/variantvalidator_api.json` |
| `benchmark_python.py` | Benchmarks ferro-hgvs Python bindings against biocommons/hgvs | N/A (prints results) |

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

# Regenerate NCBI fixture (rate-limited, takes several minutes)
python scripts/fetch_ncbi_variation.py

# Regenerate VariantValidator fixture (rate-limited, takes several minutes)
python scripts/fetch_variantvalidator.py

# Run Python binding benchmarks
python scripts/benchmark_python.py --iterations 10000
```
