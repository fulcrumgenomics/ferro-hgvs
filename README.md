[![CI](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml)
[![Nightly reference-aware tests](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml/badge.svg?branch=main)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml)
[![Codecov](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs)
[![Crates.io](https://img.shields.io/crates/v/ferro-hgvs.svg)](https://crates.io/crates/ferro-hgvs)
[![PyPI](https://img.shields.io/pypi/v/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Python versions](https://img.shields.io/pypi/pyversions/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Documentation](https://readthedocs.org/projects/ferro-hgvs/badge/?version=stable)](https://ferro-hgvs.readthedocs.io/en/stable/)
[![API docs](https://docs.rs/ferro-hgvs/badge.svg)](https://docs.rs/ferro-hgvs)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ferro-hgvs/README.html)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18703103-blue.svg)](https://doi.org/10.5281/zenodo.18703103)

# ferro-hgvs

A high-performance HGVS variant nomenclature parser and normalizer written in Rust.

**WARNING: ALPHA SOFTWARE - USE AT YOUR OWN RISK**

This software is currently in **ALPHA**. While we have extensively tested it
across a wide variety of HGVS patterns, **no guarantees are made** regarding
correctness or stability.

<p>
<a href="https://fulcrumgenomics.com">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-light.svg">
    <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-light.svg" height="100">
  </picture>
</a>
</p>

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Features

- **Full HGVS Parsing**: All coordinate systems (g/c/n/r/p/m/o) and edit types
- **Variant Normalization**: 3' shifting per HGVS specification
- **High Performance**: ~5M variants/sec single-threaded parsing (>12M/s parallel), zero-copy with nom
- **Type-Safe**: Leverages Rust's type system for correctness

## Installation

### Python

```bash
pip install ferro-hgvs
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, Apple Silicon), and Windows (x86_64) on Python 3.10+.

### Rust

Add to your `Cargo.toml`:

```toml
[dependencies]
ferro-hgvs = "0.1"
```

Or install the CLI:

```bash
cargo install ferro-hgvs
```

## Quick Start

### CLI

```bash
# Parse a variant
ferro parse "NM_000088.3:c.459A>G"

# Parse from file
ferro parse -i variants.txt -f json

# Prepare reference data (downloads RefSeq, genome, cdot — RefSeq-only by default)
ferro prepare --output-dir ferro-reference

# Verify reference data is ready
ferro check --reference ferro-reference

# (Optional) pre-build the on-disk cdot cache as a setup step, so the one-time
# cache build doesn't slow the start of a real (or timed/benchmarked) run.
ferro check --reference ferro-reference --build-cache

# Normalize with reference
ferro normalize "NM_000088.3:c.459del" --reference ferro-reference/
```

> **Read the warnings.** Normalization sometimes *repairs* a description in a way
> the normalized string does not record — separately reported cis members merged
> into one `delins` (`MEMBERS_COALESCED_FROM_REPORTED_FORM`), a `ins[100_110]`
> reference-range payload replaced by the bases it denotes
> (`INSERTED_SEQUENCE_EXPANDED`), a stated reference base that contradicted the
> reference and was accepted anyway (`REFSEQ_MISMATCH`). Those are reported as
> `warning[CODE]: message` on **stderr** (and in the `warnings` array under
> `--format json`, the `detail` column under `--format tsv`), so a pipeline
> reading only stdout will not see them. `--error-mode strict` is not a
> substitute: it rejects a specific ladder of conditions and reports the rest
> exactly as lenient does.

> **Throughput tip:** when normalizing many variants, feed them **sorted by transcript accession (or by genomic position)**. ferro caches each resolved transcript, so consecutive variants on the same transcript skip the (dominant) cost of re-reading and re-building it from the reference. Sorted input keeps the relevant transcripts resident in the cache and is markedly faster on large batches — see [Performance Comparison](#performance-comparison).

### Optional reference data

A bare `ferro prepare` builds a **RefSeq-only** reference (accessions `NM_`/`NR_`/`NP_`/`NG_`). Two opt-in flags provision additional data — pass them at prepare time; they are what a fully-provisioned ("blessed") reference is built with:

```bash
# Add Ensembl support (accessions ENST/ENSG/ENSP). Downloads the Ensembl cdot
# metadata and cDNA FASTAs (~1 GB+); off by default. Without it, an ENST/ENSG/ENSP
# input reports "Reference not found" and the message points back at this flag.
ferro prepare --output-dir ferro-reference --ensembl

# Derive version-independent NG_ placements and the NG_→transcript-version map
# (ng_hosted_transcripts) for a curated list of RefSeqGene accessions. Required to
# resolve legacy gene-symbol selectors (NG_(GENE):c.…) and bare-NG_ hosted lookups.
ferro prepare --output-dir ferro-reference \
  --derive-ng-placements path/to/ng_accessions.txt

# A fully-provisioned reference combines both in one run:
ferro prepare --output-dir ferro-reference --ensembl \
  --derive-ng-placements path/to/ng_accessions.txt
```

Both flags are incremental: re-running `ferro prepare` over an existing reference adds the requested data and preserves already-provisioned artifacts.

### Library

```rust
use ferro_hgvs::{parse_hgvs, HgvsVariant};

fn main() -> Result<(), ferro_hgvs::FerroError> {
    let variant = parse_hgvs("NM_000088.3:c.459A>G")?;

    match &variant {
        HgvsVariant::Cds(v) => println!("CDS variant: {}", v),
        HgvsVariant::Genome(v) => println!("Genomic variant: {}", v),
        _ => println!("Other: {}", variant),
    }

    Ok(())
}
```

### Python

```python
import ferro_hgvs

# Parse a variant
variant = ferro_hgvs.parse("NM_000088.3:c.459A>G")
print(variant.variant_type)  # "coding"
print(variant.reference)     # "NM_000088.3"
print(str(variant))          # "NM_000088.3:c.459A>G"

# Normalize with reference data
normalizer = ferro_hgvs.Normalizer(reference_json="ferro-reference/cdot.json")
normalized = normalizer.normalize("NM_000088.3:c.459del")

# `normalize` returns only the string, so it cannot tell you that normalization
# repaired something. `normalize_with_warnings` returns the same string plus the
# diagnostics — as a free function, or as a Normalizer method.
result = normalizer.normalize_with_warnings("NM_000088.3:c.459del")
print(str(result.result))                      # the same normalized string
print([(w.code, w.message) for w in result.warnings])
```

## Documentation

Full guides live in **[the documentation site](https://ferro-hgvs.readthedocs.io/)** (source under [`docs/src/`](docs/src/)):

- [Normalization rules](docs/src/reference/normalization-rules.md) — ferro's output contract: what a normalized description is allowed to be
- [Deriving a description from sequences](docs/src/guide/deriving-from-sequences.md) — turn a window of bases into one canonical description, no reference needed
- [Normalize variants](docs/src/guide/normalize-variants.md) · [Reference data](docs/src/guide/reference-data.md)
- [Project to another axis](docs/src/guide/project-variants.md) — re-express a variant on a transcript; the projected `c.`/`r.` axis applies the coding-axis rules a genomic axis cannot
- [Error handling](docs/src/guide/error-handling.md) — strict / lenient / silent modes and warning codes
- [Comparing normalization rules (`FERRO_PARTITION`)](docs/src/guide/comparing-rules.md)
- [Benchmarking](docs/src/guide/benchmarking.md)
- [CLI reference](docs/src/cli/overview.md) · [Supported HGVS syntax](docs/src/reference/hgvs-syntax.md)
- [Why ferro-hgvs? — tool comparison](docs/src/reference/comparison.md) — capability matrix, cross-tool benchmarks, and what `ferro prepare` builds
- [Interpreting the HGVS recommendations](docs/src/shadow-spec/index.md) — how ferro reads each clause


## Why ferro-hgvs?

ferro-hgvs provides the most comprehensive HGVS variant normalization across all pattern types, with performance orders of magnitude faster than alternatives. For the full capability matrix against mutalyzer / biocommons / hgvs-rs, the cross-tool parse/normalize benchmarks, and what `ferro prepare` builds, see [**Why ferro-hgvs? — tool comparison**](docs/src/reference/comparison.md).

## Benchmark: Reference Data & Tool Comparison

The main `ferro` binary includes commands to prepare reference data (`ferro prepare`) and check its status (`ferro check`). The `ferro-benchmark` tool (build with `--features benchmark`) extends this for tool comparison benchmarks.

| Command | Description |
|---------|-------------|
| `prepare <tool>` | Prepare reference data for a tool |
| `check <tool>` | Verify tool configuration and dependencies |
| `parse <tool>` | Parse HGVS patterns with specified tool |
| `normalize <tool>` | Normalize HGVS patterns with specified tool |
| `compare results` | Compare parse/normalize results between tools |
| `extract` | Extract patterns from ClinVar, VCFs, or create samples |
| `setup` | Set up UTA database, SeqRepo, and other services |
| `generate` | Generate summary reports and configs |
| `collate` | Aggregate sharded results |

### Quick Start

```bash
# Prepare ferro reference (main binary - no special features needed)
ferro prepare --output-dir data/ferro

# Check reference data
ferro check --reference data/ferro

# Normalize with ferro
ferro normalize -i patterns.txt --reference data/ferro

# For tool comparison, build with benchmark support
cargo build --release --features benchmark

# Prepare other tools (uses ferro reference for transcript data)
ferro-benchmark prepare mutalyzer --ferro-reference data/ferro --output-dir data/mutalyzer
ferro-benchmark prepare biocommons --seqrepo-dir data/seqrepo --uta-dump uta_20210129b.pgd.gz --ferro-reference data/ferro

# Compare results between tools
ferro-benchmark normalize mutalyzer -i patterns.txt -o mutalyzer.json --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf
ferro-benchmark compare results normalize ferro.json mutalyzer.json -o comparison.json
```

**Supported tools**: ferro-hgvs, mutalyzer, biocommons/hgvs, hgvs-rs

> **Note**: The `pixi.toml` and `pixi.lock` files in this repository define a [pixi](https://pixi.sh) environment for the Python-based external tools (mutalyzer, biocommons/hgvs, seqrepo) used in benchmarking. Run `pixi shell` to activate it.

See [docs/BENCHMARK_GUIDE.md](docs/BENCHMARK_GUIDE.md) for detailed usage.

## Development

```bash
cargo build
cargo test                        # default features
cargo clippy -- -D warnings
```

The commands above use the default feature set, and CI keeps them compiling (see
the `build` job). They do not cover the whole suite — the feature-gated tests and
the integration tree need `dev`, which is what CI runs and what you want before
opening a PR:

```bash
cargo nextest run --features dev
cargo clippy --features dev --all-targets -- -D warnings
```

## License

Licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/ferro-hgvs/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/ferro-hgvs/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
