[![CI](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml)
[![Nightly reference-aware tests](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml/badge.svg?branch=main)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml)
[![Codecov](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs)
[![Crates.io](https://img.shields.io/crates/v/ferro-hgvs.svg)](https://crates.io/crates/ferro-hgvs)
[![PyPI](https://img.shields.io/pypi/v/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Python versions](https://img.shields.io/pypi/pyversions/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Documentation](https://docs.rs/ferro-hgvs/badge.svg)](https://docs.rs/ferro-hgvs)
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
- **Variant Normalization**: 3'/5' shifting per HGVS specification
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
```

## Supported HGVS Syntax

| Type | Prefix | Example |
|------|--------|---------|
| Genomic | `g.` | `NC_000001.11:g.12345A>G` |
| Coding DNA | `c.` | `NM_000088.3:c.459A>G` |
| Non-coding | `n.` | `NR_000001.1:n.100A>G` |
| RNA | `r.` | `NM_000088.3:r.459a>g` |
| Protein | `p.` | `NP_000079.2:p.Val600Glu` |
| Mitochondrial | `m.` | `NC_012920.1:m.3243A>G` |

### Edit Types

- Substitution: `A>G`, `Val600Glu`
- Deletion: `del`, `100_200del`
- Insertion: `100_101insATG`
- Deletion-Insertion: `100_102delinsATG`
- Duplication: `100_102dup`
- Inversion: `100_200inv`
- Repeat: `100CAG[20]`

## CLI Commands

The `ferro` CLI provides commands beyond parsing and normalization:

| Command | Description |
|---------|-------------|
| `prepare` | Download and prepare reference data for normalization |
| `check` | Verify reference data setup |
| `parse` | Parse and validate HGVS variants |
| `normalize` | Normalize HGVS variants (3'/5' shifting) |
| `explain` | Explain error/warning codes (e.g., `ferro explain W1001`) |
| `annotate-vcf` | Annotate VCF files with HGVS notation |
| `vcf-to-hgvs` | Convert VCF records to HGVS |
| `hgvs-to-vcf` | Convert HGVS to VCF format |
| `liftover` | Liftover coordinates between genome builds |
| `describe` | Generate HGVS from reference/observed sequences |
| `effect` | Predict protein effect from variant |
| `backtranslate` | Reverse translate protein to DNA variants |
| `convert-gff` | Convert GFF3/GTF to transcripts.json |
| `generate` | Generate HGVS descriptions from components |
| `extract-hgvs` | Extract HGVS from VEP-annotated VCFs |

## Error Handling

ferro-hgvs provides configurable error handling with three modes:

| Mode | Behavior |
|------|----------|
| `strict` | Reject non-conformant input (default) |
| `lenient` | Auto-correct with warnings |
| `silent` | Auto-correct silently |

```bash
# Use lenient mode to auto-correct common issues
ferro parse --error-mode lenient "p.val600glu"  # Corrects to p.Val600Glu

# Ignore specific warnings
ferro parse --ignore W1001,W2001 "p.val600glu"

# Get help on any error/warning code
ferro explain W1001
ferro explain --list
```

### Configuration File

Create `.ferro.toml` in your project directory:

```toml
[error-handling]
mode = "lenient"
ignore = ["W1001", "W2001"]  # Silently correct these
reject = ["W3003"]           # Always reject these
```

## Why ferro-hgvs?

ferro-hgvs provides the most comprehensive HGVS variant normalization across all pattern types, with performance orders of magnitude faster than alternatives.

### Normalization Capabilities Comparison

<!-- DO NOT EDIT — generated from docs/tool_support_matrix.json by `generate_tool_support_tables`. -->
<!-- BEGIN tool-support:normalization_capabilities -->
| Pattern Type | ferro | mutalyzer | biocommons | hgvs-rs |
|--------------|:-----:|:---------:|:----------:|:-------:|
| Genomic (g.) | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) exonic | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) intronic | ✓ | ✓** | ✗ | ✗ |
| Non-coding (n.) | ✓ | ✓ | ✓ | ✓ |
| RNA (r.) | ✓ | ✓ | ✓ | ✓ |
| Protein (p.) | ✓ | Net* | ✗ | ✗ |

\* mutalyzer protein normalization requires network access for NP_→NM_ lookups (cannot be cached locally).
\*\* mutalyzer intronic support is enabled by default via genomic-context rewriting; disable with --no-rewrite-intronic.
<!-- END tool-support:normalization_capabilities -->

### Performance Comparison

**All tools are benchmarked in ferro's offline configuration — best case for every tool.** Reference data is preloaded locally (a local UTA database and SeqRepo) and the network is disabled, so the figures below measure parse/normalize compute, not I/O. Out of the box, hgvs-rs, biocommons/hgvs, and mutalyzer resolve each variant against a remote UTA/SeqRepo or the Mutalyzer web API — a network round-trip per variant (~100–1000 ms), i.e. roughly **1–10 variants/sec, hundreds to thousands of times slower than shown here** (an order-of-magnitude estimate from per-call network latency, not separately benchmarked). That local, offline setup is exactly what ferro's `prepare` command builds; ferro needs no external service.

<!-- DO NOT EDIT — generated from data/benchmark/perf_results.json by `generate_perf_tables`. -->

_Median patterns/sec over 5 reps on an Apple M2 Max, local/offline. All tools draw from one stratified ClinVar population; per-tool sample sizes are calibrated so each tool is measured over a meaningful interval — fast cells (e.g. ferro/hgvs-rs parse) draw from millions of patterns, while slower cells (e.g. the per-tool normalize columns) draw from as few as tens to thousands. All tools exclude process/interpreter startup from the timed region — the mutalyzer/biocommons Python subprocesses are timed by their own internal startup-excluded timer, matching ferro/hgvs-rs. Only ferro parallelizes natively (rayon); the other tools are single-threaded libraries, so their normalize *@8 workers* figures come from the benchmark harness running 8 independent instances in parallel, while parsing is not sharded for them — hence the `single-threaded` label in their parse *@8 workers* column (mutalyzer normalize likewise shows no gain at 8 workers: per-call cache and IPC overhead dominate, so sharding does not help). Every tool runs fully offline against local reference data — a local UTA database and SeqRepo, with mutalyzer's network lookups disabled — the configuration ferro's `prepare` command enables; the figures therefore reflect compute throughput, not per-variant network latency. Reference-data load is excluded for all tools. ferro full-population peak: parse 20.0M/s, normalize 77.0k/s. See `docs/BENCHMARK_RUNBOOK.md` for the full method._

**Parse**

<!-- BEGIN perf:parse -->
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 5.1M/s | 12.2M/s | — |
| mutalyzer | 352/s | single-threaded | 35,000× |
| biocommons | 3.9k/s | single-threaded | 3,100× |
| hgvs-rs | 3.6M/s | single-threaded | 3× |
<!-- END perf:parse -->

**Normalize**

<!-- BEGIN perf:normalize -->
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 78.1k/s | 260.2k/s | — |
| mutalyzer | 4/s | 4/s | 73,000× |
| biocommons | 368/s | 818/s | 320× |
| hgvs-rs | 195/s | 1.3k/s | 200× |
<!-- END perf:normalize -->

**ferro thread scaling**

<!-- BEGIN perf:ferro_scaling -->
| Threads | 1 | 2 | 4 | 8 |
|---------|--:|--:|--:|--:|
| ferro parse | 5.1M/s | 9.4M/s | 16.0M/s | 12.0M/s |
<!-- END perf:ferro_scaling -->

**Input ordering matters for batch throughput.** Resolving a transcript (reading its full sequence from the reference and rebuilding its CDS/exon metadata) dominates per-variant cost. ferro memoizes resolved transcripts in a bounded in-memory cache, so repeated lookups of the same transcript are near-free. Providing variants **sorted by transcript accession — or by genomic position, which clusters variants onto the same transcripts** — maximizes the cache hit rate and can speed up large batches by an order of magnitude versus randomly-ordered input. Ordering matters most when the number of distinct transcripts in the run exceeds the cache capacity (very large or genome-wide inputs); below that, the working set stays resident regardless of order.

### Reference Data: What ferro Prepares

The `ferro prepare` command downloads and organizes all reference data needed for comprehensive normalization. This data is then shared with other tools (mutalyzer, biocommons, hgvs-rs) to enable their local operation.

| Data Type | Source | Size | Enables |
|-----------|--------|------|---------|
| **RefSeq transcripts** | NCBI | ~1GB | NM_/NR_/XM_ normalization |
| **cdot metadata** | MANE | ~200MB | Transcript-to-genome mappings |
| **GRCh38 + GRCh37 genomes** | NCBI | ~4GB | NC_ genomic normalization |
| **RefSeqGene** (sequences + genome alignments) | NCBI | ~600MB | NG_ gene-region normalization; projecting c./n. variants into an NG_ parent's own g. frame (via the RefSeqGene→genome alignment GFF3) |
| **LRG sequences + XML** | EBI | ~50MB | LRG_ stable-reference normalization; projecting c./n. variants into an LRG_ parent's own g. frame (via the LRG XML genomic mapping) |
| **Protein sequences** | Derived from CDS | ~200MB | NP_/XP_ protein normalization |
| **Legacy transcript versions** | NCBI | ~50MB | Historical ClinVar variants |

**Key insight**: Without ferro's reference preparation, other tools require network access for each variant lookup (adding 100-1000ms latency per variant). With ferro's cached reference data, all tools can operate fully offline with consistent, reproducible results.

#### Deriving version-independent NG_ placements (#728)

`ferro prepare --derive-ng-placements <accessions.txt>` derives genomic placements for the listed `NG_` versions (one exact accession per line, e.g. `NG_012337.3`; blank lines and `#` comments ignored), writing `derived_refseqgene_placements.json` into the reference directory and wiring the manifest's `derived_refseqgene_placements` field. This fills version gaps the archived RefSeqGene→genome GFF3 snapshots do not cover. It needs cdot + the genome in the same prepare run and uses NCBI EFetch per accession; accessions that cannot be validated are skipped with a warning. The field is preserved across subsequent `prepare` runs.

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
cargo test
cargo clippy -- -D warnings
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
