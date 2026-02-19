[![CI](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs)
[![Crates.io](https://img.shields.io/crates/v/ferro-hgvs.svg)](https://crates.io/crates/ferro-hgvs)
[![Documentation](https://docs.rs/ferro-hgvs/badge.svg)](https://docs.rs/ferro-hgvs)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

# ferro-hgvs

A high-performance HGVS variant nomenclature parser and normalizer written in Rust.

**WARNING: ALPHA SOFTWARE - USE AT YOUR OWN RISK**

This software is currently in **ALPHA**. While we have extensively tested it
across a wide variety of HGVS patterns, **no guarantees are made** regarding
correctness or stability.

<p>
<a href="https://fulcrumgenomics.com"><img src="https://raw.githubusercontent.com/fulcrumgenomics/fgumi/main/.github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Features

- **Full HGVS Parsing**: All coordinate systems (g/c/n/r/p/m/o) and edit types
- **Variant Normalization**: 3'/5' shifting per HGVS specification
- **High Performance**: ~2.5M variants/sec parsing, zero-copy with nom
- **Type-Safe**: Leverages Rust's type system for correctness

## Installation

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

# Prepare reference data (downloads RefSeq, genome, cdot)
ferro prepare --output-dir ferro-reference

# Verify reference data is ready
ferro check --reference ferro-reference

# Normalize with reference
ferro normalize "NM_000088.3:c.459del" --reference ferro-reference/
```

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
reject = ["W4002"]           # Always reject these
```

## Why ferro-hgvs?

ferro-hgvs provides the most comprehensive HGVS variant normalization across all pattern types, with performance orders of magnitude faster than alternatives.

### Normalization Capabilities Comparison

| Pattern Type | ferro | mutalyzer | biocommons | hgvs-rs |
|--------------|:-----:|:---------:|:----------:|:-------:|
| Genomic (g.) | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) exonic | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) intronic | ✓ | ✓* | ✗ | ✗ |
| Non-coding (n.) | ✓ | ✓ | ✓ | ✓ |
| RNA (r.) | ✓ | ✓ | ✓ | ✓ |
| Protein (p.) | ✓ | Net** | ✗ | ✓ |

\* mutalyzer intronic support requires genomic context rewriting (enabled by default)
\** mutalyzer protein normalization requires network access for NP_→NM_ lookups

### Performance Comparison

| Tool | Speed (local) | Speed (network) | ferro Speedup |
|------|---------------|-----------------|---------------|
| **ferro-hgvs** | ~4M patterns/sec | N/A (offline) | — |
| mutalyzer | ~20 patterns/sec | ~1 pattern/sec | **200,000x** |
| biocommons/hgvs | ~20 patterns/sec | ~0.2 patterns/sec | **200,000x** |
| hgvs-rs | ~2 patterns/sec | ~0.2 patterns/sec | **2,000,000x** |

### Reference Data: What ferro Prepares

The `ferro prepare` command downloads and organizes all reference data needed for comprehensive normalization. This data is then shared with other tools (mutalyzer, biocommons, hgvs-rs) to enable their local operation.

| Data Type | Source | Size | Enables |
|-----------|--------|------|---------|
| **RefSeq transcripts** | NCBI | ~1GB | NM_/NR_/XM_ normalization |
| **cdot metadata** | MANE | ~200MB | Transcript-to-genome mappings |
| **GRCh38 + GRCh37 genomes** | NCBI | ~4GB | NC_ genomic normalization |
| **RefSeqGene** | NCBI | ~600MB | NG_ gene region normalization |
| **LRG sequences** | EBI | ~50MB | LRG_ stable reference normalization |
| **Protein sequences** | Derived from CDS | ~200MB | NP_/XP_ protein normalization |
| **Legacy transcript versions** | NCBI | ~50MB | Historical ClinVar variants |

**Key insight**: Without ferro's reference preparation, other tools require network access for each variant lookup (adding 100-1000ms latency per variant). With ferro's cached reference data, all tools can operate fully offline with consistent, reproducible results.

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
