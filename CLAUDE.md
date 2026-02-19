# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ferro-hgvs is a high-performance HGVS variant nomenclature parser and normalizer written in Rust with Python bindings. It supports all HGVS coordinate systems (g/c/n/r/p/m/o) and edit types (substitution, deletion, insertion, duplication, inversion, repeat).

## Build Commands

### Rust
```bash
cargo build                          # Debug build
cargo build --release                # Release build
cargo build --features dev           # Build with all testable features
cargo build --features python        # Build with Python bindings (use maturin instead)
cargo build --features benchmark     # Build ferro-benchmark binary
```

### Python Bindings
```bash
maturin develop --features python    # Build and install locally for development
maturin build --features python      # Build wheel
```

## Testing

### Rust Tests
```bash
cargo nextest run --features dev     # Run all tests (preferred)
cargo test --features dev            # Alternative
cargo nextest run -E 'test(parse)'   # Run specific tests by name pattern
```

### Python Tests
```bash
pytest tests/python/ -v              # Run all Python tests
pytest tests/python/test_core.py -v  # Run specific test file
pytest -k "test_parse" -v            # Run tests matching pattern
```

## Linting and Formatting

### Rust
```bash
cargo fmt --all                      # Format code
cargo clippy --features dev -- -D warnings  # Lint
```

### Python
```bash
ruff check python/ tests/python/     # Lint
ruff format python/ tests/python/    # Format
mypy python/ferro_hgvs/              # Type check
```

### Pre-commit (runs all checks)
```bash
pre-commit install                   # Install hooks (one-time)
pre-commit run --all-files           # Run manually
```

## Architecture

### Core Modules (`src/`)

- **`hgvs/`** - HGVS parsing and variant types
  - `parser/` - nom-based parser (fast_path.rs for hot paths)
  - `variant.rs` - HgvsVariant enum (Cds, Genome, NonCoding, Rna, Protein, Mito)
  - `edit.rs` - Edit types (Substitution, Deletion, Insertion, etc.)
  - `location.rs` - Position types with offset support

- **`normalize/`** - Variant normalization (3'/5' shifting)
  - `shuffle.rs` - Core shuffling algorithm
  - `rules.rs` - HGVS-specific normalization rules

- **`reference/`** - Reference sequence providers
  - `provider.rs` - ReferenceProvider trait
  - `fasta.rs` - FASTA file provider
  - `mock.rs` - MockProvider for testing

- **`convert/`** - Coordinate system conversions (c. ↔ g. ↔ n. ↔ p.)

- **`spdi/`** - SPDI format support and HGVS↔SPDI conversion

- **`vcf/`** - VCF parsing and HGVS annotation

- **`error_handling/`** - Configurable error modes (strict/lenient/silent)

- **`python.rs`** - PyO3 bindings exposing Rust API to Python

### CLI Binaries (`src/bin/`)

- `ferro.rs` - Main CLI (`ferro parse`, `ferro normalize`, `ferro prepare`, etc.)
- `benchmark.rs` - Tool comparison benchmark (requires `--features benchmark`)
- `ferro-web.rs` - Web service (requires `--features web-service`)

### Python Package (`python/ferro_hgvs/`)

- `__init__.py` - Re-exports from native extension
- `__init__.pyi` - Type stubs for IDE support

## Feature Flags

- `dev` - All testable features combined
- `python` - PyO3 Python bindings (build with maturin)
- `benchmark` - ferro-benchmark binary and tool comparison
- `validation` - Hash-based validation with ahash
- `parallel` - Rayon-based parallelism
- `web-service` - HTTP API server

## Key Patterns

### Parsing
```rust
use ferro_hgvs::{parse_hgvs, HgvsVariant};
let variant = parse_hgvs("NM_000088.3:c.459A>G")?;
```

### Python API
```python
import ferro_hgvs
variant = ferro_hgvs.parse("NM_000088.3:c.459A>G")
```

### Reference Data
The `ferro prepare` command downloads RefSeq transcripts, genome FASTAs, and cdot metadata needed for normalization. Data goes to a reference directory checked with `ferro check --reference <dir>`.
