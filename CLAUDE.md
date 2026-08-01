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

#### Normalization idempotency oracle

`FERRO_ASSERT_IDEMPOTENT=1` turns every `Normalizer::normalize` call into an
assertion that `norm(norm(x)) == norm(x)`, so any test that normalizes becomes an
idempotency check:

```bash
FERRO_ASSERT_IDEMPOTENT=1 cargo nextest run --features dev
```

It is compiled out in release builds (`#[cfg(debug_assertions)]`) and read once
into a `OnceLock`, so a disabled run pays only one atomic load. CI runs the suite
a second time with it set; the nightly reference-aware job sets it too, which is
where the manifest-backed conformance corpora are actually covered.

The check sits at the single shared exit of `normalize_core_checked`, so it covers
both `normalize()` and every `VariantProjector` path (the projection-driven
genomic/coding/protein axes). Its verification pass re-enters that same method, so
a thread-local `IN_IDEMPOTENCY_CHECK` guard breaks the recursion — the inner call
skips its own check.

#### Normalization re-parse oracle

`FERRO_ASSERT_REPARSE=1` is the second half of the same seam: it asserts that a
normalized description is one `parse_hgvs` accepts — *when normalization is what
broke it*. Two exemptions keep that scope honest: `0` and `?` are legal
whole-allele outputs that `parse_hgvs` rejects standalone because it wants an
accession, and an input that does not itself re-parse is passed over, since a
description the parser accepts inbound but cannot re-read outbound is a
parse/display round-trip asymmetry normalize merely carries through. Both are a
different bug; reporting them here would bury the one this oracle is for.

```bash
FERRO_ASSERT_REPARSE=1 cargo nextest run --features dev
```

Kept separate from `FERRO_ASSERT_IDEMPOTENT` because neither subsumes the other,
and the idempotency oracle has a blind spot this one covers: it verifies by
*re-normalizing* its own output, which it cannot do for an output that fails to
parse, so an unparseable result is invisible to it. Both sit at
`normalize_core_checked`'s single exit and CI sets both together.

### Generated spec fixture (not committed)

`tests/fixtures/grammar/hgvs_spec_normalization.json` is a **generated build artifact** — it is `.gitignore`d, not committed. It is produced by the `generate_spec_fixture` binary from the HGVS spec submodule, the parser's behavior, and the curated `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`. (Committing it made every parser PR a merge-conflict magnet, since each PR regenerated the whole file.)

```bash
cargo run --features dev --bin generate_spec_fixture            # (re)generate it
cargo run --features dev --bin generate_spec_fixture -- --output <path>   # write elsewhere
```

Replace `<path>` with the destination you want; the first form writes to the default
location above. Do not substitute a machine-specific absolute path into a committed file
(see [Repository Hygiene](#no-machine-local-paths-in-committed-files)).

The tests that read it (`tests/it/hgvs_spec_normalization_tests.rs`, `tests/it/idempotency_tests.rs`, `tests/it/mito_circular_audit.rs`, `tests/it/protein_silent_eq.rs`, `tests/it/protein_unknown_roundtrip.rs`) call a shared helper (`tests/it/common/spec_fixture.rs`) that regenerates it on demand if missing, so a fresh `cargo test` "just works" (one-time generation). CI and the pre-push hook regenerate it explicitly before the test run — plain generation, not `--check`, because generation is what validates the committed overrides against the spec checkout.

Requires the `assets/hgvs-nomenclature` submodule (`git submodule update --init assets/hgvs-nomenclature`); without it the generator fails with `no HGVS strings harvested from …`, naming that command.

`--check` answers a different question — "is my local artifact current?" — and is not a gate: an absent artifact is generated rather than reported as drift, since a gitignored file has no committed baseline to drift from. Use it when you want to know whether a code change moved the fixture.

Because the fixture is no longer in git, per-PR `parse-error → preserved` status transitions are reviewed via the PR description and the accompanying test/parser changes, not a committed diff.

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

## Repository Hygiene

### No machine-local paths in committed files

Never commit machine-specific absolute paths (e.g. `/Volumes/...`, `/Users/...`, `/home/...`) to source, tests, examples, docs, or fixtures. They break on every other machine and leak local/client directory layout into a public-org repo.

- **Source / examples / tests:** take the location from a CLI argument, an environment variable (e.g. `FERRO_MANIFEST`), or a repo-relative default (e.g. `benchmark-output/manifest.json`) — never a hardcoded absolute path.
- **Docs / runbooks:** use a placeholder (e.g. `/path/to/ferro-bench-data`) or a shell variable the reader sets, and say so explicitly.
- The prepared-reference manifest itself lives outside the repo and is `.gitignore`d; reference manifest entries are stored as relative paths so the reference is portable across hosts.

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
  - `mock.rs` - `JsonProvider` (JSON/`transcripts.json`-backed; also backs `Normalizer(reference_json=...)`). Former name `MockProvider` retained as a compatibility alias.

- **`convert/`** - Coordinate system conversions (c. ↔ g. ↔ n. ↔ p.)

- **`spdi/`** - SPDI format support and HGVS↔SPDI conversion

- **`vcf/`** - VCF parsing and HGVS annotation

- **`error_handling/`** - Configurable error modes (strict/lenient/silent)

- **`python.rs`** - PyO3 bindings exposing Rust API to Python

### CLI Binaries (`src/bin/`)

- `ferro.rs` - Main CLI (`ferro parse`, `ferro normalize`, `ferro project`, `ferro prepare`, etc.)
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

It also downloads genomic-parent placement data used to project transcript-coordinate variants into a genomic parent's own frame (#480):
- **RefSeqGene→genome alignments** (`GCF_*_refseqgene_alignments.gff3`, from NCBI's archived `alignments/ARCHIVE/all/` — the feed stopped updating in 2024; latest GRCh38 is the RS-109/p13 snapshot, valid for p14 since primary `NC_` accessions are unchanged) → manifest `refseqgene_alignments`; parsed by `MultiFastaProvider` into per-`NG_` `GenomicPlacement`. The corresponding GRCh37 snapshot (`GCF_000001405.25_105.*`) → manifest `refseqgene_alignments_grch37`, merged with the GRCh38 file so an `NG_` resolves to its build-appropriate placement (#653/#713).
- **LRG XML** genomic mapping → per-`LRG_` `GenomicPlacement` (parsed on demand).

`ReferenceProvider::genomic_placement` exposes these; `VariantProjector::project_to_genomic` composes the `NM_`→`NC_` (cdot) step with the `NC_`→`NG_`/`LRG_` affine transform to re-express coordinates in the parent frame. With no placement it declines rather than emit chromosome coordinates under the parent accession.
