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
cargo build --features benchmark     # Build ferro-benchmark binary
```

### Python Bindings
```bash
cargo check --features python        # Typecheck the bindings (fast; does NOT link)
cargo clippy --features python       # Lint the bindings (also does NOT link)
maturin develop --features python    # Build and install locally for development
maturin build --features python      # Build wheel
```

**Do not use `cargo build --features python`** — on macOS it cannot link. `pyo3` is declared with
`features = ["abi3-py310", "extension-module"]`, and `extension-module` deliberately does *not*
link against libpython: the CPython symbols are meant to resolve against the host interpreter
when the module is imported. Maturin supplies the linker allowance that makes that legal
(`-undefined dynamic_lookup`); a plain `cargo build` does not, so it compiles the whole crate and
only then dies linking the cdylib:

```
error: linking with `cc` failed: exit status: 1
  = note: Undefined symbols for architecture arm64:
            "_PyBaseObject_Type", referenced from: ...
          ld: symbol(s) not found for architecture arm64
error: could not compile `ferro-hgvs` (lib) due to 1 previous error
```

Compilation had already succeeded at that point, so the failure says nothing about the code —
it just costs a full build cycle to find out. Use `cargo check --features python` (and
`cargo clippy --features python`) to verify the bindings compile, and `maturin develop
--features python` whenever you need a module you can actually import (e.g. to run
`pytest tests/python/`).

The hard failure above is macOS-specific: `ld` there refuses undefined symbols in a dylib unless
`-undefined dynamic_lookup` is passed. Two things could pass it and neither does under a plain
`cargo build`. PyO3 ships the flags in
`pyo3_build_config::add_extension_module_link_args()` — Darwin-only, by design — but a crate has
to call that from its own `build.rs`, and this crate has no `build.rs`. Maturin passes them
itself, which is why the `maturin` commands above work. So `cargo build` gets them from nowhere
and the link dies.

Linkers that do permit undefined symbols in a shared object get no such error, but the advice is
unchanged everywhere: `cargo build` still produces no importable extension module, so reach for
`maturin` rather than treating a quiet link as success.

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

The check runs from `Normalizer::assert_seam_oracles`, which all three oracles
share, at the single exit of `normalize_core_checked`. That covers every public
normalization: `normalize()`, `normalize_with_diagnostics()`, and every
`VariantProjector` path (the projection-driven genomic/coding/protein axes).

It is one call site again as of #1382. `normalize_with_diagnostics` used to reach
`normalize_core_canonical` directly, which bypassed both the oracles and — the
actual defect — the strict-mode rejection ladder, so #1366 had to call
`assert_seam_oracles` from it separately. Routing it through
`normalize_core_checked` fixes both at once, and the extra call would now just
re-run all three oracles on that path.

The idempotency oracle's verification pass re-enters normalization, so a
thread-local `IN_IDEMPOTENCY_CHECK` guard breaks the recursion — the inner call
skips its own check.

#### Normalization re-parse oracle

`FERRO_ASSERT_REPARSE=1` is the second half of the same seam: it asserts that a
normalized description is one `parse_hgvs` accepts — *when normalization is what
broke it*. Its exemptions are a **closed list**, deliberately (#1264):

- `0` and `?` are legal whole-allele outputs that `parse_hgvs` rejects standalone
  because it wants an accession — a limit of the entry point, not a malformed output.
- An **empty allele** (`[]`), reachable only by direct construction; the projector's
  own tests build one on purpose to pin that it declines.
- A **non-flanking genomic insertion**, which is the projection *pivot*: its
  coordinates are sound and the downstream cdot derivation of the c./n./p. axes needs
  it, but its spelling is not one HGVS admits, so the projector withholds the
  *reported* genomic axis instead (`non_flanking_genomic_insertion_anchor`).

**This used to be a blanket** — "skip when the input does not itself re-parse" — and
that blanket was hiding live defects rather than scoping the oracle. Instrumenting it
to report instead of return silently, then re-running the suite, found 18 hits in four
shapes, two of which were real projector bugs (the RNA-only `spl` edit carried onto the
`g.`/`c.` axes, and an insertion whose anchors straddle a splice junction). Both are
fixed; the exemption is now narrow enough that the same class of defect fails loudly.
If you find yourself widening it, that is the signal to fix the producer instead.

Like the idempotency oracle it is compiled out in release builds
(`#[cfg(debug_assertions)]`) and read once into a `OnceLock`, so it adds no
release-build cost and a disabled debug run pays only one atomic load.

```bash
FERRO_ASSERT_REPARSE=1 cargo nextest run --features dev
```

Kept separate from `FERRO_ASSERT_IDEMPOTENT` because neither subsumes the other,
and the idempotency oracle has a blind spot this one covers: it verifies by
*re-normalizing* its own output, which it cannot do for an output that fails to
parse, so an unparseable result is invisible to it. Both run from
`assert_seam_oracles` alongside the in-bounds check, and CI sets all three together.

#### Normalization in-bounds oracle

`FERRO_ASSERT_IN_BOUNDS=1` is the third at the same seam: no coordinate a
normalized description names may be past the end of its own sequence.

```bash
FERRO_ASSERT_IN_BOUNDS=1 cargo nextest run --features dev
```

It exists because the class kept being found by hand, one shape at a time —
**#1274** (`T:g.[8_9insA;10del]` -> `g.10_11=` on a 10-base contig), **#1343**
(`c.[*10dup;*11dup]` -> `c.*11_*12insAA`) and **#1307**
(`g.[24dup;24C>G]` -> `g.[24C>G;24_25insC]`) — each filed, fixed and
regression-tested separately. Three instances of one defect class is the argument
for an invariant at the seam rather than a fourth per-shape test.

Neither of the other two asks the question. `FERRO_ASSERT_REPARSE` cannot:
`parse_hgvs` holds no provider, so `g.24_25insC` is well-formed to it.
`FERRO_ASSERT_IDEMPOTENT` catches some instances incidentally — the #1307 output
re-normalizes to `g.23_24insG`, so it is not a fixed point — but an out-of-range
output that *is* a fixed point passes it, which is the #1327 shape
(`m.16569_16570insAA`).

What makes it non-trivial, and what it deliberately does not check, is documented
on `merge::first_out_of_bounds_coordinate`. In brief:

- **Positions are converted onto the served sequence's axis first.** A `c.`
  position may legitimately exceed the CDS length (`*n` into the 3'UTR, `-n` into
  the 5'UTR), so `region_sequence_delta` runs before any comparison.
- **A reversed range is not an error.** SVD-WG006 admits `<high>_<low>` for a
  circular deletion or duplication, so the two endpoints are checked
  independently and their order is never compared.
- **An authored overrun is exempt** — W4004 `PositionPastEnd` is what reports
  those. Detecting it requires reading each endpoint independently, because a
  special position (`pter`) on one end otherwise hides a past-end coordinate on
  the other.
- **Not covered:** protein axes, and an inserted-range payload
  (`g.10_11ins[20_30]`).

Like the two above it, compiled out in release (`#[cfg(debug_assertions)]`) and
read once into a `OnceLock`. CI sets all three together, in the sharded
`test-oracle` job and the `sweeps` job (those two and nowhere else — the plain
`test` and `soak` jobs run without the flags); the nightly sets it too, where it
is the only place the check runs against true transcript and contig lengths.

**An oracle failure is visible in PR CI and silent in the nightly.** In PR CI,
`test-oracle` and `sweeps` carry no `continue-on-error`, so a fire turns the job
red. The nightly reference-aware job does carry it, deliberately — its purpose is
to surface drift in the xfail report rather than to gate, and the corpus runner
wraps normalization in `catch_panics`, so an oracle panic there lands in the
uploaded xfail artifact as a FAILing case instead of failing the workflow. That
applies equally to all three flags, and it is why a red oracle in the nightly must
be read out of the xfail report rather than from the workflow conclusion.

Red is not the same as blocked, though: `main`'s ruleset requires only eight
checks — Build, Test, Clippy, Clippy (all features), Format, Python Lint, Python
Wheel Test, abi3 floor — and neither `Test oracle` nor `Exhaustive sweeps` is
among them. So an oracle fire shows up as a failed job that a human must not
merge past, not as a ruleset-enforced merge block.

#### Exhaustive cis sweeps: `FERRO_SWEEP_SEEDS`

The three exhaustive cis sweeps — `cis_junction_crossing_shift`,
`repeat_span_sibling_overlap` and `issue_1254_sibling_crossing_shift` — draw a
deterministic corpus of 20-mers from `sweep_sequences`. They dominated a local
run: measured at 79.6s of an 86.6s suite for `cis_junction_crossing_shift`
alone, so any `cargo nextest run` cost ~90s regardless of what changed (#1295).

They now take a **4-seed prefix of that corpus by default**, and CI's
`Exhaustive sweeps` job asks for all of it:

```bash
cargo nextest run --features dev                        # 4-seed prefix (fast)
FERRO_SWEEP_SEEDS=full cargo nextest run --features dev # the full corpus, as CI runs it
FERRO_SWEEP_SEEDS=12   cargo nextest run --features dev # an explicit count, for bisecting
```

Local suite: **86.6s → 27.5s**, with no test left over nextest's 60s SLOW
threshold.

**Only the sequence diversity is reduced — never the shapes.** The loop bounds
over member spellings, positions and directions are untouched, which is the
axis worth keeping: each sweep's notes record that every blocking defect found
so far lived in a shape the generator could not emit, not in a sequence it did
not draw.

**Why the pinned counts did not become profile-dependent**, which is the hazard
#1295 flags against exactly this change. `sweep_sequences` is prefix-stable
(`sweep_sequences_is_prefix_stable`), so a reduced run enumerates a strict
*subset* of the full run's cases. An `is_empty()` assertion therefore cannot
pass at the prefix while failing at the full corpus, and
`FIVE_PRIME_DUP_DEL_SEQUENCE_CHANGES` is pinned at **zero**, which is
subset-stable for the same reason — had that residual still been 74, this knob
would have had to make it profile-dependent. The one genuinely seed-dependent
assertion in each sweep is its case-count floor, which is why those are written
per-seed rather than as absolutes.

The one thing to watch: `sweeps` is the only job that sets the variable, and it
selects on `SWEEP_FILTER`. Moving a sweep out of that filter silently drops it
to the prefix in CI, so edit the two together.

### Generated spec fixture (not committed)

`tests/fixtures/grammar/hgvs_spec_normalization.json` is a **generated build artifact** — it is `.gitignore`d, not committed. It is produced by the `generate_spec_fixture` binary from the HGVS spec submodule, the parser's behavior, and the curated `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`. (Committing it made every parser PR a merge-conflict magnet, since each PR regenerated the whole file.)

```bash
cargo run --features dev --bin generate_spec_fixture            # (re)generate it
cargo run --features dev --bin generate_spec_fixture -- --output <path>   # write elsewhere
```

Replace `<path>` with the destination you want; the first form writes to the default
location above. Do not substitute a machine-specific absolute path into a committed file
(see [Repository Hygiene](#no-machine-local-paths-in-committed-files)).

Two test modules actually read it — `tests/it/hgvs_spec_normalization_tests.rs` (the driver) and `tests/it/idempotency_tests.rs` (which uses it purely as an input corpus, so its assertions are not circular). They call a shared helper (`tests/it/common/spec_fixture.rs`) that regenerates it on demand if missing, so a fresh `cargo test` "just works" (one-time generation). `mito_circular_audit.rs`, `protein_silent_eq.rs` and `protein_unknown_roundtrip.rs` only *mention* the fixture in doc comments as a cross-link; they do not read it. (This list previously named all five as consumers — corrected 2026-08-02.) CI and the pre-push hook regenerate it explicitly before the test run — plain generation, not `--check`, because generation is what validates the committed overrides against the spec checkout.

**The fixture cannot detect a regression on its own, and neither can the enumeration's (#1272).** Both artifacts are generated from the code under test and regenerated by CI before the run, so `pinned_v21_normalization_behavior` and `enumeration_replays_recorded_behavior` compare ferro against itself — they are stale-local-artifact detectors. The guards that actually judge behaviour are the **committed** ones, and they are the ones to read when either drifts:

| guard | file | what it pins |
|---|---|---|
| `ferro_produces_the_form_the_spec_states` | `hgvs_spec_normalization_tests.rs` | ferro's output vs the **spec's** stated form, plus `KNOWN_DIVERGENT_INPUTS` |
| `status_census_is_unchanged` | `hgvs_spec_normalization_tests.rs` | per-status row counts |
| `DIVERGENCE_BUDGET` / `PASSING_CENSUS` | `spec_enumeration_tests.rs` | per-status row counts, together totalling every row |

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
