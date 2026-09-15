# CLAUDE.md

Guidance for Claude Code when working in this repository. Keep this file short: it is
operational. Design rationale, measurements and decision history live in the places
this file points to, and nowhere else.

## Project overview

ferro-hgvs is an HGVS variant nomenclature parser and normalizer in Rust with Python
bindings. It supports all HGVS coordinate systems (g/c/n/r/p/m/o) and edit types.

## Build

```bash
cargo build                          # debug
cargo build --release
cargo build --features dev           # everything testable
cargo build --features benchmark     # ferro-benchmark binary
```

Python bindings (PyO3, built with maturin, never with plain `cargo build`):

```bash
cargo check --features python                          # typecheck only
maturin develop --features python,extension-module     # importable module
maturin build   --features python,extension-module     # wheel
```

`extension-module` is a separate feature from `python`. Always name it explicitly on
maturin commands; `[tool.maturin] features` in `pyproject.toml` is a default that a
`--features` flag replaces. See `CONTRIBUTING.md` for why.

## Test

```bash
cargo t                              # alias: nextest, --features dev, --lib --test it
cargo ta                             # alias: the full suite. Run before pushing.
cargo nextest run -E 'test(parse)'   # by name
pytest tests/python/ -v              # Python bindings (needs maturin develop first)
```

Rules that will bite you:

- Use `cargo nextest`, not `cargo test`, for the lib suite. Two lib tests are not safe
  under libtest's shared thread pool (`normalize::merge::tests`, `parallel::tests`).
- Do not export `CARGO_INCREMENTAL=0`. It disables incremental builds for the very large
  `it` test crate and makes every one-line test edit a full rebuild. Fast-iteration
  recipe: `docs/TESTING.md`.
- `FERRO_ASSERT_IDEMPOTENT=1 cargo nextest run` is always red on `main`. Run the oracles
  through `scripts/run_oracle_suite.sh`, which mirrors CI's selection and flags.
  Details: `docs/ORACLES.md`.
- The exhaustive sweeps run a short seed prefix by default. `FERRO_SWEEP_SEEDS=full` is
  what CI runs.
- `tests/fixtures/grammar/hgvs_spec_normalization.json` is generated and gitignored.
  Tests regenerate it on demand; it needs the `assets/hgvs-nomenclature` submodule.
- Manifest-backed and bulk-corpus suites return early and report PASS when their data is
  absent. Use `scripts/run_conformance_axis.sh` and `FERRO_REQUIRE_BULK_FIXTURES=1`.
  See `CONTRIBUTING.md`.

## Lint and format

```bash
cargo fmt --all
cargo clippy --features dev -- -D warnings
ruff check python/ tests/python/ && ruff format python/ tests/python/
mypy python/ferro_hgvs/
pre-commit run --all-files
```

## Required PR checks that are easy to miss

- A PR touching `src/normalize/`, `src/hgvs/`, `src/spdi/`, `src/project/`,
  `src/reference/` or `src/error_handling/` must carry a `Representation-Change:`
  trailer **in the PR description**. `Representation-Change: none` is a valid answer.
  How to measure and how to phrase it: `CONTRIBUTING.md`, "Declaring a representation
  change". Keep the six paths above backticked; a Python test reads them from this file.
- The required `Test` context is a rollup over several jobs. If it is red while every
  `Test (n/N)` shard is green, read the rollup's log to find the upstream job. Read the
  `needs:` list from `.github/workflows/ci.yml`, not from prose.

## Architecture

- `src/hgvs/` parsing and variant types (`parser/`, `variant.rs`, `edit.rs`, `location.rs`)
- `src/normalize/` normalization: shifting, merging, canonical partition (`merge.rs`)
- `src/reference/` reference providers (`fasta.rs`, `mock.rs` = `JsonProvider`)
- `src/convert/` c./g./n./p. coordinate mapping
- `src/project/` cross-axis projection (`VariantProjector`)
- `src/spdi/` SPDI conversion, also the applier the oracles use
- `src/vcf/` VCF parsing and annotation
- `src/conformance/` conformance corpora and the `CaptureLedger`
- `src/error_handling/` strict / lenient / silent modes
- `src/python.rs` PyO3 bindings; `python/ferro_hgvs/` package and `.pyi` stubs
- `src/bin/ferro.rs` CLI; `benchmark.rs` and `ferro-web.rs` are feature-gated

Feature flags are declared in `Cargo.toml`. `dev` turns on everything testable. `parallel`
is a no-op kept for compatibility.

## Conventions

- **No machine-local absolute paths in committed files.** Take locations from a CLI
  argument, an env var (`FERRO_MANIFEST`), or a repo-relative default.
- **Never hand-edit `tests/it/clause_ruling_index.rs`.** It embeds a generated block; the
  test that fails tells you how to regenerate it.
- **Record every adjudication in the same change.** A decision about what the *correct*
  normalization is goes into a ruling record in
  `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, an equivalence class,
  `KNOWN_DIVERGENT_INPUTS`, or an ordinary pinned test, with the spec clause cited as
  `file:line`. `undecided` and `house_choice` are first-class. A `decided` record carries a
  one-sentence `summary`; the generator refuses one without it. The rendered ledger is
  `docs/NORMALIZATION_CONTRACT.md` (generated, do not edit); the ruleset it sits under
  is `docs/src/reference/normalization-rules.md`. Read both before arguing from spec text.
- **Assert the property, not the number.** Import the constant you are guarding; do not
  restate it. A pinned count is a change detector, not a guard. Before quoting a `0` from
  a corpus, show the generator can build the shape you changed.
- **A generator must account for what it dropped.** Route the last fallible step before a
  write through `CaptureLedger` and call `finish()`. See `CONTRIBUTING.md`.
- **Do not put measurements, dates, issue narratives or line numbers in this file.**
  They go stale within days here. Put them in the PR, the issue, or a `docs/*.md` page.

## Reading the HGVS spec

Almost no clause in `assets/hgvs-nomenclature/docs/recommendations/` is RFC 2119
normative, the spec states no minimality principle, and several forward-looking notes
describe proposals that were later rejected. Do not argue from keyword strength or from
"the direction of travel". The house guide to reading it, with the checks behind each
claim, is `docs/READING_THE_SPEC.md`.
