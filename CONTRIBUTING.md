# Contributing to ferro-hgvs

Thank you for contributing to ferro-hgvs. This page covers what you need for a first PR. The
detail lives in the pages listed under [Where things are documented](#where-things-are-documented).

## Setup

You need Rust (rustup installs the toolchain pinned in `rust-toolchain.toml`),
[cargo-nextest](https://nexte.st/), [pre-commit](https://pre-commit.com/), and, for the Python
bindings, [uv](https://docs.astral.sh/uv/) with Python 3.10+.

```bash
git clone https://github.com/fulcrumgenomics/ferro-hgvs.git
cd ferro-hgvs
git submodule update --init assets/hgvs-nomenclature   # the pinned HGVS spec
scripts/fetch-test-fixtures.sh                         # the bulk test corpora (optional)
pre-commit install                                     # commit and pre-push hooks
cargo nextest run --features dev
```

`dev` is the Cargo feature that turns on the feature-gated tests, examples and generators;
the integration tests do not build without it. Skipping the corpora fetch is safe, but a few
suites then pass without testing anything; `docs/TESTING.md` explains which, and how to make
them fail instead.

For the Python bindings, build with maturin, never with plain `cargo build`, and always name
both features (the reason is on the `extension-module` feature in `Cargo.toml`):

```bash
uv sync --group dev
uv run maturin develop --features python,extension-module
uv run pytest
```

If you change the Python dependencies in `pyproject.toml`, run `uv lock` and commit `uv.lock`.

## Before you push

```bash
cargo fmt --all
cargo clippy --workspace --features dev --all-targets -- -D warnings
cargo clippy --all-features --all-targets -- -D warnings
cargo clippy --release -- -D warnings
cargo nextest run --features dev
cargo test --doc --features dev        # nextest does not run doctests
```

These are the checks CI runs; `.github/workflows/ci.yml` is the source if they drift. A PR
that changes Python code is also checked with `uv run poe check-lint`, `check-format` and
`check-typing`, and one that edits a workflow with `zizmor` and `actionlint`. Document every
public API with a doc comment.

## Opening a PR

- **Give the PR a conventional-commit title**, such as `fix(normalize): correct boundary
  detection for UTR regions`. PRs are squash-merged, so the title becomes the commit on `main`
  and the changelog entry. Types: `feat`, `fix`, `perf`, `refactor`, `test`, `docs`, `ci`,
  `chore`. Branch names and branch commits do not matter.
- **Put `Closes #N` in the description** when the PR resolves an issue.
- **Add a `Representation-Change:` line to the description** if the PR touches a watched
  directory. See the next section.
- **A maintainer merges through the merge queue**, which re-runs the required checks against
  the latest `main`. If the queue drops a green PR without saying why, look at the
  `Merge commit signature` check: a description line that starts with a git commit-header word
  such as `committer` breaks the squash commit's signature. The check's log lists those lines.

## Declaring a representation change

Downstream users store data keyed on ferro's normalized output, so a change to that output
matters even when it is a bug fix. The rule behind this section is rule 7 of the
[normalization rules](docs/src/reference/normalization-rules.md).

If your PR touches a watched directory, the `Representation change declared` check requires a
`Representation-Change:` line in the PR description. The directories are `WATCHED_PREFIXES` in
`scripts/check_representation_change.py`, and the check's message for a missing trailer
lists them, along with the other words that decline. If the change cannot move output,
decline, optionally with a reason:

```
Representation-Change: none. Tests only.
```

If output can move, say so, even when the change is a fix, and give four facts:

1. which forms moved, old and new;
2. in which direction, toward or away from the form that already ships;
3. roughly how many rows moved, and in which corpus;
4. whether the affected inputs were previously rejected or previously accepted. Only a
   previously accepted input has stored data to migrate.

Do not name a release version ("to reproduce pre-vX.Y.Z output"); say "to reproduce output
from before this change". Which release will carry the change is not known when you write it.

```
Representation-Change: c.<cds_end>_*1ins<A> now renders as
  c.<cds_end>delins<ref ++ A>. ~57 of 7,296 rows of the junction-spanning corpus.
  Toward the already-shipped form. Previously-accepted inputs.
```

Write the line at the start of a line with no formatting around it (no backticks, quote marker
or bullet), and indent any continuation lines. Do not decline and then describe a move. To
measure what moved, the header of `examples/dump_normalized_corpus.rs` shows how to compare
two revisions. The full parsing rules are in the docstrings of
`scripts/check_representation_change.py`, and `docs/RELEASE.md` describes how a trailer reaches
the changelog.

## Tests

`docs/TESTING.md` covers running, writing and organizing tests. The short version:

- Integration tests go in `tests/it/`, one module per file, declared as a `mod` in
  `tests/it/main.rs`, or they never run.
- Generate test data in the test. Do not add new committed fixtures.
- Use `cargo nextest`, not `cargo test`, except for doctests.
- Assert the property you care about, not a number that happens to hold today, and check that
  your test can fail.

## Normalization decisions

A decision about what the *correct* normalization output is must be recorded in the same PR,
as a ruling record or a pinned test that cites the spec clause. `docs/TESTING.md` explains
where each kind of decision goes and what to regenerate after editing the ruling ledger.

## Where things are documented

Each topic has one home. Link to it rather than copying it.

| topic | where |
|---|---|
| normalization rules | `docs/src/reference/normalization-rules.md` |
| ruling records | `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`, rendered into `docs/NORMALIZATION_CONTRACT.md` |
| reading the HGVS spec | `docs/READING_THE_SPEC.md` |
| running and writing tests, fixtures, recording decisions, editing the ledger | `docs/TESTING.md` |
| oracle suites | `docs/ORACLES.md` |
| what CI runs | `.github/workflows/ci.yml` |
| representation-change check | `scripts/check_representation_change.py` |
| changelog and releases | `docs/RELEASE.md` |
| tool-support tables | `docs/tool_support_matrix.json` |
| generator rules (`CaptureLedger`) | `tests/it/generator_completeness.rs`, `src/conformance/completeness.rs` |
| agent operating rules | `CLAUDE.md` |

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
