# Contributing to ferro-hgvs

Thank you for considering contributing to ferro-hgvs!

## Getting Started

### Prerequisites

- Rust (stable)
- Git
- [uv](https://docs.astral.sh/uv/) (for Python development)

### Setup

```bash
git clone https://github.com/fulcrumgenomics/ferro-hgvs.git
cd ferro-hgvs
cargo build
cargo test --features dev
```

### Python Bindings Setup

```bash
uv sync --group dev
uv run maturin develop --features python
uv run pytest
```

After modifying Python dependencies in `pyproject.toml`, run `uv lock` and commit
the updated `uv.lock`. CI uses `--locked` and will fail if the lockfile is out of
sync with `pyproject.toml`.

## Development Workflow

### Making Changes

1. Create a branch: `git checkout -b feature/your-feature-name`
2. Write your code and add tests
3. Ensure all tests pass: `cargo nextest run --features dev` (or `cargo test --features dev`)
4. Run lints: `cargo clippy --features dev -- -D warnings && cargo fmt --check`

### Commit Messages

Follow conventional commit format:

```
type(scope): description
```

Types: `feat`, `fix`, `docs`, `refactor`, `test`, `chore`

Examples:
```
feat(parser): add support for repeat variants
fix(normalize): correct boundary detection for UTR regions
```

### Declaring a representation change

**Representation stability is a shipped guarantee.** The downstream consumer keys
its read counts on the normalized HGVS string, so what matters is not which
representation ferro picks but whether the pick *moves*. A canonical form that
churns between releases silently re-buckets every stored design — a
re-normalization of the consumer's whole library, not a bugfix they can ignore.

So if your change can alter any normalized output string, say so in the commit,
with a `Representation-Change:` trailer:

```
fix(normalize): clamp a derived insertion at the CDS/3'UTR junction

Representation-Change: c.<cds_end>_*1ins<A> now renders as
  c.<cds_end>delins<ref ++ A>. ~57 of 7,296 rows of the junction-spanning
  corpus slice (48 seeds x both directions x three shape families).
  Toward the already-shipped form: it is what a second normalize already
  produced. Previously-accepted inputs, so a real migration on paper — but
  anyone normalizing to a fixed point sees no change.
```

**If your change moves nothing, say that** — `Representation-Change: none` is a
first-class answer, not a way of opting out:

```
ci: pin the workflow's action hashes

Representation-Change: none
```

`none`, `no`, `n/a` and `na` all count, in any casing. **Say why, if it is not
obvious** — the verdict is the first word, and a full stop, semicolon, colon or
dash may introduce a reason:

```
test(conformance): gate the harvested unguarded cases on every PR

Representation-Change: none. Test-only plus two new `src/conformance/`
  modules; nothing under `src/normalize/`, `src/hgvs/`, `src/spdi/` or
  `src/project/` is touched, and no existing expectation was re-blessed.
```

A comma is *not* a terminator, because it usually introduces a qualification
that changes the verdict ("none, except two rows") — that reads as a description
of a move, and is filed as one. Same for `no rows move`: no terminator, so the
verdict is not `no`. Both err toward listing a change rather than hiding one.

**Do not decline and then describe a move.** `no. 3 rows move` is filed as a
decline, because the verdict is the first word — so the disclosure never reaches
the changelog and nobody is told. CI rejects that combination: a declining
trailer that also claims `<n> rows move` (or merge/split/respell, for any `n`
other than zero) fails the check, and the message asks you to say which it is.
Quantifying a zero is fine and encouraged — `none. 0 of 950 rows move` passes.

A declining trailer is
excluded from the changelog's **Representation changes** section, so declaring it
costs the reader nothing while leaving the judgement on the record. CI enforces
the distinction: a change touching `src/normalize/`, `src/hgvs/`, `src/spdi/` or
`src/project/` with no trailer at all fails the `Representation change declared`
check, because an absent trailer is indistinguishable from an unconsidered one.

**Indent continuation lines**, as above. Git only folds a multi-line trailer
value into the trailer when the continuations are whitespace-prefixed: measured
against git-cliff 2.13.1, the unindented form still groups correctly (its footer
matcher regexes the whole footer paragraph), but `git interpret-trailers --parse`
reports *no trailers at all* for it. So an unindented value happens to work for
the changelog today while not actually being a trailer, which is the kind of
thing that holds until something else reads it.

Four things belong in it, and the last is the one a consumer actually acts on:

1. **which forms moved**, old and new;
2. **in which direction** — toward or away from the already-shipped form;
3. **roughly how many corpus rows**, with what you measured over (an order of
   magnitude is fine; a number with no corpus named is not);
4. **whether the affected inputs were previously rejected or previously
   accepted.** Previously *rejected* means no consumer has the old string
   stored, so the change is free. Previously *accepted* means someone does, and
   it is a migration. This distinction cannot be recovered later — nobody
   reading the diff in six months can tell — so it has to be written down now.

Two related habits, from `CLAUDE.md`:

- **"Fixes non-confluence" is not a sufficient description.** Confluence (two
  spellings of one variant agree) and stability (a variant normalizes to what it
  did last release) are different properties, and fixing the first does not give
  the second: every confluence fix picks a winner, and the losing side is a
  representation change for whoever stored it.
- **Break spec-permitted ties toward the form that already ships.** Free when you
  do; expensive when you do not. Where the spec mandates the other form, change
  it and flag the migration loudly.

Commits carrying the trailer are collected into a **Representation changes**
section at the top of the changelog, and that section is surfaced on the release
PR — see [Changelog](#changelog) below. The release PR states the section even
when it is empty, so that a cycle in which nothing was declared is visible rather
than silent.

Two limits of that section are worth knowing before you rely on it. It carries
the commit **subject** only — the trailer's own text does not reach the changelog
([#1556]), so the four facts above live in your PR description and a reader has
to follow the link for them. And a declining commit is filed under **Other**
whatever its type, so a `fix:` that correctly declined is not listed under
**Fixed** ([#1557]).

[#1556]: https://github.com/fulcrumgenomics/ferro-hgvs/issues/1556
[#1557]: https://github.com/fulcrumgenomics/ferro-hgvs/issues/1557

### Changelog

`CHANGELOG.md` is generated by [release-plz](https://release-plz.dev/) from the conventional-commit history when a release PR is opened — do not hand-edit it in feature PRs. Two PRs both adding bullets under `## [Unreleased]` will always conflict on the same line range, so leaving the changelog to tooling keeps feature PRs mergeable until release time. Write a clean conventional-commit subject (`feat(scope): …`, `fix(scope): …`, etc.) and release-plz will pick it up.

That it is generated is also why a `Representation-Change:` trailer is the way to get such a note into the changelog *at release time*: release-plz regenerates the pending section of an open release PR on every run, so an entry hand-written there is lost the next time the PR updates. `release-plz.toml`'s `commit_parsers` route trailered commits into a **Representation changes** group ahead of the usual ones, and its `pr_body` template surfaces that group on the release PR itself. The parsers match the commit **footer**, so the trailer must survive squashing — put it in the PR description as well, since that is what GitHub uses for the squash commit body.

A section that has already been released is a different matter: a release only prepends to `CHANGELOG.md`, leaving earlier sections untouched, and the GitHub release body is never revisited once its tag exists. So a retrospective note — a consolidated representation summary for a release whose commits predate this convention, say — can be added afterwards, by editing the released section and `gh release edit <tag>`. That is the escape hatch for a cycle that shipped without trailers; it is not a substitute for declaring changes as you make them, because nobody can reconstruct the rejected-vs-accepted distinction after the fact.

### Submitting a Pull Request

1. Push your branch and open a PR on GitHub
2. Fill out the PR template
3. Wait for CI to pass
4. Request review from maintainers

## Testing

```bash
cargo nextest run --features dev          # All tests (preferred)
cargo test --features dev                 # Alternative
cargo nextest run -E 'test(test_name)'    # Specific test
cargo test -- --nocapture                 # With output
cargo bench --features dev               # Benchmarks (`dev` is required by seqfirst_align)
```

### Conformance axis (manifest-backed)

`cargo nextest run --features dev -E 'test(axis_)'` selects the manifest-backed
conformance tests, but each of those silently returns early (reported as
PASSED, not skipped) when `FERRO_MANIFEST` is unset or points at a missing
file — so running that filter bare can look like a clean pass while testing
nothing. Use `scripts/run_conformance_axis.sh` instead: it validates the
manifest *before* invoking nextest, so a missing one is a hard failure rather
than a silent green. See `scripts/README.md` for details, including why the
test count the script prints is a weaker signal than the manifest check.

```bash
FERRO_MANIFEST=/path/to/manifest.json scripts/run_conformance_axis.sh
```

## Code Style

- Format with `cargo fmt`
- Lint with `cargo clippy --features dev -- -D warnings`
- Document all public APIs with doc comments

## Updating the v21.0 Normalization Fixture

`tests/fixtures/grammar/hgvs_spec_normalization.json` pins ferro's current
`normalize()` output for every variant string in the HGVS v21.0 spec
(vendored at `assets/hgvs-nomenclature/`). The companion test
`tests/it/hgvs_spec_normalization_tests.rs` fails any time a row's observed
output drifts from the recorded `current`.

If your PR changes any normalization output:

1. Make sure the vendored spec submodule is checked out, snapshot the current
   fixture, then regenerate it:

   ```bash
   git submodule update --init assets/hgvs-nomenclature
   # The fixture is a gitignored build artifact, so `git diff` will never show
   # your change. Snapshot the existing copy first — it is the baseline step 2
   # compares against. (On a fresh worktree there is nothing to snapshot;
   # regenerate on `main` first if you want a baseline.)
   cp tests/fixtures/grammar/hgvs_spec_normalization.json /tmp/spec-fixture-before.json
   cargo run --features dev --bin generate_spec_fixture
   ```

2. Inspect the change against that snapshot — not `git diff`, which shows
   nothing for an untracked file:

   ```bash
   diff -u /tmp/spec-fixture-before.json tests/fixtures/grammar/hgvs_spec_normalization.json
   ```

   For each changed row, verify the new `current` against the v21.0 spec text
   under `assets/hgvs-nomenclature/docs/recommendations/`. When ferro's new
   output now matches the spec's canonical form, that row should also have
   `current == spec_expected` after regen.

3. If a row's spec-canonical form differs from the input string (a "pair" — e.g.
   spec input `c.79GC>TT` is canonicalized to `c.79_80delinsTT`), record the
   divergence in `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`
   keyed on the input. Override entry shape:

   ```jsonc
   {
     "by_input": {
       "<exact input string>": {
         "status": "diverges",                       // optional: preserved | diverges | parse-error |
                                                     //           correctly-rejected | false-acceptance |
                                                     //           needs-reference
         "spec_expected": "<spec's canonical form>", // optional: string for canonical output, null for
                                                     //           "spec rejects this", absent for default
         "input_prefixed": "<accession:c.…>",        // optional: force a specific accession for bare
                                                     //           fragments (overrides default-prefix)
         "requires_reference": true,                 // optional: skip this row at test time until #82
                                                     //           lands (3'-rule shifting examples)
         "todo": "<https://… link>"                  // optional: defaults to a #83 link for any row that
                                                     //           lands as diverges / parse-error /
                                                     //           false-acceptance / needs-reference
       }
     }
   }
   ```

   Override keys must match a real fixture input — typos are caught by the
   generator.

### Status taxonomy

| status               | meaning                                                                             |
|----------------------|-------------------------------------------------------------------------------------|
| `preserved`          | ferro accepts the input and round-trips it (`current == spec_expected`)             |
| `diverges`           | ferro accepts the input but rewrites it (`current != spec_expected`)                |
| `correctly-rejected` | spec marks invalid (via `<code class="invalid">…</code>`), ferro also rejects       |
| `false-acceptance`   | spec marks invalid, **ferro accepts** — these are bug candidates                    |
| `parse-error`        | spec mentions the input as a canonical form, ferro can't parse it                   |
| `needs-reference`    | parse succeeds, normalization needs reference data ferro can't run today (see #82)  |

`spec_expected: null` carries the spec's "this string is invalid" intent. It
is set automatically for inputs the spec marks via `<code class="invalid">…</code>`
and can be set manually via the override file.

### Default-prefix table

The spec routinely writes bare fragments like `c.1083A>C` whose accession is
implied by the surrounding paragraph. The generator prepends a default
accession per coord system before feeding the input to ferro, recording the
result in `input_prefixed`:

| coord | default accession |
|-------|-------------------|
| `c`   | `NM_004006.2`     |
| `n`   | `NR_002196.1`     |
| `r`   | `NM_004006.3`     |
| `g`   | `NC_000023.11`    |
| `p`   | `NP_003997.1`     |
| `m`   | `NC_012920.1`     |
| `o`   | `NC_000023.11`    |

These are load-bearing constants that match the most-cited accessions in the
v21.0 spec. Re-validate them whenever bumping the spec submodule (see step 5).
For a row that needs a different accession, set `input_prefixed` in the
override file.

4. CI regenerates the fixture rather than diffing it, because the fixture is a
   gitignored build artifact with no committed baseline. The generation run is
   itself the guard: it harvests the spec checkout and resolves every curated
   override against it, so a stale override fails the build.

   ```bash
   cargo run --features dev --bin generate_spec_fixture
   ```

   The pre-push hook runs the same command for the same reason. `--check` is a
   local convenience answering a different question — "is my artifact current?"
   — and is never a gate:

   ```bash
   cargo run --features dev --bin generate_spec_fixture -- --check
   ```

   (Contrast the tool-support tables below, whose outputs *are* committed. There
   `--check` is the right gate, and CI uses it.)

5. To bump the spec to a newer upstream version:
   - update the submodule pointer under `assets/hgvs-nomenclature/`
     (`git -C assets/hgvs-nomenclature checkout <new-tag>`)
   - re-validate the default-prefix table above against the new spec corpus
   - regenerate the fixture and review the change against a snapshot (step 2 —
     the artifact is untracked, so `git diff` shows nothing)

## Tool-support comparison tables

The ferro/mutalyzer/biocommons/hgvs-rs support tables in `README.md`,
`docs/BENCHMARK_GUIDE.md`, and the web help tab are **generated**. To change a
cell, edit `docs/tool_support_matrix.json` and run:

```bash
cargo run --features dev --example generate_tool_support_tables
```

Do not edit the tables in those files directly — CI (`… -- --check`) will fail.
See `docs/tool_support_matrix.json` for the full schema and inline `_comment` field for authoring notes.

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
