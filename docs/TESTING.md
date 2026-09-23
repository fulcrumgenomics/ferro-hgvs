# Testing ferro-hgvs

This page tells you how to build and run the test suites, how the tests are organized, and how
to write one. It also tells you what CI runs and which settings keep a local edit-test loop fast.

## Rust tests

Run the suite with nextest:

```bash
cargo nextest run --features dev              # the whole suite
cargo nextest run -E 'test(parse)'            # one test, or a name pattern
cargo nextest run --features dev --no-capture # with test output
cargo bench --features dev                    # benchmarks; `seqfirst_align` requires `dev`
```

Use `cargo nextest`, not `cargo test`. `cargo test` runs each test target in its own binary, and
the test cases within a binary share that one process. Two lib tests, `normalize::merge::tests` and
`parallel::tests`, are not safe in that shared pool. nextest runs each test in its own process, so
they pass.

The one exception is doctests, which nextest cannot run: use `cargo test --doc --features dev`.

### Layout: one `it` binary

All integration tests are under `tests/it/`. They compile into one binary, `it`, not one binary
per file. The header of `tests/it/main.rs` gives the reason. Declare each new test file as a
`mod` in `tests/it/main.rs`. A file that is not declared does not compile and never runs. Name
filters still work, because nextest names each test `it::<module>::<test>`. Shared helpers are in
`tests/it/common/`. The `it` binary needs `--features dev`. Do not add another standalone Rust
test file directly under `tests/`.

Generate test data in the test. The committed corpora under `tests/fixtures/` are deliberate; do
not add new ones.

### Property tests and fuzz tests

- **Property test**: a [proptest](https://docs.rs/proptest) test in `tests/it/`. Name a new
  module `*_proptest` and also `#[path]`-include it in `tests-soak/tests/soak/main.rs`. CI leaves
  `*proptest*` modules out of the normal test shards and runs them only in its optimized soak
  job, which sees just the modules listed there; miss the second step and the test runs in no
  required job. Seeds of past failures are kept under `tests/proptest-regressions/`.
- **Fuzz test**: a [cargo-fuzz](https://rust-fuzz.github.io/book/cargo-fuzz.html) target in
  `fuzz/fuzz_targets/`. It needs the nightly toolchain and runs weekly in CI
  (`.github/workflows/fuzz.yml`), not in the normal suite.

Use each word only for its own kind.

### The `tests-soak/` member

`tests-soak/` is intentionally not in `default-members`. The default commands do not change, and
no module runs twice. The member includes its modules from `tests/it/` with `#[path]`, so the
`it` binary already runs all of them. You lose no coverage when you skip it.

The member adds one check. Its working directory is `tests-soak/`, not the workspace root, so a
helper that resolves a fixture path against the current directory fails there. The `it` binary
cannot see that failure. For that reason, fixture paths go through
`common::fixture_gen::fixture_path`. Build the member only when you change something it depends
on: a module it includes, a `tests/it/common/` helper that module uses, or its build profile.

```bash
cargo nextest run -p ferro-hgvs-soak-tests
```

## Suites that pass when their data is absent

**Bulk corpora.** Four large HGVS corpora are attached to the `test-fixtures-v1` GitHub release
rather than kept in the git tree. `scripts/fetch-test-fixtures.sh` downloads them and checks each
against `tests/fixtures/CHECKSUMS.sha256`; `tests/fixtures/README.md` covers provenance. Without
them, `clinvar_hgvs_tests`, `cmrg_exhaustive_tests`, `paraphase_exhaustive_tests` and
`normalize_axis_preserving` return early and report PASSED, not skipped.
`FERRO_REQUIRE_BULK_FIXTURES=1` turns that skip into a failure. CI sets it wherever it fetches
the corpora; leave it unset locally unless you have fetched them.

```bash
scripts/fetch-test-fixtures.sh --verify   # are they present and correct?
FERRO_REQUIRE_BULK_FIXTURES=1 cargo nextest run --features dev --lib --test it
```

**Manifest-backed axis tests.** The `axis_*` tests return early and report PASSED when
`FERRO_MANIFEST` is unset or points at a missing file, so a bare `-E 'test(axis_)'` filter can
pass while testing nothing. Run them through `scripts/run_conformance_axis.sh`, which validates
the manifest first. Build a manifest with `ferro prepare`; `scripts/README.md` has the details.

```bash
FERRO_MANIFEST=/path/to/manifest.json scripts/run_conformance_axis.sh
```

## Writing a test

**A test that pins a number is a change detector for that number.** It guards the property only
while the two agree, and nothing makes them agree.

- **Assert the property, and import the constant you guard.** `longest > 1024` is a number;
  `longest_block > <the cap the normalizer applies>` is a property. A comment telling the reader
  to update a literal when the constant moves is the defect, not a mitigation.
- **Prefer a unit test next to the item over widening the public API.** When an example or
  `tests/it` must reach an item, re-export that one item with `#[doc(hidden)]`, as
  `ShuffleDirection` is in `src/lib.rs`; never make a whole module public for it. Where a
  consumer needs a debug-only item, it refuses to run without it rather than reporting a zero,
  as `measure_spec_conformance_per_arm.rs` does for `partition_blocks_cut`.
- **If a count is the right assertion, say what it counts and against which denominator.**
- **Prove the test can fail.** Break the code once, watch the test go red, restore.
- **Before you quote a zero from a generated corpus, show the generator can produce the shape
  you changed.** Otherwise it is a structural zero: it describes the corpus, not your change.
- **One reproducer proves a defect exists, not its extent.** Measure the extent and scope the
  fix against it.
- **A passing case is evidence only once you know what made it pass.** A provider that rejects
  the input passes the test for you.

**Generators under `examples/` must account for what they drop.** A generator that writes an
artifact routes its records through `CaptureLedger` (`src/conformance/completeness.rs`), and a
generator with `#[cfg(test)]` sets `test = true` on its Cargo target.
`tests/it/generator_completeness.rs` enforces both, and its header explains what each check can
and cannot see.

## Python tests

Build the extension module first, then run pytest. The tests import the built module. If the
build is stale or missing, the tests fail:

```bash
uv run maturin develop --features python,extension-module   # build first
uv run pytest tests/python/                                  # then test
```

The Rust unit tests in `src/python.rs` run under nextest with `--features python` alone. Do not
add `extension-module` there: it stops the test binary linking libpython. The binary needs a
Python built with a shared library; the `Python Wheel Test` job in `ci.yml` shows the loader
setup. Scope the run and pass `--no-tests=fail`, so a filter typo fails instead of passing with
zero tests:

```bash
cargo nextest run --features python --lib --no-tests=fail -E 'test(python::tests::)'
```

## Fast local iteration

A one-line test edit should cost seconds, not minutes. Two settings decide that.

**Leave `CARGO_INCREMENTAL` unset.** Some sccache setups export `CARGO_INCREMENTAL=0`. With that
setting, a one-line edit to any file in `tests/it/` recompiles the whole crate. Do not set it to
`1` either, because sccache stops with an error on that value. Unset it:

```bash
env -u CARGO_INCREMENTAL cargo t -E 'test(my_test)'
```

Incremental and non-incremental builds cannot share one `target/`, because a change to the flag
re-fingerprints every dependency. If you use both modes, give each mode its own
`CARGO_TARGET_DIR`.

**Build less with `cargo t`.** A bare `cargo nextest run` also builds every binary and example.
The `cargo t` alias, `nextest run --features dev --lib --test it`, builds only the library and
the `it` suite. That is about half the work. It skips the standalone integration targets, so run
the full suite with `cargo ta` before you push.

Measured on an M2 Max with the repo-pinned toolchain (`rust-toolchain.toml`), warm `target/`, no
other load: after a one-line edit to a `tests/it/` module, rebuilding a single test with
`env -u CARGO_INCREMENTAL cargo t -E 'test(<one_test>)'` fell from about 17 s to about 5 s wall, and
from about 30 s to about 4 s CPU, once both settings above are in place. The numbers are approximate
and vary with the edit and the machine. These did not help, so you do not need to try them: `lld`, a
split of the `it` binary, and `[profile.test]` debug settings.

## Exhaustive cis sweeps: `FERRO_SWEEP_SEEDS`

Three exhaustive sweeps draw a large deterministic corpus of sequences. By default they run a
4-seed prefix of that corpus. Set `FERRO_SWEEP_SEEDS` when you need more:

```bash
cargo nextest run --features dev                        # 4-seed prefix (default, fast)
FERRO_SWEEP_SEEDS=full cargo nextest run --features dev # the full corpus, as CI runs it
FERRO_SWEEP_SEEDS=12   cargo nextest run --features dev # a fixed count, for bisecting
```

If you shrink a sweep, cut sequence diversity, never shapes. The prefix does exactly that, and
that axis is safe to cut. Every blocking defect these sweeps found was in a shape the generator
could not emit, never in a sequence it did not draw.

Only the `sweeps` CI job sets `FERRO_SWEEP_SEEDS=full`, and it selects tests with `SWEEP_FILTER`.
A sweep that reads the variable but is not in that filter runs at the prefix everywhere, and CI
stays green. When you add a sweep, name it in `SWEEP_FILTER` (so `test` and `test-oracle` negate
it) and in the `sweeps` job's selection, or it runs at the prefix everywhere with nothing to
flag it.

## What CI runs

The required `Test` context is a rollup job, `test-required` in `.github/workflows/ci.yml`. It
needs every test job to succeed. If it is red while every shard is green, read the rollup's log
to find the upstream job. Read its `needs:` list from the file, not from this page. The jobs, and
how to run each one locally:

- `test`: the default suite, with `FERRO_REQUIRE_BULK_FIXTURES=1`. `test-oracle`, `censuses`, and
  the fetching jobs in `coverage.yml` and `external-validation.yml` set it too. Some suites skip
  and report PASS when their bulk data is absent. That variable turns the skip into a failure.
  Set it locally when you have the data. `scripts/run_conformance_axis.sh` runs one
  manifest-backed axis.
- `test-oracle`: the suite with the four `FERRO_ASSERT_*` seam oracles armed, over the suite minus
  the modules those oracles would silence. See `docs/ORACLES.md`, section "What CI arms, and
  where".
- `sweeps`: the three exhaustive sweeps at `FERRO_SWEEP_SEEDS=full`.
- `censuses` and `censuses-plain`: the slow census modules, from the optimized soak archive.
  `censuses-plain` runs the modules that build their corpus in code. Those modules refuse to run
  when an oracle is set.
- `soak`: the idempotency property tests, the `*proptest*` modules, from the same archive.
- `hgvs-rs-tests`: the unit tests behind the `hgvs-rs` feature. No other job enables that
  feature.

## Generated spec fixture (not committed)

`tests/fixtures/grammar/hgvs_spec_normalization.json` is a generated artifact, not a committed
file. It is gitignored. When it was committed, every parser PR caused merge conflicts in it.
Regenerate it with the `generate_spec_fixture` binary, which reads the `assets/hgvs-nomenclature`
submodule:

```bash
git submodule update --init assets/hgvs-nomenclature                     # once, if missing
cargo run --features dev --bin generate_spec_fixture                      # (re)generate
cargo run --features dev --bin generate_spec_fixture -- --output <path>   # write elsewhere
```

Never commit the file. Never put a machine-specific output path in a committed file.

**A replay test cannot catch a regression. It compares ferro against its own output.** The
fixture is generated from the code under test, and CI and the pre-push hook regenerate it before
every run. A failing replay test is therefore a stale local artifact, not a regression. The tests
that judge behaviour are the committed guards in `hgvs_spec_normalization_tests.rs` and
`spec_enumeration_tests.rs`. When a replay test drifts, read the headers of those guards to see
what they pin.

**A status records whether ferro rewrote a string, not whether two rows mean the same variant.**
For example, the spec gives `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` and a `c.[850_869del;…]`
split as two descriptions of one variant. They sit in the fixture as two untouched rows, and both
pass.

Identity is a separate assertion. The `equivalence_classes` section of
`hgvs_spec_normalization_overrides.json` declares which inputs denote one variant. Its check
fails if a class gives more than one output. The disagreements it finds today are expected and
pinned. To converge them is a downstream representation change, not a fix for a test PR. The
`rulings` section of the same file is described in "Recording a normalization decision" below.

Regenerate the fixture with plain generation, not `--check`. Generation is what validates the
committed overrides against the spec checkout, and CI does the same before every run. `--check`
answers only one question: is my local artifact current? It is not a gate, and CI does not gate
on it.

### Updating the fixture when your PR changes normalization output

1. **Snapshot the fixture, then regenerate it.** The fixture is gitignored, so `git diff` shows
   nothing and the snapshot is your only baseline. On a fresh worktree there is no fixture to
   snapshot; regenerate it on `main` first and snapshot that.

   ```bash
   cp tests/fixtures/grammar/hgvs_spec_normalization.json /tmp/spec-fixture-before.json
   cargo run --features dev --bin generate_spec_fixture
   ```

2. **Diff against the snapshot**, and check each changed row against the spec text under
   `assets/hgvs-nomenclature/docs/recommendations/`. A row whose new output matches the spec's
   canonical form shows `current == spec_expected`.

   ```bash
   diff -u /tmp/spec-fixture-before.json tests/fixtures/grammar/hgvs_spec_normalization.json
   ```

3. **Tell the generator where the spec's canonical form differs from the input; otherwise it
   expects the input unchanged.** For a row such as `c.79GC>TT`, which the spec canonicalizes to
   `c.79_80delinsTT`, add an entry to `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`
   keyed on `c.79GC>TT` with `spec_expected` set to `c.79_80delinsTT`. A key that matches no
   fixture input fails generation.

4. **To bump the spec:** move the submodule pointer with
   `git -C assets/hgvs-nomenclature checkout <new-tag>`, re-validate the default accessions in
   `prefix::DEFAULTS` (`examples/common/spec_harvest.rs`) against the new spec, then regenerate
   and review against a snapshot as in steps 1 and 2.

An override entry has this shape:

```jsonc
{
  "by_input": {
    "<exact input string>": {
      "status": "diverges",                 // optional; see the table below
      "spec_expected": "<canonical form>",  // optional; null means the spec rejects the input
      "input_prefixed": "<accession:c.…>",  // optional; accession to force for a bare fragment
      "requires_reference": true,           // optional; skip at test time, needs reference bases
      "note": "<why>"                       // optional; why this override exists
    }
  }
}
```

Bare fragments such as `c.1083A>C` get a default accession from `prefix::DEFAULTS` in
`examples/common/spec_harvest.rs`, recorded in the row's `input_prefixed`. Each row gets one of
these statuses, assigned by `classify` in the same file:

| status               | meaning                                                                             |
|----------------------|-------------------------------------------------------------------------------------|
| `preserved`          | ferro accepts the input and round-trips it (`current == spec_expected`)             |
| `diverges`           | ferro accepts the input but rewrites it (`current != spec_expected`)                |
| `correctly-rejected` | spec marks invalid (via `<code class="invalid">…</code>`), ferro also rejects       |
| `false-acceptance`   | spec marks invalid, **ferro accepts**. These are bug candidates                     |
| `parse-error`        | spec mentions the input as a canonical form, ferro cannot parse it                  |
| `needs-reference`    | parse succeeds, normalization needs reference data the test cannot supply           |

The generator sets `spec_expected: null` itself for inputs the spec marks
`<code class="invalid">…</code>`.

## Recording a normalization decision

**Record every decision about the *correct* normalization output in a committed test or ruling
record, in the same PR.** The record states the question, the ruling and the authority. A
decision that lives only in a PR description, an issue comment or a working document is lost,
and the next person re-derives it differently. The rules themselves are in
[`docs/src/reference/normalization-rules.md`](src/reference/normalization-rules.md).

A decision needs a record when it is a ruling that one clause governs another, a determination
that ferro's output is right or wrong against a cited clause, a decision to follow or deviate
from Mutalyzer, a choice between two competing representations of one variant, or a question
deliberately left open. Implementation choices, refactors and performance work need none.

| the decision is about | record it in |
|---|---|
| two spec clauses in tension | a `rulings` record in `hgvs_spec_normalization_overrides.json` |
| two spellings that must converge on one output | an `equivalence_classes` entry + `EQUIVALENCE_CLASS_VERDICTS` |
| a deliberate, known deviation from the spec's stated form | `KNOWN_DIVERGENT_INPUTS`, which pins it *as* a deviation |
| a deliberate deviation from Mutalyzer where the spec is silent | a `rulings` record. There is no spec form to diverge *from*, so `KNOWN_DIVERGENT_INPUTS` is the wrong home |
| a concrete input whose correct output is now settled | an ordinary `tests/it/*` test pinning the exact string |

All of these live under `tests/fixtures/grammar/` or `tests/it/`. For the shape of a `rulings`
record, read `adjudication-precedence-order` in the ledger; it also records why Mutalyzer is not
a spec oracle and how the spec, Mutalyzer and house choices rank.

**State which kind of record it is.** The kinds are not interchangeable:

- **adjudicated-correct**: pin the exact expected output and cite the clause.
- **adjudicated-deviation**: pin it in `KNOWN_DIVERGENT_INPUTS`, so a fixed deviation cannot
  stay in the list unnoticed.
- **undecided**: better than no record. It names no governing or deviated-from clause and cites
  at least two clauses in conflict; the generator refuses one that does not.
- **house-choice**: decided, and ours rather than the spec's (rule 5's silent limb or rule 6). It
  names no governing clause, is never citable as conformance, and says what was rejected.

Every record except an `undecided` one is `"status": "decided"` and carries a one-sentence
`summary`.

**A test that pins today's output is not an adjudication record.** It is a change detector.
What makes a record an adjudication is its authority: an exact `file:line` into
`assets/hgvs-nomenclature`, a named Mutalyzer measurement, or an explicit "undecided, and here is
why". Without one, you have frozen the current behaviour, including whatever is wrong with it.
A deviation from Mutalyzer needs a record too, or the next person measures Mutalyzer, finds the
mismatch, and "fixes" a decision.

**Record what was refuted, not only what was decided.** A measurement that kills a plausible
belief is worth as much as the ruling, because the belief will recur. Put it in the ruling
record, the issue or the PR, not in a code comment; a test's doc comment may name the belief it
guards against in one line.

### Citing the spec

Cite the clause exactly, and quote it: `general.md:33`, not "the separation rule". A clause's
directory is its jurisdiction: a claim about `r.` needs a clause under `RNA/`, because a `DNA/`
clause cannot scope `r.`. `docs/READING_THE_SPEC.md` is the guide to reading the spec.

A green citation check does not prove a quote is exact. It is a whitespace-collapsed substring
match, so a quote that differs from the spec in spacing or line breaks passes. If a claim rests
on exactness, compare against the spec file.

### After you edit the ledger

Five things are built from or checked against the ledger. Rewording a record's prose affects the
first two; changing its cited clauses affects the first three; adding or removing a record, or
changing its status, affects all five.

1. `docs/NORMALIZATION_CONTRACT.md`. Its failing test prints the bless command.
2. The shadow-spec `why` blocks and `tests/fixtures/shadow_spec/corpus.jsonl`. Their failing test
   prints the bless command. A bless updates only the generated `why` block, so check the
   hand-written text beside it.
3. The clause index in `tests/it/clause_ruling_index.rs`. There is no bless variable: the failing
   `the_rendered_index_is_current` prints the replacement block. Save it to a file and paste it
   between the `BEGIN`/`END` markers rather than retyping it. For a record count or status
   change, also update the `(records, decided, undecided)` count in `the_index_is_not_vacuous`.
4. `RULING_STATUSES` in `tests/it/hgvs_spec_normalization_tests.rs`, edited by hand and **never**
   regenerated, so that a record cannot be added, removed or decided without an edit a reviewer
   sees.
5. The gitignored `tests/fixtures/grammar/hgvs_spec_normalization.json`: regenerate it (see
   "Generated spec fixture" above).

Commit the updated files by explicit path, then re-run the tests without any bless variable.

## Generated documentation

Some committed docs are generated, and CI fails when a hand edit makes them disagree with their
source:

- **Tool-support tables** in `docs/src/reference/comparison.md`, `docs/BENCHMARK_GUIDE.md` and
  `src/service/web/static/data/tool_support_matrix.json`. Edit `docs/tool_support_matrix.json`
  and run `cargo run --features dev --example generate_tool_support_tables`. The JSON file's
  `_comment` field describes the schema.
- **`docs/NORMALIZATION_CONTRACT.md` and the shadow-spec `why` blocks**, both built from the
  ruling ledger. "After you edit the ledger" above lists what to update.
