# Testing ferro-hgvs

This page tells you how to build and run the test suites and how the tests are organized. It
also tells you what CI runs and which settings keep a local edit-test loop fast.

## Rust tests

Run the suite with nextest:

```bash
cargo nextest run --features dev              # the whole suite
cargo nextest run -E 'test(parse)'            # one test, or a name pattern
```

Use `cargo nextest`, not `cargo test`. `cargo test` runs all tests as threads in one process.
Two lib tests, `normalize::merge::tests` and `parallel::tests`, are not safe in that shared pool.
nextest runs each test in its own process, so they pass.

### Layout: one `it` binary

All integration tests are under `tests/it/`. They compile into one binary, `it`, not one binary
per file. The header of `tests/it/main.rs` gives the reason. Declare each new test file as a
`mod` in `tests/it/main.rs`. A file that is not declared does not compile and never runs. Name
filters still work, because nextest names each test `it::<module>::<test>`. Shared helpers are in
`tests/it/common/`.

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

## Python tests

Build the extension module first, then run pytest. The tests import the built module. If the
build is stale or missing, the tests fail:

```bash
maturin develop --features python,extension-module   # build first
pytest tests/python/                                  # then test
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

On an M2 Max, the two settings together cut a rebuild after a `tests/it/` edit from about 17 s to
about 5 s. CPU time fell from about 30 s to about 4 s. These did not help, so you do not need to
try them: `lld`, a split of the `it` binary, and `[profile.test]` debug settings.

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
stays green. `tests/it/sweep_filter_invariant.rs` checks that every module that reads the
variable is in the filter, and the reverse. When you add a sweep, that test tells you what to
edit.

## What CI runs

The required `Test` context is a rollup job, `test-required` in `.github/workflows/ci.yml`. It
needs every test job to succeed. If it is red while every shard is green, read the rollup's log
to find the upstream job. Read its `needs:` list from the file, not from this page. The jobs, and
how to run each one locally:

- `test`: the default suite, with `FERRO_REQUIRE_BULK_FIXTURES=1`. Some suites skip and report
  PASS when their bulk data is absent. That variable turns the skip into a failure. Set it
  locally when you have the data. `scripts/run_conformance_axis.sh` runs one manifest-backed
  axis.
- `test-oracle`: the suite with `FERRO_ASSERT_IDEMPOTENT`, `FERRO_ASSERT_REPARSE` and
  `FERRO_ASSERT_IN_BOUNDS` set, without the spec-corpus census modules that those oracles would
  silence. `FERRO_ASSERT_IDEMPOTENT=1` on a bare `cargo nextest run` is always red on `main`.
  Run `scripts/run_oracle_suite.sh` instead. It uses the same selection and flags as the job.
  `tests/it/oracle_exclude_invariant.rs` keeps the exclusion list and the job the same.
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
`rulings` section of the same file is documented in `CONTRIBUTING.md`, Adjudications.

Regenerate the fixture with plain generation, not `--check`. Generation is what validates the
committed overrides against the spec checkout, and CI does the same before every run. `--check`
answers only one question: is my local artifact current? It is not a gate, and CI does not gate
on it.
