# Testing ferro-hgvs

How to build and run the test suites, and the settings that make a local loop fast.
Moved out of `CLAUDE.md` on 2026-09-02; the content is as it stood there, and a
follow-on pass will edit it.

## Rust tests

```bash
cargo nextest run --features dev     # Run all tests (preferred)
cargo test --features dev            # Alternative
cargo nextest run -E 'test(parse)'   # Run specific tests by name pattern
```

**Those commands do not build `ferro-hgvs-soak-tests`, and they do not need to.**
That workspace member (`tests-soak/`) is the driver CI's optimized-archive jobs
run; it `#[path]`-includes the modules those jobs select, which remain under
`tests/it/` and remain
declared in `tests/it/main.rs`, so the commands above still run every one of
them, from the `it` binary. It is out of `default-members` precisely so nothing
above changes and no module is run twice locally. Build the driver itself only
when you have touched one of the modules it includes, a `tests/it/common/` helper they
reach, or the profile:

```bash
cargo nextest run -p ferro-hgvs-soak-tests
```

The one thing it exercises that `it` cannot: that binary's working directory is
`tests-soak/`, not the workspace root, so a helper resolving a fixture against
the cwd goes red there and passes in `it`. Fixture paths go through
`common::fixture_gen::fixture_path` for that reason.

## Python tests

```bash
pytest tests/python/ -v              # Run all Python tests
pytest tests/python/test_core.py -v  # Run specific test file
pytest -k "test_parse" -v            # Run tests matching pattern
```

## Fast local iteration

A one-line test edit costing two minutes is not the crate being big; it is two
settings. Measured on an M2 Max (12 core), interleaved A/B reps, `user` CPU
seconds quoted because wall clock on a shared machine is not reproducible:

| after editing… | default setup | both knobs | speedup |
|---|---|---|---|
| one `tests/it/*.rs` file | 16.8–18.0 s wall / 30.7–32.4 s CPU | 5.2–5.9 s / 4.0–4.4 s | **3.1× wall, 7.4× CPU** |
| `src/lib.rs` | 38.0–38.7 s wall / 133–136 s CPU | 9.5 s / 15.1–15.2 s | **4.0× wall, 8.9× CPU** |

**Knob 1 — do not export `CARGO_INCREMENTAL=0`.** This is the big one. If your
shell sets it (a common sccache incantation), incremental compilation is off, so
a one-line edit to any of the 428 files in `tests/it/` recompiles all ~139 k
lines of that crate from scratch. Confirm with `du -sh target/debug/incremental`
— `0B` means it is off, and rustc emitting 16 codegen units instead of 256 is
the same tell.

The premise the setting rests on does not hold for this crate anyway. **sccache
never caches any `ferro-hgvs` unit**, incremental or not: a rebuild's 35 compile
requests were classified `Non-cacheable calls`, reason **`crate-type`**, with
zero hits and zero misses (`sccache --show-stats` diffed either side of one
run). sccache is earning its keep on the ~400 dependency crates only, and those
are unaffected by this knob. So there is no cache benefit being traded away.

You cannot simply set it to `1`: sccache 0.16 **hard-errors**
(`sccache: incremental compilation is prohibited: Unset CARGO_INCREMENTAL to
continue.`) and the build dies on a dependency build script. **Unset** it
instead — cargo passes `-C incremental` only to workspace members, so sccache
never sees the flag on a dependency and keeps caching the whole dependency tree;
it soft-skips only this crate's own units, which are the ones you are
invalidating anyway.

```bash
env -u CARGO_INCREMENTAL CARGO_TARGET_DIR=$PWD/target-incr cargo t -E 'test(my_test)'
```

Two costs, both real: `target-incr/debug/incremental` reaches **2.0 GB**, and
flipping `CARGO_INCREMENTAL` **re-fingerprints every dependency** (measured: 318
crates recompiled), so the incremental and non-incremental builds cannot share
one `target/`. Use a separate `CARGO_TARGET_DIR` (`/target-*/` is gitignored) or
pick one mode and stay in it.

**Knob 2 — narrow the target set: `cargo t`.** A bare `cargo nextest run` builds
19 test binaries plus every `[[example]]` and `[[bin]]`, so a `src/` edit relinks
`ferro`, `ferro-web`, `ferro-benchmark`, both spec generators, ~20 examples and
11 example-test harnesses. The `cargo t` alias (`.cargo/config.toml`) is
`nextest run --features dev --lib --test it`, which halves the work
(37.3 → 18.1 CPU s on a `src/lib.rs` edit). It skips the two standalone
integration targets and the examples' own unit tests, so run `cargo ta` (the
full suite) before pushing.

**Measured and NOT worth doing** — recorded so nobody repeats them:

- `[profile.test] debug = 0` / `debug = "line-tables-only"`, with or without the
  same on `[profile.dev]` — no effect outside noise.
- `split-debuginfo = "unpacked"` — already the default here.
- **lld.** Not a win — consistently worse wall clock than Apple's `ld-1267`,
  which is already fast.
- **Splitting the `it` binary.** Unnecessary once knob 1 is on.

## Exhaustive cis sweeps: `FERRO_SWEEP_SEEDS`

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

The one thing to watch: `sweeps` is the only job that sets the variable, and it
selects on `SWEEP_FILTER`. Moving a sweep out of that filter silently drops it
to the prefix in CI, so edit the two together.

## Generated spec fixture (not committed)

`tests/fixtures/grammar/hgvs_spec_normalization.json` is a **generated build artifact** — it is `.gitignore`d, not committed. It is produced by the `generate_spec_fixture` binary from the HGVS spec submodule, the parser's behavior, and the curated `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`. (Committing it made every parser PR a merge-conflict magnet, since each PR regenerated the whole file.)

```bash
cargo run --features dev --bin generate_spec_fixture            # (re)generate it
cargo run --features dev --bin generate_spec_fixture -- --output <path>   # write elsewhere
```

Replace `<path>` with the destination you want; the first form writes to the default
location above. Do not substitute a machine-specific absolute path into a committed file
(see [Repository Hygiene](#no-machine-local-paths-in-committed-files)).

**Three** test modules actually read it — `tests/it/hgvs_spec_normalization_tests.rs` (the driver), `tests/it/idempotency_tests.rs` (which uses it purely as an input corpus, so its assertions are not circular) and `tests/it/normalize_axis_preserving.rs` (whose `axis_is_preserved_over_the_spec_fixture` walks it as a corpus; `ci.yml`'s `censuses` job downloads the `spec-fixtures` artifact precisely for it). They call a shared helper (`tests/it/common/spec_fixture.rs`) that regenerates it on demand if missing, so a fresh `cargo test` "just works" (one-time generation). A fourth module, `tests/it/spec_enumeration_tests.rs`, reads the *enumeration* fixture through `common/spec_enumeration.rs`, which ensures this one first as its generator's prerequisite — so **four** modules in total depend on a spec fixture being present. `mito_circular_audit.rs`, `protein_silent_eq.rs` and `protein_unknown_roundtrip.rs` only *mention* the fixture in doc comments as a cross-link; they do not read it. CI and the pre-push hook regenerate it explicitly before the test run — plain generation, not `--check`, because generation is what validates the committed overrides against the spec checkout.

**The fixture cannot detect a regression on its own, and neither can the enumeration's (#1272).** Both artifacts are generated from the code under test and regenerated by CI before the run, so `pinned_v21_normalization_behavior` and `enumeration_replays_recorded_behavior` compare ferro against itself — they are stale-local-artifact detectors. The guards that actually judge behaviour are the **committed** ones, and they are the ones to read when either drifts:

| guard | file | what it pins |
|---|---|---|
| `ferro_produces_the_form_the_spec_states` | `hgvs_spec_normalization_tests.rs` | ferro's output vs the **spec's** stated form, plus `KNOWN_DIVERGENT_INPUTS` |
| `status_census_is_unchanged` | `hgvs_spec_normalization_tests.rs` | per-status row counts |
| `spec_equivalence_classes_converge` | `hgvs_spec_normalization_tests.rs` | one normalized output per curated equivalence class, plus `EQUIVALENCE_CLASS_VERDICTS` |
| `ruling_records_are_intact` | `hgvs_spec_normalization_tests.rs` | the ruling records' ids and decided/undecided statuses (`RULING_STATUSES`) |
| `DIVERGENCE_BUDGET` / `PASSING_CENSUS` | `spec_enumeration_tests.rs` | per-status row counts, together totalling every row |

**Statuses only ask whether ferro *changed* a string — never whether two rows are the same
variant.** That is why the spec's own non-confluent pair
(`LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` and the `c.[850_869del;…]` split it names as that
variant's "alternative description", `DNA/delins.md:44-47`) sat in the fixture as two `preserved`
rows and passed. The `equivalence_classes` section of
`hgvs_spec_normalization_overrides.json` closes that: it declares which harvested inputs denote
one variant, and `spec_equivalence_classes_converge` fails when a class yields more than one
output. The disagreements it currently finds are **expected and pinned** — converging them is
#1235's job and a representation change for a downstream consumer, not something to fix in a test
PR. Members that cannot be evaluated hermetically (`needs-reference`, or a description ferro
refuses) are skipped rather than compared, so an error string never counts as a representation.

The sibling `rulings` section is the decision log for clauses that **conflict while pointing at the
same description**. (This used to say "conflict at equal RFC 2119 strength", citing `style.md:9`.
That framing is withdrawn — uppercase RFC 2119 keywords appear only **twice** outside `style.md`
itself, both in one file, so essentially every clause in play is lowercase prose and keyword
strength cannot rank them. The count is stated once, under
[Reading the spec](#reading-the-spec), with the grep it came from; do not restate it here.)
Each record names the clauses in tension, which
one governs, which is deviated from, and why. `undecided` is a first-class state and the generator
**refuses** to build a record that is `undecided` yet names a governing *or* a deviated-from clause
— either implies a side was chosen, and an unsettled question must not be able to smuggle in a
ruling nobody made. It equally refuses an `undecided` record citing fewer than two clauses: such a
record states no conflict, only a position it declines to name as one. Every citation carries a
quote that the generator checks against the spec checkout, so a submodule bump that moves a clause
fails the build instead of leaving the citation pointing at unrelated prose. That check is a
whitespace-collapsed substring match rather than a byte-for-byte one; what it does and does not
guarantee is spelled out under **Cite the clause exactly, and quote it** below.

Requires the `assets/hgvs-nomenclature` submodule (`git submodule update --init assets/hgvs-nomenclature`); without it the generator fails with `no HGVS strings harvested from …`, naming that command.

`--check` answers a different question — "is my local artifact current?" — and is not a gate: an absent artifact is generated rather than reported as drift, since a gitignored file has no committed baseline to drift from. Use it when you want to know whether a code change moved the fixture.

Because the fixture is no longer in git, per-PR `parse-error → preserved` status transitions are reviewed via the PR description and the accompanying test/parser changes, not a committed diff.

