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

#### Fast local iteration — two knobs, both measured

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
  same on `[profile.dev]`. Three-way interleaved, no effect outside noise (4.95 /
  5.08 / 4.49 CPU s). On macOS `split-debuginfo` already defaults to `unpacked`,
  so debuginfo lives in the `.o` files and is never copied into the binary — the
  linker only writes a debug map. `line-tables-only` did shrink objects 270 MB →
  164 MB and bought nothing. It is also not free to try: it re-fingerprints the
  whole dependency tree, which is what `ci.yml` means when it says
  `CARGO_PROFILE_TEST_DEBUG=0` "re-fingerprints every dependency" (confirmed:
  318 crates).
- `split-debuginfo = "unpacked"` — already the default here. No `.dSYM` bundle is
  produced, so there is no `dsymutil` step to remove.
- **lld.** `ld64.lld` ships inside the rustc sysroot (no Homebrew needed) and does
  link this crate on macOS arm64, via
  `RUSTFLAGS="-C link-arg=-fuse-ld=$(rustc --print sysroot)/lib/rustlib/aarch64-apple-darwin/bin/gcc-ld/ld64.lld"`.
  It is **not** a win: ~14 % less CPU but consistently *worse* wall clock than
  Apple's `ld-1267`, which is already fast.
- **Splitting the `it` binary.** Unnecessary once knob 1 is on: with incremental
  compilation, editing one of 428 test modules costs 3.9 s of cargo time across a
  139 k-line crate. A split would add a second full link for every change to
  `tests/it/common/`.

#### Running the oracles locally — a bare armed run CANNOT pass

Read this before running any of the four flags below over the whole suite. The
obvious command is documented nowhere else in this file any more, because it is
a trap:

```bash
FERRO_ASSERT_IDEMPOTENT=1 cargo nextest run --features dev   # ALWAYS RED on main
```

It fails **7** tests on a clean `origin/main`, and has for as long as the spec
corpus has existed. Use the local runner instead, which mirrors `ci.yml`'s
`test-oracle` job — its flags, and its selection:

```bash
scripts/run_oracle_suite.sh                     # as test-oracle runs it
scripts/run_oracle_suite.sh --print-selection   # what it would run, without running it
scripts/run_oracle_suite.sh -E 'test(my_test)'  # extra args reach nextest
```

**Why the bare form is red, and why that is not a coverage gap to close.** All 7
failures come from `FERRO_ASSERT_IDEMPOTENT` alone — measured: `FERRO_ASSERT_REPARSE`
and `FERRO_ASSERT_IN_BOUNDS` each contribute **0**, and a prepared reference
changes nothing (the same 7 fail with and without `FERRO_MANIFEST`; every one of
these modules is `MockProvider`-backed). All 7 live in modules `ORACLE_EXCLUDE`
already names, for two reasons that are worth keeping apart:

- **5 of them ASSERT the defect the oracle PANICS on.** The four
  `defect_non_idempotent_outputs` tests and
  `spec_corpus_regressions::an_insertion_at_the_cds_end_is_not_a_fixed_point`
  pin `c.*1delinsCTT` → `c.72_*1insCT` → `c.72delinsCCT` as a known non-fixed
  point. A test that pins a defect and an oracle that aborts on it cannot both
  run; there is no version of this the flag could pass.
- **2 of them COUNT it**, and arming the oracle makes the count read *better*
  than the truth. `spec_conformance_axis`'s censuses wrap normalization in
  `catch_unwind`, so a panicking row is filed `declined` and never reaches its
  family's output set. See `ORACLE_EXCLUDE`'s comment in `ci.yml` for the
  measured figures.

So the answer is not "run it in CI anyway" — that is the one thing that would
destroy the evidence. The answer is that CI already excludes it, and the local
path now excludes it the same way, from the same source of truth.

**The runner reads the whole `-E` selection and the flag set out of `ci.yml`**,
so neither can drift into a second copy. "Whole" is load-bearing: `test-oracle`
negates `test(proptest)` and `SWEEP_FILTER` as well as `ORACLE_EXCLUDE`, and a
runner that negated only the last of them would run the proptest modules and the
three exhaustive sweeps while reporting itself as mirroring that job. So the
shape of the expression is read too, and only its variable references are
expanded — rebuilding it here from a hardcoded `not test(proptest) and not (…)`
would trade today's drift for tomorrow's.

`tests/it/oracle_exclude_invariant.rs` re-derives all of it in Rust — different
anchors, deliberately — and compares; separately asserts that the agreed
selection actually negates all three (a spelling-independent check, since both
sides otherwise read the same line of `ci.yml` and a matching pair of wrong ones
would satisfy the equality); and separately forbids the script from hardcoding a
module name. Each of those guards is checked against a deliberate sabotage
rather than assumed.

#### Normalization idempotency oracle

`FERRO_ASSERT_IDEMPOTENT=1` turns every `Normalizer::normalize` call into an
assertion that `norm(norm(x)) == norm(x)`, so any test that normalizes becomes an
idempotency check. Run it through `scripts/run_oracle_suite.sh` — see
[Running the oracles locally](#running-the-oracles-locally--a-bare-armed-run-cannot-pass)
for why the bare full-suite form cannot pass.

It is compiled out in release builds (`#[cfg(debug_assertions)]`) and read once
into a `OnceLock`, so a disabled run pays only one atomic load. CI runs the suite
a second time with it set; the nightly reference-aware job sets it too, which is
where the manifest-backed conformance corpora are actually covered.

The check runs from `Normalizer::assert_seam_oracles`, which all four oracles
share, at the single exit of `normalize_core_checked`. That covers every public
normalization: `normalize()`, `normalize_with_diagnostics()`, and every
`VariantProjector` path (the projection-driven genomic/coding/protein axes).

It is one call site again as of #1382. `normalize_with_diagnostics` used to reach
`normalize_core_canonical` directly, which bypassed both the oracles and — the
actual defect — the strict-mode rejection ladder, so #1366 had to call
`assert_seam_oracles` from it separately. Routing it through
`normalize_core_checked` fixes both at once, and the extra call would now just
re-run all four oracles on that path.

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
- A **non-coding downstream position** (`n.*N`), which #1748 refuses at parse in every
  mode while `TxPos::downstream` stays public API — so the spelling is gone and the AST
  is not. Reachable only by a Rust library caller constructing one; no string entry
  point reaches it, and normalization only ever preserves the flag, never mints it.
  Keyed on the AST through `noncoding_zone_marker`, which is `Tx`-only by an exhaustive
  match, so `c.*N` — anchored to the CDS, still legal — is not exempted.

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

`scripts/run_oracle_suite.sh` arms it, alongside the two flags beside it — the
bare full-suite form is red for a reason this oracle has no part in, so it tells
you nothing about re-parsing. Measured: `FERRO_ASSERT_REPARSE` on its own fails
**0** tests over `ORACLE_EXCLUDE`'s modules.

Kept separate from `FERRO_ASSERT_IDEMPOTENT` because neither subsumes the other,
and the idempotency oracle has a blind spot this one covers: it verifies by
*re-normalizing* its own output, which it cannot do for an output that fails to
parse, so an unparseable result is invisible to it. Both run from
`assert_seam_oracles` alongside the in-bounds and denoted-sequence checks, and CI
sets these two and the in-bounds check together. The **denoted-sequence** flag is
the exception and is not set everywhere they are — see
[Where it runs, and the one place it deliberately does not](#normalization-denoted-sequence-oracle).

#### Normalization in-bounds oracle

`FERRO_ASSERT_IN_BOUNDS=1` is the third at the same seam: no coordinate a
normalized description names may be past the end of its own sequence.

`scripts/run_oracle_suite.sh` arms it too. As with the re-parse oracle,
`FERRO_ASSERT_IN_BOUNDS` on its own fails **0** tests over `ORACLE_EXCLUDE`'s
modules — the 7 failures a bare armed run shows are the idempotency oracle's,
and reading them as this one's would send you after the wrong invariant.

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
read once into a `OnceLock`. CI sets **these three** together, in the sharded
`test-oracle` job and the `sweeps` job (those two and nowhere else — the plain
`test` and `soak` jobs run without the flags); the nightly sets it too, where it
is the only place the check runs against true transcript and contig lengths.

**The fourth flag is not set in all of those places, so do not read "all four"
off this paragraph.** `FERRO_ASSERT_SEQUENCE` runs in `sweeps` and in the
nightly, and deliberately not in `test-oracle` — the reason, and what currently
keeps it out, are under
[Where it runs, and the one place it deliberately does not](#normalization-denoted-sequence-oracle)
below and in `ci.yml`'s comment on the `test-oracle` job.

#### Normalization denoted-sequence oracle

`FERRO_ASSERT_SEQUENCE=1` is the fourth at the same seam, and the only one that
asks what the output **means** rather than how it is written: apply the input to
the reference, apply the output, and assert the bases agree.

```bash
FERRO_ASSERT_SEQUENCE=1 cargo nextest run --features dev -E "$SWEEP_SELECTION"
```

**`scripts/run_oracle_suite.sh` deliberately does NOT arm this one**, because
`test-oracle` does not — the runner mirrors that job, and mirroring it faithfully
means inheriting the gap. The job where this flag is green is `sweeps`, so scope
a local run to that job's selection rather than to the whole suite. Added to
`ORACLE_EXCLUDE`'s modules it raises **5** further failures beyond the
idempotency oracle's 7, all in `spec_corpus_regressions` and all the same shape
as those: a test pinning the CDS-end defect that this oracle aborts on.

**That 5 is about `ORACLE_EXCLUDE`'s modules, and there is a second, unrelated 5
below.** `test-oracle`'s own selection *excludes* those modules, and arming the
flag over it raises 5 failures too — none of them in `spec_corpus_regressions`
and none of them the same defect. The two figures share a digit and nothing
else, so never carry one over to the other's scope; say which selection a count
was taken on, every time.

**The other three are all form questions, and a wrong sequence passes all of
them.** It is a fixed point, so `FERRO_ASSERT_IDEMPOTENT` is satisfied; it
parses, and `parse_hgvs` holds no provider so it could not know either way; and
its coordinates exist, so `FERRO_ASSERT_IN_BOUNDS` is satisfied. #1592 and #1600
each record all three passing on their own reproducer. The class had been found
by hand six times before those two — **#1254** (`g.[3_4del;9del]` → `g.12_14del`),
**#1281** (→ `g.[1del;1del]`), **#1290**, **#1304**, **#1308**, **#1312** — each
filed, fixed and regression-tested separately. Eight instances is the same
argument #1353 made for the in-bounds oracle after three.

**The applier is not the normalizer.** `spdi::compare_denoted_sequences` reaches
the bases through `hgvs_to_spdi` and an SPDI splice — the same walk
`apply_to_reference` and `tests/it/common/cis_apply_oracle.rs` use — so nothing
here can agree with the output merely because normalization produced it.
`EquivalenceChecker` is **not** usable for this: it normalizes both sides, which
is circular (`tests/it/issue_1234_sibling_clamped_shift.rs:198`).

Both descriptions are applied over the **union** of their spans, in one fetch. A
per-description window would give each its own frame and report every 3'-shift as
a difference; over a window containing both, `g.3_4del` and `g.7_8del` in a tract
denote the same bases, which is what makes shifting sequence-preserving in the
first place.

**A side that cannot be applied is counted, not silently passed** — a skip that
reads as a pass is the exact failure mode this oracle exists to remove:

| case | verdict |
|---|---|
| both apply, bases agree | pass, counted in `compared` |
| both apply, bases differ | **fire** |
| the **output** denotes no sequence while the input does | **fire** — #1281's `g.[1del;1del]`, two members claiming one base. Worse than a wrong sequence, so it must not be filed under the skip exemption |
| the **input** denotes no sequence (trans allele, `REFSEQ_MISMATCH`, an edit SPDI cannot carry) | skip, counted in `skipped` — there is no baseline, the same discipline the other two apply to already-broken inputs |
| the two name different accessions (#785 version substitution) | skip, counted |
| the union window exceeds `MAX_APPLY_WINDOW`, or the provider cannot serve it | skip, counted |

`normalize::denoted_sequence_oracle_counts()` returns `(compared, skipped)`
process-wide, so a run can say how much it actually compared. **Read it before
trusting a green oracle run**: zero comparisons and zero faults look identical
from the outside.

`tests/it/issue_1615_denoted_sequence_oracle.rs` is the oracle's own regression
guard. It pins the recorded `(reference, input, wrong output)` triple from each
of the eight issues and asserts the predicate fires — deliberately *not* by
re-normalizing, since a test that ran the normalizer would go green as each
defect is fixed and stop saying anything about the oracle. Its other half is the
negative control: legitimate re-spellings (a 3'-shift, a merge, a decomposition,
a dup/ins equivalence) must stay silent.

Compiled out in release like the three beside it, and the most expensive of the
four when on — one provider fetch and two splices per normalization — so it runs
last in `assert_seam_oracles`, after the cheaper and more specific checks have
had their chance to name the fault more precisely.

**Where it runs, and the one place it deliberately does not.** `sweeps` sets it
(gating, measured green over the full corpus) and the nightly sets it
(non-gating). `test-oracle` does **not**, which makes it the only job where the
set of four is incomplete. Two rows in that job's selection used to fire, and
both were real disagreements inside ferro rather than noise:

| row | what disagreed | state |
|---|---|---|
| #1618 — `NC_TEST.1:g.262TG[6]` → `g.259_262GT[6]` | `hgvs_to_spdi` read the anchored spelling as 6 copies replacing a **1**-copy tract, the normalizer's output as 6 copies replacing a **2**-copy tract — 14 bases against 12 | closed before `6116f84a` |
| #1619 — `NM_033517.1:c.4818dupC` → `c.4818dup` | `hgvs_to_spdi` resolved the `c.` position by **walking** the exon list while the normalizer indexes the **flat** transcript, so the two disagreed across any transcript-coordinate gap: the input applied `C`, the output `T`, at transcript position 4877. `NM_033517.1` carries a real 39-base cdot hole between exons 10 and 11 — see below, because this row **replaced** an earlier one on the same issue | closed by the flat-frame fix: `cds_to_tx`/`tx_to_cds` no longer read the exon table |

**Both are now green, and the flag is STILL not set here — the selection-wide
run those two closures were blocking on has since been done, and it is RED.**
The two-row measurement, on the #1619 branch, was this:

```bash
FERRO_ASSERT_SEQUENCE=1 cargo nextest run --features dev --test it \
  -E 'test(test_explicit_base_removal)'          # 2 passed
FERRO_ASSERT_SEQUENCE=1 cargo nextest run --features dev --test it \
  -E 'test(repeat_tract_maximization)'           # 3 passed
```

That is the two named rows, not `test-oracle`'s selection, which is the whole
suite minus `ORACLE_EXCLUDE`, `test(proptest)` and `SWEEP_FILTER`. Run over
**that** selection — `origin/main` @ 674e9c8b, this job's own env plus
`FERRO_ASSERT_SEQUENCE=1`, no `--partition` — it reads
`10383 tests run: 10378 passed, 5 failed, 289 skipped`. The five are in two
classes, **neither of them #1618 or #1619**, so neither was reachable by the
two-row check above:

- **3x `issue_1487_canonical_window_overflow`** — not an oracle fire at all.
  `attempt to add with overflow` at `src/convert/mapper.rs:114`: arming the flag
  routes the output through `hgvs_to_spdi`, which reaches `cds_to_tx`'s unchecked
  `cds_start as i64 + pos.base - 1` on an `i64::MAX`-adjacent `c.` coordinate.
  The open issue is **#1690**.
- **2x `stranded_identity_member`** — a genuine fire, of #1281's "denotes no
  sequence" shape, on a module that exists to PIN a defect (#1655's stranded
  derived identity). That is the class `ORACLE_EXCLUDE` already documents: a test
  pinning a defect and an oracle aborting on it cannot both run.

So arming it is blocked on **#1690** and on finding a home for
`stranded_identity_member` — **not** on another measurement. Do not read "both
rows are green" as "it can be armed". `ci.yml`'s comment on the `test-oracle`
job carries the full triage, including what a run with `FERRO_MANIFEST` set adds.
Suppressing a row would still hollow out the oracle; that has not changed.

**#1619's row was replaced before it was closed, and the swap is the whole
point.** The
old row (`NM_000492.4:c.1520_1522del` → `c.1521_1523del`) rode on a **fixture
defect**, not on annotation: `normalization_transcripts.json` was the real cdot
exon table with every exon's `end` decremented by one, so all 22 of its
multi-exon records carried a one-base hole at every junction — a shape that
occurs nowhere in 474,818 real builds. Repairing the fixture removed that row.

The row that stands in its place rides on **cdot's own annotation**.
`NM_033517.1` is one of the 58 gapped builds, with a genuine 39-base hole
between exon 10 (ends 1302) and exon 11 (starts 1342), and the fixture now
carries its real table. So #1619 is demonstrated on real transcript geometry
rather than on an artifact — which is a *stronger* reason to keep
`FERRO_ASSERT_SEQUENCE` out of `test-oracle`, and it hands the issue the
reproducer the repair would otherwise have destroyed. The gap is exempted by
name in `CDOT_GAP_JUNCTIONS`
(`tests/it/normalization_transcripts_exon_contract.rs`), never by loosening the
contiguity rule.

**No CI job arms the oracle where it fires** — checked, not assumed.
`FERRO_ASSERT_SEQUENCE` is set in exactly two places: `ci.yml`'s `sweeps`, which
selects `SWEEP_FILTER + test(issue_1615_denoted_sequence_oracle)`, and
`nightly-mutalyzer.yml`, which selects every reference-aware module and is
`continue-on-error`. `normalize_tests` is in neither selection, so the fire is
reachable only by setting the flag by hand. (That nightly selection named
**three** modules until the orphaned-guard sweep; `normalize_tests` was not one
of the modules added, so this paragraph's conclusion is unchanged — but do not
re-derive "the three modules" from it.)

**Here is the command the before/after comes from**, because a figure with no
command is one nobody can re-derive. Swap only the fixture, and exclude the guard
this change adds — that guard is *about* the fixture, so it cannot be a term in a
before/after on it:

```bash
# "before" = the fixture as it stands on the PR's base; "after" = the fixed one
git show <base>:tests/fixtures/sequences/normalization_transcripts.json \
  > tests/fixtures/sequences/normalization_transcripts.json
FERRO_ASSERT_IDEMPOTENT=1 FERRO_ASSERT_REPARSE=1 FERRO_ASSERT_IN_BOUNDS=1 \
FERRO_ASSERT_SEQUENCE=1 \
  cargo nextest run --features dev --lib --test it \
  -E 'not test(normalization_transcripts_exon_contract)'
```

Against `origin/main` at `d5f26fcb`: **25 failures before, 22 after, zero new.**
The three that went green are `idempotency_tests::test_normalization_idempotency`,
`normalize_tests::normalization_transformations::test_3prime_shifting::case_1` and
`normalize_tests::normalization_transformations::test_deletion_shifts_in_real_homopolymer::case_1`.

**That is two files, not three.** Three files read the fixture —
`normalize_property_tests`, `normalize_tests` and `idempotency_tests` — and the
three green tests come from the latter two; `normalize_property_tests` reads it
and contributes none. This paragraph used to call them "exactly the three suites
that read that fixture", which is the claim being corrected.

**Quote the command, not the pair.** The *difference* is stable — three tests,
always the same three — but the totals are a property of `main` on the day: the
same measurement read **20 before, 17 after** against `5f22abed` and 25/22
against `d5f26fcb` a few hours later, because unrelated PRs were landing into the
pre-existing oracle residue the whole time. A bare "25 / 22" quoted without its
base is not reproducible, which is what this paragraph originally shipped.
`tests/it/normalization_transcripts_exon_contract.rs` now guards the shape.

**The defect the issue named was live until the flat-frame fix**, and losing the
only corpus that exercised it was the hazard. `hgvs_to_spdi` resolved a `c.`
position by *walking* the exon list while the normalizer indexes the *flat*
transcript, and those two readings disagree whenever the exon table has a
transcript-coordinate gap. Real cdot has such gaps — measured over
cdot-0.2.32.refseq.GRCh38, **58** of 474,818 multi-exon builds, sizes 23–2718
bases (none of them one base, which is why the fixture's shape was diagnostic of
a generator bug rather than of real annotation).

**And the suite does exercise it, on real geometry.** `NM_033517.1`'s restored
cdot table is that corpus: one record, one genuine gap, one existing case
(`c.4818dupC`) that used to fire the denoted-sequence oracle. It is why the
issue could be closed against measured geometry instead of a corpus built from
scratch. **It is still the regression guard, so keep it that way** — do not
"simplify" it by re-flattening the record or by widening `CDOT_GAP_JUNCTIONS`:
the first destroys the reproducer, the second re-admits the synthetic one-base
holes the contract guard exists to keep out.

**Which side was wrong, and on what authority — the record is
`c-and-n-positions-are-flat-transcript-offsets` in
`hgvs_spec_normalization_overrides.json`.** In one line: the accession's own
numbering is flat (`NP_277052.1` is 1731 aa, GenBank annotates
`NM_033517.1` `CDS 1..5196`, and RefSeq's own exon table for it is contiguous),
so cdot's hole is a fact about its **genome alignment**, not about the
coordinate space a `c.` number is written in. The exon walk left
`CoordinateMapper`; genome↔transcript mapping is untouched and stays exon- and
CIGAR-aware. The second half of the defect — cdot's `start_codon`/`stop_codon`
being gap-collapsed, so the live provider serves `NM_033517.1` `cds_end = 5157`
where RefSeq says `5196` — is **not** fixed and is not what that record decides.

**Measured false-positive classes, and what closed each.** The first run of this
oracle over the suite raised **344** fires, all but 16 of them false. Each fix
below is a class, not a case, and the reasoning is on the code.

**Do not sum the `fires` column, and do not read it as a partition of the 344.**
Only the 344 and the 16 come from one run. Each row's count was measured on the
run that isolated *that* class, after the rows above it had been closed — closing
one class changes which normalizations reach the comparison at all, so the same
normalization can be counted in more than one row and the column adds up to well
over 344. Read each figure as the size of the class at the point it was diagnosed,
which is what makes it an argument for the fix beside it:

| class | fires | why it was not a defect |
|---|---|---|
| output cannot be transliterated | 328 | the input states its own deleted bases (`g.1000delC`) and converts with no provider; the output (`g.1000del`) must read the reference, which the fixture does not hold |
| insertion flush against a deletion | 233 | `triples_are_disjoint` called it an overlap; `cis_apply_oracle`'s tie-break does not, and it is well defined. **No longer true of the applier as of #1831** — `apply_triples` now sorts longer deletions first among triples sharing a position, the key `splice_denoted_sequence` and `conformance::spec_corpus`'s applier already used, so `triples_are_disjoint` sees correctly-ordered input and no longer reports this shape. The 233 is the size of the class when it was diagnosed, and is kept as that |
| overlap-conflicting input | ~12 | an insertion *interior* to a deletion — the input denotes nothing, so there is no baseline |
| `pter`/`qter` | 6 | carry no numeric coordinate (`base: 0`); `hgvs_to_spdi` reads `base` and silently resolves `c.pterdel` to the sequence's **last** base |
| corrected `REFSEQ_MISMATCH` | 5 | `c.10dupA` where the reference reads `C`: normalization is *supposed* to change the denoted sequence here |
| `r.` payload vs DNA reference | 1 | `UGC` against `TGC` — the same bases in two alphabets |
| uncertain allele `[(…)]` | 1 | the normalizer deliberately does not clamp those |

The lesson generalizes past this oracle: **the two sides of a comparison do not
need the reference equally**, and any check that treats "I could not derive it"
as "it is wrong" will fire hardest on the shapes where that asymmetry is largest.

**An oracle failure is visible in PR CI and silent in the nightly.** In PR CI,
`test-oracle` and `sweeps` carry no `continue-on-error`, so a fire turns the job
red. The nightly reference-aware job does carry it, deliberately — its purpose is
to surface drift in the xfail report rather than to gate, and the corpus runner
wraps normalization in `catch_panics`, so an oracle panic there lands in the
uploaded xfail artifact as a FAILing case instead of failing the workflow. That
applies equally to all four flags, and it is why a red oracle in the nightly must
be read out of the xfail report rather than from the workflow conclusion.

**So the nightly's `report-failure` issue is never about an oracle.** `if:
failure()` cannot fire on a step whose failure `continue-on-error` has already
absorbed, so every issue that job opens is an infrastructure one — checkout, the
`ferro` build, `ferro prepare`, a missing manifest, the artifact upload — and its
`hint:` says so. The armed reproduction command lives where the oracle failure
actually surfaces instead: an `outcome == 'failure'` step writes it to the run's
**job summary**, alongside the pointer to the `ferro-xfail-*` artifact.

**An oracle fire DOES block the merge — corrected 2026-08-09.** This section used
to say the opposite: that the ruleset requires only a subset of the jobs CI runs,
that neither `Test oracle` nor `Exhaustive sweeps` is among them, and so a fire
"shows up as a failed job that a human must not merge past, not as a
ruleset-enforced merge block." The premise is right and the conclusion does not
follow, because it reads the ruleset's context names as job names.

`Test` is not the sharded test job. It is the `test-required` job in
`.github/workflows/ci.yml` — cited by name because a line number goes stale on
the next insertion above it — a **rollup**:

```yaml
  test-required:
    name: Test
    needs: [test, test-oracle, soak, sweeps]
```

whose script echoes each result and exits 1 if any is not `success`. So the
required `Test` context transitively gates on `test-oracle`, `soak` and `sweeps`.

**`generated-docs` is NOT in that list, and still gates** — it is one hop further
out, so do not read its absence as a gap. The doc generators run inside
`spec-fixtures`, and `test`/`test-oracle` both `needs:` that job, so a generator
failure *skips* them; `test-required` demands `result == "success"` from each and
a skip is not a success, so the required `Test` context still goes red. The
comment at `ci.yml`'s "Gating is unchanged" step records the move — and it is the
reason to read the `needs:` array from the workflow rather than from prose, since
a job named there can be retired without the gate it provided going away.

Both halves of the claim are checkable from the tree rather than from prose, and
should be re-checked rather than trusted — the `needs:` list above and the
required-context list below have each changed at least once already:

```bash
# is `Test` really a rollup, and over which jobs?
grep -A4 '^  test-required:' .github/workflows/ci.yml

# is `Test` really required on main?  (the set grows — `Representation change
# declared` and `zizmor + actionlint` were both added after this was written,
# so any count quoted in prose is stale by the next contexts change)
gh api --paginate repos/fulcrumgenomics/ferro-hgvs/rulesets \
  --jq '.[] | select(.target=="branch") | .id' \
| xargs -I{} gh api repos/fulcrumgenomics/ferro-hgvs/rulesets/{} \
  --jq 'select(.enforcement=="active")
        | select(.conditions.ref_name.include
                 | any(. == "refs/heads/main" or . == "~DEFAULT_BRANCH"))
        | .rules[] | select(.type=="required_status_checks")
        | .parameters.required_status_checks[].context'
```

Two practical consequences. A red oracle is a hard merge block, not an advisory
one. And `gh pr checks` showing `Test` red while every `Test (n/6)` row is green
is not a contradiction — look at the rollup's log to find which upstream job
actually failed, rather than re-running the shards. Note the second consequence
is the *mechanism*, not a claim that any particular PR is in that state: an
earlier draft of this paragraph cited PR #1572 as having six green shards and two
red oracle shards, and by 2026-08-10 that PR read three red shards and three red
oracle shards. Re-running CI moves the anecdote; it does not move the rollup.

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

### Declaring a representation change: a REQUIRED check, not a convention

A PR touching `src/normalize/`, `src/hgvs/`, `src/spdi/`, `src/project/`, `src/reference/`
or `src/error_handling/` must carry a `Representation-Change:` trailer **in its PR
description**, or the `Representation change declared` check fails. That context is in
`main`'s ruleset, so a missing trailer **blocks the merge**.

The last two were added by #1853 on measurement rather than on reasoning, and the reasons
they are the *only* two added are on `WATCHED_PREFIXES` in
`scripts/check_representation_change.py` — read that docstring before proposing a seventh.
In short: the v0.13.0 cycle is the one place output movement was measured without relying
on anyone declaring it, and its 326,404 newly-normalizing rows attribute to one PR touching
only `src/reference/` and one touching only `src/error_handling/`. `src/conformance/` stays
out because it cannot move ferro's output, and `src/data/` stays out because every measured
disclosure that touches it also touches `src/reference/`.

```
Representation-Change: 577 rows move, 360 merge / 205 split / 12 respell.
  Previously-accepted inputs, so a real migration for the consumer.
```

**`Representation-Change: none` is a first-class answer, not an opt-out.** `none`, `no`, `n/a` and
`na` all count, in any casing, and a declining trailer is excluded from the changelog's
**Representation changes** section — so declaring it costs the reader nothing. Most changes under
those directories move no output at all (a comment, a doc, a new unit test), and that is the
expected case. What the check rejects is *silence*, because an absent trailer is
indistinguishable from an unconsidered one.

**The verdict is the first word, and a reason may follow it** — `none. Tests only; no watched
file is touched.` A full stop, semicolon, colon or dash introduces the reason; a comma does not,
because it usually introduces a qualification that changes the answer ("none, except two rows"),
which is filed as a real change. This was the #1555 defect: the value used to be anchored, so
declining *with a reason* was read as a disclosure. It is not a corner case — in the v0.13.1
cycle 8 of the 14 declines gave a reason, and the section listed 10 entries of which 2 were real.

**One limit of the section itself, and one that is now closed.** A declining commit is filed
under **Other** whatever its type, so a declining `fix:` is missing from **Fixed** (#1557);
git-cliff evaluates a parser's fields as OR, not AND, so this cannot be repaired by adding
`message` to the exclusion rule.

The other — the bullet carrying the commit subject and nothing else, so the
rejected-vs-accepted fact had to be read out of every linked PR — is #1556, and it is closed by
`scripts/inject_representation_disclosure.py`, which attaches each trailer's own text under its
bullet. **It is deliberately not a git-cliff template**, and the three template routes measured
dead in #1556 should not be re-walked: `commit.footers` is mis-parsed by git_conventional (any
`word:` line is a footer token, so a prose body shreds it), `commit.body` has the trailer split
out of it entirely (0 occurrences for both real v0.13.1 disclosures), and splitting on the
unanchored literal finds only the prose *mentions* — which is worse than nothing, because it
produces plausible output for some commits and quoted prose for others.

The route that works is not through the template: read the trailer from the raw commit message
with `git log`, and attach it to the rendered bullet afterwards — which is what
`check_changelog_grouping.py` already does to attribute trailers, so the script imports its
`trailer_value` rather than keeping a fourth copy of the column-0 rule.

Two things it measured that contradict what `CONTRIBUTING.md` asks for, and that any future work
here will hit again. **The trailer is not the last thing in the squash body**: GitHub appends
CodeRabbit's auto-generated summary after it, so reading to end-of-message quotes three
paragraphs of "Summary by CodeRabbit" as the author's disclosure. And **the trailers are not
well-formed**, their continuation lines being unindented, so indentation cannot delimit the
value — which is also why `git interpret-trailers --parse` reports *no* trailers at all for
#1537, #1535 or #1547. The script stops at the next column-0 trailer token, Markdown heading or
HTML comment instead.

CI runs it as `--check` in the `Changelog grouping audit` job: a no-op on an ordinary PR, and the
gate on the release PR, where the fix is the one command it names. It is idempotent, so re-run it
after any push that makes release-plz regenerate the pending section.

**It must be in the PR description, not only in a commit message.** GitHub builds the squash
commit body from the description, and `release-plz.toml`'s parsers match the commit *footer*; a
trailer on an intermediate commit is discarded by the squash and never reaches the changelog.

**Why this is enforced rather than encouraged.** The mechanism to route these trailers into the
changelog landed mid-cycle and was then used **zero** times: across the v0.13.0 cycle 87 commits
landed carrying no trailer while 46 touched those four directories, so the changelog's
Representation changes section rendered empty. An empty section looks identical whether it is
empty because nothing moved or because nobody declared anything, so it went unnoticed for ~90
commits and the disclosure had to be reconstructed after the fact — normalizing ClinVar, Paraphase
and CMRG (5,761,302 expressions) through both releases and diffing row by row, which found 577
stored strings moved and 326,404 inputs newly normalized. None of that was visible from the
changelog as generated.

**How to measure, if the answer is not obviously `none`.** `examples/dump_normalized_corpus.rs`
diffs two revisions over a synthetic shape-family corpus (`--out` on each side, then `--compare`);
read its `--verify-spdi` report to separate "moved" from "denotes different bases". That corpus is
deliberately enriched for churn-prone shapes, so quote its rate as *of the affected family*, never
as a repo-wide figure — for consumer impact, normalize a real corpus through both revisions with
`ferro normalize -i <inputs> --reference <dir> -f tsv` instead.

**A decline that then describes a move is rejected.** `no. 3 rows move` reads as a decline —
first word wins — so the disclosure would vanish silently, which is the worst direction this
mechanism can fail in. The checker fails a declining trailer that also claims `<n> rows
move|merge|split|respell` for non-zero `n`. Quantifying a *zero* is encouraged and passes
(`none. 0 of 950 rows move`). **The number that counts is the count, never the denominator** —
`3 rows of 500,004 move`, `3 of 500,004 rows move` and `3 of 500,004 corpus rows move` all disclose
three moved rows, and the rule reads three from each. It is a tripwire for the phrasing this repo
actually uses, not a general contradiction detector — the enumerated forms live on
`MOVEMENT_CLAIM_RE`, and that docstring is the one to extend when a new phrasing appears.

This sentence used to say "the pattern requires the count to sit immediately before `rows`", and
that was the #1647 defect rather than a description of it: in `0 of 950 rows move` the number
immediately before `rows` is the **denominator**, so the very example offered here as passing was
rejected, and the documented form was a merge blocker. Do not restate the rule positionally.

**The section itself is audited, not just each PR's trailer.** `Changelog grouping audit` renders
the changelog with the real `release-plz.toml` over `<latest tag>..HEAD` and checks the result
against the commits: nothing filed under **Representation changes** may open with a decline word,
nothing that declares a move may be filed elsewhere, and no heading may still carry its `<!-- N -->`
ordering prefix. It deliberately shares **no code** with the checker — #1555 passed the test written
to catch it precisely because that test compared the two halves against each other and they were
wrong together, so a second opinion is only worth having when it is derived differently.

**"Derived differently" means the derivation, not the verdict — and reading it the other way shipped
a defect.** The audit's rule was originally made deliberately *coarser* than the checker's (first
word, punctuation stripped, no terminator logic) on that reasoning. The result was a rule that
disagreed with `release-plz.toml` on trailer forms `CONTRIBUTING.md` documents as **correct**:
`no rows move` and `none, except two rows that merge` are filed as real changes by design, and the
audit called them declines, so **no trailer text could satisfy both halves** — and once such a
commit is on `main` the job is red for every open PR until the next release tag. The reverse also
bit: `none.Tests only.` (a missing space) split to `none.Tests`, which the checker called a decline
and the audit called a move. The two now agree on the verdict, pinned by
`test_the_audit_and_the_checker_agree_on_every_documented_form`; what stays independent is that
this check *renders the changelog through real git-cliff and reads the result*, which is the
question no vocabulary-comparison test was asking. A disagreement is not a second opinion, it is an
unsatisfiable build.

Contributor-facing guidance is in `CONTRIBUTING.md`. The checker is
`scripts/check_representation_change.py` and the audit is
`scripts/check_changelog_grouping.py`; the decline vocabulary is pinned against
`release-plz.toml` and `CONTRIBUTING.md` by tests in
`tests/python/test_representation_change_trailer.py`, because three copies of one list drift
silently.

### Never hand-edit `tests/it/clause_ruling_index.rs` — it embeds a GENERATED block

The file is committed, so it looks hand-maintainable. It is not: it carries a rendered index of
every clause the ruling records cite, and `the_rendered_index_is_current` compares that block
against what the code generates. Editing it by hand is transcription, and transcription of a
~100-row block is lossy in a way that reads as a small diff — two attempts at it were made here
and both dropped rows silently.

Capture what the test prints **programmatically** instead. The pattern that works: add a throwaway
test that writes `render_index()` to a file, run it, restore the source byte-for-byte from a
pre-edit copy, then splice the generated text in on its `BEGIN`/`END` markers. Nothing is retyped.

Note this is a *different* trap from the gitignored spec artifacts above. Those regenerate on
demand and a stale one fails loudly. This one is committed, so a lossy edit survives review as an
ordinary diff and only fails on the one assertion that re-renders it.

### Record every adjudication as a committed test or ruling record, in the same change

**Policy.** When you decide what the *correct* normalization behaviour is — not how to
implement it — that decision lands in a committed test or ruling record **in the same PR**, with
a comment stating the question, the ruling, and the authority. An adjudication that lives only in
a PR description, an issue comment, or a working document is lost: the next person re-derives it
from scratch, and often re-derives it differently.

This is cheap because the machinery already exists (the five committed guards are tabulated
directly above). The policy is about using it consistently, not about building anything.

**What counts as an adjudication:** a ruling that one clause governs another where they conflict;
a determination that ferro's output is right or wrong against a cited clause; a decision to follow
or deviate from Mutalyzer; a choice between two competing representations of one variant; or a
question deliberately left open. Implementation choices, refactors and performance work are not
adjudications and do not need records.

**Where it goes:**

| the adjudication is about | record it in |
|---|---|
| two spec clauses in tension | a `rulings` record in `hgvs_spec_normalization_overrides.json` |
| two spellings that must converge on one output | an `equivalence_classes` entry + `EQUIVALENCE_CLASS_VERDICTS` |
| a deliberate, known deviation from the spec's stated form | `KNOWN_DIVERGENT_INPUTS`, which pins it *as* a deviation |
| a deliberate deviation from Mutalyzer where the spec is silent | a `rulings` record — there is no spec form to diverge *from*, so `KNOWN_DIVERGENT_INPUTS` is the wrong home (see `adjudication-precedence-order`) |
| a concrete input whose correct output is now settled | an ordinary `tests/it/*` test pinning the exact string |

**State which kind of record it is, because they are not interchangeable:**

- **adjudicated-correct** — pin the exact expected output and cite the clause. This is a real
  guard: it fails when behaviour regresses away from a decided answer.
- **adjudicated-deviation** — pin it as a deviation so it stays visible. `KNOWN_DIVERGENT_INPUTS`
  exists for this, and `ferro_produces_the_form_the_spec_states` deliberately fails when a listed
  input *starts* matching the spec, so a fixed deviation cannot rot in the list unnoticed.
- **undecided** — a first-class state. The generator refuses an `undecided` record that names a
  governing clause, so an open question cannot smuggle in a ruling nobody made. Prefer an honest
  `undecided` record to no record.

**A test that merely pins today's output is not an adjudication record.** It is a change
detector, and this repo has already been bitten by the difference —
`pinned_v21_normalization_behavior` compares ferro against itself (see the note above about
stale-local-artifact detectors). What makes a record an adjudication is the *authority*: an exact
`file:line` into `assets/hgvs-nomenclature`, a named Mutalyzer measurement, or an explicit
"undecided, and here is why". Without one, you have frozen the current behaviour including
whatever is wrong with it.

**Cite the clause exactly, and quote it.** `rulings` citations carry a quote the generator checks
against the spec checkout, so a submodule bump that moves a clause fails the build instead of
leaving the citation pointing at unrelated prose. Do the same in prose comments: `general.md:34`,
not "the separation rule". A clause's directory is its jurisdiction: a claim about an `r.` axis
needs a clause under `RNA/` because a `DNA/` one cannot scope `r.`, while `general.md` is not
molecule-specific.

**That check is a whitespace-collapsed substring match, not a byte-for-byte one.** The generator
joins the cited line range with spaces and collapses runs of whitespace on both sides before
testing containment. Two consequences follow, and only the first is usually noticed. A quote may
**span** the cited lines — which is what makes a multi-line range work at all — and re-wrapping the
spec's prose leaves a citation valid while the clause moving out from under it still fails the
build. That much is the deliberate trade. What the guard does **not** buy is that a quote is
reproduced byte-for-byte: one whose only difference from the spec is spacing or a line break
passes, and citations in this ledger do exactly that. So do not offer the guard as evidence that a
quote is exact — if a claim rests on exactness, measure it against the spec file rather than
inferring it from a green build.

**Record what was refuted, not only what was decided.** Measurements that killed a plausible
belief are worth as much as the ruling itself, because the belief will recur. Existing examples
worth imitating: `MIN_SEPARATION_NO_FRAME`'s doc comment explains why applying the coding
exception everywhere is wrong; `apply_coding_codon_exception` records the worked
`c.10_13delinsTCAG` case showing why the rule tests the triplet rather than the hull;
`split_codon_incompatible_triplets` names the two tests that pin both sides of its boundary.

**Deviating from the reference implementation needs a record, not just a rationale.** Precedence
is **spec-explicit > Mutalyzer > our judgement**. Where the spec is explicit and Mutalyzer differs,
that is a deliberate divergence and it gets a record saying so — otherwise the next person
measures Mutalyzer, finds a mismatch, and "fixes" a conformance decision.

### Assert the property. Measure the count. Never let a count BE the property

The section below is the mechanized half of this one, and this one had never been written down —
so it kept being re-derived, one instance at a time, by whoever it had just cost a day.

**A guard that pins a number is a change detector for that number.** It is a guard for the
property only for as long as the two agree, and nothing makes them agree. When they drift, the
guard follows the number, stays green, and reads as coverage.

Three shapes, in ascending order of how convincing the wrong answer looks.

**1. A count restated instead of imported.** `the_corpus_emits_a_block_past_the_split_cap`
(`examples/dump_normalized_corpus.rs:2381`) opens `const SPLIT_CAP: u64 = 1024;` and asserts the
corpus builds a block wider than it. `MAX_SPLIT_BLOCK` is **4096** (`src/normalize/merge.rs:2000`,
derived from `MAX_CANONICAL_WINDOW`), raised from `1024` by #1899. The corpus's long cores are
`[1024, 1100]` and were chosen — the comment says so — as "the pair a change to the cap moves
between". The cap moved, both cores landed on the same side, and `1100 > 1024` kept the guard
green. It even carries a comment instructing the reader to update the literal if the constant
moves; **that comment is the defect, not a mitigation.** A guard that restates the value it guards
cannot observe that value changing. Import the constant. If it is private, make it `pub(crate)` —
that is cheaper than the class of failure it buys out of. Filed as #1925.

**2. A zero whose instrument cannot vary the property.** `examples/dump_normalized_corpus.rs`
has now been blind six times — geometry (#1456), scale (#1460), sequence structure (#1517),
transcript geometry (#1478), molecule type (#1606), and reversed ranges (#1917) — and its own
header records the governing fact: **fixing one blindness does not reveal the next.** Each time,
a `0` was available to quote as safety. So before quoting one, **name the property your change
keys on and show the generator can vary it.** A zero you cannot attribute to the change is a
claim about the instrument. Report a structural zero *as* structural, in those words.

**3. A rule generalised from the one case you reproduced.** #1917 reported an underflow and stated
the rule as "any reversed range whose span exceeds `window_size`". Measured, the failure is three
bands in `start - end` against the window: `<= w` fails somewhere else entirely (a backwards
slice), `(w, 2w]` is the reported underflow, and `> 2w` never failed at all — the fetch bounds
themselves invert and the provider rejects the read. The rule was wrong in **both** directions,
and the issue's own two real-corpus rows were in the band that never panicked. One reproducer
establishes that a defect exists. It does not establish its extent, and the extent is what a fix
is scoped against.

**Two more ways a number arrives already wrong**, both seen here:

- **Right for a reason that is not the code.** #1917's widest band passed before the fix *and*
  after it — because `JsonProvider` rejects an inverted read. That is a property of the installed
  provider, not of the thing under test. A passing case is only evidence once you know what made
  it pass.
- **Two confirmations sharing one stale source are one observation.** Two agents "independently
  confirmed" a set of `merge.rs` line numbers; both had read `main/`'s working tree, ~950 lines
  stale. Independence is a property of the derivation, not of the number of people who reported it.

**So, when you write a guard:**

- Assert the **property**, in the words you would use to explain it. `longest_block > <the cap the
  normalizer actually applies>` is a property; `longest > 1024` is a number.
- If a count is genuinely the right assertion — a census, a ledger cardinality — say what it counts
  and against which denominator, and pin it somewhere a change to the thing it counts must touch.
- **Prove the guard can fail.** Sabotage it once, watch it go red, restore. An assertion never
  observed failing is indistinguishable from one that cannot. The `ORACLE_EXCLUDE` invariants and
  `issue_1615_denoted_sequence_oracle` are both built this way; imitate them.
- When you report the number, report what it is made of. "46 checks" and "42 pass, 1 skipping, 3
  cancelled corpses" describe one run, and only the second is usable.

### A generator must account for what it dropped

This is the mechanized half of **"a corpus zero is a claim about the corpus, not about the
change"** — the doctrine that a generator's `0` reads as *"this is safe"* when it usually means
*"the corpus could not build the thing you changed"*. That doctrine covers a zero the generator
*reports*. The same trap one level down is a population the generator **silently drops**: a
fallible step whose failure is representable as a legitimate value —
`spans_of(..).unwrap_or_default()` classifying as "no gaps", `else { continue }` skipping exactly
the record that had a problem, a discarded `Result`. Nothing counts it, so a partial run and a
clean run write indistinguishable artifacts and the only thing separating them is whether a
reviewer read the right five lines.

**Route the population through `CaptureLedger`** (`src/conformance/completeness.rs`, #1550):

```rust
use ferro_hgvs::conformance::completeness::{Allowance, CaptureLedger};

let mut ledger = CaptureLedger::new("transcript windows");
for accession in &accessions {
    // Replaces `let Ok(t) = resolve(a) else { continue };` — same control flow,
    // but the drop is now accounted for.
    // `record`'s subject is `impl Into<String>`, which `&String` does not
    // satisfy and does not auto-deref to — hence `as_str()`.
    let Some(transcript) = ledger.record(accession.as_str(), resolve(accession)) else {
        continue;
    };
    capture(transcript);
}
// `finish` refuses on any shortfall, and on a pass that attempted nothing.
// Waiving one means naming it: `finish_with(Allowance::at_most(1, "why"))`.
let counts = ledger.finish()?;
```

Three things follow. **Only one of them is enforced by a test; the other two are conventions, and
saying otherwise here would be the exact defect this section is about** — a claim that reads as
coverage while providing none:

- **A generator carrying `#[cfg(test)]` must set `test = true` on its cargo target.** Enforced by
  `it::generator_completeness::examples_with_unit_tests_opt_into_running_them`, whose two halves
  are not equally strong: the manifest is *parsed* as TOML, so `test = true` is read rather than
  guessed, but `#[cfg(test)]` is found by scanning the compiled sources for that literal token —
  which matches it inside a doc comment and misses `#[cfg(all(test, feature = "dev"))]`. Say
  "parsed manifest, token scan", not "exact". Cargo does not build a target's tests unless it
  opts in, so without the flag they never run — and a committed test that has never executed is
  worse than none, because it reads as coverage. `report_conformance_reference_gaps` shipped one
  for months.
- **Refuse, do not report.** A convention, backed by the *type* rather than by a guard over your
  code. `finish()` returns a `Result` you must handle before reaching your `fs::write`, and
  `CaptureCounts` has no `Default`, so `finish().unwrap_or_default()` does not compile. Nothing
  stops `let _ = ledger.finish();`.
- **Stamp the counts into the artifact.** A convention, enforced by **nothing**. `CaptureCounts`
  is serde-serializable so a fixture *can* carry its own completeness claim and a consuming test
  *can* assert on it — no test anywhere checks that any artifact does.

`tests/it/generator_completeness.rs` holds the `test = true` invariant above and a second guard:
a generator that looks like it writes an artifact must look like it routes its population through
the ledger, or name itself in `LEDGER_EXEMPT`. **Read that second one as a floor, not a proof.**
Both of its sides are substring scans over source text — `fs::write(`, `File::create(`,
`OpenOptions::new(`, `from_path(` on one side; on the other, a `use` line naming
`conformance::completeness`, the word `CaptureLedger`, and a `finish` call, all three in the
generator's own entry point. A writer using an idiom outside that list is invisible to it, and
"routes its population through" is a semantic question a grep cannot answer. What it does buy: a
generator written the way every generator in this repo is currently written cannot be added in
silence, and the question gets asked. The semantic guarantee is carried by the ledger and by
review.

The allowlist is **shrink-only** — a row whose generator has been deleted, stopped matching the
write idioms, or now imports the ledger and calls `finish` fails the test — so it cannot rot into
a blanket exemption. The "now uses the ledger" direction is deliberately hard to trip by accident,
because acting on it *deletes* a row: an earlier revision keyed on the bare word `CaptureLedger`
anywhere in the compiled sources, and one comment appended to the shared
`examples/common/recording.rs` marked three unrelated generators as migrated and demanded their
rows be removed. Adding a row is a legitimate answer, in the same way
`Representation-Change: none` is; what neither check tolerates is silence.

**And note what the ledger itself does not claim.** It accounts for the step you route through
it, not for the artifact: `succeeded` is never compared against the rows that reach the file, so a
pass with a second unaccounted fallible step can stamp
`{"attempted":10,"succeeded":10,"dropped":0}` onto an artifact holding five rows. Route the
*last* fallible step before the write, not the first.

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
- `parallel` - **no-op alias, retained for compatibility.** Rayon-based parallelism
  (`ferro_hgvs::parallel`, and `BatchProcessor`'s `*_parallel`/`*_streaming` methods)
  is unconditional — rayon is a non-optional dependency, so there is nothing to enable.
  It formerly gated that API's visibility, which read as a capability gate and was
  misread as one (#1504, #1507).
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

`ReferenceProvider::genomic_placement` exposes these; `VariantProjector::project_to_genomic` composes the `NM_`→`NC_` (cdot) step with the `NC_`→`NG_`/`LRG_` transform to re-express coordinates in the parent frame. With no placement it declines rather than emit chromosome coordinates under the parent accession. That transform is **affine only when the placement carries no alignment gaps** — an LRG record lists `<diff>` elements alongside its `<mapping_span>`, and since #1833 the indel ones are carried on the placement and honoured, so 46 of the prepared reference's 1294 LRG records map through a gap list rather than a pure offset (#1499).

## Reading the spec

Four findings that change how a clause argument should be framed. Each is checkable from the
pinned `assets/hgvs-nomenclature` checkout; re-check rather than trust this list.

### Almost nothing in the recommendations is normative

`style.md:9` binds the spec to RFC 2119, which invites the habit of arguing that clause A is a
SHOULD and clause B a MUST. That argument almost never works here, because **uppercase RFC 2119
keywords appear exactly twice outside `style.md` itself**, and both are in one file:

```bash
grep -rnoE '\b(MUST|SHOULD|RECOMMENDED|MAY|SHALL|REQUIRED|REQUIRES|OPTIONAL)( NOT)?\b' \
  assets/hgvs-nomenclature/docs/recommendations/ | grep -v '/style\.md'
# -> docs/recommendations/RNA/adjoined_transcript.md:20:REQUIRES
# -> docs/recommendations/RNA/adjoined_transcript.md:21:SHOULD    (and nothing else)
```

**The count was recorded here as "exactly once" and that was wrong** — corrected 2026-08-10.
`REQUIRES` is missing from the RFC 2119 alternation everyone copies (RFC 2119 lists `REQUIRED`),
so `:20` — "This syntax REQUIRES the use of a range (not a single position)" — was invisible to
the grep the count was taken from. Neither hit generalises: `:20` constrains one syntax, and
`:21`'s SHOULD is assay-scoped, conditioned on the junction rather than the whole transcript
being analyzed. Do not restate "exactly once" anywhere.

Every clause this project has litigated — `general.md:34`, `:35`, `:56`, `:58`,
`DNA/delins.md:17`, `:18`, `:47`, `DNA/inversion.md:20` — is lowercase prose. Read strictly,
none of them is normative. That does not make them ignorable: it makes most of these questions
house-style choices the spec leaves open, which ferro still has to answer. Argue them on the
spec's worked examples and on downstream cost, not on keyword strength. (`general.md:58` and
`DNA/duplication.md:18` are read as strong because of their wording — "are not allowed",
"**must**" — not because RFC 2119 makes them so; say which you mean.)

The `rulings[delins-merge-vs-individual-gap-two-or-more]` record used to argue from "the two
clauses are the same RFC 2119 strength". It does not any more, for this reason.

### There is no minimality principle, and stability is the spec's first stated value

`background/basics.md:38`: "The recommendations for the description of sequence variants are
designed to be **stable**, **meaningful**, **memorable**, and **unequivocal**." Minimality is
absent, and stability leads. Ferro's column-minimal objective is therefore **our policy**, not
compliance — never cite the spec for it. `DNA/delins.md:44-47` in fact recommends a
*non-minimal* description in its own worked example (66 columns where 40 suffice), and `:47`
recommends the **spanning delins**; the split is what `:46` calls the "alternative description".
That direction has been retold backwards at least once, so state it explicitly whenever you
cite the passage.

### A split is rarely unique, which is a stability argument by itself

Exact enumeration over 40 rows found **27 that admit more than one equally-compliant split**,
median 2 and **max 125**; the spec's own `:44-47` example admits five. So adopting a split
trades one stable canonical form for an arbitrary pick out of a family — and the arbitrariness
is invisible in any single-row before/after comparison, which is how it keeps getting missed.

### Read the Q&A, not only the Notes and Examples

`DNA/delins.md`'s trailing `!!! note` blocks are adjudicative, not colour. A claim that "no
sequence-local rule separates these two cases" was withdrawn once `:86-89` was read: the Q&A
carries the worked answer (`NM_007294.3:c.2077delinsATA`) and `:89` records that the passage
permitting the two-member spelling was *removed* by the committee. `:79-84` likewise carries the
spec's own discriminator — "the two variants may have been reported (or might occur)
individually" — which is **provenance**, recoverable only from the input's spelling.

### A forward-looking note is a suggestion, and may describe a proposal that FAILED

`general.md:36-39` reads as forward guidance — "the SVD-WG is preparing a proposal… The new
recommendation will be: two variants separated by less than two nucleotides should be described
as a `delins`" — and it is tempting to treat it as the direction of travel and discount the
current rule accordingly.

**That proposal is SVD-WG010, and it was rejected.** Word for word, including its rationale.
So the note is stale text describing a change that never happened, and three conclusions drawn
from reading it as guidance are all withdrawn: the codon exception is *not* going away, the
proposal to replace it having failed; it *strengthens* `codon-carve-out-shape-restriction`
rather than undercutting it; and the "spec-admitted instability" argument built on it does not
stand.

So a rejected proposal earns a **negative guard**, not an expectation. The spec corpus builds
210 rows whose only purpose is to catch a frameless separation floor of two — what implementing
SVD-WG010 looks like from the outside — and asserts `guard_violations == 0` over them, with the
denominator asserted non-zero so `0 of 0` cannot pass as a result. See
`spec_conformance_axis.rs`.

**So: cross-check every forward-looking statement in `recommendations/` against the
`consultation/` dispositions before citing it.** "is preparing", "will be", "the new
recommendation" are all flags. The disposition table is inventoried in the consultation slice:
9 accepted, 3 rejected, 8 open, 3 unclear — a rejected proposal is generated as a NEGATIVE
guard, never as an expectation.

### Comparing `c.` positions across numbering zones is OUR policy, not compliance

`c.` has three numbering zones — `-n` upstream of the ATG, plain `n` in the CDS, `*n` downstream
of the stop — and a cis allele's members can sit in different ones. Ferro therefore needs an
endpoint ordering that spans zones, to sort members, detect overlap and decide separation.

**The spec does not give one.** `background/numbering.md` defines each zone but states no rule
for comparing positions across them, and it never speaks about alleles at all: it contains zero
occurrences of "allele" and zero of "member". `refseq.md` contains "allele" twice (`:264-265`),
but both are the population-genetics sense — the wild-type "major allele present in the human
population" — not cis-allele membership. (An earlier version of this note claimed *neither* file
contained the word; that overstated it, and the two hits are worth knowing about so the check is
not re-run and read as a contradiction.)

So ferro's cross-zone ordering is a house rule. **Never cite the spec for it.** Argue it from the
underlying transcript coordinate, which is unambiguous, and say plainly that the `c.` spelling is
a presentation of that.

This is not theoretical. The sequence-changing and denotes-no-sequence classes the spec corpus
found both localise to the `c.72`/`c.*1` transition and to nothing else: the same flush
deletion-plus-insertion shape collapses correctly at all 22 other homopolymer runs in the test
transcript and mis-normalizes only where the pair straddles the CDS/3'UTR zone boundary. See
`the_cds_end_flush_pair_is_its_two_members_normalized_separately` and
`the_five_prime_boundary_masks_the_same_per_member_defect` in `tests/it/spec_corpus_regressions.rs`
— the second is the reminder that the 5' boundary is not a working case, only a masked one.

## Adjudication records: where the open questions live

The `rulings` section of `tests/fixtures/grammar/hgvs_spec_normalization_overrides.json` is the
decision log, pinned by `ruling_records_are_intact` (ids **and** statuses; keep `RULING_STATUSES`
in sync).

**Read the ledger before adjudicating anything from spec text.** The records answer most of what
gets re-argued, and re-deriving a settled decision is this project's most repeated failure. Two
traps make it worse than simply forgetting to look:

- **A record's id states its QUESTION, not its ruling, and the two can be opposites.**
  `separation-is-a-property-of-the-spelling-not-of-the-variant` rules that separation is read off
  the partition **re-derived from the resulting sequence** — the negation of its own name. Reading
  the id and stopping gets the answer exactly backwards.
- **A record's own instructions can be stale.** `delins-merge-vs-individual-gap-two-or-more` says
  "Re-run the census before citing it" and then publishes a grep whose alternation omits
  `REQUIRES`, so following it reproduces the wrong count — one uppercase RFC 2119 keyword outside
  `style.md` when there are two (`RNA/adjoined_transcript.md:20` and `:21`, at spec checkout
  `6f85311`). Re-run the command *and* sanity-check the command.

**Do not trust a count of decided/undecided records written down here, and do not trust which
side of the line a record is on either.** Both failures have happened: the sentence this replaces
said "five of its eight records are `undecided`" long after the ledger had grown and four of those
five had been ruled on, and both tables below have carried a record on the wrong side for days
after it was decided. A stale denominator reads as a claim about coverage; a stale status reads as
an open question that is in fact settled, which is the more expensive of the two because someone
then re-derives it.

The statuses are pinned in exactly one place, `RULING_STATUSES` in
`tests/it/hgvs_spec_normalization_tests.rs`, and `ruling_records_are_intact` fails when the ledger
and that list disagree. Read them from there, or straight off the ledger:

```bash
python3 -c "import json;d=json.load(open('tests/fixtures/grammar/hgvs_spec_normalization_overrides.json'));\
[print(f\"{r['status']:<10} {r['id']}\") for r in d['rulings']]"
```

**Both tables below are now enforced against the ledger**, by
`claude_md_adjudication_tables::the_two_tables_partition_the_ledger_by_status`. Moving a record
between them is part of deciding it, and the build fails otherwise. That guard exists because this
section had drifted three records deep and the existing currency scan could not see it: that scan is
line-based over `.rs` files under `src/` and `tests/`, so a Markdown table at the repository root was
outside it twice over, and the word "open" sits in this paragraph while the ids sit in table rows —
so no line carries both. Widening the scan would not have caught it; the check had to be
section-aware.

The table below lists the records that are **open**, with what each leaves unsettled. Read them
before re-deriving the same argument:

| record | what is open |
|---|---|
| `rna-repeat-range-plus-unit-redundancy` | `RNA/repeated.md:22` calls range-plus-unit invalid, `:27` publishes exactly that shape as valid. Upstream's conflict (#466). Since #1631 every path that answers at all emits `:27`'s form and strict refuses the input — but that is **not** the record being decided for `:27`: the lenient repair emits the anchored `r.-6g[6]` and the *normalizer's* tract maximization widens it back, so the agreement is an artifact of two independent passes |
| `ring-telomere-anchoring` | Must a ring's first `::` segment start at `pter` and its last end at `qter`, so a ring naming only interior coordinates is refused? Ferro currently **accepts** an unanchored ring (pinned by `a_ring_with_no_telomere_anchor_is_still_accepted`) — that is the status quo, not a ruling. Enforcing is not simply "is the first endpoint `pter`": `:28`'s `(pter)_#` / `#_(qter)` forms and `cen` complicate the predicate, and `cen`-anchored segments do not parse at all today |
| `junction-exit-wrapper-scope-in-a-mixed-allele` | An allele mixing a member whose intronic offset ferro MANUFACTURED with one the author spelled intronic. `checklist.md:20` wants a genomic reference on the first; `DNA/alleles.md:16` factors one reference out of the whole allele, and the compact form admits only one accession — so lifting re-spells a member `bare-transcript-intronic-position` says to leave as authored, and expanding abandons `:16`'s form. Ferro **declines** today, which is the status quo #1704 shipped, not a ruling; #1723 made both answers expressible and picked neither |

The `decided` records, and the scope each was decided **at** — read the record before citing it,
because several of these rulings are narrower than their one-line summary:

| record | ruling |
|---|---|
| `adjudication-precedence-order` | **A pointer, not an order.** The canonical ruleset lives in `README.md`; this record cites it, carries the evidence base for it, and holds the escalation register. It deliberately states no ordering of its own — a rule restated in four places is how `CLAUDE.md`, the ledger and a PR body came to disagree, and how a ruling came to cite a record that had never merged |
| `canonical-form-choice-when-both-legal` | **Re-derivation governs.** Where two legal descriptions of one variant compete and no clause selects, ferro derives from the resulting sequence and emits what falls out, subject to every explicit spec tie-break — not the input's spelling, not the previously-shipped string. `general.md:157-160` is the spec's own strongest statement of the method. The counter-evidence is recorded in the record rather than left to be re-discovered as a refutation: `delins.md:83` grounds the split-by-default rule in **provenance**, which no normalizer can recover, and `protein/delins.md:50` keeps two descriptions apart while stating that the resulting proteins are identical |
| `delins-codon-carve-out-gap-one` | `delins.md:18` governs |
| `delins-adjacent-members-when-both-consume-reference` | `substitution.md:32` governs: at separation **zero** the split spelling is marked `class="invalid"` and called "not correct" by name, so two adjacent members are one `delins`. **Scoped** to both members consuming reference bases (`sub`/`del`/`delins`/`inv`). The `ins` half of the carve-out is now **closed** too (2026-08-11): the `c.2077delinsATA` triple — `DNA/delins.md:86-89`, `DNA/substitution.md:93-96`, `RNA/delins.md:68-71` — states the merge in three registers and `delins.md:89` records that the passage permitting the split was **removed**; ferro already emits it, so it moves no row. A `dup` on either side is still open, because `duplication.md:18` would be destroyed by merging, as is an adjacent repeat expansion |
| `self-cancelling-across-ring-junctions` | `DNA/complex.md:130` governs: `general.md:58` does **not** reach the `::`-joined segments of a ring. `:58` is a sub-bullet of the prioritisation bullet `:56`, so it inherits its antecedent — it redirects to a preferred single-type description — and `complex.md:5` defines complex as precisely the case where no such description exists. `:130` settles the operators on identical member content: `::` asserts a join that `;` does not, so a `::` composite is not reducible to its member set |
| `projection-codon-exception-is-decided-by-the-rendered-axis` | **2026-08-11.** `DNA/delins.md:42` governs, and reaches **only an axis that declares a reading frame**: a projection does not re-merge its derived genomic axis to match the coding one. `:42` is a conditional whose second conjunct — "together affecting one amino acid" — is unstatable on a genomic reference, so `:17` governs there unopposed and the members stay individual. `general.md:23`/`:26` make the prefix a claim about the **type of reference sequence**, so `LRG_199:g.…` is genomic however gene-scoped its accession; `general.md:44` is precedent for the spec scoping a rule by naming prefixes. Corrects an over-reading of `delins-codon-carve-out-gap-one`, which settles *when* two variants merge and is silent on which axes. Scoped to **DNA axes**: the `r.` axis merges on `RNA/delins.md:18`'s own authority and is not ruled on here; the `n.` axis is left open |
| `absolute-prohibition-enforcement-stage` | **2026-08-10.** The enforcement stage is **mode-dependent**, uniformly: strict fails at PARSE (strict validates input conformance, not merely parseability); lenient does not validate input conformance and fails only when it cannot NORMALIZE; silent is lenient without messages. A third option the record did not enumerate. It answers the record's own objection to unconditional parse refusal, because rule 1 of the README ruleset is about **output** conformance — accepting a non-conformant input and normalizing it to a conformant output trades nothing. **Being implemented clause by clause**, not all at once: `standards.md:39` landed with #1684 (`W3028`), `numbering.md:52` with #1751 (`W4008`), and `checklist.md:16`/`:45`'s genomic offsets with #1628's parse/normalize half (`W4009`), which also took `prohibition_violating_outputs` to **0**. Still open: `checklist.md:33`'s `ins6` (#1627), the remaining clauses of #1630, and #1629 for the silent arm |
| `alignment-only-symbol-in-a-description` | `standards.md:39` governs: neither `X` nor `-` may appear in a description. The dagger sits inside the symbol cell, so `:39` annotates the row rather than competing with `:36`/`:37`, and the strength grading is moot either way — `general.md:48` admits only IUPAC-IUBMB symbols |
| `rna-axis-alignment-only-symbol-reach` | **REFUSE, 2026-08-12.** `standards.md:47-61` governs — the RNA table publishes exactly fifteen symbols (the DNA table's non-daggered rows, lowercased with `u` for `t`) and neither daggered row has an RNA counterpart, while `general.md:50` names that IUPAC-IUBMB assignment as the source of an `r.` description's symbols. So `x` is not one, and a non-leading lower-case `x` in an `r.` insert is refused. The jurisdiction objection is met rather than waived: `background/standards.md` carries BOTH tables, so its RNA half is an RNA-jurisdiction citation and not `:39`'s dagger stretched across axes. **Symbol-specific, not an alphabet rule** — `Named("alu")` is untouched, and lower-case `x` on the DNA axes is unchanged. A lone or leading `x` never reaches the AST, so only a non-leading one is mode-dependent. Corpus movement is a **structural zero** — `RefShape::all()` has no `r.` shape, so 0 of 58,552 spellings are on the axis |
| `bare-transcript-intronic-position` | `checklist.md:20` governs, **as a conditional clause**: strict input hygiene refuses a bare `NM_…:c.20+2del` (`W4007`), lenient accepts. It did not excuse the 371 bare-`NM_` intronic descriptions ferro emitted itself, and **#1704 closed those**: #670's junction-crossing output is now re-parented onto the genomic reference the crossing already resolved (`NC_…(NM_…):c.20+2del`), coordinate untouched. An input that already names one is still left as authored |
| `confluence-gate-is-apply-equality-on-every-determined-axis` | **The gate is restated, not retired.** "Equivalent inputs produce one canonical output" stands, stated non-circularly as "`normalize` is constant on each equivalence class", where equivalence is apply-equality on every *determined* axis (`EquivalenceLevel::CrossAxisSequenceMatch`) and never `NormalizedMatch` — a relation over `apply`, whose codomain is bases, so it cannot collapse into the normalizer it gates. `general.md:43` governs as the spec's statement that one variant's several descriptions are one object; `duplication.md:148` supplies the counterexample making single-axis apply-equality insufficient. Protein is excluded — translation is many-to-one and `p.` states a consequence, not a denotation. Confluence is asserted only over *decided* pairs; the record also states what it does **not** solve, citing `RNA/repeated.md:20-21` and `delins.md:83`: two descriptions can be apply-equal on every axis and still make different epistemic claims, which is a canonical-form question and must not be encoded as a rung |
| `conflicting-member-geometry-refusal-scope` | `DNA/alleles.md:5` governs — the *definition*, not `general.md:58`, is what reaches nested and coincident-insertion geometries; `:58`'s stated ground literally describes only its own `del`+`dup` example. `general.md:56` is cited to record that it does **not** reach a multi-member allele |
| `delins-merge-vs-individual-gap-two-or-more` | `DNA/delins.md:44-47` governs `:17`, **scoped** to the alignment-coincidence shape `:44-47` describes. Where a separation of two or more arises from anything else, `general.md:34` still governs; unscoped the reading reaches roughly fifteen times the row set the argument was made on. **Scoped by DIRECTION too** (2026-08-11): it reaches the **net-deletion** case only — payload shorter than the span — and does **not** reach net insertions, where the split form stays canonical. That keeps `merge.rs`'s `payload.len() < span.len()` gate, whose grounds are the rejection of SVD-WG010 (`:5`), whose `:16` example is itself a net-insertion merge; `DNA/duplication.md:90-92`, which publishes a net insertion as a `[dup;ins]` split; and the fact that the ruling's whole evidence base is net-deletion |
| `delins-payload-coincidence-carve-out-is-coding-dna-scoped` | **Scoped to the coding DNA axis.** `delins.md:47`'s payload-coincidence carve-out reaches `c.` and nothing else; on `g.`/`m.`/`o.`/`n.`/`r.`, `general.md:34` governs and members stay individual. `:47`'s stated reason is preventing "incorrect predictions for the consequences on protein level", which has nothing to bite on elsewhere, and `general.md:43` shows the spec says "ALL descriptions" when it means universal. **`r.` is out in both directions** — a DNA document has no jurisdiction over the RNA axis, and `RNA/delins.md` states no `:47` counterpart, so `RNA/delins.md:17` stands unqualified. Do not gate this on `AxisFrame::reading_frame`: that is true for a coding `r.` too, which is how the defect arose. Counterweight in the record: `:47`'s other ground, simplicity, is axis-neutral. **Honoured on the SHIPPING path as of #1616** — it was honoured only by `payload_coalesce_applies`, which gates a non-shipping partition rule, while `partition_block` applied the carve-out on every axis. What it costs is `g.` rows in the payload-coincidence stratum that the ungated pass would have merged — 380 over the designed cis corpus. It does **not** cost the stored 723-base `g.` `delins` its 193-member re-derivation: that claim was withdrawn as measured false (the row is byte-identical on both arms), and its cause is the `coalesce_whole_block_inversion` family |
| `delins-recommendation-reach-when-the-input-arrives-split` | **2026-08-12. `:46` governs — an INSERTED SEQUENCE is what `:47` reaches.** `:46` states the mechanism the carve-out rests on, and states it as "parts of the **inserted sequence** 'align'", so `:47` reaches a payload-coincidence split only where some derived member supplies bases while consuming a *different* number of reference bases. A split of pure deletions inserts nothing, so nothing re-aligned, so `general.md:34`/`:17` govern it unqualified. It is the only reading under which both worked examples in the passage come out as the spec writes them — `:44-47`'s own example merges (its `895_896delinsC` is gap-bearing) and W58 stays split. A further scope on `delins-merge-vs-individual-gap-two-or-more`, sitting under that record's shape and net-deletion scopes and under the `c.`-only axis scope; it widens none of them. **Not a confluence ruling** — the pass runs last on already-derived pieces, so each class converges either way; this picks *which form*. Already implemented by `split_carries_a_gap_bearing_insert` (#1698) on the `CanonicalCoalesced` arm only, so it moves no shipped row and the real-corpus disclosure is owed by the default flip. Counterweight in the record: `:47`'s other ground, simplicity, is mechanism-neutral and is given no weight — the same second ground the axis scope also declines, so both narrow together if it is ever accepted |
| `inversion-vs-two-delins-76-83` | `inversion.md:5` governs: a whole-block reverse complement is one `inv` when the members it competes with are `delins`, since `general.md:56` ranks substitution but not `delins`. #1230's substitution case is untouched and still splits |
| `inversion-vs-a-mixed-member-competitor` | **2026-08-12.** The same question where the competitor is a MIX — lone substitutions beside multi-column members. `inversion.md:5` governs, as a **permission**: the `inv` is conformant, the competitor is conformant too, so this is a `README.md` rule 6 choice and must never be cited as conformance. It deliberately does **not** rest on `general.md:56`, which read plainly argues the other way (every competitor contains a substitution, rank (1), above inversion); the reading that saves it is `conflicting-member-geometry-refusal-scope`'s, that `:56` does not reach a multi-member allele. `general.md:58` was measured and does not reach it either — no competitor contains a `del`, `dup` or `ins`. The ground held is that notation must not turn on base coincidence |
| `contiguous-insertion-split-by-a-blocked-derivation` | **NO** — `general.md:34` is stated over "two variants", and this locus carries one. Two spellings of one contiguous insertion are not two variants separated by unchanged nucleotides, so `:34` does not reach them; the split survives because a derivation was blocked, not because a clause requires it. Verified by applying both spellings through `hgvs_to_spdi`, independently of the normalizer |
| `separation-is-a-property-of-the-spelling-not-of-the-variant` | **Read the ruling, not the id** — the title states the *question*, and the answer is the opposite of it. The separation `general.md:34` keys on is read off the partition **re-derived from the resulting sequence**, never off the input's spelling. On its own case both spellings converge on `g.[1001009_1001010del;1001013del]`. Deviates from `delins.md:79-84`, whose provenance rationale ferro cannot honour |
| `separation-rule-force-modal-or-negation` | **RULE 2, whole clause.** `general.md:34`'s "and not" names the excluded alternative; the **modal** grades the clause. Decisive internally: "described individually" and "described as a delins" are complements, so forbidding one would make the other mandatory and leave "should" doing no work — a clause cannot grade its positive half as a preference while forbidding the sole alternative. Confirmed by `DNA/delins.md:81` restating this very rule with "**preferably**", and by the spec pairing "and not" with "must" only where the alternatives are *not* exhaustive (`duplication.md:18`'s "and not as, **e.g.**, an insertion"). The prohibition line the spec does draw is at separation **zero** (`substitution.md:32`, `class="invalid"` + "not correct"). Changes classification only: rule 2 still binds and outranks maintainer judgment, but such an output is a deviation to **disclose and pin with a tripwire**, not a rule-7 bug, so it does not by itself block a release |
| `codon-carve-out-shape-restriction` | **WIDEN.** The `delins.md:18` / `general.md:35` exception applies wherever its stated precondition holds — two variants separated by one nucleotide, together affecting one amino acid — **regardless of edit type**. Ferro's restriction to a sub/unchanged/sub shape is dropped, on the ground that edit type is a property of the spelling while "together affecting one amino acid" is a property of the resulting sequence |
| `derivation-may-not-be-bounded-by-the-inputs-spelling` | **DELETE THE BOUND.** No clause compares a description to the *input*; the recommendations compare it to the sequence, and `background/basics.md:38`'s design values omit minimality. So `canonicalize_from_sequence` may not refuse a derivation for weighing more than the member list the caller happened to write — that comparand contradicts `canonical-form-choice-when-both-legal` in terms. **Read the scope**: the ruling removes an input-relative comparand, it does **not** license the merges the comparand happened to be blocking. Those are governed by `general.md:34` and by `delins-merge-vs-individual-gap-two-or-more`. **Do not attribute the masking to `best_alignment`'s single-gap restriction (#1617)** — the record measures the producers as `partition_block`'s two payload-coincidence exits plus the `MAX_SPLIT_BLOCK` length cap, and withdraws that earlier attribution by name. **All three are now closed, on the SHIPPING path** (#1616): the carve-out is scoped to `c.` inside `partition_block`, and a block the length cap hands back **unexamined** may not be emitted as one spanning member off that axis — a block nothing looked at is the absence of a finding, not a finding of no separation. `MAX_SPLIT_BLOCK` itself is untouched and raising it stays open. Costs 18 (3') / 14 (5') `spec_conformance_axis` classes their confluence, which rule 2 outranks |
| `unchanged-is-read-over-every-minimal-alignment` | **COLUMN-BASED, OVER EVERY MINIMAL ALIGNMENT, 2026-08-07 — recorded 2026-08-14.** A reference base is unchanged iff **every minimal alignment** matches it; `general.md:34`'s "separated by one or more nucleotides" and `delins.md:16`'s "two or more consecutive nucleotides" are both read against that notion. Decided while building the forced-unchanged-column detector (#1539/#1540) and left in a session log for a week, during which **two agents re-derived it differently within one hour** — this row exists because that is the exact cost the record-every-adjudication policy names. The distinction is cell-based vs column-based: `Dominators::matched_ref()` asks whether one `(ref, alt)` CELL matches on every minimal path, the clause asks whether a reference BASE is unchanged in every minimal alignment; on `GACA -> AGAT` both cost-3 alignments match ref offset 1 but to different alt offsets, so the cell notion sees nothing while the base is unchanged either way. **Equal length does NOT imply no indel** — `CAG -> AGA` has edit distance 2 (del + ins), so its two middle bases are unchanged and the changes are described individually; `delins.md:16` has no antecedent because nothing is two CONSECUTIVE CHANGED nucleotides, and `substitution.md:32` does not reach it because its own `GC -> TT` example is one where two substitutions IS minimal. **Falsifies a premise stated three times in `merge.rs`** — that an equal-length block "has no gap to place" so every matched base is a coordinate-wise identity; all three guards key on `reference.len() != result.len()` and so exempt precisely the class that breaks them, which is why `FERRO_PARTITION=live` gives the position-wise answer here and the canonical family does not. **Scope is UNMEASURED**: the affected class is every equal-length block whose minimal alignment is cheaper than position-wise, and a generator must key on **edit distance against block length**, not length alone. Not a minimality principle — minimal alignment decides which BASES ARE UNCHANGED, a fact about the sequences, never which description is preferred |
| `exon-junction-dup-converge-from-the-far-side` | **CONVERGE.** `LRG_199t1:c.3922dup` normalizes to `c.3921dup`. The canonical position is the most 3' position that does **not** cross an exon/exon junction, reached from **either** side — a description approaching the junction stops at it, one already spelled past it is pulled back. Grounded three times in `DNA/duplication.md` (`:26`, `:60`, `:148`), which is what carries it past the RFC 2119 census's lesson about single lowercase-prose clauses |
| `c-and-n-positions-are-flat-transcript-offsets` | **FLAT, 2026-08-12.** A `c.`/`n.` number is an offset on the transcript reference sequence, resolved as `cds_start + N - 1` and never by walking an exon table. `background/numbering.md:52` governs — `n.` numbering runs "from the first to the last nucleotide of the reference sequence", a count over the sequence's own bases with no term for an alignment — with `:21` the second authority, tying the `c.` axis to those same nucleotides by anchoring `c.1` on the sequence's own start codon. `:40` and `:44` are **supporting context, not authority**, and the record says so: `:40` is about introns and flanking regions rather than about bases that fail to align, and `:44` is an out-of-bounds rule `c.4818` does not engage. **The gaps are REAL and the ruling does not say otherwise**: 58 of 474,818 GRCh38 multi-exon cdot builds carry one, 23–2718 bases, which is why #1665 was closed for calling them malformed. What settles it is whose coordinate space they are in — RefSeq gives `NM_033517.1` a contiguous exon table, a GenBank `CDS 1..5196` and a 1731-aa `NP_277052.1`, all 39 bases longer than cdot's exon-covered span, so the hole is a fact about cdot's **genome alignment**. **Scoped to the sequence frame**: genome↔transcript mapping is untouched and stays exon- and CIGAR-aware. **Moves rows on gapped transcripts only.** Explicitly does **not** settle the sibling provider defect — cdot's gap-collapsed `start_codon`/`stop_codon` reaching a flat-space `cds_end`, so the live provider serves `NM_033517.1` `cds_end = 5157` against RefSeq's `5196` |
| `duplication-must-ranks-the-label-not-the-partition` | **THE LABEL, NOT THE PARTITION, 2026-08-13.** `DNA/duplication.md:17` governs — it is the clause that scopes when `dup` "may only be used" at all, and `:18`'s MUST is its sub-bullet, so the MUST ranks a *label* for a change and never requires that the edit set be partitioned so as to expose one. Applied per piece of the partition re-derived from the resulting sequence: a piece whose payload is not a copy of the reference bases immediately 5' of its insertion point is not a duplication, so `:18` never fires on it. **Two independent grounds, both cited**, because the calibration bar from `exon-junction-dup-converge-from-the-far-side` is that one lowercase-prose clause is weak authority — the second is definitional and in another file: `background/glossary.md:310-311` makes `:18`'s subject, "a variant", a difference between a reference and an observed sequence, i.e. an object prior to any spelling. **Read the scope.** It does *not* reach the `dup` half of `delins-adjacent-members-when-both-consume-reference`'s separation-zero carve-out, which is a different mechanism (merging two adjacent members, destroying a `dup` they carry) and is left exactly as that record left it. **Moves no row on the shipped default**; under `FERRO_PARTITION=canonical-coalesced` it authorises six re-pins, two of which are `contiguous-insertion-split-by-a-blocked-derivation`'s own gap closing and two more of which that record's REPRESENTATION EFFECT paragraph predicted by description. Its falsifier — `DNA/insertion.md:5`, which makes an `ins` that copies its 5' flank a rank-1 defect — is enforced by `tests/it/duplication_label_not_partition.rs` rather than asserted. The counterweight is recorded **unanswered**: the reading makes a spec MUST contingent on how coarsely the derivation cut, and the cut is set by `MAX_CANONICAL_WINDOW` and by an unrelated third member's distance |
| `c-description-against-an-unresolvable-cds-is-refused` | **REFUSE, 2026-08-13.** A `c.` description against a transcript whose CDS the reference cannot resolve is refused, not normalized. It used to come back as its own output — `status: ok`, `changed: false`, `warnings: []`, byte-identically in strict, lenient and silent — so on such an accession `normalize` was the IDENTITY and all five of `NM_000546.3:c.528del`…`c.532del`, one deletion from one `CCCCC` run, survived as five distinct "normalized" strings. `background/numbering.md:21` governs: a `c.` number is defined by the position of the `ATG` start codon, so with nothing locating that codon `c.528` names no base ferro can identify; `refseq.md:129` is the second authority on what a `c.` reference IS. **`:52` and `:44` are supporting context, not authority** — `:52` is why the `n.` axis survives untouched, `:44` reaches only the accept-`c.999999del` half. **Refused at NORMALIZE, not parse**, and this is where it parts company with `absolute-prohibition-enforcement-stage`: the input's spelling is conformant and the fact is the REFERENCE's, which `parse_hgvs` holds no provider to see; the same string against a reference carrying the CDS normalizes correctly. **Unconditional across strict/lenient/silent**, which adds no mode policy — lenient's own contract under that record is that it "fails only if it cannot NORMALIZE". Reuses `project`'s own `ConversionError` / `has no CDS start`, closing the contradiction where one binary answered `ok` on the input its own `project --axis c` exited 1 on. **The `n.` axis is untouched and must stay so.** **The spec's own worked example is in the moving set** — `refseq.md:145` publishes `NM_004006.1:c.5697del`, and that version's CDS ferro's committed slice does not resolve; `project` already refused it before this change. **A real, one-directional migration**: ~68,100 of 225,662 served `c.`-addressable accessions move from answered to refused (~66,900 `XM_` predicted models; 1,242 curated `NM_`, of which 607 are superseded versions — a legacy-input exposure, not "ferro mis-normalizes TP53"). Pinned by `tests/it/issue_1870_cds_less_transcript_refusal.rs` |
| `coding-axis-merges-are-a-disclosed-general-34-deviation` | **A DISCLOSED RULE-2 DEVIATION, 2026-08-13.** Ferro deviates from `general.md:34` and its axis-local twin `DNA/delins.md:17`; `DNA/delins.md:81` governs, as the clause whose "**preferably**" carries the force (the roles are split because the schema refuses one clause being both authority and departure — `:34` reaches these rows, and that is the ruling). #1616 merges eight coding-axis multi-member alleles whose members are authored two or more unchanged nucleotides apart, and `:34` says they are described individually. Stated plainly rather than reasoned away. The chain is that no member is gap-bearing — all eight are pure deletions, `retained == payload length` in every row — so `DNA/delins.md:47` does not reach them on `delins-recommendation-reach-when-the-input-arrives-split`'s reading (#1760), `DNA/delins.md:17` and `general.md:34` govern unqualified, and `general.md:35`'s codon exception cannot reach a separation of two by its own antecedent. `derivation-may-not-be-bounded-by-the-inputs-spelling` says in terms that it does **not** license the merges its deleted comparand was blocking, so #1616 is no authority for its own rows. `separation-rule-force-modal-or-negation` grades `:34` rule 2, whose remedy is disclosure plus a tripwire — and **the tripwire is #1710's `coding_axis_separation_two_or_more_merges` re-pinned 0 -> 8** in both directions; a rise above 8 must name the clause that carried it. The record carries the reading that was put to the operator and **declined** — that the separation is read off the derived partition, which is one member, so `:34` has no antecedent — and its falsifier: a row that gains a gap-bearing member leaves the set, so 8 is a set to re-derive, not a constant to defend |

Two things follow for representation-stability work specifically. The repository's doctrine once
read as "do not move a normalized output", while the downstream filer has said instability is
acceptable **provided it is declared a breaking change** — different bars, and the README ruleset
chooses between them: disclosure is the obligation, so the filer's bar is the project's bar. And
Mutalyzer is not a spec oracle:
its normalizer re-derives a description by minimizing a **weighted description length** with
constants dated 2014, and has no separation rule at all, so a separation disagreement with it is
two objectives meeting rather than evidence it knows something the spec does not. The full
forensics are in that record and on `MIN_PIECE_SEPARATION` in `src/normalize/merge.rs`.
