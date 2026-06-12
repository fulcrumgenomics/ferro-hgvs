# CI perf-regression gate — design

Date: 2026-06-10
Author: nilshomer
Status: design (pending implementation)

## Motivation

PR #585 fixed a perf regression that went unnoticed for months. The cdot
transcript cache (`.bin`) silently fell back to re-parsing a 512 MB JSON on
every `ferro` invocation, inflating `ferro normalize` startup from ~1.0 s to
~4.8 s. It was silent because of three compounding failures:

1. The cache layout drifted; the stale `.bin` failed to load.
2. The fallback warning went through `log`, which the CLI never initializes —
   so nothing was printed.
3. `prepare` skipped regenerating a stale `.bin` when one already existed.

No test and no CI job exercised real-data startup, so a ~5× regression shipped
undetected. This gate exists so that this *class* of bug — a fast path silently
degrading to a slow path — fails loudly, both per-PR (deterministically) and
nightly (against real reference data).

## Goals

- Catch a silent fast-path → slow-path fallback in the cdot cache **per-PR,
  deterministically**, with no reference data and no timing thresholds.
- Catch a gross real-data startup regression (the #585 ~5× class) **nightly**,
  against the real prepared manifest, with a fixed wall-time budget.
- Reuse existing CI infrastructure (the `nightly-mutalyzer.yml` prepared-manifest
  cache) rather than provisioning new multi-GB reference data.

## Non-goals (v1)

- No per-PR wall-time benchmarking (too noisy on shared GitHub runners).
- No baseline database, trend storage, or criterion-in-CI.
- No annotate-vcf / TranscriptDb (#593) budget or structural test (deferred to v2).
- No parse-throughput / FASTA / convert budgets — not the regression class #585
  represents.
- Determinism regressions stay covered by the existing conformance gate and the
  #583/#588/#594 fixes; out of scope here.

## Prerequisite — sequencing (IMPORTANT)

This gate protects behavior introduced by the **#585 → #591 → #593 cache stack**,
which is **not yet merged to `main`** as of this writing. Specifically:

- The structural test's *self-heal* assertion depends on #585 (current `main`
  `load` does not regenerate a stale cache).
- The `ferro check` cache-path co-gate depends on `ferro check --build-cache`,
  added in #585.
- The fast-path representation moves from **bincode** (current `main`) to **rkyv**
  in #591. This spec writes "archive" to mean *whatever fast binary cache is in
  `main` when the gate lands*, so it survives that swap.

**Therefore: implement this gate on top of the merged stack.** The
`ci/perf-regression-gate` branch must be rebased onto `origin/main` *after*
#585/#591/#593 land, before the implementation (not this design doc) is written.
The design doc itself has no code dependency and can be committed now.

## Architecture — two layers

| Layer | When | Mechanism | Catches | Reference data |
|-------|------|-----------|---------|----------------|
| A. Structural test | per-PR | deterministic Rust unit tests | silent fallback (the #585 root cause) | none |
| B. Startup budget | nightly | wall-time ceiling + cache-path co-gate | gross real-data regression (the #585 symptom) | shared prepared manifest |

Layer A is the primary guard: it makes the *root cause* (a silent fallback)
impossible to merge. Layer B is the backstop: it catches real-data regressions
A cannot see (e.g. an algorithmic slowdown at scale, or an eager init added
elsewhere), and re-checks the cache path against real data.

---

## Layer A — per-PR deterministic structural test

### Production change: make the load path observable

The root API hole is that `CdotMapper::load` returns `Result<Self>` with **no
signal of which path it took**. A caller cannot tell archive-load from
JSON-fallback, which is exactly why #585 was silent. Close the hole:

In `src/data/cdot.rs`:

```rust
/// Which path `load` actually took. Makes a silent JSON fallback
/// impossible to ignore at the API level (root cause of #585).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CdotLoadSource {
    /// Loaded from the fast binary archive (bincode today, rkyv after #591).
    Archive,
    /// Fell back to parsing the source JSON.
    JsonFallback,
}

impl CdotMapper {
    /// Like `load`, but reports which path produced the mapper.
    pub fn load_with_source<P: AsRef<Path>>(
        path: P,
    ) -> Result<(Self, CdotLoadSource), FerroError> {
        // current `load` body, returning the tagged source on each arm
    }

    /// Unchanged public behavior.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, FerroError> {
        Self::load_with_source(path).map(|(mapper, _)| mapper)
    }
}
```

- The successful archive arm returns `CdotLoadSource::Archive`.
- The fallback arm (currently the silent `log::warn!` at `cdot.rs:1109`) returns
  `CdotLoadSource::JsonFallback` instead of being swallowed.
- `load`'s signature and behavior are unchanged — this is purely additive.

This change is independently useful: it is also what powers the `ferro check`
co-gate in Layer B.

### Tests

Three `#[cfg(test)]` unit tests in `src/data/cdot.rs`, placed alongside the
existing `test_bincode_roundtrip_*` block (same module already does temp-file
archive round-trips, so they share helpers and idioms). All use the existing
inline `multi_build_cdot_json()` fixture written into a `tempfile::tempdir()` —
no committed data files (per project rule: generate test data programmatically).

1. **Fast path is taken (round-trip).** Build a mapper from the inline cdot JSON,
   write the sibling archive (`to_bincode_file` / the archive writer), then
   `load_with_source(json_path)` → assert source `== Archive` **and** key data
   round-trips (transcript count, a `base_to_versioned` entry, `primary_build`).
   Catches "layout drifted but load still claims success."

2. **Stale/corrupt cache surfaces, not silent (negative).** Write a valid
   archive, truncate/corrupt it on disk, then `load_with_source(json_path)` →
   assert it still returns `Ok` (graceful) **but** source `== JsonFallback`.
   Encodes "a broken cache must be observable, not silently absorbed."

3. **Self-heal regenerates (the `prepare`-skipped-stale mode; depends on #585).**
   After the corrupt-fallback load, drive the regeneration path (whatever
   `load`/`prepare` uses to rewrite a fresh archive), then `load_with_source`
   again → assert source `== Archive`. Proves a stale cache heals rather than
   re-parsing JSON forever.

No timing, no thresholds, no reference data — fully deterministic, runs in the
existing per-PR `Test` job (`cargo nextest run --features dev`).

---

## Layer B — nightly real-startup budget

### Workflow file

New file `.github/workflows/nightly-perf.yml` — **not** a job appended to
`nightly-mutalyzer.yml`, because the two have different failure semantics: the
mutalyzer nightly is `continue-on-error: true` (surfaces drift, never gates),
whereas the perf gate **must fail** on a budget breach. Separate files keep
ownership clean and allow independent `workflow_dispatch`.

Conventions (match the rest of the repo):
- `runs-on: ubuntu-latest`, `permissions: contents: read`, `timeout-minutes: 30`,
  a `concurrency` group mirroring the mutalyzer file.
- All actions pinned by commit SHA with a `# vX.Y.Z` comment — reuse the exact
  pinned SHAs already in `nightly-mutalyzer.yml` (checkout, `dtolnay/rust-toolchain`,
  `actions/cache`).
- The timing/gate steps use `shell: bash` with `set -euo pipefail`.
- **No `continue-on-error` on the gate steps.**

### Shared prepared-manifest cache (the load-bearing detail)

Copy the `Cache prepared reference manifest` step **byte-identically** from
`nightly-mutalyzer.yml`:

- `path: benchmark-output/`
- key: `${{ runner.os }}-ferro-manifest-${{ env.MANIFEST_CACHE_VERSION }}-${{ hashFiles('Cargo.lock', 'src/prepare/**', 'src/bin/ferro.rs') }}`
- `env.MANIFEST_CACHE_VERSION: v1`

GitHub Actions cache is repo-scoped, so an identical key restores the
already-prepared `benchmark-output/` that the mutalyzer nightly populated — no
second multi-GB `ferro prepare`. Keep the `if: steps.cache-manifest.outputs.cache-hit
!= 'true'` guarded `ferro prepare --output-dir benchmark-output --genome grch38`
fallback so the perf job self-heals if it ever runs on a cold cache.

**Coordination:** add a cross-referencing comment in **both** nightly files —
"manifest cache key is shared with the other nightly; keep this key and
`MANIFEST_CACHE_VERSION` identical, bump them together."

**Schedule:** mutalyzer runs `0 4 * * *`. Run perf at `0 5 * * *` (05:00 UTC) so
the mutalyzer job populates the manifest cache first and perf gets a warm
restore, avoiding two concurrent cold `prepare`s. (Also clear of the Sunday
00:00 external-validation run.) Plus `workflow_dispatch`.

### Measurement protocol

- **Command measured:** `ferro normalize --reference benchmark-output` (the
  `ferro` binary takes a reference *directory* via `--reference`; it does not
  read the `FERRO_MANIFEST` env var, which is a test-harness-only convention),
  fed a tiny fixed input written to `/tmp` from a heredoc (3–5 `c.`/`g.` variants).
  Startup is ~90% of wall time, so input size is irrelevant — we are deliberately
  measuring **startup**, which is where #585 lived.
- **Protocol:** 1 untimed warmup run (primes OS page cache → measure code, not
  disk), then **N = 7** timed runs. Capture wall time via `date +%s.%N` deltas
  around the binary; write each run's seconds to a temp file, then compute min/median
  (do not pipe a long-running command straight into a parser).
- **Statistic + ceiling:** gate value = **min of 7** (startup is fixed-cost work;
  min strips upward scheduler noise and exposes the floor a silent fallback shifts).
  **Hard-fail if min > 2.5 s.** Today's floor ~1.0 s; the bug was 4.8 s. 2.5 s sits
  ~2.5× above healthy and ~2× below the bug — comfortably outside runner noise.
- **Soft warn (summary only):** also compute median; if median > 1.5 s, emit a
  `::warning::` and a bold row in `$GITHUB_STEP_SUMMARY`, but **do not fail**.

### Deterministic co-gate (the true #585 catch)

A fast runner could squeak under 2.5 s while still silently re-parsing JSON, so
timing alone is not a reliable #585 catch. Add a timing-free assertion:

- Surface `CdotLoadSource` through `ferro check` (and `ferro check --build-cache`):
  print `cache: archive` vs `cache: json-fallback`.
- After the timed runs, the perf job runs `ferro check` against the prepared
  manifest and **hard-fails if the output is not the archive path.**

This catches the real #585 failure (silent fallback) independent of runner speed.
It is the single most important part of Layer B.

### Job summary

Write a table to `$GITHUB_STEP_SUMMARY`: min / median / all 7 timings + the
cache-path line. No artifact upload needed for v1 (the summary suffices).

---

## v1 scope summary

- **Per-PR (deterministic, no reference data):** add `CdotLoadSource{Archive,
  JsonFallback}` + `load_with_source` in `src/data/cdot.rs`; `load` wraps it.
  Three unit tests beside `test_bincode_roundtrip_*`: (1) round-trip uses
  `Archive`, (2) corrupt cache → `JsonFallback` not silent, (3) self-heal
  restores `Archive`. Fixture = existing inline `multi_build_cdot_json()` +
  `tempfile::tempdir()`.
- **Nightly (`.github/workflows/nightly-perf.yml`, new, `0 5 * * *`):** restore
  the shared `benchmark-output/` manifest via the byte-identical cache key +
  `MANIFEST_CACHE_VERSION: v1`; `prepare` only on cache miss.
- **Budget:** measure `ferro normalize` startup, 1 warmup + 7 timed runs;
  hard-fail if min > 2.5 s; median > 1.5 s → warn-only in summary.
- **Deterministic co-gate:** `ferro check` must report `cache: archive`;
  hard-fail on `json-fallback` — the true #585 catch, independent of runner speed.
- **Conventions:** actions pinned by SHA + version comment, `shell: bash` /
  `set -euo pipefail`, `permissions: contents: read`, no `continue-on-error` on
  the gate; keep the shared cache key in sync via cross-referencing comments.

## Deferred to v2 (documented, not built)

- annotate-vcf / TranscriptDb (#593) wall-time budget — add once #593 is in
  `main` and a small GFF/VCF pair is cacheable; ceiling ~1.5 s (today ~0.6 s).
- A TranscriptDb structural test — #593's binary cache is the same silent-fallback
  class and deserves the identical `LoadSource` treatment as a fast-follow.

## Risks / open questions

- **Stack-merge sequencing** (see Prerequisite): if the gate is implemented before
  #585/#591/#593 merge, test #3 and the `ferro check` co-gate will fail because
  the self-heal / `--build-cache` behavior doesn't exist yet. Mitigation: build on
  the rebased post-stack branch.
- **rkyv archive-writer API (#591):** the exact "write the sibling archive" and
  "regenerate" calls in tests #1/#3 must match whatever #591 lands; confirm the
  symbol names when rebasing onto the stack.
- **Runner-floor drift:** if `ubuntu-latest` hardware gets slower, the ~1.0 s
  floor could rise. 2.5 s has wide headroom, but revisit the ceiling if the
  warn-threshold (1.5 s median) starts firing routinely on healthy runs.
