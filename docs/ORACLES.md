# Normalization oracles

This repo uses a set of "oracles" as an extra correctness check. When armed, they validate every
normalization output the test suite produces, in the CI jobs that arm them and in local armed runs.
Each oracle is switched on by setting its `FERRO_ASSERT_*` flag; an unarmed debug build or the plain
`test` job runs none of them, and release builds carry no oracle code at all. Read this page when an
oracle fires on your PR, when you add a normalization path, or when you change the oracle jobs in
`ci.yml`.

Each oracle asks a question of an output: is it a fixed point, does it re-parse, are its
coordinates in range, does it denote the same bases as the input? They catch a wrong output that
passes every test written against an expected spelling. The sections below say how to run the
oracles, what each catches and misses, where CI arms them, and how a fire blocks a merge.

Terms used on this page:

- **oracle**: a check that judges an output by a property, with no expected spelling to compare
  against.
- **seam**: the single exit every normalization passes through. The oracles run there.
- **arm**: set a `FERRO_ASSERT_*` flag, which switches on its oracle.
- **fire**: an oracle fails. It panics, so the test that produced the output fails.

## The seam

All oracles run from a single call site. `Normalizer::assert_seam_oracles` runs at the exit
of `normalize_core_checked`, and every public normalization path uses that exit: `normalize()`,
`normalize_with_diagnostics()`, and every `VariantProjector` axis. Route a new normalization path
through `normalize_core_checked`. Do not call an oracle from anywhere else.

The checks run in a fixed order, putting more expensive checks last. Currently that order is
in-bounds, re-parse, idempotency, and denoted-sequence.

Every oracle carries `#[cfg(debug_assertions)]`, so a release build contains no oracle code. Each
flag is read once into a `OnceLock`, so a debug run with the flag unset pays one atomic load. The
idempotency oracle re-enters normalization to check its own output. A thread-local
`IN_IDEMPOTENCY_CHECK` guard makes the inner call skip its check.

## Running the oracles locally

```bash
scripts/run_oracle_suite.sh                     # run the armed suite
scripts/run_oracle_suite.sh --print-selection   # list what it would run, then stop
scripts/run_oracle_suite.sh -E 'test(my_test)'  # other arguments pass through to nextest
```

The runner selects the `oracle` profile in `.config/nextest.toml`. That profile excludes the
modules that cannot run armed and binds `scripts/arm-oracles.sh`, which sets the `FERRO_ASSERT_*`
flags. CI's `test-oracle` job selects the same profile and adds only scheduling: it shards the
run and leaves the sweeps, censuses and proptests to other jobs. A local run therefore includes
tests that job does not.
[Which oracles the CI arms, and where](#which-oracles-the-ci-arms-and-where) lists every job
and the flags it sets.

To change what the armed run excludes, edit the profile. Neither `ci.yml` nor the runner needs a
matching edit. `tests/it/oracle_exclude_invariant.rs` reads the profile and fails when it
disagrees with the test tree.

Do not set the flags by hand over the whole suite. Some modules pin a wrong output on purpose or
count conformance rows. An armed oracle fails on them by design, so that run is red on `main`.
The red is not a coverage gap. Those modules run unarmed in the plain `test` job.

## The oracles

### Idempotency oracle

`FERRO_ASSERT_IDEMPOTENT=1` asserts `norm(norm(x)) == norm(x)` for every normalized output:
normalizing again must change nothing. The check re-normalizes the output to verify this, so it
cannot judge an output that fails to parse.

### Re-parse oracle

`FERRO_ASSERT_REPARSE=1` asserts that `parse_hgvs` accepts the normalized description. The oracle
fires only when the input parsed and the output does not. `parse_hgvs` holds no provider, so this
oracle accepts a well-formed spelling that denotes the wrong bases.

A fire on any shape not listed below is a defect in the producer. Fix the producer. Do not add an
exemption. The exemptions are a closed list:

- `0` and `?`. They are legal whole-allele outputs. `parse_hgvs` rejects them standalone because it
  requires an accession.
- An empty allele, `[]`. Only direct construction reaches it. The projector's own tests build one to
  pin that the projector declines it.
- A non-flanking genomic insertion, the projection pivot. Its coordinates are sound, but HGVS admits
  no spelling for it, so the projector withholds the reported genomic axis
  (`non_flanking_genomic_insertion_anchor`).
- A non-coding downstream position, `n.*N`. `parse_hgvs` rejects it in every mode, but
  `TxPos::downstream` is public API, so a Rust caller can build one. `noncoding_zone_marker` keys
  the exemption on the AST.

### In-bounds oracle

`FERRO_ASSERT_IN_BOUNDS=1` asserts that no coordinate in a normalized description is past the end of
its sequence. The rules are on the doc comment of `merge::first_out_of_bounds_coordinate`. This page
does not repeat them. The oracle does not cover protein axes or an inserted-range payload such as
`g.10_11ins[20_30]`. Idempotency does not stand in for this check: an out-of-range coordinate that
is a fixed point passes it.

### Denoted-sequence oracle

`FERRO_ASSERT_SEQUENCE=1` checks that normalization did not change what a description does. A
description is an edit to a reference sequence. The oracle applies the input and the output to one
stretch of reference, wide enough to cover both, and asserts that the results are the same bases.
One window for both is what lets a 3'-shift inside a repeat yield the same bases instead of a
difference. In other words,
this oracle detects a wrong edit that looks right: one that parses, whose positions exist, and that
normalizing again leaves unchanged.

Read `normalize::denoted_sequence_oracle_counts()` before you trust a green run. It returns
`(compared, skipped)` for the process, and zero comparisons and zero faults look the same. A side
that cannot be applied is counted as a skip, never as a pass:

| case | verdict |
|---|---|
| both apply, bases agree | pass, counted in `compared` |
| both apply, bases differ | fire |
| the output denotes no sequence and the input does | fire. The output lost the input's meaning, for example when two members claim one base |
| the input denotes no sequence (a trans allele, a `REFSEQ_MISMATCH`, an edit SPDI cannot carry) | skip, counted in `skipped`. The input gives no baseline |
| the two name different accessions | skip, counted |
| the union window exceeds `MAX_APPLY_WINDOW`, or the provider cannot serve it | skip, counted |

The oracle stays silent on the shapes below by design. A fire on one of them is a regression in
the oracle: file it against the oracle, not the normalizer.

| class | why the oracle does not fire |
|---|---|
| output cannot be transliterated | the input states its own deleted bases and converts with no provider; the output must read a reference the fixture does not hold |
| insertion flush against a deletion | the applier's tie-break defines the order, so it is not an overlap. See #1749 and #1831 |
| overlap-conflicting input | an insertion interior to a deletion; the input denotes nothing, so there is no baseline |
| `pter`/`qter` | they carry no numeric coordinate, so the applier declines the row before it converts either side |
| corrected `REFSEQ_MISMATCH` | the input names a reference base the reference does not hold; normalization corrects it, so the denoted sequence changes on purpose; the seam's denoted-sequence check skips the row by its warning |
| `r.` payload against a DNA reference | the same bases in two alphabets |
| uncertain allele `[(…)]` | the members are uncertain, so the applier skips the row |

The applier does not call the normalizer or use its result, or the check would agree with whatever
normalization produced. `spdi::compare_denoted_sequences` reaches the bases through `hgvs_to_spdi`
and an SPDI splice, a walk that agrees with `apply_to_reference`. Do not use `EquivalenceChecker`
here: it normalizes both sides. The applier and the normalizer read a `c.` position on the same flat
transcript axis; the ruling record `c-and-n-positions-are-flat-transcript-offsets` says why, and
`CDOT_GAP_JUNCTIONS` guards it.

To add a regression test for the oracle, pin the recorded wrong output and assert that the oracle
fires on it. Do not run the normalizer: a test that re-normalizes goes green when the defect is
fixed and stops testing the oracle. `tests/it/issue_1615_denoted_sequence_oracle.rs` holds these
rows, plus a negative control that a legitimate re-spelling stays silent.

## The oracles in CI

### Which oracles the CI arms, and where

CI splits the test suite across several jobs, and each job chooses which oracles to switch on for
its tests. The `test-oracle` job takes its flags from the `oracle` profile. `sweeps`, `censuses`
and the nightly set theirs in their own workflow steps.

| job or step | IDEMPOTENT | REPARSE | IN_BOUNDS | SEQUENCE |
|---|---|---|---|---|
| `test-oracle` armed step | yes | yes | yes | yes |
| `test-oracle` re-run step | yes | yes | yes | no |
| `sweeps` | yes | yes | yes | yes |
| `censuses` armed step | yes | yes | yes | no |
| `censuses-plain` | no | no | no | no |
| `test` | no | no | no | no |
| `soak` | no | no | no | no |
| nightly | yes | yes | yes | yes |

To turn an oracle on in a CI job, first run every test that job runs with the flag set. Passing
the few tests you have in mind proves nothing about the others. `censuses` leaves
`FERRO_ASSERT_SEQUENCE` off by design until that run is done. One caveat when reading results:
`test-oracle` provisions no `FERRO_MANIFEST`, so the tests that need one return early and pass,
and a green `test-oracle` says nothing about them.

### What the oracle profile excludes

The `oracle` profile excludes two lists. The first is modules that cannot run armed at all. The
second is single rows the denoted-sequence oracle fires on for a bug not yet fixed.

The module list. These modules still run unarmed in the plain `test` job, so the list is not a
coverage exemption. A census counts spec rows. Armed, `conformance::census::measure` catches the
oracle's panic and files the row as `declined`, so the count reads better than the truth.
`conformance::census::run_census` refuses to run instead, returning `CensusError::OracleArmed` when
any `FERRO_ASSERT_*` flag is set. A pinned-defect module records a known wrong output, so an oracle
fires on it by design. Every module that reads the spec corpus
must be on this list, whether or not it fires, and every name on the list must read the corpus.
`tests/it/oracle_exclude_invariant.rs` fails on either miss. When you add a test module that
imports the corpus, add it here.

The debt list. These rows are excluded from the denoted-sequence oracle only. The `test-oracle`
re-run step runs them under the `oracle-rerun` profile with the other oracles armed. Each row sits
beside the open issue that retires it. Read the rows off `.config/nextest.toml`. To exclude a new
fire, add the row to both profiles beside its open issue. Never add an exemption inside
`assert_seam_oracles`: a row in the profile is visible and names its issue, and an exemption in
code hides the fire from every run. When the issue closes, remove the row from both profiles.

### The merge gate and the nightly test

For PRs, CI runs `test-oracle` and `sweeps` with the oracles armed. A fire fails the job, and the
required `Test` rollup needs both jobs, so the fire blocks the merge. Read the rollup's `needs:`
list from `.github/workflows/ci.yml`, not from this page.

The nightly reference-aware run tolerates known failures and alarms on new ones. Its test step
carries `continue-on-error`, so a failing test does not fail the job. A later step diffs the run's
failing set against the committed baseline and fails the job on any difference, which opens a
`report-failure` issue. A new oracle fire is a new failure, so it opens the issue. A fire on a test
already in the baseline does not, so a green nightly does not mean no oracle fired. To see what
fired, read the uploaded xfail artifact; the job summary prints the armed reproduction recipe.
