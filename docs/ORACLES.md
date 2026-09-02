# The normalization seam oracles

The four `FERRO_ASSERT_*` flags that turn every normalization in the test suite into an
invariant check, where each runs in CI, and how to run them locally. Moved out of
`CLAUDE.md` on 2026-09-02; the content is as it stood there, and a follow-on pass will
edit it.

## Running the oracles locally — a bare armed run CANNOT pass

Read this before running any of the four flags below over the whole suite. The
obvious command is documented nowhere else in this file any more, because it is
a trap:

```bash
FERRO_ASSERT_IDEMPOTENT=1 cargo nextest run --features dev   # ALWAYS RED on main
```

It has been red for as long as the spec corpus has existed. The count was **7**
before #1650 closed the idempotency half; it is deliberately not restated here as
a number, because nothing checks it and a stale one reads as a measurement — read
it off a run. Use the local runner instead, which mirrors `ci.yml`'s
`test-oracle` job — its flags, and its selection:

```bash
scripts/run_oracle_suite.sh                     # as test-oracle runs it
scripts/run_oracle_suite.sh --print-selection   # what it would run, without running it
scripts/run_oracle_suite.sh -E 'test(my_test)'  # extra args reach nextest
```

**Why the bare form is red, and why that is not a coverage gap to close.** The
failures come from `FERRO_ASSERT_IDEMPOTENT` alone — measured: `FERRO_ASSERT_REPARSE`
and `FERRO_ASSERT_IN_BOUNDS` each contribute **0**, and a prepared reference
changes nothing (the same set fails with and without `FERRO_MANIFEST`; every one of
these modules is `MockProvider`-backed). All of them live in modules `ORACLE_EXCLUDE`
already names, for two reasons that are worth keeping apart:

- **Some of them USED TO ASSERT the defect the oracle PANICS on.** The
  `defect_non_idempotent_outputs` tests and
  `spec_corpus_regressions::an_insertion_at_the_cds_end_is_a_fixed_point`
  pinned `c.*1delinsCTT` → `c.72_*1insCT` → `c.72delinsCCT` as a known non-fixed
  point, and a test that pins a defect and an oracle that aborts on it cannot both
  run. **#1650 closed that class** — the chain collapses to
  `c.*1delinsCTT` → `c.72delinsCCT` in one pass, `non_idempotent_outputs` reads 0
  in both directions, and those tests now assert the fixed point (which is why the
  last one is named `..._is_a_fixed_point`). So this reason no longer holds for
  `FERRO_ASSERT_IDEMPOTENT`. The exclusion is kept anyway, deliberately: see
  `scripts/run_oracle_suite.sh`'s header for why, and note that narrowing it needs
  its own measurement with `tests/it/oracle_exclude_invariant.rs` updated in the
  same change.
- **Others COUNT it**, and arming the oracle makes the count read *better*
  than the truth. `spec_conformance_axis`'s censuses wrap normalization in
  `catch_unwind`, so a panicking row is filed `declined` and never reaches its
  family's output set. See `ORACLE_EXCLUDE`'s comment in `ci.yml` for the
  measured figures.

So the answer is not "run it in CI anyway" — that is the one thing that would
destroy the evidence. The answer is that CI already excludes it, and the local
path now excludes it the same way, from the same source of truth.

**The runner reads the whole `-E` selection and the flag set out of `ci.yml`**,
so neither can drift into a second copy. "Whole" is load-bearing: `test-oracle`
negates `test(proptest)`, `SWEEP_FILTER`, `CENSUS_FILTER` and — since #1815 —
`SEQUENCE_ORACLE_EXCLUDE` as well as `ORACLE_EXCLUDE`, and a runner that negated
only `ORACLE_EXCLUDE` would run the proptest modules, the three exhaustive sweeps
and the three slow censuses while reporting itself as mirroring that job. So the
shape of the expression is read too, and only its variable references are
expanded — rebuilding it here from a hardcoded `not test(proptest) and not (…)`
would trade today's drift for tomorrow's. **Do not count the negated terms in
prose**; that number has changed twice and is read off the file.

`tests/it/oracle_exclude_invariant.rs` re-derives all of it in Rust — different
anchors, deliberately — and compares; separately asserts that the agreed
selection actually negates each of the named filters (a spelling-independent
check, since both sides otherwise read the same line of `ci.yml` and a matching
pair of wrong ones would satisfy the equality); separately forbids the script from
hardcoding a module name; and separately asserts that the rows withheld from the
denoted-sequence oracle still run under the other three, which is what makes
#1815 a superset rather than a trade. Each of those guards is checked against a
deliberate sabotage rather than assumed.

One anchor is worth knowing about because it is the non-obvious consequence of
#1815: `test-oracle` now has **two** nextest steps, so "the flags of
`test-oracle`" is no longer well-formed. The Rust side identifies the armed step
by its `--partition`, and the runner's awk by the step's `name:` — deliberately
different discriminators, so a rename breaks one and not the other.

## Normalization idempotency oracle

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

## Normalization re-parse oracle

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

## Normalization in-bounds oracle

`FERRO_ASSERT_IN_BOUNDS=1` is the third at the same seam: no coordinate a
normalized description names may be past the end of its own sequence.

`scripts/run_oracle_suite.sh` arms it too. As with the re-parse oracle,
`FERRO_ASSERT_IN_BOUNDS` on its own fails **0** tests over `ORACLE_EXCLUDE`'s
modules — the failures a bare armed run shows are the idempotency oracle's,
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
`test-oracle` job, the `sweeps` job, and the **armed step** of the `censuses`
job (those three and nowhere else — the plain `test` and `soak` jobs run without
the flags, and the `censuses-plain` job runs `spec_conformance_axis`
un-armed); the nightly sets it too, where it is the only place the check runs
against true transcript and contig lengths.

**The fourth flag is not set in all of those places, so do not read "all four"
off this paragraph.** `FERRO_ASSERT_SEQUENCE` runs in `sweeps`, in the nightly,
and — as of #1815 — in `test-oracle`, where it is armed over that job's selection
minus a two-row debt list (`SEQUENCE_ORACLE_EXCLUDE`). It is still deliberately
**not** set in `censuses`, whose armed step is now a strict subset of
`test-oracle`'s flag set rather than a copy of it. The reason, and what the debt
list holds, are under
[Where it runs, and the one place it deliberately does not](#normalization-denoted-sequence-oracle)
below and in `ci.yml`'s comment on the `test-oracle` job.

## Normalization denoted-sequence oracle

`FERRO_ASSERT_SEQUENCE=1` is the fourth at the same seam, and the only one that
asks what the output **means** rather than how it is written: apply the input to
the reference, apply the output, and assert the bases agree.

```bash
scripts/run_oracle_suite.sh    # arms all four, over test-oracle's own selection
```

**`scripts/run_oracle_suite.sh` DOES arm this one as of #1815**, because
`test-oracle` now does — the runner reads that job's flag set rather than
carrying its own, so it started arming the fourth oracle on the day the job did,
with no change to its flag handling. What it did need teaching is the third
negated filter, `SEQUENCE_ORACLE_EXCLUDE`; a local armed run that omits that
negation is red by construction on the rows named there. (This paragraph used to
say the runner deliberately does *not* arm it, "mirroring the gap" — that is the
statement #1815 changed.) Added to
`ORACLE_EXCLUDE`'s modules it raises **5** further failures beyond the
idempotency oracle's own, all in `spec_corpus_regressions` and all the same shape
as those: a test pinning the CDS-end defect that this oracle aborts on. (That 5
was measured before #1650; it is the denoted-sequence class, which #1650 did not
close — `sequence_changed` is unchanged at 4/0 — but re-measure rather than quote
it, for the reason given above.)

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
(gating, measured green over the full corpus), the nightly sets it (non-gating),
and **`test-oracle` sets it as of #1815** — over that job's whole selection minus
`SEQUENCE_ORACLE_EXCLUDE`, a debt list of two rows each carrying the open issue
that retires it.

`censuses`' armed step still does **not** set it. That used to be a *derived*
reason — the step exists to run what `test-oracle` ran, on a faster archive, so it
inherited that job's flag set exactly, and arming a fourth oracle there would make
the move a change of instrument. Since #1815 the derivation no longer holds: its
flag set is now a strict **subset** of `test-oracle`'s, deliberately, and arming
the fourth there has never been measured over `ORACLE_ONLY_FILTER`'s modules. So
"restore parity with `test-oracle`" is now the plausible-looking wrong edit, and
that job's header says so.

**The debt list is what to read before believing any count in this section.** It
has changed composition twice in three days, in both directions: it held
`issue_1487_canonical_window_overflow` until #1990 closed #1690 (2026-08-15), and
gained `dump_normalized_corpus`'s render-time gate when #2051 landed (2026-08-16).
Read it off a run.

Two rows in `test-oracle`'s selection used to fire before any of that, and both
were real disagreements inside ferro rather than noise:

| row | what disagreed | state |
|---|---|---|
| #1618 — `NC_TEST.1:g.262TG[6]` → `g.259_262GT[6]` | `hgvs_to_spdi` read the anchored spelling as 6 copies replacing a **1**-copy tract, the normalizer's output as 6 copies replacing a **2**-copy tract — 14 bases against 12 | closed before `6116f84a` |
| #1619 — `NM_033517.1:c.4818dupC` → `c.4818dup` | `hgvs_to_spdi` resolved the `c.` position by **walking** the exon list while the normalizer indexes the **flat** transcript, so the two disagreed across any transcript-coordinate gap: the input applied `C`, the output `T`, at transcript position 4877. `NM_033517.1` carries a real 39-base cdot hole between exons 10 and 11 (the row it replaced is recorded in #1619) | closed by the flat-frame fix: `cds_to_tx`/`tx_to_cds` no longer read the exon table |

**Both are now green, and that was still not enough — the selection-wide run
those two closures were blocking on came back RED, three times, on a different
set of rows each time.** The history below is why the flag is armed behind a debt
list rather than armed outright, and why "the rows I know about are green" has
never once been sufficient here. The two-row measurement, on the #1619 branch, was
this:

```bash
FERRO_ASSERT_SEQUENCE=1 cargo nextest run --features dev --test it \
  -E 'test(test_explicit_base_removal)'          # 2 passed
FERRO_ASSERT_SEQUENCE=1 cargo nextest run --features dev --test it \
  -E 'test(repeat_tract_maximization)'           # 3 passed
```

That is the two named rows, not `test-oracle`'s selection, which is the whole
suite minus `ORACLE_EXCLUDE`, `test(proptest)` and `SWEEP_FILTER`. The failures are in two
classes, **neither of them #1618 or #1619**, so neither was reachable by the
two-row check above:

- **3x `issue_1487_canonical_window_overflow`** — not an oracle fire at all.
  `attempt to add with overflow` at `src/convert/mapper.rs:114`: arming the flag
  routes the output through `hgvs_to_spdi`, which reaches `cds_to_tx`'s unchecked
  `cds_start as i64 + pos.base - 1` on an `i64::MAX`-adjacent `c.` coordinate.
  The issue was **#1690**, **since CLOSED** — the two arms that can overflow now
  `checked_add` and return a `ConversionError`.
- **2x `stranded_identity_member`** — a genuine fire, of #1281's "denotes no
  sequence" shape, on a module that exists to PIN a defect (#1655's stranded
  derived identity). That is the class `ORACLE_EXCLUDE` already documents: a test
  pinning a defect and an oracle aborting on it cannot both run.

**Re-measured after #1690 closed (#1990), and only one of the two classes
survives.** The three `issue_1487_canonical_window_overflow` rows are **gone** — all 8 tests
in that module pass — and **nothing new fired**. The 2 survivors are both
`stranded_identity_member`, still the recorded row:

| test | the row it fires on |
|---|---|
| `every_stranded_identity_is_non_confluent_with_its_surviving_member` | input `TEMPLATE:g.[10_11insAT;11_12inv;12_13insA]`, output `g.[11_12insTA;11_12=;12_13insA]` |
| `the_wider_multi_member_census_is_what_the_header_says` | the same row, reached through the wider census |

**The `FERRO_MANIFEST` half of that measurement is stale too.** All six
`mutalyzer_normalize_tests` axes now **pass**, and the remaining failure is
`multi_member_cis_axis::axis_normalized`, which is **not an oracle fire**: it
fails identically with `FERRO_ASSERT_SEQUENCE` **unset** (measured as a control),
on a census pin reading `unwindowed: 89` against a pinned `90` while
`sequence_changed` stays `0`. That is the HEAD-drift this file already documents
for that pin, and it is orthogonal to arming the flag.

Three failures, not two — and the third could not have been named by
either reading above, because the test did not exist when they were taken. #2051
(2026-08-16) added
`dump_normalized_corpus::the_render_time_reference_matches_what_the_pipeline_was_given`,
which re-normalizes each corpus row **itself, outside** the `catch_unwind` that
`dump`'s own pass uses, so 47 corpus inputs abort it — **#2140** (38), **#1983**
(8) and **#2139** (1). It passes unarmed in the same tree, so it is an oracle fire
and not a pre-existing red.

So arming it is **no longer blocked** — it is armed, behind those two rows.
`SEQUENCE_ORACLE_EXCLUDE` names them, each with its retiring issue, and
`test-oracle` gains a second un-partitioned step that re-runs exactly those rows
under the other three oracles, so the file's oracle coverage is a strict superset
of what it was rather than a trade. Suppressing a row **at the seam** would still
hollow out the oracle; a visible, measured, issue-numbered selection term is a
different thing, and that distinction is why the debt list is its own variable
rather than two more entries on `ORACLE_EXCLUDE`.

**Do not read "the rows I know about are green" as "it can be armed" — that
reading has been wrong three times out of three.** `ci.yml`'s comment on the
`test-oracle` job carries the full triage, including what a run with
`FERRO_MANIFEST` set adds, and `scripts/run_oracle_suite.sh` reproduces the whole
measurement locally in one command.

**What #1815 still does NOT close is the manifest half.** `test-oracle` arms the
flag but provisions no `FERRO_MANIFEST`, so the `mutalyzer_normalize_tests` /
`biocommons_normalize_tests` axes early-return there exactly as they do in the
plain `test` job, and nextest reports that skip path as **PASS**. A
denoted-sequence violation on those axes still cannot redden a required check.
Provisioning a reference per PR is not a small change — the nightly spends ~30 min
in `ferro prepare` and the result is deliberately uncached (#1218) — so #1815 stays
open on that half.

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

| class | why it was not a defect |
|---|---|
| output cannot be transliterated | the input states its own deleted bases (`g.1000delC`) and converts with no provider; the output (`g.1000del`) must read the reference, which the fixture does not hold |
| insertion flush against a deletion | `triples_are_disjoint` called it an overlap; `cis_apply_oracle`'s tie-break does not, and it is well defined. **Closed by #1749** — the cause was an order dependence, not a policy: the predicate carried a running 3'-most `reach`, and a **stable** descending sort left a zero-width triple and a span triple sharing one position in *input* order, so `g.[5_9del;4_5insA]` was disjoint and `g.[4_5insA;5_9del]` was not. It is now pairwise over `WriteFootprint`, and the application order breaks the tie by deletion length. #1831 reached the same class from the applier side and added the `apply_triples` sort key that puts longer deletions first among triples sharing a position — the key `splice_denoted_sequence` and `conformance::spec_corpus`'s applier already used; that key remains, and now carries the splice order alone rather than the predicate's verdict. The 233 is the size of the class when it was diagnosed, and is kept as that |
| overlap-conflicting input | an insertion *interior* to a deletion — the input denotes nothing, so there is no baseline |
| `pter`/`qter` | carry no numeric coordinate (`base: 0`); `hgvs_to_spdi` reads `base` and silently resolves `c.pterdel` to the sequence's **last** base |
| corrected `REFSEQ_MISMATCH` | `c.10dupA` where the reference reads `C`: normalization is *supposed* to change the denoted sequence here |
| `r.` payload vs DNA reference | `UGC` against `TGC` — the same bases in two alphabets |
| uncertain allele `[(…)]` | the normalizer deliberately does not clamp those |

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

An oracle fire blocks the merge: the required `Test` context is a rollup whose
`needs:` list includes `test-oracle` and `sweeps`; read that list from
`.github/workflows/ci.yml`. `generated-docs` is not in that list and still gates: the
generators run inside `spec-fixtures`, which `test` and `test-oracle` depend on, so a
generator failure skips them, and the rollup treats a skip as a failure.
