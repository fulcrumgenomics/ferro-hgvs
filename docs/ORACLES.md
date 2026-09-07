# The normalization seam oracles

Four `FERRO_ASSERT_*` flags turn every normalization in the test suite into an invariant check.
Each asks one question of every normalized description: is it a fixed point, does it re-parse, are
its coordinates in range, and does it denote the same bases the input did. This page is the home
for that story. The other three copies point here: the `ci.yml` job comments, the runner header,
and the invariant test header.

## The seam

The four flags share one call site. `Normalizer::assert_seam_oracles` runs at the single exit of
`normalize_core_checked`, so it covers every public normalization path: `normalize()`,
`normalize_with_diagnostics()`, and every `VariantProjector` axis. An oracle sees every normalized
description ferro hands back.

The four checks run in a fixed order: in-bounds, re-parse, idempotency, denoted-sequence.
Denoted-sequence runs last because it is the expensive one and the only one that reads the
reference. An out-of-bounds or unparseable output must be named as that, not as a sequence the
oracle could not apply.

Each flag is compiled out in release builds. Every call carries `#[cfg(debug_assertions)]`, so the
whole body disappears in release. Each flag is read once into a `OnceLock`, so a disabled run pays
only one atomic load. The idempotency oracle re-enters normalization to verify its own output, so a
thread-local `IN_IDEMPOTENCY_CHECK` guard breaks that recursion: the inner call skips its check.

## Running the oracles locally

Use the runner. It reproduces `ci.yml`'s `test-oracle` *armed step*: that step's four flags, over
that step's selection. It does **not** reproduce the job's compensating step, which re-runs the
`SEQUENCE_ORACLE_EXCLUDE` debt rows under the other three oracles — that step runs only in CI, so a
local run leaves those rows uncovered. See [What CI arms, and where](#what-ci-arms-and-where).

```bash
scripts/run_oracle_suite.sh                     # arm the flags and run the selection
scripts/run_oracle_suite.sh --print-selection   # print what it would run, then stop
scripts/run_oracle_suite.sh -E 'test(my_test)'  # extra args go through to nextest
```

Do not arm any of the four flags by hand over the whole suite. This command is red on `main`:

```bash
FERRO_ASSERT_IDEMPOTENT=1 FERRO_ASSERT_REPARSE=1 \
  FERRO_ASSERT_IN_BOUNDS=1 FERRO_ASSERT_SEQUENCE=1 \
  cargo nextest run --features dev   # red on main
```

It is red for two reasons, and neither is a coverage gap to close:

- A spec-corpus census counts the defect an oracle would panic on. `conformance::census::measure`
  catches the panic and files that row as `declined`, so the count reads better than the truth, and
  `conformance::census::run_census` refuses to run at all with a flag set. This applies to every
  flag: `run_census` refuses when any `FERRO_ASSERT_*` variable is set, and `measure` flatters
  under whichever oracle fires.
- Some rows pin a defect at the seam. `spec_corpus_regressions` pins the CDS-end flush-pair rows
  the denoted-sequence oracle fires on. A test that pins a defect and an oracle that fires on it
  cannot both run. This applies only to the denoted-sequence flag; the other three oracles fire on
  neither pinned-defect module.

`ORACLE_EXCLUDE` names those modules, and the armed CI job negates it. See
[What CI arms, and where](#what-ci-arms-and-where).

The runner reads the whole `-E` selection and the flag set out of `ci.yml`. It copies nothing, so
neither can drift into a second copy. `tests/it/oracle_exclude_invariant.rs` re-derives the same
selection in Rust from different anchors and compares the two, so a hand-built copy that drifts from
the file fails that test.

## Idempotency oracle

`FERRO_ASSERT_IDEMPOTENT=1` asserts that `norm(norm(x)) == norm(x)`. Every test that normalizes
becomes an idempotency check. Its blind spot: it verifies by re-normalizing its own output, so it
cannot judge an output that fails to parse.

## Re-parse oracle

`FERRO_ASSERT_REPARSE=1` asserts that `parse_hgvs` accepts a normalized description, when
normalization is what broke it. The exemptions are a closed list of four:

- `0` and `?` are legal whole-allele outputs that `parse_hgvs` rejects standalone because it wants
  an accession.
- An empty allele (`[]`), which only direct construction reaches; the projector's own tests build
  one to pin that it declines.
- A non-flanking genomic insertion, the projection pivot: its coordinates are sound but its spelling
  is not one HGVS admits, so the projector withholds the reported genomic axis instead
  (`non_flanking_genomic_insertion_anchor`).
- A non-coding downstream position (`n.*N`), which parse refuses in every mode while
  `TxPos::downstream` stays public API, so only a Rust caller reaches it; `noncoding_zone_marker`
  keys the exemption on the AST.

Keep the list whole. If you widen it, that is the signal to fix the producer. Its blind spot:
`parse_hgvs` holds no provider, so a well-formed spelling that denotes the wrong bases is valid to it.

## In-bounds oracle

`FERRO_ASSERT_IN_BOUNDS=1` asserts that no coordinate a normalized description names is past the end
of its own sequence. The rules live on the doc comment of `merge::first_out_of_bounds_coordinate`.
Read them there. Not covered: protein axes, and an inserted-range payload (`g.10_11ins[20_30]`).

The oracle exists because this defect class was found by hand, one shape at a time, in #1274, #1343
and #1307 before #1353 asserted it at the seam.

Its blind spot: an out-of-range coordinate that is a fixed point passes idempotency, so idempotency
does not catch it.

## Denoted-sequence oracle

`FERRO_ASSERT_SEQUENCE=1` is the only oracle that asks what the output means. It applies the input to
the reference, applies the output, and asserts the bases agree. It applies both descriptions over
the union of their spans, in one fetch, so a shared frame does not report a 3'-shift as a difference.

The other three are all form questions, and a wrong sequence passes all of them. It is a fixed point,
so idempotency is satisfied. It parses, and `parse_hgvs` holds no provider, so re-parse cannot know.
Its coordinates exist, so in-bounds is satisfied.

The class was found by hand in #1254, #1281, #1290, #1304, #1308, #1312, #1592 and #1600 before #1615
asserted it at the seam, and #1592 and #1600 each record the other three oracles passing on their
reproducer.

The applier is not the normalizer. `spdi::compare_denoted_sequences` reaches the bases through
`hgvs_to_spdi` and an SPDI splice, the same walk `apply_to_reference` and
`tests/it/common/cis_apply_oracle.rs` use, so nothing here agrees with the output merely because
normalization produced it. `EquivalenceChecker` is not usable for this: it normalizes both sides,
which is circular. For why `hgvs_to_spdi` and the normalizer read a `c.` position on the same flat
transcript axis, see the ruling record `c-and-n-positions-are-flat-transcript-offsets` in
`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`. The regression guard for that defect
is `CDOT_GAP_JUNCTIONS` in `tests/it/normalization_transcripts_exon_contract.rs`, whose doc comment
says why the NM_033517.1 record must keep its gap.

A side that cannot be applied is counted, not silently passed. A skip that reads as a pass is the
exact failure mode this oracle exists to remove:

| case | verdict |
|---|---|
| both apply, bases agree | pass, counted in `compared` |
| both apply, bases differ | fire |
| the output denotes no sequence while the input does | fire. Two members can claim one base, which is worse than a wrong sequence, so it is never a skip |
| the input denotes no sequence (a trans allele, a `REFSEQ_MISMATCH`, an edit SPDI cannot carry) | skip, counted in `skipped`. There is no baseline |
| the two name different accessions | skip, counted |
| the union window exceeds `MAX_APPLY_WINDOW`, or the provider cannot serve it | skip, counted |

`normalize::denoted_sequence_oracle_counts()` returns `(compared, skipped)` process-wide. Read it
before you trust a green oracle run, because zero comparisons and zero faults look the same from the
outside.

`tests/it/issue_1615_denoted_sequence_oracle.rs` is the oracle's own regression guard. It pins each
recorded wrong output and asserts the predicate fires without re-normalizing, so it keeps testing the
oracle rather than the fix. Its other half is a negative control where a legitimate re-spelling must
stay silent.

The first run of this oracle over the suite raised many fires, and all but a few were false. Each row
below is a class, and the reasoning is on the code:

| class | why it was not a defect |
|---|---|
| output cannot be transliterated | the input states its own deleted bases and converts with no provider; the output must read a reference the fixture does not hold |
| insertion flush against a deletion | the disjointness predicate called it an overlap; the applier's tie-break does not, and it is well defined. See #1749 and #1831. The 233 is the size of the class when it was diagnosed |
| overlap-conflicting input | an insertion interior to a deletion; the input denotes nothing, so there is no baseline |
| `pter`/`qter` | they carry no numeric coordinate, so `hgvs_to_spdi` resolves the position to the last base |
| corrected `REFSEQ_MISMATCH` | normalization is supposed to change the denoted sequence here |
| `r.` payload against a DNA reference | the same bases in two alphabets |
| uncertain allele `[(…)]` | the normalizer deliberately does not clamp those |

The two sides of a comparison do not need the reference equally. Any check that reads "I could not
derive it" as "it is wrong" fires hardest where that asymmetry is largest.

## What CI arms, and where

Each column is a flag. Each row is a job or a step. `yes` means the flag is set.

| job or step | IDEMPOTENT | REPARSE | IN_BOUNDS | SEQUENCE |
|---|---|---|---|---|
| `test-oracle` armed step | yes | yes | yes | yes |
| `test-oracle` compensating step | yes | yes | yes | no |
| `sweeps` | yes | yes | yes | yes |
| `censuses` armed step | yes | yes | yes | no |
| `censuses-plain` | no | no | no | no |
| `test` | no | no | no | no |
| `soak` | no | no | no | no |
| nightly | yes | yes | yes | yes |

`ORACLE_EXCLUDE` states which instrument may be armed while another is measuring. It is not a
coverage exemption: every module it names runs unarmed in the plain `test` job. Two census mechanisms
sit behind it. `conformance::census::measure` catches a panic and files that row as `declined`, which
flatters the count; `spec_conformance_axis` calls it. `conformance::census::run_census` refuses
outright when any `FERRO_ASSERT_*` is set and returns `CensusError::OracleArmed`;
`conformance_census_runs` calls it, and `conformance_census_instrument` reaches it through the
`ferro-benchmark` binary.

Two pinned-defect modules sit on the same list. The denoted-sequence oracle fires on rows that
`spec_corpus_regressions` pins, the CDS-end flush-pair class, and a fire reddens that test rather than
emptying any sweep or count. The idempotency, re-parse and in-bounds oracles fire on neither
`spec_corpus_regressions` nor `defect_non_idempotent_outputs`. `defect_non_idempotent_outputs` fires
under no flag; it stays on the list because `tests/it/oracle_exclude_invariant.rs` requires every
module that reads the spec corpus to be there.

`SEQUENCE_ORACLE_EXCLUDE` is a debt list. It withholds rows from the denoted-sequence oracle and
from nothing else: a module whose purpose is to pin a defect, and a gate whose corpus rows fire. Read
the rows off `ci.yml`, each beside the open issue that retires it.
The `test-oracle` compensating step re-runs exactly those rows under the other three oracles, so this
file's oracle coverage is a strict superset of what it was, not a trade. Suppressing a row at the seam
would hollow out the oracle; a visible, issue-numbered selection term does not, which is why the debt
list is its own variable.

Only a selection-wide armed run over a job's own `-E` says whether a flag can be armed. "The rows I
know about are green" has never been sufficient here. `censuses` arms three flags and not the fourth,
deliberately: its flag set is a strict subset of `test-oracle`'s. Do not restore parity by copying
`FERRO_ASSERT_SEQUENCE` down without first measuring it over `ORACLE_ONLY_FILTER`'s modules. And
`test-oracle` arms the fourth flag but provisions no `FERRO_MANIFEST`, so the reference-aware axes
early-return there and nextest reports that skip path as PASS: a denoted-sequence violation on those
axes cannot redden the required check, and that manifest half stays open (#1815).

## The nightly and the merge gate

In PR CI, `test-oracle` and `sweeps` carry no `continue-on-error`, so an oracle fire turns the job red.

The nightly reference-aware job carries `continue-on-error`. Its purpose is to surface drift in the
xfail report, not to gate. The corpus runner wraps normalization in `catch_panics`, so an oracle fire
lands in the uploaded xfail artifact as a failing case. Read a nightly oracle fire out of that report,
not from the workflow conclusion.

A nightly issue can still be about an oracle. The #1998 diff step compares the run's failing set
against the committed baseline. It carries no `continue-on-error`, so any drift fails the job and
`report-failure` opens a tracking issue, and a new oracle fire changes the failing set, so it reaches
that gate. A separate job-summary step prints the armed reproduction recipe, keyed on the test step's
`outcome`.

An oracle fire blocks the merge. The required `Test` context is a rollup. Its `needs:` list includes
`test-oracle` and `sweeps`. Read that list from `.github/workflows/ci.yml`, not from this page.
