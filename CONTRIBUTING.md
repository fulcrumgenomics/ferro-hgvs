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
git submodule update --init assets/hgvs-nomenclature   # the pinned HGVS spec
scripts/fetch-test-fixtures.sh                         # the 156 MB bulk corpora
cargo build
cargo test --features dev
```

`scripts/fetch-test-fixtures.sh` is not optional if you want the full suite. Four
large HGVS corpora (156 MB) are **not in the git tree** — they are assets on the
`test-fixtures-v1` release, because Git LFS bandwidth is metered even for a
public repository and ten CI jobs per run exhausted the budget. The script
downloads them over plain HTTPS (no token) and verifies every one against
`tests/fixtures/CHECKSUMS.sha256`, refusing anything that does not match.

Skipping it is safe but lossy: the four suites that read those corpora
(`clinvar_hgvs_tests`, `cmrg_exhaustive_tests`, `paraphase_exhaustive_tests`,
`normalize_axis_preserving`) return early and **report PASS** rather than
failing. See [Testing](#testing) for the flag that turns that into a failure, and
`tests/fixtures/README.md` for provenance and regeneration.

### Python Bindings Setup

```bash
uv sync --group dev
uv run maturin develop --features python,extension-module
uv run pytest
```

After modifying Python dependencies in `pyproject.toml`, run `uv lock` and commit
the updated `uv.lock`. CI uses `--locked` and will fail if the lockfile is out of
sync with `pyproject.toml`.

#### Why `extension-module` is a separate feature

**`extension-module` is a SEPARATE Cargo feature, not part of `python` (#2046).** `pyo3` is
declared with `features = ["abi3-py310"]`, and `pyo3/extension-module` is pulled in by the
`extension-module` feature instead.

The split exists because `extension-module` deliberately does *not* link against libpython (the
CPython symbols resolve against the host interpreter at import time). With it enabled, PyO3 emits
no libpython link line on Linux/macOS (`is_linking_libpython_for_target` in `pyo3-build-config`),
so a plain test/bin executable fails to link its undefined Python symbols — which is exactly why
`src/python.rs`'s pure-Rust `#[cfg(test)]` tests ran nowhere in CI. Keeping `extension-module`
out of `python` lets `cargo test/nextest --features python --lib` link libpython and actually run
them; CI does this in the `Python Wheel Test` job (which has a shared-library Python available).

**So every wheel/develop build must ask for `extension-module` BY NAME — `[tool.maturin]
features` does not supply it for free.** That list in `pyproject.toml` is a **default**, not an
addition: measured with maturin 1.13.3, a `--features` flag on the maturin command line
**replaces** it rather than merging with it, so `maturin build --features python` hands cargo
`--features python` alone and the pyproject entry is never consulted. The entry is therefore
load-bearing only for builds that pass no `--features` of their own (the PEP 517 backend — `pip
install .`, `python -m build`, `uv`), and every explicit invocation in this repo names the feature
itself: `ci.yml`'s "Build release wheel", `release-wheels.yml`'s three `maturin-action` blocks,
`pyproject.toml`'s two `build-*` tasks, and `CONTRIBUTING.md`.

Do not "simplify" those call sites back to `--features python` on the grounds that the wheel
still works without it. It does still work — maturin sets `PYO3_BUILD_EXTENSION_MODULE=1`, and
pyo3 0.28's `is_extension_module()` honours that env var independently of the Cargo feature,
whose doc comment there calls it *"(deprecated)"*. Measured: wheels built both ways are
indistinguishable (no libpython link line, 99 undefined `_Py` symbols each). But that makes the
linkage a property of **which maturin ran**, and `ci.yml` installs an unpinned
`maturin>=1.0,<2.0` while `release-wheels.yml` takes `maturin-action`'s default — so a maturin
predating that env var would build a libpython-linked wheel silently. Naming the feature restores
exactly the cargo feature resolution this repo had before the split, and costs no version floor.

#### Why the binding tests are scoped

The `Python Wheel Test` job reports a **required** status check, and the lib suite is not safe to
run in one process. `cargo test --features python --lib` runs all 4240 lib tests under libtest's
thread pool, where `normalize::merge::tests::a_declined_sequence_first_partition_is_counted`
asserts deltas on process-global `AtomicU64` counters (`SEQFIRST_*` in `src/normalize/merge.rs`)
that any concurrent test reaching `partition_block_for_rule` corrupts, and
`parallel::tests::test_stress_concurrent_throughput` is a timing test sensitive to
oversubscription. Measured at `--test-threads=32`: **4 of 6 runs failed** — i.e. an unrelated PR's
required check flaking. nextest is process-per-test, so the corruption is not expressible at any
thread count (10/10 green at the same thread count), and the selection keeps the job to the tests
it exists to run: 10 tests in 0.02s against 4240 in ~14s.

`--no-tests=fail` is the non-vacuity guard, and the other reason this is nextest. A selection
matching nothing exits **4** here (`error: no tests to run`); `cargo test -- python::tests`
reports `ok. 0 passed; 4240 filtered out` and exits **0**, so a typo in the filter would silently
delete the coverage without reddening anything.

**Still use `maturin`, not `cargo build --features python`, for an importable module.** A plain
`cargo build` produces a cdylib that is not a wheel-ready abi3 extension, and historically dies
linking on macOS unless `-undefined dynamic_lookup` is passed — PyO3 ships those flags in
`pyo3_build_config::add_extension_module_link_args()` (Darwin-only) but a crate must call them
from its own `build.rs`, and this crate has none. Maturin passes them itself. Use
`cargo check`/`cargo clippy --features python` to verify the bindings compile, and `maturin
develop --features python` whenever you need a module you can import (e.g. to run
`pytest tests/python/`).

The commands, verbatim from the former CLAUDE.md:

```bash
cargo check --features python        # Typecheck the bindings (fast; does NOT link)
cargo clippy --features python       # Lint the bindings (also does NOT link)
# Run src/python.rs's own Rust #[test]s (#2046). SCOPED, and nextest not
# `cargo test` — see "Why the binding tests are scoped" below.
cargo nextest run --features python --lib --no-tests=fail -E 'test(python::tests::)'
```

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

This section is the mechanism for rule 7 of the project's
[normalization rules](docs/src/reference/normalization-rules.md); read that section for what has to be
disclosed, and this one for how to declare it.

**Disclosure is the shipped guarantee — stability is not.** A downstream consumer
keys its read counts on the normalized HGVS string, so a canonical form that
churns between releases silently re-buckets everything already stored: a
re-normalization of that consumer's whole library, not a bugfix they can ignore.
That is why a move must never be silent — and *silent* is the part the project
promises to avoid.

It does **not** promise the string stays put. Every fix to rule 3 (*confluent*)
picks a winner between two spellings, and whichever side loses is a
representation change for someone; a project that refused to move a string could
not converge on one. So the ranking is **confluence, then disclosure, then
stability**. Where the spec permits either form *and confluence is unaffected
either way*, the already-shipped representation is a reasonable tiebreaker — but
it is a tiebreaker of last resort, never an argument against a confluence fix.

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
  modules; nothing under a watched directory is touched, and no existing
  expectation was re-blessed.
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
**The number that counts is the count, never the denominator** —
`3 rows of 500,004 move`, `3 of 500,004 rows move` and `3 of 500,004 corpus rows move` all disclose
three moved rows, and the rule reads three from each.

This sentence used to say "the pattern requires the count to sit immediately before `rows`", and
that was the #1647 defect rather than a description of it: in `0 of 950 rows move` the number
immediately before `rows` is the **denominator**, so the very example offered here as passing was
rejected, and the documented form was a merge blocker. Do not restate the rule positionally.

A declining trailer is
excluded from the changelog's **Representation changes** section, so declaring it
costs the reader nothing while leaving the judgement on the record. CI enforces
the distinction: a change touching `src/normalize/`, `src/hgvs/`, `src/spdi/`,
`src/project/`, `src/reference/` or `src/error_handling/` with no trailer at all
fails the `Representation change declared` check, because an absent trailer is
indistinguishable from an unconsidered one.

The last two joined that list on measurement (#1853). `src/reference/` decides
which bases a description resolves against and `src/error_handling/` decides the
accept/reject boundary, and both have moved consumer-visible output while the
gate stayed silent — between them they account for every newly-normalizing row
in the v0.13.0 cycle, whose disclosure had to be reconstructed by hand. Two
neighbours are deliberately **not** watched: `src/conformance/`, which is
measurement and adjudication code that cannot move ferro's output, and
`src/data/`, whose every measured disclosure also touches `src/reference/` and is
therefore already covered.

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

**Do not name a release version in the trailer.** Which release carries your
change is not knowable when you write it: `release-plz.toml` sets
`features_always_increment_minor`, so what a `feat:` cuts depends on the base
version at release time, a `feat:` landing before yours moves it, and a cycle can
be re-cut. Say **"to reproduce output from before this change"**, not "to
reproduce pre-vX.Y.Z output". The version-free form is correct under every cut
and loses nothing, because the changelog heading the trailer is quoted under
already names the release the reader is holding.

This is the house form because it has already been got wrong the other way:
#1835's trailer told consumers to pin `FERRO_PARTITION=live` to reproduce
**pre-v0.15.0** output, for a change shipping in **v0.14.0** ([#1886]) —
which reads to a consumer as "this lands next cycle" rather than "this is the
release you are reading", and could not be fixed in the trailer once merged.

Two related habits (see [Adjudications](#adjudications) below):

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

One limit of that section is worth knowing before you rely on it: it carries the
commit **subject** only — the trailer's own text does not reach the changelog
([#1556]), so the four facts above live in your PR description and a reader has
to follow the link for them.

A declining commit (`Representation-Change: none`) is grouped by its conventional
type, so a `fix:` that correctly declined is listed under **Fixed** rather than
buried in **Other** — a `commit_preprocessor` in `release-plz.toml` neutralizes
the declining trailer so the commit falls through to its type ([#1557]).

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

The route that works is not through the template: read the trailer from the raw commit message
with `git log`, and attach it to the rendered bullet afterwards — which is what
`check_changelog_grouping.py` already does to attribute trailers, so the script imports its
`trailer_value` rather than keeping a fourth copy of the column-0 rule.

The checker is
`scripts/check_representation_change.py` and the audit is
`scripts/check_changelog_grouping.py`; the decline vocabulary is pinned against
`release-plz.toml` and `CONTRIBUTING.md` by tests in
`tests/python/test_representation_change_trailer.py`, because three copies of one list drift
silently.

[#1556]: https://github.com/fulcrumgenomics/ferro-hgvs/issues/1556
[#1557]: https://github.com/fulcrumgenomics/ferro-hgvs/issues/1557
[#1886]: https://github.com/fulcrumgenomics/ferro-hgvs/issues/1886

Two things follow for representation-stability work specifically. The repository's doctrine once
read as "do not move a normalized output", while the downstream filer has said instability is
acceptable **provided it is declared a breaking change** — different bars, and `docs/src/reference/normalization-rules.md`
chooses between them: disclosure is the obligation, so the filer's bar is the project's bar. And
Mutalyzer is not a spec oracle:
its normalizer re-derives a description by minimizing a **weighted description length** with
constants dated 2014, and has no separation rule at all, so a separation disagreement with it is
two objectives meeting rather than evidence it knows something the spec does not. The full
forensics are in the `adjudication-precedence-order` record and on `MIN_PIECE_SEPARATION` in `src/normalize/merge.rs`.

### Adjudications

**Policy.** When you decide what the *correct* normalization behaviour is — not how to
implement it — that decision lands in a committed test or ruling record **in the same PR**, with
a comment stating the question, the ruling, and the authority. An adjudication that lives only in
a PR description, an issue comment, or a working document is lost: the next person re-derives it
from scratch, and often re-derives it differently.

This is cheap because the machinery already exists (the five committed guards are tabulated
in `docs/TESTING.md`, under "Generated spec fixture"). The policy is about using it consistently, not about building anything.

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

- **adjudicated-correct** — pin the exact expected output and cite the clause; a real guard that
  fails when behaviour regresses away from a decided answer.
- **adjudicated-deviation** — pin it as a deviation via `KNOWN_DIVERGENT_INPUTS`, so a fixed
  deviation cannot rot in the list unnoticed.
- **undecided** — a first-class state; the generator refuses one that names a governing clause, so
  an open question cannot smuggle in a ruling nobody made. Prefer an honest `undecided` record to
  no record.
- **house-choice** — decided, and ours rather than the spec's: a `rulings` record with a
  `house_choice` object, made under rule 5's silent limb or rule 6 of
  `docs/src/reference/normalization-rules.md`, naming no governing and no deviated-from clause;
  never citable as conformance. It must say what was considered and rejected.

**A test that merely pins today's output is not an adjudication record.** It is a change
detector, and this repo has already been bitten by the difference —
`pinned_v21_normalization_behavior` compares ferro against itself (see `docs/TESTING.md`, Generated spec fixture,
on stale-local-artifact detectors). What makes a record an adjudication is the *authority*: an exact
`file:line` into `assets/hgvs-nomenclature`, a named Mutalyzer measurement, or an explicit
"undecided, and here is why". Without one, you have frozen the current behaviour including
whatever is wrong with it.

**Deviating from the reference implementation needs a record, not just a rationale.** Precedence
is **spec-explicit > Mutalyzer > our judgement**. Where the spec is explicit and Mutalyzer differs,
that is a deliberate divergence and it gets a record saying so — otherwise the next person
measures Mutalyzer, finds a mismatch, and "fixes" a conformance decision.

**Record what was refuted, not only what was decided.** A measurement that kills a plausible
belief is worth as much as the ruling itself, because the belief will recur —
`MIN_SEPARATION_NO_FRAME`'s doc comment in `src/normalize/merge.rs` is a worked example of
recording such a refutation.

**Cite the clause exactly, and quote it.** Do this in prose comments too: `general.md:33`, not
"the separation rule". A clause's directory is its jurisdiction — a claim about an `r.` axis needs
a clause under `RNA/`, since a `DNA/` clause cannot scope `r.`.

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


#### Never hand-edit `tests/it/clause_ruling_index.rs`

The file is committed, so it looks hand-maintainable. It is not: it carries a rendered index of
every clause the ruling records cite, and `the_rendered_index_is_current` compares that block
against what the code generates. Editing it by hand is transcription, and transcription of a
~100-row block is lossy in a way that reads as a small diff — two attempts at it were made here
and both dropped rows silently.

Capture what the test prints **programmatically** instead. The pattern that works: add a throwaway
test that writes `render_index()` to a file, run it, restore the source byte-for-byte from a
pre-edit copy, then splice the generated text in on its `BEGIN`/`END` markers. Nothing is retyped.

Note this is a *different* trap from the gitignored spec artifacts (`docs/TESTING.md`, Generated spec fixture). Those regenerate on
demand and a stale one fails loudly. This one is committed, so a lossy edit survives review as an
ordinary diff and only fails on the one assertion that re-renders it.

### Changelog

`CHANGELOG.md` is generated by [release-plz](https://release-plz.dev/) from the conventional-commit history when a release PR is opened — do not hand-edit it in feature PRs. Two PRs both adding bullets under `## [Unreleased]` will always conflict on the same line range, so leaving the changelog to tooling keeps feature PRs mergeable until release time. Write a clean conventional-commit subject (`feat(scope): …`, `fix(scope): …`, etc.) and release-plz will pick it up.

That it is generated is also why a `Representation-Change:` trailer is the way to get such a note into the changelog *at release time*: release-plz regenerates the pending section of an open release PR on every run, so an entry hand-written there is lost the next time the PR updates. `release-plz.toml`'s `commit_parsers` route trailered commits into a **Representation changes** group ahead of the usual ones, and its `pr_body` template surfaces that group on the release PR itself. The parsers match the commit **footer**, so the trailer must survive squashing — put it in the PR description as well, since that is what GitHub uses for the squash commit body.

**The bullet release-plz writes is the commit subject and nothing else, so the trailer's own text is attached afterwards** (#1556). `scripts/inject_representation_disclosure.py` reads each trailer out of the raw commit message with `git log` and quotes it under its bullet, which is what makes the release checklist's "previously rejected or previously accepted?" question answerable from the changelog rather than from every linked PR. Run it on the release PR:

```bash
python scripts/inject_representation_disclosure.py          # attach; idempotent
python scripts/inject_representation_disclosure.py --check   # what CI asks
```

CI runs the `--check` form in the `Changelog grouping audit` job, so a release PR stays red until it has been run. Re-run it after any push that makes release-plz regenerate the pending section — a second run over an already-attached section changes nothing.

Two consequences for how you write the trailer. Keep it **at column 0 and undecorated** — no backticks, no `>`, no `-`/`*`/`+`, no `#`, no emphasis, and the separator is a hyphen. Anything else and the line is prose: an indented `Representation-Change:` is a continuation of the line above it to git, to release-plz and to both checkers, and a code-spanned one is invisible to all three at once. This is not pedantry about formatting — you would have declared something and nobody would be told, which is the exact failure the trailer exists to prevent — so the checker now refuses a line that looks like a declaration it cannot read, names the character responsible, and does so on every PR rather than only on the ones touching a watched directory. To show the form as an *example* rather than declare it, indent the line — and keep your real, unindented declaration above it, because that is what tells the checker the indented one is a quotation. An indented line on its own is reported as a near miss rather than read as an example, which is the right answer: nothing distinguishes it from a declaration that will never be read. And know that "put it last" is advice the tooling cannot rely on: GitHub appends CodeRabbit's summary after your description, so the script stops the value at the next column-0 trailer token, Markdown heading or HTML comment rather than at the end of the message. Do not put a heading inside the disclosure.

A section that has already been released is a different matter: a release only prepends to `CHANGELOG.md`, leaving earlier sections untouched, and the GitHub release body is never revisited once its tag exists. So a retrospective note — a consolidated representation summary for a release whose commits predate this convention, say — can be added afterwards, by editing the released section and `gh release edit <tag>`. That is the escape hatch for a cycle that shipped without trailers; it is not a substitute for declaring changes as you make them, because nobody can reconstruct the rejected-vs-accepted distinction after the fact.

**Correcting a disclosure that is already merged is a third case, and hand-editing `CHANGELOG.md` is the wrong tool for it.** A trailer cannot be corrected once its commit is on `main` — that would mean rewriting released history — and an edited changelog copy fails twice over: release-plz discards it the next time it regenerates the pending section, and `inject_representation_disclosure.py --check` reports it as drift against the trailer it no longer matches, which is red CI. The durable surface is the injector itself, because it re-renders the block on every run. Register the correction in its `EDITORIAL_CORRECTIONS` table, keyed by PR number, and re-run the script:

```bash
python scripts/inject_representation_disclosure.py
```

It is appended under the bullet after the trailer's own words, inside the same blockquote, and `--check` then demands it. One limit to know before you reach for it: the injector *judges* a block that is already on disk rather than rewriting it, so registering a correction against an already-rendered bullet is reported as drift. That resolves itself in the pending section, which release-plz regenerates on every push to `main`; a **released** section is never regenerated, so there you delete the stale block from `CHANGELOG.md` and let the re-run write it back with the correction attached. Two constraints, both load-bearing. Each paragraph must **open by naming itself editorial and citing the issue that raised it** — the release checklist turns on whose words a disclosure carries, and inside the block the label is the only thing that answers it. And a correction states a **fact the trailer got wrong or left out**; it is not a place to re-word, summarise or improve a disclosure, which is what keeps "the trailer's own text" a claim a reader can rely on. Append a paragraph rather than rewriting one, so two corrections to one disclosure conflict visibly instead of silently overwriting each other.

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

### Bulk corpora: a skip that reads as a pass

The same trap as the conformance axis below, from a different direction. The four
large corpora fetched by `scripts/fetch-test-fixtures.sh` are release assets
rather than git objects, and the suites reading them return early when a fixture
is absent — reported as PASSED, not skipped. So a run without them looks clean
*and faster*, which is the worst way for coverage to disappear.

`FERRO_REQUIRE_BULK_FIXTURES=1` turns every such skip into a failure naming the
missing path and the command that fixes it:

```bash
scripts/fetch-test-fixtures.sh --verify   # are they present and correct?
FERRO_REQUIRE_BULK_FIXTURES=1 cargo nextest run --features dev --lib --test it
```

CI sets it in every job that fetches the fixtures (`test`, `test-oracle`,
`coverage`, `external-validation`), so a fetch that failed reddens the job instead
of quietly deleting the ClinVar, CMRG and Paraphase coverage. Leave it unset
locally unless you have fetched the corpora.

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

### Assert the property. Measure the count. Never let a count BE the property

**A guard that pins a number is a change detector for that number.** It is a guard for the
property only for as long as the two agree, and nothing makes them agree. When they drift, the
guard follows the number, stays green, and reads as coverage.

Three shapes, in ascending order of how convincing the wrong answer looks.

**1. A count restated instead of imported.** `the_corpus_emits_a_block_past_the_split_cap` in
`examples/dump_normalized_corpus.rs` used to open `const SPLIT_CAP: u64 = 1024;` and assert the
corpus builds a block wider than it, while the cap the normalizer applies
(`MAX_SPLIT_BLOCK` in `src/normalize/merge.rs`, equal to `MAX_CANONICAL_WINDOW`) had been raised to 4096 by #1899. The
corpus's long cores are `[1024, 1100]` and were chosen — the comment says so — as "the pair a
change to the cap moves between". The cap moved, both cores landed on the same side, and
`1100 > 1024` kept the guard green. It even carries a comment instructing the reader to update
the literal if the constant moves; **that comment is the defect, not a mitigation.** A guard
that restates the value it guards cannot observe that value changing. Import the constant.
`dump_normalized_corpus` is a genuine `[[example]]`, so it links the library as an **external
crate** (`use ferro_hgvs::…`) and can see only the library's `pub` API — **`pub(crate)` is
unreachable from there.** A private item must be made `pub`, with `#[doc(hidden)]` on the
re-export to keep it off the public documentation surface; `src/lib.rs`'s `ShuffleDirection`
re-export is the worked pattern.

**Do not gate such a re-export behind `dev`.** The tree has already decided that, for this very
constant, and `src/normalize/mod.rs` gives the reason in as many words: "**NOT gated on `dev`**:
… a constant present in only some builds re-creates that gap for anyone building without the
feature." A feature gate buys back the documentation surface at the price of reintroducing the
restated literal wherever the feature is off — which is the failure this shape is about, moved
rather than removed. Where a gate genuinely is unavoidable the consumer must **refuse** in a build
that lacks the item rather than report a zero: `partition_blocks_cut` is `#[cfg(debug_assertions)]`
and `measure_spec_conformance_per_arm.rs` does exactly that.

Making the item `pub` is cheaper than the class of failure it buys out of. Filed as #1925, which
made `MAX_CANONICAL_BLOCK` `pub` and imported it in that example.

**2. A zero whose instrument cannot vary the property.** `examples/dump_normalized_corpus.rs` has
been blind in each of the ways its own header now catalogs, and that header records the
governing fact: **fixing one blindness does not reveal the next.**

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

### Adding or changing a generator under `examples/`

In one sentence: a generator must account for what it dropped.

Generators build the corpora, fixtures and tables the rest of the suite adjudicates
against, and they share one failure mode: a fallible step whose failure is
representable as a legitimate value — `unwrap_or_default()`, `else { continue }`, a
discarded `Result`. The dropped population is never counted, so a partial run and a
clean run write indistinguishable artifacts. This is the same trap as the
"a corpus zero is a claim about the corpus, not about the change" rule in
[Assert the property. Measure the count. Never let a count BE the property](#assert-the-property-measure-the-count-never-let-a-count-be-the-property),
one level down: there the misleading number is *reported*, here it is never computed.

Two rules, with two very different amounts of machinery behind them —
`tests/it/generator_completeness.rs` (#1550) has a guard for each, and **neither guard is
exact end to end**. Rule 2 parses the manifest but only token-scans the source; rule 1 is a
substring heuristic on both sides. It matters that you know which half you are leaning on:

1. **A generator that writes an artifact routes its population through
   `CaptureLedger`** (`src/conformance/completeness.rs`). Record a success or a drop
   for every record — `ledger.record(id, fallible())` is a drop-in for
   `let Ok(v) = fallible() else { continue };` — and call `finish()` before writing.
   It returns a `Result` and refuses on any shortfall, including a pass that attempted
   nothing. An expected shortfall is waived by naming it:
   `finish_with(Allowance::at_most(2, "bare LRG_ ids have no versioned index entry"))`.
   The counts are serializable, so stamp them into the artifact and let the consuming
   test assert on them.

   Route the **last** fallible step before the write, not the first. The ledger accounts
   for the step you hand it and nothing else — it never sees the artifact, so a second
   unaccounted `else { continue }` downstream will happily let you stamp
   `{"attempted":10,"succeeded":10,"dropped":0}` onto a file holding five rows.

   If you are not adopting the ledger, add your generator to `LEDGER_EXEMPT` in that
   test with a reason. That is a legitimate answer — what the check rejects is silence.
   The list is shrink-only: a row that is no longer needed fails the test.

   **The guard behind this rule is a heuristic floor.** It asks whether your source
   contains one of a few write idioms (`fs::write(`, `File::create(`,
   `OpenOptions::new(`, `from_path(`) and, if so, whether your entry point carries a
   `use` line naming `conformance::completeness`, the word `CaptureLedger`, and a
   `finish` call — all three in that one file. Both sides are substring scans: write a
   file some other way and the ratchet never notices, import the ledger through a
   fully-qualified path or a rustfmt-wrapped `use` group and it reads as un-adopted,
   and "routes its population through" is not something a grep can check. So treat a green test as "the
   question got asked", not as "this generator is accounted for" — the accounting is
   yours and your reviewer's to get right. Nothing checks that any artifact actually
   carries its `CaptureCounts`, either; stamping them is a convention.

2. **A generator carrying `#[cfg(test)]` sets `test = true` on its `[[example]]` or
   `[[bin]]` entry in `Cargo.toml`.** Cargo does not build a target's tests unless it
   opts in, so without it they silently never run — which reads as coverage while
   providing none. **The manifest side is exact and the source side is a token
   scan**: `test = true` is read out of parsed TOML, but `#[cfg(test)]` is
   detected by looking for that literal token, so gating your tests as
   `#[cfg(all(test, feature = "dev"))]` slips past it. `#[path]` includes are
   followed transitively (`test = true` is a per-target flag, so a `#[cfg(test)]`
   in a shared `examples/common/` module is dead in every target that does not opt
   in). A plain `mod helper;` include is the one boundary — the repo uses `#[path]`
   throughout.

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
