# Interpreting the HGVS Recommendations

> **Under construction.** Pages are added one HGVS spec page at a time; this section will fill in
> as they land.

## What it is

The [HGVS nomenclature recommendations](https://hgvs-nomenclature.org/) are the standard for naming
sequence variants. They are large and carefully written, but — like any natural-language standard —
in places they leave a question open, or two passages can be read as pointing different ways. A
normalizer still has to emit exactly one string. Where ferro has had to settle such a question, we
record the decision (and, where a passage looks inconsistent, file it upstream with the SVD Working
Group).

This section mirrors the recommendations page for page and records, per clause, how ferro reads
each one. Each page gives, for a span of the spec:

- the **spec text** it interprets, quoted verbatim;
- **ferro's reading** of that clause;
- **example spellings**, each with the form ferro **normalizes** it to and a verdict describing that
  output against the recommendations (see below);
- which spellings **converge** on one output;
- the **reason** ferro reads the clause the way it does, and a link to the governing decision.

## How to read the verdicts

The recommendations are HGVS's, not ferro's. ferro's job is to take any input and produce the form
the recommendations prefer. The verdict describes **ferro's output** against that goal:

- **recommended** — ferro's output is the form the HGVS recommendations prefer. This is the target,
  whether the input was already that form or ferro normalized it there.
- **conformant** — ferro's output is a valid HGVS string but not *yet* the recommended form. This is
  a known ferro limitation, and it links to the issue tracking it.
- **refused** — the input is not valid HGVS; ferro rejects it in strict mode. Refusing is the
  correct behavior.
- **bug** — ferro's output is not valid HGVS. This should never happen; where a page shows one, it
  links to the open defect.

## What it is not

It is not a second copy of the rules. The canonical ruleset lives in the
[normalization rules](../reference/normalization-rules.md) page, and
every adjudication is recorded once in the **ruling ledger**
([`hgvs_spec_normalization_overrides.json`](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/tests/fixtures/grammar/hgvs_spec_normalization_overrides.json),
rendered for reading as
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)).
These pages **draw their explanations from the ledger** rather than restating them — a rule
written in two places is a rule that drifts.

## How it stays honest

Every example is executed through ferro by a build check: each spec quote is verified against the
pinned spec checkout, each spelling is parsed (and strict-refused where marked **refused**), and each
normalized output is asserted against a real GRCh38 reference. A wrong verdict, a moved quote, or a
changed normalization fails the build — so a page cannot claim behavior ferro does not have.

## Pages

### DNA

- [Substitution](recommendations/DNA/substitution.md)
