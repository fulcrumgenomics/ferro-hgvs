# How to read a page

Every interpretation page has the same shape. This page explains that shape and the handful of
recurring terms, once — the per-page write-ups link back here instead of repeating it.

## The layout of a page

Each page mirrors one section of the official HGVS recommendations. Under a heading naming the
**exact spec lines** it covers (e.g. `` `substitution.md:5` ``), you'll find up to four things:

1. the **spec text** being interpreted, quoted verbatim;
2. **ferro's reading** of that clause, in one or two sentences;
3. a table of **example spellings**, each with ferro's verdict and the form ferro normalizes it to;
4. a **Why** block linking the governing decision in the ruling ledger (see [the ledger](#the-ruling-ledger-and-ruling-ids) below).

## The verdicts

The recommendations are HGVS's, not ferro's. ferro's job is to take any input and produce the form
the recommendations prefer. The verdict describes **ferro's output** against that goal:

<span class="verdict-badge verdict-recommended">recommended</span> &nbsp; ferro's output is the
form the recommendations prefer — whether the input was already that form, or ferro normalized it
there. This is the target.

<span class="verdict-badge verdict-conformant">conformant</span> &nbsp; ferro's output is valid
HGVS, but not *yet* the recommended form — a known limitation, or a deliberate house choice among
forms the spec treats as equally legal. Where an issue tracks it, the row links to it.

<span class="verdict-badge verdict-refused">refused</span> &nbsp; the input is not valid HGVS;
ferro rejects it in strict mode. Refusing is the correct behavior here.

<span class="verdict-badge verdict-bug">bug</span> &nbsp; ferro's output is not valid HGVS — a
defect. This should never happen; where a page shows one, it links to the open issue. Rare by design.

## The table columns

| Column | What it means |
|---|---|
| **Input** | the HGVS string fed to ferro, verbatim |
| **Verdict** | one of the four above, describing ferro's output |
| **Normalizes to** | the form ferro produces — see the two special tokens below |
| **Notes** | a short plain-language explanation of that row |

In the **Normalizes to** column:

- <span class="norm-chip norm-self">= input</span> (written `self` in the source) — the input is a
  **fixed point**: ferro leaves it exactly as given, because it is already the recommended form.
- <span class="norm-chip norm-na">not run here</span> (written `—` in the source) — this spelling
  is **not run against the test reference on this page**, usually because it uses an accession
  outside the small committed slice ferro tests against. The row is still checked for *parsing*
  (and for strict-refusal where marked <span class="verdict-badge verdict-refused">refused</span>),
  just not for its normalized output.
- anything else is the literal string ferro rewrites the input to.

## Recurring terms

### parse-only here

The example uses a real accession that isn't in the committed reference slice, so ferro can *parse*
it but can't look up its bases to normalize it here. The row proves the spelling is well-formed HGVS;
the normalization is exercised elsewhere against a full reference.

### confluence

Several different input spellings that all **normalize to the same output**. A page notes confluence
when, for example, a split form and its merged delins both land on the one recommended string — the
point being that ferro maps a family of equivalent inputs onto a single canonical answer.

### strict vs. lenient

Two error-handling modes. **Strict** rejects anything not valid HGVS at parse time (this is what a
<span class="verdict-badge verdict-refused">refused</span> verdict checks). **Lenient** instead
repairs the input where it can and emits a warning, failing only if it cannot normalize at all. Some
rows describe behavior that only the lenient tier can express.

### warning codes (`W####`)

A machine-readable code ferro attaches to a normalization event — for instance a merge or a repair —
alongside the human-readable message. They identify the exact rule that fired; you don't need to
memorize them.

### The ruling ledger, and ruling ids

Where ferro has had to settle a question the spec leaves open, the decision is recorded **once** in
the ruling ledger, keyed by a short **ruling id** like `canonical-form-choice-when-both-legal`. A
**Why** block on a page transcludes that record's one-line summary and links to its full entry.
The reasoning lives in the ledger, never re-typed on the pages — a rule written in two places is a
rule that drifts. The ledger is rendered for reading as
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md).

## How it stays honest

Every example is executed through ferro by a build check: each spec quote is verified against the
pinned spec checkout, each spelling is parsed (and strict-refused where marked
<span class="verdict-badge verdict-refused">refused</span>), and each normalized output is asserted
against a real GRCh38 reference. A wrong verdict, a moved quote, or a changed normalization fails the
build — so a page cannot claim behavior ferro does not have.
