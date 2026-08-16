# Shadow Spec

> **Under construction.** Pages are added one HGVS spec page at a time; this section will fill in
> as they land.

## What it is

The HGVS recommendations are, in places, silent, ambiguous, or self-contradictory — but a
normalizer still has to emit one string. The **shadow spec** mirrors the recommendations page for
page and records, per clause, how ferro reads it.

Each page gives, for a span of the spec:

- the **spec text** it shadows, quoted verbatim;
- **ferro's reading** of that clause;
- **example spellings** with a verdict — `recommended` (the form ferro emits), `conformant`
  (legal, but ferro may rewrite it), `invalid` (refused in strict mode), or `gap` (a known
  ferro↔spec divergence, pinned to current behavior);
- the form each input **normalizes to**, and which spellings **converge** on one output;
- links to the governing **rulings** and README rules, so the *why* is one click away.

## What it is not

It is not a second copy of the rules. The canonical ruleset lives in the project README, and the
adjudication decisions live in the ruling ledger. The shadow spec **links** to both and never
restates them — a rule written in two places is a rule that drifts.

## How it stays honest

Every example is executed through ferro by a build check: each spec quote is verified against the
pinned spec checkout, each spelling is parsed (and strict-refused where marked `invalid`), and each
`normalizes to` is asserted against a real GRCh38 reference. A wrong verdict, a moved quote, or a
changed normalization fails the build — so a page cannot claim behavior ferro does not have.

## Pages

_None published yet._ The first page (DNA substitution) is in review.
