# Interpreting the HGVS Recommendations

The [HGVS nomenclature recommendations](https://hgvs-nomenclature.org/) are the standard for naming
sequence variants. They are large and carefully written, but — like any natural-language standard —
in places they leave a question open, or two passages can be read as pointing different ways. A
normalizer still has to emit exactly one string. Where ferro has had to settle such a question, this
section records the decision, mirrored page-for-page against the official recommendations.

<div class="use-cards">
  <div class="use-card">
    <h3>32 pages, one per spec section</h3>
    <p>Each page shows the spec's own text, ferro's reading of it, and a table of example spellings
    with ferro's verdict on each.</p>
  </div>
  <div class="use-card">
    <h3>Every example is executed</h3>
    <p>A build check runs each spelling through ferro against a real GRCh38 reference. A wrong
    verdict, a moved quote, or a changed normalization fails the build.</p>
  </div>
  <div class="use-card">
    <h3>New to the vocabulary?</h3>
    <p>Start with <a href="reading-guide.md">How to read a page</a> — the verdicts, the table
    conventions, and the recurring terms, explained once.</p>
  </div>
</div>

## The verdicts, at a glance

The recommendations are HGVS's, not ferro's. ferro's job is to take any input and produce the form
the recommendations prefer; the verdict describes **ferro's output** against that goal.

<span class="verdict-badge verdict-recommended">recommended</span> &nbsp; the form the
recommendations prefer — the target. &nbsp;&nbsp;
<span class="verdict-badge verdict-conformant">conformant</span> &nbsp; valid HGVS, but not *yet*
the recommended form. &nbsp;&nbsp;
<span class="verdict-badge verdict-refused">refused</span> &nbsp; not valid HGVS; ferro correctly
rejects it. &nbsp;&nbsp;
<span class="verdict-badge verdict-bug">bug</span> &nbsp; ferro produced invalid HGVS — a defect.

Full definitions, plus the table conventions (`self`, `—`, "parse-only") and the recurring terms,
are on the [How to read a page](reading-guide.md) reference.

## Pages

### DNA

| Page | Covers |
|---|---|
| [Substitution](recommendations/DNA/substitution.md) | one base replaced by one other base |
| [Deletion](recommendations/DNA/deletion.md) | removing one or more bases |
| [Duplication](recommendations/DNA/duplication.md) | a tandem copy of existing bases |
| [Insertion](recommendations/DNA/insertion.md) | new bases between two adjacent positions |
| [Inversion](recommendations/DNA/inversion.md) | a reverse-complemented span |
| [Deletion-Insertion](recommendations/DNA/delins.md) | one or more bases replaced by one or more others |
| [Repeated Sequences](recommendations/DNA/repeated.md) | short tandem repeat notation |
| [Alleles](recommendations/DNA/alleles.md) | combining variants in cis / trans / unknown phase |
| [Complex](recommendations/DNA/complex.md) | translocations and multi-part rearrangements |
| [Other](recommendations/DNA/other.md) | cases the rest of the list doesn't cover |

### RNA

| Page | Covers |
|---|---|
| [Substitution](recommendations/RNA/substitution.md) | one base for one base, on the transcript |
| [Deletion](recommendations/RNA/deletion.md) | removing one or more transcript bases |
| [Duplication](recommendations/RNA/duplication.md) | a tandem copy on the transcript |
| [Insertion](recommendations/RNA/insertion.md) | new transcript bases between two positions |
| [Inversion](recommendations/RNA/inversion.md) | a reverse-complemented transcript span |
| [Deletion-insertion](recommendations/RNA/delins.md) | transcript bases replaced by others |
| [Repeated Sequences](recommendations/RNA/repeated.md) | tandem repeat notation, transcript coordinates |
| [Alleles](recommendations/RNA/alleles.md) | cis / trans / unknown / product combinations |
| [Splicing](recommendations/RNA/splicing.md) | splice-site consequences written as del/ins |
| [Adjoined transcript](recommendations/RNA/adjoined_transcript.md) | gene-fusion transcripts |

### Protein

Protein examples run against `NP_003997.1` (dystrophin, the translation of `NM_004006.3`) in the
committed slice. There are no genomic sequence alignments on this axis — a `p.` description has no
SPDI — so the pages rely on the executed verdict and "normalizes to" columns.

| Page | Covers |
|---|---|
| [Substitution](recommendations/protein/substitution.md) | one amino acid for one (missense, nonsense, silent) |
| [Deletion](recommendations/protein/deletion.md) | removing one or more residues, 3′-placed |
| [Duplication](recommendations/protein/duplication.md) | a tandem copy of residues |
| [Insertion](recommendations/protein/insertion.md) | residues added between two flanking positions |
| [Deletion-insertion](recommendations/protein/delins.md) | residues replaced by others |
| [Repeated Sequences](recommendations/protein/repeated.md) | amino-acid repeat notation |
| [Frameshift](recommendations/protein/frameshift.md) | translation shifted to another frame |
| [Extension](recommendations/protein/extension.md) | the sequence extended past a terminus |
| [Alleles](recommendations/protein/alleles.md) | cis / trans / unknown phase on the protein |

### Cross-cutting

These three shadow the spec's non-axis pages — the rules and notation that apply across DNA, RNA,
and protein.

| Page | Covers |
|---|---|
| [General notation](recommendations/general.md) | the DNA-level-primary rule, mandatory prefixes, prioritisation, and the special characters (`=`, `/`, `_`, …) |
| [Uncertain descriptions](recommendations/uncertain.md) | `( )` / `?` / `^` notation for predicted or incompletely-known variants, across all axes |
| [Publication checklist](recommendations/checklist.md) | the most frequently offended rules, as correct-vs-refused pairs |

## What this section is *not*

It is not a second copy of the rules. The canonical ruleset lives in the
[normalization rules](../reference/normalization-rules.md) page, and every adjudication is recorded
once in the **ruling ledger**
([`hgvs_spec_normalization_overrides.json`](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/tests/fixtures/grammar/hgvs_spec_normalization_overrides.json),
rendered for reading as
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)).
These pages **draw their explanations from the ledger** rather than restating them — a rule written
in two places is a rule that drifts.

## How it stays honest

Every example is executed through ferro by a build check: each spec quote is verified against the
pinned spec checkout, each spelling is parsed (and strict-refused where marked
<span class="verdict-badge verdict-refused">refused</span>), and each normalized output is asserted
against a real GRCh38 reference. A wrong verdict, a moved quote, or a changed normalization fails the
build — so a page cannot claim behavior ferro does not have.
