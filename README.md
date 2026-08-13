[![CI](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml/badge.svg)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/ci.yml)
[![Nightly reference-aware tests](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml/badge.svg?branch=main)](https://github.com/fulcrumgenomics/ferro-hgvs/actions/workflows/nightly-mutalyzer.yml)
[![Codecov](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/ferro-hgvs)
[![Crates.io](https://img.shields.io/crates/v/ferro-hgvs.svg)](https://crates.io/crates/ferro-hgvs)
[![PyPI](https://img.shields.io/pypi/v/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Python versions](https://img.shields.io/pypi/pyversions/ferro-hgvs.svg)](https://pypi.org/project/ferro-hgvs/)
[![Documentation](https://docs.rs/ferro-hgvs/badge.svg)](https://docs.rs/ferro-hgvs)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ferro-hgvs/README.html)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.18703103-blue.svg)](https://doi.org/10.5281/zenodo.18703103)

# ferro-hgvs

A high-performance HGVS variant nomenclature parser and normalizer written in Rust.

**WARNING: ALPHA SOFTWARE - USE AT YOUR OWN RISK**

This software is currently in **ALPHA**. While we have extensively tested it
across a wide variety of HGVS patterns, **no guarantees are made** regarding
correctness or stability.

<p>
<a href="https://fulcrumgenomics.com">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-light.svg">
    <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/ferro-hgvs/main/.github/logos/fulcrumgenomics-light.svg" height="100">
  </picture>
</a>
</p>

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## Features

- **Full HGVS Parsing**: All coordinate systems (g/c/n/r/p/m/o) and edit types
- **Variant Normalization**: 3'/5' shifting per HGVS specification
- **High Performance**: ~5M variants/sec single-threaded parsing (>12M/s parallel), zero-copy with nom
- **Type-Safe**: Leverages Rust's type system for correctness

## Installation

### Python

```bash
pip install ferro-hgvs
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, Apple Silicon), and Windows (x86_64) on Python 3.10+.

### Rust

Add to your `Cargo.toml`:

```toml
[dependencies]
ferro-hgvs = "0.1"
```

Or install the CLI:

```bash
cargo install ferro-hgvs
```

## Quick Start

### CLI

```bash
# Parse a variant
ferro parse "NM_000088.3:c.459A>G"

# Parse from file
ferro parse -i variants.txt -f json

# Prepare reference data (downloads RefSeq, genome, cdot — RefSeq-only by default)
ferro prepare --output-dir ferro-reference

# Verify reference data is ready
ferro check --reference ferro-reference

# (Optional) pre-build the on-disk cdot cache as a setup step, so the one-time
# cache build doesn't slow the start of a real (or timed/benchmarked) run.
ferro check --reference ferro-reference --build-cache

# Normalize with reference
ferro normalize "NM_000088.3:c.459del" --reference ferro-reference/
```

> **Read the warnings.** Normalization sometimes *repairs* a description in a way
> the normalized string does not record — separately reported cis members merged
> into one `delins` (`MEMBERS_COALESCED_FROM_REPORTED_FORM`), a `ins[100_110]`
> reference-range payload replaced by the bases it denotes
> (`INSERTED_SEQUENCE_EXPANDED`), a stated reference base that contradicted the
> reference and was accepted anyway (`REFSEQ_MISMATCH`). Those are reported as
> `warning[CODE]: message` on **stderr** (and in the `warnings` array under
> `--format json`, the `detail` column under `--format tsv`), so a pipeline
> reading only stdout will not see them. `--error-mode strict` is not a
> substitute: it rejects a specific ladder of conditions and reports the rest
> exactly as lenient does.

> **Throughput tip:** when normalizing many variants, feed them **sorted by transcript accession (or by genomic position)**. ferro caches each resolved transcript, so consecutive variants on the same transcript skip the (dominant) cost of re-reading and re-building it from the reference. Sorted input keeps the relevant transcripts resident in the cache and is markedly faster on large batches — see [Performance Comparison](#performance-comparison).

### Optional reference data

A bare `ferro prepare` builds a **RefSeq-only** reference (accessions `NM_`/`NR_`/`NP_`/`NG_`). Two opt-in flags provision additional data — pass them at prepare time; they are what a fully-provisioned ("blessed") reference is built with:

```bash
# Add Ensembl support (accessions ENST/ENSG/ENSP). Downloads the Ensembl cdot
# metadata and cDNA FASTAs (~1 GB+); off by default. Without it, an ENST/ENSG/ENSP
# input reports "Reference not found" and the message points back at this flag.
ferro prepare --output-dir ferro-reference --ensembl

# Derive version-independent NG_ placements and the NG_→transcript-version map
# (ng_hosted_transcripts) for a curated list of RefSeqGene accessions. Required to
# resolve legacy gene-symbol selectors (NG_(GENE):c.…) and bare-NG_ hosted lookups.
ferro prepare --output-dir ferro-reference \
  --derive-ng-placements path/to/ng_accessions.txt

# A fully-provisioned reference combines both in one run:
ferro prepare --output-dir ferro-reference --ensembl \
  --derive-ng-placements path/to/ng_accessions.txt
```

Both flags are incremental: re-running `ferro prepare` over an existing reference adds the requested data and preserves already-provisioned artifacts.

### Library

```rust
use ferro_hgvs::{parse_hgvs, HgvsVariant};

fn main() -> Result<(), ferro_hgvs::FerroError> {
    let variant = parse_hgvs("NM_000088.3:c.459A>G")?;

    match &variant {
        HgvsVariant::Cds(v) => println!("CDS variant: {}", v),
        HgvsVariant::Genome(v) => println!("Genomic variant: {}", v),
        _ => println!("Other: {}", variant),
    }

    Ok(())
}
```

### Python

```python
import ferro_hgvs

# Parse a variant
variant = ferro_hgvs.parse("NM_000088.3:c.459A>G")
print(variant.variant_type)  # "coding"
print(variant.reference)     # "NM_000088.3"
print(str(variant))          # "NM_000088.3:c.459A>G"

# Normalize with reference data
normalizer = ferro_hgvs.Normalizer(reference_json="ferro-reference/cdot.json")
normalized = normalizer.normalize("NM_000088.3:c.459del")

# `normalize` returns only the string, so it cannot tell you that normalization
# repaired something. `normalize_with_warnings` returns the same string plus the
# diagnostics — as a free function, or as a Normalizer method.
result = normalizer.normalize_with_warnings("NM_000088.3:c.459del")
print(str(result.result))                      # the same normalized string
print([(w.code, w.message) for w in result.warnings])
```

## Normalization rules

ferro's normalizer follows four rules about its output, and three about how it handles the gaps.

### The output contract

1. **Conformant.** Output follows the HGVS recommendations. *Absolute — never traded.*
   Scope: the `syntax.yaml` grammar, plus every prohibition, read from **prose force** rather than
   keyword casing — "not allowed", "not correct", "can not be used", "by definition", the
   `class="invalid"` markup, and the set `checklist.md` enumerates.

2. **Recommended form.** Where the spec prefers among conformant forms, ferro produces it.
   *Best effort.*
   Scope: "recommended", "preferred", lowercase "should", the 3' rule. A preference clause outranks
   maintainer judgment, but not rule 3: where it is not evaluable on what a normalizer holds, rule 3
   governs.

3. **Confluent.** Inputs denoting one variant produce one output. *Best effort.*
   Every rule — rule 2 included — is evaluated over the **resulting sequence**, never over the
   input's spelling. Reference context counts as sequence: transcript model, exon boundaries,
   reading frame, strand and topology are functions of the accession, not of the description.

4. **Deterministic.** Same input, same output. *Absolute.*
   Note that 4 does not imply 3 — a deterministic normalizer can be arbitrarily non-confluent.

### The procedure

5. **Where the spec is silent, ambiguous, or self-contradictory:** the issue is filed upstream
   **first** and cited, and only then does ferro ship a provisional choice.
    - **Self-contradictory** — two clauses that cannot both hold — is a **defect**. Every conforming
      tool must pick a side and none of them can be right, so filing is a bug report and is not
      optional.
    - **Silent** is merely incomplete. Ferro decides under rule 6 and violates nothing; filing is a
      feature request, worth making but never a reason to hold a release.

6. **Among multiple conformant forms:** the maintainers choose. **There are no user options for
   normalization form.** Error mode is an orthogonal axis and stays available. The 3'/5' knob is
   **not** orthogonal — it selects the frame every rule is evaluated in, so rules 2 and 3 are
   claimed *per direction* rather than across the two.

7. **Disclosure.** Any change to these rules, and any different choice made under 5 or 6, is
   disclosed: in the changelog before v1, by a major version bump after. Output that *violates*
   rules 1-4 is a **bug**, not a disclosure.

### Why 2 and 3 are best effort, and 1 and 4 are not

Rules 1 and 4 are always achievable. Conformance is checkable against the spec text, and
determinism is a property of ferro's own code; nothing external can prevent us honouring them.

Rules 2 and 3 depend on the spec **determining an answer**, and sometimes it does not:

- **No preference exists.** The spec ranks substitution, deletion, inversion, duplication and
  insertion, but says nothing about competing `delins` forms.
- **Two preferences disagree.** `general.md` ranks duplication above insertion; `DNA/inversion.md`
  prefers an insertion for inverted copies. No single output satisfies both.
- **The same clause has two versions.** `general.md`'s current text and its forthcoming NOTE give
  opposite answers for variants separated by one nucleotide.
- **The variant's decomposition is not recoverable.** For an equal-length block the column
  correspondence is unique, so "two variants separated by N nucleotides" is a fact about the
  variant. For an unequal-length block there is no correspondence — recovering one means *choosing*
  an alignment, and the spec does not say which. There is no derivable form to converge on.
- **The preference keys on information a normalizer does not hold.** A "frequently occurring
  variant" (`RNA/delins.md:41`); a repeat "variable in the population" (`RNA/repeated.md:33`,
  `protein/repeated.md:22`); two variants "reported (or might occur) individually" (`DNA/delins.md:83`).
  The spec determines an answer; ferro cannot see what it keys on.

**"Best effort" is bounded by the spec's determinacy and by what a normalizer's inputs can decide —
not by ferro's implementation quality.** A failure of rule 2 or 3 caused by ferro's code is a
**bug** under rule 7; one caused by the spec not determining an answer triggers rule 5; one caused
by a clause keying on provenance or population data is a **declared deviation**, and a permanent
one — no upstream answer can put that information into a normalizer's hands.

Permanent is the only thing the last case adds, and it does **not** exempt the question from rule
5. The grain matters, because the two grains give opposite answers. Read **on its own**, such a
clause is not silent, ambiguous or self-contradictory — it determines an answer perfectly well, and
ferro is the one that cannot see what the answer keys on. Read **against the spec's own
re-derivation mandate** — `general.md:157-160`, which has a protein description derived by
comparing the variant and reference protein sequences and says knowledge of the underlying DNA
change "should not be used", with `general.md:13` extending the method to RNA — it is one half of a
pair that cannot both hold. That is rule 5's self-contradiction limb exactly, and the ruling ledger
says so in those words. What the carve-out states is that rule 5's escalation cannot *end*
the matter here: a choice made under rule 5 is provisional pending an upstream ruling, and no
ruling upstream makes provenance visible, so the deviation outlives the filing.

Where it is recorded is worth stating precisely, because the obvious pointer does not resolve.
`canonical-form-choice-when-both-legal` has **no `deviates_from` field**: it carries all four
clauses in its rationale as the recorded counter-evidence to re-derivation, opens that paragraph
with "The spec contradicts itself here", and adopts re-derivation over them deliberately. The one
of the four that is recorded *as* a deviation is `DNA/delins.md:83`, through the
`deviates_from: ["docs/recommendations/DNA/delins.md:79-84"]` on
`separation-is-a-property-of-the-spelling-not-of-the-variant` — the record whose ruling the rule 3
example below turns on.

### A worked example of reading force from prose

`DNA/duplication.md` says a variant that *can* be described as a duplication **must** be — but the
"must" is scoped by the preceding clause, which defines when a duplication *can* be used at all:
only when the additional copy is directly 3'-flanking the original. So the rule ranks the *label*
for one span; it does not require that a partition be chosen so as to produce a duplication.
Reading the force without the scope inverts the rule.

### What rule 3 excludes

"Never over the input's spelling" is narrower than it sounds. Most of a description is context.

| Carried by the description | Treatment |
|---|---|
| Accession, axis (`g.`/`c.`/`n.`), version | **Used** — it *is* the reference context |
| Which bases end up different | **Used** — this is the variant |
| Cis against trans (`[a;b]` against `[a];[b]`) | **Used** — different variants, not two spellings of one |
| Type label (`dup` against `ins`, `inv` against `delins`) | **Re-derived**, then ranked by `general.md` |
| How the edit set is cut into members | **Excluded** — a property of who wrote the string |
| Which copy in a run of identical residues a member names | **Excluded** — the 3' rule assigns this *"arbitrarily"* (`general.md:41`) |
| Repeat unit and phase, where several are equivalent | **Excluded** |

So three rows read **Excluded**, and they are three spellings of the one thing the input does not
get to decide: the **partition**, the run-position choice that feeds it, and — where several
unit-and-phase pairs describe one tract equally well — which pair a repeat member names.

`NC_000001.11:g.1001002_1001016` reads `ATGAGGGGCCACTGT`: a `GGGG` run at `1001006-1001009`, a `CC`
run at `1001010-1001011`, a lone `C` at `1001013`. Two spellings, one denoted sequence
(`ATGAGGGCATGT`), because `1001010` and `1001011` are both `C`:

```
g.[1001009del;1001010del;1001013del]     written gaps of 0 and 2
g.[1001009del;1001011del;1001013del]     written gaps of 1 and 1
```

`general.md:34`'s "two variants separated by one or more nucleotides should be described individually"
reads those gaps, so it answers twice for one variant. Rule 3 reads them off the partition ferro
derives instead: both give `g.[1001009_1001010del;1001013del]`. Rule 2 then keeps that over the
spanning `g.1001009_1001013delinsCA`, which merges across two unchanged nucleotides. Pinned in
[`tests/it/cis_confluence_adjudication.rs`](tests/it/cis_confluence_adjudication.rs).

### Known limitation

ferro cannot today guarantee that every input form is normalized according to these rules. They are
enforced intent, not a claim of current completeness.

Rule 7's disclosure mechanism — the `Representation-Change:` trailer and how it reaches the
changelog — is documented in
[CONTRIBUTING.md](CONTRIBUTING.md#declaring-a-representation-change).

### Where the individual decisions live

Rules 5 and 6 above say *how* a question is decided where the recommendations are silent,
ambiguous or self-contradictory. **What was decided, case by case**, is recorded separately, as
adjudication records — each naming the clauses in tension, which one governs, which is deviated
from, and why. Those records are published in full, with their clause quotes, in
[docs/NORMALIZATION_CONTRACT.md](docs/NORMALIZATION_CONTRACT.md).

That document is generated from the records and gated against them, so it cannot drift; it
deliberately does not restate the seven rules above, which are stated only here.

## Deriving a description from sequences

If what you have is **bases** rather than a description — a window out of a BAM, a VCF row, an
aligner's output — `from_sequences` derives the description instead of asking you to spell one:

```python
import ferro_hgvs

ferro_hgvs.from_sequences("NC_000001.11", 1000, "AGCGT", "AGT")
# NC_000001.11:g.1002_1003del
```

```rust
use ferro_hgvs::{from_sequences, FromSequencesOptions};

let variant = from_sequences("NC_000001.11", 1000, "AGCGT", "AGT",
                             &FromSequencesOptions::default())?;
```

The axis follows the accession: `g.`, or `m.` on `NC_012920` / `NC_001807`, which HGVS requires the
`m.` coordinate system for. Every other accession class — transcript, protein, UniProt — is refused
with a message naming it.

**It reads no reference sequence.** The output is a pure function of its arguments — the
accession, the position, the two sequences and the options — so the same bases give the same
description on any machine, against any reference build, with no hidden input. `position` is
1-based, and `reference` is taken on trust: verifying it would need the reference and would make
the provider a hidden input, costing exactly the determinism the function exists to provide.

That is five values, not four. The options are not inert — `direction` moves a placement within
the window and `max_grid_cells` decides whether an answer is produced at all — so the
"four arguments" this section and the Rust docs both used to claim is withdrawn. Purity is the
property; the count was wrong.

### How to: one canonical description per variant

The job this exists for — a pipeline that aligns reads, post-processes a BAM, and wants one
description per variant, decided by the observed bases and nothing else.

**1. Get a window.** From a pileup, a VCF row, or an aligner's output: the reference bases over
some interval, the observed bases over the same interval, and the 1-based position of the window's
first base. If what you hold is a *description* rather than bases, `Normalizer::to_sequences`
produces the same window — see [Going the other way](#going-the-other-way).

**2. Derive, and read the flag.**

```python
d = ferro_hgvs.from_sequences_detailed("NC_000001.11", 1000, "AGCGT", "AGT")
d.variant                      # NC_000001.11:g.1002_1003del
d.placement_bounded_by_window  # False — the bases settled the placement, not the window

# The same variant, read through a window that stops at it:
d = ferro_hgvs.from_sequences_detailed("NC_000001.11", 1000, "AGCG", "AG")
d.variant                      # NC_000001.11:g.1002_1003del — the same answer
d.placement_bounded_by_window  # True — the deletion is flush with the window's 3' edge
```

The two rows are the flag's whole meaning: **same description, different confidence that a wider
read would agree.** The first window has a base to spare 3' of the deletion, so the placement is
settled by the bases; the second ends exactly where the deletion does, so the flag fires even
though the answer is unchanged. This example used to show only the second window while claiming
`False`, which inverted both the value and its explanation.

Prefer `from_sequences_detailed` to `from_sequences` in a pipeline. The flag is the only thing that
tells you the window may have decided the answer, and a bare `from_sequences` discards it.

**3. Store it, or normalize first.** What you have is already conformant (rule 1) and deterministic
(rule 4), so it is safe to store and to compare between runs and between machines. It is *not*
necessarily the recommended form, and it is not guaranteed to agree with a description derived from
a different window — those are rules 2 and 3, and both need the reference. If you want them, run
`normalize` on the result. That is the whole offer: **derive now, normalize later, or never.**

**4. If the flag is set and you need the reference-anchored answer**, either re-derive from a wider
window or run `normalize`. Do not treat a flagged result as wrong — see below.

#### How wide should the window be?

Wide enough to contain the whole interval over which the change could legally be placed. That is
the exact condition, and it is what the [cost section](#the-cost-stated-plainly) pins.

Operationally: **pad on both sides by at least the length of the longest ambiguous run — a
homopolymer or tandem repeat — the variant might sit in.** An insertion or duplication needs one
further base 5', since that is where a 5'-most insertion anchors. `to_sequences` pads by 128 on
**each** side, so its window is `span + 2 * pad`; that covers ordinary repeats comfortably, and is
worth raising if you work with long tracts.

Two things follow that are easy to get backwards:

- **A window that cuts the interval does not give a wrong answer, it gives a bounded one.** The
  description still denotes the same bases and carries the same canonical SPDI; it is simply placed
  at the window's edge rather than the run's. What you lose is the *preferred spelling*, not
  correctness.
- **`placement_bounded_by_window` is conservative on purpose.** It reports "this could have moved",
  not "this is wrong" — a window flush with a tract is flagged and is nonetheless the same answer a
  whole-sequence derivation gives. Distinguishing the two needs the reference, which this function
  does not read.

#### The refusals, and what to do about each

The policy is to refuse rather than quietly answer with a weaker rule, so each of these is
actionable rather than fatal:

| refusal | why | what to do |
|---|---|---|
| the accession is not genomic | this surface emits `g.` (and `m.` on the two rCRS mitochondrial accessions) and nothing else; `NM_…:g.9_10del` would be well-formed and denote nothing | pass the genomic accession, and project afterwards if you need a transcript axis |
| an inserted payload sits against the window's 5' edge | HGVS writes an insertion *between* two positions, so it would have to anchor at `position - 1` — outside the window | re-fetch with more 5' flank. How much is direction-dependent, so widen rather than adding a fixed one base |
| the alignment grid exceeds `max_grid_cells` | a cost bound — a cell is roughly 18 bytes, and the default admits a window of about 4 096 bases | raise `max_grid_cells` if you have the memory, or narrow the window. Real structural alleles are far past any sane budget: `LRG_542:g.[101177_102434delins36;107248_127198delins21]` spans 26 kb |
| a symbol outside the IUPAC-IUBMB set, a zero `position`, an empty `reference` | the input cannot denote anything | fix the input. `X` and `-` are refused deliberately (`standards.md:39`) — they are alignment symbols, not bases |
| `U` in either sequence | this surface's axis is DNA; a `g.`/`m.` description naming `U` would be well-formed and wrong | pass `T`, and project onto an `r.` axis afterwards if you need RNA |

Case is **not** a refusal: a soft-masked (lower-case) window derives exactly as its upper-case twin
does, and both sequences are folded before anything reads them.

### Which rules it delivers

This is the whole design, and it falls straight out of the four rules above:

| | rules delivered | force | needs |
|---|---|---|---|
| `from_sequences` | **1 (conformant)**, **4 (deterministic)** | both *absolute* | the caller's four arguments |
| `normalize`, afterwards | **2 (recommended form)**, **3 (confluent)** | both *best effort* | the reference |

Rules 1 and 4 are the two the section above calls always achievable, so a function that has only
the caller's arguments can still deliver both in full. Rules 2 and 3 need the reference: rule 2's
scope names the 3' rule explicitly, and a reference-anchored shift is precisely what a
window-local function cannot perform.

So an output may be 3'-shiftable further than the window allowed, and **that is not a defect**.
Run `normalize` afterwards if you want it — `Normalizer::from_sequences(..., normalize = true)`
does both in one call.

**Run it unless you have a reason not to.** Over a 6,000-shape sweep, `normalize` moved **8.6%** of
derived descriptions — in three classes: repeat notation (`g.27_28insAAA` → `g.27A[4]`),
reference-anchored member re-derivation, and an inversion spread across several members
(`g.[17C>A;19T>A;21T>G]` → `g.17_21inv`, which the alignment DAG partitions before anything can see
it, since it minimises edit distance and an inversion is not in that cost model). All three are
rule 2 and rule 3 — the recommended form and agreement with a wider view — which is exactly the
pair this design assigns to `normalize`. Rules 1 and 4 hold either way.

That figure is a claim about that sweep, not about the world: one synthetic contig, six shape
generators, genomic axis only.

### What you get for it

Two spellings of one variant, **over one window**, reach one description — because the derivation
never sees a spelling. Over the cis confluence corpus that is 5 636 classes with no divergence —
**its genomic half, and all of what this surface can reach**: the corpus is generated `--axes g,c`
at 11 272 classes, and the 5 636 `c.` classes are drawn against `NM_TEST.1`, which the `g.`-only
gate refuses, so not one of them enters the comparison. The exclusion is structural and is asserted
as such in `tests/it/from_sequences_corpus.rs`.

Over the nine externally-reported confluence pairs (#1419 / #1420 / #1421) it is nine of nine, in
both shuffle directions — where `normalize`, handed the same pairs as descriptions, currently
converges none of them.

**Read "over one window" as load-bearing, not as hedging.** It is what makes the claim arithmetic:
`to_sequences` computes its window from the *denoted bases*, so both spellings of a variant get
byte-identical `(position, reference, alternate)` triples, and a pure function of that triple can
only give one answer. It is also the exact limit — the claim is confluence over **spellings**, and
rule 3's scope is confluence over **inputs**. Two reads covering one variant differently are two
inputs, and the section below is what happens then.

Read the comparison with `normalize` as two functions answering different questions rather than as
one beating the other. A caller who *has* a description and wants it normalized still needs
`normalize` to converge; the pairs are simply the case where being handed a description is itself
the problem.

### The cost, stated plainly

A window-local derivation is **read-dependent**. Nothing may shift outside the bases you supplied,
so a read that stops partway through an ambiguous run places the change at the end of the *read*
rather than the end of the run. One deletion from the `AAAA` at 12–15 of a test contig, seen
through three windows:

| window | reference | alternate | derived | `placement_bounded_by_window` |
|---|---|---|---|---|
| 10–16 | `GCAAAAG` | `GCAAAG` | `g.15del` | `false` |
| 12–15 | `AAAA` | `AAA` | `g.15del` | `true` |
| 10–14 | `GCAAA` | `GCAA` | **`g.14del`** | `true` |

**None of these is wrong.** All three carry the same canonical SPDI (`14:A:`) and denote the same
bases — `g.14del` is a conformant description of exactly the same variant, it is simply not the
3'-most spelling. So what a truncating read costs you is the **recommended form** (rule 2) and
**agreement with a wider read** (rule 3): precisely the two rules this function never claimed,
because both need the reference. Rules 1 and 4 hold in every row. `normalize` closes the gap,
shifting to 15 regardless of what the read covered.

The boundary is exact, and pinned in `tests/it/from_sequences_window_condition.rs`:

> **Two windows that both contain the whole interval over which the change can be placed derive
> the same description. A window that cuts that interval places the change at its own edge
> instead.**

Note that for an insertion or duplication that interval already reaches one base 5' of the tract,
since that is where a 5'-most insertion anchors — so "contains the interval" subsumes the flank
requirement rather than needing a separate clause.

**`placement_bounded_by_window` is a "could move" flag, not a "is wrong" flag**, and it is
conservative in that direction on purpose. Row 2 is flagged and already correct; row 3 is flagged
and merely non-preferred. Telling those apart requires knowing what lies outside the window — the
reference — which this function does not read, so it reports the uncertainty rather than resolving
it. Treat a `true` as "re-derive from a wider window, or run `normalize`, if you need a
reference-anchored answer".

Two refusals are worth knowing about in advance, both deliberate — the policy is to refuse rather
than degrade to a weaker rule:

- **An alignment grid over budget.** The default admits a window of about 4 096 bases; a cell costs
  roughly 18 bytes. `max_grid_cells` is the knob, and the refusal names it. Real *structural*
  alleles are well past it — `LRG_542:g.[101177_102434delins36;107248_127198delins21]` spans 26 kb.

  **Do not read the multi-member census as a measure of this refusal.** Of the 592 multi-member
  alleles harvested from ClinVar, CMRG and Paraphase, 443 windows were captured and 59 derive —
  but splitting the 384 refusals by message gives **384 accession refusals and 0 grid refusals**,
  every one of them an `NM_` transcript hitting the `g.`-only gate above. The structural rows are
  filtered out of the capture before the grid is ever consulted. That figure was quoted here as
  evidence for the grid bound and is not.
- **An inserted payload against the window's 5' edge.** HGVS writes an insertion between two
  positions, so such a payload can only be anchored at `position - 1` — outside the window, and
  non-existent when `position` is 1. Supply more 5' flank; `Normalizer::to_sequences` pads both
  sides for you.

### Going the other way

`Normalizer::to_sequences` is the inverse, so a caller who already holds descriptions needs no new
plumbing to reach the derivation:

```python
pair = normalizer.to_sequences(variant, pad=128)
derived = normalizer.from_sequences(pair.accession, pair.position, pair.reference, pair.alternate)
```

The pad is not decoration: `dup` typing reads the reference bases immediately 5' of an insertion
point (`DNA/duplication.md:18`), so a member flush with the window's 5' edge comes back as an `ins`
instead of a `dup`. It is applied to **both** sides — the window is `span + 2 * pad` — and the
bases come back upper-cased, so a soft-masked region does not produce a mixed-case pair.

### Bounding a derivation to a region it must not leave

When a variant must stay inside a target region, an amplicon or a tiling window, anchor every raw
pair to that region first. The derivation is a pure function of the window it is handed, so one
window gives one answer:

```python
pair = ferro_hgvs.SequencePair("chr1", 10, "GCAAAAG", "GCAAAG")   # straight from a BAM

str(pair.derive().variant)                  # 'chr1:g.15del' — rolls to the run's end
bounded = pair.trim_to(end=14)              # hold it at 14
str(bounded.derive().variant)               # 'chr1:g.14del'
```

`trim_to` needs no reference and can only *narrow*. To widen, use `Normalizer::reanchor`, which
reads the padding bases from the reference:

```python
anchored = normalizer.reanchor(pair, start=5, end=25)   # 5' widened, 3' widened
str(anchored.derive().variant)                          # 'chr1:g.15del'

both = normalizer.reanchor(pair, start=5, end=14)       # 5' widened, 3' narrowed
str(both.derive().variant)                              # 'chr1:g.14del'
```

**`reanchor` moves a window's edges; it does not relocate the window.** Each edge may go outwards
(padded from the reference) or inwards (trimmed), in any combination — but the window you ask for
must **overlap the pair's own**, and the overlap must still hold the bases the two sequences
disagree on. `reanchor(pair, start=1000, end=1200)` on the pair above is refused, not fetched: the
changed bases exist only in the pair, so there is nothing to carry to a region the pair does not
cover. So "anchor every raw pair to my target region" works exactly when every raw pair overlaps
that region — which is the case the feature is for, and is worth checking rather than assuming.

Prefer `pair.derive()` to re-spreading the four fields: a pair returned by `trim_to` or
`reanchor` carries its **own** position, and pairing a pre-trim position with post-trim bases is
the mistake the method exists to prevent.

Both take 1-based inclusive bounds, and `None` leaves that edge where it is.

**Both refuse rather than clamp**, in every case: a bound that would cut a base the two sequences
disagree on (naming the coordinate), a bound that would empty the reference, `start` past `end`,
and — for `reanchor` — a bound outside the sequence, or a window disjoint from the pair's. A window
silently pulled back to the contig would hide a bug upstream of the call.

Case is not a disagreement: a soft-masked reference against an upper-case alternate trims normally.
`trim_to` fetches nothing and so leaves your bases as you passed them; `reanchor` reads flank from
the provider and therefore returns the whole window **upper-cased**, exactly as `to_sequences`
does, rather than splicing provider bases onto caller bases and handing back a mixed-case pair.

> **Reach for this when the bound is a requirement, not to make heterogeneous inputs agree.** For
> that, `Normalizer::from_sequences(..., normalize = true)` and a `to_sequences` round trip both
> already converge, and both reach the *reference-anchored* placement — which can shift as far as
> the sequence allows rather than as far as your window allows. Anchoring to a window that cuts an
> ambiguous run makes every caller using that window agree with each other and disagree with the
> reference. That is a legitimate contract and a poor default; `placement_bounded_by_window`
> reports it either way.

## Supported HGVS Syntax

| Type | Prefix | Example |
|------|--------|---------|
| Genomic | `g.` | `NC_000001.11:g.12345A>G` |
| Coding DNA | `c.` | `NM_000088.3:c.459A>G` |
| Non-coding | `n.` | `NR_000001.1:n.100A>G` |
| RNA | `r.` | `NM_000088.3:r.459a>g` |
| Protein | `p.` | `NP_000079.2:p.Val600Glu` |
| Mitochondrial | `m.` | `NC_012920.1:m.3243A>G` |

### Edit Types

- Substitution: `A>G`, `Val600Glu`
- Deletion: `del`, `100_200del`
- Insertion: `100_101insATG`
- Deletion-Insertion: `100_102delinsATG`
- Duplication: `100_102dup`
- Inversion: `100_200inv`
- Repeat: `100CAG[20]`

## CLI Commands

The `ferro` CLI provides commands beyond parsing and normalization:

| Command | Description |
|---------|-------------|
| `prepare` | Download and prepare reference data for normalization |
| `check` | Verify reference data setup |
| `parse` | Parse and validate HGVS variants |
| `normalize` | Normalize HGVS variants (3'/5' shifting) |
| `explain` | Explain error/warning codes (e.g., `ferro explain W1001`) |
| `annotate-vcf` | Annotate VCF files with HGVS notation |
| `vcf-to-hgvs` | Convert VCF records to HGVS |
| `hgvs-to-vcf` | Convert HGVS to VCF format |
| `liftover` | Liftover coordinates between genome builds |
| `describe` | Generate HGVS from reference/observed sequences |
| `effect` | Predict protein effect from variant |
| `backtranslate` | Reverse translate protein to DNA variants |
| `convert-gff` | Convert GFF3/GTF to transcripts.json |
| `generate` | Generate HGVS descriptions from components |
| `extract-hgvs` | Extract HGVS from VEP-annotated VCFs |

## Error Handling

ferro-hgvs provides configurable error handling with three modes:

| Mode | Behavior |
|------|----------|
| `strict` | Reject non-conformant input (default) |
| `lenient` | Auto-correct with warnings |
| `silent` | Auto-correct silently |

```bash
# Use lenient mode to auto-correct common issues
ferro parse --error-mode lenient "p.val600glu"  # Corrects to p.Val600Glu

# Ignore specific warnings
ferro parse --ignore W1001,W2001 "p.val600glu"

# Get help on any error/warning code
ferro explain W1001
ferro explain --list
```

### Configuration File

Create `.ferro.toml` in your project directory:

```toml
[error-handling]
mode = "lenient"
ignore = ["W1001", "W2001"]  # Silently correct these
reject = ["W3003"]           # Always reject these
```

## Comparing normalization rules (`FERRO_PARTITION`)

> **Unstable.** `FERRO_PARTITION` is an evaluation switch, **not** a supported feature. It is not covered by semantic versioning, its values may change, and it is expected to be removed once the normalization rule is settled. Do not depend on it in production pipelines.

Normalization cuts each changed region of sequence into allele members. `FERRO_PARTITION` selects which rule does that cutting, so a candidate rule can be measured against the shipped one over a real corpus before anything changes for users.

| Value | Rule |
|-------|------|
| unset / empty / `live` | The shipped rule. This is what every normal invocation uses. |
| `shadow` | Cut only at alignment steps common to *every* minimal alignment. |
| `canonical` | The member-count-minimal minimal alignment. |
| `canonical-coalesced` | `canonical`, plus the `delins.md:44-47` merge: a split whose payload realigns as one block is re-spelled as a single `delins`. Applied after the downstream passes rather than at partition time, so it cuts identically to `canonical` and differs only in what survives. |

With the variable unset — or set to the empty string — output is byte-identical to a build with no switch at all.

### Running an A/B comparison

Run the same input twice and diff the normalized descriptions. In `tsv` format the columns are `line, input, normalized, changed, status, detail`, so column 3 is the normalized string:

```bash
# 1. the shipped behaviour
ferro normalize --input variants.txt --reference /path/to/reference \
  --format tsv --error-mode lenient -j 10 > shipped.tsv

# 2. the candidate rule
FERRO_PARTITION=canonical ferro normalize --input variants.txt --reference /path/to/reference \
  --format tsv --error-mode lenient -j 10 > candidate.tsv

# 3. what moved
diff <(cut -f3 shipped.tsv) <(cut -f3 candidate.tsv) | head

# 4. how many moved, counting only rows that succeeded on both sides
#    (columns 3/5 are `normalized`/`status`; a row that failed carries no
#     normalized string, so comparing column 3 alone would score two different
#     failures as identical)
paste <(cut -f3,5 shipped.tsv) <(cut -f3,5 candidate.tsv) \
  | awk -F'\t' '$2 == "ok" && $4 == "ok" && $1 != $3' | wc -l

# and, separately, rows whose status itself changed
paste <(cut -f5 shipped.tsv) <(cut -f5 candidate.tsv) | awk -F'\t' '$1 != $2' | wc -l
```

This works on a **stock release build** — the switch is not behind a build feature, so no special binary or wheel is needed.

### Two traps worth knowing

**A misspelled value is refused, loudly.** `FERRO_PARTITION=canonicl` makes `ferro` exit with an error before it reads any input, naming the value you gave and the arms this build has:

```
FERRO_PARTITION="canonicl" is not a partitioner this build has. This build's arms
are: live, shadow, canonical, canonical-coalesced. Refusing rather than falling
back to `live`, because a bake-off served the shipped rule under a candidate's
name reports that the candidate changes nothing.
```

Naming *this build's* arms is the point: a value that exists on some other branch, or that used to exist, is reported as absent here rather than quietly answered as `live`.

The refusal comes from the CLI, not from the normalizer. Up to and including v0.14.0 it was a panic raised deep inside normalization, which meant a development-only switch could abort any process — a long-running service, or a Python caller, across the FFI boundary — that merely happened to have the variable set. A release build of the *library* now falls safe to `live` for a value that names no arm and keeps the refusal for its caller to report. Every binary in this repository reports it — `ferro`, `ferro-web`, `ferro-benchmark`, both spec generators, and the `dump_normalized_corpus` bake-off harness — and that is pinned by a test rather than by this sentence, so a binary added later cannot quietly skip it. If you embed ferro as a library and run bake-offs through it, read `ferro_hgvs::normalize::partition_switch_startup_error()` at startup and do the same.

Builds up to and including v0.13.1 fell back to `live` instead, and produced a clean, empty diff that read as "the candidate changes nothing". The fallback emitted a warning through the `log` facade, but the `ferro` CLI installs no logger, so that warning reached no stream and `RUST_LOG` could not surface it — there was no signal at all. If you are on an older build, treat every empty diff as unproven.

A positive control on an input **known** to differ is still worth running, because it catches the other way a comparison can be vacuous: a variable that never reached the process at all (a lost `export`, a `sudo` that scrubbed the environment, a container that did not forward it). That case is indistinguishable from `unset`, which is legitimately `live`, so no amount of validation inside ferro can catch it. On your own corpus a zero remains ambiguous between "the switch is not taking effect" and "this corpus has no affected variants".

These three inputs run against the built-in test data, so they need no `--reference` and no prepared reference directory. Between them they separate all four arms — every pair of arms disagrees on at least one row, so the control tells you *which* arm you got, not merely that something changed:

```bash
printf 'NM_001234.1:c.[5_6insAC;9del]\nNM_001234.1:c.[2del;5del]\nNM_001234.1:c.[2del;9del]\n' > control.txt

for arm in live shadow canonical canonical-coalesced; do
  echo "== $arm"
  if FERRO_PARTITION=$arm ferro normalize --input control.txt --format tsv \
       --error-mode lenient > "control.$arm.tsv"; then
    tail -n +2 "control.$arm.tsv" | cut -f3
  else
    echo "   FAILED (exit $?) -- ferro did not produce this arm's column"
  fi
done
```

The status check is not boilerplate. An unrecognised arm now aborts the process, and a pipeline reports the exit status of its *last* command — so `ferro … | tail | cut` reports `cut`'s success and the loop prints an empty column under the arm's heading. An empty column and an aborted run look identical, which is the same "a broken measurement reads as a result" failure this whole section is about. Redirecting first and testing the status makes the abort say so.

| input | `live` | `shadow` | `canonical` | `canonical-coalesced` |
|---|---|---|---|---|
| `c.[5_6insAC;9del]` | `c.6_9delinsACCAA` | `c.[5_6insAC;11del]` | `c.[5_6insAC;11del]` | `c.[5_6insAC;11del]` |
| `c.[2del;5del]` | `c.[2del;6del]` | `c.[2del;6del]` | `c.2_4delinsG` | `c.2_4delinsG` |
| `c.[2del;9del]` | `c.[2del;11del]` | `c.[2del;11del]` | `c.[2del;11del]` | `c.2_9delinsGCCCAA` |

Row 1 separates `live` from the other three, row 2 separates `{live, shadow}` from `{canonical, canonical-coalesced}`, and row 3 separates `canonical` from `canonical-coalesced`.

If the arm you selected does not produce its column above, the variable is not reaching ferro and any comparison you run is meaningless. Only once the control behaves is a zero on your own corpus informative.

> The control this section used to give did not work. It offered `c.[2del;9del]`, `c.[3del;9del]`, `c.[2del;9dup]` and claimed `canonical` answers `c.[2del;33del]`; on those three inputs **no arm differs from `live`**, so the documented check failed for a correct setup and would have been read as "the variable is not reaching ferro". The table above is verified against both a debug and a release build.

**From Python, it must be set before the first normalization.** The value is read once per process and cached, so:

```python
import os
os.environ["FERRO_PARTITION"] = "canonical"   # must precede the first normalize call
import ferro_hgvs
```

Setting it after any variant has been normalized silently does nothing, and it cannot be changed within a running process. Comparing two rules therefore means two separate processes (or two CLI runs, as above), not two calls in one script.

## Why ferro-hgvs?

ferro-hgvs provides the most comprehensive HGVS variant normalization across all pattern types, with performance orders of magnitude faster than alternatives.

### Normalization Capabilities Comparison

<!-- DO NOT EDIT — generated from docs/tool_support_matrix.json by `generate_tool_support_tables`. -->
<!-- BEGIN tool-support:normalization_capabilities -->
| Pattern Type | ferro | mutalyzer | biocommons | hgvs-rs |
|--------------|:-----:|:---------:|:----------:|:-------:|
| Genomic (g.) | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) exonic | ✓ | ✓ | ✓ | ✓ |
| Coding (c.) intronic | ✓ | ✓** | ✗ | ✗ |
| Non-coding (n.) | ✓ | ✓ | ✓ | ✓ |
| RNA (r.) | ✓ | ✓ | ✓ | ✓ |
| Protein (p.) | ✓ | Net* | ✗ | ✗ |

\* mutalyzer protein normalization requires network access for NP_→NM_ lookups (cannot be cached locally).
\*\* mutalyzer intronic support is enabled by default via genomic-context rewriting; disable with --no-rewrite-intronic.
<!-- END tool-support:normalization_capabilities -->

### Performance Comparison

**All tools are benchmarked in ferro's offline configuration — best case for every tool.** Reference data is preloaded locally (a local UTA database and SeqRepo) and the network is disabled, so the figures below measure parse/normalize compute, not I/O. Out of the box, hgvs-rs, biocommons/hgvs, and mutalyzer resolve each variant against a remote UTA/SeqRepo or the Mutalyzer web API — a network round-trip per variant (~100–1000 ms), i.e. roughly **1–10 variants/sec, hundreds to thousands of times slower than shown here** (an order-of-magnitude estimate from per-call network latency, not separately benchmarked). That local, offline setup is exactly what ferro's `prepare` command builds; ferro needs no external service.

<!-- DO NOT EDIT — generated from data/benchmark/perf_results.json by `generate_perf_tables`. -->

_Median patterns/sec over 5 reps on an Apple M2 Max, local/offline. All tools draw from one stratified ClinVar population; per-tool sample sizes are calibrated so each tool is measured over a meaningful interval — fast cells (e.g. ferro/hgvs-rs parse) draw from millions of patterns, while slower cells (e.g. the per-tool normalize columns) draw from as few as tens to thousands. All tools exclude process/interpreter startup from the timed region — the mutalyzer/biocommons Python subprocesses are timed by their own internal startup-excluded timer, matching ferro/hgvs-rs. Only ferro parallelizes natively (rayon); the other tools are single-threaded libraries, so their normalize *@8 workers* figures come from the benchmark harness running 8 independent instances in parallel, while parsing is not sharded for them — hence the `single-threaded` label in their parse *@8 workers* column (mutalyzer normalize likewise shows no gain at 8 workers: per-call cache and IPC overhead dominate, so sharding does not help). Every tool runs fully offline against local reference data — a local UTA database and SeqRepo, with mutalyzer's network lookups disabled — the configuration ferro's `prepare` command enables; the figures therefore reflect compute throughput, not per-variant network latency. Reference-data load is excluded for all tools. ferro full-population peak: parse 20.0M/s, normalize 77.0k/s. See `docs/BENCHMARK_RUNBOOK.md` for the full method._

**Parse**

<!-- BEGIN perf:parse -->
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 5.1M/s | 12.2M/s | — |
| mutalyzer | 352/s | single-threaded | 35,000× |
| biocommons | 3.9k/s | single-threaded | 3,100× |
| hgvs-rs | 3.6M/s | single-threaded | 3× |
<!-- END perf:parse -->

**Normalize**

<!-- BEGIN perf:normalize -->
| Tool | Throughput @ 1 worker | Throughput @ 8 workers | ferro speedup @ 8w |
|------|----------------------:|-----------------------:|-------------------:|
| ferro | 78.1k/s | 260.2k/s | — |
| mutalyzer | 4/s | 4/s | 73,000× |
| biocommons | 368/s | 818/s | 320× |
| hgvs-rs | 195/s | 1.3k/s | 200× |
<!-- END perf:normalize -->

**ferro thread scaling**

<!-- BEGIN perf:ferro_scaling -->
| Threads | 1 | 2 | 4 | 8 |
|---------|--:|--:|--:|--:|
| ferro parse | 5.1M/s | 9.4M/s | 16.0M/s | 12.0M/s |
<!-- END perf:ferro_scaling -->

**Input ordering matters for batch throughput.** Resolving a transcript (reading its full sequence from the reference and rebuilding its CDS/exon metadata) dominates per-variant cost. ferro memoizes resolved transcripts in a bounded in-memory cache, so repeated lookups of the same transcript are near-free. Providing variants **sorted by transcript accession — or by genomic position, which clusters variants onto the same transcripts** — maximizes the cache hit rate and can speed up large batches by an order of magnitude versus randomly-ordered input. Ordering matters most when the number of distinct transcripts in the run exceeds the cache capacity (very large or genome-wide inputs); below that, the working set stays resident regardless of order.

### Reference Data: What ferro Prepares

The `ferro prepare` command downloads and organizes all reference data needed for comprehensive normalization. This data is then shared with other tools (mutalyzer, biocommons, hgvs-rs) to enable their local operation.

| Data Type | Source | Size | Enables |
|-----------|--------|------|---------|
| **RefSeq transcripts** | NCBI | ~1GB | NM_/NR_/XM_ normalization |
| **cdot metadata** | MANE | ~200MB | Transcript-to-genome mappings |
| **GRCh38 + GRCh37 genomes** | NCBI | ~4GB | NC_ genomic normalization |
| **RefSeqGene** (sequences + genome alignments) | NCBI | ~600MB | NG_ gene-region normalization; projecting c./n. variants into an NG_ parent's own g. frame (via the RefSeqGene→genome alignment GFF3) |
| **LRG sequences + XML** | EBI | ~50MB | LRG_ stable-reference normalization; projecting c./n. variants into an LRG_ parent's own g. frame (via the LRG XML genomic mapping) |
| **Protein sequences** | Derived from CDS | ~200MB | NP_/XP_ protein normalization |
| **Legacy transcript versions** | NCBI | ~50MB | Historical ClinVar variants |

**Key insight**: Without ferro's reference preparation, other tools require network access for each variant lookup (adding 100-1000ms latency per variant). With ferro's cached reference data, all tools can operate fully offline with consistent, reproducible results.

#### Deriving version-independent NG_ placements (#728)

`ferro prepare --derive-ng-placements <accessions.txt>` derives genomic placements for the listed `NG_` versions (one exact accession per line, e.g. `NG_012337.3`; blank lines and `#` comments ignored), writing `derived_refseqgene_placements.json` into the reference directory and wiring the manifest's `derived_refseqgene_placements` field. This fills version gaps the archived RefSeqGene→genome GFF3 snapshots do not cover. It needs cdot + the genome in the same prepare run and uses NCBI EFetch per accession; accessions that cannot be validated are skipped with a warning. The field is preserved across subsequent `prepare` runs.

## Benchmark: Reference Data & Tool Comparison

The main `ferro` binary includes commands to prepare reference data (`ferro prepare`) and check its status (`ferro check`). The `ferro-benchmark` tool (build with `--features benchmark`) extends this for tool comparison benchmarks.

| Command | Description |
|---------|-------------|
| `prepare <tool>` | Prepare reference data for a tool |
| `check <tool>` | Verify tool configuration and dependencies |
| `parse <tool>` | Parse HGVS patterns with specified tool |
| `normalize <tool>` | Normalize HGVS patterns with specified tool |
| `compare results` | Compare parse/normalize results between tools |
| `extract` | Extract patterns from ClinVar, VCFs, or create samples |
| `setup` | Set up UTA database, SeqRepo, and other services |
| `generate` | Generate summary reports and configs |
| `collate` | Aggregate sharded results |

### Quick Start

```bash
# Prepare ferro reference (main binary - no special features needed)
ferro prepare --output-dir data/ferro

# Check reference data
ferro check --reference data/ferro

# Normalize with ferro
ferro normalize -i patterns.txt --reference data/ferro

# For tool comparison, build with benchmark support
cargo build --release --features benchmark

# Prepare other tools (uses ferro reference for transcript data)
ferro-benchmark prepare mutalyzer --ferro-reference data/ferro --output-dir data/mutalyzer
ferro-benchmark prepare biocommons --seqrepo-dir data/seqrepo --uta-dump uta_20210129b.pgd.gz --ferro-reference data/ferro

# Compare results between tools
ferro-benchmark normalize mutalyzer -i patterns.txt -o mutalyzer.json --mutalyzer-settings data/mutalyzer/mutalyzer_settings.conf
ferro-benchmark compare results normalize ferro.json mutalyzer.json -o comparison.json
```

**Supported tools**: ferro-hgvs, mutalyzer, biocommons/hgvs, hgvs-rs

> **Note**: The `pixi.toml` and `pixi.lock` files in this repository define a [pixi](https://pixi.sh) environment for the Python-based external tools (mutalyzer, biocommons/hgvs, seqrepo) used in benchmarking. Run `pixi shell` to activate it.

See [docs/BENCHMARK_GUIDE.md](docs/BENCHMARK_GUIDE.md) for detailed usage.

## Development

```bash
cargo build
cargo test                        # default features
cargo clippy -- -D warnings
```

The commands above use the default feature set, and CI keeps them compiling (see
the `build` job). They do not cover the whole suite — the feature-gated tests and
the integration tree need `dev`, which is what CI runs and what you want before
opening a PR:

```bash
cargo nextest run --features dev
cargo clippy --features dev --all-targets -- -D warnings
```

## License

Licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Disclaimer

This software is under active development.
While we make a best effort to test this software and to fix issues as they are reported, this software is provided as-is without any warranty (see the [license](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/LICENSE) for details).
Please submit an [issue](https://github.com/fulcrumgenomics/ferro-hgvs/issues), and better yet a [pull request](https://github.com/fulcrumgenomics/ferro-hgvs/pulls) as well, if you discover a bug or identify a missing feature.
Please contact [Fulcrum Genomics](https://www.fulcrumgenomics.com) if you are considering using this software or are interested in sponsoring its development.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
