# Deriving a description from sequences

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

That is five values, not four. `max_grid_cells` is not inert — it decides whether an answer is
produced at all — so the "four arguments" this section and the Rust docs both used to claim is
withdrawn. Purity is the property; the count was wrong. (`FromSequencesOptions` carries a
`direction` too, but it is not a caller-facing knob: it is `#[doc(hidden)]`, always 3' on every
shipped path, and exists for the internal differential oracle described under rule 6.)

## How to: one canonical description per variant

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

### How wide should the window be?

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

### The refusals, and what to do about each

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

## Which rules it delivers

This is the whole design, and it falls straight out of the four [normalization rules](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/README.md#normalization-rules):

| | rules delivered | force | needs |
|---|---|---|---|
| `from_sequences` | **1 (conformant)**, **4 (deterministic)** | both *absolute* | the caller's four arguments |
| `normalize`, afterwards | **2 (recommended form)**, **3 (confluent)** | both *best effort* | the reference |

Rules 1 and 4 are the two the [normalization rules](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/README.md#normalization-rules) call always achievable, so a function that has only
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

## What you get for it

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

## The cost, stated plainly

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

## Going the other way

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

## Bounding a derivation to a region it must not leave

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
