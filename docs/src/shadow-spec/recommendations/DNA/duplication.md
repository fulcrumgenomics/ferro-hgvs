# Duplication — ferro's reading

ferro's reading of the HGVS **duplication** recommendations, clause by clause — each spelling with
the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Duplication (`r.`)](../RNA/duplication.md).*

On the coding axis (`c.`/`n.`) the exon/exon-junction 3'rule exception applies and a duplication is
not shifted across the junction. On the `r.` axis the same NOTE switches the exception off, which
ferro does not yet honour
([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211), on the RNA duplication page).

## `duplication.md:5` — definition

> Duplication: a sequence change where, compared to a reference sequence, a copy of one or more nucleotides is inserted **directly 3'** of the original copy of that sequence.

Ferro: a duplication inserts a copy of one or more nucleotides directly 3' of the original copy; the
"directly 3'-flanking" test is what separates a `dup` from an `ins` (the tandem-only rule at
`duplication.md:17-20`). The description carries only the duplicated range, never the duplicated
bases (the payload-spelling prohibition recurs at `duplication.md:35-36`, `:50-51`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5697dup` | recommended | self | single-nucleotide duplication — the last `A` of the 8-nucleotide A-stretch `c.5690_5697` (`c.5698` is `T`), so already a fixed point |
| `NC_000023.10:g.33229410dup` | recommended | — | the spec's own one-nucleotide genomic example, paired with `c.20dup` (foreign accession — parse-only here) |

## `duplication.md:15` — the range must name two different positions

> - `positions_duplicated` should contain **two different positions**, e.g., `123_126`, not `123_123`.

Ferro: a same-position range is repaired to the single-position form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5695_5697dup` | recommended | self | a correctly ranged three-nucleotide duplication — the 3'-most three of the A-stretch, so already a fixed point |
| `NM_004006.3:c.5697_5697dup` | refused | — | same-position range — rejected at parse in strict mode; lenient repairs to the single-position `c.5697dup` |

## `duplication.md:16` — the range is listed 5' to 3'

> - the `positions_duplicated` should be listed from **5' to 3'**, e.g., `123_126`, not `126_123`.

Ferro: the ordinary range is listed 5'→3'; a reversed, non-circular range is refused rather than
reordered. (The reversed `<high>_<low>` spelling is admitted only on an `o.`/`m.` circular reference
for an origin-spanning change; a coding reference is never circular.)

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5695_5697dup` | recommended | self | correctly ordered range (`c.5695` is 5' of `c.5697`) |
| `NM_004006.3:c.5697_5695dup` | refused | — | reversed range on a non-circular reference — rejected at parse |

## `duplication.md:17-20` — tandem-only; MUST be a dup not an ins; an inverted duplication is an insertion

> - by definition, duplication may only be used when the additional copy is **directly 3'-flanking** of the original copy (a "tandem duplication").
>     - when a variant can be described as a duplication, it **must** be described as a duplication and not as, e.g., an insertion (see [_Prioritization_](../general.md)).
>     - when there is no evidence that the extra copy of a sequence detected is in tandem (directly 3'-flanking the original copy), the change can not be described as a duplication; it should be described as **an insertion** (see [Insertion](insertion.md) and [proposal SVD-WG003](../../consultation/SVD-WG003.md)).
>     - **inverted duplications** are described as an insertion (`g.234_235ins123_234inv`), not as a duplication (see [Inversion](inversion.md)).

Ferro: four points. (i) a `dup` is only for a copy directly 3'-flanking the original. (ii) a
duplicating insertion **must** be relabeled `dup` — ferro derives the label from the resulting
sequence, so `c.19_20insT` where `c.19` is the copied base collapses to the `dup`. (iii) a copy
that is *not* directly 3'-flanking is an insertion, not a `dup` (worked at `duplication.md:150-154`).
(iv) an inverted duplication is spelled `ins<range>inv`, naming the source span, never `dup` — and
the `dupinv` shorthand is refused at parse. The "no evidence … in tandem" clause is an
assay-level evidentiary question a string-to-string normalizer cannot adjudicate; ferro never
manufactures tandem evidence the input did not assert.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[duplication-must-ranks-the-label-not-the-partition](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The rule that a duplication must be labelled 'dup' ranks the label of each piece ferro derives, not the partition; the one exception is a net-longer tandem copy of a multi-base motif, where the derivation is cut to expose the dup rather than merged into a delins.
>
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum — and it is wired today only in normalize_genome, so c./n./r./m. still emit literals.
<!-- why:END:duplication-must-ranks-the-label-not-the-partition,inverted-duplication-is-derived-as-ins-range-inv -->

</details>

**The inverted-duplication `ins<range>inv` form is derived only on the genomic axis** — `c.`/`n.`/`m.`
have no equivalent, so on the coding axis ferro **expands** an `ins<range>inv` range payload to
literal reverse-complement bases rather than preserving the range-citation spelling. That literal
output is valid HGVS that re-parses and denotes the same sequence, so the verdict is `conformant`,
not `bug` — a known limitation, with #2236 tracking the render-stage fix.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5690dup` | recommended | `NM_004006.3:c.5697dup` | the first `A` of the run relabeled and 3'-shifted to the last `A` — a duplicating change written as a `dup`, landing 3'-most |
| `NM_004006.3:c.234_235ins123_234inv` | conformant | `NM_004006.3:c.234_235insGTTGACATTGTTCAGGGCATGAACTCTTGTGGATCCTTTTTCTTTTGGCAGTTTTTGCCCTGTCAGGCCTTCGAGGAGGTCTAGGAGGCGCCTCCCATCCTGTAGGTCACTG` | the spec's own inverted-duplication spelling — an `ins<range>inv`. On the `c.` axis ferro **expands** the range payload to the 112 literal reverse-complement bases of `c.123_234` rather than preserving the `ins<range>inv` citation (`inverted-duplication-is-derived-as-ins-range-inv`; #2236 tracks the render-stage fix). Valid HGVS that re-parses and denotes the same sequence, so conformant — not the recommended spelling |
| `NM_004006.3:c.123_456dupinv` | refused | — | the `dupinv` shorthand — rejected at parse in strict mode; `:20` requires the `ins<range>inv` form |

## `duplication.md:21` — triplication and quadruplication use the repeated-sequence format

> - when more than one additional copies are inserted directly 3' of the original copy, the change is indicated using the format for [Repeated sequences](repeated.md), like `[3]` (triplication), `[4]` (quadruplication), etc.

Ferro: more than one additional copy is a repeat, spelled `[N]` (see `repeated.md`).

## `duplication.md:22-23` — individual vs. delins, and the one-codon exception

> - two variants separated by one or more nucleotides should be described individually and **not** as a "delins".
>     - **exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a "delins".<br>

Ferro: the separation rule (ruleset rule 2 — a preference, not a ban); the exception folds two
changes into one delins only when they sit one nucleotide apart and together change one amino acid.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[codon-carve-out-shape-restriction](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together alter a single amino acid are ruled one delins whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled — though the coding-codon pass currently implements only the substitution/unchanged-base/substitution triplet.
<!-- why:END:separation-rule-force-modal-or-negation,codon-carve-out-shape-restriction -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | two changes more than one nucleotide apart — described individually, kept as a cis allele, not merged |
| `NM_004006.3:c.[145C>T;147C>G]` | recommended | `NM_004006.3:c.145_147delinsTGG` | one nucleotide apart, together altering the `c.145_147` codon (`CGC`→`TGG`) — folded to one delins, whatever the member edit types |

## `duplication.md:24-26` — the 3'rule, and the exon/exon-junction exception (where ferro is correct on DNA)

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).
>     - **exception**: duplications around exon/exon junctions when identical nucleotides flank the junction (see [Numbering](../../background/numbering.md#DNAc));<br>
>       when `..GA`<code class="del">T</code>`gta..//..cagTCA..` changes to `..GA`<code class="ins">TT</code>`gta..//..cagTCA..`, based on a coding DNA reference sequence, the variant is described as `LRG_199t1:c.3921dup` (`NC_000023.10:g.32459297dup`) and not as `c.3922dup` (which would translate to `g.32456507dup`).

Ferro: the 3'rule shifts a duplication in a single-residue stretch or tandem repeat to its most-3'
position. The **exception** clamps that shift at an exon/exon junction: a `c.`/`n.` duplication whose
identical flanking bases would otherwise let it shift across the junction is held at the last
position that does not cross, reached from **either** side — a description approaching the junction
stops at it (`c.3921dup` is a fixed point), and one already spelled on the far side (`c.3922dup`) is
pulled back to `c.3921dup`.

**This is where the coding axis and the RNA axis diverge, and ferro is correct here.** On `c.`/`n.`
the exception applies, so `c.3922dup` does not translate to the wrong nucleotide in the wrong exon.
On `r.` the same NOTE *switches the exception off* and the 3'rule must cross the junction — which
ferro does not yet do
([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211), on the RNA duplication page).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[exon-junction-dup-converge-from-the-far-side](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A duplication is placed at the most 3' position that does not cross an exon/exon junction, reached from either side, so a copy spelled past the junction is pulled back to it.
<!-- why:END:exon-junction-dup-converge-from-the-far-side -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5690dup` | recommended | `NM_004006.3:c.5697dup` | 3'rule over a single-residue stretch: `c.5690_5697` is `AAAAAAAA`, so the duplication lands on the last `A` — the spec's own `duplication.md:38-39` example, on the slice's transcript version |
| `NM_004006.3:c.5697dup` | recommended | self | already the 3'-most base of that run — a fixed point |
| `NM_004006.3:c.5690_5692dup` | recommended | `NM_004006.3:c.5695_5697dup` | the same rule over a multi-nucleotide duplication inside the run — shifted to the 3'-most three |
| `LRG_199t1:c.3921dup` | recommended | — | the spec's own worked example — the duplication held at the exon/exon junction (parse-only here) |
| `LRG_199t1:c.3922dup` | recommended | — | the far-side spelling, which the exception pulls back to `c.3921dup` on the coding axis (parse-only here) |

## `duplication.md:27` — uncertain duplications

> - † = see [Uncertain](../uncertain.md); when the position and/or the sequence of a duplication has not been defined.

Ferro: an uncertain duplication's parenthesised range is preserved verbatim — the 3'rule has nothing
determinate to shift when neither the exact position nor the exact length is asserted.

## `duplication.md:31-46` — one nucleotide (three worked examples)

>     - **`NM_004006.2:c.5697dup` (3'rule)**<br>

Ferro: the duplicated base is never spelled (`c.20dupT`/`c.20dupG` are `class="invalid"`); the
3'rule shifts a duplication within a run of identical bases to its most-3' position, including
under minus-strand mapping. Describing the variant as an insertion (`c.19_20insT`) is well-formed
HGVS and parses, but the `:18` MUST is honoured by *relabeling* it to `dup` at normalization rather
than rejecting it at parse.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.20dup` | recommended | — | the plain one-nucleotide example; the spec's transcript version (parse-only here) |
| `NM_004006.2:c.19_20insT` | recommended | — | the insertion spelling of a duplicating change, which `:18`'s MUST forbids — but the string is well-formed HGVS, so it **parses** in every mode; the MUST is a normalization-level relabel to `c.20dup`, not a parse rejection (parse-only here; the executable twin two rows down shows the relabel) |
| `NM_004006.3:c.5697_5698insA` | recommended | `NM_004006.3:c.5697dup` | executable twin of the `:18` MUST: an `A` inserted directly 3' of `c.5697` (the last `A` of the run, `c.5698` is `T`) is a tandem copy, so ferro relabels the insertion as the duplication it is |
| `NM_004006.2:c.20dupT` | refused | — | the disallowed payload-bearing spelling (`class="invalid"`) (parse-only here) |
| `NM_004006.3:c.5697dup` | recommended | self | executable twin on the slice's version — the last `A` of `c.5690_5697`, a fixed point |
| `NM_004006.3:c.5697dupA` | refused | — | executable twin — strict rejects the stated duplicated base; lenient drops the payload and normalizes to `c.5697dup` |
| `NC_000023.11:g.32343183dup` | recommended | — | the minus-strand cross-check (foreign accession — parse-only here) |

## `duplication.md:47-55` — several nucleotides

>     - **`NM_004006.2:c.20_23dup` (`NC_000023.11:g.33211290_33211293dup`)**<br>

Ferro: a multi-nucleotide range never carries the duplicated payload sequence either
(`c.20_23dupTAGA` is `class="invalid"`), and an ordinary range crossing an exon/intron border
(`c.260_264+48dup`) needs no special handling while it stays within the transcript.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.20_23dup` | recommended | — | plain range, no duplicated sequence spelled; the spec's transcript version (parse-only here) |
| `NM_004006.2:c.20_23dupTAGA` | refused | — | the disallowed payload-bearing range spelling (`class="invalid"`) (parse-only here) |
| `NM_004006.3:c.5695_5697dup` | recommended | self | executable several-nucleotide twin — the 3'-most three of the A-stretch, nothing left to shift |
| `NM_004006.3:c.5695_5697dupAAA` | refused | — | executable twin — strict rejects the stated duplicated sequence; lenient drops it and normalizes to `c.5695_5697dup` |
| `NC_000023.11(NM_004006.2):c.260_264+48dup` | recommended | — | range crossing an exon/intron border (parse-only here) |

## `duplication.md:56-71` — exon/intron/exon

>         - **`NC_000023.11(NM_004006.2):c.3921dup`**<br>

Ferro: three border cases. The exon/exon case (`c.3921dup`) is the exon/exon-junction exception
already adjudicated at `duplication.md:24-26`. The exon/intron (`c.1704+1dup`, not `c.1704dup`) and
intron/exon (`c.1813dup`, not `c.1813-1dup`) cases are ordinary "spell the duplicated base where it
unambiguously sits" — no exception fires, because no run of identical bases straddles the border. The
"does not depend on the effect observed on RNA level" NOTE on both is a scope statement: DNA-level
normalization is independent of any predicted splicing consequence.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):c.3921dup` | recommended | — | exon/exon border — the exception holds the duplication at `c.3921`, not `c.3922` (parse-only here; see `duplication.md:24-26`) |
| `NC_000023.11(NM_004006.2):c.1704+1dup` | recommended | — | exon/intron border — the intronic `G` spelled as `c.1704+1`, not `c.1704dup` (parse-only here) |
| `NC_000023.11(NM_004006.2):c.1813dup` | recommended | — | intron/exon border — the exonic `G` spelled as `c.1813`, not `c.1813-1dup` (parse-only here) |

## `duplication.md:72-103` — exons, the [dup;ins] split, and the beyond-transcript prohibition

>     - **`NC_000023.11(NM_004006.2):c.4072-1234_5155-246dup`**<br>

Ferro: multi-exon duplications use concrete or uncertain (parenthesized) positions that stay within
the transcript, and the size is never appended. A `dup` immediately followed by an `ins` at the
very next position is spelled as **two separate members** (`c.[…dup;…ins…]`), never merged into a
`dupins` composite — a format the spec calls "not used in HGVS nomenclature". A duplication
extending beyond the transcribed region **cannot** be described on a coding reference — a rank-1
prohibition; genomic coordinates are required. Ferro enforces this for a concrete `c.` position
(rejected in strict mode), but the spec's `class="invalid"` forms reach past the boundary through
an open `?` end, which ferro accepts.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):c.4072-1234_5155-246dup` | recommended | — | concrete multi-exon duplication (parse-only here) |
| `NM_004006.2:c.4072-1234_5155-246dup` | refused | — | the coding reference does not include the intron sequences the range spans (`class="invalid"`) (parse-only here) |
| `NC_000023.11(NM_004006.2):c.720_991dup` | recommended | — | plain multi-exon duplication (parse-only here) |
| `NC_000023.11(NM_004006.2):c.(4071+1_4072-1)_(5154+1_5155-1)dup` | recommended | — | uncertain-breakpoint form — the spec tags this as part of the still-undecided proposal SVD-WG003, so no settled ferro behaviour is asserted (parse-only here) |
| `NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]` | recommended | — | a `dup` beside an `ins`, spelled as two members — never `dupins` (parse-only here) |

## `duplication.md:104-126` — gene, chromosome, and mosaic/chimeric

> **NOTE**: the array analysis detects an extra copy of the sequences, and it has to be determined whether it is a duplication.
> When it is not sure the variant is a duplication, the variant should be described as an insertion; `g.?_?ins[NC_000023.11:g.(31060227_31100351)_(33274278_33417151)]`.

Ferro: the gene-, chromosome- and mosaic/chimeric-level rows are cross-document applications, not
duplication-local rules. The whole-gene SNP-array/MLPA forms
(`g.(31060227_31100351)_(33274278_33417151)dup`) are answered by the general uncertain-range parser;
the "when it is not sure it is a duplication, describe it as an insertion" NOTE is assay-evidentiary
language a sequence-string normalizer cannot adjudicate (there is no "sureness" input). The
chromosome-level `sup` (supernumerary) and `[2]` forms belong to `complex.md`/`repeated.md`, and the
mosaic (`=/dup`) and chimeric (`=//dup`) operators belong to the general mosaicism syntax — the `dup`
member inside carries the same 3'rule as a standalone duplication.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.33344590_33344592=/dup` | recommended | — | mosaic form, spec's own accession (parse-only here) |
| `NC_000023.11:g.33344590_33344592=//dup` | recommended | — | chimeric form, spec's own accession (parse-only here) |

## `duplication.md:129-139` — why not describe a duplication as an insertion

> !!! note "Why do we not describe a duplication as an insertion?"

Ferro: rationale prose for the dup-over-ins preference (brevity, position clarity, avoiding
intron/exon insertion-site debates, the local-slippage mechanism). It states no rule beyond the
`:18` MUST already covered at `duplication.md:17-20`.

## `duplication.md:140-144` — can I use `g.123dup6`?

> !!! note "Can I use <code class="invalid">g.123dup6</code> to describe a 6 nucleotide duplication?"

Ferro: `class="invalid"` — a duplication of more than one nucleotide gives the first and last
position separated by `_`, never a length suffix, and `g.123dup6` is additionally ambiguous about
whether the duplication starts *at* or *after* position 123. Strict rejects the length suffix.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NG_012232.1:g.123_128dup` | recommended | — | the correct range form for a six-nucleotide duplication (parse-only here) |
| `NG_012232.1:g.123dup6` | refused | — | the disallowed length-suffix spelling (`class="invalid"`) (parse-only here) |
| `NM_004006.3:c.5692_5697dup` | recommended | self | executable twin of the correct form — the 3'-most six of the A-stretch (`c.5698` is `T`), so a fixed point |

## `duplication.md:145-149` — why `c.3921dup` and not `c.3922dup`

> !!! note "In the example above, **`c.3921dup`**, should the description based on a coding DNA reference sequence not be `c.3922dup`?"

Ferro: this restates the exon/exon-junction exception (`duplication.md:24-26`) as a Q&A — the
far-side `c.3922dup` is pulled back to `c.3921dup` so that translating the position back does not
land in the wrong nucleotide in the wrong exon (`exon-junction-dup-converge-from-the-far-side`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.3922dup` | recommended | — | the far-side spelling, pulled back to `c.3921dup` by the clamp (parse-only here; see `duplication.md:24-26`) |

## `duplication.md:150-154` — a coincidental copy is an insertion, not a duplication

> !!! note "How should I describe the change `ATCG`<code class="spot1">ATCGATCGATCG</code><code class="spot2">A</code>`GGGTCCC` to `ATCG`<code class="spot1">ATCGATCGATCG</code><code class="spot2">A</code><code class="ins">ATCGATCGATCG</code>`GGGTCCC`? The fact that the inserted sequence (<code class="ins">ATCGATCGATCG</code>) is present in the original sequence, suggests it derives from a duplicative event."

Ferro: the copied span (`g.5_16`) is offset one base from the insertion point (`g.17|18`), so the
payload is not directly 3'-flanking and is correctly not mislabeled `dup`. The spec's recommended
form is `g.17_18ins5_16`, a range-citation naming the likely source span. Ferro does not derive
that range-citation form for a plain (non-inverted) coincidental insertion — it expands the range
to literal bases instead. That literal insertion is valid HGVS denoting the same sequence, so the
verdict is `conformant`, not the recommended spelling.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.17_18ins5_16` | conformant | `NM_004006.3:c.17_18insTTTGGTGGGAAG` | executable twin of the spec's own answer `g.17_18ins5_16` — an insertion naming its likely source range, **not** a `dup`, because the copy is not directly 3'-flanking. Ferro does not relabel the non-3'-flanking copy as `dup` (correct), and on the `c.` axis it **expands** the range payload to the 12 literal bases of `c.5_16` rather than preserving the `ins5_16` citation — valid HGVS denoting the same sequence, so conformant, not the spec's recommended range-citation form |
