# Inversion — ferro's reading

ferro's reading of the HGVS **inversion** recommendations (DNA axes), clause by clause — each
spelling with the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Inversion (`r.`)](../RNA/inversion.md).*

## `inversion.md:5` — definition: more than one nucleotide, reverse complement

> Inversion: a sequence change where, compared to a reference sequence, **more than one nucleotide** replacing the original sequence is the reverse complement of the original sequence.

Ferro: an inversion replaces a span of **more than one** nucleotide with its exact reverse
complement. When a span's whole content is its reverse complement, that span types as one `inv`
uniformly — whatever the competing split partition's members would type as (substitutions, delins,
or a mix), and however much of the interior happens to coincide with the reference. The property is
a fact about the whole span, with no term for its interior columns. Ferro mints this typing only on
the genomic axis; on `c.`/`n.` ferro preserves an already-spelled inversion but does not turn a
reverse-complement split into one.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[whole-span-reverse-complement-types-as-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A span whose whole content is replaced by its exact reverse complement is written as one inv, however much of its interior coincides with the reference and whatever the competing partition is made of; this is a project choice among conformant forms, not a conformance requirement.
<!-- why:END:whole-span-reverse-complement-types-as-inv -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5_8inv` | recommended | `NM_004006.3:c.5_8inv` | A four-nucleotide inversion with every column changing: `c.5_8` is `TTTG`, replaced by its reverse complement `CAAA`, so the competing split is four substitutions (`c.[5T>C;6T>A;7T>A;8G>A]`) — and `:5` types it as one `inv` regardless |
| `NM_004006.3:c.123_124inv` | recommended | `NM_004006.3:c.123_124inv` | The minimal >1-nt inversion, and a whole-span case *with* a competitor: `c.123_124` is `CA`, reverse complement `TG`, so both columns are substitutions (`C>T`, `A>G`) — the competing split is two subs, and `:5` types it as one `inv` regardless |
| `NM_004006.3:c.5690_5697inv` | recommended | `NM_004006.3:c.5690_5697inv` | The whole-span rule at its starkest, on the slice's own 8-nucleotide A-stretch: `AAAAAAAA` → reverse complement `TTTTTTTT`, every column an `A>T` substitution — the competing split is eight subs, and it still types as one `inv` |
| `NC_000023.10:g.32361330_32361333inv` | recommended | — | the spec's own four-nucleotide genomic example (`:27-28`), `..CA`TCAG`CCT..` → `..CA`CTGA`CCT..` — reverse complement of `TCAG` is `CTGA`, all four columns change (foreign accession — parse-only here) |
| `NM_004006.2:c.5657_5660inv` | recommended | — | the spec's coding example (`:30-31`), `CTGA` → `TCAG` — the wrong transcript version for the slice, so parse-only here |
| `NM_004006.2:c.4145_4160inv` | recommended | — | the spec's load-bearing 16-nucleotide example (`:33-34`) with a large unchanged interior — re-derived by the ledger as 10 of 16 columns changed, two three-base unchanged interior runs, yet typed as one `inv`. The single strongest textual evidence that the spec does not decompose an inversion by interior coincidence (wrong version — parse-only here) |
| `NC_000023.10:g.111754331_111966764inv` | recommended | — | the spec's 212,434-nucleotide example (`:36-37`) — the definition states no length ceiling; a 212-kb span types identically to a 4-nt one (foreign accession — parse-only here) |

## `inversion.md:15-16` — more than one nucleotide; a one-nucleotide inversion is a substitution

> - by definition, the region inverted (`positions_inverted`) contains **more than one nucleotide**.
>   The description <code class="invalid">g.234inv</code> is therefore not allowed; a one-nucleotide inversion should be described as a [substitution](substitution.md).

Ferro: two points, both mechanical. (i) An inversion spans more than one nucleotide, so the grammar
admits only a *range* (`range "inv"`) — a lone single-position `inv` has no production and is
refused at parse. (ii) A one-nucleotide inversion should be written as a substitution, which is
definitional: the reverse complement of a single base is its complement, i.e. a `sub`, so there is
no lone-`inv` spelling for ferro to relabel. The spec's own `g.234inv` is `class="invalid"`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.234inv` | refused | — | the executable twin of the spec's `class="invalid"` single-position `g.234inv` — rejected at parse: an inversion requires a range of more than one nucleotide, so a lone position carries no span to invert |

## `inversion.md:17-18` — the 3'rule

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).
>     - the 3'rule applies to ALL descriptions (genome, gene, transcript, and protein) of a given variant.

Ferro: the general 3'-shift rule applies to inversions exactly as to the other edit types — an
inversion writable at more than one equivalent position is assigned its most 3' placement. An
inversion 3'-shifts only inside a larger reverse-complement-symmetric context, which this
transcript's homopolymer runs do not supply, so the executable row below is a fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5_8inv` | recommended | `NM_004006.3:c.5_8inv` | `c.5_8` is `TTTG`; no equivalent 3' placement exists for this span, so the 3'rule has nothing to shift — already 3'-most |

## `inversion.md:19` — an inverted duplication is `ins<range>inv`, not `dupinv`

> - **inverted duplications** are described as an insertion (`g.234_235ins123_234inv`), not as <code class="invalid">g.123_456dupinv</code> (see [Q&A](#dupinv)).

Ferro: an inverted duplication is spelled `ins<range>inv`, naming the span the inverted copy came
from — never `dup`, and never the `dupinv` shorthand, which the grammar has no production for and
is refused at parse. That `ins<range>inv` form is the ruling's target, and ferro derives it **only
on the genomic axis**. On `c.`/`n.` ferro does not yet mint or preserve it — it **expands** an
already-spelled range payload to literal reverse-complement bases, which is valid HGVS that
re-parses and denotes the same sequence, so `conformant`, not the recommended spelling. This is an
implementation gap, not a disagreement with the ruling; #2236 tracks the render-stage fix.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum — and it is wired today only in normalize_genome, so c./n./r./m. still emit literals.
<!-- why:END:inverted-duplication-is-derived-as-ins-range-inv -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.16_17ins5_16inv` | conformant | `NM_004006.3:c.16_17insCTTCCCACCAAA` | An inverted duplication of `c.5_16` (`TTTGGTGGGAAG`) inserted 3' of the original: on `c.` ferro expands the range-`inv` payload to its literal reverse complement (`CTTCCCACCAAA`). Valid HGVS that re-parses and denotes the same sequence, so conformant, not the recommended `ins5_16inv` spelling |
| `NM_004006.3:c.123_456dupinv` | refused | — | the `dupinv` shorthand `:19` marks `class="invalid"` — rejected at parse: `dup` has no `inv` suffix production, and the clause requires the `ins<range>inv` form instead |
| `NM_004006.2:c.849_850ins850_900inv` | conformant | — | the spec's inverted-duplication example (`:39-40`), copy inserted 5' of the original. On a `c.` reference ferro expands the range-`inv` payload to literal reverse-complement bases (as on the executable `c.16_17ins5_16inv` row above), so its output is not the recommended `ins<range>inv` form — conformant (wrong version — parse-only here) |
| `NM_004006.2:c.900_901ins850_900inv` | conformant | — | the spec's mirror example (`:42-43`), copy inserted 3' of the original — the two orientations are why the genomic derivation must produce both spellings rather than normalizing to one; on `c.` ferro instead expands to literals, so conformant (wrong version — parse-only here) |
| `LRG_199t1:c.940_941ins[885_940inv;A;851_883inv]` | conformant | — | the spec's composite example (`:45-46`): a bracketed insertion payload of two inverted ranges and one literal base (`c.884` `G>A`) inside the copied interior. A payload shape wider than the atomic `ins<range>inv` derivation covers, so on `c.` ferro cannot keep the range form — conformant (foreign accession — parse-only here) |
| `NM_004006.2:c.940_941ins[903_940inv;851_885inv]` | conformant | — | the spec's second composite example (`:48-49`): two non-adjacent inverted ranges, representing a deletion inside the copied-and-inverted interior. A payload shape wider than the atomic `ins<range>inv` derivation covers, so on `c.` ferro cannot keep the range form — conformant (wrong version — parse-only here) |

## `inversion.md:20-21` — individual vs. delins, and the one-codon exception

> - two variants separated by one or more nucleotides should be described individually and **not** as a "delins".<br>
>   **exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a "delins".<br>

Ferro: this restates `general.md:33`/`:34`'s separation-and-codon clause inside the inversion doc.
Its only inversion-specific function is naming the competitor when a single whole-span inversion
competes with a multi-member split — and there the typing rule at `:5` governs, because a
whole-span inversion is *one* variant, so the separation rule has no antecedent (no "two variants")
to bite on until a split is already chosen. Two competitor-shape records are cited below; only the
pure-`delins` competitor has an executable row here (`c.76_83inv`), while the mixed
substitution/multi-column competitor is cited for completeness with no example on the committed
slice. Both outcomes stand — a whole-span reverse complement types as one `inv` — but their
competitor-type reasoning is **superseded** by the uniform `:5` rule above.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[inversion-vs-two-delins-76-83](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A span replaced by its exact reverse complement is written as a single inv even where its interior columns coincide with the reference, not as the two delins those columns would separate.
>
> **[inversion-vs-a-mixed-member-competitor](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When a span replaced by its reverse complement competes with a description mixing lone substitutions and multi-column members, ferro writes it as one inv; both forms are conformant, so this is the project's choice among them.
<!-- why:END:inversion-vs-two-delins-76-83,inversion-vs-a-mixed-member-competitor -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76_83inv` | recommended | `NM_004006.3:c.76_83inv` | The `inversion-vs-two-delins-76-83` record's own worked span: an 8-nucleotide whole-span reverse complement whose competing split would be two `delins` members (`c.[76_77delins…;82_83delins…]`) — typed and preserved as one `inv`, not split. `c.76` is `A` on the slice; the remaining bases are read from the reference at bless time |

## `inversion.md:22-23` — inversions are not used on the protein level

> - inversions are not used on protein level.
>   Depending on the (predicted) consequences of an inversion on protein level, changes are usually described as either a **delins** or a **frameshift**.

Ferro: an absolute prohibition, and it is grammar-enforced rather than merely a preferred
description — the protein (`aa:`) grammar has no `inv` edit type at all (only `alleles`, `del`,
`delins`, `dup`, `ext`, `fs`, `ins`, `rpt`, `sub`), so no production a `p.` inversion could parse
as exists. A protein-level consequence of an inversion is expressed as a `delins` or a frameshift
instead.

## `inversion.md:53-56` — Discussion: the complement is not an inversion

> !!! note "Is the change `AAGC` to `TTCG` an inversion?"

Ferro: a clarification, not a rule ferro applies. It reinforces the `:5` definition — an inversion
is the **reverse** complement, so `AAGC` inverts to `GCTT`, whereas `TTCG` is only the complement
(bases flipped, order kept).

## `inversion.md:59-62` — Discussion: the reverse is not an inversion

> !!! note "Is the change `AAGC` to `CGAA` an inversion?"

Ferro: the mirror-image clarification of the one above — `CGAA` is only the reverse (order flipped,
bases kept), not the reverse **complement** (`AAGC`'s reverse complement is `GCTT`). Both Q&A blocks
fence the definition against its two near-misses (complement-only and reverse-only).

## `inversion.md:64-69` — Discussion: why `ins<range>inv` and not `dupinv`

> !!! note "Is it not better to describe the variant `g.234_235ins123_234inv` as <code class="invalid">g.123_234dupinv</code>?"

Ferro: the spec's fullest statement of *why* `dupinv` is rejected — an inverted copy is not "a copy
inserted directly 3' of the original," so by `duplication.md`'s definition it is not a duplication
but an insertion — and of the two-spelling asymmetry (`g.122_123ins123_234inv` or
`g.234_235ins123_234inv` depending on whether the inverted copy is 5' or 3' of the original). Same
governing record as `:19` (`inverted-duplication-is-derived-as-ins-range-inv`); the verdict rows
are there.
