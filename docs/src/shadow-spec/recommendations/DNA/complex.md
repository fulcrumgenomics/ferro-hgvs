# Complex (HGVS/ISCN) — ferro's reading

ferro's reading of the HGVS **complex (HGVS/ISCN)** recommendations, clause by clause — each
spelling with the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

This page is almost entirely the ISCN named extension — translocations, ring chromosomes,
pericentric inversions, supernumerary/marker chromosomes and the symbol conventions those need.
Those descriptions are **genomic**: they sit on chromosome (`NC_`), `chr`, `NG_` or foreign
accessions, none of which the committed reference slice carries (it holds one transcript,
`NM_004006.3`). Complex descriptions on `NM_004006.3` do not occur — the shape is chromosomal — so
this page fabricates none, and **every verdict row here is parse-only** (`Normalizes to` = `—`): the
verdict is what ferro's strict parse does with the string (`refused` for a form it rejects,
`recommended`/`conformant` for one it accepts), read without reference bases. The refusals that key
on specific ISCN spellings — comma-thousands coordinates, the removed cross-chromosome `::`, the
`dupinv` form — are the parser's own grammar rejections, in every mode. No formal `rulings` record
covers them; what does exist is ferro's per-input override table (`by_input`, #546), which records
for the spec-fixture census that the spec *itself* expects its harvested ISCN / SVD-WG004 strings
to be rejected, so they count as `correctly-rejected` rather than `parse-error`. The spec text is
unambiguous enough that a rejection is plainly conformant, so a dedicated adjudication record is
optional and none is invented here.

## `complex.md:5` — definition: the residual "complex" bucket

> Complex: a sequence change where, compared to a reference sequence, a range of changes occur that can not be described as one of the basic variant types (substitution, deletion, duplication, insertion, inversion, deletion-insertion, or repeated sequence).

Ferro: "complex" is the residual category — what remains when no basic variant type fits. It states
no behavioural obligation of its own, but it is load-bearing elsewhere: it is precisely because a
ring chromosome *has* no single-type description that the self-cancelling prohibition does not reach
across its `::` junctions (see `complex.md:130` below). Descriptive.

## `complex.md:12-13` — ISCN describes structure, HGVS describes the variants

> It should be noted there is a basic difference between ISCN and HGVS: while ISCN describes the structure of the resulting chromosome(s), HGVS describes the **variant(s) detected**.
> It should be noted that the description of complex changes can become rather complicated and at some point, although literally correct, becomes effectively meaningless.

Ferro: background framing, and a caution the ledger leans on. ISCN describes the resulting
chromosome; HGVS describes the variants. The `:13` warning — that a complex description can be
"literally correct" yet "effectively meaningless" — is why a structural rejection rule is not minted
from a pair of worked examples plus biology (see the open anchoring question at `complex.md:28`).
Descriptive.

## `complex.md:46-47` — nucleotide positions with comma-thousands are not HGVS

> - in ISCN, it is allowed to describe nucleotide positions using commas to indicate thousands and millions (e.g., "108,111,982").
>   In HGVS, this is not allowed.

Ferro: a flat prohibition — one of the few strong statements in the file. A coordinate written with
comma-thousands separators is an ISCN rendering, not a valid HGVS string, and is refused — the
parser fails at the first comma, in every mode. `by_input` (#546) records the spec's own ISCN
example strings as expected-rejected, cited to this clause, so the fixture census reads them as
`correctly-rejected`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `chr2:g.pter_8,247,756::chr11:g.15,825,273_cen_qter` | refused | — | the ISCN rendering of the balanced `t(2;11)` (`complex.md:75`): comma-thousands coordinates on a `chr` accession; not an HGVS string. Rejected at parse (the first comma), in every mode; recorded as expected-rejected in `by_input` (#546), cited to `complex.md:46-47` |

## `complex.md:60-64` — translocations are described as `delins`, not with `::`

> HGVS recommendations try to avoid such conflicts wherever possible.
> HGVS, therefore, recommends to describe translocations exclusively using a "delins" format.

Ferro: ISCN2020 narrowed `::` to the ring-junction use only; the earlier cross-chromosome `::`
(translocation/transposition) use was removed. The spec therefore recommends translocations
"exclusively" as `delins[ACC:g…]` — the adjacent `delins` row is the HGVS-proper form, and the
cross-chromosome `::` string is not. The spec spells the removed form out itself and marks it
`class="invalid"` (`complex.md:72`). Ferro accepts the `delins` form (strict and lenient alike) and
refuses the cross-chromosome `::` composite at parse, in every mode — `::` is admitted only between
same-accession segments (the ring form at `complex.md:130`), so a second accession after it is a
grammar error. The ISCN `chr`-form renderings of the same translocations are recorded as
expected-rejected in `by_input` (#546), cited to `complex.md:38-41` and `complex.md:60-64`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825272]` | recommended | — | the balanced `t(2;11)` as the spec's `delins` form (`complex.md:70`); foreign genomic accessions — parse-only here |
| `NC_000002.12:g.pter_8247756::NC_000011.10:g.15825273_cen_qter` | refused | — | the removed cross-chromosome `::` form (`class="invalid"` at `complex.md:72`) — ISCN2020 dropped `::` for translocations, so this is refused and the `delins` form above is used instead |

## `complex.md:130` — `::` joins a ring; it is not `;`

> **NOTE**: "::" is used to indicate the join, instead of ";" to describe two not connected deletions.

Ferro: `::` carries connectivity that `;` does not. Two `del` segments joined by `::` form a ring
(a fused-break topology); the same two joined by `;` are two independent deletions. A `::` composite
is therefore **not** reducible to its member set, so `general.md:57`'s self-cancelling prohibition —
a sub-bullet of the type-prioritisation rule at `general.md:55`, which can only redirect to a
competing single-type description — cannot reach across a ring's `::` junctions, because a ring (being
complex by definition, `complex.md:5`) has no such single-type description. Ferro accepts ring
composites; the `ring_segment_wellformedness` and `issue_1578_followup_self_cancelling_rings` guards
exercise this.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[self-cancelling-across-ring-junctions](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The prohibition on replacing part of a sequence with part of itself does not reach the '::'-joined segments of a ring chromosome, so ferro keeps a ring's members rather than collapsing them to a linear delins.
<!-- why:END:self-cancelling-across-ring-junctions -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000022.11:g.pter_(12200001_14700000)del::(37600001_41000000)_qterdel` | recommended | — | the ring-chromosome form (`complex.md:127`): two `del` segments joined by `::`, which `;` cannot express; foreign genomic accession — parse-only here |
| `NC_000022.11:g.[pter_(12200001_14700000)del::(37600001_41000000)_qterdel]sup` | recommended | — | the supernumerary ring (`sup`, `complex.md:161`) — an additional ring not attached to other chromosomal material; parse-only here |

## `complex.md:28` — ring telomere anchoring is undecided

> the start of the chromosome is described as `pter`, the end as `qter`, and the centromere as `cen`.

Ferro: the clause defines what `pter`/`qter`/`cen` *mean*; it does not command their presence. The
open question is whether a ring's first `::` segment must start at `pter` and its last end at `qter`,
so a ring naming only interior coordinates is refused. The spec publishes only anchored rings
(`complex.md:127`, `complex.md:161`) and biology agrees a ring loses both telomeres — but no clause
states the requirement, and `complex.md:13` counsels against reading a rejection rule out of two
worked examples. `complex.md:30`'s `(pter)_#` / `#_(qter)` partial-anchor forms and `cen` further
complicate any simple "first endpoint must be `pter`" predicate. **This is undecided.** Ferro
**accepts** an unanchored ring today — the status quo, pinned by
`a_ring_with_no_telomere_anchor_is_still_accepted`, not a ruling either way. No verdict row: the
status quo is acceptance, and no fix direction is taken.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[ring-telomere-anchoring](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Whether a ring chromosome's '::'-joined segments must be anchored at 'pter' and 'qter' is undecided — no clause states the requirement, only biological reasoning and two worked examples support it — and ferro's current acceptance of an unanchored ring is the unresolved status quo, not a ruling either way.
<!-- why:END:ring-telomere-anchoring -->

</details>

## `complex.md:49-55` — the 3' rule and breakpoint ordering

> to determine the location of the break point, the general HGVS rule of maintaining the longest unchanged sequence applies (the 3' rule).
> Break point location is determined by the first break point encountered, i.e. `pter` of the chromosome is to be listed first.

Ferro: the ISCN-structural placement/ordering conventions. `complex.md:50` restates the general 3'
rule for breakpoint placement (consistent with ferro's general shifting); `complex.md:51`,
`complex.md:53` and `complex.md:55` fix segment order (`pter` to `qter`) and forward orientation.
These are consistent with ferro's existing behaviour and add no new obligation — normative-supporting
but effectively descriptive on this page.

## `complex.md:134-140` — an orientation-reversed tandem duplication is `ins(range)inv`

> - **`NC_000008.11:g.(131500001_136400000)ins(127300001_131500000)_(131500001_136400000)inv`**<br>
>   for ISCN `dup(8)(q24.22q24.21)`<br>
>   _(within a chromosome, orientation reversed relative to original sequence, breakpoint not sequenced)_

Ferro: a duplication whose copy is inserted in reversed orientation is spelled `ins<range>inv`,
naming the source span — never as `dupinv`. `dupinv` is `class="invalid"` in `inversion.md:19`, and
the primary authority for the derived form is `inversion.md:69`. Ferro derives the `ins<range>inv`
form rather than expanding the payload to literal reverse-complement bases, and refuses `dupinv`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum, so a short chance reverse-complement match is not misread as one.
<!-- why:END:inverted-duplication-is-derived-as-ins-range-inv -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000008.11:g.(127300001_131500000)_(131500001_136400000)dup` | recommended | — | the same-orientation tandem duplication (`complex.md:134`); foreign genomic accession — parse-only here |
| `NC_000008.11:g.(131500001_136400000)ins(127300001_131500000)_(131500001_136400000)inv` | recommended | — | the orientation-reversed duplication as `ins<range>inv`, naming the source span (`complex.md:138`); parse-only here |
| `NC_000008.11:g.(127300001_131500000)_(131500001_136400000)dupinv` | refused | — | the `dupinv` spelling for the same event — `class="invalid"` at `inversion.md:19`; rejected at parse (trailing `inv` after `dup`), in every mode. The SVD-WG004 harvest of this event (`chr8:…dupinv`) is recorded as expected-rejected in `by_input` (#546); the `ins<range>inv` form above is used instead |

## `complex.md:113-117` — pericentric inversions as `[;]` cis alleles

> - **`NC_000006.12:g.[776788_93191545inv;93191546T>C]`**<br>

Ferro: worked examples spelling a complex rearrangement as a `[;]` cis allele of basic edit types —
an `inv` beside a substitution (`complex.md:113`), and `del` + `inv` + `ins` (`complex.md:117`).
Their multi-megabase members are genuinely distinct events, so there is no separation/merge tension
of their own; they corroborate that connected complex changes are `;`-joined cis alleles of basic
types. Descriptive.

## `complex.md:169-175` — Discussion: ad-hoc labels are refused; `N[count]` is legal

> No, not really, it is not exact.

Ferro: three adjudicative points from the Q&A. First, ad-hoc labels (`TSD`, `insL1.603bp`) are not
valid HGVS — a description must be exact — and are refused. Second, the discussion restates the
duplication definition (`complex.md:173`): `dup` is used only when the duplicate sits directly 3′ of
the original copy, which corroborates that the `dup` label ranks the change itself, not how the
partition was cut. Third, `N[count]` is the spec's own endorsed spelling for an unknown-sequence
insert of known size (`complex.md:175`) — a legal HGVS description; the only place ferro refuses an
`N`-unit repeat is the HGVS→SPDI conversion path, where it declines rather than emit a
storable-but-wrong triple, leaving the HGVS spelling untouched.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[duplication-must-ranks-the-label-not-the-partition](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The rule that a duplication must be labelled 'dup' ranks the label of each piece ferro derives, not the partition; the one exception is a net-longer tandem copy of a multi-base motif, where the derivation is cut to expose the dup rather than merged into a delins.
>
> **[spdi-n-unit-repeat-refusal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — On the HGVS-to-SPDI path ferro refuses an 'N'-unit or 'N'-containing repeat rather than expand it to literal 'N' bases — a project choice, since 'N' states a length, not identified bases, and the recommendations do not reach SPDI conversion.
<!-- why:END:duplication-must-ranks-the-label-not-the-partition,spdi-n-unit-repeat-refusal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.1:c.123+45_123+51TSDinsL1.603bp` | refused | — | the ad-hoc `TSD` / `insL1.603bp` labels — `class="invalid"` at `complex.md:169`; a description must be exact, so this is rejected at parse |
| `NG_012232.1(NM_004006.2):c.123+51_123+52ins[N[603];123+45_123+51]` | conformant | — | the spec's endorsed unknown-insert form (`complex.md:175`): `N[603]` gives the size, the range describes the target-site duplication as an insertion. The `N[count]` HGVS spelling is legal (only the SPDI path refuses an `N`-unit repeat), but on a real `c.` reference ferro would expand the `123+45_123+51` same-reference range citation to literal bases (as the insertion page's same-reference range payloads do), so its output is not the recommended range-citation form — conformant. Foreign `NG_`/`NM_` accessions — parse-only here |

## Coverage note (structural)

The normalization corpus exercises essentially none of `complex.md`: multi-megabase genomic
breakpoints, cross-chromosome and ring `::` composites, `sup`/marker chromosomes, ISCN comma-thousands
renderings and orientation-reversed dup-as-`ins(range)inv` are all outside what the shape-family
generator (`RefShape` / `dump_normalized_corpus`) can emit — the generator has no ISCN/complex axis.
A `0` from that corpus here is a claim about the instrument, not about this file's behaviour. The
spec-fixture census (whose ISCN rows `by_input`, #546, marks expected-rejected) and the hand-authored `ring_segment_wellformedness`,
`issue_1578_followup_self_cancelling_rings` and `a_ring_with_no_telomere_anchor_is_still_accepted`
guards are the only exercise this file's behaviour gets; any future confidence about complex/ISCN
normalization must come from targeted fixtures, not the corpus.

See also → `inversion.md:69` (the primary authority for the derived `ins<range>inv` form),
`inversion.md:19` (`dupinv` marked invalid), `general.md:57` (the self-cancelling prohibition the
ring `::` join escapes).
