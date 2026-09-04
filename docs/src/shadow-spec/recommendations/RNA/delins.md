# Deletion-insertion — ferro's reading

ferro's reading of the HGVS **delins** recommendations on the transcript (`r.`) axis, clause by
clause — each spelling with the form ferro normalizes it to and a verdict on that output. New
here? See [How to read a page](../../reading-guide.md) for the verdicts, the table conventions,
and the recurring terms.

*DNA twin: [Deletion-insertion (`c.`/`g.`)](../DNA/delins.md).*

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's worked
examples at `delins.md:29-42` carry **no accession**, and their stated flanks are not
`NM_004006.3`'s (`r.775` is `a` here, not `u`; `r.142_144` is `agg`, not `cga`) — so each of those
examples is shown on a stated-flank twin: a locus on `NM_004006.3` whose bases have the same shape
the spec describes. Rows on accessions the slice does not carry (`NM_004006.2`, `NM_007294.3`, the
fusion partners) are parse-only (`—`).

Two rows below carry `refused` for a reason the legend does not cover: the spec's bracketed
cross-reference payloads (`delins[ACC:g.…]`, `delins[ACC:r.…]`) are forms ferro's parser does
**not yet accept**. That is a ferro parser gap on a spec-recommended form, not an invalid input;
each such row's Note says so. It is recorded as `refused` only because that is the mechanically
true statement about ferro today.

## `delins.md:5` — definition: not a substitution, not an inversion

> Deletion-Insertion (delins): a sequence change where, compared to a reference sequence, one or
> more nucleotides are replaced by one or more other nucleotides **and which is not** a
> substitution or inversion.

Ferro: a delins is one-or-more-for-one-or-more nucleotides that is neither a 1-for-1 substitution
nor a whole-span reverse complement (typed `inv`). The inversion half of this boundary on `r.`
would rest on `RNA/inversion.md`, which this page does not adjudicate — the DNA-side ruling
(`whole-span-reverse-complement-types-as-inv`) has not been extended to this axis; open.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.902_909delinsuuu` | recommended | self | a genuine delins (`acacacag` → `uuu`) — see `:35-37` |
| `NM_004006.3:r.5delinsa` | recommended | `NM_004006.3:r.5u>a` | 1-for-1 net effect is a substitution, not a delins — see `:16` |

## `delins.md:16` — one nucleotide for one is a substitution, not a delins

> by definition, when **one** nucleotide is replaced by **one** other nucleotide, the change is a
> [substitution](substitution.md).

Ferro: a `delins` whose net effect is exactly one nucleotide for one is normalized to the
substitution form, on rule-1 ("by definition") force.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.5delinsa` | recommended | `NM_004006.3:r.5u>a` | `r.5` is `u`; the net 1-for-1 change collapses to the recommended substitution |

## `delins.md:17-19` — separation, and the one-nucleotide codon exception

> - two variants separated by one or more nucleotides should be described individually and
>   **not** as a "delins".
>     - **exception**: two variants separated by one nucleotide, together affecting one amino
>       acid, should be described as a "delins" (e.g., `r.142_144delinsugg` `p.(Arg48Trp)`).<br>
>       **NOTE**: this prevents tools predicting the consequences of a variant to make
>       conflicting and incorrect predictions of two different substitutions at one position.<br>

Ferro: rule 2 stands unqualified on `r.` — members a nucleotide or more apart are described
individually — except for the one-nucleotide codon exception, which reaches `r.` on this
document's own authority. The DNA-side payload-coincidence widening (`DNA/delins.md:47`) has no
`RNA/delins.md` counterpart, so it does not reach this axis at all.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[delins-codon-carve-out-gap-one](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually; the merged span must stay within the reading frame's codon boundary regardless of the members' edit types.
>
> **[projection-codon-exception-is-decided-by-the-rendered-axis](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The codon merge fires only on an axis that declares a reading frame, so when a coding description merges under it ferro leaves the members individual on the derived genomic axis rather than re-merging it to match.
>
> **[delins-payload-coincidence-carve-out-is-coding-dna-scoped](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Where a split exists only because payload bases coincide with the reference, ferro writes it as one spanning delins on every DNA axis (c./g./m./n., but not r.); on the frameless axes this is a disclosed rule-2 deviation and the project's choice among conformant forms.
<!-- why:END:separation-rule-force-modal-or-negation,delins-codon-carve-out-gap-one,projection-codon-exception-is-decided-by-the-rendered-axis,delins-payload-coincidence-carve-out-is-coding-dna-scoped -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[100a>u;104a>g]` | recommended | self | separated by 3 nt and in different codons (34 and 35) — no codon-exception antecedent; stays split, `r.`'s recommended form under rule 2 |
| `NM_004006.3:r.[142a>u;144g>c]` | recommended | `NM_004006.3:r.142_144delinsugc` | codon 48 (`agg` → `ugc`, `p.(Arg48Cys)`); one nucleotide apart, one amino acid — the `:18` exception reaches `r.` on this document's own authority. Worked in full at `:38-42` below |

## `delins.md:20-21` — conversions are a delins; `con` is retired

> **conversions**, a sequence change where a **range of nucleotides** are replaced by a sequence
> from elsewhere in the genome, are described as a "delins".
> The previous format "con" is no longer used (see [Community Consultation
> SVD-WG009](../../consultation/SVD-WG009.md)).

Ferro: SVD-WG009 is accepted; the retired `con` syntax is never emitted — refused in strict,
repaired to `delins` in lenient.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:r.2623_2803delins2804_2949` | conformant | — | the delins form for a conversion (SVD-WG009 accepted); on a real reference ferro would expand the same-accession range to literal bases (see the executable probe below), so not the recommended range form — conformant (the spec's own accession, not in the slice — parse-only here) |
| `NM_004006.3:r.5_16delins123_134` | conformant | `NM_004006.3:r.[5_8delinscagu;10_13delinsaccu;15_16delinsca]` | executable conversion probe on the slice: a same-accession numeric-range payload — ferro expands the `123_134` range to literal bases and re-derives the resulting change into a three-member split, so not the recommended range form — conformant |
| `NM_004006.2:r.2623_2803con2804_2949` | refused | — | the retired `con` syntax — rejected in strict mode |

## `delins.md:22-24` — adjoined transcripts: a special-case delins

> Adjoined transcripts from gene fusions represent a special case of deletion-insertion variant.
> The fusion break point is described using **`::`**.<br>
> **NOTE**: to avoid confusion, HGVS recommends to follow the [HGNC
> guidelines](https://www.genenames.org/about/guidelines/) to describe products of gene
> translocations or fusions (format `GENESYMBOL1::GENESYMBOL2`) and readthrough transcripts
> (format `GENESYMBOL1-GENESYMBOL2`).

Ferro: descriptive — this line delegates syntax entirely to `adjoined_transcript.md`. "Special
case of delins" is not read as "may be rewritten as a delins": with two reference sequences,
nothing is re-derivable and nothing 3'-shifts, so this is not adjudicated here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | — | the spec's own `TPM3::PDGFRB` example (`adjoined_transcript.md:29`); neither partner is in the slice — parse-only. The full rules live in `adjoined_transcript.md`, not on this page |
| `NM_152263.2:c.-115_775::NM_002609.3:c.1580_*1924` | refused | — | `c.` partners inside `::` are rejected — the syntax is RNA-only, per `adjoined_transcript.md:18` |

## `delins.md:25` — the 3' rule

> for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily
> assigned to have been changed (**3'rule**).

Ferro: canonicalization derives from the resulting sequence and applies the 3' rule as part of
that re-derivation; on `r.` the exon-junction exception to it lives at `general.md:43`. For a
delins the rule bites only after re-derivation trims payload bases that coincide with the
reference into residual pieces, which then shift.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.775_776delinsa` | recommended | `NM_004006.3:r.777del` | `r.775_777` is `aaa`; the net effect is one `a` deleted from the run, and the 3' rule places it on `r.777` |
| `NM_004006.3:r.149delinsga` | recommended | self | `c[u]c` → `c[ga]c`: neither payload base coincides with the reference, nothing residual to shift — a fixed point |

## `delins.md:29-31` — example: `r.775delinsga`

> **`r.775delinsga`**<br>
> a deletion of nucleotide `r.775` (a `u`, not described), replaced by nucleotides `ga`, changing
> `..aggc`<code class="del">u</code>`cauu..` to `..aggc`<code class="ins">ga</code>`cauu..`.

Ferro: on the stated flank neither payload base equals the deleted `u`; both minimal two-member
partitions are separation zero with an `ins` on one side, so both merge to one delins under the
`ins` half of `delins-adjacent-members-when-both-consume-reference` — a fixed point. The stated
flank is not `NM_004006.3`'s (`r.775` is `a` here), so the example is shown on its stated-flank
twin `r.149` (`c[u]c`).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.149delinsga` | recommended | self | stated-flank twin (`cuc` → `cgac`); fixed point |
| `NM_004006.3:r.[149u>g;149_150insa]` | recommended | `NM_004006.3:r.149delinsga` | one minimal split, separation zero — converges to the delins |
| `NM_004006.3:r.[148_149insg;149u>a]` | recommended | `NM_004006.3:r.149delinsga` | the other minimal split — also converges |
| `NM_004006.3:r.775delinsga` | recommended | `NM_004006.3:r.774_775insg` | the spec's literal spelling on this transcript's bases: `r.775` is `a`, so the payload's `a` coincides with the reference and the variant is an insertion of `g`, not a delins — ferro re-derives from the resulting sequence rather than preserving the spelling |

**Confluence.** `{r.[149u>g;149_150insa], r.[148_149insg;149u>a], r.149delinsga} → r.149delinsga`.

## `delins.md:32-34` — example: `r.775_777delinsc`, and the `r.`-exclusion

> **`r.775_777delinsc`**<br>
> a deletion of nucleotides `r.775` to `r.777` (`uca`, not described), replaced by nucleotide `c`,
> changing `..aggc`<code class="del">uca</code>`uu..` to `..aggc`<code class="ins">c</code>`uu..`.

Ferro: on the stated flank the payload `c` coincides with reference base `776`, so the minimal
partition is `[775del;777del]` — one unchanged base apart, `:17`'s antecedent. On a DNA axis this
exact shape now re-merges under the payload-coincidence carve-out (`DNA/delins.md:47`); that
carve-out is scoped to the DNA axes and states no `r.` counterpart, so on a non-coding `r.`
transcript rule 2 governs unqualified and the block **stays split** — a rule-2-conformant `r.`
split, not a ferro limitation. On a coding transcript, whether the one-nucleotide codon exception
(`:18`) separately licenses keeping this specific block whole is a genuinely open question (see
Notes). The stated flank is not `NM_004006.3`'s (`r.775_777` is `aaa` here, so the literal
spelling has no coincidence at all); the executable rows use stated-flank twins.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-payload-coincidence-carve-out-is-coding-dna-scoped](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Where a split exists only because payload bases coincide with the reference, ferro writes it as one spanning delins on every DNA axis (c./g./m./n., but not r.); on the frameless axes this is a disclosed rule-2 deviation and the project's choice among conformant forms.
>
> **[unequal-length-block-a-placed-gap-is-not-a-separation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A lone unequal-length net-deletion delins whose payload merely coincides with the reference, with no higher-priority competing type, is kept whole on every DNA axis (c./g./m./n., but not r.) rather than split at that coincidence into a separate residual member, since a placed alignment gap is not itself a fact about the variant for a separation rule to key on.
>
> **[codon-carve-out-shape-restriction](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together alter a single amino acid are written as one delins, whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled.
<!-- why:END:delins-payload-coincidence-carve-out-is-coding-dna-scoped,unequal-length-block-a-placed-gap-is-not-a-separation,codon-carve-out-shape-restriction -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NR_x.1:r.775_777delinsc` | recommended | — | non-coding twin, illustrative accession (not resolvable — parse-only). The decided reading: payload coincides with reference `776c`, separation one, no reading frame for `:18`, `:17` unqualified — `r.[775del;777del]` is `r.`'s recommended form, since the DNA-scoped carve-out does not reach this axis (`delins-payload-coincidence-carve-out-is-coding-dna-scoped`, `unequal-length-block-a-placed-gap-is-not-a-separation`) |
| `NM_004006.3:r.123_125delinsa` | recommended | `NM_004006.3:r.[123del;125del]` | coding transcript, stated-flank twin **spanning** codons 41/42 (`u[cag]u` → `u[a]u`): payload coincides with `124a`, separation one, but the pair does not sit within one codon — `:18` has no antecedent, `:17` governs, split |
| `NM_004006.3:r.103_105delinsa` | conformant | self | coding transcript, stated-flank twin **within** codon 35 (`g[cag]c` → `g[a]c`): same shape, one codon; whether `:18`'s "together affecting one amino acid" also covers this net −2 (frameshift) pair is open — provisional, see Notes |

**Notes.** The spec publishes this example whole, and it is accession-less. On a non-coding `r.`
reading, no clause licenses re-merging the split — that is the locked reading on this page, and it
holds even though the sibling DNA axes now merge the identical shape under the payload-coincidence
carve-out. On a coding `r.` reading, the codon exception's antecedent is stated without regard to
edit type or net length, and whether a net-length-changing (frameshift) pair counts as "together
affecting one amino acid" is not yet decided in the ledger (`codon-carve-out-shape-restriction`
leaves it open; ferro's implementation merges on span-within-one-codon alone). Until it is
decided, ferro's whole-block output on a coding transcript is marked `conformant` rather than
`recommended`; this is a provisional verdict pending that ledger decision, not a known ferro bug,
and no issue is filed for it.

## `delins.md:35-37` — example: `r.902_909delinsuuu`

> **`r.902_909delinsuuu`**<br>
> a deletion of nucleotides `r.902` to `r.909`, replaced by nucleotides `uuu`.

Ferro: no flank is stated; on `NM_004006.3` (`r.902_909` is `acacacag`) no payload base coincides
with any reference base in the deleted span, and the span crosses codons 301-303. Nothing to
split, nothing to shift — fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.902_909delinsuuu` | recommended | self | no payload/reference coincidence; span crosses three codons — fixed point |

## `delins.md:38-42` — example: `r.142_144delinsugg`, and the polymorphism note

> **`r.142_144delinsugg` (`p.Arg48Trp`)**<br>
> a deletion replacing nucleotides `r.142` to `r.144` (`cga`, not described) with `ugg`.<br>
> **NOTE**: the variant can also be described as `r.[142c>u;144a>g]`, i.e. two substitutions.
> This format is preferred when either of the two variants is known as a frequently occurring
> variant ("polymorphism").

Ferro: on the stated flank, separation is one (`143g` unchanged) and codon 48 — the `:18`
exception merges the split to one delins. The spec's own preference for the split, on population
frequency, is provenance no normalizer holds; ferro re-derives from the resulting sequence
instead. `NM_004006.3`'s codon 48 is `agg`, not the stated `cga` (the literal `r.142_144delinsugg`
would re-derive to the single substitution `r.142a>u` here), so the rows use the same-shaped
`agg` → `ugc`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-codon-carve-out-gap-one](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually; the merged span must stay within the reading frame's codon boundary regardless of the members' edit types.
>
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:delins-codon-carve-out-gap-one,canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.142_144delinsugc` | recommended | self | codon 48; `p.(Arg48Cys)` — fixed point |
| `NM_004006.3:r.[142a>u;144g>c]` | recommended | `NM_004006.3:r.142_144delinsugc` | the polymorphism-preferred split shape; ferro re-derives from the resulting sequence rather than preserving that provenance |

**Confluence.** `{r.[142a>u;144g>c], r.142_144delinsugc} → r.142_144delinsugc`.

## `delins.md:43-51` — RNA conversion examples

> - RNA conversion (based on [SVD-WG009](../../consultation/SVD-WG009.md))
>     - **`NM_004006.2:r.2623_2803delins2804_2949`**<br>
>       conversion replacing nucleotides `r.2623` to `r.2803` (exon 21) with nucleotides `r.2804`
>       to `r.2949` (exon 22) as found in the _DMD_ coding RNA sequence file `NM_004006.2`.
>     - **`r.415_1655delins[AC096506.5:g.409_1649]`**<br>
>       conversion replacing nucleotides `r.415` to `r.1655` with nucleotides `g.409` to `g.1649`
>       as found in the genomic reference sequence `AC096506.5`.
>     - **`r.1401_1446delins[NR_002570.3:r.1513_1558]`**<br>
>       conversion in exon 9 of the _CYP2D6_ gene, replacing exon 9 nucleotides `r.1401` to
>       `r.1446` with those of the 3' flanking _CYP2D7P1_ gene, nucleotides `r.1513` to `r.1558`.

Ferro: three payload forms — a same-accession range, a bracketed foreign genomic range, a
bracketed foreign RNA range. Ferro parses only the same-accession range form, and does **not** keep
it as a range: on a real reference it expands the range to literal bases and normalizes the result
(the executable probe at `:20-21` shows this on the slice), so its output is not the recommended
`delins<range>` citation form — `conformant`. The two bracketed cross-reference forms it does **not
yet parse** at all — a parser gap, recorded as `refused` below.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:r.2623_2803delins2804_2949` | conformant | — | same-accession conversion range; ferro would expand it to literal bases on a real reference (see the `:20-21` probe), not preserve it — so not the recommended range form (the spec's own accession, not in the slice — parse-only) |
| `r.415_1655delins[AC096506.5:g.409_1649]` | refused | — | the spec's own foreign-genomic-range payload. **Ferro parser gap**: `parse_hgvs` does not yet accept a bracketed cross-reference payload, so this spec-recommended form cannot be built or normalized today; `refused` records that fact, not a judgment that the input is invalid |
| `r.1401_1446delins[NR_002570.3:r.1513_1558]` | refused | — | the spec's own foreign-RNA-range payload (_CYP2D6_/_CYP2D7P1_). Same ferro parser gap as the row above — a limitation on a spec-recommended form, not an invalid input |

**Notes.** `r.2623_2803delins2804_2949` is sequence-identical to
`r.[2623_2803del;2804_2949dup]` — a `del` adjacent to a `dup` at separation zero. The `dup` half of
the separation-zero merge is still open (`delins-adjacent-members-when-both-consume-reference`
decides only the `del`/`ins`/`sub`/`inv` sides), so this equivalence class is not yet in the
ledger — open, provisional.

## `delins.md:62-66` — Q&A: `gc>ug` is not a di-nucleotide substitution

> !!! note "Can I describe a `gc` to `ug` variant as a di-nucleotide substitution
> (<code class="invalid">r.4gc>ug</code>)?"
>
> No, this is not allowed.
> By definition, a substitution changes **one** nucleotide into **one** other nucleotide (see
> [Substitution](substitution.md)). The change `augu`<code class="del">gc</code>`ca` to
> `augu`<code class="ins">ug</code>`ca` should be described as `r.5_6delinsug`, i.e. a
> deletion/insertion (delins).

Ferro: the di-nucleotide substitution spelling is refused; two adjacent changed nucleotides
(separation zero) are described as one delins — RNA's own authority for the sub+sub merge.
`NM_004006.3:r.5_6` is `uu` (not the Q&A's `gc`), so the executable rows change both to `ag`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.4gc>ug` | refused | — | the disallowed di-nucleotide substitution spelling |
| `NM_004006.3:r.[5u>a;6u>g]` | recommended | `NM_004006.3:r.5_6delinsag` | adjacent, separation zero — merges to the recommended delins |
| `NM_004006.3:r.5_6delinsag` | recommended | self | the merged form is a fixed point |

**Confluence.** `{r.[5u>a;6u>g], r.5_6delinsag} → r.5_6delinsag`.

**Note.** The governing record does not yet cite this clause or its `RNA/substitution.md:16`
twin directly — an open ledger citation gap, not a behavior question.

## `delins.md:68-71` — Q&A: the BRCA1 sub+ins merge

> !!! note "The _BRCA1_ coding RNA reference sequence `NM_007294.3` from position `r.2074` to `r.2080` is
> `cau`<code class="sub">g</code>`aca`. A variant frequently found in the population is
> `cau`<code class="sub">a</code>`aca` (`NM_007294.3:r.2077g>a`). In a patient I found the sequence
> `cau`<code class="sub">a</code><code class="ins">ua</code>`aca`. Can I describe this variant as
> <code class="invalid">NM_007294.3:r.[2077g>a;2077_2078insua]</code>?"
>
> The correct description of this variant is `NM_007294.3:r.2077delinsaua`.<br>
> **NOTE**: the answer was modified, i.e. the addition "However, since the variant is likely a
> combination of two other variants, it is acceptable to describe it as
> <code class="invalid">NM_007294.3:r.[2077g>a;2077_2078insua]</code>." was removed.

Ferro: separation-zero sub+ins merges to one delins; the committee withdrew a provenance-based
permission for the split — the strongest register this corpus has for this shape. The spec's own
`NM_007294.3` rows are parse-only here (not in the slice; the stated-flank equivalence is pinned
on a synthetic transcript in `tests/it/spec_worked_example_rules.rs`); the same shape is executed
on `NM_004006.3` at `r.149` (`c[u]c` → `c[aga]c`).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_007294.3:r.2077delinsaua` | recommended | — | the spec's own answer; fixed point by construction (parse-only here) |
| `NM_007294.3:r.[2077g>a;2077_2078insua]` | recommended | — | the withdrawn split — the committee removed the passage permitting it (parse-only here) |
| `NM_004006.3:r.149delinsaga` | recommended | self | executable twin of the sub+ins shape; fixed point |
| `NM_004006.3:r.[149u>a;149_150insga]` | recommended | `NM_004006.3:r.149delinsaga` | executable twin of the withdrawn split: separation zero, sub + ins — merges to the delins |

**Confluence.** `{r.[149u>a;149_150insga], r.149delinsaga} → r.149delinsaga`.
