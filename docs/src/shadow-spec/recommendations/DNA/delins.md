# Deletion-Insertion — ferro's reading

ferro's reading of the HGVS **deletion-insertion (delins)** recommendations, clause by clause —
each spelling with the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Deletion-Insertion (`r.`)](../RNA/delins.md).*

Executable rows use `NM_004006.3`, the one transcript in the committed slice: `cds_start = 238`
(so `c.-237` is the first transcript base), the codon `c.145_147` is `CGC`, `c.76` is `A`, and
`c.123` is `C`. The spec's own worked examples sit on accessions the slice does not carry, so
those rows are parse-only (`—`).

## `delins.md:5` — definition: not a substitution, not an inversion

> Deletion-Insertion (delins): a sequence change where, compared to a reference sequence, one or more nucleotides are replaced by one or more other nucleotides **and which is not** a substitution or inversion.

Ferro: a delins is one-or-more-for-one-or-more nucleotides that is neither a 1-for-1 substitution
nor a whole-span reverse complement (`inv`, governed on `inversion.md:5`) — whatever remains once
those two are excluded.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.902_909delinsuuu` | recommended | self | illustrative only — a genuine multi-for-multi change is a delins; adjudicated on the RNA page |

## `delins.md:15` — one nucleotide for one is a substitution

> by definition, when **one** nucleotide is replaced by **one** other nucleotide, the change is a [substitution](substitution.md).

Ferro: a delins whose net effect is one nucleotide for one is re-derived as the substitution
form — the resulting sequence is a single-base change, so `>` is used, not `delins` (rule-1, "by
definition").

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76delinsG` | recommended | `NM_004006.3:c.76A>G` | `c.76` is `A`; the net 1-for-1 change collapses to the substitution form |

## `delins.md:16` — two or more consecutive nucleotides are a delins

> changes involving two or more consecutive nucleotides are described as deletion/insertion (delins) variants.

Ferro: two adjacent changed nucleotides (separation zero) merge into one delins — the spec marks
the split spelling `class="invalid"` at separation zero (`substitution.md:32`), the one place it
states an outright merge requirement rather than a preference.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[145C>A;146G>T]` | recommended | `NM_004006.3:c.145_146delinsAT` | two adjacent substitutions (`c.145_146` is `CG`), separation zero — merged to one delins |
| `NM_004006.3:c.145_146delinsAT` | recommended | self | the merged form is a fixed point — payload `AT` does not coincide with reference `CG`, nothing residual to shift |

## `delins.md:17` — the separation rule (individual, not delins)

> two variants separated by one or more nucleotides should be described individually and **not** as a "delins".

Ferro: rule 2, a preference not a ban ("should", not "not allowed"). Changes a nucleotide or more
apart are kept as a cis allele, described individually. Separation is read off the partition
re-derived from the resulting sequence, not the input's member boundaries, so two spellings of one
variant converge on one output.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[separation-is-a-property-of-the-spelling-not-of-the-variant](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Ferro reads the separation between changes off the partition it re-derives from the resulting sequence, not off the input's spelling, so two spellings of one variant converge on one output.
<!-- why:END:separation-rule-force-modal-or-negation,separation-is-a-property-of-the-spelling-not-of-the-variant -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | two changes far more than one nucleotide apart — described individually, kept as a cis allele, not merged |

## `delins.md:18-19` — the one-codon exception

> - **exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a "delins".<br>
>   **NOTE**: this prevents tools predicting the consequences of a variant to make conflicting and incorrect predictions of two different substitutions at one position (e.g., `c.235_237delinsTAT` (`p.Lys79Tyr`) versus `c.[235A>T;237G>T]` (`p.[Lys79*;Lys79Asn]`).<br>

Ferro: the exception merges two changes one nucleotide apart into one delins when together they
alter one amino acid — a fact about the resulting sequence, not the input's spelling, so it applies
whatever the members' edit types. It reaches only an axis that declares a reading frame: a
coding-axis merge under it does not propagate to the derived genomic projection, which stays split
under `delins.md:17`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-codon-carve-out-gap-one](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually.
>
> **[codon-carve-out-shape-restriction](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together alter a single amino acid are ruled one delins whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled — though the coding-codon pass currently implements only the substitution/unchanged-base/substitution triplet.
>
> **[projection-codon-exception-is-decided-by-the-rendered-axis](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The codon merge fires only on an axis that declares a reading frame, so when a coding description merges under it ferro leaves the members individual on the derived genomic axis rather than re-merging it to match.
<!-- why:END:delins-codon-carve-out-gap-one,codon-carve-out-shape-restriction,projection-codon-exception-is-decided-by-the-rendered-axis -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[145C>T;147C>G]` | recommended | `NM_004006.3:c.145_147delinsTGG` | one nucleotide apart, together altering the `c.145_147` codon (`CGC`→`TGG`) — folded to one delins regardless of member edit types |
| `LRG_199t1:c.9002_9009delinsTTT` | recommended | — | the spec's own example; its NOTE names the split `c.[145C>T;147C>G]` as *not* correct because the two together affect one amino acid (parse-only) |

## `delins.md:20-21` — conversions are a delins; `con` is retired

> - **conversions**, a sequence change where a **range of nucleotides** are replaced by a sequence from elsewhere in the genome, are described as a "delins".
>   The previous format "con" is no longer used (see [Community Consultation SVD-WG009](../../consultation/SVD-WG009.md)).

Ferro: SVD-WG009 is accepted — a conversion is a delins whose payload is a coordinate range. A
same-reference numeric-range payload is expanded to literal bases and re-derived (conformant, not
the recommended range form — the same limitation the inverted-duplication payload carries). A
cross-accession bracket payload (a translocation) is preserved verbatim, since ferro cannot fetch
the partner reference to expand it, so that citation stays the recommended form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_000797.3:c.812_829delins908_925` | conformant | — | conversion with a same-reference numeric-range payload; ferro would expand it to literal bases on a real reference, so the output is not the recommended range form (parse-only) |
| `NM_004006.3:c.5_16delins123_134` | conformant | `NM_004006.3:c.[5_8delinsCAGT;10_13delinsACCT;15_16delinsCA]` | executable: ferro expands the `123_134` range to literal bases and re-derives it into a three-member split, so the output is not the recommended range form (same limitation as the inverted-duplication payload) |
| `NM_004006.2:c.812_829delinsN[12]` | recommended | — | unknown-length insert `N[12]` (parse-only) |

## `delins.md:22` — the 3'rule

> for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).

Ferro: canonicalization derives from the resulting sequence and applies the 3'rule as part of that
re-derivation. For a delins the rule bites once re-derivation trims payload bases that coincide
with the reference into residual pieces, which then shift; where two legal spellings compete and
no clause selects, ferro emits what falls out of the resulting sequence rather than preserving the
input's spelling.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro applies the edit set to the reference and re-derives the description from the resulting sequence — subject to every explicit spec tie-break, and breaking any remaining tie toward the already-shipped form — rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32386323delinsGA` | recommended | — | canonical one-position delins example: deletes `g.32386323` (`T`), inserts `GA` (parse-only) |
| `NM_004006.2:c.6775_6777delinsC` | recommended | — | coding example: `GAG` replaced by `C` (parse-only) |
| `LRG_199t1:c.145_147delinsTGG` | recommended | — | `p.Arg49Trp` example (parse-only) |

## `delins.md:26-35` — the deleted sequence is never spelled between `del` and `ins`

> **NOTE**: the recommendation is not to describe the variant as <code class="invalid">NC_000023.11:g.32386323delTinsGA</code>, i.e. describe the deleted nucleotide sequence.

Ferro: the DNA delins grammar admits a `del<seq>ins<seq>` production, so a spelling naming the
deleted bases between `del` and `ins` parses — and, unlike the payload-bearing plain `del`/`dup`
forms (rejected in strict mode, see `deletion.md`), it is accepted in strict mode. Ferro preserves
the explicit deleted sequence for a faithful round-trip and strips it to the short `delins` form on
canonicalization — a legal-but-discouraged form ferro tolerates rather than a grammar rejection.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32386323delTinsGA` | conformant | — | the spec's `class="invalid"` deleted-sequence spelling; ferro accepts it in strict mode and canonicalizes to the short `delins` form (parse-only) |
| `NM_004006.2:c.6775_6777delGAGinsC` | conformant | — | same explicit-deleted spelling on the coding example — accepted, canonicalizes to short form (parse-only) |

## `delins.md:44-47` — when part of the inserted sequence happens to match the reference

> - **`LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`**<br>
>   a deletion of nucleotides `c.850` to `c.901`, replaced by `TTCCTCGATGCCTG`.<br>
>   **NOTE**: parts of the inserted sequence "align" with the reference sequence, giving an alternative description like `c.[850_869del;874_881del;887_897del;901_902insG]`.
>   **The "delins" format is recommended**: it is simpler and prevents software tools making incorrect predictions for the consequences on protein level.

Ferro: this passage decides the merge geometry the sibling pages defer here. When a minimal delins
would split only because some member's payload coincides with the reference bases it replaces,
ferro merges it back to the spanning delins `:47` recommends. The reach is narrow: net-deletion
only (payload shorter than the span), and only where some member is gap-bearing — a split of pure
deletions inserts nothing, so nothing realigned, and it stays governed by `delins.md:17`. A placed
gap needed to align an unequal-length block is an artifact of the alignment, not the "separated
variants" `:17` describes.

On `c.` coding positions this merged form is recommended. In the coding-transcript UTR, as on
`g.`/`m.`/`n.`, it is conformant, not recommended:
`:47`'s stated rationale (avoiding incorrect protein-consequence predictions) has nothing to bite
on where there is no reading frame, so the wider reach rests on `:47`'s axis-neutral simplicity
alone — a maintainer house choice among conformant forms, not a conformance requirement. `r.` is
out of jurisdiction (a DNA document cannot scope it, and `RNA/delins.md` states no `:47`
counterpart).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-merge-vs-individual-gap-two-or-more](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two changes two or more nucleotides apart arise only because part of the inserted sequence coincides with the reference and the block is a net deletion, ferro writes them as one spanning delins rather than individually.
>
> **[delins-recommendation-reach-when-the-input-arrives-split](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Ferro merges a re-derived split into one delins only when some member supplies inserted bases while consuming a different number of reference bases; a split of pure deletions inserts nothing and stays individual.
>
> **[unequal-length-block-a-placed-gap-is-not-a-separation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A lone unequal-length net-deletion delins whose payload merely coincides with the reference, and whose every member other than the placed gap would itself render as a delins, is kept whole on every DNA axis (c./g./m./n., but not r.) rather than split at that coincidence into a separate residual member.
>
> **[delins-payload-coincidence-carve-out-is-coding-dna-scoped](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Where a split exists only because payload bases coincide with the reference, ferro writes it as one spanning delins on every DNA axis (c./g./m./n., but not r.); on the frameless axes this is a disclosed rule-2 deviation and the project's choice among conformant forms.
<!-- why:END:delins-merge-vs-individual-gap-two-or-more,delins-recommendation-reach-when-the-input-arrives-split,unequal-length-block-a-placed-gap-is-not-a-separation,delins-payload-coincidence-carve-out-is-coding-dna-scoped -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` | recommended | — | the spanning delins is the recommended `c.` form, not the `c.[850_869del;…]` split named as the "alternative description" (parse-only) |
| `LRG_199t1:c.[850_869del;874_881del;887_897del;901_902insG]` | recommended | — | the split spelling of the same variant; some member is gap-bearing (`901_902ins`), so ferro merges it to the spanning delins above (parse-only) |

## `delins.md:49-64` — translocations and conversions

> - **`NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]`**<br>

Ferro: these translocation and conversion examples sit on foreign accessions with bracketed
cross-reference or numeric-range payloads (full rules in `complex.md`); all are parse-only. A
cross-accession bracket payload is preserved verbatim (recommended); a same-reference numeric-range
payload would be expanded to literal bases on a real reference (conformant, not the recommended
range form — see `:20-21`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]` | recommended | — | unbalanced-translocation delins, see `complex.md` (parse-only) |
| `NC_000022.10:g.42522624_42522669delins42536337_42536382` | conformant | — | conversion with a same-reference numeric-range payload — would expand to literal bases on a real reference, so not the recommended range form (parse-only) |
| `NC_000012.11:g.6128892_6128954delins[NC_000022.10:g.17179029_17179091]` | recommended | — | conversion with a bracketed cross-reference payload (parse-only) |

## `delins.md:73-77` — a `GC` to `TG` change is not a di-nucleotide substitution

> No, this is not allowed.
> By definition, a substitution changes **one** nucleotide into **one** other nucleotide (see [Substitution](substitution.md)).
> The change `TGT`<code class="del">GC</code>`CA` to `TGT`<code class="ins">TG</code>`CA` should be described as `g.4_5delinsTG`, i.e. a deletion/insertion (delins).

Ferro: the di-nucleotide substitution spelling is refused — the substitution grammar admits one
reference base and one alternate base only. Two adjacent changed nucleotides are described as one
delins, the same separation-zero merge as `delins.md:16`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.145CG>AT` | refused | — | the disallowed di-nucleotide substitution shape (`c.145_146` is `CG`) — rejected at grammar; see `delins.md:16` for the recommended delins form |

## `delins.md:79-84` — the maximum-unchanged-nucleotides Q&A

> Yes, two variants separated by one or more nucleotides should preferably be described individually and not as a "delins" (unless they together affect one amino acid).

Ferro: this restates the `:17` separation rule as a preference ("preferably"), on two grounds —
provenance (the variants may have been reported individually) and an annotation-overlap concern.
Provenance is not recoverable by a normalizer, so ferro re-derives from the resulting sequence
rather than preserving the input's spelling; the separation it measures is read off the re-derived
partition, so two spellings of one variant converge.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-is-a-property-of-the-spelling-not-of-the-variant](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Ferro reads the separation between changes off the partition it re-derives from the resulting sequence, not off the input's spelling, so two spellings of one variant converge on one output.
<!-- why:END:separation-is-a-property-of-the-spelling-not-of-the-variant -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | separated by far more than one nucleotide — described individually under rule 2 (same row as `delins.md:17`) |

## `delins.md:86-89` — the BRCA1 `c.2077` Q&A (separation-zero sub+ins merge)

> The correct description of this variant is `NM_007294.3:c.2077delinsATA`.<br>
> **NOTE**: the answer was modified, i.e. the addition "However, since the variant is likely a combination of two other variants, it is acceptable to describe it as <code class="invalid">NM_007294.3:c.[2077G>A;2077_2078insTA]</code>." was removed.

Ferro: a substitution immediately followed by an insertion, separation zero, coalesces at one
locus — both members consume reference bases there. The committee itself withdrew an earlier
passage permitting the split spelling, correcting toward the merged form. This is separation-zero,
so it says nothing about the ≥2 axis-scope question `:44-47` settles.

The executable twin does not reproduce the BRCA1 delins shape. The spec's `c.2077G>A` is a genuine
substitution, so combined with the flush insertion it is a true delins. On the slice's transcript,
`c.76` (`A`) sits inside an A-run (`c.75_77` is `AAA`), so `c.76A>G` plus a flush `insTA` re-derives
to a pure insertion — the substituted base is recovered by the run, and the minimal edit is
`c.75_76insGT`. That is the correct form for a variant that is, on this reference, purely
insertional; it is not the delins the BRCA1 locus produces.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_007294.3:c.2077delinsATA` | recommended | — | the spec's own answer; a fixed point by construction (parse-only) |
| `NM_007294.3:c.[2077G>A;2077_2078insTA]` | recommended | — | the withdrawn split — the committee removed the passage permitting it; ferro merges it to the delins above (parse-only) |
| `NM_004006.3:c.[76A>G;76_77insTA]` | recommended | `NM_004006.3:c.75_76insGT` | executable twin of the sub+ins shape: `c.76` sits in the `c.75_77` A-run, so the sub+ins re-derives to a pure insertion — the minimal edit `c.75_76insGT`, not the BRCA1 delins |
| `NM_004006.3:c.76delinsGTA` | recommended | `NM_004006.3:c.75_76insGT` | the delins of `GTA` at `c.76` (`A`, in the A-run) is sequence-identical to inserting `GT`; ferro collapses it to the minimal insertion form |

## `delins.md:68-71` — the "indel" terminology Q&A

> The term "indel" is not used in HGVS nomenclature (see [Glossary](../../background/glossary.md)).

Descriptive — the term "indel" is not used in HGVS nomenclature (it is ambiguous across
disciplines).
