# Deletion-Insertion — ferro's reading

ferro's reading of the HGVS **deletion-insertion (delins)** recommendations, clause by clause —
each spelling with the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Deletion-Insertion (`r.`)](../RNA/delins.md).*

Most of this page is CONFIRM-by-inspection against the spec text and the shipped code: the delins
definition, the substitution boundary, the conversions housekeeping, and the 3'rule are mechanical
parser or re-derivation behaviour with no clause in tension. The units that are actually
adjudicated — and carry a Why block — are the separation rule at `delins.md:17`, the one-codon
exception at `delins.md:18`, the two separation-zero merges (`delins.md:16` sub+sub and
`delins.md:86-89` sub+ins), and the central example at `delins.md:44-47` of when part of the
inserted sequence happens to match the reference, which is where the merge geometry that the
sibling `deletion.md` and `substitution.md` pages defer to this file is actually decided.

**The payload-coincidence carve-out reaches every DNA axis (`c.`, `g.`, `m.`, `n.`), not `c.` alone.**
On `c.` the merged spanning delins is the **recommended** form; on `g.`/`m.`/`n.` — the frameless
axes, where `delins.md:47`'s protein-consequence rationale has nothing to bite on — it is a
**conformant** maintainer house choice made on `:47`'s axis-neutral simplicity ground (README rule 6),
not a conformance requirement. The `r.` axis is out of jurisdiction entirely (`RNA/delins.md`
states no `:47` counterpart) and is out of scope for this DNA page.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The CDS base facts
the rows rely on are the ones the sibling pages establish: `cds_start = 238` (so `c.-237` is the
first transcript base), the codon `c.145_147` is `CGC`, `c.76` is `A`, `c.123` is `C`, and the
8-nucleotide A-stretch is `c.5690_5697` (`c.5698` is `T`). The spec's own worked examples sit on
`NM_004006.2` (a different transcript version), `LRG_199t1`, and foreign genomic accessions
(`NC_000023.11`, `NC_000002.12`, `NC_000022.10`, `NC_000012.11`, `NC_000011.10`, `NM_000797.3`,
`NM_007294.3`), none of which the slice carries, so those rows are parse-only (`—`) — ferro cannot
read their bases here. Executable `NM_004006.3` rows carry ferro's actual normalized output, as the
bless harness produced it against the committed slice.

## `delins.md:5` — definition: not a substitution, not an inversion

> Deletion-Insertion (delins): a sequence change where, compared to a reference sequence, one or more nucleotides are replaced by one or more other nucleotides **and which is not** a substitution or inversion.

Ferro: a delins is one-or-more-for-one-or-more nucleotides that is neither a 1-for-1 substitution
nor a whole-span reverse complement (which types as `inv`, `whole-span-reverse-complement-types-as-inv`,
governed on `inversion.md:5` and out of this page's scope). The delins label is what remains after
substitution and inversion are excluded, so it never competes with either for the same span.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.902_909delinsuuu` | recommended | self | shown only to name the boundary — a genuine multi-for-multi change is a delins. `r.` is adjudicated on the RNA page; this row is illustrative |

## `delins.md:15` — one nucleotide for one is a substitution

> by definition, when **one** nucleotide is replaced by **one** other nucleotide, the change is a [substitution](substitution.md).

Ferro: a `delins` whose net effect is exactly one nucleotide for one is re-derived to the
substitution form, on rule-1 ("by definition") force — the resulting sequence is a single-base
change, so the recommended spelling is `>`, not `delins`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76delinsG` | recommended | `NM_004006.3:c.76A>G` | `c.76` is `A`; the net 1-for-1 change collapses to the substitution form |

## `delins.md:16` — two or more consecutive nucleotides are a delins

> changes involving two or more consecutive nucleotides are described as deletion/insertion (delins) variants.

Ferro: two adjacent changed nucleotides (separation **zero**) are one delins — this is the one place
the spec states an outright merge requirement rather than a preference: at separation zero the split
spelling is marked `class="invalid"` and called "not correct" by name (`substitution.md:32`). Both
members consume reference bases, so they coalesce.

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

Ferro: ruleset rule 2 (a preference, not a ban) — the "and not" names the excluded alternative, the
modal "should" grades the whole clause. Two changes a nucleotide or more apart are kept as a cis
allele, described individually. The separation that this clause keys on is read off the partition
**re-derived from the resulting sequence**, never off the input's member boundaries, so two spellings
of one variant converge on one separation count and one output.

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

Ferro: the exception folds two changes into one delins when they sit exactly one nucleotide apart
and together change one amino acid — applied triplet-precisely, so the merged span stays within the
reading frame, whatever the members' edit types ("together affecting one amino acid" is a fact about
the resulting sequence, not about how the input was spelled). Because the antecedent needs an
amino-acid consequence, the exception reaches only an axis that declares a reading frame: a
coding-axis merge under it does **not** propagate to its derived genomic projection, which stays
split under `delins.md:17`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-codon-carve-out-gap-one](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually; the merged span must stay within the reading frame's codon boundary regardless of the members' edit types.
>
> **[codon-carve-out-shape-restriction](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together alter a single amino acid are written as one delins, whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled.
>
> **[projection-codon-exception-is-decided-by-the-rendered-axis](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The codon merge fires only on an axis that declares a reading frame, so when a coding description merges under it ferro leaves the members individual on the derived genomic axis rather than re-merging it to match.
<!-- why:END:delins-codon-carve-out-gap-one,codon-carve-out-shape-restriction,projection-codon-exception-is-decided-by-the-rendered-axis -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[145C>T;147C>G]` | recommended | `NM_004006.3:c.145_147delinsTGG` | one nucleotide apart, together altering the `c.145_147` codon (`CGC`→`TGG`) — folded to one delins, whatever the member edit types (the sibling `substitution.md` pins this same row) |
| `LRG_199t1:c.9002_9009delinsTTT` | recommended | — | the spec's own example whose NOTE names the split `c.[145C>T;147C>G]` as *not* correct because the two together affect one amino acid (parse-only — `LRG_199t1` not in the slice) |

## `delins.md:20-21` — conversions are a delins; `con` is retired

> - **conversions**, a sequence change where a **range of nucleotides** are replaced by a sequence from elsewhere in the genome, are described as a "delins".
>   The previous format "con" is no longer used (see [Community Consultation SVD-WG009](../../consultation/SVD-WG009.md)).

Ferro: SVD-WG009 is accepted — a conversion is a delins whose payload is a coordinate range. Ferro
does **not** keep a **same-reference** numeric-range payload: on a real reference it expands the range
to literal bases and normalizes the result (the executable probe below shows `c.5_16delins123_134`
re-derived to a three-member split), so its output is not the recommended `delins<range>` citation
form — `conformant`, the same limitation the inverted-duplication payload carries. A
**cross-accession** bracket payload (`delins[OTHER:g...]`, a translocation) is preserved instead,
because ferro cannot fetch the partner reference to expand it, so that citation *is* the recommended
form it keeps — `recommended`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_000797.3:c.812_829delins908_925` | conformant | — | conversion with a same-reference numeric-range payload; on a real reference ferro would expand it to literal bases (see the executable row below), so its output is not the recommended range form — conformant (parse-only — `NM_000797.3` not in the slice) |
| `NM_004006.3:c.5_16delins123_134` | conformant | `NM_004006.3:c.[5_8delinsCAGT;10_13delinsACCT;15_16delinsCA]` | executable conversion with a same-reference numeric-range payload: ferro expands the `123_134` range to its literal bases and re-derives the resulting change into a three-member split — so the output is not the recommended `delins<range>` citation form (a limitation, like the inverted-duplication payload) — conformant |
| `NM_004006.2:c.812_829delinsN[12]` | recommended | — | unknown-length insert `N[12]` (parse-only — `NM_004006.2` is a different transcript version, not in the slice) |

## `delins.md:22` — the 3'rule

> for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).

Ferro: canonicalization derives from the resulting sequence and applies the 3'rule as part of that
re-derivation. For a delins the rule bites only after re-derivation trims any payload bases that
coincide with the reference into residual pieces, which then shift; where two legal spellings compete
and no clause selects, ferro emits what falls out of the resulting sequence rather than preserving the
input's spelling.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32386323delinsGA` | recommended | — | the spec's canonical one-position delins example, deleting `g.32386323` (a `T`) and inserting `GA` (parse-only — foreign genomic accession) |
| `NM_004006.2:c.6775_6777delinsC` | recommended | — | the spec's own coding example, `GAG` replaced by `C` (parse-only — `NM_004006.2` version not in the slice) |
| `LRG_199t1:c.145_147delinsTGG` | recommended | — | the spec's `p.Arg49Trp` example (parse-only — `LRG_199t1` not in the slice) |

## `delins.md:26-35` — the deleted sequence is never spelled between `del` and `ins`

> **NOTE**: the recommendation is not to describe the variant as <code class="invalid">NC_000023.11:g.32386323delTinsGA</code>, i.e. describe the deleted nucleotide sequence.

Ferro: the DNA delins grammar **does** admit a `del<seq>ins<seq>` production
(`parse_delins_with_deleted_seq`), so a spelling that names the deleted bases between `del` and
`ins` parses — and, unlike the payload-bearing plain `del`/`dup` forms (which strict rejects with
`W3025`, see the sibling `deletion.md`), it is **accepted in strict mode**. Ferro preserves the
explicit deleted sequence for a faithful round-trip and strips it back to the recommended short
`delins` form on canonicalization. So this is a legal-but-discouraged form ferro tolerates, not a
grammar rejection; the `class="invalid"` markup is a spec stylistic preference ferro does not
currently enforce at parse. No Why block — this is parser/strict behaviour, not an adjudication.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32386323delTinsGA` | conformant | — | the spec's `class="invalid"` deleted-sequence spelling; ferro accepts it (strict does **not** reject the explicit-deleted `delins`) and canonicalizes to the short `delins` form (parse-only — foreign genomic) |
| `NM_004006.2:c.6775_6777delGAGinsC` | conformant | — | the same explicit-deleted spelling on the coding example — accepted, canonicalizes to short form (parse-only — `NM_004006.2` version not in the slice) |

## `delins.md:44-47` — when part of the inserted sequence happens to match the reference

> - **`LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`**<br>
>   a deletion of nucleotides `c.850` to `c.901`, replaced by `TTCCTCGATGCCTG`.<br>
>   **NOTE**: parts of the inserted sequence "align" with the reference sequence, giving an alternative description like `c.[850_869del;874_881del;887_897del;901_902insG]`.
>   **The "delins" format is recommended**: it is simpler and prevents software tools making incorrect predictions for the consequences on protein level.

Ferro: this is the passage that decides the merge geometry the sibling pages defer here. When a
minimal delins would be split **only** because some derived member's payload coincides with the
reference bases it replaces — a gap-bearing member, an inserted sequence that genuinely realigned —
ferro merges it back to the spanning delins that `:47` recommends. The reach is narrow and directional:
**net-deletion only** (payload shorter than the span), and **only** where some member is gap-bearing;
a split of pure deletions inserts nothing, so nothing realigned, and `general.md:33`/`delins.md:17`
govern it unqualified. A placed gap needed to align an unequal-length block is an artifact of the
chosen alignment, not the "separated variants" `:17` was written to describe.

**On `c.` this merged spanning delins is the recommended form.** On `g.`/`m.`/`n.` it is
**conformant, not recommended**: `:47`'s stated rationale — preventing incorrect protein-consequence
predictions — has nothing to bite on where there is no reading frame, so the widened scope rests on
`:47`'s axis-neutral *simplicity* ground alone, which makes it a README rule-6 maintainer house
choice among conformant forms rather than a conformance requirement. The `r.` axis is out of
jurisdiction (a DNA document cannot scope `r.`, and `RNA/delins.md` states no `:47` counterpart).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-merge-vs-individual-gap-two-or-more](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two changes two or more nucleotides apart arise only because part of the inserted sequence coincides with the reference and the block is a net deletion, ferro writes them as one spanning delins rather than individually.
>
> **[delins-recommendation-reach-when-the-input-arrives-split](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Ferro merges a re-derived split into one delins only when some member supplies inserted bases while consuming a different number of reference bases; a split of pure deletions inserts nothing and stays individual.
>
> **[unequal-length-block-a-placed-gap-is-not-a-separation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A lone unequal-length net-deletion delins whose payload merely coincides with the reference, with no higher-priority competing type, is kept whole on every DNA axis (c./g./m./n., but not r.) rather than split at that coincidence into a separate residual member, since a placed alignment gap is not itself a fact about the variant for a separation rule to key on.
>
> **[delins-payload-coincidence-carve-out-is-coding-dna-scoped](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Where a split exists only because payload bases coincide with the reference, ferro writes it as one spanning delins on every DNA axis (c./g./m./n., but not r.); on the frameless axes this is a disclosed rule-2 deviation and the project's choice among conformant forms.
<!-- why:END:delins-merge-vs-individual-gap-two-or-more,delins-recommendation-reach-when-the-input-arrives-split,unequal-length-block-a-placed-gap-is-not-a-separation,delins-payload-coincidence-carve-out-is-coding-dna-scoped -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.850_901delinsTTCCTCGATGCCTG` | recommended | — | the spec's own worked example — the spanning delins is the recommended `c.` form, not the `c.[850_869del;…]` split it names as the "alternative description" (parse-only — `LRG_199t1` not in the slice) |
| `LRG_199t1:c.[850_869del;874_881del;887_897del;901_902insG]` | recommended | — | the split spelling of the same variant; on `c.`, some member is gap-bearing (`901_902ins`), so ferro merges it to the spanning delins above (parse-only — the merge target is the recommended `c.` form) |

## `delins.md:49-64` — translocations and conversions

> - **`NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]`**<br>

Ferro: the translocation and conversion examples are all foreign genomic (or foreign-transcript)
accessions with bracketed cross-reference or numeric-range payloads; the full rules for translocations
live in `complex.md`. None is in the slice, so all are parse-only. A **cross-accession** bracket
payload is preserved verbatim (ferro cannot fetch the partner reference to expand it), so it is
`recommended`; a **same-reference** numeric-range payload would be expanded to literal bases on a real
reference (as `:20-21`'s probe shows), so it is `conformant`, not the recommended range form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]` | recommended | — | unbalanced-translocation delins, cross-references `complex.md` (parse-only — foreign genomic) |
| `NC_000022.10:g.42522624_42522669delins42536337_42536382` | conformant | — | conversion, same-reference numeric-range payload — ferro would expand it to literal bases on a real reference, so not the recommended range form (parse-only — foreign genomic) |
| `NC_000012.11:g.6128892_6128954delins[NC_000022.10:g.17179029_17179091]` | recommended | — | conversion, bracketed cross-reference payload (parse-only — foreign genomic) |

## `delins.md:73-77` — a `GC` to `TG` change is not a di-nucleotide substitution

> No, this is not allowed.
> By definition, a substitution changes **one** nucleotide into **one** other nucleotide (see [Substitution](substitution.md)).
> The change `TGT`<code class="del">GC</code>`CA` to `TGT`<code class="ins">TG</code>`CA` should be described as `g.4_5delinsTG`, i.e. a deletion/insertion (delins).

Ferro: the di-nucleotide substitution spelling is refused (the substitution grammar admits one
reference base and one alternate base only); two adjacent changed nucleotides are described as one
delins — the same separation-zero merge as `delins.md:16`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.145CG>AT` | refused | — | the disallowed di-nucleotide substitution shape (mirrors the spec's `class="invalid"` `g.4GC>TG`; `c.145_146` is `CG`) — rejected at grammar. See `delins.md:16` for the recommended delins form of the same change |

## `delins.md:79-84` — the maximum-unchanged-nucleotides Q&A

> Yes, two variants separated by one or more nucleotides should preferably be described individually and not as a "delins" (unless they together affect one amino acid).

Ferro: this restates the `:17` separation rule as a preference ("**preferably**") and gives two
grounds — first, provenance ("the two variants may have been reported (or might occur) individually");
second, an axis-neutral annotation-overlap concern. Provenance is not recoverable by a normalizer, so
ferro re-derives from the resulting sequence and discloses the deviation rather than preserving the
input's spelling. The separation ferro measures is read off the re-derived partition, so two spellings
of one variant converge.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-is-a-property-of-the-spelling-not-of-the-variant](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Ferro reads the separation between changes off the partition it re-derives from the resulting sequence, not off the input's spelling, so two spellings of one variant converge on one output.
<!-- why:END:separation-is-a-property-of-the-spelling-not-of-the-variant -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | separated by far more than one nucleotide — described individually, `:17`'s recommended form under rule 2 (same row as `delins.md:17`) |

## `delins.md:86-89` — the BRCA1 `c.2077` Q&A (separation-zero sub+ins merge)

> The correct description of this variant is `NM_007294.3:c.2077delinsATA`.<br>
> **NOTE**: the answer was modified, i.e. the addition "However, since the variant is likely a combination of two other variants, it is acceptable to describe it as <code class="invalid">NM_007294.3:c.[2077G>A;2077_2078insTA]</code>." was removed.

Ferro: a substitution immediately followed by an insertion, separation **zero**, coalesces at one
locus — both members consume reference bases there. This is the strongest register the corpus has
for the merge: the committee had **permitted** the split spelling and then **withdrew** that
permission, correcting itself toward the merged form. `:86-89` is separation-zero, so it says nothing
about the `≥2` axis-scope question the `:44-47` carve-out settles.

**The executable twin does not reproduce the BRCA1 delins shape**, and the reason is worth stating.
The spec's `c.2077G>A` is a *genuine* substitution (reference `G` becomes `A`), so its combination
with the flush insertion is a true `delins`. The one transcript in the slice has `c.76 = A` sitting
inside an A-run (`c.75_77` is `AAA`), so `c.76A>G` plus a flush `insTA` re-derives from the resulting
sequence to a **pure insertion** — the deleted/substituted base is recovered by the run, nothing net
changes, and the minimal edit is `c.75_76insGT`, which ferro emits. That is the correct recommended
form for a variant that is, on this reference, purely insertional; it simply is not the delins the
BRCA1 locus produces. The parse-only BRCA1 rows below carry the merge; the executable rows carry the
collapse.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_007294.3:c.2077delinsATA` | recommended | — | the spec's own answer; a fixed point by construction (parse-only — `NM_007294.3` not in the slice) |
| `NM_007294.3:c.[2077G>A;2077_2078insTA]` | recommended | — | the withdrawn split — the committee removed the passage permitting it; ferro merges it to the delins above (parse-only — merge target is the recommended form) |
| `NM_004006.3:c.[76A>G;76_77insTA]` | recommended | `NM_004006.3:c.75_76insGT` | executable twin of the sub+ins shape: `c.76` sits in the `c.75_77` A-run, so the sub+ins re-derives to a pure insertion — the minimal edit `c.75_76insGT`, not the BRCA1 delins |
| `NM_004006.3:c.76delinsGTA` | recommended | `NM_004006.3:c.75_76insGT` | the delins of `GTA` at `c.76` (`A`, in the A-run) is sequence-identical to inserting `GT`; ferro collapses it to the minimal insertion form |

## `delins.md:68-71` — the "indel" terminology Q&A

> The term "indel" is not used in HGVS nomenclature (see [Glossary](../../background/glossary.md)).

Descriptive — the term "indel" is not used in HGVS nomenclature (it is ambiguous across
disciplines). No verdict row; nothing to enforce.
