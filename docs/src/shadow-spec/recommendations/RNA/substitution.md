# Substitution — ferro's reading

ferro's reading of the HGVS **substitution** recommendations on the transcript (`r.`) axis,
clause by clause — each spelling with the form ferro normalizes it to and a verdict on that
output. New here? See [How to read a page](../../reading-guide.md) for the verdicts, the table
conventions, and the recurring terms.

*DNA twin: [Substitution (`c.`/`g.`)](../DNA/substitution.md).*

The RNA axis mirrors DNA substitution — the same 1→1 definition, the same
separation/codon-exception rule, the same ban on the multi-base spelling and the "polymorphism"
slash — and diverges only in surface: lowercase letters, `u` for `t` (`general.md:49`), and
predicted-consequence parentheses.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. A few spec examples
cite accessions or bases the slice doesn't carry; those rows are parse-only (`—`), with an
executable twin shown alongside where one exists.

## `substitution.md:5` — definition: one nucleotide for one

> Substitution: a sequence change where, compared to a reference sequence, **one** nucleotide is
> replaced by **one** other nucleotide.

Ferro: a substitution is exactly 1→1 on the `r.` axis, lower-case with `u` for `t`
(`general.md:49`). IUPAC ambiguity codes may stand in for the replacing base, but only the fifteen
symbols the RNA table assigns — `x` is refused, since the RNA table has no alignment-only symbol.
Applies uniformly across 5'UTR, CDS, and 3'UTR numbering.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[rna-axis-alignment-only-symbol-reach](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A non-leading lower-case 'x', which is not an assigned RNA nucleotide symbol, is refused in an 'r.' description, mirroring the refusal of the alignment-only 'X' on the DNA axes.
<!-- why:END:rna-axis-alignment-only-symbol-reach -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>c` | recommended | self | the spec's own example — canonical single-base substitution |
| `NM_004006.3:r.54g>h` | recommended | self | replacement is IUPAC `h` (`a`/`c`/`u`) — one of the fifteen assigned RNA symbols |
| `NM_004006.3:r.54g>x` | refused | — | `x` is not an assigned RNA symbol — rejected at parse |
| `NM_004006.1:r.-14a>c` | recommended | — | the spec's own 5'UTR example, on a transcript version the slice lacks — parse-only |
| `NM_004006.3:r.-14a>c` | recommended | self | executable twin on the slice's version, 14 nt 5' of the start codon |
| `NM_004006.3:r.*41u>a` | recommended | — | the spec's own 3'UTR example; the slice's base differs, so this row is parse-only — see the executable twin below |
| `NM_004006.3:r.*40u>a` | recommended | self | stated-flank twin on the slice, 40 nt 3' of the stop codon |
| `NM_004006.3:r.76aa>ug` | refused | — | two nt — violates 1→1; rejected as deprecated multi-base substitution syntax (see `:16`) |

## `substitution.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein
> level may be given in addition.

Ferro: guidance for authors on which level(s) to report — not a constraint on normalizing an `r.`
description, which stands on its own.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>c` | recommended | self | a lone RNA-level description is valid; no accompanying `c.`/`g.` form required |

## `substitution.md:16` — two or more nucleotides is a delins

> substitutions involving two or more consecutive nucleotides are described as deletion/insertions
> (delins)

Ferro: a multi-base replacement is never a substitution on `r.`; it is written as a delins (worked
further at `:28-29`).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76_77delinsug` | recommended | self | the two-base change written as a delins — both bases change, a genuine delins |
| `NM_004006.3:r.76_77aa>ug` | refused | — | rejected: "deprecated multi-base substitution syntax" — use delins |
| `NM_004006.3:r.76aa>ug` | refused | — | same rejection, the single-position variant of the invalid spelling |

## `substitution.md:17-19` — separation, and the one-codon exception

> - two substitutions separated by one or more nucleotides should be described individually and
>   not as a "delins".
>     - **exception**: two variants separated by one nucleotide, together affecting one amino
>       acid, should be described as a "delins" (e.g., `r.142_144delinsugg` (`p.Arg48Trp`)).<br>
>       **NOTE**: this prevents tools predicting the consequences of a variant to make
>       conflicting and incorrect predictions of two different substitutions at one position.

Ferro: the separation rule (rule 2), as on the DNA axis; the codon exception folds two
substitutions into one delins when they sit one nucleotide apart and together change one amino
acid, applied triplet-precisely within the reading frame.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[delins-codon-carve-out-gap-one](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually.
>
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro applies the edit set to the reference and re-derives the description from the resulting sequence — subject to every explicit spec tie-break, and breaking any remaining tie toward the already-shipped form — rather than preserving the input's spelling.
<!-- why:END:separation-rule-force-modal-or-negation,delins-codon-carve-out-gap-one,canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.142_144delinsugc` | recommended | self | stated-flank twin of the spec's codon-48 example (`agg`→`ugc`, `p.(Arg48Cys)`) — one codon, fixed point |
| `NM_004006.3:r.[142a>u;144g>c]` | recommended | `NM_004006.3:r.142_144delinsugc` | one nucleotide apart, one amino acid — ferro normalizes the split to the recommended delins |
| `NM_004006.3:r.142_144delinsugg` | recommended | `NM_004006.3:r.142a>u` | the spec's literal spelling; on this transcript only `r.142` actually changes, so ferro re-derives the single substitution rather than preserving the spelling |
| `NM_004006.3:r.[100a>u;104a>g]` | recommended | self | separated by three nucleotides in different codons — no codon-exception antecedent; stays split under rule 2 |

**Confluence.** `{r.[142a>u;144g>c], r.142_144delinsugc} → r.142_144delinsugc`.

See also → `RNA/delins.md:17-19`, where the same exception is adjudicated on its own document.

## `substitution.md:20` — no change is `=`, not a substitution

> nucleotides that have been tested and found **not changed** are described as

Ferro: a tested, unchanged position is `=`; a no-change substitution normalizes to it. Both the
single-position and the range form are fixed points.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123=` | recommended | self | the spec's own example — position screened and found unchanged |
| `NM_004006.3:r.123c>c` | recommended | `NM_004006.3:r.123=` | a no-change substitution normalizes to the recommended `=` |
| `NM_004006.3:r.4567_4569=` | recommended | self | the range form the clause itself gives; three positions screened and unchanged |
| `NM_004006.3:r.109a=` | recommended | self | the clause's base-bearing form (generic in the spec) |

## `substitution.md:21` — polymorphisms are not `a/g`

> it is not correct to describe "_polymorphisms_" as

Ferro: the slash form for "polymorphism" is not valid HGVS on the `r.` axis either; ferro rejects
it at parse in strict mode.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>g` | recommended | self | the correct form (`:73`) |
| `NM_004006.3:r.76a/g` | refused | — | the disallowed "polymorphism" slash — `class="invalid"`; rejected at parse |

## `substitution.md:28-29` — the worked delins example, and its two invalid spellings

> - **`NM_004006.3:r.76_77delinsug`**<br>
>   **NOTE**: based on the definition of a substitution, i.e. **one** nucleotide replaced by
>   **one** other nucleotide, this change can not be described as a substitution like
>   <code class="invalid">r.76_77aa>ug</code> or <code class="invalid">r.76aa>ug</code>.

Ferro: two adjacent changed nucleotides (separation zero) merge into one delins — the DNA axis
makes the same merge. The substitution-style spellings are `class="invalid"` and refused.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76_77delinsug` | recommended | self | the spec's own delins — both bases change, a fixed point |
| `NM_004006.3:r.[76a>u;77a>g]` | recommended | `NM_004006.3:r.76_77delinsug` | the split spelling of the same change — adjacent, separation zero — merges to the recommended delins |
| `NM_004006.3:r.76_77aa>ug` | refused | — | `class="invalid"` — the range-form multi-base substitution |
| `NM_004006.3:r.76aa>ug` | refused | — | `class="invalid"` — the single-position multi-base substitution |

**Confluence.** `{r.[76a>u;77a>g], r.76_77delinsug} → r.76_77delinsug`.

## `substitution.md:31-32` — predicted consequence in parentheses

> the predicted consequences on RNA level is a substitution of the `g` nucleotide at `r.1388` with
> an `a`.

Ferro: a predicted (not directly observed) RNA-level substitution is wrapped in uncertainty
parentheses and preserved as written — the parentheses assert provenance, which a normalizer can't
strengthen or weaken.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.(1388g>a)` | recommended | self | the spec's own predicted-consequence form, preserved |
| `NM_004006.3:r.1388g>a` | recommended | self | the same change asserted as observed — also a fixed point; ferro doesn't add or remove parentheses |

## `substitution.md:43-45` — allele: two transcripts from one DNA variant

> two different transcripts, `r.897u>g` and `r.832_960del`, derive from one variant

Ferro: an allele bracket may group independently-derived `r.` descriptions from splice variants of
one underlying DNA change (the `,` separator marks distinct transcripts, not cis members); each
member spells its own `r.` edit.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[897u>g,832_960del]` | recommended | self | the spec's own cross-transcript grouping; both members are already fixed points, so the allele is one |

## `substitution.md:47-54` — reserved notations: undetected, splicing affected, unpredictable

> no RNA from the variant allele could be detected

Ferro: `r.0`, `r.spl` and `r.?` are reserved tokens — no RNA detected, splicing likely affected,
and an expected-but-unpredictable effect, respectively. Each is preserved verbatim.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.1:r.0` | recommended | — | the spec's own accession, a version the slice lacks — parse-only |
| `NM_004006.3:r.0` | recommended | self | executable twin: no RNA from the variant allele detected — reserved notation, preserved |
| `NM_004006.3:r.spl` | recommended | self | splicing likely affected, RNA not analysed — reserved notation, preserved |
| `NM_004006.3:r.?` | recommended | self | an effect is expected but not reliably predictable — reserved notation, preserved |

## `substitution.md:56-62` — mosaic and chimeric

> - **`NM_004006.3:r.85=/u>c`**<br>
>   a mosaic case where at position 85, besides the normal sequence (a `u`, described as `=`),
>   also transcripts are found containing a `c` (`r.85u>c`).<br>
>   **NOTE**: irrespective of the frequency in which each nucleotide was found, the reference is
>   always described first.

Ferro: mosaic (`/`) and chimeric (`//`) mixtures are valid on the `r.` axis; the recommendations
write the reference allele first, and ferro reorders a variant-first substitution to that form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.85=/u>c` | recommended | self | mosaic: reference `=` written first, then `u>c` |
| `NM_004006.3:r.85=//u>c` | recommended | self | chimeric: a mix of `r.85=` and `r.85u>c` cells |
| `NM_004006.3:r.85u>c/=` | recommended | `NM_004006.3:r.85=/u>c` | variant-first spelling; ferro reorders it to the recommended reference-first form, since `substitution.md:58` writes the reference first ([#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)) |

## `substitution.md:71-73` — Q&A: polymorphisms are described as `r.76a>g`

> No, all substitutions are described as `r.76a>g`.

Ferro: the Discussion restates `:21` — the slash spelling is historical and invalid; every
substitution, "polymorphic" or not, is written with `>`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>g` | recommended | self | the Q&A's own answer |
| `NM_004006.3:r.76a/g` | refused | — | the historical "polymorphism" form — rejected at parse |
