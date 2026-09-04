# Substitution — ferro's reading

ferro's reading of the HGVS **substitution** recommendations on the transcript (`r.`) axis,
clause by clause — each spelling with the form ferro normalizes it to and a verdict on that
output. New here? See [How to read a page](../../reading-guide.md) for the verdicts, the table
conventions, and the recurring terms.

*DNA twin: [Substitution (`c.`/`g.`)](../DNA/substitution.md).*

The RNA axis mirrors DNA substitution in substance — the same 1→1 definition, the same
separation/codon-exception rule, the same ban on the multi-base spelling and on the "polymorphism"
slash — and diverges only in surface: lowercase letters, `u` for `t` (`general.md:49`), and
predicted-consequence parentheses. Where a section below governs the same ground as its DNA
sibling, it cites the same ledger record.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. Most of this page's
spec examples are written on that very accession, and on the slice their stated bases hold
(`r.76` is `a`, `r.85` is `u`, `r.123` is `c`, `r.897` is `u`, `r.1388` is `g`) — so they run as
written. Two do not: the spec's `NM_004006.1` rows are on a version the slice does not carry
(parse-only, `—`, each with an executable `.3` twin), and `r.*41` reads `g` on the slice rather
than the spec's `u`, so that example is shown on a stated-flank twin. The `:15` Notes bullet and
the Discussion's cross-level authoring guidance (`:66-68`, `:77-80`) are reporting advice rather
than normalization rules over a single `r.` description; `:15` is kept below as a one-row
section, the two Q&A items are omitted as descriptive.

## `substitution.md:5` — definition: one nucleotide for one

> Substitution: a sequence change where, compared to a reference sequence, **one** nucleotide is
> replaced by **one** other nucleotide.

Ferro: a substitution is exactly 1→1 on the `r.` axis, spelled lower-case with `u` for `t`
(`general.md:49`). IUPAC ambiguity codes stand in for the single replacing base exactly as on the
DNA axis — but only the fifteen symbols the RNA table assigns; `x`, which the DNA table lists as
alignment-only and the RNA table does not list at all, is refused. The rule applies uniformly
across the transcript's numbering zones — 5'UTR (`-n`), CDS, and 3'UTR (`*n`).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[rna-axis-alignment-only-symbol-reach](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A non-leading lower-case 'x', which is not an assigned RNA nucleotide symbol, is refused in an 'r.' description, mirroring the refusal of the alignment-only 'X' on the DNA axes.
<!-- why:END:rna-axis-alignment-only-symbol-reach -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>c` | recommended | self | the spec's own example (`:25`); `r.76` is `a` on the slice — canonical single-base substitution |
| `NM_004006.3:r.54g>h` | recommended | self | replacement is IUPAC `h` = `a`, `c`, or `u` — one of the fifteen assigned RNA symbols (`r.54` is `g`) |
| `NM_004006.3:r.54g>x` | refused | — | `x` is not an assigned RNA symbol — rejected at parse (the RNA table has no counterpart to the DNA table's daggered `X`) |
| `NM_004006.1:r.-14a>c` | recommended | — | the spec's own 5'UTR example (`:37`), on a transcript version the slice does not carry — parse-only here |
| `NM_004006.3:r.-14a>c` | recommended | self | executable twin on the slice's version: `r.-14` is `a` on `NM_004006.3` too, 14 nt 5' of the `aug` |
| `NM_004006.3:r.*41u>a` | recommended | — | the spec's own 3'UTR example (`:40`). On the committed slice `r.*41` is `g`, not the `u` the example states, so the row is parse-only; the executable twin is one base 5' |
| `NM_004006.3:r.*40u>a` | recommended | self | stated-flank twin: `r.*40` is `u` on the slice (`r.*39_*43` is `augga`), 40 nt 3' of the `uag` stop |
| `NM_004006.3:r.76aa>ug` | refused | — | two nt — violates 1→1; rejected as "deprecated multi-base substitution syntax" (see `:16`) |

## `substitution.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein
> level may be given in addition.

Ferro: guidance to authors about which level(s) to report, not a constraint on how a given `r.`
description, once written, is normalized. An `r.` description stands on its own.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>c` | recommended | self | a lone RNA-level description is fully valid; ferro does not require an accompanying `c.`/`g.` form |

## `substitution.md:16` — two or more nucleotides is a delins

> substitutions involving two or more consecutive nucleotides are described as deletion/insertions
> (delins)

Ferro: a multi-base replacement is never a substitution on the `r.` axis either; write it as a
delins. The substitution-style spellings are marked `class="invalid"` in the spec's own worked
example (`:28-29`, below).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76_77delinsug` | recommended | self | the two-base change written as a delins (`r.76_77` is `aa`; both bases change, so it is a genuine delins on the slice) |
| `NM_004006.3:r.76_77aa>ug` | refused | — | rejected: "deprecated multi-base substitution syntax" — use delins |
| `NM_004006.3:r.76aa>ug` | refused | — | same rejection, the single-position variant of the invalid spelling |

## `substitution.md:17-19` — separation, and the one-codon exception

> - two substitutions separated by one or more nucleotides should be described individually and
>   not as a "delins".
>     - **exception**: two variants separated by one nucleotide, together affecting one amino
>       acid, should be described as a "delins" (e.g., `r.142_144delinsugg` (`p.Arg48Trp`)).<br>
>       **NOTE**: this prevents tools predicting the consequences of a variant to make
>       conflicting and incorrect predictions of two different substitutions at one position.

Ferro: the separation rule (ruleset rule 2), exactly as on the DNA axis; the codon exception folds
two substitutions into one delins when they sit one nucleotide apart and together change one amino
acid, applied triplet-precisely within the reading frame — on this document's own authority for
the `r.` axis. The spec's worked example states codon 48 as `cga`; on `NM_004006.3` it is `agg`,
so the literal `r.142_144delinsugg` changes only `r.142` there and ferro re-derives it to the
single substitution it actually is. The executable rows use the same-shaped `agg` → `ugc`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[delins-codon-carve-out-gap-one](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually.
>
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:separation-rule-force-modal-or-negation,delins-codon-carve-out-gap-one,canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.142_144delinsugc` | recommended | self | stated-flank twin of the spec's codon-48 example: `agg` → `ugc`, `p.(Arg48Cys)` — one codon, fixed point |
| `NM_004006.3:r.[142a>u;144g>c]` | recommended | `NM_004006.3:r.142_144delinsugc` | one nucleotide apart, one amino acid — ferro normalizes the split to the recommended delins |
| `NM_004006.3:r.142_144delinsugg` | recommended | `NM_004006.3:r.142a>u` | the spec's literal spelling on this transcript's bases: codon 48 is `agg`, so only `r.142` changes and the variant is one substitution — ferro re-derives from the resulting sequence rather than preserving the spelling |
| `NM_004006.3:r.[100a>u;104a>g]` | recommended | self | separated by three nucleotides and in different codons (34 and 35) — no codon-exception antecedent; stays split, the recommended form under rule 2 |

**Confluence.** `{r.[142a>u;144g>c], r.142_144delinsugc} → r.142_144delinsugc`.

See also → `RNA/delins.md:17-19`, where the same exception is adjudicated on its own document.

## `substitution.md:20` — no change is `=`, not a substitution

> nucleotides that have been tested and found **not changed** are described as

Ferro: a tested, unchanged position is `=`; a no-change substitution normalizes to it. Both the
single-position and the range form are fixed points.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123=` | recommended | self | the spec's own example (`:34`): `r.123` screened, unchanged — and it is `c` on the slice, as `:35` states |
| `NM_004006.3:r.123c>c` | recommended | `NM_004006.3:r.123=` | a no-change substitution normalized to the recommended `=` (constructed by analogy to the DNA axis) |
| `NM_004006.3:r.4567_4569=` | recommended | self | the range form the clause itself gives; three positions screened and unchanged |
| `NM_004006.3:r.109a=` | recommended | self | the clause's base-bearing form (`r.109u=` in the spec, generic; `r.109` is `a` on the slice) |

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

Ferro: two adjacent changed nucleotides (separation zero) are one delins — the same
separation-zero merge the DNA axis makes at `DNA/substitution.md:32`, here on RNA's own worked
example. The two substitution-style spellings are `class="invalid"` and refused.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76_77delinsug` | recommended | self | the spec's own delins; `r.76_77` is `aa` on the slice, both bases change — fixed point |
| `NM_004006.3:r.[76a>u;77a>g]` | recommended | `NM_004006.3:r.76_77delinsug` | the split spelling of the same change — adjacent, separation zero — merges to the recommended delins |
| `NM_004006.3:r.76_77aa>ug` | refused | — | `class="invalid"` — the range-form multi-base substitution |
| `NM_004006.3:r.76aa>ug` | refused | — | `class="invalid"` — the single-position multi-base substitution |

**Confluence.** `{r.[76a>u;77a>g], r.76_77delinsug} → r.76_77delinsug`.

## `substitution.md:31-32` — predicted consequence in parentheses

> the predicted consequences on RNA level is a substitution of the `g` nucleotide at `r.1388` with
> an `a`.

Ferro: a predicted (not directly observed) RNA-level substitution is wrapped in uncertainty
parentheses and preserved as written — the parentheses are a claim about provenance that a
normalizer cannot strengthen or weaken. `r.1388` is `g` on the slice, as the example states.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.(1388g>a)` | recommended | self | the spec's own predicted-consequence form, preserved |
| `NM_004006.3:r.1388g>a` | recommended | self | the same change asserted as observed — also a fixed point; ferro neither adds nor removes the parentheses |

## `substitution.md:43-45` — allele: two transcripts from one DNA variant

> two different transcripts, `r.897u>g` and `r.832_960del`, derive from one variant

Ferro: an allele bracket may group independently-derived `r.` descriptions from splice variants of
one underlying DNA change (the `,` separator: distinct transcripts, not cis members); each member
spells its own `r.` substitution or deletion. `r.897` is `u` on the slice, as the example states.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[897u>g,832_960del]` | recommended | self | the spec's own cross-transcript grouping; both members are already fixed points, so the allele is one |

## `substitution.md:47-54` — reserved notations: undetected, splicing affected, unpredictable

> no RNA from the variant allele could be detected

Ferro: `r.0`, `r.spl` and `r.?` are reserved literal tokens for, respectively, no RNA detected,
splicing likely affected, and an expected-but-unpredictable effect — not substitution spellings
themselves, but part of the same axis vocabulary. Each is preserved verbatim; there is nothing to
re-derive.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.1:r.0` | recommended | — | the spec's own accession (`:47`), a version the slice does not carry — parse-only here |
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
write the reference allele first (`:58`, `:62`), and ferro's output is reference-first when the
input is. `r.85` is `u` on the slice, as the example states.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.85=/u>c` | recommended | self | mosaic: reference `=` written first, then `u>c` |
| `NM_004006.3:r.85=//u>c` | recommended | self | chimeric: a mix of `r.85=` and `r.85u>c` cells |
| `NM_004006.3:r.85u>c/=` | conformant | self | valid, but `:58` writes the reference first; ferro does not yet reorder this to the recommended `r.85=/u>c`. Same limitation as the DNA axis — tracked by [#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034), filed on `DNA/substitution.md:49`, whose NOTE `:58` restates word for word |

## `substitution.md:71-73` — Q&A: polymorphisms are described as `r.76a>g`

> No, all substitutions are described as `r.76a>g`.

Ferro: the Discussion restates `:21` — the slash spelling is historical and invalid; every
substitution, "polymorphic" or not, is written with `>`. A description is neutral about
frequency, which is provenance no normalizer holds.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76a>g` | recommended | self | the Q&A's own answer; `r.76` is `a` on the slice |
| `NM_004006.3:r.76a/g` | refused | — | the historical "polymorphism" form — rejected at parse |
