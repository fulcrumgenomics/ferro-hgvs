# Deletion — ferro's reading

ferro's reading of the HGVS **deletion** recommendations on the transcript (`r.`) axis, clause
by clause — each spelling with the form ferro normalizes it to and a verdict on that output.
New here? See [How to read a page](../../reading-guide.md) for the verdicts, the table
conventions, and the recurring terms.

*DNA twin: [Deletion (`c.`/`g.`)](../DNA/deletion.md).*

No ledger record cites an `RNA/deletion.md` clause directly (the exception is `:47-49`), so most
sections here are CONFIRM-by-inspection against the spec text and the shipped code, not an
adjudicated ruling.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. Spec examples on
other accessions the slice doesn't carry are parse-only (`—`).

## `deletion.md:5` — definition

> Deletion: a sequence change where, compared to a reference sequence, one or more nucleotides are
> not present (deleted).

Ferro: a deletion removes one or more nucleotides; the description itself carries only the deleted
range, never the deleted bases (see `:30-39` below for the payload-spelling NOTE).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.77del` | recommended | self | single-nucleotide RNA deletion; already 3'-most within its homopolymer run, a fixed point |
| `LRG_199t1:r.10del` | recommended | — | the spec's own one-nucleotide example (parse-only here) |

## `deletion.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein
> level may be given in addition.

Ferro: guidance for authors on which level(s) to report — not a constraint on normalizing an `r.`
description, which stands on its own.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.77del` | recommended | self | a lone RNA-level description is valid; no accompanying `c.`/`g.` form required |

## `deletion.md:16` — the range must name two different positions

> `position(s)_deleted` should contain **two different positions**, e.g., `123_126`, not `123_123`.

Ferro: rule 2 (lowercase "should"); a same-position range is repaired to the single-position form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123del` | recommended | self | the correct single-position form — no run to shift along |
| `NM_004006.3:r.123_123del` | refused | — | same-position range — rejected at parse (`W4003 SinglePositionRange`) in strict mode; lenient repairs to `r.123del` |

## `deletion.md:17` — the range is listed 5' to 3'

> the `position(s)_deleted` should be listed from **5' to 3'**, e.g., `123_126`, not `126_123`.

Ferro: a reversed, non-circular range has no exception on the `r.` axis (the `<high>_<low>` form
is admitted only for `o.`/`m.` circular references), so it is refused rather than reordered.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_126del` | recommended | self | correctly ordered range — already 3'-most |
| `NM_004006.3:r.126_123del` | refused | — | reversed range on a non-circular reference — rejected at parse |

## `deletion.md:18-20` — the 3'rule, and the exon/exon-junction exception that does not apply here

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily
>   assigned to have been changed (**3'rule**).
>     - the 3'rule also applies for changes in single residue stretches and tandem repeats (nucleotide
>       or amino acid).<br>
>       **NOTE**: the exception to the 3'rule for deletions around exon/exon junctions (see
>       [Deletions](../DNA/deletion.md)) does not apply when describing variants based on an RNA
>       reference sequence.

Ferro: the 3'rule applies on `r.` as elsewhere, including across single-residue runs and tandem
repeats. The NOTE is a separate instruction: it switches off the DNA exon/exon-junction 3'rule
exception on the RNA axis, so an `r.` deletion must still shift to its 3'-most position even
across an intron.

**Ferro currently gets this wrong.** It clamps at the exon/exon junction instead of crossing it,
matching the DNA behavior the NOTE says should not apply here. The output is still valid HGVS —
a wrong preference, not a malformed string — so the verdict is `conformant`, not `bug`. Filed as
[#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76del` | recommended | `NM_004006.3:r.77del` | 3'rule over a single-residue stretch — the deletion lands on the run's 3'-most base |
| `NM_004006.3:r.1470_1472del` | recommended | `NM_004006.3:r.1476_1478del` | 3'rule over a tandem repeat — the deleted unit lands on the 3'-most copy |
| `LRG_2t1:r.1034_1036del` | recommended | — | the spec's own 3'-most, junction-crossing form; no committed fixture for this accession — parse-only |
| `LRG_2t1:r.1033_1035del` | conformant | — | the spec marks this `class="invalid"` — the 3'rule should cross the junction to `r.1034_1036del`, but ferro wrongly clamps and keeps this spelling as a fixed point ([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211)); recommended is the 3'-most-crossing form. Code-derived — no committed fixture for this accession exercises the row |

No executable row on the committed slice can show the wrong behavior — it's built as a single flat
exon, so no `r.` deletion on it reaches an exon/exon junction.

**Non-confluence.** `{r.1033_1035del, r.1034_1036del}` describe one variant but do not converge
under the current defect — ferro treats `r.1033_1035del` as its own fixed point instead of
normalizing it to `r.1034_1036del`.

See also → `RNA/duplication.md:24` (identical NOTE for duplications, same defect shape).

## `deletion.md:21-22` — uncertain deletions

> see [Uncertain](../uncertain.md); when the position and/or the sequence of a deletion has not
> been defined, a description may have a format like `r.(100_150)delN[15]`.

Ferro: an uncertain deletion's parenthesised range is preserved verbatim — the 3'rule has nothing
determinate to shift, since neither the exact position nor the exact length is asserted.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.(100_150)delN[15]` | recommended | self | uncertain form, preserved — the spec's example is generic |

## `deletion.md:26-27` — one nucleotide

> a deletion of the `u` at position `r.10` in the reference sequence `LRG_199t1`.

Ferro: the spec's own single-nucleotide worked example.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.10del` | recommended | — | spec's own example; already 3'-most by construction (parse-only here) |

## `deletion.md:30-39` — several nucleotides: three worked examples

> - **`NM_004006.2:r.6_8del`**<br>
>   a deletion of nucleotides `r.6` to `r.8` in the reference sequence `NM_004006.2`.<br>
>   **NOTE**: the recommendation is not to describe the variant as
>   <code class="invalid">r.6_8deluug</code>, i.e., describe the deleted nucleotide sequence.
> - **`LRG_2t1:r.1034_1036del`**<br>
>   a deletion of nucleotides `r.1034` to `r.1036` (`uug`) in the reference sequence `LRG_2t1`.<br>
>   **NOTE**: since the 3'rule has to be applied, the variant, crossing the intron between
>   nucleotides `r.1035` and `r.1036`, is **not** described as
>   <code class="invalid">r.1033_1035del</code> (deletion `guu`).
> - **`LRG_199t1:r.(4072_5145del)`**<br>
>   the predicted deletion of exon 30 (starting at position `r.4072`) to exon 36 (ending at position
>   `r.5145`) of the _DMD_ gene; RNA has **not been analysed**.

Ferro: three independent points bundled in one Examples block — (i) a plain multi-nucleotide range
never carries the deleted payload sequence, (ii) the exon/exon-junction case already covered at
`:18-20`, and (iii) a predicted, unanalysed multi-exon deletion is preserved as an uncertain wrap
rather than resolved to concrete RNA coordinates.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:r.6_8del` | recommended | — | correct range-only form; spec's own accession (parse-only here) |
| `NM_004006.3:r.6_8del` | recommended | self | executable twin on the slice's version — already 3'-most |
| `NM_004006.2:r.6_8deluug` | refused | — | the disallowed payload-bearing spelling — `class="invalid"`; rejected at parse in strict mode |
| `NM_004006.3:r.6_8deluug` | refused | — | executable twin: strict rejects the payload; lenient drops it and normalizes to `r.6_8del` |
| `LRG_199t1:r.(4072_5145del)` | recommended | — | predicted multi-exon deletion, preserved as an uncertain wrap (parse-only here) |

**Exon/exon junction.** The `LRG_2t1` pair `{r.1034_1036del, r.1033_1035del}` is adjudicated at
`deletion.md:18-20` above; it is not repeated here.

## `deletion.md:41-43` — mosaic case

> - **`LRG_199t1:r.=/6_8del`**<br>
>   a mosaic case where from position `r.6` to `r.8`, besides the normal sequence, also transcripts
>   are found containing a deletion of this sequence.<br>
>   **NOTE**: for the predicted consequences of a variant, the description is
>   `LRG_199t1:r.(=/6_8del)`.

Ferro: mosaic (`=/`) syntax is valid on `r.`; the `6_8del` member inside it is subject to the same
3'rule as a standalone deletion. The predicted-consequence wrap (`(=/6_8del)`) is a separate,
uncertain form, preserved as written.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.=/6_8del` | recommended | — | mosaic form, spec's own accession (parse-only here) |
| `LRG_199t1:r.(=/6_8del)` | recommended | — | predicted-consequence wrap, spec's own accession (parse-only here) |
| `NM_004006.3:r.=/6_8del` | recommended | self | executable twin: the `6_8del` member is already 3'-most on this transcript (see `:30-39`) |
| `NM_004006.3:r.(=/6_8del)` | recommended | self | executable twin: the predicted-consequence wrap is preserved |

## `deletion.md:47-49` — a deletion is never written with a length suffix

> No, a deletion of more than one residue should mention the first and last residue deleted,
> separated using the range symbol ("_", underscore), e.g., `r.123_128del` and not
> <code class="invalid">r.123del6</code>.

Ferro: `class="invalid"` — rule 1. The repair is determinate: a length-`N` suffix starting at
position `p` becomes the range `p_(p+N-1)`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.124_129del` | recommended | self | the correct range form for a six-nucleotide deletion; already 3'-most. The spec's literal `r.123_128del` sits one base 5' on this transcript and would shift here under the 3'rule, so this row uses the fixed point to show the clause itself |
| `NM_004006.3:r.124del6` | refused | — | the disallowed length-suffix spelling — rejected at parse in strict mode (`W3011`); lenient repairs to `r.124_129del` |

## `deletion.md:51-54` — a deletion is never written by exon label

> A description like <code class="invalid">r.EX17del</code> has never been allowed. Descriptions
> should be specific and indicate the nucleotides affected by the change.

Ferro: `EX17` is not a position token in the grammar on any axis — an exon label is a fact about
annotation, not a coordinate the description can carry. No repair is possible; this is a parse
rejection, not a normalization.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.EX17del` | refused | — | not a valid position token — rejected at parse |
