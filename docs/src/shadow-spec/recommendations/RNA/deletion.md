# Deletion — ferro's reading

ferro's reading of `RNA/deletion.md`. The rules are HGVS's; ferro's job is to produce the form the
recommendations prefer. Verdicts describe **ferro's output**:

- **recommended** — ferro's output is the form the recommendations prefer (whether the input was
  already that form, or ferro normalized it there).
- **conformant** — ferro's output is valid HGVS but not *yet* the recommended form — a ferro
  limitation or a deliberate maintainer house choice among conformant forms, with a tracking
  issue where one exists.
- **refused** — the input is not valid HGVS; ferro rejects it in strict mode (correct behavior).
- **bug** — ferro's output is not valid HGVS (a defect). None on this page.

Each **Why** block is transcluded from the ruling ledger — the record's own one-line summary,
rendered here and linked to its full entry in
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md).
The reasoning lives once, in the ledger; it is never re-typed here.

No ledger record currently cites an `RNA/deletion.md` clause directly, so most sections below carry
no Why block — the reading is CONFIRM-by-inspection against the spec text and the shipped code, not
an adjudicated ruling. The one exception is `:47-49`, governed by a general mode-policy record.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
examples sit on `LRG_199t1`, `LRG_2t1` and `NM_004006.2`, none of which the slice carries, so those
rows are parse-only (`—`) — ferro cannot read their bases here.

## `deletion.md:5` — definition

> Deletion: a sequence change where, compared to a reference sequence, one or more nucleotides are
> not present (deleted).

Ferro: a deletion removes one or more nucleotides; the description itself carries only the deleted
range, never the deleted bases (see `:30-39` below for the payload-spelling NOTE).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.77del` | recommended | self | single-nucleotide RNA deletion; the 3'-most `a` of the `aaa` run at `r.75_77`, so already a fixed point |
| `LRG_199t1:r.10del` | recommended | — | the spec's own one-nucleotide example (parse-only here) |

## `deletion.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein
> level may be given in addition.

Ferro: this is guidance to authors about which level(s) to report, not a constraint on how a given
`r.` description, once written, is normalized. An `r.` description stands on its own.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.77del` | recommended | self | a lone RNA-level description is fully valid; ferro does not require an accompanying `c.`/`g.` form |

## `deletion.md:16` — the range must name two different positions

> `position(s)_deleted` should contain **two different positions**, e.g., `123_126`, not `123_123`.

Ferro: rule 2 (lowercase "should"); a same-position range is repaired to the single-position form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123del` | recommended | self | the correct single-position form (`r.123` is `c`, `r.124` is `a` — no run to shift along) |
| `NM_004006.3:r.123_123del` | refused | — | same-position range — rejected at parse (`W4003 SinglePositionRange`) in strict mode; lenient repairs to `r.123del` |

## `deletion.md:17` — the range is listed 5' to 3'

> the `position(s)_deleted` should be listed from **5' to 3'**, e.g., `123_126`, not `126_123`.

Ferro: a reversed, non-circular range has no exception on the `r.` axis (the `<high>_<low>` form
is admitted only for `o.`/`m.` circular references), so it is refused rather than reordered.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_126del` | recommended | self | correctly ordered range (`cagu`; `r.127` is `g`, so it is already 3'-most) |
| `NM_004006.3:r.126_123del` | refused | — | reversed range on a non-circular reference — rejected at parse |

## `deletion.md:18-20` — the 3'rule, and the exon/exon-junction exception that does not apply here

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily
>   assigned to have been changed (**3'rule**).
>     - the 3'rule also applies for changes in single residue stretches and tandem repeats (nucleotide
>       or amino acid).<br>
>       **NOTE**: the exception to the 3'rule for deletions around exon/exon junctions (see
>       [Deletions](../DNA/deletion.md)) does not apply when describing variants based on an RNA
>       reference sequence.

Ferro: the 3'rule itself applies on `r.` exactly as elsewhere, including across single-residue
runs and tandem repeats — the two executable rows below show each. The NOTE is a second, separate
instruction — it **switches off** the DNA exon/exon-junction 3'rule exception on the RNA axis, so
an `r.` deletion (or duplication) must still shift to its 3'-most position even when doing so
crosses an intron.

**Ferro currently gets the second half wrong.** `normalize_rna` (#334) applies the DNA exon/exon
clamp to `r.` del/dup edits — it never lets a deletion 3'-shift across an exon/exon junction on the
RNA axis — which is exactly what this NOTE says must not happen. No ledger record governs this
clause; the nearest record, `exon-junction-dup-converge-from-the-far-side`, decides only the
far-side pull-back on `c.`/`n.` and rests its own `r.` reach on `general.md:43` alone — the token
`RNA/deletion.md:20` (and `RNA/duplication.md:24`) explicitly withdraws. Filed as
[#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211). The output is valid HGVS —
it is the wrong *preference*, not a malformed string — so the verdict is `conformant`, not `bug`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76del` | recommended | `NM_004006.3:r.77del` | 3'rule over a single-residue stretch: `r.75_77` is `aaa`, so the deletion lands on `r.77` |
| `NM_004006.3:r.1470_1472del` | recommended | `NM_004006.3:r.1476_1478del` | 3'rule over a tandem repeat: `r.1470_1478` is `aca` x 3 (`r.1479` is `u`), so the deleted unit lands on the 3'-most copy |
| `LRG_2t1:r.1034_1036del` | recommended | — | the spec's own 3'-most, junction-crossing form (`del` `uug`); no committed `LRG_2t1` fixture, so parse-only here |
| `LRG_2t1:r.1033_1035del` | conformant | — | the spec marks this `class="invalid"` — on `r.` the exon/exon-junction exception should not apply, so the 3'rule requires crossing to `r.1034_1036del`. Ferro's `normalize_rna` wrongly clamps at the junction and keeps this spelling as a fixed point ([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211)); recommended is the 3'-most-crossing form. Code-derived: no committed `LRG_2t1` fixture exercises this row, and the `c.`/`n.` sibling clamp is pinned in `tests/it/issue_334_exon_junction_exception.rs`, which carries no `r.` del/dup case |

No executable `NM_004006.3` row can show the wrong behavior: the committed slice is built as a
single flat exon (`tests/fixtures/shadow_spec/build_fixture.py`), so no `r.` deletion on it can
reach an exon/exon junction. The clamp is enforced by code and, so far, by no test.

**Non-confluence.** The pair `{r.1033_1035del, r.1034_1036del}` describes one variant (deletion of
`guu`/one 3'-shifted frame of it) but, under the current defect, does not converge: ferro treats
`r.1033_1035del` as its own fixed point rather than normalizing it to the recommended
`r.1034_1036del`.

See also → `RNA/duplication.md:24` (identical NOTE for duplications, same defect shape).

## `deletion.md:21-22` — uncertain deletions

> see [Uncertain](../uncertain.md); when the position and/or the sequence of a deletion has not
> been defined, a description may have a format like `r.(100_150)delN[15]`.

Ferro: an uncertain deletion's parenthesised range is preserved verbatim — the 3'rule has nothing
determinate to shift, since neither the exact position nor the exact length is asserted.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.(100_150)delN[15]` | recommended | self | uncertain form, preserved (illustrative position; the spec's own example is generic) |

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
| `NM_004006.3:r.6_8del` | recommended | self | executable twin on the slice's version: `r.6_8` is the same `uug` on `NM_004006.3`, and `r.9` is `g`, so it is already 3'-most |
| `NM_004006.2:r.6_8deluug` | refused | — | the disallowed payload-bearing spelling — `class="invalid"`; rejected at parse in strict mode |
| `NM_004006.3:r.6_8deluug` | refused | — | executable twin: strict rejects the payload; lenient drops it and normalizes to `r.6_8del` (the spec's own spelling is pinned that way in `tests/it/error_mode_tests.rs::rna_lowercase_corrected`) |
| `LRG_199t1:r.(4072_5145del)` | recommended | — | predicted multi-exon deletion, preserved as an uncertain wrap (`RNA has not been analysed`); round-tripped by `rna_coding_consistency.rs::rna_predicted_multi_exon_del_round_trips` (parse-only here) |

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

**Why.**
<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.124_129del` | recommended | self | the correct range form for a six-nucleotide deletion (`agugac`); already 3'-most, since `r.124` is `a` and `r.130` is `c`. (The spec's literal `r.123_128del` is one base 5' of this on `NM_004006.3` — `r.123` and `r.129` are both `c` — so it would shift here under the 3'rule; the row uses the fixed point so the clause, not the shift, is what it shows) |
| `NM_004006.3:r.124del6` | refused | — | the disallowed length-suffix spelling — rejected at parse in strict mode (`W3011`); lenient repairs to `r.124_129del`. Pinned by `tests/it/issue_1079_point_size_suffix.rs` |

## `deletion.md:51-54` — a deletion is never written by exon label

> A description like <code class="invalid">r.EX17del</code> has never been allowed. Descriptions
> should be specific and indicate the nucleotides affected by the change.

Ferro: `EX17` is not a position token in the grammar on any axis — an exon label is a fact about
annotation, not a coordinate the description can carry. No repair is possible; this is a parse
rejection, not a normalization.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.EX17del` | refused | — | not a valid position token — rejected at parse. No dedicated negative-test pin found (minor gap; not a ledger matter) |
