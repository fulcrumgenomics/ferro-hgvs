# Repeated Sequences — ferro's reading

ferro's reading of the HGVS **repeated sequences** recommendations, clause by clause — each
spelling with the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Repeated Sequences (`r.`)](../RNA/repeated.md).*

Most of this page is CONFIRM-by-inspection against the spec text and the shipped code: the 3'rule
over a repeat tract, the mixed/composite listing form, the sequenced-versus-sized dichotomy, and
the strand-plus-3'rule unit selection are mechanical parser or normalizer behaviour with no clause
in tension. Two things carry a Why block, and both deserve prominence. The **multiple-of-3 coding
restriction** at `repeated.md:21-24` is the one clause with real teeth — inside the coding sequence
a short repeat may **not** be spelled as a repeat and is repaired to a `dup`/`ins` — and it is a
straight reading of the text, so it carries **no** Why block, only executable rows. The **range +
unit shape** is the open question: on the RNA page it is a live self-contradiction
(`rna-repeat-range-plus-unit-redundancy`, **undecided**), and while the DNA page's own text never
raises that conflict, the same shape recurs here (`c.-6_-3G[6]`), so it is discussed as **OPEN** at
`repeated.md:9`. The one genuinely adjudicated DNA choice is repeat-notation-versus-`dup` where a
tract-length change spells both ways (`canonical-form-choice-when-both-legal`, at
`repeated.md:36-44`).

**The grammar fact under the whole page.** `docs/syntax.yaml`'s `dna.rpt` production declares two
canonical spellings — a **single position plus a unit** ("Unique Repeat", `g.123CAG[23]`) and a
chained **Mixed Repeat** (`g.123_191CAG[19]CAA[4]`). It declares **no** range-plus-unit form
(`g.123_191CAG[23]`) even though that is the doc's dominant example shape, and **no**
positions-only form (`g.123_191[23]`, no unit) even though the parser accepts one
(`parse_repeat`'s optional unit). So the shipped parser is wider than the declared grammar; ferro
accepts the range+unit, single-position+unit, and range-only shapes alike, and the mismatch between
`syntax.yaml` and the doc's own worked examples is a documentation gap, not a ferro obligation.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. It carries two real
repeat tracts the spec never names: the 8-nucleotide A-stretch `c.5690_5697` (`AAAAAAAA`, the
coding-sequence 3'rule example `deletion.md` also leans on) and an `AC` dinucleotide run of 8
copies at `c.*460_*475` in the 3'UTR. The A-stretch anchors the **coding-restriction** row (a
mononucleotide repeat inside the CDS, where the restriction bites) and the `AC` tract anchors the
**UTR-carve-out** row (a dinucleotide repeat exempt because it sits in the 3'UTR). Every executable
output below is ferro's actual normalized output, verified against the bless harness. The
spec's own examples all sit on foreign accessions the slice does not carry — `NM_024312.4`,
`NM_000044.3`, `NM_023035.2`, `NM_000333.3`, `NM_002024.5`, `LRG_763t1` / `NM_002111.8`,
`NM_000492.3`, `NM_021080.3`, and the `NC_` / `NC_000014.8` genomic accessions — so those rows are
parse-only (`—`); ferro cannot read their bases here.

## `repeated.md:5` — definition: a repeat unit present several times

> Repeated sequence: a sequence where, compared to a reference sequence, a segment of **one or more** nucleotides (the repeat unit) is present several times, one after the other.

Ferro: the definition of the object — a unit of one or more nucleotides present in several tandem
copies ("one after the other", so strictly adjacent). The floor is unit length ≥ 1, so a
mononucleotide homopolymer run is a valid repeated sequence. It states no normalization obligation,
so there is no verdict row.

## `repeated.md:9` — Community Consultation NOTE, and the OPEN range + unit shape

> **NOTE**: a Community Consultation proposal is being prepared which will suggest to allow only the format where the **entire range** of the repeated sequence is indicated; so `g.123_191CAG[23]`, **not** `g.123CAG[23]`.

Ferro: a forward-looking NOTE describing a proposal that is being *prepared*, not a current rule.
Checked against `docs/consultation/` — there is **no numbered, dispositioned SVD-WG proposal** for
this; the range-versus-position question appears only under "Ongoing discussions", with no
committee decision for or against. So both the single-position form (`g.123CAG[23]`) and the
range+unit form (`g.123_191CAG[23]`) remain fully legal today, and a parser must **not** enforce
entire-range-only. Ferro accepts both. Both spec spellings are accession-less, so they are
discussed here in prose rather than forced into a verdict table (a bare spelling does not parse).

The related open question is the **range + unit** shape itself. On the RNA page this is a live
self-contradiction — `RNA/repeated.md:22` calls range-plus-unit invalid as redundant while `:27`
publishes exactly that shape as valid. The DNA page's own text never raises that conflict (it never
calls range+unit invalid), yet it carries the byte-parallel `c.-6_-3G[6]` at `:24` — the same
upstream text copied. So on the DNA axis the shape is presented as valid and ferro accepts it, but
whether it *should be* the canonical form is held **open**, tracked on the RNA-scoped record below
and owed upstream (`#466`).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[rna-repeat-range-plus-unit-redundancy](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Whether an RNA repeat may combine a position range with its own repeat unit is undecided — the recommendations both call that shape invalid as redundant and publish it as a worked example — and ferro's current lenient-mode agreement with the published form is an incidental effect of tract maximization rather than a resolution.
<!-- why:END:rna-repeat-range-plus-unit-redundancy -->

</details>

The record is RNA-scoped; it is cited here because the DNA page inherits the same shape by copy
(`c.-6_-3G[6]` at `:24`) while raising no conflict of its own — so the DNA reading is OPEN by
analogy, not decided. See `repeated.md:21-24` for where the carve-out that string demonstrates is
in fact correct (only the *spelling* is the open question).

## `repeated.md:17` — the size range of repeats

> - repeated sequences include both small (mono-, di-, tri-, etc., nucleotide) and larger (kilobase-sized) repeats.

Ferro: a scope statement — no ceiling or floor beyond the unit-length-≥1 of the definition. No
normalization obligation, so no verdict row.

## `repeated.md:18` — mixed repeats: range anchor, then a per-unit listing

> - for **mixed repeats**, the range of the repeat sequence is given followed by a listing of each repeat unit and the number of repeats in each unit; `NC_000012.11:g.112036755_112036823CTG[9]TTG[1]CTG[13]`.

Ferro: states, as a definite rule, that a mixed/composite repeat's anchor is a **range** (never a
single position), followed by a chained `unit[n]` listing. Ferro parses that shape. The clause does
**not** state whether the segmentation of a heterogeneous tract into named unit-runs is *unique* —
the same open shape the project has recorded for `delins` splits — but no ferro obligation is
adjudicated from it here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000012.11:g.112036755_112036823CTG[9]TTG[1]CTG[13]` | recommended | — | the spec's own mixed-repeat illustration — a range anchor followed by three unit-runs (parse-only here). Restated as a standalone Example at `repeated.md:85-86` |

## `repeated.md:19-20` — sequenced versus not-sequenced

> - `NM_000044.3:c.171_239GCA[34]` describes a repeated sequence containing 34 `GCA` units (sequenced, the reference sequence contains 23 `GCA` units).
>   `NM_000044.3:c.(92_331)insN[33]` describes an insertion of 33 nucleotides in the amplified region from position `c.92` to `c.331` (**not sequenced**), containing a repeated sequence of 24 `GCA` units in the reference sequence.

Ferro: the sequenced/not-sequenced dichotomy the page uses throughout — a range + concrete unit
when the tract was actually sequenced (the exact bases are known), and a parenthesized range +
`insN[#]` when only a size difference was measured (fragment sizing). The two are **not**
apply-equal: `c.171_239GCA[34]` denotes an exact tract, while `c.(92_331)insN[33]` denotes "33
nucleotides longer, bases undetermined" — different claims from different evidence — so there is no
confluence question to resolve and no Why block. This page is the *authority* for the `insN[#]` /
`delN[#]` shape; `uncertain.md` borrows it back.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_000044.3:c.171_239GCA[34]` | recommended | — | the sequenced form — a range + concrete `GCA` unit (parse-only here) |
| `NM_000044.3:c.(92_331)insN[33]` | recommended | — | the not-sequenced form — a parenthesized-range `insN[#]` sized insertion (parse-only here) |

## `repeated.md:21-24` — the multiple-of-3 coding restriction, and the UTR carve-out

> - **exception**: using a coding DNA reference sequence ("c." description), a repeated sequence variant description can be used only for repeat units with a length which is a multiple of 3, i.e. which can not affect the reading frame.
>   Consequently, use `NM_024312.4:c.2692_2693dup` and **not** <code class="invalid">NM_024312.4:c.2686A[10]</code>; use `NM_024312.4:c.1741_1742insTATATATA` and **not** <code class="invalid">NM_024312.4:c.1738TA[6]</code>.
>   This restriction only applies to the coding sequence, which does not include the introns or the UTR sequence.
>   As such, `NM_024312.4:c.-6_-3G[6]` is valid as the reading frame is not affected.

Ferro: **inside the coding sequence** (not introns, not UTR), a repeat may be described as a repeat
only when its unit length is a multiple of 3, so it cannot shift the reading frame. A shorter unit
(length 1 or 2) is `class="invalid"` there and must be re-expressed as a `dup`/`ins` describing the
same change. Outside the coding sequence — the 5'UTR (`c.-N`), the 3'UTR (`c.*N`), and introns —
the restriction is switched off and any unit length is legal. This is a self-contained clause: no
other recommendation backs or generalizes it, and the RNA sibling states the identical rule on its
own axis. Because nothing contests it, there is **no ledger record** and no Why block.

The coding-restriction gate (`unit_len % 3`) and the UTR exemption (`span_is_coding`) both live in
`src/normalize/`; the parse-only rows below are the spec's own worked repairs, and the two
`NM_004006.3` rows exercise both halves on the slice's real bases.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_024312.4:c.2686A[10]` | recommended | — | the spec's `class="invalid"` coding mononucleotide repeat (unit length 1) — ferro repairs it away from the barred repeat form to the spec's own recommended `c.2692_2693dup`, the only legal form in the CDS (parse-only, foreign accession) |
| `NM_024312.4:c.2692_2693dup` | recommended | — | the spec's own recommended rewrite of the row above — a `dup` describing the same change in-frame (parse-only here) |
| `NM_024312.4:c.1738TA[6]` | recommended | — | the spec's second `class="invalid"` coding repeat (di-nucleotide unit, length 2) — ferro repairs it to the spec's own recommended `c.1741_1742insTATATATA` in default/lenient (strict declines the barred spelling), mirroring the RNA sibling's di-nucleotide row (parse-only here) |
| `NM_024312.4:c.1741_1742insTATATATA` | recommended | — | the spec's recommended rewrite of the row above — an `ins` describing the same change (parse-only here) |
| `NM_024312.4:c.-6_-3G[6]` | conformant | — | the 5'UTR carve-out: a mononucleotide repeat (unit length 1) is exempt because `c.-6` sits in the 5'UTR and the reading frame is not affected, so `:24` calls it valid. This is **also** the range+unit shape the OPEN question at `repeated.md:9` concerns — the carve-out is correct, the range+unit *spelling* is what the undecided RNA record disputes. Parse-only here |
| `NM_004006.3:c.5690A[10]` | recommended | `NM_004006.3:c.5696_5697dup` | executable coding-restriction row on the slice's own A-stretch: `c.5690_5697` is `AAAAAAAA` (8 copies) in the coding sequence, so `A[10]` asserts 2 extra `A`s. Unit length 1 is not a multiple of 3 and `c.5690` is in the CDS, so the repeat form is barred; `dup` is the only legal form there, so ferro's `c.5696_5697dup` (the +2 mononucleotide expansion read as a duplication and 3'-shifted to the most-3' pair, `c.5698` is `T`) is the recommended output. Direct parallel to the spec's `c.2686A[10]`→`c.2692_2693dup` |
| `NM_004006.3:c.*460AC[9]` | recommended | `NM_004006.3:c.*475_*476dup` | executable UTR-carve-out row: `NM_004006.3` carries an `AC`×8 tract at `c.*460_*475` in the 3'UTR, so `c.*460AC[9]` asserts one extra `AC` copy. Unit length 2 is **not** a multiple of 3, yet the repeat is legal because `c.*460` is in the 3'UTR — the carve-out this section states. Ferro reads the +1-copy expansion as a duplication and applies the 3'rule, emitting the 3'-most `dup`; value mirrors the RNA sibling's measured `r.*475_*476dup` on the identical tract. |

No executable row can exercise the coding *di*-nucleotide case: the slice's short coding tandems
(`c.5_7` `TTT`, `c.11_13` `GGG`, both mononucleotide) are all unit length 1, so the length-2 branch
(`c.1738TA[6]`) stays parse-only on the spec's accession.

## `repeated.md:28-34` — unique repeat, sequenced

> - **`NC_000014.8:g.101179660_101179695TG[14]`**<br>
>   a repeated `TG` di-nucleotide sequence, starting at position `g.101179660` on human chromosome 14, with 14 `TG` copies.

Ferro: the plain sequenced unique-repeat form (a range + concrete unit), and its two-allele
cis-bracket, which composes with the standard `dna.alleles` syntax. No new rule beyond the mixed
form already covered; both spec examples sit on `NC_000014.8`, so both are parse-only.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000014.8:g.101179660_101179695TG[14]` | recommended | — | the sequenced `TG` di-nucleotide repeat, 14 copies (parse-only here) |
| `NC_000014.8:g.[101179660_101179695TG[14]];[101179660_101179695TG[18]]` | recommended | — | the di-allelic cis-bracket form — 14 copies on one allele, 18 on the other (parse-only here) |

## `repeated.md:36-44` — repeat expansion disorders: the repeat-versus-`dup` choice, and the 3'rule

> - **`NM_023035.2:c.6955_6993CAG[26]` (or `c.6955_6993dup`)**<br>
>   a repeated CAG tri-nucleotide sequence, starting at position `c.6955` in the _CACNA1A_ gene with 26 `CAG` copies (`p.(Gln2319[26])` or `p.(Gln2319_Gln2331dup)`).

Ferro: the CACNA1A example juxtaposes two legal spellings of one variant — the repeat form
`c.6955_6993CAG[26]` and `c.6955_6993dup` — with a bare "(or …)" and states **no** preference. A
+1-unit tandem expansion of the tract's terminal 9 nt *is* structurally a duplication of those 9
nt, so the two denote the same variant. No DNA-side clause states which is preferred (`general.md`
and `DNA/duplication.md` are both silent), so which one ferro emits is a `canonical-form` choice: it
derives the form from the resulting sequence rather than preserving the input's spelling. (The RNA
axis states the opposite of silence for the analogous case, preferring repeat notation over
`dup`/`del`; the DNA axis does not, which is why this is the DNA page's one genuinely adjudicated
representation choice.)

The ATXN7 row is the 3'rule over a repeat unit: literature calls it a `CAG` repeat, but on the
coding DNA reference the 3'rule selects the `AGC` rotation, so it is described `AGC[13]` — the same
strand-plus-3'rule unit selection `general.md` states generally. No dispute.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_023035.2:c.6955_6993CAG[26]` | recommended | — | the CACNA1A repeat form (parse-only here); an alternate legal spelling of the same variant as the `dup` below, and which one ferro emits is the `canonical-form` choice above |
| `NM_023035.2:c.6955_6993dup` | recommended | — | the `dup` spelling of the same +1-unit expansion (parse-only here) — both legal, one variant |
| `NC_000003.12:g.63912687_63912716AGC[13]` | recommended | — | the ATXN7 repeat, genomic axis — the `AGC` unit the 3'rule selects (parse-only here) |
| `NM_000333.3:c.89_118AGC[13]` | recommended | — | the same allele on the coding axis; `AGC`, not the literature `CAG` (parse-only here) |

## `repeated.md:46-52` — not sequenced: `insN[#]` and `delN[#]`

> - **`NC_000003.12:g.(63912602_63912844)insN[9]` / `NM_000333.3:c.(4_246)insN[9]`**<br>
>   a fragment containing the `AGC` repeat in the _ATXN7_ gene was amplified (from nucleotide `g.63912602` / `c.4` to `g.63912844` / `c.246`) and its size determined to be 9 nucleotides larger (`insN[9]`) compared to that of the reference sequence.<br>

Ferro: the sized-but-unsequenced forms — a parenthesized (uncertain) range with `insN[#]` for a
tract measured longer, or `delN[#]` for one measured shorter. The `insN[#]`/`delN[#]` shape parses;
the parenthesized range is preserved verbatim (nothing determinate to shift). No new rule.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000003.12:g.(63912602_63912844)insN[9]` | recommended | — | the "9 nt larger, not sequenced" genomic form (parse-only here) |
| `NM_000333.3:c.(4_246)insN[9]` | recommended | — | the same on the coding axis (parse-only here) |
| `NC_000003.12:g.(63912602_63912844)delN[15]` | recommended | — | the "15 nt smaller, not sequenced" form (parse-only here) |
| `NM_000333.3:c.(4_246)delN[15]` | recommended | — | the same on the coding axis (parse-only here) |

## `repeated.md:54-71` — mixed repeat reference sequence: FMR1

> - **`NM_002024.5:c.-128_-69GGC[10]GGA[1]GGC[9]GGA[1]GGC[10]`**<br>
>   a sequenced `GGC` tri-nucleotide repeat from position `c.-128` to `c.-69` contains 10 `GGC`, 1 `GGA`, 9 `GGC`, 1 `GGA`, and 10 `GGC` units (31 repeat units).

Ferro: the FMR1 block is the worked mixed-repeat form (a range anchor + per-unit listing), plus
three additional notational options and one negative example. All sit on `NM_002024.5`, so all are
parse-only.

- The explicit per-run listing (`GGC[10]GGA[1]…`) and the `GGC[68]GGA[1]GGC[10]` schematic are two
  descriptions of different-length alleles of the same tract.
- `c.-128_-69GGM[108]` uses the IUPAC ambiguity code `M` (A or C) as the repeat *unit*, standing in
  for `GGA`/`GGC` per iteration — a third way to spell a partially-determined tract. `general.md`
  admits IUPAC-IUBMB symbols, so the shape is not itself contested (whether the parser accepts an
  ambiguity code as a DNA repeat unit is a build-time check, not adjudicated here).
- `c.(-144_-16)insN[(1800_2400)]` composes the sized-insertion form with an uncertain copy-number
  range — the same conventions already established, composed.
- `c.-129CGG[79]` is marked `class="invalid"` by the spec: it "would cover only the sequence up to
  the first `AGG` interruption", so the tract cannot be described as a single pure `CGG` run. That
  invalidity is a fact about the reference bases (an interruption exists), which ferro cannot see on
  a foreign accession — so ferro parses the grammatically-valid string and cannot reject it here;
  the row is flagged for the fixup as behavior ferro cannot adjudicate on the slice.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_002024.5:c.-128_-69GGC[10]GGA[1]GGC[9]GGA[1]GGC[10]` | recommended | — | the fully-listed mixed repeat, 31 units (parse-only here) |
| `NM_002024.5:c.-128_-69GGC[68]GGA[1]GGC[10]` | recommended | — | the 79-unit `GGC`-framed schematic of a longer allele (parse-only here) |
| `NM_002024.5:c.-128_-69GGM[108]` | recommended | — | the IUPAC `M`-unit spelling of the partially-determined tract (parse-only here) |
| `NM_002024.5:c.(-144_-16)insN[(1800_2400)]` | recommended | — | sized insertion with an uncertain copy-number range (parse-only here) |
| `NM_002024.5:c.-129CGG[79]` | conformant | — | the spec's `class="invalid"` pure-`CGG` spelling — invalid because the reference tract is interrupted by `AGG`, a fact about the bases ferro cannot check on a foreign accession, so ferro parses it as valid HGVS (parse-only here; flagged for fixup — ferro cannot adjudicate the interruption without the sequence) |

## `repeated.md:73-78` — mixed repeat reference sequence: HTT

> - **`LRG_763t1:c.54_110GCA[23]`**<br>
>   a sequenced `GCA` tri-nucleotide repeat starting at position `c.54` contains 23 units, on protein level described as `NP_002102.4:p.(Gln18)[25]`.<br>

Ferro: the HTT `GCA` repeat — the `GCA` unit is again the 3'rule's rotation of the literature `CAG`.
Sequenced range + concrete unit; on `LRG_763t1`, so parse-only.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_763t1:c.54_110GCA[23]` | recommended | — | the sequenced `GCA` repeat, 23 units — the 3'rule rotation of the literature `CAG` (parse-only here) |

## `repeated.md:80-83` — mixed repeat reference sequence: CFTR intron 9

> **`NM_000492.3:c.1210-33_1210-6GT[11]T[6]`**<br>
> the mixed repeat sequence from position `c.1210-33` to `c.1210-6` contains 11 `GT` and 6 `T` copies.<br>

Ferro: a two-unit mixed repeat in intron 9 (a mixed repeat may be spelled `GT[11]T[6]`), and its
sub-range T-stretch-only spelling (`c.1210-12_1210-6T[7]`), when only that portion's variability was
reported. Consistent with the sequenced pattern; no new rule. On `NM_000492.3`, so parse-only.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_000492.3:c.1210-33_1210-6GT[11]T[6]` | recommended | — | the full intron-9 mixed repeat (parse-only here) |
| `NM_000492.3:c.1210-12_1210-6T[7]` | recommended | — | the T-stretch-only sub-range spelling (parse-only here) |

## `repeated.md:85-86` — a standalone mixed-repeat Example

> - **`NC_000012.11:g.112036755_112036823CTG[9]TTG[1]CTG[13]`**<br>
>   a complex repeated sequence from position `g.112036755` to `g.112036823` on chromosome 12 with first a `CTG` unit present in 9 copies, then a `TTG` unit present in 1 copy and then a `CTG` unit present in 13 copies.

Ferro: this restates the mixed-repeat rule of `repeated.md:18`, on `NC_000012.11` — the row is shown
there. Same open partition-uniqueness question, no ferro obligation adjudicated, no duplicate row.

## `repeated.md:88-90` — differing genomic and coding descriptions, minus strand

> `NC_000001.11:g.57367047_57367121ATAAA[15]` and `NM_021080.3:c.-136-75952_-136-75878ATTTT[15]` describe the same repeat allele in intron 3 of the _DAB1_ gene.<br>
> **NOTE**: based on the **3' rule** and the transcriptional orientation of the gene (minus strand), the description of the repeat units differs.

Ferro: `ATAAA` and `ATTTT` are reverse complements, so this is the same 3'rule-plus-strand
mechanism, demonstrated genome-versus-transcript rather than literature-versus-transcript. The two
spellings must denote the same underlying allele. On foreign accessions, so parse-only.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000001.11:g.57367047_57367121ATAAA[15]` | recommended | — | the genomic (plus-strand) unit `ATAAA` (parse-only here) |
| `NM_021080.3:c.-136-75952_-136-75878ATTTT[15]` | recommended | — | the coding (minus-strand) unit `ATTTT` — the reverse complement, the 3'rule rotation on the transcript strand (parse-only here); denotes the same allele as the row above |

## `repeated.md:94-101` — Discussion: CFTR intron 9, the `GT`-versus-`TG` phase

> !!! note "Intron 9 of the _CFTR_ gene ends with the sequence `...tgtgtgtgtgtttttttaacag`. Both the `TG` and `T` stretches are variable in length (from 9 to 13 and 5 to 9, respectively). The reference sequence has 11 `TG` copies and 7 `T`s. Is it correct to describe an allele as `c.1210-14TG[13]T[5]` or for the T stretch as `c.1210-6T[5]`?"

Ferro: two points from the worked answer. The 3'rule selects the unit's **phase/register**, not just
a position: the raw run reads as a `TG` dinucleotide by naive grouping, but the 3'rule forces the
`GT` rotation (whose boundary sits maximally 3' against the adjacent `T`-stretch). And the uncertain
copy-number bracket `[(9_13)]` / `[(4_8)]` applies the general parenthetical-uncertainty convention
to the copy-number element — a shape shown only here. Both are compositions of already-general
rules; no new obligation, no Why block. On `NM_000492.3`, so parse-only.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_000492.3:c.1210-33_1210-6GT[11]T[6]` | recommended | — | the reference allele in the `GT` phase the 3'rule selects (parse-only here) |
| `NM_000492.3:c.1210-12_1210-6T[7]` | recommended | — | the T-stretch-only reference spelling (parse-only here) |
| `NM_000492.3:c.1210-33_1210-6GT[(9_13)]T[(4_8)]` | recommended | — | the population-variability form — uncertain copy-number ranges on both units (parse-only here) |
| `NM_000492.3:c.1210-12_1210-6T[(5_9)]` | recommended | — | the T-stretch population-variability form (parse-only here) |
