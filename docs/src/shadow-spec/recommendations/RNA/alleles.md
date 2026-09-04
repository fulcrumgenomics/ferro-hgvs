# Alleles — ferro's reading

ferro's reading of `RNA/alleles.md`. The rules are HGVS's; ferro's job is to produce the form the
recommendations prefer. Verdicts describe **ferro's output**:

- **recommended** — ferro's output is the form the recommendations prefer (whether the input was
  already that form, or ferro normalized it there).
- **conformant** — ferro's output is valid HGVS but not *yet* the recommended form — a ferro
  limitation or a deliberate maintainer house choice among conformant forms, with a tracking
  issue where one exists.
- **refused** — the input is not valid HGVS; ferro rejects it in strict mode (correct behavior).
- **bug** — ferro's output is not valid HGVS (a defect). None on this page — but two kinds of
  `conformant` row here are not mere limitations, and the verdict word undersells both. Three rows
  at `:5` are inputs ferro **refuses in strict mode at normalize** (`W5002`), which the checked
  `refused` verdict — a strict *parse* check — cannot carry, so they are recorded on what the
  lenient tier emits. Two rows at `:36-52` / `:82-85` re-spell the spec's own compact form to an
  expanded one; the output re-parses, which is what the verdict measures, and the defect is filed as
  [#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214).

Each **Why** block is transcluded from the ruling ledger — the record's own one-line summary,
rendered here and linked to its full entry in
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md).
The reasoning lives once, in the ledger; it is never re-typed here.

**No ledger record cites any `RNA/alleles.md` clause** — the whole document is a record-level gap,
so every section below carries **no Why block**. The reading is CONFIRM-by-inspection against the
spec text and the shipped code, not an adjudicated ruling. The one cross-axis record that touches
allele member geometry, `conflicting-member-geometry-refusal-scope`, governs **`DNA/alleles.md:5`**
and its guard (`the_three_conflicting_geometries_are_refused_in_strict_mode_on_every_axis`)
enumerates only `c.`/`g.` shapes — so by the repo's own jurisdiction rule (an `r.` claim needs a
clause under `RNA/`) it carries no RNA authority. Where it is relevant it is documented in a Note,
not a Why block; see `alleles.md:5`.

How ferro treats an allele, which every section below rests on: the members are parsed into an
`AllelePhase` — `Cis` (`[a;b]`), `Trans` (`[a];[b]`), `Unknown` (`a(;)b`), `Products` (`,`) — and
**each member is normalized in its own frame**. Every merge/collapse pass is gated on
`phase != Cis`, so a **cis** allele's members may merge when they sit adjacent, while **trans**,
**unknown** and **Products** members are *only* individually normalized — never merged, never
reordered. An `r.` allele therefore round-trips member-wise. Measured on the executable rows
below: a shiftable member 3'-shifts inside every phase, including a predicted `[(…)]` wrap and a
`,` Products allele, and the join is preserved in each.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
examples sit on `LRG_199t1` and `NM_004006.2`, neither of which the slice carries, so those rows
are parse-only (`—`) — ferro cannot read their bases here. The `NM_004006.3` base facts the
executable rows rely on (`r.75_77` is `aaa`, so `r.76del` → `r.77del` and `r.77del` is a fixed
point; `r.6_8` is `uug` with `r.9` = `g`, a fixed point; `r.124_129del` is a fixed point) are the
same ones established on `RNA/deletion.md`. One further measured fact: `r.345del` 3'-shifts to
`r.346del`, which is what makes the syntax table's own `NM_004006.3` example executable at `:17`.

## `alleles.md:5` — definition, and the member-geometry refusal

> Allele: a series of variants in a transcript from one chromosome.

Ferro: an allele is a series of variants on one transcript from one chromosome. Two members that
claim **intersecting reference territory** — nested, overlapping, or two insertions at one point —
contradict "a series of variants" and are refused. The refusal is coordinate-based (it compares
member spans, no reference bases needed) and **is axis-general in behaviour**: on the committed
slice each `r.` row below fails strict-mode normalization with the same `W5002
OverlapConflictingEdits` message as its `c.` twin (`NM_004006.3:c.[9_12del;10_11del]`,
`c.[9_12del;11_14del]`, `c.[8_9insAT;8_9insTC]` — measured, byte-identical apart from the axis).

**The stage matters for how the rows are recorded.** The refusal fires at **normalize**, not at
parse: strict *parse* accepts all three, and the check runs on the parsed allele's member spans
before any merge. This page's `refused` verdict is checked at strict **parse**, so it cannot carry
a normalize-time refusal; the three rows are therefore recorded on what the lenient tier — the tier
`Normalizes to` is checked against — emits, which is the input unchanged, carrying the `W5002`
warning. Read `conformant` there as "strict refuses; lenient passes through", not as "ferro accepts
intersecting members".

**No ledger record scopes this to `r.`.** `conflicting-member-geometry-refusal-scope` decides the
refusal but cites only `DNA/alleles.md:5`, and its guard enumerates `RefShape::CodingSingleExon`,
`Genomic`, `CodingMultiExon` — c. and g. only, no `r.` case. `RNA/alleles.md:5` is a verbatim twin
of that definition, but a `DNA/` citation cannot carry an `r.` ruling, so the `r.` refusal is
measured behaviour that is currently un-authorised and un-guarded (worksheet AL2 proposes co-citing
`RNA/alleles.md:5` and adding an `r.` guard row — a proposal, not yet in the ledger; the extension
is not invented here).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[6_8del;77del]` | recommended | self | a well-formed cis allele — two disjoint members, each a fixed point; kept as a series of variants |
| `NM_004006.3:r.[9_12del;10_11del]` | conformant | self | nested members — the second is inside the first; intersecting territory violates the definition. **Strict refuses at normalize**: `NM_004006.3:r.9_12 has 2 coincident cis-allele edits (del, del); HGVS spec defines no canonical form for this case (OverlapConflictingEdits / W5002)`. Lenient emits the input unchanged with that warning, which is what this row records. No `r.`-jurisdiction ledger record or guard covers the refusal (worksheet AL2) |
| `NM_004006.3:r.[9_12del;11_14del]` | conformant | self | overlapping members — same `W5002` refusal in strict, same lenient pass-through |
| `NM_004006.3:r.[8_9insau;8_9insuc]` | conformant | self | two insertions at one point (`r.8_9`) — coincident insertions; strict refuses (`… has 2 coincident cis-allele edits (ins, ins) … W5002`), lenient passes through. The RNA payloads (`au`, `uc`) are valid symbols; the geometry is what is refused |

## `alleles.md:15-16` — DNA-level reporting, and diploidy

> - all variants **should be** described on the DNA level; descriptions on the RNA and/or protein level may be given in addition.
> - humans are diploid organisms and have **two alleles** at each genetic locus, with one allele inherited from each parent.

Ferro: both notes are descriptive. `:15` is authoring guidance about which level(s) to report on,
not a constraint on how a given `r.` allele, once written, is normalized — a lone `r.` description
stands on its own. `:16` is background biology. Neither is anything a normalizer enforces.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[6_8del;77del]` | recommended | self | a lone RNA-level allele is fully valid; ferro does not require an accompanying `c.`/`g.` form |

## `alleles.md:17` — cis: `r.[variant1;variant2]`

> - when two variants are identified in a transcript that derive from **one chromosome** (in cis), this should be described as `r.[variant1`<code class="spot1">;</code>`variant2]`.

Ferro: parses the `;`-joined bracket to `AllelePhase::Cis`. Members are normalized individually and,
because cis is the one phase **not** gated out of merging, adjacent members may collapse to a single
`delins` (that member-level merge behaviour is governed on `RNA/delins.md`, out of this file's
scope). Members two or more nucleotides apart stay individual, and the `;` join is preserved.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[6_8del;77del]` | recommended | self | two disjoint members, each a fixed point — individually normalized, no merge, `;` preserved |
| `NM_004006.3:r.[76del;124_129del]` | recommended | `NM_004006.3:r.[77del;124_129del]` | cis members individually normalized: `r.76del` 3'-shifts to `r.77del` over the `aaa` run, the far member is untouched, no merge |
| `NM_004006.3:r.[123c>a;345del]` | recommended | `NM_004006.3:r.[123c>a;346del]` | the syntax table's own cis example (`syntax.yaml`, `rna.alleles`), which happens to sit on the slice's accession: the substitution is a fixed point (`r.123` is `c`), the deletion 3'-shifts one base |
| `LRG_199t1:r.[76a>u;103del]` | recommended | — | the spec's own cis example (`:28-30`); parse-only here |

## `alleles.md:18` — trans: `r.[variant1];[variant2]`

> - when two variants are identified in transcripts that derive from **different chromosomes** (in trans), this should be described as `r.[variant1]`<code class="spot1">;</code>`[variant2]`.

Ferro: parses the `];[` join to `AllelePhase::Trans`. Trans is two chromosomes — two distinct
variants — so it is gated out of every merge pass: members are only individually normalized, never
merged against each other and never reordered. A distinct object from cis, and preserved as such.

The recommended spelling lists each allele's own variant, `r.[a];[b]`. Padding each allele's bracket
with the *other* allele's positions as `=` fillers — `r.[a;p=];[p=;b]` — is the cross-spelled form
the DNA twin marks `class="invalid"` by name (`DNA/alleles.md:19`, restated at `:49`: the "not
changed" indication is used only when one variant was identified). Ferro rejects it **at parse**,
in every mode: `trans_has_cross_spelled_identity` (`src/hgvs/parser/variant.rs`, on the trans
production) refuses a position that carries a concrete edit in one bracket and an `=` inside a
*multi-member* bracket elsewhere. A single-member `[77=]` beside `[77del]` is the sanctioned form
and is accepted (`:47-49` below); so is an `=` at a position nothing else touches
(`r.[77del;78=];[78=;6_8del]` parses). This is not the `:5` geometry refusal — an earlier draft
of this page attributed it there — and the clause it enforces is a `DNA/` one, so on `r.` the
rejection rests on the DNA twin's NOTE with no `RNA/alleles.md` counterpart. The error surfaces as a
generic parse failure (`Parse error at position 12 … input: "]"`) that names neither the clause nor
the cross-spelled position (minor gap).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[77del];[6_8del]` | recommended | self | the recommended trans form — each allele carries its own variant, both fixed points, preserved unmerged |
| `NM_004006.3:r.[76del];[77del]` | recommended | `NM_004006.3:r.[77del];[77del]` | each allele's member individually normalized (`r.76del` → `r.77del`); the two brackets are **not** merged into one (trans is not cis) |
| `NM_004006.3:r.[77del;77=];[77=;6_8del]` | refused | — | the cross-spelled form — the recommended spelling is `r.[77del];[6_8del]`. Rejected at parse by `trans_has_cross_spelled_identity`: `r.77` carries a concrete edit in one bracket and an `=` filler inside a multi-member bracket in the other, the shape `DNA/alleles.md:19` calls not correct. Authority is the DNA twin's; `RNA/alleles.md` states no counterpart |
| `LRG_199t1:r.[76a>u];[103del]` | recommended | — | the spec's own heterozygous trans example (`:37-40`); parse-only here |

## `alleles.md:19-20` — unknown phase: `variant1(;)variant2`

> - when two variants are identified in a transcript, but when it is **not known** whether these derive from one chromosome (in cis) or from different chromosomes (in trans), this should be described as `variant1`<code class="spot1">(;)</code>`variant2`, i.e. without using `[ ]`.<br>
>   **NOTE**: it is recommended to determine whether the changes are in the same transcript or not.

Ferro: parses the bracket-less `(;)` join to `AllelePhase::Unknown` and keeps the compact, no-`[ ]`
form. Unknown phase is gated out of merging, so members are individually normalized and never
merged or reordered. The NOTE is advisory — a recommendation to the author to resolve the phase,
not something a normalizer acts on.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.6_8del(;)77del` | recommended | self | unknown phase, no brackets — two fixed-point members, individually normalized, `(;)` and order preserved |
| `NM_004006.3:r.76del(;)124_129del` | recommended | `NM_004006.3:r.77del(;)124_129del` | the first member 3'-shifts (`r.76del` → `r.77del`); no brackets minted, members not merged |
| `NM_004006.2:r.76a>u(;)103del` | recommended | — | the spec's own unknown-phase example (`:54-57`); parse-only here |

## `alleles.md:21` — comma / Products: `r.[variant1,variant2]`  (RNA/protein-unique)

> - when two variants are identified in two different transcripts that derive from **one variant** on the DNA level, the variants are separated using a `,`; `r.[variant1`<code class="spot1">,</code>`variant2]`.

Ferro: the `,` (Products) allele has **no analogue on the DNA axis**. Its members are downstream
**products of a single DNA event** (alternative splicing) — different transcripts, not co-occurrence
phases — so they must never be partition-merged against each other. Ferro models it as
`AllelePhase::Products`, restricts the form to the `r.`/`p.` axes (`c.[a,b]` / `g.[a,b]` are
rejected), forbids mixing `,` with `;` in one bracket, normalizes each member in its own frame, and
— because Products is not Cis — never merges or reorders them. The cross-reference `general.md:80`
gives the bare form `r.[123a>u,122_154del]`.

**There is no ledger record for this form.** No record cites `RNA/alleles.md:21` or `general.md:80`,
and none adjudicates Products-allele normalization semantics, so this section carries **no Why
block**; the behaviour is a defensible but unrecorded ferro reading (worksheet AL1 proposes a
`decided` house-choice record pinning: `,` members are distinct-transcript products, never merged or
re-partitioned; restricted to `r.`/`p.`; each member normalized in its own frame — a proposal, not
yet filed).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.[897u>g,832_960del]` | recommended | — | the spec's own Products example (`:60-61`); the two members derive from one DNA variant (`LRG_199t1:c.897T>G`); parse-only here |
| `NM_004006.3:r.[6_8del,124_129del]` | recommended | self | executable Products twin — two distinct-transcript products, each a fixed point; `AllelePhase::Products`, kept unmerged and unreordered (matches `comma_products_allele.rs`). Illustrative positions; `general.md:80`'s own string is the bare `r.[123a>u,122_154del]` |
| `NM_004006.3:r.[76del,124_129del]` | recommended | `NM_004006.3:r.[77del,124_129del]` | each product normalized in its own frame: the first 3'-shifts to `r.77del`, the second is untouched, the `,` join and order are preserved |
| `NM_004006.3:c.[123del,124del]` | refused | — | the `,` form is `r.`/`p.`-only — a `c.` Products allele is rejected at parse (`find_products_bracket` restricts the comma form to the RNA/protein axes) |

## `alleles.md:27-34` — examples: variants on one allele (cis)

> - **variants on one allele**
>     - **`LRG_199t1:r.[76a>u;103del]`**<br>
>       one transcript contains two different changes, `r.76a>u` and `r.103del`.
>       The variants are found in _cis_.
>
>     - **`LRG_199t1:r.[(578c>u;1339a>g;1680del)]`**<br>
>       one transcript contains three different predicted changes, `r.(578c>u)`, `r.(1339a>g)`, and `r.(1680del)`.
>       The variants are found in _cis_.

Ferro: the plain cis allele (ex 1) is `AllelePhase::Cis` with members individually normalized and
preserved. The second example wraps the whole allele in `[(...)]` — a **predicted** (uncertain) cis
allele; ferro preserves the uncertainty flag and normalizes members inside it, without collapsing
the prediction to a concrete assertion.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.[76a>u;103del]` | recommended | — | example 1, a plain cis allele; parse-only here |
| `LRG_199t1:r.[(578c>u;1339a>g;1680del)]` | recommended | — | example 2, a predicted `[(…)]` cis allele — uncertainty flag preserved; parse-only here |
| `NM_004006.3:r.[(6_8del;77del)]` | recommended | self | executable twin of the predicted form — the `[(…)]` wrap is preserved, members individually normalized (both fixed points) |
| `NM_004006.3:r.[(76del;124_129del)]` | recommended | `NM_004006.3:r.[(77del;124_129del)]` | the wrap is preserved **and** the members inside it are normalized: `r.76del` 3'-shifts to `r.77del` without the prediction being collapsed |

## `alleles.md:36-52` — examples: variants on two alleles (trans), and `[76=]` vs `[=]`

> - **variants on two alleles**
>     - **`LRG_199t1:r.[76a>u];[103del]`**<br>
>       the two transcript alleles each contain a different change, `r.76a>u` and `r.103del`.
>       A **heterozygous** case (compound heterozygote, e.g., in a recessive disease).
>       The variants are found in _trans_.
>
>     - **`NM_004006.2:r.[76a>u];[76a>u]`**<br>
>       both transcript alleles contain the same variant, `r.76a>u`.
>       A **homozygous** case (e.g., in a recessive disease).<br>
>       **NOTE**: `LRG_199t1:r.76a>u(;)(76a>u)` indicates analysis detects one variant (`r.76a>u`), suggesting both transcript alleles contain this variant, but it can not be excluded the other allele is deleted or not expressed.
>
>     - **`LRG_199t1:r.[76a>u];[76=]`**<br>
>       one transcript allele contains a variant, `r.76a>u`, the other transcript allele contains at position `r.76` the reference sequence, `r.76=` (is **wild-type**).<br>
>       **NOTE**: the description `r.[76a>u];[=]`, containing `r.76a>u` and `r.=`, is different since it indicates the entire coding RNA reference sequence was analysed and the only variant identified was `r.76a>u` (on one allele).
>
>     - **`NM_004006.2:r.[76a>u];[?]`**<br>
>       one transcript allele contains a variant, `r.76a>u`, while a variant in the other transcript allele is expected but not yet identified (`r.?`) (e.g., in individuals affected by a recessive disease).

Ferro: four trans shapes, all preserved unmerged — heterozygous (two different members), homozygous
(the same member on both alleles), one allele wild-type at a named position (`[76=]`), and one
allele's variant not yet identified (`[?]`). The `:47-49` NOTE draws a **meaning** distinction the
parser must not erase: `[76=]` asserts wild-type **at position 76**, while `[=]` asserts the **whole
sequence** was analysed and unchanged. Ferro preserves each exactly and never collapses one into the
other.

**The `[?]` shape is re-spelled.** Ferro keeps the `[?]` member (it lowers to the dedicated
`UnknownAllele` sentinel) but renders the whole allele in the **expanded** per-bracket-accession
form, `[NM_004006.3:r.77del];[?]` — the accession moves inside the first bracket. That happens at
parse → `Display`, before normalization touches it, in every mode, on `c.` as well as `r.`, and for
`[0]` too (`:82-85`); a `[76=]` or `[=]` second member stays compact. The output is valid HGVS and
re-parses — the expanded form is the one `DNA/alleles.md:110` recommends for a *two-accession*
trans allele — but it is not the shape of the spec's own same-accession example at `:51`, and
`trans_compact_anchor` (`src/hgvs/variant.rs`) states in its own doc that `[0]`/`[?]` members are
meant to keep the compact form. Filed as
[#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214); the existing pins
(`allele_trans_phase.rs::test_trans_n_with_unknown_allele_round_trip`,
`issue_396_trans_allele_consolidation.rs::rna_trans_with_unknown_allele_roundtrip`) assert only a
round-trip or a `contains(";[?]")`, which both renderings satisfy.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.[76a>u];[103del]` | recommended | — | heterozygous trans; parse-only here |
| `NM_004006.2:r.[76a>u];[76a>u]` | recommended | — | homozygous trans — both alleles carry the same member, preserved (not merged); parse-only here |
| `NM_004006.3:r.[77del];[77del]` | recommended | self | executable homozygous twin — identical members on two alleles, kept as two brackets |
| `LRG_199t1:r.[76a>u];[76=]` | recommended | — | `[76=]` = position-wild-type member, preserved; parse-only here |
| `NM_004006.3:r.[77del];[76=]` | recommended | self | executable twin — the `[76=]` position-reference member is preserved, not collapsed to `[=]` |
| `NM_004006.3:r.[76del];[76=]` | conformant | `NM_004006.3:r.[77del];[76=]` | the variant 3'-shifts to `r.77` while the `[76=]` filler stays where it was authored — an identity member has nothing to shift. Valid, and still true of the second allele; but the spec's examples always pair the filler with the variant's own position, so `r.[77del];[76=]` and the authored-at-77 `r.[77del];[77=]` are two spellings of one observation that do not converge. Not adjudicated and no issue filed — flagged for a ruling on whether the filler should be re-anchored to the shifted member |
| `NM_004006.3:r.[77del];[=]` | recommended | self | the whole-sequence `[=]` member — a *different* meaning from `[76=]`, preserved as written |
| `NM_004006.2:r.[76a>u];[?]` | recommended | — | `[?]` = unknown second allele, preserved; parse-only here |
| `NM_004006.3:r.[77del];[?]` | conformant | `[NM_004006.3:r.77del];[?]` | executable twin — the `[?]` member is kept, but the allele is re-spelled from the compact form the spec's example uses to the expanded per-bracket-accession form. Valid, re-parses; not the recommended shape ([#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214)) |

## `alleles.md:54-57` — alleles not certain (unknown phase)

> - **alleles not certain**
>     - **`NM_004006.2:r.76a>u(;)103del`**<br>
>       two variants are found in a transcript, `r.76a>u` and `r.103del`, but it is not known whether they derive from the same or from different transcript alleles (chromosomes).<br>
>       **NOTE**: when it is not known on which allele a variant is, allele brackets should not be used.

Ferro: the worked example of the `:19-20` unknown-phase rule — the bracket-less `(;)` form is kept
bracket-less, members individually normalized, phase not inferred.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:r.76a>u(;)103del` | recommended | — | the spec's own unknown-phase example; parse-only here |
| `NM_004006.3:r.6_8del(;)77del` | recommended | self | executable twin — no brackets, members preserved |

## `alleles.md:59-61` — one allele, two transcripts (Products)

> - **one allele, two transcripts**
>     - **`LRG_199t1:r.[897u>g,832_960del]`**<br>
>       two different transcripts, `r.897u>g` and `r.832_960del`, derive from one variant (`LRG_199t1:c.897T>G` on the DNA level).

Ferro: the worked example of the `:21` Products rule — the two `,`-separated members are different
transcripts derived from **one** DNA variant, so they are `AllelePhase::Products` and are never
merged or reordered against each other, each normalized in its own frame. (Adjudication is unrecorded
— see `alleles.md:21`; no Why block.)

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.[897u>g,832_960del]` | recommended | — | the spec's own Products example; members from one DNA event (`c.897T>G`), kept distinct; parse-only here |
| `NM_004006.3:r.[6_8del,124_129del]` | recommended | self | executable Products twin — two distinct-transcript products, each a fixed point, unmerged |

## `alleles.md:65-69` — Discussion: the retracted `+` separator

> !!! note "Was originally the recommendation to use the format <code class="invalid">[r.76a>c+r.83g>c]</code>?"
>
>     Indeed, originally [den Dunnen and Antonarakis, 2000](http://dx.doi.org/10.1002/%28SICI%291098-1004%28200001%2915:1%3c7::AID-HUMU4%3e3.0.CO;2-N) the suggestion was to describe two changes in a transcript from one chromosome as <code class="invalid">[r.76a>c+r.83g>c]</code>, i.e. using a "+"-character to separate the two changes, while an earlier publication suggested to use a ";" (<code class="invalid">[r.76a>c;r.83g>c]</code> [(Antonarakis and the Nomenclature Working Group, 1998](http://dx.doi.org/10.1002/%28SICI%291098-1004%281998%2911:1%3c1::AID-HUMU1%3e3.0.CO;2-O)).
>     To prevent confusion with older publications, to improve overall consistency, and to keep descriptions as short as possible, the 2000 proposal was retracted.
>     The recommended format is `r.[76a>c;83g>c]`.

Ferro: the `+` separator inside allele brackets is `class="invalid"` — a retracted 2000 proposal.
The `;` is the recommended cis separator; ferro does not emit `+` and rejects it at parse.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[6_8del+77del]` | refused | — | the retracted `+` separator — not valid HGVS; rejected at parse (the grammar stops at `+77del`). The recommended cis form uses `;` (`r.[6_8del;77del]`) |

## `alleles.md:71-76` — Discussion: recording variant combinations in recessive disease

> !!! note "In recessive diseases, is it important I show in which combination variants were found?"
>
>     When in one individual you find more than one variant, it is essential that you clearly indicate on which transcript allele(s) variant(s) were found.
>
>     - disease severity will depend on the combination of variants found;
>     - in recessive disease, when two variants are in one transcript, an individual is a carrier or you might not have found the variant on transcripts from the second allele.

Ferro: advisory guidance to the author to record the cis/trans combination clearly. It picks no
canonical form and constrains no normalization — descriptive, nothing for the normalizer to enforce.
The cis/trans/unknown machinery that lets the author express the combination is covered at
`alleles.md:17`, `:18` and `:19-20`.

## `alleles.md:78-80` — Discussion: the empty allele `[]`

> !!! note "I find the notation <code class="invalid">r.[76a>c]</code> without describing the second transcript allele misleading; not enough researchers know this refers to only one of the two transcripts present. Would using <code class="invalid">r.[76a>c];[]</code> be OK?"
>
>     No, the recommended description is `r.76[a>c];[=]`, i.e. `r.76=` for "no change" at position `r.76` on the second transcript.

Ferro: the empty second allele `[]` is `class="invalid"` — refused at parse. The spec's recommended
answer is the position-factored trans form `r.76[a>c];[=]`, which shares the position out of the
brackets. Ferro neither produces nor **accepts** that factored spelling — it rejects it at parse
and emits the fully-bracketed `r.[76a>c];[76=]`-style form. That is the DNA twin's reading, which
forbids the analogous factoring by name: `DNA/alleles.md:54` marks `c.2376[G>C];[G>C]`
`class="invalid"` ("it is not allowed to shorten this"), and the RNA syntax table (`syntax.yaml`,
`rna.alleles`) publishes no factored form. So the RNA Q&A's answer contradicts both its DNA
parallel and the grammar — a spec-internal inconsistency for the upstream-filing lane, not
ferro-actionable (worksheet marks it DESCRIPTIVE). Ferro's fully-bracketed output is the
recommended form on the reading the rest of the spec supports.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[77del];[]` | refused | — | the empty allele `[]` is `class="invalid"` — rejected at parse. Recommended: name the second allele's state, e.g. `r.[77del];[77=]` |
| `NM_004006.3:r.76[a>u];[=]` | refused | — | the Q&A's own factored answer — rejected at parse (`Parse error at position 12 … "[a>u];[=]"`). The factoring is not in the RNA syntax table and its DNA analogue is `class="invalid"` at `DNA/alleles.md:54`; the bracketed `r.[76a>u];[76=]` is what ferro emits and accepts |
| `NM_004006.3:r.[77del];[77=]` | recommended | self | the fully-bracketed form ferro emits for "variant on one allele, wild-type at that position on the other" |

## `alleles.md:82-85` — Discussion: X-chromosome male, `r.0` for an absent transcript

> !!! note "How should I describe the variants detected in males and females for a transcript from the X-chromosome?"
>
>     In **females**, the description is straightforward, like `LRG_199t1:r.[76a>c];[76=]`.
>     In **males**, there is no transcript from the second allele (X-chromosome), which can be described as `LRG_199t1:r.[76a>c];[0]`, i.e. using `r.0` to indicate the absence of a transcript from the second X-chromosome.

Ferro: `r.0` denotes an absent transcript (no expression from the second X allele in a male). Ferro
accepts the `[0]` member and lowers it to the dedicated `NullAllele` sentinel on the `r.` axis
(pinned by `issue_277_trans_allele_protein_zero.rs::non_protein_compact_trans_allele_zero_stays_null_allele`
on `NR_002196.2:r.[100a>g];[0]`, and by
`issue_396_trans_allele_consolidation.rs::rna_trans_with_null_allele_roundtrip`). The female case is
the ordinary `[variant];[76=]` trans form from `:47-49`. The male case is subject to the same
re-spelling as `[?]` at `:36-52`: the member is preserved, but the allele is rendered in the
expanded per-bracket-accession form ([#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214)).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.[76a>c];[76=]` | recommended | — | the female case — an ordinary position-wild-type trans allele; parse-only here |
| `LRG_199t1:r.[76a>c];[0]` | recommended | — | the male case — `[0]` = absent transcript, preserved; parse-only here |
| `NM_004006.3:r.[77del];[0]` | conformant | `[NM_004006.3:r.77del];[0]` | executable twin — the `[0]` absent-transcript member is accepted and kept unmerged, but the allele is re-spelled from the compact form the spec writes to the expanded one. Valid, re-parses; not the recommended shape ([#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214)) |

See also → `DNA/alleles.md:5` (the twin allele definition, the geometry-refusal record's actual
jurisdiction), `DNA/alleles.md:19` and `:49` (the cross-spelled `=` form the `:18` parse guard
enforces), `DNA/alleles.md:54` (the factoring the `:78-80` Q&A recommends and the DNA twin forbids),
`DNA/alleles.md:110` (the two-accession expanded trans form ferro emits for `[?]`/`[0]`),
`RNA/deletion.md:18-20` (the `NM_004006.3` 3'-shift facts the executable members rely on),
`general.md:80` (the `,` Products separator and its bare example `r.[123a>u,122_154del]`).
