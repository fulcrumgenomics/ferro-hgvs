# Alleles — ferro's reading

ferro's reading of the HGVS **alleles** recommendations, clause by clause — each spelling with the
form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Alleles (`r.`)](../RNA/alleles.md).*

Two `conformant` verdicts on this page are not mere limitations. Three rows at `alleles.md:5` are
inputs ferro **refuses in strict mode at normalize** (`W5002`) — a normalize-time refusal the
checked `refused` verdict (a strict *parse* check) cannot carry, so they are recorded on what the
lenient tier emits instead. The `[?]` and `[0]` rows re-spell the spec's own compact form to an
expanded one; the output re-parses, but this is filed as
[#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214).

How ferro treats an allele, which every section below rests on: the members are parsed into an
`AllelePhase` — `Cis` (`[a;b]`), `Trans` (`[a];[b]`), `Unknown` (`a(;)b`) — and **each member is
normalized in its own frame**. Every merge/collapse pass early-returns unless `phase == Cis`, so a
**cis** allele's members may merge when they sit adjacent, while **trans** and **unknown** members
are *only* individually normalized — never merged, never reordered.

Executable rows use `NM_004006.3`, the one transcript in the committed slice, over the base facts
established on `substitution.md` (`c.76` is `A`, `c.79` is `G`, `c.80` is `C`, `c.123` is `C`).
Substitution members are fixed points, which keeps these rows decidable on the hermetic slice. The
spec's own examples sit on accessions the slice does not carry, so those rows are parse-only (`—`).

## `alleles.md:5` — definition, and the member-geometry refusal

> Allele: a series of variants on one chromosome.

Ferro: an allele is a series of variants on one chromosome. Two members claiming intersecting
reference territory — nested, overlapping, or two insertions at one point — contradict that and
are refused. The refusal is coordinate-based (member spans only, no reference bases needed).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[conflicting-member-geometry-refusal-scope](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two members of one allele that claim intersecting reference territory — nested, overlapping, or two insertions at one interbase — are refused, whatever edit types they render as.
>
> **[inversion-vs-a-mixed-member-competitor](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When a span replaced by its reverse complement competes with a description mixing lone substitutions and multi-column members, ferro writes it as one inv; both forms are conformant, so this is the project's choice among them.
<!-- why:END:conflicting-member-geometry-refusal-scope,inversion-vs-a-mixed-member-competitor -->

</details>

The geometry-refusal record's jurisdiction is DNA, governing this clause directly: the definition
("a series of variants on **one chromosome**") is what makes a nested/overlapping/coincident-
insertion pair undenotable, since there is no single chromosome sequence such a pair can be read
off. The inversion record cites this clause too, but only as background — its substantive ruling
rests on `inversion.md:5`.

The refusal fires at **normalize**, not at parse — strict parse accepts all three, and the check
runs on the parsed allele's member spans before any merge. This page's `refused` verdict is checked
at strict parse, so it cannot carry a normalize-time refusal; the three rows are recorded on what
the lenient tier emits instead — the input unchanged, carrying the `W5002` warning. Read
`conformant` there as "strict refuses; lenient passes through", not as ferro accepting intersecting
members.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | a well-formed cis allele — two disjoint members, each a fixed-point substitution; kept as a series of variants |
| `NM_004006.3:c.[9_12del;10_11del]` | conformant | self | nested members — the second is inside the first; intersecting territory violates the definition. Strict refuses at normalize (`OverlapConflictingEdits` / `W5002`); lenient emits the input unchanged with that warning, which is what this row records |
| `NM_004006.3:c.[9_12del;11_14del]` | conformant | self | overlapping members — same `W5002` refusal in strict, same lenient pass-through |
| `NM_004006.3:c.[8_9insAT;8_9insTC]` | conformant | self | two insertions at one point (`c.8_9`) — coincident insertions; strict refuses (`W5002`), lenient passes through. The payloads are valid; the geometry is what is refused |

## `alleles.md:15` — diploidy

> - humans are diploid organisms and have **two alleles** at each genetic locus, with one allele inherited from each parent.

Ferro: background biology — no syntactic or normalization rule to enforce. Descriptive.

## `alleles.md:16` — cis: `g.[variant1;variant2]`

> - when two variants are identified in a gene that are on **one chromosome** (in cis), this should be described as `g.[variant1`<code class="spot1">;</code>`variant2]`.

Ferro: parses the `;`-joined bracket to `AllelePhase::Cis`. Members are normalized individually,
and because cis is the one phase not gated out of merging, adjacent members may collapse to a
single delins (governed on `delins.md`). Members two or more nucleotides apart stay individual, and
the `;` join is preserved. The compact, single-prefix form is used only when every member shares
accession and coordinate type — the shape of every worked cis example here.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[junction-exit-wrapper-scope-in-a-mixed-allele](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An allele holds one member whose intronic offset ferro MANUFACTURED and one the author spelled intronic, both on one bare transcript. `checklist.md:20` wants a genomic reference for the manufactured member; `DNA/alleles.md:16` wants one reference factored out of the whole allele. Lift the wrapper to the whole description, or expand to per-member accessions?
<!-- why:END:junction-exit-wrapper-scope-in-a-mixed-allele -->

</details>

**The mixed-allele wrapper-scope question is open (undecided).** `alleles.md:16` describes the cis
shape only when a shared accession+axis prefix already exists. It says nothing about a member that
is an author-spelled bare intronic position paired with a sibling carrying an intronic offset ferro
manufactured, which needs a genomic wrapper — the compact form admits only one accession, so
lifting the wrapper would re-spell the authored sibling, while expanding to per-member accessions
abandons `:16`'s factored form. Ferro declines today (the offset ships bare), which is the status
quo, not a decision. No executable row: the shape needs a cdot the hermetic slice does not carry.

### Examples — cis (`alleles.md:30-41`)

> - **`NM_004006.2:c.[2376G>C;3103del]`**<br>

Ferro: ordinary same-accession cis alleles (plus two single-variant repeat notations at `:37` and
`:40`, cis framing rather than the cis-bracket grammar). All parse; the foreign-accession ones are
parse-only, except `alleles.md:40`'s open-ended repeat, whose disclosed render divergence needs no
reference bases and so executes.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.[2376G>C;3103del]` | recommended | — | the two-change cis allele; foreign accession — parse-only here |
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | executable cis twin — two disjoint fixed-point substitutions, `;` preserved, no merge (two nucleotides-plus apart) |
| `NC_000023.10:g.[30683643A>G;33038273T>G]` | recommended | — | one accession, two genes (`GK`, `DMD`); a gene boundary is metadata, not grammar (parse-only) |
| `NC_000003.12:g.63912687AGC[(50_60)]` | recommended | — | a single-variant repeat with a bounded uncertain copy-number range (parse-only) |
| `NC_000003.12:g.63912687AGC[(60_?)]` | conformant | `NC_000003.12:g.63912687AGC[60_?]` | open-ended uncertain repeat count. Ferro **drops the uncertainty parentheses** from the range — a deliberate render choice. The re-spelling needs no reference bases, so this executes despite the foreign accession; it re-parses as valid HGVS but is not the spec's stated form. Open finding: no ledger record or issue filed |

## `alleles.md:17` — trans: `g.[variant1];[variant2]`

> - when two variants are identified in a gene that are **on different chromosomes** (in trans), this should be described as `g.[variant1]`<code class="spot1">;</code>`[variant2]`.

Ferro: parses the `];[` join to `AllelePhase::Trans`. Trans is two chromosomes, so it is gated out
of every merge pass — members are only individually normalized, never merged or reordered. The
compact single-prefix form is used only when every member shares accession and coordinate type; a
`[?]` or `[0]` marker member forces the expanded per-bracket-accession form instead (see below).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.[2376G>C];[3103del]` | recommended | — | the recommended trans form — each allele carries its own variant. Foreign accession — parse-only here |
| `NM_004006.3:c.[76A>G];[123C>T]` | recommended | self | executable heterozygous twin — two fixed-point substitutions on two alleles, kept as two brackets (trans is not cis) |

### Examples — trans (`alleles.md:45-66`)

> - **`NM_004006.2:c.[2376G>C];[2376G>C]`**<br>

Ferro: the trans example shapes — heterozygous (two different members), homozygous (the same
member on both alleles, not shortened), a reference-allele arm (`[2376=]`, position-wild-type, kept
distinct from the whole-sequence `[=]`), and an unknown second allele (`[?]`). The two invalid
shortenings of the homozygous case are `class="invalid"` and refused at parse.

The `alleles.md:59` NOTE draws a meaning distinction the parser must not erase: `c.[2376=]` asserts
wild-type **at position 2376**, while `c.[=]` asserts the **whole** coding DNA reference sequence
was analysed and unchanged. Ferro preserves each exactly and never collapses one into the other.

**The `[?]` shape is re-spelled.** Ferro keeps the `[?]` member but renders the whole allele in the
expanded per-bracket-accession form — the accession moves inside the first bracket, in every mode,
on the `c.` axis as documented on the RNA twin, and for `[0]` too. A `[2376=]` or `[=]` second
member stays compact. The output is valid HGVS and re-parses — the expanded form is what
`alleles.md:110` recommends for a *two-accession* trans allele — but it is not the shape of the
spec's same-accession example. Filed as
[#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G];[76A>G]` | recommended | self | executable homozygous twin — identical members on two alleles, kept as two brackets, not shortened |
| `NM_004006.3:c.76[A>G];[A>G]` | refused | — | the factored-prefix shortening (`class="invalid"`) — rejected at parse |
| `NM_004006.3:c.76A>G[];[]` | refused | — | the empty-arm shortening (`class="invalid"`) — rejected at parse |
| `NM_004006.2:c.[296T>G;476T>C;1083A>C];[296T>G;1083A>C]` | recommended | — | three-variant-cis vs two-variant-cis trans compound (parse-only) |
| `NM_004006.3:c.[76A>G];[76=]` | recommended | self | `c.76=` = position-wild-type member (the other allele carries the reference at `c.76`), preserved, not collapsed to `[=]` |
| `NM_004006.3:c.[76A>G];[=]` | recommended | self | the whole-sequence `[=]` member — a *different* meaning from `[76=]`, preserved as written |
| `NM_004006.3:c.[76A>G;123=];[76=;123del]` | refused | — | the cross-spelled form the `:18-19` prohibition names — rejected at parse; see `alleles.md:18-19` |
| `NM_004006.2:c.[2376G>C];[?]` | recommended | — | `[?]` = unknown second allele, preserved (parse-only) |
| `NM_004006.3:c.[76A>G];[?]` | conformant | `[NM_004006.3:c.76A>G];[?]` | executable twin — the `[?]` member is kept, but the allele is re-spelled from the compact form to the expanded per-bracket-accession form. Valid, re-parses; not the recommended shape ([#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214)) |

## `alleles.md:18-19` — the cross-spelled `=` filler is not correct

> - using allele descriptions involving two or more different variants: you do not indicate the other allele does not contain the variant (unless one allele contains no variants at all).
>   `LRG_199t1:c.[2376G>C];[3103del]` is correct, <code class="invalid">LRG_199t1:c.[2376G>C;3103=];[2376=;3103del]</code> is not correct.

Ferro: this is the one clause on this page carrying real prose force ("is not correct",
`class="invalid"`) rather than a soft "should". The recommended spelling lists each allele's own
variant, `c.[a];[b]`. Padding each allele's bracket with the *other* allele's positions as `=`
fillers is the cross-spelled form marked invalid — the "not changed" indication is used only when
one variant was identified. Ferro rejects it at parse, in every mode, regardless of spelling
variant. A single-member `[3103=]` beside `[3103del]` is the sanctioned form (`alleles.md:59`) and
is accepted. This is an absolute prohibition, enforced unconditionally.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.[2376G>C];[3103del]` | recommended | — | the correct form — each allele carries its own variant (parse-only) |
| `NM_004006.2:c.[2376G>C;3103=];[2376=;3103del]` | refused | — | the `class="invalid"` cross-spelled form — rejected at parse |
| `NM_004006.3:c.[76A>G;123=];[76=;123del]` | refused | — | executable twin — `c.76` carries a concrete edit in one bracket and an `=` filler inside a multi-member bracket in the other; rejected at parse |

## `alleles.md:20-22` — unknown / uncertain phase: `variant1(;)variant2`

> - when two variants are identified in a gene, but when it is **not known** whether these are on one chromosome (in cis) or on different chromosomes (in trans), this should be described as `variant1`<code class="spot1">(;)</code>`variant2`, i.e. without using `[ ]`.<br>
>   **NOTE**: in the latest publication of the recommendations ([Den Dunnen et al. (2016)](http://onlinelibrary.wiley.com/doi/10.1002/humu.22981/pdf)), the example given is not correct.<br>
>   **NOTE**: it is recommended to determine whether the changes are on the same chromosome or not.

Ferro: parses the bracket-less `(;)` join to `AllelePhase::Unknown` and keeps the compact, no-`[ ]`
form. Unknown phase is gated out of merging — members are individually normalized and never merged
or reordered; brackets are refused around an unknown-phase pair. Both NOTEs are descriptive, not
something a normalizer acts on.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.2376G>C(;)3103del` | recommended | — | the spec's own unknown-phase example (parse-only) |
| `NM_004006.3:c.76A>G(;)123C>T` | recommended | self | executable twin — no brackets, two fixed-point members individually normalized, `(;)` and order preserved |
| `NM_004006.2:c.[296T>G;476T>C];[476T>C](;)1083A>C` | recommended | — | a three-way mixed cis/trans/unknown-phase compound (parse-only); ferro's phase-mixing grammar composes cis groups joined by `(;)` |

Two example shapes here are soft gaps. `alleles.md:75`'s `c.2376G>C(;)(2376G>C)` (the same variant
repeated in parens, "the other allele possibly deleted") is unverified to round-trip as its own
composite — the parenthesized/predicted member is a known DNA construct, so it is likely accepted
structurally, but unconfirmed. `alleles.md:81`'s four-variant compound has no string-for-string pin
either. Neither is asserted as a verdict row.

## `alleles.md:23` — different reference sequence types

> - descriptions combining variants based on different reference sequence types (e.g., <code class="invalid">c.[76A>C];g.[10091C>G]</code>) should not be used.

Ferro: the prose reads "should not be used" — a preference (rule 2), not an absolute prohibition. A
compact form with a mixed prefix inside one bracket is a genuine grammar violation and is refused
at parse, but the fully-expanded, syntactically distinct trans shape (each member carrying its own
accession and coordinate type) is accepted verbatim as spec-discouraged. The spec's own
illustration is not itself well-formed under the compact-form grammar, so ferro's split — compact
mixed prefix rejected, expanded cross-reference accepted-but-discouraged — is a sharper distinction
than the spec draws.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;g.10091C>G]` | refused | — | a mixed prefix *inside* one compact bracket — a genuine grammar violation; rejected at parse |
| `[NM_004006.3:c.76A>G];[NC_000023.10:g.33038255C>A]` | conformant | — | the fully-expanded cross-reference trans form — different reference sequence types across brackets. Ferro accepts it as spec-discouraged (rule-2 preference, not an absolute ban); parse-only. Deliberate, disclosed policy, no issue filed |

## `alleles.md:86-90` — Discussion: the retracted `+` separator

> !!! note "Was originally the recommendation to use the format <code class="invalid">[c.76A>C+c.83G>C]</code>?"

Ferro: the `+` separator inside allele brackets is `class="invalid"` — a retracted proposal. The
`;` is the recommended cis separator; ferro rejects `+` at parse.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G+123C>T]` | refused | — | the retracted `+` separator — not valid HGVS; rejected at parse. The recommended cis form uses `;` (`c.[76A>G;123C>T]`) |

## `alleles.md:92-97` — Discussion: recording variant combinations

> !!! note "In recessive diseases, is it important I show in which combination variants were found?"

Ferro: advisory guidance to record the cis/trans combination clearly — picks no canonical form,
constrains no normalization. Descriptive; the machinery is covered at `alleles.md:16`, `:17` and
`:20-22`.

## `alleles.md:99-101` — Discussion: the empty allele `[]`

> !!! note "I find the notation <code class="invalid">c.[76A>C]</code> without describing the second allele misleading; not enough researchers know this refers to only one of the two alleles present. Would using <code class="invalid">c.[76A>C];[]</code> be OK?"

Ferro: the empty second allele `[]` is `class="invalid"` — refused at parse. The spec's recommended
answer names the second allele's state instead, `c.[76A>C];[76=]` — the fully-bracketed
position-wild-type form ferro emits and accepts.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G];[]` | refused | — | the empty allele `[]` (`class="invalid"`) — rejected at parse. Recommended: name the second allele's state, e.g. `c.[76A>G];[76=]` |
| `NM_004006.3:c.[76A>G];[76=]` | recommended | self | the fully-bracketed form ferro emits for "variant on one allele, wild-type at that position on the other" |

## `alleles.md:103-106` — Discussion: X-chromosome male, `c.0` for an absent allele

> !!! note "How should I describe the variants detected in males and females for a gene on the X-chromosome?"

Ferro: `c.0` denotes an absent second allele (no second X-chromosome in a male). Ferro accepts the
`[0]` member as a distinct construct from a bare `c.0` coordinate, which is separately and
correctly rejected as an invalid position. The female case is the ordinary `[variant];[76=]` trans
form. The male case is subject to the same re-spelling as `[?]`: the member is preserved, but the
allele is rendered in the expanded per-bracket-accession form
([#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214)).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.[76A>C];[76=]` | recommended | — | the female case — an ordinary position-wild-type trans allele (parse-only) |
| `LRG_199t1:c.[76A>C];[0]` | recommended | — | the male case — `[0]` = absent second allele, preserved (parse-only) |
| `NM_004006.3:c.[76A>G];[0]` | conformant | `[NM_004006.3:c.76A>G];[0]` | executable twin — the `[0]` absent-allele member is accepted and kept unmerged, but the allele is re-spelled from the compact form the spec writes to the expanded one. Valid, re-parses; not the recommended shape ([#2214](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2214)) |

## `alleles.md:108-113` — Discussion: gene symbols in trans

> !!! note "I have a patient with hearing loss and variants in the _GJB2_ (`c.35delG`) and _GJB6_ (`c.689_690insT`) genes _in trans_, how should I describe this? (Nancy Carson, Ottawa, Canada)"

Ferro: the recommended standard-accession trans form is an ordinary expanded cross-reference trans
allele — accepted (parse-only). The alternative gene-symbol-only form (HGNC symbols as the sole
reference identifier) is a soft gap: ferro covers gene-selector-on-transcript forms
(`NM_X(GENE):c…`), but the gene-symbol-as-sole-identifier shape at this exact spelling is
unconfirmed.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `[NM_004004.2:c.35del];[NM_006783.1:c.689_690insT]` | recommended | — | the standard-accession trans form; foreign accessions — parse-only here |

## `alleles.md:115-125` — Discussion: haplotype short format

> !!! note "When for a haplotype I have a range of changes to report, is there a suggested short format to use?"

Ferro: the short haplotype notation is out of scope by design — it is interpretable only once the
reference sequences and variant list are "clearly specified" elsewhere, side-channel context a
context-free HGVS validator has no place to carry. Not a parser gap; intentionally out of scope.

## `alleles.md:127-140` — Discussion: Peter Taschner's suggested shorthand

> !!! note "Suggestion (Peter Taschner, Leiden)"

Ferro: the informal genotype shorthand is explicitly labelled a **suggestion**, not adopted into
the recommendations. Descriptive, no implementation expectation.

See also → `substitution.md:5` (the `NM_004006.3` base facts the executable substitution members
rely on), `RNA/alleles.md:36-52` (the `c.`-axis `[?]`/`[0]` re-spelling documented on the RNA twin),
`inversion.md:5` (where the substantive inversion ruling that only cites `alleles.md:5` as
background actually rests).
