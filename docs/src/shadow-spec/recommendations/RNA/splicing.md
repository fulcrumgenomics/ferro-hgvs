# Splicing — ferro's reading

ferro's reading of the HGVS **splicing** recommendations on the transcript (`r.`) axis, clause
by clause, with each spelling's normalized form and a verdict. See
[How to read a page](../../reading-guide.md) for verdicts, table conventions, and recurring terms.

**No ledger record cites any `RNA/splicing.md` clause** — the whole document is a record-level gap,
so **every section below carries no Why block**; the reading is CONFIRM-by-inspection against the
spec text and the shipped code. Three findings on this page are live and unrecorded, written up as
Notes rather than rulings:

- the intron-retention `r.` insert (`splicing.md:34`, `:45`, `:50-57`) — the page's central form —
  is spec-**recommended** (five bold examples), yet ferro declines it at normalize with a
  *spec-undefined* reason that is a **misstatement on the `r.` axis**. The proposed reading is that
  the form is spec-defined and the house choice should be coordinate-space-resolvability: flatten a
  range the rendering reference can resolve, **preserve** an intronic-offset range it cannot.
  Documented at `splicing.md:27-36`;
- the same shape is wrapped in `<code class="invalid">` at `RNA/insertion.md:44` — a
  formatting accident that contradicts this page's five bold examples; `splicing.md` governs, and
  the upstream inconsistency is owed a filing (noted at `splicing.md:27-36`);
- `r.(spl?)` (`splicing.md:79`) versus `r.spl?` (`uncertain.md:194`) — one statement, two spellings,
  confluence unadjudicated; and `r.(spl)`, which appears nowhere in the spec yet ferro's
  **projector** manufactures it. Both open; documented at `splicing.md:75-80`.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
splicing examples sit on accessions outside the slice, so those rows are parse-only (`—`), and the
intron-retention forms among them additionally decline at normalize (see `splicing.md:27-36`). The
four splice markers (`r.spl`, `r.spl?`, `r.(spl?)`, `r.(spl)`) are reserved tokens reachable only on
the `r.` parse path — `c.spl`/`g.spl` do not parse.

## `splicing.md:5` — definition

> Splicing: a sequence change where, compared to a reference sequence, the normal RNA splicing pattern is altered.

Ferro: this defines a **biological outcome**, not a description form. Splicing has no edit type of
its own (see the Syntax below), so there is nothing here a normalizer applies or violates.
**Descriptive.**

## `splicing.md:8-11` — syntax

> Variants affecting RNA splicing result in either a [deletion](deletion.md) or [insertion](insertion.md) on the RNA level and should be described as such.
> Recommendations on representing adjoined transcripts formed by gene fusions are discussed in the Notes and Examples below.

Ferro: a splicing event is written as a `del` or an `ins` on the `r.` axis — it has no edit type of
its own. A normalizer receives those `del`/`ins` descriptions and treats them as any other RNA deletion
or insertion (governed by `RNA/deletion.md` and `RNA/insertion.md`, not here). Ferro does not derive
splicing outcomes from a DNA variant, so there is nothing to measure at this clause. **Descriptive.**

## `splicing.md:14` — describe on the DNA level

> - all variants **should be** described on the DNA level; descriptions on the RNA and/or protein level may be given in addition.

Ferro: authoring guidance about which molecular level(s) a submitter should report on — a lowercase
"should" over a *level* choice, not a constraint on how a given `r.` description, once written, is
normalized. A lone `r.` description stands on its own; ferro does not require an accompanying `c.`/`g.`
form. Nothing for the normalizer to enforce. **Descriptive.**

## `splicing.md:15` — the comma / Products separator

> - a `,` (comma) is used to separate different transcripts/proteins derived from one allele; `r.[123a>u,122_154del]`.

Ferro: the `,` separates **different transcripts** — downstream products of one DNA event — a
fourth allele relationship distinct from cis `;`, trans `];[`, and unknown `(;)`. Ferro models it
as `AllelePhase::Products`, restricts the form to the `r.`/`p.` axes, normalizes each member in
its own frame, and — because Products is not `Cis` — **never merges or reorders** the members.
Members may **overlap** in coordinates (`123a>u` sits inside `122_154del`) without that being a
conflict, and are **not position-ordered** (`r.[897u>g,832_960del]` lists 897 before 832). The
full treatment lives on `RNA/alleles.md:21`, whose `general.md:80` cross-reference gives the bare
form `r.[123a>u,122_154del]`.

**No ledger record adjudicates Products-allele normalization semantics** — so this section carries
no Why block. The behaviour is a defensible but unrecorded ferro reading.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.[6_8del,124_129del]` | recommended | self | executable Products allele — two distinct-transcript products, each a fixed point; `AllelePhase::Products`, kept unmerged and unreordered (matches `comma_products_allele.rs`). Illustrative positions; `general.md:80`'s own string is the bare `r.[123a>u,122_154del]` |
| `NM_004006.3:r.[76del,124_129del]` | recommended | `NM_004006.3:r.[77del,124_129del]` | each product normalized in its own frame — the first 3'-shifts to `r.77del`, the second is untouched, the `,` join and order are preserved |
| `NM_004006.3:c.[123del,124del]` | refused | — | the `,` form is `r.`/`p.`-only — a `c.` Products allele is rejected at parse (`find_products_bracket` restricts the comma form to the RNA/protein axes) |

See also → `RNA/alleles.md:21` (the full Products treatment), `general.md:80` (the `,` separator and
its bare example).

## `splicing.md:16-19` — HGNC / VICC gene-fusion nomenclature

<details class="ss-spec"><summary>Spec discussion (verbatim)</summary>

> - HGVS recommends following the [HGNC guidelines](https://www.genenames.org/about/guidelines/) and the [VICC Gene Fusion Specification](https://fusions.cancervariants.org/en/latest) nomenclature to describe products of gene fusions.
>     - The HGNC recommendations include using a `GENESYMBOL1::GENESYMBOL2` syntax for gene-level fusion descriptions, and a `GENESYMBOL1-GENESYMBOL2` syntax for read-through transcripts.
>     - The VICC nomenclature extends the HGNC recommendations to include a terminology, information model, and nomenclature for gene-level and exon-level representation, with components for disambiguating regulatory fusions from chimeric transcript fusions.
>     - HGVS also recommends the use of adjoined transcripts (see examples) for precise and unambiguous characterization of chimeric transcripts at the sequence level.

</details>

Ferro: gene-symbol-level fusion nomenclature is delegated to HGNC/VICC and is **not an HGVS
description**. The only HGVS-level form here is the adjoined transcript, which is governed by its own
page (`RNA/adjoined_transcript.md`) and shown at `splicing.md:59-61` below. Nothing on these lines is
a normalization rule. **Descriptive.**

## `splicing.md:23-25` — example: one variant, several transcripts

> - **one variant, several transcripts**
>     - **`NC_000023.11(NM_004006.2):r.[897u>g,832_960del]`**<br>
>       two different transcripts, `r.897u>g` and `r.832_960del`, derive from one variant (`c.897T>G` on the DNA level).

Ferro: the worked example of the `splicing.md:15` Products rule — two members from **one** DNA variant,
listed out of position order (897 before 832). Ferro keeps them as `AllelePhase::Products`, normalizes
each in its own frame, and neither merges nor reorders them. (Adjudication is unrecorded — see
`splicing.md:15`; no Why block.)

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):r.[897u>g,832_960del]` | recommended | — | the spec's own Products example; two transcripts from one DNA event (`c.897T>G`), kept distinct and unordered; parse-only here (accession outside the slice) |
| `NM_004006.3:r.[6_8del,124_129del]` | recommended | self | executable Products twin — two distinct-transcript products, each a fixed point, kept unmerged |

## `splicing.md:27-36` — example: splice acceptor site

<details class="ss-spec"><summary>Spec discussion (verbatim)</summary>

> - **splice acceptor site**
>     - **`NC_000023.11(NM_004006.2):r.650_831del`**<br>
>       as a consequence of a variant destroying a splice acceptor site, the sequence from nucleotide `r.650` to `r.831` (exon 8) is deleted from the transcript.
>
>     - **`NC_000023.11(NM_004006.2):r.650_712del`**<br>
>       as a consequence of a variant destroying a splice acceptor site, a new acceptor site in exon 8 (positions `712` / `713`) is activated and the sequence from nucleotide `r.650` to `r.712` is deleted from the transcript.
>
>     - **`NC_000023.11(NM_004006.2):r.649_650ins[650-52_650-2;c]`**<br>
>       as a consequence of a variant destroying a splice acceptor site (`c.650-1G>C`), a new acceptor site in intron 7 is activated and the intron 7 sequence from positions `650-52` to `650-1` is inserted in the transcript.<br>
>       **NOTE**: nucleotide `650-1` changed from `g` to `c`.

</details>

Ferro: two of these are ordinary RNA deletions (governed by `RNA/deletion.md`) and preserve; the
third is the page's **central form** — intron retention written as an `r.` insertion whose payload
is a **range in intronic (`c.`-style offset) coordinates**, with the mutated base spelled literally
in a `[…;…]` composite. It is bold with no `invalid` markup: the spec's recommended way to write
intron retention. Ferro parses it, then **declines at normalize** as spec-undefined, producing no
output. Nothing malformed is emitted, so the row is `conformant`, not `bug` and not `refused`
(strict parse accepts this input).

**Note — the decline reason misstates the spec on `r.`.** The spec **defines** this form five
times on this page as the recommended way to describe intron retention, so labelling the decline
"spec-undefined" is wrong, even though declining may be legitimate: expanding the range to literal
bases needs the pre-mRNA sequence the `r.` reference does not hold. The proposed reading is
**coordinate-space-resolvability**: flatten a range whose coordinates the rendering reference
resolves, and **preserve** one whose coordinates lie outside it (an intronic-offset range →
preserve, not decline). No ledger record exists for this yet. Ferro's flatten-vs-preserve policy is
separately inconsistent across payload shapes on this axis — a resolvable same-reference `r.` range
payload flattens to a literal elsewhere, with no ledger record for that choice either — so the
intronic-range shape inherits a decline instead of either half of an unrecorded policy.

**Note — an upstream inconsistency, owed a filing.** `RNA/insertion.md:44` wraps the identical
shape in `<code class="invalid">` while describing it as a valid intron-retention insertion and
cross-linking back to this page — a formatting accident that contradicts this page's five bold
examples. Until answered, **`splicing.md` (five bold examples) governs over one formatting-commit
class attribute.**

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):r.650_831del` | recommended | — | exon-8 deletion — an ordinary RNA deletion (governed by `RNA/deletion.md`), preserved; parse-only here |
| `NC_000023.11(NM_004006.2):r.650_712del` | recommended | — | new-acceptor-site deletion — preserved; parse-only here |
| `NC_000023.11(NM_004006.2):r.649_650ins[650-52_650-2;c]` | conformant | — | the spec's central intron-retention form — a range in intronic offset coordinates plus the mutated base `c`. Ferro parses it (`InsertedPart::CdsPositionRange`) then **declines at normalize** with the *spec-undefined* message. Nothing wrong is emitted, so `conformant`; parse-only here. The message misstates the spec on `r.` and the decided preserve policy is unrecorded — see the Notes above |

## `splicing.md:38-47` — example: splice donor site

<details class="ss-spec"><summary>Spec discussion (verbatim)</summary>

> - **splice donor site** (`c.831+2T>A`)
>     - **`NC_000023.11(NM_004006.2):r.650_831del`**<br>
>       as a consequence of a variant destroying the exon 8 splice donor site, the sequence from nucleotide `r.650` to `r.831` (exon 8) is deleted from the transcript.
>
>     - **`NC_000023.11(NM_004006.2):r.778_831del`**<br>
>       as a consequence of a variant destroying the exon 8 splice donor site, a new donor site in exon 8 (positions `777` / `778`) is activated and the sequence from nucleotide `r.778` to `r.831` is deleted from the transcript.
>
>     - **`NC_000023.11(NM_004006.2):r.831_832ins[ga;831+3_831+60]`**<br>
>       as a consequence of a variant destroying the exon 8 splice donor site, a new donor site in intron 8 (positions `831+60` / `831+61`) is activated and the intron 8 sequence from positions `831+1` to `831+60` is inserted in the transcript.<br>
>       **NOTE**: nucleotide `831+2` changed from `u` to `a`.

</details>

Ferro: the acceptor-site shape at `splicing.md:27-36`, mirrored to the donor side — the literal `ga`
spells `831+1` and the mutated `831+2`, and the range resumes at `831+3`. Two ordinary deletions
preserve; the intron-retention insert declines at normalize with the same *spec-undefined* message.
Same finding as `splicing.md:27-36`; the false wording and the decided preserve policy are documented
there.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):r.650_831del` | recommended | — | exon-8 deletion (governed by `RNA/deletion.md`), preserved; parse-only here |
| `NC_000023.11(NM_004006.2):r.778_831del` | recommended | — | new-donor-site deletion — preserved; parse-only here |
| `NC_000023.11(NM_004006.2):r.831_832ins[ga;831+3_831+60]` | conformant | — | the donor-side intron-retention form — a literal `ga` plus an intronic-offset range. Parses, then declines at normalize with the same *spec-undefined* message; nothing wrong emitted, so `conformant`; parse-only here. See `splicing.md:27-36` |

## `splicing.md:49-57` — example: intron variant

<details class="ss-spec"><summary>Spec discussion (verbatim)</summary>

> - **intron variant**
>     - **`NC_000023.11(NM_004006.2):r.649_650ins650-50_650-1`**<br>
>       as a consequence of an intron 7 variant (`c.650-52_650-51del`), a new stronger exon 8 splice acceptor site is created (positions `650-51` / `650-50`) and the intron 7 sequence from positions `650-50` to `650-1` is inserted in the transcript.
>
>     - **`NC_000023.11(NM_004006.2):r.831_832ins831+1_831+67`**<br>
>       as a consequence of an intron 8 variant (`c.831+71C>A`), a new stronger exon 8 splice donor site is created (positions `831+67` / `831+68`) and the intron 8 sequence from positions `831+1` to `831+67` is inserted in the transcript.
>
>     - **`NC_000023.11(NM_004006.2):r.649_650ins650-1400_650-1268`**<br>
>       as a consequence of an intron 7 variant (`c.650-1401T>G`), a new exon is created and its sequence (positions `650-1400` to `650-1268`) is inserted in the transcript.

</details>

Ferro: the **unbracketed** intronic-range payload — the plainest form of the intron-retention shape,
with no literal component at all. All three decline at normalize with the *spec-undefined* message.
This is the strongest evidence that "spec-undefined" is the wrong word: there is no composite, no
literal, nothing but the spec's own canonical spelling of intron retention. Same finding as
`splicing.md:27-36`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):r.649_650ins650-50_650-1` | conformant | — | a new-acceptor-site intron-retention insert, unbracketed intronic range. Parses, declines at normalize (*spec-undefined*); nothing wrong emitted; parse-only here |
| `NC_000023.11(NM_004006.2):r.831_832ins831+1_831+67` | conformant | — | the donor-side unbracketed intronic range — same shape, same decline |
| `NC_000023.11(NM_004006.2):r.649_650ins650-1400_650-1268` | conformant | — | a cryptic-exon (pseudo-exon) inclusion — same shape, same decline |

## `splicing.md:59-61` — example: adjoined transcript

> - **adjoined transcript** (based on [SVD-WG007](../../consultation/SVD-WG007.md))
>     - **`NM_002354.2:r.-358_555::NM_000251.2:r.212_*279`**<br>
>       describes an adjoined transcript from an `EPCAM::MSH2` gene fusion, where nucleotides `r.-358` to `r.555` (_EPCAM_ gene, reference transcript `NM_002354.2`) are spliced to nucleotides `r.212` to `r.*279` (_MSH2_ gene, reference transcript `NM_000251.2`).

Ferro: a cross-reference to `RNA/adjoined_transcript.md` — nothing on this page adds a rule. Ferro
parses the `::`-joined adjoined form and preserves it. The reading is owned by the adjoined-transcript
page, not here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279` | recommended | — | the spec's own `EPCAM::MSH2` adjoined transcript; preserved. Governed by `RNA/adjoined_transcript.md`; both accessions are outside the slice, so parse-only here |

## `splicing.md:63-71` — example: uncertain (RNA not analysed)

<details class="ss-spec"><summary>Spec discussion (verbatim)</summary>

> - **uncertain** (RNA not analysed)
>     - **`NC_000023.11(NM_004006.2):r.(76a>c)`**<br>
>       RNA was not analysed, but a substitution of the `a` nucleotide at `r.76` by a `c` is predicted.
>
>     - **`NC_000023.11(NM_004006.2):r.?`**<br>
>       an effect on the RNA level is expected, but it is not possible to give a reliable prediction of the consequences (RNA not analysed).
>
>     - **`NC_000023.11(NM_004006.2):r.spl`**<br>
>       RNA has not been analysed, but it is very likely that splicing is affected.

</details>

Ferro: three whole-transcript or predicted statements. `r.(76a>c)` is an ordinary predicted
substitution (the parentheses mark it predicted); `r.?` and `r.spl` carry no position and no sequence
change, so normalization can only be identity. All three preserve. `spl` is reserved and reachable only
on the `r.` parse path (`parser/variant.rs`; `c.spl`/`g.spl` do not parse), and `r.spl` is kept
distinct from `r.?` (`rna_spl_marker.rs`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):r.(76a>c)` | recommended | — | the spec's predicted substitution; parse-only here |
| `NM_004006.3:r.(76a>c)` | recommended | self | executable twin — `r.76` reads `a` (from the `aaa` run at `r.75_77`), matching the spec's `a>c`; the predicted substitution is preserved, the parentheses kept |
| `NC_000023.11(NM_004006.2):r.?` | recommended | — | whole-transcript unknown RNA effect; parse-only here |
| `NM_004006.3:r.?` | recommended | self | executable twin — the whole-transcript `?` is an identity, preserved |
| `NC_000023.11(NM_004006.2):r.spl` | recommended | — | the `r.spl` splice marker; parse-only here |
| `NM_004006.3:r.spl` | recommended | self | executable — `spl` is reachable only on the `r.` parse path; preserved as itself, kept distinct from `r.?` (`rna_spl_marker.rs`) |

## `splicing.md:75-80` — discussion: `r.spl` and `r.(spl?)`

<details class="ss-spec"><summary>Spec discussion (verbatim)</summary>

> !!! note "A variant changes the +1 intron sequence (`GT` to `AT`). Although I did not analyse RNA, I am quite sure that normal splicing is affected. How can I best indicate this?"
>
>     HGVS recommends to use the format `r.spl` to indicate that RNA was not analysed, but splicing is most probably affected.
>     In general, the format is used for variants changing the +1, +2, -2 and -1 position of an intron, i.e. affecting the `GT` splice donor and `AG` splice acceptor site (excluding `GT` to `GC` and `GC` to `GT` variants).
>     `r.(spl?)` is frequently used to indicate that normal splicing might be affected as a consequence of variants in the first or last nucleotide of an exon, the +3 to +5 intron position (splice donor site), and variants generating a new `AG` di-nucleotide close to the normal splice acceptor site (`AG`).
>     See [Uncertain](../uncertain.md).

</details>

Ferro: three things live in this Q&A, two of them unrecorded open questions.

- *When to write `r.spl` versus the weaker marker* (±1/±2 versus exon-edge/+3..+5/new AG, `GT`↔`GC`
  excluded) is guidance to an author or effect predictor. Ferro's effect module never constructs
  the splice edit outside the parser, so it neither applies nor violates this. Descriptive.

- *The weaker marker's spelling differs across two spec pages*: this page writes `r.(spl?)`;
  `uncertain.md:194` writes `r.spl?` for the identical statement. The spec never says whether they
  are equivalent, and since `(…)` means "predicted / RNA not analysed" while `r.spl?` already
  asserts RNA-not-analysed, the parens arguably add nothing. Ferro keeps `r.spl?` and `r.(spl?)` as
  **distinct, non-canonicalized** values — a defensible but ungoverned choice; no ledger record
  exists.

- *`r.(spl)` appears nowhere in the spec.* Ferro accepts it as the symmetric predicted-but-certain
  form and, on the **projector** (not the normalizer), manufactures it: projecting `r.spl` onto the
  rna axis returns the uncertain-wrapped splice edit, rendered `r.(spl)`, on a stability argument.
  Emitting a spelling the spec never shows, for an input the spec does show, is open on the
  projection axis; `adjudication-precedence-order` ranks stability as a tiebreaker only. The
  provisional preference is for the projector to emit the input's own spelling instead.

**Note.** The projector also once carried the RNA-only `spl` token verbatim onto other axes
(`g.spl`, `c.1spl`, `n.245spl`); it now refuses those axes because an RNA-level effect has no
genomic representation. That principle is stated only in code, with no ledger record. The
normalizer is unaffected — on the normalize path all four markers are identity fixed points.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.spl` | recommended | self | the strong marker (±1/±2 intron positions); preserved as itself |
| `NM_004006.3:r.(spl?)` | recommended | self | `splicing.md:79`'s weaker marker; preserved |
| `NM_004006.3:r.spl?` | recommended | self | the same weaker marker as `uncertain.md:194` spells it; preserved, and kept **distinct** from `r.(spl?)` (`rna_spl_marker.rs`). Whether the two should converge is unadjudicated — see the Note above |
| `NM_004006.3:r.(spl)` | conformant | self | a spelling the spec never publishes; ferro accepts and preserves it on the normalize path, and its **projector** manufactures it from `r.spl` on a stability argument. Valid HGVS, re-parses; open on the watched projection axis (provisional house choice: emit the input's own spelling), unrecorded — see the Note above |

## `splicing.md:82-84` — discussion: the protein level

> !!! note "How can I best describe the predicted consequences on the protein level of a variant that most probably affects splicing?"
>
>     The best format seems to use `p.?`, meaning "I do not know what to expect on the protein level".

Ferro: hedged ("seems"), about the protein axis, authored-description advice — not a normalization
rule. One observation for the projection axis: the projector **declines** the protein axis for
`r.spl` ("has no representation on another axis") rather than emitting `p.?`, while the spec's soft
preference here would be `p.?`. That is weaker than rule-2 force and belongs to any future
projection-axis record, not a ledger change of its own. Nothing for the normalizer to enforce.
**Descriptive.**
