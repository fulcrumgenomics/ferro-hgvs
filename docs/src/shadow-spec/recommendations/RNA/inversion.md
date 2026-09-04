# Inversion — ferro's reading

ferro's reading of `RNA/inversion.md`. The rules are HGVS's; ferro's job is to produce the form the
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

Two ledger records reach this page's clauses, and neither is an `r.`-jurisdiction ruling on its own
terms. `whole-span-reverse-complement-types-as-inv` governs the whole-span reverse-complement →
`inv` typing at `:5`, but it is grounded entirely on DNA clauses (`DNA/inversion.md:5`,
`general.md:55`) and guarded only on the coding and genomic axes — it carries **no `r.` test**.
`RNA/inversion.md:5` states the identical whole-span definition on the RNA axis's own authority, so
the record *could* be grounded on `r.` without stretching a DNA clause; adding that authority and an
`r.` guard is a **proposed** ledger change, not one made, and is surfaced here for the typing rule
rather than claimed as settled `r.` conformance. `inverted-duplication-is-derived-as-ins-range-inv`
cites `RNA/insertion.md:18` and `RNA/duplication.md:21` directly, so it carries `r.` authority on
its own molecule-native basis and governs `:19`. The `:15`, `:16-17`, `:18`, `:20`, `:21-22` and the
Example/Discussion clauses carry **no Why block** — the reading is CONFIRM-by-inspection against the
spec text and the shipped code, not an adjudicated ruling.

**One derivation gap sits under two of this page's clauses.** Ferro runs full inversion and
inverted-duplication derivation **only via `normalize_genome`**; the `c.`/`n.`/`r.`/`m.` axes still
emit literals. So on the `r.` axis ferro **preserves** an already-spelled `r.NNN_NNNinv` (or
`ins<range>inv`) input, but it does not **mint** an `inv` from a reverse-complement `delins`, nor
keep an `ins<range>inv` range payload — it expands the range to literal bases. This is an
implementation gap, not a spec question, and it does not change any verdict below: an already-`inv`
input is a fixed point, and a literal-expanded `ins<range>inv` output is valid HGVS that re-parses.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
Notes and Examples are written on bare, illustrative sequences (`r.177_180inv`, `r.203_506inv`,
`r.234inv`) whose bases are not the slice's, and the `:44-46` Q&A cites the foreign accession
`AB053210.2`, which the slice does not carry — so those spec spellings are shown on the slice
accession or left parse-only (`—`). The `NM_004006.3` base facts the executable rows rely on, all
measured on the slice: `r.177_180` is `gcaa`, whose reverse complement `uugc` changes all four
columns, and no other four-base inversion placement in the flank (`r.170_187` is
`ugacagggcaaaaacugc`) reproduces the same sequence; `r.203_506` is a 304-nucleotide span of which
222 columns change under inversion; and the reverse complement of `r.123_234` (112 bases) is
exactly the literal payload ferro emits for `r.234_235ins123_234inv` — the same expansion measured
on `RNA/insertion.md` and `RNA/duplication.md`. Both well-formed inversions (`r.177_180inv`,
`r.203_506inv`) are measured fixed points on the slice, not merely code-derived ones.

## `inversion.md:5` — definition: more than one nucleotide, reverse complement

> Inversion: a sequence change where, compared to a reference sequence, **more than one nucleotide** replacing the original sequence is the reverse complement of the original sequence.

Ferro: an inversion replaces a span of **more than one** nucleotide with its exact reverse
complement. When a span's whole content is its reverse complement, that span types as one `inv`,
whatever the competing split partition's members would type as — the property is a fact about the
whole span, with no term for its interior columns.

**Why.**
<!-- why:START -->
> **[whole-span-reverse-complement-types-as-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A span whose whole content is replaced by its exact reverse complement is written as one inv, however much of its interior coincides with the reference and whatever the competing partition is made of — a project choice among conformant forms, not a conformance requirement, since the spec's type-ranking rule cannot settle a merge-versus-split question at all.
<!-- why:END:whole-span-reverse-complement-types-as-inv -->

`whole-span-reverse-complement-types-as-inv` is grounded on `DNA/inversion.md:5` and `general.md:55`
and guarded on the coding and genomic axes only — no `r.` guard exists. `RNA/inversion.md:5` here
states the identical whole-span definition on the RNA axis's own authority, so grounding the rule on
`r.` needs no DNA-clause stretch; the record simply does not yet cite it, and the extension is a
worksheet proposal, not an implemented one. And per the derivation gap above, the `inv`-from-`delins`
typing runs only in `normalize_genome` — so on `r.` ferro preserves an already-spelled inversion but
does not mint one from a reverse-complement span.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | a well-formed inversion of a four-nucleotide span (>1 nt): on the slice `r.177_180` is `gcaa`, replaced by its reverse complement `uugc`, every column changing. Already an `inv`, so a measured fixed point. The executable twin of the spec's illustrative `r.177_180inv` (`:26-27`) — the spec's reverse-complement math is worked on its own bases (`cuga`), not the slice's |

## `inversion.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein level may be given in addition.

Ferro: this is guidance to authors about which level(s) to report, not a constraint on how a given
`r.` description, once written, is normalized. An `r.` inversion stands on its own — ferro does not
require an accompanying `c.`/`g.` form. This is the twin of `RNA/duplication.md:15`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | a lone RNA-level inversion is fully valid; ferro normalizes it without demanding a DNA-level companion |

## `inversion.md:16-17` — more than one nucleotide; a one-nucleotide inversion is a substitution

> - by definition, the region inverted (`positions_inverted`) contains **more than one nucleotide**.
>   The description <code class="invalid">r.234inv</code> is therefore not allowed; a one-nucleotide inversion should be described as a [substitution](substitution.md)

Ferro: two points. (i) An inversion spans more than one nucleotide, so the grammar admits only a
*range* — a lone `r.234inv` is refused at parse. (ii) A one-nucleotide inversion should be written
as a substitution, which is definitional: the reverse complement of a single base is its complement,
i.e. a `sub`, so there is no lone-`inv` spelling for ferro to relabel.

**The `r.`-axis one-nucleotide redirect is un-adjudicated.** No ledger record scopes the "one
nucleotide → substitution" redirect on the RNA axis; the DNA analogue at `DNA/inversion.md:16` is cited
only inside the DNA whole-span record, and `whole-span-reverse-complement-types-as-inv`'s own "more
than one nucleotide" floor is DNA-grounded (see `:5`). The behaviour is nonetheless CONFIRM by
inspection: the grammar forecloses the lone-`inv` spelling before any redirect is needed, and ferro
mints no `inv` from a one-nucleotide change that would require the redirect. Treat the `r.`-axis
redirect as **open/provisional** — conformant by construction, but not recorded on its own authority.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | a valid inversion of more than one nucleotide — the shape `:16` requires |
| `NM_004006.3:r.234inv` | refused | — | the spec's own `class="invalid"` single-position `inv` — rejected at parse: an inversion requires a range of more than one nucleotide, so a lone position carries no span to invert |

## `inversion.md:18` — the 3'rule

> for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).

Ferro: the general 3'-shift rule applies to inversions exactly as to the other edit types — an
inversion that could be written at more than one equivalent position is assigned its most 3'
placement. A shiftable inversion is hard to exhibit on this transcript: an inversion 3'-shifts only
when it sits inside a larger reverse-complement-symmetric context, which the slice's homopolymer runs
do not supply, so the executable row below is a measured fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | no equivalent 3' placement exists for this span — measured: over the flank `r.170_187` (`ugacagggcaaaaacugc`), `r.177_180` is the only four-base inversion placement that yields the same sequence — so the 3'rule has nothing to shift; already 3'-most |

## `inversion.md:19` — an inverted duplication is `ins<range>inv`, not `dupinv`

> **inverted duplications** are described as an insertion (`r.234_235ins123_234inv`), not as <code class="invalid">r.123_456dupinv</code>.

Ferro: an inverted duplication is spelled `ins<range>inv`, naming the span the inverted copy came
from — never `dup`, and never the `dupinv` shorthand, which is refused at parse. This clause is the
`r.`-native basis for the `ins<range>inv` form worked from the insertion and duplication sides.

**Why.**
<!-- why:START -->
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum, so a short chance reverse-complement match is not misread as one.
<!-- why:END:inverted-duplication-is-derived-as-ins-range-inv -->

`inverted-duplication-is-derived-as-ins-range-inv` reaches the `r.` axis on RNA's own authority
(`RNA/insertion.md:18`, `RNA/duplication.md:21`), not by stretching a `DNA/` clause. Per the
derivation gap above, the mint/keep of the `ins<range>inv` spelling is wired only in
`normalize_genome`, so on `r.` ferro **expands** an already-spelled range payload to literal
reverse-complement bases via `try_expand_rna_ins` rather than preserving the range-inv form —
valid HGVS that re-parses and denotes the same sequence, so `conformant`, not the recommended
spelling. #1946's render stage is the named fix.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.234_235ins123_234inv` | conformant | `NM_004006.3:r.234_235insguugacauuguucagggcaugaacucuuguggauccuuuuucuuuuggcaguuuuugcccugucaggccuucgaggaggucuaggaggcgccucccauccuguaggucacug` | the spec's own inverted-duplication spelling — an `ins<range>inv`, on the slice accession. Ferro expands the range payload to literal reverse-complement bases (`try_expand_rna_ins`) because the range-inv derivation is wired only in `normalize_genome`; the same output and gap as `RNA/insertion.md:18` and `RNA/duplication.md:18-21` (`inverted-duplication-is-derived-as-ins-range-inv`; #1946 names the render-stage fix). Valid HGVS that re-parses and denotes the same sequence, so conformant, not the recommended spelling |
| `NM_004006.3:r.123_456dupinv` | refused | — | the `dupinv` shorthand `:19` marks `class="invalid"` — rejected at parse ("Unexpected trailing characters: 'inv'"); the clause requires the `ins<range>inv` form. Pinned by `reject_dupins`/`recommended_form_pins` |

## `inversion.md:20` — large genomic inversions on RNA give deletion or delins

> since exon splice signals will be inverted, large genomic inversions on the RNA level usually give [deletion](deletion.md) or [deletion-insertion (delins)](delins.md) variants

Ferro: descriptive biology, not a normalization operation. It observes that a large genomic
inversion, seen on the RNA level, typically surfaces as a deletion or delins because the inverted
region disrupts splice signals — a statement about what the observed RNA looks like, not a rewrite
ferro applies to a given `r.` inversion. Ferro emits no automatic inversion → deletion/delins
recasting, so there is no verdict row.

## `inversion.md:21-22` — inversions are not used on the protein level

> - inversions are not used on protein level.
>   Depending on the (predicted) consequences of an inversion on protein level, changes are usually described as either a **delins** or a **frameshift**.

Ferro: descriptive protein-level redirect. The `p.` axis has no `inv` edit type at all; a
protein-level consequence of an inversion is expressed as a delins or a frameshift instead. Ferro
emits no protein `inv`, so this is CONFIRM-by-inspection with no verdict row.

## `inversion.md:26-27` — Example: a four-nucleotide inversion

> - **`r.177_180inv`**<br>
>   inversion of nucleotides `r.177` to `r.180`, changing `..agg`<code class="del">cuga</code>`uu..` to `..agg`<code class="ins">ucag</code>`uu..`.

Ferro: the spec's own worked example, and its reverse-complement math checks out — the reverse of
`cuga` is `aguc`, whose complement is `ucag` (`u`↔`a`, `c`↔`g`), exactly the `<code class="ins">`
result. A four-nucleotide whole-span reverse complement types as one `r.177_180inv`
(`whole-span-reverse-complement-types-as-inv`, cited at `:5`, with the DNA-only-guard caveat noted
there). The spec spelling is accession-less and worked on illustrative bases; the executable twin on
the slice accession is a well-formed inversion and a fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | executable twin of the spec's bare `r.177_180inv` — a four-nucleotide inversion (`gcaa` → `uugc` on the slice, against the spec's illustrative `cuga` → `ucag`), already `inv` and 3'-most, so a measured fixed point (the accession-less spec spelling is not table-executable) |

## `inversion.md:29-30` — Example: a 304-nucleotide inversion

> - **`r.203_506inv`**<br>
>   inversion of the 304 nucleotides from position `r.203` to `r.506`.

Ferro: the spec's larger worked example — a 304-nucleotide range inversion. No bases are given, so
there is nothing to recompute; it is a well-formed `inv` over a span far larger than one nucleotide,
which satisfies `:16`'s lower bound. Shown on the slice accession as a measured fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.203_506inv` | recommended | self | executable twin of the spec's bare `r.203_506inv` — a well-formed 304-nucleotide range inversion (222 of its 304 columns change on the slice), a measured fixed point (the accession-less spec spelling is not table-executable) |

## `inversion.md:34-37` — Discussion: the complement is not an inversion

> !!! note "Is the change `aagc` to `uucg` an inversion?"
>
>     No, an inversion would change `aagc` to `gcuu`, its **reverse-complement**.
>     `uucg` is only the **complement** of `aagc`.

Ferro: a clarification, not a rule ferro applies. It reinforces the `:5` definition — an inversion
is the **reverse** complement, so `aagc` inverts to `gcuu`, whereas `uucg` is only the complement
(bases flipped, order kept). Ferro's inversion handling already uses the reverse complement, so this
is CONFIRM-by-inspection with no verdict row.

## `inversion.md:39-42` — Discussion: the reverse is not an inversion

> !!! note "Is the change `aagc` to `cgaa` an inversion?"
>
>     No, an inversion would change `aagc` to `gcuu`, its **reverse-complement**.
>     `cgaa` is only the **reverse** of `aagc`.

Ferro: the mirror-image clarification of the one above — `cgaa` is only the reverse (order flipped,
bases kept), not the reverse **complement**. Both Q&A blocks fence the definition against its two
near-misses (complement-only and reverse-only). Descriptive; no verdict row.

## `inversion.md:44-46` — Discussion: the old "o" strand marker is not used

> !!! note "On the [old nomenclature website](http://www.HGVS.org/mutnomen/examplesRNA.html) (bottom), you had the example <code class="invalid">r.124_500delinsoAB053210.2:r.1289-365_1289-73</code>, i.e. the "o" indicating the inserted sequence `AB053210.2:r.1289-365_1289-73` was from the opposite transcriptional strand. Is the "o" still used?"
>
>     No, the "o" is not used, the insertion is in an inverted orientation, so "inv" should be used; <code class="invalid">r.124_500delins[AB053210.2:r.1289-365_1289-73inv]</code>.

Ferro: descriptive history. The deprecated `o` strand marker is gone; an inserted sequence taken
from the opposite strand is spelled by an `inv`-suffixed insertion payload
(`delins[…inv]`) instead. This is about how an **insertion payload** carries orientation, not about
inversion typing, so it touches the insertion recommendations rather than this page's normalization
behaviour. The cited spelling is on the foreign accession `AB053210.2` and is not in the slice —
there is no executable verdict row.

See also → `RNA/insertion.md:18` and `RNA/duplication.md:18-21` (the `ins<range>inv`
inverted-duplication form, worked from the insertion and duplication sides), `RNA/substitution.md`
(where a one-nucleotide "inversion" lands, per `:16-17`), `DNA/inversion.md` (where the whole-span
reverse-complement → `inv` typing derivation is wired and guarded today).
