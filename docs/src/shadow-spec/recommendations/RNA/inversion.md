# Inversion — ferro's reading

ferro's reading of the HGVS **inversion** recommendations on the transcript (`r.`) axis, clause
by clause — each spelling with the form ferro normalizes it to and a verdict on that output.
New here? See [How to read a page](../../reading-guide.md) for the verdicts, the table
conventions, and the recurring terms.

*DNA twin: [Inversion (`c.`/`g.`)](../DNA/inversion.md).*

<details class="ss-notes"><summary>Page notes — ledger grounding, the <code>r.</code>-axis derivation gap, and the slice base facts</summary>

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

</details>

## `inversion.md:5` — definition: a reverse-complemented span of more than one base

> Inversion: a sequence change where, compared to a reference sequence, **more than one nucleotide** replacing the original sequence is the reverse complement of the original sequence.

Ferro: an inversion replaces a span of **more than one** nucleotide with its exact reverse
complement. A whole span that equals its own reverse complement types as a single `inv`, whatever a
competing split partition would call its interior.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[whole-span-reverse-complement-types-as-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A span whose whole content is replaced by its exact reverse complement is written as one inv, however much of its interior coincides with the reference and whatever the competing partition is made of; this is a project choice among conformant forms, not a conformance requirement.
<!-- why:END:whole-span-reverse-complement-types-as-inv -->

</details>

<details class="ss-notes"><summary>On the <code>r.</code> axis: grounding &amp; the derivation gap</summary>

`whole-span-reverse-complement-types-as-inv` is grounded on `DNA/inversion.md:5` and `general.md:55`
and guarded on the coding and genomic axes only — no `r.` guard exists. `RNA/inversion.md:5` here
states the identical whole-span definition on the RNA axis's own authority, so grounding the rule on
`r.` needs no DNA-clause stretch; the record simply does not yet cite it, and the extension is a
worksheet proposal, not an implemented one. And per the derivation gap above, the `inv`-from-`delins`
typing runs only in `normalize_genome` — so on `r.` ferro preserves an already-spelled inversion but
does not mint one from a reverse-complement span.

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | well-formed 4-nt inversion; on the slice `r.177_180` = `gcaa` → `uugc`. Already `inv`, so a fixed point (the spec's own math uses illustrative bases `cuga`) |

## `inversion.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein level may be given in addition.

Ferro: author guidance about which level(s) to report, not a normalizer rule. An `r.` inversion
stands on its own — ferro requires no accompanying `c.`/`g.` form. Twin of `RNA/duplication.md:15`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | a lone RNA-level inversion is valid; no DNA companion required |

## `inversion.md:16-17` — more than one nucleotide; a one-nucleotide inversion is a substitution

> - by definition, the region inverted (`positions_inverted`) contains **more than one nucleotide**.
>   The description <code class="invalid">r.234inv</code> is therefore not allowed; a one-nucleotide inversion should be described as a [substitution](substitution.md)

Ferro: (i) the grammar admits only a *range*, so a lone `r.234inv` is refused at parse; (ii) a
one-nucleotide inversion is by definition a substitution (the reverse complement of one base is its
complement), so there is no lone-`inv` spelling for ferro to relabel.

<details class="ss-notes"><summary>On the <code>r.</code> axis: the one-nucleotide redirect is un-adjudicated</summary>

No ledger record scopes the "one nucleotide → substitution" redirect on the RNA axis; the DNA
analogue at `DNA/inversion.md:16` is cited only inside the DNA whole-span record, and
`whole-span-reverse-complement-types-as-inv`'s own "more than one nucleotide" floor is DNA-grounded
(see `:5`). The behaviour is nonetheless CONFIRM by inspection: the grammar forecloses the lone-`inv`
spelling before any redirect is needed, and ferro mints no `inv` from a one-nucleotide change that
would require the redirect. Treat the `r.`-axis redirect as **open/provisional** — conformant by
construction, but not recorded on its own authority.

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | a valid inversion over >1 nt — the shape `:16` requires |
| `NM_004006.3:r.234inv` | refused | — | the spec's `class="invalid"` single-position `inv` — refused at parse; an inversion needs a range |

## `inversion.md:18` — the 3'rule

> for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).

Ferro: the 3'-shift rule applies to inversions as to any edit type — an inversion writable at
several equivalent positions takes its most 3' one. A shiftable inversion is hard to exhibit here
(the slice's homopolymers lack the needed reverse-complement symmetry), so the row below is a fixed
point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | no equivalent 3' placement over the flank `r.170_187`, so the 3'rule has nothing to shift; already 3'-most |

## `inversion.md:19` — an inverted duplication is `ins<range>inv`, not `dupinv`

> **inverted duplications** are described as an insertion (`r.234_235ins123_234inv`), not as <code class="invalid">r.123_456dupinv</code>.

Ferro: an inverted duplication is spelled `ins<range>inv`, naming the span the inverted copy came
from — never `dup`, never the `dupinv` shorthand (refused at parse). This is the `r.`-native basis
for the `ins<range>inv` form.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum.
<!-- why:END:inverted-duplication-is-derived-as-ins-range-inv -->

</details>

<details class="ss-notes"><summary>On the <code>r.</code> axis: why this is <code>conformant</code>, not recommended</summary>

`inverted-duplication-is-derived-as-ins-range-inv` reaches the `r.` axis on RNA's own authority
(`RNA/insertion.md:18`, `RNA/duplication.md:21`), not by stretching a `DNA/` clause. Per the
derivation gap above, the mint/keep of the `ins<range>inv` spelling is wired only in
`normalize_genome`, so on `r.` ferro **expands** an already-spelled range payload to literal
reverse-complement bases via `try_expand_rna_ins` rather than preserving the range-inv form —
valid HGVS that re-parses and denotes the same sequence, so `conformant`, not the recommended
spelling. #1946's render stage is the named fix.

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.234_235ins123_234inv` | conformant | `NM_004006.3:r.234_235insguugacauuguucagggcaugaacucuuguggauccuuuuucuuuuggcaguuuuugcccugucaggccuucgaggaggucuaggaggcgccucccauccuguaggucacug` | the spec's inverted-duplication spelling. On `r.` ferro expands the range payload to literal reverse-complement bases (`try_expand_rna_ins`) — the derivation is wired only in `normalize_genome` — so valid HGVS, same sequence, **conformant** not recommended. Same gap as `RNA/insertion.md:18`; #1946 names the fix |
| `NM_004006.3:r.123_456dupinv` | refused | — | the `dupinv` shorthand, marked `class="invalid"` — refused at parse; the `ins<range>inv` form is required |

## `inversion.md:20` — large genomic inversions on RNA give deletion or delins

> since exon splice signals will be inverted, large genomic inversions on the RNA level usually give [deletion](deletion.md) or [deletion-insertion (delins)](delins.md) variants

Ferro: descriptive biology, not a rewrite. A large genomic inversion seen on RNA usually surfaces as
a deletion or delins because it disrupts splice signals. Ferro applies no automatic recasting, so
there's no verdict row.

## `inversion.md:21-22` — inversions are not used on the protein level

> - inversions are not used on protein level.
>   Depending on the (predicted) consequences of an inversion on protein level, changes are usually described as either a **delins** or a **frameshift**.

Ferro: the `p.` axis has no `inv` edit type; a protein consequence of an inversion is written as a
delins or frameshift. Ferro emits no protein `inv` — CONFIRM-by-inspection, no verdict row.

## `inversion.md:26-27` — Example: a four-nucleotide inversion

> - **`r.177_180inv`**<br>
>   inversion of nucleotides `r.177` to `r.180`, changing `..agg`<code class="del">cuga</code>`uu..` to `..agg`<code class="ins">ucag</code>`uu..`.

Ferro: the spec's worked example; its reverse-complement math checks out (reverse of `cuga` is
`aguc`, complement `ucag`). A four-nucleotide whole-span reverse complement types as one `inv`
(`whole-span-reverse-complement-types-as-inv`, cited at `:5`). The accession-less spec spelling isn't
table-executable; the executable twin on the slice is a fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.177_180inv` | recommended | self | executable twin of the spec's bare `r.177_180inv` (`gcaa` → `uugc` on the slice); already `inv` and 3'-most, a fixed point |

## `inversion.md:29-30` — Example: a 304-nucleotide inversion

> - **`r.203_506inv`**<br>
>   inversion of the 304 nucleotides from position `r.203` to `r.506`.

Ferro: the spec's larger worked example — a 304-nucleotide range inversion. No bases are given, so
nothing to recompute; well-formed and well over the one-nucleotide floor (`:16`). Shown on the slice
as a fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.203_506inv` | recommended | self | executable twin of the spec's bare `r.203_506inv` (222 of 304 columns change on the slice); a fixed point |

## `inversion.md:34-37` — Discussion: the complement is not an inversion

<details class="ss-spec"><summary>Spec discussion — is the complement an inversion?</summary>

> !!! note "Is the change `aagc` to `uucg` an inversion?"
>
>     No, an inversion would change `aagc` to `gcuu`, its **reverse-complement**.
>     `uucg` is only the **complement** of `aagc`.

</details>

Ferro: an inversion is the **reverse** complement, so `aagc` → `gcuu`; `uucg` is only the complement
(bases flipped, order kept). Ferro already uses the reverse complement — CONFIRM-by-inspection.

## `inversion.md:39-42` — Discussion: the reverse is not an inversion

<details class="ss-spec"><summary>Spec discussion — is the reverse an inversion?</summary>

> !!! note "Is the change `aagc` to `cgaa` an inversion?"
>
>     No, an inversion would change `aagc` to `gcuu`, its **reverse-complement**.
>     `cgaa` is only the **reverse** of `aagc`.

</details>

Ferro: the mirror image — `cgaa` is only the **reverse** (order flipped, bases kept), not the reverse
complement. Together these two fence the definition against its near-misses. Descriptive; no verdict row.

## `inversion.md:44-46` — Discussion: the old "o" strand marker is not used

<details class="ss-spec"><summary>Spec discussion — is the old "o" strand marker still used?</summary>

> !!! note "On the [old nomenclature website](http://www.HGVS.org/mutnomen/examplesRNA.html) (bottom), you had the example <code class="invalid">r.124_500delinsoAB053210.2:r.1289-365_1289-73</code>, i.e. the "o" indicating the inserted sequence `AB053210.2:r.1289-365_1289-73` was from the opposite transcriptional strand. Is the "o" still used?"
>
>     No, the "o" is not used, the insertion is in an inverted orientation, so "inv" should be used; <code class="invalid">r.124_500delins[AB053210.2:r.1289-365_1289-73inv]</code>.

</details>

Ferro: descriptive history. The deprecated `o` strand marker is gone; opposite-strand insert sequence
is spelled with an `inv`-suffixed payload (`delins[…inv]`). This concerns insertion-payload
orientation, not inversion typing; the cited spelling is on the foreign accession `AB053210.2`, so
there's no executable row.

See also → `RNA/insertion.md:18` and `RNA/duplication.md:18-21` (the `ins<range>inv`
inverted-duplication form, worked from the insertion and duplication sides), `RNA/substitution.md`
(where a one-nucleotide "inversion" lands, per `:16-17`), `DNA/inversion.md` (where the whole-span
reverse-complement → `inv` typing derivation is wired and guarded today).
