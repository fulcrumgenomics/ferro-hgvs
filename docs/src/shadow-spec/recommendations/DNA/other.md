# Other — ferro's reading

ferro's reading of the HGVS **other** recommendations — the "no change" (`=`) assertion,
methylation (`|gom`/`|lom`/`|met=`), and the mosaic (`=/`) and chimeric (`=//`) mixtures — clause
by clause, each spelling with the form ferro normalizes it to and a verdict on that output. New
here? See [How to read a page](../../reading-guide.md) for the verdicts, the table conventions,
and the recurring terms.

No ruling record in the ledger cites `other.md`: the notation families here are either mechanical
(the `=` identity and the methylation tokens) or governed — for mosaicism and chimerism — by
clauses `general.md:91-92` attributes to `substitution.md`, adjudicated on that page.

The one interesting behaviour is a mosaic whose member 3'-shifts: when a `del` member inside `=/`
shifts under the 3'rule while the `=` range does not, ferro abandons the compact form and repeats
the accession — the same finding the deletion page carries at `deletion.md:105-110`. It gets its
own section below.

Executable rows use `NM_004006.3`, the one transcript in the committed slice, over known bases:
`c.123` is `C`, `c.124` is `A`, and `c.5690_5697` is an A-stretch. The spec's own worked examples
sit on accessions the slice does not carry, so those rows are parse-only (`—`).

## `other.md:5-7` — definition: no change, mosaic, chimeric

> No change: a sequence was analysed, but no variant was detected.
> Mosaic: the occurrence in one individual of two or more cell populations, derived from a single zygote, with different sequences.
> Chimeric: the occurrence in one individual of two or more cell populations, derived from different zygotes, with different sequences.

Ferro: pure terminology — biological definitions of three terms, with no notation, position or
edit a normalizer can act on. Descriptive.

## `other.md:23-28` — no change: one nucleotide, and the self-substitution prohibition

> **NOTE**: the description <code class="invalid">c.123C>C</code> is not allowed.

Ferro: a tested, unchanged position is spelled `<pos>=` — there is exactly one legal spelling per
(position, "no change") pair. The self-substitution `c.123C>C` is `class="invalid"` and is not a
synonym for `c.123=` — strict rejects it, lenient repairs it to `c.123=` (the same behaviour
adjudicated on `substitution.md:19`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32867908=` | recommended | — | the spec's own genomic "no change" example (parse-only) |
| `LRG_199t1:c.123=` | recommended | — | the spec's own coding example (parse-only) |
| `NM_004006.3:c.123=` | recommended | self | executable twin — `c.123` is reference `C`, so the identity assertion is already a fixed point |
| `NM_004006.3:c.123C>C` | recommended | `NM_004006.3:c.123=` | the `class="invalid"` self-substitution repaired to the recommended `=` (matches `substitution.md:19`) |

## `other.md:31-38` — no change: several nucleotides, and the identity allele

> a screen was performed showing nucleotides `c.123`, `c.456`, and `c.789` (all on the same allele) were identical to the coding DNA reference (the nucleotides were not changed).

Ferro: the identity form extends to a range (`<start>_<end>=`) and to a multi-member allele of
independent "no change" assertions (`[<pos>=;<pos>=;…]`). Each is syntactically unambiguous, so
each is a fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32867908_32867923=` | recommended | — | the spec's own genomic range example (parse-only) |
| `LRG_199t1:c.123_145=` | recommended | — | the spec's own coding range example (parse-only) |
| `NM_004006.2:c.[123=;456=;789=]` | recommended | — | the spec's own multi-member identity allele (parse-only) |
| `NM_004006.3:c.123_124=` | recommended | self | executable range twin — identity over `c.123` (`C`) and `c.124` (`A`), a fixed point |
| `NM_004006.3:c.[123=;124=]` | recommended | self | executable multi-member identity allele on the slice's transcript, preserved |

## `other.md:42-46` — methylation (`|gom`/`|lom`/`|met=`)

> the sequence from position `g.1999904` to `g.1999946` showed a gain of methylation (`|gom`).

Ferro: the three methylation-state tokens — gain (`|gom`), loss (`|lom`) and normal (`|met=`) —
parse and round-trip, so ferro accepts and preserves each. A methylation marker carries no
reference bases to shift, so each is a fixed point. This is the one notation family this doc
genuinely originates (`general.md:93` names `other.md` as its authority), yet it has no grammar
entry in `syntax.yaml` and no committed worked-example pin, so its round-trip is unverified by any
guard. The spec's three examples — `NC_000011.10:g.1999904_1999946|gom` (gain), `…|lom` (loss) and
`…|met=` (normal) — sit on a foreign accession, so they are parse-only and shown inline rather than
in a table (`|` is not a markdown-table-safe character).

## `other.md:50-52` — mosaicism (`=/`): the compact reference-first form

> **NOTE**: irrespective of the frequency in which each nucleotide was found, the reference is always described first.

Ferro: the mosaic operator `/` joins two members of a single mixed sample; the recommendations
write the reference member first, in the compact form `<pos>=/<edit>` where the `=` identity
carries the position and the bare edit inherits it. For a two-member substitution mixture, ferro
reorders a variant-first input to that form
([#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)), governed on
`substitution.md:47-49`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.85=/T>C` | recommended | — | the spec's own mosaic example, reference-first (parse-only) |
| `NM_004006.3:c.123=/C>T` | recommended | self | executable mosaic-substitution twin — reference `=` written first per the NOTE, single-base member so no shift |
| `NM_004006.3:c.123C>T/=` | recommended | `NM_004006.3:c.123=/C>T` | variant-first shorthand; ferro reorders it reference-first, as the NOTE and `substitution.md:49` write it ([#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)) |

## `other.md:49` — a mosaic whose `del` member 3'-shifts (the finding)

> **mosaicism**

The compact mosaic form `<pos>=/<edit>` shares one range between the `=` identity and the edit
member. When that edit is a `del` sitting in a shiftable run, the 3'rule moves the `del` member but
not the `=` range — so the shared range no longer describes both members, and ferro abandons the
compact form and repeats the accession after `=/` rather than shifting the shared range with the
deletion. Same finding as `deletion.md:105-110`; no tracking issue filed yet.

The two fixed-point rows show the compact form is preserved when the member is already 3'-most;
the last row is the finding — the `del` member shifts off the shared range.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5695_5697=/del` | recommended | self | executable mosaic twin — the `del` member is already 3'-most within the A-stretch, a fixed point, so the compact form holds |
| `NM_004006.3:c.5695_5697=//del` | recommended | self | executable chimeric twin — same 3'-most member, compact form preserved |
| `NM_004006.3:c.5690_5692=/del` | conformant | `NM_004006.3:c.5690_5692=/NM_004006.3:c.5695_5697del` | **the finding** — the `del` member 3'-shifts to `c.5695_5697` while the `=` range stays at `c.5690_5692`, so ferro abandons the compact form and repeats the accession after `=/`. Re-parses (conformant, not a bug), but the recommended spelling is the compact `c.5695_5697=/del`. No tracking issue yet |

## `other.md:55-57` — chimerism (`=//`)

> **NOTE**: irrespective of the frequency in which each nucleotide was found, the reference is always described first.

Ferro: identical in structure and jurisdiction to mosaicism — `//` joins two members of a mix of
cell populations, governed by the same `substitution.md` authority as `/`. Ferro accepts and
preserves the compact chimeric form, and reorders the variant-first substitution shorthand to it
([#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.85=//T>C` | recommended | — | the spec's own chimeric example (parse-only) |
| `NM_004006.3:c.123=//C>T` | recommended | self | executable chimeric-substitution twin — reference `=` written first per the NOTE |
| `NM_004006.3:c.123C>T//=` | recommended | `NM_004006.3:c.123=//C>T` | variant-first shorthand; ferro reorders it reference-first, as the NOTE writes it ([#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)) |
