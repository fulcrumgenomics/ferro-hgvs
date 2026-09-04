# Other — ferro's reading

ferro's reading of the HGVS **other** recommendations — the "no change" (`=`) assertion,
methylation (`|gom`/`|lom`/`|met=`), and the mosaic (`=/`) and chimeric (`=//`) mixtures — clause
by clause, each spelling with the form ferro normalizes it to and a verdict on that output. New
here? See [How to read a page](../../reading-guide.md) for the verdicts, the table conventions,
and the recurring terms.

**No Why block appears on this page.** No ruling record in the ledger cites `other.md`, and the
notation families this doc carries are either mechanical (the `=` identity and the methylation
tokens) or governed — for mosaicism and chimerism — by clauses `general.md:91-92` attributes to
`substitution.md`, adjudicated on that page. This doc republishes the same two examples under its
own headers; the reasoning lives once, on the substitution page.

Most of this page is CONFIRM-by-inspection against the spec text and the shipped code: the `=`
identity form and its range and multi-member-allele variants, the self-substitution prohibition
(`c.123C>C` is not a synonym for `c.123=`), and the three methylation tokens are all mechanical
parser or round-trip behaviour with no clause in tension. The **one interesting behaviour** is a
mosaic whose member 3'-shifts: when a `del` member inside `=/` shifts under the 3'rule while the
`=` range does not, ferro abandons the compact form and repeats the accession. It is the same
finding the deletion page carries at `deletion.md:105-110`; no tracking issue is filed yet, and it
gets its own section below.

Executable rows use `NM_004006.3`, the one transcript in the committed slice, over known bases:
`c.123` is `C`, `c.124` is `A`, and `c.5690_5697` is the 8-nucleotide A-stretch (`AAAAAAAA`). The
spec's own worked examples sit on `LRG_199t1`, `NM_004006.2` and foreign genomic accessions
(`NC_000023.11`, `NC_000011.10`), none of which the slice carries, so those rows are parse-only
(`—`) — ferro cannot read their bases here. Rows built on constructed mosaic/chimeric edits over
`NM_004006.3` show ferro's actual normalized output, verified against the bless harness.

## `other.md:5-7` — definition: no change, mosaic, chimeric

> No change: a sequence was analysed, but no variant was detected.
> Mosaic: the occurrence in one individual of two or more cell populations, derived from a single zygote, with different sequences.
> Chimeric: the occurrence in one individual of two or more cell populations, derived from different zygotes, with different sequences.

Ferro: pure terminology — biological definitions of three terms, with no notation, position or
edit a normalizer can act on. Descriptive; no verdict row.

## `other.md:23-28` — no change: one nucleotide, and the self-substitution prohibition

> **NOTE**: the description <code class="invalid">c.123C>C</code> is not allowed.

Ferro: a tested, unchanged position is spelled `<pos>=`; there is exactly one legal spelling per
(position, "no change") pair, so nothing for the normalizer to choose. The self-substitution
`c.123C>C` is `class="invalid"` and is not a synonym for `c.123=` — strict rejects it, lenient
repairs it to `c.123=` (the same behaviour adjudicated on `substitution.md:19`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32867908=` | recommended | — | the spec's own genomic "no change" example (foreign accession — parse-only here) |
| `LRG_199t1:c.123=` | recommended | — | the spec's own coding example (`LRG_199t1` — parse-only here) |
| `NM_004006.3:c.123=` | recommended | self | executable twin — `c.123` is reference `C`, so the identity assertion is already a fixed point |
| `NM_004006.3:c.123C>C` | recommended | `NM_004006.3:c.123=` | the `class="invalid"` self-substitution repaired to the recommended `=` (matches `substitution.md:19`) |

## `other.md:31-38` — no change: several nucleotides, and the identity allele

> a screen was performed showing nucleotides `c.123`, `c.456`, and `c.789` (all on the same allele) were identical to the coding DNA reference (the nucleotides were not changed).

Ferro: the identity form extends to a range (`<start>_<end>=`) and to a multi-member allele of
independent "no change" assertions (`[<pos>=;<pos>=;…]`). Each is syntactically unambiguous, so
each is a fixed point.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.32867908_32867923=` | recommended | — | the spec's own genomic range example (foreign accession — parse-only here) |
| `LRG_199t1:c.123_145=` | recommended | — | the spec's own coding range example (`LRG_199t1` — parse-only here) |
| `NM_004006.2:c.[123=;456=;789=]` | recommended | — | the spec's own multi-member identity allele (`NM_004006.2`, wrong transcript version — parse-only here) |
| `NM_004006.3:c.123_124=` | recommended | self | executable range twin — identity over `c.123` (`C`) and `c.124` (`A`), a fixed point |
| `NM_004006.3:c.[123=;124=]` | recommended | self | executable multi-member identity allele on the slice's transcript, preserved |

## `other.md:42-46` — methylation (`|gom`/`|lom`/`|met=`)

> the sequence from position `g.1999904` to `g.1999946` showed a gain of methylation (`|gom`).

Ferro: the three methylation-state tokens — gain (`|gom`), loss (`|lom`) and normal (`|met=`) —
parse and round-trip (`src/hgvs/edit.rs`'s `MethylationStatus`), so ferro accepts each and
preserves it. A methylation marker carries no reference bases to shift, so normalization has
nothing to move; each is a fixed point. Note the coverage gap the clean-room review recorded:
this is the one notation family this doc genuinely originates (`general.md:93` names `other.md`
as its authority), yet it has no `syntax.yaml` grammar entry and no committed worked-example pin,
so its round-trip is unverified by any guard. The spec's three examples —
`NC_000011.10:g.1999904_1999946|gom` (gain), `…|lom` (loss) and `…|met=` (normal) — are all on
that foreign genomic accession, so they are parse-only here and are shown inline rather than in an
executable table (the `|` methylation delimiter is not a markdown-table-safe character).

## `other.md:50-52` — mosaicism (`=/`): the compact reference-first form

> **NOTE**: irrespective of the frequency in which each nucleotide was found, the reference is always described first.

Ferro: the mosaic operator `/` joins two members of a single mixed sample; the recommendations
write the reference member first, in the compact form `<pos>=/<edit>` where the `=` identity
carries the position and the bare edit inherits it. The reference-first ordering the NOTE requires
is not enforced against a variant-first input: `substitution.md:47-49` is the governing text, and
ferro's non-reordering of the variant-first shorthand (`<edit>/=`) is tracked there by
[#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.85=/T>C` | recommended | — | the spec's own mosaic example, reference-first (`LRG_199t1` — parse-only here) |
| `NM_004006.3:c.123=/C>T` | recommended | self | executable mosaic-substitution twin — reference `=` written first per the NOTE, single-base member so no shift |
| `NM_004006.3:c.123C>T/=` | conformant | self | valid, but the reference should be written first (NOTE); ferro does not reorder the variant-first shorthand to `c.123=/C>T` (same limitation as `substitution.md:49`, tracked by [#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)) |

## `other.md:49` — a mosaic whose `del` member 3'-shifts (the finding)

> **mosaicism**

The one behaviour on this page worth flagging. The compact mosaic form `<pos>=/<edit>` shares one
range between the `=` identity and the edit member. When that edit is a `del` sitting in a
shiftable run, the 3'rule moves the `del` member but **not** the `=` range — so the single shared
range no longer describes both members, and ferro abandons the compact form and **repeats the
accession** after `=/` rather than shifting the shared range with the deletion. This is the same
finding the deletion page carries at `deletion.md:105-110`; no tracking issue is filed for
the mosaic/chimeric case yet.

The two fixed-point rows show the compact form is preserved when the member is already 3'-most;
the last row is the finding — the `del` member shifts off the shared range.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5695_5697=/del` | recommended | self | executable mosaic twin — the `del` member is the 3'-most three of the A-stretch, already a fixed point, so the compact form holds |
| `NM_004006.3:c.5695_5697=//del` | recommended | self | executable chimeric twin — same 3'-most member, compact form preserved |
| `NM_004006.3:c.5690_5692=/del` | conformant | `NM_004006.3:c.5690_5692=/NM_004006.3:c.5695_5697del` | **the finding** — the `del` member 3'-shifts to `c.5695_5697` while the `=` range stays at `c.5690_5692`, so ferro abandons the compact form and repeats the accession after `=/`. Re-parses (so conformant, not a bug), but the recommended spelling is the compact `c.5695_5697=/del` (both ranges name the same change). No tracking issue yet |

## `other.md:55-57` — chimerism (`=//`)

> **NOTE**: irrespective of the frequency in which each nucleotide was found, the reference is always described first.

Ferro: identical in structure and jurisdiction to mosaicism — `//` joins two members of a mix of
cell populations, and `general.md:92` names `substitution.md` as the authority for `//` exactly as
it does for `/`. The spec's spelling is reproduced byte-for-byte at `substitution.md:51-53`, and
the ordering NOTE is the same. Ferro accepts and preserves the compact chimeric form; the
variant-first shorthand is again not reordered (tracked by
[#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.85=//T>C` | recommended | — | the spec's own chimeric example (`NM_004006.2`, wrong transcript version — parse-only here) |
| `NM_004006.3:c.123=//C>T` | recommended | self | executable chimeric-substitution twin — reference `=` written first per the NOTE |
| `NM_004006.3:c.123C>T//=` | conformant | self | valid, but the reference should be written first (NOTE); ferro does not reorder the variant-first shorthand to `c.123=//C>T` (tracked by [#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)) |
