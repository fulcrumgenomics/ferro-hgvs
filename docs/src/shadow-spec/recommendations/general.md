# General — ferro's reading

ferro's reading of the HGVS **general recommendations** — the cross-cutting rules and the special
characters. New here? See [How to read a page](../reading-guide.md). Examples run against
`NM_004006.3` / `NP_003997.1` in the committed slice. The per-edit-type rules live on the
[DNA](DNA/substitution.md), [RNA](RNA/substitution.md), and [protein](protein/substitution.md) pages.

## `general.md:10-11` — describe at the DNA level; other levels in addition

> all variants should be described at the most basic level, **the DNA level**

Ferro: the DNA-level description is primary; RNA and protein descriptions are given in addition, and
a predicted one is parenthesized. ferro parses and preserves each.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76A>G` | recommended | self | the DNA-level description — the basic level |
| `NP_003997.1:p.(Trp24Cys)` | recommended | self | a predicted protein consequence, given in addition |

## `general.md:22-30` — a letter prefix is mandatory

> a **letter prefix** is mandatory to indicate the type of reference sequence used

Ferro: every description carries a coordinate prefix — `c` / `g` / `m` / `n` / `o` / `p` / `r`. An
unknown prefix is refused at parse.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76A>G` | recommended | self | a valid `c.` (coding DNA) prefix |
| `NM_004006.3:x.76A>G` | refused | — | `x.` is not an accepted prefix |

## `general.md:55-57` — prioritisation, and overlapping allele members

> the preferred description is: (1) substitution, (2) deletion, (3) inversion, (4) duplication, (5) insertion

Ferro: when several types could describe one change, the earlier type wins — a tandem-copy insertion
is written as a duplication (`DNA/insertion.md:17`). A *separate* rule bars two members of one allele
from claiming **overlapping reference territory**: `c.[762_768del;767_774dup]` has `767`–`768` in
both members, so it is refused. That refusal is governed by the allele definition (`DNA/alleles.md:5`),
**not** by the prioritisation order.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[conflicting-member-geometry-refusal-scope](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two members of one allele that claim intersecting reference territory — nested, overlapping, or two insertions at one interbase — are refused, whatever edit types they render as.
<!-- why:END:conflicting-member-geometry-refusal-scope -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[762_768del;767_774dup]` | refused | — | members share reference bases `767`–`768` — overlapping territory (`DNA/alleles.md:5`) |

## `general.md:90-91` — the `=` (unchanged) and `/` (mosaic) characters

> `=` (equals) is used to indicate a sequence was tested but found unchanged

Ferro: `=` marks a tested-unchanged position; `/` marks a mosaic mixture, written reference-first.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76=` | recommended | self | `c.76` tested, unchanged |
| `NM_004006.3:c.76=/A>G` | recommended | self | mosaic, reference (`=`) written first |
| `NM_004006.3:c.76A>G/=` | conformant | self | valid, but the reference should be written first; ferro does not yet reorder it ([#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034)) |

## `general.md:95` — spaces are not permitted

> spaces are *not* permitted in any HGVS description

Ferro: a space anywhere in a description makes it invalid; ferro refuses it at parse.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76A>G` | recommended | self | no spaces |
| `NM_004006.3:c.76 A>G` | refused | — | an internal space — not a valid HGVS description |

See also → [DNA insertion](DNA/insertion.md) (dup-over-insertion prioritisation), [uncertain descriptions](uncertain.md), [publication checklist](checklist.md).
