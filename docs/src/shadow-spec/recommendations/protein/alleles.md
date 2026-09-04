# Alleles — ferro's reading

ferro's reading of the HGVS protein **allele** recommendations, clause by clause. New here? See
[How to read a page](../../reading-guide.md). Examples run against `NP_003997.1` (dystrophin) in
the committed slice. For the shared allele grammar, see also [DNA alleles](../DNA/alleles.md).

## `alleles.md:17` — two variants in cis: `p.[a;b]`

> derive from **one chromosome** (in cis)

Ferro: two variants on one chromosome are one bracketed allele, `;`-separated; each member is
normalized in place.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[Ser68Arg;Asn594del]` | recommended | self | `Ser68Arg` and `Asn594del` in cis |

## `alleles.md:18` — two variants in trans: `p.[a];[b]`

> derive from **different chromosomes** (in trans)

Ferro: variants on different chromosomes are separate bracketed alleles; a member that can shift is
placed by the 3′ rule like any other.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[Ser68Arg];[Asn594del]` | recommended | self | compound heterozygote, in trans |
| `NP_003997.1:p.[Trp3del];[Ser68Arg]` | recommended | `NP_003997.1:p.[Trp4del];[Ser68Arg]` | the first member's deletion is 3′-shifted to `Trp4` |

## `alleles.md:16` + `alleles.md:34` — predicted parentheses go *inside* the brackets

> should be given in parentheses inside the square brackets

Ferro: for a predicted allele the parentheses sit *inside* the square brackets, around each member;
the outside-parentheses form is not valid HGVS and is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[(Ser68Arg;Asn594del)]` | recommended | self | predicted, parentheses inside the brackets |
| `NP_003997.1:p.([Ser68Arg;Asn594del])` | refused | — | parentheses outside the brackets — rejected at parse |

## `alleles.md:19` — unknown phase: `a(;)b`, no brackets

> it is **not known** whether these derive from one chromosome (in cis) or from different chromosomes (in trans)

Ferro: when the phase is unknown, the `(;)` separator is used *without* brackets; ferro preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.(Ser68Arg)(;)(Asn594del)` | recommended | self | phase unknown — `(;)`, no allele brackets |

## `alleles.md:58-61` — a wild-type second allele is not `p.=`

> is different since it indicates the entire protein reference sequence was analysed

Ferro: `[Ser68=]` (the second allele carries the reference residue at that position) is a different
claim from `[=]` (the whole reference sequence was analysed, nothing found). Both are preserved,
distinctly. A male X-linked case with no second allele uses `[0]`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[Ser68Arg];[Ser68=]` | recommended | self | second allele wild-type at `Ser68` |
| `NP_003997.1:p.[Ser68Arg];[=]` | recommended | self | whole reference analysed — a distinct claim from `[Ser68=]` |
| `NP_003997.1:p.[Ser68Arg];[0]` | recommended | self | no protein from the second allele (e.g. X-linked, male) — `alleles.md:104-108` |

## `alleles.md:21` — one allele encoding two proteins: the `,` separator

> the variants are separated using a `,`

Ferro: two proteins from *one* DNA-level variant (via two transcripts) are `,`-separated inside one
bracket; ferro preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[Lys31Asn,Val25_Lys31del]` | recommended | self | two proteins from one allele — `,`-separated |

See also → `substitution.md:24`, `deletion.md:19`, and the [DNA allele grammar](../DNA/alleles.md).
