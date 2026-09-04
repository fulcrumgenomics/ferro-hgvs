# Repeated Sequences — ferro's reading

ferro's reading of the HGVS protein **repeated sequence** recommendations, clause by clause. New
here? See [How to read a page](../../reading-guide.md). Examples run against `NP_003997.1`
(dystrophin) in the committed slice, whose N-terminus carries a `TrpTrp` (positions 3–4) and a
`LysLys` (positions 18–19) tandem pair.

## `repeated.md:5` — definition: a unit repeated in tandem

> a segment of **one or more** amino acids (the repeat unit) is present several times, one after the other

Ferro: a repeat names the first residue of the unit and the copy count in `[n]`; ferro accepts and
preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp3[2]` | recommended | self | the reference `TrpTrp` stated as two copies from position 3 |

## `repeated.md:20-23` — the repeat count and the del / dup form are both legal

> when the repeat is variable in the population and the reference sequence has 10 units, the description `p.Ala2[9]` is preferred over `p.Ala11del`

Ferro: a change to a repeat count can be spelled as a repeat (`[n]`) or as a `del` / `dup` — both are
legal, and the spec prefers the repeat form **only when the repeat is variable in the population**, a
fact ferro cannot derive from sequence. So ferro preserves whichever form it is given rather than
rewriting one into the other. The `del` / `dup` spellings of the same change, which ferro *does*
treat as recommended, are on the [deletion](deletion.md) and [duplication](duplication.md) pages.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp3[1]` | conformant | self | a contraction to one copy — the same protein as `p.Trp4del`; ferro preserves the repeat form |
| `NP_003997.1:p.Lys18[3]` | conformant | self | an expansion to three copies — the same protein as `p.Lys19dup`; ferro preserves the repeat form |

## `repeated.md:25-26` — two alleles

> present in 10 copies on one allele and 11 copies on the other allele

Ferro: per-allele repeat counts use the `[m];[n]` form; ferro preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp3[2];[3]` | recommended | self | two copies on one allele, three on the other |

See also → `deletion.md:19`, `duplication.md:22`, `alleles.md:5`.
