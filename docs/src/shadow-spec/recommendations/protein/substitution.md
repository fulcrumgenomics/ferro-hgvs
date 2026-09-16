# Substitution — ferro's reading

ferro's reading of the HGVS protein **substitution** recommendations, clause by clause — each
spelling with the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

Protein examples run against `NP_003997.1` (dystrophin, the translation of `NM_004006.3`) in the
committed slice, so residues and positions are real. Per the spec, a protein change **should also**
be described on the DNA level — these pages read the `p.` axis in isolation.

## `substitution.md:5` — definition: one amino acid for one

> **one** amino acid is replaced by **one** other amino acid

Ferro: a substitution is exactly 1→1; a single residue replaced by two or more is a delins.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp24Cys` | recommended | self | canonical missense, `Trp24`→`Cys` |
| `NP_003997.1:p.TrpVal24CysArg` | refused | — | two-for-two is not a substitution — write a delins (`delins.md:5`) |

## `substitution.md:16` — predicted consequences go in parentheses

> should be given in parentheses

Ferro: when the protein change is inferred (no protein sequenced), it is parenthesized; ferro
accepts and preserves the predicted form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.(Trp24Cys)` | recommended | self | predicted missense (e.g. inferred from DNA-level data) |

## `substitution.md:17` — a nonsense variant is a substitution

> changing an amino acid to a translation termination (stop) codon, is described as a **substitution**

Ferro: a residue-to-stop change is a substitution to `Ter`, not a C-terminal deletion. The `*`
shorthand for the stop codon is normalized to `Ter`. An *immediate* stop is a nonsense
substitution, never a frameshift.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp24Ter` | recommended | self | nonsense: `Trp24` to a stop codon |
| `NP_003997.1:p.Trp24*` | recommended | `NP_003997.1:p.Trp24Ter` | the `*` stop shorthand normalized to `Ter` |
| `NP_003997.1:p.Trp4TerfsTer1` | refused | — | an immediate stop is a nonsense substitution, not a frameshift (`substitution.md:20`, `frameshift.md:22`) |

## `substitution.md:24` — a tested, unchanged residue is `=`

> found **not changed** (silent) are described as `p.Cys123=`

Ferro: a screened, unchanged position is `=`; `p.=` is the distinct claim that the *entire* coding
region was analysed and nothing changed.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Cys188=` | recommended | self | `Cys188` screened, unchanged |
| `NP_003997.1:p.=` | recommended | self | whole-protein "no change" — a distinct claim from a position `=` |

Ferro also accepts the older `p.Cys188Cys` spelling (a residue written as unchanged into itself) and
**preserves** it — it does not collapse it to `p.Cys188=`. The DNA analogue `c.…C>C` is disallowed
and the DNA axis rewrites it to `=`; the protein axis does not yet, so `p.Cys188Cys` stays as
written (a known limitation, not the recommended form —
[#2220](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2220)).

## `substitution.md:45-49` — the translation initiation codon

> no protein is produced

Ferro: a variant in the start codon is written `p.0` (no protein), `p.0?` (predicted none), or
`p.(Met1?)` (consequence unknown) — never as a substitution of `Met1`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.0` | recommended | self | no protein produced |
| `NP_003997.1:p.0?` | recommended | self | predicted no protein |
| `NP_003997.1:p.(Met1?)` | recommended | self | consequence on the protein unknown |
| `NP_003997.1:p.Met1Thr` | refused | — | an initiation-codon change is not a substitution (`substitution.md:49`); use `p.0` / `p.0?` / `p.(Met1?)` |

## `substitution.md:76-78` — an uncertain set of replacing residues

> changed to an `Ala`, `Ser`, or `Cys`

Ferro: the `^` "one of" set is a valid uncertainty; ferro parses and preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.(Gly56Ala^Ser^Cys)` | recommended | self | `Gly56` changed to one of `Ala`, `Ser`, `Cys` |

## `substitution.md:80-84` — mosaic, reference first

> the reference is always described first

Ferro: a mosaic mixture uses `=/`, with the reference allele written first; ferro accepts and
preserves the order.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp24=/Cys` | recommended | self | mosaic: reference `Trp24=` first, then `Cys` |
| `NP_003997.1:p.(Trp24=/Cys)` | recommended | self | the predicted-consequence form of the same mosaic |

## `substitution.md:88-92` — the "polymorphism" slash is not HGVS

> all substitutions are described as `NP_003997.1:p.Gln2366Lys`

Ferro: the historical `p.2366Gln/Lys` "polymorphism" form is not valid HGVS; ferro rejects it at
parse. The only conformant form names the reference and the replacing residue.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Lys2366Gln` | recommended | self | the correct substitution form |
| `NP_003997.1:p.2366Gln/Lys` | refused | — | the disallowed "polymorphism" slash — rejected at parse |

## `substitution.md:100-104` — "any amino acid" is `Xaa`

> we suggest to use `Xaa` three-letter amino acid code only

Ferro: the IUPAC symbol for "any amino acid" is `Xaa` (three-letter), not `X` — `X` historically
meant a stop codon. The three-letter code is the preferred spelling generally: a one-letter code is
valid HGVS but not preferred, and ferro rewrites it to three-letter (strict mode refuses one-letter).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Arg782Xaa` | recommended | self | `Arg782` to any amino acid, written `Xaa` |
| `NP_003997.1:p.W24C` | recommended | `NP_003997.1:p.Trp24Cys` | one-letter input rewritten to the preferred three-letter form (strict mode refuses one-letter) |

See also → `delins.md:5`, `extension.md:5`, `frameshift.md:22`.
