# Frameshift — ferro's reading

ferro's reading of the HGVS protein **frameshift** recommendations, clause by clause. New here? See
[How to read a page](../../reading-guide.md). Examples run against `NP_003997.1` (dystrophin) in
the committed slice.

## `frameshift.md:5` — definition: translation shifts frame

> translation shifts to another reading frame

Ferro: a frameshift names the first residue changed, the first new residue, and the position of the
new stop (`Ter#` / `*#`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Thr97ProfsTer23` | recommended | self | `Thr97`→`Pro`, new stop 23 codons downstream |
| `NP_003997.1:p.Glu5ValfsTer5` | recommended | self | `Glu5`→`Val`, new stop after five codons |

## `frameshift.md:20-22` — `*` is `Ter`, and an immediate stop is nonsense

> the shortest frameshift variant possible contains `fsTer2`

Ferro: the `*` stop shorthand is normalized to `Ter`. A frameshift needs at least `fsTer2`; an
*immediate* stop (`fsTer1`) is a nonsense substitution (`substitution.md:17`), not a frameshift.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Thr97Profs*23` | recommended | `NP_003997.1:p.Thr97ProfsTer23` | the `*` stop shorthand normalized to `Ter` |
| `NP_003997.1:p.Trp4TerfsTer1` | refused | — | an immediate stop is a nonsense substitution, not a frameshift |

## `frameshift.md:23` — the short format

> frameshifts can also be described using a **short format**; `p.Arg123fs`

Ferro: the short form gives only the first residue changed and `fs`; ferro accepts and preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Thr97fs` | recommended | self | the short form of the frameshift at `Thr97` |

## `frameshift.md:44-45` — no new stop codon is `Ter?`

> the new reading frame does not encounter a new translation termination (stop) codon

Ferro: when the shifted frame reaches no new stop, the position is `Ter?` / `*?`; the `*` form is
normalized to `Ter`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Asp327Argfs*?` | recommended | `NP_003997.1:p.Asp327ArgfsTer?` | no new stop reached; `*?` normalized to `Ter?` |

## `frameshift.md:47-49` — anchored at the first residue changed

> Since frameshift variants start with the **first amino acid changed**, the description

Ferro: the anchor is the first residue the shifted frame *alters* — not the first codon touched on
the DNA level. A frameshift anchored at an unchanged residue is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Asn150Asnfs*10` | refused | — | `Asn150` names an unchanged residue — anchor at the first changed residue |

See also → `substitution.md:17`, `delins.md:5`, `extension.md:5`.
