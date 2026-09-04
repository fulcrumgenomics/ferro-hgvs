# Deletion-Insertion — ferro's reading

ferro's reading of the HGVS protein **deletion-insertion (delins)** recommendations, clause by
clause. New here? See [How to read a page](../../reading-guide.md). Examples run against
`NP_003997.1` (dystrophin) in the committed slice.

## `delins.md:5` — definition: residues replaced by others

> one or more amino acids are replaced by one or more other amino acids **and which is not** a substitution or frameshift

Ferro: a delins replaces one or more residues with one or more others, where the change is neither a
substitution nor a frameshift.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Cys10delinsTrpVal` | recommended | self | `Cys10` replaced by `TrpVal` |
| `NP_003997.1:p.Trp24_Val25delinsCysArg` | recommended | self | `Trp24Val25` replaced by `CysArg` |

## `delins.md:17` — one-for-one is a substitution, not a delins

> when **one** amino acid is replaced by **one** other amino acid, the change is a

Ferro: a single residue replaced by a single residue is a substitution; ferro rewrites the 1→1
delins spelling to it. A single residue replaced by *two or more* is a genuine delins.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp24delinsCys` | recommended | `NP_003997.1:p.Trp24Cys` | 1→1 delins normalized to the substitution |
| `NP_003997.1:p.Trp24delinsCysArg` | recommended | self | 1→2 is a genuine delins |

## `delins.md:18-20` — two or more consecutive changes are one delins

> the description `p.Arg76_Cys77delinsSerTrp` is correct

Ferro: adjacent per-residue changes are written as one delins, derived from the resulting sequence —
the split spelling is marked "not correct" at separation zero.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
<!-- why:END:delins-adjacent-members-when-both-consume-reference -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[Trp24Ser;Val25Arg]` | recommended | `NP_003997.1:p.Trp24_Val25delinsSerArg` | two adjacent substitutions merged into one delins |
| `NP_003997.1:p.Trp24_Val25delinsSerArg` | recommended | self | the recommended merged form |
| `NP_003997.1:p.TrpVal24CysArg` | refused | — | the multi-residue substitution spelling — write a delins |

## `delins.md:62-64` — a separated pair is described individually

> the variant is not described as `p.Ser44_Trp46delinsArgLeuArg`

Ferro: two changes separated by an unchanged residue stay separate members; they are *not* fused
into a delins that would falsely claim the residue between them changed.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.[Leu44Arg;Asp46Arg]` | recommended | self | separated by `Gln45` — kept as two members, not one delins |

## `delins.md:22` — the 3′ rule reduces a delins to a deletion

> the **most C-terminal position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**)

Ferro: when a delins's payload coincides with part of the reference, the net change is a plain
deletion; ferro reduces it and applies the 3′ rule.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp3_Trp4delinsTrp` | recommended | `NP_003997.1:p.Trp4del` | deleting `TrpTrp` and inserting `Trp` is a one-`Trp` deletion, 3′-placed |

## `delins.md:43-45` — a stop in the payload

> amino acids after the translation termination codon are **not** listed

Ferro: an inserted peptide ends at its first stop. A delins whose payload *ends* at `Ter` is
accepted; residues written *after* the `Ter` make the description invalid, and ferro refuses it. (An
*immediate* stop is a nonsense substitution — `substitution.md:17`.)

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp24_Val25delinsLeuTer` | recommended | self | payload ends at `Ter` — accepted |
| `NP_003997.1:p.Glu5_Glu6delinsTerGluAsp` | refused | — | residues listed after `Ter` — rejected at parse (`delins.md:45`) |

See also → `substitution.md:17`, `frameshift.md:5`, `extension.md:5`.
