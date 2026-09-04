# Deletion — ferro's reading

ferro's reading of the HGVS protein **deletion** recommendations, clause by clause. New here? See
[How to read a page](../../reading-guide.md) for the verdicts and table conventions. Examples run
against `NP_003997.1` (dystrophin) in the committed slice.

## `deletion.md:5` — definition: one or more amino acids removed

> one or more amino acids are not present (deleted)

Ferro: a deletion removes one or more residues between the start and stop codon.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Val7del` | recommended | self | a single-residue deletion |
| `NP_003997.1:p.(Val7del)` | recommended | self | the predicted-consequence form |
| `NP_003997.1:p.Lys23_Val25del` | recommended | self | a three-residue deletion, `Lys23` to `Val25` |

## `deletion.md:17-18` — two different positions, listed 5′ to 3′

> should contain **two different** positions, i.e. `Cys76_Glu79`, not `Cys76_Cys76`

Ferro: a range names two *different* endpoints from 5′ to 3′; a range whose endpoints coincide is
normalized to the single-residue form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Val7_Val7del` | recommended | `NP_003997.1:p.Val7del` | identical endpoints collapse to the single-residue spelling |

## `deletion.md:19` — the 3′ rule

> the **most C-terminal position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**)

Ferro: in a run of identical residues (or a tandem repeat) the deletion is placed at the **most
C-terminal** position. In the reference `MetLeuTrpTrpGlu`, deleting either `Trp` yields the same
protein, so ferro writes the 3′-most.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp3del` | recommended | `NP_003997.1:p.Trp4del` | 3′ rule: shifted to the C-terminal `Trp4` |
| `NP_003997.1:p.Trp4del` | recommended | self | already the 3′-most residue of the `TrpTrp` run |

## `deletion.md:20-21` — a nonsense variant is not a C-terminal deletion

> is described as a **substitution** (`p.Trp26Ter` or `p.Trp26*`

Ferro: a residue-to-stop change is a nonsense **substitution** (`substitution.md:17`), never a
deletion of the C-terminal end (`p.Trp24_Lys1623del`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp24Ter` | recommended | self | the nonsense form — not a C-terminal deletion |

## `deletion.md:81-83` — mosaic deletion

> a mosaic case where at amino acid position `7`

Ferro: a mosaic deletion uses `=/del`, reference first; ferro preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Val7=/del` | recommended | self | mosaic: reference `Val7=` first, then the deletion |
| `NP_003997.1:p.(Val7=/del)` | recommended | self | the predicted-consequence form |

## `deletion.md:87-89` — a size suffix is not a deletion length

> a deletion of more than one residue should mention the first and last residue deleted

Ferro: a multi-residue deletion names its first and last residue with a range; a `del<size>` suffix
on a single anchor names only one endpoint and is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Lys23_Val25del` | recommended | self | the correct range form |
| `NP_003997.1:p.Lys23del6` | refused | — | a size suffix names only one endpoint — write `<start>_<end>del` |

## `deletion.md:91-94` — an exon label is not a residue range

> has never been allowed

Ferro: a description must name the amino acids affected; an `EX17del`-style exon label is not HGVS
and is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.EX17del` | refused | — | an exon label, not a residue range — never allowed |

See also → `substitution.md:17`, `duplication.md:5`, `frameshift.md:5`.
