# Duplication — ferro's reading

ferro's reading of the HGVS protein **duplication** recommendations, clause by clause. New here?
See [How to read a page](../../reading-guide.md). Examples run against `NP_003997.1` (dystrophin)
in the committed slice.

## `duplication.md:5` — definition: a tandem copy, directly C-terminal

> a copy of one or more amino acids is inserted **directly C-terminal** of the original copy

Ferro: a duplication is a copy placed immediately C-terminal of the original — a tandem copy.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Val7dup` | recommended | self | a single-residue tandem duplication |
| `NP_003997.1:p.(Val7dup)` | recommended | self | the predicted-consequence form |
| `NP_003997.1:p.Lys23_Val25dup` | recommended | self | a three-residue tandem duplication |

## `duplication.md:19-20` — tandem only; otherwise an insertion

> duplication may only be used when the additional copy is **directly C-terminal**

Ferro: `dup` is only for a copy directly C-terminal of the original. A copy that is *not* in tandem
is an insertion (`insertion.md:20`) — and, conversely, an insertion whose payload is a tandem copy
of its N-terminal neighbour is rewritten to `dup`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp4_Glu5insTrp` | recommended | `NP_003997.1:p.Trp4dup` | a tandem-copy insertion is normalized to the `dup` form |

## `duplication.md:22` — the 3′ rule

> the **most C-terminal position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**)

Ferro: in a run of identical residues the duplication is placed at the **most C-terminal** position.
In `MetLeuTrpTrpGlu`, duplicating either `Trp` gives the same protein, so ferro writes the 3′-most.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp3dup` | recommended | `NP_003997.1:p.Trp4dup` | 3′ rule: shifted to the C-terminal `Trp4` |
| `NP_003997.1:p.Trp4dup` | recommended | self | already the 3′-most residue of the `TrpTrp` run |

## `duplication.md:71-74` — mosaic duplication

> a mosaic case where at amino acid position 7

Ferro: a mosaic duplication uses `=/dup`, reference first; ferro preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Val7=/dup` | recommended | self | mosaic: reference `Val7=` first, then the duplication |

See also → `insertion.md:20`, `deletion.md:19`, `extension.md:5`.
