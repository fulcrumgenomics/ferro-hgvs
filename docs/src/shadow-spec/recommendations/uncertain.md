# Uncertain — ferro's reading

ferro's reading of the HGVS **uncertain** recommendations — how a variant is described when not all
details are known. New here? See [How to read a page](../reading-guide.md). Executable examples run
against `NM_004006.3` / `NP_003997.1` in the committed slice; genomic-breakpoint examples use foreign
accessions and are parse-only.

## `uncertain.md:9-16` — the uncertainty characters

> `( )` (parentheses) are used to indicate uncertainties

Ferro: `( )` marks a predicted or uncertain form, `^` "or" (a set of possible residues), and
`Xaa[n]` / `N[n]` a run of unknown residues / bases; `?` marks an unknown position (shown in the
breakpoint examples below). ferro parses and preserves all of them.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.(Trp24Cys)` | recommended | self | `( )` — a predicted consequence |
| `NM_004006.3:c.(71_72)G>A` | recommended | self | `( )` — an uncertain position |
| `NP_003997.1:p.(Gly56Ala^Ser^Cys)` | recommended | self | `^` — one of `Ala`, `Ser`, `Cys` |
| `NP_003997.1:p.Glu719(Ala^Ser)fsTer23` | recommended | self | `^` inside a frameshift — the first new residue is `Ala` or `Ser` |

## `uncertain.md:181-182` — RNA when RNA was not analysed

> an effect on the RNA level is expected, but that it is not possible to give a reliable prediction of the consequences

Ferro: the RNA "not analysed" forms — `r.?` (effect expected, unpredictable), `r.(=)` (no change
expected), `r.spl` / `r.spl?` (splicing very likely / might be affected), `r.(76a>c)` (expected but
unverified), `r.(?)` (no change beyond the DNA-level consequence) — are all parsed and preserved.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.?` | recommended | self | an RNA effect is expected but not predictable |
| `NM_004006.3:r.spl` | recommended | self | splicing very likely affected (RNA not analysed) |
| `NM_004006.3:r.(=)` | recommended | self | no RNA change expected |
| `NM_004006.3:r.(76a>c)` | recommended | self | expected but RNA not analysed |

The spec also lists `r.0?` (possibly no transcript); ferro does not yet parse the bare `r.0?` form —
though the protein counterpart `p.0?` does parse (an axis asymmetry).

## `uncertain.md:206-207` — protein when the consequence is unknown

> an effect on the protein level is expected

Ferro: `p.?` (effect expected, unpredictable), `p.0?` (possibly no protein) and `p.(=)` (no change
expected) are all parsed and preserved.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.?` | recommended | self | a protein effect is expected but not predictable |
| `NP_003997.1:p.0?` | recommended | self | possibly no protein produced |
| `NP_003997.1:p.(=)` | recommended | self | no protein change expected |

## `uncertain.md:215-222` — an unknown position between two residues

> at an unknown position between amino acids `Ala123` and `Pro131`

Ferro: a change whose exact *position* is unknown is anchored to a residue range `(start_end)` — for
a nonsense (`)Ter`), a frameshift (`)fs`), or an insertion (`)insXaa[n]`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.(Val123_Leu131)Ter` | recommended | self | a stop at an unknown position in `Val123`–`Leu131` |
| `NP_003997.1:p.(Val123_Leu131)fs` | recommended | self | a frameshift starting at an unknown position |
| `NP_003997.1:p.(Val123_Leu131)insXaa[4]` | recommended | self | four unknown residues inserted at an unknown position |

## `uncertain.md:32-33` — a deletion whose breakpoints were not sequenced

> `(A_B)_(C_D)del`, where `B_C` describes the **minimal** extent and `A_D` the **maximal** extent

Ferro: an unsequenced breakpoint is written as nested `(min)_(max)` position ranges, with `?` for an
open end. These examples use genomic coordinates on a foreign build, so they are parse-only here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.(31729716_31774235)_(32216847_32287541)del` | recommended | — | nested min/max breakpoints (parse-only here) |
| `NC_000023.10:g.(?_32238146)_(32984039_?)del` | recommended | — | open-ended breakpoints with `?` (parse-only here) |

See also → [predicted substitutions](DNA/substitution.md), [protein frameshift](protein/frameshift.md), [general notation](general.md).
