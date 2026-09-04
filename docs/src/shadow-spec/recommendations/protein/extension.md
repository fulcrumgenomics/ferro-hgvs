# Extension — ferro's reading

ferro's reading of the HGVS protein **extension** recommendations, clause by clause. New here? See
[How to read a page](../../reading-guide.md). Examples run against `NP_003997.1` (dystrophin) in
the committed slice, whose stop codon sits at position 3686.

## `extension.md:5` — definition: the sequence extends past a terminus

> extending the reference amino acid sequence at the N- or C-terminal end with one or more amino acids

Ferro: an extension lengthens the protein at the N-terminus (a new upstream start, `ext-n`) or the
C-terminus (a lost stop, `extTer#`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Met1ext-5` | recommended | self | N-terminal: a new upstream start at `Met-5` |
| `NP_003997.1:p.Ter3686GlnextTer17` | recommended | self | C-terminal: the stop becomes `Gln`, a new stop 17 codons on |

## `extension.md:26-28` — a start-lost change is an insertion, not an extension

> this variant is **not** described as an extension

Ferro: when the change *alters* `Met1` itself (activating an upstream start), it is an insertion
(`insertion.md:5`) or deletion-insertion, not an extension — because part of the normal sequence
changes. Writing it as an `ext` of `Met1` is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Met1_Leu2insArgSerThrVal` | recommended | self | start-lost with an upstream start → insertion |
| `NP_003997.1:p.Met1Valext-4` | refused | — | `Met1` itself changes — not an extension (`extension.md:28`) |

## `extension.md:30-31` — the C-terminal no-stop extension, `*` is `Ter`

> a variant in the stop codon (`Ter` / `*`) at position 110, changing it to a `Gln`-codon

Ferro: a lost stop (no-stop variant) extends the C-terminus to the next stop; the `*` shorthand for
both the old and new stop is normalized to `Ter`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Ter3686GlnextTer17` | recommended | self | stop→`Gln`, new stop at `+17` |
| `NP_003997.1:p.*3686Glnext*17` | recommended | `NP_003997.1:p.Ter3686GlnextTer17` | the `*` stop shorthand normalized to `Ter` |

## `extension.md:33-34` — no new stop is `extTer?`

> adding a tail of new amino acids of unknown length (position `Ter?`)

Ferro: when the extended frame reaches no new stop, the length is unknown and the position is
`extTer?`; ferro accepts and preserves it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Ter3686ArgextTer?` | recommended | self | stop→`Arg`, no new stop encountered |

See also → `insertion.md:5`, `substitution.md:45-49`, `frameshift.md:5`.
