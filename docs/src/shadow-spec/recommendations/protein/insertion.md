# Insertion — ferro's reading

ferro's reading of the HGVS protein **insertion** recommendations, clause by clause. New here? See
[How to read a page](../../reading-guide.md). Examples run against `NP_003997.1` (dystrophin) in
the committed slice.

## `insertion.md:5` — definition: residues added between two others

> one or more amino acids are inserted, which is not a frameshift

Ferro: an insertion adds residues between two flanking positions, and is neither a frameshift nor a
tandem copy of the N-terminal neighbour.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Lys23_Trp24insAla` | recommended | self | `Ala` inserted between `Lys23` and `Trp24` |
| `NP_003997.1:p.Leu2_Trp3insGlnSerLys` | recommended | self | a three-residue insertion |

## `insertion.md:17-18` — two flanking positions, never one

> an insertion can not be described using **one** amino acid position

Ferro: the anchor is two *flanking* residues; a single-position insertion is ambiguous (at, or
after?) and is refused. The parser is deliberately lenient about *adjacency*, though — it accepts a
non-adjacent pair such as `p.Lys23_Val25insAla`, which round-trips, even though the spec asks for two
flanking residues ([#1264](https://github.com/fulcrumgenomics/ferro-hgvs/issues/1264)).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Lys23_Trp24insAla` | recommended | self | the correct two-residue anchor |
| `NP_003997.1:p.Lys23insAsp` | refused | — | a single-position anchor — refused (`insertion.md:18`) |

## `insertion.md:20` — a duplicating insertion is a duplication

> duplicating insertions should be described as duplications

Ferro: an insertion whose payload is a tandem copy of its N-terminal neighbour is rewritten to the
`dup` form — the description follows from the resulting sequence, not the input's spelling.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Trp4_Glu5insTrp` | recommended | `NP_003997.1:p.Trp4dup` | tandem-copy insertion → `dup` (`duplication.md:19`) |

## `insertion.md:21` — a large insertion may be sized

> the insertion may be described by its length

Ferro: a long insertion may be written by length (`insXaa[n]`) or by a stop position (`insTer#`);
ferro accepts and preserves both.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.Lys23_Trp24insXaa[23]` | recommended | self | a 23-residue in-frame insertion, given by length |
| `NP_003997.1:p.Lys23_Trp24insTer12` | recommended | self | an insertion sized by a stop at position 12 |
| `NP_003997.1:p.(Ser332_Ser333insXaa)` | recommended | self | a single unknown inserted residue (predicted) |
| `NP_003997.1:p.(Val582_Asn583insXaa[5])` | recommended | self | five unknown inserted residues (predicted) |

## `insertion.md:70-74` — the caret form is not allowed

> insertions can not be described using the format

Ferro: the `^` "between" caret is not valid HGVS; ferro rejects it at parse.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.123^124Ala` | refused | — | the disallowed `^` caret form — rejected at parse |

## `insertion.md:76-80` — a copy from elsewhere is still an insertion

> The variant should be described as an insertion; `p.His7_Gln8insGly4_Ser6`

Ferro: when the inserted residues are a copy of a sequence elsewhere in the protein (not a tandem
neighbour), the spec writes an insertion that *names the source range* — `insGly4_Ser6`. Ferro does
**not** yet accept this `ins<range>` back-reference on the `p.` axis: inserted residues must be given
literally, so this form is rejected at parse. This is a ferro limitation, not an invalid input — the
spelling is valid HGVS. Tracked by [#2219](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2219).

See also → `duplication.md:19`, `delins.md:5`, `frameshift.md:5`.
