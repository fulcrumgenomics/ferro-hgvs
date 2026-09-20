# Checklist — ferro's reading

ferro's reading of the HGVS **checklist** — the most frequently offended rules, as ferro
enforces them. New here? See [How to read a page](../reading-guide.md). Examples run against
`NM_004006.3` / `NP_003997.1` in the committed slice. The spec's checklist has 14 numbered items;
the sections below shadow the ones with a mechanically-checkable form — the rest are
publication-practice advice.

## `checklist.md:19-20` — an intronic `c.` position needs a genomic reference

> can only be used to describe variants in introns using a `c.` prefix when a genomic reference sequence is given

Ferro: an intronic `c.` position (`c.94-2A>G`) names an intron base a bare `NM_` transcript does not
contain. It is valid HGVS only against a genomic-anchored reference — so ferro refuses it on a bare
transcript in **strict** mode (`W4007`) and accepts it in **lenient** mode. Because the input is
valid, not malformed, this is a house input-hygiene policy, not a `refused`-as-invalid verdict.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[bare-transcript-intronic-position](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A bare transcript intronic position such as 'NM_...:c.20+2del' is refused in strict input-hygiene mode and accepted in lenient mode, since the spec allows it only when a genomic reference is given.
<!-- why:END:bare-transcript-intronic-position -->

</details>

The genomic-anchored form `NC_000023.11(NM_004006.3):c.94-2A>G` carries the intron bases and is
accepted in both modes.

## `checklist.md:24-26` — intronic descriptions: no labels, no incomplete ranges

> descriptions referring to exon or intron numbers instead of nucleotide positions

Ferro: an `IVS`/exon-number label is not a position and is refused; an intronic range must name both
endpoints in full (`c.123-65_123-50`), not abbreviate the second (`c.123-65_-50`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.IVS4-2A>G` | refused | — | an `IVS` intron label, not a nucleotide position |
| `NM_004006.3:c.123-65_-50` | refused | — | an incomplete intronic range — the second endpoint is abbreviated |

## `checklist.md:30-33` — insertions: two positions, and give the sequence

> insertions should be reported using the format `c.51_52insT`

Ferro: an insertion names two flanking positions and the inserted sequence; a single-position anchor
or a bare count is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.51_52insT` | recommended | self | the correct form: two positions, inserted base given |
| `NM_004006.3:c.52insT` | refused | — | a single-position anchor — ambiguous |
| `NM_004006.3:c.5439_5440ins6` | refused | — | a bare count — the inserted sequence must be given (or `insN[6]`) |

## `checklist.md:42-45` — the range sign is `_`, not `-`

> The sign used to indicate a range is the `_` (underscore), not a `-` (hyphen-minus)

Ferro: `_` joins a range; a hyphen is a coding-DNA intronic offset. `c.12-14del` is **not** "12 to
14" — it is a deletion at the intronic position `c.12-14` (14 nt into the intron 5′ of `c.12`), a
*valid but entirely different* description. That silent change of meaning is exactly why the hyphen
is a trap; on a bare transcript ferro refuses that intronic form in strict mode (see above).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.12_14del` | recommended | self | a deletion of `c.12` to `c.14` — the range you meant |

## `checklist.md:47-49` — a deletion names its first and last residue

> correct is `g.123_125del`

Ferro: a multi-residue deletion names both endpoints with a range; a `del<size>` suffix is refused.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_125del` | recommended | self | the correct range form |
| `NM_004006.3:c.123del3` | refused | — | a size suffix names only one endpoint |

## `checklist.md:60-65` — protein: `Ter` not `X`, silent is `=`, no `Met1` substitution

> predicted "**silent**" protein level variants are described as **p.(Leu54=)**

Ferro: a stop codon is `Ter` (or `*`), never `X`; a silent change is `=`; and an initiation-codon
change is not a substitution.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NP_003997.1:p.(Leu54=)` | recommended | self | the recommended silent form |
| `NP_003997.1:p.Trp24Ter` | recommended | self | a stop written `Ter` |
| `NP_003997.1:p.Leu54Leu` | conformant | self | silent (`Leu`→`Leu`); the recommended form is `p.Leu54=` — ferro preserves it ([#2220](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2220)) |
| `NP_003997.1:p.(Met1Val)` | refused | — | an initiation-codon change is not a substitution |

## `checklist.md:67-69` — no "polymorphism" slash

> should not be described using the `/` (slash)

Ferro: the `/` "polymorphism" form is not valid HGVS; describe the change as an ordinary variant.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.127G>A` | recommended | self | the ordinary substitution form |
| `NM_004006.3:c.127G/A` | refused | — | the disallowed "polymorphism" slash |

See also → [general notation](general.md), [DNA insertion](DNA/insertion.md), [protein substitution](protein/substitution.md).
