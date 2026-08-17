# Substitution — ferro's reading

ferro's reading of `substitution.md`. The rules are HGVS's; ferro's job is to produce the form the
recommendations prefer. Verdicts describe **ferro's output**:

- **recommended** — ferro's output is the form the recommendations prefer (whether the input was
  already that form, or ferro normalized it there).
- **conformant** — ferro's output is valid HGVS but not *yet* the recommended form — a ferro
  limitation, with a tracking issue.
- **refused** — the input is not valid HGVS; ferro rejects it in strict mode (correct behavior).
- **bug** — ferro's output is not valid HGVS (a defect). None on this page.

Each **Why** block is transcluded from the ruling ledger — the record's own one-line summary,
rendered here and linked to its full entry in
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md).
The reasoning lives once, in the ledger; it is never re-typed here.

## `substitution.md:5` — definition: one nucleotide for one

> **one** nucleotide is replaced by **one** other nucleotide

Ferro: a substitution is exactly 1→1; IUPAC ambiguity codes stand in for the single replacing base.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76A>G` | recommended | self | canonical single-base substitution |
| `NM_004006.3:c.54G>H` | recommended | self | replacement is IUPAC `H` = `A`, `C`, or `T` |
| `NC_000023.10:g.33038255C>A` | recommended | — | the spec's genomic example (GRCh37 accession — parse-only here) |
| `NM_004006.3:c.79GC>TT` | refused | — | two nt — violates 1→1; rejected as "deprecated multi-base substitution syntax" |

## `substitution.md:15` — two or more nucleotides is a delins

> two or more consecutive nucleotides are described as deletion/insertion (delins) variants

Ferro: a multi-base replacement is never a substitution; write it as a delins.

**Why.**
<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.79_80delinsTT` | recommended | self | the two-base change written as a delins |
| `NM_004006.3:c.79_80GC>TT` | refused | — | rejected: "deprecated multi-base substitution syntax" — use delins |

## `substitution.md:16-18` — separation, and the one-codon exception

> two variants separated by one or more nucleotides should be described individually

Ferro: the separation rule (ruleset rule 2); the exception folds two subs into one delins only when
they sit one nucleotide apart and together change one amino acid.

**Why.**
<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[codon-carve-out-shape-restriction](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together alter a single amino acid are written as one delins, whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled.
<!-- why:END:separation-rule-force-modal-or-negation,codon-carve-out-shape-restriction -->

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.145_147delinsTGG` | recommended | self | one-codon delins, `CGC`→`TGG` |
| `NM_004006.3:c.[145C>T;147C>G]` | recommended | `NM_004006.3:c.145_147delinsTGG` | ferro normalizes the split to the recommended delins (warns `MEMBERS_COALESCED_FROM_REPORTED_FORM`) |
| `NM_004006.3:c.235_237delinsTAT` | recommended | self | the Lys79 codon; the split predicts `p.[Lys79*;Lys79Asn]`, the delins `p.Lys79Tyr` |
| `NM_004006.3:c.[235A>T;237G>T]` | recommended | `NM_004006.3:c.235_237delinsTAT` | split → the recommended delins (same codon exception) |

**Confluence.** Each split and its delins converge on one recommended output —
`{c.[145C>T;147C>G], c.145_147delinsTGG} → c.145_147delinsTGG`, and the Lys79 pair →
`c.235_237delinsTAT`.

## `substitution.md:19` — no change is `=`, not a substitution

> found **not changed** are described as

Ferro: a tested, unchanged position is `=`; a no-change substitution is normalized to it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123=` | recommended | self | `c.123` screened, unchanged (reference `C`) |
| `NM_004006.3:c.123C>C` | recommended | `NM_004006.3:c.123=` | a no-change substitution is normalized to the recommended `=` |

## `substitution.md:20` — polymorphisms are not `A/G`

> it is not correct to describe

Ferro: the slash form for "polymorphism" is not valid HGVS; ferro rejects it.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.76A>G` | recommended | self | the correct form |
| `NM_004006.3:c.76A/G` | refused | — | the disallowed "polymorphism" slash — rejected at parse |

## `substitution.md:32` — adjacent substitutions are one delins

> changes involving two or more consecutive nucleotides are described as deletion-insertion (delins) so the description

Ferro: at separation zero the split is not correct; two adjacent members that both consume reference
bases are one delins. The same rule licenses re-merging adjacency that per-member 3′ shifting
creates (protein `p.[Gly16Ala;Gly17del]` → `p.Gly16_Gly17delinsAla`), keeping normalization
idempotent.

**Why.**
<!-- why:START -->
> **[delins-adjacent-members-when-both-consume-reference](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero.
>
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:delins-adjacent-members-when-both-consume-reference,canonical-form-choice-when-both-legal -->

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.79_80delinsTT` | recommended | self | adjacent `c.79`,`c.80` written as one delins |
| `NM_004006.3:c.[79G>T;80C>T]` | recommended | `NM_004006.3:c.79_80delinsTT` | ferro normalizes the adjacent split to the recommended delins |
| `NM_004006.3:c.79_80GC>TT` | refused | — | the multi-base substitution spelling — rejected |

**Confluence.** `{c.[79G>T;80C>T], c.79_80delinsTT} → c.79_80delinsTT`.

See also → `substitution.md:15`, `substitution.md:16-18`.

## `substitution.md:47-49` — mosaic and chimeric

> a mosaic case where at position

Ferro: mosaic (`/`) and chimeric (`//`) mixtures are valid; the recommendations write the reference
allele first.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.85=/T>C` | recommended | self | mosaic: reference `=` written first, then `T>C` |
| `NM_004006.3:c.85=//T>C` | recommended | self | chimeric: a mix of `c.85=` and `c.85T>C` cells |
| `NM_004006.3:c.85T>C/=` | conformant | self | valid, but `substitution.md:49` writes the reference first; ferro does not yet reorder this to the recommended `c.85=/T>C`. Tracked by [#2034](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2034). |
