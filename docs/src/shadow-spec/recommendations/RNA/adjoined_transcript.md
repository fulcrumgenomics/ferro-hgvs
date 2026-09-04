# Adjoined transcript — ferro's reading

ferro's reading of `RNA/adjoined_transcript.md`. The rules are HGVS's; ferro's job is to produce the
form the recommendations prefer. Verdicts describe **ferro's output**:

- **recommended** — ferro's output is the form the recommendations prefer (whether the input was
  already that form, or ferro normalized it there).
- **conformant** — ferro's output is valid HGVS but not *yet* the recommended form — a ferro
  limitation or a deliberate maintainer house choice among conformant forms, with a tracking
  issue where one exists.
- **refused** — the input is not valid HGVS; ferro rejects it in strict mode (correct behavior).
- **bug** — ferro's output is not valid HGVS (a defect). None on this page — but note two `conformant`
  rows that are *not* mere limitations: ferro's output re-parses in both, which is what the verdict
  measures, while the parse behind it is wrong or an RFC 2119 MUST goes unenforced. Those are filed as
  [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212) (single-position breakpoints
  accepted) and [#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213) (linker silently
  absorbed into the 3' accession).

Each **Why** block is transcluded from the ruling ledger — the record's own one-line summary,
rendered here and linked to its full entry in
[NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md).
The reasoning lives once, in the ledger; it is never re-typed here.

**No ledger record cites any `RNA/adjoined_transcript.md` clause** — the whole document is a
record-level gap, so every section below carries **no Why block**. The reading is CONFIRM- (or
DISPUTE-) by-inspection against the spec text and the shipped code, not an adjudicated ruling. This
includes the document's two RFC 2119 clauses (`:20` REQUIRES, `:21` SHOULD) — the only two uppercase
keyword uses in all of `recommendations/` outside `style.md`.

An adjoined transcript is a two-partner fusion (`5'partner :: 3'partner`); ferro parses it to an
`RnaFusion` variant and normalization is a **pass-through** (`normalize/mod.rs:5388`), so the
round-trip is fully executable **without any reference** — an accepted fusion normalizes to itself,
and no partner accession has to sit in the committed slice for that to be observable. So the
`Normalizes to self` rows below are real round-trips, not parse-only stand-ins. (The example
accessions — `NM_152263.2`, `NM_002609.3`, `NM_002354.2`, `NM_000251.2` — are the spec's own; none
is in the slice, but nothing here reads their bases.)

## `adjoined_transcript.md:5` — definition

> Adjoined Transcript: a transcript (RNA molecule) composed of adjoined RNA from two or more contributing transcripts.

Ferro: an adjoined transcript is modelled as `RnaFusionVariant`, a 5' partner joined to a 3' partner
by `::`. The *definition* admits "two or more" contributing transcripts, but the *syntax* is
two-partner only — a deliberate, stated limitation (see `:17`), not a contradiction. `:14`
(`Adjoined transcripts are a product of some gene fusions.`) is descriptive framing.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | self | the canonical two-partner form (`TPM3::PDGFRB`); parsed to `RnaFusion`, preserved byte-identically through the pass-through normalizer |

## `adjoined_transcript.md:17` — two-partner adjoined transcripts only

> - This syntax is for two-partner adjoined transcripts only.

Ferro: exactly two partners are admitted; a third `::`-joined partner is rejected. Enforcement is
structural — the 3' interval parser rejects the trailing `::…` — not a named check.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | self | two partners — accepted; pinned by `f15_adjoined_transcript_two_partner` |
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924::NM_000251.2:r.212_*279` | refused | — | three partners — rejected ("Unexpected content after RNA interval"); pinned by `f15_three_partner_adjoined_invalid` |

## `adjoined_transcript.md:18` — RNA sequence only (no `c.` / `n.`)

> - This syntax is for RNA sequence only (no use of coding (`c.`) / non-coding DNA (`n.`) reference sequences).

Ferro: the fusion production is `r.`-only. The dispatcher routes to `parse_rna_fusion` only when
both halves carry `:r.`; a `c.`/`n.` `::` form falls through to the ordinary parser and dies on the
trailing `::`. `docs/NORMALIZATION_STAGE_AUDIT.md` (line 282) is correct that `:18` is a parser
constraint, not a normalization rule.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:c.-115_775::NM_002609.3:c.1580_*1924` | refused | — | `c.` axis — rejected on the trailing `::…`; pinned by `f17_adjoined_c_prefix_forbidden` |
| `NR_000001.1:n.1_100::NR_000002.1:n.200_300` | refused | — | `n.` axis — rejected the same way. The rejection message does not name the `:18` clause (minor gap) |

## `adjoined_transcript.md:19` — linker sequences use RNA character codes

> - Linker sequences are specified using [General Recommendations](../general.md) for RNA sequence character codes, e.g. `aggcucccuugg`

Ferro: the grammar element `[ linker_sequence "::" ]` (`syntax.yaml:336`) is **not modelled** —
`RnaFusionVariant` has no linker field. A linker-bearing input is accepted only **by accident**:
`parse_rna_fusion` splits at the *first* `::`, and the 3' half (linker + second `::` + partner) goes
to `parse_simple_accession`, which scans forward for the next `X.`-prefixed accession and swallows
every SAM-refname char — `:` included — before it. The linker is therefore absorbed into the 3'
partner accession, which becomes the custom accession `aggcucccuugg::NM_002609.3`. `Display` renders
custom accessions verbatim, so the string round-trips while the AST is wrong.

**This is silent corruption, filed as [#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213).**
The row below is `conformant` rather than `bug` only because the verdict measures ferro's *output*,
and that output re-parses byte-identically — the defect is in the AST, which the string cannot show.
What looks like alphabet enforcement — W3020 refusing `t`/`T` — is the whole-string thymine pre-scan
for any `r.` description (`error_handling/preprocessor.rs:1275`), not linker-aware: `::hello::`,
`::ACGU::`, an empty linker (`::::`) and a doubled linker (`::x::x::`) all pass and round-trip. The
alphabet's authority (`background/standards.md:45-61`, fifteen lowercase RNA symbols, via
`general.md:47`) is exactly what the existing `rna-axis-alignment-only-symbol-reach` ruling already
establishes for the `r.` axis, so a real linker validator has its authority ready-made — it simply
does not exist yet.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::aggcucccuugg::NM_002609.3:r.1580_*1924` | conformant | self | the `syntax.yaml` linker example. Round-trips byte-identically — **but the `::linker::` element is not modeled**: ferro silently absorbs it into the 3' partner accession (`three_prime.accession == "aggcucccuugg::NM_002609.3"`), so the parse is wrong and the byte-identical round-trip masks it. Silent AST corruption ([#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213)). Recommended behaviour: model `linker` as a first-class element and re-emit `::{linker}::` |

## `adjoined_transcript.md:20` — a range is REQUIRED, not a single position

> - This syntax REQUIRES the use of a range (not a single position) for `five_prime_range` / `three_prime_range`.

Ferro: `REQUIRES` is RFC 2119 MUST-strength (`style.md:9-12`) — a well-formed adjoined transcript has
a `_`-range on **both** sides. Ferro **accepts single-position breakpoints** on either or both sides
in every mode: `parse_rna_interval`'s last alternative builds a bare `RnaInterval::point`,
`parse_rna_fusion_breakpoint` performs no range check, and the pass-through normalizer emits the
input verbatim. Since ferro's output for a fusion is its input, an accepted single position becomes an
output that violates the `:20` REQUIRES — a ferro **enforcement gap** on the document's only
MUST-strength clause.

**Filed as [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212).** Worse than
unguarded: `test_parse_rna_fusion_single_positions` (`variant.rs:10042`) **pins the violation as
correct behaviour**. There is no recommended single-position spelling — a single position cannot
express the required range, so such a variant is simply not expressible as a conformant adjoined
transcript; the fix is to refuse it at parse in strict mode (new W-code naming `:20`). The rows
below are `conformant` in the verdict's sense — the output is well-formed HGVS that re-parses — not
`bug`; the recommended form uses `_`-ranges on both sides.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.775::NM_002609.3:r.1580` | conformant | self | single position on **both** sides — accepted and round-tripped; violates `adjoined_transcript.md:20` (RFC 2119 REQUIRES a range). Ferro should refuse in strict mode — [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212); pinned as accepted by `test_parse_rna_fusion_single_positions`. Recommended form uses `_`-ranges on both sides |
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580` | conformant | self | one-sided violation — a range on the 5' side, a single position on the 3' side; violates `adjoined_transcript.md:20` (RFC 2119 REQUIRES a range). Ferro should refuse in strict mode — [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212); recommended form uses `_`-ranges on both sides |
| `NM_152263.2:r.(-115_775)::NM_002609.3:r.1580_*1924` | conformant | self | `(a_b)` is **one uncertain position** (`uncertain.md`), not a range, so the 5' side violates `adjoined_transcript.md:20` (RFC 2119 REQUIRES a range); accepted and round-tripped. Ferro should refuse in strict mode — [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212); recommended form uses `_`-ranges on both sides |

## `adjoined_transcript.md:21-22` — `?` outer bounds when only the junction is analyzed

> - When the adjoined transcript junction but not the entire transcript is analyzed, the outer range bounds SHOULD be
>   specified with `?`, e.g. `NM_152263.2:r.?_775::NM_002609.3:r.1580_?`

Ferro: this SHOULD (`style.md:18-19`, recommended, ignorable with cause) is conditioned on an
**assay fact** — "the junction but not the entire transcript is analyzed" — which is invisible to a
normalizer (it keys on information the normalizer does not hold, the `RNA/delins.md:41` class). A
conformant tool must therefore (i) *accept* `?` outer bounds, (ii) *not rewrite* concrete bounds to
`?` or `?` to concrete (either direction asserts an assay fact), and (iii) not try to *enforce* the
SHOULD. Ferro does all three: `?_775` / `1580_?` is a range (unlike a bare single position, `:20`),
it is accepted, and the pass-through normalizer preserves it verbatim. Honoured to the full extent
evaluable.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.?_775::NM_002609.3:r.1580_?` | recommended | self | the spec's own `?`-outer-bound form — accepted and preserved verbatim (both `?` retained on both sides). No test pins this today (minor gap) |
| `NM_152263.2:r.-115_?::NM_002609.3:r.?_*1924` | recommended | self | `?` on the **inner** (junction) bounds — describes an unknown junction, which `:21` does not address; the grammar admits it and the spec is silent, so accept-and-preserve is the only conformant behaviour |

## `adjoined_transcript.md:23-25` — one format for all mechanisms

> - All adjoined transcripts are described using the same format, irrespective of whether they derive
>   from inter-chromosomal or intra-chromosomal DNA rearrangements (translocation, deletion, inversion)
>   or other mechanisms (trans-splicing).

Ferro: one format for every mechanism; the description must not encode the mechanism. `RnaFusionVariant`
has no mechanism field and the output is the input, so this is trivially conformant. Note the contrast
with DNA: HGVS *removed* `::` for DNA translocations (`DNA/complex.md`, delins only), so ferro's
`by_input` entries correctly reject DNA `::`; the RNA `::` here is a different production and is
rightly accepted — no tension.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | self | one format, whatever the underlying mechanism; ferro records no mechanism and preserves the form |

## `adjoined_transcript.md:29-41` — worked examples

> - **translocation-derived adjoined transcript**<br>
>     - **`NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924`**<br>
>         describes an adjoined transcript from a `TPM3::PDGFRB` gene fusion, where nucleotides `r.-115` to `r.775` (reference transcript
>         `NM_152263.2`, _TPM3_ gene) are coupled to nucleotides `r.1580` to `r.*1924` (reference transcript `NM_002609.3`, _PDGFRB_ gene).
>
> - **deletion-derived adjoined transcripts**
>     - **`NM_002354.2:r.-358_555::NM_000251.2:r.212_*279`**<br>
>         describes an adjoined transcript from an `EPCAM::MSH2` gene fusion, where nucleotides `r.-358` to `r.555` (reference transcript
>         `NM_002354.2`, _EPCAM_ gene) are coupled to nucleotides `r.212` to `r.*279` (reference transcript `NM_000251.2`, _MSH2_ gene).
>
>     - **`NM_002354.2:r.?_555::guaugauuuuuuaataa::NM_000251.2:r.212_?`**<br>
>         describes an adjoined transcript from an `EPCAM::MSH2` gene fusion, where only the fusion break point has been characterised,
>         showing the insertion of a 17 nucleotide sequence (`guaugauuuuuuaataa`) between two adjoined transcripts.

Ferro: the two range-only examples (ex 1, ex 2a) parse to `RnaFusion` with both accessions correct
and round-trip byte-identically — honest greens. The linker-bearing ex 2b is disputed **twice
over**, once by the spec and once by ferro.

**The spec's own example violates its own `:19`.** The linker `guaugauuuuuuaa**t**aa` contains `t`,
which is not an RNA character code (`background/standards.md:45-61`; `u` stands for `t`), contradicting `:19`
four lines above. Two clauses that cannot both hold — an open **upstream** inconsistency in
`assets/hgvs-nomenclature` that should be filed against the spec. The intended string is almost
certainly `guaugauuuuuuaauaa` (the 12-mer at `:19` and in `syntax.yaml` is spelled correctly); only
the alphabet is wrong, the length is 17 nt as stated.

**Ferro's parse is wrong in every mode.** Strict *refuses* ex 2b via W3020 (thymine in an `r.`
description) — defensible, but W3020 is a whole-string pre-scan, not linker validation. Lenient
rewrites `t→u` and *then* mis-parses the linker into the 3' accession, the `:19` / #2213 corruption.
So no mode produces a correct AST.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | self | example 1 (`TPM3::PDGFRB`) — both accessions correct (`test_parse_rna_fusion_basic`), preserved through normalize |
| `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279` | recommended | self | example 2a (`EPCAM::MSH2`) — range-only, parsed and preserved correctly |
| `NM_002354.2:r.?_555::guaugauuuuuuaataa::NM_000251.2:r.212_?` | refused | — | example 2b, as the spec spells it. The `t` in the linker violates the page's own `:19` (upstream spec inconsistency — file against `assets/hgvs-nomenclature`; corrected form `guaugauuuuuuaauaa`). Ferro strict refuses via the W3020 thymine pre-scan; lenient rewrites `t→u` then absorbs the linker into the 3' accession ([#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213)) — no mode produces a correct AST |

See also → `adjoined_transcript.md:19` (the linker-corruption defect, #2213), `adjoined_transcript.md:20`
(single-position breakpoints accepted, #2212).
