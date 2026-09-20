# Adjoined transcript — ferro's reading

ferro's reading of the HGVS **adjoined transcript** recommendations on the transcript (`r.`)
axis, clause by clause, with each spelling's normalized form and a verdict. See
[How to read a page](../../reading-guide.md) for verdicts, table conventions, and recurring terms.

Four `conformant` rows on this page are not mere limitations: ferro's output re-parses in each,
which is what the verdict measures, while the parse behind it is wrong or a requirement the spec
states plainly goes unenforced. Three are single-position breakpoints ferro accepts
([#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212)); one is the linker silently
absorbed into the 3' accession ([#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213)).
The worked linker example that ferro strict-refuses is separate, and is not one of these four.

**No ledger record cites any `RNA/adjoined_transcript.md` clause** — the whole document is a
record-level gap, so every section below carries **no Why block**. The reading is
CONFIRM-by-inspection (or DISPUTE-by-inspection) against the spec text and the shipped code, not an
adjudicated ruling. This
includes the document's two uppercase-keyword clauses (`:20` REQUIRES, `:21` SHOULD) — the only
uppercase keyword uses in the spec's own quoted `recommendations/` clauses outside `style.md` (a
scope that excludes these shadow pages' own prose and quoted error messages, which use `MUST`
freely). Note only `:21`'s **SHOULD**
is an actual RFC 2119 keyword; `:20`'s **REQUIRES** is not in RFC 2119 (which lists REQUIRED), so it
reads as emphatic spec prose, not a formal keyword — and per `READING_THE_SPEC.md`, clauses are not
ranked by keyword strength regardless.

An adjoined transcript is a two-partner fusion (`5'partner :: 3'partner`); ferro parses it to an
`RnaFusion` variant and normalization is a **pass-through**, so the round-trip is fully executable
**without any reference** — an accepted fusion normalizes to itself, and no partner accession has
to sit in the committed slice for that to be observable. So the `Normalizes to self` rows below are
real round-trips, not parse-only stand-ins. (The example accessions are the spec's own; none is in
the slice, but nothing here reads their bases.)

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
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | self | two partners — accepted |
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924::NM_000251.2:r.212_*279` | refused | — | three partners — rejected ("Unexpected content after RNA interval") |

## `adjoined_transcript.md:18` — RNA sequence only (no `c.` / `n.`)

> - This syntax is for RNA sequence only (no use of coding (`c.`) / non-coding DNA (`n.`) reference sequences).

Ferro: the fusion production is `r.`-only. The dispatcher routes to the fusion parser only when
both halves carry `:r.`; a `c.`/`n.` `::` form falls through to the ordinary parser and dies on the
trailing `::`. This is a parser constraint, not a normalization rule.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:c.-115_775::NM_002609.3:c.1580_*1924` | refused | — | `c.` axis — rejected on the trailing `::…` |
| `NR_000001.1:n.1_100::NR_000002.1:n.200_300` | refused | — | `n.` axis — rejected the same way. The rejection message does not name the `:18` clause (minor gap) |

## `adjoined_transcript.md:19` — linker sequences use RNA character codes

> - Linker sequences are specified using [General Recommendations](../general.md) for RNA sequence character codes, e.g. `aggcucccuugg`

Ferro: the grammar's `linker_sequence` element between the two `::` is **not modelled** —
`RnaFusionVariant` has no linker field. A linker-bearing input is accepted only **by accident**:
the parser splits at the *first* `::`, and the 3' half (linker + second `::` + partner) is scanned
forward for the next accession, swallowing every character — `:` included — before it. The linker
is therefore absorbed into the 3' partner accession, becoming the custom accession
`aggcucccuugg::NM_002609.3`. Custom accessions render verbatim, so the string round-trips while the
AST is wrong.

**This is silent corruption, filed as [#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213).**
The row below is `conformant` rather than `bug` only because the verdict measures ferro's *output*,
which re-parses byte-identically — the defect is in the AST, which the string cannot show. What
looks like alphabet enforcement — W3020 refusing `t`/`T` — is a whole-string thymine pre-scan for
any `r.` description, not linker-aware: `::hello::`, `::ACGU::`, an empty linker (`::::`) and a
doubled linker (`::x::x::`) all pass and round-trip. The RNA alphabet's authority is already
established for the `r.` axis by `rna-axis-alignment-only-symbol-reach`, so a real linker validator
has its authority ready-made — it does not exist yet.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::aggcucccuugg::NM_002609.3:r.1580_*1924` | conformant | self | the grammar's linker example. Round-trips byte-identically — **but the `::linker::` element is not modeled**: ferro silently absorbs it into the 3' partner accession (`aggcucccuugg::NM_002609.3`), so the parse is wrong and the byte-identical round-trip masks it. Silent AST corruption ([#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213)). Recommended behaviour: model `linker` as a first-class element and re-emit `::{linker}::` |

## `adjoined_transcript.md:20` — a range is REQUIRED, not a single position

> - This syntax REQUIRES the use of a range (not a single position) for `five_prime_range` / `three_prime_range`.

Ferro: `:20` states this plainly — "REQUIRES the use of a range (not a single position)". (`REQUIRES`
is not itself an RFC 2119 keyword — RFC 2119 lists REQUIRED — so this reads as emphatic spec prose,
not a formal MUST, though the wording here is unambiguous.) A well-formed adjoined transcript has a
`_`-range on **both** sides. Ferro **accepts single-position breakpoints** on either or both sides in
every mode, with no range check, and the pass-through normalizer emits the input verbatim. Since
ferro's output for a fusion is its input, an accepted single position becomes an output that
violates `:20`'s stated requirement — a ferro **enforcement gap** on the one clause in this document
that reads as a hard requirement.

**Filed as [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212).** Worse than
unguarded: an existing test pins the violation as correct behaviour. There is no recommended
single-position spelling — a single position cannot express the required range, so such a variant is
not expressible as a conformant adjoined transcript; the fix is to refuse it at parse in
strict mode. The rows below are `conformant` in the verdict's sense — the output is well-formed
HGVS that re-parses — not `bug`; the recommended form uses `_`-ranges on both sides.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.775::NM_002609.3:r.1580` | conformant | self | single position on **both** sides — accepted and round-tripped; violates `:20`'s range requirement. Ferro should refuse in strict mode — [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212) |
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580` | conformant | self | one-sided violation — a range on the 5' side, a single position on the 3' side; same `:20` violation and [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212) |
| `NM_152263.2:r.(-115_775)::NM_002609.3:r.1580_*1924` | conformant | self | `(a_b)` is **one uncertain position** (`uncertain.md`), not a range, so the 5' side violates `:20`; accepted and round-tripped. Same [#2212](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2212) |

## `adjoined_transcript.md:21-22` — `?` outer bounds when only the junction is analyzed

> - When the adjoined transcript junction but not the entire transcript is analyzed, the outer range bounds SHOULD be
>   specified with `?`, e.g. `NM_152263.2:r.?_775::NM_002609.3:r.1580_?`

Ferro: this SHOULD is conditioned on an **assay fact** — "the junction but not the entire transcript
is analyzed" — invisible to a normalizer. A conformant tool must therefore (i) *accept* `?` outer
bounds, (ii) *not rewrite* concrete bounds to `?` or vice versa (either direction asserts an assay
fact), and (iii) not try to *enforce* the SHOULD. Ferro does all three: `?_775` / `1580_?` is a
range (unlike a bare single position, `:20`), it is accepted, and the pass-through normalizer
preserves it verbatim.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.?_775::NM_002609.3:r.1580_?` | recommended | self | the spec's own `?`-outer-bound form — accepted and preserved verbatim (both `?` retained on both sides). No test pins this today (minor gap) |
| `NM_152263.2:r.-115_?::NM_002609.3:r.?_*1924` | conformant | self | `?` on the **inner** (junction) bounds — ferro accepts and preserves this unspecified form; `:21` addresses outer bounds only |

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
which is not an RNA character code (`u` stands for `t`), contradicting `:19` four lines above — an
open **upstream** inconsistency that should be filed against the spec. The intended string is
almost certainly `guaugauuuuuuaauaa`; only the alphabet is wrong.

**Ferro's parse is wrong in every mode.** Strict *refuses* ex 2b via W3020 (thymine in an `r.`
description) — defensible, but W3020 is a whole-string pre-scan, not linker validation. Lenient
rewrites `t→u` and *then* mis-parses the linker into the 3' accession, the same #2213 corruption.
No mode produces a correct AST.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_152263.2:r.-115_775::NM_002609.3:r.1580_*1924` | recommended | self | example 1 (`TPM3::PDGFRB`) — both accessions correct (`test_parse_rna_fusion_basic`), preserved through normalize |
| `NM_002354.2:r.-358_555::NM_000251.2:r.212_*279` | recommended | self | example 2a (`EPCAM::MSH2`) — range-only, parsed and preserved correctly |
| `NM_002354.2:r.?_555::guaugauuuuuuaataa::NM_000251.2:r.212_?` | refused | — | example 2b, as the spec spells it. The `t` in the linker violates the page's own `:19` (upstream spec inconsistency — file against `assets/hgvs-nomenclature`; corrected form `guaugauuuuuuaauaa`). Ferro strict refuses via the W3020 thymine pre-scan; lenient rewrites `t→u` then absorbs the linker into the 3' accession ([#2213](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2213)) — no mode produces a correct AST |

See also → `adjoined_transcript.md:19` (the linker-corruption defect, #2213), `adjoined_transcript.md:20`
(single-position breakpoints accepted, #2212).
