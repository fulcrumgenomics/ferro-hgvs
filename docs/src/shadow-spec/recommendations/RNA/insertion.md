# Insertion — ferro's reading

ferro's reading of the HGVS **insertion** recommendations on the transcript (`r.`) axis, clause
by clause — each spelling with the form ferro normalizes it to and a verdict on that output.
New here? See [How to read a page](../../reading-guide.md) for the verdicts, the table
conventions, and the recurring terms.

*DNA twin: [Insertion (`c.`/`g.`)](../DNA/insertion.md).*

Five ledger records reach this page's clauses.
`inverted-duplication-is-derived-as-ins-range-inv` cites `RNA/insertion.md:18` directly, so it
carries `r.` authority on its own molecule-native basis. `duplication-must-ranks-the-label-not-the-partition`
cites the DNA twin `DNA/duplication.md:18`, so on the `r.` axis the MUST that separates a `dup`
from an `ins` here rests on the verbatim RNA twin `RNA/duplication.md:19` — noted where it applies
rather than claimed as an `r.`-jurisdiction ruling. `canonical-form-choice-when-both-legal` and
`absolute-prohibition-enforcement-stage` are molecule-neutral (a general mechanism, a general
mode-policy) and reach the `r.` axis unchanged. `spdi-n-unit-repeat-refusal` is a house choice on
the HGVS→SPDI conversion path — not a conformance ruling and not axis-scoped — and reaches the
`n`-repeat examples at `:38-42`. The definition, the flanking/order Notes, the 3'rule Note and the
literal-base Examples carry **no Why block** — the reading is CONFIRM-by-inspection against the
spec text and the shipped code, not an adjudicated ruling.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
examples are written on `LRG_199t1`, `NM_004006.2` and on an `NG_012232.1(NM_004006.2)` splicing
form, none of which the slice carries, so those rows are parse-only (`—`) — ferro cannot read their
bases here. The `NM_004006.3` base facts the executable rows rely on (`r.123` is `c`, `r.124` is
`a`; `r.6_7` is `uu`; `r.5_16` is `uuuggugggaag`; `r.123_125` and `r.457_459` both read `cag`, with
`r.460` = `g`) are the same ones established on `RNA/deletion.md` and `RNA/duplication.md`, plus
the `r.457_460` facts measured for the tandem-copy rows at `:17`.

## `insertion.md:5` — definition: inserted, and not a copy of the 5' flank

> Insertion: a sequence change where, compared to the reference sequence, one or more nucleotides are inserted **and** where the insertion is not a copy of a sequence immediately 5'.

Ferro: an insertion adds one or more nucleotides between two flanking positions, and the "not a
copy of a sequence immediately 5'" clause is exactly the `dup`/`ins` boundary — a copy directly
3'-flanking of its original is a duplication (see `:17` and `RNA/duplication.md:18-21`), so ferro
relabels a duplicating insertion to `dup` rather than leaving it as an `ins`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_124insg` | recommended | self | a genuine single-nucleotide insertion between `r.123` (`c`) and `r.124` (`a`); the inserted `g` is not a copy of either flank, so it is neither a `dup` nor 3'-shiftable — a fixed point |
| `NM_004006.3:r.6_7insu` | recommended | `NM_004006.3:r.7dup` | the boundary case: the inserted `u` **is** a copy of the `u` at `r.6` immediately 5', so by the definition it is not an insertion — ferro relabels it to the recommended `dup` (`rules::insertion_is_duplication`, reached by `normalize_rna`), landing on the 3'-most `u` (`r.7`). Mirrors `issue_736::rna_insertion_with_u_collapses_to_dup` |

## `insertion.md:15` — the flanking positions must be two *adjacent* nucleotides

> - `positions flanking` should contain **two flanking nucleotides**, e.g., `123_124`, not `123_125`.

Ferro: an insertion occupies the zero-width junction between two adjacent positions, so the anchor
must name a pair that is exactly one apart. `validate_no_point_insertion` refuses a non-adjacent
anchor with a structured `InvalidEdit` error citing the DNA twin of this Note — there is no way to
recover which flanking pair the author meant, so this is a rejection, not a repair.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_124insg` | recommended | self | correctly flanking (adjacent) pair |
| `NM_004006.3:r.123_125insg` | refused | — | non-adjacent flanks — rejected at parse ("insertion anchor positions are not flanking … the two positions MUST be adjacent"); the axis-general check in `validate_no_point_insertion` handles `r.` identically to its `g.`/`c.` twins |

## `insertion.md:16` — the flanking positions are listed 5' to 3'

> - the `positions_flanking` should be listed from **5' to 3'**, e.g., `123_124`, not `124_123`.

Ferro: a reversed, non-circular flank order has no exception on the `r.` axis (the `<high>_<low>`
form is admitted only for `o.`/`m.` circular references), so it is refused rather than reordered —
the twin of `RNA/deletion.md:17` and `RNA/duplication.md:17`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_124insg` | recommended | self | correctly ordered flank |
| `NM_004006.3:r.124_123insg` | refused | — | reversed flank order on a non-circular reference — rejected at parse |

## `insertion.md:17` — a tandem duplication is a `dup`, not an insertion

> - tandem duplications are described as a duplication (`r.123_456`<code class="spot1">dup</code>), not an insertion (<code class="invalid">r.456_457ins123_456</code>, see [Prioritization](../general.md)).

Ferro: the `dup`-over-`ins` MUST seen from the insertion side — a copy directly 3'-flanking of its
original is spelled `dup`, so a duplicating insertion is relabeled rather than kept as an `ins`.
`rules::insertion_is_duplication` (reached by `normalize_rna`) performs the relabel; the
inverted-duplication sub-bullet (`:18`) is worked separately below. The spec's own invalid
tandem-as-insertion spelling `r.456_457ins123_456` is executable on the slice accession: ferro
expands the range payload to literal bases (`try_expand_rna_ins`, `INSERTED_SEQUENCE_EXPANDED`),
finds the payload is a copy of the 334 bases immediately 5', and relabels it `dup`.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[duplication-must-ranks-the-label-not-the-partition](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The rule that a duplication must be labelled 'dup' ranks the label of each piece ferro derives, not the partition; the one exception is a net-longer tandem copy of a multi-base motif, where the derivation is cut to expose the dup rather than merged into a delins.
<!-- why:END:duplication-must-ranks-the-label-not-the-partition -->

</details>

`duplication-must-ranks-the-label-not-the-partition` cites `DNA/duplication.md:18`; on the `r.`
axis the same MUST is carried verbatim by `RNA/duplication.md:19`, and a `DNA/` citation cannot
scope `r.` on its own — the record is surfaced here for its mechanism (the MUST ranks a *label*,
not a partition), with the `r.`-jurisdiction authority being the RNA twin.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.6_7insu` | recommended | `NM_004006.3:r.7dup` | the duplicating insertion `r.6_7insu` copies the `u` at `r.6` directly 3', so the MUST relabels it to the recommended `dup` (landing 3'-most, `r.7`) rather than keeping the `ins`. This is `:17`'s "not an insertion" case, executable |
| `NM_004006.3:r.456_457ins123_456` | recommended | `NM_004006.3:r.126_459dup` | the spec's own `class="invalid"` tandem-as-insertion spelling, on the slice accession. Ferro expands the `123_456` range payload to literal bases, relabels the directly-3'-flanking copy to `dup`, and the 3'rule (`:19`) then carries the 334-base duplication three positions 3': `r.457_459` reads `cag`, the same as `r.123_125`, and `r.460` (`g`) is where the match stops. Output is the 3'-most `dup`, as `:17` and `:19` together require. This placement agrees on both the committed slice and the full current-hg38 reference |
| `NM_004006.3:r.123_456dup` | recommended | — | the `dup` form the spec says the row above must be rewritten to. It parses and stays a valid `dup`, but its exact 3' placement is **reference-tier-dependent**, so no fixed placement is pinned: the committed slice 3'-shifts it to `r.126_459dup`, while the full-length current-hg38 transcript leaves it a fixed point at `r.123_456dup`. Ferro applies the long-span 3'-shift on the insertion→`dup` relabel path (the row above, which shifts on **both** tiers) but not on a plain long-span `dup` input over the full transcript — a real inconsistency in the long-span shift gate, flagged for the operator, no issue filed here |

**Confluence, with a caveat.** The pair `{r.456_457ins123_456, r.123_456dup}` is one variant by
`:17`'s own statement. On the committed slice both spellings converge on the 3'-most
`r.126_459dup`. On the full current-hg38 reference they **diverge**: the insertion still relabels
and shifts to `r.126_459dup`, but the plain `r.123_456dup` input is left unshifted at
`r.123_456dup` — ferro's long-span shift gate treats the two entry paths differently over a
full-length transcript. This is the flagged finding above; it is a `dup`/`dup` placement question,
not a validity one, so both spellings stay valid HGVS.

## `insertion.md:18` — an inverted duplication is `ins<range>inv`, not `dup`

> - **inverted duplications** are described as an insertion (`r.234_235ins123_234inv`), not as a duplication (see [Inversion](inversion.md)).

Ferro: an inverted duplication is spelled `ins<range>inv`, naming the span the inverted copy came
from — never `dup`, and never the `dupinv` shorthand, which is refused at parse. `:18` is the one
clause on this page the ledger cites directly on its own `r.` authority.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum, so a short chance reverse-complement match is not misread as one.
<!-- why:END:inverted-duplication-is-derived-as-ins-range-inv -->

</details>

**The inverted-duplication form is not yet derived on `r.`** The `ins<range>inv` re-spell is wired
only in `normalize_genome` (the sole production caller of `rules::inverted_adjacent_copy_span`);
`normalize_rna` has no equivalent, so ferro does not **mint** it from a reverse-complement literal
payload on the `r.` axis, and it **expands** an already-spelled range payload to literal
reverse-complement bases via `try_expand_rna_ins`. The ledger record discloses this gap itself
("wired only in `normalize_genome`") and names #1946's render stage as the intended fix. The output
stays literal on `r.` — valid HGVS that re-parses and denotes the same sequence, so the verdict is
`conformant`, not `bug`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.234_235ins123_234inv` | conformant | `NM_004006.3:r.234_235insguugacauuguucagggcaugaacucuuguggauccuuuuucuuuuggcaguuuuugcccugucaggccuucgaggaggucuaggaggcgccucccauccuguaggucacug` | the spec's own inverted-duplication spelling. On a real reference ferro **expands** the range payload to literal reverse-complement bases via `try_expand_rna_ins` rather than deriving/preserving the `ins<range>inv` form: the mint/keep of the range-inv spelling is wired only in `normalize_genome`, so the `r.` axis emits literals (`inverted-duplication-is-derived-as-ins-range-inv`; #1946 names the render-stage fix). Valid HGVS that re-parses and denotes the same sequence, so conformant, not the recommended spelling |
| `NM_004006.3:r.123_456dupinv` | refused | — | the `dupinv` shorthand — rejected at parse ("Unexpected trailing characters: 'inv'"); `:18` requires the `ins<range>inv` form. Pinned by `reject_dupins`/`recommended_form_pins` |

See also → `RNA/duplication.md:18-21` (the same `ins<range>inv` limitation, worked from the
duplication side); the `ins<range>inv` spelling itself is named by the inversion recommendations,
which `inverted-duplication-is-derived-as-ins-range-inv` cites at `DNA/inversion.md:69`.

## `insertion.md:19` — the 3'rule

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).

Ferro: the 3'rule applies to insertions exactly as to deletions and duplications — an inserted
sequence that could be written at more than one equivalent position is assigned its most 3'
placement. There is no exon/exon-junction NOTE on this page (contrast `RNA/deletion.md:18-20` and
`RNA/duplication.md:22-25`, whose NOTE and #2211 defect are specific to del/dup), so this clause is
CONFIRM-by-inspection.

A pure-insertion 3'-shift is hard to exhibit on this transcript's homopolymer runs: an inserted
base that matches an adjacent run base is by definition a copy of its 5' flank, so it is a
duplication (`:17`), not an insertion, and is relabeled `dup` before any insertion-level shift can
be shown. The executable row below is therefore a fixed-point insertion; the shiftable cases show
up on the duplication page, and the one long-span case where ferro misses the 3'rule is worked at
`:17` above.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_124insg` | recommended | self | the inserted `g` matches neither flank, so the 3'rule has nothing to shift it into — already 3'-most, a fixed point |

## `insertion.md:20-21` — the inserted sequence: literal bases, a reference range, or another reference

> - the **"inserted_sequence"** can be given as the nucleotides inserted (e.g., `insagc`) or, for larger insert sequences, by referring to the sequence in the reference sequence (e.g., `r.849_850ins858_895`) or another reference (see Examples).
>   When the inserted sequence is not present in a reference sequence, it should be submitted to a database (e.g., [GenBank](http://www.ncbi.nlm.nih.gov/genbank/submit/)); the accession.version number obtained can then be used to describe the variant, like `r.123_124ins[L37425.1:r.23_361]`.

Ferro: the spec admits three payload spellings — literal bases (`insagc`), a range-citation into
the same reference (`ins858_895`), and a cross-reference into another accession
(`ins[L37425.1:r.23_361]`). Ferro parses all three; on `r.` it **expands** a range-citation payload
to literal bases via `try_expand_rna_ins`, and resolves a cross-reference payload against the
named accession when the provider serves it (#422, `issue_422_cross_reference_ins`), declining
with `Reference not found` when it does not. Since `:20` states in as many words that both the
literal and the range-citation forms are legal, and no clause selects between them, which spelling
ferro emits is a representation choice governed by the ledger — not a defect.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.756_757insuacu` | recommended | self | a literal-base insert (the executable form of the spec's `LRG_199t1:r.756_757insuacu` at `:29-30`); the literal payload is the recommended spelling for a short insert. Measured fixed point on the slice |
| `NM_004006.3:r.849_850ins858_895` | conformant | `NM_004006.3:r.849_850insugagagaacuucuuccccuaagccucgauucaagagcu` | the spec's own range-citation payload, on the slice accession: ferro expands `858_895` to the 38 literal bases it names (`INSERTED_SEQUENCE_EXPANDED`). Both spellings are legal by `:20`; the expansion is a representation choice (`canonical-form-choice-when-both-legal`), and the range form is the one #1946's render stage would restore. Valid HGVS that re-parses |
| `LRG_199t1:r.849_850ins858_895` | conformant | — | the spec's own accession for the same range-citation form; not in the slice, so parse-only |
| `LRG_199t1:r.123_124ins[L37425.1:r.23_361]` | recommended | — | the spec's cross-reference payload; ferro parses `ins[ACC:...]` forms and expands them against the named reference (#422), but neither `L37425.1` nor `LRG_199t1` is in the slice — parse-only. On `NM_004006.3` the same payload parses and normalization declines with `Reference not found: L37425.1`, which is the correct answer for an unserved accession |

## `insertion.md:22` — uncertain: the dagger form

> - † = see [Uncertain](../uncertain.md); when the position and/or the sequence of an inserted sequence has not been defined, a description may have a format like `r.(100_150)insn[25]`.

Ferro: an uncertain insertion's parenthesised range and unspecified-length payload are preserved
verbatim — the 3'rule has nothing determinate to shift when neither the exact position nor the
exact sequence is asserted. `issue_1079_flanking_insertion_anchor.rs` pins the analogous
`g.(100_150)insN[25]` as accepted.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.(100_150)insn[25]` | recommended | self | the uncertain form, preserved (illustrative position; the spec's own example is generic). The parenthesised range and the `n[25]` unspecified payload are both left as written — the bracketed-count rewrite at `:38-42` does not fire under an uncertain position |

## `insertion.md:26-30` — literal-base examples

> - **`LRG_199t1:r.426_427insa`**<br>
>   the insertion of an `a` nucleotide between nucleotides `r.426` and `r.427`.
>
> - **`LRG_199t1:r.756_757insuacu`**<br>
>   the insertion of nucleotides `uacu` between nucleotides `r.756` and `r.757`.

Ferro: the spec's two literal-base worked examples — a single nucleotide and a short run, each
written with the literal payload, which is the recommended spelling for a short insert. Both sit on
`LRG_199t1`, so they are parse-only here; the executable twins on the slice accession are measured
fixed points.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:r.426_427insa` | recommended | — | spec's own single-nucleotide literal insertion (parse-only here) |
| `NM_004006.3:r.426_427insa` | recommended | self | executable twin on the slice accession — a fixed point (the inserted `a` copies neither flank) |
| `LRG_199t1:r.756_757insuacu` | recommended | — | spec's own multi-nucleotide literal insertion (parse-only here) |
| `NM_004006.3:r.756_757insuacu` | recommended | self | executable twin on the slice accession — a fixed point |

## `insertion.md:32-42` — unspecified nucleotides and an uncertain position

> - **`NM_004006.2:r.(222_226)insg` (`p.Asn75fs`)**<br>
>   the insertion of a `g` at an unknown position in the sequence encoding amino acid 75.
>
> - **`NM_004006.2:r.549_550insn`**<br>
>   the insertion of one not specified nucleotide (`n`) between positions `r.549` and `r.550`.
>
> - **`NM_004006.2:r.761_762insnnnnn` (alternatively `r.761_762insn[5]`)**<br>
>   the insertion of 5 not specified nucleotides (`nnnnn`) between positions `r.761` and `r.762`.
>
> - **`LRG_199t1:r.1149_1150insn[100]`**<br>
>   the insertion of 100 not specified nucleotides between positions `r.1149` and `r.1150`.

Ferro: `n` is a legal HGVS symbol (`general.md:49` admits IUPAC-IUBMB symbols), so an unspecified
insert parses and normalizes as itself, and the parenthesised `(222_226)` uncertain position is
preserved. The spelled-out run `insnnnnn` and its repeat-count alternative `insn[5]` are two
spellings of the same unspecified 5-mer; the spec offers both, and ferro rewrites a bracketed
exact count to the literal run (the #920 rewrite in `normalize_core`, which is what lets a
one-copy insert of the 5' flank become a `dup`). For the 5-mer that lands on the spec's primary
spelling; for the 100-mer it lands on a 100-`n` literal where the spec's own example uses the
count form — both legal by `:38`, and no clause selects.

**On the HGVS→SPDI path these `n`-repeat forms are refused, not expanded.** `insn[100]` states a
*length*, not bases (`n` is any of `a`/`c`/`g`/`u`), so emitting 100 literal `N`s would be a
storable-but-wrong SPDI triple; `hgvs_to_spdi` declines with `UnrepresentableInSpdi` at
`expand_repeat_unit` (`src/spdi/convert.rs`). That is a conversion-path house choice — it moves no
normalized `r.` string, and the rows below stay valid HGVS — so it is recorded here rather than
executed: the harness runs `normalize`, not the SPDI conversion.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[spdi-n-unit-repeat-refusal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — On the HGVS-to-SPDI path ferro refuses an 'N'-unit or 'N'-containing repeat rather than expand it to literal 'N' bases — a project choice, since 'N' states a length, not identified bases, and the recommendations do not reach SPDI conversion.
<!-- why:END:spdi-n-unit-repeat-refusal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:r.(222_226)insg` | recommended | — | insertion at an uncertain position inside a codon; the parenthesised range is preserved (parse-only here) |
| `NM_004006.3:r.(222_226)insg` | recommended | self | executable twin on the slice accession — the uncertain position is preserved, nothing is clamped or shifted |
| `NM_004006.2:r.549_550insn` | recommended | — | one unspecified nucleotide; `n` preserved (parse-only here) |
| `NM_004006.3:r.549_550insn` | recommended | self | executable twin — a single `n` is a fixed point |
| `NM_004006.2:r.761_762insnnnnn` | recommended | — | spelled-out unspecified 5-mer, the spec's primary spelling (parse-only here) |
| `NM_004006.3:r.761_762insnnnnn` | recommended | self | executable twin — the literal run is a fixed point |
| `NM_004006.3:r.761_762insn[5]` | recommended | `NM_004006.3:r.761_762insnnnnn` | the spec's "alternatively" spelling, rewritten to the literal run (#920) — which is `:38`'s primary spelling, so recommended |
| `LRG_199t1:r.1149_1150insn[100]` | recommended | — | 100 unspecified nucleotides via the repeat-count spelling, the spec's own accession (parse-only here) |
| `NM_004006.3:r.1149_1150insn[100]` | conformant | `NM_004006.3:r.1149_1150insnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn` | executable twin — the same #920 rewrite expands the count to a 100-`n` literal. Valid HGVS that re-parses and states the same thing; the spec's own spelling for this case is the count form and `:38` makes the two alternatives, so this is a representation choice (`canonical-form-choice-when-both-legal`) rather than a defect — recorded as `conformant` because the emitted spelling is not the one the spec chose for a 100-mer. On the SPDI path this input refuses (`spdi-n-unit-repeat-refusal`) |

## `insertion.md:44-47` — an intronic insertion from a splicing event

> - <code class="invalid">NG_012232.1(NM_004006.2):r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]</code><br>
>   the insertion of intronic nucleotides `r.2950-30` to `r.2950-12` and `r.2950-4` to `r.2950-1` between nucleotides `r.2949` and `r.2950` (caused by the deletion `NC_000023.10(NM_004006.2):c.2950-11_2950-5del`).
>   Alternative description <code class="invalid">r.2949_2950ins[2950-30_2950-12;uuag]</code>.<br>
>   **NOTE**: for more examples of variants affecting splicing, see [Splicing](splicing.md).

Ferro: a splicing-driven insertion whose payload cites *intronic* coordinates on a genomic-parent
accession (`NG_012232.1(NM_004006.2)`). The spec marks both spellings `class="invalid"` and
redirects to the splicing recommendations. **Ferro's strict parser accepts both** — the
`ins[…]` bracket grammar admits an offset-bearing range part — and normalization then declines
with `UnsupportedVariant` ("CDS-offset range (intronic or UTR-marker) is spec-undefined and not
yet supported by ferro; see follow-up", `rules.rs`), so no output is produced. Under
`absolute-prohibition-enforcement-stage` a `class="invalid"` spelling should be refused at parse
in strict mode; here it is accepted at parse and refused at normalize, with no repair in lenient
either. That is an unenforced prohibition, not a malformed output, so the rows are recorded as
`conformant`/parse-only rather than `refused` (which the harness checks mechanically and which
strict does not do). `cross_doc_compliance::xdoc_rna_invalid_intronic_anchor_insertion` pins the
round-trip parse, not a rejection. **No tracking issue found** for the parse-stage half — flagged
for the operator, not filed here.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NG_012232.1(NM_004006.2):r.2949_2950ins[2950-30_2950-12;2950-4_2950-1]` | conformant | — | the spec's own `class="invalid"` splicing example. Strict parse **accepts** it; normalization refuses at the `CdsPositionRange` payload part with `UnsupportedVariant`, so there is no output to grade. Parse-only here (the accession is outside the slice); the governing treatment lives with the splicing recommendations |
| `NG_012232.1(NM_004006.2):r.2949_2950ins[2950-30_2950-12;uuag]` | conformant | — | the spec's `class="invalid"` alternative spelling (one intronic range part plus a literal `uuag`). Same shape: strict parse accepts, normalization refuses at the offset-bearing part; parse-only here |

## `insertion.md:51-57` — Discussion: a single-position insertion is not allowed

> !!! note "Can I describe a variant as <code class="invalid">r.123insg</code>?"
>
>     No, since the description is not unequivocal, it is not allowed.
>     What does the description mean, the insertion of a `g` **at** position `r.123` or the insertion of a `g` **after** position `r.123`?
>
>     The situation becomes even more complex when, using a coding RNA reference sequence, a "-" character is used; e.g., <code class="invalid">r.-14insG</code>.
>     When the insertion is **after** nucleotide `r.-14`, is this position `r.-13` or `r.-15`?

Ferro: a single-position insertion anchor is ambiguous ("at" vs "after") and cannot be repaired
into a flanking pair, so ferro refuses it in every mode — the strict rejection is at parse (with a
structured `InvalidEdit` code citing the DNA twin of this Note: "single-position insertion not
allowed on r.; the anchor MUST be two adjacent positions"), and lenient/silent also fail because
there is nothing to normalize the ambiguity into. `reject_single_position_insertion.rs` pins the
refusal across `g.`/`c.`/`n.`/`r.`/`m.`/`o.` "by construct-symmetry".

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_124insg` | recommended | self | the correct two-position flanking form |
| `NM_004006.3:r.123insg` | refused | — | the disallowed single-position spelling (`class="invalid"`) — rejected at parse with a structured `InvalidEdit` code; unrepairable, so refused in lenient and silent too (the `cannot normalize` branch) |
| `NM_004006.3:r.-14insG` | refused | — | the spec's 5'UTR variant of the same question — rejected at parse by the same single-position check, before the `-13`/`-15` ambiguity the Note describes is ever reached |

## `insertion.md:59-63` — Discussion: the `^` character is not allowed

> !!! note "Can I use the "^" character to describe an insertion?"
>
>     No, insertions can not be described using the format <code class="invalid">r.123ˆ124insu</code> or <code class="invalid">r.123ˆ124u</code>.
>     The recommendations try to restrict the number of different characters used to a minimum.
>     Since a character was already used to indicate a range (the *underscore*), a new character was not required.

Ferro: the `^` separator is rejected at parse — the range underscore already spells the flanking
pair, so `^` is a character the grammar does not admit. `input_hygiene_rejections.rs::rejects_caret_separator_in_insertion`
pins both `g.123^124insG` and the bare `g.123^124G` form; the `r.` twins below fail the same way.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_124insu` | recommended | self | the correct underscore-flanked form |
| `NM_004006.3:r.123^124insu` | refused | — | the disallowed `^` separator (`class="invalid"`) — rejected at parse; use the underscore range |
| `NM_004006.3:r.123^124u` | refused | — | the spec's second `class="invalid"` form, with the `ins` keyword dropped as well — rejected at parse at the same character |

## `insertion.md:65-69` — Discussion: a non-tandem duplicative copy is an insertion

> !!! note "How should I describe the change `aucg`<code class="spot1">aucgaucgaucg</code><code class="spot2">a</code>`ggguccc` to `aucg`<code class="spot1">aucgaucgaucg</code><code class="spot2">a</code><code class="ins">aucgaucgaucg</code>`ggguccc`? The fact that the inserted sequence (<code class="ins">aucgaucgaucg</code>) is present in the original sequence, suggests it derives from a duplicative event."
>
>     The variant should be described as an insertion; `r.17_18ins5_16`.
>     A description using "dup" is not correct since, by definition, a duplication should be **directly 3'-flanking of the original copy** (in tandem).
>     Note that the description given still makes it clear that the sequence inserted between `r.17` and `r.18` is probably derived from nearby, i.e. positions `r.5` to `r.16`, and thus likely derived from a duplicative event.

Ferro: the worked case of `:5`/`:17` — a copy present in the reference but **not** directly
3'-flanking is an **insertion**, not a duplication, even when it plainly derives from a duplicative
event. Ferro will not relabel a non-3'-flanking copy as `dup` (`rules::insertion_is_duplication`
checks only the immediate 5' flank). The by-range payload (`ins5_16`) versus a literal payload is
the representation choice `:20` leaves open, governed by `canonical-form-choice-when-both-legal`;
ferro expands the range payload to literal bases on `r.` via `try_expand_rna_ins`, the same
`r.`-axis gap as the `ins<range>inv` form at `:18` (#1946's render stage is the shared home). This
is the verbatim twin of `RNA/duplication.md:47-51`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.17_18ins5_16` | conformant | `NM_004006.3:r.17_18insuuuggugggaag` | executable twin of the spec's own answer `r.17_18ins5_16` (accession-less, not table-executable) — an insertion naming its likely source range, **not** a `dup`, because the copy is not directly 3'-flanking. ferro will not relabel a non-3'-flanking copy as `dup`; it does expand the range payload to literal bases via `try_expand_rna_ins`, a representation choice left open by the spec (`insertion.md:20` admits both), governed by `canonical-form-choice-when-both-legal` — valid HGVS that re-parses, so conformant |

See also → `RNA/duplication.md:18-21` (the `dup`/`ins` boundary and the same `ins<range>inv`
limitation), `RNA/duplication.md:47-51` (the identical duplicative-event worked case),
`general.md` (the prioritisation the `:17` tandem-duplication rule points at).
