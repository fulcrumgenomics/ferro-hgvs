# Insertion — ferro's reading

ferro's reading of the HGVS **insertion** recommendations, clause by clause — each spelling with
the form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Insertion (`r.`)](../RNA/insertion.md).*

Three ledger records reach this page's clauses. `duplication-must-ranks-the-label-not-the-partition`
governs the `dup`/`ins` boundary the definition (`insertion.md:5`) and the tandem-duplication Note
(`insertion.md:17`) draw. `inverted-duplication-is-derived-as-ins-range-inv` governs the
inverted-duplication Note (`insertion.md:18`). `canonical-form-choice-when-both-legal` governs
which spelling ferro emits when a payload may be written as literal bases or as a reference-range
citation (`insertion.md:22-23`). The flanking/order Notes, the 3'rule Note, the separation
boilerplate and the literal-base Examples carry **no Why block** — the reading is
CONFIRM-by-inspection against the spec text and the shipped code, not an adjudicated ruling.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
examples are written on `LRG_199t1`, `NM_004006.2` and on foreign genomic accessions
(`NC_000023.10`, `NC_000002.11`, `NC_000004.11`, `NC_000006.11`), none of which the slice carries,
so those rows are parse-only (`—`) — ferro cannot read their bases here. The `NM_004006.3` base
facts the executable rows rely on are the ones `deletion.md`, `substitution.md` and the RNA twin
establish: `c.123` is `C` and `c.124` is `A`; `c.5_7` is the run `TTT` (so `c.6` and `c.7` are both
`T`); `c.5_16` reads `TTTGGTGGGAAG`; the transcript spans `c.-237` … `c.*2697`. Every executable
`NM_004006.3` row below shows ferro's actual normalized output, verified against the bless
harness.

## `insertion.md:5` — definition: inserted, and not a copy of the 5' flank

> Insertion: a sequence change where, compared to the reference sequence, one or more nucleotides are inserted **and** where the insertion is not a copy of a sequence immediately 5'.

Ferro: an insertion adds one or more nucleotides between two flanking positions, and the "not a
copy of a sequence immediately 5'" clause is exactly the `dup`/`ins` boundary — a copy directly
3'-flanking of its original is a duplication (see `insertion.md:17` and `DNA/duplication.md:5`),
so ferro relabels a duplicating insertion to `dup` rather than leaving it as an `ins`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | a genuine single-nucleotide insertion between `c.123` (`C`) and `c.124` (`A`); the inserted `G` copies neither flank, so it is neither a `dup` nor 3'-shiftable — a fixed point |
| `NM_004006.3:c.6_7insT` | recommended | `NM_004006.3:c.7dup` | the boundary case: the inserted `T` **is** a copy of the `T` at `c.6` immediately 5', so by the definition it is not an insertion — ferro relabels it to the recommended `dup` (`rules::insertion_is_duplication`, reached by `normalize_cds`), landing on the 3'-most `T` of the `c.5_7` run (`c.7`) |

## `insertion.md:15` — the flanking positions must be two *adjacent* nucleotides

> - `positions flanking` should contain **two flanking nucleotides**, e.g., `123_124`, not `123_125`.

Ferro: an insertion occupies the zero-width junction between two adjacent positions, so the anchor
must name a pair that is exactly one apart. A non-adjacent anchor is refused at parse with a
structured `InvalidEdit` error — there is no way to recover which flanking pair the author meant,
so this is a rejection, not a repair.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | correctly flanking (adjacent) pair |
| `NM_004006.3:c.123_125insG` | refused | — | non-adjacent flanks — rejected at parse ("the two positions MUST be adjacent"); unrepairable |

## `insertion.md:16` — the flanking positions are listed 5' to 3'

> - the `positions_flanking` should be listed from **5' to 3'**, e.g., `123_124`, not `124_123`.

Ferro: a reversed, non-circular flank order has no exception on the coding axis (the
`<high>_<low>` form is admitted only for `o.`/`m.` circular references), so it is refused rather
than reordered.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | correctly ordered flank |
| `NM_004006.3:c.124_123insG` | refused | — | reversed flank order on a non-circular reference — rejected at parse |

## `insertion.md:17` — a tandem duplication is a `dup`, not an insertion

> - tandem duplications are described as a duplication (`g.123_456`<code class="spot1">dup</code>), not an insertion (<code class="invalid">g.456_457ins123_456</code>, see [Prioritization](../general.md)).

Ferro: the `dup`-over-`ins` preference seen from the insertion side — a copy directly
3'-flanking of its original is spelled `dup`, so a duplicating insertion is relabeled rather than
kept as an `ins`. `rules::insertion_is_duplication` performs the relabel. The spec's own invalid
tandem-as-insertion spelling `g.456_457ins123_456` is accession-less here, so it is not
table-executable, but on a real reference ferro expands the range payload to literal bases, finds
the payload is a copy of the bases immediately 5', and relabels it `dup`. The inverted-duplication
sub-bullet (`insertion.md:18`) is worked separately below.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[duplication-must-ranks-the-label-not-the-partition](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The rule that a duplication must be labelled 'dup' ranks the label of each piece ferro derives, not the partition; the one exception is a net-longer tandem copy of a multi-base motif, where the derivation is cut to expose the dup rather than merged into a delins.
<!-- why:END:duplication-must-ranks-the-label-not-the-partition -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.6_7insT` | recommended | `NM_004006.3:c.7dup` | the duplicating insertion `c.6_7insT` copies the `T` at `c.6` immediately 5', so the relabel makes it the recommended `dup` (landing 3'-most on the `c.5_7` run, `c.7`) rather than keeping the `ins`. This is `:17`'s "not an insertion" case, executable |

## `insertion.md:18` — an inverted duplication is `ins<range>inv`, not `dup`

> - **inverted duplications** are described as an insertion (`g.234_235ins123_234inv`), not as a duplication (see [Inversion](inversion.md)).

Ferro: an inverted duplication is spelled `ins<range>inv`, naming the span the inverted copy came
from — never `dup`, and never the `dupinv` shorthand, which is refused at parse.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum, so a short chance reverse-complement match is not misread as one.
<!-- why:END:inverted-duplication-is-derived-as-ins-range-inv -->

</details>

**The inverted-duplication form is not yet derived on the coding axis.** The `ins<range>inv`
re-spell is wired only in `normalize_genome` (the sole caller of
`rules::inverted_adjacent_copy_span`); `normalize_cds` has no equivalent, so on the `c.` axis ferro
**expands** an already-spelled `ins<range>inv` payload to literal reverse-complement bases rather
than deriving or preserving the range-inv spelling. The ledger record discloses this gap itself
("wired only in `normalize_genome`") and names #1946's render stage as the intended fix. The output
stays literal on `c.` — valid HGVS that re-parses and denotes the same sequence, so the verdict is
`conformant`, not `bug`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.849_850ins850_900inv` | conformant | `NM_004006.3:c.849_850insGGCATAGCTCTTGAATCGAGGCTTAGGGGAAGAAGTTCTCTCATATCCCTG` | executable twin of the spec's inverted-duplication spelling: ferro expands the `850_900inv` payload to literal reverse-complement bases via the `INSERTED_SEQUENCE_EXPANDED` path rather than deriving/preserving the `ins<range>inv` form, because the mint/keep is wired only in `normalize_genome`. Valid HGVS that re-parses (`inverted-duplication-is-derived-as-ins-range-inv`; #1946) |
| `NM_004006.2:c.849_850ins850_900inv` | conformant | — | the spec's own accession (inverted copy inserted 5' of the original) — parse-only |
| `NM_004006.2:c.900_901ins850_900inv` | conformant | — | the spec's own accession (inverted copy inserted 3' of the original) — parse-only |
| `NM_004006.3:c.850_900dupinv` | refused | — | the `dupinv` shorthand — rejected at parse; `:18` requires the `ins<range>inv` form |

## `insertion.md:19-20` — individual variants, not a delins

> - two variants separated by one or more nucleotides should be described individually and **not** as a "delins".<br>
>   **exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a "delins".<br>

Ferro: this is the shared separation/codon clause (`general.md:33`/`:34`), reproduced verbatim
across nine files. A single `ins` variant has no internal separation, so this Note is boilerplate
about a *second, separate* variant appearing near the insertion — not about the insertion's own
grammar. The clause is adjudicated where the merge geometry actually lives (`deletion.md:18-19`,
`delins.md`), out of this file's scope; nothing in `insertion.md`'s own text bears on it
differently, so there is no insertion-specific question to re-litigate and no verdict row here.

## `insertion.md:21` — the 3'rule

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).

Ferro: the 3'rule applies to insertions exactly as to deletions and duplications — an inserted
sequence that could be written at more than one equivalent position is assigned its most 3'
placement. There is no exon/exon-junction NOTE on this page (contrast the deletion page's
`deletion.md:20-22`), so this clause is CONFIRM-by-inspection.

A pure-insertion 3'-shift is hard to exhibit on this transcript's homopolymer runs: an inserted
base that matches an adjacent run base is by definition a copy of its 5' flank, so it is a
duplication (`insertion.md:17`), not an insertion, and is relabeled `dup` before any
insertion-level shift can be shown. The executable row below is therefore a fixed-point insertion.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | the inserted `G` matches neither flank, so the 3'rule has nothing to shift it into — already 3'-most, a fixed point |

## `insertion.md:22-23` — the inserted sequence: literal bases, a reference range, or another reference

> - the **"inserted_sequence"** can be given as the nucleotides inserted (e.g., `insAGC`) or, for larger insert sequences, by referring to the sequence in the reference sequence (e.g., `c.849_850ins858_895`) or another reference (e.g., `NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]`).
>   When the inserted sequence is not present in the reference genome, it should be submitted to a database (e.g., [GenBank](http://www.ncbi.nlm.nih.gov/genbank/submit/)); the accession.version number obtained can then be used to describe the variant.

Ferro: the spec admits three payload spellings — literal bases (`insAGC`), a range-citation into
the same reference (`ins858_895`), and a cross-reference into another accession
(`ins[NC_…:g.A_B]`). Ferro parses all three; on the `c.` axis it **expands** a range-citation
payload to literal bases rather than preserving it. Since `:22` states in as many words that both
the literal and the range-citation forms are legal, and no clause selects between them, which
spelling ferro emits is a representation choice governed by the ledger — not a defect. The final
sentence is authoring workflow advice (submit novel sequence to a database, cite the accession),
not a normalization obligation.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | a literal-base insert — the recommended spelling for a short insert |
| `NM_004006.2:c.849_850ins858_895` | conformant | — | the spec's own same-reference range-citation payload — parse-only; ferro expands the range to literal bases on a real reference (`canonical-form-choice-when-both-legal`; the range form is what #1946's render stage would restore) |
| `NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]` | conformant | — | the spec's cross-reference citation payload — parse-only |

## `insertion.md:24` — uncertain: the dagger form

> - † = see [Uncertain](../uncertain.md); when the position and/or the sequence of an inserted sequence has not been defined, a description may have a format like `g.(100_150)insN[25]`.

Ferro: an uncertain insertion's parenthesised range and unspecified-length payload are preserved
verbatim — the 3'rule has nothing determinate to shift when neither the exact position nor the
exact sequence is asserted.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.(100_150)insN[25]` | recommended | self | the uncertain form, preserved (illustrative position; the spec's own example is generic genomic). The parenthesised range and the `N[25]` unspecified payload are both left as written |

## `insertion.md:28-39` — Examples: simple insertions

> - **`NM_004006.2:c.849_850ins858_895`**<br>
>   the insertion of a copy of nucleotides `c.858` to `c.895` between nucleotides `c.849` and `c.850`.

Ferro: the spec's four simple worked forms span the three payload kinds — literal bases, a
same-reference range citation, and a cross-reference range citation — at the simplest complexity
level. Each is the canonical worked form; all sit on foreign accessions or a foreign transcript
version, so they are parse-only here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.10:g.32867861_32867862insT` | recommended | — | spec's single-nucleotide literal insertion — parse-only |
| `NC_000023.10:g.32862923_32862924insCCT` | recommended | — | spec's multi-nucleotide literal insertion — parse-only |
| `NM_004006.2:c.169_170insA` | recommended | — | spec's coding literal insertion — parse-only |
| `NM_004006.2:c.849_850ins858_895` | conformant | — | spec's same-reference range citation — parse-only (expands to literal bases on a real reference) |
| `NC_000002.11:g.47643464_47643465ins[NC_000022.10:g.35788169_35788352]` | conformant | — | spec's cross-reference citation — parse-only |
| `NM_004006.3:c.123_124insG` | recommended | self | executable literal-insertion twin on the slice accession — a fixed point |

## `insertion.md:41-49` — Examples: complex insertions

> - **`NM_004006.2:c.419_420ins[T;401_419]`**<br>
>   the insertion of `T` followed by a copy of the sequence from `c.401` to `c.419` (a duplication not directly flanking the original sequence).

Ferro: bracketed composite payloads — a semicolon-joined list mixing literal bases, same-reference
range citations, and cross-reference-plus-repeat parts in one payload. Nothing here is ambiguous in
meaning; each part is independently well-defined, and concatenation follows 5'→3' order of the
resulting inserted sequence. Ferro parses all three; a range part expands to literal bases on a
real reference (`canonical-form-choice-when-both-legal`). All sit on foreign accessions, so they
are parse-only here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.419_420ins[T;401_419]` | conformant | — | literal + same-reference range composite — parse-only |
| `LRG_199t1:c.419_420ins[T;450_470;AGGG]` | conformant | — | literal + range + literal composite — parse-only |
| `NC_000006.11:g.10791926_10791927ins[NC_000004.11:g.106370094_106370420;A[26]]` | conformant | — | cross-reference Alu range + `A[26]` repeat composite — parse-only |

## `insertion.md:51-62` — Examples: insertion of inverted duplicated copies

> - **`NM_004006.2:c.849_850ins850_900inv`**<br>
>   a copy of nucleotides `c.850` to `c.900` is inserted, in inverted orientation, 5' of the original sequence, between nucleotides `c.849` and `c.850`.

Ferro: the direct worked examples for `insertion.md:18`'s ruling (5' vs 3' of the original copy),
plus two intricate composites — an inverted-copy part combined with a substitution (`A`) or with a
second inverted-copy part carrying an implied internal deletion (a "hole" in the cited range). The
`ins<range>inv` part-type is the derived form on the genomic axis; on the `c.` axis ferro expands
it to literal reverse-complement bases (see `insertion.md:18`). All spec examples sit on foreign
accessions, so they are parse-only; the executable twin below is on the slice accession.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.849_850ins850_900inv` | conformant | — | inverted copy inserted 5' of the original — parse-only |
| `NM_004006.2:c.900_901ins850_900inv` | conformant | — | inverted copy inserted 3' of the original — parse-only |
| `LRG_199t1:c.940_941ins[885_940inv;A;851_883inv]` | conformant | — | inverted copy with an internal `G>A` substitution — parse-only |
| `NM_004006.2:c.940_941ins[903_940inv;851_885inv]` | conformant | — | inverted copy with an internal deletion (a "hole" in the cited range) — parse-only |
| `NM_004006.3:c.849_850ins850_900inv` | conformant | `NM_004006.3:c.849_850insGGCATAGCTCTTGAATCGAGGCTTAGGGGAAGAAGTTCTCTCATATCCCTG` | executable twin — ferro expands the range-inv payload to literal reverse-complement bases on the `c.` axis (`inverted-duplication-is-derived-as-ins-range-inv`; #1946) |

## `insertion.md:64-87` — Examples: incomplete descriptions

>     - **`NM_004006.2:c.761_762insNNNNN` (alternatively `NM_004006.2:c.761_762insN[5]`)**<br>
>     the insertion of 5 not specified nucleotides (`NNNNN`) between positions `c.761` and `c.762`.

Ferro: the uncertainty forms — a `G` at a parenthesised uncertain position, and runs of
unspecified `N` bases given either spelled out (`insNNNNN`) or as a repeat count (`insN[5]`,
`insN[100]`, `insN[(80_120)]`, `insN[?]`). `N` is a legal HGVS symbol (`general.md:47` admits
IUPAC-IUBMB symbols), so an unspecified insert parses and normalizes as itself, and a parenthesised
uncertain position is preserved. The spec states `insNNNNN` and `insN[5]` as equivalent
("alternatively"); ferro rewrites the bracketed exact count to the literal run (the #920 rewrite),
so `insN[5]` lands on the spec's primary spelling `insNNNNN`. This `NNNNN`≡`N[5]` equivalence is
stated by the spec but not currently pinned as an `equivalence_classes` entry in the ledger — a
candidate for one (flagged for the operator). Most spec examples are genomic/foreign, so parse-only;
the executable twins are on the slice accession.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.2:c.(222_226)insG` | recommended | — | insertion at an uncertain position inside a codon; the parenthesised range is preserved — parse-only |
| `NM_004006.3:c.(222_226)insG` | recommended | self | executable twin — the uncertain position is preserved, nothing clamped or shifted |
| `NC_000004.11:g.(3076562_3076732)insN[12]` | recommended | — | 12 unspecified nucleotides at an uncertain position — parse-only |
| `NC_000023.10:g.32717298_32717299insN` | recommended | — | one unspecified nucleotide — parse-only |
| `NM_004006.2:c.761_762insNNNNN` | recommended | — | spelled-out unspecified 5-mer, the spec's primary spelling — parse-only |
| `NM_004006.3:c.761_762insNNNNN` | recommended | self | executable twin — the literal run is a fixed point |
| `NM_004006.3:c.761_762insN[5]` | recommended | `NM_004006.3:c.761_762insNNNNN` | the spec's "alternatively" spelling, rewritten to the literal run (#920) — which is the spec's primary spelling, so recommended |
| `NC_000023.10:g.32717298_32717299insN[100]` | recommended | — | 100 unspecified nucleotides via the repeat-count spelling — parse-only |
| `NC_000023.10:g.32717298_32717299insN[(80_120)]` | recommended | — | 80–120 unspecified nucleotides (uncertain count) — parse-only |
| `NC_000023.10:g.32717298_32717299insN[?]` | recommended | — | an unknown number of unspecified nucleotides — parse-only |
| `NC_000006.11:g.8897754_8897755ins[N[543];8897743_8897754]` | conformant | — | an undefined 543-nt run plus a 12-nt target-site-duplication range citation; like the `ins[T;401_419]` composite above, ferro would expand the same-reference range citation to literal bases on a real reference rather than preserve it, so it is not the recommended range-citation form (parse-only here) |

## `insertion.md:89-91` — Examples: other

> - **`g.?_?ins[NC_000023.10:g.(12345_23456)_(34567_45678)]`**<br>
>   the insertion of a sequence from the X-chromosome (`NC_000023.10`), maximally involving nucleotides `12345_45678` but certainly nucleotides `23456_34567`, at an unknown position (`g.?_?`) in the genome (see [Uncertain](../uncertain.md)).

Ferro: the doc's single most complex worked form — an entirely unknown insertion site (`g.?_?`)
combined with a cross-reference payload whose own span is doubly uncertain. It is legal by
composition of two already-legal conventions (the uncertain-range and cross-reference citations),
and it **parses**. It carries no accession of its own, so it is not table-executable here; and on a
real reference ferro's normalizer currently **declines** to expand this exact combination
(cross-reference payload plus uncertain range), naming its own scope limit ("Expansion currently
covers g./m./o./c./n./r. axes with simple positive-integer positions or ranges, optionally with a
trailing `inv` on a range") rather than mis-normalizing. This is the one worked example on the page
that does not round-trip through the normalizer — a real, locatable implementation gap (no clause
is in tension; ferro has simply not built the combination). No ledger record names it; flagged for
the operator as an issue candidate, not a ledger question.

## `insertion.md:95-101` — Discussion: a single-position insertion is not allowed

> !!! note "Can I describe a variant as <code class="invalid">g.123insG</code>?"
>
>     No, since the description is not unequivocal, it is not allowed.
>     What does the description mean, the insertion of a `G` **at** position `g.123` or the insertion of a `G` **after** position `g.123`?
>
>     The situation becomes even more complex when, using a coding DNA reference sequence, a "-" character is used; e.g., <code class="invalid">c.-14insG</code> or <code class="invalid">c.456-13insG</code>.
>     In the description <code class="invalid">c.456-13insG</code>, when the insertion is **after** intronic nucleotide `c.456-13`, is this position `c.456-12` or `c.456-14`?

Ferro: a single-position insertion anchor is ambiguous ("at" vs "after") and cannot be repaired
into a flanking pair, so ferro refuses it in every mode — the strict rejection is at parse (with a
structured `InvalidEdit` code). The intronic-offset variant (`c.456-13insG`) is refused by the same
single-position check, before the `-12`/`-14` ambiguity the Note describes is ever reached.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | the correct two-position flanking form |
| `NM_004006.3:c.123insG` | refused | — | the disallowed single-position spelling (`class="invalid"`) — rejected at parse; the "at"/"after" ambiguity is unrepairable |
| `NM_004006.3:c.-14insG` | refused | — | the spec's 5'UTR variant of the same question — rejected at parse by the same single-position check |
| `NM_004006.3:c.456-13insG` | refused | — | the spec's intronic variant — rejected at parse (single-position); the `-12`/`-14` ambiguity is never reached |

## `insertion.md:103-108` — Discussion: the `^` character is not allowed

> !!! note "Can I use the "^" character to describe an insertion?"
>
>     No, insertions can not be described using the format <code class="invalid">g.123ˆ124insG</code> or <code class="invalid">g.123ˆ124G</code>.
>     The recommendations try to restrict the number of different characters used to a minimum.
>     Since a character was already used to indicate a range (the *underscore*), a new character was not required.

Ferro: the `^` separator is rejected at parse — the range underscore already spells the flanking
pair, so `^` is a character the grammar does not admit. The bare form with the `ins` keyword
dropped fails the same way.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.123_124insG` | recommended | self | the correct underscore-flanked form |
| `NM_004006.3:c.123^124insG` | refused | — | the disallowed `^` separator (`class="invalid"`) — rejected at parse; use the underscore range |
| `NM_004006.3:c.123^124G` | refused | — | the spec's second `class="invalid"` form, with the `ins` keyword dropped too — rejected at parse |

## `insertion.md:109-113` — Discussion: a non-tandem duplicative copy is an insertion

> !!! note "How should I describe the change `ATCG`<code class="spot1">ATCGATCGATCG</code><code class="spot2">A</code>`GGGTCCC` to `ATCG`<code class="spot1">ATCGATCGATCG</code><code class="spot2">A</code><code class="ins">ATCGATCGATCG</code>`GGGTCCC`? The fact that the inserted sequence (<code class="ins">ATCGATCGATCG</code>) is present in the original sequence suggests it derives from a duplicative event."
>
>     The variant should be described as an insertion; `g.17_18ins5_16`.
>     A description using "dup" is not correct since, by definition, a duplication should be **directly 3'-flanking of the original copy** (in tandem).
>     Note that the description given still makes it clear that the sequence inserted between `g.17` and `g.18` is probably derived from nearby, i.e. positions `g.5` to `g.16`, and thus likely derived from a duplicative event.

Ferro: the worked case of the definition (`insertion.md:5`) and `:17` — a copy present in the
reference but **not** directly 3'-flanking is an **insertion**, not a duplication, even when it
plainly derives from a duplicative event. Ferro will not relabel a non-3'-flanking copy as `dup`
(`rules::insertion_is_duplication` checks only the immediate 5' flank). The final sentence is about
provenance ("probably derived from"), which no normalizer can recover — the correct description is
stated flatly as `ins` regardless. The by-range payload (`ins5_16`) versus a literal payload is the
representation choice `:22-23` leaves open, governed by `canonical-form-choice-when-both-legal`;
ferro expands the range payload to literal bases on the `c.` axis.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.17_18ins5_16` | conformant | `NM_004006.3:c.17_18insTTTGGTGGGAAG` | executable twin of the spec's own answer `g.17_18ins5_16` (accession-less, not table-executable) — an insertion naming its likely source range, **not** a `dup`, because the copy is not directly 3'-flanking. ferro expands the range payload (`c.5_16` = `TTTGGTGGGAAG`) to literal bases, a representation choice `insertion.md:22` leaves open (`canonical-form-choice-when-both-legal`) — valid HGVS that re-parses |

## `insertion.md:115-119` — Discussion: `c.23ins24` and the `c.9_32dup` correction

> !!! note "A variant in the _CDKN2A_ gene, duplicating the first 24 nucleotides of the coding DNA reference sequence, has been described as <code class="invalid">c.23ins24</code>. My interpretation is it should be described as `c.1_24dup`, is this correct?"
>
>     Since the sequence in that region is <span class="sequence">cagc<code class="spot1">ATGGAGCC</code>GGCGGCGGGGAGCAGC<code class="spot1">ATGGAGCC</code>TTCG</span>, the correct description is `c.9_32dup` (`p.(Ala4_Pro11dup)`).
>     `c.1_24dup` seems correct but neglects the **3'rule** (3' shift possible for the highlighted region).
>     `c.23ins24` is not correct since the position of the insertion is not described properly and because "ins24" does not define the sequence inserted.

Ferro: a compound worked example touching three rules at once — the single-position/adjacency
requirement (`c.23` is one position, not a flanking pair, so "the position … is not described
properly"), the absolute prohibition on a bare count (`checklist.md:33` — `ins24` names no bases), and the 3'rule's
interaction with duplication placement (`c.1_24dup` must be 3'-shifted within the repeated
`ATGGAGCC…ATGGAGCC` unit to `c.9_32dup`). The disallowed spelling `c.23ins24` is refused at parse
on both grounds. The `c.1_24dup`→`c.9_32dup` correction is on the _CDKN2A_ coding sequence, not on
`NM_004006.3`'s bases, so the `dup` half is conceptual/parse-only here and its exact 3'-shift is
not reproduced on the slice.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.23ins24` | refused | — | the disallowed spelling (`class="invalid"`) — a single-position anchor (`c.23` is not a flanking pair) and a bare count (`ins24` defines no bases); rejected at parse |
