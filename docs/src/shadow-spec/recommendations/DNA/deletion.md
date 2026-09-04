# Deletion — ferro's reading

ferro's reading of the HGVS **deletion** recommendations, clause by clause — each spelling with the
form ferro normalizes it to and a verdict on that output. New here? See
[How to read a page](../../reading-guide.md) for the verdicts, the table conventions, and the
recurring terms.

*RNA twin: [Deletion (`r.`)](../RNA/deletion.md).*

Most of this page is CONFIRM-by-inspection against the spec text and the shipped code: the plain
3'rule, the range-order and two-position rules, and the payload/length-suffix prohibitions are
mechanical parser or W-code behaviour with no clause in tension. Three units are adjudicated and
carry a Why block: the separation-and-codon clause at `deletion.md:18-19`, the exon/exon-junction
3'rule exception at `deletion.md:20-22`, and the length-suffix prohibition at `deletion.md:114-117`.
The exon/exon exception is the one place worth watching against the RNA twin — **on the coding axis
the exception applies and ferro is correct** (a `c.`/`n.` deletion is *not* shifted across the
junction), the opposite of the `r.` axis, where the same NOTE switches the exception off and ferro
does not yet honour it ([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211), on the
RNA deletion page).

Executable rows use `NM_004006.3`, the one transcript in the committed slice. It carries the spec's
own worked 3'rule example: the 8-nucleotide A-stretch `c.5690_5697` (`ATTG`·`AAAAAAAA`·`TT` at
`c.5686_5699`) that `deletion.md:33-34` describes on `NM_004006.2`, so the `c.5690del`/`c.5697del`
pair and the length-suffix repair below are the spec's own bases, one transcript version later.
The slice's transcript begins at `c.-237` and ends at `c.*2697` (the boundary rows at
`deletion.md:67-96` use those, where the spec's `NM_004006.2` runs `c.-244` to `c.*2691`). The
CDS base facts the separation rows rely on (`c.76` is `A`; the codon `c.145_147` is `CGC`;
`c.123` is `C`) are the ones `substitution.md` and `alleles.md` establish. The spec's other examples
sit on `LRG_199t1`, `NM_004006.2`, `NM_000492.3`, `NG_012232.1` and foreign genomic accessions
(`NC_000023.11`, `NC_000023.10`, `NC_000003.12`), none of which the slice carries, so those rows
are parse-only (`—`) — ferro cannot read their bases here.

## `deletion.md:5` — definition

> Deletion: a sequence change where, compared to a reference sequence, one or more nucleotides are not present (deleted).

Ferro: a deletion removes one or more nucleotides; the description carries only the deleted range,
never the deleted bases (the payload-spelling prohibition is at `deletion.md:27-40` below).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5697del` | recommended | self | single-nucleotide deletion — the last `A` of the 8-nucleotide A-stretch `c.5690_5697` (`c.5698` is `T`), so already a fixed point |
| `NC_000023.11:g.33344591del` | recommended | — | the spec's own one-nucleotide genomic example (foreign accession — parse-only here) |

## `deletion.md:15` — the range must name two different positions

> - `position(s)_deleted` should contain **two different positions**, e.g., `123_126`, not `123_123`.

Ferro: rule 2 (lowercase "should"); a same-position range is repaired to the single-position form.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5695_5697del` | recommended | self | a correctly ranged three-nucleotide deletion — the 3'-most three of the A-stretch, so already a fixed point |
| `NM_004006.3:c.5697_5697del` | refused | — | same-position range — rejected at parse (`W3016 SinglePositionRange`) in strict mode; lenient repairs to `c.5697del` (`correct_single_position_range`) |

## `deletion.md:16-17` — the range is listed 5' to 3', except on a circular reference

> - the `position(s)_deleted` should be listed from **5' to 3'**, e.g., `123_126`, not `126_123`.
>     - **exception**: when a circular genomic reference sequence is used ("o" and "m" prefix) nucleotide positions may be listed from 3' to 5' when the deletion includes both the last and first nucleotides of the reference sequence.

Ferro: the ordinary range is listed 5'→3'; a reversed, non-circular range is refused rather than
reordered. The circular carve-out (SVD-WG006) admits the reversed `<high>_<low>` spelling only on an
`o.`/`m.` reference for an origin-spanning deletion, and ferro honours it there (the reversed-range
allowance in `src/normalize/footprint.rs`); the two endpoints of such a deletion are checked
independently rather than order-compared.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5695_5697del` | recommended | self | correctly ordered range (`c.5695` is 5' of `c.5697`) |
| `NM_004006.3:c.5697_5695del` | refused | — | reversed range on a non-circular reference — rejected at parse |

## `deletion.md:18-19` — individual vs. delins, and the one-codon exception

> - two variants separated by one or more nucleotides should be described individually and **not** as a "delins".
>     - **exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a "delins".<br>

Ferro: the separation rule (ruleset rule 2 — a preference, not a ban); the exception folds two
changes into one delins only when they sit one nucleotide apart and together change one amino acid.
This is the shared separation clause (`general.md:33`/`:34`), reproduced verbatim across nine files;
the deletion-specific merge geometry is adjudicated on `delins.md`, out of this file's scope.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[separation-rule-force-modal-or-negation](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero.
>
> **[codon-carve-out-shape-restriction](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Two changes one nucleotide apart that together alter a single amino acid are written as one delins, whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled.
<!-- why:END:separation-rule-force-modal-or-negation,codon-carve-out-shape-restriction -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.[76A>G;123C>T]` | recommended | self | two changes more than one nucleotide apart — described individually, kept as a cis allele, not merged |
| `NM_004006.3:c.[145C>T;147C>G]` | recommended | `NM_004006.3:c.145_147delinsTGG` | one nucleotide apart, together altering the `c.145_147` codon (`CGC`→`TGG`) — folded to one delins, whatever the member edit types |

## `deletion.md:20-22` — the 3'rule, and the exon/exon-junction exception (where ferro is correct on DNA)

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily assigned to have been changed (**3'rule**).
>     - **exception**: deletions around exon/exon junctions when identical nucleotides flank the junction (see [Numbering](../../background/numbering.md#DNAc));<br>
>       when `..GA`<code class="del">T</code>`gta..//..cagTCA..` changes to `..GAgta..//..cagTCA..`, based on a coding DNA reference sequence, the variant is described as `LRG_199t1:c.3921del` (`NC_000023.10:g.32459297del`) and not as `c.3922del` (which would translate to `g.32456507del`).

Ferro: the 3'rule shifts a deletion in a single-residue stretch or tandem repeat to its most-3'
position. The **exception** clamps that shift at an exon/exon junction: a `c.`/`n.` deletion whose
identical flanking bases would otherwise let it shift across the junction is held at the last
position that does not cross, reached from **either** side — a description approaching the junction
stops at it, and one already spelled on the far side (`c.3922del`) is pulled back to `c.3921del`.

**This is where the coding axis and the RNA axis diverge, and ferro is correct here.** On `c.`/`n.`
the exception applies (`src/normalize/mod.rs`'s `edit_is_del_or_dup` clamp), so `c.3922del` does not
translate to the wrong nucleotide in the wrong exon. On `r.` the same NOTE *switches the exception
off* and the 3'rule must cross the junction — which ferro does not yet do
([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211), on the RNA deletion page).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[exon-junction-dup-converge-from-the-far-side](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — A duplication is placed at the most 3' position that does not cross an exon/exon junction, reached from either side, so a copy spelled past the junction is pulled back to it.
<!-- why:END:exon-junction-dup-converge-from-the-far-side -->

</details>

The governing record's title and worked example read as duplication-only (`c.3922dup`→`c.3921dup`),
but its SCOPE paragraph and the shipped clamp both cover deletion on the `c.`/`n.` axes — the code
gates on `edit_is_del_or_dup`, not on `dup` alone. Broadening the record's title/question to name
deletion is a proposed cosmetic strengthen for discoverability, not a ruling change, and is not done.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.5690del` | recommended | `NM_004006.3:c.5697del` | 3'rule over a single-residue stretch: `c.5690_5697` is `AAAAAAAA`, so the deletion lands on the last `A` — the spec's own `deletion.md:33-34` example, on the slice's transcript version |
| `NM_004006.3:c.5697del` | recommended | self | already the 3'-most base of that run — a fixed point |
| `NM_004006.3:c.5690_5692del` | recommended | `NM_004006.3:c.5695_5697del` | the same rule over a multi-nucleotide deletion inside the run — shifted to the 3'-most three |
| `LRG_199t1:c.3921del` | recommended | — | the spec's own worked example — the deletion held at the exon/exon junction (parse-only here) |
| `LRG_199t1:c.3922del` | recommended | — | the far-side spelling, which the exception pulls back to `c.3921del` on the coding axis (parse-only here — no `LRG_199t1` in the slice; both directions of the `c.`/`n.` deletion clamp are pinned on a synthetic two-exon transcript in `tests/it/issue_1621_exon_junction_far_side.rs`, and the near-side halt in `tests/it/issue_334_exon_junction_exception.rs`) |

No executable `NM_004006.3` row can exercise the exon/exon clamp: the committed slice is built as a
single flat exon, so no deletion on it reaches an exon/exon junction.

See also → `deletion.md:119-122` (the same exception restated as a Q&A).

## `deletion.md:23` — uncertain deletions

> - † = see [Uncertain](../uncertain.md); when the position and/or the sequence of a deletion has not been defined, a description may have a format like `g.(100_150)delN[15]`.

Ferro: an uncertain deletion's parenthesised range is preserved verbatim — the 3'rule has nothing
determinate to shift when neither the exact position nor the exact length is asserted. The
`delN[15]` shape parses (`src/hgvs/parser/edit.rs`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.(100_150)delN[15]` | recommended | self | uncertain form, preserved (illustrative position; the spec's own example is generic genomic) |

## `deletion.md:27-40` — one nucleotide (three worked examples)

>     - **`NM_004006.2:c.5697del` (3'rule)**<br>

Ferro: three sub-examples — a plain deletion (`g.33344591del`), the 3'rule over an A-stretch
(`c.5697del`, the last A of an 8-nucleotide run), and the same underlying variant expressed on the
minus strand (`g.32343183del` ↔ `c.5697del`, confirming that minus-strand mapping and the 3'rule
compose). Each carries a NOTE that the deleted base is **never** spelled (`delA` is `class="invalid"`).
The A-stretch is present unchanged on the slice's `NM_004006.3`, so the middle example — including
the `class="invalid"` `c.5690del` spelling the minus-strand NOTE at `deletion.md:40` names — runs
here on the spec's own bases.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.33344591del` | recommended | — | the plain one-nucleotide example (parse-only here) |
| `NM_004006.2:c.5697del` | recommended | — | the 3'rule over the A-stretch; the spec's transcript version (parse-only here) |
| `NM_004006.3:c.5697del` | recommended | self | executable twin on the slice's version — the last `A` of `c.5690_5697`, a fixed point |
| `NM_004006.3:c.5690del` | recommended | `NM_004006.3:c.5697del` | the `class="invalid"` first-`A` spelling of `deletion.md:40` — the 3'rule carries it to the last `A` |
| `NC_000023.11:g.32343183del` | recommended | — | the minus-strand cross-check (parse-only here) |
| `NC_000023.11:g.33344591delA` | refused | — | the disallowed payload-bearing spelling (`class="invalid"`) (parse-only here) |
| `NM_004006.3:c.5697delA` | refused | — | executable twin — strict rejects the stated deleted base; lenient drops the payload and normalizes to `c.5697del` |

## `deletion.md:42-49` — several nucleotides

>     - **`NC_000023.11:g.33344590_33344592del`**<br>

Ferro: a multi-nucleotide range never carries the deleted payload sequence either (`delGAT` is
`class="invalid"`), and an ordinary range crossing an exon/intron border needs no special handling
while it stays within the transcript.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.33344590_33344592del` | recommended | — | plain range, no deleted sequence spelled (parse-only here) |
| `NM_004006.3:c.5695_5697del` | recommended | self | executable several-nucleotide twin — the 3'-most three of the A-stretch, nothing left to shift |
| `NC_000023.11:g.33344590_33344592delGAT` | refused | — | the disallowed payload-bearing range spelling (`class="invalid"`) (parse-only here) |
| `NM_004006.3:c.5695_5697delAAA` | refused | — | executable twin — strict rejects the stated deleted sequence; lenient drops it and normalizes to `c.5695_5697del` |
| `NC_000023.11(NM_004006.2):c.183_186+48del` | recommended | — | range crossing an exon/intron border (parse-only here) |

## `deletion.md:51-65` — exon/intron/exon

>         - **`LRG_199t1:c.1704+1del`**<br>

Ferro: three border cases. The exon/exon case (`c.3921del`) is the exon/exon-junction exception
already adjudicated at `deletion.md:20-22`. The exon/intron (`c.1704+1del`, not `c.1704del`) and
intron/exon (`c.1813del`, not `c.1813-1del`) cases are ordinary "spell the deleted base where it
unambiguously sits" — no exception fires, because no run of identical bases straddles the border.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.1704+1del` | recommended | — | exon/intron border — the intronic `G` spelled as `c.1704+1`, not `c.1704del` (parse-only here) |
| `LRG_199t1:c.1813del` | recommended | — | intron/exon border — the exonic `G` spelled as `c.1813`, not `c.1813-1del` (parse-only here) |

## `deletion.md:67-96` — exons, and the beyond-transcript prohibition

>       following current recommendations (see [Numbering](../../background/numbering.md)), it is not allowed to describe variants in nucleotides beyond the boundaries of a reference sequence.

Ferro: multi-exon deletions are written with concrete or uncertain (parenthesized) positions that
stay within the transcript. The two `class="invalid"` forms describe a deletion running **beyond the
transcript's own boundaries** (5' of `c.-244` or 3' of `c.*2691` on `NM_004006.2`); that is a rank-1
prohibition — genomic coordinates are required instead. Ferro enforces the boundary for a
*concrete* transcript position: a `c.` coordinate past the CDS or transcript end is `W4004
PositionPastEnd` (refused in strict), and an `n.` position past either end is `W4008`
(`src/hgvs/noncoding_zones.rs`, which cites the `background/numbering.md` clause forbidding a
position beyond a transcript's boundaries). The two spec forms are different: each
names only in-bounds concrete positions and reaches past the boundary through an **open** `?`
end, and ferro accepts them — the fact that `c.-244` *is* the first nucleotide is the reference's,
invisible at parse, and no normalize-time check reads the `?`-open end against the transcript's
extent. The size of the deletion is never appended (`delXXXXX`; the length-suffix prohibition is
adjudicated at `deletion.md:114-117`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11(NM_004006.2):c.4072-1234_5155-246del` | recommended | — | concrete multi-exon deletion (parse-only here) |
| `NC_000023.11(NM_004006.2):c.(4071+1_4072-1)_(5154+1_5155-1)del` | recommended | — | the "exon-based" uncertain-position form (parse-only here) |
| `NC_000023.11(NM_004006.2):c.(3996_4196)_(5090_5284)del` | recommended | — | the probe-based uncertain-position form (parse-only here) |
| `NC_000023.11(NM_004006.2):c.(?_-244)_(31+1_32-1)del` | conformant | — | extends 5' beyond the transcript boundary (`class="invalid"`) — genomic coordinates required. Ferro **accepts** it: strict parses it, since the boundary is a reference fact, and nothing at normalize refuses the `?`-open end (parse-only here; no tracking issue yet) |
| `NC_000023.11(NM_004006.2):c.(10086+1_10087-1)_(*2691_?)del` | conformant | — | extends 3' beyond the transcript boundary (`class="invalid"`) — genomic coordinates required. Same acceptance as the row above (parse-only here; no tracking issue yet) |

## `deletion.md:98-103` — gene-level deletion

>     - **`NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del`**<br>

Ferro: whole-gene deletions in uncertain-range genomic syntax (SNP-array, MLPA) — a "does the parser
accept this shape" question, answered by the general uncertain-range machinery; no special
normalization.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del` | recommended | — | SNP-array uncertain-range whole-gene deletion (parse-only here) |
| `NC_000023.11:g.(?_31120496)_(33339477_?)del` | recommended | — | MLPA half-uncertain (`?`) whole-gene deletion (parse-only here) |

## `deletion.md:105-110` — mosaic and chimeric

>     - **`NC_000023.11:g.33344590_33344592=/del`**<br>

Ferro: the deletion-specific instances of the general mosaic (`=/`) and chimeric (`=//`) compact
forms — not deletion-specific rules; the `del` member inside carries the same 3'rule as a standalone
deletion (`src/hgvs/parser/variant.rs`). When that 3'rule moves the `del` member, the compact form's
shared range no longer holds, and ferro re-spells the description with the accession repeated on
each side of the `=/` rather than shifting the shared range with the deletion (see the last row).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NC_000023.11:g.33344590_33344592=/del` | recommended | — | mosaic form, spec's own accession (parse-only here) |
| `NC_000023.11:g.33344590_33344592=//del` | recommended | — | chimeric form, spec's own accession (parse-only here) |
| `NM_004006.3:c.5695_5697=/del` | recommended | self | executable mosaic twin — the `del` member is the 3'-most three of the A-stretch, already a fixed point |
| `NM_004006.3:c.5695_5697=//del` | recommended | self | executable chimeric twin — same member, preserved |
| `NM_004006.3:c.5690_5692=/del` | conformant | `NM_004006.3:c.5690_5692=/NM_004006.3:c.5695_5697del` | the `del` member 3'-shifts to `c.5695_5697` but the `=` range stays put, so ferro abandons the compact form and repeats the accession after `=/`. Re-parses, but the recommended spelling is the compact `c.5695_5697=/del` (both ranges name the same change). No tracking issue yet |

## `deletion.md:114-117` — a deletion is never written with a length suffix

> !!! note "Can I use <code class="invalid">NG_012232.1:g.123del6</code> to describe a 6 nucleotide deletion?"

Ferro: `class="invalid"` — rule 1. The repair is determinate: a length-`N` suffix starting at
position `p` becomes the range `p_(p+N-1)` (`detect_del_size_suffix`, whose fixture is this clause's
own `NG_012232.1:g.123del6`).

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[absolute-prohibition-enforcement-stage](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize.
<!-- why:END:absolute-prohibition-enforcement-stage -->

</details>

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NG_012232.1:g.123_128del` | recommended | — | the correct range form for a six-nucleotide deletion (parse-only here) |
| `NG_012232.1:g.123del6` | refused | — | the disallowed length-suffix spelling (`class="invalid"`) — strict rejects (`DelSizeSuffix`); lenient repairs to `g.123_128del` (parse-only here) |
| `NM_004006.3:c.5692_5697del` | recommended | self | executable twin of the correct form — the 3'-most six of the A-stretch (`c.5698` is `T`), so a fixed point |
| `NM_004006.3:c.5692del6` | refused | — | executable twin — strict rejects the length suffix; lenient repairs it to the range form `c.5692_5697del`, which is already 3'-most |

## `deletion.md:119-122` — why `c.3921del` and not `c.3922del`

> !!! note "In the example above, `LRG_199t1:c.3921del`, should the description based on a coding DNA reference sequence not be `LRG_199t1:c.3922del`?"

Ferro: this restates the exon/exon-junction exception adjudicated at `deletion.md:20-22` as a Q&A —
the far-side `c.3922del` is pulled back to `c.3921del` so that translating the position back does not
land in the wrong exon. Same governing record (`exon-junction-dup-converge-from-the-far-side`); ferro
is correct here on the coding axis.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `LRG_199t1:c.3922del` | recommended | — | the far-side spelling, pulled back to `c.3921del` by the clamp (parse-only here — see `deletion.md:20-22` for where the pull-back is pinned) |

## `deletion.md:124-127` — a deletion is never written by exon label

> !!! note "Is the description of a deletion of exon 17 as <code class="invalid">c.EX17del</code> still allowed?"

Ferro: `EX17` is not a legal position token on any axis — an exon label is a fact about annotation,
not a coordinate a description can carry. No conformant grammar accepts it; this is a parse
rejection, not a normalization.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:c.EX17del` | refused | — | `EX17` is not a valid position token — rejected at parse |

## `deletion.md:129-132` — BRCA1 Alu breakpoint, 3'rule reaffirmed

> !!! note "Deletions in the _BRCA1_ gene are usually mediated by Alu sequences having a very high homology, reaching 100% in the breakpoint region. In such cases, what nucleotide should be used to describe the deletion breakpoint?"

Ferro: a plain restatement of the 3'rule mechanics — shift the alignment as far 3' as possible, and
the first nucleotide that then differs is the first nucleotide deleted. No new content beyond the
3'rule already covered at `deletion.md:20-22`; it is independent textual support for the general
shuffle algorithm's correctness criterion. Descriptive — nothing deletion-specific to adjudicate,
no verdict row.

## `deletion.md:134-141` — PCR breakpoint uncertainty (Kamsteeg letter)

> !!! note "PCR analysis of a gene on the X-chromosome shows products for exons 1_3, no product is detected for exons 4_14 (exon 14 is the last exon of the gene). Since PCR fails already when one primer is not hybridising, we are not sure whether exon 4 and 14 are completely absent, or only partially. To describe the deletion I would therefore like to use the last base of exon 3 with "+?" and the last base of exon 13 with a "+?". What are your recommendations? (Erik-Jan Kamsteeg, Nijmegen, Nederland)"

Ferro: the spec's own answer ends in an open rhetorical question ("Is this really more
informative...") and does not commit to one of `c.(987+123_?)del` or `c.(987+1_?)del` as canonical —
both are legal uncertain-position deletions denoting *different breakpoint claims from different
evidence*, so there is no confluence question to resolve. Ferro accepts both; there is nothing to
enforce beyond that. Descriptive / open — no verdict row, and not a ledger-change candidate.

## `deletion.md:143-150` — CFTR deltaF508

> !!! note "In literature I often see the description "deltaF508" for a variant in the _CFTR_ gene in patients with Cystic Fibrosis. Is the variant detected in these patients <code class="invalid">NM_000492.3:c.1522_1524delTTT</code>?"

Ferro: a worked 3'rule example inside a repeat-like run (`..ATCTTTGGT..`): three literal deletions
all give `p.Phe508del`, and the 3'rule alone reduces them only to one of *two* forms
(`c.1521_1523del` or `c.1522_1524del`) — it cannot distinguish which triplet was truly altered
without external evidence. That is a statement about the limits of sequence-alone re-derivation, not
a normalization defect. The payload-bearing `delTTT` is `class="invalid"`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_000492.3:c.1521_1523del` | recommended | — | the 3'-most form, mostly the change found in patients (parse-only here) |
| `NM_000492.3:c.1522_1524delTTT` | refused | — | the payload-bearing spelling (`class="invalid"`) (parse-only here) |

## `deletion.md:152-163` — the "los"/"dec" naming suggestion

> !!! note "Suggestion to use "los" for a loss from a mono-nucleotide stretch."

Ferro: historical correspondence — the proposal to add `los`/`dec` terms alongside `del` was not
adopted; `del` stands, and no `los`/`dec` token exists in the grammar. Nothing to enforce. The
closing line — "a description should be clear/unequivocal and it is not intended to contain other
information" — is quoted in the ledger's `canonical-form-choice-when-both-legal` rationale as spec
support for re-derivation-from-sequence over provenance-preservation, but the passage itself states
no rule. Descriptive — no verdict row.
