# Repeated Sequences — ferro's reading

ferro's reading of the HGVS **repeated sequence** recommendations on the transcript (`r.`)
axis, clause by clause — each spelling with the form ferro normalizes it to and a verdict on
that output. New here? See [How to read a page](../../reading-guide.md) for the verdicts, the
table conventions, and the recurring terms.

*DNA twin: [Repeated Sequences (`c.`/`g.`)](../DNA/repeated.md).*

Three ledger records reach this page's clauses. `rna-repeat-range-plus-unit-redundancy` is the one
record that governs the central `:22`/`:27` conflict, and it is **undecided** on purpose — its
summary states the conflict rather than resolving it. `confluence-gate-is-apply-equality-on-every-determined-axis`
cites `RNA/repeated.md:20-21` directly, so it carries `r.` authority on its own molecule-native
basis and governs the sequenced-versus-interruptible split at `:19-21`.
`canonical-form-choice-when-both-legal` is molecule-neutral (a general mechanism) and reaches the
range-only-versus-start-plus-unit choice at `:31-34` unchanged. The definition, the Community
Consultation NOTE, the DNA-level and size Notes, the multiple-of-3 restriction with its UTR
carve-out, the composite-repeat syntax, and the FMR1/HD 3'rule Examples carry **no Why block** —
the reading is CONFIRM-by-inspection (or, for `:22`/`:27`, an open finding) against the spec text
and the shipped code, not an adjudicated ruling.

The spec's own examples all sit on a foreign accession the committed slice does not carry —
`NM_024312.4` (the coding-restriction examples), `NM_002024.5` (FMR1), `NM_002111.6` (HD) — or on a
bare, accession-less spec spelling (`r.-124_-123[14]`, `r.123cug[23]`). None of those is executable
here, so every row shown on a spec accession is parse-only (`—`), and an accession-less spelling is
discussed in prose rather than forced into a verdict table because it does not parse. The one slice
transcript is `NM_004006.3`, and it does carry a real repeat tract the spec never names: an `ac`
dinucleotide run of 8 copies at `r.*460_*475` in the 3'UTR. That anchors the single **executable**
row on this page (`:24-27`), where ferro's `Normalizes to` is its measured output on the committed
slice rather than a `—`.

The grammar fact under the whole page: `docs/syntax.yaml`'s `rna.rpt` production admits **exactly
two** canonical repeat spellings — a **range only** (`r.9495_9497[4]`, no unit) or a **single start
position plus a unit** (`r.9495caa[4]`, no range). A description carrying *both* a range and a unit
(`r.-125_-123cug[4]`, `r.-6_-3g[6]`) matches neither production. That backbone is what makes `:22`
and `:27` a genuine self-contradiction rather than a matter of taste, and it is the independent
ground the `:22`/`:27` section leans on.

## `repeated.md:5` — definition: a repeat unit present several times

> Repeated sequence: a sequence where, compared to a reference sequence, a segment of **one or more** nucleotides (the repeat unit) is present several times, one after the other.

Ferro: the definition of the object. A repeat is a unit (one or more nucleotides) tandemly present
in several copies. It states no normalization obligation, so there is no verdict row.

## `repeated.md:9` — Community Consultation: entire-range-only is a *proposal*, not yet a rule

> **NOTE**: a Community Consultation proposal is being prepared which will suggest to allow only the format where the **entire range** of the repeated sequence is indicated; so `r.123_191cag[23]`, **not** `r.123cag[23]`.

Ferro: a forward-looking NOTE describing a proposal that is being *prepared*, not a current rule.
It is the inverse of the page's own worked Examples, which use single-anchor forms
(`r.53agc[19]`, `r.-124ug[14]`), so a parser must **not** enforce entire-range-only today. Ferro
accepts single-anchor start-plus-unit forms, which is correct against the current text. A
"is being prepared" / "will suggest" NOTE is a suggestion whose disposition has to be checked
before it is treated as guidance, so no obligation is adjudicated from it here.

## `repeated.md:17-18` — DNA-level reporting, and the size range of repeats

> - all variants **should be** described on the DNA level; descriptions on the RNA and/or protein level may be given in addition.
> - repeated sequences include both small (mono-, di-, tri-, etc., nucleotide) and larger (kilobase-sized) repeats.

Ferro: both are reporting-practice statements, not constraints on how a written `r.` description is
normalized. The DNA-level recommendation ("may be given in addition") does not forbid a lone `r.`
repeat — ferro normalizes one without demanding a `c.`/`g.` companion — and the size range is
descriptive scope. This is the twin of `RNA/inversion.md:15`. That RNA descriptions are secondary
to and derived from the DNA level is load-bearing below: it is one reason the DNA page's parallel
treatment weighs in when the RNA page contradicts itself at `:22`/`:27`. No verdict row.

## `repeated.md:19-21` — repeat-position form preferred; the sequenced-versus-interruptible split

> - the format based on **repeat position** is preferred, descriptions of the repeat sequence quickly become too lengthy.<br>
>   **NOTE**: while `r.123cug[23]` describes a repeat of 23 `cug` units, `r.123_125[23]` describes a tri-nucleotide repeat of 23 units which **could be interrupted** with other units (e.g., a rare `cua`).
>   The description `r.123cug[23]` can thus only be used when the repeat was sequenced.

Ferro: the range-only form (`r.123_125[23]`) and the start-plus-unit form (`r.123cug[23]`) make
different *epistemic* claims — the range-only admits interruptions, the unit-spelled asserts a pure
tract "only... when the repeat was sequenced" — yet on a pure tract they denote the same sequence.
Ferro's confluence gate is apply-equality on the denoted bases, so it deliberately treats the two
as one equivalence class and encodes no rung for the provenance difference; encoding it would
re-introduce spelling-dependence the gate exists to remove. `:21` grounds the distinction in
whether the tract "was sequenced", which is provenance, not sequence, and is not recoverable from
the bases.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[confluence-gate-is-apply-equality-on-every-determined-axis](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Ferro's release gate asserts that inputs denoting the same sequence on every determined axis, the protein axis excluded, normalize to one output; it is asserted over decided equivalence classes only, and equivalence is judged by applying the descriptions to the reference rather than by comparing normalized strings.
<!-- why:END:confluence-gate-is-apply-equality-on-every-determined-axis -->

</details>

The two forms are accession-less in the spec, so there is no executable verdict row; the record's
own summary names `RNA/repeated.md:20-21` as the clause it does **not** try to collapse — the
epistemic gap is a canonical-form question, not a confluence one. The concrete instance of this
split appears at the HD Example (`:49-56`), where `r.53agc[19]` (pure, needs sequencing) and
`r.53_55[31]` (range, admits interruptions) describe the same locus.

## `repeated.md:22` — range-plus-unit is called invalid, and `:27` publishes it as valid

This is the page's central conflict, and it is **upstream's**, not ferro's. `:22` states a rule:

> - the format <code class="invalid">r.-125_-123cug[4]</code> should not be used; it contains redundant information (`-125_-123` and `cug`).

Yet four lines later, `repeated.md:27` publishes exactly that range-plus-unit shape as a valid
worked example:

>   As such, `NM_024312.4:r.-6_-3g[6]` is valid as the reading frame is not affected.

Ferro's reading — set out because the governing record is deliberately **undecided**, and because
the reading matters more than the verdict here:

1. **Both strings are the same ill-formed shape: a range plus a unit.** Neither matches either
   `rna.rpt` production (see the grammar fact in the intro). `:22` names that shape and gives its
   ground — the range already fixes the unit — which agrees with the grammar; `:22` is a **rule
   with a stated reason**. `:27` is an **example** that violates both the grammar and `:22`.

2. **`:27`'s instance is not merely redundant — read literally it is self-contradictory.** In
   `:22`'s example the range length equals the unit length (`-125_-123` is 3 nt, `cug` is 3 nt), so
   the two pieces are consistent-but-redundant. In `:27` the range `-6_-3` spans **4 positions**
   while the unit `g` is **1 nt** repeated **6×** — the range says a 4-nt unit, the `g[6]` says a
   1-nt unit in six copies, and these cannot both be true. `:27`'s string is worse than the shape
   `:22` forbids.

3. **What `:27` is trying to demonstrate is sound and untouched by this.** Its job is the UTR
   carve-out to the multiple-of-3 restriction (`:24-27`): a mononucleotide repeat is legal in the
   5'UTR because the reading frame is not affected. That carve-out is correct and implemented (see
   `:24-27` below); only the *spelling* is wrong. The grammar-conformant spelling of the intended
   variant is the start-plus-unit `r.-6g[6]` — a single start position does not fix the unit, so it
   is **not** barred by `:22`.

4. **The DNA sibling does not corroborate `:27`.** `DNA/repeated.md:24` carries the byte-parallel
   `c.-6_-3G[6]` with no redundancy rule of its own, but it is plainly the same upstream text
   copied — the same error propagated, not a second witness.

So the campaign position leans toward `:22` (the rule generalises; the example does not), with the
canonical valid form being the start-plus-unit `r.-6g[6]` — but it is held **open pending an
upstream filing**, since the conflict is genuinely upstream's to correct (the ledger notes it as
owed upstream, `#466`). The record stays undecided.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[rna-repeat-range-plus-unit-redundancy](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — Whether an RNA repeat may combine a position range with its own repeat unit is undecided — the recommendations both call that shape invalid as redundant and publish it as a worked example — and ferro's current lenient-mode agreement with the published form is an incidental effect of tract maximization rather than a resolution.
<!-- why:END:rna-repeat-range-plus-unit-redundancy -->

</details>

**Ferro's behaviour is what makes the decision non-free, and it hides an open finding.** Since the
`#1631` lenient repair (`tests/it/issue_1631_repeat_label_repair.rs`), a repair hands the
normalizer the anchored `r.-6g[6]` — and the normalizer's tract maximization then **re-widens** the
single-nucleotide anchor back into a range, re-emitting `r.-6_-3g[6]`: the very `:22`-invalid,
`:27`-contradictory string. So ferro's agreement with `:27` is an artifact of two independent
passes, exactly as the record's summary says, not a resolution. If the ledger ever decides for
`:22`, ferro's lenient/normalizer output on this input becomes non-conformant, and the fix is to
stop the repeat renderer widening a single-nucleotide anchor into a range. That widening is a real
open finding tied to this conflict; **no issue is filed** for it — flagged for the operator, not
filed here.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_024312.4:r.-6g[6]` | conformant | — | the grammar-valid start-plus-unit anchored form (a single start position does not fix the unit, so it is **not** barred by `:22`) — the intended variant's conformant spelling. The `#1631` lenient repair emits this, but the normalizer's tract maximization re-widens it to the range-plus-unit `r.-6_-3g[6]`. Parse-only (foreign accession, outside the slice); the widening is the open finding above, no issue filed |
| `NM_024312.4:r.-6_-3g[6]` | conformant | — | `:27`'s published-valid range-plus-unit shape, and the string ferro's normalizer manufactures by widening the row above. Ferro round-trips it, so it re-parses as valid HGVS; whether it is the *right* form is what the undecided record leaves open — and read literally it is internally contradictory (4-position range against a 1-nt unit ×6). Parse-only |
| `NM_024312.4:r.-125_-123cug[4]` | conformant | — | `:22`'s own `class="invalid"` range-plus-unit example, shown on the page's foreign accession (`:22` writes it bare). Here range length equals unit length (3 nt each), so it is redundant rather than contradictory — the milder of the two `:22`-shape instances. Same undecided treatment; the grammar-conformant re-spellings are `r.-125cug[4]` (start-plus-unit) or `r.-125_-123[4]` (range-only). Parse-only |

## `repeated.md:23` — composite repeats: successive per-unit listing

> - for **composite repeats**, the basic format can be used, successively listing each different repeat unit; <code class="invalid">r.456_465[4]466_489[9]490_499[3]</code>.

Ferro: descriptive syntax for composite (mixed-unit) repeats. The example is wrapped
`<code class="invalid">` in the source — it illustrates the *shape*, not a valid string, and it uses
a range-per-unit listing where the DNA page's mixed-repeat canonical form is unit-per-segment. That
page-internal inconsistency is noted but out of scope; no ferro obligation is adjudicated here, so
there is no verdict row.

## `repeated.md:24-27` — the multiple-of-3 coding restriction, and the UTR carve-out

> - **exception**: using a coding RNA reference sequence, a repeated sequence variant description can be used only for repeat units with a length which is a multiple of 3, i.e. which can not affect the reading frame.
>   Consequently, use `NM_024312.4:r.2692_2693dup` and **not** <code class="invalid">NM_024312.4:r.2686a[10]</code>; use `NM_024312.4:r.1741_1742insuauauaua` and **not** <code class="invalid">NM_024312.4:r.1738ua[6]</code>.
>   This restriction only applies to the coding sequence, which does not include the UTR sequence.
>   As such, `NM_024312.4:r.-6_-3g[6]` is valid as the reading frame is not affected.

Ferro: on a coding RNA reference a repeat unit must be a multiple of 3 (so it cannot shift the
reading frame); an out-of-frame coding repeat is repaired away to a `dup`/`ins` describing the same
change. The restriction is gated off in the UTR, where any unit length is legal. Both halves are
implemented: the `unit_len % 3` gate lives in the repeat rules (`normalize_repeat`,
`deletion_to_repeat`, `duplication_to_repeat` in `rules.rs`), and the UTR exemption via
`span_is_coding` (`src/normalize/mod.rs`); the worked repairs are pinned in
`tests/fixtures/spec-worked-examples/cases.json`. There is **no dedicated ledger record** — the
behaviour is a straight reading of the text, not an adjudicated conflict — so no Why block.

The carve-out verdict does **not** depend on `:27`'s bad string: the mononucleotide-in-5'UTR
exemption is correct, and only the range-plus-unit spelling it is shown with is what `:22` disputes
(see `:22` above).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_024312.4:r.2686a[10]` | recommended | — | `:25`'s `class="invalid"` coding mononucleotide repeat (unit length 1, not a multiple of 3); ferro repairs it away from the barred repeat form to the spec's own recommended `r.2692_2693dup`, the only legal form in the CDS (all modes, `cases.json`). Parse-only (foreign accession) |
| `NM_024312.4:r.2692_2693dup` | recommended | — | the spec's own recommended rewrite of the row above — a `dup` describing the same change in-frame |
| `NM_024312.4:r.1738ua[6]` | recommended | — | `:25`'s second `class="invalid"` coding repeat (di-nucleotide unit, length 2); ferro repairs it to the spec's own recommended `r.1741_1742insuauauaua` in default/lenient (strict declines the barred spelling), per `cases.json` |
| `NM_024312.4:r.1741_1742insuauauaua` | recommended | — | the spec's recommended rewrite of the row above — an `ins` describing the same change |
| `NM_024312.4:r.-6_-3g[6]` | conformant | — | the UTR carve-out itself: a 5'UTR mononucleotide repeat is exempt from the multiple-of-3 rule (reading frame not affected), so `:27` calls it valid. This is **also** the `:22`-shape string the section above disputes — the carve-out is correct, the range-plus-unit *spelling* is the open conflict. Parse-only |
| `NM_004006.3:r.*460ac[9]` | recommended | `NM_004006.3:r.*475_*476dup` | the one executable row on this page: a real 3'UTR `ac` dinucleotide tract on the slice accession. `NM_004006.3` carries an `ac`×8 repeat at `r.*460_*475` (3'UTR, `cds_end` = 11295), so `r.*460ac[9]` asserts one extra `ac` copy. Unit length 2 is **not** a multiple of 3, yet the repeat is legal because `*460` sits in the 3'UTR — the carve-out this section states. Ferro reads the single-copy expansion as a duplication and applies the 3'rule, emitting the 3'-most `r.*475_*476dup` (the same output `r.*460_*461dup` and `r.*474_*475dup` both converge on — one variant, 3'-shifted). Measured on the committed slice |

## `repeated.md:31-34` — Example #1: `r.-124_-123[14]`, and the population-preference NOTE

> - **`r.-124_-123[14]` (alternatively `r.-124ug[14]`)**<br>
>   a repeated di-nucleotide sequence, with the first unit located from position `r.-124` to `r.-123`, is present in 14 copies.<br>
>   **NOTE**: when the repeat is variable in the population and the reference sequence has 15 units, the description `r.-123ug[14]` is preferred over `r.-97_-96del`.<br>
>   **NOTE**: when the repeat is variable in the population and the reference sequence has 15 units, the description `r.-123ug[17]` is preferred over `r.-99_-96dup`.

Ferro: two observations. First, the range-only `r.-124_-123[14]` (2-nt unit) and the start-plus-unit
`r.-124ug[14]` are the two canonical spellings the grammar admits, offered as equals — and this pair
is *legitimately* both-legal (range length 2 equals unit length 2), unlike the `:22`/`:27` hybrid.
Which of the two ferro emits is the representation choice `canonical-form-choice-when-both-legal`
governs: ferro derives the form from the resulting sequence rather than preserving the input's
spelling.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[canonical-form-choice-when-both-legal](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — When two descriptions of one variant are both legal and no clause chooses between them, ferro derives the form from the resulting sequence rather than preserving the input's spelling.
<!-- why:END:canonical-form-choice-when-both-legal -->

</details>

Second, the "preferred when variable in the population" NOTEs are **provenance-conditioned**: which
description is preferred depends on the repeat being variable in the population and on the reference
having 15 units — facts about population frequency, not about the sequence in front of the
normalizer. Ferro cannot recover them, so it cannot honour the preference. That is a **disclosed
GAP** (upstream/unimplementable), not a defect: the spec asks for a choice ferro has no input to
make. Both spec forms are accession-less, so there is no executable verdict row; the both-legal
choice and the GAP are the adjudicated content.

## `repeated.md:36-37` — Example #2: the two-allele form `r.-124_-123[14];[18]`

> - **`r.-124_-123[14];[18]` (alternatively `r.-124ug[14];[18]`)**<br>
>   a repeated di-nucleotide sequence, with the first unit located from position `r.-124` to `r.-123`, is present in 14 copies on one allele and 18 copies on the other allele.

Ferro: the same range-only-versus-start-plus-unit both-legal choice from Example #1, carried onto a
`;`-separated two-allele description. The canonical-form choice is governed identically
(`canonical-form-choice-when-both-legal`, cited at `:31-34`); nothing new is adjudicated. Both
forms are accession-less, so there is no verdict row.

## `repeated.md:39-47` — Example #3: the FMR1 `ggc` repeat and the 3'rule

> - _FMR1_ `GGC`-repeat: in literature, the Fragile-X tri-nucleotide repeat is known as the `CGG`-repeat.
>   However, based on a coding RNA reference sequence (GenBank `NM_002024.5`) and applying the **3'rule**, on the RNA level, the repeat has to be described as a `ggc`-repeat (see [Recommendations](../general.md)).
>     - **`r.-128_-126[79]`**<br>
>       an extended repeat of exactly 79 units.<br>
>       **NOTE**: `r.-128ggc[79]` can only be used when the repeat has been sequenced, excluding it is interrupted by one or more `gga`-triplets.
>
>     - **`r.-128_-126[(600_800)]`**<br>
>       the repeated tri-nucleotide sequence, starting at position `c.-128`, has an estimated size of between 600 and 800 copies.<br>
>       **NOTE**: the repeat can be pure or a mix of `ggc` and `gga` triplets.

Ferro: the `ggc` unit (not `cgg`) is a consequence of the 3'rule — the most 3' equivalent placement
of the tract — consistent with `general.md:40-42` and with `DNA/repeated.md:57-58`. The extended
and estimated-size forms below are descriptive syntax. The `r.-128ggc[79]` NOTE is the same
epistemic point as `:19-21` — the unit-spelled form can only be used when the tract was sequenced,
excluding a `gga`-triplet interruption. No independent obligation is adjudicated, so these rows are
CONFIRM-by-inspection.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_002024.5:r.-128_-126[79]` | recommended | — | the spec's own range-only form for an extended repeat of exactly 79 units, on its coding accession. The `ggc` unit follows from the 3'rule (`general.md:40-42`). Parse-only (foreign accession, outside the slice) |
| `NM_002024.5:r.-128_-126[(600_800)]` | recommended | — | the uncertain-count form `[(600_800)]` for an estimated 600–800 copies — descriptive syntax, preserved as written. Parse-only |

## `repeated.md:49-56` — Example #4: the HD `agc` repeat, and the epistemic pair made concrete

> - HD `AGC`-repeat: based on the _HTT_ (huntingtin) coding DNA reference sequence (GenBank `NM_002111.6`), applying the **3'rule**, on the RNA level, the Huntington's Disease tri-nucleotide repeat is described as an `agc` (not `cag`) repeat.
>     - **`r.53agc[19]`**<br>
>       **NOTE**: the coding RNA reference sequence (`NM_002111.6`) contains an allele of 21 `agc` repeats.<br>
>       **NOTE**: on protein level, the reference allele contains 21 `Gln`s, described as `p.Gln[21]` (alternatively `p.Q[21]`).
>       The difference derives from the fact that the `agc` repeat is interrupted by a `aac`-triplet (`caa` coding) at position 20.
>
>     - **`r.53_55[31]`**<br>
>       the coding RNA reference sequence (`NM_002111.6`) contains a tri-nucleotide allele of 32 repeats (`agc`-19, `aac`, `agc`, `cgc`, `cac`, `cgc`-7, `cuc`-2) encoding 21 `Gln` and 11 `Pro`-residues.

Ferro: the 3'rule gives `agc` (not `cag`), the same shift consequence as the FMR1 example. This
locus is where the sequenced-versus-interruptible split of `:19-21` becomes concrete: `r.53agc[19]`
(the pure, unit-spelled form, usable only when sequenced) and `r.53_55[31]` (the range-only form,
which admits the `aac` interruption at position 20) describe the same tract as different epistemic
claims — apply-equal, and one class to ferro
(`confluence-gate-is-apply-equality-on-every-determined-axis`, cited at `:19-21`). Faithful; no
independent obligation missed.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_002111.6:r.53agc[19]` | recommended | — | the spec's start-plus-unit form for the pure `agc` tract, on its coding accession — usable only when the tract is sequenced (it asserts no `aac` interruption). Parse-only (foreign accession, outside the slice) |
| `NM_002111.6:r.53_55[31]` | recommended | — | the spec's range-only form for the same locus (a tri-nucleotide allele of 32 repeats, interrupted: `agc`-19, `aac`, `agc`, `cgc`, `cac`, `cgc`-7, `cuc`-2), which admits the interruptions the unit-spelled form excludes. Apply-equal to the row above on a pure tract; the epistemic gap is a canonical-form question, not a rung. Parse-only |

See also → `RNA/insertion.md:20-21` (the by-range-versus-literal payload choice, the same
`canonical-form-choice-when-both-legal` mechanism), `general.md:40-42` (the 3'rule the FMR1/HD
examples apply), `DNA/repeated.md:57-58` (the parallel `ggc`/`agc` 3'rule treatment on the DNA
axis, and — at `DNA/repeated.md:24` — the copied `c.-6_-3G[6]` that mirrors `:27`'s conflict).
