# Duplication — ferro's reading

ferro's reading of the HGVS **duplication** recommendations on the transcript (`r.`) axis,
clause by clause — each spelling with the form ferro normalizes it to and a verdict on that
output. New here? See [How to read a page](../../reading-guide.md) for the verdicts, the table
conventions, and the recurring terms.

*DNA twin: [Duplication (`c.`/`g.`)](../DNA/duplication.md).*

Two ledger records reach this page's clauses. `inverted-duplication-is-derived-as-ins-range-inv`
cites `RNA/duplication.md:21` (and `RNA/insertion.md:18`) directly, so it carries `r.` authority on
its own molecule-native basis; `duplication-must-ranks-the-label-not-the-partition` cites the DNA
twin `DNA/duplication.md:18`, so on the `r.` axis the MUST it grades rests on the verbatim RNA twin
`RNA/duplication.md:19` — noted where it applies rather than claimed as an `r.`-jurisdiction ruling.
The `:16`, `:17`, `:22-25` and the definition/Discussion clauses carry **no Why block** — the
reading is CONFIRM-by-inspection against the spec text and the shipped code, not an adjudicated
ruling.

Executable rows use `NM_004006.3`, the one transcript in the committed slice. The spec's own
examples are written on bare illustrative sequences (`r.7dup`, `r.6_8dup`) and on `LRG_2t1`, which
the slice does not carry, so those rows are parse-only (`—`) — ferro cannot read their bases here.
The `NM_004006.3` base facts the executable rows rely on (`r.75_77` is `aaa`; `r.6_8` is `uug` with
`r.9` = `g`; `r.1470_1478` is `aca` x 3 with `r.1479` = `u`; `r.123` is `c`) are the same ones
established on `RNA/deletion.md`.

## `duplication.md:5` — definition: a copy inserted directly 3'

> Duplication: a sequence change where, compared to a reference sequence, a copy of one or more
> nucleotides is inserted **directly 3'** of the original copy of that sequence.

Ferro: a duplication inserts a copy of one or more nucleotides directly 3' of the original copy;
the "directly 3'-flanking" test is what separates a `dup` from an `ins` (see `:18-21`). The
description carries only the duplicated range, never the duplicated bases (see `:32-35`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.7dup` | recommended | self | single-nucleotide RNA duplication; the 3'-most `u` of the `r.6_7` `uu` pair, so already a fixed point. The spec's own bare example is the accession-less `r.7dup` (`:28-31`) — not table-executable, so shown here on the slice accession |

## `duplication.md:15` — DNA-level reporting is authoring practice, not a normalizer rule

> all variants **should be** described on the DNA level; descriptions on the RNA and/or protein
> level may be given in addition.

Ferro: this is guidance to authors about which level(s) to report, not a constraint on how a given
`r.` description, once written, is normalized. An `r.` description stands on its own.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.7dup` | recommended | self | a lone RNA-level description is fully valid; ferro does not require an accompanying `c.`/`g.` form |

## `duplication.md:16` — the range must name two different positions

> `positions_duplicated` should contain **two different positions**, e.g., `123_126`, not `123_123`.

Ferro: rule 2 (lowercase "should"); a same-position range is repaired to the single-position form.
The `W4003 SinglePositionRange` machinery is axis-general, so it handles an `r.` range identically
to its `c.` twin.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123dup` | recommended | self | the correct single-position form (`r.123` is `c`, with no `c` run to shift along) |
| `NM_004006.3:r.123_123dup` | refused | — | same-position range — rejected at parse (`W4003 SinglePositionRange`) in strict mode; lenient repairs to `r.123dup` (`error_mode_tests::w4003_single_position_range_emission` pins the `c.` twin, the machinery is axis-general) |

## `duplication.md:17` — the range is listed 5' to 3'

> the `positions_duplicated` should be listed from **5' to 3'**, e.g., `123_126`, not `126_123`.

Ferro: a reversed, non-circular range has no exception on the `r.` axis (the `<high>_<low>` form is
admitted only for `o.`/`m.` circular references), so it is refused rather than reordered.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.123_126dup` | recommended | self | correctly ordered range; no 3' tandem copy to shift into, so it is already 3'-most |
| `NM_004006.3:r.126_123dup` | refused | — | reversed range on a non-circular reference — rejected at parse |

## `duplication.md:18-21` — tandem-only; MUST be a dup not an ins; an inverted duplication is an insertion

> - by definition, duplication may only be used when the additional copy is **directly 3'-flanking**
>   of the original copy (a "tandem duplication").
>     - when a variant can be described as a duplication, it **must** be described as a duplication
>       and not as e.g., an insertion (see [Prioritization](../general.md)).
>     - when there is no evidence that the extra copy of a sequence detected is in tandem (directly
>       3'-flanking the original copy), the change can not be described as a duplication, it should
>       be described as **an insertion** (see [Insertion](insertion.md)).
>     - **inverted duplications** are described as an insertion (`r.234_235ins123_234inv`), not as a
>       duplication (see [Inversion](inversion.md)).

Ferro: four points. (i) a `dup` is only for a copy directly 3'-flanking the original; (ii) a
duplicating insertion **must** be relabeled `dup` — ferro collapses `r.X_Yins<5'-flank>` to `dup`
through the shared canonical-form path (`rules::insertion_is_duplication`, reached by
`normalize_rna`); (iii) a copy that is *not* directly 3'-flanking is an insertion, not a `dup`
(worked at `:47-52`); (iv) an inverted duplication is spelled `ins<range>inv`, never `dup` — the
`dupinv` shorthand is refused at parse.

<details class="ss-why"><summary>Why ferro reads it this way</summary>

<!-- why:START -->
> **[duplication-must-ranks-the-label-not-the-partition](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — The rule that a duplication must be labelled 'dup' ranks the label of each piece ferro derives, not the partition; the one exception is a net-longer tandem copy of a multi-base motif, where the derivation is cut to expose the dup rather than merged into a delins.
>
> **[inverted-duplication-is-derived-as-ins-range-inv](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md)** — An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum.
<!-- why:END:duplication-must-ranks-the-label-not-the-partition,inverted-duplication-is-derived-as-ins-range-inv -->

</details>

`duplication-must-ranks-the-label-not-the-partition` cites `DNA/duplication.md:18`; on the `r.`
axis the same MUST is carried verbatim by `RNA/duplication.md:19`, and a `DNA/` citation cannot
scope `r.` on its own — the record is surfaced here for its mechanism (the MUST ranks a *label*,
not a partition), with the `r.`-jurisdiction authority being the RNA twin.

**The inverted-duplication form is not yet derived on `r.`** The `ins<range>inv` re-spell is wired
only in `normalize_genome`; `normalize_rna` has no equivalent, so ferro **preserves** an already
spelled `r.234_235ins123_234inv` input but does not **mint** it from a reverse-complement literal
payload on the `r.` axis. The ledger record calls this its "largest known gap" and names #1946's
render stage as the intended fix. An inverted-dup input whose payload is spelled as literal
reverse-complement bases therefore stays literal on `r.` — valid HGVS that re-parses, so the
verdict is `conformant`, not `bug`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.6_7insu` | recommended | `NM_004006.3:r.7dup` | the spec example-1 NOTE: the duplicating insertion `r.6_7insu` copies the `u` at `r.6` directly 3', so `:19`'s MUST relabels it to the `dup`, landing on the 3'-most `u` (`r.7`). Mirrors `issue_736::rna_insertion_with_u_collapses_to_dup` |
| `NM_004006.3:r.234_235ins123_234inv` | conformant | `NM_004006.3:r.234_235insguugacauuguucagggcaugaacucuuguggauccuuuuucuuuuggcaguuuuugcccugucaggccuucgaggaggucuaggaggcgccucccauccuguaggucacug` | the spec's own inverted-duplication spelling — an `ins<range>inv`. On a real reference ferro **expands** the range payload to literal reverse-complement bases via `try_expand_rna_ins` rather than preserving the `ins<range>inv` spelling: the derivation that mints/keeps the range-inv form is wired only in `normalize_genome`, so the `r.` axis emits literals (`inverted-duplication-is-derived-as-ins-range-inv`; #1946 names the render-stage fix). Valid HGVS that re-parses and denotes the same sequence, so conformant, not the recommended spelling |
| `NM_004006.3:r.123_456dupinv` | refused | — | the `dupinv` shorthand — rejected at parse ("Unexpected trailing characters: 'inv'"); `:21` requires the `ins<range>inv` form. Pinned by `reject_dupins`/`recommended_form_pins` |

## `duplication.md:22-25` — the 3'rule, and the exon/exon-junction exception that does not apply here

> - for all descriptions, the **most 3' position** possible of the reference sequence is arbitrarily
>   assigned to have been changed (**3'rule**).
>     - the 3'rule also applies for changes in single residue stretches and tandem repeats.
>     - **NOTE**: the exception to the 3'rule for duplications around exon/exon junctions (see
>       [Duplications](../DNA/duplication.md)) does not apply when describing variants based on an
>       RNA reference sequence.

Ferro: the 3'rule itself applies on `r.` exactly as elsewhere, including across single-residue runs
and tandem repeats — the two executable rows below show each. The NOTE is a second, separate
instruction — it **switches off** the DNA exon/exon-junction 3'rule exception on the RNA axis, so
an `r.` duplication must still shift to its 3'-most position even when doing so crosses an intron.
This is the twin of `RNA/deletion.md:20`.

**Ferro currently gets the second half wrong.** `normalize_rna` (#334) applies the DNA exon/exon
clamp to `r.` del/dup edits — it never lets a duplication 3'-shift across an exon/exon junction on
the RNA axis — which is exactly what this NOTE says must not happen. No ledger record governs this
clause; the nearest record, `exon-junction-dup-converge-from-the-far-side`, decides only the
far-side pull-back on `c.`/`n.` and rests its own `r.` reach on `general.md:43`, which lists `r.`
as in-scope and is in tension with this NOTE's flat exclusion of the RNA axis. Filed as
[#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211). The output is valid HGVS — it
is the wrong *preference*, not a malformed string — so the verdict is `conformant`, not `bug`.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.76dup` | recommended | `NM_004006.3:r.77dup` | 3'rule over a single-residue stretch: `r.75_77` is `aaa`, so the duplicated `a` lands on `r.77` |
| `NM_004006.3:r.1470_1472dup` | recommended | `NM_004006.3:r.1476_1478dup` | 3'rule over a tandem repeat: `r.1470_1478` is `aca` x 3 (`r.1479` is `u`), so the duplicated unit lands on the 3'-most copy |
| `LRG_2t1:r.1034_1036dup` | recommended | — | illustrative junction-crossing form (the exon/exon boundary sits between `r.1035` and `r.1036`, as in `RNA/deletion.md`'s `LRG_2t1` example); on `r.` the 3'rule must cross it. No committed `LRG_2t1` fixture, so parse-only here |
| `LRG_2t1:r.1033_1035dup` | conformant | — | on `r.` the exon/exon-junction exception should not apply, so the 3'rule requires crossing to `r.1034_1036dup`. Ferro's `normalize_rna` wrongly clamps at the junction and keeps this spelling as a fixed point ([#2211](https://github.com/fulcrumgenomics/ferro-hgvs/issues/2211)); recommended is the 3'-most-crossing form. Code-derived: no committed `LRG_2t1` fixture exercises this row, and the `c.`/`n.` sibling clamp is pinned in `tests/it/issue_334_exon_junction_exception.rs`, which carries no `r.` del/dup case |

No executable `NM_004006.3` row can show the wrong behavior: the committed slice is built as a
single flat exon (`tests/fixtures/shadow_spec/build_fixture.py`), so no `r.` duplication on it can
reach an exon/exon junction. The clamp is enforced by code and, so far, by no test.

**Non-confluence.** The pair `{r.1033_1035dup, r.1034_1036dup}` describes one variant but, under the
current defect, does not converge: ferro treats `r.1033_1035dup` as its own fixed point rather than
normalizing it to the recommended `r.1034_1036dup`.

See also → `RNA/deletion.md:18-20` (the identical NOTE for deletions, same defect shape, #2211).

## `duplication.md:28-31` — one nucleotide

> - **`r.7dup` (one nucleotide)**<br>
>   the duplication of a `u` at position `r.7` in the sequence `acuuacu`<code class="ins">u</code>`gcc`.<br>
>   **NOTE**: it is **not** allowed to describe the variant as <code class="invalid">r.6_7insu</code>
>   (see [prioritisation](../general.md)).

Ferro: the spec's own single-nucleotide worked example, plus the NOTE that the duplicating
insertion `r.6_7insu` must be spelled `dup` (the MUST of `:19`).

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.7dup` | recommended | self | executable twin of the spec's bare `r.7dup`: `r.7` is the 3'-most `u` of the `r.6_7` `uu` pair, so it is a fixed point (the accession-less spec spelling is not table-executable) |
| `NM_004006.3:r.6_7insu` | recommended | `NM_004006.3:r.7dup` | the NOTE's forbidden `insu` spelling — relabeled to the recommended `dup` (see `:18-21`); a duplicating insertion is written as a duplication |

## `duplication.md:32-35` — several nucleotides, and no payload-bearing spelling

> - **`r.6_8dup` (several nucleotides)**<br>
>   a duplication from position `r.6` to `r.8` in the sequence `acaauugc`<code class="ins">ugc</code>`c`.<br>
>   **NOTE**: the recommendation is not to describe the variant as
>   <code class="invalid">r.6_8dupugc</code>, i.e., describe the deleted nucleotide sequence.

Ferro: a plain multi-nucleotide duplication carries only the range, never the duplicated payload
sequence. `canonicalize_edit` drops an explicit `dup` sequence/length even without a reference.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.6_8dup` | recommended | self | executable twin of the spec's bare `r.6_8dup`: `r.6_8` is `uug`, with `r.9` = `g` and no 3' tandem copy, so it is already 3'-most (the accession-less spec spelling is not table-executable) |
| `NM_004006.3:r.6_8dupugc` | refused | — | the disallowed payload-bearing spelling — `class="invalid"`; rejected at parse in strict mode; lenient drops the payload and normalizes to `r.6_8dup` (`issue_1615_denoted_sequence_oracle.rs` pins the analogous `NR_TEST.1:r.6_8dupugc` → `r.6_8dup`) |

## `duplication.md:38-46` — Discussion: why not describe a duplication as an insertion

> !!! note "Why do we not describe a duplication as an insertion?"
>
>     Although duplications are basically a special type of insertion, there are several reasons why
>     the recommendation is to describe duplications separately.
>
>     - the description is simple and shorter;
>     - it is clear and prevents confusion regarding the position when an insertion is incorrectly
>       reported, like <code class="invalid">c.22insG</code>;
>     - insertion more or less means "coming from elsewhere".
>       Mechanistically, a duplication is most likely caused by a local event, DNA polymerase
>       slippage, duplicating a local sequence.

Ferro: rationale prose for the dup-over-ins preference. It states no rule beyond the `:19` MUST
already covered at `:18-21` (a duplicating insertion is written as a `dup`); nothing here is a
constraint the normalizer applies separately.

## `duplication.md:47-51` — Discussion: a non-tandem duplicative copy is an insertion

> !!! note "How should I describe the change `aucg`<code class="spot1">aucgaucgaucg</code><code class="spot2">a</code>`ggguccc` to `aucg`<code class="spot1">aucgaucgaucg</code><code class="spot2">a</code><code class="ins">aucgaucgaucg</code>`ggguccc`? The fact that the inserted sequence (<code class="ins">aucgaucgaucg</code>) is present in the original sequence, suggests it derives from a duplicative event."
>
>     The variant should be described as an insertion; `r.17_18ins5_16`.
>     A description using "dup" is not correct since, by definition, a duplication should be
>     **directly 3'-flanking of the original copy** (in tandem).
>     Note that the description given still makes it clear that the sequence inserted between `r.17`
>     and `r.18` is probably derived from nearby, i.e. positions `r.5` to `r.16`, and thus likely
>     derived from a duplicative event.

Ferro: the worked case of `:20` — a copy present in the reference but **not** directly 3'-flanking
is an **insertion**, not a duplication, even when it plainly derives from a duplicative event. The
classification is what the spec fixes here; the by-range payload (`ins5_16`) versus a literal
payload is a representation choice the spec leaves open (`insertion.md:20` admits both), governed
by `canonical-form-choice-when-both-legal`. Ferro will not relabel a non-3'-flanking copy as `dup`
(`rules::insertion_is_duplication` checks only the immediate 5' flank) and expands range payloads
to literal bases on `r.` (`try_expand_rna_ins`) — both conformant.

| Input | Verdict | Normalizes to | Notes |
|---|---|---|---|
| `NM_004006.3:r.17_18ins5_16` | conformant | `NM_004006.3:r.17_18insuuuggugggaag` | executable twin of the spec's own answer `r.17_18ins5_16` (accession-less, not table-executable) — an insertion naming its likely source range, **not** a `dup`, because the copy is not directly 3'-flanking. ferro will not relabel a non-3'-flanking copy as `dup`; it does expand the range payload to literal bases via `try_expand_rna_ins`, a representation choice left open by the spec (`insertion.md:20` admits both), governed by `canonical-form-choice-when-both-legal` — valid HGVS that re-parses, so conformant |

See also → `RNA/deletion.md:18-20` (the twin exon/exon-junction NOTE and the #2211 defect shape),
`RNA/insertion.md` (the `ins`/`ins<range>inv` forms a non-tandem or inverted copy takes),
`general.md` (the prioritisation the `:19` MUST points at).
