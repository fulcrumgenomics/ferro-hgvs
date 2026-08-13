<!--
GENERATED FILE — do not edit by hand.

Rendered from tests/fixtures/grammar/normalization_stage_audit.json and the
adjudication ledger, by tests/it/normalization_stage_audit_doc.rs. Edit the
fixture, then regenerate:

    BLESS_STAGE_AUDIT_DOC=1 cargo nextest run --features dev --test it \
      -E 'test(normalization_stage_audit_doc)'
-->

# What governs each normalization stage

Ferro's normalizer makes a decision at every stage of its pipeline. Some of those decisions are
backed by a clause of the HGVS recommendations, some by an adjudication record, and some by
nothing more than what the implementation happens to do. This document is the inventory of
which is which, so that a deliberate choice can be told from an accident — by a contributor
deciding a new case, and by a reviewer reading a diff.

**7 stages. 3 of the records cited below are still open. 2 decisions are
governed by nothing.**

## How to read it

Each stage names the decision it makes, the records that rule on it, and the clauses it is
decided under. A record's status is read from the ledger when this file is generated and is
never copied into the fixture, so a stage cannot be described as settled by a record that is
not.

The rows worth reading first are the **ungoverned** ones. Each is classified as either a
*genuine gap* — a decision with real consequences that no clause and no record settles — or as
*no record warranted*, meaning the spec settles it so plainly that writing a record would be
ceremony. Only the first kind is a finding.

## What this document is not

**It is not the ruleset.** What ferro's output is allowed to be, and what happens where the
spec determines no answer, is stated once, in
[README.md, *Normalization rules*](../README.md#normalization-rules), and is deliberately not
restated here.

**It is not the records.** It says which record governs a stage; it does not reproduce the
record. The records live in
[`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json).

**It is not exhaustive about clauses.** A stage's clause list carries the clauses that decide
it, not every clause that touches it.

## The ungoverned decisions

### Axis projection — genuine-gap

*How two `c.` positions in different numbering zones (`-n` upstream of the ATG, plain `n` in the CDS, `*n` downstream of the stop) are ordered against each other.*

Ferro needs a total order on `c.` endpoints to sort an allele's members, detect overlap between them, and compute the separation that the merge/split stage keys on. `background/numbering.md` defines each zone and states no rule for comparing positions across them; it never speaks about alleles at all. So the ordering is a house rule with no clause behind it — and, unlike every other house rule at this grain, it carries no ruling record either. It is currently written down only as repository prose. This is not theoretical: the sequence-changing and denotes-no-sequence classes the spec corpus found localise to the `c.72`/`c.*1` transition and to nothing else.

### Shifting — genuine-gap

*What the 5' direction does — at an exon/exon junction, in a repeat tract, and at a reference boundary.*

The recommendations state a 3' rule and only a 3' rule; `5' rule` and `5'-rule` occur nowhere in the checkout. Ferro nevertheless offers a 5' direction, and the direction is not an orthogonal knob — it selects the frame every other rule is evaluated in, which is why the README claims its best-effort rules per direction rather than across the two. So every clause and every ruling at this stage is authority for the 3' arm; the 5' arm's answers are ferro's own, mirrored by construction rather than by a decision anyone recorded. That is defensible and may well be the only possible answer, but it is nowhere adjudicated, and the mirroring is not obviously right at a junction, where the 3' answer is itself an exception.

## The stages

### Parse / validate (`parse-validate`)

**Decides.** Which authored forms are refused, which are repaired, and at which stage the refusal happens.

**Ruled by.**

- `absolute-prohibition-enforcement-stage` — decided
- `alignment-only-symbol-in-a-description` — decided
- `rna-axis-alignment-only-symbol-reach` — decided
- `bare-transcript-intronic-position` — decided
- `conflicting-member-geometry-refusal-scope` — decided
- `ring-telomere-anchoring` — undecided
- `rna-repeat-range-plus-unit-redundancy` — undecided

**Governed by nothing.** Nothing recorded at this stage.

The stage question itself is settled — `absolute-prohibition-enforcement-stage` rules that enforcement is mode-dependent and uniform: strict fails at parse, lenient fails only when it cannot normalize. What that record leaves open is coverage, clause by clause, and it names the remaining clauses itself. Two further questions at this stage carry open records rather than no record, which is the state this audit is looking for and does not need to add to.

### Axis projection (`axis-projection`)

**Decides.** How a coordinate maps across axes, and what happens at boundaries.

**Ruled by.**

- `projection-codon-exception-is-decided-by-the-rendered-axis` — decided
- `c-and-n-positions-are-flat-transcript-offsets` — decided
- `junction-exit-wrapper-scope-in-a-mixed-allele` — undecided

**Decided under.**

- `docs/background/numbering.md:52`
  > nucleotide numbering is `n.1`, `n.2`, `n.3`, ..., etc., from the first to the last nucleotide of the reference sequence

**Governed by nothing.**

- **genuine-gap** — How two `c.` positions in different numbering zones (`-n` upstream of the ATG, plain `n` in the CDS, `*n` downstream of the stop) are ordered against each other.

### Sequence derivation (`sequence-derivation`)

**Decides.** How the resulting sequence is computed from the authored edits, and what may bound that computation.

**Ruled by.**

- `canonical-form-choice-when-both-legal` — decided
- `derivation-may-not-be-bounded-by-the-inputs-spelling` — decided
- `confluence-gate-is-apply-equality-on-every-determined-axis` — decided

**Decided under.**

- `docs/background/basics.md:38`
  > The recommendations for the description of sequence variants are designed to be **stable**, **meaningful**, **memorable**, and **unequivocal**.

**Governed by nothing.** Nothing recorded at this stage.

The most heavily ruled stage in the pipeline, and the one the rest depend on: `canonical-form-choice-when-both-legal` makes the resulting sequence the thing every other rule is evaluated over, and `derivation-may-not-be-bounded-by-the-inputs-spelling` removes the one comparand that contradicted it. Note what `background/basics.md:38` is cited for here — it lists the design values and minimality is absent from them, so ferro's column-minimal objective is policy and never compliance.

### Partition + typing (`partition-and-typing`)

**Decides.** Which of several legal descriptions is emitted: how the edit set is cut into members, and what each piece is labelled.

**Ruled by.**

- `canonical-form-choice-when-both-legal` — decided
- `separation-is-a-property-of-the-spelling-not-of-the-variant` — decided
- `duplication-must-ranks-the-label-not-the-partition` — decided
- `inversion-vs-two-delins-76-83` — decided
- `inversion-vs-a-mixed-member-competitor` — decided
- `contiguous-insertion-split-by-a-blocked-derivation` — decided

**Decided under.**

- `docs/recommendations/general.md:55`
  > when a description is possible according to several types, the preferred description is: (1) substitution, (2) deletion, (3) inversion, (4) duplication, (5) insertion

**Governed by nothing.** Nothing recorded at this stage.

Settled at the level of principle and in flight at the level of shipped behaviour: several of these rulings are implemented only under a candidate partitioner arm, and the default is being flipped separately. #1553 places this stage out of scope for re-litigation, and this row is an index into the records rather than an invitation to reopen them.

### Shifting (`shifting`)

**Decides.** Where the 3' rule places a change within a run of equivalent positions, and where that placement is clamped.

**Ruled by.**

- `exon-junction-dup-converge-from-the-far-side` — decided

**Decided under.**

- `docs/recommendations/general.md:40`
  > **3'rule**: for all descriptions, the most 3' position possible of the reference sequence is arbitrarily assigned to have been changed.
- `docs/recommendations/general.md:43`
  > **exception**: deletions/duplications around exon/exon junctions using **c.**, **r.** or **n.** reference sequences
- `docs/recommendations/DNA/complex.md:50`
  > the general HGVS rule of maintaining the longest unchanged sequence applies (the 3' rule)

**Governed by nothing.**

- **genuine-gap** — What the 5' direction does — at an exon/exon junction, in a repeat tract, and at a reference boundary.

### Merge / split (`merge-split`)

**Decides.** When consecutive or separated changes combine into one description rather than staying individual.

**Ruled by.**

- `delins-merge-vs-individual-gap-two-or-more` — decided
- `delins-codon-carve-out-gap-one` — decided
- `delins-adjacent-members-when-both-consume-reference` — decided
- `codon-carve-out-shape-restriction` — decided
- `separation-rule-force-modal-or-negation` — decided
- `separation-is-a-property-of-the-spelling-not-of-the-variant` — decided
- `self-cancelling-across-ring-junctions` — decided

**Decided under.**

- `docs/recommendations/general.md:33`
  > two variants separated by one or more nucleotides should be described individually and **not** as a "delins"
- `docs/recommendations/DNA/substitution.md:32`
  > changes involving two or more consecutive nucleotides are described as deletion-insertion (delins)

**Governed by nothing.** Nothing recorded at this stage.

Seven records, more than any other stage, and the residue is named inside them rather than missing from them: `delins-adjacent-members-when-both-consume-reference` is explicitly scoped to members that consume reference bases and records the `dup` and repeat-expansion cases as still open. An audit row that reported this stage as fully settled would be reading the ruling and not its scope.

### Re-spelling (`re-spelling`)

**Decides.** Terminal coalescing and type re-selection, applied to pieces that are already derived.

**Ruled by.**

- `delins-payload-coincidence-carve-out-is-coding-dna-scoped` — decided
- `delins-recommendation-reach-when-the-input-arrives-split` — decided

**Governed by nothing.** Nothing recorded at this stage.

Both records here are axis- or mechanism-scoped rather than general, and both were decided in the direction of narrowing what the carve-out reaches. Neither is live on the shipped default, so this stage's rulings and this stage's behaviour are currently different questions.

## The axis dimension

The recommendations are not uniform across axes, and the non-uniformity is deliberate rather
than accidental. `general.md:162` settles that outright: asked whether the `g.` and `c.`
descriptions of one variant may name different nucleotides, the spec answers **yes**, for a
gene on the minus strand inside a repeated sequence.

That matters because a normalizer's axes disagreeing has, in this repository, always turned out
to be a bug — a dropped offset, a skipped pass, a missing gate. Every such case so far was an
*accidental* divergence, and none was a deliberate axis-scoped rule. Without a table recording
which is which, the next deliberate divergence is indistinguishable from the next defect.

Two things shape the table. **Rules are keyed on the property, not on the axis label**: what
matters for the codon carve-out is having a reading frame, which correctly groups `n.` with
`g.` rather than with `c.`, a pairing an axis-name-keyed table would get wrong. And **some
rules cannot be modelled at all** — a clause keyed on provenance or on population frequency is
not a function of the sequences a normalizer holds — so those are recorded as `unmodellable`
rather than left looking unimplemented.

| clause | scope | normalization-relevant? |
|---|---|---|
| `docs/recommendations/general.md:34` | Coding only by construction — the condition names amino acids, which an axis with no reading frame cannot state. | modelled |
| `docs/recommendations/general.md:43` | `c.`, `r.` and `n.` only — the spec scopes the exception by naming prefixes. | modelled |
| `docs/recommendations/DNA/repeated.md:21` | `c.` only. | modelled |
| `docs/recommendations/general.md:161-166` | The `g.`/`c.` pair, on the minus strand, inside a repeated sequence. | modelled |
| `docs/recommendations/general.md:13` | RNA and protein axes. | modelled |
| `docs/recommendations/general.md:135` | `c.` only. | not-a-normalization-rule |
| `docs/recommendations/RNA/adjoined_transcript.md:18` | `r.` only. | not-a-normalization-rule |
| `docs/background/refseq.md:47` | Intronic positions on transcript axes. | not-a-normalization-rule |
| `docs/recommendations/DNA/delins.md:83` | Axis-neutral. Listed here because it is the clearest case of a rule this table must record as unmodellable rather than leave looking unimplemented. | unmodellable |
| `docs/recommendations/RNA/delins.md:41` | `r.`. | unmodellable |

### The clauses, with their text

#### `docs/recommendations/general.md:34` — modelled

> **exception**: two variants separated by one nucleotide, together affecting one amino acid, should be described as a "delins"

**Scope.** Coding only by construction — the condition names amino acids, which an axis with no reading frame cannot state.

Keyed on the property, not on the axis label: `axis_min_separation(reading_frame: bool)` selects the floor, which correctly groups `n.` — a transcript axis with no frame — with `g.` rather than with `c.`. An axis-name-keyed table would get that pair wrong.

#### `docs/recommendations/general.md:43` — modelled

> **exception**: deletions/duplications around exon/exon junctions using **c.**, **r.** or **n.** reference sequences

**Scope.** `c.`, `r.` and `n.` only — the spec scopes the exception by naming prefixes.

Also precedent, cited as such by `projection-codon-exception-is-decided-by-the-rendered-axis`: the spec does scope a rule by naming prefixes, so an axis-scoped reading of another clause is not an invention.

#### `docs/recommendations/DNA/repeated.md:21` — modelled

> using a coding DNA reference sequence ("c." description), a repeated sequence variant description can be used only for repeat units with a length which is a multiple of 3

**Scope.** `c.` only.

#### `docs/recommendations/general.md:161-166` — modelled

> Yes, when a gene is on the minus strand of a chromosome (opposite transcriptional orientation) and the change is located in a repeated sequence (mono-, di-, tri-, etc. nucleotide stretches), the 3'rule has this as a consequence.

**Scope.** The `g.`/`c.` pair, on the minus strand, inside a repeated sequence.

The load-bearing row of this table, and the reason the table exists. The Q&A answers **yes**: `g.` and `c.` descriptions of one variant may legitimately name different nucleotides and different positions. So cross-axis divergence is endorsed by the spec, not a defect — while this repository has historically filed axis disagreement as a bug, correctly in every case so far, because every one of those was an accidental divergence rather than a deliberate axis-scoped rule. Without a table recording which is which, the next deliberate divergence is indistinguishable from the next bug.

#### `docs/recommendations/general.md:13` — modelled

> descriptions on RNA/protein level should describe the changes observed on that level (RNA/protein) and not try to incorporate any knowledge regarding the change on

**Scope.** RNA and protein axes.

The method half of `canonical-form-choice-when-both-legal`, whose governing clause is `general.md:156` — the same statement made for protein and extended here to RNA.

#### `docs/recommendations/general.md:135` — not-a-normalization-rule

> The minus sign should only be used as a minus in the description of variants based on a coding DNA reference sequence.

**Scope.** `c.` only.

A parser constraint. It decides which strings are well formed, not which of several well-formed descriptions is emitted.

#### `docs/recommendations/RNA/adjoined_transcript.md:18` — not-a-normalization-rule

> This syntax is for RNA sequence only (no use of coding (`c.`) / non-coding DNA (`n.`) reference sequences).

**Scope.** `r.` only.

Expressibility. It says a syntax is unavailable on an axis, which is upstream of any choice between forms.

#### `docs/background/refseq.md:47` — not-a-normalization-rule

> intronic sequences are considered to be **within the boundaries** of a transcript reference sequence and may only be used to describe a variant when a genomic r

**Scope.** Intronic positions on transcript axes.

Reference selection rather than form selection — though the sibling clause `checklist.md:20` is what `bare-transcript-intronic-position` rules under, and that ruling does move output.

#### `docs/recommendations/DNA/delins.md:83` — unmodellable

> First, the two variants may have been reported (or might occur) individually.

**Scope.** Axis-neutral. Listed here because it is the clearest case of a rule this table must record as unmodellable rather than leave looking unimplemented.

The discriminator is **provenance** — whether the two variants were reported individually — which is not a function of (reference, observed) and is unrecoverable from any description. `canonical-form-choice-when-both-legal` records this among its counter-evidence and adopts re-derivation over it deliberately; `separation-is-a-property-of-the-spelling-not-of-the-variant` records the departure formally, in its `deviates_from`.

#### `docs/recommendations/RNA/delins.md:41` — unmodellable

> This format is preferred when either of the two variants is known as a frequently occurring variant ("polymorphism").

**Scope.** `r.`.

Population frequency, like the row above, is not a function of the sequences a normalizer holds. Recording it as unmodellable is the point: an audit that left it blank would read as a gap someone could close.
