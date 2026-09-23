<!--
GENERATED FILE — do not edit by hand.

Rendered from the `rulings` section of
tests/fixtures/grammar/hgvs_spec_normalization_overrides.json
by tests/it/normalization_contract_doc.rs. Edit the ledger, then regenerate:

    BLESS_CONTRACT_DOC=1 cargo nextest run --features dev --test it \
      -E 'test(normalization_contract_doc)'

An edit made here instead is reverted by the next regeneration, and fails CI
before that.
-->

# Ferro's normalization contract

**38 adjudication records — 35 decided, 3 open.**

## What this document is

The HGVS recommendations are, in places, silent, ambiguous, or self-contradictory. A
normalizer still has to emit one string. Where ferro has had to decide such a question, the
decision is recorded as a **ruling record**, and this document is an index of every one of
those records: its status, the clause that governs it where one applies (a decided record names
its governing clause or, as a house choice, cites none; an open record names no authority), and
the record's own one-sentence statement of the ruling.

It is generated, not written. Each ruling summary is the record's own `summary` field, copied
verbatim from
[`tests/fixtures/grammar/hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json);
nothing here is paraphrased by the renderer. The full reasoning behind each ruling, the spec
clauses it quotes, the text each clause was quoted against, and any scope on the ruling are
**not** reproduced here — they live in the ledger, which is what the build enforces.

## What this document is not

**It is not the ruleset.** What ferro's output is allowed to be — which properties are
absolute and which are best effort, what happens where the spec determines no answer, and what
must be disclosed when a choice changes — is stated once, in
[the normalization rules](src/reference/normalization-rules.md). That statement is
deliberately not reproduced here, and single-sourcing it is itself part of one of the rulings
below. Where this document and the ruleset page appear to disagree, the ruleset page governs and this
document has a bug.

**It is not a substitute for the records.** It is a reading of them. The records are what the
build enforces.

**It is not a spec.** The HGVS recommendations are upstream — rendered at
[hgvs-nomenclature.org](https://hgvs-nomenclature.org/) and sourced from the
`HGVSnomenclature/hgvs-nomenclature` repository that `assets/hgvs-nomenclature` vendors. Each
governing clause below is spelled as `path:line` and links to that exact line in the **pinned**
commit of that repository, so the citation resolves to the spec version ferro actually pins
rather than to whatever the site currently shows.

## How to read the index

- **The id states the QUESTION, not the ruling, and the two can be opposites.** Read the
  ruling column, not the id. One record below is titled as the position it *rejects*.
- **`undecided` is a first-class state**, not an oversight. An open record states a conflict
  and declines to settle it; whatever ferro does with that case today is the status quo, and
  citing the behaviour as a decision is the error the record exists to prevent.
- **A one-line ruling states the decision, not how to carry it out.** Many records add an
  explicit scope — one axis, one direction, one shape — and the mechanism or tie-break that
  turns the decision into an output; this index shows neither. Read it to scan every ruling,
  and open the ledger record before *acting* on one. A summary that reads like a complete rule
  may still be narrower, or more conditional, than it looks.

## Which build this describes

The rulings are decisions about what ferro's output *should* be, and they do not depend on
which block partitioner is selected at run time. Several records are nevertheless explicit
about whether their ruling is already live in shipped output or is implemented only under a
candidate arm, so the answer depends on the default — and that default is currently in motion.

**As generated, the shipped default — what `FERRO_PARTITION` unset selects — is
`canonical-coalesced`.** That sentence is not written by hand: the generator reads the arm out of
`src/normalize/merge.rs`, so a change to the default fails this document's own test until it
is regenerated, and cannot leave a stale claim behind. See
[Comparing normalization rules](src/guide/comparing-rules.md#comparing-normalization-rules-ferro_partition)
for the knob and its traps.

## The records

Every ruling record on one screen: its status, the clause that governs it where one applies (a decided record names its governing clause or, as a house choice, cites none; an open record names no authority), and a one-sentence statement of the ruling. Read the ruling, not the id — an id states the record's *question*, and at least one of the questions below is answered in the negative. An `open` record names a conflict ferro has **not** settled; whatever ferro does with that case today is the status quo, not a ruling, and must not be cited as one. The full reasoning, the clauses each record quotes, and any scope on the ruling live in the ledger itself, [`hgvs_spec_normalization_overrides.json`](../tests/fixtures/grammar/hgvs_spec_normalization_overrides.json); this document is an index into it, not a copy of it.

| Record | Status | Governing clause | Ruling |
|---|---|---|---|
| `junction-exit-wrapper-scope-in-a-mixed-allele` | undecided | — | *Undecided.* An allele holds one member whose intronic offset ferro MANUFACTURED and one the author spelled intronic, both on one bare transcript. `checklist.md:20` wants a genomic reference for the manufactured member; `DNA/alleles.md:16` wants one reference factored out of the whole allele. Lift the wrapper to the whole description, or expand to per-member accessions? |
| `ring-telomere-anchoring` | undecided | — | *Undecided.* Must a ring chromosome's first `::` segment start at `pter` and its last end at `qter`, so that a ring whose segments name only interior coordinates is refused? |
| `rna-repeat-range-plus-unit-redundancy` | undecided | — | *Undecided.* On an RNA reference, may a repeat be written with BOTH a position range and the repeat unit (`r.-6_-3g[6]`)? `RNA/repeated.md:22` marks that shape invalid as redundant; `:27` publishes exactly that shape as valid. |
| `absolute-prohibition-enforcement-stage` | decided | [`docs/recommendations/checklist.md:5`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/checklist.md?plain=1#L5) | Spellings the spec prohibits are rejected — at parse in strict mode; lenient mode instead repairs the input where it can and fails only if it cannot normalize. |
| `adjudication-precedence-order` | decided | [`docs/recommendations/DNA/complex.md:63`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/complex.md?plain=1#L63) | The canonical normalization ruleset is stated once, in ferro's reference documentation, and this record deliberately does not restate it — it holds only the escalation register for the rare case where the spec prefers one conformant form but supplies nothing to choose between the candidates that satisfy it. |
| `alignment-only-symbol-in-a-description` | decided | [`docs/background/standards.md:39`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/background/standards.md?plain=1#L39) | The alignment-only symbols 'X' and '-' are refused in a description, since they are not among the IUPAC-IUBMB nucleotide symbols a description may use. |
| `bare-transcript-intronic-position` | decided | [`docs/recommendations/checklist.md:20`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/checklist.md?plain=1#L20) | A bare transcript intronic position such as 'NM_...:c.20+2del' is refused in strict input-hygiene mode and accepted in lenient mode, since the spec allows it only when a genomic reference is given. |
| `c-and-n-positions-are-flat-transcript-offsets` | decided | [`docs/background/numbering.md:52`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/background/numbering.md?plain=1#L52) | A 'c.'/'n.' position names the base at that offset in the flat transcript sequence and is never resolved by walking the exon table, so it lands correctly even on a transcript whose cdot alignment carries a gap. |
| `c-description-against-an-unresolvable-cds-is-refused` | decided | [`docs/background/numbering.md:21`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/background/numbering.md?plain=1#L21) | A 'c.' description against a transcript whose CDS start the reference cannot resolve is refused rather than returned unchanged, because ferro has no origin to number the coordinate from. |
| `canonical-form-choice-when-both-legal` | decided | [`docs/recommendations/general.md:156`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/general.md?plain=1#L156) | When two descriptions of one variant are both legal and no clause chooses between them, ferro applies the edit set to the reference and re-derives the description from the resulting sequence — subject to every explicit spec tie-break, and breaking any remaining tie toward the already-shipped form — rather than preserving the input's spelling. |
| `cds-unknown-position-is-refused-at-conversion` | decided | house choice (cites none) | Converting a 'c.?' (unknown) position toward a transcript coordinate is refused rather than coerced into a concrete position — a project choice, since fabricating a coordinate for an explicitly unknown one yields a storable-but-wrong result. |
| `coding-axis-merges-are-a-disclosed-general-34-deviation` | decided | [`docs/recommendations/DNA/delins.md:81`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L81) | A disclosure that currently fires on nothing: it records a rule-2 deviation from the spec's preference to describe members individually — a coding-axis merge of allele members two or more unchanged nucleotides apart with no gap-bearing member — but on the current partition rule that set is empty, so the record is a candidate for narrowing or retirement rather than a live rule. |
| `codon-carve-out-shape-restriction` | decided | [`docs/recommendations/DNA/delins.md:18`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L18) | Two changes one nucleotide apart that together alter a single amino acid are ruled one delins whatever the edit types — because "together affecting one amino acid" is a fact about the resulting sequence, not about how the input was spelled — though the coding-codon pass currently implements only the substitution/unchanged-base/substitution triplet. |
| `conflicting-member-geometry-refusal-scope` | decided | [`docs/recommendations/DNA/alleles.md:5`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/alleles.md?plain=1#L5) | Two members of one allele that claim intersecting reference territory — nested, overlapping, or two insertions at one interbase — are refused, whatever edit types they render as. |
| `confluence-gate-is-apply-equality-on-every-determined-axis` | decided | [`docs/recommendations/general.md:42`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/general.md?plain=1#L42) | Ferro's release gate asserts that inputs denoting the same sequence on every determined axis, the protein axis excluded, normalize to one output; it is asserted over decided equivalence classes only, and equivalence is judged by applying the descriptions to the reference rather than by comparing normalized strings. |
| `contiguous-insertion-split-by-a-blocked-derivation` | decided | [`docs/recommendations/general.md:33`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/general.md?plain=1#L33) | A single contiguous insertion spelled as an insertion plus a duplication one base apart is ruled to be the one insertion, because the separation rule keys on two variants and this locus carries only one — a decided target not yet implemented, so ferro still emits the multi-member form when a distant third member blocks the derivation. |
| `cross-zone-c-positions-order-by-transcript-coordinate` | decided | house choice (cites none) | Two `c.` positions in different numbering zones (`c.-n`/`c.n`/`c.*n`) are ordered by their flat transcript-sequence offset, not by the `c.` spelling — a house choice under `normalization-rules.md` rule 5's silent limb, because the recommendations define the zones but state no rule for comparing across them. |
| `delins-adjacent-members-when-both-consume-reference` | decided | [`docs/recommendations/DNA/substitution.md:32`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/substitution.md?plain=1#L32) | Two adjacent changes that both consume reference bases are written as a single delins; the spec marks the split spelling "not correct" at separation zero. |
| `delins-codon-carve-out-gap-one` | decided | [`docs/recommendations/DNA/delins.md:18`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L18) | Two changes one nucleotide apart that together affect a single amino acid are written as one delins on the coding sequence, the explicit exception the spec makes to describing them individually. |
| `delins-merge-vs-individual-gap-two-or-more` | decided | [`docs/recommendations/DNA/delins.md:47`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L47) | When two changes two or more nucleotides apart arise only because part of the inserted sequence coincides with the reference and the block is a net deletion, ferro writes them as one spanning delins rather than individually. |
| `delins-payload-coincidence-carve-out-is-coding-dna-scoped` | decided | [`docs/recommendations/DNA/delins.md:47`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L47) | Where a split exists only because payload bases coincide with the reference, ferro writes it as one spanning delins on every DNA axis (c./g./m./n., but not r.); on the frameless axes this is a disclosed rule-2 deviation and the project's choice among conformant forms. |
| `delins-recommendation-reach-when-the-input-arrives-split` | decided | [`docs/recommendations/DNA/delins.md:46`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L46) | Ferro merges a re-derived split into one delins only when some member supplies inserted bases while consuming a different number of reference bases; a split of pure deletions inserts nothing and stays individual. |
| `derivation-may-not-be-bounded-by-the-inputs-spelling` | decided | [`docs/recommendations/DNA/delins.md:47`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L47) | Ferro never refuses a re-derived description for naming more change than the input's own spelling did; the derived form is chosen from the sequence, not weighed against how the input was written. |
| `duplication-must-ranks-the-label-not-the-partition` | decided | [`docs/recommendations/DNA/duplication.md:17`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/duplication.md?plain=1#L17) | The rule that a duplication must be labelled 'dup' ranks the label of each piece ferro derives, not the partition; the one exception is a net-longer tandem copy of a multi-base motif, where the derivation is cut to expose the dup rather than merged into a delins. |
| `exon-junction-dup-converge-from-the-far-side` | decided | [`docs/recommendations/DNA/duplication.md:26`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/duplication.md?plain=1#L26) | A duplication is placed at the most 3' position that does not cross an exon/exon junction, reached from either side, so a copy spelled past the junction is pulled back to it. |
| `inversion-vs-a-mixed-member-competitor` | decided | [`docs/recommendations/DNA/inversion.md:5`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/inversion.md?plain=1#L5) | When a span replaced by its reverse complement competes with a description mixing lone substitutions and multi-column members, ferro writes it as one inv; both forms are conformant, so this is the project's choice among them. |
| `inversion-vs-two-delins-76-83` | decided | [`docs/recommendations/DNA/inversion.md:5`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/inversion.md?plain=1#L5) | A span replaced by its exact reverse complement is written as a single inv even where its interior columns coincide with the reference, not as the two delins those columns would separate. |
| `inverted-duplication-is-derived-as-ins-range-inv` | decided | [`docs/recommendations/DNA/inversion.md:69`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/inversion.md?plain=1#L69) | An inverted duplication is written as 'ins<range>inv', naming the span the inverted copy came from, rather than expanded to reverse-complemented literal bases; whether a payload counts as an inverted copy at all is gated by a house coincidence-probability floor, not any spec-stated minimum — and it is wired today only in normalize_genome, so c./n./r./m. still emit literals. |
| `past-cds-end-coordinate-is-non-conformant` | decided | [`docs/background/numbering.md:21`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/background/numbering.md?plain=1#L21) | A 'c.N' coordinate past the CDS end is non-conformant: strict mode refuses it, while lenient and silent modes repair it to the equivalent 'c.*' position — applied the same way to a lone position and to a member of a cis allele. |
| `projection-codon-exception-is-decided-by-the-rendered-axis` | decided | [`docs/recommendations/DNA/delins.md:42`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L42) | The codon merge fires only on an axis that declares a reading frame, so when a coding description merges under it ferro leaves the members individual on the derived genomic axis rather than re-merging it to match. |
| `rna-axis-alignment-only-symbol-reach` | decided | [`docs/background/standards.md:47-61`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/background/standards.md?plain=1#L47-L61) | A non-leading lower-case 'x', which is not an assigned RNA nucleotide symbol, is refused in an 'r.' description, mirroring the refusal of the alignment-only 'X' on the DNA axes. |
| `self-cancelling-across-ring-junctions` | decided | [`docs/recommendations/DNA/complex.md:130`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/complex.md?plain=1#L130) | The prohibition on replacing part of a sequence with part of itself does not reach the '::'-joined segments of a ring chromosome, so ferro keeps a ring's members rather than collapsing them to a linear delins. |
| `separation-is-a-property-of-the-spelling-not-of-the-variant` | decided | [`docs/recommendations/general.md:33`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/general.md?plain=1#L33) | Ferro reads the separation between changes off the partition it re-derives from the resulting sequence, not off the input's spelling, so two spellings of one variant converge on one output. |
| `separation-rule-force-modal-or-negation` | decided | [`docs/recommendations/DNA/delins.md:81`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L81) | Two changes a nucleotide or more apart are described individually — this is the spec's preference (ruleset rule 2), not an outright ban; the only spelling the recommendations forbid is the split at separation zero. |
| `spdi-n-unit-repeat-refusal` | decided | house choice (cites none) | On the HGVS-to-SPDI path ferro refuses an 'N'-unit or 'N'-containing repeat rather than expand it to literal 'N' bases — a project choice, since 'N' states a length, not identified bases, and the recommendations do not reach SPDI conversion. |
| `unchanged-is-read-over-every-minimal-alignment` | decided | house choice (cites none) | Ferro treats a reference base as unchanged only when every minimum-edit alignment of the block matches it — a project choice the recommendations neither require nor forbid — except that an equal-length run with no base surviving at its own coordinate is rendered as one spanning delins. |
| `unequal-length-block-a-placed-gap-is-not-a-separation` | decided | [`docs/recommendations/DNA/delins.md:47`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/delins.md?plain=1#L47) | A lone unequal-length net-deletion delins whose payload merely coincides with the reference, and whose every member other than the placed gap would itself render as a delins, is kept whole on every DNA axis (c./g./m./n., but not r.) rather than split at that coincidence into a separate residual member. |
| `whole-span-reverse-complement-types-as-inv` | decided | [`docs/recommendations/DNA/inversion.md:5`](https://github.com/HGVSnomenclature/hgvs-nomenclature/blob/565b9734a4a44fa9bfaa6264e3a3ff092e5e82af/docs/recommendations/DNA/inversion.md?plain=1#L5) | A span whose whole content is replaced by its exact reverse complement is written as one inv, however much of its interior coincides with the reference and whatever the competing partition is made of; this is a project choice among conformant forms, not a conformance requirement. |
