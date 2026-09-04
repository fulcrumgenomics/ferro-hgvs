# Reading the HGVS recommendations

Findings about the spec text in `assets/hgvs-nomenclature` that change how a clause argument
should be framed. The submodule is the HGVS nomenclature repository, pinned at one commit. Every
`file:line` below cites it. Open the cited lines; do not trust this summary. Where a reading is a
house choice, the section names the ruling record in `docs/NORMALIZATION_CONTRACT.md`. The record
holds the reasoning.

## Almost nothing in the recommendations is normative

Do not rank clauses by RFC 2119 keyword strength. `style.md:9` binds the document to RFC 2119,
but the recommendations almost never use the keywords. At the pinned commit, uppercase RFC 2119
keywords appear exactly twice outside `style.md`, both in one file:

```bash
grep -rnoE '\b(MUST|SHOULD|RECOMMENDED|MAY|SHALL|REQUIRED|REQUIRES|OPTIONAL)( NOT)?\b' \
  assets/hgvs-nomenclature/docs/recommendations/ | grep -v '/style\.md'
# -> assets/hgvs-nomenclature/docs/recommendations/RNA/adjoined_transcript.md:20:REQUIRES
# -> assets/hgvs-nomenclature/docs/recommendations/RNA/adjoined_transcript.md:21:SHOULD
# (and nothing else)
```

`REQUIRES` is not in `style.md:9`'s list. The census counts it as the authors' intent: a
capitalised form of REQUIRED in a document bound to RFC 2119. Under the strict list the count is
one, at `:21`.

Every merge and split clause this repository argues over (`general.md:33`, `:34`, `:55`, `:57`,
`DNA/delins.md:16`, `:17`, `:18`, `:47`, `DNA/inversion.md:20`) is lowercase prose. Read strictly,
none is normative. They are not ignorable: each raises a house choice the spec leaves open, which
ferro must answer and record. Argue them from the spec's worked examples and from downstream cost.

Some clauses read as prohibitions because of their wording: `general.md:57` says such
descriptions "are not allowed"; `DNA/duplication.md:18` says a duplication "**must**" be
described as one. They are strong because of their words, not because of RFC 2119. Say which you
mean.

## There is no minimality principle, and stability is the spec's first stated value

Never cite the spec for minimality. `background/basics.md:38` states the design values: "stable,
meaningful, memorable, and unequivocal". Minimality is not among them; stability is first. An
argument from length, or from bases touched, is a house argument.

The spec's own worked example prefers the larger form. `DNA/delins.md:44` gives
`LRG_199t1:c.850_901delinsTTCCTCGATGCCTG`, a delins that replaces 52 bases. `:46` gives "an
alternative description" as a split, `c.[850_869del;874_881del;887_897del;901_902insG]`, which
deletes 39 bases and inserts one. `:47` recommends the delins: "it is simpler and prevents
software tools making incorrect predictions for the consequences on protein level". State that
direction when you cite the passage, and cite the ruling that scopes it,
`rulings[delins-merge-vs-individual-gap-two-or-more]`. That is the spec preferring the larger form *for this example's stated reasons*, not a general "prefer the larger delins" — the scope is that ruling, nothing wider.

## A split is rarely unique, which is a stability argument by itself

The spec has two explicit tie-breaks: type prioritisation at `general.md:55` and the 3' rule at
`general.md:40`. Beyond those it gives no rule for choosing among equally compliant splits of one
change, and it never claims a split is unique. Before you adopt a split, count the alternatives.
In the rows ferro enumerated, most admit more than one compliant set of members. Picking one
trades a stable canonical form for an arbitrary member of a family, and one before-and-after
comparison cannot show that. Check the family, not the row. An argument for a split is an
argument against `basics.md:38`'s stability; say so.

There are two kinds of stability. `basics.md:38` means a variant keeps one description across
time and tools. Ferro also requires determinism, rule 4 of
`docs/src/reference/normalization-rules.md`: same input, same output, on every run. That is not
confluence, rule 3, which is about inputs denoting one variant producing one output. The spec
says nothing about determinism; do not cite it for that. Where ferro must pick, the pick depends
on the input alone, never on iteration order or scheduling, and it is a house choice to record.
For splits the record is `rulings[canonical-form-choice-when-both-legal]`: re-derive from the
resulting sequence, then make a deterministic choice among equally-minimal partitions. The record
also holds the enumeration behind "most": 27 of 40 rows admit more than one split, the worst
admits 125, and the spec's own split at `DNA/delins.md:46` is one of five.

## Read the Q&A, not only the Notes and Examples

Read each recommendations page to its end. The `!!! note` blocks at the foot are community
questions and suggestions with the committee's answers. They can settle what the body leaves
open, and they can record a changed decision.

`DNA/delins.md:79-84` answers whether two nearby substitutions are one delins or two variants:
two, unless together they affect one amino acid. The first reason is that "the two variants may
have been reported (or might occur) individually". That reason is provenance: how the changes
were observed. No rule over the sequence alone can recover it. So any sequence-local merge or
split rule in ferro is a house approximation of the spec's criterion; argue it as one. The
governing records are `rulings[delins-merge-vs-individual-gap-two-or-more]` and
`rulings[canonical-form-choice-when-both-legal]`. Provenance is how the change was observed, not "preserve the input's spelling": ferro re-derives from the resulting sequence, so the input's spelling gets no weight.

`DNA/delins.md:86-89` answers a BRCA1 case: a substitution plus an adjacent insertion is
`NM_007294.3:c.2077delinsATA`. `:89` records that a sentence permitting the two-member spelling
was removed. Do not cite the removed text.

## Cross-check forward-looking statements against `consultation/`

A sentence that describes a future rule is not a rule. The sign is future tense about the
recommendations themselves: "is preparing", "is being prepared", "will suggest", "the new
recommendation will be". Before you cite or build on such a sentence, find its proposal in
`consultation/` and read the status. If the proposal is rejected or undecided, the current text
stands. If it is rejected and you write a test about it, the test must fail when ferro produces
what the proposal asked for; it must never expect that output. No forward-looking note at the
pinned commit points at an accepted proposal. If one appears, record an adjudication in the
ledger; do not infer a reading from this page. A reference to an adopted proposal, such as
"based on Community Consultation proposal SVD-WG005" at `DNA/other.md:41`, is not
forward-looking. The text around it is current.

The notes found at the pinned commit:

- `general.md:35-38`: "the SVD-WG is preparing a proposal to modify this recommendation ... The
  new recommendation will be: two variants separated by less than two nucleotides should be
  described as a delins." Read alone, it says the rule at `:33-34` will change. The proposal is
  SVD-WG010: `consultation/SVD-WG010.md:12` states the same rule, `:27` the same reason, and
  `:30-33` names `general.md:33-34` as the rule to replace. It was rejected (`SVD-WG010.md:5`).
  The note was not updated. The rule at `:33-34` stands, codon exception included.
- `DNA/duplication.md:86` marks the worked example above it as "part of proposal SVD-WG003
  (undecided)". `SVD-WG003.md:5` reads "new proposal to be made"; `:10` says the example's format
  "does follow current recommendations". The example stands.
- `DNA/repeated.md:9` and `RNA/repeated.md:9` say a proposal "is being prepared" to require the
  full range in repeat descriptions. No such proposal exists in `consultation/`; the topic is in
  `consultation/open-issues.md:82-110`. The syntax after the note stands.

## Comparing `c.` positions across numbering zones is OUR policy, not compliance

Never cite the spec for how two `c.` positions in different numbering zones compare. The `c.`
axis has three zones: `-n` upstream of the ATG, plain `n` in the CDS, `*n` downstream of the
stop. `background/numbering.md` defines each zone and the flat `n.` numbering under them. It
states no rule for comparing positions across zones, and it never mentions alleles or members:
zero occurrences of "allele", zero of "member". `refseq.md`'s only uses of "allele" (`:264-265`)
are the population-genetics sense.

Ferro needs such an order: a cis allele's members can sit in different zones, and normalization
must sort them, detect overlap and decide separation. The order is a house choice,
`rulings[cross-zone-c-positions-order-by-transcript-coordinate]`: positions are ordered by their
flat transcript coordinate, and the `c.` spelling is a presentation of that coordinate. Argue
from the coordinate. The failure point is the CDS/3'UTR boundary, where the last coding base and
`c.*1` are adjacent on the transcript and far apart in zone-local arithmetic. The defect class is
pinned in `tests/it/spec_corpus_regressions.rs`.
