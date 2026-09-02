# Reading the HGVS recommendations

Findings about the spec text in `assets/hgvs-nomenclature` that change how a clause
argument should be framed. Each is checkable against the pinned checkout; re-check
rather than trust. Moved out of `CLAUDE.md` on 2026-09-02.

## Almost nothing in the recommendations is normative

`style.md:9` binds the spec to RFC 2119, which invites the habit of arguing that clause A is a
SHOULD and clause B a MUST. That argument almost never works here, because **uppercase RFC 2119
keywords appear exactly twice outside `style.md` itself**, and both are in one file:

```bash
grep -rnoE '\b(MUST|SHOULD|RECOMMENDED|MAY|SHALL|REQUIRED|REQUIRES|OPTIONAL)( NOT)?\b' \
  assets/hgvs-nomenclature/docs/recommendations/ | grep -v '/style\.md'
# -> docs/recommendations/RNA/adjoined_transcript.md:20:REQUIRES
# -> docs/recommendations/RNA/adjoined_transcript.md:21:SHOULD    (and nothing else)
```

Every clause this project has litigated — `general.md:33`, `:34`, `:55`, `:57`,
`DNA/delins.md:17`, `:18`, `:47`, `DNA/inversion.md:20` — is lowercase prose. Read strictly,
none of them is normative. That does not make them ignorable: it makes most of these questions
house-style choices the spec leaves open, which ferro still has to answer. Argue them on the
spec's worked examples and on downstream cost, not on keyword strength. (`general.md:57` and
`DNA/duplication.md:18` are read as strong because of their wording — "are not allowed",
"**must**" — not because RFC 2119 makes them so; say which you mean.)

The `rulings[delins-merge-vs-individual-gap-two-or-more]` record used to argue from "the two
clauses are the same RFC 2119 strength". It does not any more, for this reason.

## There is no minimality principle, and stability is the spec's first stated value

`background/basics.md:38`: "The recommendations for the description of sequence variants are
designed to be **stable**, **meaningful**, **memorable**, and **unequivocal**." Minimality is
absent, and stability leads. Ferro's column-minimal objective is therefore **our policy**, not
compliance — never cite the spec for it. `DNA/delins.md:44-47` in fact recommends a
*non-minimal* description in its own worked example (66 columns where 40 suffice), and `:47`
recommends the **spanning delins**; the split is what `:46` calls the "alternative description".
That direction has been retold backwards at least once, so state it explicitly whenever you
cite the passage.

## A split is rarely unique, which is a stability argument by itself

Exact enumeration over 40 rows found **27 that admit more than one equally-compliant split**,
median 2 and **max 125**; the spec's own `:44-47` example admits five. So adopting a split
trades one stable canonical form for an arbitrary pick out of a family — and the arbitrariness
is invisible in any single-row before/after comparison, which is how it keeps getting missed.

## Read the Q&A, not only the Notes and Examples

`DNA/delins.md`'s trailing `!!! note` blocks are adjudicative, not colour. A claim that "no
sequence-local rule separates these two cases" was withdrawn once `:86-89` was read: the Q&A
carries the worked answer (`NM_007294.3:c.2077delinsATA`) and `:89` records that the passage
permitting the two-member spelling was *removed* by the committee. `:79-84` likewise carries the
spec's own discriminator — "the two variants may have been reported (or might occur)
individually" — which is **provenance**, recoverable only from the input's spelling.

## A forward-looking note is a suggestion, and may describe a proposal that FAILED

`general.md:35-38` reads as forward guidance — "the SVD-WG is preparing a proposal… The new
recommendation will be: two variants separated by less than two nucleotides should be described
as a `delins`" — and it is tempting to treat it as the direction of travel and discount the
current rule accordingly.

**That proposal is SVD-WG010, and it was rejected.** Word for word, including its rationale.
So the note is stale text describing a change that never happened, and three conclusions drawn
from reading it as guidance are all withdrawn: the codon exception is *not* going away, the
proposal to replace it having failed; it *strengthens* `codon-carve-out-shape-restriction`
rather than undercutting it; and the "spec-admitted instability" argument built on it does not
stand.

So a rejected proposal earns a **negative guard**, not an expectation. The spec corpus builds
210 rows whose only purpose is to catch a frameless separation floor of two — what implementing
SVD-WG010 looks like from the outside — and asserts `guard_violations == 0` over them, with the
denominator asserted non-zero so `0 of 0` cannot pass as a result. See
`spec_conformance_axis.rs`.

**So: cross-check every forward-looking statement in `recommendations/` against the
`consultation/` dispositions before citing it.** "is preparing", "will be", "the new
recommendation" are all flags. The disposition table is inventoried in the consultation slice:
9 accepted, 3 rejected, 8 open, 3 unclear — a rejected proposal is generated as a NEGATIVE
guard, never as an expectation.

## Comparing `c.` positions across numbering zones is OUR policy, not compliance

`c.` has three numbering zones — `-n` upstream of the ATG, plain `n` in the CDS, `*n` downstream
of the stop — and a cis allele's members can sit in different ones. Ferro therefore needs an
endpoint ordering that spans zones, to sort members, detect overlap and decide separation.

**The spec does not give one.** `background/numbering.md` defines each zone but states no rule
for comparing positions across them, and it never speaks about alleles at all: it contains zero
occurrences of "allele" and zero of "member". `refseq.md` contains "allele" twice (`:264-265`),
but both are the population-genetics sense — the wild-type "major allele present in the human
population" — not cis-allele membership.

So ferro's cross-zone ordering is a house rule. **Never cite the spec for it.** Argue it from the
underlying transcript coordinate, which is unambiguous, and say plainly that the `c.` spelling is
a presentation of that.

This is not theoretical. The sequence-changing and denotes-no-sequence classes the spec corpus
found both localise to the `c.72`/`c.*1` transition and to nothing else: the same flush
deletion-plus-insertion shape collapses correctly at all 22 other homopolymer runs in the test
transcript and mis-normalizes only where the pair straddles the CDS/3'UTR zone boundary. See
`the_cds_end_flush_pair_is_its_two_members_normalized_separately` and
`the_five_prime_boundary_masks_the_same_per_member_defect` in `tests/it/spec_corpus_regressions.rs`
— the second is the reminder that the 5' boundary is not a working case, only a masked one.

