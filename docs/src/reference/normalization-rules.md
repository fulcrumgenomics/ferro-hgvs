# Normalization rules

ferro's normalizer follows four rules about its output, and three about how it handles the gaps.

## The output contract

1. **Conformant.** Output follows the HGVS recommendations. *Absolute — never traded.*
   Scope: the `syntax.yaml` grammar, plus every prohibition, read from **prose force** rather than
   keyword casing — "not allowed", "not correct", "can not be used", "by definition", the
   `class="invalid"` markup, and the set `checklist.md` enumerates.

2. **Recommended form.** Where the spec prefers among conformant forms, ferro produces it.
   *Best effort.*
   Scope: "recommended", "preferred", lowercase "should", the 3' rule. A preference clause outranks
   maintainer judgment, but not rule 3: where it is not evaluable on what a normalizer holds, rule 3
   governs.

3. **Confluent.** Inputs denoting one variant produce one output. *Best effort.*
   Every rule — rule 2 included — is evaluated over the **resulting sequence**, never over the
   input's spelling. Reference context counts as sequence: transcript model, exon boundaries,
   reading frame, strand and topology are functions of the accession, not of the description.

4. **Deterministic.** Same input, same output. *Absolute.*
   Note that 4 does not imply 3 — a deterministic normalizer can be arbitrarily non-confluent.

## The procedure

5. **Where the spec is silent, ambiguous, or self-contradictory:** the issue is filed upstream
   **first** and cited, and only then does ferro ship a provisional choice.
    - **Self-contradictory** — two clauses that cannot both hold — is a **defect**. Every conforming
      tool must pick a side and none of them can be right, so filing is a bug report and is not
      optional.
    - **Silent** is merely incomplete. Ferro decides under rule 6 and violates nothing; filing is a
      feature request, worth making but never a reason to hold a release.

6. **Among multiple conformant forms:** the maintainers choose. **There are no user options for
   normalization form.** Error mode is an orthogonal axis and stays available. The 3'/5' shuffle
   direction was the one exception, and it is now removed from the public surface rather than
   excused: it is **not** orthogonal — it selects the frame every rule is evaluated in — so it was
   a user option for normalization form sitting inside this rule. ferro shifts 3', the only
   direction the HGVS recommendations describe, and no CLI flag, Python keyword or service config
   key selects otherwise. The 5' arm survives internally as a differential oracle over ferro's own
   test suite; an instrument is not a user option, and rules 2 and 3 are now claimed once rather
   than per direction.

7. **Disclosure.** Any change to these rules, and any different choice made under 5 or 6, is
   disclosed: in the changelog before v1, by a major version bump after. Output that *violates*
   rules 1-4 is a **bug**, not a disclosure.

## Why 2 and 3 are best effort, and 1 and 4 are not

Rules 1 and 4 are always achievable. Conformance is checkable against the spec text, and
determinism is a property of ferro's own code; nothing external can prevent us honouring them.

Rules 2 and 3 depend on the spec **determining an answer**, and sometimes it does not:

- **No preference exists.** The spec ranks substitution, deletion, inversion, duplication and
  insertion, but says nothing about competing `delins` forms.
- **Two preferences disagree.** `general.md` ranks duplication above insertion; `DNA/inversion.md`
  prefers an insertion for inverted copies. No single output satisfies both.
- **The same clause has two versions.** `general.md`'s current text and its forthcoming NOTE give
  opposite answers for variants separated by one nucleotide.
- **The variant's decomposition is not recoverable.** Recovering one means *choosing* an alignment,
  and the spec does not say which, so there is no derivable form to converge on. **Block length does
  not settle it**: an equal-length block can still carry a balanced del+ins pair, so its column
  correspondence need not be unique — `CAG -> AGA` is equal length with edit distance 2, not the
  position-wise 3. What decides whether a reference base is unchanged is whether **every minimal
  alignment** matches it, so the property to key on is edit distance against block length, never
  length alone. See `rulings[unchanged-is-read-over-every-minimal-alignment]`.
- **The preference keys on information a normalizer does not hold.** A "frequently occurring
  variant" (`RNA/delins.md:41`); a repeat "variable in the population" (`RNA/repeated.md:33`,
  `protein/repeated.md:22`); two variants "reported (or might occur) individually" (`DNA/delins.md:83`).
  The spec determines an answer; ferro cannot see what it keys on.

**"Best effort" is bounded by the spec's determinacy and by what a normalizer's inputs can decide —
not by ferro's implementation quality.** A failure of rule 2 or 3 caused by ferro's code is a
**bug** under rule 7; one caused by the spec not determining an answer triggers rule 5; one caused
by a clause keying on provenance or population data is a **declared deviation**, and a permanent
one — no upstream answer can put that information into a normalizer's hands.

Permanent is the only thing the last case adds, and it does **not** exempt the question from rule
5. The grain matters, because the two grains give opposite answers. Read **on its own**, such a
clause is not silent, ambiguous or self-contradictory — it determines an answer perfectly well, and
ferro is the one that cannot see what the answer keys on. Read **against the spec's own
re-derivation mandate** — `general.md:157-160`, which has a protein description derived by
comparing the variant and reference protein sequences and says knowledge of the underlying DNA
change "should not be used", with `general.md:13` extending the method to RNA — it is one half of a
pair that cannot both hold. That is rule 5's self-contradiction limb exactly, and the ruling ledger
says so in those words. What the carve-out states is that rule 5's escalation cannot *end*
the matter here: a choice made under rule 5 is provisional pending an upstream ruling, and no
ruling upstream makes provenance visible, so the deviation outlives the filing.

Where it is recorded is worth stating precisely, because the obvious pointer does not resolve.
`canonical-form-choice-when-both-legal` has **no `deviates_from` field**: it carries all four
clauses in its rationale as the recorded counter-evidence to re-derivation, opens that paragraph
with "The spec contradicts itself here", and adopts re-derivation over them deliberately. The one
of the four that is recorded *as* a deviation is `DNA/delins.md:83`, through the
`deviates_from: ["docs/recommendations/DNA/delins.md:79-84"]` on
`separation-is-a-property-of-the-spelling-not-of-the-variant` — the record whose ruling the rule 3
example below turns on.

## A worked example of reading force from prose

`DNA/duplication.md` says a variant that *can* be described as a duplication **must** be — but the
"must" is scoped by the preceding clause, which defines when a duplication *can* be used at all:
only when the additional copy is directly 3'-flanking the original. So the rule ranks the *label*
for one span; it does not require that a partition be chosen so as to produce a duplication.
Reading the force without the scope inverts the rule.

## What rule 3 excludes

"Never over the input's spelling" is narrower than it sounds. Most of a description is context.

| Carried by the description | Treatment |
|---|---|
| Accession, axis (`g.`/`c.`/`n.`), version | **Used** — it *is* the reference context |
| Which bases end up different | **Used** — this is the variant |
| Cis against trans (`[a;b]` against `[a];[b]`) | **Used** — different variants, not two spellings of one |
| Type label (`dup` against `ins`, `inv` against `delins`) | **Re-derived**, then ranked by `general.md` |
| How the edit set is cut into members | **Excluded** — a property of who wrote the string |
| Which copy in a run of identical residues a member names | **Excluded** — the 3' rule assigns this *"arbitrarily"* (`general.md:41`) |
| Repeat unit and phase, where several are equivalent | **Excluded** |

So three rows read **Excluded**, and they are three spellings of the one thing the input does not
get to decide: the **partition**, the run-position choice that feeds it, and — where several
unit-and-phase pairs describe one tract equally well — which pair a repeat member names.

`NC_000001.11:g.1001002_1001016` reads `ATGAGGGGCCACTGT`: a `GGGG` run at `1001006-1001009`, a `CC`
run at `1001010-1001011`, a lone `C` at `1001013`. Two spellings, one denoted sequence
(`ATGAGGGCATGT`), because `1001010` and `1001011` are both `C`:

```
g.[1001009del;1001010del;1001013del]     written gaps of 0 and 2
g.[1001009del;1001011del;1001013del]     written gaps of 1 and 1
```

`general.md:34`'s "two variants separated by one or more nucleotides should be described individually"
reads those gaps, so it answers twice for one variant. Rule 3 reads them off the partition ferro
derives instead: both give `g.[1001009_1001010del;1001013del]`. Rule 2 then keeps that over the
spanning `g.1001009_1001013delinsCA`, which merges across two unchanged nucleotides. Pinned in
[`tests/it/cis_confluence_adjudication.rs`](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/tests/it/cis_confluence_adjudication.rs).

## Known limitation

ferro cannot today guarantee that every input form is normalized according to these rules. They are
enforced intent, not a claim of current completeness.

Rule 7's disclosure mechanism — the `Representation-Change:` trailer and how it reaches the
changelog — is documented in
[CONTRIBUTING.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/CONTRIBUTING.md#declaring-a-representation-change).

## Where the individual decisions live

Rules 5 and 6 above say *how* a question is decided where the recommendations are silent,
ambiguous or self-contradictory. **What was decided, case by case**, is recorded separately, as
adjudication records — each naming the clauses in tension, which one governs, which is deviated
from, and why. Those records are published in full, with their clause quotes, in
[docs/NORMALIZATION_CONTRACT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_CONTRACT.md).

That document is generated from the records and gated against them, so it cannot drift; it
deliberately does not restate the seven rules above, which are stated only here.

The inverse index — for each stage of the normalizer, *which* record or clause governs it, and
which decisions are governed by **nothing** — is
[docs/NORMALIZATION_STAGE_AUDIT.md](https://github.com/fulcrumgenomics/ferro-hgvs/blob/main/docs/NORMALIZATION_STAGE_AUDIT.md). Read that one if what
you want to know is whether a behaviour you are looking at was chosen or merely happened.
