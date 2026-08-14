//! The two worked examples that explain the weight bound, side by side.
//!
//! # Why this file exists
//!
//! `canonicalize_from_sequence` carries a bound — `derived_columns >
//! changed_columns_of_edits(&edits)` in `merge.rs` — whose comment states it as
//! *"a canonicalization may re-partition and re-type the change; it may not
//! describe **more** change than the input already did."*
//!
//! It is the single largest confluence blocker measured in this codebase, and
//! also the thing standing between several inputs and an answer the spec
//! recommends. Whether it should exist at all is an open question. What kept
//! making that question expensive to ask is that its two faces were only ever
//! discussed one at a time, in prose, in documents outside the repository. This
//! file pins both, against the same reference, so the next reader gets the
//! contrast in one screen and does not have to reconstruct it.
//!
//! **Nothing here asserts that the bound is right or wrong.** Each row pins what
//! ferro does today, with the arithmetic that produced it and the clause that
//! bears on it. If the bound is removed or narrowed, these tests fail, and their
//! failure messages say which way to read the change.
//!
//! # The measure
//!
//! A description's weight is `sum over members of max(ref_len, alt_len)`. The
//! bound compares the weight of the *derived* partition against the weight of
//! the partition the *input* spelled, and refuses the derivation when it is
//! heavier — handing the variant back to the per-member pipeline, which never
//! re-aligns across members.
//!
//! Two consequences follow from the measure alone, and they are what the two
//! examples below show:
//!
//! * a re-derivation that marks *more columns changed* than the input is
//!   refused — which is the protection;
//! * **a span outweighs a split whenever the split keeps reference bases the
//!   span must cover**, so a merge across retained bases is refused however
//!   small the block — which is the cost.
//!
//! Stating the second one conditionally is not hedging; the unconditional form
//! is false. Writing `g` for the reference bases a split keeps but a span must
//! span, the two weights are `span = max(Σrᵢ, Σaᵢ) + g` and
//! `split = Σ max(rᵢ, aᵢ)`, so
//!
//! ```text
//! span − split  =  g  −  ( Σ max(rᵢ, aᵢ)  −  max(Σrᵢ, Σaᵢ) )
//!                  ^                    ^
//!            retained bases      slack: how much the members' own
//!            the span pays for   maxima exceed the block's
//! ```
//!
//! The refusal needs `g` to exceed that slack. Example 2 below has `g = 2` and
//! zero slack (`9 − 9`), so it is refused. A gap-free split does not survive the
//! comparison: two adjacent members weighing `max(3,1) + max(1,3) = 6` are
//! spanned by one member weighing `max(4,4) = 4`, and **that merge is accepted**.
//! So the bound does not refuse merges on principle — it refuses them on the
//! retained bases, which is a narrower and more interesting claim.
//!
//! # Example 1 — the protection
//!
//! [`INSERT_AND_DELETE`]. An insertion and a deletion whose lengths cancel, so
//! the block is net zero. The input spells two members, each weighing 1. A
//! re-derivation across them can only place one contiguous gap, so every base
//! between the two edits is offset by one and the derived form marks columns
//! changed that the input left alone. The bound refuses it, and the per-member
//! pipeline answers instead — which is why the output is the input.
//!
//! **This row is NOT protection — adjudicated 2026-08-10, by measurement.** The
//! columns below are the **position-wise** ones, which is what a net-zero block
//! lets you write down; they are not, by themselves, an answer to which bases
//! `general.md:34` counts as unchanged:
//!
//! ```text
//! pos      2  3  4  5  6  7  8  9 10 11 12
//! ref      T  G  C  A  C  C  A  G  T  C  A
//! observed T  A  G  C  A  C  C  A  G  T  C
//!             \___________/  ^  \___________/
//!              g.3_6           |   g.8_12
//!              delinsAGCA      |   delinsCAGTC
//!                        unchanged (pos 7)
//! ```
//!
//! Equal length, 11 against 11, with **9 of the 11 columns changed**.
//! `DNA/delins.md:16` types each consecutive run of changed columns as a
//! `delins`, and `:17` keeps the two runs separate — they are separated by one
//! unchanged base, and `:18`'s exception ("together affecting one amino acid")
//! cannot apply to a `g.` description, which has no amino acid. **So the derived
//! form is exactly what the spec's stated rules produce, and the input is a
//! spelling they never generate.**
//!
//! **This passage used to argue that from the wrong premise, and the premise is
//! settled elsewhere.** It read "for a net-zero block the column correspondence
//! is unique, so *which columns changed* is a fact rather than a choice", stated
//! as a general truth about equal-length blocks. It is not one:
//! `rulings[unchanged-is-read-over-every-minimal-alignment]` records `CAG ->
//! AGA` as an equal-length block with edit distance 2, where the position-wise
//! reading is *not* minimal and two of the three bases are unchanged. That
//! record is the only statement of the rule; do not restate it here. For **this**
//! block the position-wise reading happens to be the one the derivation takes,
//! which is all the worked example needs and all it now claims.
//!
//! The weights make the refusal vivid: the input spells `1 + 1 = 2`, the
//! derivation `max(4,4) + max(5,5) = 9`. `9 > 2`, so the bound refuses — even
//! though the heavier description is the conformant one.
//!
//! A survey of all twenty rows the bound's removal moves found **four** of this
//! shape — every one an equal-length block whose no-bound output the clauses
//! produce verbatim — **zero** rows the clauses contradict, and fifteen on which
//! `:15`/`:16`/`:17` have no defined input at all because the block is
//! unequal-length.
//!
//! **Corrected 2026-08-10, and worth stating because the wrong version was
//! persuasive.** This passage previously described the block as `TTTTTTTAAT` →
//! `ATTTTTTTAA` with changed columns `{0, 7, 9}`, predicted the no-bound output
//! to be a run of *substitutions*, and drew a conclusion about repetitive
//! context — that shifting by one base still matches at most positions **by
//! coincidence**, so an insertion-plus-deletion is column-wise indistinguishable
//! from a few scattered substitutions. Every one of those numbers belongs to a
//! **different fixture's** reference (`cis_junction_crossing_shift`'s poly-T
//! `TRACT`); the string `TTTTTTTAAT` does not occur in [`SEQ`] at all. On the
//! reference this file actually uses the region is *not* repetitive, so a
//! one-base shift changes almost everything rather than almost nothing, and the
//! derived form is two `delins` runs rather than scattered substitutions.
//!
//! The coincidence argument is still a real phenomenon — it is simply about
//! repetitive contexts, which this fixture does not build, so this row is not
//! evidence for it. Whether such an output is *meaningful* — one of the four
//! design values at `background/basics.md:38`, alongside stable, memorable and
//! unequivocal, in a list from which minimality is absent — remains a question
//! for the committee, not for a bound in this crate.
//!
//! # Example 2 — the cost
//!
//! [`SPANNING_DELINS`] and [`SEPARATED_MEMBERS`] are two spellings of one
//! variant — asserted below against an applier independent of the normalizer,
//! not against `EquivalenceChecker`, which is documented elsewhere in this suite
//! as circular for this purpose.
//!
//! ```text
//! pos        29  30  31  32  33
//! ref         C   A   C   G   T          (asserted below)
//! observed    A   A   C   ACATACTG
//!             ^   \_____/   ^
//!          changed   kept  changed
//! ```
//!
//! The `AC` at 30–31 survives in the split spelling **only because the payload
//! happens to carry `AC` at that offset** — with five reference bases becoming
//! eleven there is no column correspondence to force it. That is precisely the
//! construction `DNA/delins.md:46` builds, and `:47` answers it: *"**The
//! "delins" format is recommended**"*.
//!
//! The weights:
//!
//! ```text
//! split   29C>A                 max(1,1) =  1
//!         32_33delinsACATACTG   max(2,8) =  8      total  9
//! span    29_33delinsAACACATACTG max(5,11) = 11    total 11
//! ```
//!
//! `11 > 9`, so the bound refuses the span. Both spellings are therefore fixed
//! points, and one variant has two normalized representations — which is the
//! non-confluence reported as #1421.
//!
//! # Assert-then-flip
//!
//! Every expectation here is today's behaviour, asserted rather than
//! `#[ignore]`d, so the day it moves this file has to be read. The failure
//! messages name the form each row is expected to move *to*.

use crate::common::cis_apply_oracle as oracle;

/// The synthetic reference, verbatim from the #1235 reproductions so this file
/// and those reports stay comparable.
const SEQ: &str = concat!(
    "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
    "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA"
);

/// Example 1 — an insertion and a deletion whose lengths cancel.
const INSERT_AND_DELETE: &str = "TEMPLATE:g.[2_3insA;12del]";

/// Example 2, the spanning spelling. The form `delins.md:47` recommends.
const SPANNING_DELINS: &str = "TEMPLATE:g.29_33delinsAACACATACTG";

/// Example 2, the separated spelling. Cheaper by the bound's measure, because
/// the payload's `AC` coincides with the reference at 30–31.
const SEPARATED_MEMBERS: &str = "TEMPLATE:g.[29C>A;32_33delinsACATACTG]";

/// The reference bases every argument in the module docs rests on.
///
/// Asserted rather than quoted: a reference that does not carry them would make
/// the rows below pin something other than what this file says they pin.
#[test]
fn the_reference_is_the_one_the_examples_were_measured_on() {
    assert_eq!(
        &SEQ[28..33],
        "CACGT",
        "TEMPLATE:g.29_33 must read CACGT for the coincidence argument to hold"
    );
    assert_eq!(SEQ.len(), 125, "the reproductions used a 125 nt TEMPLATE");
}

/// Example 1 is preserved: the bound refuses a derivation heavier than the input.
///
/// **Assert-then-flip.** If the bound is removed or narrowed this becomes the
/// re-derived form `TEMPLATE:g.[3_6delinsAGCA;8_12delinsCAGTC]` — **measured**,
/// not predicted (see the module docs' 2026-08-10 correction; an earlier version
/// of this comment predicted a substitution run, from another fixture's
/// reference).
#[test]
fn a_net_zero_insert_delete_pair_is_left_to_the_per_member_pipeline() {
    let out = oracle::normalize(SEQ, INSERT_AND_DELETE);
    assert_eq!(
        out, INSERT_AND_DELETE,
        "the weight bound refuses a re-derivation that marks more columns changed \
         than the input did, so the per-member pipeline answers and the input is \
         preserved. If this just failed with \
         `TEMPLATE:g.[3_6delinsAGCA;8_12delinsCAGTC]`, the bound moved and that is \
         the expected new form: the block is equal-length (11/11) with 9 columns \
         changed, so `DNA/delins.md:16` types each run as a delins and `:17` keeps \
         them separate. Flip this expectation and record the adjudication. Any \
         OTHER form is a finding — re-derive the columns before accepting it"
    );
}

/// The two Example 2 spellings denote one sequence — established independently.
///
/// Uses the module's applier, which converts each member to its SPDI triple and
/// splices the reference, rather than `EquivalenceChecker`. That distinction is
/// deliberate: a sequence claim checked with the normalizer's own equivalence
/// machinery cannot separate "these agree" from "these agree because the same
/// code decided both".
#[test]
fn the_two_spellings_denote_the_same_sequence() {
    let span = oracle::apply(SEQ, SPANNING_DELINS);
    let split = oracle::apply(SEQ, SEPARATED_MEMBERS);
    assert!(
        span.is_some(),
        "the applier must resolve {SPANNING_DELINS}; if it declines, this file \
         proves nothing"
    );
    assert_eq!(
        span, split,
        "{SPANNING_DELINS} and {SEPARATED_MEMBERS} must denote one sequence — \
         that is what makes their two normalized forms a confluence defect \
         rather than two different variants"
    );
}

/// Example 2 CONVERGES (#1835) — one variant, one normalized string.
///
/// # The assert-then-flip fired, and it flipped the OTHER WAY
///
/// This row used to pin a non-confluence and carried an instruction for closing
/// it: "confirm they agree on [`SPANNING_DELINS`] (`DNA/delins.md:47`, and the
/// decided ruling `delins-merge-vs-individual-gap-two-or-more`)". **Do not follow
/// that instruction.** The two spellings now agree, so the defect is indeed
/// fixed — but they agree on the SPLIT, and two decided rulings postdating that
/// sentence say that is the correct side. Each is sufficient on its own:
///
/// * **`delins-payload-coincidence-carve-out-is-coding-dna-scoped`** scopes
///   `delins.md:47` to the coding DNA axis. `TEMPLATE:g.` is genomic, so `:47`
///   never reaches this row and `general.md:34` governs unqualified — members
///   separated by unchanged nucleotides are described individually. `:47`'s own
///   stated ground is preventing "incorrect predictions for the consequences on
///   protein level", which has nothing to bite on where no protein is predicted.
/// * **`delins-merge-vs-individual-gap-two-or-more`** — the very record the old
///   instruction cited — was later scoped BY DIRECTION to the net-deletion case
///   and "does **not** reach net insertions, where the split form stays
///   canonical". This row replaces a 5 nt span with an 11 nt payload: a net
///   insertion, so the record does not reach it either.
///
/// So the old sentence named the right authorities and the wrong side of both,
/// because both were scoped after it was written. Left in place above rather than
/// deleted, because the mistake is instructive: a failure message that tells you
/// which way to re-pin is only as current as the ledger it was copied from.
///
/// **Mutalyzer merging it is not authority here.** `adjudication-precedence-order`
/// records why Mutalyzer is not a spec oracle — it minimizes a weighted
/// description length with constants dated 2014 and has no separation rule at
/// all, so a separation disagreement is two objectives meeting rather than
/// evidence.
///
/// The converged form is not the input's spelling either, which is what
/// `canonical-form-choice-when-both-legal` requires: the answer is derived from
/// the resulting sequence, and neither authored spelling survives it.
#[test]
fn one_variant_normalizes_to_two_strings_because_a_span_outweighs_its_split() {
    /// What both spellings now reach. Neither authored form.
    const CONVERGED: &str = "TEMPLATE:g.[28_29insAA;32G>A;33_34insACTG]";

    let from_span = oracle::normalize(SEQ, SPANNING_DELINS);
    let from_split = oracle::normalize(SEQ, SEPARATED_MEMBERS);

    assert_eq!(
        from_span, CONVERGED,
        "the spanning spelling is re-derived rather than preserved"
    );
    assert_eq!(
        from_split, CONVERGED,
        "and the separated spelling reaches the same form"
    );
    assert_eq!(
        from_span, from_split,
        "the two spellings of one variant must converge. If they have diverged \
         again this is a confluence regression, not a re-pin: see \
         `confluence-gate-is-apply-equality-on-every-determined-axis`"
    );
}

/// The bound's measure, stated as arithmetic rather than as prose.
///
/// Pinned because every explanation of this defect has had to re-derive these
/// four numbers, and because the retained bases are what make the cost
/// structural: the span must pay for the `AC` at 30–31 that the split keeps, and
/// no threshold on the bound can be tuned to let that through without also
/// admitting the derivations Example 1 exists to refuse.
///
/// **Scoped deliberately.** An earlier revision of this test was named
/// `a_span_always_outweighs_a_split_of_the_same_block` and its message claimed
/// "any input-relative bound refuses every merge of this shape". That is false:
/// with `g` retained bases, `span − split = g − (Σ max(rᵢ,aᵢ) − max(Σrᵢ,Σaᵢ))`,
/// so a gap-free split can be *cheaper* than its span and the merge is accepted
/// — `max(3,1) + max(1,3) = 6` against a span of `max(4,4) = 4`. The module docs
/// carry the derivation. What this row pins is the retained-gap case, which is
/// the one the fixture builds and the one #1421 reports.
#[test]
fn a_span_outweighs_a_split_that_keeps_reference_bases() {
    // max(ref_len, alt_len) per member.
    let split_weight = 1 + 8; // 29C>A ; 32_33delinsACATACTG
    let span_weight = 11; // 29_33delinsAACACATACTG

    assert_eq!(split_weight, 9);
    assert_eq!(span_weight, 11);

    // The gap is the whole difference: 30-31 are kept by the split and spanned
    // by the merge, and the members' maxima already sum to the block's own
    // (9 == max(3, 9)), leaving no slack to absorb it.
    let retained_bases = 2; // 30-31, the `AC` the payload coincides with
    assert_eq!(
        span_weight - split_weight,
        retained_bases,
        "the refusal is priced entirely by the bases the split keeps and the \
         span must cover; if this difference ever stops equalling the gap, the \
         measure has changed and the module docs need re-deriving"
    );
    assert!(
        span_weight > split_weight,
        "a merge across retained reference bases is refused: the span pays for \
         the kept bases while the split does not, so no threshold separates \
         this from the derivations Example 1 refuses. Note this is the \
         retained-gap case specifically — a gap-free split can be heavier than \
         its span, and that merge is accepted"
    );
}
