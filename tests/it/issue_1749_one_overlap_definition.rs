//! #1749 — "do two cis members overlap?" was decided in four places, on three
//! different geometries, and they disagreed.
//!
//! The correct notion was already established — key on what a member **writes**,
//! not on what it **reads** (#1411, #1437, #1448) — but it lived as prose in a
//! comment plus two partial implementations that shared no code, so `repeat`
//! ended up with three different answers inside one file and nothing could see
//! the disagreement.
//!
//! The geometry now lives once, in `normalize::footprint::WriteFootprint`, and
//! every consumer derives its question from it. This file pins the **end-to-end**
//! consequences; `normalize::overlap`'s own `FOOTPRINT_MATRIX` pins the detector
//! agreement at unit level, and `normalize::footprint`'s tests pin the geometry
//! itself.
//!
//! ## Measured on real hg38, `origin/main` at `b49c5ce3`
//!
//! `ferro normalize --reference <prepared dir> --error-mode strict`:
//!
//! ```text
//! NC_000017.11:g.[43045721dup;43045721_43045722dup]          -> g.43045721_43045722insAGA
//! NC_000017.11:g.[43045721_43045724del;43045722_43045723del] -> REFUSED, W5002
//! NC_000017.11:g.[43045721_43045725del;43045723_43045727dup] -> REFUSED at PARSE (E3006)
//! ```
//!
//! Those coordinates need a prepared reference, so the same three geometries are
//! pinned below against the committed fixture, where they run unconditionally. A
//! reference-aware test that skips when the manifest is absent reports as PASS
//! rather than SKIP, which would make this file look like coverage it is not.
//!
//! ## The fourth decider is deliberately untouched
//!
//! Row 3 above is refused by `HgvsVariant::detect_self_cancelling_pair` (E3006),
//! which keys on **read** spans. Under the write-footprint geometry the `del`
//! writes `[43045721, 43045725]` and the `dup` writes at the junction
//! `43045727|43045728`, so the pair is disjoint and denotes one sequence.
//!
//! It is still refused, and that is correct: `general.md:58` names this
//! edit-type pair explicitly and its own worked example
//! (`NM_004006.2:c.[762_768del;767_774dup]`) overlaps at 767-768 without
//! cancelling — so the clause keys on read spans by construction. Precedence is
//! spec-explicit over our judgement, so the clause wins over the geometry. What
//! #1749 changes is that this is now *stated* as a scoped spec prohibition
//! rather than sitting unremarked as a fourth opinion on the word "overlap".

use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer};

/// 39 bp, with an `AA` at 10-11 and a `GTGTGT` tandem tract at 20-25.
/// (9 T + `AA` + 8 T + `GTGTGT` + 14 T.)
const SEQUENCE: &str = "TTTTTTTTTAATTTTTTTTGTGTGTTTTTTTTTTTTTTT";

fn provider() -> MockProvider {
    let mut provider = MockProvider::new();
    provider.add_genomic_sequence("NC_TEST.1", SEQUENCE.to_string());
    provider
}

fn strict_accepts(input: &str) -> bool {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::with_config(provider(), NormalizeConfig::strict())
        .normalize(&variant)
        .is_ok()
}

/// What strict mode settles `input` on. Panics if it refuses, so a caller that
/// means to pin a *form* cannot silently degrade into pinning a refusal.
fn strict_settles_on(input: &str) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::with_config(provider(), NormalizeConfig::strict())
        .normalize(&variant)
        .expect("strict must accept this input")
        .to_string()
}

/// Two `dup`s over overlapping read spans compose, because they write at two
/// different junctions — the row Mutalyzer refuses ("Variants overlap;
/// overlapping variant locations are not allowed") and ferro accepts on purpose.
#[test]
fn two_duplications_over_overlapping_read_spans_still_merge() {
    assert!(
        strict_accepts("NC_TEST.1:g.[10dup;10_11dup]"),
        "the two `dup`s write at junctions 10 and 11, so their overlapping read \
         spans are not a collision and strict must accept"
    );
}

/// Two deletions whose spans intersect are still refused — the #1448 half, which
/// the unification must not loosen.
#[test]
fn two_intersecting_deletions_are_still_refused() {
    for input in [
        "NC_TEST.1:g.[10_14del;11_12del]",
        "NC_TEST.1:g.[10_14del;12_16del]",
        "NC_TEST.1:g.[10_14del;10_16del]",
    ] {
        assert!(
            !strict_accepts(input),
            "`{input}`: two deletions rewriting an intersecting range denote no \
             single sequence, so strict must refuse"
        );
    }
}

/// `general.md:58`'s `del`+`dup` prohibition is enforced at parse and is
/// **not** derived from the write-footprint geometry.
///
/// Pinned as its own row so that a later change to the geometry cannot silently
/// take this refusal with it: the two are independent, and the clause is the
/// authority for this one.
#[test]
fn the_spec_del_dup_prohibition_is_still_enforced_at_parse() {
    let refused = parse_hgvs("NC_TEST.1:g.[10_14del;12_16dup]");
    assert!(
        refused.is_err(),
        "`general.md:58` names a `del` combined with a `dup` of part of the same \
         sequence and disallows it, keying on read spans by its own worked \
         example — so this is refused at parse whatever the write footprints say"
    );
}

/// A `repeat` is no longer answered three ways.
///
/// The three passes used to disagree: the coincident-bounds pass gave it its
/// whole span, the intersection pass dropped it entirely, and the junction pass
/// made it junction-*only*. Which answer you got depended only on which pass
/// happened to look.
///
/// It is now decided by arithmetic the description already carries — `unit x
/// count` against the tract length — so it needs no provider and no guess.
mod a_repeat_is_decided_by_its_own_arithmetic {
    use super::*;

    /// **Grows.** `GT[4]` over the 6-base tract 20-25 asserts those bases
    /// unchanged and lands one extra copy at the junction 3' of 25 — a `dup`'s
    /// geometry, which is what keeps the two spellings of one shape in
    /// agreement (#1437).
    ///
    /// The second row asserted a **refusal** until the dup-plus-insertion
    /// ruling, on the ground that "both write at junction 25, with nothing in
    /// the description saying which goes first". The spec says otherwise, and
    /// this is the clause that settles it — `duplication.md:90` publishes a
    /// pair at one junction and orders it in prose:
    ///
    /// ```text
    /// NC_000001.11(NM_206933.2):c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]
    ///   a duplication of the sequence from c.675-542 to c.1211-703,
    ///   followed by the insertion of the sequence GTAAA.
    /// ```
    ///
    /// The insertion's junction *is* the duplication's own write junction, and
    /// "followed by" supplies the composition order — which is also the genomic
    /// order `alleles.md:118` already requires members to be listed in. So the
    /// pair is ordered, not ambiguous.
    ///
    /// Keeping the refusal here while the `dup` spelling of the same shape is
    /// accepted would be precisely the #1437 disagreement this file exists to
    /// close, so the row is flipped rather than carved out.
    ///
    /// It pins the **composed form**, not merely acceptance: ordering was the
    /// whole of the original objection, so a test that only checked the input
    /// were accepted would no longer be about anything. `insGTA` is the extra
    /// `GT` followed by the `A`; the other order would spell `insAGT`.
    #[test]
    fn a_growing_repeat_writes_only_at_its_junction() {
        assert!(
            strict_accepts("NC_TEST.1:g.[20_25GT[4];21_22insA]"),
            "the repeat grows, so it writes at junction 25 and the insertion at \
             21 collides with nothing"
        );
        assert_eq!(
            strict_settles_on("NC_TEST.1:g.[20_25GT[4];25_26insA]"),
            "NC_TEST.1:g.25_26insGTA",
            "both write at junction 25, and `duplication.md:90` orders that pair \
             by the members' genomic order — the repeat's extra copy first, then \
             the insertion"
        );
    }

    /// **Shrinks.** `GT[1]` over the same tract keeps 20-21 and removes 22-25,
    /// so its footprint is that trailing range — which a sibling deletion can
    /// collide with, and which the intersection pass used to be blind to
    /// because it dropped `repeat` outright.
    #[test]
    fn a_shrinking_repeat_collides_over_the_bases_it_removes() {
        assert!(
            !strict_accepts("NC_TEST.1:g.[20_25GT[1];23_27del]"),
            "the repeat removes 22-25 and the deletion claims 23-27, so the two \
             rewrite an intersecting range"
        );
        // The kept prefix is untouched reference, so a sibling there collides
        // with nothing the repeat writes. The sibling must be one that cannot
        // 3'-shift, though: inside a tandem tract a `20_21del` is the same
        // variant as a deletion of the tract's last unit, so it shifts onto the
        // very bases the repeat removes and the refusal there is correct.
        assert!(
            strict_accepts("NC_TEST.1:g.[20_25GT[1];20G>C]"),
            "20 is in the KEPT prefix and a substitution does not shift, so it \
             collides with nothing the repeat writes"
        );
    }

    /// The two spellings of one shape get one verdict, whichever the per-member
    /// pipeline picked. This is #1437's argument and the reason the unification
    /// could not simply give `repeat` its whole span.
    ///
    /// **The pairs must denote the same edit, or the test is about nothing.**
    /// They read `A[4]` against `10_11dup` because the tract at 10-11 is `AA`,
    /// so a fourth `A`-unit is exactly the copy the `dup` adds. An earlier
    /// version paired `A[3]`, which adds one base where the `dup` adds two —
    /// invisible while the assertion compared only verdicts.
    #[test]
    fn the_dup_and_repeat_spellings_of_one_shape_agree() {
        for (repeat, dup) in [
            (
                "NC_TEST.1:g.[10_11A[4];10_11insT]",
                "NC_TEST.1:g.[10_11dup;10_11insT]",
            ),
            (
                "NC_TEST.1:g.[10_11A[4];11_12insT]",
                "NC_TEST.1:g.[10_11dup;11_12insT]",
            ),
        ] {
            assert_eq!(
                strict_accepts(repeat),
                strict_accepts(dup),
                "`{repeat}` and `{dup}` are one shape in two spellings; the \
                 per-member pipeline picks between them BY AXIS, so a detector \
                 that answered differently would be sensitive to a choice the \
                 author never made"
            );
        }
    }

    /// Agreeing on the *verdict* is not the whole contract once the verdict is
    /// ACCEPT — two spellings that both refuse agree trivially, and two that
    /// both accept can still settle on different forms. So the shared-junction
    /// pair is pinned on the form as well.
    ///
    /// `insAAT` is the duplicated `AA` followed by the inserted `T`, which is
    /// `duplication.md:90`'s "followed by" order. Both spellings reach it.
    #[test]
    fn the_two_spellings_settle_on_one_form_at_a_shared_junction() {
        assert_eq!(
            strict_settles_on("NC_TEST.1:g.[10_11A[4];11_12insT]"),
            strict_settles_on("NC_TEST.1:g.[10_11dup;11_12insT]"),
            "one shape, two spellings, one junction — the forms must agree"
        );
        assert_eq!(
            strict_settles_on("NC_TEST.1:g.[10_11dup;11_12insT]"),
            "NC_TEST.1:g.11_12insAAT",
            "the duplicated `AA` followed by the inserted `T`"
        );
    }

    /// **Known divergence, #1996.** When the sibling insertion is *interior* to
    /// the duplicated span rather than at its write junction, the two spellings
    /// are both accepted and settle on different forms: the `dup` merges into
    /// one member, the `repeat` stays as two.
    ///
    /// Pinned as measured rather than left to the agreement test above, which
    /// compares verdicts and so cannot see it. This divergence is **older than
    /// the dup-plus-insertion ruling** — it reproduces byte-for-byte on
    /// `origin/main` at `de8532ee`, where both rows are likewise accepted.
    /// Recorded here so that fixing #1996 has a place to land, and so that the
    /// difference is not mistaken for fallout of the junction ruling.
    #[test]
    fn an_interior_insertion_still_splits_the_two_spellings() {
        assert_eq!(
            strict_settles_on("NC_TEST.1:g.[10_11A[4];10_11insT]"),
            "NC_TEST.1:g.[10_11insT;10_11A[4]]",
            "#1996: the repeat spelling stays as two members"
        );
        assert_eq!(
            strict_settles_on("NC_TEST.1:g.[10_11dup;10_11insT]"),
            "NC_TEST.1:g.10_11insTAA",
            "#1996: the dup spelling merges into one"
        );
    }
}

/// An insertion flush against either edge of a deletion composes uniquely, and
/// says so **in either member order**.
///
/// `CLAUDE.md` records this as a 233-fire false-positive class of the
/// denoted-sequence oracle. The cause was that `spdi::apply::triples_are_disjoint`
/// carried a running 3'-most `reach` compared against each triple's position,
/// which read a tie between a zero-width triple and a span triple differently
/// depending on which the **stable** descending sort happened to visit first.
#[test]
fn an_insertion_flush_against_a_deletion_composes_in_either_order() {
    for input in [
        "NC_TEST.1:g.[10_14del;14_15insA]",
        "NC_TEST.1:g.[14_15insA;10_14del]",
        "NC_TEST.1:g.[10_14del;9_10insA]",
        "NC_TEST.1:g.[9_10insA;10_14del]",
    ] {
        assert!(
            strict_accepts(input),
            "`{input}`: an insertion at the edge of a deletion claims no base \
             the deletion claims, so it is well defined — and the verdict must \
             not depend on which member was written first"
        );
    }
}

/// An insertion *interior* to a deletion is still accepted (#1406) and an
/// insertion interior to a span that KEEPS its bases is still refused.
///
/// The pair is the discriminator: a rule that simply stopped looking at
/// junctions would pass the first and break the second.
#[test]
fn the_deletion_carve_out_still_has_a_boundary() {
    assert!(
        strict_accepts("NC_TEST.1:g.[10_14del;12_13insA]"),
        "the deletion removes every base it spans, so an interior junction has \
         nothing left to be positioned against (#1406)"
    );
    assert!(
        !strict_accepts("NC_TEST.1:g.[10_14delinsGG;12_13insA]"),
        "a `delins` payload survives, so an interior insertion has a position \
         relative to it and the pair is genuinely ambiguous"
    );
    assert!(
        !strict_accepts("NC_TEST.1:g.[10_14inv;12_13insA]"),
        "an inversion keeps the bases it spans, so the interior junction stays \
         meaningful"
    );
}
