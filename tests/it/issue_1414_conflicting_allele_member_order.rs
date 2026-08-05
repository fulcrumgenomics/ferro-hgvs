//! #1414 — a conflicting cis allele skipped the genomic-order sort even when
//! every one of its members had been rewritten.
//!
//! `normalize_allele` gated the sort on whether an `OverlapConflict` was
//! reported. The skip itself is right and long-standing: an overlap-conflicting
//! allele is emitted **verbatim** (#395), and reordering an allele you are
//! handing back untouched would contradict that.
//!
//! The premise is what failed. *Conflict was reported* and *allele was emitted
//! verbatim* are different predicates, and they come apart exactly when the
//! repair passes succeed: `11_12inv` cancels to `11_12=` because the span is its
//! own reverse complement, `11_12insAA` respells as `10_11A[4]`, and the
//! pipeline then declines to order the two on the strength of a contract about
//! *not* rewriting them.
//!
//! Two defects followed, and the first is the more serious:
//!
//! 1. **Out-of-genomic-order members** — a direct violation of #1235's
//!    criterion 2 ("no normalized cis allele contains overlapping or
//!    out-of-order members"), on a shape its confluence proptest cannot
//!    generate (its indel model draws no inversion beside a co-located
//!    insertion).
//! 2. **Non-idempotence**, for the subset whose reordered output no longer
//!    conflicts: pass 2 then takes the *other* branch of the same gate and
//!    sorts. `FERRO_ASSERT_IDEMPOTENT` fires on those.
//!
//! **#1423 is not this fix, and did not make this unreachable.** It repaired the
//! single input the issue cited (`g.[11_12inv;11_12insAA]`, now the single
//! member `g.10_11A[4]`) by dropping an identity a sibling *overlaps* — and a
//! zero-width junction insertion overlaps nothing, so the whole `[inv;ins]`
//! family survived it. Measured on `origin/main` at `94817bf6` over four
//! sequences and eleven member shapes: **41** inputs in four families still
//! reached the divergence.
//!
//! The gate now asks the question the contract is actually about — was this
//! allele handed back untouched? — so an allele whose members were rewritten is
//! ordered like any other, and one authored exactly as it settles is not.

use crate::common::cis_apply_oracle::{normalize, normalize_in};
use ferro_hgvs::ShuffleDirection;

/// The `A`/`T` corpus sequence the sweep draws, where the `[inv;ins]` and
/// `[dup;inv]` families both settle out of order on `main`.
const AT_SEQUENCE: &str = "TTTTTTTTTAATATATTTTA";

/// The issue's own reference: 18 bases, `10` and `11` are `A`, `12` is `T`.
const CG_SEQUENCE: &str = "CGCGCGCGCAATCGCGCG";

/// The headline: a rewritten conflicting allele comes back in genomic order.
///
/// **Rewritten, not necessarily wholly rewritten.** In both rows below the
/// inversion survives *as an inversion* and only its sibling moves — the
/// duplication respells as an insertion at its own junction, the insertion
/// 3'-shifts within the inverted span. That is enough: the trigger is that the
/// allele was not emitted verbatim, so there is no authored form left for #395
/// to preserve, and one moved member is as disqualifying as all of them.
///
/// (An earlier revision of this comment said "every member here is rewritten —
/// the inversion cancels to an identity". That described the rows this test
/// carried before #1445, which are now preserved rather than repaired and have
/// moved to `issue_1406_lenient_output_keeps_the_conflict`.)
///
/// Pinned as exact strings rather than checked for orderedness. A "members are
/// sorted" predicate is satisfied by a refusal, by a dropped member, and by any
/// other ordered spelling; what must come back is these members, in this order.
///
/// **Re-blessed against #1406/#1445.** Two of the rows this test originally
/// carried — `g.[9_10insA;9_10inv]` and `g.[11_12dup;12_13inv]` — are no longer
/// rewritten at all. #1445 hands a laundered conflict back as authored, because
/// the independent applier declines such an input (`ApplyFailure::Overlapping`)
/// and so there is no denoted sequence for a repair to preserve. They are
/// therefore not examples of "a rewritten conflicting allele" any more, and
/// `issue_1406_lenient_output_keeps_the_conflict` is where they are now pinned.
///
/// What still reaches this gate is an allele whose members are rewritten and
/// which **still conflicts afterwards**, so #1445's restore does not fire. Both
/// rows below are that shape, and the gate is load-bearing for them: measured by
/// reverting it to `!has_overlap_conflict` on top of #1445, this test and
/// `issue_1235_transcript_axes::overlap_conflicting_allele_is_not_canonicalized`
/// both fail.
#[test]
fn a_rewritten_conflicting_allele_is_returned_in_genomic_order() {
    for (input, expected) in [
        // `[inv;dup]` — the duplication respells as an insertion at its own
        // junction; the inversion survives as an inversion.
        (
            "TEMPLATE:g.[11_13inv;11dup]",
            "TEMPLATE:g.[11_12insA;11_13inv]",
        ),
        // `[ins;inv]` where the insertion 3'-shifts inside the inversion's span
        // rather than out of it, so the pair still conflicts and is still
        // ordered by this gate.
        (
            "TEMPLATE:g.[12_13insA;12_14inv]",
            "TEMPLATE:g.[12_14inv;13_14insA]",
        ),
    ] {
        assert_eq!(
            normalize(AT_SEQUENCE, input),
            expected,
            "`{input}` rewrote every member, so the genomic-order sort must run"
        );
    }
}

/// …and the same on the 5' axis, which carries the larger share of the moved
/// rows (1178 of 2037 measured). The gate is direction-agnostic, so a fix that
/// only reached the default direction would be half a fix.
#[test]
fn the_five_prime_axis_gets_the_same_answer() {
    let input = "TEMPLATE:g.[10_11inv;10dup]";
    let five = normalize_in(AT_SEQUENCE, input, ShuffleDirection::FivePrime);
    assert_eq!(
        five, "TEMPLATE:g.[10_11insA;10_11inv]",
        "the 5' axis must order a rewritten conflicting allele too"
    );
    assert_eq!(
        five,
        normalize_in(AT_SEQUENCE, &five, ShuffleDirection::FivePrime),
        "the 5' answer must be a fixed point too"
    );
}

/// The idempotency half, on the subset where the reordered output no longer
/// conflicts — the shape that fires `FERRO_ASSERT_IDEMPOTENT`.
///
/// 209 of the 2037 moved rows were non-idempotent on `main` this way; the rest
/// were stable but out of order. Both are criterion-2 violations, and this pins
/// the one that is also a live oracle failure.
#[test]
fn the_rewritten_output_is_a_fixed_point_in_one_pass() {
    // The same two rows the headline test uses, and for the same reason: since
    // #1445 the inputs this test originally named are handed back as authored,
    // so they are fixed points trivially and would no longer exercise the
    // second-pass branch this test exists to reach.
    for input in [
        "TEMPLATE:g.[11_13inv;11dup]",
        "TEMPLATE:g.[12_13insA;12_14inv]",
    ] {
        let once = normalize(AT_SEQUENCE, input);
        let twice = normalize(AT_SEQUENCE, &once);
        assert_eq!(
            once, twice,
            "`{input}` must settle in one pass; got `{once}` then `{twice}`"
        );
    }
}

/// The contract the gate exists to honour is unchanged: an allele authored
/// exactly as its members settle is still handed back verbatim, in its authored
/// order, however un-genomic that order is.
///
/// This is the discriminating half of the fix. Replacing the conflict test with
/// an unconditional sort would pass every assertion above and break #395 here.
#[test]
fn an_allele_emitted_verbatim_keeps_its_authored_order() {
    // `main`'s own out-of-order output for `g.[12_13insA;12_14inv]`. Both
    // members are already settled — the insertion sits at its 3'-most junction,
    // the inversion at its own span — so nothing is rewritten and the verbatim
    // contract genuinely applies, un-genomic order and all.
    let verbatim = "TEMPLATE:g.[13_14insA;12_14inv]";
    assert_eq!(
        normalize(AT_SEQUENCE, verbatim),
        verbatim,
        "an allele whose members are already settled must be returned untouched"
    );
    // Deliberately the reverse of genomic order, and it must survive.
    assert_ne!(
        members_of(verbatim),
        sorted_members_of(verbatim),
        "the fixture must actually be out of genomic order, or it guards nothing"
    );
}

/// #1423's own input stays fixed: the fix must not disturb a single-member
/// output, which has no ordering question at all.
#[test]
fn the_single_member_output_of_1423_is_untouched() {
    assert_eq!(
        normalize(CG_SEQUENCE, "TEMPLATE:g.[11_12inv;11_12insAA]"),
        "TEMPLATE:g.10_11A[4]",
        "#1423's collapse to one member must survive this change"
    );
}

/// A disjoint, non-conflicting allele is still sorted, as it was before. The
/// gate's new term only ever *adds* sorting, and this pins that the ordinary
/// path did not move.
#[test]
fn a_disjoint_allele_is_still_sorted() {
    let out = normalize(CG_SEQUENCE, "TEMPLATE:g.[12T>G;10A>G]");
    assert_eq!(
        out, "TEMPLATE:g.[10A>G;12T>G]",
        "a disjoint allele must still be rendered in genomic order"
    );
}

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/// The member list of a bracketed allele, in the order it is rendered.
fn members_of(description: &str) -> Vec<String> {
    description
        .rsplit_once('[')
        .and_then(|(_, rest)| rest.strip_suffix(']'))
        .map(|inner| inner.split(';').map(str::to_string).collect())
        .unwrap_or_default()
}

/// The same members in genomic order, by leading coordinate then by text.
///
/// A presentational sort for the assertions above, not a re-implementation of
/// `cis_member_order_key` — the exact-string pins are what judge the order; this
/// only has to be right enough to tell "ordered" from "reversed".
fn sorted_members_of(description: &str) -> Vec<String> {
    let mut members = members_of(description);
    members.sort_by_key(|member| {
        let digits: String = member.chars().take_while(|c| c.is_ascii_digit()).collect();
        (digits.parse::<u32>().unwrap_or(u32::MAX), member.clone())
    });
    members
}
