//! #1831 — the applier's overlap verdict must not depend on member order.
//!
//! An insertion at interbase `p` and a `sub`/`del` starting at `p` share an SPDI
//! position. `apply_triples` sorted on position alone with a **stable** sort, so
//! the two kept whatever order the author wrote them in, and
//! `triples_are_disjoint`'s running-reach walk then called the pair overlapping
//! or disjoint accordingly. `g.[9_10insT;10G>A]` was refused; `g.[10G>A;9_10insT]`
//! — the same variant — applied.
//!
//! Both halves of that are pinned here, and pinned as **equality between the two
//! orderings** rather than as "the first one now applies": the property is that
//! member order carries no meaning (#1098), so a future change that made the
//! *second* ordering refuse would be just as wrong and must fail this too.
//!
//! # What must NOT change
//!
//! The fix is a sort key, not a loosening of the disjointness rule, and the two
//! shapes below are the discriminating cases. A tempting wider fix — dropping or
//! weakening `triples_are_disjoint` — passes every assertion above and fails
//! these:
//!
//! * members that genuinely claim the same base still refuse;
//! * two pure insertions at one interbase still refuse, because their order is
//!   unstated and so their result is undefined — that refusal comes from
//!   `variant_edit_triples_reason`, not from this sort, and it is a different
//!   message.
//!
//! # One direction, deliberately
//!
//! Every normalizer here is built `ThreePrime`. That is not an omission:
//! `apply_triples` takes no `ShuffleDirection` and reads none, so there is no
//! second arm to cover — the sort key is a function of the triples alone. A
//! `FivePrime` copy of each test would re-run identical code and read as
//! coverage of a distinction that does not exist.

use crate::common::cis_apply_oracle::normalizer_for;
use ferro_hgvs::{MockProvider, Normalizer, ShuffleDirection};

/// 1-based:      1234567890
///
/// Position 9 is `G` and 10 is `G`, so an insertion at `9_10` sits flush against
/// a member acting on 10. The `A` run at 12-15 is not used here; the contig is
/// shared with the rest of the cis tests so this module cannot drift onto a
/// fixture of its own.
const SEQUENCE: &str = "GGATTACAGGCAAAAGCCTGAGGATTACAGGCATTAGCCT";

/// The window `to_sequences` returns for `spelling`, or the refusal.
fn applied(n: &Normalizer<MockProvider>, spelling: &str) -> Result<(u64, String, String), String> {
    let variant = ferro_hgvs::parse_hgvs(spelling).expect("committed row parses");
    n.to_sequences(&variant, 3)
        .map(|p| (p.position, p.reference, p.alternate))
        .map_err(|e| e.to_string())
}

/// **An insertion flush against a sibling applies, and applies the same way
/// whichever order the two are written in.**
///
/// The resulting bases are pinned exactly, not merely compared between the two
/// orderings: two orderings agreeing on a *wrong* sequence would satisfy an
/// equality-only test, and the payload's placement relative to the substituted
/// base is precisely what the sort decides.
#[test]
fn an_insertion_flush_against_a_sibling_applies_in_either_member_order() {
    let n = normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime);

    for (label, insertion_first, sibling_first, expected) in [
        (
            "insertion flush against a substitution",
            "TEMPLATE:g.[9_10insT;10G>A]",
            "TEMPLATE:g.[10G>A;9_10insT]",
            // 7-13 is `CAGGCAA`; the payload lands 5' of the substituted base.
            (7u64, "CAGGCAA", "CAGTACAA"),
        ),
        (
            "insertion flush against a deletion",
            "TEMPLATE:g.[9_10insT;10del]",
            "TEMPLATE:g.[10del;9_10insT]",
            (7, "CAGGCAA", "CAGTCAA"),
        ),
    ] {
        let first = applied(&n, insertion_first);
        let second = applied(&n, sibling_first);

        assert_eq!(
            first,
            Ok((expected.0, expected.1.to_string(), expected.2.to_string())),
            "{label}: the insertion-first spelling did not apply to the expected bases"
        );
        assert_eq!(
            first, second,
            "{label}: the two member orders denote different things, so the applier is \
             still reading order as meaning (#1098 says it carries none)"
        );
    }
}

/// **The same property over a grid, because two rows cannot carry it.**
///
/// The test above pins the two shapes the defect was reported on, exactly. That
/// is the right test for those bases and the wrong one for the *claim*: member
/// order carries no meaning (#1098) is a property of every allele, and a pair of
/// rows can only witness it at one locus. A sort key is exactly the kind of fix
/// that can be right where it was measured and wrong one base over — the old one
/// was, since `Reverse(position)` decides correctly wherever no two members
/// share a position.
///
/// So the grid walks the contig: at every interior site, an insertion at the 5'
/// flank, at the 3' flank, and one base clear, against four sibling edits — and
/// asserts the two authorings agree, whatever the verdict. A refusal is a
/// verdict too, and both orders must reach it, which is what makes this a
/// symmetry test rather than a "these now apply" test.
///
/// The two sides are compared on the **applied window**, not on the error
/// string: a refusal message may quote the description it refused, which
/// differs between the orders by construction. `Ok(w) == Ok(w)` is agreement on
/// the bases, `None == None` is agreement that there is no single resulting
/// sequence, and the two floors below stop either from being reached vacuously.
///
/// Measured on this contig: **336 pairs, 280 of which apply and 56 refuse**.
/// The floors sit well under both, so a legitimate change to what
/// `to_sequences` accepts does not force a re-bless, while a grid that stopped
/// exercising either verdict still fails.
#[test]
fn the_verdict_is_order_independent_over_a_grid_of_flush_and_clear_pairs() {
    let n = normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime);
    let bases = SEQUENCE.as_bytes();

    let (mut pairs, mut applied_pairs, mut refused_pairs) = (0usize, 0usize, 0usize);
    for site in 6..=33usize {
        let base = bases[site - 1] as char;
        let alt = if base == 'A' { 'C' } else { 'A' };
        // Flush against the site's 5' flank, flush against its 3' flank, and one
        // base clear of it — so the grid holds both sides of the shared-position
        // case and a control that never shared one.
        for insertion_point in [site - 1, site, site - 2] {
            let insertion = format!("{}_{}insT", insertion_point, insertion_point + 1);
            for sibling in [
                format!("{site}{base}>{alt}"),
                format!("{site}del"),
                format!("{}_{}del", site, site + 1),
                format!("{}_{}delinsTT", site, site + 1),
            ] {
                let insertion_first = format!("TEMPLATE:g.[{insertion};{sibling}]");
                let sibling_first = format!("TEMPLATE:g.[{sibling};{insertion}]");
                let first = applied(&n, &insertion_first);
                let second = applied(&n, &sibling_first);

                assert_eq!(
                    first.as_ref().ok(),
                    second.as_ref().ok(),
                    "`{insertion_first}` and `{sibling_first}` are the same variant written in \
                     the two orders, and they do not denote the same thing — the applier is \
                     still reading member order as meaning (#1098 says it carries none)"
                );

                pairs += 1;
                if first.is_ok() {
                    applied_pairs += 1;
                } else {
                    refused_pairs += 1;
                }
            }
        }
    }

    // A symmetry assertion is satisfied by a grid that refuses everything, and
    // one that applies everything never reaches the disjointness walk this
    // change orders the input for. Both floors are therefore part of the test.
    assert_eq!(
        pairs, 336,
        "the grid changed shape; the floors below were measured on it"
    );
    assert!(
        applied_pairs >= 200,
        "only {applied_pairs} of {pairs} grid pairs apply, so the agreement above is mostly \
         two refusals agreeing"
    );
    assert!(
        refused_pairs >= 40,
        "only {refused_pairs} of {pairs} grid pairs refuse, so the grid no longer exercises the \
         disjointness verdict this change orders the input for"
    );
}

/// **A sibling one base clear of the insertion is unaffected**, and so is a pair
/// of plain substitutions.
///
/// The negative control for the assertion above: these applied before the fix
/// and must apply identically after it, or the sort has moved something it was
/// not supposed to reach.
#[test]
fn members_that_already_applied_are_byte_identical() {
    let n = normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime);

    assert_eq!(
        applied(&n, "TEMPLATE:g.[9_10insT;11C>A]"),
        Ok((7, "CAGGCAAA".to_string(), "CAGTGAAAA".to_string())),
        "an insertion with one base of clearance moved"
    );
    assert_eq!(
        applied(&n, "TEMPLATE:g.[10G>A;11C>A]"),
        Ok((7, "CAGGCAAA".to_string(), "CAGAAAAA".to_string())),
        "two plain substitutions moved"
    );
}

/// **Members that genuinely claim the same base still refuse.**
///
/// The discriminating case. `10G>A` sits *inside* `9_11del`, so the two claim
/// base 10 between them and there is no single resulting sequence — which is
/// what `triples_are_disjoint` is for. A fix that reached this would have
/// deleted the guard rather than ordered its input.
///
/// # Why the message alone cannot carry this
///
/// `src/spdi/apply.rs` emits one string for two causes — "its members overlap,
/// **or** a stated reference base disagrees with the reference" — so
/// `contains("members overlap")` is satisfied by a refusal that was really a
/// base mismatch, and this test would stay green if a future change broke base
/// resolution instead of the overlap guard.
///
/// The positive control closes that: the same members with the substitution
/// moved *outside* the deletion's span apply cleanly, which can only happen if
/// every stated reference base here is correct. So the refusal above is the
/// overlap, by elimination rather than by trusting a disjunction.
#[test]
fn a_member_interior_to_a_deletion_still_refuses() {
    let n = normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime);
    let error = applied(&n, "TEMPLATE:g.[9_11del;10G>A]")
        .expect_err("a substitution interior to a deletion claims a base the deletion also claims");
    assert!(
        error.contains("members overlap"),
        "the refusal no longer names the overlap: {error}"
    );

    // 9-11 is `GGC`; 12 is the first base of the `A` run, and `12A>T` claims no
    // base the deletion claims.
    assert!(
        applied(&n, "TEMPLATE:g.[9_11del;12A>T]").is_ok(),
        "the control does not apply, so the refusal above cannot be attributed \
         to the overlap rather than to a stated base disagreeing"
    );
}

/// **Two pure insertions at one interbase still refuse**, and still refuse for
/// their own reason.
///
/// They tie on both sort components, so the stable sort leaves them in authored
/// order and this module's change cannot reach them. The refusal comes from
/// `variant_edit_triples_reason` — an unstated order between two payloads at one
/// point has no defined result — and the message is asserted so that a future
/// change cannot silently reroute it through the overlap path instead.
#[test]
fn coincident_pure_insertions_still_refuse_for_their_own_reason() {
    let n = normalizer_for(SEQUENCE, ShuffleDirection::ThreePrime);
    let error = applied(&n, "TEMPLATE:g.[9_10insT;9_10insA]")
        .expect_err("two payloads at one interbase have no stated order");
    assert!(
        error.contains("no single resulting sequence"),
        "the refusal is no longer the undefined-result one: {error}"
    );
    // That message is itself a disjunction of four causes, so the positive
    // content above is weak on its own. What distinguishes this path from the
    // one this change touches is that it is NOT the overlap path — asserted
    // directly, so a future edit cannot reroute it there unnoticed.
    assert!(
        !error.contains("members overlap"),
        "coincident insertions are now refused as an overlap, which is the path \
         this change orders the input for — they claim no base: {error}"
    );
}

/// **#1420's alignment rows reach the applier at all.**
///
/// The reason this defect mattered. Comment 2 of #1420 states one read against
/// one reference under four gap placements — the only alignment-varying material
/// in the reported corpus — and two of the four carry `129_130insT` flush
/// against `130G>A`, so neither could be applied, and the shape most worth
/// testing was the one that could not be reached.
///
/// Run here against the contig those rows are drawn on, which is **139 bases**
/// — `reported_confluence_pairs::TEMPLATE` is the same sequence truncated to
/// 125, and positions 129-136 do not exist in it.
#[test]
fn the_1420_alignment_rows_now_apply_and_denote_one_sequence() {
    let n = normalizer_for(
        concat!(
            "ATGCACCAGTCACCAGTCTGATGCGGATCACGTGCAATTGCACGTGCAATTGGATCCGATCG",
            "TACGTACGATCGGCATGCATGCTAGCTAGCATCGATCGTAGCTAGCTAGCATCGATCGATCGA",
            "TTTAGGCGACATTT",
        ),
        ShuffleDirection::ThreePrime,
    );

    // All four spellings of the one read, as #1420 comment 2 states them.
    let spellings = [
        "TEMPLATE:g.[130_131del;131_132insTAC;134_135insG;135_136insTT]",
        "TEMPLATE:g.[129_130insT;130G>A;131G>C;134_135insG;135_136insTT]",
        "TEMPLATE:g.[129_130insT;130G>A;131G>C;134_135insGCT;135C>T]",
        "TEMPLATE:g.[130_132del;133G>T;135_136insCGAGCTT]",
    ];

    for spelling in spellings {
        assert_eq!(
            applied(&n, spelling),
            // Derived from #1420's own read, not echoed back from ferro. The
            // issue states the alignment over 126-139 as
            //
            //     ref : TTTAGGCGACATTT
            //     read: TTTATACCGAGCTTATTT
            //
            // and `applied` asks for `pad = 3`, giving 127-138 — one base of
            // flank trimmed from each end of both strings. `TTTAGGCGACATTT`
            // minus its first and last `T` is `TTAGGCGACATT`; `TTTATACCGAGCTTATTT`
            // minus its first and last `T` is `TTATACCGAGCTTATT`. So the
            // expected value is the issue's, transformed by an arithmetic a
            // reader can check, which is what keeps this an oracle rather than
            // a change detector.
            Ok((
                127,
                "TTAGGCGACATT".to_string(),
                "TTATACCGAGCTTATT".to_string()
            )),
            "{spelling} does not apply to the read the issue draws it from"
        );
    }
}
