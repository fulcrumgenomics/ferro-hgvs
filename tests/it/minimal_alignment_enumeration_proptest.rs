//! Property tests for `common::minimal_alignment`.
//!
//! `minimal_alignment_enumeration` pins the record's own worked cases, which
//! checks the instrument against the six blocks whose answers are already
//! written down. These properties check the other direction: that the
//! enumerator is internally coherent on inputs nobody chose, so a pinned case
//! cannot be right by accident.
//!
//! # What is asserted, and what deliberately is not
//!
//! | property | why it is checkable |
//! |---|---|
//! | every emitted alignment is well formed and costs the reported minimum | recomputes the cost from the columns, independently of the DP grid |
//! | the emitted alignments are pairwise distinct | a duplicate would inflate the count without changing the intersection, which is the one error the count cannot reveal |
//! | `unchanged` is exactly the intersection of the per-alignment matched sets | the record's definition, restated against the emitted set |
//! | the enumerator agrees with #1539's closed-form detector | a second implementation of the same notion, by a route that never enumerates |
//! | the cap refuses iff the set does not fit, and otherwise changes nothing | a truncated intersection is too *permissive*, so this is the safety property |
//! | reversing both blocks mirrors the answer | catches a DP or backtrack asymmetry, which a symmetric corpus would hide |
//! | `cost(sub = 1) <= cost(sub = 2) <= 2 * cost(sub = 1)` | the two models are related even where their unchanged sets are not |
//! | on equal-length blocks the minimal cost never exceeds the position-wise cost | the record's premise, stated as an inequality |
//!
//! **The brief for this module proposed one further property — that `unchanged`
//! is always a subset of the position-wise matched set — and it is FALSE.**
//! Brute force over every `{A, C}` pair up to length 5 under substitution cost 1
//! finds 178 counterexamples to that containment and 250 to its converse; the
//! record's own `CAG -> AGA` is one of the first, with `unchanged = {1, 2}` and
//! no position-wise match at all. The two notions are incomparable, which is a
//! stronger statement of the record's headline claim than "they can differ", and
//! it is pinned with witnesses in
//! `minimal_alignment_enumeration::neither_notion_contains_the_other`. A
//! proptest cannot assert that a difference *exists*, so it is not attempted
//! here — asserting the false containment would have gone red on roughly one
//! generated pair in four.

use crate::common::minimal_alignment::{
    minimal_alignments, minimal_alignments_capped, position_wise_matches, CostModel, Step,
    DEFAULT_ALIGNMENT_CAP,
};
use proptest::prelude::*;

/// The longest block either side may be.
///
/// Six, so that the worst case stays comfortably inside
/// [`DEFAULT_ALIGNMENT_CAP`] and the properties can assert enumeration
/// *succeeds* rather than tolerating a refusal — a tolerated refusal is how a
/// property test goes vacuous. The worst case at length `n` is an all-mismatch
/// block under substitution cost 2, which has `D(n, n)` optimal alignments (the
/// central Delannoy numbers): 8989 at `n = 6`, against a cap of 65536. Chosen
/// for runtime rather than for the cap, which would admit `n = 7`.
const MAX_BLOCK_LEN: usize = 6;

/// Bases to draw from.
///
/// Two alphabets rather than one. `ACGT` is the realistic draw, but at length
/// six a four-letter alphabet rarely produces a block whose minimal alignment
/// is interesting — most pairs share nothing and the answer is trivial. The
/// two-letter draw concentrates the generator on repeat-like blocks, which is
/// where several minimal alignments compete and the intersection actually does
/// work.
fn block_alphabet() -> impl Strategy<Value = Vec<u8>> {
    prop_oneof![Just(b"AC".to_vec()), Just(b"ACGT".to_vec()),]
}

/// A `(reference, alternate)` pair with independently drawn lengths.
fn block_pair() -> impl Strategy<Value = (Vec<u8>, Vec<u8>)> {
    block_alphabet().prop_flat_map(|alphabet| {
        let base = prop::sample::select(alphabet);
        (
            prop::collection::vec(base.clone(), 0..=MAX_BLOCK_LEN),
            prop::collection::vec(base, 0..=MAX_BLOCK_LEN),
        )
    })
}

/// A `(reference, alternate)` pair of equal length — the shape the record's
/// headline claim is about, since equal length is what makes the position-wise
/// reading look available.
fn equal_length_pair() -> impl Strategy<Value = (Vec<u8>, Vec<u8>)> {
    block_alphabet().prop_flat_map(|alphabet| {
        (0..=MAX_BLOCK_LEN).prop_flat_map(move |len| {
            let base = prop::sample::select(alphabet.clone());
            (
                prop::collection::vec(base.clone(), len..=len),
                prop::collection::vec(base, len..=len),
            )
        })
    })
}

/// Every model, so each property sweeps both rather than picking one.
fn cost_model() -> impl Strategy<Value = CostModel> {
    prop::sample::select(CostModel::ALL.to_vec())
}

/// Reference offsets an alignment matches, recomputed here rather than read
/// back from the report, so the properties do not check the instrument against
/// itself.
fn matched_reference_offsets(alignment: &[Step]) -> Vec<u32> {
    let mut matched = Vec::new();
    let mut offset = 0u32;
    for &step in alignment {
        match step {
            Step::Match => {
                matched.push(offset);
                offset += 1;
            }
            Step::Sub | Step::Del => offset += 1,
            Step::Ins => {}
        }
    }
    matched
}

proptest! {
    /// **Every emitted alignment is a well-formed alignment of exactly these two
    /// blocks, and costs the reported minimum.**
    ///
    /// The cost is recomputed from the columns, so this checks the enumerated
    /// paths against the dynamic-programming grid rather than against
    /// themselves. `Match` and `Sub` are checked to agree with the bases they
    /// pair, which is what makes the matched-offset bookkeeping downstream
    /// meaningful.
    #[test]
    fn every_emitted_alignment_is_well_formed_and_minimal(
        (reference, alternate) in block_pair(),
        model in cost_model(),
    ) {
        let report = minimal_alignments(&reference, &alternate, model)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;

        for alignment in report.alignments() {
            let mut i = 0usize;
            let mut j = 0usize;
            let mut cost = 0u32;
            for &step in alignment {
                cost += model.step_cost(step);
                match step {
                    Step::Match => {
                        prop_assert_eq!(
                            reference[i], alternate[j],
                            "a Match column pairs unequal bases at ({}, {})", i, j
                        );
                        i += 1;
                        j += 1;
                    }
                    Step::Sub => {
                        prop_assert_ne!(
                            reference[i], alternate[j],
                            "a Sub column pairs equal bases at ({}, {})", i, j
                        );
                        i += 1;
                        j += 1;
                    }
                    Step::Del => i += 1,
                    Step::Ins => j += 1,
                }
            }
            prop_assert_eq!(i, reference.len(), "the alignment did not consume the reference");
            prop_assert_eq!(j, alternate.len(), "the alignment did not consume the alternate");
            prop_assert_eq!(cost, report.cost(), "an emitted alignment is not minimal");
        }
    }

    /// **No alignment is emitted twice.**
    ///
    /// A duplicate would inflate `count()` while leaving `unchanged()` correct,
    /// so it is invisible to every other property here — and `count()` is what
    /// separates "one alignment pins this" from "several agree", which is the
    /// distinction the record turns on.
    #[test]
    fn the_emitted_alignments_are_pairwise_distinct(
        (reference, alternate) in block_pair(),
        model in cost_model(),
    ) {
        let report = minimal_alignments(&reference, &alternate, model)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        let mut seen: Vec<&Vec<Step>> = report.alignments().iter().collect();
        seen.sort_unstable();
        let before = seen.len();
        seen.dedup();
        prop_assert_eq!(seen.len(), before, "an alignment was enumerated more than once");
        prop_assert_eq!(report.count(), before, "count() disagrees with alignments()");
    }

    /// **`unchanged` is exactly the intersection of the per-alignment matched
    /// reference offsets** — the record's definition, restated.
    #[test]
    fn unchanged_is_the_intersection_over_every_minimal_alignment(
        (reference, alternate) in block_pair(),
        model in cost_model(),
    ) {
        let report = minimal_alignments(&reference, &alternate, model)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;

        let mut intersection: Option<Vec<u32>> = None;
        for alignment in report.alignments() {
            let matched = matched_reference_offsets(alignment);
            intersection = Some(match intersection {
                None => matched,
                Some(so_far) => so_far.into_iter().filter(|o| matched.contains(o)).collect(),
            });
        }
        let expected = intersection.unwrap_or_default();
        prop_assert_eq!(report.unchanged(), expected.as_slice());

        // And every unchanged offset really is a reference offset.
        for &offset in report.unchanged() {
            prop_assert!((offset as usize) < reference.len());
        }
    }

    /// **The cap refuses exactly when the set does not fit, and otherwise
    /// changes nothing.**
    ///
    /// This is the safety property. Because `unchanged` is an intersection,
    /// enumerating fewer alignments can only *add* offsets to it, so a
    /// silently truncated answer would report bases as unchanged that are not.
    /// A cap one below the true count must therefore refuse, and a cap at the
    /// true count must give bit-identical output to the uncapped call.
    #[test]
    fn the_cap_refuses_rather_than_truncating(
        (reference, alternate) in block_pair(),
        model in cost_model(),
    ) {
        let report = minimal_alignments(&reference, &alternate, model)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        let count = report.count();

        let exact = minimal_alignments_capped(&reference, &alternate, model, count)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        prop_assert_eq!(&exact, &report, "a cap at the exact count changed the answer");

        let tight = minimal_alignments_capped(&reference, &alternate, model, count - 1);
        let exceeded = tight.expect_err("a cap below the true count must refuse");
        prop_assert_eq!(exceeded.cap, count - 1);
        prop_assert_eq!(exceeded.cost, report.cost(), "the cost survives a refusal");
        prop_assert_eq!(exceeded.model, model);
    }

    /// **Reversing both blocks mirrors the answer.**
    ///
    /// Edit distance is invariant under reversal, so the alignment count must be
    /// too and each unchanged offset `i` must become `len - 1 - i`. A backtrack
    /// that preferred one end — a 5'/3' bias smuggled into the enumerator —
    /// would break this while leaving every pinned case, all of which are
    /// checked in one direction only, entirely green.
    #[test]
    fn reversing_both_blocks_mirrors_the_answer(
        (reference, alternate) in block_pair(),
        model in cost_model(),
    ) {
        let forward = minimal_alignments(&reference, &alternate, model)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;

        let reversed_reference: Vec<u8> = reference.iter().rev().copied().collect();
        let reversed_alternate: Vec<u8> = alternate.iter().rev().copied().collect();
        let backward = minimal_alignments(&reversed_reference, &reversed_alternate, model)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;

        prop_assert_eq!(forward.cost(), backward.cost());
        prop_assert_eq!(forward.count(), backward.count());

        let last = reference.len() as u32;
        let mut mirrored: Vec<u32> = backward
            .unchanged()
            .iter()
            .map(|&offset| last - 1 - offset)
            .collect();
        mirrored.sort_unstable();
        prop_assert_eq!(forward.unchanged(), mirrored.as_slice());
    }

    /// **`cost(sub = 1) <= cost(sub = 2) <= 2 * cost(sub = 1)`.**
    ///
    /// The two models are incomparable on their unchanged sets — which is the
    /// gap `the_two_cost_models_disagree_on_1420_v4` records — but not
    /// unrelated: every column costs at least as much under the second model,
    /// and any alignment can trade each substitution for a deletion and an
    /// insertion at exactly double the price.
    #[test]
    fn the_two_cost_models_bracket_each_other(
        (reference, alternate) in block_pair(),
    ) {
        let unit = minimal_alignments(&reference, &alternate, CostModel::Levenshtein)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        let paired = minimal_alignments(&reference, &alternate, CostModel::SubstitutionCostsTwo)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        prop_assert!(unit.cost() <= paired.cost());
        prop_assert!(paired.cost() <= 2 * unit.cost());
    }

    /// **The enumerator and the closed form agree on every block.**
    ///
    /// `issue_1539_split_member_separation::forced_unchanged_columns` computes
    /// the same notion in `Θ(n·m)` and without materialising an alignment: it
    /// marks a column changeable when *some* minimum-cost path deletes or
    /// mismatches it. That detector is the one
    /// `rulings[unchanged-is-read-over-every-minimal-alignment]` names as where
    /// the ruling came from, so agreement with it is agreement with the record's
    /// working definition — reached by a genuinely different route, since one
    /// side never enumerates and the other never inspects the grid.
    ///
    /// It is hardcoded to substitution cost 1, so only
    /// [`CostModel::Levenshtein`] is comparable. That is itself the strongest
    /// evidence available about which model the record intends, and it is still
    /// not a ruling — see
    /// `minimal_alignment_enumeration::the_two_cost_models_disagree_on_1420_v4`.
    #[test]
    fn the_closed_form_and_the_enumerator_agree(
        (reference, alternate) in block_pair(),
    ) {
        let report = minimal_alignments(&reference, &alternate, CostModel::Levenshtein)
            .map_err(|e| TestCaseError::fail(e.to_string()))?;
        let closed_form = crate::issue_1539_split_member_separation::forced_unchanged_columns(
            &reference,
            &alternate,
        )
        .expect("these blocks are far below the closed form's own matrix cap");

        let from_closed_form: Vec<u32> = closed_form
            .iter()
            .enumerate()
            .filter(|(_, forced)| **forced)
            .map(|(i, _)| i as u32)
            .collect();
        prop_assert_eq!(
            report.unchanged(),
            from_closed_form.as_slice(),
            "{:?} -> {:?}: the enumerator and the closed form disagree",
            String::from_utf8_lossy(&reference),
            String::from_utf8_lossy(&alternate),
        );
    }

    /// **On an equal-length block the minimal cost never exceeds the
    /// position-wise cost** — the record's premise, as an inequality.
    ///
    /// The position-wise correspondence is always *an* alignment of two blocks
    /// of the same length, so it can never beat the minimum. When it is strictly
    /// beaten, the block is one where the record's rule and the position-wise
    /// reading may disagree — `CAG -> AGA` being the smallest instance. The
    /// inequality holds either way, which is why this is the assertable half.
    #[test]
    fn on_equal_length_blocks_the_minimum_never_exceeds_the_position_wise_cost(
        (reference, alternate) in equal_length_pair(),
    ) {
        let position_wise_cost = (reference.len() - position_wise_matches(&reference, &alternate).len()) as u32;
        for model in CostModel::ALL {
            let report = minimal_alignments(&reference, &alternate, model)
                .map_err(|e| TestCaseError::fail(e.to_string()))?;
            prop_assert!(
                report.cost() <= position_wise_cost * model.substitution_cost(),
                "{:?} -> {:?} under {model:?}: minimal cost {} exceeds the position-wise cost {}",
                String::from_utf8_lossy(&reference),
                String::from_utf8_lossy(&alternate),
                report.cost(),
                position_wise_cost * model.substitution_cost(),
            );
        }
    }
}

/// The `MAX_BLOCK_LEN` bound above is only sound if the worst case at that
/// length fits the default cap. Checked rather than asserted in prose, since a
/// later change to either constant would otherwise make every property above
/// tolerate a refusal it currently rules out.
#[test]
fn the_generator_bound_stays_inside_the_default_cap() {
    let worst_reference = vec![b'A'; MAX_BLOCK_LEN];
    let worst_alternate = vec![b'C'; MAX_BLOCK_LEN];
    let report = minimal_alignments(
        &worst_reference,
        &worst_alternate,
        CostModel::SubstitutionCostsTwo,
    )
    .expect("the worst case at MAX_BLOCK_LEN must fit the default cap");
    assert_eq!(report.count(), 8989, "the central Delannoy number D(6, 6)");
    assert!(report.count() < DEFAULT_ALIGNMENT_CAP);
}
