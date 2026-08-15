//! The one definition of what a cis-allele member **writes**, and of when two
//! members collide.
//!
//! # Why this module exists
//!
//! "Do two cis members overlap?" used to be decided in four places by three
//! different geometries, and they disagreed (#1749). The correct notion was
//! already established — key on what a member *writes*, not on what it *reads*
//! (#1411, #1437, #1448) — but it lived as prose in a comment plus two partial,
//! non-shared implementations, and a third that had never been updated to it.
//!
//! Two disagreements were internal to `overlap.rs` alone. A `repeat` was a span
//! to the coincident-bounds pass, invisible to the intersection pass, and a
//! junction-*only* occupant to the junction pass — three ladders, three answers,
//! none reachable from a comment on the others.
//!
//! So the geometry is defined **once**, here, and every consumer derives its
//! question from it. Adding a fourth ladder is now a visible act rather than an
//! oversight.
//!
//! # The model
//!
//! A member writes in at most two places, and the two are independent:
//!
//! - over a **closed range of reference bases**, which it may or may not leave
//!   standing; and
//! - at a **zero-width junction**, the slot immediately 3' of some position.
//!
//! | edit | bases | junction | bases survive |
//! |---|---|---|---|
//! | `sub` / `delins` / `inv` | its span | — | yes |
//! | `del` | its span | — | **no** |
//! | `dup` (certain extent) | — | 3' of its span (`duplication.md:5`) | n/a |
//! | `dup` (uncertain extent) | its span | — | yes |
//! | `ins` | — | between its two flanks | n/a |
//! | `repeat` (grows) | — | 3' of its span | n/a |
//! | `repeat` (shrinks) | the trailing bases it removes | — | **no** |
//!
//! Two entries in that table are the ones that were being answered three ways,
//! and each is a deliberate choice rather than a simplification:
//!
//! **`repeat` is decided by arithmetic, not by a provider.** A repeat states
//! both the tract it replaces and what that tract becomes, so `unit x count`
//! against the tract length says whether it grows (the extra copies land at the
//! junction 3' of the tract — a `dup`'s geometry) or shrinks (the trailing bases
//! are removed and the leading ones are untouched reference). The three ladders
//! this module replaces each guessed instead: one kept the whole span, one
//! dropped the member, one made it junction-only. The construction is on
//! `super::overlap::repeat_footprint`, along with the shapes it declines.
//!
//! **An uncertain-extent `dup` keeps its span.** A certain `dup`'s write
//! footprint is reference-*independent* — the copy goes directly 3' of the
//! original, whatever the bases are. An uncertain one's is not known to be that
//! junction at all, so it is not entitled to the narrow footprint.
//!
//! # What this module is not
//!
//! It is not the spec's `general.md:58` prohibition on a `del` combined with a
//! `dup` of part of the same sequence. That clause names an **edit-type pair**
//! and keys on *read* spans by its own worked example
//! (`NM_004006.2:c.[762_768del;767_774dup]`, which overlaps at 767–768 and does
//! not cancel), so it is not a statement of this geometry and cannot be derived
//! from it. It stays where it is, in `HgvsVariant::detect_self_cancelling_pair`,
//! and is documented there as a spec-explicit exception rather than as a fourth
//! opinion on overlap. See #1749 and
//! `rulings[conflicting-member-geometry-refusal-scope]`.

/// Where a member writes, in whatever position space the caller is working in.
///
/// `P` is the position type: [`super::overlap`] instantiates it at ranked
/// HGVS axis positions, [`crate::spdi::apply`] at 0-based interbase
/// coordinates. The geometry is identical in both; only the mapping onto the
/// axis differs, which is exactly the part that should not be shared.
///
/// Ranges are **closed** — `bases: Some((4, 4))` is one base. An interbase
/// caller converts its half-open `[position, position + len)` on the way in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct WriteFootprint<P> {
    /// The closed range of reference bases this member rewrites, or `None` when
    /// it rewrites none.
    pub(crate) bases: Option<(P, P)>,
    /// The zero-width junction immediately 3' of this position, when the member
    /// writes there.
    pub(crate) junction: Option<P>,
    /// Whether the bases named by `bases` are still there afterwards.
    ///
    /// Only a pure deletion sets this `false`, and it is what makes
    /// `g.[2_3del;2_3insAA]` compose uniquely: the deletion removes every base
    /// it spans, so an interior junction has nothing left to be positioned
    /// against, whichever interior junction the insertion named (#1406).
    /// Meaningless when `bases` is `None`.
    pub(crate) preserves_bases: bool,
}

impl<P: Ord + Copy> WriteFootprint<P> {
    /// A member that rewrites the closed range `start..=end` and nothing else.
    pub(crate) fn spanning(start: P, end: P, preserves_bases: bool) -> Self {
        Self {
            bases: Some((start, end)),
            junction: None,
            preserves_bases,
        }
    }

    /// A member that writes only at the junction immediately 3' of `after`.
    pub(crate) fn at_junction(after: P) -> Self {
        Self {
            bases: None,
            junction: Some(after),
            preserves_bases: true,
        }
    }

    /// The closed base range, when it is one this module is willing to read.
    ///
    /// A reversed range is declined rather than judged. SVD-WG006 admits
    /// `<high>_<low>` for a circular deletion or duplication, and nothing
    /// upstream reorders those, so read as an ordinary interval `m.16569_5del`
    /// claims to end before it begins and every test below would decide it by
    /// accident. A missed report is a verdict this module simply does not make;
    /// a wrong one is a rejection the caller cannot explain (#1423).
    fn readable_bases(&self) -> Option<(P, P)> {
        self.bases.filter(|(start, end)| start <= end)
    }

    /// Whether this member's junction falls **strictly inside** `other`'s
    /// surviving bases.
    ///
    /// Strict at the 3' end (`gap < end`), so a junction flush against either
    /// edge of a span is not interior: `g.[10_14del;14_15insT]` and
    /// `g.[10_14del;9_10insT]` each compose uniquely. It is also what stops a
    /// `repeat` conflicting with itself — its own junction is its own `end`.
    fn junction_interior_to(&self, other: &Self) -> bool {
        let (Some(gap), Some((start, end))) = (self.junction, other.readable_bases()) else {
            return false;
        };
        other.preserves_bases && start <= gap && gap < end
    }

    /// Whether two members write to a place that cannot hold both.
    ///
    /// **Never call this with a member against itself** — a member's bases
    /// trivially intersect their own. Callers pair distinct indices.
    pub(crate) fn conflicts_with(&self, other: &Self) -> bool {
        // (1) Two members rewriting an intersecting range of reference bases.
        //     Closed intervals intersect when neither ends before the other
        //     begins.
        if let (Some((a_start, a_end)), Some((b_start, b_end))) =
            (self.readable_bases(), other.readable_bases())
        {
            if a_start <= b_end && b_start <= a_end {
                return true;
            }
        }

        // (2) One junction, two writers, and nothing in the description saying
        //     which goes first. HGVS spells "insert both, in this order" as a
        //     single ordered compound payload (`ins[A;C]`, general.md:79), so an
        //     author who meant an order had a way to say it.
        if let (Some(a), Some(b)) = (self.junction, other.junction) {
            if a == b {
                return true;
            }
        }

        // (3) A junction strictly interior to a span whose bases survive, where
        //     the junction's position within the span is meaningful and the
        //     combination therefore has no single answer.
        self.junction_interior_to(other) || other.junction_interior_to(self)
    }

    /// [`Self::conflicts_with`], but treating every span as if its bases
    /// survive.
    ///
    /// This is the **splice algorithm's precondition**, which is a different
    /// question from whether the allele denotes one sequence, and conflating
    /// the two is what [`crate::spdi::apply`] used to do.
    ///
    /// `apply_triples` walks 3' -> 5' so that each stated deletion can be
    /// validated against the pristine reference. An insertion sitting inside a
    /// deletion breaks that walk — the deletion's splice would consume bases the
    /// insertion had already put into the buffer — even though the pair denotes
    /// one sequence perfectly well (#1406). So the applier must decline a shape
    /// it cannot *execute*, without that decline being read as a verdict on the
    /// description.
    ///
    /// **Rule (2) is deliberately omitted.** Two zero-width writes at one
    /// junction are not a splice hazard: both splice at the same offset, and
    /// the buffer survives. The ambiguity they *do* carry — that nothing in the
    /// description says which goes first — is caught one level up, in
    /// `spdi::apply::variant_edit_triples_reason`, and it has to stay there.
    /// This predicate is shared with `EquivalenceChecker`, where the permissive
    /// reading is correct because a decline only forgoes *upgrading* a
    /// `NotEquivalent` verdict and so can never invent an equivalence, and two
    /// pinned tests depend on that. One predicate cannot serve both a checker
    /// that may guess and a published key that may not.
    pub(crate) fn obstructs_splice_of(&self, other: &Self) -> bool {
        if let (Some((a_start, a_end)), Some((b_start, b_end))) =
            (self.readable_bases(), other.readable_bases())
        {
            if a_start <= b_end && b_start <= a_end {
                return true;
            }
        }
        let survive = |f: &Self| Self {
            preserves_bases: true,
            ..*f
        };
        survive(self).junction_interior_to(&survive(other))
            || survive(other).junction_interior_to(&survive(self))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span(start: u64, end: u64) -> WriteFootprint<u64> {
        WriteFootprint::spanning(start, end, true)
    }

    fn deletion(start: u64, end: u64) -> WriteFootprint<u64> {
        WriteFootprint::spanning(start, end, false)
    }

    #[test]
    fn intersecting_base_ranges_conflict() {
        assert!(span(10, 14).conflicts_with(&span(12, 16)));
        assert!(span(10, 14).conflicts_with(&span(10, 16)));
        assert!(span(10, 14).conflicts_with(&span(10, 14)));
        assert!(!span(10, 14).conflicts_with(&span(15, 19)));
    }

    #[test]
    fn a_junction_is_disjoint_from_the_span_it_sits_beside() {
        // Flush at the 3' edge and at the 5' edge.
        assert!(!WriteFootprint::at_junction(14).conflicts_with(&span(10, 14)));
        assert!(!WriteFootprint::at_junction(9).conflicts_with(&span(10, 14)));
    }

    #[test]
    fn an_interior_junction_conflicts_only_where_the_bases_survive() {
        assert!(WriteFootprint::at_junction(14).conflicts_with(&span(10, 20)));
        assert!(!WriteFootprint::at_junction(14).conflicts_with(&deletion(10, 20)));
    }

    #[test]
    fn the_splice_precondition_ignores_whether_bases_survive() {
        // The #1406 pair denotes one sequence but cannot be spliced in one walk.
        assert!(!WriteFootprint::at_junction(14).conflicts_with(&deletion(10, 20)));
        assert!(WriteFootprint::at_junction(14).obstructs_splice_of(&deletion(10, 20)));
        // A flush junction obstructs nothing either way.
        assert!(!WriteFootprint::at_junction(9).obstructs_splice_of(&deletion(10, 20)));
    }

    #[test]
    fn two_writers_at_one_junction_conflict() {
        assert!(WriteFootprint::at_junction(14).conflicts_with(&WriteFootprint::at_junction(14)));
        assert!(!WriteFootprint::at_junction(14).conflicts_with(&WriteFootprint::at_junction(15)));
    }

    /// A `<high>_<low>` range is declined by [`WriteFootprint::readable_bases`]
    /// rather than read as an interval.
    ///
    /// # The comparand has to SPAN both stated endpoints
    ///
    /// The first two cases below cannot see the decline at all, and this test
    /// shipped with only those: against `span(1, 10)` the intersection arm is
    /// `16569 <= 10 && 1 <= 5`, whose *first* conjunct is already false, so it
    /// answers "no conflict" whether or not the reversed range was declined.
    /// Deleting the `start <= end` filter passed the entire suite (10,487
    /// tests) including this test — the guard had no test that could fail.
    ///
    /// `span(1, 20000)` is the discriminating shape: read as an ordinary
    /// interval, `(16569, 5)` intersects it at both ends (`16569 <= 20000` and
    /// `1 <= 5`), so an undeclined reversed range reports a conflict here. Any
    /// comparand asserted against must contain both `16569` and `5`, which is
    /// exactly the geometry SVD-WG006's circular `m.16569_5del` has against a
    /// whole-molecule sibling.
    ///
    /// Asserted on both predicates, since both read their spans through
    /// `readable_bases`.
    #[test]
    fn a_reversed_range_is_declined_rather_than_judged() {
        let reversed = WriteFootprint::spanning(16569, 5, true);
        assert!(!reversed.conflicts_with(&span(1, 10)));
        assert!(!WriteFootprint::at_junction(7).conflicts_with(&reversed));
        assert!(
            !reversed.conflicts_with(&span(1, 20000)),
            "a reversed range must be declined, not read as an interval that \
             happens to intersect a comparand spanning both its endpoints"
        );
        assert!(
            !reversed.obstructs_splice_of(&span(1, 20000)),
            "the splice precondition reads its spans through the same \
             `readable_bases`, so it declines a reversed range too"
        );
    }
}
