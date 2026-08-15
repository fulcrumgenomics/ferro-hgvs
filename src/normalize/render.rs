//! The single stage that chooses how a settled variant is **spelled** (#1946).
//!
//! Normalization interleaves two kinds of rewrite. **Semantic** rewrites — the
//! 3'/5' shuffle, the cis collapse, the sequence-first partition — change *which*
//! variant is described, or how it is cut into members; every downstream stage has
//! to be able to evaluate their output. **Presentational** rewrites — choosing
//! `dup` over `ins`, a copy count over `dup`, `ins<range>inv` over a literal
//! payload — change only how the same change is written.
//!
//! Today the presentational choices are made **per member**, inside
//! `normalize_na_edit`, before the stages that must compute on them, and the
//! codebase pays to undo them afterwards: `lower_repeat_edits`,
//! `demote_repeats_spanning_siblings`, `demote_coincident_tract_repeats` and
//! `respell_colliding_duplications` exist for that and nothing else. Measured by
//! ablation over the 86,170-row synthetic corpus, those four move **1,279**,
//! **1,011**, **0** and **505** rows respectively, for a union of **2,761 (3.2%)**
//! — the three live ones are near-disjoint, overlapping on 34 rows in one pair.
//! (`demote_coincident_tract_repeats`'s zero is structural, not safety: it builds a
//! candidate group 7,413 times and every one is of size 1.)
//!
//! This module is where those choices are meant to move, once the member set is
//! final. **It does nothing yet, deliberately** — see
//! [`render_canonical_spelling`].

use crate::hgvs::variant::HgvsVariant;
use crate::reference::ReferenceProvider;

/// Apply the `general.md:56` prioritisation ladder to a settled member list.
///
/// **A strict no-op today.** The seam is landed on its own so that the placement
/// and the argument for it are reviewable separately from any behaviour change;
/// when a form does start moving, the diff shows only the form. Every existing row
/// of the synthetic corpus is byte-identical either side of this commit, which is
/// the same proof an added corpus family carries.
///
/// # The placement is pinned, not chosen
///
/// This runs between `detect_overlap_conflicts` and
/// `sort_cis_members_by_genomic_order` in `normalize_allele`, and both bounds are
/// forced rather than convenient:
///
/// - **After `detect_overlap_conflicts`.** A tandem-prefix split emits two pieces
///   sharing a junction. Run before the detector, that reads as an overlap conflict
///   and the allele is preserved verbatim and rejected by strict mode — the stage
///   would manufacture the very fault the detector exists to report.
/// - **Before `sort_cis_members_by_genomic_order`.** A stage that may emit new
///   pieces must emit them before the pass that puts members in coordinate order,
///   or the output is unordered and two spellings of one allele stop converging
///   (#1098).
/// - **After `emitted_verbatim` is computed**, which is the finer point inside that
///   window. That flag asks whether the allele came through the *derivation*
///   unchanged, and it gates the #1406 rule that a real conflict must survive into
///   the output. A presentational re-spelling must not be able to flip a claim
///   about the derivation, so the question is asked before this stage answers a
///   different one.
///
/// # The invariant is apply-equality, and it is not member-count stability
///
/// This stage **may change the number of members.** A tandem-prefix split turns one
/// long literal insertion into `[dup;ins]`, and those are different member sets
/// denoting identical bases.
///
/// So the property it preserves is denotational: applying the member list before
/// and after must produce the same sequence. That is exactly the relation the
/// release gate is stated over — the `confluence-gate-is-apply-equality-on-every-determined-axis`
/// ruling restates "equivalent inputs produce one canonical output" as *normalize is
/// constant on each equivalence class, where equivalence is apply-equality on every
/// determined axis* (`EquivalenceLevel::CrossAxisSequenceMatch`), never
/// `NormalizedMatch`. Because that relation is over `apply`, whose codomain is
/// bases, it survives a member-count change by construction. Anything keying on
/// member identity or member count instead will not, and that is the thing to check
/// when this stage gains behaviour.
///
/// # What the spec does and does not require here
///
/// Stated precisely, because the obvious reading is wrong and would otherwise
/// propagate. `DNA/duplication.md:18` — "when a variant can be described as a
/// duplication, it **must** be described as a duplication" — does **not** require a
/// tandem-prefix split. The `duplication-must-ranks-the-label-not-the-partition`
/// ruling holds that `:18` requires a change which *is* a tandem duplication to be
/// **labelled** `dup`, and does not require that the edit set be partitioned so as
/// to produce one; a split undertaken in order to expose a `dup` is precisely that
/// partitioning.
///
/// What the split does have is `:17`, which *permits* the label once the piece
/// exists ("directly 3'-flanking"), and `:90`, which publishes the geometry
/// (`c.[675-542_1211-703dup;1211-703_1211-702insGTAAA]`, with a note that it is
/// deliberately not spelled `dupins`). `delins-merge-vs-individual-gap-two-or-more`
/// cites that same passage for the same reading — it publishes a net insertion as a
/// `[dup;ins]` split. Legal and idiomatic, therefore, but **ours to choose**.
///
/// So the justification is a **house choice** in this repository's own sense, and it
/// owes a `rulings` record carrying a `house_choice` object — `decided`, naming no
/// governing and no deviated-from clause, saying what was considered and rejected,
/// and never citable as conformance. **No such record is on `origin/main` today**
/// (checked against the ledger, and against all 78 remote branches), so nothing here
/// cites one by id: a doc comment pointing at a ruling that does not exist is worse
/// than one that says the decision is still owed. Write the record before the
/// tandem-prefix consumer ships, not after.
///
/// # The coincidence floor is house policy too, and its constant must not be a corpus size
///
/// Both consumers need a floor below which a match is chance rather than structure —
/// a short tandem prefix matches by accident constantly. Two cautions, because the
/// obvious shortcuts are both wrong.
///
/// `changed_columns_dominate_the_span` is **not** that floor and must not be
/// borrowed as one. Its own doc says it "is **not** evidence about whether a given
/// reverse-complement relation is structural", and its analysis says why: relative
/// spread is `sqrt(3/m)`, so a density reading carries almost no information about a
/// short block — which is exactly the regime a tandem *prefix* lives in.
///
/// And whatever constant the floor ends up with **must read the corpus size or be
/// stated as an order of magnitude, never hardcode a row count.** The synthetic
/// corpus has moved twice in this campaign alone — 85,642 to 85,930 to 86,170 — so a
/// literal pinned to one of those is stale by the next family. That is the failure
/// `the_corpus_emits_a_block_past_the_split_cap` already shipped once, where a guard
/// restating `MAX_SPLIT_BLOCK` as `1024` stayed green after the constant moved to
/// 4096.
///
/// # Why a provider is taken now
///
/// Both waiting consumers need reference context — an `ins<range>inv` payload has to
/// be checked against the bases the range names, and a tandem prefix against the
/// reference immediately 5' of the insertion point — so the borrow is established
/// with the seam rather than threaded through later.
pub(crate) fn render_canonical_spelling<P: ReferenceProvider>(
    members: &mut Vec<HgvsVariant>,
    provider: &P,
) {
    // Deliberately empty. Named rather than absent so the call site, the ordering
    // argument above and the invariant are reviewable on their own.
    let _ = (members, provider);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_hgvs;
    use crate::reference::MockProvider;
    use crate::spdi::apply::{compare_denoted_sequences, DenotedSequenceComparison};

    /// A 40-base contig with a tandem `ACG` tract, enough for the shapes below.
    fn provider() -> MockProvider {
        let mut p = MockProvider::new();
        p.add_genomic_sequence("NC_REND.1", "GTCAACGACGACGTGCATTGCATTGCATTGCATTGCATTG");
        p
    }

    fn members(descriptions: &[&str]) -> Vec<HgvsVariant> {
        descriptions
            .iter()
            .map(|d| parse_hgvs(d).expect("test descriptions parse"))
            .collect()
    }

    /// The stage preserves the bases its member list denotes.
    ///
    /// **This is the invariant, not a no-op check, and the distinction is the whole
    /// reason it is written this way.** A test asserting "this function changes
    /// nothing" is a change detector for the current implementation: it must be
    /// *deleted* the moment the stage does its job, so it guards nothing past today.
    /// Apply-equality is the property the stage must hold forever — including once it
    /// splits one member into two — so this test keeps working rather than being
    /// thrown away.
    ///
    /// # `Agree`, never merely "not a fault"
    ///
    /// `compare_denoted_sequences` declines by default, returning `NotComparable` for
    /// a whole family of reasons, so asserting the *absence* of `Differ` is satisfied
    /// by a comparison that never happened — a skip reading as a pass. That is not a
    /// hypothetical: an earlier revision of this test asserted the weaker form, and a
    /// stage body mutated to **drop a member** passed it. The sabotage is what found
    /// it, which is why the sabotage is worth running rather than assuming.
    ///
    /// # What this harness can and cannot ask, stated rather than hidden
    ///
    /// Every case below is a **single member**, and that is a limit of the instrument
    /// rather than a choice about coverage. Measured in this harness, a hand-built
    /// two-member `AlleleVariant` returns `NotComparable(InputDenotesNoSequence)` —
    /// on `[g.14_15insACG;g.7C>A]`, `[g.6_14dup;g.20G>A]` and `[g.7_9del;g.20G>A]`
    /// alike — while the single-member forms of the same edits compare and agree. A
    /// `Repeat` member declines too, for a documented reason: SPDI cannot carry that
    /// edit.
    ///
    /// So the multi-member half of the invariant is **not** asserted here. It is
    /// covered where the comparison does run over real alleles: the corpus harness's
    /// `--verify-spdi` pass, which compared 86,170 rows and reported the seam moving
    /// none. Saying that plainly is the point — listing a case the oracle cannot
    /// answer would put the silence straight back.
    #[test]
    fn the_stage_preserves_the_bases_its_members_denote() {
        let provider = provider();
        let cases = [
            "NC_REND.1:g.18_19insTTTTTT",
            "NC_REND.1:g.6_14dup",
            "NC_REND.1:g.7_9del",
            "NC_REND.1:g.20G>A",
        ];
        for case in cases {
            let before = members(&[case]);
            let mut after = before.clone();
            render_canonical_spelling(&mut after, &provider);

            // A member-count change has no member-wise counterpart, so the
            // comparison is made over the whole list rather than pairwise — and a
            // list that no longer holds exactly one member cannot be wrapped, which
            // is what makes a dropped member visible here.
            let wrap = |v: Vec<HgvsVariant>| match v.len() {
                1 => Some(v.into_iter().next().expect("length checked")),
                _ => None,
            };
            let (Some(lhs), Some(rhs)) = (wrap(before), wrap(after)) else {
                panic!("[{case}] the stage changed the member count on a single-member list");
            };
            let verdict = compare_denoted_sequences(&lhs, &rhs, &provider);
            assert!(
                matches!(verdict, DenotedSequenceComparison::Agree),
                "[{case}] the comparison did not agree (and must not merely have \
                 declined): {verdict:?}"
            );
        }
    }

    /// The stage is reachable and does not panic on an empty member list.
    ///
    /// Cheap, and it pins the one shape the call site can hand it that the cases
    /// above do not build: `normalize_allele` reaches this seam for every cis allele
    /// it settles, including degenerate ones.
    #[test]
    fn an_empty_member_list_is_accepted() {
        let mut empty: Vec<HgvsVariant> = Vec::new();
        render_canonical_spelling(&mut empty, &provider());
        assert!(empty.is_empty());
    }
}
