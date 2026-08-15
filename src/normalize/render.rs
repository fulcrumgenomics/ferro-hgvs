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
//! ablation over an earlier **86,170-row** revision of the synthetic corpus — the
//! proportion is the argument, and the corpus is 96,182 rows on this branch at base
//! `5567412f`, so read these as a share rather than as current counts — those four
//! move **1,279**, **1,011**, **0** and **505** rows respectively, a union of
//! **2,761 (3.2%)**
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
/// corpus has moved four times in this campaign alone — 85,642 to 85,930 to 86,170,
/// and then to 96,182 when a rebase brought in the families other branches had been
/// adding the whole time — so a literal pinned to one of those is stale by the next
/// family, and this comment has had to be re-pointed once already. That is the failure
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

/// Which of the per-member pipeline's two reference shapes a [`MintReference`]
/// rebuilds (#1946).
///
/// This is what the reference gate keys its comparison on. It is derived from the
/// member's `CisKind` and from nothing else —
/// an earlier revision inferred the same distinction from the *contents* of the
/// rebuilt reference (`complete && offset == 0`), which happened to agree only
/// because a genomic window's padding is narrower than the transcript-flanking
/// pad, and would have started disagreeing the moment either constant moved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MintFrame {
    /// `c.`/`r.`/`n.` — the entire transcript, exactly as `normalize_cds`,
    /// `normalize_rna` and `normalize_tx` serve it. There is no window to choose.
    WholeTranscript,
    /// `g.`/`m.` — a window around the member's own span.
    GenomicWindow,
}

/// The reference a render-time mint is evaluated against, and the frame its
/// positions live in (#1946).
///
/// `ref_seq` is **not one thing** in the per-member pipeline — it is the entire
/// spliced transcript on `c.`/`r.`/`n.`, a growing window on `g.`/`m.`, and a
/// reverse-complemented genomic window on the intronic paths — so a stage that
/// receives only a settled member list has to rebuild it deliberately rather than
/// substitute one uniform window. Doing the latter silently truncates tract
/// discovery near an edge: `insertion_to_repeat` walks outward bounded only by
/// `ref_seq.len()`.
#[derive(Debug, Clone)]
pub struct MintReference {
    /// Reference bases, **exactly as the pipeline's own fetch produces them**.
    ///
    /// Deliberately not case-normalized here: byte-equality with what
    /// `normalize_na_edit` was handed is the property the gate asserts, so this
    /// has to inherit each path's own convention rather than impose one.
    /// `fetch_canonical_window` upper-cases, so [`MintFrame::GenomicWindow`]
    /// bases are upper-case; a transcript's bases are served verbatim, so
    /// [`MintFrame::WholeTranscript`] bases are whatever the provider holds.
    /// (This field's doc claimed "upper-cased" flatly, which was true of neither
    /// the transcript branch nor of the growth loop it is compared against —
    /// `normalize_in_grown_window` passes `provider.get_sequence(..).as_bytes()`
    /// with no case fold at all.)
    pub bases: Vec<u8>,
    /// Which shape this is. Ask this, not the flags below, when the question is
    /// "which pipeline call does this model".
    pub frame: MintFrame,
    /// A 1-based **sequence** position `p` indexes `bases[(p - 1 - offset)]`.
    /// Zero when `bases` starts at the sequence's own first base.
    pub offset: i64,
    /// The 5' edge of `bases` is the sequence's own start, so there is nothing
    /// further 5' for a tract walk to reach.
    pub at_sequence_start: bool,
    /// The 3' edge of `bases` is the sequence's own end, so there is nothing
    /// further 3' for a tract walk to reach.
    ///
    /// Separate from [`Self::at_sequence_start`] on purpose: the predicate this
    /// replaces (`complete`) was set from a 3'-clamped fetch alone and then
    /// documented as "`bases` is the complete reference", which is a strictly
    /// stronger claim than the code checked — a window clamped only at the
    /// contig end reported itself complete while its 5' side sat in the interior.
    pub at_sequence_end: bool,
}

impl MintReference {
    /// Whether `bases` is the entire reference sequence, so no edge can truncate
    /// a tract walk and there is nothing to grow into on either side.
    #[must_use]
    pub fn is_whole_sequence(&self) -> bool {
        self.at_sequence_start && self.at_sequence_end
    }
}

/// Half-width of the first genomic window, matching `NormalizeConfig`'s default.
pub const MINT_WINDOW: i64 = 100;

/// Rebuild the reference a settled `member`'s mint should be evaluated against.
///
/// **Unused by `render_canonical_spelling` on purpose.** This lands with its
/// equivalence test and nothing else: the gate on relocating the repeat mint is
/// whether a reference rebuilt here matches what the per-member pipeline was given,
/// and that question is answered before any mint moves.
///
/// Mirrors the existing callers rather than unifying them:
///
/// - **`c.`/`r.`/`n.`** — the whole transcript, exactly as `normalize_cds`,
///   `normalize_rna` and `normalize_tx` pass it. [`MintReference::is_whole_sequence`]
///   holds, so no growth is possible or needed.
/// - **`g.`/`m.`** — a window, whose edges are flagged individually. Growth is the
///   caller's business and is **triggered differently from
///   `normalize_in_grown_window`**: that loop grows when the *shuffle* runs to an
///   edge, and at render time there is no shuffle, so a render-time mint must grow
///   when the *tract walk* reaches one.
///
/// # What this does NOT model, and why the gate has to say so
///
/// A `c.` member whose span crosses an exon/exon junction, or which sits in an
/// intron, is **not** normalized against the spliced transcript: `normalize_cds`
/// hands it to `normalize_boundary_spanning_cds` or `normalize_intronic_cds`, which
/// fetch a *genomic* window and shuffle there. This function returns the whole
/// transcript for such a member, because that is what its `CisKind` says — so for
/// those members it rebuilds a reference the pipeline's junction-crossing call was
/// never given.
///
/// That is a real gap in this function, and the gate reports it as its own
/// population rather than letting it hide inside a match: see
/// `the_render_time_reference_matches_what_the_pipeline_was_given`, which selects
/// the recorded calls by [`NaEditCallSite`](crate::normalize::NaEditCallSite) and
/// counts members whose depth-0 calls were all at a site this function does not
/// model. Closing it means minting the same genomic window those two paths build,
/// and is the next step rather than a widening to do in passing.
///
/// # Measured over `dump(1)`
///
/// Comparing what this rebuilds against what `normalize_na_edit` was actually
/// handed (recorded by `normalize::take_na_edit_references`), selecting by call
/// site — the figures are printed by the gate and are quoted in its doc comment
/// rather than here, so that one measurement has one home. What is worth stating
/// on the function is the shape of the answer: the whole-transcript frame is
/// asserted **byte-identical**, the genomic frame is **reported** and cannot be
/// byte-equal by construction (the pipeline's window is centred on the *input*
/// position, and it starts one base 3' of this one — `normalize_in_grown_window`
/// fetches from `start - w` where this fetches from `start - 1 - w`).
pub fn mint_reference_for<P: ReferenceProvider>(
    member: &HgvsVariant,
    provider: &P,
    half_width: i64,
) -> Option<MintReference> {
    use crate::normalize::merge::{
        cis_axis_parts, cis_kind_of, fetch_canonical_window, region_sequence_delta, CisKind,
    };

    let kind = cis_kind_of(member)?;
    let (accession, region, axis_start, axis_end, _edit) = cis_axis_parts(member, kind)?;
    let key = accession.transcript_accession();

    // Transcript axes: the reference is served whole, so there is no window to
    // choose. Keyed on the kind, which is what decides the frame — never on the
    // shape of what came back.
    if !matches!(kind, CisKind::Genome | CisKind::Mt) {
        let transcript = provider.get_transcript_for_variant(member).ok()?;
        let bases = transcript.sequence.as_deref()?.as_bytes().to_vec();
        return Some(MintReference {
            bases,
            frame: MintFrame::WholeTranscript,
            offset: 0,
            at_sequence_start: true,
            at_sequence_end: true,
        });
    }

    // Genomic: a window around the member's own span, in sequence coordinates.
    let delta = region_sequence_delta(region, &key, provider)?;
    let (lo, hi) = (
        axis_start.min(axis_end) + delta,
        axis_start.max(axis_end) + delta,
    );
    let start0 = (lo - 1 - half_width).max(0);
    let end0 = hi + half_width;
    if end0 <= start0 {
        return None;
    }
    let bases = fetch_canonical_window(provider, &key, start0, end0)?;
    // A short read means the fetch ran off the contig, which is the only way to
    // learn the 3' edge without a second provider round-trip. The 5' edge is
    // decidable from `start0` alone.
    let at_sequence_end = (bases.len() as i64) < end0 - start0;
    Some(MintReference {
        bases,
        frame: MintFrame::GenomicWindow,
        offset: start0,
        at_sequence_start: start0 == 0,
        at_sequence_end,
    })
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
    /// `--verify-spdi` pass, which compared 86,170 rows at the time and reported the
    /// seam moving none. The seam's no-op-ness is re-established at each base rather
    /// than carried: at `5567412f` the branch's dump with its three added families
    /// stripped is `diff`-identical to the base's own, over all 95,614 shared rows.
    /// Saying that plainly is the point — listing a case the oracle cannot
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
