//! #1482 — a cis member spanning a region boundary was invisible to every
//! sibling-awareness pass, and the allele that fell out is one ferro's own
//! parser rejects.
//!
//! ```text
//! NM_TEST.1:c.[15_*1del;15_*1insCC]  ->  NM_TEST.1:c.[14_15dup;15_*1del]
//!   Self-cancelling allele: variants at index 1 (del) and 0 (dup) describe
//!   overlapping reference positions
//! ```
//!
//! # The issue was wrong about the cause, and that matters for what is pinned
//!
//! It reported the repair as "choosing a dup spelling whose implied insertion
//! point falls inside the sibling's deleted span". There is no such choice.
//! Measured, each member alone:
//!
//! ```text
//! c.15_*1insCC  ALONE -> c.14_15dup    (3')
//! c.15_*1del    ALONE -> c.15_*1del    unchanged
//! ```
//!
//! `c.14_15dup` is the *correct* normalization of that insertion in isolation —
//! `duplication.md:18` makes `dup` mandatory when a variant can be described as
//! one — and the allele output is simply the two per-member results side by
//! side. Input order makes no difference.
//!
//! The real cause is one line in `merge::join_pos`, which refuses an interval
//! whose endpoints are in different regions on the stated grounds that
//! cross-region ranges "have no valid HGVS syntax". **They do.** `c.15_*1del`
//! deletes across a stop codon, a common real variant, and the comment's own
//! example of the impossible thing — `c.-1_1` — is written out in
//! `consultation/SVD-WG001.md:37`. `member_span` inherited that refusal and
//! returned `None`, and every sibling-aware pass drops a `None` through
//! `filter_map` rather than declining. So the deletion was invisible **as a
//! sibling**, and `respell_colliding_duplications` — the pass that exists to
//! turn exactly that duplication back into an insertion, "a form that claims
//! nothing and so cannot collide" (#1320/#1323) — never saw the collision.
//!
//! The root cause is pinned where it lives, as unit tests on `member_span` in
//! `src/normalize/merge.rs`. This file pins the consequence end to end.
//!
//! # Every fixture here was checked to fail on the pre-fix revision
//!
//! Worth stating because the first cut of this file did not. Three of its five
//! tests passed unchanged on `a530cdaf` — the boundary shapes they used simply
//! never reached the collision on the core they were written against, so they
//! asserted a property that was already true and would have guarded nothing.
//! That is #1435's failure mode, and it is invisible unless you go and run the
//! test against the old code.
//!
//! What the sweep that replaced them found: the defect needs the insertion to
//! canonicalise to a duplication that *lands on* the sibling's span, which
//! depends on the bases either side of the boundary. So the fixtures below name
//! a specific core **and** a specific CDS for each case, and every one of them
//! produced an unparseable output before the fix. They are recorded next to each
//! assertion rather than left implicit.
//!
//! # None of the seam oracles saw it
//!
//! `FERRO_ASSERT_IN_BOUNDS` passes — every coordinate exists.
//! `FERRO_ASSERT_IDEMPOTENT` passes — the bad output is a fixed point.
//! `FERRO_ASSERT_REPARSE` is the one that *should* have fired, and the
//! rejection is exactly what it tests for; the shape simply was not in any
//! suite until #1482 filed it.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

/// The issue's own core. Under a `1..=15` CDS, `c.15` is `C` and `c.*1` is `T`.
const CORE: &str = "GCGCTAGTCTCGCCCTGTTA";

/// A non-coding record, for the `n.` axis cases at the bottom of this file.
const NONCODING_ACCESSION: &str = "NR_TEST.1";

/// A second core, needed for the 5'UTR/CDS boundary: on [`CORE`] the insertion
/// there never canonicalises onto the deletion's span, so that boundary cannot
/// reproduce the defect on it at all.
const AT_RICH_CORE: &str = "TTATTTAAATAAAAATAAAA";

fn provider(core: &str, cds_start: u64, cds_end: u64) -> MockProvider {
    let mut provider = MockProvider::new();
    let length = core.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        core.to_string(),
        Some(cds_start),
        Some(cds_end),
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

fn normalize(provider: &MockProvider, input: &str, direction: ShuffleDirection) -> String {
    let variant = parse_hgvs(input).expect("input must parse");
    Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::default().with_direction(direction),
    )
    .normalize(&variant)
    .expect("normalization must succeed")
    .to_string()
}

/// The property the issue is about: whatever spelling normalization settles on,
/// ferro must be able to read it back. Asserted rather than a pinned string, so
/// a later canonicalisation change cannot quietly re-break it while looking like
/// a re-blessing.
fn assert_output_reparses(provider: &MockProvider, input: &str, direction: ShuffleDirection) {
    let output = normalize(provider, input, direction);
    assert!(
        parse_hgvs(&output).is_ok(),
        "{input} normalized to {output}, which ferro's own parser rejects: {:?}",
        parse_hgvs(&output).err()
    );
}

/// Both members of the reported pair, normalized alone.
///
/// Pinned because the issue's stated cause hangs on it: if `c.15_*1insCC` alone
/// did *not* give `c.14_15dup`, the defect really would be a choice made by the
/// repair, and the fix would belong somewhere else entirely. Both values are
/// unchanged by this change — they were the same before it.
#[test]
fn each_member_alone_normalizes_as_the_spec_requires() {
    let provider = provider(CORE, 1, 15);
    assert_eq!(
        normalize(
            &provider,
            "NM_TEST.1:c.15_*1insCC",
            ShuffleDirection::ThreePrime
        ),
        "NM_TEST.1:c.14_15dup",
        "`duplication.md:18` makes `dup` mandatory here; the repair does not choose it"
    );
    assert_eq!(
        normalize(
            &provider,
            "NM_TEST.1:c.15_*1del",
            ShuffleDirection::ThreePrime
        ),
        "NM_TEST.1:c.15_*1del",
        "the deletion is already its own normal form"
    );
}

/// The reported allele, in both input orders.
///
/// Pre-fix, both orders gave `c.[14_15dup;15_*1del]` under a 3' shift — the
/// rejected output in the issue. The 5' direction was already fine on this core
/// (`c.[13_14dup;15_*1del]`), which is why the issue reports the defect as
/// 3'-only; it is asserted here anyway, since the fix must not cost the
/// direction that worked.
#[test]
fn the_reported_allele_reparses_in_either_input_order() {
    let provider = provider(CORE, 1, 15);
    for input in [
        "NM_TEST.1:c.[15_*1del;15_*1insCC]",
        "NM_TEST.1:c.[15_*1insCC;15_*1del]",
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert_output_reparses(&provider, input, direction);
        }
    }
}

/// The same collision at the **5'UTR/CDS** boundary, where the axis has no zero
/// and the crossing is `c.-1_1`.
///
/// Not merely a second instance: `region_sequence_delta` gives the 5'UTR its own
/// offset (`cds_start`, against the CDS body's `cds_start - 1`), so this
/// boundary exercises different arithmetic from the one above.
///
/// Pre-fix on this core with a `4..=15` CDS, both directions gave
/// `c.[-1dup;-1_1del]` — rejected as self-cancelling.
#[test]
fn the_five_prime_boundary_reparses() {
    let provider = provider(AT_RICH_CORE, 4, 15);
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_output_reparses(&provider, "NM_TEST.1:c.[-1_1del;-1_1insA]", direction);
    }
}

/// The stop codon moved, so the crossing member is `c.12_*1` rather than
/// `c.15_*1`.
///
/// Guards against a fix that happens to work only where the CDS ends at the
/// position the first fixture uses. Pre-fix, this gave `c.[11_12dup;12_*1del]`
/// under a 3' shift.
#[test]
fn the_defect_is_not_pinned_to_one_cds_end() {
    let provider = provider(CORE, 4, 15);
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_output_reparses(&provider, "NM_TEST.1:c.[12_*1del;12_*1insCC]", direction);
    }
}

/// A single-base insertion beside the same deletion.
///
/// A distinct shape rather than more of the same: a one-base `dup` sorts and
/// acts at the same position, so it reaches the collision test by the other of
/// the two routes `respell_colliding_duplications` checks. Pre-fix on this core
/// with a `1..=15` CDS it gave `c.[15dup;15_*1del]` (3') and
/// `c.[15_*1del;*1dup]` (5') — both rejected.
#[test]
fn a_single_base_insertion_reaches_the_same_collision() {
    let provider = provider(AT_RICH_CORE, 1, 15);
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_output_reparses(&provider, "NM_TEST.1:c.[15_*1del;15_*1insA]", direction);
    }
}

/// The overlap the parser **tolerates**, on the same boundary.
///
/// `parse_hgvs` rejects a `dup` overlapping a `del` and accepts one overlapping
/// an `inv`, so the `inv` form of the reported shape produced an overlapping
/// allele that no oracle complained about: pre-fix,
/// `c.[15_*1inv;15_*1insCC]` gave `c.[14_15dup;15_*1inv]`, two members claiming
/// `c.15`. It is fixed by the same change, and pinned here because the
/// reparse-based assertions above are structurally blind to it — that
/// inconsistency in the parser is worth its own issue and is not addressed here.
/// **Asserted as the contract, not as one exact string.** Both members sit at
/// `15_*1`, so the order between them is not determined by position and is an
/// artifact of which pass last touched the set. Pinning the rendered allele
/// verbatim made this test depend on that artifact: with #1511 (cis overlap
/// detection across a region boundary) in `main`, the same input renders
/// `c.[15_*1inv;15_*1insCC]` rather than `c.[15_*1insCC;15_*1inv]` — same two
/// members, opposite order, and the property this test exists for holds in
/// both. Landing order alone would otherwise turn it red.
///
/// What must hold is what the doc above states: the `insCC` is still an
/// **insertion**. Pre-fix it became `14_15dup`, and a `dup` claims `c.15`,
/// which the `inv` also claims — that is the overlap.
#[test]
fn the_inversion_form_no_longer_overlaps_its_sibling() {
    let provider = provider(CORE, 1, 15);
    let out = normalize(
        &provider,
        "NM_TEST.1:c.[15_*1inv;15_*1insCC]",
        ShuffleDirection::ThreePrime,
    );
    assert!(
        out.contains("15_*1insCC"),
        "the insertion must stay an insertion, which claims no bases and so \
         cannot overlap; got {out}"
    );
    assert!(
        out.contains("15_*1inv"),
        "the inversion member must survive unchanged; got {out}"
    );
    assert!(
        !out.contains("dup"),
        "a `dup` here claims `c.15` beside the `inv` that also claims it — the \
         overlap this test exists to rule out; got {out}"
    );
    // And the whole allele must still be one ferro's own parser accepts, which
    // is what an overlapping pair would fail.
    assert!(
        parse_hgvs(&out).is_ok(),
        "the normalized allele must re-parse; got {out}"
    );
}

/// The same two shapes wholly **inside** the CDS, pinned to the strings they
/// already produced.
///
/// The control for the whole file. This change makes a previously invisible span
/// visible, and the thing to rule out is that it also moved the spans that were
/// visible all along — so both values were measured on the pre-fix revision as
/// well as this one, and are identical there.
#[test]
fn an_in_region_pair_keeps_the_spelling_it_already_had() {
    let provider = provider(CORE, 1, 15);
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        assert_eq!(
            normalize(&provider, "NM_TEST.1:c.[11_12del;11_12insCC]", direction),
            "NM_TEST.1:c.12G>C"
        );
        assert_eq!(
            normalize(&provider, "NM_TEST.1:c.[11_12inv;11_12insCC]", direction),
            "NM_TEST.1:c.[11_12inv;11_12insCC]"
        );
    }
}

fn noncoding_provider() -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        NONCODING_ACCESSION.to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        CORE.to_string(),
        None,
        None,
        vec![Exon::new(1, 1, length)],
        None,
        None,
        None,
        Default::default(),
        ManeStatus::None,
        None,
        None,
    ));
    provider
}

/// `n.-N` and `n.*N` members stay visible to their siblings.
///
/// The counterweight to the whole change, and it is here because the first cut
/// of it got this wrong. Making `member_span` convert onto the sequence axis
/// means it now needs a conversion to exist, and `region_sequence_delta`
/// **refuses** the two `n.` regions outside the transcript — there are no bases
/// to read there. Routing `member_span` through it therefore *removed*
/// visibility for `n.-N` / `n.*N` members while adding it for boundary-crossing
/// ones: the exact defect this change is about, in a new place.
///
/// Both values below are `main`'s. Without `region_span_delta` they came out as
/// `n.[-6_-4del;-5=]` and `n.[-5=;-5del]` — the second being two members on one
/// position, which is what `drop_identity_members_covered_by_siblings` exists to
/// remove.
///
/// **The corpus cannot see any of this.** `dump_normalized_corpus` emits `c`,
/// `cx` and `g` rows and no `n.` rows at all, so its "0 newly wrong" is a
/// statement about the axes it builds. That is why this is pinned by hand.
#[test]
fn n_axis_members_outside_the_transcript_stay_visible_to_siblings() {
    let provider = noncoding_provider();
    for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
        // The identity is covered by the deletion, so it is dropped.
        assert_eq!(
            normalize(&provider, "NR_TEST.1:n.[-5=;-6_-4del]", direction),
            "NR_TEST.1:n.-6_-4del"
        );
        // Coincident, which is the shape that would otherwise have survived as
        // an allele claiming one position twice.
        assert_eq!(
            normalize(&provider, "NR_TEST.1:n.[-5=;-5del]", direction),
            "NR_TEST.1:n.-5del"
        );
        // The 3' side of the transcript, whose offset is the length rather than
        // a constant.
        assert_eq!(
            normalize(&provider, "NR_TEST.1:n.[*5=;*4_*6del]", direction),
            "NR_TEST.1:n.*4_*6del"
        );
        // In-transcript control.
        assert_eq!(
            normalize(&provider, "NR_TEST.1:n.[5=;4_6del]", direction),
            "NR_TEST.1:n.4_6del"
        );
    }
}
