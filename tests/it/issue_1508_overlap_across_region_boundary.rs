//! #1508 — both cis overlap detectors were blind to a member whose span
//! crosses a region boundary, so strict mode silently accepted overlaps it
//! rejects when the same pair sits inside one region.
//!
//! Measured on `a530cdaf`, a 20-base transcript with CDS `1..=15`:
//!
//! ```text
//! c.[11_12dup;12_13inv]      REJECTED  W5002
//! c.[14_15dup;15_*1inv]      ACCEPTED            <-- same shape, across the boundary
//! c.[11_12dup;12_13delinsAA] REJECTED  W5002
//! c.[14_15dup;15_*1delinsAA] ACCEPTED
//! c.[11_12del;11_12inv]      REJECTED  W5002
//! c.[15_*1del;15_*1inv]      ACCEPTED
//! c.[11_12insAA;11_12insCC]  REJECTED  W5002
//! c.[15_*1insAA;15_*1insCC]  ACCEPTED  -> c.[15delinsCAA;15_*1insCC]
//! ```
//!
//! Nothing about a pair changes when it moves across the boundary except that
//! one member's span stops being readable. `cds_range` / `tx_range` /
//! `rna_range` each ended with `if rs != re { return None }`, the same false
//! premise #1482 found in `merge::join_pos`, and `None` there is not a decline:
//! it reaches `let Some(..) = .. else { continue }`, so the member was dropped
//! from the analysis rather than treated conservatively.
//!
//! # Why the parser is not the layer at fault
//!
//! `parse_hgvs` rejects `c.[14_15dup;15_*1del]` with E3006 across the boundary,
//! correctly, because `detect_self_cancelling_pair` compares a region-ranked
//! `CanonicalPos`. E3006 is deliberately narrow — it implements one spec rule
//! (`general.md:58`, example `c.[762_768del;767_774dup]`) — and it implements it
//! on both sides of a boundary. The general "two members claim the same
//! positions" case is W5002's, and W5002 is the one that could not see across.
//!
//! # Blast radius
//!
//! This makes strict mode reject inputs it used to accept, so it was measured
//! rather than assumed. Over 76,638 two-member cis alleles — 3 cores × 3 CDS
//! layouts × every adjacent position pair spanning both boundaries × 5 edit
//! spellings × both shuffle directions:
//!
//! | | count |
//! |---|---|
//! | strict-mode acceptances, before | 68,624 |
//! | strict-mode acceptances, after | 67,497 |
//! | **newly rejected** | **1,127** |
//! | newly accepted | 0 |
//!
//! All 1,127 are alleles that denote no sequence: 1,062 have inputs the
//! independent applier already declines, and for the remaining 65 the input
//! applies but `main`'s own lenient *output* does not — checked one by one with
//! `canonical_spdi`, 65 of 65. So the change reports real overlaps and invents
//! none.

use ferro_hgvs::reference::transcript::{Exon, ManeStatus, Strand, Transcript};
use ferro_hgvs::{parse_hgvs, MockProvider, NormalizeConfig, Normalizer, ShuffleDirection};

const CORE: &str = "GCGCTAGTCTCGCCCTGTTA";

/// A transcript with a real 5'UTR *and* a real 3'UTR, so both boundaries exist.
fn provider(cds_start: u64, cds_end: u64) -> MockProvider {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        "NM_TEST.1".to_string(),
        Some("TESTGENE".to_string()),
        Strand::Plus,
        CORE.to_string(),
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

/// Whether strict mode accepts `input` — where "rejects" means **specifically
/// the overlap**, not merely that something went wrong.
///
/// A bare `.is_ok()` would collapse every error to `false`, so a row asserting
/// "this must be rejected" would be satisfied by an unrelated failure — a
/// provider miss, a parse-level refusal, a coordinate check. That is exactly
/// the outcome these tests must not accept, because the claim they defend is
/// that the *overlap detector* now sees across the boundary. Any other error
/// panics rather than quietly reading as a rejection.
fn strict_accepts(provider: &MockProvider, input: &str, direction: ShuffleDirection) -> bool {
    let variant = parse_hgvs(input).expect("input must parse");
    match Normalizer::with_config(
        provider.clone(),
        NormalizeConfig::strict().with_direction(direction),
    )
    .normalize(&variant)
    {
        Ok(_) => true,
        Err(e) => {
            let msg = e.to_string();
            assert!(
                msg.contains("W5002") || msg.contains("overlap"),
                "`{input}` was rejected, but not for the overlap this test is about: {msg}"
            );
            false
        }
    }
}

/// The four shapes from the issue, each paired with the in-region control that
/// was already rejected.
///
/// Written as pairs rather than as two lists because the control is the whole
/// argument: these are the *same* geometries, and the only thing that changed
/// was which side of the stop codon they sit on.
#[test]
fn an_overlap_is_reported_on_both_sides_of_the_cds_boundary() {
    let provider = provider(1, 15);
    for (across, inside) in [
        (
            "NM_TEST.1:c.[14_15dup;15_*1inv]",
            "NM_TEST.1:c.[11_12dup;12_13inv]",
        ),
        (
            "NM_TEST.1:c.[14_15dup;15_*1delinsAA]",
            "NM_TEST.1:c.[11_12dup;12_13delinsAA]",
        ),
        (
            "NM_TEST.1:c.[15_*1del;15_*1inv]",
            "NM_TEST.1:c.[11_12del;11_12inv]",
        ),
        (
            "NM_TEST.1:c.[15_*1insAA;15_*1insCC]",
            "NM_TEST.1:c.[11_12insAA;11_12insCC]",
        ),
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert!(
                !strict_accepts(&provider, inside, direction),
                "{inside} was already rejected before this change and must stay rejected"
            );
            assert!(
                !strict_accepts(&provider, across, direction),
                "{across} is the same geometry across the boundary and must be rejected too"
            );
        }
    }
}

/// The 5'UTR/CDS boundary, which has its own arithmetic: the axis has no zero,
/// so `c.-1` and `c.1` are adjacent while their coordinates differ by 2.
#[test]
fn the_five_prime_boundary_is_covered_too() {
    let provider = provider(4, 15);
    for input in [
        "NM_TEST.1:c.[-2_-1dup;-1_1inv]",
        "NM_TEST.1:c.[-1_1del;-1_1inv]",
        "NM_TEST.1:c.[-1_1insAA;-1_1insCC]",
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert!(
                !strict_accepts(&provider, input, direction),
                "{input} overlaps across the 5'UTR/CDS boundary"
            );
        }
    }
}

/// Members that merely *neighbour* a boundary, and members on opposite sides of
/// one, are still accepted.
///
/// The control that matters most for a change that makes strict mode stricter:
/// ranking the regions puts every member on one axis, and the risk is that
/// pairs which were invisible to each other become *falsely* visible. Each of
/// these is disjoint, and each is accepted before and after.
#[test]
fn disjoint_members_across_a_boundary_are_still_accepted() {
    let provider = provider(1, 15);
    for input in [
        // A CDS member and a 3'UTR member with a gap between them.
        "NM_TEST.1:c.[11_12del;*3_*4inv]",
        // Flush on either side of the boundary, but not overlapping.
        "NM_TEST.1:c.[14_15del;*1_*2dup]",
        // A span ending at the boundary beside an insertion well 5' of it.
        "NM_TEST.1:c.[13_14insAA;15_*1del]",
        // Two members entirely inside the 3'UTR.
        "NM_TEST.1:c.[*1_*2del;*4_*5inv]",
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert!(
                strict_accepts(&provider, input, direction),
                "{input} has disjoint members and must stay accepted"
            );
        }
    }
}

/// A non-coding `n.` record, where the ranked regions are
/// `TxUpstream` / `Tx` / `TxDownstream` rather than the UTRs.
///
/// `tx_range` carried the same refusal as `cds_range`, and the corpus harness
/// emits no `n.` rows at all, so this axis is only covered by hand.
#[test]
fn the_transcript_axis_boundaries_are_covered() {
    let mut provider = MockProvider::new();
    let length = CORE.len() as u64;
    provider.add_transcript(Transcript::new(
        "NR_TEST.1".to_string(),
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
    for input in [
        "NR_TEST.1:n.[19_20dup;20_*1inv]",
        "NR_TEST.1:n.[20_*1del;20_*1inv]",
        "NR_TEST.1:n.[-1_1del;-1_1inv]",
    ] {
        for direction in [ShuffleDirection::ThreePrime, ShuffleDirection::FivePrime] {
            assert!(
                !strict_accepts(&provider, input, direction),
                "{input} overlaps across an `n.` region boundary"
            );
        }
    }
}
